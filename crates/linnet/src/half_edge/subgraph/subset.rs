use std::{
    cmp::Ordering,
    hash::{Hash, Hasher},
    marker::PhantomData,
    ops::{Index, Not, Range, RangeFrom, RangeInclusive, RangeTo, RangeToInclusive},
};

use bitvec::{bitvec, order::Lsb0, vec::BitVec};
#[cfg(feature = "rkyv")]
use rkyv::{
    with::{ArchiveWith, DeserializeWith, SerializeWith},
    Archive, Archived, Deserialize, Fallible, Resolver, Serialize,
};

use crate::half_edge::{
    involution::{Flow, Hedge, HedgePair},
    nodestore::{NodeStorage, NodeStorageOps},
    subgraph::{
        BaseSubgraph, Inclusion, ModifySubSet, SubGraphLike, SubGraphOps, SubSetLike, SubSetOps,
        BASE62_ALPHABET,
    },
    swap::Swap,
    typed_vec::IndexLike,
    HedgeGraph,
};

pub type SuBitGraph = SubSet<Hedge>;

#[cfg(feature = "rkyv")]
#[derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)]
#[archive(check_bytes)]
pub(crate) struct BitVecArchiveRepr {
    raw: Vec<usize>,
    bit_len: usize,
}

#[cfg(feature = "rkyv")]
impl From<&BitVec> for BitVecArchiveRepr {
    fn from(value: &BitVec) -> Self {
        Self {
            raw: value.as_raw_slice().to_vec(),
            bit_len: value.len(),
        }
    }
}

#[cfg(feature = "rkyv")]
struct BitVecRkyv;

#[cfg(feature = "rkyv")]
impl ArchiveWith<BitVec> for BitVecRkyv {
    type Archived = Archived<BitVecArchiveRepr>;
    type Resolver = Resolver<BitVecArchiveRepr>;

    unsafe fn resolve_with(
        field: &BitVec,
        pos: usize,
        resolver: Self::Resolver,
        out: *mut Self::Archived,
    ) {
        BitVecArchiveRepr::from(field).resolve(pos, resolver, out);
    }
}

#[cfg(feature = "rkyv")]
impl<S: Fallible + ?Sized> SerializeWith<BitVec, S> for BitVecRkyv
where
    BitVecArchiveRepr: Serialize<S>,
{
    fn serialize_with(field: &BitVec, serializer: &mut S) -> Result<Self::Resolver, S::Error> {
        BitVecArchiveRepr::from(field).serialize(serializer)
    }
}

#[cfg(feature = "rkyv")]
impl<D: Fallible + ?Sized> DeserializeWith<Archived<BitVecArchiveRepr>, BitVec, D> for BitVecRkyv
where
    Archived<BitVecArchiveRepr>: Deserialize<BitVecArchiveRepr, D>,
{
    fn deserialize_with(
        field: &Archived<BitVecArchiveRepr>,
        deserializer: &mut D,
    ) -> Result<BitVec, D::Error> {
        let repr = field.deserialize(deserializer)?;
        let mut bitvec = BitVec::from_vec(repr.raw);
        bitvec.truncate(repr.bit_len);
        Ok(bitvec)
    }
}

const WORD_BITS: usize = usize::BITS as usize;

/// Mask of the live bits in the last raw word of a `len`-bit set.
///
/// `bitvec` leaves the dead bits past `len` unspecified (`!` and `truncate`
/// leave them set), so every word-level read of a [`SubSet`] masks its last word.
fn tail_mask(len: usize) -> usize {
    usize::MAX >> ((WORD_BITS - len % WORD_BITS) % WORD_BITS)
}

/// A set of IDs `0..size` stored as a bitset.
///
/// Equality, hashing, ordering and the set algebra work on whole words. The
/// order is lexicographic from ID 0 (the first ID in only one set decides,
/// and the set that has it is greater), with a strict prefix before its
/// extensions, i.e. the order of `bitvec`'s `BitSlice`.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
#[cfg_attr(feature = "rkyv", archive(check_bytes))]
pub struct SubSet<ID> {
    #[cfg_attr(feature = "bincode", bincode(with_serde))]
    #[cfg_attr(feature = "rkyv", with(BitVecRkyv))]
    set: BitVec,
    id: std::marker::PhantomData<ID>,
}

impl<ID> PartialEq for SubSet<ID> {
    fn eq(&self, other: &Self) -> bool {
        self.set.len() == other.set.len() && self.words().eq(other.words())
    }
}

impl<ID> Eq for SubSet<ID> {}

impl<ID> Hash for SubSet<ID> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_usize(self.set.len());
        self.words().for_each(|word| state.write_usize(word));
    }
}

impl<ID> PartialOrd for SubSet<ID> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<ID> Ord for SubSet<ID> {
    fn cmp(&self, other: &Self) -> Ordering {
        let common = self.set.len().min(other.set.len());
        let n_words = common.div_ceil(WORD_BITS);
        let words = self.raw()[..n_words].iter().zip(&other.raw()[..n_words]);
        for (i, (&a, &b)) in words.enumerate() {
            let diff = if i + 1 == n_words {
                (a ^ b) & tail_mask(common)
            } else {
                a ^ b
            };
            if diff != 0 {
                let first = diff.isolate_lowest_one();
                return if a & first != 0 {
                    Ordering::Greater
                } else {
                    Ordering::Less
                };
            }
        }
        self.set.len().cmp(&other.set.len())
    }
}

impl<ID> Not for &SubSet<ID> {
    type Output = SubSet<ID>;
    fn not(self) -> Self::Output {
        Self::Output {
            set: self.set.clone().not(),
            id: std::marker::PhantomData,
        }
    }
}

impl<ID> FromIterator<bool> for SubSet<ID> {
    fn from_iter<T: IntoIterator<Item = bool>>(iter: T) -> Self {
        let set: BitVec = iter.into_iter().collect();
        Self {
            set,
            id: std::marker::PhantomData,
        }
    }
}

impl<ID: IndexLike> Swap<ID> for SubSet<ID> {
    fn is_zero_length(&self) -> bool {
        self.set.is_empty()
    }

    fn len(&self) -> ID {
        self.set.len().into()
    }

    fn swap(&mut self, i: ID, j: ID) {
        self.set.swap(i.into(), j.into());
    }
}

impl<ID> SubSet<ID> {
    pub fn is_full(&self) -> bool {
        self.set.all()
    }

    pub fn iter(&self) -> impl Iterator<Item = bool> + '_ {
        self.set.iter().by_vals()
    }
    pub fn split_off(&mut self, at: ID) -> Self
    where
        ID: IndexLike,
    {
        let mut set = self.set.split_off(at.into());
        // `split_off` keeps the offset of `at` within its word.
        set.force_align();
        Self {
            set,
            id: PhantomData,
        }
    }
    pub fn full(size: usize) -> Self {
        Self {
            set: bitvec![usize, Lsb0; 1; size],
            id: std::marker::PhantomData,
        }
    }

    pub fn push(&mut self, val: bool) {
        self.set.push(val);
    }

    pub fn from_usize(num: usize, size: usize) -> Self {
        let mut set = BitVec::from_element(num);
        set.truncate(size);
        Self {
            set,
            id: std::marker::PhantomData,
        }
    }

    pub fn pop(&mut self) -> Option<bool> {
        self.set.pop()
    }

    pub fn clear(&mut self) {
        self.set.fill(false);
    }

    /// The raw words. A set starts at bit 0 of its first word, so bit `i` of
    /// word `w` is ID `w * WORD_BITS + i`; the dead bits past the last ID are
    /// unspecified.
    fn raw(&self) -> &[usize] {
        debug_assert_eq!(
            self.set.as_bitptr().bit().into_inner(),
            0,
            "unaligned SubSet"
        );
        self.set.as_raw_slice()
    }

    fn raw_mut(&mut self) -> &mut [usize] {
        debug_assert_eq!(
            self.set.as_bitptr().bit().into_inner(),
            0,
            "unaligned SubSet"
        );
        self.set.as_raw_mut_slice()
    }

    /// The raw words, with the dead bits past the last ID cleared.
    fn words(&self) -> impl Iterator<Item = usize> + '_ {
        let mask = tail_mask(self.set.len());
        let (last, body) = match self.raw().split_last() {
            Some((&last, body)) => (Some(last & mask), body),
            None => (None, &[][..]),
        };
        body.iter().copied().chain(last)
    }

    /// Replaces each word `a` of `self` by `op(a, b)` with `b` the matching word
    /// of `other`, zero-extended or truncated to the length of `self`.
    fn zip_words_with(&mut self, other: &Self, op: impl Fn(usize, usize) -> usize) {
        let mut other = other.words();
        for word in self.raw_mut() {
            *word = op(*word, other.next().unwrap_or(0));
        }
    }
}

/// An iterator over the IDs included in a [`SubSet`], ascending from the front
/// and descending from the back.
///
/// It consumes the lowest set bit of the front word or the highest set bit of
/// the back word, and refills either from the untouched words in between.
#[derive(Clone, Debug)]
pub struct SubSetIter<'a, ID = Hedge> {
    /// Remaining bits of the front word, whose bit 0 is ID `front_base`.
    front: usize,
    front_base: usize,
    /// Untouched words strictly between the front and the back word.
    middle: &'a [usize],
    /// Remaining bits of the back word, whose bit 0 is ID `back_base`.
    back: usize,
    back_base: usize,
    id: PhantomData<fn() -> ID>,
}

impl<'a, ID> SubSetIter<'a, ID> {
    fn new(set: &'a SubSet<ID>) -> Self {
        let mask = tail_mask(set.set.len());
        let words = set.raw();
        let (front, middle, back) = match words {
            [] => (0, &[][..], 0),
            [only] => (0, &[][..], only & mask),
            [first, middle @ .., last] => (*first, middle, last & mask),
        };
        Self {
            front,
            front_base: 0,
            middle,
            back,
            back_base: words.len().saturating_sub(1) * WORD_BITS,
            id: PhantomData,
        }
    }
}

impl<ID: From<usize>> Iterator for SubSetIter<'_, ID> {
    type Item = ID;
    fn next(&mut self) -> Option<Self::Item> {
        while self.front == 0 {
            if let Some((&word, rest)) = self.middle.split_first() {
                self.front = word;
                self.front_base += WORD_BITS;
                self.middle = rest;
            } else if self.back != 0 {
                self.front = std::mem::take(&mut self.back);
                self.front_base = self.back_base;
            } else {
                return None;
            }
        }
        let bit = self.front.trailing_zeros() as usize;
        self.front &= self.front - 1;
        Some(ID::from(self.front_base + bit))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.len();
        (len, Some(len))
    }
}

impl<ID: From<usize>> DoubleEndedIterator for SubSetIter<'_, ID> {
    fn next_back(&mut self) -> Option<Self::Item> {
        while self.back == 0 {
            if let Some((&word, rest)) = self.middle.split_last() {
                self.back = word;
                self.back_base -= WORD_BITS;
                self.middle = rest;
            } else if self.front != 0 {
                self.back = std::mem::take(&mut self.front);
                self.back_base = self.front_base;
            } else {
                return None;
            }
        }
        let bit = WORD_BITS - 1 - self.back.leading_zeros() as usize;
        self.back ^= 1 << bit;
        Some(ID::from(self.back_base + bit))
    }
}

impl<ID: From<usize>> ExactSizeIterator for SubSetIter<'_, ID> {
    fn len(&self) -> usize {
        let middle: u32 = self.middle.iter().map(|word| word.count_ones()).sum();
        (self.front.count_ones() + middle + self.back.count_ones()) as usize
    }
}

#[cfg(feature = "rkyv")]
impl<ID: IndexLike> ArchivedSubSet<ID> {
    pub fn contains(&self, index: ID) -> bool {
        self.contains_index(index.into())
    }

    pub fn contains_index(&self, index: usize) -> bool {
        if index >= self.size() {
            return false;
        }

        let bits_per_word = usize::BITS as usize;
        let word_index = index / bits_per_word;
        let bit_index = index % bits_per_word;
        self.set
            .raw
            .as_slice()
            .get(word_index)
            .is_some_and(|word| ((*word >> bit_index) & 1) == 1)
    }

    pub fn size(&self) -> usize {
        self.set.bit_len.try_into().unwrap()
    }

    pub fn intersects_owned(&self, other: &SubSet<ID>) -> bool {
        other.included_iter().any(|index| self.contains(index))
    }
}

impl<ID> Not for SubSet<ID> {
    type Output = Self;
    fn not(self) -> Self::Output {
        Self {
            set: self.set.not(),
            id: std::marker::PhantomData,
        }
    }
}
impl<ID: Into<usize>> Index<ID> for SubSet<ID> {
    type Output = bool;

    fn index(&self, index: ID) -> &Self::Output {
        self.set.index(index.into())
    }
}
impl Inclusion<HedgePair> for SuBitGraph {
    fn includes(&self, other: &HedgePair) -> bool {
        match other {
            HedgePair::Unpaired { hedge, .. } => self.includes(hedge),
            HedgePair::Split {
                source,
                sink,
                split,
            } => match split {
                Flow::Sink => self.includes(sink),
                Flow::Source => self.includes(source),
            },
            HedgePair::Paired { source, sink } => self.includes(source) && self.includes(sink),
        }
    }

    fn intersects(&self, other: &HedgePair) -> bool {
        match other {
            HedgePair::Unpaired { hedge, .. } => self.includes(hedge),
            HedgePair::Split {
                source,
                sink,
                split,
            } => match split {
                Flow::Sink => self.includes(sink),
                Flow::Source => self.includes(source),
            },
            HedgePair::Paired { source, sink } => self.includes(source) || self.includes(sink),
        }
    }
}

// impl<ID> Inclusion<SubSet<ID>> for SubSet<ID> {
//     fn includes(&self, other: &SubSet<ID>) -> bool {
//         &self.intersection(other) == other
//     }

//     fn intersects(&self, other: &SubSet<ID>) -> bool {
//         self.intersection(other).is_empty()
//     }
// }

impl<ID: IndexLike> Inclusion<Range<ID>> for SubSet<ID> {
    fn includes(&self, other: &Range<ID>) -> bool {
        (other.start.into()..other.end.into()).all(|a| self.includes(&ID::from(a)))
    }

    fn intersects(&self, other: &Range<ID>) -> bool {
        (other.start.into()..other.end.into()).any(|a| self.includes(&ID::from(a)))
    }
}

impl<ID: IndexLike> Inclusion<SubSet<ID>> for SubSet<ID> {
    fn includes(&self, other: &SubSet<ID>) -> bool {
        self.set.len() == other.set.len()
            && self.words().zip(other.words()).all(|(a, b)| b & !a == 0)
    }

    fn intersects(&self, other: &SubSet<ID>) -> bool {
        self.words().zip(other.words()).any(|(a, b)| a & b != 0)
    }
}

impl<ID: IndexLike> Inclusion<RangeTo<ID>> for SubSet<ID> {
    fn includes(&self, other: &RangeTo<ID>) -> bool {
        (0..other.end.into()).all(|a| self.includes(&ID::from(a)))
    }

    fn intersects(&self, other: &RangeTo<ID>) -> bool {
        (0..other.end.into()).any(|a| self.includes(&ID::from(a)))
    }
}

impl<ID: IndexLike> Inclusion<RangeToInclusive<ID>> for SubSet<ID> {
    fn includes(&self, other: &RangeToInclusive<ID>) -> bool {
        (0..=other.end.into()).all(|a| self.includes(&ID::from(a)))
    }

    fn intersects(&self, other: &RangeToInclusive<ID>) -> bool {
        (0..=other.end.into()).any(|a| self.includes(&ID::from(a)))
    }
}

impl<ID: IndexLike> Inclusion<RangeFrom<ID>> for SubSet<ID> {
    fn includes(&self, other: &RangeFrom<ID>) -> bool {
        (other.start.into()..).all(|a| self.includes(&ID::from(a)))
    }

    fn intersects(&self, other: &RangeFrom<ID>) -> bool {
        (other.start.into()..).any(|a| self.includes(&ID::from(a)))
    }
}

impl<ID: IndexLike> Inclusion<RangeInclusive<ID>> for SubSet<ID> {
    fn includes(&self, other: &RangeInclusive<ID>) -> bool {
        ((*other.start()).into()..=(*other.end()).into()).all(|a| self.includes(&ID::from(a)))
    }

    fn intersects(&self, other: &RangeInclusive<ID>) -> bool {
        ((*other.start()).into()..=(*other.end()).into()).any(|a| self.includes(&ID::from(a)))
    }
}

impl<ID: IndexLike> Inclusion<ID> for SubSet<ID> {
    fn includes(&self, other: &ID) -> bool {
        self[*other]
    }

    fn intersects(&self, other: &ID) -> bool {
        self.includes(other)
    }
}

impl<ID: IndexLike> ModifySubSet<ID> for SubSet<ID> {
    fn add(&mut self, hedge: ID) {
        self.set.set(hedge.into(), true);
    }

    fn sub(&mut self, hedge: ID) {
        self.set.set(hedge.into(), false);
    }
}

impl ModifySubSet<HedgePair> for SuBitGraph {
    fn add(&mut self, index: HedgePair) {
        match index {
            HedgePair::Paired { source, sink } => {
                self.add(source);
                self.add(sink);
            }
            HedgePair::Split {
                source,
                sink,
                split,
            } => match split {
                Flow::Source => {
                    self.add(source);
                    self.sub(sink);
                }
                Flow::Sink => {
                    self.add(sink);
                    self.sub(source);
                }
            },
            HedgePair::Unpaired { hedge, .. } => {
                self.add(hedge);
            }
        }
    }

    fn sub(&mut self, index: HedgePair) {
        match index {
            HedgePair::Paired { source, sink } => {
                self.sub(source);
                self.sub(sink);
            }
            HedgePair::Split {
                source,
                sink,
                split,
            } => match split {
                Flow::Source => {
                    self.sub(source);
                    // self.sub(sink);
                }
                Flow::Sink => {
                    self.sub(sink);
                    // self.sub(source);
                }
            },
            HedgePair::Unpaired { hedge, .. } => {
                self.sub(hedge);
            }
        }
    }
}

impl BaseSubgraph for SuBitGraph {
    fn from_filter<E, V, H, N: NodeStorageOps<NodeData = V>, F: FnMut(&E) -> bool>(
        graph: &HedgeGraph<E, V, H, N>,
        mut filter: F,
    ) -> Self {
        let mut empty: SuBitGraph = graph.empty_subgraph();

        for (p, _, d) in graph.iter_edges() {
            if filter(d.data) {
                empty.add(p);
            }
        }

        empty
    }

    fn from_hedge_iter<I: Iterator<Item = Hedge>>(iter: I, len: usize) -> Self {
        let mut subgraph = SuBitGraph::empty(len);

        for h in iter {
            subgraph.add(h);
        }

        subgraph
    }
}

impl SubGraphLike for SuBitGraph {
    fn hairs(&self, node: impl Iterator<Item = Hedge>) -> SuBitGraph {
        let mut hairs = SuBitGraph::empty(self.size());

        for h in node {
            if self.includes(&h) {
                hairs.add(h);
            }
        }

        hairs
    }

    fn nedges<E, V, H, N: NodeStorageOps<NodeData = V>>(
        &self,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> usize {
        let mut count = 0;
        for i in self.included_iter() {
            if i != graph.inv(i) && self.includes(&graph.inv(i)) {
                count += 1;
            }
        }
        count / 2
    }
}

impl<ID: IndexLike> SubSetLike<ID> for SubSet<ID> {
    type Base = SubSet<ID>;
    type BaseIter<'a>
        = SubSetIter<'a, ID>
    where
        ID: 'a;
    fn included(&self) -> &Self::Base {
        self
    }

    fn size(&self) -> usize {
        self.set.len()
    }

    fn has_greater(&self, hedge: ID) -> bool {
        (hedge.into()..self.size()).any(|h| self.includes(&ID::from(h)))
    }

    fn join_mut(&mut self, other: Self) {
        self.set.extend(other.set);
    }
    fn included_iter(&self) -> Self::BaseIter<'_> {
        SubSetIter::new(self)
    }

    fn n_included(&self) -> usize {
        self.words().map(|word| word.count_ones() as usize).sum()
    }

    fn empty(size: usize) -> Self {
        Self {
            set: bitvec![usize, Lsb0; 0; size],
            id: PhantomData,
        }
    }

    fn string_label(&self) -> String {
        if self.is_empty() {
            return "0".to_string();
        }

        let mut digits = vec![0u8]; // Initialize with a single zero digit

        // Iterate over the bits from MSB to LSB
        for bit in self.set.iter().by_vals().rev() {
            let mut carry = 0u8;

            // Multiply existing digits by 2 (shift left)
            for digit in &mut digits {
                let temp = (*digit as u16) * 2 + carry as u16;
                *digit = (temp % 62) as u8;
                carry = (temp / 62) as u8;
            }

            if carry > 0 {
                digits.push(carry);
            }

            // Add the current bit (if it's 1)
            if bit {
                let mut carry = 1u8;
                for digit in &mut digits {
                    let temp = *digit as u16 + carry as u16;
                    *digit = (temp % 62) as u8;
                    carry = (temp / 62) as u8;

                    if carry == 0 {
                        break;
                    }
                }
                if carry > 0 {
                    digits.push(carry);
                }
            }
        }

        // Map digits to base62 characters and reverse the result
        let base62_string: String = digits
            .iter()
            .rev()
            .map(|&d| BASE62_ALPHABET[d as usize] as char)
            .collect();

        base62_string
    }

    fn from_base62(label: &str, size: usize) -> Option<Self> {
        let label = label.trim();
        if label.is_empty() {
            return None;
        }

        let mut digits: Vec<u8> = Vec::with_capacity(label.len());
        for byte in label.as_bytes() {
            let digit = BASE62_ALPHABET.iter().position(|&b| b == *byte)? as u8;
            digits.push(digit);
        }

        while digits.len() > 1 && digits.first() == Some(&0) {
            digits.remove(0);
        }

        let mut bits = Vec::new();
        while digits.iter().any(|&d| d != 0) {
            let mut carry = 0u8;
            let mut next = Vec::with_capacity(digits.len());
            for digit in digits {
                let value = (carry as u16) * 62 + digit as u16;
                let next_digit = (value / 2) as u8;
                carry = (value % 2) as u8;
                if !next.is_empty() || next_digit != 0 {
                    next.push(next_digit);
                }
            }
            bits.push(carry == 1);
            digits = next;
        }

        if bits.len() > size {
            return None;
        }

        let mut set = bitvec![usize, Lsb0; 0; size];
        for (index, bit) in bits.into_iter().enumerate() {
            if bit {
                set.set(index, true);
            }
        }

        Some(Self {
            set,
            id: PhantomData,
        })
    }

    fn is_empty(&self) -> bool {
        self.words().all(|word| word == 0)
    }
}

impl<ID: IndexLike> SubSetOps<ID> for SubSet<ID> {
    fn intersect_with(&mut self, other: &Self) {
        self.zip_words_with(other, |a, b| a & b)
    }

    fn union_with(&mut self, other: &Self) {
        self.zip_words_with(other, |a, b| a | b)
    }

    fn sym_diff_with(&mut self, other: &Self) {
        self.zip_words_with(other, |a, b| a ^ b)
    }

    fn empty_union(&self, other: &Self) -> bool {
        self.is_empty() && other.is_empty()
    }

    fn empty_intersection(&self, other: &Self) -> bool {
        !self.intersects(other)
    }

    fn subtract_with(&mut self, other: &Self) {
        self.zip_words_with(other, |a, b| a & !b)
    }
}

impl SubGraphOps for SuBitGraph {
    fn complement<E, V, H, N: NodeStorage<NodeData = V>>(
        &self,
        _graph: &HedgeGraph<E, V, H, N>,
    ) -> Self {
        !self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::half_edge::involution::Hedge;
    use crate::half_edge::subgraph::SubSetLike;
    use std::{
        collections::VecDeque,
        iter::{repeat, repeat_n},
    };

    use proptest::{prelude::*, sample::Index as SampleIndex};

    #[test]
    fn base62_roundtrip() {
        let size = 128;
        let mut subset = SuBitGraph::empty(size);
        subset.add(Hedge(0));
        subset.add(Hedge(1));
        subset.add(Hedge(2));
        subset.add(Hedge(31));
        subset.add(Hedge(63));
        subset.add(Hedge(64));
        subset.add(Hedge(127));

        let label = subset.string_label();
        let parsed = <SuBitGraph as SubSetLike>::from_base62(&label, size).unwrap();
        assert_eq!(subset, parsed);

        let empty_label = SuBitGraph::empty(size).string_label();
        let empty_parsed = <SuBitGraph as SubSetLike>::from_base62(&empty_label, size).unwrap();
        assert!(empty_parsed.is_empty());

        let a = "rt12UAag4";
        assert_eq!(SuBitGraph::from_base62(a, 1000).unwrap().n_included(), 27);
        assert_eq!(
            SuBitGraph::from_base62(a, 1000)
                .unwrap()
                .string_label()
                .as_str(),
            a
        )
    }

    /// Lengths around the word boundaries, including the empty set.
    const LENS: [usize; 9] = [0, 1, 63, 64, 65, 127, 128, 129, 200];

    /// A set with the live bits `bits` and the dead bits of its last word set
    /// from `garbage`, as `!` or `truncate` leave them.
    fn subset(bits: &[bool], garbage: usize) -> SuBitGraph {
        let mut set = reference(bits);
        if let Some(last) = set.as_raw_mut_slice().last_mut() {
            *last |= garbage & !tail_mask(bits.len());
        }
        SubSet {
            set,
            id: PhantomData,
        }
    }

    /// A bitvec with zeroed dead bits, on which `bitvec`'s own ops are exact.
    fn reference(bits: &[bool]) -> BitVec {
        let mut set = BitVec::repeat(false, bits.len());
        for (i, &bit) in bits.iter().enumerate() {
            set.set(i, bit);
        }
        set
    }

    fn live(set: &SuBitGraph) -> BitVec {
        set.set.iter().by_vals().collect()
    }

    fn bits() -> impl Strategy<Value = Vec<bool>> {
        prop::sample::select(LENS.as_slice())
            .prop_flat_map(|n| prop::collection::vec(any::<bool>(), n))
    }

    /// Pairs that are independent (any lengths), of equal length, subsets,
    /// equal, one bit apart, or prefix and extension of each other.
    fn pair() -> impl Strategy<Value = (Vec<bool>, Vec<bool>)> {
        (bits(), bits(), 0..6u8, any::<SampleIndex>()).prop_map(|(a, r, relation, index)| {
            let padded = || r.iter().copied().chain(repeat(false));
            let b = match relation {
                0 => r.clone(),
                1 => a.iter().zip(padded()).map(|(_, y)| y).collect(),
                2 => a.iter().zip(padded()).map(|(x, y)| *x && y).collect(),
                3 => a.clone(),
                4 => {
                    let mut b = a.clone();
                    if !b.is_empty() {
                        let i = index.index(b.len());
                        b[i] = !b[i];
                    }
                    b
                }
                _ => a.iter().chain(&r).take(r.len()).copied().collect(),
            };
            (a, b)
        })
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(2000))]

        #[test]
        fn set_algebra_matches_bitvec(
            (a, b) in pair(),
            (ga, gb) in (any::<usize>(), any::<usize>()),
        ) {
            let (sa, sb) = (subset(&a, ga), subset(&b, gb));
            let (ra, rb) = (reference(&a), reference(&b));
            let other = |i: usize| b.get(i).copied().unwrap_or(false);

            prop_assert_eq!(live(&sa.intersection(&sb)), ra.clone() & &rb);
            prop_assert_eq!(live(&sa.union(&sb)), ra.clone() | &rb);
            prop_assert_eq!(live(&sa.sym_diff(&sb)), ra.clone() ^ &rb);
            let difference: BitVec = a.iter().enumerate().map(|(i, &x)| x && !other(i)).collect();
            prop_assert_eq!(live(&sa.subtract(&sb)), difference);
            prop_assert_eq!(live(&!&sa), !ra.clone());
            prop_assert_eq!(live(&!sa.clone()), !ra.clone());

            let common = a.iter().zip(&b).any(|(x, y)| *x && *y);
            let included = a.len() == b.len() && a.iter().zip(&b).all(|(x, y)| *x || !*y);
            prop_assert_eq!(sa.intersects(&sb), common);
            prop_assert_eq!(sa.empty_intersection(&sb), !common);
            prop_assert_eq!(sa.includes(&sb), included);
            prop_assert_eq!(sa.empty_union(&sb), ra.not_any() && rb.not_any());
            prop_assert_eq!(sa.is_empty(), ra.not_any());
            prop_assert_eq!(sa.n_included(), ra.count_ones());
        }

        #[test]
        fn eq_hash_and_ord_match_bitvec(
            (a, b) in pair(),
            (ga, gb, ga2) in (any::<usize>(), any::<usize>(), any::<usize>()),
        ) {
            let (sa, sb) = (subset(&a, ga), subset(&b, gb));
            let (ra, rb) = (reference(&a), reference(&b));

            prop_assert_eq!(sa == sb, ra == rb);
            prop_assert_eq!(sa.cmp(&sb), ra.cmp(&rb));
            prop_assert_eq!(sa.partial_cmp(&sb), Some(ra.cmp(&rb)));

            // The same live bits with other dead bits: equal, equal hash.
            let sa2 = subset(&a, ga2);
            let state = ahash::RandomState::with_seeds(1, 2, 3, 4);
            prop_assert_eq!(&sa, &sa2);
            prop_assert_eq!(sa.cmp(&sa2), Ordering::Equal);
            prop_assert_eq!(state.hash_one(&sa), state.hash_one(&sa2));
            if sa == sb {
                prop_assert_eq!(state.hash_one(&sa), state.hash_one(&sb));
            }
        }

        #[test]
        fn iteration_matches_bitvec(
            a in bits(),
            garbage in any::<usize>(),
            from_back in prop::collection::vec(any::<bool>(), 0..210),
        ) {
            let set = subset(&a, garbage);
            let ones: Vec<usize> = reference(&a).iter_ones().collect();

            prop_assert_eq!(set.included_iter().map(|h| h.0).collect::<Vec<_>>(), ones.clone());
            prop_assert_eq!(
                set.included_iter().rev().map(|h| h.0).collect::<Vec<_>>(),
                ones.iter().rev().copied().collect::<Vec<_>>()
            );

            let mut expected: VecDeque<usize> = ones.into();
            let mut iter = set.included_iter();
            for back in from_back.into_iter().chain(repeat_n(false, expected.len() + 1)) {
                prop_assert_eq!(iter.len(), expected.len());
                prop_assert_eq!(iter.size_hint(), (expected.len(), Some(expected.len())));
                let (got, want) = if back {
                    (iter.next_back(), expected.pop_back())
                } else {
                    (iter.next(), expected.pop_front())
                };
                prop_assert_eq!(got.map(|h| h.0), want);
            }
        }

        #[test]
        fn split_off_halves_are_word_aligned(
            a in bits(),
            garbage in any::<usize>(),
            at in any::<SampleIndex>(),
        ) {
            // `bitvec`'s `split_off` keeps the in-word offset of `at`.
            let at = at.index(a.len() + 1);
            let mut head = subset(&a, garbage);
            let tail = head.split_off(Hedge(at));
            for (half, bits) in [(&head, &a[..at]), (&tail, &a[at..])] {
                let expected = subset(bits, 0);
                prop_assert_eq!(half, &expected);
                prop_assert_eq!(half.n_included(), expected.n_included());
                prop_assert!(half.included_iter().eq(expected.included_iter()));
                prop_assert_eq!(live(&half.union(&expected)), reference(bits));
            }
        }
    }

    #[test]
    fn complement_ignores_dead_bits() {
        // `!` sets the dead bits of the last raw word.
        let full = !SuBitGraph::empty(34);
        assert_ne!(full.set.as_raw_slice()[0], tail_mask(34));
        assert_eq!(full, SuBitGraph::full(34));
        assert_eq!(full.n_included(), 34);
        assert_eq!(full.included_iter().count(), 34);
        assert_eq!(full.included_iter().next_back(), Some(Hedge(33)));
        assert!((!full.clone()).is_empty());
        assert!(full.includes(&SuBitGraph::full(34)));
        let state = ahash::RandomState::with_seeds(1, 2, 3, 4);
        assert_eq!(state.hash_one(&full), state.hash_one(SuBitGraph::full(34)));
    }
}
