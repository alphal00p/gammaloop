use std::{
    hash::Hash,
    marker::PhantomData,
    ops::{
        BitAndAssign, BitOrAssign, BitXorAssign, Index, Not, Range, RangeFrom, RangeInclusive,
        RangeTo, RangeToInclusive,
    },
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
        BaseSubgraph, Inclusion, ModifySubSet, SubGraphLike, SubGraphOps, SubSetIter, SubSetLike,
        SubSetOps, BASE62_ALPHABET,
    },
    swap::Swap,
    typed_vec::IndexLike,
    HedgeGraph,
};

pub type SuBitGraph = SubSet<Hedge>;

#[cfg(feature = "rkyv")]
#[derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)]
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
#[derive(Clone, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct SubSet<ID> {
    #[cfg_attr(feature = "bincode", bincode(with_serde))]
    #[cfg_attr(feature = "rkyv", with(BitVecRkyv))]
    set: BitVec,
    id: std::marker::PhantomData<ID>,
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
        Self {
            set: self.set.split_off(at.into()),
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
        &self.intersection(other) == other
    }

    fn intersects(&self, other: &SubSet<ID>) -> bool {
        !(self.intersection(other).is_empty())
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
        SubSetIter {
            iter: self.set.iter_ones().map(ID::from),
            // len: self.n_included(),
        }
    }

    fn n_included(&self) -> usize {
        self.set.count_ones()
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
        self.n_included() == 0
    }
}

impl<ID: IndexLike> SubSetOps<ID> for SubSet<ID> {
    fn intersect_with(&mut self, other: &Self) {
        self.set.bitand_assign(&other.set)
    }

    fn union_with(&mut self, other: &Self) {
        self.set.bitor_assign(&other.set)
    }

    fn sym_diff_with(&mut self, other: &Self) {
        self.set.bitxor_assign(&other.set)
    }

    fn empty_union(&self, other: &Self) -> bool {
        self.union(other).is_empty()
    }

    fn empty_intersection(&self, other: &Self) -> bool {
        self.intersection(other).is_empty()
    }

    fn subtract_with(&mut self, other: &Self) {
        self.set.bitand_assign(!other.set.clone());
    }
}

impl SubGraphOps for SuBitGraph {
    fn complement<E, V, H, N: NodeStorage<NodeData = V>>(
        &self,
        _graph: &HedgeGraph<E, V, H, N>,
    ) -> Self {
        Self {
            set: !self.set.clone(),
            id: PhantomData,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::half_edge::involution::Hedge;
    use crate::half_edge::subgraph::SubSetLike;

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
}
