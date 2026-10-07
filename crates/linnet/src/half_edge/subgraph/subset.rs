use std::{
    cmp::Ordering,
    fmt,
    hash::{Hash, Hasher},
    marker::PhantomData,
    ops::{Index, Not, Range, RangeFrom, RangeInclusive, RangeTo, RangeToInclusive},
};

#[cfg(feature = "rkyv")]
use rkyv::{
    bytecheck::CheckBytes,
    validation::ArchiveContext,
    vec::{ArchivedVec, VecResolver},
    Archive, Archived, Fallible,
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

/// Words stored inline: sets of up to `64 * INLINE_WORDS` IDs never allocate.
const INLINE_WORDS: usize = 2;

const WORD_BITS: usize = u64::BITS as usize;

/// Mask of the live bits in the last word of a `len`-bit set.
fn tail_mask(len: usize) -> u64 {
    u64::MAX >> ((WORD_BITS - len % WORD_BITS) % WORD_BITS)
}

fn n_words(len: usize) -> usize {
    len.div_ceil(WORD_BITS)
}

/// Whether `op(a, b)` is non-zero for some pair of words of `a` and `b`
/// (zipped to the shorter), OR-reducing blocks of words so the scan vectorises.
#[inline]
fn any_word(a: &[u64], b: &[u64], op: impl Fn(u64, u64) -> u64) -> bool {
    const BLOCK: usize = 8;
    let n = a.len().min(b.len());
    let (a, b) = (&a[..n], &b[..n]);
    if n < BLOCK {
        return a.iter().zip(b).any(|(&x, &y)| op(x, y) != 0);
    }
    let ((a_blocks, a_rest), (b_blocks, b_rest)) = (a.as_chunks::<BLOCK>(), b.as_chunks::<BLOCK>());
    a_blocks
        .iter()
        .zip(b_blocks)
        .any(|(x, y)| (0..BLOCK).fold(0, |acc, i| acc | op(x[i], y[i])) != 0)
        || a_rest.iter().zip(b_rest).any(|(&x, &y)| op(x, y) != 0)
}

/// Why stored words do not form a set of the stored length: a `len`-ID set
/// takes exactly `len.div_ceil(bits)` words of `bits` bits, with the bits past
/// `len` zero.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum SubSetWordsError {
    #[error("a set of {len} IDs takes {expected} words, not {found}")]
    WordCount {
        len: usize,
        expected: usize,
        found: usize,
    },
    #[error("a set of {len} IDs has bits set past its last ID")]
    DeadBits { len: usize },
}

impl SubSetWordsError {
    /// Checks that `words`, with ID `w * bits + i` at bit `i` of word `w` for
    /// `bits` the width of `W`, hold exactly a set of `len` IDs.
    fn check<W: Copy + Into<u64>>(len: usize, words: &[W]) -> Result<(), Self> {
        let bits = 8 * size_of::<W>();
        let expected = len.div_ceil(bits);
        if words.len() != expected {
            return Err(Self::WordCount {
                len,
                expected,
                found: words.len(),
            });
        }
        let live = len - expected.saturating_sub(1) * bits;
        match words.last() {
            Some(&last) if last.into().checked_shr(live as u32).unwrap_or(0) != 0 => {
                Err(Self::DeadBits { len })
            }
            _ => Ok(()),
        }
    }
}

/// The words of a set: in `inline` (unused words zero) when they fit in
/// `INLINE_WORDS`, otherwise exactly `n_words(len)` words in `heap` (and
/// `inline` all zero). `repr(C)` keeps `inline` at offset 0, so a whole-array
/// store of it forwards to the 16-byte loads that copy the set.
#[derive(Clone)]
#[repr(C)]
struct Words {
    inline: [u64; INLINE_WORDS],
    heap: Option<Box<[u64]>>,
}

impl Words {
    fn zeroed(n: usize) -> Self {
        Words {
            inline: [0; INLINE_WORDS],
            heap: (n > INLINE_WORDS).then(|| vec![0; n].into_boxed_slice()),
        }
    }

    /// `words` (the exact word count of the set) in canonical form.
    fn from_vec(words: Vec<u64>) -> Self {
        let mut inline = [0; INLINE_WORDS];
        if words.len() <= INLINE_WORDS {
            inline[..words.len()].copy_from_slice(&words);
            Words { inline, heap: None }
        } else {
            Words {
                inline,
                heap: Some(words.into_boxed_slice()),
            }
        }
    }
}

/// A set of IDs `0..size` stored as a bitset of 64-bit words.
///
/// Sets of up to 128 IDs live inline and never allocate. The bits past `size`
/// are always zero, so equality, hashing and the set algebra work on whole
/// words. The order is lexicographic from ID 0 (the first ID in only one set
/// decides, and the set that has it is greater), with a strict prefix before
/// its extensions.
///
/// Serde writes `{len, words}` with 32-bit words, ID `32 * w + i` at bit `i` of
/// word `w`, so every word survives decoders whose integers are doubles
/// (JavaScript) or signed 64-bit (Typst). Decoding rejects a word count other
/// than `len.div_ceil(32)` and bits set past `len`; archive validation checks
/// the same for the archived 64-bit words.
#[repr(C)]
pub struct SubSet<ID> {
    words: Words,
    len: usize,
    id: PhantomData<ID>,
}

impl<ID> fmt::Debug for SubSet<ID> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SubSet")
            .field("len", &self.len)
            .field("words", &self.words())
            .finish()
    }
}

/// The serialized 32-bit words of a set.
#[cfg(feature = "serde")]
struct SerialWords<'a, ID>(&'a SubSet<ID>);

#[cfg(feature = "serde")]
impl<ID> serde::Serialize for SerialWords<'_, ID> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let words = self.0.words();
        serializer.collect_seq(
            (0..self.0.len.div_ceil(32)).map(|i| (words[i / 2] >> (32 * (i % 2))) as u32),
        )
    }
}

#[cfg(feature = "serde")]
impl<ID> serde::Serialize for SubSet<ID> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;
        let mut state = serializer.serialize_struct("SubSet", 2)?;
        state.serialize_field("len", &self.len)?;
        state.serialize_field("words", &SerialWords(self))?;
        state.end()
    }
}

#[cfg(feature = "serde")]
#[derive(serde::Deserialize)]
#[serde(rename = "SubSet", deny_unknown_fields)]
struct SerialSubSet {
    len: usize,
    words: Vec<u32>,
}

#[cfg(feature = "serde")]
impl<'de, ID> serde::Deserialize<'de> for SubSet<ID> {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let SerialSubSet { len, words } = SerialSubSet::deserialize(deserializer)?;
        SubSetWordsError::check(len, &words).map_err(serde::de::Error::custom)?;
        let words = words
            .chunks(2)
            .map(|pair| {
                let high = pair.get(1).map_or(0, |&high| u64::from(high) << 32);
                u64::from(pair[0]) | high
            })
            .collect();
        Ok(Self::from_words(len, words))
    }
}

#[cfg(feature = "bincode")]
impl<ID, Context> bincode::Decode<Context> for SubSet<ID> {
    fn decode<D: bincode::de::Decoder<Context = Context>>(
        decoder: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        <bincode::serde::Compat<Self> as bincode::Decode<Context>>::decode(decoder)
            .map(|compat| compat.0)
    }
}

#[cfg(feature = "bincode")]
impl<'de, ID, Context> bincode::BorrowDecode<'de, Context> for SubSet<ID> {
    fn borrow_decode<D: bincode::de::BorrowDecoder<'de, Context = Context>>(
        decoder: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        <bincode::serde::Compat<Self> as bincode::Decode<Context>>::decode(decoder)
            .map(|compat| compat.0)
    }
}

#[cfg(feature = "bincode")]
impl<ID> bincode::Encode for SubSet<ID> {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        bincode::serde::BorrowCompat(self).encode(encoder)
    }
}

/// The archived form of a [`SubSet`]: its length and exactly
/// `len.div_ceil(64)` words with the bits past `len` zero, which validation
/// checks.
#[cfg(feature = "rkyv")]
#[repr(C)]
pub struct ArchivedSubSet<ID> {
    words: ArchivedVec<u64>,
    len: Archived<usize>,
    id: PhantomData<ID>,
}

#[cfg(feature = "rkyv")]
/// The resolver of an archived [`SubSet`]: where its words were written.
pub struct SubSetResolver {
    words: VecResolver,
}

#[cfg(feature = "rkyv")]
impl<ID> Archive for SubSet<ID> {
    type Archived = ArchivedSubSet<ID>;
    type Resolver = SubSetResolver;

    unsafe fn resolve(&self, pos: usize, resolver: Self::Resolver, out: *mut Self::Archived) {
        let (offset, words) = rkyv::out_field!(out.words);
        ArchivedVec::resolve_from_slice(self.words(), pos + offset, resolver.words, words);
        let (offset, len) = rkyv::out_field!(out.len);
        self.len.resolve(pos + offset, (), len);
    }
}

#[cfg(feature = "rkyv")]
impl<ID, S: rkyv::ser::ScratchSpace + rkyv::ser::Serializer + ?Sized> rkyv::Serialize<S>
    for SubSet<ID>
{
    fn serialize(&self, serializer: &mut S) -> Result<Self::Resolver, S::Error> {
        Ok(SubSetResolver {
            words: ArchivedVec::serialize_from_slice(self.words(), serializer)?,
        })
    }
}

#[cfg(feature = "rkyv")]
impl<ID, D: Fallible + ?Sized> rkyv::Deserialize<SubSet<ID>, D> for ArchivedSubSet<ID> {
    /// Rebuilds the canonical form, so an archive that skipped validation
    /// still yields a well-formed set.
    fn deserialize(&self, _: &mut D) -> Result<SubSet<ID>, D::Error> {
        Ok(SubSet::from_words(self.size(), self.words.to_vec()))
    }
}

/// Why an archived [`SubSet`] is invalid.
#[cfg(feature = "rkyv")]
#[derive(Debug, thiserror::Error)]
pub enum ArchivedSubSetError<E> {
    #[error("invalid archived subset words: {0}")]
    Words(E),
    #[error(transparent)]
    Invalid(#[from] SubSetWordsError),
}

#[cfg(feature = "rkyv")]
impl<ID, C: ArchiveContext + ?Sized> CheckBytes<C> for ArchivedSubSet<ID>
where
    C::Error: rkyv::bytecheck::Error,
{
    type Error = ArchivedSubSetError<<ArchivedVec<u64> as CheckBytes<C>>::Error>;

    unsafe fn check_bytes<'a>(
        value: *const Self,
        context: &mut C,
    ) -> Result<&'a Self, Self::Error> {
        let words = ArchivedVec::<u64>::check_bytes(std::ptr::addr_of!((*value).words), context)
            .map_err(ArchivedSubSetError::Words)?;
        let len = Archived::<usize>::check_bytes(std::ptr::addr_of!((*value).len), context)
            .unwrap_or_else(|never| match never {});
        SubSetWordsError::check(*len as usize, words.as_slice())?;
        Ok(&*value)
    }
}

#[cfg(feature = "rkyv")]
impl<ID> ArchivedSubSet<ID> {
    pub fn size(&self) -> usize {
        self.len as usize
    }
}

#[cfg(feature = "rkyv")]
impl<ID: IndexLike> ArchivedSubSet<ID> {
    pub fn contains(&self, index: ID) -> bool {
        self.contains_index(index.into())
    }

    pub fn contains_index(&self, index: usize) -> bool {
        index < self.size()
            && self
                .words
                .as_slice()
                .get(index / WORD_BITS)
                .is_some_and(|word| (*word >> (index % WORD_BITS)) & 1 == 1)
    }

    pub fn intersects_owned(&self, other: &SubSet<ID>) -> bool {
        any_word(self.words.as_slice(), other.words(), |a, b| a & b)
    }
}

impl<ID> Clone for SubSet<ID> {
    fn clone(&self) -> Self {
        Self {
            words: self.words.clone(),
            len: self.len,
            id: PhantomData,
        }
    }
}

impl<ID> PartialEq for SubSet<ID> {
    fn eq(&self, other: &Self) -> bool {
        self.len == other.len
            && self.words.inline == other.words.inline
            && self.words.heap == other.words.heap
    }
}

impl<ID> Eq for SubSet<ID> {}

impl<ID> Hash for SubSet<ID> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_usize(self.len);
        for &word in self.words() {
            state.write_u64(word);
        }
    }
}

impl<ID> PartialOrd for SubSet<ID> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<ID> Ord for SubSet<ID> {
    fn cmp(&self, other: &Self) -> Ordering {
        // The zip ends with the last word of the shorter set, whose dead bits
        // are zero: an ID of the longer set past the shorter length decides
        // there as the greater length would.
        for (&a, &b) in self.words().iter().zip(other.words()) {
            let diff = a ^ b;
            if diff != 0 {
                let first = diff.isolate_lowest_one();
                return if a & first != 0 {
                    Ordering::Greater
                } else {
                    Ordering::Less
                };
            }
        }
        self.len.cmp(&other.len)
    }
}

impl<ID> Not for &SubSet<ID> {
    type Output = SubSet<ID>;
    fn not(self) -> Self::Output {
        self.clone().not()
    }
}

impl<ID> FromIterator<bool> for SubSet<ID> {
    fn from_iter<T: IntoIterator<Item = bool>>(iter: T) -> Self {
        let (mut words, mut len) = (Vec::new(), 0);
        for bit in iter {
            if len % WORD_BITS == 0 {
                words.push(0);
            }
            *words.last_mut().unwrap() |= u64::from(bit) << (len % WORD_BITS);
            len += 1;
        }
        Self::from_words(len, words)
    }
}

impl<ID: IndexLike> Swap<ID> for SubSet<ID> {
    fn is_zero_length(&self) -> bool {
        self.len == 0
    }

    fn len(&self) -> ID {
        self.len.into()
    }

    fn swap(&mut self, i: ID, j: ID) {
        let (i, j) = (i.into(), j.into());
        let (a, b) = (self.get(i), self.get(j));
        if a != b {
            self.set_bit(i, b);
            self.set_bit(j, a);
        }
    }
}

impl<ID> SubSet<ID> {
    fn with_len(len: usize) -> Self {
        Self {
            len,
            words: Words::zeroed(n_words(len)),
            id: PhantomData,
        }
    }

    /// The set of `len` IDs whose IDs are the set bits of `words` below `len`,
    /// with missing words zero.
    fn from_words(len: usize, mut words: Vec<u64>) -> Self {
        words.resize(n_words(len), 0);
        if let Some(last) = words.last_mut() {
            *last &= tail_mask(len);
        }
        Self {
            len,
            words: Words::from_vec(words),
            id: PhantomData,
        }
    }

    /// The words: bit `i` of word `w` is ID `w * 64 + i`; dead bits are zero.
    #[inline]
    pub(crate) fn words(&self) -> &[u64] {
        match &self.words.heap {
            None => &self.words.inline[..n_words(self.len)],
            Some(words) => words,
        }
    }

    #[inline]
    fn words_mut(&mut self) -> &mut [u64] {
        let n = n_words(self.len);
        match &mut self.words.heap {
            None => &mut self.words.inline[..n],
            Some(words) => words,
        }
    }

    #[inline]
    fn get(&self, index: usize) -> bool {
        assert!(
            index < self.len,
            "index {index} out of range for SubSet of size {}",
            self.len
        );
        let word = match &self.words.heap {
            None => self.words.inline[index / WORD_BITS],
            Some(words) => words[index / WORD_BITS],
        };
        word >> (index % WORD_BITS) & 1 != 0
    }

    #[inline]
    fn set_bit(&mut self, index: usize, value: bool) {
        assert!(
            index < self.len,
            "index {index} out of range for SubSet of size {}",
            self.len
        );
        let bit = 1 << (index % WORD_BITS);
        let word = match &mut self.words.heap {
            None => &mut self.words.inline[index / WORD_BITS],
            Some(words) => &mut words[index / WORD_BITS],
        };
        if value {
            *word |= bit;
        } else {
            *word &= !bit;
        }
    }

    /// Clears the bits past `len`, including unused inline words.
    #[inline]
    fn clear_tail(&mut self) {
        let (n, mask) = (n_words(self.len), tail_mask(self.len));
        match &mut self.words.heap {
            None => {
                let words = &mut self.words.inline;
                if let Some(last) = n.checked_sub(1) {
                    words[last] &= mask;
                }
                words[n..].fill(0);
            }
            Some(words) => {
                if let Some(last) = words.last_mut() {
                    *last &= mask;
                }
            }
        }
    }

    pub fn is_full(&self) -> bool {
        match self.words().split_last() {
            None => true,
            Some((&last, body)) => {
                body.iter().all(|&w| w == u64::MAX) && last == tail_mask(self.len)
            }
        }
    }

    pub fn iter(&self) -> impl DoubleEndedIterator<Item = bool> + '_ {
        (0..self.len).map(|i| self.get(i))
    }

    pub fn split_off(&mut self, at: ID) -> Self
    where
        ID: IndexLike,
    {
        let at: usize = at.into();
        assert!(
            at <= self.len,
            "split_off index {at} past the end {}",
            self.len
        );
        let (first, shift) = (at / WORD_BITS, at % WORD_BITS);
        let words = self.words();
        let tail = (first..words.len())
            .map(|i| {
                let high = match (shift, words.get(i + 1)) {
                    (0, _) | (_, None) => 0,
                    (_, Some(&next)) => next << (WORD_BITS - shift),
                };
                words[i] >> shift | high
            })
            .collect();
        let tail = Self::from_words(self.len - at, tail);
        *self = Self::from_words(at, words[..n_words(at)].to_vec());
        tail
    }

    pub fn full(size: usize) -> Self {
        let mut set = Self::with_len(size);
        set.words_mut().fill(u64::MAX);
        set.clear_tail();
        set
    }

    pub fn push(&mut self, val: bool) {
        if self.len.is_multiple_of(WORD_BITS) && n_words(self.len) >= INLINE_WORDS {
            let mut words = self.words().to_vec();
            words.push(0);
            self.words = Words::from_vec(words);
        }
        self.len += 1;
        if val {
            self.set_bit(self.len - 1, true);
        }
    }

    /// The set of size `min(size, 64)` whose IDs are the set bits of `num`.
    pub fn from_usize(num: usize, size: usize) -> Self {
        let mut set = Self::with_len(size.min(usize::BITS as usize));
        if let Some(word) = set.words_mut().first_mut() {
            *word = num as u64;
        }
        set.clear_tail();
        set
    }

    pub fn pop(&mut self) -> Option<bool> {
        let last = self.len.checked_sub(1)?;
        let bit = self.get(last);
        if self.words.heap.is_some() && last.is_multiple_of(WORD_BITS) {
            *self = Self::from_words(last, self.words().to_vec());
        } else {
            self.set_bit(last, false);
            self.len = last;
        }
        Some(bit)
    }

    pub fn clear(&mut self) {
        self.words_mut().fill(0);
    }

    /// Replaces each word `a` of `self` by `op(a, b)` with `b` the matching word
    /// of `other`, zero-extended or truncated to the length of `self`.
    /// `op(0, 0)` must be 0.
    #[inline]
    fn zip_words_with(&mut self, other: &Self, op: impl Fn(u64, u64) -> u64) {
        match (&mut self.words.heap, &other.words.heap) {
            // Whole-array op on the inline buffers (unused words are zero on
            // both sides): straight-line SIMD and one store per 16 bytes.
            (None, None) => {
                let (a, b) = (&mut self.words.inline, &other.words.inline);
                *a = std::array::from_fn(|i| op(a[i], b[i]));
            }
            _ if self.len == other.len => {
                for (a, &b) in self.words_mut().iter_mut().zip(other.words()) {
                    *a = op(*a, b);
                }
            }
            _ => {
                let mut other = other.words().iter().copied();
                for word in self.words_mut() {
                    *word = op(*word, other.next().unwrap_or(0));
                }
            }
        }
        if self.len != other.len {
            self.clear_tail();
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
    front: u64,
    front_base: usize,
    /// Untouched words strictly between the front and the back word.
    middle: &'a [u64],
    /// Remaining bits of the back word, whose bit 0 is ID `back_base`.
    back: u64,
    back_base: usize,
    id: PhantomData<fn() -> ID>,
}

impl<'a, ID> SubSetIter<'a, ID> {
    fn new(set: &'a SubSet<ID>) -> Self {
        let words = set.words();
        let (front, middle, back) = match words {
            [] => (0, &[][..], 0),
            [only] => (0, &[][..], *only),
            [first, middle @ .., last] => (*first, middle, *last),
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
    #[inline]
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

impl<ID> Not for SubSet<ID> {
    type Output = Self;
    fn not(mut self) -> Self::Output {
        for word in self.words_mut().iter_mut() {
            *word = !*word;
        }
        self.clear_tail();
        self
    }
}

static TRUE: bool = true;
static FALSE: bool = false;

impl<ID: Into<usize>> Index<ID> for SubSet<ID> {
    type Output = bool;

    fn index(&self, index: ID) -> &Self::Output {
        if self.get(index.into()) {
            &TRUE
        } else {
            &FALSE
        }
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
        self.len == other.len && !any_word(self.words(), other.words(), |a, b| b & !a)
    }

    fn intersects(&self, other: &SubSet<ID>) -> bool {
        any_word(self.words(), other.words(), |a, b| a & b)
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
        self.set_bit(hedge.into(), true);
    }

    fn sub(&mut self, hedge: ID) {
        self.set_bit(hedge.into(), false);
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
        self.len
    }

    fn has_greater(&self, hedge: ID) -> bool {
        let from: usize = hedge.into();
        if from >= self.len {
            return false;
        }
        let first = from / WORD_BITS;
        let words = self.words();
        words[first] >> (from % WORD_BITS) != 0 || words[first + 1..].iter().any(|&w| w != 0)
    }

    fn join_mut(&mut self, other: Self) {
        let shift = self.len % WORD_BITS;
        let mut words = self.words().to_vec();
        if shift == 0 {
            words.extend_from_slice(other.words());
        } else {
            let mut carry_in = other.words().iter().copied();
            if let Some(first) = carry_in.next() {
                *words.last_mut().unwrap() |= first << shift;
                let mut prev = first;
                for word in carry_in {
                    words.push((prev >> (WORD_BITS - shift)) | (word << shift));
                    prev = word;
                }
                words.push(prev >> (WORD_BITS - shift));
            }
        }
        *self = Self::from_words(self.len + other.len, words);
    }

    fn included_iter(&self) -> Self::BaseIter<'_> {
        SubSetIter::new(self)
    }

    fn n_included(&self) -> usize {
        self.words()
            .iter()
            .map(|word| word.count_ones() as usize)
            .sum()
    }

    fn empty(size: usize) -> Self {
        Self::with_len(size)
    }

    fn string_label(&self) -> String {
        if self.is_empty() {
            return "0".to_string();
        }

        let mut digits = vec![0u8]; // Initialize with a single zero digit

        // Iterate over the bits from MSB to LSB
        for bit in self.iter().rev() {
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

        let mut set = Self::with_len(size);
        for (index, bit) in bits.into_iter().enumerate() {
            if bit {
                set.set_bit(index, true);
            }
        }

        Some(set)
    }

    fn is_empty(&self) -> bool {
        self.words().iter().all(|&word| word == 0)
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
    use bitvec::vec::BitVec;
    use std::{
        collections::VecDeque,
        hash::BuildHasher,
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

    /// Lengths around the word, inline-capacity and 8-word `any_word` block
    /// boundaries (sets of more than 448 IDs take the block path).
    const LENS: [usize; 18] = [
        0, 1, 63, 64, 65, 127, 128, 129, 191, 200, 256, 257, 300, 448, 449, 512, 600, 1000,
    ];

    fn subset(bits: &[bool]) -> SuBitGraph {
        bits.iter().copied().collect()
    }

    /// The `bitvec` reference on which its own ops define the expected result;
    /// its order is the derived order of the former `bitvec`-backed set.
    fn reference(bits: &[bool]) -> BitVec {
        bits.iter().copied().collect()
    }

    /// The bits of `set`, after checking its canonical form: inline exactly
    /// when it fits, the exact word count and every unused bit zero.
    fn live(set: &SuBitGraph) -> BitVec {
        match &set.words.heap {
            None => {
                let words = &set.words.inline;
                assert!(n_words(set.len) <= INLINE_WORDS, "large set inline");
                assert!(
                    words[n_words(set.len)..].iter().all(|&w| w == 0),
                    "unused inline words set"
                );
            }
            Some(words) => {
                assert!(words.len() > INLINE_WORDS, "small set on the heap");
                assert_eq!(set.words.inline, [0; INLINE_WORDS], "inline words set");
            }
        }
        assert_eq!(set.words().len(), n_words(set.len));
        if let Some(last) = set.words().last() {
            assert_eq!(last & !tail_mask(set.len), 0, "dead bits set");
        }
        set.iter().collect()
    }

    fn bits() -> impl Strategy<Value = Vec<bool>> {
        prop::sample::select(LENS.as_slice())
            .prop_flat_map(|n| prop::collection::vec(any::<bool>(), n))
    }

    /// Pairs that are independent (any lengths), of equal length, subsets,
    /// equal, one bit apart (inclusion hinges on one word), disjoint but for at
    /// most one ID (overlap hinges on one word), or prefix and extension of
    /// each other.
    fn pair() -> impl Strategy<Value = (Vec<bool>, Vec<bool>)> {
        (bits(), bits(), 0..7u8, any::<SampleIndex>()).prop_map(|(a, r, relation, index)| {
            let padded = || r.iter().copied().chain(repeat(false));
            let flip_one = |mut b: Vec<bool>| {
                if !b.is_empty() {
                    let i = index.index(b.len());
                    b[i] = !b[i];
                }
                b
            };
            let b = match relation {
                0 => r.clone(),
                1 => a.iter().zip(padded()).map(|(_, y)| y).collect(),
                2 => a.iter().zip(padded()).map(|(x, y)| *x && y).collect(),
                3 => a.clone(),
                4 => flip_one(a.clone()),
                5 => flip_one(a.iter().map(|x| !x).collect()),
                _ => a.iter().chain(&r).take(r.len()).copied().collect(),
            };
            (a, b)
        })
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(2000))]

        #[test]
        fn set_algebra_matches_bitvec((a, b) in pair()) {
            let (sa, sb) = (subset(&a), subset(&b));
            let (ra, rb) = (reference(&a), reference(&b));
            let other = |i: usize| b.get(i).copied().unwrap_or(false);
            let and: BitVec = a.iter().enumerate().map(|(i, &x)| x && other(i)).collect();
            let or: BitVec = a.iter().enumerate().map(|(i, &x)| x || other(i)).collect();
            let xor: BitVec = a.iter().enumerate().map(|(i, &x)| x ^ other(i)).collect();
            let difference: BitVec = a.iter().enumerate().map(|(i, &x)| x && !other(i)).collect();

            prop_assert_eq!(live(&sa.intersection(&sb)), and.clone());
            prop_assert_eq!(live(&sa.union(&sb)), or.clone());
            prop_assert_eq!(live(&sa.sym_diff(&sb)), xor.clone());
            prop_assert_eq!(live(&sa.subtract(&sb)), difference.clone());
            for (op, expected) in [(0, and), (1, or), (2, xor), (3, difference)] {
                let mut s = sa.clone();
                match op {
                    0 => s.intersect_with(&sb),
                    1 => s.union_with(&sb),
                    2 => s.sym_diff_with(&sb),
                    _ => s.subtract_with(&sb),
                }
                prop_assert_eq!(live(&s), expected);
            }
            prop_assert_eq!(live(&!&sa), !ra.clone());
            prop_assert_eq!(live(&!sa.clone()), !ra.clone());

            let common = a.iter().zip(&b).any(|(x, y)| *x && *y);
            let included = a.len() == b.len() && a.iter().zip(&b).all(|(x, y)| *x || !*y);
            prop_assert_eq!(sa.intersects(&sb), common);
            prop_assert_eq!(sa.empty_intersection(&sb), !common);
            prop_assert_eq!(sa.includes(&sb), included);
            prop_assert_eq!(sa.empty_union(&sb), ra.not_any() && rb.not_any());
            prop_assert_eq!(sa.is_empty(), ra.not_any());
            prop_assert_eq!(sa.is_full(), ra.all());
            prop_assert_eq!(sa.n_included(), ra.count_ones());

            let joined = sa.clone().join(sb.clone());
            let mut rjoined = ra.clone();
            rjoined.extend_from_bitslice(&rb);
            prop_assert_eq!(live(&joined), rjoined);
        }

        #[test]
        fn eq_hash_and_ord_match_bitvec((a, b) in pair()) {
            let (sa, sb) = (subset(&a), subset(&b));
            let (ra, rb) = (reference(&a), reference(&b));

            prop_assert_eq!(sa == sb, ra == rb);
            prop_assert_eq!(sa.cmp(&sb), ra.cmp(&rb));
            prop_assert_eq!(sa.partial_cmp(&sb), Some(ra.cmp(&rb)));
            // The hash is the length and then the live words, as the former
            // `bitvec`-backed set wrote them (dead bits zero).
            let state = ahash::RandomState::with_seeds(1, 2, 3, 4);
            let mut hasher = state.build_hasher();
            hasher.write_usize(a.len());
            for &word in ra.as_raw_slice() {
                hasher.write_usize(word);
            }
            prop_assert_eq!(state.hash_one(&sa), hasher.finish());
            if sa == sb {
                prop_assert_eq!(state.hash_one(&sa), state.hash_one(&sb));
            }
        }

        #[test]
        fn constructions_agree(a in bits(), at in any::<SampleIndex>()) {
            let set = subset(&a);
            prop_assert_eq!(live(&set), reference(&a));

            let mut pushed = SuBitGraph::empty(0);
            let mut added = SuBitGraph::empty(a.len());
            for (i, &bit) in a.iter().enumerate() {
                pushed.push(bit);
                if bit {
                    added.add(Hedge(i));
                }
            }
            let ones = a.iter().enumerate().filter(|(_, &x)| x).map(|(i, _)| Hedge(i));
            let from_hedges = SuBitGraph::from_hedge_iter(ones, a.len());
            for other in [&pushed, &added, &from_hedges] {
                prop_assert_eq!(live(other), reference(&a));
                prop_assert_eq!(other, &set);
            }
            prop_assert_eq!(
                SuBitGraph::from_base62(&set.string_label(), a.len()),
                Some(set.clone())
            );
            prop_assert_eq!(live(&SuBitGraph::full(a.len())), reference(&vec![true; a.len()]));
            prop_assert_eq!(SuBitGraph::full(a.len()), !SuBitGraph::empty(a.len()));

            let mut cleared = set.clone();
            cleared.clear();
            prop_assert_eq!(live(&cleared), reference(&vec![false; a.len()]));

            let num = a.iter().take(64).enumerate().fold(0, |acc, (i, &x)| acc | usize::from(x) << i);
            let size = at.index(a.len() + 1);
            let low: Vec<bool> = a.iter().take(size.min(64)).copied().collect();
            prop_assert_eq!(live(&SuBitGraph::from_usize(num, size)), reference(&low));
        }

        #[test]
        fn single_id_ops_match_bitvec(
            a in bits(),
            ops in prop::collection::vec((0..3u8, any::<SampleIndex>(), any::<SampleIndex>()), 0..20),
        ) {
            prop_assume!(!a.is_empty());
            let mut set = subset(&a);
            let mut expected = reference(&a);
            for (op, i, j) in ops {
                let (i, j) = (i.index(a.len()), j.index(a.len()));
                match op {
                    0 => {
                        set.add(Hedge(i));
                        expected.set(i, true);
                    }
                    1 => {
                        set.sub(Hedge(i));
                        expected.set(i, false);
                    }
                    _ => {
                        set.swap(Hedge(i), Hedge(j));
                        expected.swap(i, j);
                    }
                }
                prop_assert_eq!(live(&set), expected.clone());
            }
            for i in 0..a.len() {
                prop_assert_eq!(set[Hedge(i)], expected[i]);
                prop_assert_eq!(set.includes(&Hedge(i)), expected[i]);
            }
        }

        #[test]
        fn iteration_matches_bitvec(
            a in bits(),
            from_back in prop::collection::vec(any::<bool>(), 0..310),
        ) {
            let set = subset(&a);
            let ones: Vec<usize> = reference(&a).iter_ones().collect();

            prop_assert_eq!(set.included_iter().map(|h| h.0).collect::<Vec<_>>(), ones.clone());
            prop_assert_eq!(
                set.included_iter().rev().map(|h| h.0).collect::<Vec<_>>(),
                ones.iter().rev().copied().collect::<Vec<_>>()
            );
            prop_assert!(set.iter().rev().eq(a.iter().rev().copied()));

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
        fn split_off_push_pop_and_has_greater(
            a in bits(),
            at in any::<SampleIndex>(),
        ) {
            let at = at.index(a.len() + 1);
            let mut head = subset(&a);
            let tail = head.split_off(Hedge(at));
            for (half, bits) in [(&head, &a[..at]), (&tail, &a[at..])] {
                prop_assert_eq!(live(half), reference(bits));
                prop_assert_eq!(half, &subset(bits));
            }
            let set = subset(&a);
            for from in [0, at, a.len()] {
                let expected = a.iter().skip(from).any(|&x| x);
                prop_assert_eq!(set.has_greater(Hedge(from)), expected);
            }
            let mut popped = set.clone();
            let mut rpopped = a.clone();
            while let Some(bit) = popped.pop() {
                prop_assert_eq!(Some(bit), rpopped.pop());
                // `==` compares the whole storage, so equality with the
                // canonical `subset` also checks the form (and stays fast for
                // the quadratic number of bits popped from long sets).
                prop_assert_eq!(&popped, &subset(&rpopped));
            }
            prop_assert_eq!(rpopped.pop(), None);
        }
    }

    #[cfg(feature = "serde")]
    proptest! {
        #[test]
        fn serde_roundtrips_with_32_bit_words(a in bits()) {
            let set = subset(&a);
            let json = serde_json::to_value(&set).unwrap();
            prop_assert_eq!(&json["len"], &serde_json::json!(a.len()));
            let words = json["words"].as_array().unwrap();
            prop_assert_eq!(words.len(), a.len().div_ceil(32));
            for (w, word) in words.iter().enumerate() {
                let word = word.as_u64().unwrap();
                prop_assert!(word <= u64::from(u32::MAX));
                for i in 0..32 {
                    let id = 32 * w + i;
                    prop_assert_eq!(word >> i & 1 == 1, a.get(id).copied().unwrap_or(false));
                }
            }
            let back: SuBitGraph = serde_json::from_value(json).unwrap();
            prop_assert_eq!(live(&back), reference(&a));
            prop_assert_eq!(back, set);
        }
    }

    #[cfg(feature = "bincode")]
    proptest! {
        #[test]
        fn bincode_roundtrips(a in bits()) {
            let set = subset(&a);
            let config = bincode::config::standard();
            let bytes = bincode::encode_to_vec(&set, config).unwrap();
            let (back, read): (SuBitGraph, usize) = bincode::decode_from_slice(&bytes, config).unwrap();
            prop_assert_eq!(read, bytes.len());
            prop_assert_eq!(live(&back), reference(&a));
        }
    }

    #[cfg(feature = "rkyv")]
    proptest! {
        #[test]
        fn rkyv_roundtrips((a, b) in pair()) {
            let (set, other) = (subset(&a), subset(&b));
            let bytes = rkyv::to_bytes::<_, 256>(&set).unwrap();
            let archived = rkyv::check_archived_root::<SuBitGraph>(&bytes).unwrap();
            prop_assert_eq!(archived.words.len(), n_words(a.len()));
            prop_assert_eq!(archived.size(), a.len());
            for i in 0..a.len() + 2 {
                prop_assert_eq!(archived.contains_index(i), a.get(i).copied().unwrap_or(false));
            }
            prop_assert_eq!(archived.intersects_owned(&other), set.intersects(&other));
            let back: SuBitGraph =
                rkyv::Deserialize::deserialize(archived, &mut rkyv::Infallible).unwrap();
            prop_assert_eq!(live(&back), reference(&a));
        }
    }

    #[test]
    fn inline_heap_transitions_at_128_and_129() {
        let mut set = SuBitGraph::empty(0);
        for i in 0..129 {
            set.push(i % 2 == 0);
            assert_eq!(set.words.heap.is_some(), set.size() > 128);
            live(&set);
        }
        let full = set.clone();
        assert_eq!(full.included_iter().next_back(), Some(Hedge(128)));
        assert_eq!(set.pop(), Some(true));
        assert!(set.words.heap.is_none(), "128 IDs fit inline");
        assert_eq!(set.pop(), Some(false));
        live(&set);

        let mut long = full.clone();
        let tail = long.split_off(Hedge(128));
        assert!(long.words.heap.is_none() && tail.words.heap.is_none());
        assert_eq!(long.size(), 128);
        assert_eq!(long.clone().join(tail), full);

        let mut wide = SuBitGraph::full(129);
        wide.intersect_with(&SuBitGraph::full(128));
        assert_eq!(
            live(&wide),
            reference(&[[true; 128].as_slice(), &[false]].concat())
        );
        let mut narrow = SuBitGraph::full(128);
        narrow.union_with(&SuBitGraph::full(129));
        assert_eq!(narrow, SuBitGraph::full(128));
    }

    #[test]
    fn clones_are_independent() {
        for len in [34, 300] {
            let set = SuBitGraph::full(len);
            let mut clone = set.clone();
            clone.sub(Hedge(5));
            assert!(set.includes(&Hedge(5)));
            assert!(!clone.includes(&Hedge(5)));
            assert_eq!(set, SuBitGraph::full(len));
        }
    }

    #[test]
    fn complement_keeps_dead_bits_clear() {
        let full = !SuBitGraph::empty(34);
        assert_eq!(full.words()[0], tail_mask(34));
        assert_eq!(full, SuBitGraph::full(34));
        assert_eq!(full.n_included(), 34);
        assert_eq!(full.included_iter().count(), 34);
        assert_eq!(full.included_iter().next_back(), Some(Hedge(33)));
        assert!((!full.clone()).is_empty());
        assert!(full.includes(&SuBitGraph::full(34)));
        let state = ahash::RandomState::with_seeds(1, 2, 3, 4);
        assert_eq!(state.hash_one(&full), state.hash_one(SuBitGraph::full(34)));
    }

    #[test]
    fn ord_puts_a_prefix_before_its_extensions() {
        // Extensions with and without IDs in the prefix's last word.
        for prefix in [vec![false; 34], vec![true; 34]] {
            for (len, id) in [(36, Some(35)), (40, None), (100, Some(70))] {
                let mut extension = prefix.clone();
                extension.resize(len, false);
                if let Some(id) = id {
                    extension[id] = true;
                }
                let (p, e) = (subset(&prefix), subset(&extension));
                assert_eq!(
                    reference(&prefix).cmp(&reference(&extension)),
                    Ordering::Less
                );
                assert_eq!(p.cmp(&e), Ordering::Less);
                assert_eq!(e.cmp(&p), Ordering::Greater);
            }
        }
    }

    #[test]
    fn from_usize_matches_powerset_semantics() {
        let set = SuBitGraph::from_usize(0b1011, 3);
        assert_eq!(set.size(), 3);
        assert_eq!(
            set.included_iter().collect::<Vec<_>>(),
            [Hedge(0), Hedge(1)]
        );
    }

    #[test]
    fn words_check_rejects_wrong_counts_and_dead_bits() {
        assert_eq!(SubSetWordsError::check::<u64>(0, &[]), Ok(()));
        assert_eq!(SubSetWordsError::check(64, &[u64::MAX]), Ok(()));
        assert_eq!(SubSetWordsError::check(32, &[u32::MAX]), Ok(()));
        assert_eq!(
            SubSetWordsError::check(0, &[0u64]),
            Err(SubSetWordsError::WordCount {
                len: 0,
                expected: 0,
                found: 1
            })
        );
        assert_eq!(
            SubSetWordsError::check(65, &[0u64]),
            Err(SubSetWordsError::WordCount {
                len: 65,
                expected: 2,
                found: 1
            })
        );
        assert_eq!(
            SubSetWordsError::check(33, &[0, 2u32]),
            Err(SubSetWordsError::DeadBits { len: 33 })
        );
        assert_eq!(
            SubSetWordsError::check(63, &[1u64 << 63]),
            Err(SubSetWordsError::DeadBits { len: 63 })
        );
    }

    #[cfg(feature = "serde")]
    #[test]
    fn serde_format_and_rejection() {
        let mut set = SuBitGraph::empty(34);
        set.add(Hedge(0));
        set.add(Hedge(33));
        assert_eq!(
            serde_json::to_string(&set).unwrap(),
            r#"{"len":34,"words":[1,2]}"#
        );
        let parse = serde_json::from_str::<SuBitGraph>;
        assert_eq!(parse(r#"{"len":34,"words":[1,2]}"#).unwrap(), set);
        assert_eq!(
            parse(r#"{"len":0,"words":[]}"#).unwrap(),
            SuBitGraph::empty(0)
        );
        for malformed in [
            // Wrong word counts.
            r#"{"len":34,"words":[1]}"#,
            r#"{"len":34,"words":[1,2,0]}"#,
            r#"{"len":0,"words":[0]}"#,
            // Bits set past the length.
            r#"{"len":3,"words":[8]}"#,
            r#"{"len":34,"words":[1,4]}"#,
            // Words wider than 32 bits.
            r#"{"len":64,"words":[4294967296,0]}"#,
            r#"{"len":64,"words":[-1,0]}"#,
            r#"{"len":34}"#,
            // Unknown fields.
            r#"{"len":34,"words":[1,2],"id":null}"#,
            // The former `bitvec` encoding.
            r#"{"set":{"order":"bitvec::order::Lsb0","head":{"width":64,"index":0},"bits":3,"data":[5]},"id":null}"#,
        ] {
            assert!(parse(malformed).is_err(), "accepted {malformed}");
        }
    }

    #[cfg(feature = "rkyv")]
    #[test]
    fn rkyv_validation_rejects_malformed_archives() {
        let mut set = SuBitGraph::empty(100);
        for h in [3, 40, 63, 64, 99] {
            set.add(Hedge(h));
        }
        let bytes = rkyv::to_bytes::<_, 256>(&set).unwrap();
        // SAFETY: `bytes` was just archived from a `SuBitGraph`.
        let archived = unsafe { rkyv::archived_root::<SuBitGraph>(&bytes) };
        let len_at = std::ptr::addr_of!(archived.len) as usize - bytes.as_ptr() as usize;
        let with_len = |len: usize| {
            let mut patched = bytes.clone();
            let len = len as rkyv::FixedUsize;
            patched[len_at..len_at + size_of_val(&len)].copy_from_slice(&len.to_ne_bytes());
            patched
        };

        for wrong_count in [0, 64, 129, 200] {
            let patched = with_len(wrong_count);
            assert!(rkyv::check_archived_root::<SuBitGraph>(&patched).is_err());
        }
        let dead_bit = with_len(99);
        assert!(rkyv::check_archived_root::<SuBitGraph>(&dead_bit).is_err());

        // Deserializing without validation still yields the canonical set.
        // SAFETY: only the length changed, so the archive is structurally valid.
        let unchecked = unsafe { rkyv::archived_root::<SuBitGraph>(&dead_bit) };
        let back: SuBitGraph =
            rkyv::Deserialize::deserialize(unchecked, &mut rkyv::Infallible).unwrap();
        assert_eq!(back.size(), 99);
        live(&back);

        // A longer length in the same words is valid.
        let longer = with_len(128);
        let archived = rkyv::check_archived_root::<SuBitGraph>(&longer).unwrap();
        assert_eq!(archived.size(), 128);
        assert!(archived.contains_index(99) && !archived.contains_index(100));
    }
}
