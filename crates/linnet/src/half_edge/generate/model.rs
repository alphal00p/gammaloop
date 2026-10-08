//! Dense form of the generation input: stub types, signature lookup tables, external-leg classes
//! and edge labels.
//!
//! A stub type is `3 * field + direction` (undirected, incoming, outgoing), with fields ranked in
//! their `Ord` order. The stubs at a vertex are a [`StubSet`]: the sorted multiset packed with the
//! vertex degree into one `u128` key, looked up in two tables built from the signatures: every
//! sub-multiset of a signature of degree `d` (partial vertices that can still be completed) and
//! every full signature (the vertex rules with exactly these stubs).
//!
//! External legs are packed by label class (interchangeable legs: equal label and stub) on complete
//! topologies, and by slot class while the topology search places them, where distinct labels with
//! a common stub can share an anonymous class.

use ahash::{AHashMap, AHashSet};

use super::GenerationError;

/// Largest vertex degree: a [`StubSet`] holds 15 stubs and its degree.
pub(crate) const MAX_DEGREE: usize = 15;
/// Largest number of fields: stub types are bytes and `StubSet` stores them shifted by one.
pub(crate) const MAX_FIELDS: usize = 84;
pub(crate) const UNDIRECTED: u8 = 0;
pub(crate) const INCOMING: u8 = 1;
pub(crate) const OUTGOING: u8 = 2;

/// Sorted multiset of stub types at a vertex of known degree; byte `15` is the degree and bytes
/// `0..len` the sorted types plus one, so the packed key is exact.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct StubSet {
    bytes: [u8; 16],
}

impl StubSet {
    pub(crate) fn new(degree: u8) -> Self {
        let mut bytes = [0; 16];
        bytes[15] = degree;
        StubSet { bytes }
    }

    pub(crate) fn insert(&mut self, stub: u8) {
        let b = stub + 1;
        let i = self.bytes[..15]
            .iter()
            .position(|&x| x == 0 || x > b)
            .expect("a vertex holds at most its degree in stubs");
        self.bytes.copy_within(i..14, i + 1);
        self.bytes[i] = b;
    }

    pub(crate) fn remove(&mut self, stub: u8) {
        let i = self.bytes[..15]
            .iter()
            .position(|&x| x == stub + 1)
            .expect("removed stub is present");
        self.bytes.copy_within(i + 1..15, i);
        self.bytes[14] = 0;
    }

    /// Holds as many stubs as its degree.
    pub(crate) fn is_full(&self) -> bool {
        let degree = self.bytes[15] as usize;
        degree == 0 || self.bytes[degree - 1] != 0
    }

    pub(crate) fn key(&self) -> u128 {
        u128::from_le_bytes(self.bytes)
    }
}

/// A paired edge between two distinct vertices `u < v` as seen from the pair: the stub it puts at
/// `u` and at `v`, its image when `u` and `v` swap, and how to materialize it.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Label {
    pub(crate) at_u: u8,
    pub(crate) at_v: u8,
    pub(crate) flip: u8,
    pub(crate) field: u8,
    pub(crate) direction: Direction,
}

/// A self-loop: the two stubs it puts at its vertex.
#[derive(Clone, Copy, Debug)]
pub(crate) struct LoopLabel {
    pub(crate) stubs: [u8; 2],
    pub(crate) field: u8,
    pub(crate) directed: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Direction {
    Undirected,
    /// From the lower to the higher vertex of the pair.
    Forward,
    Backward,
}

/// External legs grouped into classes of interchangeable legs (equal label and stub). A multiset
/// of legs is packed into a `u64`: class `c` owns `width[c]` bits at `offset[c]`, enough for its
/// multiplicity, so sums and differences of packed multisets never carry.
#[derive(Clone, Debug, Default)]
pub(crate) struct Legs {
    pub(crate) stub: Vec<u8>,
    pub(crate) mult: Vec<u8>,
    offset: Vec<u8>,
    width: Vec<u8>,
    /// All legs, packed.
    pub(crate) all: u64,
    pub(crate) total: u32,
}

impl Legs {
    pub(crate) fn new(classes: &[(u8, u8)]) -> Result<Self, GenerationError> {
        let mut legs = Legs::default();
        let mut offset = 0u32;
        for &(stub, mult) in classes {
            let width = u8::BITS - mult.leading_zeros();
            if offset + width > u64::BITS {
                return Err(GenerationError::TooManyLegs);
            }
            legs.stub.push(stub);
            legs.mult.push(mult);
            legs.offset.push(offset as u8);
            legs.width.push(width as u8);
            legs.all |= (mult as u64) << offset;
            legs.total += mult as u32;
            offset += width;
        }
        Ok(legs)
    }

    pub(crate) fn n_classes(&self) -> usize {
        self.stub.len()
    }

    pub(crate) fn count(&self, packed: u64, class: usize) -> u8 {
        ((packed >> self.offset[class]) & ((1 << self.width[class]) - 1)) as u8
    }

    pub(crate) fn one(&self, class: usize) -> u64 {
        1 << self.offset[class]
    }

    /// Adding the legs `add` to `used` stays within `total`, class by class.
    pub(crate) fn fits(&self, used: u64, add: u64, total: u64) -> bool {
        (0..self.n_classes()).all(|c| {
            self.count(used, c) as u32 + self.count(add, c) as u32 <= self.count(total, c) as u32
        })
    }

    /// Number of legs in `packed`.
    pub(crate) fn size(&self, packed: u64) -> u32 {
        (0..self.n_classes())
            .map(|c| self.count(packed, c) as u32)
            .sum()
    }

    /// Every sub-multiset of `packed` with at most `max_len` legs, with its size.
    pub(crate) fn subsets(&self, packed: u64, max_len: u32, out: &mut Vec<(u64, u32)>) {
        out.clear();
        out.push((0, 0));
        for c in 0..self.n_classes() {
            let available = self.count(packed, c) as u32;
            for i in 0..out.len() {
                let (base, len) = out[i];
                for k in 1..=available.min(max_len.saturating_sub(len)) {
                    out.push((base + k as u64 * self.one(c), len + k));
                }
            }
        }
    }
}

/// The dense model shared by both stages.
#[derive(Debug)]
pub(crate) struct Model {
    /// Distinct signature degrees, ascending.
    pub(crate) degrees: Vec<u8>,
    partial: AHashSet<u128>,
    exact: AHashMap<u128, Vec<u32>>,
    pub(crate) labels: Vec<Label>,
    pub(crate) loop_labels: Vec<LoopLabel>,
    /// The external legs by label class, as complete topologies carry them.
    pub(crate) legs: Legs,
    /// The external legs as the topology search places them: three or more label classes of
    /// multiplicity one with a common stub form one anonymous slot class, whose labels are assigned
    /// once the topology is shaped; every other label class is its own slot class. With two such
    /// legs, labelling costs about what the symmetric topologies save.
    pub(crate) slots: Legs,
    /// Per slot class: its label class, `None` for an anonymous one.
    slot_labels: Vec<Option<usize>>,
    /// The labels of anonymous slots in assignment order: `(slot class, label class)`.
    pub(crate) labelling: Vec<(usize, usize)>,
    /// Every vertex has at most one rule for its stubs.
    pub(crate) unique_rules: bool,
}

impl Model {
    /// `signatures` and `legs` are given as stub types; legs as classes `(stub, multiplicity)`.
    pub(crate) fn new(
        n_fields: usize,
        signatures: &[Vec<u8>],
        legs: &[(u8, u8)],
    ) -> Result<Self, GenerationError> {
        if n_fields > MAX_FIELDS {
            return Err(GenerationError::TooManyFields);
        }
        let mut partial = AHashSet::new();
        let mut exact: AHashMap<u128, Vec<u32>> = AHashMap::new();
        let mut degrees = vec![];
        let mut present = vec![false; 3 * n_fields];
        for (i, signature) in signatures.iter().enumerate() {
            let degree = signature.len();
            if degree == 0 || degree > MAX_DEGREE {
                return Err(GenerationError::UnsupportedDegree(degree));
            }
            degrees.push(degree as u8);
            let mut sorted = signature.clone();
            sorted.sort_unstable();
            for &s in &sorted {
                present[s as usize] = true;
            }
            // every sub-multiset: choose how many of each run of equal stubs to keep
            let runs: Vec<(u8, usize)> = sorted
                .chunk_by(|a, b| a == b)
                .map(|run| (run[0], run.len()))
                .collect();
            let mut keep = vec![0usize; runs.len()];
            loop {
                let mut set = StubSet::new(degree as u8);
                for (&(stub, _), &k) in runs.iter().zip(&keep) {
                    for _ in 0..k {
                        set.insert(stub);
                    }
                }
                partial.insert(set.key());
                let Some(r) = (0..runs.len()).find(|&r| keep[r] < runs[r].1) else {
                    exact.entry(set.key()).or_default().push(i as u32);
                    break;
                };
                keep[r] += 1;
                keep[..r].fill(0);
            }
        }
        degrees.sort_unstable();
        degrees.dedup();

        let mut labels = vec![];
        let mut loop_labels = vec![];
        for field in 0..n_fields as u8 {
            let stub = |direction: u8| 3 * field + direction;
            if present[stub(UNDIRECTED) as usize] {
                let s = stub(UNDIRECTED);
                labels.push(Label {
                    at_u: s,
                    at_v: s,
                    flip: labels.len() as u8,
                    field,
                    direction: Direction::Undirected,
                });
                loop_labels.push(LoopLabel {
                    stubs: [s, s],
                    field,
                    directed: false,
                });
            }
            if present[stub(INCOMING) as usize] && present[stub(OUTGOING) as usize] {
                let (i, o) = (stub(INCOMING), stub(OUTGOING));
                let forward = labels.len() as u8;
                labels.push(Label {
                    at_u: o,
                    at_v: i,
                    flip: forward + 1,
                    field,
                    direction: Direction::Forward,
                });
                labels.push(Label {
                    at_u: i,
                    at_v: o,
                    flip: forward,
                    field,
                    direction: Direction::Backward,
                });
                loop_labels.push(LoopLabel {
                    stubs: [i, o],
                    field,
                    directed: true,
                });
            }
        }
        let unique_rules = exact.values().all(|rules| rules.len() == 1);
        // slot classes in the order of their first label class
        let anonymous = |&(stub, mult): &(u8, u8)| {
            mult == 1 && legs.iter().filter(|&&other| other == (stub, 1)).count() > 2
        };
        let mut slots: Vec<(u8, u8)> = vec![];
        let mut slot_labels = vec![];
        let mut labelling = vec![];
        for (c, class) in legs.iter().enumerate() {
            if !anonymous(class) {
                slots.push(*class);
                slot_labels.push(Some(c));
                continue;
            }
            let s = match (0..slots.len())
                .find(|&s| slot_labels[s].is_none() && slots[s].0 == class.0)
            {
                Some(s) => {
                    slots[s].1 += 1;
                    s
                }
                None => {
                    slots.push(*class);
                    slot_labels.push(None);
                    slots.len() - 1
                }
            };
            labelling.push((s, c));
        }
        Ok(Model {
            degrees,
            partial,
            exact,
            labels,
            loop_labels,
            legs: Legs::new(legs)?,
            slots: Legs::new(&slots)?,
            slot_labels,
            labelling,
            unique_rules,
        })
    }

    /// The legs of the labelled slot classes in `slots` (packed by slot class), packed by label
    /// class.
    pub(crate) fn labelled(&self, slots: u64) -> u64 {
        self.slot_labels
            .iter()
            .enumerate()
            .filter_map(|(s, c)| Some(self.slots.count(slots, s) as u64 * self.legs.one((*c)?)))
            .sum()
    }

    /// The legs of the anonymous slot classes in `slots`.
    pub(crate) fn anonymous(&self, slots: u64) -> u64 {
        self.slot_labels
            .iter()
            .enumerate()
            .filter(|(_, c)| c.is_none())
            .map(|(s, _)| self.slots.count(slots, s) as u64 * self.slots.one(s))
            .sum()
    }

    /// The stubs fit into some signature of their degree.
    pub(crate) fn feasible(&self, stubs: &StubSet) -> bool {
        self.partial.contains(&stubs.key())
    }

    /// The signatures with exactly these stubs.
    pub(crate) fn rules(&self, stubs: &StubSet) -> &[u32] {
        self.exact.get(&stubs.key()).map_or(&[], Vec::as_slice)
    }

    /// The legs in `slots` (packed by slot class) fit at a vertex of this degree.
    pub(crate) fn legs_fit(&self, degree: u8, slots: u64) -> bool {
        slots == 0 || self.feasible(&self.leg_stubs(degree, slots))
    }

    /// The stubs of the legs in `slots` at a vertex of this degree.
    pub(crate) fn leg_stubs(&self, degree: u8, slots: u64) -> StubSet {
        let mut stubs = StubSet::new(degree);
        for c in 0..self.slots.n_classes() {
            for _ in 0..self.slots.count(slots, c) {
                stubs.insert(self.slots.stub[c]);
            }
        }
        stubs
    }

    /// Edges can take only one label and vertices one rule: a topology has one diagram.
    pub(crate) fn single_diagram(&self) -> bool {
        self.labels.len() <= 1 && self.loop_labels.len() <= 1 && self.unique_rules
    }
}

#[cfg(test)]
mod tests {
    use super::{Legs, StubSet};

    #[test]
    fn stub_sets_are_sorted_multisets() {
        let mut a = StubSet::new(4);
        for s in [5, 1, 5, 0] {
            a.insert(s);
        }
        let mut b = StubSet::new(4);
        for s in [0, 5, 1, 5] {
            b.insert(s);
        }
        assert_eq!(a.key(), b.key());
        assert!(a.is_full());
        a.remove(5);
        assert!(!a.is_full());
        b.remove(5);
        assert_eq!(a.key(), b.key());
        assert_ne!(a.key(), StubSet::new(3).key());
    }

    #[test]
    fn packed_legs_count_and_enumerate_sub_multisets() {
        // classes of multiplicity 2, 1 and 3
        let legs = Legs::new(&[(0, 2), (1, 1), (3, 3)]).unwrap();
        assert_eq!(legs.total, 6);
        assert_eq!(
            (0..3).map(|c| legs.count(legs.all, c)).collect::<Vec<_>>(),
            [2, 1, 3]
        );
        let mut subsets = vec![];
        legs.subsets(legs.all, u32::MAX, &mut subsets);
        assert_eq!(subsets.len(), 3 * 2 * 4);
        legs.subsets(legs.all, 1, &mut subsets);
        assert_eq!(subsets.len(), 4);
        assert!(subsets
            .iter()
            .all(|&(s, n)| (0..3).map(|c| legs.count(s, c) as u32).sum::<u32>() == n));
    }
}
