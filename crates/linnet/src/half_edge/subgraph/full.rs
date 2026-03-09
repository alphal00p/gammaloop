use std::ops::{Range, RangeFrom, RangeInclusive, RangeTo, RangeToInclusive};

use crate::half_edge::{
    involution::{Hedge, HedgePair},
    nodestore::NodeStorageOps,
    subgraph::{SuBitGraph, SubGraphLike},
    HedgeGraph,
};

use super::{Inclusion, SubSetLike};

#[derive(Clone, Debug, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
/// Represents a subgraph that is either completely full (contains all possible half-edges
/// in the graph context) or entirely empty.
///
/// This struct provides a lightweight way to specify these two common boundary conditions
/// for subgraph operations. The state (full or empty) is determined by the sign
/// of the `size` field.
pub struct FullOrEmpty {
    /// Stores the size of the graph context.
    /// If `size` is positive, it represents a full subgraph of that many half-edges.
    /// If `size` is negative or zero, it represents an empty subgraph; the absolute
    /// value can still indicate the size of the parent graph context.
    size: isize,
}

impl FullOrEmpty {
    pub fn full(size: usize) -> FullOrEmpty {
        FullOrEmpty {
            size: size as isize,
        }
    }
    pub fn empty(size: usize) -> FullOrEmpty {
        FullOrEmpty {
            size: -(size as isize),
        }
    }
}

impl Inclusion<Hedge> for FullOrEmpty {
    fn includes(&self, hedge_id: &Hedge) -> bool {
        (hedge_id.0 as isize) < self.size
    }

    fn intersects(&self, other: &Hedge) -> bool {
        (other.0 as isize) < self.size
    }
}

impl Inclusion<HedgePair> for FullOrEmpty {
    fn includes(&self, hedge_id: &HedgePair) -> bool {
        match hedge_id {
            HedgePair::Split { source, sink, .. } | HedgePair::Paired { source, sink } => {
                (source.0 as isize) < self.size && (sink.0 as isize) < self.size
            }
            HedgePair::Unpaired { hedge, .. } => (hedge.0 as isize) < self.size,
        }
    }

    fn intersects(&self, other: &HedgePair) -> bool {
        match other {
            HedgePair::Split { source, sink, .. } | HedgePair::Paired { source, sink } => {
                (source.0 as isize) < self.size && (sink.0 as isize) < self.size
            }
            HedgePair::Unpaired { hedge, .. } => (hedge.0 as isize) < self.size,
        }
    }
}

impl Inclusion<FullOrEmpty> for FullOrEmpty {
    fn includes(&self, other: &FullOrEmpty) -> bool {
        // true
        self.size == other.size
    }

    fn intersects(&self, other: &FullOrEmpty) -> bool {
        // true
        self.size == other.size
    }
}

impl Inclusion<SuBitGraph> for FullOrEmpty {
    fn includes(&self, other: &SuBitGraph) -> bool {
        self.size == other.size() as isize
    }

    fn intersects(&self, other: &SuBitGraph) -> bool {
        self.size == other.size() as isize
    }
}

impl Inclusion<Range<Hedge>> for FullOrEmpty {
    fn includes(&self, other: &Range<Hedge>) -> bool {
        if let Ok(size) = usize::try_from(self.size) {
            size >= other.end.0
        } else {
            other.start >= other.end
        }
    }

    fn intersects(&self, other: &Range<Hedge>) -> bool {
        if let Ok(size) = usize::try_from(self.size) {
            size >= other.end.0
        } else {
            other.start >= other.end
        }
    }
}

impl Inclusion<RangeTo<Hedge>> for FullOrEmpty {
    fn includes(&self, other: &RangeTo<Hedge>) -> bool {
        if let Ok(size) = usize::try_from(self.size) {
            size >= other.end.0
        } else {
            false
        }
    }

    fn intersects(&self, other: &RangeTo<Hedge>) -> bool {
        if let Ok(size) = usize::try_from(self.size) {
            size >= other.end.0
        } else {
            false
        }
    }
}

impl Inclusion<RangeToInclusive<Hedge>> for FullOrEmpty {
    fn includes(&self, other: &RangeToInclusive<Hedge>) -> bool {
        if let Ok(size) = usize::try_from(self.size) {
            size > other.end.0
        } else {
            false
        }
    }

    fn intersects(&self, other: &RangeToInclusive<Hedge>) -> bool {
        if let Ok(size) = usize::try_from(self.size) {
            size > other.end.0
        } else {
            false
        }
    }
}

impl Inclusion<RangeFrom<Hedge>> for FullOrEmpty {
    fn includes(&self, other: &RangeFrom<Hedge>) -> bool {
        if let Ok(size) = usize::try_from(self.size) {
            size > other.start.0
        } else {
            false
        }
    }

    fn intersects(&self, other: &RangeFrom<Hedge>) -> bool {
        if let Ok(size) = usize::try_from(self.size) {
            size > other.start.0
        } else {
            false
        }
    }
}

impl Inclusion<RangeInclusive<Hedge>> for FullOrEmpty {
    fn includes(&self, other: &RangeInclusive<Hedge>) -> bool {
        if let Ok(size) = usize::try_from(self.size) {
            size > other.end().0
        } else {
            other.start() > other.end()
        }
    }

    fn intersects(&self, other: &RangeInclusive<Hedge>) -> bool {
        if let Ok(size) = usize::try_from(self.size) {
            size > other.end().0
        } else {
            other.start() > other.end()
        }
    }
}

/// An iterator that yields [`Hedge`] identifiers over a given numerical range.
///
/// This is used as the `BaseIter` for the [`FullOrEmpty`] subgraph type when it
/// represents a full graph, allowing iteration over all hedges from `0` up to `size-1`.
pub struct RangeHedgeIter {
    /// The underlying `Range<usize>` that this iterator traverses.
    iter: Range<usize>,
}

impl From<Range<Hedge>> for RangeHedgeIter {
    fn from(value: Range<Hedge>) -> Self {
        RangeHedgeIter {
            iter: value.start.0..value.end.0,
        }
    }
}

impl DoubleEndedIterator for RangeHedgeIter {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.next_back().map(Hedge)
    }
}

impl Iterator for RangeHedgeIter {
    type Item = Hedge;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(Hedge)
    }
}

impl SubGraphLike for FullOrEmpty {
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
impl SubSetLike<Hedge> for FullOrEmpty {
    type Base = FullOrEmpty;
    type BaseIter<'a> = RangeHedgeIter;

    fn has_greater(&self, _hedge: Hedge) -> bool {
        false
    }

    fn size(&self) -> usize {
        self.size.unsigned_abs()
    }

    fn join_mut(&mut self, other: Self) {
        if self.is_empty() && other.is_empty() {
            self.size -= other.size
        } else if !self.is_empty() && !other.is_empty() {
            self.size += other.size
        }
    }

    fn included_iter(&self) -> Self::BaseIter<'_> {
        if self.is_empty() {
            (Hedge(1)..Hedge(0)).into()
        } else {
            (Hedge(0)..Hedge(self.size as usize)).into()
        }
    }

    fn n_included(&self) -> usize {
        self.size.try_into().unwrap_or(0)
    }

    fn empty(size: usize) -> Self {
        FullOrEmpty {
            size: -(size as isize),
        }
    }

    fn included(&self) -> &FullOrEmpty {
        self
    }

    fn string_label(&self) -> String {
        if self.is_empty() {
            "∅".into()
        } else {
            "full".into()
        }
    }
    fn from_base62(label: &str, size: usize) -> Option<Self> {
        if label.is_empty() {
            Some(FullOrEmpty::empty(size))
        } else {
            Some(FullOrEmpty {
                size: size as isize,
            })
        }
    }
    fn is_empty(&self) -> bool {
        self.size <= 0
    }
}
