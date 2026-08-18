use std::collections::HashMap;
use std::fmt::{Display, Formatter};

use super::abstract_index::AbstractIndex;
use super::representation::LibraryRep;
use super::slot::{AbsInd, DummyAind, IsAbstractSlot, Slot};
use super::{OrderedStructure, PermutedStructure, TensorStructure};

/// An occurrence-local identifier for an unresolved tensor port.
///
/// Open port identifiers are metadata only. They are canonicalized in logical
/// port order and are never serialized into a Symbolica atom.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct OpenPortId(pub usize);

impl Display for OpenPortId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "open_{}", self.0)
    }
}

/// An explicit abstract index or an unresolved occurrence-local port.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum PartialIndex<A> {
    Explicit(A),
    Open(OpenPortId),
}

impl<A: Display> Display for PartialIndex<A> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Explicit(index) => index.fmt(f),
            Self::Open(id) => id.fmt(f),
        }
    }
}

impl<A: AbsInd> AbsInd for PartialIndex<A> {}

pub type PartialSlot = Slot<LibraryRep, PartialIndex<AbstractIndex>>;
pub type PartialStructure =
    PermutedStructure<OrderedStructure<LibraryRep, PartialIndex<AbstractIndex>>>;

/// Operations that preserve the logical order hidden by `OrderedStructure`'s
/// canonical sorting permutations.
pub trait PartialStructureExt: Sized {
    fn from_logical_slots(slots: impl IntoIterator<Item = PartialSlot>) -> Self;
    fn logical_slots(&self) -> Vec<PartialSlot>;
    fn canonicalize_open_ports(&self) -> Self;
    fn open_positions(&self) -> Vec<usize>;
    fn materialize_open_ports(
        &self,
        replacements: &HashMap<OpenPortId, AbstractIndex>,
    ) -> PermutedStructure<OrderedStructure<LibraryRep, AbstractIndex>>;
    fn materialize_all_open_ports(
        &self,
    ) -> (
        PermutedStructure<OrderedStructure<LibraryRep, AbstractIndex>>,
        HashMap<OpenPortId, AbstractIndex>,
    );
}

impl PartialStructureExt for PartialStructure {
    fn from_logical_slots(slots: impl IntoIterator<Item = PartialSlot>) -> Self {
        let slots = slots.into_iter().collect::<Vec<_>>();
        let mut next = 0;
        let slots = slots.into_iter().map(|mut slot| {
            if matches!(slot.aind, PartialIndex::Open(_)) {
                slot.set_aind(PartialIndex::Open(OpenPortId(next)));
                next += 1;
            }
            slot
        });
        OrderedStructure::new(slots.collect())
    }

    fn logical_slots(&self) -> Vec<PartialSlot> {
        let canonical = self.structure.external_structure();
        let representation_sorted = self.index_permutation.apply_slice_inv(&canonical);
        self.rep_permutation.apply_slice_inv(&representation_sorted)
    }

    fn canonicalize_open_ports(&self) -> Self {
        Self::from_logical_slots(self.logical_slots())
    }

    fn open_positions(&self) -> Vec<usize> {
        self.logical_slots()
            .iter()
            .enumerate()
            .filter_map(|(position, slot)| {
                matches!(slot.aind, PartialIndex::Open(_)).then_some(position)
            })
            .collect()
    }

    fn materialize_open_ports(
        &self,
        replacements: &HashMap<OpenPortId, AbstractIndex>,
    ) -> PermutedStructure<OrderedStructure<LibraryRep, AbstractIndex>> {
        let logical = self.logical_slots().into_iter().map(|slot| {
            let index = match slot.aind {
                PartialIndex::Explicit(index) => index,
                PartialIndex::Open(id) => replacements[&id],
            };
            slot.rep().slot(index)
        });
        OrderedStructure::new(logical.collect())
    }

    fn materialize_all_open_ports(
        &self,
    ) -> (
        PermutedStructure<OrderedStructure<LibraryRep, AbstractIndex>>,
        HashMap<OpenPortId, AbstractIndex>,
    ) {
        let replacements = self
            .logical_slots()
            .into_iter()
            .filter_map(|slot| match slot.aind {
                PartialIndex::Open(id) => Some((id, AbstractIndex::new_dummy())),
                PartialIndex::Explicit(_) => None,
            })
            .collect::<HashMap<_, _>>();
        (self.materialize_open_ports(&replacements), replacements)
    }
}

impl PartialIndex<AbstractIndex> {
    pub fn explicit(index: AbstractIndex) -> Self {
        Self::Explicit(index)
    }

    pub fn open(position: usize) -> Self {
        Self::Open(OpenPortId(position))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::dimension::Dimension;
    use crate::structure::representation::{ExtendibleReps, RepName};

    #[test]
    fn canonicalizes_open_ids_in_logical_order() {
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let structure = PartialStructure::from_logical_slots([
            rep.slot(PartialIndex::Open(OpenPortId(9))),
            rep.slot(PartialIndex::Open(OpenPortId(2))),
        ]);

        assert_eq!(
            structure
                .logical_slots()
                .into_iter()
                .map(|slot| slot.aind)
                .collect::<Vec<_>>(),
            vec![PartialIndex::open(0), PartialIndex::open(1)]
        );
    }

    #[test]
    fn logical_slots_reverse_both_sorting_permutations() {
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let logical = vec![
            mink.slot(PartialIndex::open(0)),
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(8))),
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(3))),
        ];
        let structure = PartialStructure::from_logical_slots(logical.clone());

        assert_eq!(structure.logical_slots(), logical);
    }
}
