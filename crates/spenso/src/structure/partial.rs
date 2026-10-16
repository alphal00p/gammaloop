use std::collections::HashMap;
use std::fmt::{Display, Formatter};

use super::abstract_index::AbstractIndex;
use super::representation::{LibraryRep, Representation};
use super::slot::{AbsInd, DummyAind, IsAbstractSlot, Slot};
use super::{Canonicalized, OrderedStructure, StructureError, TensorStructure};

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
    Canonicalized<OrderedStructure<LibraryRep, PartialIndex<AbstractIndex>>>;

/// Canonical signature access and operations on the stored logical axis order.
pub trait PartialStructureExt: Sized {
    fn from_logical_slots(slots: impl IntoIterator<Item = PartialSlot>) -> Self;
    fn logical_slots(&self) -> Vec<PartialSlot>;
    /// Return all explicitly indexed slots in canonical external-slot order.
    ///
    /// Fails at the first unresolved axis, without allocating an index or
    /// dropping axes. Scalars return an empty vector.
    fn slots(&self) -> Result<Vec<Slot<LibraryRep, AbstractIndex>>, StructureError>;
    /// Return all unresolved representations in canonical representation order.
    ///
    /// Fails at the first explicitly indexed axis; indices are never discarded.
    /// Different spaces and dimensions are allowed. Scalars return an empty vector.
    fn representations(&self) -> Result<Vec<Representation<LibraryRep>>, StructureError>;
    fn canonicalize_open_ports(&self) -> Self;
    fn open_positions(&self) -> Vec<usize>;
    fn materialize_open_ports(
        &self,
        replacements: &HashMap<OpenPortId, AbstractIndex>,
    ) -> Canonicalized<OrderedStructure<LibraryRep, AbstractIndex>>;
    fn materialize_all_open_ports(
        &self,
    ) -> (
        Canonicalized<OrderedStructure<LibraryRep, AbstractIndex>>,
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
        let canonical = self.canonical().external_structure();
        self.layout().canonical_to_logical(&canonical)
    }

    fn slots(&self) -> Result<Vec<Slot<LibraryRep, AbstractIndex>>, StructureError> {
        self.canonical().slots()
    }

    fn representations(&self) -> Result<Vec<Representation<LibraryRep>>, StructureError> {
        self.canonical().representations()
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
    ) -> Canonicalized<OrderedStructure<LibraryRep, AbstractIndex>> {
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
        Canonicalized<OrderedStructure<LibraryRep, AbstractIndex>>,
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

impl OrderedStructure<LibraryRep, PartialIndex<AbstractIndex>> {
    /// Return every external slot in canonical order, or fail at the first unresolved axis.
    pub fn slots(&self) -> Result<Vec<Slot<LibraryRep, AbstractIndex>>, StructureError> {
        self.external_structure()
            .into_iter()
            .enumerate()
            .map(|(axis, slot)| match slot.aind {
                PartialIndex::Explicit(index) => Ok(slot.rep().slot(index)),
                PartialIndex::Open(_) => Err(StructureError::ExpectedExplicitSlot { axis }),
            })
            .collect()
    }

    /// Return every unresolved representation in canonical order, or fail at the first indexed axis.
    pub fn representations(&self) -> Result<Vec<Representation<LibraryRep>>, StructureError> {
        self.external_structure()
            .into_iter()
            .enumerate()
            .map(|(axis, slot)| match slot.aind {
                PartialIndex::Open(_) => Ok(slot.rep()),
                PartialIndex::Explicit(_) => Err(StructureError::ExpectedRepresentation { axis }),
            })
            .collect()
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
    fn homogeneous_accessors_use_canonical_order() {
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let slots = [
            mink.slot(AbstractIndex::Normal(9)),
            euc.slot(AbstractIndex::Normal(2)),
        ];
        let explicit = PartialStructure::from_logical_slots(
            slots
                .iter()
                .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind))),
        );
        assert_eq!(explicit.slots().unwrap(), [slots[1], slots[0]]);
        assert!(matches!(
            explicit.representations(),
            Err(StructureError::ExpectedRepresentation { axis: 0 })
        ));
        let open = PartialStructure::from_logical_slots([
            mink.slot(PartialIndex::open(9)),
            euc.slot(PartialIndex::open(2)),
        ]);
        assert_eq!(open.representations().unwrap(), [euc, mink]);
        assert!(matches!(
            open.slots(),
            Err(StructureError::ExpectedExplicitSlot { axis: 0 })
        ));
    }

    #[test]
    fn homogeneous_accessors_reject_mixed_axes_without_losing_information() {
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let structure = PartialStructure::from_logical_slots([
            rep.slot(PartialIndex::Explicit(AbstractIndex::Normal(7))),
            rep.slot(PartialIndex::open(8)),
        ]);
        let before = structure.clone();
        assert!(matches!(
            structure.slots(),
            Err(StructureError::ExpectedExplicitSlot { axis: 1 })
        ));
        assert!(matches!(
            structure.representations(),
            Err(StructureError::ExpectedRepresentation { axis: 0 })
        ));
        assert_eq!(structure, before);
        let reversed =
            PartialStructure::from_logical_slots(structure.logical_slots().into_iter().rev());
        assert!(matches!(
            reversed.slots(),
            Err(StructureError::ExpectedExplicitSlot { axis: 1 })
        ));
        assert!(matches!(
            reversed.representations(),
            Err(StructureError::ExpectedRepresentation { axis: 0 })
        ));
    }

    #[test]
    fn scalar_structure_supports_both_homogeneous_accessors() {
        let scalar = PartialStructure::from_logical_slots([]);
        assert!(scalar.slots().unwrap().is_empty());
        assert!(scalar.representations().unwrap().is_empty());
    }

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
