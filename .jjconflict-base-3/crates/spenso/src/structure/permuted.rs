use std::fmt::Display;

use bincode_trait_derive::{Decode, Encode};
use linnet::permutation::Permutation;
use serde::{Deserialize, Serialize};

use super::{HasName, StructureError, TensorStructure, slot::IsAbstractSlot};

/// The axis order that moves logical input into representation-grouped storage.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Deserialize, Serialize, Encode, Decode)]
pub struct RepresentationOrder(Permutation);

/// The axis order that moves representation-grouped storage into canonical storage.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Deserialize, Serialize, Encode, Decode)]
pub struct IndexOrder(Permutation);

/// Immutable history of the two stages used to canonicalize a tensor structure.
///
/// Representation ordering belongs to the input boundary. Both fields are history
/// only; later index-only movement is represented by [`Reindexed`].
#[derive(Debug, Clone, PartialEq, Eq, Hash, Deserialize, Serialize, Encode, Decode)]
pub struct CanonicalLayout {
    representation: RepresentationOrder,
    index: IndexOrder,
}

impl CanonicalLayout {
    pub fn identity(order: usize) -> Self {
        Self::from_permutations(Permutation::id(order), Permutation::id(order))
    }

    pub(crate) fn from_permutations(representation: Permutation, index: Permutation) -> Self {
        Self {
            representation: RepresentationOrder(representation),
            index: IndexOrder(index),
        }
    }

    pub(crate) fn representation_permutation(&self) -> &Permutation {
        &self.representation.0
    }

    pub(crate) fn index_permutation(&self) -> &Permutation {
        &self.index.0
    }

    pub fn order(&self) -> usize {
        self.representation.0.map().len()
    }

    pub fn logical_to_representation<T: Clone>(&self, values: &[T]) -> Vec<T> {
        self.representation.0.apply_slice(values)
    }

    pub fn representation_to_logical<T: Clone>(&self, values: &[T]) -> Vec<T> {
        self.representation.0.apply_slice_inv(values)
    }

    pub fn representation_to_canonical<T: Clone>(&self, values: &[T]) -> Vec<T> {
        self.index.0.apply_slice(values)
    }

    pub fn canonical_to_representation<T: Clone>(&self, values: &[T]) -> Vec<T> {
        self.index.0.apply_slice_inv(values)
    }

    pub fn logical_to_canonical<T: Clone>(&self, values: &[T]) -> Vec<T> {
        let representation_order = self.logical_to_representation(values);
        self.representation_to_canonical(&representation_order)
    }

    pub fn canonical_to_logical<T: Clone>(&self, values: &[T]) -> Vec<T> {
        let representation_order = self.canonical_to_representation(values);
        self.representation_to_logical(&representation_order)
    }
}

/// A value in canonical order together with the immutable history of that order.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Deserialize, Serialize, Encode, Decode)]
pub struct Canonicalized<S> {
    canonical: S,
    layout: CanonicalLayout,
}

impl<S: Display> Display for Canonicalized<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.canonical.fmt(f)
    }
}

impl<S> Canonicalized<S> {
    pub fn identity(canonical: S) -> Self
    where
        S: TensorStructure,
    {
        let layout = CanonicalLayout::identity(canonical.order());
        Self { canonical, layout }
    }

    pub(crate) fn from_parts(canonical: S, layout: CanonicalLayout) -> Self {
        Self { canonical, layout }
    }

    pub fn canonical(&self) -> &S {
        &self.canonical
    }

    pub fn into_canonical(self) -> S {
        self.canonical
    }

    pub fn layout(&self) -> &CanonicalLayout {
        &self.layout
    }

    pub fn into_parts(self) -> (S, CanonicalLayout) {
        (self.canonical, self.layout)
    }

    pub fn map_canonical<U>(self, f: impl FnOnce(S) -> U) -> Canonicalized<U> {
        Canonicalized {
            canonical: f(self.canonical),
            layout: self.layout,
        }
    }
}

impl<S: HasName> HasName for Canonicalized<S> {
    type Name = S::Name;
    type Args = S::Args;

    fn name(&self) -> Option<Self::Name> {
        self.canonical.name()
    }

    fn args(&self) -> Option<Self::Args> {
        self.canonical.args()
    }

    fn set_name(&mut self, name: Self::Name) {
        self.canonical.set_name(name);
    }
}

/// A checked index-only storage permutation.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct PendingIndexPermutation(Permutation);

impl PendingIndexPermutation {
    pub(crate) fn checked<S: IsAbstractSlot>(
        source: &[S],
        target: &[S],
        permutation: Permutation,
    ) -> Result<Self, StructureError> {
        if source.len() != target.len() || permutation.map().len() != source.len() {
            return Err(StructureError::InvalidIndexPermutation(
                "source, target, and permutation orders differ".to_string(),
            ));
        }

        for (source_position, target_position) in permutation.map().iter().copied().enumerate() {
            let source_slot = source[source_position];
            let target_slot = target[target_position];
            if source_slot.rep() != target_slot.rep() || source_slot.dim() != target_slot.dim() {
                return Err(StructureError::InvalidIndexPermutation(format!(
                    "axis {source_position} maps across representation or dimension",
                )));
            }
        }

        Ok(Self(permutation))
    }

    pub fn apply_slice<T: Clone>(&self, values: &[T]) -> Vec<T> {
        self.0.apply_slice(values)
    }

    pub fn apply_slice_inverse<T: Clone>(&self, values: &[T]) -> Vec<T> {
        self.0.apply_slice_inv(values)
    }

    pub fn is_identity(&self) -> bool {
        self.0.is_identity()
    }

    pub(crate) fn permutation(&self) -> &Permutation {
        &self.0
    }
}

/// A value whose metadata is canonical but whose payload still uses the source index order.
#[must_use = "apply the pending index permutation before using the reindexed value"]
#[derive(Debug, PartialEq, Eq, Hash)]
pub struct Reindexed<T> {
    target: T,
    pending: PendingIndexPermutation,
}

impl<T> Reindexed<T> {
    pub(crate) fn from_parts(target: T, pending: PendingIndexPermutation) -> Self {
        Self { target, pending }
    }

    pub fn map_target<U>(self, f: impl FnOnce(T) -> U) -> Reindexed<U> {
        Reindexed {
            target: f(self.target),
            pending: self.pending,
        }
    }

    pub fn apply(self) -> <T as ApplyPendingIndexPermutation>::Output
    where
        T: ApplyPendingIndexPermutation,
    {
        self.target.apply_pending_index_permutation(&self.pending)
    }
}

/// Moves a payload from its previous index order into its canonical target order.
pub trait ApplyPendingIndexPermutation {
    type Output;

    fn apply_pending_index_permutation(self, pending: &PendingIndexPermutation) -> Self::Output;
}

impl<T: ApplyPendingIndexPermutation> ApplyPendingIndexPermutation for Canonicalized<T> {
    type Output = Canonicalized<T::Output>;

    fn apply_pending_index_permutation(self, pending: &PendingIndexPermutation) -> Self::Output {
        self.map_canonical(|canonical| canonical.apply_pending_index_permutation(pending))
    }
}

/// Builds the identity tensor associated with a pair of compatible slots.
pub trait TensorIdentity {
    type Id;
    type IdSlot;

    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id;
}

// The previous transition wrapper carried these debugging probes while composing
// representation and index permutations:
// println!("Layout:{layout:?}");
// println!("Pending:{pending:?}");
