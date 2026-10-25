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

    /// Returns the canonical storage position of a logical axis in constant time.
    ///
    /// # Panics
    ///
    /// Panics if `position >= self.order()`.
    ///
    /// # Examples
    ///
    /// ```
    /// use spenso::structure::CanonicalLayout;
    ///
    /// let layout = CanonicalLayout::identity(3);
    /// assert_eq!(layout.logical_to_canonical_position(1), 1);
    /// ```
    pub fn logical_to_canonical_position(&self, position: usize) -> usize {
        let representation_position = self.representation.0.map()[position];
        self.index.0.map()[representation_position]
    }

    /// Returns the logical position of a canonical storage axis in constant time.
    ///
    /// # Panics
    ///
    /// Panics if `position >= self.order()`.
    ///
    /// # Examples
    ///
    /// ```
    /// use spenso::structure::CanonicalLayout;
    ///
    /// let layout = CanonicalLayout::identity(3);
    /// assert_eq!(layout.canonical_to_logical_position(1), 1);
    /// ```
    pub fn canonical_to_logical_position(&self, position: usize) -> usize {
        let representation_position = self.index.0.inv()[position];
        self.representation.0.inv()[representation_position]
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
        (0..self.order())
            .map(|position| values[self.logical_to_canonical_position(position)].clone())
            .collect()
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_layout_preserves_noncommuting_axis_orders() {
        let layout = CanonicalLayout::from_permutations(
            Permutation::from_map(vec![1, 2, 0, 3]),
            Permutation::from_map(vec![0, 3, 2, 1]),
        );
        let logical = [
            String::from("a"),
            String::from("b"),
            String::from("c"),
            String::from("d"),
        ];
        let canonical = [
            String::from("c"),
            String::from("d"),
            String::from("b"),
            String::from("a"),
        ];

        assert_eq!(layout.logical_to_canonical(&logical), canonical);
        assert_eq!(layout.canonical_to_logical(&canonical), logical);
        for (logical_position, canonical_position) in [3, 2, 0, 1].into_iter().enumerate() {
            assert_eq!(
                layout.logical_to_canonical_position(logical_position),
                canonical_position,
            );
            assert_eq!(
                layout.canonical_to_logical_position(canonical_position),
                logical_position,
            );
        }
    }

    #[test]
    fn identity_layout_preserves_empty_single_and_multiple_axes() {
        for order in 0..=3 {
            let layout = CanonicalLayout::identity(order);
            let values: Vec<_> = (0..order).map(|position| position.to_string()).collect();
            assert_eq!(layout.canonical_to_logical(&values), values);
            assert_eq!(layout.logical_to_canonical(&values), values);
            for position in 0..order {
                assert_eq!(layout.logical_to_canonical_position(position), position);
                assert_eq!(layout.canonical_to_logical_position(position), position);
            }
        }
    }

    #[test]
    #[should_panic]
    fn logical_position_rejects_out_of_bounds_axis() {
        CanonicalLayout::identity(3).logical_to_canonical_position(3);
    }

    #[test]
    #[should_panic]
    fn canonical_position_rejects_out_of_bounds_axis() {
        CanonicalLayout::identity(0).canonical_to_logical_position(0);
    }
}
