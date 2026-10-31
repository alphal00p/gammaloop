use std::fmt::Display;

use bincode::de::Decoder;
use bincode::error::DecodeError;
use bincode::{Decode as BincodeDecode, impl_borrow_decode};
use bincode_trait_derive::{Decode, Encode};
use linnet::permutation::Permutation;
use serde::{Deserialize, Serialize, Serializer};

use super::slot::IsAbstractSlot;
use super::{HasName, StructureError, TensorStructure};

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
/// Canonical results store only their order, without allocating identity maps.
/// Direct Serde and Bincode encodings distinguish this form from mapped history.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Deserialize, Encode)]
#[serde(try_from = "CanonicalStages")]
pub struct CanonicalLayout(CanonicalStages);

#[derive(Debug, Clone, PartialEq, Eq, Hash, Deserialize, Serialize, Encode, Decode)]
enum CanonicalStages {
    Identity {
        order: usize,
    },
    Mapped {
        representation: RepresentationOrder,
        index: IndexOrder,
    },
}

impl TryFrom<CanonicalStages> for CanonicalLayout {
    type Error = &'static str;

    fn try_from(stages: CanonicalStages) -> Result<Self, Self::Error> {
        match stages {
            CanonicalStages::Identity { order } => Ok(Self::identity(order)),
            CanonicalStages::Mapped {
                representation,
                index,
            } => {
                Self::validate_permutation(&representation.0)?;
                Self::validate_permutation(&index.0)?;
                if representation.0.map().len() != index.0.map().len() {
                    return Err("representation and index orders differ");
                }
                // Cancelling nonidentity stages still describe a meaningful
                // intermediate order. Validate both directions before checking
                // identity so corrupt inverse maps cannot be normalized away.
                if representation.0.is_identity() && index.0.is_identity() {
                    return Ok(Self::identity(representation.0.map().len()));
                }
                Ok(Self(CanonicalStages::Mapped {
                    representation,
                    index,
                }))
            }
        }
    }
}

impl Serialize for CanonicalLayout {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.0.serialize(serializer)
    }
}

impl<Context> BincodeDecode<Context> for CanonicalLayout {
    fn decode<D: Decoder<Context = Context>>(decoder: &mut D) -> Result<Self, DecodeError> {
        Self::try_from(CanonicalStages::decode(decoder)?).map_err(DecodeError::Other)
    }
}

impl_borrow_decode!(CanonicalLayout);

impl CanonicalLayout {
    fn validate_permutation(permutation: &Permutation) -> Result<(), &'static str> {
        let map = permutation.map();
        let inverse = permutation.inv();
        if map.len() != inverse.len() {
            return Err("permutation forward and inverse orders differ");
        }
        for (source, &target) in map.iter().enumerate() {
            if inverse.get(target) != Some(&source) {
                return Err("permutation forward map has an invalid inverse");
            }
        }
        for (target, &source) in inverse.iter().enumerate() {
            if map.get(source) != Some(&target) {
                return Err("permutation inverse map has an invalid forward map");
            }
        }
        Ok(())
    }

    /// Records an already canonical interface without allocating permutation maps.
    ///
    /// ```
    /// use spenso::structure::CanonicalLayout;
    ///
    /// let layout = CanonicalLayout::identity(2);
    /// assert_eq!(layout.canonical_to_logical(&["i", "j"]), ["i", "j"]);
    /// ```
    pub fn identity(order: usize) -> Self {
        Self(CanonicalStages::Identity { order })
    }

    pub(crate) fn from_permutations(representation: Permutation, index: Permutation) -> Self {
        Self::try_from(CanonicalStages::Mapped {
            representation: RepresentationOrder(representation),
            index: IndexOrder(index),
        })
        .expect("canonicalization produces valid permutations of the same order")
    }

    pub(crate) fn representation_is_identity(&self) -> bool {
        match &self.0 {
            CanonicalStages::Identity { .. } => true,
            CanonicalStages::Mapped { representation, .. } => representation.0.is_identity(),
        }
    }

    pub(crate) fn into_index_permutation(self) -> Permutation {
        match self.0 {
            CanonicalStages::Identity { order } => Permutation::id(order),
            CanonicalStages::Mapped { index, .. } => index.0,
        }
    }

    pub fn order(&self) -> usize {
        match &self.0 {
            CanonicalStages::Identity { order } => *order,
            CanonicalStages::Mapped { representation, .. } => representation.0.map().len(),
        }
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
        assert!(position < self.order(), "logical axis is out of bounds");
        match &self.0 {
            CanonicalStages::Identity { .. } => position,
            CanonicalStages::Mapped {
                representation,
                index,
            } => index.0.map()[representation.0.map()[position]],
        }
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
        assert!(position < self.order(), "canonical axis is out of bounds");
        match &self.0 {
            CanonicalStages::Identity { .. } => position,
            CanonicalStages::Mapped {
                representation,
                index,
            } => representation.0.inv()[index.0.inv()[position]],
        }
    }

    pub fn logical_to_representation<T: Clone>(&self, values: &[T]) -> Vec<T> {
        match &self.0 {
            CanonicalStages::Identity { order } => values[..*order].to_vec(),
            CanonicalStages::Mapped { representation, .. } => representation.0.apply_slice(values),
        }
    }

    pub fn representation_to_logical<T: Clone>(&self, values: &[T]) -> Vec<T> {
        match &self.0 {
            CanonicalStages::Identity { order } => values[..*order].to_vec(),
            CanonicalStages::Mapped { representation, .. } => {
                representation.0.apply_slice_inv(values)
            }
        }
    }

    pub fn representation_to_canonical<T: Clone>(&self, values: &[T]) -> Vec<T> {
        match &self.0 {
            CanonicalStages::Identity { order } => values[..*order].to_vec(),
            CanonicalStages::Mapped { index, .. } => index.0.apply_slice(values),
        }
    }

    pub fn canonical_to_representation<T: Clone>(&self, values: &[T]) -> Vec<T> {
        match &self.0 {
            CanonicalStages::Identity { order } => values[..*order].to_vec(),
            CanonicalStages::Mapped { index, .. } => index.0.apply_slice_inv(values),
        }
    }

    pub fn logical_to_canonical<T: Clone>(&self, values: &[T]) -> Vec<T> {
        match &self.0 {
            CanonicalStages::Identity { order } => values[..*order].to_vec(),
            CanonicalStages::Mapped {
                representation,
                index,
            } => index.0.apply_slice(representation.0.apply_slice(values)),
        }
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
    use std::collections::HashSet;
    use std::panic;

    use bincode as binary;
    use bincode::config;
    use ron as ron_codec;
    use serde_json as json;

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
        let representation = [
            String::from("c"),
            String::from("a"),
            String::from("b"),
            String::from("d"),
        ];

        assert_eq!(layout.logical_to_representation(&logical), representation);
        assert_eq!(layout.representation_to_logical(&representation), logical);
        assert_eq!(
            layout.representation_to_canonical(&representation),
            canonical
        );
        assert_eq!(
            layout.canonical_to_representation(&canonical),
            representation
        );
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
            assert!(matches!(layout.0, CanonicalStages::Identity { .. }));
            assert!(matches!(layout.clone().0, CanonicalStages::Identity { .. }));
            let values: Vec<_> = (0..order).map(|position| position.to_string()).collect();
            assert_eq!(layout.logical_to_representation(&values), values);
            assert_eq!(layout.representation_to_logical(&values), values);
            assert_eq!(layout.representation_to_canonical(&values), values);
            assert_eq!(layout.canonical_to_representation(&values), values);
            assert_eq!(layout.canonical_to_logical(&values), values);
            assert_eq!(layout.logical_to_canonical(&values), values);
            for position in 0..order {
                assert_eq!(layout.logical_to_canonical_position(position), position);
                assert_eq!(layout.canonical_to_logical_position(position), position);
            }
        }
    }

    #[test]
    fn identity_stages_normalize_with_consistent_equality_and_hashing() {
        for order in 0..=3 {
            let compact = CanonicalLayout::identity(order);
            let normalized =
                CanonicalLayout::from_permutations(Permutation::id(order), Permutation::id(order));
            assert!(matches!(normalized.0, CanonicalStages::Identity { .. }));
            assert_eq!(compact, normalized);
            assert_eq!(HashSet::from([compact, normalized]).len(), 1);
        }
        assert_ne!(CanonicalLayout::identity(0), CanonicalLayout::identity(1));
    }

    #[test]
    fn cancelling_stages_retain_the_intermediate_representation_order() {
        let layout = CanonicalLayout::from_permutations(
            Permutation::from_map(vec![1, 2, 0]),
            Permutation::from_map(vec![2, 0, 1]),
        );
        let values = ["a", "b", "c"];
        let representation = ["c", "a", "b"];
        assert!(matches!(layout.0, CanonicalStages::Mapped { .. }));
        assert_eq!(layout.logical_to_canonical(&values), values);
        assert_eq!(layout.canonical_to_logical(&values), values);
        assert_eq!(layout.logical_to_representation(&values), representation);
        assert_eq!(layout.canonical_to_representation(&values), representation);
        assert_eq!(layout.representation_to_logical(&representation), values);
        assert_eq!(layout.representation_to_canonical(&representation), values);
        assert_ne!(layout, CanonicalLayout::identity(3));
        assert_eq!(
            HashSet::from([layout, CanonicalLayout::identity(3)]).len(),
            2
        );
    }

    #[test]
    fn either_nonidentity_stage_remains_mapped_and_can_transfer_its_index_order() {
        let representation_only = CanonicalLayout::from_permutations(
            Permutation::from_map(vec![1, 2, 0]),
            Permutation::id(3),
        );
        let index_only = CanonicalLayout::from_permutations(
            Permutation::id(3),
            Permutation::from_map(vec![2, 0, 1]),
        );
        assert!(matches!(
            representation_only.0,
            CanonicalStages::Mapped { .. }
        ));
        assert!(!representation_only.representation_is_identity());
        assert!(representation_only.into_index_permutation().is_identity());
        assert!(matches!(index_only.0, CanonicalStages::Mapped { .. }));
        assert!(index_only.representation_is_identity());
        assert_eq!(index_only.into_index_permutation().map(), &[2, 0, 1]);
        assert_eq!(
            CanonicalLayout::identity(3).into_index_permutation(),
            Permutation::id(3),
        );
    }

    #[test]
    fn layout_codecs_preserve_identity_and_intermediate_stages() {
        for layout in [
            CanonicalLayout::identity(0),
            CanonicalLayout::identity(3),
            CanonicalLayout::from_permutations(
                Permutation::from_map(vec![1, 2, 0]),
                Permutation::from_map(vec![2, 0, 1]),
            ),
            CanonicalLayout::from_permutations(
                Permutation::from_map(vec![1, 2, 0, 3]),
                Permutation::from_map(vec![0, 3, 2, 1]),
            ),
        ] {
            let encoded = json::to_string(&layout).unwrap();
            assert_eq!(json::from_str::<CanonicalLayout>(&encoded).unwrap(), layout);
            let encoded = ron_codec::to_string(&layout).unwrap();
            assert_eq!(
                ron_codec::from_str::<CanonicalLayout>(&encoded).unwrap(),
                layout
            );
            let encoded = binary::encode_to_vec(&layout, config::standard()).unwrap();
            let (decoded, consumed): (CanonicalLayout, _) =
                binary::decode_from_slice(&encoded, config::standard()).unwrap();
            assert_eq!(consumed, encoded.len());
            assert_eq!(decoded, layout);
            let (borrowed, consumed): (CanonicalLayout, _) =
                binary::borrow_decode_from_slice(&encoded, config::standard()).unwrap();
            assert_eq!(consumed, encoded.len());
            assert_eq!(borrowed, layout);
        }
    }

    #[test]
    fn layout_codecs_normalize_identity_stages_and_reject_mismatched_orders() {
        for index_order in [2, 3] {
            let stages = CanonicalStages::Mapped {
                representation: RepresentationOrder(Permutation::id(2)),
                index: IndexOrder(Permutation::id(index_order)),
            };
            let encoded = json::to_string(&stages).unwrap();
            let decoded = json::from_str::<CanonicalLayout>(&encoded);
            let bytes = binary::encode_to_vec(&stages, config::standard()).unwrap();
            let decoded_binary =
                binary::decode_from_slice::<CanonicalLayout, _>(&bytes, config::standard());
            if index_order == 2 {
                let decoded = decoded.unwrap();
                assert!(matches!(decoded.0, CanonicalStages::Identity { order: 2 }));
                assert_eq!(decoded, CanonicalLayout::identity(2));
                assert_eq!(decoded_binary.unwrap().0, decoded);
            } else {
                assert!(decoded.is_err());
                assert!(decoded_binary.is_err());
            }
        }
    }

    #[test]
    fn layout_codecs_reject_malformed_stage_permutations_before_normalizing() {
        let malformed_cases = [
            (vec![1, 0], vec![0, 1]),
            (vec![0, 1], vec![1, 0]),
            (vec![0, 1], vec![0]),
            (vec![0], vec![0, 1]),
            (vec![0, 2], vec![0, 1]),
            (vec![0, 1], vec![0, 2]),
            (vec![0, 0], vec![0, 1]),
            (vec![0, 1], vec![0, 0]),
            (vec![], vec![0]),
            (vec![0], vec![]),
        ];
        for (map, inverse) in malformed_cases {
            // Linnet's raw codec permits these values; CanonicalLayout owns the
            // validation required before stage lookups or identity normalization.
            let malformed: Permutation =
                json::from_str(&format!(r#"{{"map":{map:?},"inv":{inverse:?}}}"#)).unwrap();
            for malformed_representation in [false, true] {
                let valid = Permutation::id(malformed.map().len());
                let (representation, index) = if malformed_representation {
                    (malformed.clone(), valid)
                } else {
                    (valid, malformed.clone())
                };
                let stages = CanonicalStages::Mapped {
                    representation: RepresentationOrder(representation),
                    index: IndexOrder(index),
                };
                let encoded = json::to_string(&stages).unwrap();
                assert!(json::from_str::<CanonicalLayout>(&encoded).is_err());
                let encoded = ron_codec::to_string(&stages).unwrap();
                assert!(ron_codec::from_str::<CanonicalLayout>(&encoded).is_err());
                let bytes = binary::encode_to_vec(&stages, config::standard()).unwrap();
                assert!(
                    binary::decode_from_slice::<CanonicalLayout, _>(&bytes, config::standard())
                        .is_err()
                );
                assert!(
                    binary::borrow_decode_from_slice::<CanonicalLayout, _>(
                        &bytes,
                        config::standard(),
                    )
                    .is_err()
                );
            }
        }
    }

    #[test]
    fn slice_mappings_use_exactly_the_layout_order_and_reject_short_inputs() {
        type SliceMapping = fn(&CanonicalLayout, &[u8]) -> Vec<u8>;
        let mappings: [SliceMapping; 6] = [
            CanonicalLayout::logical_to_representation,
            CanonicalLayout::representation_to_logical,
            CanonicalLayout::representation_to_canonical,
            CanonicalLayout::canonical_to_representation,
            CanonicalLayout::logical_to_canonical,
            CanonicalLayout::canonical_to_logical,
        ];
        for layout in [
            CanonicalLayout::identity(2),
            CanonicalLayout::from_permutations(
                Permutation::from_map(vec![1, 0]),
                Permutation::from_map(vec![1, 0]),
            ),
        ] {
            for mapping in mappings {
                assert_eq!(mapping(&layout, &[0, 1, 2]), mapping(&layout, &[0, 1]));
                assert!(panic::catch_unwind(|| mapping(&layout, &[0])).is_err());
            }
        }
    }

    #[test]
    #[should_panic(expected = "logical axis is out of bounds")]
    fn mapped_logical_position_rejects_out_of_bounds_axis() {
        CanonicalLayout::from_permutations(Permutation::from_map(vec![1, 0]), Permutation::id(2))
            .logical_to_canonical_position(2);
    }

    #[test]
    #[should_panic(expected = "canonical axis is out of bounds")]
    fn mapped_canonical_position_rejects_out_of_bounds_axis() {
        CanonicalLayout::from_permutations(Permutation::id(2), Permutation::from_map(vec![1, 0]))
            .canonical_to_logical_position(2);
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
