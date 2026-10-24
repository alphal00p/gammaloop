//! Persist authoritative tensor metadata and operation flags, never planner caches.
use super::SymbolicTensor;
use super::composition::ExplicitIndexOccurrences;
use super::inference::BorrowedTensorParts;
use crate::CookSettings;
use crate::tensor::TensorProofs;
use bincode::de::Decoder;
use bincode::enc::Encoder;
use bincode::error::{DecodeError, EncodeError};
use bincode::{Decode, Encode};
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{
    PartialIndex, PartialSlot, PartialStructure, PartialStructureExt,
};
use spenso::structure::representation::{LibrarySlot, RepName};
use spenso::structure::slot::IsAbstractSlot;
use std::sync::Arc;
use symbolica::atom::Atom;
use symbolica::state::HasStateMap;

pub(super) const TENSOR_ARCHIVE_VERSION: u8 = 1;

impl Encode for ExplicitIndexOccurrences {
    fn encode<E: Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
        self.validate()
            .map_err(|error| EncodeError::OtherString(error.to_string()))?;
        let mut indices = Vec::new();
        let mut dummies = Vec::new();
        let mut orientations = Vec::new();
        for (slot, total, base, dual) in self.oriented_entries() {
            indices.push(slot);
            dummies.push(total == 2);
            orientations.push((base as u8) | ((dual as u8) << 2));
        }
        indices.encode(encoder)?;
        dummies.encode(encoder)?;
        orientations.encode(encoder)
    }
}

impl<C: HasStateMap> Decode<C> for ExplicitIndexOccurrences {
    fn decode<D: Decoder<Context = C>>(decoder: &mut D) -> Result<Self, DecodeError> {
        let indices = Vec::<LibrarySlot<AbstractIndex>>::decode(decoder)?;
        let dummies = Vec::<bool>::decode(decoder)?;
        let orientations = Vec::<u8>::decode(decoder)?;
        if indices.len() != dummies.len() || indices.len() != orientations.len() {
            return Err(DecodeError::OtherString(
                "serialized tensor index identities, dummy mask and orientation bounds have different lengths".into(),
            ));
        }
        let mut occurrences = Self::default();
        for ((slot, dummy), orientations) in indices.into_iter().zip(dummies).zip(orientations) {
            if occurrences.count(slot) != 0 {
                return Err(DecodeError::OtherString(
                    "serialized tensor index identities must have unique collision keys".into(),
                ));
            }
            let total = 1 + usize::from(dummy);
            let base = usize::from(orientations & 3);
            let dual = usize::from((orientations >> 2) & 3);
            let representative = if !slot.rep().rep.is_base() {
                dual
            } else {
                base
            };
            if orientations & !15 != 0
                || base > 2
                || dual > 2
                || base > total
                || dual > total
                || base + dual < total
                || representative == 0
                || (slot.rep().rep.is_self_dual() && dual != 0)
            {
                return Err(DecodeError::OtherString(
                    "serialized tensor orientation bounds disagree with its index multiplicity or representative".into(),
                ));
            }
            // Imported symbols can have a different ordering. The metadata owner
            // restores sorted identities after applying the decoder's state map.
            occurrences.add_oriented(slot, total, base, dual);
        }
        Ok(occurrences)
    }
}

impl Encode for BorrowedTensorParts<'_> {
    fn encode<E: Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
        let occurrences = self
            .explicit_index_occurrences()
            .map_err(|error| EncodeError::OtherString(error.to_string()))?;
        let _ = self.product_is_normalized();
        TENSOR_ARCHIVE_VERSION.encode(encoder)?;
        self.expression.encode(encoder)?;
        // Canonical storage order alone cannot reconstruct the public axis order.
        self.structure.logical_slots().encode(encoder)?;
        // Operation flags need not match root syntax: negating a metric preserves
        // its metric classification.
        self.is_metric.encode(encoder)?;
        self.is_composite.encode(encoder)?;
        // Readiness controls the next product's notation, so import retains it
        // rather than classifying the normalized payload again.
        self.proofs
            .normalized_product
            .get()
            .copied()
            .encode(encoder)?;
        occurrences.encode(encoder)
    }
}

impl Encode for SymbolicTensor<PartialStructure> {
    fn encode<E: Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
        self.borrowed_parts().encode(encoder)
    }
}

impl<C: HasStateMap> Decode<C> for SymbolicTensor<PartialStructure> {
    fn decode<D: Decoder<Context = C>>(decoder: &mut D) -> Result<Self, DecodeError> {
        if u8::decode(decoder)? != TENSOR_ARCHIVE_VERSION {
            return Err(DecodeError::OtherString(
                "unsupported symbolic tensor archive version".into(),
            ));
        }
        let expression = Atom::decode(decoder)?;
        let slots = Vec::<PartialSlot>::decode(decoder)?;
        let is_metric = bool::decode(decoder)?;
        let is_composite = bool::decode(decoder)?;
        let normalized_product = Option::<bool>::decode(decoder)?;
        let occurrences = ExplicitIndexOccurrences::decode(decoder)?;
        let structure = PartialStructure::from_logical_slots(slots.iter().copied());
        if structure.logical_slots() != slots {
            return Err(DecodeError::OtherString(
                "serialized tensor interface has noncanonical open-port identities".into(),
            ));
        }
        if normalized_product == Some(true) {
            for slot in &slots {
                let open = match slot.aind {
                    PartialIndex::Open(_) => true,
                    PartialIndex::Explicit(index) => CookSettings::index_has_open_port(index)
                        .map_err(|error| DecodeError::OtherString(error.to_string()))?,
                };
                if open {
                    return Err(DecodeError::OtherString(
                        "serialized tensor product readiness cannot include unresolved public ports".into(),
                    ));
                }
            }
        }
        let boundary = ExplicitIndexOccurrences::from_interface(&structure);
        boundary
            .validate()
            .map_err(|error| DecodeError::OtherString(error.to_string()))?;
        if !boundary.is_bounded_by(&occurrences) {
            return Err(DecodeError::OtherString(
                "serialized tensor index metadata does not cover its free interface".into(),
            ));
        }
        // This is an admitted tensor archive, not raw expression/parts intake.
        // Its construction history is authoritative even after cancellation or
        // power coalescence changes the normalized payload's apparent indices.
        // Retained uses and independent variance maxima are construction
        // history, so they need not equal the exposed free-port inventory.
        let proofs = TensorProofs::default();
        let _ = proofs.validated.set(true);
        let _ = proofs.occurrences.set(Arc::new(occurrences));
        if let Some(normalized) = normalized_product {
            let _ = proofs.normalized_product.set(normalized);
        }
        Ok(Self {
            expression,
            structure,
            is_metric,
            is_composite,
            proofs,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::super::inference::tests::{
        INFERENCE_CALLS, OCCURRENCE_FUNCTION_VISITS, SCOPE_VALIDATIONS,
    };
    use super::{ExplicitIndexOccurrences, SymbolicTensor, TENSOR_ARCHIVE_VERSION};
    use crate::representations::ColorFundamental;
    use crate::tensor::composition::TensorCompositionError;
    use crate::test_support;
    use base64::Engine;
    use base64::engine::general_purpose::STANDARD;
    use bincode::config;
    use spenso::structure::abstract_index::AbstractIndex;
    use spenso::structure::partial::{
        OpenPortId, PartialIndex, PartialStructure, PartialStructureExt,
    };
    use spenso::structure::representation::{Euclidean, LibrarySlot, Minkowski, RepName};
    use spenso::structure::slot::IsAbstractSlot;
    use spenso::{bracket, tensor_symbol};
    use std::env;
    use std::io::Cursor;
    use std::process::Command;
    use std::sync::Arc;
    use symbolica::atom::{Atom, AtomCore, Symbol};
    use symbolica::state::{State, StateMap};
    use symbolica::{function, symbol};

    #[test]
    fn tensor_archive_preserves_ranked_zero_and_logical_open_port_order() {
        let structure = PartialStructure::from_logical_slots([
            Minkowski {}
                .new_rep(4)
                .to_lib()
                .slot(PartialIndex::Open(OpenPortId(0))),
            Euclidean {}
                .new_rep(3)
                .to_lib()
                .slot(PartialIndex::Open(OpenPortId(1))),
        ]);
        let source = SymbolicTensor::from_produced_parts(Atom::Zero, structure).unwrap();
        let bytes = bincode::encode_to_vec(&source, config::standard()).unwrap();
        let (decoded, consumed): (SymbolicTensor<PartialStructure>, _) =
            bincode::decode_from_slice_with_context(
                &bytes,
                config::standard(),
                StateMap::default(),
            )
            .unwrap();
        assert_eq!(consumed, bytes.len());
        assert!(decoded.same_value(&source));
        assert_eq!(decoded, source);
        assert_eq!(decoded.rank(), 2);
        assert_eq!(
            decoded.structure.logical_slots(),
            source.structure.logical_slots()
        );
        assert_eq!(decoded.proofs.validated.get(), Some(&true));
        assert_eq!(
            decoded.proofs.normalized_product.get(),
            source.proofs.normalized_product.get()
        );
        assert!(decoded.proofs.observations.get().is_none());
        assert!(decoded.proofs.shallow_graph.get().is_none());
    }

    #[test]
    fn tensor_archive_admits_raw_parts_before_encoding() {
        let structure = PartialStructure::from_logical_slots([Minkowski {}
            .new_rep(4)
            .to_lib()
            .slot(PartialIndex::Open(OpenPortId(0)))]);
        let raw = SymbolicTensor::from_normalized_parts(Atom::one(), structure);
        assert!(bincode::encode_to_vec(raw, config::standard()).is_err());
    }

    #[test]
    fn tensor_archive_rejects_malformed_index_metadata() {
        let slot = Minkowski {}
            .new_rep(4)
            .to_lib()
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(97101));
        let free = vec![
            slot.rep()
                .slot::<PartialIndex<AbstractIndex>, _>(PartialIndex::Explicit(slot.aind)),
        ];
        let dual = ColorFundamental {}
            .new_rep(3)
            .to_lib()
            .dual()
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(97103));
        for (version, slots, indices, dummies, expected) in [
            (
                TENSOR_ARCHIVE_VERSION + 1,
                vec![],
                vec![],
                vec![],
                "version",
            ),
            (
                TENSOR_ARCHIVE_VERSION,
                vec![],
                vec![slot],
                vec![],
                "lengths",
            ),
            (
                TENSOR_ARCHIVE_VERSION,
                vec![],
                vec![slot, slot],
                vec![true, true],
                "unique",
            ),
            (
                TENSOR_ARCHIVE_VERSION,
                free.clone(),
                vec![],
                vec![],
                "cover",
            ),
            (
                TENSOR_ARCHIVE_VERSION,
                vec![free[0], free[0], free[0]],
                vec![slot],
                vec![true],
                "more than two",
            ),
            (
                TENSOR_ARCHIVE_VERSION,
                vec![],
                vec![dual, dual.rep().base().slot(dual.aind)],
                vec![false, true],
                "unique",
            ),
            (
                TENSOR_ARCHIVE_VERSION,
                vec![slot.rep().slot(PartialIndex::Open(OpenPortId(7)))],
                vec![],
                vec![],
                "noncanonical",
            ),
        ] {
            let orientations = indices
                .iter()
                .zip(&dummies)
                .map(|(slot, dummy)| {
                    let total = 1 + u8::from(*dummy);
                    if !slot.rep().rep.is_base() {
                        total << 2
                    } else {
                        total
                    }
                })
                .collect::<Vec<_>>();
            let bytes = bincode::encode_to_vec(
                (
                    version,
                    Atom::Zero,
                    slots,
                    false,
                    true,
                    None::<bool>,
                    indices,
                    dummies,
                    orientations,
                ),
                config::standard(),
            )
            .unwrap();
            let error = bincode::decode_from_slice_with_context::<
                _,
                SymbolicTensor<PartialStructure>,
                _,
            >(&bytes, config::standard(), StateMap::default())
            .unwrap_err()
            .to_string();
            assert!(error.contains(expected), "{error}");
        }
    }

    #[test]
    fn tensor_archive_preserves_independent_orientation_maxima() {
        test_support::test_initialize();
        let rep = ColorFundamental {}.new_rep(3).to_lib();
        let base = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(97151));
        let dual = rep
            .dual()
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(97153));
        let self_dual = Minkowski {}
            .new_rep(4)
            .to_lib()
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(97155));
        let mut source = ExplicitIndexOccurrences::default();
        // Distinct alternatives can have a base maximum of two and a dual
        // maximum of one while no branch exceeds two total occurrences.
        source.add_oriented(base, 2, 2, 1);
        source.add_oriented(dual, 1, 0, 1);
        source.add_oriented(self_dual, 2, 2, 0);
        let bytes = bincode::encode_to_vec(&source, config::standard()).unwrap();
        let (restored, used): (ExplicitIndexOccurrences, _) =
            bincode::decode_from_slice_with_context(
                &bytes,
                config::standard(),
                StateMap::default(),
            )
            .unwrap();
        assert_eq!(used, bytes.len());
        assert_eq!(
            restored.oriented_entries().collect::<Vec<_>>(),
            source.oriented_entries().collect::<Vec<_>>(),
        );
        assert_eq!(restored.count(base), 2);
        assert_eq!(restored.count(dual), 1);
        assert_eq!(restored.count(self_dual), 2);
    }

    #[test]
    fn tensor_archive_rejects_malformed_orientation_bounds() {
        test_support::test_initialize();
        let base = ColorFundamental {}
            .new_rep(3)
            .to_lib()
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(97157));
        let dual = base.rep().dual().slot::<AbstractIndex, _>(base.aind);
        let self_dual = Minkowski {}
            .new_rep(4)
            .to_lib()
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(97159));
        for (slot, dummy, bounds) in [
            (base, false, vec![]),
            (base, false, vec![1u8, 1]),
            (base, false, vec![0]),
            (base, true, vec![3]),
            (base, true, vec![12]),
            (base, false, vec![2]),
            (base, true, vec![1]),
            (base, false, vec![4]),
            (dual, false, vec![1]),
            (self_dual, true, vec![5]),
            (base, false, vec![17]),
        ] {
            let bytes =
                bincode::encode_to_vec((vec![slot], vec![dummy], bounds), config::standard())
                    .unwrap();
            assert!(
                bincode::decode_from_slice_with_context::<_, ExplicitIndexOccurrences, _>(
                    &bytes,
                    config::standard(),
                    StateMap::default(),
                )
                .is_err()
            );
        }
        // The representative can be base-oriented, but its bounds must still
        // cover a declared dual free port with the same collision identity.
        let bytes = bincode::encode_to_vec(
            (
                TENSOR_ARCHIVE_VERSION,
                Atom::Zero,
                vec![
                    dual.rep()
                        .slot::<PartialIndex<AbstractIndex>, _>(PartialIndex::Explicit(dual.aind)),
                ],
                false,
                true,
                Some(true),
                vec![base],
                vec![false],
                vec![1u8],
            ),
            config::standard(),
        )
        .unwrap();
        assert!(
            bincode::decode_from_slice_with_context::<_, SymbolicTensor<PartialStructure>, _>(
                &bytes,
                config::standard(),
                StateMap::default(),
            )
            .is_err()
        );
    }

    #[test]
    fn tensor_archive_retains_product_readiness_without_finishing_the_payload() {
        test_support::test_initialize();
        let slot = Minkowski {}
            .new_rep(4)
            .to_lib()
            .slot(AbstractIndex::Normal(97143));
        let expression = function!(tensor_symbol!("archive_unfinished_product"), slot.to_atom());
        let source = SymbolicTensor::from_declared_parts(
            bracket!(expression),
            PartialStructure::from_logical_slots([slot
                .rep()
                .slot(PartialIndex::Explicit(slot.aind))]),
        )
        .unwrap();
        assert_eq!(source.proofs.normalized_product.get(), Some(&false));
        let bytes = bincode::encode_to_vec(&source, config::standard()).unwrap();
        let (restored, used): (SymbolicTensor<PartialStructure>, _) =
            bincode::decode_from_slice_with_context(
                &bytes,
                config::standard(),
                StateMap::default(),
            )
            .unwrap();
        assert_eq!(used, bytes.len());
        assert!(restored.same_value(&source));
        assert_eq!(restored.proofs.normalized_product.get(), Some(&false));
    }

    #[test]
    fn tensor_archive_rejects_product_readiness_with_unresolved_public_ports() {
        let representation = Minkowski {}.new_rep(4).to_lib();
        for index in [
            PartialIndex::Open(OpenPortId(0)),
            PartialIndex::Explicit(AbstractIndex::Open {
                owner: 97145,
                axis: 0,
            }),
        ] {
            let bytes = bincode::encode_to_vec(
                (
                    TENSOR_ARCHIVE_VERSION,
                    Atom::Zero,
                    vec![representation.slot::<PartialIndex<AbstractIndex>, _>(index)],
                    false,
                    true,
                    Some(true),
                    Vec::<LibrarySlot<AbstractIndex>>::new(),
                    Vec::<bool>::new(),
                    Vec::<u8>::new(),
                ),
                config::standard(),
            )
            .unwrap();
            let error = bincode::decode_from_slice_with_context::<
                _,
                SymbolicTensor<PartialStructure>,
                _,
            >(&bytes, config::standard(), StateMap::default())
            .unwrap_err()
            .to_string();
            assert!(error.contains("readiness"), "{error}");
        }
    }

    #[test]
    fn tensor_archive_preserves_retained_dual_orientation() {
        test_support::test_initialize();
        let slot = ColorFundamental {}
            .new_rep(3)
            .to_lib()
            .dual()
            .slot(AbstractIndex::Normal(97141));
        let source = SymbolicTensor::<PartialStructure>::infer(function!(
            tensor_symbol!("archive_dual_orientation_tensor"),
            slot.to_atom()
        ))
        .unwrap();
        assert_eq!(source.rank(), 1);
        assert_eq!(
            source.structure.logical_slots(),
            [slot.rep().slot(PartialIndex::Explicit(slot.aind))],
        );
        assert_eq!(
            source
                .explicit_index_occurrences()
                .unwrap()
                .entries()
                .collect::<Vec<_>>(),
            [(slot, 1)],
        );
        let bytes = bincode::encode_to_vec(&source, config::standard()).unwrap();
        let (restored, used): (SymbolicTensor<PartialStructure>, _) =
            bincode::decode_from_slice_with_context(
                &bytes,
                config::standard(),
                StateMap::default(),
            )
            .unwrap();
        assert_eq!(used, bytes.len());
        assert!(restored.same_value(&source));
        let metadata = restored.explicit_index_occurrences().unwrap();
        assert_eq!(metadata.entries().collect::<Vec<_>>(), [(slot, 1)]);
        assert_eq!(metadata.count(slot), 1);
        assert_eq!(metadata.count(slot.rep().base().slot(slot.aind)), 1);
    }

    #[test]
    fn tensor_archive_preserves_scalar_singletons_and_conservative_free_port_counts() {
        test_support::test_initialize();
        let slot = Minkowski {}
            .new_rep(4)
            .to_lib()
            .slot(AbstractIndex::Normal(97131));
        let scalar = SymbolicTensor::<PartialStructure>::infer(slot.to_atom()).unwrap();
        assert!(scalar.is_scalar());
        assert_eq!(scalar.explicit_index_occurrences().unwrap().count(slot), 1);
        let vector = SymbolicTensor::<PartialStructure>::infer(function!(
            tensor_symbol!("archive_scalar_singleton_vector"),
            slot.to_atom()
        ))
        .unwrap();
        let product = scalar.multiply(&vector).unwrap();
        assert_eq!(product.rank(), 1);
        assert_eq!(product.explicit_index_occurrences().unwrap().count(slot), 2);
        for source in [scalar, product] {
            let expected = source.explicit_index_occurrences().unwrap().count(slot);
            INFERENCE_CALLS.with(|count| count.set(0));
            OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
            SCOPE_VALIDATIONS.with(|count| count.set(0));
            let bytes = bincode::encode_to_vec(&source, config::standard()).unwrap();
            let (restored, used): (SymbolicTensor<PartialStructure>, _) =
                bincode::decode_from_slice_with_context(
                    &bytes,
                    config::standard(),
                    StateMap::default(),
                )
                .unwrap();
            assert_eq!(used, bytes.len());
            assert!(restored.same_value(&source));
            assert_eq!(
                restored.explicit_index_occurrences().unwrap().count(slot),
                expected
            );
            assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
            assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()), 0);
            assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
        }
    }

    #[test]
    fn tensor_archive_and_clones_retain_cancelled_dummy_collisions_without_observation() {
        test_support::test_initialize();
        let slot = Minkowski {}
            .new_rep(4)
            .to_lib()
            .slot(AbstractIndex::Normal(97111));
        let pair = function!(tensor_symbol!("archive_cancel_left"), slot.to_atom())
            * function!(tensor_symbol!("archive_cancel_right"), slot.to_atom());
        let x = Atom::var(symbol!("archive_cancel_x"));
        let polynomial = (&x + 1).pow(2) - &x * &x - Atom::num(2) * &x - 1;
        let source = SymbolicTensor::<PartialStructure>::infer(pair * polynomial)
            .unwrap()
            .expanded(None, false)
            .unwrap();
        assert!(source.is_zero());
        assert_eq!(source.explicit_index_occurrences().unwrap().count(slot), 2);
        let extra = SymbolicTensor::<PartialStructure>::infer(function!(
            tensor_symbol!("archive_cancel_extra"),
            slot.to_atom()
        ))
        .unwrap();
        INFERENCE_CALLS.with(|count| count.set(0));
        OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
        SCOPE_VALIDATIONS.with(|count| count.set(0));
        let bytes = bincode::encode_to_vec(&source, config::standard()).unwrap();
        let (restored, used): (SymbolicTensor<PartialStructure>, _) =
            bincode::decode_from_slice_with_context(
                &bytes,
                config::standard(),
                StateMap::default(),
            )
            .unwrap();
        assert_eq!(used, bytes.len());
        assert!(restored.same_value(&source));
        assert_eq!(restored.proofs.normalized_product.get(), Some(&true));
        let cloned = restored.clone();
        assert!(Arc::ptr_eq(
            restored.proofs.occurrences.get().unwrap(),
            cloned.proofs.occurrences.get().unwrap(),
        ));
        for value in [source, restored, cloned] {
            assert!(matches!(
                value.multiply(&extra),
                Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
            ));
        }
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()), 0);
        assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
    }

    #[test]
    fn tensor_archive_preserves_operation_flags_after_negating_a_metric() {
        crate::test_support::test_initialize();
        let metric = SymbolicTensor::infer(function!(
            spenso::network::library::symbolic::ETS.metric,
            Minkowski {}.new_rep(4).to_symbolic([Atom::num(91)]),
            Minkowski {}.new_rep(4).to_symbolic([Atom::num(93)])
        ))
        .unwrap();
        let source = -metric;
        assert!(source.is_metric());
        assert!(source.is_composite());
        let bytes = bincode::encode_to_vec(&source, bincode::config::standard()).unwrap();
        let (decoded, consumed): (SymbolicTensor<PartialStructure>, _) =
            bincode::decode_from_slice_with_context(
                &bytes,
                bincode::config::standard(),
                StateMap::default(),
            )
            .unwrap();
        assert_eq!(consumed, bytes.len());
        assert_eq!(decoded, source);
    }

    fn indexed_tensor(dimension: Symbol, index: Symbol) -> SymbolicTensor<PartialStructure> {
        SymbolicTensor::infer(function!(
            tensor_symbol!("idenso::persistence_test::vector"),
            Minkowski {}
                .new_rep(dimension)
                .to_symbolic([Atom::var(index)])
        ))
        .unwrap()
    }

    #[test]
    fn tensor_archive_remaps_symbols_in_both_payload_and_interface() {
        const CHILD: &str = "IDENSO_TENSOR_ARCHIVE_CHILD";
        if env::var_os(CHILD).is_some() {
            let dimension = symbol!("idenso::persistence_test::dimension");
            let index = symbol!("idenso::persistence_test::index");
            let source = indexed_tensor(dimension, index);
            let other_index = symbol!("idenso::persistence_test::other_index");
            let other = indexed_tensor(dimension, other_index);
            let bytes = bincode::encode_to_vec(&source, config::standard()).unwrap();
            let x = Atom::var(symbol!("idenso::persistence_test::cancel_x"));
            let polynomial = (&x + 1).pow(2) - &x * &x - Atom::num(2) * &x - 1;
            let coefficient = SymbolicTensor::<PartialStructure>::infer(polynomial).unwrap();
            let hidden = source
                .multiply(&source)
                .unwrap()
                .multiply(&other.multiply(&other).unwrap())
                .unwrap()
                .multiply(&coefficient)
                .unwrap()
                .expanded(None, false)
                .unwrap();
            assert!(hidden.is_zero());
            assert!(hidden.is_scalar());
            let hidden = bincode::encode_to_vec(hidden, config::standard()).unwrap();
            let mut state = Vec::new();
            State::export(&mut state).unwrap();
            println!("TENSOR_INDEX={}", index.get_id());
            println!("TENSOR_STATE={}", STANDARD.encode(state));
            println!("TENSOR_BYTES={}", STANDARD.encode(bytes));
            println!("TENSOR_HIDDEN_BYTES={}", STANDARD.encode(hidden));
            return;
        }

        let output = Command::new(env::current_exe().unwrap())
            .args([
                "--exact",
                "tensor::persistence::tests::tensor_archive_remaps_symbols_in_both_payload_and_interface",
                "--nocapture",
            ])
            .env(CHILD, "1")
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let output = String::from_utf8(output.stdout).unwrap();
        let field = |prefix| {
            output
                .lines()
                .find_map(|line| line.strip_prefix(prefix))
                .unwrap()
        };
        let old_index: u32 = field("TENSOR_INDEX=").parse().unwrap();
        // Register unrelated names first so remapping is exercised even when
        // this is the only selected test in the parent process.
        for position in 0..=old_index {
            let _ = symbol!(&format!("idenso::persistence_test::padding_{position}"));
        }
        // Reverse the registration order of the two hidden metadata identities.
        // Their imported ordering must be restored before binary-search lookup.
        let other_index = symbol!("idenso::persistence_test::other_index");
        let index = symbol!("idenso::persistence_test::index");
        let dimension = symbol!("idenso::persistence_test::dimension");
        assert_ne!(index.get_id(), old_index);
        let expected = indexed_tensor(dimension, index);
        let state = STANDARD.decode(field("TENSOR_STATE=")).unwrap();
        let map = State::import(&mut Cursor::new(&state), None).unwrap();
        let bytes = STANDARD.decode(field("TENSOR_BYTES=")).unwrap();
        let (decoded, consumed): (SymbolicTensor<PartialStructure>, _) =
            bincode::decode_from_slice_with_context(&bytes, config::standard(), map).unwrap();
        assert_eq!(consumed, bytes.len());
        assert!(decoded.same_value(&expected));
        assert_eq!(decoded, expected);
        assert_eq!(
            decoded.structure.logical_slots(),
            expected.structure.logical_slots()
        );
        let map = State::import(&mut Cursor::new(state), None).unwrap();
        let bytes = STANDARD.decode(field("TENSOR_HIDDEN_BYTES=")).unwrap();
        let (hidden, used): (SymbolicTensor<PartialStructure>, _) =
            bincode::decode_from_slice_with_context(&bytes, config::standard(), map).unwrap();
        assert_eq!(used, bytes.len());
        assert!(hidden.is_zero());
        assert!(hidden.is_scalar());
        let index = AbstractIndex::Symbol(index.into());
        let slot = Minkowski {}.new_rep(dimension).to_lib().slot(index);
        assert_eq!(hidden.explicit_index_occurrences().unwrap().count(slot), 2);
        let other_slot = slot.rep().slot(AbstractIndex::Symbol(other_index.into()));
        assert_eq!(
            hidden
                .explicit_index_occurrences()
                .unwrap()
                .count(other_slot),
            2
        );
        assert!(hidden.multiply(&expected).is_err());
    }
}
