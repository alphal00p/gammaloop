//! Persist tensor data and operation flags, never proof or planner caches.

use bincode::{Decode, Encode};
use spenso::structure::partial::{PartialSlot, PartialStructure, PartialStructureExt};
use symbolica::{atom::Atom, state::HasStateMap};

use super::SymbolicTensor;

impl Encode for SymbolicTensor<PartialStructure> {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        self.expression.encode(encoder)?;
        // Canonical storage order alone cannot reconstruct the public axis order.
        self.structure.logical_slots().encode(encoder)?;
        // These operation flags participate in value equality and need not match
        // the root syntax: negating a metric preserves its metric classification.
        self.is_metric.encode(encoder)?;
        self.is_composite.encode(encoder)
    }
}

impl<C: HasStateMap> Decode<C> for SymbolicTensor<PartialStructure> {
    fn decode<D: bincode::de::Decoder<Context = C>>(
        decoder: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        let expression = Atom::decode(decoder)?;
        let slots = Vec::<PartialSlot>::decode(decoder)?;
        let is_metric = bool::decode(decoder)?;
        let is_composite = bool::decode(decoder)?;
        // Atom, dimension and index decoding all use the same symbol-state map.
        // Imported parts establish a fresh certificate only after their payload,
        // interface, placeholder scopes and index multiplicities agree.
        let mut value =
            Self::from_declared_parts(expression, PartialStructure::from_logical_slots(slots))
                .map_err(|error| {
                    bincode::error::DecodeError::OtherString(format!(
                        "serialized symbolic tensor has an invalid payload or interface: {error}"
                    ))
                })?;
        value.is_metric = is_metric;
        value.is_composite = is_composite;
        Ok(value)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use base64::{Engine, engine::general_purpose::STANDARD};
    use spenso::{
        structure::{
            partial::{OpenPortId, PartialIndex, PartialStructure, PartialStructureExt},
            representation::{Euclidean, Minkowski, RepName},
        },
        tensor_symbol,
    };
    use symbolica::{
        atom::{Atom, Symbol},
        function,
        state::{State, StateMap},
        symbol,
    };

    use super::SymbolicTensor;

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
        let source = SymbolicTensor::checked_parts(Atom::Zero, structure).unwrap();
        let bytes = bincode::encode_to_vec(&source, bincode::config::standard()).unwrap();
        let (decoded, consumed): (SymbolicTensor<PartialStructure>, _) =
            bincode::decode_from_slice_with_context(
                &bytes,
                bincode::config::standard(),
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
        assert!(decoded.proofs.observations.get().is_none());
        assert!(decoded.proofs.shallow_graph.get().is_none());
    }

    #[test]
    fn tensor_archive_rejects_a_nonzero_payload_with_an_unrelated_interface() {
        let slots = PartialStructure::from_logical_slots([Minkowski {}
            .new_rep(4)
            .to_lib()
            .slot(PartialIndex::Open(OpenPortId(0)))])
        .logical_slots();
        let bytes = bincode::encode_to_vec(
            (Atom::one(), slots, false, false),
            bincode::config::standard(),
        )
        .unwrap();
        assert!(
            bincode::decode_from_slice_with_context::<_, SymbolicTensor<PartialStructure>, _>(
                &bytes,
                bincode::config::standard(),
                StateMap::default(),
            )
            .is_err()
        );
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
        if std::env::var_os(CHILD).is_some() {
            let dimension = symbol!("idenso::persistence_test::dimension");
            let index = symbol!("idenso::persistence_test::index");
            let source = indexed_tensor(dimension, index);
            let bytes = bincode::encode_to_vec(source, bincode::config::standard()).unwrap();
            let mut state = Vec::new();
            State::export(&mut state).unwrap();
            println!("TENSOR_INDEX={}", index.get_id());
            println!("TENSOR_STATE={}", STANDARD.encode(state));
            println!("TENSOR_BYTES={}", STANDARD.encode(bytes));
            return;
        }

        let output = std::process::Command::new(std::env::current_exe().unwrap())
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
        let index = symbol!("idenso::persistence_test::index");
        let dimension = symbol!("idenso::persistence_test::dimension");
        assert_ne!(index.get_id(), old_index);
        let expected = indexed_tensor(dimension, index);
        let state = STANDARD.decode(field("TENSOR_STATE=")).unwrap();
        let map = State::import(&mut Cursor::new(state), None).unwrap();
        let bytes = STANDARD.decode(field("TENSOR_BYTES=")).unwrap();
        let (decoded, consumed): (SymbolicTensor<PartialStructure>, _) =
            bincode::decode_from_slice_with_context(&bytes, bincode::config::standard(), map)
                .unwrap();
        assert_eq!(consumed, bytes.len());
        assert!(decoded.same_value(&expected));
        assert_eq!(decoded, expected);
        assert_eq!(
            decoded.structure.logical_slots(),
            expected.structure.logical_slots()
        );
    }
}
