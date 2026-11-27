//! Persist authoritative tensor metadata and operation flags, never planner caches.
//!
//! These archives are trusted GammaLoop state files. Decoding trusts stored
//! metadata, restores stored regional facts without inference, and delegates
//! payload decoding to Symbolica, whose Atom decoder is
//! not safe for untrusted bytes. The fixed magic and version reject unsupported
//! formats before the Atom decoder can interpret their bytes as an allocation length.
use super::SymbolicTensor;
use super::composition::ExplicitIndexOccurrences;
use super::inference::BorrowedTensorParts;
use super::regions::AdmittedInterfaces;
use crate::CookSettings;
use crate::tensor::TensorProofs;
use ahash::{AHashMap, AHashSet};
use bincode::config;
use bincode::de::Decoder;
use bincode::enc::Encoder;
use bincode::error::{AllowedEnumVariants, DecodeError, EncodeError};
use bincode::{Decode, Encode};
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{
    PartialIndex, PartialSlot, PartialStructure, PartialStructureExt,
};
use spenso::structure::representation::{LibrarySlot, RepName};
use spenso::structure::slot::IsAbstractSlot;
use std::collections::HashMap;
use std::io::{Cursor, Error as IoError, ErrorKind, Read, Result as IoResult, Write};
use std::sync::{Arc, Mutex};
use symbolica::atom::{Atom, AtomView};
use symbolica::state::HasStateMap;

pub(super) const TENSOR_ARCHIVE_MAGIC: [u8; 8] = [0xd9, b'I', b'D', b'E', b'N', b'S', b'O', 0];
pub(super) const TENSOR_ARCHIVE_VERSION: u8 = 2;

/// Public Symbolica stream fragments. Payload bytes occur once in the archive;
/// every fact refers to the exact borrowed subtree bytes written by Symbolica.
/// Its framing is captured verbatim instead of reproducing its encoding format.
#[derive(Debug, Encode, Decode)]
enum AtomFragment {
    Bytes(Vec<u8>),
    Payload { offset: usize, length: usize },
}

struct AtomFragmentWriter<'a> {
    payload: &'a [u8],
    prefix_length: usize,
    fragments: Vec<AtomFragment>,
}

impl Write for AtomFragmentWriter<'_> {
    fn write(&mut self, bytes: &[u8]) -> IoResult<usize> {
        let offset = (bytes.as_ptr() as usize).checked_sub(self.payload.as_ptr() as usize);
        if let Some(offset) = offset
            && offset
                .checked_add(bytes.len())
                .is_some_and(|end| end <= self.payload.len())
        {
            self.fragments.push(AtomFragment::Payload {
                offset: self.prefix_length + offset,
                length: bytes.len(),
            });
        } else {
            self.fragments.push(AtomFragment::Bytes(bytes.to_vec()));
        }
        Ok(bytes.len())
    }

    fn flush(&mut self) -> IoResult<()> {
        Ok(())
    }
}

struct AtomFragmentReader<'a> {
    payload: &'a [u8],
    fragments: &'a [AtomFragment],
    fragment: usize,
    offset: usize,
}

impl Read for AtomFragmentReader<'_> {
    fn read(&mut self, output: &mut [u8]) -> IoResult<usize> {
        while let Some(fragment) = self.fragments.get(self.fragment) {
            let bytes = match fragment {
                AtomFragment::Bytes(bytes) => bytes.as_slice(),
                AtomFragment::Payload { offset, length } => offset
                    .checked_add(*length)
                    .and_then(|end| self.payload.get(*offset..end))
                    .ok_or_else(|| {
                        IoError::new(
                            ErrorKind::InvalidData,
                            "tensor fact payload range is outside its atom",
                        )
                    })?,
            };
            if self.offset == bytes.len() {
                self.fragment += 1;
                self.offset = 0;
                continue;
            }
            let count = output.len().min(bytes.len() - self.offset);
            output[..count].copy_from_slice(&bytes[self.offset..self.offset + count]);
            self.offset += count;
            return Ok(count);
        }
        Ok(0)
    }
}

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
        TENSOR_ARCHIVE_MAGIC.encode(encoder)?;
        TENSOR_ARCHIVE_VERSION.encode(encoder)?;
        let payload = bincode::encode_to_vec(self.expression, config::standard())?;
        let raw = self.expression.as_view().get_data();
        let prefix_length = payload
            .len()
            .checked_sub(raw.len())
            .ok_or(EncodeError::Other(
                "Symbolica atom stream does not contain its borrowed payload",
            ))?;
        if payload[prefix_length..] != *raw {
            return Err(EncodeError::Other(
                "Symbolica atom stream does not end with its borrowed payload",
            ));
        }
        payload.encode(encoder)?;
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
        occurrences.encode(encoder)?;
        let read_leaf = self.proofs.leaf_interfaces.reader();
        let mut pending = vec![self.expression.as_view()];
        let mut seen = AHashSet::new();
        let mut records = Vec::new();
        while let Some(value) = pending.pop() {
            if !seen.insert(value.get_data()) {
                continue;
            }
            match value {
                AtomView::Fun(function) => pending.extend(function.iter()),
                AtomView::Add(sum) => pending.extend(sum.iter()),
                AtomView::Mul(product) => pending.extend(product.iter()),
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    pending.extend([base, exponent]);
                }
                AtomView::Num(_) | AtomView::Var(_) => {}
            }
            let interface = self.proofs.retained_interface(value);
            let leaf = read_leaf(value);
            if interface.is_none() && leaf.is_none() {
                continue;
            }
            let mut writer = AtomFragmentWriter {
                payload: raw,
                prefix_length,
                fragments: Vec::new(),
            };
            value
                .write(&mut writer)
                .map_err(|error| EncodeError::OtherString(error.to_string()))?;
            records.push((
                writer.fragments,
                interface.map(|value| value.logical_slots()),
                leaf,
            ));
        }
        records.encode(encoder)
    }
}

impl Encode for SymbolicTensor<PartialStructure> {
    fn encode<E: Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
        self.borrowed_parts().encode(encoder)
    }
}

impl<C: HasStateMap> Decode<C> for SymbolicTensor<PartialStructure> {
    fn decode<D: Decoder<Context = C>>(decoder: &mut D) -> Result<Self, DecodeError> {
        if <[u8; 8]>::decode(decoder)? != TENSOR_ARCHIVE_MAGIC {
            return Err(DecodeError::Other(
                "unrecognized symbolic tensor archive magic; legacy archives are unsupported",
            ));
        }
        let version = u8::decode(decoder)?;
        if version != TENSOR_ARCHIVE_VERSION {
            return Err(DecodeError::UnexpectedVariant {
                type_name: "symbolic tensor archive version",
                allowed: &AllowedEnumVariants::Allowed(&[2]),
                found: u32::from(version),
            });
        }
        let payload = Vec::<u8>::decode(decoder)?;
        let mut reader = Cursor::new(&payload);
        let expression = Atom::import_with_map(&mut reader, decoder.context().get_state_map())
            .map_err(|error| DecodeError::OtherString(error.to_string()))?;
        if reader.position() != payload.len() as u64 {
            return Err(DecodeError::Other("trailing bytes in tensor payload"));
        }
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
        let mut proofs = TensorProofs::default();
        let _ = proofs.occurrences.set(Arc::new(occurrences));
        if let Some(normalized) = normalized_product {
            let _ = proofs.normalized_product.set(normalized);
        }
        let records = Vec::<(
            Vec<AtomFragment>,
            Option<Vec<PartialSlot>>,
            Option<(bool, Vec<Atom>)>,
        )>::decode(decoder)?;
        // Same-state decoding uses verified subtree ranges directly. If symbol
        // remapping changes layout or normalization, import each stored key
        // independently, so its facts follow the correct normalized subtree.
        // Storage is linear in payload plus fact records. Remapped decoding can
        // visit nested subtree bytes repeatedly; it never infers their interface.
        let raw = expression.as_view().get_data();
        let prefix_length = payload.len().checked_sub(raw.len());
        let mut borrowed = AHashMap::new();
        if let Some(prefix_length) = prefix_length
            && payload[prefix_length..] == *raw
        {
            let mut pending = vec![expression.as_view()];
            while let Some(value) = pending.pop() {
                let bytes = value.get_data();
                if let Some(offset) = (bytes.as_ptr() as usize).checked_sub(raw.as_ptr() as usize)
                    && offset
                        .checked_add(bytes.len())
                        .is_some_and(|end| end <= raw.len())
                {
                    borrowed.insert((prefix_length + offset, bytes.len()), value);
                }
                match value {
                    AtomView::Fun(function) => pending.extend(function.iter()),
                    AtomView::Add(sum) => pending.extend(sum.iter()),
                    AtomView::Mul(product) => pending.extend(product.iter()),
                    AtomView::Pow(power) => {
                        let (base, exponent) = power.get_base_exp();
                        pending.extend([base, exponent]);
                    }
                    AtomView::Num(_) | AtomView::Var(_) => {}
                }
            }
        }
        let mut interfaces = AdmittedInterfaces::for_expression(expression.as_view());
        let mut leaves = HashMap::new();
        for (fragments, interface, leaf) in records {
            let cached = match fragments.as_slice() {
                [prefix @ .., AtomFragment::Payload { offset, length }]
                    if prefix
                        .iter()
                        .all(|part| matches!(part, AtomFragment::Bytes(_))) =>
                {
                    borrowed.get(&(*offset, *length)).copied()
                }
                _ => None,
            };
            let imported;
            let value = if let Some(value) = cached {
                value
            } else {
                let mut reader = AtomFragmentReader {
                    payload: &payload,
                    fragments: &fragments,
                    fragment: 0,
                    offset: 0,
                };
                imported = Atom::import_with_map(&mut reader, decoder.context().get_state_map())
                    .map_err(|error| DecodeError::OtherString(error.to_string()))?;
                let mut extra = [0];
                if reader
                    .read(&mut extra)
                    .map_err(|error| DecodeError::OtherString(error.to_string()))?
                    != 0
                {
                    return Err(DecodeError::Other("trailing bytes in tensor fact key"));
                }
                imported.as_view()
            };
            if let Some(slots) = interface {
                let interface = PartialStructure::from_logical_slots(slots.iter().copied());
                if interface.logical_slots() != slots {
                    return Err(DecodeError::Other(
                        "serialized regional interface has noncanonical open-port identities",
                    ));
                }
                interfaces.insert(value, interface);
            }
            if let Some(leaf) = leaf {
                leaves.insert(value.to_owned(), leaf);
            }
        }
        interfaces.insert(expression.as_view(), structure.clone());
        let _ = proofs.admitted_interfaces.set(Arc::new(interfaces));
        proofs.leaf_interfaces = Arc::new(Mutex::new(leaves)).into();
        Ok(Self {
            expression,
            structure,
            is_metric,
            is_composite,
            work: Default::default(),
            proofs,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::super::inference::tests::{
        INFERENCE_CALLS, INFERENCE_VISITS, OCCURRENCE_FUNCTION_VISITS, SCOPE_VALIDATIONS,
    };
    use super::{
        ExplicitIndexOccurrences, SymbolicTensor, TENSOR_ARCHIVE_MAGIC, TENSOR_ARCHIVE_VERSION,
    };
    use crate::representations::ColorFundamental;
    use crate::tensor::composition::{MatrixChannel, TensorCompositionError};
    use crate::tensor::network_metadata_tests::RAW_PARSER_ADMISSIONS;
    use crate::tensor::{SymbolicNetExt, SymbolicNetParse};
    use crate::test_support;
    use base64::Engine;
    use base64::engine::general_purpose::STANDARD;
    use bincode::config;
    use bincode::error::DecodeError;
    use spenso::network::library::symbolic::ETS;
    use spenso::network::parsing::ParseSettings;
    use spenso::structure::abstract_index::AbstractIndex;
    use spenso::structure::partial::{
        OpenPortId, PartialIndex, PartialStructure, PartialStructureExt,
    };
    use spenso::structure::representation::{
        Euclidean, LibraryRep, LibrarySlot, Minkowski, RepName, Representation,
    };
    use spenso::structure::slot::IsAbstractSlot;
    use spenso::{bracket, tensor_symbol};
    use std::env;
    use std::io::Cursor;
    use std::process::Command;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
    use symbolica::state::{State, StateMap};
    use symbolica::{function, symbol};

    #[test]
    fn tensor_archive_shares_nested_fact_payload_bytes() {
        test_support::test_initialize();
        let wrapper = symbol!("archive_scalar_wrapper");
        let mut expression = Atom::var(symbol!("archive_scalar_argument"));
        for _ in 0..256 {
            expression = function!(wrapper, expression);
        }
        let source = SymbolicTensor::infer(expression).unwrap();
        let mut pending = vec![source.expression().as_view()];
        while let Some(value) = pending.pop() {
            source
                .proofs
                .retain_interface(value, PartialStructure::from_logical_slots([]));
            match value {
                AtomView::Fun(function) => pending.extend(function.iter()),
                AtomView::Add(sum) => pending.extend(sum.iter()),
                AtomView::Mul(product) => pending.extend(product.iter()),
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    pending.extend([base, exponent]);
                }
                AtomView::Num(_) | AtomView::Var(_) => {}
            }
        }
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let bytes = bincode::encode_to_vec(&source, config::standard()).unwrap();
        assert!(
            bytes.len() < 32 * source.expression().as_view().get_byte_size(),
            "nested facts must refer to shared payload bytes, not serialize every full subtree",
        );
        let (restored, used): (SymbolicTensor<PartialStructure>, _) =
            bincode::decode_from_slice_with_context(
                &bytes,
                config::standard(),
                StateMap::default(),
            )
            .unwrap();
        assert_eq!(used, bytes.len());
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        let mut pending = vec![restored.expression().as_view()];
        while let Some(value) = pending.pop() {
            assert_eq!(
                restored.proofs.retained_interface(value),
                Some(PartialStructure::from_logical_slots([])),
            );
            match value {
                AtomView::Fun(function) => pending.extend(function.iter()),
                AtomView::Add(sum) => pending.extend(sum.iter()),
                AtomView::Mul(product) => pending.extend(product.iter()),
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    pending.extend([base, exponent]);
                }
                AtomView::Num(_) | AtomView::Var(_) => {}
            }
        }
    }

    #[test]
    fn tensor_archive_preserves_republished_subtree_facts() {
        test_support::test_initialize();
        let rep = Minkowski {}.new_rep(4).to_lib();
        let ports = [97215, 97211].map(|index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let leaf = function!(tensor_symbol!("archive_kept_leaf"), &ports[0], &ports[1]);
        let sibling = function!(tensor_symbol!("archive_omitted_leaf"), &ports[0], &ports[1]);
        let term = SymbolicTensor::infer(Atom::num(2) * &leaf).unwrap();
        let other = SymbolicTensor::infer(sibling.clone()).unwrap();
        let leaf_interface = term.proofs.retained_interface(leaf.as_view()).unwrap();
        let sum = term.try_add(&other).unwrap();
        assert_eq!(
            sum.proofs.retained_interface(term.expression().as_view()),
            Some(term.structure.clone())
        );
        assert!(sum.proofs.retained_interface(sibling.as_view()).is_some());
        assert!(
            sum.admitted_leaf_interfaces()
                .lock()
                .unwrap()
                .contains_key(&sibling)
        );

        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        RAW_PARSER_ADMISSIONS.with(|count| count.set(0));
        let published = sum
            .with_produced_expression(term.expression().clone())
            .unwrap();
        let bytes = bincode::encode_to_vec(&published, config::standard()).unwrap();
        let (restored, used): (SymbolicTensor<PartialStructure>, _) =
            bincode::decode_from_slice_with_context(
                &bytes,
                config::standard(),
                StateMap::default(),
            )
            .unwrap();
        assert_eq!(used, bytes.len());
        for value in [&published, &restored] {
            assert_eq!(value.expression(), term.expression());
            assert_eq!(value.structure, sum.structure);
            assert_eq!(
                value.proofs.retained_interface(term.expression().as_view()),
                Some(sum.structure.clone())
            );
            assert_eq!(
                value.proofs.retained_interface(leaf.as_view()),
                Some(leaf_interface.clone())
            );
            assert!(value.proofs.retained_interface(sibling.as_view()).is_none());
            let leaves = value.admitted_leaf_interfaces();
            let leaves = leaves.lock().unwrap();
            let (tensorial, written_ports) = &leaves[&leaf];
            assert!(*tensorial);
            assert_eq!(written_ports.as_slice(), ports.as_slice());
            assert!(!leaves.contains_key(&sibling));
        }
        let network = restored
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(network.simple_execute::<()>().unwrap(), *term.expression());
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        assert_eq!(RAW_PARSER_ADMISSIONS.with(|count| count.get()), 0);
    }

    #[test]
    fn decoded_nested_regions_remain_contractible_after_scalar_rewrites() {
        test_support::test_initialize();
        let representation = Minkowski {}.new_rep(4).to_lib();
        let [free, left, right] = [97201, 97203, 97205].map(|index| {
            representation
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let open = tensor_symbol!(
            "archive_nested_open",
            norm = move |_node, _out| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let pair = function!(ETS.metric, &left, &right)
            * function!(tensor_symbol!("archive_nested_left"), &left)
            * function!(tensor_symbol!("archive_nested_right"), &right);
        let weight = Atom::var(symbol!("archive_nested_weight"));
        for nested in [pair.clone(), (&pair + &weight).pow(2)] {
            let source = SymbolicTensor::infer(function!(open, &free) * nested).unwrap();
            let rewritten = Atom::num(2) * source.expression();
            let expected = source
                .with_produced_expression(rewritten.clone())
                .unwrap()
                .contract(Default::default())
                .unwrap();
            let bytes = bincode::encode_to_vec(&source, config::standard()).unwrap();
            calls.store(0, Ordering::Relaxed);
            INFERENCE_CALLS.with(|count| count.set(0));
            INFERENCE_VISITS.with(|count| count.set(0));
            let (restored, used): (SymbolicTensor<PartialStructure>, _) =
                bincode::decode_from_slice_with_context(
                    &bytes,
                    config::standard(),
                    StateMap::default(),
                )
                .unwrap();
            assert_eq!(used, bytes.len());
            assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
            assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
            assert_eq!(
                calls.load(Ordering::Relaxed),
                0,
                "decode must not replay normalizers"
            );
            assert_eq!(restored.structure, source.structure);
            assert_eq!(
                restored
                    .explicit_index_occurrences()
                    .unwrap()
                    .oriented_entries()
                    .collect::<Vec<_>>(),
                source
                    .explicit_index_occurrences()
                    .unwrap()
                    .oriented_entries()
                    .collect::<Vec<_>>()
            );
            INFERENCE_CALLS.with(|count| count.set(0));
            INFERENCE_VISITS.with(|count| count.set(0));
            let actual = restored
                .with_produced_expression(rewritten)
                .unwrap()
                .contract(Default::default())
                .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
            assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        }
    }

    #[test]
    fn tensor_archive_preserves_detached_zero_channel_inference_convention() {
        test_support::test_initialize();
        let minkowski = Minkowski {}.new_rep(4).to_lib();
        let euclidean = Euclidean {}.new_rep(4).to_lib();
        let fundamental = ColorFundamental {}.new_rep(3).to_lib();
        let port = |representation: Representation<LibraryRep>, index| {
            representation.slot(PartialIndex::Explicit(AbstractIndex::Normal(index)))
        };
        let cases = [
            (
                SymbolicTensor::composite_interface([port(minkowski, 90), port(minkowski, 30)]),
                Some(MatrixChannel {
                    input: 0,
                    output: 1,
                }),
            ),
            (
                PartialStructure::from_logical_slots([
                    port(fundamental.dual(), 30),
                    port(fundamental, 90),
                ]),
                Some(MatrixChannel {
                    input: 1,
                    output: 0,
                }),
            ),
            (
                SymbolicTensor::composite_interface([
                    port(minkowski, 90),
                    port(minkowski, 30),
                    port(euclidean, 9),
                    port(euclidean, 8),
                ]),
                None,
            ),
        ];
        for (interface, expected) in cases {
            let source = SymbolicTensor::from_declared_parts(Atom::Zero, interface).unwrap();
            assert_eq!(source.matrix_channel(), expected);
            let bytes = bincode::encode_to_vec(&source, config::standard()).unwrap();
            let (decoded, consumed): (SymbolicTensor<PartialStructure>, _) =
                bincode::decode_from_slice_with_context(
                    &bytes,
                    config::standard(),
                    StateMap::default(),
                )
                .unwrap();
            assert_eq!(consumed, bytes.len());
            assert_eq!(decoded, source);
            assert_eq!(decoded.matrix_channel(), expected);
        }
    }

    #[test]
    fn tensor_archive_preserves_canonical_composite_axes_and_zero_history() {
        test_support::test_initialize();
        let representation = Minkowski {}.new_rep(4).to_lib();
        let ports = [98197, 98191]
            .map(|index| representation.slot(PartialIndex::Explicit(AbstractIndex::Normal(index))));
        let atom = function!(
            tensor_symbol!("archive_composite_axes"),
            representation
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(98197))
                .to_atom(),
            representation
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(98191))
                .to_atom()
        );
        let leaf = SymbolicTensor::from_declared_parts(
            atom.clone(),
            PartialStructure::from_logical_slots(ports),
        )
        .unwrap();
        let scalar = SymbolicTensor::infer(Atom::num(2)).unwrap();
        let probe = SymbolicTensor::infer(function!(
            tensor_symbol!("archive_composite_probe"),
            representation
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(98191))
                .to_atom()
        ))
        .unwrap();
        for source in [
            leaf.multiply(&scalar).unwrap(),
            leaf.try_sub(&leaf).unwrap(),
        ] {
            assert_eq!(source.structure.logical_slots(), [ports[1], ports[0]]);
            let bytes = bincode::encode_to_vec(&source, config::standard()).unwrap();
            let (decoded, consumed): (SymbolicTensor<PartialStructure>, _) =
                bincode::decode_from_slice_with_context(
                    &bytes,
                    config::standard(),
                    StateMap::default(),
                )
                .unwrap();
            assert_eq!(consumed, bytes.len());
            assert_eq!(decoded, source);
            assert_eq!(decoded.presentation_atom(), source.expression);
            let contracted = decoded.multiply(&probe).unwrap();
            assert_eq!(contracted.structure.logical_slots(), [ports[0]]);
            assert!(matches!(
                contracted.multiply(&probe),
                Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
            ));
        }
        assert_eq!(leaf.structure.logical_slots(), ports);
        assert_eq!(leaf.expression, atom);
    }

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
        assert!(decoded.proofs.is_admitted());
        assert_eq!(
            decoded.proofs.normalized_product.get(),
            source.proofs.normalized_product.get()
        );
        assert!(decoded.work.observations().is_none());
        assert!(decoded.work.graph().is_none());
    }

    #[test]
    fn tensor_intake_rejects_invalid_parts_before_encoding() {
        let structure = PartialStructure::from_logical_slots([Minkowski {}
            .new_rep(4)
            .to_lib()
            .slot(PartialIndex::Open(OpenPortId(0)))]);
        assert!(SymbolicTensor::from_declared_parts(Atom::one(), structure).is_err());
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
                    TENSOR_ARCHIVE_MAGIC,
                    version,
                    bincode::encode_to_vec(Atom::Zero, config::standard()).unwrap(),
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
                TENSOR_ARCHIVE_MAGIC,
                TENSOR_ARCHIVE_VERSION,
                bincode::encode_to_vec(Atom::Zero, config::standard()).unwrap(),
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
                    TENSOR_ARCHIVE_MAGIC,
                    TENSOR_ARCHIVE_VERSION,
                    bincode::encode_to_vec(Atom::Zero, config::standard()).unwrap(),
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
        let opaque = SymbolicTensor::<PartialStructure>::infer(slot.to_atom()).unwrap();
        assert!(opaque.is_scalar());
        assert!(opaque.explicit_index_occurrences().unwrap().is_empty());
        // Declared intake retains hidden index uses; opaque inference does not.
        let scalar = SymbolicTensor::from_declared_parts(
            slot.to_atom(),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
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
        for source in [opaque, scalar, product] {
            let expected = source.explicit_index_occurrences().unwrap().count(slot);
            INFERENCE_CALLS.with(|count| count.set(0));
            INFERENCE_VISITS.with(|count| count.set(0));
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
            assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
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
        INFERENCE_VISITS.with(|count| count.set(0));
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
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
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

    fn indexed_pair(dimension: Symbol, index: Symbol) -> SymbolicTensor<PartialStructure> {
        let slot = Minkowski {}
            .new_rep(dimension)
            .to_symbolic([Atom::var(index)]);
        SymbolicTensor::infer(
            function!(tensor_symbol!("idenso::persistence_test::pair_left"), &slot)
                * function!(
                    tensor_symbol!("idenso::persistence_test::pair_right"),
                    &slot
                ),
        )
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
            let nested = indexed_pair(dimension, index)
                .try_add(&indexed_pair(dimension, other_index))
                .unwrap();
            let nested = bincode::encode_to_vec(&nested, config::standard()).unwrap();
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
            println!("TENSOR_NESTED_BYTES={}", STANDARD.encode(nested));
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
        // The two sum branches reverse order when their symbolic indices are
        // remapped. Their retained leaf and regional interfaces must follow
        // the imported atoms, rather than their old preorder positions.
        let map = State::import(&mut Cursor::new(&state), None).unwrap();
        let bytes = STANDARD.decode(field("TENSOR_NESTED_BYTES=")).unwrap();
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let (nested, used): (SymbolicTensor<PartialStructure>, _) =
            bincode::decode_from_slice_with_context(&bytes, config::standard(), map).unwrap();
        assert_eq!(used, bytes.len());
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        for index in [index, other_index] {
            let expected = indexed_pair(dimension, index);
            assert_eq!(
                nested
                    .proofs
                    .retained_interface(expected.expression().as_view()),
                Some(expected.structure.clone()),
            );
            for factor in expected.expression().as_view().children() {
                assert!(nested.proofs.retained_interface(factor).is_some());
                assert!(
                    nested
                        .admitted_leaf_interfaces()
                        .lock()
                        .unwrap()
                        .contains_key(factor.get_data())
                );
            }
        }
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
    #[test]
    fn tensor_archive_rejects_legacy_and_crafted_lengths_before_atom_decode() {
        let legacy_atom = bincode::encode_to_vec(Atom::one(), config::standard()).unwrap();
        let mut old_tensor = vec![1u8];
        old_tensor.extend(&legacy_atom);
        let mut crafted = vec![1u8, 1u8];
        crafted.extend(u64::MAX.to_le_bytes());
        for bytes in [legacy_atom, old_tensor, crafted] {
            let error = bincode::decode_from_slice_with_context::<
                _,
                SymbolicTensor<PartialStructure>,
                _,
            >(&bytes, config::standard(), StateMap::default())
            .unwrap_err();
            assert!(matches!(
                error,
                DecodeError::Other(
                    "unrecognized symbolic tensor archive magic; legacy archives are unsupported"
                )
            ));
        }
    }

    #[test]
    fn tensor_archive_rejects_unknown_version_before_atom_decode() {
        let bytes = bincode::encode_to_vec(
            (TENSOR_ARCHIVE_MAGIC, TENSOR_ARCHIVE_VERSION + 1),
            config::standard(),
        )
        .unwrap();
        let error =
            bincode::decode_from_slice_with_context::<_, SymbolicTensor<PartialStructure>, _>(
                &bytes,
                config::standard(),
                StateMap::default(),
            )
            .unwrap_err();
        assert!(matches!(
            error,
            DecodeError::UnexpectedVariant {
                type_name: "symbolic tensor archive version",
                ..
            }
        ));
    }
}
