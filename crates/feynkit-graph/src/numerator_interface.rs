//! Source notation for retained numerator payloads and their logical ports.

use idenso::{CookMode, CookSettings, tensor::SymbolicTensor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use spenso::structure::{
    Canonicalized, OrderedStructure, TensorStructure,
    abstract_index::AbstractIndex,
    partial::{OpenPortId, PartialIndex, PartialStructure, PartialStructureExt},
    representation::LibraryRep,
    slot::{AbsInd, IsAbstractSlot, SlotMatcher},
};
use symbolica::atom::{Atom, AtomCore, AtomView};

use crate::{DiagramError, FeynmanDiagram};

/// Ordered source-artifact ports for non-scalar numerators, including computed zeros.
/// Declared logical axes survive different symbol-registration orders when source text
/// is loaded elsewhere. Compact representations denote distinct unresolved ports.
#[derive(Clone, Debug)]
pub struct NumeratorInterface(
    /// Source port declarations in logical order. `admit` validates them before
    /// they become the interface of a retained tensor.
    pub Vec<Atom>,
);

impl NumeratorInterface {
    pub fn parse(source: &str) -> Result<Self, DiagramError> {
        let ports: Vec<String> = serde_json::from_str(source)?;
        Ok(Self(
            ports
                .into_iter()
                .map(|port| FeynmanDiagram::parse_expression("numerator interface", port))
                .collect::<Result<_, _>>()?,
        ))
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string(self).expect("strings always serialize to JSON")
    }

    /// Admit authored source notation once, checking any declared logical ports.
    pub fn admit(
        expression: Atom,
        interface: Option<Self>,
    ) -> Result<SymbolicTensor<PartialStructure>, DiagramError> {
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let cook = |value: &Atom| {
            cooking
                .try_cook(value.as_view())
                .map_err(|error| DiagramError::Invariant {
                    operation: "admitting tensor source notation",
                    message: error.to_string(),
                })
        };
        let expression = cook(&expression)?;
        let Some(interface) = interface else {
            return Ok(SymbolicTensor::infer(expression)?);
        };
        let mut matcher = SlotMatcher::default();
        let ports = interface
            .0
            .into_iter()
            .map(|port| {
                let port = cook(&port)?;
                let invalid = |error: String| DiagramError::Invariant {
                    operation: "admitting tensor source interface",
                    message: error,
                };
                if matcher.compact_representation(port.as_view()).is_some() {
                    Ok(matcher
                        .parse_representation::<LibraryRep>(port.as_view())
                        .map_err(|error| invalid(error.to_string()))?
                        .slot(PartialIndex::Open(OpenPortId(0))))
                } else {
                    let slot = matcher
                        .parse::<LibraryRep, AbstractIndex>(port.as_view())
                        .map_err(|error| invalid(error.to_string()))?;
                    Ok(slot.rep().slot(PartialIndex::Explicit(slot.aind)))
                }
            })
            .collect::<Result<Vec<_>, DiagramError>>()?;
        Ok(SymbolicTensor::from_declared_parts(
            expression,
            PartialStructure::from_logical_slots(ports),
        )?)
    }

    /// Export source notation with its complete nonempty logical interface.
    pub fn export<A>(
        tensor: &SymbolicTensor<Canonicalized<OrderedStructure<LibraryRep, A>>>,
    ) -> (Atom, Option<Self>)
    where
        A: AbsInd + Into<PartialIndex<AbstractIndex>>,
    {
        let canonical = tensor.structure().canonical().external_structure();
        let logical = tensor
            .structure()
            .layout()
            .canonical_to_logical(&canonical)
            .into_iter()
            .map(|slot| slot.rep().slot(slot.aind.into()))
            .collect::<Vec<_>>();
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let expression = cooking.uncook(tensor.expression().as_view());
        if logical.is_empty() {
            return (expression, None);
        }
        let ports = logical
            .into_iter()
            .map(|slot| {
                let port = match slot.aind {
                    PartialIndex::Explicit(index) => {
                        slot.rep().slot::<AbstractIndex, _>(index).to_atom()
                    }
                    PartialIndex::Open(_) => slot.rep().to_symbolic([]),
                };
                cooking.uncook(port.as_view())
            })
            .collect();
        (expression, Some(Self(ports)))
    }

    /// Transport source-owned labels while retaining the live tensor carrier.
    /// Only individual encoded index payloads are decoded. The native rename
    /// owner maps their occurrences and logical ports simultaneously, including
    /// ports carried only by a zero; remaining momentum/scalar heads are mapped
    /// in the native expression after index transport.
    pub(crate) fn translate(
        tensor: &SymbolicTensor<PartialStructure>,
        translate: impl Fn(&Atom) -> Result<Atom, DiagramError>,
    ) -> Result<SymbolicTensor<PartialStructure>, DiagramError> {
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let tensor = tensor.with_cooked_indices(&cooking)?;
        let mut symbols = std::collections::HashSet::new();
        let mut collect = |value: &Atom| {
            value.visitor(&mut |node| {
                if let AtomView::Var(variable) = node {
                    symbols.insert(variable.get_symbol());
                }
                true
            });
        };
        collect(tensor.expression());
        for slot in tensor.structure().logical_slots() {
            if let PartialIndex::Explicit(index) = slot.aind {
                collect(&slot.rep().slot::<AbstractIndex, _>(index).to_atom());
            }
        }
        let mut renaming = std::collections::HashMap::new();
        for symbol in symbols {
            let original = Atom::var(symbol);
            let decoded = cooking.uncook(original.as_view());
            if decoded == original {
                continue;
            }
            let mapped = translate(&decoded)?;
            if mapped == decoded {
                continue;
            }
            let Some(function) = mapped.as_fun_view() else {
                return Err(DiagramError::Invariant {
                    operation: "transporting encoded diagram indices",
                    message: "an index-owner map must retain a function index".into(),
                });
            };
            let mapped =
                cooking
                    .cook_function(function)
                    .map_err(|error| DiagramError::Invariant {
                        operation: "transporting encoded diagram indices",
                        message: error.to_string(),
                    })?;
            let Some(target) = mapped.get_symbol() else {
                return Err(DiagramError::Invariant {
                    operation: "transporting encoded diagram indices",
                    message: "a transported index was not reversibly encoded".into(),
                });
            };
            renaming.insert(symbol, target);
        }
        let tensor = tensor.rename_symbols(&renaming)?;
        Ok(tensor.with_checked_expression(translate(tensor.expression())?)?)
    }
}

impl Serialize for NumeratorInterface {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.0
            .iter()
            .map(AtomCore::to_canonical_string)
            .collect::<Vec<_>>()
            .serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for NumeratorInterface {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let ports = Vec::<String>::deserialize(deserializer)?;
        ports
            .into_iter()
            .map(|port| FeynmanDiagram::parse_expression("numerator interface", port))
            .collect::<Result<Vec<_>, _>>()
            .map(Self)
            .map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json as json;
    use spenso::network::tags::SPENSO_TAG;
    use spenso::structure::representation::{Euclidean, Minkowski, RepName};
    use std::env;
    use std::process::Command;
    use symbolica::function;
    use symbolica::id::Replacement;

    #[test]
    fn source_codec_omits_scalar_interfaces() {
        for input in [Atom::zero(), Atom::one(), Atom::num((2, 3))] {
            let source = NumeratorInterface::admit(input.clone(), None).unwrap();
            let (expression, metadata) = NumeratorInterface::export(&source);
            assert_eq!(expression, input);
            assert!(metadata.is_none());
            assert!(
                NumeratorInterface::admit(expression, metadata)
                    .unwrap()
                    .same_value(&source)
            );
        }
    }

    #[test]
    fn source_codec_keeps_ports_across_registration_order() {
        const PAYLOAD: &str = "FEYNKIT_INTERFACE_REGISTRATION_ORDER_TEST_PAYLOAD";
        let payload = env::var(PAYLOAD).ok();
        let names = [
            "feynkit_graph::portable_rank_one_a",
            "feynkit_graph::portable_rank_one_b",
        ];
        let (a, b) = if payload.is_some() {
            let b = SPENSO_TAG.rank_one_tensor_symbol(names[1]);
            let a = SPENSO_TAG.rank_one_tensor_symbol(names[0]);
            (a, b)
        } else {
            let a = SPENSO_TAG.rank_one_tensor_symbol(names[0]);
            let b = SPENSO_TAG.rank_one_tensor_symbol(names[1]);
            (a, b)
        };

        if let Some(payload) = payload {
            let (expression, interface): (String, String) = json::from_str(&payload).unwrap();
            let expression =
                FeynmanDiagram::parse_expression("test numerator", expression).unwrap();
            let inferred = NumeratorInterface::admit(expression.clone(), None).unwrap();
            let metadata = NumeratorInterface::parse(&interface).unwrap();
            let restored = NumeratorInterface::admit(expression, Some(metadata)).unwrap();
            assert_eq!(restored.expression(), inferred.expression());
            assert_ne!(
                restored.structure().logical_slots(),
                inferred.structure().logical_slots(),
                "opposite registration must expose the inference-order difference"
            );
            let (_, metadata) = NumeratorInterface::export(&restored);
            assert_eq!(
                metadata.unwrap().to_json(),
                interface,
                "the reader must preserve the writer's logical axis order"
            );
            return;
        }

        let ports = [1, 2].map(|index| Minkowski {}.new_rep(4).to_symbolic([Atom::num(index)]));
        let input = function!(a, &ports[0]) * function!(b, &ports[1]);
        let source = NumeratorInterface::admit(input.clone(), None).unwrap();
        let (expression, metadata) = NumeratorInterface::export(&source);
        assert_eq!(expression, input);
        assert_eq!(
            NumeratorInterface::admit(expression.clone(), None)
                .unwrap()
                .structure()
                .logical_slots(),
            source.structure().logical_slots()
        );
        let metadata = metadata.expect("portable source text requires the writer's ordered ports");
        let payload =
            json::to_string(&(expression.to_canonical_string(), metadata.to_json())).unwrap();
        let reader = Command::new(env::current_exe().unwrap())
            .args([
                "--exact",
                "numerator_interface::tests::source_codec_keeps_ports_across_registration_order",
                "--test-threads=1",
                "--nocapture",
            ])
            .env(PAYLOAD, payload)
            .output()
            .unwrap();
        let stdout = String::from_utf8_lossy(&reader.stdout);
        assert!(
            reader.status.success() && stdout.contains("1 passed; 0 failed;"),
            "fresh reader failed:\n{}\n{}",
            stdout,
            String::from_utf8_lossy(&reader.stderr)
        );
    }

    #[test]
    fn source_codec_preserves_anonymous_axes_and_ranked_zero() {
        let matrix = spenso::tensor_symbol!("feynkit_graph::source_matrix");
        let input = function!(
            matrix,
            Minkowski {}.new_rep(4).to_symbolic([]),
            Euclidean {}.new_rep(2).to_symbolic([])
        );
        let source = NumeratorInterface::admit(input, None).unwrap();
        let zero = source.with_checked_expression(Atom::Zero).unwrap();
        let (expression, metadata) = NumeratorInterface::export(&zero);
        let metadata = metadata.expect("a ranked zero requires declared ports");
        let imported = NumeratorInterface::admit(
            expression,
            Some(NumeratorInterface::parse(&metadata.to_json()).unwrap()),
        )
        .unwrap();
        assert!(imported.same_value(&zero));
        assert_eq!(
            imported.structure().logical_slots(),
            source.structure().logical_slots()
        );
        assert!(NumeratorInterface::admit(Atom::one(), Some(metadata)).is_err());
    }

    #[test]
    fn source_label_permutations_retain_zero_ports_without_decoding_payloads() {
        let matrix = spenso::tensor_symbol!("feynkit_graph::indexed_source_matrix");
        let indices = [1, 2].map(|owner| function!(crate::symbols::hedge_index(), owner, 1));
        let input = function!(
            matrix,
            Minkowski {}.new_rep(4).to_symbolic([indices[0].clone()]),
            Euclidean {}.new_rep(2).to_symbolic([indices[1].clone()])
        );
        let source = NumeratorInterface::admit(input.clone(), None).unwrap();
        let translate = |atom: &Atom| {
            Ok(atom.replace_multiple([
                Replacement::new(indices[0].to_pattern(), indices[1].clone()),
                Replacement::new(indices[1].to_pattern(), indices[0].clone()),
            ]))
        };
        let expected = NumeratorInterface::admit(translate(&input).unwrap(), None).unwrap();
        let mapped = NumeratorInterface::translate(&source, translate).unwrap();
        assert!(mapped.same_value(&expected));
        let zero = source.with_checked_expression(Atom::Zero).unwrap();
        let mapped_zero = NumeratorInterface::translate(&zero, translate).unwrap();
        assert!(mapped_zero.is_zero());
        assert_eq!(
            mapped_zero.structure().logical_slots(),
            expected.structure().logical_slots()
        );
        assert!(
            NumeratorInterface::translate(&source, |atom| Ok(atom
                .replace(indices[0].to_pattern())
                .with(indices[1].clone())))
            .is_err()
        );
    }
}
