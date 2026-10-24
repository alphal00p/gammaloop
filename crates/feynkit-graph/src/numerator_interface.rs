//! Source notation for retained numerator payloads and their logical ports.

use std::collections::{HashMap, HashSet};

use idenso::tensor::SymbolicTensor;
use idenso::tensor::inference::InterfaceInference;
use idenso::{CookMode, CookSettings};
use serde::de::Error as DeserializeError;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{OpenPortId, PartialIndex, PartialStructure, PartialStructureExt};
use spenso::structure::representation::LibraryRep;
use spenso::structure::slot::{AbsInd, IsAbstractSlot, ParseableAind, SlotMatch, SlotMatcher};
use spenso::structure::{Canonicalized, OrderedStructure, TensorStructure};
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};

use crate::{DiagramError, FeynmanDiagram};

/// Optional source-artifact ports that cannot be recovered from the expression,
/// notably a computed zero's rank. Entries retain logical axis order; compact
/// representations denote distinct unresolved ports.
///
/// # Example
///
/// ```
/// use feynkit_graph::NumeratorInterface;
/// use symbolica::atom::Atom;
///
/// let ports = NumeratorInterface::parse(r#"["spenso::mink(4)"]"#)?;
/// let zero = NumeratorInterface::admit(Atom::Zero, Some(ports))?;
/// let (expression, interface) = NumeratorInterface::export(&zero);
/// let interface = interface.expect("zero syntax cannot carry its axis");
/// let restored = NumeratorInterface::admit(
///     expression, Some(NumeratorInterface::parse(&interface.to_json())?),
/// )?;
/// assert!(restored.same_value(&zero));
/// # Ok::<(), feynkit_graph::DiagramError>(())
/// ```
#[derive(Clone, Debug)]
pub struct NumeratorInterface(
    /// Source port declarations in logical order. `admit` validates them before
    /// they become the interface of a retained tensor.
    pub Vec<Atom>,
);

impl NumeratorInterface {
    /// Decode source port declarations from their JSON string-array representation.
    /// See the round trip in [`NumeratorInterface`].
    ///
    /// # Errors
    /// Returns a JSON or symbolic parse error for malformed declarations.
    pub fn parse(source: &str) -> Result<Self, DiagramError> {
        let ports: Vec<String> = serde_json::from_str(source)?;
        Ok(Self(
            ports
                .into_iter()
                .map(|port| FeynmanDiagram::parse_expression("numerator interface", port))
                .collect::<Result<_, _>>()?,
        ))
    }

    /// Encode logical ports in canonical source notation; see [`NumeratorInterface`].
    #[must_use]
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).expect("strings always serialize to JSON")
    }

    /// Admit authored source notation once, checking any declared logical ports.
    /// See [`NumeratorInterface`] for a zero with a retained axis.
    ///
    /// # Errors
    /// Returns the cooking, slot, or tensor-inference error when notation or
    /// declared ports are invalid or inconsistent with the expression.
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
                .map_err(|error| DiagramError::TensorCooking {
                    operation: "admitting tensor source notation",
                    source: Box::new(error),
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
                let invalid = |source| DiagramError::TensorSlot {
                    operation: "admitting tensor source interface",
                    source: Box::new(source),
                };
                if matcher.compact_representation(port.as_view()).is_some() {
                    Ok(matcher
                        .parse_representation::<LibraryRep>(port.as_view())
                        .map_err(invalid)?
                        .slot(PartialIndex::Open(OpenPortId(0))))
                } else {
                    let slot = matcher
                        .parse::<LibraryRep, AbstractIndex>(port.as_view())
                        .map_err(invalid)?;
                    Ok(slot.rep().slot(PartialIndex::Explicit(slot.aind)))
                }
            })
            .collect::<Result<Vec<_>, DiagramError>>()?;
        Ok(SymbolicTensor::from_declared_parts(
            expression,
            PartialStructure::from_logical_slots(ports),
        )?)
    }

    /// Export source notation, retaining metadata whenever syntax loses ports.
    /// See the round trip in [`NumeratorInterface`].
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
        if Self::admit(expression.clone(), None)
            .is_ok_and(|restored| restored.structure().logical_slots() == logical)
        {
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
    /// Native indices and individual reversible labels share one injective map,
    /// including indices carried only by a zero. The original expression is
    /// translated once; checked tensor parts retain the established port order.
    pub(crate) fn translate(
        tensor: &SymbolicTensor<PartialStructure>,
        translate: impl Fn(&Atom) -> Result<Atom, DiagramError>,
    ) -> Result<SymbolicTensor<PartialStructure>, DiagramError> {
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let invalid = |message: String| DiagramError::Invariant {
            operation: "transporting diagram indices",
            message,
        };
        let mut indices = HashSet::new();
        let mut slots = SlotMatcher::default();
        let mut collect = |value: &Atom| {
            let mut failure = None;
            value.visitor(&mut |node| {
                if let SlotMatch::Explicit(slot) = slots.classify(node) {
                    match AbstractIndex::try_from(slot.index()) {
                        Ok(index) => {
                            indices.insert(index);
                        }
                        Err(source) => {
                            failure = Some(DiagramError::TensorIndex {
                                operation: "collecting transported tensor indices",
                                source: Box::new(source),
                            })
                        }
                    }
                }
                if let AtomView::Var(variable) = node {
                    // Include scalar symbols too, so an encoded label cannot
                    // capture an existing coefficient or a scoped identity.
                    indices.insert(AbstractIndex::from(variable.get_symbol()));
                }
                failure.is_none()
            });
            failure.map_or(Ok(()), Err)
        };
        collect(tensor.expression())?;
        for slot in tensor.structure().logical_slots() {
            if let PartialIndex::Explicit(index) = slot.aind {
                collect(&slot.rep().slot::<AbstractIndex, _>(index).to_atom())?;
            }
        }
        let mut renaming = HashMap::new();
        let mut encoded = HashMap::<Symbol, Atom>::new();
        let mut targets = HashSet::new();
        for index in indices {
            let original = index.to_atom();
            let decoded = cooking.uncook(original.as_view());
            let mut mapped = translate(&decoded)?;
            let target = if mapped == decoded {
                index
            } else {
                // Scopes are already native index structure. Re-encode only
                // their terminal label after decoding and mapping its owner.
                let mut scopes = Vec::new();
                while let Some(function) = mapped.as_fun_view()
                    && function.get_symbol() == AbstractIndex::scope_symbol()
                {
                    let arguments = function.iter().collect::<Vec<_>>();
                    let [AtomView::Var(scope), inner] = arguments.as_slice() else {
                        return Err(invalid(
                            "owner mapping produced an invalid index scope".into(),
                        ));
                    };
                    scopes.push(scope.get_symbol());
                    mapped = (*inner).to_owned();
                }
                let mapped = if let Some(function) = mapped.as_fun_view() {
                    cooking.cook_function(function).map_err(|source| {
                        DiagramError::TensorCooking {
                            operation: "encoding transported tensor indices",
                            source: Box::new(source),
                        }
                    })?
                } else {
                    mapped
                };
                scopes.into_iter().rev().fold(
                    AbstractIndex::try_from(mapped.as_view()).map_err(|source| {
                        DiagramError::TensorIndex {
                            operation: "parsing transported tensor indices",
                            source: Box::new(source),
                        }
                    })?,
                    |index, scope| index.scoped(scope),
                )
            };
            if !targets.insert(target) {
                return Err(invalid("index-owner transport must be injective".into()));
            }
            renaming.insert(index, target);
            if decoded != original && target != index {
                // Native functions are mapped by the expression translator.
                // Encoded variables remain opaque until that pass has finished,
                // otherwise a newly native target would be transported twice.
                if let AbstractIndex::Symbol(symbol) = index {
                    encoded.insert(symbol.into(), target.to_atom());
                }
            }
        }
        let original_slots = tensor.structure().logical_slots();
        let logical = original_slots
            .iter()
            .copied()
            .map(|mut slot| {
                if let PartialIndex::Explicit(index) = slot.aind {
                    slot.aind = PartialIndex::Explicit(renaming[&index]);
                }
                slot
            })
            .collect::<Vec<_>>();
        let expression = translate(tensor.expression())?;
        let expression = expression.replace_map(|node, _, output| {
            if let AtomView::Var(variable) = node
                && let Some(target) = encoded.get(&variable.get_symbol())
            {
                **output = target.clone();
            }
        });
        if !tensor.structure().open_positions().is_empty() && !expression.is_zero() {
            if logical == original_slots {
                return Ok(tensor.with_checked_expression(expression)?);
            }
            if !InterfaceInference::normalization_is_intrinsic(tensor.expression().as_view()) {
                return Err(invalid(
                    "cannot certify unresolved port identities while transporting explicit indices through callbacks"
                        .into(),
                ));
            }
            let replacements = original_slots
                .iter()
                .zip(&logical)
                .enumerate()
                .filter_map(|(position, (source, target))| match target.aind {
                    PartialIndex::Explicit(index) if source.aind != target.aind => {
                        Some((position, index))
                    }
                    _ => None,
                })
                .collect();
            // The native positional rewrite certifies the mapped boundary.
            // Its checked-expression guard then protects anonymous occurrences
            // from changes to the remaining momentum and scalar heads.
            let boundary = tensor.reindex_interface_ports(&replacements)?;
            return Ok(boundary.with_checked_expression(expression)?);
        }
        // Translation is an existing checked boundary: callback output must
        // return a typed error in release and debug builds alike. Admit the
        // declared result directly instead of claiming producer trust first.
        Ok(SymbolicTensor::from_declared_parts(
            expression,
            PartialStructure::from_logical_slots(logical),
        )?)
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
            .map_err(DeserializeError::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symbols;
    use spenso::structure::representation::{Euclidean, Minkowski, RepName};
    use spenso::tensor_symbol;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use symbolica::function;
    use symbolica::id::Replacement;
    use symbolica::symbol;

    #[test]
    fn source_codec_preserves_anonymous_axes_and_ranked_zero() {
        let matrix = tensor_symbol!("feynkit_graph::source_matrix");
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
    fn native_label_permutations_retain_zero_ports_and_reject_collisions() {
        let matrix = tensor_symbol!("feynkit_graph::indexed_source_matrix");
        let indices = [1, 2].map(|owner| function!(symbols::hedge_index(), owner, 1));
        let input = function!(
            matrix,
            Minkowski {}.new_rep(4).to_symbolic([indices[0].clone()]),
            Euclidean {}.new_rep(2).to_symbolic([indices[1].clone()])
        );
        let source = NumeratorInterface::admit(input.clone(), None).unwrap();
        assert_eq!(source.expression(), &input);
        assert!(
            source
                .structure()
                .logical_slots()
                .iter()
                .all(|slot| matches!(slot.aind, PartialIndex::Explicit(AbstractIndex::Named(..))))
        );
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

    #[test]
    fn native_owner_transport_preserves_scopes_and_declared_logical_order() {
        let head = tensor_symbol!("feynkit_graph::native_owner_transport");
        let scope = symbol!("feynkit_graph::source_scope");
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let source_indices = [
            AbstractIndex::Named(symbols::hedge_index().into(), 1, 0),
            AbstractIndex::Named(symbols::edge_index().into(), 2, 0).scoped(scope),
        ];
        let target_indices = [
            AbstractIndex::Named(symbols::hedge_index().into(), 7, 0),
            AbstractIndex::Named(symbols::edge_index().into(), 8, 0).scoped(scope),
        ];
        let expression = function!(
            head,
            rep.slot::<AbstractIndex, _>(source_indices[0]).to_atom(),
            rep.slot::<AbstractIndex, _>(source_indices[1]).to_atom()
        );
        let source = SymbolicTensor::from_declared_parts(
            expression,
            PartialStructure::from_logical_slots(
                source_indices
                    .into_iter()
                    .rev()
                    .map(|index| rep.slot(PartialIndex::Explicit(index))),
            ),
        )
        .unwrap();
        let translate = |atom: &Atom| {
            Ok(atom.replace_multiple([
                Replacement::new(
                    source_indices[0].to_atom().to_pattern(),
                    target_indices[0].to_atom(),
                ),
                Replacement::new(
                    source_indices[1].to_atom().to_pattern(),
                    target_indices[1].to_atom(),
                ),
            ]))
        };
        for value in [
            source.clone(),
            source.with_checked_expression(Atom::Zero).unwrap(),
        ] {
            let expected = value
                .reindex_interface_ports(&HashMap::from([
                    (0, target_indices[1]),
                    (1, target_indices[0]),
                ]))
                .unwrap();
            let mapped = NumeratorInterface::translate(&value, translate).unwrap();
            assert!(mapped.same_value(&expected));
            assert_eq!(mapped.is_zero(), value.is_zero());
            assert_eq!(
                mapped.structure().logical_slots(),
                target_indices
                    .into_iter()
                    .rev()
                    .map(|index| rep.slot(PartialIndex::Explicit(index)))
                    .collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn native_transport_keeps_anonymous_and_interface_only_ports() {
        let rep = Minkowski {}.new_rep(4);
        let source_index = AbstractIndex::Named(symbols::hedge_index().into(), 1, 0);
        let target_index = AbstractIndex::Named(symbols::hedge_index().into(), 9, 0);
        let source = NumeratorInterface::admit(
            function!(
                tensor_symbol!("feynkit_graph::mixed_source_ports"),
                rep.to_symbolic([]),
                rep.slot::<AbstractIndex, _>(source_index).to_atom(),
                Euclidean {}.new_rep(2).to_symbolic([])
            ),
            None,
        )
        .unwrap();
        let translate = |atom: &Atom| {
            Ok(atom
                .replace(source_index.to_atom().to_pattern())
                .with(target_index.to_atom()))
        };
        let mut expected = source.structure().logical_slots();
        expected[1].aind = PartialIndex::Explicit(target_index);
        for value in [
            source.clone(),
            source.with_checked_expression(Atom::Zero).unwrap(),
        ] {
            let mapped = NumeratorInterface::translate(&value, translate).unwrap();
            assert_eq!(mapped.structure().logical_slots(), expected);
            assert_eq!(mapped.structure().open_positions(), vec![0, 2]);
            assert_eq!(mapped.is_zero(), value.is_zero());
            if !mapped.is_zero() {
                assert_eq!(
                    mapped.expression(),
                    &translate(source.expression()).unwrap()
                );
            }
        }
    }

    #[test]
    fn anonymous_owner_transport_rejects_ambiguous_factor_reordering() {
        let rep = Minkowski {}.new_rep(4);
        let momentum = symbols::momentum();
        let source_momenta = [0, 1].map(|owner| function!(momentum, owner, rep.to_symbolic([])));
        let target_momenta = [3, 2].map(|owner| function!(momentum, owner, rep.to_symbolic([])));
        let source =
            NumeratorInterface::admit(&source_momenta[0] * &source_momenta[1], None).unwrap();
        let translate = |atom: &Atom| {
            Ok(atom.replace_multiple([
                Replacement::new(source_momenta[0].to_pattern(), target_momenta[0].clone()),
                Replacement::new(source_momenta[1].to_pattern(), target_momenta[1].clone()),
            ]))
        };
        let expression = translate(source.expression()).unwrap();
        assert!(source.with_checked_expression(expression).is_err());
        assert!(NumeratorInterface::translate(&source, translate).is_err());
    }

    #[test]
    fn reversible_labels_remain_supported_and_can_transport_to_native_indices() {
        let head = tensor_symbol!("feynkit_graph::encoded_source_ports");
        let rep = Minkowski {}.new_rep(4);
        let source_index = function!(symbols::hedge_index(), 1);
        let target_index = function!(symbols::hedge_index(), 9);
        let opaque = symbol!("feynkit_graph::opaque_source_label");
        let translate = |atom: &Atom| {
            Ok(atom
                .replace(source_index.to_pattern())
                .with(target_index.clone()))
        };
        for (label, encoding) in [
            (
                function!(opaque, &source_index),
                CookSettings::indices().with_mode(CookMode::ReversibleEncoding),
            ),
            (source_index.clone(), CookSettings::reversible()),
        ] {
            let encoded = encoding
                .cook_function(label.as_fun_view().unwrap())
                .unwrap();
            assert!(encoded.get_symbol().is_some());
            let source =
                SymbolicTensor::infer(function!(head, rep.to_symbolic([encoded]))).unwrap();
            let expected = NumeratorInterface::admit(
                function!(head, rep.to_symbolic([translate(&label).unwrap()])),
                None,
            )
            .unwrap();
            for value in [
                source.clone(),
                source.with_checked_expression(Atom::Zero).unwrap(),
            ] {
                let mapped = NumeratorInterface::translate(&value, translate).unwrap();
                assert_eq!(
                    mapped.structure().logical_slots(),
                    expected.structure().logical_slots()
                );
                if !value.is_zero() {
                    assert!(mapped.same_value(&expected));
                } else {
                    assert!(mapped.is_zero());
                }
            }
            let slot = expected.structure().logical_slots()[0];
            assert_eq!(
                matches!(slot.aind, PartialIndex::Explicit(AbstractIndex::Named(..))),
                label == source_index
            );
        }
    }

    #[test]
    fn scoped_encoded_labels_transport_their_inner_owner_and_zero_ports() {
        let rep = Minkowski {}.new_rep(4);
        let head = tensor_symbol!("feynkit_graph::scoped_encoded_transport");
        let scope = symbol!("feynkit_graph::encoded_transport_scope");
        let original = function!(
            symbol!("feynkit_graph::scoped_opaque_label"),
            function!(symbols::hedge_index(), 1)
        );
        let translated = original
            .replace(function!(symbols::hedge_index(), 1).to_pattern())
            .with(function!(symbols::hedge_index(), 9));
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let encoded = cooking
            .cook_function(original.as_fun_view().unwrap())
            .unwrap();
        let target = cooking
            .cook_function(translated.as_fun_view().unwrap())
            .unwrap();
        let target = AbstractIndex::try_from(target.as_view())
            .unwrap()
            .scoped(scope);
        let scoped = AbstractIndex::try_from(encoded.as_view())
            .unwrap()
            .scoped(scope);
        let legacy = CookSettings::reversible()
            .cook_function(
                function!(AbstractIndex::scope_symbol(), Atom::var(scope), original)
                    .as_fun_view()
                    .unwrap(),
            )
            .unwrap();
        for source_index in [scoped.to_atom(), legacy] {
            let source =
                SymbolicTensor::infer(function!(head, rep.to_symbolic([source_index]))).unwrap();
            for value in [
                source.clone(),
                source.with_checked_expression(Atom::Zero).unwrap(),
            ] {
                let expected = value
                    .reindex_interface_ports(&HashMap::from([(0, target)]))
                    .unwrap();
                let mapped = NumeratorInterface::translate(&value, |atom| {
                    Ok(atom
                        .replace(function!(symbols::hedge_index(), 1).to_pattern())
                        .with(function!(symbols::hedge_index(), 9)))
                })
                .unwrap();
                assert!(mapped.same_value(&expected));
                assert_eq!(
                    mapped.structure().logical_slots(),
                    expected.structure().logical_slots()
                );
            }
        }
    }

    #[test]
    fn native_and_encoded_indices_cannot_capture_each_other() {
        let head = tensor_symbol!("feynkit_graph::mixed_identity_collision");
        let rep = Minkowski {}.new_rep(4);
        let native = function!(symbols::edge_index(), 2);
        let opaque = function!(symbol!("feynkit_graph::collision_label"), 7);
        let source = NumeratorInterface::admit(
            function!(
                head,
                rep.to_symbolic([native.clone()]),
                rep.to_symbolic([opaque.clone()])
            ),
            None,
        )
        .unwrap();
        for value in [
            source.clone(),
            source.with_checked_expression(Atom::Zero).unwrap(),
        ] {
            for (from, to) in [(&native, &opaque), (&opaque, &native)] {
                assert!(
                    NumeratorInterface::translate(&value, |atom| {
                        Ok(atom.replace(from.to_pattern()).with(to.clone()))
                    })
                    .is_err()
                );
            }
        }
    }

    #[test]
    fn native_transport_keeps_internal_contractions_and_rejects_capture() {
        let rep = Minkowski {}.new_rep(4);
        let [dummy, external, new_dummy, new_external] =
            [1, 2, 7, 8].map(|owner| function!(symbols::hedge_index(), owner));
        let left = tensor_symbol!("feynkit_graph::contracted_left");
        let right = tensor_symbol!("feynkit_graph::contracted_right");
        let input = function!(left, rep.to_symbolic([dummy.clone()]))
            * function!(
                right,
                rep.to_symbolic([dummy.clone()]),
                rep.to_symbolic([external.clone()])
            );
        let source = NumeratorInterface::admit(input.clone(), None).unwrap();
        assert_eq!(source.rank(), 1);
        let translate = |atom: &Atom| {
            Ok(atom.replace_multiple([
                Replacement::new(dummy.to_pattern(), new_dummy.clone()),
                Replacement::new(external.to_pattern(), new_external.clone()),
            ]))
        };
        let expected = NumeratorInterface::admit(translate(&input).unwrap(), None).unwrap();
        let mapped = NumeratorInterface::translate(&source, translate).unwrap();
        assert!(mapped.same_value(&expected));
        assert_eq!(mapped.rank(), 1);
        assert!(mapped.expression().contains(&new_dummy));
        assert!(!mapped.expression().contains(&dummy));
        assert!(
            NumeratorInterface::translate(&source, |atom| {
                Ok(atom.replace(dummy.to_pattern()).with(external.clone()))
            })
            .is_err()
        );
    }

    #[test]
    fn native_transport_checks_callback_port_loss_without_replaying_normalization() {
        let source_index = function!(symbols::edge_index(), 1);
        let target_index = function!(symbols::edge_index(), 2);
        let rep = Minkowski {}.new_rep(4);
        let target = rep.to_symbolic([target_index.clone()]);
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let head = tensor_symbol!(
            "feynkit_graph::transport_callback_port_loss",
            norm = move |node, output| {
                observed.fetch_add(1, Ordering::Relaxed);
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let source = NumeratorInterface::admit(
            function!(head, rep.to_symbolic([source_index.clone()])),
            None,
        )
        .unwrap();
        calls.store(0, Ordering::Relaxed);
        assert!(
            NumeratorInterface::translate(&source, |atom| {
                Ok(atom
                    .replace(source_index.to_pattern())
                    .with(target_index.clone()))
            })
            .is_err()
        );
        assert_eq!(calls.load(Ordering::Relaxed), 1);
    }
}
