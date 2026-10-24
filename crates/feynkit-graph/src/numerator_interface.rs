//! Source notation for retained numerator payloads and their logical ports.

use std::{cell::RefCell, collections::HashMap, rc::Rc};

use idenso::{CookMode, CookSettings, tensor::SymbolicTensor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use spenso::structure::{
    Canonicalized, OrderedStructure, TensorStructure,
    abstract_index::AbstractIndex,
    partial::{OpenPortId, PartialIndex, PartialSlot, PartialStructure, PartialStructureExt},
    representation::LibraryRep,
    slot::{AbsInd, IsAbstractSlot, SlotMatcher},
};
use symbolica::atom::{Atom, AtomCore, AtomView};

use crate::{DiagramError, FeynmanDiagram};

/// Optional source-artifact ports that cannot be recovered from the expression,
/// notably a computed zero's rank. Entries retain logical axis order; compact
/// representations denote distinct unresolved ports.
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
            return Ok(Self::canonical_ports(SymbolicTensor::infer(expression)?));
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

    /// Order explicit logical ports by their source notation. Inference follows
    /// the in-memory term order of a product, which depends on the order in
    /// which Symbolica registered its symbols and so differs between processes;
    /// admitted and serialized interfaces must not. Open ports are identified by
    /// position, so a tensor with any keeps its order.
    pub(crate) fn canonical_ports(
        tensor: SymbolicTensor<PartialStructure>,
    ) -> SymbolicTensor<PartialStructure> {
        let Some(mut ports) = Self::keyed_ports(&tensor) else {
            return tensor;
        };
        if Self::in_canonical_order(&ports) {
            return tensor;
        }
        ports.sort_by(|(left, _), (right, _)| left.cmp(right));
        // Reordering retains every port, so the established interface stays valid.
        // `SymbolicTensor::permuted` would first rewrite and re-check the whole
        // expression to name open ports, which explicit ports do not need.
        SymbolicTensor::from_normalized_parts(
            tensor.into_expression(),
            PartialStructure::from_logical_slots(ports.into_iter().map(|(_, slot)| slot)),
        )
    }

    /// Logical ports with their source notation, or `None` if any port is open.
    fn keyed_ports(
        tensor: &SymbolicTensor<PartialStructure>,
    ) -> Option<Vec<(Rc<str>, PartialSlot)>> {
        tensor
            .structure()
            .logical_slots()
            .into_iter()
            .map(|slot| Self::port_key(slot).map(|key| (key, slot)))
            .collect()
    }

    fn in_canonical_order(ports: &[(Rc<str>, PartialSlot)]) -> bool {
        ports.is_sorted_by(|(left, _), (right, _)| left <= right)
    }

    /// The source notation of an explicit port, or `None` for an open port.
    fn port_key(slot: PartialSlot) -> Option<Rc<str>> {
        thread_local! {
            // A key depends only on its slot, and the same few external and
            // half-edge ports recur in every fragment and export of a process.
            static KEYS: RefCell<HashMap<PartialSlot, Rc<str>>> = RefCell::default();
        }
        let PartialIndex::Explicit(index) = slot.aind else {
            return None;
        };
        Some(KEYS.with_borrow_mut(|keys| {
            keys.entry(slot)
                .or_insert_with(|| {
                    let port = slot.rep().slot::<AbstractIndex, _>(index).to_atom();
                    CookSettings::indices()
                        .with_mode(CookMode::ReversibleEncoding)
                        .uncook(port.as_view())
                        .to_canonical_string()
                        .into()
                })
                .clone()
        }))
    }

    /// Export source notation, retaining metadata whenever syntax loses ports.
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
    /// Only individual encoded index payloads are decoded. The native rename
    /// owner maps their occurrences and logical ports simultaneously, including
    /// ports carried only by a zero; remaining momentum/scalar heads are mapped
    /// in the native expression after index transport. A canonical port order,
    /// declared or inferred, stays canonical for the new labels; any other
    /// order, which only a declaration produces, moves positionally.
    pub(crate) fn translate(
        tensor: &SymbolicTensor<PartialStructure>,
        translate: impl Fn(&Atom) -> Result<Atom, DiagramError>,
    ) -> Result<SymbolicTensor<PartialStructure>, DiagramError> {
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let tensor = tensor.with_cooked_indices(&cooking)?;
        let canonical =
            Self::keyed_ports(&tensor).is_some_and(|ports| Self::in_canonical_order(&ports));
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
        let tensor = tensor.with_checked_expression(translate(tensor.expression())?)?;
        Ok(if canonical {
            Self::canonical_ports(tensor)
        } else {
            tensor
        })
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
    use spenso::structure::representation::{Euclidean, Minkowski, RepName};
    use symbolica::{function, id::Replacement};

    #[test]
    fn inferred_ports_are_independent_of_symbol_registration_order() {
        let vector = |head, owner: usize| {
            let index = function!(crate::symbols::hedge_index(), owner, 1);
            function!(head, Minkowski {}.new_rep(4).to_symbolic([index]))
        };
        let ports = |first, second| {
            NumeratorInterface::admit(vector(first, 1) * vector(second, 2), None)
                .unwrap()
                .structure()
                .logical_slots()
        };
        // Symbolica orders product factors by the registration of their heads:
        // `a` precedes `b`, whereas `d` precedes `c`.
        let a = spenso::tensor_symbol!("feynkit_graph::port_order_a");
        let b = spenso::tensor_symbol!("feynkit_graph::port_order_b");
        let d = spenso::tensor_symbol!("feynkit_graph::port_order_d");
        let c = spenso::tensor_symbol!("feynkit_graph::port_order_c");
        assert_eq!(ports(a, b), ports(c, d));
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

    #[test]
    fn transport_keeps_canonical_orders_canonical_and_moves_other_orders_positionally() {
        let index = |owner: usize| function!(crate::symbols::hedge_index(), owner, 1);
        let port = |owner| Minkowski {}.new_rep(4).to_symbolic([index(owner)]);
        let matrix = spenso::tensor_symbol!("feynkit_graph::transported_matrix");
        let source = function!(matrix, port(1), port(2));
        let relabel = |map: [(usize, usize); 2]| {
            move |atom: &Atom| {
                Ok::<_, DiagramError>(atom.replace_multiple(
                    map.map(|(from, to)| Replacement::new(index(from).to_pattern(), index(to))),
                ))
            }
        };
        let declare = |owners: [usize; 2]| Some(NumeratorInterface(owners.map(port).to_vec()));

        // A canonical order, inferred or declared, stays canonical for the
        // swapped labels.
        let swap = relabel([(1, 2), (2, 1)]);
        let expected = NumeratorInterface::admit(swap(&source).unwrap(), None).unwrap();
        for interface in [None, declare([1, 2])] {
            let canonical = NumeratorInterface::admit(source.clone(), interface).unwrap();
            assert!(
                NumeratorInterface::translate(&canonical, swap)
                    .unwrap()
                    .same_value(&expected)
            );
        }

        // Any other order, which only a declaration produces, moves positionally.
        let shift = relabel([(1, 3), (2, 4)]);
        let declared = NumeratorInterface::admit(source.clone(), declare([2, 1])).unwrap();
        let expected = NumeratorInterface::admit(shift(&source).unwrap(), declare([4, 3])).unwrap();
        assert!(
            NumeratorInterface::translate(&declared, shift)
                .unwrap()
                .same_value(&expected)
        );
    }
}
