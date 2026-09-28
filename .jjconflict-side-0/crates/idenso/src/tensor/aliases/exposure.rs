//! Occurrence-local exposure of an already parsed alias definition.

use super::*;
use crate::tensor::{
    SymbolicNet,
    composition::{matching_interface_position, rewrite_interface_ports},
    inference::LeafInference,
};
use linnet::half_edge::involution::Hedge;
use spenso::{
    network::{
        TensorNetworkError, graph::NetworkEdge, parsing::ParseState,
        store::TensorScalarStoreMapping,
    },
    structure::{
        OrderedStructure, TensorStructure,
        representation::LibraryRep,
        slot::{DualSlotTo, Slot},
    },
};

impl SymbolicTensor<AliasInterfaces, AliasedAtom> {
    /// Expose one selected literal use without reparsing its definition.
    /// `encoded_body` is the cached logical boundary of `network`, before any
    /// occurrence relabelling. Internal incidences receive operation-fresh names;
    /// boundary ports retain the use's literal labels.
    ///
    /// The collector must first certify intrinsic normalization for every original
    /// definition in `definitions`, including bodies reached by nested literal
    /// registration, and reserve its complete input in `state`. Refused exposure
    /// leaves the literal opaque. No callback is speculatively executed here.
    /// `None` means this valid template cannot be exposed with unchanged topology
    /// (for example, metric normalization consumes a port). Invalid registry or
    /// layout evidence is an error. Observations are published only on success.
    pub(crate) fn expose_literal_graph<E: AtomCore + Clone + From<E::Output>>(
        network: &SymbolicNet<AbstractIndex, E>,
        definition: &Definition,
        encoded_body: &PartialStructure,
        current_use: &SymbolicTensor<PartialStructure>,
        state: &ParseState<AbstractIndex>,
        definitions: &mut HashMap<Atom, Definition>,
        observe: &mut impl FnMut(&Atom, &Atom),
    ) -> Result<Option<SymbolicNet<AbstractIndex, Atom>>> {
        let (handle, body) = definition;
        if network.graph.has_bound_ports() || body.expression.is_zero() {
            return Ok(None);
        }
        if !InterfaceInference::normalization_is_intrinsic(body.expression.as_view()) {
            return Ok(None);
        }
        let (AtomView::Fun(before), AtomView::Fun(after)) = (
            handle.expression.as_view(),
            current_use.expression.as_view(),
        ) else {
            return Err(TensorInferenceError::invalid(
                "exposure requires a tensor alias literal",
            ));
        };
        if (before.get_symbol() != spenso::tensor_symbol!("idenso::tensor_alias")
            && before.get_symbol() != symbolica::symbol!("idenso::scalar_alias"; Scalar))
            || after.get_symbol() != before.get_symbol()
            || after.get_nargs() != before.get_nargs()
        {
            return Err(TensorInferenceError::invalid(
                "alias exposure changed its head or arity",
            ));
        }
        let logical = body.structure.logical_slots();
        let handle_slots = handle.structure.logical_slots();
        let mut claimed = vec![false; logical.len()];
        let mut source_ports = vec![None; logical.len()];
        for slot in encoded_body.logical_slots() {
            let PartialIndex::Explicit(index) = slot.aind else {
                return Ok(None);
            };
            let source = slot.rep().slot::<AbstractIndex, _>(index);
            let position =
                matching_interface_position(source.to_atom().as_view(), &logical, &claimed)
                    .ok_or_else(|| {
                        TensorInferenceError::invalid("alias template has a different boundary")
                    })?;
            claimed[position] = true;
            source_ports[position] = Some(source);
        }
        if claimed.contains(&false) {
            return Err(TensorInferenceError::invalid(
                "alias template lost a boundary port",
            ));
        }
        claimed.fill(false);
        let mut handle_claimed = vec![false; handle_slots.len()];
        let mut replacements = Vec::new();
        for (old, new) in before.iter().zip(after.iter()) {
            let position = matching_interface_position(old, &logical, &claimed);
            let handle_position = matching_interface_position(old, &handle_slots, &handle_claimed);
            let (Some(position), Some(handle_position)) = (position, handle_position) else {
                if old != new {
                    return Err(TensorInferenceError::invalid(
                        "alias exposure changed non-port metadata",
                    ));
                }
                continue;
            };
            claimed[position] = true;
            handle_claimed[handle_position] = true;
            let Ok(target) = Slot::<LibraryRep, AbstractIndex>::try_from(new) else {
                // Compact vector bindings belong to the graph contraction owner.
                return Ok(None);
            };
            let source = source_ports[position].unwrap();
            if target.rep() != source.rep() || handle_slots[handle_position].rep() != source.rep() {
                return Err(TensorInferenceError::invalid(
                    "alias exposure changed a port representation",
                ));
            }
            state.reserve_index(target.aind());
            replacements.push((source, target));
        }
        if claimed.contains(&false) || handle_claimed.contains(&false) {
            return Err(TensorInferenceError::invalid(
                "alias handle lost a logical port",
            ));
        }
        let mut graph = network.map_ref(Clone::clone, Clone::clone);
        let mut boundary = Vec::new();
        for index in 0..graph.graph.graph.n_hedges() {
            let hedge = Hedge(index);
            if let NetworkEdge::Slot(slot) = graph.graph.graph[[&hedge]] {
                state.reserve_index(slot.aind());
                if graph.graph.graph.inv(hedge) == hedge {
                    let target = replacements
                        .iter()
                        .find_map(|(source, target)| {
                            // Sewing retains one descriptor for both endpoints.
                            // Keep that descriptor's variance while changing its
                            // index; map_occurrences restores each stored leaf's
                            // endpoint representation from its structural ordinal.
                            (*source == slot || source.matches(&slot))
                                .then(|| slot.rep().slot(target.aind()))
                        })
                        .ok_or_else(|| {
                            TensorInferenceError::invalid(
                                "parsed alias has an unrecorded boundary port",
                            )
                        })?;
                    boundary.push((hedge, target));
                }
            }
        }
        if boundary.len() != replacements.len() {
            return Err(TensorInferenceError::invalid(
                "parsed alias boundary is incomplete",
            ));
        }
        let mut visited = graph
            .graph
            .relabel_slot_components(&boundary)
            .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
        for index in 0..graph.graph.graph.n_hedges() {
            let hedge = Hedge(index);
            if visited.contains(&hedge) {
                continue;
            }
            if let NetworkEdge::Slot(slot) = graph.graph.graph[[&hedge]] {
                visited.extend(
                    graph
                        .graph
                        .relabel_slot_components(&[(hedge, slot.rep().slot(state.fresh_index()))])
                        .map_err(|error| TensorInferenceError::invalid(error.to_string()))?,
                );
            }
        }

        let mut unavailable = false;
        let mut literal_uses = Vec::new();
        let mapped = graph.map_occurrences(
            |scalar| Ok(scalar.as_atom_view().to_owned()),
            |tensor, ports, logical_order| {
                let fail = |message: &str| TensorNetworkError::Other(eyre::eyre!("{message}"));
                let Some(logical_order) = logical_order else {
                    return Err(fail("parsed alias occurrence has no logical layout"));
                };
                if ports.iter().any(|(_, _, bound)| bound.is_some())
                    || !InterfaceInference::normalization_is_intrinsic(tensor.expression.as_atom_view())
                {
                    unavailable = true;
                    return Err(fail("alias exposure requires an unbound intrinsic occurrence"));
                }
                let source = tensor.structure.external_structure();
                let source_interface = PartialStructure::from_logical_slots(logical_order.iter().map(|&position| {
                    let slot = source[position];
                    let index = match slot.aind() {
                        AbstractIndex::Open { axis, .. } => PartialIndex::open(axis),
                        index => PartialIndex::Explicit(index),
                    };
                    slot.rep().slot(index)
                }));
                let targets = ports.iter().map(|(_, slot, _)| *slot).collect::<Vec<_>>();
                let replacements = logical_order.iter().enumerate().map(|(logical, &storage)| {
                    (logical, targets[storage].to_atom())
                }).collect::<HashMap<_, _>>();
                let source_tensor = SymbolicTensor::from_normalized_parts(tensor.expression.as_atom_view().to_owned(), source_interface);
                let expression = rewrite_interface_ports(&source_tensor, &replacements, &mut |old, new| {
                    if matches!(old, AtomView::Fun(fun) if fun.get_symbol() == spenso::tensor_symbol!("idenso::tensor_alias")) {
                        literal_uses.push((old.to_owned(), new.clone()));
                    }
                }).map_err(|error| fail(&error.to_string()))?;
                let expected = PartialStructure::from_logical_slots(logical_order.iter().map(|&position| {
                    let slot = targets[position];
                    slot.rep().slot(PartialIndex::Explicit(slot.aind()))
                }));
                if SymbolicTensor::validate_observed_interface(
                    &expression,
                    &expected,
                    LeafInference::ObserveStorage,
                ).is_err() {
                    unavailable = true;
                    return Err(fail("normalized alias occurrence changed its interface"));
                }
                Ok(SymbolicTensor {
                    expression,
                    structure: OrderedStructure::new(targets).into_canonical(),
                    is_metric: tensor.is_metric,
                    is_composite: tensor.is_composite,
                    proofs: Default::default(),
                })
            },
        );
        let mapped = match mapped {
            Ok(mapped) => mapped,
            Err(_) if unavailable => return Ok(None),
            Err(error) => return Err(TensorInferenceError::invalid(error.to_string())),
        };
        for (old, new) in &literal_uses {
            Self::register_literal_use(old, new, definitions)?;
        }
        for (old, new) in &literal_uses {
            observe(old, new);
        }
        Ok(Some(mapped))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::{
        network::{library::DummyLibrary, parsing::ParseSettings},
        structure::representation::{Minkowski, RepName},
    };
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };

    type Tensor = SymbolicTensor<PartialStructure>;
    type Aliased = SymbolicTensor<AliasInterfaces, AliasedAtom>;

    fn parse(body: &Tensor) -> SymbolicNet<AbstractIndex, Atom> {
        SymbolicNet::try_from_view::<OrderedStructure, _>(
            body.expression.as_view(),
            &DummyLibrary::<SymbolicTensor<OrderedStructure>>::new(),
            &ParseSettings {
                precontract_scalars: false,
                parse_composite_scalars_as_tensors: true,
                ..ParseSettings::default()
            },
        )
        .unwrap()
    }

    fn leaf(symbol: symbolica::atom::Symbol, indices: &[usize]) -> Atom {
        FunctionBuilder::new(symbol)
            .add_args(indices.iter().map(|&index| {
                Minkowski {}
                    .new_slot::<AbstractIndex, _, _>(4, AbstractIndex::Normal(index))
                    .to_lib()
                    .to_atom()
            }))
            .finish()
    }

    fn internal_indices(network: &SymbolicNet<AbstractIndex, Atom>) -> HashSet<AbstractIndex> {
        let boundary = network.graph.dangling_indices();
        (0..network.graph.graph.n_hedges())
            .filter_map(|i| match network.graph.graph[[&Hedge(i)]] {
                NetworkEdge::Slot(slot) if !boundary.contains(&slot) => Some(slot.aind()),
                _ => None,
            })
            .collect()
    }

    #[test]
    fn exposed_occurrences_freshen_internal_incidence_without_renaming_boundary() {
        crate::test_support::test_initialize();
        let p = spenso::tensor_symbol!("exposure_scope_p");
        let q = spenso::tensor_symbol!("exposure_scope_q");
        let body = Tensor::infer(leaf(p, &[811, 812]) * leaf(q, &[812])).unwrap();
        let handle = body.alias_handle().unwrap();
        let current = handle
            .reindex_interface_ports(&HashMap::from([(0, AbstractIndex::Normal(812))]))
            .unwrap();
        let definition = (handle.clone(), body.clone());
        let network = parse(&body);
        let original = format!("{:?}", network.graph);
        let encoded = InterfaceInference::replacement_interface(body.expression.as_view()).unwrap();
        let state = ParseState::default();
        state.reserve_indices(body.expression.as_view());
        let mut definitions = HashMap::from([(handle.expression.clone(), definition.clone())]);
        let mut previous = HashSet::new();
        for _ in 0..2 {
            let exposed = Aliased::expose_literal_graph(
                &network,
                &definition,
                &encoded,
                &current,
                &state,
                &mut definitions,
                &mut |_, _| panic!("no nested alias"),
            )
            .unwrap()
            .unwrap();
            assert_eq!(
                exposed.graph.dangling_indices(),
                vec![
                    Slot::<LibraryRep, AbstractIndex>::try_from(spenso::mink!(4, 812).as_view())
                        .unwrap()
                ]
            );
            let internal = internal_indices(&exposed);
            assert_eq!(internal.len(), 1);
            assert!(!internal.contains(&AbstractIndex::Normal(811)));
            assert!(!internal.contains(&AbstractIndex::Normal(812)));
            assert!(previous.is_disjoint(&internal));
            previous.extend(internal);
            for tensor in &exposed.store.tensors {
                tensor
                    .validate_rewritten_interface(&tensor.expression)
                    .unwrap();
            }
            assert_eq!(exposed.graph.graph.n_nodes(), network.graph.graph.n_nodes());
            assert_eq!(
                exposed.graph.graph.n_hedges(),
                network.graph.graph.n_hedges()
            );
        }
        assert_eq!(format!("{:?}", network.graph), original);
        assert_eq!(definitions.len(), 1);
    }

    #[test]
    fn exposed_scalar_aliases_and_sum_seams_keep_copy_local_contractions() {
        crate::test_support::test_initialize();
        let p = spenso::tensor_symbol!("exposure_sum_p");
        let q = spenso::tensor_symbol!("exposure_sum_q");
        let r = spenso::tensor_symbol!("exposure_sum_r");
        // Both alternatives share one external sum seam before meeting q.
        let body = Tensor::infer((leaf(p, &[821]) + leaf(r, &[821])) * leaf(q, &[821])).unwrap();
        let handle = body.alias_handle().unwrap();
        assert!(handle.is_scalar());
        let definition = (handle.clone(), body.clone());
        let network = parse(&body);
        let encoded = InterfaceInference::replacement_interface(body.expression.as_view()).unwrap();
        let mut definitions = HashMap::from([(handle.expression.clone(), definition.clone())]);
        let exposed = Aliased::expose_literal_graph(
            &network,
            &definition,
            &encoded,
            &handle,
            &ParseState::default(),
            &mut definitions,
            &mut |_, _| panic!("no nested alias"),
        )
        .unwrap()
        .unwrap();
        assert!(exposed.graph.dangling_indices().is_empty());
        assert_eq!(
            internal_indices(&exposed).len(),
            1,
            "sum alternatives use the same sewn incidence"
        );
        assert!(!internal_indices(&exposed).contains(&AbstractIndex::Normal(821)));
        for tensor in &exposed.store.tensors {
            tensor
                .validate_rewritten_interface(&tensor.expression)
                .unwrap();
        }
    }

    #[test]
    fn exposed_nested_literals_are_registered_before_observation() {
        crate::test_support::test_initialize();
        let p = spenso::tensor_symbol!("exposure_nested_p");
        let q = spenso::tensor_symbol!("exposure_nested_q");
        let inner = Tensor::infer(leaf(p, &[831, 832])).unwrap();
        let inner_handle = inner.alias_handle().unwrap();
        let body = Tensor::infer(&inner_handle.expression * leaf(q, &[832])).unwrap();
        let handle = body.alias_handle().unwrap();
        let definition = (handle.clone(), body.clone());
        let network = parse(&body);
        let encoded = InterfaceInference::replacement_interface(body.expression.as_view()).unwrap();
        let mut definitions = HashMap::from([
            (handle.expression.clone(), definition.clone()),
            (
                inner_handle.expression.clone(),
                (inner_handle.clone(), inner),
            ),
        ]);
        let mut observations = Vec::new();
        let exposed = Aliased::expose_literal_graph(
            &network,
            &definition,
            &encoded,
            &handle,
            &ParseState::default(),
            &mut definitions,
            &mut |old, new| observations.push((old.clone(), new.clone())),
        )
        .unwrap()
        .unwrap();
        let [(old, new)] = observations.as_slice() else {
            panic!("one rewritten nested literal");
        };
        assert_eq!(old, &inner_handle.expression);
        assert_ne!(old, new);
        let (next_handle, next_body) = &definitions[new];
        Tensor::validate_encoded_interface(&next_body.expression, &next_handle.structure).unwrap();
        assert!(
            exposed
                .store
                .tensors
                .iter()
                .any(|tensor| tensor.expression == new)
        );
    }

    #[test]
    fn exposure_declines_intrinsic_rank_collapse_and_callback_templates() {
        crate::test_support::test_initialize();
        let body = Tensor::infer(leaf(
            spenso::network::library::symbolic::ETS.metric,
            &[841, 842],
        ))
        .unwrap();
        let handle = body.alias_handle().unwrap();
        let AtomView::Fun(function) = handle.expression.as_view() else {
            unreachable!()
        };
        let target = Tensor::infer(
            FunctionBuilder::new(function.get_symbol())
                .add_arg(function.iter().next().unwrap())
                .add_arg(spenso::mink!(4, 841))
                .add_arg(spenso::mink!(4, 841))
                .finish(),
        )
        .unwrap();
        let definition = (handle.clone(), body.clone());
        let mut definitions = HashMap::from([(handle.expression.clone(), definition.clone())]);
        assert!(
            Aliased::expose_literal_graph(
                &parse(&body),
                &definition,
                &InterfaceInference::replacement_interface(body.expression.as_view()).unwrap(),
                &target,
                &ParseState::default(),
                &mut definitions,
                &mut |_, _| panic!("refusal must publish nothing")
            )
            .unwrap()
            .is_none()
        );
        assert_eq!(definitions.len(), 1);

        let calls = Arc::new(AtomicUsize::new(0));
        let captured = Arc::clone(&calls);
        let callback = spenso::tensor_symbol!(
            "exposure_callback",
            norm = move |_, _| {
                captured.fetch_add(1, Ordering::Relaxed);
            }
        );
        let body = Tensor::infer(leaf(callback, &[843])).unwrap();
        let handle = body.alias_handle().unwrap();
        let definition = (handle.clone(), body.clone());
        let network = parse(&body);
        let encoded = InterfaceInference::replacement_interface(body.expression.as_view()).unwrap();
        calls.store(0, Ordering::Relaxed);
        let mut definitions = HashMap::from([(handle.expression.clone(), definition.clone())]);
        assert!(
            Aliased::expose_literal_graph(
                &network,
                &definition,
                &encoded,
                &handle,
                &ParseState::default(),
                &mut definitions,
                &mut |_, _| panic!("refusal must publish nothing")
            )
            .unwrap()
            .is_none()
        );
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn exposure_preserves_explicit_identity_and_mixed_auto_axis_order() {
        crate::test_support::test_initialize();
        use spenso::structure::representation::Euclidean;
        let head = spenso::tensor_symbol!("exposure_axis_matrix");
        let slot = |index| {
            Minkowski {}
                .new_slot::<AbstractIndex, _, _>(4, index)
                .to_lib()
        };
        let a = slot(AbstractIndex::Normal(851));
        let b = slot(AbstractIndex::Normal(852));
        let expression = FunctionBuilder::new(head)
            .add_arg(a.to_atom())
            .add_arg(b.to_atom())
            .finish();
        let body = Tensor::checked_parts(
            expression,
            PartialStructure::from_logical_slots(
                [b, a].map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
            ),
        )
        .unwrap();
        let handle = body.alias_handle().unwrap();
        let current = handle
            .reindex_interface_ports(&HashMap::from([
                (0, AbstractIndex::Normal(862)),
                (1, AbstractIndex::Normal(861)),
            ]))
            .unwrap();
        let expected = FunctionBuilder::new(head)
            .add_arg(slot(AbstractIndex::Normal(861)).to_atom())
            .add_arg(slot(AbstractIndex::Normal(862)).to_atom())
            .finish();
        let open = |axis| {
            Euclidean {}
                .new_slot::<AbstractIndex, _, _>(2, AbstractIndex::Open { owner: 853, axis })
                .to_lib()
                .to_atom()
        };
        let auto = Tensor::infer(
            FunctionBuilder::new(head)
                .add_arg(open(0))
                .add_arg(a.to_atom())
                .add_arg(open(1))
                .finish(),
        )
        .unwrap();
        let auto_handle = auto.alias_handle().unwrap();
        let auto_current = auto_handle
            .reindex_interface_ports(&HashMap::from([
                (0, AbstractIndex::Normal(871)),
                (2, AbstractIndex::Normal(872)),
            ]))
            .unwrap();
        let auto_expected = FunctionBuilder::new(head)
            .add_arg(
                Euclidean {}
                    .new_slot::<AbstractIndex, _, _>(2, AbstractIndex::Normal(871))
                    .to_lib()
                    .to_atom(),
            )
            .add_arg(a.to_atom())
            .add_arg(
                Euclidean {}
                    .new_slot::<AbstractIndex, _, _>(2, AbstractIndex::Normal(872))
                    .to_lib()
                    .to_atom(),
            )
            .finish();
        for (body, handle, current, expected) in [
            (body, handle, current, expected),
            (auto, auto_handle, auto_current, auto_expected),
        ] {
            let definition = (handle.clone(), body.clone());
            let mut definitions = HashMap::from([(handle.expression.clone(), definition.clone())]);
            let exposed = Aliased::expose_literal_graph(
                &parse(&body),
                &definition,
                &InterfaceInference::replacement_interface(body.expression.as_view()).unwrap(),
                &current,
                &ParseState::default(),
                &mut definitions,
                &mut |_, _| panic!("no nested alias"),
            )
            .unwrap()
            .unwrap();
            assert_eq!(exposed.store.tensors.len(), 1);
            assert_eq!(exposed.store.tensors[0].expression, expected);
        }
    }
    #[test]
    fn extracted_dual_endpoints_keep_their_stored_representation() {
        use linnet::half_edge::subgraph::{ModifySubSet, SuBitGraph};
        use spenso::network::graph::{NetworkLeaf, NetworkNode};
        use spenso::structure::slot::IsAbstractSlot;

        crate::test_support::test_initialize();
        let rep = LibraryRep::from(crate::representations::ColorFundamental {}).new_rep(3);
        let index = AbstractIndex::Normal(881);
        let target = AbstractIndex::Normal(882);
        let heads = [
            spenso::tensor_symbol!("exposure_endpoint_fundamental"),
            spenso::tensor_symbol!("exposure_endpoint_antifundamental"),
        ];
        let slots = [rep.slot::<AbstractIndex, _>(index), rep.dual().slot(index)];
        let bodies = heads
            .into_iter()
            .zip(slots)
            .map(|(head, slot)| {
                Tensor::infer(FunctionBuilder::new(head).add_arg(slot.to_atom()).finish()).unwrap()
            })
            .collect::<Vec<_>>();
        let input = Tensor::infer(&bodies[0].expression * &bodies[1].expression).unwrap();
        let network = parse(&input);
        let original = format!("{:?}", network.graph);
        let mut reversed = 0;
        for (body, slot) in bodies.into_iter().zip(slots) {
            let node = network
                .graph
                .graph
                .iter_nodes()
                .find_map(|(node, _, value)| match value {
                    NetworkNode::Leaf(NetworkLeaf::LocalTensor(index))
                        if network.store.tensors[*index].expression == body.expression =>
                    {
                        Some(node)
                    }
                    _ => None,
                })
                .unwrap();
            let mut subset: SuBitGraph = network.graph.graph.empty_subgraph();
            for hedge in network.graph.graph.iter_crown(node) {
                subset.add(hedge);
            }
            let mut template = network.map_ref(Clone::clone, Clone::clone);
            template.graph = template.graph.extract(&subset);
            let stored = template.graph.dangling_indices();
            assert_eq!(stored.len(), 1);
            assert!(stored[0] == slot || stored[0].matches(&slot));
            reversed += usize::from(stored[0] != slot);

            let handle = body.alias_handle().unwrap();
            let current = handle
                .reindex_interface_ports(&HashMap::from([(0, target)]))
                .unwrap();
            let encoded =
                InterfaceInference::replacement_interface(body.expression.as_view()).unwrap();
            let definition = (handle.clone(), body.clone());
            let mut definitions = HashMap::from([(handle.expression.clone(), definition.clone())]);
            let state = ParseState::default();
            state.reserve_indices(input.expression.as_view());
            let exposed = Aliased::expose_literal_graph(
                &template,
                &definition,
                &encoded,
                &current,
                &state,
                &mut definitions,
                &mut |_, _| panic!("no nested alias"),
            )
            .unwrap()
            .unwrap();
            let AtomView::Fun(function) = body.expression.as_view() else {
                unreachable!()
            };
            let expected = FunctionBuilder::new(function.get_symbol())
                .add_arg(slot.rep().slot::<AbstractIndex, _>(target).to_atom())
                .finish();
            assert_eq!(exposed.store.tensors.len(), 1);
            let tensor = &exposed.store.tensors[0];
            assert_eq!(tensor.expression, expected);
            assert_eq!(
                tensor.structure.external_structure(),
                vec![slot.rep().slot(target)]
            );
            tensor
                .validate_rewritten_interface(&tensor.expression)
                .unwrap();
            assert_eq!(definitions.len(), 1);
        }
        assert_eq!(
            reversed, 1,
            "exactly one endpoint inherits its dual's shared descriptor"
        );
        assert_eq!(format!("{:?}", network.graph), original);
    }
    #[test]
    fn exposed_closed_trace_retains_its_stored_local_incidences() {
        crate::test_support::test_initialize();
        let body = Tensor::infer(symbolica::parse!(
            "trace(cof(2),cyclic(t(coad(3,trace_axis),in,out),t(coad(3,trace_axis),in,out)))",
            default_namespace = "spenso"
        ))
        .unwrap();
        assert!(body.is_scalar());
        let handle = body.alias_handle().unwrap();
        let definition = (handle.clone(), body.clone());
        let network = SymbolicNet::try_from_view::<OrderedStructure, _>(
            body.expression.as_view(),
            &DummyLibrary::<SymbolicTensor<OrderedStructure>>::new(),
            &ParseSettings {
                precontract_scalars: false,
                parse_composite_scalars_as_tensors: true,
                shorthand_parsing: spenso::network::parsing::ShorthandParsing::Opaque,
                ..ParseSettings::default()
            },
        )
        .unwrap();
        let encoded = InterfaceInference::replacement_interface(body.expression.as_view()).unwrap();
        let mut definitions = HashMap::from([(handle.expression.clone(), definition.clone())]);
        let exposed = Aliased::expose_literal_graph(
            &network,
            &definition,
            &encoded,
            &handle,
            &ParseState::default(),
            &mut definitions,
            &mut |_, _| panic!("no nested alias"),
        )
        .unwrap()
        .unwrap();
        assert!(exposed.graph.dangling_indices().is_empty());
        assert_eq!(exposed.graph.graph.n_nodes(), network.graph.graph.n_nodes());
        assert_eq!(
            exposed.graph.graph.n_hedges(),
            network.graph.graph.n_hedges()
        );
        assert_eq!(exposed.store.tensors.len(), network.store.tensors.len());
        let trace_position = network
            .store
            .tensors
            .iter()
            .position(|tensor| {
                tensor.expression == body.expression
                    && tensor.structure.external_structure().len() == 2
            })
            .expect("the parsed trace retains both locally sewn axes");
        assert_eq!(
            exposed.store.tensors[trace_position]
                .structure
                .external_structure()
                .len(),
            2
        );
        for (original, rewritten) in network.store.tensors.iter().zip(&exposed.store.tensors) {
            assert_eq!(
                original.structure.external_structure().len(),
                rewritten.structure.external_structure().len()
            );
        }
        assert!(internal_indices(&network).is_disjoint(&internal_indices(&exposed)));
    }
}
