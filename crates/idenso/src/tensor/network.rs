//! Checked lowering from a symbolic carrier to finite component execution.

use super::TensorProofs;
use super::composition;
use super::inference::{InterfaceInference, TensorInferenceError};
use super::{SymbolicNet, SymbolicNetParse, SymbolicTensor};
use crate::CookSettings;
use spenso::network::library::{DummyKey, DummyLibrary, FunctionLibrary};
use spenso::network::parsing::{
    LeafInterfaceCache, ParseSettings, StrictTensorFilter, StructureFromAtom, TensorFromExpression,
    TensorLibraryFor,
};
use spenso::network::store::NetworkStore;
use spenso::network::tags::SPENSO_TAG;
use spenso::network::{Network, TensorNetworkError};
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialIndex, PartialStructure, PartialStructureExt};
use spenso::structure::representation::{LibraryRep, LibrarySlot};
use spenso::structure::slot::{
    AbsInd, DummyAind, IsAbstractSlot, ParseableAind, SlotMatch, SlotMatcher,
};
use spenso::structure::{
    Canonicalized, HasStructure, IndexedStructure, OrderedStructure, ScalarStructure,
    TensorStructure,
};
use spenso::tensors::data::SparseTensor;
use std::fmt::{Debug, Display};
use std::sync::Arc;
use symbolica::atom::{Atom, AtomView, Symbol};

impl TensorProofs {
    pub(crate) fn parser_leaf_interfaces(&self) -> LeafInterfaceCache {
        Arc::clone(self.leaf_interfaces.get_or_init(Default::default))
    }

    /// Opaque arithmetic and shorthand regions expose their retained free
    /// boundary. Ordinary functions keep the separate raw-port witness, which
    /// includes repeated internal indices and their written order.
    fn parser_interfaces(
        &self,
        expression: AtomView<'_>,
        decode_port: impl Fn(&Atom) -> Atom,
    ) -> Result<LeafInterfaceCache, TensorInferenceError> {
        let cache = Arc::clone(self.leaf_interfaces.get_or_init(Default::default));
        let mut pending = vec![expression];
        let mut slots = SlotMatcher::with_leaf_interfaces([Arc::clone(&cache)]);
        while let Some(value) = pending.pop() {
            if !matches!(slots.classify(value), SlotMatch::Other) {
                continue;
            }
            let composite = match value {
                AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_) => true,
                AtomView::Fun(function) => {
                    let head = function.get_symbol();
                    if head.is_scalar() {
                        continue;
                    }
                    composition::is_composite_head(head, &SPENSO_TAG)
                        || head == SPENSO_TAG.dot
                        || InterfaceInference::is_compound_metric(function, &mut slots)
                }
                _ => false,
            };
            if composite {
                let retained = cache
                    .lock()
                    .expect("leaf interface cache is not poisoned")
                    .get(value.get_data())
                    .is_some_and(|(tensorial, _)| *tensorial);
                if !retained {
                    let interface = self.admitted_interface(value).ok_or_else(|| {
                        TensorInferenceError::invalid(format!(
                            "component lowering is missing a retained regional interface for {value}"
                        ))
                    })?;
                    let ports = interface
                        .logical_slots()
                        .into_iter()
                        .map(|slot| decode_port(&composition::port_atom(slot)))
                        .collect();
                    cache
                        .lock()
                        .expect("leaf interface cache is not poisoned")
                        .insert(value.to_owned(), (true, ports));
                }
            }
            match value {
                AtomView::Fun(function) => pending.extend(function.iter()),
                AtomView::Add(sum) => pending.extend(sum.iter()),
                AtomView::Mul(product) => pending.extend(product.iter()),
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    pending.extend([base, exponent]);
                }
                _ => {}
            }
        }
        Ok(cache)
    }

    /// Publish the ports selected by a tensor constructor before its expression
    /// loses provenance. Raw slots retain internal pairs for graph lowering;
    /// the region boundary contracts those pairs without inspecting the payload.
    pub(crate) fn retain_leaf_ports(
        &self,
        expression: AtomView<'_>,
        ports: Vec<Atom>,
    ) -> Result<(), TensorInferenceError> {
        self.retain_boundary(
            expression,
            &Self::boundary_from_leaf_ports(ports.iter().map(Atom::as_view))?,
        );
        if matches!(expression, AtomView::Fun(_)) {
            self.leaf_interfaces
                .get_or_init(Default::default)
                .lock()
                .expect("leaf interface cache is not poisoned")
                .insert(expression.to_owned(), (true, ports));
        }
        Ok(())
    }

    pub(super) fn boundary_from_leaf_ports<'a>(
        ports: impl IntoIterator<Item = AtomView<'a>>,
    ) -> Result<PartialStructure, TensorInferenceError> {
        let mut matcher = SlotMatcher::default();
        let interface = ports
            .into_iter()
            .enumerate()
            .filter(|(_, port)| !matches!(port, AtomView::Var(var)
                if var.get_symbol() == SPENSO_TAG.chain_in || var.get_symbol() == SPENSO_TAG.chain_out))
            .map(|(position, port)| {
                if let Ok(slot) = matcher.parse::<LibraryRep, AbstractIndex>(port) {
                    Ok(slot.rep().slot(PartialIndex::Explicit(slot.aind())))
                } else {
                    matcher
                        .parse_representation::<LibraryRep>(port)
                        .map(|rep| rep.slot(PartialIndex::open(position)))
                        .map_err(|reason| TensorInferenceError::invalid(reason.to_string()))
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        InterfaceInference::merge_explicit_interface_sequence(&[
            PartialStructure::from_logical_slots(interface),
        ])
    }
}

impl SymbolicTensor<PartialStructure> {
    /// Isolate facts emitted by one producer. Unchanged expression keys remain
    /// shared, while new records cannot accumulate on a retained source tensor.
    pub(crate) fn producer_context(&self) -> (TensorProofs, LeafInterfaceCache) {
        let mut proofs = self.proofs.clone();
        proofs.publish_interfaces(self.expression.as_view(), self.structure.clone(), []);
        let leaves = Arc::clone(proofs.leaf_interfaces.get_or_init(Default::default));
        (proofs, leaves)
    }

    /// Lower this admitted tensor into a component network, retaining its boundary.
    ///
    /// The parser still owns shorthand topology, fresh dummy allocation and library
    /// realization. Optional index decoding is confined to this lowering boundary;
    /// the symbolic carrier keeps its representable indices and validation facts.
    /// Open ports are materialized by the existing positional composition owner.
    /// A ranked zero becomes an empty component tensor with that same interface.
    /// The returned layout maps canonical component axes back to the declared
    /// logical axes and remains valid throughout contraction of internal ports.
    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn into_network<T, S, K, Aind, Lib, FunLib>(
        self,
        library: &Lib,
        functions: &FunLib,
        settings: &ParseSettings,
        index_decoding: Option<&CookSettings>,
    ) -> Result<
        Canonicalized<Network<NetworkStore<T, Atom>, K, Symbol, Aind>>,
        TensorNetworkError<K, Symbol>,
    >
    where
        Aind: AbsInd + DummyAind + ParseableAind,
        K: Clone + Debug + Display,
        S: TensorStructure<Slot = LibrarySlot<Aind>>
            + ScalarStructure
            + Clone
            + StructureFromAtom
            + From<OrderedStructure<LibraryRep, Aind>>,
        T: HasStructure<Structure = S>
            + TensorStructure<Slot = LibrarySlot<Aind>>
            + Clone
            + From<SparseTensor<Atom, S>>
            + for<'src> TensorFromExpression<'src, S, Atom, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Atom, Key = Symbol>,
    {
        let invalid = |message: String| TensorNetworkError::Other(eyre::eyre!(message));
        if index_decoding.is_some_and(|settings| {
            settings.mode() != crate::CookMode::ReversibleEncoding
                || !matches!(
                    settings.source_filter(),
                    crate::CookSourceFilter::RepresentationPayload {
                        indices: true,
                        dimensions: false,
                        ..
                    }
                )
        }) {
            return Err(invalid(
                "component lowering only decodes reversibly encoded index payloads".into(),
            ));
        }
        self.explicit_index_occurrences()
            .map_err(|error| invalid(error.to_string()))?;
        let mut source = if self.structure.open_positions().is_empty() {
            self
        } else {
            let source = self
                .materialized()
                .map_err(|error| invalid(error.to_string()))?;
            drop(self);
            source
        };
        // Component construction needs admitted facts, not the symbolic
        // execution plan. Release that graph before decoding can run callbacks.
        source.replace_admission(source.producer_context().0);
        let decode = |expression: &Atom| match index_decoding {
            Some(settings) => settings.uncook(expression.as_view()),
            None => expression.clone(),
        };
        // The mapper runs each rebuilt function's normalizer once and carries
        // exact raw-port witnesses alongside it. Temporary regional facts keep
        // the source labels until their boundary ports are decoded below;
        // they are never published as a symbolic tensor with stale indices.
        let (expression, leaf_interfaces) = if let Some(decoding) = index_decoding {
            let (expression, mapped) = source
                .map_admitted_expression(
                    |value| {
                        if matches!(value, AtomView::Fun(function)
                            if function.get_symbol().is_scalar())
                        {
                            return Ok(Some(value.to_owned()));
                        }
                        if matches!(value, AtomView::Var(_)) {
                            let decoded = decoding.uncook(value);
                            return Ok((decoded.as_view() != value).then_some(decoded));
                        }
                        Ok(None)
                    },
                    |interface| Ok(interface.clone()),
                )
                .map_err(|error| invalid(error.to_string()))?;
            let leaves = mapped
                .parser_interfaces(expression.as_view(), decode)
                .map_err(|error| invalid(error.to_string()))?;
            (Some(expression), leaves)
        } else {
            let leaves = source
                .proofs
                .parser_interfaces(source.expression.as_view(), decode)
                .map_err(|error| invalid(error.to_string()))?;
            (None, leaves)
        };
        let (original, structure) = source.into_parts();
        let expression = expression.unwrap_or(original);
        let expected = structure
            .logical_slots()
            .into_iter()
            .map(|slot| {
                let PartialIndex::Explicit(index) = slot.aind else {
                    return Err(invalid("component lowering left an unresolved port".into()));
                };
                let encoded = decode(&index.to_atom());
                let index = Aind::from_view(encoded.as_view())
                    .map_err(|_| invalid(format!("component index cannot represent {encoded}")))?;
                Ok(slot.rep().slot(index))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let interface: Canonicalized<OrderedStructure<LibraryRep, Aind>> =
            expected.iter().copied().collect();

        if expression.is_zero() && !expected.is_empty() {
            let tensor = SparseTensor::empty(S::from(interface.canonical().clone()), Atom::Zero);
            return Ok(interface.map_canonical(|_| Network::from_tensor(T::from(tensor))));
        }

        let network =
            Network::<NetworkStore<T, Atom>, K, Symbol, Aind>::try_from_admitted_view_with_function_library(
                expression.as_view(),
                library,
                functions,
                settings,
                leaf_interfaces,
                None,
            )?;
        let mut actual = network.graph.dangling_indices();
        for slot in expected {
            let Some(position) = actual.iter().position(|candidate| *candidate == slot) else {
                return Err(invalid(
                    "component lowering changed the admitted tensor interface".into(),
                ));
            };
            actual.swap_remove(position);
        }
        if !actual.is_empty() {
            return Err(invalid(
                "component lowering introduced an unadmitted tensor port".into(),
            ));
        }
        Ok(interface.map_canonical(|_| network))
    }
}

impl SymbolicNetParse for SymbolicTensor<PartialStructure> {
    fn parse_to_symbolic_net<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<SymbolicNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        let materialized;
        let source = if self.structure.open_positions().is_empty() {
            self
        } else {
            materialized = self
                .materialized()
                .map_err(|error| TensorNetworkError::Other(eyre::eyre!(error.to_string())))?;
            &materialized
        };
        let slots = source.structure.logical_slots().into_iter().map(|slot| {
            let PartialIndex::Explicit(index) = slot.aind else {
                unreachable!("materialized symbolic-network ports are explicit")
            };
            slot.rep().slot(index)
        });
        source.parse_admitted_symbolic_net(settings, slots)
    }
}

impl SymbolicNetParse for SymbolicTensor<IndexedStructure> {
    fn parse_to_symbolic_net<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<SymbolicNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        self.parse_admitted_symbolic_net(settings, self.structure.logical_slots())
    }
}

impl<S> SymbolicTensor<S> {
    fn parse_admitted_symbolic_net<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        settings: &ParseSettings,
        slots: impl IntoIterator<Item = LibrarySlot<AbstractIndex>>,
    ) -> Result<SymbolicNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        let slots = slots.into_iter().collect::<Vec<_>>();
        if self.expression.is_zero() {
            let slots = slots
                .into_iter()
                .map(|slot| {
                    let index = slot.aind().to_atom();
                    Aind::from_view(index.as_view())
                        .map(|index| slot.rep().slot(index))
                        .map_err(|_| {
                            TensorNetworkError::Other(eyre::eyre!(
                                "symbolic network index cannot represent {index}"
                            ))
                        })
                })
                .collect::<Result<Vec<_>, _>>()?;
            return Ok(SymbolicNet::from_tensor(SymbolicTensor {
                expression: Atom::Zero,
                structure: OrderedStructure::new(slots).into_canonical(),
                work: self.work.clone(),
                proofs: self.proofs.clone(),
                is_metric: false,
                is_composite: false,
            }));
        }
        let library = DummyLibrary::<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>::new();
        let settings = settings
            .clone()
            .with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
        let mut proofs = self.proofs.clone();
        let interface = PartialStructure::from_logical_slots(
            slots
                .into_iter()
                .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
        );
        proofs.publish_interfaces(self.expression.as_view(), interface, []);
        let leaves = proofs.parser_leaf_interfaces();
        SymbolicNet::try_from_admitted_view_with_context::<OrderedStructure<LibraryRep, Aind>, _>(
            self.expression.as_view(),
            &library,
            &settings,
            leaves,
            Some(&proofs),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tensor::inference::tests::{INFERENCE_CALLS, INFERENCE_VISITS};
    use crate::tensor::planning::PlanningNetwork;
    use crate::{CookMode, CookSettings, test_support};
    use spenso::iterators::IteratableTensor;
    use spenso::network::library::panicing::ErroringLibrary;
    use spenso::network::library::{DummyKey, DummyLibrary};
    use spenso::network::parsing::{ShadowedStructure, ShorthandParsing};
    use spenso::network::{ExecutionResult, Sequential, SmallestDegree};
    use spenso::structure::abstract_index::AbstractIndex;
    use spenso::structure::representation::{Euclidean, RepName};
    use spenso::tensors::parametric::ParamTensor;
    use spenso::{index_symbol, tensor_symbol};
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::sync::{Arc, Mutex, Weak};
    use symbolica::atom::{AtomCore, AtomView};
    use symbolica::function;

    type Structure = ShadowedStructure<AbstractIndex>;
    type Tensor = ParamTensor<Structure>;
    type Net = Network<NetworkStore<Tensor, Atom>, DummyKey, Symbol, AbstractIndex>;

    #[test]
    fn typed_symbolic_network_retains_leaf_loops_and_ranked_zero() {
        test_support::test_initialize();
        let rep = Euclidean {}.new_rep(2);
        let slot = rep.slot::<AbstractIndex, _>(83).to_atom();
        let source = SymbolicTensor::infer(function!(
            tensor_symbol!("typed_symbolic_network_loop"),
            &slot,
            &slot
        ))
        .unwrap();
        let typed = source
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        let raw = source
            .expression()
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert!(typed.graph.dangling_indices().is_empty());
        assert_eq!(
            typed.graph.graph.iter_edges().count(),
            raw.graph.graph.iter_edges().count(),
        );

        let vector = SymbolicTensor::infer(function!(
            tensor_symbol!("typed_symbolic_network_vector"),
            &slot
        ))
        .unwrap();
        let zero = vector.with_produced_expression(Atom::Zero).unwrap();
        let network = zero
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(network.graph.dangling_indices().len(), 1);
        assert_eq!(network.graph.dangling_indices()[0].to_atom(), slot);
    }

    #[test]
    fn typed_component_lowering_matches_direct_component_execution() {
        crate::test_support::test_initialize();
        let rep = Euclidean {}.new_rep(2);
        let mu = rep.slot::<AbstractIndex, _>(71).to_atom();
        let nu = rep.slot::<AbstractIndex, _>(73).to_atom();
        let left = spenso::tensor_symbol!("typed_component_lowering_left");
        let right = spenso::tensor_symbol!("typed_component_lowering_right");
        let expression = symbolica::parse!("(typed_component_x+typed_component_y)^7")
            * function!(left, &mu)
            * spenso::network::library::symbolic::ETS.metric_literal(mu, nu.clone())
            * function!(right, nu);
        let source = SymbolicTensor::infer(expression.clone()).unwrap();
        let library = DummyLibrary::<Tensor>::new();
        let functions = ErroringLibrary::<Symbol>::new();
        let settings = ParseSettings {
            shorthand_parsing: ShorthandParsing::expand_all(),
            ..Default::default()
        };
        let mut typed: Net = source
            .clone()
            .into_network(&library, &functions, &settings, None)
            .unwrap()
            .into_canonical();
        let mut direct = Net::try_from_view_with_function_library(
            expression.as_view(),
            &library,
            &functions,
            &settings,
        )
        .unwrap();
        typed
            .execute::<Sequential, SmallestDegree, _, _, _>(&library, &functions)
            .unwrap();
        direct
            .execute::<Sequential, SmallestDegree, _, _, _>(&library, &functions)
            .unwrap();
        let typed: Atom = typed.result_scalar().unwrap().into();
        let direct: Atom = direct.result_scalar().unwrap().into();
        assert_eq!(typed, direct);
        assert_eq!(source.expression(), &expression);
    }

    #[test]
    fn component_lowering_transports_decoded_loops_for_every_filter() {
        test_support::test_initialize();
        let decoding = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let label = function!(index_symbol!("component_decoded_loop_index"), 18, 2);
        let encoded = CookSettings::reversible()
            .cook_function(label.as_fun_view().unwrap())
            .unwrap();
        let rep = Euclidean {}.new_rep(2);
        let slot = rep.to_symbolic([encoded]);
        let head = tensor_symbol!("component_decoded_loop");
        let expression = function!(head, &slot, &slot);
        let source = SymbolicTensor::infer(expression.clone()).unwrap();
        let decoded = decoding.uncook(expression.as_view());
        let library = DummyLibrary::<Tensor>::new();
        let functions = ErroringLibrary::<Symbol>::new();
        for filter in [
            StrictTensorFilter::Tagged,
            StrictTensorFilter::TaggedChecked,
            StrictTensorFilter::ContainsReps,
        ] {
            for composite_scalars in [false, true] {
                let settings = ParseSettings {
                    strict_tensor_filter: filter,
                    parse_composite_scalars_as_tensors: composite_scalars,
                    ..Default::default()
                };
                let mut expected = Net::try_from_view_with_function_library(
                    decoded.as_view(),
                    &library,
                    &functions,
                    &settings,
                )
                .unwrap();
                expected
                    .execute::<Sequential, SmallestDegree, _, _, _>(&library, &functions)
                    .unwrap();
                let expected: Atom = expected.result_scalar().unwrap().into();
                INFERENCE_CALLS.with(|count| count.set(0));
                INFERENCE_VISITS.with(|count| count.set(0));
                let mut network: Net = source
                    .clone()
                    .into_network(&library, &functions, &settings, Some(&decoding))
                    .unwrap()
                    .into_canonical();
                assert!(network.graph.dangling_indices().is_empty());
                network
                    .execute::<Sequential, SmallestDegree, _, _, _>(&library, &functions)
                    .unwrap();
                let actual: Atom = network.result_scalar().unwrap().into();
                assert_eq!(actual, expected);
                assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
                assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
            }
        }
    }

    #[test]
    fn typed_component_lowering_materializes_open_ports_and_preserves_ranked_zero() {
        crate::test_support::test_initialize();
        let rep = Euclidean {}.new_rep(2).to_symbolic([]);
        let source = SymbolicTensor::infer(function!(
            spenso::tensor_symbol!("typed_component_open_matrix"),
            &rep,
            &rep
        ))
        .unwrap();
        assert_eq!(source.structure.open_positions().len(), 2);
        let library = DummyLibrary::<Tensor>::new();
        let functions = ErroringLibrary::<Symbol>::new();
        let network: Net = source
            .clone()
            .into_network(&library, &functions, &ParseSettings::default(), None)
            .unwrap()
            .into_canonical();
        let slots = network.graph.dangling_indices();
        assert_eq!(slots.len(), 2);
        assert_ne!(slots[0].aind(), slots[1].aind());
        assert_eq!(source.structure.open_positions().len(), 2);

        let zero = source.with_algebra_result(Atom::Zero).unwrap();
        let network: Net = zero
            .into_network(&library, &functions, &ParseSettings::default(), None)
            .unwrap()
            .into_canonical();
        assert_eq!(network.graph.dangling_indices().len(), 2);
        let ExecutionResult::Val(tensor) = network.result_tensor(&library).unwrap() else {
            panic!("ranked zero must retain its tensor interface");
        };
        assert_eq!(tensor.order(), 2);
        assert!(tensor.iter_flat().all(|(_, component)| component.is_zero()));
    }

    #[test]
    fn component_lowering_retains_mixed_logical_axes_for_nonzero_and_zero() {
        use spenso::{
            network::library::symbolic::{ExplicitKey, TensorLibrary},
            structure::representation::Minkowski,
        };

        crate::test_support::test_initialize();
        type Lib = TensorLibrary<ParamTensor<ExplicitKey<AbstractIndex>>, AbstractIndex>;
        type ComponentNet =
            Network<NetworkStore<Tensor, Atom>, ExplicitKey<AbstractIndex>, Symbol, AbstractIndex>;
        let logical_reps = [
            Minkowski {}.new_rep(4).to_lib(),
            Euclidean {}.new_rep(2).to_lib(),
        ];
        let head = spenso::tensor_symbol!("component_lowering_mixed_matrix");
        let signature = ExplicitKey::from_iter(logical_reps, head, None);
        let source = SymbolicTensor::infer(function!(
            head,
            logical_reps[0].to_symbolic([]),
            logical_reps[1].to_symbolic([])
        ))
        .unwrap()
        .multiply(&SymbolicTensor::infer(Atom::num(2)).unwrap())
        .unwrap();
        let mut library = Lib::new();
        library
            .insert_explicit_sparse(
                signature,
                (0..4).flat_map(|row| {
                    (0..2).map(move |column| (vec![row, column], Atom::num(10 * row + column + 1)))
                }),
                Atom::Zero,
            )
            .unwrap();
        let functions = ErroringLibrary::<Symbol>::new();
        let zero = source.with_algebra_result(Atom::Zero).unwrap();
        for (input, is_zero) in [(source, false), (zero, true)] {
            let lowered: Canonicalized<ComponentNet> = input
                .into_network(&library, &functions, &ParseSettings::default(), None)
                .unwrap();
            let (mut network, layout) = lowered.into_parts();
            network
                .execute::<Sequential, SmallestDegree, _, _, _>(&library, &functions)
                .unwrap();
            let ExecutionResult::Val(tensor) = network.result_tensor(&library).unwrap() else {
                panic!("a matrix result must retain both logical axes, including zero");
            };
            let canonical_reps = tensor.external_reps_iter().collect::<Vec<_>>();
            assert_ne!(canonical_reps, logical_reps);
            assert_eq!(layout.canonical_to_logical(&canonical_reps), logical_reps);
            let mut entries = 0;
            for (indices, value) in tensor.iter_expanded() {
                let logical = layout.canonical_to_logical(&indices);
                let expected = if is_zero {
                    Atom::Zero
                } else {
                    Atom::num(2 * (10 * logical[0] + logical[1] + 1))
                };
                assert_eq!(value.to_owned(), expected);
                entries += 1;
            }
            assert!(is_zero || entries == 8);
        }
    }

    #[test]
    fn typed_component_lowering_checks_callbacks_after_index_decoding() {
        test_support::test_initialize();
        let settings = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let raw_index = function!(index_symbol!("typed_component_owned_index"), 7, 2);
        let cooked_index = CookSettings::reversible()
            .cook_function(raw_index.as_fun_view().unwrap())
            .unwrap();
        let rep = Euclidean {}.new_rep(2);
        let raw_slot = rep.to_symbolic([raw_index]);
        let cooked_slot = rep.to_symbolic([cooked_index]);
        assert_ne!(cooked_slot, raw_slot);
        let ordinary = SymbolicTensor::infer(function!(
            tensor_symbol!("typed_component_decoded_vector"),
            &cooked_slot
        ))
        .unwrap();
        let library = DummyLibrary::<Tensor>::new();
        let functions = ErroringLibrary::<Symbol>::new();
        let network: Net = ordinary
            .into_network(
                &library,
                &functions,
                &ParseSettings::default(),
                Some(&settings),
            )
            .unwrap()
            .into_canonical();
        assert_eq!(network.graph.dangling_indices().len(), 1);
        assert_eq!(network.graph.dangling_indices()[0].to_atom(), raw_slot);
        let callback = tensor_symbol!(
            "typed_component_index_callback",
            norm = move |node, output| {
                if let AtomView::Fun(fun) = node
                    && fun.iter().any(|argument| argument == raw_slot.as_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let source = SymbolicTensor::infer(function!(callback, cooked_slot)).unwrap();
        let network: Net = source
            .clone()
            .into_network(&library, &functions, &ParseSettings::default(), None)
            .unwrap()
            .into_canonical();
        assert_eq!(network.graph.dangling_indices().len(), 1);
        let decoded: Result<Canonicalized<Net>, _> = source.clone().into_network(
            &library,
            &functions,
            &ParseSettings::default(),
            Some(&settings),
        );
        assert!(
            decoded
                .unwrap_err()
                .to_string()
                .contains("admitted tensor interface")
        );
        assert_eq!(source.rank(), 1);
    }

    #[test]
    fn component_lowering_releases_symbolic_graph_before_decoding() {
        test_support::test_initialize();
        let settings = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let raw_index = function!(index_symbol!("component_lifetime_index"), 11, 3);
        let cooked_index = CookSettings::reversible()
            .cook_function(raw_index.as_fun_view().unwrap())
            .unwrap();
        let rep = Euclidean {}.new_rep(2);
        let raw_slot = rep.to_symbolic([raw_index]);
        let cooked_slot = rep.to_symbolic([cooked_index]);
        assert_ne!(cooked_slot, raw_slot);
        let graph = Arc::new(Mutex::new(None::<Weak<PlanningNetwork<Atom>>>));
        let observed = Arc::clone(&graph);
        let decoded = Arc::new(AtomicBool::new(false));
        let did_decode = Arc::clone(&decoded);
        let callback = tensor_symbol!(
            "component_lifetime_tensor",
            norm = move |node, _output| {
                if let AtomView::Fun(function) = node
                    && function
                        .iter()
                        .any(|argument| argument == raw_slot.as_view())
                {
                    let graph = observed.lock().unwrap();
                    assert!(graph.as_ref().unwrap().upgrade().is_none());
                    did_decode.store(true, Ordering::Relaxed);
                }
            }
        );
        let source = SymbolicTensor::infer(function!(callback, cooked_slot)).unwrap();
        let planned = source.shallow_graph().unwrap();
        *graph.lock().unwrap() = Some(Arc::downgrade(&planned));
        drop(planned);
        let network: Net = source
            .into_network(
                &DummyLibrary::<Tensor>::new(),
                &ErroringLibrary::<Symbol>::new(),
                &ParseSettings::default(),
                Some(&settings),
            )
            .unwrap()
            .into_canonical();
        assert_eq!(network.graph.dangling_indices().len(), 1);
        assert!(decoded.load(Ordering::Relaxed));
    }
}
