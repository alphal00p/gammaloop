//! Checked lowering from a symbolic carrier to finite component execution.

use std::fmt::{Debug, Display};

use spenso::{
    network::{
        Network, TensorNetworkError,
        library::FunctionLibrary,
        parsing::{ParseSettings, StructureFromAtom, TensorFromExpression, TensorLibraryFor},
        store::NetworkStore,
    },
    structure::{
        Canonicalized, HasStructure, OrderedStructure, ScalarStructure, TensorStructure,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
        representation::{LibraryRep, LibrarySlot},
        slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind},
    },
    tensors::data::SparseTensor,
};
use symbolica::atom::{Atom, Symbol};

use super::SymbolicTensor;
use crate::CookSettings;

impl SymbolicTensor<PartialStructure> {
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
        self.ensure_validated()
            .map_err(|error| invalid(error.to_string()))?;
        let source = if self.structure.open_positions().is_empty() {
            self
        } else {
            let source = self
                .materialized()
                .map_err(|error| invalid(error.to_string()))?;
            drop(self);
            source
        };
        // Retain the admitted boundary and payload, but release symbolic planning
        // graphs before allocating the independent component-execution network.
        let (expression, structure) = source.into_parts();
        let decode = |expression: &Atom| match index_decoding {
            Some(settings) => settings.uncook(expression.as_view()),
            None => expression.clone(),
        };
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

        let expression = match index_decoding {
            Some(_) => {
                let decoded = decode(&expression);
                drop(expression);
                decoded
            }
            None => expression,
        };
        let network =
            Network::<NetworkStore<T, Atom>, K, Symbol, Aind>::try_from_view_with_function_library(
                expression.as_view(),
                library,
                functions,
                settings,
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

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::{
        iterators::IteratableTensor,
        network::{
            ExecutionResult, Sequential, SmallestDegree,
            library::{DummyKey, DummyLibrary, panicing::ErroringLibrary},
            parsing::{ShadowedStructure, ShorthandParsing},
        },
        structure::{
            abstract_index::AbstractIndex,
            representation::{Euclidean, RepName},
        },
        tensors::parametric::ParamTensor,
    };
    use symbolica::{
        atom::{AtomCore, AtomView},
        function,
    };

    type Structure = ShadowedStructure<AbstractIndex>;
    type Tensor = ParamTensor<Structure>;
    type Net = Network<NetworkStore<Tensor, Atom>, DummyKey, Symbol, AbstractIndex>;

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
        crate::test_support::test_initialize();
        let settings =
            crate::CookSettings::indices().with_mode(crate::CookMode::ReversibleEncoding);
        let raw_slot = function!(
            Euclidean {}
                .to_symbolic(std::iter::empty::<Atom>())
                .as_fun_view()
                .unwrap()
                .get_symbol(),
            2,
            function!(spenso::index_symbol!("typed_component_owned_index"), 7, 2)
        );
        let cooked_slot = settings.try_cook_indices(raw_slot.as_view()).unwrap();
        assert_ne!(cooked_slot, raw_slot);
        let ordinary = SymbolicTensor::infer(function!(
            spenso::tensor_symbol!("typed_component_decoded_vector"),
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
        let callback = spenso::tensor_symbol!(
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
        use std::sync::{
            Arc, Mutex, Weak,
            atomic::{AtomicBool, Ordering},
        };

        crate::test_support::test_initialize();
        let settings =
            crate::CookSettings::indices().with_mode(crate::CookMode::ReversibleEncoding);
        let raw_slot = Euclidean {}.new_rep(2).to_symbolic([function!(
            spenso::index_symbol!("component_lifetime_index"),
            11,
            3
        )]);
        let cooked_slot = settings.try_cook_indices(raw_slot.as_view()).unwrap();
        let graph = Arc::new(Mutex::new(
            None::<Weak<crate::tensor::SymbolicNet<AbstractIndex>>>,
        ));
        let observed = Arc::clone(&graph);
        let decoded = Arc::new(AtomicBool::new(false));
        let did_decode = Arc::clone(&decoded);
        let callback = spenso::tensor_symbol!(
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
