//! Finite component materialization of the existing literal alias registry.

use super::*;
use crate::{CookSettings, tensor::composition::matching_interface_position};
use spenso::{
    algebra::ScalarMul,
    network::{
        ExecutionResult, FastTensorSum, Network, Ref, TensorNetworkError,
        graph::NMul,
        library::{DummyLibrary, FunctionLibrary, function_lib::Wrap},
        parsing::{
            ParseSettings, StructureFromAtom, StructureInferenceMode, TensorFromExpression,
            TensorLibraryFor,
        },
        store::{NetworkStore, TensorScalarStoreMapping},
    },
    structure::{
        ApplyPendingIndexPermutation, HasStructure, OrderedStructure, ScalarStructure,
        ScalarTensor, TensorStructure,
        representation::{LibraryRep, LibrarySlot},
        slot::{AbsInd, DummyAind, ParseableAind, Slot},
    },
};
use std::{
    fmt::{Debug, Display},
    ops::AddAssign,
};
use symbolica::atom::Symbol;

impl SymbolicTensor<AliasInterfaces, AliasedAtom> {
    /// Materialize literal tensor definitions through the existing component network.
    ///
    /// Tensor values are cached by the complete literal handle, before a component
    /// library can erase index labels. Scalar definitions remain registered; the
    /// caller merges this registry with its network result before scoping handles.
    /// Rejected scalar aliases are substituted at their literal uses, in dependency
    /// order, so selectors remain visible without resolving the tensor DAG.
    /// Disconnected definitions are validated by the registry but are not executed.
    ///
    /// `execute` and `tensor_map` are the caller's existing finite-component policy.
    /// Cooking is reversed only at this boundary; the shared carrier is unchanged.
    #[allow(
        clippy::type_complexity,
        clippy::too_many_arguments,
        clippy::result_large_err
    )]
    pub fn to_network<T, S, K, Aind, Lib, FunLib>(
        &self,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        cooking: Option<&CookSettings>,
        mut execute: impl FnMut(
            &mut Network<NetworkStore<T, Atom>, K, Symbol, Aind>,
        ) -> std::result::Result<(), TensorNetworkError<K, Symbol>>,
        mut tensor_map: impl FnMut(T) -> std::result::Result<T, TensorNetworkError<K, Symbol>>,
        mut retain_scalar_alias: impl FnMut(&Atom) -> bool,
    ) -> std::result::Result<
        (Network<NetworkStore<T, Atom>, K, Symbol, Aind>, AliasedAtom),
        TensorNetworkError<K, Symbol>,
    >
    where
        Aind: AbsInd + DummyAind + ParseableAind,
        K: Clone + Debug + Display,
        T: Clone
            + HasStructure<Structure = S>
            + TensorStructure<Slot = LibrarySlot<Aind>, Indexed = T>
            + ScalarTensor<Scalar = Atom>
            + ApplyPendingIndexPermutation<Output = T>
            + Ref
            + FastTensorSum
            + ScalarMul<Atom, Output = T>
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
        S: Clone + ScalarStructure + StructureFromAtom + TensorStructure<Slot = LibrarySlot<Aind>>,
        for<'a> T: TensorFromExpression<'a, S, Atom, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        Lib::LibraryTensor: TensorStructure<Indexed = T>,
        FunLib: FunctionLibrary<T, Atom, Key = Symbol>,
    {
        let error = |message: &str| TensorNetworkError::Other(eyre::eyre!("{message}"));
        let definitions = self.aliases().map_err(|e| error(&e.to_string()))?;
        let uncook =
            |value: &Atom| cooking.map_or_else(|| value.clone(), |c| c.uncook(value.as_view()));
        let zero = |handle: &Atom| {
            let structure = S::structure_from_atom(handle.as_view(), StructureInferenceMode::Fast)?;
            let tensor = <T as TensorFromExpression<
                '_,
                S,
                Atom,
                K,
                Symbol,
                Aind,
                Lib,
                FunLib,
            >>::tensor_from_leaf(handle.as_view().into(), structure)?;
            tensor
                .scalar_mul(&Atom::Zero)
                .ok_or(TensorNetworkError::FailedScalarMul)
        };
        let mut values: HashMap<Atom, T> = HashMap::new();
        let mut visible: HashMap<Atom, Atom> = HashMap::new();
        let mut registry = AliasedAtom::from(Atom::one());
        let substitute = |value: &Atom, visible: &HashMap<Atom, Atom>| {
            value.replace_map(|node, _, output| {
                if let Some(body) = visible.get(node.get_data()) {
                    output.set_from_view(&body.as_view());
                }
            })
        };
        let build =
            |value: &Atom,
             values: &HashMap<Atom, T>,
             execute: &mut dyn FnMut(
                &mut Network<NetworkStore<T, Atom>, K, Symbol, Aind>,
            )
                -> std::result::Result<(), TensorNetworkError<K, Symbol>>,
             tensor_map: &mut dyn FnMut(
                T,
            ) -> std::result::Result<
                T,
                TensorNetworkError<K, Symbol>,
            >| {
                // A dummy library retains every ordinary tensor's exact source literal.
                // Component-library lookup is deferred until map_occurrences below.
                let graph = Network::<
                    NetworkStore<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, Atom>,
                    K,
                    Symbol,
                    Aind,
                >::try_from_view_with_function_library(
                    value.as_view(),
                    &DummyLibrary::new(),
                    &Wrap,
                    settings,
                )?;
                graph.map_occurrences(
                |scalar| Ok(scalar.clone()),
                |source, ports, _logical| {
                    let mut component = if let Some(value) = values.get(source.expression()) {
                        value.clone()
                    } else {
                        let mut leaf = Network::<NetworkStore<T, Atom>, K, Symbol, Aind>
                            ::try_from_view_with_function_library(
                                source.expression().as_view(), library, function_library, settings,
                            )?;
                        leaf = leaf.map_result(Ok, &mut *tensor_map)?;
                        execute(&mut leaf)?;
                        match leaf.result_tensor(library)? {
                            ExecutionResult::Val(value) => value.into_owned(),
                            ExecutionResult::Zero => T::new_scalar(Atom::Zero),
                            ExecutionResult::One => T::new_scalar(Atom::one()),
                        }
                    };
                    let original = source.structure().external_structure();
                    let indices =
                        component
                            .external_structure_iter()
                            .map(|slot| {
                                let position =
                                    original.iter().position(|old| *old == slot).ok_or_else(
                                        || error("component alias changed its encoded interface"),
                                    )?;
                                ports
                                    .iter()
                                    .find(|(storage, _, _)| *storage == position)
                                    .map(|(_, slot, _)| slot.aind())
                                    .ok_or_else(|| error("component alias occurrence lost a port"))
                            })
                            .collect::<std::result::Result<Vec<_>, _>>()?;
                    component = component.reindex_storage(&indices)?.apply();
                    let mut network = Network::from_tensor(component);
                    let mut bound = false;
                    for (_, slot, value) in ports {
                        if let Some(value) = value {
                            // The ordinary component parser owns compact-vector
                            // recognition. Bind its one surviving axis using the
                            // graph's already admitted port representation.
                            let mut vector = Network::<NetworkStore<T, Atom>, K, Symbol, Aind>
                                ::try_from_view_with_function_library(
                                    value.as_view(), library, function_library, settings,
                                )?;
                            vector = vector.map_result(Ok, &mut *tensor_map)?;
                            execute(&mut vector)?;
                            let ExecutionResult::Val(vector) = vector.result_tensor(library)?
                            else {
                                return Err(error("bound component vector has no tensor value"));
                            };
                            if vector.order() != 1
                                || vector.external_structure()[0].rep() != slot.rep()
                            {
                                return Err(error(
                                    "bound component vector changed its admitted representation",
                                ));
                            }
                            let vector =
                                vector.into_owned().reindex_storage(&[slot.aind()])?.apply();
                            network = network.n_mul([Network::from_tensor(vector)]);
                            bound = true;
                        }
                    }
                    if bound {
                        execute(&mut network)?;
                    }
                    match network.result_tensor(library)? {
                        ExecutionResult::Val(value) => Ok(value.into_owned()),
                        ExecutionResult::Zero => Ok(T::new_scalar(Atom::Zero)),
                        ExecutionResult::One => Ok(T::new_scalar(Atom::one())),
                    }
                },
            )
            };
        let reachable = self.reachable_definitions();
        for (handle, body) in definitions {
            if !reachable.contains(handle.expression()) {
                continue;
            }
            let scalar = handle.is_scalar();
            let handle_atom = uncook(handle.expression());
            let body_atom = substitute(&uncook(body.expression()), &visible);
            if body_atom.is_zero() && !scalar {
                values.insert(handle_atom.clone(), zero(&handle_atom)?);
                continue;
            }
            let mut network = build(&body_atom, &values, &mut execute, &mut tensor_map)?;
            execute(&mut network)?;
            let component = match network.result_tensor(library)? {
                ExecutionResult::Val(value) => value.into_owned(),
                ExecutionResult::Zero => T::new_scalar(Atom::Zero),
                ExecutionResult::One => T::new_scalar(Atom::one()),
            };
            if component.structure().is_scalar() != scalar {
                return Err(error("component alias changed its declared rank"));
            }
            if scalar {
                let body = component
                    .scalar()
                    .ok_or_else(|| error("scalar component has no scalar payload"))?;
                if retain_scalar_alias(&body) {
                    registry.register_alias(handle_atom, body);
                } else {
                    visible.insert(handle_atom, body);
                }
            } else {
                // Explicit labels keep their identity even if handle and body
                // declare different logical presentations. AUTO axes instead use
                // the composition owner's positional policy, never owner equality.
                let logical = body.structure().logical_slots();
                let mut claimed = vec![false; logical.len()];
                let encoded =
                    InterfaceInference::replacement_interface(body.expression().as_view())
                        .map_err(|e| error(&e.to_string()))?;
                let mut body_ports = vec![None; logical.len()];
                for slot in encoded.logical_slots() {
                    let PartialIndex::Explicit(index) = slot.aind else {
                        return Err(error("component alias body has an unresolved encoded port"));
                    };
                    let atom = slot.rep().slot::<AbstractIndex, _>(index).to_atom();
                    let position = matching_interface_position(atom.as_view(), &logical, &claimed)
                        .ok_or_else(|| {
                            error("component alias body changed its logical interface")
                        })?;
                    claimed[position] = true;
                    body_ports[position] = Some(
                        Slot::<LibraryRep, Aind>::try_from(uncook(&atom).as_view())
                            .map_err(spenso::structure::StructureError::from)?,
                    );
                }
                let AtomView::Fun(handle_function) = handle.expression().as_view() else {
                    return Err(error(
                        "tensor component alias has no literal function handle",
                    ));
                };
                claimed.fill(false);
                let mut replacements = HashMap::new();
                for argument in handle_function.iter() {
                    if let Some(position) =
                        matching_interface_position(argument, &logical, &claimed)
                    {
                        claimed[position] = true;
                        let source = body_ports[position]
                            .ok_or_else(|| error("component alias body lost an encoded port"))?;
                        let target = Slot::<LibraryRep, Aind>::try_from(
                            uncook(&argument.to_owned()).as_view(),
                        )
                        .map_err(spenso::structure::StructureError::from)?;
                        replacements.insert(source, target.aind());
                    }
                }
                let indices = component
                    .external_structure_iter()
                    .map(|slot| {
                        replacements
                            .get(&slot)
                            .copied()
                            .ok_or_else(|| error("component alias changed its encoded interface"))
                    })
                    .collect::<std::result::Result<Vec<_>, _>>()?;
                if indices.len() != logical.len() {
                    return Err(error("component alias changed its declared interface"));
                }
                values.insert(handle_atom, component.reindex_storage(&indices)?.apply());
            }
        }
        let root = substitute(&uncook(self.expression.get_root()), &visible);
        let network = if root.is_zero() && !self.structure.root.canonical().is_scalar() {
            let handle = self
                .root()
                .alias_handle()
                .map_err(|e| error(&e.to_string()))?;
            Network::from_tensor(zero(&uncook(handle.expression()))?)
        } else {
            build(&root, &values, &mut execute, &mut tensor_map)?
        };
        Ok((network, registry))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::{
        iterators::IteratableTensor,
        network::{
            Sequential, SmallestDegree,
            library::{
                panicing::ErroringLibrary,
                symbolic::{ExplicitKey, TensorLibrary},
            },
            parsing::ShadowedStructure,
        },
        structure::{
            representation::{Euclidean, RepName},
            slot::IsAbstractSlot,
        },
        tensors::parametric::ParamTensor,
    };
    use symbolica::{function, symbol};

    type Tensor = ParamTensor<ShadowedStructure<AbstractIndex>>;
    type Lib = TensorLibrary<ParamTensor<ExplicitKey<AbstractIndex>>, AbstractIndex>;
    type Net = Network<NetworkStore<Tensor, Atom>, ExplicitKey<AbstractIndex>, Symbol>;

    fn execute(
        net: &mut Net,
        lib: &Lib,
    ) -> std::result::Result<(), TensorNetworkError<ExplicitKey<AbstractIndex>, Symbol>> {
        net.execute::<Sequential, SmallestDegree, _, _, _>(lib, &ErroringLibrary::new())?;
        Ok(())
    }

    fn scalar(expression: Atom) -> SymbolicTensor<PartialStructure> {
        SymbolicTensor::checked_parts(expression, PartialStructure::from_logical_slots([])).unwrap()
    }

    fn result(net: &Net) -> Atom {
        match net.result_scalar().unwrap() {
            ExecutionResult::Val(value) => value.into_owned(),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::One => Atom::one(),
        }
    }

    fn fixture() -> (Lib, Symbol, Symbol) {
        crate::test_support::test_initialize();
        let p = spenso::vector_symbol!("alias_components::p");
        let q = spenso::vector_symbol!("alias_components::q");
        let mut lib = Lib::new();
        for (symbol, components) in [(p, [1, 2]), (q, [3, 4])] {
            lib.insert_explicit_sparse(
                ExplicitKey::from_iter([LibraryRep::from(Euclidean {}).new_rep(2)], symbol, None),
                components
                    .into_iter()
                    .enumerate()
                    .map(|(i, v)| (vec![i], Atom::num(v))),
                Atom::Zero,
            )
            .unwrap();
        }
        (lib, p, q)
    }

    fn vector(symbol: Symbol, index: usize) -> SymbolicTensor<PartialStructure> {
        let slot = LibraryRep::from(Euclidean {})
            .new_rep(2)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
            .to_atom();
        SymbolicTensor::infer(function!(symbol, slot)).unwrap()
    }

    #[test]
    fn component_aliases_use_exact_literal_labels_before_library_keys() {
        let (lib, p, q) = fixture();
        let pa = vector(p, 97001);
        let qb = vector(q, 97003);
        let first = pa.alias_handle().unwrap();
        let second = first
            .reindex_interface_ports(&HashMap::from([(0, AbstractIndex::Normal(97003))]))
            .unwrap();
        let root = SymbolicTensor::infer(
            first.expression() * second.expression() * pa.expression() * qb.expression(),
        )
        .unwrap();
        let value = root.with_aliases([(first, pa), (second, qb)]).unwrap();
        let before = value.clone();
        let (mut net, registry) = value
            .to_network(
                &lib,
                &ErroringLibrary::new(),
                &ParseSettings::default(),
                None,
                |net| execute(net, &lib),
                |_| panic!("a library tensor must not be remapped as a parsed local tensor"),
                |_| true,
            )
            .unwrap();
        execute(&mut net, &lib).unwrap();
        assert_eq!(result(&net), Atom::num(125));
        assert!(registry.get_aliases().is_empty());
        assert_eq!(value.root(), before.root());
        assert_eq!(value.aliases().unwrap(), before.aliases().unwrap());
    }

    #[test]
    fn component_scalar_aliases_evaluate_tensor_bodies_and_expose_transitive_guards() {
        let (lib, p, _) = fixture();
        let pa = vector(p, 97101);
        let norm = SymbolicTensor::infer(pa.expression().pow(2)).unwrap();
        let handle = norm.alias_handle().unwrap();
        let first = scalar(Atom::var(symbol!("alias_components::guarded")));
        let second = scalar(Atom::var(symbol!("alias_components::transitive")));
        let guard = symbol!("alias_components::selector");
        let x = Atom::var(symbol!("alias_components::x"));
        let guarded = function!(guard, &x) * handle.expression();
        let body = (first.expression() + 1).pow(2);
        let value = second
            .clone()
            .with_aliases([
                (second, scalar(body)),
                (first, scalar(guarded)),
                (handle.clone(), norm),
            ])
            .unwrap();
        let (mut net, registry) = value
            .to_network(
                &lib,
                &ErroringLibrary::new(),
                &ParseSettings::default(),
                None,
                |net| execute(net, &lib),
                Ok,
                |body| !body.contains_symbol(guard),
            )
            .unwrap();
        execute(&mut net, &lib).unwrap();
        assert_eq!(registry.get_aliases().len(), 1);
        assert_eq!(registry.get_aliases()[handle.expression()], Atom::num(5));
        assert_eq!(
            result(&net),
            (function!(guard, x) * handle.expression() + 1).pow(2)
        );
        assert!(result(&net).contains_symbol(guard));
    }

    #[test]
    fn component_aliases_keep_open_typed_zero_and_skip_disconnected_definitions() {
        let (lib, p, q) = fixture();
        let body = vector(p, 97201);
        let handle = body.alias_handle().unwrap();
        let unused = vector(q, 97203);
        let zero = body.with_checked_expression(Atom::Zero).unwrap();
        let value = handle
            .clone()
            .with_aliases([
                (handle.clone(), zero),
                (unused.alias_handle().unwrap(), unused),
            ])
            .unwrap();
        let mut executions = 0;
        let (mut net, registry) = value
            .to_network(
                &lib,
                &ErroringLibrary::new(),
                &ParseSettings::default(),
                None,
                |net| {
                    executions += 1;
                    execute(net, &lib)
                },
                Ok,
                |_| true,
            )
            .unwrap();
        assert_eq!(executions, 0);
        execute(&mut net, &lib).unwrap();
        let ExecutionResult::Val(tensor) = net.result_tensor(&lib).unwrap() else {
            panic!("open zero must retain its component interface")
        };
        assert_eq!(tensor.order(), 1);
        assert_eq!(
            tensor.external_structure(),
            vec![
                LibraryRep::from(Euclidean {})
                    .new_rep(2)
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(97201))
            ]
        );
        assert!(registry.get_aliases().is_empty());
    }

    #[test]
    fn component_aliases_align_auto_owners_without_swapping_explicit_labels() {
        let (mut lib, p, q) = fixture();
        let matrix = spenso::tensor_symbol!("alias_components::matrix");
        let rep = LibraryRep::from(Euclidean {}).new_rep(2);
        lib.insert_explicit_sparse(
            ExplicitKey::from_iter([rep, rep], matrix, None),
            [1, 2, 3, 4]
                .into_iter()
                .enumerate()
                .map(|(i, v)| (vec![i / 2, i % 2], Atom::num(v))),
            Atom::Zero,
        )
        .unwrap();
        let owner = AbstractIndex::fresh_open_owner();
        let auto = SymbolicTensor::infer(
            FunctionBuilder::new(matrix)
                .add_args((0..2).map(|axis| {
                    rep.slot::<AbstractIndex, _>(AbstractIndex::Open { owner, axis })
                        .to_atom()
                }))
                .finish(),
        )
        .unwrap();
        let auto_handle = auto.alias_handle().unwrap();
        assert_ne!(auto.expression(), auto_handle.expression());
        let value = auto_handle
            .clone()
            .with_aliases([(auto_handle, auto)])
            .unwrap();
        let (mut net, _) = value
            .to_network(
                &lib,
                &ErroringLibrary::new(),
                &ParseSettings::default(),
                None,
                |net| execute(net, &lib),
                Ok,
                |_| true,
            )
            .unwrap();
        execute(&mut net, &lib).unwrap();
        let ExecutionResult::Val(tensor) = net.result_tensor(&lib).unwrap() else {
            panic!("matrix missing")
        };
        let mut entries = tensor
            .iter_flat()
            .map(|(index, value)| (index, value.to_owned()))
            .collect::<Vec<_>>();
        entries.sort_by_key(|(index, _)| *index);
        assert_eq!(
            entries
                .into_iter()
                .map(|(_, value)| value)
                .collect::<Vec<_>>(),
            [1, 2, 3, 4].map(Atom::num)
        );

        let compact = function!(p, rep.to_symbolic([]));
        let body = SymbolicTensor::infer(function!(
            matrix,
            compact,
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(97303))
                .to_atom()
        ))
        .unwrap();
        let handle = body.alias_handle().unwrap();
        let root =
            SymbolicTensor::infer(handle.expression() * vector(q, 97303).expression()).unwrap();
        let compact_value = root.with_aliases([(handle, body)]).unwrap();
        let (mut net, _) = compact_value
            .to_network(
                &lib,
                &ErroringLibrary::new(),
                &ParseSettings::default(),
                None,
                |net| execute(net, &lib),
                Ok,
                |_| true,
            )
            .unwrap();
        execute(&mut net, &lib).unwrap();
        assert_eq!(result(&net), Atom::num(61));

        let indices = [97301, 97303].map(AbstractIndex::Normal);
        let body = SymbolicTensor::infer(
            FunctionBuilder::new(matrix)
                .add_args(indices.map(|i| rep.slot::<AbstractIndex, _>(i).to_atom()))
                .finish(),
        )
        .unwrap();
        let mut handle = body.alias_handle().unwrap();
        handle.structure = PartialStructure::from_logical_slots(
            handle.structure.logical_slots().into_iter().rev(),
        );
        let root = SymbolicTensor::infer(
            handle.expression() * vector(p, 97301).expression() * vector(q, 97303).expression(),
        )
        .unwrap();
        let value = root.with_aliases([(handle, body)]).unwrap();
        let (mut net, _) = value
            .to_network(
                &lib,
                &ErroringLibrary::new(),
                &ParseSettings::default(),
                None,
                |net| execute(net, &lib),
                Ok,
                |_| true,
            )
            .unwrap();
        execute(&mut net, &lib).unwrap();
        assert_eq!(result(&net), Atom::num(61));
    }

    #[test]
    fn component_aliases_preserve_open_zero_from_finite_cancellation() {
        let (mut lib, p, _) = fixture();
        let same = spenso::vector_symbol!("alias_components::same_as_p");
        let rep = LibraryRep::from(Euclidean {}).new_rep(2);
        lib.insert_explicit_sparse(
            ExplicitKey::from_iter([rep], same, None),
            [(vec![0], Atom::num(1)), (vec![1], Atom::num(2))],
            Atom::Zero,
        )
        .unwrap();
        let body = vector(p, 97401).try_sub(&vector(same, 97401)).unwrap();
        assert!(!body.expression().is_zero());
        let handle = body.alias_handle().unwrap();
        let value = handle.clone().with_aliases([(handle, body)]).unwrap();
        let (mut net, _) = value
            .to_network(
                &lib,
                &ErroringLibrary::new(),
                &ParseSettings::default(),
                None,
                |net| execute(net, &lib),
                Ok,
                |_| true,
            )
            .unwrap();
        execute(&mut net, &lib).unwrap();
        let ExecutionResult::Val(tensor) = net.result_tensor(&lib).unwrap() else {
            panic!("computed open zero lost its interface")
        };
        assert_eq!(
            tensor.external_structure(),
            vec![rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(97401))]
        );
        assert!(tensor.iter_flat().all(|(_, value)| value.is_zero()));
    }
}
