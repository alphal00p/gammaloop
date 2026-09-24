use crate::network::library::DummyKey;

use super::TensorNetworkError;

#[test]
fn display() {
    let a = TensorNetworkError::<i8, DummyKey>::Infallible;

    println!("{a}")
}

#[test]
fn display_labels_match_network_errors() {
    let invalid_dot = TensorNetworkError::<i8, DummyKey>::InvalidDotFunction("dot(p)".to_owned());
    let invalid_power =
        TensorNetworkError::<i8, DummyKey>::NonSelfDualTensorPower("T^2".to_owned());

    assert_eq!(invalid_dot.to_string(), "Invalid dot function: dot(p)");
    assert_eq!(invalid_power.to_string(), "Non self-dual tensor power: T^2");
}

#[cfg(feature = "shadowing")]
#[test]
fn scalar_alias_refs_resolve_to_original_atom() {
    use symbolica::{
        atom::{Atom, AtomCore},
        function, parse, symbol,
    };

    use super::{
        Network,
        store::{NetworkStore, TensorScalarStore},
        tags::scalar_store_alias,
    };

    let original = Atom::var(symbol!("x"));
    let mut net: Network<NetworkStore<(), Atom>, i8, i8> = Network::from_scalar(original.clone());
    let aliases = net.alias_scalar_refs(|_, _| true);

    assert_eq!(
        net.resolve_scalar_aliases(&aliases, scalar_store_alias(0)),
        original
    );

    // Network results may mix opaque references with newly formed scalar bodies.
    // Definitions can share a value or refer to earlier aliases without requiring
    // a second common-subexpression extraction over that result.
    let original = parse!("x+y");
    let mut net: Network<NetworkStore<(), Atom>, i8, i8> = Network::from_scalar(original.clone());
    let nested = net.store.add_scalar(scalar_store_alias(0) + parse!("z"));
    let repeated = net.store.add_scalar(original.clone());
    net.store.add_scalar(parse!("unused"));
    let aliases = net.alias_scalar_refs(|_, _| true);

    let root = scalar_store_alias(nested).pow(2)
        + function!(symbol!("f"), scalar_store_alias(0))
        + scalar_store_alias(repeated).pow(-1)
        + original.pow(3)
        + scalar_store_alias(99);
    let expected = parse!("(x+y+z)^2+f(x+y)+1/(x+y)+(x+y)^3") + scalar_store_alias(99);
    let aliased = net.aliased_atom(&aliases, root.clone());
    assert_eq!(aliased.get_root(), &root);
    assert_eq!(aliased.clone().into_inner(), expected);
    assert_eq!(net.resolve_scalar_aliases(&aliases, root.clone()), expected);

    // The former recompression path could rewrite the constructed root.
    // Registration must preserve that factorization and resolve every alias;
    // alias-map layout is not part of the mathematical result.
    for root in [Atom::Zero, Atom::num(1), scalar_store_alias(99)] {
        assert_eq!(net.resolve_scalar_aliases(&aliases, root.clone()), root);
    }
    let aliases = net.alias_scalar_refs(|_, _| false);
    assert!(aliases.is_empty());
    assert_eq!(
        net.resolve_scalar_aliases(&aliases, expected.clone()),
        expected
    );
}

#[cfg(feature = "shadowing")]
#[test]
fn scalar_alias_resolution_preserves_nested_and_unregistered_handles() {
    use symbolica::{
        atom::{Atom, AtomCore},
        function, parse, symbol,
    };

    use super::{
        Network,
        store::{NetworkStore, TensorScalarStore},
        tags::{SPENSO_TAG, scalar_store_alias},
    };

    let mut net: Network<NetworkStore<(), Atom>, i8, i8> = Network::from_scalar(Atom::num(1));
    net.store.add_scalar(scalar_store_alias(3));
    net.store.add_scalar(parse!("unregistered"));
    net.store.add_scalar(parse!("(x+y)^7*(z+w)^5"));
    net.store.add_scalar(scalar_store_alias(4));
    let aliases = net.alias_scalar_refs(|index, _| index != 2);

    // Forward definitions and handles produced by normalization require the
    // same fixed point as the generic alias map. Spectator powers stay intact.
    let nested = function!(SPENSO_TAG.scalar, scalar_store_alias(0));
    let expected = parse!("(x+y)^7*(z+w)^5");
    assert_eq!(
        net.resolve_scalar_aliases(&aliases, nested.clone()),
        expected
    );
    for root in [
        nested,
        scalar_store_alias(1).pow(3) * parse!("(a+b)^4"),
        function!(symbol!("f"), scalar_store_alias(1)),
        scalar_store_alias(2),
        scalar_store_alias(4),
        scalar_store_alias(99),
        Atom::var(SPENSO_TAG.scalar),
        function!(SPENSO_TAG.scalar, -1),
        function!(SPENSO_TAG.scalar, parse!("1/2")),
        function!(SPENSO_TAG.scalar, 0, 1),
        function!(SPENSO_TAG.scalar, scalar_store_alias(0), parse!("x")),
        function!(SPENSO_TAG.scalar, parse!("x")),
        expected,
    ] {
        let expected = net.aliased_atom(&aliases, root.clone()).into_inner();
        assert_eq!(net.resolve_scalar_aliases(&aliases, root), expected);
    }
}

#[cfg(feature = "shadowing")]
#[test]
fn auto_serializes_unlicensed_symbolic_fast_tensor_sum() {
    use std::collections::HashMap;

    use symbolica::{atom::Atom, parse};

    use crate::{
        network::FastTensorSum,
        structure::{
            OrderedStructure,
            concrete_index::FlatIndex,
            representation::{Euclidean, RepName},
        },
        symbolic_parallelism::{
            SymbolicParallelism, scoped_symbolica_rayon_setting_for_test, symbolica_rayon_enabled,
        },
        tensors::{
            data::{DataTensor, SparseTensor},
            parametric::ParamTensor,
        },
    };

    // The repository's Symbolica build may itself be licensed. Injecting the
    // unlicensed result makes this failure mode deterministic while exercising
    // the exact same Auto resolution and cached global used in production.
    let _guard = scoped_symbolica_rayon_setting_for_test(SymbolicParallelism::Auto, || false);
    assert!(!symbolica_rayon_enabled());

    let structure: OrderedStructure<Euclidean> =
        OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).into_canonical();
    let tensor = |first: Atom, second: Atom| {
        ParamTensor::composite(DataTensor::Sparse(SparseTensor {
            elements: HashMap::from([(FlatIndex::from(0), first), (FlatIndex::from(1), second)]),
            zero: Atom::Zero,
            structure: structure.clone(),
        }))
    };
    let left = tensor(parse!("x"), parse!("y"));
    let right = tensor(parse!("a"), parse!("b"));

    // Each output entry contains two atoms, so this executes the guarded
    // Atom::add_many path that previously always ran on a Rayon worker.
    let result = <ParamTensor<_> as FastTensorSum>::fast_tensor_sum(&[&left, &right], None)
        .expect("two compatible sparse tensors should use the fast symbolic sum");
    let DataTensor::Sparse(result) = result.tensor else {
        panic!("the fast sparse sum should remain sparse");
    };
    assert_eq!(result.elements[&FlatIndex::from(0)], parse!("a+x"));
    assert_eq!(result.elements[&FlatIndex::from(1)], parse!("b+y"));
}

#[cfg(feature = "shadowing")]
#[test]
fn fast_tensor_sum_parallel_candidate_rejects_light_or_unbalanced_work() {
    use super::{AtomSumShapeStats, FastTensorSumWorkload};

    let candidate =
        |logical_entries, entries, total_terms, total_bytes, max_group_terms, max_group_bytes| {
            FastTensorSumWorkload {
                shape: AtomSumShapeStats {
                    logical_entries,
                    entries,
                    total_terms,
                    total_bytes,
                    ..Default::default()
                },
                max_group_terms,
                max_group_bytes,
            }
        };

    // Many one-term atoms still lose to Rayon because each individual merge is
    // too cheap, even though their aggregate size is substantial.
    let trivial = candidate(256, 1_024, 1_024, 64 * 1024, 4, 256);
    assert!(!trivial.meets_parallel_shape_floor());

    // This balanced shape represents the first conservatively profitable
    // benchmark bucket: each merge combines expression-heavy atoms.
    let balanced = candidate(4, 16, 512, 32 * 1024, 128, 8 * 1024);
    assert!(balanced.meets_parallel_shape_floor());
    assert!(balanced.parallel_is_profitable_for(8));

    let too_few_groups = candidate(2, 8, 512, 32 * 1024, 256, 16 * 1024);
    assert!(!too_few_groups.meets_parallel_shape_floor());

    // Disjoint sparse supports only clone one Atom per output group; they do
    // not exercise the add-many work calibrated by the current benchmark.
    let disjoint = candidate(8, 8, 768, 32 * 1024, 96, 4 * 1024);
    assert!(!disjoint.meets_parallel_shape_floor());

    assert!(!balanced.parallel_is_profitable_for(1));
    assert!(!balanced.parallel_is_profitable_for(16));

    let unbalanced = candidate(4, 16, 512, 32 * 1024, 400, 20 * 1024);
    assert!(!unbalanced.meets_parallel_shape_floor());
}

#[test]
fn executed_scaled_tensors_add_distinct_tensors() {
    use crate::{
        network::{
            ExecutionResult, Network, NetworkLeaf, NetworkNode, Sequential, SmallestDegree,
            library::{DummyLibrary, DummyLibraryTensor, panicing::ErroringLibrary},
            store::NetworkStore,
        },
        structure::{
            OrderedStructure,
            representation::{Euclidean, RepName},
        },
        tensors::data::DenseTensor,
    };

    type Tensor = DenseTensor<f64, OrderedStructure<Euclidean>>;
    type Store = NetworkStore<Tensor, f64>;
    type Net = Network<Store, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;
    type FnLib = ErroringLibrary<DummyKey>;

    fn execute(net: &mut Net, lib: &Lib, fn_lib: &FnLib) {
        net.execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(lib, fn_lib)
            .unwrap();
    }

    fn assert_scaled_tensor(net: &Net) {
        let (node, _, _) = net.graph.result().unwrap();
        assert!(
            matches!(node, NetworkNode::Leaf(NetworkLeaf::ScaledTensor(_))),
            "expected executed scalar-tensor product to stay a ScaledTensor, found {node:?}",
        );
    }

    let structure = OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).into_canonical();
    let a = DenseTensor::from_storage_data(vec![1.0, 2.0], structure.clone()).unwrap();
    let b = DenseTensor::from_storage_data(vec![3.0, 4.0], structure).unwrap();
    let lib = Lib::new();
    let fn_lib = FnLib::new();

    let mut left = Net::from_scalar(2.0) * Net::from_tensor(a);
    execute(&mut left, &lib, &fn_lib);
    assert_scaled_tensor(&left);

    let mut right = Net::from_scalar(3.0) * Net::from_tensor(b);
    execute(&mut right, &lib, &fn_lib);
    assert_scaled_tensor(&right);

    let mut sum = left + right;
    execute(&mut sum, &lib, &fn_lib);

    let ExecutionResult::Val(result) = sum.result_tensor::<LibTensor, Lib>(&lib).unwrap() else {
        panic!("expected tensor result");
    };
    assert_eq!(result.data, vec![11.0, 16.0]);
}

#[test]
fn lazy_scalar_tensor_leaves_add_scalars_in_either_order() {
    use crate::{
        network::{
            ExecutionResult, Network, NetworkLeaf, NetworkNode, Sequential, SmallestDegree,
            graph::{ScalarRef, ScaledTensorRef},
            library::{DummyLibrary, DummyLibraryTensor, panicing::ErroringLibrary},
            store::{NetworkStore, TensorScalarStore},
        },
        structure::{OrderedStructure, representation::Euclidean},
        tensors::data::DenseTensor,
    };

    type Tensor = DenseTensor<f64, OrderedStructure<Euclidean>>;
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    let lib = DummyLibrary::<Tensor, DummyKey>::new();
    let functions = ErroringLibrary::<DummyKey>::new();
    let structure = OrderedStructure::new(vec![]).into_canonical();
    let mut tensor =
        Net::from_tensor(DenseTensor::from_storage_data(vec![3.0], structure.clone()).unwrap());
    let second = tensor
        .store
        .add_tensor(DenseTensor::from_storage_data(vec![7.0], structure).unwrap());
    let scale = tensor.store.add_scalar(2.0);
    let other_scale = tensor.store.add_scalar(-5.0);
    tensor.store.scalar_aliases = vec![Some(11.0), Some(-13.0)];

    // Distinct entries and coefficients detect reused references; aliases must
    // resolve through the store without broadcasting scales into tensors.
    for (leaf, expected) in [
        (NetworkLeaf::Scalar(scale.into()), 2.0),
        (NetworkLeaf::Scalar(ScalarRef::Alias(scale)), 11.0),
        (NetworkLeaf::LocalTensor(0), 3.0),
        (NetworkLeaf::TensorSum(vec![0, second]), 10.0),
        (
            NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(0, scale)),
            6.0,
        ),
        (
            NetworkLeaf::ScaledTensor(ScaledTensorRef::tensor(second)),
            7.0,
        ),
        (
            NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled_ref(0, ScalarRef::Alias(scale))),
            33.0,
        ),
        (
            NetworkLeaf::ScaledTensorSum(vec![
                ScaledTensorRef::scaled(0, scale),
                ScaledTensorRef::scaled(second, other_scale),
                ScaledTensorRef::tensor(0),
            ]),
            -26.0,
        ),
        (
            NetworkLeaf::ScaledTensorSum(vec![
                ScaledTensorRef::scaled_ref(0, ScalarRef::Alias(scale)),
                ScaledTensorRef::scaled_ref(second, ScalarRef::Alias(other_scale)),
                ScaledTensorRef::tensor(0),
            ]),
            -55.0,
        ),
    ] {
        let mut left = tensor.clone();
        let root = left.graph.result().unwrap().1;
        left.graph.graph[root] = NetworkNode::Leaf(leaf.clone());
        for mut sum in [
            left.clone() + Net::from_scalar(1.0),
            Net::from_scalar(1.0) + left,
        ] {
            sum.execute::<Sequential, SmallestDegree, DummyLibraryTensor<Tensor>, _, _>(
                &lib, &functions,
            )
            .unwrap();
            let ExecutionResult::Val(actual) = sum.result_scalar().unwrap() else {
                panic!("expected a scalar result for {leaf:?}");
            };
            assert_eq!(*actual, expected + 1.0, "{leaf:?}");
            assert!(matches!(
                sum.graph.result().unwrap().0,
                NetworkNode::Leaf(NetworkLeaf::Scalar(_))
            ));
            assert!(sum.store.tensors.is_empty());
        }
    }
}

#[test]
fn lazy_scalar_tensor_sum_rejects_free_indices() {
    use crate::{
        network::{
            Network, NetworkLeaf,
            graph::ScaledTensorRef,
            store::{NetworkStore, TensorScalarStore},
        },
        structure::{
            OrderedStructure,
            representation::{Euclidean, RepName},
        },
        tensors::data::DenseTensor,
    };

    type Tensor = DenseTensor<f64, OrderedStructure<Euclidean>>;
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    let structure = OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).into_canonical();
    let mut tensor =
        Net::from_tensor(DenseTensor::from_storage_data(vec![1.0, 2.0], structure).unwrap());
    let scale = tensor.store.add_scalar(2.0);
    let root = tensor.graph.result().unwrap().1;
    let scalar = NetworkLeaf::Scalar(scale.into());
    for leaf in [
        NetworkLeaf::LocalTensor(0),
        NetworkLeaf::TensorSum(vec![0, 0]),
        NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(0, scale)),
        NetworkLeaf::ScaledTensorSum(vec![
            ScaledTensorRef::scaled(0, scale),
            ScaledTensorRef::tensor(0),
        ]),
    ] {
        for targets in [
            [(root, &leaf), (root, &scalar)],
            [(root, &scalar), (root, &leaf)],
        ] {
            assert!(
                super::try_balanced_scalar_sum::<DummyKey, super::AbstractIndex, _>(
                    &mut tensor.store,
                    &targets,
                    None,
                )
                .is_none(),
                "a tensor with a free index must not be treated as a scalar: {leaf:?}"
            );
        }
    }

    // Even a malformed lazy sum whose first tensor is scalar must fail the
    // scalar eligibility check when a later tensor carries a free index.
    let scalar_tensor = tensor.store.add_tensor(
        DenseTensor::from_storage_data(vec![3.0], OrderedStructure::new(vec![]).into_canonical())
            .unwrap(),
    );
    for leaf in [
        NetworkLeaf::TensorSum(vec![scalar_tensor, 0]),
        NetworkLeaf::ScaledTensorSum(vec![
            ScaledTensorRef::scaled(scalar_tensor, scale),
            ScaledTensorRef::tensor(0),
        ]),
    ] {
        for targets in [
            [(root, &leaf), (root, &scalar)],
            [(root, &scalar), (root, &leaf)],
        ] {
            assert!(
                super::try_balanced_scalar_sum::<DummyKey, super::AbstractIndex, _>(
                    &mut tensor.store,
                    &targets,
                    None,
                )
                .is_none()
            );
        }
    }
}

#[test]
fn lazy_scalar_tensor_contraction_adds_scalar_in_either_order() {
    use crate::{
        network::{
            ExecutionResult, Network, NetworkLeaf, NetworkNode, Sequential, SmallestDegree,
            graph::ScaledTensorRef,
            library::{DummyLibrary, DummyLibraryTensor, panicing::ErroringLibrary},
            store::{NetworkStore, TensorScalarStore},
        },
        structure::{
            OrderedStructure,
            representation::{Euclidean, RepName},
        },
        tensors::data::DenseTensor,
    };

    type Tensor = DenseTensor<f64, OrderedStructure<Euclidean>>;
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    let lib = DummyLibrary::<Tensor, DummyKey>::new();
    let functions = ErroringLibrary::<DummyKey>::new();
    let structure = OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).into_canonical();
    let mut sum = Net::from_tensor(
        DenseTensor::from_storage_data(vec![1.0, 2.0], structure.clone()).unwrap(),
    );
    let second = sum
        .store
        .add_tensor(DenseTensor::from_storage_data(vec![3.0, 4.0], structure.clone()).unwrap());
    let two = sum.store.add_scalar(2.0);
    let three = sum.store.add_scalar(3.0);
    let root = sum.graph.result().unwrap().1;
    // Construct the deferred 2A + 3B explicitly so the regression does not
    // depend on the eager-sum size threshold or process-global environment.
    sum.graph.graph[root] = NetworkNode::Leaf(NetworkLeaf::ScaledTensorSum(vec![
        ScaledTensorRef::scaled(0, two),
        ScaledTensorRef::scaled(second, three),
    ]));
    let c = Net::from_tensor(DenseTensor::from_storage_data(vec![5.0, 6.0], structure).unwrap());
    let contracted = sum * c;
    for mut sum in [
        contracted.clone() + Net::from_scalar(1.0),
        Net::from_scalar(1.0) + contracted,
    ] {
        sum.execute::<Sequential, SmallestDegree, DummyLibraryTensor<Tensor>, _, _>(
            &lib, &functions,
        )
        .unwrap();
        let ExecutionResult::Val(actual) = sum.result_scalar().unwrap() else {
            panic!("expected a fully contracted scalar");
        };
        assert_eq!(*actual, 152.0);
    }
}

#[test]
fn odd_tensor_powers_keep_the_base_square_fixed() {
    use super::{
        ExecutionResult, Network, NetworkGraph, NetworkLeaf, ScaledTensorRef, Sequential,
        SmallestDegree,
        library::{DummyLibrary, DummyLibraryTensor, panicing::ErroringLibrary},
        store::{NetworkStore, TensorScalarStore},
    };
    use crate::{
        structure::{
            OrderedStructure, ScalarTensor,
            representation::{Euclidean, RepName},
        },
        tensors::data::DenseTensor,
    };

    type Tensor = DenseTensor<f64, OrderedStructure<Euclidean>>;
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;
    type FnLib = ErroringLibrary<DummyKey>;
    let lib = Lib::new();
    let fn_lib = FnLib::new();

    for open_index in [false, true] {
        let mut net = Net::from_tensor(if open_index {
            let structure =
                OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).into_canonical();
            DenseTensor::from_storage_data(vec![1.0, 2.0], structure).unwrap()
        } else {
            Tensor::new_scalar(3.0)
        });
        net.store.add_tensor(net.store.tensors[0].clone());
        let twice = net.store.add_scalar(2.0);
        let thrice = net.store.add_scalar(3.0);
        for (leaf, scale) in [
            (NetworkLeaf::LocalTensor(0), 1.0_f64),
            (NetworkLeaf::TensorSum(vec![0, 1]), 2.0),
            (
                NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(0, twice)),
                2.0,
            ),
            (
                NetworkLeaf::ScaledTensorSum(vec![
                    ScaledTensorRef::scaled(0, twice),
                    ScaledTensorRef::scaled(1, thrice),
                ]),
                5.0,
            ),
        ] {
            net.graph = NetworkGraph::tensor(&net.store.tensors[0], leaf);
            for power in [1, 3, 5] {
                let mut powered = net.clone().pow(power);
                powered
                    .execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
                    .unwrap();
                let ExecutionResult::Val(actual) =
                    powered.result_tensor::<LibTensor, Lib>(&lib).unwrap()
                else {
                    panic!("expected tensor result");
                };
                let expected = if open_index {
                    // Each pair contracts to the squared norm; the last factor
                    // retains the original free index and tensor components.
                    let factor = scale * (5.0 * scale * scale).powi(i32::from(power / 2));
                    vec![factor, 2.0 * factor]
                } else {
                    vec![(3.0 * scale).powi(i32::from(power))]
                };
                assert_eq!(actual.data, expected, "open={open_index}, power={power}");
            }
        }
    }
}

#[cfg(feature = "shadowing")]
#[test]
fn odd_library_tensor_powers_keep_the_base_square_fixed() {
    use symbolica::{
        atom::{Atom, Symbol},
        symbol,
    };

    use super::{
        ExecutionResult, Network, Sequential, SmallestDegree,
        library::{
            panicing::ErroringLibrary,
            symbolic::{ExplicitKey, TensorLibrary},
        },
        store::NetworkStore,
    };
    use crate::{
        structure::{
            NamedStructure, TensorStructure,
            abstract_index::AbstractIndex,
            representation::{Euclidean, LibraryRep, RepName, Representation},
        },
        tensors::data::DataTensor,
    };

    type Tensor = DataTensor<f64, NamedStructure<Symbol, Vec<Atom>, LibraryRep>>;
    type LibTensor = DataTensor<f64, ExplicitKey<AbstractIndex>>;
    type Lib = TensorLibrary<LibTensor, AbstractIndex>;
    type FnLib = ErroringLibrary<DummyKey>;
    type Net = Network<NetworkStore<Tensor, f64>, ExplicitKey<AbstractIndex>, DummyKey>;
    let mut lib = Lib::new();
    let fn_lib = FnLib::new();
    let key = ExplicitKey::from_iter(
        [] as [Representation<LibraryRep>; 0],
        symbol!("odd_power_scalar_tensor"),
        None,
    );
    lib.insert_explicit_dense(key.clone(), vec![3.0]).unwrap();
    let structure = key
        .canonical()
        .clone()
        .reindex_storage(&([] as [AbstractIndex; 0]))
        .unwrap()
        .apply();

    for power in [1, 3, 5] {
        let mut net = Net::library_tensor(&structure, key.clone()).pow(power);
        net.execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
            .unwrap();
        let ExecutionResult::Val(actual) = net.result_tensor::<LibTensor, Lib>(&lib).unwrap()
        else {
            panic!("expected tensor result");
        };
        assert_eq!(
            actual.into_owned().to_bare_dense().data,
            vec![3.0_f64.powi(i32::from(power))]
        );
    }

    let key = ExplicitKey::from_iter(
        [Euclidean {}.new_rep(2).to_lib()],
        symbol!("sum_boundary_library_vector"),
        None,
    );
    lib.insert_explicit_dense(key.clone(), vec![1.0, 2.0])
        .unwrap();
    let structure = key
        .canonical()
        .clone()
        .reindex_storage(&[AbstractIndex::from(1)])
        .unwrap()
        .apply();
    let vector = Net::library_tensor(&structure, key);
    let original = vector.clone() * (vector.clone() + vector);
    let mut prepared = original.clone();
    prepared.graph.contract_ready_sum_boundaries();
    for mut network in [original, prepared] {
        network
            .execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
            .unwrap();
        let ExecutionResult::Val(actual) = network.result_scalar().unwrap() else {
            panic!("expected the contracted library-vector sum")
        };
        assert_eq!(*actual, 10.0);
    }
}

#[test]
fn lazy_scalar_tensors_add_scalars_in_either_order() {
    use super::{
        ExecutionResult, Network, NetworkGraph, NetworkLeaf, ScaledTensorRef, Sequential,
        SmallestDegree,
        library::{DummyLibrary, DummyLibraryTensor, panicing::ErroringLibrary},
        store::{NetworkStore, TensorScalarStore},
    };
    use crate::{
        structure::{
            OrderedStructure, ScalarTensor,
            representation::{Euclidean, RepName},
        },
        tensors::data::DenseTensor,
    };

    type Tensor = DenseTensor<f64, OrderedStructure<Euclidean>>;
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;
    type FnLib = ErroringLibrary<DummyKey>;
    let lib = Lib::new();
    let fn_lib = FnLib::new();

    for open_index in [false, true] {
        let mut lazy = Net::from_tensor(if open_index {
            let structure =
                OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).into_canonical();
            DenseTensor::from_storage_data(vec![3.0, 5.0], structure).unwrap()
        } else {
            Tensor::new_scalar(3.0)
        });
        lazy.store.add_tensor(if open_index {
            lazy.store.tensors[0].clone()
        } else {
            Tensor::new_scalar(5.0)
        });
        let negative = lazy.store.add_scalar(-2.0);
        let positive = lazy.store.add_scalar(4.0);
        for (leaf, expected) in [
            (NetworkLeaf::TensorSum(vec![0, 1]), 8.0),
            (
                NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(0, negative)),
                -6.0,
            ),
            (
                NetworkLeaf::ScaledTensorSum(vec![
                    ScaledTensorRef::scaled(0, negative),
                    ScaledTensorRef::scaled(1, positive),
                ]),
                14.0,
            ),
        ] {
            // These are the lazy leaves produced by contraction. Keep them
            // lazy at the next sum boundary, where operand order used to
            // decide whether a valid scalar-valued tensor was accepted.
            lazy.graph = NetworkGraph::tensor(&lazy.store.tensors[0], leaf.clone());
            if !open_index {
                let ExecutionResult::Val(actual) = lazy.result_scalar().unwrap() else {
                    panic!("expected scalar-valued lazy tensor");
                };
                assert_eq!(*actual, expected);
            }
            for scalar_first in [false, true] {
                if open_index {
                    // A well-formed network rejects incompatible free indices
                    // before execution. The public scalar accessor must likewise
                    // reject an open lazy tensor without altering its value.
                    assert!(matches!(
                        lazy.result_scalar(),
                        Err(super::TensorNetworkError::NoScalar)
                    ));
                    let ExecutionResult::Val(actual) =
                        lazy.result_tensor::<LibTensor, Lib>(&lib).unwrap()
                    else {
                        panic!("expected an unchanged open tensor")
                    };
                    let factor = match &leaf {
                        NetworkLeaf::TensorSum(_) => 2.0,
                        NetworkLeaf::ScaledTensor(_) => -2.0,
                        NetworkLeaf::ScaledTensorSum(_) => 2.0,
                        _ => unreachable!(),
                    };
                    assert_eq!(actual.data, vec![3.0 * factor, 5.0 * factor]);
                    continue;
                }
                let scalar = Net::from_scalar(7.0);
                let mut sum = if scalar_first {
                    scalar + lazy.clone()
                } else {
                    lazy.clone() + scalar
                };
                sum.execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
                    .unwrap();
                let ExecutionResult::Val(actual) = sum.result_scalar().unwrap() else {
                    panic!("expected scalar sum");
                };
                assert_eq!(*actual, expected + 7.0);
            }
        }
    }
}

#[cfg(feature = "shadowing")]
#[test]
fn bulk_atom_sum_accepts_closed_tensor_leaves_and_preserves_aliases() {
    use linnet::half_edge::{builder::HedgeGraphBuilder, involution::Flow};
    use std::collections::HashMap;
    use symbolica::{
        atom::{Atom, AtomCore},
        parse,
    };

    use super::{
        ExecutionResult, Network, NetworkGraph, NetworkLeaf, ScaledTensorRef, Sequential,
        SmallestDegree,
        graph::{NetworkEdge, NetworkNode, NetworkOp, ScalarRef},
        library::{DummyLibrary, DummyLibraryTensor, panicing::ErroringLibrary},
        store::{NetworkStore, TensorScalarStore},
        tags::scalar_store_alias,
    };
    use crate::{
        structure::{
            OrderedStructure, ScalarTensor,
            representation::{Euclidean, RepName},
        },
        tensors::data::{DataTensor, DenseTensor, SparseTensor},
    };

    type Tensor = DataTensor<Atom, OrderedStructure<Euclidean>>;
    type Net = Network<NetworkStore<Tensor, Atom>, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;
    type FnLib = ErroringLibrary<DummyKey>;
    let lib = Lib::new();
    let fn_lib = FnLib::new();
    let mut initial: NetworkStore<Tensor, Atom> = NetworkStore::default();
    let value = parse!("(x+y)^3*(z+w)^2");
    let positive = initial.add_scalar(value.clone());
    let negative = initial.add_scalar(-&value);
    let factor = initial.add_scalar(parse!("(a+b)^2"));
    let zero = initial.add_scalar(Atom::Zero);
    initial
        .scalar_aliases
        .resize_with(initial.scalar.len(), || None);
    initial.scalar_aliases[factor] = Some(scalar_store_alias(factor));
    let tensor = initial.add_tensor(Tensor::new_scalar(value.clone()));
    let opposite = initial.add_tensor(Tensor::new_scalar(-value.clone()));
    let sparse_zero = initial.add_tensor(DataTensor::Sparse(SparseTensor {
        elements: HashMap::new(),
        zero: Atom::Zero,
        structure: OrderedStructure::new(vec![]).into_canonical(),
    }));
    let aliased_scale = ScaledTensorRef {
        tensor,
        scale: Some(ScalarRef::Alias(factor)),
    };
    let leaves: Vec<NetworkLeaf<DummyKey>> = vec![
        NetworkLeaf::Scalar(positive.into()),
        NetworkLeaf::Scalar(negative.into()),
        NetworkLeaf::LocalTensor(tensor),
        NetworkLeaf::LocalTensor(sparse_zero),
        NetworkLeaf::TensorSum(vec![tensor, opposite]),
        NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(tensor, factor)),
        NetworkLeaf::ScaledTensor(aliased_scale.clone()),
        NetworkLeaf::ScaledTensorSum(vec![
            aliased_scale,
            ScaledTensorRef::scaled(opposite, factor),
        ]),
        NetworkLeaf::ScaledTensorSum(vec![]),
        NetworkLeaf::Scalar(ScalarRef::Alias(factor)),
    ];

    // Small and large complete sums cover scalar aliases, exact cancellation,
    // sparse zeros and factored scales through ordinary network execution.
    let scale = parse!("(a+b)^2");
    let alias = scalar_store_alias(factor);
    let expected_terms = [
        value.clone(),
        -&value,
        value.clone(),
        Atom::Zero,
        Atom::Zero,
        &value * &scale,
        &value * &alias,
        &value * (&alias - &scale),
        Atom::Zero,
        alias,
    ];
    for count in [3, 4, 31, 32] {
        for reversed in [false, true] {
            let mut terms = (0..count)
                .map(|index| index % leaves.len())
                .collect::<Vec<_>>();
            if reversed {
                terms.reverse();
            }
            let expected = Atom::add_many(
                terms
                    .iter()
                    .map(|&index| expected_terms[index].clone())
                    .collect::<Vec<_>>(),
            );
            let mut builder = HedgeGraphBuilder::new();
            let root = builder.add_node(NetworkNode::Op(NetworkOp::Sum));
            builder.add_external_edge(root, NetworkEdge::Head, true, Flow::Source);
            for index in terms {
                let child = builder.add_node(NetworkNode::Leaf(leaves[index].clone()));
                builder.add_edge(child, root, NetworkEdge::Head, true);
            }
            let mut network = Net::from_scalar(Atom::Zero);
            network.store = initial.clone();
            network.graph = builder.into();
            network
                .execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
                .unwrap();
            let ExecutionResult::Val(actual) = network.result_scalar().unwrap() else {
                panic!("expected the complete scalar sum")
            };
            assert!((actual.into_owned() - expected).expand().is_zero());
        }
    }

    // Even a zero-scaled open tensor remains ineligible. The public scalar
    // result must reject the tensor rather than expose a partial conversion.
    let open = initial.add_tensor(DataTensor::Dense(
        DenseTensor::from_storage_data(
            vec![Atom::Zero; 2],
            OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).into_canonical(),
        )
        .unwrap(),
    ));
    let open_index = open;
    let open = NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(open_index, zero));
    let mut network = Net::from_scalar(Atom::Zero);
    network.store = initial.clone();
    network.graph = NetworkGraph::tensor(&initial.tensors[open_index], open);
    assert!(matches!(
        network.result_scalar(),
        Err(super::TensorNetworkError::NoScalar)
    ));

    let empty = NetworkLeaf::<DummyKey>::TensorSum(vec![]);
    for count in [4, 32] {
        for zero_leaves in [[&empty, &empty], [&leaves[0], &leaves[1]]] {
            let mut builder = HedgeGraphBuilder::new();
            let root = builder.add_node(NetworkNode::Op(NetworkOp::Sum));
            builder.add_external_edge(root, NetworkEdge::Head, true, Flow::Source);
            for index in 0..count {
                let child = builder.add_node(NetworkNode::Leaf(zero_leaves[index % 2].clone()));
                builder.add_edge(child, root, NetworkEdge::Head, true);
            }
            let mut network = Net::from_scalar(Atom::Zero);
            network.store = initial.clone();
            network.graph = builder.into();
            network
                .execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
                .unwrap();
            match network.result_scalar().unwrap() {
                ExecutionResult::Zero => {}
                ExecutionResult::Val(actual) => assert_eq!(*actual, Atom::Zero),
                ExecutionResult::One => panic!("expected zero"),
            }
        }
    }
}

#[cfg(feature = "shadowing")]
#[test]
fn fused_numeric_contraction_matches_dual_indices_in_storage_order() {
    use symbolica::{atom::Atom, parse};

    use super::{FastTensorSumContract, fused_numeric_tensor_sum_contract};
    use crate::{
        structure::{
            HasStructure, OrderedStructure, TensorDataLayout,
            representation::{LibraryRep, Lorentz, RepName},
            slot::{DualSlotTo, IsAbstractSlot, Slot},
        },
        tensors::{
            data::{SetTensorData, SparseTensor},
            parametric::ParamTensor,
        },
    };

    let tensor = |slots: Vec<Slot<LibraryRep>>, coefficient: Atom| {
        let coordinates = (0..slots.len()).collect::<Vec<_>>();
        let structure = OrderedStructure::new(slots);
        let layout = TensorDataLayout::from_canonicalized(&structure).unwrap();
        let mut tensor = SparseTensor::empty(structure.into_canonical(), Atom::Zero);
        tensor
            .set(
                &layout
                    .logical_expanded_to_storage_expanded(&coordinates)
                    .unwrap(),
                coefficient,
            )
            .unwrap();
        ParamTensor::composite(tensor.into())
    };

    // Canonical storage groups base and dual slots separately. Closing these
    // indices requires a transposition or a three-cycle, not equal positions.
    for rank in [2, 3] {
        let slots = (0..rank)
            .map(|i| {
                if i + 1 == rank {
                    Lorentz {}.dual().new_slot(3, i).to_lib()
                } else {
                    Lorentz {}.new_slot(3, i).to_lib()
                }
            })
            .collect::<Vec<_>>();
        let numeric = tensor(slots.iter().map(|slot| slot.dual()).collect(), Atom::num(3));
        for factor in [Atom::num(1), parse!("x+y")] {
            let term = tensor(slots.clone(), Atom::num(2) * &factor);
            for count in [1, 2, 8] {
                let terms = vec![&term; count];
                for terms_on_left in [false, true] {
                    let result = fused_numeric_tensor_sum_contract(&terms, &numeric, terms_on_left)
                        .expect("sparse numeric contraction uses the fused path")
                        .unwrap();
                    let actual = match result {
                        FastTensorSumContract::Materialized(tensor) => tensor.scalar().unwrap(),
                        FastTensorSumContract::Terms(terms) => {
                            Atom::add_many(terms.into_iter().map(|tensor| tensor.scalar().unwrap()))
                        }
                        FastTensorSumContract::ScaledTerms(terms) => {
                            Atom::add_many(terms.into_iter().map(|term| {
                                term.tensor.scalar().unwrap() * term.scale.unwrap_or(Atom::num(1))
                            }))
                        }
                    };
                    assert_eq!(
                        actual,
                        Atom::num(6 * count as i64) * &factor,
                        "rank={rank}, count={count}, terms_on_left={terms_on_left}"
                    );
                }
            }
        }
    }
}

#[cfg(feature = "shadowing")]
#[test]
fn sparse_pair_estimate_counts_shared_output_coordinates() {
    use linnet::permutation::Permutation;
    use symbolica::atom::Atom;

    use super::FastTensorSumContractible;
    use crate::{
        structure::{
            OrderedStructure, TensorStructure,
            representation::{Euclidean, RepName},
        },
        tensors::{
            data::{DataTensor, SparseTensor},
            parametric::ParamTensor,
        },
    };

    let structure: OrderedStructure<Euclidean> = OrderedStructure::new(vec![
        Euclidean {}.new_slot(2, 1),
        Euclidean {}.new_slot(2, 2),
    ])
    .into_canonical();
    let tensor = |support: &[[usize; 2]]| {
        ParamTensor::composite(DataTensor::Sparse(SparseTensor {
            elements: support
                .iter()
                .map(|coordinates| (structure.flat_index(*coordinates).unwrap(), Atom::num(1)))
                .collect(),
            zero: Atom::Zero,
            structure: structure.clone(),
        }))
    };
    let left = tensor(&[[0, 0], [0, 1], [1, 0]]);
    let right = tensor(&[[0, 0], [1, 0], [1, 1]]);

    // Contracting the first axis produces (0,0) twice, from distinct contracted
    // groups, and (1,0)/(0,1) once each. Full contraction joins two coordinates;
    // both products have the same pair of empty free-coordinate keys.
    for (matches, limit, products, entries, max_products) in [
        ([true, false], 4, 4, 3, 2),
        ([true, false], 3, 4, 4, 2),
        ([true, true], 2, 2, 1, 2),
    ] {
        let matched_axes = matches.iter().filter(|matched| **matched).count();
        let output_dense_size = if matched_axes == 2 { 1 } else { 4 };
        let actual = left.contraction_pair_estimate(
            &right,
            &Permutation::id(matched_axes),
            &matches,
            &matches,
            left.contraction_profile(),
            right.contraction_profile(),
            output_dense_size,
            limit,
        );
        assert_eq!(actual.estimated_products, products);
        assert_eq!(actual.estimated_output_entries, entries);
        assert_eq!(actual.output_dense_size, output_dense_size);
        assert_eq!(actual.max_output_entry_products, max_products);
    }

    // These supports overlap only after reordering the matched coordinates.
    let left = tensor(&[[0, 1]]);
    let right = tensor(&[[1, 0]]);
    for (permutation, products) in [
        (Permutation::id(2), 0),
        (Permutation::from_map(vec![1, 0]), 1),
    ] {
        let estimate = left.contraction_pair_estimate(
            &right,
            &permutation,
            &[true, true],
            &[true, true],
            left.contraction_profile(),
            right.contraction_profile(),
            1,
            1,
        );
        assert_eq!(estimate.estimated_products, products);
    }
}

#[cfg(feature = "shadowing")]
#[test]
fn large_scaled_tensor_sum_preserves_results_across_strategies() {
    use std::collections::HashMap;

    use symbolica::{
        atom::{Atom, AtomCore},
        function, symbol,
    };

    use crate::{
        network::{
            ExecutionResult, MinIntermediateCost, MinResultRank, Network, Sequential,
            SequentialRef, TensorOrScalarOrKey,
            library::{DummyLibrary, DummyLibraryTensor, panicing::ErroringLibrary},
            store::NetworkStore,
        },
        structure::{
            OrderedStructure, TensorStructure,
            concrete_index::FlatIndex,
            representation::{Euclidean, RepName},
        },
        tensors::{
            data::{DataTensor, GetTensorData, SparseTensor},
            parametric::ParamTensor,
        },
    };

    type Tensor = ParamTensor<OrderedStructure<Euclidean>>;
    type Store = NetworkStore<Tensor, Atom>;
    type Net = Network<Store, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;
    type FnLib = ErroringLibrary<DummyKey>;

    macro_rules! check_method {
        ($method:ty) => {{
            let lib = Lib::new();
            let fn_lib = FnLib::new();
            let coefficient = Atom::add_many(
                (0..2048)
                    .map(|i| function!(symbol!("coefficient"), Atom::num(i)))
                    .collect::<Vec<_>>(),
            );
            let other_coefficient = function!(symbol!("other_coefficient"), &coefficient);

            for width in [2usize, 4096] {
                let structure =
                    OrderedStructure::new(vec![Euclidean {}.new_slot(width, 1)]).into_canonical();
                let a = (0..width).map(|i| i as i64 % 3 - 1).collect::<Vec<_>>();
                let b = (0..width).map(|i| i as i64 % 5 - 2).collect::<Vec<_>>();
                let weights = (0..width).map(|i| i as i64 % 7 - 3).collect::<Vec<_>>();
                let tensor = |entries: &[i64]| {
                    ParamTensor::composite(DataTensor::Sparse(SparseTensor {
                        elements: entries
                            .iter()
                            .enumerate()
                            .filter(|(_, value)| **value != 0)
                            .map(|(i, &value)| (FlatIndex::from(i), Atom::num(value)))
                            .collect::<HashMap<_, _>>(),
                        zero: Atom::Zero,
                        structure: structure.clone(),
                    }))
                };
                let sum = Net::from_scalar(coefficient.clone()) * Net::from_tensor(tensor(&a))
                    + Net::from_scalar(other_coefficient.clone()) * Net::from_tensor(tensor(&b));
                // Tiny and large sums share the same outward contract even when their
                // intermediate execution chooses different materialization strategies.

                // The independent reference contracts integer vectors first and only
                // then multiplies their dot products by the symbolic coefficients.
                let dot_a: i64 = a.iter().zip(&weights).map(|(a, w)| a * w).sum();
                let dot_b: i64 = b.iter().zip(&weights).map(|(b, w)| b * w).sum();
                let expected =
                    &coefficient * Atom::num(dot_a) + &other_coefficient * Atom::num(dot_b);
                let mut contracted = sum.clone() * Net::from_tensor(tensor(&weights));
                contracted
                    .execute::<Sequential, $method, LibTensor, Lib, FnLib>(&lib, &fn_lib)
                    .unwrap();
                let ExecutionResult::Val(actual) = contracted.result_scalar().unwrap() else {
                    panic!("expected a scalar contraction");
                };
                assert!((actual.as_ref() - &expected).expand().is_zero());
                drop(contracted);

                // A terminal sum must still expose its tensor through result(), even
                // when a negation wraps it. These are intentionally complete graphs;
                // callers do not have to request materialization with result_tensor().
                for (mut terminal, sign) in [(sum.clone(), 1), (-sum.clone(), -1)] {
                    terminal
                        .execute::<SequentialRef, $method, LibTensor, Lib, FnLib>(&lib, &fn_lib)
                        .unwrap();
                    let ExecutionResult::Val(TensorOrScalarOrKey::Tensor { tensor, .. }) =
                        terminal.result().unwrap()
                    else {
                        panic!("expected a materialized terminal tensor");
                    };
                    for i in [0, width / 2, width - 1] {
                        let actual = tensor.tensor.get_ref_linear(FlatIndex::from(i)).unwrap();
                        let expected = &coefficient * Atom::num(sign * a[i])
                            + &other_coefficient * Atom::num(sign * b[i]);
                        assert!((actual - expected).expand().is_zero());
                    }
                }

                // Leave one index open after contraction to cover a terminal product,
                // including the parallel execution boundary.
                let matrix_structure = OrderedStructure::new(vec![
                    Euclidean {}.new_slot(width, 1),
                    Euclidean {}.new_slot(2, 2),
                ]);
                let matrix = ParamTensor::composite(DataTensor::Sparse(SparseTensor {
                    elements: weights
                        .iter()
                        .enumerate()
                        .filter(|(_, weight)| **weight != 0)
                        .flat_map(|(i, &weight)| {
                            let structure = &matrix_structure;
                            (0..2).map(move |j| {
                                // OrderedStructure can reorder unequal dimensions;
                                // apply the same permutations to the data coordinates.
                                let mut indices = [i, j];
                                let canonical_indices =
                                    structure.layout().logical_to_canonical(&indices);
                                indices.copy_from_slice(&canonical_indices);
                                (
                                    structure.canonical().flat_index(indices).unwrap(),
                                    Atom::num(weight * (j as i64 + 1)),
                                )
                            })
                        })
                        .collect::<HashMap<_, _>>(),
                    zero: Atom::Zero,
                    structure: matrix_structure.into_canonical(),
                }));
                let mut product = sum * Net::from_tensor(matrix);
                product
                    .execute_parallel::<$method, LibTensor, Lib, FnLib>(&lib, &fn_lib)
                    .unwrap();
                let ExecutionResult::Val(TensorOrScalarOrKey::Tensor { tensor, .. }) =
                    product.result().unwrap()
                else {
                    panic!("expected a materialized nonscalar product");
                };
                for j in 0..2 {
                    let actual = tensor.tensor.get_ref_linear(FlatIndex::from(j)).unwrap();
                    assert!(
                        (actual - &expected * Atom::num(j as i64 + 1))
                            .expand()
                            .is_zero()
                    );
                }
            }
        }};
    }
    check_method!(MinResultRank);
    check_method!(MinIntermediateCost);
}

#[test]
fn contraction_strategies_match_independent_coordinate_sum() {
    use crate::{
        network::{
            ExecutionResult, MinIntermediateCost, MinResultRank, Network, Sequential,
            library::{DummyLibrary, DummyLibraryTensor, panicing::ErroringLibrary},
            store::NetworkStore,
        },
        structure::{
            OrderedStructure, TensorStructure,
            representation::{Euclidean, RepName},
        },
        tensors::data::DenseTensor,
    };

    type Tensor = DenseTensor<f64, OrderedStructure<Euclidean>>;
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;
    type FnLib = ErroringLibrary<DummyKey>;

    // Supply values in the written index order, independently of the storage
    // permutations imposed by OrderedStructure for unequal dimensions.
    fn tensor(axes: &[(usize, usize)], value: impl Fn(&[usize]) -> i64) -> Tensor {
        let ordered = OrderedStructure::new(
            axes.iter()
                .map(|&(dim, id)| Euclidean {}.new_slot(dim, id))
                .collect(),
        );
        let size = axes.iter().map(|(dim, _)| dim).product();
        let mut data = vec![0.0; size];
        for linear in 0..size {
            let mut remainder = linear;
            let mut indices = vec![0; axes.len()];
            for (axis, &(dimension, _)) in axes.iter().enumerate().rev() {
                indices[axis] = remainder % dimension;
                remainder /= dimension;
            }
            let entry = value(&indices);
            let canonical_indices = ordered.layout().logical_to_canonical(&indices);
            indices.copy_from_slice(&canonical_indices);
            let flat: usize = ordered.canonical().flat_index(indices).unwrap().into();
            data[flat] = entry as f64;
        }
        DenseTensor::from_storage_data(data, ordered.into_canonical()).unwrap()
    }

    let a = |i: usize, j: usize, x: usize| (1 + i + 2 * j + x % 3) as i64;
    let b = |i: usize, j: usize, y: usize| (2 * i + j + y % 5) as i64 - 3;
    let c = |x: usize, z: usize| (x % 7 + z) as i64 - 2;
    let original = Net::from_tensor(tensor(&[(2, 1), (3, 2), (16, 3)], |v| a(v[0], v[1], v[2])))
        * Net::from_tensor(tensor(&[(2, 1), (3, 2), (16, 4)], |v| b(v[0], v[1], v[2])))
        * Net::from_tensor(tensor(&[(16, 3), (1, 5)], |v| c(v[0], v[1])));

    // Independent coordinate sum: R_yz = sum_ijx A_ijx B_ijy C_xz.
    // Build the oracle in integers; every input and intermediate is exactly
    // representable by this network's f64 backend. No floating tolerance or
    // contraction kernel is used to manufacture the expected result.
    let expected = tensor(&[(16, 4), (1, 5)], |v| {
        let mut sum = 0;
        for i in 0..2 {
            for j in 0..3 {
                for x in 0..16 {
                    sum += a(i, j, x) * b(i, j, v[0]) * c(x, v[1]);
                }
            }
        }
        sum
    });

    let lib = Lib::new();
    let fn_lib = FnLib::new();
    let mut rank_first = original.clone();
    rank_first
        .execute::<Sequential, MinResultRank, LibTensor, Lib, FnLib>(&lib, &fn_lib)
        .unwrap();
    let mut intermediate_cost = original;
    intermediate_cost
        .execute::<Sequential, MinIntermediateCost, LibTensor, Lib, FnLib>(&lib, &fn_lib)
        .unwrap();

    for net in [&rank_first, &intermediate_cost] {
        let ExecutionResult::Val(actual) = net.result_tensor::<LibTensor, Lib>(&lib).unwrap()
        else {
            panic!("expected the tensor with open y,z indices");
        };
        assert_eq!(actual.structure, expected.structure);
        assert_eq!(actual.data, expected.data);
    }

    // AB has rank two but 16*16 entries. AC has rank three with only
    // 2*3*1 entries; its final contraction has 16 entries. Both possible
    // contraction orders must produce the complete coordinate sum above.
}

#[test]
fn ready_sum_boundary_closure_matches_coordinates_across_execution_strategies() {
    use crate::{
        network::{
            ExecutionResult, Network, Sequential, SequentialExtract, SequentialRef, SmallestDegree,
            graph::{NAdd, NMul},
            library::{DummyLibrary, DummyLibraryTensor, panicing::ErroringLibrary},
            store::NetworkStore,
        },
        structure::{
            OrderedStructure, TensorStructure,
            representation::{Euclidean, Lorentz, RepName},
            slot::IsAbstractSlot,
        },
        tensors::data::DenseTensor,
    };

    type Tensor = DenseTensor<f64, OrderedStructure>;
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;
    type FnLib = ErroringLibrary<DummyKey>;

    let tensor = |axes: &[usize], offset: usize| {
        let ordered = OrderedStructure::new(
            axes.iter()
                .map(|axis| Euclidean {}.new_slot(2, *axis).to_lib())
                .collect(),
        );
        let mut data = vec![0.0; 1 << axes.len()];
        for linear in 0..data.len() {
            let mut indices: Vec<_> = (0..axes.len())
                .map(|axis| (linear >> (axes.len() - axis - 1)) & 1)
                .collect();
            let value = (offset + linear) as f64;
            let canonical_indices = ordered.layout().logical_to_canonical(&indices);
            indices.copy_from_slice(&canonical_indices);
            let flat: usize = ordered.canonical().flat_index(indices).unwrap().into();
            data[flat] = value;
        }
        Net::from_tensor(DenseTensor::from_storage_data(data, ordered.into_canonical()).unwrap())
    };
    // Each native arm is a product with four exposed indices. Its matrix factors
    // stay intact while the four external vectors close it. Reversed written
    // matrix indices exercise slot permutations while closing tensor boundaries.
    let left = tensor(&[2, 1], 1).n_mul([tensor(&[3, 4], 2)]);
    let right = tensor(&[1, 2], 3).n_mul([tensor(&[4, 3], 4)]);
    let sum = left.n_add([right]);
    let original = sum.clone().n_mul([
        tensor(&[1], 1),
        tensor(&[2], 2),
        tensor(&[3], 3),
        tensor(&[4], 4),
        Net::from_scalar(7.0),
    ]);
    let reordered = tensor(&[1], 1).n_mul([
        tensor(&[2], 2),
        tensor(&[3], 3),
        tensor(&[4], 4),
        Net::from_scalar(7.0),
        sum,
    ]);
    let mut expected = 0usize;
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                for l in 0..2 {
                    expected += 7
                        * (1 + i)
                        * (2 + j)
                        * (3 + k)
                        * (4 + l)
                        * ((1 + 2 * j + i) * (2 + 2 * k + l) + (3 + 2 * i + j) * (4 + 2 * l + k));
                }
            }
        }
    }
    let lib = Lib::new();
    let fn_lib = FnLib::new();
    let mut prepared = original.clone();
    prepared.graph.contract_ready_sum_boundaries();
    prepared.graph.contract_ready_sum_boundaries();
    // With closing factors first, the original Sum endpoints have opposite
    // graph flows. Their residual seam must survive every partial closure.
    let mut reordered_prepared = reordered.clone();
    reordered_prepared.graph.contract_ready_sum_boundaries();
    assert_eq!(reordered_prepared.graph.n_dangling(), 0);

    macro_rules! check {
        ($strategy:ty) => {{
            for mut net in [
                original.clone(),
                prepared.clone(),
                reordered.clone(),
                reordered_prepared.clone(),
            ] {
                net.execute::<$strategy, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
                    .unwrap();
                let ExecutionResult::Val(actual) = net.result_scalar().unwrap() else {
                    panic!("expected a closed scalar");
                };
                assert_eq!(*actual, expected as f64);
            }
        }};
    }
    check!(Sequential);
    check!(SequentialRef);
    check!(SequentialExtract);
    prepared
        .execute_parallel::<SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
        .unwrap();
    let ExecutionResult::Val(actual) = prepared.result_scalar().unwrap() else {
        panic!("expected a closed scalar");
    };
    assert_eq!(*actual, expected as f64);

    // Equal slot labels in independent scopes must never identify their seams.
    let mut disjoint = (0..8).map(|offset| {
        tensor(&[1], offset + 1)
            .n_add([tensor(&[1], offset + 3)])
            .n_mul([tensor(&[1], 2)])
    });
    let disjoint = disjoint.next().unwrap().n_add(disjoint);
    let disjoint_expected: usize = (0..8)
        .flat_map(|offset| (0..2).map(move |i| (2 * offset + 4 + 2 * i) * (2 + i)))
        .sum();

    // Two independent closures keep two independent sums; neither additive
    // operand becomes a Cartesian product of the other sum's arms.
    // All four Sum targets share one Product and should be prepared together.
    let siblings = Net::from_scalar(1.0).n_mul((1..=4).flat_map(|axis| {
        [
            tensor(&[axis], 4).n_add([tensor(&[axis], 6)]),
            tensor(&[axis], 2),
        ]
    }));
    let sibling_expected: usize = (0..2).map(|i| (10 + 2 * i) * (2 + i)).sum();

    // The two ready vectors close together, while the remaining index
    // still connects the original two sums, without a Cartesian product.
    // The b seam couples two sums while a and c close independently. Neither
    // the shared residual edge nor its axis order may be reconstructed by name.
    let coupled = tensor(&[2, 1], 1).n_add([tensor(&[1, 2], 3)]).n_mul([
        tensor(&[2, 3], 4).n_add([tensor(&[3, 2], 6)]),
        tensor(&[1], 1),
        tensor(&[3], 2),
    ]);
    let coupled_reordered = tensor(&[2, 3], 4).n_add([tensor(&[3, 2], 6)]).n_mul([
        tensor(&[1], 1),
        tensor(&[3], 2),
        tensor(&[2, 1], 1).n_add([tensor(&[1, 2], 3)]),
    ]);
    let mut coupled_expected = 0usize;
    for a in 0..2 {
        for b in 0..2 {
            for c in 0..2 {
                coupled_expected += (1 + a) * (2 + c) * (4 + 3 * a + 3 * b) * (10 + 3 * b + 3 * c);
            }
        }
    }

    // Closing the outer b boundary exposes another ready Sum in the next wave.
    let nested = tensor(&[2, 1], 1)
        .n_add([tensor(&[1, 2], 3)])
        .n_mul([tensor(&[1], 1)])
        .n_add([tensor(&[2], 5)])
        .n_mul([tensor(&[2], 2)]);
    let nested_expected: usize = (0..2)
        .map(|b| (2 + b) * ((0..2).map(|a| (1 + a) * (4 + 3 * a + 3 * b)).sum::<usize>() + 5 + b))
        .sum();

    // A ready leaf can also carry an internal trace; clone its two half-edges
    // and their axis positions, not just its free boundary to the Sum.
    let traced = tensor(&[2], 4)
        .n_add([tensor(&[2], 6)])
        .n_mul([tensor(&[9, 9, 2], 3)]);
    let traced_expected: usize = (0..2).map(|i| (10 + 2 * i) * (12 + 2 * i)).sum();

    for (name, original, expected) in [
        ("disjoint", disjoint, disjoint_expected),
        ("siblings", siblings, sibling_expected.pow(4)),
        ("coupled", coupled, coupled_expected),
        ("coupled reordered", coupled_reordered, coupled_expected),
        ("nested", nested, nested_expected),
        ("traced", traced, traced_expected),
    ] {
        let mut prepared = original.clone();
        prepared.graph.contract_ready_sum_boundaries();
        prepared.graph.graph.check().unwrap();
        assert_eq!(prepared.graph.n_dangling(), 0, "{name}");
        prepared.graph.contract_ready_sum_boundaries();
        for mut net in [original, prepared] {
            net.execute::<SequentialRef, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
                .unwrap();
            let ExecutionResult::Val(actual) = net.result_scalar().unwrap() else {
                panic!("expected a closed scalar for {name}");
            };
            assert_eq!(*actual, expected as f64, "{name}");
        }
    }

    let lorentz = OrderedStructure::new(vec![Lorentz {}.new_slot(2, 9).to_lib()]).into_canonical();
    let dual =
        OrderedStructure::new(vec![Lorentz {}.dual().new_slot(2, 9).to_lib()]).into_canonical();
    let left_dual =
        Net::from_tensor(DenseTensor::from_storage_data(vec![2.0, 3.0], lorentz.clone()).unwrap());
    let right_dual =
        Net::from_tensor(DenseTensor::from_storage_data(vec![5.0, 7.0], dual.clone()).unwrap());
    let dual_spectators = left_dual.clone() * right_dual.clone();
    let open_sum = tensor(&[1, 2], 1).n_add([tensor(&[1, 2], 3)]);
    let partially_closed =
        open_sum.n_mul([tensor(&[1], 1), Net::from_scalar(7.0), dual_spectators]);
    let open_product = tensor(&[1, 2], 1).n_mul([tensor(&[1], 1).n_add([tensor(&[1], 3)])]);
    let scalar_spectator = Net::from_scalar(7.0).n_mul([tensor(&[2], 1).n_add([tensor(&[2], 3)])]);
    let expected_structure =
        OrderedStructure::new(vec![Euclidean {}.new_slot(2, 2).to_lib()]).into_canonical();
    for (name, original, expected) in [
        (
            "dual spectators and partial closure",
            partially_closed,
            vec![4340.0, 5642.0],
        ),
        ("partly open closing factor", open_product, vec![22.0, 32.0]),
        ("scalar spectator", scalar_spectator, vec![28.0, 42.0]),
    ] {
        let mut prepared = original.clone();
        prepared.graph.contract_ready_sum_boundaries();
        prepared.graph.contract_ready_sum_boundaries();
        for mut network in [original, prepared] {
            network
                .execute::<SequentialRef, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
                .unwrap();
            let ExecutionResult::Val(actual) =
                network.result_tensor::<LibTensor, Lib>(&lib).unwrap()
            else {
                panic!("expected the complete open tensor for {name}")
            };
            assert_eq!(actual.structure, expected_structure, "{name}");
            assert_eq!(actual.data, expected, "{name}");
        }
    }

    // A moved dualizable boundary needs endpoint metadata beyond the shared
    // edge descriptor; retain this valid contraction for ordinary execution.
    let dual_sum = left_dual.clone()
        * (right_dual
            + Net::from_tensor(
                DenseTensor::from_storage_data(vec![11.0, 13.0], dual.clone()).unwrap(),
            ));
    let additive =
        tensor(&[1], 1).n_add([tensor(&[1], 3)]) * tensor(&[1], 1).n_add([tensor(&[1], 3)]);
    for (name, original, expected) in [
        ("dual boundary", dual_sum, 92.0),
        ("additive closing factor", additive, 52.0),
    ] {
        let mut prepared = original.clone();
        prepared.graph.contract_ready_sum_boundaries();
        for mut network in [original, prepared] {
            network
                .execute::<SequentialRef, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
                .unwrap();
            let ExecutionResult::Val(actual) = network.result_scalar().unwrap() else {
                panic!("expected the complete scalar for {name}")
            };
            assert_eq!(*actual, expected, "{name}");
        }
    }

    // An unconsumed dualizable boundary also needs endpoint restoration.
    // Retain both original operand orders without cutting their shared edge.
    let right = Net::from_tensor(DenseTensor::from_storage_data(vec![5.0, 6.0], lorentz).unwrap());
    let mixed_sum = (tensor(&[1], 1) * left_dual).n_add([tensor(&[1], 3) * right]);
    let closing = tensor(&[1], 7);
    let outside = Net::from_tensor(DenseTensor::from_storage_data(vec![11.0, 13.0], dual).unwrap());
    for original in [
        outside.clone().n_mul([closing.clone(), mixed_sum.clone()]),
        mixed_sum.n_mul([closing, outside]),
    ] {
        let mut prepared = original.clone();
        prepared.graph.contract_ready_sum_boundaries();
        for mut network in [original, prepared] {
            network
                .execute::<SequentialRef, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
                .unwrap();
            let ExecutionResult::Val(actual) = network.result_scalar().unwrap() else {
                panic!("expected the mixed-boundary scalar")
            };
            assert_eq!(*actual, 8452.0);
        }
    }
}

#[cfg(feature = "shadowing")]
#[test]
fn sparse_contraction_preserves_overlapping_and_disjoint_support() {
    use crate::{
        contraction::Contract,
        structure::{
            HasStructure, OrderedStructure, TensorStructure,
            representation::{Euclidean, RepName},
        },
        tensors::{
            data::{DataTensor, GetTensorData, SparseTensor},
            parametric::ParamTensor,
        },
    };
    use symbolica::atom::Atom;

    type Tensor = ParamTensor<OrderedStructure<Euclidean>>;
    fn tensor(axes: [(usize, usize); 2], values: &[([usize; 2], i64)]) -> Tensor {
        let ordered = OrderedStructure::new(
            axes.into_iter()
                .map(|(dim, id)| Euclidean {}.new_slot(dim, id))
                .collect(),
        );
        let elements = values
            .iter()
            .map(|&(mut indices, value)| {
                let canonical_indices = ordered.layout().logical_to_canonical(&indices);
                indices.copy_from_slice(&canonical_indices);
                (
                    ordered.canonical().flat_index(indices).unwrap(),
                    Atom::num(value),
                )
            })
            .collect();
        ParamTensor::composite(DataTensor::Sparse(SparseTensor {
            elements,
            zero: Atom::Zero,
            structure: ordered.into_canonical(),
        }))
    }

    // Three entries on either side would suggest nine Cartesian products.
    // Only k=0 overlaps: two left entries times two right entries gives four.
    let left = tensor([(2, 1), (3, 2)], &[([0, 0], 2), ([1, 0], 3), ([0, 1], 5)]);
    let right = tensor([(3, 2), (2, 3)], &[([0, 0], 7), ([0, 1], 11), ([2, 1], 13)]);
    let actual = left.contract(&right).unwrap();
    let expected = tensor(
        [(2, 1), (2, 3)],
        &[([0, 0], 14), ([0, 1], 22), ([1, 0], 21), ([1, 1], 33)],
    );
    assert_eq!(actual.structure(), expected.structure());
    for i in 0..2 {
        for j in 0..2 {
            assert_eq!(
                actual.tensor.get_ref([i, j]).unwrap(),
                expected.tensor.get_ref([i, j]).unwrap()
            );
        }
    }

    // Stored zero entries and disjoint nonzero supports both have an exactly
    // empty join. Neither an exact join nor a conservative planning estimate
    // changes the zero result or its correct open-index shape.
    for right in [
        tensor([(3, 2), (2, 3)], &[]),
        tensor([(3, 2), (2, 3)], &[([0, 0], 0)]),
        tensor([(3, 2), (2, 3)], &[([2, 0], 7), ([2, 1], 11)]),
    ] {
        let actual = left.contract(&right).unwrap();
        assert_eq!(actual.structure(), expected.structure());
        // Sparse zero entries are absent; materialize this four-entry result
        // with its stored zero before checking every coordinate.
        let actual = actual.tensor.to_bare_dense();
        for i in 0..2 {
            for j in 0..2 {
                assert_eq!(actual.get_ref([i, j]).unwrap(), &Atom::Zero);
            }
        }
    }
}
