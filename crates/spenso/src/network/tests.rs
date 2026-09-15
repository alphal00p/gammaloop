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
    use symbolica::{atom::Atom, symbol};

    use super::{
        Network, NetworkLeaf, NetworkNode,
        store::{NetworkStore, TensorScalarStore},
        tags::scalar_store_alias,
    };

    let original = Atom::var(symbol!("x"));
    let mut net: Network<NetworkStore<(), Atom>, i8, i8> = Network::from_scalar(original.clone());
    let aliases = net.alias_scalar_refs(|_, _| true);

    assert_eq!(aliases.aliases_created(), 1);
    let (node, _, _) = net.graph.result().unwrap();
    let NetworkNode::Leaf(NetworkLeaf::Scalar(scalar)) = node else {
        panic!("expected scalar result node");
    };

    assert_eq!(net.store.get_scalar_ref(*scalar), &scalar_store_alias(0));
    assert_eq!(
        net.resolve_scalar_aliases(&aliases, scalar_store_alias(0)),
        original
    );
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
            assert_eq!(sum.store.tensors.len(), tensor.store.tensors.len());
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
