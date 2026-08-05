use crate::network::library::DummyKey;

use super::TensorNetworkError;

#[test]
fn display() {
    let a = TensorNetworkError::<i8, DummyKey>::Infallible;

    println!("{a}")
}

#[test]
fn leaf_reference_mapping_covers_every_store_backed_variant() {
    use linnet::permutation::Permutation;

    use crate::structure::PermutedStructure;

    use super::{NetworkLeaf, ScaledTensorRef};
    use crate::network::graph::ScalarRef;

    let mut leaves: Vec<NetworkLeaf<DummyKey>> = vec![
        NetworkLeaf::LocalTensor(0),
        NetworkLeaf::TensorSum(vec![1, 2]),
        NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled_ref(3, ScalarRef::Store(4))),
        NetworkLeaf::ScaledTensorSum(vec![
            ScaledTensorRef::tensor(5),
            ScaledTensorRef::scaled_ref(6, ScalarRef::Alias(7)),
        ]),
        NetworkLeaf::LibraryKey {
            key: PermutedStructure {
                structure: DummyKey::default(),
                rep_permutation: Permutation::id(0),
                index_permutation: Permutation::id(0),
            },
            indices: vec![],
        },
        NetworkLeaf::Scalar(ScalarRef::Store(8)),
    ];

    for leaf in &mut leaves {
        leaf.map_tensor_refs(|tensor| tensor + 10);
        leaf.map_scalar_refs(|scalar| scalar.map_index(|index| index + 20));
    }

    assert_eq!(leaves[0], NetworkLeaf::LocalTensor(10));
    assert_eq!(leaves[1], NetworkLeaf::TensorSum(vec![11, 12]));
    assert_eq!(
        leaves[2],
        NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled_ref(13, ScalarRef::Store(24)))
    );
    assert_eq!(
        leaves[3],
        NetworkLeaf::ScaledTensorSum(vec![
            ScaledTensorRef::tensor(15),
            ScaledTensorRef::scaled_ref(16, ScalarRef::Alias(27)),
        ])
    );
    assert!(matches!(leaves[4], NetworkLeaf::LibraryKey { .. }));
    assert_eq!(leaves[5], NetworkLeaf::Scalar(ScalarRef::Store(28)));
}

#[test]
fn binary_and_nary_assembly_rebase_store_references_once() {
    use linnet::permutation::Permutation;

    use super::{NAdd, NMul, Network, NetworkLeaf, NetworkNode, NetworkState, ScaledTensorRef};
    use crate::network::graph::ScalarRef;
    use crate::network::{graph::NetworkGraph, store::NetworkStore};
    use crate::structure::PermutedStructure;

    type Store = NetworkStore<i32, i32>;
    type Net = Network<Store, DummyKey, DummyKey>;

    fn network(leaf: NetworkLeaf<DummyKey>, tensors: Vec<i32>, scalars: Vec<i32>) -> Net {
        let mut graph = NetworkGraph::scalar(0);
        let (_, _, node) = graph.graph.iter_nodes_mut().next().unwrap();
        *node = NetworkNode::Leaf(leaf);
        let scalar_aliases = vec![None; scalars.len()];
        Network {
            graph,
            store: NetworkStore {
                tensors,
                scalar: scalars,
                scalar_aliases,
            },
            state: NetworkState::PureScalar,
        }
    }

    fn assert_leaf(net: &Net, expected: &NetworkLeaf<DummyKey>) {
        assert!(
            net.graph
                .graph
                .iter_nodes()
                .any(|(_, _, node)| matches!(node, NetworkNode::Leaf(leaf) if leaf == expected)),
            "missing rebased leaf {expected:?}",
        );
    }

    let base = network(NetworkLeaf::Scalar(0.into()), vec![10], vec![100]);
    let variants = [
        (NetworkLeaf::LocalTensor(0), NetworkLeaf::LocalTensor(1)),
        (
            NetworkLeaf::TensorSum(vec![0, 1]),
            NetworkLeaf::TensorSum(vec![1, 2]),
        ),
        (
            NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(0, 0)),
            NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(1, 1)),
        ),
        (
            NetworkLeaf::ScaledTensorSum(vec![
                ScaledTensorRef::tensor(0),
                ScaledTensorRef::scaled_ref(1, ScalarRef::Alias(1)),
            ]),
            NetworkLeaf::ScaledTensorSum(vec![
                ScaledTensorRef::tensor(1),
                ScaledTensorRef::scaled_ref(2, ScalarRef::Alias(2)),
            ]),
        ),
        (
            NetworkLeaf::Scalar(ScalarRef::Store(0)),
            NetworkLeaf::Scalar(ScalarRef::Store(1)),
        ),
        (
            NetworkLeaf::LibraryKey {
                key: PermutedStructure {
                    structure: DummyKey::default(),
                    rep_permutation: Permutation::id(0),
                    index_permutation: Permutation::id(0),
                },
                indices: vec![],
            },
            NetworkLeaf::LibraryKey {
                key: PermutedStructure {
                    structure: DummyKey::default(),
                    rep_permutation: Permutation::id(0),
                    index_permutation: Permutation::id(0),
                },
                indices: vec![],
            },
        ),
    ];

    for (source, expected) in variants {
        let source = network(source, vec![20, 21], vec![200, 201]);
        for combined in [
            base.clone() + source.clone(),
            base.clone() * source.clone(),
            base.clone().n_add([source.clone()]),
            base.clone().n_mul([source.clone()]),
        ] {
            assert_eq!(combined.store.tensors, vec![10, 20, 21]);
            assert_eq!(combined.store.scalar, vec![100, 200, 201]);
            assert_leaf(&combined, &expected);
        }
    }

    let second = network(
        NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(0, 0)),
        vec![20],
        vec![200],
    );
    let third = network(
        NetworkLeaf::ScaledTensorSum(vec![
            ScaledTensorRef::tensor(0),
            ScaledTensorRef::scaled(1, 0),
        ]),
        vec![30, 31],
        vec![300],
    );
    let expected_second = NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(1, 1));
    let expected_third = NetworkLeaf::ScaledTensorSum(vec![
        ScaledTensorRef::tensor(2),
        ScaledTensorRef::scaled(3, 2),
    ]);

    for combined in [
        base.clone().n_add([second.clone(), third.clone()]),
        base.n_mul([second, third]),
    ] {
        assert_eq!(combined.store.tensors, vec![10, 20, 30, 31]);
        assert_eq!(combined.store.scalar, vec![100, 200, 300]);
        assert_leaf(&combined, &expected_second);
        assert_leaf(&combined, &expected_third);
    }
}

#[test]
fn negation_preserves_lazy_tensor_leaf_variants() {
    use linnet::permutation::Permutation;

    use super::{Network, NetworkLeaf, NetworkNode, NetworkState, ScaledTensorRef};
    use crate::{
        network::{graph::NetworkGraph, store::NetworkStore},
        structure::PermutedStructure,
    };

    type Net = Network<NetworkStore<i32, i32>, DummyKey, DummyKey>;

    let variants = [
        NetworkLeaf::LibraryKey {
            key: PermutedStructure {
                structure: DummyKey::default(),
                rep_permutation: Permutation::id(0),
                index_permutation: Permutation::id(0),
            },
            indices: vec![],
        },
        NetworkLeaf::LocalTensor(0),
        NetworkLeaf::TensorSum(vec![0, 1]),
        NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(0, 0)),
        NetworkLeaf::ScaledTensorSum(vec![
            ScaledTensorRef::tensor(0),
            ScaledTensorRef::scaled(1, 0),
        ]),
    ];

    for leaf in variants {
        let mut graph = NetworkGraph::scalar(0);
        let (_, _, node) = graph.graph.iter_nodes_mut().next().unwrap();
        *node = NetworkNode::Leaf(leaf.clone());
        let negated = -Net {
            graph,
            store: NetworkStore {
                tensors: vec![1, 2],
                scalar: vec![3],
                scalar_aliases: vec![None],
            },
            state: NetworkState::SelfDualTensor,
        };

        assert!(negated.graph.graph.iter_nodes().any(|(_, _, node)| {
            matches!(node, NetworkNode::Leaf(candidate) if candidate == &leaf)
        }));
    }
}

#[test]
fn executed_negation_preserves_lazy_tensor_leaf_variants() {
    use crate::{
        network::{
            ExecutionResult, Network, NetworkLeaf, NetworkNode, ScaledTensorRef, Sequential,
            SmallestDegree,
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
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;
    type FnLib = ErroringLibrary<DummyKey>;

    fn leaf_network(leaf: NetworkLeaf<DummyKey>, tensors: Vec<Tensor>, scalars: Vec<f64>) -> Net {
        let mut network = Net::from_tensor(tensors[0].clone());
        let (_, _, node) = network.graph.graph.iter_nodes_mut().next().unwrap();
        *node = NetworkNode::Leaf(leaf);
        network.store = NetworkStore {
            tensors,
            scalar_aliases: vec![None; scalars.len()],
            scalar: scalars,
        };
        network
    }

    let structure = OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).structure;
    let tensor = |data| DenseTensor::from_data(data, structure.clone()).unwrap();
    let variants = [
        (
            "TensorSum",
            leaf_network(
                NetworkLeaf::TensorSum(vec![0, 1]),
                vec![tensor(vec![1.0, 2.0]), tensor(vec![3.0, 4.0])],
                vec![],
            ),
        ),
        (
            "ScaledTensor",
            leaf_network(
                NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(0, 0)),
                vec![tensor(vec![2.0, 3.0])],
                vec![2.0],
            ),
        ),
        (
            "ScaledTensorSum",
            leaf_network(
                NetworkLeaf::ScaledTensorSum(vec![
                    ScaledTensorRef::scaled(0, 0),
                    ScaledTensorRef::tensor(1),
                ]),
                vec![tensor(vec![0.5, 1.0]), tensor(vec![3.0, 4.0])],
                vec![2.0],
            ),
        ),
    ];
    let library = Lib::new();
    let functions = FnLib::new();

    for (kind, network) in variants {
        let mut negated = -network;
        negated
            .execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(&library, &functions)
            .unwrap();
        let (node, _, _) = negated.graph.result().unwrap();
        assert!(
            matches!(
                (kind, node),
                ("TensorSum", NetworkNode::Leaf(NetworkLeaf::TensorSum(_)))
                    | (
                        "ScaledTensor",
                        NetworkNode::Leaf(NetworkLeaf::ScaledTensor(_))
                    )
                    | (
                        "ScaledTensorSum",
                        NetworkNode::Leaf(NetworkLeaf::ScaledTensorSum(_))
                    )
            ),
            "executed Neg changed the {kind} representation to {node:?}",
        );
        let ExecutionResult::Val(result) =
            negated.result_tensor::<LibTensor, Lib>(&library).unwrap()
        else {
            panic!("expected a tensor result for {kind}");
        };
        assert_eq!(result.data, vec![-4.0, -6.0], "wrong Neg result for {kind}");
    }
}

#[test]
fn sum_execution_preserves_configured_lazy_eager_policy() {
    const TEST_MODE: &str = "SPENSO_TEST_LAZY_TENSOR_SUM_MODE";
    const LAZY_MODE: &str = "lazy";
    const CHILD_OK: &str = "spenso-sum-policy-child-ok";

    let Ok(mode) = std::env::var(TEST_MODE) else {
        let executable = std::env::current_exe().unwrap();
        let test_name = std::thread::current().name().unwrap().to_owned();
        for mode in ["eager", LAZY_MODE] {
            let mut command = std::process::Command::new(&executable);
            command
                .arg("--exact")
                .arg(&test_name)
                .arg("--nocapture")
                .env(TEST_MODE, mode);
            if mode == LAZY_MODE {
                command.env("SPENSO_NETWORK_LAZY_TENSOR_SUMS", "1");
            } else {
                command.env_remove("SPENSO_NETWORK_LAZY_TENSOR_SUMS");
            }
            let output = command.output().unwrap();
            let stdout = String::from_utf8_lossy(&output.stdout);
            assert!(
                output.status.success() && stdout.contains(&format!("{CHILD_OK}:{mode}")),
                "{mode} Sum child failed:\nstdout:\n{}\nstderr:\n{}",
                stdout,
                String::from_utf8_lossy(&output.stderr),
            );
        }
        return;
    };

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
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;
    type FnLib = ErroringLibrary<DummyKey>;

    let structure = OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).structure;
    let left = DenseTensor::from_data(vec![1.0, 2.0], structure.clone()).unwrap();
    let right = DenseTensor::from_data(vec![3.0, 4.0], structure).unwrap();
    let mut sum = Net::from_tensor(left) + Net::from_tensor(right);
    let library = Lib::new();
    sum.execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(&library, &FnLib::new())
        .unwrap();

    let (node, _, _) = sum.graph.result().unwrap();
    if mode == LAZY_MODE {
        assert!(
            matches!(node, NetworkNode::Leaf(NetworkLeaf::TensorSum(_))),
            "lazy Sum produced {node:?}",
        );
    } else {
        assert!(
            matches!(node, NetworkNode::Leaf(NetworkLeaf::LocalTensor(_))),
            "eager Sum produced {node:?}",
        );
    }
    let ExecutionResult::Val(result) = sum.result_tensor::<LibTensor, Lib>(&library).unwrap()
    else {
        panic!("expected a tensor Sum result");
    };
    assert_eq!(result.data, vec![4.0, 6.0]);
    println!("{CHILD_OK}:{mode}");
}

#[test]
fn add_i8_rebases_the_appended_scalar_reference() {
    use crate::{
        network::{Network, NetworkLeaf, NetworkNode, store::NetworkStore},
        structure::OrderedStructure,
        tensors::data::DenseTensor,
    };

    type Tensor = DenseTensor<i8, OrderedStructure>;
    type Net = Network<NetworkStore<Tensor, i8>, DummyKey, DummyKey>;

    let combined = Net::from_scalar(11) + 7;
    let mut scalar_indices = combined
        .graph
        .graph
        .iter_nodes()
        .filter_map(|(_, _, node)| match node {
            NetworkNode::Leaf(NetworkLeaf::Scalar(scalar)) => Some(scalar.index()),
            _ => None,
        })
        .collect::<Vec<_>>();
    scalar_indices.sort_unstable();

    assert_eq!(scalar_indices, vec![0, 1]);
    assert_eq!(combined.store.scalar, vec![11, 7]);
}

#[test]
fn parallel_replacement_rebases_only_overlay_references() {
    use super::{NetworkLeaf, ScaledTensorRef, remap_parallel_replacement};
    use crate::network::graph::ScalarRef;

    let mut replacement: NetworkLeaf<DummyKey> = NetworkLeaf::ScaledTensorSum(vec![
        ScaledTensorRef::scaled_ref(4, ScalarRef::Alias(6)),
        ScaledTensorRef::scaled_ref(5, ScalarRef::Store(7)),
        ScaledTensorRef::scaled_ref(7, ScalarRef::Alias(9)),
    ]);

    remap_parallel_replacement(&mut replacement, 5, 7, 11, 13);

    assert_eq!(
        replacement,
        NetworkLeaf::ScaledTensorSum(vec![
            ScaledTensorRef::scaled_ref(4, ScalarRef::Alias(6)),
            ScaledTensorRef::scaled_ref(11, ScalarRef::Store(13)),
            ScaledTensorRef::scaled_ref(13, ScalarRef::Alias(15)),
        ])
    );
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

#[test]
fn appended_populated_scalar_alias_keeps_its_value_and_variant() {
    use super::{
        Network, NetworkLeaf, NetworkNode,
        graph::ScalarRef,
        store::{NetworkStore, TensorScalarStore},
    };

    let mut right: Network<NetworkStore<(), i32>, i8, i8> = Network::from_scalar(2);
    right.store.scalar_aliases[0] = Some(7);
    right
        .graph
        .map_scalar_refs(|scalar| scalar.with_alias_state(|_| true));

    let combined = Network::from_scalar(11) * right;
    assert!(combined.graph.graph.iter_nodes().any(|(_, _, node)| {
        matches!(
            node,
            NetworkNode::Leaf(NetworkLeaf::Scalar(ScalarRef::Alias(1)))
        )
    }));
    assert_eq!(combined.store.get_scalar_ref(ScalarRef::Alias(1)), &7,);
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
        OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).structure;
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

    let structure = OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).structure;
    let a = DenseTensor::from_data(vec![1.0, 2.0], structure.clone()).unwrap();
    let b = DenseTensor::from_data(vec![3.0, 4.0], structure).unwrap();
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
fn mixed_scalar_sum_is_independent_of_leaf_variant_and_child_order() {
    use itertools::Itertools;

    use crate::{
        network::{
            ExecutionResult, NAdd, Network, NetworkLeaf, NetworkNode, NetworkState,
            ScaledTensorRef, Sequential, SmallestDegree,
            graph::NetworkGraph,
            library::{DummyLibrary, DummyLibraryTensor, panicing::ErroringLibrary},
            store::NetworkStore,
        },
        structure::{OrderedStructure, ScalarTensor, representation::Euclidean},
        tensors::data::DenseTensor,
    };

    type Tensor = DenseTensor<f64, OrderedStructure<Euclidean>>;
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;
    type FnLib = ErroringLibrary<DummyKey>;

    fn leaf_network(leaf: NetworkLeaf<DummyKey>, tensors: Vec<Tensor>, scalars: Vec<f64>) -> Net {
        let mut graph = NetworkGraph::scalar(0);
        let (_, _, node) = graph.graph.iter_nodes_mut().next().unwrap();
        *node = NetworkNode::Leaf(leaf);
        let scalar_aliases = vec![None; scalars.len()];
        Network {
            graph,
            store: NetworkStore {
                tensors,
                scalar: scalars,
                scalar_aliases,
            },
            state: NetworkState::PureScalar,
        }
    }

    fn execute_scalar(
        mut net: Net,
        lib: &Lib,
        fn_lib: &FnLib,
    ) -> Result<f64, TensorNetworkError<DummyKey, DummyKey>> {
        net.execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(lib, fn_lib)?;
        Ok(match net.result_scalar()? {
            ExecutionResult::Val(value) => *value,
            ExecutionResult::Zero => 0.0,
            ExecutionResult::One => 1.0,
        })
    }

    fn sum(children: &[Net], order: &[usize]) -> Net {
        let mut ordered = order.iter().map(|index| children[*index].clone());
        ordered
            .next()
            .expect("a test permutation has at least one child")
            .n_add(ordered)
    }

    let scalar_tensor = |value| Tensor::new_scalar(value);
    let children = vec![
        Net::from_scalar(5.0),
        leaf_network(
            NetworkLeaf::LocalTensor(0),
            vec![scalar_tensor(3.0)],
            vec![],
        ),
        leaf_network(
            NetworkLeaf::TensorSum(vec![0, 1]),
            vec![scalar_tensor(3.0), scalar_tensor(4.0)],
            vec![],
        ),
        leaf_network(
            NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(0, 0)),
            vec![scalar_tensor(3.0)],
            vec![2.0],
        ),
        leaf_network(
            NetworkLeaf::ScaledTensorSum(vec![
                ScaledTensorRef::scaled(0, 0),
                ScaledTensorRef::tensor(1),
            ]),
            vec![scalar_tensor(3.0), scalar_tensor(4.0)],
            vec![2.0],
        ),
    ];
    let lib = Lib::new();
    let fn_lib = FnLib::new();

    for order in (0..children.len()).permutations(children.len()) {
        assert_eq!(
            execute_scalar(sum(&children, &order), &lib, &fn_lib).unwrap(),
            31.0,
            "mixed scalar Sum changed for child order {order:?}",
        );
    }

    for invalid in [
        NetworkLeaf::TensorSum(Vec::new()),
        NetworkLeaf::ScaledTensorSum(Vec::new()),
    ] {
        let invalid_children = vec![
            Net::from_scalar(5.0),
            leaf_network(
                NetworkLeaf::LocalTensor(0),
                vec![scalar_tensor(3.0)],
                vec![],
            ),
            leaf_network(invalid, vec![], vec![]),
        ];
        for order in (0..invalid_children.len()).permutations(invalid_children.len()) {
            assert!(
                matches!(
                    execute_scalar(sum(&invalid_children, &order), &lib, &fn_lib),
                    Err(TensorNetworkError::EmptyTensorSumLeaf)
                ),
                "mixed scalar Sum error changed for child order {order:?}",
            );
        }
    }
}

#[test]
fn empty_stored_tensor_sum_leaves_are_rejected_at_use_boundaries() {
    use crate::{
        network::{
            Network, NetworkLeaf, NetworkNode, NetworkState, Sequential, SmallestDegree,
            graph::NetworkGraph,
            library::{DummyLibrary, DummyLibraryTensor, panicing::ErroringLibrary},
            store::NetworkStore,
        },
        structure::{OrderedStructure, representation::Euclidean},
        tensors::data::DenseTensor,
    };

    type Tensor = DenseTensor<f64, OrderedStructure<Euclidean>>;
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;
    type FnLib = ErroringLibrary<DummyKey>;

    fn leaf_network(leaf: NetworkLeaf<DummyKey>) -> Net {
        let mut graph = NetworkGraph::scalar(0);
        let (_, _, node) = graph.graph.iter_nodes_mut().next().unwrap();
        *node = NetworkNode::Leaf(leaf);
        Network {
            graph,
            store: NetworkStore::default(),
            state: NetworkState::PureScalar,
        }
    }

    let lib = Lib::new();
    let fn_lib = FnLib::new();
    for leaf in [
        NetworkLeaf::TensorSum(vec![]),
        NetworkLeaf::ScaledTensorSum(vec![]),
    ] {
        assert!(!super::sum_target_is_scalar(
            &NetworkStore::<Tensor, f64>::default(),
            &leaf,
        ));

        let net = leaf_network(leaf);
        assert!(matches!(
            net.result_scalar(),
            Err(TensorNetworkError::EmptyTensorSumLeaf)
        ));
        assert!(matches!(
            net.result_tensor::<LibTensor, Lib>(&lib),
            Err(TensorNetworkError::EmptyTensorSumLeaf)
        ));

        let mut wrapped = net.fun(DummyKey::default());
        assert!(matches!(
            wrapped.execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib),
            Err(TensorNetworkError::EmptyTensorSumLeaf)
        ));
    }
}

#[test]
fn function_materializes_every_stored_tensor_leaf_variant() {
    use crate::{
        network::{
            ExecutionResult, Network, NetworkLeaf, NetworkNode, NetworkState, ScaledTensorRef,
            Sequential, SmallestDegree,
            graph::NetworkGraph,
            library::{DummyLibrary, DummyLibraryTensor, FunctionLibrary, FunctionLibraryError},
            store::NetworkStore,
        },
        structure::{
            OrderedStructure,
            representation::{Euclidean, RepName},
        },
        tensors::data::DenseTensor,
    };

    type Tensor = DenseTensor<f64, OrderedStructure<Euclidean>>;
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;

    struct IdentityFunctions;

    impl FunctionLibrary<Tensor, f64> for IdentityFunctions {
        type Key = DummyKey;

        fn apply(
            &self,
            _key: &Self::Key,
            tensor: Tensor,
        ) -> Result<Tensor, FunctionLibraryError<Self::Key>> {
            Ok(tensor)
        }

        fn apply_scalar(
            &self,
            _key: &Self::Key,
            scalar: f64,
        ) -> Result<f64, FunctionLibraryError<Self::Key>> {
            Ok(scalar)
        }
    }

    fn function_network(
        leaf: NetworkLeaf<DummyKey>,
        tensors: Vec<Tensor>,
        scalars: Vec<f64>,
    ) -> Net {
        let mut graph = NetworkGraph::scalar(0);
        let (_, _, node) = graph.graph.iter_nodes_mut().next().unwrap();
        *node = NetworkNode::Leaf(leaf);
        let graph = graph.function(DummyKey::default());
        let scalar_aliases = vec![None; scalars.len()];
        Network {
            graph,
            store: NetworkStore {
                tensors,
                scalar: scalars,
                scalar_aliases,
            },
            state: NetworkState::SelfDualTensor,
        }
    }

    let structure = OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).structure;
    let tensor = |data| DenseTensor::from_data(data, structure.clone()).unwrap();
    let variants = [
        (
            function_network(
                NetworkLeaf::LocalTensor(0),
                vec![tensor(vec![1.0, 2.0])],
                vec![],
            ),
            vec![1.0, 2.0],
        ),
        (
            function_network(
                NetworkLeaf::TensorSum(vec![0, 1]),
                vec![tensor(vec![1.0, 2.0]), tensor(vec![3.0, 4.0])],
                vec![],
            ),
            vec![4.0, 6.0],
        ),
        (
            function_network(
                NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(0, 0)),
                vec![tensor(vec![1.0, 2.0])],
                vec![2.0],
            ),
            vec![2.0, 4.0],
        ),
        (
            function_network(
                NetworkLeaf::ScaledTensorSum(vec![
                    ScaledTensorRef::scaled(0, 0),
                    ScaledTensorRef::tensor(1),
                ]),
                vec![tensor(vec![1.0, 2.0]), tensor(vec![3.0, 4.0])],
                vec![2.0],
            ),
            vec![5.0, 8.0],
        ),
    ];
    let lib = Lib::new();
    let fn_lib = IdentityFunctions;

    for (mut net, expected) in variants {
        net.execute::<Sequential, SmallestDegree, LibTensor, Lib, IdentityFunctions>(&lib, &fn_lib)
            .unwrap();
        let ExecutionResult::Val(result) = net.result_tensor::<LibTensor, Lib>(&lib).unwrap()
        else {
            panic!("expected materialized function tensor result");
        };
        assert_eq!(result.data, expected);
    }
}

#[test]
fn self_loop_trace_materializes_stored_tensor_leaf_variants() {
    use crate::{
        network::{
            ExecutionResult, Network, NetworkLeaf, NetworkNode, ScaledTensorRef, Sequential,
            SmallestDegree,
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
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;
    type FnLib = ErroringLibrary<DummyKey>;

    fn trace_network(leaf: NetworkLeaf<DummyKey>, tensors: Vec<Tensor>, scalars: Vec<f64>) -> Net {
        let mut net = Net::from_tensor(tensors[0].clone());
        let (_, _, node) = net.graph.graph.iter_nodes_mut().next().unwrap();
        *node = NetworkNode::Leaf(leaf);
        let scalar_aliases = vec![None; scalars.len()];
        net.store = NetworkStore {
            tensors,
            scalar: scalars,
            scalar_aliases,
        };
        net
    }

    let slot = Euclidean {}.new_slot(2, 1);
    let structure = OrderedStructure::new(vec![slot, slot]).structure;
    let tensor = |diagonal: (f64, f64)| {
        DenseTensor::from_data(vec![diagonal.0, 0.0, 0.0, diagonal.1], structure.clone()).unwrap()
    };
    let variants = [
        (
            trace_network(
                NetworkLeaf::LocalTensor(0),
                vec![tensor((1.0, 4.0))],
                vec![],
            ),
            5.0,
        ),
        (
            trace_network(
                NetworkLeaf::TensorSum(vec![0, 1]),
                vec![tensor((1.0, 4.0)), tensor((10.0, 40.0))],
                vec![],
            ),
            55.0,
        ),
        (
            trace_network(
                NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(0, 0)),
                vec![tensor((1.0, 4.0))],
                vec![2.0],
            ),
            10.0,
        ),
        (
            trace_network(
                NetworkLeaf::ScaledTensorSum(vec![
                    ScaledTensorRef::scaled(0, 0),
                    ScaledTensorRef::tensor(1),
                ]),
                vec![tensor((1.0, 4.0)), tensor((10.0, 40.0))],
                vec![2.0],
            ),
            60.0,
        ),
    ];
    let lib = Lib::new();
    let fn_lib = FnLib::new();

    for (mut net, expected) in variants {
        net.execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
            .unwrap();
        let ExecutionResult::Val(result) = net.result_scalar().unwrap() else {
            panic!("expected traced scalar result");
        };
        assert_eq!(*result, expected);
    }
}

#[test]
fn zero_and_identity_powers_are_normalized_before_execution() {
    use crate::network::store::NetworkStore;
    use crate::network::{Network, NetworkNode, NetworkOp, NetworkState};

    type Net = Network<NetworkStore<(), f64>, DummyKey, DummyKey>;

    assert_eq!(NetworkState::Tensor.pow(0), NetworkState::PureScalar);

    let original = Net::from_scalar(3.0);
    assert_eq!(original.clone().pow(1), original);

    let mut stale_tensor_state = original;
    stale_tensor_state.state = NetworkState::Tensor;
    let zero = stale_tensor_state.pow(0);
    assert_eq!(zero.state, NetworkState::PureScalar);
    assert!(zero.store.scalar.is_empty());
    let (result, _, _) = zero.graph.result().unwrap();
    assert!(matches!(result, NetworkNode::Op(NetworkOp::Product)));
    assert!(
        !zero
            .graph
            .graph
            .iter_nodes()
            .any(|(_, _, node)| { matches!(node, NetworkNode::Op(NetworkOp::Power(_))) })
    );
}

#[test]
fn power_error_variants_have_matching_messages() {
    assert_eq!(
        TensorNetworkError::<DummyKey, DummyKey>::InvalidDotFunction(": bad".into()).to_string(),
        "Invalid dot function: bad",
    );
    assert_eq!(
        TensorNetworkError::<DummyKey, DummyKey>::NonSelfDualTensorPower(": bad".into())
            .to_string(),
        "Non self-dual tensor power: bad",
    );
}

#[test]
fn scalar_and_tensor_powers_use_the_exact_unsigned_magnitude() {
    use crate::{
        network::{
            ExecutionResult, Network, NetworkState, Sequential, SmallestDegree,
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
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
    type LibTensor = DummyLibraryTensor<Tensor>;
    type Lib = DummyLibrary<Tensor, DummyKey>;
    type FnLib = ErroringLibrary<DummyKey>;

    fn execute(net: &mut Net, lib: &Lib, fn_lib: &FnLib) {
        net.execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(lib, fn_lib)
            .unwrap();
    }

    fn scalar_result(net: &Net) -> f64 {
        match net.result_scalar().unwrap() {
            ExecutionResult::One => 1.0,
            ExecutionResult::Zero => 0.0,
            ExecutionResult::Val(value) => *value,
        }
    }

    let lib = Lib::new();
    let fn_lib = FnLib::new();
    for exponent in -5_i8..=5 {
        let mut power = Net::from_scalar(2.0).pow(exponent);
        execute(&mut power, &lib, &fn_lib);
        assert_eq!(scalar_result(&power), 2.0_f64.powi(i32::from(exponent)));
    }

    let mut minimum_power = Net::from_scalar(2.0).pow(i8::MIN);
    execute(&mut minimum_power, &lib, &fn_lib);
    assert_eq!(scalar_result(&minimum_power), 2.0_f64.powi(-128));

    let structure = OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).structure;
    let vector = DenseTensor::from_data(vec![1.0, 2.0], structure).unwrap();

    let binary_closed = Net::from_tensor(vector.clone()) * Net::from_tensor(vector.clone());
    assert_eq!(binary_closed.state, NetworkState::Scalar);
    let mut negative_odd_closed = binary_closed.pow(-3);
    assert_eq!(negative_odd_closed.state, NetworkState::Scalar);
    execute(&mut negative_odd_closed, &lib, &fn_lib);
    assert_eq!(scalar_result(&negative_odd_closed), 5.0_f64.powi(-3));

    let convenience_closed = Net::from_tensor(vector.clone()) * vector.clone();
    assert_eq!(convenience_closed.state, NetworkState::Scalar);

    let identity = Net::from_tensor(vector.clone());
    assert_eq!(identity.clone().pow(1), identity);

    let zero = Net::from_tensor(vector.clone()).pow(0);
    assert!(zero.store.scalar.is_empty());
    assert!(zero.store.tensors.is_empty());
    assert!(matches!(zero.result().unwrap(), ExecutionResult::One));

    for (exponent, expected) in [(3, vec![5.0, 10.0]), (5, vec![25.0, 50.0])] {
        let mut power = Net::from_tensor(vector.clone()).pow(exponent);
        execute(&mut power, &lib, &fn_lib);
        let ExecutionResult::Val(result) = power.result_tensor::<LibTensor, Lib>(&lib).unwrap()
        else {
            panic!("expected a tensor result for exponent {exponent}");
        };
        assert_eq!(result.data, expected);
    }

    for (exponent, expected) in [(2, 5.0), (4, 25.0), (-2, 0.2), (-4, 0.04)] {
        let mut power = Net::from_tensor(vector.clone()).pow(exponent);
        execute(&mut power, &lib, &fn_lib);
        assert_eq!(scalar_result(&power), expected);
    }

    let mut negative_odd = Net::from_tensor(vector).pow(-3);
    assert!(matches!(
        negative_odd.execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib),
        Err(TensorNetworkError::NegativeExponentNonScalar(_)),
    ));
}

#[test]
fn equivalent_tensor_leaf_variants_share_eager_execution_semantics() {
    use std::{
        borrow::Cow,
        sync::atomic::{AtomicUsize, Ordering},
    };

    use linnet::permutation::Permutation;

    use crate::{
        network::{
            ExecutionResult, Network, NetworkLeaf, ScaledTensorRef, Sequential, SmallestDegree,
            graph::NetworkGraph,
            library::{FunctionLibrary, FunctionLibraryError, Library, LibraryError},
            store::NetworkStore,
        },
        structure::{
            OrderedStructure, PermutedStructure, ScalarTensor, TensorStructure,
            representation::{Euclidean, RepName},
        },
        tensors::data::{DataTensor, DenseTensor, HasTensorData},
    };

    type Structure = OrderedStructure<Euclidean>;
    type Tensor = DataTensor<f64, Structure>;
    type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;

    struct ValueLibrary {
        value: PermutedStructure<Tensor>,
        get_calls: AtomicUsize,
        fail_get: bool,
    }

    impl Library<Structure> for ValueLibrary {
        type Key = DummyKey;
        type Value = PermutedStructure<Tensor>;

        fn key_for_structure(
            &self,
            _structure: &PermutedStructure<Structure>,
        ) -> Result<Self::Key, LibraryError<Self::Key>> {
            Ok(DummyKey::default())
        }

        fn get<'a>(
            &'a self,
            key: &Self::Key,
        ) -> Result<Cow<'a, Self::Value>, LibraryError<Self::Key>> {
            self.get_calls.fetch_add(1, Ordering::Relaxed);
            if self.fail_get {
                Err(LibraryError::NotFound(*key))
            } else {
                Ok(Cow::Borrowed(&self.value))
            }
        }
    }

    struct IdentityFunctions;

    impl FunctionLibrary<Tensor, f64> for IdentityFunctions {
        type Key = DummyKey;

        fn apply(
            &self,
            _key: &Self::Key,
            tensor: Tensor,
        ) -> Result<Tensor, FunctionLibraryError<Self::Key>> {
            Ok(tensor)
        }

        fn apply_scalar(
            &self,
            _key: &Self::Key,
            scalar: f64,
        ) -> Result<f64, FunctionLibraryError<Self::Key>> {
            Ok(scalar)
        }
    }

    struct EquivalentLeaves {
        library: ValueLibrary,
        variants: Vec<(&'static str, Net)>,
    }

    impl EquivalentLeaves {
        fn new(value: Tensor, half: Tensor) -> Self {
            fn leaf_network(
                value: &Tensor,
                leaf: NetworkLeaf<DummyKey>,
                tensors: Vec<Tensor>,
                scalars: Vec<f64>,
            ) -> Net {
                let graph = NetworkGraph::tensor(value, leaf);
                let state = graph.state();
                let scalar_aliases = vec![None; scalars.len()];
                Network {
                    graph,
                    store: NetworkStore {
                        tensors,
                        scalar: scalars,
                        scalar_aliases,
                    },
                    state,
                }
            }

            let order = value.order();
            let indices = value.external_indices_iter().collect();
            let key = PermutedStructure {
                structure: DummyKey::default(),
                rep_permutation: Permutation::id(order),
                index_permutation: Permutation::id(order),
            };
            let variants = vec![
                (
                    "LibraryKey",
                    leaf_network(
                        &value,
                        NetworkLeaf::LibraryKey { key, indices },
                        vec![],
                        vec![],
                    ),
                ),
                (
                    "LocalTensor",
                    leaf_network(
                        &value,
                        NetworkLeaf::LocalTensor(0),
                        vec![value.clone()],
                        vec![],
                    ),
                ),
                (
                    "TensorSum",
                    leaf_network(
                        &value,
                        NetworkLeaf::TensorSum(vec![0, 1]),
                        vec![half.clone(), half.clone()],
                        vec![],
                    ),
                ),
                (
                    "ScaledTensor",
                    leaf_network(
                        &value,
                        NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled(0, 0)),
                        vec![half.clone()],
                        vec![2.0],
                    ),
                ),
                (
                    "ScaledTensorSum",
                    leaf_network(
                        &value,
                        NetworkLeaf::ScaledTensorSum(vec![
                            ScaledTensorRef::scaled(0, 0),
                            ScaledTensorRef::tensor(1),
                        ]),
                        vec![value.clone(), half],
                        vec![0.5],
                    ),
                ),
            ];

            Self {
                library: ValueLibrary {
                    value: PermutedStructure::identity(value),
                    get_calls: AtomicUsize::new(0),
                    fail_get: false,
                },
                variants,
            }
        }
    }

    fn dense(data: Vec<f64>, structure: Structure) -> Tensor {
        DenseTensor::from_data(data, structure).unwrap().into()
    }

    fn execute(net: &mut Net, library: &ValueLibrary) {
        net.execute::<Sequential, SmallestDegree, Tensor, ValueLibrary, IdentityFunctions>(
            library,
            &IdentityFunctions,
        )
        .unwrap();
    }

    fn tensor_result(net: &Net, library: &ValueLibrary) -> Vec<f64> {
        let ExecutionResult::Val(result) =
            net.result_tensor::<Tensor, ValueLibrary>(library).unwrap()
        else {
            panic!("expected tensor result");
        };
        result.as_ref().data()
    }

    fn scalar_result(net: &Net) -> f64 {
        match net.result_scalar().unwrap() {
            ExecutionResult::Val(value) => *value,
            ExecutionResult::Zero => 0.0,
            ExecutionResult::One => 1.0,
        }
    }

    let slot = Euclidean {}.new_slot(2, 1);
    let vector_structure = OrderedStructure::new(vec![slot]).structure;
    let vector = dense(vec![1.0, 0.0], vector_structure.clone());
    let half_vector = dense(vec![0.5, 0.0], vector_structure);

    let identity_fixture = EquivalentLeaves::new(vector.clone(), half_vector.clone());
    for (kind, base) in &identity_fixture.variants {
        identity_fixture
            .library
            .get_calls
            .store(0, Ordering::Relaxed);
        let mut zero = base.clone().pow(0);
        execute(&mut zero, &identity_fixture.library);
        assert!(
            matches!(zero.result().unwrap(), ExecutionResult::One),
            "{kind} Power(0) did not become one before materialization",
        );
        assert!(zero.store.tensors.is_empty());
        assert!(zero.store.scalar.is_empty());
        assert_eq!(
            identity_fixture.library.get_calls.load(Ordering::Relaxed),
            0,
            "{kind} Power(0) read the tensor library",
        );

        let mut identity = base.clone().pow(1);
        execute(&mut identity, &identity_fixture.library);
        assert_eq!(
            identity, *base,
            "{kind} Power(1) changed the leaf before materialization",
        );
        assert_eq!(
            identity_fixture.library.get_calls.load(Ordering::Relaxed),
            0,
            "{kind} Power(1) read the tensor library",
        );

        let mut function = base.clone().fun(DummyKey::default());
        execute(&mut function, &identity_fixture.library);
        assert_eq!(
            tensor_result(&function, &identity_fixture.library),
            vec![1.0, 0.0],
            "{kind} Function result differs",
        );
    }

    for exponent in [1, 2, 3, 4, 5, -2, -4, i8::MIN] {
        for (kind, base) in &identity_fixture.variants {
            let mut power = (-base.clone()).pow(exponent);
            execute(&mut power, &identity_fixture.library);
            if exponent > 0 && exponent % 2 != 0 {
                assert_eq!(
                    tensor_result(&power, &identity_fixture.library),
                    vec![-1.0, 0.0],
                    "{kind} Power({exponent}) differs",
                );
            } else {
                assert_eq!(
                    scalar_result(&power),
                    1.0,
                    "{kind} Power({exponent}) differs",
                );
            }
        }
    }

    for exponent in [-5, -3, -1] {
        for (kind, base) in &identity_fixture.variants {
            let mut power = (-base.clone()).pow(exponent);
            assert!(
                matches!(
                    power.execute::<
                        Sequential,
                        SmallestDegree,
                        Tensor,
                        ValueLibrary,
                        IdentityFunctions,
                    >(&identity_fixture.library, &IdentityFunctions),
                    Err(TensorNetworkError::NegativeExponentNonScalar(_)),
                ),
                "{kind} Power({exponent}) did not return the shared typed error",
            );
        }
    }

    let scalar_fixture = EquivalentLeaves::new(Tensor::new_scalar(2.0), Tensor::new_scalar(1.0));
    for exponent in (-5_i8..=5).chain([i8::MIN]) {
        for (kind, base) in &scalar_fixture.variants {
            let mut power = (-base.clone()).pow(exponent);
            execute(&mut power, &scalar_fixture.library);
            assert_eq!(
                scalar_result(&power),
                (-2.0_f64).powi(i32::from(exponent)),
                "{kind} scalar Power({exponent}) differs",
            );
        }
    }

    let trace_structure = OrderedStructure::new(vec![slot, slot]).structure;
    let trace_fixture = EquivalentLeaves::new(
        dense(vec![1.0, 0.0, 0.0, 3.0], trace_structure.clone()),
        dense(vec![0.5, 0.0, 0.0, 1.5], trace_structure),
    );

    let failing_library = ValueLibrary {
        value: identity_fixture.library.value.clone(),
        get_calls: AtomicUsize::new(0),
        fail_get: true,
    };
    let assert_lookup_failure = |kind: &str, mut network: Net| {
        failing_library.get_calls.store(0, Ordering::Relaxed);
        let error = network
            .execute::<Sequential, SmallestDegree, Tensor, ValueLibrary, IdentityFunctions>(
                &failing_library,
                &IdentityFunctions,
            )
            .unwrap_err();
        assert!(
            matches!(
                &error,
                TensorNetworkError::LibErr(LibraryError::NotFound(_))
            ),
            "{kind} returned {error}",
        );
        assert_eq!(
            failing_library.get_calls.load(Ordering::Relaxed),
            1,
            "{kind} retried a failed tensor-library lookup",
        );
    };
    assert_lookup_failure(
        "Function",
        identity_fixture.variants[0]
            .1
            .clone()
            .fun(DummyKey::default()),
    );
    assert_lookup_failure("Power", identity_fixture.variants[0].1.clone().pow(2));
    assert_lookup_failure("self-loop trace", trace_fixture.variants[0].1.clone());

    for (kind, mut trace) in trace_fixture.variants.clone() {
        execute(&mut trace, &trace_fixture.library);
        assert_eq!(scalar_result(&trace), 4.0, "{kind} self-loop trace differs",);
    }
}
