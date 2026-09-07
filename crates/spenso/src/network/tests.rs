use crate::network::library::DummyKey;

use super::TensorNetworkError;

#[test]
fn display() {
    let a = TensorNetworkError::<i8, DummyKey>::Infallible;

    println!("{a}")
}

#[cfg(feature = "shadowing")]
#[test]
fn scalar_alias_refs_resolve_to_original_atom() {
    use symbolica::{
        atom::{AliasedAtom, Atom, AtomCore},
        function, parse, symbol,
    };

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

    // Network results may mix opaque references with newly formed scalar bodies.
    // Definitions can share a value or refer to earlier aliases without requiring
    // a second common-subexpression extraction over that result.
    let original = parse!("x+y");
    let mut net: Network<NetworkStore<(), Atom>, i8, i8> = Network::from_scalar(original.clone());
    assert_eq!(net.store.add_scalar(scalar_store_alias(0) + parse!("z")), 1);
    assert_eq!(net.store.add_scalar(original.clone()), 2);
    assert_eq!(net.store.add_scalar(parse!("unused")), 3);
    let aliases = net.alias_scalar_refs(|_, _| true);
    assert_eq!(aliases.aliases_created(), 4);

    let root = scalar_store_alias(1).pow(2)
        + function!(symbol!("f"), scalar_store_alias(0))
        + scalar_store_alias(2).pow(-1)
        + original.pow(3)
        + scalar_store_alias(99);
    let expected = parse!("(x+y+z)^2+f(x+y)+1/(x+y)+(x+y)^3") + scalar_store_alias(99);
    let aliased = net.aliased_atom(&aliases, root.clone());
    assert_eq!(aliased.get_root(), &root);
    assert_eq!(aliased.get_aliases().len(), 4);
    for index in aliases.aliased_indices() {
        assert_eq!(
            aliased.get_aliases()[&scalar_store_alias(index)],
            net.store.scalar[index]
        );
    }
    assert_eq!(aliased.clone().into_inner(), expected);
    assert_eq!(net.resolve_scalar_aliases(&aliases, root.clone()), expected);

    // Compare the former recompression path on this small exact fixture, while
    // requiring registration to preserve the already constructed root verbatim.
    let mut recompressed = AliasedAtom::from(root);
    for index in aliases.aliased_indices() {
        recompressed =
            recompressed.add_alias(scalar_store_alias(index), net.store.scalar[index].clone());
    }
    assert_eq!(aliased.get_aliases(), recompressed.get_aliases());
    assert_eq!(recompressed.into_inner(), expected);
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

#[cfg(feature = "shadowing")]
#[test]
fn large_scaled_tensor_sum_contracts_without_broadcasting_coefficients() {
    use std::collections::HashMap;

    use symbolica::{
        atom::{Atom, AtomCore},
        function, symbol,
    };

    use super::AtomSumShapeDiagnostics;
    use crate::{
        network::{
            ExecutionResult, ExecutionStrategy, MinIntermediateCost, MinResultRank, Network,
            NetworkLeaf, NetworkNode, Sequential, SequentialRef, TensorOrScalarOrKey,
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
        let structure = OrderedStructure::new(vec![Euclidean {}.new_slot(width, 1)]).structure;
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
        if width == 2 {
            let mut small = sum.clone();
            // Inspect the internal strategy before terminal finalization: tiny
            // sums should still be eager during intermediate execution.
            <Sequential as ExecutionStrategy<Store, FnLib, Lib, DummyKey, DummyKey, _>>::execute_all::<$method>(
                &mut small.store, &mut small.graph, &lib, &fn_lib,
            ).unwrap();
            assert!(matches!(
                small.graph.result().unwrap().0,
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(_))
            ));
        }

        // The independent reference contracts integer vectors first and only
        // then multiplies their dot products by the symbolic coefficients.
        let dot_a: i64 = a.iter().zip(&weights).map(|(a, w)| a * w).sum();
        let dot_b: i64 = b.iter().zip(&weights).map(|(b, w)| b * w).sum();
        let expected = &coefficient * Atom::num(dot_a) + &other_coefficient * Atom::num(dot_b);
        let mut contracted = sum.clone() * Net::from_tensor(tensor(&weights));
        contracted
            .execute::<Sequential, $method, LibTensor, Lib, FnLib>(&lib, &fn_lib)
            .unwrap();
        let ExecutionResult::Val(actual) = contracted.result_scalar().unwrap() else {
            panic!("expected a scalar contraction");
        };
        assert!((actual.as_ref() - &expected).expand().is_zero());
        assert!(
            contracted.store.tensors.iter().all(|tensor| {
                tensor.atom_sum_shape_stats(true).total_bytes < super::MAX_EAGER_TENSOR_SUM_BYTES
            }),
            "contraction should not broadcast the large coefficients into an intermediate tensor"
        );
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
                        structure.rep_permutation.apply_slice_in_place(&mut indices);
                        structure
                            .index_permutation
                            .apply_slice_in_place(&mut indices);
                        (
                            structure.structure.flat_index(indices).unwrap(),
                            Atom::num(weight * (j as i64 + 1)),
                        )
                    })
                })
                .collect::<HashMap<_, _>>(),
            zero: Atom::Zero,
            structure: matrix_structure.structure,
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
fn intermediate_cost_avoids_a_lower_rank_larger_numeric_intermediate() {
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
            ordered.rep_permutation.apply_slice_in_place(&mut indices);
            ordered.index_permutation.apply_slice_in_place(&mut indices);
            let flat: usize = ordered.structure.flat_index(indices).unwrap().into();
            data[flat] = entry as f64;
        }
        DenseTensor::from_data(data, ordered.structure).unwrap()
    }

    let a = |i: usize, j: usize, x: usize| (1 + i + 2 * j + x % 3) as i64;
    let b = |i: usize, j: usize, y: usize| (2 * i + j + y % 5) as i64 - 3;
    let c = |x: usize, z: usize| (x % 7 + z) as i64 - 2;
    let original = Net::from_tensor(tensor(&[(2, 1), (3, 2), (16, 3)], |v| a(v[0], v[1], v[2])))
        * Net::from_tensor(tensor(&[(2, 1), (3, 2), (16, 4)], |v| b(v[0], v[1], v[2])))
        * Net::from_tensor(tensor(&[(16, 3), (1, 5)], |v| c(v[0], v[1])));
    let input_tensors = original.store.tensors.len();

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

    let largest_generated = |net: &Net| {
        net.store
            .tensors
            .iter()
            .skip(input_tensors)
            .map(|tensor| tensor.data.len())
            .max()
            .unwrap()
    };
    // AB has rank two but 16*16 entries. AC has rank three with only
    // 2*3*1 entries; its final contraction has 16 entries. The append-only
    // store retains intermediates, so this inspects actual generated tensors.
    assert_eq!(largest_generated(&rank_first), 256);
    assert_eq!(largest_generated(&intermediate_cost), 16);
}

#[cfg(feature = "shadowing")]
#[test]
fn sparse_pair_cost_counts_only_matching_contracted_coordinates() {
    use crate::{
        contraction::Contract,
        network::FastTensorSumContractible,
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
                ordered.rep_permutation.apply_slice_in_place(&mut indices);
                ordered.index_permutation.apply_slice_in_place(&mut indices);
                (
                    ordered.structure.flat_index(indices).unwrap(),
                    Atom::num(value),
                )
            })
            .collect();
        ParamTensor::composite(DataTensor::Sparse(SparseTensor {
            elements,
            zero: Atom::Zero,
            structure: ordered.structure,
        }))
    }

    // Three entries on either side would suggest nine Cartesian products.
    // Only k=0 overlaps: two left entries times two right entries gives four.
    let left = tensor([(2, 1), (3, 2)], &[([0, 0], 2), ([1, 0], 3), ([0, 1], 5)]);
    let right = tensor([(3, 2), (2, 3)], &[([0, 0], 7), ([0, 1], 11), ([2, 1], 13)]);
    let (permutation, left_matches, right_matches) =
        left.structure().match_indices(right.structure()).unwrap();
    for exact_join_limit in [0, 20_000] {
        let estimate = left.contraction_pair_estimate(
            &right,
            &permutation,
            &left_matches,
            &right_matches,
            left.contraction_profile(),
            right.contraction_profile(),
            4,
            exact_join_limit,
        );
        assert_eq!(estimate.estimated_products, 4);
        assert_eq!(estimate.estimated_output_entries, 4);
    }
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
    // empty join. The planner must not replace that known zero work by the
    // profile fallback; contraction must retain the correct open-index shape.
    for right in [
        tensor([(3, 2), (2, 3)], &[]),
        tensor([(3, 2), (2, 3)], &[([0, 0], 0)]),
        tensor([(3, 2), (2, 3)], &[([2, 0], 7), ([2, 1], 11)]),
    ] {
        let (permutation, left_matches, right_matches) =
            left.structure().match_indices(right.structure()).unwrap();
        for exact_join_limit in [0, 20_000] {
            let estimate = left.contraction_pair_estimate(
                &right,
                &permutation,
                &left_matches,
                &right_matches,
                left.contraction_profile(),
                right.contraction_profile(),
                4,
                exact_join_limit,
            );
            assert_eq!(estimate.estimated_products, 0);
        }
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
