use std::{fmt::Display, ops::AddAssign};

use linnet::half_edge::NodeIndex;

use crate::{
    algebra::ScalarMul,
    structure::{HasStructure, PermutedStructure, TensorStructure, permuted::PermuteTensor},
};

use super::{
    FastTensorSum, Ref, TensorNetworkError, balanced_ref_sum,
    graph::{NetworkGraph, NetworkLeaf, ScalarRef, ScaledTensorRef},
    library::{Library, LibraryTensor},
    store::NetworkStoreAccess,
};
use crate::structure::slot::{AbsInd, IsAbstractSlot};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum TensorSumOrder {
    Balanced,
    Linear,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum TensorSumFastPath {
    Enabled,
    Disabled,
}

pub(super) fn reduce_tensor_sum<T>(
    terms: Vec<T>,
    order: TensorSumOrder,
    sum_start: Option<std::time::Instant>,
) -> T
where
    T: Ref + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
{
    debug_assert!(!terms.is_empty());

    match order {
        TensorSumOrder::Balanced => balanced_ref_sum(terms, sum_start),
        TensorSumOrder::Linear => {
            let mut terms = terms.into_iter();
            let mut result = terms.next().expect("tensor sum has at least one term");
            for term in terms {
                result += term.refer();
            }
            result
        }
    }
}

pub(super) fn materialize_tensor_sum_owned<'a, K, FK, T>(
    indices: &[usize],
    order: TensorSumOrder,
    fast_path: TensorSumFastPath,
    sum_start: Option<std::time::Instant>,
    tensor: impl Fn(usize) -> &'a T,
) -> Result<T, TensorNetworkError<K, FK>>
where
    T: 'a + Clone + Ref + FastTensorSum + for<'b> AddAssign<<T as Ref>::Ref<'b>>,
    K: Display,
    FK: Display,
{
    if indices.is_empty() {
        return Err(TensorNetworkError::EmptyTensorSumLeaf);
    }

    if fast_path == TensorSumFastPath::Enabled
        && let Some(result) = {
            let terms = indices
                .iter()
                .map(|index| tensor(*index))
                .collect::<Vec<_>>();
            T::fast_tensor_sum(&terms, sum_start)
        }
    {
        return Ok(result);
    }

    if order == TensorSumOrder::Linear {
        let mut indices = indices.iter();
        let first = *indices
            .next()
            .expect("a nonempty tensor sum has a first term");
        let mut result = tensor(first).clone();
        for index in indices {
            result += tensor(*index).refer();
        }
        return Ok(result);
    }

    let terms = indices.iter().map(|index| tensor(*index).clone()).collect();
    Ok(reduce_tensor_sum(terms, order, sum_start))
}

fn materialize_tensor_sum<K, FK, T, Store>(
    store: &mut Store,
    indices: &[usize],
    order: TensorSumOrder,
    sum_start: Option<std::time::Instant>,
) -> Result<usize, TensorNetworkError<K, FK>>
where
    Store: NetworkStoreAccess<Tensor = T>,
    T: Clone + Ref + FastTensorSum + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    K: Display,
    FK: Display,
{
    let result = materialize_tensor_sum_owned::<K, FK, _>(
        indices,
        order,
        TensorSumFastPath::Enabled,
        sum_start,
        |index| store.tensor(index),
    )?;
    Ok(store.push_tensor(result))
}

pub(super) fn materialize_scaled_tensor_owned<'a, K, FK, T, Sc>(
    term: &ScaledTensorRef,
    tensor: impl Fn(usize) -> &'a T,
    scalar: impl Fn(ScalarRef) -> &'a Sc,
) -> Result<T, TensorNetworkError<K, FK>>
where
    T: 'a + Clone + ScalarMul<Sc, Output = T>,
    Sc: 'a,
    K: Display,
    FK: Display,
{
    match term.scale {
        Some(scale) => tensor(term.tensor)
            .scalar_mul(scalar(scale))
            .ok_or(TensorNetworkError::FailedScalarMul),
        None => Ok(tensor(term.tensor).clone()),
    }
}

pub(super) fn materialize_scaled_tensors_owned<'a, K, FK, T, Sc>(
    terms: &[ScaledTensorRef],
    order: TensorSumOrder,
    sum_start: Option<std::time::Instant>,
    tensor: impl Fn(usize) -> &'a T,
    scalar: impl Fn(ScalarRef) -> &'a Sc,
) -> Result<T, TensorNetworkError<K, FK>>
where
    T: 'a
        + Clone
        + Ref
        + FastTensorSum
        + ScalarMul<Sc, Output = T>
        + for<'b> AddAssign<<T as Ref>::Ref<'b>>,
    Sc: 'a,
    K: Display,
    FK: Display,
{
    if terms.is_empty() {
        return Err(TensorNetworkError::EmptyTensorSumLeaf);
    }

    if let [term] = terms {
        return materialize_scaled_tensor_owned::<K, FK, _, _>(term, tensor, scalar);
    }

    if terms.iter().all(|term| term.scale.is_none()) {
        let indices = terms.iter().map(|term| term.tensor).collect::<Vec<_>>();
        return materialize_tensor_sum_owned::<K, FK, _>(
            &indices,
            order,
            TensorSumFastPath::Enabled,
            sum_start,
            tensor,
        );
    }

    let materialized = terms
        .iter()
        .map(|term| materialize_scaled_tensor_owned::<K, FK, _, _>(term, &tensor, &scalar))
        .collect::<Result<Vec<_>, _>>()?;
    let result = {
        let terms = materialized.iter().collect::<Vec<_>>();
        T::fast_tensor_sum(&terms, sum_start)
    }
    .unwrap_or_else(|| reduce_tensor_sum(materialized, order, sum_start));
    Ok(result)
}

pub(super) fn materialize_scaled_tensors<K, FK, T, Sc, Store>(
    store: &mut Store,
    terms: &[ScaledTensorRef],
    order: TensorSumOrder,
    sum_start: Option<std::time::Instant>,
) -> Result<usize, TensorNetworkError<K, FK>>
where
    Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
    T: Clone
        + Ref
        + FastTensorSum
        + ScalarMul<Sc, Output = T>
        + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    K: Display,
    FK: Display,
{
    let result = materialize_scaled_tensors_owned::<K, FK, _, _>(
        terms,
        order,
        sum_start,
        |index| store.tensor(index),
        |scalar| store.scalar_ref(scalar),
    )?;
    Ok(store.push_tensor(result))
}

#[allow(clippy::too_many_arguments, clippy::result_large_err)]
pub(super) fn materialize_tensor_leaf<LT, L, K, FK, Aind, T, Sc, Store>(
    store: &mut Store,
    graph: &NetworkGraph<K, FK, Aind>,
    node: NodeIndex,
    leaf: &NetworkLeaf<K, Aind>,
    lib: &L,
    order: TensorSumOrder,
    sum_start: Option<std::time::Instant>,
) -> Result<usize, TensorNetworkError<K, FK>>
where
    Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
    T: HasStructure
        + Clone
        + Ref
        + FastTensorSum
        + ScalarMul<Sc, Output = T>
        + for<'a> AddAssign<<T as Ref>::Ref<'a>>
        + From<LT::WithIndices>,
    Sc: Clone,
    LT: LibraryTensor + Clone,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    K: Display + std::fmt::Debug,
    FK: Display + std::fmt::Debug,
    Aind: AbsInd,
    LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    let tensor = match leaf {
        NetworkLeaf::Scalar(_) => return Err(TensorNetworkError::SlotEdgeToScalarNode),
        NetworkLeaf::LibraryKey { .. } => {
            let tensor = T::from(graph.get_lib_data::<_, LT, L>(lib, node)?);
            store.push_tensor(tensor)
        }
        NetworkLeaf::LocalTensor(tensor) => *tensor,
        NetworkLeaf::TensorSum(indices) => {
            materialize_tensor_sum::<K, FK, _, _>(store, indices, order, sum_start)?
        }
        NetworkLeaf::ScaledTensor(term) => materialize_scaled_tensors::<K, FK, _, _, _>(
            store,
            std::slice::from_ref(term),
            order,
            sum_start,
        )?,
        NetworkLeaf::ScaledTensorSum(terms) => {
            materialize_scaled_tensors::<K, FK, _, _, _>(store, terms, order, sum_start)?
        }
    };
    Ok(tensor)
}

#[cfg(test)]
mod tests {
    use std::{
        ops::AddAssign,
        sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        },
    };

    use super::{
        TensorSumFastPath, TensorSumOrder, materialize_scaled_tensor_owned,
        materialize_scaled_tensors, materialize_scaled_tensors_owned, materialize_tensor_sum,
        materialize_tensor_sum_owned, reduce_tensor_sum,
    };
    use crate::{
        algebra::ScalarMul,
        network::{FastTensorSum, Ref, ScaledTensorRef, TensorNetworkError, store::NetworkStore},
    };

    #[derive(Clone, Copy, Debug, PartialEq)]
    struct TestFloat(f64);

    impl Ref for TestFloat {
        type Ref<'a> = &'a Self;

        fn refer(&self) -> Self::Ref<'_> {
            self
        }
    }

    impl AddAssign<&Self> for TestFloat {
        fn add_assign(&mut self, rhs: &Self) {
            self.0 += rhs.0;
        }
    }

    impl FastTensorSum for TestFloat {}

    impl ScalarMul<()> for TestFloat {
        type Output = Self;

        fn scalar_mul(&self, _rhs: &()) -> Option<Self::Output> {
            Some(*self)
        }
    }

    #[test]
    fn reduction_order_is_explicit() {
        let terms = vec![
            TestFloat(1.0e16),
            TestFloat(1.0),
            TestFloat(-1.0e16),
            TestFloat(1.0),
        ];

        assert_eq!(
            reduce_tensor_sum(terms.clone(), TensorSumOrder::Balanced, None),
            TestFloat(0.0),
        );
        assert_eq!(
            reduce_tensor_sum(terms, TensorSumOrder::Linear, None),
            TestFloat(1.0),
        );
    }

    #[cfg(feature = "shadowing")]
    #[test]
    fn symbolic_reduction_orders_preserve_the_exact_sum() {
        use symbolica::parse;

        let terms = ["a", "b", "c", "d"]
            .map(|term| parse!(term, default_namespace = "materialize_reduction"))
            .to_vec();
        let expected = parse!("a+b+c+d", default_namespace = "materialize_reduction");

        assert_eq!(
            reduce_tensor_sum(terms.clone(), TensorSumOrder::Balanced, None),
            expected,
        );
        assert_eq!(
            reduce_tensor_sum(terms, TensorSumOrder::Linear, None),
            expected,
        );
    }

    #[test]
    fn result_extraction_and_mutating_materialization_share_owned_core() {
        use crate::{
            network::{
                ExecutionResult, Network, NetworkLeaf, NetworkNode,
                library::{DummyKey, DummyLibrary, DummyLibraryTensor},
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

        let structure = OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).structure;
        let tensor = |data| DenseTensor::from_data(data, structure.clone()).unwrap();
        let store = Store {
            tensors: vec![tensor(vec![1.0, 2.0]), tensor(vec![3.0, 4.0])],
            scalar: vec![2.0],
            scalar_aliases: vec![None],
        };
        let library = Lib::new();

        let tensor_sum = materialize_tensor_sum_owned::<DummyKey, DummyKey, _>(
            &[0, 1],
            TensorSumOrder::Balanced,
            TensorSumFastPath::Enabled,
            None,
            |index| &store.tensors[index],
        )
        .unwrap();
        assert_eq!(tensor_sum.data, vec![4.0, 6.0]);
        let mut pushed = store.clone();
        let pushed_tensor_sum = materialize_tensor_sum::<DummyKey, DummyKey, _, _>(
            &mut pushed,
            &[0, 1],
            TensorSumOrder::Balanced,
            None,
        )
        .unwrap();
        assert_eq!(pushed.tensors[pushed_tensor_sum].data, tensor_sum.data);

        let mut network = Net::from_tensor(store.tensors[0].clone());
        let (_, _, node) = network.graph.graph.iter_nodes_mut().next().unwrap();
        *node = NetworkNode::Leaf(NetworkLeaf::TensorSum(vec![0, 1]));
        network.store = store.clone();
        let ExecutionResult::Val(extracted) =
            network.result_tensor::<LibTensor, Lib>(&library).unwrap()
        else {
            panic!("expected an owned TensorSum result");
        };
        assert_eq!(extracted.data, tensor_sum.data);

        let scaled_terms = [ScaledTensorRef::scaled(0, 0), ScaledTensorRef::tensor(1)];
        let scaled_sum = materialize_scaled_tensors_owned::<DummyKey, DummyKey, _, _>(
            &scaled_terms,
            TensorSumOrder::Balanced,
            None,
            |index| &store.tensors[index],
            |scalar| &store.scalar[scalar.index()],
        )
        .unwrap();
        let mut pushed = store.clone();
        let pushed_scaled_sum = materialize_scaled_tensors::<DummyKey, DummyKey, _, _, _>(
            &mut pushed,
            &scaled_terms,
            TensorSumOrder::Balanced,
            None,
        )
        .unwrap();
        assert_eq!(pushed.tensors[pushed_scaled_sum].data, scaled_sum.data);

        let mut network = Net::from_tensor(store.tensors[0].clone());
        let (_, _, node) = network.graph.graph.iter_nodes_mut().next().unwrap();
        *node = NetworkNode::Leaf(NetworkLeaf::ScaledTensorSum(scaled_terms.to_vec()));
        network.store = store;
        let ExecutionResult::Val(extracted) =
            network.result_tensor::<LibTensor, Lib>(&library).unwrap()
        else {
            panic!("expected an owned ScaledTensorSum result");
        };
        assert_eq!(extracted.data, scaled_sum.data);
        assert_eq!(scaled_sum.data, vec![5.0, 8.0]);
    }

    #[test]
    fn empty_stored_tensor_sums_return_the_typed_error() {
        let mut store = NetworkStore::<TestFloat, ()> {
            tensors: Vec::new(),
            scalar: Vec::new(),
            scalar_aliases: Vec::new(),
        };

        assert!(matches!(
            materialize_tensor_sum::<&str, &str, _, _>(
                &mut store,
                &[],
                TensorSumOrder::Balanced,
                None,
            ),
            Err(TensorNetworkError::EmptyTensorSumLeaf)
        ));
        assert!(matches!(
            materialize_scaled_tensors::<&str, &str, _, _, _>(
                &mut store,
                &[],
                TensorSumOrder::Balanced,
                None,
            ),
            Err(TensorNetworkError::EmptyTensorSumLeaf)
        ));
    }

    #[derive(Debug)]
    struct FastValue {
        value: i32,
        clones: Arc<AtomicUsize>,
        probes: Arc<AtomicUsize>,
        use_fast_sum: bool,
    }

    impl Clone for FastValue {
        fn clone(&self) -> Self {
            self.clones.fetch_add(1, Ordering::Relaxed);
            Self {
                value: self.value,
                clones: Arc::clone(&self.clones),
                probes: Arc::clone(&self.probes),
                use_fast_sum: self.use_fast_sum,
            }
        }
    }

    impl Ref for FastValue {
        type Ref<'a> = &'a Self;

        fn refer(&self) -> Self::Ref<'_> {
            self
        }
    }

    impl AddAssign<&Self> for FastValue {
        fn add_assign(&mut self, rhs: &Self) {
            self.value += rhs.value;
        }
    }

    impl FastTensorSum for FastValue {
        fn fast_tensor_sum(
            terms: &[&Self],
            _sum_start: Option<std::time::Instant>,
        ) -> Option<Self> {
            let first = terms.first()?;
            first.probes.fetch_add(1, Ordering::Relaxed);
            if !first.use_fast_sum {
                return None;
            }
            Some(Self {
                value: terms.iter().map(|term| term.value).sum(),
                clones: Arc::clone(&first.clones),
                probes: Arc::clone(&first.probes),
                use_fast_sum: true,
            })
        }
    }

    impl ScalarMul<()> for FastValue {
        type Output = Self;

        fn scalar_mul(&self, _rhs: &()) -> Option<Self::Output> {
            Some(self.clone())
        }
    }

    #[test]
    fn unscaled_materialization_uses_borrowed_fast_sum() {
        let clones = Arc::new(AtomicUsize::new(0));
        let probes = Arc::new(AtomicUsize::new(0));
        let mut store = NetworkStore::<FastValue, ()> {
            tensors: vec![
                FastValue {
                    value: 1,
                    clones: Arc::clone(&clones),
                    probes: Arc::clone(&probes),
                    use_fast_sum: true,
                },
                FastValue {
                    value: 2,
                    clones: Arc::clone(&clones),
                    probes: Arc::clone(&probes),
                    use_fast_sum: true,
                },
            ],
            scalar: Vec::new(),
            scalar_aliases: Vec::new(),
        };

        let result = materialize_scaled_tensors::<&str, &str, _, _, _>(
            &mut store,
            &[ScaledTensorRef::tensor(0), ScaledTensorRef::tensor(1)],
            TensorSumOrder::Balanced,
            None,
        )
        .unwrap();

        assert_eq!(store.tensors[result].value, 3);
        assert_eq!(clones.load(Ordering::Relaxed), 0);
        assert_eq!(probes.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn unscaled_materialization_probes_fast_sum_once_before_fallback() {
        let clones = Arc::new(AtomicUsize::new(0));
        let probes = Arc::new(AtomicUsize::new(0));
        let mut store = NetworkStore::<FastValue, ()> {
            tensors: vec![
                FastValue {
                    value: 1,
                    clones: Arc::clone(&clones),
                    probes: Arc::clone(&probes),
                    use_fast_sum: false,
                },
                FastValue {
                    value: 2,
                    clones: Arc::clone(&clones),
                    probes: Arc::clone(&probes),
                    use_fast_sum: false,
                },
            ],
            scalar: Vec::new(),
            scalar_aliases: Vec::new(),
        };

        let result = materialize_scaled_tensors::<&str, &str, _, _, _>(
            &mut store,
            &[ScaledTensorRef::tensor(0), ScaledTensorRef::tensor(1)],
            TensorSumOrder::Balanced,
            None,
        )
        .unwrap();

        assert_eq!(store.tensors[result].value, 3);
        assert_eq!(clones.load(Ordering::Relaxed), 2);
        assert_eq!(probes.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn one_scaled_tensor_has_the_same_boundary_as_direct_materialization() {
        let clones = Arc::new(AtomicUsize::new(0));
        let probes = Arc::new(AtomicUsize::new(0));
        let tensor = FastValue {
            value: 7,
            clones: Arc::clone(&clones),
            probes: Arc::clone(&probes),
            use_fast_sum: true,
        };
        let term = ScaledTensorRef::tensor(0);

        let direct =
            materialize_scaled_tensor_owned::<&str, &str, _, _>(&term, |_| &tensor, |_| &())
                .unwrap();
        let through_sum = materialize_scaled_tensors_owned::<&str, &str, _, _>(
            std::slice::from_ref(&term),
            TensorSumOrder::Balanced,
            None,
            |_| &tensor,
            |_| &(),
        )
        .unwrap();

        assert_eq!(through_sum.value, direct.value);
        assert_eq!(clones.load(Ordering::Relaxed), 2);
        assert_eq!(probes.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn linear_unscaled_fallback_clones_only_the_accumulator() {
        let clones = Arc::new(AtomicUsize::new(0));
        let probes = Arc::new(AtomicUsize::new(0));
        let tensors = [1, 2, 3].map(|value| FastValue {
            value,
            clones: Arc::clone(&clones),
            probes: Arc::clone(&probes),
            use_fast_sum: false,
        });

        let result = materialize_tensor_sum_owned::<&str, &str, _>(
            &[0, 1, 2],
            TensorSumOrder::Linear,
            TensorSumFastPath::Enabled,
            None,
            |index| &tensors[index],
        )
        .unwrap();

        assert_eq!(result.value, 6);
        assert_eq!(clones.load(Ordering::Relaxed), 1);
        assert_eq!(probes.load(Ordering::Relaxed), 1);

        clones.store(0, Ordering::Relaxed);
        probes.store(0, Ordering::Relaxed);
        let strict = materialize_tensor_sum_owned::<&str, &str, _>(
            &[0, 1, 2],
            TensorSumOrder::Linear,
            TensorSumFastPath::Disabled,
            None,
            |index| &tensors[index],
        )
        .unwrap();
        assert_eq!(strict.value, 6);
        assert_eq!(clones.load(Ordering::Relaxed), 1);
        assert_eq!(probes.load(Ordering::Relaxed), 0);
    }
}
