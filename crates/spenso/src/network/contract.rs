use std::{
    fmt::{Debug, Display},
    ops::{AddAssign, MulAssign},
};

use eyre::eyre;
use linnet::half_edge::NodeIndex;

use crate::{
    algebra::ScalarMul,
    contraction::Contract,
    network::graph::{NetworkLeaf, NetworkNode, NetworkOp, NetworkOperation},
    structure::{
        HasStructure, PermutedStructure, TensorStructure,
        permuted::PermuteTensor,
        slot::{AbsInd, IsAbstractSlot},
    },
};

use super::{
    Ref, TensorNetworkError,
    graph::NetworkGraph,
    library::{Library, LibraryTensor},
    profile,
    store::NetworkStoreAccess,
};

const MAX_LAZY_TENSOR_SUM_DISTRIBUTED_TERMS: usize = 96;

pub struct SmallestDegree<CStrat = ()> {
    phantom: std::marker::PhantomData<CStrat>,
}

impl<CStrat> Default for SmallestDegree<CStrat> {
    fn default() -> Self {
        Self {
            phantom: std::marker::PhantomData,
        }
    }
}

pub struct SmallestDegreeIter<const N: usize, CStrat = ()> {
    phantom: std::marker::PhantomData<CStrat>,
}

pub struct MinResultRank<CStrat = ()> {
    phantom: std::marker::PhantomData<CStrat>,
}

impl<CStrat> Default for MinResultRank<CStrat> {
    fn default() -> Self {
        Self {
            phantom: std::marker::PhantomData,
        }
    }
}

pub struct ContractScalars<CStrat = ()> {
    phantom: std::marker::PhantomData<CStrat>,
}

impl<CStrat> Default for ContractScalars<CStrat> {
    fn default() -> Self {
        Self {
            phantom: std::marker::PhantomData,
        }
    }
}

pub struct SingleSmallestDegree<const D: bool, CStrat = ()> {
    phantom: std::marker::PhantomData<CStrat>,
}

pub struct SingleLargestDegree<const D: bool, CStrat = ()> {
    phantom: std::marker::PhantomData<CStrat>,
}

pub trait ContractionStrategy<E, L, K, FK, Aind>: Sized {
    #[allow(clippy::result_large_err)]
    fn contract(
        executor: &mut E,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display;
}

#[derive(Debug, Clone)]
pub struct ProductOperand<K, Aind> {
    leaf: NetworkLeaf<K, Aind>,
    source: Option<NodeIndex>,
}

impl<K, Aind> ProductOperand<K, Aind> {
    pub fn leaf(&self) -> &NetworkLeaf<K, Aind> {
        &self.leaf
    }

    pub fn source(&self) -> Option<NodeIndex> {
        self.source
    }
}

#[derive(Debug, Clone)]
pub struct ProductContraction<K, Aind> {
    operands: Vec<ProductOperand<K, Aind>>,
}

impl<K, Aind: AbsInd> ProductContraction<K, Aind> {
    #[allow(clippy::result_large_err)]
    pub fn from_operation<FK>(
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
    ) -> Result<Self, TensorNetworkError<K, FK>>
    where
        K: Clone + Debug + Display,
        FK: Debug + Display,
        Aind: AbsInd,
    {
        if !matches!(operation.op(), NetworkOp::Product) {
            return Err(TensorNetworkError::Other(eyre!(
                "contraction strategy received non-product operation: {}",
                operation.op()
            )));
        }

        let mut operands = Vec::with_capacity(operation.children().len());
        for child in operation.children() {
            match &graph.graph[*child] {
                NetworkNode::Leaf(leaf) => operands.push(ProductOperand {
                    leaf: leaf.clone(),
                    source: Some(*child),
                }),
                NetworkNode::Op(op) => {
                    return Err(TensorNetworkError::Other(eyre!(
                        "ready product child was still an operation: {op}"
                    )));
                }
            }
        }

        Ok(Self { operands })
    }

    pub fn operands(&self) -> &[ProductOperand<K, Aind>] {
        &self.operands
    }

    pub fn local_tensor_index(&self, operand: usize) -> Option<usize> {
        match self.operands.get(operand)?.leaf {
            NetworkLeaf::LocalTensor(index) => Some(index),
            NetworkLeaf::TensorSum(_) | NetworkLeaf::LibraryKey { .. } | NetworkLeaf::Scalar(_) => {
                None
            }
        }
    }

    pub fn tensor_structure<'a, Store>(
        &self,
        executor: &'a Store,
        operand: usize,
    ) -> Option<&'a <Store::Tensor as HasStructure>::Structure>
    where
        Store: NetworkStoreAccess,
        Store::Tensor: HasStructure,
    {
        match &self.operands.get(operand)?.leaf {
            NetworkLeaf::LocalTensor(index) => Some(executor.tensor(*index).structure()),
            NetworkLeaf::TensorSum(indices) => Some(executor.tensor(*indices.first()?).structure()),
            NetworkLeaf::LibraryKey { .. } | NetworkLeaf::Scalar(_) => None,
        }
    }

    pub fn local_tensor_structure<'a, Store>(
        &self,
        executor: &'a Store,
        operand: usize,
    ) -> Option<&'a <Store::Tensor as HasStructure>::Structure>
    where
        Store: NetworkStoreAccess,
        Store::Tensor: HasStructure,
    {
        self.tensor_structure(executor, operand)
    }

    #[allow(clippy::result_large_err)]
    pub fn materialize_libraries<LT, T, L, Sc, FK, Store>(
        &mut self,
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
        K: Clone + Debug + Display,
        FK: Debug + Display,
        Aind: AbsInd,
        LT: LibraryTensor + Clone,
        T: HasStructure + From<LT::WithIndices>,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        for operand in &mut self.operands {
            if matches!(operand.leaf, NetworkLeaf::LibraryKey { .. }) {
                let source = operand.source.ok_or_else(|| {
                    TensorNetworkError::Other(eyre!("library product operand lost source node"))
                })?;
                let tensor = graph.get_lib_data(lib, source).ok_or_else(|| {
                    TensorNetworkError::Other(eyre!(
                        "failed to materialize library product operand"
                    ))
                })?;
                let index = executor.push_tensor(T::from(tensor));
                operand.leaf = NetworkLeaf::LocalTensor(index);
                operand.source = None;
            }
        }
        Ok(())
    }

    #[allow(clippy::result_large_err)]
    pub fn contract_scalars<LT, T, L, Sc, FK, Store>(
        &mut self,
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<bool, TensorNetworkError<K, FK>>
    where
        Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
        K: Clone + Debug + Display,
        FK: Debug + Display,
        Aind: AbsInd,
        LT: LibraryTensor + Clone,
        T: HasStructure + Clone + ScalarMul<Sc, Output = T>,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        Sc: for<'a> MulAssign<Sc::Ref<'a>>
            + Clone
            + for<'a> MulAssign<T::ScalarRef<'a>>
            + From<T::Scalar>
            + Ref,
        LT::WithIndices: ScalarMul<Sc, Output = T> + PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        let mut accumulator = None;
        let mut scalar_positions = Vec::new();
        let mut tensor_positions = Vec::new();

        for (position, operand) in self.operands.iter().enumerate() {
            let is_scalar = match operand.leaf {
                NetworkLeaf::Scalar(index) => {
                    if let Some(accumulator) = &mut accumulator {
                        *accumulator *= executor.scalar(index).refer();
                    } else {
                        accumulator = Some(executor.scalar(index).clone());
                    }
                    true
                }
                NetworkLeaf::LocalTensor(index) => {
                    if let Some(scalar) = executor.tensor(index).scalar_ref() {
                        if let Some(accumulator) = &mut accumulator {
                            *accumulator *= scalar;
                        } else {
                            accumulator =
                                Some(Sc::from(executor.tensor(index).clone().scalar().unwrap()));
                        }
                        true
                    } else {
                        false
                    }
                }
                NetworkLeaf::TensorSum(_) | NetworkLeaf::LibraryKey { .. } => false,
            };

            if is_scalar {
                scalar_positions.push(position);
            } else {
                tensor_positions.push(position);
            }
        }

        let Some(accumulator) = accumulator else {
            return Ok(false);
        };

        match tensor_positions.len() {
            0 => {
                if scalar_positions.len() == 1
                    && matches!(
                        self.operands[scalar_positions[0]].leaf,
                        NetworkLeaf::Scalar(_)
                    )
                {
                    return Ok(false);
                }
                let index = executor.push_scalar(accumulator);
                self.replace_operands(
                    &scalar_positions,
                    ProductOperand {
                        leaf: NetworkLeaf::Scalar(index),
                        source: None,
                    },
                );
                Ok(true)
            }
            1 => {
                let tensor_position = tensor_positions[0];
                let leaf = self.scalar_multiply_operand::<LT, T, L, Sc, FK, Store>(
                    tensor_position,
                    accumulator,
                    executor,
                    graph,
                    lib,
                )?;
                let mut positions = scalar_positions;
                positions.push(tensor_position);
                positions.sort_unstable();
                self.replace_operands(&positions, ProductOperand { leaf, source: None });
                Ok(true)
            }
            _ if scalar_positions.len() > 1 => {
                let index = executor.push_scalar(accumulator);
                self.replace_operands(
                    &scalar_positions,
                    ProductOperand {
                        leaf: NetworkLeaf::Scalar(index),
                        source: None,
                    },
                );
                Ok(true)
            }
            _ => Ok(false),
        }
    }

    #[allow(clippy::result_large_err)]
    fn scalar_multiply_operand<LT, T, L, Sc, FK, Store>(
        &self,
        operand: usize,
        scalar: Sc,
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
        K: Clone + Debug + Display,
        FK: Debug + Display,
        Aind: AbsInd,
        LT: LibraryTensor + Clone,
        T: HasStructure + ScalarMul<Sc, Output = T>,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        LT::WithIndices: ScalarMul<Sc, Output = T> + PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        let leaf = match &self.operands[operand].leaf {
            NetworkLeaf::LocalTensor(index) => executor
                .tensor(*index)
                .scalar_mul(&scalar)
                .ok_or(TensorNetworkError::FailedScalarMul)?,
            NetworkLeaf::TensorSum(indices) => {
                let mut scaled = Vec::with_capacity(indices.len());
                for index in indices {
                    let tensor = executor
                        .tensor(*index)
                        .scalar_mul(&scalar)
                        .ok_or(TensorNetworkError::FailedScalarMul)?;
                    scaled.push(executor.push_tensor(tensor));
                }
                return Ok(NetworkLeaf::TensorSum(scaled));
            }
            NetworkLeaf::LibraryKey { .. } => {
                let source = self.operands[operand].source.ok_or_else(|| {
                    TensorNetworkError::Other(eyre!("library product operand lost source node"))
                })?;
                let tensor = graph.get_lib_data(lib, source).ok_or_else(|| {
                    TensorNetworkError::Other(eyre!(
                        "failed to materialize library product operand"
                    ))
                })?;
                tensor
                    .scalar_mul(&scalar)
                    .ok_or(TensorNetworkError::FailedScalarMul)?
            }
            NetworkLeaf::Scalar(_) => return Err(TensorNetworkError::SlotEdgeToScalarNode),
        };

        let index = executor.push_tensor(leaf);
        Ok(NetworkLeaf::LocalTensor(index))
    }

    fn tensor_term_indices(&self, operand: usize) -> Option<Vec<usize>> {
        match &self.operands.get(operand)?.leaf {
            NetworkLeaf::LocalTensor(index) => Some(vec![*index]),
            NetworkLeaf::TensorSum(indices) => Some(indices.clone()),
            NetworkLeaf::LibraryKey { .. } | NetworkLeaf::Scalar(_) => None,
        }
    }

    fn tensor_sum_leaf<T, Store>(executor: &mut Store, indices: Vec<usize>) -> NetworkLeaf<K, Aind>
    where
        Store: NetworkStoreAccess<Tensor = T>,
        T: HasStructure + Clone + Ref + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    {
        debug_assert!(!indices.is_empty());

        if indices.len() == 1 {
            return NetworkLeaf::LocalTensor(indices[0]);
        }

        if indices
            .iter()
            .all(|index| executor.tensor(*index).scalar_ref().is_some())
        {
            let mut iter = indices.into_iter();
            let first = iter
                .next()
                .expect("tensor sum with at least one term has first term");
            let mut materialized = executor.tensor(first).clone();
            for index in iter {
                materialized += executor.tensor(index).refer();
            }
            NetworkLeaf::LocalTensor(executor.push_tensor(materialized))
        } else {
            NetworkLeaf::TensorSum(indices)
        }
    }

    fn materialize_tensor_sum<T, Store>(executor: &mut Store, indices: &[usize]) -> usize
    where
        Store: NetworkStoreAccess<Tensor = T>,
        T: Clone + Ref + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    {
        debug_assert!(!indices.is_empty());

        let mut iter = indices.iter();
        let first = *iter
            .next()
            .expect("tensor sum with at least one term has first term");
        let mut materialized = executor.tensor(first).clone();
        for index in iter {
            materialized += executor.tensor(*index).refer();
        }
        executor.push_tensor(materialized)
    }

    #[allow(clippy::result_large_err)]
    pub fn contract_pair<LT, T, L, Sc, CStrat, FK, Store>(
        &mut self,
        left: usize,
        right: usize,
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
        K: Clone + Debug + Display,
        FK: Debug + Display,
        Aind: AbsInd,
        LT: LibraryTensor + Clone,
        T: HasStructure
            + From<LT::WithIndices>
            + Contract<T, CStrat, LCM = T>
            + Clone
            + Ref
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        self.materialize_libraries::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;

        let mut left_terms = self
            .tensor_term_indices(left)
            .ok_or(TensorNetworkError::SlotEdgeToScalarNode)?;
        let mut right_terms = self
            .tensor_term_indices(right)
            .ok_or(TensorNetworkError::SlotEdgeToScalarNode)?;
        let distributed_terms = left_terms.len() * right_terms.len();

        if left_terms.len() > 1
            && right_terms.len() > 1
            && distributed_terms > MAX_LAZY_TENSOR_SUM_DISTRIBUTED_TERMS
        {
            if profile::enabled() {
                eprintln!(
                    "spenso_profile product.lazy_tensor_sum_materialize left_terms={} right_terms={} distributed_terms={} max_distributed_terms={}",
                    left_terms.len(),
                    right_terms.len(),
                    distributed_terms,
                    MAX_LAZY_TENSOR_SUM_DISTRIBUTED_TERMS,
                );
            }
            if left_terms.len() > 1 {
                left_terms = vec![Self::materialize_tensor_sum(executor, &left_terms)];
            }
            if right_terms.len() > 1 {
                right_terms = vec![Self::materialize_tensor_sum(executor, &right_terms)];
            }
        }

        if profile::enabled() && (left_terms.len() > 1 || right_terms.len() > 1) {
            eprintln!(
                "spenso_profile product.lazy_tensor_sum left_terms={} right_terms={} distributed_terms={}",
                left_terms.len(),
                right_terms.len(),
                left_terms.len() * right_terms.len(),
            );
        }
        let mut contracted_indices = Vec::with_capacity(left_terms.len() * right_terms.len());
        for left_tensor in &left_terms {
            for right_tensor in &right_terms {
                let contracted = executor
                    .tensor(*left_tensor)
                    .contract(executor.tensor(*right_tensor))?;
                contracted_indices.push(executor.push_tensor(contracted));
            }
        }
        let leaf = Self::tensor_sum_leaf(executor, contracted_indices);

        let mut positions = [left, right];
        positions.sort_unstable();
        self.replace_operands(&positions, ProductOperand { leaf, source: None });
        Ok(())
    }

    #[allow(clippy::result_large_err)]
    pub fn contract_one_by_degree<
        const DEBUG: bool,
        const LARGEST: bool,
        LT,
        T,
        L,
        Sc,
        CStrat,
        FK,
        Store,
    >(
        &mut self,
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<bool, TensorNetworkError<K, FK>>
    where
        Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
        K: Clone + Debug + Display,
        FK: Debug + Display,
        Aind: AbsInd,
        LT: LibraryTensor + Clone,
        T: HasStructure
            + From<LT::WithIndices>
            + Contract<T, CStrat, LCM = T>
            + Clone
            + Ref
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
        T::Structure: Display,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        self.materialize_libraries::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        let Some((left, right, _degree)) = self.best_degree_pair::<Store, LARGEST>(executor) else {
            return Ok(false);
        };

        if DEBUG {
            println!(
                "Contracting {} with {}",
                self.tensor_structure(executor, left).unwrap(),
                self.tensor_structure(executor, right).unwrap()
            );
        }
        let profile_pair = profile::enabled();
        let log_pair = profile::verbose() && (self.operands.len() > 4 || _degree > 1);
        let pair_start = if profile_pair {
            if log_pair {
                eprintln!(
                    "spenso_profile product.pair_start operands={} left_operand={} right_operand={} degree={} left={} right={}",
                    self.operands.len(),
                    left,
                    right,
                    _degree,
                    self.tensor_structure(executor, left).unwrap(),
                    self.tensor_structure(executor, right).unwrap(),
                );
            }
            Some(std::time::Instant::now())
        } else {
            None
        };

        self.contract_pair::<LT, T, L, Sc, CStrat, FK, Store>(left, right, executor, graph, lib)?;
        if let Some(pair_start) = pair_start {
            let elapsed = pair_start.elapsed();
            if profile::verbose() || elapsed.as_millis() >= 100 {
                eprintln!(
                    "spenso_profile product.slow_pair elapsed_ms={:.3}",
                    elapsed.as_secs_f64() * 1000.0,
                );
            }
        }

        if DEBUG && let Some(structure) = self.tensor_structure(executor, left.min(right)) {
            println!("Obtained {structure}");
        }

        Ok(true)
    }

    pub fn best_degree_pair<Store, const LARGEST: bool>(
        &self,
        executor: &Store,
    ) -> Option<(usize, usize, u32)>
    where
        Store: NetworkStoreAccess,
        Store::Tensor: HasStructure,
    {
        self.best_tensor_pair_structure_by::<Store, _, LARGEST>(executor, |_, _, degree, _, _| {
            if degree == 0 { u32::MAX } else { degree }
        })
        .map(|(left, right, degree, _)| (left, right, degree))
    }

    pub fn best_result_rank_pair<Store>(&self, executor: &Store) -> Option<(usize, usize, u32)>
    where
        Store: NetworkStoreAccess,
        Store::Tensor: HasStructure,
    {
        self.best_tensor_pair_structure_by::<Store, _, false>(
            executor,
            |_, _, degree, left, right| {
                if degree == 0 {
                    return (u32::MAX, u32::MAX);
                }

                let result_rank = left.order() + right.order() - 2 * degree as usize;

                (result_rank as u32, u32::MAX - degree)
            },
        )
        .map(|(left, right, degree, _)| (left, right, degree))
    }

    #[allow(clippy::result_large_err)]
    pub fn contract_one_by_result_rank<LT, T, L, Sc, CStrat, FK, Store>(
        &mut self,
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<bool, TensorNetworkError<K, FK>>
    where
        Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
        K: Clone + Debug + Display,
        FK: Debug + Display,
        Aind: AbsInd,
        LT: LibraryTensor + Clone,
        T: HasStructure
            + From<LT::WithIndices>
            + Contract<T, CStrat, LCM = T>
            + Clone
            + Ref
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
        T::Structure: Display,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        self.materialize_libraries::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        let Some((left, right, degree)) = self.best_result_rank_pair::<Store>(executor) else {
            return Ok(false);
        };

        let profile_pair = profile::enabled();
        let log_pair = profile::verbose() && (self.operands.len() > 4 || degree > 1);
        let pair_start = if profile_pair {
            if log_pair {
                eprintln!(
                    "spenso_profile product.result_rank_pair_start operands={} left_operand={} right_operand={} degree={} left={} right={}",
                    self.operands.len(),
                    left,
                    right,
                    degree,
                    self.tensor_structure(executor, left).unwrap(),
                    self.tensor_structure(executor, right).unwrap(),
                );
            }
            Some(std::time::Instant::now())
        } else {
            None
        };

        self.contract_pair::<LT, T, L, Sc, CStrat, FK, Store>(left, right, executor, graph, lib)?;

        if let Some(pair_start) = pair_start {
            let elapsed = pair_start.elapsed();
            if profile::verbose() || elapsed.as_millis() >= 100 {
                eprintln!(
                    "spenso_profile product.result_rank_slow_pair elapsed_ms={:.3}",
                    elapsed.as_secs_f64() * 1000.0,
                );
            }
        }

        Ok(true)
    }

    pub fn best_tensor_pair_by<Store, Score: Ord, const LARGEST: bool>(
        &self,
        executor: &Store,
        mut score: impl FnMut(usize, usize, u32, &Store::Tensor, &Store::Tensor) -> Score,
    ) -> Option<(usize, usize, u32, Score)>
    where
        Store: NetworkStoreAccess,
        Store::Tensor: HasStructure,
    {
        let tensors = self
            .operands
            .iter()
            .enumerate()
            .filter_map(|(operand, _)| {
                self.local_tensor_index(operand)
                    .map(|index| (operand, index))
            })
            .collect::<Vec<_>>();

        let mut best = None;
        for (left_position, (left_operand, left_tensor)) in tensors.iter().enumerate() {
            for (right_operand, right_tensor) in &tensors[left_position + 1..] {
                let left = executor.tensor(*left_tensor);
                let right = executor.tensor(*right_tensor);
                let degree = matching_degree(left.structure(), right.structure());
                let candidate_score = score(*left_operand, *right_operand, degree, left, right);

                let replace = match &best {
                    None => true,
                    Some((_, _, _, best_score)) => {
                        if LARGEST {
                            candidate_score > *best_score
                        } else {
                            candidate_score < *best_score
                        }
                    }
                };

                if replace {
                    best = Some((*left_operand, *right_operand, degree, candidate_score));
                }
            }
        }

        best
    }

    pub fn best_tensor_pair_structure_by<Store, Score: Ord, const LARGEST: bool>(
        &self,
        executor: &Store,
        mut score: impl FnMut(
            usize,
            usize,
            u32,
            &<Store::Tensor as HasStructure>::Structure,
            &<Store::Tensor as HasStructure>::Structure,
        ) -> Score,
    ) -> Option<(usize, usize, u32, Score)>
    where
        Store: NetworkStoreAccess,
        Store::Tensor: HasStructure,
    {
        let tensors = self
            .operands
            .iter()
            .enumerate()
            .filter_map(|(operand, _)| {
                self.tensor_structure(executor, operand)
                    .map(|structure| (operand, structure))
            })
            .collect::<Vec<_>>();

        let mut best = None;
        for (left_position, (left_operand, left)) in tensors.iter().enumerate() {
            for (right_operand, right) in &tensors[left_position + 1..] {
                let degree = matching_degree(*left, *right);
                let candidate_score = score(*left_operand, *right_operand, degree, left, right);

                let replace = match &best {
                    None => true,
                    Some((_, _, _, best_score)) => {
                        if LARGEST {
                            candidate_score > *best_score
                        } else {
                            candidate_score < *best_score
                        }
                    }
                };

                if replace {
                    best = Some((*left_operand, *right_operand, degree, candidate_score));
                }
            }
        }

        best
    }

    #[allow(clippy::result_large_err)]
    pub fn finish<LT, T, L, Sc, FK, Store>(
        mut self,
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
        K: Clone + Debug + Display,
        FK: Debug + Display,
        Aind: AbsInd,
        LT: LibraryTensor + Clone,
        T: HasStructure + Clone + ScalarMul<Sc, Output = T>,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        Sc: for<'a> MulAssign<Sc::Ref<'a>>
            + Clone
            + for<'a> MulAssign<T::ScalarRef<'a>>
            + From<T::Scalar>
            + Ref,
        LT::WithIndices: ScalarMul<Sc, Output = T> + PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        while self.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)? {}

        if self.operands.len() == 1 {
            Ok(self.operands.remove(0).leaf)
        } else {
            Err(TensorNetworkError::Other(eyre!(
                "product contraction did not collapse to one leaf; {} operands remain",
                self.operands.len()
            )))
        }
    }

    fn replace_operands(&mut self, positions: &[usize], replacement: ProductOperand<K, Aind>) {
        let insert_at = positions[0];
        for position in positions.iter().rev() {
            self.operands.remove(*position);
        }
        self.operands
            .insert(insert_at.min(self.operands.len()), replacement);
    }
}

fn matching_degree<S>(left: &S, right: &S) -> u32
where
    S: TensorStructure,
{
    left.match_indices(right)
        .map(|(_, matched, _)| matched.into_iter().filter(|matched| *matched).count() as u32)
        .unwrap_or(0)
}

impl<
    CStrat,
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<T, CStrat, LCM = T>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>
        + Ref
        + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
> ContractionStrategy<Store, L, K, FK, Aind> for ContractScalars<CStrat>
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    fn contract(
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        product.finish::<LT, T, L, Sc, FK, Store>(executor, graph, lib)
    }
}

impl<
    CStrat,
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<T, CStrat, LCM = T>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>
        + Ref
        + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
> ContractionStrategy<Store, L, K, FK, Aind> for SmallestDegree<CStrat>
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <LT::WithIndices as HasStructure>::Structure: Display,
    T::Structure: Display,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    fn contract(
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        while product.contract_one_by_degree::<false, false, LT, T, L, Sc, CStrat, FK, Store>(
            executor, graph, lib,
        )? {
            product.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        }
        product.finish::<LT, T, L, Sc, FK, Store>(executor, graph, lib)
    }
}

impl<
    CStrat,
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<T, CStrat, LCM = T>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>
        + Ref
        + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
    const N: usize,
> ContractionStrategy<Store, L, K, FK, Aind> for SmallestDegreeIter<N, CStrat>
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <LT::WithIndices as HasStructure>::Structure: Display,
    T::Structure: Display,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    fn contract(
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        for _ in 0..N {
            if !product.contract_one_by_degree::<false, false, LT, T, L, Sc, CStrat, FK, Store>(
                executor, graph, lib,
            )? {
                break;
            }
            product.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        }
        product.finish::<LT, T, L, Sc, FK, Store>(executor, graph, lib)
    }
}

impl<
    CStrat,
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<T, CStrat, LCM = T>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>
        + Ref
        + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
> ContractionStrategy<Store, L, K, FK, Aind> for MinResultRank<CStrat>
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <LT::WithIndices as HasStructure>::Structure: Display,
    T::Structure: Display,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    fn contract(
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        while product
            .contract_one_by_result_rank::<LT, T, L, Sc, CStrat, FK, Store>(executor, graph, lib)?
        {
            product.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        }
        product.finish::<LT, T, L, Sc, FK, Store>(executor, graph, lib)
    }
}

impl<
    CStrat,
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<T, CStrat, LCM = T>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>
        + Ref
        + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
    const D: bool,
> ContractionStrategy<Store, L, K, FK, Aind> for SingleSmallestDegree<D, CStrat>
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <LT::WithIndices as HasStructure>::Structure: Display,
    T::Structure: Display,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    fn contract(
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        product.contract_one_by_degree::<D, false, LT, T, L, Sc, CStrat, FK, Store>(
            executor, graph, lib,
        )?;
        product.finish::<LT, T, L, Sc, FK, Store>(executor, graph, lib)
    }
}

impl<
    CStrat,
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<T, CStrat, LCM = T>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>
        + Ref
        + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
    const D: bool,
> ContractionStrategy<Store, L, K, FK, Aind> for SingleLargestDegree<D, CStrat>
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <LT::WithIndices as HasStructure>::Structure: Display,
    T::Structure: Display,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    fn contract(
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        product.contract_one_by_degree::<D, true, LT, T, L, Sc, CStrat, FK, Store>(
            executor, graph, lib,
        )?;
        product.finish::<LT, T, L, Sc, FK, Store>(executor, graph, lib)
    }
}
