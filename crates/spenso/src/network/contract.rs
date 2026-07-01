use std::{
    env,
    fmt::{Debug, Display},
    ops::{AddAssign, MulAssign},
    sync::OnceLock,
};

use eyre::eyre;
use linnet::half_edge::{
    NodeIndex,
    subgraph::{SuBitGraph, SubSetLike},
};

use crate::{
    algebra::ScalarMul,
    contraction::Contract,
    network::graph::{NetworkLeaf, NetworkNode, NetworkOp, NetworkOperation, TensorTerm},
    structure::{
        HasStructure, PermutedStructure, StructureContract, TensorStructure,
        permuted::PermuteTensor,
        slot::{AbsInd, IsAbstractSlot},
    },
};

use super::{
    AtomComponentOptimizable, FastTensorSum, FastTensorSumContract, FastTensorSumContractible, Ref,
    TensorCommonFactor, TensorContractionProfile, TensorNetworkError,
    graph::NetworkGraph,
    library::{Library, LibraryTensor},
    profile,
    store::NetworkStoreAccess,
};

const MAX_LAZY_TENSOR_SUM_DISTRIBUTED_TERMS: usize = 96;

fn max_lazy_tensor_sum_distributed_terms() -> usize {
    static VALUE: OnceLock<usize> = OnceLock::new();
    *VALUE.get_or_init(|| {
        env::var("SPENSO_NETWORK_MAX_LAZY_DISTRIBUTED_TERMS")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(MAX_LAZY_TENSOR_SUM_DISTRIBUTED_TERMS)
    })
}

const PAIR_SCORE_ORDER_LEN: usize = 12;
const PAIR_SCORE_ORDER_END: u8 = u8::MAX;

pub const PAIR_SCORE_RESULT_RANK: u8 = 0;
pub const PAIR_SCORE_ESTIMATED_OUTPUT_ENTRIES: u8 = 1;
pub const PAIR_SCORE_OUTPUT_DENSE_SIZE: u8 = 2;
pub const PAIR_SCORE_MAX_OUTPUT_ENTRY_PRODUCTS: u8 = 3;
pub const PAIR_SCORE_SIMPLE_TENSOR_PENALTY: u8 = 4;
pub const PAIR_SCORE_COMMON_FACTOR_PENALTY: u8 = 5;
pub const PAIR_SCORE_ATOM_WORK: u8 = 6;
pub const PAIR_SCORE_MAX_ENTRY_GROWTH: u8 = 7;
pub const PAIR_SCORE_ENTRY_WORK: u8 = 8;
pub const PAIR_SCORE_INVERSE_DEGREE: u8 = 9;

pub const DEFAULT_EXACT_JOIN_LIMIT: usize = 20_000;

pub const fn pack_pair_score_order(components: [u8; PAIR_SCORE_ORDER_LEN]) -> u128 {
    let mut packed = 0u128;
    let mut index = 0;
    while index < PAIR_SCORE_ORDER_LEN {
        packed |= (components[index] as u128) << (index * 8);
        index += 1;
    }
    packed
}

pub const fn pair_score_order_contains(score_order: u128, component: u8) -> bool {
    let mut index = 0;
    while index < PAIR_SCORE_ORDER_LEN {
        let current = ((score_order >> (index * 8)) & 0xff) as u8;
        if current == PAIR_SCORE_ORDER_END {
            return false;
        }
        if current == component {
            return true;
        }
        index += 1;
    }
    false
}

const fn pair_score_order_needs_sparse_estimate(score_order: u128) -> bool {
    pair_score_order_contains(score_order, PAIR_SCORE_ESTIMATED_OUTPUT_ENTRIES)
        || pair_score_order_contains(score_order, PAIR_SCORE_MAX_OUTPUT_ENTRY_PRODUCTS)
}

const fn pair_score_order_needs_output_dense_size(score_order: u128) -> bool {
    pair_score_order_needs_sparse_estimate(score_order)
        || pair_score_order_contains(score_order, PAIR_SCORE_OUTPUT_DENSE_SIZE)
}

pub const PAIR_SCORE_SPARSE_ATOM_AWARE: u128 = pack_pair_score_order([
    PAIR_SCORE_RESULT_RANK,
    PAIR_SCORE_ESTIMATED_OUTPUT_ENTRIES,
    PAIR_SCORE_OUTPUT_DENSE_SIZE,
    PAIR_SCORE_MAX_OUTPUT_ENTRY_PRODUCTS,
    PAIR_SCORE_SIMPLE_TENSOR_PENALTY,
    PAIR_SCORE_COMMON_FACTOR_PENALTY,
    PAIR_SCORE_ATOM_WORK,
    PAIR_SCORE_MAX_ENTRY_GROWTH,
    PAIR_SCORE_ENTRY_WORK,
    PAIR_SCORE_INVERSE_DEGREE,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
]);

pub const PAIR_SCORE_ATOM_AWARE: u128 = pack_pair_score_order([
    PAIR_SCORE_RESULT_RANK,
    PAIR_SCORE_ATOM_WORK,
    PAIR_SCORE_MAX_ENTRY_GROWTH,
    PAIR_SCORE_ENTRY_WORK,
    PAIR_SCORE_INVERSE_DEGREE,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
]);

pub const PAIR_SCORE_RESULT_RANK_ONLY: u128 = pack_pair_score_order([
    PAIR_SCORE_RESULT_RANK,
    PAIR_SCORE_INVERSE_DEGREE,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
]);

pub const PAIR_SCORE_ENTRY_AWARE: u128 = pack_pair_score_order([
    PAIR_SCORE_RESULT_RANK,
    PAIR_SCORE_ENTRY_WORK,
    PAIR_SCORE_ATOM_WORK,
    PAIR_SCORE_INVERSE_DEGREE,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
    PAIR_SCORE_ORDER_END,
]);

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ResultRankPairScore {
    result_rank: u128,
    estimated_output_entries: u128,
    output_dense_size: u128,
    max_output_entry_products: u128,
    simple_tensor_penalty: u128,
    common_factor_penalty: u128,
    atom_work: u128,
    max_entry_growth: u128,
    entry_work: u128,
    inverse_degree: u128,
}

impl ResultRankPairScore {
    fn new(
        result_rank: u32,
        degree: u32,
        left: TensorContractionProfile,
        right: TensorContractionProfile,
        pair: super::TensorContractionPairEstimate,
    ) -> Self {
        let left_entries = left.entries.max(1) as u128;
        let right_entries = right.entries.max(1) as u128;
        let left_terms = left.total_terms.max(1) as u128;
        let right_terms = right.total_terms.max(1) as u128;
        let left_bytes = left.total_bytes.max(1) as u128;
        let right_bytes = right.total_bytes.max(1) as u128;
        let left_max_terms = left.max_terms.max(1) as u128;
        let right_max_terms = right.max_terms.max(1) as u128;
        let left_max_bytes = left.max_bytes.max(1) as u128;
        let right_max_bytes = right.max_bytes.max(1) as u128;

        Self {
            result_rank: result_rank as u128,
            estimated_output_entries: pair.estimated_output_entries,
            output_dense_size: pair.output_dense_size,
            max_output_entry_products: pair.max_output_entry_products,
            simple_tensor_penalty: pair.simple_tensor_penalty,
            common_factor_penalty: pair.common_factor_penalty,
            atom_work: left_bytes
                .saturating_mul(right_terms)
                .saturating_add(right_bytes.saturating_mul(left_terms)),
            max_entry_growth: left_max_bytes
                .saturating_mul(right_max_terms)
                .saturating_add(right_max_bytes.saturating_mul(left_max_terms)),
            entry_work: left_entries.saturating_mul(right_entries),
            inverse_degree: (u32::MAX - degree) as u128,
        }
    }

    fn less_by<const SCORE_ORDER: u128>(&self, other: &Self) -> bool {
        let mut index = 0;
        while index < PAIR_SCORE_ORDER_LEN {
            let component = ((SCORE_ORDER >> (index * 8)) & 0xff) as u8;
            if component == PAIR_SCORE_ORDER_END {
                return false;
            }
            let left = self.component(component);
            let right = other.component(component);
            if left != right {
                return left < right;
            }
            index += 1;
        }
        false
    }

    fn component(&self, component: u8) -> u128 {
        match component {
            PAIR_SCORE_RESULT_RANK => self.result_rank,
            PAIR_SCORE_ESTIMATED_OUTPUT_ENTRIES => self.estimated_output_entries,
            PAIR_SCORE_OUTPUT_DENSE_SIZE => self.output_dense_size,
            PAIR_SCORE_MAX_OUTPUT_ENTRY_PRODUCTS => self.max_output_entry_products,
            PAIR_SCORE_SIMPLE_TENSOR_PENALTY => self.simple_tensor_penalty,
            PAIR_SCORE_COMMON_FACTOR_PENALTY => self.common_factor_penalty,
            PAIR_SCORE_ATOM_WORK => self.atom_work,
            PAIR_SCORE_MAX_ENTRY_GROWTH => self.max_entry_growth,
            PAIR_SCORE_ENTRY_WORK => self.entry_work,
            PAIR_SCORE_INVERSE_DEGREE => self.inverse_degree,
            _ => u128::MAX,
        }
    }
}

pub struct SmallestDegree<CStrat = (), COpt = ()> {
    phantom: std::marker::PhantomData<(CStrat, COpt)>,
}

impl<CStrat, COpt> Default for SmallestDegree<CStrat, COpt> {
    fn default() -> Self {
        Self {
            phantom: std::marker::PhantomData,
        }
    }
}

pub struct SmallestDegreeIter<const N: usize, CStrat = (), COpt = ()> {
    phantom: std::marker::PhantomData<(CStrat, COpt)>,
}

pub struct MinResultRankWith<
    const SCORE_ORDER: u128,
    const EXACT_JOIN_LIMIT: usize = DEFAULT_EXACT_JOIN_LIMIT,
    CStrat = (),
    COpt = (),
> {
    phantom: std::marker::PhantomData<(CStrat, COpt)>,
}

pub type MinResultRank<COpt = ()> =
    MinResultRankWith<PAIR_SCORE_SPARSE_ATOM_AWARE, DEFAULT_EXACT_JOIN_LIMIT, (), COpt>;

impl<const SCORE_ORDER: u128, const EXACT_JOIN_LIMIT: usize, CStrat, COpt> Default
    for MinResultRankWith<SCORE_ORDER, EXACT_JOIN_LIMIT, CStrat, COpt>
{
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

pub struct SingleSmallestDegree<const D: bool, CStrat = (), COpt = ()> {
    phantom: std::marker::PhantomData<(CStrat, COpt)>,
}

pub struct SingleLargestDegree<const D: bool, CStrat = (), COpt = ()> {
    phantom: std::marker::PhantomData<(CStrat, COpt)>,
}

pub enum ContractionMode {
    SmallestDegree,
    MinResultRank,
}

pub trait ContractionStrategy<E, L, K, FK, Aind>: Sized {
    const SUPPORTS_PARTIAL_GRAPH_REWRITE: bool = false;

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

    #[allow(clippy::result_large_err)]
    fn contract_product_in_place(
        executor: &mut E,
        graph: &mut NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
        ignored: &mut SuBitGraph,
    ) -> Result<bool, TensorNetworkError<K, FK>>
    where
        K: Display + Debug,
        FK: Display + Debug,
        Aind: AbsInd,
    {
        let replacement = Self::contract(executor, graph, operation, lib)?;
        graph
            .identify_subgraph_nodes_without_deleting_self_edges(
                operation.subgraph(),
                NetworkNode::Leaf(replacement),
                ignored,
            )
            .ok_or_else(|| {
                TensorNetworkError::Other(eyre!(
                    "ready operation subgraph did not contain any nodes"
                ))
            })?;
        graph.finish_deferred_node_identifications();
        Ok(true)
    }
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

struct ProductPairReplacement<K, Aind> {
    nodes: [NodeIndex; 2],
    position: usize,
    leaf: NetworkLeaf<K, Aind>,
}

enum ProductRewriteProgress {
    NoProgress,
    UpdatedProduct,
    CollapsedProduct,
}

impl ProductRewriteProgress {
    fn made_progress(&self) -> bool {
        !matches!(self, ProductRewriteProgress::NoProgress)
    }
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
        match &self.operands.get(operand)?.leaf {
            NetworkLeaf::LocalTensor(index) => Some(*index),
            NetworkLeaf::TensorTerm(term) if term.scalar.is_none() => Some(term.tensor),
            NetworkLeaf::TensorSum(_)
            | NetworkLeaf::TensorTerm(_)
            | NetworkLeaf::TensorTermSum(_)
            | NetworkLeaf::LibraryKey { .. }
            | NetworkLeaf::Scalar(_) => None,
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
            NetworkLeaf::TensorTerm(term) => Some(executor.tensor(term.tensor).structure()),
            NetworkLeaf::TensorTermSum(terms) => {
                Some(executor.tensor(terms.first()?.tensor).structure())
            }
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

    fn tensor_operand_profile<Sc, Store>(
        &self,
        executor: &Store,
        operand: usize,
    ) -> Option<(Option<usize>, TensorContractionProfile)>
    where
        Store: NetworkStoreAccess,
        Store::Tensor: FastTensorSumContractible<Sc>,
    {
        let terms = self.tensor_terms(operand)?;
        let single_tensor = (terms.len() == 1).then_some(terms[0].tensor);
        let mut profile = TensorContractionProfile {
            entries: 0,
            total_terms: 0,
            max_terms: 0,
            total_bytes: 0,
            max_bytes: 0,
            common_factor_count: 0,
            simple_tensor: true,
        };

        for term in terms {
            let term_profile = executor.tensor(term.tensor).contraction_profile();
            profile.entries = profile.entries.saturating_add(term_profile.entries);
            profile.total_terms = profile.total_terms.saturating_add(term_profile.total_terms);
            profile.max_terms = profile.max_terms.max(term_profile.max_terms);
            profile.total_bytes = profile.total_bytes.saturating_add(term_profile.total_bytes);
            profile.max_bytes = profile.max_bytes.max(term_profile.max_bytes);
            profile.common_factor_count = profile
                .common_factor_count
                .saturating_add(term_profile.common_factor_count);
            profile.simple_tensor &= term_profile.simple_tensor && term.scalar.is_none();
        }

        profile.entries = profile.entries.max(1);
        profile.total_terms = profile.total_terms.max(1);
        profile.max_terms = profile.max_terms.max(1);
        profile.total_bytes = profile.total_bytes.max(1);
        profile.max_bytes = profile.max_bytes.max(1);

        Some((single_tensor, profile))
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
                let tensor = graph.get_lib_data(lib, source)?;
                let index = executor.push_tensor(T::from(tensor));
                operand.leaf = NetworkLeaf::LocalTensor(index);
            }
        }
        Ok(())
    }

    fn multiply_scalar_indices<Sc, Store>(
        executor: &mut Store,
        left: Option<super::graph::ScalarRef>,
        right: Option<super::graph::ScalarRef>,
    ) -> Option<super::graph::ScalarRef>
    where
        Store: NetworkStoreAccess<Scalar = Sc>,
        Sc: for<'a> MulAssign<Sc::Ref<'a>> + Clone + Ref,
    {
        match (left, right) {
            (None, None) => None,
            (Some(index), None) | (None, Some(index)) => Some(index),
            (Some(left), Some(right)) => {
                let mut scalar = executor.scalar_ref(left).clone();
                scalar *= executor.scalar_ref(right).refer();
                Some(executor.push_scalar(scalar).into())
            }
        }
    }

    fn multiply_scalar_value_with_term<Sc, Store>(
        executor: &mut Store,
        scalar: &Sc,
        term: &TensorTerm,
    ) -> TensorTerm
    where
        Store: NetworkStoreAccess<Scalar = Sc>,
        Sc: for<'a> MulAssign<Sc::Ref<'a>> + Clone + Ref,
    {
        let scalar = match term.scalar {
            Some(index) => {
                let mut combined = scalar.clone();
                combined *= executor.scalar_ref(index).refer();
                executor.push_scalar(combined)
            }
            None => executor.push_scalar(scalar.clone()),
        };

        TensorTerm::scaled(term.tensor, scalar)
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
            let is_scalar = match &operand.leaf {
                NetworkLeaf::Scalar(index) => {
                    if let Some(accumulator) = &mut accumulator {
                        *accumulator *= executor.scalar_ref(*index).refer();
                    } else {
                        accumulator = Some(executor.scalar_ref(*index).clone());
                    }
                    true
                }
                NetworkLeaf::LocalTensor(index) => {
                    if let Some(scalar) = executor.tensor(*index).scalar_ref() {
                        if let Some(accumulator) = &mut accumulator {
                            *accumulator *= scalar;
                        } else {
                            accumulator =
                                Some(Sc::from(executor.tensor(*index).clone().scalar().unwrap()));
                        }
                        true
                    } else {
                        false
                    }
                }
                NetworkLeaf::TensorTerm(term) => {
                    if let Some(tensor_scalar) = executor.tensor(term.tensor).scalar_ref() {
                        if let Some(term_scalar) = term.scalar {
                            if let Some(accumulator) = &mut accumulator {
                                *accumulator *= executor.scalar_ref(term_scalar).refer();
                            } else {
                                accumulator = Some(executor.scalar_ref(term_scalar).clone());
                            }
                        }

                        if let Some(accumulator) = &mut accumulator {
                            *accumulator *= tensor_scalar;
                        } else {
                            accumulator = Some(Sc::from(
                                executor.tensor(term.tensor).clone().scalar().unwrap(),
                            ));
                        }
                        true
                    } else {
                        false
                    }
                }
                NetworkLeaf::TensorSum(_)
                | NetworkLeaf::TensorTermSum(_)
                | NetworkLeaf::LibraryKey { .. } => false,
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
                let index = executor.push_scalar(accumulator).into();
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
                let index = executor.push_scalar(accumulator).into();
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
    fn contract_scalars_in_place<LT, T, L, Sc, FK, Store>(
        &mut self,
        executor: &mut Store,
        graph: &mut NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
        ignored: &mut SuBitGraph,
    ) -> Result<ProductRewriteProgress, TensorNetworkError<K, FK>>
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
            let is_scalar = match &operand.leaf {
                NetworkLeaf::Scalar(index) => {
                    if let Some(accumulator) = &mut accumulator {
                        *accumulator *= executor.scalar_ref(*index).refer();
                    } else {
                        accumulator = Some(executor.scalar_ref(*index).clone());
                    }
                    true
                }
                NetworkLeaf::LocalTensor(index) => {
                    if let Some(scalar) = executor.tensor(*index).scalar_ref() {
                        if let Some(accumulator) = &mut accumulator {
                            *accumulator *= scalar;
                        } else {
                            accumulator =
                                Some(Sc::from(executor.tensor(*index).clone().scalar().unwrap()));
                        }
                        true
                    } else {
                        false
                    }
                }
                NetworkLeaf::TensorTerm(term) => {
                    if let Some(tensor_scalar) = executor.tensor(term.tensor).scalar_ref() {
                        if let Some(term_scalar) = term.scalar {
                            if let Some(accumulator) = &mut accumulator {
                                *accumulator *= executor.scalar_ref(term_scalar).refer();
                            } else {
                                accumulator = Some(executor.scalar_ref(term_scalar).clone());
                            }
                        }

                        if let Some(accumulator) = &mut accumulator {
                            *accumulator *= tensor_scalar;
                        } else {
                            accumulator = Some(Sc::from(
                                executor.tensor(term.tensor).clone().scalar().unwrap(),
                            ));
                        }
                        true
                    } else {
                        false
                    }
                }
                NetworkLeaf::TensorSum(_)
                | NetworkLeaf::TensorTermSum(_)
                | NetworkLeaf::LibraryKey { .. } => false,
            };

            if is_scalar {
                scalar_positions.push(position);
            } else {
                tensor_positions.push(position);
            }
        }

        let collapse_product = |graph: &mut NetworkGraph<K, FK, Aind>,
                                ignored: &mut SuBitGraph,
                                leaf: NetworkLeaf<K, Aind>|
         -> Result<(), TensorNetworkError<K, FK>> {
            graph
                .identify_subgraph_nodes_without_deleting_self_edges(
                    operation.subgraph(),
                    NetworkNode::Leaf(leaf),
                    ignored,
                )
                .ok_or_else(|| {
                    TensorNetworkError::Other(eyre!(
                        "ready operation subgraph did not contain any nodes"
                    ))
                })?;
            graph.finish_deferred_node_identifications();
            Ok(())
        };

        let Some(accumulator) = accumulator else {
            if self.operands.len() == 1 {
                collapse_product(graph, ignored, self.operands[0].leaf.clone())?;
                return Ok(ProductRewriteProgress::CollapsedProduct);
            }
            return Ok(ProductRewriteProgress::NoProgress);
        };

        match tensor_positions.len() {
            0 => {
                let leaf = if scalar_positions.len() == 1
                    && matches!(
                        self.operands[scalar_positions[0]].leaf,
                        NetworkLeaf::Scalar(_)
                    ) {
                    self.operands[scalar_positions[0]].leaf.clone()
                } else {
                    NetworkLeaf::Scalar(executor.push_scalar(accumulator).into())
                };
                collapse_product(graph, ignored, leaf)?;
                Ok(ProductRewriteProgress::CollapsedProduct)
            }
            1 => {
                let leaf = self.scalar_multiply_operand::<LT, T, L, Sc, FK, Store>(
                    tensor_positions[0],
                    accumulator,
                    executor,
                    graph,
                    lib,
                )?;
                collapse_product(graph, ignored, leaf)?;
                Ok(ProductRewriteProgress::CollapsedProduct)
            }
            _ if scalar_positions.len() > 1 => {
                let scalar = executor.push_scalar(accumulator).into();
                let scalar_nodes = scalar_positions
                    .iter()
                    .map(|position| {
                        self.operands[*position].source.ok_or_else(|| {
                            TensorNetworkError::Other(eyre!(
                                "scalar product operand lost source node"
                            ))
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let node = graph.identify_nodes_marking_self_edges_and_duplicate_heads(
                    &scalar_nodes,
                    NetworkNode::Leaf(NetworkLeaf::Scalar(scalar)),
                    ignored,
                );
                graph.finish_deferred_node_identifications();
                self.replace_operands(
                    &scalar_positions,
                    ProductOperand {
                        leaf: NetworkLeaf::Scalar(scalar),
                        source: Some(node),
                    },
                );
                Ok(ProductRewriteProgress::UpdatedProduct)
            }
            _ => Ok(ProductRewriteProgress::NoProgress),
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
        Sc: for<'a> MulAssign<Sc::Ref<'a>> + Clone + Ref,
    {
        let leaf = match &self.operands[operand].leaf {
            NetworkLeaf::LocalTensor(index) => {
                let scalar = executor.push_scalar(scalar);
                return Ok(NetworkLeaf::TensorTerm(TensorTerm::scaled(*index, scalar)));
            }
            NetworkLeaf::TensorSum(indices) => {
                let scalar = executor.push_scalar(scalar);
                let mut scaled = Vec::with_capacity(indices.len());
                for index in indices {
                    scaled.push(TensorTerm::scaled(*index, scalar));
                }
                return Ok(NetworkLeaf::TensorTermSum(scaled));
            }
            NetworkLeaf::TensorTerm(term) => {
                let term = Self::multiply_scalar_value_with_term(executor, &scalar, term);
                return Ok(NetworkLeaf::TensorTerm(term));
            }
            NetworkLeaf::TensorTermSum(terms) => {
                let scaled = terms
                    .iter()
                    .map(|term| Self::multiply_scalar_value_with_term(executor, &scalar, term))
                    .collect();
                return Ok(NetworkLeaf::TensorTermSum(scaled));
            }
            NetworkLeaf::LibraryKey { .. } => {
                let source = self.operands[operand].source.ok_or_else(|| {
                    TensorNetworkError::Other(eyre!("library product operand lost source node"))
                })?;
                let tensor = graph.get_lib_data(lib, source)?;
                tensor
                    .scalar_mul(&scalar)
                    .ok_or(TensorNetworkError::FailedScalarMul)?
            }
            NetworkLeaf::Scalar(_) => return Err(TensorNetworkError::SlotEdgeToScalarNode),
        };

        let index = executor.push_tensor(leaf);
        Ok(NetworkLeaf::LocalTensor(index))
    }

    fn tensor_terms(&self, operand: usize) -> Option<Vec<TensorTerm>> {
        match &self.operands.get(operand)?.leaf {
            NetworkLeaf::LocalTensor(index) => Some(vec![TensorTerm::tensor(*index)]),
            NetworkLeaf::TensorSum(indices) => {
                Some(indices.iter().copied().map(TensorTerm::tensor).collect())
            }
            NetworkLeaf::TensorTerm(term) => Some(vec![term.clone()]),
            NetworkLeaf::TensorTermSum(terms) => Some(terms.clone()),
            NetworkLeaf::LibraryKey { .. } | NetworkLeaf::Scalar(_) => None,
        }
    }

    fn tensor_sum_leaf<T, Store>(
        executor: &mut Store,
        terms: Vec<TensorTerm>,
    ) -> NetworkLeaf<K, Aind>
    where
        Store: NetworkStoreAccess<Tensor = T>,
        T: HasStructure + Clone + Ref + FastTensorSum + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    {
        debug_assert!(!terms.is_empty());

        if terms.len() == 1 {
            let term = terms.into_iter().next().expect("single tensor term");
            return match term.scalar {
                Some(_) => NetworkLeaf::TensorTerm(term),
                None => NetworkLeaf::LocalTensor(term.tensor),
            };
        }

        if terms.iter().all(|term| term.scalar.is_none())
            && terms
                .iter()
                .all(|term| executor.tensor(term.tensor).scalar_ref().is_some())
        {
            let mut iter = terms.into_iter();
            let first = iter
                .next()
                .expect("tensor sum with at least one term has first term")
                .tensor;
            let mut materialized = executor.tensor(first).clone();
            for term in iter {
                materialized += executor.tensor(term.tensor).refer();
            }
            NetworkLeaf::LocalTensor(executor.push_tensor(materialized))
        } else if terms.iter().all(|term| term.scalar.is_none()) {
            NetworkLeaf::TensorSum(terms.into_iter().map(|term| term.tensor).collect())
        } else {
            NetworkLeaf::TensorTermSum(terms)
        }
    }

    fn materialize_tensor_terms<T, Sc, Store>(executor: &mut Store, terms: &[TensorTerm]) -> usize
    where
        Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
        T: Clone
            + Ref
            + FastTensorSum
            + ScalarMul<Sc, Output = T>
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
        Sc: Clone,
    {
        debug_assert!(!terms.is_empty());

        if terms.iter().all(|term| term.scalar.is_none()) {
            let indices = terms.iter().map(|term| term.tensor).collect::<Vec<_>>();
            return Self::materialize_tensor_sum(executor, &indices);
        }

        let mut materialized_terms = Vec::with_capacity(terms.len());
        for term in terms {
            let tensor = match term.scalar {
                Some(scalar) => executor
                    .tensor(term.tensor)
                    .scalar_mul(executor.scalar_ref(scalar))
                    .expect("scaled tensor term should support scalar multiplication"),
                None => executor.tensor(term.tensor).clone(),
            };
            materialized_terms.push(executor.push_tensor(tensor));
        }

        Self::materialize_tensor_sum(executor, &materialized_terms)
    }

    fn materialize_tensor_sum<T, Store>(executor: &mut Store, indices: &[usize]) -> usize
    where
        Store: NetworkStoreAccess<Tensor = T>,
        T: Clone + Ref + FastTensorSum + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    {
        debug_assert!(!indices.is_empty());

        if let Some(materialized) = {
            let terms = indices
                .iter()
                .map(|index| executor.tensor(*index))
                .collect::<Vec<_>>();
            T::fast_tensor_sum(&terms, None)
        } {
            return executor.push_tensor(materialized);
        }

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

    fn split_common_tensor_factor<T, Sc, Store>(
        executor: &mut Store,
        term: TensorTerm,
    ) -> TensorTerm
    where
        Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
        T: TensorCommonFactor<Sc>,
        Sc: for<'a> MulAssign<Sc::Ref<'a>> + Clone + Ref,
    {
        let Some((reduced, factor)) = executor.tensor(term.tensor).split_common_factor() else {
            return term;
        };

        let tensor = executor.push_tensor(reduced);
        let factor = executor.push_scalar(factor);
        let scalar = Self::multiply_scalar_indices(executor, term.scalar, Some(factor.into()));

        TensorTerm { tensor, scalar }
    }

    fn try_fast_tensor_sum_contract<T, Sc, Store, CStrat, COpt, FK>(
        executor: &mut Store,
        left_terms: &[TensorTerm],
        right_terms: &[TensorTerm],
    ) -> Result<Option<NetworkLeaf<K, Aind>>, TensorNetworkError<K, FK>>
    where
        Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
        T: HasStructure
            + Clone
            + Ref
            + FastTensorSum
            + FastTensorSumContractible<Sc>
            + AtomComponentOptimizable<COpt>
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
        Sc: Clone,
        K: Display,
        FK: Display,
    {
        if left_terms.len() > 1
            && right_terms.len() == 1
            && left_terms.iter().all(|term| term.scalar.is_none())
            && right_terms[0].scalar.is_none()
        {
            let terms = left_terms
                .iter()
                .map(|term| executor.tensor(term.tensor))
                .collect::<Vec<_>>();
            if let Some(result) = T::fast_tensor_sum_contract::<CStrat>(
                &terms,
                executor.tensor(right_terms[0].tensor),
                true,
            ) {
                return result
                    .map(|result| {
                        Some(Self::fast_tensor_sum_contract_leaf::<T, Sc, Store, COpt>(
                            executor, result,
                        ))
                    })
                    .map_err(|error| TensorNetworkError::Other(error.into()));
            }
        }

        if right_terms.len() > 1
            && left_terms.len() == 1
            && right_terms.iter().all(|term| term.scalar.is_none())
            && left_terms[0].scalar.is_none()
        {
            let terms = right_terms
                .iter()
                .map(|term| executor.tensor(term.tensor))
                .collect::<Vec<_>>();
            if let Some(result) = T::fast_tensor_sum_contract::<CStrat>(
                &terms,
                executor.tensor(left_terms[0].tensor),
                false,
            ) {
                return result
                    .map(|result| {
                        Some(Self::fast_tensor_sum_contract_leaf::<T, Sc, Store, COpt>(
                            executor, result,
                        ))
                    })
                    .map_err(|error| TensorNetworkError::Other(error.into()));
            }
        }

        Ok(None)
    }

    fn fast_tensor_sum_contract_leaf<T, Sc, Store, COpt>(
        executor: &mut Store,
        result: FastTensorSumContract<T, Sc>,
    ) -> NetworkLeaf<K, Aind>
    where
        Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
        T: HasStructure
            + Clone
            + Ref
            + FastTensorSum
            + AtomComponentOptimizable<COpt>
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    {
        match result {
            FastTensorSumContract::Materialized(tensor) => {
                NetworkLeaf::LocalTensor(executor.push_tensor(tensor.optimize_atom_components()))
            }
            FastTensorSumContract::Terms(terms) => {
                let terms = terms
                    .into_iter()
                    .map(|tensor| {
                        TensorTerm::tensor(executor.push_tensor(tensor.optimize_atom_components()))
                    })
                    .collect::<Vec<_>>();
                Self::tensor_sum_leaf(executor, terms)
            }
            FastTensorSumContract::ScaledTerms(terms) => {
                let terms = terms
                    .into_iter()
                    .map(|term| TensorTerm {
                        tensor: executor.push_tensor(term.tensor.optimize_atom_components()),
                        scalar: term
                            .scalar
                            .map(|scalar| executor.push_scalar(scalar).into()),
                    })
                    .collect::<Vec<_>>();
                Self::tensor_sum_leaf(executor, terms)
            }
        }
    }

    #[allow(clippy::result_large_err)]
    pub fn contract_pair<LT, T, L, Sc, CStrat, COpt, FK, Store>(
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
            + AtomComponentOptimizable<COpt>
            + Clone
            + Ref
            + FastTensorSum
            + FastTensorSumContractible<Sc>
            + TensorCommonFactor<Sc>
            + ScalarMul<Sc, Output = T>
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        Sc: for<'a> MulAssign<Sc::Ref<'a>> + Clone + Ref,
        T::Structure: Display,
        LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        self.materialize_libraries::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;

        let mut left_terms = self
            .tensor_terms(left)
            .ok_or(TensorNetworkError::SlotEdgeToScalarNode)?;
        let mut right_terms = self
            .tensor_terms(right)
            .ok_or(TensorNetworkError::SlotEdgeToScalarNode)?;
        left_terms = left_terms
            .into_iter()
            .map(|term| Self::split_common_tensor_factor::<T, Sc, Store>(executor, term))
            .collect();
        right_terms = right_terms
            .into_iter()
            .map(|term| Self::split_common_tensor_factor::<T, Sc, Store>(executor, term))
            .collect();
        let distributed_terms = left_terms.len() * right_terms.len();

        let max_distributed_terms = max_lazy_tensor_sum_distributed_terms();
        if distributed_terms > max_distributed_terms
            && (left_terms.len() > 1 || right_terms.len() > 1)
        {
            if profile::enabled() {
                eprintln!(
                    "spenso_profile product.lazy_tensor_sum_materialize left_terms={} right_terms={} distributed_terms={} max_distributed_terms={}",
                    left_terms.len(),
                    right_terms.len(),
                    distributed_terms,
                    max_distributed_terms,
                );
            }
            if left_terms.len() > 1 {
                left_terms = vec![TensorTerm::tensor(Self::materialize_tensor_terms(
                    executor,
                    &left_terms,
                ))];
            }
            if right_terms.len() > 1 {
                right_terms = vec![TensorTerm::tensor(Self::materialize_tensor_terms(
                    executor,
                    &right_terms,
                ))];
            }
        }

        if profile::enabled() && (left_terms.len() > 1 || right_terms.len() > 1) {
            eprintln!(
                "spenso_profile product.lazy_tensor_sum left_terms={} right_terms={} distributed_terms={} left_structure={} right_structure={}",
                left_terms.len(),
                right_terms.len(),
                left_terms.len() * right_terms.len(),
                self.tensor_structure(executor, left)
                    .map(|structure| structure.to_string())
                    .unwrap_or_else(|| "<none>".to_string()),
                self.tensor_structure(executor, right)
                    .map(|structure| structure.to_string())
                    .unwrap_or_else(|| "<none>".to_string()),
            );
        }

        if let Some(leaf) = Self::try_fast_tensor_sum_contract::<T, Sc, Store, CStrat, COpt, FK>(
            executor,
            &left_terms,
            &right_terms,
        )? {
            let mut positions = [left, right];
            positions.sort_unstable();
            self.replace_operands(&positions, ProductOperand { leaf, source: None });
            return Ok(());
        }

        if left_terms.len() == 1
            && right_terms.len() == 1
            && left_terms[0].scalar.is_none()
            && right_terms[0].scalar.is_none()
        {
            if let Some(result) = T::fast_tensor_sum_contract::<CStrat>(
                &[executor.tensor(left_terms[0].tensor)],
                executor.tensor(right_terms[0].tensor),
                true,
            ) {
                let leaf = result
                    .map(|result| {
                        Self::fast_tensor_sum_contract_leaf::<T, Sc, Store, COpt>(executor, result)
                    })
                    .map_err(|error| TensorNetworkError::Other(error.into()))?;
                let mut positions = [left, right];
                positions.sort_unstable();
                self.replace_operands(&positions, ProductOperand { leaf, source: None });
                return Ok(());
            }

            if let Some(result) = T::fast_tensor_sum_contract::<CStrat>(
                &[executor.tensor(right_terms[0].tensor)],
                executor.tensor(left_terms[0].tensor),
                false,
            ) {
                let leaf = result
                    .map(|result| {
                        Self::fast_tensor_sum_contract_leaf::<T, Sc, Store, COpt>(executor, result)
                    })
                    .map_err(|error| TensorNetworkError::Other(error.into()))?;
                let mut positions = [left, right];
                positions.sort_unstable();
                self.replace_operands(&positions, ProductOperand { leaf, source: None });
                return Ok(());
            }
        }

        let mut contracted_terms = Vec::with_capacity(left_terms.len() * right_terms.len());
        for left_tensor in left_terms {
            for right_tensor in &right_terms {
                let contracted = executor
                    .tensor(left_tensor.tensor)
                    .contract(executor.tensor(right_tensor.tensor))?;
                let (contracted, extra_scalar) =
                    if let Some((reduced, factor)) = contracted.split_common_factor() {
                        (reduced, Some(executor.push_scalar(factor).into()))
                    } else {
                        (contracted, None)
                    };
                let contracted = contracted.optimize_atom_components();
                let tensor = executor.push_tensor(contracted);
                let scalar = Self::multiply_scalar_indices(
                    executor,
                    left_tensor.scalar,
                    right_tensor.scalar,
                );
                let scalar = Self::multiply_scalar_indices(executor, scalar, extra_scalar);
                contracted_terms.push(TensorTerm { tensor, scalar });
            }
        }
        let leaf = Self::tensor_sum_leaf(executor, contracted_terms);

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
        COpt,
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
            + AtomComponentOptimizable<COpt>
            + Clone
            + Ref
            + FastTensorSum
            + FastTensorSumContractible<Sc>
            + TensorCommonFactor<Sc>
            + ScalarMul<Sc, Output = T>
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
        T::Structure: Display,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        Sc: for<'a> MulAssign<Sc::Ref<'a>> + Clone + Ref,
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

        self.contract_pair::<LT, T, L, Sc, CStrat, COpt, FK, Store>(
            left, right, executor, graph, lib,
        )?;
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

    #[allow(clippy::result_large_err)]
    fn contract_one_by_degree_replacement<
        const DEBUG: bool,
        const LARGEST: bool,
        LT,
        T,
        L,
        Sc,
        CStrat,
        COpt,
        FK,
        Store,
    >(
        &mut self,
        executor: &mut Store,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<Option<ProductPairReplacement<K, Aind>>, TensorNetworkError<K, FK>>
    where
        Store: NetworkStoreAccess<Tensor = T, Scalar = Sc>,
        K: Clone + Debug + Display,
        FK: Debug + Display,
        Aind: AbsInd,
        LT: LibraryTensor + Clone,
        T: HasStructure
            + From<LT::WithIndices>
            + Contract<T, CStrat, LCM = T>
            + AtomComponentOptimizable<COpt>
            + Clone
            + Ref
            + FastTensorSum
            + FastTensorSumContractible<Sc>
            + TensorCommonFactor<Sc>
            + ScalarMul<Sc, Output = T>
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
        T::Structure: Display,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        Sc: for<'a> MulAssign<Sc::Ref<'a>> + Clone + Ref,
        LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        self.materialize_libraries::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        let Some((left, right, degree)) = self.best_degree_pair::<Store, LARGEST>(executor) else {
            return Ok(None);
        };
        let left_node = self.operands[left].source.ok_or_else(|| {
            TensorNetworkError::Other(eyre!("tensor product operand lost source node"))
        })?;
        let right_node = self.operands[right].source.ok_or_else(|| {
            TensorNetworkError::Other(eyre!("tensor product operand lost source node"))
        })?;

        if DEBUG {
            println!(
                "Contracting {} with {}",
                self.tensor_structure(executor, left).unwrap(),
                self.tensor_structure(executor, right).unwrap()
            );
        }
        let profile_pair = profile::enabled();
        let log_pair = profile::verbose() && (self.operands.len() > 4 || degree > 1);
        let pair_start = if profile_pair {
            if log_pair {
                eprintln!(
                    "spenso_profile product.pair_start operands={} left_operand={} right_operand={} degree={} left={} right={}",
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

        self.contract_pair::<LT, T, L, Sc, CStrat, COpt, FK, Store>(
            left, right, executor, graph, lib,
        )?;
        if let Some(pair_start) = pair_start {
            let elapsed = pair_start.elapsed();
            if profile::verbose() || elapsed.as_millis() >= 100 {
                eprintln!(
                    "spenso_profile product.slow_pair elapsed_ms={:.3}",
                    elapsed.as_secs_f64() * 1000.0,
                );
            }
        }

        let replacement_position = left.min(right);
        if DEBUG && let Some(structure) = self.tensor_structure(executor, replacement_position) {
            println!("Obtained {structure}");
        }

        Ok(Some(ProductPairReplacement {
            nodes: [left_node, right_node],
            position: replacement_position,
            leaf: self.operands[replacement_position].leaf.clone(),
        }))
    }

    #[allow(clippy::result_large_err)]
    fn contract_one_by_degree_in_place<
        const DEBUG: bool,
        const LARGEST: bool,
        LT,
        T,
        L,
        Sc,
        CStrat,
        COpt,
        FK,
        Store,
    >(
        &mut self,
        executor: &mut Store,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        ignored: &mut SuBitGraph,
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
            + AtomComponentOptimizable<COpt>
            + Clone
            + Ref
            + FastTensorSum
            + FastTensorSumContractible<Sc>
            + TensorCommonFactor<Sc>
            + ScalarMul<Sc, Output = T>
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
        T::Structure: Display,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        Sc: for<'a> MulAssign<Sc::Ref<'a>> + Clone + Ref,
        LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        let Some(replacement) = self
            .contract_one_by_degree_replacement::<
                DEBUG,
                LARGEST,
                LT,
                T,
                L,
                Sc,
                CStrat,
                COpt,
                FK,
                Store,
            >(executor, graph, lib)?
        else {
            return Ok(false);
        };

        let node = graph.identify_nodes_marking_self_edges_and_duplicate_heads(
            &replacement.nodes,
            NetworkNode::Leaf(replacement.leaf),
            ignored,
        );
        graph.finish_deferred_node_identifications();
        self.operands[replacement.position].source = Some(node);
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

    fn best_result_rank_pair_by_degree<
        Sc,
        Store,
        const SCORE_ORDER: u128,
        const EXACT_JOIN_LIMIT: usize,
    >(
        &self,
        executor: &Store,
        mut accept_degree: impl FnMut(u32) -> bool,
    ) -> Option<(usize, usize, u32)>
    where
        Store: NetworkStoreAccess,
        Store::Tensor: HasStructure + FastTensorSumContractible<Sc>,
    {
        let tensors = self
            .operands
            .iter()
            .enumerate()
            .filter_map(|(operand, _)| {
                let structure = self.tensor_structure(executor, operand)?;
                let (single_tensor, profile) =
                    self.tensor_operand_profile::<Sc, Store>(executor, operand)?;
                Some((operand, single_tensor, structure, profile))
            })
            .collect::<Vec<_>>();

        let mut best = None;
        for (left_position, (left_operand, left_tensor, left_structure, left_profile)) in
            tensors.iter().enumerate()
        {
            for (right_operand, right_tensor, right_structure, right_profile) in
                &tensors[left_position + 1..]
            {
                let Some((left_match_permutation, left_matches, right_matches)) =
                    left_structure.match_indices(right_structure)
                else {
                    continue;
                };
                let degree = left_matches.iter().filter(|matched| **matched).count() as u32;
                if !accept_degree(degree) {
                    continue;
                }

                let result_rank = left_matches.iter().filter(|matched| !**matched).count()
                    + right_matches.iter().filter(|matched| !**matched).count();
                let output_dense_size = if pair_score_order_needs_output_dense_size(SCORE_ORDER) {
                    output_dense_size(
                        *left_structure,
                        *right_structure,
                        &left_matches,
                        &right_matches,
                    )
                } else {
                    1
                };
                let pair_estimate = if pair_score_order_needs_sparse_estimate(SCORE_ORDER) {
                    match (left_tensor, right_tensor) {
                        (Some(left_tensor), Some(right_tensor)) => {
                            executor.tensor(*left_tensor).contraction_pair_estimate(
                                executor.tensor(*right_tensor),
                                &left_match_permutation,
                                &left_matches,
                                &right_matches,
                                *left_profile,
                                *right_profile,
                                output_dense_size,
                                EXACT_JOIN_LIMIT,
                            )
                        }
                        _ => super::TensorContractionPairEstimate::from_profiles(
                            *left_profile,
                            *right_profile,
                            output_dense_size,
                        ),
                    }
                } else {
                    super::TensorContractionPairEstimate::from_profiles(
                        *left_profile,
                        *right_profile,
                        output_dense_size,
                    )
                };
                let score = ResultRankPairScore::new(
                    result_rank as u32,
                    degree,
                    *left_profile,
                    *right_profile,
                    pair_estimate,
                );

                let replace = match &best {
                    None => true,
                    Some((_, _, _, best_score)) => score.less_by::<SCORE_ORDER>(best_score),
                };

                if replace {
                    best = Some((*left_operand, *right_operand, degree, score));
                }
            }
        }

        best.map(|(left, right, degree, _)| (left, right, degree))
    }

    pub fn best_result_rank_pair<
        Sc,
        Store,
        const SCORE_ORDER: u128,
        const EXACT_JOIN_LIMIT: usize,
    >(
        &self,
        executor: &Store,
    ) -> Option<(usize, usize, u32)>
    where
        Store: NetworkStoreAccess,
        Store::Tensor: HasStructure + FastTensorSumContractible<Sc>,
    {
        self.best_result_rank_pair_by_degree::<Sc, Store, SCORE_ORDER, EXACT_JOIN_LIMIT>(
            executor,
            |degree| degree > 0,
        )
    }

    pub fn best_exterior_result_rank_pair<
        Sc,
        Store,
        const SCORE_ORDER: u128,
        const EXACT_JOIN_LIMIT: usize,
    >(
        &self,
        executor: &Store,
    ) -> Option<(usize, usize, u32)>
    where
        Store: NetworkStoreAccess,
        Store::Tensor: HasStructure + FastTensorSumContractible<Sc>,
        <Store::Tensor as HasStructure>::Structure: StructureContract,
    {
        let tensors = self
            .operands
            .iter()
            .enumerate()
            .filter_map(|(operand, _)| {
                let structure = self.tensor_structure(executor, operand)?;
                let (_, profile) = self.tensor_operand_profile::<Sc, Store>(executor, operand)?;
                Some((operand, structure, profile))
            })
            .collect::<Vec<_>>();

        let mut best = None;
        for (left_position, (left_operand, left_structure, left_profile)) in
            tensors.iter().enumerate()
        {
            for (right_operand, right_structure, right_profile) in &tensors[left_position + 1..] {
                let Ok((result_structure, left_matches, _, _)) =
                    left_structure.merge(right_structure)
                else {
                    continue;
                };
                let degree = left_matches.n_included() as u32;
                if degree != 0 {
                    continue;
                }

                let output_dense_size = if pair_score_order_needs_output_dense_size(SCORE_ORDER) {
                    result_structure
                        .size()
                        .map(|size| size as u128)
                        .unwrap_or(u128::MAX)
                        .max(1)
                } else {
                    1
                };
                let pair_estimate = super::TensorContractionPairEstimate::from_profiles(
                    *left_profile,
                    *right_profile,
                    output_dense_size,
                );
                let score = ResultRankPairScore::new(
                    result_structure.order() as u32,
                    degree,
                    *left_profile,
                    *right_profile,
                    pair_estimate,
                );

                let replace = match &best {
                    None => true,
                    Some((_, _, _, best_score)) => score.less_by::<SCORE_ORDER>(best_score),
                };

                if replace {
                    best = Some((*left_operand, *right_operand, degree, score));
                }
            }
        }

        best.map(|(left, right, degree, _)| (left, right, degree))
    }

    #[allow(clippy::result_large_err)]
    pub fn contract_one_by_result_rank<
        LT,
        T,
        L,
        Sc,
        CStrat,
        COpt,
        FK,
        Store,
        const SCORE_ORDER: u128,
        const EXACT_JOIN_LIMIT: usize,
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
            + AtomComponentOptimizable<COpt>
            + Clone
            + Ref
            + FastTensorSum
            + FastTensorSumContractible<Sc>
            + TensorCommonFactor<Sc>
            + ScalarMul<Sc, Output = T>
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
        T::Structure: Display,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        Sc: for<'a> MulAssign<Sc::Ref<'a>> + Clone + Ref,
        LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        self.materialize_libraries::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        let Some((left, right, degree)) =
            self.best_result_rank_pair::<Sc, Store, SCORE_ORDER, EXACT_JOIN_LIMIT>(executor)
        else {
            return Ok(false);
        };

        self.contract_result_rank_pair::<LT, T, L, Sc, CStrat, COpt, FK, Store>(
            left,
            right,
            degree,
            "result_rank_pair_start",
            "result_rank_slow_pair",
            executor,
            graph,
            lib,
        )?;

        Ok(true)
    }

    #[allow(clippy::result_large_err)]
    pub fn contract_one_exterior_by_result_rank<
        LT,
        T,
        L,
        Sc,
        CStrat,
        COpt,
        FK,
        Store,
        const SCORE_ORDER: u128,
        const EXACT_JOIN_LIMIT: usize,
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
            + AtomComponentOptimizable<COpt>
            + Clone
            + Ref
            + FastTensorSum
            + FastTensorSumContractible<Sc>
            + TensorCommonFactor<Sc>
            + ScalarMul<Sc, Output = T>
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
        T::Structure: Display + StructureContract,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        Sc: for<'a> MulAssign<Sc::Ref<'a>> + Clone + Ref,
        LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        self.materialize_libraries::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        let Some((left, right, degree)) = self
            .best_exterior_result_rank_pair::<Sc, Store, SCORE_ORDER, EXACT_JOIN_LIMIT>(executor)
        else {
            return Ok(false);
        };

        self.contract_result_rank_pair::<LT, T, L, Sc, CStrat, COpt, FK, Store>(
            left,
            right,
            degree,
            "result_rank_exterior_pair_start",
            "result_rank_exterior_slow_pair",
            executor,
            graph,
            lib,
        )?;

        Ok(true)
    }

    #[allow(clippy::result_large_err, clippy::too_many_arguments)]
    fn contract_result_rank_pair<LT, T, L, Sc, CStrat, COpt, FK, Store>(
        &mut self,
        left: usize,
        right: usize,
        degree: u32,
        start_event: &str,
        slow_event: &str,
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
            + AtomComponentOptimizable<COpt>
            + Clone
            + Ref
            + FastTensorSum
            + FastTensorSumContractible<Sc>
            + TensorCommonFactor<Sc>
            + ScalarMul<Sc, Output = T>
            + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
        T::Structure: Display,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        Sc: for<'a> MulAssign<Sc::Ref<'a>> + Clone + Ref,
        LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        let profile_pair = profile::enabled();
        let left_profile = profile_pair
            .then(|| self.tensor_operand_profile::<Sc, Store>(executor, left))
            .flatten()
            .map(|(_, profile)| profile);
        let right_profile = profile_pair
            .then(|| self.tensor_operand_profile::<Sc, Store>(executor, right))
            .flatten()
            .map(|(_, profile)| profile);
        let large_profile_pair = match (left_profile, right_profile) {
            (Some(left_profile), Some(right_profile)) => {
                left_profile.entries.saturating_mul(right_profile.entries) > 100_000
                    || left_profile
                        .total_bytes
                        .saturating_add(right_profile.total_bytes)
                        > 100_000_000
            }
            _ => false,
        };
        let log_pair = (profile_pair && large_profile_pair)
            || (profile::verbose() && (self.operands.len() > 4 || degree > 1));
        let pair_start = if profile_pair {
            if log_pair {
                eprintln!(
                    "spenso_profile product.{start_event} operands={} left_operand={} right_operand={} degree={} left={} right={} left_profile={:?} right_profile={:?}",
                    self.operands.len(),
                    left,
                    right,
                    degree,
                    self.tensor_structure(executor, left).unwrap(),
                    self.tensor_structure(executor, right).unwrap(),
                    left_profile,
                    right_profile,
                );
            }
            Some(std::time::Instant::now())
        } else {
            None
        };

        self.contract_pair::<LT, T, L, Sc, CStrat, COpt, FK, Store>(
            left, right, executor, graph, lib,
        )?;

        if let Some(pair_start) = pair_start {
            let elapsed = pair_start.elapsed();
            if profile::verbose() || elapsed.as_millis() >= 100 {
                eprintln!(
                    "spenso_profile product.{slow_event} elapsed_ms={:.3}",
                    elapsed.as_secs_f64() * 1000.0,
                );
            }
        }

        Ok(())
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
            let operands = self
                .operands
                .iter()
                .enumerate()
                .map(|(position, operand)| {
                    let kind = match &operand.leaf {
                        NetworkLeaf::LocalTensor(_) => "LocalTensor".to_string(),
                        NetworkLeaf::TensorSum(indices) => format!("TensorSum({})", indices.len()),
                        NetworkLeaf::TensorTerm(term) => {
                            format!("TensorTerm(scaled={})", term.scalar.is_some())
                        }
                        NetworkLeaf::TensorTermSum(terms) => {
                            let scaled = terms.iter().filter(|term| term.scalar.is_some()).count();
                            format!("TensorTermSum(terms={},scaled={scaled})", terms.len())
                        }
                        NetworkLeaf::LibraryKey { .. } => "LibraryKey".to_string(),
                        NetworkLeaf::Scalar(_) => "Scalar".to_string(),
                    };
                    let structure = self
                        .tensor_structure(executor, position)
                        .map(|structure| {
                            format!(
                                "order={},dense_size={}",
                                structure.order(),
                                structure.size().unwrap_or(0)
                            )
                        })
                        .unwrap_or_else(|| "no_tensor_structure".to_string());
                    format!("#{position}:{kind}:{structure}:source={:?}", operand.source)
                })
                .collect::<Vec<_>>();
            let mut matching_pairs = Vec::new();
            for left in 0..self.operands.len() {
                let Some(left_structure) = self.tensor_structure(executor, left) else {
                    continue;
                };
                for right in left + 1..self.operands.len() {
                    let Some(right_structure) = self.tensor_structure(executor, right) else {
                        continue;
                    };
                    let Some((_, left_matches, _)) = left_structure.match_indices(right_structure)
                    else {
                        continue;
                    };
                    let degree = left_matches.iter().filter(|matched| **matched).count();
                    if degree > 0 {
                        matching_pairs.push(format!("#{left}-#{right}:degree={degree}"));
                    }
                }
            }

            Err(TensorNetworkError::Other(eyre!(
                "product contraction did not collapse to one leaf; {} operands remain; operands=[{}]; matching_pairs=[{}]",
                self.operands.len(),
                operands.join(", "),
                matching_pairs.join(", ")
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

fn output_dense_size<S>(left: &S, right: &S, left_matches: &[bool], right_matches: &[bool]) -> u128
where
    S: TensorStructure,
{
    let left_size = left
        .external_dims_iter()
        .enumerate()
        .filter(|(axis, _)| !left_matches.get(*axis).copied().unwrap_or(false))
        .map(|(_, dim)| usize::try_from(dim).unwrap_or(usize::MAX) as u128)
        .fold(1u128, u128::saturating_mul);
    let right_size = right
        .external_dims_iter()
        .enumerate()
        .filter(|(axis, _)| !right_matches.get(*axis).copied().unwrap_or(false))
        .map(|(_, dim)| usize::try_from(dim).unwrap_or(usize::MAX) as u128)
        .fold(1u128, u128::saturating_mul);
    left_size.saturating_mul(right_size).max(1)
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
        + FastTensorSum
        + FastTensorSumContractible<Sc>
        + TensorCommonFactor<Sc>
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
    const SUPPORTS_PARTIAL_GRAPH_REWRITE: bool = true;

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

    fn contract_product_in_place(
        executor: &mut Store,
        graph: &mut NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
        ignored: &mut SuBitGraph,
    ) -> Result<bool, TensorNetworkError<K, FK>>
    where
        K: Display + Debug,
        FK: Display + Debug,
        Aind: AbsInd,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product
            .contract_scalars_in_place::<LT, T, L, Sc, FK, Store>(
                executor, graph, operation, lib, ignored,
            )
            .map(|progress| progress.made_progress())
    }
}

impl<
    CStrat,
    COpt,
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<T, CStrat, LCM = T>
        + AtomComponentOptimizable<COpt>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>
        + Ref
        + FastTensorSum
        + FastTensorSumContractible<Sc>
        + TensorCommonFactor<Sc>
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
> ContractionStrategy<Store, L, K, FK, Aind> for SmallestDegree<CStrat, COpt>
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
            .contract_one_by_degree::<false, false, LT, T, L, Sc, CStrat, COpt, FK, Store>(
                executor, graph, lib,
            )?
        {
            product.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        }
        product.finish::<LT, T, L, Sc, FK, Store>(executor, graph, lib)
    }
}

impl<
    CStrat,
    COpt,
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<T, CStrat, LCM = T>
        + AtomComponentOptimizable<COpt>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>
        + Ref
        + FastTensorSum
        + FastTensorSumContractible<Sc>
        + TensorCommonFactor<Sc>
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
> ContractionStrategy<Store, L, K, FK, Aind> for SmallestDegreeIter<N, CStrat, COpt>
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <LT::WithIndices as HasStructure>::Structure: Display,
    T::Structure: Display,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    const SUPPORTS_PARTIAL_GRAPH_REWRITE: bool = true;

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
            if !product
                .contract_one_by_degree::<false, false, LT, T, L, Sc, CStrat, COpt, FK, Store>(
                    executor, graph, lib,
                )?
            {
                break;
            }
            product.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        }
        product.finish::<LT, T, L, Sc, FK, Store>(executor, graph, lib)
    }

    fn contract_product_in_place(
        executor: &mut Store,
        graph: &mut NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
        ignored: &mut SuBitGraph,
    ) -> Result<bool, TensorNetworkError<K, FK>>
    where
        K: Display + Debug,
        FK: Display + Debug,
        Aind: AbsInd,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        let mut did_progress = false;

        match product.contract_scalars_in_place::<LT, T, L, Sc, FK, Store>(
            executor, graph, operation, lib, ignored,
        )? {
            ProductRewriteProgress::NoProgress => {}
            ProductRewriteProgress::UpdatedProduct => did_progress = true,
            ProductRewriteProgress::CollapsedProduct => return Ok(true),
        }

        for _ in 0..N {
            if !product
                .contract_one_by_degree_in_place::<
                    false,
                    false,
                    LT,
                    T,
                    L,
                    Sc,
                    CStrat,
                    COpt,
                    FK,
                    Store,
                >(
                    executor, graph, lib, ignored,
                )?
            {
                break;
            }
            did_progress = true;

            match product.contract_scalars_in_place::<LT, T, L, Sc, FK, Store>(
                executor, graph, operation, lib, ignored,
            )? {
                ProductRewriteProgress::NoProgress => {}
                ProductRewriteProgress::UpdatedProduct => did_progress = true,
                ProductRewriteProgress::CollapsedProduct => return Ok(true),
            }
        }

        Ok(did_progress)
    }
}

impl<
    CStrat,
    COpt,
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<T, CStrat, LCM = T>
        + AtomComponentOptimizable<COpt>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>
        + Ref
        + FastTensorSum
        + FastTensorSumContractible<Sc>
        + TensorCommonFactor<Sc>
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
    const SCORE_ORDER: u128,
    const EXACT_JOIN_LIMIT: usize,
> ContractionStrategy<Store, L, K, FK, Aind>
    for MinResultRankWith<SCORE_ORDER, EXACT_JOIN_LIMIT, CStrat, COpt>
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <LT::WithIndices as HasStructure>::Structure: Display,
    T::Structure: Display + StructureContract,
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
            .contract_one_by_result_rank::<
                LT,
                T,
                L,
                Sc,
                CStrat,
                COpt,
                FK,
                Store,
                SCORE_ORDER,
                EXACT_JOIN_LIMIT,
            >(executor, graph, lib)?
        {
            product.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        }
        while product
            .contract_one_exterior_by_result_rank::<
                LT,
                T,
                L,
                Sc,
                CStrat,
                COpt,
                FK,
                Store,
                SCORE_ORDER,
                EXACT_JOIN_LIMIT,
            >(executor, graph, lib)?
        {
            product.contract_scalars::<LT, T, L, Sc, FK, Store>(executor, graph, lib)?;
        }
        product.finish::<LT, T, L, Sc, FK, Store>(executor, graph, lib)
    }
}

impl<
    CStrat,
    COpt,
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<T, CStrat, LCM = T>
        + AtomComponentOptimizable<COpt>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>
        + Ref
        + FastTensorSum
        + FastTensorSumContractible<Sc>
        + TensorCommonFactor<Sc>
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
> ContractionStrategy<Store, L, K, FK, Aind> for SingleSmallestDegree<D, CStrat, COpt>
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <LT::WithIndices as HasStructure>::Structure: Display,
    T::Structure: Display,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    const SUPPORTS_PARTIAL_GRAPH_REWRITE: bool = true;

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
        product.contract_one_by_degree::<D, false, LT, T, L, Sc, CStrat, COpt, FK, Store>(
            executor, graph, lib,
        )?;
        product.finish::<LT, T, L, Sc, FK, Store>(executor, graph, lib)
    }

    fn contract_product_in_place(
        executor: &mut Store,
        graph: &mut NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
        ignored: &mut SuBitGraph,
    ) -> Result<bool, TensorNetworkError<K, FK>>
    where
        K: Display + Debug,
        FK: Display + Debug,
        Aind: AbsInd,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_one_by_degree_in_place::<D, false, LT, T, L, Sc, CStrat, COpt, FK, Store>(
            executor, graph, lib, ignored,
        )
    }
}

impl<
    CStrat,
    COpt,
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<T, CStrat, LCM = T>
        + AtomComponentOptimizable<COpt>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>
        + Ref
        + FastTensorSum
        + FastTensorSumContractible<Sc>
        + TensorCommonFactor<Sc>
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
> ContractionStrategy<Store, L, K, FK, Aind> for SingleLargestDegree<D, CStrat, COpt>
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <LT::WithIndices as HasStructure>::Structure: Display,
    T::Structure: Display,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    const SUPPORTS_PARTIAL_GRAPH_REWRITE: bool = true;

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
        product.contract_one_by_degree::<D, true, LT, T, L, Sc, CStrat, COpt, FK, Store>(
            executor, graph, lib,
        )?;
        product.finish::<LT, T, L, Sc, FK, Store>(executor, graph, lib)
    }

    fn contract_product_in_place(
        executor: &mut Store,
        graph: &mut NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
        ignored: &mut SuBitGraph,
    ) -> Result<bool, TensorNetworkError<K, FK>>
    where
        K: Display + Debug,
        FK: Display + Debug,
        Aind: AbsInd,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_one_by_degree_in_place::<D, true, LT, T, L, Sc, CStrat, COpt, FK, Store>(
            executor, graph, lib, ignored,
        )
    }
}
