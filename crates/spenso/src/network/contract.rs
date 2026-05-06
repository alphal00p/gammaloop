use std::{
    fmt::{Debug, Display},
    ops::MulAssign,
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
    store::NetworkStore,
};

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
            NetworkLeaf::LibraryKey { .. } | NetworkLeaf::Scalar(_) => None,
        }
    }

    pub fn local_tensor_structure<'a, T, Sc>(
        &self,
        executor: &'a NetworkStore<T, Sc>,
        operand: usize,
    ) -> Option<&'a T::Structure>
    where
        T: HasStructure,
    {
        self.local_tensor_index(operand)
            .map(|index| executor.tensors[index].structure())
    }

    #[allow(clippy::result_large_err)]
    pub fn materialize_libraries<LT, T, L, Sc, FK>(
        &mut self,
        executor: &mut NetworkStore<T, Sc>,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
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
                let index = executor.tensors.len();
                executor.tensors.push(T::from(tensor));
                operand.leaf = NetworkLeaf::LocalTensor(index);
                operand.source = None;
            }
        }
        Ok(())
    }

    #[allow(clippy::result_large_err)]
    pub fn contract_scalars<LT, T, L, Sc, FK>(
        &mut self,
        executor: &mut NetworkStore<T, Sc>,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<bool, TensorNetworkError<K, FK>>
    where
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
                        *accumulator *= executor.scalar[index].refer();
                    } else {
                        accumulator = Some(executor.scalar[index].clone());
                    }
                    true
                }
                NetworkLeaf::LocalTensor(index) => {
                    if let Some(scalar) = executor.tensors[index].scalar_ref() {
                        if let Some(accumulator) = &mut accumulator {
                            *accumulator *= scalar;
                        } else {
                            accumulator =
                                Some(Sc::from(executor.tensors[index].clone().scalar().unwrap()));
                        }
                        true
                    } else {
                        false
                    }
                }
                NetworkLeaf::LibraryKey { .. } => false,
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
                let index = executor.scalar.len();
                executor.scalar.push(accumulator);
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
                let leaf = self.scalar_multiply_operand(
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
                let index = executor.scalar.len();
                executor.scalar.push(accumulator);
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
    fn scalar_multiply_operand<LT, T, L, Sc, FK>(
        &self,
        operand: usize,
        scalar: Sc,
        executor: &mut NetworkStore<T, Sc>,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
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
            NetworkLeaf::LocalTensor(index) => executor.tensors[*index]
                .scalar_mul(&scalar)
                .ok_or(TensorNetworkError::FailedScalarMul)?,
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

        let index = executor.tensors.len();
        executor.tensors.push(leaf);
        Ok(NetworkLeaf::LocalTensor(index))
    }

    #[allow(clippy::result_large_err)]
    pub fn contract_pair<LT, T, L, Sc, CStrat, FK>(
        &mut self,
        left: usize,
        right: usize,
        executor: &mut NetworkStore<T, Sc>,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Clone + Debug + Display,
        FK: Debug + Display,
        Aind: AbsInd,
        LT: LibraryTensor + Clone,
        T: HasStructure + From<LT::WithIndices> + Contract<T, CStrat, LCM = T>,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        self.materialize_libraries::<LT, T, L, Sc, FK>(executor, graph, lib)?;

        let left_tensor = self
            .local_tensor_index(left)
            .ok_or(TensorNetworkError::SlotEdgeToScalarNode)?;
        let right_tensor = self
            .local_tensor_index(right)
            .ok_or(TensorNetworkError::SlotEdgeToScalarNode)?;
        let contracted = executor.tensors[left_tensor].contract(&executor.tensors[right_tensor])?;
        let index = executor.tensors.len();
        executor.tensors.push(contracted);

        let mut positions = [left, right];
        positions.sort_unstable();
        self.replace_operands(
            &positions,
            ProductOperand {
                leaf: NetworkLeaf::LocalTensor(index),
                source: None,
            },
        );
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
    >(
        &mut self,
        executor: &mut NetworkStore<T, Sc>,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<bool, TensorNetworkError<K, FK>>
    where
        K: Clone + Debug + Display,
        FK: Debug + Display,
        Aind: AbsInd,
        LT: LibraryTensor + Clone,
        T: HasStructure + From<LT::WithIndices> + Contract<T, CStrat, LCM = T>,
        T::Structure: Display,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
        LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        self.materialize_libraries::<LT, T, L, Sc, FK>(executor, graph, lib)?;
        let Some((left, right, _degree)) = self.best_degree_pair::<T, Sc, LARGEST>(executor) else {
            return Ok(false);
        };

        if DEBUG {
            let left_tensor = self.local_tensor_index(left).unwrap();
            let right_tensor = self.local_tensor_index(right).unwrap();
            println!(
                "Contracting {} with {}",
                executor.tensors[left_tensor].structure(),
                executor.tensors[right_tensor].structure()
            );
        }

        self.contract_pair::<LT, T, L, Sc, CStrat, FK>(left, right, executor, graph, lib)?;

        if DEBUG && let NetworkLeaf::LocalTensor(index) = self.operands[left.min(right)].leaf {
            println!("Obtained {}", executor.tensors[index].structure());
        }

        Ok(true)
    }

    pub fn best_degree_pair<T, Sc, const LARGEST: bool>(
        &self,
        executor: &NetworkStore<T, Sc>,
    ) -> Option<(usize, usize, u32)>
    where
        T: HasStructure,
    {
        self.best_tensor_pair_by::<T, Sc, _, LARGEST>(executor, |_, _, degree, _, _| {
            if degree == 0 { u32::MAX } else { degree }
        })
        .map(|(left, right, degree, _)| (left, right, degree))
    }

    pub fn best_tensor_pair_by<T, Sc, Score: Ord, const LARGEST: bool>(
        &self,
        executor: &NetworkStore<T, Sc>,
        mut score: impl FnMut(usize, usize, u32, &T, &T) -> Score,
    ) -> Option<(usize, usize, u32, Score)>
    where
        T: HasStructure,
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
                let left = &executor.tensors[*left_tensor];
                let right = &executor.tensors[*right_tensor];
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

    #[allow(clippy::result_large_err)]
    pub fn finish<LT, T, L, Sc, FK>(
        mut self,
        executor: &mut NetworkStore<T, Sc>,
        graph: &NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
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
        while self.contract_scalars::<LT, T, L, Sc, FK>(executor, graph, lib)? {}

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
        + From<LT::WithIndices>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
> ContractionStrategy<NetworkStore<T, Sc>, L, K, FK, Aind> for ContractScalars<CStrat>
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    fn contract(
        executor: &mut NetworkStore<T, Sc>,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_scalars::<LT, T, L, Sc, FK>(executor, graph, lib)?;
        product.finish::<LT, T, L, Sc, FK>(executor, graph, lib)
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
        + From<LT::WithIndices>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
> ContractionStrategy<NetworkStore<T, Sc>, L, K, FK, Aind> for SmallestDegree<CStrat>
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
        executor: &mut NetworkStore<T, Sc>,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_scalars::<LT, T, L, Sc, FK>(executor, graph, lib)?;
        while product.contract_one_by_degree::<false, false, LT, T, L, Sc, CStrat, FK>(
            executor, graph, lib,
        )? {
            product.contract_scalars::<LT, T, L, Sc, FK>(executor, graph, lib)?;
        }
        product.finish::<LT, T, L, Sc, FK>(executor, graph, lib)
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
        + From<LT::WithIndices>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
    const N: usize,
> ContractionStrategy<NetworkStore<T, Sc>, L, K, FK, Aind> for SmallestDegreeIter<N, CStrat>
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
        executor: &mut NetworkStore<T, Sc>,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_scalars::<LT, T, L, Sc, FK>(executor, graph, lib)?;
        for _ in 0..N {
            if !product.contract_one_by_degree::<false, false, LT, T, L, Sc, CStrat, FK>(
                executor, graph, lib,
            )? {
                break;
            }
            product.contract_scalars::<LT, T, L, Sc, FK>(executor, graph, lib)?;
        }
        product.finish::<LT, T, L, Sc, FK>(executor, graph, lib)
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
        + From<LT::WithIndices>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
    const D: bool,
> ContractionStrategy<NetworkStore<T, Sc>, L, K, FK, Aind> for SingleSmallestDegree<D, CStrat>
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
        executor: &mut NetworkStore<T, Sc>,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_scalars::<LT, T, L, Sc, FK>(executor, graph, lib)?;
        product
            .contract_one_by_degree::<D, false, LT, T, L, Sc, CStrat, FK>(executor, graph, lib)?;
        product.finish::<LT, T, L, Sc, FK>(executor, graph, lib)
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
        + From<LT::WithIndices>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
    const D: bool,
> ContractionStrategy<NetworkStore<T, Sc>, L, K, FK, Aind> for SingleLargestDegree<D, CStrat>
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
        executor: &mut NetworkStore<T, Sc>,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_scalars::<LT, T, L, Sc, FK>(executor, graph, lib)?;
        product
            .contract_one_by_degree::<D, true, LT, T, L, Sc, CStrat, FK>(executor, graph, lib)?;
        product.finish::<LT, T, L, Sc, FK>(executor, graph, lib)
    }
}
