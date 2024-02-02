use indexmap::IndexMap;

use intmap::IntMap;
use itertools::Itertools;

use petgraph::{
    dot::{Config, Dot},
    graph::EdgeIndex,
    graph::Graph,
    graph::NodeIndex,
    visit::EdgeRef,
    Undirected,
};
use symbolica::{
    representations::{Identifier, Num},
    state::{State, Workspace},
};

use crate::HasIntegrand;
use nohash_hasher::BuildNoHashHasher;
type NHIndexMap<K, V> = IndexMap<K, V, BuildNoHashHasher<K>>;

use self::mixed_tensor::{MixedTensor, MixedTensors, SymbolicContract};

use super::*;
use smartstring::alias::String;
use std::{
    fmt::{Debug, Display},
    ops::Neg,
};
impl<T, I> DenseTensor<T, I>
where
    T: for<'a> std::ops::AddAssign<&'a T>
        + for<'b> std::ops::SubAssign<&'b T>
        + Neg<Output = T>
        + Clone
        + std::fmt::Debug,
    I: Clone,
{
    pub fn internal_contract(&self) -> Self {
        let mut result: DenseTensor<T, I> = self.clone();
        for trace in self.traces() {
            let mut new_structure = self.structure().clone();
            new_structure.trace(trace[0], trace[1]);

            let mut new_result = DenseTensor::from_data_coerced(&self.data, new_structure).unwrap();
            for (idx, t) in result.iter_trace(trace) {
                new_result.set(&idx, t);
            }
            result = new_result;
        }
        result
    }
}

impl<T, I> SparseTensor<T, I>
where
    T: for<'a> std::ops::AddAssign<&'a T>
        + for<'b> std::ops::SubAssign<&'b T>
        + std::ops::Neg<Output = T>
        + Clone
        + Default
        + PartialEq,
    I: Clone,
{
    pub fn internal_contract(&self) -> Self {
        let trace = self.traces()[0];

        // println!("trace {:?}", trace);
        let mut new_structure = self.structure().clone();
        new_structure.trace(trace[0], trace[1]);

        let mut new_result = SparseTensor::empty(new_structure);
        for (idx, t) in self.iter_trace(trace).filter(|(_, t)| *t != T::default()) {
            new_result.set(&idx, t).unwrap();
        }

        if new_result.traces().is_empty() {
            new_result
        } else {
            new_result.internal_contract()
        }
    }
}

// pub trait NumTensor {}

// impl<T> NumTensor for DenseTensor<T> {}
// impl<T> NumTensor for SparseTensor<T> {}
// pub trait Densest<A: NumTensor, B: NumTensor> {}
// pub trait Contractable<T> {
//     fn contract<C: NumTensor + Densest<SparseTensor<T>, Self>>(
//         &self,
//         other: &SparseTensor<T>,
//     ) -> Option<C>;
// }

// pub fn mul<'a, 'b, T, U>(right: &'a T, left: &'b U) -> <&'b U as SmallestUpgrade<&'b T>>::LCM
// where
//     &'b U: SmallestUpgrade<&'b T>,
//     &'a T: SmallestUpgrade<&'b U, LCM = <&'b U as SmallestUpgrade<&'b T>>::LCM>,
//     <&'b U as SmallestUpgrade<&'b T>>::LCM: std::ops::Mul<
//         <&'b U as SmallestUpgrade<&'b T>>::LCM,
//         Output = <&'b U as SmallestUpgrade<&'b T>>::LCM,
//     >,
// {
//     left.upgrade() * right.upgrade()
// }

// pub fn mul<T, U>(right: T, left: U) -> U::LCM
// where
//     U: SmallestUpgrade<T>,
//     T: SmallestUpgrade<U>,
//     U::LCM: std::ops::Mul<T::LCM, Output = U::LCM>,
// {
//     left.upgrade() * right.upgrade()
// }
pub trait Contract<T> {
    type LCM;
    fn contract(&self, other: &T) -> Option<Self::LCM>;
}

impl<T, U, Out, I> Contract<DenseTensor<T, I>> for DenseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = Out>,
    Out: std::ops::AddAssign<Out>
        + std::ops::SubAssign<Out>
        + Neg<Output = Out>
        + Clone
        + num::traits::Zero
        + for<'a> std::ops::AddAssign<&'a Out>
        + for<'b> std::ops::SubAssign<&'b Out>
        + std::fmt::Debug,
    I: Clone,
{
    type LCM = DenseTensor<Out, I>;
    fn contract(&self, other: &DenseTensor<T, I>) -> Option<Self::LCM> {
        if let Some((i, j)) = self.structure().match_index(other.structure()) {
            // println!("{},{}", i, j);
            let self_shape = self.shape();

            let dimension_of_contraction = self_shape[i];
            let metric = self.get_ith_metric(i).unwrap();

            let final_structure = self.structure().merge_at(other.structure(), (i, j));

            // Initialize result tensor with default values
            let mut result_data = vec![Out::zero(); final_structure.size()];
            let mut result_index = 0;

            for fiber_a in self.iter_fibers(i) {
                for fiber_b in other.iter_fibers(j) {
                    // final_structure
                    //     .flat_index(
                    //         &index_a[..i]
                    //             .iter()
                    //             .chain(&index_a[i + 1..])
                    //             .chain(&index_b[..j])
                    //             .chain(&index_b[j + 1..])
                    //             .cloned()
                    //             .collect::<Vec<_>>(),
                    //     )
                    //     .unwrap();

                    for k in 0..dimension_of_contraction {
                        // Adjust indices for fetching from the other tensor
                        if metric[k] {
                            result_data[result_index] -= fiber_a[k] * fiber_b[k];
                        } else {
                            result_data[result_index] += fiber_a[k] * fiber_b[k];
                        }
                    }
                    result_index += 1;
                }
            }

            let result: DenseTensor<Out, I> = DenseTensor {
                data: result_data,
                structure: final_structure,
            };

            if result.traces().is_empty() {
                return Some(result);
            } else {
                return Some(result.internal_contract());
            }
        }
        None
    }
}

impl<T, U, I, Out> Contract<DenseTensor<T, I>> for SparseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = Out>,
    Out: std::ops::AddAssign<Out>
        + std::ops::SubAssign<Out>
        + Neg<Output = Out>
        + Clone
        + num::traits::Zero
        + for<'a> std::ops::AddAssign<&'a Out>
        + for<'b> std::ops::SubAssign<&'b Out>
        + std::fmt::Debug,
    I: Clone,
    U: Clone,
{
    type LCM = DenseTensor<Out, I>;
    fn contract(&self, other: &DenseTensor<T, I>) -> Option<Self::LCM> {
        // println!("Contracting SparseTensor DenseTensor");
        if let Some((i, j)) = self.structure().match_index(other.structure()) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = vec![Out::zero(); final_structure.size()];

            let metric = self.get_ith_metric(i).unwrap();
            let mut result_index = 0;

            for (skipped, nonzeros, fiber_a) in self.iter_fibers(i) {
                result_index += skipped;
                for fiber_b in other.iter_fibers(j) {
                    for (i, k) in nonzeros.iter().enumerate() {
                        // Adjust indices for fetching from the other tensor
                        if metric[*k] {
                            result_data[result_index] -= fiber_a[i] * fiber_b[*k];
                        } else {
                            result_data[result_index] += fiber_a[i] * fiber_b[*k];
                        }
                    }
                    result_index += 1;
                }
            }

            let result = DenseTensor {
                data: result_data,
                structure: final_structure,
            };

            if result.traces().is_empty() {
                return Some(result);
            } else {
                return Some(result.internal_contract());
            }
        }
        None
    }
}

impl<T, U, Out, I> Contract<SparseTensor<T, I>> for SparseTensor<U, I>
where
    // T: SmallestUpgrade<U> + Copy,
    // U: SmallestUpgrade<T, LCM = T::LCM> + Copy,
    // T::LCM: std::ops::AddAssign<T::LCM>
    //     + std::ops::SubAssign<T::LCM>
    //     + for<'a> std::ops::AddAssign<&'a T::LCM>
    //     + for<'b> std::ops::SubAssign<&'b T::LCM>
    //     + std::fmt::Debug
    //     + Neg<Output = T::LCM>
    //     + std::cmp::PartialOrd
    //     + Default
    //     + Clone
    //     + std::ops::Mul<T::LCM, Output = T::LCM>,
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = Out>,
    Out: std::ops::AddAssign<Out>
        + std::ops::SubAssign<Out>
        + Neg<Output = Out>
        + Clone
        + Default
        + num::traits::Zero
        + for<'a> std::ops::AddAssign<&'a Out>
        + for<'b> std::ops::SubAssign<&'b Out>
        + std::fmt::Debug
        + std::cmp::PartialEq,
    I: Clone,
    U: Clone,
    T: Clone,
{
    type LCM = SparseTensor<Out, I>;
    fn contract(&self, other: &SparseTensor<T, I>) -> Option<Self::LCM> {
        // println!("Contracting SparseTensor SparseTensor");
        if let Some((i, j)) = self.structure().match_index(other.structure()) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = IntMap::default();
            let one = 1;
            let stride_other = *final_structure
                .strides()
                .get(self.structure().order() - 2)
                .unwrap_or(&one);

            let metric = self.get_ith_metric(i).unwrap();
            let mut result_index = 0;

            for (skipped_a, nonzeros_a, fiber_a) in self.iter_fibers(i) {
                result_index += skipped_a;
                for (skipped_b, nonzeros_b, fiber_b) in other.iter_fibers(j) {
                    result_index += skipped_b * stride_other;
                    let mut value = Out::zero();
                    let mut nonzero = false;
                    for (i, j, x) in nonzeros_a.iter().enumerate().filter_map(|(i, &x)| {
                        nonzeros_b.binary_search(&x).ok().map(|j| (i, j, x)) // Only store the positions
                    }) {
                        // Adjust indices for fetching from the other tensor
                        if metric[x] {
                            value -= fiber_a[i] * fiber_b[j];
                        } else {
                            value += fiber_a[i] * fiber_b[j];
                        }

                        nonzero = true;
                    }

                    if nonzero && value != Out::zero() {
                        result_data.insert(result_index as u64, value);
                    }
                    result_index += 1;
                }
            }

            let result = SparseTensor {
                elements: result_data,
                structure: final_structure,
            };

            if result.traces().is_empty() {
                return Some(result);
            } else {
                return Some(result.internal_contract());
            }
        }
        None
    }
}

impl<T, U, Out, I> Contract<SparseTensor<T, I>> for DenseTensor<U, I>
where
    // T: SmallestUpgrade<U> + Copy,
    // U: SmallestUpgrade<T, LCM = T::LCM> + Copy,
    // T::LCM: std::ops::AddAssign<T::LCM>
    //     + std::ops::SubAssign<T::LCM>
    //     + for<'a> std::ops::AddAssign<&'a T::LCM>
    //     + for<'b> std::ops::SubAssign<&'b T::LCM>
    //     + std::fmt::Debug
    //     + Neg<Output = T::LCM>
    //     + Default
    //     + Clone
    //     + std::ops::Mul<T::LCM, Output = T::LCM>,
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = Out>,
    Out: std::ops::AddAssign<Out>
        + std::ops::SubAssign<Out>
        + Neg<Output = Out>
        + Clone
        + num::traits::Zero
        + for<'a> std::ops::AddAssign<&'a Out>
        + for<'b> std::ops::SubAssign<&'b Out>
        + std::fmt::Debug,
    I: Clone,
    T: Clone,
{
    type LCM = DenseTensor<Out, I>;
    fn contract(&self, other: &SparseTensor<T, I>) -> Option<Self::LCM> {
        // self is dense U, other is sparse T
        if let Some((i, j)) = self.structure().match_index(other.structure()) {
            // final structure is the structure of self appended with the structure of other, with the contracted indices removed
            let final_structure = self.structure().merge_at(other.structure(), (i, j));

            let one = 1;
            let stride_other = *final_structure
                .strides()
                .get(self.structure().order() - 2)
                .unwrap_or(&one);
            let mut result_data = vec![Out::zero(); final_structure.size()];

            let metric = other.get_ith_metric(j).unwrap();
            let mut result_index = 0;

            for fiber_a in self.iter_fibers(i) {
                for (skipped, nonzeros, fiber_b) in other.iter_fibers(j) {
                    //nonzeros is the indices of the non-zero elements of the other tensor
                    result_index += skipped * stride_other;

                    for (i, k) in nonzeros.iter().enumerate() {
                        // Adjust indices for fetching from the other tensor
                        if metric[*k] {
                            result_data[result_index] -= fiber_a[*k] * fiber_b[i];
                        } else {
                            result_data[result_index] += fiber_a[*k] * fiber_b[i];
                        }
                    }
                    result_index += 1;
                }
            }

            let result = DenseTensor {
                data: result_data,
                structure: final_structure,
            };

            if result.traces().is_empty() {
                return Some(result);
            } else {
                return Some(result.internal_contract());
            }
        }
        None
    }
}

impl<T, U, I, Out> Contract<NumTensor<T, I>> for NumTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = Out>,
    Out: std::ops::AddAssign<Out>
        + std::ops::SubAssign<Out>
        + Neg<Output = Out>
        + Clone
        + Default
        + num::traits::Zero
        + for<'a> std::ops::AddAssign<&'a Out>
        + for<'b> std::ops::SubAssign<&'b Out>
        + std::fmt::Debug
        + std::cmp::PartialEq,
    I: Clone,
    T: Clone,
    U: Clone,
{
    type LCM = NumTensor<Out, I>;
    fn contract(&self, other: &NumTensor<T, I>) -> Option<NumTensor<Out, I>> {
        match (self, other) {
            (NumTensor::Dense(s), NumTensor::Dense(o)) => Some(NumTensor::Dense(s.contract(o)?)),
            (NumTensor::Dense(s), NumTensor::Sparse(o)) => Some(NumTensor::Dense(s.contract(o)?)),
            (NumTensor::Sparse(s), NumTensor::Dense(o)) => Some(NumTensor::Dense(s.contract(o)?)),
            (NumTensor::Sparse(s), NumTensor::Sparse(o)) => Some(NumTensor::Sparse(s.contract(o)?)),
        }
    }
}

impl<I> Contract<NumTensors<I>> for NumTensors<I>
where
    I: Clone,
{
    type LCM = NumTensors<I>;
    fn contract(&self, other: &NumTensors<I>) -> Option<Self::LCM> {
        match (self, other) {
            (NumTensors::Float(a), NumTensors::Float(b)) => Some(NumTensors::Float(a.contract(b)?)),
            (NumTensors::Float(a), NumTensors::Complex(b)) => {
                Some(NumTensors::Complex(a.contract(b)?))
            }
            (NumTensors::Complex(a), NumTensors::Float(b)) => {
                Some(NumTensors::Complex(a.contract(b)?))
            }
            (NumTensors::Complex(a), NumTensors::Complex(b)) => {
                Some(NumTensors::Complex(a.contract(b)?))
            }
        }
    }
}

// trait ContractableTensorInt: HasTensorStructure {}

// impl<T> ContractableTensorInt for DenseTensor<T> {}
// impl<T> ContractableTensorInt for SparseTensor<T> {}

// trait ContractableTensor: ContractableTensorInt {
//     fn contract<T: Contract<Self>>(&self, other: &T) -> Option<T::LCM> {
//         other.contract(self)
//     }
// }

// impl<T,U> ContractableTensor for DenseTensor<T> {
//     fn contract<DenseTensor<U>>(&self, other: &DenseTensor<U>) -> Option<Self> {
//         self.contract(other)
//     }
// }

// enum UpgradableType {
//     F64(f64),
//     Complex(Complex<f64>),
// }

// impl Mul for UpgradableType {
//     type Output = Self;
//     fn mul(self, rhs: Self) -> Self::Output {
//         match (self, rhs) {
//             (UpgradableType::F64(a), UpgradableType::F64(b)) => UpgradableType::F64(a * b),
//             (UpgradableType::Complex(a), UpgradableType::Complex(b)) => {
//                 UpgradableType::Complex(a * b)
//             }
//             (UpgradableType::F64(a), UpgradableType::Complex(b)) => {
//                 UpgradableType::Complex(Complex::new(a, 0.0) * b)
//             }
//             (UpgradableType::Complex(a), UpgradableType::F64(b)) => {
//                 UpgradableType::Complex(a * Complex::new(b, 0.0))
//             }
//         }
//     }
// }

// impl From<f64> for UpgradableType {
//     fn from(other: f64) -> Self {
//         UpgradableType::F64(other)
//     }
// }

// impl From<Complex<f64>> for UpgradableType {
//     fn from(other: Complex<f64>) -> Self {
//         UpgradableType::Complex(other)
//     }
// }

// #[enum_dispatch(HasTensorStructure)]
// enum ContractableTensor {
//     Dense(DenseTensor<UpgradableType>),
//     Sparse(SparseTensor<UpgradableType>),
// }

// impl ContractableTensor {
//     fn contract(&self, other: &Self) -> Option<Self> {
//         match (self, other) {
//             (ContractableTensor::Dense(a), ContractableTensor::Dense(b)) => {
//                 a.contract(b).map(ContractableTensor::Dense)
//             }
//             (ContractableTensor::Dense(a), ContractableTensor::Sparse(b)) => {
//                 a.contract(b).map(ContractableTensor::Dense)
//             }
//             (ContractableTensor::Sparse(a), ContractableTensor::Dense(b)) => {
//                 a.contract(b).map(ContractableTensor::Dense)
//             }
//             (ContractableTensor::Sparse(a), ContractableTensor::Sparse(b)) => {
//                 a.contract(b).map(ContractableTensor::Dense)
//             }
//         }
//     }
// }

// // struct TensorNetwork {
// //     tensors: Vec<ContractableTensor>,
// // }
#[derive(Debug, Clone)]
pub struct TensorNetwork<T> {
    graph: Graph<T, Slot, Undirected>,
}

impl<T> TensorNetwork<T> {
    fn edge_to_min_degree_node(&self) -> Option<EdgeIndex> {
        let mut max_degree = 0;
        let mut min_degree = usize::MAX;
        let mut edge_to_min_degree_node = None;
        for node in self.graph.node_indices() {
            if self.graph.edges(node).count() > max_degree {
                max_degree = self.graph.edges(node).count();
            }
            if self.graph.edges(node).count() < min_degree {
                min_degree = self.graph.edges(node).count();
                if min_degree > 0 {
                    edge_to_min_degree_node = Some(self.graph.edges(node).next().unwrap().id());
                }
            }
        }
        edge_to_min_degree_node
    }

    pub fn to_vec(self) -> Vec<T> {
        self.graph
            .into_nodes_edges()
            .0
            .into_iter()
            .map(|n| n.weight)
            .collect()
    }
}

impl<N> TensorNetwork<MixedTensor<N>>
where
    N: Debug,
{
    pub fn to_symbolic_tensor_vec(self) -> Vec<NumTensor<Atom, N>> {
        self.graph
            .into_nodes_edges()
            .0
            .into_iter()
            .filter(|n| n.weight.is_symbolic())
            .map(|n| n.weight.try_into_symbolic().unwrap())
            .collect()
    }
}
impl<T> TensorNetwork<T>
where
    T: HasTensorStructure,
{
    pub fn new(tensors: Vec<T>) -> Self {
        TensorNetwork {
            graph: Self::generate_network_graph(tensors),
        }
    }

    fn generate_network_graph(tensors: Vec<T>) -> Graph<T, Slot, Undirected> {
        let mut graph = Graph::<T, Slot, Undirected>::new_undirected();

        for tensor in tensors {
            graph.add_node(tensor);
        }

        for (ni, n) in graph.node_indices().enumerate() {
            for (_, m) in graph.node_indices().enumerate().filter(|(i, _)| *i > ni) {
                let a = graph.node_weight(n).unwrap();
                let b = graph.node_weight(m).unwrap();

                if let Some((i, _)) = a.structure().match_index(b.structure()) {
                    graph.add_edge(n, m, a.structure()[i]);
                }
            }
        }
        graph
    }

    fn merge_nodes(&mut self, a: NodeIndex, b: NodeIndex, weight: T) {
        let neighsa = self.graph.neighbors(a).count();
        let neighsb = self.graph.neighbors(b).count();

        let (m, d) = if neighsa < neighsb { (b, a) } else { (a, b) };

        let mut new_edges = vec![];

        for e in self.graph.edges(d) {
            let (a, b) = (e.target(), e.source());
            let n = if a == d { b } else { a };
            if n == m || n == d {
                continue;
            }
            new_edges.push((m, n, *e.weight()));
        }

        for (a, b, w) in new_edges {
            self.graph.add_edge(a, b, w);
        }

        if let Some(w) = self.graph.node_weight_mut(m) {
            *w = weight;
        }

        self.graph.remove_node(d);
    }

    pub fn edge_to_min_degree_node_with_depth(&self, depth: usize) -> Option<EdgeIndex> {
        let mut min_degree = usize::MAX;
        let mut edge_to_min_degree_node = None;
        for edge in self
            .graph
            .edge_indices()
            .filter(|&e| {
                let (nl, nr) = self.graph.edge_endpoints(e).unwrap();
                self.graph.node_weight(nl).unwrap().contractions_num()
                    + self.graph.node_weight(nr).unwrap().contractions_num()
                    < depth
            })
            .sorted_by(|&e1, &e2| {
                let (nl1, nr1) = self.graph.edge_endpoints(e1).unwrap();
                let (nl2, nr2) = self.graph.edge_endpoints(e2).unwrap();
                (self.graph.node_weight(nl1).unwrap().contractions_num()
                    + self.graph.node_weight(nr1).unwrap().contractions_num())
                .cmp(
                    &(self.graph.node_weight(nl2).unwrap().contractions_num()
                        + self.graph.node_weight(nr2).unwrap().contractions_num()),
                )
            })
        {
            let (nl, nr) = self.graph.edge_endpoints(edge).unwrap();
            let degrees = self.graph.edges(nl).count() + self.graph.edges(nr).count();
            if degrees < min_degree {
                min_degree = degrees;

                if min_degree > 0 {
                    edge_to_min_degree_node = Some(edge);
                }
            }
        }
        edge_to_min_degree_node
    }
}
impl<T> TensorNetwork<T>
where
    T: Clone,
{
    pub fn result(&self) -> T {
        self.graph
            .node_weight(self.graph.node_indices().next().unwrap())
            .unwrap()
            .clone()
    }
}

impl<T> TensorNetwork<T>
where
    T: Debug + HasTensorStructure<Name = String>,
{
    pub fn dot(&self) -> String {
        format!(
            "{:?}",
            Dot::with_attr_getters(
                &self.graph,
                &[Config::EdgeNoLabel, Config::NodeNoLabel],
                &|_, e| { format!("label=\"{}\"", e.weight()).into() },
                &|_, n| { format!("label=\"{}\"", n.1.structure()).into() }
            )
        )
        .into()
    }
}

impl<T> TensorNetwork<T>
where
    T: Debug + HasTensorStructure<Name = Identifier>,
{
    pub fn dotsym(&self, state: &State) -> String {
        format!(
            "{:?}",
            Dot::with_attr_getters(
                &self.graph,
                &[Config::EdgeNoLabel, Config::NodeNoLabel],
                &|_, e| { format!("label=\"{}\"", e.weight()).into() },
                &|_, n| { format!("label=\"{}\"", n.1.structure().to_string(state)).into() }
            )
        )
        .into()
    }
}

impl<T> TensorNetwork<T>
where
    T: HasTensorStructure<Name = Identifier> + Clone,
{
    pub fn symbolic_shadow(
        &mut self,
        name: &str,
        state: &mut State,
        ws: &Workspace,
    ) -> TensorNetwork<MixedTensors> {
        self.graph.node_indices().for_each(|n| {
            self.graph
                .node_weight_mut(n)
                .unwrap()
                .mut_structure()
                .set_global_name(
                    state
                        .get_or_insert_fn(format!("{}{}", name, n.index()), None)
                        .unwrap(),
                );
        });
        let g: Graph<MixedTensors, Slot, Undirected> = Graph::map(
            &self.graph,
            |_, nw| {
                MixedTensor::<Identifier>::from(nw.structure().clone().to_dense(state, ws).unwrap())
            },
            |_, &w| w,
        );
        TensorNetwork { graph: g }
    }
}

impl<T> TensorNetwork<T>
where
    T: HasTensorStructure<Name = String>,
{
    pub fn name(&mut self, name: T::Name) {
        self.graph.node_indices().for_each(|n| {
            self.graph
                .node_weight_mut(n)
                .unwrap()
                .mut_structure()
                .set_global_name(format!("{}{}", name, n.index()).into());
        });
    }
}

impl<T> TensorNetwork<T>
where
    T: HasTensorStructure<Name = Identifier>,
{
    pub fn namesym(&mut self, name: &str, state: &mut State) {
        self.graph.node_indices().for_each(|n| {
            self.graph
                .node_weight_mut(n)
                .unwrap()
                .mut_structure()
                .set_global_name(
                    state
                        .get_or_insert_fn(format!("{}{}", name, n.index()), None)
                        .unwrap(),
                );
        });
    }
}

impl<T> TensorNetwork<T>
where
    T: Contract<T, LCM = T>,
    T: HasTensorStructure,
{
    fn contract_edge(&mut self, edge_idx: EdgeIndex) {
        let (a, b) = self.graph.edge_endpoints(edge_idx).unwrap();

        let ai = self.graph.node_weight(a).unwrap();
        let bi = self.graph.node_weight(b).unwrap();

        let f = ai.contract(bi).unwrap();

        self.merge_nodes(a, b, f);
    }
    pub fn contract(&mut self) {
        if let Some(e) = self.edge_to_min_degree_node() {
            self.contract_edge(e);
            self.contract();
        }
    }
}

impl<T> TensorNetwork<T>
where
    T: HasTensorStructure + SymbolicContract<T, LCM = T> + Debug,
{
    fn contract_edge_sym(&mut self, edge_idx: EdgeIndex, state: &State, ws: &Workspace) {
        let (a, b) = self.graph.edge_endpoints(edge_idx).unwrap();

        let ai = self.graph.node_weight(a).unwrap();
        let bi = self.graph.node_weight(b).unwrap();

        let f = ai.contract_sym(bi, state, ws).unwrap();

        self.merge_nodes(a, b, f);
    }
    pub fn contract_sym(&mut self, state: &State, ws: &Workspace) {
        if let Some(e) = self.edge_to_min_degree_node() {
            self.contract_edge_sym(e, state, ws);
            self.contract_sym(state, ws)
        }
    }
    pub fn contract_sym_depth(&mut self, depth: usize, state: &State, ws: &Workspace) {
        if let Some(e) = self.edge_to_min_degree_node_with_depth(depth) {
            self.contract_edge_sym(e, state, ws);
            self.contract_sym_depth(depth, state, ws)
        }
    }
}
