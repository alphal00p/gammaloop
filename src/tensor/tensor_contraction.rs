use ahash::AHashMap;
use enum_dispatch::enum_dispatch;
use hyperdual::Zero;
use num::Complex;
use petgraph::{
    dot::{Config, Dot},
    graph::Edge,
    graph::EdgeIndex,
    graph::NodeIndex,
    visit::{EdgeRef, IntoEdges},
    Graph, Undirected,
};

use super::*;
use std::{
    borrow::Borrow,
    ops::{Mul, Neg},
    process::Output,
};
impl<T> DenseTensor<T>
where
    T: for<'a> std::ops::AddAssign<&'a T>
        + for<'b> std::ops::SubAssign<&'b T>
        + Neg<Output = T>
        + Clone
        + std::fmt::Debug,
{
    pub fn internal_contract(&self) -> Self {
        let mut result = self.clone();
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

impl<T> SparseTensor<T>
where
    T: for<'a> std::ops::AddAssign<&'a T>
        + for<'b> std::ops::SubAssign<&'b T>
        + std::ops::Neg<Output = T>
        + Clone
        + Default
        + PartialEq,
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
pub trait Contract<T>: HasTensorStructure {
    type LCM;
    fn contract(&self, other: &T) -> Option<Self::LCM>;
}

impl<T, U, Out> Contract<DenseTensor<T>> for DenseTensor<U>
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
    // T::LCM: std::ops::AddAssign<T::LCM>
    //     + std::ops::SubAssign<T::LCM>
    //     + for<'a> std::ops::AddAssign<&'a T::LCM>
    //     + for<'b> std::ops::SubAssign<&'b T::LCM>
    //     + std::fmt::Debug
    //     + Neg<Output = T::LCM>
    //     + Default
    //     + Clone
    //     + std::ops::Mul<T::LCM, Output = T::LCM>,
{
    type LCM = DenseTensor<Out>;
    fn contract(&self, other: &DenseTensor<T>) -> Option<Self::LCM> {
        if let Some((i, j)) = self.match_index(other) {
            // println!("{},{}", i, j);
            let self_shape = self.shape();

            let dimension_of_contraction = self_shape[i];
            let metric = self.get_ith_metric(i).unwrap();

            let final_structure = self.structure().merge_at(other.structure(), (i, j));

            // Initialize result tensor with default values
            let mut result_data = vec![Out::zero(); final_structure.size()];

            for (index_a, fiber_a) in self.iter_fibers(i) {
                for (index_b, fiber_b) in other.iter_fibers(j) {
                    let result_index = final_structure
                        .flat_index(
                            &index_a[..i]
                                .iter()
                                .chain(&index_a[i + 1..])
                                .chain(&index_b[..j])
                                .chain(&index_b[j + 1..])
                                .cloned()
                                .collect::<Vec<_>>(),
                        )
                        .unwrap();

                    for k in 0..dimension_of_contraction {
                        // Adjust indices for fetching from the other tensor
                        if metric[k] {
                            result_data[result_index] -= fiber_a[k] * fiber_b[k];
                        } else {
                            result_data[result_index] += fiber_a[k] * fiber_b[k];
                        }
                    }
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

impl<T, U, Out> Contract<DenseTensor<T>> for SparseTensor<U>
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
{
    type LCM = DenseTensor<Out>;
    fn contract(&self, other: &DenseTensor<T>) -> Option<Self::LCM> {
        // println!("Contracting SparseTensor DenseTensor");
        if let Some((i, j)) = self.match_index(other) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = vec![Out::zero(); final_structure.size()];

            let metric = self.get_ith_metric(i).unwrap();

            for (index_a, nonzeros, fiber_a) in self.iter_fibers(i) {
                for (index_b, fiber_b) in other.iter_fibers(j) {
                    let result_index = final_structure
                        .flat_index(
                            &index_a[..i]
                                .iter()
                                .chain(&index_a[i + 1..])
                                .chain(&index_b[..j])
                                .chain(&index_b[j + 1..])
                                .cloned()
                                .collect::<Vec<_>>(),
                        )
                        .unwrap();
                    for (i, k) in nonzeros.iter().enumerate() {
                        // Adjust indices for fetching from the other tensor
                        if metric[*k] {
                            result_data[result_index] -= fiber_a[i] * fiber_b[*k];
                        } else {
                            result_data[result_index] += fiber_a[i] * fiber_b[*k];
                        }
                    }
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

impl<T, U, Out> Contract<SparseTensor<T>> for SparseTensor<U>
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
{
    type LCM = SparseTensor<Out>;
    fn contract(&self, other: &SparseTensor<T>) -> Option<Self::LCM> {
        // println!("Contracting SparseTensor SparseTensor");
        if let Some((i, j)) = self.match_index(other) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = AHashMap::new();

            let metric = self.get_ith_metric(i).unwrap();

            for (index_a, nonzeros_a, fiber_a) in self.iter_fibers(i) {
                for (index_b, nonzeros_b, fiber_b) in other.iter_fibers(j) {
                    let result_index = index_a[..i]
                        .iter()
                        .chain(&index_a[i + 1..])
                        .chain(&index_b[..j])
                        .chain(&index_b[j + 1..])
                        .cloned()
                        .collect::<Vec<_>>();

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
                        result_data.insert(result_index, value);
                    }
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

impl<T, U, Out> Contract<SparseTensor<T>> for DenseTensor<U>
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
{
    type LCM = DenseTensor<Out>;
    fn contract(&self, other: &SparseTensor<T>) -> Option<Self::LCM> {
        // self is dense U, other is sparse T
        if let Some((i, j)) = self.match_index(other) {
            // final structure is the structure of self appended with the structure of other, with the contracted indices removed
            let final_structure = self.structure().merge_at(other.structure(), (i, j));

            let mut result_data = vec![Out::zero(); final_structure.size()];

            let metric = other.get_ith_metric(j).unwrap();

            for (index_a, fiber_a) in self.iter_fibers(i) {
                for (index_b, nonzeros, fiber_b) in other.iter_fibers(j) {
                    //nonzeros is the indices of the non-zero elements of the other tensor
                    let result_index = final_structure
                        .flat_index(
                            &index_a[..i]
                                .iter()
                                .chain(&index_a[i + 1..])
                                .chain(&index_b[..j])
                                .chain(&index_b[j + 1..])
                                .cloned()
                                .collect::<Vec<_>>(),
                        )
                        .unwrap();

                    for (i, k) in nonzeros.iter().enumerate() {
                        // Adjust indices for fetching from the other tensor
                        if metric[*k] {
                            result_data[result_index] -= fiber_a[*k] * fiber_b[i];
                        } else {
                            result_data[result_index] += fiber_a[*k] * fiber_b[i];
                        }
                    }
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

impl<T, U, Out> Contract<NumTensor<T>> for NumTensor<U>
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
{
    type LCM = NumTensor<Out>;
    fn contract(&self, other: &NumTensor<T>) -> Option<NumTensor<Out>> {
        match (self, other) {
            (NumTensor::Dense(s), NumTensor::Dense(o)) => Some(NumTensor::Dense(s.contract(o)?)),
            (NumTensor::Dense(s), NumTensor::Sparse(o)) => Some(NumTensor::Dense(s.contract(o)?)),
            (NumTensor::Sparse(s), NumTensor::Dense(o)) => Some(NumTensor::Dense(s.contract(o)?)),
            (NumTensor::Sparse(s), NumTensor::Sparse(o)) => Some(NumTensor::Sparse(s.contract(o)?)),
        }
    }
}

impl Contract<NumTensors> for NumTensors {
    type LCM = NumTensors;
    fn contract(&self, other: &NumTensors) -> Option<Self::LCM> {
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

pub struct TensorNetwork {
    graph: Graph<NumTensors, Slot, Undirected>,
}

impl TensorNetwork {
    pub fn new(tensors: Vec<NumTensors>) -> Self {
        TensorNetwork {
            graph: Self::generate_network_graph(tensors),
        }
    }

    fn generate_network_graph(tensors: Vec<NumTensors>) -> Graph<NumTensors, Slot, Undirected> {
        let mut graph = Graph::<NumTensors, Slot, Undirected>::new_undirected();

        for tensor in tensors {
            graph.add_node(tensor);
        }

        for (ni, n) in graph.node_indices().enumerate() {
            for (mi, m) in graph.node_indices().enumerate().filter(|(i, _)| *i > ni) {
                let a = graph.node_weight(n).unwrap();
                let b = graph.node_weight(m).unwrap();

                if let Some((i, j)) = a.match_index(b) {
                    graph.add_edge(n, m, a.structure()[i]);
                }
            }
        }
        graph
    }

    fn contract_edge(&mut self, edge_idx: EdgeIndex) {
        let (a, b) = self.graph.edge_endpoints(edge_idx).unwrap();

        let ai = self.graph.node_weight(a).unwrap();
        let bi = self.graph.node_weight(b).unwrap();

        let f = ai.contract(bi).unwrap();

        self.merge_nodes(a, b, f);
    }

    fn merge_nodes(&mut self, a: NodeIndex, b: NodeIndex, weight: NumTensors) {
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

    pub fn dot(&self) -> String {
        format!(
            "{:?}",
            Dot::with_attr_getters(
                &self.graph,
                &[Config::EdgeNoLabel, Config::NodeNoLabel],
                &|_, e| { format!("{}", e.weight()) },
                &|_, n| { format!("{}", n.1.structure()) }
            )
        )
    }

    pub fn contract(&mut self) {
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

        if let Some(e) = edge_to_min_degree_node {
            self.contract_edge(e);
            self.contract();
        }
    }

    pub fn result(&self) -> NumTensors {
        self.graph
            .node_weight(self.graph.node_indices().next().unwrap())
            .unwrap()
            .clone()
    }
}
