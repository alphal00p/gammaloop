use ahash::AHashMap;
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
    representations::Identifier,
    state::{State, Workspace},
};

use self::mixed_tensor::{MixedTensor, MixedTensors, SymbolicContract};

use super::*;
use smartstring::alias::String;
use std::{fmt::Debug, ops::Neg};
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
pub trait Contract<T> {
    type LCM;
    fn contract(&self, other: &T) -> Option<Self::LCM>;
}

impl<T, U, Out, I> Contract<DenseTensor<T, I>> for DenseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = Out>,
    Out: std::ops::AddAssign<Out>
        + std::ops::SubAssign<Out>
        + Neg<Output = Out>
        + Clone
        + num::traits::Zero
        + for<'a> std::ops::AddAssign<&'a Out>
        + for<'b> std::ops::SubAssign<&'b Out>
        + std::fmt::Debug,
    I: Clone + Debug,
    T: Debug,
    U: Debug,
{
    type LCM = DenseTensor<Out, I>;
    fn contract(&self, other: &DenseTensor<T, I>) -> Option<Self::LCM> {
        if let Some((single, i, j)) = self.structure().match_index(other.structure()) {
            if i >= j {
                if single {
                    let final_structure = self.structure().merge_at(other.structure(), (i, j));
                    let mut result_data = vec![Out::zero(); final_structure.size()];
                    let mut result_index = 0;

                    let mut self_iter = self.iter_fiber(i);
                    let mut other_iter = other.iter_fiber(j);
                    let dimension = self_iter.fiber_dimension;
                    let metric = self.structure()[i].representation.negative();

                    for fiber_a in self_iter.by_ref() {
                        for fiber_b in other_iter.by_ref() {
                            for i in 0..dimension {
                                if metric[i] {
                                    result_data[result_index] -= fiber_a[i] * fiber_b[i];
                                } else {
                                    result_data[result_index] += fiber_a[i] * fiber_b[i];
                                }
                            }
                            result_index += 1;
                        }
                        other_iter.reset();
                    }

                    let result = DenseTensor {
                        data: result_data,
                        structure: final_structure,
                    };

                    return Some(result);
                } else {
                    // println!("SparseTensor DenseTensor");
                    let (permutation, self_matches, other_matches) =
                        self.structure().match_indices(other.structure()).unwrap();

                    let mut final_structure = self.structure().clone();
                    final_structure.merge(other.structure());

                    // Initialize result tensor with default values
                    let mut result_data = vec![Out::zero(); final_structure.size()];
                    let mut result_index = 0;
                    let mut selfiter = self.iter_multi_fibers_metric(&self_matches, permutation);

                    let mut other_iter = other.iter_multi_fibers(&other_matches);
                    while let Some(fiber_a) = selfiter.next() {
                        for fiber_b in other_iter.by_ref() {
                            for k in 0..fiber_a.len() {
                                if fiber_a[k].1 {
                                    result_data[result_index] -=
                                        fiber_b[selfiter.map[k]] * fiber_a[k].0;
                                } else {
                                    result_data[result_index] +=
                                        fiber_b[selfiter.map[k]] * fiber_a[k].0;
                                }
                            }
                            result_index += 1;
                        }
                        other_iter.reset();
                    }
                    let result: DenseTensor<Out, I> = DenseTensor {
                        data: result_data,
                        structure: final_structure,
                    };

                    return Some(result);
                }
            } else {
                return other.contract(self);
            }
        }
        None
    }
}

impl<T, U, I, Out> Contract<DenseTensor<T, I>> for SparseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = Out>,
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
        if let Some((single, i, j)) = self.structure().match_index(other.structure()) {
            if i >= j {
                if single {
                    let final_structure = self.structure().merge_at(other.structure(), (i, j));
                    let mut result_data = vec![Out::zero(); final_structure.size()];
                    let mut result_index = 0;

                    let mut self_iter = self.iter_fiber(i);
                    let mut other_iter = other.iter_fiber(j);

                    let metric = self.structure()[i].representation.negative();

                    for (skipped, nonzeros, fiber_a) in self_iter.by_ref() {
                        result_index += skipped;
                        for fiber_b in other_iter.by_ref() {
                            for (i, k) in nonzeros.iter().enumerate() {
                                if metric[*k] {
                                    result_data[result_index] -= fiber_a[i] * fiber_b[*k];
                                } else {
                                    result_data[result_index] += fiber_a[i] * fiber_b[*k];
                                }
                            }
                            result_index += 1;
                        }
                        other_iter.reset();
                    }

                    let result = DenseTensor {
                        data: result_data,
                        structure: final_structure,
                    };

                    Some(result)
                } else {
                    // println!("SparseTensor DenseTensor");
                    let (permutation, self_matches, other_matches) =
                        self.structure().match_indices(other.structure()).unwrap();

                    let mut final_structure = self.structure().clone();
                    final_structure.merge(other.structure());

                    let mut result_data = vec![Out::zero(); final_structure.size()];
                    let mut result_index = 0;

                    let mut selfiter = self.iter_multi_fibers_metric(&self_matches, permutation);
                    let mut other_iter = other.iter_multi_fibers(&other_matches);

                    while let Some((skipped, nonzeros, fiber_a)) = selfiter.next() {
                        result_index += skipped;
                        for fiber_b in other_iter.by_ref() {
                            for (i, k) in nonzeros.iter().enumerate() {
                                if fiber_a[i].1 {
                                    result_data[result_index] -=
                                        fiber_a[i].0 * fiber_b[selfiter.map[*k]];
                                } else {
                                    result_data[result_index] +=
                                        fiber_a[i].0 * fiber_b[selfiter.map[*k]];
                                }
                            }
                            result_index += 1;
                        }
                        other_iter.reset();
                    }
                    let result = DenseTensor {
                        data: result_data,
                        structure: final_structure,
                    };

                    Some(result)
                }
            } else {
                other.contract(self)
            }
        } else {
            None
        }
    }
}

impl<T, U, Out, I> Contract<SparseTensor<T, I>> for SparseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = Out>,
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
        if let Some((single, i, j)) = self.structure().match_index(other.structure()) {
            if i >= j {
                if single {
                    let final_structure = self.structure().merge_at(other.structure(), (i, j));
                    let mut result_data = AHashMap::default();
                    let mut result_index = 0;

                    let mut self_iter = self.iter_fiber(i);
                    let mut other_iter = other.iter_fiber(j);

                    let metric = self.structure()[i].representation.negative();

                    let stride = *final_structure
                        .strides()
                        .get(i.saturating_sub(1))
                        .unwrap_or(&1);

                    for (skipped_a, nonzeros_a, fiber_a) in self_iter.by_ref() {
                        result_index += skipped_a;
                        for (skipped_b, nonzeros_b, fiber_b) in other_iter.by_ref() {
                            result_index += skipped_b * stride;
                            let mut value = Out::zero();
                            let mut nonzero = false;

                            for (i, j, x) in nonzeros_a.iter().enumerate().filter_map(|(i, &x)| {
                                nonzeros_b.binary_search(&x).ok().map(|j| (i, j, x))
                            }) {
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
                            result_index += 1;
                        }
                        other_iter.reset();
                    }

                    let result = SparseTensor {
                        elements: result_data,
                        structure: final_structure,
                    };

                    return Some(result);
                } else {
                    let (permutation, self_matches, other_matches) =
                        self.structure().match_indices(other.structure()).unwrap();

                    let mut final_structure = self.structure().clone();
                    final_structure.merge(other.structure());
                    let mut result_data = AHashMap::default();
                    let one = if let Some(o) = final_structure.strides().first() {
                        *o
                    } else {
                        1
                    };
                    let stride = *final_structure
                        .strides()
                        .get(i.checked_sub(1)?)
                        .unwrap_or(&one);

                    let mut result_index = 0;

                    let mut self_iter = self.iter_multi_fibers_metric(&self_matches, permutation);
                    let mut other_iter = other.iter_multi_fibers(&other_matches);
                    while let Some((skipped_a, nonzeros_a, fiber_a)) = self_iter.next() {
                        result_index += skipped_a;
                        for (skipped_b, nonzeros_b, fiber_b) in other_iter.by_ref() {
                            result_index += skipped_b * stride;
                            let mut value = Out::zero();
                            let mut nonzero = false;

                            for (i, j, _x) in nonzeros_a.iter().enumerate().filter_map(|(i, &x)| {
                                nonzeros_b
                                    .binary_search(&self_iter.map[x])
                                    .ok()
                                    .map(|j| (i, j, x))
                            }) {
                                if fiber_a[i].1 {
                                    value -= fiber_a[i].0 * fiber_b[j];
                                } else {
                                    value += fiber_a[i].0 * fiber_b[j];
                                }
                                nonzero = true;
                            }

                            if nonzero && value != Out::zero() {
                                result_data.insert(result_index, value);
                            }
                            result_index += 1;
                        }
                        other_iter.reset();
                    }

                    let result = SparseTensor {
                        elements: result_data,
                        structure: final_structure,
                    };

                    return Some(result);
                }
            } else {
                return other.contract(self);
            }
        }
        None
    }
}

impl<T, U, Out, I> Contract<SparseTensor<T, I>> for DenseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = Out>,
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
        if let Some((single, i, j)) = self.structure().match_index(other.structure()) {
            if i >= j {
                if single {
                    let final_structure = self.structure().merge_at(other.structure(), (i, j));
                    let mut result_data = vec![Out::zero(); final_structure.size()];
                    let mut result_index = 0;

                    let mut self_iter = self.iter_fiber(i);
                    let mut other_iter = other.iter_fiber(j);

                    let metric = self.structure()[i].representation.negative();

                    let stride = *final_structure
                        .strides()
                        .get(i.saturating_sub(1))
                        .unwrap_or(&1);

                    for fiber_a in self_iter.by_ref() {
                        for (skipped, nonzeros, fiber_b) in other_iter.by_ref() {
                            result_index += skipped * stride;
                            for (i, k) in nonzeros.iter().enumerate() {
                                if metric[*k] {
                                    result_data[result_index] -= fiber_a[*k] * fiber_b[i];
                                } else {
                                    result_data[result_index] += fiber_a[*k] * fiber_b[i];
                                }
                            }
                            result_index += 1;
                        }
                        other_iter.reset();
                    }

                    let result = DenseTensor {
                        data: result_data,
                        structure: final_structure,
                    };

                    Some(result)
                } else {
                    let (permutation, self_matches, other_matches) =
                        self.structure().match_indices(other.structure()).unwrap();

                    let mut final_structure = self.structure().clone();
                    final_structure.merge(other.structure());

                    let mut result_data = vec![Out::zero(); final_structure.size()];
                    let mut result_index = 0;

                    let one = if let Some(o) = final_structure.strides().first() {
                        *o
                    } else {
                        1
                    };
                    let stride = *final_structure
                        .strides()
                        .get(i.saturating_sub(1))
                        .unwrap_or(&one);

                    let mut selfiter = self.iter_multi_fibers_metric(&self_matches, permutation);

                    while let Some(fiber_a) = selfiter.next() {
                        for (skipped, nonzeros, fiber_b) in other.iter_multi_fibers(&other_matches)
                        {
                            result_index += skipped * stride;
                            for (i, k) in nonzeros.iter().enumerate() {
                                if fiber_a[*k].1 {
                                    result_data[result_index] -=
                                        fiber_a[selfiter.map[*k]].0 * fiber_b[i];
                                } else {
                                    result_data[result_index] +=
                                        fiber_a[selfiter.map[*k]].0 * fiber_b[i];
                                }
                            }
                            result_index += 1;
                        }
                    }

                    let result = DenseTensor {
                        data: result_data,
                        structure: final_structure,
                    };

                    Some(result)
                }
            } else {
                other.contract(self)
            }
        } else {
            None
        }
    }
}

impl<T, U, I, Out> Contract<NumTensor<T, I>> for NumTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = Out>,
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
    I: Clone + Debug,
    T: Clone + Debug,
    U: Clone + Debug,
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
    I: Clone + Debug,
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

                if let Some((_, i, _)) = a.structure().match_index(b.structure()) {
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
                &|_, e| { format!("label=\"{}\"", e.weight()) },
                &|_, n| { format!("label=\"{}\"", n.1.structure()) }
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
                &|_, e| { format!("label=\"{}\"", e.weight()) },
                &|_, n| { format!("label=\"{}\"", n.1.structure().to_string(state)) }
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
