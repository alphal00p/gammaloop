use ahash::AHashMap;
use itertools::Itertools;

use num::traits::Inv;
use petgraph::{
    dot::{Config, Dot},
    graph::{self, EdgeIndex, Graph, NodeIndex},
    visit::EdgeRef,
    Direction::Outgoing,
    Undirected,
};

use serde::{Deserialize, Serialize};
use slotmap::{new_key_type, DenseSlotMap, Key, SecondaryMap, SlotMap};
use symbolica::{
    representations::Identifier,
    state::{State, Workspace},
};

use self::mixed_tensor::{MixedTensor, MixedTensors, SymbolicContract};
use self::tensor_structure::HistoryStructure;

use super::{
    mixed_tensor, tensor_structure, Atom, DataTensor, DenseTensor, GetTensorData, HasName,
    HasTensorData, NumTensor, SetTensorData, Shadowable, Slot, SparseTensor, StructureContract,
    SymbolicZero, TensorStructure, TracksCount,
};
use smartstring::alias::String;
use std::{
    collections::BTreeMap,
    fmt::{Debug, Display},
    // intrinsics::needs_drop,
    ops::Neg,
    result,
};

impl<T, I> DenseTensor<T, I>
where
    T: for<'a> std::ops::AddAssign<&'a T>
        + for<'b> std::ops::SubAssign<&'b T>
        + Neg<Output = T>
        + Clone
        + std::fmt::Debug,
    I: TensorStructure + Clone + StructureContract,
{
    #[must_use]

    /// Contract the tensor with itself, i.e. trace over all matching indices.
    pub fn internal_contract(&self) -> Self {
        let mut result: DenseTensor<T, I> = self.clone();
        for trace in self.traces() {
            let mut new_structure = self.structure.clone();
            new_structure.trace(trace[0], trace[1]);

            let mut new_result = DenseTensor::from_data_coerced(&self.data, new_structure)
                .unwrap_or_else(|_| unreachable!());
            for (idx, t) in result.iter_trace(trace) {
                new_result.set(&idx, t).unwrap_or_else(|_| unreachable!());
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
    I: TensorStructure + Clone + StructureContract,
{
    #[must_use]
    /// Contract the tensor with itself, i.e. trace over all matching indices.
    pub fn internal_contract(&self) -> Self {
        let trace = self.traces()[0];

        // println!("trace {:?}", trace);
        let mut new_structure = self.structure.clone();
        new_structure.trace(trace[0], trace[1]);

        let mut new_result = SparseTensor::empty(new_structure);
        for (idx, t) in self.iter_trace(trace).filter(|(_, t)| *t != T::default()) {
            new_result.set(&idx, t).unwrap_or_else(|_| unreachable!());
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
    I: TensorStructure + Clone + StructureContract,
    T: Debug,
    U: Debug,
{
    type LCM = DenseTensor<Out, I>;
    fn contract(&self, other: &DenseTensor<T, I>) -> Option<Self::LCM> {
        if let Some((single, i, j)) = self.structure.match_index(&other.structure) {
            if i >= j {
                if single {
                    let final_structure = self.structure.merge_at(&other.structure, (i, j));
                    let mut result_data = vec![Out::zero(); final_structure.size()];
                    let mut result_index = 0;

                    let mut self_iter = self.iter_fiber(i);
                    let mut other_iter = other.iter_fiber(j);
                    let dimension = self_iter.fiber_dimension;
                    let metric = self.external_structure()[i].representation.negative();

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
                }
                // println!("SparseTensor DenseTensor");
                let (permutation, self_matches, other_matches) =
                    self.structure().match_indices(other.structure()).unwrap();

                let mut final_structure = self.structure.clone();
                final_structure.merge(&other.structure);

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
            return other.contract(self);
        }
        let mut final_structure = self.structure.clone();
        final_structure.merge(&other.structure);

        let mut result_data = vec![Out::zero(); final_structure.size()];
        let stride = other.size();

        for (i, u) in self.iter_flat() {
            for (j, t) in other.iter_flat() {
                result_data[i * stride + j] = u * t;
            }
        }
        let result = DenseTensor {
            data: result_data,
            structure: final_structure,
        };
        Some(result)
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
    I: TensorStructure + Clone + StructureContract,
    U: Clone,
{
    type LCM = DenseTensor<Out, I>;
    fn contract(&self, other: &DenseTensor<T, I>) -> Option<Self::LCM> {
        if let Some((single, i, j)) = self.structure.match_index(&other.structure) {
            if i >= j {
                if single {
                    let final_structure = self.structure.merge_at(&other.structure, (i, j));
                    let mut result_data = vec![Out::zero(); final_structure.size()];
                    let mut result_index = 0;

                    let mut self_iter = self.iter_fiber(i);
                    let mut other_iter = other.iter_fiber(j);

                    let metric = self.external_structure()[i].representation.negative();

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

                    let mut final_structure = self.structure.clone();
                    final_structure.merge(&other.structure);

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
            let mut final_structure = self.structure.clone();
            final_structure.merge(&other.structure);

            let mut result_data = vec![Out::zero(); final_structure.size()];
            let stride = other.size();

            for (i, u) in self.iter_flat() {
                for (j, t) in other.iter_flat() {
                    result_data[i * stride + j] = u * t;
                }
            }
            let result = DenseTensor {
                data: result_data,
                structure: final_structure,
            };
            Some(result)
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
    I: Clone + TensorStructure + StructureContract,
    U: Clone,
    T: Clone,
{
    type LCM = SparseTensor<Out, I>;
    #[allow(clippy::too_many_lines)]
    fn contract(&self, other: &SparseTensor<T, I>) -> Option<Self::LCM> {
        if let Some((single, i, j)) = self.structure.match_index(&other.structure) {
            if i >= j {
                if single {
                    let final_structure = self.structure.merge_at(&other.structure, (i, j));
                    let mut result_data = AHashMap::default();
                    let mut result_index = 0;

                    let mut self_iter = self.iter_fiber(i);
                    let mut other_iter = other.iter_fiber(j);

                    let metric = self.external_structure()[i].representation.negative();

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
                }
                let (permutation, self_matches, other_matches) =
                    self.structure().match_indices(other.structure()).unwrap();

                let mut final_structure = self.structure.clone();
                final_structure.merge(&other.structure);
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
            return other.contract(self);
        }
        let mut final_structure = self.structure.clone();
        final_structure.merge(&other.structure);

        let mut result_data = AHashMap::default();
        let stride = other.size();

        for (i, u) in self.iter_flat() {
            for (j, t) in other.iter_flat() {
                result_data.insert(i * stride + j, u * t);
            }
        }
        let result = SparseTensor {
            elements: result_data,
            structure: final_structure,
        };
        Some(result)
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
    I: Clone + TensorStructure + StructureContract,
    T: Clone,
{
    type LCM = DenseTensor<Out, I>;
    fn contract(&self, other: &SparseTensor<T, I>) -> Option<Self::LCM> {
        // self is dense U, other is sparse T
        if let Some((single, i, j)) = self.structure.match_index(&other.structure) {
            if i >= j {
                if single {
                    let final_structure = self.structure.merge_at(&other.structure, (i, j));
                    let mut result_data = vec![Out::zero(); final_structure.size()];
                    let mut result_index = 0;

                    let mut self_iter = self.iter_fiber(i);
                    let mut other_iter = other.iter_fiber(j);

                    let metric = self.external_structure()[i].representation.negative();

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

                    let mut final_structure = self.structure.clone();
                    final_structure.merge(&other.structure);

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
            let mut final_structure = self.structure.clone();
            final_structure.merge(&other.structure);

            let mut result_data = vec![Out::zero(); final_structure.size()];
            let stride = other.size();

            for (i, u) in self.iter_flat() {
                for (j, t) in other.iter_flat() {
                    result_data[i * stride + j] = u * t;
                }
            }
            let result = DenseTensor {
                data: result_data,
                structure: final_structure,
            };
            Some(result)
        }
    }
}

impl<T, U, I, Out> Contract<DataTensor<T, I>> for DataTensor<U, I>
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
    I: Clone + Debug + TensorStructure + StructureContract,
    T: Clone + Debug,
    U: Clone + Debug,
{
    type LCM = DataTensor<Out, I>;
    fn contract(&self, other: &DataTensor<T, I>) -> Option<DataTensor<Out, I>> {
        match (self, other) {
            (DataTensor::Dense(s), DataTensor::Dense(o)) => Some(DataTensor::Dense(s.contract(o)?)),
            (DataTensor::Dense(s), DataTensor::Sparse(o)) => {
                Some(DataTensor::Dense(s.contract(o)?))
            }
            (DataTensor::Sparse(s), DataTensor::Dense(o)) => {
                Some(DataTensor::Dense(s.contract(o)?))
            }
            (DataTensor::Sparse(s), DataTensor::Sparse(o)) => {
                Some(DataTensor::Sparse(s.contract(o)?))
            }
        }
    }
}

impl<I> Contract<NumTensor<I>> for NumTensor<I>
where
    I: Clone + Debug + TensorStructure + StructureContract,
{
    type LCM = NumTensor<I>;
    fn contract(&self, other: &NumTensor<I>) -> Option<Self::LCM> {
        match (self, other) {
            (NumTensor::Float(a), NumTensor::Float(b)) => Some(NumTensor::Float(a.contract(b)?)),
            (NumTensor::Float(a), NumTensor::Complex(b)) => {
                Some(NumTensor::Complex(a.contract(b)?))
            }
            (NumTensor::Complex(a), NumTensor::Float(b)) => {
                Some(NumTensor::Complex(a.contract(b)?))
            }
            (NumTensor::Complex(a), NumTensor::Complex(b)) => {
                Some(NumTensor::Complex(a.contract(b)?))
            }
        }
    }
}

new_key_type! {
    pub struct NodeId;
    pub struct HedgeId;
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HalfEdgeGraph<N, E> {
    pub edges: DenseSlotMap<HedgeId, E>,
    pub involution: SecondaryMap<HedgeId, HedgeId>,
    pub neighbors: SecondaryMap<HedgeId, HedgeId>,
    pub nodes: DenseSlotMap<NodeId, N>,
    pub nodemap: SecondaryMap<HedgeId, NodeId>,
    pub reverse_nodemap: SecondaryMap<NodeId, HedgeId>,
}

struct IncidentIterator<'a> {
    neighbors: &'a SecondaryMap<HedgeId, HedgeId>,
    current: Option<HedgeId>,
    start: HedgeId,
}

impl<'a> Iterator for IncidentIterator<'a> {
    type Item = HedgeId;
    fn next(&mut self) -> Option<HedgeId> {
        let current = self.current?;

        self.current = Some(self.neighbors[current]);

        if self.current == Some(self.start) {
            self.current = None;
        }

        Some(current)
    }
}

impl<'a> IncidentIterator<'a> {
    fn new<N, E>(graph: &'a HalfEdgeGraph<N, E>, initial: HedgeId) -> Self {
        IncidentIterator {
            neighbors: &graph.neighbors,
            current: Some(initial),
            start: initial,
        }
    }
}

impl<N, E> HalfEdgeGraph<N, E> {
    fn new() -> Self {
        HalfEdgeGraph {
            involution: SecondaryMap::new(),
            nodemap: SecondaryMap::new(),
            neighbors: SecondaryMap::new(),
            reverse_nodemap: SecondaryMap::new(),
            nodes: DenseSlotMap::with_key(),
            edges: DenseSlotMap::with_key(),
        }
    }

    fn dot(&self) -> std::string::String {
        let mut out = "graph {".to_string();
        out.push_str("  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\";");

        // for (i, n) in &self.nodes {
        //     out.push_str(&format!("\n {}", i.data().as_ffi()));
        // }
        for (i, _) in &self.neighbors {
            if i > self.involution[i] {
                out.push_str(&format!(
                    "\n {} -- {}",
                    self.nodemap[i].data().as_ffi(),
                    self.nodemap[self.involution[i]].data().as_ffi()
                ));
            } else if i == self.involution[i] {
                out.push_str(&format!(
                    "ext{} [shape=none, label=\"\"];",
                    i.data().as_ffi()
                ));
                out.push_str(&format!(
                    "\n {} -- ext{}",
                    self.nodemap[i].data().as_ffi(),
                    i.data().as_ffi()
                ));
            }
        }

        out += "}";
        out
    }

    fn add_node(&mut self, data: N) -> NodeId {
        self.nodes.insert(data)
    }

    fn node_indices(&self) -> slotmap::dense::Keys<'_, NodeId, N> {
        self.nodes.keys()
    }

    /// Add a node with a list of edget with associated data. Matches edges by equality.
    fn add_node_with_edges(&mut self, data: N, edges: &[E]) -> NodeId
    where
        E: Eq + Clone,
    {
        let idx = self.add_node(data);
        for e in edges {
            let mut found_match = false;
            for (i, other_e) in &self.edges {
                if *e == *other_e && self.involution[i] == i {
                    found_match = true;
                    let eid = self.edges.insert(e.clone());
                    self.involution.insert(eid, i);
                    self.involution.insert(i, eid);
                    self.nodemap.insert(eid, idx);
                    if let Some(prev_eid) = self.reverse_nodemap.insert(idx, eid) {
                        let next_eid = self.neighbors.insert(prev_eid, eid).unwrap();
                        self.neighbors.insert(eid, next_eid);
                    } else {
                        self.neighbors.insert(eid, eid);
                    }
                    break;
                }
            }
            if !found_match {
                let eid = self.edges.insert(e.clone());
                self.involution.insert(eid, eid);
                self.nodemap.insert(eid, idx);
                if let Some(prev_eid) = self.reverse_nodemap.insert(idx, eid) {
                    let next_eid = self.neighbors.insert(prev_eid, eid).unwrap();
                    self.neighbors.insert(eid, next_eid);
                } else {
                    self.neighbors.insert(eid, eid);
                }
            }
        }

        idx
    }

    pub fn validate_neighbors(&self) -> bool {
        for (i, n) in &self.reverse_nodemap {
            for j in IncidentIterator::new(self, *n) {
                if self.nodemap[j] != i {
                    return false;
                }
            }
        }
        true
    }

    fn node_labels(&self) -> String
    where
        N: Display,
    {
        let mut out = String::new();
        for (i, n) in &self.nodes {
            out.push_str(&format!("{}[label= \"{}\"]\n", i.data().as_ffi(), n));
        }
        out
    }

    #[allow(clippy::too_many_lines)]
    fn merge_nodes(&mut self, a: NodeId, b: NodeId, data: N) -> NodeId {
        let c = self.nodes.insert(data);

        // New initial edge for reverse_nodemap, that does not link to b
        // if none is found, all incident edges are link to b and must be removed from the neighbors list
        let mut new_initial_a = self
            .edges_incident(a)
            .find(|x| self.nodemap[self.involution[*x]] != b && self.involution[*x] != *x);

        if new_initial_a.is_none() {
            new_initial_a = self
                .edges_incident(a)
                .find(|x| self.nodemap[self.involution[*x]] != b);
        }

        let mut a_edge = new_initial_a;

        if a_edge.is_none() {
            // all edges link to b, and must be removed
            let initial = self.reverse_nodemap[a];
            let mut current = Some(initial);
            loop {
                if current.is_none() {
                    break;
                }
                let next = self.neighbors.remove(current.unwrap());

                if next == Some(initial) {
                    current = None;
                } else {
                    current = next;
                }
            }
        } else {
            loop {
                let mut next = self.neighbors[a_edge.unwrap()];

                while self.nodemap[self.involution[next]] == b {
                    next = self.neighbors.remove(next).unwrap();
                }

                self.nodemap.insert(a_edge.unwrap(), c);
                self.neighbors.insert(a_edge.unwrap(), next);

                if new_initial_a == Some(next) {
                    break;
                }

                a_edge = Some(next);
            }
        }

        let mut new_initial_b = self
            .edges_incident(b)
            .find(|x| self.nodemap[self.involution[*x]] != a && self.involution[*x] != *x);

        if new_initial_b.is_none() {
            new_initial_b = self
                .edges_incident(b)
                .find(|x| self.nodemap[self.involution[*x]] != a);
        }
        let mut b_edge = new_initial_b;

        if b_edge.is_none() {
            let initial = self.reverse_nodemap[b];
            let mut current = Some(initial);
            loop {
                if current.is_none() {
                    break;
                }
                let next = self.neighbors.remove(current.unwrap());

                if next == Some(initial) {
                    current = None;
                } else {
                    current = next;
                }
            }
        } else {
            loop {
                let mut next = self.neighbors[b_edge.unwrap()];

                while self.nodemap[self.involution[next]] == a {
                    next = self.neighbors.remove(next).unwrap();
                }

                self.nodemap.insert(b_edge.unwrap(), c);
                self.neighbors.insert(b_edge.unwrap(), next);

                if new_initial_b == Some(next) {
                    break;
                }

                b_edge = Some(next);
            }
        }

        match (new_initial_a, new_initial_b) {
            (Some(new_edge_a), Some(new_edge_b)) => {
                self.reverse_nodemap.insert(c, new_edge_a);
                self.reverse_nodemap.remove(a);
                self.reverse_nodemap.remove(b);
                let old_neig = self.neighbors.insert(new_edge_a, new_edge_b).unwrap();
                self.neighbors.insert(b_edge.unwrap(), old_neig);
            }
            (Some(new_edge_a), None) => {
                self.reverse_nodemap.insert(c, new_edge_a);
                self.reverse_nodemap.remove(a);
                self.reverse_nodemap.remove(b);
            }
            (None, Some(new_edge_b)) => {
                self.reverse_nodemap.insert(c, new_edge_b);
                self.reverse_nodemap.remove(a);
                self.reverse_nodemap.remove(b);
            }
            (None, None) => {
                self.reverse_nodemap.remove(b);
                self.reverse_nodemap.remove(a);
            }
        }
        self.nodes.remove(a);
        self.nodes.remove(b);
        c
    }

    /// Add an internal edge between two nodes.
    fn add_edge(&mut self, a: NodeId, b: NodeId, data: E) -> HedgeId
    where
        E: Clone,
    {
        let hedge_id_a = self.edges.insert(data.clone());
        let hedge_id_b = self.edges.insert(data);
        self.involution.insert(hedge_id_a, hedge_id_b);
        self.involution.insert(hedge_id_b, hedge_id_a);
        self.nodemap.insert(hedge_id_a, a);
        if let Some(prev_eid) = self.reverse_nodemap.insert(a, hedge_id_a) {
            let next_eid = self.neighbors.insert(prev_eid, hedge_id_a).unwrap();
            self.neighbors.insert(hedge_id_a, next_eid).unwrap();
        } else {
            self.neighbors.insert(hedge_id_a, hedge_id_a);
        }
        self.nodemap.insert(hedge_id_b, b);
        if let Some(prev_eid) = self.reverse_nodemap.insert(b, hedge_id_b) {
            let next_eid = self.neighbors.insert(prev_eid, hedge_id_b).unwrap();
            self.neighbors.insert(hedge_id_b, next_eid).unwrap();
        } else {
            self.neighbors.insert(hedge_id_b, hedge_id_b);
        }
        hedge_id_a
    }

    /// Add external, as a fixed point involution half edge.
    fn add_external(&mut self, a: NodeId, data: E) -> HedgeId {
        let id = self.edges.insert(data);
        self.involution.insert(id, id);
        self.nodemap.insert(id, a);
        if let Some(prev_eid) = self.reverse_nodemap.insert(a, id) {
            let next_eid = self.neighbors.insert(prev_eid, id).unwrap();
            self.neighbors.insert(id, next_eid).unwrap();
        } else {
            self.neighbors.insert(id, id);
        }
        id
    }

    fn edges_incident(&self, node: NodeId) -> impl Iterator<Item = HedgeId> + '_ {
        IncidentIterator::new(self, self.reverse_nodemap[node])
    }

    fn edges_between(&self, a: NodeId, b: NodeId) -> impl Iterator<Item = HedgeId> + '_ {
        self.edges_incident(a)
            .filter(move |&i| self.nodemap[self.involution[i]] == b)
    }

    fn internal_edges_incident(&self, node: NodeId) -> impl Iterator<Item = HedgeId> + '_ {
        self.edges_incident(node)
            .filter(move |&i| self.nodemap[self.involution[i]] != node)
    }

    fn external_edges_incident(&self, node: NodeId) -> impl Iterator<Item = HedgeId> + '_ {
        self.edges_incident(node)
            .filter(move |&i| self.nodemap[self.involution[i]] == node)
    }

    fn degree(&self, node: NodeId) -> usize {
        self.edges_incident(node).collect::<Vec<_>>().len()
    }

    fn neighbors(&self, node: NodeId) -> impl Iterator<Item = NodeId> + '_ {
        self.edges_incident(node)
            .map(move |i| self.nodemap[self.involution[i]])
    }

    // fn map_nodes<F, U>(&self, f: F) -> HalfEdgeGraph<U, E>
    // where
    //     F: Fn(&N) -> U,
    //     E: Clone,
    // {
    //     let edges = self.edges.clone();
    //     let involution = self.involution.clone();

    //     let mut nodes = DenseSlotMap::with_key();
    //     let mut nodemap = SecondaryMap::new();

    //     for n in &self.nodes {
    //         let nid = nodes.insert(f(n.1));
    //         for e in self.edges_incident(n.0) {
    //             nodemap.insert(e, nid);
    //         }
    //     }

    //     HalfEdgeGraph {
    //         edges,
    //         involution,
    //         nodes,
    //         nodemap,
    //     }
    // }
}

#[test]
fn merge() {
    let mut graph = HalfEdgeGraph::new();
    let a = graph.add_node_with_edges(1, &[1, 2, 3, 4, 5]);
    let b = graph.add_node_with_edges(2, &[1, 2, 6, 7, 8]);
    let c = graph.add_node_with_edges(4, &[4, 6, 9, 10, 11]);

    println!("{}", graph.dot());
    println!("{}", graph.degree(a));
    println!("{}", graph.degree(b));

    for (i, n) in &graph.neighbors {
        println!("{} {}", graph.edges[i], graph.edges[*n]);
    }

    let d = graph.merge_nodes(a, b, 3);

    // for (i, n) in &graph.neighbors {
    //     println!("{} {}", graph.edges[i], graph.edges[*n]);
    // }

    println!("{}", graph.dot());
    println!("{}", graph.degree(c));
    println!("{}", graph.neighbors.len());

    let e = graph.merge_nodes(c, d, 5);

    println!("{}", graph.dot());
    println!("{}", graph.degree(e));
    println!("{}", graph.neighbors.len());

    let mut graph = HalfEdgeGraph::new();
    let a = graph.add_node_with_edges("a", &[10, 2, 3]);
    let b = graph.add_node_with_edges("b", &[20, 3, 4]);
    let c = graph.add_node_with_edges("c", &[30, 4, 2]);
    let d = graph.add_node_with_edges("d", &[20]);
    let e = graph.add_node_with_edges("e", &[30]);

    println!("Test {}", graph.dot());
    println!("{}", graph.degree(a));
    println!("{}", graph.degree(b));

    for (i, n) in &graph.neighbors {
        println!("{} {}", graph.edges[i], graph.edges[*n]);
    }

    let d = graph.merge_nodes(d, b, "bd");

    // for (i, n) in &graph.neighbors {
    //     println!("{} {}", graph.edges[i], graph.edges[*n]);
    // }

    println!("{}", graph.degree(c));
    println!("{}", graph.neighbors.len());

    println!("{}", graph.dot());

    let e = graph.merge_nodes(c, e, "ce");

    if graph.validate_neighbors() {
        println!("valid");
    } else {
        println!("invalid");
    }

    println!("{}", graph.dot());
    let f = graph.merge_nodes(d, e, "de");

    if graph.validate_neighbors() {
        println!("valid");
    } else {
        println!("invalid");
    }

    println!("{}", graph.dot());
    println!("{}", graph.node_labels());
    println!("{}", graph.degree(a));
    println!("{}", graph.neighbors.len());

    let g = graph.merge_nodes(a, f, "af");

    if graph.validate_neighbors() {
        println!("valid");
    } else {
        println!("invalid");
    }

    println!("{}", graph.dot());
    println!("{}", graph.neighbors.len());
    println!("{}", graph.degree(g));

    // println!("{}", graph.degree(b));
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TensorNetwork<T> {
    pub graph: HalfEdgeGraph<T, Slot>,
}

impl<T> TensorNetwork<T> {
    fn edge_to_min_degree_node(&self) -> Option<HedgeId> {
        let mut neighs = self.graph.reverse_nodemap.clone();
        if neighs.is_empty() {
            return None;
        }

        loop {
            let mut all_ext = true;
            for (node, initial) in &mut neighs {
                *initial = self.graph.neighbors[*initial];
                let start = self.graph.reverse_nodemap[node];

                if self.graph.involution[start] != start {
                    all_ext = false;
                    if *initial == start {
                        return Some(start);
                    }
                }
            }
            if all_ext {
                return None;
            }
        }
    }

    pub fn to_vec(&self) -> Vec<&T> {
        self.graph.nodes.values().collect()
    }
}

impl<N> TensorNetwork<MixedTensor<N>>
where
    N: Debug + TensorStructure,
{
    pub fn to_symbolic_tensor_vec(mut self) -> Vec<DataTensor<Atom, N>> {
        self.graph
            .nodes
            .drain()
            .filter(|(_, n)| n.is_symbolic())
            .map(|(_, n)| n.try_into_symbolic().unwrap())
            .collect()
    }
}

impl<T> From<Vec<T>> for TensorNetwork<T>
where
    T: TensorStructure,
{
    fn from(tensors: Vec<T>) -> Self {
        TensorNetwork {
            graph: Self::generate_network_graph(tensors),
        }
    }
}

impl<T> TensorNetwork<T>
where
    T: TensorStructure,
{
    pub fn new() -> Self {
        TensorNetwork {
            graph: HalfEdgeGraph::new(),
        }
    }

    pub fn push(&mut self, tensor: T) -> NodeId {
        let slots = tensor.external_structure().to_vec();
        self.graph.add_node_with_edges(tensor, &slots)
    }

    fn generate_network_graph(tensors: Vec<T>) -> HalfEdgeGraph<T, Slot> {
        let mut graph = HalfEdgeGraph::<T, Slot>::new();

        for tensor in tensors {
            let slots = tensor.external_structure().to_vec();
            graph.add_node_with_edges(tensor, &slots);
        }

        graph
    }

    pub fn edge_to_min_degree_node_with_depth(&self, depth: usize) -> Option<HedgeId>
    where
        T: TracksCount,
    {
        let mut min_degree = usize::MAX;
        let mut edge_to_min_degree_node = None;
        for (e, n) in &self.graph.nodemap {
            let edge_depth = self.graph.nodes[*n].contractions_num()
                + self.graph.nodes[self.graph.nodemap[self.graph.involution[e]]].contractions_num();

            if edge_depth < depth {
                let degree = self.graph.degree(*n);
                if degree < min_degree {
                    min_degree = degree;
                    if min_degree > 0 {
                        edge_to_min_degree_node = Some(self.graph.reverse_nodemap[*n]);
                    }
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
        self.graph.nodes.iter().next().unwrap().1.clone()
    }
}

impl<T> TensorNetwork<T> {
    pub fn dot(&self) -> std::string::String {
        // format!(
        //     "{:?}",
        //     Dot::with_attr_getters(
        //         &self.graph,
        //         &[Config::EdgeNoLabel, Config::NodeNoLabel],
        //         &|_, e| { format!("label=\"{}\"", e.weight()) },
        //         &|_, n| { format!("label=\"{}\"", n.1.structure()) }
        //     )
        // )
        // .into()
        self.graph.dot()
    }
}

impl<T> TensorNetwork<T>
where
    T: Debug + TensorStructure<Structure = HistoryStructure<Identifier>>,
{
    pub fn dotsym(&self, state: &State) -> std::string::String {
        // format!(
        //     "{:?}",
        //     Dot::with_attr_getters(
        //         &self.graph,
        //         &[Config::EdgeNoLabel, Config::NodeNoLabel],
        //         &|_, e| { format!("label=\"{}\"", e.weight()) },
        //         &|_, n| { format!("label=\"{}\"", n.1.structure().to_string(state)) }
        //     )
        // )
        // .into()
        self.graph.dot()
    }
}

impl<T> TensorNetwork<T>
where
    T: TensorStructure<Structure = HistoryStructure<Identifier>> + Clone,
{
    pub fn symbolic_shadow(
        &mut self,
        name: &str,
        state: &mut State,
        ws: &Workspace,
    ) -> TensorNetwork<MixedTensors> {
        for (i, n) in &mut self.graph.nodes {
            n.mut_structure().set_name(
                &state
                    .get_or_insert_fn(format!("{}{}", name, i.data().as_ffi()), None)
                    .unwrap(),
            );
        }

        let edges = self.graph.edges.clone();
        let involution = self.graph.involution.clone();
        let neighbors = self.graph.neighbors.clone();

        let mut nodes = DenseSlotMap::with_key();
        let mut nodemap = SecondaryMap::new();
        let mut reverse_nodemap = SecondaryMap::new();

        for (i, n) in &self.graph.nodes {
            let nid = nodes.insert(MixedTensor::<HistoryStructure<Identifier>>::from(
                n.structure().clone().shadow(state, ws).unwrap(),
            ));
            let mut first = true;
            for e in self.graph.edges_incident(i) {
                if first {
                    reverse_nodemap.insert(nid, e);
                    first = false;
                }
                nodemap.insert(e, nid);
            }
        }

        let g = HalfEdgeGraph {
            edges,
            involution,
            reverse_nodemap,
            neighbors,
            nodes,
            nodemap,
        };

        // let g: Graph<MixedTensors, Slot, Undirected> = Graph::map(
        //     &self.graph,
        //     |_, nw| {
        //         MixedTensor::<HistoryStructure<Identifier>>::from(
        //             nw.structure().clone().shadow(state, ws).unwrap(),
        //         )
        //     },
        //     |_, &w| w,
        // );
        TensorNetwork { graph: g }
    }
}

impl<T> TensorNetwork<T>
where
    T: HasName,
{
    pub fn name(&mut self, name: T::Name)
    where
        T::Name: From<std::string::String> + Display,
    {
        for (id, n) in &mut self.graph.nodes {
            n.set_name(&format!("{}{}", name, id.data().as_ffi()).into());
        }
    }
}

impl<T> TensorNetwork<T>
where
    T: HasName<Name = Identifier>,
{
    pub fn namesym(&mut self, name: &str, state: &mut State) {
        for (id, n) in &mut self.graph.nodes {
            n.set_name(
                &state
                    .get_or_insert_fn(format!("{}{}", name, id.data().as_ffi()), None)
                    .unwrap(),
            );
        }
    }
}

impl<T> TensorNetwork<T>
where
    T: Contract<T, LCM = T> + TensorStructure,
{
    pub fn contract_algo(&mut self, edge_choice: fn(&TensorNetwork<T>) -> Option<HedgeId>) {
        if let Some(e) = edge_choice(self) {
            self.contract_edge(e);
            // println!("{}", self.dot());
            self.contract_algo(edge_choice);
        }
    }
    fn contract_edge(&mut self, edge_idx: HedgeId) {
        let a = self.graph.nodemap[edge_idx];
        let b = self.graph.nodemap[self.graph.involution[edge_idx]];

        let ai = self.graph.nodes.get(a).unwrap();
        let bi = self.graph.nodes.get(b).unwrap();

        let f = ai.contract(bi).unwrap();

        self.graph.merge_nodes(a, b, f);
    }

    pub fn contract(&mut self) {
        self.contract_algo(|tn| tn.edge_to_min_degree_node())
    }
}

impl<T> TensorNetwork<T>
where
    T: SymbolicContract<T, LCM = T>
        + Debug
        + TensorStructure<Structure = HistoryStructure<Identifier>>
        + TracksCount,
{
    fn contract_edge_sym(&mut self, edge_idx: HedgeId, state: &State, ws: &Workspace) {
        let a = self.graph.nodemap[edge_idx];
        let b = self.graph.nodemap[self.graph.involution[edge_idx]];

        let ai = self.graph.nodes.get(a).unwrap();
        let bi = self.graph.nodes.get(b).unwrap();
        let f = ai.contract_sym(bi, state, ws).unwrap();

        self.graph.merge_nodes(a, b, f);
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
