use std::ops::{Index, IndexMut};

use linnet::half_edge::builder::HedgeGraphBuilder;
use linnet::half_edge::involution::{EdgeData, EdgeIndex, Flow, Hedge, HedgePair, Orientation};
use linnet::half_edge::nodestore::{NodeStorage, NodeStorageVec};
use linnet::half_edge::subgraph::{
    cut::OrientedCut, cycle::Cycle, HedgeNode, SuBitGraph, SubGraphOps, SubSetLike,
};
use linnet::half_edge::swap::Swap;
use linnet::half_edge::tree::SimpleTraversalTree;
use linnet::half_edge::{EdgeAccessors, HedgeGraph, HedgeGraphError, NodeIndex};
use linnet::permutation::Permutation;
use linnet::tree::{child_vec::ChildVecStore, Forest};
use pyo3::prelude::*;
use pyo3_stub_gen::derive::gen_stub_pyclass_enum;

use crate::graph::{EdgeRecord, HalfEdgeRecord, NodeRecord};

pub(crate) type VecNodeStore<V> = NodeStorageVec<V>;
pub(crate) type ForestNodeStore<V> = Forest<V, ChildVecStore<()>>;
type VecGraph = HedgeGraph<EdgeRecord, NodeRecord, HalfEdgeRecord, VecNodeStore<NodeRecord>>;
type ForestGraph = HedgeGraph<EdgeRecord, NodeRecord, HalfEdgeRecord, ForestNodeStore<NodeRecord>>;

/// Storage strategy used for graph nodes.
#[gen_stub_pyclass_enum]
#[pyclass(from_py_object, eq, eq_int, name = "NodeStore")]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PyNodeStore {
    /// Dense node records with one incidence bitset per node.
    #[default]
    Vec,
    /// Identification-preserving forest node records.
    Forest,
}

/// The two concrete Linnet graph representations exposed through one Python class.
pub(crate) enum PyHedgeGraph {
    Vec(VecGraph),
    Forest(ForestGraph),
}

#[derive(Clone)]
pub(crate) enum NativeNeighbors<'a> {
    Vec(<VecNodeStore<NodeRecord> as NodeStorage>::NeighborsIter<'a>),
    Forest(<ForestNodeStore<NodeRecord> as NodeStorage>::NeighborsIter<'a>),
}

impl Iterator for NativeNeighbors<'_> {
    type Item = Hedge;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Vec(iter) => iter.next(),
            Self::Forest(iter) => iter.next(),
        }
    }
}

impl ExactSizeIterator for NativeNeighbors<'_> {
    fn len(&self) -> usize {
        match self {
            Self::Vec(iter) => iter.len(),
            Self::Forest(iter) => iter.len(),
        }
    }
}

impl PyHedgeGraph {
    pub(crate) fn build(
        builder: HedgeGraphBuilder<EdgeRecord, NodeRecord, HalfEdgeRecord>,
        node_store: PyNodeStore,
    ) -> Self {
        match node_store {
            PyNodeStore::Vec => Self::Vec(builder.build()),
            PyNodeStore::Forest => Self::Forest(builder.build()),
        }
    }

    pub(crate) fn empty(node_store: PyNodeStore) -> Self {
        Self::build(HedgeGraphBuilder::new(), node_store)
    }

    pub(crate) fn node_store(&self) -> PyNodeStore {
        match self {
            Self::Vec(_) => PyNodeStore::Vec,
            Self::Forest(_) => PyNodeStore::Forest,
        }
    }

    pub(crate) fn into_node_store(self, node_store: PyNodeStore) -> Result<Self, HedgeGraphError> {
        match (self, node_store) {
            (graph @ Self::Vec(_), PyNodeStore::Vec)
            | (graph @ Self::Forest(_), PyNodeStore::Forest) => Ok(graph),
            (Self::Vec(graph), PyNodeStore::Forest) => graph
                .into_node_store::<ForestNodeStore<NodeRecord>>()
                .map(Self::Forest),
            (Self::Forest(graph), PyNodeStore::Vec) => graph
                .into_node_store::<VecNodeStore<NodeRecord>>()
                .map(Self::Vec),
        }
    }

    pub(crate) fn copy(&self, py: Python<'_>) -> PyResult<Self> {
        match self {
            Self::Vec(graph) => graph
                .map_data_ref_result(
                    |_, _, node| node.copy(py),
                    |_, _, _, edge| edge.map_result(|edge| edge.copy(py)),
                    |(_, half_edge)| half_edge.copy(py),
                )
                .map(Self::Vec),
            Self::Forest(graph) => graph
                .map_data_ref_result(
                    |_, _, node| node.copy(py),
                    |_, _, _, edge| edge.map_result(|edge| edge.copy(py)),
                    |(_, half_edge)| half_edge.copy(py),
                )
                .map(Self::Forest),
        }
    }

    pub(crate) fn n_nodes(&self) -> usize {
        match self {
            Self::Vec(graph) => graph.n_nodes(),
            Self::Forest(graph) => graph.n_nodes(),
        }
    }

    pub(crate) fn n_edges(&self) -> usize {
        match self {
            Self::Vec(graph) => graph.n_edges(),
            Self::Forest(graph) => graph.n_edges(),
        }
    }

    pub(crate) fn n_hedges(&self) -> usize {
        match self {
            Self::Vec(graph) => graph.n_hedges(),
            Self::Forest(graph) => graph.n_hedges(),
        }
    }

    pub(crate) fn iter_crown(&self, node: NodeIndex) -> NativeNeighbors<'_> {
        match self {
            Self::Vec(graph) => NativeNeighbors::Vec(graph.iter_crown(node)),
            Self::Forest(graph) => NativeNeighbors::Forest(graph.iter_crown(node)),
        }
    }

    pub(crate) fn iter_nodes(
        &self,
    ) -> impl Iterator<Item = (NodeIndex, NativeNeighbors<'_>, &NodeRecord)> {
        (0..self.n_nodes()).map(|index| {
            let node = NodeIndex(index);
            (node, self.iter_crown(node), &self[node])
        })
    }

    pub(crate) fn iter_node_ids(&self) -> impl Iterator<Item = NodeIndex> + '_ {
        (0..self.n_nodes()).map(NodeIndex)
    }

    pub(crate) fn iter_hedges(&self) -> impl Iterator<Item = (Hedge, &HalfEdgeRecord)> {
        (0..self.n_hedges()).map(|index| {
            let hedge = Hedge(index);
            (hedge, &self[hedge])
        })
    }

    pub(crate) fn iter_edges(
        &self,
    ) -> impl Iterator<Item = (HedgePair, EdgeIndex, EdgeData<&EdgeRecord>)> {
        (0..self.n_edges()).map(|index| {
            let edge = EdgeIndex(index);
            let pair = self[&edge].1;
            (
                pair,
                edge,
                EdgeData::new(&self[edge], self.orientation(edge)),
            )
        })
    }

    pub(crate) fn iter_edges_of<'a, S: SubSetLike + 'a>(
        &'a self,
        subgraph: &'a S,
    ) -> Box<dyn Iterator<Item = (HedgePair, EdgeIndex, EdgeData<&'a EdgeRecord>)> + 'a> {
        match self {
            Self::Vec(graph) => Box::new(graph.iter_edges_of(subgraph)),
            Self::Forest(graph) => Box::new(graph.iter_edges_of(subgraph)),
        }
    }

    pub(crate) fn iter_nodes_of<'a, S: SubSetLike + 'a>(
        &'a self,
        subgraph: &'a S,
    ) -> Box<dyn Iterator<Item = (NodeIndex, NativeNeighbors<'a>, &'a NodeRecord)> + 'a> {
        match self {
            Self::Vec(graph) => Box::new(
                graph
                    .iter_nodes_of(subgraph)
                    .map(|(node, crown, data)| (node, NativeNeighbors::Vec(crown), data)),
            ),
            Self::Forest(graph) => Box::new(
                graph
                    .iter_nodes_of(subgraph)
                    .map(|(node, crown, data)| (node, NativeNeighbors::Forest(crown), data)),
            ),
        }
    }

    pub(crate) fn flow(&self, hedge: Hedge) -> Flow {
        match self {
            Self::Vec(graph) => graph.flow(hedge),
            Self::Forest(graph) => graph.flow(hedge),
        }
    }

    pub(crate) fn set_flow(&mut self, hedge: Hedge, flow: Flow) {
        match self {
            Self::Vec(graph) => graph.set_flow(hedge, flow),
            Self::Forest(graph) => graph.set_flow(hedge, flow),
        }
    }

    pub(crate) fn orientation(&self, edge: EdgeIndex) -> Orientation {
        match self {
            Self::Vec(graph) => graph.orientation(edge),
            Self::Forest(graph) => graph.orientation(edge),
        }
    }

    pub(crate) fn set_orientation(&mut self, edge: EdgeIndex, orientation: Orientation) {
        match self {
            Self::Vec(graph) => graph.set_orientation(edge, orientation),
            Self::Forest(graph) => graph.set_orientation(edge, orientation),
        }
    }

    pub(crate) fn node_id(&self, hedge: Hedge) -> NodeIndex {
        match self {
            Self::Vec(graph) => graph.node_id(hedge),
            Self::Forest(graph) => graph.node_id(hedge),
        }
    }

    pub(crate) fn inv(&self, hedge: Hedge) -> Hedge {
        match self {
            Self::Vec(graph) => graph.inv(hedge),
            Self::Forest(graph) => graph.inv(hedge),
        }
    }

    pub(crate) fn full_filter(&self) -> SuBitGraph {
        match self {
            Self::Vec(graph) => graph.full_filter(),
            Self::Forest(graph) => graph.full_filter(),
        }
    }

    pub(crate) fn complement(&self, subgraph: &SuBitGraph) -> SuBitGraph {
        match self {
            Self::Vec(graph) => subgraph.complement(graph),
            Self::Forest(graph) => subgraph.complement(graph),
        }
    }

    pub(crate) fn connected_components(&self, subgraph: &SuBitGraph) -> Vec<SuBitGraph> {
        match self {
            Self::Vec(graph) => graph.connected_components(subgraph),
            Self::Forest(graph) => graph.connected_components(subgraph),
        }
    }

    pub(crate) fn count_connected_components(&self, subgraph: &SuBitGraph) -> usize {
        match self {
            Self::Vec(graph) => graph.count_connected_components(subgraph),
            Self::Forest(graph) => graph.count_connected_components(subgraph),
        }
    }

    pub(crate) fn cyclotomatic_number(&self, subgraph: &SuBitGraph) -> usize {
        match self {
            Self::Vec(graph) => graph.cyclotomatic_number(subgraph),
            Self::Forest(graph) => graph.cyclotomatic_number(subgraph),
        }
    }

    pub(crate) fn bridges_of(&self, subgraph: &SuBitGraph) -> SuBitGraph {
        match self {
            Self::Vec(graph) => graph.bridges_of(subgraph),
            Self::Forest(graph) => graph.bridges_of(subgraph),
        }
    }

    pub(crate) fn cycle_basis_of(&self, subgraph: &SuBitGraph) -> (Vec<Cycle>, SuBitGraph) {
        match self {
            Self::Vec(graph) => graph.cycle_basis_of(subgraph),
            Self::Forest(graph) => graph.cycle_basis_of(subgraph),
        }
    }

    pub(crate) fn all_spanning_forests_of(&self, subgraph: &SuBitGraph) -> Vec<SuBitGraph> {
        match self {
            Self::Vec(graph) => graph.all_spanning_forests_of(subgraph),
            Self::Forest(graph) => graph.all_spanning_forests_of(subgraph),
        }
    }

    pub(crate) fn combine_to_single_hedgenode(&self, nodes: &[NodeIndex]) -> HedgeNode {
        match self {
            Self::Vec(graph) => graph.combine_to_single_hedgenode(nodes),
            Self::Forest(graph) => graph.combine_to_single_hedgenode(nodes),
        }
    }

    pub(crate) fn all_cuts(
        &self,
        source: HedgeNode,
        target: HedgeNode,
    ) -> Vec<(SuBitGraph, OrientedCut, SuBitGraph)> {
        match self {
            Self::Vec(graph) => graph.all_cuts(source, target),
            Self::Forest(graph) => graph.all_cuts(source, target),
        }
    }

    pub(crate) fn empty_traversal_tree(&self) -> SimpleTraversalTree {
        match self {
            Self::Vec(graph) => SimpleTraversalTree::empty(graph),
            Self::Forest(graph) => SimpleTraversalTree::empty(graph),
        }
    }

    pub(crate) fn depth_first_traverse(
        &self,
        subgraph: &SuBitGraph,
        root: &NodeIndex,
        include: Option<Hedge>,
    ) -> Result<SimpleTraversalTree, HedgeGraphError> {
        match self {
            Self::Vec(graph) => {
                SimpleTraversalTree::depth_first_traverse(graph, subgraph, root, include)
            }
            Self::Forest(graph) => {
                SimpleTraversalTree::depth_first_traverse(graph, subgraph, root, include)
            }
        }
    }

    pub(crate) fn breadth_first_traverse(
        &self,
        subgraph: &SuBitGraph,
        root: &NodeIndex,
        include: Option<Hedge>,
    ) -> Result<SimpleTraversalTree, HedgeGraphError> {
        match self {
            Self::Vec(graph) => {
                SimpleTraversalTree::breadth_first_traverse(graph, subgraph, root, include)
            }
            Self::Forest(graph) => {
                SimpleTraversalTree::breadth_first_traverse(graph, subgraph, root, include)
            }
        }
    }

    pub(crate) fn extract_nodes(
        &mut self,
        nodes: impl IntoIterator<Item = NodeIndex>,
        split_edge: impl FnMut(EdgeData<&EdgeRecord>) -> EdgeData<EdgeRecord>,
        internal_edge: impl FnMut(EdgeData<EdgeRecord>) -> EdgeData<EdgeRecord>,
    ) -> Self {
        let nodes = nodes.into_iter().collect::<Vec<_>>();
        match self {
            Self::Vec(graph) => Self::Vec(graph.extract_nodes(nodes, split_edge, internal_edge)),
            Self::Forest(graph) => {
                Self::Forest(graph.extract_nodes(nodes, split_edge, internal_edge))
            }
        }
    }

    pub(crate) fn extract(
        &mut self,
        subgraph: &SuBitGraph,
        split_edge: impl FnMut(EdgeData<&EdgeRecord>) -> EdgeData<EdgeRecord>,
        internal_edge: impl FnMut(EdgeData<EdgeRecord>) -> EdgeData<EdgeRecord>,
        split_node: impl FnMut(&NodeRecord) -> NodeRecord,
        internal_node: impl FnMut(NodeRecord) -> NodeRecord,
    ) -> Self {
        match self {
            Self::Vec(graph) => Self::Vec(graph.extract(
                subgraph,
                split_edge,
                internal_edge,
                split_node,
                internal_node,
            )),
            Self::Forest(graph) => Self::Forest(graph.extract(
                subgraph,
                split_edge,
                internal_edge,
                split_node,
                internal_node,
            )),
        }
    }

    pub(crate) fn delete_hedges(&mut self, subgraph: &SuBitGraph) {
        match self {
            Self::Vec(graph) => graph.delete_hedges(subgraph),
            Self::Forest(graph) => graph.delete_hedges(subgraph),
        }
    }

    pub(crate) fn contract_subgraph(&mut self, subgraph: &SuBitGraph, node: NodeRecord) {
        match self {
            Self::Vec(graph) => graph.contract_subgraph(subgraph, node),
            Self::Forest(graph) => graph.contract_subgraph(subgraph, node),
        }
    }

    pub(crate) fn append_disconnected_mut(
        &mut self,
        other: Self,
    ) -> Result<Hedge, HedgeGraphError> {
        let other = other.into_node_store(self.node_store())?;
        match (self, other) {
            (Self::Vec(left), Self::Vec(right)) => left.append_disconnected_mut(right),
            (Self::Forest(left), Self::Forest(right)) => left.append_disconnected_mut(right),
            _ => unreachable!("node stores were normalized"),
        }
    }

    pub(crate) fn join_with_hedge_data(
        self,
        other: Self,
        matching: impl Fn(
            Hedge,
            Flow,
            EdgeData<&EdgeRecord>,
            &HalfEdgeRecord,
            Hedge,
            Flow,
            EdgeData<&EdgeRecord>,
            &HalfEdgeRecord,
        ) -> bool,
        merge: impl Fn(
            Flow,
            EdgeData<EdgeRecord>,
            Flow,
            EdgeData<EdgeRecord>,
        ) -> (Flow, EdgeData<EdgeRecord>),
    ) -> Result<Self, HedgeGraphError> {
        let node_store = self.node_store();
        let other = other.into_node_store(node_store)?;
        match (self, other) {
            (Self::Vec(left), Self::Vec(right)) => left
                .join_with_hedge_data(right, matching, merge)
                .map(Self::Vec),
            (Self::Forest(left), Self::Forest(right)) => left
                .join_with_hedge_data(right, matching, merge)
                .map(Self::Forest),
            _ => unreachable!("node stores were normalized"),
        }
    }
}

macro_rules! index_native {
    ($index:ty, $output:ty) => {
        impl Index<$index> for PyHedgeGraph {
            type Output = $output;

            fn index(&self, index: $index) -> &Self::Output {
                match self {
                    Self::Vec(graph) => &graph[index],
                    Self::Forest(graph) => &graph[index],
                }
            }
        }

        impl IndexMut<$index> for PyHedgeGraph {
            fn index_mut(&mut self, index: $index) -> &mut Self::Output {
                match self {
                    Self::Vec(graph) => &mut graph[index],
                    Self::Forest(graph) => &mut graph[index],
                }
            }
        }
    };
}

index_native!(NodeIndex, NodeRecord);
index_native!(EdgeIndex, EdgeRecord);
index_native!(Hedge, HalfEdgeRecord);

impl Index<&EdgeIndex> for PyHedgeGraph {
    type Output = (EdgeRecord, HedgePair);

    fn index(&self, index: &EdgeIndex) -> &Self::Output {
        match self {
            Self::Vec(graph) => &graph[index],
            Self::Forest(graph) => &graph[index],
        }
    }
}

impl Index<&Hedge> for PyHedgeGraph {
    type Output = EdgeIndex;

    fn index(&self, index: &Hedge) -> &Self::Output {
        match self {
            Self::Vec(graph) => &graph[index],
            Self::Forest(graph) => &graph[index],
        }
    }
}

macro_rules! swap_native {
    ($index:ty) => {
        impl Swap<$index> for PyHedgeGraph {
            fn swap(&mut self, left: $index, right: $index) {
                match self {
                    Self::Vec(graph) => graph.swap(left, right),
                    Self::Forest(graph) => graph.swap(left, right),
                }
            }

            fn len(&self) -> $index {
                match self {
                    Self::Vec(graph) => <VecGraph as Swap<$index>>::len(graph),
                    Self::Forest(graph) => <ForestGraph as Swap<$index>>::len(graph),
                }
            }

            fn is_zero_length(&self) -> bool {
                match self {
                    Self::Vec(graph) => <VecGraph as Swap<$index>>::is_zero_length(graph),
                    Self::Forest(graph) => <ForestGraph as Swap<$index>>::is_zero_length(graph),
                }
            }

            fn permute(&mut self, permutation: &Permutation) {
                match self {
                    Self::Vec(graph) => <VecGraph as Swap<$index>>::permute(graph, permutation),
                    Self::Forest(graph) => {
                        <ForestGraph as Swap<$index>>::permute(graph, permutation)
                    }
                }
            }
        }
    };
}

swap_native!(NodeIndex);
swap_native!(EdgeIndex);
swap_native!(Hedge);
