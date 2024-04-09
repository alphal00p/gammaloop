use std::{hash::Hash, path::Display};

use ahash::{AHashMap, AHashSet};
use bitvec::vec::BitVec;
use indexmap::IndexMap;
use num::traits::{ops::inv, Inv};
use serde::{Deserialize, Serialize};

use crate::{
    graph::{Edge, EdgeType, Graph, Vertex},
    StabilityLevelSetting,
};

use bitvec::prelude::*;
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
struct NodeIndex(usize);

impl std::fmt::Display for NodeIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
struct LeafIndex(usize);

#[derive(Clone, Debug, Serialize, Deserialize)]
struct ForestNode<N, E> {
    parent: Option<NodeIndex>,
    label: N,
    content: NodeContent<E>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
enum NodeContent<L> {
    Internal(Vec<NodeIndex>),
    Leaves(Vec<L>),
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct Forest<N, E = LeafIndex> {
    nodes: AHashMap<NodeIndex, ForestNode<N, E>>,
    base_nodes: usize,
    base_leaves: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
enum ForestError {
    NodeNotFound,
    NodeIsInternal,
}

impl<N, E> Forest<N, E> {
    fn new() -> Self {
        Forest {
            nodes: AHashMap::new(),
            base_nodes: 0,
            base_leaves: 0,
        }
    }

    // Add a new leaf node to the forest
    fn add_leaf_node(&mut self, mut leaf_data: Vec<E>, node_data: N) -> NodeIndex {
        let data_index = self.base_leaves;

        let new_node = ForestNode {
            parent: None,
            label: node_data,
            content: NodeContent::Leaves(leaf_data),
        };
        let new_node_index = Self::leaf_node_index_from_usize(self.base_nodes);
        self.nodes.insert(new_node_index, new_node);
        self.base_nodes += 1;
        new_node_index
    }

    pub fn leaf_node_index_from_usize(i: usize) -> NodeIndex {
        NodeIndex(1 << i)
    }

    fn insert_leaf_at(&mut self, node_index: usize, leaf_data: E) -> Result<(), ForestError> {
        let node_id = Self::leaf_node_index_from_usize(node_index);
        if let Some(v) = self.nodes.get_mut(&node_id) {
            match &mut v.content {
                NodeContent::Leaves(leaf_index) => leaf_index.push(leaf_data),
                NodeContent::Internal(_) => {
                    return Err(ForestError::NodeIsInternal);
                }
            }
        } else {
            return Err(ForestError::NodeNotFound);
        }
        Ok(())
    }

    fn add_leaf_to_node(&mut self, node_index: NodeIndex, leaf_data: E) -> Result<(), ForestError> {
        if let Some(node) = self.nodes.get_mut(&node_index) {
            match &mut node.content {
                NodeContent::Leaves(leaf_index) => {
                    leaf_index.push(leaf_data);
                    Ok(())
                }
                NodeContent::Internal(_) => Err(ForestError::NodeIsInternal),
            }
        } else {
            Err(ForestError::NodeNotFound)
        }
    }

    fn collect_leaf_data(&self, node_index: NodeIndex) -> Vec<&E> {
        let mut leaf_data = Vec::new();
        self.dfs_collect_leaf_data(node_index, &mut leaf_data);
        leaf_data
    }

    // Helper for recursively collecting leaf indices.
    fn dfs_collect_leaf_data<'a>(&'a self, node_index: NodeIndex, leaf_data: &mut Vec<&'a E>) {
        match &self.nodes[&node_index].content {
            NodeContent::Leaves(leaf_index) => leaf_data.extend(leaf_index),
            NodeContent::Internal(children) => {
                for &child_index in children {
                    self.dfs_collect_leaf_data(child_index, leaf_data);
                }
            }
        }
    }

    // Access all leaves' data starting from a specific node.
    pub fn get_leaves_from_node(&self, node_index: NodeIndex) -> Vec<&E> {
        let leaf_indices = self.collect_leaf_data(node_index);
        leaf_indices
    }

    // Generalized union function for n-ary unions, works for any node indices.
    pub fn union(&mut self, node_indices: Vec<NodeIndex>, node_data: N) -> NodeIndex {
        let new_node = ForestNode {
            label: node_data,
            parent: None,
            content: NodeContent::Internal(node_indices.clone()),
        };

        let mut new_node_id = 0;
        for n in &node_indices {
            new_node_id += n.0;
        }

        let new_node_id = NodeIndex(new_node_id);

        self.nodes.insert(new_node_id, new_node);

        // Update parent references for each child node.
        for node_index in node_indices {
            if let Some(node) = self.nodes.get_mut(&node_index) {
                node.parent = Some(new_node_id);
            }
        }

        new_node_id
    }

    fn node_depth(&self, node_index: NodeIndex, current_depth: usize) -> usize {
        if let Some(node) = self.nodes.get(&node_index) {
            if let Some(parent_index) = node.parent {
                return self.node_depth(parent_index, current_depth + 1);
            }
        }
        current_depth
    }

    // Other methods as previously defined...
}

#[derive(Clone, Debug, Serialize, Deserialize)]
enum InvolutiveMapping<E> {
    Identity(Option<E>),
    Source((Option<E>, usize)),
    Sink(usize),
}

impl<E> InvolutiveMapping<E> {
    pub fn make_source(&mut self, sink: usize) -> Option<InvolutiveMapping<E>> {
        let data = match self {
            InvolutiveMapping::Identity(d) => d.take(),
            _ => None,
        };
        Some(InvolutiveMapping::Source((data, sink)))
    }

    pub fn new_identity(data: E) -> Self {
        InvolutiveMapping::Identity(Some(data))
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct Involution<N, E> {
    inv: Vec<(N, InvolutiveMapping<E>)>,
}

impl<N, E> Involution<N, E> {
    fn new() -> Self {
        Involution { inv: Vec::new() }
    }

    fn add_identity(&mut self, data: E, id: N) -> usize {
        let index = self.inv.len();
        self.inv.push((id, InvolutiveMapping::new_identity(data)));
        index
    }

    fn get_data(&self, index: usize) -> &E {
        match &self.inv[index].1 {
            InvolutiveMapping::Identity(data) => data.as_ref().unwrap(),
            InvolutiveMapping::Source((data, _)) => data.as_ref().unwrap(),
            InvolutiveMapping::Sink(i) => self.get_data(*i),
        }
    }

    pub fn get_node_id(&self, index: usize) -> &N {
        &self.inv[index].0
    }

    pub fn get_connected_node_id(&self, index: usize) -> &N {
        &self.inv[self.inv(index)].0
    }
    fn connect_to_identity(&mut self, source: usize, id: N) -> usize {
        let sink = self.inv.len();
        self.inv.push((id, InvolutiveMapping::Sink(source)));
        self.inv[source].1 = self.inv[source]
            .1
            .make_source(sink)
            .unwrap_or_else(|| panic!("Source must be an identity mapping"));
        sink
    }

    fn connect(&mut self, source: usize, sink: usize)
    where
        E: Clone,
    {
        self.inv[source].1 = self.inv[source].1.make_source(sink).unwrap();

        if let InvolutiveMapping::Identity(_) = &self.inv[sink].1 {
            self.inv[sink].1 = InvolutiveMapping::Sink(source);
        } else {
            panic!("Sink must be an identity mapping")
        }
    }

    fn add_pair(&mut self, data: E, source_id: N, sink_id: N) -> (usize, usize) {
        let source = self.add_identity(data, source_id);
        let sink = self.connect_to_identity(source, sink_id);
        (source, sink)
    }

    fn len(&self) -> usize {
        self.inv.len()
    }

    pub fn inv(&self, index: usize) -> usize {
        match &self.inv[index].1 {
            InvolutiveMapping::Identity(_) => index,
            InvolutiveMapping::Source((_, i)) => *i,
            InvolutiveMapping::Sink(i) => *i,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct NestingGraph<E, V> {
    nodes: IndexMap<NestingNode, V>, // Forest of nodes, that contain half-edges, with data E
    involution: Involution<NestingNode, E>, // Involution of half-edges
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
struct NestingNode {
    internalhedges: BitVec,
    externalhedges: BitVec,
}

impl NestingNode {
    pub fn new(len: usize) -> Self {
        NestingNode {
            internalhedges: bitvec![usize, Lsb0; 0; len],
            externalhedges: bitvec![usize, Lsb0; 0; len],
        }
    }

    pub fn from_builder<V>(builder: &NestingNodeBuilder<V>, len: usize) -> Self {
        let internalhedges = bitvec![usize, Lsb0; 0; len];
        let mut externalhedges = bitvec![usize, Lsb0; 0; len];

        for hedge in &builder.hedges {
            let mut bit = externalhedges.get_mut(*hedge).unwrap();
            *bit = true;
        }

        NestingNode {
            internalhedges,
            externalhedges,
        }
    }

    pub fn is_node(&self) -> bool {
        self.internalhedges.len() == 0
    }

    pub fn is_subgraph(&self) -> bool {
        self.internalhedges.len() > 0
    }
}

#[derive(Clone, Debug)]
struct NestingNodeBuilder<V> {
    data: V,
    hedges: Vec<usize>,
}

struct NestingGraphBuilder<E, V> {
    nodes: Vec<NestingNodeBuilder<V>>,
    involution: Involution<NodeIndex, E>,
}

impl<E, V> NestingGraphBuilder<E, V> {
    pub fn new() -> Self {
        NestingGraphBuilder {
            nodes: Vec::new(),
            involution: Involution::new(),
        }
    }

    pub fn add_node(&mut self, data: V) -> NodeIndex {
        let index = self.nodes.len();
        self.nodes.push(NestingNodeBuilder {
            data,
            hedges: Vec::new(),
        });
        NodeIndex(index)
    }

    pub fn add_edge(&mut self, source: NodeIndex, sink: NodeIndex, data: E) {
        let id = self.involution.add_pair(data, source, sink);
        self.nodes[source.0].hedges.push(id.0);

        self.nodes[sink.0].hedges.push(id.1);
    }

    pub fn add_external_edge(&mut self, source: NodeIndex, data: E) {
        let id = self.involution.add_identity(data, source);
        self.nodes[source.0].hedges.push(id);
    }
}

impl<E, V> From<NestingGraphBuilder<E, V>> for NestingGraph<E, V> {
    fn from(builder: NestingGraphBuilder<E, V>) -> Self {
        let len = builder.involution.len();
        let involution = Involution {
            inv: builder
                .involution
                .inv
                .into_iter()
                .map(|(n, i)| (NestingNode::from_builder(&builder.nodes[n.0], len), i))
                .collect(),
        };
        let nodes = builder
            .nodes
            .into_iter()
            .map(|x| (NestingNode::from_builder(&x, len), x.data))
            .collect();
        NestingGraph { nodes, involution }
    }
}

impl<E, V> NestingGraph<E, V> {
    // fn node_out_edges(&self, node_index: NodeIndex) -> Vec<usize> {
    //     let leaves = AHashSet::from_iter(
    //         self.nodes
    //             .get_leaves_from_node(node_index)
    //             .into_iter()
    //             .map(|x| *x),
    //     );

    //     let inv_leaves = leaves
    //         .clone()
    //         .into_iter()
    //         .map(|x| self.involutions.inv(x))
    //         .collect::<AHashSet<_>>();

    //     let externals: Vec<usize> = leaves.difference(&inv_leaves).map(|x| *x).collect();

    //     externals
    // }

    // fn out_degree(&self, node_index: NodeIndex) -> usize {
    //     self.node_out_edges(node_index).len()
    // }

    pub fn base_dot(&self) -> String {
        let mut out = "graph {".to_string();
        out.push_str("  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\";");
        for (i, (n, e)) in self.involution.inv.iter().enumerate() {
            match e {
                InvolutiveMapping::Identity(_) => {
                    out.push_str(&format!("ext{} [shape=none, label=\"\"];", i));
                    out.push_str(&format!(
                        "\n {} -- ext{};",
                        self.nodes.get_index_of(n).unwrap(),
                        i
                    ));
                }
                InvolutiveMapping::Source((_, _)) => {
                    out.push_str(&format!(
                        "  {} -- {};\n",
                        self.nodes.get_index_of(n).unwrap(),
                        self.nodes
                            .get_index_of(self.involution.get_connected_node_id(i))
                            .unwrap()
                    ));
                }
                InvolutiveMapping::Sink(_) => {}
            }
        }
        out += "}";
        out
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct UVEdge {}

impl UVEdge {
    pub fn from_edge(edge: &Edge) -> Self {
        UVEdge {}
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct UVNode {}

impl UVNode {
    pub fn from_vertex(vertex: &Vertex) -> Self {
        UVNode {}
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct UVGraph(NestingGraph<UVEdge, UVNode>);

impl UVGraph {
    pub fn from_graph(graph: &Graph) -> Self {
        let mut uv_graph = NestingGraphBuilder::new();

        for n in &graph.vertices {
            uv_graph.add_node(UVNode::from_vertex(n));
        }

        for edge in &graph.edges {
            let ves = edge.vertices;

            let sink = NodeIndex(ves[1]);
            let source = NodeIndex(ves[0]);

            match edge.edge_type {
                EdgeType::Virtual => {
                    uv_graph.add_edge(source, sink, UVEdge::from_edge(edge));
                }
                EdgeType::Outgoing => {
                    uv_graph.add_external_edge(source, UVEdge::from_edge(edge));
                }
                EdgeType::Incoming => {
                    uv_graph.add_external_edge(sink, UVEdge::from_edge(edge));
                }
            }
        }

        UVGraph(uv_graph.into())
    }
}

#[cfg(test)]
mod tests;
