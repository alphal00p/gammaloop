use ahash::{AHashMap, AHashSet};
use num::traits::Inv;

use crate::graph::{Edge, Graph, Vertex};

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct NodeIndex(usize);

#[derive(Clone, Copy, Debug)]
struct LeafIndex(usize);

#[derive(Clone)]
struct ForestNode<N, E> {
    parent: Option<NodeIndex>,
    label: N,
    content: NodeContent<E>,
}

#[derive(Clone)]
enum NodeContent<L> {
    Internal(Vec<NodeIndex>),
    Leaves(Vec<L>),
}

struct Forest<N, E = LeafIndex> {
    nodes: AHashMap<NodeIndex, ForestNode<N, E>>,
    base_nodes: usize,
    base_leaves: usize,
}

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

struct Involution<E> {
    inv: Vec<InvolutiveMapping<E>>,
}

impl<E> Involution<E> {
    fn new() -> Self {
        Involution { inv: Vec::new() }
    }

    fn add_identity(&mut self, data: E) -> usize {
        let index = self.inv.len();
        self.inv.push(InvolutiveMapping::new_identity(data));
        index
    }

    fn get_data(&self, index: usize) -> &E {
        match &self.inv[index] {
            InvolutiveMapping::Identity(data) => data.as_ref().unwrap(),
            InvolutiveMapping::Source((data, _)) => data.as_ref().unwrap(),
            InvolutiveMapping::Sink(i) => self.get_data(*i),
        }
    }

    fn connect_to_identity(&mut self, source: usize) -> usize {
        let sink = self.inv.len();
        self.inv.push(InvolutiveMapping::Sink(source));
        self.inv[source] = self.inv[source]
            .make_source(sink)
            .unwrap_or_else(|| panic!("Source must be an identity mapping"));
        sink
    }

    fn connect(&mut self, source: usize, sink: usize)
    where
        E: Clone,
    {
        self.inv[source] = self.inv[source].make_source(sink).unwrap();

        if let InvolutiveMapping::Identity(_) = &self.inv[sink] {
            self.inv[sink] = InvolutiveMapping::Sink(source);
        } else {
            panic!("Sink must be an identity mapping")
        }
    }

    fn add_pair(&mut self, data: E) -> (usize, usize) {
        let source = self.add_identity(data);
        let sink = self.connect_to_identity(source);
        (source, sink)
    }

    fn len(&self) -> usize {
        self.inv.len()
    }

    pub fn inv(&self, index: usize) -> usize {
        match &self.inv[index] {
            InvolutiveMapping::Identity(_) => index,
            InvolutiveMapping::Source((_, i)) => *i,
            InvolutiveMapping::Sink(i) => *i,
        }
    }
}

struct NestingGraph<E, V> {
    nodes: Forest<V, usize>, // Forest of nodes, that contain half-edges, with data E
    involutions: Involution<E>, // Involution of half-edges
}

enum NestingGraphError {
    NodeNotFound,
    NodeIsNotLeaf,
    SourceNodeNotFound,
    SinkNodeNotFound,
}

impl From<ForestError> for NestingGraphError {
    fn from(e: ForestError) -> Self {
        match e {
            ForestError::NodeNotFound => NestingGraphError::NodeNotFound,
            ForestError::NodeIsInternal => NestingGraphError::NodeIsNotLeaf,
        }
    }
}

impl<E, V> NestingGraph<E, V> {
    fn new() -> Self {
        NestingGraph {
            nodes: Forest::new(),
            involutions: Involution::new(),
        }
    }

    fn add_node(&mut self, data: V) -> NodeIndex {
        self.nodes.add_leaf_node(Vec::new(), data)
    }

    fn add_edge(
        &mut self,
        source: NodeIndex,
        sink: NodeIndex,
        data: E,
    ) -> Result<(), NestingGraphError> {
        let id = (self.involutions.len(), self.involutions.len() + 1);
        self.nodes
            .add_leaf_to_node(source, id.0)
            .map_err(|e| match e {
                ForestError::NodeNotFound => NestingGraphError::SourceNodeNotFound,
                ForestError::NodeIsInternal => NestingGraphError::NodeIsNotLeaf,
            })?;

        self.nodes
            .add_leaf_to_node(sink, id.1)
            .map_err(|e| match e {
                ForestError::NodeNotFound => NestingGraphError::SourceNodeNotFound,
                ForestError::NodeIsInternal => NestingGraphError::NodeIsNotLeaf,
            })?;

        self.involutions.add_pair(data);
        Ok(())
    }

    fn add_external_edge(&mut self, source: NodeIndex, data: E) {
        let id = self.involutions.add_identity(data);
        self.nodes.add_leaf_to_node(source, id);
    }

    fn node_out_edges(&self, node_index: NodeIndex) -> Vec<usize> {
        let leaves = AHashSet::from_iter(
            self.nodes
                .get_leaves_from_node(node_index)
                .into_iter()
                .map(|x| *x),
        );

        let inv_leaves = leaves
            .clone()
            .into_iter()
            .map(|x| self.involutions.inv(x))
            .collect::<AHashSet<_>>();

        let externals: Vec<usize> = leaves.difference(&inv_leaves).map(|x| *x).collect();

        externals
    }

    fn out_degree(&self, node_index: NodeIndex) -> usize {
        self.node_out_edges(node_index).len()
    }
}

struct UVEdge {}

impl UVEdge {
    pub fn from_edge(edge: &Edge) -> Self {
        UVEdge {}
    }
}

struct UVNode {}

impl UVNode {
    pub fn from_vertex(vertex: &Vertex) -> Self {
        UVNode {}
    }
}

struct UVGraph(NestingGraph<UVEdge, UVNode>);

impl UVGraph {
    fn from_graph(graph: &Graph) -> Self {
        let mut uv_graph = UVGraph(NestingGraph::new());

        for edge in &graph.edges {
            let ves = edge.vertices;

            let sink = Forest::<(), ()>::leaf_node_index_from_usize(ves[0]);
            let source = Forest::<(), ()>::leaf_node_index_from_usize(ves[1]);

            let mut edge_exists = false;

            while !edge_exists {
                match uv_graph.0.add_edge(source, sink, UVEdge::from_edge(edge)) {
                    Ok(_) => edge_exists = true,
                    Err(NestingGraphError::SourceNodeNotFound) => {
                        uv_graph
                            .0
                            .add_node(UVNode::from_vertex(&graph.vertices[ves[1]]));
                    }
                    Err(NestingGraphError::SinkNodeNotFound) => {
                        uv_graph
                            .0
                            .add_node(UVNode::from_vertex(&graph.vertices[ves[0]]));
                    }
                    Err(_) => panic!("Unexpected error"),
                }
            }
        }

        uv_graph
    }
}
