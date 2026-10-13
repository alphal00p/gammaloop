use super::{
    hedgevec::SmartEdgeVec,
    involution::{
        EdgeData, Flow, Hedge, HedgePair, HedgeVec, Involution, InvolutionError, Orientation,
    },
    nodestore::{NodeStore, Topology},
    subgraph::{ModifySubSet, SuBitGraph, SubSetLike},
    swap::Swap,
    HedgeGraph, NoData, NodeIndex,
};

#[derive(Clone, Debug)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct HedgeData<H> {
    pub data: H,
    pub is_in_subgraph: bool,
    pub node: NodeIndex,
}

impl<H: Default> From<NodeIndex> for HedgeData<H> {
    fn from(node: NodeIndex) -> Self {
        HedgeData {
            data: H::default(),
            is_in_subgraph: false,
            node,
        }
    }
}

impl<H> HedgeData<H> {
    pub fn set_in_subgraph(mut self, is_in_subgraph: bool) -> Self {
        self.is_in_subgraph = is_in_subgraph;
        self
    }

    pub fn map<F, T>(self, f: F) -> HedgeData<T>
    where
        F: FnOnce(H) -> T,
    {
        HedgeData {
            data: f(self.data),
            is_in_subgraph: self.is_in_subgraph,
            node: self.node,
        }
    }

    pub fn map_result<F, T, E>(self, f: F) -> Result<HedgeData<T>, E>
    where
        F: FnOnce(H) -> Result<T, E>,
    {
        match f(self.data) {
            Ok(data) => Ok(HedgeData {
                data,
                is_in_subgraph: self.is_in_subgraph,
                node: self.node,
            }),
            Err(err) => Err(err),
        }
    }
}

#[derive(Clone, Debug)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
/// A builder for programmatically constructing [`HedgeGraph`] instances.
///
/// This builder allows for the incremental addition of nodes and edges (both
/// paired and external/dangling) before finalizing the graph structure.
///
/// # Type Parameters
///
/// - `E`: The type of data to be associated with edges.
/// - `V`: The type of data to be associated with nodes.
///
/// # Example
///
/// ```rust,ignore
/// use linnet::half_edge::builder::HedgeGraphBuilder;
/// use linnet::half_edge::involution::{Flow, Orientation};
/// use linnet::half_edge::HedgeGraph;
///
/// let mut builder = HedgeGraphBuilder::<&str, &str>::new();
///
/// let node0 = builder.add_node("Node_0_Data");
/// let node1 = builder.add_node("Node_1_Data");
/// let node2 = builder.add_node("Node_2_Data");
///
/// builder.add_edge(node0, node1, "Edge_01_Data", Orientation::Default);
/// builder.add_external_edge(node2, "External_Edge_Data", Orientation::Undirected, Flow::Source);
///
/// let graph: HedgeGraph<&str, &str> = builder.build();
/// ```
pub struct HedgeGraphBuilder<E, V, H = NoData> {
    hedge_data: HedgeVec<H>,
    nodes: Vec<V>,
    /// Record the partition directly in half-edge order, without a temporary
    /// incidence vector for each node and a second pass to flatten it.
    node_labels: Vec<NodeIndex>,
    /// The [`Involution`] structure managing the half-edges being added to the graph.
    pub(crate) involution: Involution<E>,
}

impl<E, V, H> HedgeGraphBuilder<E, V, H> {
    pub fn new() -> Self {
        HedgeGraphBuilder {
            hedge_data: HedgeVec::new(),
            nodes: Vec::new(),
            node_labels: Vec::new(),
            involution: Involution::new(),
        }
    }

    /// An empty builder with room for `nodes` nodes and `hedges` half-edges.
    pub fn with_capacity(nodes: usize, hedges: usize) -> Self {
        HedgeGraphBuilder {
            hedge_data: HedgeVec::with_capacity(hedges),
            nodes: Vec::with_capacity(nodes),
            node_labels: Vec::with_capacity(hedges),
            involution: Involution {
                inv: HedgeVec::with_capacity(hedges),
            },
        }
    }

    /// Retain nodes in the requested order, preserving each node's hedge order.
    ///
    /// This compacts construction state before a node store is built. Edges
    /// crossing the retained boundary become dangling with their existing flow
    /// and orientation, as in [`Involution::delete`]. The returned map translates
    /// old hedge IDs to new IDs; omitted hedges map to `None`.
    ///
    /// # Panics
    /// Panics if a requested node is out of range or appears more than once.
    pub fn retain_nodes_in_order(mut self, order: &[NodeIndex]) -> (Self, HedgeVec<Option<Hedge>>) {
        let mut used = vec![false; self.nodes.len()];
        let mut incident = vec![Vec::new(); self.nodes.len()];
        for (hedge, node) in self.node_labels.iter().enumerate() {
            incident[node.0].push(Hedge(hedge));
        }
        let mut hedges = Vec::new();
        for &node in order {
            assert!(!used[node.0], "retained nodes must be unique");
            used[node.0] = true;
            hedges.extend_from_slice(&incident[node.0]);
        }
        let len = self.involution.len().0;
        let mut mapping: HedgeVec<_> = std::iter::repeat_n(None, len).collect();
        let mut at = (0..len).collect::<Vec<_>>();
        let mut positions = at.clone();
        for (position, &hedge) in hedges.iter().enumerate() {
            mapping[hedge] = Some(Hedge(position));
            let from = positions[hedge.0];
            if from != position {
                self.involution.swap(Hedge(position), Hedge(from));
                self.hedge_data.swap(Hedge(position), Hedge(from));
                self.node_labels.swap(position, from);
                positions.swap(at[position], at[from]);
                at.swap(position, from);
            }
        }
        let mut discarded = SuBitGraph::empty(len);
        for index in hedges.len()..len {
            discarded.add(Hedge(index));
        }
        self.involution.delete(&discarded);
        self.hedge_data.truncate(hedges.len());
        self.node_labels.truncate(hedges.len());
        let mut node_map = vec![NodeIndex(0); self.nodes.len()];
        for (new, old) in order.iter().enumerate() {
            node_map[old.0] = NodeIndex(new);
        }
        for node in &mut self.node_labels {
            *node = node_map[node.0];
        }
        let mut nodes = self.nodes.into_iter().map(Some).collect::<Vec<_>>();
        self.nodes = order
            .iter()
            .map(|node| nodes[node.0].take().expect("unique retained node"))
            .collect();
        (self, mapping)
    }

    pub fn build<N: Topology>(self) -> HedgeGraph<E, V, H, N> {
        self.into()
    }

    pub fn build_with_map<V2, N: Topology>(
        self,
        map: impl FnMut(V) -> V2,
    ) -> HedgeGraph<E, V2, H, N> {
        HedgeGraph {
            node_store: NodeStore::new(
                N::from_labels(self.node_labels.into_iter(), self.nodes.len()),
                self.nodes.into_iter().map(map).collect(),
            ),
            edge_store: SmartEdgeVec::new(self.involution),
            hedge_data: self.hedge_data,
        }
    }

    /// Connect two dangling hedges before materializing edge storage.
    ///
    /// Node incidence and hedge payloads keep their order. Merge ordering and
    /// invalid-input behavior are owned by [`Involution::connect_identities`].
    pub fn connect_identities(
        &mut self,
        source: Hedge,
        sink: Hedge,
        merge: impl FnOnce(Flow, EdgeData<E>, Flow, EdgeData<E>) -> (Flow, EdgeData<E>),
    ) -> Result<HedgePair, InvolutionError>
    where
        E: Clone,
    {
        self.involution.connect_identities(source, sink, merge)
    }

    pub fn add_node(&mut self, data: V) -> NodeIndex {
        let index = self.nodes.len();
        self.nodes.push(data);
        NodeIndex(index)
    }

    pub fn add_edge(
        &mut self,
        source: impl Into<HedgeData<H>>,
        sink: impl Into<HedgeData<H>>,
        data: E,
        directed: impl Into<Orientation>,
    ) {
        let source = source.into();
        let sink = sink.into();
        assert!(source.node.0 < self.nodes.len() && sink.node.0 < self.nodes.len());
        self.involution.add_pair(data, directed);
        self.hedge_data.push(source.data);
        self.hedge_data.push(sink.data);

        self.node_labels.extend([source.node, sink.node]);
    }

    pub fn add_external_edge(
        &mut self,
        source: impl Into<HedgeData<H>>,
        data: E,
        orientation: impl Into<Orientation>,
        underlying: Flow,
    ) {
        let source = source.into();
        assert!(source.node.0 < self.nodes.len());
        self.involution.add_identity(data, orientation, underlying);
        self.node_labels.push(source.node);
        self.hedge_data.push(source.data);
    }
}

impl<E, V, H> Default for HedgeGraphBuilder<E, V, H> {
    fn default() -> Self {
        Self::new()
    }
}

impl<E, V, H, N: Topology> From<HedgeGraphBuilder<E, V, H>> for HedgeGraph<E, V, H, N> {
    fn from(builder: HedgeGraphBuilder<E, V, H>) -> Self {
        builder.build_with_map(std::convert::identity)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::half_edge::nodestore::{reference::RefTopology, NodeForest};

    fn retained<N: Topology>() {
        let mut builder = HedgeGraphBuilder::<usize, usize>::new();
        let nodes = (0..4).map(|i| builder.add_node(i)).collect::<Vec<_>>();
        builder.add_edge(nodes[0], nodes[1], 10, Orientation::Default);
        builder.add_edge(nodes[1], nodes[1], 11, Orientation::Reversed);
        builder.add_external_edge(nodes[2], 12, Orientation::Undirected, Flow::Sink);
        let (builder, first) = builder.retain_nodes_in_order(&[nodes[1], nodes[0], nodes[3]]);
        assert_eq!(
            first
                .into_iter()
                .map(|(_, value)| value)
                .collect::<Vec<_>>(),
            vec![
                Some(Hedge(3)),
                Some(Hedge(0)),
                Some(Hedge(1)),
                Some(Hedge(2)),
                None
            ]
        );
        let (builder, second) = builder.retain_nodes_in_order(&[NodeIndex(0), NodeIndex(2)]);
        assert_eq!(
            second
                .into_iter()
                .map(|(_, value)| value)
                .collect::<Vec<_>>(),
            vec![Some(Hedge(0)), Some(Hedge(1)), Some(Hedge(2)), None]
        );
        let graph: HedgeGraph<usize, usize, NoData, N> = builder.build_with_map(|i| i + 100);
        assert_eq!(graph.node_bound(), 2);
        assert_eq!(
            graph
                .iter_nodes()
                .map(|(_, crown, data)| (crown.collect::<Vec<_>>(), *data))
                .collect::<Vec<_>>(),
            vec![(vec![Hedge(0), Hedge(1), Hedge(2)], 101), (vec![], 103)]
        );
        let edges = graph
            .iter_edges()
            .map(|(pair, _, data)| (pair, *data.data, data.orientation))
            .collect::<Vec<_>>();
        assert_eq!(
            edges,
            vec![
                (
                    HedgePair::Unpaired {
                        hedge: Hedge(0),
                        flow: Flow::Sink
                    },
                    10,
                    Orientation::Default
                ),
                (
                    HedgePair::Paired {
                        source: Hedge(1),
                        sink: Hedge(2)
                    },
                    11,
                    Orientation::Reversed
                )
            ]
        );
    }

    #[test]
    fn flat_incidence_retains_exact_frames_in_each_topology() {
        retained::<NodeForest>();
        retained::<RefTopology>();
    }
}
