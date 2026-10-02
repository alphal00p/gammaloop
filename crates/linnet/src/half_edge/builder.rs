use super::{
    hedgevec::SmartEdgeVec,
    involution::{
        EdgeData, Flow, Hedge, HedgePair, HedgeVec, Involution, InvolutionError, Orientation,
    },
    nodestore::NodeStorageOps,
    subgraph::{BaseSubgraph, ModifySubSet, SuBitGraph, SubSetLike},
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
/// A temporary structure used during the construction of a [`HedgeGraph`].
///
/// It holds the data for a node and a list of half-edges that are incident to it
/// before the full graph topology (e.g., specific `NodeStorage` format) is finalized.
pub struct HedgeNodeBuilder<V> {
    /// The data associated with the node being built.
    pub(crate) data: V,
    /// A list of [`Hedge`] identifiers that are incident to this node.
    pub(crate) hedges: Vec<Hedge>,
}

impl<V> HedgeNodeBuilder<V> {
    pub fn to_base<S: BaseSubgraph>(&self, len: usize) -> S {
        let mut subgraph = S::empty(len);

        for hedge in &self.hedges {
            subgraph.add(*hedge);
        }

        subgraph
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
/// // Assuming a NodeStorage type MyNodeStore is defined and implements NodeStorageOps
/// // let graph: HedgeGraph<&str, &str, MyNodeStore> = builder.build();
/// ```
pub struct HedgeGraphBuilder<E, V, H = NoData> {
    hedge_data: HedgeVec<H>,
    /// A list of nodes currently being built, stored as [`HedgeNodeBuilder`] instances.
    nodes: Vec<HedgeNodeBuilder<V>>,
    /// The [`Involution`] structure managing the half-edges being added to the graph.
    pub(crate) involution: Involution<E>,
}

impl<E, V, H> HedgeGraphBuilder<E, V, H> {
    pub fn new() -> Self {
        HedgeGraphBuilder {
            hedge_data: HedgeVec::new(),
            nodes: Vec::new(),
            involution: Involution::new(),
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
        let mut hedges = Vec::new();
        for &node in order {
            assert!(!used[node.0], "retained nodes must be unique");
            used[node.0] = true;
            hedges.extend_from_slice(&self.nodes[node.0].hedges);
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
        let mut nodes = self.nodes.into_iter().map(Some).collect::<Vec<_>>();
        self.nodes = order
            .iter()
            .map(|node| {
                let mut node = nodes[node.0].take().expect("unique retained node");
                for hedge in &mut node.hedges {
                    *hedge = mapping[*hedge].expect("retained node hedge");
                }
                node
            })
            .collect();
        (self, mapping)
    }

    pub fn build<N: NodeStorageOps<NodeData = V>>(self) -> HedgeGraph<E, V, H, N> {
        self.into()
    }

    pub fn build_with_map<V2, N: NodeStorageOps<NodeData = V2>>(
        self,
        map: impl FnMut(V) -> V2,
    ) -> HedgeGraph<E, V2, H, N> {
        let len: Hedge = self.involution.len();

        HedgeGraph {
            node_store: N::build_with_mapping(self.nodes, len.0, map),
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
        self.nodes.push(HedgeNodeBuilder {
            data,
            hedges: Vec::new(),
        });
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
        let (sourceh, sinkh) = self.involution.add_pair(data, directed);
        self.hedge_data.push(source.data);
        self.hedge_data.push(sink.data);

        self.nodes[source.node.0].hedges.push(sourceh);
        self.nodes[sink.node.0].hedges.push(sinkh);
    }

    pub fn add_external_edge(
        &mut self,
        source: impl Into<HedgeData<H>>,
        data: E,
        orientation: impl Into<Orientation>,
        underlying: Flow,
    ) {
        let source = source.into();
        let id = self.involution.add_identity(data, orientation, underlying);
        self.nodes[source.node.0].hedges.push(id);
        self.hedge_data.push(source.data);
    }
}

impl<E, V, H> Default for HedgeGraphBuilder<E, V, H> {
    fn default() -> Self {
        Self::new()
    }
}

impl<E, V, H, N: NodeStorageOps<NodeData = V>> From<HedgeGraphBuilder<E, V, H>>
    for HedgeGraph<E, V, H, N>
{
    fn from(builder: HedgeGraphBuilder<E, V, H>) -> Self {
        let len: Hedge = builder.involution.len();

        HedgeGraph {
            node_store: N::build(builder.nodes, len.0),
            edge_store: SmartEdgeVec::new(builder.involution),
            hedge_data: builder.hedge_data,
        }
    }
}
