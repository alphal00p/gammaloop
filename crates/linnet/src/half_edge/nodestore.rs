//! Node storage: a [`Topology`] partitions the half-edges into nodes and a
//! [`NodeStore`] pairs it with one payload per node.

use std::collections::BTreeSet;

use crate::tree::{
    parent_pointer::{PPNode, ParentPointerStore},
    Forest, RootData, RootId, TreeNodeId,
};

use super::{
    involution::Hedge,
    remap::{HedgeRemap, NodeRemap, Renumbering, SplitRenumbering},
    subgraph::{BaseSubgraph, Inclusion, SuBitGraph, SubSetLike},
    HedgeGraphError, NodeIndex, NodeVec,
};

/// The node partition of a graph's half-edges, without node payloads.
///
/// Nodes are numbered `0..node_bound()` and each one is *live* or *absorbed*.
/// [`Self::identify`] merges nodes without renumbering anything: the merged
/// nodes stay in place, absorbed and without half-edges, until
/// [`Self::compact`] drops them. Implementations agree on every observable
/// result, numbering included, so one store can replace another:
///
/// - Crowns ascend in every state.
/// - The survivor of an identification is the requested node with the largest
///   crown, the first such in request order on ties; repeated and absorbed nodes
///   resolve to their [`Self::representative`]. The survivor keeps its anchor or
///   adopts the first absorbed one in request order.
/// - [`Self::push_node`], [`Self::attach`] and [`Self::extend`] append behind the
///   existing half-edges and nodes, which keep their indices. Half-edges and
///   node slots stay fewer than `2^31`, so that remaps can flag merges in the
///   top bit; growing beyond panics.
/// - [`Self::compact`] fills the absorbed slots from the tail with the shared
///   kernel [`super::swap::partition_with`].
/// - [`Self::partition_hedges`] moves half-edges with the same kernel, so the
///   edge and hedge stores stay aligned, and orders nodes as
///   `[Kept | Split | Moved]` with one kernel pass per class.
/// - A node's anchor follows its half-edge and falls back to the smallest
///   remaining crown half-edge when its own half-edge leaves the node.
///
/// Renumbering operations ([`Self::compact`], [`Self::partition_hedges`],
/// [`Self::permute`], [`Self::extend`]) fail with
/// [`HedgeGraphError::CheckpointOpen`] while a [`Checkpoint`] is open;
/// [`Self::partition_hedges`] and [`Self::permute`] also fail with
/// [`HedgeGraphError::PendingIdentifications`] until absorbed nodes are compacted.
pub trait Topology: Clone + Sized {
    /// Iterates over a crown in ascending half-edge order.
    type Crown<'a>: ExactSizeIterator<Item = Hedge> + Clone
    where
        Self: 'a;

    /// Number of half-edges.
    fn n_hedges(&self) -> usize;
    /// Number of node slots, live or absorbed: the bound of every [`NodeIndex`]
    /// and the size of dense side tables.
    fn node_bound(&self) -> usize;
    /// Number of live nodes.
    fn n_live(&self) -> usize;
    /// The live node of `hedge`.
    fn node_of(&self, hedge: Hedge) -> NodeIndex;
    /// The half-edges of `node`, ascending; empty for absorbed nodes.
    fn crown(&self, node: NodeIndex) -> Self::Crown<'_>;
    /// Whether `node` has not been absorbed by an identification.
    fn is_live(&self, node: NodeIndex) -> bool;
    /// The live node that `node` was, possibly transitively, identified into;
    /// `node` itself while it is live.
    fn representative(&self, node: NodeIndex) -> NodeIndex;
    /// The anchor half-edge of `node`, if it has one.
    fn anchor(&self, node: NodeIndex) -> Option<Hedge>;
    /// Sets or clears the anchor of the live `node`.
    ///
    /// # Panics
    ///
    /// If `anchor` is not in the crown of `node`.
    fn set_anchor(&mut self, node: NodeIndex, anchor: Option<Hedge>);

    /// Builds `node_bound` live nodes without anchors from each half-edge's node.
    ///
    /// # Panics
    ///
    /// If a label is not below `node_bound`.
    fn from_labels(labels: impl ExactSizeIterator<Item = NodeIndex>, node_bound: usize) -> Self;
    /// Appends a live node without half-edges.
    fn push_node(&mut self) -> NodeIndex;
    /// Appends an unpaired half-edge to the live `node`.
    ///
    /// Fails with [`HedgeGraphError::NoNode`] if `node` is out of range and with
    /// [`HedgeGraphError::InvalidNode`] if it was absorbed, both before changing
    /// anything.
    fn attach(&mut self, node: NodeIndex) -> Result<Hedge, HedgeGraphError>;
    /// Appends `other` behind `self`, shifting its half-edges and nodes; pending
    /// identifications on either side are kept.
    fn extend(&mut self, other: Self) -> Result<(), HedgeGraphError>;
    /// Merges `nodes` into one survivor and returns it; renumbers nothing.
    ///
    /// # Panics
    ///
    /// If `nodes` is empty.
    fn identify(&mut self, nodes: &[NodeIndex]) -> NodeIndex;
    /// Drops the absorbed nodes. Absorbed nodes map to their representative's
    /// new index and are marked merged; live nodes are the origins.
    fn compact(&mut self) -> Result<NodeRemap, HedgeGraphError>;
    /// Moves the half-edges selected by `moved` behind the others and splits the
    /// nodes by `fate`, which sees each node with its numbers of staying and
    /// moved half-edges. `self` keeps the staying half-edges and the
    /// [`NodeFate::Kept`] and [`NodeFate::Split`] nodes; with `build_moved`, the
    /// moved half-edges with the split and moved nodes become a second topology.
    fn partition_hedges(
        &mut self,
        moved: impl Fn(Hedge) -> bool,
        fate: impl Fn(NodeIndex, usize, usize) -> NodeFate,
        build_moved: bool,
    ) -> Result<Partition<Self>, HedgeGraphError>;
    /// Renumbers half-edges and nodes by the bijections `hedges` and `nodes`.
    fn permute(&mut self, hedges: &HedgeRemap, nodes: &NodeRemap) -> Result<(), HedgeGraphError>;

    /// Checks the contract through the interface: crowns ascend, agree with
    /// [`Self::node_of`] and cover every half-edge once; absorbed nodes have no
    /// crown or anchor and lead to live nodes; anchors lie in their crowns; and
    /// [`Self::n_live`] counts the live nodes.
    fn check(&self) -> Result<(), HedgeGraphError> {
        let partition = HedgeGraphError::NodesDoNotPartition;
        let (n_hedges, mut covered, mut n_live) = (self.n_hedges(), 0, 0);
        for node in (0..self.node_bound()).map(NodeIndex) {
            let crown = self.crown(node);
            if !self.is_live(node) {
                if crown.len() != 0 || self.anchor(node).is_some() {
                    return Err(partition(format!("absorbed node {node} has a crown")));
                }
                if !self.is_live(self.representative(node)) {
                    return Err(partition(format!("absorbed {node} reaches no live node")));
                }
                continue;
            }
            n_live += 1;
            covered += crown.len();
            let mut previous = None;
            for hedge in crown {
                if hedge.0 >= n_hedges {
                    return Err(HedgeGraphError::InvalidHedge(hedge));
                }
                if self.node_of(hedge) != node {
                    return Err(HedgeGraphError::NodeDoesNotContainHedge(node, hedge));
                }
                if previous.is_some_and(|previous| previous >= hedge) {
                    return Err(partition(format!("crown of {node} does not ascend")));
                }
                previous = Some(hedge);
            }
            if self
                .anchor(node)
                .is_some_and(|anchor| anchor.0 >= n_hedges || self.node_of(anchor) != node)
            {
                return Err(partition(format!("anchor of {node} is outside its crown")));
            }
        }
        if covered != n_hedges {
            return Err(partition("crowns do not cover the half-edges".to_owned()));
        }
        if n_live != self.n_live() {
            return Err(partition("live node count is stale".to_owned()));
        }
        Ok(())
    }
}

/// Undo support for the edits that renumber nothing.
///
/// Between [`Self::checkpoint`] and the matching [`Self::rollback`] or
/// [`Self::commit`], [`Topology::push_node`], [`Topology::attach`],
/// [`Topology::identify`] and [`Topology::set_anchor`] are recorded, while the
/// renumbering operations fail with [`HedgeGraphError::CheckpointOpen`].
/// Checkpoints nest. They cover the topology only: payloads that a
/// [`NodeStore`] overwrites are not restored.
pub trait Checkpoint: Topology {
    /// Opens a checkpoint.
    fn checkpoint(&mut self);
    /// Undoes the edits since the innermost open checkpoint and closes it.
    ///
    /// # Panics
    ///
    /// If no checkpoint is open.
    fn rollback(&mut self);
    /// Keeps the edits since the innermost open checkpoint and closes it.
    ///
    /// # Panics
    ///
    /// If no checkpoint is open.
    fn commit(&mut self);
    /// Whether a checkpoint is open.
    fn has_checkpoint(&self) -> bool;
}

/// Where [`Topology::partition_hedges`] puts a node.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NodeFate {
    /// Stays, with all its half-edges.
    Kept,
    /// Stays and also appears on the moved side, each copy with its own
    /// half-edges.
    Split,
    /// Leaves with all its half-edges.
    Moved,
}

impl NodeFate {
    /// The fate in a removal of half-edges: untouched nodes stay, nodes that
    /// keep some half-edges split, and nodes that lose all of them leave;
    /// isolated nodes stay.
    pub(crate) fn of_removal(_: NodeIndex, staying: usize, moved: usize) -> Self {
        match (staying, moved) {
            (_, 0) => NodeFate::Kept,
            (0, _) => NodeFate::Moved,
            _ => NodeFate::Split,
        }
    }
}

/// The result of [`Topology::partition_hedges`].
#[derive(Clone, Debug)]
pub struct Partition<T> {
    /// The moved half-edges with the split and moved nodes, if requested.
    pub moved: Option<T>,
    /// Each old half-edge's position in `[staying | moved]`.
    pub hedges: HedgeRemap,
    /// The number of staying half-edges.
    pub n_staying: usize,
    /// Each old node's position in `[Kept | Split | Moved]`.
    pub nodes: NodeRemap,
    /// The first split node.
    pub split: usize,
    /// The first moved node.
    pub inside: usize,
}

impl<T> Partition<T> {
    /// Where the staying half-edges and the kept and split nodes went.
    fn into_kept(self) -> Renumbering {
        Renumbering {
            hedges: self.hedges.window(0..self.n_staying),
            nodes: self.nodes.window(0..self.inside),
        }
    }

    /// The moved side, and where the half-edges and nodes of both sides went.
    fn into_split(self) -> (Option<T>, SplitRenumbering) {
        let (n_hedges, n_nodes) = (self.hedges.new_bound(), self.nodes.new_bound());
        let kept = Renumbering {
            hedges: self.hedges.clone().window(0..self.n_staying),
            nodes: self.nodes.clone().window(0..self.inside),
        };
        let moved = Renumbering {
            hedges: self.hedges.window(self.n_staying..n_hedges),
            nodes: self.nodes.window(self.split..n_nodes),
        };
        (self.moved, SplitRenumbering { kept, moved })
    }
}

/// A [`Topology`] with one payload per node slot, absorbed nodes included.
///
/// Everything that does not depend on how the topology is stored lives here:
/// payload access and mapping, iteration, building, consistency checks and
/// moving payloads along the topology's renumberings.
///
/// A store lends its topology out only by shared reference, so it never holds
/// an open [`Checkpoint`]: its renumberings fail only on pending
/// identifications, and compacting and appending cannot fail.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "serde",
    serde(
        try_from = "serialization::Parts<T, V>",
        bound(deserialize = "T: Topology + serde::Deserialize<'de>, V: serde::Deserialize<'de>")
    )
)]
pub struct NodeStore<T, V> {
    topology: T,
    data: NodeVec<V>,
}

const NO_CHECKPOINT: &str = "node stores never hold an open checkpoint";

impl<T: Topology, V> NodeStore<T, V> {
    /// Pairs `topology` with one payload per node slot.
    ///
    /// # Panics
    ///
    /// If the payload count differs from the node bound.
    pub(crate) fn new(topology: T, data: NodeVec<V>) -> Self {
        assert_eq!(
            data.0.len(),
            topology.node_bound(),
            "one payload per node slot"
        );
        NodeStore { topology, data }
    }

    /// The node partition.
    pub fn topology(&self) -> &T {
        &self.topology
    }

    /// Builds nodes from source and sink neighbor sets with default payloads.
    pub(crate) fn random(sources: &[SuBitGraph], sinks: &[SuBitGraph]) -> Self
    where
        V: Default,
    {
        let n_hedges = sources.first().map_or(0, |source| source.size());
        let mut label = vec![NodeIndex(0); n_hedges];
        for (node, crown) in sources.iter().chain(sinks).enumerate() {
            for hedge in crown.included_iter() {
                label[hedge.0] = NodeIndex(node);
            }
        }
        let n_nodes = sources.len() + sinks.len();
        let topology = T::from_labels(label.into_iter(), n_nodes);
        let data = (0..n_nodes).map(|_| V::default()).collect();
        NodeStore { topology, data }
    }

    /// Number of node slots, live or absorbed: the size of dense side tables.
    pub fn node_bound(&self) -> usize {
        self.topology.node_bound()
    }

    /// Number of live nodes.
    pub fn n_live(&self) -> usize {
        self.topology.n_live()
    }

    /// Number of half-edges.
    pub fn n_hedges(&self) -> usize {
        self.topology.n_hedges()
    }

    /// The node of `hedge`.
    #[inline]
    pub fn node_id(&self, hedge: Hedge) -> NodeIndex {
        self.topology.node_of(hedge)
    }

    /// The half-edges of `node`, ascending.
    #[inline]
    pub fn crown(&self, node: NodeIndex) -> T::Crown<'_> {
        self.topology.crown(node)
    }

    /// Whether `node` has not been absorbed by an identification.
    pub fn is_live(&self, node: NodeIndex) -> bool {
        self.topology.is_live(node)
    }

    /// The anchor half-edge of `node`, if it has one.
    pub fn anchor(&self, node: NodeIndex) -> Option<Hedge> {
        self.topology.anchor(node)
    }

    /// Sets or clears the anchor of `node`; see [`Topology::set_anchor`].
    pub fn set_anchor(&mut self, node: NodeIndex, anchor: Option<Hedge>) {
        self.topology.set_anchor(node, anchor);
    }

    /// The payload of `node`.
    pub fn get(&self, node: NodeIndex) -> &V {
        &self.data[node]
    }

    /// The payload of `node`, mutably.
    pub fn get_mut(&mut self, node: NodeIndex) -> &mut V {
        &mut self.data[node]
    }

    /// Every node slot's payload, absorbed nodes included.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (NodeIndex, &V)> {
        self.data.iter()
    }

    /// The live nodes with their crowns and payloads.
    pub(crate) fn iter_nodes(&self) -> impl Iterator<Item = (NodeIndex, T::Crown<'_>, &V)> {
        self.data
            .iter()
            .filter(|(node, _)| self.topology.is_live(*node))
            .map(|(node, data)| (node, self.topology.crown(node), data))
    }

    /// The live nodes with their crowns and mutable payloads.
    pub(crate) fn iter_nodes_mut(
        &mut self,
    ) -> impl Iterator<Item = (NodeIndex, T::Crown<'_>, &mut V)> {
        let topology = &self.topology;
        self.data
            .iter_mut()
            .filter(|(node, _)| topology.is_live(*node))
            .map(move |(node, data)| (node, topology.crown(node), data))
    }

    /// The live node ids.
    pub(crate) fn iter_node_ids(&self) -> impl Iterator<Item = NodeIndex> + '_ {
        (0..self.node_bound())
            .map(NodeIndex)
            .filter(|&node| self.topology.is_live(node))
    }

    /// A dense table over every node slot.
    pub(crate) fn new_nodevec<'a, U>(
        &'a self,
        mut f: impl FnMut(NodeIndex, T::Crown<'a>, &'a V) -> U,
    ) -> NodeVec<U> {
        self.data
            .iter()
            .map(|(node, data)| f(node, self.topology.crown(node), data))
            .collect()
    }

    /// Maps borrowed payloads, keeping the topology.
    pub(crate) fn map_ref<'a, U>(
        &'a self,
        f: impl FnMut(NodeIndex, T::Crown<'a>, &'a V) -> U,
    ) -> NodeStore<T, U> {
        NodeStore {
            data: self.new_nodevec(f),
            topology: self.topology.clone(),
        }
    }

    /// Fallible counterpart of [`Self::map_ref`].
    pub(crate) fn try_map_ref<'a, U, E>(
        &'a self,
        mut f: impl FnMut(NodeIndex, T::Crown<'a>, &'a V) -> Result<U, E>,
    ) -> Result<NodeStore<T, U>, E> {
        let data = (self.data.iter())
            .map(|(node, data)| f(node, self.topology.crown(node), data))
            .collect::<Result<_, E>>()?;
        Ok(NodeStore {
            topology: self.topology.clone(),
            data,
        })
    }

    /// Maps payloads through mutable references, keeping the topology.
    pub(crate) fn map_mut<'a, U>(
        &'a mut self,
        mut f: impl FnMut(T::Crown<'a>, &'a mut V) -> U,
    ) -> NodeStore<T, U> {
        let topology = &self.topology;
        let data = (self.data.iter_mut())
            .map(|(node, data)| f(topology.crown(node), data))
            .collect();
        NodeStore {
            topology: topology.clone(),
            data,
        }
    }

    /// Consumes the store and maps owned payloads.
    pub(crate) fn map<U>(self, mut f: impl FnMut(NodeIndex, V) -> U) -> NodeStore<T, U> {
        NodeStore {
            data: (self.data.into_iter())
                .map(|(node, data)| f(node, data))
                .collect(),
            topology: self.topology,
        }
    }

    /// Fallible counterpart of [`Self::map`].
    pub(crate) fn try_map<U, E>(
        self,
        mut f: impl FnMut(NodeIndex, V) -> Result<U, E>,
    ) -> Result<NodeStore<T, U>, E> {
        Ok(NodeStore {
            data: (self.data.into_iter())
                .map(|(node, data)| f(node, data))
                .collect::<Result<_, E>>()?,
            topology: self.topology,
        })
    }

    /// Appends a node without half-edges and returns its index.
    pub(crate) fn add_node(&mut self, data: V) -> NodeIndex {
        let node = self.topology.push_node();
        self.data.push(data);
        node
    }

    /// Appends an unpaired half-edge to `node`; see [`Topology::attach`].
    pub(crate) fn attach(&mut self, node: NodeIndex) -> Result<Hedge, HedgeGraphError> {
        self.topology.attach(node)
    }

    /// Appends `other` behind `self`; see [`Topology::extend`].
    pub(crate) fn extend(&mut self, other: Self) {
        (self.topology.extend(other.topology)).expect(NO_CHECKPOINT);
        self.data.extend(other.data);
    }

    /// Merges `nodes` into one survivor carrying `data`; see [`Topology::identify`].
    pub(crate) fn identify(&mut self, nodes: &[NodeIndex], data: V) -> NodeIndex {
        let survivor = self.topology.identify(nodes);
        self.data[survivor] = data;
        survivor
    }

    /// Drops the absorbed nodes with their payloads; see [`Topology::compact`].
    pub(crate) fn compact(&mut self) -> NodeRemap {
        let remap = self.topology.compact().expect(NO_CHECKPOINT);
        remap.fill_holes(&mut self.data.0);
        remap
    }

    /// Removes the half-edges of `subgraph`, and the nodes that lose all of
    /// their half-edges, with their payloads.
    pub(crate) fn delete(
        &mut self,
        subgraph: &impl SubSetLike,
    ) -> Result<Renumbering, HedgeGraphError> {
        let mut partition = self.topology.partition_hedges(
            |hedge| subgraph.includes(&hedge),
            NodeFate::of_removal,
            false,
        )?;
        partition.nodes.permute_in_place(&mut self.data.0);
        self.data.0.truncate(partition.inside);
        Ok(partition.into_kept())
    }

    /// Moves the half-edges of `subgraph` into a new store. Nodes crossed by the
    /// selection boundary stay and leave a copy mapped by `split_node`; nodes
    /// wholly inside leave with their payload mapped by `owned_node`.
    pub(crate) fn extract<U>(
        &mut self,
        subgraph: &impl SubSetLike,
        split_node: impl FnMut(&V) -> U,
        owned_node: impl FnMut(V) -> U,
    ) -> Result<(NodeStore<T, U>, SplitRenumbering), HedgeGraphError> {
        let partition = self.topology.partition_hedges(
            |hedge| subgraph.includes(&hedge),
            NodeFate::of_removal,
            true,
        )?;
        Ok(self.split_off(partition, split_node, owned_node))
    }

    /// Moves the live `nodes`, isolated ones included, into a new store with
    /// their half-edges, which are returned as a subset of the old half-edges.
    pub(crate) fn extract_nodes(
        &mut self,
        nodes: impl IntoIterator<Item = NodeIndex>,
    ) -> Result<(SuBitGraph, Self, SplitRenumbering), HedgeGraphError> {
        let nodes: BTreeSet<NodeIndex> = nodes.into_iter().collect();
        let selected = SuBitGraph::from_hedge_iter(
            nodes.iter().flat_map(|&node| self.topology.crown(node)),
            self.n_hedges(),
        );
        let partition = self.topology.partition_hedges(
            |hedge| selected.includes(&hedge),
            |node, _, _| {
                if nodes.contains(&node) {
                    NodeFate::Moved
                } else {
                    NodeFate::Kept
                }
            },
            true,
        )?;
        let (store, renumbering) =
            self.split_off(partition, |_| unreachable!("no node splits"), |data| data);
        Ok((selected, store, renumbering))
    }

    /// Moves the payloads along a partition built with its moved side.
    fn split_off<U>(
        &mut self,
        mut partition: Partition<T>,
        split_node: impl FnMut(&V) -> U,
        owned_node: impl FnMut(V) -> U,
    ) -> (NodeStore<T, U>, SplitRenumbering) {
        partition.nodes.permute_in_place(&mut self.data.0);
        let mut data: Vec<U> = (self.data.0[partition.split..partition.inside].iter())
            .map(split_node)
            .collect();
        data.extend(self.data.0.drain(partition.inside..).map(owned_node));
        let (topology, renumbering) = partition.into_split();
        let topology = topology.expect("partitions for extraction build their moved side");
        (NodeStore::new(topology, NodeVec(data)), renumbering)
    }

    /// Renumbers half-edges and nodes; see [`Topology::permute`].
    pub(crate) fn permute(
        &mut self,
        hedges: &HedgeRemap,
        nodes: &NodeRemap,
    ) -> Result<(), HedgeGraphError> {
        self.topology.permute(hedges, nodes)?;
        nodes.apply(&mut self.data.0);
        Ok(())
    }

    /// Checks the topology contract and that every node slot has a payload.
    pub fn check_nodes(&self) -> Result<(), HedgeGraphError> {
        if self.data.0.len() != self.node_bound() {
            return Err(HedgeGraphError::DataLengthMismatch);
        }
        self.topology.check()
    }

    /// Copies the topology into a parent-pointer forest, one tree per node slot
    /// rooted at its smallest half-edge with the others as children, while
    /// mapping borrowed payloads; absorbed and isolated nodes become empty trees.
    pub(crate) fn to_forest<U, H>(
        &self,
        map_data: impl Fn(&V) -> U,
    ) -> Forest<U, ParentPointerStore<H>> {
        let mut nodes: Vec<Option<PPNode<H>>> = std::iter::repeat_with(|| None)
            .take(self.n_hedges())
            .collect();
        let mut roots = Vec::with_capacity(self.node_bound());
        for (node, data) in self.data.iter() {
            let mut crown = self.topology.crown(node);
            let root_id = match crown.next() {
                Some(head) => {
                    nodes[head.0] = Some(PPNode::dataless_root(RootId(node.0)));
                    for hedge in crown {
                        nodes[hedge.0] = Some(PPNode::dataless_child(TreeNodeId(head.0)));
                    }
                    TreeNodeId(head.0)
                }
                None => TreeNodeId::EMPTY,
            };
            roots.push(RootData {
                root_id,
                data: map_data(data),
            });
        }
        Forest {
            nodes: nodes
                .into_iter()
                .map(|node| node.expect("crowns partition the half-edges"))
                .collect(),
            roots,
        }
    }
}

#[cfg(feature = "serde")]
mod serialization {
    use super::{NodeStore, NodeVec, Topology};

    #[derive(serde::Deserialize)]
    #[serde(rename = "NodeStore")]
    pub(super) struct Parts<T, V> {
        topology: T,
        data: NodeVec<V>,
    }

    impl<T: Topology, V> TryFrom<Parts<T, V>> for NodeStore<T, V> {
        type Error = String;

        fn try_from(parts: Parts<T, V>) -> Result<Self, String> {
            if parts.data.0.len() != parts.topology.node_bound() {
                return Err("one payload per node slot expected".to_owned());
            }
            Ok(NodeStore {
                topology: parts.topology,
                data: parts.data,
            })
        }
    }
}

#[cfg(feature = "bincode")]
mod binary {
    use bincode::{
        de::{BorrowDecoder, Decoder},
        enc::Encoder,
        error::{DecodeError, EncodeError},
        BorrowDecode, Decode, Encode,
    };

    use super::{NodeStore, NodeVec, Topology};

    impl<T: Encode, V: Encode> Encode for NodeStore<T, V> {
        fn encode<E: Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
            self.topology.encode(encoder)?;
            self.data.encode(encoder)
        }
    }

    fn checked<T: Topology, V>(
        topology: T,
        data: NodeVec<V>,
    ) -> Result<NodeStore<T, V>, DecodeError> {
        if data.0.len() != topology.node_bound() {
            return Err(DecodeError::Other("one payload per node slot expected"));
        }
        Ok(NodeStore { topology, data })
    }

    impl<Context, T: Topology + Decode<Context>, V: Decode<Context>> Decode<Context>
        for NodeStore<T, V>
    {
        fn decode<D: Decoder<Context = Context>>(decoder: &mut D) -> Result<Self, DecodeError> {
            checked(Decode::decode(decoder)?, Decode::decode(decoder)?)
        }
    }

    impl<'de, Context, T, V> BorrowDecode<'de, Context> for NodeStore<T, V>
    where
        T: Topology + BorrowDecode<'de, Context>,
        V: BorrowDecode<'de, Context>,
    {
        fn borrow_decode<D: BorrowDecoder<'de, Context = Context>>(
            decoder: &mut D,
        ) -> Result<Self, DecodeError> {
            checked(
                BorrowDecode::borrow_decode(decoder)?,
                BorrowDecode::borrow_decode(decoder)?,
            )
        }
    }
}

mod node_forest;

#[cfg(feature = "rkyv")]
pub use node_forest::{ArchivedNodeForest, ArchivedNodeStore};
pub use node_forest::{NodeForest, RingNeighborIter};

#[cfg(test)]
pub(crate) mod reference;
#[cfg(test)]
mod test;
