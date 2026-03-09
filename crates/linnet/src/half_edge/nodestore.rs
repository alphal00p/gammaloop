#[cfg(feature = "nodestore-forest-parent-child")]
use crate::tree::child_pointer::ParentChildStore;
#[cfg(feature = "nodestore-forest-child-vec")]
use crate::tree::child_vec::ChildVecStore;
use crate::tree::{parent_pointer::ParentPointerStore, Forest};

use super::{
    builder::HedgeNodeBuilder,
    involution::{EdgeIndex, Hedge, Involution},
    subgraph::{BaseSubgraph, SubSetLike},
    swap::Swap,
    HedgeGraph, HedgeGraphError, NodeIndex, NodeVec,
};

/// Defines a comprehensive set of operations for a node storage backend within a [`HedgeGraph`].
///
/// This trait extends [`NodeStorage`] by adding methods for modifying the storage (e.g.,
/// extending, extracting, deleting, identifying nodes), building it, iterating over nodes,
/// and mapping node data. It abstracts over the specific data structures used to store
/// node information and their incident half-edges.
pub trait NodeStorageOps: NodeStorage + Swap<Hedge> + Swap<NodeIndex> {
    /// The type of storage returned after mapping operations that might change the node data type.
    type OpStorage<N>: NodeStorageOps<NodeData = N>;
    /// The base type of subgraph used by this node storage (e.g., a `BitVec` to represent a set of hedges).
    type Base: BaseSubgraph;
    // where
    // Self: 'a;
    fn extend(self, other: Self) -> Self;

    fn new_nodevec<'a, V2>(
        &'a self,
        node_map: impl FnMut(NodeIndex, Self::NeighborsIter<'a>, &'a Self::NodeData) -> V2,
    ) -> NodeVec<V2>;

    fn extend_mut(&mut self, other: Self);
    fn extract<S: SubSetLike<Base = Self::Base>, V2>(
        &mut self,
        subgraph: &S,
        split_node: impl FnMut(&Self::NodeData) -> V2,
        owned_node: impl FnMut(Self::NodeData) -> V2,
    ) -> Self::OpStorage<V2>;

    fn extract_nodes(&mut self, nodes: impl IntoIterator<Item = NodeIndex>) -> (Self::Base, Self);

    fn delete<S: SubSetLike<Base = Self::Base>>(&mut self, subgraph: &S);

    // fn add_node(&mut self, node_data: Self::NodeData) -> NodeIndex;

    /// Identifies nodes, essentially turning them into a single node
    /// Invalidates all previous NodeIndex values if using VecNodeStore.
    /// Does not invalidate if using Forest as NodeStore
    fn identify_nodes(&mut self, nodes: &[NodeIndex], node_data_merge: Self::NodeData)
        -> NodeIndex;

    fn forget_identification_history(&mut self) -> NodeVec<Self::NodeData>;

    fn to_forest<U, H>(
        &self,
        map_data: impl Fn(&Self::NodeData) -> U,
    ) -> Forest<U, ParentPointerStore<H>>;

    fn build<I: IntoIterator<Item = HedgeNodeBuilder<Self::NodeData>>>(
        nodes: I,
        n_hedges: usize,
    ) -> Self;

    fn build_with_mapping<I: IntoIterator<Item = HedgeNodeBuilder<ND>>, ND>(
        nodes: I,
        n_hedges: usize,
        map_data: impl FnMut(ND) -> Self::NodeData,
    ) -> Self;

    fn add_dangling_edge(self, source: NodeIndex) -> Result<Self, HedgeGraphError>;

    fn random(sources: &[Self::Neighbors], sinks: &[Self::Neighbors]) -> Self
    where
        Self::NodeData: Default;

    fn drain(self) -> impl Iterator<Item = (NodeIndex, Self::NodeData)>;
    fn iter(&self) -> impl Iterator<Item = (NodeIndex, &Self::NodeData)>;

    fn check_and_set_nodes(&mut self) -> Result<(), HedgeGraphError>;

    fn check_nodes(&self) -> Result<(), HedgeGraphError>;

    fn map_data_ref_graph<'a, E, V2, H>(
        &'a self,
        graph: &'a HedgeGraph<E, Self::NodeData, H, Self>,
        node_map: impl FnMut(
            &'a HedgeGraph<E, Self::NodeData, H, Self>,
            Self::NeighborsIter<'a>,
            &'a Self::NodeData,
        ) -> V2,
    ) -> Self::OpStorage<V2>;

    fn map_data_ref_mut_graph<'a, V2>(
        &'a mut self,
        node_map: impl FnMut(Self::NeighborsIter<'a>, &'a mut Self::NodeData) -> V2,
    ) -> Self::OpStorage<V2>;

    fn map_data_ref_graph_result<'a, E, V2, H, Er>(
        &'a self,
        graph: &'a HedgeGraph<E, Self::NodeData, H, Self>,
        node_map: impl FnMut(
            &'a HedgeGraph<E, Self::NodeData, H, Self>,
            Self::NeighborsIter<'a>,
            &'a Self::NodeData,
        ) -> Result<V2, Er>,
    ) -> Result<Self::OpStorage<V2>, Er>;

    fn map_data_graph<'a, V2>(
        self,
        involution: &'a Involution<EdgeIndex>,
        f: impl FnMut(&'a Involution<EdgeIndex>, NodeIndex, Self::NodeData) -> V2,
    ) -> Self::OpStorage<V2>;

    fn map_data_graph_result<'a, V2, Err>(
        self,
        involution: &'a Involution<EdgeIndex>,
        f: impl FnMut(&'a Involution<EdgeIndex>, NodeIndex, Self::NodeData) -> Result<V2, Err>,
    ) -> Result<Self::OpStorage<V2>, Err>;

    fn iter_node_id(&self) -> impl Iterator<Item = NodeIndex> {
        (0..<Self as Swap<NodeIndex>>::len(self).0).map(NodeIndex)
    }
    fn iter_nodes(
        &self,
    ) -> impl Iterator<Item = (NodeIndex, Self::NeighborsIter<'_>, &Self::NodeData)>;
    fn iter_nodes_mut(
        &mut self,
    ) -> impl Iterator<Item = (NodeIndex, Self::NeighborsIter<'_>, &mut Self::NodeData)>;
    /// Retrieves the [`NodeIndex`] that the given [`Hedge`] is incident to (or originates from).
    fn node_id_ref(&self, hedge: Hedge) -> NodeIndex;
    /// Returns an iterator over the half-edges incident to the specified `node_id`.
    fn get_neighbor_iterator(&self, node_id: NodeIndex) -> Self::NeighborsIter<'_>;
    /// Gets a reference to the data associated with `node_id`.
    fn get_node_data(&self, node_id: NodeIndex) -> &Self::NodeData;
    /// Gets a mutable reference to the data associated with `node_id`.
    fn get_node_data_mut(&mut self, node_id: NodeIndex) -> &mut Self::NodeData;
}

/// The fundamental trait defining the interface for a node storage backend in a [`HedgeGraph`].
///
/// This trait abstracts the way nodes and their incident half-edges are stored. Different
/// implementations can offer various performance trade-offs for accessing and
/// manipulating node data and topology.
pub trait NodeStorage: Sized {
    /// The type of data stored for each node in the graph.
    type NodeData;
    /// A generic way to refer to this storage type when it holds node data of type `N`.
    /// This is often used in `NodeStorageOps` for methods that transform node data.
    type Storage<N>; // TODO: This might need better definition or usage clarification.

    /// The type used to represent the collection of half-edges incident to a single node.
    /// This must implement the [`SubGraph`] trait.
    type Neighbors: SubSetLike;
    /// An iterator that yields the [`Hedge`] identifiers incident to a node.
    /// It must be cloneable and provide an exact size hint.
    type NeighborsIter<'a>: ExactSizeIterator<Item = Hedge> + Clone
    where
        Self: 'a;
}

mod bitvec_find;
mod forest;
mod vec;

pub use vec::BitVecNeighborIter;
pub use vec::NodeStorageVec;

#[cfg(all(
    feature = "nodestore-forest-child-vec",
    feature = "nodestore-forest-parent-child"
))]
compile_error!(
    "Select only one forest nodestore feature: nodestore-forest-child-vec or nodestore-forest-parent-child."
);

#[cfg(feature = "nodestore-forest-child-vec")]
pub type DefaultNodeStore<V> = Forest<V, ChildVecStore<()>>;

#[cfg(all(
    not(feature = "nodestore-forest-child-vec"),
    feature = "nodestore-forest-parent-child"
))]
pub type DefaultNodeStore<V> = Forest<V, ParentChildStore<()>>;

#[cfg(all(
    not(feature = "nodestore-forest-child-vec"),
    not(feature = "nodestore-forest-parent-child")
))]
pub type DefaultNodeStore<V> = NodeStorageVec<V>;

#[cfg(test)]
mod test;
