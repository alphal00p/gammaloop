use crate::{
    half_edge::{
        involution::Hedge,
        nodestore::{BitVecNeighborIter, NodeStorage, NodeStorageOps, NodeStorageVec},
        subgraph::SuBitGraph,
        swap::Swap,
        NodeIndex, NodeVec,
    },
    tree::{parent_pointer::PPNode, Forest, RootData, RootId},
    union_find::{SetIndex, UFNode, UnionFind},
};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Stores the data for a single graph node when using [`UnionFindNodeStore`].
///
/// This struct combines the custom node data (`V`) with a `BitVec` representing
/// the set of half-edges incident to this node. It serves as the data type `U`
/// within the `UnionFind<U>` structure used by `UnionFindNodeStore`.
///
/// # Type Parameters
///
/// - `V`: The type of custom data associated with the graph node.
pub struct HedgeNodeStore<V> {
    /// The custom data for the graph node.
    data: V,
    /// A bitmask representing the set of half-edges incident to this node.
    /// Each bit position corresponds to a [`Hedge`] index.
    node: SuBitGraph,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// An implementation of [`NodeStorage`] and [`NodeStorageOps`] that uses a
/// [`UnionFind`] data structure to manage nodes.
///
/// This storage strategy is particularly useful when frequent node identification
/// (merging of nodes) is required, as `UnionFind` is highly optimized for such
/// operations. Each set in the `UnionFind` structure corresponds to a graph node,
/// and the data associated with each set is a [`HedgeNodeStore<V>`], which holds
/// the node's custom data and its incident half-edges.
///
/// Note: Several operations in the `NodeStorageOps` trait are marked as `todo!()`
/// in this implementation, suggesting it might be specialized or partially implemented.
///
/// # Type Parameters
///
/// - `V`: The type of custom data associated with each graph node.
pub struct UnionFindNodeStore<V> {
    /// The underlying [`UnionFind`] data structure. Each disjoint set in this
    /// structure represents a node in the graph, and the data associated with
    /// each set is a [`HedgeNodeStore<V>`].
    pub nodes: UnionFind<HedgeNodeStore<V>>,
}

impl<V> From<NodeStorageVec<V>> for UnionFindNodeStore<V> {
    fn from(vecnodes: NodeStorageVec<V>) -> Self {
        let node_data: Vec<_> = vecnodes
            .node_data
            .into_iter()
            .zip(vecnodes.nodes)
            .map(|((_, v), (_, h))| HedgeNodeStore { data: v, node: h })
            .collect();

        let hairs: Vec<_> = node_data.iter().map(|a| a.node.clone()).collect();

        UnionFindNodeStore {
            nodes: UnionFind::from_partition(node_data.into_iter().zip(hairs).collect()).unwrap(),
        }
    }
}

impl SetIndex {
    pub fn node_ref(&self) -> &NodeIndex {
        unsafe { std::mem::transmute(self) }
    }
}

impl<V> NodeStorage for UnionFindNodeStore<V> {
    type Storage<N> = UnionFindNodeStore<N>;
    type NodeData = V;
    type Neighbors = SuBitGraph;
    type NeighborsIter<'a>
        = BitVecNeighborIter<'a>
    where
        Self: 'a;
}

impl<V> UnionFindNodeStore<V> {
    pub fn root_hedge(&self, node: NodeIndex) -> Hedge {
        self.nodes[&SetIndex(node.0)].root_pointer
    }
}

impl<V> Swap<Hedge> for UnionFindNodeStore<V> {
    fn is_zero_length(&self) -> bool {
        todo!()
    }
    fn len(&self) -> Hedge {
        todo!()
    }
    fn swap(&mut self, _i: Hedge, _j: Hedge) {
        todo!()
    }
}
impl<V> Swap<NodeIndex> for UnionFindNodeStore<V> {
    fn is_zero_length(&self) -> bool {
        todo!()
    }
    fn len(&self) -> NodeIndex {
        todo!()
    }
    fn swap(&mut self, _i: NodeIndex, _j: NodeIndex) {
        todo!()
    }
}

impl<V> NodeStorageOps for UnionFindNodeStore<V> {
    type OpStorage<A> = Self::Storage<A>;
    type Base = SuBitGraph;

    fn check_nodes(&self) -> Result<(), crate::half_edge::HedgeGraphError> {
        todo!()
    }
    fn extract_nodes(&mut self, _nodes: impl IntoIterator<Item = NodeIndex>) -> (SuBitGraph, Self) {
        todo!()
    }
    fn forget_identification_history(&mut self) -> NodeVec<Self::NodeData> {
        todo!()
    }

    fn delete<S: crate::half_edge::subgraph::SubSetLike<Base = Self::Base>>(
        &mut self,
        _subgraph: &S,
    ) {
        todo!()
    }

    fn new_nodevec<'a, V2>(
        &'a self,
        _node_map: impl FnMut(NodeIndex, Self::NeighborsIter<'a>, &'a Self::NodeData) -> V2,
    ) -> NodeVec<V2> {
        todo!()
    }

    fn extract<S: crate::half_edge::subgraph::SubSetLike<Base = Self::Base>, V2>(
        &mut self,
        _subgraph: &S,
        _spit_node: impl FnMut(&Self::NodeData) -> V2,
        _owned_node: impl FnMut(Self::NodeData) -> V2,
    ) -> Self::OpStorage<V2> {
        todo!()
    }

    fn get_neighbor_iterator(&self, node_id: NodeIndex) -> Self::NeighborsIter<'_> {
        (&self.nodes[SetIndex(node_id.0)].node).into()
    }

    fn identify_nodes(
        &mut self,
        nodes: &[NodeIndex],
        node_data_merge: Self::NodeData,
    ) -> NodeIndex {
        let hedges = nodes
            .iter()
            .map(|n| self.root_hedge(*n))
            .collect::<Vec<_>>();

        let n_init_root = hedges[0];
        for n in hedges.iter().skip(1) {
            self.nodes.union(n_init_root, *n, |a, _| a);
        }
        self.nodes.replace_set_data_of(n_init_root, |mut a| {
            a.data = node_data_merge;
            a
        });
        NodeIndex(self.nodes.find_data_index(n_init_root).0)
    }

    fn to_forest<U, H>(
        &self,
        map_data: impl Fn(&Self::NodeData) -> U,
    ) -> crate::tree::Forest<U, crate::tree::parent_pointer::ParentPointerStore<H>> {
        Forest {
            roots: self
                .nodes
                .set_data
                .iter()
                .map(|a| RootData {
                    root_id: a.root_pointer.into(),
                    data: map_data(&a.data.as_ref().unwrap().data),
                })
                .collect(),
            nodes: self
                .nodes
                .nodes
                .iter()
                .map(|a| {
                    let ad = a.get();
                    let n = match &ad {
                        UFNode::Child(c) => PPNode::dataless_child((*c).into()),
                        UFNode::Root { set_data_idx, .. } => {
                            PPNode::dataless_root(RootId(set_data_idx.0))
                        }
                    };

                    a.set(ad);

                    n
                })
                .collect(),
        }
    }

    fn check_and_set_nodes(&mut self) -> Result<(), crate::half_edge::HedgeGraphError> {
        Ok(())
    }

    fn map_data_ref_graph<'a, E, V2, H>(
        &'a self,
        graph: &'a crate::half_edge::HedgeGraph<E, Self::NodeData, H, Self>,
        mut node_map: impl FnMut(
            &'a crate::half_edge::HedgeGraph<E, Self::NodeData, H, Self>,
            Self::NeighborsIter<'a>,
            &'a Self::NodeData,
        ) -> V2,
    ) -> Self::OpStorage<V2> {
        UnionFindNodeStore {
            nodes: self.nodes.map_set_data_ref(|n| HedgeNodeStore {
                data: node_map(graph, (&n.node).into(), &n.data),
                node: n.node.clone(),
            }),
        }
    }

    fn map_data_ref_mut_graph<'a, V2>(
        &'a mut self,
        mut node_map: impl FnMut(Self::NeighborsIter<'a>, &'a mut Self::NodeData) -> V2,
    ) -> Self::OpStorage<V2> {
        UnionFindNodeStore {
            nodes: self.nodes.map_set_data_ref_mut(|n| HedgeNodeStore {
                data: node_map((&n.node).into(), &mut n.data),
                node: n.node.clone(),
            }),
        }
    }

    fn map_data_ref_graph_result<'a, E, V2, H, Er>(
        &'a self,
        graph: &'a crate::half_edge::HedgeGraph<E, Self::NodeData, H, Self>,
        mut node_map: impl FnMut(
            &'a crate::half_edge::HedgeGraph<E, Self::NodeData, H, Self>,
            Self::NeighborsIter<'a>,
            &'a Self::NodeData,
        ) -> Result<V2, Er>,
    ) -> Result<Self::OpStorage<V2>, Er> {
        Ok(UnionFindNodeStore {
            nodes: self.nodes.map_set_data_ref_result(|n| {
                match node_map(graph, (&n.node).into(), &n.data) {
                    Ok(data) => Ok(HedgeNodeStore {
                        data,
                        node: n.node.clone(),
                    }),
                    Err(err) => Err(err),
                }
            })?,
        })
    }

    fn iter_nodes(
        &self,
    ) -> impl Iterator<Item = (NodeIndex, Self::NeighborsIter<'_>, &Self::NodeData)> {
        self.nodes
            .iter_set_data()
            .map(|(i, d)| (NodeIndex(i.0), (&d.node).into(), &d.data))
    }
    fn iter_nodes_mut(
        &mut self,
    ) -> impl Iterator<Item = (NodeIndex, Self::NeighborsIter<'_>, &mut Self::NodeData)> {
        self.nodes
            .iter_set_data_mut()
            .map(|(i, d)| (NodeIndex(i.0), (&d.node).into(), &mut d.data))
    }

    fn node_id_ref(&self, hedge: crate::half_edge::involution::Hedge) -> NodeIndex {
        NodeIndex(self.nodes.find_data_index(hedge).0)
    }

    fn get_node_data_mut(&mut self, node_id: NodeIndex) -> &mut Self::NodeData {
        &mut self.nodes[SetIndex(node_id.0)].data
    }

    fn iter_node_id(&self) -> impl Iterator<Item = NodeIndex> {
        self.nodes.iter_set_data().map(|(i, _)| NodeIndex(i.0))
    }

    fn get_node_data(&self, node_id: NodeIndex) -> &Self::NodeData {
        &self.nodes[SetIndex(node_id.0)].data
    }

    fn map_data_graph<'a, V2>(
        self,
        involution: &'a crate::half_edge::involution::Involution<
            crate::half_edge::involution::EdgeIndex,
        >,
        mut f: impl FnMut(
            &'a crate::half_edge::involution::Involution<crate::half_edge::involution::EdgeIndex>,
            NodeIndex,
            Self::NodeData,
        ) -> V2,
    ) -> Self::OpStorage<V2> {
        UnionFindNodeStore {
            nodes: self.nodes.map_set_data(|i, n| HedgeNodeStore {
                data: f(involution, NodeIndex(i.0), n.data),
                node: n.node.clone(),
            }),
        }
    }
    fn map_data_graph_result<'a, V2, Err>(
        self,
        _involution: &'a crate::half_edge::involution::Involution<
            crate::half_edge::involution::EdgeIndex,
        >,
        _f: impl FnMut(
            &'a crate::half_edge::involution::Involution<crate::half_edge::involution::EdgeIndex>,
            NodeIndex,
            Self::NodeData,
        ) -> Result<V2, Err>,
    ) -> Result<Self::OpStorage<V2>, Err> {
        todo!()
    }

    fn extend(mut self, other: Self) -> Self {
        self.nodes.extend(other.nodes);
        self
    }

    fn extend_mut(&mut self, other: Self) {
        self.nodes.extend(other.nodes);
    }

    fn iter(&self) -> impl Iterator<Item = (NodeIndex, &Self::NodeData)> {
        self.nodes
            .iter_set_data()
            .map(|(i, d)| (NodeIndex(i.0), &d.data))
    }

    fn build<
        I: IntoIterator<Item = crate::half_edge::builder::HedgeNodeBuilder<Self::NodeData>>,
    >(
        nodes: I,
        n_hedges: usize,
    ) -> Self {
        NodeStorageVec::build(nodes, n_hedges).into()
    }

    fn build_with_mapping<
        I: IntoIterator<Item = crate::half_edge::builder::HedgeNodeBuilder<ND>>,
        ND,
    >(
        nodes: I,
        n_hedges: usize,
        map_data: impl FnMut(ND) -> Self::NodeData,
    ) -> Self {
        NodeStorageVec::build_with_mapping(nodes, n_hedges, map_data).into()
    }

    fn drain(self) -> impl Iterator<Item = (NodeIndex, Self::NodeData)> {
        self.nodes
            .drain_set_data()
            .map(|(s, d)| (NodeIndex(s.0), d.data))
    }

    fn random(sources: &[Self::Neighbors], sinks: &[Self::Neighbors]) -> Self
    where
        Self::NodeData: Default,
    {
        NodeStorageVec::random(sources, sinks).into()
    }

    fn add_dangling_edge(
        mut self,
        source: NodeIndex,
    ) -> Result<Self, crate::half_edge::HedgeGraphError> {
        let setid = SetIndex(source.0);
        let _ = self.nodes.add_child(setid);
        self.nodes.iter_set_data_mut().for_each(|(s, n)| {
            if s == setid {
                n.node.push(true);
            } else {
                n.node.push(false);
            }
        });
        Ok(self)
    }
}
