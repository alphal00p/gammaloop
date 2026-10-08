use std::{collections::HashSet, fmt};

use super::{NodeStorage, NodeStorageOps, NodeStorageVec};
use crate::{
    half_edge::{
        involution::Hedge,
        subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
        swap::{partition_with, Swap},
        NodeIndex, NodeVec,
    },
    tree::{
        parent_pointer::{PPNode, ParentPointerStore},
        Forest, ForestNodeStore, ForestNodeStorePreorder, RootData, RootId, TreeNodeId,
    },
};

// pub struct ForestNodeStore<Tree, N> {
//     forest: Forest<N, Tree>,
// }

#[derive(Clone)]
/// An iterator over the half-edges incident to a graph node when the node
/// storage is implemented using a [`Forest`].
///
/// In this context, a graph node corresponds to a tree in the `Forest`, and the
/// half-edges incident to that graph node are represented as nodes within that tree.
/// This iterator typically traverses the nodes of a specific tree in pre-order,
/// yielding [`Hedge`] identifiers.
///
/// # Type Parameters
///
/// - `'a`: The lifetime of the borrow from the underlying forest node store.
/// - `P`: The type of the [`ForestNodeStore`] used by the `Forest`, which must
///   also implement [`ForestNodeStorePreorder`] to allow pre-order traversal
///   of tree nodes (representing half-edges).
pub enum ForestNeighborIter<'a, P: ForestNodeStorePreorder + 'a> {
    Empty,
    Iter(P::Iterator<'a>),
}

impl<'a, P: ForestNodeStorePreorder + 'a> ForestNeighborIter<'a, P> {
    fn from_root(nodes: &'a P, root_id: TreeNodeId) -> Self {
        if root_id.is_empty() {
            Self::Empty
        } else {
            Self::Iter(nodes.iter_preorder(root_id))
        }
    }
}

impl<P: ForestNodeStorePreorder> fmt::Debug for ForestNeighborIter<'_, P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("ForestNeighborIter::Empty"),
            Self::Iter(_) => f.write_str("ForestNeighborIter::Iter(_)"),
        }
    }
}

impl<'a, P: ForestNodeStorePreorder + 'a> Iterator for ForestNeighborIter<'a, P> {
    type Item = Hedge;
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Empty => None,
            Self::Iter(iter) => iter.next().map(|a| a.into()),
        }
    }
}

impl<'a, P: ForestNodeStorePreorder + 'a> ExactSizeIterator for ForestNeighborIter<'a, P> {
    fn len(&self) -> usize {
        match self {
            Self::Empty => 0,
            Self::Iter(iter) => iter.clone().count(),
        }
    }
}

impl<V, P: ForestNodeStore + ForestNodeStorePreorder + Clone> NodeStorage for Forest<V, P> {
    type Storage<N> = Forest<N, P>;
    type NodeData = V;
    type Neighbors = SuBitGraph;
    type NeighborsIter<'a>
        = ForestNeighborIter<'a, P>
    where
        Self: 'a;
}

impl From<RootId> for NodeIndex {
    fn from(value: RootId) -> Self {
        NodeIndex(value.0)
    }
}

impl<V, P: ForestNodeStore> Swap<Hedge> for Forest<V, P> {
    fn len(&self) -> Hedge {
        Hedge(self.nodes.n_nodes())
    }

    fn is_zero_length(&self) -> bool {
        self.nodes.n_nodes() == 0
    }

    fn swap(&mut self, a: Hedge, b: Hedge) {
        assert!(!self.has_identification_history(), "{HISTORY_RETAINED}");
        if a != b {
            let a = a.into();
            let b = b.into();
            self.nodes.swap(a, b);
            let rid = self.nodes.root_node(a);
            let r = self.root(rid);
            self.roots[r.0].root_id = rid;
            let rid = self.nodes.root_node(b);
            let r = self.root(rid);
            self.roots[r.0].root_id = rid;
        }

        #[cfg(test)]
        self.validate_structure().unwrap()
    }
}
impl<V, P: ForestNodeStore> Swap<NodeIndex> for Forest<V, P> {
    fn len(&self) -> NodeIndex {
        NodeIndex(self.roots.len())
    }

    fn is_zero_length(&self) -> bool {
        self.roots.is_empty()
    }
    fn swap(&mut self, i: NodeIndex, j: NodeIndex) {
        assert!(!self.has_identification_history(), "{HISTORY_RETAINED}");
        self.swap_roots(RootId(i.0), RootId(j.0));
    }
}

/// Identification history is retained only until it is forgotten: absorbed
/// roots point at tree nodes by half-edge position and at their survivors by
/// node index, so moving half-edges or renumbering nodes would corrupt it.
const HISTORY_RETAINED: &str =
    "moving half-edges or renumbering nodes of a Forest node store requires \
     forget_identification_history() after identify_nodes()";

impl<V, P: ForestNodeStore + ForestNodeStorePreorder + Clone> Forest<V, P> {
    /// Moves the half-edges of `subgraph` to the end with the same swaps as the
    /// edge store, then orders the nodes as `[untouched | split | inside]`.
    ///
    /// Returns the first moved half-edge and the starts of the split and inside
    /// node ranges.
    fn partition_for_removal<S: SubSetLike>(
        &mut self,
        subgraph: &S,
    ) -> (Hedge, NodeIndex, NodeIndex) {
        assert!(!self.has_identification_history(), "{HISTORY_RETAINED}");
        let left = <Self as Swap<Hedge>>::partition(self, |h| !subgraph.includes(h));
        let split = <Self as Swap<NodeIndex>>::partition_from(self, NodeIndex(0), |s, n| {
            s.get_neighbor_iterator(n).all(|h| h < left)
        });
        let inside = <Self as Swap<NodeIndex>>::partition_from(self, split, |s, n| {
            s.get_neighbor_iterator(n).any(|h| h < left)
        });
        (left, split, inside)
    }
}

impl<V, P: ForestNodeStore + ForestNodeStorePreorder + Clone> NodeStorageOps for Forest<V, P> {
    type Base = SuBitGraph;
    type OpStorage<N> = Forest<N, P>;

    /// Checks that the live roots' trees partition the half-edges. An absorbed
    /// root's tree hangs below its survivor, so it is covered there.
    fn check_nodes(&self) -> Result<(), crate::half_edge::HedgeGraphError> {
        let n_hedges: Hedge = <Self as Swap<Hedge>>::len(self);
        let mut cover = SuBitGraph::empty(n_hedges.0);

        for i in 0..self.roots.len() {
            let node_id = NodeIndex(i);
            if self.is_absorbed(node_id) {
                continue;
            }
            for h in self.get_neighbor_iterator(node_id) {
                if h >= n_hedges {
                    return Err(crate::half_edge::HedgeGraphError::InvalidHedge(h));
                }
                if cover.includes(&h) {
                    return Err(crate::half_edge::HedgeGraphError::NodesDoNotPartition(
                        format!("They overlap for node {node_id:?}: Cover:{cover:?}, crown:{h:?}"),
                    ));
                }
                cover.add(h);
            }
        }

        let full = !SuBitGraph::empty(n_hedges.0);
        if !(cover.sym_diff(&full).is_empty()) {
            return Err(crate::half_edge::HedgeGraphError::NodesDoNotPartition(
                format!("They do not cover all hedges: Cover:{cover:?}"),
            ));
        }

        Ok(())
    }
    fn extract_nodes(&mut self, nodes: impl IntoIterator<Item = NodeIndex>) -> (SuBitGraph, Self) {
        assert!(!self.has_identification_history(), "{HISTORY_RETAINED}");
        let nodes: HashSet<NodeIndex> = nodes.into_iter().collect();
        let n_hedges: Hedge = <Self as Swap<Hedge>>::len(self);
        let mut extracted = SuBitGraph::empty(n_hedges.0);
        let mut isolated = vec![false; self.roots.len()];
        for node in 0..self.roots.len() {
            let node = NodeIndex(node);
            if nodes.contains(&node) {
                if self.roots[node.0].root_id.is_empty() {
                    isolated[node.0] = true;
                } else {
                    extracted.union_with_iter(self.get_neighbor_iterator(node));
                }
            }
        }

        // Isolated roots have no half-edge representation, so move them out before the
        // crown-based extraction.
        let left = <Self as Swap<NodeIndex>>::partition(self, |n| !isolated[n.0]);
        let isolated_roots = self.roots.split_off(left.0);

        let mut extracted_store = self.extract(
            &extracted,
            |_| panic!("Forest::extract_nodes encountered a split node; expected full nodes only."),
            |d| d,
        );
        extracted_store.roots.extend(isolated_roots);

        (extracted, extracted_store)
    }

    fn iter(&self) -> impl Iterator<Item = (NodeIndex, &Self::NodeData)> {
        self.iter_roots()
            .enumerate()
            .map(|(rid, (d, _))| (NodeIndex(rid), d))
    }

    fn build<
        I: IntoIterator<Item = crate::half_edge::builder::HedgeNodeBuilder<Self::NodeData>>,
    >(
        nodes: I,
        n_hedges: usize,
    ) -> Self {
        Self::build_with_mapping(nodes, n_hedges, std::convert::identity)
    }

    fn build_with_mapping<
        I: IntoIterator<Item = crate::half_edge::builder::HedgeNodeBuilder<ND>>,
        ND,
    >(
        nodes: I,
        n_hedges: usize,
        mut map_data: impl FnMut(ND) -> Self::NodeData,
    ) -> Self {
        // The builder already owns each node's incidence list. Construct the
        // same parent-pointer forest directly instead of allocating a full
        // hedge bitset for every node and decoding those bitsets again.
        let mut entries = Vec::new();
        let mut roots = Vec::new();
        let mut covered = Vec::new();
        for node in nodes {
            let data = map_data(node.data);
            if roots.is_empty() {
                entries.resize_with(n_hedges, || None);
                covered.resize(n_hedges, false);
            } else {
                // Retain the existing partition constructor's admission: its
                // coverage starts after the first partition. Builder-produced
                // incidence is disjoint, but direct trait calls keep the same
                // behavior for repeated or malformed incidence.
                assert!(
                    node.hedges.iter().all(|h| !covered[h.0]),
                    "hedge incidence must form a partition"
                );
                for h in &node.hedges {
                    covered[h.0] = true;
                }
            }
            // Bitset iteration picked the lowest hedge, independent of the
            // order or duplicate entries in the supplied incidence list.
            let first = node.hedges.iter().min().copied().map(TreeNodeId::from);
            let root = RootId(roots.len());
            for h in node.hedges {
                entries[h.0] = Some(if Some(TreeNodeId::from(h)) == first {
                    PPNode::dataless_root(root)
                } else {
                    PPNode::dataless_child(first.unwrap())
                });
            }
            roots.push(RootData {
                root_id: first.unwrap_or(TreeNodeId::EMPTY),
                data,
            });
        }
        Self {
            nodes: P::from_parents(
                entries
                    .into_iter()
                    .map(|entry| entry.expect("hedge incidence must cover all half-edges")),
            ),
            roots,
            merged_into: Vec::new(),
        }
    }

    fn drain(self) -> impl Iterator<Item = (NodeIndex, Self::NodeData)> {
        self.roots
            .into_iter()
            .enumerate()
            .map(|(id, a)| (NodeIndex(id), a.data))
    }

    fn forget_identification_history(&mut self) -> NodeVec<Self::NodeData> {
        // Swapping roots directly: `Swap<NodeIndex>` refuses while history is retained.
        let active = partition_with(
            self,
            NodeIndex(0),
            NodeIndex(self.roots.len()),
            |s, n| s.is_active_root(n.into()),
            |s, a, b| s.swap_roots(a.into(), b.into()),
        );
        if active.0 < self.roots.len() {
            // Once history is discarded, its hierarchy has no remaining owner.
            // Each forgotten root's former tree now hangs below an active root;
            // flattening just those branches keeps future hedge-to-node lookups
            // within one parent hop at a cost proportional to the merged crowns.
            let absorbed = self.roots[active.0..]
                .iter()
                .map(|root| root.root_id)
                .filter(|root_id| !root_id.is_empty());
            self.nodes.flatten_branches(absorbed);
        }
        self.merged_into.clear();
        self.roots
            .split_off(active.0)
            .into_iter()
            .map(|r| r.data)
            .collect()
    }

    fn extend_mut(&mut self, other: Self) {
        let nodeshift = TreeNodeId(self.nodes.n_nodes());
        let shift_roots_by = RootId(self.roots.len());
        if self.has_identification_history() || other.has_identification_history() {
            let mut links = other.merged_into;
            links.resize(other.roots.len(), None);
            self.merged_into.resize(shift_roots_by.0, None);
            self.merged_into.extend(
                links
                    .into_iter()
                    .map(|merged| merged.map(|r| RootId(r.0 + shift_roots_by.0))),
            );
        }
        self.roots.extend(other.roots.into_iter().map(|mut a| {
            if !a.root_id.is_empty() {
                a.root_id.0 += nodeshift.0;
            }
            a
        }));

        self.nodes.extend(other.nodes, shift_roots_by);
    }

    fn extend(mut self, other: Self) -> Self {
        self.extend_mut(other);
        self
    }

    fn random(sources: &[Self::Neighbors], sinks: &[Self::Neighbors]) -> Self
    where
        Self::NodeData: Default,
    {
        let a = NodeStorageVec::<()>::random(sources, sinks)
            .to_forest::<Self::NodeData, P::NodeData>(|_| Default::default());

        Forest {
            roots: a.roots,
            nodes: P::from_pp(a.nodes),
            merged_into: a.merged_into,
        }
    }

    fn delete<S: crate::half_edge::subgraph::SubSetLike<Base = Self::Base>>(
        &mut self,
        subgraph: &S,
    ) {
        let (left, _, overlapping_nodes) = self.partition_for_removal(subgraph);

        // Now all hedges in the subgraph are at the end of the storage (swapped)
        // and all nodes in the subgraph are also at the end of the nodestore
        //
        // We can safely split off the roots after the ones with overlap.
        let _ = self.roots.split_off(overlapping_nodes.0);

        // We now need to adjust the pointer structure that associates a hedge with its node.
        // Splitting it off aswell.
        let _ = self.nodes.split_off(left.into());

        // we need to adjust the root_node_pointer that the root stores, and also shift the root_id that the root node stores.
        //

        for n in self.nodes.iter_node_id() {
            let rid = self.nodes.root_node(n);
            let r = self.nodes.root(rid);

            self.roots[r.0].root_id = rid;
        }

        #[cfg(test)]
        self.nodes.validate().unwrap();
    }

    fn extract<S: crate::half_edge::subgraph::SubSetLike<Base = Self::Base>, V2>(
        &mut self,
        subgraph: &S,
        mut split_node: impl FnMut(&Self::NodeData) -> V2,
        mut owned_node: impl FnMut(Self::NodeData) -> V2,
    ) -> Self::OpStorage<V2> {
        let (left, left_nodes, overlapping_nodes) = self.partition_for_removal(subgraph);

        // Now all hedges in the subgraph are at the end of the storage (swapped)
        // and all nodes in the subgraph are also at the end of the nodestore
        //
        // We can safely split off the roots after the ones with overlap.
        let extracted_roots: Vec<_> = self
            .roots
            .split_off(overlapping_nodes.0)
            .into_iter()
            .map(|a| RootData {
                data: owned_node(a.data),
                root_id: a.root_id,
            })
            .collect();

        let mut overlapping_roots = vec![];

        // For those with overlap, duplicate them with the closure.
        for i in (left_nodes.0)..(overlapping_nodes.0) {
            let data = split_node(&self.roots[i].data);
            // println!("overlapping node found");

            overlapping_roots.push(RootData {
                data,
                root_id: self.roots[i].root_id,
            });
        }

        // We start with the overlapping and then the fully included, so that we do not need to swap
        overlapping_roots.extend(extracted_roots);

        // We now need to adjust the pointer structure that associates a hedge with its node.
        // Splitting it off aswell.
        let mut extracted_nodes = self.nodes.split_off(left.into());

        // we need to adjust the root_node_pointer that the root stores, and also shift the root_id that the root node stores.
        //

        #[cfg(test)]
        extracted_nodes.validate().unwrap();
        // println!("{}", extracted_nodes.debug_draw(|_| None));
        for n in extracted_nodes.iter_node_id() {
            let rid = extracted_nodes.root_node(n);
            let r = extracted_nodes.root(rid);

            overlapping_roots[r.0 - left_nodes.0].root_id = rid;
        }

        for n in self.nodes.iter_node_id() {
            let rid = self.nodes.root_node(n);
            let r = self.nodes.root(rid);

            self.roots[r.0].root_id = rid;
        }

        for (i, r) in overlapping_roots.iter().enumerate() {
            if !r.root_id.is_empty() {
                extracted_nodes.set_root(r.root_id, RootId(i))
            }
        }
        // println!("{}", extracted_nodes.debug_draw(|_| None));
        // println!("{}", self.nodes.debug_draw(|_| None));
        #[cfg(test)]
        self.nodes.validate().unwrap();

        Forest {
            nodes: extracted_nodes,
            roots: overlapping_roots,
            merged_into: Vec::new(),
        }
    }

    fn to_forest<U, H>(
        &self,
        map_data: impl Fn(&Self::NodeData) -> U,
    ) -> Forest<U, ParentPointerStore<H>> {
        Forest {
            nodes: self.nodes.forgetful_map().to_store().into(),
            roots: self
                .roots
                .iter()
                .map(|a| crate::tree::RootData {
                    data: map_data(&a.data),
                    root_id: a.root_id,
                })
                .collect(),
            merged_into: self.merged_into.clone(),
        }
    }

    fn get_neighbor_iterator(&self, node_id: NodeIndex) -> Self::NeighborsIter<'_> {
        let root_id = self.roots[node_id.0].root_id;
        ForestNeighborIter::from_root(&self.nodes, root_id)
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
        Ok(Forest {
            nodes: self.nodes.clone(),
            roots: self
                .roots
                .iter()
                .enumerate()
                .map(|(i, r)| {
                    match node_map(graph, self.get_neighbor_iterator(NodeIndex(i)), &r.data) {
                        Err(err) => Err(err),
                        Ok(data) => Ok(RootData {
                            data,
                            root_id: r.root_id,
                        }),
                    }
                })
                .collect::<Result<Vec<_>, Er>>()?,
            merged_into: self.merged_into.clone(),
        })
    }

    fn check_and_set_nodes(&mut self) -> Result<(), crate::half_edge::HedgeGraphError> {
        Ok(())
    }

    fn iter_nodes(
        &self,
    ) -> impl Iterator<Item = (NodeIndex, Self::NeighborsIter<'_>, &Self::NodeData)> {
        self.roots.iter().enumerate().map(|(i, d)| {
            (
                NodeIndex(i),
                self.get_neighbor_iterator(NodeIndex(i)),
                &d.data,
            )
        })
    }

    fn node_id_ref(&self, hedge: Hedge) -> NodeIndex {
        self.root(hedge.into()).into()
    }

    fn get_node_data(&self, node_id: NodeIndex) -> &Self::NodeData {
        &self.roots[node_id.0].data
    }

    fn identify_nodes(
        &mut self,
        nodes: &[NodeIndex],
        node_data_merge: Self::NodeData,
    ) -> NodeIndex {
        let survivor = self.active_root(nodes[0].into());
        for &node in &nodes[1..] {
            // Resolving after the previous merges also skips repeated nodes.
            let absorbed = self.active_root(node.into());
            if absorbed == survivor {
                continue;
            }
            self.merge_into(absorbed, survivor);
            let tree = self.roots[absorbed.0].root_id;
            if tree.is_empty() {
                continue;
            }
            let survivor_tree = self.roots[survivor.0].root_id;
            if survivor_tree.is_empty() {
                // The survivor takes over the whole tree, so the absorbed root
                // keeps no crown of its own.
                self.nodes.set_root(tree, survivor);
                self.roots[survivor.0].root_id = tree;
                self.roots[absorbed.0].root_id = TreeNodeId::EMPTY;
            } else {
                self.nodes.make_root_of(survivor_tree, tree);
            }
        }
        self.roots[survivor.0].data = node_data_merge;
        survivor.into()
    }

    fn is_absorbed(&self, node: NodeIndex) -> bool {
        !self.is_active_root(node.into())
    }

    fn iter_nodes_mut(
        &mut self,
    ) -> impl Iterator<Item = (NodeIndex, Self::NeighborsIter<'_>, &mut Self::NodeData)> {
        self.roots.iter_mut().enumerate().map(|(i, d)| {
            let root_id = d.root_id;
            (
                NodeIndex(i),
                ForestNeighborIter::from_root(&self.nodes, root_id),
                &mut d.data,
            )
        })
    }

    fn add_dangling_edge_mut(
        &mut self,
        source: NodeIndex,
    ) -> Result<(), crate::half_edge::HedgeGraphError> {
        if source.0 >= self.roots.len() {
            return Err(crate::half_edge::HedgeGraphError::NoNode);
        }
        if !self.is_active_root(source.into()) {
            return Err(crate::half_edge::HedgeGraphError::InvalidNode(source));
        }
        let root = &mut self.roots[source.0];
        if root.root_id.is_empty() {
            root.root_id = self.nodes.add_dataless_root(source.into());
        } else {
            self.nodes.add_dataless_child(root.root_id);
        }
        Ok(())
    }

    fn get_node_data_mut(&mut self, node_id: NodeIndex) -> &mut Self::NodeData {
        &mut self.roots[node_id.0].data
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
        Forest {
            roots: self
                .roots
                .into_iter()
                .enumerate()
                .map(|(i, r)| RootData {
                    data: f(involution, NodeIndex(i), r.data),
                    root_id: r.root_id,
                })
                .collect(),
            nodes: self.nodes,
            merged_into: self.merged_into,
        }
    }

    fn map_data_graph_result<'a, V2, Err>(
        self,
        involution: &'a crate::half_edge::involution::Involution<
            crate::half_edge::involution::EdgeIndex,
        >,
        mut f: impl FnMut(
            &'a crate::half_edge::involution::Involution<crate::half_edge::involution::EdgeIndex>,
            NodeIndex,
            Self::NodeData,
        ) -> Result<V2, Err>,
    ) -> Result<Self::OpStorage<V2>, Err> {
        Ok(Forest {
            roots: self
                .roots
                .into_iter()
                .enumerate()
                .map(|(i, r)| match f(involution, NodeIndex(i), r.data) {
                    Ok(data) => Ok(RootData {
                        data,
                        root_id: r.root_id,
                    }),
                    Err(e) => Err(e),
                })
                .collect::<Result<_, Err>>()?,
            nodes: self.nodes,
            merged_into: self.merged_into,
        })
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
        Forest {
            nodes: self.nodes.clone(),
            roots: self
                .roots
                .iter()
                .enumerate()
                .map(|(i, r)| RootData {
                    data: node_map(graph, self.get_neighbor_iterator(NodeIndex(i)), &r.data),
                    root_id: r.root_id,
                })
                .collect(),
            merged_into: self.merged_into.clone(),
        }
    }

    fn new_nodevec<'a, V2>(
        &'a self,
        mut node_map: impl FnMut(NodeIndex, Self::NeighborsIter<'a>, &'a Self::NodeData) -> V2,
    ) -> NodeVec<V2> {
        self.roots
            .iter()
            .enumerate()
            .map(|(i, r)| {
                node_map(
                    NodeIndex(i),
                    self.get_neighbor_iterator(NodeIndex(i)),
                    &r.data,
                )
            })
            .collect()
    }

    fn map_data_ref_mut_graph<'a, V2>(
        &'a mut self,
        mut node_map: impl FnMut(Self::NeighborsIter<'a>, &'a mut Self::NodeData) -> V2,
    ) -> Self::OpStorage<V2> {
        Forest {
            nodes: self.nodes.clone(),
            roots: self
                .roots
                .iter_mut()
                // .enumerate()
                .map(|r| {
                    let root_id = r.root_id;

                    RootData {
                        data: node_map(
                            ForestNeighborIter::from_root(&self.nodes, root_id),
                            &mut r.data,
                        ),
                        root_id: r.root_id,
                    }
                })
                .collect(),
            merged_into: self.merged_into.clone(),
        }
    }
}

// impl<V,P:ForestNodeStoreDown> Forest<V,P>{
//     fn hedge_node(&self,root:RootId)->HedgeNode{

//     }
// }

// impl<V, P: ForestNodeStoreDown> NodeStorageOps for Forest<V, P> {
//     type Base = BitVec;
//     type OpStorage<N> = Self::Storage<N>;

//     fn iter(&self) -> impl Iterator<Item = (crate::half_edge::NodeIndex, &Self::NodeData)> {
//         self.iter_roots()
//             .enumerate()
//             .map(|(id, (n, hedge))| (NodeIndex::from(id), n))
//     }

//     fn iter_nodes(
//         &self,
//     ) -> impl Iterator<
//         Item = (
//             &crate::half_edge::subgraph::HedgeNode,
//             NodeIndex,
//             &Self::NodeData,
//         ),
//     > {
//     }
// }

#[cfg(test)]
mod builder_tests {
    use super::*;
    use crate::{
        half_edge::{builder::HedgeNodeBuilder, subgraph::BaseSubgraph},
        tree::{child_pointer::ParentChildStore, child_vec::ChildVecStore},
    };

    fn builders(parts: &[Vec<usize>]) -> Vec<HedgeNodeBuilder<usize>> {
        parts
            .iter()
            .enumerate()
            .map(|(data, hedges)| HedgeNodeBuilder {
                data,
                hedges: hedges.iter().copied().map(Hedge).collect(),
            })
            .collect()
    }

    fn reference<P: ForestNodeStore>(parts: &[Vec<usize>], n_hedges: usize) -> Forest<usize, P> {
        Forest::from_bitvec_partition(builders(parts).into_iter().map(|node| {
            (
                node.data,
                SuBitGraph::from_hedge_iter(node.hedges.into_iter(), n_hedges),
            )
        }))
        .unwrap()
    }

    #[test]
    fn discarding_history_flattens_crowns_without_changing_graph_data() {
        use crate::{
            half_edge::{builder::HedgeGraphBuilder, involution::Flow, HedgeGraph},
            tree::parent_pointer::ParentId,
        };

        fn check<P: ForestNodeStore<NodeData = ()> + ForestNodeStorePreorder + Clone>() {
            let mut builder = HedgeGraphBuilder::<usize, usize, usize>::new();
            let nodes: Vec<_> = (0..34).map(|i| builder.add_node(i)).collect();
            for i in 0..32 {
                builder.add_edge(
                    nodes[i].add_data(2 * i),
                    nodes[(i + 1) % 32].add_data(2 * i + 1),
                    i,
                    i % 2 == 0,
                );
            }
            builder.add_edge(nodes[0].add_data(64), nodes[0].add_data(65), 32, false);
            builder.add_edge(nodes[0].add_data(66), nodes[1].add_data(67), 33, true);
            builder.add_external_edge(nodes[1].add_data(68), 34, true, Flow::Sink);
            let mut graph: HedgeGraph<_, _, _, Forest<_, P>> = builder.build();
            for group in 0..2 {
                let mut previous = nodes[group];
                for layer in 1..16 {
                    previous = graph.identify_nodes(&[nodes[group + 2 * layer], previous], group);
                }
            }
            // Include an identified isolated node and a surviving isolated node.
            graph.identify_nodes(&[nodes[30], nodes[32]], 0);
            let snapshot = |graph: &HedgeGraph<usize, usize, usize, Forest<usize, P>>| {
                graph
                    .iter_hedges()
                    .map(|(hedge, data)| {
                        (
                            hedge,
                            *data,
                            graph[graph.node_id(hedge)],
                            graph.inv(hedge),
                            graph.flow(hedge),
                            graph[&hedge],
                            graph[graph[&hedge]],
                        )
                    })
                    .collect::<Vec<_>>()
            };
            let before = snapshot(&graph);
            let crowns: Vec<_> = [nodes[30], nodes[31], nodes[33]]
                .into_iter()
                .map(|node| (graph[node], graph.iter_crown(node).collect::<Vec<_>>()))
                .collect();
            let removed = graph.forget_identification_history();
            assert_eq!(removed.iter().count(), 31);
            assert_eq!(snapshot(&graph), before);
            assert_eq!(graph.n_nodes(), 3);
            graph.check().unwrap();
            for (node, crown, payload) in graph.iter_nodes() {
                assert_eq!(
                    crown.collect::<Vec<_>>(),
                    crowns.iter().find(|(p, _)| p == payload).unwrap().1
                );
                for hedge in graph.iter_crown(node) {
                    if let ParentId::Node(parent) = graph.node_store.nodes.parent(hedge.into()) {
                        assert!(matches!(
                            graph.node_store.nodes.parent(parent),
                            ParentId::Root(_)
                        ));
                    }
                }
            }
            assert_eq!(graph.forget_identification_history().iter().count(), 0);
            assert_eq!(snapshot(&graph), before);
            // Mutations and subgraph extraction still work after flattening.
            let a = graph.node_id(Hedge(0));
            let b = graph.node_id(Hedge(1));
            graph = graph
                .add_pair(a.add_data(69), b.add_data(70), 35, true)
                .unwrap()
                .2;
            graph.check().unwrap();
            let (_, extracted) = graph.node_store.extract_nodes([a]);
            extracted.check_nodes().unwrap();
            graph.node_store.check_nodes().unwrap();
        }
        check::<ParentChildStore<()>>();
        check::<ChildVecStore<()>>();
    }

    #[test]
    fn batched_forgetting_flattens_only_merged_trees() {
        use crate::{
            half_edge::{builder::HedgeGraphBuilder, HedgeGraph, NoData},
            tree::parent_pointer::ParentId,
        };
        use std::collections::HashMap;

        type Graph<P> = HedgeGraph<usize, usize, NoData, Forest<usize, P>>;
        type Crowns = HashMap<usize, Vec<(Hedge, Option<TreeNodeId>)>>;

        // Active crowns in order with each half-edge's parent node, by payload.
        // Root labels are omitted because compaction renumbers nodes.
        fn crowns<P: ForestNodeStorePreorder<NodeData = ()> + Clone>(graph: &Graph<P>) -> Crowns {
            let store = &graph.node_store;
            graph
                .iter_nodes()
                .filter(|(node, _, _)| store.is_active_root((*node).into()))
                .map(|(_, crown, &payload)| {
                    let parents = crown.map(|h| match store.nodes.parent(h.into()) {
                        ParentId::Node(parent) => (h, Some(parent)),
                        _ => (h, None),
                    });
                    (payload, parents.collect())
                })
                .collect()
        }

        fn check<P: ForestNodeStorePreorder<NodeData = ()> + Clone>() {
            let n = 24;
            let mut builder = HedgeGraphBuilder::new();
            let nodes: Vec<_> = (0..n).map(|i| builder.add_node(i)).collect();
            for i in 0..n {
                builder.add_edge(nodes[i], nodes[(i + 1) % n], i, false);
                builder.add_edge(nodes[i], nodes[(i + 5) % n], n + i, false);
            }
            let mut graph: Graph<P> = builder.build();
            let reps: Vec<_> = (0..n)
                .map(|i| graph.iter_crown(NodeIndex(i)).next().unwrap())
                .collect();
            let mut payload = n;
            let mut merge = |graph: &mut Graph<P>, a: usize, b: usize| {
                payload += 1;
                let nodes = [graph.node_id(reps[a]), graph.node_id(reps[b])];
                graph.identify_nodes(&nodes, payload);
                payload
            };
            // Node 23 is never merged. Rerooting deepens its tree, which
            // forgetting merge history leaves as it is.
            let rerooted = 23;
            let new_root = graph.iter_crown(NodeIndex(rerooted)).nth(2).unwrap();
            let store = &mut graph.node_store;
            let mut parents: ParentPointerStore<()> = store.nodes.clone().to_store().into();
            let root = parents.change_root(new_root.into());
            store.nodes = P::from_pp(parents);
            store.roots[root.0].root_id = new_root.into();
            let mut previous = crowns(&graph);
            let is_flat = |crown: &[(Hedge, Option<TreeNodeId>)]| {
                let root = Some(crown[0].0.into());
                crown[0].1.is_none() && crown[1..].iter().all(|(_, parent)| *parent == root)
            };
            assert!(!is_flat(&previous[&rerooted]));
            for batch in 0..4 {
                let [a, b, c, d, e] = [1, 2, 3, 4, 5].map(|i| 5 * batch + i);
                merge(&mut graph, 0, a);
                // Nest absorbed trees three levels deep below the hub.
                merge(&mut graph, b, c);
                merge(&mut graph, d, b);
                let mut hub = merge(&mut graph, 0, d);
                if batch % 2 == 1 {
                    // Move the whole hub below a fresh root.
                    hub = merge(&mut graph, e, 0);
                }
                let mut merged = vec![hub];
                if batch == 0 {
                    // A separate merged node, left untouched by later batches.
                    merged.push(merge(&mut graph, 21, 22));
                }
                let before = crowns(&graph);
                graph.forget_identification_history();
                assert_eq!(graph.n_nodes(), before.len());
                graph.check().unwrap();
                graph.node_store.check_nodes().unwrap();
                let after = crowns(&graph);
                assert_eq!(after.len(), before.len());
                for (payload, crown) in &after {
                    let hedges =
                        |crown: &[(Hedge, _)]| crown.iter().map(|(h, _)| *h).collect::<Vec<_>>();
                    assert_eq!(hedges(crown), hedges(&before[payload]));
                    // Every merged crown is flat: preorder starts at the root.
                    assert_eq!(is_flat(crown), *payload != rerooted);
                    if !merged.contains(payload) {
                        assert_eq!(crown, &previous[payload]);
                    }
                }
                previous = after;
            }
            // The hub, the separate merge, node 23 and the unused nodes 5 and 15.
            assert_eq!(graph.n_nodes(), 5);
        }
        check::<ParentChildStore<()>>();
        check::<ChildVecStore<()>>();
    }

    #[test]
    fn direct_forest_builder_preserves_parent_root_and_crown_order() {
        for (parts, n_hedges) in [
            (vec![], 0),
            (vec![], 3),
            (vec![vec![]], 0),
            (vec![vec![], vec![], vec![]], 0),
            (vec![vec![3, 0, 1, 2]], 4),
            (vec![vec![3, 1], vec![], vec![2, 0]], 4),
            (vec![vec![2, 2, 0], vec![3, 1, 3]], 4),
            (vec![vec![], vec![0], vec![]], 1),
        ] {
            macro_rules! compare {
                ($store:ty) => {{
                    let actual = <Forest<usize, $store> as NodeStorageOps>::build(
                        builders(&parts),
                        n_hedges,
                    );
                    let expected = reference::<$store>(&parts, n_hedges);
                    assert_eq!(actual, expected);
                    let actual_parents: ParentPointerStore<()> = actual.nodes.clone().into();
                    let expected_parents =
                        reference::<ParentPointerStore<()>>(&parts, n_hedges).nodes;
                    assert_eq!(actual_parents, expected_parents);
                    assert_eq!(actual.nodes.n_nodes(), expected.nodes.n_nodes());
                    for i in 0..parts.len() {
                        assert_eq!(
                            actual
                                .get_neighbor_iterator(NodeIndex(i))
                                .collect::<Vec<_>>(),
                            expected
                                .get_neighbor_iterator(NodeIndex(i))
                                .collect::<Vec<_>>()
                        );
                    }
                }};
            }
            compare!(ParentChildStore<()>);
            compare!(ChildVecStore<()>);
        }
    }

    #[test]
    fn direct_forest_builder_retains_partition_failure_boundary() {
        for (parts, n_hedges) in [
            (vec![vec![1]], 1),
            (vec![vec![0]], 2),
            (vec![vec![]], 1),
            (vec![vec![0], vec![1], vec![1]], 2),
            // The old partition constructor does not add the first partition
            // to its overlap mask. Preserve, rather than fix, that behavior.
            (vec![vec![0], vec![0]], 1),
        ] {
            let actual = std::panic::catch_unwind(|| {
                <Forest<usize, ParentChildStore<()>> as NodeStorageOps>::build(
                    builders(&parts),
                    n_hedges,
                )
            });
            let expected =
                std::panic::catch_unwind(|| reference::<ParentChildStore<()>>(&parts, n_hedges));
            assert_eq!(actual.is_ok(), expected.is_ok());
            if let (Ok(actual), Ok(expected)) = (actual, expected) {
                assert_eq!(actual, expected);
            }
        }
    }

    #[test]
    fn direct_forest_builder_maps_nonclone_payloads_once_in_node_order() {
        use std::sync::Mutex;
        let nodes = (0..3).map(|i| HedgeNodeBuilder {
            data: Mutex::new(i),
            hedges: vec![Hedge(i)],
        });
        let mut calls = Vec::new();
        let actual =
            <Forest<Mutex<usize>, ParentChildStore<()>> as NodeStorageOps>::build_with_mapping(
                nodes,
                3,
                |data| {
                    let value = data.into_inner().unwrap();
                    calls.push(value);
                    Mutex::new(value + 10)
                },
            );
        assert_eq!(calls, vec![0, 1, 2]);
        for (i, root) in actual.roots.iter().enumerate() {
            assert_eq!(*root.data.lock().unwrap(), i + 10);
            assert_eq!(root.root_id, TreeNodeId(i));
        }
    }
}

#[cfg(test)]
mod identification_tests {
    use std::{
        collections::BTreeSet,
        panic::{catch_unwind, AssertUnwindSafe},
    };

    use super::*;
    use crate::{
        half_edge::{builder::HedgeGraphBuilder, involution::Flow, HedgeGraph, HedgeGraphError},
        tree::{child_pointer::ParentChildStore, child_vec::ChildVecStore},
    };

    trait Store: ForestNodeStore<NodeData = ()> + ForestNodeStorePreorder + Clone {}
    impl<P: ForestNodeStore<NodeData = ()> + ForestNodeStorePreorder + Clone> Store for P {}

    type Graph<P> = HedgeGraph<(), usize, (), Forest<usize, P>>;

    /// A four-cycle on nodes 0..4 and two isolated nodes 4 and 5.
    fn graph<P: Store>() -> Graph<P> {
        let mut builder = HedgeGraphBuilder::new();
        let nodes: Vec<_> = (0..6).map(|i| builder.add_node(i)).collect();
        for i in 0..4 {
            builder.add_edge(nodes[i], nodes[(i + 1) % 4], (), true);
        }
        builder.build()
    }

    fn crown<P: Store>(graph: &Graph<P>, node: usize) -> BTreeSet<Hedge> {
        graph.iter_crown(NodeIndex(node)).collect()
    }

    macro_rules! for_both_stores {
        ($check:ident) => {
            $check::<ChildVecStore<()>>();
            $check::<ParentChildStore<()>>();
        };
    }

    #[test]
    fn identifying_absorbed_nodes_merges_their_survivors() {
        fn check<P: Store>() {
            let mut graph = graph::<P>();
            let original: Vec<_> = (0..6).map(|node| crown(&graph, node)).collect();
            let n = NodeIndex;
            assert_eq!(graph.identify_nodes(&[n(0), n(1)], 10), n(0));
            // Repeated and absorbed nodes resolve to the node that absorbed them.
            assert_eq!(graph.identify_nodes(&[n(1), n(0), n(1)], 11), n(0));
            assert_eq!(graph.identify_nodes(&[n(3), n(3)], 12), n(3));
            assert_eq!(graph.identify_nodes(&[n(2), n(1)], 13), n(2));
            // Nodes without half-edges keep their survivor as well.
            assert_eq!(graph.identify_nodes(&[n(4), n(5)], 14), n(4));
            assert_eq!(graph.identify_nodes(&[n(5), n(3)], 15), n(4));
            assert_eq!(graph.identify_nodes(&[n(3), n(1), n(5)], 16), n(4));
            for hedge in original.iter().flatten() {
                assert_eq!(graph.node_id(*hedge), n(4));
            }
            assert_eq!(graph[n(4)], 16);
            assert_eq!(
                crown(&graph, 4),
                original.iter().flatten().copied().collect::<BTreeSet<_>>()
            );

            let mut removed: Vec<_> = graph
                .forget_identification_history()
                .into_iter()
                .map(|(_, data)| data)
                .collect();
            removed.sort();
            assert_eq!(removed, [1, 5, 11, 12, 13]);
            assert_eq!(graph.n_nodes(), 1);
            assert_eq!(graph[n(0)], 16);
            graph.check().unwrap();
        }
        for_both_stores!(check);
    }

    #[test]
    fn an_isolated_survivor_does_not_share_its_crown() {
        fn check<P: Store>() {
            let mut graph = graph::<P>();
            let absorbed = crown(&graph, 0);
            assert_eq!(
                graph.identify_nodes(&[NodeIndex(4), NodeIndex(0)], 10),
                NodeIndex(4)
            );
            assert_eq!(crown(&graph, 4), absorbed);
            assert!(crown(&graph, 0).is_empty());
            for hedge in &absorbed {
                assert_eq!(graph.node_id(*hedge), NodeIndex(4));
            }
            assert_eq!(
                graph
                    .forget_identification_history()
                    .into_iter()
                    .map(|(_, data)| data)
                    .collect::<Vec<_>>(),
                [0]
            );
            graph.check().unwrap();
            assert_eq!(graph.n_nodes(), 5);
        }
        for_both_stores!(check);
    }

    #[test]
    fn dangling_edges_need_an_existing_node() {
        fn check<P: Store>() {
            assert!(matches!(
                graph::<P>().node_store.add_dangling_edge(NodeIndex(6)),
                Err(HedgeGraphError::NoNode)
            ));
        }
        for_both_stores!(check);
    }

    #[test]
    fn dangling_edges_need_an_active_node() {
        fn check<P: Store>() {
            let mut graph = graph::<P>();
            graph.identify_nodes(&[NodeIndex(0), NodeIndex(1)], 10);
            graph.identify_nodes(&[NodeIndex(4), NodeIndex(5)], 11);
            for absorbed in [1, 5] {
                assert!(matches!(
                    graph.node_store.clone().add_dangling_edge(NodeIndex(absorbed)),
                    Err(HedgeGraphError::InvalidNode(NodeIndex(node))) if node == absorbed
                ));
            }
            for survivor in [0, 4] {
                let store = graph
                    .node_store
                    .clone()
                    .add_dangling_edge(NodeIndex(survivor));
                let store = store.unwrap();
                assert_eq!(store.node_id_ref(Hedge(8)), NodeIndex(survivor));
            }
        }
        for_both_stores!(check);
    }

    #[test]
    fn moving_half_edges_or_nodes_requires_forgotten_history() {
        fn check<P: Store>() {
            let mut graph = graph::<P>();
            graph.identify_nodes(&[NodeIndex(1), NodeIndex(0)], 10);
            let mut subgraph: SuBitGraph = graph.empty_subgraph();
            subgraph.add(Hedge(0));
            for operation in 0..5 {
                let run = |store: &mut Forest<usize, P>| match operation {
                    0 => store.delete(&subgraph),
                    1 => drop(store.extract(&subgraph, |d| *d, |d| d)),
                    2 => drop(store.extract_nodes([NodeIndex(2)])),
                    3 => Swap::<Hedge>::swap(store, Hedge(0), Hedge(7)),
                    _ => Swap::<NodeIndex>::swap(store, NodeIndex(2), NodeIndex(3)),
                };
                let mut store = graph.node_store.clone();
                assert!(catch_unwind(AssertUnwindSafe(|| run(&mut store))).is_err());
                let mut store = graph.node_store.clone();
                store.forget_identification_history();
                run(&mut store);
                store.check_nodes().unwrap();
            }
        }
        for_both_stores!(check);
    }

    #[test]
    fn appending_a_store_shifts_its_identification_links() {
        fn check<P: Store>() {
            let mut appended = graph::<P>();
            appended.identify_nodes(&[NodeIndex(1), NodeIndex(0)], 10);
            for self_history in [false, true] {
                let mut store = graph::<P>().node_store;
                if self_history {
                    store.identify_nodes(&[NodeIndex(2), NodeIndex(3)], 11);
                }
                store.extend_mut(appended.node_store.clone());
                store.validate_structure().unwrap();
                assert_eq!(
                    store.identify_nodes(&[NodeIndex(6), NodeIndex(8)], 12),
                    NodeIndex(7)
                );
                assert!(store.is_active_root(RootId(2)));
                assert_eq!(store.is_active_root(RootId(3)), !self_history);
                assert_eq!(
                    store.forget_identification_history().len(),
                    NodeIndex(2 + usize::from(self_history))
                );
                assert!(!store.has_identification_history());
                store.check_nodes().unwrap();
            }
        }
        for_both_stores!(check);
    }

    /// Random identifications of arbitrary, possibly absorbed or repeated nodes,
    /// dangling edges and forgetting, against a union-find over node indices.
    #[test]
    fn random_identifications_match_a_union_find() {
        fn check<P: Store>() {
            let mut state = 0x9e37_79b9_7f4a_7c15_u64;
            let mut random = |bound: usize| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                (state % bound as u64) as usize
            };
            let find = |parent: &[usize], mut node: usize| {
                while parent[node] != node {
                    node = parent[node];
                }
                node
            };
            for _ in 0..300 {
                let mut builder = HedgeGraphBuilder::new();
                let n_nodes = 1 + random(7);
                let nodes: Vec<_> = (0..n_nodes).map(|i| builder.add_node(i)).collect();
                for _ in 0..random(10) {
                    builder.add_edge(nodes[random(n_nodes)], nodes[random(n_nodes)], (), true);
                }
                let mut graph: Graph<P> = builder.build();
                let node_ids = |graph: &Graph<P>| -> Vec<usize> {
                    (0..graph.n_hedges())
                        .map(|hedge| graph.node_id(Hedge(hedge)).0)
                        .collect()
                };
                // Node payloads stay unique, so they identify nodes across forgetting.
                let mut parent: Vec<usize> = (0..n_nodes).collect();
                let mut owner = node_ids(&graph);
                for step in 0..12 {
                    match random(4) {
                        0 | 1 => {
                            let requested: Vec<_> =
                                (0..=random(4)).map(|_| random(parent.len())).collect();
                            let survivor = find(&parent, requested[0]);
                            for &node in &requested {
                                let node = find(&parent, node);
                                parent[node] = survivor;
                            }
                            let requested: Vec<_> = requested.into_iter().map(NodeIndex).collect();
                            assert_eq!(
                                graph.identify_nodes(&requested, 100 * (step + 1)),
                                NodeIndex(survivor)
                            );
                        }
                        2 => {
                            let node = random(parent.len());
                            match graph.clone().add_dangling_edge(
                                NodeIndex(node),
                                (),
                                Flow::Sink,
                                true,
                            ) {
                                Ok((_, extended)) => {
                                    assert_eq!(find(&parent, node), node);
                                    graph = extended;
                                    owner.push(node);
                                }
                                Err(HedgeGraphError::InvalidNode(_)) => {
                                    assert_ne!(find(&parent, node), node)
                                }
                                Err(error) => panic!("{error}"),
                            }
                        }
                        _ => {
                            let payload =
                                |graph: &Graph<P>, node| graph[NodeIndex(find(&parent, node))];
                            let survivors: BTreeSet<_> = (0..parent.len())
                                .map(|node| payload(&graph, node))
                                .collect();
                            let expected: Vec<_> =
                                owner.iter().map(|&node| payload(&graph, node)).collect();
                            graph.forget_identification_history();
                            graph.check().unwrap();
                            let nodes = graph.iter_nodes().map(|(_, _, data)| *data);
                            assert_eq!(nodes.collect::<BTreeSet<_>>(), survivors);
                            parent = (0..graph.n_nodes()).collect();
                            owner = node_ids(&graph);
                            let actual: Vec<_> =
                                owner.iter().map(|&node| graph[NodeIndex(node)]).collect();
                            assert_eq!(actual, expected);
                        }
                    }
                    for (hedge, &node) in owner.iter().enumerate() {
                        assert_eq!(graph.node_id(Hedge(hedge)), NodeIndex(find(&parent, node)));
                    }
                }
            }
        }
        for_both_stores!(check);
    }
}
