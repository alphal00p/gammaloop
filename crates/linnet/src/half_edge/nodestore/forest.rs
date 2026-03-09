use super::{NodeStorage, NodeStorageOps, NodeStorageVec};
use crate::{
    half_edge::{
        involution::Hedge,
        subgraph::{BaseSubgraph, Inclusion, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
        swap::Swap,
        NodeIndex, NodeVec,
    },
    tree::{
        parent_pointer::ParentPointerStore, Forest, ForestNodeStore, ForestNodeStorePreorder,
        RootData, RootId, TreeNodeId,
    },
};

// pub struct ForestNodeStore<Tree, N> {
//     forest: Forest<N, Tree>,
// }

#[derive(Debug, Clone)]
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
pub struct ForestNeighborIter<'a, P: ForestNodeStorePreorder + 'a> {
    // root: Option<TreeNodeId>, // Potentially for future use or alternative iteration strategies
    /// The underlying pre-order iterator from the forest's node store.
    iter: P::Iterator<'a>,
}

impl<'a, P: ForestNodeStorePreorder + 'a> Iterator for ForestNeighborIter<'a, P> {
    type Item = Hedge;
    fn next(&mut self) -> Option<Self::Item> {
        // if let Some(root) = self.root.take() {
        //     Some(root.into())
        // } else {
        self.iter.next().map(|a| a.into())
        // }
    }
}

impl<'a, P: ForestNodeStorePreorder + 'a> ExactSizeIterator for ForestNeighborIter<'a, P> {
    fn len(&self) -> usize {
        self.iter.clone().count()
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
        if i == j {
            return;
        }
        self.swap_roots(RootId(i.0), RootId(j.0));
    }
}

impl<V, P: ForestNodeStore + ForestNodeStorePreorder + Clone> NodeStorageOps for Forest<V, P> {
    type Base = SuBitGraph;
    type OpStorage<N> = Forest<N, P>;

    fn check_nodes(&self) -> Result<(), crate::half_edge::HedgeGraphError> {
        let n_hedges: Hedge = <Self as Swap<Hedge>>::len(self);
        let mut cover = SuBitGraph::empty(n_hedges.0);

        for (i, _) in self.roots.iter().enumerate() {
            let node_id = NodeIndex(i);
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
        let nodes: Vec<NodeIndex> = nodes.into_iter().collect();
        let n_hedges: Hedge = <Self as Swap<Hedge>>::len(self);
        let mut extracted = SuBitGraph::empty(n_hedges.0);
        for node in nodes {
            extracted.union_with_iter(self.get_neighbor_iterator(node));
        }

        let extracted_store = self.extract(
            &extracted,
            |_| panic!("Forest::extract_nodes encountered a split node; expected full nodes only."),
            |d| d,
        );

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
        Forest::from_bitvec_partition(nodes.into_iter().map(|n| {
            (
                n.data,
                SuBitGraph::from_hedge_iter(n.hedges.into_iter(), n_hedges),
            )
        }))
        .unwrap()
    }

    fn build_with_mapping<
        I: IntoIterator<Item = crate::half_edge::builder::HedgeNodeBuilder<ND>>,
        ND,
    >(
        nodes: I,
        n_hedges: usize,
        mut map_data: impl FnMut(ND) -> Self::NodeData,
    ) -> Self {
        Forest::from_bitvec_partition(nodes.into_iter().map(|n| {
            (
                map_data(n.data),
                SuBitGraph::from_hedge_iter(n.hedges.into_iter(), n_hedges),
            )
        }))
        .unwrap()
    }

    fn drain(self) -> impl Iterator<Item = (NodeIndex, Self::NodeData)> {
        self.roots
            .into_iter()
            .enumerate()
            .map(|(id, a)| (NodeIndex(id), a.data))
    }

    fn forget_identification_history(&mut self) -> NodeVec<Self::NodeData> {
        let mut active_nodes_upper_bound = NodeIndex(0);
        let mut historical_nodes_lower_bound = self.len();

        // first we swap to the front all nodes that have correct pointers.
        while active_nodes_upper_bound < historical_nodes_lower_bound {
            if self
                .nodes
                .root_node(self.roots[active_nodes_upper_bound.0].root_id)
                == self.roots[active_nodes_upper_bound.0].root_id
            {
                //left is in the right place

                active_nodes_upper_bound.0 += 1;
            } else {
                //left needs to be swapped
                historical_nodes_lower_bound.0 -= 1;
                if self
                    .nodes
                    .root_node(self.roots[historical_nodes_lower_bound.0].root_id)
                    == self.roots[historical_nodes_lower_bound.0].root_id
                {
                    //only with an extracted that is in the wrong spot
                    self.swap_roots(
                        active_nodes_upper_bound.into(),
                        historical_nodes_lower_bound.into(),
                    );
                    active_nodes_upper_bound.0 += 1;
                }
            }
        }
        self.roots
            .split_off(active_nodes_upper_bound.0)
            .into_iter()
            .map(|r| r.data)
            .collect()
    }

    fn extend_mut(&mut self, other: Self) {
        let nodeshift = TreeNodeId(self.nodes.n_nodes());
        let shift_roots_by = RootId(self.roots.len());
        self.roots.extend(other.roots.into_iter().map(|mut a| {
            a.root_id.0 += nodeshift.0;
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
        }
    }

    fn delete<S: crate::half_edge::subgraph::SubSetLike<Base = Self::Base>>(
        &mut self,
        subgraph: &S,
    ) {
        let mut left = Hedge(0);
        let mut extracted = self.len();
        // println!("{}", self.nodes.debug_draw(|_| None));
        // Do the same swapping as for the edge store, so that they line up
        while left < extracted {
            // println!("{left},{extracted}");
            if !subgraph.includes(&left) {
                //left is in the right place
                left.0 += 1;
            } else {
                //left needs to be swapped
                extracted.0 -= 1;
                if !subgraph.includes(&extracted) {
                    //only with an extracted that is in the wrong spot
                    // println!("{left}<->{extracted}");
                    self.swap(left, extracted);
                    left.0 += 1;
                }
            }
        }
        // println!("{}", self.nodes.debug_draw(|_| None));

        // Now need to partition the nodes into 3 ranges,
        // - 0..left_nodes do not contain any half edges in the subgraph,
        // - left_nodes..overlapping_nodes are all nodes who have incident half-edges inside the subgraph, but not all of them. These ones need to be split
        // - overlapping_nodes..len are the nodes containing only half-edges in the subgraph
        let mut left_nodes = NodeIndex(0);
        let mut extracted_nodes = self.len();

        // first we swap to the front all nodes that contain no edges from subgraph. Since we have swapped all edges in the subgraph to be after left, in the above, now we check if the neighbor iter contains no hedge>=left
        while left_nodes < extracted_nodes {
            if self.get_neighbor_iterator(left_nodes).all(|h| h < left) {
                //left is in the right place

                left_nodes.0 += 1;
            } else {
                //left needs to be swapped
                extracted_nodes.0 -= 1;
                if self
                    .get_neighbor_iterator(extracted_nodes)
                    .all(|h| h < left)
                {
                    //only with an extracted that is in the wrong spot
                    self.swap_roots(left_nodes.into(), extracted_nodes.into());
                    left_nodes.0 += 1;
                }
            }
        }

        let mut overlapping_nodes = left_nodes;
        let mut non_overlapping_extracted = self.len();

        // println!(
        //     "RootsplitnonOverlapping:{}",
        //     self.nodes.debug_draw(|_| None)
        // );

        // now we place all overlapping nodes after all nodes totally outside. We can check whether a node is overlapping if it now contains a hedge<left.
        while overlapping_nodes < non_overlapping_extracted {
            // println!("looking at possibly overlapping:{overlapping_nodes}");
            if self
                .get_neighbor_iterator(overlapping_nodes)
                .any(|h| h < left)
            {
                //overlapping is in the right place, as it intersects (is after left_nodes) but isn't fully included
                overlapping_nodes.0 += 1;
                // println!("found overlap");
            } else {
                //overlapping needs to be swapped
                non_overlapping_extracted.0 -= 1;
                if self
                    .get_neighbor_iterator(non_overlapping_extracted)
                    .any(|h| h < left)
                {
                    //only with an extracted that is in the wrong spot
                    self.swap_roots(overlapping_nodes.into(), non_overlapping_extracted.into());
                    // println!("swap overlapping");
                    overlapping_nodes.0 += 1;
                }
            }
        }

        // println!("Rootsplit:{}", self.nodes.debug_draw(|_| None));

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
        let mut left = Hedge(0);
        let mut extracted = self.len();
        // println!("{}", self.nodes.debug_draw(|_| None));
        // Do the same swapping as for the edge store, so that they line up
        while left < extracted {
            // println!("{left},{extracted}");
            if !subgraph.includes(&left) {
                //left is in the right place
                left.0 += 1;
            } else {
                //left needs to be swapped
                extracted.0 -= 1;
                if !subgraph.includes(&extracted) {
                    //only with an extracted that is in the wrong spot
                    // println!("{left}<->{extracted}");
                    self.swap(left, extracted);
                    left.0 += 1;
                }
            }
        }
        // println!("{}", self.nodes.debug_draw(|_| None));

        // Now need to partition the nodes into 3 ranges,
        // - 0..left_nodes do not contain any half edges in the subgraph,
        // - left_nodes..overlapping_nodes are all nodes who have incident half-edges inside the subgraph, but not all of them. These ones need to be split
        // - overlapping_nodes..len are the nodes containing only half-edges in the subgraph
        let mut left_nodes = NodeIndex(0);
        let mut extracted_nodes = self.len();

        // first we swap to the front all nodes that contain no edges from subgraph. Since we have swapped all edges in the subgraph to be after left, in the above, now we check if the neighbor iter contains no hedge>=left
        while left_nodes < extracted_nodes {
            if self.get_neighbor_iterator(left_nodes).all(|h| h < left) {
                //left is in the right place

                left_nodes.0 += 1;
            } else {
                //left needs to be swapped
                extracted_nodes.0 -= 1;
                if self
                    .get_neighbor_iterator(extracted_nodes)
                    .all(|h| h < left)
                {
                    //only with an extracted that is in the wrong spot
                    self.swap_roots(left_nodes.into(), extracted_nodes.into());
                    left_nodes.0 += 1;
                }
            }
        }

        let mut overlapping_nodes = left_nodes;
        let mut non_overlapping_extracted = self.len();

        // println!(
        //     "RootsplitnonOverlapping:{}",
        //     self.nodes.debug_draw(|_| None)
        // );

        // now we place all overlapping nodes after all nodes totally outside. We can check whether a node is overlapping if it now contains a hedge<left.
        while overlapping_nodes < non_overlapping_extracted {
            // println!("looking at possibly overlapping:{overlapping_nodes}");
            if self
                .get_neighbor_iterator(overlapping_nodes)
                .any(|h| h < left)
            {
                //overlapping is in the right place, as it intersects (is after left_nodes) but isn't fully included
                overlapping_nodes.0 += 1;
                // println!("found overlap");
            } else {
                //overlapping needs to be swapped
                non_overlapping_extracted.0 -= 1;
                if self
                    .get_neighbor_iterator(non_overlapping_extracted)
                    .any(|h| h < left)
                {
                    //only with an extracted that is in the wrong spot
                    self.swap_roots(overlapping_nodes.into(), non_overlapping_extracted.into());
                    // println!("swap overlapping");
                    overlapping_nodes.0 += 1;
                }
            }
        }

        // println!("Rootsplit:{}", self.nodes.debug_draw(|_| None));

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
            extracted_nodes.set_root(r.root_id, RootId(i))
        }
        // println!("{}", extracted_nodes.debug_draw(|_| None));
        // println!("{}", self.nodes.debug_draw(|_| None));
        #[cfg(test)]
        self.nodes.validate().unwrap();

        Forest {
            nodes: extracted_nodes,
            roots: overlapping_roots,
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
        }
    }

    fn get_neighbor_iterator(&self, node_id: NodeIndex) -> Self::NeighborsIter<'_> {
        let root_id = self.roots[node_id.0].root_id;
        ForestNeighborIter {
            // root: Some(root_id),
            iter: self.nodes.iter_preorder(root_id),
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
        let first = nodes[0];
        let x = self.roots[first.0].root_id;

        for i in nodes.iter().skip(1) {
            let y = self.roots[i.0].root_id;
            self.nodes.make_root_of(x, y);
        }

        self.roots[first.0].data = node_data_merge;

        first
    }

    fn iter_nodes_mut(
        &mut self,
    ) -> impl Iterator<Item = (NodeIndex, Self::NeighborsIter<'_>, &mut Self::NodeData)> {
        self.roots.iter_mut().enumerate().map(|(i, d)| {
            let root_id = d.root_id;
            (
                NodeIndex(i),
                ForestNeighborIter {
                    // root: Some(root_id),
                    iter: self.nodes.iter_preorder(root_id),
                },
                &mut d.data,
            )
        })
    }

    fn add_dangling_edge(
        mut self,
        source: NodeIndex,
    ) -> Result<Self, crate::half_edge::HedgeGraphError> {
        let parent = self.roots[source.0].root_id;
        self.nodes.add_dataless_child(parent);
        Ok(self)
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
                            ForestNeighborIter {
                                // root: Some(root_id),
                                iter: self.nodes.iter_preorder(root_id),
                            },
                            &mut r.data,
                        ),
                        root_id: r.root_id,
                    }
                })
                .collect(),
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
