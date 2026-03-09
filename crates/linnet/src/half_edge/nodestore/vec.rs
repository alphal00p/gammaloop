use std::collections::HashSet;

use crate::{
    half_edge::{
        builder::HedgeNodeBuilder,
        involution::{EdgeIndex, Hedge, HedgeVec, Involution},
        subgraph::{
            subset::SubSet, BaseSubgraph, HedgeNode, Inclusion, InternalSubGraph, ModifySubSet,
            SuBitGraph, SubSetIter, SubSetLike, SubSetOps,
        },
        swap::Swap,
        HedgeGraph, HedgeGraphError, NodeIndex, NodeVec,
    },
    tree::{
        parent_pointer::{PPNode, ParentPointerStore},
        Forest, RootData, RootId,
    },
};

use super::{NodeStorage, NodeStorageOps};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// An implementation of [`NodeStorage`] and [`NodeStorageOps`] that uses `Vec`s
/// and `BitVec`s to store node information and their incident half-edges.
///
/// This strategy is often straightforward but can have performance implications
/// for certain operations like node deletion or identification if it involves
/// frequent re-indexing or large data movements.
///
/// # Type Parameters
///
/// - `N`: The type of custom data associated with each node.
///
/// # Fields
///
/// - `node_data`: A `Vec<N>` storing the custom data for each node. The index
///   in this vector corresponds to a `NodeIndex`.
/// - `hedge_data`: A `Vec<NodeIndex>` where the index of the vector is a `Hedge`'s
///   underlying `usize` value, and the element at that index is the `NodeIndex`
///   to which the hedge is incident.
/// - `nodes`: A `Vec<BitVec>`. Each index in the outer `Vec` corresponds to a
///   `NodeIndex`. The `BitVec` at that index is a bitmask representing the set
///   of half-edges incident to that node.
pub struct NodeStorageVec<N> {
    /// Stores the custom data for each node. Indexed by `NodeIndex.0`.
    pub(crate) node_data: NodeVec<N>,
    /// Maps each half-edge index (`Hedge.0`) to the `NodeIndex` it belongs to.
    pub(crate) hedge_data: HedgeVec<NodeIndex>,
    /// For each node (indexed by `NodeIndex.0`), stores a `BitVec` representing
    /// the set of half-edges incident to it.
    pub(crate) nodes: NodeVec<SuBitGraph>, // Nodes
}

#[derive(Clone, Debug)]
/// An iterator that yields the [`Hedge`] identifiers incident to a node,
/// based on iterating over the set bits in a `BitVec`.
///
/// This is typically used by [`NodeStorageVec`] to provide an iterator
/// for its `NeighborsIter` associated type.
pub struct BitVecNeighborIter<'a> {
    /// The underlying iterator over set bits in the `BitVec`.
    iter_ones: SubSetIter<'a>,
    /// The total number of possible hedges (size of the `BitVec`), used for `ExactSizeIterator`.
    len: Hedge,
}

impl<'a> From<&'a SuBitGraph> for BitVecNeighborIter<'a> {
    fn from(value: &'a SuBitGraph) -> Self {
        Self {
            iter_ones: value.included_iter(),
            len: Hedge(value.size()),
        }
    }
}

impl<'a> From<BitVecNeighborIter<'a>> for SuBitGraph {
    fn from(value: BitVecNeighborIter<'a>) -> Self {
        let len = value.len;
        SuBitGraph::from_hedge_iter(value, len.0)
    }
}

impl<'a> From<BitVecNeighborIter<'a>> for HedgeNode {
    fn from(value: BitVecNeighborIter<'a>) -> Self {
        HedgeNode {
            internal_graph: InternalSubGraph::empty(value.len.0),
            hairs: value.into(),
        }
    }
}

impl Iterator for BitVecNeighborIter<'_> {
    type Item = Hedge;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter_ones.next()
    }
}

impl ExactSizeIterator for BitVecNeighborIter<'_> {
    fn len(&self) -> usize {
        self.iter_ones.len()
    }
}

impl<N> NodeStorage for NodeStorageVec<N> {
    type NodeData = N;
    type Neighbors = SuBitGraph;
    type NeighborsIter<'a>
        = BitVecNeighborIter<'a>
    where
        Self: 'a;
    type Storage<M> = NodeStorageVec<M>;
}

impl<N> NodeStorageVec<N> {
    // fn swap_nodes(&mut self, a: NodeIndex, b: NodeIndex) {
    //     if a != b {
    //         for i in self.nodes[a].included_iter() {
    //             self.hedge_data[i] = b;
    //         }
    //         for i in self.nodes[b].included_iter() {
    //             self.hedge_data[i] = a;
    //         }
    //         self.node_data.swap(a, b);
    //         self.nodes.swap(a, b);
    //     }
    // }

    fn from_hairs_and_data(
        node_data: impl Into<NodeVec<N>>,
        nodes: impl Into<NodeVec<SuBitGraph>>,
    ) -> Option<Self> {
        let nodes = nodes.into();
        let node_data = node_data.into();
        let n_hedges = nodes[NodeIndex(0)].size();
        let mut hedge_data: HedgeVec<_> = vec![None; n_hedges].into();

        for (i, n) in nodes.iter() {
            // println!("{:?}", n);
            for h in n.included_iter() {
                hedge_data[h] = Some(i);
            }
        }
        Some(Self {
            node_data,
            hedge_data: hedge_data
                .into_iter()
                .map(|(_, i)| i)
                .collect::<Option<HedgeVec<_>>>()?,
            nodes,
        })
    }
}

impl<N> Swap<Hedge> for NodeStorageVec<N> {
    fn is_zero_length(&self) -> bool {
        self.hedge_data.is_zero_length()
    }
    fn len(&self) -> Hedge {
        self.hedge_data.len()
    }

    fn swap(&mut self, a: Hedge, b: Hedge) {
        if a != b {
            let node_a = self.hedge_data[a];
            let node_b = self.hedge_data[b];

            self.hedge_data.swap(a, b);
            self.hedge_data[a] = node_b;
            self.hedge_data[b] = node_a;

            self.nodes[node_a].swap(a, b);
            self.nodes[node_b].swap(a, b);
        }
    }
}

impl<N> Swap<NodeIndex> for NodeStorageVec<N> {
    fn is_zero_length(&self) -> bool {
        self.nodes.is_zero_length()
    }

    fn len(&self) -> NodeIndex {
        self.nodes.len()
    }

    fn swap(&mut self, a: NodeIndex, b: NodeIndex) {
        if a != b {
            for i in self.nodes[a].included_iter() {
                self.hedge_data[i] = b;
            }
            for i in self.nodes[b].included_iter() {
                self.hedge_data[i] = a;
            }
            self.node_data.swap(a, b);
            self.nodes.swap(a, b);
        }
    }
}

impl<N> NodeStorageOps for NodeStorageVec<N> {
    type OpStorage<A> = Self::Storage<A>;
    type Base = SuBitGraph;

    fn check_nodes(&self) -> Result<(), HedgeGraphError> {
        for (h, n) in &self.hedge_data {
            if !self.nodes[*n].includes(&h) {
                Err(HedgeGraphError::NodeDoesNotContainHedge(*n, h))?;
            }
        }

        let self_n_h: Hedge = self.len();
        let mut cover = Self::Base::empty(self_n_h.0);
        for (i, node) in self.nodes.iter() {
            for h in node.included_iter() {
                if cover.includes(&h) {
                    return Err(HedgeGraphError::NodesDoNotPartition(format!(
                        "They overlap for node {i}: Cover:{cover:?}, crown: {h:?}"
                    )));
                } else {
                    cover.add(h);
                }
            }
        }

        let full = !Self::Base::empty(self_n_h.0);

        if !(cover.sym_diff(&full).is_empty()) {
            return Err(HedgeGraphError::NodesDoNotPartition(format!(
                "They do not cover the whole graph: cover {cover:?}"
            )));
        }

        Ok(())
    }

    fn delete<S: SubSetLike<Base = Self::Base>>(&mut self, subgraph: &S) {
        // println!("Deleting subgraph");
        let mut left = Hedge(0);
        let mut extracted = self.len();
        while left < extracted {
            if !subgraph.includes(&left) {
                //left is in the right place
                left.0 += 1;
            } else {
                //left needs to be swapped
                extracted.0 -= 1;
                if !subgraph.includes(&extracted) {
                    //only with an extracted that is in the wrong spot
                    self.swap(left, extracted);
                    left.0 += 1;
                }
            }
        }

        // println!("left{}", left);

        let mut left_nodes = NodeIndex(0);
        let mut extracted_nodes = self.len();
        while left_nodes < extracted_nodes {
            if !self.nodes[left_nodes].has_greater(left) {
                //left is in the right place
                left_nodes.0 += 1;
            } else {
                //left needs to be swapped
                extracted_nodes.0 -= 1;
                if !self.nodes[extracted_nodes].has_greater(left) {
                    //only with an extracted that is in the wrong spot
                    self.swap(left_nodes, extracted_nodes);
                    left_nodes.0 += 1;
                }
            }
        }

        let mut overlapping_nodes = left_nodes;
        let mut non_overlapping_extracted = self.len();

        while overlapping_nodes < non_overlapping_extracted {
            if self.nodes[overlapping_nodes].intersects(&(..left)) {
                //overlapping is in the right place, as it intersects (is after left_nodes) but isn't fully included
                overlapping_nodes.0 += 1;
            } else {
                //overlapping needs to be swapped
                non_overlapping_extracted.0 -= 1;
                if self.nodes[non_overlapping_extracted].intersects(&(..left)) {
                    //only with an extracted that is in the wrong spot
                    self.swap(overlapping_nodes, non_overlapping_extracted);
                    overlapping_nodes.0 += 1;
                }
            }
        }

        let _ = self.nodes.split_off(overlapping_nodes);
        let _ = self.node_data.split_off(overlapping_nodes);
        let _ = self.hedge_data.split_off(left);

        for i in 0..(left_nodes.0) {
            let _ = self.nodes[NodeIndex(i)].split_off(left);
            // self.nodes[i].internal_graph.filter.split_off(left.0);

            // split == 0;
        }
        for i in (left_nodes.0)..(overlapping_nodes.0) {
            let _ = self.nodes[NodeIndex(i)].split_off(left);
        }
    }

    fn extract_nodes(&mut self, nodes: impl IntoIterator<Item = NodeIndex>) -> (SuBitGraph, Self) {
        let nodes: HashSet<NodeIndex> = nodes.into_iter().collect();
        let left_nodes = <Self as Swap<NodeIndex>>::partition(self, |n| !nodes.contains(n));

        let mut extracted = SuBitGraph::empty(self.hedge_data.len().0);

        for i in left_nodes.0..self.node_data.len().0 {
            extracted.union_with(&self.nodes[NodeIndex(i)]);
        }

        let left = <Self as Swap<Hedge>>::partition(self, |h| !extracted.includes(h));

        let mut extracted_neighbors = self.nodes.split_off(left_nodes);
        for (_, s) in &mut self.nodes {
            s.split_off(left);
        }

        for (_, s) in &mut extracted_neighbors {
            *s = s.split_off(left);
        }

        let extracted_data = self.node_data.split_off(left_nodes);
        let extracted_hedges = self
            .hedge_data
            .split_off(left)
            .into_iter()
            .map(|(h, mut n)| {
                n -= left_nodes;
                // println!("Extracted hedge: {:?}{n:?}", h);
                (h, n)
            })
            .collect();

        (
            extracted,
            Self {
                nodes: extracted_neighbors,
                node_data: extracted_data,
                hedge_data: extracted_hedges,
            },
        )
    }

    fn extract<S: SubSetLike<Base = SuBitGraph>, V2>(
        &mut self,
        subgraph: &S,
        mut split_node: impl FnMut(&Self::NodeData) -> V2,
        mut owned_node: impl FnMut(Self::NodeData) -> V2,
    ) -> Self::OpStorage<V2> {
        let mut left = Hedge(0);
        let mut extracted = self.len();
        // println!("HOI");
        while left < extracted {
            if !subgraph.includes(&left) {
                //left is in the right place
                left.0 += 1;
            } else {
                //left needs to be swapped
                extracted.0 -= 1;
                if !subgraph.includes(&extracted) {
                    //only with an extracted that is in the wrong spot
                    self.swap(left, extracted);
                    left.0 += 1;
                }
            }
        }

        // println!("left{}", left);

        let mut left_nodes = NodeIndex(0);
        let mut extracted_nodes = self.len();
        while left_nodes < extracted_nodes {
            if !self.nodes[left_nodes].has_greater(left) {
                //left is in the right place
                left_nodes.0 += 1;
            } else {
                // println!(
                //     "Needs swapping left {} extracted {}",
                //     left_nodes.0, extracted_nodes.0
                // );
                //left needs to be swapped
                extracted_nodes.0 -= 1;
                // println!("{}", self.nodes[extracted_nodes].nhedges());
                // println!("{:?}", self.nodes[extracted_nodes]);
                if !self.nodes[extracted_nodes].has_greater(left) {
                    //only with an extracted that is in the wrong spot
                    self.swap(left_nodes, extracted_nodes);
                    left_nodes.0 += 1;
                }
            }
        }
        // println!("left: {}", left_nodes.0);
        // println!("extracted: {}", extracted_nodes.0);
        let mut overlapping_nodes = left_nodes;
        let mut non_overlapping_extracted = self.len();

        while overlapping_nodes < non_overlapping_extracted {
            if self.nodes[overlapping_nodes].intersects(&(..left)) {
                //overlapping is in the right place, as it intersects (is after left_nodes) but isn't fully included
                overlapping_nodes.0 += 1;
            } else {
                //overlapping needs to be swapped
                non_overlapping_extracted.0 -= 1;
                if self.nodes[non_overlapping_extracted].intersects(&(..left)) {
                    //only with an extracted that is in the wrong spot
                    self.swap(overlapping_nodes, non_overlapping_extracted);
                    overlapping_nodes.0 += 1;
                }
            }
        }

        // println!("overlapping_nodes: {}", overlapping_nodes.0);
        // println!("non_overlapping_extracted: {}", non_overlapping_extracted.0);

        let mut extracted_nodes = self.nodes.split_off(overlapping_nodes);
        let mut extracted_data: NodeVec<_> = self
            .node_data
            .split_off(overlapping_nodes)
            .into_iter()
            .map(|(_nid, n)| owned_node(n))
            .collect();

        let _ = self.hedge_data.split_off(left);

        let mut overlapping_node_hairs = NodeVec::new();
        let mut overlapping_data = NodeVec::new();

        for i in 0..(left_nodes.0) {
            let _ = self.nodes[NodeIndex(i)].split_off(left);
            // self.nodes[i].internal_graph.filter.split_off(left.0);

            // split == 0;
        }
        for i in (left_nodes.0)..(overlapping_nodes.0) {
            overlapping_data.push(split_node(&self.node_data[NodeIndex(i)]));
            // println!("og {}", self.nodes[i].nhedges());

            let overlapped = self.nodes[NodeIndex(i)].split_off(left);
            // println!("overlapped {}", overlapped.nhedges());
            overlapping_node_hairs.push(overlapped);
        }

        for (_, h) in &mut extracted_nodes {
            // println!("Init nhedges {}", h.nhedges());
            *h = h.split_off(left);

            // println!("After nhedges {}", h.nhedges());
        }

        extracted_nodes.extend(overlapping_node_hairs);
        extracted_data.extend(overlapping_data);

        NodeStorageVec::from_hairs_and_data(extracted_data, extracted_nodes)
            .expect("Extracted nodes should cover extracted hedges")
    }

    // fn add_node(&mut self, node_data: Self::NodeData) -> NodeIndex {
    //     let empty = HedgeNode::empty(self.hedge_len());
    //     self.nodes.push(empty);
    //     self.node_data.push(node_data);
    //     NodeIndex(self.nodes.len() - 1)
    // }

    fn identify_nodes(
        &mut self,
        nodes: &[NodeIndex],
        node_data_merge: Self::NodeData,
    ) -> NodeIndex {
        let n_nodes: NodeIndex = self.len();
        let n_hedgs: Hedge = self.len();
        let mut removed = SubSet::<NodeIndex>::empty(n_nodes.0);
        let mut full_node = SuBitGraph::empty(n_hedgs.0);

        for n in nodes {
            removed.add(*n);
            full_node.union_with(&self.nodes[*n]);
        }

        let replacement = removed.included_iter().next().unwrap();

        for r in removed.included_iter().skip(1).rev() {
            // let last_index = self.nodes.len() - 1;

            // Before doing anything, update any hedge pointers that point to the node being removed.
            for (_, hedge) in self.hedge_data.iter_mut() {
                if *hedge == r {
                    *hedge = replacement;
                }
            }

            // if r != last_index {
            //     // Swap the target with the last element in both vectors.
            //     self.nodes.swap(r, last_index);
            //     self.node_data.swap(r, last_index);

            //     // After swapping, update any hedge pointer that pointed to the moved element.
            //     // It used to be at last_index, now it is at r.
            //     for hedge in self.hedge_data.iter_mut() {
            //         if *hedge == NodeIndex(last_index) {
            //             *hedge = NodeIndex(r);
            //         }
            //     }
            // }
            // // Remove the (now last) element.

            // self.nodes.pop();
            // self.node_data.pop();
        }

        self.nodes[replacement] = full_node;
        self.node_data[replacement] = node_data_merge;

        replacement
    }

    fn forget_identification_history(&mut self) -> NodeVec<Self::NodeData> {
        let n_nodes: NodeIndex = self.len();
        let mut to_keep: SubSet<NodeIndex> = SubSet::empty(n_nodes.0);

        for (_, h) in &self.hedge_data {
            to_keep.add(*h);
        }

        let notk = !to_keep.clone();
        let n_hedges: Hedge = self.len();
        for n in notk.included_iter() {
            self.nodes[n] = SuBitGraph::empty(n_hedges.0);
        }

        let mut left_nodes = NodeIndex(0);
        let mut extracted_nodes = n_nodes;
        while left_nodes < extracted_nodes {
            if to_keep[left_nodes] {
                //left is in the right place
                left_nodes.0 += 1;
            } else {
                //left needs to be swapped
                extracted_nodes.0 -= 1;
                if to_keep[extracted_nodes] {
                    //only with an extracted that is in the wrong spot
                    self.swap(left_nodes, extracted_nodes);
                    // self.nodes.swap(left_nodes.0, extracted_nodes.0);
                    left_nodes.0 += 1;
                }
            }
        }

        let _ = self.nodes.split_off(left_nodes);
        self.node_data.split_off(left_nodes)
    }

    fn to_forest<U, H>(
        &self,
        map_data: impl Fn(&Self::NodeData) -> U,
    ) -> Forest<U, ParentPointerStore<H>> {
        let n_hedges: Hedge = self.len();
        let mut nodes: Vec<_> = std::iter::repeat_with(|| None).take(n_hedges.0).collect();

        let mut roots = vec![];

        for ((_, set), (_, d)) in self.nodes.iter().zip(&self.node_data) {
            let mut first = None;
            for i in set.included_iter() {
                if let Some(root) = first {
                    nodes[i.0] = Some(PPNode::dataless_child(root))
                } else {
                    first = Some(i.into());
                    nodes[i.0] = Some(PPNode::dataless_root(RootId(roots.len())));
                }
            }
            roots.push(RootData {
                root_id: first.unwrap(),
                data: map_data(d),
            });
        }
        Forest {
            nodes: nodes
                .into_iter()
                .collect::<Option<Vec<_>>>()
                .unwrap()
                .into_iter()
                .collect(),
            roots,
        }
    }

    fn iter(&self) -> impl Iterator<Item = (NodeIndex, &Self::NodeData)> {
        self.node_data.iter()
    }

    fn drain(self) -> impl Iterator<Item = (NodeIndex, Self::NodeData)> {
        self.node_data.into_iter()
    }
    fn build<I: IntoIterator<Item = HedgeNodeBuilder<N>>>(node_iter: I, n_hedges: usize) -> Self {
        let mut nodes: NodeVec<SuBitGraph> = NodeVec::new();
        let mut node_data = NodeVec::new();
        let mut hedgedata: HedgeVec<_> = vec![None; n_hedges].into();

        for (i, n) in node_iter.into_iter().enumerate() {
            for h in &n.hedges {
                hedgedata[*h] = Some(NodeIndex(i));
            }
            nodes.push(n.to_base(n_hedges));
            node_data.push(n.data);
        }

        let hedge_data = hedgedata.into_iter().map(|(_, x)| x.unwrap()).collect();

        NodeStorageVec {
            node_data,
            hedge_data,
            nodes,
        }
    }

    fn build_with_mapping<I: IntoIterator<Item = HedgeNodeBuilder<ND>>, ND>(
        node_iter: I,
        n_hedges: usize,
        mut map_data: impl FnMut(ND) -> Self::NodeData,
    ) -> Self {
        let mut nodes: NodeVec<SuBitGraph> = NodeVec::new();
        let mut node_data = NodeVec::new();
        let mut hedgedata: HedgeVec<_> = vec![None; n_hedges].into();

        for (i, n) in node_iter.into_iter().enumerate() {
            for h in &n.hedges {
                hedgedata[*h] = Some(NodeIndex(i));
            }
            nodes.push(n.to_base(n_hedges));
            node_data.push(map_data(n.data));
        }

        let hedge_data = hedgedata.into_iter().map(|(_, x)| x.unwrap()).collect();

        NodeStorageVec {
            node_data,
            hedge_data,
            nodes,
        }
    }

    fn iter_nodes(
        &self,
    ) -> impl Iterator<Item = (NodeIndex, Self::NeighborsIter<'_>, &Self::NodeData)> {
        self.nodes
            .iter()
            .map(|(_, b)| b.into())
            .zip(self.node_data.iter())
            .map(|(node, (id, data))| (id, node, data))
    }

    fn iter_nodes_mut(
        &mut self,
    ) -> impl Iterator<Item = (NodeIndex, Self::NeighborsIter<'_>, &mut Self::NodeData)> {
        self.nodes
            .iter()
            .map(|(_, b)| b.into())
            .zip(self.node_data.iter_mut())
            .map(|(node, (id, data))| (id, node, data))
    }

    fn node_id_ref(&self, hedge: Hedge) -> NodeIndex {
        self.hedge_data[hedge]
    }

    // fn get_node(&self, node_id: NodeIndex) -> Self:: {
    //     &self.nodes[node_id.0]
    // }

    fn get_neighbor_iterator(&self, node_id: NodeIndex) -> Self::NeighborsIter<'_> {
        BitVecNeighborIter {
            iter_ones: self.nodes[node_id].included_iter(),
            len: self.len(),
        }
    }

    fn get_node_data(&self, node_id: NodeIndex) -> &N {
        &self.node_data[node_id]
    }

    fn get_node_data_mut(&mut self, node_id: NodeIndex) -> &mut Self::NodeData {
        &mut self.node_data[node_id]
    }

    fn extend(self, other: Self) -> Self {
        let self_n_h: Hedge = self.len();
        let other_n_h: Hedge = other.len();
        let self_empty_filter = SuBitGraph::empty(self_n_h.0);
        let other_empty_filter = SuBitGraph::empty(other_n_h.0);
        let mut node_data = self.node_data;
        node_data.extend(other.node_data);

        let nodes: NodeVec<_> = self
            .nodes
            .into_iter()
            .map(|(_, k)| k.join(other_empty_filter.clone()))
            .chain(other.nodes.into_iter().map(|(_, k)| {
                let new_hairs = self_empty_filter.clone();
                new_hairs.join(k.clone())
            }))
            .collect();

        let mut hedge_data = self.hedge_data;
        hedge_data.extend(other.hedge_data);

        NodeStorageVec {
            node_data,
            hedge_data,
            nodes,
        }
    }

    fn extend_mut(&mut self, other: Self) {
        let self_n_h: Hedge = self.len();
        let other_n_h: Hedge = other.len();
        let self_empty_filter = SuBitGraph::empty(self_n_h.0);
        let other_empty_filter = SuBitGraph::empty(other_n_h.0);
        let node_data = &mut self.node_data;
        node_data.extend(other.node_data);

        for (_, n) in self.nodes.iter_mut() {
            n.join_mut(other_empty_filter.clone());
        }

        let nodes: NodeVec<_> = other
            .nodes
            .into_iter()
            .map(|(_, k)| {
                let new_hairs = self_empty_filter.clone();
                new_hairs.join(k.clone())
            })
            .collect();

        self.nodes.extend(nodes);

        self.hedge_data.extend(other.hedge_data);
    }
    fn add_dangling_edge(self, source: NodeIndex) -> Result<Self, HedgeGraphError> {
        if self.nodes.len() <= source {
            return Err(HedgeGraphError::NoNode);
        }
        let nodes: NodeVec<_> = self
            .nodes
            .into_iter()
            .map(|(i, mut k)| {
                if i == source {
                    k.push(true);
                } else {
                    k.push(false);
                }
                k
            })
            .collect();
        let mut hedge_data = self.hedge_data;
        hedge_data.push(source);

        Ok(NodeStorageVec {
            node_data: self.node_data,
            hedge_data,
            nodes,
        })
    }

    fn random(sources: &[Self::Neighbors], sinks: &[Self::Neighbors]) -> Self
    where
        N: Default,
    {
        let mut nodes = NodeVec::new();
        let mut node_data: NodeVec<N> = NodeVec::new();

        let mut hedge_data: HedgeVec<_> = vec![NodeIndex(0); sources[0].n_included()].into();

        for (nid, n) in sources.iter().enumerate() {
            nodes.push(n.clone());
            node_data.push(N::default());
            for i in n.included_iter() {
                hedge_data[i] = NodeIndex(nid);
            }
        }

        let len = nodes.len();

        for (nid, n) in sinks.iter().enumerate() {
            nodes.push(n.clone());
            node_data.push(N::default());

            for i in n.included_iter() {
                hedge_data[i] = NodeIndex(nid) + len;
            }
        }

        NodeStorageVec {
            node_data,
            hedge_data,
            nodes,
        }
    }

    fn check_and_set_nodes(&mut self) -> Result<(), HedgeGraphError> {
        let self_n_h: Hedge = self.len();
        let mut cover = SuBitGraph::empty(self_n_h.0);
        for (i, node) in self.nodes.iter() {
            for h in node.included_iter() {
                if cover.includes(&h) {
                    return Err(HedgeGraphError::NodesDoNotPartition(format!(
                        "They overlap. Cover:{cover:?}, crown: {h:?}"
                    )));
                } else {
                    cover.add(h);
                    self.hedge_data[h] = i;
                }
            }
        }

        let full = !SuBitGraph::empty(self_n_h.0);

        if !(cover.sym_diff(&full).is_empty()) {
            return Err(HedgeGraphError::NodesDoNotPartition(format!(
                "They do not cover the whole graph: cover {cover:?}"
            )));
        }

        Ok(())
    }

    fn map_data_ref_graph<'a, E, V2, H>(
        &'a self,
        graph: &'a HedgeGraph<E, Self::NodeData, H, Self>,
        mut node_map: impl FnMut(
            &'a HedgeGraph<E, Self::NodeData, H, Self>,
            Self::NeighborsIter<'a>,
            &'a Self::NodeData,
        ) -> V2,
    ) -> Self::OpStorage<V2> {
        let node_data = self
            .node_data
            .iter()
            .zip(self.nodes.iter())
            .map(|((_, v), (_, h))| node_map(graph, h.into(), v))
            .collect();

        NodeStorageVec {
            node_data,
            hedge_data: self.hedge_data.clone(),
            nodes: self.nodes.clone(),
        }
    }

    // fn map_data_ref_graph<'a, E, V2>(
    //     &'a self,
    //     graph: &'a HedgeGraph<E, Self::NodeData, Self>,
    //     mut node_map: impl FnMut(
    //         &'a HedgeGraph<E, Self::NodeData, Self>,
    //         &'a HedgeNode,
    //         &'a Self::NodeData,
    //     ) -> V2,
    // ) -> Self::Storage<V2> {
    //     let node_data = self
    //         .node_data
    //         .iter()
    //         .zip(self.nodes.iter())
    //         .map(|(v, h)| node_map(graph, h, v))
    //         .collect();

    //     NodeStorageVec {
    //         node_data,
    //         hedge_data: self.hedge_data.clone(),
    //         nodes: self.nodes.clone(),
    //     }
    // }

    // fn map_data_ref_mut_graph<'a, V2>(
    //     &'a mut self,
    //     mut node_map: impl FnMut(&'a HedgeNode, &'a mut Self::NodeData) -> V2,
    // ) -> Self::Storage<V2> {
    //     let node_data = self
    //         .node_data
    //         .iter_mut()
    //         .zip(self.nodes.iter())
    //         .map(|(v, h)| node_map(h, v))
    //         .collect();

    //     NodeStorageVec {
    //         node_data,
    //         hedge_data: self.hedge_data.clone(),
    //         nodes: self.nodes.clone(),
    //     }
    // }
    fn map_data_ref_mut_graph<'a, V2>(
        &'a mut self,
        mut node_map: impl FnMut(Self::NeighborsIter<'a>, &'a mut Self::NodeData) -> V2,
    ) -> Self::OpStorage<V2> {
        let node_data = self
            .node_data
            .iter_mut()
            .zip(self.nodes.iter())
            .map(|((_, v), (_, h))| node_map(h.into(), v))
            .collect();

        NodeStorageVec {
            node_data,
            hedge_data: self.hedge_data.clone(),
            nodes: self.nodes.clone(),
        }
    }

    fn map_data_ref_graph_result<'a, E, V2, H, Er>(
        &'a self,
        graph: &'a HedgeGraph<E, Self::NodeData, H, Self>,
        mut node_map: impl FnMut(
            &'a HedgeGraph<E, Self::NodeData, H, Self>,
            Self::NeighborsIter<'a>,
            &'a Self::NodeData,
        ) -> Result<V2, Er>,
    ) -> Result<Self::OpStorage<V2>, Er> {
        let node_data: Result<NodeVec<_>, Er> = self
            .node_data
            .iter()
            .zip(self.nodes.iter())
            .map(|((_, v), (_, h))| node_map(graph, h.into(), v))
            .collect();

        Ok(NodeStorageVec {
            node_data: node_data?,
            hedge_data: self.hedge_data.clone(),
            nodes: self.nodes.clone(),
        })
    }

    fn map_data_graph<'a, V2>(
        self,
        involution: &'a Involution<EdgeIndex>,
        mut f: impl FnMut(&'a Involution<EdgeIndex>, NodeIndex, Self::NodeData) -> V2,
    ) -> Self::OpStorage<V2> {
        let node_data = self
            .node_data
            .into_iter()
            .map(|(i, v)| f(involution, i, v))
            .collect();

        NodeStorageVec {
            node_data,
            hedge_data: self.hedge_data,
            nodes: self.nodes,
        }
    }

    fn map_data_graph_result<'a, V2, Err>(
        self,
        involution: &'a Involution<EdgeIndex>,
        mut f: impl FnMut(&'a Involution<EdgeIndex>, NodeIndex, Self::NodeData) -> Result<V2, Err>,
    ) -> Result<Self::OpStorage<V2>, Err> {
        let node_data = self
            .node_data
            .into_iter()
            .map(|(i, v)| f(involution, i, v))
            .collect::<Result<NodeVec<_>, Err>>()?;
        Ok(NodeStorageVec {
            node_data,
            hedge_data: self.hedge_data,
            nodes: self.nodes,
        })
    }

    fn new_nodevec<'a, V2>(
        &'a self,
        mut node_map: impl FnMut(NodeIndex, Self::NeighborsIter<'a>, &'a Self::NodeData) -> V2,
    ) -> NodeVec<V2> {
        self.node_data
            .iter()
            .zip(self.nodes.iter())
            .map(|((i, v), (_, h))| node_map(i, h.into(), v))
            .collect()
    }
}
