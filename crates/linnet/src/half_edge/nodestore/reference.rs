//! A naive [`Topology`] that states the contract as plainly as possible: each
//! half-edge stores its live node and crowns are found by scanning. It orders
//! partitions and compactions by the tail-fill rule stated directly rather
//! than by the shared kernel. The differential tests check
//! [`super::NodeForest`] against it.

use std::ops::Range;

use super::{NodeFate, Partition, Topology};
use crate::half_edge::{
    involution::Hedge,
    remap::{positions, HedgeRemap, NodeRemap},
    HedgeGraphError, NodeIndex,
};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct RefTopology {
    /// Each half-edge's live node.
    label: Vec<usize>,
    anchor: Vec<Option<usize>>,
    /// The node each absorbed node was merged into.
    absorbed_into: Vec<Option<usize>>,
}

/// The tail-fill order stated directly, independently of the shared kernel:
/// the positions `start..` of `order` accepted by `keep` come first, the
/// rejected ones below the split, ascending, trading places with the accepted
/// ones above it, descending. Returns the split.
fn tail_fill(order: &mut [u32], start: usize, keep: impl Fn(u32) -> bool) -> usize {
    let split = start + order[start..].iter().filter(|&&old| keep(old)).count();
    let holes: Vec<usize> = (start..split).filter(|&at| !keep(order[at])).collect();
    let fillers: Vec<usize> = (split..order.len())
        .rev()
        .filter(|&at| keep(order[at]))
        .collect();
    for (hole, filler) in holes.into_iter().zip(fillers) {
        order.swap(hole, filler);
    }
    split
}

impl RefTopology {
    fn members(&self, node: usize) -> Vec<Hedge> {
        (0..self.label.len())
            .filter(|&hedge| self.label[hedge] == node)
            .map(Hedge)
            .collect()
    }

    /// The topology of the half-edges and nodes at positions `hedges` and
    /// `nodes` of the orders `order` and `node_order`, whose inverses are `at`
    /// and `new_id`.
    fn side(
        &self,
        (order, at): (&[u32], &[u32]),
        (node_order, new_id): (&[u32], &[u32]),
        hedges: Range<usize>,
        nodes: Range<usize>,
    ) -> Self {
        let label: Vec<usize> = (order[hedges.clone()].iter())
            .map(|&old| new_id[self.label[old as usize]] as usize - nodes.start)
            .collect();
        let anchor = (node_order[nodes.clone()].iter().enumerate())
            .map(|(new, &old)| {
                let anchor = at[self.anchor[old as usize]?] as usize;
                match hedges.contains(&anchor) {
                    true => Some(anchor - hedges.start),
                    false => label.iter().position(|&node| node == new),
                }
            })
            .collect();
        RefTopology {
            label,
            anchor,
            absorbed_into: vec![None; nodes.len()],
        }
    }
}

impl Topology for RefTopology {
    type Crown<'a> = std::vec::IntoIter<Hedge>;

    fn n_hedges(&self) -> usize {
        self.label.len()
    }

    fn node_bound(&self) -> usize {
        self.anchor.len()
    }

    fn n_live(&self) -> usize {
        self.absorbed_into
            .iter()
            .filter(|into| into.is_none())
            .count()
    }

    fn node_of(&self, hedge: Hedge) -> NodeIndex {
        NodeIndex(self.label[hedge.0])
    }

    fn crown(&self, node: NodeIndex) -> Self::Crown<'_> {
        self.members(node.0).into_iter()
    }

    fn is_live(&self, node: NodeIndex) -> bool {
        self.absorbed_into[node.0].is_none()
    }

    fn representative(&self, mut node: NodeIndex) -> NodeIndex {
        while let Some(into) = self.absorbed_into[node.0] {
            node = NodeIndex(into);
        }
        node
    }

    fn anchor(&self, node: NodeIndex) -> Option<Hedge> {
        self.anchor[node.0].map(Hedge)
    }

    fn set_anchor(&mut self, node: NodeIndex, anchor: Option<Hedge>) {
        assert!(
            anchor.is_none_or(|hedge| self.label[hedge.0] == node.0),
            "not in the crown"
        );
        self.anchor[node.0] = anchor.map(|hedge| hedge.0);
    }

    fn from_labels(labels: impl ExactSizeIterator<Item = NodeIndex>, node_bound: usize) -> Self {
        let label: Vec<usize> = labels.map(|node| node.0).collect();
        assert!(label.iter().all(|&node| node < node_bound));
        RefTopology {
            label,
            anchor: vec![None; node_bound],
            absorbed_into: vec![None; node_bound],
        }
    }

    fn push_node(&mut self) -> NodeIndex {
        self.anchor.push(None);
        self.absorbed_into.push(None);
        NodeIndex(self.anchor.len() - 1)
    }

    fn attach(&mut self, node: NodeIndex) -> Result<Hedge, HedgeGraphError> {
        match self.absorbed_into.get(node.0) {
            None => Err(HedgeGraphError::NoNode),
            Some(Some(_)) => Err(HedgeGraphError::InvalidNode(node)),
            Some(None) => {
                self.label.push(node.0);
                Ok(Hedge(self.label.len() - 1))
            }
        }
    }

    fn extend(&mut self, other: Self) -> Result<(), HedgeGraphError> {
        let (hedges, nodes) = (self.label.len(), self.anchor.len());
        self.label
            .extend(other.label.iter().map(|node| node + nodes));
        self.anchor
            .extend(other.anchor.iter().map(|a| a.map(|a| a + hedges)));
        self.absorbed_into
            .extend(other.absorbed_into.iter().map(|n| n.map(|n| n + nodes)));
        Ok(())
    }

    fn identify(&mut self, nodes: &[NodeIndex]) -> NodeIndex {
        let size = |store: &Self, node: NodeIndex| store.members(node.0).len();
        let mut survivor = self.representative(nodes[0]);
        for &node in &nodes[1..] {
            let node = self.representative(node);
            if size(self, node) > size(self, survivor) {
                survivor = node;
            }
        }
        for &node in nodes {
            let node = self.representative(node).0;
            if node != survivor.0 {
                for label in self.label.iter_mut().filter(|label| **label == node) {
                    *label = survivor.0;
                }
                self.absorbed_into[node] = Some(survivor.0);
                let anchor = self.anchor[node].take();
                self.anchor[survivor.0] = self.anchor[survivor.0].or(anchor);
            }
        }
        survivor
    }

    fn compact(&mut self) -> Result<NodeRemap, HedgeGraphError> {
        let n = self.anchor.len();
        let mut order: Vec<u32> = (0..n as u32).collect();
        let live = tail_fill(&mut order, 0, |node| self.is_live(NodeIndex(node as usize)));
        let new_id = positions(&order);
        let remap = NodeRemap::from_entries(
            (0..n).map(|node| {
                let root = self.representative(NodeIndex(node)).0;
                Some((new_id[root] as usize, root == node))
            }),
            live,
        );
        for label in &mut self.label {
            *label = new_id[*label] as usize;
        }
        self.anchor = order[..live]
            .iter()
            .map(|&old| self.anchor[old as usize])
            .collect();
        self.absorbed_into = vec![None; live];
        Ok(remap)
    }

    fn partition_hedges(
        &mut self,
        moved: impl Fn(Hedge) -> bool,
        fate: impl Fn(NodeIndex, usize, usize) -> NodeFate,
        build_moved: bool,
    ) -> Result<Partition<Self>, HedgeGraphError> {
        if self.n_live() != self.node_bound() {
            return Err(HedgeGraphError::PendingIdentifications);
        }
        let (n_hedges, n_nodes) = (self.label.len(), self.anchor.len());
        let mut order: Vec<u32> = (0..n_hedges as u32).collect();
        let n_staying = tail_fill(&mut order, 0, |hedge| !moved(Hedge(hedge as usize)));
        let fates: Vec<NodeFate> = (0..n_nodes)
            .map(|node| {
                let members = self.members(node);
                let n_moved = members.iter().filter(|&&hedge| moved(hedge)).count();
                fate(NodeIndex(node), members.len() - n_moved, n_moved)
            })
            .collect();
        let mut node_order: Vec<u32> = (0..n_nodes as u32).collect();
        let split = tail_fill(&mut node_order, 0, |node| {
            fates[node as usize] == NodeFate::Kept
        });
        let inside = tail_fill(&mut node_order, split, |node| {
            fates[node as usize] == NodeFate::Split
        });
        let (at, new_id) = (positions(&order), positions(&node_order));
        let side = |hedges, nodes| self.side((&order, &at), (&node_order, &new_id), hedges, nodes);
        let moved_side = build_moved.then(|| side(n_staying..n_hedges, split..n_nodes));
        let kept = side(0..n_staying, 0..inside);
        *self = kept;
        Ok(Partition {
            moved: moved_side,
            hedges: HedgeRemap::bijection(at),
            n_staying,
            nodes: NodeRemap::bijection(new_id),
            split,
            inside,
        })
    }

    fn permute(&mut self, hedges: &HedgeRemap, nodes: &NodeRemap) -> Result<(), HedgeGraphError> {
        if self.n_live() != self.node_bound() {
            return Err(HedgeGraphError::PendingIdentifications);
        }
        let new_hedge = |hedge: usize| hedges.get(Hedge(hedge)).unwrap().0;
        let new_node = |node: usize| nodes.get(NodeIndex(node)).unwrap().0;
        let mut label = vec![0; self.label.len()];
        for (hedge, &node) in self.label.iter().enumerate() {
            label[new_hedge(hedge)] = new_node(node);
        }
        let mut anchor = vec![None; self.anchor.len()];
        for (node, &old) in self.anchor.iter().enumerate() {
            anchor[new_node(node)] = old.map(new_hedge);
        }
        (self.label, self.anchor) = (label, anchor);
        Ok(())
    }
}
