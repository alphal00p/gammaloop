use std::ops::{Range, RangeFrom, RangeInclusive, RangeTo, RangeToInclusive};

use crate::half_edge::builder::HedgeNodeBuilder;
use crate::half_edge::involution::HedgePair;
use crate::half_edge::subgraph::{SuBitGraph, SubGraphLike, SubGraphOps};
use crate::half_edge::{Hedge, HedgeGraph, NodeStorageOps};

use super::contracted::ContractedSubGraph;
use super::{internal::InternalSubGraph, SubSetLike, SubSetOps};
use super::{Inclusion, ModifySubSet, SubSetIter};

#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Represents a "node" in a graph that might itself have internal structure,
/// often used in graph contraction or hierarchical graph representations.
///
/// A `HedgeNode` consists of two main parts:
/// 1.  An [`InternalSubGraph`]: A set of edges that are fully contained within this node.
/// 2.  A set of "hairs": These are half-edges that connect the `internal_graph` (or
///     other hairs) to the rest of the larger graph, acting as the connection points
///     for this `HedgeNode`.
///
/// This structure allows a complex part of a graph to be treated as a single logical node
/// while retaining information about its internal topology and external connections.
pub struct HedgeNode {
    /// The subgraph of edges fully contained within this node.
    /// This internal graph itself has no dangling edges to the outside;
    /// all such connections are represented by `hairs`.
    pub internal_graph: InternalSubGraph,
    /// A bitmask representing the set of half-edges that are "hairs" of this node.
    /// Hairs are external connections: one end is implicitly connected to the
    /// `internal_graph` or other hairs of this node, and the other end connects
    /// to other parts of the main graph.
    pub hairs: SuBitGraph,
}

impl Inclusion<Hedge> for HedgeNode {
    fn includes(&self, hedge_id: &Hedge) -> bool {
        self.internal_graph.includes(hedge_id) || self.hairs.includes(hedge_id)
    }

    fn intersects(&self, other: &Hedge) -> bool {
        self.includes(other)
    }
}

impl Inclusion<HedgePair> for HedgeNode {
    fn includes(&self, hedge_id: &HedgePair) -> bool {
        self.internal_graph.includes(hedge_id) || self.hairs.includes(hedge_id)
    }

    fn intersects(&self, other: &HedgePair) -> bool {
        self.includes(other)
    }
}

impl Inclusion<HedgeNode> for HedgeNode {
    fn includes(&self, other: &HedgeNode) -> bool {
        self.internal_graph.includes(&other.internal_graph)
    }

    fn intersects(&self, other: &HedgeNode) -> bool {
        self.hairs.intersects(&other.hairs)
    }
}

impl Inclusion<SuBitGraph> for HedgeNode {
    fn includes(&self, other: &SuBitGraph) -> bool {
        self.internal_graph.includes(other) || self.hairs.includes(other)
    }

    fn intersects(&self, other: &SuBitGraph) -> bool {
        self.hairs.intersects(other)
    }
}

impl SubGraphLike for HedgeNode {
    fn nedges<E, V, H, N: NodeStorageOps<NodeData = V>>(
        &self,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> usize {
        self.internal_graph.nedges(graph)
    }

    fn hairs(&self, node: impl Iterator<Item = Hedge>) -> SuBitGraph {
        let mut hairs = SuBitGraph::empty(self.size());
        for h in node {
            if self.includes(&h) {
                hairs.add(h)
            }
        }
        hairs
    }
}

impl SubSetLike<Hedge> for HedgeNode {
    type Base = SuBitGraph;
    type BaseIter<'a> = SubSetIter<'a, Hedge>;

    fn included_iter(&self) -> Self::BaseIter<'_> {
        self.hairs.included_iter()
    }

    fn has_greater(&self, hedge: Hedge) -> bool {
        self.hairs.has_greater(hedge)
    }

    fn size(&self) -> usize {
        self.hairs.size()
    }

    fn join_mut(&mut self, other: Self) {
        self.hairs.join_mut(other.hairs);
        self.internal_graph.join_mut(other.internal_graph);
    }

    fn n_included(&self) -> usize {
        self.hairs.n_included()
    }

    fn included(&self) -> &SuBitGraph {
        self.hairs.included()
    }

    fn string_label(&self) -> String {
        (self.hairs.string_label() + "⦻") + self.internal_graph.string_label().as_str()
    }

    fn is_empty(&self) -> bool {
        self.hairs.is_empty() && self.internal_graph.is_empty()
    }

    fn empty(size: usize) -> Self {
        Self {
            internal_graph: InternalSubGraph::empty(size),
            hairs: SuBitGraph::empty(size),
        }
    }

    fn from_base62(label: &str, size: usize) -> Option<Self> {
        let hairs = SuBitGraph::from_base62(label, size)?;
        Some(HedgeNode {
            internal_graph: InternalSubGraph::empty(size),
            hairs,
        })
    }
}
impl Inclusion<Range<Hedge>> for HedgeNode {
    fn includes(&self, other: &Range<Hedge>) -> bool {
        (other.start.0..other.end.0).all(|a| self.includes(&Hedge(a)))
    }

    fn intersects(&self, other: &Range<Hedge>) -> bool {
        (other.start.0..other.end.0).any(|a| self.includes(&Hedge(a)))
    }
}

impl Inclusion<RangeTo<Hedge>> for HedgeNode {
    fn includes(&self, other: &RangeTo<Hedge>) -> bool {
        (0..other.end.0).all(|a| self.includes(&Hedge(a)))
    }

    fn intersects(&self, other: &RangeTo<Hedge>) -> bool {
        (0..other.end.0).any(|a| self.includes(&Hedge(a)))
    }
}

impl Inclusion<RangeToInclusive<Hedge>> for HedgeNode {
    fn includes(&self, other: &RangeToInclusive<Hedge>) -> bool {
        (0..=other.end.0).all(|a| self.includes(&Hedge(a)))
    }

    fn intersects(&self, other: &RangeToInclusive<Hedge>) -> bool {
        (0..=other.end.0).any(|a| self.includes(&Hedge(a)))
    }
}

impl Inclusion<RangeFrom<Hedge>> for HedgeNode {
    fn includes(&self, other: &RangeFrom<Hedge>) -> bool {
        (other.start.0..).all(|a| self.includes(&Hedge(a)))
    }

    fn intersects(&self, other: &RangeFrom<Hedge>) -> bool {
        (other.start.0..).any(|a| self.includes(&Hedge(a)))
    }
}

impl Inclusion<RangeInclusive<Hedge>> for HedgeNode {
    fn includes(&self, other: &RangeInclusive<Hedge>) -> bool {
        (other.start().0..=other.end().0).all(|a| self.includes(&Hedge(a)))
    }

    fn intersects(&self, other: &RangeInclusive<Hedge>) -> bool {
        (other.start().0..=other.end().0).any(|a| self.includes(&Hedge(a)))
    }
}
impl SubSetOps<Hedge> for HedgeNode {
    fn union_with_iter(&mut self, other: impl Iterator<Item = Hedge>) {
        self.hairs.union_with_iter(other);
    }

    fn union_with(&mut self, other: &Self) {
        // union is the intersection of the internal graphs, and the union of the external graph.
        self.internal_graph.intersect_with(&other.internal_graph);
        self.hairs.union_with(&other.hairs);
    }

    fn intersect_with(&mut self, other: &Self) {
        // intersection is the union of the internal graphs, and the intersection of the external graph.
        self.internal_graph.union_with(&other.internal_graph);
        self.hairs.intersect_with(&other.hairs);
    }

    fn sym_diff_with(&mut self, other: &Self) {
        // external hedges that are only present in one of the two graphs.
        // contracted parts unioned
        self.internal_graph.union_with(&other.internal_graph);
        self.hairs.sym_diff_with(&other.hairs);
    }

    fn empty_intersection(&self, other: &Self) -> bool {
        self.hairs.empty_intersection(&other.hairs)
    }

    fn empty_union(&self, other: &Self) -> bool {
        self.hairs.empty_union(&other.hairs)
    }

    fn subtract_with(&mut self, other: &Self) {
        self.internal_graph.union_with(&other.internal_graph);
        self.hairs.subtract_with(&other.hairs);
    }
}

impl SubGraphOps for HedgeNode {
    fn complement<E, V, H, N: NodeStorageOps<NodeData = V>>(
        &self,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> Self {
        Self::from_internal_graph(self.internal_graph.complement(graph), graph)
    }
}
impl HedgeNode {
    pub fn from_internal_graph<E, V, H, N: NodeStorageOps<NodeData = V>>(
        subgraph: InternalSubGraph,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> Self {
        graph.nesting_node_from_subgraph(subgraph)
    }

    pub fn weakly_disjoint(&self, other: &HedgeNode) -> bool {
        let internals = self
            .internal_graph
            .filter
            .intersection(&other.internal_graph.filter);

        internals.is_empty()
    }

    pub fn strongly_disjoint(&self, other: &HedgeNode) -> bool {
        let internals = self
            .internal_graph
            .filter
            .intersection(&other.internal_graph.filter);

        let externals_in_self = self.internal_graph.filter.intersection(&other.hairs);
        let externals_in_other = self.hairs.intersection(&other.internal_graph.filter);

        internals.is_empty() && externals_in_self.is_empty() && externals_in_other.is_empty()
    }

    pub fn internal_graph_union(&self, other: &HedgeNode) -> InternalSubGraph {
        InternalSubGraph {
            filter: self
                .internal_graph
                .filter
                .union(&other.internal_graph.filter),
            loopcount: None,
        }
    }

    pub fn from_builder<V>(builder: &HedgeNodeBuilder<V>, len: usize) -> Self {
        let internal_graph = InternalSubGraph::empty(len);
        let mut externalhedges = SuBitGraph::empty(len);

        for hedge in &builder.hedges {
            externalhedges.add(*hedge);
        }

        HedgeNode {
            internal_graph,
            hairs: externalhedges,
        }
    }

    pub fn is_node(&self) -> bool {
        self.internal_graph.is_empty()
    }

    /// Fixes the node by ensuring that all hairs are true hairs and not internal, and if they are move them to the internal graph. ! Also moves dangling edges to the internal graph.// not sure if this is kosher
    pub fn fix<E, V, H, N: NodeStorageOps<NodeData = V>>(
        &mut self,
        graph: &HedgeGraph<E, V, H, N>,
    ) {
        for i in self.hairs.included_iter() {
            let invh = graph.inv(i);
            if self.hairs.includes(&invh) {
                self.internal_graph.filter.add(i);
                self.internal_graph.filter.add(invh);
            }
        }
        self.hairs.subtract_with(&self.internal_graph.filter);
    }

    /// adds all hairs possible
    pub fn add_all_hairs<E, V, H, N: NodeStorageOps<NodeData = V>>(
        &mut self,
        graph: &HedgeGraph<E, V, H, N>,
    ) {
        let mut hairs: SuBitGraph = graph.empty_subgraph();

        for i in self.internal_graph.included_iter() {
            hairs.union_with_iter(graph.neighbors(i));
        }

        for i in self.hairs.included_iter() {
            hairs.union_with_iter(graph.neighbors(i));
        }

        let mut nh = !hairs;
        nh.union_with(&self.internal_graph.filter);
        self.hairs = !nh;
    }

    pub fn valid<E, V, H, N: NodeStorageOps<NodeData = V>>(
        &self,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> bool {
        for i in self.hairs.included_iter() {
            let invh = graph.inv(i);
            if self.hairs.includes(&invh) {
                return false;
            }
        }
        true
    }

    pub fn internal_and_hairs(&self) -> SuBitGraph {
        self.internal_graph.filter.union(&self.hairs)
    }

    pub fn is_subgraph(&self) -> bool {
        !self.is_node()
    }
}

impl From<HedgeNode> for ContractedSubGraph {
    fn from(value: HedgeNode) -> Self {
        ContractedSubGraph {
            internal_graph: value.internal_graph,
            allhedges: value.hairs,
        }
    }
}

impl From<SuBitGraph> for HedgeNode {
    fn from(value: SuBitGraph) -> Self {
        HedgeNode {
            internal_graph: InternalSubGraph::empty(value.size()),
            hairs: value,
        }
    }
}
