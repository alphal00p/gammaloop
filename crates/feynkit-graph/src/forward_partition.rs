//! Cut partitions of unsewn forward graphs with whole vertex sides.
//!
//! A forward graph is unsewn when each external leg ends at its own
//! degree-one vertex: the incoming legs belong on the left of every cut and
//! the outgoing legs on the right. Linnet's cut search returns each side as the
//! half-edges it reached, together with the far half of each bridge it
//! skipped, so a vertex behind such a bridge can appear split between the
//! sides. [`ForwardGraph::partition`] rebuilds whole vertex sides from a cut
//! boundary, and [`ForwardGraph::leg_seeded_cuts`] enumerates cuts as main's
//! generator does, with only the leg vertices fixed.
//!
//! ```
//! use feynkit_graph::{ForwardGraph, LegSeededCut};
//! use linnet::half_edge::HedgeGraph;
//! use linnet::half_edge::builder::HedgeGraphBuilder;
//! use linnet::half_edge::involution::Orientation;
//! use linnet::half_edge::subgraph::SubSetLike;
//!
//! // in -> a, a = b (two lines), b -> out.
//! let mut builder = HedgeGraphBuilder::new();
//! let [incoming, a, b, outgoing] = [(); 4].map(|()| builder.add_node(()));
//! builder.add_edge(incoming, a, (), Orientation::Default);
//! builder.add_edge(a, b, (), Orientation::Default);
//! builder.add_edge(a, b, (), Orientation::Default);
//! builder.add_edge(b, outgoing, (), Orientation::Default);
//! let graph: HedgeGraph<(), ()> = builder.build();
//!
//! let (incoming, outgoing) = ([incoming], [outgoing]);
//! let forward = ForwardGraph::new(&graph, &incoming, &outgoing).unwrap();
//! let partitions: Vec<_> = forward
//!     .leg_seeded_cuts(|_| false)
//!     .into_iter()
//!     .filter_map(|cut| match cut {
//!         LegSeededCut::Partition(partition) => Some(partition),
//!         LegSeededCut::CrossesLeg { .. } => None,
//!     })
//!     .collect();
//! // Only the two parallel lines separate the legs without crossing one.
//! assert_eq!(partitions.len(), 1);
//! assert_eq!(partitions[0].cut().n_included(), 2);
//! ```

use std::collections::BTreeSet;

use linnet::half_edge::HedgeGraph;
use linnet::half_edge::NodeIndex;
use linnet::half_edge::involution::Hedge;
use linnet::half_edge::nodestore::NodeStorageOps;
use linnet::half_edge::subgraph::{
    Inclusion, ModifySubSet, OrientedCut, SuBitGraph, SubGraphOps, SubSetLike,
};
use thiserror::Error;

/// An unsewn forward graph with its incoming and outgoing leg vertices.
pub struct ForwardGraph<'a, E, V, H, N: NodeStorageOps<NodeData = V>> {
    graph: &'a HedgeGraph<E, V, H, N>,
    incoming: &'a [NodeIndex],
    outgoing: &'a [NodeIndex],
    legs: BTreeSet<NodeIndex>,
}

/// Why leg vertices do not describe an unsewn forward graph.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum ForwardGraphError {
    #[error("a forward graph needs at least one incoming and one outgoing leg")]
    MissingLeg,
    #[error("leg vertex {0:?} is not a vertex of the graph")]
    UnknownVertex(NodeIndex),
    #[error("leg vertex {0:?} is not a vertex of degree one")]
    NotDegreeOne(NodeIndex),
    #[error("leg vertex {0:?} is listed more than once")]
    RepeatedLeg(NodeIndex),
}

/// A partition of an unsewn forward graph whose sides consist of whole vertex
/// crowns, with every leg on the side of its attachment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ForwardPartition {
    left: SuBitGraph,
    right: SuBitGraph,
    cut: SuBitGraph,
}

/// One cut found with only the leg vertices fixed.
pub enum LegSeededCut {
    /// A cut that crosses no leg, as a partition with whole vertex sides.
    Partition(ForwardPartition),
    /// A cut that crosses a leg, as the cut search returned it.
    CrossesLeg {
        left: SuBitGraph,
        cut: OrientedCut,
        right: SuBitGraph,
    },
}

impl ForwardPartition {
    /// The half-edges of the vertices on the incoming side.
    pub fn left(&self) -> &SuBitGraph {
        &self.left
    }

    /// The half-edges of the vertices on the outgoing side.
    pub fn right(&self) -> &SuBitGraph {
        &self.right
    }

    /// The incoming-side half of every edge between the sides.
    pub fn cut(&self) -> &SuBitGraph {
        &self.cut
    }
}

impl<'a, E, V, H, N: NodeStorageOps<NodeData = V>> ForwardGraph<'a, E, V, H, N> {
    /// Describe `graph` as a forward graph whose `incoming` legs belong on the
    /// left and `outgoing` legs on the right. Both lists must be non-empty and
    /// hold distinct vertices of degree one.
    ///
    /// # Errors
    ///
    /// Returns [`ForwardGraphError::MissingLeg`] when either list is empty,
    /// [`ForwardGraphError::UnknownVertex`] for an index outside the graph,
    /// [`ForwardGraphError::NotDegreeOne`] for a vertex with more or fewer
    /// than one half-edge, and [`ForwardGraphError::RepeatedLeg`] for a vertex
    /// listed twice.
    pub fn new(
        graph: &'a HedgeGraph<E, V, H, N>,
        incoming: &'a [NodeIndex],
        outgoing: &'a [NodeIndex],
    ) -> Result<Self, ForwardGraphError> {
        if incoming.is_empty() || outgoing.is_empty() {
            return Err(ForwardGraphError::MissingLeg);
        }
        let mut legs = BTreeSet::new();
        for &leg in incoming.iter().chain(outgoing) {
            if leg.0 >= graph.n_nodes() {
                return Err(ForwardGraphError::UnknownVertex(leg));
            }
            if graph.iter_crown(leg).count() != 1 {
                return Err(ForwardGraphError::NotDegreeOne(leg));
            }
            if !legs.insert(leg) {
                return Err(ForwardGraphError::RepeatedLeg(leg));
            }
        }
        Ok(Self {
            graph,
            incoming,
            outgoing,
            legs,
        })
    }

    /// Whether `boundary`, holding one half of each crossing edge, crosses a
    /// leg.
    pub fn crosses_leg(&self, boundary: &SuBitGraph) -> bool {
        boundary.included_iter().any(|hedge| {
            self.legs.contains(&self.graph.node_id(hedge))
                || self
                    .legs
                    .contains(&self.graph.node_id(self.graph.inv(hedge)))
        })
    }

    /// The partition whose left side holds every vertex reached from the
    /// incoming legs without crossing `boundary`, given as one half of each
    /// crossing edge. Returns `None` when the boundary crosses a leg, and
    /// unless every crossing edge joins the two sides and every outgoing leg
    /// is on the right.
    pub fn partition(&self, boundary: &SuBitGraph) -> Option<ForwardPartition> {
        if self.crosses_leg(boundary) {
            return None;
        }
        let graph = self.graph;
        let crosses =
            |hedge: Hedge| boundary.includes(&hedge) || boundary.includes(&graph.inv(hedge));
        let mut vertices: BTreeSet<NodeIndex> = self.incoming.iter().copied().collect();
        let mut pending = self.incoming.to_vec();
        while let Some(vertex) = pending.pop() {
            for hedge in graph.iter_crown(vertex) {
                let partner = graph.inv(hedge);
                if partner != hedge && !crosses(hedge) {
                    let next = graph.node_id(partner);
                    if vertices.insert(next) {
                        pending.push(next);
                    }
                }
            }
        }
        if self.outgoing.iter().any(|leg| vertices.contains(leg)) {
            return None;
        }
        let mut left: SuBitGraph = graph.empty_subgraph();
        for &vertex in &vertices {
            for hedge in graph.iter_crown(vertex) {
                left.add(hedge);
            }
        }
        let mut cut: SuBitGraph = graph.empty_subgraph();
        for hedge in boundary.included_iter() {
            let (near, far) = (hedge, graph.inv(hedge));
            match (left.includes(&near), left.includes(&far)) {
                (true, false) => cut.add(near),
                (false, true) => cut.add(far),
                _ => return None,
            }
        }
        let right = left.complement(graph);
        Some(ForwardPartition { left, right, cut })
    }

    /// Enumerate the cuts with only the leg vertices fixed. A cut that crosses
    /// a leg is returned as the search found it. Any other cut is rebuilt by
    /// [`Self::partition`] unless `known` accepts its boundary, so callers can
    /// skip partitions they already hold; each partition is returned once and
    /// invalid boundaries are dropped.
    pub fn leg_seeded_cuts(&self, known: impl Fn(&SuBitGraph) -> bool) -> Vec<LegSeededCut> {
        let mut seen = BTreeSet::new();
        let mut found = Vec::new();
        for (left, cut, right) in self.graph.all_cuts_from_ids(self.incoming, self.outgoing) {
            if self.crosses_leg(&cut.left) {
                found.push(LegSeededCut::CrossesLeg { left, cut, right });
            } else if !known(&cut.left)
                && !seen.contains(&cut.left)
                && let Some(partition) = self.partition(&cut.left)
            {
                seen.insert(cut.left);
                found.push(LegSeededCut::Partition(partition));
            }
        }
        found
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use linnet::half_edge::builder::HedgeGraphBuilder;
    use linnet::half_edge::involution::Orientation;

    /// in -> a -> out, whose only edges are the two legs.
    fn chain() -> (HedgeGraph<(), ()>, [NodeIndex; 3]) {
        let mut builder = HedgeGraphBuilder::new();
        let nodes = [(); 3].map(|()| builder.add_node(()));
        builder.add_edge(nodes[0], nodes[1], (), Orientation::Default);
        builder.add_edge(nodes[1], nodes[2], (), Orientation::Default);
        (builder.build(), nodes)
    }

    #[test]
    fn leg_vertices_are_checked_before_use() {
        let (graph, [incoming, middle, outgoing]) = chain();
        let outside = NodeIndex(graph.n_nodes());
        for (incoming, outgoing, expected) in [
            (vec![], vec![outgoing], ForwardGraphError::MissingLeg),
            (vec![incoming], vec![], ForwardGraphError::MissingLeg),
            (
                vec![incoming],
                vec![outside],
                ForwardGraphError::UnknownVertex(outside),
            ),
            (
                vec![incoming],
                vec![middle],
                ForwardGraphError::NotDegreeOne(middle),
            ),
            (
                vec![incoming],
                vec![incoming],
                ForwardGraphError::RepeatedLeg(incoming),
            ),
        ] {
            assert_eq!(
                ForwardGraph::new(&graph, &incoming, &outgoing).err(),
                Some(expected)
            );
        }
    }

    #[test]
    fn partitions_never_separate_a_leg_from_its_attachment() {
        let (graph, [incoming, _, outgoing]) = chain();
        let (incoming, outgoing) = ([incoming], [outgoing]);
        let forward = ForwardGraph::new(&graph, &incoming, &outgoing).unwrap();
        // The boundary holding only the incoming leg puts `in` alone on the left.
        let mut boundary: SuBitGraph = graph.empty_subgraph();
        boundary.add(graph.iter_crown(incoming[0]).next().unwrap());
        assert!(forward.crosses_leg(&boundary));
        assert_eq!(forward.partition(&boundary), None);
        // Every cut of this graph crosses a leg.
        assert!(
            forward
                .leg_seeded_cuts(|_| false)
                .iter()
                .all(|cut| matches!(cut, LegSeededCut::CrossesLeg { .. }))
        );
    }
}
