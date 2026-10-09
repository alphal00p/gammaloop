//! The generated diagram before finalization.

use std::hash::{Hash, Hasher};

use feynkit_graph::ExternalState;
use feynkit_model::{ParticleId, VertexRuleId};
use linnet::half_edge::{
    HedgeGraph, NodeIndex,
    builder::HedgeGraphBuilder,
    involution::{EdgeIndex, Flow, Hedge, HedgePair, Orientation},
};

/// An external leg of a generated diagram with its momentum label.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ExternalNode {
    pub index: usize,
    pub state: ExternalState,
    /// The physical particle of the external state.
    pub particle: ParticleId,
}

/// An end of a [`GeneratedEdge`]: an interaction vertex or an external leg.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum Endpoint {
    Vertex(usize),
    External(ExternalNode),
}

impl Endpoint {
    pub(crate) fn external(self) -> Option<ExternalNode> {
        match self {
            Self::Vertex(_) => None,
            Self::External(external) => Some(external),
        }
    }
}

/// An edge from `source` to `sink`. A directed edge carries `particle` from its source to its
/// sink; an undirected edge carries a self-conjugate particle and keeps its endpoint order only
/// as a frame for half-edge references.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct GeneratedEdge {
    pub(crate) source: Endpoint,
    pub(crate) sink: Endpoint,
    pub(crate) particle: ParticleId,
    pub(crate) directed: bool,
}

impl GeneratedEdge {
    pub(crate) fn is_external(&self) -> bool {
        self.source.external().is_some() || self.sink.external().is_some()
    }
}

/// A generated diagram: interaction vertices labelled by their rule, internal edges between
/// them and external legs as dangling half-edges labelled by their external state.
///
/// Edge identifiers are the positions in [`Self::edges`], and vertex identifiers the positions in
/// [`Self::vertices`]. The source/sink half-edges of the underlying graph are the edge's
/// [`GeneratedEdge::source`] and [`GeneratedEdge::sink`]; a leg is the half-edge at its vertex.
#[derive(Debug, Clone)]
pub(crate) struct GeneratedGraph {
    graph: HedgeGraph<ParticleId, VertexRuleId, Option<ExternalNode>>,
}

/// Equal vertex rules and edges in identifier order.
impl PartialEq for GeneratedGraph {
    fn eq(&self, other: &Self) -> bool {
        self.vertices().eq(other.vertices()) && self.edges().eq(other.edges())
    }
}

impl Eq for GeneratedGraph {}

impl Hash for GeneratedGraph {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.n_vertices().hash(state);
        self.vertices().for_each(|rule| rule.hash(state));
        self.edges().for_each(|edge| edge.hash(state));
    }
}

impl GeneratedGraph {
    /// Build a graph from its vertex rules and its edges, in identifier order.
    ///
    /// # Panics
    /// Panics if an edge has no vertex endpoint.
    pub(crate) fn new(
        vertices: impl IntoIterator<Item = VertexRuleId>,
        edges: impl IntoIterator<Item = GeneratedEdge>,
    ) -> Self {
        let mut builder = HedgeGraphBuilder::new();
        let nodes: Vec<NodeIndex> = vertices
            .into_iter()
            .map(|rule| builder.add_node(rule))
            .collect();
        for edge in edges {
            let orientation = if edge.directed {
                Orientation::Default
            } else {
                Orientation::Undirected
            };
            match (edge.source, edge.sink) {
                (Endpoint::Vertex(source), Endpoint::Vertex(sink)) => {
                    builder.add_edge(nodes[source], nodes[sink], edge.particle, orientation)
                }
                (Endpoint::Vertex(vertex), Endpoint::External(external)) => builder
                    .add_external_edge(
                        nodes[vertex].add_data(Some(external)),
                        edge.particle,
                        orientation,
                        Flow::Source,
                    ),
                (Endpoint::External(external), Endpoint::Vertex(vertex)) => builder
                    .add_external_edge(
                        nodes[vertex].add_data(Some(external)),
                        edge.particle,
                        orientation,
                        Flow::Sink,
                    ),
                (Endpoint::External(_), Endpoint::External(_)) => {
                    panic!("an edge needs an interaction vertex")
                }
            }
        }
        Self {
            graph: builder.build(),
        }
    }

    pub(crate) fn n_vertices(&self) -> usize {
        self.graph.node_bound()
    }

    pub(crate) fn n_edges(&self) -> usize {
        self.graph.n_edges()
    }

    /// The vertex rules in vertex order.
    pub(crate) fn vertices(&self) -> impl Iterator<Item = VertexRuleId> + '_ {
        (0..self.n_vertices()).map(|vertex| self.graph[NodeIndex(vertex)])
    }

    pub(crate) fn edges(&self) -> impl Iterator<Item = GeneratedEdge> + '_ {
        (0..self.n_edges()).map(|edge| self.edge(edge))
    }

    pub(crate) fn edge(&self, edge: usize) -> GeneratedEdge {
        let (particle, pair) = &self.graph[&EdgeIndex(edge)];
        let any_hedge = match *pair {
            HedgePair::Paired { source, .. } | HedgePair::Split { source, .. } => source,
            HedgePair::Unpaired { hedge, .. } => hedge,
        };
        let directed =
            self.graph.get_edge_data_full(any_hedge).orientation != Orientation::Undirected;
        let vertex = |hedge: Hedge| Endpoint::Vertex(self.graph.node_id(hedge).0);
        let (source, sink) = match *pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                (vertex(source), vertex(sink))
            }
            HedgePair::Unpaired { hedge, flow } => {
                let external = Endpoint::External(
                    self.graph[hedge].expect("a dangling half-edge is an external leg"),
                );
                match flow {
                    Flow::Source => (vertex(hedge), external),
                    Flow::Sink => (external, vertex(hedge)),
                }
            }
        };
        GeneratedEdge {
            source,
            sink,
            particle: *particle,
            directed,
        }
    }

    /// The edge ends at `vertex` in half-edge order; a self-loop contributes both ends.
    pub(crate) fn incident(&self, vertex: usize) -> Vec<(usize, Flow)> {
        let mut hedges: Vec<Hedge> = self.graph.iter_crown(NodeIndex(vertex)).collect();
        hedges.sort();
        hedges
            .into_iter()
            .map(|hedge| (self.graph[&hedge].0, self.graph.flow(hedge)))
            .collect()
    }

    /// The external legs as `(edge, external, attachment vertex)`, in edge order.
    pub(crate) fn externals(&self) -> impl Iterator<Item = (usize, ExternalNode, usize)> + '_ {
        self.edges()
            .enumerate()
            .filter_map(|(id, edge)| match (edge.source, edge.sink) {
                (Endpoint::External(external), Endpoint::Vertex(vertex))
                | (Endpoint::Vertex(vertex), Endpoint::External(external)) => {
                    Some((id, external, vertex))
                }
                _ => None,
            })
    }

    /// Number edge endpoints as plain-graph nodes: vertices first, then one node per external
    /// leg in edge order. Algorithms written for graphs with degree-one external vertices run on
    /// this numbering unchanged.
    pub(crate) fn plain_nodes(&self) -> PlainNodes {
        let mut legs = vec![None; self.n_edges()];
        let mut externals = vec![];
        for (edge, external, _) in self.externals() {
            legs[edge] = Some(self.n_vertices() + externals.len());
            externals.push(external);
        }
        PlainNodes {
            vertices: self.n_vertices(),
            legs,
            externals,
        }
    }

    /// The graph with every external leg ending on a degree-one vertex, after the interaction
    /// vertices, keeping the edge identifiers and endpoint order. Sewing and cut enumeration act
    /// on these terminals.
    pub(crate) fn with_external_vertices(&self) -> HedgeGraph<ParticleId, PlainVertex> {
        let mut builder = HedgeGraphBuilder::new();
        let mut nodes: Vec<NodeIndex> = self
            .vertices()
            .map(|rule| builder.add_node(PlainVertex::Interaction(rule)))
            .collect();
        let mut node = |endpoint: Endpoint| match endpoint {
            Endpoint::Vertex(vertex) => nodes[vertex],
            Endpoint::External(external) => {
                nodes.push(builder.add_node(PlainVertex::External(external)));
                *nodes.last().expect("just added")
            }
        };
        let edges: Vec<_> = self
            .edges()
            .map(|edge| (node(edge.source), node(edge.sink), edge))
            .collect();
        for (source, sink, edge) in edges {
            let orientation = if edge.directed {
                Orientation::Default
            } else {
                Orientation::Undirected
            };
            builder.add_edge(source, sink, edge.particle, orientation);
        }
        builder.build()
    }
}

/// A vertex of [`GeneratedGraph::with_external_vertices`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PlainVertex {
    Interaction(VertexRuleId),
    External(ExternalNode),
}

/// See [`GeneratedGraph::plain_nodes`].
pub(crate) struct PlainNodes {
    vertices: usize,
    /// The node of each leg, by edge.
    legs: Vec<Option<usize>>,
    /// The external of each leg node, by node after the vertices.
    externals: Vec<ExternalNode>,
}

impl PlainNodes {
    pub(crate) fn len(&self) -> usize {
        self.vertices + self.externals.len()
    }

    /// The plain-graph nodes of the source and sink of `edge`.
    pub(crate) fn ends(&self, edge_id: usize, edge: &GeneratedEdge) -> (usize, usize) {
        let node = |endpoint: Endpoint| match endpoint {
            Endpoint::Vertex(vertex) => vertex,
            Endpoint::External(_) => self.legs[edge_id].expect("an external leg"),
        };
        (node(edge.source), node(edge.sink))
    }

    /// The external leg of `node`, if it is a leg node.
    pub(crate) fn external(&self, node: usize) -> Option<ExternalNode> {
        node.checked_sub(self.vertices)
            .map(|leg| self.externals[leg])
    }

    pub(crate) fn endpoint(&self, node: usize) -> Endpoint {
        self.external(node)
            .map_or(Endpoint::Vertex(node), Endpoint::External)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edges_round_trip_with_legs_on_either_end() {
        let external = |index, state| ExternalNode {
            index,
            state,
            particle: ParticleId::from_index(3),
        };
        let edges = [
            GeneratedEdge {
                source: Endpoint::External(external(0, ExternalState::Incoming)),
                sink: Endpoint::Vertex(1),
                particle: ParticleId::from_index(3),
                directed: true,
            },
            GeneratedEdge {
                source: Endpoint::Vertex(1),
                sink: Endpoint::Vertex(0),
                particle: ParticleId::from_index(4),
                directed: false,
            },
            GeneratedEdge {
                source: Endpoint::Vertex(0),
                sink: Endpoint::Vertex(0),
                particle: ParticleId::from_index(5),
                directed: true,
            },
            GeneratedEdge {
                source: Endpoint::Vertex(0),
                sink: Endpoint::External(external(1, ExternalState::Outgoing)),
                particle: ParticleId::from_index(3),
                directed: true,
            },
        ];
        let graph = GeneratedGraph::new(
            [VertexRuleId::from_index(0), VertexRuleId::from_index(2)],
            edges,
        );
        assert_eq!(graph.edges().collect::<Vec<_>>(), edges);
        assert_eq!(
            graph.vertices().collect::<Vec<_>>(),
            [VertexRuleId::from_index(0), VertexRuleId::from_index(2)]
        );
        assert_eq!(
            graph.incident(0),
            [
                (1, Flow::Sink),
                (2, Flow::Source),
                (2, Flow::Sink),
                (3, Flow::Source)
            ]
        );
        let plain = graph.plain_nodes();
        assert_eq!(plain.len(), 4);
        assert_eq!(plain.ends(0, &edges[0]), (2, 1));
        assert_eq!(plain.ends(3, &edges[3]), (0, 3));
    }
}
