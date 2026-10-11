//! The generated diagram before finalization.

use std::{
    hash::{Hash, Hasher},
    sync::Arc,
};

use feynkit_graph::ExternalState;
use feynkit_model::{ParticleId, VertexRuleId};
use linnet::{
    half_edge::{
        HedgeGraph, NodeIndex,
        builder::HedgeGraphBuilder,
        involution::{Flow, Orientation},
    },
    permutation::Permutation,
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

/// A generated diagram with shared incidence and a composed vertex/edge frame.
///
/// Vertex and edge identifiers refer to the current frame. Reordering, flow normalization and
/// CP rule repair change that frame without reconstructing the topology. A self-loop retains
/// both incident ends, and each external leg retains its attached interaction vertex.
#[derive(Debug, Clone)]
pub(crate) struct GeneratedGraph {
    topology: Arc<GeneratedTopology>,
    frame: Arc<GeneratedFrame>,
}

#[derive(Debug)]
struct GeneratedTopology {
    edges: Box<[GeneratedEdge]>,
    offsets: Box<[usize]>,
    incident: Box<[(usize, Flow)]>,
}

#[derive(Debug, Clone)]
struct GeneratedFrame {
    /// Base indices to current indices, composed directly rather than through nested views.
    vertices: Permutation,
    edges: Permutation,
    rules: Vec<VertexRuleId>,
    edge_data: Vec<FramedEdge>,
}

#[derive(Debug, Clone, Copy)]
struct FramedEdge {
    /// Relative to the immutable base edge, including for undirected edges and self-loops.
    reversed: bool,
    particle: ParticleId,
    external: Option<ExternalNode>,
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
    /// Retain a topology from its vertex rules and edges, in identifier order.
    ///
    /// # Panics
    /// Panics if an edge has no vertex endpoint or references an unknown vertex.
    pub(crate) fn new(
        vertices: impl IntoIterator<Item = VertexRuleId>,
        edges: impl IntoIterator<Item = GeneratedEdge>,
    ) -> Self {
        let rules = vertices.into_iter().collect::<Vec<_>>();
        let edges = edges.into_iter().collect::<Box<[_]>>();
        let mut offsets = vec![0; rules.len() + 1];
        for edge in &edges {
            assert!(
                matches!(edge.source, Endpoint::Vertex(_))
                    || matches!(edge.sink, Endpoint::Vertex(_)),
                "an edge needs an interaction vertex"
            );
            for endpoint in [edge.source, edge.sink] {
                if let Endpoint::Vertex(vertex) = endpoint {
                    offsets[vertex + 1] += 1;
                }
            }
        }
        for vertex in 0..rules.len() {
            offsets[vertex + 1] += offsets[vertex];
        }
        let mut incident = vec![(0, Flow::Source); offsets[rules.len()]];
        let mut next = offsets[..rules.len()].to_vec();
        for (id, edge) in edges.iter().enumerate() {
            for (endpoint, flow) in [(edge.source, Flow::Source), (edge.sink, Flow::Sink)] {
                if let Endpoint::Vertex(vertex) = endpoint {
                    incident[next[vertex]] = (id, flow);
                    next[vertex] += 1;
                }
            }
        }
        let frame = GeneratedFrame {
            vertices: Permutation::id(rules.len()),
            edges: Permutation::id(edges.len()),
            rules,
            edge_data: edges
                .iter()
                .map(|edge| FramedEdge {
                    reversed: false,
                    particle: edge.particle,
                    external: edge.source.external().or_else(|| edge.sink.external()),
                })
                .collect(),
        };
        Self {
            topology: Arc::new(GeneratedTopology {
                edges,
                offsets: offsets.into_boxed_slice(),
                incident: incident.into_boxed_slice(),
            }),
            frame: Arc::new(frame),
        }
    }

    /// Adopt a new frame without rebuilding incidence. `vertices` maps current vertices to
    /// their new indices; `edges` lists current edge IDs in new order, with an endpoint reversal
    /// and the particle in that orientation. `external` maps the current external labels.
    ///
    /// # Panics
    /// Panics unless the vertex and edge maps are permutations of the current identifiers.
    pub(crate) fn reframed(
        &self,
        vertices: &[usize],
        edges: impl IntoIterator<Item = (usize, bool, ParticleId)>,
        external: impl Fn(ExternalNode) -> ExternalNode,
    ) -> Self {
        assert_eq!(vertices.len(), self.n_vertices());
        let edges = edges.into_iter().collect::<Vec<_>>();
        assert_eq!(edges.len(), self.n_edges());
        if vertices.iter().enumerate().all(|(old, &new)| old == new)
            && edges
                .iter()
                .enumerate()
                .all(|(new, &(old, reversed, particle))| {
                    old == new && !reversed && particle == self.frame.edge_data[old].particle
                })
            && self
                .frame
                .edge_data
                .iter()
                .filter_map(|edge| edge.external)
                .all(|leg| external(leg) == leg)
        {
            return self.clone();
        }
        let mut rules = vec![None; vertices.len()];
        for (old, &new) in vertices.iter().enumerate() {
            assert!(rules[new].replace(self.frame.rules[old]).is_none());
        }
        let vertices = Permutation::from_map(vertices.to_vec()).compose(&self.frame.vertices);
        let mut positions = vec![usize::MAX; self.n_edges()];
        let mut edge_data = Vec::with_capacity(self.n_edges());
        for (old, reversed, particle) in edges {
            let base = self.frame.edges.inv()[old];
            assert_eq!(
                std::mem::replace(&mut positions[base], edge_data.len()),
                usize::MAX,
                "an edge appears twice in a frame"
            );
            let old = self.frame.edge_data[old];
            edge_data.push(FramedEdge {
                reversed: old.reversed ^ reversed,
                particle,
                external: old.external.map(&external),
            });
        }
        assert_eq!(edge_data.len(), self.n_edges());
        Self {
            topology: Arc::clone(&self.topology),
            frame: Arc::new(GeneratedFrame {
                vertices,
                edges: Permutation::from_map(positions),
                rules: rules
                    .into_iter()
                    .map(|rule| rule.expect("vertex positions are a permutation"))
                    .collect(),
                edge_data,
            }),
        }
    }

    /// Replace rules in current vertex order, retaining this topology and coordinate frame.
    pub(crate) fn with_vertex_rules(
        mut self,
        rules: impl IntoIterator<Item = VertexRuleId>,
    ) -> Self {
        let rules = rules.into_iter().collect::<Vec<_>>();
        assert_eq!(rules.len(), self.n_vertices());
        Arc::make_mut(&mut self.frame).rules = rules;
        self
    }

    pub(crate) fn n_vertices(&self) -> usize {
        self.frame.rules.len()
    }

    pub(crate) fn n_edges(&self) -> usize {
        self.topology.edges.len()
    }

    /// The vertex rules in current vertex order.
    pub(crate) fn vertices(&self) -> impl Iterator<Item = VertexRuleId> + '_ {
        self.frame.rules.iter().copied()
    }

    pub(crate) fn edges(&self) -> impl Iterator<Item = GeneratedEdge> + '_ {
        (0..self.n_edges()).map(|edge| self.edge(edge))
    }

    pub(crate) fn edge(&self, edge: usize) -> GeneratedEdge {
        let base = self.topology.edges[self.frame.edges.inv()[edge]];
        let data = self.frame.edge_data[edge];
        let endpoint = |endpoint| match endpoint {
            Endpoint::Vertex(vertex) => Endpoint::Vertex(self.frame.vertices.map()[vertex]),
            Endpoint::External(_) => Endpoint::External(data.external.expect("external leg")),
        };
        let mut source = endpoint(base.source);
        let mut sink = endpoint(base.sink);
        if data.reversed {
            std::mem::swap(&mut source, &mut sink);
        }
        GeneratedEdge {
            source,
            sink,
            particle: data.particle,
            directed: base.directed,
        }
    }

    /// Incident ends in the half-edge order of a graph built in the current frame: increasing
    /// edge IDs, source before sink. Reversing a self-loop must not reverse its UFO slot order.
    pub(crate) fn incident(&self, vertex: usize) -> Vec<(usize, Flow)> {
        let base = self.frame.vertices.inv()[vertex];
        let mut incident = self.topology.incident
            [self.topology.offsets[base]..self.topology.offsets[base + 1]]
            .iter()
            .map(|&(edge, flow)| {
                let edge = self.frame.edges.map()[edge];
                let flow = if self.frame.edge_data[edge].reversed {
                    -flow
                } else {
                    flow
                };
                (edge, flow)
            })
            .collect::<Vec<_>>();
        incident.sort_unstable_by_key(|&(edge, flow)| (edge, flow.is_sink()));
        incident
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
    #[test]
    fn composed_frames_keep_incidence_and_current_external_labels() {
        use std::collections::hash_map::DefaultHasher;

        let a = ExternalNode {
            index: 0,
            state: ExternalState::Incoming,
            particle: ParticleId::from_index(3),
        };
        let b = ExternalNode {
            index: 1,
            state: ExternalState::Outgoing,
            particle: ParticleId::from_index(7),
        };
        let edge = |source, sink, particle, directed| GeneratedEdge {
            source,
            sink,
            particle: ParticleId::from_index(particle),
            directed,
        };
        let rules = [0, 1, 2].map(VertexRuleId::from_index);
        let graph = GeneratedGraph::new(
            rules,
            [
                edge(Endpoint::External(a), Endpoint::Vertex(0), 3, true),
                edge(Endpoint::Vertex(0), Endpoint::Vertex(1), 4, false),
                edge(Endpoint::Vertex(1), Endpoint::Vertex(1), 5, true),
                edge(Endpoint::Vertex(1), Endpoint::Vertex(2), 6, true),
                edge(Endpoint::Vertex(2), Endpoint::External(b), 7, true),
                edge(Endpoint::Vertex(0), Endpoint::Vertex(1), 8, true),
            ],
        );
        let first = graph.reframed(
            &[2, 0, 1],
            [
                (4, true, 17),
                (2, true, 15),
                (0, false, 13),
                (5, true, 18),
                (3, false, 16),
                (1, true, 14),
            ]
            .map(|(edge, reversed, particle)| (edge, reversed, ParticleId::from_index(particle))),
            |leg| if leg == a { b } else { a },
        );
        // The reversed self-loop still occupies its current source slot before its sink slot.
        assert_eq!(
            first.incident(0),
            [
                (1, Flow::Source),
                (1, Flow::Sink),
                (3, Flow::Source),
                (4, Flow::Source),
                (5, Flow::Source)
            ]
        );
        let second = first.reframed(
            &[1, 2, 0],
            [
                (2, true, 23),
                (5, false, 24),
                (1, true, 25),
                (4, true, 26),
                (0, false, 27),
                (3, true, 28),
            ]
            .map(|(edge, reversed, particle)| (edge, reversed, ParticleId::from_index(particle))),
            |leg| ExternalNode {
                index: leg.index + 10,
                ..leg
            },
        );
        // Literal final edges independently specify permutation composition, orientation,
        // particle replacement, parallel edges and the two composed external-label changes.
        let expected_edges = [
            edge(
                Endpoint::Vertex(0),
                Endpoint::External(ExternalNode { index: 11, ..b }),
                23,
                true,
            ),
            edge(Endpoint::Vertex(1), Endpoint::Vertex(0), 24, false),
            edge(Endpoint::Vertex(1), Endpoint::Vertex(1), 25, true),
            edge(Endpoint::Vertex(2), Endpoint::Vertex(1), 26, true),
            edge(
                Endpoint::External(ExternalNode { index: 10, ..a }),
                Endpoint::Vertex(2),
                27,
                true,
            ),
            edge(Endpoint::Vertex(0), Endpoint::Vertex(1), 28, true),
        ];
        let expected = GeneratedGraph::new(rules, expected_edges);
        assert_eq!(second.edges().collect::<Vec<_>>(), expected_edges);
        assert_eq!(second.vertices().collect::<Vec<_>>(), rules);
        assert_eq!(second, expected);
        let hash = |graph: &GeneratedGraph| {
            let mut state = DefaultHasher::new();
            graph.hash(&mut state);
            state.finish()
        };
        assert_eq!(hash(&second), hash(&expected));
        assert!(Arc::ptr_eq(&graph.topology, &first.topology));
        assert!(Arc::ptr_eq(&graph.topology, &second.topology));

        // Compare the interpreted incidence with native half-edge insertion order, including
        // source/sink flow at each interaction endpoint of a dangling leg.
        let materialized = expected.with_external_vertices();
        for vertex in 0..second.n_vertices() {
            let mut hedges = materialized
                .iter_crown(NodeIndex(vertex))
                .collect::<Vec<_>>();
            hedges.sort_unstable();
            let native = hedges
                .into_iter()
                .map(|hedge| (materialized[&hedge].0, materialized.flow(hedge)))
                .collect::<Vec<_>>();
            assert_eq!(second.incident(vertex), native);
        }
        assert_eq!(
            second.externals().collect::<Vec<_>>(),
            expected.externals().collect::<Vec<_>>()
        );
        let plain = second.plain_nodes();
        assert_eq!(plain.len(), 5);
        assert_eq!(plain.ends(0, &second.edge(0)), (0, 3));
        assert_eq!(plain.ends(4, &second.edge(4)), (4, 2));
        assert_eq!(plain.external(3), Some(ExternalNode { index: 11, ..b }));
        assert_eq!(plain.external(4), Some(ExternalNode { index: 10, ..a }));

        let unchanged = second.reframed(
            &[0, 1, 2],
            second
                .edges()
                .enumerate()
                .map(|(id, edge)| (id, false, edge.particle)),
            |leg| leg,
        );
        assert!(Arc::ptr_eq(&second.frame, &unchanged.frame));
        let changed_rules = [12, 13, 14].map(VertexRuleId::from_index);
        let repaired = unchanged.with_vertex_rules(changed_rules);
        assert!(Arc::ptr_eq(&second.topology, &repaired.topology));
        assert_eq!(repaired.vertices().collect::<Vec<_>>(), changed_rules);
        assert_eq!(repaired.edges().collect::<Vec<_>>(), expected_edges);
        assert_eq!(second.vertices().collect::<Vec<_>>(), rules);
    }

    #[test]
    fn noncommuting_vertex_frames_preserve_intermediate_rule_repairs() {
        let particle = ParticleId::from_index(1);
        let graph = GeneratedGraph::new(
            [10, 20, 30].map(VertexRuleId::from_index),
            [(0, 1), (1, 2)].map(|(source, sink)| GeneratedEdge {
                source: Endpoint::Vertex(source),
                sink: Endpoint::Vertex(sink),
                particle,
                directed: true,
            }),
        );
        let edges = [(0, false, particle), (1, false, particle)];
        let first = graph.reframed(&[2, 0, 1], edges, |leg| leg);
        let repaired = first
            .clone()
            .with_vertex_rules([21, 31, 11].map(VertexRuleId::from_index));
        let second = repaired.reframed(&[1, 0, 2], edges, |leg| leg);
        assert_eq!(
            second.vertices().collect::<Vec<_>>(),
            [31, 21, 11].map(VertexRuleId::from_index)
        );
        assert_eq!(second.edge(0).source, Endpoint::Vertex(2));
        assert_eq!(second.edge(0).sink, Endpoint::Vertex(1));
        assert_eq!(second.edge(1).source, Endpoint::Vertex(1));
        assert_eq!(second.edge(1).sink, Endpoint::Vertex(0));
        assert_eq!(second.incident(0), [(1, Flow::Sink)]);
        assert_eq!(second.incident(1), [(0, Flow::Sink), (1, Flow::Source)]);
        assert_eq!(second.incident(2), [(0, Flow::Source)]);
        assert_eq!(
            first.vertices().collect::<Vec<_>>(),
            [20, 30, 10].map(VertexRuleId::from_index)
        );
        assert!(Arc::ptr_eq(&graph.topology, &second.topology));
    }

    #[test]
    fn frames_keep_isolated_vertices_and_empty_graphs() {
        let rules = [1, 4].map(VertexRuleId::from_index);
        let graph = GeneratedGraph::new(rules, []);
        let swapped = graph.reframed(&[1, 0], [], |leg| leg);
        assert_eq!(swapped.vertices().collect::<Vec<_>>(), [rules[1], rules[0]]);
        assert_eq!(swapped.n_edges(), 0);
        assert!(swapped.incident(0).is_empty());
        assert!(swapped.incident(1).is_empty());
        assert!(Arc::ptr_eq(&graph.topology, &swapped.topology));
        let empty = GeneratedGraph::new([], []);
        assert_eq!(empty.reframed(&[], [], |leg| leg), empty);
    }
}
