//! Bounded request-local reuse of exact reduced-canonizer inputs.

use std::sync::{Arc, Mutex};

use linnet::half_edge::{NoData, canon::CanonicalForm, generate::DiagramTopology};

use super::graph::{Endpoint, GeneratedEdge, GeneratedGraph};

type ColoredEdge = (GeneratedEdge, usize, Option<[i64; 2]>);
const PLANS_PER_TOPOLOGY: usize = 8;

pub(super) struct CanonicalNodes {
    pub(super) form: CanonicalForm<NoData, usize, Option<i64>>,
    pub(super) positions: Vec<usize>,
}

/// The complete ordered `KeyedEdge` input, including the current coordinate frame.
#[derive(PartialEq, Eq)]
enum InputEdge {
    Paired {
        source: usize,
        sink: usize,
        color: usize,
    },
    Dangling {
        vertex: usize,
        color: usize,
        label: Option<i64>,
    },
}

impl InputEdge {
    fn new((edge, color, labels): &ColoredEdge, swapped: bool) -> Self {
        match (edge.source, edge.sink) {
            (Endpoint::Vertex(source), Endpoint::Vertex(sink)) => Self::Paired {
                source,
                sink,
                color: *color,
            },
            (Endpoint::Vertex(vertex), Endpoint::External(_))
            | (Endpoint::External(_), Endpoint::Vertex(vertex)) => Self::Dangling {
                vertex,
                color: *color,
                label: labels.map(|labels| labels[usize::from(swapped)]),
            },
            _ => unreachable!("generated edges touch an interaction vertex"),
        }
    }
}

struct Entry {
    flags: [bool; 4],
    input: Vec<InputEdge>,
    nodes: Arc<CanonicalNodes>,
}

impl Entry {
    fn matches(&self, edges: &[ColoredEdge], flags: [bool; 4]) -> bool {
        self.flags == flags
            && self.input.len() == edges.len()
            && self
                .input
                .iter()
                .zip(edges)
                .all(|(input, edge)| *input == InputEdge::new(edge, flags[3]))
    }
}

struct Family {
    shape: Arc<DiagramTopology>,
    entries: Vec<Entry>,
}

/// Retain only the latest shape and at most eight exact input plans. Sibling
/// decorations stream together; parallel interleaving may lose hits, never correctness.
#[derive(Default)]
pub(super) struct TopologyReuse {
    family: Mutex<Option<Family>>,
}

impl TopologyReuse {
    pub(super) fn get(
        &self,
        graph: &GeneratedGraph,
        edges: &[ColoredEdge],
        flags: [bool; 4],
    ) -> Option<Arc<CanonicalNodes>> {
        if !graph.topology().may_have_multiple_decorations() {
            return None;
        }
        let mut family = self.family.lock().expect("topology reuse poisoned");
        if !family
            .as_ref()
            .is_some_and(|family| Arc::ptr_eq(&family.shape, graph.topology()))
        {
            *family = Some(Family {
                shape: Arc::clone(graph.topology()),
                entries: Vec::new(),
            });
        }
        family
            .as_ref()
            .expect("current shape")
            .entries
            .iter()
            .find(|entry| entry.matches(edges, flags))
            .map(|entry| Arc::clone(&entry.nodes))
    }

    pub(super) fn insert(
        &self,
        graph: &GeneratedGraph,
        edges: &[ColoredEdge],
        flags: [bool; 4],
        nodes: &Arc<CanonicalNodes>,
    ) {
        if !graph.topology().may_have_multiple_decorations() {
            return;
        }
        let mut family = self.family.lock().expect("topology reuse poisoned");
        let Some(family) = family
            .as_mut()
            .filter(|family| Arc::ptr_eq(&family.shape, graph.topology()))
        else {
            return;
        };
        if family.entries.len() < PLANS_PER_TOPOLOGY
            && !family
                .entries
                .iter()
                .any(|entry| entry.matches(edges, flags))
        {
            family.entries.push(Entry {
                flags,
                input: edges
                    .iter()
                    .map(|edge| InputEdge::new(edge, flags[3]))
                    .collect(),
                nodes: Arc::clone(nodes),
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use feynkit_model::{ParticleId, VertexRuleId};
    use linnet::half_edge::{NodeIndex, canon::KeyedEdge, involution::Orientation};

    fn graph() -> GeneratedGraph {
        let particle = ParticleId::from_index(0);
        GeneratedGraph::new(
            [VertexRuleId::from_index(0); 2],
            [GeneratedEdge {
                source: Endpoint::Vertex(0),
                sink: Endpoint::Vertex(1),
                particle,
                directed: false,
            }],
        )
    }

    fn colored(graph: &GeneratedGraph) -> Vec<ColoredEdge> {
        graph.edges().map(|edge| (edge, 0, None)).collect()
    }

    fn plan() -> Arc<CanonicalNodes> {
        let (form, nodes) = CanonicalForm::label_edges(
            [NoData {}; 2],
            [KeyedEdge::Paired {
                source: (NodeIndex(0), None),
                sink: (NodeIndex(1), None),
                key: 0,
                orientation: Orientation::Undirected,
            }],
        )
        .unwrap();
        Arc::new(CanonicalNodes {
            form,
            positions: nodes.into_iter().map(|(_, node)| node.0).collect(),
        })
    }

    #[test]
    fn reuse_requires_exact_current_frame_colors_and_flags() {
        let graph = graph();
        let edges = colored(&graph);
        let cache = TopologyReuse::default();
        let flags = [false; 4];
        assert!(cache.get(&graph, &edges, flags).is_none());
        let nodes = plan();
        cache.insert(&graph, &edges, flags, &nodes);
        assert!(Arc::ptr_eq(
            &nodes,
            &cache.get(&graph, &edges, flags).unwrap()
        ));
        let mut changed = edges.clone();
        changed[0].1 = 1;
        assert!(cache.get(&graph, &changed, flags).is_none());
        assert!(
            cache
                .get(&graph, &edges, [true, false, false, false])
                .is_none()
        );
        let reframed = graph.reframed(&[1, 0], [(0, false, edges[0].0.particle)], |external| {
            external
        });
        assert!(Arc::ptr_eq(graph.topology(), reframed.topology()));
        assert!(cache.get(&reframed, &colored(&reframed), flags).is_none());
    }

    #[test]
    fn moving_to_the_next_topology_releases_the_previous_shape() {
        let cache = TopologyReuse::default();
        let first = graph();
        let old = Arc::downgrade(first.topology());
        let edges = colored(&first);
        let flags = [false; 4];
        assert!(cache.get(&first, &edges, flags).is_none());
        cache.insert(&first, &edges, flags, &plan());
        drop(first);
        assert!(old.upgrade().is_some());
        let next = graph();
        assert!(cache.get(&next, &colored(&next), flags).is_none());
        assert!(old.upgrade().is_none());
    }
}
