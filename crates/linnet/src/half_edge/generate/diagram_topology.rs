//! Immutable output incidence shared by decorations of one completed topology.

use std::sync::OnceLock;

use crate::half_edge::involution::Flow;

/// An owned graph shape in native edge order, independently of fields and rules.
///
/// A dangling edge has one absent endpoint. Its stored source/target roles are
/// coordinates only; a decoration can reverse them without rebuilding incidence.
#[derive(Debug)]
pub struct DiagramTopology {
    vertices: usize,
    endpoints: Box<[[Option<usize>; 2]]>,
    incidence: OnceLock<Incidence>,
    pub(super) single_diagram: bool,
}

#[derive(Debug)]
struct Incidence {
    offsets: Box<[usize]>,
    ends: Box<[(usize, Flow)]>,
}

impl DiagramTopology {
    pub fn new(vertices: usize, endpoints: impl IntoIterator<Item = [Option<usize>; 2]>) -> Self {
        let endpoints = endpoints.into_iter().collect::<Box<[_]>>();
        for edge in &endpoints {
            assert!(
                edge.iter().any(Option::is_some),
                "an edge needs an attached endpoint"
            );
            assert!(
                edge.iter().flatten().all(|&vertex| vertex < vertices),
                "unknown endpoint vertex"
            );
        }
        Self {
            vertices,
            endpoints,
            incidence: OnceLock::new(),
            single_diagram: false,
        }
    }

    pub fn vertex_count(&self) -> usize {
        self.vertices
    }

    pub fn edge_count(&self) -> usize {
        self.endpoints.len()
    }

    pub fn endpoints(&self, edge: usize) -> [Option<usize>; 2] {
        self.endpoints[edge]
    }

    /// Both ends of a self-loop occur, in edge order with source before sink.
    pub fn incident(&self, vertex: usize) -> &[(usize, Flow)] {
        let incidence = self.incidence.get_or_init(|| {
            let mut offsets = vec![0; self.vertices + 1];
            for vertex in self.endpoints.iter().flatten().flatten() {
                offsets[*vertex + 1] += 1;
            }
            for vertex in 0..self.vertices {
                offsets[vertex + 1] += offsets[vertex];
            }
            let mut ends = vec![(0, Flow::Source); offsets[self.vertices]];
            let mut next = offsets[..self.vertices].to_vec();
            for (edge, endpoints) in self.endpoints.iter().enumerate() {
                for (vertex, flow) in endpoints.iter().zip([Flow::Source, Flow::Sink]) {
                    if let Some(vertex) = vertex {
                        ends[next[*vertex]] = (edge, flow);
                        next[*vertex] += 1;
                    }
                }
            }
            Incidence {
                offsets: offsets.into_boxed_slice(),
                ends: ends.into_boxed_slice(),
            }
        });
        &incidence.ends[incidence.offsets[vertex]..incidence.offsets[vertex + 1]]
    }

    /// Whether the generation model can decorate this shape more than once.
    pub fn may_have_multiple_decorations(&self) -> bool {
        !self.single_diagram
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn incidence_retains_self_loop_twins_and_dangling_roles() {
        let topology =
            DiagramTopology::new(2, [[Some(0), Some(0)], [None, Some(0)], [Some(1), None]]);
        assert_eq!(
            topology.incident(0),
            [(0, Flow::Source), (0, Flow::Sink), (1, Flow::Sink)]
        );
        assert_eq!(topology.incident(1), [(2, Flow::Source)]);
        assert_eq!(topology.incident(0).as_ptr(), topology.incident(0).as_ptr());
    }

    #[test]
    fn native_sibling_decorations_share_shape_and_keep_distinct_fields() {
        use crate::half_edge::generate::{GraphGenerator, Stub};
        use std::sync::Arc;
        let a = Stub::undirected(0_u8);
        let b = Stub::incoming(1_u8);
        for (legs, vertices) in [(2, 2), (1, 1)] {
            let diagrams = GraphGenerator::new(
                (0..legs).map(|label| (label, a.clone())),
                [
                    vec![a.clone(); 3],
                    vec![b.clone(), b.clone().flipped(), a.clone()],
                ],
            )
            .loops(1..=1)
            .max_vertices(vertices)
            .allow_self_loops(vertices == 1)
            .iter()
            .unwrap()
            .collect::<Vec<_>>();
            assert!(diagrams.len() >= 2);
            let first = &diagrams[0];
            for diagram in &diagrams {
                assert!(Arc::ptr_eq(first.topology(), diagram.topology()));
                assert_eq!(diagram.rules().len(), vertices);
            }
            assert!(diagrams
                .iter()
                .any(|diagram| diagram.edges().any(|edge| *edge.stub.unwrap().field == 1)));
            assert!(diagrams
                .iter()
                .any(|diagram| diagram.edges().all(|edge| *edge.stub.unwrap().field == 0)));
            assert!(first.topology().may_have_multiple_decorations());
        }
    }
}
