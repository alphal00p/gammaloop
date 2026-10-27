//! Immutable output incidence shared by the decorations of one complete skeleton.

use std::sync::OnceLock;

use crate::half_edge::involution::Flow;

/// The edge endpoints of a skeleton's diagrams in native edge order (external legs first, then
/// internal edges in decorator order), independently of fields and rules.
///
/// A dangling edge has one absent endpoint. Its stored source/sink roles are coordinates only;
/// a decoration can reverse them without rebuilding incidence.
#[derive(Debug)]
pub struct SkeletonIncidence {
    vertices: usize,
    endpoints: Box<[[Option<usize>; 2]]>,
    adjacency: OnceLock<Adjacency>,
}

/// Compressed vertex adjacency, built on first demand.
#[derive(Debug)]
struct Adjacency {
    offsets: Box<[usize]>,
    ends: Box<[(usize, Flow)]>,
}

impl SkeletonIncidence {
    /// # Panics
    /// Panics if an edge has no endpoint or references an unknown vertex.
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
            adjacency: OnceLock::new(),
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
        let adjacency = self.adjacency.get_or_init(|| {
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
            Adjacency {
                offsets: offsets.into_boxed_slice(),
                ends: ends.into_boxed_slice(),
            }
        });
        &adjacency.ends[adjacency.offsets[vertex]..adjacency.offsets[vertex + 1]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn incidence_retains_self_loop_twins_and_dangling_roles() {
        let incidence =
            SkeletonIncidence::new(2, [[Some(0), Some(0)], [None, Some(0)], [Some(1), None]]);
        assert_eq!(
            incidence.incident(0),
            [(0, Flow::Source), (0, Flow::Sink), (1, Flow::Sink)]
        );
        assert_eq!(incidence.incident(1), [(2, Flow::Source)]);
        assert_eq!(
            incidence.incident(0).as_ptr(),
            incidence.incident(0).as_ptr()
        );
    }
}
