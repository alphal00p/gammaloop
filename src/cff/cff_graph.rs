use ahash::{HashMap, HashSet, HashSetExt};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::hash::Hash;

use crate::{
    cff::hsurface::Hsurface,
    graph::{BareGraph, EdgeType},
};

use super::{
    esurface::{Esurface, ExternalShift},
    surface::HybridSurface,
};

const MAX_VERTEX_COUNT: usize = 64;

#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
struct CFFVertex {
    vertex_set: VertexSet,
    pub incoming_edges: Vec<CFFEdge>,
    pub outgoing_edges: Vec<CFFEdge>,
}

impl CFFVertex {
    fn new(id: usize) -> Self {
        let vertex_set = VertexSet::from_usize(id);

        CFFVertex {
            vertex_set,
            incoming_edges: Vec::new(),
            outgoing_edges: Vec::new(),
        }
    }

    fn get_vertex_type(&self) -> VertexType {
        let is_sink = self
            .outgoing_edges
            .iter()
            .all(|edge| edge.edge_type != CFFEdgeType::Virtual);

        let is_source = self
            .incoming_edges
            .iter()
            .all(|edge| edge.edge_type != CFFEdgeType::Virtual);

        if is_sink {
            VertexType::Sink
        } else if is_source {
            VertexType::Source
        } else {
            VertexType::None
        }
    }

    fn generates_esurface(&self) -> bool {
        let vertex_type = self.get_vertex_type();

        match vertex_type {
            VertexType::None => false,
            VertexType::Sink => self
                .outgoing_edges
                .iter()
                .all(|edge| edge.edge_type == CFFEdgeType::External),
            VertexType::Source => self
                .incoming_edges
                .iter()
                .all(|edge| edge.edge_type == CFFEdgeType::External),
        }
    }

    fn join(&self, other: &Self) -> Self {
        let new_vertex_set = self.vertex_set.join(&other.vertex_set);
        // find all parrallel edges between the vertices

        let incoming_edges_of_new = self
            .incoming_edges
            .iter()
            .filter(|edge| !other.outgoing_edges.contains(edge))
            .chain(
                other
                    .incoming_edges
                    .iter()
                    .filter(|edge| !self.outgoing_edges.contains(edge)),
            )
            .copied()
            .sorted_by(|edge_1, edge_2| edge_1.edge_id.cmp(&edge_2.edge_id))
            .collect_vec();

        let outgoing_edges_of_new = self
            .outgoing_edges
            .iter()
            .filter(|edge| !other.incoming_edges.contains(edge))
            .chain(
                other
                    .outgoing_edges
                    .iter()
                    .filter(|edge| !self.incoming_edges.contains(edge)),
            )
            .copied()
            .sorted_by(|edge_1, edge_2| edge_1.edge_id.cmp(&edge_2.edge_id))
            .collect_vec();

        CFFVertex {
            vertex_set: new_vertex_set,
            incoming_edges: incoming_edges_of_new,
            outgoing_edges: outgoing_edges_of_new,
        }
    }

    fn iter_all_edges(&self) -> impl Iterator<Item = &CFFEdge> {
        self.incoming_edges.iter().chain(self.outgoing_edges.iter())
    }

    fn iter_all_edges_mut(&mut self) -> impl Iterator<Item = &mut CFFEdge> {
        self.incoming_edges
            .iter_mut()
            .chain(self.outgoing_edges.iter_mut())
    }

    fn has_edge(&self, edge_id: usize) -> bool {
        self.iter_all_edges().any(|edge| edge.edge_id == edge_id)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct VertexSet {
    vertex_set: u64,
}

impl VertexSet {
    pub fn from_usize(id: usize) -> Self {
        assert!(id < MAX_VERTEX_COUNT, "Vertex ID out of bounds");

        VertexSet {
            vertex_set: 1 << id,
        }
    }

    fn join(&self, other: &VertexSet) -> VertexSet {
        assert!(
            self.vertex_set & other.vertex_set == 0,
            "Vertex sets overlap"
        );

        VertexSet {
            vertex_set: self.vertex_set | other.vertex_set,
        }
    }

    fn contains_vertices(&self) -> Vec<VertexSet> {
        (0..MAX_VERTEX_COUNT)
            .filter(|id| self.vertex_set & (1 << id) != 0)
            .map(VertexSet::from_usize)
            .collect()
    }

    pub fn dummy() -> Self {
        VertexSet { vertex_set: 0 }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
struct CFFEdge {
    edge_id: usize,
    edge_type: CFFEdgeType,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CFFEdgeType {
    External,
    Virtual,
    VirtualExternal,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum VertexType {
    Source,
    Sink,
    None,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CFFGenerationGraph {
    vertices: Vec<CFFVertex>,
    pub global_orientation: Vec<bool>,
}

impl PartialEq for CFFGenerationGraph {
    fn eq(&self, other: &Self) -> bool {
        self.vertices == other.vertices
    }
}

impl Eq for CFFGenerationGraph {}

impl Hash for CFFGenerationGraph {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.vertices.hash(state);
    }
}

impl CFFGenerationGraph {
    fn get_vertex(&self, vertex_set: &VertexSet) -> &CFFVertex {
        self.vertices
            .iter()
            .find(|node| node.vertex_set == *vertex_set)
            .unwrap_or_else(|| panic!("Vertex not found"))
    }

    fn are_directed_adjacent(&self, left: &VertexSet, right: &VertexSet) -> bool {
        let left_vertex = self.get_vertex(left);
        let right_vertex = self.get_vertex(right);

        left_vertex.outgoing_edges.iter().any(|outgoing_edge| {
            right_vertex
                .incoming_edges
                .iter()
                .any(|right_incoming| right_incoming.edge_id == outgoing_edge.edge_id)
        })
    }

    fn are_adjacent(&self, vertex_1: &VertexSet, vertex_2: &VertexSet) -> bool {
        self.are_directed_adjacent(vertex_1, vertex_2)
            || self.are_directed_adjacent(vertex_2, vertex_1)
    }

    // helper function for tests
    #[cfg(test)]
    pub fn from_vec(
        edges: Vec<(usize, usize)>,
        incoming_vertices: Vec<(usize, CFFEdgeType)>,
        virtual_orientation: Vec<bool>,
    ) -> Self {
        let edges = edges
            .into_iter()
            .map(|(from, to)| (VertexSet::from_usize(from), VertexSet::from_usize(to)))
            .collect_vec();

        let mut unique_vertex_ids = HashSet::new();
        let mut global_orientation = vec![];

        for edge in edges.iter() {
            unique_vertex_ids.insert(edge.0);
            unique_vertex_ids.insert(edge.1);
        }

        let mut unique_vertices = HashMap::default();

        for vertex in unique_vertex_ids.iter() {
            let cff_vertex = CFFVertex {
                vertex_set: *vertex,
                incoming_edges: Vec::new(),
                outgoing_edges: Vec::new(),
            };

            unique_vertices.insert(*vertex, cff_vertex);
        }

        for (edge_id, (incoming_vertex, edge_type)) in incoming_vertices.iter().enumerate() {
            assert_ne!(*edge_type, CFFEdgeType::Virtual);

            let vertex_set = VertexSet::from_usize(*incoming_vertex);
            let cff_vertex = unique_vertices.get_mut(&vertex_set).unwrap();
            let incoming_edge = CFFEdge {
                edge_id,
                edge_type: *edge_type,
            };
            cff_vertex.incoming_edges.push(incoming_edge);
            global_orientation.push(true);
        }

        for (edge_id, (left, right)) in edges.iter().enumerate() {
            let edge_id = edge_id + incoming_vertices.len();
            let cff_edge = CFFEdge {
                edge_id,
                edge_type: CFFEdgeType::Virtual,
            };

            let left_vertex = unique_vertices.get_mut(left).unwrap();
            left_vertex.outgoing_edges.push(cff_edge);

            let right_vertex = unique_vertices.get_mut(right).unwrap();
            right_vertex.incoming_edges.push(cff_edge);
        }

        global_orientation.extend(virtual_orientation);
        let nodes = unique_vertices.into_values().collect_vec();

        CFFGenerationGraph {
            vertices: nodes,
            global_orientation,
        }
    }

    fn get_vertex_that_is_not_v(&self, vertex: &VertexSet) -> &CFFVertex {
        self.vertices
            .iter()
            .find(|node| node.vertex_set != *vertex)
            .unwrap_or_else(|| panic!("Could not find vertex that is not v"))
    }

    fn get_directed_neighbours(&self, vertex: &VertexSet) -> Vec<&CFFVertex> {
        let outgoing_of_vertex = self.get_vertex(vertex).outgoing_edges.iter();

        outgoing_of_vertex
            .filter_map(|edge| {
                self.vertices
                    .iter()
                    .find(|node| node.incoming_edges.contains(edge))
            })
            .collect()
    }

    fn get_undirected_neighbours(&self, vertex: &VertexSet) -> Vec<&CFFVertex> {
        self.vertices
            .iter()
            .filter(|other_vertex| {
                self.are_adjacent(vertex, &other_vertex.vertex_set)
                    && other_vertex.vertex_set != *vertex
            })
            .collect()
    }

    fn depth_first_search(
        &self,
        vertex: &VertexSet,
        visited: &mut HashSet<VertexSet>,
        stack: &mut Vec<VertexSet>,
    ) -> bool {
        if visited.contains(vertex) {
            return stack.contains(vertex);
        }

        visited.insert(*vertex);
        stack.push(*vertex);

        let neighbours = self.get_directed_neighbours(vertex);
        for neighbour in neighbours.iter() {
            if self.depth_first_search(&neighbour.vertex_set, visited, stack) {
                return true;
            }
        }

        stack.pop();
        false
    }

    fn has_directed_cycle(&self, seed_vertex: &VertexSet) -> bool {
        let mut visited = HashSet::new();
        let mut stack = Vec::new();

        self.depth_first_search(seed_vertex, &mut visited, &mut stack)
    }

    pub fn has_directed_cycle_initial(&self) -> bool {
        self.vertices
            .iter()
            .any(|vertex| self.has_directed_cycle(&vertex.vertex_set))
    }

    fn has_connected_complement(&self, vertex: &VertexSet) -> bool {
        if self.vertices.len() == 1 {
            return true;
        }

        let vertex_that_is_not_v = self.get_vertex_that_is_not_v(vertex);

        let mut current_vertices = HashSet::default();
        current_vertices.insert(vertex_that_is_not_v.vertex_set);

        let mut visited_vertices = HashSet::default();
        visited_vertices.insert(&vertex_that_is_not_v.vertex_set);

        let mut vertices_found_in_previous_iteration = HashSet::default();
        let mut delta = 1;

        while delta > 0 {
            delta = 0; // reset delta

            for current_vertex in current_vertices.iter() {
                for adjacent_vertex in self
                    .get_undirected_neighbours(current_vertex)
                    .iter()
                    .filter(|v| v.vertex_set != *vertex)
                {
                    if visited_vertices.insert(&adjacent_vertex.vertex_set) {
                        delta += 1;
                        vertices_found_in_previous_iteration.insert(adjacent_vertex.vertex_set);
                    }
                }
            }

            // current_vertices = vertices_found_in_previous_iteration.clone();
            // vertices_found_in_previous_iteration.clear();
            std::mem::swap(
                &mut current_vertices,
                &mut vertices_found_in_previous_iteration,
            );
            vertices_found_in_previous_iteration.clear();
        }

        visited_vertices.len() == self.vertices.len() - 1
    }

    #[cfg(test)]
    fn get_vertex_type(&self, vertex: &VertexSet) -> VertexType {
        let vertex = self.get_vertex(vertex);
        vertex.get_vertex_type()
    }

    fn is_valid_source_or_sink(&self, vertex: VertexSet) -> bool {
        let vertex = self.get_vertex(&vertex);
        let vertex_type = vertex.get_vertex_type();
        if vertex_type == VertexType::None {
            return false;
        }

        self.has_connected_complement(&vertex.vertex_set)
    }

    #[allow(unused)]
    fn get_source_sink_candidate_list(&self) -> Vec<&CFFVertex> {
        self.vertices
            .iter()
            .filter(|vertex| {
                let vertex_type = vertex.get_vertex_type();
                let has_connected_complement = self.has_connected_complement(&vertex.vertex_set);
                (vertex_type == VertexType::Sink || vertex_type == VertexType::Source)
                    && has_connected_complement
            })
            .collect()
    }

    #[allow(unused)]
    fn get_source_sink_greedy(&self) -> Option<&CFFVertex> {
        self.vertices.iter().find(|vertex| {
            let vertex_type = vertex.get_vertex_type();
            if vertex_type != VertexType::Sink && vertex_type != VertexType::Source {
                false
            } else {
                self.has_connected_complement(&vertex.vertex_set)
            }
        })
    }

    fn contract_vertices(&self, vertex_1: &VertexSet, vertex_2: &VertexSet) -> Self {
        let vertex_1 = self.get_vertex(vertex_1);
        let vertex_2 = self.get_vertex(vertex_2);

        let new_vertex = vertex_1.join(vertex_2);

        let mut new_vertices = self.vertices.clone();

        new_vertices.retain(|vertex| vertex.vertex_set != vertex_1.vertex_set);
        new_vertices.retain(|vertex| vertex.vertex_set != vertex_2.vertex_set);
        new_vertices.push(new_vertex);
        new_vertices.shrink_to_fit();
        new_vertices.sort_by(|a, b| a.vertex_set.vertex_set.cmp(&b.vertex_set.vertex_set));

        CFFGenerationGraph {
            vertices: new_vertices,
            global_orientation: self.global_orientation.clone(),
        }
    }

    #[allow(unused)]
    fn get_source_or_sink_slow(&self) -> &CFFVertex {
        let mut source_sink_candidates = self.get_source_sink_candidate_list();

        source_sink_candidates
            .sort_by(|a, b| a.vertex_set.vertex_set.cmp(&b.vertex_set.vertex_set));

        if source_sink_candidates.is_empty() {
            panic!("No source or sink candidates found for graph {:#?}", self);
        }

        source_sink_candidates[0]
    }

    #[allow(unused)]
    fn get_source_or_sink_smart(&self, vertices_used: &mut Vec<VertexSet>) -> &CFFVertex {
        let mut vertices_checked = vec![];

        for vertex in vertices_used.iter() {
            if !self.vertices.iter().any(|v| v.vertex_set == *vertex) {
                continue;
            }

            let neighbours_of_vertex = self.get_undirected_neighbours(vertex);

            for candidate in neighbours_of_vertex.iter() {
                if vertices_checked.contains(&candidate.vertex_set) {
                    continue;
                } else {
                    let is_valid = self.is_valid_source_or_sink(candidate.vertex_set);
                    if is_valid {
                        vertices_used.push(candidate.vertex_set);
                        return candidate;
                    } else {
                        vertices_checked.push(candidate.vertex_set);
                    }
                }
            }
        }

        for vertex in self.vertices.iter() {
            if vertices_used.contains(&vertex.vertex_set) {
                continue;
            }

            let is_valid = self.is_valid_source_or_sink(vertex.vertex_set);
            if is_valid {
                return vertex;
            }
        }

        panic!("No source or sink candidates found for graph {:#?}", self);
    }

    #[allow(unused)]
    fn get_vertex_with_conn_complement(&self) -> &CFFVertex {
        self.vertices
            .iter()
            .find(|vertex| self.has_connected_complement(&vertex.vertex_set))
            .unwrap_or_else(|| panic!("Could not find vertex with connected complement"))
    }

    #[allow(unused_variables)]
    #[allow(clippy::ptr_arg)]
    pub fn generate_children(
        &self,
        vertices_used: &mut Vec<VertexSet>,
    ) -> (Option<Vec<Self>>, HybridSurface) {
        if self.vertices.len() < 2 {
            return (None, HybridSurface::Unit(()));
        }

        let vertex = if let Some(vertex) = self.get_source_sink_greedy() {
            vertex
        } else {
            panic!("could not find vertex to contract from")
        };

        let vertex_type = vertex.get_vertex_type();

        let external_shift: ExternalShift = vertex
            .incoming_edges
            .iter()
            .filter(|edge| edge.edge_type == CFFEdgeType::External)
            .map(|edge| {
                let edge_id = edge.edge_id;
                let shift_sign = match vertex_type {
                    VertexType::Source => -1,
                    VertexType::Sink => 1,
                    VertexType::None => panic!("vertex is not a source or a sink"),
                };
                (edge_id, shift_sign)
            })
            .chain(
                vertex
                    .outgoing_edges
                    .iter()
                    .filter(|edge| edge.edge_type == CFFEdgeType::External)
                    .map(|edge| {
                        let edge_id = edge.edge_id;
                        let shift_sign = match vertex_type {
                            VertexType::Source => 1,
                            VertexType::Sink => -1,
                            VertexType::None => panic!("vertex is not a source or a sink"),
                        };
                        (edge_id, shift_sign)
                    }),
            )
            .sorted_by(|(edge_1, _), (edge_2, _)| edge_1.cmp(edge_2))
            .collect_vec();

        let positive_energies = vertex
            .incoming_edges
            .iter()
            .chain(vertex.outgoing_edges.iter())
            .filter(|edge| edge.edge_type == CFFEdgeType::Virtual)
            .map(|edge| edge.edge_id)
            .sorted()
            .collect_vec();

        let surface = if vertex.generates_esurface() {
            let mut extra_positive_energies = vertex
                .iter_all_edges()
                .filter(|edge| edge.edge_type == CFFEdgeType::VirtualExternal)
                .map(|edge| edge.edge_id)
                .collect();

            let mut positive_energies = positive_energies;
            positive_energies.append(&mut extra_positive_energies);
            positive_energies.sort();

            let sub_orientation = positive_energies
                .iter()
                .map(|id| self.global_orientation[*id])
                .collect();

            let esurface = Esurface {
                energies: positive_energies,
                sub_orientation,
                external_shift,
                circled_vertices: vertex.vertex_set,
            };

            HybridSurface::Esurface(esurface)
        } else {
            let (mut extra_positive_energies, mut negative_energies) = match vertex_type {
                VertexType::Sink => (
                    vertex
                        .incoming_edges
                        .iter()
                        .filter(|edge| edge.edge_type == CFFEdgeType::VirtualExternal)
                        .map(|edge| edge.edge_id)
                        .collect(),
                    vertex
                        .outgoing_edges
                        .iter()
                        .filter(|edge| edge.edge_type == CFFEdgeType::VirtualExternal)
                        .map(|edge| edge.edge_id)
                        .collect_vec(),
                ),
                VertexType::Source => (
                    vertex
                        .outgoing_edges
                        .iter()
                        .filter(|edge| edge.edge_type == CFFEdgeType::VirtualExternal)
                        .map(|edge| edge.edge_id)
                        .collect(),
                    vertex
                        .incoming_edges
                        .iter()
                        .filter(|edge| edge.edge_type == CFFEdgeType::VirtualExternal)
                        .map(|edge| edge.edge_id)
                        .collect_vec(),
                ),
                VertexType::None => {
                    unreachable!()
                }
            };

            let mut positive_energies = positive_energies;
            positive_energies.append(&mut extra_positive_energies);
            positive_energies.sort();

            negative_energies.sort();

            let hsurface = Hsurface {
                positive_energies,
                negative_energies,
                external_shift,
            };

            HybridSurface::Hsurface(hsurface)
        };

        if self.vertices.len() > 2 {
            let adjacent_vertices = self.get_undirected_neighbours(&vertex.vertex_set);

            let children = adjacent_vertices
                .iter()
                .map(|adjacent_vertex| {
                    let contracted =
                        self.contract_vertices(&vertex.vertex_set, &adjacent_vertex.vertex_set);
                    let new_vertex_set = vertex.vertex_set.join(&adjacent_vertex.vertex_set);
                    (contracted, new_vertex_set)
                })
                .filter(|(graph, new_vertex_set)| !graph.has_directed_cycle(new_vertex_set))
                .map(|(graph, _)| graph)
                .collect_vec();

            (Some(children), surface)
        } else {
            (None, surface)
        }
    }

    pub fn new(graph: &BareGraph, virtual_orientation: Vec<bool>) -> Self {
        let virtual_index_to_edge_index = graph
            .get_virtual_edges_iterator()
            .map(|(_virtual_id, edge)| graph.get_edge_position(&edge.name).unwrap())
            .collect_vec();

        let edge_index_to_virtual_index = virtual_index_to_edge_index
            .iter()
            .enumerate()
            .map(|(index, &value)| (value, index))
            .collect::<HashMap<_, _>>();

        let global_orientation = graph
            .edges
            .iter()
            .map(|edge| match edge.edge_type {
                EdgeType::Virtual => {
                    let edge_id = graph.get_edge_position(&edge.name).unwrap();
                    virtual_orientation[edge_index_to_virtual_index[&edge_id]]
                }
                EdgeType::Incoming => true,
                EdgeType::Outgoing => true,
            })
            .collect_vec();

        let mut vertices = vec![];
        let mut cff_vertex_to_graph_vertex = vec![];
        let mut graph_vertex_to_cff_vertex = HashMap::default();

        graph
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, vertex)| vertex.edges.len() > 1)
            .enumerate()
            .for_each(|(cff_vertex_id, (graph_vertex_id, _vertex))| {
                let cff_vertex = CFFVertex::new(cff_vertex_id);
                vertices.push(cff_vertex);
                cff_vertex_to_graph_vertex.push(graph_vertex_id);
                graph_vertex_to_cff_vertex.insert(graph_vertex_id, cff_vertex_id);
            });

        for edge in graph.edges.iter() {
            let edge_id = graph.get_edge_position(&edge.name).unwrap();
            match edge.edge_type {
                EdgeType::Virtual => {
                    let from = edge.vertices[0];
                    let to = edge.vertices[1];
                    let orientation = global_orientation[edge_id];

                    if orientation {
                        vertices[graph_vertex_to_cff_vertex[&from]]
                            .outgoing_edges
                            .push(CFFEdge {
                                edge_id,
                                edge_type: CFFEdgeType::Virtual,
                            });
                        vertices[graph_vertex_to_cff_vertex[&to]]
                            .incoming_edges
                            .push(CFFEdge {
                                edge_id,
                                edge_type: CFFEdgeType::Virtual,
                            });
                    } else {
                        vertices[graph_vertex_to_cff_vertex[&from]]
                            .incoming_edges
                            .push(CFFEdge {
                                edge_id,
                                edge_type: CFFEdgeType::Virtual,
                            });
                        vertices[graph_vertex_to_cff_vertex[&to]]
                            .outgoing_edges
                            .push(CFFEdge {
                                edge_id,
                                edge_type: CFFEdgeType::Virtual,
                            });
                    }
                }
                EdgeType::Incoming => {
                    let vertex_id = edge.vertices[1];
                    vertices[graph_vertex_to_cff_vertex[&vertex_id]]
                        .incoming_edges
                        .push(CFFEdge {
                            edge_id,
                            edge_type: CFFEdgeType::External,
                        });
                }
                EdgeType::Outgoing => {
                    let vertex_id = edge.vertices[0];
                    vertices[graph_vertex_to_cff_vertex[&vertex_id]]
                        .outgoing_edges
                        .push(CFFEdge {
                            edge_id,
                            edge_type: CFFEdgeType::External,
                        });
                }
            }
        }

        Self {
            vertices,
            global_orientation,
        }
    }

    pub fn generate_cut(&self, circled_vertices: VertexSet) -> (Self, Self) {
        let vertices_in_cut = circled_vertices.contains_vertices();
        let mut vertices = self.vertices.clone();

        let mut left = vec![];

        for vertex_in_cut in vertices_in_cut.iter() {
            let vertex_position = vertices
                .iter()
                .position(|vertex| vertex.vertex_set == *vertex_in_cut)
                .unwrap();

            left.push(vertices.remove(vertex_position));
        }

        let mut left_graph = CFFGenerationGraph {
            vertices: left,
            global_orientation: self.global_orientation.clone(),
        };

        let mut right_graph = CFFGenerationGraph {
            vertices,
            global_orientation: self.global_orientation.clone(),
        };

        let edges_of_left_graph = left_graph.get_edges();

        let cut_edges = edges_of_left_graph
            .iter()
            .filter(|&&edge_id| right_graph.has_edge(edge_id))
            .copied()
            .collect_vec();

        left_graph
            .iter_all_edges_mut()
            .chain(right_graph.iter_all_edges_mut())
            .filter(|edge| cut_edges.contains(&edge.edge_id))
            .for_each(|edge_to_edit| edge_to_edit.edge_type = CFFEdgeType::VirtualExternal);

        (left_graph, right_graph)
    }

    fn get_edges(&self) -> Vec<usize> {
        let mut unique_edges = vec![];

        for vertex in self.vertices.iter() {
            for edge in vertex.iter_all_edges() {
                if !unique_edges.contains(&edge.edge_id) {
                    unique_edges.push(edge.edge_id)
                }
            }
        }

        // sort is always good
        unique_edges.sort();
        unique_edges
    }

    #[cfg(test)]
    /// has duplicates
    fn iter_all_edges(&self) -> impl Iterator<Item = &CFFEdge> {
        self.vertices
            .iter()
            .flat_map(|vertex| vertex.iter_all_edges())
    }

    fn iter_all_edges_mut(&mut self) -> impl Iterator<Item = &mut CFFEdge> {
        self.vertices
            .iter_mut()
            .flat_map(|vertex| vertex.iter_all_edges_mut())
    }

    fn has_edge(&self, edge_id: usize) -> bool {
        self.vertices.iter().any(|vertex| vertex.has_edge(edge_id))
    }
}

#[cfg(test)]
mod test {
    use super::CFFGenerationGraph;
    use crate::cff::cff_graph::{CFFEdgeType, VertexSet};
    use itertools::Itertools;

    #[test]
    fn test_graph_struct_triangle() {
        let triangle = vec![(0, 1), (1, 2), (2, 0)];
        let incoming_vertices = vec![
            (0, CFFEdgeType::External),
            (1, CFFEdgeType::External),
            (2, CFFEdgeType::External),
        ];

        let vertex_sets = [
            VertexSet::from_usize(0),
            VertexSet::from_usize(1),
            VertexSet::from_usize(2),
        ];

        let dummy_orientation = vec![true; 6];
        let cff_triangle =
            CFFGenerationGraph::from_vec(triangle, incoming_vertices, dummy_orientation);

        assert_eq!(cff_triangle.vertices.len(), 3);
        println!("node count test passed");

        assert!(cff_triangle.are_adjacent(&vertex_sets[0], &vertex_sets[1]));
        assert!(cff_triangle.are_adjacent(&vertex_sets[1], &vertex_sets[2]));
        assert!(cff_triangle.are_adjacent(&vertex_sets[2], &vertex_sets[0]));

        assert!(cff_triangle.are_adjacent(&vertex_sets[1], &vertex_sets[0]));
        assert!(cff_triangle.are_adjacent(&vertex_sets[2], &vertex_sets[1]));
        assert!(cff_triangle.are_adjacent(&vertex_sets[0], &vertex_sets[2]));

        println!("Adjacency test passed");

        assert!(cff_triangle.are_directed_adjacent(&vertex_sets[0], &vertex_sets[1]));
        assert!(cff_triangle.are_directed_adjacent(&vertex_sets[1], &vertex_sets[2]));
        assert!(cff_triangle.are_directed_adjacent(&vertex_sets[2], &vertex_sets[0]));

        assert!(!cff_triangle.are_directed_adjacent(&vertex_sets[1], &vertex_sets[0]));
        assert!(!cff_triangle.are_directed_adjacent(&vertex_sets[2], &vertex_sets[1]));
        assert!(!cff_triangle.are_directed_adjacent(&vertex_sets[0], &vertex_sets[2]));

        println!("Directed adjacency test passed");
    }

    #[test]
    fn test_graph_struct_with_virtext() {
        let line = vec![(0, 1)];
        let incoming_vertices = vec![
            (0, CFFEdgeType::VirtualExternal),
            (0, CFFEdgeType::External),
            (1, CFFEdgeType::VirtualExternal),
            (1, CFFEdgeType::External),
        ];

        let vertex_sets = [VertexSet::from_usize(0), VertexSet::from_usize(1)];

        let dummy_orientation = vec![true; 5];

        let cff_line = CFFGenerationGraph::from_vec(line, incoming_vertices, dummy_orientation);
        assert_eq!(cff_line.vertices.len(), 2);

        assert!(cff_line.are_adjacent(&vertex_sets[0], &vertex_sets[1]));
        assert!(cff_line.are_directed_adjacent(&vertex_sets[0], &vertex_sets[1]));
        assert!(!cff_line.are_directed_adjacent(&vertex_sets[1], &vertex_sets[0]));

        assert_eq!(cff_line.get_edges().len(), 5);
    }

    #[test]
    fn test_graph_struct_double_box() {
        let double_box = vec![(0, 1), (4, 5), (2, 3), (0, 4), (4, 2), (1, 5), (5, 3)];
        let incoming_vertices = (0..4).map(|i| (i, CFFEdgeType::External)).collect_vec();

        let dummy_orienatation = vec![true; double_box.len() + incoming_vertices.len()];

        let v = (0..=5).map(VertexSet::from_usize).collect::<Vec<_>>();

        let cff_double_box =
            CFFGenerationGraph::from_vec(double_box, incoming_vertices, dummy_orienatation);

        assert_eq!(cff_double_box.vertices.len(), 6);
        println!("node count test passed");

        assert!(cff_double_box.are_adjacent(&v[0], &v[1]));
        assert!(cff_double_box.are_adjacent(&v[1], &v[0]));
        assert!(cff_double_box.are_adjacent(&v[4], &v[5]));
        assert!(cff_double_box.are_adjacent(&v[5], &v[4]));
        assert!(cff_double_box.are_adjacent(&v[2], &v[3]));
        assert!(cff_double_box.are_adjacent(&v[3], &v[2]));
        assert!(cff_double_box.are_adjacent(&v[0], &v[4]));
        assert!(cff_double_box.are_adjacent(&v[4], &v[0]));
        assert!(cff_double_box.are_adjacent(&v[4], &v[2]));
        assert!(cff_double_box.are_adjacent(&v[2], &v[4]));
        assert!(cff_double_box.are_adjacent(&v[1], &v[5]));
        assert!(cff_double_box.are_adjacent(&v[5], &v[1]));
        assert!(cff_double_box.are_adjacent(&v[5], &v[3]));
        assert!(cff_double_box.are_adjacent(&v[3], &v[5]));

        assert!(!cff_double_box.are_adjacent(&v[0], &v[2]));
        assert!(!cff_double_box.are_adjacent(&v[0], &v[3]));
        assert!(!cff_double_box.are_adjacent(&v[0], &v[5]));
        assert!(!cff_double_box.are_adjacent(&v[1], &v[2]));
        assert!(!cff_double_box.are_adjacent(&v[1], &v[3]));
        assert!(!cff_double_box.are_adjacent(&v[1], &v[4]));
        // etc.

        println!("Adjacency test passed");

        assert!(cff_double_box.are_directed_adjacent(&v[0], &v[1]));
        assert!(cff_double_box.are_directed_adjacent(&v[4], &v[5]));
        assert!(cff_double_box.are_directed_adjacent(&v[2], &v[3]));
        assert!(cff_double_box.are_directed_adjacent(&v[0], &v[4]));
        assert!(cff_double_box.are_directed_adjacent(&v[4], &v[2]));
        assert!(cff_double_box.are_directed_adjacent(&v[1], &v[5]));
        assert!(cff_double_box.are_directed_adjacent(&v[5], &v[3]));

        assert!(!cff_double_box.are_directed_adjacent(&v[1], &v[0]));
        assert!(!cff_double_box.are_directed_adjacent(&v[5], &v[4]));
        assert!(!cff_double_box.are_directed_adjacent(&v[3], &v[2]));
        assert!(!cff_double_box.are_directed_adjacent(&v[4], &v[0]));
        assert!(!cff_double_box.are_directed_adjacent(&v[2], &v[4]));
        assert!(!cff_double_box.are_directed_adjacent(&v[5], &v[1]));
        assert!(!cff_double_box.are_directed_adjacent(&v[3], &v[5]));

        println!("Directed adjacency test passed");
    }

    #[test]
    fn test_vertex_that_is_not_v() {
        let triangle = vec![(0, 1), (1, 2), (2, 0)];
        let incoming_vertices = (0..3).map(|i| (i, CFFEdgeType::External)).collect_vec();

        let vertex_sets = [
            VertexSet::from_usize(0),
            VertexSet::from_usize(1),
            VertexSet::from_usize(2),
        ];

        let dummy_orientation = vec![true; triangle.len() + incoming_vertices.len()];

        let cff_triangle =
            CFFGenerationGraph::from_vec(triangle, incoming_vertices, dummy_orientation);
        let other_vertex = cff_triangle.get_vertex_that_is_not_v(&vertex_sets[0]);
        assert_ne!(other_vertex.vertex_set, vertex_sets[0]);
    }

    #[test]
    fn test_get_directed_neighbours() {
        let triangle = vec![(0, 1), (1, 2), (2, 0)];
        let incoming_vertices = (0..3).map(|i| (i, CFFEdgeType::External)).collect_vec();
        let vertex_sets = [
            VertexSet::from_usize(0),
            VertexSet::from_usize(1),
            VertexSet::from_usize(2),
        ];

        let dummy_orientation = vec![true; triangle.len() + incoming_vertices.len()];

        let cff_triangle =
            CFFGenerationGraph::from_vec(triangle, incoming_vertices, dummy_orientation);

        let neighbours = cff_triangle.get_directed_neighbours(&vertex_sets[0]);
        assert_eq!(neighbours.len(), 1);
        assert_eq!(neighbours[0].vertex_set, vertex_sets[1]);
    }

    #[test]

    fn test_has_directed_cycle() {
        let triangle = vec![(0, 1), (1, 2), (2, 0)];
        let incoming_vertices = (0..3).map(|i| (i, CFFEdgeType::External)).collect_vec();

        let vertex_sets = [
            VertexSet::from_usize(0),
            VertexSet::from_usize(1),
            VertexSet::from_usize(2),
        ];

        let dummy_orientation = vec![true; triangle.len() + incoming_vertices.len()];

        let cff_triangle =
            CFFGenerationGraph::from_vec(triangle, incoming_vertices, dummy_orientation.clone());
        assert!(cff_triangle.has_directed_cycle(&vertex_sets[0]));

        let triangle = vec![(0, 1), (1, 2), (0, 2)];
        let incoming_vertices = (0..3).map(|i| (i, CFFEdgeType::External)).collect_vec();

        let vertex_sets = [
            VertexSet::from_usize(0),
            VertexSet::from_usize(1),
            VertexSet::from_usize(2),
        ];

        let cff_triangle =
            CFFGenerationGraph::from_vec(triangle, incoming_vertices, dummy_orientation);
        assert!(!cff_triangle.has_directed_cycle(&vertex_sets[0]));
    }

    #[test]
    fn test_connected_complement() {
        let triangle = vec![(0, 1), (1, 2), (2, 0)];
        let incoming_vertices = (0..3).map(|i| (i, CFFEdgeType::External)).collect_vec();

        let vertex_sets = [
            VertexSet::from_usize(0),
            VertexSet::from_usize(1),
            VertexSet::from_usize(2),
        ];

        let dummy_orientation = vec![true; triangle.len() + incoming_vertices.len()];

        let cff_triangle =
            CFFGenerationGraph::from_vec(triangle, incoming_vertices, dummy_orientation);
        assert!(cff_triangle.has_connected_complement(&vertex_sets[0]));
        assert!(cff_triangle.has_connected_complement(&vertex_sets[1]));
        assert!(cff_triangle.has_connected_complement(&vertex_sets[2]));

        println!("triangle passed");

        let double_bubble = vec![(0, 1), (0, 1), (1, 2), (1, 2)];
        let incoming_vertices = [0, 2]
            .into_iter()
            .map(|i| (i, CFFEdgeType::External))
            .collect_vec();

        let vertex_sets = [
            VertexSet::from_usize(0),
            VertexSet::from_usize(1),
            VertexSet::from_usize(2),
        ];

        let dummy_orientation = vec![true; double_bubble.len() + incoming_vertices.len()];

        let cff_double_bubble =
            CFFGenerationGraph::from_vec(double_bubble, incoming_vertices, dummy_orientation);

        assert!(cff_double_bubble.has_connected_complement(&vertex_sets[0]));
        assert!(!cff_double_bubble.has_connected_complement(&vertex_sets[1]));
        assert!(cff_double_bubble.has_connected_complement(&vertex_sets[2]));
        println!("double bubble passed");

        let single_bubble = vec![(0, 1), (0, 1)];
        let incoming_vertices = (0..2).map(|i| (i, CFFEdgeType::External)).collect_vec();

        let vertex_sets = [VertexSet::from_usize(0), VertexSet::from_usize(1)];

        let dummy_orientation = vec![true; single_bubble.len() + incoming_vertices.len()];

        let cff_single_bubble =
            CFFGenerationGraph::from_vec(single_bubble, incoming_vertices, dummy_orientation);

        assert!(cff_single_bubble.has_connected_complement(&vertex_sets[0]));
        assert!(cff_single_bubble.has_connected_complement(&vertex_sets[1]));

        println!("single bubble passed");
    }

    #[test]
    fn test_get_vertex_type() {
        let triangle = vec![(0, 1), (2, 1), (0, 2)];
        let incoming_vertices = (0..3).map(|i| (i, CFFEdgeType::External)).collect_vec();
        let vertex_sets = [
            VertexSet::from_usize(0),
            VertexSet::from_usize(1),
            VertexSet::from_usize(2),
        ];

        let dummy_orientation = vec![true; triangle.len() + incoming_vertices.len()];

        let cff_triangle =
            CFFGenerationGraph::from_vec(triangle, incoming_vertices, dummy_orientation);

        assert_eq!(
            cff_triangle.get_vertex_type(&vertex_sets[0]),
            super::VertexType::Source
        );
        assert_eq!(
            cff_triangle.get_vertex_type(&vertex_sets[1]),
            super::VertexType::Sink
        );
        assert_eq!(
            cff_triangle.get_vertex_type(&vertex_sets[2]),
            super::VertexType::None
        );
    }

    #[test]
    fn test_source_sink_candidate_list() {
        let triangle = vec![(0, 1), (2, 1), (0, 2)];
        let incoming_vertices = (0..3).map(|i| (i, CFFEdgeType::External)).collect_vec();
        let vertex_set = [
            VertexSet::from_usize(0),
            VertexSet::from_usize(1),
            VertexSet::from_usize(2),
        ];

        let dummy_orientation = vec![true; triangle.len() + incoming_vertices.len()];

        let cff_triangle =
            CFFGenerationGraph::from_vec(triangle, incoming_vertices, dummy_orientation);

        let source_sink_candidates = cff_triangle.get_source_sink_candidate_list();
        assert_eq!(source_sink_candidates.len(), 2);

        let source_sink_candidates = source_sink_candidates
            .iter()
            .map(|vertex| vertex.vertex_set)
            .collect_vec();

        assert!(source_sink_candidates.contains(&vertex_set[0]));
        assert!(source_sink_candidates.contains(&vertex_set[1]));

        println!("triangle passed");

        let bubble = vec![(0, 1), (0, 1)];
        let incoming_vertices = (0..2).map(|i| (i, CFFEdgeType::External)).collect_vec();
        let vertex_set = [VertexSet::from_usize(0), VertexSet::from_usize(1)];

        let dummy_orientation = vec![true; bubble.len() + incoming_vertices.len()];

        let cff_bubble = CFFGenerationGraph::from_vec(bubble, incoming_vertices, dummy_orientation);

        let source_sink_candidates = cff_bubble.get_source_sink_candidate_list();
        let source_sink_candidates = source_sink_candidates
            .iter()
            .map(|vertex| vertex.vertex_set)
            .collect_vec();

        assert_eq!(source_sink_candidates.len(), 2);
        assert!(source_sink_candidates.contains(&vertex_set[0]));
        assert!(source_sink_candidates.contains(&vertex_set[1]));

        println!("bubble passed");
    }

    #[test]
    fn test_contract_from_vertex() {
        let triangle = vec![(0, 1), (2, 1), (0, 2)];
        let incoming_vertices = (0..3).map(|i| (i, CFFEdgeType::External)).collect_vec();

        let vertex_set = [
            VertexSet::from_usize(0),
            VertexSet::from_usize(1),
            VertexSet::from_usize(2),
        ];

        let joined_vertex = vertex_set[0].join(&vertex_set[1]);
        let dummy_orientation = vec![true; triangle.len() + incoming_vertices.len()];

        let cff_triangle =
            CFFGenerationGraph::from_vec(triangle, incoming_vertices, dummy_orientation);

        let contracted_graph = cff_triangle.contract_vertices(&vertex_set[0], &vertex_set[1]);

        assert_eq!(contracted_graph.vertices.len(), 2);

        let contracted_vertex = contracted_graph.get_vertex(&joined_vertex);
        assert_eq!(contracted_vertex.incoming_edges.len(), 3);
        assert_eq!(contracted_vertex.outgoing_edges.len(), 1);

        let box_like = vec![(1, 0), (3, 0), (1, 2), (2, 3)];
        let incoming_vertices = vec![
            (0, CFFEdgeType::VirtualExternal),
            (1, CFFEdgeType::VirtualExternal),
            (2, CFFEdgeType::External),
            (3, CFFEdgeType::External),
        ];

        let vertex_set = [
            VertexSet::from_usize(0),
            VertexSet::from_usize(1),
            VertexSet::from_usize(2),
            VertexSet::from_usize(3),
        ];

        let dummy_orientation = vec![true; box_like.len() + incoming_vertices.len()];
        let cff_box = CFFGenerationGraph::from_vec(box_like, incoming_vertices, dummy_orientation);

        let joined_vertex = vertex_set[0].join(&vertex_set[1]);
        let contracted_box = cff_box.contract_vertices(&vertex_set[0], &vertex_set[1]);

        let contracted_vertex = contracted_box.get_vertex(&joined_vertex);

        assert_eq!(contracted_vertex.incoming_edges.len(), 3);
        assert_eq!(contracted_vertex.outgoing_edges.len(), 1);
    }

    #[test]
    fn test_get_edges() {
        let triangle = vec![(0, 1), (2, 1), (0, 2)];
        let incoming_vertices = (0..3).map(|i| (i, CFFEdgeType::External)).collect_vec();

        let dummy_orientation = vec![true; triangle.len() + incoming_vertices.len()];

        let cff_triangle =
            CFFGenerationGraph::from_vec(triangle, incoming_vertices, dummy_orientation);

        let edges = cff_triangle.get_edges();
        assert_eq!(edges, vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_has_edge() {
        let triangle = vec![(0, 1), (2, 1), (0, 2)];
        let incoming_vertices = (0..3).map(|i| (i, CFFEdgeType::External)).collect_vec();

        let dummy_orientation = vec![true; triangle.len() + incoming_vertices.len()];

        let cff_triangle =
            CFFGenerationGraph::from_vec(triangle, incoming_vertices, dummy_orientation);

        for edge in 0..6 {
            assert!(cff_triangle.has_edge(edge))
        }

        for edge in 6..12 {
            assert!(!cff_triangle.has_edge(edge))
        }
    }

    #[test]
    fn test_generate_cut() {
        let box_edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];
        let incoming_vertices = (0..4).map(|i| (i, CFFEdgeType::External)).collect_vec();
        let incoming_vertices_len = incoming_vertices.len();

        let dummy_orientation = vec![true; box_edges.len() + incoming_vertices.len()];

        let cff_box = CFFGenerationGraph::from_vec(box_edges, incoming_vertices, dummy_orientation);

        let vertex_sets = [
            VertexSet::from_usize(0),
            VertexSet::from_usize(1),
            VertexSet::from_usize(2),
            VertexSet::from_usize(3),
        ];

        let circled = vertex_sets[0].join(&vertex_sets[1]);

        let (left_cut, right_cut) = cff_box.generate_cut(circled);

        assert_eq!(left_cut.vertices.len(), 2);
        assert_eq!(right_cut.vertices.len(), 2);

        assert!(left_cut.has_edge(incoming_vertices_len));
        assert!(left_cut.has_edge(1 + incoming_vertices_len));
        assert!(left_cut.has_edge(3 + incoming_vertices_len));
        assert!(!left_cut.has_edge(2 + incoming_vertices_len));

        assert!(right_cut.has_edge(1 + incoming_vertices_len));
        assert!(right_cut.has_edge(2 + incoming_vertices_len));
        assert!(right_cut.has_edge(3 + incoming_vertices_len));
        assert!(!right_cut.has_edge(incoming_vertices_len));

        #[allow(clippy::if_same_then_else)]
        for edge in left_cut.iter_all_edges() {
            if edge.edge_id == incoming_vertices_len {
                assert_eq!(edge.edge_type, CFFEdgeType::Virtual)
            } else if edge.edge_id == 1 + incoming_vertices_len {
                assert_eq!(edge.edge_type, CFFEdgeType::VirtualExternal)
            } else if edge.edge_id == 3 + incoming_vertices_len {
                assert_eq!(edge.edge_type, CFFEdgeType::VirtualExternal)
            } else {
                assert_eq!(edge.edge_type, CFFEdgeType::External)
            }
        }

        #[allow(clippy::if_same_then_else)]
        for edge in right_cut.iter_all_edges() {
            if edge.edge_id == 2 + incoming_vertices_len {
                assert_eq!(edge.edge_type, CFFEdgeType::Virtual)
            } else if edge.edge_id == 1 + incoming_vertices_len {
                assert_eq!(edge.edge_type, CFFEdgeType::VirtualExternal)
            } else if edge.edge_id == 3 + incoming_vertices_len {
                assert_eq!(edge.edge_type, CFFEdgeType::VirtualExternal)
            } else {
                assert_eq!(edge.edge_type, CFFEdgeType::External)
            }
        }

        //println!("left {:#?}", left_cut);
        //println!("right {:#?}", right_cut);
    }
}
