use crate::cff::hsurface::Hsurface;
use ahash::{HashMap, HashSet, HashSetExt};
use bincode::Encode;
use bincode_trait_derive::Decode;
use bitvec::vec::BitVec;
use color_eyre::Result;
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::{
    hedgevec::EdgeVec,
    involution::{EdgeIndex, Flow, HedgePair, Orientation},
    subgraph::{ModifySubgraph, SubGraph},
    HedgeGraph, NodeIndex,
};
use serde::{Deserialize, Serialize};
use std::hash::Hash;

use super::{
    esurface::{Esurface, ExternalShift},
    surface::{HybridSurface, UnitSurface},
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

    fn contract(&self, other: &Self, remove_single_edge: Option<EdgeIndex>) -> Self {
        let new_vertex_set = self.vertex_set.join(&other.vertex_set);

        let incoming_edges_of_new = self
            .incoming_edges
            .iter()
            .filter(|edge| {
                if let Some(edge_to_be_removed) = remove_single_edge {
                    edge.edge_id != edge_to_be_removed
                } else {
                    !other.outgoing_edges.contains(edge)
                }
            })
            .chain(other.incoming_edges.iter().filter(|edge| {
                if let Some(edge_to_be_removed) = remove_single_edge {
                    edge.edge_id != edge_to_be_removed
                } else {
                    !self.outgoing_edges.contains(edge)
                }
            }))
            .copied()
            .sorted_by(|edge_1, edge_2| edge_1.edge_id.cmp(&edge_2.edge_id))
            .collect_vec();

        let outgoing_edges_of_new = self
            .outgoing_edges
            .iter()
            .filter(|edge| {
                if let Some(edge_to_be_removed) = remove_single_edge {
                    edge.edge_id != edge_to_be_removed
                } else {
                    !other.incoming_edges.contains(edge)
                }
            })
            .chain(other.outgoing_edges.iter().filter(|edge| {
                if let Some(edge_to_be_removed) = remove_single_edge {
                    edge.edge_id != edge_to_be_removed
                } else {
                    !self.incoming_edges.contains(edge)
                }
            }))
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

    fn has_impossible_edge(&self) -> bool {
        let mut dedup_outgoing = self.outgoing_edges.clone();
        dedup_outgoing.dedup();
        let mut dedup_incoming = self.incoming_edges.clone();
        dedup_incoming.dedup();

        dedup_incoming.len() != self.incoming_edges.len()
            || dedup_outgoing.len() != self.outgoing_edges.len()
    }

    fn has_edge(&self, edge_id: EdgeIndex) -> bool {
        self.iter_all_edges().any(|edge| edge.edge_id == edge_id)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize, Encode, Decode)]
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

    pub fn join(&self, other: &VertexSet) -> VertexSet {
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

    fn get_nodes(&self) -> Vec<NodeIndex> {
        (0..MAX_VERTEX_COUNT)
            .filter(|id| self.vertex_set & (1 << id) != 0)
            .map(NodeIndex::from)
            .collect()
    }

    pub fn subgraph<E, V>(&self, graph: &HedgeGraph<E, V>) -> BitVec {
        let mut result: BitVec = graph.empty_subgraph();
        for hedge in self
            .get_nodes()
            .iter()
            .flat_map(|node_id| graph.iter_crown(*node_id))
        {
            result.add(hedge);
        }
        result
    }

    pub fn dummy() -> Self {
        VertexSet { vertex_set: 0 }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
struct CFFEdge {
    edge_id: EdgeIndex,
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
    pub global_orientation: EdgeVec<Orientation>,
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
    fn has_impossible_edge(&self) -> bool {
        self.vertices
            .iter()
            .any(|vertex| vertex.has_impossible_edge())
    }

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
        orientation: Option<EdgeVec<Orientation>>,
    ) -> Self {
        use crate::utils;

        let total_num_edges = edges.len() + incoming_vertices.len();

        let edges = edges
            .into_iter()
            .map(|(from, to)| (VertexSet::from_usize(from), VertexSet::from_usize(to)))
            .collect_vec();

        let mut unique_vertex_ids = HashSet::new();

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
                edge_id: EdgeIndex::from(edge_id),
                edge_type: *edge_type,
            };
            cff_vertex.incoming_edges.push(incoming_edge);
        }

        for (edge_id, (left, right)) in edges.iter().enumerate() {
            let edge_id = edge_id + incoming_vertices.len();
            let cff_edge = CFFEdge {
                edge_id: EdgeIndex::from(edge_id),
                edge_type: CFFEdgeType::Virtual,
            };

            let left_vertex = unique_vertices.get_mut(left).unwrap();
            left_vertex.outgoing_edges.push(cff_edge);

            let right_vertex = unique_vertices.get_mut(right).unwrap();
            right_vertex.incoming_edges.push(cff_edge);
        }

        let nodes = unique_vertices.into_values().collect_vec();
        let global_orientation = match orientation {
            Some(orientation) => orientation,
            None => utils::dummy_hedge_graph(total_num_edges)
                .new_edgevec_from_iter(vec![Orientation::Default; total_num_edges])
                .unwrap(),
        };

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

    fn remove_edge(&mut self, edge_id: EdgeIndex) {
        for vertex in self.vertices.iter_mut() {
            vertex.incoming_edges.retain(|edge| edge.edge_id != edge_id);
            vertex.outgoing_edges.retain(|edge| edge.edge_id != edge_id);
        }
    }

    pub fn remove_self_edges(&mut self) {
        let self_edges = self.get_self_edges();

        for self_edge in self_edges.iter() {
            self.remove_edge(*self_edge);
        }
    }

    pub fn get_self_edges(&self) -> Vec<EdgeIndex> {
        let mut self_edges = vec![];

        for vertex in self.vertices.iter() {
            for edge in vertex.incoming_edges.iter() {
                if vertex.outgoing_edges.contains(edge) {
                    self_edges.push(edge.edge_id);
                }
            }
        }

        self_edges
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
        self.contract_vertices_impl(vertex_1, vertex_2, None)
    }

    pub fn contract_edge(&self, edge_id: EdgeIndex) -> Self {
        let (source, sink) = self.get_source_sink_of_edge(edge_id);
        let vertex_1 = &source.vertex_set;
        let vertex_2 = &sink.vertex_set;

        // self edges need a special treatment
        if vertex_1 == vertex_2 {
            let mut new_graph = self.clone();
            new_graph.remove_edge(edge_id);
            return new_graph;
        }

        self.contract_vertices_impl(vertex_1, vertex_2, Some(edge_id))
    }

    fn contract_vertices_impl(
        &self,
        vertex_1: &VertexSet,
        vertex_2: &VertexSet,
        remove_single_edge: Option<EdgeIndex>,
    ) -> Self {
        let vertex_1 = self.get_vertex(vertex_1);
        let vertex_2 = self.get_vertex(vertex_2);

        let new_vertex = vertex_1.contract(vertex_2, remove_single_edge);

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

    fn get_source_sink_of_edge(&self, edge_id: EdgeIndex) -> (&CFFVertex, &CFFVertex) {
        let source = self
            .vertices
            .iter()
            .find(|vertex| {
                vertex
                    .outgoing_edges
                    .iter()
                    .any(|edge| edge.edge_id == edge_id)
            })
            .expect("not a virtual edge");

        let sink = self
            .vertices
            .iter()
            .find(|vertex| {
                vertex
                    .incoming_edges
                    .iter()
                    .any(|edge| edge.edge_id == edge_id)
            })
            .expect("not a virtual edge");

        (source, sink)
    }

    pub fn generate_children(&self) -> (Option<Vec<Self>>, HybridSurface) {
        if self.vertices.len() < 2 {
            return (None, HybridSurface::Unit(UnitSurface {}));
        }

        let vertex = if let Some(vertex) = self.get_source_sink_greedy() {
            vertex
        } else {
            panic!(
                "could not find vertex to contract from for graph {:#?}",
                self
            );
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

            let esurface = Esurface {
                energies: positive_energies,
                external_shift,
                vertex_set: vertex.vertex_set,
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

    /// for now only non-cut graphs are supported
    pub fn new<E, V, H>(
        graph: &HedgeGraph<E, V, H>,
        global_orientation: EdgeVec<Orientation>,
    ) -> Self {
        let mut vertices = (0..graph.n_nodes()).map(CFFVertex::new).collect_vec();

        for (hedge_pair, edge_id, _) in graph.iter_edges() {
            match hedge_pair {
                HedgePair::Unpaired { hedge, flow } => {
                    let vertex = Into::<usize>::into(graph.node_id(hedge));
                    let edge_type = CFFEdgeType::External;
                    let edge = CFFEdge { edge_id, edge_type };
                    match flow {
                        Flow::Source => {
                            vertices[vertex].outgoing_edges.push(edge);
                        }
                        Flow::Sink => {
                            vertices[vertex].incoming_edges.push(edge);
                        }
                    }
                }
                HedgePair::Paired { source, sink } => {
                    let source_vertex = Into::<usize>::into(graph.node_id(source));
                    let sink_vertex = Into::<usize>::into(graph.node_id(sink));
                    let edge_type = CFFEdgeType::Virtual;
                    let edge = CFFEdge { edge_id, edge_type };

                    match global_orientation[edge_id] {
                        Orientation::Default => {
                            vertices[source_vertex].outgoing_edges.push(edge);
                            vertices[sink_vertex].incoming_edges.push(edge);
                        }
                        Orientation::Reversed => {
                            vertices[source_vertex].incoming_edges.push(edge);
                            vertices[sink_vertex].outgoing_edges.push(edge);
                        }
                        Orientation::Undirected => {
                            panic!("Can not generate CFF with undirected edges")
                        }
                    }
                }
                HedgePair::Split { .. } => unreachable!(),
            }
        }

        Self {
            vertices,
            global_orientation,
        }
    }

    pub fn new_from_subgraph<E, V, H, S: SubGraph>(
        graph: &HedgeGraph<E, V, H>,
        global_orientation: EdgeVec<Orientation>,
        subgraph: &S,
    ) -> Result<Self> {
        let mut vertices = HashMap::default();

        for (node_id, _, _data) in graph.iter_nodes_of(subgraph) {
            let vertex = CFFVertex::new(node_id.into());
            vertices.insert(node_id, vertex);
        }

        for (hedge_pair, edge_index, _data) in graph.iter_edges_of(subgraph) {
            match hedge_pair {
                HedgePair::Unpaired { hedge, flow } => {
                    let vertex = graph.node_id(hedge);
                    let edge_type = CFFEdgeType::External;
                    let edge = CFFEdge {
                        edge_id: edge_index,
                        edge_type,
                    };
                    match flow {
                        Flow::Source => {
                            vertices.get_mut(&vertex).unwrap().outgoing_edges.push(edge);
                        }
                        Flow::Sink => {
                            vertices.get_mut(&vertex).unwrap().incoming_edges.push(edge);
                        }
                    }
                }
                HedgePair::Paired { source, sink } => {
                    let source_vertex = graph.node_id(source);
                    let sink_vertex = graph.node_id(sink);
                    let edge_type = CFFEdgeType::Virtual;
                    let edge = CFFEdge {
                        edge_id: edge_index,
                        edge_type,
                    };

                    match global_orientation[edge_index] {
                        Orientation::Default => {
                            vertices
                                .get_mut(&source_vertex)
                                .unwrap()
                                .outgoing_edges
                                .push(edge);
                            vertices
                                .get_mut(&sink_vertex)
                                .unwrap()
                                .incoming_edges
                                .push(edge);
                        }
                        Orientation::Reversed => {
                            vertices
                                .get_mut(&source_vertex)
                                .unwrap()
                                .incoming_edges
                                .push(edge);
                            vertices
                                .get_mut(&sink_vertex)
                                .unwrap()
                                .outgoing_edges
                                .push(edge);
                        }
                        Orientation::Undirected => {
                            return Err(eyre!(
                                "undirected edge found, edge_id: {}, subgraph: \n {}",
                                edge_index,
                                graph.dot(subgraph)
                            ));
                        }
                    }
                }
                HedgePair::Split {
                    source,
                    sink,
                    split,
                } => {
                    let edge_type = CFFEdgeType::VirtualExternal;
                    let edge = CFFEdge {
                        edge_id: edge_index,
                        edge_type,
                    };
                    match split {
                        Flow::Source => {
                            let vertex = graph.node_id(source);
                            match global_orientation[edge_index] {
                                Orientation::Default => {
                                    vertices.get_mut(&vertex).unwrap().outgoing_edges.push(edge)
                                }
                                Orientation::Reversed => {
                                    vertices.get_mut(&vertex).unwrap().incoming_edges.push(edge)
                                }
                                Orientation::Undirected => {
                                    return Err(eyre!(
                                        "undirected edge found for source split, edge_id: {}, subgraph: \n {}",
                                        edge_index,
                                        graph.dot(subgraph)
                                    ))
                                }
                            }
                        }
                        Flow::Sink => {
                            let vertex = graph.node_id(sink);
                            match global_orientation[edge_index] {
                                Orientation::Default => {
                                    vertices.get_mut(&vertex).unwrap().incoming_edges.push(edge)
                                }
                                Orientation::Reversed => {
                                    vertices.get_mut(&vertex).unwrap().outgoing_edges.push(edge)
                                }
                                Orientation::Undirected => {
                                    return Err(eyre!(
                                        "undirected edge found for sink split, edge_id: {}, subgraph: \n {}",
                                        edge_index,
                                        graph.dot(subgraph)
                                    ))
                                }
                            }
                        }
                    }
                }
            }
        }

        let vertices = vertices
            .into_values()
            .sorted_by(|a, b| a.vertex_set.vertex_set.cmp(&b.vertex_set.vertex_set))
            .collect_vec();

        let res = Self {
            vertices,
            global_orientation,
        };

        Ok(res)
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

    fn get_edges(&self) -> Vec<EdgeIndex> {
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

    fn has_edge(&self, edge_id: EdgeIndex) -> bool {
        self.vertices.iter().any(|vertex| vertex.has_edge(edge_id))
    }
}

#[cfg(test)]
mod test {
    use super::CFFGenerationGraph;
    use crate::cff::cff_graph::{CFFEdge, CFFEdgeType, VertexSet};
    use bitvec::vec::BitVec;
    use itertools::Itertools;
    use linnet::half_edge::{
        builder::HedgeGraphBuilder,
        involution::{EdgeIndex, Flow, Orientation},
        nodestore::NodeStorageVec,
        subgraph::SubGraphOps,
        HedgeGraph,
    };

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

        let cff_triangle = CFFGenerationGraph::from_vec(triangle, incoming_vertices, None);

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

        let cff_line = CFFGenerationGraph::from_vec(line, incoming_vertices, None);
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

        let v = (0..=5).map(VertexSet::from_usize).collect::<Vec<_>>();

        let cff_double_box = CFFGenerationGraph::from_vec(double_box, incoming_vertices, None);

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

        let cff_triangle = CFFGenerationGraph::from_vec(triangle, incoming_vertices, None);
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

        let cff_triangle = CFFGenerationGraph::from_vec(triangle, incoming_vertices, None);

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

        let cff_triangle = CFFGenerationGraph::from_vec(triangle, incoming_vertices, None);
        assert!(cff_triangle.has_directed_cycle(&vertex_sets[0]));

        let triangle = vec![(0, 1), (1, 2), (0, 2)];
        let incoming_vertices = (0..3).map(|i| (i, CFFEdgeType::External)).collect_vec();

        let vertex_sets = [
            VertexSet::from_usize(0),
            VertexSet::from_usize(1),
            VertexSet::from_usize(2),
        ];

        let cff_triangle = CFFGenerationGraph::from_vec(triangle, incoming_vertices, None);
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

        let cff_triangle = CFFGenerationGraph::from_vec(triangle, incoming_vertices, None);
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

        let cff_double_bubble =
            CFFGenerationGraph::from_vec(double_bubble, incoming_vertices, None);

        assert!(cff_double_bubble.has_connected_complement(&vertex_sets[0]));
        assert!(!cff_double_bubble.has_connected_complement(&vertex_sets[1]));
        assert!(cff_double_bubble.has_connected_complement(&vertex_sets[2]));
        println!("double bubble passed");

        let single_bubble = vec![(0, 1), (0, 1)];
        let incoming_vertices = (0..2).map(|i| (i, CFFEdgeType::External)).collect_vec();

        let vertex_sets = [VertexSet::from_usize(0), VertexSet::from_usize(1)];

        let cff_single_bubble =
            CFFGenerationGraph::from_vec(single_bubble, incoming_vertices, None);

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

        let cff_triangle = CFFGenerationGraph::from_vec(triangle, incoming_vertices, None);

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

        let cff_triangle = CFFGenerationGraph::from_vec(triangle, incoming_vertices, None);

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

        let cff_bubble = CFFGenerationGraph::from_vec(bubble, incoming_vertices, None);

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

        let cff_triangle = CFFGenerationGraph::from_vec(triangle, incoming_vertices, None);

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

        let cff_box = CFFGenerationGraph::from_vec(box_like, incoming_vertices, None);

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

        let cff_triangle = CFFGenerationGraph::from_vec(triangle, incoming_vertices, None);

        let edges = cff_triangle.get_edges();
        let comp_edges = (0..6).map(EdgeIndex::from).collect::<Vec<EdgeIndex>>();

        //assert_eq!(edges, vec![0, 1, 2, 3, 4, 5]);
        assert_eq!(edges, comp_edges);
    }

    #[test]
    fn test_has_edge() {
        let triangle = vec![(0, 1), (2, 1), (0, 2)];
        let incoming_vertices = (0..3).map(|i| (i, CFFEdgeType::External)).collect_vec();

        let cff_triangle = CFFGenerationGraph::from_vec(triangle, incoming_vertices, None);

        for edge in (0..6).map(EdgeIndex::from) {
            assert!(cff_triangle.has_edge(edge))
        }

        for edge in (6..12).map(EdgeIndex::from) {
            assert!(!cff_triangle.has_edge(edge))
        }
    }

    #[test]
    fn test_generate_cut() {
        let box_edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];
        let incoming_vertices = (0..4).map(|i| (i, CFFEdgeType::External)).collect_vec();
        let incoming_vertices_len = incoming_vertices.len();

        let cff_box = CFFGenerationGraph::from_vec(box_edges, incoming_vertices, None);

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

        assert!(left_cut.has_edge(EdgeIndex::from(incoming_vertices_len)));
        assert!(left_cut.has_edge(EdgeIndex::from(1 + incoming_vertices_len)));
        assert!(left_cut.has_edge(EdgeIndex::from(3 + incoming_vertices_len)));
        assert!(!left_cut.has_edge(EdgeIndex::from(2 + incoming_vertices_len)));

        assert!(right_cut.has_edge(EdgeIndex::from(1 + incoming_vertices_len)));
        assert!(right_cut.has_edge(EdgeIndex::from(2 + incoming_vertices_len)));
        assert!(right_cut.has_edge(EdgeIndex::from(3 + incoming_vertices_len)));
        assert!(!right_cut.has_edge(EdgeIndex::from(incoming_vertices_len)));

        #[allow(clippy::if_same_then_else)]
        for edge in left_cut.iter_all_edges() {
            if edge.edge_id == incoming_vertices_len.into() {
                assert_eq!(edge.edge_type, CFFEdgeType::Virtual)
            } else if edge.edge_id == (1 + incoming_vertices_len).into() {
                assert_eq!(edge.edge_type, CFFEdgeType::VirtualExternal)
            } else if edge.edge_id == (3 + incoming_vertices_len).into() {
                assert_eq!(edge.edge_type, CFFEdgeType::VirtualExternal)
            } else {
                assert_eq!(edge.edge_type, CFFEdgeType::External)
            }
        }

        #[allow(clippy::if_same_then_else)]
        for edge in right_cut.iter_all_edges() {
            if edge.edge_id == (2 + incoming_vertices_len).into() {
                assert_eq!(edge.edge_type, CFFEdgeType::Virtual)
            } else if edge.edge_id == (1 + incoming_vertices_len).into() {
                assert_eq!(edge.edge_type, CFFEdgeType::VirtualExternal)
            } else if edge.edge_id == (3 + incoming_vertices_len).into() {
                assert_eq!(edge.edge_type, CFFEdgeType::VirtualExternal)
            } else {
                assert_eq!(edge.edge_type, CFFEdgeType::External)
            }
        }

        //println!("left {:#?}", left_cut);
        //println!("right {:#?}", right_cut);
    }

    #[test]
    fn test_contract_edge() {
        let mut hedge_graph_builder = HedgeGraphBuilder::new();
        let nodes = (0..2)
            .map(|_| hedge_graph_builder.add_node(()))
            .collect_vec();

        hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);

        hedge_graph_builder.add_external_edge(nodes[0], (), Orientation::Undirected, Flow::Sink);
        hedge_graph_builder.add_external_edge(nodes[1], (), Orientation::Undirected, Flow::Source);

        let hedge_graph: HedgeGraph<_, _, ()> = hedge_graph_builder.build::<NodeStorageVec<_>>();
        let global_orientation = hedge_graph.new_edgevec(|_, _, _| Orientation::Default);

        let cff_graph = CFFGenerationGraph::new(&hedge_graph, global_orientation);
        let contracted = cff_graph.contract_edge(EdgeIndex::from(0));

        assert!(!contracted.has_edge(EdgeIndex::from(0)));
        assert!(contracted.has_edge(EdgeIndex::from(1)));
        assert!(contracted.has_edge(EdgeIndex::from(2)));
        assert!(contracted.has_edge(EdgeIndex::from(3)));
        assert!(contracted.has_edge(EdgeIndex::from(4)));

        assert_eq!(contracted.vertices.len(), 1);

        let vertex = contracted.vertices[0].clone();
        assert_eq!(
            vertex.incoming_edges,
            vec![
                CFFEdge {
                    edge_id: EdgeIndex::from(1),
                    edge_type: CFFEdgeType::Virtual,
                },
                CFFEdge {
                    edge_id: EdgeIndex::from(2),
                    edge_type: CFFEdgeType::Virtual,
                },
                CFFEdge {
                    edge_id: EdgeIndex::from(3),
                    edge_type: CFFEdgeType::External
                }
            ]
        );

        assert_eq!(
            vertex.outgoing_edges,
            vec![
                CFFEdge {
                    edge_id: EdgeIndex::from(1),
                    edge_type: CFFEdgeType::Virtual,
                },
                CFFEdge {
                    edge_id: EdgeIndex::from(2),
                    edge_type: CFFEdgeType::Virtual,
                },
                CFFEdge {
                    edge_id: EdgeIndex::from(4),
                    edge_type: CFFEdgeType::External
                }
            ]
        );

        let mut tri_box_builder = HedgeGraphBuilder::new();

        let nodes = (0..5).map(|_| tri_box_builder.add_node(())).collect_vec();

        tri_box_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
        tri_box_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);
        tri_box_builder.add_edge(nodes[0], nodes[2], (), Orientation::Undirected);

        tri_box_builder.add_edge(nodes[1], nodes[3], (), Orientation::Undirected);
        tri_box_builder.add_edge(nodes[2], nodes[4], (), Orientation::Undirected);
        tri_box_builder.add_edge(nodes[3], nodes[4], (), Orientation::Undirected);

        tri_box_builder.add_external_edge(nodes[0], (), Orientation::Undirected, Flow::Sink);
        tri_box_builder.add_external_edge(nodes[3], (), Orientation::Undirected, Flow::Source);
        tri_box_builder.add_external_edge(nodes[4], (), Orientation::Undirected, Flow::Source);

        let tri_box: HedgeGraph<(), (), ()> = tri_box_builder.build::<NodeStorageVec<_>>();
        let global_orientation = tri_box.new_edgevec(|_, _, _| Orientation::Default);
        let mut tri_box_cff_graph = CFFGenerationGraph::new(&tri_box, global_orientation);

        assert!(!tri_box_cff_graph.has_impossible_edge());
        tri_box_cff_graph = tri_box_cff_graph.contract_edge(EdgeIndex::from(0));
        assert!(!tri_box_cff_graph.has_impossible_edge());
        tri_box_cff_graph = tri_box_cff_graph.contract_edge(EdgeIndex::from(1));
        assert!(!tri_box_cff_graph.has_impossible_edge());
        tri_box_cff_graph = tri_box_cff_graph.contract_edge(EdgeIndex::from(2));
        assert!(!tri_box_cff_graph.has_impossible_edge());

        tri_box_cff_graph.remove_self_edges();
        assert!(!tri_box_cff_graph.has_impossible_edge());

        for vertex in tri_box_cff_graph.vertices.iter() {
            let all_edges = vertex.iter_all_edges().collect_vec();
            assert_eq!(all_edges.len(), 3);
            let num_external = all_edges
                .iter()
                .filter(|edge| edge.edge_type == CFFEdgeType::External)
                .count();
            assert_eq!(num_external, 1);

            let num_virtual = all_edges
                .iter()
                .filter(|edge| edge.edge_type == CFFEdgeType::Virtual)
                .count();
            assert_eq!(num_virtual, 2);

            let num_virtual_external = all_edges
                .iter()
                .filter(|edge| edge.edge_type == CFFEdgeType::VirtualExternal)
                .count();
            assert_eq!(num_virtual_external, 0);
        }
    }

    #[test]
    fn test_new_from_subgraph() {
        let mut hedge_graph_builder = HedgeGraphBuilder::new();

        let nodes = (0..4)
            .map(|_| hedge_graph_builder.add_node(()))
            .collect_vec();

        hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[0], nodes[2], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[1], nodes[3], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[2], nodes[3], (), Orientation::Undirected);

        hedge_graph_builder.add_external_edge(nodes[0], (), Orientation::Undirected, Flow::Sink);
        hedge_graph_builder.add_external_edge(nodes[3], (), Orientation::Undirected, Flow::Source);

        let hedge_graph: HedgeGraph<(), (), ()> = hedge_graph_builder.build::<NodeStorageVec<_>>();
        let node_0: BitVec = hedge_graph.iter_crown(nodes[0]).into();
        let node_1: BitVec = hedge_graph.iter_crown(nodes[1]).into();
        let node_2: BitVec = hedge_graph.iter_crown(nodes[2]).into();

        let left_triangle = node_0.union(&node_1).union(&node_2);

        let global_orientation = hedge_graph.new_edgevec(|_, _, _| Orientation::Default);

        let cff_graph =
            CFFGenerationGraph::new_from_subgraph(&hedge_graph, global_orientation, &left_triangle)
                .unwrap();

        assert!(cff_graph.has_edge(EdgeIndex::from(0)));
        assert!(cff_graph.has_edge(EdgeIndex::from(1)));
        assert!(cff_graph.has_edge(EdgeIndex::from(2)));
        assert!(cff_graph.has_edge(EdgeIndex::from(3)));
        assert!(cff_graph.has_edge(EdgeIndex::from(4)));
        assert!(cff_graph.has_edge(EdgeIndex::from(5)));
        assert!(!cff_graph.has_edge(EdgeIndex::from(6)));

        for edge in cff_graph.iter_all_edges() {
            match edge.edge_id.into() {
                0 => assert_eq!(edge.edge_type, CFFEdgeType::Virtual),
                1 => assert_eq!(edge.edge_type, CFFEdgeType::Virtual),
                2 => assert_eq!(edge.edge_type, CFFEdgeType::Virtual),
                3 => assert_eq!(edge.edge_type, CFFEdgeType::VirtualExternal),
                4 => assert_eq!(edge.edge_type, CFFEdgeType::VirtualExternal),
                5 => assert_eq!(edge.edge_type, CFFEdgeType::External),
                _ => unreachable!(),
            }
        }
    }
}
