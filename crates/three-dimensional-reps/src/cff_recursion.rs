use std::collections::{BTreeSet, HashSet};

use serde::{Deserialize, Serialize};

use crate::{LinearEnergyExpr, ParsedGraph};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
struct EdgeRef {
    edge_id: usize,
    edge_type: EdgeType,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
enum EdgeType {
    Virtual,
    External,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
struct CffVertex {
    nodes: BTreeSet<usize>,
    incoming: Vec<EdgeRef>,
    outgoing: Vec<EdgeRef>,
}

impl CffVertex {
    fn vertex_type(&self) -> VertexType {
        let is_sink = self
            .outgoing
            .iter()
            .all(|edge| edge.edge_type != EdgeType::Virtual);
        let is_source = self
            .incoming
            .iter()
            .all(|edge| edge.edge_type != EdgeType::Virtual);
        if is_sink {
            VertexType::Sink
        } else if is_source {
            VertexType::Source
        } else {
            VertexType::None
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VertexType {
    Source,
    Sink,
    None,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
struct CffGenerationGraph {
    vertices: Vec<CffVertex>,
}

impl CffGenerationGraph {
    fn new(mut vertices: Vec<CffVertex>) -> Self {
        vertices.sort_by_key(|vertex| vertex.nodes.iter().copied().collect::<Vec<_>>());
        Self { vertices }
    }

    fn vertex(&self, nodes: &BTreeSet<usize>) -> &CffVertex {
        self.vertices
            .iter()
            .find(|vertex| &vertex.nodes == nodes)
            .expect("vertex must exist")
    }

    fn virtual_adjacency(&self) -> Vec<(BTreeSet<usize>, BTreeSet<usize>)> {
        let mut adjacency = Vec::new();
        for vertex in &self.vertices {
            for edge in vertex
                .outgoing
                .iter()
                .filter(|edge| edge.edge_type == EdgeType::Virtual)
            {
                if let Some(other) = self
                    .vertices
                    .iter()
                    .find(|other| other.nodes != vertex.nodes && other.incoming.contains(edge))
                {
                    adjacency.push((vertex.nodes.clone(), other.nodes.clone()));
                }
            }
        }
        adjacency
    }

    fn has_directed_cycle(&self) -> bool {
        let adjacency_edges = self.virtual_adjacency();
        let mut visited = HashSet::new();
        let mut stack = HashSet::new();

        fn dfs(
            node: &BTreeSet<usize>,
            adjacency_edges: &[(BTreeSet<usize>, BTreeSet<usize>)],
            visited: &mut HashSet<BTreeSet<usize>>,
            stack: &mut HashSet<BTreeSet<usize>>,
        ) -> bool {
            if stack.contains(node) {
                return true;
            }
            if !visited.insert(node.clone()) {
                return false;
            }
            stack.insert(node.clone());
            for (_, next) in adjacency_edges.iter().filter(|(source, _)| source == node) {
                if dfs(next, adjacency_edges, visited, stack) {
                    return true;
                }
            }
            stack.remove(node);
            false
        }

        self.vertices
            .iter()
            .any(|vertex| dfs(&vertex.nodes, &adjacency_edges, &mut visited, &mut stack))
    }

    fn are_directed_adjacent(&self, left: &BTreeSet<usize>, right: &BTreeSet<usize>) -> bool {
        let left_vertex = self.vertex(left);
        let right_vertex = self.vertex(right);
        left_vertex
            .outgoing
            .iter()
            .any(|edge| edge.edge_type == EdgeType::Virtual && right_vertex.incoming.contains(edge))
    }

    fn are_adjacent(&self, left: &BTreeSet<usize>, right: &BTreeSet<usize>) -> bool {
        self.are_directed_adjacent(left, right) || self.are_directed_adjacent(right, left)
    }

    fn undirected_neighbours(&self, nodes: &BTreeSet<usize>) -> Vec<&CffVertex> {
        self.vertices
            .iter()
            .filter(|vertex| &vertex.nodes != nodes && self.are_adjacent(nodes, &vertex.nodes))
            .collect()
    }

    fn has_connected_complement(&self, removed_nodes: &BTreeSet<usize>) -> bool {
        let others = self
            .vertices
            .iter()
            .filter(|vertex| &vertex.nodes != removed_nodes)
            .map(|vertex| vertex.nodes.clone())
            .collect::<Vec<_>>();
        let Some(start) = others.first() else {
            return true;
        };

        let mut visited = HashSet::from([start.clone()]);
        let mut frontier = vec![start.clone()];
        while let Some(current) = frontier.pop() {
            for neighbour in self.undirected_neighbours(&current) {
                if &neighbour.nodes == removed_nodes {
                    continue;
                }
                if visited.insert(neighbour.nodes.clone()) {
                    frontier.push(neighbour.nodes.clone());
                }
            }
        }
        visited.len() == others.len()
    }

    fn source_sink_greedy(&self) -> Option<&CffVertex> {
        self.vertices.iter().find(|vertex| {
            matches!(vertex.vertex_type(), VertexType::Source | VertexType::Sink)
                && self.has_connected_complement(&vertex.nodes)
        })
    }

    fn contract_vertices(
        &self,
        left_nodes: &BTreeSet<usize>,
        right_nodes: &BTreeSet<usize>,
    ) -> Self {
        let left = self.vertex(left_nodes);
        let right = self.vertex(right_nodes);
        let right_outgoing = right.outgoing.iter().copied().collect::<HashSet<_>>();
        let right_incoming = right.incoming.iter().copied().collect::<HashSet<_>>();
        let left_outgoing = left.outgoing.iter().copied().collect::<HashSet<_>>();
        let left_incoming = left.incoming.iter().copied().collect::<HashSet<_>>();

        let mut incoming = left
            .incoming
            .iter()
            .filter(|edge| !right_outgoing.contains(edge))
            .chain(
                right
                    .incoming
                    .iter()
                    .filter(|edge| !left_outgoing.contains(edge)),
            )
            .copied()
            .collect::<Vec<_>>();
        let mut outgoing = left
            .outgoing
            .iter()
            .filter(|edge| !right_incoming.contains(edge))
            .chain(
                right
                    .outgoing
                    .iter()
                    .filter(|edge| !left_incoming.contains(edge)),
            )
            .copied()
            .collect::<Vec<_>>();
        incoming.sort();
        outgoing.sort();

        let mut nodes = left.nodes.clone();
        nodes.extend(right.nodes.iter().copied());
        let mut vertices = self
            .vertices
            .iter()
            .filter(|vertex| &vertex.nodes != left_nodes && &vertex.nodes != right_nodes)
            .cloned()
            .collect::<Vec<_>>();
        vertices.push(CffVertex {
            nodes,
            incoming,
            outgoing,
        });
        Self::new(vertices)
    }

    fn reverse_virtual_edge(&mut self, edge_id: usize) {
        let edge = EdgeRef {
            edge_id,
            edge_type: EdgeType::Virtual,
        };
        for vertex in &mut self.vertices {
            if let Some(position) = vertex.outgoing.iter().position(|item| *item == edge) {
                vertex.outgoing.remove(position);
                vertex.incoming.push(edge);
                vertex.incoming.sort();
            } else if let Some(position) = vertex.incoming.iter().position(|item| *item == edge) {
                vertex.incoming.remove(position);
                vertex.outgoing.push(edge);
                vertex.outgoing.sort();
            }
        }
    }
}

pub(crate) fn enumerate_cff_surface_chains(
    parsed: &ParsedGraph,
    edge_signs: &[i32],
) -> Vec<Vec<LinearEnergyExpr>> {
    let mut graph = build_base_graph_from_parsed(parsed);
    for (edge_id, sign) in edge_signs.iter().enumerate() {
        if *sign < 0 {
            graph.reverse_virtual_edge(edge_id);
        }
    }
    if graph.has_directed_cycle() {
        return Vec::new();
    }

    let mut branches = Vec::new();
    enumerate_cff_branches(&graph, parsed, &mut branches);
    branches
}

fn build_base_graph_from_parsed(parsed: &ParsedGraph) -> CffGenerationGraph {
    let n_vertices = parsed
        .node_name_to_internal
        .values()
        .copied()
        .max()
        .map(|value| value + 1)
        .unwrap_or(0);
    let mut vertices = (0..n_vertices)
        .map(|vertex_id| CffVertex {
            nodes: BTreeSet::from([vertex_id]),
            incoming: Vec::new(),
            outgoing: Vec::new(),
        })
        .collect::<Vec<_>>();

    for edge in &parsed.internal_edges {
        let edge_ref = EdgeRef {
            edge_id: edge.edge_id,
            edge_type: EdgeType::Virtual,
        };
        vertices[edge.tail].outgoing.push(edge_ref);
        vertices[edge.head].incoming.push(edge_ref);
    }
    for edge in &parsed.external_edges {
        let edge_ref = EdgeRef {
            edge_id: edge.edge_id,
            edge_type: EdgeType::External,
        };
        if let Some(source) = edge.source {
            vertices[source].outgoing.push(edge_ref);
        }
        if let Some(destination) = edge.destination {
            vertices[destination].incoming.push(edge_ref);
        }
    }
    for vertex in &mut vertices {
        vertex.incoming.sort();
        vertex.outgoing.sort();
    }
    CffGenerationGraph::new(vertices)
}

fn enumerate_cff_branches(
    graph: &CffGenerationGraph,
    parsed: &ParsedGraph,
    branch_acc: &mut Vec<Vec<LinearEnergyExpr>>,
) {
    if graph.vertices.len() < 2 {
        branch_acc.push(Vec::new());
        return;
    }
    let Some(vertex) = graph.source_sink_greedy() else {
        return;
    };
    let surface = cff_surface_for_vertex(parsed, vertex);
    if graph.vertices.len() == 2 {
        branch_acc.push(vec![surface]);
        return;
    }

    let mut emitted = false;
    for neighbour in graph.undirected_neighbours(&vertex.nodes) {
        let child = graph.contract_vertices(&vertex.nodes, &neighbour.nodes);
        if child.has_directed_cycle() {
            continue;
        }
        let mut sub = Vec::new();
        enumerate_cff_branches(&child, parsed, &mut sub);
        for mut chain in sub {
            let mut full_chain = vec![surface.clone()];
            full_chain.append(&mut chain);
            branch_acc.push(full_chain);
            emitted = true;
        }
    }
    if !emitted {
        branch_acc.push(vec![surface]);
    }
}

fn cff_surface_for_vertex(parsed: &ParsedGraph, vertex: &CffVertex) -> LinearEnergyExpr {
    let mut expr = LinearEnergyExpr::zero();
    for edge in vertex
        .incoming
        .iter()
        .chain(vertex.outgoing.iter())
        .filter(|edge| edge.edge_type == EdgeType::Virtual)
    {
        expr =
            expr + LinearEnergyExpr::ose(linnet::half_edge::involution::EdgeIndex(edge.edge_id), 1);
    }

    let mut external_shift = boundary_external_shift_from_internal_labels(parsed, &vertex.nodes);
    if vertex.vertex_type() == VertexType::Source {
        for coeff in external_shift.values_mut() {
            *coeff *= -1;
        }
    }
    for (external_id, coeff) in external_shift {
        expr = expr
            + LinearEnergyExpr::external(
                linnet::half_edge::involution::EdgeIndex(external_id),
                i64::from(coeff),
            );
    }
    expr.canonical()
}

fn boundary_external_shift_from_internal_labels(
    parsed: &ParsedGraph,
    node_set: &BTreeSet<usize>,
) -> std::collections::BTreeMap<usize, i32> {
    let mut acc = std::collections::BTreeMap::<usize, i32>::new();
    for edge in &parsed.internal_edges {
        let sign = if node_set.contains(&edge.tail) && !node_set.contains(&edge.head) {
            1
        } else if node_set.contains(&edge.head) && !node_set.contains(&edge.tail) {
            -1
        } else {
            continue;
        };
        for (external_id, coeff) in edge.signature.external_signature.iter().enumerate() {
            if *coeff != 0 {
                *acc.entry(external_id).or_default() += sign * *coeff;
            }
        }
    }
    acc.retain(|_, coeff| *coeff != 0);
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cff_recursion_finds_nontrivial_box_pow3_contraction_branch() {
        let parsed = crate::graph_io::test_graphs::box_pow3_graph();
        let has_nontrivial = (0..(1usize << parsed.internal_edges.len())).any(|bitmask| {
            let signs = (0..parsed.internal_edges.len())
                .map(|edge_index| {
                    if bitmask & (1usize << edge_index) == 0 {
                        1
                    } else {
                        -1
                    }
                })
                .collect::<Vec<_>>();
            enumerate_cff_surface_chains(&parsed, &signs)
                .iter()
                .any(|branch| branch.len() > 1)
        });

        assert!(has_nontrivial);
    }
}
