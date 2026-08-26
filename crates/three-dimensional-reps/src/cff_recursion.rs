use std::collections::{BTreeMap, BTreeSet, HashSet};

use serde::{Deserialize, Serialize};

use crate::{
    LinearEnergyExpr, MediumMode, ParsedGraph, ThermalDistributionFactor, ThermalNumerator,
    ThermalWeight,
};
use linnet::half_edge::involution::EdgeIndex;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
struct EdgeRef {
    edge_id: usize,
    edge_type: EdgeType,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
enum EdgeType {
    Virtual,
    External,
    InitialStateCut,
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

    fn strip_thermal_distribution_factors(
        &mut self,
        signs: &[i32],
    ) -> Vec<ThermalDistributionFactor> {
        let mut factors = Vec::new();
        loop {
            let self_edges = self
                .vertices
                .iter()
                .flat_map(|vertex| {
                    vertex
                        .incoming
                        .iter()
                        .filter(|edge| {
                            edge.edge_type == EdgeType::Virtual && vertex.outgoing.contains(edge)
                        })
                        .map(|edge| edge.edge_id)
                })
                .collect::<BTreeSet<_>>();
            let mut removed = self_edges.clone();
            factors.extend(
                self_edges
                    .into_iter()
                    .map(|edge_id| ThermalDistributionFactor {
                        edge_id: EdgeIndex(edge_id),
                        sign: signs[edge_id],
                        derivative_order: 0,
                    }),
            );
            self.remove_virtual_edges(&removed);
            if let Some(cycle) = self.vertices.iter().find_map(|vertex| {
                self.detachable_cycle(
                    &vertex.nodes,
                    &vertex.nodes,
                    &mut vec![vertex.nodes.clone()],
                    &mut Vec::new(),
                )
            }) {
                factors.push(ThermalDistributionFactor {
                    edge_id: EdgeIndex(*cycle.iter().min().expect("cycle has edges")),
                    sign: 1,
                    derivative_order: cycle.len() - 1,
                });
                removed.extend(cycle.iter().copied());
                self.remove_virtual_edges(&cycle.into_iter().collect());
            }
            if removed.is_empty() {
                break;
            }
        }
        factors
    }

    fn remove_virtual_edges(&mut self, edges: &BTreeSet<usize>) {
        for vertex in &mut self.vertices {
            vertex.incoming.retain(|edge| {
                edge.edge_type != EdgeType::Virtual || !edges.contains(&edge.edge_id)
            });
            vertex.outgoing.retain(|edge| {
                edge.edge_type != EdgeType::Virtual || !edges.contains(&edge.edge_id)
            });
        }
        self.vertices
            .retain(|vertex| !vertex.incoming.is_empty() || !vertex.outgoing.is_empty());
    }

    fn detachable_cycle(
        &self,
        start: &BTreeSet<usize>,
        current: &BTreeSet<usize>,
        visited: &mut Vec<BTreeSet<usize>>,
        path: &mut Vec<usize>,
    ) -> Option<Vec<usize>> {
        for edge in self
            .vertex(current)
            .outgoing
            .iter()
            .filter(|edge| edge.edge_type == EdgeType::Virtual)
        {
            let Some(next) = self
                .vertices
                .iter()
                .find(|vertex| vertex.incoming.contains(edge))
            else {
                continue;
            };
            if &next.nodes == start {
                if path.is_empty() {
                    continue;
                }
                let mut cycle = path.clone();
                cycle.push(edge.edge_id);
                let attachments = visited
                    .iter()
                    .filter(|nodes| {
                        let vertex = self.vertex(nodes);
                        vertex.incoming.iter().chain(&vertex.outgoing).any(|edge| {
                            edge.edge_type != EdgeType::Virtual || !cycle.contains(&edge.edge_id)
                        })
                    })
                    .count();
                if attachments <= 1 {
                    return Some(cycle);
                }
            } else if !visited.contains(&next.nodes) {
                visited.push(next.nodes.clone());
                path.push(edge.edge_id);
                if let Some(cycle) = self.detachable_cycle(start, &next.nodes, visited, path) {
                    return Some(cycle);
                }
                path.pop();
                visited.pop();
            }
        }
        None
    }

    fn reverse_virtual_edge(&mut self, edge_id: usize) {
        let edge = EdgeRef {
            edge_id,
            edge_type: EdgeType::Virtual,
        };
        for vertex in &mut self.vertices {
            if vertex.incoming.contains(&edge) && vertex.outgoing.contains(&edge) {
                continue;
            }
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

pub(crate) struct CffSurfaceChain {
    pub surfaces: Vec<LinearEnergyExpr>,
    pub thermal_weight: ThermalWeight,
    pub sign: i32,
}

pub(crate) fn enumerate_cff_surface_chains(
    parsed: &ParsedGraph,
    edge_signs: &[i32],
    medium_mode: MediumMode,
) -> Vec<CffSurfaceChain> {
    let mut graph = build_base_graph_from_parsed(parsed);
    for (edge_id, sign) in edge_signs.iter().enumerate() {
        if *sign < 0 {
            graph.reverse_virtual_edge(edge_id);
        }
    }
    if medium_mode == MediumMode::Vacuum && graph.has_directed_cycle() {
        return Vec::new();
    }
    let mut branches = Vec::new();
    enumerate_cff_branches(&graph, parsed, edge_signs, medium_mode, &mut branches);
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
        let edge_type = if parsed.is_initial_state_cut_edge(edge.edge_id) {
            EdgeType::InitialStateCut
        } else {
            EdgeType::Virtual
        };
        let edge_ref = EdgeRef {
            edge_id: edge.edge_id,
            edge_type,
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
    // Isolated zero-edge components contribute the multiplicative identity;
    // retaining them would manufacture a causal surface with value zero.
    vertices.retain(|vertex| !vertex.incoming.is_empty() || !vertex.outgoing.is_empty());
    CffGenerationGraph::new(vertices)
}

fn enumerate_cff_branches(
    graph: &CffGenerationGraph,
    parsed: &ParsedGraph,
    edge_signs: &[i32],
    medium_mode: MediumMode,
    branch_acc: &mut Vec<CffSurfaceChain>,
) {
    let mut graph = graph.clone();
    let thermal = medium_mode != MediumMode::Vacuum;
    let weight = ThermalWeight {
        medium_mode,
        distributions: if thermal {
            graph.strip_thermal_distribution_factors(edge_signs)
        } else {
            Vec::new()
        },
        numerators: Vec::new(),
    };
    if graph.vertices.len() < 2 {
        branch_acc.push(CffSurfaceChain {
            surfaces: Vec::new(),
            thermal_weight: weight,
            sign: 1,
        });
        return;
    }
    // First, try to find a source or sink with connected complement.
    // Thermal orientations can have no source or sink; then contract a vertex
    // of degree at least three with connected complement.
    let Some(vertex) = graph.source_sink_greedy().or_else(|| {
        thermal
            .then(|| {
                graph.vertices.iter().find(|vertex| {
                    vertex.incoming.len() + vertex.outgoing.len() >= 3
                        && graph.has_connected_complement(&vertex.nodes)
                })
            })
            .flatten()
    }) else {
        return;
    };
    let (surface, surface_sign) = if thermal {
        cff_surface_for_vertex_thermal(parsed, vertex)
    } else {
        (cff_surface_for_vertex(parsed, vertex), 1)
    };
    if !thermal && graph.vertices.len() == 2 {
        branch_acc.push(CffSurfaceChain {
            surfaces: vec![surface],
            thermal_weight: weight,
            sign: 1,
        });
        return;
    }

    let mut emitted = false;
    for neighbour in graph.undirected_neighbours(&vertex.nodes) {
        let child = graph.contract_vertices(&vertex.nodes, &neighbour.nodes);
        if !thermal && child.has_directed_cycle() {
            continue;
        }
        let mut branch_weight = weight.clone();
        let mut sign = 1;
        if thermal {
            // Edges joining the contracted vertices carry the thermal numerator
            // for that contraction, with outgoing minus incoming ordering.
            let outgoing = vertex
                .outgoing
                .iter()
                .filter(|edge| neighbour.incoming.contains(edge))
                .map(|edge| EdgeIndex(edge.edge_id))
                .collect();
            let incoming = vertex
                .incoming
                .iter()
                .filter(|edge| neighbour.outgoing.contains(edge))
                .map(|edge| EdgeIndex(edge.edge_id))
                .collect();
            let (numerator, numerator_sign) =
                ThermalNumerator::from_edge_lists_canonicalized(outgoing, incoming);
            sign = surface_sign * numerator_sign;
            if !numerator.is_trivial() {
                branch_weight.numerators.push(numerator);
            }
        }
        let mut sub = Vec::new();
        enumerate_cff_branches(&child, parsed, edge_signs, medium_mode, &mut sub);
        for mut chain in sub {
            chain.surfaces.insert(0, surface.clone());
            chain.thermal_weight = branch_weight.product(&chain.thermal_weight);
            chain.sign *= sign;
            branch_acc.push(chain);
            emitted = true;
        }
    }
    if !emitted && !thermal {
        branch_acc.push(CffSurfaceChain {
            surfaces: vec![surface],
            thermal_weight: weight,
            sign: 1,
        });
    }
}

fn cff_surface_for_vertex_thermal(
    parsed: &ParsedGraph,
    vertex: &CffVertex,
) -> (LinearEnergyExpr, i32) {
    let virtual_ids = |edges: &[EdgeRef]| {
        edges
            .iter()
            .filter(|edge| edge.edge_type == EdgeType::Virtual)
            .map(|edge| edge.edge_id)
            .collect::<Vec<_>>()
    };
    let incoming = virtual_ids(&vertex.incoming);
    let outgoing = virtual_ids(&vertex.outgoing);
    // Canonicalize an H-surface with the larger positive side, using the
    // smaller edge IDs as the tie break when the two sides have equal size.
    let flip = match incoming.len().cmp(&outgoing.len()) {
        std::cmp::Ordering::Greater => true,
        std::cmp::Ordering::Less => false,
        std::cmp::Ordering::Equal => outgoing > incoming,
    };
    let sign = if flip { -1 } else { 1 };
    let mut expr = LinearEnergyExpr::zero();
    for edge in outgoing {
        expr = expr + LinearEnergyExpr::ose(EdgeIndex(edge), i64::from(sign));
    }
    for edge in incoming {
        expr = expr + LinearEnergyExpr::ose(EdgeIndex(edge), -i64::from(sign));
    }
    let mut shift = boundary_external_shift_from_internal_labels(parsed, &vertex.nodes);
    add_initial_state_cut_external_shift(parsed, vertex, &mut shift);
    for (edge, coeff) in shift {
        expr = expr + LinearEnergyExpr::external(EdgeIndex(edge), -i64::from(sign * coeff));
    }
    (expr.canonical(), sign)
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
    add_initial_state_cut_external_shift(parsed, vertex, &mut external_shift);
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

fn add_initial_state_cut_external_shift(
    parsed: &ParsedGraph,
    vertex: &CffVertex,
    external_shift: &mut BTreeMap<usize, i32>,
) {
    for (edge, incidence_sign) in vertex
        .incoming
        .iter()
        .filter(|edge| edge.edge_type == EdgeType::InitialStateCut)
        .map(|edge| (edge, -1))
        .chain(
            vertex
                .outgoing
                .iter()
                .filter(|edge| edge.edge_type == EdgeType::InitialStateCut)
                .map(|edge| (edge, 1)),
        )
    {
        let Some(cut_edge) = parsed.initial_state_cut_edge(edge.edge_id) else {
            continue;
        };
        *external_shift.entry(cut_edge.external_id).or_default() +=
            -incidence_sign * cut_edge.external_sign;
    }
    external_shift.retain(|_, coeff| *coeff != 0);
}

fn boundary_external_shift_from_internal_labels(
    parsed: &ParsedGraph,
    node_set: &BTreeSet<usize>,
) -> std::collections::BTreeMap<usize, i32> {
    let mut acc = std::collections::BTreeMap::<usize, i32>::new();
    let initial_state_external_ids = parsed
        .initial_state_cut_edges
        .iter()
        .map(|edge| edge.external_id)
        .collect::<BTreeSet<_>>();
    for edge in &parsed.internal_edges {
        if parsed.is_initial_state_cut_edge(edge.edge_id) {
            continue;
        }
        let sign = if node_set.contains(&edge.tail) && !node_set.contains(&edge.head) {
            1
        } else if node_set.contains(&edge.head) && !node_set.contains(&edge.tail) {
            -1
        } else {
            continue;
        };
        for (external_id, coeff) in edge.signature.external_signature.iter().enumerate() {
            if *coeff != 0 && !initial_state_external_ids.contains(&external_id) {
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
            enumerate_cff_surface_chains(&parsed, &signs, MediumMode::Vacuum)
                .iter()
                .any(|branch| branch.surfaces.len() > 1)
        });

        assert!(has_nontrivial);
    }
}

#[cfg(test)]
mod thermal_tests {
    use super::*;

    #[test]
    fn thermal_detachable_cycles_keep_distribution_derivatives() {
        let mut parsed = crate::graph_io::test_graphs::box_graph();
        parsed.external_edges.clear();
        parsed.external_names.clear();
        for edge in &mut parsed.internal_edges {
            edge.signature.external_signature.clear();
        }
        for (count, derivative_order) in [(2, 1), (3, 2)] {
            let mut cycle = parsed.clone();
            cycle.internal_edges.truncate(count);
            cycle.internal_edges[count - 1].head = 0;
            let chains = enumerate_cff_surface_chains(
                &cycle,
                &vec![1; count],
                MediumMode::ThermodynamicEquilibrium,
            );
            assert_eq!(chains.len(), 1);
            assert!(chains[0].surfaces.is_empty());
            assert_eq!(
                chains[0].thermal_weight.distributions,
                vec![ThermalDistributionFactor {
                    edge_id: EdgeIndex(0),
                    sign: 1,
                    derivative_order,
                }]
            );
        }
    }

    #[test]
    fn thermal_cycles_with_multiple_attachments_are_not_stripped() {
        let parsed = crate::graph_io::test_graphs::box_graph();
        let mut graph = build_base_graph_from_parsed(&parsed);
        assert!(graph.strip_thermal_distribution_factors(&[1; 4]).is_empty());
        let chains =
            enumerate_cff_surface_chains(&parsed, &[1; 4], MediumMode::ThermodynamicEquilibrium);
        assert!(!chains.is_empty());
        assert!(chains.iter().all(|chain| !chain.surfaces.is_empty()));
    }
}
