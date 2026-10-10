use std::collections::{BTreeMap, BTreeSet, HashSet};

use serde::{Deserialize, Serialize};
use symbolica::domains::integer::Integer;

use crate::{
    GenerationError, LinearEnergyExpr, MediumMode, ParsedGraph, ThermalDistributionFactor,
    ThermalNumerator, ThermalWeight, generation::Result, utils::Rational,
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
                && vertex
                    .incoming
                    .iter()
                    .chain(&vertex.outgoing)
                    .any(|edge| edge.edge_type == EdgeType::Virtual)
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

    fn strip_tadpoles_and_cyclic_chains(
        &mut self,
        parsed: &ParsedGraph,
        medium_mode: MediumMode,
        edge_signs: &[i32],
        preferred_tadpole_poles: &BTreeMap<usize, i32>,
    ) -> (Vec<ThermalDistributionFactor>, Rational) {
        let mut factors = Vec::new();
        let mut prefactor = Rational::one();
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
            for edge_id in self_edges {
                // Free tadpoles average the two closures. A contact keeps the
                // closure inherited from its parent, independently of the medium.
                // The physical pole sign still controls its numerator and mu.
                let contour_sign = preferred_tadpole_poles
                    .get(&edge_id)
                    .map_or(0, |pole| edge_signs[edge_id] * pole);
                if medium_mode == MediumMode::Vacuum {
                    prefactor *= Rational::new(1 + i64::from(contour_sign), 2);
                } else {
                    factors.push(ThermalDistributionFactor {
                        edge_id: EdgeIndex(edge_id),
                        sign: contour_sign,
                        derivative_order: 0,
                    });
                }
            }
            self.remove_virtual_edges(&removed);
            if medium_mode != MediumMode::Vacuum
                && self.has_directed_cycle()
                && let Some(cycle) = self.vertices.iter().find_map(|vertex| {
                    self.detachable_cycle(
                        &|cycle| {
                            let first = &parsed.internal_edges[cycle[0]];
                            let (signature, _) = first.signature.canonical_up_to_sign();
                            cycle.iter().all(|&edge_id| {
                                let edge = &parsed.internal_edges[edge_id];
                                edge.mass_key == first.mass_key
                                    && edge.signature.canonical_up_to_sign().0 == signature
                            })
                        },
                        &vertex.nodes,
                        &vertex.nodes,
                        &mut vec![vertex.nodes.clone()],
                        &mut Vec::new(),
                    )
                })
            {
                // An m-edge cyclic chain with a common pole contributes
                // (-1)^(m-1)/(m-1)! times the ordinary (m-1)th energy derivative.
                // Unequal poles stay in the recursion as divided differences;
                // self-loops retain their contour prescription above.
                let derivative_order = cycle.len() - 1;
                let sign = if derivative_order.is_multiple_of(2) {
                    1
                } else {
                    -1
                };
                prefactor *= Rational::from(sign)
                    / Rational::from(&Integer::factorial(derivative_order as u32));
                factors.push(ThermalDistributionFactor {
                    edge_id: EdgeIndex(*cycle.iter().min().expect("cycle has edges")),
                    sign: 1,
                    derivative_order,
                });
                removed.extend(cycle.iter().copied());
                self.remove_virtual_edges(&cycle.into_iter().collect());
            }
            if removed.is_empty() {
                break;
            }
        }
        (factors, prefactor)
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
        accept: &impl Fn(&[usize]) -> bool,
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
                if attachments <= 1 && accept(&cycle) {
                    return Some(cycle);
                }
            } else if !visited.contains(&next.nodes) {
                visited.push(next.nodes.clone());
                path.push(edge.edge_id);
                if let Some(cycle) =
                    self.detachable_cycle(accept, start, &next.nodes, visited, path)
                {
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
    pub prefactor: Rational,
}

pub(crate) fn enumerate_cff_surface_chains(
    parsed: &ParsedGraph,
    edge_signs: &[i32],
    medium_mode: MediumMode,
    preferred_tadpole_poles: &BTreeMap<usize, i32>,
) -> Result<Vec<CffSurfaceChain>> {
    let mut graph = build_base_graph_from_parsed(parsed);
    for (edge_id, sign) in edge_signs.iter().enumerate() {
        if *sign < 0 {
            graph.reverse_virtual_edge(edge_id);
        }
    }
    if medium_mode == MediumMode::Vacuum && graph.has_directed_cycle() {
        return Ok(Vec::new());
    }
    let mut branches = Vec::new();
    enumerate_cff_branches(
        graph,
        parsed,
        edge_signs,
        medium_mode,
        preferred_tadpole_poles,
        &mut branches,
    )?;
    Ok(branches)
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
        // A cut with both endpoints at the same vertex carries no boundary
        // momentum and must not act as an external attachment to a cycle.
        if edge_type == EdgeType::InitialStateCut && edge.tail == edge.head {
            continue;
        }
        let edge_ref = EdgeRef {
            edge_id: edge.edge_id,
            edge_type,
        };
        vertices[edge.tail].outgoing.push(edge_ref);
        vertices[edge.head].incoming.push(edge_ref);
    }
    for edge in &parsed.external_edges {
        if edge.source.is_some() && edge.source == edge.destination {
            continue;
        }
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
    mut graph: CffGenerationGraph,
    parsed: &ParsedGraph,
    edge_signs: &[i32],
    medium_mode: MediumMode,
    preferred_tadpole_poles: &BTreeMap<usize, i32>,
    branch_acc: &mut Vec<CffSurfaceChain>,
) -> Result<()> {
    let thermal = medium_mode != MediumMode::Vacuum;
    let (distributions, reduction_prefactor) = graph.strip_tadpoles_and_cyclic_chains(
        parsed,
        medium_mode,
        edge_signs,
        preferred_tadpole_poles,
    );
    if reduction_prefactor.is_zero() {
        return Ok(());
    }
    let weight = ThermalWeight {
        medium_mode,
        distributions,
        numerators: Vec::new(),
    };
    if graph.vertices.len() < 2 {
        branch_acc.push(CffSurfaceChain {
            surfaces: Vec::new(),
            thermal_weight: weight,
            prefactor: reduction_prefactor,
        });
        return Ok(());
    }
    // First, try to find a source or sink with connected complement.
    // For a mixed vertex, connected complement and virtual degree >= 3
    // provide independent loop flow through the boundary. A bivalent mixed
    // vertex instead needs external attachments on both sides: no attachments,
    // or all attachments on one side, allow an identically zero E1-E2 surface.
    let vertex = graph
        .source_sink_greedy()
        .or_else(|| {
            let external_vertices = graph
                .vertices
                .iter()
                .filter(|vertex| {
                    vertex
                        .incoming
                        .iter()
                        .chain(&vertex.outgoing)
                        .any(|edge| edge.edge_type != EdgeType::Virtual)
                })
                .map(|vertex| &vertex.nodes)
                .collect::<BTreeSet<_>>();
            graph.vertices.iter().find(|vertex| {
                let incoming = vertex
                    .incoming
                    .iter()
                    .filter(|edge| edge.edge_type == EdgeType::Virtual)
                    .count();
                let outgoing = vertex
                    .outgoing
                    .iter()
                    .filter(|edge| edge.edge_type == EdgeType::Virtual)
                    .count();
                (incoming + outgoing >= 3
                    || (incoming == 1
                        && outgoing == 1
                        && external_vertices.len() >= 2
                        && external_vertices.contains(&vertex.nodes)))
                    && graph.has_connected_complement(&vertex.nodes)
            })
        })
        .or_else(|| {
            // A detachable cyclic chain with unequal poles is a divided
            // difference. Contract a mass transition away from its attachment;
            // common-pole child cycles then use the derivative reduction above.
            // This pole-dependent operation is confined to certified chains.
            graph.vertices.iter().find_map(|start| {
                let cycle = graph.detachable_cycle(
                    &|cycle| {
                        let (signature, _) = parsed.internal_edges[cycle[0]]
                            .signature
                            .canonical_up_to_sign();
                        cycle.iter().all(|&edge_id| {
                            parsed.internal_edges[edge_id]
                                .signature
                                .canonical_up_to_sign()
                                .0
                                == signature
                        })
                    },
                    &start.nodes,
                    &start.nodes,
                    &mut vec![start.nodes.clone()],
                    &mut Vec::new(),
                )?;
                graph.vertices.iter().find(|vertex| {
                    match (vertex.incoming.as_slice(), vertex.outgoing.as_slice()) {
                        ([incoming], [outgoing])
                            if incoming.edge_type == EdgeType::Virtual
                                && outgoing.edge_type == EdgeType::Virtual
                                && cycle.contains(&incoming.edge_id)
                                && cycle.contains(&outgoing.edge_id) =>
                        {
                            // Distinct symbolic mass keys must represent distinct
                            // masses when evaluating this divided difference.
                            parsed.internal_edges[incoming.edge_id].mass_key
                                != parsed.internal_edges[outgoing.edge_id].mass_key
                        }
                        _ => false,
                    }
                })
            })
        })
        .ok_or_else(|| GenerationError::NoAdmissibleCffVertex {
            medium_mode,
            edge_signs: edge_signs.to_vec(),
            vertices: graph
                .vertices
                .iter()
                .map(|vertex| vertex.nodes.iter().copied().collect())
                .collect(),
        })?;
    let (surface, surface_sign) = cff_surface_for_vertex(parsed, vertex);

    for neighbour in graph.undirected_neighbours(&vertex.nodes) {
        let child = graph.contract_vertices(&vertex.nodes, &neighbour.nodes);
        if !thermal && child.has_directed_cycle() {
            continue;
        }
        let mut branch_weight = weight.clone();
        let mut prefactor = reduction_prefactor.clone();
        // Edges joining the contracted vertices carry the thermal numerator
        // for that contraction, with outgoing minus incoming ordering.
        let outgoing = vertex
            .outgoing
            .iter()
            .filter(|edge| edge.edge_type == EdgeType::Virtual && neighbour.incoming.contains(edge))
            .map(|edge| EdgeIndex(edge.edge_id))
            .collect();
        let incoming = vertex
            .incoming
            .iter()
            .filter(|edge| edge.edge_type == EdgeType::Virtual && neighbour.outgoing.contains(edge))
            .map(|edge| EdgeIndex(edge.edge_id))
            .collect();
        let (numerator, numerator_sign) =
            ThermalNumerator::from_edge_lists_canonicalized(outgoing, incoming);
        prefactor *= Rational::from(surface_sign * numerator_sign);
        if thermal && !numerator.is_trivial() {
            branch_weight.numerators.push(numerator);
        }
        let mut sub = Vec::new();
        enumerate_cff_branches(
            child,
            parsed,
            edge_signs,
            medium_mode,
            preferred_tadpole_poles,
            &mut sub,
        )?;
        for mut chain in sub {
            chain.surfaces.insert(0, surface.clone());
            chain.thermal_weight = branch_weight.product(&chain.thermal_weight);
            chain.prefactor *= &prefactor;
            branch_acc.push(chain);
        }
    }
    Ok(())
}

fn cff_surface_for_vertex(parsed: &ParsedGraph, vertex: &CffVertex) -> (LinearEnergyExpr, i32) {
    let virtual_ids = |edges: &[EdgeRef]| {
        edges
            .iter()
            .filter(|edge| edge.edge_type == EdgeType::Virtual)
            .map(|edge| EdgeIndex(edge.edge_id))
            .collect::<Vec<_>>()
    };
    let (sides, sign) = ThermalNumerator::from_edge_lists_canonicalized(
        virtual_ids(&vertex.outgoing),
        virtual_ids(&vertex.incoming),
    );
    let mut expr = LinearEnergyExpr::zero();
    for edge in sides.positive_energies {
        expr = expr + LinearEnergyExpr::ose(edge, 1);
    }
    for edge in sides.negative_energies {
        expr = expr + LinearEnergyExpr::ose(edge, -1);
    }
    let shift = boundary_external_shift_from_internal_labels(parsed, &vertex.nodes);
    for (edge, coeff) in shift {
        expr = expr + LinearEnergyExpr::external(EdgeIndex(edge), -i64::from(sign * coeff));
    }
    (expr.canonical(), sign)
}

fn boundary_external_shift_from_internal_labels(
    parsed: &ParsedGraph,
    node_set: &BTreeSet<usize>,
) -> std::collections::BTreeMap<usize, i32> {
    let mut acc = std::collections::BTreeMap::<usize, i32>::new();
    // Virtual routing already includes the flow from cut aliases and ordinary
    // external legs. Keep their combined boundary shift even when the same
    // external coordinate also labels a cut elsewhere in the connected graph.
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
    use symbolica::{
        atom::{Atom, AtomCore},
        parse,
    };

    #[test]
    fn tadpole_pruning_averages_both_poles_in_every_medium() {
        let mut parsed = crate::graph_io::test_graphs::box_graph();
        parsed.external_edges.clear();
        parsed.external_names.clear();
        parsed.loop_names = vec!["k".to_string(), "l".to_string()];
        parsed.internal_edges.truncate(2);
        for edge in &mut parsed.internal_edges {
            edge.tail = 0;
            edge.head = 0;
            edge.signature.loop_signature = (0..2)
                .map(|index| i32::from(index == edge.edge_id))
                .collect();
            edge.signature.external_signature.clear();
        }
        for medium in [
            MediumMode::Vacuum,
            MediumMode::ThermodynamicEquilibrium,
            MediumMode::ZeroTemperatureEquilibrium,
        ] {
            for signs in [[1, 1], [1, -1], [-1, 1], [-1, -1]] {
                let chains =
                    enumerate_cff_surface_chains(&parsed, &signs, medium, &BTreeMap::new())
                        .unwrap();
                assert_eq!(chains.len(), 1);
                let chain = &chains[0];
                assert!(chain.surfaces.is_empty());
                assert!(chain.thermal_weight.numerators.is_empty());
                if medium == MediumMode::Vacuum {
                    assert_eq!(chain.prefactor, Rational::new(1, 4));
                    assert!(chain.thermal_weight.distributions.is_empty());
                } else {
                    assert_eq!(chain.prefactor, Rational::one());
                    assert_eq!(
                        chain.thermal_weight.distributions,
                        (0..2)
                            .map(|edge_id| ThermalDistributionFactor {
                                edge_id: EdgeIndex(edge_id),
                                sign: 0,
                                derivative_order: 0,
                            })
                            .collect::<Vec<_>>()
                    );
                }
            }
        }
    }

    #[test]
    fn attached_tadpole_weights_survive_vertex_contractions() {
        for count in [2, 3] {
            let mut cycle = crate::graph_io::test_graphs::box_graph();
            cycle.external_edges.clear();
            cycle.external_names.clear();
            cycle.internal_edges.truncate(count);
            cycle.internal_edges[count - 1].head = 0;
            for edge in &mut cycle.internal_edges {
                edge.signature.external_signature.clear();
            }
            let mut attached = cycle.clone();
            attached.loop_names.push("tadpole".to_string());
            for edge in &mut attached.internal_edges {
                edge.signature.loop_signature.push(0);
            }
            let mut tadpole = attached.internal_edges[0].clone();
            tadpole.edge_id = count;
            tadpole.head = tadpole.tail;
            tadpole.signature.loop_signature = vec![0, 1];
            attached.internal_edges.push(tadpole);
            let mut signs = vec![1; count];
            signs[count - 1] = -1;

            for medium in [
                MediumMode::Vacuum,
                MediumMode::ThermodynamicEquilibrium,
                MediumMode::ZeroTemperatureEquilibrium,
            ] {
                let base =
                    enumerate_cff_surface_chains(&cycle, &signs, medium, &BTreeMap::new()).unwrap();
                assert!(!base.is_empty());
                for sign in [-1, 1] {
                    let mut attached_signs = signs.clone();
                    attached_signs.push(sign);
                    let chains = enumerate_cff_surface_chains(
                        &attached,
                        &attached_signs,
                        medium,
                        &BTreeMap::new(),
                    )
                    .unwrap();
                    assert_eq!(chains.len(), base.len());
                    for (chain, base) in chains.iter().zip(&base) {
                        assert_eq!(chain.surfaces, base.surfaces);
                        let mut expected_weight = base.thermal_weight.clone();
                        let expected_prefactor = if medium == MediumMode::Vacuum {
                            &base.prefactor * &Rational::new(1, 2)
                        } else {
                            expected_weight
                                .distributions
                                .push(ThermalDistributionFactor {
                                    edge_id: EdgeIndex(count),
                                    sign: 0,
                                    derivative_order: 0,
                                });
                            expected_weight.canonicalize();
                            base.prefactor.clone()
                        };
                        assert_eq!(chain.prefactor, expected_prefactor);
                        assert_eq!(chain.thermal_weight, expected_weight);
                    }
                }
            }
        }
    }

    #[test]
    fn acyclic_contractions_match_the_thermal_vacuum_projection() {
        for tadpole_vertex in [None, Some(0), Some(2)] {
            let mut parsed = crate::graph_io::test_graphs::box_graph();
            if let Some(vertex) = tadpole_vertex {
                parsed.loop_names.push("tadpole".to_string());
                for edge in &mut parsed.internal_edges {
                    edge.signature.loop_signature.push(0);
                }
                let mut tadpole = parsed.internal_edges[0].clone();
                tadpole.edge_id = parsed.internal_edges.len();
                tadpole.tail = vertex;
                tadpole.head = vertex;
                tadpole.signature.loop_signature = vec![0, 1];
                parsed.internal_edges.push(tadpole);
            }
            for orientation in 0..1 << parsed.internal_edges.len() {
                let signs = (0..parsed.internal_edges.len())
                    .map(|edge| {
                        if orientation & (1 << edge) == 0 {
                            1
                        } else {
                            -1
                        }
                    })
                    .collect::<Vec<_>>();
                let mut graph = build_base_graph_from_parsed(&parsed);
                for (edge, sign) in signs.iter().enumerate() {
                    if *sign < 0 {
                        graph.reverse_virtual_edge(edge);
                    }
                }
                if graph.has_directed_cycle() {
                    continue;
                }
                let vacuum = enumerate_cff_surface_chains(
                    &parsed,
                    &signs,
                    MediumMode::Vacuum,
                    &BTreeMap::new(),
                )
                .unwrap()
                .into_iter()
                .map(|chain| (chain.surfaces, Atom::num(chain.prefactor)))
                .collect::<Vec<_>>();
                assert!(!vacuum.is_empty());
                assert!(vacuum.iter().all(|(surfaces, _)| surfaces.len() == 3));
                for medium in [
                    MediumMode::ThermodynamicEquilibrium,
                    MediumMode::ZeroTemperatureEquilibrium,
                ] {
                    let finite = medium.is_finite_temperature();
                    let projected =
                        enumerate_cff_surface_chains(&parsed, &signs, medium, &BTreeMap::new())
                            .unwrap()
                            .into_iter()
                            .filter_map(|chain| {
                                let mut coefficient = Atom::num(chain.prefactor);
                                for numerator in &chain.thermal_weight.numerators {
                                    let mut body = numerator.to_atom(finite);
                                    for edge in numerator
                                        .positive_energies
                                        .iter()
                                        .chain(&numerator.negative_energies)
                                    {
                                        for sign in [-1, 1] {
                                            body = body
                                                .replace(
                                                    ThermalDistributionFactor {
                                                        edge_id: *edge,
                                                        sign,
                                                        derivative_order: 0,
                                                    }
                                                    .to_atom(finite),
                                                )
                                                .with(ThermalDistributionFactor::vacuum_atom(
                                                    0,
                                                    Atom::num(sign),
                                                ));
                                        }
                                    }
                                    coefficient *= body;
                                }
                                for factor in chain.thermal_weight.distributions {
                                    coefficient *= ThermalDistributionFactor::vacuum_atom(
                                        factor.derivative_order,
                                        Atom::num(factor.sign),
                                    );
                                }
                                (!coefficient.is_zero()).then_some((chain.surfaces, coefficient))
                            })
                            .collect::<Vec<_>>();
                    // Compare every ordered surface chain and its exact weight,
                    // including the terminal contraction and tadpole half-weight.
                    assert_eq!(
                        projected, vacuum,
                        "{tadpole_vertex:?}, {signs:?}, {medium:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn thermal_surface_canonicalization_preserves_external_shift_sign() {
        let parsed = crate::graph_io::test_graphs::box_graph();
        for (outgoing, incoming, positive, negative, expected_sign) in [
            (vec![2, 0], vec![3, 1], vec![0, 2], vec![1, 3], 1),
            (vec![3, 1], vec![2, 0], vec![0, 2], vec![1, 3], -1),
            (vec![3], vec![2, 0, 1], vec![0, 1, 2], vec![3], -1),
        ] {
            let vertex = CffVertex {
                nodes: BTreeSet::from([0, 2]),
                outgoing: outgoing
                    .into_iter()
                    .map(|edge_id| EdgeRef {
                        edge_id,
                        edge_type: EdgeType::Virtual,
                    })
                    .collect(),
                incoming: incoming
                    .into_iter()
                    .map(|edge_id| EdgeRef {
                        edge_id,
                        edge_type: EdgeType::Virtual,
                    })
                    .collect(),
            };
            let (surface, sign) = cff_surface_for_vertex(&parsed, &vertex);
            assert_eq!(sign, expected_sign);
            let expected = LinearEnergyExpr {
                internal_terms: positive
                    .into_iter()
                    .map(|edge| (EdgeIndex(edge), Rational::from(1)))
                    .chain(
                        negative
                            .into_iter()
                            .map(|edge| (EdgeIndex(edge), Rational::from(-1))),
                    )
                    .collect(),
                external_terms: vec![
                    (EdgeIndex(0), Rational::from(expected_sign)),
                    (EdgeIndex(2), Rational::from(expected_sign)),
                ],
                ..LinearEnergyExpr::zero()
            };
            assert_eq!(surface, expected.canonical());
        }
    }

    #[test]
    fn thermal_detachable_cycles_respect_attachment_vertices() {
        for (name, edges, incoming, expected_cycles) in [
            (
                "one attachment among multiple candidate cycles",
                vec![
                    (0, 1),
                    (1, 2),
                    (2, 0), // Vertex 2 prevents detaching this cycle.
                    (0, 3),
                    (3, 0), // Only vertex 0 attaches this cycle to the remainder.
                    (0, 4),
                    (5, 6),
                    (6, 5), // Both vertices have external attachments.
                ],
                vec![
                    (0, EdgeType::External),
                    (0, EdgeType::InitialStateCut),
                    (2, EdgeType::InitialStateCut),
                    (5, EdgeType::External),
                    (6, EdgeType::External),
                ],
                vec![(0, Some(vec![8, 9])), (5, None)],
            ),
            (
                "two attachment vertices",
                vec![(0, 1), (1, 2), (2, 0), (0, 3)],
                vec![(2, EdgeType::InitialStateCut)],
                vec![(0, None)],
            ),
            (
                "external attachments on one vertex",
                vec![(0, 1), (1, 2), (2, 0)],
                vec![(0, EdgeType::External), (0, EdgeType::InitialStateCut)],
                vec![(0, Some(vec![2, 3, 4]))],
            ),
            (
                "detachable cycle inside a larger directed component",
                vec![(0, 1), (1, 2), (2, 0), (0, 3), (3, 0), (0, 4)],
                vec![
                    (0, EdgeType::External),
                    (0, EdgeType::InitialStateCut),
                    (2, EdgeType::InitialStateCut),
                ],
                vec![(0, Some(vec![6, 7]))],
            ),
            (
                "branching internal edges",
                vec![(0, 1), (1, 2), (2, 0), (0, 2), (2, 1)],
                vec![],
                vec![(0, None)],
            ),
        ] {
            let vertex_count = edges.iter().flat_map(|&(a, b)| [a, b]).max().unwrap() + 1;
            let mut vertices = (0..vertex_count)
                .map(|node| CffVertex {
                    nodes: BTreeSet::from([node]),
                    incoming: Vec::new(),
                    outgoing: Vec::new(),
                })
                .collect::<Vec<_>>();
            // Keep nonvirtual edges first, preserving the original fixture edge IDs.
            for (edge_id, &(node, edge_type)) in incoming.iter().enumerate() {
                vertices[node].incoming.push(EdgeRef { edge_id, edge_type });
            }
            for (index, &(tail, head)) in edges.iter().enumerate() {
                let edge = EdgeRef {
                    edge_id: incoming.len() + index,
                    edge_type: EdgeType::Virtual,
                };
                vertices[tail].outgoing.push(edge);
                vertices[head].incoming.push(edge);
            }
            let graph = CffGenerationGraph::new(vertices);
            for (start, expected) in expected_cycles {
                let start = BTreeSet::from([start]);
                assert_eq!(
                    graph.detachable_cycle(
                        &|_| true,
                        &start,
                        &start,
                        &mut vec![start.clone()],
                        &mut Vec::new(),
                    ),
                    expected,
                    "{name}",
                );
            }
        }
    }

    #[test]
    fn thermal_cycle_stripping_preserves_remaining_edges_and_vertices() {
        for (name, edges, expected_factors, remaining_edges, remaining_nodes) in [
            (
                "self edge with attached tail",
                vec![(0, 0), (0, 1)],
                vec![(0, 0)],
                vec![1],
                vec![0, 1],
            ),
            (
                "two-edge cycle with attached tail",
                vec![(0, 1), (1, 0), (0, 2)],
                vec![(0, 1)],
                vec![2],
                vec![0, 2],
            ),
            (
                "repeated stripping reaches a fixed point",
                vec![(0, 1), (1, 0), (2, 3), (3, 2), (0, 4)],
                vec![(0, 1), (2, 1)],
                vec![4],
                vec![0, 4],
            ),
            (
                "two attachments prevent stripping",
                vec![(0, 1), (1, 0), (0, 2), (1, 3)],
                vec![],
                vec![0, 1, 2, 3],
                vec![0, 1, 2, 3],
            ),
            (
                "only isolated cycle vertices are removed",
                vec![(0, 1), (1, 2), (2, 0), (0, 3)],
                vec![(0, 2)],
                vec![3],
                vec![0, 3],
            ),
        ] {
            let vertex_count = edges.iter().flat_map(|&(a, b)| [a, b]).max().unwrap() + 1;
            let mut vertices = (0..vertex_count)
                .map(|node| CffVertex {
                    nodes: BTreeSet::from([node]),
                    incoming: Vec::new(),
                    outgoing: Vec::new(),
                })
                .collect::<Vec<_>>();
            for (edge_id, &(tail, head)) in edges.iter().enumerate() {
                let edge = EdgeRef {
                    edge_id,
                    edge_type: EdgeType::Virtual,
                };
                vertices[tail].outgoing.push(edge);
                vertices[head].incoming.push(edge);
            }
            let mut graph = CffGenerationGraph::new(vertices);
            let mut parsed = crate::graph_io::test_graphs::box_graph();
            parsed.internal_edges = edges
                .iter()
                .enumerate()
                .map(|(edge_id, &(tail, head))| {
                    let mut edge = parsed.internal_edges[0].clone();
                    edge.edge_id = edge_id;
                    edge.tail = tail;
                    edge.head = head;
                    edge
                })
                .collect();
            assert_eq!(
                graph
                    .strip_tadpoles_and_cyclic_chains(
                        &parsed,
                        MediumMode::ThermodynamicEquilibrium,
                        &[],
                        &BTreeMap::new()
                    )
                    .0,
                expected_factors
                    .into_iter()
                    .map(|(edge, derivative_order)| ThermalDistributionFactor {
                        edge_id: EdgeIndex(edge),
                        sign: i32::from(derivative_order != 0),
                        derivative_order,
                    })
                    .collect::<Vec<_>>(),
                "{name}",
            );
            assert_eq!(
                graph
                    .vertices
                    .iter()
                    .map(|vertex| vertex.nodes.clone())
                    .collect::<Vec<_>>(),
                remaining_nodes
                    .iter()
                    .map(|&node| BTreeSet::from([node]))
                    .collect::<Vec<_>>(),
                "{name}",
            );
            for vertex in &graph.vertices {
                let node = *vertex.nodes.first().unwrap();
                for (actual, incoming) in [(&vertex.incoming, true), (&vertex.outgoing, false)] {
                    let expected = remaining_edges
                        .iter()
                        .filter(|&&edge| {
                            if incoming {
                                edges[edge].1 == node
                            } else {
                                edges[edge].0 == node
                            }
                        })
                        .map(|&edge_id| EdgeRef {
                            edge_id,
                            edge_type: EdgeType::Virtual,
                        })
                        .collect::<Vec<_>>();
                    assert_eq!(*actual, expected, "{name}");
                }
            }
            let fixed_point = graph.clone();
            assert!(
                graph
                    .strip_tadpoles_and_cyclic_chains(
                        &parsed,
                        MediumMode::ThermodynamicEquilibrium,
                        &[],
                        &BTreeMap::new()
                    )
                    .0
                    .is_empty(),
                "{name}"
            );
            assert_eq!(graph, fixed_point, "{name}");
        }
    }

    #[test]
    fn thermal_reversed_self_edge_keeps_adjacency_after_contraction() {
        let edge = |edge_id| EdgeRef {
            edge_id,
            edge_type: EdgeType::Virtual,
        };
        let graph = CffGenerationGraph::new(vec![
            CffVertex {
                nodes: BTreeSet::from([0]),
                incoming: vec![edge(0)],
                outgoing: vec![edge(0), edge(1)],
            },
            CffVertex {
                nodes: BTreeSet::from([1]),
                incoming: vec![edge(1)],
                outgoing: vec![edge(2)],
            },
            CffVertex {
                nodes: BTreeSet::from([2]),
                incoming: vec![edge(2)],
                outgoing: vec![],
            },
        ]);
        // Vertex contraction removes every joining edge; the self edge and tail survive.
        let mut graph = graph.contract_vertices(&BTreeSet::from([0]), &BTreeSet::from([1]));
        let contracted = graph.clone();
        graph.reverse_virtual_edge(0);
        assert_eq!(graph, contracted);
        let vertex = graph.vertex(&BTreeSet::from([0, 1]));
        assert_eq!(vertex.incoming, vec![edge(0)]);
        assert_eq!(vertex.outgoing, vec![edge(0), edge(2)]);
        let parsed = crate::graph_io::test_graphs::box_graph();
        assert_eq!(
            graph
                .strip_tadpoles_and_cyclic_chains(
                    &parsed,
                    MediumMode::ThermodynamicEquilibrium,
                    &[],
                    &BTreeMap::new()
                )
                .0,
            vec![ThermalDistributionFactor {
                edge_id: EdgeIndex(0),
                sign: 0,
                derivative_order: 0,
            }],
        );
        assert!(graph.vertex(&BTreeSet::from([0, 1])).incoming.is_empty());
        assert_eq!(
            graph.vertex(&BTreeSet::from([0, 1])).outgoing,
            vec![edge(2)]
        );
        assert_eq!(graph.vertex(&BTreeSet::from([2])).incoming, vec![edge(2)]);
    }

    #[test]
    fn thermal_detachable_cycles_keep_distribution_derivatives() {
        let mut parsed = crate::graph_io::test_graphs::box_graph();
        parsed.external_edges.clear();
        parsed.external_names.clear();
        for edge in &mut parsed.internal_edges {
            edge.signature.external_signature.clear();
            edge.mass_key = Some("m".to_string());
        }
        for (count, derivative_order, numerator, denominator) in
            [(1, 0, 1, 1), (2, 1, -1, 1), (3, 2, 1, 2), (4, 3, -1, 6)]
        {
            let mut cycle = parsed.clone();
            cycle.internal_edges.truncate(count);
            cycle.internal_edges[count - 1].head = 0;
            let chains = enumerate_cff_surface_chains(
                &cycle,
                &vec![1; count],
                MediumMode::ThermodynamicEquilibrium,
                &BTreeMap::new(),
            )
            .unwrap();
            assert_eq!(chains.len(), 1);
            assert!(chains[0].surfaces.is_empty());
            assert_eq!(chains[0].prefactor, Rational::new(numerator, denominator));
            assert_eq!(
                chains[0].thermal_weight.distributions,
                vec![ThermalDistributionFactor {
                    edge_id: EdgeIndex(0),
                    sign: i32::from(derivative_order != 0),
                    derivative_order,
                }]
            );
        }
    }

    #[test]
    fn thermal_cycles_with_distinct_and_repeated_masses_are_divided_differences() {
        for masses in [
            vec![0, 1],
            vec![0, 1, 2],
            vec![0, 1, 2, 3],
            vec![0, 0, 1],
            vec![0, 1, 0],
            vec![0, 0, 0],
        ] {
            for external_count in [0, 2, 4] {
                let count = masses.len();
                let mut parsed = crate::graph_io::test_graphs::box_graph();
                // Any number of external legs at the sole attachment must leave
                // the same cyclic-chain divided difference.
                parsed.external_edges = (0..external_count)
                    .map(|edge_id| crate::graph_io::ParsedGraphExternalEdge {
                        edge_id,
                        source: (edge_id % 2 == 0).then_some(0),
                        destination: (edge_id % 2 == 1).then_some(0),
                        label: format!("p{edge_id}"),
                        external_coefficients: Vec::new(),
                    })
                    .collect();
                parsed.external_names.clear();
                parsed.internal_edges.truncate(count);
                for (edge, mass) in parsed.internal_edges.iter_mut().zip(&masses) {
                    edge.head = (edge.edge_id + 1) % count;
                    edge.mass_key = Some(format!("m{mass}"));
                    edge.signature.external_signature.clear();
                }
                let signs = vec![1; count];
                assert!(
                    enumerate_cff_surface_chains(
                        &parsed,
                        &signs,
                        MediumMode::Vacuum,
                        &BTreeMap::new()
                    )
                    .unwrap()
                    .is_empty()
                );
                let expected = if masses == [0, 0, 0] {
                    parse!("d2n0/2")
                } else if masses.iter().filter(|&&mass| mass == 0).count() == 2 {
                    parse!("(n1-n0+(e0-e1)*dn0)/(e0-e1)^2")
                } else {
                    // The ordinary residue sum is independent of contraction order.
                    masses.iter().fold(Atom::Zero, |sum, mass| {
                        let energy = parse!(format!("e{mass}"));
                        let denominator = masses
                            .iter()
                            .filter(|other| *other != mass)
                            .fold(Atom::one(), |product, other| {
                                product * (&energy - parse!(format!("e{other}")))
                            });
                        let sign = if count.is_multiple_of(2) { -1 } else { 1 };
                        sum + Atom::num(sign) * parse!(format!("n{mass}")) / denominator
                    })
                };
                for mode in [
                    MediumMode::ThermodynamicEquilibrium,
                    MediumMode::ZeroTemperatureEquilibrium,
                ] {
                    let finite = mode.is_finite_temperature();
                    let chains =
                        enumerate_cff_surface_chains(&parsed, &signs, mode, &BTreeMap::new())
                            .unwrap();
                    assert!(!chains.is_empty(), "{masses:?}, {mode:?}");
                    let mut actual = chains.iter().fold(Atom::Zero, |sum, chain| {
                        let numerator = chain
                            .thermal_weight
                            .numerators
                            .iter()
                            .map(|numerator| numerator.to_atom(finite))
                            .chain(
                                chain
                                    .thermal_weight
                                    .distributions
                                    .iter()
                                    .map(|factor| factor.to_atom(finite)),
                            )
                            .fold(Atom::num(chain.prefactor.clone()), |product, factor| {
                                product * factor
                            });
                        let denominator =
                            chain.surfaces.iter().fold(Atom::one(), |product, surface| {
                                product * surface.to_atom(&[])
                            });
                        sum + numerator / denominator
                    });
                    for (edge_id, mass) in masses.iter().enumerate() {
                        for sign in [-1, 1] {
                            for derivative_order in 0..count {
                                let value = match derivative_order {
                                    0 => Atom::num((1 + sign) / 2) + parse!(format!("n{mass}")),
                                    1 => parse!(format!("dn{mass}")),
                                    order => parse!(format!("d{order}n{mass}")),
                                };
                                actual = actual
                                    .replace(
                                        ThermalDistributionFactor {
                                            edge_id: EdgeIndex(edge_id),
                                            sign,
                                            derivative_order,
                                        }
                                        .to_atom(finite),
                                    )
                                    .with(value);
                            }
                        }
                        actual = actual
                            .replace(LinearEnergyExpr::ose(EdgeIndex(edge_id), 1).to_atom(&[]))
                            .with(parse!(format!("e{mass}")));
                    }
                    assert_eq!(
                        (&actual - &expected).together(),
                        Atom::Zero,
                        "{masses:?}, {mode:?}: {actual} != {expected}",
                    );
                }
            }
        }
    }

    #[test]
    fn thermal_cut_self_edges_do_not_obstruct_cycle_reduction() {
        for masses in [["m", "m"], ["m0", "m1"]] {
            let mut parsed = crate::graph_io::test_graphs::box_graph();
            parsed.external_edges.clear();
            parsed.external_names = vec!["p".to_string()];
            parsed.internal_edges.truncate(2);
            for (edge, mass) in parsed.internal_edges.iter_mut().zip(masses) {
                edge.head = (edge.edge_id + 1) % 2;
                edge.mass_key = Some(mass.to_string());
                edge.signature.external_signature = vec![0];
            }
            let mode = MediumMode::ThermodynamicEquilibrium;
            let expected =
                enumerate_cff_surface_chains(&parsed, &[1; 2], mode, &BTreeMap::new()).unwrap();
            for node in 0..2 {
                let mut cut = parsed.internal_edges[0].clone();
                cut.edge_id = parsed.internal_edges.len();
                cut.tail = node;
                cut.head = node;
                cut.signature.loop_signature.fill(0);
                cut.signature.external_signature = vec![1];
                parsed.initial_state_cut_edges.push(
                    crate::graph_io::ParsedGraphInitialStateCutEdge {
                        edge_id: cut.edge_id,
                        external_id: 0,
                        external_sign: 1,
                    },
                );
                parsed.internal_edges.push(cut);
            }
            let actual =
                enumerate_cff_surface_chains(&parsed, &[1; 4], mode, &BTreeMap::new()).unwrap();
            assert_eq!(actual.len(), expected.len());
            for (actual, expected) in actual.iter().zip(expected) {
                assert_eq!(actual.surfaces, expected.surfaces);
                assert_eq!(actual.thermal_weight, expected.thermal_weight);
                assert_eq!(actual.prefactor, expected.prefactor);
            }
        }
    }

    #[test]
    fn thermal_cycle_stripping_skips_distinct_poles_and_finds_later_equal_cycles() {
        for (first_masses, first_shifts, expected_edges) in [
            (["m0", "m1"], [0, 0], vec![2]),
            (["m0", "m0"], [0, 1], vec![2]),
            (["m0", "m0"], [0, 0], vec![0, 2]),
        ] {
            let mut parsed = crate::graph_io::test_graphs::box_graph();
            parsed.external_edges.clear();
            for (edge, (tail, head)) in
                parsed
                    .internal_edges
                    .iter_mut()
                    .zip([(0, 1), (1, 0), (0, 2), (2, 0)])
            {
                edge.tail = tail;
                edge.head = head;
                edge.mass_key = Some("m2".to_string());
                edge.signature.external_signature.fill(0);
            }
            for (edge, (mass, shift)) in parsed
                .internal_edges
                .iter_mut()
                .zip(first_masses.into_iter().zip(first_shifts))
            {
                edge.mass_key = Some(mass.to_string());
                edge.signature.external_signature[0] = shift;
            }
            let mut graph = build_base_graph_from_parsed(&parsed);
            let (distributions, prefactor) = graph.strip_tadpoles_and_cyclic_chains(
                &parsed,
                MediumMode::ThermodynamicEquilibrium,
                &[],
                &BTreeMap::new(),
            );
            assert_eq!(
                distributions,
                expected_edges
                    .iter()
                    .map(|&edge_id| ThermalDistributionFactor {
                        edge_id: EdgeIndex(edge_id),
                        sign: 1,
                        derivative_order: 1,
                    })
                    .collect::<Vec<_>>(),
            );
            assert_eq!(
                prefactor,
                Rational::from(if expected_edges.len() == 1 { -1 } else { 1 }),
            );
            if expected_edges.len() == 1 {
                assert_eq!(graph.vertices.len(), 2);
                assert!(graph.has_directed_cycle());
            } else {
                assert!(graph.vertices.is_empty());
            }
        }
    }

    #[test]
    fn thermal_cycle_prefactors_compose_across_stripping_and_recursion() {
        for (edges, expected) in [
            (
                vec![(0, 1), (1, 0), (0, 2), (2, 0)],
                vec![(0, Rational::one(), vec![(0, 1), (2, 1)])],
            ),
            (
                vec![(0, 1), (1, 0), (0, 2)],
                vec![(1, Rational::from(-1), vec![(0, 1)])],
            ),
            (
                // Contracting 0 with 2 exposes a two-edge cycle below recursion.
                vec![(0, 1), (0, 2), (2, 0), (1, 2)],
                vec![
                    (1, Rational::from(-1), vec![(0, 1)]),
                    (2, Rational::one(), vec![]),
                ],
            ),
        ] {
            let mut parsed = crate::graph_io::test_graphs::box_graph();
            parsed.external_edges.clear();
            parsed.external_names.clear();
            parsed.internal_edges.truncate(edges.len());
            for (edge, &(tail, head)) in parsed.internal_edges.iter_mut().zip(&edges) {
                edge.tail = tail;
                edge.head = head;
                edge.signature.external_signature.clear();
                edge.mass_key = Some("m".to_string());
            }
            for mode in [
                MediumMode::ThermodynamicEquilibrium,
                MediumMode::ZeroTemperatureEquilibrium,
            ] {
                let mut actual = enumerate_cff_surface_chains(
                    &parsed,
                    &vec![1; edges.len()],
                    mode,
                    &BTreeMap::new(),
                )
                .unwrap()
                .into_iter()
                .map(|chain| {
                    (
                        chain.surfaces.len(),
                        chain.prefactor,
                        chain
                            .thermal_weight
                            .distributions
                            .iter()
                            .map(|factor| (factor.edge_id.0, factor.derivative_order))
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>();
                actual.sort();
                assert_eq!(actual, expected, "{edges:?}, {mode:?}");
            }
        }
    }

    #[test]
    fn thermal_cycles_with_multiple_attachments_are_not_stripped() {
        let parsed = crate::graph_io::test_graphs::box_graph();
        let mut graph = build_base_graph_from_parsed(&parsed);
        assert!(
            graph
                .strip_tadpoles_and_cyclic_chains(
                    &parsed,
                    MediumMode::ThermodynamicEquilibrium,
                    &[],
                    &BTreeMap::new()
                )
                .0
                .is_empty()
        );
        let chains = enumerate_cff_surface_chains(
            &parsed,
            &[1; 4],
            MediumMode::ThermodynamicEquilibrium,
            &BTreeMap::new(),
        )
        .unwrap();
        assert!(!chains.is_empty());
        assert!(chains.iter().all(|chain| !chain.surfaces.is_empty()));
    }

    #[test]
    fn thermal_unreducible_orientation_is_an_error() {
        let mut parsed = crate::graph_io::test_graphs::box_graph();
        // Inconsistent routing without the external attachments needed to
        // support it cannot be contracted or reduced as a common-momentum cycle.
        // It must not silently become an empty orientation contribution.
        parsed.external_edges.clear();
        let mode = MediumMode::ThermodynamicEquilibrium;
        let error = enumerate_cff_surface_chains(&parsed, &[1; 4], mode, &BTreeMap::new())
            .err()
            .expect("unreducible orientation must fail generation");
        assert!(matches!(
            error,
            GenerationError::NoAdmissibleCffVertex {
                medium_mode,
                edge_signs,
                vertices,
            } if medium_mode == mode
                && edge_signs == [1; 4]
                && vertices == vec![vec![0], vec![1], vec![2], vec![3]]
        ));

        // A cut-only endpoint has no virtual neighbour to contract. Selecting
        // it as a source/sink would otherwise silently emit no thermal chains.
        let mut parsed = crate::graph_io::test_graphs::initial_state_cut_line_graph(1);
        parsed
            .node_name_to_internal
            .insert("cut_only".to_string(), 2);
        parsed.internal_edges[0].tail = 2;
        assert!(matches!(
            enumerate_cff_surface_chains(&parsed, &[1; 2], mode, &BTreeMap::new()),
            Err(GenerationError::NoAdmissibleCffVertex { .. })
        ));
    }
}
