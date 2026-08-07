use std::collections::{BTreeMap, BTreeSet};

use ahash::AHashSet;
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeIndex, Flow, HedgePair},
    subgraph::{InternalSubGraph, ModifySubSet, SuBitGraph},
};
use spenso::structure::concrete_index::ExpandedIndex;
use symbolica::domains::rational::Rational;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    function,
    id::Replacement,
};
use three_dimensional_reps::{
    EnergyEdgeIndexMap, LinearSurface, MomentumSignature, ParsedGraph, ThreeDGraphSource,
    graph_io::{
        GraphIoError, ParsedGraphExternalEdge, ParsedGraphInitialStateCutEdge,
        ParsedGraphInternalEdge, initial_state_cut_external_alias,
    },
    utils::{rank_i64, solve_rational_system},
};

use crate::{
    graph::{Graph, LMBext},
    momentum::SignOrZero,
    utils::{GS, W_},
};

impl ThreeDGraphSource for Graph {
    fn to_three_d_parsed_graph(&self) -> three_dimensional_reps::graph_io::Result<ParsedGraph> {
        GraphThreeDSource::new(self, &[])?.to_three_d_parsed_graph()
    }

    fn energy_edge_index_map(&self, parsed: &ParsedGraph) -> Option<EnergyEdgeIndexMap> {
        GraphThreeDSource::new(self, &[])
            .ok()?
            .energy_edge_index_map(parsed)
    }
}

pub(crate) struct GraphThreeDSource<'a> {
    graph: &'a Graph,
    contract_edges: AHashSet<EdgeIndex>,
    initial_state_cut_edges: AHashSet<EdgeIndex>,
    inner_loop_count: usize,
    outer_loop_edges: Vec<EdgeIndex>,
    edge_loop_coordinates: BTreeMap<EdgeIndex, Vec<Rational>>,
    parent_loop_coordinates: Vec<Vec<Rational>>,
    exact_denominators: Option<&'a [FourDDenominator]>,
    exact_signatures: Vec<MomentumSignature>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct FourDDenominator {
    pub(crate) source_edge: EdgeIndex,
    pub(crate) momentum: Atom,
    pub(crate) mass_squared: Atom,
    pub(crate) full_expr: Atom,
}

impl<'a> GraphThreeDSource<'a> {
    pub(crate) fn new(
        graph: &'a Graph,
        contract_edges: &[EdgeIndex],
    ) -> three_dimensional_reps::graph_io::Result<Self> {
        let contract_edges = contract_edges.iter().copied().collect::<AHashSet<_>>();
        let initial_state_cut_edges = graph
            .iter_edges_of(&graph.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect::<AHashSet<_>>();
        let contracts_edge = |edge_id: EdgeIndex| {
            contract_edges.contains(&edge_id) && !initial_state_cut_edges.contains(&edge_id)
        };

        let mut contracted_filter: SuBitGraph = graph.empty_subgraph();
        for (pair, edge_id, edge_data) in graph.underlying.iter_edges() {
            if pair.is_paired() && !edge_data.data.is_dummy && contracts_edge(edge_id) {
                contracted_filter.add(pair);
            }
        }
        let contracted_subgraph =
            InternalSubGraph::cleaned_filter_optimist(contracted_filter, &graph.underlying);
        let inner_lmb = graph
            .try_compatible_sub_lmb(
                &contracted_subgraph,
                graph.dummy_stripped_external_flows_of(&contracted_subgraph),
                &graph.loop_momentum_basis,
            )
            .map_err(|error| {
                GraphIoError::Source(format!(
                    "contracted edges do not admit a parent-compatible inner loop basis: {error}"
                ))
            })?;

        let parent_loop_count = graph.loop_momentum_basis.loop_edges.len();
        let mut basis_rows = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .filter(|edge_id| inner_lmb.loop_edges.contains(edge_id))
            .map(|edge_id| {
                graph.loop_momentum_basis.edge_signatures[*edge_id]
                    .internal
                    .iter()
                    .map(|sign| sign_to_i32(*sign))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let inner_loop_count = basis_rows.len();
        if inner_loop_count != inner_lmb.loop_edges.len()
            || exact_integer_rank(&basis_rows) != inner_loop_count
        {
            return Err(GraphIoError::Source(format!(
                "parent-compatible inner loop basis has {} carriers but rank {} in the {}-dimensional production basis",
                inner_lmb.loop_edges.len(),
                exact_integer_rank(&basis_rows),
                parent_loop_count,
            )));
        }

        let mut outer_loop_edges = Vec::new();
        let mut basis_rank = inner_loop_count;
        for (pair, edge_id, edge_data) in graph
            .underlying
            .iter_edges()
            .sorted_by_key(|(_, edge_id, _)| *edge_id)
        {
            if !pair.is_paired()
                || edge_data.data.is_dummy
                || initial_state_cut_edges.contains(&edge_id)
                || contracts_edge(edge_id)
            {
                continue;
            }
            let row = graph.loop_momentum_basis.edge_signatures[edge_id]
                .internal
                .iter()
                .map(|sign| sign_to_i32(*sign))
                .collect::<Vec<_>>();
            let mut trial = basis_rows.clone();
            trial.push(row.clone());
            let trial_rank = exact_integer_rank(&trial);
            if trial_rank > basis_rank {
                basis_rows.push(row);
                outer_loop_edges.push(edge_id);
                basis_rank = trial_rank;
            }
            if basis_rank == parent_loop_count {
                break;
            }
        }
        if basis_rank != parent_loop_count {
            return Err(GraphIoError::Source(format!(
                "contracted source supplies loop-signature rank {basis_rank}, expected the production rank {parent_loop_count}"
            )));
        }

        let basis_matrix = (0..parent_loop_count)
            .map(|column| {
                basis_rows
                    .iter()
                    .map(|row| Rational::from(row[column]))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let solve_coordinates = |row: Vec<SignOrZero>| {
            if row.is_empty() {
                return Ok(Vec::new());
            }
            let rhs = row
                .iter()
                .map(|sign| Rational::from(sign_to_i32(*sign)))
                .collect::<Vec<_>>();
            solve_rational_system(basis_matrix.clone(), rhs).ok_or_else(|| {
                GraphIoError::Source(
                    "failed to solve exact contracted-source loop coordinates".to_string(),
                )
            })
        };
        let edge_loop_coordinates = graph
            .loop_momentum_basis
            .edge_signatures
            .iter()
            .map(|(edge_id, signature)| {
                Ok((
                    edge_id,
                    solve_coordinates(signature.internal.iter().copied().collect())?,
                ))
            })
            .collect::<three_dimensional_reps::graph_io::Result<BTreeMap<_, _>>>()?;
        for (pair, edge_id, edge_data) in graph.underlying.iter_edges() {
            if pair.is_paired()
                && !edge_data.data.is_dummy
                && !contracts_edge(edge_id)
                && edge_loop_coordinates[&edge_id]
                    .iter()
                    .take(inner_loop_count)
                    .any(|coordinate| *coordinate != 0)
            {
                return Err(GraphIoError::Source(format!(
                    "surviving edge {} retains an inner-loop coordinate after contraction",
                    usize::from(edge_id)
                )));
            }
        }
        let parent_loop_coordinates = (0..parent_loop_count)
            .map(|loop_index| {
                let row = (0..parent_loop_count)
                    .map(|column| {
                        if column == loop_index {
                            SignOrZero::Plus
                        } else {
                            SignOrZero::Zero
                        }
                    })
                    .collect::<Vec<_>>();
                solve_coordinates(row)
            })
            .collect::<three_dimensional_reps::graph_io::Result<Vec<_>>>()?;

        Ok(Self {
            graph,
            contract_edges,
            initial_state_cut_edges,
            inner_loop_count,
            outer_loop_edges,
            edge_loop_coordinates,
            parent_loop_coordinates,
            exact_denominators: None,
            exact_signatures: Vec::new(),
        })
    }

    pub(crate) fn from_exact_denominators(
        graph: &'a Graph,
        denominators: &'a [FourDDenominator],
    ) -> three_dimensional_reps::graph_io::Result<Self> {
        let parent_signatures = denominators
            .iter()
            .map(|denominator| denominator.momentum_signature(graph))
            .collect::<three_dimensional_reps::graph_io::Result<Vec<_>>>()?;
        // Exact 4D sources must use the same coordinate frame as the ordinary
        // production CFF. The configured graph LMB is not necessarily the
        // canonical carrier basis selected by `GraphThreeDSource::new`.
        let production_source = Self::new(graph, &[])?;
        let parent_loop_count = graph.loop_momentum_basis.loop_edges.len();
        let integral_coordinate = |coordinate: &Rational| {
            (coordinate.denominator_ref().to_i64() == Some(1))
                .then(|| coordinate.numerator_ref().to_i64())
                .flatten()
                .and_then(|numerator| i32::try_from(numerator).ok())
        };
        let production_loop_rows = parent_signatures
            .iter()
            .map(|signature| {
                (0..parent_loop_count)
                    .map(|production_column| {
                        let coordinate = signature
                            .loop_signature
                            .iter()
                            .zip(&production_source.parent_loop_coordinates)
                            .fold(Rational::from(0), |sum, (coefficient, parent_coordinates)| {
                                sum + Rational::from(*coefficient)
                                    * parent_coordinates[production_column].clone()
                            });
                        integral_coordinate(&coordinate).ok_or_else(|| {
                            GraphIoError::Source(format!(
                                "exact 4D denominator has non-integral production CFF coordinate {coordinate}"
                            ))
                        })
                    })
                    .collect()
            })
            .collect::<three_dimensional_reps::graph_io::Result<Vec<Vec<_>>>>()?;
        let (outer_basis_rows, active_loop_columns) =
            select_exact_loop_basis(&production_loop_rows, parent_loop_count)?;
        let inactive_loop_columns = (0..parent_loop_count)
            .filter(|column| !active_loop_columns.contains(column))
            .collect::<Vec<_>>();
        let mut source_basis_rows = inactive_loop_columns
            .iter()
            .map(|inactive_column| {
                (0..parent_loop_count)
                    .map(|column| i32::from(column == *inactive_column))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        source_basis_rows.extend(outer_basis_rows.iter().cloned());
        let basis_matrix = (0..parent_loop_count)
            .map(|column| {
                source_basis_rows
                    .iter()
                    .map(|row| Rational::from(row[column]))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let solve_coordinates = |row: &[Rational]| {
            if row.is_empty() {
                return Ok(Vec::new());
            }
            solve_rational_system(basis_matrix.clone(), row.to_vec()).ok_or_else(|| {
                GraphIoError::Source(
                    "failed to solve canonical exact-source loop coordinates".to_string(),
                )
            })
        };
        let edge_loop_coordinates = graph
            .loop_momentum_basis
            .edge_signatures
            .iter()
            .map(|(edge_id, _)| {
                let production_coordinates = production_source
                    .edge_loop_coordinates
                    .get(&edge_id)
                    .expect("production source contains every graph edge");
                Ok((edge_id, solve_coordinates(production_coordinates)?))
            })
            .collect::<three_dimensional_reps::graph_io::Result<BTreeMap<_, _>>>()?;
        let parent_loop_coordinates = production_source
            .parent_loop_coordinates
            .iter()
            .map(|coordinates| solve_coordinates(coordinates))
            .collect::<three_dimensional_reps::graph_io::Result<Vec<_>>>()?;
        let inner_loop_count = inactive_loop_columns.len();
        let exact_signatures = parent_signatures
            .iter()
            .zip(&production_loop_rows)
            .map(|(signature, production_row)| {
                let coordinates = solve_coordinates(
                    &production_row
                        .iter()
                        .map(|coefficient| Rational::from(*coefficient))
                        .collect::<Vec<_>>(),
                )?;
                if coordinates[..inner_loop_count]
                    .iter()
                    .any(|coordinate| *coordinate != 0)
                {
                    return Err(GraphIoError::Source(
                        "exact 4D denominator retains an inactive-loop coordinate".to_string(),
                    ));
                }
                let loop_signature = coordinates[inner_loop_count..]
                    .iter()
                    .map(|coordinate| {
                        integral_coordinate(coordinate).ok_or_else(|| {
                            GraphIoError::Source(format!(
                                "exact 4D denominator has non-integral source coordinate {coordinate}"
                            ))
                        })
                    })
                    .collect::<three_dimensional_reps::graph_io::Result<Vec<_>>>()?;
                Ok(MomentumSignature {
                    loop_signature,
                    external_signature: signature.external_signature.clone(),
                })
            })
            .collect::<three_dimensional_reps::graph_io::Result<Vec<_>>>()?;
        let outer_loop_edges = outer_basis_rows
            .iter()
            .map(|basis_row| {
                production_loop_rows
                    .iter()
                    .zip(denominators)
                    .filter_map(|(row, denominator)| {
                        (row == basis_row
                            || row
                                .iter()
                                .zip(basis_row)
                                .all(|(coefficient, basis_coefficient)| {
                                    coefficient.checked_neg() == Some(*basis_coefficient)
                                }))
                        .then_some(denominator.source_edge)
                    })
                    .min()
                    .expect("every exact source basis row comes from a denominator")
            })
            .collect();
        let active_cograph_edges = denominators
            .iter()
            .filter(|denominator| !denominator.uses_uv_loop_basis())
            .map(|denominator| denominator.source_edge)
            .collect::<AHashSet<_>>();
        let initial_state_cut_edges = graph
            .iter_edges_of(&graph.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect::<AHashSet<_>>();
        for edge_id in &initial_state_cut_edges {
            if edge_loop_coordinates[edge_id][..inner_loop_count]
                .iter()
                .any(|coordinate| *coordinate != 0)
            {
                return Err(GraphIoError::Source(format!(
                    "initial-state cut carrier edge {} retains an inactive-loop coordinate",
                    usize::from(*edge_id)
                )));
            }
        }
        let contract_edges = graph
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge_id, edge)| {
                (pair.is_paired()
                    && !edge.data.is_dummy
                    && !active_cograph_edges.contains(&edge_id)
                    && !initial_state_cut_edges.contains(&edge_id))
                .then_some(edge_id)
            })
            .collect();

        Ok(Self {
            graph,
            contract_edges,
            initial_state_cut_edges,
            inner_loop_count,
            outer_loop_edges,
            edge_loop_coordinates,
            parent_loop_coordinates,
            exact_denominators: Some(denominators),
            exact_signatures,
        })
    }

    pub(crate) fn edge_loop_coordinates(&self, edge_id: EdgeIndex) -> Option<&[Rational]> {
        self.edge_loop_coordinates.get(&edge_id).map(Vec::as_slice)
    }

    pub(crate) fn reconstructible_outer_loop_coordinates(
        &self,
        edge_id: EdgeIndex,
    ) -> Option<&[Rational]> {
        let coordinates = self.edge_loop_coordinates(edge_id)?;
        coordinates[..self.inner_loop_count]
            .iter()
            .all(|coordinate| *coordinate == 0)
            .then_some(&coordinates[self.inner_loop_count..])
    }

    fn outer_loop_signature(
        &self,
        edge_id: EdgeIndex,
    ) -> three_dimensional_reps::graph_io::Result<Vec<i32>> {
        self.edge_loop_coordinates(edge_id)
            .ok_or_else(|| {
                GraphIoError::Source(format!(
                    "missing contracted-source coordinates for edge {}",
                    usize::from(edge_id)
                ))
            })?
            .iter()
            .skip(self.inner_loop_count)
            .map(|coordinate| {
                let numerator = coordinate.numerator_ref().to_i64().ok_or_else(|| {
                    GraphIoError::Source(format!(
                        "contracted-source coordinate {coordinate} is outside the i64 range"
                    ))
                })?;
                let denominator = coordinate.denominator_ref().to_i64().ok_or_else(|| {
                    GraphIoError::Source(format!(
                        "contracted-source coordinate {coordinate} has an out-of-range denominator"
                    ))
                })?;
                if denominator != 1 {
                    return Err(GraphIoError::Source(format!(
                        "contracted-source coordinate {coordinate} is not integral"
                    )));
                }
                i32::try_from(numerator).map_err(|_| {
                    GraphIoError::Source(format!(
                        "contracted-source coordinate {coordinate} is outside the i32 range"
                    ))
                })
            })
            .collect()
    }

    fn contracts_edge(&self, edge_id: EdgeIndex) -> bool {
        self.contract_edges.contains(&edge_id) && !self.initial_state_cut_edges.contains(&edge_id)
    }

    pub(crate) fn physical_energy_edge_index_map(&self) -> Option<EnergyEdgeIndexMap> {
        let denominators = self.exact_denominators?;
        let parsed = self.exact_parsed_graph().ok()?;
        let offset = self.graph.underlying.n_edges();
        let mut internal = denominators
            .iter()
            .enumerate()
            .map(|(local, denominator)| (offset + local, usize::from(denominator.source_edge)))
            .collect::<BTreeMap<_, _>>();
        for (cut_edge, source_edge) in parsed
            .initial_state_cut_edges
            .iter()
            .zip(self.initial_state_cut_edges.iter().copied().sorted())
        {
            internal.insert(offset + cut_edge.edge_id, usize::from(source_edge));
        }
        Some(EnergyEdgeIndexMap {
            internal,
            // Exact-source generation already remaps external coordinates to
            // physical graph edge IDs. This second map changes only the
            // occurrence-local internal-energy namespace.
            external: BTreeMap::new(),
            orientation_edge_count: self.graph.underlying.n_edges(),
        })
    }

    pub(crate) fn physical_linear_surface(&self, surface: &LinearSurface) -> Option<LinearSurface> {
        let denominators = self.exact_denominators?;
        let edge_map = self.physical_energy_edge_index_map()?;
        let offset = self.graph.underlying.n_edges();
        let all_physical = surface
            .expression
            .internal_terms
            .iter()
            .all(|(edge_id, _)| {
                let local = usize::from(*edge_id);
                if let Some(denominator) = local
                    .checked_sub(offset)
                    .and_then(|local| denominators.get(local))
                {
                    denominator.is_original_graph_denominator(self.graph)
                } else {
                    edge_map.internal.contains_key(&local)
                }
            });
        all_physical.then(|| {
            let mut surface = surface.clone();
            surface.expression = surface
                .expression
                .remap_energy_edges(&edge_map.internal, &BTreeMap::new());
            surface
        })
    }

    pub(crate) fn exact_ose_replacements(&self) -> Option<Vec<Replacement>> {
        let offset = self.graph.underlying.n_edges();
        let mut replacements = self
            .exact_denominators?
            .iter()
            .enumerate()
            .map(|(local, denominator)| {
                Replacement::new(
                    crate::utils::ose_atom_from_index(EdgeIndex(offset + local)).to_pattern(),
                    denominator.on_shell_energy().to_pattern(),
                )
            })
            .collect::<Vec<_>>();
        let parsed = self.exact_parsed_graph().ok()?;
        replacements.extend(
            parsed
                .initial_state_cut_edges
                .iter()
                .zip(self.initial_state_cut_edges.iter().copied().sorted())
                .map(|(cut_edge, source_edge)| {
                    Replacement::new(
                        crate::utils::ose_atom_from_index(EdgeIndex(offset + cut_edge.edge_id))
                            .to_pattern(),
                        crate::utils::ose_atom_from_index(source_edge).to_pattern(),
                    )
                }),
        );
        Some(replacements)
    }

    pub(crate) fn exact_inverse_energy_product(&self) -> Option<Atom> {
        Some(
            self.exact_denominators?
                .iter()
                .map(|denominator| -Atom::num(2) * denominator.on_shell_energy())
                .reduce(|product, energy| product * energy)
                .map(|product| Atom::one() / product)
                .unwrap_or_else(Atom::one),
        )
    }

    pub(crate) fn active_loop_count(&self) -> usize {
        self.outer_loop_edges.len()
    }

    pub(crate) fn contract_subgraph(&self) -> SuBitGraph {
        let mut subgraph: SuBitGraph = self.graph.empty_subgraph();
        for (pair, edge_id, _) in self.graph.underlying.iter_edges() {
            if pair.is_paired() && self.contracts_edge(edge_id) {
                subgraph.add(pair);
            }
        }
        subgraph
    }

    fn exact_parsed_graph(&self) -> three_dimensional_reps::graph_io::Result<ParsedGraph> {
        let denominators = self
            .exact_denominators
            .expect("exact parsed source has exact denominators");

        // Denominators absent from this additive 4D term are contracted only
        // in the source topology. The typed atom remains untouched and keeps
        // exact multiplicity and denominator expressions.
        let mut parent = (0..self.graph.n_nodes()).collect_vec();
        for (pair, edge_id, edge_data) in self.graph.underlying.iter_edges() {
            if edge_data.data.is_dummy || !self.contracts_edge(edge_id) {
                continue;
            }
            if let HedgePair::Paired { source, sink } = pair {
                union_parent(
                    &mut parent,
                    usize::from(self.graph.node_id(source)),
                    usize::from(self.graph.node_id(sink)),
                );
            }
        }

        let node_ids = self
            .graph
            .underlying
            .iter_nodes()
            .map(|(node_id, _, _)| node_id)
            .sorted()
            .collect::<Vec<_>>();
        let mut root_to_internal = BTreeMap::<usize, usize>::new();
        for node_id in &node_ids {
            let root = find_parent(&mut parent, usize::from(*node_id));
            if !root_to_internal.contains_key(&root) {
                root_to_internal.insert(root, root_to_internal.len());
            }
        }
        let node_to_internal = node_ids
            .iter()
            .map(|node_id| {
                let root = find_parent(&mut parent, usize::from(*node_id));
                (*node_id, root_to_internal[&root])
            })
            .collect::<BTreeMap<_, _>>();
        let uv_nodes = denominators
            .iter()
            .filter(|denominator| denominator.uses_uv_loop_basis())
            .flat_map(|denominator| {
                let (_, pair) = self.graph[&denominator.source_edge];
                match pair {
                    HedgePair::Paired { source, sink } => {
                        vec![self.graph.node_id(source), self.graph.node_id(sink)]
                    }
                    _ => Vec::new(),
                }
            })
            .collect::<BTreeSet<_>>();
        let uv_node_to_internal = uv_nodes
            .into_iter()
            .enumerate()
            .map(|(offset, node)| (node, root_to_internal.len() + offset))
            .collect::<BTreeMap<_, _>>();

        let mut internal_edges = Vec::new();
        for (denominator, signature) in denominators.iter().zip(&self.exact_signatures) {
            let (_, pair) = self.graph[&denominator.source_edge];
            let HedgePair::Paired { source, sink } = pair else {
                return Err(GraphIoError::Source(format!(
                    "4D denominator edge {} is not a paired internal edge",
                    usize::from(denominator.source_edge)
                )));
            };
            let node_map = if denominator.uses_uv_loop_basis() {
                &uv_node_to_internal
            } else {
                &node_to_internal
            };
            let local_edge_id = internal_edges.len();
            internal_edges.push(ParsedGraphInternalEdge {
                edge_id: local_edge_id,
                tail: node_map[&self.graph.node_id(source)],
                head: node_map[&self.graph.node_id(sink)],
                label: self.graph.underlying[denominator.source_edge]
                    .name
                    .value
                    .clone(),
                mass_key: Some(denominator.mass_squared.to_canonical_string()),
                signature: signature.clone(),
                had_pow: false,
            });
        }

        let mut initial_state_cut_edges = Vec::new();
        for edge_index in self.initial_state_cut_edges.iter().copied().sorted() {
            let (_, pair) = self.graph[&edge_index];
            let HedgePair::Paired { source, sink } = pair else {
                return Err(GraphIoError::Source(format!(
                    "initial-state cut carrier edge {} is not paired",
                    usize::from(edge_index)
                )));
            };
            let signature = &self.graph.loop_momentum_basis.edge_signatures[edge_index];
            let local_edge_id = internal_edges.len();
            internal_edges.push(ParsedGraphInternalEdge {
                edge_id: local_edge_id,
                tail: node_to_internal[&self.graph.node_id(source)],
                head: node_to_internal[&self.graph.node_id(sink)],
                label: self.graph.underlying[edge_index].name.value.clone(),
                mass_key: Some(
                    self.graph.underlying[edge_index]
                        .particle
                        .mass_atom()
                        .to_canonical_string(),
                ),
                signature: MomentumSignature {
                    loop_signature: self.outer_loop_signature(edge_index)?,
                    external_signature: (&signature.external)
                        .into_iter()
                        .map(sign_to_i32)
                        .collect(),
                },
                had_pow: false,
            });
            let (external_id, external_sign) = initial_state_cut_external_alias(
                usize::from(edge_index),
                &internal_edges[local_edge_id].signature,
            )?;
            initial_state_cut_edges.push(ParsedGraphInitialStateCutEdge {
                edge_id: local_edge_id,
                external_id,
                external_sign,
            });
        }

        let mut external_edges = self
            .graph
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge_id, edge)| match pair {
                HedgePair::Unpaired { hedge, flow } if !edge.data.is_dummy => {
                    let node = node_to_internal[&self.graph.node_id(hedge)];
                    let (source, destination) = match flow {
                        Flow::Source => (Some(node), None),
                        Flow::Sink => (None, Some(node)),
                    };
                    Some(ParsedGraphExternalEdge {
                        edge_id: usize::from(edge_id),
                        source,
                        destination,
                        label: edge.data.name.value.clone(),
                        external_coefficients: (&self.graph.loop_momentum_basis.edge_signatures
                            [edge_id]
                            .external)
                            .into_iter()
                            .map(sign_to_i32)
                            .collect(),
                    })
                }
                _ => None,
            })
            .collect::<Vec<_>>();
        complete_external_balance_edges(
            &mut external_edges,
            &internal_edges,
            &initial_state_cut_edges,
            &self
                .graph
                .loop_momentum_basis
                .ext_edges
                .iter()
                .map(|edge_id| self.graph.underlying[*edge_id].name.value.clone())
                .collect::<Vec<_>>(),
        );

        Ok(ParsedGraph {
            internal_edges,
            external_edges,
            initial_state_cut_edges,
            loop_names: self
                .outer_loop_edges
                .iter()
                .enumerate()
                .map(|(source_index, edge_id)| {
                    format!(
                        "{}#exact{source_index}",
                        self.graph.underlying[*edge_id].name.value
                    )
                })
                .collect(),
            external_names: self
                .graph
                .loop_momentum_basis
                .ext_edges
                .iter()
                .map(|edge_id| self.graph.underlying[*edge_id].name.value.clone())
                .collect(),
            node_name_to_internal: root_to_internal
                .into_iter()
                .map(|(root, node)| (format!("n{root}"), node))
                .chain(
                    uv_node_to_internal
                        .into_iter()
                        .map(|(node, internal)| (format!("uv{}", usize::from(node)), internal)),
                )
                .collect(),
        })
    }
}

impl FourDDenominator {
    fn uses_uv_loop_basis(&self) -> bool {
        self.mass_squared == Atom::var(GS.m_uv_expansion).pow(2)
    }

    fn is_original_graph_denominator(&self, graph: &Graph) -> bool {
        self.momentum
            == FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(self.source_edge))
                .finish()
            && self.mass_squared
                == graph.underlying[self.source_edge]
                    .particle
                    .mass_atom()
                    .pow(2)
    }

    fn on_shell_energy(&self) -> Atom {
        let spatial_square = (1..=3).fold(Atom::Zero, |sum, spatial_index| {
            let component_index = Atom::from(ExpandedIndex::from_iter([spatial_index]));
            let component = self
                .momentum
                .replace(function!(GS.emr_mom, W_.i_))
                .with(function!(GS.emr_mom, W_.i_, component_index));
            sum + component.clone() * component
        });
        (&self.mass_squared + spatial_square).sqrt()
    }

    pub(crate) fn momentum_signature(
        &self,
        graph: &Graph,
    ) -> three_dimensional_reps::graph_io::Result<MomentumSignature> {
        let mut loop_signature = vec![0; graph.loop_momentum_basis.loop_edges.len()];
        let mut external_signature = vec![0; graph.loop_momentum_basis.ext_edges.len()];
        accumulate_momentum_signature(
            graph,
            self.momentum.as_view(),
            1,
            self.uses_uv_loop_basis(),
            &mut loop_signature,
            &mut external_signature,
        )
        .map_err(|error| {
            GraphIoError::Source(format!(
                "could not extract 4D denominator momentum signature for edge {} from `{}`: {error}",
                usize::from(self.source_edge),
                self.momentum
            ))
        })?;
        Ok(MomentumSignature {
            loop_signature,
            external_signature,
        })
    }

    pub(crate) fn depends_on_loop(
        &self,
        graph: &Graph,
    ) -> three_dimensional_reps::graph_io::Result<bool> {
        Ok(self
            .momentum_signature(graph)?
            .loop_signature
            .into_iter()
            .any(|coefficient| coefficient != 0))
    }
}

impl ThreeDGraphSource for GraphThreeDSource<'_> {
    fn to_three_d_parsed_graph(&self) -> three_dimensional_reps::graph_io::Result<ParsedGraph> {
        if self.exact_denominators.is_some() {
            return self.exact_parsed_graph();
        }
        // Keep contraction virtual here: linnet graph contraction mutates/deletes
        // edges, while the CFF expression must still map back to the original
        // EdgeIndex values used by GammaLoop's energy and surface caches.
        let mut parent = (0..self.graph.n_nodes()).collect_vec();
        for (pair, edge_id, edge_data) in self.graph.underlying.iter_edges() {
            if edge_data.data.is_dummy || !self.contracts_edge(edge_id) {
                continue;
            }
            if let HedgePair::Paired { source, sink } = pair {
                union_parent(
                    &mut parent,
                    usize::from(self.graph.node_id(source)),
                    usize::from(self.graph.node_id(sink)),
                );
            }
        }

        let node_ids = self
            .graph
            .underlying
            .iter_nodes()
            .map(|(node_id, _, _)| node_id)
            .sorted()
            .collect::<Vec<_>>();
        let mut root_to_internal = BTreeMap::<usize, usize>::new();
        for node_id in &node_ids {
            let root = find_parent(&mut parent, usize::from(*node_id));
            if !root_to_internal.contains_key(&root) {
                root_to_internal.insert(root, root_to_internal.len());
            }
        }
        let mut node_to_internal = BTreeMap::new();
        for node_id in &node_ids {
            let root = find_parent(&mut parent, usize::from(*node_id));
            node_to_internal.insert(*node_id, root_to_internal[&root]);
        }

        let loop_names = self
            .outer_loop_edges
            .iter()
            .map(|edge_id| self.graph.underlying[*edge_id].name.value.clone())
            .collect::<Vec<_>>();
        let external_names = self
            .graph
            .loop_momentum_basis
            .ext_edges
            .iter()
            .map(|edge_id| self.graph.underlying[*edge_id].name.value.clone())
            .collect::<Vec<_>>();

        let mut internal_edges = Vec::new();
        let mut external_edges = Vec::new();
        let mut initial_state_cut_edges = Vec::new();
        let mut next_external_id = 10_000_000usize;

        for (pair, edge_index, edge_data) in self
            .graph
            .underlying
            .iter_edges()
            .sorted_by_key(|(_, edge_index, _)| *edge_index)
        {
            if edge_data.data.is_dummy {
                continue;
            }
            let signature = &self.graph.loop_momentum_basis.edge_signatures[edge_index];
            let momentum_signature = MomentumSignature {
                loop_signature: self.outer_loop_signature(edge_index)?,
                external_signature: (&signature.external).into_iter().map(sign_to_i32).collect(),
            };
            let label = edge_data.data.name.value.clone();
            match pair {
                HedgePair::Paired { source, sink } => {
                    if self.contracts_edge(edge_index) {
                        continue;
                    }
                    let local_edge_id = internal_edges.len();
                    let tail = *node_to_internal
                        .get(&self.graph.node_id(source))
                        .ok_or_else(|| {
                            GraphIoError::Source(format!(
                                "missing contracted source node mapping for edge {}",
                                usize::from(edge_index)
                            ))
                        })?;
                    let head =
                        *node_to_internal
                            .get(&self.graph.node_id(sink))
                            .ok_or_else(|| {
                                GraphIoError::Source(format!(
                                    "missing contracted sink node mapping for edge {}",
                                    usize::from(edge_index)
                                ))
                            })?;
                    internal_edges.push(ParsedGraphInternalEdge {
                        edge_id: local_edge_id,
                        tail,
                        head,
                        label,
                        mass_key: Some(edge_data.data.particle.mass_atom().to_canonical_string()),
                        signature: momentum_signature,
                        had_pow: false,
                    });
                    if self.initial_state_cut_edges.contains(&edge_index) {
                        let (external_id, external_sign) = initial_state_cut_external_alias(
                            usize::from(edge_index),
                            &internal_edges[local_edge_id].signature,
                        )?;
                        initial_state_cut_edges.push(ParsedGraphInitialStateCutEdge {
                            edge_id: local_edge_id,
                            external_id,
                            external_sign,
                        });
                    }
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let node = *node_to_internal
                        .get(&self.graph.node_id(hedge))
                        .ok_or_else(|| {
                            GraphIoError::Source(format!(
                                "missing contracted external node mapping for edge {}",
                                usize::from(edge_index)
                            ))
                        })?;
                    let (source, destination) = match flow {
                        Flow::Source => (Some(node), None),
                        Flow::Sink => (None, Some(node)),
                    };
                    external_edges.push(ParsedGraphExternalEdge {
                        edge_id: next_external_id,
                        source,
                        destination,
                        label,
                        external_coefficients: momentum_signature.external_signature,
                    });
                    next_external_id += 1;
                }
                HedgePair::Split { .. } => {
                    return Err(GraphIoError::Source(
                        "split edges are not supported when extracting GammaLoop Graph input"
                            .to_string(),
                    ));
                }
            }
        }

        Ok(ParsedGraph {
            internal_edges,
            external_edges,
            initial_state_cut_edges,
            loop_names,
            external_names,
            node_name_to_internal: root_to_internal
                .into_iter()
                .map(|(root, node)| (format!("n{root}"), node))
                .collect(),
        })
    }

    fn energy_edge_index_map(&self, _parsed: &ParsedGraph) -> Option<EnergyEdgeIndexMap> {
        if self.exact_denominators.is_some() {
            let internal_count = self.exact_parsed_graph().ok()?.internal_edges.len();
            let offset = self.graph.underlying.n_edges();
            return Some(EnergyEdgeIndexMap {
                internal: (0..internal_count)
                    .map(|edge_id| (edge_id, offset + edge_id))
                    .collect(),
                external: self
                    .graph
                    .loop_momentum_basis
                    .ext_edges
                    .iter()
                    .enumerate()
                    .map(|(external_id, edge_id)| (external_id, usize::from(*edge_id)))
                    .collect(),
                orientation_edge_count: offset + internal_count,
            });
        }
        let internal = self
            .graph
            .underlying
            .iter_edges()
            .sorted_by_key(|(_, edge_index, _)| *edge_index)
            .filter_map(|(pair, edge_index, edge_data)| {
                (pair.is_paired() && !edge_data.data.is_dummy && !self.contracts_edge(edge_index))
                    .then_some(usize::from(edge_index))
            })
            .enumerate()
            .collect::<BTreeMap<_, _>>();

        let external = self
            .graph
            .loop_momentum_basis
            .ext_edges
            .iter()
            .enumerate()
            .map(|(external_id, edge_id)| (external_id, usize::from(*edge_id)))
            .collect::<BTreeMap<_, _>>();

        Some(EnergyEdgeIndexMap {
            internal,
            external,
            orientation_edge_count: self.graph.underlying.n_edges(),
        })
    }
}

fn find_parent(parent: &mut [usize], node: usize) -> usize {
    let parent_node = parent[node];
    if parent_node == node {
        node
    } else {
        let root = find_parent(parent, parent_node);
        parent[node] = root;
        root
    }
}

fn union_parent(parent: &mut [usize], left: usize, right: usize) {
    let left_root = find_parent(parent, left);
    let right_root = find_parent(parent, right);
    if left_root != right_root {
        parent[right_root] = left_root;
    }
}

fn complete_external_balance_edges(
    external_edges: &mut Vec<ParsedGraphExternalEdge>,
    internal_edges: &[ParsedGraphInternalEdge],
    initial_state_cut_edges: &[ParsedGraphInitialStateCutEdge],
    external_names: &[String],
) {
    let external_count = external_names.len();
    let initial_state_cut_edge_ids = initial_state_cut_edges
        .iter()
        .map(|edge| edge.edge_id)
        .collect::<BTreeSet<_>>();
    let initial_state_external_ids = initial_state_cut_edges
        .iter()
        .map(|edge| edge.external_id)
        .collect::<BTreeSet<_>>();
    let mut balances = BTreeMap::<usize, Vec<i32>>::new();
    for edge in internal_edges {
        if initial_state_cut_edge_ids.contains(&edge.edge_id) {
            continue;
        }
        balances
            .entry(edge.tail)
            .or_insert_with(|| vec![0; external_count]);
        balances
            .entry(edge.head)
            .or_insert_with(|| vec![0; external_count]);
        for (external_id, coefficient) in edge.signature.external_signature.iter().enumerate() {
            if initial_state_external_ids.contains(&external_id) {
                continue;
            }
            balances.get_mut(&edge.tail).unwrap()[external_id] -= coefficient;
            balances.get_mut(&edge.head).unwrap()[external_id] += coefficient;
        }
    }
    for edge in external_edges.iter() {
        for node in [edge.source, edge.destination].into_iter().flatten() {
            balances
                .entry(node)
                .or_insert_with(|| vec![0; external_count]);
            for (external_id, coefficient) in edge.external_coefficients.iter().enumerate() {
                if !initial_state_external_ids.contains(&external_id) {
                    balances.get_mut(&node).unwrap()[external_id] += coefficient;
                }
            }
        }
    }

    let mut next_external_id = external_edges
        .iter()
        .map(|edge| edge.edge_id)
        .max()
        .map(|edge_id| edge_id + 1)
        .unwrap_or(0);
    for (node, balance) in balances {
        for (external_id, coefficient) in balance.into_iter().enumerate() {
            if initial_state_external_ids.contains(&external_id) {
                continue;
            }
            let name = external_names
                .get(external_id)
                .cloned()
                .unwrap_or_else(|| format!("p{external_id}"));
            for _ in 0..coefficient.max(0) {
                let mut external_coefficients = vec![0; external_count];
                external_coefficients[external_id] = -1;
                external_edges.push(ParsedGraphExternalEdge {
                    edge_id: next_external_id,
                    source: Some(node),
                    destination: None,
                    label: format!("-{name}"),
                    external_coefficients,
                });
                next_external_id += 1;
            }
            for _ in 0..(-coefficient).max(0) {
                let mut external_coefficients = vec![0; external_count];
                external_coefficients[external_id] = 1;
                external_edges.push(ParsedGraphExternalEdge {
                    edge_id: next_external_id,
                    source: None,
                    destination: Some(node),
                    label: name.clone(),
                    external_coefficients,
                });
                next_external_id += 1;
            }
        }
    }
}

fn select_exact_loop_basis(
    rows: &[Vec<i32>],
    loop_count: usize,
) -> three_dimensional_reps::graph_io::Result<(Vec<Vec<i32>>, Vec<usize>)> {
    if rows.iter().any(|row| row.len() != loop_count) {
        return Err(GraphIoError::Source(
            "exact 4D denominator loop signatures have inconsistent dimensions".to_string(),
        ));
    }
    let active_rank = exact_integer_rank(rows);
    if active_rank == 0 {
        return Ok((Vec::new(), Vec::new()));
    }

    // Prefer the production carriers, but only accept denominator rows whose
    // active minor keeps all rewritten momentum signatures integral.
    let mut candidates = rows
        .iter()
        .filter(|row| row.iter().any(|coefficient| *coefficient != 0))
        .map(|row| {
            let mut canonical = row.clone();
            if canonical
                .iter()
                .find(|coefficient| **coefficient != 0)
                .is_some_and(|coefficient| *coefficient < 0)
            {
                for coefficient in &mut canonical {
                    *coefficient = coefficient.checked_neg().ok_or_else(|| {
                        GraphIoError::Source(
                            "exact 4D denominator loop coefficient cannot be canonicalized"
                                .to_string(),
                        )
                    })?;
                }
            }
            Ok(canonical)
        })
        .collect::<three_dimensional_reps::graph_io::Result<BTreeSet<_>>>()?
        .into_iter()
        .collect::<Vec<_>>();
    let unit_axis = |row: &[i32]| {
        row.iter().enumerate().find_map(|(axis, coefficient)| {
            (*coefficient == 1
                && row
                    .iter()
                    .enumerate()
                    .all(|(column, value)| column == axis || *value == 0))
            .then_some(axis)
        })
    };
    candidates.sort_by(|left, right| match (unit_axis(left), unit_axis(right)) {
        (Some(left), Some(right)) => left.cmp(&right),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => left.cmp(right),
    });

    for selected in candidates.iter().cloned().combinations(active_rank) {
        for columns in (0..loop_count).combinations(active_rank) {
            let matrix = selected
                .iter()
                .map(|row| {
                    columns
                        .iter()
                        .map(|column| Rational::from(row[*column]))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let inverse_is_integral = (0..active_rank).all(|rhs_index| {
                let rhs = (0..active_rank)
                    .map(|index| Rational::from(i32::from(index == rhs_index)))
                    .collect::<Vec<_>>();
                solve_rational_system(matrix.clone(), rhs).is_some_and(|solution| {
                    solution
                        .iter()
                        .all(|coordinate| coordinate.denominator_ref().to_i64() == Some(1))
                })
            });
            if inverse_is_integral {
                return Ok((selected, columns));
            }
        }
    }

    Err(GraphIoError::Source(format!(
        "exact 4D source rank {active_rank} has no unimodular denominator-carrier basis in the {loop_count}-dimensional production CFF coordinates"
    )))
}

fn accumulate_momentum_signature(
    graph: &Graph,
    view: AtomView<'_>,
    coefficient: i32,
    use_uv_loop_basis: bool,
    loop_signature: &mut [i32],
    external_signature: &mut [i32],
) -> std::result::Result<(), String> {
    match view {
        AtomView::Add(add) => {
            for term in add.iter() {
                accumulate_momentum_signature(
                    graph,
                    term,
                    coefficient,
                    use_uv_loop_basis,
                    loop_signature,
                    external_signature,
                )?;
            }
        }
        AtomView::Mul(mul) => {
            let mut scalar = coefficient;
            let mut momentum_factor = None;
            for factor in mul.iter() {
                if let Ok(integer) = i64::try_from(factor) {
                    scalar = scalar
                        .checked_mul(i32::try_from(integer).map_err(|error| error.to_string())?)
                        .ok_or_else(|| {
                            "integer coefficient overflow in 4D denominator momentum".to_string()
                        })?;
                } else if momentum_factor.replace(factor).is_some() {
                    return Err(format!(
                        "expected a linear momentum expression, found product `{}`",
                        view.to_owned()
                    ));
                }
            }
            if let Some(momentum_factor) = momentum_factor {
                accumulate_momentum_signature(
                    graph,
                    momentum_factor,
                    scalar,
                    use_uv_loop_basis,
                    loop_signature,
                    external_signature,
                )?;
            } else if scalar != 0 {
                return Err(format!(
                    "expected a momentum factor, found scalar `{}`",
                    view.to_owned()
                ));
            }
        }
        AtomView::Fun(function) => {
            if function.get_symbol() != GS.emr_mom || function.get_nargs() != 1 {
                return Err(format!(
                    "expected an EMR momentum, found `{}`",
                    view.to_owned()
                ));
            }
            let edge = EdgeIndex(usize::try_from(function.get(0)).map_err(|_| {
                format!(
                    "EMR momentum has non-integer edge id `{}`",
                    function.get(0).to_owned()
                )
            })?);
            if use_uv_loop_basis
                && let Some(loop_index) = graph
                    .loop_momentum_basis
                    .loop_edges
                    .iter()
                    .position(|loop_edge| *loop_edge == edge)
            {
                loop_signature[loop_index] += coefficient;
                return Ok(());
            }
            let signature = graph
                .loop_momentum_basis
                .edge_signatures
                .get(edge)
                .ok_or_else(|| format!("edge {edge} is not present in the loop-momentum basis"))?;
            for (loop_index, sign) in signature.internal.iter_enumerated() {
                loop_signature[usize::from(loop_index)] += coefficient * sign_to_i32(*sign);
            }
            for (external_index, sign) in signature.external.iter_enumerated() {
                external_signature[usize::from(external_index)] += coefficient * sign_to_i32(*sign);
            }
        }
        AtomView::Num(_) if view.to_owned().is_zero() => {}
        AtomView::Num(_) => {
            return Err(format!(
                "expected a momentum expression, found scalar `{}`",
                view.to_owned()
            ));
        }
        _ => {
            return Err(format!(
                "expected a linear momentum expression, found `{}`",
                view.to_owned()
            ));
        }
    }
    Ok(())
}

fn exact_integer_rank(rows: &[Vec<i32>]) -> usize {
    rank_i64(
        &rows
            .iter()
            .map(|row| {
                row.iter()
                    .map(|coefficient| i64::from(*coefficient))
                    .collect()
            })
            .collect::<Vec<Vec<_>>>(),
    )
}

fn sign_to_i32(sign: SignOrZero) -> i32 {
    match sign {
        SignOrZero::Minus => -1,
        SignOrZero::Zero => 0,
        SignOrZero::Plus => 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        dot, graph::parse::IntoGraph, initialisation::test_initialise, momentum::sample::LoopIndex,
    };
    use color_eyre::Result;

    #[test]
    fn contracted_source_preserves_emr_ids_without_denominator_aliases() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
            d -> e [id=3 lmb_id=1]
            e -> f [id=4]
            f -> d [id=5]
            c -> d [id=6]
        })?;
        let contracted_edges = [EdgeIndex::from(0), EdgeIndex::from(1), EdgeIndex::from(2)];
        let source = GraphThreeDSource::new(&graph, &contracted_edges)?;
        assert_eq!(source.inner_loop_count, 1);
        assert_eq!(source.outer_loop_edges.len(), 1);

        let inner_carrier = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .find(|edge_id| contracted_edges.contains(edge_id))
            .expect("contracted loop should retain a parent-basis carrier");
        assert!(
            source
                .reconstructible_outer_loop_coordinates(inner_carrier)
                .is_none()
        );
        let parsed = source.to_three_d_parsed_graph()?;
        let edge_map = source
            .energy_edge_index_map(&parsed)
            .expect("GammaLoop sources provide an edge-index map");
        assert!(
            !edge_map
                .internal
                .values()
                .any(|edge| *edge == usize::from(inner_carrier))
        );
        let outer_carrier = source.outer_loop_edges[0];
        assert!(
            edge_map
                .internal
                .values()
                .any(|edge| *edge == usize::from(outer_carrier))
        );

        let contracted_tree_edge = EdgeIndex::from(6);
        let tree_source = GraphThreeDSource::new(&graph, &[contracted_tree_edge])?;
        assert_eq!(tree_source.inner_loop_count, 0);
        assert_eq!(
            tree_source.reconstructible_outer_loop_coordinates(contracted_tree_edge),
            tree_source.edge_loop_coordinates(contracted_tree_edge),
        );
        let parsed = tree_source.to_three_d_parsed_graph()?;
        let edge_map = tree_source
            .energy_edge_index_map(&parsed)
            .expect("GammaLoop sources provide an edge-index map");
        assert!(
            !edge_map
                .internal
                .values()
                .any(|edge| *edge == usize::from(contracted_tree_edge)),
            "a shrunken edge must not be aliased to a surviving denominator"
        );
        Ok(())
    }

    #[test]
    fn exact_source_preserves_dotted_same_edge_occurrences() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
            c -> d [id=3]
        })?;
        let edge = EdgeIndex::from(0);
        let momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(edge))
            .finish();
        let mass_squared = graph.underlying[edge].particle.mass_atom().pow(2);
        let first = FourDDenominator {
            source_edge: edge,
            momentum: momentum.clone(),
            mass_squared: mass_squared.clone(),
            full_expr: Atom::var(symbolica::symbol!("three_d_source_test::first")),
        };
        let second = FourDDenominator {
            source_edge: edge,
            momentum,
            mass_squared,
            full_expr: Atom::var(symbolica::symbol!("three_d_source_test::second")),
        };
        let denominators = [first.clone(), first, second];
        let source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        let parsed = source.to_three_d_parsed_graph()?;

        assert_eq!(parsed.internal_edges.len(), 3);
        assert_eq!(
            parsed
                .internal_edges
                .iter()
                .map(|edge| edge.edge_id)
                .collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
        let offset = graph.underlying.n_edges();
        assert_eq!(
            source.energy_edge_index_map(&parsed).unwrap().internal,
            BTreeMap::from([(0, offset), (1, offset + 1), (2, offset + 2)])
        );
        assert_eq!(
            source.physical_energy_edge_index_map().unwrap().internal,
            BTreeMap::from([
                (offset, usize::from(edge)),
                (offset + 1, usize::from(edge)),
                (offset + 2, usize::from(edge)),
            ])
        );
        Ok(())
    }

    #[test]
    fn exact_full_source_uses_canonical_production_carriers() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_production_basis {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0]
            b -> c [id=1]
            c -> a [id=2 lmb_id=0]
            c -> d [id=3]
            d -> e [id=4]
            e -> c [id=5 lmb_id=1]
        })?;
        let production = GraphThreeDSource::new(&graph, &[])?;
        let configured = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .collect::<Vec<_>>();
        assert_ne!(production.outer_loop_edges, configured);

        let denominators = (0..6)
            .map(|edge| {
                let edge = EdgeIndex(edge);
                FourDDenominator {
                    source_edge: edge,
                    momentum: FunctionBuilder::new(GS.emr_mom)
                        .add_arg(usize::from(edge))
                        .finish(),
                    mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
                    full_expr: Atom::one(),
                }
            })
            .collect::<Vec<_>>();
        let exact = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        assert_eq!(exact.outer_loop_edges, production.outer_loop_edges);

        let production_parsed = production.to_three_d_parsed_graph()?;
        let exact_parsed = exact.to_three_d_parsed_graph()?;
        assert_eq!(exact_parsed.internal_edges.len(), denominators.len());
        for (exact_edge, denominator) in exact_parsed.internal_edges.iter().zip(&denominators) {
            assert_eq!(
                exact_edge.signature,
                production_parsed.internal_edges[usize::from(denominator.source_edge)].signature
            );
        }

        let mut reversed_denominators = denominators.clone();
        reversed_denominators.reverse();
        let reversed = GraphThreeDSource::from_exact_denominators(&graph, &reversed_denominators)?;
        let reversed_parsed = reversed.to_three_d_parsed_graph()?;
        assert_eq!(reversed.outer_loop_edges, exact.outer_loop_edges);
        assert_eq!(reversed_parsed.loop_names, exact_parsed.loop_names);
        for (exact_edge, denominator) in reversed_parsed
            .internal_edges
            .iter()
            .zip(&reversed_denominators)
        {
            assert_eq!(
                exact_edge.signature,
                production_parsed.internal_edges[usize::from(denominator.source_edge)].signature
            );
        }
        Ok(())
    }

    #[test]
    fn exact_rank_deficient_source_keeps_the_complete_active_direction() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_rank_deficient {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0]
            b -> c [id=1]
            c -> a [id=2 lmb_id=0]
            c -> d [id=3]
            d -> e [id=4]
            e -> c [id=5 lmb_id=1]
        })?;
        let parent_carriers = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .collect::<Vec<_>>();
        let denominator = FourDDenominator {
            source_edge: EdgeIndex(0),
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(parent_carriers[0]))
                .finish()
                + FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(parent_carriers[1]))
                    .finish(),
            mass_squared: graph.underlying[EdgeIndex(0)].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        };
        let raw_signature = denominator.momentum_signature(&graph)?.loop_signature;
        let denominators = [denominator];
        let source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        assert_eq!(source.inner_loop_count, 1);
        assert_eq!(source.active_loop_count(), 1);

        let source_coordinates = (0..2)
            .map(|source_column| {
                raw_signature.iter().enumerate().fold(
                    Rational::from(0),
                    |sum, (parent_column, coefficient)| {
                        sum + Rational::from(*coefficient)
                            * source.parent_loop_coordinates[parent_column][source_column].clone()
                    },
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(source_coordinates[0], Rational::from(0));
        assert_eq!(
            source_coordinates[1]
                .numerator_ref()
                .to_i64()
                .unwrap()
                .abs(),
            1
        );
        assert_eq!(source_coordinates[1].denominator_ref().to_i64(), Some(1));
        assert_eq!(
            source.exact_signatures[0].loop_signature,
            vec![source_coordinates[1].numerator_ref().to_i64().unwrap() as i32]
        );
        Ok(())
    }

    #[test]
    fn exact_uv_carriers_from_one_source_edge_have_unique_loop_names() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_uv_names {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0]
            b -> c [id=1]
            c -> a [id=2 lmb_id=0]
            c -> d [id=3]
            d -> e [id=4]
            e -> c [id=5 lmb_id=1]
        })?;
        let denominators = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .map(|carrier| FourDDenominator {
                source_edge: EdgeIndex(0),
                momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(*carrier))
                    .finish(),
                mass_squared: Atom::var(GS.m_uv_expansion).pow(2),
                full_expr: Atom::one(),
            })
            .collect::<Vec<_>>();
        let source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        let parsed = source.to_three_d_parsed_graph()?;

        assert_eq!(parsed.loop_names.len(), 2);
        assert_eq!(
            parsed.loop_names.iter().collect::<BTreeSet<_>>().len(),
            parsed.loop_names.len()
        );
        Ok(())
    }

    #[test]
    fn exact_source_separates_loop_denominators_from_external_tree_factors() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph shifted_loop {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=4]
            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
            c -> d [id=3]
            d -> outgoing [id=5]
        })?;
        let loop_edge = graph.loop_momentum_basis.loop_edges[LoopIndex::from(0)];
        let tree_edge = EdgeIndex::from(3);
        let loop_denominator = FourDDenominator {
            source_edge: loop_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(loop_edge))
                .finish(),
            mass_squared: graph.underlying[loop_edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        };
        let tree_denominator = FourDDenominator {
            source_edge: tree_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(tree_edge))
                .finish(),
            mass_squared: graph.underlying[tree_edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        };

        assert!(loop_denominator.depends_on_loop(&graph)?);
        assert!(!tree_denominator.depends_on_loop(&graph)?);

        let external_edge = *graph
            .loop_momentum_basis
            .ext_edges
            .first()
            .expect("shifted loop graph has an external momentum carrier");
        let shifted = FourDDenominator {
            source_edge: loop_edge,
            momentum: &loop_denominator.momentum
                + FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(external_edge))
                    .finish(),
            mass_squared: Atom::var(GS.m_uv_expansion).pow(2),
            full_expr: Atom::one(),
        };
        let shifted_signature = shifted.momentum_signature(&graph)?;
        assert_eq!(shifted_signature.loop_signature, vec![1]);
        assert_eq!(
            shifted_signature.external_signature,
            (&graph.loop_momentum_basis.edge_signatures[external_edge].external)
                .into_iter()
                .map(sign_to_i32)
                .collect::<Vec<_>>()
        );
        Ok(())
    }

    #[test]
    fn pure_tree_exact_denominators_have_no_active_loop_source() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0]
            b -> c [id=1]
        })?;
        for edge in [EdgeIndex::from(0), EdgeIndex::from(1)] {
            let denominator = FourDDenominator {
                source_edge: edge,
                momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(edge))
                    .finish(),
                mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
                full_expr: Atom::one(),
            };
            assert!(!denominator.depends_on_loop(&graph)?);
        }
        Ok(())
    }
}
