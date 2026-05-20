use std::collections::BTreeMap;

use ahash::AHashSet;
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, Flow, HedgePair};
use symbolica::atom::AtomCore;
use symbolica::domains::rational::Rational;
use three_dimensional_reps::{
    EnergyEdgeIndexMap, MomentumSignature, ParsedGraph, ThreeDGraphSource,
    graph_io::{
        GraphIoError, ParsedGraphExternalEdge, ParsedGraphInitialStateCutEdge,
        ParsedGraphInternalEdge, initial_state_cut_external_alias,
    },
};

use crate::{
    graph::Graph,
    momentum::{SignOrZero, sample::LoopIndex},
};

impl ThreeDGraphSource for Graph {
    fn to_three_d_parsed_graph(&self) -> three_dimensional_reps::graph_io::Result<ParsedGraph> {
        GraphThreeDSource::new(self, &[]).to_three_d_parsed_graph()
    }

    fn energy_edge_index_map(&self, parsed: &ParsedGraph) -> Option<EnergyEdgeIndexMap> {
        GraphThreeDSource::new(self, &[]).energy_edge_index_map(parsed)
    }
}

pub(crate) struct GraphThreeDSource<'a> {
    graph: &'a Graph,
    contract_edges: AHashSet<EdgeIndex>,
    initial_state_cut_edges: AHashSet<EdgeIndex>,
}

impl<'a> GraphThreeDSource<'a> {
    pub(crate) fn new(graph: &'a Graph, contract_edges: &[EdgeIndex]) -> Self {
        Self {
            graph,
            contract_edges: contract_edges.iter().copied().collect(),
            initial_state_cut_edges: graph
                .iter_edges_of(&graph.initial_state_cut)
                .map(|(_, edge_id, _)| edge_id)
                .collect(),
        }
    }

    fn contracts_edge(&self, edge_id: EdgeIndex) -> bool {
        self.contract_edges.contains(&edge_id) && !self.initial_state_cut_edges.contains(&edge_id)
    }

    fn active_loop_signature_rows(&self) -> Vec<Vec<i32>> {
        self.graph
            .underlying
            .iter_edges()
            .sorted_by_key(|(_, edge_index, _)| *edge_index)
            .filter_map(|(pair, edge_index, edge_data)| {
                (pair.is_paired()
                    && !edge_data.data.is_dummy
                    && !self.initial_state_cut_edges.contains(&edge_index)
                    && !self.contracts_edge(edge_index))
                .then(|| {
                    self.graph.loop_momentum_basis.edge_signatures[edge_index]
                        .internal
                        .iter()
                        .map(|sign| sign_to_i32(*sign))
                        .collect::<Vec<_>>()
                })
            })
            .collect_vec()
    }
}

impl ThreeDGraphSource for GraphThreeDSource<'_> {
    fn to_three_d_parsed_graph(&self) -> three_dimensional_reps::graph_io::Result<ParsedGraph> {
        // Keep contraction virtual here: linnet graph contraction mutates/deletes
        // edges, while the 3D expression must still map back to original EdgeIndex
        // values used by GammaLoop energy caches and residual denominators.
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

        let active_loop_signature_rows = self.active_loop_signature_rows();
        let active_loop_columns = independent_loop_columns(&active_loop_signature_rows);

        let loop_names = self
            .graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .enumerate()
            .filter(|(loop_index, _)| active_loop_columns.contains(loop_index))
            .map(|(_, edge_id)| self.graph.underlying[*edge_id].name.value.clone())
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
                loop_signature: active_loop_columns
                    .iter()
                    .map(|loop_index| sign_to_i32(signature.internal[LoopIndex::from(*loop_index)]))
                    .collect(),
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

fn independent_loop_columns(rows: &[Vec<i32>]) -> Vec<usize> {
    let Some(column_count) = rows.iter().map(Vec::len).max() else {
        return Vec::new();
    };
    let mut selected = Vec::new();
    let mut rank = 0;
    for column in 0..column_count {
        let mut candidate = selected.clone();
        candidate.push(column);
        let candidate_rank = integer_matrix_column_rank(rows, &candidate);
        if candidate_rank > rank {
            selected.push(column);
            rank = candidate_rank;
        }
    }
    selected
}

fn integer_matrix_column_rank(rows: &[Vec<i32>], columns: &[usize]) -> usize {
    let mut matrix = rows
        .iter()
        .map(|row| {
            columns
                .iter()
                .map(|column| Rational::from(row.get(*column).copied().unwrap_or_default()))
                .collect::<Vec<_>>()
        })
        .filter(|row| row.iter().any(|entry| *entry != 0))
        .collect::<Vec<_>>();

    let mut rank = 0;
    let mut pivot_column = 0;
    while rank < matrix.len() && pivot_column < columns.len() {
        let Some(pivot_row) = (rank..matrix.len()).find(|row| matrix[*row][pivot_column] != 0)
        else {
            pivot_column += 1;
            continue;
        };
        matrix.swap(rank, pivot_row);
        let pivot = matrix[rank][pivot_column].clone();
        for entry in matrix[rank].iter_mut().skip(pivot_column) {
            *entry /= pivot.clone();
        }
        let pivot_row = matrix[rank].clone();
        for (row, row_values) in matrix.iter_mut().enumerate() {
            if row == rank {
                continue;
            }
            let factor = row_values[pivot_column].clone();
            if factor == 0 {
                continue;
            }
            for (entry, pivot_entry) in row_values
                .iter_mut()
                .zip(pivot_row.iter())
                .skip(pivot_column)
            {
                *entry -= factor.clone() * pivot_entry.clone();
            }
        }
        rank += 1;
        pivot_column += 1;
    }
    rank
}

fn sign_to_i32(sign: SignOrZero) -> i32 {
    match sign {
        SignOrZero::Minus => -1,
        SignOrZero::Zero => 0,
        SignOrZero::Plus => 1,
    }
}
