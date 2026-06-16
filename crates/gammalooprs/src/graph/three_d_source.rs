use std::collections::BTreeMap;

use ahash::AHashSet;
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeIndex, Flow, HedgePair},
    subgraph::{InternalSubGraph, ModifySubSet, SuBitGraph},
};
use symbolica::atom::AtomCore;
use symbolica::domains::rational::Rational;
use three_dimensional_reps::{
    EnergyEdgeIndexMap, MomentumSignature, ParsedGraph, ThreeDGraphSource,
    graph_io::{
        GraphIoError, ParsedGraphExternalEdge, ParsedGraphInitialStateCutEdge,
        ParsedGraphInternalEdge, initial_state_cut_external_alias,
    },
    utils::{rank_i64, solve_rational_system},
};

use crate::{
    graph::{Graph, LMBext},
    momentum::SignOrZero,
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
        Ok(Self {
            graph,
            contract_edges,
            initial_state_cut_edges,
            inner_loop_count,
            outer_loop_edges,
            edge_loop_coordinates,
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
}

impl ThreeDGraphSource for GraphThreeDSource<'_> {
    fn to_three_d_parsed_graph(&self) -> three_dimensional_reps::graph_io::Result<ParsedGraph> {
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
    use crate::{dot, graph::parse::IntoGraph, initialisation::test_initialise};
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
}
