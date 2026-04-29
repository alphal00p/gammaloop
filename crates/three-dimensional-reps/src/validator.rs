use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::graph_io::{ParsedGraph, repeated_groups};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ParallelDuplicateViolation {
    pub edge_pair: (usize, usize),
    pub mass: Option<String>,
    pub edge_ids: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GraphValidation {
    pub ok: bool,
    pub n_internal_edges: usize,
    pub n_external_edges: usize,
    pub n_loops_from_labels: usize,
    pub n_external_symbols: usize,
    pub vertex_balance_violations: BTreeMap<usize, Vec<i32>>,
    pub vertex_external_balance_info: BTreeMap<usize, Vec<i32>>,
    pub parallel_duplicate_violations: Vec<ParallelDuplicateViolation>,
    pub pow_attribute_violations: Vec<usize>,
    pub repeated_groups: Vec<Vec<usize>>,
}

pub fn validate_parsed_graph(parsed: &ParsedGraph) -> GraphValidation {
    let loop_count = parsed.loop_names.len();
    let external_count = parsed.external_names.len();
    let mut balances = BTreeMap::<usize, Vec<i32>>::new();

    let mut pow_attribute_violations = Vec::new();
    for edge in &parsed.internal_edges {
        if edge.had_pow {
            pow_attribute_violations.push(edge.edge_id);
        }
        balances
            .entry(edge.tail)
            .or_insert_with(|| vec![0; loop_count + external_count]);
        balances
            .entry(edge.head)
            .or_insert_with(|| vec![0; loop_count + external_count]);
        let coeffs = edge
            .signature
            .loop_signature
            .iter()
            .chain(edge.signature.external_signature.iter())
            .copied()
            .collect::<Vec<_>>();
        for (index, coeff) in coeffs.iter().enumerate() {
            balances.get_mut(&edge.tail).unwrap()[index] -= coeff;
            balances.get_mut(&edge.head).unwrap()[index] += coeff;
        }
    }

    for edge in &parsed.external_edges {
        let mut coeffs = vec![0; loop_count];
        coeffs.extend(edge.external_coefficients.iter().copied());
        if let Some(source) = edge.source {
            balances
                .entry(source)
                .or_insert_with(|| vec![0; loop_count + external_count]);
            for (index, coeff) in coeffs.iter().enumerate() {
                balances.get_mut(&source).unwrap()[index] += coeff;
            }
        }
        if let Some(destination) = edge.destination {
            balances
                .entry(destination)
                .or_insert_with(|| vec![0; loop_count + external_count]);
            for (index, coeff) in coeffs.iter().enumerate() {
                balances.get_mut(&destination).unwrap()[index] += coeff;
            }
        }
    }

    let vertex_balance_violations: BTreeMap<usize, Vec<i32>> = balances
        .iter()
        .filter_map(|(vertex, balance)| {
            let loop_balance = balance[..loop_count].to_vec();
            if loop_balance.iter().any(|coeff| *coeff != 0) {
                Some((*vertex, loop_balance))
            } else {
                None
            }
        })
        .collect();

    let vertex_external_balance_info: BTreeMap<usize, Vec<i32>> = balances
        .iter()
        .filter_map(|(vertex, balance)| {
            let external_balance = balance[loop_count..].to_vec();
            if external_balance.iter().any(|coeff| *coeff != 0) {
                Some((*vertex, external_balance))
            } else {
                None
            }
        })
        .collect();

    let mut duplicate_map = BTreeMap::<_, Vec<usize>>::new();
    for edge in &parsed.internal_edges {
        duplicate_map
            .entry((
                (edge.tail, edge.head),
                edge.signature.clone(),
                edge.mass_key.clone(),
            ))
            .or_default()
            .push(edge.edge_id);
    }

    let parallel_duplicate_violations = duplicate_map
        .into_iter()
        .filter_map(|((edge_pair, _, mass), edge_ids)| {
            if edge_ids.len() > 1 {
                Some(ParallelDuplicateViolation {
                    edge_pair,
                    mass,
                    edge_ids,
                })
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    let repeated_groups = repeated_groups(parsed)
        .into_iter()
        .map(|group| group.edge_ids)
        .collect::<Vec<_>>();

    GraphValidation {
        ok: vertex_balance_violations.is_empty()
            && parallel_duplicate_violations.is_empty()
            && pow_attribute_violations.is_empty(),
        n_internal_edges: parsed.internal_edges.len(),
        n_external_edges: parsed.external_edges.len(),
        n_loops_from_labels: loop_count,
        n_external_symbols: external_count,
        vertex_balance_violations,
        vertex_external_balance_info,
        parallel_duplicate_violations,
        pow_attribute_violations,
        repeated_groups,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validates_local_parsed_graph_fixtures() {
        assert!(validate_parsed_graph(&crate::graph_io::test_graphs::box_graph()).ok);
        assert!(validate_parsed_graph(&crate::graph_io::test_graphs::box_pow3_graph()).ok);
        assert!(validate_parsed_graph(&crate::graph_io::test_graphs::sunrise_pow4_graph()).ok);
    }

    #[test]
    fn rejects_imbalanced_graph() {
        assert!(!validate_parsed_graph(&crate::graph_io::test_graphs::imbalanced_graph()).ok);
    }

    #[test]
    fn one_loop_external_shape_is_validated_locally() {
        let graph = crate::graph_io::test_graphs::box_graph();
        let validation = validate_parsed_graph(&graph);
        assert!(validation.ok);
        assert_eq!(validation.n_internal_edges, 4);
        assert_eq!(validation.n_loops_from_labels, 1);
        assert_eq!(validation.n_external_symbols, 3);
        assert_eq!(validation.n_external_edges, 6);
        assert!(validation.repeated_groups.is_empty());
    }
}
