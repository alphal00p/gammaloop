use std::collections::BTreeMap;

use itertools::Itertools;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::graph_signatures::MomentumSignature;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ParsedGraphInternalEdge {
    pub edge_id: usize,
    pub tail: usize,
    pub head: usize,
    pub label: String,
    pub mass_key: Option<String>,
    pub signature: MomentumSignature,
    pub had_pow: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ParsedGraphExternalEdge {
    pub edge_id: usize,
    pub source: Option<usize>,
    pub destination: Option<usize>,
    pub label: String,
    pub external_coefficients: Vec<i32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedGraph {
    pub internal_edges: Vec<ParsedGraphInternalEdge>,
    pub external_edges: Vec<ParsedGraphExternalEdge>,
    pub loop_names: Vec<String>,
    pub external_names: Vec<String>,
    pub node_name_to_internal: BTreeMap<String, usize>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct RepeatedGroup {
    pub key: RepeatedGroupKey,
    pub edge_ids: Vec<usize>,
    pub relative_signs: Vec<i32>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct RepeatedGroupKey {
    pub signature: MomentumSignature,
    pub mass_key: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GraphInfo {
    pub loop_names: Vec<String>,
    pub ext_names: Vec<String>,
    pub n_internal_edges: usize,
    pub n_external_edges: usize,
    pub repeated_groups: Vec<Vec<usize>>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EnergyEdgeIndexMap {
    pub internal: BTreeMap<usize, usize>,
    pub external: BTreeMap<usize, usize>,
    pub orientation_edge_count: usize,
}

impl EnergyEdgeIndexMap {
    pub fn identity(n_internal_edges: usize, n_external_edges: usize) -> Self {
        Self {
            internal: (0..n_internal_edges)
                .map(|edge_id| (edge_id, edge_id))
                .collect(),
            external: (0..n_external_edges)
                .map(|edge_id| (edge_id, edge_id))
                .collect(),
            orientation_edge_count: n_internal_edges,
        }
    }

    pub fn internal_to_local(&self) -> BTreeMap<usize, usize> {
        self.internal
            .iter()
            .map(|(local, target)| (*target, *local))
            .collect()
    }

    pub fn remap_bounds_to_local(
        &self,
        bounds: &[(usize, usize)],
    ) -> std::result::Result<Vec<(usize, usize)>, usize> {
        let reverse = self.internal_to_local();
        bounds
            .iter()
            .map(|(edge_id, degree)| {
                reverse
                    .get(edge_id)
                    .copied()
                    .map(|local_edge_id| (local_edge_id, *degree))
                    .ok_or(*edge_id)
            })
            .collect()
    }
}

#[derive(Debug, Error)]
pub enum GraphIoError {
    #[error("failed to extract three-dimensional representation graph input: {0}")]
    Source(String),
}

pub type Result<T> = std::result::Result<T, GraphIoError>;

pub trait ThreeDGraphSource {
    fn to_three_d_parsed_graph(&self) -> Result<ParsedGraph>;

    fn energy_edge_index_map(&self, parsed: &ParsedGraph) -> Option<EnergyEdgeIndexMap> {
        let _ = parsed;
        None
    }
}

impl ThreeDGraphSource for ParsedGraph {
    fn to_three_d_parsed_graph(&self) -> Result<ParsedGraph> {
        Ok(self.clone())
    }
}

pub fn repeated_groups(parsed: &ParsedGraph) -> Vec<RepeatedGroup> {
    let mut groups = BTreeMap::<RepeatedGroupKey, Vec<(usize, i32)>>::new();
    for edge in &parsed.internal_edges {
        let (signature, relative_sign) = edge.signature.canonical_up_to_sign();
        groups
            .entry(RepeatedGroupKey {
                signature,
                mass_key: edge.mass_key.clone(),
            })
            .or_default()
            .push((edge.edge_id, relative_sign));
    }

    groups
        .into_iter()
        .filter_map(|(key, mut members)| {
            if members.len() <= 1 {
                return None;
            }
            members.sort_by_key(|(edge_id, _)| *edge_id);
            Some(RepeatedGroup {
                key,
                edge_ids: members.iter().map(|(edge_id, _)| *edge_id).collect(),
                relative_signs: members.iter().map(|(_, sign)| *sign).collect(),
            })
        })
        .sorted_by_key(|group| group.edge_ids.clone())
        .collect()
}

pub fn graph_info(parsed: &ParsedGraph) -> GraphInfo {
    GraphInfo {
        loop_names: parsed.loop_names.clone(),
        ext_names: parsed.external_names.clone(),
        n_internal_edges: parsed.internal_edges.len(),
        n_external_edges: parsed.external_edges.len(),
        repeated_groups: repeated_groups(parsed)
            .into_iter()
            .map(|group| group.edge_ids)
            .collect(),
    }
}

#[cfg(test)]
pub(crate) mod test_graphs {
    use super::*;

    pub(crate) fn box_graph() -> ParsedGraph {
        ParsedGraph {
            internal_edges: vec![
                internal(0, 0, 1, "q0", [1], [0, 0, 0], "m1"),
                internal(1, 1, 2, "q1", [1], [1, 0, 0], "m2"),
                internal(2, 2, 3, "q2", [1], [1, 1, 0], "m3"),
                internal(3, 3, 0, "q3", [1], [1, 1, 1], "m4"),
            ],
            external_edges: box_external_edges(4, 3),
            loop_names: vec!["k1".to_string()],
            external_names: vec!["p1".to_string(), "p2".to_string(), "p3".to_string()],
            node_name_to_internal: node_map(4),
        }
    }

    pub(crate) fn box_pow3_graph() -> ParsedGraph {
        ParsedGraph {
            internal_edges: vec![
                internal(0, 0, 1, "q0", [1], [0, 0, 0], "m1"),
                internal(1, 1, 2, "q1", [1], [1, 0, 0], "m2"),
                internal(2, 2, 3, "q2", [1], [1, 1, 0], "m3"),
                internal(3, 3, 4, "q3", [1], [1, 1, 1], "m4"),
                internal(4, 4, 5, "q4", [1], [1, 1, 1], "m4"),
                internal(5, 5, 0, "q5", [1], [1, 1, 1], "m4"),
            ],
            external_edges: box_external_edges(6, 3),
            loop_names: vec!["k1".to_string()],
            external_names: vec!["p1".to_string(), "p2".to_string(), "p3".to_string()],
            node_name_to_internal: node_map(6),
        }
    }

    pub(crate) fn sunrise_pow4_graph() -> ParsedGraph {
        ParsedGraph {
            internal_edges: vec![
                internal(0, 0, 1, "q0", [-1, -1], [1], "m1"),
                internal(1, 0, 1, "q1", [1, 0], [0], "m2"),
                internal(2, 0, 2, "q2", [0, 1], [0], "m3"),
                internal(3, 2, 3, "q3", [0, 1], [0], "m3"),
                internal(4, 3, 4, "q4", [0, 1], [0], "m3"),
                internal(5, 4, 1, "q5", [0, 1], [0], "m3"),
            ],
            external_edges: vec![
                external(10_000_000, None, Some(0), "p1", [1]),
                external(10_000_001, Some(1), None, "-p1", [-1]),
            ],
            loop_names: vec!["k1".to_string(), "k2".to_string()],
            external_names: vec!["p1".to_string()],
            node_name_to_internal: node_map(5),
        }
    }

    pub(crate) fn imbalanced_graph() -> ParsedGraph {
        let mut parsed = box_graph();
        parsed.internal_edges[0].signature.loop_signature = vec![0];
        parsed
    }

    fn internal<const L: usize, const E: usize>(
        edge_id: usize,
        tail: usize,
        head: usize,
        label: &str,
        loop_signature: [i32; L],
        external_signature: [i32; E],
        mass_key: &str,
    ) -> ParsedGraphInternalEdge {
        ParsedGraphInternalEdge {
            edge_id,
            tail,
            head,
            label: label.to_string(),
            mass_key: Some(mass_key.to_string()),
            signature: MomentumSignature {
                loop_signature: loop_signature.to_vec(),
                external_signature: external_signature.to_vec(),
            },
            had_pow: false,
        }
    }

    fn external<const E: usize>(
        edge_id: usize,
        source: Option<usize>,
        destination: Option<usize>,
        label: &str,
        external_coefficients: [i32; E],
    ) -> ParsedGraphExternalEdge {
        ParsedGraphExternalEdge {
            edge_id,
            source,
            destination,
            label: label.to_string(),
            external_coefficients: external_coefficients.to_vec(),
        }
    }

    fn box_external_edges(n_nodes: usize, n_external: usize) -> Vec<ParsedGraphExternalEdge> {
        assert!(n_external <= n_nodes);
        let mut edges = Vec::new();
        let mut next_edge_id = 10_000_000;
        for external_id in 0..n_external {
            let mut coeffs = vec![0; n_external];
            coeffs[external_id] = 1;
            edges.push(ParsedGraphExternalEdge {
                edge_id: next_edge_id,
                source: None,
                destination: Some(external_id + 1),
                label: format!("p{}", external_id + 1),
                external_coefficients: coeffs.clone(),
            });
            next_edge_id += 1;

            coeffs[external_id] = -1;
            edges.push(ParsedGraphExternalEdge {
                edge_id: next_edge_id,
                source: Some(0),
                destination: None,
                label: format!("-p{}", external_id + 1),
                external_coefficients: coeffs,
            });
            next_edge_id += 1;
        }
        assert!(n_nodes > n_external);
        edges
    }

    fn node_map(n_nodes: usize) -> BTreeMap<String, usize> {
        (0..n_nodes)
            .map(|node| (format!("v{node}"), node))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::{test_graphs, *};

    #[test]
    fn repeated_groups_track_mass_and_signature_sign() {
        let parsed = test_graphs::box_pow3_graph();
        let groups = repeated_groups(&parsed);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].edge_ids, vec![3, 4, 5]);
        assert_eq!(groups[0].relative_signs, vec![-1, -1, -1]);
    }

    #[test]
    fn graph_info_tracks_box_shape() {
        let parsed = test_graphs::box_graph();
        let info = graph_info(&parsed);

        assert_eq!(info.loop_names, vec!["k1"]);
        assert_eq!(info.ext_names, vec!["p1", "p2", "p3"]);
        assert_eq!(info.n_internal_edges, 4);
        assert_eq!(info.n_external_edges, 6);
        assert!(info.repeated_groups.is_empty());
    }
}
