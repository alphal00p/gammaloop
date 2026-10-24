//! Edge-marked incidence graphs let native automorphism generators include
//! parallel physical propagators. No adapter graph-isomorphism search occurs.
use super::{charge, invalid};
use crate::VakintError;
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use symbolica::atom::Atom;
use symbolica::graph::{CanonicalForm, Empty, Graph};

pub(super) struct EdgeGraph {
    pub canonical: CanonicalForm<u8, Empty>,
    /// Canonical graph vertex for each physical propagator, in physical order.
    edges: Vec<usize>,
}
impl EdgeGraph {
    pub fn new(graph: &Graph<Empty, Atom>) -> Result<Self, VakintError> {
        let mut marked = Graph::new();
        for _ in graph.nodes() {
            marked.add_node(0u8);
        }
        let mut edges = Vec::new();
        for edge in graph.edges() {
            let vertex = marked.add_node(1u8);
            // Two incidences are retained even for a self-loop. Undirected
            // physical orientation remains a momentum-map proof obligation.
            marked
                .add_edge(vertex, edge.vertices.0, false, Empty)
                .map_err(invalid)?;
            marked
                .add_edge(vertex, edge.vertices.1, false, Empty)
                .map_err(invalid)?;
            edges.push(vertex);
        }
        let canonical = marked.canonize();
        let edges = edges
            .into_iter()
            .map(|vertex| canonical.vertex_map[vertex])
            .collect();
        Ok(Self { canonical, edges })
    }
}

pub(super) fn visit(
    source: &EdgeGraph,
    target: &EdgeGraph,
    limit: Option<usize>,
    work: &mut usize,
    mut visit: impl FnMut(&[usize]) -> Result<(), VakintError>,
) -> Result<(), VakintError> {
    if source.canonical.graph != target.canonical.graph {
        return Ok(());
    }
    let lookup = target
        .edges
        .iter()
        .enumerate()
        .map(|(edge, &vertex)| (vertex, edge))
        .collect::<BTreeMap<_, _>>();
    let base = source
        .edges
        .iter()
        .map(|vertex| {
            lookup
                .get(vertex)
                .copied()
                .ok_or_else(|| invalid("native canonical edge colors disagree"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut generators = BTreeSet::new();
    for cycles in &target.canonical.orbit_generators {
        let mut permutation = (0..target.canonical.graph.nodes().len()).collect::<Vec<_>>();
        for cycle in cycles {
            if cycle.is_empty() {
                continue;
            }
            for (&from, &to) in cycle
                .iter()
                .zip(cycle.iter().cycle().skip(1))
                .take(cycle.len())
            {
                if from >= permutation.len() || to >= permutation.len() {
                    return Err(invalid("native automorphism vertex is out of range"));
                }
                permutation[from] = to;
            }
        }
        generators.insert(
            target
                .edges
                .iter()
                .map(|&vertex| {
                    lookup
                        .get(&permutation[vertex])
                        .copied()
                        .ok_or_else(|| invalid("native automorphism changed an edge color"))
                })
                .collect::<Result<Vec<_>, _>>()?,
        );
    }
    let mut pending = VecDeque::from([base.clone()]);
    let mut seen = BTreeSet::from([base]);
    while let Some(permutation) = pending.pop_front() {
        charge(work, limit, "native edge automorphism witnesses")?;
        visit(&permutation)?;
        for generator in &generators {
            let next = permutation
                .iter()
                .map(|&edge| generator[edge])
                .collect::<Vec<_>>();
            if !seen.contains(&next) {
                // Explicit diagnostics also account for queued witnesses.
                if let Some(limit) = limit
                    && work.saturating_add(pending.len()).saturating_add(1) > limit
                {
                    return Err(VakintError::MomentumMatchResourceLimit {
                        resource: "native edge automorphism witnesses",
                        requested: work.saturating_add(pending.len()).saturating_add(1),
                        limit,
                    });
                }
                seen.insert(next.clone());
                pending.push_back(next);
            }
        }
    }
    Ok(())
}
