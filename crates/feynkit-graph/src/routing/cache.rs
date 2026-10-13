//! Share canonical routing on an exact ordered amplitude skeleton.

use std::{
    collections::HashMap,
    sync::{Arc, RwLock},
};

use linnet::half_edge::involution::HedgePair;

use crate::{DiagramError, EdgeId, ExternalState, FeynmanDiagram, LoopMomentumBasis};

/// Request-owned canonical routing shared across particle decorations.
///
/// Keys retain exact endpoints, half-edge order and external roots. They omit
/// particle/rule payloads and physical orientation, which routing does not read.
/// No isomorphism search or change of momentum coordinates is performed here.
#[derive(Debug, Default)]
pub struct CanonicalRoutingCache {
    plans: RwLock<HashMap<Box<[usize]>, Arc<LoopMomentumBasis>>>,
}

impl CanonicalRoutingCache {
    /// One packed word sequence avoids hashing each enum field separately.
    /// The node/edge counts, edge tags and fixed-width records make its encoding
    /// unambiguous; the final length-prefixed sequence retains canonical DFS order.
    /// Paired externals, dummy attachments, cuts and vacuum routing stay on the
    /// general forest path. Only connected results are inserted into the cache.
    fn key(diagram: &FeynmanDiagram, ordered: &[EdgeId]) -> Option<Vec<usize>> {
        let graph = &diagram.graph;
        if !diagram.cuts.is_empty()
            || graph.node_bound() < 2
            || graph.n_live() != graph.node_bound()
        {
            return None;
        }
        let dangling = 2 * graph.n_edges() - graph.n_hedges();
        let mut key = Vec::with_capacity(3 + 5 * graph.n_edges() + 2 * dangling + ordered.len());
        key.extend([graph.node_bound(), graph.n_edges()]);
        let mut external_count = 0;
        for (pair, _, data) in graph.iter_edges() {
            if data.data.is_dummy {
                return None;
            }
            match (pair, &data.data.external) {
                (HedgePair::Paired { source, sink }, None) => {
                    // Keep both endpoint roles, including self-loop twins and
                    // nonadjacent half-edges. The same key fixes every crown.
                    key.extend([
                        0,
                        graph.node_id(source).0,
                        graph.node_id(sink).0,
                        source.0,
                        sink.0,
                    ]);
                }
                (HedgePair::Unpaired { hedge, flow }, Some(external)) => {
                    external_count += 1;
                    key.extend([
                        1,
                        graph.node_id(hedge).0,
                        hedge.0,
                        usize::from(flow.is_source()),
                        external.index,
                        external.connection,
                        usize::from(external.state == ExternalState::Outgoing),
                    ]);
                }
                _ => return None,
            }
        }
        key.push(ordered.len());
        key.extend(ordered.iter().map(|edge| edge.0));
        (external_count > 0).then_some(key)
    }

    pub(crate) fn route(
        &self,
        diagram: FeynmanDiagram,
        ordered: &[EdgeId],
    ) -> Result<FeynmanDiagram, DiagramError> {
        let Some(key) = Self::key(&diagram, ordered) else {
            return diagram.with_canonical_loop_momentum_basis_finalized(ordered);
        };
        if let Some(basis) = self
            .plans
            .read()
            .expect("routing cache poisoned")
            .get(key.as_slice())
        {
            return Ok(diagram.with_routing(Arc::clone(basis)));
        }
        let diagram = diagram.with_canonical_loop_momentum_basis_finalized(ordered)?;
        // A complete spanning tree certifies the connected amplitude domain.
        // Disconnected amplitudes retain their component-specific forest path.
        if diagram.loop_momentum_basis.tree_edges.len() + 1 == diagram.graph.node_bound() {
            self.plans
                .write()
                .expect("routing cache poisoned")
                .entry(key.into_boxed_slice())
                .or_insert_with(|| Arc::clone(&diagram.loop_momentum_basis));
        }
        Ok(diagram)
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use linnet::{
        half_edge::{involution::Flow, remap::HedgeRemap},
        permutation::Permutation,
    };
    use symbolica::atom::Atom;

    use super::*;
    use crate::{DiagramEdge, DiagramVertex, ExternalLeg, FeynmanDiagramBuilder, VertexSlot};

    fn amplitude(
        reversals: usize,
        root: usize,
        isolated: bool,
        cache: Option<Arc<CanonicalRoutingCache>>,
    ) -> FeynmanDiagram {
        let model = crate::tests::scalar_model();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagramBuilder::new(model, None).canonical_loop_momentum_basis();
        if let Some(cache) = cache {
            builder = builder.canonical_routing_cache(cache);
        }
        let nodes = (0..3 + usize::from(isolated))
            .map(|i| {
                builder.add_vertex(DiagramVertex {
                    name: format!("v{i}"),
                    interaction: None,
                    numerator: Atom::one(),
                })
            })
            .collect::<Vec<_>>();
        for (index, state, vertex, flow) in [
            (0, ExternalState::Incoming, root, Flow::Sink),
            (1, ExternalState::Outgoing, 2, Flow::Source),
        ] {
            builder
                .add_generation_external(
                    nodes[vertex],
                    flow,
                    DiagramEdge::new(particle, false),
                    VertexSlot(0),
                    ExternalLeg {
                        name: format!("ext{index}"),
                        index,
                        state,
                        connection: index,
                    },
                )
                .unwrap();
        }
        for (i, (a, b)) in [(0, 1), (0, 1), (1, 2), (2, 0), (1, 1)]
            .into_iter()
            .enumerate()
        {
            let (a, b) = if reversals & (1 << i) == 0 {
                (a, b)
            } else {
                (b, a)
            };
            builder
                .add_edge(nodes[a], nodes[b], DiagramEdge::new(particle, false))
                .unwrap();
        }
        builder.build().unwrap()
    }

    #[test]
    fn exact_skeleton_cache_preserves_all_internal_directions_and_external_roots() {
        let cache = Arc::new(CanonicalRoutingCache::default());
        for root in [0, 1] {
            for reversed in 0..32 {
                let expected = amplitude(reversed, root, false, None);
                let actual = amplitude(reversed, root, false, Some(Arc::clone(&cache)));
                assert_eq!(actual.loop_momentum_basis, expected.loop_momentum_basis);
                actual.loop_momentum_basis.validate(&actual).unwrap();
                let repeated = amplitude(reversed, root, false, Some(Arc::clone(&cache)));
                assert!(Arc::ptr_eq(
                    &actual.loop_momentum_basis,
                    &repeated.loop_momentum_basis
                ));
            }
        }
        // Reversing a self-loop's equal node endpoints does not alter its
        // builder half-edge order; the other four edges give sixteen frames.
        assert_eq!(cache.plans.read().unwrap().len(), 32);
    }

    #[test]
    fn skeleton_plan_preserves_self_loop_half_edge_coordinates() {
        let cache = CanonicalRoutingCache::default();
        let diagram = amplitude(0, 0, false, None);
        let ordered = diagram.canonical_internal_edge_order().unwrap();
        cache.route(diagram.clone(), &ordered).unwrap();
        let mut reversed = diagram;
        let (source, sink) = reversed
            .graph
            .iter_edges()
            .find_map(|(pair, _, _)| {
                if let HedgePair::Paired { source, sink } = pair
                    && reversed.graph.node_id(source) == reversed.graph.node_id(sink)
                {
                    Some((source, sink))
                } else {
                    None
                }
            })
            .unwrap();
        let mut map = (0..reversed.graph.n_hedges()).collect::<Vec<_>>();
        map.swap(source.0, sink.0);
        reversed
            .graph
            .permute_hedges(&HedgeRemap::from_permutation(&Permutation::from_map(map)))
            .unwrap();
        let expected = reversed
            .clone()
            .with_canonical_loop_momentum_basis_finalized(&ordered)
            .unwrap();
        let actual = cache.route(reversed, &ordered).unwrap();
        assert_eq!(actual.loop_momentum_basis, expected.loop_momentum_basis);
        actual.loop_momentum_basis.validate(&actual).unwrap();
        assert_eq!(cache.plans.read().unwrap().len(), 2);
    }

    #[test]
    fn disconnected_and_cut_diagrams_keep_the_general_routing_path() {
        let cache = Arc::new(CanonicalRoutingCache::default());
        let expected = amplitude(5, 0, true, None);
        let actual = amplitude(5, 0, true, Some(Arc::clone(&cache)));
        assert_eq!(actual.loop_momentum_basis, expected.loop_momentum_basis);
        assert!(cache.plans.read().unwrap().is_empty());

        let cut = crate::tests::cut_scalar_bubble();
        let ordered = cut.canonical_internal_edge_order().unwrap();
        assert!(CanonicalRoutingCache::key(&cut, &ordered).is_none());
        let expected = cut
            .clone()
            .with_canonical_loop_momentum_basis_finalized(&ordered)
            .unwrap();
        let actual = cache.route(cut, &ordered).unwrap();
        assert_eq!(actual.loop_momentum_basis, expected.loop_momentum_basis);
        assert!(cache.plans.read().unwrap().is_empty());
    }
    #[test]
    fn interleaved_internal_pairs_keep_their_distinct_incidence_order() {
        let cache = CanonicalRoutingCache::default();
        let mut diagram = amplitude(0, 0, false, None);
        let ordered = diagram.canonical_internal_edge_order().unwrap();
        let original_key = CanonicalRoutingCache::key(&diagram, &ordered).unwrap();
        cache.route(diagram.clone(), &ordered).unwrap();
        let paired = diagram
            .graph
            .iter_edges()
            .filter_map(|(pair, _, _)| {
                if let HedgePair::Paired { source, sink } = pair {
                    Some((source, sink))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        let mut map = (0..diagram.graph.n_hedges()).collect::<Vec<_>>();
        map.swap(paired[0].1.0, paired[1].1.0);
        diagram
            .graph
            .permute_hedges(&HedgeRemap::from_permutation(&Permutation::from_map(map)))
            .unwrap();
        assert_ne!(
            CanonicalRoutingCache::key(&diagram, &ordered).unwrap(),
            original_key
        );
        let expected = diagram
            .clone()
            .with_canonical_loop_momentum_basis_finalized(&ordered)
            .unwrap();
        let actual = cache.route(diagram, &ordered).unwrap();
        assert_eq!(actual.loop_momentum_basis, expected.loop_momentum_basis);
        actual.loop_momentum_basis.validate(&actual).unwrap();
        assert_eq!(cache.plans.read().unwrap().len(), 2);
    }

    #[test]
    fn shared_basis_outlives_cache_and_payload_edits_cannot_mutate_other_diagrams() {
        let cache = CanonicalRoutingCache::default();
        let first = amplitude(0, 0, false, None);
        let ordered = first.canonical_internal_edge_order().unwrap();
        let first = cache.route(first, &ordered).unwrap();
        let mut decorated = amplitude(0, 0, false, None);
        decorated.graph[linnet::half_edge::NodeIndex(0)].numerator = Atom::num(37);
        let second = cache.route(decorated, &ordered).unwrap();
        assert!(Arc::ptr_eq(
            &first.loop_momentum_basis,
            &second.loop_momentum_basis
        ));
        drop(cache);
        let original = first.loop_momentum_basis().clone();
        let mut loops = original.loop_edges.clone();
        loops.reverse();
        let changed = second.with_loop_momentum_edges(&loops).unwrap();
        assert_eq!(first.loop_momentum_basis(), &original);
        assert_ne!(
            first.loop_momentum_basis().loop_edges,
            changed.loop_momentum_basis().loop_edges
        );
        first.loop_momentum_basis().validate(&first).unwrap();
        changed.loop_momentum_basis().validate(&changed).unwrap();
    }
}
