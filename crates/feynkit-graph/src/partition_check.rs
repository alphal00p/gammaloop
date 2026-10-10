//! Validation of the stored cut and threshold-candidate partitions.
//!
//! [`PartitionCheck`] re-derives a stored partition on dense arrays indexed by
//! endpoint and by vertex. It accepts a partition only when the set-based
//! checks ([`FeynmanDiagram::partition_crossing`] and
//! [`FeynmanDiagram::summarize_cut_side`]) would accept it with the same
//! crossing edges, coupling orders and loop counts. Every other partition goes
//! through the set-based checks, which report its first error, if any.
use linnet::half_edge::{involution::Hedge, subgraph::ModifySubSet};

use super::*;

const LEFT: u8 = 1;
const RIGHT: u8 = 2;

impl FeynmanDiagram {
    /// Check that every stored cut and threshold candidate partitions the
    /// diagram, and that its oriented crossing edges (and for cuts, the side
    /// coupling orders and loop counts) are the ones the partition determines.
    pub(crate) fn validate_partitions(
        &self,
        paired_external_edges: &BTreeSet<EdgeId>,
    ) -> Result<(), DiagramError> {
        if self.cuts.is_empty() && self.topology_threshold_candidates.is_empty() {
            return Ok(());
        }
        let mut dense = PartitionCheck::new(self);
        // Only the set-based checks use the vertex crowns.
        let mut vertex_half_edges = None;
        for (index, cut) in self.cuts.iter().enumerate() {
            if !dense.accepts_cut(cut) {
                let vertex_half_edges =
                    vertex_half_edges.get_or_insert_with(|| self.vertex_half_edges());
                self.check_cut(index, cut, vertex_half_edges)?;
            }
        }
        if self.topology_threshold_candidates.is_empty() {
            return Ok(());
        }
        let ignored = self.initial_state_tree().0;
        let ignored_mask = dense.ignored_mask(&ignored);
        for (index, candidate) in self.topology_threshold_candidates.iter().enumerate() {
            if !dense.accepts_threshold_candidate(candidate, &ignored_mask, paired_external_edges) {
                let vertex_half_edges =
                    vertex_half_edges.get_or_insert_with(|| self.vertex_half_edges());
                self.check_threshold_candidate(
                    index,
                    candidate,
                    &ignored,
                    paired_external_edges,
                    vertex_half_edges,
                )?;
            }
        }
        Ok(())
    }

    /// The set-based check of one stored cut, which names its first error.
    fn check_cut(
        &self,
        index: usize,
        cut: &DiagramCut,
        vertex_half_edges: &[Vec<DiagramHalfEdge>],
    ) -> Result<(), DiagramError> {
        let invalid = |message: &str| DiagramError::InvalidCut {
            cut: index,
            message: message.to_owned(),
        };
        let oriented = cut.cut.iter().copied().collect::<BTreeSet<_>>();
        if oriented.len() != cut.cut.len() {
            return Err(invalid("half-edge sets contain duplicates"));
        }
        let normalized = self.cut_from_partitions(
            index,
            &cut.left.half_edges,
            &cut.right.half_edges,
            vertex_half_edges,
        )?;
        if oriented != normalized.cut.iter().copied().collect::<BTreeSet<_>>() {
            return Err(invalid(
                "oriented cut endpoints do not equal the crossing edges",
            ));
        }
        if normalized.left.coupling_orders != cut.left.coupling_orders
            || normalized.left.loop_count != cut.left.loop_count
            || normalized.right.coupling_orders != cut.right.coupling_orders
            || normalized.right.loop_count != cut.right.loop_count
        {
            return Err(invalid(
                "stored side coupling orders or loop counts do not match the finalized topology",
            ));
        }
        Ok(())
    }

    /// The set-based check of one stored threshold candidate, which names its
    /// first error.
    fn check_threshold_candidate(
        &self,
        index: usize,
        candidate: &DiagramThresholdCandidate,
        ignored: &SuBitGraph,
        paired_external_edges: &BTreeSet<EdgeId>,
        vertex_half_edges: &[Vec<DiagramHalfEdge>],
    ) -> Result<(), DiagramError> {
        let invalid = |message: String| DiagramError::InvalidThresholdCandidate {
            candidate: index,
            message,
        };
        let oriented = candidate.cut.iter().copied().collect::<BTreeSet<_>>();
        if oriented.len() != candidate.cut.len() {
            return Err(invalid("half-edge sets contain duplicates".to_owned()));
        }
        // Threshold candidates store no side summaries: only the crossing edges
        // are compared.
        let partition = self
            .partition_crossing(
                index,
                &candidate.left,
                &candidate.right,
                vertex_half_edges,
                Some(ignored),
            )
            .map_err(|error| match error {
                DiagramError::InvalidCut { message, .. } => invalid(message),
                error => error,
            })?;
        if candidate.non_initial_crossings(paired_external_edges) <= 1 {
            return Err(invalid(
                "threshold candidates must cross at least two non-initial-state edges".to_owned(),
            ));
        }
        if oriented != partition.crossing {
            return Err(invalid(
                "oriented cut endpoints do not equal the crossing edges".to_owned(),
            ));
        }
        Ok(())
    }
}

impl DiagramThresholdCandidate {
    /// Crossing edges that are not sewn initial states.
    fn non_initial_crossings(&self, paired_external_edges: &BTreeSet<EdgeId>) -> usize {
        self.cut
            .iter()
            .filter(|half_edge| !paired_external_edges.contains(&half_edge.edge))
            .count()
    }
}

/// Per-endpoint and per-vertex arrays of one diagram, reused for every stored
/// partition. Endpoint `(edge, endpoint)` has the dense index
/// `2 * edge + endpoint`, with the source endpoint first.
struct PartitionCheck<'a> {
    diagram: &'a FeynmanDiagram,
    /// Native half-edge and vertex of every existing endpoint, by dense index.
    existing: Vec<Option<(Hedge, usize)>>,
    /// Internal edges and whether they carry loop momentum, built at first use.
    internal: Option<Vec<(EdgeId, bool)>>,
    /// Side of every endpoint in the current partition: [`LEFT`], [`RIGHT`] or 0.
    marks: Vec<u8>,
    /// Sides touching every vertex in the current partition.
    seen: Vec<u8>,
    /// Coupling orders of the vertex rules met so far.
    rule_orders: Vec<(VertexRuleId, BTreeMap<String, usize>)>,
}

impl<'a> PartitionCheck<'a> {
    fn new(diagram: &'a FeynmanDiagram) -> Self {
        let n_edges = diagram.graph.n_edges();
        let mut existing = vec![None; 2 * n_edges];
        for half_edge in diagram.half_edges() {
            let hedge = diagram.half_edge_id(half_edge).expect("existing half-edge");
            let index = Self::index(half_edge, n_edges).expect("edge of the graph");
            let previous = existing[index].replace((hedge, diagram.graph.node_id(hedge).0));
            // The two half-edges of a paired edge have opposite flows.
            debug_assert!(previous.is_none(), "{half_edge:?} names two half-edges");
        }
        Self {
            diagram,
            marks: vec![0; existing.len()],
            seen: vec![0; diagram.graph.node_bound()],
            existing,
            internal: None,
            rule_orders: Vec::new(),
        }
    }

    /// Dense index of an endpoint reference, or `None` for an edge outside the
    /// graph. Deserialized partitions can name any edge id, and doubling one of
    /// at least 2^63 would wrap onto an existing endpoint.
    fn index(half_edge: DiagramHalfEdge, n_edges: usize) -> Option<usize> {
        (half_edge.edge.0 < n_edges).then(|| {
            2 * half_edge.edge.0 + usize::from(half_edge.endpoint == DiagramEndpoint::Target)
        })
    }

    /// Whether the set-based check accepts `cut`. `false` means it might not,
    /// for instance because the stored crossing list is not sorted.
    fn accepts_cut(&mut self, cut: &DiagramCut) -> bool {
        self.crossing(&cut.left.half_edges, &cut.right.half_edges, None)
            .is_some_and(|crossing| crossing == cut.cut)
            && self.side_matches(LEFT, &cut.left)
            && self.side_matches(RIGHT, &cut.right)
    }

    /// Whether the set-based check accepts `candidate`. `false` means it might
    /// not.
    fn accepts_threshold_candidate(
        &mut self,
        candidate: &DiagramThresholdCandidate,
        ignored_mask: &[bool],
        paired_external_edges: &BTreeSet<EdgeId>,
    ) -> bool {
        candidate.non_initial_crossings(paired_external_edges) > 1
            && self
                .crossing(&candidate.left, &candidate.right, Some(ignored_mask))
                .is_some_and(|crossing| crossing == candidate.cut)
    }

    /// Endpoints excluded from threshold sides, by dense index.
    fn ignored_mask(&self, ignored: &SuBitGraph) -> Vec<bool> {
        self.existing
            .iter()
            .map(|existing| existing.is_some_and(|(hedge, _)| ignored.includes(&hedge)))
            .collect()
    }

    /// Sorted oriented crossing endpoints of a valid partition. `None` when the
    /// sides repeat or overlap, miss or exceed the existing endpoints outside
    /// `ignored`, or split a vertex.
    fn crossing(
        &mut self,
        left: &[DiagramHalfEdge],
        right: &[DiagramHalfEdge],
        ignored: Option<&[bool]>,
    ) -> Option<Vec<DiagramHalfEdge>> {
        let n_edges = self.existing.len() / 2;
        self.marks.fill(0);
        for (side, half_edges) in [(LEFT, left), (RIGHT, right)] {
            for half_edge in half_edges {
                let mark = &mut self.marks[Self::index(*half_edge, n_edges)?];
                if *mark != 0 {
                    return None;
                }
                *mark = side;
            }
        }
        let is_ignored = |index: usize| ignored.is_some_and(|ignored| ignored[index]);
        for (index, mark) in self.marks.iter().enumerate() {
            if (*mark != 0) != (self.existing[index].is_some() && !is_ignored(index)) {
                return None;
            }
        }
        let diagram = self.diagram;
        let internal = self.internal.get_or_insert_with(|| {
            diagram
                .edges()
                .filter(|(_, endpoints, data)| {
                    data.external.is_none()
                        && endpoints.source.is_some()
                        && endpoints.target.is_some()
                })
                .map(|(edge, _, _)| {
                    let carries_loop_momentum = diagram.loop_momentum_basis.edge_signatures[&edge]
                        .loops
                        .iter()
                        .any(|sign| !sign.is_zero());
                    (edge, carries_loop_momentum)
                })
                .collect()
        });
        let marks = &self.marks;
        let left_of = |index: usize| marks[index] == LEFT;
        let right_of = |index: usize| marks[index] == RIGHT;
        let mut oriented = Vec::new();
        // The orientation rule of `FeynmanDiagram::partition_crossing`.
        for &(edge, crosses_ignored) in internal.iter() {
            let (source, target) = (2 * edge.0, 2 * edge.0 + 1);
            let source_to_target = left_of(source) && right_of(target)
                || crosses_ignored
                    && (left_of(source) && is_ignored(target)
                        || is_ignored(source) && right_of(target));
            let target_to_source = left_of(target) && right_of(source)
                || crosses_ignored
                    && (left_of(target) && is_ignored(source)
                        || is_ignored(target) && right_of(source));
            let endpoint = match (source_to_target, target_to_source) {
                (true, false) => DiagramEndpoint::Source,
                (false, true) => DiagramEndpoint::Target,
                _ => continue,
            };
            oriented.push(DiagramHalfEdge { edge, endpoint });
        }
        self.seen.fill(0);
        for (index, mark) in self.marks.iter().enumerate() {
            if *mark != 0 {
                let (_, vertex) = self.existing[index]?;
                self.seen[vertex] |= *mark;
                if self.seen[vertex] == LEFT | RIGHT {
                    return None;
                }
            }
        }
        oriented.sort();
        Some(oriented)
    }

    /// Whether one side of the last partition accepted by [`Self::crossing`]
    /// has the stored coupling orders and loop count.
    fn side_matches(&mut self, side: u8, stored: &DiagramCutSide) -> bool {
        let diagram = self.diagram;
        let rules = self
            .seen
            .iter()
            .enumerate()
            .filter(|(_, seen)| *seen & side != 0)
            .filter_map(|(vertex, _)| {
                diagram
                    .vertex(VertexId(vertex))
                    .and_then(|vertex| vertex.interaction)
            });
        let rule_orders = &mut self.rule_orders;
        for rule in rules.clone() {
            if !rule_orders.iter().any(|(id, _)| *id == rule) {
                let Ok(found) = diagram.model.vertex_rule_by_id(rule) else {
                    return false;
                };
                rule_orders.push((rule, found.coupling_orders(&diagram.model)));
            }
        }
        let mut orders = BTreeMap::<&str, usize>::new();
        for rule in rules {
            let (_, rule_orders) = rule_orders
                .iter()
                .find(|(id, _)| *id == rule)
                .expect("memoized rule");
            for (name, order) in rule_orders {
                *orders.entry(name.as_str()).or_insert(0) += order;
            }
        }
        if !stored
            .coupling_orders
            .iter()
            .map(|(name, order)| (name.as_str(), *order))
            .eq(orders)
        {
            return false;
        }
        let mut selected = diagram.graph.empty_subgraph::<SuBitGraph>();
        for (mark, existing) in self.marks.iter().zip(&self.existing) {
            if *mark == side
                && let Some((hedge, _)) = existing
            {
                selected.add(*hedge);
            }
        }
        diagram.graph.cyclotomatic_number(&selected) == stored.loop_count
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use serde_json::json;

    use super::*;
    use crate::tests::{cut_scalar_bubble, cut_scalar_line, scalar_model};

    /// Xorshift64, so that every run checks the same partitions.
    struct Rng(u64);

    impl Rng {
        fn below(&mut self, bound: usize) -> usize {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            (self.0 % bound.max(1) as u64) as usize
        }
    }

    /// The scalar model with a QED-ordered coupling on its cubic rule.
    fn coupled_scalar_model() -> Arc<Model> {
        let mut definition = serde_json::to_value(&*scalar_model()).unwrap();
        definition["orders"] = json!([{"name":"QED","expansion_order":99,"hierarchy":1}]);
        definition["couplings"] =
            json!([{"name":"GC_1","expression":"M","orders":[["QED",1]],"value":[1.0,0.0]}]);
        definition["vertex_rules"][0]["couplings"] = json!([["GC_1"]]);
        Arc::new(Model::from_json(&definition.to_string()).unwrap())
    }

    /// A three-loop scalar forward scattering: two physical cuts, four
    /// threshold candidates, and a side with a loop.
    fn three_loop_cross_section() -> FeynmanDiagram {
        FeynmanDiagram::from_dot(
            coupled_scalar_model(),
            r#"digraph three_loop {
                graph [final_state="phi,phi"];
                ext [style=invis];
                ext -> a [particle="phi", is_cut=0];
                b -> ext [particle="phi", is_cut=0];
                a -> c [particle="phi"];
                a -> d [particle="phi"];
                c -> d [particle="phi"];
                c -> b [particle="phi"];
                d -> b [particle="phi"];
            }"#,
        )
        .unwrap()
    }

    /// A bubble between two loop-independent bridges, which attach it to the
    /// initial state, and a second bubble on one attachment vertex. The
    /// initial-state tree holds the crowns of both attachment vertices, which
    /// threshold sides leave out, so the second bubble crosses into the tree.
    fn bridged_cross_section() -> FeynmanDiagram {
        let model = coupled_scalar_model();
        let rule = model.vertex_rule_id("V_1").unwrap();
        let phi = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "bridged-bubble");
        let untyped = |name: &str| DiagramVertex {
            name: name.into(),
            interaction: None,
            numerator: Atom::one(),
        };
        let a = builder.add_vertex(untyped("a"));
        let b = builder.add_vertex(untyped("b"));
        let c = builder.add_vertex(DiagramVertex::interaction("c", rule));
        let d = builder.add_vertex(DiagramVertex::interaction("d", rule));
        let e = builder.add_vertex(untyped("e"));
        let mut initial = DiagramEdge::new(phi, false);
        initial.external = Some(ExternalLeg {
            name: "p".into(),
            index: 0,
            state: ExternalState::Incoming,
            connection: 0,
        });
        builder.add_edge(a, b, initial).unwrap();
        for (source, target) in [(b, c), (c, d), (c, d), (d, a), (b, e), (e, b)] {
            builder
                .add_edge(source, target, DiagramEdge::new(phi, false))
                .unwrap();
        }
        let diagram = builder.build().unwrap();
        let crowns = diagram.vertex_half_edges();
        let sides = |left: &[VertexId]| {
            let (left, right): (Vec<_>, Vec<_>) =
                (0..crowns.len()).partition(|vertex| left.iter().any(|left| left.0 == *vertex));
            [left, right].map(|vertices| {
                vertices
                    .into_iter()
                    .flat_map(|vertex| crowns[vertex].clone())
                    .collect::<Vec<_>>()
            })
        };
        let [left, right] = sides(&[b, c, e]);
        let diagram = diagram
            .with_cut_partitions(vec![(left.clone(), right.clone())])
            .unwrap()
            .with_topology_threshold_partitions(vec![(left, right)])
            .unwrap();
        assert!(!diagram.initial_state_tree().1.is_empty());
        diagram
    }

    fn flip(endpoint: DiagramEndpoint) -> DiagramEndpoint {
        match endpoint {
            DiagramEndpoint::Source => DiagramEndpoint::Target,
            DiagramEndpoint::Target => DiagramEndpoint::Source,
        }
    }

    /// Remove, flip, repeat, swap or push one endpoint reference, or move it to
    /// an edge id whose double wraps onto its own dense index.
    fn corrupt(rng: &mut Rng, half_edges: &mut Vec<DiagramHalfEdge>, n_edges: usize) {
        let (position, other) = (rng.below(half_edges.len()), rng.below(half_edges.len()));
        match (rng.below(6), half_edges.is_empty()) {
            (0, false) => {
                half_edges.remove(position);
            }
            (1, false) => half_edges[position].endpoint = flip(half_edges[position].endpoint),
            (2, false) => half_edges.push(half_edges[position]),
            (3, false) => half_edges.swap(position, other),
            (4, false) => half_edges[position].edge.0 += 1 << (usize::BITS - 1),
            _ => half_edges.push(DiagramHalfEdge {
                edge: EdgeId(rng.below(n_edges + 2)),
                endpoint: [DiagramEndpoint::Source, DiagramEndpoint::Target][rng.below(2)],
            }),
        }
    }

    /// Move one half-edge, or every half-edge of one vertex, to the other side.
    fn regroup(
        rng: &mut Rng,
        crowns: &[Vec<DiagramHalfEdge>],
        left: &mut Vec<DiagramHalfEdge>,
        right: &mut Vec<DiagramHalfEdge>,
    ) {
        let (from, to) = if rng.below(2) == 0 {
            (left, right)
        } else {
            (right, left)
        };
        if rng.below(2) == 0 {
            if !from.is_empty() {
                to.push(from.remove(rng.below(from.len())));
            }
        } else {
            let crown = &crowns[rng.below(crowns.len())];
            to.extend(from.extract_if(.., |half_edge| crown.contains(half_edge)));
        }
    }

    #[test]
    fn dense_checks_accept_exactly_the_sorted_partitions_the_set_based_checks_accept() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        let strictly_sorted = |cut: &[DiagramHalfEdge]| cut.is_sorted_by(|a, b| a < b);
        for diagram in [
            cut_scalar_line(),
            cut_scalar_bubble(),
            three_loop_cross_section(),
            bridged_cross_section(),
        ] {
            let n_edges = diagram.graph.n_edges();
            let paired_external_edges = diagram
                .edges()
                .filter(|(_, endpoints, edge)| {
                    edge.external.is_some()
                        && endpoints.source.is_some()
                        && endpoints.target.is_some()
                })
                .map(|(edge, _, _)| edge)
                .collect::<BTreeSet<_>>();
            let crowns = diagram.vertex_half_edges();
            let ignored = diagram.initial_state_tree().0;
            let mut dense = PartitionCheck::new(&diagram);
            let ignored_mask = dense.ignored_mask(&ignored);
            // Every stored partition takes the dense path.
            assert!(diagram.cuts().iter().all(|cut| dense.accepts_cut(cut)));
            assert!(
                diagram
                    .topology_threshold_candidates()
                    .iter()
                    .all(|candidate| dense.accepts_threshold_candidate(
                        candidate,
                        &ignored_mask,
                        &paired_external_edges
                    ))
            );

            let mut outcomes = BTreeMap::<(&str, bool), usize>::new();
            let mut cuts = diagram.cuts().to_vec();
            for _ in 0..3000 {
                let mut cut = cuts[rng.below(cuts.len())].clone();
                let side = if rng.below(2) == 0 {
                    &mut cut.left
                } else {
                    &mut cut.right
                };
                match rng.below(6) {
                    0 => corrupt(&mut rng, &mut cut.cut, n_edges),
                    1 => corrupt(&mut rng, &mut side.half_edges, n_edges),
                    2 => side.loop_count = rng.below(3),
                    3 => *side.coupling_orders.entry("QED".to_owned()).or_default() += 1,
                    4 => side.coupling_orders.clear(),
                    _ => {
                        regroup(
                            &mut rng,
                            &crowns,
                            &mut cut.left.half_edges,
                            &mut cut.right.half_edges,
                        );
                        // Usually a valid cut with other crossings and sides.
                        if let Ok(recomputed) = diagram.cut_from_partitions(
                            0,
                            &cut.left.half_edges,
                            &cut.right.half_edges,
                            &crowns,
                        ) {
                            cut = recomputed;
                        }
                    }
                }
                let accepted = dense.accepts_cut(&cut);
                let valid = diagram.check_cut(0, &cut, &crowns).is_ok();
                assert_eq!(accepted, valid && strictly_sorted(&cut.cut), "{cut:?}");
                *outcomes.entry(("cut", accepted)).or_default() += 1;
                // Corrupt valid variants further, too.
                if valid {
                    cuts.push(cut);
                }
            }

            let mut candidates = diagram.topology_threshold_candidates().to_vec();
            for _ in 0..if candidates.is_empty() { 0 } else { 3000 } {
                let mut candidate = candidates[rng.below(candidates.len())].clone();
                match rng.below(4) {
                    0 => corrupt(&mut rng, &mut candidate.cut, n_edges),
                    1 => corrupt(&mut rng, &mut candidate.left, n_edges),
                    2 => corrupt(&mut rng, &mut candidate.right, n_edges),
                    _ => {
                        regroup(&mut rng, &crowns, &mut candidate.left, &mut candidate.right);
                        // A valid partition, possibly crossing too few edges.
                        if let Ok(partition) = diagram.partition_crossing(
                            0,
                            &candidate.left,
                            &candidate.right,
                            &crowns,
                            Some(&ignored),
                        ) {
                            candidate.cut = partition.crossing.into_iter().collect();
                        }
                    }
                }
                let accepted = dense.accepts_threshold_candidate(
                    &candidate,
                    &ignored_mask,
                    &paired_external_edges,
                );
                let valid = diagram
                    .check_threshold_candidate(
                        0,
                        &candidate,
                        &ignored,
                        &paired_external_edges,
                        &crowns,
                    )
                    .is_ok();
                assert_eq!(
                    accepted,
                    valid && strictly_sorted(&candidate.cut),
                    "{candidate:?}"
                );
                *outcomes.entry(("threshold", accepted)).or_default() += 1;
                if valid {
                    candidates.push(candidate);
                }
            }
            // Both outcomes occur for every kind of stored partition.
            let kinds = 1 + usize::from(!diagram.topology_threshold_candidates().is_empty());
            assert_eq!(
                outcomes.len(),
                2 * kinds,
                "{}: {outcomes:?}",
                diagram.name()
            );
        }
    }

    #[test]
    fn partitions_naming_edges_beyond_the_graph_are_rejected() {
        // Doubling an edge id of at least 2^63 wraps onto the endpoint of the
        // edge it was derived from.
        let alias = |half_edge: &mut DiagramHalfEdge| half_edge.edge.0 += 1 << (usize::BITS - 1);
        let diagram = cut_scalar_bubble();
        let mut cut = diagram.cuts()[0].clone();
        alias(&mut cut.left.half_edges[0]);
        assert!(matches!(
            diagram.clone().with_cuts(vec![cut]),
            Err(DiagramError::InvalidCut { cut: 0, message })
                if message.contains("disjoint and complementary")
        ));
        let mut candidate = diagram.topology_threshold_candidates()[0].clone();
        alias(&mut candidate.right[0]);
        assert!(matches!(
            diagram.clone().with_topology_threshold_candidates(vec![candidate]),
            Err(DiagramError::InvalidThresholdCandidate { candidate: 0, message })
                if message.contains("disjoint and complementary")
        ));

        // A deserialized diagram is validated without normalizing its cuts.
        let mut serialized = diagram.serde_view();
        alias(&mut serialized.cuts[0].left.half_edges[0]);
        assert!(matches!(
            FeynmanDiagram::from_json(
                diagram.model_arc(),
                &serde_json::to_string(&serialized).unwrap()
            ),
            Err(DiagramError::InvalidCut { cut: 0, message })
                if message.contains("disjoint and complementary")
        ));
    }
}
