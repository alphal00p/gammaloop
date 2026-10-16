//! Physics invariants retained from GammaLoop's completed soft foundation.

use super::*;

#[test]
fn paired_incoming_three_cycles_preserve_the_canonical_connection_frame() {
    let model = Model::phi4();
    let particle = model.particle_id("phi").unwrap();
    let interaction = model.vertex_rule_id(&model.vertex_rules()[0].name).unwrap();
    // An asymmetric attachment graph exercises the canonical coordinate map
    // independently of Feynman-rule evaluation. A three-cycle is not its own inverse.
    let mut graph = Graph::new();
    for _ in 0..6 {
        graph.add_node(ColoredNode::Interaction(interaction));
    }
    for index in 0..6 {
        let external = graph.add_node(ColoredNode::External(ExternalNode {
            index,
            state: if index < 3 {
                ExternalState::Incoming
            } else {
                ExternalState::Outgoing
            },
            particle,
        }));
        let (source, target) = if index < 3 {
            (external, index)
        } else {
            (index, external)
        };
        graph
            .add_edge(
                source,
                target,
                false,
                EdgeColor {
                    particle,
                    direction: None,
                },
            )
            .unwrap();
    }
    for (source, target) in [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (0, 2),
        (0, 3),
        (1, 4),
        (3, 5),
    ] {
        graph
            .add_edge(
                source,
                target,
                false,
                EdgeColor {
                    particle,
                    direction: None,
                },
            )
            .unwrap();
    }
    let relabel = |paired: bool| {
        let mut output = Graph::new();
        for node in graph.nodes() {
            let mut data = node.data.clone();
            if let ColoredNode::External(external) = &mut data
                && (paired || external.index < 3)
            {
                external.index = (external.index / 3) * 3 + (external.index + 1) % 3;
            }
            output.add_node(data);
        }
        for edge in graph.edges() {
            output
                .add_edge(edge.vertices.0, edge.vertices.1, edge.directed, edge.data)
                .unwrap();
        }
        output
    };
    let process = Process::new(["phi", "phi", "phi"], ["phi", "phi"]);
    for left_right in [false, true] {
        let options = GenerationOptions::default()
            .symmetrize_initial(true)
            .symmetrize_left_right(left_right);
        let resolved =
            ResolvedProcess::new(&model, &process, GenerationType::CrossSection, &options).unwrap();
        let canonical = resolved
            .canonicalize_numerator_graph(&graph, &model, &NumeratorGrouping::None)
            .unwrap();
        assert_eq!(
            canonical,
            resolved
                .canonicalize_numerator_graph(&relabel(true), &model, &NumeratorGrouping::None)
                .unwrap()
        );
        assert_eq!(
            canonical,
            resolved
                .canonicalize_numerator_graph(&canonical, &model, &NumeratorGrouping::None)
                .unwrap()
        );
        assert_ne!(
            canonical,
            resolved
                .canonicalize_numerator_graph(&relabel(false), &model, &NumeratorGrouping::None)
                .unwrap(),
            "one-sided permutations must not exchange physical incoming pairs"
        );
    }
}

#[test]
fn identical_incoming_scalar_born_weights_only_divide_by_final_state_symmetry() {
    use feynkit_graph::expressions::evaluate_overall_factor;

    let model = Arc::new(Model::phi4());
    let process = Process::new(["phi", "phi"], ["phi", "phi"]);
    for initial in [false, true] {
        for final_state in [false, true] {
            for left_right in [false, true] {
                let options = GenerationOptions::default()
                    .with_loop_count(1, 1)
                    .unwrap()
                    .max_vertices(2)
                    .symmetrize_initial(initial)
                    .symmetrize_final(final_state)
                    .symmetrize_left_right(left_right)
                    .threads(1);
                let result = process
                    .generate_cross_section(Arc::clone(&model), &options)
                    .unwrap();
                assert_eq!(result.diagrams.len(), 1);
                assert_eq!(result.diagrams[0].cuts().len(), 1);
                assert_eq!(
                    evaluate_overall_factor(result.diagrams[0].overall_factor().as_view()),
                    Atom::num((1, 2)),
                    "initial={initial}, final={final_state}, left_right={left_right}: a labelled two-particle flux has no extra initial-state factorial"
                );
            }
        }
    }
}

#[test]
fn directed_self_loop_symmetry_factors_match_wick_contractions() {
    let neutral = HalfEdge::undirected("a");
    let charged = HalfEdge::incoming("q");
    let other = HalfEdge::incoming("r");
    // Production consumes these graph automorphism denominators directly.
    // Directed loop ends are distinct particle/antiparticle fields; neutral
    // flips and permutations of identical propagators still count.
    for (signature, vertices, expected) in [
        (vec![charged, charged.flip()], 1, vec![(1, 1)]),
        (
            vec![charged, charged.flip(), other, other.flip()],
            1,
            vec![(2, 1)],
        ),
        (
            vec![charged, charged.flip(), charged, charged.flip()],
            1,
            vec![(2, 2)],
        ),
        (vec![neutral; 4], 1, vec![(0, 8)]),
        (
            vec![neutral, neutral, charged, charged.flip()],
            1,
            vec![(1, 2)],
        ),
        (
            vec![neutral, charged, charged.flip()],
            2,
            vec![(0, 2), (2, 2)],
        ),
    ] {
        let graphs = Graph::<usize, &str>::generate(
            &[],
            std::slice::from_ref(&signature),
            GenerationSettings::new()
                .max_vertices(vertices)
                .max_loops(2)
                .allow_self_loops(true)
                .allow_zero_flow_edges(true),
        )
        .unwrap();
        let mut observed = graphs
            .iter()
            .map(|(graph, symmetry)| {
                (
                    graph
                        .edges()
                        .iter()
                        .filter(|edge| edge.directed && edge.vertices.0 == edge.vertices.1)
                        .count(),
                    u64::try_from(symmetry.clone()).unwrap(),
                )
            })
            .collect::<Vec<_>>();
        observed.sort();
        assert_eq!(observed, expected, "vertex signature {signature:?}");
    }
}

struct GhostChart {
    model: Model,
}

impl GhostChart {
    fn new() -> Self {
        Self {
            model: Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json"))
                .unwrap(),
        }
    }

    fn process(&self, cross_section: bool, cp: bool) -> ResolvedProcess {
        ResolvedProcess::new(
            &self.model,
            &Process::new(
                [22_i64],
                if cross_section {
                    vec![1_i64, -1]
                } else {
                    vec![22]
                },
            ),
            if cross_section {
                GenerationType::CrossSection
            } else {
                GenerationType::Amplitude
            },
            &GenerationOptions::default().symmetrize_left_right(cp),
        )
        .unwrap()
    }

    fn graph(
        &self,
        triangle: bool,
        reverse: bool,
        reverse_storage: bool,
    ) -> Graph<ColoredNode, EdgeColor> {
        let mut graph = Graph::new();
        for node in 0..8 {
            let rule = match (triangle, node) {
                (true, 0..=2) | (false, 0..=3) => "V_35",
                (true, 3..=5) | (false, 4..=5) => "V_74",
                _ => "V_71",
            };
            graph.add_node(ColoredNode::Interaction(
                self.model.vertex_rule_id(rule).unwrap(),
            ));
        }
        for (index, state) in [(0, ExternalState::Incoming), (1, ExternalState::Outgoing)] {
            graph.add_node(ColoredNode::External(ExternalNode {
                index,
                state,
                particle: self.model.particle_id_by_pdg(22).unwrap(),
            }));
        }
        let edges = if triangle {
            [
                (0, 1, 9000005),
                (2, 0, 9000005),
                (0, 3, 21),
                (1, 2, 9000005),
                (1, 5, 21),
                (2, 4, 21),
                (3, 4, 1),
                (5, 3, 1),
                (4, 6, 1),
                (7, 5, 1),
                (6, 7, 1),
            ]
        } else {
            [
                (0, 2, 9000005),
                (3, 0, 9000005),
                (0, 4, 21),
                (2, 1, 9000005),
                (1, 3, 9000005),
                (1, 5, 21),
                (2, 3, 21),
                (4, 5, 1),
                (6, 4, 1),
                (5, 7, 1),
                (7, 6, 1),
            ]
        };
        for (a, b, pdg) in edges.into_iter().chain([(8, 6, 22), (7, 9, 22)]) {
            let ghost = pdg == 9000005;
            let (a, b) = if reverse
                && (if triangle {
                    ghost
                } else {
                    pdg != 21 && pdg != 22
                }) {
                (b, a)
            } else {
                (a, b)
            };
            let (a, b, pdg) = if reverse_storage && ghost {
                (b, a, -pdg)
            } else {
                (a, b, pdg)
            };
            graph
                .add_edge(
                    a,
                    b,
                    true,
                    EdgeColor {
                        particle: self.model.particle_id_by_pdg(pdg).unwrap(),
                        direction: None,
                    },
                )
                .unwrap();
        }
        graph
    }

    fn chart(
        &self,
        graph: &Graph<ColoredNode, EdgeColor>,
        reverse_quarks: bool,
    ) -> Vec<((usize, usize), bool, EdgeColor)> {
        let mut edges = graph
            .edges()
            .iter()
            .map(|edge| {
                let particle = self.model.particle_by_id(edge.data.particle).unwrap();
                let endpoints = if reverse_quarks && particle.is_fermion() && !particle.is_ghost() {
                    (edge.vertices.1, edge.vertices.0)
                } else {
                    edge.vertices
                };
                (endpoints, edge.directed, edge.data)
            })
            .collect::<Vec<_>>();
        edges.sort();
        edges
    }
}

#[test]
fn ghost_box_reflection_preserves_the_derivative_momentum_chart() {
    let fixture = GhostChart::new();
    let process = fixture.process(false, false);
    let grouping = NumeratorGrouping::UpToScalar(crate::GraphGroupingOptions::default());
    let canonical =
        [(false, false), (false, true), (true, false), (true, true)].map(|(reverse, storage)| {
            let source = fixture.graph(false, reverse, storage);
            let sorted = process
                .canonicalize_numerator_graph(&source, &fixture.model, &grouping)
                .unwrap();
            let repeated = process
                .canonicalize_numerator_graph(&sorted, &fixture.model, &grouping)
                .unwrap();
            assert_eq!(sorted.nodes(), repeated.nodes());
            assert_eq!(
                fixture.chart(&sorted, false),
                fixture.chart(&repeated, false)
            );
            process
                .normalize_fermion_flows(&sorted, &fixture.model)
                .unwrap()
                .graph
        });
    for graph in &canonical[1..] {
        assert_eq!(
            canonical[0]
                .nodes()
                .iter()
                .map(|node| &node.data)
                .collect::<Vec<_>>(),
            graph
                .nodes()
                .iter()
                .map(|node| &node.data)
                .collect::<Vec<_>>()
        );
    }
    assert_eq!(
        fixture.chart(&canonical[0], false),
        fixture.chart(&canonical[2], true)
    );
    for index in [0, 2] {
        assert_eq!(
            fixture.chart(&canonical[index], false),
            fixture.chart(&canonical[index + 1], false)
        );
    }
}

#[test]
fn cp_preserves_neutral_adjoint_ghost_derivative_flow() {
    let fixture = GhostChart::new();
    let grouping = NumeratorGrouping::UpToScalar(crate::GraphGroupingOptions::default());
    let canonical = |reverse, storage, cp| {
        let process = fixture.process(true, cp);
        let source = fixture.graph(true, reverse, storage);
        process
            .canonicalize_numerator_graph(&source, &fixture.model, &grouping)
            .unwrap()
    };
    assert_ne!(
        fixture.chart(&canonical(false, false, false), false),
        fixture.chart(&canonical(true, false, false), false)
    );
    let process = fixture.process(true, true);
    let normalized =
        [(false, false), (false, true), (true, false), (true, true)].map(|(reverse, storage)| {
            let sorted = canonical(reverse, storage, true);
            let normalized = process
                .normalize_fermion_flows(&sorted, &fixture.model)
                .unwrap();
            let repeated = process
                .normalize_fermion_flows(&normalized.graph, &fixture.model)
                .unwrap();
            assert!(!normalized.signs.external_ordering_negative);
            assert!(!repeated.signs.external_ordering_negative);
            assert_eq!(normalized.graph.nodes(), repeated.graph.nodes());
            assert_eq!(
                fixture.chart(&normalized.graph, false),
                fixture.chart(&repeated.graph, false)
            );
            normalized.graph
        });
    // CP reverses the quark pentagon twice, while retaining the derivative
    // momentum carried by each ghost edge.
    for graph in &normalized[1..] {
        assert_eq!(
            normalized[0]
                .nodes()
                .iter()
                .map(|node| &node.data)
                .collect::<Vec<_>>(),
            graph
                .nodes()
                .iter()
                .map(|node| &node.data)
                .collect::<Vec<_>>()
        );
        assert_eq!(
            fixture.chart(&normalized[0], false),
            fixture.chart(graph, false)
        );
    }
}
