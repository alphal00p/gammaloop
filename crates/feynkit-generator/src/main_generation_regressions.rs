//! Physics invariants retained from GammaLoop's completed soft foundation.

use super::tests::{ColoredNode, Graph};
use super::*;
use crate::{Cuts, Loops, Symmetrize, TopologyRules};

/// A cross-section request with the beam symmetrization and the left/right identification.
fn cross_section_request(process: &Process, initial: bool, left_right: bool) -> Request {
    process
        .clone()
        .cross_section(Loops::forward(0..=3))
        .with_symmetrize(Symmetrize {
            initial,
            ..Symmetrize::default()
        })
        .with_cuts(Cuts::new().with_symmetrize_left_right(left_right))
}

#[test]
fn paired_incoming_three_cycles_preserve_the_canonical_connection_frame() {
    let model = Model::phi4();
    let particle = model.particle_id("phi").unwrap();
    let interaction = model.vertex_rule_id(&model.vertex_rules()[0].name).unwrap();
    // An asymmetric attachment graph exercises the canonical coordinate map
    // independently of Feynman-rule evaluation. A three-cycle is not its own inverse.
    let relabel = |cycled: Option<bool>| {
        let mut graph = Graph::new();
        for _ in 0..6 {
            graph.add_node(ColoredNode::Interaction(interaction));
        }
        for label in 0..6 {
            // A paired relabelling cycles the incoming legs and their sewn copies together;
            // a one-sided one only the incoming legs.
            let index = match cycled {
                Some(paired) if paired || label < 3 => (label / 3) * 3 + (label + 1) % 3,
                _ => label,
            };
            let external = graph.add_node(ColoredNode::External(ExternalNode {
                index,
                state: if index < 3 {
                    ExternalState::Incoming
                } else {
                    ExternalState::Outgoing
                },
                particle,
            }));
            let (source, target) = if label < 3 {
                (external, label)
            } else {
                (label, external)
            };
            graph.add_edge(source, target, false, particle).unwrap();
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
            graph.add_edge(source, target, false, particle).unwrap();
        }
        graph.build()
    };
    let process = Process::new(["phi", "phi", "phi"], ["phi", "phi"]);
    // Exchanging identical beams is a rotation only in a 2 -> X process. The sewn-pair
    // quotient that identifies paired relabellings therefore rejects three identical beams.
    let symmetrized = ResolvedProcess::new(
        &model,
        &process,
        GenerationType::CrossSection,
        &cross_section_request(&process, true, false),
    )
    .and_then(|resolved| resolved.externals(&model).map(|_| ()));
    assert!(matches!(
        symmetrized,
        Err(GenerationError::UnsupportedBeamSymmetrization { incoming: 3 })
    ));
    for left_right in [false, true] {
        let resolved = ResolvedProcess::new(
            &model,
            &process,
            GenerationType::CrossSection,
            &cross_section_request(&process, false, left_right),
        )
        .unwrap();
        let canonical = resolved
            .canonicalize_numerator_graph(&relabel(None), &model, &Grouping::None)
            .unwrap();
        assert_eq!(
            canonical,
            resolved
                .canonicalize_numerator_graph(&canonical, &model, &Grouping::None)
                .unwrap()
        );
        assert_ne!(
            canonical,
            resolved
                .canonicalize_numerator_graph(&relabel(Some(false)), &model, &Grouping::None)
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
    let run = Run::default().threads(1);
    let request = |initial, left_right| {
        process
            .clone()
            .cross_section(Loops::forward(1..=1))
            .with_topology(TopologyRules::new().with_max_vertices(2))
            .with_symmetrize(Symmetrize {
                initial,
                ..Symmetrize::default()
            })
            .with_cuts(Cuts::new().with_symmetrize_left_right(left_right))
    };
    for initial in [false, true] {
        for left_right in [false, true] {
            let result = request(initial, left_right)
                .generate(Arc::clone(&model), &run)
                .unwrap()
                .build(
                    &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                    &run,
                )
                .unwrap();
            assert_eq!(result.diagrams().len(), 1);
            assert_eq!(result.diagrams()[0].cuts().len(), 1);
            assert_eq!(
                evaluate_overall_factor(result.diagrams()[0].overall_factor().as_view()),
                Atom::num((1, 2)),
                "initial={initial}, left_right={left_right}: a labelled two-particle flux has no extra initial-state factorial"
            );
        }
    }
    // Final-state symmetrization is an amplitude option.
    assert!(matches!(
        request(false, false)
            .with_symmetrize(Symmetrize::FINAL)
            .generate(Arc::clone(&model), &run),
        Err(GenerationError::CrossSectionFinalStateSymmetrization)
    ));
}

#[test]
fn directed_self_loop_symmetry_factors_match_wick_contractions() {
    use linnet::half_edge::generate::{
        Legs, Plan, Query, Run as EngineRun, Structure, Stub, Theory,
    };

    const NEUTRAL: u8 = 0;
    const CHARGED: u8 = 1;
    const OTHER: u8 = 2;
    let neutral = Stub::undirected(NEUTRAL);
    let charged = Stub::incoming(CHARGED);
    let other = Stub::incoming(OTHER);
    // Production consumes these graph automorphism denominators directly.
    // Directed loop ends are distinct particle/antiparticle fields; neutral
    // flips and permutations of identical propagators still count.
    for (signature, vertices, expected) in [
        (
            vec![charged.clone(), charged.clone().flipped()],
            1,
            vec![(1, 1)],
        ),
        (
            vec![
                charged.clone(),
                charged.clone().flipped(),
                other.clone(),
                other.clone().flipped(),
            ],
            1,
            vec![(2, 1)],
        ),
        (
            vec![
                charged.clone(),
                charged.clone().flipped(),
                charged.clone(),
                charged.clone().flipped(),
            ],
            1,
            vec![(2, 2)],
        ),
        (vec![neutral.clone(); 4], 1, vec![(0, 8)]),
        (
            vec![
                neutral.clone(),
                neutral.clone(),
                charged.clone(),
                charged.clone().flipped(),
            ],
            1,
            vec![(1, 2)],
        ),
        (
            vec![neutral.clone(), charged.clone(), charged.clone().flipped()],
            2,
            vec![(0, 2), (2, 2)],
        ),
    ] {
        let theory = Theory::builder()
            .vertex(signature.clone(), [])
            .build()
            .unwrap();
        let query = Query::new(&Legs::<u8, u8>::new())
            .structure(Structure {
                max_vertices: Some(vertices),
                self_loops: true,
                zero_flow_edges: true,
                ..Structure::default()
            })
            .loops(0..=2);
        let mut observed = Plan::new(&theory, &query)
            .unwrap()
            .diagrams(&EngineRun::default())
            .unwrap()
            .map(|diagram| {
                (
                    diagram
                        .edges()
                        .filter(|edge| edge.flow.is_some() && edge.vertices[0] == edge.vertices[1])
                        .count(),
                    diagram.symmetry_factor().unwrap(),
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
        let process = Process::new(
            [22_i64],
            if cross_section {
                vec![1_i64, -1]
            } else {
                vec![22]
            },
        );
        let request = if cross_section {
            process
                .clone()
                .cross_section(Loops::forward(3..=3))
                .with_cuts(Cuts::new().with_symmetrize_left_right(cp))
        } else {
            process.clone().amplitude(Loops::total(2..=2))
        };
        let kind = request.kind();
        ResolvedProcess::new(&self.model, &process, kind, &request).unwrap()
    }

    fn graph(&self, triangle: bool, reverse: bool, reverse_storage: bool) -> GeneratedGraph {
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
                .add_edge(a, b, true, self.model.particle_id_by_pdg(pdg).unwrap())
                .unwrap();
        }
        graph.build()
    }

    /// The edges as sorted `(ends, directed, particle)`, optionally with every quark reversed.
    fn chart(
        &self,
        graph: &GeneratedGraph,
        reverse_quarks: bool,
    ) -> Vec<((Endpoint, Endpoint), bool, ParticleId)> {
        let mut edges = graph
            .edges()
            .map(|edge| {
                let particle = self.model.particle_by_id(edge.particle).unwrap();
                let endpoints = if reverse_quarks && particle.is_fermion() && !particle.is_ghost() {
                    (edge.sink, edge.source)
                } else {
                    (edge.source, edge.sink)
                };
                (endpoints, edge.directed, edge.particle)
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
    let grouping = Grouping::UpToScalar(crate::GraphGroupingOptions::default());
    let canonical =
        [(false, false), (false, true), (true, false), (true, true)].map(|(reverse, storage)| {
            let source = fixture.graph(false, reverse, storage);
            let sorted = process
                .canonicalize_numerator_graph(&source, &fixture.model, &grouping)
                .unwrap();
            let repeated = process
                .canonicalize_numerator_graph(&sorted, &fixture.model, &grouping)
                .unwrap();
            assert!(sorted.vertices().eq(repeated.vertices()));
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
        assert!(canonical[0].vertices().eq(graph.vertices()));
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
    let grouping = Grouping::UpToScalar(crate::GraphGroupingOptions::default());
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
            assert!(normalized.graph.vertices().eq(repeated.graph.vertices()));
            assert_eq!(
                fixture.chart(&normalized.graph, false),
                fixture.chart(&repeated.graph, false)
            );
            normalized.graph
        });
    // CP reverses the quark pentagon twice, while retaining the derivative
    // momentum carried by each ghost edge.
    for graph in &normalized[1..] {
        assert!(normalized[0].vertices().eq(graph.vertices()));
        assert_eq!(
            fixture.chart(&normalized[0], false),
            fixture.chart(graph, false)
        );
    }
}
