//! Physics invariants retained from GammaLoop's completed soft foundation.

use feynkit_graph::expressions::evaluate_overall_factor;
use serde_json::Value;
use symbolica::function;

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
fn grouped_amplitude_orbits_share_the_canonical_symmetrized_external_frame() {
    let model = Model::phi4();
    let particle = model.particle_id("phi").unwrap();
    let interaction = model.vertex_rule_id(&model.vertex_rules()[0].name).unwrap();
    let edge = EdgeColor {
        particle,
        direction: None,
    };
    // Exchanging the two incoming attachments of this asymmetric graph gives
    // a distinct labelled graph, as Symbolica may return for another orbit.
    let graph = |incoming: [usize; 2]| {
        let mut graph = Graph::new();
        for _ in 0..4 {
            graph.add_node(ColoredNode::Interaction(interaction));
        }
        for (index, vertex) in incoming.into_iter().chain([2, 3]).enumerate() {
            let state = if index < 2 {
                ExternalState::Incoming
            } else {
                ExternalState::Outgoing
            };
            let external = graph.add_node(ColoredNode::External(ExternalNode {
                index,
                state,
                particle,
            }));
            let (source, target) = if index < 2 {
                (external, vertex)
            } else {
                (vertex, external)
            };
            graph.add_edge(source, target, false, edge).unwrap();
        }
        for (source, target) in [(0, 1), (1, 2), (2, 3), (0, 3), (0, 2)] {
            graph.add_edge(source, target, false, edge).unwrap();
        }
        graph
    };
    let (direct, exchanged) = (graph([0, 1]), graph([1, 0]));
    let process = Process::new(["phi", "phi"], ["phi", "phi"]);
    let grouping = NumeratorGrouping::UpToScalar(Default::default());
    for (symmetrized, grouping) in [
        (true, &grouping),
        (false, &grouping),
        (true, &NumeratorGrouping::None),
    ] {
        let options = GenerationOptions::default().symmetrize_initial(symmetrized);
        let resolved =
            ResolvedProcess::new(&model, &process, GenerationType::Amplitude, &options).unwrap();
        let canonical = |graph: &Graph<ColoredNode, EdgeColor>| {
            resolved
                .canonicalize_numerator_graph(graph, &model, grouping)
                .unwrap()
        };
        let frame = canonical(&direct);
        assert_eq!(frame, canonical(&frame));
        if symmetrized && !matches!(grouping, NumeratorGrouping::None) {
            assert_eq!(frame, canonical(&exchanged));
        } else {
            // Distinguishable legs, and diagrams that are not numerator
            // grouped, keep the labelling they were generated with.
            assert_ne!(frame, canonical(&exchanged));
        }
    }
}

/// Canonicalize `graph` for a numerator-grouped amplitude of `process`.
fn grouped_amplitude_frame(
    model: &Model,
    process: &Process,
    options: &GenerationOptions,
    graph: &Graph<ColoredNode, EdgeColor>,
) -> Graph<ColoredNode, EdgeColor> {
    ResolvedProcess::new(model, process, GenerationType::Amplitude, options)
        .unwrap()
        .canonicalize_numerator_graph(
            graph,
            model,
            &NumeratorGrouping::UpToScalar(Default::default()),
        )
        .unwrap()
}

/// An amplitude graph whose external leg `index` attaches to `attachments[index]`.
fn attached_amplitude_graph(
    rules: &[VertexRuleId],
    internal: &[(usize, usize)],
    legs: &[(ExternalState, ParticleId)],
    attachments: &[usize],
) -> Graph<ColoredNode, EdgeColor> {
    let mut graph = Graph::new();
    for rule in rules {
        graph.add_node(ColoredNode::Interaction(*rule));
    }
    for (index, ((state, particle), vertex)) in legs.iter().zip(attachments).enumerate() {
        let external = graph.add_node(ColoredNode::External(ExternalNode {
            index,
            state: *state,
            particle: *particle,
        }));
        let (source, target) = match state {
            ExternalState::Incoming => (external, *vertex),
            ExternalState::Outgoing => (*vertex, external),
        };
        let edge = EdgeColor {
            particle: *particle,
            direction: None,
        };
        graph.add_edge(source, target, false, edge).unwrap();
    }
    for (source, target) in internal {
        let particle = legs[0].1;
        let edge = EdgeColor {
            particle,
            direction: None,
        };
        graph.add_edge(*source, *target, false, edge).unwrap();
    }
    graph
}

#[test]
fn reduced_skeleton_ties_are_resolved_by_the_complete_colored_graph() {
    // Two separately named phi^3 interactions are invisible to the reduced
    // skeleton, whose automorphism exchanges the two incoming branches.
    let mut definition: Value = serde_json::from_str(&Model::phi3().to_json().unwrap()).unwrap();
    let mut second = definition["vertex_rules"][0].clone();
    second["name"] = "second_phi3_rule".into();
    definition["vertex_rules"]
        .as_array_mut()
        .unwrap()
        .push(second);
    let model = Model::from_json(&definition.to_string()).unwrap();
    let phi = model.particle_id("phi").unwrap();
    let rules = [
        model.vertex_rule_id_at(0).unwrap(),
        model.vertex_rule_id_at(1).unwrap(),
    ];
    let legs = [
        (ExternalState::Incoming, phi),
        (ExternalState::Incoming, phi),
        (ExternalState::Outgoing, phi),
    ];
    let graph = |attachments: [usize; 3]| {
        attached_amplitude_graph(
            &[rules[0], rules[1], rules[0]],
            &[(0, 1), (1, 2), (0, 2)],
            &legs,
            &attachments,
        )
    };
    let process = Process::new(["phi", "phi"], ["phi"]);
    let options = GenerationOptions::default().symmetrize_initial(true);
    let frame = |graph: &Graph<ColoredNode, EdgeColor>| {
        grouped_amplitude_frame(&model, &process, &options, graph)
    };
    let direct = frame(&graph([0, 1, 2]));
    assert_eq!(direct, frame(&graph([1, 0, 2])));
    assert_eq!(direct, frame(&direct));
}

#[test]
fn interchangeable_three_cycles_reach_one_external_frame() {
    let model = Model::phi3();
    let phi = model.particle_id("phi").unwrap();
    let rule = model.vertex_rule_id_at(0).unwrap();
    let legs = [
        (ExternalState::Incoming, phi),
        (ExternalState::Incoming, phi),
        (ExternalState::Incoming, phi),
        (ExternalState::Outgoing, phi),
    ];
    let graph = |attachments: [usize; 4]| {
        attached_amplitude_graph(
            &[rule; 5],
            &[(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (2, 4)],
            &legs,
            &attachments,
        )
    };
    let process = Process::new(["phi", "phi", "phi"], ["phi"]);
    let options = GenerationOptions::default().symmetrize_initial(true);
    let frame = |graph: &Graph<ColoredNode, EdgeColor>| {
        grouped_amplitude_frame(&model, &process, &options, graph)
    };
    let direct = frame(&graph([0, 1, 2, 3]));
    for attachments in [[1, 2, 0, 3], [2, 0, 1, 3]] {
        let cycled = frame(&graph(attachments));
        assert_eq!(direct, cycled);
        assert_eq!(cycled, frame(&cycled));
    }
}

#[test]
fn mixed_external_classes_keep_their_distinguishable_legs() {
    // Two photon classes are interchangeable; the unsymmetrized electrons
    // are not, so moving one of them is a different diagram.
    let model = Model::qed();
    let photon = model.particle_id("a").unwrap();
    let electron = model.particle_id("e-").unwrap();
    let rule = model.vertex_rule_id_at(0).unwrap();
    let legs = [
        (ExternalState::Incoming, photon),
        (ExternalState::Incoming, photon),
        (ExternalState::Incoming, electron),
        (ExternalState::Outgoing, photon),
        (ExternalState::Outgoing, photon),
        (ExternalState::Outgoing, electron),
    ];
    let graph = |attachments: [usize; 6]| {
        attached_amplitude_graph(
            &[rule; 5],
            &[(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (0, 2)],
            &legs,
            &attachments,
        )
    };
    let process = Process::new(["a", "a", "e-"], ["a", "a", "e-"]);
    let options = GenerationOptions::default()
        .symmetrize_initial(true)
        .symmetrize_final(true);
    let frame = |graph: &Graph<ColoredNode, EdgeColor>| {
        grouped_amplitude_frame(&model, &process, &options, graph)
    };
    let direct = frame(&graph([0, 1, 2, 3, 4, 1]));
    assert_eq!(direct, frame(&direct));
    for attachments in [[1, 0, 2, 3, 4, 1], [0, 1, 2, 4, 3, 1], [1, 0, 2, 4, 3, 1]] {
        assert_eq!(direct, frame(&graph(attachments)));
    }
    assert_ne!(direct, frame(&graph([0, 1, 2, 3, 4, 3])));
}

#[test]
fn interchangeable_leg_frames_need_no_permutation_enumeration() {
    // Ten interchangeable legs, two on each vertex of a path. Enumerating
    // their relabellings would canonicalize 10! = 3,628,800 skeletons.
    let model = Model::phi4();
    let phi = model.particle_id("phi").unwrap();
    let rule = model.vertex_rule_id_at(0).unwrap();
    let mut legs = vec![(ExternalState::Incoming, phi); 10];
    legs.push((ExternalState::Outgoing, phi));
    let graph = |attachments: &[usize]| {
        attached_amplitude_graph(
            &[rule; 5],
            &[(0, 1), (1, 2), (2, 3), (3, 4)],
            &legs,
            attachments,
        )
    };
    let process = Process::new(["phi"; 10], ["phi"]);
    let options = GenerationOptions::default().symmetrize_initial(true);
    let frame = |graph: &Graph<ColoredNode, EdgeColor>| {
        grouped_amplitude_frame(&model, &process, &options, graph)
    };
    let attachments = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0];
    let direct = frame(&graph(&attachments));
    assert_eq!(direct, frame(&direct));
    for shift in [1, 3, 7] {
        let mut relabelled = attachments;
        relabelled[..10].rotate_left(shift);
        assert_eq!(direct, frame(&graph(&relabelled)));
    }
}

#[test]
fn symmetrized_fermion_relabelling_carries_no_ordering_sign() {
    // With every leg symmetrized, the t- and u-channel exchanges are one
    // orbit and its external ordering sign is omitted by construction. The
    // grouped representative is a singleton, so its weight is its own.
    let model = Arc::new(Model::qed());
    let process = Process::new(["e-", "e-"], ["e-", "e-"]);
    let options = GenerationOptions::default()
        .symmetrize_initial(true)
        .symmetrize_final(true)
        .symmetrize_external_fermions(true)
        .threads(1);
    let weights = |options: &GenerationOptions| {
        process
            .generate_diagrams(Arc::clone(&model), options)
            .unwrap()
            .diagrams
            .iter()
            .map(|diagram| evaluate_overall_factor(diagram.overall_factor().as_view()))
            .collect::<Vec<_>>()
    };
    let ungrouped = weights(&options);
    assert_eq!(ungrouped.len(), 1);
    assert_eq!(
        ungrouped,
        weights(
            &options
                .clone()
                .numerator_grouping(NumeratorGrouping::UpToScalar(Default::default()))
        )
    );
}

#[test]
fn symmetrized_light_by_light_orientations_group_in_one_frame() {
    // An electron box with four labelled photons has six oriented cyclic
    // orders. Exchanging the two incoming photons leaves three orbits: one
    // closed under loop reversal, and two that are each other's reversal.
    // Charge conjugation gives a reversed four-photon loop the factor
    // (-1)^4, so those two orbits match with ratio one, but only once their
    // photon labels, momenta and polarizations share one frame.
    let model = Arc::new(Model::qed());
    let veto = model
        .particles()
        .iter()
        .map(|particle| particle.pdg_code)
        .filter(|pdg| ![11, 22].contains(&pdg.abs()))
        .map(ParticleSelector::Pdg)
        .collect();
    let process = Process::new(["a", "a"], ["a", "a"]);
    let options = GenerationOptions::default()
        .with_loop_count(1, 1)
        .unwrap()
        .with_graph_filter(GenerationFilter::ParticleVeto(veto))
        .threads(1);
    let grouping = NumeratorGrouping::UpToScalar(Default::default());
    for (symmetrized, orbits, groups) in [(false, 6, 3), (true, 3, 2)] {
        let options = options.clone().symmetrize_initial(symmetrized);
        let ungrouped = process
            .generate_diagrams(Arc::clone(&model), &options)
            .unwrap();
        assert_eq!(
            ungrouped.diagrams.len(),
            orbits,
            "symmetrized={symmetrized}"
        );
        let grouped = process
            .generate_diagrams(
                Arc::clone(&model),
                &options.clone().numerator_grouping(grouping.clone()),
            )
            .unwrap();
        assert_eq!(grouped.diagrams.len(), groups, "symmetrized={symmetrized}");
        let members = grouped
            .groups
            .iter()
            .flat_map(|group| &group.members)
            .collect::<Vec<_>>();
        assert_eq!(members.len(), orbits, "symmetrized={symmetrized}");
        assert!(
            members.iter().all(|member| member.ratio == Atom::num(1)),
            "symmetrized={symmetrized}"
        );
    }
}

#[test]
fn identical_incoming_scalar_born_weights_only_divide_by_final_state_symmetry() {
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

fn amplitude_weight(
    model: &Arc<Model>,
    process: &Process,
    options: &GenerationOptions,
) -> (usize, Atom) {
    let result = process
        .generate_diagrams(Arc::clone(model), options)
        .unwrap();
    let weight = result.diagrams.iter().fold(Atom::num(0), |sum, diagram| {
        sum + diagram.overall_factor().clone()
    });
    (result.diagrams.len(), weight)
}

#[test]
fn symmetrized_amplitude_orbits_carry_their_permutation_multiplicity() {
    let model = Arc::new(Model::qed());
    let process = Process::new(["a", "a"], ["e-", "e+"]);
    let options = GenerationOptions::default()
        .symmetrize_initial(true)
        .threads(1);
    // The t- and u-channel exchanges form one orbit of the photon swap.
    let (count, weight) = amplitude_weight(&model, &process, &options);
    assert_eq!(count, 1);
    assert_eq!(evaluate_overall_factor(weight.as_view()), Atom::num(2));
}

#[test]
fn symmetrized_amplitude_weights_sum_to_the_labelled_weights() {
    // Symmetrized external fermions drop the ordering sign, which differs
    // between the labelled diagrams of one orbit; compare without it there.
    let ordering_sign = function!(
        FACTOR_HEADS.external_fermion_ordering_sign,
        Atom::var(symbol!("x_"))
    )
    .to_pattern();
    let qed = Arc::new(Model::qed());
    let cases = [
        (
            Arc::new(Model::phi4()),
            Process::new(["phi", "phi"], ["phi", "phi"]),
            2,
            false,
        ),
        (
            Arc::clone(&qed),
            Process::new(["a", "a"], ["e-", "e+"]),
            1,
            false,
        ),
        // Two species per state: both photons and electrons are permuted
        // only when external fermions are symmetrized.
        (
            Arc::clone(&qed),
            Process::new(["e-", "e-"], ["e-", "e-", "a", "a"]),
            0,
            false,
        ),
        (
            Arc::clone(&qed),
            Process::new(["e-", "e-"], ["e-", "e-", "a", "a"]),
            0,
            true,
        ),
    ];
    for (model, process, loops, fermions) in cases {
        let options = GenerationOptions::default()
            .with_loop_count(0, loops)
            .unwrap()
            .symmetrize_external_fermions(fermions)
            .threads(1);
        let (labelled_count, labelled) = amplitude_weight(&model, &process, &options);
        for (initial, final_state) in [(true, false), (false, true), (true, true)] {
            let symmetrized_options = options
                .clone()
                .symmetrize_initial(initial)
                .symmetrize_final(final_state);
            let (count, symmetrized) = amplitude_weight(&model, &process, &symmetrized_options);
            let (labelled, symmetrized) = if fermions {
                (
                    labelled
                        .replace(ordering_sign.clone())
                        .with(Atom::num(1).to_pattern()),
                    symmetrized
                        .replace(ordering_sign.clone())
                        .with(Atom::num(1).to_pattern()),
                )
            } else {
                (labelled.clone(), symmetrized)
            };
            assert!(count <= labelled_count);
            assert_eq!(
                evaluate_overall_factor(symmetrized.as_view()).expand(),
                evaluate_overall_factor(labelled.as_view()).expand(),
                "{process:?} up to {loops} loops, fermions={fermions}, initial={initial}, final={final_state}"
            );
        }
    }
}
#[test]
fn external_orbits_merge_mixed_cubic_rules_only_when_the_states_permit() {
    let mut definition = serde_json::to_value(Model::phi3()).unwrap();
    let mut rule_b = definition["vertex_rules"][0].clone();
    rule_b["name"] = "V_B".into();
    definition["vertex_rules"]
        .as_array_mut()
        .unwrap()
        .push(rule_b);
    let model = Arc::new(Model::from_json(&definition.to_string()).unwrap());
    let rule_a = model.vertex_rule_id("V_SCALAR").unwrap();
    let process = Process::new(["phi", "phi"], ["phi", "phi"]);
    let inventory = |symmetrize| {
        let options = GenerationOptions::default()
            .symmetrize_initial(symmetrize)
            .symmetrize_final(symmetrize)
            .numerator_grouping(NumeratorGrouping::None)
            .threads(1);
        let generated = process
            .generate_diagrams(Arc::clone(&model), &options)
            .unwrap();
        let mut orbits = BTreeMap::<usize, (usize, Atom)>::new();
        for diagram in generated.diagrams {
            assert_eq!(diagram.vertices().count(), 2);
            let mut a_vertices = diagram
                .vertices()
                .filter(|(_, v)| v.interaction == Some(rule_a));
            let Some((vertex_a, _)) = a_vertices.next() else {
                continue;
            };
            if a_vertices.next().is_some() {
                continue;
            }
            // With exactly one interaction of each kind, this count classifies
            // the tree up to permutations within the incoming and outgoing pairs.
            let incoming_at_a = diagram
                .edges()
                .filter(|(_, endpoints, e)| {
                    e.external
                        .as_ref()
                        .is_some_and(|leg| leg.state == ExternalState::Incoming)
                        && (endpoints.source == Some(vertex_a)
                            || endpoints.target == Some(vertex_a))
                })
                .count();
            let (count, weight) = orbits.entry(incoming_at_a).or_insert((0, Atom::num(0)));
            *count += 1;
            *weight += evaluate_overall_factor(diagram.overall_factor().as_view());
        }
        orbits
    };
    let labelled = inventory(false);
    let symmetrized = inventory(true);
    for (incoming_at_a, labelled_count) in [(0, 1), (1, 4), (2, 1)] {
        assert_eq!(labelled[&incoming_at_a].0, labelled_count);
        assert_eq!(
            symmetrized[&incoming_at_a].0, 1,
            "incoming legs at A: {incoming_at_a}"
        );
        assert_eq!(symmetrized[&incoming_at_a].1, labelled[&incoming_at_a].1);
    }
}

#[test]
fn generated_ghost_loop_has_one_statistics_minus() {
    let model = Arc::new(
        Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap(),
    );
    let process = Process::new(["g"], ["g"]);
    let options = |fermion_loops| {
        GenerationOptions::default()
            .with_loop_count(1, 1)
            .unwrap()
            .threads(1)
            .with_graph_filter(GenerationFilter::VertexAllow(vec![VertexSelector::Name(
                "V_35".to_owned(),
            )]))
            // The public filter excludes ghost loops without changing their statistics.
            .with_graph_filter(GenerationFilter::FermionLoopCountRange((
                fermion_loops,
                fermion_loops,
            )))
    };
    let generated = Generator::new(Arc::clone(&model))
        .generate(&process, &options(0), GenerationType::Amplitude)
        .unwrap();
    let [diagram] = generated.diagrams.as_slice() else {
        panic!("expected exactly one ghost vacuum-polarization graph");
    };
    assert_eq!(diagram.loop_count(), 1);
    assert_eq!(
        evaluate_overall_factor(diagram.overall_factor().as_view()),
        Atom::num(-1)
    );
    assert!(
        Generator::new(model)
            .generate(&process, &options(1), GenerationType::Amplitude)
            .unwrap()
            .diagrams
            .is_empty()
    );
}

/// A cubic charged cycle, using the shared fixture for its external legs.
fn charged_amplitude_cycle(
    model: &Model,
    particle: ParticleId,
    external: ParticleId,
    attachments: &[usize],
    reverse_first: bool,
    conjugate_first: bool,
) -> Graph<ColoredNode, EdgeColor> {
    let antiparticle = model.particle_by_id(particle).unwrap().antiparticle;
    let rule = model
        .vertex_rules()
        .iter()
        .position(|rule| {
            rule.particles.len() == 3
                && [particle, antiparticle, external]
                    .iter()
                    .all(|particle| rule.particles.contains(particle))
        })
        .map(|index| model.vertex_rule_id_at(index).unwrap())
        .unwrap();
    let legs = (0..attachments.len())
        .map(|index| {
            (
                if index < 2 {
                    ExternalState::Incoming
                } else {
                    ExternalState::Outgoing
                },
                external,
            )
        })
        .collect::<Vec<_>>();
    let mut graph = attached_amplitude_graph(&vec![rule; legs.len()], &[], &legs, attachments);
    for source in 0..legs.len() {
        let target = (source + 1) % legs.len();
        let (source, target, particle) = if reverse_first && source == 0 {
            (
                target,
                source,
                if conjugate_first {
                    antiparticle
                } else {
                    particle
                },
            )
        } else {
            (source, target, particle)
        };
        graph
            .add_edge(
                source,
                target,
                true,
                EdgeColor {
                    particle,
                    direction: None,
                },
            )
            .unwrap();
    }
    graph
}

#[test]
fn directed_charged_orbit_frames_preserve_physical_edge_storage() {
    let model = Model::qed();
    let photon = model.particle_id("a").unwrap();
    // Triangles and boxes, with adjacent or opposite incoming attachments.
    for attachments in [vec![0, 1, 2], vec![0, 1, 2, 3], vec![0, 2, 1, 3]] {
        let process = Process::new(["a", "a"], vec!["a"; attachments.len() - 2]);
        let options = GenerationOptions::default().symmetrize_initial(true);
        let resolved =
            ResolvedProcess::new(&model, &process, GenerationType::Amplitude, &options).unwrap();
        for name in ["e-", "e+"] {
            let particle = model.particle_id(name).unwrap();
            for exact_particles in [false, true] {
                let grouping = NumeratorGrouping::UpToScalar(crate::GraphGroupingOptions {
                    differentiate_particle_masses_only: !exact_particles,
                    ..Default::default()
                });
                let frame = |graph: &Graph<ColoredNode, EdgeColor>| {
                    resolved
                        .canonicalize_numerator_graph(graph, &model, &grouping)
                        .unwrap()
                };
                let graph = |attachments: &[usize], reverse| {
                    charged_amplitude_cycle(&model, particle, photon, attachments, reverse, reverse)
                };
                let canonical = frame(&graph(&attachments, false));
                let mut exchanged = attachments.clone();
                exchanged.swap(0, 1);
                for (case, actual) in [
                    ("idempotence", frame(&canonical)),
                    ("incoming exchange", frame(&graph(&exchanged, false))),
                    ("reverse and conjugate", frame(&graph(&attachments, true))),
                ] {
                    assert_eq!(
                        canonical, actual,
                        "{case}, {attachments:?}, {name}, exact={exact_particles}"
                    );
                }
            }
        }
    }
}

#[test]
fn ghost_orbit_frames_preserve_derivative_arrow_direction() {
    let model = Model::standard_model();
    let ghost = model.particle_id("ghG").unwrap();
    let gluon = model.particle_id("g").unwrap();
    let process = Process::new(["g", "g"], ["g"]);
    let options = GenerationOptions::default().symmetrize_initial(true);
    let frame = |graph: &Graph<ColoredNode, EdgeColor>| {
        grouped_amplitude_frame(&model, &process, &options, graph)
    };
    let graph = |reverse, conjugate| {
        charged_amplitude_cycle(&model, ghost, gluon, &[0, 1, 2], reverse, conjugate)
    };
    let canonical = frame(&graph(false, false));
    assert_eq!(canonical, frame(&canonical));
    assert_eq!(canonical, frame(&graph(true, true)));
    // Reversing the derivative arrow alone changes the physical ghost flow.
    assert_ne!(canonical, frame(&graph(true, false)));
}
