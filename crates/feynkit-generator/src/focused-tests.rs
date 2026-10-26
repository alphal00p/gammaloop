// Included inside generation.rs's existing #[cfg(test)] mod tests.
mod paired_external_momenta_regressions {
    use super::*;
    use std::sync::Arc;

    fn resolved_pair_process(
        model: &Model,
        kind: GenerationType,
        left_right: bool,
    ) -> ResolvedProcess {
        ResolvedProcess::new(
            model,
            &Process::new([-1_i64, 1_i64], [22_i64]),
            kind,
            &GenerationOptions::default().symmetrize_left_right(left_right),
        )
        .unwrap()
    }

    fn charged_pair_frame(
        model: &Model,
        fixed: bool,
        mirror: bool,
    ) -> Graph<ColoredNode, EdgeColor> {
        let external_rows = if fixed {
            [
                (1, ExternalState::Incoming, 1),
                (0, ExternalState::Incoming, -1),
                (3, ExternalState::Outgoing, 1),
                (2, ExternalState::Outgoing, -1),
            ]
        } else {
            [
                (3, ExternalState::Outgoing, 1),
                (2, ExternalState::Outgoing, -1),
                (1, ExternalState::Incoming, 1),
                (0, ExternalState::Incoming, -1),
            ]
        };
        let edge_pdgs = if fixed {
            [1, -1, -1, 1]
        } else {
            [-1, 1, 1, -1]
        };
        let mut graph = Graph::new();
        for (index, state, pdg) in external_rows {
            graph.add_node(ColoredNode::External(ExternalNode {
                index,
                state,
                particle: particle(model, pdg),
            }));
        }
        graph.add_node(ColoredNode::Interaction(vertex(model, "V")));
        graph.add_node(ColoredNode::Interaction(vertex(model, "V")));
        let (first_attachment, second_attachment) = if mirror { (5, 4) } else { (4, 5) };
        for (source, pdg) in edge_pdgs.into_iter().enumerate() {
            let target = if source < 2 {
                first_attachment
            } else {
                second_attachment
            };
            graph
                .add_edge(source, target, true, edge(model, pdg))
                .unwrap();
        }
        graph.add_edge(4, 5, false, edge(model, 22)).unwrap();
        graph
    }

    fn assert_non_assignment_data_unchanged(
        before: &Graph<ColoredNode, EdgeColor>,
        after: &Graph<ColoredNode, EdgeColor>,
    ) {
        assert_eq!(before.nodes().len(), after.nodes().len());
        assert_eq!(before.edges().len(), after.edges().len());
        for (old, new) in before.nodes().iter().zip(after.nodes()) {
            assert_eq!(
                old.edges, new.edges,
                "incident-edge order must remain fixed"
            );
            assert_eq!(old.valence, new.valence);
            match (&old.data, &new.data) {
                (ColoredNode::External(old), ColoredNode::External(new)) => {
                    assert_eq!(
                        old.particle, new.particle,
                        "physical external species must not be conjugated"
                    );
                }
                _ => assert_eq!(old.data, new.data, "interaction vertices must remain fixed"),
            }
        }
        for (old, new) in before.edges().iter().zip(after.edges()) {
            assert_eq!(old.vertices, new.vertices, "no endpoint or edge reordering");
            assert_eq!(old.directed, new.directed);
            assert_eq!(old.data.direction, new.data.direction);
            if !matches!(
                before.nodes()[old.vertices.0].data,
                ColoredNode::External(_)
            ) {
                assert_eq!(
                    old.data, new.data,
                    "internal particle/color data must remain fixed"
                );
            }
        }
    }

    fn assert_native_half_edge_and_ufo_slot_consistency(
        model: &Arc<Model>,
        graph: &Graph<ColoredNode, EdgeColor>,
        mirror: bool,
    ) {
        // Invoke the existing owner slot matcher, not numerator construction or
        // numerical evaluation. V orders its fields [f~, f, a]. Conjugating the
        // external all-incoming species must exchange the two charged slots.
        let generator = Generator::new(Arc::clone(model));
        let rule = model.vertex_rule_by_id(vertex(model, "V")).unwrap();
        let charged_slots: [[usize; 2]; 2] = if mirror {
            [[2, 3], [1, 0]]
        } else {
            [[1, 0], [2, 3]]
        };
        for (offset, vertex) in [4, 5].into_iter().enumerate() {
            let legs = generator
                .interaction_half_edges(graph, vertex, rule)
                .unwrap();
            assert_eq!(
                legs.iter().map(|leg| leg.edge).collect::<Vec<_>>(),
                vec![charged_slots[offset][0], charged_slots[offset][1], 4]
            );
            for (slot, leg) in legs.iter().enumerate() {
                let edge = &graph.edges()[leg.edge];
                let expected_flow = if slot == 2 && vertex == 4 {
                    Flow::Source
                } else {
                    Flow::Sink
                };
                assert_eq!(leg.flow, expected_flow);
                assert_eq!(
                    leg.hedge,
                    2 * leg.edge + usize::from(leg.flow == Flow::Sink)
                );
                let attachment = if leg.flow == Flow::Source {
                    edge.vertices.0
                } else {
                    edge.vertices.1
                };
                assert_eq!(attachment, vertex);
                let particle = model.particle_by_id(rule.particles[slot]).unwrap();
                assert_eq!(leg.spin, particle.spin);
                assert_eq!(leg.color, particle.color);
            }
        }
    }

    #[test]
    fn paired_external_momenta_sorted_roles_conjugate_only_external_edge_species() {
        let model = Arc::new(fermion_model());
        let resolved = resolved_pair_process(&model, GenerationType::CrossSection, true);
        let mut graph = charged_pair_frame(&model, false, false);
        let original = graph.clone();
        resolved
            .canonicalize_paired_external_momenta_assignment(&mut graph, &model)
            .unwrap();
        assert_non_assignment_data_unchanged(&original, &graph);
        assert_eq!(graph, charged_pair_frame(&model, true, false));
        assert_native_half_edge_and_ufo_slot_consistency(&model, &graph, false);

        let tags = graph
            .nodes()
            .iter()
            .filter_map(|node| match &node.data {
                ColoredNode::External(external) => Some(external.index + 1),
                ColoredNode::Interaction(_) => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(
            tags,
            vec![2, 1, 4, 3],
            "legacy [4,3,2,1] paired labels follow rebuilt edge order"
        );
        let pdgs = graph
            .edges()
            .iter()
            .map(|edge| model.particle_by_id(edge.data.particle).unwrap().pdg_code)
            .collect::<Vec<_>>();
        assert_eq!(pdgs, vec![1, -1, -1, 1, 22]);
        for edge in &graph.edges()[..4] {
            let ColoredNode::External(external) = &graph.nodes()[edge.vertices.0].data else {
                panic!("fixture external edges must be external-first");
            };
            assert_eq!(
                model.particle_by_id(edge.data.particle).unwrap().pdg_code,
                all_incoming_particle(external, &model).unwrap().pdg_code
            );
        }
    }

    #[test]
    fn paired_external_momenta_mirrored_side_counterparts_reach_a_fixed_assignment() {
        let model = Arc::new(fermion_model());
        let resolved = resolved_pair_process(&model, GenerationType::CrossSection, true);
        for mirror in [false, true] {
            let mut side_reversed = charged_pair_frame(&model, false, mirror);
            let mut already_fixed = charged_pair_frame(&model, true, mirror);
            let reversed_input = side_reversed.clone();
            let fixed_input = already_fixed.clone();
            resolved
                .canonicalize_paired_external_momenta_assignment(&mut side_reversed, &model)
                .unwrap();
            resolved
                .canonicalize_paired_external_momenta_assignment(&mut already_fixed, &model)
                .unwrap();
            assert_non_assignment_data_unchanged(&reversed_input, &side_reversed);
            assert_non_assignment_data_unchanged(&fixed_input, &already_fixed);
            assert_eq!(side_reversed, already_fixed);
            assert_eq!(
                already_fixed, fixed_input,
                "already assigned counterpart is unchanged"
            );
            assert_native_half_edge_and_ufo_slot_consistency(&model, &side_reversed, mirror);
            let once = side_reversed.clone();
            resolved
                .canonicalize_paired_external_momenta_assignment(&mut side_reversed, &model)
                .unwrap();
            assert_eq!(
                side_reversed, once,
                "fixed edge order makes role restoration idempotent"
            );
        }
    }

    #[test]
    fn paired_external_momenta_scope_gates_and_invalid_metadata_do_not_mutate_graphs() {
        let model = Arc::new(fermion_model());
        for (kind, left_right) in [
            (GenerationType::Amplitude, true),
            (GenerationType::CrossSection, false),
        ] {
            let resolved = resolved_pair_process(&model, kind, left_right);
            let mut graph = charged_pair_frame(&model, false, false);
            let before = graph.clone();
            resolved
                .canonicalize_paired_external_momenta_assignment(&mut graph, &model)
                .unwrap();
            assert_eq!(graph, before);
        }
        let vacuum = ResolvedProcess::new(
            &model,
            &Process::new(Vec::<i64>::new(), Vec::<i64>::new()),
            GenerationType::CrossSection,
            &GenerationOptions::default().symmetrize_left_right(true),
        )
        .unwrap();
        let mut vacuum_graph = Graph::new();
        vacuum_graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        vacuum_graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        vacuum_graph.add_edge(0, 1, true, edge(&model, 1)).unwrap();
        vacuum_graph.add_edge(1, 0, true, edge(&model, 1)).unwrap();
        vacuum_graph
            .add_edge(0, 1, false, edge(&model, 22))
            .unwrap();
        let before = vacuum_graph.clone();
        vacuum
            .canonicalize_paired_external_momenta_assignment(&mut vacuum_graph, &model)
            .unwrap();
        assert_eq!(vacuum_graph, before);

        let resolved = resolved_pair_process(&model, GenerationType::CrossSection, true);
        let mut malformed = charged_pair_frame(&model, false, false);
        let mut duplicate = match &malformed.nodes()[0].data {
            ColoredNode::External(external) => external.clone(),
            ColoredNode::Interaction(_) => unreachable!(),
        };
        duplicate.index = 2;
        malformed.set_node_data(0, ColoredNode::External(duplicate));
        let before = malformed.clone();
        assert!(matches!(
            resolved.canonicalize_paired_external_momenta_assignment(&mut malformed, &model),
            Err(GenerationError::GraphConstruction(_))
        ));
        assert_eq!(
            malformed, before,
            "failed validation must not commit partial assignments"
        );
    }
}
