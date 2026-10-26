mod physical_graph_order_regressions {
    use super::*;

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum ClosedFermionFlow {
        Forward,
        Reversed,
    }

    fn physical_flow_counterpart(
        model: &Model,
        flow: ClosedFermionFlow,
    ) -> Graph<ColoredNode, EdgeColor> {
        // Retained raw physical frame of the two historical flow counterparts.
        // These are semantic rules/PDGs, never model IDs or generated GL labels.
        let mut graph = Graph::new();
        for (index, state, pdg) in [
            (1, ExternalState::Incoming, 11),
            (0, ExternalState::Incoming, -11),
            (3, ExternalState::Outgoing, 11),
            (2, ExternalState::Outgoing, -11),
        ] {
            graph.add_node(ColoredNode::External(ExternalNode {
                index,
                state,
                particle: particle(model, pdg),
            }));
        }
        for name in [
            "V_141", "V_141", "V_134", "V_134", "V_137", "V_137", "V_98", "V_98",
        ] {
            graph.add_node(ColoredNode::Interaction(vertex(model, name)));
        }
        for (source, target, directed, pdg) in [
            (11, 0, true, -11),
            (1, 11, true, -11),
            (10, 2, true, 11),
            (3, 10, true, 11),
            (4, 5, false, 25),
            (4, 7, true, 6),
            (9, 4, true, 6),
            (6, 5, true, 6),
            (5, 9, true, 6),
            (8, 6, true, 6),
            (6, 11, false, 22),
            (7, 8, true, 6),
            (7, 10, false, 22),
            (8, 9, false, 21),
        ] {
            let (source, target) = if pdg == 6 && flow == ClosedFermionFlow::Reversed {
                (target, source)
            } else {
                (source, target)
            };
            graph
                .add_edge(source, target, directed, edge(model, pdg))
                .unwrap();
        }
        graph
    }

    #[test]
    fn physical_graph_order_selects_the_lower_raw_closed_fermion_flow() {
        let model = standard_model();
        let forward = physical_flow_counterpart(&model, ClosedFermionFlow::Forward);
        let reversed = physical_flow_counterpart(&model, ClosedFermionFlow::Reversed);
        let forward_before = forward.clone();
        let reversed_before = reversed.clone();
        let forward_key = LegacyPhysicalGraphOrdKey::from_graph(&forward, &model).unwrap();
        let reversed_key = LegacyPhysicalGraphOrdKey::from_graph(&reversed, &model).unwrap();

        assert_eq!(forward_key.nodes, reversed_key.nodes);
        assert_eq!(forward_key.edges[..5], reversed_key.edges[..5]);
        assert_eq!(forward_key.edges[5].vertices, (4, 7));
        assert_eq!(reversed_key.edges[5].vertices, (7, 4));
        assert!(forward_key < reversed_key);
        let mut ordered = [
            (reversed_key, ClosedFermionFlow::Reversed),
            (forward_key, ClosedFermionFlow::Forward),
        ];
        ordered.sort_by(|left, right| left.0.cmp(&right.0));
        assert_eq!(ordered[0].1, ClosedFermionFlow::Forward);
        assert_eq!(
            forward, forward_before,
            "projection must not change a graph"
        );
        assert_eq!(
            reversed, reversed_before,
            "projection must not change a graph"
        );
    }

    #[test]
    fn physical_graph_order_preserves_storage_metadata_and_legacy_key_direction() {
        let model = standard_model();
        let mut graph = Graph::new();
        graph.add_node(ColoredNode::External(ExternalNode {
            index: 0,
            state: ExternalState::Incoming,
            particle: particle(&model, -11),
        }));
        graph.add_node(ColoredNode::Interaction(vertex(&model, "V_98")));
        // The stored self-loop precedes the lexicographically smaller edge.
        graph.add_edge(1, 1, false, edge(&model, 22)).unwrap();
        graph.add_edge(0, 1, true, edge(&model, -11)).unwrap();
        let before = graph.clone();
        let key = LegacyPhysicalGraphOrdKey::from_graph(&graph, &model).unwrap();
        assert_eq!(key.nodes[0].data, (1, "external".to_owned()));
        assert_eq!(key.nodes[0].incident_edges, vec![1]);
        assert_eq!(key.nodes[0].valence, 1);
        assert_eq!(key.nodes[1].data, (0, "V_98".to_owned()));
        assert_eq!(key.nodes[1].incident_edges, vec![0, 1]);
        assert_eq!(key.nodes[1].valence, 3, "self-loops count twice in valence");
        assert_eq!(key.edges[0].vertices, (1, 1));
        assert_eq!(key.edges[0].signed_pdg, 22);
        assert_eq!(key.edges[1].vertices, (0, 1));
        assert_eq!(key.edges[1].signed_pdg, -11);
        assert!(key.edges.iter().all(|edge| edge.directed));
        assert!(
            !graph.edges()[0].directed,
            "legacy flag belongs only to the key"
        );
        assert_eq!(graph, before);

        let mut longer: Graph<ColoredNode, EdgeColor> = Graph::new();
        for _ in 0..3 {
            longer.add_node(ColoredNode::Interaction(vertex(&model, "V_98")));
        }
        assert!(
            LegacyPhysicalGraphOrdKey::from_graph(&longer, &model).unwrap() < key,
            "raw node data precede vector length; counts must not lead the key"
        );

        let mut assignment_hint = graph.edges()[0].data;
        assignment_hint.direction = Some(true);
        graph.set_edge_data(0, assignment_hint);
        assert_eq!(
            key,
            LegacyPhysicalGraphOrdKey::from_graph(&graph, &model).unwrap(),
            "the modern half-edge assignment hint is absent from legacy EdgeColor"
        );
    }

    #[test]
    fn physical_graph_order_resolves_rule_names_and_signed_pdgs_instead_of_model_ids() {
        let model = standard_model();
        let earlier_name = vertex(&model, "V_134");
        let later_name = vertex(&model, "V_98");
        assert!(earlier_name.index() > later_name.index());
        let mut earlier: Graph<ColoredNode, EdgeColor> = Graph::new();
        earlier.add_node(ColoredNode::Interaction(earlier_name));
        let mut later: Graph<ColoredNode, EdgeColor> = Graph::new();
        later.add_node(ColoredNode::Interaction(later_name));
        assert!(
            LegacyPhysicalGraphOrdKey::from_graph(&earlier, &model).unwrap()
                < LegacyPhysicalGraphOrdKey::from_graph(&later, &model).unwrap()
        );

        let negative_pdg = particle(&model, -11);
        let positive_pdg = particle(&model, 11);
        assert!(negative_pdg.index() > positive_pdg.index());
        for graph in [&mut earlier, &mut later] {
            graph.set_node_data(0, ColoredNode::Interaction(later_name));
            graph.add_node(ColoredNode::Interaction(later_name));
        }
        earlier.add_edge(0, 1, true, edge(&model, -11)).unwrap();
        later.add_edge(0, 1, true, edge(&model, 11)).unwrap();
        assert!(
            LegacyPhysicalGraphOrdKey::from_graph(&earlier, &model).unwrap()
                < LegacyPhysicalGraphOrdKey::from_graph(&later, &model).unwrap()
        );
    }

    #[test]
    fn physical_graph_order_scope_excludes_amplitudes_unsymmetrized_xs_and_vacuum() {
        let model = standard_model();
        for (kind, left_right, expected) in [
            (GenerationType::CrossSection, true, true),
            (GenerationType::CrossSection, false, false),
            (GenerationType::Amplitude, true, false),
            (GenerationType::Amplitude, false, false),
        ] {
            let resolved = ResolvedProcess::new(
                &model,
                &Process::new([-11_i64, 11_i64], [6_i64, -6_i64, 25_i64]),
                kind,
                &GenerationOptions::default().symmetrize_left_right(left_right),
            )
            .unwrap();
            assert_eq!(resolved.uses_legacy_physical_graph_order(), expected);
        }
        let vacuum = ResolvedProcess::new(
            &model,
            &Process::new(Vec::<i64>::new(), Vec::<i64>::new()),
            GenerationType::CrossSection,
            &GenerationOptions::default().symmetrize_left_right(true),
        )
        .unwrap();
        assert!(!vacuum.uses_legacy_physical_graph_order());

        let mut invalid: Graph<ColoredNode, EdgeColor> = Graph::new();
        invalid.add_node(ColoredNode::External(ExternalNode {
            index: usize::MAX,
            state: ExternalState::Incoming,
            particle: particle(&model, -11),
        }));
        assert!(matches!(
            LegacyPhysicalGraphOrdKey::from_graph(&invalid, &model),
            Err(GenerationError::ArithmeticOverflow(_))
        ));
    }
}
