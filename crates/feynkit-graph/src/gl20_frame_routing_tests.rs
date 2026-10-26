mod gl20_frame_routing_tests {
    use super::*;

    #[derive(Clone, Copy, Debug)]
    enum StoredFlow {
        Forward,
        Backward,
    }

    #[derive(Clone, Copy, Debug)]
    enum ExternalFrame {
        Paired {
            incoming: usize,
            outgoing: usize,
        },
        OneSided {
            state: ExternalState,
            attachment: usize,
        },
    }

    fn routing_fixture(frame: ExternalFrame, flow: StoredFlow) -> FeynmanDiagram {
        let model = scalar_model();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "native-root-routing");
        let vertices = (0..8)
            .map(|vertex| {
                builder.add_vertex(DiagramVertex {
                    name: format!("v{vertex}"),
                    interaction: None,
                    numerator: SymbolicTensor::infer(Atom::one()).unwrap(),
                })
            })
            .collect::<Vec<_>>();
        for connection in 0..2 {
            match frame {
                ExternalFrame::Paired { incoming, outgoing } => {
                    let initial = builder.add_generation_external(
                        format!("in{connection}"),
                        connection,
                        ExternalState::Incoming,
                        connection,
                    );
                    let final_state = builder.add_generation_external(
                        format!("out{connection}"),
                        connection + 2,
                        ExternalState::Outgoing,
                        connection,
                    );
                    builder
                        .add_edge(
                            initial,
                            vertices[incoming],
                            DiagramEdge::new(particle, true),
                        )
                        .unwrap();
                    builder
                        .add_edge(
                            vertices[outgoing],
                            final_state,
                            DiagramEdge::new(particle, true),
                        )
                        .unwrap();
                }
                ExternalFrame::OneSided { state, attachment } => {
                    let external = builder.add_generation_external(
                        format!("ext{connection}"),
                        connection,
                        state,
                        connection,
                    );
                    let (source, target) = match state {
                        ExternalState::Incoming => (external, vertices[attachment]),
                        ExternalState::Outgoing => (vertices[attachment], external),
                    };
                    builder
                        .add_edge(source, target, DiagramEdge::new(particle, true))
                        .unwrap();
                }
            }
        }
        // The retained eight-vertex topology, in its stored internal-edge
        // order. Reversing scalar carriers changes their native endpoints
        // while preserving the unoriented ordering and spanning forests.
        for (source, target) in [
            (0, 1),
            (0, 3),
            (5, 0),
            (2, 1),
            (1, 5),
            (4, 2),
            (2, 7),
            (3, 4),
            (3, 6),
            (4, 5),
        ] {
            let (source, target) = match flow {
                StoredFlow::Forward => (source, target),
                StoredFlow::Backward => (target, source),
            };
            builder
                .add_edge(
                    vertices[source],
                    vertices[target],
                    DiagramEdge::new(particle, true),
                )
                .unwrap();
        }
        let diagram = builder.canonical_loop_momentum_basis().build().unwrap();
        if let ExternalFrame::Paired { incoming, .. } = frame {
            // Retained old-master cut 0 has internal vertices [1, 2, 4] on
            // the left; native complete partitions also include attachment 7.
            // Reverse the two sides when the mirrored fixture enters at 6.
            let historical_left = [1, 2, 4, 7];
            let incoming_on_historical_left = historical_left.contains(&incoming);
            let (left, right): (Vec<DiagramHalfEdge>, Vec<DiagramHalfEdge>) =
                diagram.half_edges().partition(|half| {
                    let hedge = diagram.half_edge_id(*half).unwrap();
                    let vertex = diagram.underlying().node_id(hedge).0;
                    historical_left.contains(&vertex) == incoming_on_historical_left
                });
            let left_set = left.iter().copied().collect::<BTreeSet<_>>();
            let right_set = right.iter().copied().collect::<BTreeSet<_>>();
            assert!(left_set.is_disjoint(&right_set));
            assert_eq!(left.len() + right.len(), diagram.half_edges().count());
            assert_eq!(
                left_set.union(&right_set).copied().collect::<BTreeSet<_>>(),
                diagram.half_edges().collect::<BTreeSet<_>>()
            );
            let diagram = diagram.with_cut_partitions(vec![(left, right)]).unwrap();
            assert_eq!(
                diagram.cuts()[0]
                    .cut
                    .iter()
                    .map(|half| half.edge)
                    .collect::<BTreeSet<_>>(),
                BTreeSet::from([EdgeId(2), EdgeId(6), EdgeId(9), EdgeId(11)])
            );
            diagram
        } else {
            diagram
        }
    }

    fn assert_native_root(diagram: &FeynmanDiagram, expected: VertexId) {
        let first = diagram.half_edges().next().unwrap();
        let hedge = diagram.half_edge_id(first).unwrap();
        assert_eq!(hedge.0, 0);
        assert_eq!(VertexId(diagram.underlying().node_id(hedge).0), expected);
    }

    #[test]
    fn gl20_frame_paired_routing_follows_native_first_half_edge_in_both_flows() {
        for flow in [StoredFlow::Forward, StoredFlow::Backward] {
            let diagram = routing_fixture(
                ExternalFrame::Paired {
                    incoming: 7,
                    outgoing: 6,
                },
                flow,
            );
            assert_native_root(&diagram, VertexId(7));
            let (_, endpoints, _) = diagram.edges().next().unwrap();
            assert_eq!(endpoints.source, Some(VertexId(6)));
            assert_eq!(endpoints.target, Some(VertexId(7)));
            assert_eq!(
                diagram.loop_momentum_basis().tree_edges,
                [
                    EdgeId(2),
                    EdgeId(3),
                    EdgeId(6),
                    EdgeId(7),
                    EdgeId(8),
                    EdgeId(10),
                    EdgeId(11)
                ]
            );
            assert_eq!(
                diagram.loop_momentum_basis().loop_edges,
                [EdgeId(4), EdgeId(5), EdgeId(9)]
            );
            diagram.loop_momentum_basis().validate(&diagram).unwrap();
            diagram.validate().unwrap();

            let native_order = diagram.half_edges().collect::<Vec<_>>();
            let rerouted = diagram
                .with_first_loop_momentum_basis_in_edge_order(
                    &(2..12).map(EdgeId).collect::<Vec<_>>(),
                )
                .unwrap();
            assert_eq!(rerouted.half_edges().collect::<Vec<_>>(), native_order);
            assert_eq!(
                rerouted.loop_momentum_basis().loop_edges,
                [EdgeId(4), EdgeId(5), EdgeId(9)]
            );
        }
    }

    #[test]
    fn gl20_frame_paired_mirror_uses_its_actual_native_owner() {
        for flow in [StoredFlow::Forward, StoredFlow::Backward] {
            let diagram = routing_fixture(
                ExternalFrame::Paired {
                    incoming: 6,
                    outgoing: 7,
                },
                flow,
            );
            assert_native_root(&diagram, VertexId(6));
            let (_, endpoints, _) = diagram.edges().next().unwrap();
            assert_eq!(endpoints.source, Some(VertexId(7)));
            assert_eq!(endpoints.target, Some(VertexId(6)));
            assert_eq!(
                diagram.loop_momentum_basis().tree_edges,
                [
                    EdgeId(2),
                    EdgeId(5),
                    EdgeId(6),
                    EdgeId(8),
                    EdgeId(9),
                    EdgeId(10),
                    EdgeId(11)
                ]
            );
            assert_eq!(
                diagram.loop_momentum_basis().loop_edges,
                [EdgeId(3), EdgeId(4), EdgeId(7)]
            );
            diagram.loop_momentum_basis().validate(&diagram).unwrap();
        }
    }

    #[test]
    fn gl20_frame_one_sided_empty_and_disconnected_fallbacks_are_retained() {
        for (state, attachment, loops) in [
            (
                ExternalState::Outgoing,
                7,
                [EdgeId(4), EdgeId(5), EdgeId(9)],
            ),
            (
                ExternalState::Incoming,
                6,
                [EdgeId(3), EdgeId(4), EdgeId(7)],
            ),
        ] {
            let diagram = routing_fixture(
                ExternalFrame::OneSided { state, attachment },
                StoredFlow::Forward,
            );
            assert_native_root(&diagram, VertexId(attachment));
            assert_eq!(diagram.loop_momentum_basis().loop_edges, loops);
            diagram.loop_momentum_basis().validate(&diagram).unwrap();
        }
        let empty = FeynmanDiagram::builder(scalar_model(), "empty-routing")
            .canonical_loop_momentum_basis()
            .build()
            .unwrap();
        assert!(empty.half_edges().next().is_none());
        assert!(empty.loop_momentum_basis().loop_edges.is_empty());
        assert!(empty.loop_momentum_basis().tree_edges.is_empty());
        empty.loop_momentum_basis().validate(&empty).unwrap();

        let model = scalar_model();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "paired-disconnected-routing");
        let vertices = (0..5)
            .map(|vertex| {
                builder.add_vertex(DiagramVertex {
                    name: format!("v{vertex}"),
                    interaction: None,
                    numerator: SymbolicTensor::infer(Atom::one()).unwrap(),
                })
            })
            .collect::<Vec<_>>();
        let incoming = builder.add_generation_external("in", 0, ExternalState::Incoming, 0);
        let outgoing = builder.add_generation_external("out", 1, ExternalState::Outgoing, 0);
        builder
            .add_edge(incoming, vertices[0], DiagramEdge::new(particle, false))
            .unwrap();
        builder
            .add_edge(vertices[1], outgoing, DiagramEdge::new(particle, false))
            .unwrap();
        // A parallel-edge component, an unrooted component with a self-loop,
        // and an isolated vertex require the existing uncovered-edge fallback.
        for (source, target) in [(0, 1), (0, 1), (2, 3), (2, 2)] {
            builder
                .add_edge(
                    vertices[source],
                    vertices[target],
                    DiagramEdge::new(particle, false),
                )
                .unwrap();
        }
        let disconnected = builder.canonical_loop_momentum_basis().build().unwrap();
        // Fixture-only cut through the parallel component: the incoming
        // attachment is on the left, and every remaining half-edge is right.
        let (left, right): (Vec<DiagramHalfEdge>, Vec<DiagramHalfEdge>) =
            disconnected.half_edges().partition(|half| {
                let hedge = disconnected.half_edge_id(*half).unwrap();
                disconnected.underlying().node_id(hedge).0 == vertices[0].0
            });
        let left_set = left.iter().copied().collect::<BTreeSet<_>>();
        let right_set = right.iter().copied().collect::<BTreeSet<_>>();
        assert!(left_set.is_disjoint(&right_set));
        assert_eq!(left.len() + right.len(), disconnected.half_edges().count());
        assert_eq!(
            left_set.union(&right_set).copied().collect::<BTreeSet<_>>(),
            disconnected.half_edges().collect::<BTreeSet<_>>()
        );
        let disconnected = disconnected
            .with_cut_partitions(vec![(left, right)])
            .unwrap();
        assert_eq!(
            disconnected.cuts()[0]
                .cut
                .iter()
                .map(|half| half.edge)
                .collect::<BTreeSet<_>>(),
            BTreeSet::from([EdgeId(1), EdgeId(2)])
        );
        assert_native_root(&disconnected, VertexId(0));
        assert_eq!(
            disconnected.loop_momentum_basis().tree_edges,
            [EdgeId(2), EdgeId(3)]
        );
        assert_eq!(
            disconnected.loop_momentum_basis().loop_edges,
            [EdgeId(1), EdgeId(4)]
        );
        disconnected
            .loop_momentum_basis()
            .validate(&disconnected)
            .unwrap();
        disconnected.validate().unwrap();
    }
}
