use feynkit_graph::{DiagramEdge, EdgeEndpoints, EdgeId, VertexId};
use linnet::half_edge::{
    involution::{EdgeIndex, HedgePair},
    subgraph::SubSetLike,
};
use typed_index_collections::ti_vec;

use super::Graph;
use crate::{
    finalized_runtime_dot,
    graph::{
        FeynmanGraph, GraphGroup,
        parse::{IntoFinalizedRuntimeGraph, complete_group_parsing},
    },
    initialisation::test_initialise,
    momentum::sample::LoopIndex,
    processes::DotExportSettings,
    uv::uv_graph::UVE,
};

fn triangle_fixture(loop_edge: usize, loop_edge_mass: Option<i64>) -> String {
    let mass = loop_edge_mass.map_or_else(String::new, |mass| format!(" mass={mass}"));
    let loop_attribute = |edge| {
        if edge == loop_edge { " lmb_id=0" } else { "" }
    };
    format!(
        r#"digraph triangle {{
            graph [projector="1"]
            node [num="1"]
            edge [pdg=1000 num="1" dir=none]
            ext [style=invis]
            ext -> v4 [id=0 sink="{{ufo_order:0}}"]
            v5 -> ext [id=1 source="{{ufo_order:0}}"]
            v6 -> ext [id=2 source="{{ufo_order:0}}"]
            v4 -> v5 [id=3{}{} source="{{ufo_order:1}}" sink="{{ufo_order:1}}"]
            v5 -> v6 [id=4{} source="{{ufo_order:2}}" sink="{{ufo_order:1}}"]
            v6 -> v4 [id=5{} source="{{ufo_order:2}}" sink="{{ufo_order:2}}"]
        }}"#,
        loop_attribute(3),
        if loop_edge == 3 { mass.as_str() } else { "" },
        loop_attribute(4),
        loop_attribute(5),
    )
}

fn parse_triangle_fixture(loop_edge: usize, loop_edge_mass: Option<i64>) -> Graph {
    triangle_fixture(loop_edge, loop_edge_mass)
        .into_finalized_runtime_graph(&crate::utils::load_generic_model("scalars"))
        .unwrap()
}

#[test]
fn legacy_internal_order_key_recovers_pre_normalized_fermion_flow() {
    let model = crate::utils::load_generic_model("sm");
    let forward = DiagramEdge::new(model.particle_id("u").unwrap(), true);
    let reversed = DiagramEdge::new(model.particle_id("u~").unwrap(), true);

    let forward_key = forward
        .canonical_order_key(
            &model,
            EdgeId(7),
            EdgeEndpoints {
                source: Some(VertexId(2)),
                target: Some(VertexId(5)),
            },
        )
        .unwrap();
    let reversed_key = reversed
        .canonical_order_key(
            &model,
            EdgeId(7),
            EdgeEndpoints {
                source: Some(VertexId(5)),
                target: Some(VertexId(2)),
            },
        )
        .unwrap();

    assert_eq!(reversed_key, forward_key);
    assert_eq!(forward_key.2, 2);
}

fn split_forward_fixture(outgoing_momentum: &str, explicit_connection: bool) -> String {
    let connection = if explicit_connection {
        r#" is_cut="0" initial_state_connection="true""#
    } else {
        ""
    };
    format!(
        r#"digraph split_forward {{
            graph [projector="1"]
            node [num="1"]
            edge [pdg=1000 num="1" dir=none]
            ext_in [style=invis]
            ext_out [style=invis]
            ext_in -> A [id=0 lmb_rep="P(0,a___)" sink="{{ufo_order:0}}"{connection}]
            A -> B [id=1 lmb_id=0 lmb_rep="K(0,a___)" source="{{ufo_order:1}}" sink="{{ufo_order:0}}"]
            A -> B [id=2 lmb_rep="-K(0,a___)+P(0,a___)" source="{{ufo_order:2}}" sink="{{ufo_order:1}}"]
            B -> ext_out [id=3 lmb_rep="{outgoing_momentum}" source="{{ufo_order:2}}"{connection}]
        }}"#
    )
}

fn two_split_forward_connections_fixture() -> &'static str {
    r#"digraph two_split_forward_connections {
        graph [projector="1"]
        node [num="1"]
        edge [pdg=1000 num="1" dir=none]
        in0 [style=invis]
        out0 [style=invis]
        in1 [style=invis]
        out1 [style=invis]
        in0 -> A [id=0 lmb_rep="P(0,a___)" sink="{ufo_order:0}" is_cut="0" initial_state_connection="true"]
        in1 -> A [id=1 lmb_rep="P(1,a___)" sink="{ufo_order:1}" is_cut="1" initial_state_connection="true"]
        A -> B [id=2 lmb_rep="P(0,a___)+P(1,a___)" source="{ufo_order:2}" sink="{ufo_order:0}"]
        B -> out0 [id=3 lmb_rep="P(0,a___)" source="{ufo_order:1}" is_cut="0" initial_state_connection="true"]
        B -> out1 [id=4 lmb_rep="P(1,a___)" source="{ufo_order:2}" is_cut="1" initial_state_connection="true"]
    }"#
}

#[test]
fn explicitly_marked_split_forward_endpoints_are_one_connection() {
    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("scalars");
    let graph: Graph = split_forward_fixture("P(0,a___)", true)
        .into_finalized_runtime_graph(&model)
        .unwrap();

    assert!(graph.get_external_signature().is_empty());
    assert_eq!(graph.loop_momentum_basis.ext_edges.len(), 1);
    assert_eq!(graph.loop_momentum_basis.loop_edges.len(), 1);

    let serialized = graph.dot_serialize(&DotExportSettings {
        split_xs_by_initial_states: true,
        ..DotExportSettings::default()
    });
    assert_eq!(serialized.matches("is_cut=").count(), 2, "{serialized}");

    let round_trip: Graph = serialized.into_finalized_runtime_graph(&model).unwrap();
    assert!(round_trip.get_external_signature().is_empty());
    assert_eq!(round_trip.loop_momentum_basis.ext_edges.len(), 1);
}

#[test]
fn equal_external_momenta_without_explicit_markers_remain_distinct_legs() {
    test_initialise().unwrap();
    let graph: Graph = split_forward_fixture("P(0,a___)", false)
        .into_finalized_runtime_graph(&crate::utils::load_generic_model("scalars"))
        .unwrap();

    assert_eq!(graph.get_external_signature().len(), 2);
    assert_eq!(graph.loop_momentum_basis.ext_edges.len(), 2);
    assert!(graph.initial_state_cut.is_empty());
}

#[test]
fn multiple_undirected_forward_connections_survive_split_round_trip() {
    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("scalars");
    let graph: Graph = two_split_forward_connections_fixture()
        .into_finalized_runtime_graph(&model)
        .unwrap();
    assert!(graph.get_external_signature().is_empty());
    assert_eq!(graph.loop_momentum_basis.ext_edges.len(), 2);

    let serialized = graph.dot_serialize(&DotExportSettings {
        split_xs_by_initial_states: true,
        ..DotExportSettings::default()
    });
    assert_eq!(serialized.matches("is_cut=").count(), 4, "{serialized}");

    let round_trip: Graph = serialized.into_finalized_runtime_graph(&model).unwrap();
    assert!(round_trip.get_external_signature().is_empty());
    assert_eq!(round_trip.loop_momentum_basis.ext_edges.len(), 2);
}

#[test]
fn test_loop_momentum_basis() {
    test_initialise().unwrap();
    let g = parse_triangle_fixture(3, None);

    assert_eq!(g.loop_momentum_basis.loop_edges.len(), 1);
    assert_eq!(
        g.loop_momentum_basis.loop_edges[LoopIndex::from(0)],
        EdgeIndex::from(3)
    );
    assert_eq!(g.loop_momentum_basis.ext_edges.len(), 3);
}

#[test]
fn dod_rescales_only_internal_edge_qs() {
    test_initialise().unwrap();
    let g: Graph = finalized_runtime_dot!(
        digraph G {
            graph [projector="1"]
            edge [pdg=1000 num="1" dir=none]
            node [num="1"]
            ext0 [style=invis]
            ext1 [style=invis]
            ext0 -> A [id=0 sink="{ufo_order:0}"]
            B -> ext1 [id=1 source="{ufo_order:0}"]
            A -> B [id=2, lmb_id=0, source="{ufo_order:1}", sink="{ufo_order:1}", num="Q(2,spenso::mink(4,edge(2,1)))*Q(0,spenso::mink(4,edge(2,1)))"]
            A -> B [id=3, source="{ufo_order:2}", sink="{ufo_order:2}"]
        },
        "scalars"
    )
    .unwrap();

    assert_eq!(g.underlying[EdgeIndex::from(2)].dod.value, -1);
}

#[test]
fn parse_and_build_edgevec() {
    test_initialise().unwrap();
    let graph_1 = parse_triangle_fixture(3, None);
    let graph_2 = parse_triangle_fixture(3, None);

    let graph_1 = &graph_1.underlying;
    let graph_2 = &graph_2.underlying;

    let test_data = vec![true, true, true];

    let graph_1_test = graph_1.new_edgevec(|_, _, p| p.is_paired());

    let graph_2_test = graph_2.new_edgevec(|_, _, p| p.is_paired());

    let mut graph_1_iter = test_data.clone().into_iter();
    let mut graph_2_iter = test_data.clone().into_iter();

    let graph_1_test_2 = graph_1
        .new_edgevec_from_iter(graph_1.iter_edges().map(|(pair, _, _)| {
            if matches!(pair, HedgePair::Paired { .. }) {
                graph_1_iter.next().unwrap()
            } else {
                false
            }
        }))
        .unwrap();

    let graph_2_test_2 = graph_2
        .new_edgevec_from_iter(graph_2.iter_edges().map(|(pair, _, _)| {
            if matches!(pair, HedgePair::Paired { .. }) {
                graph_2_iter.next().unwrap()
            } else {
                false
            }
        }))
        .unwrap();

    assert_eq!(graph_1_test, graph_1_test_2);
    assert_eq!(graph_2_test, graph_2_test_2);
}

fn grouped_graph_fixture(name: &str, group_id: Option<usize>, master: Option<bool>) -> String {
    let group_id = group_id.map_or_else(String::new, |id| format!(" group_id={id}"));
    let master = master.map_or_else(String::new, |value| format!(" is_group_master={value}"));
    format!(
        r#"digraph {name} {{
            graph [projector="1"{group_id}{master}]
            node [num="1"]
            edge [pdg=1000 num="1" dir=none]
            ext [style=invis]
            ext -> A [sink="{{ufo_order:0}}"]
            A -> B [lmb_id=0 source="{{ufo_order:1}}" sink="{{ufo_order:0}}"]
            A -> B [source="{{ufo_order:2}}" sink="{{ufo_order:1}}"]
            B -> ext [source="{{ufo_order:2}}"]
        }}"#
    )
}

fn parse_grouped_graphs(specs: &[(&str, Option<usize>, Option<bool>)]) -> Vec<Graph> {
    let input = specs
        .iter()
        .map(|(name, group_id, master)| grouped_graph_fixture(name, *group_id, *master))
        .collect::<Vec<_>>()
        .join("\n");
    input
        .into_finalized_runtime_graph(&crate::utils::load_generic_model("scalars"))
        .unwrap()
}

#[test]
fn test_group_parsing_1() {
    test_initialise().unwrap();
    let mut graphs = parse_grouped_graphs(&[
        ("G1", Some(0), Some(true)),
        ("G2", Some(0), Some(false)),
        ("G3", Some(1), Some(true)),
        ("G4", Some(1), Some(false)),
    ]);

    let groups = complete_group_parsing(&mut graphs).unwrap();
    let expected_result = ti_vec![
        GraphGroup {
            master: 0,
            remaining: vec![1],
        },
        GraphGroup {
            master: 2,
            remaining: vec![3],
        },
    ];

    assert_eq!(groups, expected_result);
}

#[test]
fn test_group_parsing_2() {
    test_initialise().unwrap();
    let mut graphs = parse_grouped_graphs(&[
        ("G1", Some(0), Some(true)),
        ("G2", Some(0), None),
        ("G3", Some(1), None),
        ("G4", Some(1), None),
    ]);

    let groups = complete_group_parsing(&mut graphs).unwrap();
    let expected_result = ti_vec![
        GraphGroup {
            master: 0,
            remaining: vec![1],
        },
        GraphGroup {
            master: 2,
            remaining: vec![3],
        },
    ];

    assert_eq!(groups, expected_result);
}

#[test]
fn test_group_parsing_3() {
    test_initialise().unwrap();
    let mut graphs = parse_grouped_graphs(&[
        ("G1", Some(0), Some(true)),
        ("G2", Some(0), Some(false)),
        ("G3", None, None),
        ("G4", None, None),
    ]);

    let groups = complete_group_parsing(&mut graphs).unwrap();
    let expected_result = ti_vec![
        GraphGroup {
            master: 0,
            remaining: vec![1],
        },
        GraphGroup {
            master: 2,
            remaining: vec![],
        },
        GraphGroup {
            master: 3,
            remaining: vec![],
        },
    ];

    assert_eq!(groups, expected_result);
}

#[test]
fn test_group_parsing_4() {
    test_initialise().unwrap();
    let mut graphs = parse_grouped_graphs(&[
        ("G1", Some(0), Some(true)),
        ("G2", Some(0), Some(true)),
        ("G3", None, None),
        ("G4", None, None),
    ]);

    let groups = complete_group_parsing(&mut graphs);
    assert!(groups.is_err());
}

#[test]
fn test_group_parsing_5() {
    test_initialise().unwrap();
    let mut graphs = parse_grouped_graphs(&[
        ("G1", Some(1), Some(true)),
        ("G2", Some(1), Some(false)),
        ("G3", Some(2), Some(true)),
        ("G4", Some(2), Some(false)),
    ]);

    assert!(complete_group_parsing(&mut graphs).is_err());
}

#[test]
fn parse_triangle_lmb() {
    test_initialise().unwrap();

    for edge in 3..=5 {
        let g = parse_triangle_fixture(edge, None);
        assert_eq!(
            g.loop_momentum_basis.loop_edges[LoopIndex::from(0)],
            EdgeIndex::from(edge)
        );
    }
}

#[test]
fn edge_mass_attribute_drives_evaluated_edge_mass() {
    test_initialise().unwrap();

    let g = parse_triangle_fixture(3, Some(7));

    let model = crate::utils::load_generic_model("scalars");
    assert_eq!(
        g.underlying[EdgeIndex::from(3)]
            .mass_atom(&model)
            .to_string(),
        "7"
    );

    let evaluated_mass = g.underlying[EdgeIndex::from(3)]
        .mass_value::<f64>(&model, &g.param_builder)
        .unwrap();
    assert_eq!(evaluated_mass.re.0, 7.0);
    assert_eq!(evaluated_mass.im.0, 0.0);
}

#[test]
fn explicit_hedge_payload_round_trips_in_dot_export() {
    test_initialise().unwrap();
    let g: Graph = finalized_runtime_dot!(
        digraph payload_graph{
            graph [projector="1"]
            ext_in [style=invis]
            ext_out [style=invis]
            A [num=1 dod=0]
            ext_in -> A [name=e_in num=1 dod=-2 sink="{ufo_order:0,dod:-2}"]
            A -> ext_out [name=e_out num=1 dod=-2 source="{ufo_order:1,dod:-2}"]
        }
    )
    .unwrap();

    let serialized = g.dot_serialize(&DotExportSettings::default());
    assert!(serialized.contains("ufo_order"));
    assert!(serialized.contains("dod"));
    assert!(serialized.contains("source=") || serialized.contains("sink="));
}

fn assert_native_runtime_parity(diagram: &feynkit_graph::FeynmanDiagram) {
    use crate::uv::UltravioletGraph;
    use linnet::half_edge::{
        EdgeAccessors,
        involution::{Flow, Hedge},
    };
    use symbolica::atom::AtomCore;

    let runtime = Graph::from_feynkit(diagram, None, true).unwrap();
    let source = diagram.underlying();
    assert_eq!(runtime.underlying.n_nodes(), source.n_nodes());
    assert_eq!(runtime.underlying.n_edges(), source.n_edges());
    assert_eq!(runtime.underlying.n_hedges(), source.n_hedges());
    for (pair, edge, data) in source.iter_edges() {
        let (runtime_data, runtime_pair) = &runtime.underlying[&edge];
        assert_eq!(*runtime_pair, pair);
        assert_eq!(runtime.underlying.orientation(edge), data.orientation);
        assert_eq!(runtime_data.num.value, data.data.numerator);
        assert_eq!(runtime_data.is_dummy, data.data.is_dummy);
    }
    for index in 0..source.n_hedges() {
        let hedge = Hedge(index);
        assert_eq!(runtime.underlying.node_id(hedge), source.node_id(hedge));
        assert_eq!(runtime.underlying.flow(hedge), source.flow(hedge));
        assert_eq!(runtime.underlying[&hedge], source[&hedge]);
    }
    assert_eq!(runtime.overall_factor, *diagram.overall_factor());
    assert_eq!(runtime.global_prefactor.num, *diagram.numerator_prefactor());
    assert_eq!(runtime.global_prefactor.projector, *diagram.projector());
    let transported: feynkit_graph::routing::MomentumBasis = (&runtime.loop_momentum_basis).into();
    assert_eq!(
        transported,
        diagram.loop_momentum_basis().to_routing(source)
    );
    let full = source.full_filter();
    let empty = source.empty_subgraph();
    assert_eq!(
        runtime
            .numerator(&full, &empty)
            .get_single_atom()
            .unwrap()
            .expand(),
        diagram.numerator().expand()
    );
    for (_, edge, _) in source.iter_edges() {
        let mut selected: linnet::half_edge::subgraph::SuBitGraph = source.empty_subgraph();
        use linnet::half_edge::subgraph::ModifySubSet;
        selected.add(source[&edge].1);
        assert_eq!(
            runtime
                .numerator(&selected, &empty)
                .get_single_atom()
                .unwrap(),
            diagram.numerator_of(&selected, &empty)
        );
        assert_eq!(
            runtime.denominator(&selected, diagram.model(), |_| -2),
            diagram
                .denominator_of(
                    &selected,
                    &std::collections::BTreeMap::from([(EdgeId(edge.0), -2)])
                )
                .unwrap()
        );
    }
    let half_edges = |filter: &linnet::half_edge::subgraph::SuBitGraph| {
        filter
            .included_iter()
            .map(|hedge| feynkit_graph::DiagramHalfEdge {
                edge: EdgeId(source[&hedge].0),
                endpoint: if source.flow(hedge) == Flow::Source {
                    feynkit_graph::DiagramEndpoint::Source
                } else {
                    feynkit_graph::DiagramEndpoint::Target
                },
            })
            .collect::<std::collections::BTreeSet<_>>()
    };
    assert_eq!(runtime.finalized_cuts.len(), diagram.cuts().len());
    for (runtime, source) in runtime.finalized_cuts.iter().zip(diagram.cuts()) {
        assert_eq!(
            half_edges(&runtime.left),
            source.left.half_edges.iter().copied().collect()
        );
        assert_eq!(
            half_edges(&runtime.right),
            source.right.half_edges.iter().copied().collect()
        );
        assert_eq!(
            half_edges(&runtime.cut.left),
            source.cut.iter().copied().collect()
        );
    }
    assert_eq!(
        runtime.finalized_topology_threshold_candidates.len(),
        diagram.topology_threshold_candidates().len()
    );
    for (runtime, source) in runtime
        .finalized_topology_threshold_candidates
        .iter()
        .zip(diagram.topology_threshold_candidates())
    {
        assert_eq!(
            half_edges(&runtime.left),
            source.left.iter().copied().collect()
        );
        assert_eq!(
            half_edges(&runtime.right),
            source.right.iter().copied().collect()
        );
        assert_eq!(
            half_edges(&runtime.cut.left),
            source.cut.iter().copied().collect()
        );
    }
}

#[test]
fn native_amplitude_runtime_enrichment_preserves_every_graph_coordinate() {
    use feynkit_generator::{GenerationOptions, Process};
    test_initialise().unwrap();
    let model = std::sync::Arc::new(crate::utils::load_generic_model("scalars"));
    let generated = Process::new(["scalar_1"], ["scalar_1"])
        .generate_diagrams(
            model,
            &GenerationOptions::default()
                .with_loop_count(1, 1)
                .unwrap()
                .threads(1)
                .max_vertices(2),
        )
        .unwrap();
    assert!(!generated.diagrams.is_empty());
    for diagram in generated.diagrams {
        assert_native_runtime_parity(&diagram);
    }
}

#[test]
fn native_cross_section_runtime_enrichment_preserves_cuts_and_positive_external_flow() {
    use feynkit_generator::{GenerationOptions, Process};
    test_initialise().unwrap();
    let model = std::sync::Arc::new(crate::utils::load_generic_model("scalars"));
    let generated = Process::new(["scalar_1"], ["scalar_1", "scalar_1"])
        .generate_cross_section(
            model,
            &GenerationOptions::default()
                .with_loop_count(1, 1)
                .unwrap()
                .threads(1)
                .max_vertices(2),
        )
        .unwrap();
    assert!(!generated.diagrams.is_empty());
    for diagram in generated.diagrams {
        assert!(!diagram.cuts().is_empty());
        for (slot, edge) in diagram
            .loop_momentum_basis()
            .external_edges
            .iter()
            .enumerate()
        {
            assert_eq!(
                diagram.loop_momentum_basis().edge_signatures[edge]
                    .external
                    .get(slot),
                Some(feynkit_kinematics::SignOrZero::Plus)
            );
        }
        assert_native_runtime_parity(&diagram);
    }
}

#[test]
fn native_compton_runtime_enrichment_preserves_sewn_fermions() {
    use feynkit_generator::{GenerationFilter, GenerationOptions, Process};
    test_initialise().unwrap();
    let model = std::sync::Arc::new(
        feynkit_model::Model::from_json(include_str!(
            "../../../../feynkit-model/tests/fixtures/sm.json"
        ))
        .unwrap(),
    );

    for pdg in [11_i64, -11] {
        let generated = Process::new([pdg, 22], [pdg, 22])
            .generate_cross_section(
                model.clone(),
                &GenerationOptions::default()
                    .with_loop_count(1, 1)
                    .unwrap()
                    .threads(1)
                    .max_vertices(4)
                    .with_graph_filter(GenerationFilter::VertexAllow(vec!["V_98".into()])),
            )
            .unwrap();
        assert_eq!(generated.diagrams.len(), 4);
        for diagram in generated.diagrams {
            assert_native_runtime_parity(&diagram);
        }
    }
}

#[test]
fn native_runtime_enrichment_preserves_explicit_loop_order() {
    use feynkit_graph::{DiagramVertex, ExternalState, FeynmanDiagram};
    use symbolica::atom::Atom;
    test_initialise().unwrap();
    let model = std::sync::Arc::new(crate::utils::load_generic_model("scalars"));
    let particle = model.particle_id("scalar_1").unwrap();
    let mut builder = FeynmanDiagram::builder(model, "ordered_self_loops");
    let node = builder.add_vertex(DiagramVertex {
        name: "vertex".into(),
        interaction: None,
        numerator: Atom::one(),
    });
    let mut incoming = DiagramEdge::new(particle, false);
    incoming.external = Some(feynkit_graph::ExternalLeg {
        name: "incoming".into(),
        index: 0,
        state: ExternalState::Incoming,
        connection: 0,
    });
    let mut outgoing = DiagramEdge::new(particle, false);
    outgoing.external = Some(feynkit_graph::ExternalLeg {
        name: "outgoing".into(),
        index: 1,
        state: ExternalState::Outgoing,
        connection: 1,
    });
    builder.add_edge(None, Some(node), incoming).unwrap();
    builder.add_edge(Some(node), None, outgoing).unwrap();
    builder
        .add_edge(node, node, DiagramEdge::new(particle, false))
        .unwrap();
    builder
        .add_edge(node, node, DiagramEdge::new(particle, false))
        .unwrap();
    let diagram = builder.build().unwrap();
    let mut selected = diagram.loop_momentum_basis().loop_edges.clone();
    assert_eq!(selected.len(), 2);
    selected.reverse();
    let diagram = diagram.with_loop_momentum_edges(&selected).unwrap();
    assert_eq!(diagram.loop_momentum_basis().loop_edges, selected);
    assert_native_runtime_parity(&diagram);
}

#[test]
fn native_subgraph_routing_matches_runtime_for_empty_and_disconnected_selections() {
    use crate::graph::LMBext;
    use feynkit_graph::{DiagramVertex, FeynmanDiagram};
    use linnet::half_edge::subgraph::{ModifySubSet, SuBitGraph};
    use symbolica::atom::Atom;
    test_initialise().unwrap();
    let model = std::sync::Arc::new(crate::utils::load_generic_model("scalars"));
    let particle = model.particle_id("scalar_1").unwrap();
    let mut builder = FeynmanDiagram::builder(model, "disconnected_selections");
    let nodes = (0..4)
        .map(|index| {
            builder.add_vertex(DiagramVertex {
                name: format!("v{index}"),
                interaction: None,
                numerator: Atom::one(),
            })
        })
        .collect::<Vec<_>>();
    for (index, state, source, target) in [
        (
            0,
            feynkit_graph::ExternalState::Incoming,
            None,
            Some(nodes[0]),
        ),
        (
            1,
            feynkit_graph::ExternalState::Outgoing,
            Some(nodes[3]),
            None,
        ),
    ] {
        let mut edge = DiagramEdge::new(particle, false);
        edge.external = Some(feynkit_graph::ExternalLeg {
            name: format!("external{index}"),
            index,
            state,
            connection: index,
        });
        builder.add_edge(source, target, edge).unwrap();
    }
    for (source, target) in [(0, 1), (0, 1), (1, 2), (2, 3), (2, 3)] {
        builder
            .add_edge(
                nodes[source],
                nodes[target],
                DiagramEdge::new(particle, false),
            )
            .unwrap();
    }
    let diagram = builder.build().unwrap();
    let runtime = Graph::from_feynkit(&diagram, None, true).unwrap();
    let empty: SuBitGraph = diagram.underlying().empty_subgraph();
    let mut disconnected = empty.clone();
    for edge in [2, 3, 5, 6] {
        disconnected.add(diagram.underlying()[&EdgeIndex(edge)].1);
    }
    let mut partial = empty.clone();
    for edge in [2, 3] {
        partial.add(diagram.underlying()[&EdgeIndex(edge)].1);
    }
    for selected in [empty, partial, disconnected] {
        let native = diagram
            .momentum_basis_of(&selected)
            .unwrap()
            .to_routing(diagram.underlying());
        let reference: feynkit_graph::routing::MomentumBasis = (&runtime.lmb_of(&selected)).into();
        assert_eq!(
            native.tree.included_iter().collect::<Vec<_>>(),
            reference.tree.included_iter().collect::<Vec<_>>()
        );
        assert_eq!(native.loop_edges, reference.loop_edges);
        assert_eq!(native.ext_edges, reference.ext_edges);
        assert_eq!(native.edge_signatures, reference.edge_signatures);
    }
}

#[path = "finalized_reference.rs"]
mod finalized_reference;
