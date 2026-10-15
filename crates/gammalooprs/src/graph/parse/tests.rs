use feynkit_graph::{DiagramEdge, EdgeEndpoints, EdgeId, VertexId};
use linnet::half_edge::NodeIndex;
use linnet::half_edge::{
    involution::{EdgeIndex, HedgePair},
    subgraph::SubSetLike,
};
use spenso::{
    network::tags::SPENSO_TAG,
    structure::{
        abstract_index::AIND_SYMBOLS,
        representation::{Minkowski, RepName},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder},
    symbol,
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
    utils::GS,
    uv::{UltravioletGraph, uv_graph::UVE},
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
            source
                .left
                .half_edges
                .iter()
                .copied()
                .collect::<std::collections::BTreeSet<_>>()
        );
        assert_eq!(
            half_edges(&runtime.right),
            source
                .right
                .half_edges
                .iter()
                .copied()
                .collect::<std::collections::BTreeSet<_>>()
        );
        assert_eq!(
            half_edges(&runtime.cut.left),
            source
                .cut
                .iter()
                .copied()
                .collect::<std::collections::BTreeSet<_>>()
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
            source
                .left
                .iter()
                .copied()
                .collect::<std::collections::BTreeSet<_>>()
        );
        assert_eq!(
            half_edges(&runtime.right),
            source
                .right
                .iter()
                .copied()
                .collect::<std::collections::BTreeSet<_>>()
        );
        assert_eq!(
            half_edges(&runtime.cut.left),
            source
                .cut
                .iter()
                .copied()
                .collect::<std::collections::BTreeSet<_>>()
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

#[test]
fn native_self_energy_tensor_ports_match_runtime_component_interfaces() {
    use crate::numerator::{aind::Aind, symbolica_ext::NumeratorAtomExt};
    use idenso::{dirac::AGS, representations::Bispinor};
    use spenso::{
        network::library::symbolic::ETS,
        structure::{TensorStructure, slot::IsAbstractSlot},
    };
    use symbolica::function;

    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("sm");
    let source = r#"digraph port_self_energy {
        projector=1
        ext[style=invis]
        ext -> a[id=0 particle=t]
        a -> b[id=1 particle=t lmb_id=0]
        a -> b[id=2 particle=g]
        b -> ext[id=3 particle=t]
    }"#;
    let graph = Graph::import_dot(source, &model).unwrap().remove(0);
    let numerator = graph
        .numerator(&graph.full_filter(), &graph.empty_subgraph())
        .get_single_atom()
        .unwrap();
    let mut source_ports = numerator.parse_into_net().unwrap().graph.dangling_indices();
    source_ports.sort();
    let external = graph.external_filter::<linnet::half_edge::subgraph::SuBitGraph>();
    let mut runtime_ports = Vec::new();
    for hedge in external.included_iter() {
        let data = &graph.underlying[hedge];
        runtime_ports.extend(
            data.num_indices
                .spin_indices
                .vertex_indices
                .external_structure_iter(),
        );
        runtime_ports.extend(
            data.num_indices
                .color_indices
                .vertex_indices
                .external_structure_iter(),
        );
        let spin = data
            .edge_spin_slots()
            .next()
            .expect("a top boundary has one spin port");
        assert_eq!(spin.aind, Aind::Hedge(hedge.0 as u16, 1));
        assert!(
            numerator.contains(spin.to_atom()),
            "the component port must occur in its source factor"
        );
        let dummy = Bispinor {}
            .new_rep(4)
            .to_symbolic([Atom::var(symbol!("port_interface::dummy"))]);
        let port = spin.to_atom();
        let projector = (function!(ETS.metric, &port, &dummy)
            + function!(AGS.gamma0, &port, &dummy))
            / Atom::num(2);
        let projected = projector
            * numerator
                .replace(port.to_pattern())
                .with(dummy.to_pattern());
        let mut projected_ports = projected.parse_into_net().unwrap().graph.dangling_indices();
        projected_ports.sort();
        assert_eq!(
            projected_ports, source_ports,
            "ordered projector insertion cannot create a phantom port"
        );
    }
    runtime_ports.sort();
    assert_eq!(
        runtime_ports, source_ports,
        "spin and oriented color interfaces must both agree"
    );

    // A literal artifact owns its supplied labels, including the zero-based
    // labels of hand-authored numerators. Reindex complete local factors only.
    let reindex = |atom: &Atom| {
        (0..graph.underlying.n_hedges()).fold(atom.clone(), |atom, hedge| {
            atom.replace(Atom::from(Aind::Hedge(hedge as u16, 1)).to_pattern())
                .with(Atom::from(Aind::Hedge(hedge as u16, 0)).to_pattern())
        })
    };
    let mut literal = graph.clone();
    literal.underlying = graph.underlying.map_data_ref(
        |_, _, vertex| {
            let mut vertex = vertex.clone();
            vertex.num.value = reindex(&vertex.num);
            vertex
        },
        |_, _, _, edge| {
            edge.map(|edge| {
                let mut edge = edge.clone();
                edge.num.value = reindex(&edge.num);
                edge
            })
        },
        |_, data| data.clone(),
    );
    let restored = Graph::import_dot(
        &literal.dot_serialize(&DotExportSettings::default()),
        &model,
    )
    .unwrap()
    .remove(0);
    for hedge in external.included_iter() {
        let port = restored.underlying[hedge].edge_spin_slots().next().unwrap();
        assert_eq!(port.aind, Aind::Hedge(hedge.0 as u16, 0));
        assert!(
            restored
                .numerator(&restored.full_filter(), &restored.empty_subgraph())
                .get_single_atom()
                .unwrap()
                .contains(port.to_atom())
        );
    }

    let hedge = external.included_iter().next().unwrap();
    let [first, second] = [0, 1].map(|local| {
        Bispinor {}
            .new_rep(4)
            .to_symbolic([Atom::from(Aind::Hedge(hedge.0 as u16, local))])
    });
    literal.global_prefactor.num *= function!(ETS.metric, first, second);
    let error = Graph::import_dot(
        &literal.dot_serialize(&DotExportSettings::default()),
        &model,
    )
    .err()
    .expect("conflicting literal labels must not select a guessed port");
    assert!(
        error
            .to_string()
            .contains("ambiguous literal tensor interface"),
        "{error:?}"
    );

    literal.global_prefactor.num = Atom::one();
    literal.underlying = literal.underlying.map_data_ref(
        |_, _, vertex| {
            let mut vertex = vertex.clone();
            vertex.num.value = Atom::one();
            vertex
        },
        |_, _, _, edge| {
            edge.map(|edge| {
                let mut edge = edge.clone();
                edge.num.value = Atom::one();
                edge
            })
        },
        |_, data| data.clone(),
    );
    let scalar = Graph::import_dot(
        &literal.dot_serialize(&DotExportSettings::default()),
        &model,
    )
    .unwrap()
    .remove(0);
    for hedge in external.included_iter() {
        assert!(
            scalar.underlying[hedge].edge_spin_slots().next().is_none(),
            "literal scalar fragments cannot manufacture a component spin port"
        );
    }
}

#[test]
fn compact_model_import_preserves_metadata_in_the_finalized_frame() {
    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("scalars");
    let mut parsed = linnet::parser::DotGraph::from_string(
        r#"digraph mapped_bubble {
        final_state="scalar_1,scalar_1"
        group_id=0
        is_group_master=true
        params="source_weight"
        ext [style=invis]
        ext -> a [id=0 particle=scalar_1 is_cut=7]
        a -> b [id=1 particle=scalar_1 lmb_id=0]
        b -> ext [id=2 particle=scalar_1 is_cut=7]
        a -> b [id=3 particle=scalar_1]
    }"#,
    )
    .unwrap();
    parsed.global_data.statements.insert(
        "threshold_counterterms".into(),
        r#"
schema_version = 1
[[cuts]]
edges = [1, 3]
[[cuts.thresholds]]
edges = [1, 3]
[[cuts.thresholds.counterterms]]
subspace = [3]
parent_lmb = [3]
"#
        .into(),
    );
    let imported = Graph::from_compact_model_dot(parsed, &model).unwrap();
    assert_eq!(imported.group_id, Some(crate::graph::GroupId(0)));
    assert!(imported.is_group_master);
    assert_eq!(
        imported.param_builder.pairs.additional_params.params.len(),
        1
    );
    let cut = &imported.threshold_counterterms.cuts[0];
    assert_eq!(cut.edges, [EdgeIndex(1), EdgeIndex(2)]);
    assert_eq!(cut.thresholds[0].edges, [EdgeIndex(1), EdgeIndex(2)]);
    let variant = &cut.thresholds[0].counterterms[0];
    assert_eq!(variant.subspace.as_ref().unwrap(), &[EdgeIndex(2)]);
    assert_eq!(variant.parent_lmb.as_ref().unwrap(), &[EdgeIndex(2)]);
    assert!(!imported.finalized_cuts.is_empty());
}

#[test]
fn compact_qcd_import_instantiates_rules_and_retains_explicit_fragments() {
    use feynkit_generator::FeynmanDiagramRulesExt;
    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("sm");
    let source = r#"digraph physical_top_self_energy {
        ext [style=invis]
        ext -> a [id=0 particle=t]
        a -> b [id=1 particle=t lmb_id=0]
        a -> b [id=2 particle=g]
        b -> ext [id=3 particle=t]
    }"#;
    let topology =
        feynkit_graph::FeynmanDiagram::from_dot(std::sync::Arc::new(model.clone()), source)
            .unwrap();
    assert_eq!(topology.numerator(), &Atom::one());
    let reference = topology.apply_feynman_rules().unwrap();
    let imported = Graph::import_dot(source, &model).unwrap().remove(0);
    assert_ne!(reference.numerator(), &Atom::one());
    assert_eq!(
        imported
            .numerator(&imported.full_filter(), &imported.empty_subgraph())
            .get_single_atom()
            .unwrap(),
        *reference.numerator()
    );
    assert_eq!(
        imported.global_prefactor.projector,
        reference.external_state_projector().unwrap()
    );
    assert!(!imported.global_prefactor.projector.is_one());
    let supplied = source.replace(
        "ext [style=invis]",
        "projector=1; a[num=7]; ext [style=invis]",
    );
    let explicit = Graph::import_dot(&supplied, &model).unwrap().remove(0);
    assert!(explicit.global_prefactor.projector.is_one());
    assert_eq!(
        explicit.underlying[linnet::half_edge::NodeIndex(0)]
            .num
            .value,
        Atom::num(7)
    );
    assert!(!explicit.underlying[EdgeIndex(1)].num.value.is_one());
}

#[test]
fn compact_gamma_loop_indices_bind_to_native_rule_ports() {
    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("sm");
    let source =
        include_str!("../../../../../tests/resources/graphs/local_ir_gluon_self_energy.dot");
    let qualified = Graph::import_dot(source, &model).unwrap().remove(0);
    let unqualified = Graph::import_dot(&source.replace("gammalooprs::", ""), &model)
        .unwrap()
        .remove(0);
    // This explicitly supplied gluon projector closes the independently
    // instantiated color and Lorentz ports; disconnected copies of either
    // interface would fail exact network closure.
    qualified.validate_full_numerator_tensor_network().unwrap();
    unqualified
        .validate_full_numerator_tensor_network()
        .unwrap();
    assert_eq!(qualified.global_prefactor, unqualified.global_prefactor);
    let projector = &qualified.global_prefactor.projector;
    for hedge in qualified
        .external_filter::<linnet::half_edge::subgraph::SuBitGraph>()
        .included_iter()
    {
        let native = Atom::from(crate::numerator::aind::Aind::Hedge(hedge.0 as u16, 1));
        let source = Atom::from(crate::numerator::aind::Aind::Hedge(hedge.0 as u16, 0));
        assert!(projector.contains(&native));
        assert!(!projector.contains(&source));
    }
}

#[test]
fn compact_sewn_projector_uses_the_source_half_edge_chart() {
    use crate::numerator::aind::Aind;
    use idenso::representations::Bispinor;
    use linnet::half_edge::involution::Flow;
    use symbolica::function;
    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("sm");
    let mut source = linnet::parser::DotGraph::from_string(
        r#"digraph sewn_top_self_energy {
        final_state="t,g"
        ext[style=invis]
        ext -> a[id=0 particle=t is_cut=7]
        a -> b[id=1 particle=t lmb_id=0]
        b -> ext[id=2 particle=t is_cut=7]
        a -> b[id=3 particle=g]
    }"#,
    )
    .unwrap();
    let automatic = Graph::from_compact_model_dot(source.clone(), &model).unwrap();
    let mut projector = Atom::one();
    for (pair, edge, _) in source.iter_edges() {
        if let HedgePair::Unpaired { hedge, flow } = pair {
            let slot = Bispinor {}
                .new_rep(4)
                .to_symbolic([Atom::from(Aind::Hedge(hedge.0 as u16, 0))]);
            projector *= function!(
                if flow == Flow::Sink { GS.u } else { GS.ubar },
                edge.0,
                slot
            );
        }
    }
    source
        .global_data
        .statements
        .insert("projector".into(), projector.to_canonical_string());
    let explicit = Graph::from_compact_model_dot(source, &model).unwrap();
    assert_eq!(
        explicit.global_prefactor.projector, automatic.global_prefactor.projector,
        "supplied wavefunctions must follow actual source incidence through sewing"
    );
    assert_eq!(explicit.loop_momentum_basis, automatic.loop_momentum_basis);
    assert_eq!(
        explicit
            .numerator(&explicit.full_filter(), &explicit.empty_subgraph())
            .get_single_atom()
            .unwrap(),
        automatic
            .numerator(&automatic.full_filter(), &automatic.empty_subgraph())
            .get_single_atom()
            .unwrap()
    );
}

#[test]
fn full_numerator_requires_one_q3_structure_argument() {
    test_initialise().unwrap();
    let mut graph: Graph = crate::graph::Graph::import_dot(
        stringify!(
        digraph q3_validation {
            edge [pdg=1000]
            ext [style=invis]
            ext -> v4
            ext -> v5
            v6 -> ext
            v5 -> v4 [lmb_id=0]
            v6 -> v5
            v4 -> v6
        }),
        &crate::utils::load_generic_model("scalars"),
    )
    .map(|mut graphs| {
        assert_eq!(graphs.len(), 1);
        graphs.remove(0)
    })
    .unwrap();

    graph.overall_factor = GS.emr_vec.call_args([Atom::num(3)]);
    let error = graph.validate_full_numerator_tensor_network().unwrap_err();
    assert!(
        error
            .to_string()
            .contains("Q3 in graph 'q3_validation' must have exactly two arguments"),
        "unexpected Q3 arity error: {error:?}"
    );

    graph.overall_factor = GS
        .emr_vec
        .call_args([Atom::var(symbol!("invalid_q3_edge")), GS.cind(1)]);
    let error = graph.validate_full_numerator_tensor_network().unwrap_err();
    assert!(
        error
            .to_string()
            .contains("found non-canonical edge argument invalid_q3_edge"),
        "unexpected Q3 edge error: {error:?}"
    );

    graph.overall_factor = GS.emr_vec(EdgeIndex(3), Atom::var(symbol!("invalid_q3_structure")));
    let error = graph.validate_full_numerator_tensor_network().unwrap_err();
    assert!(
        error
            .to_string()
            .contains("found non-canonical second argument invalid_q3_structure"),
        "unexpected Q3 structure error: {error:?}"
    );

    graph.overall_factor = GS.emr_vec(EdgeIndex(3), Minkowski {}.new_rep(3).to_symbolic([]));
    let error = graph.validate_full_numerator_tensor_network().unwrap_err();
    assert!(
        error
            .to_string()
            .contains("found non-canonical second argument"),
        "unexpected Q3 dimension error: {error:?}"
    );

    let out_of_range_edge = EdgeIndex(graph.underlying.n_edges());
    graph.overall_factor = GS.emr_vec(out_of_range_edge, Minkowski {}.new_rep(4).to_symbolic([]));
    let error = graph.validate_full_numerator_tensor_network().unwrap_err();
    assert!(
        error.to_string().contains(&format!(
            "found edge ID {} outside the valid range 0..{}",
            out_of_range_edge.0,
            graph.underlying.n_edges(),
        )),
        "unexpected out-of-range Q3 edge error: {error:?}"
    );

    graph.overall_factor = GS.emr_mom(out_of_range_edge, GS.cind(1));
    let error = graph.validate_full_numerator_tensor_network().unwrap_err();
    assert!(
        error.to_string().contains(&format!(
            "found edge ID {} outside the valid range 0..{}",
            out_of_range_edge.0,
            graph.underlying.n_edges(),
        )),
        "unexpected out-of-range Q edge error: {error:?}"
    );

    let minkowski = Minkowski {}.new_rep(4).to_symbolic([]);
    let represented = GS.emr_vec(EdgeIndex(3), minkowski);
    graph.overall_factor = SPENSO_TAG.dot.call_args([represented.clone(), represented]);
    graph.validate_full_numerator_tensor_network().unwrap();

    let indexed_minkowski = Minkowski {}
        .new_rep(4)
        .to_symbolic([Atom::var(symbol!("q3_validator_mu"))]);
    let indexed = GS.emr_vec(EdgeIndex(3), indexed_minkowski);
    graph.overall_factor = &indexed * &indexed;
    graph.validate_full_numerator_tensor_network().unwrap();

    graph.overall_factor = GS.emr_vec(EdgeIndex(3), GS.cind(1));
    assert_eq!(
        graph.overall_factor,
        GS.emr_mom(EdgeIndex(3), GS.cind(1)),
        "a transient tensor-network Q3 component must normalize immediately",
    );

    let abstract_slot = Minkowski {}.new_rep(4).to_symbolic([Atom::num(1)]);
    graph.overall_factor = GS.emr_vec(EdgeIndex(3), abstract_slot);
    let error = graph.validate_full_numerator_tensor_network().unwrap_err();
    assert!(
        error
            .to_string()
            .contains("Full numerator still has dangling tensor indices"),
        "a numeric Minkowski index must be validated as an abstract index: {error:?}",
    );
    assert!(
        graph.overall_factor.contains_symbol(GS.emr_vec)
            && !graph.overall_factor.contains_symbol(GS.emr_mom),
        "a numeric Minkowski index must remain an abstract Q3 label, not normalize to a Q component",
    );

    for malformed_component in [
        AIND_SYMBOLS.cind.call_args([Atom::num(-1)]),
        AIND_SYMBOLS.cind.call_args([Atom::num(4)]),
        AIND_SYMBOLS.cind.call_args([Atom::num(1), Atom::num(2)]),
    ] {
        graph.overall_factor = GS.emr_vec(EdgeIndex(3), malformed_component);
        let error = graph.validate_full_numerator_tensor_network().unwrap_err();
        assert!(
            error
                .to_string()
                .contains("found non-canonical second argument"),
            "unexpected malformed Q3 component error: {error:?}",
        );
    }

    let wrapper = symbol!("q3_validation_wrapper");
    let scalar = Atom::var(symbol!("q3_validation_scalar"));
    let malformed_q3 = GS.emr_vec.call_args([Atom::num(3)]);
    for nested in [
        wrapper.call_args([malformed_q3.clone()]),
        wrapper.call_args([malformed_q3.pow(2)]),
        wrapper.call_args([&malformed_q3 + &scalar]),
        wrapper.call_args([&malformed_q3 * &scalar]),
    ] {
        graph.overall_factor = nested;
        let error = graph.validate_spatial_momentum_convention().unwrap_err();
        assert!(
            error.to_string().contains("Q3 in graph 'q3_validation'"),
            "nested factorized payload must retain Q3 validation: {error:?}",
        );
    }

    let malformed_q = FunctionBuilder::new(GS.emr_mom).finish();
    graph.overall_factor = wrapper.call_args([&malformed_q3 + &malformed_q.pow(2)]);
    let error = graph.validate_spatial_momentum_convention().unwrap_err();
    assert!(
        error.to_string().contains("Spatial momentum Q3"),
        "the Q3 pass must retain priority over malformed Q: {error:?}",
    );

    let nested_q = wrapper.call_args([scalar.pow(malformed_q)]);
    graph.overall_factor = nested_q.clone();
    let error = graph.validate_spatial_momentum_convention().unwrap_err();
    assert!(
        error.to_string().contains("found no edge argument"),
        "a malformed Q in a power exponent must remain visible: {error:?}",
    );

    graph.overall_factor = Atom::one();
    let original_global = graph.global_prefactor.clone();
    graph.global_prefactor.num = nested_q.clone();
    assert!(graph.validate_spatial_momentum_convention().is_err());
    graph.global_prefactor = original_global.clone();
    graph.global_prefactor.projector = nested_q;
    assert!(graph.validate_spatial_momentum_convention().is_err());
    graph.global_prefactor = original_global;

    // Component normalization happens before either traversal sees the atom.
    let component = GS.emr_vec(out_of_range_edge, GS.cind(1));
    assert_eq!(component, GS.emr_mom(out_of_range_edge, GS.cind(1)));
    graph.overall_factor = wrapper.call_args([component]);
    let error = graph.validate_spatial_momentum_convention().unwrap_err();
    assert!(
        error.to_string().contains("Momentum Q in graph"),
        "a normalized Q3 component must be checked as Q: {error:?}",
    );

    // Validate the normalized product, including cancellation and annihilation.
    graph.overall_factor = &malformed_q3 - &malformed_q3;
    graph.validate_spatial_momentum_convention().unwrap();
    graph.overall_factor = malformed_q3;
    graph.global_prefactor.num = Atom::Zero;
    graph.validate_spatial_momentum_convention().unwrap();
}

#[test]
fn graph_set_import_enforces_the_q3_convention_without_tensor_closure() {
    test_initialise().unwrap();
    for numerator in [
        r#"num = "gammalooprs::Q3(0)""#,
        r#"edge [num = "gammalooprs::Q3(0)"]"#,
        r#"node [num = "gammalooprs::Q3(0)"]"#,
    ] {
        let malformed = format!(
            r#"
            digraph q3_graph_set_validation {{
                {numerator}
                edge [pdg=1000]
                ext [style=invis]
                ext -> v4
                ext -> v5
                v6 -> ext
                v5 -> v4 [lmb_id=0]
                v6 -> v5
                v4 -> v6
            }}
            "#
        );

        let error = crate::graph::Graph::import_dot(
            &malformed,
            &crate::utils::load_generic_model("scalars"),
        )
        .err()
        .expect("multi-graph imports must reject the obsolete one-argument Q3 form");
        assert!(
            format!("{error:#}")
                .contains("Q3 in graph 'q3_graph_set_validation' must have exactly two arguments"),
            "unexpected graph-set Q3 validation error: {error:?}"
        );
    }
}

#[test]
fn hedge_payload_does_not_require_dod() {
    test_initialise().unwrap();
    let graph: Graph = crate::finalized_runtime_dot!(digraph missing_hedge_dod{
    graph [projector=1]
                   ext_in [style=invis]
                   ext_out [style=invis]
                   A [num=1 dod=0]
                   ext_in -> A [name=e_in num=1 dod=-2 sink="{ufo_order:0}"]
                   A -> ext_out [name=e_out num=1 dod=-2 source="{ufo_order:1}"]
               })
    .unwrap();

    let mut orders = graph
        .underlying
        .iter_hedges()
        .map(|(_, hedge)| {
            assert!(!hedge.ufo_order.autogenerated);
            hedge.ufo_order.value
        })
        .collect::<Vec<_>>();
    orders.sort();
    assert_eq!(orders, vec![0, 1]);
    assert_eq!(graph.underlying[NodeIndex(0)].dod.value, 0);
    assert!(!graph.underlying[NodeIndex(0)].dod.autogenerated);
    for edge in [EdgeIndex(0), EdgeIndex(1)] {
        assert_eq!(graph.underlying[edge].dod.value, -2);
        assert!(!graph.underlying[edge].dod.autogenerated);
    }
}

#[test]
fn initial_state_cut_permutation_preserves_specified_loop_edges() {
    test_initialise().unwrap();
    let graph: Graph = crate::finalized_runtime_dot!(digraph cut_lmb_permutation {
 graph [projector=1]
        edge [num=1 mass=2]
        node [num=1]
        a -> b [id=0 name=first_loop lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:0}"]
        a -> b [id=1 name=second_loop lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:1}"]
        a -> b [id=2 source="{ufo_order:2}" sink="{ufo_order:2}"]
        b -> a [id=3 is_cut=0 source="{ufo_order:3}" sink="{ufo_order:3}"]
        b -> a [id=4 is_cut=1 source="{ufo_order:4}" sink="{ufo_order:4}"]
    }, &crate::utils::load_generic_model("sm"))
    .unwrap();

    assert_eq!(
        graph.loop_momentum_basis.ext_edges.raw,
        [EdgeIndex(0), EdgeIndex(1)]
    );
    for (loop_id, name) in ["first_loop", "second_loop"].into_iter().enumerate() {
        assert_eq!(
            graph.loop_momentum_basis.loop_edges[LoopIndex(loop_id)],
            graph.edge_name_to_index(name).unwrap()
        );
    }
}

#[test]
fn missing_vertex_dod_is_autogenerated() {
    test_initialise().unwrap();
    let graph: Graph = crate::finalized_runtime_dot!(digraph missing_vertex_dod{
    graph [projector=1]
                   ext_in [style=invis]
                   ext_out [style=invis]
                   A [num=1]
                   ext_in -> A [name=e_in num=1 dod=-2 sink="{ufo_order:0,dod:-2}"]
                   A -> ext_out [name=e_out num=1 dod=-2 source="{ufo_order:1,dod:-2}"]
               })
    .unwrap();

    let dod = &graph.underlying[NodeIndex(0)].dod;
    assert_eq!(dod.value, 0);
    assert!(dod.autogenerated);
}

#[test]
fn open_fundamental_color_requires_an_explicit_projector() {
    test_initialise().unwrap();
    let mut graph: Graph = crate::graph::Graph::import_dot(stringify!(
        digraph open_fundamental_color {
            num = "1/3*spenso::g(spenso::dind(spenso::cof(3,gammalooprs::hedge(0))),spenso::cof(3,gammalooprs::hedge(1)))"
            ext [style=invis]

            ext -> v1:0 [id=0 particle=d]
            v2:1 -> ext [id=1 particle=d]
            v1:2 -> v2:3 [id=2 particle=d lmb_id=0]
            v1:4 -> v2:5 [id=3 particle=g]
        }
    ), &crate::utils::load_generic_model("sm")).map(|mut graphs| { assert_eq!(graphs.len(),1); graphs.remove(0) })
    .unwrap();

    graph.validate_full_numerator_tensor_network().unwrap();
    graph.global_prefactor.num = Atom::num(1);
    let error = graph.validate_full_numerator_tensor_network().unwrap_err();
    assert!(
        error
            .to_string()
            .contains("Full numerator still has dangling tensor indices"),
        "removing delta_i^j/N_c must expose the open fundamental color string: {error:?}"
    );
}

#[test]
fn vacuum_loop_edges_must_leave_a_spanning_forest() {
    test_initialise().unwrap();
    let graph: eyre::Result<Graph> = crate::finalized_runtime_dot!(digraph invalid_vacuum_basis {
 graph [projector=1]
        edge [num=1 mass=2]
        node [num=1]
        a -> b [id=0 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:0}"]
        b -> a [id=1 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:1}"]
    }, &crate::utils::load_generic_model("sm"));
    let error = graph
        .err()
        .expect("an overcomplete loop basis must be rejected");
    assert!(format!("{error:?}").contains("loop edges specified are not actual loop edges"));
}

#[test]
fn vacuum_tadpoles_keep_all_loop_generators_without_external_legs() {
    use linnet::half_edge::subgraph::SubSetLike;

    test_initialise().unwrap();
    let graphs: Vec<Graph> = crate::finalized_runtime_dot!(
        digraph canonical_vacuum {
 graph [projector=1]
            edge [num=1 mass=2]
            node [num=1]
            a -> a [id=0 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:1}"]
            a -> a [id=1 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:3}"]
        }
        digraph specified_vacuum {
 graph [projector=1]
            edge [num=1 mass=2]
            node [num=1]
            a -> a [id=0 lmb_id=1 source="{ufo_order:4}" sink="{ufo_order:5}"]
            a -> a [id=1 lmb_id=0 source="{ufo_order:6}" sink="{ufo_order:7}"]
        }
        digraph disconnected_vacuum {
 graph [projector=1]
            edge [num=1 mass=2]
            node [num=1]
            a -> a [id=0 lmb_id=1 source="{ufo_order:8}" sink="{ufo_order:9}"]
            b -> b [id=1 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:1}"]
        }
    , &crate::utils::load_generic_model("sm"))
    .unwrap();

    for (graph, expected_edges) in graphs.iter().zip([[0, 1], [1, 0], [1, 0]]) {
        let lmb = &graph.loop_momentum_basis;
        assert!(graph.initial_state_cut.is_empty());
        assert!(lmb.tree.is_empty());
        assert!(lmb.ext_edges.is_empty());
        assert_eq!(lmb.loop_edges.raw, expected_edges.map(EdgeIndex));
        for (loop_id, edge) in lmb.loop_edges.iter().enumerate() {
            let signature = &lmb.edge_signatures[*edge];
            assert!(signature.external.is_empty());
            assert_eq!(
                signature
                    .internal
                    .iter()
                    .map(|sign| *sign * 1_i32)
                    .collect::<Vec<_>>(),
                (0..2)
                    .map(|index| i32::from(index == loop_id))
                    .collect::<Vec<_>>()
            );
        }
        let exported = graph.dot_serialize(&DotExportSettings::default());
        let reparsed: Graph = exported
            .into_finalized_runtime_graph(&crate::utils::load_generic_model("sm"))
            .unwrap();
        assert_eq!(reparsed.loop_momentum_basis, *lmb);
    }
}
