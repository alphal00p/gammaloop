use feynkit_graph::{DiagramEdge, EdgeEndpoints, EdgeId, VertexId};
use linnet::half_edge::{
    involution::{EdgeIndex, HedgePair},
    subgraph::SubSetLike,
};
use typed_index_collections::ti_vec;

use super::Graph;
use crate::{
    dot,
    graph::{
        FeynmanGraph, GraphGroup,
        parse::{IntoGraph, complete_group_parsing},
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
            ext -> v4 [id=0]
            v5 -> ext [id=1]
            v6 -> ext [id=2]
            v4 -> v5 [id=3{}{}]
            v5 -> v6 [id=4{}]
            v6 -> v4 [id=5{}]
        }}"#,
        loop_attribute(3),
        if loop_edge == 3 { mass.as_str() } else { "" },
        loop_attribute(4),
        loop_attribute(5),
    )
}

fn parse_triangle_fixture(loop_edge: usize, loop_edge_mass: Option<i64>) -> Graph {
    triangle_fixture(loop_edge, loop_edge_mass)
        .into_graph(&crate::utils::load_generic_model("scalars"))
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
        r#" is_cut="0""#
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
            ext_in -> A [id=0 lmb_rep="P(0,a___)"{connection}]
            A -> B [id=1 lmb_id=0 lmb_rep="K(0,a___)"]
            A -> B [id=2 lmb_rep="-K(0,a___)+P(0,a___)"]
            B -> ext_out [id=3 lmb_rep="{outgoing_momentum}"{connection}]
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
        in0 -> A [id=0 lmb_rep="P(0,a___)" is_cut="0"]
        in1 -> A [id=1 lmb_rep="P(1,a___)" is_cut="1"]
        A -> B [id=2 lmb_rep="P(0,a___)+P(1,a___)"]
        B -> out0 [id=3 lmb_rep="P(0,a___)" is_cut="0"]
        B -> out1 [id=4 lmb_rep="P(1,a___)" is_cut="1"]
    }"#
}

#[test]
fn explicitly_marked_split_forward_endpoints_are_one_connection() {
    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("scalars");
    let graph: Graph = split_forward_fixture("P(0,a___)", true)
        .into_graph(&model)
        .unwrap();

    assert!(graph.get_external_signature().is_empty());
    assert_eq!(graph.loop_momentum_basis.ext_edges.len(), 1);
    assert_eq!(graph.loop_momentum_basis.loop_edges.len(), 1);

    let serialized = graph.dot_serialize(&DotExportSettings {
        split_xs_by_initial_states: true,
        ..DotExportSettings::default()
    });
    assert_eq!(serialized.matches("is_cut=").count(), 2, "{serialized}");

    let round_trip: Graph = serialized.into_graph(&model).unwrap();
    assert!(round_trip.get_external_signature().is_empty());
    assert_eq!(round_trip.loop_momentum_basis.ext_edges.len(), 1);
}

#[test]
fn equal_external_momenta_without_explicit_markers_remain_distinct_legs() {
    test_initialise().unwrap();
    let graph: Graph = split_forward_fixture("P(0,a___)", false)
        .into_graph(&crate::utils::load_generic_model("scalars"))
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
        .into_graph(&model)
        .unwrap();
    assert!(graph.get_external_signature().is_empty());
    assert_eq!(graph.loop_momentum_basis.ext_edges.len(), 2);

    let serialized = graph.dot_serialize(&DotExportSettings {
        split_xs_by_initial_states: true,
        ..DotExportSettings::default()
    });
    assert_eq!(serialized.matches("is_cut=").count(), 4, "{serialized}");

    let round_trip: Graph = serialized.into_graph(&model).unwrap();
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
    let g: Graph = dot!(
        digraph G {
            graph [projector="1"]
            edge [pdg=1000 num="1" dir=none]
            node [num="1"]
            ext0 [style=invis]
            ext1 [style=invis]
            ext0 -> A [id=0]
            B -> ext1 [id=1]
            A -> B [id=2, lmb_id=0 num="Q(2,spenso::mink(4,edge(2,1)))*Q(0,spenso::mink(4,edge(2,1)))"]
            A -> B [id=3]
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
            ext -> A
            A -> B [lmb_id=0]
            A -> B
            B -> ext
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
        .into_graph(&crate::utils::load_generic_model("scalars"))
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
fn gammaloop_attributes_round_trip_through_dot_export() {
    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("sm");
    let g: Graph = dot!(
        digraph attribute_graph {
            graph [projector="1" params="x;y"]
            ext_in [style=invis]
            ext_out [style=invis]
            A [num=1 dod=1]
            ext_in -> A [id=0 particle=H mass=0 name=e_in dod=-3 vakint_edge_power=2]
            A -> ext_out [id=1 particle=H mass=0 name=e_out momtrop_edge_power="1/2"]
        }
    )
    .unwrap();

    let serialized = g.dot_serialize(&DotExportSettings::default());
    let round_trip: Graph = serialized.as_str().into_graph(&model).unwrap();
    for graph in [&g, &round_trip] {
        let edges = graph
            .underlying
            .iter_edges()
            .map(|(_, _, edge)| edge.data)
            .collect::<Vec<_>>();
        assert_eq!(edges[0].name.value, "e_in");
        assert_eq!(edges[1].name.value, "e_out");
        assert!(edges.iter().all(|edge| edge.mass_atom(&model).is_zero()));
        assert_eq!(edges[0].dod.value, -3);
        assert_eq!(edges[0].extra_data.vakint_edge_power, Some(2));
        assert_eq!(
            edges[1].extra_data.momtrop_edge_power,
            Some(symbolica::parse!("1/2"))
        );
        assert_eq!(
            graph.underlying[linnet::half_edge::NodeIndex(0)].dod.value,
            1
        );
        assert_eq!(graph.param_builder.pairs.additional_params.params.len(), 2);
    }
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
    assert_eq!(runtime.underlying.n_live(), source.n_live());
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

fn native_fragment_fixture(
    vertex_numerator: symbolica::atom::Atom,
    dangling: &[symbolica::atom::Atom],
    paired: bool,
) -> feynkit_graph::FeynmanDiagram {
    use feynkit_graph::{DiagramVertex, FeynmanDiagram};
    let model = std::sync::Arc::new(crate::utils::load_generic_model("scalars"));
    let particle = model.particle_id("scalar_1").unwrap();
    let mut builder = FeynmanDiagram::builder(model, "native_fragments");
    let vertex = builder.add_vertex(DiagramVertex {
        name: "v".into(),
        interaction: None,
        numerator: vertex_numerator,
    });
    if paired {
        builder
            .add_edge(vertex, vertex, DiagramEdge::new(particle, false))
            .unwrap();
    }
    for numerator in dangling {
        let mut edge = DiagramEdge::new(particle, false);
        edge.is_dummy = true;
        edge.numerator = numerator.clone();
        builder.add_edge(Some(vertex), None, edge).unwrap();
    }
    let diagram = builder.build().unwrap();
    diagram.validate().unwrap();
    diagram
}

#[test]
fn native_runtime_conversion_certifies_local_fragments() {
    use symbolica::{atom::Atom, parse};
    test_initialise().unwrap();
    let diagram = native_fragment_fixture(Atom::one(), &[], true);
    let parsed = super::ParseGraph::from_feynkit_diagram(&diagram).unwrap();
    let vertex = parsed.graph.iter_node_ids().next().unwrap();
    let edge = EdgeIndex(0);

    for changed_vertex in [true, false] {
        let mut changed = parsed.clone();
        if changed_vertex {
            changed.graph[vertex].num = Some(Atom::num(2));
        } else {
            changed.graph[&edge].0.num = Some(Atom::num(2));
        }
        let error = Graph::from_feynkit_parse_graph(&diagram, changed)
            .err()
            .expect("runtime conversion must reject the changed diagram");
        assert!(
            error
                .to_string()
                .contains(if changed_vertex { "vertex" } else { "edge" }),
            "{error:?}"
        );
    }
    let mut compensated = parsed;
    compensated.graph[vertex].num = Some(Atom::num(2));
    compensated.graph[&edge].0.num = Some(Atom::one() / 2);
    assert!(
        Graph::from_feynkit_parse_graph(&diagram, compensated)
            .err()
            .expect("runtime conversion must reject the changed diagram")
            .to_string()
            .contains("vertex")
    );

    let factored = parse!("x*(y+z)", default_namespace = "runtime_fragment");
    let expanded = parse!("x*y+x*z", default_namespace = "runtime_fragment");
    // The conservative coefficient proof cancels factored minus distributed,
    // while the reverse remains inconclusive. Retain the former aggregate
    // comparison's runtime-minus-finalized operand order in both cases.
    let diagram = native_fragment_fixture(expanded.clone(), &[], true);
    let mut parsed = super::ParseGraph::from_feynkit_diagram(&diagram).unwrap();
    parsed.graph[vertex].num = Some(factored.clone());
    Graph::from_feynkit_parse_graph(&diagram, parsed).unwrap();

    let diagram = native_fragment_fixture(factored, &[], true);
    let mut parsed = super::ParseGraph::from_feynkit_diagram(&diagram).unwrap();
    parsed.graph[vertex].num = Some(expanded);
    assert!(
        Graph::from_feynkit_parse_graph(&diagram, parsed)
            .err()
            .expect("runtime conversion must reject an unproved fragment identity")
            .to_string()
            .contains("vertex")
    );

    let diagram = native_fragment_fixture(
        parse!("(x+y)^2", default_namespace = "runtime_fragment"),
        &[],
        true,
    );
    let mut parsed = super::ParseGraph::from_feynkit_diagram(&diagram).unwrap();
    parsed.graph[vertex].num = Some(parse!(
        "x^2+2*x*y+y^2",
        default_namespace = "runtime_fragment"
    ));
    assert!(
        Graph::from_feynkit_parse_graph(&diagram, parsed)
            .err()
            .expect("runtime conversion must reject the changed diagram")
            .to_string()
            .contains("vertex")
    );
}

#[test]
fn native_runtime_conversion_rejects_changed_frame_before_cut_indexing() {
    use linnet::half_edge::{EdgeAccessors, involution::Orientation};
    use symbolica::atom::Atom;
    test_initialise().unwrap();
    let diagram = native_fragment_fixture(Atom::one(), &[], true);
    let parsed = super::ParseGraph::from_feynkit_diagram(&diagram).unwrap();
    let mut changed = parsed.clone();
    let vertex = changed.graph.iter_nodes().next().unwrap().2.clone();
    changed.graph.add_node(vertex);
    assert!(
        Graph::from_feynkit_parse_graph(&diagram, changed)
            .err()
            .expect("runtime conversion must reject the changed diagram")
            .to_string()
            .contains("changed its finalized topology")
    );
    let mut changed = parsed;
    changed
        .graph
        .set_orientation(EdgeIndex(0), Orientation::Default);
    assert!(
        Graph::from_feynkit_parse_graph(&diagram, changed)
            .err()
            .expect("runtime conversion must reject the changed diagram")
            .to_string()
            .contains("changed its finalized topology")
    );
}

#[test]
fn native_runtime_conversion_preserves_dangling_dummy_product_acceptance() {
    use symbolica::atom::Atom;
    test_initialise().unwrap();
    for paired in [false, true] {
        for (vertex, omitted, accepted) in [
            (Atom::one(), vec![Atom::one()], true),
            (Atom::one(), vec![Atom::num(2)], false),
            (Atom::one(), vec![Atom::num(2), Atom::one() / 2], true),
            (Atom::Zero, vec![Atom::num(2)], true),
        ] {
            let diagram = native_fragment_fixture(vertex, &omitted, paired);
            let runtime = Graph::from_feynkit(&diagram, None, true);
            assert_eq!(
                runtime.is_ok(),
                accepted,
                "paired={paired}, omitted={omitted:?}, error={:?}",
                runtime.as_ref().err()
            );
        }
    }
}

#[test]
fn native_amplitude_runtime_enrichment_preserves_every_graph_coordinate() {
    use feynkit_generator::{Loops, Naming, NumeratorOptions, Process, Run, TopologyRules};
    test_initialise().unwrap();
    let model = std::sync::Arc::new(crate::utils::load_generic_model("scalars"));
    let run = Run::default().threads(1);
    let generated = Process::new(["scalar_1"], ["scalar_1"])
        .amplitude(Loops::total(1..=1))
        .with_topology(
            TopologyRules::new()
                .with_self_loops(false)
                .with_max_vertices(2),
        )
        .generate(model, &run)
        .unwrap()
        .build(
            &NumeratorOptions::default().with_naming(Naming::prefix("D")),
            &run,
        )
        .unwrap();
    assert!(!generated.is_empty());
    for diagram in generated.iter() {
        assert_native_runtime_parity(diagram);
    }
}

#[test]
fn native_cross_section_runtime_enrichment_preserves_cuts_and_positive_external_flow() {
    use feynkit_generator::{Cuts, Loops, Naming, NumeratorOptions, Process, Run, TopologyRules};
    test_initialise().unwrap();
    let model = std::sync::Arc::new(crate::utils::load_generic_model("scalars"));
    let run = Run::default().threads(1);
    let generated = Process::new(["scalar_1"], ["scalar_1", "scalar_1"])
        .cross_section(Loops::forward(1..=1))
        .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
        .with_topology(
            TopologyRules::new()
                .with_self_loops(false)
                .with_max_vertices(2),
        )
        .generate(model, &run)
        .unwrap()
        .build(
            &NumeratorOptions::default().with_naming(Naming::prefix("D")),
            &run,
        )
        .unwrap();
    assert!(!generated.is_empty());
    for diagram in generated.iter() {
        assert!(!diagram.cuts().is_empty());
        let mut changed = super::ParseGraph::from_feynkit_diagram(diagram).unwrap();
        let vertex = changed.graph.iter_nodes().next().unwrap().2.clone();
        changed.graph.add_node(vertex);
        assert!(
            Graph::from_feynkit_parse_graph(diagram, changed)
                .err()
                .expect("runtime conversion must reject the changed diagram")
                .to_string()
                .contains("changed its finalized topology")
        );
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
        assert_native_runtime_parity(diagram);
    }
}

#[test]
fn native_compton_runtime_enrichment_preserves_sewn_fermions() {
    use feynkit_generator::{
        Cuts, Loops, Naming, NumeratorOptions, Process, Restrict, Run, TopologyRules,
    };
    test_initialise().unwrap();
    let model = std::sync::Arc::new(
        feynkit_model::Model::from_json(include_str!(
            "../../../../feynkit-model/tests/fixtures/sm.json"
        ))
        .unwrap(),
    );

    let run = Run::default().threads(1);
    for pdg in [11_i64, -11] {
        let generated = Process::new([pdg, 22], [pdg, 22])
            .cross_section(Loops::forward(1..=1))
            .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
            .with_topology(
                TopologyRules::new()
                    .with_self_loops(false)
                    .with_max_vertices(4),
            )
            .with_restrict(Restrict::new().only_vertices(["V_98"]))
            .generate(model.clone(), &run)
            .unwrap()
            .build(
                &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                &run,
            )
            .unwrap();
        assert_eq!(generated.len(), 4);
        for diagram in generated.iter() {
            assert_native_runtime_parity(diagram);
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
