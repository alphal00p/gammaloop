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
fn numerator_interface_source_roundtrip_keeps_two_identical_open_axes() {
    use super::NumeratorInterface;
    use idenso::tensor::SymbolicTensor;
    use spenso::structure::partial::{
        OpenPortId, PartialIndex, PartialStructure, PartialStructureExt,
    };
    test_initialise().unwrap();
    let representation = Minkowski {}.new_rep(4).to_lib();
    let interface = PartialStructure::from_logical_slots([
        representation.slot(PartialIndex::Open(OpenPortId(0))),
        representation.slot(PartialIndex::Open(OpenPortId(1))),
    ]);
    let zero = SymbolicTensor::from_declared_parts(Atom::Zero, interface).unwrap();
    let (expression, metadata) = NumeratorInterface::export(&zero);
    let metadata = metadata.expect("a ranked zero requires interface metadata");
    let restored = NumeratorInterface::admit(
        expression,
        Some(NumeratorInterface::parse(&metadata.to_json()).unwrap()),
    )
    .unwrap();
    assert!(restored.same_value(&zero));
    let slots = restored.structure().logical_slots();
    assert_eq!(slots.len(), 2);
    assert_eq!(slots[0].aind, PartialIndex::Open(OpenPortId(0)));
    assert_eq!(slots[1].aind, PartialIndex::Open(OpenPortId(1)));
    assert!(
        NumeratorInterface::admit(Atom::one(), None)
            .unwrap()
            .is_scalar()
    );
    assert!(
        NumeratorInterface::admit(Atom::Zero, None)
            .unwrap()
            .is_scalar()
    );
    assert!(NumeratorInterface::admit(Atom::one(), Some(metadata)).is_err());
}

#[test]
fn runtime_graph_rejects_anonymous_ports_without_assigning_indices() {
    use super::NumeratorInterface;
    use spenso::structure::{abstract_index::AbstractIndex, slot::IsAbstractSlot};

    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("scalars");
    let rep = Minkowski {}.new_rep(4).to_lib();
    let anonymous = AbstractIndex::Open {
        owner: 801,
        axis: 0,
    };
    for port in [
        rep.to_symbolic([]),
        rep.slot::<AbstractIndex, _>(anonymous).to_atom(),
        rep.slot::<AbstractIndex, _>(anonymous.scoped(symbol!("runtime_open_scope")))
            .to_atom(),
    ] {
        let interface = NumeratorInterface(vec![port]);
        let metadata = serde_json::to_string(&interface.to_json()).unwrap();
        // A portable ranked zero is valid source data. GammaLoop requires an
        // explicit graph interface even when its expression vanishes.
        assert_eq!(
            NumeratorInterface::admit(Atom::Zero, Some(interface))
                .unwrap()
                .rank(),
            1
        );
        for source in [
            triangle_fixture(3, None).replace(
                "graph [projector=\"1\"]",
                &format!(
                    "graph [projector=\"1\" num=\"0\" numerator_prefactor_interface={metadata}]"
                ),
            ),
            triangle_fixture(3, None).replace(
                "node [num=\"1\"]",
                &format!("node [num=\"0\" num_interface={metadata}]"),
            ),
            triangle_fixture(3, None).replace(
                "edge [pdg=1000 num=\"1\" dir=none]",
                &format!("edge [pdg=1000 num=\"0\" num_interface={metadata} dir=none]"),
            ),
        ] {
            let result: color_eyre::Result<Graph> = source.into_finalized_runtime_graph(&model);
            let error = result
                .err()
                .expect("anonymous runtime ports must be rejected");
            assert!(
                format!("{error:#}").contains("requires an explicit index"),
                "{error:?}"
            );
        }
    }
    // Compact internal contractions are scalar; fully indexed interfaces do
    // not require expanding tensor shorthand or changing the source expression.
    let vector = GS.emr_vec(EdgeIndex(3), rep.to_symbolic([]));
    let scalar = SPENSO_TAG.dot.call_args([vector.clone(), vector]);
    let source = triangle_fixture(3, None).replace(
        "graph [projector=\"1\"]",
        &format!(
            "graph [projector=\"1\" num={}]",
            serde_json::to_string(&scalar.to_canonical_string()).unwrap()
        ),
    );
    let graph: Graph = source.into_finalized_runtime_graph(&model).unwrap();
    assert!(graph.global_prefactor.num.is_scalar());
    assert_eq!(
        NumeratorInterface::export(&graph.global_prefactor.num).0,
        scalar
    );
}

#[test]
fn runtime_and_canonical_dot_preserve_ranked_zero_fragments() {
    use super::NumeratorInterface;
    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("scalars");
    let canonical = Graph::import_dot(
        include_str!("../../../../../tests/resources/graphs/ir_safe_thresholds/no_directive_cross_section.dot"),
        &model,
    ).unwrap().remove(0);
    let metric = symbolica::function!(
        spenso::network::library::symbolic::ETS.metric,
        Minkowski {}.new_rep(4).to_symbolic([Atom::num(37)]),
        Minkowski {}.new_rep(4).to_symbolic([Atom::num(41)])
    );
    let metric = NumeratorInterface::admit(metric, None)
        .unwrap()
        .try_into_indexed()
        .unwrap();
    let zero = metric.with_algebra_result(Atom::Zero).unwrap();
    for (mut graph, edge) in [
        (parse_triangle_fixture(3, None), EdgeIndex(3)),
        (canonical, EdgeIndex(2)),
    ] {
        assert!(
            !graph
                .dot_serialize(&DotExportSettings::default())
                .contains("num_interface")
        );
        let vertex = graph.underlying.iter_nodes().next().unwrap().0;
        graph.underlying[vertex].num.value = zero.clone();
        graph.underlying[edge].num.value = metric.clone();
        let serialized = graph.dot_serialize(&DotExportSettings::default());
        assert!(serialized.contains("num_interface"));
        let restored = Graph::import_dot(&serialized, &model).unwrap().remove(0);
        assert!(restored.underlying[vertex].num.value.same_value(&zero));
        assert!(restored.underlying[edge].num.value.same_value(&metric));
        assert_eq!(restored.finalized_cuts.len(), graph.finalized_cuts.len());
        assert_eq!(restored.overall_factor, graph.overall_factor);
        let numerator = restored
            .numerator(&restored.full_filter(), &restored.empty_subgraph())
            .unwrap()
            .state
            .expr;
        assert!(numerator.is_zero());
        assert!(
            numerator.is_scalar(),
            "the zero fragment must still close its paired ports"
        );
    }
}

#[test]
fn global_tensor_json_and_dot_roundtrips_preserve_zero_ports_and_logical_order() {
    use super::NumeratorInterface;
    use crate::numerator::GlobalPrefactor;
    use crate::uv::UltravioletGraph;
    use idenso::tensor::SymbolicTensor;
    use spenso::structure::{IndexedStructure, abstract_index::AbstractIndex};
    test_initialise().unwrap();
    let representation = Minkowski {}.new_rep(4).to_lib();
    let interface: IndexedStructure = [41, 37]
        .map(|index| representation.slot(AbstractIndex::Normal(index)))
        .into_iter()
        .collect();
    let zero = SymbolicTensor::from_indexed_parts(Atom::Zero, interface.clone()).unwrap();
    let metric = NumeratorInterface::admit(
        symbolica::function!(
            spenso::network::library::symbolic::ETS.metric,
            representation.to_symbolic([Atom::num(37)]),
            representation.to_symbolic([Atom::num(41)])
        ),
        None,
    )
    .unwrap();
    let metric =
        SymbolicTensor::from_indexed_parts(metric.expression().clone(), interface.clone()).unwrap();
    let model = crate::utils::load_generic_model("scalars");
    let canonical = Graph::import_dot(
        include_str!("../../../../../tests/resources/graphs/ir_safe_thresholds/no_directive_cross_section.dot"),
        &model,
    ).unwrap().remove(0);
    for mut graph in [parse_triangle_fixture(3, None), canonical] {
        graph.global_prefactor.num = zero.clone();
        graph.global_prefactor.projector = metric.clone();
        let json = serde_json::to_value(&graph.global_prefactor).unwrap();
        assert!(json.get("num_interface").is_some());
        assert_eq!(
            symbolica::try_parse!(json["projector"].as_str().unwrap()).unwrap(),
            NumeratorInterface::export(&metric).0,
            "serialized tensor payloads must retain source syntax independently of display printers"
        );
        assert_eq!(
            symbolica::try_parse!(json["num"].as_str().unwrap()).unwrap(),
            Atom::Zero
        );
        let restored_global: GlobalPrefactor = serde_json::from_value(json).unwrap();
        assert!(restored_global.num.same_value(&zero));
        assert!(restored_global.projector.same_value(&metric));
        let product = restored_global
            .num
            .multiply(restored_global.projector)
            .unwrap();
        assert!(product.is_zero() && product.is_scalar());

        let dot = graph.dot_serialize(&DotExportSettings::default());
        assert!(dot.contains("numerator_prefactor_interface"));
        let restored = Graph::import_dot(&dot, &model).unwrap().remove(0);
        assert!(restored.global_prefactor.num.same_value(&zero));
        assert!(restored.global_prefactor.projector.same_value(&metric));
        assert_eq!(restored.global_prefactor.num.rank(), 2);
        assert_eq!(
            restored.global_prefactor.num.structure().logical_slots(),
            interface.logical_slots()
        );
        assert_eq!(
            restored
                .global_prefactor
                .projector
                .structure()
                .logical_slots(),
            interface.logical_slots()
        );
        assert_eq!(restored.finalized_cuts.len(), graph.finalized_cuts.len());
        let full = restored
            .numerator(&restored.full_filter(), &restored.empty_subgraph())
            .unwrap()
            .state
            .expr
            .multiply(restored.global_numerator().unwrap())
            .unwrap();
        assert!(full.is_zero() && full.is_scalar());
        restored.validate_full_numerator_tensor_network().unwrap();
    }
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
        assert_eq!(
            runtime_data.num.value,
            data.data.numerator.clone().try_into_indexed().unwrap()
        );
        assert_eq!(runtime_data.is_dummy, data.data.is_dummy);
    }
    for index in 0..source.n_hedges() {
        let hedge = Hedge(index);
        assert_eq!(runtime.underlying.node_id(hedge), source.node_id(hedge));
        assert_eq!(runtime.underlying.flow(hedge), source.flow(hedge));
        assert_eq!(runtime.underlying[&hedge], source[&hedge]);
    }
    assert_eq!(runtime.overall_factor, *diagram.overall_factor());
    assert_eq!(
        runtime.global_prefactor.num,
        diagram
            .numerator_prefactor()
            .clone()
            .try_into_indexed()
            .unwrap()
    );
    assert_eq!(
        runtime.global_prefactor.projector,
        diagram.projector().clone().try_into_indexed().unwrap()
    );
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
            .unwrap()
            .get_single_atom()
            .unwrap()
            .expand(),
        super::NumeratorInterface::export(diagram.numerator())
            .0
            .expand()
    );
    for (_, edge, _) in source.iter_edges() {
        let mut selected: linnet::half_edge::subgraph::SuBitGraph = source.empty_subgraph();
        use linnet::half_edge::subgraph::ModifySubSet;
        selected.add(source[&edge].1);
        assert_eq!(
            runtime
                .numerator(&selected, &empty)
                .unwrap()
                .get_single_atom()
                .unwrap(),
            super::NumeratorInterface::export(&diagram.numerator_of(&selected, &empty).unwrap()).0
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

impl Graph {
    fn assert_fixture_hedge_routing(&self, expected: &[(u8, Option<u16>)]) {
        use spenso::structure::TensorStructure;

        // Cut permutations change graph positions while the tensor slots
        // retain their parser indices. Generated UFO orders are omitted
        // from DOT payloads, so check both pieces of routing directly.
        assert_eq!(self.underlying.n_hedges(), expected.len());
        for ((hedge, data), (order, index)) in self.underlying.iter_hedges().zip(expected) {
            assert_eq!(data.ufo_order.value, *order, "UFO order at {hedge}");
            let expected_index = index.map(|index| Aind::Hedge(index, 0));
            let mut has_indices = false;
            for structure in [
                &data.num_indices.color_indices.edge_indices,
                &data.num_indices.color_indices.vertex_indices,
                &data.num_indices.spin_indices.edge_indices,
                &data.num_indices.spin_indices.vertex_indices,
            ] {
                for slot in structure.external_structure_iter() {
                    has_indices = true;
                    assert_eq!(Some(slot.aind), expected_index, "tensor slot at {hedge}");
                }
            }
            assert_eq!(
                has_indices,
                expected_index.is_some(),
                "slot presence at {hedge}"
            );
        }
    }
}

#[test]
fn xs_parsing() {
    test_initialise().unwrap();
    let mut model = load_generic_model("sm");
    let mut restriction =
        InputParamCard::from_str(LIGHT_QUARKS_MASSLESS_RESTRICTION.into(), "json").unwrap();
    model.simplify(&mut restriction).unwrap();
    let g: Graph = dot!(
        digraph{
          ext0	 [style=invis is_cut=0];
          ext1	 [style=invis is_cut=1];
          ext0	-> 0	 [dir=back particle="e+"];
          ext1	-> 0	 [particle="e-"];
            0	-> 1	 [particle="a"];
            1 -> 2	 [particle="d"];
            2 -> 4   [particle="g"];
            2 -> 3	 [particle="d"];
            3 -> 4	 [particle="d"];
            4 -> 1	 [particle="d"];
            3->5	 [particle="a"];
            5	-> ext0	 [dir=back particle="e+"];
            5	-> ext1	 [particle="e-"];
        }, &model
    )
    .unwrap();
    g.assert_fixture_hedge_routing(&[
        (0, Some(16)),
        (1, Some(17)),
        (0, Some(2)),
        (1, Some(3)),
        (0, Some(4)),
        (1, Some(5)),
        (2, Some(6)),
        (2, Some(7)),
        (0, Some(8)),
        (1, Some(9)),
        (2, Some(10)),
        (2, Some(11)),
        (0, Some(12)),
        (1, Some(13)),
        (1, Some(14)),
        (0, Some(15)),
        (2, Some(0)),
        (2, Some(1)),
    ]);
    assert_snapshot!(g.debug_dot_with_settings(&DotExportSettings {
            include_autogenerated_fields: true,
            ..DotExportSettings::default()
        }),@r#"
        digraph {
          num = "1";
          overall_factor = "1";
          overall_factor_evaluated = "1";
          projector = "u(1,spenso::bis(4,hedge(17)))*ubar(1,spenso::bis(4,hedge(15)))*v(0,spenso::bis(4,hedge(14)))*vbar(0,spenso::bis(4,hedge(16)))";
          threshold_counterterms = "
        schema_version = 1
        cuts = []
        ";
          "0"[dod="0" dod_autogen="true" int_id="V_98" name_autogen="true" num="UFO::GC_3*spenso::gamma(spenso::bis(4,hedge(16)),spenso::bis(4,hedge(17)),spenso::mink(4,hedge(0)))" num_autogen="true"];
          "1"[dod="0" dod_autogen="true" int_id="V_71" name_autogen="true" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(13)),spenso::dind(spenso::cof(3,hedge(2))))*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(13)),spenso::mink(4,hedge(1)))" num_autogen="true"];
          "2"[dod="0" dod_autogen="true" int_id="V_74" name_autogen="true" num="UFO::GC_11*spenso::gamma(spenso::bis(4,hedge(4)),spenso::bis(4,hedge(3)),spenso::mink(4,hedge(6)))*spenso::t(spenso::coad(8,hedge(6)),spenso::cof(3,hedge(3)),spenso::dind(spenso::cof(3,hedge(4))))" num_autogen="true"];
          "3"[dod="0" dod_autogen="true" int_id="V_71" name_autogen="true" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(5)),spenso::dind(spenso::cof(3,hedge(8))))*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(5)),spenso::mink(4,hedge(10)))" num_autogen="true"];
          "4"[dod="0" dod_autogen="true" int_id="V_74" name_autogen="true" num="UFO::GC_11*spenso::gamma(spenso::bis(4,hedge(12)),spenso::bis(4,hedge(9)),spenso::mink(4,hedge(7)))*spenso::t(spenso::coad(8,hedge(7)),spenso::cof(3,hedge(9)),spenso::dind(spenso::cof(3,hedge(12))))" num_autogen="true"];
          "5"[dod="0" dod_autogen="true" int_id="V_98" name_autogen="true" num="UFO::GC_3*spenso::gamma(spenso::bis(4,hedge(15)),spenso::bis(4,hedge(14)),spenso::mink(4,hedge(11)))" num_autogen="true"];

          "5":14	-> "0":0	 [id=0 dir=back  dod="-1" dod_autogen="true" is_cut="0" lmb_rep="P(0,a___)" name="e0" name_autogen="true" num="1" num_autogen="true" particle="e+"];
          "5":15	-> "0":1	 [id=1  dod="-1" dod_autogen="true" is_cut="1" lmb_rep="P(1,a___)" name="e1" name_autogen="true" num="1" num_autogen="true" particle="e-"];
          "2":4	-> "3":5	 [id=2  dod="-1" dod_autogen="true" lmb_rep="-1*K(1,a___)+K(0,a___)" name="e2" name_autogen="true" num="Q(2,spenso::mink(4,edge(2,1)))*spenso::g(spenso::cof(3,hedge(4)),spenso::dind(spenso::cof(3,hedge(5))))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(4)),spenso::mink(4,edge(2,1)))*𝑖" num_autogen="true" particle="d"];
          "2":6	-> "4":7	 [id=3 dir=none  dod="-2" dod_autogen="true" lmb_id="1" lmb_rep="K(1,a___)" name="e3" name_autogen="true" num="-𝑖*spenso::g(spenso::coad(8,hedge(6)),spenso::coad(8,hedge(7)))*spenso::g(spenso::mink(4,hedge(6)),spenso::mink(4,hedge(7)))" num_autogen="true" particle="g"];
          "3":8	-> "4":9	 [id=4  dod="-1" dod_autogen="true" lmb_rep="-1*K(1,a___)+-1*P(0,a___)+-1*P(1,a___)+K(0,a___)" name="e4" name_autogen="true" num="Q(4,spenso::mink(4,edge(4,1)))*spenso::g(spenso::cof(3,hedge(8)),spenso::dind(spenso::cof(3,hedge(9))))*spenso::gamma(spenso::bis(4,hedge(9)),spenso::bis(4,hedge(8)),spenso::mink(4,edge(4,1)))*𝑖" num_autogen="true" particle="d"];
          "3":10	-> "5":11	 [id=5 dir=none  dod="-2" dod_autogen="true" lmb_rep="P(0,a___)+P(1,a___)" name="e5" name_autogen="true" num="-𝑖*spenso::g(spenso::mink(4,hedge(10)),spenso::mink(4,hedge(11)))" num_autogen="true" particle="a"];
          "4":12	-> "1":13	 [id=6  dod="-1" dod_autogen="true" lmb_rep="-1*P(0,a___)+-1*P(1,a___)+K(0,a___)" name="e6" name_autogen="true" num="Q(6,spenso::mink(4,edge(6,1)))*spenso::g(spenso::cof(3,hedge(12)),spenso::dind(spenso::cof(3,hedge(13))))*spenso::gamma(spenso::bis(4,hedge(13)),spenso::bis(4,hedge(12)),spenso::mink(4,edge(6,1)))*𝑖" num_autogen="true" particle="d"];
          "0":16	-> "1":17	 [id=7 dir=none  dod="-2" dod_autogen="true" lmb_rep="P(0,a___)+P(1,a___)" name="e7" name_autogen="true" num="-𝑖*spenso::g(spenso::mink(4,hedge(0)),spenso::mink(4,hedge(1)))" num_autogen="true" particle="a"];
          "1":2	-> "2":3	 [id=8  dod="-1" dod_autogen="true" lmb_id="0" lmb_rep="K(0,a___)" name="e8" name_autogen="true" num="Q(8,spenso::mink(4,edge(8,1)))*spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(3))))*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(2)),spenso::mink(4,edge(8,1)))*𝑖" num_autogen="true" particle="d"];
        }
        "#);

    let g: Graph = dot!(
            digraph{
                edge [
                    pdg=1000
                ]
                ext0[style=invis is_cut=0];
                ext0->a
                a->b->c->d->a
                b->e
                c->e
                d->e
                c->ext0


            }
            ,"scalars")
    .unwrap();
    g.assert_fixture_hedge_routing(&[
        (2, None),
        (0, None),
        (1, None),
        (0, None),
        (2, None),
        (0, None),
        (1, None),
        (0, None),
        (2, None),
        (1, None),
        (3, None),
        (1, None),
        (1, None),
        (2, None),
        (2, None),
        (0, None),
    ]);
    assert_snapshot!(g.debug_dot_with_settings(&DotExportSettings {
            include_autogenerated_fields: true,
            ..DotExportSettings::default()
        }),@r#"
        digraph {
          num = "1";
          overall_factor = "1";
          overall_factor_evaluated = "1";
          projector = "1";
          threshold_counterterms = "
        schema_version = 1
        cuts = []
        ";
          a[dod="0" dod_autogen="true" int_id="V_3_SCALAR_000" num="UFO::SCALAR_COUPLING" num_autogen="true"];
          b[dod="0" dod_autogen="true" int_id="V_3_SCALAR_000" num="UFO::SCALAR_COUPLING" num_autogen="true"];
          c[dod="0" dod_autogen="true" int_id="V_4_SCALAR_0000" num="UFO::SCALAR_COUPLING" num_autogen="true"];
          d[dod="0" dod_autogen="true" int_id="V_3_SCALAR_000" num="UFO::SCALAR_COUPLING" num_autogen="true"];
          e[dod="0" dod_autogen="true" int_id="V_3_SCALAR_000" num="UFO::SCALAR_COUPLING" num_autogen="true"];

          c:10	-> a:0	 [id=0 dir=none  dod="-2" dod_autogen="true" is_cut="0" lmb_rep="P(0,a___)" name="e0" name_autogen="true" num="1" num_autogen="true" particle="scalar_0"];
          b:2	-> c:3	 [id=1 dir=none  dod="-2" dod_autogen="true" lmb_id="0" lmb_rep="K(0,a___)" name="e1" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
          b:4	-> e:5	 [id=2 dir=none  dod="-2" dod_autogen="true" lmb_rep="-1*K(0,a___)+K(2,a___)+P(0,a___)" name="e2" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
          c:6	-> d:7	 [id=3 dir=none  dod="-2" dod_autogen="true" lmb_rep="-1*K(1,a___)+-1*P(0,a___)+K(0,a___)" name="e3" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
          c:8	-> e:9	 [id=4 dir=none  dod="-2" dod_autogen="true" lmb_id="1" lmb_rep="K(1,a___)" name="e4" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
          a:15	-> b:1	 [id=5 dir=none  dod="-2" dod_autogen="true" lmb_rep="K(2,a___)+P(0,a___)" name="e5" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
          d:11	-> a:12	 [id=6 dir=none  dod="-2" dod_autogen="true" lmb_id="2" lmb_rep="K(2,a___)" name="e6" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
          d:13	-> e:14	 [id=7 dir=none  dod="-2" dod_autogen="true" lmb_rep="-1*K(1,a___)+-1*K(2,a___)+-1*P(0,a___)+K(0,a___)" name="e7" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
        }
        "#);

    let g: Graph = dot!(
            digraph GL06{
              num = "1";
              overall_factor = "(AutG(1))^(-1)*AntiFermionSpinSumSign(-1)*ExternalFermionOrderingSign(-1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2)";
              0[int_id=V_141];
              1[int_id=V_141];
              2[int_id=V_134];
              3[int_id=V_134];
              4[int_id=V_137];
              5[int_id=V_137];
              6[int_id=V_98];
              7[int_id=V_98];
              ext0 [style=invis];
              ext0-> 7:0 [id=0 dir=back is_cut=0 particle="e+"];
              ext1 [style=invis];
              ext1-> 7:1 [id=1 is_cut=1 particle="e-"];
              0:2-> 1:3 [id=2  particle="t"];
              0:4-> 1:5 [id=3 dir=none  particle="H"];
              5:6-> 0:7 [id=4  particle="t"];
              1:8-> 3:9 [id=5  particle="t"];
              4:10-> 2:11 [id=6  particle="t"];
              2:12-> 5:13 [id=7  particle="t"];
              2:14-> 7:15 [id=8 dir=none  particle="a"];
              3:16-> 4:17 [id=9  particle="t"];
              3:18-> 6:19 [id=10 dir=none  particle="a"];
              4:20-> 5:21 [id=11 dir=none  particle="g"];
              ext12 [style=invis];
              6:22-> ext12 [id=12 dir=back is_cut=0 particle="e+"];
              ext13 [style=invis];
              6:23-> ext13 [id=13 is_cut=1 particle="e-"];
            }, &model
        )
        .unwrap();

    println!("{}", g.global_prefactor.projector);

    let g: Graph = dot!(
    digraph GL06{
        num = "-2";
            exte0    [style=invis, is_cut=0];
            exte0   -> 7:0   [id=0 particle="e+"];
            exte1    [style=invis, is_cut=1];
            exte1   -> 7:1   [id=1 particle="e-"];
            4:22    -> 5:23  [id=2 particle="g"];
            3:20    -> 6:21  [id=3 particle="a"];
            0:4     -> 1:5   [id=4 particle="t" lmb_id=0];
            0:6     -> 1:7   [id=5 particle="H" lmb_id=1];
            5:8     -> 0:9   [id=6 particle="t"];
            1:10    -> 3:11  [id=7 particle="t"];
            4:12    -> 2:13  [id=8 particle="t"];
            2:14    -> 5:15  [id=9 particle="t" lmb_id=2];
            2:16    -> 7:17  [id=10 particle="a"];
            3:18    -> 4:19  [id=11 particle="t"];

            6:3     -> exte0        [id=12 particle="e+"];

            6:2     -> exte1        [id=13 particle="e-"];
    }, &model)
    .unwrap();

    println!("{}", g.global_prefactor.projector);
}

#[test]
fn xs_glueing() {
    test_initialise().unwrap();
    let mut model = load_generic_model("sm");
    let mut restriction =
        InputParamCard::from_str(LIGHT_QUARKS_MASSLESS_RESTRICTION.into(), "json").unwrap();
    // This fixture retains its massive strange propagator while u and d stay massless.
    restriction.remove(&"MS".into());
    restriction.remove(&"yms".into());
    model.simplify(&mut restriction).unwrap();
    let g: Graph = dot!(
        digraph{
          num = "1";
          overall_factor = "1";
          0[int_id=V_89];
          1[int_id=V_127];
          2[int_id=V_123];
          3[int_id=V_93];

          ext0	 [style=invis is_cut=0];
          ext1	 [style=invis is_cut=1];
          2	-> ext0	     [id=0  particle="d"];
          ext0	-> 0	 [id=1  particle="d"];
          3:2	-> ext1	 [id=2  particle="c"];
          ext1	-> 1:3	 [id=3  particle="c"];
          0:4	-> 2:5	 [id=4  particle="u"];
          1:6	-> 3:7	 [id=5  particle="s"];
          0:8	-> 1:9	 [id=6  particle="W-"];
          2:10	-> 3:11	 [id=7  particle="W+"];
        }, &model
    )
    .unwrap();
    g.assert_fixture_hedge_routing(&[
        (1, Some(1)),
        (1, Some(3)),
        (0, Some(2)),
        (0, Some(0)),
        (0, Some(4)),
        (1, Some(5)),
        (0, Some(6)),
        (1, Some(7)),
        (2, Some(8)),
        (2, Some(9)),
        (2, Some(10)),
        (2, Some(11)),
    ]);
    assert_snapshot!(g.debug_dot_with_settings(&DotExportSettings {
            include_autogenerated_fields: true,
            ..DotExportSettings::default()
        }),@r#"
        digraph {
          num = "1";
          overall_factor = "1";
          overall_factor_evaluated = "1";
          projector = "u(0,spenso::bis(4,hedge(1)))*u(1,spenso::bis(4,hedge(3)))*ubar(0,spenso::bis(4,hedge(0)))*ubar(1,spenso::bis(4,hedge(2)))";
          threshold_counterterms = "
        schema_version = 1
        cuts = []
        ";
          "0"[dod="0" dod_autogen="true" int_id="V_89" name_autogen="true" num="(-1/2*spenso::gamma(spenso::bis(4,hedge(4)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(8)))*spenso::gamma5(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(1)))+1/2*spenso::gamma(spenso::bis(4,hedge(4)),spenso::bis(4,hedge(1)),spenso::mink(4,hedge(8))))*UFO::GC_100*spenso::g(spenso::cof(3,hedge(1)),spenso::dind(spenso::cof(3,hedge(4))))" num_autogen="true"];
          "1"[dod="0" dod_autogen="true" int_id="V_127" name_autogen="true" num="(-1/2*spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(9)))*spenso::gamma5(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3)))+1/2*spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,hedge(3)),spenso::mink(4,hedge(9))))*UFO::GC_45*spenso::g(spenso::cof(3,hedge(3)),spenso::dind(spenso::cof(3,hedge(6))))" num_autogen="true"];
          "2"[dod="0" dod_autogen="true" int_id="V_123" name_autogen="true" num="(-1/2*spenso::gamma(spenso::bis(4,hedge(0)),spenso::bis(4,vertex(2,1)),spenso::mink(4,hedge(10)))*spenso::gamma5(spenso::bis(4,vertex(2,1)),spenso::bis(4,hedge(5)))+1/2*spenso::gamma(spenso::bis(4,hedge(0)),spenso::bis(4,hedge(5)),spenso::mink(4,hedge(10))))*UFO::GC_41*spenso::g(spenso::cof(3,hedge(5)),spenso::dind(spenso::cof(3,hedge(0))))" num_autogen="true"];
          "3"[dod="0" dod_autogen="true" int_id="V_93" name_autogen="true" num="(-1/2*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(3,1)),spenso::mink(4,hedge(11)))*spenso::gamma5(spenso::bis(4,vertex(3,1)),spenso::bis(4,hedge(7)))+1/2*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(7)),spenso::mink(4,hedge(11))))*UFO::GC_104*spenso::g(spenso::cof(3,hedge(7)),spenso::dind(spenso::cof(3,hedge(2))))" num_autogen="true"];

          "2":3	-> "0":0	 [id=0  dod="-1" dod_autogen="true" is_cut="0" lmb_rep="P(0,a___)" name="e0" name_autogen="true" num="spenso::g(spenso::cof(3,hedge(0)),spenso::dind(spenso::cof(3,hedge(1))))" num_autogen="true" particle="d"];
          "3":2	-> "1":1	 [id=1  dod="-1" dod_autogen="true" is_cut="1" lmb_rep="P(1,a___)" name="e1" name_autogen="true" num="spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(3))))" num_autogen="true" particle="c"];
          "0":8	-> "1":9	 [id=2 dir=none  dod="-2" dod_autogen="true" lmb_rep="-1*K(0,a___)+P(0,a___)" name="e2" name_autogen="true" num="-𝑖*spenso::g(spenso::mink(4,hedge(8)),spenso::mink(4,hedge(9)))" num_autogen="true" particle="W-"];
          "2":10	-> "3":11	 [id=3 dir=none  dod="-2" dod_autogen="true" lmb_rep="-1*P(0,a___)+K(0,a___)" name="e3" name_autogen="true" num="-𝑖*spenso::g(spenso::mink(4,hedge(10)),spenso::mink(4,hedge(11)))" num_autogen="true" particle="W+"];
          "0":4	-> "2":5	 [id=4  dod="-1" dod_autogen="true" lmb_id="0" lmb_rep="K(0,a___)" name="e4" name_autogen="true" num="Q(4,spenso::mink(4,edge(4,1)))*spenso::g(spenso::cof(3,hedge(4)),spenso::dind(spenso::cof(3,hedge(5))))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(4)),spenso::mink(4,edge(4,1)))*𝑖" num_autogen="true" particle="u"];
          "1":6	-> "3":7	 [id=5  dod="-1" dod_autogen="true" lmb_rep="-1*K(0,a___)+P(0,a___)+P(1,a___)" name="e5" name_autogen="true" num="(Q(5,spenso::mink(4,edge(5,1)))*spenso::gamma(spenso::bis(4,hedge(7)),spenso::bis(4,hedge(6)),spenso::mink(4,edge(5,1)))+UFO::MS*spenso::g(spenso::bis(4,hedge(6)),spenso::bis(4,hedge(7))))*spenso::g(spenso::cof(3,hedge(6)),spenso::dind(spenso::cof(3,hedge(7))))*𝑖" num_autogen="true" particle="s"];
        }
        "#);

    let g: Graph = dot!(
        digraph{
            a->b [particle= c];
            b->a [particle= d];
            b->a [particle= "W+"];
            a[int_id=V_90]
        }, &model
    )
    .unwrap();
    g.assert_fixture_hedge_routing(&[
        (0, Some(0)),
        (1, Some(1)),
        (0, Some(2)),
        (1, Some(3)),
        (2, Some(4)),
        (2, Some(5)),
    ]);
    assert_snapshot!(g.debug_dot_with_settings(&DotExportSettings {
            include_autogenerated_fields: true,
            ..DotExportSettings::default()
        }),@r#"
        digraph {
          num = "1";
          overall_factor = "1";
          overall_factor_evaluated = "1";
          projector = "1";
          threshold_counterterms = "
        schema_version = 1
        cuts = []
        ";
          a[dod="0" dod_autogen="true" int_id="V_90" num="(-1/2*spenso::gamma(spenso::bis(4,hedge(0)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(5)))*spenso::gamma5(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(3)))+1/2*spenso::gamma(spenso::bis(4,hedge(0)),spenso::bis(4,hedge(3)),spenso::mink(4,hedge(5))))*UFO::GC_103*spenso::g(spenso::cof(3,hedge(3)),spenso::dind(spenso::cof(3,hedge(0))))" num_autogen="true"];
          b[dod="0" dod_autogen="true" int_id="V_126" num="(-1/2*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(4)))*spenso::gamma5(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(1)))+1/2*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(1)),spenso::mink(4,hedge(4))))*UFO::GC_44*spenso::g(spenso::cof(3,hedge(1)),spenso::dind(spenso::cof(3,hedge(2))))" num_autogen="true"];

          a:0	-> b:1	 [id=0  dod="-1" dod_autogen="true" lmb_id="0" lmb_rep="K(0,a___)" name="e0" name_autogen="true" num="(Q(0,spenso::mink(4,edge(0,1)))*spenso::gamma(spenso::bis(4,hedge(1)),spenso::bis(4,hedge(0)),spenso::mink(4,edge(0,1)))+UFO::MC*spenso::g(spenso::bis(4,hedge(0)),spenso::bis(4,hedge(1))))*spenso::g(spenso::cof(3,hedge(0)),spenso::dind(spenso::cof(3,hedge(1))))*𝑖" num_autogen="true" particle="c"];
          b:2	-> a:3	 [id=1  dod="-1" dod_autogen="true" lmb_id="1" lmb_rep="K(1,a___)" name="e1" name_autogen="true" num="Q(1,spenso::mink(4,edge(1,1)))*spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(3))))*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(2)),spenso::mink(4,edge(1,1)))*𝑖" num_autogen="true" particle="d"];
          b:4	-> a:5	 [id=2 dir=none  dod="-2" dod_autogen="true" lmb_rep="-1*K(1,a___)+K(0,a___)" name="e2" name_autogen="true" num="-𝑖*spenso::g(spenso::mink(4,hedge(4)),spenso::mink(4,hedge(5)))" num_autogen="true" particle="W+"];
        }
        "#);
}

#[test]
fn vertex_normalization() {
    test_initialise().unwrap();
    let mut model = load_generic_model("sm");
    let mut restriction =
        InputParamCard::from_str(LIGHT_QUARKS_MASSLESS_RESTRICTION.into(), "json").unwrap();
    model.simplify(&mut restriction).unwrap();
    let g: Graph = dot!(
            digraph GL1{
                num = "1";
                overall_factor = "AutG(1)^-1*InternalFermionLoopSign(-1)*ExternalFermionOrderingSign(1)*AntiFermionSpinSumSign(1)";
                0[int_id=V_117];
                1[int_id=V_82];
                2[int_id=V_71];
                3[int_id=V_11];

                ext0 [style=invis];
                2:0 -> ext0 [id=0 dir=back is_cut=0 particle="a"];
                ext1 [style=invis];
                ext1 -> 3:1 [id=1 is_cut=0 particle="a"];
                0:2 -> 1:3 [id=2 particle="c~"];
                0:4 -> 2:5 [id=3 particle="d"];
                0:6 -> 3:7 [id=4 particle="G+"];
                1:8 -> 2:9 [id=5 particle="d~"];
                1:10 -> 3:11 [id=6 particle="G-"];
            }, &model
        ).unwrap();
    g.assert_fixture_hedge_routing(&[
        (0, Some(1)),
        (2, Some(0)),
        (1, Some(2)),
        (0, Some(3)),
        (0, Some(4)),
        (1, Some(5)),
        (2, None),
        (2, None),
        (1, Some(8)),
        (0, Some(9)),
        (2, None),
        (1, None),
    ]);
    assert_snapshot!(g.debug_dot_with_settings(&DotExportSettings {
            include_autogenerated_fields: true,
            ..DotExportSettings::default()
        }),@r#"
        digraph GL1 {
          num = "1";
          overall_factor = "(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)";
          overall_factor_evaluated = "-1";
          projector = "ϵ(0,spenso::mink(4,hedge(1)))*ϵbar(0,spenso::mink(4,hedge(0)))";
          threshold_counterterms = "
        schema_version = 1
        cuts = []
        ";
          "0"[dod="0" dod_autogen="true" int_id="V_117" name_autogen="true" num="(1/2*spenso::g(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(4)))+1/2*spenso::gamma5(spenso::bis(4,hedge(4)),spenso::bis(4,hedge(2))))*UFO::GC_22*spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(4))))" num_autogen="true"];
          "1"[dod="0" dod_autogen="true" int_id="V_82" name_autogen="true" num="(-1/2*spenso::gamma5(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(8)))+1/2*spenso::g(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(8))))*UFO::GC_16*spenso::g(spenso::cof(3,hedge(8)),spenso::dind(spenso::cof(3,hedge(3))))" num_autogen="true"];
          "2"[dod="0" dod_autogen="true" int_id="V_71" name_autogen="true" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(5)),spenso::dind(spenso::cof(3,hedge(9))))*spenso::gamma(spenso::bis(4,hedge(9)),spenso::bis(4,hedge(5)),spenso::mink(4,hedge(0)))" num_autogen="true"];
          "3"[dod="1" dod_autogen="true" int_id="V_11" name_autogen="true" num="(-1*Q(4,spenso::mink(4,hedge(1)))+Q(1,spenso::mink(4,hedge(1))))*UFO::GC_3" num_autogen="true"];

          "2":1	-> "3":0	 [id=0 dir=none  dod="-2" dod_autogen="true" is_cut="0" lmb_rep="P(0,a___)" name="e0" name_autogen="true" num="1" num_autogen="true" particle="a"];
          "1":10	-> "3":11	 [id=1 dir=none  dod="-2" dod_autogen="true" lmb_rep="-1*K(1,a___)+-1*P(0,a___)" name="e1" name_autogen="true" num="𝑖" num_autogen="true" particle="G-"];
          "0":2	-> "1":3	 [id=2 dir=back  dod="-1" dod_autogen="true" lmb_id="0" lmb_rep="K(0,a___)" name="e2" name_autogen="true" num="(-1*Q(2,spenso::mink(4,edge(2,1)))*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(3)),spenso::mink(4,edge(2,1)))+UFO::MC*spenso::g(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(3))))*spenso::g(spenso::cof(3,hedge(3)),spenso::dind(spenso::cof(3,hedge(2))))*𝑖" num_autogen="true" particle="c~"];
          "0":4	-> "2":5	 [id=3  dod="-1" dod_autogen="true" lmb_rep="-1*K(0,a___)+-1*K(1,a___)" name="e3" name_autogen="true" num="Q(3,spenso::mink(4,edge(3,1)))*spenso::g(spenso::cof(3,hedge(4)),spenso::dind(spenso::cof(3,hedge(5))))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(4)),spenso::mink(4,edge(3,1)))*𝑖" num_autogen="true" particle="d"];
          "0":6	-> "3":7	 [id=4 dir=none  dod="-2" dod_autogen="true" lmb_id="1" lmb_rep="K(1,a___)" name="e4" name_autogen="true" num="𝑖" num_autogen="true" particle="G+"];
          "1":8	-> "2":9	 [id=5 dir=back  dod="-1" dod_autogen="true" lmb_rep="K(0,a___)+K(1,a___)+P(0,a___)" name="e5" name_autogen="true" num="-𝑖*Q(5,spenso::mink(4,edge(5,1)))*spenso::g(spenso::cof(3,hedge(9)),spenso::dind(spenso::cof(3,hedge(8))))*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(9)),spenso::mink(4,edge(5,1)))" num_autogen="true" particle="d~"];
        }
        "#);
}

#[test]
fn validate_verterx() {
    test_initialise().unwrap();
    let mut model = load_generic_model("sm");
    let mut restriction =
        InputParamCard::from_str(LIGHT_QUARKS_MASSLESS_RESTRICTION.into(), "json").unwrap();
    model.simplify(&mut restriction).unwrap();
    let g: Graph = dot!(digraph GL1{
      //  0[dod=0 int_id=V_82 num="UFO::GC_16*spenso::g(spenso::dind(spenso::cof(3,hedge(2))),spenso::dind(spenso::cof(3,hedge(4))))*spenso::projp(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(4)))"];
      // 1[dod=0 int_id=V_117 num="UFO::GC_22*spenso::g(spenso::dind(spenso::cof(3,hedge(8))),spenso::cof(3,hedge(3)))*spenso::projm(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(3)))"];
      // 2[dod=0 int_id=V_71 num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(5)),spenso::cof(3,hedge(9)))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(9)),spenso::mink(4,hedge(0)))"];
      // 3[dod=0 int_id=V_11 num="UFO::GC_3*(UFO::P(spenso::mink(4,hedge(1)),2)-UFO::P(spenso::mink(4,hedge(1)),3))"];

      2:1   -> 3:0   [id=0 is_cut=0 particle="a"];
      1:8   -> 2:9   [id=1  particle="d"];
      0:2   -> 1:3   [id=2 particle="c"];
      0:4   -> 2:5   [id=3  particle="d~"];
      0:6   -> 3:7   [id=4  particle="G-"];
      1:10  -> 3:11  [id=5 particle="G+"];
    }, &model)
    .unwrap();

    g.assert_fixture_hedge_routing(&[
        (0, Some(0)),
        (2, Some(1)),
        (0, Some(2)),
        (1, Some(3)),
        (1, Some(4)),
        (0, Some(5)),
        (2, None),
        (1, None),
        (0, Some(8)),
        (1, Some(9)),
        (2, None),
        (2, None),
    ]);
    assert_snapshot!(g.debug_dot_with_settings(&DotExportSettings {
            include_autogenerated_fields: true,
            ..DotExportSettings::default()
        }),@r#"
        digraph GL1 {
          num = "1";
          overall_factor = "1";
          overall_factor_evaluated = "1";
          projector = "ϵ(0,spenso::mink(4,hedge(0)))*ϵbar(0,spenso::mink(4,hedge(1)))";
          threshold_counterterms = "
        schema_version = 1
        cuts = []
        ";
          "0"[dod="0" dod_autogen="true" int_id="V_82" name_autogen="true" num="(-1/2*spenso::gamma5(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(4)))+1/2*spenso::g(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(4))))*UFO::GC_16*spenso::g(spenso::cof(3,hedge(4)),spenso::dind(spenso::cof(3,hedge(2))))" num_autogen="true"];
          "1"[dod="0" dod_autogen="true" int_id="V_117" name_autogen="true" num="(1/2*spenso::g(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(8)))+1/2*spenso::gamma5(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(3))))*UFO::GC_22*spenso::g(spenso::cof(3,hedge(3)),spenso::dind(spenso::cof(3,hedge(8))))" num_autogen="true"];
          "2"[dod="0" dod_autogen="true" int_id="V_71" name_autogen="true" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(9)),spenso::dind(spenso::cof(3,hedge(5))))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(9)),spenso::mink(4,hedge(1)))" num_autogen="true"];
          "3"[dod="1" dod_autogen="true" int_id="V_11" name_autogen="true" num="(-1*Q(5,spenso::mink(4,hedge(0)))+Q(4,spenso::mink(4,hedge(0))))*UFO::GC_3" num_autogen="true"];

          "2":1	-> "3":0	 [id=0 dir=none  dod="-2" dod_autogen="true" is_cut="0" lmb_rep="P(0,a___)" name="e0" name_autogen="true" num="1" num_autogen="true" particle="a"];
          "1":8	-> "2":9	 [id=1  dod="-1" dod_autogen="true" lmb_rep="K(0,a___)+K(1,a___)+P(0,a___)" name="e1" name_autogen="true" num="Q(1,spenso::mink(4,edge(1,1)))*spenso::g(spenso::cof(3,hedge(8)),spenso::dind(spenso::cof(3,hedge(9))))*spenso::gamma(spenso::bis(4,hedge(9)),spenso::bis(4,hedge(8)),spenso::mink(4,edge(1,1)))*𝑖" num_autogen="true" particle="d"];
          "0":2	-> "1":3	 [id=2  dod="-1" dod_autogen="true" lmb_id="0" lmb_rep="K(0,a___)" name="e2" name_autogen="true" num="(Q(2,spenso::mink(4,edge(2,1)))*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(2)),spenso::mink(4,edge(2,1)))+UFO::MC*spenso::g(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(3))))*spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(3))))*𝑖" num_autogen="true" particle="c"];
          "0":4	-> "2":5	 [id=3 dir=back  dod="-1" dod_autogen="true" lmb_rep="-1*K(0,a___)+-1*K(1,a___)" name="e3" name_autogen="true" num="-𝑖*Q(3,spenso::mink(4,edge(3,1)))*spenso::g(spenso::cof(3,hedge(5)),spenso::dind(spenso::cof(3,hedge(4))))*spenso::gamma(spenso::bis(4,hedge(4)),spenso::bis(4,hedge(5)),spenso::mink(4,edge(3,1)))" num_autogen="true" particle="d~"];
          "0":6	-> "3":7	 [id=4 dir=none  dod="-2" dod_autogen="true" lmb_id="1" lmb_rep="K(1,a___)" name="e4" name_autogen="true" num="𝑖" num_autogen="true" particle="G-"];
          "1":10	-> "3":11	 [id=5 dir=none  dod="-2" dod_autogen="true" lmb_rep="-1*K(1,a___)+-1*P(0,a___)" name="e5" name_autogen="true" num="𝑖" num_autogen="true" particle="G+"];
        }
        "#);
}

#[test]
fn test_ufo() {
    test_initialise().unwrap();
    let mut model = load_generic_model("sm");
    let mut restriction =
        InputParamCard::from_str(LIGHT_QUARKS_MASSLESS_RESTRICTION.into(), "json").unwrap();
    model.simplify(&mut restriction).unwrap();
    let g: Graph = dot!(
            digraph GL0{
                num = "1";
                overall_factor = "AutG(1)^-1*InternalFermionLoopSign(-1)*ExternalFermionOrderingSign(1)*AntiFermionSpinSumSign(1)";
                0[int_id=V_79];
            1[int_id=V_79];
            2[int_id=V_71];
            3[int_id=V_71];
            ext0 [style=invis];
            2:0-> ext0 [id=0 dir=back is_cut=0  particle="a"];
            ext1 [style=invis];
            ext1-> 3:1 [id=1 is_cut=0  particle="a"];
            0:2-> 1:3 [id=2   particle="d"];
            0:4-> 1:5 [id=3   particle="Z"];
            0:6-> 3:7 [id=4   particle="d~"];
            1:8-> 2:9 [id=5   particle="d"];
            2:10-> 3:11 [id=6   particle="d"];
            }, &model
        ).unwrap();

    g.assert_fixture_hedge_routing(&[
        (2, Some(1)),
        (2, Some(0)),
        (0, Some(2)),
        (1, Some(3)),
        (2, Some(4)),
        (2, Some(5)),
        (1, Some(6)),
        (0, Some(7)),
        (0, Some(8)),
        (1, Some(9)),
        (0, Some(10)),
        (1, Some(11)),
    ]);
    assert_snapshot!(g.debug_dot_with_settings(&DotExportSettings {
            include_autogenerated_fields: true,
            ..DotExportSettings::default()
        }),@r#"
        digraph GL0 {
          num = "1";
          overall_factor = "(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)";
          overall_factor_evaluated = "-1";
          projector = "ϵ(0,spenso::mink(4,hedge(1)))*ϵbar(0,spenso::mink(4,hedge(0)))";
          threshold_counterterms = "
        schema_version = 1
        cuts = []
        ";
          "0"[dod="0" dod_autogen="true" int_id="V_79" name_autogen="true" num="((-1/2*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(6)),spenso::mink(4,hedge(4)))+-3/2*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::gamma5(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(6))))*UFO::GC_58+(-1/2*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::gamma5(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(6)))+1/2*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(6)),spenso::mink(4,hedge(4))))*UFO::GC_50)*spenso::g(spenso::cof(3,hedge(6)),spenso::dind(spenso::cof(3,hedge(2))))" num_autogen="true"];
          "1"[dod="0" dod_autogen="true" int_id="V_79" name_autogen="true" num="((-1/2*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(3)),spenso::mink(4,hedge(5)))+-3/2*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::gamma5(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3))))*UFO::GC_58+(-1/2*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::gamma5(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3)))+1/2*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(3)),spenso::mink(4,hedge(5))))*UFO::GC_50)*spenso::g(spenso::cof(3,hedge(3)),spenso::dind(spenso::cof(3,hedge(8))))" num_autogen="true"];
          "2"[dod="0" dod_autogen="true" int_id="V_71" name_autogen="true" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(9)),spenso::dind(spenso::cof(3,hedge(10))))*spenso::gamma(spenso::bis(4,hedge(10)),spenso::bis(4,hedge(9)),spenso::mink(4,hedge(0)))" num_autogen="true"];
          "3"[dod="0" dod_autogen="true" int_id="V_71" name_autogen="true" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(11)),spenso::dind(spenso::cof(3,hedge(7))))*spenso::gamma(spenso::bis(4,hedge(7)),spenso::bis(4,hedge(11)),spenso::mink(4,hedge(1)))" num_autogen="true"];

          "2":1	-> "3":0	 [id=0 dir=none  dod="-2" dod_autogen="true" is_cut="0" lmb_rep="P(0,a___)" name="e0" name_autogen="true" num="1" num_autogen="true" particle="a"];
          "2":10	-> "3":11	 [id=1  dod="-1" dod_autogen="true" lmb_rep="-1*K(1,a___)+-1*P(0,a___)" name="e1" name_autogen="true" num="Q(1,spenso::mink(4,edge(1,1)))*spenso::g(spenso::cof(3,hedge(10)),spenso::dind(spenso::cof(3,hedge(11))))*spenso::gamma(spenso::bis(4,hedge(11)),spenso::bis(4,hedge(10)),spenso::mink(4,edge(1,1)))*𝑖" num_autogen="true" particle="d"];
          "0":2	-> "1":3	 [id=2  dod="-1" dod_autogen="true" lmb_id="0" lmb_rep="K(0,a___)" name="e2" name_autogen="true" num="Q(2,spenso::mink(4,edge(2,1)))*spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(3))))*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(2)),spenso::mink(4,edge(2,1)))*𝑖" num_autogen="true" particle="d"];
          "0":4	-> "1":5	 [id=3 dir=none  dod="-2" dod_autogen="true" lmb_rep="-1*K(0,a___)+-1*K(1,a___)" name="e3" name_autogen="true" num="-𝑖*spenso::g(spenso::mink(4,hedge(4)),spenso::mink(4,hedge(5)))" num_autogen="true" particle="Z"];
          "0":6	-> "3":7	 [id=4 dir=back  dod="-1" dod_autogen="true" lmb_id="1" lmb_rep="K(1,a___)" name="e4" name_autogen="true" num="-𝑖*Q(4,spenso::mink(4,edge(4,1)))*spenso::g(spenso::cof(3,hedge(7)),spenso::dind(spenso::cof(3,hedge(6))))*spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,hedge(7)),spenso::mink(4,edge(4,1)))" num_autogen="true" particle="d~"];
          "1":8	-> "2":9	 [id=5  dod="-1" dod_autogen="true" lmb_rep="-1*K(1,a___)" name="e5" name_autogen="true" num="Q(5,spenso::mink(4,edge(5,1)))*spenso::g(spenso::cof(3,hedge(8)),spenso::dind(spenso::cof(3,hedge(9))))*spenso::gamma(spenso::bis(4,hedge(9)),spenso::bis(4,hedge(8)),spenso::mink(4,edge(5,1)))*𝑖" num_autogen="true" particle="d"];
        }
        "#);
}

#[test]
fn parse_lmbsetting() {
    test_initialise().unwrap();
    let g: Graph = dot!(
        digraph G{
            graph [
                multiplicity_factor = "1/2"
                overall_factor = "2*x"
                color = "spenso::t"
                colorless = "spenso::gamma"
            ]
            ext [style=invis]
            A [num="1" color_num="a"]
            B [num="1"]
            C [num="1"]
            D [num="1"]
            E [num="1"]
            ext -> A  [pdg=1000]
            ext -> C  [pdg=1000]
            C -> A    [pdg=1000, num="P(eid,spenso::mink(4,0))"]
            A -> D    [pdg=1000, num="P(eid,spenso::mink(4,0))"]
            D -> B    [pdg=1000, num="1"]
            B -> ext   [pdg=1000]
            E -> B    [pdg=1000,lmb_id=0]
            E -> A    [pdg=1000,lmb_id=1]
            E -> C   [pdg=1000]
            C -> D    [pdg=1000, num="1"]
            B -> ext   [pdg=1000]
        },"scalars"
    )
    .unwrap();

    g.assert_fixture_hedge_routing(&[
        (0, None),
        (0, None),
        (0, None),
        (1, None),
        (0, None),
        (1, None),
        (1, None),
        (1, None),
        (2, None),
        (2, None),
        (0, None),
        (2, None),
        (1, None),
        (3, None),
        (2, None),
        (2, None),
        (3, None),
        (3, None),
    ]);
    assert_snapshot!(g.debug_dot_with_settings(&DotExportSettings {
            include_autogenerated_fields: true,
            ..DotExportSettings::default()
        }),@r#"
        digraph G {
          num = "1";
          overall_factor = "2*x";
          overall_factor_evaluated = "2*x";
          projector = "1";
          threshold_counterterms = "
        schema_version = 1
        cuts = []
        ";
          A[dod="0" dod_autogen="true" num="1"];
          B[dod="0" dod_autogen="true" num="1"];
          C[dod="0" dod_autogen="true" num="1"];
          D[dod="0" dod_autogen="true" num="1"];
          E[dod="0" dod_autogen="true" num="1"];

          A:0	-> D:1	 [id=0 dir=none  dod="-2" dod_autogen="true" lmb_rep="K(1,a___)+K(2,a___)+P(2,a___)" name="e0" name_autogen="true" num="P(0,spenso::mink(4,0))" particle="scalar_0"];
          ext1	 [style=invis];
          B:2	-> ext1	 [id=1 dir=none dod="-2" dod_autogen="true" lmb_rep="P(0,a___)" name="e1" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
          ext2	 [style=invis];
          B:3	-> ext2	 [id=2 dir=none dod="-2" dod_autogen="true" lmb_rep="P(1,a___)" name="e2" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
          C:4	-> A:5	 [id=3 dir=none  dod="-2" dod_autogen="true" lmb_id="2" lmb_rep="K(2,a___)" name="e3" name_autogen="true" num="P(3,spenso::mink(4,0))" particle="scalar_0"];
          C:6	-> D:7	 [id=4 dir=none  dod="-2" dod_autogen="true" lmb_rep="-1*K(0,a___)+-1*K(1,a___)+-1*K(2,a___)+-1*P(2,a___)+P(0,a___)+P(1,a___)" name="e4" name_autogen="true" num="1" particle="scalar_0"];
          D:8	-> B:9	 [id=5 dir=none  dod="-2" dod_autogen="true" lmb_rep="-1*K(0,a___)+P(0,a___)+P(1,a___)" name="e5" name_autogen="true" num="1" particle="scalar_0"];
          E:10	-> A:11	 [id=6 dir=none  dod="-2" dod_autogen="true" lmb_id="1" lmb_rep="K(1,a___)" name="e6" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
          E:12	-> B:13	 [id=7 dir=none  dod="-2" dod_autogen="true" lmb_id="0" lmb_rep="K(0,a___)" name="e7" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
          E:14	-> C:15	 [id=8 dir=none  dod="-2" dod_autogen="true" lmb_rep="-1*K(0,a___)+-1*K(1,a___)" name="e8" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
          ext9	 [style=invis];
          ext9	-> A:16	 [id=9 dir=none dod="-2" dod_autogen="true" lmb_rep="P(2,a___)" name="e9" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
          ext10	 [style=invis];
          ext10	-> C:17	 [id=10 dir=none dod="-2" dod_autogen="true" lmb_rep="-1*P(2,a___)+P(0,a___)+P(1,a___)" name="e10" name_autogen="true" num="𝑖" num_autogen="true" particle="scalar_0"];
        }
        "#);

    println!(
        "{}",
        g.underlying
            .dot_lmb_of(&g.underlying.full(), &g.loop_momentum_basis)
    );

    let num = Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter());

    let expr = num.state.expr.as_view();

    let lib: DummyLibrary<SymbolicTensor<Aind>> = DummyLibrary::<_>::new();
    let net =
        Network::<NetworkStore<SymbolicTensor<Aind>, Atom>, _, Symbol, Aind>::try_from_view::<
            SymbolicTensor<Aind>,
            _,
        >(expr, &lib, &ParseSettings::default())
        .unwrap();

    println!("{}", expr);
    println!(
        "{}",
        net.dot_display_impl(
            ToString::to_string,
            |_| None,
            |a| {
                if let Ok(a) = ShadowedStructure::<Aind>::parse(a.expression.as_view()) {
                    a.name()
                        .map(|s| {
                            if let Some(a) = a.args() {
                                FunctionBuilder::new(s).add_args(&a).finish().to_string()
                            } else {
                                s.to_string()
                            }
                        })
                        .unwrap_or("".to_string())
                } else {
                    "".to_string()
                }
            },
            ToString::to_string
        )
    );

    let num = num.color_simplify();

    println!("{}", num.state.expr);
    let num = num.gamma_simplify();

    println!("{}", num.state.expr);
}

#[test]
fn massive_gluon() {
    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("sm");
    let mut graph = Graph::from_dot(
        linnet::dot!(
        digraph ddxaaapentagon{
            // num = "-2";
            ext   [style=invis];
            ext	-> 0:0	 [id=0 particle="d"];
            ext	-> 1:1	 [id=1 particle="d~"];
            2:2	-> ext	 [id=2 particle="a"];
            3:3	-> ext	 [id=3 particle="a"];
            4:4	-> ext   [id=4 particle="a"];

            0	-> 2	 [id=5 particle="d"];
            2	-> 3	 [id=6 particle="d"];
            3	-> 4	 [id=7 particle="d"];
            4 -> 1     [id=8 particle="d"];
            // 1 -> 0     [id=9 lmb_id="0" particle="g"];
            // The option below would introduce an IR regulator
            1 -> 0     [id=9 lmb_id="0" particle="g" mass="1"];
        })
        .unwrap(),
        &model,
    )
    .unwrap();
    graph.global_prefactor.num *=
        graph.underlying[Hedge(0)].color_kronekers(&graph.underlying[Hedge(1)]);
    graph.validate_full_numerator_tensor_network().unwrap();
}

mod failing {
    use super::*;

    #[test]
    fn polarizations() {
        test_initialise().unwrap();
        let gs: Vec<Graph> = dot!(
            digraph g1{
               ext [style=invis];
                n [num=1]
                ext -> n [particle="d"];
            }
            digraph g2{
               ext [style=invis];
                n [num=1]
                n -> ext [particle="d"];
            }
            digraph g3{
               ext [style=invis];
                n [num=1]
                ext -> n [dir=back particle="d~"];
            }
            digraph g4{
               ext [style=invis];
                n [num=1]
                n -> ext [dir=back particle="d~"];
            }

            digraph g3p{
               ext [style=invis];
                n [num=1]
                ext -> n [ particle="d~"];
            }
            digraph g4p{
               ext [style=invis];
                n [num=1]
                n -> ext [particle="d~"];
            }
        )
        .unwrap();

        insta::assert_snapshot!(gs[0].global_prefactor.projector.to_bare_ordered_string(),@"u(0,bis(4,hedge(0)))");
        insta::assert_snapshot!(gs[1].global_prefactor.projector.to_bare_ordered_string(),@"ubar(0,bis(4,hedge(0)))");
        insta::assert_snapshot!(gs[2].global_prefactor.projector.to_bare_ordered_string(),@"vbar(0,bis(4,hedge(0)))");
        insta::assert_snapshot!(gs[3].global_prefactor.projector.to_bare_ordered_string(),@"v(0,bis(4,hedge(0)))");
        insta::assert_snapshot!(gs[4].global_prefactor.projector.to_bare_ordered_string(),@"vbar(0,bis(4,hedge(0)))");
        insta::assert_snapshot!(gs[5].global_prefactor.projector.to_bare_ordered_string(),@"v(0,bis(4,hedge(0)))");
    }

    #[test]
    fn vertex_rules() {
        use idenso::IndexTooling;
        use itertools::Itertools;
        use symbolica::{atom::AtomCore, parse};

        test_initialise().unwrap();
        let gs: Vec<Graph> = dot!(
            digraph g1{
               ext [style=invis];
                ext -> n [particle="g"];
                ext -> n [particle="g"];
                ext -> n [particle="g"];
            }

            digraph g2{
               ext [style=invis];
                ext -> n [particle="d"];
                ext -> n [particle="d~"];
                ext -> n [particle="g"];
            }

            digraph g3{
               ext [style=invis];
                ext -> n [id = 0 particle="ghG"];
                ext -> n [id = 1 particle="ghG~"];
                ext -> n [id = 2 particle="g"];
            }

            digraph g4{
               ext [style=invis];
                ext -> n [id = 0 particle="ghG"];
                n -> ext [id = 1 particle="ghG"];
                ext -> n [id = 2 particle="g"];
            }

            digraph GL55{
              ext   [style=invis];
              0:0   -> ext   [id=0 dir=none particle="g"];
              0:1   -> ext   [id=1 dir=none particle="ghG~"];
              0:2   -> ext   [id=2 dir=none particle="ghG"];
            }

            digraph GL54{
              ext   [style=invis];
              0:0   -> ext  [id=0 dir=none  particle="g"];
              0:1   -> ext  [id=1 dir=none  particle="ghG"];
              0:2   -> ext  [id=2 dir=none  particle="ghG~"];
            }

        )
        .unwrap();

        // Compare complete vertex tensors, preserving physical hedge/EMR labels while
        // allowing equivalent algebra and contracted dummy names to vary.
        let expected = [
            parse!(
                "\
                (-1*Q(0,spenso::mink(4,hedge(1)))*spenso::g(spenso::mink(4,hedge(0)),spenso::mink(4,hedge(2)))+
                -1*Q(1,spenso::mink(4,hedge(2)))*spenso::g(spenso::mink(4,hedge(0)),spenso::mink(4,hedge(1)))+
                -1*Q(2,spenso::mink(4,hedge(0)))*spenso::g(spenso::mink(4,hedge(1)),spenso::mink(4,hedge(2)))+
                Q(0,spenso::mink(4,hedge(2)))*spenso::g(spenso::mink(4,hedge(0)),spenso::mink(4,hedge(1)))+
                Q(1,spenso::mink(4,hedge(0)))*spenso::g(spenso::mink(4,hedge(1)),spenso::mink(4,hedge(2)))+
                Q(2,spenso::mink(4,hedge(1)))*spenso::g(spenso::mink(4,hedge(0)),spenso::mink(4,hedge(2))))*UFO::GC_10*spenso::f(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)),spenso::coad(8,hedge(2)))",
                default_namespace = "gammalooprs"
            ),
            parse!(
                "\
                UFO::GC_11*spenso::gamma(spenso::bis(4,hedge(1)),spenso::bis(4,hedge(0)),spenso::mink(4,hedge(2)))*spenso::t(spenso::coad(8,hedge(2)),spenso::cof(3,hedge(0)),spenso::dind(spenso::cof(3,hedge(1))))",
                default_namespace = "gammalooprs"
            ),
            // The corrected UUV1 rule reverses these four ghost-vector orientations.
            parse!(
                "\
                UFO::GC_10*Q(1,spenso::mink(4,hedge(2)))*spenso::f(spenso::coad(8,hedge(2)),spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))",
                default_namespace = "gammalooprs"
            ),
            parse!(
                "\
                -1*UFO::GC_10*Q(1,spenso::mink(4,hedge(1)))*spenso::f(spenso::coad(8,hedge(1)),spenso::coad(8,hedge(0)),spenso::coad(8,hedge(2)))",
                default_namespace = "gammalooprs"
            ),
            parse!(
                "\
                -1*UFO::GC_10*Q(2,spenso::mink(4,hedge(0)))*spenso::f(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)),spenso::coad(8,hedge(2)))",
                default_namespace = "gammalooprs"
            ),
            parse!(
                "\
                -1*UFO::GC_10*Q(1,spenso::mink(4,hedge(0)))*spenso::f(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(2)),spenso::coad(8,hedge(1)))",
                default_namespace = "gammalooprs"
            ),
        ];
        assert_eq!(
            gs.len(),
            expected.len(),
            "every input vertex needs an oracle"
        );
        for (graph, expected) in gs.iter().zip(expected) {
            let actual = graph
                .underlying
                .iter_nodes()
                .map(|(_, _, vertex)| vertex.get_num())
                .exactly_one()
                .unwrap_or_else(|_| panic!("input contains one physical vertex"))
                .canonize(Aind::Dummy)
                .unwrap();
            let expected = expected.canonize(Aind::Dummy).unwrap();
            assert!(
                (actual - expected).expand().is_zero(),
                "vertex rule for {} changed",
                graph.name,
            );
        }

        let gs: Vec<Graph> = dot!(
            digraph g1{
               ext [style=invis];
                ext -> n [particle="scalar_0"];
                ext -> n [particle="scalar_0"];
                ext -> n [particle="graviton"];
            },"scalar_gravity"
        )
        .unwrap();
        let expected = parse!(
            "\
            (-1*Q(0,spenso::mink(4,vertex(0,1)))*Q(1,spenso::mink(4,vertex(0,1)))*spenso::g(spenso::mink(4,hedge(2)),spenso::mink(4,hedge(2,1)))+
            Q(0,spenso::mink(4,hedge(2)))*Q(1,spenso::mink(4,hedge(2,1)))+
            Q(0,spenso::mink(4,hedge(2,1)))*Q(1,spenso::mink(4,hedge(2))))*UFO::SST+
            -1*UFO::SSTmpart0*spenso::g(spenso::mink(4,hedge(2)),spenso::mink(4,hedge(2,1)))",
            default_namespace = "gammalooprs"
        );
        assert_eq!(gs.len(), 1, "the scalar-gravity input needs one oracle");
        let actual = gs[0]
            .underlying
            .iter_nodes()
            .map(|(_, _, vertex)| vertex.get_num())
            .exactly_one()
            .unwrap_or_else(|_| panic!("input contains one physical vertex"))
            .canonize(Aind::Dummy)
            .unwrap();
        let expected = expected.canonize(Aind::Dummy).unwrap();
        assert!(
            (actual - expected).expand().is_zero(),
            "scalar-gravity vertex rule changed",
        );
    }

    #[test]
    fn dot_export_can_include_autogenerated_fields() {
        test_initialise().unwrap();
        let g: Graph = dot!(
            digraph GL1{
                num = "1";
                overall_factor = "AutG(1)^-1*InternalFermionLoopSign(-1)*ExternalFermionOrderingSign(1)*AntiFermionSpinSumSign(1)";
                0[int_id=V_117];
                1[int_id=V_82];
                2[int_id=V_71];
                3[int_id=V_11];

                ext0 [style=invis];
                2:0 -> ext0 [id=0 dir=back is_cut=0 particle="a"];
                ext1 [style=invis];
                ext1 -> 3:1 [id=1 is_cut=0 particle="a"];
                0:2 -> 1:3 [id=2 particle="c~"];
                0:4 -> 2:5 [id=3 particle="d"];
                0:6 -> 3:7 [id=4 particle="G+"];
                1:8 -> 2:9 [id=5 particle="d~"];
                1:10 -> 3:11 [id=6 particle="G-"];
            }
        )
        .unwrap();

        assert_eq!(g.underlying[NodeIndex(3)].dod.value, 1);
        assert_eq!(g.underlying[EdgeIndex(2)].dod.value, -1);
        assert_eq!(g.underlying[EdgeIndex(4)].dod.value, -2);

        let default_dot = g.debug_dot();
        assert!(!default_dot.contains("dod_autogen"));
        assert!(!default_dot.contains("num_autogen"));
        assert!(!default_dot.contains("dod=\"1\""));

        let with_autogenerated_fields = g.debug_dot_with_settings(&DotExportSettings {
            include_autogenerated_fields: true,
            ..DotExportSettings::default()
        });

        assert!(with_autogenerated_fields.contains("dod=\"1\""));
        assert!(with_autogenerated_fields.contains("dod_autogen"));
        assert!(with_autogenerated_fields.contains("num_autogen"));
        assert!(with_autogenerated_fields.contains("name_autogen"));
    }

    #[test]
    fn test_pols_tree() {
        test_initialise().unwrap();
        let g: Graph = dot!(digraph qqx_aaa_tree_1 {
                    ext    [style=invis]
                    ext -> v1:1 [particle="d" id=1];
                    ext -> v3:2 [dir=back particle="d~" id=2];
                    v1:3 -> ext [particle="a" id=3];
                    v2:4 -> ext [particle="a" id=4];
                    v3:0 -> ext [particle="a" id=0];
                    v1:5 -> v2:6 [particle="d" id=5];
                    v2 -> v3 [particle="d" id=6];
        })
        .unwrap();

        assert_snapshot!(g.debug_dot(),@r#"
        digraph qqx_aaa_tree_1{
          num = "1";
          overall_factor = "1";
          overall_factor_evaluated = "1";
          projector = "u(1,spenso::bis(4,hedge(1)))*vbar(2,spenso::bis(4,hedge(2)))*ϵbar(0,spenso::mink(4,hedge(0)))*ϵbar(3,spenso::mink(4,hedge(3)))*ϵbar(4,spenso::mink(4,hedge(4)))";
          v1[dod="0" int_id="V_71" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(1)),spenso::dind(spenso::cof(3,hedge(5))))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(1)),spenso::mink(4,hedge(3)))"];
          v2[dod="0" int_id="V_71" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(6)),spenso::dind(spenso::cof(3,hedge(7))))*spenso::gamma(spenso::bis(4,hedge(7)),spenso::bis(4,hedge(6)),spenso::mink(4,hedge(4)))"];
          v3[dod="0" int_id="V_71" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(8)),spenso::dind(spenso::cof(3,hedge(2))))*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(8)),spenso::mink(4,hedge(0)))"];

          ext0	 [style=invis];
          v3:0	-> ext0	 [id=0 dir=none source="{ufo_order:2}" dod="-2"  lmb_rep="P(0,a___)" name="e0" num="1𝑖" particle="a"];
          ext1	 [style=invis];
          ext1	-> v1:1	 [id=1 sink="{ufo_order:1}" dod="-2"  lmb_rep="P(1,a___)" name="e1" num="1𝑖" particle="d"];
          ext2	 [style=invis];
          ext2	-> v3:2	 [id=2 dir=back sink="{ufo_order:0}" dod="-2"  lmb_rep="P(2,a___)" name="e2" num="1𝑖" particle="d~"];
          ext3	 [style=invis];
          v1:3	-> ext3	 [id=3 dir=none source="{ufo_order:2}" dod="-2"  lmb_rep="P(3,a___)" name="e3" num="1𝑖" particle="a"];
          ext4	 [style=invis];
          v2:4	-> ext4	 [id=4 dir=none source="{ufo_order:2}" dod="-2"  lmb_rep="-1*P(0,a___)+-1*P(3,a___)+P(1,a___)+P(2,a___)" name="e4" num="1𝑖" particle="a"];
          v1:5	-> v2:6	 [id=5 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_rep="-1*P(3,a___)+P(1,a___)" name="e5" num="Q(5,spenso::mink(4,edge(5,1)))*spenso::g(spenso::cof(3,hedge(5)),spenso::dind(spenso::cof(3,hedge(6))))*spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,hedge(5)),spenso::mink(4,edge(5,1)))" particle="d"];
          v2:7	-> v3:8	 [id=6 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_rep="-1*P(2,a___)+P(0,a___)" name="e6" num="Q(6,spenso::mink(4,edge(6,1)))*spenso::g(spenso::cof(3,hedge(7)),spenso::dind(spenso::cof(3,hedge(8))))*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(7)),spenso::mink(4,edge(6,1)))" particle="d"];
        }
        "#);
    }

    #[test]
    fn parse_local() {
        test_initialise().unwrap();
        let graphs: Vec<Graph> = dot!(
            digraph G{
                edge [num="1" dod=-1000 pdg=1000]
                node [num="1" dod=-1000]
                ext     [style=invis]
                A [id=1]
                B [id=0]
                C
                D [id=2]
                ext -> A  [id=1]
                ext -> C
                C -> A    [num="P(eid,spenso::mink(4,eid(0)))*Q(eid,spenso::mink(4,eid(0)))"]
                A -> D    [num="P(eid,spenso::mink(4,0))"]
                D:1 -> B:2    [id=2]
                B -> ext
                C -> D
                B -> ext
            }

            digraph G2{
                ext [style=invis]
                A [num="1"]
                B [num="1"]
                C [num="1"]
                D [num="1"]
                ext -> A  [pdg=1000]
                ext -> C  [pdg=1000]
                C -> A    [pdg=1000, num="P(eid,spenso::mink(4,eid(0)))*Q(eid,spenso::mink(4,eid(0)))"]
                A -> D    [pdg=1000, num="K(eid,spenso::mink(4,0))"]
                D -> B    [pdg=1000, num="1"]
                B -> ext   [pdg=1000]
                C -> D    [pdg=1000, num="1"]
                B -> ext   [pdg=1000]
            }
            ,
            "scalars"
        )
        .unwrap();

        let g = &graphs[1];

        println!(
            "{}",
            g.underlying
                .dot_lmb_of(&g.underlying.full(), &g.loop_momentum_basis)
        );

        let num = Numerator::<UnInit>::default().from_new_graph(g, &g.underlying.full_filter());

        let expr = num.state.expr.as_view();

        let lib: DummyLibrary<SymbolicTensor<Aind>> = DummyLibrary::<_>::new();
        let net =
            Network::<NetworkStore<SymbolicTensor<Aind>, Atom>, _, Symbol, Aind>::try_from_view::<
                SymbolicTensor<Aind>,
                _,
            >(expr, &lib, &ParseSettings::default())
            .unwrap();

        println!("{}", expr);
        println!(
            "{}",
            net.dot_display_impl(
                ToString::to_string,
                |_| None,
                |a| {
                    if let Ok(a) = ShadowedStructure::<Aind>::parse(a.expression.as_view()) {
                        a.name()
                            .map(|s| {
                                if let Some(a) = a.args() {
                                    FunctionBuilder::new(s).add_args(&a).finish().to_string()
                                } else {
                                    s.to_string()
                                }
                            })
                            .unwrap_or("".to_string())
                    } else {
                        "".to_string()
                    }
                },
                ToString::to_string
            )
        );

        let num = num.color_simplify();

        println!("{}", num.state.expr);
        // let num = num.gamma_simplify();

        // println!(
        //     "{}",
        //     num.state
        //         .colorless
        //         .get_owned_linear(FlatIndex::from(0))
        //         .unwrap()
        //         .to_dots()
        // );
    }

    #[test]
    fn explicit_hedge_payload_round_trips_in_dot_export() {
        test_initialise().unwrap();
        let g: Graph = dot!(
            digraph payload_graph{
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

    #[test]
    fn split_dot_export_places_initial_state_cut_inputs_left_and_outputs_right() {
        test_initialise().unwrap();
        let g: Graph = dot!(
            digraph cut_external_flow_graph{
                num = "1";
                overall_factor = "AutG(1)^-1";
                0[num=1 dod=0]
                1[num=1 dod=0]
                ext_in [style=invis]
                ext_out [style=invis]
                ext_in -> 0:0 [id=0 is_cut=0 name=e_in num=1 dod=-2 sink="{ufo_order:0,dod:-2}"]
                1:1 -> ext_out [id=1 is_cut=0 name=e_out num=1 dod=-2 source="{ufo_order:1,dod:-2}"]
                0:2 -> 1:3 [id=2 name=prop num=1 dod=-2]
            }
        )
        .unwrap();

        let split = g.dot_serialize(&DotExportSettings {
            split_xs_by_initial_states: true,
            ..DotExportSettings::default()
        });

        assert!(
            g.initial_state_cut
                .iter_edges_idx_relative()
                .next()
                .is_some()
        );

        let incoming_edge = split
            .lines()
            .find(|line| line.contains("sink=") && line.contains("is_cut="))
            .unwrap_or_else(|| panic!("missing incoming cut edge in split DOT:\n{split}"));
        assert!(
            incoming_edge.contains("pos=\"x:@-left!,y:@edge"),
            "{incoming_edge}"
        );
        assert!(incoming_edge.contains("!\""), "{incoming_edge}");

        let outgoing_edge = split
            .lines()
            .find(|line| line.contains("source=") && line.contains("is_cut="))
            .unwrap_or_else(|| panic!("missing outgoing cut edge in split DOT:\n{split}"));
        assert!(
            outgoing_edge.contains("pos=\"x:@+right!,y:@edge"),
            "{outgoing_edge}"
        );
        assert!(outgoing_edge.contains("!\""), "{outgoing_edge}");
    }
}

#[test]
fn missing_vertex_dod_is_autogenerated() {
    test_initialise().unwrap();
    let graph: Graph = dot!(
        digraph missing_vertex_dod{
            ext_in [style=invis]
            ext_out [style=invis]
            A [num=1]
            ext_in -> A [name=e_in num=1 dod=-2 sink="{ufo_order:0,dod:-2}"]
            A -> ext_out [name=e_out num=1 dod=-2 source="{ufo_order:1,dod:-2}"]
        }
    )
    .unwrap();

    let dod = &graph.underlying[NodeIndex(0)].dod;
    assert_eq!(dod.value, 0);
    assert!(dod.autogenerated);
}

#[test]
fn hedge_payload_does_not_require_dod() {
    test_initialise().unwrap();
    let graph: Graph = dot!(
        digraph missing_hedge_dod{
            ext_in [style=invis]
            ext_out [style=invis]
            A [num=1 dod=0]
            ext_in -> A [name=e_in num=1 dod=-2 sink="{ufo_order:0}"]
            A -> ext_out [name=e_out num=1 dod=-2 source="{ufo_order:1}"]
        }
    )
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
