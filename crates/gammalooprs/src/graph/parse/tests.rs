use color_eyre::Result;
use feynkit_graph::{DiagramEdge, EdgeEndpoints, EdgeId, VertexId};
use idenso::tensor::inference::TensorInferenceError;
use linnet::half_edge::NodeIndex;
use linnet::half_edge::involution::{EdgeIndex, HedgePair};
use linnet::half_edge::subgraph::SubSetLike;
use serde_json::Error as JsonError;
use serde_json::ser;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::abstract_index::{AIND_SYMBOLS, AbstractIndex};
use spenso::structure::representation::{Minkowski, RepName};
use spenso::structure::slot::IsAbstractSlot;
use spenso::tensor_symbol;
use symbolica::atom::{Atom, AtomCore, FunctionBuilder};
use symbolica::function;
use symbolica::parse;
use symbolica::symbol;
use typed_index_collections::ti_vec;

use super::{Graph, NumeratorInterface};
use crate::dot;
use crate::graph::parse::{IntoGraph, complete_group_parsing};
use crate::graph::{FeynmanGraph, GraphGroup};
use crate::initialisation::test_initialise;
use crate::momentum::sample::LoopIndex;
use crate::processes::DotExportSettings;
use crate::utils::{self, GS};
use crate::uv::UltravioletGraph;
use crate::uv::uv_graph::UVE;

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
fn numerator_interface_source_roundtrip_keeps_two_identical_open_axes() {
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
    test_initialise().unwrap();
    let model = utils::load_generic_model("scalars");
    let rep = Minkowski {}.new_rep(4).to_lib();
    let anonymous = AbstractIndex::Open {
        owner: 801,
        axis: 0,
    };
    let vector_head = tensor_symbol!("gammalooprs::runtime_anonymous_vector");
    for port in [
        rep.to_symbolic([]),
        rep.slot::<AbstractIndex, _>(anonymous).to_atom(),
        rep.slot::<AbstractIndex, _>(anonymous.scoped(symbol!("runtime_open_scope")))
            .to_atom(),
    ] {
        let vector = function!(vector_head, &port);
        let interface = NumeratorInterface(vec![port]);
        let metadata = ser::to_string(&interface.to_json()).unwrap();
        // A portable ranked zero is valid source data. GammaLoop requires an
        // explicit graph interface even when its expression vanishes.
        assert_eq!(
            NumeratorInterface::admit(Atom::Zero, Some(interface))
                .unwrap()
                .rank(),
            1
        );
        for expression in [Atom::Zero, vector] {
            let ranked_zero = expression.is_zero();
            let expression = ser::to_string(&expression.to_canonical_string()).unwrap();
            for source in [
                triangle_fixture(3, None).replace(
                    "graph [projector=\"1\"]",
                    &format!("graph [projector={expression} projector_interface={metadata}]"),
                ),
                triangle_fixture(3, None).replace(
                    "graph [projector=\"1\"]",
                    &format!(
                        "graph [projector=\"1\" num={expression} \
                         numerator_prefactor_interface={metadata}]"
                    ),
                ),
                triangle_fixture(3, None).replace(
                    "node [num=\"1\"]",
                    &format!("node [num={expression} num_interface={metadata}]"),
                ),
                triangle_fixture(3, None).replace(
                    "edge [pdg=1000 num=\"1\" dir=none]",
                    &format!("edge [pdg=1000 num={expression} num_interface={metadata} dir=none]"),
                ),
            ] {
                let result: Result<Graph> = source.into_graph(&model);
                let error = result
                    .err()
                    .expect("anonymous runtime ports must be rejected");
                let diagnostic = format!("{error:#}");
                let runtime_refusal = diagnostic.contains("requires an explicit index");
                // An authored open index in a nonzero payload can already fail
                // native interface transport. Preserve that typed source error;
                // ranked zeros must reach the runtime explicit-index boundary.
                let native_refusal = !ranked_zero
                    && error
                        .chain()
                        .any(|cause| cause.is::<TensorInferenceError>());
                assert!(runtime_refusal || native_refusal, "{error:?}");
            }
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
            ser::to_string(&scalar.to_canonical_string()).unwrap()
        ),
    );
    let graph: Graph = source.into_graph(&model).unwrap();
    assert!(graph.global_prefactor.num.is_scalar());
    assert_eq!(
        NumeratorInterface::export(&graph.global_prefactor.num).0,
        scalar
    );
}

#[test]
fn runtime_interface_validation_preserves_native_source_errors() {
    test_initialise().unwrap();
    let model = utils::load_generic_model("scalars");
    let interface = NumeratorInterface(vec![Minkowski {}.new_rep(4).to_symbolic([])]);
    let metadata = ser::to_string(&interface.to_json()).unwrap();
    let source = triangle_fixture(3, None)
        .replace(
            "graph [projector=\"1\"]",
            &format!("graph [projector=\"1\" num=\"0\" numerator_prefactor_interface={metadata}]"),
        )
        .replace(
            "edge [pdg=1000 num=\"1\" dir=none]",
            "edge [pdg=1000 num=\"1\" num_interface=\"not-json\" dir=none]",
        );
    // Source admission must finish before the runtime-only interface policy.
    let error = Graph::from_dot(&source, &model)
        .err()
        .expect("malformed edge interface metadata must fail native admission");
    assert!(
        error.chain().any(|cause| cause.is::<JsonError>()),
        "{error:?}"
    );
}

#[test]
fn runtime_and_canonical_dot_preserve_ranked_zero_fragments() {
    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("scalars");
    let canonical = Graph::from_dot(
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
        let restored = Graph::from_dot(&serialized, &model).unwrap().remove(0);
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
    use crate::numerator::GlobalPrefactor;
    use crate::uv::UltravioletGraph;
    use idenso::tensor::SymbolicTensor;
    use spenso::structure::IndexedStructure;
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
    let canonical = Graph::from_dot(
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
        let restored = Graph::from_dot(&dot, &model).unwrap().remove(0);
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
        assert_eq!(edges[1].extra_data.momtrop_edge_power, Some(parse!("1/2")));
        assert_eq!(
            graph.underlying[linnet::half_edge::NodeIndex(0)].dod.value,
            1
        );
        assert_eq!(graph.param_builder.pairs.additional_params.params.len(), 2);
    }
}

fn assert_native_runtime_parity(diagram: &feynkit_graph::FeynmanDiagram) {
    use crate::uv::UltravioletGraph;
    use linnet::half_edge::EdgeAccessors;
    use linnet::half_edge::involution::{Flow, Hedge};
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
        super::NumeratorInterface::export(diagram.numerator().unwrap())
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
        numerator: super::NumeratorInterface::admit(vertex_numerator, None).unwrap(),
    });
    if paired {
        builder
            .add_edge(vertex, vertex, DiagramEdge::new(particle, false))
            .unwrap();
    }
    for numerator in dangling {
        let mut edge = DiagramEdge::new(particle, false);
        edge.is_dummy = true;
        edge.numerator = super::NumeratorInterface::admit(numerator.clone(), None).unwrap();
        builder.add_edge(Some(vertex), None, edge).unwrap();
    }
    let diagram = builder.build().unwrap();
    diagram.validate().unwrap();
    diagram
}

#[test]
fn native_runtime_conversion_rejects_changed_frame_before_cut_indexing() {
    use linnet::half_edge::EdgeAccessors;
    use linnet::half_edge::involution::Orientation;
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
        numerator: super::NumeratorInterface::admit(Atom::one(), None).unwrap(),
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
                numerator: super::NumeratorInterface::admit(Atom::one(), None).unwrap(),
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
    use crate::numerator::aind::Aind;
    use crate::numerator::symbolica_ext::NumeratorAtomExt;
    use idenso::dirac::AGS;
    use idenso::representations::Bispinor;
    use spenso::network::library::symbolic::ETS;
    use spenso::structure::TensorStructure;
    use symbolica::function;

    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("sm");
    // The model's projector closes the external spinors; the numerator keeps them open.
    let source = r#"digraph port_self_energy {
        ext[style=invis]
        ext -> a[id=0 particle=t]
        a -> b[id=1 particle=t lmb_id=0]
        a -> b[id=2 particle=g]
        b -> ext[id=3 particle=t]
    }"#;
    let graph = Graph::from_dot(source, &model).unwrap().remove(0);
    let numerator = graph
        .numerator(&graph.full_filter(), &graph.empty_subgraph())
        .unwrap()
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
}

#[test]
fn compact_model_import_preserves_metadata_in_the_finalized_frame() {
    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("scalars");
    let mut parsed: linnet::parser::DotGraph = linnet::parser::DotGraph::from_string(
        r#"digraph mapped_bubble {
        final_state="scalar_1,scalar_1"
        group_id=0
        is_group_master=true
        params="source_weight;source_mass"
        ext [style=invis]
        ext -> a [id=0 particle=scalar_1 is_cut=7]
        a -> b [id=1 particle=scalar_1 lmb_id=0]
        b -> ext [id=2 particle=scalar_1 is_cut=7]
        a -> b [id=3 particle=scalar_1 mass="source_mass"]
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
    let imported = Graph::from_dot(&parsed.debug_dot(), &model)
        .unwrap()
        .remove(0);
    assert_eq!(imported.group_id, Some(crate::graph::GroupId(0)));
    assert!(imported.is_group_master);
    assert_eq!(
        imported.param_builder.pairs.additional_params.params.len(),
        2
    );
    let mass = symbol!("gammalooprs::source_mass").to_atom();
    assert_eq!(imported.underlying[EdgeIndex(2)].mass_atom(&model), mass);
    let energy = imported
        .param_builder
        .reps
        .iter()
        .find(|entry| entry.lhs == GS.ose(EdgeIndex(2)))
        .expect("the remapped loop edge has a registered energy");
    let spatial_square = Atom::add_many(
        (1..4).map(|axis| symbolica::function!(GS.emr_mom, 2, GS.cind(axis)).pow(2)),
    );
    assert_eq!(energy.rhs, (mass.pow(2) + spatial_square).pow((1, 2)));
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
    assert_eq!(topology.numerator().unwrap().expression(), &Atom::one());
    let reference = topology.apply_feynman_rules().unwrap();
    let imported = Graph::from_dot(source, &model).unwrap().remove(0);
    assert_ne!(reference.numerator().unwrap().expression(), &Atom::one());
    assert_eq!(
        imported
            .numerator(&imported.full_filter(), &imported.empty_subgraph())
            .unwrap()
            .get_single_atom()
            .unwrap(),
        super::NumeratorInterface::export(reference.numerator().unwrap()).0
    );
    assert_eq!(
        imported.global_prefactor.projector,
        reference
            .external_state_projector()
            .unwrap()
            .try_into_indexed()
            .unwrap()
    );
    assert!(!imported.global_prefactor.projector.expression().is_one());
    // Explicit factors override the model's. Import requires a contracted network,
    // so supply the model's projector and seven times the model's factor of `a`.
    let vertex = super::NumeratorInterface::export(
        &imported.underlying[linnet::half_edge::NodeIndex(0)]
            .num
            .value,
    )
    .0;
    let projector = super::NumeratorInterface::export(&imported.global_prefactor.projector).0;
    let supplied = source.replace(
        "ext [style=invis]",
        &format!(
            r#"projector="{}"; a[num="7*({})"]; ext [style=invis]"#,
            projector.to_canonical_string(),
            vertex.to_canonical_string()
        ),
    );
    let explicit = Graph::from_dot(&supplied, &model).unwrap().remove(0);
    assert!(
        explicit
            .global_prefactor
            .projector
            .same_value(&imported.global_prefactor.projector)
    );
    assert_eq!(
        super::NumeratorInterface::export(
            &explicit.underlying[linnet::half_edge::NodeIndex(0)]
                .num
                .value
        )
        .0,
        vertex * 7
    );
    assert!(
        !explicit.underlying[EdgeIndex(1)]
            .num
            .value
            .expression()
            .is_one()
    );
}

#[test]
fn compact_gamma_loop_indices_bind_to_native_rule_ports() {
    test_initialise().unwrap();
    let model = crate::utils::load_generic_model("sm");
    let source =
        include_str!("../../../../../tests/resources/graphs/local_ir_gluon_self_energy.dot");
    let qualified = Graph::from_dot(source, &model).unwrap().remove(0);
    let unqualified = Graph::from_dot(&source.replace("gammalooprs::", ""), &model)
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
    let projector = super::NumeratorInterface::export(&qualified.global_prefactor.projector).0;
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
    let mut source: linnet::parser::DotGraph = linnet::parser::DotGraph::from_string(
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
    let automatic = Graph::from_dot(&source.debug_dot(), &model)
        .unwrap()
        .remove(0);
    let mut projector = Atom::one();
    for (pair, edge, _) in source.iter_edges() {
        if let HedgePair::Unpaired { hedge, flow } = pair {
            let slot = Bispinor {}
                .new_rep(4)
                .to_symbolic([Atom::from(Aind::Hedge(hedge.0 as u16, 1))]);
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
    let explicit = Graph::from_dot(&source.debug_dot(), &model)
        .unwrap()
        .remove(0);
    // An admitted projector orders its ports canonically, the generated one by
    // construction; the wavefunctions and their slots must agree.
    assert_eq!(
        explicit.global_prefactor.projector.expression(),
        automatic.global_prefactor.projector.expression(),
        "supplied wavefunctions must follow actual source incidence through sewing"
    );
    assert_eq!(explicit.loop_momentum_basis, automatic.loop_momentum_basis);
    assert_eq!(
        explicit
            .numerator(&explicit.full_filter(), &explicit.empty_subgraph())
            .unwrap()
            .get_single_atom()
            .unwrap(),
        automatic
            .numerator(&automatic.full_filter(), &automatic.empty_subgraph())
            .unwrap()
            .get_single_atom()
            .unwrap()
    );
}

#[test]
fn full_numerator_requires_one_q3_structure_argument() {
    test_initialise().unwrap();
    let mut graph: Graph = crate::graph::Graph::from_dot(
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
    graph.overall_factor = Atom::one();
    graph.global_prefactor.num =
        super::NumeratorInterface::admit(GS.emr_vec(EdgeIndex(3), abstract_slot), None)
            .unwrap()
            .try_into_indexed()
            .unwrap();
    let error = graph.validate_full_numerator_tensor_network().unwrap_err();
    assert!(
        error
            .to_string()
            .contains("Full numerator still has dangling tensor indices"),
        "a numeric Minkowski index must be validated as an abstract index: {error:?}",
    );
    assert!(
        graph
            .global_prefactor
            .num
            .expression()
            .contains_symbol(GS.emr_vec)
            && !graph
                .global_prefactor
                .num
                .expression()
                .contains_symbol(GS.emr_mom),
        "a numeric Minkowski index must remain an abstract Q3 label, not normalize to a Q component",
    );

    graph.global_prefactor.num = super::NumeratorInterface::admit(Atom::one(), None)
        .unwrap()
        .try_into_indexed()
        .unwrap();

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

    graph.global_prefactor.num =
        super::NumeratorInterface::admit(GS.emr_mom(out_of_range_edge, GS.cind(1)), None)
            .unwrap()
            .try_into_indexed()
            .unwrap();
    graph.overall_factor = malformed_q3.clone();
    let error = graph.validate_spatial_momentum_convention().unwrap_err();
    assert!(
        error.to_string().contains("Spatial momentum Q3"),
        "Q3 source validation must precede Q bounds in a retained tensor: {error:?}",
    );
    graph.global_prefactor.num = super::NumeratorInterface::admit(Atom::one(), None)
        .unwrap()
        .try_into_indexed()
        .unwrap();

    let nested_q = wrapper.call_args([scalar.pow(malformed_q)]);
    graph.overall_factor = nested_q.clone();
    let error = graph.validate_spatial_momentum_convention().unwrap_err();
    assert!(
        error.to_string().contains("found no edge argument"),
        "a malformed Q in a power exponent must remain visible: {error:?}",
    );

    graph.overall_factor = Atom::one();
    for attribute in ["num", "projector"] {
        let mut source: linnet::parser::DotGraph =
            linnet::parser::DotGraph::from_string(graph.dot_serialize(&Default::default()))
                .unwrap();
        source
            .global_data
            .statements
            .insert(attribute.into(), nested_q.to_canonical_string());
        let error = Graph::from_dot(&source.debug_dot(), &graph.model)
            .err()
            .expect("malformed raw global momentum must be rejected");
        assert!(
            error.to_string().contains("found no edge argument"),
            "malformed raw global {attribute} must retain source validation: {error:?}"
        );
    }

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
    graph.global_prefactor.num = super::NumeratorInterface::admit(Atom::Zero, None)
        .unwrap()
        .try_into_indexed()
        .unwrap();
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
                ext [style=invis num=1]
                ext -> v4 [num=1]
                ext -> v5 [num=1]
                v6 -> ext [num=1]
                v5 -> v4 [lmb_id=0]
                v6 -> v5
                v4 -> v6
            }}
            "#
        );

        let error =
            crate::graph::Graph::from_dot(&malformed, &crate::utils::load_generic_model("scalars"))
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
    let graph: Graph = dot!(digraph missing_hedge_dod{
    graph [projector=1]
                   ext_in [style=invis]
                   ext_out [style=invis]
                   A [num=1 dod=0]
                   ext_in -> A [particle=H mass=0 name=e_in num=1 dod=-2]
                   A -> ext_out [particle=H mass=0 name=e_out num=1 dod=-2]
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
    let graph = Graph::from_dot(
        stringify!(digraph cut_lmb_permutation {
            graph [projector=1 final_state="scalar_1,scalar_1,scalar_1"]
            edge [num=1 particle=scalar_1]
            node [num=1]
            a -> b [id=0 name=first_loop lmb_id=0]
            a -> b [id=1 name=second_loop lmb_id=1]
            a -> b [id=2]
            ext [style=invis]
            ext -> a [id=3 is_cut=0]
            b -> ext [id=4 is_cut=0]
            ext -> a [id=5 is_cut=1]
            b -> ext [id=6 is_cut=1]
        }),
        &crate::utils::load_generic_model("scalars"),
    )
    .unwrap()
    .remove(0);

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
    let graph: Graph = dot!(digraph missing_vertex_dod{
    graph [projector=1]
                   ext_in [style=invis]
                   ext_out [style=invis]
                   A [num=1]
                   ext_in -> A [particle=H mass=0 name=e_in num=1 dod=-2]
                   A -> ext_out [particle=H mass=0 name=e_out num=1 dod=-2]
               })
    .unwrap();

    let dod = &graph.underlying[NodeIndex(0)].dod;
    assert_eq!(dod.value, 0);
    assert!(dod.autogenerated);
}

#[test]
fn open_fundamental_color_requires_an_explicit_projector() {
    test_initialise().unwrap();
    let mut graph: Graph = crate::graph::Graph::from_dot(stringify!(
        digraph open_fundamental_color {
            num = "1/3*spenso::g(spenso::dind(spenso::cof(3,gammalooprs::hedge(0,1))),spenso::cof(3,gammalooprs::hedge(1,1)))"
            ext [style=invis]

            ext -> v1:0 [id=0 particle=d]
            v2:1 -> ext [id=1 particle=d]
            v1:2 -> v2:3 [id=2 particle=d lmb_id=0]
            v1:4 -> v2:5 [id=3 particle=g]
        }
    ), &crate::utils::load_generic_model("sm")).map(|mut graphs| { assert_eq!(graphs.len(),1); graphs.remove(0) })
    .unwrap();

    graph.validate_full_numerator_tensor_network().unwrap();
    graph.global_prefactor.num = super::NumeratorInterface::admit(Atom::one(), None)
        .unwrap()
        .try_into_indexed()
        .unwrap();
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
    let graph: eyre::Result<Graph> = dot!(digraph invalid_vacuum_basis {
 graph [projector=1]
        edge [num=1 mass=2]
        node [num=1]
        a -> b [particle=H id=0 lmb_id=0]
        b -> a [particle=H id=1 lmb_id=1]
    }, &crate::utils::load_generic_model("sm"));
    let error = graph
        .err()
        .expect("an overcomplete loop basis must be rejected");
    assert!(
        format!("{error:?}").contains("do not form a spanning forest"),
        "unexpected basis-validation error: {error:?}"
    );
}

#[test]
fn vacuum_tadpoles_keep_all_loop_generators_without_external_legs() {
    use linnet::half_edge::subgraph::SubSetLike;

    test_initialise().unwrap();
    let graphs: Vec<Graph> = dot!(
        digraph canonical_vacuum {
 graph [projector=1]
            edge [num=1 mass=2]
            node [num=1]
            a -> a [particle=H id=0 lmb_id=0]
            a -> a [particle=H id=1 lmb_id=1]
        }
        digraph specified_vacuum {
 graph [projector=1]
            edge [num=1 mass=2]
            node [num=1]
            a -> a [particle=H id=0 lmb_id=1]
            a -> a [particle=H id=1 lmb_id=0]
        }
        digraph disconnected_vacuum {
 graph [projector=1]
            edge [num=1 mass=2]
            node [num=1]
            a -> a [particle=H id=0 lmb_id=1]
            b -> b [particle=H id=1 lmb_id=0]
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
            .into_graph(&crate::utils::load_generic_model("sm"))
            .unwrap();
        assert_eq!(reparsed.loop_momentum_basis, *lmb);
    }
}
