use std::{collections::BTreeMap, fs};

use figment::providers::Serialized;
use figment::{Figment, Profile};
use linnet::half_edge::involution::EdgeIndex;
use linnet::half_edge::layout::spring::{Constraint, ShiftDirection};
use linnet::half_edge::swap::Swap;
use linnet::{
    dot,
    parser::{set::DotGraphSet, DotGraph},
};
use serde::{Deserialize, Serialize};

use crate::{
    expand_template, graph_archived_compass_subgraph_bytes, graph_compass_subgraph_bytes,
    graph_cycle_basis_bytes, graph_dot_bytes, graph_edges_bytes,
    graph_edges_of_archived_subgraph_bytes, graph_edges_of_bytes, graph_from_spec_bytes,
    graph_info_bytes, graph_nodes_bytes, graph_nodes_of_archived_subgraph_bytes,
    graph_nodes_of_bytes, graph_spanning_forests_bytes, graph_subgraph_bytes, layout_graph_bytes,
    layout_parsed_graph_bytes, layout_parsed_graphs_bytes, parse_dot_graphs_bytes,
    subgraph_contains_hedge_bytes, subgraph_hedges_bytes, subgraph_label_bytes, CBORTypstGraph,
    DotPlacementExpr, PinConstraint, TreeInitCfg, TypstDotEdge, TypstDotGraphInfo, TypstDotNode,
    TypstGraph, TypstPoint,
};

fn test_figment() -> Figment {
    Figment::from(Serialized::from(
        BTreeMap::<String, String>::new(),
        Profile::Default,
    ))
}

fn empty_config_bytes() -> Vec<u8> {
    let mut bytes = Vec::new();
    ciborium::ser::into_writer(&BTreeMap::<String, String>::new(), &mut bytes).unwrap();
    bytes
}

fn decode_graphs(bytes: &[u8]) -> Vec<Vec<u8>> {
    ciborium::de::from_reader(bytes).unwrap()
}

fn encode_cbor<T: Serialize>(value: &T) -> Vec<u8> {
    let mut bytes = Vec::new();
    ciborium::ser::into_writer(value, &mut bytes).unwrap();
    bytes
}

fn decode_cbor<T: for<'de> Deserialize<'de>>(bytes: &[u8]) -> T {
    ciborium::de::from_reader(bytes).unwrap()
}

fn one_statement(key: &str, value: &str) -> BTreeMap<String, String> {
    BTreeMap::from([(key.to_string(), value.to_string())])
}

fn assert_point_close(left: &TypstPoint, right: &TypstPoint) {
    assert!((left.x - right.x).abs() < 1e-9, "{left:?} != {right:?}");
    assert!((left.y - right.y).abs() < 1e-9, "{left:?} != {right:?}");
}

#[test]
fn test_pin_direction_syntaxes_parse() {
    assert!(matches!(
        PinConstraint::parse("x:@+right"),
        Some(PinConstraint::LinkX(group)) if group == "+right"
    ));
    assert!(matches!(
        PinConstraint::parse("x:+@right"),
        Some(PinConstraint::LinkX(group)) if group == "+right"
    ));
    assert!(matches!(
        PinConstraint::parse("y:@-edge0"),
        Some(PinConstraint::LinkY(group)) if group == "-edge0"
    ));
    assert!(matches!(
        PinConstraint::parse("y:-@edge0"),
        Some(PinConstraint::LinkY(group)) if group == "-edge0"
    ));
}

#[test]
fn dot_pos_refs_resolve_by_explicit_node_and_edge_ids() {
    let figment = test_figment();
    let typst_graph = TypstGraph::from_dot(
        dot!(digraph {
            a [id=0 pos="1,2!"]
            b [id=1 pos="ref(node:0)+3,-1!"]
            c [id=2]
            a -> b [id=0 pos="ref(node:1)+0,2!"]
            b -> c [id=1 pos="ref(edge:0)+1,0!"]
        })
        .unwrap(),
        &figment,
    );

    let mut node_positions = BTreeMap::new();
    for (_, _, node) in typst_graph.iter_nodes() {
        node_positions.insert(node.index.unwrap().0, node.pos);
    }

    assert_eq!(node_positions[&0].x, 1.0);
    assert_eq!(node_positions[&0].y, 2.0);
    assert_eq!(node_positions[&1].x, 4.0);
    assert_eq!(node_positions[&1].y, 1.0);

    let mut edge_positions = BTreeMap::new();
    for (_, edge_index, edge) in typst_graph.iter_edges() {
        edge_positions.insert(edge_index.0, edge.data.pos);
    }

    assert_eq!(edge_positions[&0].x, 4.0);
    assert_eq!(edge_positions[&0].y, 3.0);
    assert_eq!(edge_positions[&1].x, 5.0);
    assert_eq!(edge_positions[&1].y, 3.0);
}

#[test]
fn dot_pos_constraint_value_replaces_external_pin_attribute() {
    let figment = test_figment();
    let typst_graph = TypstGraph::from_dot(
        dot!(digraph {
            ext [style=invis]
            a
            ext -> a [id=0 pos="x:@-left!,y:@edge0!"]
        })
        .unwrap(),
        &figment,
    );

    let (_, _, edge) = typst_graph.iter_edges().next().unwrap();
    assert!(matches!(
        edge.data.constraints.x,
        Constraint::Grouped(_, ShiftDirection::NegativeOnly)
    ));
    assert!(matches!(
        edge.data.constraints.y,
        Constraint::Grouped(_, ShiftDirection::Any)
    ));
    assert_eq!(
        edge.data.statements.get("pos").map(String::as_str),
        Some("0,0")
    );
    assert!(!edge.data.statements.contains_key("pin"));
}

#[test]
fn dot_pos_axis_constraints_support_numeric_partial_pins() {
    let figment = test_figment();
    let typst_graph = TypstGraph::from_dot(
        dot!(digraph {
            a [id=0 pos="x:2!"]
            b [id=1 pos="y:-1!"]
            a -> b
        })
        .unwrap(),
        &figment,
    );

    let mut nodes = BTreeMap::new();
    for (_, _, node) in typst_graph.iter_nodes() {
        nodes.insert(node.index.unwrap().0, node);
    }

    assert_eq!(nodes[&0].pos.x, 2.0);
    assert_eq!(nodes[&0].pos.y, 0.0);
    assert!(matches!(nodes[&0].constraints.x, Constraint::Fixed));
    assert!(matches!(nodes[&0].constraints.y, Constraint::Free));

    assert_eq!(nodes[&1].pos.x, 0.0);
    assert_eq!(nodes[&1].pos.y, -1.0);
    assert!(matches!(nodes[&1].constraints.x, Constraint::Free));
    assert!(matches!(nodes[&1].constraints.y, Constraint::Fixed));
}

#[test]
fn dot_pos_axis_bangs_apply_per_axis() {
    let figment = test_figment();
    let typst_graph = TypstGraph::from_dot(
        dot!(digraph {
            a [id=0 pos="x:2!,y:3"]
            b [id=1 pos="x:4,y:-1!"]
            a -> b
        })
        .unwrap(),
        &figment,
    );

    let mut nodes = BTreeMap::new();
    for (_, _, node) in typst_graph.iter_nodes() {
        nodes.insert(node.index.unwrap().0, node);
    }

    assert_eq!(nodes[&0].pos.x, 2.0);
    assert_eq!(nodes[&0].pos.y, 3.0);
    assert!(matches!(nodes[&0].constraints.x, Constraint::Fixed));
    assert!(matches!(nodes[&0].constraints.y, Constraint::Free));

    assert_eq!(nodes[&1].pos.x, 4.0);
    assert_eq!(nodes[&1].pos.y, -1.0);
    assert!(matches!(nodes[&1].constraints.x, Constraint::Free));
    assert!(matches!(nodes[&1].constraints.y, Constraint::Fixed));
}

#[test]
fn dot_pos_group_axis_requires_axis_bang() {
    let err = DotPlacementExpr::parse("x:@left,y:0").unwrap_err();
    assert!(err.contains("must be pinned with !"));
}

#[test]
fn test_template_expansion_replaces_known_keys_and_escapes_braces() {
    let statements = BTreeMap::from([("label".to_string(), "\"a-c\"".to_string())]);

    assert_eq!(expand_template("[{label}]", &statements), "[a-c]");
    assert_eq!(expand_template("{{label}}", &statements), "{label}");
    assert_eq!(expand_template("{missing}", &statements), "{missing}");
}

#[derive(Serialize)]
struct TestGraphSpec {
    name: String,
    #[serde(default)]
    statements: BTreeMap<String, String>,
    nodes: Vec<TestNodeSpec>,
    edges: Vec<TestEdgeSpec>,
}

#[derive(Serialize)]
#[serde(rename_all = "kebab-case")]
struct TestTemplatedGraphSpec {
    name: String,
    #[serde(default)]
    statements: BTreeMap<String, String>,
    #[serde(default)]
    edge_statements: BTreeMap<String, String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    source_style_eval: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    sink_style_eval: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    label_eval: Option<String>,
    #[serde(default)]
    node_statements: BTreeMap<String, String>,
    nodes: Vec<TestNodeSpec>,
    edges: Vec<TestEdgeSpec>,
}

#[derive(Serialize)]
struct TestNodeSpec {
    name: String,
    #[serde(default)]
    statements: BTreeMap<String, String>,
}

#[derive(Serialize)]
struct TestEdgeSpec {
    #[serde(skip_serializing_if = "Option::is_none")]
    source: Option<TestEndpointSpec>,
    #[serde(skip_serializing_if = "Option::is_none")]
    sink: Option<TestEndpointSpec>,
    #[serde(default)]
    statements: BTreeMap<String, String>,
}

#[derive(Serialize)]
struct TestPlacementGraphSpec {
    name: String,
    nodes: Vec<TestPlacedNodeSpec>,
    edges: Vec<TestPlacedEdgeSpec>,
}

#[derive(Serialize)]
struct TestPlacedNodeSpec {
    name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pos: Option<TestPlacementSpec>,
    #[serde(default)]
    statements: BTreeMap<String, String>,
}

#[derive(Serialize)]
struct TestPlacedEdgeSpec {
    #[serde(skip_serializing_if = "Option::is_none")]
    source: Option<TestEndpointSpec>,
    #[serde(skip_serializing_if = "Option::is_none")]
    sink: Option<TestEndpointSpec>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pos: Option<TestPlacementSpec>,
    #[serde(default)]
    statements: BTreeMap<String, String>,
}

#[derive(Serialize)]
#[serde(rename_all = "kebab-case")]
struct TestPlacementSpec {
    #[serde(skip_serializing_if = "Option::is_none")]
    mode: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    x: Option<TestPlacementCoord>,
    #[serde(skip_serializing_if = "Option::is_none")]
    y: Option<TestPlacementCoord>,
    #[serde(skip_serializing_if = "Option::is_none", rename = "ref")]
    reference: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    dx: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    dy: Option<f64>,
}

#[derive(Serialize)]
#[serde(untagged)]
enum TestPlacementCoord {
    Number(f64),
    Group(TestPlacementGroup),
}

#[derive(Serialize)]
struct TestPlacementGroup {
    kind: &'static str,
    name: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    side: Option<&'static str>,
}

#[derive(Serialize)]
struct TestEndpointSpec {
    node: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    compass: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    statement: Option<String>,
}

#[test]
fn test_parse_pass_returns_archived_dot_graphs() {
    let parsed =
        parse_dot_graphs_bytes(br#"digraph first { a -> b } digraph second { c -> d }"#).unwrap();

    let archived = decode_graphs(&parsed);

    assert_eq!(archived.len(), 2);
    let first = DotGraph::archived_view(&archived[0]);
    let second = DotGraph::archived_view(&archived[1]);
    assert_eq!(first.global_data().name.as_str(), "first");
    assert_eq!(second.global_data().name.as_str(), "second");
}

#[test]
fn test_graph_query_api_reads_nodes_edges_and_subgraphs() {
    let parsed = parse_dot_graphs_bytes(
        br#"digraph first {
            graph [full_num="x+y"];
            a [eval="left"];
            b [eval="right"];
            a:out:n -> b:in:s [label="ab"];
            ext -> a [dir=back, source="src"];
        }"#,
    )
    .unwrap();

    let graphs = decode_graphs(&parsed);
    let graph = &graphs[0];

    let info: TypstDotGraphInfo =
        ciborium::de::from_reader(graph_info_bytes(graph).unwrap().as_slice()).unwrap();
    assert_eq!(info.name, "first");
    assert_eq!(
        info.global_statements.get("full_num").map(String::as_str),
        Some("x+y")
    );

    let nodes: Vec<TypstDotNode> =
        ciborium::de::from_reader(graph_nodes_bytes(graph).unwrap().as_slice()).unwrap();
    assert_eq!(nodes.len(), 3);
    assert_eq!(nodes[0].name.as_deref(), Some("a"));

    let edges: Vec<TypstDotEdge> =
        ciborium::de::from_reader(graph_edges_bytes(graph).unwrap().as_slice()).unwrap();
    assert_eq!(edges.len(), 2);
    assert_eq!(edges[0].orientation, "default");

    let north = graph_compass_subgraph_bytes(graph, &encode_cbor(&"n")).unwrap();
    let north_label: String = ciborium::de::from_reader(north.as_slice()).unwrap();

    let north_edges: Vec<TypstDotEdge> = ciborium::de::from_reader(
        graph_edges_of_bytes(graph, &encode_cbor(&north_label))
            .unwrap()
            .as_slice(),
    )
    .unwrap();
    assert_eq!(north_edges.len(), 1);

    let view = DotGraph::archived_view(graph);
    let mut south_bits = vec![false; view.n_hedges()];
    south_bits[view.n_hedges() - 1] = true;
    let subgraph = graph_subgraph_bytes(graph, &encode_cbor(&south_bits)).unwrap();
    let south_label: String = ciborium::de::from_reader(subgraph.as_slice()).unwrap();
    let south_nodes: Vec<TypstDotNode> = ciborium::de::from_reader(
        graph_nodes_of_bytes(graph, &encode_cbor(&south_label))
            .unwrap()
            .as_slice(),
    )
    .unwrap();
    assert_eq!(south_nodes.len(), 1);
}

#[test]
fn test_graph_spec_constructor_reads_nodes_edges_and_subgraphs() {
    let spec = TestGraphSpec {
        name: "constructed".to_string(),
        statements: one_statement("full_num", "x+y"),
        nodes: vec![
            TestNodeSpec {
                name: "a".to_string(),
                statements: BTreeMap::new(),
            },
            TestNodeSpec {
                name: "b".to_string(),
                statements: BTreeMap::new(),
            },
        ],
        edges: vec![
            TestEdgeSpec {
                source: Some(TestEndpointSpec {
                    node: 0,
                    compass: Some("e".to_string()),
                    statement: Some("out".to_string()),
                }),
                sink: Some(TestEndpointSpec {
                    node: 1,
                    compass: Some("w".to_string()),
                    statement: Some("in".to_string()),
                }),
                statements: one_statement("label", "a-to-b"),
            },
            TestEdgeSpec {
                source: None,
                sink: Some(TestEndpointSpec {
                    node: 0,
                    compass: Some("n".to_string()),
                    statement: None,
                }),
                statements: one_statement("label", "incoming"),
            },
            TestEdgeSpec {
                source: Some(TestEndpointSpec {
                    node: 1,
                    compass: Some("s".to_string()),
                    statement: None,
                }),
                sink: None,
                statements: one_statement("label", "outgoing"),
            },
        ],
    };

    let graph = graph_from_spec_bytes(&encode_cbor(&spec)).unwrap();
    let info: TypstDotGraphInfo =
        ciborium::de::from_reader(graph_info_bytes(&graph).unwrap().as_slice()).unwrap();
    assert_eq!(info.name, "constructed");
    assert_eq!(
        info.global_statements.get("full_num").map(String::as_str),
        Some("x+y")
    );

    let nodes: Vec<TypstDotNode> =
        ciborium::de::from_reader(graph_nodes_bytes(&graph).unwrap().as_slice()).unwrap();
    assert_eq!(nodes.len(), 2);

    let edges: Vec<TypstDotEdge> =
        ciborium::de::from_reader(graph_edges_bytes(&graph).unwrap().as_slice()).unwrap();
    assert_eq!(edges.len(), 3);

    let dot: String =
        ciborium::de::from_reader(graph_dot_bytes(&graph).unwrap().as_slice()).unwrap();
    assert!(dot.contains("digraph constructed"));
    assert!(dot.contains("a"));
    assert!(dot.contains("b"));

    let north = graph_compass_subgraph_bytes(&graph, &encode_cbor(&"n")).unwrap();
    let north_label: String = ciborium::de::from_reader(north.as_slice()).unwrap();
    let north_edges: Vec<TypstDotEdge> = ciborium::de::from_reader(
        graph_edges_of_bytes(&graph, &encode_cbor(&north_label))
            .unwrap()
            .as_slice(),
    )
    .unwrap();
    assert_eq!(north_edges.len(), 1);

    let internal =
        graph_subgraph_bytes(&graph, &encode_cbor(&vec![true, true, false, false])).unwrap();
    let internal_label: String = ciborium::de::from_reader(internal.as_slice()).unwrap();
    let internal_nodes: Vec<TypstDotNode> = ciborium::de::from_reader(
        graph_nodes_of_bytes(&graph, &encode_cbor(&internal_label))
            .unwrap()
            .as_slice(),
    )
    .unwrap();
    assert_eq!(internal_nodes.len(), 2);
}

#[test]
fn test_graph_spec_payloads_are_opaque_and_survive_layout() {
    #[derive(Serialize)]
    struct PayloadGraphSpec {
        name: String,
        payload: Vec<u8>,
        nodes: Vec<PayloadNodeSpec>,
        edges: Vec<PayloadEdgeSpec>,
    }

    #[derive(Serialize)]
    struct PayloadNodeSpec {
        name: String,
        payload: Vec<u8>,
    }

    #[derive(Serialize)]
    struct PayloadEdgeSpec {
        payload: Vec<u8>,
        source: PayloadEndpointSpec,
        sink: PayloadEndpointSpec,
    }

    #[derive(Serialize)]
    struct PayloadEndpointSpec {
        node: usize,
        payload: Vec<u8>,
    }

    fn assert_payloads(graph: &[u8]) {
        let info: TypstDotGraphInfo = decode_cbor(&graph_info_bytes(graph).unwrap());
        assert_eq!(info.payload.as_deref(), Some(&b"graph"[..]));

        let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(graph).unwrap());
        assert_eq!(nodes[0].payload.as_deref(), Some(&b"node-a"[..]));
        assert_eq!(nodes[1].payload.as_deref(), Some(&b"node-b"[..]));

        let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(graph).unwrap());
        assert_eq!(edges[0].payload.as_deref(), Some(&b"edge"[..]));
        assert_eq!(
            edges[0]
                .source
                .as_ref()
                .and_then(|source| source.payload.as_deref()),
            Some(&b"source"[..])
        );
        assert_eq!(
            edges[0]
                .sink
                .as_ref()
                .and_then(|sink| sink.payload.as_deref()),
            Some(&b"sink"[..])
        );
    }

    let spec = PayloadGraphSpec {
        name: "payloads".to_string(),
        payload: b"graph".to_vec(),
        nodes: vec![
            PayloadNodeSpec {
                name: "a".to_string(),
                payload: b"node-a".to_vec(),
            },
            PayloadNodeSpec {
                name: "b".to_string(),
                payload: b"node-b".to_vec(),
            },
        ],
        edges: vec![PayloadEdgeSpec {
            payload: b"edge".to_vec(),
            source: PayloadEndpointSpec {
                node: 0,
                payload: b"source".to_vec(),
            },
            sink: PayloadEndpointSpec {
                node: 1,
                payload: b"sink".to_vec(),
            },
        }],
    };

    let graph = graph_from_spec_bytes(&encode_cbor(&spec)).unwrap();
    assert_payloads(&graph);

    let laid_out = layout_parsed_graph_bytes(&graph, &empty_config_bytes()).unwrap();
    assert_payloads(&laid_out);
}

#[test]
fn test_graph_spec_exposes_half_edge_ids() {
    #[derive(Serialize)]
    struct HalfIdGraphSpec {
        name: String,
        nodes: Vec<TestNodeSpec>,
        edges: Vec<HalfIdEdgeSpec>,
    }

    #[derive(Serialize)]
    struct HalfIdEdgeSpec {
        id: usize,
        source: HalfIdEndpointSpec,
        sink: HalfIdEndpointSpec,
    }

    #[derive(Serialize)]
    struct HalfIdEndpointSpec {
        node: usize,
        id: usize,
    }

    let graph = graph_from_spec_bytes(&encode_cbor(&HalfIdGraphSpec {
        name: "half-ids".to_string(),
        nodes: vec![
            TestNodeSpec {
                name: "a".to_string(),
                statements: BTreeMap::new(),
            },
            TestNodeSpec {
                name: "b".to_string(),
                statements: BTreeMap::new(),
            },
        ],
        edges: vec![HalfIdEdgeSpec {
            id: 5,
            source: HalfIdEndpointSpec { node: 0, id: 7 },
            sink: HalfIdEndpointSpec { node: 1, id: 11 },
        }],
    }))
    .unwrap();

    let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(&graph).unwrap());
    assert_eq!(edges[0].id, Some(5));
    assert_eq!(edges[0].source.as_ref().unwrap().id, Some(7));
    assert_eq!(edges[0].sink.as_ref().unwrap().id, Some(11));
}

#[test]
fn test_graph_spec_uses_first_class_placements() {
    let spec = TestPlacementGraphSpec {
        name: "placed".to_string(),
        nodes: vec![
            TestPlacedNodeSpec {
                name: "a".to_string(),
                pos: Some(TestPlacementSpec {
                    mode: Some("pin"),
                    x: Some(TestPlacementCoord::Number(1.0)),
                    y: Some(TestPlacementCoord::Number(2.0)),
                    reference: None,
                    dx: None,
                    dy: None,
                }),
                statements: one_statement("label", "left"),
            },
            TestPlacedNodeSpec {
                name: "b".to_string(),
                pos: Some(TestPlacementSpec {
                    mode: Some("pin"),
                    x: None,
                    y: None,
                    reference: Some(0),
                    dx: Some(3.0),
                    dy: Some(-1.0),
                }),
                statements: one_statement("label", "right"),
            },
        ],
        edges: vec![TestPlacedEdgeSpec {
            source: Some(TestEndpointSpec {
                node: 0,
                compass: None,
                statement: None,
            }),
            sink: Some(TestEndpointSpec {
                node: 1,
                compass: None,
                statement: None,
            }),
            pos: Some(TestPlacementSpec {
                mode: Some("pin"),
                x: Some(TestPlacementCoord::Group(TestPlacementGroup {
                    kind: "group",
                    name: "edge-column",
                    side: Some("+"),
                })),
                y: Some(TestPlacementCoord::Number(0.5)),
                reference: None,
                dx: None,
                dy: None,
            }),
            statements: BTreeMap::new(),
        }],
    };

    let graph = graph_from_spec_bytes(&encode_cbor(&spec)).unwrap();
    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&graph).unwrap());
    let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(&graph).unwrap());

    assert_eq!(nodes[0].pos, Some(crate::TypstPoint { x: 1.0, y: 2.0 }));
    assert_eq!(nodes[1].pos, Some(crate::TypstPoint { x: 4.0, y: 1.0 }));
    assert_eq!(
        nodes[0].statements.get("label").map(String::as_str),
        Some("left")
    );
    assert!(!nodes[0].statements.contains_key("pos"));
    assert!(!nodes[0].statements.contains_key("pin"));
    assert_eq!(edges[0].pos, Some(crate::TypstPoint { x: 0.0, y: 0.5 }));
    assert!(!edges[0].statements.contains_key("pos"));
    assert!(!edges[0].statements.contains_key("pin"));

    let laid_out = layout_parsed_graph_bytes(&graph, &empty_config_bytes()).unwrap();
    let laid_out_nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&laid_out).unwrap());
    assert_eq!(
        laid_out_nodes[0].pos,
        Some(crate::TypstPoint { x: 1.0, y: 2.0 })
    );
    assert_eq!(
        laid_out_nodes[1].pos,
        Some(crate::TypstPoint { x: 4.0, y: 1.0 })
    );
}

#[test]
fn test_graph_spec_placement_defaults_to_pin() {
    let spec = TestPlacementGraphSpec {
        name: "placed-default".to_string(),
        nodes: vec![TestPlacedNodeSpec {
            name: "a".to_string(),
            pos: Some(TestPlacementSpec {
                mode: None,
                x: Some(TestPlacementCoord::Number(1.0)),
                y: Some(TestPlacementCoord::Number(2.0)),
                reference: None,
                dx: None,
                dy: None,
            }),
            statements: BTreeMap::new(),
        }],
        edges: Vec::new(),
    };

    let graph = graph_from_spec_bytes(&encode_cbor(&spec)).unwrap();
    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&graph).unwrap());
    assert_eq!(nodes[0].pos, Some(crate::TypstPoint { x: 1.0, y: 2.0 }));

    let laid_out = layout_parsed_graph_bytes(&graph, &empty_config_bytes()).unwrap();
    let laid_out_nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&laid_out).unwrap());
    assert_eq!(
        laid_out_nodes[0].pos,
        Some(crate::TypstPoint { x: 1.0, y: 2.0 })
    );
}

#[test]
fn test_graph_spec_constructor_expands_default_statement_templates() {
    let spec = TestTemplatedGraphSpec {
        name: "templated".to_string(),
        statements: BTreeMap::new(),
        edge_statements: BTreeMap::new(),
        source_style_eval: Some("(stroke: red + 0.5pt)".to_string()),
        sink_style_eval: Some("(stroke: blue + 0.5pt)".to_string()),
        label_eval: Some("(text(fill: rgb(\"#{color}\"))[{label}])".to_string()),
        node_statements: BTreeMap::from([(
            "eval".to_string(),
            "(fill: rgb(\"#{color}\"))".to_string(),
        )]),
        nodes: vec![
            TestNodeSpec {
                name: "a".to_string(),
                statements: one_statement("color", "ff0000"),
            },
            TestNodeSpec {
                name: "b".to_string(),
                statements: one_statement("color", "00aa00"),
            },
        ],
        edges: vec![TestEdgeSpec {
            source: Some(TestEndpointSpec {
                node: 0,
                compass: Some("e".to_string()),
                statement: None,
            }),
            sink: Some(TestEndpointSpec {
                node: 1,
                compass: Some("w".to_string()),
                statement: None,
            }),
            statements: BTreeMap::from([
                ("color".to_string(), "0055ff".to_string()),
                ("label".to_string(), "a-c".to_string()),
            ]),
        }],
    };

    let graph = graph_from_spec_bytes(&encode_cbor(&spec)).unwrap();
    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&graph).unwrap());
    assert_eq!(nodes[0].eval.as_deref(), Some("(fill: rgb(\"#ff0000\"))"));
    assert_eq!(nodes[1].eval.as_deref(), Some("(fill: rgb(\"#00aa00\"))"));

    let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(&graph).unwrap());
    assert_eq!(
        edges[0].label_eval.as_deref(),
        Some("(text(fill: rgb(\"#0055ff\"))[a-c])")
    );
    assert_eq!(
        edges[0].source_style_eval.as_deref(),
        Some("(stroke: red + 0.5pt)")
    );
    assert_eq!(
        edges[0].sink_style_eval.as_deref(),
        Some("(stroke: blue + 0.5pt)")
    );
    assert_eq!(
        edges[0].statements.get("label-eval").map(String::as_str),
        Some("(text(fill: rgb(\"#0055ff\"))[a-c])")
    );
}

#[test]
fn test_archived_graph_and_subgraph_api() {
    let graph = graph_from_spec_bytes(&encode_cbor(&TestGraphSpec {
        name: "constructed".to_string(),
        statements: BTreeMap::new(),
        nodes: vec![
            TestNodeSpec {
                name: "a".to_string(),
                statements: BTreeMap::new(),
            },
            TestNodeSpec {
                name: "b".to_string(),
                statements: BTreeMap::new(),
            },
        ],
        edges: vec![TestEdgeSpec {
            source: Some(TestEndpointSpec {
                node: 0,
                compass: Some("e".to_string()),
                statement: None,
            }),
            sink: Some(TestEndpointSpec {
                node: 1,
                compass: Some("w".to_string()),
                statement: None,
            }),
            statements: one_statement("join", "ab"),
        }],
    }))
    .unwrap();

    let east = graph_archived_compass_subgraph_bytes(&graph, &encode_cbor(&"e")).unwrap();
    let east_label: String = decode_cbor(&subgraph_label_bytes(&east).unwrap());
    assert!(!east_label.is_empty());
    let east_hedges: Vec<usize> = decode_cbor(&subgraph_hedges_bytes(&east).unwrap());
    assert_eq!(east_hedges, vec![0]);
    let contains: bool =
        decode_cbor(&subgraph_contains_hedge_bytes(&east, &encode_cbor(&0)).unwrap());
    assert!(contains);

    let east_edges: Vec<TypstDotEdge> =
        decode_cbor(&graph_edges_of_archived_subgraph_bytes(&graph, &east).unwrap());
    assert_eq!(east_edges.len(), 1);
    let east_nodes: Vec<TypstDotNode> =
        decode_cbor(&graph_nodes_of_archived_subgraph_bytes(&graph, &east).unwrap());
    assert_eq!(east_nodes.len(), 1);

    let cycles: Vec<Vec<u8>> = decode_cbor(&graph_cycle_basis_bytes(&graph).unwrap());
    assert!(cycles.is_empty());
    let forests: Vec<Vec<u8>> = decode_cbor(&graph_spanning_forests_bytes(&graph).unwrap());
    assert_eq!(forests.len(), 1);
}

#[test]
fn test_graph_spec_expands_default_statement_templates() {
    let graph = graph_from_spec_bytes(&encode_cbor(&TestTemplatedGraphSpec {
        name: "constructed".to_string(),
        statements: BTreeMap::new(),
        edge_statements: BTreeMap::new(),
        source_style_eval: Some("(stroke: red + 0.5pt)".to_string()),
        sink_style_eval: Some("(stroke: blue + 0.5pt)".to_string()),
        label_eval: Some("(text(fill: rgb(\"#{color}\"))[{label}])".to_string()),
        node_statements: BTreeMap::from([(
            "eval".to_string(),
            "(fill: rgb(\"#{color}\"))".to_string(),
        )]),
        nodes: vec![
            TestNodeSpec {
                name: "a".to_string(),
                statements: one_statement("color", "ff0000"),
            },
            TestNodeSpec {
                name: "b".to_string(),
                statements: one_statement("color", "00aa00"),
            },
        ],
        edges: vec![TestEdgeSpec {
            source: Some(TestEndpointSpec {
                node: 0,
                compass: None,
                statement: None,
            }),
            sink: Some(TestEndpointSpec {
                node: 1,
                compass: None,
                statement: None,
            }),
            statements: BTreeMap::from([
                ("color".to_string(), "0055ff".to_string()),
                ("label".to_string(), "a-c".to_string()),
            ]),
        }],
    }))
    .unwrap();

    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&graph).unwrap());
    assert_eq!(nodes[0].eval.as_deref(), Some("(fill: rgb(\"#ff0000\"))"));
    let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(&graph).unwrap());
    assert_eq!(
        edges[0].label_eval.as_deref(),
        Some("(text(fill: rgb(\"#0055ff\"))[a-c])")
    );
    assert_eq!(
        edges[0].source_style_eval.as_deref(),
        Some("(stroke: red + 0.5pt)")
    );
    assert_eq!(
        edges[0].sink_style_eval.as_deref(),
        Some("(stroke: blue + 0.5pt)")
    );
}

#[test]
fn test_graph_join_matches_half_edge_statement() {
    let left = graph_from_spec_bytes(&encode_cbor(&TestGraphSpec {
        name: "left".to_string(),
        statements: BTreeMap::new(),
        nodes: vec![TestNodeSpec {
            name: "a".to_string(),
            statements: BTreeMap::new(),
        }],
        edges: vec![TestEdgeSpec {
            source: None,
            sink: Some(TestEndpointSpec {
                node: 0,
                compass: None,
                statement: Some("join".to_string()),
            }),
            statements: BTreeMap::new(),
        }],
    }))
    .unwrap();
    let right = graph_from_spec_bytes(&encode_cbor(&TestGraphSpec {
        name: "right".to_string(),
        statements: BTreeMap::new(),
        nodes: vec![TestNodeSpec {
            name: "b".to_string(),
            statements: BTreeMap::new(),
        }],
        edges: vec![TestEdgeSpec {
            source: Some(TestEndpointSpec {
                node: 0,
                compass: None,
                statement: Some("join".to_string()),
            }),
            sink: None,
            statements: BTreeMap::new(),
        }],
    }))
    .unwrap();

    let joined = crate::graph_join_by_hedge_key_bytes(
        &left,
        &right,
        &encode_cbor(&one_statement("key", "statement")),
    )
    .unwrap();
    let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(&joined).unwrap());
    assert_eq!(edges.len(), 1);
    assert!(edges[0].source.is_some());
    assert!(edges[0].sink.is_some());
}

#[test]
fn test_split_parse_layout_matches_wrapper() {
    let input = br#"digraph archived { a [pin="x:1.0"]; a -> b }"#;
    let config = empty_config_bytes();

    let parsed = parse_dot_graphs_bytes(input).unwrap();
    let split = layout_parsed_graphs_bytes(&parsed, &config).unwrap();
    let wrapped = layout_graph_bytes(input, &config).unwrap();

    assert_eq!(split, wrapped);

    let output = decode_graphs(&split);
    assert_eq!(output.len(), 1);
}

#[test]
fn test_single_graph_layout_mutates_graph_bytes() {
    let parsed =
        parse_dot_graphs_bytes(br#"digraph archived { a [pin="x:1.0"]; a -> b }"#).unwrap();
    let graphs = decode_graphs(&parsed);
    let laid_out = layout_parsed_graph_bytes(&graphs[0], &empty_config_bytes()).unwrap();

    let nodes: Vec<TypstDotNode> =
        ciborium::de::from_reader(graph_nodes_bytes(&laid_out).unwrap().as_slice()).unwrap();
    assert_eq!(nodes.len(), 2);
    assert!(nodes.iter().all(|node| node.pos.is_some()));
}

#[test]
fn test_dot_layout_layers_directed_graph() {
    let graph = graph_from_spec_bytes(&encode_cbor(&TestGraphSpec {
        name: "dot".to_string(),
        statements: BTreeMap::new(),
        nodes: vec![
            TestNodeSpec {
                name: "a".to_string(),
                statements: BTreeMap::new(),
            },
            TestNodeSpec {
                name: "b".to_string(),
                statements: BTreeMap::new(),
            },
            TestNodeSpec {
                name: "c".to_string(),
                statements: BTreeMap::new(),
            },
            TestNodeSpec {
                name: "d".to_string(),
                statements: BTreeMap::new(),
            },
        ],
        edges: vec![
            TestEdgeSpec {
                source: Some(TestEndpointSpec {
                    node: 0,
                    compass: None,
                    statement: None,
                }),
                sink: Some(TestEndpointSpec {
                    node: 1,
                    compass: None,
                    statement: None,
                }),
                statements: BTreeMap::new(),
            },
            TestEdgeSpec {
                source: Some(TestEndpointSpec {
                    node: 0,
                    compass: None,
                    statement: None,
                }),
                sink: Some(TestEndpointSpec {
                    node: 2,
                    compass: None,
                    statement: None,
                }),
                statements: BTreeMap::new(),
            },
            TestEdgeSpec {
                source: Some(TestEndpointSpec {
                    node: 1,
                    compass: None,
                    statement: None,
                }),
                sink: Some(TestEndpointSpec {
                    node: 3,
                    compass: None,
                    statement: None,
                }),
                statements: BTreeMap::new(),
            },
            TestEdgeSpec {
                source: Some(TestEndpointSpec {
                    node: 2,
                    compass: None,
                    statement: None,
                }),
                sink: Some(TestEndpointSpec {
                    node: 3,
                    compass: None,
                    statement: None,
                }),
                statements: BTreeMap::new(),
            },
        ],
    }))
    .unwrap();

    let laid_out = layout_parsed_graph_bytes(
        &graph,
        &encode_cbor(&BTreeMap::from([
            ("layout-algo".to_string(), "dot".to_string()),
            ("label-steps".to_string(), "0".to_string()),
        ])),
    )
    .unwrap();
    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&laid_out).unwrap());
    let by_name = nodes
        .iter()
        .map(|node| (node.name.as_deref().unwrap(), node.pos.as_ref().unwrap()))
        .collect::<BTreeMap<_, _>>();

    assert!(by_name["a"].y < by_name["b"].y);
    assert_eq!(by_name["b"].y, by_name["c"].y);
    assert!(by_name["b"].y < by_name["d"].y);
}

#[test]
fn test_tree_layout_uses_explicit_root_order() {
    let parsed = parse_dot_graphs_bytes(
        br#"digraph roots {
            a -> b
            b -> c
            d -> e
        }"#,
    )
    .unwrap();
    let graph = decode_graphs(&parsed).remove(0);

    let laid_out = layout_parsed_graph_bytes(
        &graph,
        &encode_cbor(&ciborium::Value::Map(vec![
            (
                ciborium::Value::Text("layout-algo".to_string()),
                ciborium::Value::Text("tree".to_string()),
            ),
            (
                ciborium::Value::Text("layout-roots".to_string()),
                ciborium::Value::Array(vec![
                    ciborium::Value::Integer(3.into()),
                    ciborium::Value::Integer(0.into()),
                ]),
            ),
            (
                ciborium::Value::Text("label-steps".to_string()),
                ciborium::Value::Text("0".to_string()),
            ),
        ])),
    )
    .unwrap();
    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&laid_out).unwrap());
    let by_name = nodes
        .iter()
        .map(|node| (node.name.as_deref().unwrap(), node.pos.as_ref().unwrap()))
        .collect::<BTreeMap<_, _>>();

    assert!(by_name["d"].x < by_name["a"].x);
    assert_eq!(by_name["d"].y, by_name["a"].y);
}

#[test]
fn test_dot_layout_uses_explicit_roots_as_first_rank() {
    let parsed = parse_dot_graphs_bytes(
        br#"digraph roots {
            a -> b
            c -> d
        }"#,
    )
    .unwrap();
    let graph = decode_graphs(&parsed).remove(0);

    let laid_out = layout_parsed_graph_bytes(
        &graph,
        &encode_cbor(&ciborium::Value::Map(vec![
            (
                ciborium::Value::Text("layout-algo".to_string()),
                ciborium::Value::Text("dot".to_string()),
            ),
            (
                ciborium::Value::Text("layout-roots".to_string()),
                ciborium::Value::Integer(2.into()),
            ),
            (
                ciborium::Value::Text("label-steps".to_string()),
                ciborium::Value::Text("0".to_string()),
            ),
        ])),
    )
    .unwrap();
    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&laid_out).unwrap());
    let by_name = nodes
        .iter()
        .map(|node| (node.name.as_deref().unwrap(), node.pos.as_ref().unwrap()))
        .collect::<BTreeMap<_, _>>();

    assert!(by_name["c"].x < by_name["a"].x);
    assert_eq!(by_name["c"].y, by_name["a"].y);
    assert!(by_name["c"].y < by_name["d"].y);
}

#[test]
fn test_partial_tree_layout_keeps_non_tree_edges_straight() {
    let parsed =
        parse_dot_graphs_bytes(br#"digraph partial { a -> b; b -> c; c -> d; d -> a; a -> c }"#)
            .unwrap();
    let graph = decode_graphs(&parsed).remove(0);
    let forests: Vec<Vec<u8>> = decode_cbor(&graph_spanning_forests_bytes(&graph).unwrap());
    let forest_label: String = decode_cbor(&subgraph_label_bytes(&forests[0]).unwrap());
    let forest_hedges: Vec<usize> = decode_cbor(&subgraph_hedges_bytes(&forests[0]).unwrap());

    let laid_out = layout_parsed_graph_bytes(
        &graph,
        &encode_cbor(&BTreeMap::from([
            ("layout-algo".to_string(), "tree".to_string()),
            ("subgraph".to_string(), forest_label),
            ("label-steps".to_string(), "0".to_string()),
        ])),
    )
    .unwrap();

    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&laid_out).unwrap());
    let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(&laid_out).unwrap());
    let non_tree_edge = edges
        .iter()
        .find(|edge| {
            let source = edge.source.as_ref().unwrap().hedge;
            let sink = edge.sink.as_ref().unwrap().hedge;
            !(forest_hedges.contains(&source) && forest_hedges.contains(&sink))
        })
        .unwrap();
    let source = nodes[non_tree_edge.source.as_ref().unwrap().node]
        .pos
        .as_ref()
        .unwrap();
    let sink = nodes[non_tree_edge.sink.as_ref().unwrap().node]
        .pos
        .as_ref()
        .unwrap();
    let expected_midpoint = TypstPoint {
        x: 0.5 * (source.x + sink.x),
        y: 0.5 * (source.y + sink.y),
    };

    assert_point_close(non_tree_edge.pos.as_ref().unwrap(), &expected_midpoint);
}

#[test]
fn test_fixed_node_subgraph_tree_layout_moves_only_selected_edges() {
    let parsed = parse_dot_graphs_bytes(
        br#"digraph fixed_edges {
            a [pos="0,0"]
            b [pos="4,0"]
            c [pos="8,0"]
            a -> b
            b -> c [pos="88,88"]
        }"#,
    )
    .unwrap();
    let graph = decode_graphs(&parsed).remove(0);
    let selected_label: String = decode_cbor(
        &graph_subgraph_bytes(&graph, &encode_cbor(&vec![true, true, false, false])).unwrap(),
    );

    let laid_out = layout_parsed_graph_bytes(
        &graph,
        &encode_cbor(&BTreeMap::from([
            ("layout-algo".to_string(), "tree".to_string()),
            ("layout-nodes".to_string(), "fixed".to_string()),
            ("subgraph".to_string(), selected_label),
            ("label-steps".to_string(), "0".to_string()),
        ])),
    )
    .unwrap();

    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&laid_out).unwrap());
    let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(&laid_out).unwrap());

    assert_point_close(
        nodes[0].pos.as_ref().unwrap(),
        &TypstPoint { x: 0.0, y: 0.0 },
    );
    assert_point_close(
        nodes[1].pos.as_ref().unwrap(),
        &TypstPoint { x: 4.0, y: 0.0 },
    );
    assert_point_close(
        nodes[2].pos.as_ref().unwrap(),
        &TypstPoint { x: 8.0, y: 0.0 },
    );
    assert_point_close(
        edges[0].pos.as_ref().unwrap(),
        &TypstPoint { x: 2.0, y: 0.0 },
    );
    assert_point_close(
        edges[1].pos.as_ref().unwrap(),
        &TypstPoint { x: 88.0, y: 88.0 },
    );
}

#[test]
fn test_partial_iterative_layout_preserves_complement_positions() {
    let parsed = parse_dot_graphs_bytes(
        br#"digraph partial {
            a [pos="0,0"]
            b [pos="1,0"]
            c [pos="4,4"]
            a -> b [pos="0.5,0"]
            b -> c [pos="4,5"]
        }"#,
    )
    .unwrap();
    let graph = decode_graphs(&parsed).remove(0);
    let selected_label: String = decode_cbor(
        &graph_subgraph_bytes(&graph, &encode_cbor(&vec![true, true, false, false])).unwrap(),
    );

    for algo in ["force", "anneal"] {
        let laid_out = layout_parsed_graph_bytes(
            &graph,
            &encode_cbor(&BTreeMap::from([
                ("layout-algo".to_string(), algo.to_string()),
                ("subgraph".to_string(), selected_label.clone()),
                ("label-steps".to_string(), "0".to_string()),
                ("steps".to_string(), "4".to_string()),
                ("epochs".to_string(), "2".to_string()),
            ])),
        )
        .unwrap();

        let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&laid_out).unwrap());
        let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(&laid_out).unwrap());

        assert_point_close(
            nodes[2].pos.as_ref().unwrap(),
            &TypstPoint { x: 4.0, y: 4.0 },
        );
        assert_point_close(
            edges[1].pos.as_ref().unwrap(),
            &TypstPoint { x: 4.0, y: 5.0 },
        );
    }
}

#[test]
fn test_fixed_node_iterative_layout_keeps_nodes_and_complement_fixed() {
    let parsed = parse_dot_graphs_bytes(
        br#"digraph fixed_force {
            a [pos="0,0"]
            b [pos="4,0"]
            c [pos="8,0"]
            a -> b [pos="2,3"]
            b -> c [pos="88,88"]
        }"#,
    )
    .unwrap();
    let graph = decode_graphs(&parsed).remove(0);
    let selected_label: String = decode_cbor(
        &graph_subgraph_bytes(&graph, &encode_cbor(&vec![true, true, false, false])).unwrap(),
    );

    let laid_out = layout_parsed_graph_bytes(
        &graph,
        &encode_cbor(&BTreeMap::from([
            ("layout-algo".to_string(), "force".to_string()),
            ("layout-nodes".to_string(), "fixed".to_string()),
            ("subgraph".to_string(), selected_label),
            ("label-steps".to_string(), "0".to_string()),
            ("steps".to_string(), "4".to_string()),
            ("epochs".to_string(), "2".to_string()),
        ])),
    )
    .unwrap();

    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&laid_out).unwrap());
    let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(&laid_out).unwrap());

    assert_point_close(
        nodes[0].pos.as_ref().unwrap(),
        &TypstPoint { x: 0.0, y: 0.0 },
    );
    assert_point_close(
        nodes[1].pos.as_ref().unwrap(),
        &TypstPoint { x: 4.0, y: 0.0 },
    );
    assert_point_close(
        nodes[2].pos.as_ref().unwrap(),
        &TypstPoint { x: 8.0, y: 0.0 },
    );
    assert_point_close(
        edges[1].pos.as_ref().unwrap(),
        &TypstPoint { x: 88.0, y: 88.0 },
    );
    let selected_edge = edges[0].pos.as_ref().unwrap();
    assert!(
        (selected_edge.x - 2.0).abs() > 1e-6 || (selected_edge.y - 3.0).abs() > 1e-6,
        "selected edge control point should remain free when layout-nodes is fixed"
    );
}

#[test]
fn test_fixed_node_second_layout_respects_previous_node_positions() {
    let parsed = parse_dot_graphs_bytes(
        br#"digraph partial { a -> b; a -> b; a -> b; b -> c; c -> d; d -> a }"#,
    )
    .unwrap();
    let graph = decode_graphs(&parsed).remove(0);

    let tree_layout = layout_parsed_graph_bytes(
        &graph,
        &encode_cbor(&BTreeMap::from([
            ("layout-algo".to_string(), "tree".to_string()),
            ("label-steps".to_string(), "0".to_string()),
        ])),
    )
    .unwrap();
    let tree_nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&tree_layout).unwrap());

    let fixed_force_layout = layout_parsed_graph_bytes(
        &tree_layout,
        &encode_cbor(&BTreeMap::from([
            ("layout-algo".to_string(), "force".to_string()),
            ("layout-nodes".to_string(), "fixed".to_string()),
            ("label-steps".to_string(), "0".to_string()),
            ("steps".to_string(), "8".to_string()),
            ("epochs".to_string(), "3".to_string()),
        ])),
    )
    .unwrap();
    let fixed_nodes: Vec<TypstDotNode> =
        decode_cbor(&graph_nodes_bytes(&fixed_force_layout).unwrap());

    for (tree_node, fixed_node) in tree_nodes.iter().zip(fixed_nodes.iter()) {
        assert_point_close(
            fixed_node.pos.as_ref().unwrap(),
            tree_node.pos.as_ref().unwrap(),
        );
    }

    let info: TypstDotGraphInfo = decode_cbor(&graph_info_bytes(&fixed_force_layout).unwrap());
    assert!(!info.global_statements.contains_key("layout-nodes"));
    assert!(!info.global_statements.contains_key("layout-algo"));
    assert!(!info.global_statements.contains_key("steps"));
}

#[test]
fn test_fixed_node_tree_layout_straightens_existing_edges() {
    let parsed = parse_dot_graphs_bytes(
        br#"digraph straight {
            a [pos="0,0"]
            b [pos="4,0"]
            a -> b [pos="0,8"]
        }"#,
    )
    .unwrap();
    let graph = decode_graphs(&parsed).remove(0);

    let fixed_tree_layout = layout_parsed_graph_bytes(
        &graph,
        &encode_cbor(&BTreeMap::from([
            ("layout-algo".to_string(), "tree".to_string()),
            ("layout-nodes".to_string(), "fixed".to_string()),
            ("label-steps".to_string(), "0".to_string()),
        ])),
    )
    .unwrap();
    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&fixed_tree_layout).unwrap());
    let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(&fixed_tree_layout).unwrap());
    let source = nodes[0].pos.as_ref().unwrap();
    let sink = nodes[1].pos.as_ref().unwrap();

    assert_point_close(
        edges[0].pos.as_ref().unwrap(),
        &TypstPoint {
            x: 0.5 * (source.x + sink.x),
            y: 0.5 * (source.y + sink.y),
        },
    );

    let force_layout = layout_parsed_graph_bytes(
        &graph,
        &encode_cbor(&BTreeMap::from([
            ("layout-algo".to_string(), "force".to_string()),
            ("label-steps".to_string(), "0".to_string()),
            ("steps".to_string(), "8".to_string()),
            ("epochs".to_string(), "2".to_string()),
        ])),
    )
    .unwrap();
    let fixed_tree_layout = layout_parsed_graph_bytes(
        &force_layout,
        &encode_cbor(&BTreeMap::from([
            ("layout-algo".to_string(), "tree".to_string()),
            ("layout-nodes".to_string(), "fixed".to_string()),
            ("label-steps".to_string(), "0".to_string()),
        ])),
    )
    .unwrap();
    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&fixed_tree_layout).unwrap());
    let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(&fixed_tree_layout).unwrap());
    let source = nodes[0].pos.as_ref().unwrap();
    let sink = nodes[1].pos.as_ref().unwrap();

    assert_point_close(
        edges[0].pos.as_ref().unwrap(),
        &TypstPoint {
            x: 0.5 * (source.x + sink.x),
            y: 0.5 * (source.y + sink.y),
        },
    );
}

#[test]
fn test_layout_preserves_hedge_compass_subgraphs() {
    let graph = graph_from_spec_bytes(&encode_cbor(&TestGraphSpec {
        name: "compass".to_string(),
        statements: BTreeMap::new(),
        nodes: vec![
            TestNodeSpec {
                name: "a".to_string(),
                statements: BTreeMap::new(),
            },
            TestNodeSpec {
                name: "b".to_string(),
                statements: BTreeMap::new(),
            },
        ],
        edges: vec![TestEdgeSpec {
            source: Some(TestEndpointSpec {
                node: 0,
                compass: None,
                statement: None,
            }),
            sink: Some(TestEndpointSpec {
                node: 1,
                compass: Some("e".to_string()),
                statement: None,
            }),
            statements: BTreeMap::new(),
        }],
    }))
    .unwrap();

    let laid_out = layout_parsed_graph_bytes(&graph, &empty_config_bytes()).unwrap();
    let east = graph_archived_compass_subgraph_bytes(&laid_out, &encode_cbor(&"e")).unwrap();
    let hedges: Vec<usize> = decode_cbor(&subgraph_hedges_bytes(&east).unwrap());

    assert_eq!(hedges.len(), 1);
    let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(&laid_out).unwrap());
    assert_eq!(
        edges[0].sink.as_ref().unwrap().compass.as_deref(),
        Some("e")
    );
    assert_eq!(hedges[0], edges[0].sink.as_ref().unwrap().hedge);
}

#[test]
fn test_layout_handles_disconnected_spec_nodes() {
    let graph = graph_from_spec_bytes(&encode_cbor(&TestGraphSpec {
        name: "demo".to_string(),
        statements: BTreeMap::new(),
        nodes: vec![
            TestNodeSpec {
                name: "a".to_string(),
                statements: BTreeMap::new(),
            },
            TestNodeSpec {
                name: "c".to_string(),
                statements: BTreeMap::new(),
            },
            TestNodeSpec {
                name: "d".to_string(),
                statements: BTreeMap::new(),
            },
            TestNodeSpec {
                name: "e".to_string(),
                statements: BTreeMap::new(),
            },
        ],
        edges: vec![
            TestEdgeSpec {
                source: Some(TestEndpointSpec {
                    node: 0,
                    compass: None,
                    statement: None,
                }),
                sink: Some(TestEndpointSpec {
                    node: 3,
                    compass: None,
                    statement: None,
                }),
                statements: BTreeMap::new(),
            },
            TestEdgeSpec {
                source: Some(TestEndpointSpec {
                    node: 3,
                    compass: None,
                    statement: None,
                }),
                sink: Some(TestEndpointSpec {
                    node: 2,
                    compass: None,
                    statement: None,
                }),
                statements: BTreeMap::new(),
            },
        ],
    }))
    .unwrap();

    let laid_out = layout_parsed_graph_bytes(&graph, &empty_config_bytes()).unwrap();
    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&laid_out).unwrap());
    assert_eq!(nodes.len(), 4);
    assert!(nodes.iter().all(|node| node.pos.is_some()));
}

#[test]
fn test_force_center_gravity_keeps_disconnected_nodes_close() {
    fn max_node_radius(g_center: &str) -> f64 {
        let graph = graph_from_spec_bytes(&encode_cbor(&TestGraphSpec {
            name: "demo".to_string(),
            statements: BTreeMap::new(),
            nodes: vec![
                TestNodeSpec {
                    name: "a".to_string(),
                    statements: BTreeMap::new(),
                },
                TestNodeSpec {
                    name: "c".to_string(),
                    statements: BTreeMap::new(),
                },
                TestNodeSpec {
                    name: "d".to_string(),
                    statements: BTreeMap::new(),
                },
            ],
            edges: vec![TestEdgeSpec {
                source: Some(TestEndpointSpec {
                    node: 0,
                    compass: None,
                    statement: None,
                }),
                sink: Some(TestEndpointSpec {
                    node: 1,
                    compass: None,
                    statement: None,
                }),
                statements: BTreeMap::new(),
            }],
        }))
        .unwrap();

        let config = BTreeMap::from([
            ("layout-algo".to_string(), "force".to_string()),
            ("seed".to_string(), "2".to_string()),
            ("steps".to_string(), "10".to_string()),
            ("epochs".to_string(), "3".to_string()),
            ("step".to_string(), "0.81".to_string()),
            ("delta".to_string(), "0.4".to_string()),
            ("beta".to_string(), "46.1".to_string()),
            ("g-center".to_string(), g_center.to_string()),
            ("length-scale".to_string(), "0.01".to_string()),
        ]);
        let laid_out = layout_parsed_graph_bytes(&graph, &encode_cbor(&config)).unwrap();
        let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&laid_out).unwrap());

        nodes
            .iter()
            .map(|node| {
                let pos = node.pos.as_ref().unwrap();
                (pos.x.powi(2) + pos.y.powi(2)).sqrt()
            })
            .fold(0.0, f64::max)
    }

    let uncentered = max_node_radius("0");
    let centered = max_node_radius("500000");

    assert!(
        centered < uncentered,
        "expected center gravity to reduce max radius, got {centered} >= {uncentered}"
    );
}

#[test]
fn dot_cbor() {
    let figment = test_figment();
    let g = TypstGraph::from_dot(dot!(digraph{ a; a->b}).unwrap(), &figment);

    let _cbor = g.to_cbor();
}

#[test]
fn test_cbor_serialization() {
    let figment = test_figment();
    let g = TypstGraph::from_dot(dot!(digraph{ a; a->b; b->c}).unwrap(), &figment);
    let cbor = g.to_cbor();

    let test_path = "test_graph.cbor";

    cbor.serialize_to_file(test_path)
        .expect("Failed to serialize to file");

    let deserialized =
        CBORTypstGraph::deserialize_from_file(test_path).expect("Failed to deserialize from file");

    assert_eq!(cbor.nodes.len(), deserialized.nodes.len());
    assert_eq!(cbor.edges.len(), deserialized.edges.len());

    fs::remove_file(test_path).ok();
}

#[test]
fn test_typst_graph_convenience_serialization() {
    let figment = test_figment();
    let g = TypstGraph::from_dot(dot!(digraph{ a->b; b->c}).unwrap(), &figment);
    let test_path = "test_convenience_graph.cbor";

    g.serialize_to_file(test_path)
        .expect("Failed to serialize using convenience method");

    let deserialized =
        CBORTypstGraph::deserialize_from_file(test_path).expect("Failed to deserialize from file");

    assert_eq!(g.to_cbor().nodes.len(), deserialized.nodes.len());
    assert_eq!(g.to_cbor().edges.len(), deserialized.edges.len());

    fs::remove_file(test_path).ok();
}

#[test]
fn typst_graph_rkyv_roundtrip() {
    let figment = test_figment();
    let mut graph = TypstGraph::from_dot(
        dot!(digraph {
            a [pin="x:1,y:-2"];
            b;
            a -> b [pin="@edge", "label-pos"="0.4,0.2", "label-angle"="0.3rad"];
        })
        .unwrap(),
        &figment,
    );
    graph.layout();

    let bytes = rkyv::to_bytes::<_, 4096>(&graph).expect("failed to archive TypstGraph");
    let restored = unsafe {
        rkyv::from_bytes_unchecked::<TypstGraph>(&bytes).expect("failed to restore TypstGraph")
    };

    assert_eq!(
        graph.to_dot_graph().debug_dot(),
        restored.to_dot_graph().debug_dot()
    );
}

#[test]
fn test_pin_parsing() {
    let figment = test_figment();
    let mut g = TypstGraph::from_dot(
        dot!( digraph dot_80_0_GL208 {

            steps=1
            step=0.4
            beta =13.1
            "k-spring"=20.3;
            "g-center"=0
            "gamma-ee"=0.3
            "gamma-ev"=0.01
            "length-scale" = 0.25
            node[
              eval="(stroke:blue,fill :black,
         radius:2pt,
         outset: -2pt)"
            ]

            edge[
              eval=top
            ]
            v0[pin="x:@initial,y:@p1", style=invis]
            v1[pin="x:@initial,y:@p2",style=invis]
            v2[pin="x:@final,y:@p1",style=invis]
            v3[pin="x:@final,y:@p2",style=invis]
            v0 -> v11 [eval=photon]
            v1 -> v10 [eval="(..photon,label:[$gamma$],label-side: left)"]
            v9 -> v2 [eval=photon]
            v8 -> v3 [eval=photon]
            v4 -> v10
            v10 -> v5
            v5 -> v11 [dir=back]
            v11 -> v4
            v4 -> v7 [eval=gluon]
            v5 -> v6 [eval=gluon]
            v6 -> v8
            v8 -> v7
            v7 -> v9
            v9 -> v6
        })
        .unwrap(),
        &figment,
    );

    g.layout();
    println!("{}", g.to_dot_graph().debug_dot())
}

#[test]
fn test_parsing() {
    let figment = test_figment();
    let g = TypstGraph::from_dot(
        dot!(digraph qqx_aaa_pentagon {
            steps=600
            step=.2
            beta =3.1
            "k-spring"=15.3;
            "g-center"=0
            "gamma-ee"=0.3
            "gamma-ev"=0.01
            "length-scale" = 0.2

             node[
              eval="(stroke:blue,fill :black,
          radius:2pt,
          outset: -2pt)"
            ]

        exte0 [style=invis pin="x:-4"];
        v3:0 -> exte0;
        })
        .unwrap(),
        &figment,
    );

    println!("{}", g.to_dot_graph().debug_dot())
}

#[test]
fn test_directional_constraints() {
    let figment = test_figment();
    let typst_graph = TypstGraph::from_dot(
        dot!(digraph {
            "tree-dy" = 2.0
            "tree-dx" = 3.0
            a [pin="x:+@group1"]
            b [pin="x:-@group1"]
            c [pin="y:+@group2"]
            d [pin="y:-@group2"]
            e -> f [pin="x:+@edge_group"]
            g -> h [pin="x:-@edge_group"]
            a -> b
            c -> d
        })
        .unwrap(),
        &figment,
    );
    let cfg = TreeInitCfg { dy: 2.0, dx: 3.0 };
    let (node_positions, edge_positions) = typst_graph.new_positions(cfg);

    let mut group1_positive = None;
    let mut group1_negative = None;
    let mut group2_positive = None;
    let mut group2_negative = None;

    for (node_idx, _, node) in typst_graph.iter_nodes() {
        match &node.constraints.x {
            Constraint::Grouped(_, ShiftDirection::PositiveOnly) => {
                group1_positive = Some(node_positions[node_idx].x);
            }
            Constraint::Grouped(_, ShiftDirection::NegativeOnly) => {
                group1_negative = Some(node_positions[node_idx].x);
            }
            _ => {}
        }
        match &node.constraints.y {
            Constraint::Grouped(_, ShiftDirection::PositiveOnly) => {
                group2_positive = Some(node_positions[node_idx].y);
            }
            Constraint::Grouped(_, ShiftDirection::NegativeOnly) => {
                group2_negative = Some(node_positions[node_idx].y);
            }
            _ => {}
        }
    }

    if let Some(pos_x) = group1_positive {
        assert!(pos_x >= 0.0);
    }

    if let Some(neg_x) = group1_negative {
        assert!(neg_x <= 0.0);
    }

    if let Some(pos_y) = group2_positive {
        assert!(pos_y >= 0.0);
    }

    if let Some(neg_y) = group2_negative {
        assert!(neg_y <= 0.0);
    }

    if let (Some(pos_x), Some(neg_x)) = (group1_positive, group1_negative) {
        assert!(pos_x > neg_x);
    }

    if let (Some(pos_y), Some(neg_y)) = (group2_positive, group2_negative) {
        assert!(pos_y > neg_y);
    }

    let mut edge_group_positive = None;
    let mut edge_group_negative = None;

    for (_, eid, edge_data) in typst_graph.iter_edges() {
        let edge = edge_data.data;
        match &edge.constraints.x {
            Constraint::Grouped(_, ShiftDirection::PositiveOnly) => {
                edge_group_positive = Some(edge_positions[eid].x);
            }
            Constraint::Grouped(_, ShiftDirection::NegativeOnly) => {
                edge_group_negative = Some(edge_positions[eid].x);
            }
            _ => {}
        }
    }

    if let Some(pos_edge_x) = edge_group_positive {
        assert!(pos_edge_x >= 0.0);
    }

    if let Some(neg_edge_x) = edge_group_negative {
        assert!(neg_edge_x <= 0.0);
    }

    if let (Some(pos_edge_x), Some(neg_edge_x)) = (edge_group_positive, edge_group_negative) {
        assert!(pos_edge_x > neg_edge_x);
    }

    assert!(group1_positive.is_some());
    assert!(group1_negative.is_some());
    assert!(group2_positive.is_some());
    assert!(group2_negative.is_some());
    assert!(edge_group_positive.is_some());
    assert!(edge_group_negative.is_some());
}

#[test]
fn grouped_directional_constraints_enforce_pin_side_after_layout() {
    let figment = test_figment();
    let typst_graph = TypstGraph::from_dot(
        dot!(digraph {
            a -> b [pin="x:@+right"]
            b -> c [pin="x:@-left"]
        })
        .unwrap(),
        &figment,
    );
    let (mut node_positions, mut edge_positions) =
        typst_graph.new_positions(TreeInitCfg { dx: 1.0, dy: 1.0 });

    edge_positions[EdgeIndex(0)].x = -2.0;
    edge_positions[EdgeIndex(1)].x = 2.0;

    typst_graph.apply_grouped_constraints(&mut node_positions, &mut edge_positions);

    assert!(edge_positions[EdgeIndex(0)].x >= 0.0);
    assert!(edge_positions[EdgeIndex(1)].x <= 0.0);
}

#[test]
fn gammaloop_external_edge_pins_align_by_group_and_side() {
    let figment = test_figment();
    let mut typst_graph = TypstGraph::from_dot(
        dot!(digraph {
            "layout-algo" = "force"
            steps = 10
            epochs = 2
            ext_r0 [style=invis]
            ext_r1 [style=invis]
            ext_l0 [style=invis]
            ext_l1 [style=invis]
            a
            b
            ext_r0 -> a [id=0 pin="x:@+right,y:@edgee0"]
            ext_r1 -> a [id=1 pin="x:@+right,y:@edgee1"]
            b -> ext_l0 [id=2 pin="x:@-left,y:@edgee0"]
            b -> ext_l1 [id=3 pin="x:@-left,y:@edgee1"]
            a -> b
        })
        .unwrap(),
        &figment,
    );

    typst_graph.layout();

    let mut edge_positions = BTreeMap::new();
    for (_, edge_index, edge_data) in typst_graph.iter_edges() {
        let edge = edge_data.data;
        edge_positions.insert(edge_index.0, edge.pos);
    }

    let right0 = edge_positions.get(&0).unwrap();
    let right1 = edge_positions.get(&1).unwrap();
    let left0 = edge_positions.get(&2).unwrap();
    let left1 = edge_positions.get(&3).unwrap();

    assert!(right0.x > 0.0, "{right0:?}");
    assert!(right1.x > 0.0, "{right1:?}");
    assert!(left0.x < 0.0, "{left0:?}");
    assert!(left1.x < 0.0, "{left1:?}");
    assert!((right0.y - left0.y).abs() < 1e-9, "{right0:?} {left0:?}");
    assert!((right1.y - left1.y).abs() < 1e-9, "{right1:?} {left1:?}");
}

#[test]
fn gammaloop_external_edge_pos_aligns_by_group_and_side() {
    let figment = test_figment();
    let mut typst_graph = TypstGraph::from_dot(
        dot!(digraph {
            "layout-algo" = "force"
            steps = 10
            epochs = 2
            ext_l0 [style=invis]
            ext_l1 [style=invis]
            ext_r0 [style=invis]
            ext_r1 [style=invis]
            a
            b
            ext_l0 -> a [id=0 pos="x:@-left!,y:@edgee0!"]
            ext_l1 -> a [id=1 pos="x:@-left!,y:@edgee1!"]
            b -> ext_r0 [id=2 pos="x:@+right!,y:@edgee0!"]
            b -> ext_r1 [id=3 pos="x:@+right!,y:@edgee1!"]
            a -> b
        })
        .unwrap(),
        &figment,
    );

    typst_graph.layout();

    let mut edge_positions = BTreeMap::new();
    for (_, edge_index, edge_data) in typst_graph.iter_edges() {
        let edge = edge_data.data;
        edge_positions.insert(edge_index.0, edge.pos);
    }

    let left0 = edge_positions.get(&0).unwrap();
    let left1 = edge_positions.get(&1).unwrap();
    let right0 = edge_positions.get(&2).unwrap();
    let right1 = edge_positions.get(&3).unwrap();

    assert!(left0.x < 0.0, "{left0:?}");
    assert!(left1.x < 0.0, "{left1:?}");
    assert!(right0.x > 0.0, "{right0:?}");
    assert!(right1.x > 0.0, "{right1:?}");
    assert!((left0.y - right0.y).abs() < 1e-9, "{left0:?} {right0:?}");
    assert!((left1.y - right1.y).abs() < 1e-9, "{left1:?} {right1:?}");
}

#[test]
fn test_full() {
    let input = stringify!(
       digraph {
           steps=600
           step=0.4
           beta =13.1
           "k-spring"=3.3;
           "g-center"=0
           "gamma-dangling"=50
           "gamma-ee"=0.3
           "gamma-ev"=0.01
           "length-scale" = 0.25
           a->b
       }
    );

    let figment = Figment::from(Serialized::from(
        BTreeMap::<String, String>::new(),
        Profile::Default,
    ));

    let dots = DotGraphSet::from_string_with_figment(input, figment.clone())
        .map_err(|a| a.to_string())
        .unwrap()
        .into_iter();

    let mut graphs = Vec::new();
    for g in dots {
        let mut typst_graph = TypstGraph::from_dot(g, &figment);

        typst_graph.layout();
        graphs.push((
            typst_graph.to_cbor(),
            typst_graph.to_dot_graph().debug_dot(),
        ));
    }
}
