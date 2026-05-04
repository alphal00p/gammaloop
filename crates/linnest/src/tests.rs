use std::{collections::BTreeMap, fs};

use figment::providers::Serialized;
use figment::{Figment, Profile};
use linnet::half_edge::layout::spring::{Constraint, ShiftDirection};
use linnet::half_edge::swap::Swap;
use linnet::{
    dot,
    parser::{set::DotGraphSet, DotGraph},
};
use serde::{Deserialize, Serialize};

use crate::{
    builder_add_edge_bytes, builder_add_node_bytes, builder_finish_bytes, builder_new_bytes,
    expand_template, graph_archived_compass_subgraph_bytes, graph_compass_subgraph_bytes,
    graph_cycle_basis_bytes, graph_dot_bytes, graph_edges_bytes,
    graph_edges_of_archived_subgraph_bytes, graph_edges_of_bytes, graph_from_spec_bytes,
    graph_info_bytes, graph_nodes_bytes, graph_nodes_of_archived_subgraph_bytes,
    graph_nodes_of_bytes, graph_spanning_forests_bytes, graph_subgraph_bytes, layout_graph_bytes,
    layout_parsed_graph_bytes, layout_parsed_graphs_bytes, parse_dot_graphs_bytes,
    subgraph_contains_hedge_bytes, subgraph_hedges_bytes, subgraph_label_bytes, CBORTypstGraph,
    TreeInitCfg, TypstDotEdge, TypstDotGraphInfo, TypstDotNode, TypstGraph,
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
struct TestTemplatedGraphSpec {
    name: String,
    #[serde(default)]
    statements: BTreeMap<String, String>,
    #[serde(default)]
    edge_statements: BTreeMap<String, String>,
    #[serde(default)]
    node_statements: BTreeMap<String, String>,
    nodes: Vec<TestNodeSpec>,
    edges: Vec<TestEdgeSpec>,
}

#[derive(Serialize)]
struct TestBuilderSpec {
    name: String,
    #[serde(default)]
    statements: BTreeMap<String, String>,
    #[serde(default)]
    edge_statements: BTreeMap<String, String>,
    #[serde(default)]
    node_statements: BTreeMap<String, String>,
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
struct TestEndpointSpec {
    node: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    compass: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    statement: Option<String>,
}

#[derive(Deserialize)]
struct TestBuilderNodeResult {
    builder: Vec<u8>,
    node: usize,
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
fn test_graph_spec_constructor_expands_default_statement_templates() {
    let spec = TestTemplatedGraphSpec {
        name: "templated".to_string(),
        statements: BTreeMap::new(),
        edge_statements: BTreeMap::from([(
            "eval_label".to_string(),
            "(text(fill: rgb(\"#{color}\"))[{label}])".to_string(),
        )]),
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
        edges[0].eval_label.as_deref(),
        Some("(text(fill: rgb(\"#0055ff\"))[a-c])")
    );
    assert_eq!(
        edges[0].statements.get("eval_label").map(String::as_str),
        Some("(text(fill: rgb(\"#0055ff\"))[a-c])")
    );
}

#[test]
fn test_archived_builder_and_subgraph_api() {
    let builder = builder_new_bytes(&encode_cbor(&BTreeMap::from([(
        "name".to_string(),
        "builder".to_string(),
    )])))
    .unwrap();

    let a: TestBuilderNodeResult = decode_cbor(
        &builder_add_node_bytes(
            &builder,
            &encode_cbor(&TestNodeSpec {
                name: "a".to_string(),
                statements: BTreeMap::new(),
            }),
        )
        .unwrap(),
    );
    let b: TestBuilderNodeResult = decode_cbor(
        &builder_add_node_bytes(
            &a.builder,
            &encode_cbor(&TestNodeSpec {
                name: "b".to_string(),
                statements: BTreeMap::new(),
            }),
        )
        .unwrap(),
    );

    let builder = builder_add_edge_bytes(
        &b.builder,
        &encode_cbor(&TestEdgeSpec {
            source: Some(TestEndpointSpec {
                node: a.node,
                compass: Some("e".to_string()),
                statement: None,
            }),
            sink: Some(TestEndpointSpec {
                node: b.node,
                compass: Some("w".to_string()),
                statement: None,
            }),
            statements: one_statement("join", "ab"),
        }),
    )
    .unwrap();
    let graph = builder_finish_bytes(&builder).unwrap();

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
fn test_builder_pattern_expands_default_statement_templates() {
    let builder = builder_new_bytes(&encode_cbor(&TestBuilderSpec {
        name: "builder".to_string(),
        statements: BTreeMap::new(),
        edge_statements: BTreeMap::from([(
            "eval_label".to_string(),
            "(text(fill: rgb(\"#{color}\"))[{label}])".to_string(),
        )]),
        node_statements: BTreeMap::from([(
            "eval".to_string(),
            "(fill: rgb(\"#{color}\"))".to_string(),
        )]),
    }))
    .unwrap();

    let a: TestBuilderNodeResult = decode_cbor(
        &builder_add_node_bytes(
            &builder,
            &encode_cbor(&TestNodeSpec {
                name: "a".to_string(),
                statements: one_statement("color", "ff0000"),
            }),
        )
        .unwrap(),
    );
    let b: TestBuilderNodeResult = decode_cbor(
        &builder_add_node_bytes(
            &a.builder,
            &encode_cbor(&TestNodeSpec {
                name: "b".to_string(),
                statements: one_statement("color", "00aa00"),
            }),
        )
        .unwrap(),
    );

    let builder = builder_add_edge_bytes(
        &b.builder,
        &encode_cbor(&TestEdgeSpec {
            source: Some(TestEndpointSpec {
                node: a.node,
                compass: None,
                statement: None,
            }),
            sink: Some(TestEndpointSpec {
                node: b.node,
                compass: None,
                statement: None,
            }),
            statements: BTreeMap::from([
                ("color".to_string(), "0055ff".to_string()),
                ("label".to_string(), "a-c".to_string()),
            ]),
        }),
    )
    .unwrap();
    let graph = builder_finish_bytes(&builder).unwrap();

    let nodes: Vec<TypstDotNode> = decode_cbor(&graph_nodes_bytes(&graph).unwrap());
    assert_eq!(nodes[0].eval.as_deref(), Some("(fill: rgb(\"#ff0000\"))"));
    let edges: Vec<TypstDotEdge> = decode_cbor(&graph_edges_bytes(&graph).unwrap());
    assert_eq!(
        edges[0].eval_label.as_deref(),
        Some("(text(fill: rgb(\"#0055ff\"))[a-c])")
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
fn test_layout_handles_disconnected_builder_nodes() {
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
            ("layout_algo".to_string(), "force".to_string()),
            ("seed".to_string(), "2".to_string()),
            ("steps".to_string(), "10".to_string()),
            ("epochs".to_string(), "3".to_string()),
            ("step".to_string(), "0.81".to_string()),
            ("delta".to_string(), "0.4".to_string()),
            ("beta".to_string(), "46.1".to_string()),
            ("g_center".to_string(), g_center.to_string()),
            ("length_scale".to_string(), "0.01".to_string()),
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
fn test_pin_parsing() {
    let figment = test_figment();
    let mut g = TypstGraph::from_dot(dot!( digraph dot_80_0_GL208 {

       steps=1
       step=0.4
       beta =13.1
       k_spring=20.3;
       g_center=0
       gamma_ee=0.3
       gamma_ev=0.01
       length_scale = 0.25
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
       v1 -> v10 [eval="(..photon,label:[$gamma$],label-side: left)", mom_eval="(label:[$p_1$],label-sep:0mm)"]
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
   }).unwrap(), &figment);

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
            k_spring=15.3;
            g_center=0
            gamma_ee=0.3
            gamma_ev=0.01
            length_scale = 0.2

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
            tree_dy = 2.0
            tree_dx = 3.0
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
}

#[test]
fn test_full() {
    let input = stringify!(
       digraph {
           steps=600
           step=0.4
           beta =13.1
           k_spring=3.3;
           g_center=0
           gamma_dangling=50
           gamma_ee=0.3
           gamma_ev=0.01
           length_scale = 0.25
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
