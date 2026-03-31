use std::collections::BTreeMap;

use dot_parser::ast::CompassPt;
use linnet::{
    half_edge::{involution::ArchivedOrientation, subgraph::SubSetLike},
    parser::{
        ArchivedDotEdgeView, ArchivedDotEndpointView, ArchivedDotGraphView, ArchivedDotVertexView,
        DotGraph,
    },
};
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct TypstDotGraphInfo {
    pub name: String,
    pub global_statements: BTreeMap<String, String>,
    pub edge_statements: BTreeMap<String, String>,
    pub node_statements: BTreeMap<String, String>,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct TypstPoint {
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct TypstDotNode {
    pub node: usize,
    pub name: Option<String>,
    pub index: Option<usize>,
    pub pos: Option<TypstPoint>,
    pub shift: Option<TypstPoint>,
    pub eval: Option<String>,
    pub statements: BTreeMap<String, String>,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct TypstDotEndpoint {
    pub node: usize,
    pub hedge: usize,
    pub statement: Option<String>,
    pub port_label: Option<String>,
    pub compass: Option<String>,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct TypstDotEdge {
    pub edge: usize,
    pub orientation: String,
    pub source: Option<TypstDotEndpoint>,
    pub sink: Option<TypstDotEndpoint>,
    pub pos: Option<TypstPoint>,
    pub shift: Option<TypstPoint>,
    pub label_pos: Option<TypstPoint>,
    pub label_angle: Option<f64>,
    pub bend: Option<f64>,
    pub eval_sink: Option<String>,
    pub eval_source: Option<String>,
    pub eval_label: Option<String>,
    pub mom_eval: Option<String>,
    pub statements: BTreeMap<String, String>,
}

pub fn encode_cbor<T: Serialize>(value: &T) -> Result<Vec<u8>, String> {
    let mut buffer = Vec::new();
    ciborium::ser::into_writer(value, &mut buffer).map_err(|err| err.to_string())?;
    Ok(buffer)
}

pub fn decode_graph_bytes_list(arg: &[u8]) -> Result<Vec<Vec<u8>>, String> {
    ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize graph bytes list: {err}"))
}

pub fn graph_info_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let graph = DotGraph::archived_view(arg);
    encode_cbor(&graph_info(&graph))
}

pub fn graph_nodes_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let graph = DotGraph::archived_view(arg);
    encode_cbor(
        &graph
            .vertex_data()
            .map(node_view_to_output)
            .collect::<Vec<_>>(),
    )
}

pub fn graph_nodes_of_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = DotGraph::archived_view(arg);
    let subgraph = decode_subgraph_spec(&graph, arg2)?;
    encode_cbor(
        &graph
            .vertex_data_of(&subgraph)
            .map(node_view_to_output)
            .collect::<Vec<_>>(),
    )
}

pub fn graph_edges_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let graph = DotGraph::archived_view(arg);
    encode_cbor(
        &graph
            .edge_data()
            .map(|edge| edge_view_to_output(&graph, edge))
            .collect::<Vec<_>>(),
    )
}

pub fn graph_edges_of_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = DotGraph::archived_view(arg);
    let subgraph = decode_subgraph_spec(&graph, arg2)?;
    encode_cbor(
        &graph
            .edge_data_of(&subgraph)
            .map(|edge| edge_view_to_output(&graph, edge))
            .collect::<Vec<_>>(),
    )
}

pub fn graph_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = DotGraph::archived_view(arg);
    let subgraph = decode_subgraph_spec(&graph, arg2)?;
    encode_cbor(&subgraph.string_label())
}

pub fn graph_compass_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = DotGraph::archived_view(arg);
    let compass = decode_compass(arg2)?;
    encode_cbor(&graph.compass_subgraph(compass).string_label())
}

fn graph_info(graph: &ArchivedDotGraphView<'_>) -> TypstDotGraphInfo {
    let global_data = graph.global_data();
    TypstDotGraphInfo {
        name: global_data.name.as_str().to_string(),
        global_statements: global_data
            .statements
            .iter()
            .map(|(key, value)| (key.as_str().to_string(), value.as_str().to_string()))
            .collect(),
        edge_statements: global_data
            .edge_statements
            .iter()
            .map(|(key, value)| (key.as_str().to_string(), value.as_str().to_string()))
            .collect(),
        node_statements: global_data
            .node_statements
            .iter()
            .map(|(key, value)| (key.as_str().to_string(), value.as_str().to_string()))
            .collect(),
    }
}

fn node_view_to_output(vertex: ArchivedDotVertexView<'_>) -> TypstDotNode {
    let statements = vertex
        .data
        .statements
        .iter()
        .map(|(key, value)| (key.as_str().to_string(), value.as_str().to_string()))
        .collect::<BTreeMap<_, _>>();

    TypstDotNode {
        node: vertex.node.0,
        name: vertex
            .data
            .name
            .as_ref()
            .map(|value| value.as_str().to_string()),
        index: vertex
            .data
            .index
            .as_ref()
            .map(|value| value.0.try_into().unwrap()),
        pos: parse_point(&statements, "pos"),
        shift: parse_point(&statements, "shift"),
        eval: statements.get("eval").cloned(),
        statements,
    }
}

fn edge_view_to_output(
    graph: &ArchivedDotGraphView<'_>,
    edge: ArchivedDotEdgeView<'_>,
) -> TypstDotEdge {
    let statements = edge
        .data
        .statements
        .iter()
        .map(|(key, value)| (key.as_str().to_string(), value.as_str().to_string()))
        .collect::<BTreeMap<_, _>>();
    let endpoints = graph.endpoints_of_edge(&edge);

    TypstDotEdge {
        edge: edge.edge.0,
        orientation: orientation_to_string(edge.orientation).to_string(),
        source: endpoints.source.map(endpoint_to_output),
        sink: endpoints.sink.map(endpoint_to_output),
        pos: parse_point(&statements, "pos"),
        shift: parse_point(&statements, "shift"),
        label_pos: parse_point(&statements, "label_pos"),
        label_angle: parse_rad(&statements, "label_angle"),
        bend: parse_rad(&statements, "bend"),
        eval_sink: statements.get("eval_sink").cloned(),
        eval_source: statements.get("eval_source").cloned(),
        eval_label: statements.get("eval_label").cloned(),
        mom_eval: statements.get("mom_eval").cloned(),
        statements,
    }
}

fn endpoint_to_output(endpoint: ArchivedDotEndpointView<'_>) -> TypstDotEndpoint {
    TypstDotEndpoint {
        node: endpoint.node.0,
        hedge: endpoint.hedge.0,
        statement: endpoint
            .data
            .statement
            .as_ref()
            .map(|value| value.as_str().to_string()),
        port_label: endpoint
            .data
            .port_label
            .as_ref()
            .map(|value| value.as_str().to_string()),
        compass: endpoint
            .data
            .compasspt
            .as_ref()
            .copied()
            .map(compass_to_string),
    }
}

fn decode_subgraph_spec(
    graph: &ArchivedDotGraphView<'_>,
    arg: &[u8],
) -> Result<linnet::half_edge::subgraph::SuBitGraph, String> {
    if let Ok(bits) = ciborium::de::from_reader::<Vec<bool>, _>(arg) {
        return graph.subgraph_from_bools(bits);
    }

    if let Ok(label) = ciborium::de::from_reader::<String, _>(arg) {
        return graph.subgraph_from_base62(&label);
    }

    let label = std::str::from_utf8(arg).map_err(|_| {
        "Subgraph must be provided as a CBOR-encoded bool list, a CBOR string, or UTF-8 bytes"
            .to_string()
    })?;
    graph.subgraph_from_base62(label)
}

fn decode_compass(arg: &[u8]) -> Result<Option<CompassPt>, String> {
    if let Ok(value) = ciborium::de::from_reader::<Option<String>, _>(arg) {
        let compass = match value {
            Some(value) => parse_compass(&value)?,
            None => None,
        };
        return Ok(compass);
    }

    if let Ok(value) = ciborium::de::from_reader::<String, _>(arg) {
        return parse_compass(&value);
    }

    let value = std::str::from_utf8(arg)
        .map_err(|_| "Compass must be provided as a CBOR string, CBOR none, or UTF-8 bytes")?;
    parse_compass(value)
}

fn parse_point(statements: &BTreeMap<String, String>, key: &str) -> Option<TypstPoint> {
    let value = statements.get(key)?;
    let cleaned = value
        .trim()
        .trim_matches('"')
        .trim_matches(|c| c == '(' || c == ')');
    let parts: Vec<_> = cleaned
        .split([',', ' '])
        .filter(|part| !part.is_empty())
        .collect();

    if parts.len() != 2 {
        return None;
    }

    let x = parts[0].trim().parse::<f64>().ok()?;
    let y = parts[1].trim().parse::<f64>().ok()?;
    Some(TypstPoint { x, y })
}

fn parse_rad(statements: &BTreeMap<String, String>, key: &str) -> Option<f64> {
    statements.get(key).and_then(|value| {
        value
            .trim()
            .trim_matches('"')
            .trim_end_matches("rad")
            .trim()
            .parse::<f64>()
            .ok()
    })
}

fn parse_compass(value: &str) -> Result<Option<CompassPt>, String> {
    match value.trim().to_ascii_lowercase().as_str() {
        "" | "none" => Ok(None),
        "n" => Ok(Some(CompassPt::N)),
        "ne" => Ok(Some(CompassPt::NE)),
        "e" => Ok(Some(CompassPt::E)),
        "se" => Ok(Some(CompassPt::SE)),
        "s" => Ok(Some(CompassPt::S)),
        "sw" => Ok(Some(CompassPt::SW)),
        "w" => Ok(Some(CompassPt::W)),
        "nw" => Ok(Some(CompassPt::NW)),
        "c" => Ok(Some(CompassPt::C)),
        other => Err(format!("Invalid compass point: {other}")),
    }
}

fn orientation_to_string(orientation: &ArchivedOrientation) -> &'static str {
    match orientation {
        ArchivedOrientation::Default => "default",
        ArchivedOrientation::Reversed => "reversed",
        ArchivedOrientation::Undirected => "undirected",
    }
}

fn compass_to_string(compass: u8) -> String {
    match compass {
        0 => "n",
        1 => "ne",
        2 => "e",
        3 => "se",
        4 => "s",
        5 => "sw",
        6 => "w",
        7 => "nw",
        8 => "c",
        _ => "?",
    }
    .to_string()
}
