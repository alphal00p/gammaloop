use std::collections::BTreeMap;

use dot_parser::ast::CompassPt;
use linnet::{
    half_edge::{
        builder::{HedgeData, HedgeGraphBuilder},
        involution::{ArchivedOrientation, EdgeData, Flow, Hedge, Orientation},
        nodestore::DefaultNodeStore,
        subgraph::{Inclusion, SuBitGraph, SubSetLike},
        NodeIndex,
    },
    parser::{
        ArchivedDotEdgeView, ArchivedDotEndpointView, ArchivedDotGraphView, ArchivedDotVertexView,
        DotEdgeData, DotGraph, DotHedgeData, DotVertexData, GlobalData,
    },
};
use serde::{Deserialize, Serialize};

type DotBuilder = HedgeGraphBuilder<DotEdgeData, DotVertexData, DotHedgeData>;

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

#[derive(Debug, Deserialize)]
pub struct TypstGraphSpec {
    #[serde(default)]
    pub name: Option<String>,
    #[serde(default)]
    pub statements: BTreeMap<String, String>,
    #[serde(default)]
    pub edge_statements: BTreeMap<String, String>,
    #[serde(default)]
    pub node_statements: BTreeMap<String, String>,
    #[serde(default)]
    pub nodes: Vec<TypstNodeSpec>,
    #[serde(default)]
    pub edges: Vec<TypstEdgeSpec>,
}

#[derive(Debug, Deserialize)]
pub struct TypstNodeSpec {
    #[serde(default)]
    pub name: Option<String>,
    #[serde(default)]
    pub index: Option<usize>,
    #[serde(default)]
    pub statements: BTreeMap<String, String>,
}

#[derive(Debug, Deserialize)]
pub struct TypstEdgeSpec {
    #[serde(default)]
    pub source: Option<TypstEndpointSpec>,
    #[serde(default)]
    pub sink: Option<TypstEndpointSpec>,
    #[serde(default)]
    pub orientation: Option<String>,
    #[serde(default)]
    pub flow: Option<String>,
    #[serde(default)]
    pub id: Option<usize>,
    #[serde(default)]
    pub statements: BTreeMap<String, String>,
}

#[derive(Debug, Deserialize)]
pub struct TypstEndpointSpec {
    pub node: usize,
    #[serde(default)]
    pub statement: Option<String>,
    #[serde(default)]
    pub id: Option<usize>,
    #[serde(default)]
    pub port_label: Option<String>,
    #[serde(default)]
    pub compass: Option<String>,
    #[serde(default)]
    pub in_subgraph: bool,
}

#[derive(Debug, Clone, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)]
pub struct TypstDotGraphBuilder {
    pub global_data: GlobalData,
    pub builder: DotBuilder,
    pub node_count: usize,
}

#[derive(Debug, Deserialize)]
pub struct TypstBuilderSpec {
    #[serde(default)]
    pub name: Option<String>,
    #[serde(default)]
    pub statements: BTreeMap<String, String>,
    #[serde(default)]
    pub edge_statements: BTreeMap<String, String>,
    #[serde(default)]
    pub node_statements: BTreeMap<String, String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct TypstBuilderNodeResult {
    pub builder: Vec<u8>,
    pub node: usize,
}

#[derive(Debug, Deserialize)]
pub struct TypstJoinSpec {
    pub key: String,
}

pub fn encode_cbor<T: Serialize>(value: &T) -> Result<Vec<u8>, String> {
    let mut buffer = Vec::new();
    ciborium::ser::into_writer(value, &mut buffer).map_err(|err| err.to_string())?;
    Ok(buffer)
}

fn to_rkyv_bytes<T, const N: usize>(value: &T) -> Result<Vec<u8>, String>
where
    T: rkyv::Serialize<rkyv::ser::serializers::AllocSerializer<N>>,
{
    rkyv::to_bytes::<_, N>(value)
        .map(|bytes| bytes.to_vec())
        .map_err(|err| err.to_string())
}

fn decode_cbor<T: for<'de> Deserialize<'de>>(arg: &[u8], name: &str) -> Result<T, String> {
    ciborium::de::from_reader(arg).map_err(|err| format!("Failed to deserialize {name}: {err}"))
}

fn decode_builder(arg: &[u8]) -> Result<TypstDotGraphBuilder, String> {
    unsafe { rkyv::from_bytes_unchecked::<TypstDotGraphBuilder>(arg) }
        .map_err(|err| format!("Failed to deserialize archived graph builder: {err}"))
}

fn encode_builder(builder: &TypstDotGraphBuilder) -> Result<Vec<u8>, String> {
    to_rkyv_bytes::<_, 4096>(builder)
}

fn decode_subgraph(arg: &[u8]) -> Result<SuBitGraph, String> {
    unsafe { rkyv::from_bytes_unchecked::<SuBitGraph>(arg) }
        .map_err(|err| format!("Failed to deserialize archived subgraph: {err}"))
}

fn encode_subgraph(subgraph: &SuBitGraph) -> Result<Vec<u8>, String> {
    to_rkyv_bytes::<_, 256>(subgraph)
}

pub fn builder_new_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec = if arg.is_empty() {
        TypstBuilderSpec {
            name: None,
            statements: BTreeMap::new(),
            edge_statements: BTreeMap::new(),
            node_statements: BTreeMap::new(),
        }
    } else {
        decode_cbor(arg, "builder spec")?
    };

    encode_builder(&TypstDotGraphBuilder {
        global_data: GlobalData {
            name: spec.name.unwrap_or_else(|| "constructed".to_string()),
            statements: spec.statements,
            edge_statements: spec.edge_statements,
            node_statements: spec.node_statements,
        },
        builder: HedgeGraphBuilder::new(),
        node_count: 0,
    })
}

pub fn builder_add_node_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let mut graph_builder = decode_builder(arg)?;
    let node: TypstNodeSpec = decode_cbor(arg2, "node spec")?;
    let node_index = graph_builder.builder.add_node(DotVertexData {
        name: node.name,
        index: node
            .index
            .map(NodeIndex)
            .or(Some(NodeIndex(graph_builder.node_count))),
        statements: node.statements,
    });
    graph_builder.node_count += 1;

    encode_cbor(&TypstBuilderNodeResult {
        builder: encode_builder(&graph_builder)?,
        node: node_index.0,
    })
}

pub fn builder_add_edge_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let mut graph_builder = decode_builder(arg)?;
    let edge: TypstEdgeSpec = decode_cbor(arg2, "edge spec")?;
    add_edge_to_builder(&mut graph_builder.builder, graph_builder.node_count, edge)?;
    encode_builder(&graph_builder)
}

pub fn builder_finish_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let graph_builder = decode_builder(arg)?;
    let graph = DotGraph {
        global_data: graph_builder.global_data,
        graph: graph_builder
            .builder
            .build::<DefaultNodeStore<DotVertexData>>(),
    };
    graph.to_rkyv_bytes::<4096>().map(|bytes| bytes.to_vec())
}

pub fn graph_from_spec_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: TypstGraphSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize graph spec: {err}"))?;
    let builder = builder_from_spec(spec)?;
    builder_finish_bytes(&encode_builder(&builder)?)
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

pub fn graph_nodes_of_archived_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = DotGraph::archived_view(arg);
    let subgraph = decode_subgraph(arg2)?;
    encode_cbor(
        &graph
            .vertex_data_of(&subgraph)
            .map(node_view_to_output)
            .collect::<Vec<_>>(),
    )
}

pub fn graph_edges_of_archived_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = DotGraph::archived_view(arg);
    let subgraph = decode_subgraph(arg2)?;
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

pub fn graph_archived_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = DotGraph::archived_view(arg);
    let subgraph = decode_subgraph_spec(&graph, arg2)?;
    encode_subgraph(&subgraph)
}

pub fn graph_compass_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = DotGraph::archived_view(arg);
    let compass = decode_compass(arg2)?;
    encode_cbor(&graph.compass_subgraph(compass).string_label())
}

pub fn graph_archived_compass_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = DotGraph::archived_view(arg);
    let compass = decode_compass(arg2)?;
    encode_subgraph(&graph.compass_subgraph(compass))
}

pub fn subgraph_label_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let subgraph = decode_subgraph(arg)?;
    encode_cbor(&subgraph.string_label())
}

pub fn subgraph_hedges_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let subgraph = decode_subgraph(arg)?;
    encode_cbor(
        &subgraph
            .included_iter()
            .map(|hedge| hedge.0)
            .collect::<Vec<_>>(),
    )
}

pub fn subgraph_contains_hedge_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let subgraph = decode_subgraph(arg)?;
    let hedge: usize = decode_cbor(arg2, "hedge index")?;
    encode_cbor(&subgraph.includes(&Hedge(hedge)))
}

pub fn graph_cycle_basis_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let graph: DotGraph = unsafe { rkyv::from_bytes_unchecked(arg) }
        .map_err(|err| format!("Failed to deserialize archived dot graph: {err}"))?;
    let (cycles, _) = graph.cycle_basis();
    let archived = cycles
        .into_iter()
        .map(|cycle| encode_subgraph(&cycle.filter))
        .collect::<Result<Vec<_>, _>>()?;
    encode_cbor(&archived)
}

pub fn graph_spanning_forests_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let graph: DotGraph = unsafe { rkyv::from_bytes_unchecked(arg) }
        .map_err(|err| format!("Failed to deserialize archived dot graph: {err}"))?;
    let forests = graph.all_spanning_forests_of(&graph.full_filter());
    let archived = forests
        .into_iter()
        .map(|forest| encode_subgraph(&forest))
        .collect::<Result<Vec<_>, _>>()?;
    encode_cbor(&archived)
}

pub fn graph_join_by_edge_key_bytes(
    arg: &[u8],
    arg2: &[u8],
    arg3: &[u8],
) -> Result<Vec<u8>, String> {
    let left: DotGraph = unsafe { rkyv::from_bytes_unchecked(arg) }
        .map_err(|err| format!("Failed to deserialize left archived dot graph: {err}"))?;
    let right: DotGraph = unsafe { rkyv::from_bytes_unchecked(arg2) }
        .map_err(|err| format!("Failed to deserialize right archived dot graph: {err}"))?;
    let spec: TypstJoinSpec = decode_cbor(arg3, "join spec")?;
    let key = spec.key;

    let global_data = left.global_data.clone();
    let graph = left
        .graph
        .join(
            right.graph,
            |left_flow, left_data, right_flow, right_data| {
                left_flow == -right_flow
                    && statement_value(left_data, &key)
                        .zip(statement_value(right_data, &key))
                        .is_some_and(|(left, right)| left == right)
            },
            |left_flow, left_data, _, _| (left_flow, left_data),
        )
        .map_err(|err| err.to_string())?;

    DotGraph { global_data, graph }
        .to_rkyv_bytes::<4096>()
        .map(|bytes| bytes.to_vec())
}

pub fn graph_join_by_hedge_key_bytes(
    arg: &[u8],
    arg2: &[u8],
    arg3: &[u8],
) -> Result<Vec<u8>, String> {
    let left: DotGraph = unsafe { rkyv::from_bytes_unchecked(arg) }
        .map_err(|err| format!("Failed to deserialize left archived dot graph: {err}"))?;
    let right: DotGraph = unsafe { rkyv::from_bytes_unchecked(arg2) }
        .map_err(|err| format!("Failed to deserialize right archived dot graph: {err}"))?;
    let spec: TypstJoinSpec = decode_cbor(arg3, "join spec")?;
    let key = spec.key;

    let global_data = left.global_data.clone();
    let graph = left
        .graph
        .join_with_hedge_data(
            right.graph,
            |_, left_flow, _, left_hedge, _, right_flow, _, right_hedge| {
                left_flow == -right_flow
                    && hedge_value(left_hedge, &key)
                        .zip(hedge_value(right_hedge, &key))
                        .is_some_and(|(left, right)| left == right)
            },
            |left_flow, left_data, _, _| (left_flow, left_data),
        )
        .map_err(|err| err.to_string())?;

    DotGraph { global_data, graph }
        .to_rkyv_bytes::<4096>()
        .map(|bytes| bytes.to_vec())
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

fn builder_from_spec(spec: TypstGraphSpec) -> Result<TypstDotGraphBuilder, String> {
    let node_count = spec.nodes.len();
    let mut builder = HedgeGraphBuilder::<DotEdgeData, DotVertexData, DotHedgeData>::new();

    for (node_index, node) in spec.nodes.into_iter().enumerate() {
        builder.add_node(DotVertexData {
            name: node.name,
            index: node.index.map(NodeIndex).or(Some(NodeIndex(node_index))),
            statements: node.statements,
        });
    }

    for edge in spec.edges {
        add_edge_to_builder(&mut builder, node_count, edge)?;
    }

    Ok(TypstDotGraphBuilder {
        global_data: GlobalData {
            name: spec.name.unwrap_or_else(|| "constructed".to_string()),
            statements: spec.statements,
            edge_statements: spec.edge_statements,
            node_statements: spec.node_statements,
        },
        builder,
        node_count,
    })
}

fn add_edge_to_builder(
    builder: &mut DotBuilder,
    node_count: usize,
    edge: TypstEdgeSpec,
) -> Result<(), String> {
    let orientation = edge
        .orientation
        .as_deref()
        .map(parse_orientation)
        .transpose()?
        .unwrap_or(Orientation::Default);
    let edge_data = DotEdgeData {
        statements: edge.statements,
        local_statements: BTreeMap::new(),
        edge_id: edge.id.map(linnet::half_edge::involution::EdgeIndex::from),
    };

    match (edge.source, edge.sink) {
        (Some(source), Some(sink)) => {
            let source = endpoint_spec_to_hedge_data(source, node_count)?;
            let sink = endpoint_spec_to_hedge_data(sink, node_count)?;
            builder.add_edge(source, sink, edge_data, orientation);
        }
        (Some(source), None) => {
            let source = endpoint_spec_to_hedge_data(source, node_count)?;
            let flow = edge
                .flow
                .as_deref()
                .map(parse_flow)
                .transpose()?
                .unwrap_or(Flow::Source);
            builder.add_external_edge(source, edge_data, orientation, flow);
        }
        (None, Some(sink)) => {
            let sink = endpoint_spec_to_hedge_data(sink, node_count)?;
            let flow = edge
                .flow
                .as_deref()
                .map(parse_flow)
                .transpose()?
                .unwrap_or(Flow::Sink);
            builder.add_external_edge(sink, edge_data, orientation, flow);
        }
        (None, None) => return Err("Graph edge must have a source or sink endpoint".into()),
    }

    Ok(())
}

fn endpoint_spec_to_hedge_data(
    endpoint: TypstEndpointSpec,
    node_count: usize,
) -> Result<HedgeData<DotHedgeData>, String> {
    if endpoint.node >= node_count {
        return Err(format!(
            "Endpoint references node {}, but graph has {node_count} nodes",
            endpoint.node
        ));
    }

    Ok(HedgeData {
        data: DotHedgeData {
            statement: endpoint.statement,
            id: endpoint.id.map(Hedge),
            port_label: endpoint.port_label,
            compasspt: endpoint
                .compass
                .as_deref()
                .map(parse_endpoint_compass)
                .transpose()?
                .flatten(),
        },
        is_in_subgraph: endpoint.in_subgraph,
        node: NodeIndex(endpoint.node),
    })
}

fn parse_orientation(value: &str) -> Result<Orientation, String> {
    match value.trim().to_ascii_lowercase().as_str() {
        "default" | "forward" | "source" => Ok(Orientation::Default),
        "reversed" | "reverse" | "back" | "sink" => Ok(Orientation::Reversed),
        "undirected" | "none" => Ok(Orientation::Undirected),
        other => Err(format!("Invalid edge orientation: {other}")),
    }
}

fn parse_flow(value: &str) -> Result<Flow, String> {
    match value.trim().to_ascii_lowercase().as_str() {
        "source" | "out" | "outgoing" => Ok(Flow::Source),
        "sink" | "in" | "incoming" => Ok(Flow::Sink),
        other => Err(format!("Invalid external edge flow: {other}")),
    }
}

fn parse_endpoint_compass(value: &str) -> Result<Option<CompassPt>, String> {
    parse_compass(value)
}

fn statement_value<'a>(data: EdgeData<&'a DotEdgeData>, key: &str) -> Option<&'a str> {
    data.data.statements.get(key).map(String::as_str)
}

fn hedge_value(data: &DotHedgeData, key: &str) -> Option<String> {
    match key {
        "statement" => data.statement.clone(),
        "port_label" => data.port_label.clone(),
        "compass" => data.compasspt.map(compass_pt_to_string),
        "id" => data.id.map(|id| id.0.to_string()),
        _ => None,
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

fn compass_pt_to_string(compass: CompassPt) -> String {
    match compass {
        CompassPt::N => "n",
        CompassPt::NE => "ne",
        CompassPt::E => "e",
        CompassPt::SE => "se",
        CompassPt::S => "s",
        CompassPt::SW => "sw",
        CompassPt::W => "w",
        CompassPt::NW => "nw",
        CompassPt::C => "c",
        CompassPt::Underscore => "_",
    }
    .to_string()
}
