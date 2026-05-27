use std::collections::BTreeMap;

use dot_parser::ast::CompassPt;
use linnet::{
    half_edge::{
        builder::{HedgeData, HedgeGraphBuilder},
        involution::{
            ArchivedOrientation, EdgeData, EdgeIndex, Flow, Hedge, HedgePair, Orientation,
        },
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

use crate::{default_figment, TypstGraph};

type DotBuilder = HedgeGraphBuilder<DotEdgeData, DotVertexData, DotHedgeData>;
const TYPST_EDGE_NAME_KEY: &str = "__linnest-edge-name";

fn normalize_statement_key(key: &str) -> String {
    key.trim().trim_matches('"').to_string()
}

fn statement_map_value<'a>(
    statements: &'a BTreeMap<String, String>,
    key: &str,
) -> Option<&'a String> {
    statements
        .get(key)
        .or_else(|| statements.get(&format!("\"{key}\"")))
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct TypstDotGraphInfo {
    pub name: String,
    pub data: Option<Vec<u8>>,
    pub global_statements: BTreeMap<String, String>,
    pub default_edge_statements: BTreeMap<String, String>,
    pub default_node_statements: BTreeMap<String, String>,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct TypstPoint {
    pub x: f64,
    pub y: f64,
}

#[derive(
    Debug,
    Deserialize,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Default,
    rkyv::Archive,
    rkyv::Serialize,
    rkyv::Deserialize,
)]
#[serde(rename_all = "lowercase")]
enum PlacementMode {
    Start,
    #[default]
    Pin,
}

#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct TypstPlacementSpec {
    #[serde(default)]
    mode: PlacementMode,
    #[serde(default)]
    x_mode: Option<PlacementMode>,
    #[serde(default)]
    y_mode: Option<PlacementMode>,
    #[serde(default)]
    x: Option<TypstPlacementCoord>,
    #[serde(default)]
    y: Option<TypstPlacementCoord>,
    #[serde(default, rename = "ref")]
    reference: Option<usize>,
    #[serde(default, deserialize_with = "deserialize_optional_f64")]
    dx: Option<f64>,
    #[serde(default, deserialize_with = "deserialize_optional_f64")]
    dy: Option<f64>,
}

#[derive(Debug, Deserialize, Clone)]
#[serde(untagged)]
enum TypstPlacementCoord {
    Number(TypstNumber),
    Group(TypstPlacementGroup),
}

#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(untagged)]
enum TypstNumber {
    Float(f64),
    Signed(i64),
    Unsigned(u64),
}

#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
struct TypstPlacementGroup {
    kind: String,
    name: String,
    #[serde(default)]
    side: Option<String>,
}

#[derive(Debug, Clone, Copy, Default, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)]
pub struct ResolvedPoint {
    x: f64,
    y: f64,
    x_set: bool,
    y_set: bool,
}

#[derive(Debug, Clone)]
struct ResolvedPlacement {
    point: ResolvedPoint,
    pin: Option<String>,
    mode: PlacementMode,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct TypstDotNode {
    pub node: usize,
    pub name: Option<String>,
    pub data: Option<Vec<u8>>,
    pub pos: Option<TypstPoint>,
    pub shift: Option<TypstPoint>,
    pub statements: BTreeMap<String, String>,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct TypstDotEndpoint {
    pub node: usize,
    pub hedge: usize,
    pub data: Option<Vec<u8>>,
    pub statement: Option<String>,
    pub port_label: Option<String>,
    pub compass: Option<String>,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct TypstDotEdge {
    pub edge: usize,
    pub name: Option<String>,
    pub data: Option<Vec<u8>>,
    pub orientation: String,
    pub source: Option<TypstDotEndpoint>,
    pub sink: Option<TypstDotEndpoint>,
    pub pos: Option<TypstPoint>,
    pub shift: Option<TypstPoint>,
    pub label_pos: Option<TypstPoint>,
    pub label_angle: Option<f64>,
    pub bend: Option<f64>,
    pub statements: BTreeMap<String, String>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub struct TypstGraphSpec {
    #[serde(default)]
    pub name: Option<String>,
    #[serde(default)]
    pub data: Option<Vec<u8>>,
    #[serde(default)]
    pub statements: BTreeMap<String, String>,
    #[serde(default)]
    pub default_edge_statements: BTreeMap<String, String>,
    #[serde(default)]
    pub default_node_statements: BTreeMap<String, String>,
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
    pub data: Option<Vec<u8>>,
    #[serde(default)]
    pub pos: Option<TypstPlacementSpec>,
    #[serde(default)]
    pub statements: BTreeMap<String, String>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub struct TypstEdgeSpec {
    #[serde(default)]
    pub source: Option<TypstEndpointSpec>,
    #[serde(default)]
    pub sink: Option<TypstEndpointSpec>,
    #[serde(default)]
    pub data: Option<Vec<u8>>,
    #[serde(default)]
    pub orientation: Option<String>,
    #[serde(default)]
    pub flow: Option<String>,
    #[serde(default)]
    pub id: Option<usize>,
    #[serde(default)]
    pub pos: Option<TypstPlacementSpec>,
    #[serde(default)]
    pub statements: BTreeMap<String, String>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub struct TypstEndpointSpec {
    pub node: usize,
    #[serde(default)]
    pub statement: Option<String>,
    #[serde(default)]
    pub id: Option<usize>,
    #[serde(default)]
    pub data: Option<Vec<u8>>,
    #[serde(default)]
    pub port_label: Option<String>,
    #[serde(default)]
    pub compass: Option<String>,
    #[serde(default)]
    pub in_subgraph: bool,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub struct TypstGraphDataPatch {
    #[serde(default)]
    pub data: Option<Vec<u8>>,
    #[serde(default)]
    pub nodes: Vec<TypstIndexedDataPatch>,
    #[serde(default)]
    pub edges: Vec<TypstEdgeDataPatch>,
    #[serde(default)]
    pub hedges: Vec<TypstIndexedDataPatch>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct TypstGraphStructuralPatch {
    #[serde(default)]
    pub nodes: Vec<TypstNodeStructuralPatch>,
    #[serde(default)]
    pub edges: Vec<TypstEdgeStructuralPatch>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct TypstNodeStructuralPatch {
    pub index: usize,
    #[serde(default)]
    pub pos: Option<TypstPlacementSpec>,
    #[serde(default)]
    pub shift: Option<TypstPointSpec>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct TypstEdgeStructuralPatch {
    pub index: usize,
    #[serde(default)]
    pub pos: Option<TypstPlacementSpec>,
    #[serde(default)]
    pub shift: Option<TypstPointSpec>,
    #[serde(default)]
    pub label_pos: Option<TypstPointSpec>,
    #[serde(default, deserialize_with = "deserialize_optional_f64")]
    pub label_angle: Option<f64>,
    #[serde(default, deserialize_with = "deserialize_optional_f64")]
    pub bend: Option<f64>,
}

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum TypstPointSpec {
    Point { x: TypstNumber, y: TypstNumber },
    Tuple(TypstNumber, TypstNumber),
    Text(String),
}

#[derive(Debug, Deserialize)]
pub struct TypstIndexedDataPatch {
    pub index: usize,
    #[serde(default)]
    pub data: Option<Vec<u8>>,
}

#[derive(Debug, Deserialize)]
pub struct TypstEdgeDataPatch {
    pub index: usize,
    #[serde(default)]
    pub data: Option<Vec<u8>>,
    #[serde(default)]
    pub source: Option<Vec<u8>>,
    #[serde(default)]
    pub sink: Option<Vec<u8>>,
}

#[derive(Debug, Deserialize)]
pub struct TypstNamedDataPatch {
    pub name: String,
    pub data: Vec<u8>,
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

fn decode_subgraph(arg: &[u8]) -> Result<SuBitGraph, String> {
    unsafe { rkyv::from_bytes_unchecked::<SuBitGraph>(arg) }
        .map_err(|err| format!("Failed to deserialize archived subgraph: {err}"))
}

fn decode_typst_graph(arg: &[u8]) -> Result<TypstGraph, String> {
    unsafe { rkyv::from_bytes_unchecked::<TypstGraph>(arg) }
        .map_err(|err| format!("Failed to deserialize archived Typst graph: {err}"))
}

fn encode_typst_graph(graph: &TypstGraph) -> Result<Vec<u8>, String> {
    to_rkyv_bytes::<_, 4096>(graph)
}

fn typst_graph_from_dot(dot: DotGraph) -> TypstGraph {
    TypstGraph::from_dot(dot, &default_figment())
}

fn dot_bytes_from_typst_graph(graph: &TypstGraph) -> Result<Vec<u8>, String> {
    graph
        .to_dot_graph()
        .to_rkyv_bytes::<4096>()
        .map(|bytes| bytes.to_vec())
}

fn with_dot_view<T>(
    graph: &TypstGraph,
    callback: impl FnOnce(&ArchivedDotGraphView<'_>) -> Result<T, String>,
) -> Result<T, String> {
    let dot_bytes = dot_bytes_from_typst_graph(graph)?;
    let dot = DotGraph::archived_view(&dot_bytes);
    callback(&dot)
}

fn encode_subgraph(subgraph: &SuBitGraph) -> Result<Vec<u8>, String> {
    to_rkyv_bytes::<_, 256>(subgraph)
}

pub fn graph_from_spec_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: TypstGraphSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize graph spec: {err}"))?;
    encode_typst_graph(&typst_graph_from_dot(graph_from_spec(spec)?))
}

pub fn graph_with_data_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let layout_config = graph.layout_config.clone();
    let mut graph = graph.to_dot_graph();
    let patch: TypstGraphDataPatch = decode_cbor(arg2, "graph data patch")?;
    apply_graph_data_patch(&mut graph, patch)?;
    encode_typst_graph(&TypstGraph::from_dot_with_layout_config(
        graph,
        layout_config,
    ))
}

pub fn graph_node_data_by_name_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let name: String = decode_cbor(arg2, "node name")?;
    let graph = decode_typst_graph(arg)?;
    let data = with_dot_view(&graph, |graph| {
        Ok(graph
            .vertex_data()
            .find(|node| {
                node.data
                    .name
                    .as_ref()
                    .is_some_and(|value| value.as_str() == name)
            })
            .ok_or_else(|| format!("No node named {name:?}"))?
            .data
            .payload
            .as_ref()
            .map(|value| value.to_vec()))
    })?;
    encode_cbor(&data)
}

pub fn graph_edge_data_by_name_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let name: String = decode_cbor(arg2, "edge name")?;
    let graph = decode_typst_graph(arg)?;
    let data = with_dot_view(&graph, |graph| {
        Ok(graph
            .edge_data()
            .find(|edge| archived_edge_name(*edge).as_deref() == Some(name.as_str()))
            .ok_or_else(|| format!("No edge named {name:?}"))?
            .data
            .payload
            .as_ref()
            .map(|value| value.to_vec()))
    })?;
    encode_cbor(&data)
}

pub fn graph_set_node_data_by_name_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let layout_config = graph.layout_config.clone();
    let mut graph = graph.to_dot_graph();
    let patch: TypstNamedDataPatch = decode_cbor(arg2, "named node data patch")?;
    let node = graph
        .graph
        .iter_nodes()
        .find_map(|(index, _, data)| {
            (data.name.as_deref() == Some(patch.name.as_str())).then_some(index)
        })
        .ok_or_else(|| format!("No node named {:?}", patch.name))?;
    graph.graph[node].payload = Some(patch.data);
    encode_typst_graph(&TypstGraph::from_dot_with_layout_config(
        graph,
        layout_config,
    ))
}

pub fn graph_set_edge_data_by_name_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let layout_config = graph.layout_config.clone();
    let mut graph = graph.to_dot_graph();
    let patch: TypstNamedDataPatch = decode_cbor(arg2, "named edge data patch")?;
    let edge = graph
        .graph
        .iter_edges()
        .find_map(|(_, index, data)| {
            (data
                .data
                .statements
                .get(TYPST_EDGE_NAME_KEY)
                .map(String::as_str)
                == Some(patch.name.as_str()))
            .then_some(index)
        })
        .ok_or_else(|| format!("No edge named {:?}", patch.name))?;
    graph.graph[edge].payload = Some(patch.data);
    encode_typst_graph(&TypstGraph::from_dot_with_layout_config(
        graph,
        layout_config,
    ))
}

pub fn graph_apply_structural_patches_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let layout_config = graph.layout_config.clone();
    let mut dot = graph.to_dot_graph();
    let patch: TypstGraphStructuralPatch = decode_cbor(arg2, "graph structural patch")?;
    apply_graph_structural_patch(&mut dot, patch)?;
    encode_typst_graph(&TypstGraph::from_dot_with_layout_config(dot, layout_config))
}

fn apply_graph_data_patch(graph: &mut DotGraph, patch: TypstGraphDataPatch) -> Result<(), String> {
    if let Some(data) = patch.data {
        graph.global_data.payload = Some(data);
    }

    for node in patch.nodes {
        if node.index >= graph.graph.n_nodes() {
            return Err(format!(
                "Node data patch index {} is out of bounds for graph with {} nodes",
                node.index,
                graph.graph.n_nodes()
            ));
        }
        if let Some(data) = node.data {
            graph.graph[NodeIndex(node.index)].payload = Some(data);
        }
    }

    for edge in patch.edges {
        if edge.index >= graph.graph.n_edges() {
            return Err(format!(
                "Edge data patch index {} is out of bounds for graph with {} edges",
                edge.index,
                graph.graph.n_edges()
            ));
        }
        let edge_index = EdgeIndex(edge.index);
        let pair = graph
            .graph
            .iter_edges()
            .find_map(|(pair, index, _)| (index == edge_index).then_some(pair))
            .ok_or_else(|| format!("Edge data patch index {} could not be resolved", edge.index))?;

        if let Some(data) = edge.data {
            graph.graph[edge_index].payload = Some(data);
        }
        if let Some(data) = edge.source {
            let hedge = endpoint_hedge(pair, Flow::Source, edge.index)?;
            graph.graph[hedge].payload = Some(data);
        }
        if let Some(data) = edge.sink {
            let hedge = endpoint_hedge(pair, Flow::Sink, edge.index)?;
            graph.graph[hedge].payload = Some(data);
        }
    }

    for hedge in patch.hedges {
        if hedge.index >= graph.graph.n_hedges() {
            return Err(format!(
                "Half-edge data patch index {} is out of bounds for graph with {} half-edges",
                hedge.index,
                graph.graph.n_hedges()
            ));
        }
        if let Some(data) = hedge.data {
            graph.graph[Hedge(hedge.index)].payload = Some(data);
        }
    }

    Ok(())
}

fn apply_graph_structural_patch(
    graph: &mut DotGraph,
    patch: TypstGraphStructuralPatch,
) -> Result<(), String> {
    let mut node_positions = graph_node_positions(graph);

    for node in patch.nodes {
        if node.index >= graph.graph.n_nodes() {
            return Err(format!(
                "Node structural patch index {} is out of bounds for graph with {} nodes",
                node.index,
                graph.graph.n_nodes()
            ));
        }

        let index = NodeIndex(node.index);
        if let Some(pos) = node.pos {
            let placement = pos.resolve(&node_positions, "node structural patch")?;
            let statements = std::mem::take(&mut graph.graph[index].statements);
            graph.graph[index].statements =
                apply_placement_statements(statements, Some(&placement));
            node_positions[node.index] = placement.point;
        }
        if let Some(shift) = node.shift {
            graph.graph[index]
                .statements
                .insert("shift".to_string(), shift.to_statement()?);
        }
    }

    for edge in patch.edges {
        if edge.index >= graph.graph.n_edges() {
            return Err(format!(
                "Edge structural patch index {} is out of bounds for graph with {} edges",
                edge.index,
                graph.graph.n_edges()
            ));
        }

        let index = EdgeIndex(edge.index);
        let mut statements = std::mem::take(&mut graph.graph[index].statements);
        if let Some(pos) = edge.pos {
            let placement = pos.resolve(&node_positions, "edge structural patch")?;
            statements = apply_placement_statements(statements, Some(&placement));
        }
        if let Some(shift) = edge.shift {
            statements.insert("shift".to_string(), shift.to_statement()?);
        }
        if let Some(label_pos) = edge.label_pos {
            statements.insert("label-pos".to_string(), label_pos.to_statement()?);
        }
        if let Some(label_angle) = edge.label_angle {
            statements.insert("label-angle".to_string(), format!("{label_angle}rad"));
        }
        if let Some(bend) = edge.bend {
            statements.insert("bend".to_string(), format!("{bend}rad"));
        }
        graph.graph[index].local_statements = statements.clone();
        graph.graph[index].statements = statements;
    }

    Ok(())
}

fn graph_node_positions(graph: &DotGraph) -> Vec<ResolvedPoint> {
    let mut positions = vec![ResolvedPoint::default(); graph.graph.n_nodes()];
    for (index, _, data) in graph.graph.iter_nodes() {
        positions[index.0] = resolved_point_from_statements(&data.statements);
    }
    positions
}

fn endpoint_hedge(pair: HedgePair, side: Flow, edge_index: usize) -> Result<Hedge, String> {
    match (pair, side) {
        (HedgePair::Paired { source, .. } | HedgePair::Split { source, .. }, Flow::Source) => {
            Ok(source)
        }
        (HedgePair::Paired { sink, .. } | HedgePair::Split { sink, .. }, Flow::Sink) => Ok(sink),
        (HedgePair::Unpaired { hedge, flow }, side) if flow == side => Ok(hedge),
        (HedgePair::Unpaired { flow, .. }, side) => Err(format!(
            "Edge {edge_index} has only a {flow:?} endpoint, cannot patch {side:?} data"
        )),
    }
}

pub fn decode_graph_bytes_list(arg: &[u8]) -> Result<Vec<Vec<u8>>, String> {
    ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize graph bytes list: {err}"))
}

pub fn graph_info_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let info = with_dot_view(&graph, |graph| Ok(graph_info(graph)))?;
    encode_cbor(&info)
}

pub fn graph_dot_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    encode_cbor(&graph.to_dot_graph().debug_dot())
}

pub fn graph_nodes_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let nodes = with_dot_view(&graph, |graph| {
        Ok(graph
            .vertex_data()
            .map(node_view_to_output)
            .collect::<Vec<_>>())
    })?;
    encode_cbor(&nodes)
}

pub fn graph_nodes_of_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let nodes = with_dot_view(&graph, |graph| {
        let subgraph = decode_subgraph_spec(graph, arg2)?;
        Ok(graph
            .vertex_data_of(&subgraph)
            .map(node_view_to_output)
            .collect::<Vec<_>>())
    })?;
    encode_cbor(&nodes)
}

pub fn graph_edges_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let edges = with_dot_view(&graph, |graph| {
        Ok(graph
            .edge_data()
            .map(|edge| edge_view_to_output(graph, edge))
            .collect::<Vec<_>>())
    })?;
    encode_cbor(&edges)
}

pub fn graph_edges_of_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let edges = with_dot_view(&graph, |graph| {
        let subgraph = decode_subgraph_spec(graph, arg2)?;
        Ok(graph
            .edge_data_of(&subgraph)
            .map(|edge| edge_view_to_output(graph, edge))
            .collect::<Vec<_>>())
    })?;
    encode_cbor(&edges)
}

pub fn graph_nodes_of_archived_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let subgraph = decode_subgraph(arg2)?;
    let nodes = with_dot_view(&graph, |graph| {
        Ok(graph
            .vertex_data_of(&subgraph)
            .map(node_view_to_output)
            .collect::<Vec<_>>())
    })?;
    encode_cbor(&nodes)
}

pub fn graph_edges_of_archived_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let subgraph = decode_subgraph(arg2)?;
    let edges = with_dot_view(&graph, |graph| {
        Ok(graph
            .edge_data_of(&subgraph)
            .map(|edge| edge_view_to_output(graph, edge))
            .collect::<Vec<_>>())
    })?;
    encode_cbor(&edges)
}

pub fn graph_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let subgraph = with_dot_view(&graph, |graph| decode_subgraph_spec(graph, arg2))?;
    encode_cbor(&subgraph.string_label())
}

pub fn graph_archived_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let subgraph = with_dot_view(&graph, |graph| decode_subgraph_spec(graph, arg2))?;
    encode_subgraph(&subgraph)
}

pub fn graph_compass_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let compass = decode_compass(arg2)?;
    let label = graph
        .to_dot_graph()
        .compass_subgraph::<SuBitGraph>(compass)
        .string_label();
    encode_cbor(&label)
}

pub fn graph_archived_compass_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let compass = decode_compass(arg2)?;
    encode_subgraph(&graph.to_dot_graph().compass_subgraph::<SuBitGraph>(compass))
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
    let graph = decode_typst_graph(arg)?.to_dot_graph();
    let (cycles, _) = graph.cycle_basis();
    let archived = cycles
        .into_iter()
        .map(|cycle| encode_subgraph(&cycle.filter))
        .collect::<Result<Vec<_>, _>>()?;
    encode_cbor(&archived)
}

pub fn graph_spanning_forests_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?.to_dot_graph();
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
    let left_typst = decode_typst_graph(arg)?;
    let left_layout_config = left_typst.layout_config.clone();
    let left = left_typst.to_dot_graph();
    let right = decode_typst_graph(arg2)?.to_dot_graph();
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

    encode_typst_graph(&TypstGraph::from_dot_with_layout_config(
        DotGraph { global_data, graph },
        left_layout_config,
    ))
}

pub fn graph_join_by_hedge_key_bytes(
    arg: &[u8],
    arg2: &[u8],
    arg3: &[u8],
) -> Result<Vec<u8>, String> {
    let left_typst = decode_typst_graph(arg)?;
    let left_layout_config = left_typst.layout_config.clone();
    let left = left_typst.to_dot_graph();
    let right = decode_typst_graph(arg2)?.to_dot_graph();
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

    encode_typst_graph(&TypstGraph::from_dot_with_layout_config(
        DotGraph { global_data, graph },
        left_layout_config,
    ))
}

fn graph_info(graph: &ArchivedDotGraphView<'_>) -> TypstDotGraphInfo {
    let global_data = graph.global_data();
    TypstDotGraphInfo {
        name: global_data.name.as_str().to_string(),
        data: global_data.payload.as_ref().map(|value| value.to_vec()),
        global_statements: global_data
            .statements
            .iter()
            .map(|(key, value)| {
                (
                    normalize_statement_key(key.as_str()),
                    value.as_str().to_string(),
                )
            })
            .collect(),
        default_edge_statements: global_data
            .edge_statements
            .iter()
            .map(|(key, value)| {
                (
                    normalize_statement_key(key.as_str()),
                    value.as_str().to_string(),
                )
            })
            .collect(),
        default_node_statements: global_data
            .node_statements
            .iter()
            .map(|(key, value)| {
                (
                    normalize_statement_key(key.as_str()),
                    value.as_str().to_string(),
                )
            })
            .collect(),
    }
}

fn archived_edge_name(edge: ArchivedDotEdgeView<'_>) -> Option<String> {
    edge.data.statements.iter().find_map(|(key, value)| {
        (key.as_str() == TYPST_EDGE_NAME_KEY).then(|| value.as_str().to_string())
    })
}

fn graph_from_spec(spec: TypstGraphSpec) -> Result<DotGraph, String> {
    let node_count = spec.nodes.len();
    let mut builder = HedgeGraphBuilder::<DotEdgeData, DotVertexData, DotHedgeData>::new();
    let mut node_positions = Vec::with_capacity(node_count);
    let global_data = global_data_from_parts(
        spec.name,
        spec.data,
        spec.statements,
        spec.default_edge_statements,
        spec.default_node_statements,
    );

    for (node_index, node) in spec.nodes.into_iter().enumerate() {
        let placement = resolve_placement(node.pos.as_ref(), &node_positions, "node")?;
        let node_data = node_data_from_spec(&global_data, node, node_index, placement.as_ref());
        node_positions.push(
            placement
                .as_ref()
                .map(|placement| placement.point)
                .unwrap_or_else(|| resolved_point_from_statements(&node_data.statements)),
        );
        builder.add_node(node_data);
    }

    for edge in spec.edges {
        add_edge_to_builder(
            &global_data,
            &mut builder,
            node_count,
            &node_positions,
            edge,
        )?;
    }

    let mut graph = DotGraph {
        global_data,
        graph: builder.build::<DefaultNodeStore<DotVertexData>>(),
    };
    graph.apply_explicit_id_ordering();
    Ok(graph)
}

fn add_edge_to_builder(
    global_data: &GlobalData,
    builder: &mut DotBuilder,
    node_count: usize,
    node_positions: &[ResolvedPoint],
    edge: TypstEdgeSpec,
) -> Result<(), String> {
    let orientation = edge
        .orientation
        .as_deref()
        .map(parse_orientation)
        .transpose()?
        .unwrap_or(Orientation::Default);
    let placement = resolve_placement(edge.pos.as_ref(), node_positions, "edge")?;
    let local_statements = apply_placement_statements(edge.statements, placement.as_ref());
    let edge_data = DotEdgeData {
        payload: edge.data,
        statements: merged_statements(&global_data.edge_statements, &local_statements),
        local_statements,
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

fn global_data_from_parts(
    name: Option<String>,
    data: Option<Vec<u8>>,
    statements: BTreeMap<String, String>,
    edge_statements: BTreeMap<String, String>,
    node_statements: BTreeMap<String, String>,
) -> GlobalData {
    GlobalData {
        name: name.unwrap_or_else(|| "constructed".to_string()),
        payload: data,
        statements,
        edge_statements,
        node_statements,
    }
}

fn resolve_placement(
    placement: Option<&TypstPlacementSpec>,
    references: &[ResolvedPoint],
    context: &str,
) -> Result<Option<ResolvedPlacement>, String> {
    placement
        .map(|placement| placement.resolve(references, context))
        .transpose()
}

fn apply_placement_statements(
    mut statements: BTreeMap<String, String>,
    placement: Option<&ResolvedPlacement>,
) -> BTreeMap<String, String> {
    let Some(placement) = placement else {
        return statements;
    };

    if placement.point.x_set || placement.point.y_set {
        statements.insert(
            "pos".to_string(),
            format!("{},{}", placement.point.x, placement.point.y),
        );
        statements.insert("pos-x-set".to_string(), placement.point.x_set.to_string());
        statements.insert("pos-y-set".to_string(), placement.point.y_set.to_string());
        statements.insert(
            "pos-mode".to_string(),
            match placement.mode {
                PlacementMode::Start => "start",
                PlacementMode::Pin => "pin",
            }
            .to_string(),
        );
    }

    if let Some(pin) = &placement.pin {
        statements.insert("pin".to_string(), pin.clone());
    }

    statements
}

fn resolved_point_from_statements(statements: &BTreeMap<String, String>) -> ResolvedPoint {
    parse_statement_point(statements, "pos").unwrap_or_default()
}

fn parse_statement_point(
    statements: &BTreeMap<String, String>,
    attr: &str,
) -> Option<ResolvedPoint> {
    let value = statements.get(attr)?;
    let unquoted = value.trim().trim_matches('"');
    let cleaned = unquoted.trim().trim_matches(|c| c == '(' || c == ')');
    let parts: Vec<&str> = cleaned
        .split([',', ' '])
        .filter(|part| !part.is_empty())
        .collect();

    if parts.len() != 2 {
        return None;
    }

    let x = parts[0].trim().parse::<f64>().ok()?;
    let y = parts[1].trim().parse::<f64>().ok()?;
    Some(ResolvedPoint {
        x,
        y,
        x_set: parse_bool_statement(statements, "pos-x-set").unwrap_or(true),
        y_set: parse_bool_statement(statements, "pos-y-set").unwrap_or(true),
    })
}

fn parse_bool_statement(statements: &BTreeMap<String, String>, key: &str) -> Option<bool> {
    statements
        .get(key)
        .and_then(|value| value.trim().trim_matches('"').parse::<bool>().ok())
}

impl TypstPlacementSpec {
    fn resolve(
        &self,
        references: &[ResolvedPoint],
        context: &str,
    ) -> Result<ResolvedPlacement, String> {
        let mut point = if let Some(reference) = self.reference {
            let reference = references.get(reference).ok_or_else(|| {
                format!("{context} placement references node {reference}, but it is not available")
            })?;
            ResolvedPoint {
                x: reference.x + self.dx.unwrap_or(0.0),
                y: reference.y + self.dy.unwrap_or(0.0),
                x_set: true,
                y_set: true,
            }
        } else {
            ResolvedPoint {
                x: self.dx.unwrap_or(0.0),
                y: self.dy.unwrap_or(0.0),
                x_set: self.dx.is_some(),
                y_set: self.dy.is_some(),
            }
        };

        let x_mode = self.x_mode.unwrap_or(self.mode);
        let y_mode = self.y_mode.unwrap_or(self.mode);
        let x_pin = self
            .x
            .as_ref()
            .map(|coord| coord.resolve("x", x_mode, &mut point))
            .transpose()?
            .flatten();
        let y_pin = self
            .y
            .as_ref()
            .map(|coord| coord.resolve("y", y_mode, &mut point))
            .transpose()?
            .flatten();

        let mut parts = Vec::new();
        if let Some(x_pin) = x_pin {
            parts.push(x_pin);
        } else if x_mode == PlacementMode::Pin && point.x_set {
            parts.push(format!("x:{}", point.x));
        }
        if let Some(y_pin) = y_pin {
            parts.push(y_pin);
        } else if y_mode == PlacementMode::Pin && point.y_set {
            parts.push(format!("y:{}", point.y));
        }

        let pin = if parts.is_empty() {
            if self.mode == PlacementMode::Pin && self.x_mode.is_none() && self.y_mode.is_none() {
                return Err(format!(
                    "{context} pin placement must constrain x, y, or ref"
                ));
            }
            None
        } else {
            Some(parts.join(","))
        };

        let mode = if pin.is_some() {
            PlacementMode::Pin
        } else {
            PlacementMode::Start
        };

        Ok(ResolvedPlacement { point, pin, mode })
    }
}

impl TypstPlacementCoord {
    fn resolve(
        &self,
        axis: &str,
        _mode: PlacementMode,
        point: &mut ResolvedPoint,
    ) -> Result<Option<String>, String> {
        match self {
            TypstPlacementCoord::Number(value) => {
                let value = value.as_f64();
                if axis == "x" {
                    point.x = value;
                    point.x_set = true;
                } else {
                    point.y = value;
                    point.y_set = true;
                }
                Ok(None)
            }
            TypstPlacementCoord::Group(group) => {
                if group.kind != "group" {
                    return Err(format!(
                        "placement {axis} coordinate expected graph.group(...), got kind {:?}",
                        group.kind
                    ));
                }
                if axis == "x" {
                    point.x_set = true;
                } else {
                    point.y_set = true;
                }
                Ok(Some(format!("{axis}:{}", group.pin_token()?)))
            }
        }
    }
}

impl TypstNumber {
    fn as_f64(self) -> f64 {
        match self {
            TypstNumber::Float(value) => value,
            TypstNumber::Signed(value) => value as f64,
            TypstNumber::Unsigned(value) => value as f64,
        }
    }
}

impl TypstPointSpec {
    fn to_statement(&self) -> Result<String, String> {
        match self {
            TypstPointSpec::Point { x, y } => Ok(format!("{},{}", x.as_f64(), y.as_f64())),
            TypstPointSpec::Tuple(x, y) => Ok(format!("{},{}", x.as_f64(), y.as_f64())),
            TypstPointSpec::Text(value) => {
                let value = value.trim();
                if parse_point_text(value).is_some() {
                    Ok(value.to_string())
                } else {
                    Err(format!(
                        "point value {value:?} must contain two numeric coordinates"
                    ))
                }
            }
        }
    }
}

fn parse_point_text(value: &str) -> Option<(f64, f64)> {
    let cleaned = value
        .trim()
        .trim_matches('"')
        .trim_matches(|c| c == '(' || c == ')');
    let parts: Vec<&str> = cleaned
        .split([',', ' '])
        .filter(|part| !part.is_empty())
        .collect();
    if parts.len() != 2 {
        return None;
    }
    Some((parts[0].parse().ok()?, parts[1].parse().ok()?))
}

impl TypstPlacementGroup {
    fn pin_token(&self) -> Result<String, String> {
        match self.side.as_deref() {
            None => Ok(format!("@{}", self.name)),
            Some("+") | Some("positive") => Ok(format!("@+{}", self.name)),
            Some("-") | Some("negative") => Ok(format!("@-{}", self.name)),
            Some(side) => Err(format!(
                "graph.group side must be none, \"+\", \"-\", \"positive\", or \"negative\", got {side:?}"
            )),
        }
    }
}

fn deserialize_optional_f64<'de, D>(deserializer: D) -> Result<Option<f64>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    Option::<TypstNumber>::deserialize(deserializer).map(|value| value.map(TypstNumber::as_f64))
}

fn node_data_from_spec(
    global_data: &GlobalData,
    node: TypstNodeSpec,
    default_index: usize,
    placement: Option<&ResolvedPlacement>,
) -> DotVertexData {
    DotVertexData {
        name: node.name,
        index: node.index.map(NodeIndex).or(Some(NodeIndex(default_index))),
        payload: node.data,
        statements: merged_statements(
            &global_data.node_statements,
            &apply_placement_statements(node.statements, placement),
        ),
    }
}

fn merged_statements(
    defaults: &BTreeMap<String, String>,
    local: &BTreeMap<String, String>,
) -> BTreeMap<String, String> {
    let mut statements = defaults.clone();
    statements.extend(local.clone());
    statements
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
            payload: endpoint.data,
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
        "port-label" => data.port_label.clone(),
        "compass" => data.compasspt.map(compass_pt_to_string),
        "id" => data.id.map(|id| id.0.to_string()),
        _ => None,
    }
}

fn node_view_to_output(vertex: ArchivedDotVertexView<'_>) -> TypstDotNode {
    let raw_statements = vertex
        .data
        .statements
        .iter()
        .map(|(key, value)| {
            (
                normalize_statement_key(key.as_str()),
                value.as_str().to_string(),
            )
        })
        .collect::<BTreeMap<_, _>>();

    TypstDotNode {
        node: vertex.node.0,
        name: vertex
            .data
            .name
            .as_ref()
            .map(|value| value.as_str().to_string()),
        data: vertex.data.payload.as_ref().map(|value| value.to_vec()),
        pos: parse_point(&raw_statements, "pos"),
        shift: parse_point(&raw_statements, "shift"),
        statements: public_statements(raw_statements),
    }
}

fn edge_view_to_output(
    graph: &ArchivedDotGraphView<'_>,
    edge: ArchivedDotEdgeView<'_>,
) -> TypstDotEdge {
    let raw_statements = edge
        .data
        .statements
        .iter()
        .map(|(key, value)| {
            (
                normalize_statement_key(key.as_str()),
                value.as_str().to_string(),
            )
        })
        .collect::<BTreeMap<_, _>>();
    let endpoints = graph.endpoints_of_edge(&edge);

    TypstDotEdge {
        edge: edge.edge.0,
        name: archived_edge_name(edge),
        data: edge.data.payload.as_ref().map(|value| value.to_vec()),
        orientation: orientation_to_string(edge.orientation).to_string(),
        source: endpoints.source.map(endpoint_to_output),
        sink: endpoints.sink.map(endpoint_to_output),
        pos: parse_point(&raw_statements, "pos"),
        shift: parse_point(&raw_statements, "shift"),
        label_pos: parse_point(&raw_statements, "label-pos"),
        label_angle: parse_rad(&raw_statements, "label-angle"),
        bend: parse_rad(&raw_statements, "bend"),
        statements: public_statements(raw_statements),
    }
}

fn public_statements(mut statements: BTreeMap<String, String>) -> BTreeMap<String, String> {
    for key in [
        "pos",
        "pos-x-set",
        "pos-y-set",
        "pos-mode",
        "pin",
        TYPST_EDGE_NAME_KEY,
    ] {
        statements.remove(key);
    }
    statements
}

fn endpoint_to_output(endpoint: ArchivedDotEndpointView<'_>) -> TypstDotEndpoint {
    TypstDotEndpoint {
        node: endpoint.node.0,
        hedge: endpoint.hedge.0,
        data: endpoint.data.payload.as_ref().map(|value| value.to_vec()),
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
    let value = statement_map_value(statements, key)?;
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
    statement_map_value(statements, key).and_then(|value| {
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
