use std::collections::{BTreeMap, HashMap, HashSet};

use cgmath::{Point2, Rad, Vector2};
use dot_parser::ast::CompassPt;
use linnet::{
    half_edge::{
        builder::{HedgeData, HedgeGraphBuilder},
        involution::{
            ArchivedOrientation, EdgeData, EdgeIndex, Flow, Hedge, HedgePair, Involution,
            InvolutiveMapping, Orientation,
        },
        layout::spring::{Constraint, LayoutPointIndex, PointConstraint},
        nodestore::{DefaultNodeStore, NodeStorageOps},
        subgraph::{Inclusion, SuBitGraph, SubSetLike},
        swap::Swap,
        NodeIndex,
    },
    parser::{
        ArchivedDotEdgeView, ArchivedDotEndpointView, ArchivedDotGraphView, ArchivedDotVertexView,
        DotEdgeData, DotGraph, DotHedgeData, DotVertexData, GlobalData,
    },
};
use serde::{Deserialize, Serialize};

use crate::{default_figment, PinConstraint, TypstEdge, TypstGraph, TypstHedge, TypstNode};

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
    Serialize,
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

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct TypstPlacementSpec {
    #[serde(default)]
    mode: PlacementMode,
    #[serde(default)]
    x_mode: Option<PlacementMode>,
    #[serde(default)]
    y_mode: Option<PlacementMode>,
    #[serde(default)]
    z_mode: Option<PlacementMode>,
    #[serde(default)]
    x: Option<TypstPlacementCoord>,
    #[serde(default)]
    y: Option<TypstPlacementCoord>,
    #[serde(default)]
    z: Option<TypstNumber>,
    #[serde(default, rename = "ref")]
    reference: Option<usize>,
    #[serde(default, deserialize_with = "deserialize_optional_f64")]
    dx: Option<f64>,
    #[serde(default, deserialize_with = "deserialize_optional_f64")]
    dy: Option<f64>,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[serde(untagged)]
enum TypstPlacementCoord {
    Number(TypstNumber),
    Group(TypstPlacementGroup),
}

#[derive(Debug, Serialize, Deserialize, Clone, Copy, PartialEq)]
#[serde(untagged)]
enum TypstNumber {
    Float(f64),
    Signed(i64),
    Unsigned(u64),
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[serde(rename_all = "kebab-case")]
struct TypstPlacementGroup {
    kind: String,
    name: String,
    #[serde(default)]
    side: Option<String>,
    #[serde(default)]
    start: Option<TypstNumber>,
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
    point: Option<ResolvedPoint>,
    z: Option<(f64, PlacementMode)>,
    pin: Option<String>,
    mode: PlacementMode,
    group_start_x: bool,
    group_start_y: bool,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct TypstDotNode {
    pub node: usize,
    pub name: Option<String>,
    pub data: Option<Vec<u8>>,
    pub pos: Option<TypstPoint>,
    #[serde(rename = "pos-x-set")]
    pub pos_x_set: bool,
    #[serde(rename = "pos-y-set")]
    pub pos_y_set: bool,
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
    pub route_points: Vec<TypstPoint>,
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
    pub pos_x_set: bool,
    pub pos_y_set: bool,
    pub shift: Option<TypstPoint>,
    pub label_pos: Option<TypstPoint>,
    pub label_angle: Option<f64>,
    pub bend: Option<f64>,
    pub statements: BTreeMap<String, String>,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
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

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
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

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct TypstEdgeSpec {
    #[serde(default)]
    pub name: Option<String>,
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

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
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

pub const GRAPH_SPEC_SCHEMA: &str = "linnest-graph-spec";
pub const GRAPH_SPEC_VERSION: u32 = 1;

/// Versioned wire envelope shared by native renderers and the Typst plugin.
#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct TypstGraphSpecEnvelope {
    pub schema: String,
    pub version: u32,
    pub graph: TypstGraphSpec,
}

#[derive(Serialize)]
struct TypstGraphSpecRefEnvelope<'a> {
    schema: &'static str,
    version: u32,
    graph: &'a TypstGraphSpec,
}

impl TypstGraphSpecEnvelope {
    pub fn new(graph: TypstGraphSpec) -> Self {
        Self {
            schema: GRAPH_SPEC_SCHEMA.to_owned(),
            version: GRAPH_SPEC_VERSION,
            graph,
        }
    }

    fn into_graph(self) -> Result<TypstGraphSpec, String> {
        if self.schema != GRAPH_SPEC_SCHEMA {
            return Err(format!(
                "Invalid graph spec schema {:?}; expected {GRAPH_SPEC_SCHEMA:?}",
                self.schema
            ));
        }
        if self.version != GRAPH_SPEC_VERSION {
            return Err(format!(
                "Unsupported graph spec version {}; expected {GRAPH_SPEC_VERSION}",
                self.version
            ));
        }
        Ok(self.graph)
    }
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
    #[serde(default)]
    pub hedges: Vec<TypstHedgeStructuralPatch>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct TypstNodeStructuralPatch {
    pub index: usize,
    #[serde(default)]
    pub pos: Option<TypstPlacementSpec>,
    #[serde(default)]
    pub shift: Option<TypstPointSpec>,
    #[serde(default)]
    pub statements: BTreeMap<String, String>,
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
    #[serde(default)]
    pub statements: BTreeMap<String, String>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct TypstHedgeStructuralPatch {
    pub index: usize,
    #[serde(default)]
    pub statement: Option<String>,
    #[serde(default)]
    pub port_label: Option<String>,
    #[serde(default)]
    pub compass: Option<String>,
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
    let bytes = aligned_copy(arg);
    rkyv::from_bytes::<SuBitGraph>(&bytes)
        .map_err(|err| format!("Failed to deserialize archived subgraph: {err}"))
}

pub(crate) fn decode_typst_graph(arg: &[u8]) -> Result<TypstGraph, String> {
    let bytes = aligned_copy(arg);
    rkyv::from_bytes::<TypstGraph>(&bytes)
        .map_err(|err| format!("Failed to deserialize archived Typst graph: {err}"))
}

fn aligned_copy(bytes: &[u8]) -> rkyv::AlignedVec {
    let mut aligned = rkyv::AlignedVec::with_capacity(bytes.len());
    aligned.extend_from_slice(bytes);
    aligned
}

fn encode_typst_graph(graph: &TypstGraph) -> Result<Vec<u8>, String> {
    to_rkyv_bytes::<_, 4096>(graph)
}

fn typst_graph_from_dot(dot: DotGraph) -> TypstGraph {
    TypstGraph::from_dot(dot, &default_figment())
}

fn dot_bytes_from_typst_graph(graph: &TypstGraph) -> Result<rkyv::AlignedVec, String> {
    graph.to_dot_graph().to_rkyv_bytes::<4096>()
}

fn with_dot_view<T>(
    graph: &TypstGraph,
    callback: impl FnOnce(&ArchivedDotGraphView<'_>) -> Result<T, String>,
) -> Result<T, String> {
    let dot_bytes = dot_bytes_from_typst_graph(graph)?;
    let dot = DotGraph::archived_view(&dot_bytes)?;
    callback(&dot)
}

fn encode_subgraph(subgraph: &SuBitGraph) -> Result<Vec<u8>, String> {
    to_rkyv_bytes::<_, 256>(subgraph)
}

pub fn graph_from_spec_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let envelope: TypstGraphSpecEnvelope = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize graph spec: {err}"))?;
    encode_typst_graph(&typst_graph_from_dot(graph_from_spec(
        envelope.into_graph()?,
    )?))
}

pub fn encode_graph_spec_bytes(spec: &TypstGraphSpec) -> Result<Vec<u8>, String> {
    encode_cbor(&TypstGraphSpecRefEnvelope {
        schema: GRAPH_SPEC_SCHEMA,
        version: GRAPH_SPEC_VERSION,
        graph: spec,
    })
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
    let mut graph = decode_typst_graph(arg)?;
    let patch: TypstGraphStructuralPatch = decode_cbor(arg2, "graph structural patch")?;
    apply_typst_graph_structural_patch(&mut graph, patch)?;
    encode_typst_graph(&graph)
}

fn apply_typst_graph_structural_patch(
    graph: &mut TypstGraph,
    patch: TypstGraphStructuralPatch,
) -> Result<(), String> {
    let mut node_positions = typst_node_positions(graph);
    let mut refresh_positions = false;

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
            if let Some(point) = placement.point {
                refresh_positions = true;
                node_positions[node.index] = point;
            }
            let statements = std::mem::take(&mut graph.graph[index].statements);
            graph.graph[index].statements =
                apply_placement_statements(statements, Some(&placement));
        }
        if let Some(shift) = node.shift {
            graph.graph[index]
                .statements
                .insert("shift".to_string(), shift.to_statement()?);
        }
        refresh_positions |= point_statements_changed(&node.statements);
        graph.graph[index].statements.extend(node.statements);
        if let Some((x, y)) = statement_point(&graph.graph[index].statements, "shift") {
            graph.graph[index].shift = Some(Vector2::new(x, y));
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
        if let Some(pos) = edge.pos {
            let placement = pos.resolve(&node_positions, "edge structural patch")?;
            refresh_positions |= placement.point.is_some();
            let statements = std::mem::take(&mut graph.graph[index].statements);
            graph.graph[index].statements =
                apply_placement_statements(statements, Some(&placement));
        }
        if let Some(shift) = edge.shift {
            graph.graph[index]
                .statements
                .insert("shift".to_string(), shift.to_statement()?);
        }
        if let Some(label_pos) = edge.label_pos {
            graph.graph[index]
                .statements
                .insert("label-pos".to_string(), label_pos.to_statement()?);
        }
        if let Some(label_angle) = edge.label_angle {
            graph.graph[index]
                .statements
                .insert("label-angle".to_string(), format!("{label_angle}rad"));
        }
        if let Some(bend) = edge.bend {
            graph.graph[index]
                .statements
                .insert("bend".to_string(), format!("{bend}rad"));
        }
        refresh_positions |= point_statements_changed(&edge.statements);
        graph.graph[index].statements.extend(edge.statements);
        let data = &mut graph.graph[index];
        if let Some((x, y)) = statement_point(&data.statements, "shift") {
            data.shift = Some(Vector2::new(x, y));
        }
        if let Some((x, y)) = statement_point(&data.statements, "label-pos") {
            data.label_pos = Some(Point2::new(x, y));
        }
        if let Some(angle) = statement_radians(&data.statements, "label-angle") {
            data.label_angle = Some(angle);
        }
        if let Some(bend) = statement_radians(&data.statements, "bend") {
            data.bend = Ok(Rad(bend));
            data.bend_explicit = true;
        }
    }

    for hedge in patch.hedges {
        if hedge.index >= graph.graph.n_hedges() {
            return Err(format!(
                "Half-edge structural patch index {} is out of bounds for graph with {} half-edges",
                hedge.index,
                graph.graph.n_hedges()
            ));
        }

        let data = &mut graph.graph[Hedge(hedge.index)];
        if let Some(statement) = hedge.statement {
            data.statement = Some(statement);
        }
        if let Some(port_label) = hedge.port_label {
            data.port_label = Some(port_label);
        }
        if let Some(compass) = hedge.compass {
            data.compasspt = parse_endpoint_compass(&compass)?.map(compass_pt_to_string);
        }
    }

    if refresh_positions {
        refresh_structural_state_from_statements(graph);
    }
    Ok(())
}

fn point_statements_changed(statements: &BTreeMap<String, String>) -> bool {
    statements.keys().any(|key| {
        matches!(
            normalize_statement_key(key).as_str(),
            "pos"
                | "pin"
                | "pos-x-set"
                | "pos-y-set"
                | "pos-mode"
                | "group-start-x"
                | "group-start-y"
        )
    })
}

fn typst_node_positions(graph: &TypstGraph) -> Vec<ResolvedPoint> {
    let mut positions = vec![ResolvedPoint::default(); graph.graph.n_nodes()];
    for (index, _, node) in graph.graph.iter_nodes() {
        positions[index.0] = ResolvedPoint {
            x: node.pos.x,
            y: node.pos.y,
            x_set: true,
            y_set: true,
        };
    }
    positions
}

fn refresh_structural_state_from_statements(graph: &mut TypstGraph) {
    let mut group_map = HashMap::new();
    for index in 0..graph.graph.n_nodes() {
        let index = NodeIndex(index);
        let node = &mut graph.graph[index];
        let statements = node.statements.clone();
        if let Some((x, y)) = statement_point(&statements, "shift") {
            node.shift = Some(Vector2::new(x, y));
        }
        refresh_point_state(
            LayoutPointIndex::Node(index),
            &statements,
            &mut node.pos,
            &mut node.constraints,
            &mut node.start_x,
            &mut node.start_y,
            &mut group_map,
        );
    }

    for index in 0..graph.graph.n_edges() {
        let index = EdgeIndex(index);
        let edge = &mut graph.graph[index];
        let statements = edge.statements.clone();
        if let Some((x, y)) = statement_point(&statements, "shift") {
            edge.shift = Some(Vector2::new(x, y));
        }
        if let Some((x, y)) = statement_point(&statements, "label-pos") {
            edge.label_pos = Some(Point2::new(x, y));
        }
        if let Some(angle) = statement_radians(&statements, "label-angle") {
            edge.label_angle = Some(angle);
        }
        if let Some(bend) = statement_radians(&statements, "bend") {
            edge.bend = Ok(Rad(bend));
            edge.bend_explicit = true;
        }
        refresh_point_state(
            LayoutPointIndex::Edge(index),
            &statements,
            &mut edge.pos,
            &mut edge.constraints,
            &mut edge.start_x,
            &mut edge.start_y,
            &mut group_map,
        );
    }
}

fn refresh_point_state(
    index: LayoutPointIndex,
    statements: &BTreeMap<String, String>,
    point: &mut Point2<f64>,
    constraints: &mut PointConstraint,
    start_x: &mut bool,
    start_y: &mut bool,
    group_map: &mut HashMap<String, LayoutPointIndex>,
) {
    let position = parse_statement_point(statements, "pos");
    let pin = statement_map_value(statements, "pin").and_then(|value| PinConstraint::parse(value));

    match pin {
        Some(pin) => {
            let (mut pin_point, pin_constraints) = pin.point_constraint(index, group_map);
            if let Some(position) = position {
                if !matches!(pin_constraints.x, Constraint::Fixed) && position.x_set {
                    pin_point.x = position.x;
                }
                if !matches!(pin_constraints.y, Constraint::Fixed) && position.y_set {
                    pin_point.y = position.y;
                }
                *start_x = position.x_set;
                *start_y = position.y_set;
            } else {
                *start_x = false;
                *start_y = false;
            }
            *point = pin_point;
            *constraints = pin_constraints;
        }
        None => {
            if let Some(position) = position {
                *point = Point2::new(position.x, position.y);
                // Parsed group placements keep their constraints in the graph
                // while exposing only the resolved point to structural maps.
                if statement_map_value(statements, "pos-mode")
                    .is_none_or(|mode| mode.trim().trim_matches('"') != "pin")
                {
                    *constraints = PointConstraint::default();
                }
                *start_x = position.x_set;
                *start_y = position.y_set;
            }
        }
    }
}

fn statement_point(statements: &BTreeMap<String, String>, key: &str) -> Option<(f64, f64)> {
    statement_map_value(statements, key).and_then(|value| parse_point_text(value))
}

fn statement_radians(statements: &BTreeMap<String, String>, key: &str) -> Option<f64> {
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
    encode_cbor(&edges_with_route_points(&graph, edges))
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
    encode_cbor(&edges_with_route_points(&graph, edges))
}

pub fn graph_nodes_of_archived_subgraph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph = decode_typst_graph(arg)?;
    let subgraph = decode_subgraph(arg2)?;
    if subgraph.size() != graph.n_hedges() {
        return Err(format!(
            "Archived subgraph has {} bits, but graph has {} half-edges; sizes must match",
            subgraph.size(),
            graph.n_hedges()
        ));
    }
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
    if subgraph.size() != graph.n_hedges() {
        return Err(format!(
            "Archived subgraph has {} bits, but graph has {} half-edges; sizes must match",
            subgraph.size(),
            graph.n_hedges()
        ));
    }
    let edges = with_dot_view(&graph, |graph| {
        Ok(graph
            .edge_data_of(&subgraph)
            .map(|edge| edge_view_to_output(graph, edge))
            .collect::<Vec<_>>())
    })?;
    encode_cbor(&edges_with_route_points(&graph, edges))
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

#[cfg(any(test, all(target_arch = "wasm32", feature = "typst-plugin")))]
pub(crate) fn subgraph_size_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let subgraph = decode_subgraph(arg)?;
    encode_cbor(&subgraph.size())
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

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct TypstCutEntry {
    left: usize,
    right: usize,
    winding: usize,
}

#[derive(Debug, Serialize, Deserialize, PartialEq, Eq)]
struct TypstCutEdgeOrigin {
    edge: usize,
    segment: Option<usize>,
    winding: usize,
}

#[derive(Debug, Serialize, Deserialize, PartialEq, Eq)]
struct TypstCutBoundary {
    edge: usize,
    node: Option<usize>,
    hedge: usize,
    side: String,
    crossing: usize,
}

#[derive(Debug, Serialize, Deserialize)]
struct TypstCutResult {
    graph: Vec<u8>,
    nodes: Vec<Option<usize>>,
    edges: Vec<TypstCutEdgeOrigin>,
    hedges: Vec<Option<usize>>,
    boundaries: Vec<TypstCutBoundary>,
}

impl TypstEdge {
    fn cut_payload(&self) -> Option<ciborium::Value> {
        let mut bytes = self.data.as_deref()?;
        let value = ciborium::de::from_reader(&mut bytes).ok()?;
        bytes.is_empty().then_some(value)
    }

    fn cut_name(&self) -> Result<Option<String>, String> {
        let statement = statement_map_value(&self.statements, TYPST_EDGE_NAME_KEY);
        let mut name = statement.cloned();
        if let Some(ciborium::Value::Map(fields)) = self.cut_payload() {
            for (key, value) in fields {
                if key.as_text() != Some("name") || value == ciborium::Value::Null {
                    continue;
                }
                let value = value
                    .as_text()
                    .ok_or("Cut edge payload name must be text")?;
                if name.as_deref().is_some_and(|name| name != value) {
                    return Err("Conflicting cut edge names in statements and payload".into());
                }
                name = Some(value.to_owned());
            }
        }
        Ok(name)
    }

    fn cut_fragment(&self, name: String) -> Result<Self, String> {
        let mut statements = self.statements.clone();
        statements.retain(|key, _| {
            let key = normalize_statement_key(key);
            key != "pos"
                && !key.starts_with("pos-")
                && !matches!(
                    key.as_str(),
                    "pin"
                        | "z"
                        | "z-mode"
                        | "shift"
                        | "bend"
                        | "label-pos"
                        | "label-angle"
                        | "route-points"
                        | "group-start-x"
                        | "group-start-y"
                        | TYPST_EDGE_NAME_KEY
                )
        });
        statements.insert(TYPST_EDGE_NAME_KEY.into(), name.clone());
        // Payloads are otherwise opaque. Only rewrite an existing CBOR name;
        // sidecar data keys and every other field retain their input meanings.
        let mut data = self.data.clone();
        if let Some(ciborium::Value::Map(mut fields)) = self.cut_payload() {
            let mut renamed = false;
            for (key, value) in &mut fields {
                if key.as_text() == Some("name") {
                    *value = ciborium::Value::Text(name.clone());
                    renamed = true;
                }
            }
            if renamed {
                data = Some(encode_cbor(&ciborium::Value::Map(fields))?);
            }
        }
        Ok(Self {
            data,
            statements,
            ..Self::default()
        })
    }
}

impl TypstCutEntry {
    fn endpoints(&self, graph: &TypstGraph) -> Result<(Hedge, Hedge, EdgeIndex), String> {
        if self.winding == 0 {
            return Err("Cut winding must be positive".into());
        }
        if self.left >= graph.n_hedges() || self.right >= graph.n_hedges() {
            return Err(format!(
                "Cut hedge pair ({}, {}) is out of bounds",
                self.left, self.right
            ));
        }
        let (left, right) = (Hedge(self.left), Hedge(self.right));
        if left == right || graph.inv(left) != right || graph.inv(right) != left {
            return Err(format!(
                "Cut hedges ({left}, {right}) must be paired inverses"
            ));
        }
        let (source, sink) = match graph.underlying_hedge_orientation(left) {
            Flow::Source => (left, right),
            Flow::Sink => (right, left),
        };
        Ok((source, sink, graph[&source]))
    }

    fn boundary(
        &self,
        edge: usize,
        node: Option<usize>,
        hedge: Hedge,
        crossing: usize,
    ) -> TypstCutBoundary {
        TypstCutBoundary {
            edge,
            node,
            hedge: hedge.0,
            side: if hedge.0 == self.left {
                "left"
            } else {
                "right"
            }
            .into(),
            crossing,
        }
    }
}

impl TypstGraph {
    fn validate_cut_topology(&self) -> Result<(), String> {
        // rkyv checks the archive's memory representation, not cross-index
        // invariants. Validate these before calling Linnet's indexing APIs.
        let invalid = || "Invalid cut input graph topology".to_string();
        let involution: &Involution = self.graph.as_ref();
        let mappings: Vec<_> = involution.iter().map(|(_, mapping)| mapping).collect();
        let node_hedges: Hedge = self.graph.node_store.len();
        if mappings.len() != self.n_hedges()
            || self.iter_hedges().count() != self.n_hedges()
            || node_hedges.0 != self.n_hedges()
            || self.graph.node_store.iter().count() != self.n_nodes()
        {
            return Err(invalid());
        }
        let mut edges = vec![false; self.n_edges()];
        for (index, mapping) in mappings.iter().enumerate() {
            let hedge = Hedge(index);
            let data = match mapping {
                InvolutiveMapping::Identity { data, .. } => data,
                InvolutiveMapping::Source { data, sink_idx } => {
                    if !matches!(mappings.get(sink_idx.0), Some(InvolutiveMapping::Sink { source_idx }) if *source_idx == hedge)
                    {
                        return Err(invalid());
                    }
                    data
                }
                InvolutiveMapping::Sink { source_idx } => {
                    if !matches!(mappings.get(source_idx.0), Some(InvolutiveMapping::Source { sink_idx, .. }) if *sink_idx == hedge)
                    {
                        return Err(invalid());
                    }
                    continue;
                }
            };
            let seen = edges.get_mut(data.data.0).ok_or_else(invalid)?;
            if std::mem::replace(seen, true) {
                return Err(invalid());
            }
        }
        if edges.contains(&false) {
            return Err(invalid());
        }
        self.check().map_err(|err| err.to_string())?;
        let mut hedges = vec![false; self.n_hedges()];
        for (node, neighbors, _) in self.iter_nodes() {
            for hedge in neighbors {
                let seen = hedges.get_mut(hedge.0).ok_or_else(invalid)?;
                if std::mem::replace(seen, true) || self.node_id(hedge) != node {
                    return Err(invalid());
                }
            }
        }
        if hedges.contains(&false) {
            return Err(invalid());
        }
        Ok(())
    }

    /// Cut each selected paired edge into source/sink stubs and `winding - 1`
    /// disconnected middle edges. All origins refer to the immediate input graph,
    /// indexed by output IDs. Segment/crossing order follows underlying flow, not
    /// superficial orientation or the caller's choice of left/right.
    /// Requests are capped at 64 MiB of graph bytes and 1 MiB of cut CBOR;
    /// nonempty cuts also bound each output ID space to 65,536 entries and
    /// estimated expanded storage (including node bitsets) to 64 MiB.
    pub fn cut_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
        const MAX_BYTES: usize = 64 * 1024 * 1024;
        const MAX_ITEMS: usize = 65_536;
        let too_large = || "Cut exceeds topology or allocation limits".to_string();
        if arg.len() > MAX_BYTES || arg2.len() > 1024 * 1024 {
            return Err(too_large());
        }
        let mut remaining = arg2;
        let entries: Vec<TypstCutEntry> = ciborium::de::from_reader(&mut remaining)
            .map_err(|err| format!("Failed to deserialize cut entries: {err}"))?;
        if !remaining.is_empty() {
            return Err("Trailing bytes after cut entries CBOR array".into());
        }
        let mut graph = decode_typst_graph(arg)?;
        graph.validate_cut_topology()?;
        let (n_nodes, n_edges, n_hedges) = (graph.n_nodes(), graph.n_edges(), graph.n_hedges());
        let mut result = TypstCutResult {
            graph: Vec::new(),
            nodes: (0..n_nodes).map(Some).collect(),
            edges: (0..n_edges)
                .map(|edge| TypstCutEdgeOrigin {
                    edge,
                    segment: None,
                    winding: 0,
                })
                .collect(),
            hedges: (0..n_hedges).map(Some).collect(),
            boundaries: Vec::new(),
        };
        if entries.is_empty() {
            result.graph = arg.to_vec();
            return encode_cbor(&result);
        }
        if entries.len() > n_edges {
            return Err("Cut selects more entries than input edges".into());
        }

        let mut selected = HashSet::new();
        let mut windings = 0usize;
        let mut estimated_bytes = arg.len();
        // Preflight before splitting or allocating in proportion to a winding.
        for entry in &entries {
            let (source, sink, edge) = entry.endpoints(&graph)?;
            if !selected.insert(edge) {
                return Err(format!("Duplicate cut selection for edge {edge}"));
            }
            windings = windings
                .checked_add(entry.winding)
                .filter(|n| *n <= MAX_ITEMS)
                .ok_or_else(too_large)?;
            let fragment_bytes = to_rkyv_bytes::<_, 256>(&graph[edge])?.len()
                + to_rkyv_bytes::<_, 256>(&graph[source])?.len()
                + to_rkyv_bytes::<_, 256>(&graph[sink])?.len()
                + 1024;
            estimated_bytes = entry
                .winding
                .checked_mul(fragment_bytes)
                .and_then(|bytes| estimated_bytes.checked_add(bytes))
                .filter(|bytes| *bytes <= MAX_BYTES)
                .ok_or_else(too_large)?;
        }
        let extra_endpoints = 2 * (windings - entries.len());
        let output_nodes = n_nodes
            .checked_add(extra_endpoints)
            .filter(|n| *n <= MAX_ITEMS)
            .ok_or_else(too_large)?;
        let output_hedges = n_hedges
            .checked_add(extra_endpoints)
            .filter(|n| *n <= MAX_ITEMS)
            .ok_or_else(too_large)?;
        n_edges
            .checked_add(windings)
            .filter(|n| *n <= MAX_ITEMS)
            .ok_or_else(too_large)?;
        // The default vector node store uses one hedge bitset per node.
        let node_storage_bytes = output_nodes * output_hedges.div_ceil(8);
        if estimated_bytes
            .checked_add(node_storage_bytes)
            .is_none_or(|bytes| bytes > MAX_BYTES)
        {
            return Err(too_large());
        }

        let mut names: HashSet<String> = graph
            .iter_nodes()
            .filter_map(|(_, _, node)| node.name.clone())
            .collect();
        let mut edge_names = Vec::with_capacity(n_edges);
        for (_, edge, data) in graph.iter_edges() {
            let name = data.data.cut_name()?;
            names.extend(name.clone());
            edge_names.push(name.unwrap_or_else(|| format!("__linnest_cut_edge_{}", edge.0)));
        }
        let mut builder = HedgeGraphBuilder::<TypstEdge, TypstNode, TypstHedge>::new();
        let mut middle_origins = Vec::with_capacity(windings - entries.len());
        for entry in entries {
            let (source, sink, edge) = entry.endpoints(&graph)?;
            let original = graph[edge].clone();
            let orientation = graph.get_edge_data_full(source).orientation;
            graph.graph[source].route_points.clear();
            graph.graph[sink].route_points.clear();
            for segment in 0..=entry.winding {
                let name = format!("{}.{segment}", edge_names[edge.0]);
                if !names.insert(name.clone()) {
                    return Err(format!("Cut generated name collision: {name:?}"));
                }
                let mut fragment = original.cut_fragment(name)?;
                let origin = TypstCutEdgeOrigin {
                    edge: edge.0,
                    segment: Some(segment),
                    winding: entry.winding,
                };
                if segment == 0 {
                    fragment.from = Some((graph.node_id(source), source));
                    graph.graph[edge] = fragment;
                    result.edges[edge.0] = origin;
                    result
                        .boundaries
                        .push(entry.boundary(edge.0, None, source, 0));
                    continue;
                }
                if segment == entry.winding {
                    let last_edge = graph.n_edges();
                    fragment.to = Some((graph.node_id(sink), sink));
                    graph
                        .graph
                        .split_edge(sink, EdgeData::new(fragment, orientation))
                        .map_err(|err| err.to_string())?;
                    result.edges.push(origin);
                    result
                        .boundaries
                        .push(entry.boundary(last_edge, None, sink, segment - 1));
                    continue;
                }
                let middle_edge = n_edges + selected.len() + middle_origins.len();
                // Crossing endpoints carry the opposite original hedge's data:
                // a middle source continues the original sink's side of the previous
                // crossing, and its sink ends at the original source's side of the next.
                let [source_endpoint, sink_endpoint] =
                    [("source", sink), ("sink", source)].map(|(role, origin)| {
                        let name = format!("__linnest_cut_{}_{segment}_{role}", edge.0);
                        if !names.insert(name.clone()) {
                            return Err(format!("Cut generated name collision: {name:?}"));
                        }
                        let node = builder.add_node(TypstNode {
                            name: Some(name),
                            ..TypstNode::default()
                        });
                        let output_node = NodeIndex(result.nodes.len());
                        let output_hedge = Hedge(result.hedges.len());
                        result.nodes.push(None);
                        result.hedges.push(Some(origin.0));
                        let endpoint = Some((output_node, output_hedge));
                        let crossing = if role == "source" {
                            fragment.from = endpoint;
                            segment - 1
                        } else {
                            fragment.to = endpoint;
                            segment
                        };
                        result.boundaries.push(entry.boundary(
                            middle_edge,
                            Some(output_node.0),
                            origin,
                            crossing,
                        ));
                        let mut data = graph[origin].clone();
                        data.id = None;
                        data.from = 0;
                        data.to = 0;
                        Ok(HedgeData {
                            node,
                            data,
                            is_in_subgraph: false,
                        })
                    });
                builder.add_edge(source_endpoint?, sink_endpoint?, fragment, orientation);
                middle_origins.push(origin);
            }
        }
        if !middle_origins.is_empty() {
            graph
                .graph
                .append_disconnected_mut(builder.build())
                .map_err(|err| err.to_string())?;
        }
        result.edges.extend(middle_origins);
        result.graph = encode_typst_graph(&graph)?;
        encode_cbor(&result)
    }
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
                .and_then(|placement| placement.point)
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
    graph
        .apply_explicit_id_ordering()
        .map_err(|error| error.to_string())?;
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
    let has_statement_position = statement_map_value(&edge.statements, "pos").is_some()
        || statement_map_value(&global_data.edge_statements, "pos").is_some();
    let placement = resolve_placement(edge.pos.as_ref(), node_positions, "edge")?;
    let mut local_statements = apply_placement_statements(edge.statements, placement.as_ref());
    let midpoint = edge
        .source
        .as_ref()
        .zip(edge.sink.as_ref())
        .and_then(|(source, sink)| {
            if has_statement_position || placement.as_ref().is_some_and(|p| p.point.is_some()) {
                return None;
            }
            let (source, sink) = (
                node_positions.get(source.node)?,
                node_positions.get(sink.node)?,
            );
            (source.x_set && source.y_set && sink.x_set && sink.y_set).then_some(
                ResolvedPlacement {
                    point: Some(ResolvedPoint {
                        x: (source.x + sink.x) / 2.0,
                        y: (source.y + sink.y) / 2.0,
                        x_set: false,
                        y_set: false,
                    }),
                    z: None,
                    pin: None,
                    mode: PlacementMode::Start,
                    group_start_x: false,
                    group_start_y: false,
                },
            )
        });
    local_statements = apply_placement_statements(local_statements, midpoint.as_ref());
    if let Some(name) = edge.name {
        local_statements.insert(TYPST_EDGE_NAME_KEY.to_owned(), name);
    }
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

    if let Some((z, mode)) = placement.z {
        statements.remove("\"pos-z\"");
        statements.remove("\"pos-z-mode\"");
        statements.insert("pos-z".to_string(), z.to_string());
        statements.insert(
            "pos-z-mode".to_string(),
            match mode {
                PlacementMode::Start => "start",
                PlacementMode::Pin => "pin",
            }
            .to_string(),
        );
    }
    let Some(point) = placement.point else {
        return statements;
    };

    // Explicit XY placement replaces old pins; inferred midpoints and Z-only
    // patches must leave them intact.
    if point.x_set || point.y_set {
        statements.remove("pin");
        statements.remove("\"pin\"");
    }
    statements.insert("pos".to_string(), format!("{},{}", point.x, point.y));
    statements.insert("pos-x-set".to_string(), point.x_set.to_string());
    statements.insert("pos-y-set".to_string(), point.y_set.to_string());
    statements.insert(
        "pos-mode".to_string(),
        match placement.mode {
            PlacementMode::Start => "start",
            PlacementMode::Pin => "pin",
        }
        .to_string(),
    );
    statements.remove("group-start-x");
    statements.remove("group-start-y");
    if placement.group_start_x {
        statements.insert("group-start-x".to_string(), "true".to_string());
    }
    if placement.group_start_y {
        statements.insert("group-start-y".to_string(), "true".to_string());
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
        let z = self.z.map(TypstNumber::as_f64);
        if z.is_some_and(|z| !z.is_finite()) {
            return Err(format!("{context} placement z must be a finite number"));
        }
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
        let group_start_x = matches!(
            &self.x,
            Some(TypstPlacementCoord::Group(group)) if group.start.is_some() || point.x_set
        );
        let group_start_y = matches!(
            &self.y,
            Some(TypstPlacementCoord::Group(group)) if group.start.is_some() || point.y_set
        );
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
            if self.mode == PlacementMode::Pin
                && self.x_mode.is_none()
                && self.y_mode.is_none()
                && z.is_none()
            {
                return Err(format!(
                    "{context} pin placement must constrain x, y, z, or ref"
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

        Ok(ResolvedPlacement {
            // Depth-only placements must not synthesize or refresh XY state.
            point: (z.is_none() || point.x_set || point.y_set).then_some(point),
            z: z.map(|z| (z, self.z_mode.unwrap_or(self.mode))),
            pin,
            mode,
            group_start_x,
            group_start_y,
        })
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
                if let Some(start) = group.start {
                    if axis == "x" {
                        point.x = start.as_f64();
                    } else {
                        point.y = start.as_f64();
                    }
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
        pos_x_set: parse_bool_statement(&raw_statements, "pos-x-set").unwrap_or(false),
        pos_y_set: parse_bool_statement(&raw_statements, "pos-y-set").unwrap_or(false),
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
        pos_x_set: parse_bool_statement(&raw_statements, "pos-x-set").unwrap_or(false),
        pos_y_set: parse_bool_statement(&raw_statements, "pos-y-set").unwrap_or(false),
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
        "group-start-x",
        "group-start-y",
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
        route_points: Vec::new(),
    }
}

fn edges_with_route_points(graph: &TypstGraph, mut edges: Vec<TypstDotEdge>) -> Vec<TypstDotEdge> {
    for edge in &mut edges {
        if let Some(source) = &mut edge.source {
            source.route_points = route_points_for_hedge(graph, source.hedge);
        }
        if let Some(sink) = &mut edge.sink {
            sink.route_points = route_points_for_hedge(graph, sink.hedge);
        }
    }
    edges
}

fn route_points_for_hedge(graph: &TypstGraph, hedge: usize) -> Vec<TypstPoint> {
    graph.graph[Hedge(hedge)]
        .route_points
        .iter()
        .map(|point| TypstPoint {
            x: point.x,
            y: point.y,
        })
        .collect()
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
        "_" => Ok(Some(CompassPt::Underscore)),
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
        9 => "_",
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

#[cfg(test)]
mod tests {
    use super::*;
    use ciborium::Value;

    fn cut_fixture(orientation: Orientation) -> TypstGraph {
        let mut builder = DotBuilder::new();
        for (index, name) in ["b", "c", "isolated"].into_iter().enumerate() {
            builder.add_node(
                TypstNode {
                    name: Some(name.into()),
                    index: Some(NodeIndex(index)),
                    data: Some(vec![index as u8, 255]),
                    pos: Point2::new(index as f64, 4.0),
                    constraints: PointConstraint {
                        x: Constraint::Fixed,
                        y: Constraint::Fixed,
                    },
                    statements: BTreeMap::from([("pin".into(), "x:1,y:4".into())]),
                    ..TypstNode::default()
                }
                .to_dot(),
            );
        }
        for (index, name) in ["k", "uncut"].into_iter().enumerate() {
            let mut endpoints = [HedgeData::from(NodeIndex(0)), HedgeData::from(NodeIndex(1))];
            for (side, endpoint) in endpoints.iter_mut().enumerate() {
                endpoint.data = TypstHedge {
                    id: Some(2 * index + side),
                    data: Some(vec![index as u8, side as u8, 255]),
                    statement: Some(format!("weight={}", side + 2)),
                    port_label: Some(format!("port-{side}")),
                    compasspt: Some(if side == 0 { "ne" } else { "sw" }.into()),
                    ..TypstHedge::default()
                }
                .to_dot();
            }
            let [source, sink] = endpoints;
            let payload = BTreeMap::from([
                ("name", Value::Text(name.into())),
                ("data-key", Value::Integer((17 + index).into())),
                (
                    "nested",
                    Value::Array(vec![Value::Bytes(vec![0, 255]), Value::Bool(true)]),
                ),
            ]);
            let data = TypstEdge {
                data: Some(encode_cbor(&payload).unwrap()),
                pos: Point2::new(2.0, 3.0),
                start_x: true,
                start_y: true,
                shift: Some(Vector2::new(0.2, 0.3)),
                bend: Ok(Rad(0.4)),
                label_pos: Some(Point2::new(5.0, 6.0)),
                label_angle: Some(0.7),
                statements: BTreeMap::from([
                    (TYPST_EDGE_NAME_KEY.into(), name.into()),
                    ("pin".into(), "x:2,y:3".into()),
                    ("\"pos-z\"".into(), "5".into()),
                    ("pos-z-mode".into(), "pin".into()),
                    ("pos-mode".into(), "pin".into()),
                    ("group-start-x".into(), "true".into()),
                    ("group-start-y".into(), "true".into()),
                    ("route-points".into(), "stale".into()),
                    ("spring-length".into(), "2.5".into()),
                    ("source-color".into(), "red".into()),
                    ("sink-color".into(), "blue".into()),
                    ("custom".into(), "keep".into()),
                ]),
                ..TypstEdge::default()
            }
            .to_dot();
            builder.add_edge(source, sink, data, orientation);
        }
        builder.add_external_edge(
            NodeIndex(0),
            TypstEdge {
                data: Some(vec![255, 0, 128]),
                pos: Point2::new(7.0, 8.0),
                ..TypstEdge::default()
            }
            .to_dot(),
            Orientation::Reversed,
            Flow::Sink,
        );
        let global_data = global_data_from_parts(
            Some("cut-test".into()),
            Some(vec![128, 0, 255]),
            BTreeMap::from([
                ("eval".into(), "global".into()),
                ("custom".into(), "graph".into()),
            ]),
            BTreeMap::from([
                ("pin".into(), "x:9".into()),
                ("spring-length".into(), "2.5".into()),
            ]),
            BTreeMap::from([
                ("pin".into(), "y:8".into()),
                ("label".into(), "default".into()),
            ]),
        );
        let mut graph = typst_graph_from_dot(DotGraph {
            global_data,
            graph: builder.build(),
        });
        graph.layout_config =
            crate::LayoutConfig::from_figment(&default_figment().merge(("viewport-w", 23.0)));
        for index in 0..graph.n_hedges() {
            graph.graph[Hedge(index)].route_points = vec![Point2::new(index as f64, 1.0)];
            graph.graph[Hedge(index)].weight = index as f64 + 1.0;
        }
        graph
    }

    fn cut_result(graph: &TypstGraph, entries: &[TypstCutEntry]) -> Result<TypstCutResult, String> {
        let bytes = TypstGraph::cut_bytes(&encode_typst_graph(graph)?, &encode_cbor(&entries)?)?;
        decode_cbor(&bytes, "cut result")
    }

    #[test]
    fn weighted_cut_windings_orientations_payloads_and_origins() {
        for orientation in [
            Orientation::Default,
            Orientation::Reversed,
            Orientation::Undirected,
        ] {
            for winding in 1..=3 {
                for source_is_left in [false, true] {
                    let graph = cut_fixture(orientation);
                    let entry = TypstCutEntry {
                        left: usize::from(!source_is_left),
                        right: usize::from(source_is_left),
                        winding,
                    };
                    let result = cut_result(&graph, &[entry]).unwrap();
                    let output = decode_typst_graph(&result.graph).unwrap();
                    output.check().unwrap();
                    assert_eq!(output.n_nodes(), 3 + 2 * (winding - 1));
                    assert_eq!(output.n_edges(), 3 + winding);
                    assert_eq!(output.n_hedges(), 5 + 2 * (winding - 1));
                    assert_eq!(result.nodes.len(), output.n_nodes());
                    assert_eq!(result.edges.len(), output.n_edges());
                    assert_eq!(result.hedges.len(), output.n_hedges());
                    assert_eq!(output.name, graph.name);
                    assert_eq!(output.data, graph.data);
                    assert_eq!(output.global_eval, graph.global_eval);
                    assert_eq!(output.global_statements, graph.global_statements);
                    assert_eq!(
                        output.default_node_statements,
                        graph.default_node_statements
                    );
                    assert_eq!(
                        output.default_edge_statements,
                        graph.default_edge_statements
                    );
                    assert_eq!(
                        encode_cbor(&output.layout_config).unwrap(),
                        encode_cbor(&graph.layout_config).unwrap()
                    );
                    assert_eq!(&result.nodes[..3], &[Some(0), Some(1), Some(2)]);
                    assert_eq!(
                        &result.hedges[..5],
                        &[Some(0), Some(1), Some(2), Some(3), Some(4)]
                    );
                    for (node, neighbors, data) in output.iter_nodes() {
                        if node.0 < 3 {
                            assert_eq!(
                                encode_cbor(data).unwrap(),
                                encode_cbor(&graph[node]).unwrap()
                            );
                            if node.0 == 2 {
                                assert_eq!(neighbors.count(), 0);
                            }
                        } else {
                            assert_eq!(result.nodes[node.0], None);
                            assert_eq!(neighbors.count(), 1);
                            assert!(data.data.is_none());
                            assert!(data.statements.is_empty());
                            assert!(data.index.is_none());
                        }
                    }
                    for (pair, edge, data) in output.iter_edges() {
                        let origin = &result.edges[edge.0];
                        if let Some(segment) = origin.segment {
                            assert_eq!(origin.edge, 0);
                            assert_eq!(origin.winding, winding);
                            assert_eq!(data.orientation, orientation);
                            let name = format!("k.{segment}");
                            assert_eq!(
                                data.data.cut_name().unwrap().as_deref(),
                                Some(name.as_str())
                            );
                            assert_eq!(data.data.statements[TYPST_EDGE_NAME_KEY], name);
                            let mut expected: BTreeMap<String, Value> =
                                decode_cbor(graph[EdgeIndex(0)].data.as_ref().unwrap(), "payload")
                                    .unwrap();
                            expected.insert("name".into(), Value::Text(name));
                            let actual: BTreeMap<String, Value> =
                                decode_cbor(data.data.data.as_ref().unwrap(), "payload").unwrap();
                            assert_eq!(actual, expected);
                            assert_eq!(data.data.statements.len(), 5);
                            for key in ["spring-length", "source-color", "sink-color", "custom"] {
                                assert_eq!(
                                    data.data.statements[key],
                                    graph[EdgeIndex(0)].statements[key]
                                );
                            }
                            assert!(matches!(data.data.constraints.x, Constraint::Free));
                            assert!(matches!(data.data.constraints.y, Constraint::Free));
                            assert_eq!(data.data.pos, Point2::new(0.0, 0.0));
                            assert!(
                                !data.data.start_x
                                    && !data.data.start_y
                                    && !data.data.bend_explicit
                            );
                            assert!(data.data.bend.is_err());
                            assert!(
                                data.data.shift.is_none()
                                    && data.data.label_pos.is_none()
                                    && data.data.label_angle.is_none()
                            );
                            if segment == 0 {
                                assert_eq!(edge.0, 0);
                                assert_eq!(
                                    pair,
                                    HedgePair::Unpaired {
                                        hedge: Hedge(0),
                                        flow: Flow::Source
                                    }
                                );
                                assert_eq!(data.data.from, Some((NodeIndex(0), Hedge(0))));
                                assert_eq!(data.data.to, None);
                            } else if segment == winding {
                                assert_eq!(edge.0, 3);
                                assert_eq!(
                                    pair,
                                    HedgePair::Unpaired {
                                        hedge: Hedge(1),
                                        flow: Flow::Sink
                                    }
                                );
                                assert_eq!(data.data.from, None);
                                assert_eq!(data.data.to, Some((NodeIndex(1), Hedge(1))));
                            } else {
                                assert_eq!(edge.0, 3 + segment);
                                let source = Hedge(5 + 2 * (segment - 1));
                                let sink = Hedge(source.0 + 1);
                                assert_eq!(pair, HedgePair::Paired { source, sink });
                                assert_eq!(data.data.from, Some((output.node_id(source), source)));
                                assert_eq!(data.data.to, Some((output.node_id(sink), sink)));
                                assert_eq!(result.hedges[source.0], Some(1));
                                assert_eq!(result.hedges[sink.0], Some(0));
                            }
                        } else {
                            assert_eq!(origin.edge, edge.0);
                            assert_eq!(origin.winding, 0);
                            assert_eq!(
                                encode_cbor(data.data).unwrap(),
                                encode_cbor(&graph[edge]).unwrap()
                            );
                        }
                    }
                    for (hedge, data) in output.iter_hedges() {
                        let mut expected = graph[Hedge(result.hedges[hedge.0].unwrap())].clone();
                        if hedge.0 < 2 || hedge.0 >= 5 {
                            expected.route_points.clear();
                        }
                        if hedge.0 >= 5 {
                            expected.id = None;
                            expected.from = 0;
                            expected.to = 0;
                        }
                        assert_eq!(encode_cbor(data).unwrap(), encode_cbor(&expected).unwrap());
                        if hedge.0 < 5 {
                            assert_eq!(output.node_id(hedge), graph.node_id(hedge));
                        }
                    }
                    let source_side = if source_is_left { "left" } else { "right" };
                    let sink_side = if source_is_left { "right" } else { "left" };
                    let mut expected = vec![(0, None, 0, source_side, 0)];
                    for segment in 1..winding {
                        expected.push((
                            3 + segment,
                            Some(3 + 2 * (segment - 1)),
                            1,
                            sink_side,
                            segment - 1,
                        ));
                        expected.push((
                            3 + segment,
                            Some(4 + 2 * (segment - 1)),
                            0,
                            source_side,
                            segment,
                        ));
                    }
                    expected.push((3, None, 1, sink_side, winding - 1));
                    let boundaries: Vec<_> = result
                        .boundaries
                        .iter()
                        .map(|b| (b.edge, b.node, b.hedge, b.side.as_str(), b.crossing))
                        .collect();
                    assert_eq!(boundaries, expected);
                    assert_eq!(boundaries.len(), 2 * winding);
                    let edges: Vec<TypstDotEdge> =
                        decode_cbor(&graph_edges_bytes(&result.graph).unwrap(), "edges").unwrap();
                    assert_eq!(edges.len(), result.edges.len());
                }
            }
        }
    }

    #[test]
    fn weighted_cut_can_be_relaid_out() {
        let graph = TypstGraph::parse(
            r#"digraph {
            a [pos="-2,0" pin="x:-2,y:0"]
            b [pos="2,0" pin="x:2,y:0"]
            a -> b [pos="0,0" pin="x:0,y:0" bend=0.3 "spring-length"=2]
        }"#,
        )
        .unwrap();
        let cut = cut_result(
            &graph,
            &[TypstCutEntry {
                left: 1,
                right: 0,
                winding: 3,
            }],
        )
        .unwrap();
        let before = decode_typst_graph(&cut.graph).unwrap();
        for algorithm in ["layered", "force"] {
            let settings = encode_cbor(&BTreeMap::from([("layout-algo", algorithm)])).unwrap();
            let bytes = crate::api::layout_parsed_graph_bytes(&cut.graph, &settings).unwrap();
            let after = decode_typst_graph(&bytes).unwrap();
            after.check().unwrap();
            assert_eq!(after.n_nodes(), before.n_nodes());
            assert_eq!(after.n_hedges(), before.n_hedges());
            assert_eq!(after.n_edges(), before.n_edges());
            for ((pair, edge, data), (original_pair, original_edge, _)) in
                after.iter_edges().zip(before.iter_edges())
            {
                assert_eq!((pair, edge), (original_pair, original_edge));
                assert!(data.data.pos.x.is_finite() && data.data.pos.y.is_finite());
                assert_eq!(data.data.statements["spring-length"], "2");
            }
        }
    }

    #[test]
    fn weighted_cut_empty_spec_and_wire_shape() {
        for graph in [
            cut_fixture(Orientation::Default),
            TypstGraph::parse("digraph { isolated }").unwrap(),
            TypstGraph::parse("digraph {}").unwrap(),
        ] {
            let bytes = encode_typst_graph(&graph).unwrap();
            let wire =
                TypstGraph::cut_bytes(&bytes, &encode_cbor(&Vec::<TypstCutEntry>::new()).unwrap())
                    .unwrap();
            let fields: BTreeMap<String, Value> = decode_cbor(&wire, "wire fields").unwrap();
            assert_eq!(
                fields.keys().map(String::as_str).collect::<Vec<_>>(),
                ["boundaries", "edges", "graph", "hedges", "nodes"]
            );
            assert!(fields["graph"].is_array());
            let result: TypstCutResult = decode_cbor(&wire, "cut").unwrap();
            assert_eq!(result.graph, bytes);
            assert_eq!(
                result.nodes,
                (0..graph.n_nodes()).map(Some).collect::<Vec<_>>()
            );
            assert_eq!(
                result.hedges,
                (0..graph.n_hedges()).map(Some).collect::<Vec<_>>()
            );
            assert!(result.boundaries.is_empty());
            for (edge, origin) in result.edges.iter().enumerate() {
                assert_eq!(
                    origin,
                    &TypstCutEdgeOrigin {
                        edge,
                        segment: None,
                        winding: 0
                    }
                );
            }
        }
    }

    #[test]
    fn weighted_cut_mixed_entries_and_immediate_origins() {
        let graph = cut_fixture(Orientation::Reversed);
        let result = cut_result(
            &graph,
            &[
                TypstCutEntry {
                    left: 2,
                    right: 3,
                    winding: 3,
                },
                TypstCutEntry {
                    left: 1,
                    right: 0,
                    winding: 2,
                },
            ],
        )
        .unwrap();
        let output = decode_typst_graph(&result.graph).unwrap();
        output.check().unwrap();
        let origins: Vec<_> = result
            .edges
            .iter()
            .map(|o| (o.edge, o.segment, o.winding))
            .collect();
        assert_eq!(
            origins,
            [
                (0, Some(0), 2),
                (1, Some(0), 3),
                (2, None, 0),
                (1, Some(3), 3),
                (0, Some(2), 2),
                (1, Some(1), 3),
                (1, Some(2), 3),
                (0, Some(1), 2),
            ]
        );
        assert_eq!(
            &result.hedges[5..],
            &[Some(3), Some(2), Some(3), Some(2), Some(1), Some(0)]
        );
        assert_eq!(result.boundaries.len(), 10);
        let pair = output[EdgeIndex(5)].from.unwrap().1;
        let recut = cut_result(
            &output,
            &[TypstCutEntry {
                left: pair.0,
                right: output.inv(pair).0,
                winding: 1,
            }],
        )
        .unwrap();
        let recut_graph = decode_typst_graph(&recut.graph).unwrap();
        assert_eq!(recut.edges[5].edge, 5);
        assert_eq!(recut.edges.last().unwrap().edge, 5);
        assert_eq!(
            recut_graph[EdgeIndex(5)].cut_name().unwrap().as_deref(),
            Some("uncut.1.0")
        );
        assert_eq!(
            recut.hedges,
            (0..output.n_hedges()).map(Some).collect::<Vec<_>>()
        );
        assert_eq!(
            recut.nodes,
            (0..output.n_nodes()).map(Some).collect::<Vec<_>>()
        );
    }

    #[test]
    fn archived_subgraph_size_requires_exact_mask_length() {
        let graph = cut_fixture(Orientation::Default);
        let graph_bytes = encode_typst_graph(&graph).unwrap();
        let n_hedges = graph.n_hedges();
        for size in [0, 1, n_hedges - 1, n_hedges, n_hedges + 1, 64, 65] {
            let mask = SuBitGraph::empty(size);
            let bytes = encode_subgraph(&mask).unwrap();
            let size_bytes = subgraph_size_bytes(&bytes).unwrap();
            assert_eq!(
                decode_cbor::<usize>(&size_bytes, "subgraph size").unwrap(),
                size
            );
            for read in [
                graph_nodes_of_archived_subgraph_bytes,
                graph_edges_of_archived_subgraph_bytes,
            ] {
                let result = read(&graph_bytes, &bytes);
                if size == n_hedges {
                    let selected: Vec<Value> = decode_cbor(&result.unwrap(), "selection").unwrap();
                    assert!(selected.is_empty());
                } else {
                    assert_eq!(result.unwrap_err(), format!(
                        "Archived subgraph has {size} bits, but graph has {n_hedges} half-edges; sizes must match"
                    ));
                }
            }
        }
        assert!(subgraph_size_bytes(&[255]).is_err());
    }

    #[test]
    fn weighted_cut_rejects_trailing_cbor() {
        let graph = encode_typst_graph(&cut_fixture(Orientation::Default)).unwrap();
        let specs = [
            encode_cbor(&Vec::<TypstCutEntry>::new()).unwrap(),
            encode_cbor(&[TypstCutEntry {
                left: 0,
                right: 1,
                winding: 2,
            }])
            .unwrap(),
            vec![0x9f, 0xff],
        ];
        for spec in specs {
            assert!(TypstGraph::cut_bytes(&graph, &spec).is_ok());
            for trailing in [0x00, 0x80, 0xf6, 0xff] {
                let mut bytes = spec.clone();
                bytes.push(trailing);
                assert_eq!(
                    TypstGraph::cut_bytes(&graph, &bytes).unwrap_err(),
                    "Trailing bytes after cut entries CBOR array"
                );
            }
        }
    }

    #[test]
    fn weighted_cut_invalid_specs_are_errors() {
        let graph = cut_fixture(Orientation::Default);
        for (left, right, winding, message) in [
            (0, 1, 0, "positive"),
            (0, 0, 1, "paired inverses"),
            (0, 3, 1, "paired inverses"),
            (4, 4, 1, "paired inverses"),
            (0, 5, 1, "out of bounds"),
            (usize::MAX, 1, 1, "out of bounds"),
            (0, 1, usize::MAX, "limits"),
            (0, 1, 65_536, "limits"),
            (0, 1, 20_000, "limits"),
        ] {
            let error = cut_result(
                &graph,
                &[TypstCutEntry {
                    left,
                    right,
                    winding,
                }],
            )
            .unwrap_err();
            assert!(error.contains(message), "{error}");
        }
        let error = cut_result(
            &graph,
            &[
                TypstCutEntry {
                    left: 0,
                    right: 1,
                    winding: 1,
                },
                TypstCutEntry {
                    left: 1,
                    right: 0,
                    winding: 2,
                },
            ],
        )
        .unwrap_err();
        assert!(error.contains("Duplicate"), "{error}");
        let bytes = encode_typst_graph(&graph).unwrap();
        for value in [
            Value::Null,
            Value::Map(vec![]),
            Value::Array(vec![Value::Null]),
            Value::Array(vec![Value::Map(vec![])]),
        ] {
            assert!(TypstGraph::cut_bytes(&bytes, &encode_cbor(&value).unwrap()).is_err());
        }
        for value in [
            Value::Integer((-1).into()),
            Value::Float(1.5),
            Value::Text("2".into()),
            Value::Null,
        ] {
            let entries = [BTreeMap::from([
                ("left", Value::Integer(0.into())),
                ("right", Value::Integer(1.into())),
                ("winding", value),
            ])];
            assert!(TypstGraph::cut_bytes(&bytes, &encode_cbor(&entries).unwrap()).is_err());
        }
        assert!(
            TypstGraph::cut_bytes(&bytes, &[0x9b, 255, 255, 255, 255, 255, 255, 255, 255]).is_err()
        );
        assert!(TypstGraph::cut_bytes(&bytes, &[255]).is_err());
        assert!(TypstGraph::cut_bytes(&[255], &[0x80]).is_err());
    }

    #[test]
    fn weighted_cut_rejects_malformed_archived_topology() {
        for field in ["source_idx", "sink_idx", "source", "sink", "edge-index"] {
            let graph = cut_fixture(Orientation::Default);
            let mut value: Value = decode_cbor(&encode_cbor(&graph).unwrap(), "graph").unwrap();
            let mut pending = vec![&mut value];
            let mut changed = false;
            while let Some(value) = pending.pop() {
                match value {
                    Value::Map(fields) => {
                        let edge_data = fields
                            .iter()
                            .any(|(key, _)| key.as_text() == Some("orientation"));
                        for (key, value) in fields {
                            if key.as_text() == Some(field)
                                || (field == "edge-index"
                                    && edge_data
                                    && key.as_text() == Some("data")
                                    && value.is_integer())
                            {
                                *value = Value::Integer(u64::MAX.into());
                                changed = true;
                            } else {
                                pending.push(value);
                            }
                        }
                    }
                    Value::Array(values) => pending.extend(values),
                    _ => {}
                }
            }
            assert!(changed, "No serialized {field} found");
            let graph: TypstGraph =
                decode_cbor(&encode_cbor(&value).unwrap(), "corrupt graph").unwrap();
            assert!(cut_result(
                &graph,
                &[TypstCutEntry {
                    left: 0,
                    right: 1,
                    winding: 2
                }]
            )
            .is_err());
        }
        let mut graph = cut_fixture(Orientation::Default);
        graph.graph.node_store = HedgeGraphBuilder::<TypstEdge, TypstNode, TypstHedge>::new()
            .build::<DefaultNodeStore<TypstNode>>()
            .node_store;
        assert!(cut_result(
            &graph,
            &[TypstCutEntry {
                left: 0,
                right: 1,
                winding: 2
            }]
        )
        .is_err());
    }

    #[test]
    fn weighted_cut_names_do_not_overwrite_user_names() {
        for collision in 0..5 {
            let mut graph = cut_fixture(Orientation::Default);
            match collision {
                0 => {
                    graph.graph[NodeIndex(2)].name = Some("k.2".into());
                }
                1 => {
                    graph.graph[NodeIndex(2)].name = Some("__linnest_cut_0_1_source".into());
                }
                2 => {
                    graph.graph[EdgeIndex(1)].data = None;
                    graph.graph[EdgeIndex(1)]
                        .statements
                        .insert(TYPST_EDGE_NAME_KEY.into(), "k.0".into());
                }
                3 => {
                    graph.graph[EdgeIndex(1)]
                        .statements
                        .remove(TYPST_EDGE_NAME_KEY);
                    graph.graph[EdgeIndex(1)].data =
                        Some(encode_cbor(&BTreeMap::from([("name", "k.1")])).unwrap());
                }
                _ => {
                    graph.graph[EdgeIndex(0)].data =
                        Some(encode_cbor(&BTreeMap::from([("name", "different")])).unwrap());
                }
            }
            let error = cut_result(
                &graph,
                &[TypstCutEntry {
                    left: 0,
                    right: 1,
                    winding: 2,
                }],
            )
            .unwrap_err();
            assert!(
                error.contains(if collision == 4 {
                    "Conflicting"
                } else {
                    "collision"
                }),
                "{error}"
            );
        }
        for payload_only in [false, true] {
            let mut graph = cut_fixture(Orientation::Default);
            if payload_only {
                graph.graph[EdgeIndex(0)]
                    .statements
                    .remove(TYPST_EDGE_NAME_KEY);
            } else {
                graph.graph[EdgeIndex(0)].data = Some(vec![255, 0, 128]);
            }
            let result = cut_result(
                &graph,
                &[TypstCutEntry {
                    left: 0,
                    right: 1,
                    winding: 2,
                }],
            )
            .unwrap();
            let output = decode_typst_graph(&result.graph).unwrap();
            for edge in [0, 3, 4] {
                assert!(output[EdgeIndex(edge)]
                    .cut_name()
                    .unwrap()
                    .unwrap()
                    .starts_with("k."));
                if !payload_only {
                    assert_eq!(output[EdgeIndex(edge)].data, graph[EdgeIndex(0)].data);
                }
            }
        }
    }

    #[test]
    fn weighted_cut_preserves_opaque_and_unnamed_payloads() {
        let mut trailing = encode_cbor(&BTreeMap::from([("name", "opaque-prefix")])).unwrap();
        trailing.push(255);
        for payload in [
            None,
            Some(vec![255]),
            Some(trailing),
            Some(encode_cbor(&vec![1, 2, 3]).unwrap()),
        ] {
            let mut graph = cut_fixture(Orientation::Undirected);
            graph.graph[EdgeIndex(0)]
                .statements
                .remove(TYPST_EDGE_NAME_KEY);
            graph.graph[EdgeIndex(0)].data = payload.clone();
            let result = cut_result(
                &graph,
                &[TypstCutEntry {
                    left: 0,
                    right: 1,
                    winding: 2,
                }],
            )
            .unwrap();
            let output = decode_typst_graph(&result.graph).unwrap();
            for edge in [0, 3, 4] {
                let segment = result.edges[edge].segment.unwrap();
                assert_eq!(
                    output[EdgeIndex(edge)].cut_name().unwrap(),
                    Some(format!("__linnest_cut_edge_0.{segment}"))
                );
                assert_eq!(output[EdgeIndex(edge)].data, payload);
            }
        }
    }

    #[test]
    fn weighted_cut_self_loop_and_reordered_hedges() {
        let mut dot = TypstGraph::parse("digraph { a -> a }")
            .unwrap()
            .to_dot_graph();
        dot.graph[Hedge(0)].id = Some(Hedge(1));
        dot.graph[Hedge(1)].id = Some(Hedge(0));
        dot.apply_explicit_id_ordering().unwrap();
        let graph = typst_graph_from_dot(dot);
        let (pair, edge, _) = graph.iter_edges().next().unwrap();
        let HedgePair::Paired { source, sink } = pair else {
            panic!("Expected self-loop");
        };
        assert_eq!((source, sink), (Hedge(1), Hedge(0)));
        let result = cut_result(
            &graph,
            &[TypstCutEntry {
                left: sink.0,
                right: source.0,
                winding: 3,
            }],
        )
        .unwrap();
        let output = decode_typst_graph(&result.graph).unwrap();
        output.check().unwrap();
        assert_eq!(output[edge].from, Some((NodeIndex(0), source)));
        assert_eq!(output[EdgeIndex(1)].to, Some((NodeIndex(0), sink)));
        assert_eq!(result.nodes[0], Some(0));
        assert!(output[edge].cut_name().unwrap().unwrap().ends_with(".0"));
    }

    #[test]
    fn auxiliary_z_modes_preserve_xy_statements() {
        let mut spec: TypstPlacementSpec = decode_cbor(
            &encode_cbor(&BTreeMap::from([("z", -2.5)])).unwrap(),
            "placement",
        )
        .unwrap();
        let xy = BTreeMap::from([
            ("pos".into(), "3,4".into()),
            ("pos-x-set".into(), "true".into()),
            ("pos-y-set".into(), "false".into()),
            ("pos-mode".into(), "pin".into()),
            ("pin".into(), "x:@column".into()),
            ("group-start-x".into(), "true".into()),
            ("label".into(), "unchanged".into()),
            ("\"pos-z\"".into(), "9".into()),
            ("\"pos-z-mode\"".into(), "start".into()),
        ]);
        for mode in [PlacementMode::Pin, PlacementMode::Start] {
            for z_mode in [None, Some(PlacementMode::Pin), Some(PlacementMode::Start)] {
                spec.mode = mode;
                spec.z_mode = z_mode;
                let decoded: TypstPlacementSpec =
                    decode_cbor(&encode_cbor(&spec).unwrap(), "placement").unwrap();
                assert_eq!(decoded, spec);
                let resolved = decoded.resolve(&[], "test").unwrap();
                assert!(resolved.point.is_none());
                assert!(resolved.pin.is_none());
                let mut statements = apply_placement_statements(xy.clone(), Some(&resolved));
                assert_eq!(statements.remove("pos-z").as_deref(), Some("-2.5"));
                assert_eq!(
                    statements.remove("pos-z-mode").as_deref(),
                    Some(match z_mode.unwrap_or(mode) {
                        PlacementMode::Pin => "pin",
                        PlacementMode::Start => "start",
                    })
                );
                let mut expected = xy.clone();
                expected.remove("\"pos-z\"");
                expected.remove("\"pos-z-mode\"");
                assert_eq!(statements, expected);
            }
        }

        let statements = apply_placement_statements(xy, Some(&spec.resolve(&[], "test").unwrap()));
        spec.z = None;
        spec.x = Some(TypstPlacementCoord::Number(TypstNumber::Signed(7)));
        let statements =
            apply_placement_statements(statements, Some(&spec.resolve(&[], "test").unwrap()));
        assert_eq!(statements["pos-z"], "-2.5");
        assert_eq!(statements["pos-z-mode"], "start");
        assert_eq!(statements["pos-x-set"], "true");
        assert_eq!(statements["pos-y-set"], "false");
    }

    #[test]
    fn auxiliary_z_requires_finite_numeric_coordinates() {
        for z in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let spec: TypstPlacementSpec = decode_cbor(
                &encode_cbor(&BTreeMap::from([("z", z)])).unwrap(),
                "placement",
            )
            .unwrap();
            assert!(spec
                .resolve(&[], "test")
                .unwrap_err()
                .contains("finite number"));
        }
        for z in [
            Value::Integer((-3).into()),
            Value::Integer(u64::MAX.into()),
            Value::Float(f64::MAX),
        ] {
            let spec: TypstPlacementSpec = decode_cbor(
                &encode_cbor(&BTreeMap::from([("z", z)])).unwrap(),
                "placement",
            )
            .unwrap();
            let statements = apply_placement_statements(
                BTreeMap::new(),
                Some(&spec.resolve(&[], "test").unwrap()),
            );
            assert!(statements["pos-z"].parse::<f64>().unwrap().is_finite());
            assert_eq!(statements["pos-z-mode"], "pin");
            assert_eq!(statements.len(), 2);
        }
        for (key, value) in [
            ("z", Value::Text("2".into())),
            (
                "z",
                Value::Map(vec![(
                    Value::Text("kind".into()),
                    Value::Text("group".into()),
                )]),
            ),
            ("dz", Value::Integer(1.into())),
            ("ref-depth", Value::Integer(0.into())),
            ("z-mode", Value::Text("group".into())),
        ] {
            assert!(decode_cbor::<TypstPlacementSpec>(
                &encode_cbor(&BTreeMap::from([(key, value)])).unwrap(),
                "placement",
            )
            .is_err());
        }
    }

    #[test]
    fn auxiliary_z_build_keeps_defaults_and_edge_midpoint_seeds() {
        let mut spec: TypstGraphSpec = decode_cbor(
            &encode_cbor(&BTreeMap::from([(
                "nodes",
                vec![
                    BTreeMap::from([("name", "a")]),
                    BTreeMap::from([("name", "b")]),
                ],
            )]))
            .unwrap(),
            "graph",
        )
        .unwrap();
        spec.nodes[1].pos = Some(
            decode_cbor(
                &encode_cbor(&BTreeMap::from([("x", 8), ("y", 10)])).unwrap(),
                "placement",
            )
            .unwrap(),
        );
        spec.edges.push(
            decode_cbor(
                &encode_cbor(&BTreeMap::from([
                    ("source", BTreeMap::from([("node", 0)])),
                    ("sink", BTreeMap::from([("node", 1)])),
                ]))
                .unwrap(),
                "edge",
            )
            .unwrap(),
        );
        let z: TypstPlacementSpec = decode_cbor(
            &encode_cbor(&BTreeMap::from([("z", 3)])).unwrap(),
            "placement",
        )
        .unwrap();

        for with_defaults in [false, true] {
            if with_defaults {
                spec.default_node_statements = BTreeMap::from([
                    ("pos".into(), "2,4".into()),
                    ("pos-z".into(), "9".into()),
                    ("pos-z-mode".into(), "start".into()),
                ]);
            }
            for edge_pos in [None, Some("11,12")] {
                spec.default_edge_statements = edge_pos
                    .map(|pos| BTreeMap::from([("pos".into(), pos.into())]))
                    .unwrap_or_default();
                let baseline = graph_from_spec(spec.clone()).unwrap();
                let mut placed = spec.clone();
                placed.nodes[0].pos = Some(z.clone());
                placed.edges[0].pos = Some(z.clone());
                let graph = graph_from_spec(placed.clone()).unwrap();
                for (actual, expected) in [
                    (
                        &graph.graph[NodeIndex(0)].statements,
                        &baseline.graph[NodeIndex(0)].statements,
                    ),
                    (
                        &graph.graph[EdgeIndex(0)].statements,
                        &baseline.graph[EdgeIndex(0)].statements,
                    ),
                ] {
                    let mut expected = expected.clone();
                    expected.insert("pos-z".into(), "3".into());
                    expected.insert("pos-z-mode".into(), "pin".into());
                    assert_eq!(actual, &expected);
                }
                if with_defaults && edge_pos.is_none() {
                    assert_eq!(graph.graph[EdgeIndex(0)].statements["pos"], "5,7");
                    assert_eq!(graph.graph[EdgeIndex(0)].statements["pos-x-set"], "false");
                    assert_eq!(graph.graph[EdgeIndex(0)].statements["pos-y-set"], "false");
                }
                assert_eq!(
                    graph.graph[NodeIndex(1)].statements,
                    baseline.graph[NodeIndex(1)].statements
                );
                let bytes =
                    graph_from_spec_bytes(&encode_graph_spec_bytes(&placed).unwrap()).unwrap();
                let nodes: Vec<TypstDotNode> =
                    decode_cbor(&graph_nodes_bytes(&bytes).unwrap(), "nodes").unwrap();
                let edges: Vec<TypstDotEdge> =
                    decode_cbor(&graph_edges_bytes(&bytes).unwrap(), "edges").unwrap();
                assert_eq!(nodes[0].statements["pos-z"], "3");
                assert_eq!(edges[0].statements["pos-z-mode"], "pin");
                if with_defaults && edge_pos.is_none() {
                    assert_eq!(edges[0].pos, Some(TypstPoint { x: 5.0, y: 7.0 }));
                    assert!(!edges[0].pos_x_set && !edges[0].pos_y_set);
                }
            }
        }
    }

    #[test]
    fn auxiliary_z_structural_patches_preserve_xy_state_and_metadata() {
        for dot in [
            "digraph { a; b; a -> b; }",
            "digraph { a [pos=\"x:4!\"]; b [pos=\"x:@column!,y:8\"]; a -> b [pos=\"x:3,y:5!\"]; }",
        ] {
            let graphs = crate::parse_dot_graphs_bytes(dot.as_bytes()).unwrap();
            let graph = decode_graph_bytes_list(&graphs).unwrap().remove(0);
            let before = decode_typst_graph(&graph).unwrap();
            let z_patch = BTreeMap::from([
                ("index", Value::Integer(0.into())),
                (
                    "pos",
                    Value::Map(vec![
                        (Value::Text("z".into()), Value::Integer(6.into())),
                        (Value::Text("z-mode".into()), Value::Text("start".into())),
                    ]),
                ),
            ]);
            let patch =
                BTreeMap::from([("nodes", vec![z_patch.clone()]), ("edges", vec![z_patch])]);
            let patched =
                graph_apply_structural_patches_bytes(&graph, &encode_cbor(&patch).unwrap())
                    .unwrap();
            let after = decode_typst_graph(&patched).unwrap();
            for (actual, expected) in [
                (
                    &after.graph[NodeIndex(0)].statements,
                    &before.graph[NodeIndex(0)].statements,
                ),
                (
                    &after.graph[EdgeIndex(0)].statements,
                    &before.graph[EdgeIndex(0)].statements,
                ),
            ] {
                let mut expected = expected.clone();
                expected.insert("pos-z".into(), "6".into());
                expected.insert("pos-z-mode".into(), "start".into());
                assert_eq!(actual, &expected);
            }
            assert_eq!(
                after.graph[NodeIndex(0)].pos,
                before.graph[NodeIndex(0)].pos
            );
            assert_eq!(
                after.graph[EdgeIndex(0)].pos,
                before.graph[EdgeIndex(0)].pos
            );
            assert_eq!(
                encode_cbor(&after.graph[NodeIndex(0)].constraints).unwrap(),
                encode_cbor(&before.graph[NodeIndex(0)].constraints).unwrap(),
            );
            assert_eq!(
                encode_cbor(&after.graph[EdgeIndex(0)].constraints).unwrap(),
                encode_cbor(&before.graph[EdgeIndex(0)].constraints).unwrap(),
            );
            let old_nodes: Vec<TypstDotNode> =
                decode_cbor(&graph_nodes_bytes(&graph).unwrap(), "nodes").unwrap();
            let mut nodes: Vec<TypstDotNode> =
                decode_cbor(&graph_nodes_bytes(&patched).unwrap(), "nodes").unwrap();
            let old_edges: Vec<TypstDotEdge> =
                decode_cbor(&graph_edges_bytes(&graph).unwrap(), "edges").unwrap();
            let mut edges: Vec<TypstDotEdge> =
                decode_cbor(&graph_edges_bytes(&patched).unwrap(), "edges").unwrap();
            for statements in [&mut nodes[0].statements, &mut edges[0].statements] {
                assert_eq!(statements.remove("pos-z").as_deref(), Some("6"));
                assert_eq!(statements.remove("pos-z-mode").as_deref(), Some("start"));
            }
            assert_eq!(nodes, old_nodes);
            assert_eq!(edges, old_edges);

            let xy_patch = BTreeMap::from([
                ("index", Value::Integer(0.into())),
                (
                    "pos",
                    Value::Map(vec![(Value::Text("x".into()), Value::Integer(7.into()))]),
                ),
            ]);
            let patch =
                BTreeMap::from([("nodes", vec![xy_patch.clone()]), ("edges", vec![xy_patch])]);
            let patched =
                graph_apply_structural_patches_bytes(&patched, &encode_cbor(&patch).unwrap())
                    .unwrap();
            let after = decode_typst_graph(&patched).unwrap();
            for statements in [
                &after.graph[NodeIndex(0)].statements,
                &after.graph[EdgeIndex(0)].statements,
            ] {
                assert_eq!(statements["pos-z"], "6");
                assert_eq!(statements["pos-z-mode"], "start");
                assert_eq!(statements["pos-x-set"], "true");
                assert_eq!(statements["pos-y-set"], "false");
            }
        }
    }

    #[test]
    fn auxiliary_z_patch_keeps_xy_references_in_the_same_batch() {
        let graphs =
            crate::parse_dot_graphs_bytes(b"digraph { a [pos=\"4,5\"]; b; a -> b; }").unwrap();
        let graph = decode_graph_bytes_list(&graphs).unwrap().remove(0);
        let patch = BTreeMap::from([(
            "nodes",
            vec![
                BTreeMap::from([
                    ("index", Value::Integer(0.into())),
                    (
                        "pos",
                        Value::Map(vec![(Value::Text("z".into()), Value::Integer(9.into()))]),
                    ),
                ]),
                BTreeMap::from([
                    ("index", Value::Integer(1.into())),
                    (
                        "pos",
                        Value::Map(vec![
                            (Value::Text("ref".into()), Value::Integer(0.into())),
                            (Value::Text("dx".into()), Value::Integer(2.into())),
                            (Value::Text("z".into()), Value::Integer((-3).into())),
                        ]),
                    ),
                ]),
            ],
        )]);
        let patched =
            graph_apply_structural_patches_bytes(&graph, &encode_cbor(&patch).unwrap()).unwrap();
        let nodes: Vec<TypstDotNode> =
            decode_cbor(&graph_nodes_bytes(&patched).unwrap(), "nodes").unwrap();
        assert_eq!(nodes[0].pos, Some(TypstPoint { x: 4.0, y: 5.0 }));
        assert_eq!(nodes[1].pos, Some(TypstPoint { x: 6.0, y: 5.0 }));
        assert_eq!(nodes[0].statements["pos-z"], "9");
        assert_eq!(nodes[1].statements["pos-z"], "-3");
    }
}
