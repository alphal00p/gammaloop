mod api;
pub mod geom;
mod graph_api;
mod pin;
#[cfg(test)]
mod tests;
#[cfg(feature = "custom")]
mod wasm_random;

pub use api::{
    layout_graph_bytes, layout_parsed_graph_bytes, layout_parsed_graphs_bytes,
    parse_dot_graphs_bytes,
};
pub use graph_api::{
    graph_apply_structural_patches_bytes, graph_archived_compass_subgraph_bytes,
    graph_archived_subgraph_bytes, graph_compass_subgraph_bytes, graph_cycle_basis_bytes,
    graph_dot_bytes, graph_edge_data_by_name_bytes, graph_edges_bytes,
    graph_edges_of_archived_subgraph_bytes, graph_edges_of_bytes, graph_from_spec_bytes,
    graph_info_bytes, graph_join_by_edge_key_bytes, graph_join_by_hedge_key_bytes,
    graph_node_data_by_name_bytes, graph_nodes_bytes, graph_nodes_of_archived_subgraph_bytes,
    graph_nodes_of_bytes, graph_set_edge_data_by_name_bytes, graph_set_node_data_by_name_bytes,
    graph_spanning_forests_bytes, graph_subgraph_bytes, graph_with_data_bytes,
    subgraph_contains_hedge_bytes, subgraph_hedges_bytes, subgraph_label_bytes, TypstDotEdge,
    TypstDotEndpoint, TypstDotGraphInfo, TypstDotNode, TypstPoint,
};
pub use pin::PinConstraint;

use cgmath::{EuclideanSpace, InnerSpace, Point2, Rad, Vector2, Zero};
use dot_parser::ast::CompassPt;
use figment::{providers::Serialized, Figment, Profile};
use linnet::half_edge::swap::Swap;
use linnet::{
    half_edge::{
        involution::{EdgeData, EdgeIndex, EdgeVec, Flow, Hedge, HedgePair, Involution},
        layout::{
            force::{force_directed_layout, ForceLayoutConfig},
            layered::{
                LayeredConfig, LayeredEdgeRoute, LayeredGeometry, LayeredOutput, LayeredProfile,
                LayeredRouteExit,
            },
            simulatedanneale::{anneal, GeoSchedule, SAConfig},
            spring::{
                Constraint, HasPointConstraint, LayoutState, ParamTuning, PinnedLayoutNeighbor,
                PointConstraint, ShiftDirection, Shiftable, SpringChargeEnergy,
            },
        },
        nodestore::{DefaultNodeStore, NodeStorageOps},
        subgraph::{SuBitGraph, SubSetLike},
        EdgeAccessors, HedgeGraph, NodeIndex, NodeVec,
    },
    parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData, GlobalData, HedgeParseError},
};

use rand::rngs::SmallRng;
use rkyv::with::{ArchiveWith, DeserializeWith, SerializeWith};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::BufWriter;
use std::ops::{Deref, IndexMut};
use std::path::Path;
use std::{f64, fs::File};

#[cfg(target_arch = "wasm32")]
use wasm_minimal_protocol::*;

#[cfg(target_arch = "wasm32")]
initiate_protocol!();

// Custom getrandom implementation for WASM
#[cfg(feature = "custom")]
use getrandom::register_custom_getrandom;
#[cfg(feature = "custom")]
use wasm_random::custom_getrandom;

use crate::geom::{tangent_angle_toward_c_side, GeomError};

#[cfg(feature = "custom")]
register_custom_getrandom!(custom_getrandom);

fn dot_statement_value<'a>(
    statements: &'a BTreeMap<String, String>,
    key: &str,
) -> Option<&'a String> {
    statements
        .get(key)
        .or_else(|| statements.get(&format!("\"{key}\"")))
}

struct Point2Rkyv;

impl ArchiveWith<Point2<f64>> for Point2Rkyv {
    type Archived = rkyv::Archived<(f64, f64)>;
    type Resolver = rkyv::Resolver<(f64, f64)>;

    unsafe fn resolve_with(
        field: &Point2<f64>,
        pos: usize,
        resolver: Self::Resolver,
        out: *mut Self::Archived,
    ) {
        rkyv::Archive::resolve(&(field.x, field.y), pos, resolver, out);
    }
}

impl<S: rkyv::Fallible + ?Sized> SerializeWith<Point2<f64>, S> for Point2Rkyv
where
    (f64, f64): rkyv::Serialize<S>,
{
    fn serialize_with(field: &Point2<f64>, serializer: &mut S) -> Result<Self::Resolver, S::Error> {
        rkyv::Serialize::serialize(&(field.x, field.y), serializer)
    }
}

impl<D: rkyv::Fallible + ?Sized> DeserializeWith<rkyv::Archived<(f64, f64)>, Point2<f64>, D>
    for Point2Rkyv
where
    rkyv::Archived<(f64, f64)>: rkyv::Deserialize<(f64, f64), D>,
{
    fn deserialize_with(
        field: &rkyv::Archived<(f64, f64)>,
        deserializer: &mut D,
    ) -> Result<Point2<f64>, D::Error> {
        let (x, y) = rkyv::Deserialize::deserialize(field, deserializer)?;
        Ok(Point2::new(x, y))
    }
}

struct Vector2Rkyv;

impl ArchiveWith<Vector2<f64>> for Vector2Rkyv {
    type Archived = rkyv::Archived<(f64, f64)>;
    type Resolver = rkyv::Resolver<(f64, f64)>;

    unsafe fn resolve_with(
        field: &Vector2<f64>,
        pos: usize,
        resolver: Self::Resolver,
        out: *mut Self::Archived,
    ) {
        rkyv::Archive::resolve(&(field.x, field.y), pos, resolver, out);
    }
}

impl<S: rkyv::Fallible + ?Sized> SerializeWith<Vector2<f64>, S> for Vector2Rkyv
where
    (f64, f64): rkyv::Serialize<S>,
{
    fn serialize_with(
        field: &Vector2<f64>,
        serializer: &mut S,
    ) -> Result<Self::Resolver, S::Error> {
        rkyv::Serialize::serialize(&(field.x, field.y), serializer)
    }
}

impl<D: rkyv::Fallible + ?Sized> DeserializeWith<rkyv::Archived<(f64, f64)>, Vector2<f64>, D>
    for Vector2Rkyv
where
    rkyv::Archived<(f64, f64)>: rkyv::Deserialize<(f64, f64), D>,
{
    fn deserialize_with(
        field: &rkyv::Archived<(f64, f64)>,
        deserializer: &mut D,
    ) -> Result<Vector2<f64>, D::Error> {
        let (x, y) = rkyv::Deserialize::deserialize(field, deserializer)?;
        Ok(Vector2::new(x, y))
    }
}

struct BendRkyv;

impl ArchiveWith<Result<Rad<f64>, GeomError>> for BendRkyv {
    type Archived = rkyv::Archived<Result<f64, GeomError>>;
    type Resolver = rkyv::Resolver<Result<f64, GeomError>>;

    unsafe fn resolve_with(
        field: &Result<Rad<f64>, GeomError>,
        pos: usize,
        resolver: Self::Resolver,
        out: *mut Self::Archived,
    ) {
        rkyv::Archive::resolve(&field.map(|angle| angle.0), pos, resolver, out);
    }
}

impl<S: rkyv::Fallible + ?Sized> SerializeWith<Result<Rad<f64>, GeomError>, S> for BendRkyv
where
    Result<f64, GeomError>: rkyv::Serialize<S>,
{
    fn serialize_with(
        field: &Result<Rad<f64>, GeomError>,
        serializer: &mut S,
    ) -> Result<Self::Resolver, S::Error> {
        rkyv::Serialize::serialize(&field.map(|angle| angle.0), serializer)
    }
}

impl<D: rkyv::Fallible + ?Sized>
    DeserializeWith<rkyv::Archived<Result<f64, GeomError>>, Result<Rad<f64>, GeomError>, D>
    for BendRkyv
where
    rkyv::Archived<Result<f64, GeomError>>: rkyv::Deserialize<Result<f64, GeomError>, D>,
{
    fn deserialize_with(
        field: &rkyv::Archived<Result<f64, GeomError>>,
        deserializer: &mut D,
    ) -> Result<Result<Rad<f64>, GeomError>, D::Error> {
        Ok(rkyv::Deserialize::deserialize(field, deserializer)?.map(Rad))
    }
}

#[derive(
    Debug, Serialize, Deserialize, Clone, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize,
)]
pub struct TypstNode {
    name: Option<String>,
    index: Option<NodeIndex>,
    data: Option<Vec<u8>>,
    #[with(Point2Rkyv)]
    pos: Point2<f64>,
    start_x: bool,
    start_y: bool,
    constraints: PointConstraint,
    #[with(rkyv::with::Map<Vector2Rkyv>)]
    shift: Option<Vector2<f64>>,
    statements: BTreeMap<String, String>,
    eval: Option<String>,
}

impl Default for TypstNode {
    fn default() -> Self {
        TypstNode {
            name: None,
            index: None,
            data: None,
            pos: Point2::origin(),
            start_x: false,
            start_y: false,
            constraints: PointConstraint::default(),
            shift: None,
            statements: BTreeMap::new(),
            eval: None,
        }
    }
}

impl Shiftable for TypstNode {
    fn shift<I: From<usize> + PartialEq + Copy, R: std::ops::IndexMut<I, Output = Point2<f64>>>(
        &self,
        shift: Vector2<f64>,
        index: I,
        values: &mut R,
    ) -> bool {
        self.constraints.shift(shift, index, values)
    }
}

impl HasPointConstraint for TypstNode {
    fn point_constraint(&self) -> &PointConstraint {
        &self.constraints
    }
}

impl TypstNode {
    /// Convert back to DotVertexData
    fn to_dot(&self) -> DotVertexData {
        let mut statements = self.statements.clone();

        // Add position as pos attribute
        statements.insert("pos".to_string(), format!("{},{}", self.pos.x, self.pos.y));

        if let Some(s) = self.shift {
            statements.insert("shift".to_string(), format!("{},{}", s.x, s.y));
        }

        if let Some(eval) = &self.eval {
            statements.insert("eval".to_string(), eval.clone());
        }

        DotVertexData {
            name: self.name.clone(),
            index: self.index,
            payload: self.data.clone(),
            statements,
        }
    }

    fn parse(
        _inv: &Involution,
        nid: NodeIndex,
        data: DotVertexData,
        init_points: &NodeVec<(Point2<f64>, PointConstraint)>,
    ) -> Self {
        let shift =
            Self::parse_position(&data.statements, "shift").map(|(x, y)| Vector2::new(x, y));

        let eval: Option<String> = data.get("eval").transpose().unwrap();

        let (pos, constraints) = init_points[nid];
        let (start_x, start_y) = Self::parse_position_flags(&data.statements);

        Self {
            name: data.name,
            index: data.index,
            data: data.payload,
            statements: data.statements,
            pos,
            start_x,
            start_y,
            constraints,
            shift,
            eval,
        }
    }

    fn parse_position(
        statements: &std::collections::BTreeMap<String, String>,
        attr: &str,
    ) -> Option<(f64, f64)> {
        if let Some(value) = dot_statement_value(statements, attr) {
            // Remove outer quotes from DOT parsing: "\"1.0,2.0\"" -> "1.0,2.0"
            let unquoted = value.trim().trim_matches('"');
            // Parse formats like "1.0,2.0" or "1.0 2.0" or "(1.0,2.0)"
            let cleaned = unquoted.trim().trim_matches(|c| c == '(' || c == ')');
            let parts: Vec<&str> = cleaned
                .split([',', ' '])
                .filter(|s| !s.is_empty())
                .collect();

            if parts.len() == 2 {
                if let (Ok(x), Ok(y)) = (
                    parts[0].trim().parse::<f64>(),
                    parts[1].trim().parse::<f64>(),
                ) {
                    return Some((x, y));
                }
            }
        }
        None
    }

    fn parse_position_flags(
        statements: &std::collections::BTreeMap<String, String>,
    ) -> (bool, bool) {
        if dot_statement_value(statements, "pos").is_none() {
            return (false, false);
        }
        (
            Self::parse_bool_statement(statements, "pos-x-set").unwrap_or(true),
            Self::parse_bool_statement(statements, "pos-y-set").unwrap_or(true),
        )
    }

    fn parse_bool_statement(
        statements: &std::collections::BTreeMap<String, String>,
        key: &str,
    ) -> Option<bool> {
        dot_statement_value(statements, key)
            .and_then(|value| value.trim().trim_matches('"').parse::<bool>().ok())
    }
}

fn parse_rad_statement(value: &str) -> Option<Rad<f64>> {
    let unquoted = value.trim().trim_matches('"');
    let cleaned = unquoted.trim().trim_end_matches("rad");
    cleaned.trim().parse::<f64>().ok().map(Rad)
}

#[derive(
    Debug, Serialize, Deserialize, Clone, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize,
)]
pub struct TypstEdge {
    from: Option<(NodeIndex, Hedge)>,
    to: Option<(NodeIndex, Hedge)>,
    data: Option<Vec<u8>>,
    #[with(BendRkyv)]
    bend: Result<Rad<f64>, GeomError>,
    #[with(Point2Rkyv)]
    pos: Point2<f64>,
    start_x: bool,
    start_y: bool,
    #[with(rkyv::with::Map<Point2Rkyv>)]
    label_pos: Option<Point2<f64>>,
    label_angle: Option<f64>,
    #[with(rkyv::with::Map<Vector2Rkyv>)]
    shift: Option<Vector2<f64>>,
    bend_explicit: bool,
    statements: BTreeMap<String, String>,
    pub constraints: PointConstraint,
}
impl Shiftable for TypstEdge {
    fn shift<I: From<usize> + PartialEq + Copy, R: std::ops::IndexMut<I, Output = Point2<f64>>>(
        &self,
        shift: Vector2<f64>,
        index: I,
        values: &mut R,
    ) -> bool {
        self.constraints.shift(shift, index, values)
    }
}

impl HasPointConstraint for TypstEdge {
    fn point_constraint(&self) -> &PointConstraint {
        &self.constraints
    }
}
impl Default for TypstEdge {
    fn default() -> Self {
        Self {
            from: None,
            to: None,
            data: None,
            bend: Err(GeomError::NotComputed),
            pos: Point2::origin(),
            start_x: false,
            start_y: false,
            label_pos: None,
            label_angle: None,
            shift: None,
            bend_explicit: false,
            statements: BTreeMap::new(),
            constraints: PointConstraint::default(),
        }
    }
}

impl TypstEdge {
    fn parse<N: NodeStorageOps>(
        _inv: &Involution,
        node_store: &N,
        p: HedgePair,
        eid: EdgeIndex,
        data: EdgeData<DotEdgeData>,
        pin_constaints: &EdgeVec<(Point2<f64>, PointConstraint)>,
    ) -> EdgeData<Self> {
        data.map(|d| {
            let shift =
                TypstNode::parse_position(&d.statements, "shift").map(|(x, y)| Vector2::new(x, y));

            let label_pos = TypstNode::parse_position(&d.statements, "label-pos")
                .map(|(x, y)| Point2::new(x, y));
            let label_angle = dot_statement_value(&d.statements, "label-angle").and_then(|value| {
                let unquoted = value.trim().trim_matches('"');
                let cleaned = unquoted.trim().trim_end_matches("rad");
                cleaned.trim().parse::<f64>().ok()
            });
            let bend = dot_statement_value(&d.statements, "bend")
                .and_then(|value| parse_rad_statement(value))
                .ok_or(GeomError::NotComputed);
            let bend_explicit = bend.is_ok();

            let mut from = None;
            let mut to = None;
            match p {
                HedgePair::Split { source, sink, .. } | HedgePair::Paired { source, sink } => {
                    from = Some((node_store.node_id_ref(source), source));
                    to = Some((node_store.node_id_ref(sink), sink));
                }
                HedgePair::Unpaired {
                    hedge,
                    flow: Flow::Source,
                } => {
                    from = Some((node_store.node_id_ref(hedge), hedge));
                }

                HedgePair::Unpaired {
                    hedge,
                    flow: Flow::Sink,
                } => {
                    to = Some((node_store.node_id_ref(hedge), hedge));
                }
            }

            let (pos, constraints) = pin_constaints[eid];
            let (start_x, start_y) = TypstNode::parse_position_flags(&d.statements);

            Self {
                from,
                to,
                pos,
                start_x,
                start_y,
                constraints,
                label_pos,
                label_angle,
                shift,
                bend,
                bend_explicit,
                data: d.payload,
                statements: d.statements,
            }
        })
    }

    /// Convert back to DotEdgeData
    fn to_dot(&self) -> DotEdgeData {
        let mut statements = self.statements.clone();

        // Add position as pos attribute
        statements.insert("pos".to_string(), format!("{},{}", self.pos.x, self.pos.y));

        // Add shift if present
        if let Some(s) = self.shift {
            statements.insert("shift".to_string(), format!("{},{}", s.x, s.y));
        }

        // Add bend if non-default
        if let Ok(b) = self.bend {
            statements.insert("bend".to_string(), format!("{}rad", b.0));
        }

        if let Some(p) = self.label_pos {
            statements.insert("label-pos".to_string(), format!("{},{}", p.x, p.y));
        }

        if let Some(a) = self.label_angle {
            statements.insert("label-angle".to_string(), format!("{a}rad"));
        }

        DotEdgeData {
            payload: self.data.clone(),
            local_statements: statements.clone(),
            statements,
            edge_id: None,
        }
    }
}

#[derive(
    Debug, Serialize, Deserialize, Clone, Default, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize,
)]
pub struct TypstHedge {
    from: usize,
    to: usize,
    weight: f64,
    #[with(rkyv::with::Map<Point2Rkyv>)]
    route_points: Vec<Point2<f64>>,
    statement: Option<String>,
    id: Option<usize>,
    data: Option<Vec<u8>>,
    port_label: Option<String>,
    compasspt: Option<String>,
}

impl TypstHedge {
    /// Convert back to DotHedgeData
    fn to_dot(&self) -> DotHedgeData {
        let statement = if self.statement.is_some() {
            self.statement.clone()
        } else if self.weight != 0.0 {
            Some(format!("weight={}", self.weight))
        } else {
            None
        };

        DotHedgeData {
            statement,
            id: self.id.map(Hedge),
            payload: self.data.clone(),
            port_label: self.port_label.clone(),
            compasspt: self
                .compasspt
                .as_deref()
                .and_then(|value| parse_hedge_compass(value).ok())
                .flatten(),
        }
    }

    fn parse(_h: Hedge, data: DotHedgeData) -> Self {
        Self {
            statement: data.statement,
            id: data.id.map(|id| id.0),
            data: data.payload,
            port_label: data.port_label,
            compasspt: data.compasspt.map(hedge_compass_to_string),
            ..Default::default()
        }
    }
}

fn parse_hedge_compass(value: &str) -> Result<Option<CompassPt>, String> {
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

fn hedge_compass_to_string(compass: CompassPt) -> String {
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

#[derive(Debug, Serialize, Deserialize, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)]
pub struct TypstGraph {
    graph: HedgeGraph<TypstEdge, TypstNode, TypstHedge>,
    global_eval: Option<String>,
    name: String,
    data: Option<Vec<u8>>,
    global_statements: BTreeMap<String, String>,
    default_edge_statements: BTreeMap<String, String>,
    default_node_statements: BTreeMap<String, String>,
    layout_config: LayoutConfig,
}

impl Deref for TypstGraph {
    type Target = HedgeGraph<TypstEdge, TypstNode, TypstHedge>;
    fn deref(&self) -> &Self::Target {
        &self.graph
    }
}

#[derive(Debug, Clone)]
enum DotPlacementExpr {
    Point {
        point: Point2<f64>,
        pinned: bool,
    },
    Ref {
        target: DotPlacementRef,
        offset: Vector2<f64>,
        pinned: bool,
    },
    Axes {
        x: Option<DotAxisPlacement>,
        y: Option<DotAxisPlacement>,
    },
}

#[derive(Debug, Clone)]
struct DotAxisPlacement {
    coord: DotAxisCoord,
    pinned: bool,
}

#[derive(Debug, Clone)]
enum DotAxisCoord {
    Number(f64),
    Group(String),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum DotPlacementRef {
    Node(usize),
    Edge(usize),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum DotPlacementTarget {
    Node(usize),
    Edge(usize),
}

#[derive(Debug, Clone)]
struct DotResolvedPlacement {
    point: Point2<f64>,
    x_set: bool,
    y_set: bool,
    pin: Option<PinConstraint>,
}

#[derive(Clone, Copy, Debug)]
struct LayoutRect {
    min_x: f64,
    max_x: f64,
    min_y: f64,
    max_y: f64,
}

impl LayoutRect {
    fn centered(center: Point2<f64>, half_width: f64, half_height: f64) -> Self {
        Self {
            min_x: center.x - half_width,
            max_x: center.x + half_width,
            min_y: center.y - half_height,
            max_y: center.y + half_height,
        }
    }

    fn overlaps(&self, other: &Self) -> bool {
        self.max_x > other.min_x
            && other.max_x > self.min_x
            && self.max_y > other.min_y
            && other.max_y > self.min_y
    }

    fn overlap_area(&self, other: &Self) -> f64 {
        let width = (self.max_x.min(other.max_x) - self.min_x.max(other.min_x)).max(0.0);
        let height = (self.max_y.min(other.max_y) - self.min_y.max(other.min_y)).max(0.0);
        width * height
    }
}

struct DotPlacementContext {
    node_exprs: HashMap<usize, DotPlacementExpr>,
    edge_exprs: HashMap<usize, DotPlacementExpr>,
    node_id_map: HashMap<usize, usize>,
    edge_id_map: HashMap<usize, usize>,
    memo: HashMap<DotPlacementTarget, Option<DotResolvedPlacement>>,
    visiting: HashSet<DotPlacementTarget>,
}

impl DotPlacementContext {
    fn new(dot: &DotGraph) -> Self {
        let mut node_exprs = HashMap::new();
        let mut edge_exprs = HashMap::new();
        let mut node_id_map = HashMap::new();
        let mut edge_id_map = HashMap::new();

        for (node_id, _, node) in dot.graph.iter_nodes() {
            if let Some(index) = node.index {
                node_id_map.insert(index.0, node_id.0);
            }
            if let Some(expr) = dot_statement_value(&node.statements, "pos").and_then(|value| {
                DotPlacementExpr::parse(value)
                    .unwrap_or_else(|err| panic!("invalid DOT node pos {value:?}: {err}"))
            }) {
                node_exprs.insert(node_id.0, expr);
            }
        }

        for (_, edge_id, edge) in dot.graph.iter_edges() {
            if let Some(index) = edge.data.edge_id {
                edge_id_map.insert(index.0, edge_id.0);
            }
            if let Some(expr) =
                dot_statement_value(&edge.data.statements, "pos").and_then(|value| {
                    DotPlacementExpr::parse(value)
                        .unwrap_or_else(|err| panic!("invalid DOT edge pos {value:?}: {err}"))
                })
            {
                edge_exprs.insert(edge_id.0, expr);
            }
        }

        Self {
            node_exprs,
            edge_exprs,
            node_id_map,
            edge_id_map,
            memo: HashMap::new(),
            visiting: HashSet::new(),
        }
    }

    fn resolve_node(&mut self, node_id: NodeIndex) -> Option<DotResolvedPlacement> {
        self.resolve_target(DotPlacementTarget::Node(node_id.0))
            .unwrap_or_else(|err| panic!("{err}"))
    }

    fn resolve_edge(&mut self, edge_id: EdgeIndex) -> Option<DotResolvedPlacement> {
        self.resolve_target(DotPlacementTarget::Edge(edge_id.0))
            .unwrap_or_else(|err| panic!("{err}"))
    }

    fn resolve_target(
        &mut self,
        target: DotPlacementTarget,
    ) -> Result<Option<DotResolvedPlacement>, String> {
        if let Some(resolved) = self.memo.get(&target) {
            return Ok(resolved.clone());
        }

        if !self.visiting.insert(target) {
            return Err(format!("cyclic DOT pos reference involving {target:?}"));
        }

        let expr = match target {
            DotPlacementTarget::Node(index) => self.node_exprs.get(&index).cloned(),
            DotPlacementTarget::Edge(index) => self.edge_exprs.get(&index).cloned(),
        };

        let resolved = expr
            .as_ref()
            .map(|expr| self.resolve_expr(expr))
            .transpose()?;

        self.visiting.remove(&target);
        self.memo.insert(target, resolved.clone());
        Ok(resolved)
    }

    fn resolve_expr(&mut self, expr: &DotPlacementExpr) -> Result<DotResolvedPlacement, String> {
        match expr {
            DotPlacementExpr::Point { point, pinned } => Ok(DotResolvedPlacement {
                point: *point,
                x_set: true,
                y_set: true,
                pin: pinned.then_some(PinConstraint::Fixed(point.x, point.y)),
            }),
            DotPlacementExpr::Ref {
                target,
                offset,
                pinned,
            } => {
                let reference_target = self.resolve_ref(*target)?;
                let reference = self
                    .resolve_target(reference_target)?
                    .ok_or_else(|| format!("DOT pos reference {target:?} has no position"))?;
                let point = reference.point + *offset;
                Ok(DotResolvedPlacement {
                    point,
                    x_set: true,
                    y_set: true,
                    pin: pinned.then_some(PinConstraint::Fixed(point.x, point.y)),
                })
            }
            DotPlacementExpr::Axes { x, y } => {
                let mut point = Point2::new(0.0, 0.0);
                let mut x_pin = None;
                let mut y_pin = None;

                if let Some(axis) = x {
                    point.x = axis.coord.point_value();
                    if axis.pinned {
                        x_pin = Some(axis.coord.x_pin());
                    }
                }

                if let Some(axis) = y {
                    point.y = axis.coord.point_value();
                    if axis.pinned {
                        y_pin = Some(axis.coord.y_pin());
                    }
                }

                let pin = combine_axis_pins(x_pin, y_pin);
                Ok(DotResolvedPlacement {
                    point,
                    x_set: x.is_some(),
                    y_set: y.is_some(),
                    pin,
                })
            }
        }
    }

    fn resolve_ref(&self, target: DotPlacementRef) -> Result<DotPlacementTarget, String> {
        match target {
            DotPlacementRef::Node(index) => self
                .node_id_map
                .get(&index)
                .copied()
                .map(DotPlacementTarget::Node)
                .ok_or_else(|| {
                    format!("DOT pos references node id {index}, but it does not exist")
                }),
            DotPlacementRef::Edge(index) => self
                .edge_id_map
                .get(&index)
                .copied()
                .map(DotPlacementTarget::Edge)
                .ok_or_else(|| {
                    format!("DOT pos references edge id {index}, but it does not exist")
                }),
        }
    }
}

impl DotPlacementExpr {
    fn parse(input: &str) -> Result<Option<Self>, String> {
        let input = strip_dot_quotes(input);
        if input.is_empty() {
            return Ok(None);
        }

        if let Some(expr) = Self::parse_axes(input)? {
            return Ok(Some(expr));
        }

        let (input, pinned) = input
            .strip_suffix('!')
            .map(|input| (input.trim(), true))
            .unwrap_or((input, false));

        if let Some(expr) = Self::parse_ref(input, pinned) {
            return Ok(Some(expr));
        }

        if let Some((x, y)) = parse_dot_point(input) {
            return Ok(Some(DotPlacementExpr::Point {
                point: Point2::new(x, y),
                pinned,
            }));
        }

        Ok(None)
    }

    fn parse_axes(input: &str) -> Result<Option<Self>, String> {
        if !input
            .split(',')
            .any(|part| part.trim().starts_with("x:") || part.trim().starts_with("y:"))
        {
            return Ok(None);
        }

        let mut x = None;
        let mut y = None;

        for part in input.split(',') {
            let part = part.trim();
            if part.is_empty() {
                return Err("empty axis entry in DOT pos".into());
            }

            let (axis, value) = part.split_once(':').ok_or_else(|| {
                format!("axis entry {part:?} must be written as x:<coord> or y:<coord>")
            })?;
            let target = match axis.trim() {
                "x" => &mut x,
                "y" => &mut y,
                other => return Err(format!("unknown DOT pos axis {other:?}")),
            };

            if target.is_some() {
                return Err(format!("duplicate DOT pos axis {axis:?}"));
            }

            *target = Some(parse_dot_axis_placement(value.trim())?);
        }

        Ok(Some(DotPlacementExpr::Axes { x, y }))
    }

    fn parse_ref(input: &str, pinned: bool) -> Option<Self> {
        let rest = input.strip_prefix("ref(")?;
        let close = rest.find(')')?;
        let target = Self::parse_ref_target(&rest[..close])?;
        let offset = parse_ref_offset(rest[close + 1..].trim())?;
        Some(DotPlacementExpr::Ref {
            target,
            offset,
            pinned,
        })
    }

    fn parse_ref_target(input: &str) -> Option<DotPlacementRef> {
        let (kind, index) = input.split_once(':')?;
        let index = index.trim().parse::<usize>().ok()?;
        match kind.trim() {
            "node" => Some(DotPlacementRef::Node(index)),
            "edge" => Some(DotPlacementRef::Edge(index)),
            _ => None,
        }
    }
}

impl DotAxisCoord {
    fn point_value(&self) -> f64 {
        match self {
            DotAxisCoord::Number(value) => *value,
            DotAxisCoord::Group(_) => 0.0,
        }
    }

    fn x_pin(&self) -> PinConstraint {
        match self {
            DotAxisCoord::Number(value) => PinConstraint::FixX(*value),
            DotAxisCoord::Group(group) => PinConstraint::LinkX(group.clone()),
        }
    }

    fn y_pin(&self) -> PinConstraint {
        match self {
            DotAxisCoord::Number(value) => PinConstraint::FixY(*value),
            DotAxisCoord::Group(group) => PinConstraint::LinkY(group.clone()),
        }
    }
}

fn parse_dot_axis_placement(input: &str) -> Result<DotAxisPlacement, String> {
    let (input, pinned) = input
        .strip_suffix('!')
        .map(|input| (input.trim(), true))
        .unwrap_or((input, false));

    if let Some(group) = parse_dot_link_group(input) {
        if !pinned {
            return Err(format!(
                "group coordinate {input:?} in DOT pos axis syntax must be pinned with !"
            ));
        }
        return Ok(DotAxisPlacement {
            coord: DotAxisCoord::Group(group),
            pinned,
        });
    }

    let value = input
        .parse::<f64>()
        .map_err(|_| format!("DOT pos axis coordinate {input:?} is not a number or group"))?;
    Ok(DotAxisPlacement {
        coord: DotAxisCoord::Number(value),
        pinned,
    })
}

fn parse_dot_link_group(input: &str) -> Option<String> {
    if let Some(group) = input.strip_prefix('@') {
        Some(group.to_string())
    } else if let Some(group) = input.strip_prefix("+@") {
        Some(format!("+{group}"))
    } else {
        input.strip_prefix("-@").map(|group| format!("-{group}"))
    }
}

fn combine_axis_pins(
    x_pin: Option<PinConstraint>,
    y_pin: Option<PinConstraint>,
) -> Option<PinConstraint> {
    match (x_pin, y_pin) {
        (Some(x_pin), Some(y_pin)) => {
            Some(PinConstraint::Combined(Box::new(x_pin), Box::new(y_pin)))
        }
        (Some(pin), None) | (None, Some(pin)) => Some(pin),
        (None, None) => None,
    }
}

fn parse_ref_offset(input: &str) -> Option<Vector2<f64>> {
    if input.is_empty() {
        return Some(Vector2::zero());
    }

    let input = input.strip_prefix('+').unwrap_or(input);
    let (dx, dy) = parse_dot_point(input)?;
    Some(Vector2::new(dx, dy))
}

fn strip_dot_quotes(input: &str) -> &str {
    input.trim().trim_matches('"').trim()
}

fn parse_dot_point(input: &str) -> Option<(f64, f64)> {
    let cleaned = strip_dot_quotes(input)
        .trim_matches(|c| c == '(' || c == ')')
        .trim();
    let parts: Vec<&str> = cleaned
        .split([',', ' '])
        .filter(|part| !part.is_empty())
        .collect();

    if parts.len() != 2 {
        return None;
    }

    Some((parts[0].parse().ok()?, parts[1].parse().ok()?))
}

fn apply_dot_placement_statements(
    statements: &mut BTreeMap<String, String>,
    placement: &DotResolvedPlacement,
) {
    statements.insert(
        "pos".to_string(),
        format!("{},{}", placement.point.x, placement.point.y),
    );
    statements.insert("pos-x-set".to_string(), placement.x_set.to_string());
    statements.insert("pos-y-set".to_string(), placement.y_set.to_string());
    statements.insert(
        "pos-mode".to_string(),
        if placement.pin.is_some() {
            "pin"
        } else {
            "start"
        }
        .to_string(),
    );
}

fn initial_point_constraint(
    index: usize,
    statements: &BTreeMap<String, String>,
    placement: Option<&DotResolvedPlacement>,
    group_map: &mut HashMap<String, usize>,
) -> (Point2<f64>, PointConstraint) {
    let explicit_pin =
        dot_statement_value(statements, "pin").and_then(|value| PinConstraint::parse(value));
    let pin = explicit_pin.or_else(|| placement.and_then(|placement| placement.pin.clone()));

    if let Some(pin) = pin {
        let (mut point, constraints) = pin.point_constraint(index, group_map);
        if let Some(placement) = placement {
            if matches!(constraints.x, Constraint::Free) && placement.x_set {
                point.x = placement.point.x;
            }
            if matches!(constraints.y, Constraint::Free) && placement.y_set {
                point.y = placement.point.y;
            }
        } else if let Some((x, y)) = TypstNode::parse_position(statements, "pos") {
            if matches!(constraints.x, Constraint::Free) {
                point.x = x;
            }
            if matches!(constraints.y, Constraint::Free) {
                point.y = y;
            }
        }
        (point, constraints)
    } else if let Some(placement) = placement {
        (placement.point, PointConstraint::default())
    } else if let Some((x, y)) = TypstNode::parse_position(statements, "pos") {
        (Point2::new(x, y), PointConstraint::default())
    } else {
        (Point2::new(0., 0.), PointConstraint::default())
    }
}

impl TypstGraph {
    pub fn from_dot(dot: DotGraph, figment: &Figment) -> Self {
        let config = LayoutConfig::from_figment(figment);
        Self::from_dot_with_layout_config(dot, config)
    }

    pub(crate) fn from_dot_with_layout_config(dot: DotGraph, config: LayoutConfig) -> Self {
        let mut placement_context = DotPlacementContext::new(&dot);
        let node_placements = dot
            .graph
            .new_nodevec(|nid, _, _| placement_context.resolve_node(nid));
        let edge_placements = dot
            .graph
            .new_edgevec(|_, eid, _| placement_context.resolve_edge(eid));

        let mut group_map = HashMap::new();
        let edge_pin_constrains: EdgeVec<(Point2<f64>, PointConstraint)> =
            dot.graph.new_edgevec(|e, eid, _| {
                initial_point_constraint(
                    eid.0,
                    &e.statements,
                    edge_placements[eid].as_ref(),
                    &mut group_map,
                )
            });

        let mut group_map = HashMap::new();

        let node_pin_constrains: NodeVec<(Point2<f64>, PointConstraint)> =
            dot.graph.new_nodevec(|nid, _, n| {
                initial_point_constraint(
                    nid.0,
                    &n.statements,
                    node_placements[nid].as_ref(),
                    &mut group_map,
                )
            });

        let graph = dot.graph.map(
            |inv, nid, mut data| {
                if let Some(placement) = &node_placements[nid] {
                    apply_dot_placement_statements(&mut data.statements, placement);
                }
                TypstNode::parse(inv, nid, data, &node_pin_constrains)
            },
            |inv, store, p, eid, data| {
                let mut data = data;
                if let Some(placement) = &edge_placements[eid] {
                    apply_dot_placement_statements(&mut data.data.statements, placement);
                }
                TypstEdge::parse(inv, store, p, eid, data, &edge_pin_constrains)
            },
            TypstHedge::parse,
        );

        let global_eval: Option<String> = dot.global_data.statements.get("eval").cloned();

        Self {
            graph,
            global_eval,
            name: dot.global_data.name,
            data: dot.global_data.payload,
            global_statements: dot
                .global_data
                .statements
                .into_iter()
                .map(|(k, v)| {
                    let clean_key = k.trim().trim_matches('"').to_string();
                    let clean_value = v.trim().trim_matches('"').to_string();
                    (clean_key, clean_value)
                })
                .collect(),
            default_edge_statements: dot
                .global_data
                .edge_statements
                .into_iter()
                .map(|(k, v)| {
                    let clean_key = k.trim().trim_matches('"').to_string();
                    let clean_value = v.trim().trim_matches('"').to_string();
                    (clean_key, clean_value)
                })
                .collect(),
            default_node_statements: dot
                .global_data
                .node_statements
                .into_iter()
                .map(|(k, v)| {
                    let clean_key = k.trim().trim_matches('"').to_string();
                    let clean_value = v.trim().trim_matches('"').to_string();
                    (clean_key, clean_value)
                })
                .collect(),
            layout_config: config,
        }
    }
}

impl From<(DotGraph, Figment)> for TypstGraph {
    fn from(value: (DotGraph, Figment)) -> Self {
        TypstGraph::from_dot(value.0, &value.1)
    }
}

#[derive(Clone, Copy)]
pub struct TreeInitCfg {
    pub dy: f64, // ≈ 1.2 * L
    pub dx: f64, // ≈ 0.9 * L
}

#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize,
)]
#[serde(rename_all = "kebab-case")]
enum LayoutAlgo {
    Anneal,
    Dot,
    Force,
    #[serde(alias = "railroad")]
    StableLayered,
    Tree,
}

fn default_layout_algo() -> LayoutAlgo {
    LayoutAlgo::Force
}

#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize,
)]
#[serde(rename_all = "lowercase")]
enum LayoutNodeMode {
    Fixed,
    Layout,
}

fn default_layout_node_mode() -> LayoutNodeMode {
    LayoutNodeMode::Layout
}

impl LayoutNodeMode {
    fn nodes_are_fixed(self) -> bool {
        matches!(self, LayoutNodeMode::Fixed)
    }
}

#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize,
)]
#[serde(rename_all = "kebab-case")]
enum LabelLayout {
    DanglingTangent,
    FixedLength,
    Normal,
}

fn default_label_layout() -> LabelLayout {
    LabelLayout::Normal
}

#[derive(
    Debug, Clone, Serialize, Deserialize, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize,
)]
#[serde(rename_all = "kebab-case")]
struct LayoutConfig {
    #[serde(default = "default_viewport_w", deserialize_with = "deserialize_f64")]
    viewport_w: f64,
    #[serde(default = "default_viewport_h", deserialize_with = "deserialize_f64")]
    viewport_h: f64,
    #[serde(default = "default_tree_dy", deserialize_with = "deserialize_f64")]
    tree_dy: f64,
    #[serde(default = "default_tree_dx", deserialize_with = "deserialize_f64")]
    tree_dx: f64,
    #[serde(default = "default_step", deserialize_with = "deserialize_f64")]
    step: f64,
    #[serde(default = "default_temp", deserialize_with = "deserialize_f64")]
    temp: f64,
    #[serde(default = "default_seed", deserialize_with = "deserialize_u64")]
    seed: u64,
    #[serde(default = "default_delta", deserialize_with = "deserialize_f64")]
    delta: f64,
    #[serde(
        default = "default_directional_force",
        deserialize_with = "deserialize_f64"
    )]
    directional_force: f64,
    #[serde(default = "default_z_spring", deserialize_with = "deserialize_f64")]
    z_spring: f64,
    #[serde(
        default = "default_z_spring_growth",
        deserialize_with = "deserialize_f64"
    )]
    z_spring_growth: f64,
    #[serde(
        default = "default_label_steps",
        deserialize_with = "deserialize_usize"
    )]
    label_steps: usize,
    #[serde(default = "default_label_layout")]
    label_layout: LabelLayout,
    #[serde(default = "default_label_step", deserialize_with = "deserialize_f64")]
    label_step: f64,
    #[serde(
        default = "default_label_length_scale",
        deserialize_with = "deserialize_f64"
    )]
    label_length_scale: f64,
    #[serde(default = "default_label_charge", deserialize_with = "deserialize_f64")]
    label_charge: f64,
    #[serde(default = "default_label_spring", deserialize_with = "deserialize_f64")]
    label_spring: f64,
    #[serde(
        default = "default_label_early_tol",
        deserialize_with = "deserialize_f64"
    )]
    label_early_tol: f64,
    #[serde(
        default = "default_label_max_delta_scale",
        deserialize_with = "deserialize_f64"
    )]
    label_max_delta_scale: f64,
    #[serde(
        default = "default_incremental_energy",
        deserialize_with = "deserialize_bool"
    )]
    incremental_energy: bool,
    #[serde(default = "default_layout_algo")]
    layout_algo: LayoutAlgo,
    #[serde(default = "default_layout_node_mode")]
    layout_nodes: LayoutNodeMode,
    #[serde(default, deserialize_with = "deserialize_usize_vec")]
    layout_roots: Vec<usize>,
    #[serde(default, deserialize_with = "deserialize_string_vec")]
    rank_same: Vec<String>,
    #[serde(
        default = "default_route_edge_weight",
        deserialize_with = "deserialize_f64"
    )]
    route_edge_weight: f64,
    #[serde(
        default = "default_route_exit_weight",
        deserialize_with = "deserialize_f64"
    )]
    route_exit_weight: f64,
    #[serde(
        default = "default_route_label_width_scale",
        deserialize_with = "deserialize_f64"
    )]
    route_label_width_scale: f64,
    #[serde(
        default = "default_route_label_width_cap",
        deserialize_with = "deserialize_f64"
    )]
    route_label_width_cap: f64,
    #[serde(default, flatten)]
    spring: SpringConfig,
    #[serde(default, flatten)]
    schedule: ScheduleConfig,
}

impl Default for LayoutConfig {
    fn default() -> Self {
        Self {
            viewport_w: default_viewport_w(),
            viewport_h: default_viewport_h(),
            tree_dy: default_tree_dy(),
            tree_dx: default_tree_dx(),
            step: default_step(),
            temp: default_temp(),
            seed: default_seed(),
            delta: default_delta(),
            directional_force: default_directional_force(),
            z_spring: default_z_spring(),
            z_spring_growth: default_z_spring_growth(),
            label_steps: default_label_steps(),
            label_layout: default_label_layout(),
            label_step: default_label_step(),
            label_length_scale: default_label_length_scale(),
            label_charge: default_label_charge(),
            label_spring: default_label_spring(),
            label_early_tol: default_label_early_tol(),
            label_max_delta_scale: default_label_max_delta_scale(),
            incremental_energy: default_incremental_energy(),
            layout_algo: default_layout_algo(),
            layout_nodes: default_layout_node_mode(),
            layout_roots: Vec::new(),
            rank_same: Vec::new(),
            route_edge_weight: default_route_edge_weight(),
            route_exit_weight: default_route_exit_weight(),
            route_label_width_scale: default_route_label_width_scale(),
            route_label_width_cap: default_route_label_width_cap(),
            spring: SpringConfig::default(),
            schedule: ScheduleConfig::default(),
        }
    }
}

impl LayoutConfig {
    fn from_figment(figment: &Figment) -> Self {
        figment.extract::<Self>().unwrap_or_default()
    }

    fn is_layout_statement_key(key: &str) -> bool {
        let normalized = key.replace('_', "-");
        matches!(
            normalized.as_str(),
            "accept-floor"
                | "beta"
                | "cool"
                | "crossing-penalty"
                | "delta"
                | "directional-force"
                | "early-tol"
                | "epochs"
                | "eps"
                | "g-center"
                | "gamma-dangling"
                | "gamma-ee"
                | "gamma-ev"
                | "incremental-energy"
                | "k-spring"
                | "label-charge"
                | "label-early-tol"
                | "label-layout"
                | "label-length-scale"
                | "label-max-delta-scale"
                | "label-spring"
                | "label-step"
                | "label-steps"
                | "layout-algo"
                | "layout-nodes"
                | "layout-roots"
                | "length-scale"
                | "rank-same"
                | "route-edge-weight"
                | "route-exit-weight"
                | "route-label-width-cap"
                | "route-label-width-scale"
                | "seed"
                | "sink-route-exit"
                | "source-route-exit"
                | "step"
                | "step-shrink"
                | "steps"
                | "temp"
                | "tree-dx"
                | "tree-dy"
                | "viewport-h"
                | "viewport-w"
                | "z-spring"
                | "z-spring-growth"
        )
    }
}

pub(crate) fn default_figment() -> Figment {
    Figment::from(Serialized::from(
        BTreeMap::<String, String>::new(),
        Profile::Default,
    ))
}

fn default_viewport_w() -> f64 {
    10.0
}

fn default_viewport_h() -> f64 {
    10.0
}

fn default_tree_dy() -> f64 {
    1.2
}

fn default_tree_dx() -> f64 {
    0.9
}

fn default_step() -> f64 {
    0.81
}

fn default_temp() -> f64 {
    0.3
}

fn default_seed() -> u64 {
    2
}

fn default_delta() -> f64 {
    0.4
}

fn default_directional_force() -> f64 {
    5.0
}

fn default_z_spring() -> f64 {
    2.0
}

fn default_z_spring_growth() -> f64 {
    1.0
}

fn default_label_steps() -> usize {
    20
}

fn default_label_step() -> f64 {
    0.15
}

fn default_label_length_scale() -> f64 {
    0.6
}

fn default_label_charge() -> f64 {
    3.0
}

fn default_label_spring() -> f64 {
    23.0
}

fn default_label_early_tol() -> f64 {
    1e-3
}

fn default_label_max_delta_scale() -> f64 {
    0.5
}

fn default_incremental_energy() -> bool {
    true
}

fn default_route_edge_weight() -> f64 {
    0.15
}

fn default_route_exit_weight() -> f64 {
    4.0
}

fn default_route_label_width_scale() -> f64 {
    1.0
}

fn default_route_label_width_cap() -> f64 {
    2.0
}

fn default_crossing_penalty() -> f64 {
    30.0
}

#[derive(
    Debug, Clone, Serialize, Deserialize, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize,
)]
#[serde(rename_all = "kebab-case")]
struct SpringConfig {
    #[serde(default = "default_length_scale", deserialize_with = "deserialize_f64")]
    length_scale: f64,
    #[serde(default = "default_k_spring", deserialize_with = "deserialize_f64")]
    k_spring: f64,
    #[serde(default = "default_beta", deserialize_with = "deserialize_f64")]
    beta: f64,
    #[serde(
        default = "default_gamma_dangling",
        deserialize_with = "deserialize_f64"
    )]
    gamma_dangling: f64,
    #[serde(default = "default_gamma_ev", deserialize_with = "deserialize_f64")]
    gamma_ev: f64,
    #[serde(default = "default_gamma_ee", deserialize_with = "deserialize_f64")]
    gamma_ee: f64,
    #[serde(default = "default_g_center", deserialize_with = "deserialize_f64")]
    g_center: f64,
    #[serde(
        default = "default_crossing_penalty",
        deserialize_with = "deserialize_f64"
    )]
    crossing_penalty: f64,
    #[serde(default = "default_eps", deserialize_with = "deserialize_f64")]
    eps: f64,
}

impl Default for SpringConfig {
    fn default() -> Self {
        Self {
            length_scale: default_length_scale(),
            k_spring: default_k_spring(),
            beta: default_beta(),
            gamma_dangling: default_gamma_dangling(),
            gamma_ev: default_gamma_ev(),
            gamma_ee: default_gamma_ee(),
            g_center: default_g_center(),
            crossing_penalty: default_crossing_penalty(),
            eps: default_eps(),
        }
    }
}

impl From<&SpringConfig> for ParamTuning {
    fn from(cfg: &SpringConfig) -> Self {
        ParamTuning {
            length_scale: cfg.length_scale,
            k_spring: cfg.k_spring,
            beta: cfg.beta,
            gamma_dangling: cfg.gamma_dangling,
            gamma_ev: cfg.gamma_ev,
            gamma_ee: cfg.gamma_ee,
            g_center: cfg.g_center,
            crossing_penalty: cfg.crossing_penalty,
            eps: cfg.eps,
        }
    }
}

#[derive(
    Debug, Clone, Serialize, Deserialize, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize,
)]
#[serde(rename_all = "kebab-case")]
struct ScheduleConfig {
    #[serde(default = "default_steps", deserialize_with = "deserialize_usize")]
    steps: usize,
    #[serde(default = "default_epochs", deserialize_with = "deserialize_usize")]
    epochs: usize,
    #[serde(default = "default_cool", deserialize_with = "deserialize_f64")]
    cool: f64,
    #[serde(default = "default_accept_floor", deserialize_with = "deserialize_f64")]
    accept_floor: f64,
    #[serde(default = "default_step_shrink", deserialize_with = "deserialize_f64")]
    step_shrink: f64,
    #[serde(default = "default_early_tol", deserialize_with = "deserialize_f64")]
    early_tol: f64,
}

impl Default for ScheduleConfig {
    fn default() -> Self {
        Self {
            steps: default_steps(),
            epochs: default_epochs(),
            cool: default_cool(),
            accept_floor: default_accept_floor(),
            step_shrink: default_step_shrink(),
            early_tol: default_early_tol(),
        }
    }
}

impl From<&ScheduleConfig> for GeoSchedule {
    fn from(cfg: &ScheduleConfig) -> Self {
        GeoSchedule::new(
            cfg.steps,
            cfg.epochs,
            cfg.cool,
            cfg.accept_floor,
            cfg.step_shrink,
            cfg.early_tol,
        )
    }
}

fn default_length_scale() -> f64 {
    0.5
}

fn default_k_spring() -> f64 {
    11.0
}

fn default_beta() -> f64 {
    50.0
}

fn default_gamma_dangling() -> f64 {
    5.0
}

fn default_gamma_ev() -> f64 {
    0.01
}

fn default_gamma_ee() -> f64 {
    0.10
}

fn default_g_center() -> f64 {
    0.005
}

fn default_eps() -> f64 {
    1e-4
}

fn default_steps() -> usize {
    30
}

fn default_epochs() -> usize {
    30
}

fn default_cool() -> f64 {
    0.85
}

fn default_accept_floor() -> f64 {
    0.15
}

fn default_step_shrink() -> f64 {
    0.21
}

fn default_early_tol() -> f64 {
    1e-6
}

fn deserialize_f64<'de, D>(deserializer: D) -> Result<f64, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    struct F64Visitor;

    impl<'de> serde::de::Visitor<'de> for F64Visitor {
        type Value = f64;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("a number or a numeric string")
        }

        fn visit_f64<E>(self, value: f64) -> Result<f64, E> {
            Ok(value)
        }

        fn visit_u64<E>(self, value: u64) -> Result<f64, E> {
            Ok(value as f64)
        }

        fn visit_i64<E>(self, value: i64) -> Result<f64, E> {
            Ok(value as f64)
        }

        fn visit_str<E>(self, value: &str) -> Result<f64, E>
        where
            E: serde::de::Error,
        {
            value.parse().map_err(E::custom)
        }

        fn visit_string<E>(self, value: String) -> Result<f64, E>
        where
            E: serde::de::Error,
        {
            self.visit_str(&value)
        }
    }

    deserializer.deserialize_any(F64Visitor)
}

fn deserialize_u64<'de, D>(deserializer: D) -> Result<u64, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    struct U64Visitor;

    impl<'de> serde::de::Visitor<'de> for U64Visitor {
        type Value = u64;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("a non-negative integer or numeric string")
        }

        fn visit_u64<E>(self, value: u64) -> Result<u64, E> {
            Ok(value)
        }

        fn visit_i64<E>(self, value: i64) -> Result<u64, E>
        where
            E: serde::de::Error,
        {
            if value < 0 {
                return Err(E::custom("expected non-negative integer"));
            }
            Ok(value as u64)
        }

        fn visit_f64<E>(self, value: f64) -> Result<u64, E>
        where
            E: serde::de::Error,
        {
            if value.is_sign_negative() {
                return Err(E::custom("expected non-negative integer"));
            }
            Ok(value as u64)
        }

        fn visit_str<E>(self, value: &str) -> Result<u64, E>
        where
            E: serde::de::Error,
        {
            value.parse().map_err(E::custom)
        }

        fn visit_string<E>(self, value: String) -> Result<u64, E>
        where
            E: serde::de::Error,
        {
            self.visit_str(&value)
        }
    }

    deserializer.deserialize_any(U64Visitor)
}

fn deserialize_usize<'de, D>(deserializer: D) -> Result<usize, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    struct UsizeVisitor;

    impl<'de> serde::de::Visitor<'de> for UsizeVisitor {
        type Value = usize;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("a non-negative integer or numeric string")
        }

        fn visit_u64<E>(self, value: u64) -> Result<usize, E> {
            Ok(value as usize)
        }

        fn visit_i64<E>(self, value: i64) -> Result<usize, E>
        where
            E: serde::de::Error,
        {
            if value < 0 {
                return Err(E::custom("expected non-negative integer"));
            }
            Ok(value as usize)
        }

        fn visit_f64<E>(self, value: f64) -> Result<usize, E>
        where
            E: serde::de::Error,
        {
            if value.is_sign_negative() {
                return Err(E::custom("expected non-negative integer"));
            }
            Ok(value as usize)
        }

        fn visit_str<E>(self, value: &str) -> Result<usize, E>
        where
            E: serde::de::Error,
        {
            value.parse().map_err(E::custom)
        }

        fn visit_string<E>(self, value: String) -> Result<usize, E>
        where
            E: serde::de::Error,
        {
            self.visit_str(&value)
        }
    }

    deserializer.deserialize_any(UsizeVisitor)
}

fn deserialize_bool<'de, D>(deserializer: D) -> Result<bool, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    struct BoolVisitor;

    impl<'de> serde::de::Visitor<'de> for BoolVisitor {
        type Value = bool;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("a boolean or boolean string")
        }

        fn visit_bool<E>(self, value: bool) -> Result<bool, E> {
            Ok(value)
        }

        fn visit_str<E>(self, value: &str) -> Result<bool, E>
        where
            E: serde::de::Error,
        {
            match value.trim().to_ascii_lowercase().as_str() {
                "true" | "on" | "yes" | "1" => Ok(true),
                "false" | "off" | "no" | "0" => Ok(false),
                other => Err(E::custom(format!("invalid boolean string {other:?}"))),
            }
        }

        fn visit_string<E>(self, value: String) -> Result<bool, E>
        where
            E: serde::de::Error,
        {
            self.visit_str(&value)
        }
    }

    deserializer.deserialize_any(BoolVisitor)
}

fn deserialize_usize_vec<'de, D>(deserializer: D) -> Result<Vec<usize>, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    struct UsizeVecVisitor;

    impl<'de> serde::de::Visitor<'de> for UsizeVecVisitor {
        type Value = Vec<usize>;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("an array of node indices or a comma-separated index string")
        }

        fn visit_seq<A>(self, mut seq: A) -> Result<Vec<usize>, A::Error>
        where
            A: serde::de::SeqAccess<'de>,
        {
            let mut values = Vec::new();
            while let Some(value) = seq.next_element::<usize>()? {
                values.push(value);
            }
            Ok(values)
        }

        fn visit_u64<E>(self, value: u64) -> Result<Vec<usize>, E>
        where
            E: serde::de::Error,
        {
            Ok(vec![usize::try_from(value).map_err(E::custom)?])
        }

        fn visit_i64<E>(self, value: i64) -> Result<Vec<usize>, E>
        where
            E: serde::de::Error,
        {
            if value < 0 {
                return Err(E::custom("layout root indices must be non-negative"));
            }
            Ok(vec![usize::try_from(value).map_err(E::custom)?])
        }

        fn visit_str<E>(self, value: &str) -> Result<Vec<usize>, E>
        where
            E: serde::de::Error,
        {
            if value.trim().is_empty() {
                return Ok(Vec::new());
            }
            value
                .split([',', ' '])
                .filter(|part| !part.trim().is_empty())
                .map(|part| part.trim().parse().map_err(E::custom))
                .collect()
        }

        fn visit_string<E>(self, value: String) -> Result<Vec<usize>, E>
        where
            E: serde::de::Error,
        {
            self.visit_str(&value)
        }
    }

    deserializer.deserialize_any(UsizeVecVisitor)
}

fn deserialize_string_vec<'de, D>(deserializer: D) -> Result<Vec<String>, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    struct StringVecVisitor;

    impl<'de> serde::de::Visitor<'de> for StringVecVisitor {
        type Value = Vec<String>;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("an array of strings or a comma-separated string")
        }

        fn visit_seq<A>(self, mut seq: A) -> Result<Vec<String>, A::Error>
        where
            A: serde::de::SeqAccess<'de>,
        {
            let mut values = Vec::new();
            while let Some(value) = seq.next_element::<String>()? {
                values.push(value);
            }
            Ok(values)
        }

        fn visit_str<E>(self, value: &str) -> Result<Vec<String>, E>
        where
            E: serde::de::Error,
        {
            if value.trim().is_empty() {
                return Ok(Vec::new());
            }
            Ok(value
                .split([',', ' '])
                .filter_map(|part| {
                    let part = part.trim();
                    (!part.is_empty()).then(|| part.to_string())
                })
                .collect())
        }

        fn visit_string<E>(self, value: String) -> Result<Vec<String>, E>
        where
            E: serde::de::Error,
        {
            self.visit_str(&value)
        }
    }

    deserializer.deserialize_any(StringVecVisitor)
}

impl TypstGraph {
    pub fn layout(&mut self) {
        self.layout_with_subgraph(None)
            .expect("full graph layout should not fail");
    }

    pub fn layout_with_subgraph(&mut self, subgraph: Option<&SuBitGraph>) -> Result<(), String> {
        let spring_params = ParamTuning::from(&self.layout_config.spring);
        self.clear_hedge_route_points();

        let (tree_cfg, energy) = self.tree_init_cfg(&spring_params);
        let fixed_nodes = self.layout_config.layout_nodes.nodes_are_fixed();
        if let Some(subgraph) = subgraph {
            let (mut vertex_points, mut edge_points) = match self.layout_config.layout_algo {
                LayoutAlgo::Tree | LayoutAlgo::Dot | LayoutAlgo::StableLayered => {
                    if fixed_nodes {
                        let (_, selected_edges) = self.selected_layout_items(subgraph);
                        self.fixed_node_partial_layout_positions(
                            tree_cfg,
                            self.layout_config.layout_algo,
                            Some(&selected_edges),
                        )
                    } else {
                        self.partial_layout_positions(
                            tree_cfg,
                            self.layout_config.layout_algo,
                            subgraph,
                        )
                    }
                }
                LayoutAlgo::Anneal | LayoutAlgo::Force => {
                    self.partial_optimized_positions(tree_cfg, subgraph, &energy)
                }
            };
            self.apply_layout_constraints(&mut vertex_points, &mut edge_points);
            self.update_positions(vertex_points, edge_points);
            if matches!(
                self.layout_config.layout_algo,
                LayoutAlgo::Dot | LayoutAlgo::StableLayered
            ) {
                self.clear_edge_label_positions();
            } else {
                self.layout_edge_labels(energy.spring_length);
            }
            return Ok(());
        }

        if matches!(
            self.layout_config.layout_algo,
            LayoutAlgo::Tree | LayoutAlgo::Dot | LayoutAlgo::StableLayered
        ) {
            let (mut vertex_points, mut edge_points) = if fixed_nodes {
                self.fixed_node_direct_layout_positions(tree_cfg, self.layout_config.layout_algo)
            } else {
                self.direct_layout_positions(tree_cfg, self.layout_config.layout_algo)
            };
            self.apply_layout_constraints(&mut vertex_points, &mut edge_points);
            self.update_positions(vertex_points, edge_points);
            if matches!(
                self.layout_config.layout_algo,
                LayoutAlgo::Dot | LayoutAlgo::StableLayered
            ) {
                self.clear_edge_label_positions();
            } else {
                self.layout_edge_labels(energy.spring_length);
            }
            return Ok(());
        }

        let (mut vertex_points, mut edge_points) = if fixed_nodes {
            let full = self.full_filter();
            self.partial_optimized_positions(tree_cfg, &full, &energy)
        } else {
            let (pos_n, pos_e) = self.new_positions(tree_cfg);
            self.optimized_positions(pos_n, pos_e, &energy)
        };

        self.apply_layout_constraints(&mut vertex_points, &mut edge_points);
        self.update_positions(vertex_points, edge_points);
        self.layout_edge_labels(energy.spring_length);
        Ok(())
    }

    pub fn layout_energy_state(
        &self,
    ) -> (
        LayoutState<'_, TypstEdge, TypstNode, TypstHedge, DefaultNodeStore<TypstNode>>,
        SpringChargeEnergy,
    ) {
        let spring_params = ParamTuning::from(&self.layout_config.spring);
        let (tree_cfg, energy) = self.tree_init_cfg(&spring_params);
        let (pos_n, pos_e) = self.new_positions(tree_cfg);
        let state = self.graph.new_layout_state(
            pos_n,
            pos_e,
            self.layout_config.delta,
            self.layout_config.directional_force,
            self.layout_config.incremental_energy,
        );

        (state, energy)
    }

    fn tree_init_cfg(&self, tune: &ParamTuning) -> (TreeInitCfg, SpringChargeEnergy) {
        let energycfg = SpringChargeEnergy::from_graph(
            self.n_nodes(),
            self.layout_config.viewport_w,
            self.layout_config.viewport_h,
            *tune,
        );

        let l = energycfg.spring_length;
        (
            TreeInitCfg {
                dy: self.layout_config.tree_dy * l,
                dx: self.layout_config.tree_dx * l,
            },
            energycfg,
        )
    }

    fn apply_grouped_constraints(
        &self,
        pos_v: &mut NodeVec<Point2<f64>>,
        pos_e: &mut EdgeVec<Point2<f64>>,
    ) {
        let node_len = pos_v.len().0;
        for i in 0..node_len {
            let idx = NodeIndex(i);
            let constraints = &self[idx].constraints;
            if let Constraint::Grouped(reference, _) = constraints.x {
                if reference < node_len && reference != i {
                    pos_v[idx].x = pos_v[NodeIndex(reference)].x;
                }
            }
            if let Constraint::Grouped(reference, _) = constraints.y {
                if reference < node_len && reference != i {
                    pos_v[idx].y = pos_v[NodeIndex(reference)].y;
                }
            }
        }

        let edge_len = pos_e.len().0;
        for i in 0..edge_len {
            let idx = EdgeIndex(i);
            let constraints = &self[idx].constraints;
            if let Constraint::Grouped(reference, _) = constraints.x {
                if reference < edge_len && reference != i {
                    pos_e[idx].x = pos_e[EdgeIndex(reference)].x;
                }
            }
            if let Constraint::Grouped(reference, _) = constraints.y {
                if reference < edge_len && reference != i {
                    pos_e[idx].y = pos_e[EdgeIndex(reference)].y;
                }
            }
        }

        for i in 0..node_len {
            let idx = NodeIndex(i);
            Self::apply_directional_constraint(&self[idx].constraints.x, &mut pos_v[idx].x);
            Self::apply_directional_constraint(&self[idx].constraints.y, &mut pos_v[idx].y);
        }

        for i in 0..edge_len {
            let idx = EdgeIndex(i);
            Self::apply_directional_constraint(&self[idx].constraints.x, &mut pos_e[idx].x);
            Self::apply_directional_constraint(&self[idx].constraints.y, &mut pos_e[idx].y);
        }
    }

    fn apply_layout_constraints(
        &self,
        pos_v: &mut NodeVec<Point2<f64>>,
        pos_e: &mut EdgeVec<Point2<f64>>,
    ) {
        if self.layout_config.layout_nodes.nodes_are_fixed() {
            self.apply_edge_grouped_constraints(pos_e);
        } else {
            self.apply_grouped_constraints(pos_v, pos_e);
        }
    }

    fn apply_edge_grouped_constraints(&self, pos_e: &mut EdgeVec<Point2<f64>>) {
        let edge_len = pos_e.len().0;
        for i in 0..edge_len {
            let idx = EdgeIndex(i);
            let constraints = &self[idx].constraints;
            if let Constraint::Grouped(reference, _) = constraints.x {
                if reference < edge_len && reference != i {
                    pos_e[idx].x = pos_e[EdgeIndex(reference)].x;
                }
            }
            if let Constraint::Grouped(reference, _) = constraints.y {
                if reference < edge_len && reference != i {
                    pos_e[idx].y = pos_e[EdgeIndex(reference)].y;
                }
            }
        }

        for i in 0..edge_len {
            let idx = EdgeIndex(i);
            Self::apply_directional_constraint(&self[idx].constraints.x, &mut pos_e[idx].x);
            Self::apply_directional_constraint(&self[idx].constraints.y, &mut pos_e[idx].y);
        }
    }

    fn apply_directional_constraint(constraint: &Constraint, value: &mut f64) {
        match constraint {
            Constraint::Grouped(_, ShiftDirection::PositiveOnly) if *value <= 0.0 => {
                *value = value.abs().max(f64::EPSILON);
            }
            Constraint::Grouped(_, ShiftDirection::NegativeOnly) if *value >= 0.0 => {
                *value = -value.abs().max(f64::EPSILON);
            }
            _ => {}
        }
    }

    fn optimized_positions(
        &self,
        pos_n: NodeVec<Point2<f64>>,
        pos_e: EdgeVec<Point2<f64>>,
        energy: &SpringChargeEnergy,
    ) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        let spring_length = energy.spring_length;

        match self.layout_config.layout_algo {
            LayoutAlgo::Anneal => {
                let state = self.graph.new_layout_state(
                    pos_n,
                    pos_e,
                    self.layout_config.delta * spring_length,
                    self.layout_config.directional_force,
                    self.layout_config.incremental_energy,
                );
                let mut schedule = GeoSchedule::from(&self.layout_config.schedule);
                let (out, _stats) = anneal::<_, _, _, _, SmallRng>(
                    state,
                    SAConfig {
                        temp: self.layout_config.temp * spring_length.powi(2),
                        step: self.layout_config.step * spring_length,
                        seed: self.layout_config.seed,
                    },
                    &PinnedLayoutNeighbor,
                    energy,
                    &mut schedule,
                );
                (out.vertex_points, out.edge_points)
            }
            LayoutAlgo::Force => {
                let mut state = self.graph.new_layout_state(
                    pos_n,
                    pos_e,
                    self.layout_config.delta * spring_length,
                    self.layout_config.directional_force * spring_length,
                    self.layout_config.incremental_energy,
                );
                force_directed_layout(
                    &mut state,
                    energy,
                    ForceLayoutConfig {
                        steps: self.layout_config.schedule.steps,
                        epochs: self.layout_config.schedule.epochs,
                        step: self.layout_config.step,
                        cool: self.layout_config.schedule.cool,
                        max_delta: self.layout_config.delta * spring_length,
                        early_tol: self.layout_config.schedule.early_tol * spring_length,
                        seed: self.layout_config.seed,
                        z_spring: self.layout_config.z_spring,
                        z_spring_growth: self.layout_config.z_spring_growth,
                    },
                );
                (state.vertex_points, state.edge_points)
            }
            LayoutAlgo::Dot | LayoutAlgo::StableLayered | LayoutAlgo::Tree => unreachable!(),
        }
    }

    fn directional_target(value: f64, direction: ShiftDirection, fallback: f64) -> f64 {
        match direction {
            ShiftDirection::Any => value,
            ShiftDirection::PositiveOnly if value > 0.0 => value,
            ShiftDirection::PositiveOnly => fallback.abs().max(f64::EPSILON),
            ShiftDirection::NegativeOnly if value < 0.0 => value,
            ShiftDirection::NegativeOnly => -fallback.abs().max(f64::EPSILON),
        }
    }

    pub fn update_positions(&mut self, node: NodeVec<Point2<f64>>, edge: EdgeVec<Point2<f64>>) {
        node.into_iter().for_each(|(i, p)| {
            let p = p + self[i].shift.unwrap_or(Vector2::zero());
            self.graph[i].pos = p;
        });

        edge.into_iter().for_each(|(i, p)| {
            let p = p + self[i].shift.unwrap_or(Vector2::zero());

            let angle = {
                let (_, pair) = &self.graph[&i];

                match pair {
                    HedgePair::Split { source, sink, .. } | HedgePair::Paired { source, sink } => {
                        let so = self.node_id(*source);
                        let a = self[so].pos;
                        let si = self.node_id(*sink);
                        let b = self[si].pos;
                        tangent_angle_toward_c_side(a, b, p)
                    }
                    _ => Ok(Rad::zero()),
                }
            };

            if !self.graph[i].bend_explicit {
                self.graph[i].bend = angle;
            }
            self.graph[i].pos = p;
        });
    }

    fn layout_edge_labels(&mut self, spring_length: f64) {
        let cfg = &self.layout_config;
        if cfg.label_steps == 0 {
            return;
        }

        let label_length = cfg.label_length_scale * spring_length;
        let label_charge = cfg.label_charge * spring_length.powi(2);
        let label_spring = cfg.label_spring;
        let step = cfg.label_step;
        let max_delta = cfg.label_max_delta_scale * spring_length;
        let eps = 1e-4;

        let axes: EdgeVec<Vector2<f64>> =
            self.new_edgevec(|_e, idx, pair| self.edge_label_axis(idx, pair));
        let label_radii = self.edge_label_radii();
        let node_radii = self.node_layout_radii();

        let base_labels: EdgeVec<Point2<f64>> = self.new_edgevec(|_e, idx, _pair| {
            let edge_pos = self.graph[idx].pos;
            edge_pos + axes[idx] * label_length
        });
        let mut labels = base_labels.clone();

        for _ in 0..cfg.label_steps {
            let mut max_move: f64 = 0.0;

            for i in 0..labels.len().0 {
                let idx = EdgeIndex(i);
                let edge_pos = self.graph[idx].pos;
                let mut force = match cfg.label_layout {
                    LabelLayout::DanglingTangent | LabelLayout::Normal => {
                        let target = edge_pos + axes[idx] * label_length;
                        if label_spring != 0.0 {
                            (target - labels[idx]) * label_spring
                        } else {
                            Vector2::zero()
                        }
                    }
                    LabelLayout::FixedLength => Vector2::zero(),
                };

                if label_charge != 0.0 {
                    for n in 0..self.n_nodes() {
                        let ni = NodeIndex(n);
                        let np = self[ni].pos;
                        let d = labels[idx] - np;
                        let dist = d.magnitude();
                        if dist > 1e-9 {
                            let dir = d / dist;
                            let repel_dist =
                                (dist - label_radii[idx] - node_radii[ni]).max(0.0) + eps;
                            force += dir * (label_charge / repel_dist.powi(2));
                        }
                    }

                    for e in 0..labels.len().0 {
                        let ej = EdgeIndex(e);
                        if ej == idx {
                            continue;
                        }
                        let ep = self.graph[ej].pos;
                        let d = labels[idx] - ep;
                        let dist = d.magnitude();
                        if dist > 1e-9 {
                            let dir = d / dist;
                            let repel_dist = (dist - label_radii[idx]).max(0.0) + eps;
                            force += dir * (label_charge / repel_dist.powi(2));
                        }
                    }

                    for e in 0..labels.len().0 {
                        if e == i {
                            continue;
                        }
                        let ej = EdgeIndex(e);
                        let d = labels[idx] - labels[ej];
                        let dist = d.magnitude();
                        if dist > 1e-9 {
                            let dir = d / dist;
                            let repel_dist =
                                (dist - label_radii[idx] - label_radii[ej]).max(0.0) + eps;
                            force += dir * (label_charge / repel_dist.powi(2));
                        }
                    }
                }

                let mut move_vec = match cfg.label_layout {
                    LabelLayout::DanglingTangent | LabelLayout::Normal => force * step,
                    LabelLayout::FixedLength => {
                        let offset = labels[idx] - edge_pos;
                        let radial = if offset.magnitude2() > 1e-12 {
                            offset.normalize()
                        } else {
                            axes[idx]
                        };
                        let tangent = Vector2::new(-radial.y, radial.x);
                        tangent * force.dot(tangent) * step
                    }
                };
                let mag = move_vec.magnitude();
                if max_delta > 0.0 && mag > max_delta {
                    move_vec *= max_delta / mag;
                }
                labels[idx] += move_vec;
                let offset = labels[idx] - edge_pos;
                match cfg.label_layout {
                    LabelLayout::DanglingTangent | LabelLayout::Normal => {
                        if offset.dot(axes[idx]) < 0.0 {
                            let dist = offset.magnitude();
                            labels[idx] = edge_pos + axes[idx] * dist;
                        }
                    }
                    LabelLayout::FixedLength => {
                        let radius = label_length.abs();
                        labels[idx] = if radius <= 1e-12 {
                            edge_pos
                        } else if offset.magnitude2() > 1e-12 {
                            edge_pos + offset.normalize() * radius
                        } else {
                            edge_pos + axes[idx] * label_length
                        };
                    }
                }
                max_move = max_move.max(move_vec.magnitude());
            }

            if max_move < cfg.label_early_tol {
                break;
            }
        }

        let label_gap = (spring_length * 0.08).max(0.04);
        let mut measured_labels = base_labels;
        self.separate_edge_label_positions_from_boxes(
            &mut measured_labels,
            &axes,
            label_length,
            label_gap,
            spring_length,
        );
        for i in 0..labels.len().0 {
            let idx = EdgeIndex(i);
            let (half_width, half_height) = self.edge_label_route_half_extents(idx, label_gap);
            if half_width > label_gap || half_height > label_gap {
                labels[idx] = measured_labels[idx];
            }
        }

        for i in 0..labels.len().0 {
            let idx = EdgeIndex(i);
            self.graph[idx].label_pos = Some(labels[idx]);
            let angle = self.edge_label_angle(idx);
            self.graph[idx].label_angle = Some(angle);
        }
    }

    fn clear_edge_label_positions(&mut self) {
        for i in 0..self.graph.n_edges() {
            self.graph[EdgeIndex(i)].label_pos = None;
        }
    }

    fn clear_hedge_route_points(&mut self) {
        for hedge in 0..self.graph.n_hedges() {
            self.graph[Hedge(hedge)].route_points.clear();
        }
    }

    fn apply_layered_edge_routes(&mut self, output: &LayeredOutput) {
        for (edge, route) in output.edge_routes.iter() {
            self.apply_layered_edge_route(edge, route);
        }
    }

    fn apply_layered_edge_route(&mut self, edge: EdgeIndex, route: &LayeredEdgeRoute) {
        let (_, pair) = self.graph[&edge];
        match pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                self.graph[source].route_points = route.source.clone();
                self.graph[sink].route_points = route.sink.clone();
            }
            HedgePair::Unpaired {
                hedge,
                flow: Flow::Source,
            } => {
                self.graph[hedge].route_points = route.source.clone();
            }
            HedgePair::Unpaired {
                hedge,
                flow: Flow::Sink,
            } => {
                self.graph[hedge].route_points = route.sink.clone();
            }
        }
    }

    fn separate_edge_label_positions_from_boxes(
        &self,
        labels: &mut EdgeVec<Point2<f64>>,
        axes: &EdgeVec<Vector2<f64>>,
        label_length: f64,
        label_gap: f64,
        spring_length: f64,
    ) {
        let pos_v = self.new_nodevec(|_, _, node| node.pos);
        let node_boxes = self.layout_node_boxes(&pos_v, label_gap);
        let mut ordered_labels = labels
            .iter()
            .filter_map(|(edge, &target)| {
                let (half_width, half_height) = self.edge_label_route_half_extents(edge, label_gap);
                (half_width > label_gap && half_height > label_gap).then_some((
                    edge,
                    target,
                    half_width,
                    half_height,
                    half_width * half_height,
                ))
            })
            .collect::<Vec<_>>();

        ordered_labels.sort_by(|left, right| {
            right
                .4
                .partial_cmp(&left.4)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    right
                        .1
                        .y
                        .partial_cmp(&left.1.y)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| {
                    left.1
                        .x
                        .partial_cmp(&right.1.x)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| left.0 .0.cmp(&right.0 .0))
        });

        let mut placed = Vec::<LayoutRect>::new();
        for (edge, base_target, half_width, half_height, _) in ordered_labels {
            let axis = Self::normalized_or(axes[edge], Vector2::unit_y());
            let tangent = Vector2::new(-axis.y, axis.x);
            let edge_pos = self.graph[edge].pos;
            let normal_step = (2.0 * half_height + label_gap).max(spring_length * 0.10);
            let tangent_step = (2.0 * half_width + label_gap).max(spring_length * 0.10);
            let min_axis_distance = (label_length.abs() * 0.20).min(label_length.abs());
            let candidates = Self::edge_label_collision_candidates(
                base_target,
                edge_pos,
                axis,
                tangent,
                normal_step,
                tangent_step,
                min_axis_distance,
            );

            let mut best_target = base_target;
            let mut best_score = f64::INFINITY;
            for candidate in candidates {
                let candidate_box = LayoutRect::centered(candidate, half_width, half_height);
                let overlap_score = node_boxes
                    .iter()
                    .chain(placed.iter())
                    .map(|other| candidate_box.overlap_area(other))
                    .sum::<f64>();
                if overlap_score == 0.0 {
                    best_target = candidate;
                    break;
                }

                let distance_score = (candidate - base_target).magnitude2() * 1e-3;
                let score = overlap_score + distance_score;
                if score < best_score {
                    best_score = score;
                    best_target = candidate;
                }
            }

            labels[edge] = best_target;
            placed.push(LayoutRect::centered(best_target, half_width, half_height));
        }
    }

    fn edge_label_collision_candidates(
        base_target: Point2<f64>,
        edge_pos: Point2<f64>,
        axis: Vector2<f64>,
        tangent: Vector2<f64>,
        normal_step: f64,
        tangent_step: f64,
        min_axis_distance: f64,
    ) -> Vec<Point2<f64>> {
        let mut candidates = vec![base_target];
        for normal_lane in 0..=12 {
            let normal_offset = normal_lane as f64 * normal_step;
            for tangent_lane in -12i32..=12 {
                if normal_lane == 0 && tangent_lane == 0 {
                    continue;
                }
                let tangent_offset = f64::from(tangent_lane) * tangent_step;
                let candidate = base_target + axis * normal_offset + tangent * tangent_offset;
                if (candidate - edge_pos).dot(axis) + 1e-9 >= min_axis_distance {
                    candidates.push(candidate);
                }
            }
        }
        candidates.sort_by(|left, right| {
            let left_dist = (*left - base_target).magnitude2();
            let right_dist = (*right - base_target).magnitude2();
            left_dist
                .partial_cmp(&right_dist)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        candidates
    }

    fn normalized_or(value: Vector2<f64>, fallback: Vector2<f64>) -> Vector2<f64> {
        if value.magnitude2() > 1e-12 {
            value.normalize()
        } else {
            fallback
        }
    }

    fn node_layout_radii(&self) -> NodeVec<f64> {
        let widths = self.node_layout_extents("layout-width");
        let heights = self.node_layout_extents("layout-height");
        self.new_nodevec(|node, _, _| 0.5 * widths[node].hypot(heights[node]))
    }

    fn edge_label_radii(&self) -> EdgeVec<f64> {
        self.new_edgevec(|edge, _, _pair| {
            let width =
                Self::positive_statement_f64(&edge.statements, "label-width").unwrap_or(0.0);
            let height =
                Self::positive_statement_f64(&edge.statements, "label-height").unwrap_or(0.0);
            0.5 * width.hypot(height)
        })
    }

    fn edge_label_axis(&self, edge: EdgeIndex, pair: &HedgePair) -> Vector2<f64> {
        if matches!(
            self.layout_config.label_layout,
            LabelLayout::DanglingTangent
        ) && matches!(pair, HedgePair::Unpaired { .. })
        {
            self.edge_label_tangent(edge, pair)
        } else {
            self.edge_label_normal(edge, pair)
        }
    }

    fn edge_label_tangent(&self, edge: EdgeIndex, pair: &HedgePair) -> Vector2<f64> {
        let edge_pos = self.graph[edge].pos;
        let mut tangent = match pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                let a = self[self.node_id(*source)].pos;
                let b = self[self.node_id(*sink)].pos;
                b - a
            }
            HedgePair::Unpaired { hedge, .. } => {
                let a = self[self.node_id(*hedge)].pos;
                edge_pos - a
            }
        };

        if tangent.magnitude2() <= 1e-12 {
            tangent = Vector2::new(1.0, 0.0);
        }

        tangent.normalize()
    }

    fn edge_label_normal(&self, edge: EdgeIndex, pair: &HedgePair) -> Vector2<f64> {
        let edge_pos = self.graph[edge].pos;
        let tangent = self.edge_label_tangent(edge, pair);
        let mut normal = Vector2::new(-tangent.y, tangent.x);
        if normal.magnitude2() <= 1e-12 {
            normal = Vector2::new(0.0, 1.0);
        }

        if let HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } = pair {
            let a = self[self.node_id(*source)].pos;
            let b = self[self.node_id(*sink)].pos;
            let mid = a + (b - a) * 0.5;
            let to_ctrl = edge_pos - mid;
            let cross = tangent.x * to_ctrl.y - tangent.y * to_ctrl.x;
            if cross < 0.0 {
                normal = -normal;
            }
        }

        normal.normalize()
    }

    fn edge_label_angle(&self, edge: EdgeIndex) -> f64 {
        let (_, pair) = &self.graph[&edge];
        let mut dir = match pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                let a = self[self.node_id(*source)].pos;
                let b = self[self.node_id(*sink)].pos;
                b - a
            }
            HedgePair::Unpaired { hedge, .. } => {
                let a = self[self.node_id(*hedge)].pos;
                self.graph[edge].pos - a
            }
        };

        if dir.magnitude2() <= 1e-12 {
            dir = Vector2::new(1.0, 0.0);
        }

        let mut angle = dir.y.atan2(dir.x);
        if angle > std::f64::consts::FRAC_PI_2 {
            angle -= std::f64::consts::PI;
        } else if angle < -std::f64::consts::FRAC_PI_2 {
            angle += std::f64::consts::PI;
        }
        angle
    }

    fn apply_target_point<I, R>(
        constraints: &PointConstraint,
        start_x: bool,
        start_y: bool,
        target: Point2<f64>,
        fallback: Vector2<f64>,
        index: I,
        points: &mut R,
    ) where
        I: From<usize> + PartialEq + Copy,
        R: IndexMut<I, Output = Point2<f64>>,
    {
        match (constraints.x, constraints.y) {
            (Constraint::Free, Constraint::Free) => {
                if !start_x {
                    points[index].x = target.x;
                }
                if !start_y {
                    points[index].y = target.y;
                }
            }
            (Constraint::Fixed, Constraint::Fixed) => {}
            (Constraint::Free, Constraint::Fixed) => {
                if !start_x {
                    points[index].x = target.x;
                }
            }
            (Constraint::Fixed, Constraint::Free) => {
                if !start_y {
                    points[index].y = target.y;
                }
            }
            (Constraint::Grouped(_, x_dir), Constraint::Free) => {
                if !start_y {
                    points[index].y = target.y;
                }
                let target_x = Self::directional_target(target.x, x_dir, fallback.x);
                constraints.shift((target_x, 0.0).into(), index, points);
            }
            (Constraint::Free, Constraint::Grouped(_, y_dir)) => {
                if !start_x {
                    points[index].x = target.x;
                }
                let target_y = Self::directional_target(target.y, y_dir, fallback.y);
                constraints.shift((0.0, target_y).into(), index, points);
            }
            (Constraint::Grouped(_, x_dir), Constraint::Grouped(_, y_dir)) => {
                let target_x = Self::directional_target(target.x, x_dir, fallback.x);
                let target_y = Self::directional_target(target.y, y_dir, fallback.y);
                constraints.shift((target_x, target_y).into(), index, points);
            }
            (Constraint::Fixed, Constraint::Grouped(_, y_dir)) => {
                let target_y = Self::directional_target(target.y, y_dir, fallback.y);
                constraints.shift((0.0, target_y).into(), index, points);
            }
            (Constraint::Grouped(_, x_dir), Constraint::Fixed) => {
                let target_x = Self::directional_target(target.x, x_dir, fallback.x);
                constraints.shift((target_x, 0.0).into(), index, points);
            }
        }
    }

    fn selected_layout_nodes(
        &self,
        subgraph: &SuBitGraph,
        include_all_nodes: bool,
    ) -> NodeVec<bool> {
        let mut included = self.new_nodevec(|_, _, _| include_all_nodes);
        if !include_all_nodes {
            for hedge in subgraph.included_iter() {
                included[self.node_id(hedge)] = true;
            }
        }
        included
    }

    fn spanning_forest_of(&self, subgraph: &SuBitGraph) -> SuBitGraph {
        self.cycle_basis_of(subgraph).1
    }

    fn direct_layout_positions(
        &mut self,
        cfg: TreeInitCfg,
        algo: LayoutAlgo,
    ) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        let full = self.full_filter();
        let layout_subgraph = match algo {
            LayoutAlgo::Tree => self.spanning_forest_of(&full),
            LayoutAlgo::Dot | LayoutAlgo::StableLayered => full.clone(),
            LayoutAlgo::Anneal | LayoutAlgo::Force => unreachable!(),
        };
        self.layout_positions_for_subgraph(cfg, algo, &layout_subgraph, &full, true)
    }

    fn partial_layout_positions(
        &mut self,
        cfg: TreeInitCfg,
        algo: LayoutAlgo,
        subgraph: &SuBitGraph,
    ) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        let layout_subgraph = match algo {
            LayoutAlgo::Tree => self.spanning_forest_of(subgraph),
            LayoutAlgo::Dot | LayoutAlgo::StableLayered => subgraph.clone(),
            LayoutAlgo::Anneal | LayoutAlgo::Force => unreachable!(),
        };
        self.layout_positions_for_subgraph(cfg, algo, &layout_subgraph, subgraph, false)
    }

    fn fixed_node_direct_layout_positions(
        &mut self,
        cfg: TreeInitCfg,
        algo: LayoutAlgo,
    ) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        self.fixed_node_partial_layout_positions(cfg, algo, None)
    }

    fn fixed_node_partial_layout_positions(
        &mut self,
        cfg: TreeInitCfg,
        algo: LayoutAlgo,
        selected_edges: Option<&EdgeVec<bool>>,
    ) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        match algo {
            LayoutAlgo::Tree => {
                self.edge_layout_positions_from_current_nodes(cfg, selected_edges, false)
            }
            LayoutAlgo::Dot | LayoutAlgo::StableLayered => {
                self.dot_edge_layout_positions_from_current_nodes(cfg, selected_edges, false)
            }
            LayoutAlgo::Anneal | LayoutAlgo::Force => unreachable!(),
        }
    }

    fn partial_optimized_positions(
        &mut self,
        cfg: TreeInitCfg,
        subgraph: &SuBitGraph,
        energy: &SpringChargeEnergy,
    ) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        let (mut selected_nodes, selected_edges) = self.selected_layout_items(subgraph);
        let (pos_n, pos_e) = if self.layout_config.layout_nodes.nodes_are_fixed() {
            selected_nodes = self.new_nodevec(|_, _, _| false);
            self.edge_layout_positions_from_current_nodes(cfg, Some(&selected_edges), true)
        } else {
            let (mut pos_n, mut pos_e) =
                self.partial_layout_positions(cfg, LayoutAlgo::Tree, subgraph);

            for (node, selected) in selected_nodes.iter() {
                if !selected {
                    pos_n[node] = self.graph[node].pos;
                }
            }
            for (edge, selected) in selected_edges.iter() {
                if !selected {
                    pos_e[edge] = self.graph[edge].pos;
                }
            }

            (pos_n, pos_e)
        };

        let saved_node_constraints = self.new_nodevec(|_, _, node| node.constraints);
        let saved_edge_constraints = self.new_edgevec(|edge, _, _| edge.constraints);
        self.freeze_unselected_layout_items(&selected_nodes, &selected_edges);
        let positions = self.optimized_positions(pos_n, pos_e, energy);
        self.restore_layout_constraints(saved_node_constraints, saved_edge_constraints);
        positions
    }

    fn selected_layout_items(&self, subgraph: &SuBitGraph) -> (NodeVec<bool>, EdgeVec<bool>) {
        let nodes = self.selected_layout_nodes(subgraph, false);
        let mut edges = self.new_edgevec(|_, _, _| false);
        for (pair, _, _) in self.iter_edges_of(subgraph) {
            edges[self[&pair.any_hedge()]] = true;
        }
        (nodes, edges)
    }

    fn freeze_unselected_layout_items(
        &mut self,
        selected_nodes: &NodeVec<bool>,
        selected_edges: &EdgeVec<bool>,
    ) {
        let fixed = PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Fixed,
        };
        for (node, selected) in selected_nodes.iter() {
            if !selected {
                self.graph[node].constraints = fixed;
            }
        }
        for (edge, selected) in selected_edges.iter() {
            if !selected {
                self.graph[edge].constraints = fixed;
            }
        }
    }

    fn restore_layout_constraints(
        &mut self,
        node_constraints: NodeVec<PointConstraint>,
        edge_constraints: EdgeVec<PointConstraint>,
    ) {
        for (node, constraints) in node_constraints.into_iter() {
            self.graph[node].constraints = constraints;
        }
        for (edge, constraints) in edge_constraints.into_iter() {
            self.graph[edge].constraints = constraints;
        }
    }

    fn layout_positions_for_subgraph(
        &mut self,
        cfg: TreeInitCfg,
        algo: LayoutAlgo,
        subgraph: &SuBitGraph,
        node_subgraph: &SuBitGraph,
        include_all_nodes: bool,
    ) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        match algo {
            LayoutAlgo::Tree => {
                let targets =
                    self.tree_layout_targets(cfg, subgraph, node_subgraph, include_all_nodes);
                self.positions_from_node_targets(cfg, targets)
            }
            LayoutAlgo::Dot => self.layered_layout_positions(
                cfg,
                subgraph,
                include_all_nodes,
                LayeredProfile::Dot,
                true,
            ),
            LayoutAlgo::StableLayered => self.layered_layout_positions(
                cfg,
                subgraph,
                include_all_nodes,
                LayeredProfile::Stable,
                true,
            ),
            LayoutAlgo::Anneal | LayoutAlgo::Force => unreachable!(),
        }
    }

    fn layered_layout_positions(
        &mut self,
        cfg: TreeInitCfg,
        subgraph: &SuBitGraph,
        include_all_nodes: bool,
        profile: LayeredProfile,
        route_unselected_edges: bool,
    ) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        let output = self.layered_layout_output(cfg, subgraph, include_all_nodes, profile);
        self.apply_layered_edge_routes(&output);
        self.positions_from_layered_output(cfg, output, true, route_unselected_edges)
    }

    fn layered_layout_output(
        &self,
        cfg: TreeInitCfg,
        subgraph: &SuBitGraph,
        include_all_nodes: bool,
        profile: LayeredProfile,
    ) -> LayeredOutput {
        let config = LayeredConfig {
            profile,
            layer_gap: cfg.dy,
            node_gap: cfg.dx,
            edge_gap: cfg.dx * 0.35,
            route_edge_weight: self.layout_config.route_edge_weight,
            route_exit_weight: self.layout_config.route_exit_weight,
            route_label_width_scale: self.layout_config.route_label_width_scale,
            route_label_width_cap: self.layout_config.route_label_width_cap,
            sweeps: self.layout_config.schedule.epochs.clamp(1, 24),
            roots: self
                .layout_config
                .layout_roots
                .iter()
                .copied()
                .filter(|&root| root < self.n_nodes())
                .map(NodeIndex)
                .collect(),
            rank_same: self.rank_same_node_groups(),
            include_all_nodes,
        };
        let geometry = self.layered_geometry();
        self.graph.layered_layout(subgraph, &config, &geometry)
    }

    fn layered_geometry(&self) -> LayeredGeometry {
        LayeredGeometry {
            node_widths: self.node_layout_extents("layout-width"),
            node_heights: self.node_layout_extents("layout-height"),
            edge_label_widths: self.edge_layout_extents("label-width"),
            edge_label_heights: self.edge_layout_extents("label-height"),
            edge_minlens: self.new_edgevec(|edge, _, _| {
                Self::positive_statement_usize(&edge.statements, "minlen").unwrap_or(1)
            }),
            edge_weights: self.new_edgevec(|edge, _, _| {
                Self::positive_statement_f64(&edge.statements, "weight").unwrap_or(1.0)
            }),
            edge_source_exits: self.new_edgevec(|edge, _, pair| {
                Self::route_exit_statement(&edge.statements, "source-route-exit")
                    .or_else(|| self.route_exit_from_compass(pair, true))
                    .unwrap_or_default()
            }),
            edge_sink_exits: self.new_edgevec(|edge, _, pair| {
                Self::route_exit_statement(&edge.statements, "sink-route-exit")
                    .or_else(|| self.route_exit_from_compass(pair, false))
                    .unwrap_or_default()
            }),
            edge_constrained: self.new_edgevec(|edge, _, _| {
                Self::statement_bool(&edge.statements, "constraint").unwrap_or(true)
            }),
        }
    }

    fn route_exit_from_compass(
        &self,
        pair: &HedgePair,
        source_side: bool,
    ) -> Option<LayeredRouteExit> {
        let hedge = match (source_side, pair) {
            (true, HedgePair::Paired { source, .. } | HedgePair::Split { source, .. }) => source,
            (false, HedgePair::Paired { sink, .. } | HedgePair::Split { sink, .. }) => sink,
            (_, HedgePair::Unpaired { .. }) => return None,
        };
        self.graph[*hedge]
            .compasspt
            .as_deref()
            .and_then(Self::parse_route_exit)
    }

    fn rank_same_node_groups(&self) -> Vec<Vec<NodeIndex>> {
        self.layout_config
            .rank_same
            .iter()
            .filter_map(|label| SuBitGraph::from_base62(label, self.n_hedges()))
            .map(|subgraph| {
                let mut seen = vec![false; self.n_nodes()];
                let mut nodes = Vec::new();
                for hedge in subgraph.included_iter() {
                    let node = self.node_id(hedge);
                    if !seen[node.0] {
                        seen[node.0] = true;
                        nodes.push(node);
                    }
                }
                nodes
            })
            .filter(|nodes| !nodes.is_empty())
            .collect()
    }

    fn positions_from_layered_output(
        &self,
        cfg: TreeInitCfg,
        output: LayeredOutput,
        respect_edge_start: bool,
        route_unselected_edges: bool,
    ) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        let mut pos_v = self.new_nodevec(|_, _, n| n.pos);
        let mut pos_e = self.new_edgevec(|e, _, _| e.pos);

        for (node, target) in output.node_positions.iter() {
            let Some(target) = *target else {
                continue;
            };
            Self::apply_target_point(
                &self[node].constraints,
                self[node].start_x,
                self[node].start_y,
                target,
                Vector2::new(cfg.dx, cfg.dy),
                node,
                &mut pos_v,
            );
        }

        for (edge, target) in output.edge_positions.iter() {
            let Some(target) = *target else {
                continue;
            };
            let (start_x, start_y) = if respect_edge_start {
                (self[edge].start_x, self[edge].start_y)
            } else {
                (false, false)
            };
            Self::apply_target_point(
                &self[edge].constraints,
                start_x,
                start_y,
                target,
                Vector2::new(cfg.dx * 0.5, cfg.dy * 0.5),
                edge,
                &mut pos_e,
            );
        }

        if route_unselected_edges {
            let route_edges = self.new_edgevec(|_, edge, _| output.edge_positions[edge].is_none());
            let (_, routed_edges) = self.dot_edge_layout_positions_from_nodes(
                cfg,
                Some(&route_edges),
                true,
                pos_v.clone(),
            );
            for (edge, route_edge) in route_edges.iter() {
                if *route_edge {
                    pos_e[edge] = routed_edges[edge];
                }
            }
        }

        self.apply_grouped_constraints(&mut pos_v, &mut pos_e);
        (pos_v, pos_e)
    }

    fn tree_layout_targets(
        &self,
        cfg: TreeInitCfg,
        subgraph: &SuBitGraph,
        node_subgraph: &SuBitGraph,
        include_all_nodes: bool,
    ) -> NodeVec<Option<Point2<f64>>> {
        let included = self.selected_layout_nodes(node_subgraph, include_all_nodes);
        let mut adjacency = self.new_nodevec(|_, _, _| Vec::<NodeIndex>::new());
        for (pair, _, _) in self.iter_edges_of(subgraph) {
            match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    let source_node = self.node_id(source);
                    let sink_node = self.node_id(sink);
                    if source_node != sink_node {
                        adjacency[source_node].push(sink_node);
                        adjacency[sink_node].push(source_node);
                    }
                }
                HedgePair::Unpaired { .. } => {}
            }
        }

        let ordered_nodes = self.ordered_layout_nodes(&included);
        let order = self.layout_order_positions(&ordered_nodes);
        for (node, children) in adjacency.iter_mut() {
            if included[node] {
                children.sort_by_key(|child| order[child.0]);
            }
        }

        let mut visited = self.new_nodevec(|_, _, _| false);
        let mut depth = self.new_nodevec(|_, _, _| 0usize);
        let mut children = self.new_nodevec(|_, _, _| Vec::<NodeIndex>::new());
        let mut roots = Vec::new();

        for root_node in ordered_nodes {
            if !included[root_node] || visited[root_node] {
                continue;
            }

            let mut q = std::collections::VecDeque::new();
            visited[root_node] = true;
            q.push_back(root_node);
            roots.push(root_node);
            while let Some(v) = q.pop_front() {
                for &u in &adjacency[v] {
                    if !included[u] || visited[u] {
                        continue;
                    }
                    visited[u] = true;
                    depth[u] = depth[v] + 1;
                    children[v].push(u);
                    q.push_back(u);
                }
            }
        }

        let widths = self.node_layout_extents("layout-width");
        let heights = self.node_layout_extents("layout-height");
        let mut subtree_widths = self.new_nodevec(|_, _, _| 0.0);
        for &root in &roots {
            Self::compute_tidy_tree_width(root, &children, &widths, cfg.dx, &mut subtree_widths);
        }

        let mut x_slots = self.new_nodevec(|_, _, _| None::<f64>);
        let mut cursor = 0.0;
        for (root_index, &root) in roots.iter().enumerate() {
            if root_index > 0 {
                cursor += cfg.dx;
            }
            Self::assign_tidy_tree_x(
                root,
                &children,
                &subtree_widths,
                cfg.dx,
                cursor,
                &mut x_slots,
            );
            cursor += subtree_widths[root];
        }

        let level_y = self.layer_y_positions(&included, &depth, &heights, cfg.dy);
        let mut targets = self.new_nodevec(|_, _, _| None);
        for (node, slot) in x_slots.iter() {
            if let Some(slot) = *slot {
                targets[node] = Some(Point2::new(slot, level_y[depth[node]]));
            }
        }
        Self::center_targets_x_by_extents(&mut targets, &widths);
        targets
    }

    fn compute_tidy_tree_width(
        node: NodeIndex,
        children: &NodeVec<Vec<NodeIndex>>,
        widths: &NodeVec<f64>,
        gap: f64,
        subtree_widths: &mut NodeVec<f64>,
    ) -> f64 {
        let children_width = if children[node].is_empty() {
            0.0
        } else {
            children[node]
                .iter()
                .map(|&child| {
                    Self::compute_tidy_tree_width(child, children, widths, gap, subtree_widths)
                })
                .sum::<f64>()
                + gap * children[node].len().saturating_sub(1) as f64
        };
        let width = widths[node].max(children_width);
        subtree_widths[node] = width;
        width
    }

    fn assign_tidy_tree_x(
        node: NodeIndex,
        children: &NodeVec<Vec<NodeIndex>>,
        subtree_widths: &NodeVec<f64>,
        gap: f64,
        left: f64,
        x_slots: &mut NodeVec<Option<f64>>,
    ) -> f64 {
        let subtree_width = subtree_widths[node];
        let x = left + 0.5 * subtree_width;
        x_slots[node] = Some(x);

        if !children[node].is_empty() {
            let children_width = children[node]
                .iter()
                .map(|&child| subtree_widths[child])
                .sum::<f64>()
                + gap * children[node].len().saturating_sub(1) as f64;
            let mut child_left = left + 0.5 * (subtree_width - children_width).max(0.0);
            for &child in &children[node] {
                Self::assign_tidy_tree_x(child, children, subtree_widths, gap, child_left, x_slots);
                child_left += subtree_widths[child] + gap;
            }
        }
        x
    }

    fn layer_y_positions(
        &self,
        included: &NodeVec<bool>,
        ranks: &NodeVec<usize>,
        heights: &NodeVec<f64>,
        gap: f64,
    ) -> Vec<f64> {
        let max_rank = ranks
            .iter()
            .filter_map(|(node, rank)| included[node].then_some(*rank))
            .max()
            .unwrap_or(0);
        let mut rank_heights: Vec<f64> = vec![0.0; max_rank + 1];
        for (node, &is_included) in included.iter() {
            if is_included {
                rank_heights[ranks[node]] = rank_heights[ranks[node]].max(heights[node]);
            }
        }
        Self::layer_positions_from_heights(&rank_heights, gap)
    }

    fn layer_positions_from_heights(heights: &[f64], gap: f64) -> Vec<f64> {
        let mut positions = vec![0.0; heights.len()];
        for index in 1..heights.len() {
            positions[index] =
                positions[index - 1] - (0.5 * heights[index - 1] + gap + 0.5 * heights[index]);
        }
        positions
    }

    fn node_layout_extents(&self, key: &str) -> NodeVec<f64> {
        self.new_nodevec(|_, _, node| {
            Self::positive_statement_f64(&node.statements, key).unwrap_or(0.0)
        })
    }

    fn edge_layout_extents(&self, key: &str) -> EdgeVec<f64> {
        self.new_edgevec(|edge, _, _| {
            Self::positive_statement_f64(&edge.statements, key).unwrap_or(0.0)
        })
    }

    fn positive_statement_f64(statements: &BTreeMap<String, String>, key: &str) -> Option<f64> {
        let value = dot_statement_value(statements, key)?
            .trim()
            .trim_matches('"')
            .parse::<f64>()
            .ok()?;
        value.is_finite().then_some(value.max(0.0))
    }

    fn positive_statement_usize(statements: &BTreeMap<String, String>, key: &str) -> Option<usize> {
        let value = dot_statement_value(statements, key)?
            .trim()
            .trim_matches('"')
            .parse::<usize>()
            .ok()?;
        Some(value)
    }

    fn statement_bool(statements: &BTreeMap<String, String>, key: &str) -> Option<bool> {
        dot_statement_value(statements, key).and_then(|value| {
            match value.trim().trim_matches('"').to_ascii_lowercase().as_str() {
                "true" | "on" | "yes" | "1" => Some(true),
                "false" | "off" | "no" | "0" => Some(false),
                _ => None,
            }
        })
    }

    fn route_exit_statement(
        statements: &BTreeMap<String, String>,
        key: &str,
    ) -> Option<LayeredRouteExit> {
        dot_statement_value(statements, key)
            .map(|value| value.trim().trim_matches('"'))
            .and_then(Self::parse_route_exit)
    }

    fn parse_route_exit(value: &str) -> Option<LayeredRouteExit> {
        match value.trim().to_ascii_lowercase().as_str() {
            "auto" | "_" | "center" | "c" => Some(LayeredRouteExit::Auto),
            "up" | "top" | "north" | "n" => Some(LayeredRouteExit::Up),
            "down" | "bottom" | "south" | "s" => Some(LayeredRouteExit::Down),
            _ => None,
        }
    }

    fn center_targets_x_by_extents(
        targets: &mut NodeVec<Option<Point2<f64>>>,
        widths: &NodeVec<f64>,
    ) {
        let mut min_x = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        for (node, target) in targets.iter() {
            if let Some(target) = target {
                min_x = min_x.min(target.x - 0.5 * widths[node]);
                max_x = max_x.max(target.x + 0.5 * widths[node]);
            }
        }
        if !min_x.is_finite() || !max_x.is_finite() {
            return;
        }

        let center = 0.5 * (min_x + max_x);
        for (_, target) in targets.iter_mut() {
            if let Some(target) = target {
                target.x -= center;
            }
        }
    }

    fn ordered_layout_nodes(&self, included: &NodeVec<bool>) -> Vec<NodeIndex> {
        let mut seen = vec![false; included.len().0];
        let mut nodes = Vec::new();

        for &root in &self.layout_config.layout_roots {
            if root < included.len().0 {
                let node = NodeIndex(root);
                if included[node] && !seen[root] {
                    seen[root] = true;
                    nodes.push(node);
                }
            }
        }

        for (node, &is_included) in included.iter() {
            if is_included && !seen[node.0] {
                seen[node.0] = true;
                nodes.push(node);
            }
        }

        nodes
    }

    fn layout_order_positions(&self, ordered_nodes: &[NodeIndex]) -> Vec<usize> {
        let mut order = vec![usize::MAX; self.n_nodes()];
        for (position, node) in ordered_nodes.iter().enumerate() {
            order[node.0] = position;
        }
        order
    }

    fn positions_from_node_targets(
        &self,
        cfg: TreeInitCfg,
        targets: NodeVec<Option<Point2<f64>>>,
    ) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        let mut pos_v = self.new_nodevec(|_, _, n| n.pos);
        let mut pos_e = self.new_edgevec(|e, _, _| e.pos);

        for (node, target) in targets.iter() {
            let Some(target) = *target else {
                continue;
            };
            Self::apply_target_point(
                &self[node].constraints,
                self[node].start_x,
                self[node].start_y,
                target,
                Vector2::new(cfg.dx, cfg.dy),
                node,
                &mut pos_v,
            );
        }

        for (pair, _, _) in self.iter_edges() {
            let eid = self[&pair.any_hedge()];
            let target = match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    let a = pos_v[self.node_id(source)];
                    let b = pos_v[self.node_id(sink)];
                    a.midpoint(b)
                }
                HedgePair::Unpaired { hedge, .. } => {
                    let node = self.node_id(hedge);
                    pos_v[node] + Vector2::new(0.0, cfg.dy * 0.5)
                }
            };
            Self::apply_target_point(
                &self[eid].constraints,
                self[eid].start_x,
                self[eid].start_y,
                target,
                Vector2::new(cfg.dx * 0.5, cfg.dy * 0.5),
                eid,
                &mut pos_e,
            );
        }

        self.apply_grouped_constraints(&mut pos_v, &mut pos_e);
        (pos_v, pos_e)
    }

    fn edge_layout_positions_from_current_nodes(
        &self,
        cfg: TreeInitCfg,
        selected_edges: Option<&EdgeVec<bool>>,
        respect_start: bool,
    ) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        let pos_v = self.new_nodevec(|_, _, n| n.pos);
        let mut pos_e = self.new_edgevec(|e, _, _| e.pos);

        for (pair, _, _) in self.iter_edges() {
            let eid = self[&pair.any_hedge()];
            if let Some(selected_edges) = selected_edges {
                if !selected_edges[eid] {
                    continue;
                }
            }

            let target = match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    let a = pos_v[self.node_id(source)];
                    let b = pos_v[self.node_id(sink)];
                    a.midpoint(b)
                }
                HedgePair::Unpaired { hedge, .. } => {
                    let node = self.node_id(hedge);
                    pos_v[node] + Vector2::new(0.0, cfg.dy * 0.5)
                }
            };
            let (start_x, start_y) = if respect_start {
                (self[eid].start_x, self[eid].start_y)
            } else {
                (false, false)
            };
            Self::apply_target_point(
                &self[eid].constraints,
                start_x,
                start_y,
                target,
                Vector2::new(cfg.dx * 0.5, cfg.dy * 0.5),
                eid,
                &mut pos_e,
            );
        }

        (pos_v, pos_e)
    }

    fn dot_edge_layout_positions_from_current_nodes(
        &self,
        cfg: TreeInitCfg,
        selected_edges: Option<&EdgeVec<bool>>,
        respect_start: bool,
    ) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        let pos_v = self.new_nodevec(|_, _, n| n.pos);
        self.dot_edge_layout_positions_from_nodes(cfg, selected_edges, respect_start, pos_v)
    }

    fn dot_edge_layout_positions_from_nodes(
        &self,
        cfg: TreeInitCfg,
        selected_edges: Option<&EdgeVec<bool>>,
        respect_start: bool,
        pos_v: NodeVec<Point2<f64>>,
    ) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        let mut pos_e = self.new_edgevec(|e, _, _| e.pos);
        let mut targets = self.new_edgevec(|_, _, _| None::<Point2<f64>>);
        let mut paired_groups = BTreeMap::<(usize, usize), Vec<EdgeIndex>>::new();
        let mut same_rank_groups = BTreeMap::<i64, Vec<EdgeIndex>>::new();
        let mut dangling_groups = BTreeMap::<(usize, Flow), Vec<EdgeIndex>>::new();

        for (pair, _, _) in self.iter_edges() {
            let eid = self[&pair.any_hedge()];
            if let Some(selected_edges) = selected_edges {
                if !selected_edges[eid] {
                    continue;
                }
            }

            match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    let source_node = self.node_id(source);
                    let sink_node = self.node_id(sink);
                    let source_pos = pos_v[source_node];
                    let sink_pos = pos_v[sink_node];
                    if let Some(rank) = Self::same_rank_key(source_pos, sink_pos) {
                        same_rank_groups.entry(rank).or_default().push(eid);
                    } else {
                        paired_groups
                            .entry((source_node.0, sink_node.0))
                            .or_default()
                            .push(eid);
                    }
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let node = self.node_id(hedge);
                    dangling_groups.entry((node.0, flow)).or_default().push(eid);
                }
            }
        }

        for edges in paired_groups.values_mut() {
            edges.sort_by_key(|edge| edge.0);
            let offsets =
                self.centered_edge_label_offsets(edges, cfg.dx * 0.35, cfg.dx * 0.12, true);
            for (&eid, offset) in edges.iter().zip(offsets) {
                let (_, pair) = &self.graph[&eid];
                let (source, sink) = match pair {
                    HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                        (*source, *sink)
                    }
                    HedgePair::Unpaired { .. } => unreachable!(),
                };
                let a = pos_v[self.node_id(source)];
                let b = pos_v[self.node_id(sink)];
                let normal = Self::edge_route_normal(a, b);
                targets[eid] = Some(a.midpoint(b) + normal * offset);
            }
        }

        for edges in same_rank_groups.values_mut() {
            edges.sort_by(|a, b| {
                let ax = self.edge_midpoint(*a, &pos_v).x;
                let bx = self.edge_midpoint(*b, &pos_v).x;
                ax.partial_cmp(&bx)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.0.cmp(&b.0))
            });
            let lane_gap = cfg.dx * 0.12;
            let default_spacing = cfg.dx * 0.35;
            let lane_step = self.same_rank_label_lane_step(edges, default_spacing, lane_gap);
            let mut lane_right_edges = Vec::<f64>::new();
            for &eid in edges.iter() {
                let midpoint = self.edge_midpoint(eid, &pos_v);
                let half_width = self
                    .edge_label_route_half_extent(eid, false)
                    .max(default_spacing * 0.5);
                let left = midpoint.x - half_width;
                let right = midpoint.x + half_width;
                let lane = lane_right_edges
                    .iter()
                    .position(|&lane_right| left >= lane_right + lane_gap)
                    .unwrap_or_else(|| {
                        lane_right_edges.push(f64::NEG_INFINITY);
                        lane_right_edges.len() - 1
                    });
                lane_right_edges[lane] = right;
                targets[eid] = Some(midpoint + Vector2::new(0.0, (lane + 1) as f64 * lane_step));
            }
        }

        for ((node, _flow), edges) in dangling_groups.iter_mut() {
            edges.sort_by_key(|edge| edge.0);
            let offsets =
                self.centered_edge_label_offsets(edges, cfg.dx * 0.35, cfg.dx * 0.12, false);
            let anchor = pos_v[NodeIndex(*node)];
            for (&eid, offset) in edges.iter().zip(offsets) {
                targets[eid] = Some(anchor + Vector2::new(offset, cfg.dy * 0.5));
            }
        }

        self.separate_edge_label_targets_from_boxes(&mut targets, &pos_v, cfg);

        for (eid, target) in targets.iter() {
            let Some(target) = *target else {
                continue;
            };
            let (start_x, start_y) = if respect_start {
                (self[eid].start_x, self[eid].start_y)
            } else {
                (false, false)
            };
            Self::apply_target_point(
                &self[eid].constraints,
                start_x,
                start_y,
                target,
                Vector2::new(cfg.dx * 0.5, cfg.dy * 0.5),
                eid,
                &mut pos_e,
            );
        }

        (pos_v, pos_e)
    }

    fn separate_edge_label_targets_from_boxes(
        &self,
        targets: &mut EdgeVec<Option<Point2<f64>>>,
        pos_v: &NodeVec<Point2<f64>>,
        cfg: TreeInitCfg,
    ) {
        let label_gap = (cfg.dx * 0.12).max(0.05);
        let node_boxes = self.layout_node_boxes(pos_v, label_gap);
        let mut placed_labels = Vec::<LayoutRect>::new();
        let mut ordered_targets = targets
            .iter()
            .filter_map(|(edge, target)| target.map(|target| (edge, target)))
            .collect::<Vec<_>>();
        ordered_targets.sort_by(|(left_edge, left), (right_edge, right)| {
            right
                .y
                .partial_cmp(&left.y)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    left.x
                        .partial_cmp(&right.x)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| left_edge.0.cmp(&right_edge.0))
        });

        for (edge, base_target) in ordered_targets {
            let (half_width, half_height) = self.edge_label_route_half_extents(edge, label_gap);
            if half_width <= label_gap && half_height <= label_gap {
                continue;
            }

            let normal = self.edge_route_normal_for_edge(edge, pos_v);
            let midpoint = self.edge_midpoint(edge, pos_v);
            let side = if (base_target - midpoint).dot(normal) < 0.0 {
                -1.0
            } else {
                1.0
            };
            let step = (2.0 * half_height + label_gap).max(cfg.dx * 0.35).max(0.1);
            let mut best_target = base_target;

            'candidates: for lane in 0..24 {
                let lane_offset = lane as f64 * step;
                let offsets = if lane == 0 {
                    vec![0.0]
                } else if side < 0.0 {
                    vec![-lane_offset, lane_offset]
                } else {
                    vec![lane_offset, -lane_offset]
                };

                for offset in offsets {
                    let candidate = base_target + normal * offset;
                    let candidate_box = LayoutRect::centered(candidate, half_width, half_height);
                    let overlaps_node = node_boxes
                        .iter()
                        .any(|node_box| candidate_box.overlaps(node_box));
                    let overlaps_label = placed_labels
                        .iter()
                        .any(|label_box| candidate_box.overlaps(label_box));
                    if !overlaps_node && !overlaps_label {
                        best_target = candidate;
                        break 'candidates;
                    }
                }
            }

            targets[edge] = Some(best_target);
            placed_labels.push(LayoutRect::centered(best_target, half_width, half_height));
        }
    }

    fn layout_node_boxes(&self, pos_v: &NodeVec<Point2<f64>>, padding: f64) -> Vec<LayoutRect> {
        let widths = self.node_layout_extents("layout-width");
        let heights = self.node_layout_extents("layout-height");
        widths
            .iter()
            .filter_map(|(node, &width)| {
                let height = heights[node];
                (width > 0.0 || height > 0.0).then(|| {
                    LayoutRect::centered(pos_v[node], 0.5 * width + padding, 0.5 * height + padding)
                })
            })
            .collect()
    }

    fn same_rank_key(a: Point2<f64>, b: Point2<f64>) -> Option<i64> {
        if (a.y - b.y).abs() <= 1e-6 {
            Some(((a.y + b.y) * 500_000.0).round() as i64)
        } else {
            None
        }
    }

    fn edge_midpoint(&self, edge: EdgeIndex, pos_v: &NodeVec<Point2<f64>>) -> Point2<f64> {
        let (_, pair) = &self.graph[&edge];
        match pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                let a = pos_v[self.node_id(*source)];
                let b = pos_v[self.node_id(*sink)];
                a.midpoint(b)
            }
            HedgePair::Unpaired { .. } => self.graph[edge].pos,
        }
    }

    fn edge_route_normal_for_edge(
        &self,
        edge: EdgeIndex,
        pos_v: &NodeVec<Point2<f64>>,
    ) -> Vector2<f64> {
        let (_, pair) = &self.graph[&edge];
        match pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                let a = pos_v[self.node_id(*source)];
                let b = pos_v[self.node_id(*sink)];
                Self::edge_route_normal(a, b)
            }
            HedgePair::Unpaired { .. } => Vector2::unit_x(),
        }
    }

    fn same_rank_label_lane_step(
        &self,
        edges: &[EdgeIndex],
        default_spacing: f64,
        label_gap: f64,
    ) -> f64 {
        let max_label_height = edges
            .iter()
            .map(|&edge| {
                Self::positive_statement_f64(&self.graph[edge].statements, "label-height")
                    .unwrap_or(0.0)
            })
            .fold(0.0, f64::max);
        default_spacing.max(max_label_height + label_gap.max(0.0))
    }

    fn centered_edge_label_offsets(
        &self,
        edges: &[EdgeIndex],
        default_spacing: f64,
        label_gap: f64,
        radial_extent: bool,
    ) -> Vec<f64> {
        if edges.is_empty() {
            return Vec::new();
        }

        let half_extents = edges
            .iter()
            .map(|&edge| self.edge_label_route_half_extent(edge, radial_extent))
            .collect::<Vec<_>>();

        let mut offsets = vec![0.0; edges.len()];
        for index in 1..edges.len() {
            let measured_spacing =
                half_extents[index - 1] + half_extents[index] + label_gap.max(0.0);
            offsets[index] = offsets[index - 1] + default_spacing.max(measured_spacing);
        }

        let min = offsets.first().copied().unwrap_or(0.0) - half_extents[0];
        let max =
            offsets.last().copied().unwrap_or(0.0) + half_extents.last().copied().unwrap_or(0.0);
        let center = 0.5 * (min + max);
        for offset in &mut offsets {
            *offset -= center;
        }
        offsets
    }

    fn edge_label_route_half_extent(&self, edge: EdgeIndex, radial_extent: bool) -> f64 {
        let edge = &self.graph[edge];
        let width = Self::positive_statement_f64(&edge.statements, "label-width").unwrap_or(0.0);
        let height = Self::positive_statement_f64(&edge.statements, "label-height").unwrap_or(0.0);
        if radial_extent {
            0.5 * width.hypot(height)
        } else {
            0.5 * width
        }
    }

    fn edge_label_route_half_extents(&self, edge: EdgeIndex, padding: f64) -> (f64, f64) {
        let edge = &self.graph[edge];
        let width = Self::positive_statement_f64(&edge.statements, "label-width").unwrap_or(0.0);
        let height = Self::positive_statement_f64(&edge.statements, "label-height").unwrap_or(0.0);
        (0.5 * width + padding, 0.5 * height + padding)
    }

    fn edge_route_normal(a: Point2<f64>, b: Point2<f64>) -> Vector2<f64> {
        let direction = b - a;
        let length = direction.magnitude();
        if length <= 1e-9 {
            Vector2::unit_x()
        } else {
            Vector2::new(-direction.y, direction.x) / length
        }
    }

    /// Generate new positions based on the default traversal-tree placement.
    pub fn new_positions(&self, cfg: TreeInitCfg) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        let full = self.full_filter();
        let forest = self.spanning_forest_of(&full);
        let targets = self.tree_layout_targets(cfg, &forest, &full, true);
        self.positions_from_node_targets(cfg, targets)
    }

    pub fn parse<'a>(dot_str: &str) -> Result<Self, HedgeParseError<'a, (), (), (), ()>> {
        let dot = DotGraph::from_string(dot_str)?;
        Ok(TypstGraph::from_dot(dot, &default_figment()))
    }

    pub fn to_cbor(&self) -> CBORTypstGraph {
        CBORTypstGraph {
            edges: self.new_edgevec(|e, i, _p| EdgeData::new(e.clone(), self.orientation(i))),
            nodes: self.new_nodevec(|_id, _h, v| v.clone()),
            global_eval: self.global_eval.clone(),
            name: self.name.clone(),
            data: self.data.clone(),
            global_statements: self.global_statements.clone(),
            default_edge_statements: self.default_edge_statements.clone(),
            default_node_statements: self.default_node_statements.clone(),
        }
    }

    pub fn serialize_to_file<P: AsRef<Path>>(
        &self,
        path: P,
    ) -> Result<(), Box<dyn std::error::Error>> {
        self.to_cbor().serialize_to_file(path)
    }

    pub fn to_dot_graph(&self) -> DotGraph {
        let graph = self.graph.map_data_ref(
            |_graph, _neighbors, node| node.to_dot(),
            |_graph, _eid, _pair, edge_data| edge_data.map(|e| e.to_dot()),
            |_hedge, hedge_data| hedge_data.to_dot(),
        );

        let mut global_data = GlobalData::from(());
        global_data.name = self.name.clone();
        global_data.payload = self.data.clone();
        global_data.statements = self
            .global_statements
            .iter()
            .filter(|(key, _)| !LayoutConfig::is_layout_statement_key(key))
            .map(|(key, value)| (key.clone(), value.clone()))
            .collect();
        global_data.edge_statements = self.default_edge_statements.clone();
        global_data.node_statements = self.default_node_statements.clone();
        if let Some(eval) = &self.global_eval {
            global_data
                .statements
                .insert("eval".to_string(), eval.clone());
        }
        DotGraph { graph, global_data }
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub struct CBORTypstGraph {
    edges: EdgeVec<EdgeData<TypstEdge>>,
    nodes: NodeVec<TypstNode>,
    global_eval: Option<String>,
    name: String,
    data: Option<Vec<u8>>,
    global_statements: BTreeMap<String, String>,
    default_edge_statements: BTreeMap<String, String>,
    default_node_statements: BTreeMap<String, String>,
}

impl CBORTypstGraph {
    pub fn serialize_to_file<P: AsRef<Path>>(
        &self,
        path: P,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::create(path)?;
        let writer = BufWriter::new(file);
        ciborium::ser::into_writer(self, writer)?;
        Ok(())
    }

    pub fn deserialize_from_file<P: AsRef<Path>>(
        path: P,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(path)?;
        let graph = ciborium::de::from_reader(file)?;
        Ok(graph)
    }
}
