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
    builder_add_edge_bytes, builder_add_node_bytes, builder_finish_bytes, builder_new_bytes,
    graph_archived_compass_subgraph_bytes, graph_archived_subgraph_bytes,
    graph_compass_subgraph_bytes, graph_cycle_basis_bytes, graph_dot_bytes, graph_edges_bytes,
    graph_edges_of_archived_subgraph_bytes, graph_edges_of_bytes, graph_from_spec_bytes,
    graph_info_bytes, graph_join_by_edge_key_bytes, graph_join_by_hedge_key_bytes,
    graph_nodes_bytes, graph_nodes_of_archived_subgraph_bytes, graph_nodes_of_bytes,
    graph_spanning_forests_bytes, graph_subgraph_bytes, subgraph_contains_hedge_bytes,
    subgraph_hedges_bytes, subgraph_label_bytes, TypstDotEdge, TypstDotEndpoint, TypstDotGraphInfo,
    TypstDotNode, TypstPoint,
};
pub use pin::{expand_template, PinConstraint};

use cgmath::{EuclideanSpace, InnerSpace, Point2, Rad, Vector2, Zero};
use dot_parser::ast::CompassPt;
use figment::{providers::Serialized, Figment, Profile};
use linnet::half_edge::swap::Swap;
use linnet::{
    half_edge::{
        involution::{EdgeData, EdgeIndex, EdgeVec, Flow, Hedge, HedgePair, Involution},
        layout::{
            force::{force_directed_layout, ForceLayoutConfig},
            simulatedanneale::{anneal, GeoSchedule, SAConfig},
            spring::{
                Constraint, HasPointConstraint, LayoutState, ParamTuning, PinnedLayoutNeighbor,
                PointConstraint, ShiftDirection, Shiftable, SpringChargeEnergy,
            },
        },
        nodestore::{DefaultNodeStore, NodeStorageOps},
        EdgeAccessors, HedgeGraph, NodeIndex, NodeVec,
    },
    parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData, GlobalData, HedgeParseError},
};

use rand::rngs::SmallRng;
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::{collections::BTreeMap, io::BufWriter};
use std::{collections::HashMap, ops::Deref};
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

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct TypstNode {
    name: Option<String>,
    index: Option<NodeIndex>,
    pos: Point2<f64>,
    constraints: PointConstraint,
    shift: Option<Vector2<f64>>,
    statements: BTreeMap<String, String>,
    eval: Option<String>,
}

impl Default for TypstNode {
    fn default() -> Self {
        TypstNode {
            name: None,
            index: None,
            pos: Point2::origin(),
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

        let mut eval: Option<String> = data.get("eval").transpose().unwrap();

        // Apply template expansion and clean quotes
        eval = eval.map(|template| {
            let clean_template = template
                .strip_prefix('"')
                .unwrap_or(&template)
                .strip_suffix('"')
                .unwrap_or(&template);
            expand_template(clean_template, &data.statements)
        });

        let (pos, constraints) = init_points[nid];

        Self {
            name: data.name,
            index: data.index,
            statements: data.statements,
            pos,
            constraints,
            shift,
            eval,
        }
    }

    fn parse_position(
        statements: &std::collections::BTreeMap<String, String>,
        attr: &str,
    ) -> Option<(f64, f64)> {
        if let Some(value) = statements.get(attr) {
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
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct TypstEdge {
    from: Option<(NodeIndex, Hedge)>,
    to: Option<(NodeIndex, Hedge)>,
    bend: Result<Rad<f64>, GeomError>,
    pos: Point2<f64>,
    label_pos: Option<Point2<f64>>,
    label_angle: Option<f64>,
    eval_sink: Option<String>,
    eval_source: Option<String>,
    eval_label: Option<String>,
    mom_eval: Option<String>,
    shift: Option<Vector2<f64>>,
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
            bend: Err(GeomError::NotComputed),
            pos: Point2::origin(),
            label_pos: None,
            label_angle: None,
            shift: None,
            eval_label: None,
            eval_sink: None,
            eval_source: None,
            mom_eval: None,
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

            let label_pos = TypstNode::parse_position(&d.statements, "label_pos")
                .map(|(x, y)| Point2::new(x, y));
            let label_angle = d.statements.get("label_angle").and_then(|value| {
                let unquoted = value.trim().trim_matches('"');
                let cleaned = unquoted.trim().trim_end_matches("rad");
                cleaned.trim().parse::<f64>().ok()
            });

            let mut eval_sink: Option<String> = d.get("eval_sink").transpose().unwrap();

            // Apply template expansion and clean quotes for eval
            eval_sink = eval_sink.map(|template| {
                let clean_template = template
                    .strip_prefix('"')
                    .unwrap_or(&template)
                    .strip_suffix('"')
                    .unwrap_or(&template);
                expand_template(clean_template, &d.statements)
            });

            let mut eval_source: Option<String> = d.get("eval_source").transpose().unwrap();

            // Apply template expansion and clean quotes for eval
            eval_source = eval_source.map(|template| {
                let clean_template = template
                    .strip_prefix('"')
                    .unwrap_or(&template)
                    .strip_suffix('"')
                    .unwrap_or(&template);
                expand_template(clean_template, &d.statements)
            });

            let mut eval_label: Option<String> = d.get("eval_label").transpose().unwrap();

            // Apply template expansion and clean quotes for eval
            eval_label = eval_label.map(|template| {
                let clean_template = template
                    .strip_prefix('"')
                    .unwrap_or(&template)
                    .strip_suffix('"')
                    .unwrap_or(&template);
                expand_template(clean_template, &d.statements)
            });

            let mut mom_eval: Option<String> = d.get("mom_eval").transpose().unwrap();

            // Apply template expansion and clean quotes for mom_eval
            mom_eval = mom_eval.map(|template| {
                let clean_template = template
                    .strip_prefix('"')
                    .unwrap_or(&template)
                    .strip_suffix('"')
                    .unwrap_or(&template);
                expand_template(clean_template, &d.statements)
            });
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

            Self {
                from,
                to,
                pos,
                constraints,
                label_pos,
                label_angle,
                eval_label,
                eval_sink,
                eval_source,
                mom_eval,
                shift,
                statements: d.statements,
                ..Default::default()
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
            statements.insert("label_pos".to_string(), format!("{},{}", p.x, p.y));
        }

        if let Some(a) = self.label_angle {
            statements.insert("label_angle".to_string(), format!("{a}rad"));
        }

        if let Some(eval) = &self.eval_sink {
            statements.insert("eval_sink".to_string(), eval.clone());
        }

        if let Some(eval) = &self.eval_source {
            statements.insert("eval_source".to_string(), eval.clone());
        }

        if let Some(eval) = &self.eval_label {
            statements.insert("eval_label".to_string(), eval.clone());
        }

        if let Some(eval) = &self.mom_eval {
            statements.insert("mom_eval".to_string(), eval.clone());
        }

        DotEdgeData {
            local_statements: statements.clone(),
            statements,
            edge_id: None,
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone, Default)]
pub struct TypstHedge {
    from: usize,
    to: usize,
    weight: f64,
    statement: Option<String>,
    id: Option<usize>,
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

#[derive(Debug, Serialize, Deserialize)]
pub struct TypstGraph {
    graph: HedgeGraph<TypstEdge, TypstNode, TypstHedge>,
    global_eval: Option<String>,
    name: String,
    global_statements: BTreeMap<String, String>,
    layout_config: LayoutConfig,
}

impl Deref for TypstGraph {
    type Target = HedgeGraph<TypstEdge, TypstNode, TypstHedge>;
    fn deref(&self) -> &Self::Target {
        &self.graph
    }
}

impl TypstGraph {
    pub fn from_dot(dot: DotGraph, _figment: &Figment) -> Self {
        let mut group_map = HashMap::new();
        let edge_pin_constrains: EdgeVec<(Point2<f64>, PointConstraint)> =
            dot.graph.new_edgevec(|e, eid, _| {
                let a = e
                    .get::<_, String>("pin")
                    .transpose()
                    .unwrap()
                    .and_then(|a| PinConstraint::parse(&a));

                if let Some(a) = a {
                    a.point_constraint(eid.0, &mut group_map)
                } else {
                    (Point2::new(0., 0.), PointConstraint::default())
                }
            });

        let mut group_map = HashMap::new();

        let node_pin_constrains: NodeVec<(Point2<f64>, PointConstraint)> =
            dot.graph.new_nodevec(|nid, _, n| {
                let a = n
                    .get::<_, String>("pin")
                    .transpose()
                    .unwrap()
                    .and_then(|a| PinConstraint::parse(&a));

                if let Some(a) = a {
                    a.point_constraint(nid.0, &mut group_map)
                } else {
                    (Point2::new(0., 0.), PointConstraint::default())
                }
            });

        let graph = dot.graph.map(
            |inv, nid, data| TypstNode::parse(inv, nid, data, &node_pin_constrains),
            |inv, store, p, eid, data| {
                TypstEdge::parse(inv, store, p, eid, data, &edge_pin_constrains)
            },
            TypstHedge::parse,
        );

        let figment = Figment::from(Serialized::from(
            dot.global_data.statements.clone(),
            Profile::Default,
        ));

        let config = LayoutConfig::from_figment(&figment);

        let mut global_eval: Option<String> = dot.global_data.statements.get("eval").cloned();

        if let Some(g) = &mut global_eval {
            let clean_template = g
                .strip_prefix('"')
                .unwrap_or(g)
                .strip_suffix('"')
                .unwrap_or(g);
            *g = expand_template(clean_template, &dot.global_data.statements);
        }

        Self {
            graph,
            global_eval,
            name: dot.global_data.name,
            global_statements: dot
                .global_data
                .statements
                .into_iter()
                .map(|(k, v)| {
                    let clean_value = v.trim().trim_matches('"').to_string();
                    (k, clean_value)
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

pub struct TreeInitCfg {
    pub dy: f64, // ≈ 1.2 * L
    pub dx: f64, // ≈ 0.9 * L
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
enum LayoutAlgo {
    Anneal,
    Force,
}

fn default_layout_algo() -> LayoutAlgo {
    LayoutAlgo::Anneal
}

#[derive(Debug, Clone, Serialize, Deserialize)]
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
    #[serde(default = "default_incremental_energy")]
    incremental_energy: bool,
    #[serde(default = "default_layout_algo")]
    layout_algo: LayoutAlgo,
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
            label_step: default_label_step(),
            label_length_scale: default_label_length_scale(),
            label_charge: default_label_charge(),
            label_spring: default_label_spring(),
            label_early_tol: default_label_early_tol(),
            label_max_delta_scale: default_label_max_delta_scale(),
            incremental_energy: default_incremental_energy(),
            layout_algo: default_layout_algo(),
            spring: SpringConfig::default(),
            schedule: ScheduleConfig::default(),
        }
    }
}

impl LayoutConfig {
    fn from_figment(figment: &Figment) -> Self {
        figment.extract::<Self>().unwrap_or_default()
    }

    fn add_to_global(&self, global_data: &mut GlobalData) {
        macro_rules! insert {
            ($key:expr, $value:expr) => {
                global_data
                    .statements
                    .insert($key.to_string(), $value.to_string());
            };
        }

        insert!("viewport_w", self.viewport_w);
        insert!("viewport_h", self.viewport_h);
        insert!("tree_dx", self.tree_dx);
        insert!("tree_dy", self.tree_dy);
        insert!("step", self.step);
        insert!("temp", self.temp);
        insert!("seed", self.seed);
        insert!("delta", self.delta);
        insert!("directional_force", self.directional_force);
        insert!("z_spring", self.z_spring);
        insert!("z_spring_growth", self.z_spring_growth);
        insert!("label_steps", self.label_steps);
        insert!("label_step", self.label_step);
        insert!("label_length_scale", self.label_length_scale);
        insert!("label_charge", self.label_charge);
        insert!("label_spring", self.label_spring);
        insert!("label_early_tol", self.label_early_tol);
        insert!("label_max_delta_scale", self.label_max_delta_scale);
        insert!(
            "layout_algo",
            match self.layout_algo {
                LayoutAlgo::Anneal => "anneal",
                LayoutAlgo::Force => "force",
            }
        );
        global_data.statements.insert(
            "incremental_energy".to_string(),
            self.incremental_energy.to_string(),
        );

        let spring_params = ParamTuning::from(&self.spring);
        spring_params.add_to_global(global_data);
        let schedule = GeoSchedule::from(&self.schedule);
        schedule.add_to_global(global_data);
    }
}

fn default_figment() -> Figment {
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
    0.2
}

fn default_temp() -> f64 {
    1.1
}

fn default_seed() -> u64 {
    42
}

fn default_delta() -> f64 {
    0.4
}

fn default_directional_force() -> f64 {
    0.0
}

fn default_z_spring() -> f64 {
    0.0
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
    0.3
}

fn default_label_charge() -> f64 {
    0.2
}

fn default_label_spring() -> f64 {
    1.0
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

fn default_crossing_penalty() -> f64 {
    0.0
}

#[derive(Debug, Clone, Serialize, Deserialize)]
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

#[derive(Debug, Clone, Serialize, Deserialize)]
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
    1.0
}

fn default_k_spring() -> f64 {
    1.0
}

fn default_beta() -> f64 {
    0.14
}

fn default_gamma_dangling() -> f64 {
    0.14
}

fn default_gamma_ev() -> f64 {
    0.20
}

fn default_gamma_ee() -> f64 {
    0.10
}

fn default_g_center() -> f64 {
    0.05
}

fn default_eps() -> f64 {
    1e-4
}

fn default_steps() -> usize {
    200
}

fn default_epochs() -> usize {
    8
}

fn default_cool() -> f64 {
    0.85
}

fn default_accept_floor() -> f64 {
    0.15
}

fn default_step_shrink() -> f64 {
    0.7
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

impl TypstGraph {
    pub fn layout(&mut self) {
        let spring_params = ParamTuning::from(&self.layout_config.spring);

        let (tree_cfg, energy) = self.tree_init_cfg(&spring_params);
        let (pos_n, pos_e) = self.new_positions(tree_cfg);
        let mut state = self.graph.new_layout_state(
            pos_n,
            pos_e,
            self.layout_config.delta,
            self.layout_config.directional_force,
            self.layout_config.incremental_energy,
        );

        let (mut vertex_points, mut edge_points) = match self.layout_config.layout_algo {
            LayoutAlgo::Anneal => {
                let mut schedule = GeoSchedule::from(&self.layout_config.schedule);
                let (out, _stats) = anneal::<_, _, _, _, SmallRng>(
                    state,
                    SAConfig {
                        temp: self.layout_config.temp,
                        step: self.layout_config.step,
                        seed: self.layout_config.seed,
                    },
                    &PinnedLayoutNeighbor,
                    &energy,
                    &mut schedule,
                );
                (out.vertex_points, out.edge_points)
            }
            LayoutAlgo::Force => {
                force_directed_layout(
                    &mut state,
                    &energy,
                    ForceLayoutConfig {
                        steps: self.layout_config.schedule.steps,
                        epochs: self.layout_config.schedule.epochs,
                        step: self.layout_config.step,
                        cool: self.layout_config.schedule.cool,
                        max_delta: self.layout_config.delta,
                        early_tol: self.layout_config.schedule.early_tol,
                        seed: self.layout_config.seed,
                        z_spring: self.layout_config.z_spring,
                        z_spring_growth: self.layout_config.z_spring_growth,
                    },
                );
                (state.vertex_points, state.edge_points)
            }
        };

        self.apply_grouped_constraints(&mut vertex_points, &mut edge_points);
        self.update_positions(vertex_points, edge_points);
        self.layout_edge_labels(energy.spring_length);
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

            self.graph[i].bend = angle;
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

        let normals: EdgeVec<Vector2<f64>> =
            self.new_edgevec(|_e, idx, pair| self.edge_label_normal(idx, pair));

        let mut labels: EdgeVec<Point2<f64>> = self.new_edgevec(|_e, idx, _pair| {
            let edge_pos = self.graph[idx].pos;
            edge_pos + normals[idx] * label_length
        });

        for _ in 0..cfg.label_steps {
            let mut max_move: f64 = 0.0;

            for i in 0..labels.len().0 {
                let idx = EdgeIndex(i);
                let mut force = Vector2::zero();
                let edge_pos = self.graph[idx].pos;
                let target = edge_pos + normals[idx] * label_length;
                if label_spring != 0.0 {
                    force += (target - labels[idx]) * label_spring;
                }

                if label_charge != 0.0 {
                    for n in 0..self.n_nodes() {
                        let ni = NodeIndex(n);
                        let np = self[ni].pos;
                        let d = labels[idx] - np;
                        let dist = d.magnitude();
                        if dist > 1e-9 {
                            let dir = d / dist;
                            force += dir * (label_charge / (dist + eps).powi(2));
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
                            force += dir * (label_charge / (dist + eps).powi(2));
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
                            force += dir * (label_charge / (dist + eps).powi(2));
                        }
                    }
                }

                let mut move_vec = force * step;
                let mag = move_vec.magnitude();
                if max_delta > 0.0 && mag > max_delta {
                    move_vec *= max_delta / mag;
                }
                labels[idx] += move_vec;
                let offset = labels[idx] - edge_pos;
                if offset.dot(normals[idx]) < 0.0 {
                    let dist = offset.magnitude();
                    labels[idx] = edge_pos + normals[idx] * dist;
                }
                max_move = max_move.max(move_vec.magnitude());
            }

            if max_move < cfg.label_early_tol {
                break;
            }
        }

        for i in 0..labels.len().0 {
            let idx = EdgeIndex(i);
            self.graph[idx].label_pos = Some(labels[idx]);
            let angle = self.edge_label_angle(idx);
            self.graph[idx].label_angle = Some(angle);
        }
    }

    fn edge_label_normal(&self, edge: EdgeIndex, pair: &HedgePair) -> Vector2<f64> {
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

    /// Generate new positions based on tree layout while respecting constraints.
    /// This method:
    /// 1. Computes a tree-based layout with proper level spacing
    /// 2. Only applies tree positions to coordinates that are free to move
    /// 3. Respects Fixed and Grouped constraints by leaving them unchanged
    /// 4. For directional constraints (+/-), ensures absolute positioning:
    ///    - `+@group` constraints result in non-negative positions
    ///    - `-@group` constraints result in non-positive positions
    ///    - This overrides tree layout if necessary to maintain directional requirements
    pub fn new_positions(&self, cfg: TreeInitCfg) -> (NodeVec<Point2<f64>>, EdgeVec<Point2<f64>>) {
        let mut pos_v = self.new_nodevec(|_, _, n| n.pos);
        let mut pos_e = self.new_edgevec(|e, _, _| e.pos);

        let forest = self
            .all_spanning_forests_of(&self.full_filter())
            .into_iter()
            .next()
            .unwrap_or_else(|| self.empty_subgraph());
        let mut adjacency = self.new_nodevec(|_, _, _| Vec::<NodeIndex>::new());
        for (pair, _, _) in self.iter_edges_of(&forest) {
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

        let mut level: NodeVec<i32> = self.new_nodevec(|_, _, _| -1);
        let mut n_per_level: Vec<usize> = vec![];
        for (root_node, _) in adjacency.iter() {
            if level[root_node] >= 0 {
                continue;
            }

            let mut q = std::collections::VecDeque::new();
            level[root_node] = 0;
            if n_per_level.is_empty() {
                n_per_level.push(1);
            } else {
                n_per_level[0] += 1;
            }
            q.push_back(root_node);
            while let Some(v) = q.pop_front() {
                let next_level = level[v] + 1;
                for &u in &adjacency[v] {
                    if level[u] < 0 {
                        level[u] = next_level;
                        if level[u] == n_per_level.len() as i32 {
                            n_per_level.push(1);
                        } else if level[u] < n_per_level.len() as i32 {
                            n_per_level[level[u] as usize] += 1;
                        } else {
                            panic!("Level out of bounds");
                        }

                        q.push_back(u);
                    }
                }
            }
        }

        let mut level_counters: Vec<usize> = vec![0; n_per_level.len()];

        for (cl, &n) in n_per_level.iter().enumerate() {
            // place on a horizontal line
            let k = n as f64;
            let width = (k - 1.0) * cfg.dx;
            let y = (cl as f64) * cfg.dy;
            for (i, &l) in level.iter() {
                if l != (cl as i32) {
                    continue;
                }
                let position_in_level = level_counters[cl];
                let x = -0.5 * width + (position_in_level as f64) * cfg.dx;
                level_counters[cl] += 1;

                // Calculate desired tree position
                let tree_position = Point2::new(x, y);

                // Handle constraints with proper absolute positioning for directional constraints
                match (&self[i].constraints.x, &self[i].constraints.y) {
                    (Constraint::Free, Constraint::Free) => {
                        pos_v[i] = tree_position;
                    }
                    (Constraint::Fixed, Constraint::Fixed) => {
                        // Keep original position - no change needed
                    }
                    (Constraint::Free, Constraint::Fixed) => {
                        pos_v[i].x = tree_position.x;
                        // y stays as original constraint position
                    }
                    (Constraint::Fixed, Constraint::Free) => {
                        pos_v[i].y = tree_position.y;
                        // x stays as original constraint position
                    }
                    (Constraint::Grouped(_, x_dir), Constraint::Free) => {
                        // For grouped x with direction, ensure absolute position respects direction
                        pos_v[i].y = tree_position.y;
                        match x_dir {
                            ShiftDirection::PositiveOnly => {
                                // Ensure x position is positive
                                let target_x = if tree_position.x >= 0.0 {
                                    tree_position.x
                                } else {
                                    cfg.dx.abs()
                                };
                                self[i]
                                    .constraints
                                    .shift((target_x, 0.0).into(), i, &mut pos_v);
                            }
                            ShiftDirection::NegativeOnly => {
                                // Ensure x position is negative
                                let target_x = if tree_position.x <= 0.0 {
                                    tree_position.x
                                } else {
                                    -cfg.dx.abs()
                                };
                                self[i]
                                    .constraints
                                    .shift((target_x, 0.0).into(), i, &mut pos_v);
                            }
                            ShiftDirection::Any => {
                                self[i].constraints.shift(
                                    (tree_position.x, 0.0).into(),
                                    i,
                                    &mut pos_v,
                                );
                            }
                        }
                    }
                    (Constraint::Free, Constraint::Grouped(_, y_dir)) => {
                        // For grouped y with direction, ensure absolute position respects direction
                        pos_v[i].x = tree_position.x;
                        match y_dir {
                            ShiftDirection::PositiveOnly => {
                                // Ensure y position is positive
                                let target_y = if tree_position.y >= 0.0 {
                                    tree_position.y
                                } else {
                                    cfg.dy.abs()
                                };
                                self[i]
                                    .constraints
                                    .shift((0.0, target_y).into(), i, &mut pos_v);
                            }
                            ShiftDirection::NegativeOnly => {
                                // Ensure y position is negative
                                let target_y = if tree_position.y <= 0.0 {
                                    tree_position.y
                                } else {
                                    -cfg.dy.abs()
                                };
                                self[i]
                                    .constraints
                                    .shift((0.0, target_y).into(), i, &mut pos_v);
                            }
                            ShiftDirection::Any => {
                                self[i].constraints.shift(
                                    (0.0, tree_position.y).into(),
                                    i,
                                    &mut pos_v,
                                );
                            }
                        }
                    }
                    (Constraint::Grouped(_, x_dir), Constraint::Grouped(_, y_dir)) => {
                        // For both coordinates grouped with direction
                        let target_x = match x_dir {
                            ShiftDirection::PositiveOnly => {
                                if tree_position.x >= 0.0 {
                                    tree_position.x
                                } else {
                                    cfg.dx.abs()
                                }
                            }
                            ShiftDirection::NegativeOnly => {
                                if tree_position.x <= 0.0 {
                                    tree_position.x
                                } else {
                                    -cfg.dx.abs()
                                }
                            }
                            ShiftDirection::Any => tree_position.x,
                        };
                        let target_y = match y_dir {
                            ShiftDirection::PositiveOnly => {
                                if tree_position.y >= 0.0 {
                                    tree_position.y
                                } else {
                                    cfg.dy.abs()
                                }
                            }
                            ShiftDirection::NegativeOnly => {
                                if tree_position.y <= 0.0 {
                                    tree_position.y
                                } else {
                                    -cfg.dy.abs()
                                }
                            }
                            ShiftDirection::Any => tree_position.y,
                        };
                        self[i]
                            .constraints
                            .shift((target_x, target_y).into(), i, &mut pos_v);
                    }
                    (Constraint::Fixed, Constraint::Grouped(_, y_dir)) => {
                        let target_y = match y_dir {
                            ShiftDirection::PositiveOnly => {
                                if tree_position.y >= 0.0 {
                                    tree_position.y
                                } else {
                                    cfg.dy.abs()
                                }
                            }
                            ShiftDirection::NegativeOnly => {
                                if tree_position.y <= 0.0 {
                                    tree_position.y
                                } else {
                                    -cfg.dy.abs()
                                }
                            }
                            ShiftDirection::Any => tree_position.y,
                        };
                        self[i]
                            .constraints
                            .shift((0.0, target_y).into(), i, &mut pos_v);
                    }
                    (Constraint::Grouped(_, x_dir), Constraint::Fixed) => {
                        let target_x = match x_dir {
                            ShiftDirection::PositiveOnly => {
                                if tree_position.x >= 0.0 {
                                    tree_position.x
                                } else {
                                    cfg.dx.abs()
                                }
                            }
                            ShiftDirection::NegativeOnly => {
                                if tree_position.x <= 0.0 {
                                    tree_position.x
                                } else {
                                    -cfg.dx.abs()
                                }
                            }
                            ShiftDirection::Any => tree_position.x,
                        };
                        self[i]
                            .constraints
                            .shift((target_x, 0.0).into(), i, &mut pos_v);
                    }
                }
            }
        }

        // 4) edge control points
        for (pair, _, _) in self.iter_edges() {
            let h = pair.any_hedge();
            let eid = self[&h];
            match pair {
                // internal edge: midpoint + perpendicular bulge
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    let a = pos_v[self.node_id(source)];
                    let b = pos_v[self.node_id(sink)];
                    let mid = a.midpoint(b);

                    // Handle edge constraints with proper absolute positioning for directional constraints
                    match (&self[eid].constraints.x, &self[eid].constraints.y) {
                        (Constraint::Free, Constraint::Free) => {
                            pos_e[eid] = mid;
                        }
                        (Constraint::Fixed, Constraint::Fixed) => {
                            // Keep original position
                        }
                        (Constraint::Free, Constraint::Fixed) => {
                            pos_e[eid].x = mid.x;
                        }
                        (Constraint::Fixed, Constraint::Free) => {
                            pos_e[eid].y = mid.y;
                        }
                        (Constraint::Grouped(_, x_dir), Constraint::Free) => {
                            pos_e[eid].y = mid.y;
                            match x_dir {
                                ShiftDirection::PositiveOnly => {
                                    // Ensure x position is positive
                                    let target_x = if mid.x >= 0.0 {
                                        mid.x
                                    } else {
                                        cfg.dx.abs() * 0.5
                                    };
                                    self[eid].constraints.shift(
                                        (target_x, 0.0).into(),
                                        eid,
                                        &mut pos_e,
                                    );
                                }
                                ShiftDirection::NegativeOnly => {
                                    // Ensure x position is negative
                                    let target_x = if mid.x <= 0.0 {
                                        mid.x
                                    } else {
                                        -cfg.dx.abs() * 0.5
                                    };
                                    self[eid].constraints.shift(
                                        (target_x, 0.0).into(),
                                        eid,
                                        &mut pos_e,
                                    );
                                }
                                ShiftDirection::Any => {
                                    self[eid].constraints.shift(
                                        (mid.x, 0.0).into(),
                                        eid,
                                        &mut pos_e,
                                    );
                                }
                            }
                        }
                        (Constraint::Free, Constraint::Grouped(_, y_dir)) => {
                            pos_e[eid].x = mid.x;
                            match y_dir {
                                ShiftDirection::PositiveOnly => {
                                    // Ensure y position is positive
                                    let target_y = if mid.y >= 0.0 {
                                        mid.y
                                    } else {
                                        cfg.dy.abs() * 0.5
                                    };
                                    self[eid].constraints.shift(
                                        (0.0, target_y).into(),
                                        eid,
                                        &mut pos_e,
                                    );
                                }
                                ShiftDirection::NegativeOnly => {
                                    // Ensure y position is negative
                                    let target_y = if mid.y <= 0.0 {
                                        mid.y
                                    } else {
                                        -cfg.dy.abs() * 0.5
                                    };
                                    self[eid].constraints.shift(
                                        (0.0, target_y).into(),
                                        eid,
                                        &mut pos_e,
                                    );
                                }
                                ShiftDirection::Any => {
                                    self[eid].constraints.shift(
                                        (0.0, mid.y).into(),
                                        eid,
                                        &mut pos_e,
                                    );
                                }
                            }
                        }
                        (Constraint::Grouped(_, x_dir), Constraint::Grouped(_, y_dir)) => {
                            let target_x = match x_dir {
                                ShiftDirection::PositiveOnly => {
                                    if mid.x >= 0.0 {
                                        mid.x
                                    } else {
                                        cfg.dx.abs() * 0.5
                                    }
                                }
                                ShiftDirection::NegativeOnly => {
                                    if mid.x <= 0.0 {
                                        mid.x
                                    } else {
                                        -cfg.dx.abs() * 0.5
                                    }
                                }
                                ShiftDirection::Any => mid.x,
                            };
                            let target_y = match y_dir {
                                ShiftDirection::PositiveOnly => {
                                    if mid.y >= 0.0 {
                                        mid.y
                                    } else {
                                        cfg.dy.abs() * 0.5
                                    }
                                }
                                ShiftDirection::NegativeOnly => {
                                    if mid.y <= 0.0 {
                                        mid.y
                                    } else {
                                        -cfg.dy.abs() * 0.5
                                    }
                                }
                                ShiftDirection::Any => mid.y,
                            };
                            self[eid].constraints.shift(
                                (target_x, target_y).into(),
                                eid,
                                &mut pos_e,
                            );
                        }
                        (Constraint::Fixed, Constraint::Grouped(_, y_dir)) => {
                            let target_y = match y_dir {
                                ShiftDirection::PositiveOnly => {
                                    if mid.y >= 0.0 {
                                        mid.y
                                    } else {
                                        cfg.dy.abs() * 0.5
                                    }
                                }
                                ShiftDirection::NegativeOnly => {
                                    if mid.y <= 0.0 {
                                        mid.y
                                    } else {
                                        -cfg.dy.abs() * 0.5
                                    }
                                }
                                ShiftDirection::Any => mid.y,
                            };
                            self[eid]
                                .constraints
                                .shift((0.0, target_y).into(), eid, &mut pos_e);
                        }
                        (Constraint::Grouped(_, x_dir), Constraint::Fixed) => {
                            let target_x = match x_dir {
                                ShiftDirection::PositiveOnly => {
                                    if mid.x >= 0.0 {
                                        mid.x
                                    } else {
                                        cfg.dx.abs() * 0.5
                                    }
                                }
                                ShiftDirection::NegativeOnly => {
                                    if mid.x <= 0.0 {
                                        mid.x
                                    } else {
                                        -cfg.dx.abs() * 0.5
                                    }
                                }
                                ShiftDirection::Any => mid.x,
                            };
                            self[eid]
                                .constraints
                                .shift((target_x, 0.0).into(), eid, &mut pos_e);
                        }
                    }
                }
                HedgePair::Unpaired { hedge, .. } => {
                    let v = self.node_id(hedge);
                    let default_pos = Point2::new(pos_v[v].x + 1., pos_v[v].y + 1.);

                    // Handle unpaired edge constraints with proper absolute positioning for directional constraints
                    match (&self[eid].constraints.x, &self[eid].constraints.y) {
                        (Constraint::Free, Constraint::Free) => {
                            pos_e[eid] = default_pos;
                        }
                        (Constraint::Fixed, Constraint::Fixed) => {
                            // Keep original position
                        }
                        (Constraint::Free, Constraint::Fixed) => {
                            pos_e[eid].x = default_pos.x;
                        }
                        (Constraint::Fixed, Constraint::Free) => {
                            pos_e[eid].y = default_pos.y;
                        }
                        (Constraint::Grouped(_, x_dir), Constraint::Free) => {
                            pos_e[eid].y = default_pos.y;
                            match x_dir {
                                ShiftDirection::PositiveOnly => {
                                    // Ensure x position is positive
                                    let target_x = if default_pos.x >= 0.0 {
                                        default_pos.x
                                    } else {
                                        cfg.dx.abs() * 0.5
                                    };
                                    self[eid].constraints.shift(
                                        (target_x, 0.0).into(),
                                        eid,
                                        &mut pos_e,
                                    );
                                }
                                ShiftDirection::NegativeOnly => {
                                    // Ensure x position is negative
                                    let target_x = if default_pos.x <= 0.0 {
                                        default_pos.x
                                    } else {
                                        -cfg.dx.abs() * 0.5
                                    };
                                    self[eid].constraints.shift(
                                        (target_x, 0.0).into(),
                                        eid,
                                        &mut pos_e,
                                    );
                                }
                                ShiftDirection::Any => {
                                    self[eid].constraints.shift(
                                        (default_pos.x, 0.0).into(),
                                        eid,
                                        &mut pos_e,
                                    );
                                }
                            }
                        }
                        (Constraint::Free, Constraint::Grouped(_, y_dir)) => {
                            pos_e[eid].x = default_pos.x;
                            match y_dir {
                                ShiftDirection::PositiveOnly => {
                                    // Ensure y position is positive
                                    let target_y = if default_pos.y >= 0.0 {
                                        default_pos.y
                                    } else {
                                        cfg.dy.abs() * 0.5
                                    };
                                    self[eid].constraints.shift(
                                        (0.0, target_y).into(),
                                        eid,
                                        &mut pos_e,
                                    );
                                }
                                ShiftDirection::NegativeOnly => {
                                    // Ensure y position is negative
                                    let target_y = if default_pos.y <= 0.0 {
                                        default_pos.y
                                    } else {
                                        -cfg.dy.abs() * 0.5
                                    };
                                    self[eid].constraints.shift(
                                        (0.0, target_y).into(),
                                        eid,
                                        &mut pos_e,
                                    );
                                }
                                ShiftDirection::Any => {
                                    self[eid].constraints.shift(
                                        (0.0, default_pos.y).into(),
                                        eid,
                                        &mut pos_e,
                                    );
                                }
                            }
                        }
                        (Constraint::Grouped(_, x_dir), Constraint::Grouped(_, y_dir)) => {
                            let target_x = match x_dir {
                                ShiftDirection::PositiveOnly => {
                                    if default_pos.x >= 0.0 {
                                        default_pos.x
                                    } else {
                                        cfg.dx.abs() * 0.5
                                    }
                                }
                                ShiftDirection::NegativeOnly => {
                                    if default_pos.x <= 0.0 {
                                        default_pos.x
                                    } else {
                                        -cfg.dx.abs() * 0.5
                                    }
                                }
                                ShiftDirection::Any => default_pos.x,
                            };
                            let target_y = match y_dir {
                                ShiftDirection::PositiveOnly => {
                                    if default_pos.y >= 0.0 {
                                        default_pos.y
                                    } else {
                                        cfg.dy.abs() * 0.5
                                    }
                                }
                                ShiftDirection::NegativeOnly => {
                                    if default_pos.y <= 0.0 {
                                        default_pos.y
                                    } else {
                                        -cfg.dy.abs() * 0.5
                                    }
                                }
                                ShiftDirection::Any => default_pos.y,
                            };
                            self[eid].constraints.shift(
                                (target_x, target_y).into(),
                                eid,
                                &mut pos_e,
                            );
                        }
                        (Constraint::Fixed, Constraint::Grouped(_, y_dir)) => {
                            let target_y = match y_dir {
                                ShiftDirection::PositiveOnly => {
                                    if default_pos.y >= 0.0 {
                                        default_pos.y
                                    } else {
                                        cfg.dy.abs() * 0.5
                                    }
                                }
                                ShiftDirection::NegativeOnly => {
                                    if default_pos.y <= 0.0 {
                                        default_pos.y
                                    } else {
                                        -cfg.dy.abs() * 0.5
                                    }
                                }
                                ShiftDirection::Any => default_pos.y,
                            };
                            self[eid]
                                .constraints
                                .shift((0.0, target_y).into(), eid, &mut pos_e);
                        }
                        (Constraint::Grouped(_, x_dir), Constraint::Fixed) => {
                            let target_x = match x_dir {
                                ShiftDirection::PositiveOnly => {
                                    if default_pos.x >= 0.0 {
                                        default_pos.x
                                    } else {
                                        cfg.dx.abs() * 0.5
                                    }
                                }
                                ShiftDirection::NegativeOnly => {
                                    if default_pos.x <= 0.0 {
                                        default_pos.x
                                    } else {
                                        -cfg.dx.abs() * 0.5
                                    }
                                }
                                ShiftDirection::Any => default_pos.x,
                            };
                            self[eid]
                                .constraints
                                .shift((target_x, 0.0).into(), eid, &mut pos_e);
                        }
                    }
                }
            }
        }

        self.apply_grouped_constraints(&mut pos_v, &mut pos_e);
        (pos_v, pos_e)
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
            global_statements: self.global_statements.clone(),
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

        // Reconstruct GlobalData from layout parameters

        let mut global_data = GlobalData::from(());
        global_data.name = self.name.clone();
        global_data.statements = self.global_statements.clone();
        if let Some(eval) = &self.global_eval {
            global_data
                .statements
                .insert("eval".to_string(), eval.clone());
        }
        self.layout_config.add_to_global(&mut global_data);
        DotGraph { graph, global_data }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CBORTypstGraph {
    edges: EdgeVec<EdgeData<TypstEdge>>,
    nodes: NodeVec<TypstNode>,
    global_eval: Option<String>,
    name: String,
    global_statements: BTreeMap<String, String>,
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
