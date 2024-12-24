use std::sync::{Arc, Mutex};

// use crate::half_edge::drawing::Decoration;
use argmin::{
    core::{CostFunction, Executor, State},
    solver::simulatedannealing::{Anneal, SimulatedAnnealing},
};
use cgmath::{Angle, InnerSpace, Rad, Vector2, Zero};
use rand::{rngs::SmallRng, Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use super::{
    drawing::{CetzEdge, CetzString, Decoration, EdgeGeometry},
    subgraph::SubGraph,
    Flow, Hedge, HedgeGraph, Involution, InvolutiveMapping, NodeIndex, Orientation,
};

pub struct LayoutVertex<V> {
    pub data: V,
    pos: Vector2<f64>,
}

impl<V> LayoutVertex<V> {
    pub fn new(data: V, x: f64, y: f64) -> Self {
        LayoutVertex {
            data,
            pos: Vector2::new(x, y),
        }
    }
}

pub struct LayoutEdge<E> {
    pub data: E,
    geometry: EdgeGeometry,
}

impl<E> LayoutEdge<E> {
    pub fn map<E2>(self, f: impl FnOnce(E) -> E2) -> LayoutEdge<E2> {
        LayoutEdge {
            data: f(self.data),
            geometry: self.geometry,
        }
    }

    pub fn new_internal(
        data: E,
        source: &Vector2<f64>,
        sink: &Vector2<f64>,
        x: f64,
        y: f64,
    ) -> Self {
        let pos = Vector2::new(x, y);
        let (center, _) = EdgeGeometry::circumcircle([source, sink, &pos]);

        let center_to_pos = pos - center;
        let center_to_source = source - center;
        let center_to_sink = sink - center;

        let sink_angle = center_to_source.angle(center_to_sink).normalize();
        let pos_angle = center_to_source.angle(center_to_pos).normalize();

        let angle = Vector2::unit_y().angle(center_to_pos)
            + if sink_angle < pos_angle {
                Rad::turn_div_2()
            } else {
                Rad::zero()
            };

        LayoutEdge {
            data,
            geometry: EdgeGeometry::Simple {
                pos: Vector2::new(x, y),
                angle,
            },
        }
    }

    pub fn new_external(
        data: E,
        source: &Vector2<f64>,
        x: f64,
        y: f64,
        orientation: Orientation,
    ) -> Self {
        let pos = Vector2::new(x, y);
        let mut pos_angle = Vector2::unit_x().angle(pos - source);

        if let Orientation::Reversed = orientation {
        } else {
            pos_angle += Rad::turn_div_2()
        }
        LayoutEdge {
            data,
            geometry: EdgeGeometry::Simple {
                pos: Vector2::new(x, y),
                angle: pos_angle,
            },
        }
    }

    pub fn to_fancy(
        &mut self,
        source: Vector2<f64>,
        sink: Option<Vector2<f64>>,
        flow: Flow,
        settings: &FancySettings,
    ) {
        match self.geometry {
            EdgeGeometry::Simple { .. } => {}
            _ => {
                return;
            }
        };
        self.geometry = self.geometry.clone().to_fancy(source, sink, flow, settings);
    }
}

impl<E> LayoutEdge<E> {
    pub fn cetz(&self, source: usize, orientation: Orientation, sink: Option<usize>) -> String
    where
        E: CetzEdge,
    {
        self.cetz_impl(
            &|e| e.label(),
            &|e| e.decoration(),
            source,
            orientation,
            sink,
        )
    }

    pub fn cetz_impl(
        &self,
        label: &impl Fn(&E) -> String,
        decoration: &impl Fn(&E) -> Decoration,
        source: usize,
        _orientation: Orientation,
        sink: Option<usize>,
    ) -> String {
        match &self.geometry {
            EdgeGeometry::Simple { pos, angle } => {
                // if let Orientation::Reversed = orientation {
                //     angle += Rad::turn_div_2();
                //     println!("reversed")
                // }
                if let Some(sink) = sink {
                    format!(
                        "edge(node{}.pos,{},node{}.pos,decoration:{},angle:{})\n",
                        source,
                        pos.to_cetz(),
                        sink,
                        decoration(&self.data).to_cetz(),
                        angle.to_cetz()
                    )
                } else {
                    format!(
                        "edge(node{}.pos,{},decoration:{},angle:{})\n",
                        source,
                        pos.to_cetz(),
                        decoration(&self.data).to_cetz(),
                        angle.to_cetz(),
                    )
                }
            }
            EdgeGeometry::Fancy {
                pos,
                label_pos,
                label_angle,
            } => {
                if let Some(sink) = sink {
                    format!(
                    "edge(node{}.pos,{},node{}.pos,decoration:{},angle:{})\ncetz.draw.content({},angle:{},[{}])\n",
                    source,
                    pos.to_cetz(),
                    sink,
                    decoration(&self.data).to_cetz(),
                    label_angle.to_cetz(),
                    label_pos.to_cetz(),
                    label_angle.to_cetz(),
                    label(&self.data)
                )
                } else {
                    format!(
                    "edge(node{}.pos,{},decoration:{},angle:{})\ncetz.draw.content({},angle:{},[{}])\n",
                    source,
                    pos.to_cetz(),
                    decoration(&self.data).to_cetz(),
                    label_angle.to_cetz(),
                    label_pos.to_cetz(),
                    label_angle.to_cetz(),
                    label(&self.data)
                )
                }
            }
            EdgeGeometry::FancyArrow {
                pos,
                arrow_arc,
                label_pos,
                label_angle,
            } => {
                if let Some(sink) = sink {
                    format!(
                "edge(node{}.pos,{},node{}.pos,decoration:{},angle:{})\ncetz.draw.content({},angle:{},[{}])\n{}\n",
                source,
                pos.to_cetz(),
                sink,
                decoration(&self.data).to_cetz(),
                label_angle.to_cetz(),
                label_pos.to_cetz(),
                label_angle.to_cetz(),
                label(&self.data),
                arrow_arc.hobby_arrow(),
            )
                } else {
                    format!(
                "edge(node{}.pos,{},decoration:{},angle:{})\ncetz.draw.content({},angle:{},[{}])\n{}\n",
                source,
                pos.to_cetz(),
                decoration(&self.data).to_cetz(),
                label_angle.to_cetz(),
                label_pos.to_cetz(),
                label_angle.to_cetz(),
                label(&self.data),
                arrow_arc.hobby_arrow(),
                )
                }
            }
        }
    }
}

pub type PositionalHedgeGraph<E, V> = HedgeGraph<LayoutEdge<E>, LayoutVertex<V>>;

const CETZ_PREAMBLE: &str = r#"
let node(pos)=cetz.draw.circle(pos,radius:0.02,fill: black)
let stroke = 0.7pt
let amplitude = 0.051
let arrow-scale = 0.8
let segment-length = 0.0521
let edge(..points,decoration:"",angle:0deg)={
    if decoration == "coil"{
    cetz.decorations.coil(cetz.draw.hobby(..points),amplitude:amplitude,stroke:stroke,align:"MID")
    } else if decoration == "wave" {
        cetz.decorations.wave(cetz.draw.hobby(..points),amplitude:amplitude,stroke:stroke)
    }  else if decoration == "arrow"{
        let points = points.pos()
        if points.len()==2{
          let center = (0.5*(points.at(0).at(0)+points.at(1).at(0)),0.5*(points.at(0).at(1)+points.at(1).at(1)))
          cetz.draw.hobby(..points,stroke:stroke)
          cetz.draw.mark(center,angle,symbol: ">", fill: black,anchor: "center",scale:arrow-scale)
        } else {
          let (first,center,..other)=points
          cetz.draw.hobby(first,center,..other,stroke:stroke)
            cetz.draw.mark(center,angle,symbol: ">", fill: black,anchor: "center",scale:arrow-scale)
        }

    }else {
            cetz.draw.hobby(..points,stroke:stroke)
    }
}
"#;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FancySettings {
    pub label_shift: f64,
    pub arrow_angle_percentage: Option<f64>,
    pub arrow_shift: f64,
}

impl Default for FancySettings {
    fn default() -> Self {
        Self {
            label_shift: 0.1,
            arrow_angle_percentage: None,
            arrow_shift: 0.1,
        }
    }
}

impl FancySettings {
    pub fn label_shift(&self) -> f64 {
        self.label_shift
            + if self.arrow_angle_percentage.is_some() {
                self.arrow_shift
            } else {
                0.0
            }
    }
}

impl<E, V> PositionalHedgeGraph<E, V> {
    pub fn to_fancy(&mut self, settings: &FancySettings) {
        for (eid, nid) in self.involution.hedge_data.iter().enumerate() {
            let i = self[*nid].pos;

            if let InvolutiveMapping::Identity { data, underlying } = &mut self.involution.inv[eid]
            {
                if let Some(d) = &mut data.data {
                    d.to_fancy(i, None, *underlying, settings);
                }
            }

            if let Some(cn) = self.involved_node_hairs(super::Hedge(eid)) {
                let j = self[cn].pos;

                if let InvolutiveMapping::Source { data, .. } = &mut self.involution.inv[eid] {
                    if let Some(d) = &mut data.data {
                        d.to_fancy(i, Some(j), Flow::Source, settings);
                    }
                }
            }
        }
    }

    fn cetz_preamble() -> String {
        let mut out = String::new();

        out.push_str(
            r#"#import "@preview/cetz:0.3.1"
"#,
        );
        out
    }

    pub fn cetz_impl_collection(
        graphs: &[(String, String, Vec<Self>)],
        edge_label: &impl Fn(&E) -> String,
        edge_decoration: &impl Fn(&E) -> Decoration,
    ) -> String {
        let mut out = Self::cetz_preamble();
        out.push_str("#{\nlet cols = (30%,30%,30%)\n");

        out.push_str(CETZ_PREAMBLE);
        for (col_name, _, gs) in graphs.iter() {
            for (j, g) in gs.iter().enumerate() {
                out.push_str(&format!(
                    "let {col_name}{j}={}",
                    g.cetz_bare(edge_label, edge_decoration)
                ))
            }
        }
        for (col_name, label, gs) in graphs.iter() {
            out.push_str(&format!("[{}]\n", label));
            out.push_str("grid(columns: cols,gutter: 20pt,");
            for i in 0..gs.len() {
                out.push_str(&format!("{col_name}{i},"));
            }
            out.push_str(")\n");
        }
        out.push('}');
        out
    }

    fn cetz_bare(
        &self,
        edge_label: &impl Fn(&E) -> String,
        edge_decoration: &impl Fn(&E) -> Decoration,
    ) -> String {
        let mut out = String::from("cetz.canvas(length:50%,{\n");
        for a in 0..self.base_nodes {
            out.push_str(&format!(
                "let node{}= (pos:{})\n",
                a,
                self.node_data[a].pos.to_cetz()
            ));
            out.push_str(&format!("node(node{}.pos)\n", a));
        }

        for (eid, (e, nid)) in self.involution.iter() {
            let i = nid.0;

            if let Some(cn) = self.involved_node_hairs(eid) {
                let j = self.id_from_hairs(cn).unwrap().0;

                if let InvolutiveMapping::Source { data, .. } = e {
                    if let Some(l) = &data.data {
                        out.push_str(&l.cetz_impl(
                            edge_label,
                            edge_decoration,
                            i,
                            data.orientation,
                            Some(j),
                        ));
                    }
                }
            }

            if let InvolutiveMapping::Identity { data, .. } = e {
                if let Some(l) = &data.data {
                    out.push_str(&l.cetz_impl(
                        edge_label,
                        edge_decoration,
                        i,
                        data.orientation,
                        None,
                    ));
                }
            }
        }

        out.push_str("})\n");
        out
    }

    pub fn cetz_impl(
        &self,
        edge_label: &impl Fn(&E) -> String,
        edge_decoration: &impl Fn(&E) -> Decoration,
    ) -> String {
        let mut out = Self::cetz_preamble();
        out.push_str("#{\n");
        out.push_str(CETZ_PREAMBLE);
        out.push_str(&self.cetz_bare(edge_label, edge_decoration));
        out.push_str("}\n");
        out
    }
}

pub struct GraphLayout<'a, E, V> {
    pub graph: &'a HedgeGraph<E, V>,
    pub positions: Positions,
    pub params: LayoutParams,
    pub rng: Arc<Mutex<SmallRng>>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct LayoutParams {
    pub spring_constant: f64,
    pub spring_length: f64,
    pub global_edge_repulsion: f64,
    pub edge_vertex_repulsion: f64,
    pub charge_constant_e: f64,
    pub charge_constant_v: f64,
    pub external_constant: f64,
    pub central_force_constant: f64,
}

impl Default for LayoutParams {
    fn default() -> Self {
        LayoutParams {
            spring_length: 1.,
            spring_constant: 0.5,
            global_edge_repulsion: 0.15,
            edge_vertex_repulsion: 1.5,
            charge_constant_e: 0.7,
            charge_constant_v: 13.2,
            external_constant: 0.0,
            central_force_constant: 0.0,
        }
    }
}

// LayoutParams {
//     spring_length: 1.,
//     spring_constant: 0.5,
//     global_edge_repulsion: 0.15,
//     edge_vertex_repulsion: 0.015,
//     charge_constant_e: 0.4,
//     charge_constant_v: 13.2,
//     external_constant: 0.0,
//     central_force_constant: 0.0,
// }

#[derive(Debug, Clone)]
#[allow(clippy::type_complexity)]
pub struct Positions {
    vertex_positions: Vec<(Option<(f64, f64)>, usize, usize)>,
    edge_positions: Involution<(), (Option<(f64, f64)>, usize, usize)>,
}

impl Positions {
    pub fn scale(&mut self, scale: f64) {
        self.vertex_positions.iter_mut().for_each(|(a, _, _)| {
            if let Some(pos) = a {
                pos.0 /= scale;
                pos.1 /= scale;
            }
        });

        for i in 0..self.edge_positions.inv.len() {
            match &mut self.edge_positions.inv[i] {
                InvolutiveMapping::Identity { data, .. } => data.data.as_mut().map(|(a, _, _)| {
                    a.as_mut().map(|pos| {
                        pos.0 /= scale;
                        pos.1 /= scale;
                    })
                }),
                InvolutiveMapping::Source { data, .. } => data.data.as_mut().map(|(a, _, _)| {
                    a.as_mut().map(|pos| {
                        pos.0 /= scale;
                        pos.1 /= scale;
                    })
                }),
                _ => None,
            };
        }
    }
    pub fn max(&self, params: &[f64]) -> Option<f64> {
        let vertex_max = self
            .vertex_positions
            .iter()
            .map(|(v, i, j)| {
                if let Some((x, y)) = v {
                    x.abs().max(y.abs())
                } else {
                    params[*i].abs().max(params[*j].abs())
                }
            })
            .reduce(f64::max);

        let mut max_edges = None;
        for i in self.edge_positions.iter_idx() {
            let data = self.edge_positions.edge_data(i).data.map(|(a, i, j)| {
                if let Some((x, y)) = a {
                    x.abs().max(y.abs())
                } else {
                    params[i].abs().max(params[j].abs())
                }
            });
            if let Some(data) = data {
                if let Some(max) = max_edges {
                    if max < data {
                        max_edges = Some(data);
                    }
                } else {
                    max_edges = Some(data);
                }
            }
        }

        match (vertex_max, max_edges) {
            (Some(a), Some(b)) => Some(a.max(b)),
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (None, None) => None,
        }
    }
    pub fn to_graph<E, V>(
        self,
        graph: HedgeGraph<E, V>,
        params: &[f64],
    ) -> PositionalHedgeGraph<E, V> {
        // let mut new_graph = graph;
        let node_data: Vec<LayoutVertex<V>> = graph
            .node_data
            .into_iter()
            .enumerate()
            .map(|(pos, v)| {
                let pos = self.vertex_positions.get(pos).unwrap();
                LayoutVertex::new(v, params[pos.1], params[pos.2])
            })
            .collect();

        let hedge_data = graph.involution.hedge_data;

        let inv: Vec<_> = graph
            .involution
            .inv
            .into_iter()
            .enumerate()
            .map(|(i, m)| {
                let pos = &self.get_edge_position(super::Hedge(i), params);

                let mut source = node_data[hedge_data.get(i).unwrap().0].pos;
                let m = match m {
                    InvolutiveMapping::Sink { source_idx } => {
                        InvolutiveMapping::Sink { source_idx }
                    }
                    InvolutiveMapping::Source { data, sink_idx } => {
                        let mut sink = node_data[hedge_data.get(sink_idx.0).unwrap().0].pos;

                        if let Orientation::Reversed = data.orientation {
                            std::mem::swap(&mut sink, &mut source);
                        }
                        InvolutiveMapping::Source {
                            data: data
                                .map(|s| LayoutEdge::new_internal(s, &source, &sink, pos.0, pos.1)),
                            sink_idx,
                        }
                    }
                    InvolutiveMapping::Identity { data, underlying } => {
                        let orientation = data.orientation;
                        InvolutiveMapping::Identity {
                            data: data.map(|s| {
                                LayoutEdge::new_external(
                                    s,
                                    &source,
                                    pos.0,
                                    pos.1,
                                    orientation.relative_to(underlying),
                                )
                            }),
                            underlying,
                        }
                    }
                };

                m
            })
            .collect();

        HedgeGraph {
            node_data,
            nodes: graph.nodes,
            involution: Involution { inv, hedge_data },
            base_nodes: graph.base_nodes,
        }
    }
    pub fn new<E, V>(graph: &HedgeGraph<E, V>, seed: u64) -> (Vec<f64>, Self) {
        let mut rng = SmallRng::seed_from_u64(seed);
        let mut vertex_positions = Vec::new();
        let mut edge_positions = graph.involution.forgetful_map_node_data_ref(|_| ());

        let range = -1.0..1.0;
        let ext_range = 2.0..4.0;

        let mut params = Vec::new();

        for (i, h) in graph.involution.inv.iter().enumerate() {
            let j = params.len();
            if h.is_identity() {
                params.push(rng.gen_range(ext_range.clone()));
                params.push(rng.gen_range(ext_range.clone()));
            } else {
                params.push(rng.gen_range(range.clone()));
                params.push(rng.gen_range(range.clone()));
            }
            edge_positions.set_edge_data(super::Hedge(i), (None, j, j + 1));
        }

        for _ in graph.base_nodes_iter() {
            let j = params.len();
            params.push(rng.gen_range(range.clone()));
            params.push(rng.gen_range(range.clone()));
            vertex_positions.push((None, j, j + 1));
        }

        (
            params,
            Positions {
                vertex_positions,
                edge_positions,
            },
        )
    }
    pub fn circle_ext<E, V>(graph: &HedgeGraph<E, V>, seed: u64, radius: f64) -> (Vec<f64>, Self) {
        let mut rng = SmallRng::seed_from_u64(seed);
        let mut vertex_positions = Vec::new();
        let mut edge_positions = graph.involution.forgetful_map_node_data_ref(|_| ());

        let range = -1.0..1.0;

        let mut params = Vec::new();

        let mut angle = 0.;
        let angle_step = 2. * std::f64::consts::PI / f64::from(graph.n_externals() as u32);

        for i in graph.involution.iter_idx() {
            let j = params.len();
            params.push(rng.gen_range(range.clone()));
            params.push(rng.gen_range(range.clone()));
            if i == graph.involution.inv(i) {
                let x = radius * f64::cos(angle);
                let y = radius * f64::sin(angle);
                edge_positions.set_edge_data(i, (Some((x, y)), j, j + 1));
                angle += angle_step;
            } else {
                edge_positions.set_edge_data(i, (None, j, j + 1));
            }
        }

        for _ in graph.base_nodes_iter() {
            let j = params.len();
            params.push(rng.gen_range(range.clone()));
            params.push(rng.gen_range(range.clone()));
            vertex_positions.push((None, j, j + 1));
        }

        (
            params,
            Positions {
                vertex_positions,
                edge_positions,
            },
        )
    }

    pub fn iter_vertex_positions<'a>(
        &'a self,
        params: &'a [f64],
    ) -> impl Iterator<Item = (NodeIndex, (f64, f64))> + 'a {
        self.vertex_positions
            .iter()
            .enumerate()
            .map(|(i, (p, x, y))| {
                if let Some(p) = p {
                    (NodeIndex(i), (p.0, p.1))
                } else {
                    (NodeIndex(i), (params[*x], params[*y]))
                }
            })
    }

    pub fn iter_edge_positions<'a>(
        &'a self,
        params: &'a [f64],
    ) -> impl Iterator<Item = (f64, f64)> + 'a {
        self.edge_positions.iter_edge_data().filter_map(|e| {
            if let Some((p, x, y)) = e.data {
                if let Some(p) = p {
                    Some((p.0, p.1))
                } else {
                    Some((params[x], params[y]))
                }
            } else {
                None
            }
        })
    }

    pub fn get_edge_position(&self, edge: super::Hedge, params: &[f64]) -> (f64, f64) {
        let (p, ix, iy) = self.edge_positions.edge_data(edge).data.unwrap();
        if let Some(p) = p {
            (p.0, p.1)
        } else {
            (params[ix], params[iy])
        }
    }
}

impl<E, V> CostFunction for GraphLayout<'_, E, V> {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        let mut cost = 0.0;

        // global edge repulsion:
        //

        for (x, y) in self.positions.iter_edge_positions(param) {
            for (ex, ey) in self.positions.iter_edge_positions(param) {
                if ex == x && ey == y {
                    continue;
                }
                let dx = x - ex;
                let dy = y - ey;

                cost += self.params.global_edge_repulsion / (dx * dx + dy * dy).sqrt();
            }
        }

        for (node, (x, y)) in self.positions.iter_vertex_positions(param) {
            for (ex, ey) in self.positions.iter_edge_positions(param) {
                let dx = x - ex;
                let dy = y - ey;

                let dist = (dx * dx + dy * dy).sqrt();
                cost += self.params.edge_vertex_repulsion / dist;
            }
            for e in self.graph.hairs_from_id(node).hairs.included_iter() {
                let (ex, ey) = self.positions.get_edge_position(e, param);

                let dx = x - ex;
                let dy = y - ey;

                let dist = (dx * dx + dy * dy).sqrt();

                cost +=
                    0.5 * self.params.spring_constant * (self.params.spring_length - dist).powi(2);
                // cost += self.params.charge_constant_e / dist;

                for othere in self.graph.hairs_from_id(node).hairs.included_iter() {
                    let a = self.positions.edge_positions.inv(othere);
                    if e > othere && a != e {
                        let (ox, oy) = self.positions.get_edge_position(othere, param);
                        let dx = ex - ox;
                        let dy = ey - oy;
                        let dist = (dx * dx + dy * dy).sqrt();
                        cost += self.params.charge_constant_e / dist;
                    }
                }
            }

            for (other, (ox, oy)) in self.positions.iter_vertex_positions(param) {
                if node != other {
                    let dx = x - ox;
                    let dy = y - oy;
                    let dist = (dx * dx + dy * dy).sqrt();
                    cost += 0.5 * self.params.charge_constant_v / dist;

                    let dist_to_center = (ox * ox + oy * oy).sqrt();
                    if dist_to_center > 1.0 {
                        cost += 0.5 * self.params.central_force_constant / dist_to_center;
                    }

                    // cost += 0.5 * self.central_force_constant / dist_to_center;
                }
            }

            // for e in 0..self.graph.involution.inv.len() {
            //     match self.positions.edge_positions.inv[e] {
            //         InvolutiveMapping::Identity(_) => {
            //             let (ox, oy) = self.positions.get_edge_position(e, param);
            //             let dx = x - ox;
            //             let dy = y - oy;
            //             let dist = (dx * dx + dy * dy).sqrt();
            //             cost += 0.5 * self.charge_constant_v / dist;
            //             let dist_to_center = (ox * ox + oy * oy).sqrt();

            //             cost += 0.5 * self.external_constant / dist_to_center;
            //         }
            //         _ => {}
            //     }
            // }
        }

        // for othere in 0..self.graph.involution.inv.len() {
        //     for e in 0..self.graph.involution.inv.len() {
        //         match self.positions.edge_positions.inv[e].1 {
        //             InvolutiveMapping::Sink(_) => {}
        //             _ => {
        //                 if e > othere {
        //                     let (ox, oy) = self.positions.get_edge_position(othere, param);
        //                     let (ex, ey) = self.positions.get_edge_position(e, param);
        //                     let dx = ex - ox;
        //                     let dy = ey - oy;
        //                     let dist = (dx * dx + dy * dy).sqrt();
        //                     cost += self.charge_constant_e / dist;
        //                 }
        //             }
        //         }
        //     }
        // }
        Ok(cost)
    }
}

impl<E, V> Anneal for GraphLayout<'_, E, V> {
    type Param = Vec<f64>;
    type Output = Vec<f64>;
    type Float = f64;

    /// Anneal a parameter vector
    fn anneal(&self, param: &Vec<f64>, temp: f64) -> Result<Vec<f64>, argmin::core::Error> {
        let mut param_n = param.clone();
        let mut rng = self.rng.lock().unwrap();
        // let distr = Uniform::from(0..param.len());
        // Perform modifications to a degree proportional to the current temperature `temp`.
        for _ in 0..(temp.floor() as u64 + 1) {
            // Compute random index of the parameter vector using the supplied random number
            // generator.
            let idx: usize = rng.gen_range(0..param.len());

            // Compute random number in [0.1, 0.1].
            let val = rng.gen_range(-0.5..0.5);

            // modify previous parameter value at random position `idx` by `val`
            param_n[idx] += val;
        }
        Ok(param_n)
    }
}

#[derive(Debug, Clone)]
pub struct LayoutSettings {
    params: LayoutParams,
    positions: Positions,
    iters: LayoutIters,
    init_params: Vec<f64>,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct LayoutIters {
    n_iters: u64,
    temp: f64,
    seed: u64,
}

impl LayoutSettings {
    pub fn new<E, V>(graph: &HedgeGraph<E, V>, params: LayoutParams, iters: LayoutIters) -> Self {
        let (init_params, positions) = Positions::new(graph, iters.seed);

        LayoutSettings {
            params,
            positions,
            iters,
            init_params,
        }
    }

    pub fn left_right_square<E, V>(
        graph: &HedgeGraph<E, V>,
        params: LayoutParams,
        iters: LayoutIters,
        edge: f64,
        left: Vec<Hedge>,
        right: Vec<Hedge>,
    ) -> Self {
        let mut rng = SmallRng::seed_from_u64(iters.seed);
        let mut vertex_positions = Vec::new();
        let mut edge_positions = graph.involution.forgetful_map_node_data_ref(|_| ());

        let range = -edge..edge;

        let mut init_params = Vec::new();

        let left_step = unsafe {
            if left.len() > 1 {
                edge / (left.len().unchecked_sub(1) as f64)
            } else {
                edge
            }
        };
        let mut left_bot_corner = if left.len() <= 1 {
            (-edge / 2., 0.)
        } else {
            (-edge / 2., -edge / 2.)
        };

        let right_step = unsafe {
            if right.len() > 1 {
                edge / (right.len().unchecked_sub(1) as f64)
            } else {
                edge
            }
        };
        let mut right_bot_corner = if right.len() <= 1 {
            (edge / 2., 0.)
        } else {
            (edge / 2., -edge / 2.)
        };

        for i in graph.involution.iter_idx() {
            let j = init_params.len();
            init_params.push(rng.gen_range(range.clone()));
            init_params.push(rng.gen_range(range.clone()));
            edge_positions.set_edge_data(i, (None, j, j + 1));
        }

        for i in left {
            let (_, a, b) = edge_positions.edge_data(i).data.unwrap();
            edge_positions.set_edge_data(i, (Some(left_bot_corner), a, b));
            left_bot_corner.1 += left_step;
        }

        for i in right {
            let (_, a, b) = edge_positions.edge_data(i).data.unwrap();
            edge_positions.set_edge_data(i, (Some(right_bot_corner), a, b));
            right_bot_corner.1 += right_step;
        }

        for _ in graph.base_nodes_iter() {
            let j = init_params.len();
            init_params.push(rng.gen_range(range.clone()));
            init_params.push(rng.gen_range(range.clone()));
            vertex_positions.push((None, j, j + 1));
        }

        LayoutSettings {
            params,
            iters,
            positions: Positions {
                vertex_positions,
                edge_positions,
            },
            init_params,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn circle_ext<E, V>(
        graph: &HedgeGraph<E, V>,
        params: LayoutParams,
        iters: LayoutIters,
        angle_factors: Vec<u32>,
        _n_div: usize,
        shift: Rad<f64>,
        radius: f64,
    ) -> Self {
        let mut rng = SmallRng::seed_from_u64(iters.seed);
        let mut vertex_positions = Vec::new();
        let mut edge_positions = graph.involution.forgetful_map_node_data_ref(|_| ());

        let range = -1.0..1.0;

        let mut init_params = Vec::new();

        let angle_step = Rad(2. * std::f64::consts::PI / f64::from(graph.n_externals() as u32));

        let mut exti = 0;

        for i in graph.involution.iter_idx() {
            let j = init_params.len();
            init_params.push(rng.gen_range(range.clone()));
            init_params.push(rng.gen_range(range.clone()));
            if i == graph.involution.inv(i) {
                let angle = shift + angle_step * (angle_factors[exti] as f64);
                exti += 1;
                let x = radius * angle.cos();
                let y = radius * angle.sin();
                edge_positions.set_edge_data(i, (Some((x, y)), j, j + 1));
            } else {
                edge_positions.set_edge_data(i, (None, j, j + 1));
            }
        }

        for _ in graph.base_nodes_iter() {
            let j = init_params.len();
            init_params.push(rng.gen_range(range.clone()));
            init_params.push(rng.gen_range(range.clone()));
            vertex_positions.push((None, j, j + 1));
        }

        LayoutSettings {
            params,

            iters,

            positions: Positions {
                vertex_positions,
                edge_positions,
            },
            init_params,
        }
    }
}

impl<E, V> HedgeGraph<E, V> {
    pub fn layout(
        self,
        mut settings: LayoutSettings,
    ) -> HedgeGraph<LayoutEdge<E>, LayoutVertex<V>> {
        let layout = GraphLayout {
            rng: Arc::new(Mutex::new(SmallRng::seed_from_u64(0))),
            positions: settings.positions.clone(),
            graph: &self,
            params: settings.params,
        };

        let solver = SimulatedAnnealing::new(settings.iters.temp).unwrap();

        let best = {
            let res = Executor::new(layout, solver)
                .configure(|state| {
                    state
                        .param(settings.init_params)
                        .max_iters(settings.iters.n_iters)
                })
                // run the solver on the defined problem
                .run()
                .unwrap();
            res.state().get_best_param().unwrap().clone()
        };

        let max = settings.positions.max(&best).unwrap();

        let best: Vec<_> = best.into_iter().map(|a| a / max).collect();
        settings.positions.scale(max);

        settings.positions.to_graph(self, &best)
    }
}
