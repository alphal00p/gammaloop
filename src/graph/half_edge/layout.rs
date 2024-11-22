use std::sync::{Arc, Mutex};

// use crate::half_edge::drawing::Decoration;
use argmin::{
    core::{CostFunction, Executor, Gradient, State},
    solver::simulatedannealing::{Anneal, SimulatedAnnealing},
};
use cgmath::{Angle, Rad, Vector2};
use indexmap::IndexMap;
use itertools::Itertools;
use rand::{rngs::SmallRng, Rng, SeedableRng};
use statrs::statistics::Statistics;

use super::{
    drawing::{CetzEdge, Decoration, EdgeGeometry},
    subgraph::node::HedgeNode,
    HedgeGraph, Involution, InvolutiveMapping,
};

pub struct LayoutVertex<V> {
    data: V,
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
    data: E,
    geometry: EdgeGeometry,
}

impl<E> LayoutEdge<E> {
    pub fn new(data: E, x: f64, y: f64) -> Self {
        LayoutEdge {
            data,
            geometry: EdgeGeometry::Simple {
                pos: Vector2::new(x, y),
            },
        }
    }

    pub fn to_fancy(
        &mut self,
        source: Vector2<f64>,
        sink: Option<Vector2<f64>>,
        label: f64, //label shift
        arrow: Option<(f64, f64)>,
    ) {
        match self.geometry {
            EdgeGeometry::Simple { .. } => {}
            _ => {
                return;
            }
        };
        self.geometry = self.geometry.clone().to_fancy(source, sink, label, arrow);
    }
}

impl<E> LayoutEdge<E> {
    pub fn cetz(&self, source: usize, sink: Option<usize>) -> String
    where
        E: CetzEdge,
    {
        self.cetz_impl(&|e| e.label(), &|e| e.decoration(), source, sink)
    }

    pub fn cetz_impl(
        &self,
        label: &impl Fn(&E) -> String,
        decoration: &impl Fn(&E) -> Decoration,
        source: usize,
        sink: Option<usize>,
    ) -> String {
        match &self.geometry {
            EdgeGeometry::Simple { pos } => {
                if let Some(sink) = sink {
                    format!(
                        "edge(node{}.pos,({:.4},{:.4}),node{}.pos,decoration:{})\n",
                        source,
                        pos.x,
                        pos.y,
                        sink,
                        decoration(&self.data).to_cetz()
                    )
                } else {
                    format!(
                        "edge(node{}.pos,({:.4},{:.4}),decoration:{})\n",
                        source,
                        pos.x,
                        pos.y,
                        decoration(&self.data).to_cetz()
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
                    "edge(node{}.pos,({:.4},{:.4}),node{}.pos,decoration:{})\ncontent({:.4},{:.4},angle:{:.2}rad,[{}])\n",
                    source,
                    pos.x,
                    pos.y,
                    sink,
                    decoration(&self.data).to_cetz(),
                    label_pos.x,
                    label_pos.y,
                    label_angle.0,
                    label(&self.data)
                )
                } else {
                    format!(
                    "edge(node{}.pos,({:.4},{:.4}),decoration:{})\ncontent({:.4},{:.4},angle:{:.2}rad,[{}])\n",
                    source,
                    pos.x,
                    pos.y,
                    decoration(&self.data).to_cetz(),
                    label_pos.x,
                    label_pos.y,
                    label_angle.0,
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
                "edge(node{}.pos,({:.4},{:.4}),node{}.pos,decoration:{})\ncontent(({:.4},{:.4}),angle:{:.2}rad,[{}])\n{}\n",
                source,
                pos.x,
                pos.y,
                sink,
                decoration(&self.data).to_cetz(),
                label_pos.x,
                label_pos.y,
                label_angle.0,
                label(&self.data),
                arrow_arc.hobby_arrow(),
            )
                } else {
                    format!(
                "edge(node{}.pos,({:.4},{:.4}),decoration:{})\ncontent(({:.4},{:.4}),angle:{:.2}rad,[{}])\n{}\n",
                source,
                pos.x,
                pos.y,
                decoration(&self.data).to_cetz(),
                label_pos.x,
                label_pos.y,
                label_angle.0,
                label(&self.data),
                arrow_arc.hobby_arrow(),
                )
                }
            }
        }
    }
}

pub type PositionalHedgeGraph<E, V> = HedgeGraph<LayoutEdge<E>, LayoutVertex<V>>;

impl<E, V> PositionalHedgeGraph<E, V> {
    pub fn to_fancy(&mut self) {
        for (eid, nid) in self.involution.hedge_data.iter().enumerate() {
            let i = self.nodes.get(nid).unwrap().pos;

            if let InvolutiveMapping::Identity(d) = &mut self.involution.inv[eid] {
                if let Some(d) = &mut d.data {
                    d.to_fancy(i, None, 0.2, Some((0.2, 0.4)));
                }
            }

            if let Some(cn) = self.involution.get_connected_node_id(eid) {
                let j = self.nodes.get(cn).unwrap().pos;

                match &mut self.involution.inv[eid] {
                    InvolutiveMapping::Source((d, _)) => {
                        if let Some(d) = &mut d.data {
                            d.to_fancy(i, Some(j), 0.2, Some((0.2, 0.4)));
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    fn cetz_preamble() -> String {
        let mut out = String::new();

        out.push_str(
            "#import \"@preview/cetz:0.3.0\"\n
            #set page(width: auto,height:auto)\n
          ",
        );
        out
    }

    pub fn cetz_impl_collection(
        graphs: &[Self],
        col_name: &str,
        edge_label: &impl Fn(&E) -> String,
        edge_decoration: &impl Fn(&E) -> Decoration,
    ) -> String {
        let mut out = Self::cetz_preamble();
        for (i, g) in graphs.iter().enumerate() {
            out.push_str(&format!(
                "#let {col_name}{i}={}",
                g.cetz_bare(edge_label, edge_decoration)
            ))
        }
        out.push_str(
            "#grid(
          columns: 3,
          gutter: 20pt,
        ",
        );
        for i in 0..graphs.len() {
            out.push_str(&format!("{col_name}{i},"));
        }
        out.push_str(")/n");
        out
    }

    fn cetz_bare(
        &self,
        edge_label: &impl Fn(&E) -> String,
        edge_decoration: &impl Fn(&E) -> Decoration,
    ) -> String {
        let mut out = String::from(
            "cetz.canvas({ \n
          import cetz.draw: * \n
          let node(pos)=circle(pos,radius:0.12,fill: black)\n
            let stroke = 1.2pt\n
          let edge(..points,decoration:\"\")={\n
            if decoration == \"coil\"{\n
              cetz.decorations.coil(hobby(..points),amplitude:0.23,segment-length:0.2,stroke:stroke,align:\"MID\")\n
            } else if decoration == \"wave\" {\n
               cetz.decorations.wave(hobby(..points),amplitude:0.24,segment-length:0.3,stroke:stroke)\n
            } else {\n
               hobby(..points,stroke:stroke)\n
            }\n
          }\n",
        );
        for a in 0..self.base_nodes {
            out.push_str(&format!(
                "let node{}= (pos:({:.6},{:.6}))\n",
                a, self.nodes[a].pos.x, self.nodes[a].pos.y
            ));
            out.push_str(&format!("node(node{}.pos)\n", a));
        }

        for (eid, (nid, e)) in self
            .involution
            .hedge_data
            .iter()
            .zip_eq(self.involution.inv.iter())
            .enumerate()
        {
            let i = self.get_node_pos(&nid);

            if let Some(cn) = self.involution.get_connected_node_id(eid) {
                let j = self.get_node_pos(&cn);

                match e {
                    InvolutiveMapping::Source((d, _)) => {
                        if let Some(d) = &d.data {
                            out.push_str(&d.cetz_impl(edge_label, edge_decoration, i, Some(j)));
                        }
                    }
                    _ => {}
                }
            }

            if let InvolutiveMapping::Identity(d) = e {
                if let Some(d) = &d.data {
                    println!("{:?}", d.geometry);
                    out.push_str(&d.cetz_impl(edge_label, edge_decoration, i, None));
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
        let mut out = String::new();

        out.push_str(
            "#import \"@preview/cetz:0.3.0\"\n
        #set page(width: auto,height:auto)\n
        #cetz.canvas({\n
          import cetz.draw: *\n",
        );

        out.push_str(
            "let node(pos)=circle(pos,radius:0.12,fill: black)\n
              let stroke = 1.2pt\n
            let edge(..points,decoration:\"\")={\n
              if decoration == \"coil\"{\n
                cetz.decorations.coil(hobby(..points),amplitude:0.23,segment-length:0.2,stroke:stroke,align:\"MID\")\n
              } else if decoration == \"wave\" {\n
                 cetz.decorations.wave(hobby(..points),amplitude:0.24,segment-length:0.3,stroke:stroke)\n
              } else {\n
                 hobby(..points,stroke:stroke)\n
              }\n
            }\n",
        );

        for a in 0..self.base_nodes {
            out.push_str(&format!(
                "let node{}= (pos:({:.6},{:.6}))\n",
                a, self.nodes[a].pos.x, self.nodes[a].pos.y
            ));
            out.push_str(&format!("node(node{}.pos)\n", a));
        }

        for (eid, (nid, e)) in self
            .involution
            .hedge_data
            .iter()
            .zip_eq(self.involution.inv.iter())
            .enumerate()
        {
            let i = self.get_node_pos(&nid);

            if let Some(cn) = self.involution.get_connected_node_id(eid) {
                let j = self.get_node_pos(&cn);

                match e {
                    InvolutiveMapping::Source((d, _)) => {
                        if let Some(d) = &d.data {
                            out.push_str(&d.cetz_impl(edge_label, edge_decoration, i, Some(j)));
                        }
                    }
                    _ => {}
                }
            }

            if let InvolutiveMapping::Identity(d) = e {
                if let Some(d) = &d.data {
                    println!("{:?}", d.geometry);
                    out.push_str(&d.cetz_impl(edge_label, edge_decoration, i, None));
                }
            }
        }

        out.push_str("})\n");
        out
    }
}

pub struct GraphLayout<'a, E, V> {
    pub graph: &'a HedgeGraph<E, V>,
    pub positions: Positions,
    pub params: LayoutParams,
    pub rng: Arc<Mutex<SmallRng>>,
}

pub struct LayoutParams {
    pub spring_constant: f64,
    pub spring_length: f64,
    pub charge_constant_e: f64,
    pub charge_constant_v: f64,
    pub external_constant: f64,
    pub central_force_constant: f64,
}

impl Default for LayoutParams {
    fn default() -> Self {
        LayoutParams {
            spring_length: 0.2,
            spring_constant: 0.2,
            charge_constant_e: 0.6,
            charge_constant_v: 0.5,
            external_constant: 0.0,
            central_force_constant: 0.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Positions {
    vertex_positions: IndexMap<HedgeNode, (Option<(f64, f64)>, usize, usize)>,
    edge_positions: Involution<(), (Option<(f64, f64)>, usize, usize)>,
}

impl Positions {
    pub fn to_graph<E, V>(
        self,
        graph: HedgeGraph<E, V>,
        params: &[f64],
    ) -> PositionalHedgeGraph<E, V> {
        // let mut new_graph = graph;
        let nodes: IndexMap<HedgeNode, LayoutVertex<V>> = graph
            .nodes
            .into_iter()
            .map(|(e, v)| {
                let pos = self.vertex_positions.get(&e).unwrap();
                (e, LayoutVertex::new(v, params[pos.1], params[pos.2]))
            })
            .collect();

        let inv: Vec<_> = graph
            .involution
            .inv
            .into_iter()
            .enumerate()
            .map(|(i, m)| {
                let pos = &self.get_edge_position(i, params);

                let m = match m {
                    InvolutiveMapping::Sink(s) => InvolutiveMapping::Sink(s),
                    InvolutiveMapping::Source((s, i)) => {
                        InvolutiveMapping::Source((s.map(|s| LayoutEdge::new(s, pos.0, pos.1)), i))
                    }
                    InvolutiveMapping::Identity(s) => {
                        InvolutiveMapping::Identity(s.map(|s| LayoutEdge::new(s, pos.0, pos.1)))
                    }
                };

                m
            })
            .collect();

        HedgeGraph {
            nodes,
            involution: Involution {
                inv,
                hedge_data: graph.involution.hedge_data,
            },
            base_nodes: graph.base_nodes,
        }
    }
    pub fn new<E, V>(graph: &HedgeGraph<E, V>, seed: u64) -> (Vec<f64>, Self) {
        let mut rng = SmallRng::seed_from_u64(seed);
        let mut vertex_positions = IndexMap::new();
        let mut edge_positions = graph.involution.forgetful_map_node_data_ref(|_| ());

        let range = -1.0..1.0;

        let mut params = Vec::new();

        for i in 0..graph.involution.inv.len() {
            let j = params.len();
            params.push(rng.gen_range(range.clone()));
            params.push(rng.gen_range(range.clone()));
            edge_positions.set_data(i, (None, j, j + 1));
        }

        for (node, data) in graph.nodes.iter() {
            let j = params.len();
            params.push(rng.gen_range(range.clone()));
            params.push(rng.gen_range(range.clone()));
            vertex_positions.insert(node.clone(), (None, j, j + 1));
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
        let mut vertex_positions = IndexMap::new();
        let mut edge_positions = graph.involution.forgetful_map_node_data_ref(|_| ());

        let range = -1.0..1.0;

        let mut params = Vec::new();

        let mut angle = 0.;
        let angle_step = 2. * std::f64::consts::PI / f64::from(graph.n_externals() as u32);
        let mut ext = 0.;

        for i in 0..graph.involution.inv.len() {
            let j = params.len();
            params.push(rng.gen_range(range.clone()));
            params.push(rng.gen_range(range.clone()));
            if i == graph.involution.inv(i) {
                let x = radius * f64::cos(angle);
                let y = radius * f64::sin(angle);
                edge_positions.set_data(i, (Some((x, y)), j, j + 1));
                angle += angle_step;
            } else {
                edge_positions.set_data(i, (None, j, j + 1));
            }
        }

        for (node, data) in graph.nodes.iter() {
            let j = params.len();
            params.push(rng.gen_range(range.clone()));
            params.push(rng.gen_range(range.clone()));
            vertex_positions.insert(node.clone(), (None, j, j + 1));
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
    ) -> impl Iterator<Item = (&'a HedgeNode, (f64, f64))> + 'a {
        self.vertex_positions.iter().map(|(node, (p, x, y))| {
            if let Some(p) = p {
                (node, (p.0, p.1))
            } else {
                (node, (params[*x], params[*y]))
            }
        })
    }

    pub fn get_edge_position(&self, edge: usize, params: &[f64]) -> (f64, f64) {
        let (p, ix, iy) = self.edge_positions.get_data(edge).data.unwrap();
        if let Some(p) = p {
            (p.0, p.1)
        } else {
            (params[*ix], params[*iy])
        }
    }
}

impl<'a, E, V> CostFunction for GraphLayout<'a, E, V> {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        let mut cost = 0.0;

        for (node, (x, y)) in self.positions.iter_vertex_positions(param) {
            for e in node.hairs.iter_ones() {
                let (ex, ey) = self.positions.get_edge_position(e, param);

                let dx = x - ex;
                let dy = y - ey;

                let dist = (dx * dx + dy * dy).sqrt();

                cost +=
                    0.5 * self.params.spring_constant * (self.params.spring_length - dist).powi(2);
                cost += self.params.charge_constant_e / dist;

                for othere in node.hairs.iter_ones() {
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

impl<'a, E, V> Anneal for GraphLayout<'a, E, V> {
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

pub struct LayoutSettings {
    params: LayoutParams,
    temperature: f64,
    // iters: usize,
    positions: Positions,
    // seed: u64,
    init_params: Vec<f64>,
}

impl LayoutSettings {
    pub fn new<E, V>(
        graph: &HedgeGraph<E, V>,
        params: LayoutParams,
        seed: u64,
        temperature: f64,
    ) -> Self {
        let (init_params, positions) = Positions::new(graph, seed);

        LayoutSettings {
            params,
            temperature,
            positions,
            init_params,
        }
    }

    pub fn left_right_square<E, V>(
        graph: &HedgeGraph<E, V>,
        params: LayoutParams,
        seed: u64,
        temperature: f64,
        edge: f64,
        left: Vec<usize>,
        right: Vec<usize>,
    ) -> Self {
        let mut rng = SmallRng::seed_from_u64(seed);
        let mut vertex_positions = IndexMap::new();
        let mut edge_positions = graph.involution.forgetful_map_node_data_ref(|_| ());

        let range = -edge..edge;

        let mut init_params = Vec::new();

        let left_step = edge / (left.len() as f64);
        let mut left_bot_corner = (-edge / 2., -edge / 2.);
        let right_step = edge / (right.len() as f64);
        let mut right_bot_corner = (edge / 2., -edge / 2.);

        for i in 0..graph.involution.inv.len() {
            let j = init_params.len();
            init_params.push(rng.gen_range(range.clone()));
            init_params.push(rng.gen_range(range.clone()));
            edge_positions.set_data(i, (None, j, j + 1));
        }

        for i in left {
            let (_, a, b) = *edge_positions.get_data(i).data.unwrap();
            edge_positions.set_data(i, (Some(left_bot_corner), a, b));
            left_bot_corner.1 += left_step;
        }

        for i in right {
            let (_, a, b) = *edge_positions.get_data(i).data.unwrap();
            edge_positions.set_data(i, (Some(right_bot_corner), a, b));
            right_bot_corner.1 += right_step;
        }

        for (node, data) in graph.nodes.iter() {
            let j = init_params.len();
            init_params.push(rng.gen_range(range.clone()));
            init_params.push(rng.gen_range(range.clone()));
            vertex_positions.insert(node.clone(), (None, j, j + 1));
        }

        LayoutSettings {
            params,
            temperature,
            positions: Positions {
                vertex_positions,
                edge_positions,
            },
            init_params,
        }
    }

    pub fn circle_ext<E, V>(
        graph: &HedgeGraph<E, V>,
        params: LayoutParams,
        seed: u64,
        temperature: f64,
        angle_factors: Vec<u32>,
        n_div: usize,
        shift: Rad<f64>,
        radius: f64,
    ) -> Self {
        let mut rng = SmallRng::seed_from_u64(seed);
        let mut vertex_positions = IndexMap::new();
        let mut edge_positions = graph.involution.forgetful_map_node_data_ref(|_| ());

        let range = -1.0..1.0;

        let mut init_params = Vec::new();

        let angle_step = Rad(2. * std::f64::consts::PI / f64::from(graph.n_externals() as u32));

        let mut exti = 0;

        for i in 0..graph.involution.inv.len() {
            let j = init_params.len();
            init_params.push(rng.gen_range(range.clone()));
            init_params.push(rng.gen_range(range.clone()));
            if i == graph.involution.inv(i) {
                let angle = shift + angle_step * (angle_factors[exti] as f64);
                exti += 1;
                let x = radius * angle.cos();
                let y = radius * angle.sin();
                edge_positions.set_data(i, (Some((x, y)), j, j + 1));
            } else {
                edge_positions.set_data(i, (None, j, j + 1));
            }
        }

        for (node, data) in graph.nodes.iter() {
            let j = init_params.len();
            init_params.push(rng.gen_range(range.clone()));
            init_params.push(rng.gen_range(range.clone()));
            vertex_positions.insert(node.clone(), (None, j, j + 1));
        }

        LayoutSettings {
            params,
            temperature,
            positions: Positions {
                vertex_positions,
                edge_positions,
            },
            init_params,
        }
    }
}

impl<E, V> HedgeGraph<E, V> {
    pub fn layout(self, settings: LayoutSettings) -> HedgeGraph<LayoutEdge<E>, LayoutVertex<V>> {
        let layout = GraphLayout {
            rng: Arc::new(Mutex::new(SmallRng::seed_from_u64(0))),
            positions: settings.positions.clone(),
            graph: &self,
            params: settings.params,
        };

        let solver = SimulatedAnnealing::new(settings.temperature).unwrap();

        let best = {
            let res = Executor::new(layout, solver)
                .configure(|state| state.param(settings.init_params).max_iters(10000))
                // run the solver on the defined problem
                .run()
                .unwrap();
            res.state().get_best_param().unwrap().clone()
        };

        let max = best.iter().map(|a| a.abs()).reduce(f64::min).unwrap();

        let best: Vec<_> = best.into_iter().map(|a| a / max).collect();

        settings.positions.to_graph(self, &best)
    }
}
