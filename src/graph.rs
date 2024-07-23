use crate::{
    cff::{
        esurface::{
            generate_esurface_data, get_existing_esurfaces, EsurfaceDerivedData, ExistingEsurfaces,
            ExternalShift,
        },
        expression::CFFExpression,
        generation::generate_cff_expression,
    },
    gammaloop_integrand::DefaultSample,
    ltd::{generate_ltd_expression, LTDExpression, SerializableLTDExpression},
    model::{self, EdgeSlots, Model, VertexSlots},
    momentum::{FourMomentum, ThreeMomentum},
    numerator::Numerator,
    subtraction::{
        overlap::{find_maximal_overlap, OverlapStructure},
        static_counterterm,
    },
    tropical::{self, TropicalSubgraphTable},
    utils::{
        compute_four_momentum_from_three, compute_three_momentum_from_four, sorted_vectorize,
        FloatLike, F,
    },
};

use ahash::RandomState;

use color_eyre::{Help, Report};
use enum_dispatch::enum_dispatch;
use eyre::eyre;
use itertools::Itertools;
use log::{debug, warn};
use nalgebra::DMatrix;
#[allow(unused_imports)]
use spenso::Contract;
use spenso::{
    ufo::{preprocess_ufo_color_wrapped, preprocess_ufo_spin_wrapped},
    Complex, *,
};

use core::panic;
use serde::{Deserialize, Serialize};
use smallvec::{smallvec, SmallVec};
use smartstring::{LazyCompact, SmartString};
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
    sync::Arc,
};

use symbolica::{
    atom::Atom,
    domains::float::NumericalFloatLike,
    id::{Pattern, Replacement},
};
//use symbolica::{atom::Symbol,state::State};

use constcat::concat;

const MAX_COLOR_INNER_CONTRACTIONS: usize = 3;

#[derive(Debug, Default, Clone, Serialize, Deserialize, PartialEq)]
pub enum EdgeType {
    #[serde(rename = "in")]
    Incoming,
    #[serde(rename = "out")]
    Outgoing,
    #[default]
    #[serde(rename = "virtual")]
    Virtual,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[allow(non_snake_case)]
#[serde(tag = "type")]
pub enum SerializableVertexInfo {
    #[serde(rename = "external_vertex_info")]
    ExternalVertexInfo(SerializableExternalVertexInfo),
    #[serde(rename = "interacton_vertex_info")]
    InteractonVertexInfo(SerializableInteractionVertexInfo),
}

#[enum_dispatch]
pub trait HasVertexInfo {
    fn get_type(&self) -> SmartString<LazyCompact>;
}

#[derive(Debug, Clone)]
#[enum_dispatch(HasVertexInfo)]
pub enum VertexInfo {
    ExternalVertexInfo(ExternalVertexInfo),
    InteractonVertexInfo(InteractionVertexInfo),
}

impl VertexInfo {
    pub fn generate_vertex_slots(
        &self,
        shifts: (usize, usize, usize),
    ) -> (VertexSlots, (usize, usize, usize)) {
        match self {
            VertexInfo::ExternalVertexInfo(e) => {
                let (e, shifts) = e.particle.slots(shifts);
                (e.into(), shifts)
            }
            VertexInfo::InteractonVertexInfo(i) => i.vertex_rule.generate_vertex_slots(shifts),
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct SerializableInteractionVertexInfo {
    vertex_rule: SmartString<LazyCompact>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct SerializableExternalVertexInfo {
    direction: EdgeType,
    particle: SmartString<LazyCompact>,
}

#[derive(Debug, Clone)]
pub struct ExternalVertexInfo {
    direction: EdgeType,
    particle: Arc<model::Particle>,
}

impl HasVertexInfo for ExternalVertexInfo {
    fn get_type(&self) -> SmartString<LazyCompact> {
        SmartString::<LazyCompact>::from("external_vertex_info")
    }
}

#[derive(Debug, Clone)]
pub struct InteractionVertexInfo {
    #[allow(unused)]
    pub vertex_rule: Arc<model::VertexRule>,
}

impl InteractionVertexInfo {
    pub fn apply_vertex_rule(
        &self,
        edges: &[usize],
        vertex_pos: usize,
        graph: &Graph,
    ) -> Option<[DataTensor<Atom>; 3]> {
        let vertex_slots = &graph.derived_data.vertex_slots.as_ref().unwrap()[vertex_pos];

        let spin_structure = self
            .vertex_rule
            .lorentz_structures
            .iter()
            .map(|ls| {
                let mut atom = ls.structure.clone();

                for (i, e) in edges.iter().enumerate() {
                    let momentum_in_pattern = Pattern::parse(&format!("P(x_,{})", i + 1)).unwrap();

                    let momentum_out_pattern =
                        Pattern::parse(&format!("Q({},aind(lor(4,indexid(x_))))", e)).unwrap();

                    atom = momentum_in_pattern.replace_all(
                        atom.as_view(),
                        &momentum_out_pattern,
                        None,
                        None,
                    );
                }

                atom = preprocess_ufo_spin_wrapped(atom);

                for (i, _) in edges.iter().enumerate() {
                    let replacements = vertex_slots[i].replacements(i + 1);

                    let reps: Vec<Replacement> = replacements
                        .iter()
                        .map(|(pat, rhs)| Replacement::new(pat, rhs))
                        .collect();

                    atom = atom.replace_all_multiple(&reps);
                }

                atom
            })
            .collect_vec();

        let color_structure: Vec<Atom> = self
            .vertex_rule
            .color_structures
            .iter()
            .map(|cs| {
                let mut atom = cs.clone();
                //a is adjoint index, i is fundamental index, ia is antifundamental index, s is sextet ,sa is antisextet

                atom = preprocess_ufo_color_wrapped(atom);
                //T(1,2,3) Fundamental representation matrix (T a1 )ı  ̄3 i2
                // f(1,2,3) Antisymmetric structure constant f a1a2a3
                // d(1,2,3) Symmetric structure constant da1 a2 a3
                // Epsilon(1,2,3) Fundamental Levi-Civita tensor εi1 i2 i3 EpsilonBar(1,2,3) Antifundamental Levi-Civita tensor εı  ̄1 ı  ̄2 ı  ̄3
                // T6(1,2,3) Sextet representation matrix (T a1 6 ) ̄ α3 α2
                // K6(1,2,3) Sextet Clebsch-Gordan coefficient (K6)ı  ̄2 ı  ̄3 α1 K6Bar(1,2,3) Antisextet Clebsch-Gordan coefficient (K6) ̄ α1 i2i3

                // First process kronkers, with respect to spin:

                let spins: Vec<isize> =
                    self.vertex_rule.particles.iter().map(|s| s.color).collect();

                for (i, s) in spins.iter().enumerate() {
                    let id1 = Pattern::parse(&format!("Identity({},x_)", i + 1)).unwrap();

                    let id2 = Pattern::parse(&format!("id(aind(x_,{}))", i + 1)).unwrap();

                    let ind = match s {
                        1 => concat!(EUCLIDEAN, "(1,"),
                        3 => concat!(COLORFUND, "(3,"),
                        -3 => concat!(COLORANTIFUND, "(3,"),
                        6 => concat!(COLORSEXT, "(6,"),
                        -6 => concat!(COLORANTISEXT, "(6,"),
                        8 => concat!(COLORADJ, "(8,"),
                        i => panic!("Color {i}not supported "),
                    };

                    atom = id1.replace_all(
                        atom.as_view(),
                        &Pattern::parse(&format!("id(aind({}indexid({})),x_))", ind, i + 1))
                            .unwrap(),
                        None,
                        None,
                    );

                    atom = id2.replace_all(
                        atom.as_view(),
                        &Pattern::parse(&format!("id(aind(x_,{}indexid({}))))", ind, i + 1))
                            .unwrap(),
                        None,
                        None,
                    );
                }

                for (i, _) in edges.iter().enumerate() {
                    let replacements = vertex_slots[i].replacements(i + 1);

                    let reps: Vec<Replacement> = replacements
                        .iter()
                        .map(|(pat, rhs)| Replacement::new(pat, rhs))
                        .collect();

                    atom = atom.replace_all_multiple(&reps);
                }

                for i in 0..MAX_COLOR_INNER_CONTRACTIONS {
                    let pat: Pattern = Atom::parse(&format!("indexid({})", -1 - i as i64))
                        .unwrap()
                        .into_pattern();

                    atom = pat.replace_all(
                        atom.as_view(),
                        &Atom::new_num(i as i64).into_pattern(),
                        None,
                        None,
                    );
                }

                atom
            })
            .collect();

        let i = Slot {
            index: AbstractIndex::try_from(format!("i{}", vertex_pos)).unwrap(),
            representation: Representation::Euclidean(color_structure.len().into()),
        };

        let color_structure = DataTensor::Dense(
            DenseTensor::from_data(&color_structure, VecStructure::from(vec![i])).unwrap(),
        );

        let j = Slot {
            index: AbstractIndex::try_from(format!("j{}", vertex_pos)).unwrap(),
            representation: Representation::Euclidean(spin_structure.len().into()),
        };
        let spin_structure = DataTensor::Dense(
            DenseTensor::from_data(&spin_structure, VecStructure::from(vec![j])).unwrap(),
        );

        let mut couplings: DataTensor<Atom> =
            DataTensor::Sparse(SparseTensor::empty(VecStructure::from(vec![i, j])));

        for (i, row) in self.vertex_rule.couplings.iter().enumerate() {
            for (j, col) in row.iter().enumerate() {
                if let Some(atom) = col {
                    couplings.set(&[i, j], atom.expression.clone()).unwrap();
                }
            }
        }

        Some([spin_structure, color_structure, couplings])
    }
}

impl HasVertexInfo for InteractionVertexInfo {
    fn get_type(&self) -> SmartString<LazyCompact> {
        SmartString::<LazyCompact>::from("interacton_vertex_info")
    }
}

fn serializable_vertex_info(vertex_info: &VertexInfo) -> SerializableVertexInfo {
    match vertex_info {
        VertexInfo::ExternalVertexInfo(external_vertex_info) => {
            SerializableVertexInfo::ExternalVertexInfo(SerializableExternalVertexInfo {
                direction: external_vertex_info.direction.clone(),
                particle: external_vertex_info.particle.name.clone(),
            })
        }
        VertexInfo::InteractonVertexInfo(interaction_vertex_info) => {
            SerializableVertexInfo::InteractonVertexInfo(SerializableInteractionVertexInfo {
                vertex_rule: interaction_vertex_info.vertex_rule.name.clone(),
            })
        }
    }
}

fn deserialize_vertex_info(
    model: &model::Model,
    serialized_vertex_info: &SerializableVertexInfo,
) -> VertexInfo {
    match serialized_vertex_info {
        SerializableVertexInfo::ExternalVertexInfo(external_vertex_info) => {
            VertexInfo::ExternalVertexInfo(
                ExternalVertexInfo {
                    direction: external_vertex_info.direction.clone(),
                    particle: model.get_particle(&external_vertex_info.particle),
                }
                .clone(),
            )
        }
        SerializableVertexInfo::InteractonVertexInfo(interaction_vertex_info) => {
            VertexInfo::InteractonVertexInfo(InteractionVertexInfo {
                vertex_rule: model
                    .get_vertex_rule(&interaction_vertex_info.vertex_rule)
                    .clone(),
            })
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct SerializableEdge {
    name: SmartString<LazyCompact>,
    edge_type: EdgeType,
    particle: SmartString<LazyCompact>,
    propagator: SmartString<LazyCompact>,
    vertices: [SmartString<LazyCompact>; 2],
}

impl SerializableEdge {
    pub fn from_edge(graph: &Graph, edge: &Edge) -> SerializableEdge {
        SerializableEdge {
            name: edge.name.clone(),
            edge_type: edge.edge_type.clone(),
            particle: edge.particle.name.clone(),
            propagator: edge.propagator.name.clone(),
            vertices: [
                graph.vertices[edge.vertices[0]].name.clone(),
                graph.vertices[edge.vertices[1]].name.clone(),
            ],
        }
    }
}

#[derive(Debug, Clone)]
pub struct Edge {
    pub name: SmartString<LazyCompact>,
    pub edge_type: EdgeType,
    pub propagator: Arc<model::Propagator>,
    pub particle: Arc<model::Particle>,
    pub vertices: [usize; 2],
}

impl Edge {
    pub fn from_serializable_edge(
        model: &model::Model,
        graph: &Graph,
        serializable_edge: &SerializableEdge,
    ) -> Edge {
        Edge {
            name: serializable_edge.name.clone(),
            edge_type: serializable_edge.edge_type.clone(),
            particle: model.get_particle(&serializable_edge.particle),
            propagator: model.get_propagator(&serializable_edge.propagator),
            vertices: [
                graph
                    .get_vertex_position(&serializable_edge.vertices[0])
                    .unwrap(),
                graph
                    .get_vertex_position(&serializable_edge.vertices[1])
                    .unwrap(),
            ],
        }
    }

    pub fn is_incoming_to(&self, vertex: usize) -> bool {
        self.vertices[1] == vertex
    }

    pub fn denominator(&self, graph: &Graph) -> (Atom, Atom) {
        let num = *graph.edge_name_to_position.get(&self.name).unwrap();
        let mom = Atom::parse(&format!("Q{num}")).unwrap();
        let mass = self
            .particle
            .mass
            .expression
            .clone()
            .unwrap_or(Atom::new_num(0));

        (mom, mass)
    }

    pub fn substitute_lmb(&self, atom: Atom, graph: &Graph, lmb: &LoopMomentumBasis) -> Atom {
        let num = *graph.edge_name_to_position.get(&self.name).unwrap();
        let mom = Pattern::parse(&format!("Q({num},x_)")).unwrap();
        let mom_rep = lmb.pattern(num);
        atom.replace_all(&mom, &mom_rep, None, None)
    }

    // pub fn edge_momentum_symbol(&self, graph: &Graph) -> Symbol {
    //     let num = *graph.edge_name_to_position.get(&self.name).unwrap();
    //     State::get_symbol(format!("Q{num}"))
    // }

    pub fn in_slot<'a>(&self, graph: &'a Graph) -> &'a EdgeSlots {
        let local_pos_in_sink_vertex =
            graph.vertices[self.vertices[1]].get_local_edge_position(self, graph);

        &graph.derived_data.vertex_slots.as_ref().unwrap()[self.vertices[1]]
            [local_pos_in_sink_vertex]
    }

    pub fn out_slot<'a>(&self, graph: &'a Graph) -> &'a EdgeSlots {
        let local_pos_in_sink_vertex =
            graph.vertices[self.vertices[0]].get_local_edge_position(self, graph);

        &graph.derived_data.vertex_slots.as_ref().unwrap()[self.vertices[0]]
            [local_pos_in_sink_vertex]
    }

    pub fn numerator(&self, graph: &Graph) -> (Atom, usize) {
        let num = *graph.edge_name_to_position.get(&self.name).unwrap();
        let in_slots = self.in_slot(graph);
        let out_slots = self.out_slot(graph);

        match self.edge_type {
            EdgeType::Incoming => (self.particle.incoming_polarization_atom(in_slots, num), 0),
            EdgeType::Outgoing => (self.particle.outgoing_polarization_atom(out_slots, num), 0),
            EdgeType::Virtual => {
                let mut atom = self.propagator.numerator.clone();

                let pfun = Pattern::parse("P(x_)").unwrap();
                atom = pfun.replace_all(
                    atom.as_view(),
                    &Pattern::parse(&format!("Q({},aind(lor(4,x_)))", num)).unwrap(),
                    None,
                    None,
                );

                let pslashfun = Pattern::parse("PSlash(i_,j_)").unwrap();
                let pindex_num = graph.shifts.0 + 1;
                atom = pslashfun.replace_all(
                    atom.as_view(),
                    &Pattern::parse(&format!(
                        "Q({},aind(lor(4,{})))*Gamma({},i_,j_)",
                        num, pindex_num, pindex_num
                    ))
                    .unwrap(),
                    None,
                    None,
                );

                atom = preprocess_ufo_spin_wrapped(atom);

                let replacements_in = in_slots.replacements(1);

                let mut replacements_out = out_slots.replacements(2);

                replacements_out.push((
                    Atom::parse("indexid(x_)").unwrap().into_pattern(),
                    Atom::parse("x_").unwrap().into_pattern(),
                ));

                for (&cin, &cout) in in_slots.color.iter().zip(out_slots.color.iter()) {
                    let id: NamedStructure<&str, ()> =
                        NamedStructure::from_iter([cin, cout], "id", None);
                    atom = atom * &id.to_symbolic().unwrap();
                }

                let reps: Vec<Replacement> = replacements_in
                    .iter()
                    .chain(replacements_out.iter())
                    .map(|(pat, rhs)| Replacement::new(pat, rhs))
                    .collect();

                (atom.replace_all_multiple(&reps), pindex_num + 1)
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableVertex {
    name: SmartString<LazyCompact>,
    vertex_info: SerializableVertexInfo,
    edges: Vec<SmartString<LazyCompact>>,
}

impl SerializableVertex {
    pub fn from_vertex(graph: &Graph, vertex: &Vertex) -> SerializableVertex {
        SerializableVertex {
            name: vertex.name.clone(),
            vertex_info: serializable_vertex_info(&vertex.vertex_info),
            edges: vertex
                .edges
                .iter()
                .map(|e| graph.edges[*e].name.clone())
                .collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Vertex {
    pub name: SmartString<LazyCompact>,
    pub vertex_info: VertexInfo,
    pub edges: Vec<usize>,
}

impl Vertex {
    pub fn get_local_edge_position(&self, edge: &Edge, graph: &Graph) -> usize {
        let global_id: usize = graph.edge_name_to_position[&edge.name];
        self.edges
            .iter()
            .enumerate()
            .find(|(_, &e)| e == global_id)
            .unwrap()
            .0
    }

    pub fn generate_vertex_slots(
        &self,
        shifts: (usize, usize, usize),
    ) -> (VertexSlots, (usize, usize, usize)) {
        self.vertex_info.generate_vertex_slots(shifts)
    }

    pub fn from_serializable_vertex(model: &model::Model, vertex: &SerializableVertex) -> Vertex {
        Vertex {
            name: vertex.name.clone(),
            vertex_info: deserialize_vertex_info(model, &vertex.vertex_info),
            // This will be filled in later during deserialization of the complete graph
            edges: vec![],
        }
    }

    pub fn is_edge_incoming(&self, edge: usize, graph: &Graph) -> bool {
        graph.is_edge_incoming(edge, graph.get_vertex_position(&self.name).unwrap())
    }

    pub fn apply_vertex_rule(&self, graph: &Graph) -> Option<[DataTensor<Atom>; 3]> {
        match &self.vertex_info {
            VertexInfo::ExternalVertexInfo(_) => None,
            VertexInfo::InteractonVertexInfo(interaction_vertex_info) => interaction_vertex_info
                .apply_vertex_rule(
                    &self.edges,
                    graph.get_vertex_position(&self.name).unwrap(),
                    graph,
                ),
        }
    }

    pub fn contracted_vertex_rule(&self, graph: &Graph) -> Option<Atom> {
        let all = self.apply_vertex_rule(graph)?;
        let scalar = all
            .into_iter()
            .reduce(|acc, tensor| acc.contract(&tensor).unwrap())
            .unwrap()
            .get(&[])
            .unwrap()
            .clone();

        Some(scalar)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableGraph {
    name: SmartString<LazyCompact>,
    vertices: Vec<SerializableVertex>,
    edges: Vec<SerializableEdge>,
    overall_factor: f64,
    #[allow(clippy::type_complexity)]
    external_connections: Vec<(
        Option<SmartString<LazyCompact>>,
        Option<SmartString<LazyCompact>>,
    )>,
    loop_momentum_basis: Vec<SmartString<LazyCompact>>,
    #[serde(with = "sorted_vectorize")]
    edge_signatures: HashMap<SmartString<LazyCompact>, (Vec<isize>, Vec<isize>)>,
}

impl SerializableGraph {
    pub fn from_graph(graph: &Graph) -> SerializableGraph {
        SerializableGraph {
            name: graph.name.clone(),
            vertices: graph
                .vertices
                .iter()
                .map(|v| SerializableVertex::from_vertex(graph, v))
                .collect(),
            edges: graph
                .edges
                .iter()
                .map(|e| SerializableEdge::from_edge(graph, e))
                .collect(),

            overall_factor: graph.overall_factor,
            external_connections: graph
                .external_connections
                .iter()
                .map(|(v1, v2)| {
                    (
                        v1.as_ref().map(|&v| graph.vertices[v].name.clone()),
                        v2.as_ref().map(|&v| graph.vertices[v].name.clone()),
                    )
                })
                .collect(),
            loop_momentum_basis: graph
                .loop_momentum_basis
                .basis
                .iter()
                .map(|&e| graph.edges[e].name.clone())
                .collect(),
            edge_signatures: graph
                .loop_momentum_basis
                .edge_signatures
                .iter()
                .enumerate()
                .map(|(i_e, sig)| (graph.edges[i_e].name.clone(), sig.clone()))
                .collect(),
        }
    }

    pub fn from_yaml_str(yaml_str: String) -> Result<SerializableGraph, Report> {
        serde_yaml::from_str(yaml_str.as_str())
            .map_err(|e| eyre!(format!("Error parsing graph yaml: {}", e)))
            .suggestion("Is it a correct yaml file")
    }
}

#[derive(Debug, Clone)]
pub struct Graph {
    pub name: SmartString<LazyCompact>,
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub external_edges: Vec<usize>,
    pub overall_factor: f64,
    pub external_connections: Vec<(Option<usize>, Option<usize>)>,
    pub loop_momentum_basis: LoopMomentumBasis,
    pub vertex_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub edge_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub derived_data: DerivedGraphData,
    pub shifts: (usize, usize, usize),
}

impl Graph {
    fn generate_vertex_slots(&mut self) {
        let (v, s) = self
            .vertices
            .iter()
            .fold((vec![], self.shifts), |(mut acc, shifts), v| {
                let (e, new_shifts) = v.generate_vertex_slots(shifts);
                acc.push(e);
                (acc, new_shifts)
            });
        self.shifts = s;
        self.derived_data.vertex_slots = Some(v);
    }

    pub fn dot(&self) -> String {
        let mut dot = String::new();
        dot.push_str("digraph G {\n");
        for edge in &self.edges {
            let from = self.vertices[edge.vertices[0]].name.clone();
            let to = self.vertices[edge.vertices[1]].name.clone();
            dot.push_str(&format!(
                "\"{}\" -> \"{}\" [label=\"{}\"];\n",
                from, to, edge.name
            ));
        }
        dot.push_str("}\n");
        dot
    }

    pub fn dot_internal(&self) -> String {
        let mut dot = String::new();
        dot.push_str("digraph G {\n");
        for edge in &self.edges {
            let from = self.vertices[edge.vertices[0]].name.clone();
            let to = self.vertices[edge.vertices[1]].name.clone();
            dot.push_str(&format!(
                "\"{}\" -> \"{}\" [label=\"{}\"];\n",
                self.vertex_name_to_position[&from],
                self.vertex_name_to_position[&to],
                self.edge_name_to_position[&edge.name]
            ));
        }
        dot.push_str("}\n");
        dot
    }

    pub fn dot_internal_vertices(&self) -> String {
        let mut dot = String::new();
        // let mut seen_edges = HashSet::new();
        dot.push_str("digraph G {\n");
        for (vi, v) in self.vertices.iter().enumerate() {
            dot.push_str(&format!("\"{}\" [label=\"{}: {:?}\"] \n", vi, vi, v.edges));
            for e in &v.edges {
                let [from, to] = self.edges[*e].vertices;
                // if seen_edges.insert(e) {
                dot.push_str(&format!("\"{}\" -> \"{}\" [label=\"{}\"];\n", from, to, e));
                // }
            }
        }
        dot.push_str("}\n");
        dot
    }
    pub fn from_serializable_graph(model: &model::Model, graph: &SerializableGraph) -> Graph {
        // First build vertices
        let mut vertices: Vec<Vertex> = vec![];
        let mut vertex_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState> =
            HashMap::default();
        for vertex in &graph.vertices {
            let vertex = Vertex::from_serializable_vertex(model, vertex);
            vertex_name_to_position.insert(vertex.name.clone(), vertices.len());
            vertices.push(vertex);
        }

        let mut g = Graph {
            name: graph.name.clone(),
            vertices,
            edges: vec![],
            external_edges: vec![],
            overall_factor: graph.overall_factor,
            external_connections: vec![],
            loop_momentum_basis: LoopMomentumBasis {
                basis: vec![],
                edge_signatures: vec![],
            },
            vertex_name_to_position,
            edge_name_to_position: HashMap::default(),
            derived_data: DerivedGraphData::new_empty(),
            shifts: (0, 0, MAX_COLOR_INNER_CONTRACTIONS),
        };

        // let mut edges: Vec<Edge> = vec![];
        let mut edge_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState> =
            HashMap::default();
        for edge in &graph.edges {
            let edge = Edge::from_serializable_edge(model, &g, edge);
            edge_name_to_position.insert(edge.name.clone(), g.edges.len());
            g.edges.push(edge);
        }

        debug!("Loaded {} edges", g.edges.len());
        debug!("Loaded graph: {}", g.dot());
        g.external_edges = g
            .edges
            .iter()
            .enumerate()
            .filter_map(|(i, e)| {
                if e.edge_type == EdgeType::Incoming || e.edge_type == EdgeType::Outgoing {
                    Some(i)
                } else {
                    None
                }
            })
            .collect();
        for (i_e, e) in g.edges.iter().enumerate() {
            edge_name_to_position.insert(e.name.clone(), i_e);
        }
        for (vertex, serializable_vertex) in g.vertices.iter_mut().zip(graph.vertices.iter()) {
            vertex.edges = serializable_vertex
                .edges
                .iter()
                .map(|e| *edge_name_to_position.get(e).unwrap())
                .collect();
        }
        g.edge_name_to_position = edge_name_to_position;

        let mut edge_signatures: Vec<(Vec<isize>, Vec<isize>)> =
            vec![(vec![], vec![]); graph.edges.len()];
        for (e_name, sig) in graph.edge_signatures.iter() {
            edge_signatures[g.get_edge_position(e_name).unwrap()] = sig.clone();
        }
        g.loop_momentum_basis.edge_signatures = edge_signatures;

        g.external_connections = graph
            .external_connections
            .iter()
            .map(|(v1, v2)| {
                (
                    v1.as_ref().map(|v| {
                        g.get_vertex_position(v).unwrap_or_else(|| {
                            panic!("Cannot find vertex named {} in this graph.", v)
                        })
                    }),
                    v2.as_ref().map(|v| {
                        g.get_vertex_position(v).unwrap_or_else(|| {
                            panic!("Cannot find vertex named {} in this graph.", v)
                        })
                    }),
                )
            })
            .collect();

        g.loop_momentum_basis.basis = graph
            .loop_momentum_basis
            .iter()
            .map(|e| g.get_edge_position(e).unwrap())
            .collect();

        g.generate_vertex_slots();

        // panic!("{:?}", g.edge_name_to_position);

        g
    }

    #[inline]
    pub fn get_vertex(&self, name: &SmartString<LazyCompact>) -> Option<&Vertex> {
        match self.vertex_name_to_position.get(name) {
            Some(position) => Some(&self.vertices[*position]),
            None => None,
        }
    }

    #[inline]
    pub fn get_vertex_position(&self, name: &SmartString<LazyCompact>) -> Option<usize> {
        self.vertex_name_to_position.get(name).copied()
    }

    #[inline]
    pub fn get_edge(&self, name: &SmartString<LazyCompact>) -> Option<&Edge> {
        match self.edge_name_to_position.get(name) {
            Some(position) => Some(&self.edges[*position]),
            None => None,
        }
    }

    #[inline]
    pub fn get_edge_position(&self, name: &SmartString<LazyCompact>) -> Option<usize> {
        self.edge_name_to_position.get(name).copied()
    }

    #[inline]
    pub fn compute_emr<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
    ) -> Vec<ThreeMomentum<F<T>>> {
        let lmb_specification = LoopMomentumBasisSpecification::Literal(&self.loop_momentum_basis);
        self.compute_emr_in_lmb(loop_moms, external_moms, &lmb_specification)
    }

    #[inline]
    pub fn compute_emr_in_lmb<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
        lmb_specification: &LoopMomentumBasisSpecification,
    ) -> Vec<ThreeMomentum<F<T>>> {
        //do we really want 4-momenta here? -lucien
        let lmb = lmb_specification.basis(self);

        lmb.edge_signatures
            .iter()
            .map(|sig| compute_three_momentum_from_four(sig, loop_moms, external_moms))
            .collect()
    }

    #[inline]
    pub fn compute_onshell_energies<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
    ) -> Vec<F<T>> {
        let lmb_sepcification = LoopMomentumBasisSpecification::Literal(&self.loop_momentum_basis);
        self.compute_onshell_energies_in_lmb(loop_moms, external_moms, &lmb_sepcification)
    }

    #[inline]
    pub fn get_mass_vector(&self) -> Vec<Option<Complex<F<f64>>>> {
        self.edges.iter().map(|e| e.particle.mass.value).collect()
    }

    #[inline]
    pub fn get_real_mass_vector(&self) -> Vec<F<f64>> {
        self.edges
            .iter()
            .map(|e| match e.particle.mass.value {
                Some(mass) => mass.re,
                None => F::from_f64(0.0),
            })
            .collect()
    }

    #[inline]
    pub fn compute_onshell_energies_in_lmb<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
        lmb_specification: &LoopMomentumBasisSpecification,
    ) -> Vec<F<T>> {
        let lmb = lmb_specification.basis(self);

        lmb.edge_signatures
            .iter()
            .map(|sig| compute_four_momentum_from_three(sig, loop_moms, external_moms))
            .zip(self.edges.iter())
            .map(|(emr_mom, edge)| match edge.edge_type {
                EdgeType::Virtual => {
                    emr_mom
                        .spatial
                        .on_shell_energy(edge.particle.mass.value.map(|m| {
                            if m.im.is_non_zero() {
                                panic!("Complex masses not yet supported in gammaLoop")
                            }
                            F::<T>::from_ff64(m.re)
                        }))
                        .value
                }
                _ => emr_mom.temporal.value, // a wierd way of just obtaining the energy of the external particles
            })
            .collect()
    }

    #[inline]
    pub fn compute_energy_product<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
    ) -> F<T> {
        let lmb_specification = LoopMomentumBasisSpecification::Literal(&self.loop_momentum_basis);
        self.compute_energy_product_in_lmb(loop_moms, external_moms, &lmb_specification)
    }

    #[inline]
    pub fn compute_energy_product_in_lmb<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
        lmb_specification: &LoopMomentumBasisSpecification,
    ) -> F<T> {
        let all_energies =
            self.compute_onshell_energies_in_lmb(loop_moms, external_moms, lmb_specification);

        self.edges
            .iter()
            .zip(all_energies.iter())
            .filter(|(e, _)| e.edge_type == EdgeType::Virtual)
            .map(|(_, val)| F::<T>::from_f64(2.) * val) // why times 2? -lucien
            .fold(F::<T>::from_f64(1.), |acc, x| acc * x)
    }

    #[inline]
    pub fn get_edge_type_list(&self) -> Vec<EdgeType> {
        self.edges.iter().map(|e| e.edge_type.clone()).collect()
    }

    pub fn generate_loop_momentum_bases(&mut self) {
        let loop_number = self.loop_momentum_basis.basis.len();
        let num_edges = self.edges.len();
        let external_signature_length = self.loop_momentum_basis.edge_signatures[0].1.len();

        // the full virtual signature matrix in the form of a flattened vector
        let signature_matrix_flattened = self
            .loop_momentum_basis
            .edge_signatures
            .iter()
            .flat_map(|sig| sig.0.iter().map(|s| *s as f64).collect_vec())
            .collect_vec();

        // convert to dmatrix
        let signature_matrix =
            DMatrix::from_row_slice(num_edges, loop_number, &signature_matrix_flattened);

        // the full external signature matrix in the form of a flattened vector
        let external_signature_matrix_flattened = self
            .loop_momentum_basis
            .edge_signatures
            .iter()
            .flat_map(|sig| sig.1.iter().map(|s| *s as f64).collect_vec())
            .collect_vec();

        // convert to dmatrix
        let external_signature_matrix = DMatrix::from_row_slice(
            num_edges,
            external_signature_length,
            &external_signature_matrix_flattened,
        );

        let possible_lmbs = self
            .edges
            .iter()
            .filter(|e| e.edge_type == EdgeType::Virtual)
            .map(|e| self.get_edge_position(&e.name).unwrap())
            .combinations(loop_number);

        let valid_lmbs = possible_lmbs
            .map(|basis| {
                let reduced_signature_matrix_flattened = basis
                    .iter()
                    .flat_map(|e| {
                        self.loop_momentum_basis.edge_signatures[*e]
                            .0
                            .iter()
                            .map(|s| *s as f64)
                    })
                    .collect_vec();

                (
                    basis,
                    DMatrix::from_row_slice(
                        loop_number,
                        loop_number,
                        &reduced_signature_matrix_flattened,
                    ),
                )
            })
            .filter(|(_basis, reduced_signature_matrix)| {
                reduced_signature_matrix.determinant() != 0. // nonzero determinant means valid lmb
            })
            .map(|(basis, reduced_signature_matrix)| {
                let mut sorted_basis = basis;
                sorted_basis.sort();
                (sorted_basis, reduced_signature_matrix)
            });

        let lmbs = valid_lmbs
            .map(|(basis, reduced_signature_matrix)| {
                let mut new_signatures = self.loop_momentum_basis.edge_signatures.clone();

                // construct the signatures

                // for this we need the reduced external signatures
                let reduced_external_signatures_vec = basis
                    .iter()
                    .flat_map(|&e| {
                        self.loop_momentum_basis.edge_signatures[e]
                            .1
                            .iter()
                            .map(|s| *s as f64)
                    })
                    .collect_vec();

                let reduced_external_signatures = DMatrix::from_row_slice(
                    loop_number,
                    external_signature_length,
                    &reduced_external_signatures_vec,
                );

                // also the inverse of the reduced_signature matrix
                let reduced_signature_matrix_inverse =
                    reduced_signature_matrix.clone().try_inverse().unwrap();

                let new_virtual_signatures =
                    signature_matrix.clone() * reduced_signature_matrix_inverse.clone();
                let new_external_signatures = external_signature_matrix.clone()
                    - signature_matrix.clone()
                        * reduced_signature_matrix_inverse.clone()
                        * reduced_external_signatures;

                // convert back to isize vecs
                for (edge_id, new_signature) in new_signatures.iter_mut().enumerate() {
                    let new_virtual_signature = new_virtual_signatures
                        .row(edge_id)
                        .iter()
                        .map(|s| s.round() as isize)
                        .collect_vec();
                    let new_external_signature = new_external_signatures
                        .row(edge_id)
                        .iter()
                        .map(|s| s.round() as isize)
                        .collect_vec();

                    *new_signature = (new_virtual_signature, new_external_signature);
                }
                LoopMomentumBasis {
                    basis,
                    edge_signatures: new_signatures,
                }
            })
            .collect_vec();

        self.derived_data.loop_momentum_bases = Some(lmbs);
    }

    pub fn generate_loop_momentum_bases_if_not_exists(&mut self) {
        if self.derived_data.loop_momentum_bases.is_none() {
            self.generate_loop_momentum_bases();
        }
    }

    pub fn generate_lmb_replacement_rules(&self) -> Vec<(Atom, Atom)> {
        self.loop_momentum_basis_replacement_rule(&self.loop_momentum_basis)
    }

    fn loop_momentum_basis_replacement_rule(&self, lmb: &LoopMomentumBasis) -> Vec<(Atom, Atom)> {
        let mut rule = vec![];

        for (i, signature) in lmb.edge_signatures.iter().enumerate() {
            rule.push((
                Atom::parse(&format!("Q{}(x{}__)", i, i)).unwrap(),
                self.replacement_rule_from_signature(i, signature),
            ));
        }

        rule
    }

    fn replacement_rule_from_signature(
        &self,
        index: usize,

        signature: &(Vec<isize>, Vec<isize>),
    ) -> Atom {
        let mut acc = Atom::new_num(0);
        for (i_l, sign) in signature.0.iter().enumerate() {
            match sign {
                1 => {
                    acc = &acc + &Atom::parse(&format!("K{}(x{}__)", i_l, index)).unwrap();
                }
                -1 => {
                    acc = &acc - &Atom::parse(&format!("K{}(x{}__)", i_l, index)).unwrap();
                }
                _ => {}
            }
        }

        for (i_e, sign) in signature.1.iter().enumerate() {
            match sign {
                1 => {
                    acc = &acc + &Atom::parse(&format!("P{}(x{}__)", i_e, index)).unwrap();
                }
                -1 => {
                    acc = &acc + &Atom::parse(&format!("P{}(x{}__)", i_e, index)).unwrap();
                }
                _ => {}
            }
        }
        acc
    }

    pub fn generate_ltd(&mut self) {
        self.derived_data.ltd_expression = Some(generate_ltd_expression(self));
    }

    pub fn denominator(self) -> Vec<(Atom, Atom)> {
        self.edges.iter().map(|e| e.denominator(&self)).collect()
    }

    pub fn generate_cff(&mut self) {
        self.derived_data.cff_expression = Some(generate_cff_expression(self).unwrap());
    }

    pub fn generate_tropical_subgraph_table(&mut self) {
        let num_virtual_loop_edges = self.get_loop_edges_iterator().count();

        let num_loops = self.loop_momentum_basis.basis.len();

        let default_weight = 0.5;
        let dod =
            num_virtual_loop_edges as f64 * default_weight - (tropical::D * num_loops) as f64 / 2.;
        let minimum_dod = 1.0;

        let weight = if dod < minimum_dod {
            (minimum_dod + (tropical::D * num_loops) as f64 / 2.) / num_virtual_loop_edges as f64
        } else {
            default_weight
        };

        let weight_guess = vec![weight; num_virtual_loop_edges];

        let table = TropicalSubgraphTable::generate_from_graph(self, &weight_guess);

        if let Ok(table) = table {
            self.derived_data.tropical_subgraph_table = Some(table);
        } else {
            warn!("Tropical subgraph table generation failed 🥥");
        }
    }

    pub fn generate_numerator(&mut self) {
        self.derived_data.numerator = Some(Numerator::generate(self));
    }

    pub fn smart_generate_numerator(&mut self) {
        if self.derived_data.numerator.is_none() {
            self.generate_numerator();
        }
    }

    #[inline]
    pub fn evaluate_ltd_expression<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
    ) -> Complex<F<T>> {
        let one = loop_moms[0].px.one();
        let zero = one.zero();
        let i = Complex::new(zero, one);
        let loop_number = self.loop_momentum_basis.basis.len();
        let prefactor = i.pow(loop_number as u64);

        prefactor
            * self.derived_data.ltd_expression.as_ref().unwrap().evaluate(
                loop_moms,
                external_moms,
                self,
            )
    }

    #[inline]
    pub fn evaluate_ltd_expression_in_lmb<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
        lmb_specification: &LoopMomentumBasisSpecification,
    ) -> Complex<F<T>> {
        let one = loop_moms[0].px.one();
        let zero = one.zero();
        let i = Complex::new(zero, one);
        let loop_number = self.loop_momentum_basis.basis.len();
        let prefactor = i.pow(loop_number as u64);

        prefactor
            * self
                .derived_data
                .ltd_expression
                .as_ref()
                .unwrap()
                .evaluate_in_lmb(loop_moms, external_moms, self, lmb_specification)
    }

    #[inline]
    /// evaluates the cff expression at the given loop momenta and external momenta. The loop momenta are assumed to be in the loop momentum basis specified, and have irrelevant energy components. The output is a vector of complex numbers, one for each orientation.
    pub fn evaluate_cff_orientations<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        lmb_specification: &LoopMomentumBasisSpecification,
        debug: usize,
    ) -> Vec<F<T>> {
        let energy_cache = self.compute_onshell_energies_in_lmb(
            &sample.loop_moms,
            &sample.external_moms,
            lmb_specification,
        );

        self.derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .evaluate_orientations(&energy_cache, debug)
    }

    pub fn numerator_substitute_model_params(&mut self, model: &Model) {
        if let Some(numerator) = self.derived_data.numerator.as_mut() {
            numerator.substitute_model_params(model);
        }
    }

    fn emr_from_lmb<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        lmb: &LoopMomentumBasis,
    ) -> Vec<FourMomentum<F<T>>> {
        let massless = lmb.to_massless_emr(sample);
        self.edges
            .iter()
            .zip(massless)
            .map(|(edge, emr)| match edge.edge_type {
                EdgeType::Virtual => {
                    emr.spatial
                        .into_on_shell_four_momentum(edge.particle.mass.value.map(|m| {
                            if m.im.is_non_zero() {
                                panic!("Complex masses not yet supported in gammaLoop")
                            }
                            F::<T>::from_ff64(m.re)
                        }))
                }
                _ => emr,
            })
            .collect()
    }

    pub fn evaluate_model_params(&mut self, _model: &Model) {}

    pub fn evaluate_numerator<T: FloatLike>(&self, emr: Vec<FourMomentum<F<T>>>) -> Complex<F<T>> {
        self.derived_data
            .numerator
            .as_ref()
            .unwrap()
            .evaluate(&emr, self)
    }

    #[inline]
    /// evaluates the numerator at the given loop momenta and external momenta. The loop momenta are assumed to be in the loop momentum basis specified, and have irrelevant energy components. The output is a vector of complex numbers, one for each orientation.
    pub fn evaluate_numerator_orientations<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        lmb_specification: &LoopMomentumBasisSpecification,
    ) -> Vec<Complex<F<T>>> {
        let mut out = vec![];
        let lmb = lmb_specification.basis(self);

        let emr = self.emr_from_lmb(sample, lmb);

        // debug!("Numerator: {}", numerator);

        for orient in self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .map(|e| e.orientation.clone())
        {
            let mut emr = emr.clone();
            for ((i, _), sign) in self.get_virtual_edges_iterator().zip(orient.into_iter()) {
                if !sign {
                    emr[i].temporal.value.negate()
                }
            }

            out.push(self.evaluate_numerator(emr));
        }

        out
    }

    #[inline]
    /// evaluates the cff expression at the given loop momenta and external momenta. The loop momenta are assumed to be in the loop momentum basis specified, and have irrelevant energy components. The output is a vector of complex numbers, one for each orientation.
    pub fn evaluate_cff_expression_in_lmb<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        lmb_specification: &LoopMomentumBasisSpecification,
        debug: usize,
    ) -> Complex<F<T>> {
        let one = sample.one();
        let zero = one.zero();
        let i = Complex::new(zero.clone(), one.clone());

        let loop_number = self.loop_momentum_basis.basis.len();
        let internal_vertex_number = self.vertices.len() - self.external_connections.len();

        let prefactor =
            i.pow(loop_number as u64) * (-i.ref_one()).pow(internal_vertex_number as u64 - 1);

        // here numerator evaluation can be weaved into the summation
        let res = prefactor
            * self
                .evaluate_cff_orientations(sample, lmb_specification, debug)
                .into_iter()
                .zip(
                    // self.evaluate_numerator_orientations(sample, lmb_specification),
                    0.., // dummy values for performance test
                )
                .map(|(cff, _num)| {
                    let zero = cff.zero();
                    Complex::new(cff, zero)
                })
                .reduce(|acc, e| acc + &e)
                .unwrap_or_else(|| panic!("no orientations to evaluate"));

        if debug > 1 {
            println!("sum over all orientations including numerator: {}", &res)
        }

        res
    }

    #[inline]
    pub fn evaluate_cff_expression<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        debug: usize,
    ) -> Complex<F<T>> {
        let lmb_specification = LoopMomentumBasisSpecification::Literal(&self.loop_momentum_basis);
        self.evaluate_cff_expression_in_lmb(sample, &lmb_specification, debug)
    }

    pub fn process_numerator(&mut self, model: &Model) {
        self.smart_generate_numerator();
        let mut numerator = self.derived_data.numerator.clone().unwrap();

        numerator.process(self, model);
        self.derived_data.numerator = Some(numerator);
    }

    pub fn load_derived_data(&mut self, path: &Path) -> Result<(), Report> {
        let derived_data_path = path.join(format!("derived_data_{}.bin", self.name.as_str()));
        let derived_data = DerivedGraphData::load_from_path(&derived_data_path)?;
        self.derived_data = derived_data;

        let loaded_compiled = self
            .derived_data
            .cff_expression
            .as_mut()
            .unwrap()
            .load_compiled(path.into());

        if let Err(e) = loaded_compiled {
            warn!("could not load compiled cff: {}", e)
        }

        // if the user has edited the lmb in amplitude.yaml, this will set the right signature.
        let lmb_indices = self.loop_momentum_basis.basis.clone();
        self.set_lmb(&lmb_indices)?;
        Ok(())
    }

    // attempt to set a new loop momentum basis
    pub fn set_lmb(&mut self, lmb: &[usize]) -> Result<(), Report> {
        let position = self.derived_data.search_lmb_position(lmb)?;
        self.loop_momentum_basis =
            self.derived_data.loop_momentum_bases.as_ref().unwrap()[position].clone();
        Ok(())
    }

    #[inline]
    pub fn get_virtual_edges_iterator(&self) -> impl Iterator<Item = (usize, &Edge)> {
        self.edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.edge_type == EdgeType::Virtual)
    }

    /// iterate over all edges which are part of a loop, (tree-level attachments removed)
    #[inline]
    pub fn get_loop_edges_iterator(&self) -> impl Iterator<Item = (usize, &Edge)> {
        self.edges.iter().enumerate().filter(|(index, e)| {
            e.edge_type == EdgeType::Virtual && {
                self.loop_momentum_basis.edge_signatures[*index]
                    .0
                    .iter()
                    .any(|x| *x != 0)
            }
        })
    }

    /// iterate over all edges which are virtual but do not carry a loop momentum.
    #[inline]
    pub fn get_tree_level_edges_iterator(&self) -> impl Iterator<Item = (usize, &Edge)> {
        self.edges.iter().enumerate().filter(|(index, e)| {
            e.edge_type == EdgeType::Virtual && {
                self.loop_momentum_basis.edge_signatures[*index]
                    .0
                    .iter()
                    .all(|x| *x == 0)
            }
        })
    }

    /// For a given edge, return the indices of the edges which have the same signature
    #[inline]
    pub fn is_edge_raised(&self, edge_index: usize) -> SmallVec<[usize; 2]> {
        let edge_signature = &self.loop_momentum_basis.edge_signatures[edge_index];
        let virtual_edges = self.get_virtual_edges_iterator();

        virtual_edges
            .filter(|(index, _)| {
                self.loop_momentum_basis.edge_signatures[*index] == *edge_signature
                    && *index != edge_index
            })
            .map(|(index, _)| index)
            .collect()
    }

    /// Returns groups of edges which all have the same signature
    pub fn group_edges_by_signature(&self) -> Vec<SmallVec<[usize; 3]>> {
        let mut edges: Vec<usize> = self
            .get_virtual_edges_iterator()
            .map(|(index, _)| index)
            .collect();

        let mut grouped_edges = Vec::with_capacity(edges.len());

        while !edges.is_empty() {
            let current_edge = edges.remove(0);

            let mut group = smallvec![current_edge];
            let mut index = 0;

            while index < edges.len() {
                if self.loop_momentum_basis.edge_signatures[current_edge]
                    == self.loop_momentum_basis.edge_signatures[edges[index]]
                {
                    group.push(edges.remove(index));
                } else {
                    index += 1;
                }
            }

            grouped_edges.push(group);
        }

        grouped_edges
    }

    pub fn generate_edge_groups(&mut self) {
        self.derived_data.edge_groups = Some(self.group_edges_by_signature());
    }

    pub fn generate_esurface_data(&mut self) -> Result<(), Report> {
        let data = generate_esurface_data(self, &self.get_cff().esurfaces)?;
        self.derived_data.esurface_derived_data = Some(data);

        Ok(())
    }

    pub fn is_edge_incoming(&self, edge: usize, vertex: usize) -> bool {
        self.edges[edge].is_incoming_to(vertex)
    }

    // helper function
    #[inline]
    pub fn get_cff(&self) -> &CFFExpression {
        self.derived_data.cff_expression.as_ref().unwrap()
    }

    #[inline]
    pub fn get_tropical_subgraph_table(&self) -> &TropicalSubgraphTable {
        self.derived_data.tropical_subgraph_table.as_ref().unwrap()
    }

    #[inline]
    pub fn get_esurface_derived_data(&self) -> &EsurfaceDerivedData {
        self.derived_data.esurface_derived_data.as_ref().unwrap()
    }

    #[inline]
    pub fn get_existing_esurfaces<T: FloatLike>(
        &self,
        externals: &[FourMomentum<F<T>>],
        e_cm: F<f64>,
        debug: usize,
    ) -> ExistingEsurfaces {
        get_existing_esurfaces(
            &self.get_cff().esurfaces,
            self.get_esurface_derived_data(),
            externals,
            &self.loop_momentum_basis,
            debug,
            e_cm,
        )
    }

    #[inline]
    pub fn get_maximal_overlap(
        &self,
        externals: &[FourMomentum<F<f64>>],
        e_cm: F<f64>,
        debug: usize,
    ) -> OverlapStructure {
        let existing_esurfaces = self.get_existing_esurfaces(externals, e_cm, debug);
        find_maximal_overlap(
            &self.loop_momentum_basis,
            &existing_esurfaces,
            &self.get_cff().esurfaces,
            &self.get_mass_vector(),
            externals,
            debug,
        )
    }

    pub fn get_dep_mom_expr(&self) -> (usize, ExternalShift) {
        let external_edges = self
            .edges
            .iter()
            .enumerate()
            .filter(|(_index, edge)| edge.edge_type != EdgeType::Virtual)
            .collect_vec();

        // find the external leg which does not appear in it's own signature
        let (_, (dep_mom, _)) = external_edges
            .iter()
            .enumerate()
            .find(|(external_index, (index, _edge))| {
                self.loop_momentum_basis.edge_signatures[*index].1[*external_index] != 1
            })
            .unwrap_or_else(|| panic!("could not determine dependent momenta"));

        let dep_mom_signature = &self.loop_momentum_basis.edge_signatures[*dep_mom].1;

        let external_shift = external_edges
            .iter()
            .zip(dep_mom_signature.iter())
            .filter(|(_external_edge, dep_mom_sign)| **dep_mom_sign != 0)
            .map(|((external_edge, _), dep_mom_sign)| (*external_edge, *dep_mom_sign as i64))
            .collect();

        (*dep_mom, external_shift)
    }

    pub fn build_compiled_expression(
        &mut self,
        export_path: PathBuf,
        compile_seperate_orientations: bool,
    ) -> Result<(), Report> {
        let params = self.build_params_for_cff();
        match self.derived_data.cff_expression.as_mut() {
            Some(cff) => cff.build_compiled_experssion::<f64>(
                &params,
                export_path,
                compile_seperate_orientations,
            ),
            None => {
                self.generate_cff();
                self.build_compiled_expression(export_path, compile_seperate_orientations)
            }
        }
    }

    pub fn build_params_for_cff(&self) -> Vec<Atom> {
        self.edges
            .iter()
            .enumerate()
            .map(|(id, edge)| match edge.edge_type {
                EdgeType::Virtual => Atom::parse(&format!("E{}", id)).unwrap(),
                _ => Atom::parse(&format!("p{}", id)).unwrap(),
            })
            .collect()
    }
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DerivedGraphData {
    pub loop_momentum_bases: Option<Vec<LoopMomentumBasis>>,
    pub cff_expression: Option<CFFExpression>,
    pub ltd_expression: Option<LTDExpression>,
    pub tropical_subgraph_table: Option<TropicalSubgraphTable>,
    pub edge_groups: Option<Vec<SmallVec<[usize; 3]>>>,
    pub esurface_derived_data: Option<EsurfaceDerivedData>,
    pub static_counterterm: Option<static_counterterm::CounterTerm>,
    pub vertex_slots: Option<Vec<VertexSlots>>,
    pub numerator: Option<Numerator>,
}

impl DerivedGraphData {
    fn new_empty() -> Self {
        DerivedGraphData {
            loop_momentum_bases: None,
            cff_expression: None,
            ltd_expression: None,
            tropical_subgraph_table: None,
            edge_groups: None,
            esurface_derived_data: None,
            numerator: None,
            static_counterterm: None,
            vertex_slots: None,
        }
    }

    pub fn to_serializable(&self) -> SerializableDerivedGraphData {
        SerializableDerivedGraphData {
            loop_momentum_bases: self
                .loop_momentum_bases
                .clone()
                .map(|lmbs| lmbs.iter().map(|lmb| lmb.to_serializable()).collect_vec()),
            cff_expression: self.cff_expression.clone(),
            ltd_expression: self.ltd_expression.clone().map(|ltd| ltd.to_serializable()),
            tropical_subgraph_table: self.tropical_subgraph_table.clone(),
            edge_groups: self.edge_groups.clone().map(|groups| {
                groups
                    .iter()
                    .map(|group| group.clone().into_iter().collect())
                    .collect()
            }),
            esurface_derived_data: self.esurface_derived_data.clone(),
            numerator: self.numerator.clone(),
            vertex_slots: self.vertex_slots.clone(),
        }
    }

    pub fn from_serializable(serializable: SerializableDerivedGraphData) -> Self {
        DerivedGraphData {
            loop_momentum_bases: serializable.loop_momentum_bases.map(|lmbs| {
                lmbs.iter()
                    .map(LoopMomentumBasis::from_serializable)
                    .collect_vec()
            }),
            cff_expression: serializable.cff_expression,
            ltd_expression: serializable
                .ltd_expression
                .map(LTDExpression::from_serializable),
            tropical_subgraph_table: serializable.tropical_subgraph_table,
            edge_groups: serializable
                .edge_groups
                .map(|groups| groups.into_iter().map(|group| group.into()).collect()),
            esurface_derived_data: serializable.esurface_derived_data,
            numerator: serializable.numerator,
            static_counterterm: None,
            vertex_slots: serializable.vertex_slots,
        }
    }

    pub fn load_from_path(path: &Path) -> Result<Self, Report> {
        match std::fs::read(path) {
            Ok(derived_data_bytes) => {
                let derived_data: SerializableDerivedGraphData =
                    bincode::deserialize(&derived_data_bytes)?;
                Ok(Self::from_serializable(derived_data))
            }
            Err(_) => {
                warn!("no derived data found");
                Ok(Self::new_empty())
            }
        }
    }

    // search the lmb position in the list of lmbs
    fn search_lmb_position(&self, potential_lmb: &[usize]) -> Result<usize, Report> {
        match &self.loop_momentum_bases {
            None => Err(eyre!("loop momentum bases not yet generated")),
            Some(lmbs) => {
                let sorted_potential_lmb = potential_lmb.iter().sorted().collect_vec();

                for (position, lmb) in lmbs.iter().enumerate() {
                    let sorted_lmb = lmb.basis.iter().sorted().collect_vec();

                    if sorted_lmb == sorted_potential_lmb {
                        return Ok(position);
                    }
                }
                Err(eyre!("lmb not found"))
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableDerivedGraphData {
    pub loop_momentum_bases: Option<Vec<SerializableLoopMomentumBasis>>,
    pub cff_expression: Option<CFFExpression>,
    pub ltd_expression: Option<SerializableLTDExpression>,
    pub tropical_subgraph_table: Option<TropicalSubgraphTable>,
    pub edge_groups: Option<Vec<Vec<usize>>>,
    pub esurface_derived_data: Option<EsurfaceDerivedData>,
    pub numerator: Option<Numerator>,
    pub vertex_slots: Option<Vec<VertexSlots>>,
}

#[derive(Debug, Clone)]
pub struct LoopMomentumBasis {
    pub basis: Vec<usize>,
    pub edge_signatures: Vec<(Vec<isize>, Vec<isize>)>,
}

impl LoopMomentumBasis {
    pub fn to_serializable(&self) -> SerializableLoopMomentumBasis {
        SerializableLoopMomentumBasis {
            basis: self.basis.clone(),
            edge_signatures: self.edge_signatures.clone(),
        }
    }

    pub fn from_serializable(serializable: &SerializableLoopMomentumBasis) -> LoopMomentumBasis {
        LoopMomentumBasis {
            basis: serializable.basis.clone(),
            edge_signatures: serializable.edge_signatures.clone(),
        }
    }

    pub fn to_massless_emr<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
    ) -> Vec<FourMomentum<F<T>>> {
        self.edge_signatures
            .iter()
            .map(|sig| {
                compute_four_momentum_from_three(sig, &sample.loop_moms, &sample.external_moms)
            })
            .collect()
    }

    pub fn pattern(&self, edge_id: usize) -> Pattern {
        let (loop_signature, external_signature) = self.edge_signatures[edge_id].clone();

        let mut atom = Atom::new_num(0);

        for (i, sign) in loop_signature.iter().enumerate() {
            match sign {
                1 => {
                    atom = &atom + &Atom::parse(&format!("K({},x_))", i)).unwrap();
                }
                -1 => {
                    atom = &atom - &Atom::parse(&format!("K({},x_)", i)).unwrap();
                }
                _ => {}
            }
        }

        for (i, sign) in external_signature.iter().enumerate() {
            match sign {
                1 => {
                    atom = &atom + &Atom::parse(&format!("P({},x_)", i)).unwrap();
                }
                -1 => {
                    atom = &atom - &Atom::parse(&format!("P({},x_)", i)).unwrap();
                }
                _ => {}
            }
        }

        atom.into_pattern()
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SerializableLoopMomentumBasis {
    pub basis: Vec<usize>,
    pub edge_signatures: Vec<(Vec<isize>, Vec<isize>)>,
}

// helper enum for different ways of specifying an lmb to compute stuff in.
pub enum LoopMomentumBasisSpecification<'a> {
    Literal(&'a LoopMomentumBasis),
    FromList(usize),
}

impl<'a> LoopMomentumBasisSpecification<'a> {
    pub fn basis(&self, graph: &'a Graph) -> &'a LoopMomentumBasis {
        match self {
            LoopMomentumBasisSpecification::Literal(basis) => basis,
            LoopMomentumBasisSpecification::FromList(idx) => &graph
                .derived_data
                .loop_momentum_bases
                .as_ref()
                .unwrap_or_else(|| panic!("Loop momentum bases not yet generated"))[*idx],
        }
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;

    use crate::cross_section::OutputMetaData;

    use super::*;

    fn model_sm() -> Model {
        let path = Path::new("./src/test_resources/lbl/");
        let output_meta_data: OutputMetaData =
            serde_yaml::from_reader(File::open(path.join("output_metadata.yaml")).unwrap())
                .unwrap();
        Model::from_file(String::from(
            path.join(format!(
                "sources/model/{}.yaml",
                output_meta_data.model_name
            ))
            .to_str()
            .unwrap(),
        ))
        .unwrap()
    }

    #[test]
    fn vertex_rule() {
        let model = model_sm();

        let v = &model.vertex_rules[model.vertex_rule_name_to_position["V_44"]];

        let _vertex = Vertex {
            name: "v".into(),
            vertex_info: VertexInfo::InteractonVertexInfo(InteractionVertexInfo {
                vertex_rule: v.clone(),
            }),
            edges: vec![2, 3],
        };

        let _lorentz = v.lorentz_structures[0].structure.clone();
        // println!("{}");
    }
}
