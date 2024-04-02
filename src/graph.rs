use crate::{
    cff::{generate_cff_expression, CFFExpression, SerializableCFFExpression},
    ltd::{generate_ltd_expression, LTDExpression, SerializableLTDExpression},
    model::{self, Model},
    numerator::generate_numerator,
    tensor::{
        AbstractIndex, Contract, DataTensor, DenseTensor, GetTensorData, Representation,
        SetTensorData, Slot, SparseTensor, VecStructure, COLORADJ, COLORANTIFUND, COLORANTISEXT,
        COLORFUND, COLORSEXT, EUCLIDEAN,
    },
    tropical::{self, TropicalSubgraphTable},
    utils::{compute_momentum, FloatLike},
};
use ahash::RandomState;
use color_eyre::{Help, Report};
use enum_dispatch::enum_dispatch;
use eyre::eyre;
use itertools::Itertools;
use log::warn;
use lorentz_vector::LorentzVector;
use nalgebra::DMatrix;
#[allow(unused_imports)]
use num::traits::Float;
use num::Complex;

use serde::{Deserialize, Serialize};
use smartstring::{LazyCompact, SmartString};
use std::{collections::HashMap, path::Path, sync::Arc};
use symbolica::{id::Pattern, representations::Atom};

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
        let spin_structure = self
            .vertex_rule
            .lorentz_structures
            .iter()
            .map(|ls| {
                let mut atom = ls.structure.clone();

                for (i, e) in edges.iter().enumerate() {
                    let momentum_in_pattern = Pattern::parse(&format!("P(x_,{})", i + 1)).unwrap();

                    let momentum_out_pattern =
                        Pattern::parse(&format!("Q{}(lor(4,x_))", e)).unwrap();

                    atom = momentum_in_pattern.replace_all(
                        atom.as_view(),
                        &momentum_out_pattern,
                        None,
                        None,
                    );
                }

                for (i, e) in edges.iter().enumerate() {
                    let pat: Pattern = Atom::new_num((i + 1) as i64).into_pattern();

                    let dir = if graph.is_edge_incoming(*e, vertex_pos) {
                        "in"
                    } else {
                        "out"
                    };

                    let index_num = AbstractIndex::try_from(format!("{}{}", dir, e)).unwrap().0;

                    atom = pat.replace_all(
                        atom.as_view(),
                        &Atom::new_num(index_num as i64).into_pattern(),
                        None,
                        None,
                    );
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

                let t_pat = Pattern::parse("T(a_,i_,ia_)").unwrap();

                let t_rep = Pattern::parse("T(coad(a_),cof(i_),coaf(ia_))").unwrap();

                atom = t_pat.replace_all(atom.as_view(), &t_rep, None, None);

                let f_pat = Pattern::parse("f(a1_,a2_,a3_)").unwrap();

                let f_rep = Pattern::parse("f(coad(a1_),coad(a2_),coad(a3_))").unwrap();

                atom = f_pat.replace_all(atom.as_view(), &f_rep, None, None);

                let d_pat = Pattern::parse("d(a1_,a2_,a3_)").unwrap();

                let d_rep = Pattern::parse("d(coad(a1_),coad(a2_),coad(a3_))").unwrap();

                atom = d_pat.replace_all(atom.as_view(), &d_rep, None, None);

                let eps_pat = Pattern::parse("Epsilon(i1_,i2_,i3_)").unwrap();

                let eps_rep = Pattern::parse("EpsilonBar(cof(i1_),cof(i2_),cof(i3_))").unwrap();

                atom = eps_pat.replace_all(atom.as_view(), &eps_rep, None, None);

                let eps_bar_pat = Pattern::parse("EpsilonBar(ia_1,ia_2,ia_3)").unwrap();

                let eps_bar_rep =
                    Pattern::parse("Epsilon(coaf(ia_1),coaf(ia_2),coaf(ia_3))").unwrap();

                atom = eps_bar_pat.replace_all(atom.as_view(), &eps_bar_rep, None, None);

                let t6_pat = Pattern::parse("T6(a_,s_,as_)").unwrap();

                let t6_rep = Pattern::parse("T6(coad(a_),cos(s_),coas(as_))").unwrap();

                atom = t6_pat.replace_all(atom.as_view(), &t6_rep, None, None);

                let k6_pat = Pattern::parse("K6(ia1_,ia2_,s_)").unwrap();

                let k6_rep = Pattern::parse("K6(coaf(ia1_),coaf(ia2_),cos(s_))").unwrap();

                atom = k6_pat.replace_all(atom.as_view(), &k6_rep, None, None);

                let k6_bar_pat = Pattern::parse("K6Bar(as_,i1_,i2_)").unwrap();

                let k6_bar_rep = Pattern::parse("K6Bar(coas(as_),cof(i1_),cof(i2_))").unwrap();

                atom = k6_bar_pat.replace_all(atom.as_view(), &k6_bar_rep, None, None);

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

                    let id2 = Pattern::parse(&format!("Identity(x_,{})", i + 1)).unwrap();

                    let ind = match s {
                        1 => concat!(EUCLIDEAN, "("),
                        3 => concat!(COLORFUND, "("),
                        -3 => concat!(COLORANTIFUND, "("),
                        6 => concat!(COLORSEXT, "("),
                        -6 => concat!(COLORANTISEXT, "("),
                        8 => concat!(COLORADJ, "("),
                        i => panic!("Color {i}not supported "),
                    };

                    atom = id1.replace_all(
                        atom.as_view(),
                        &Pattern::parse(&format!("Identity({}{}),x_)", ind, i + 1)).unwrap(),
                        None,
                        None,
                    );

                    atom = id2.replace_all(
                        atom.as_view(),
                        &Pattern::parse(&format!("Identity(x_,{}{}))", ind, i + 1)).unwrap(),
                        None,
                        None,
                    );
                }

                let id3 = Pattern::parse("Identity(x__)").unwrap();

                atom = id3.replace_all(
                    atom.as_view(),
                    &Pattern::parse("id(x__)").unwrap(),
                    None,
                    None,
                );

                for (i, e) in edges.iter().enumerate() {
                    let pat: Pattern = Atom::new_num((i + 1) as i64).into_pattern();

                    let dir = if graph.is_edge_incoming(*e, vertex_pos) {
                        "in"
                    } else {
                        "out"
                    };

                    let index_num = AbstractIndex::try_from(format!("{}{}", dir, e)).unwrap().0;

                    atom = pat.replace_all(
                        atom.as_view(),
                        &Atom::new_num(index_num as i64).into_pattern(),
                        None,
                        None,
                    );
                }

                for i in 0..MAX_COLOR_INNER_CONTRACTIONS {
                    let pat: Pattern = Atom::new_num(-1 - i as i64).into_pattern();

                    let index_num = AbstractIndex::try_from(format!("inner{}", i)).unwrap().0;
                    atom = pat.replace_all(
                        atom.as_view(),
                        &Atom::new_num(index_num as i64).into_pattern(),
                        None,
                        None,
                    );
                }

                let euc = Pattern::parse(concat!(EUCLIDEAN, "(x_)")).unwrap();
                let euc_rep = Pattern::parse(concat!(EUCLIDEAN, "(1,x_)")).unwrap();

                atom = euc.replace_all(atom.as_view(), &euc_rep, None, None);

                let cof = Pattern::parse(concat!(COLORFUND, "(x_)")).unwrap();
                let cof_rep = Pattern::parse(concat!(COLORFUND, "(3,x_)")).unwrap();

                atom = cof.replace_all(atom.as_view(), &cof_rep, None, None);

                let coaf = Pattern::parse(concat!(COLORANTIFUND, "(x_)")).unwrap();
                let coaf_rep = Pattern::parse(concat!(COLORANTIFUND, "(3,x_)")).unwrap();

                atom = coaf.replace_all(atom.as_view(), &coaf_rep, None, None);

                let cos = Pattern::parse(concat!(COLORSEXT, "(x_)")).unwrap();
                let cos_rep = Pattern::parse(concat!(COLORSEXT, "(6,x_)")).unwrap();

                atom = cos.replace_all(atom.as_view(), &cos_rep, None, None);

                let coas = Pattern::parse(concat!(COLORANTISEXT, "(x_)")).unwrap();
                let coas_rep = Pattern::parse(concat!(COLORANTISEXT, "(6,x_)")).unwrap();

                atom = coas.replace_all(atom.as_view(), &coas_rep, None, None);

                let coa = Pattern::parse(concat!(COLORADJ, "(x_)")).unwrap();
                let coa_rep = Pattern::parse(concat!(COLORADJ, "(8,x_)")).unwrap();

                atom = coa.replace_all(atom.as_view(), &coa_rep, None, None);

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

    pub fn numerator(&self, graph: &Graph) -> Atom {
        let num = *graph.edge_name_to_position.get(&self.name).unwrap();
        match self.edge_type {
            EdgeType::Incoming => self.particle.incoming_polarization_atom(num),
            EdgeType::Outgoing => self.particle.outgoing_polarization_atom(num),
            EdgeType::Virtual => {
                let mut atom = self.propagator.numerator.clone();

                let pindex_num = AbstractIndex::try_from(format!("p{}", num)).unwrap().0;
                let pfun = Pattern::parse("P(x_)").unwrap();
                atom = pfun.replace_all(
                    atom.as_view(),
                    &Pattern::parse(&format!("Q{}(lor(4,x_))", num)).unwrap(),
                    None,
                    None,
                );

                let pslashfun = Pattern::parse("PSlash(x__)").unwrap();
                atom = pslashfun.replace_all(
                    atom.as_view(),
                    &Pattern::parse(&format!(
                        "Q{}(lor(4,{}))Gamma({},x__)",
                        num, pindex_num, pindex_num
                    ))
                    .unwrap(),
                    None,
                    None,
                );

                let pat: Pattern = Atom::new_num(1).into_pattern();

                let in_index_num = AbstractIndex::try_from(format!("in{}", num)).unwrap().0;

                atom = pat.replace_all(
                    atom.as_view(),
                    &Pattern::parse(&format!("{in_index_num}")).unwrap(),
                    None,
                    None,
                );

                let pat: Pattern = Atom::new_num(2).into_pattern();

                let out_index_num = AbstractIndex::try_from(format!("out{}", num)).unwrap().0;
                atom = pat.replace_all(
                    atom.as_view(),
                    &Pattern::parse(&format!("{out_index_num}")).unwrap(),
                    None,
                    None,
                );

                atom
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
}

impl Graph {
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
        };

        let mut edge_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState> =
            HashMap::default();
        // Then build edges
        g.edges = graph
            .edges
            .iter()
            .map(|e| Edge::from_serializable_edge(model, &g, e))
            .collect();

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
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
    ) -> Vec<LorentzVector<T>> {
        let lmb_specification = LoopMomentumBasisSpecification::Literal(&self.loop_momentum_basis);
        self.compute_emr_in_lmb(loop_moms, external_moms, &lmb_specification)
    }

    #[inline]
    pub fn compute_emr_in_lmb<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
        lmb_specification: &LoopMomentumBasisSpecification,
    ) -> Vec<LorentzVector<T>> {
        let lmb = match lmb_specification {
            LoopMomentumBasisSpecification::FromList(lmb_idx) => &self
                .derived_data
                .loop_momentum_bases
                .as_ref()
                .unwrap_or_else(|| panic!("Loop momentum bases not yet generated"))[*lmb_idx],
            LoopMomentumBasisSpecification::Literal(basis) => basis,
        };

        lmb.edge_signatures
            .iter()
            .map(|sig| compute_momentum(sig, loop_moms, external_moms))
            .collect()
    }

    #[inline]
    pub fn compute_onshell_energies<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
    ) -> Vec<T> {
        let lmb_sepcification = LoopMomentumBasisSpecification::Literal(&self.loop_momentum_basis);
        self.compute_onshell_energies_in_lmb(loop_moms, external_moms, &lmb_sepcification)
    }

    #[inline]
    pub fn compute_onshell_energies_in_lmb<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
        lmb_specification: &LoopMomentumBasisSpecification,
    ) -> Vec<T> {
        let lmb = match lmb_specification {
            LoopMomentumBasisSpecification::FromList(lmb_idx) => &self
                .derived_data
                .loop_momentum_bases
                .as_ref()
                .unwrap_or_else(|| panic!("Loop momentum bases not yet generated"))[*lmb_idx],
            LoopMomentumBasisSpecification::Literal(basis) => basis,
        };

        lmb.edge_signatures
            .iter()
            .map(|sig| compute_momentum(sig, loop_moms, external_moms))
            .zip(self.edges.iter())
            .map(|(emr_mom, edge)| match edge.edge_type {
                EdgeType::Virtual => {
                    if let Some(mass_value) = edge.particle.mass.value {
                        if mass_value.im != 0. {
                            panic!("Complex masses not yet supported in gammaLoop")
                        }
                        let energy_squared = emr_mom.spatial_squared()
                            + Into::<T>::into(mass_value.re * mass_value.re);
                        energy_squared.sqrt()
                    } else {
                        emr_mom.spatial_distance()
                    }
                }
                _ => emr_mom.t,
            })
            .collect()
    }

    #[inline]
    pub fn compute_energy_product<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
    ) -> T {
        let lmb_specification = LoopMomentumBasisSpecification::Literal(&self.loop_momentum_basis);
        self.compute_energy_product_in_lmb(loop_moms, external_moms, &lmb_specification)
    }

    #[inline]
    pub fn compute_energy_product_in_lmb<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
        lmb_specification: &LoopMomentumBasisSpecification,
    ) -> T {
        let all_energies =
            self.compute_onshell_energies_in_lmb(loop_moms, external_moms, lmb_specification);

        self.edges
            .iter()
            .zip(all_energies.iter())
            .filter(|(e, _)| e.edge_type == EdgeType::Virtual)
            .map(|(_, val)| Into::<T>::into(2.) * val)
            .fold(Into::<T>::into(1.), |acc, x| acc * x)
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

    pub fn generate_numerator(&mut self, model: &Model) {
        self.derived_data.numerator = Some(generate_numerator(self, model));
    }

    pub fn smart_generate_numerator(&mut self, model: &Model) {
        if self.derived_data.numerator.is_none() {
            self.derived_data.numerator = Some(generate_numerator(self, model));
        }
    }

    #[inline]
    pub fn evaluate_ltd_expression<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
    ) -> Complex<T> {
        let loop_number = self.loop_momentum_basis.basis.len();
        let prefactor = Complex::new(T::zero(), T::one()).powi(loop_number as i32);

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
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
        lmb_specification: &LoopMomentumBasisSpecification,
    ) -> Complex<T> {
        let loop_number = self.loop_momentum_basis.basis.len();
        let prefactor = Complex::new(T::zero(), T::one()).powi(loop_number as i32);

        prefactor
            * self
                .derived_data
                .ltd_expression
                .as_ref()
                .unwrap()
                .evaluate_in_lmb(loop_moms, external_moms, self, lmb_specification)
    }

    #[inline]
    pub fn evaluate_cff_orientations<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        independent_external_momenta: &[LorentzVector<T>],
        lmb_specification: &LoopMomentumBasisSpecification,
    ) -> Vec<T> {
        let lmb = match lmb_specification {
            LoopMomentumBasisSpecification::FromList(lmb_idx) => &self
                .derived_data
                .loop_momentum_bases
                .as_ref()
                .unwrap_or_else(|| panic!("Loop momentum bases not yet generated"))[*lmb_idx],
            LoopMomentumBasisSpecification::Literal(basis) => basis,
        };

        let mut energy_cache = vec![T::zero(); self.edges.len()];
        // some gymnastics to account for the sign of outgoing momenta

        let mut flipped_externals = Vec::with_capacity(independent_external_momenta.len() + 1);

        for (index, edge) in self
            .edges
            .iter()
            .filter(|edge| edge.edge_type != EdgeType::Virtual)
            .enumerate()
        {
            if index < independent_external_momenta.len() {
                match edge.edge_type {
                    EdgeType::Incoming => {
                        flipped_externals.push(independent_external_momenta[index].t);
                    }
                    EdgeType::Outgoing => {
                        flipped_externals.push(-independent_external_momenta[index].t);
                    }
                    _ => unreachable!(),
                }
            } else {
                flipped_externals.push(-flipped_externals.iter().fold(T::zero(), |acc, x| acc + x));
            }
        }

        // here we still use the non_flipped_externals for the virtual edges, since otherwise we would have to change the signature matrix as well.
        for (index, edge) in self.edges.iter().enumerate() {
            energy_cache[index] = match edge.edge_type {
                EdgeType::Virtual => {
                    if let Some(mass_value) = edge.particle.mass.value {
                        let energy_squared = compute_momentum(
                            &lmb.edge_signatures[index],
                            loop_moms,
                            independent_external_momenta,
                        )
                        .spatial_squared()
                            + Into::<T>::into(mass_value.re * mass_value.re);
                        energy_squared.sqrt()
                    } else {
                        compute_momentum(
                            &lmb.edge_signatures[index],
                            loop_moms,
                            independent_external_momenta,
                        )
                        .spatial_distance()
                    }
                }
                _ => flipped_externals[index],
            };
        }

        self.derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .evaluate_orientations(&energy_cache)
    }

    #[inline]
    pub fn evaluate_cff_expression_in_lmb<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
        lmb_specification: &LoopMomentumBasisSpecification,
    ) -> Complex<T> {
        let loop_number = self.loop_momentum_basis.basis.len();
        let internal_vertex_number = self.vertices.len() - self.external_connections.len();

        let prefactor = Complex::new(T::zero(), T::one()).powi(loop_number as i32)
            * Complex::new(-T::one(), T::zero()).powi(internal_vertex_number as i32 - 1);

        // here numerator evaluation can be weaved into the summation
        prefactor
            * self
                .evaluate_cff_orientations(loop_moms, external_moms, lmb_specification)
                .into_iter()
                .sum::<T>()
    }

    #[inline]
    pub fn evaluate_cff_expression<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
    ) -> Complex<T> {
        let lmb_specification = LoopMomentumBasisSpecification::Literal(&self.loop_momentum_basis);
        self.evaluate_cff_expression_in_lmb(loop_moms, external_moms, &lmb_specification)
    }

    pub fn load_derived_data(&mut self, path: &Path) -> Result<(), Report> {
        let derived_data = DerivedGraphData::load_from_path(path)?;
        self.derived_data = derived_data;

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

    pub fn is_edge_incoming(&self, edge: usize, vertex: usize) -> bool {
        self.edges[edge].is_incoming_to(vertex)
    }
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DerivedGraphData {
    pub loop_momentum_bases: Option<Vec<LoopMomentumBasis>>,
    pub cff_expression: Option<CFFExpression>,
    pub ltd_expression: Option<LTDExpression>,
    pub tropical_subgraph_table: Option<TropicalSubgraphTable>,
    pub numerator: Option<Atom>,
}

impl DerivedGraphData {
    fn new_empty() -> Self {
        DerivedGraphData {
            loop_momentum_bases: None,
            cff_expression: None,
            ltd_expression: None,
            tropical_subgraph_table: None,
            numerator: None,
        }
    }

    pub fn to_serializable(&self) -> SerializableDerivedGraphData {
        SerializableDerivedGraphData {
            loop_momentum_bases: self
                .loop_momentum_bases
                .clone()
                .map(|lmbs| lmbs.iter().map(|lmb| lmb.to_serializable()).collect_vec()),
            cff_expression: self.cff_expression.clone().map(|cff| cff.to_serializable()),
            ltd_expression: self.ltd_expression.clone().map(|ltd| ltd.to_serializable()),
            tropical_subgraph_table: self.tropical_subgraph_table.clone(),
        }
    }

    pub fn from_serializable(serializable: SerializableDerivedGraphData) -> Self {
        DerivedGraphData {
            loop_momentum_bases: serializable.loop_momentum_bases.map(|lmbs| {
                lmbs.iter()
                    .map(LoopMomentumBasis::from_serializable)
                    .collect_vec()
            }),
            cff_expression: serializable
                .cff_expression
                .map(CFFExpression::from_serializable),
            ltd_expression: serializable
                .ltd_expression
                .map(LTDExpression::from_serializable),
            tropical_subgraph_table: serializable.tropical_subgraph_table,
            numerator: None,
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
    pub cff_expression: Option<SerializableCFFExpression>,
    pub ltd_expression: Option<SerializableLTDExpression>,
    pub tropical_subgraph_table: Option<TropicalSubgraphTable>,
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
