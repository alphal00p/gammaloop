use crate::{
    cff::{generate_cff_expression, CFFExpression, SerializableCFFExpression},
    ltd::{generate_ltd_expression, LTDExpression, SerializableLTDExpression},
    model::{self, Model},
    numerator::generate_numerator,
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
use symbolica::{
    id::Pattern,
    representations::{Atom, AtomView},
};

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

    pub fn numerator(&self, graph: &Graph) -> Atom {
        let mut atom = self.propagator.numerator.clone();

        let pindex_atom = Atom::parse(&format!(
            "p{}",
            graph.edge_name_to_position.get(&self.name).unwrap()
        ))
        .unwrap();

        let pindex_num = if let AtomView::Var(v) = pindex_atom.as_view() {
            v.get_symbol().get_id()
        } else {
            unreachable!()
        };

        let pfun = Pattern::parse("P(x_)").unwrap();
        atom = pfun.replace_all(
            atom.as_view(),
            &Pattern::parse(&format!(
                "P{}(lor(4,x_))",
                graph.edge_name_to_position.get(&self.name).unwrap()
            ))
            .unwrap(),
            None,
            None,
        );

        let pslashfun = Pattern::parse("PSlash(x__)").unwrap();
        atom = pslashfun.replace_all(
            atom.as_view(),
            &Pattern::parse(&format!(
                "P{}(lor(4,{}))Gamma({},x__)",
                graph.edge_name_to_position.get(&self.name).unwrap(),
                pindex_num,
                pindex_num
            ))
            .unwrap(),
            None,
            None,
        );

        let pat: Pattern = Atom::new_num(1).into_pattern();
        let index_atom = Atom::parse(&format!(
            "in{}",
            graph.edge_name_to_position.get(&self.name).unwrap()
        ))
        .unwrap();

        let index_num = if let AtomView::Var(v) = index_atom.as_view() {
            v.get_symbol().get_id()
        } else {
            unreachable!()
        };

        atom = pat.replace_all(
            atom.as_view(),
            &Pattern::parse(&format!("{index_num}")).unwrap(),
            None,
            None,
        );

        let pat: Pattern = Atom::new_num(2).into_pattern();
        let index_atom = Atom::parse(&format!(
            "out{}",
            graph.edge_name_to_position.get(&self.name).unwrap()
        ))
        .unwrap();

        let index_num = if let AtomView::Var(v) = index_atom.as_view() {
            v.get_symbol().get_id()
        } else {
            unreachable!()
        };

        atom = pat.replace_all(
            atom.as_view(),
            &Pattern::parse(&format!("{index_num}")).unwrap(),
            None,
            None,
        );

        atom
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

    pub fn apply_vertex_rule(&self, graph: &Graph) -> Vec<Atom> {
        match &self.vertex_info {
            VertexInfo::ExternalVertexInfo(_) => vec![],
            VertexInfo::InteractonVertexInfo(interaction_vertex_info) => {
                let info = interaction_vertex_info;
                info.vertex_rule
                    .lorentz_structures
                    .iter()
                    .map(|ls| {
                        let mut atom = ls.structure.clone();
                        for (i, e) in self.edges.iter().enumerate() {
                            let pat: Pattern = Atom::new_num((i + 1) as i64).into_pattern();
                            let dir = if self.is_edge_incoming(*e, graph) {
                                "in"
                            } else {
                                "out"
                            };
                            let index_atom = Atom::parse(&format!("{}{}", dir, e)).unwrap();
                            let index_num = if let AtomView::Var(v) = index_atom.as_view() {
                                v.get_symbol().get_id()
                            } else {
                                unreachable!()
                            };

                            atom = pat.replace_all(
                                atom.as_view(),
                                &Atom::new_num(index_num as i64).into_pattern(),
                                None,
                                None,
                            );
                        }
                        atom
                    })
                    .collect_vec()
            }
        }
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
        self.loop_momentum_basis
            .edge_signatures
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
        self.loop_momentum_basis
            .edge_signatures
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
        let all_energies = self.compute_onshell_energies(loop_moms, external_moms);

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

    pub fn generate_ltd(&mut self) {
        self.derived_data.ltd_expression = Some(generate_ltd_expression(self));
    }

    pub fn generate_cff(&mut self) {
        self.derived_data.cff_expression = Some(generate_cff_expression(self).unwrap());
    }

    pub fn generate_numerator(&mut self, model: &Model) {
        self.derived_data.numerator = Some(generate_numerator(self, model));
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
    pub fn evaluate_cff_orientations<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        independent_external_momenta: &[LorentzVector<T>],
    ) -> Vec<T> {
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
                            &self.loop_momentum_basis.edge_signatures[index],
                            loop_moms,
                            independent_external_momenta,
                        )
                        .spatial_squared()
                            + Into::<T>::into(mass_value.re * mass_value.re);
                        energy_squared.sqrt()
                    } else {
                        compute_momentum(
                            &self.loop_momentum_basis.edge_signatures[index],
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
    pub fn evaluate_cff_expression<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        independent_external_momenta: &[LorentzVector<T>],
    ) -> Complex<T> {
        let loop_number = self.loop_momentum_basis.basis.len();
        let internal_vertex_numer = self.vertices.len() - self.external_connections.len();

        let prefactor = Complex::new(T::zero(), T::one()).powi(loop_number as i32)
            * Complex::new(-T::one(), T::zero()).powi(internal_vertex_numer as i32 - 1);

        // here numerator evaluation can be weaved into the summation
        prefactor
            * self
                .evaluate_cff_orientations(loop_moms, independent_external_momenta)
                .into_iter()
                .sum::<T>()
    }

    pub fn load_derived_data(&mut self, path: &Path) -> Result<(), Report> {
        let derived_data = DerivedGraphData::load_from_path(path)?;
        self.derived_data = derived_data;
        Ok(())
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
    pub numerator: Option<Atom>,
}

impl DerivedGraphData {
    fn new_empty() -> Self {
        DerivedGraphData {
            loop_momentum_bases: None,
            cff_expression: None,
            ltd_expression: None,
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
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableDerivedGraphData {
    pub loop_momentum_bases: Option<Vec<SerializableLoopMomentumBasis>>,
    pub cff_expression: Option<SerializableCFFExpression>,
    pub ltd_expression: Option<SerializableLTDExpression>,
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
