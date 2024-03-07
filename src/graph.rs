use crate::{
    cff::{generate_cff_expression, CFFExpression, SerializableCFFExpression},
    ltd::{generate_ltd_expression, LTDExpression, SerializableLTDExpression},
    model,
    subtraction::esurface_data::{self, EsurfaceDerivedData},
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
use num::Complex;
#[allow(unused_imports)]
use num_traits::Float;
use serde::{Deserialize, Serialize};
use smallvec::{smallvec, SmallVec};
use smartstring::{LazyCompact, SmartString};
use std::{collections::HashMap, path::Path, sync::Arc};

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
    vertices: [SmartString<LazyCompact>; 2],
}

impl SerializableEdge {
    pub fn from_edge(graph: &Graph, edge: &Edge) -> SerializableEdge {
        SerializableEdge {
            name: edge.name.clone(),
            edge_type: edge.edge_type.clone(),
            particle: edge.particle.name.clone(),
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
        self.compute_onshell_energies_in_lmb(loop_moms, external_moms, lmb_sepcification)
    }

    #[inline]
    pub fn compute_onshell_energies_in_lmb<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
        lmb_specification: LoopMomentumBasisSpecification,
    ) -> Vec<T> {
        let lmb = match lmb_specification {
            LoopMomentumBasisSpecification::FromList(lmb_idx) => &self
                .derived_data
                .loop_momentum_bases
                .as_ref()
                .unwrap_or_else(|| panic!("Loop momentum bases not yet generated"))[lmb_idx],
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
        self.compute_energy_product_in_lmb(loop_moms, external_moms, lmb_specification)
    }

    #[inline]
    pub fn compute_energy_product_in_lmb<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
        lmb_specification: LoopMomentumBasisSpecification,
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

    pub fn generate_ltd(&mut self) {
        self.derived_data.ltd_expression = Some(generate_ltd_expression(self));
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
        let data = esurface_data::generate_esurface_data(self)?;
        self.derived_data.esurface_derived_data = Some(data);

        Ok(())
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
            edge_groups: self.edge_groups.clone().map(|groups| {
                groups
                    .iter()
                    .map(|group| group.clone().into_iter().collect())
                    .collect()
            }),
            esurface_derived_data: self.esurface_derived_data.clone(),
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
            edge_groups: serializable
                .edge_groups
                .map(|groups| groups.into_iter().map(|group| group.into()).collect()),
            esurface_derived_data: serializable.esurface_derived_data,
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
    pub edge_groups: Option<Vec<Vec<usize>>>,
    pub esurface_derived_data: Option<EsurfaceDerivedData>,
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
