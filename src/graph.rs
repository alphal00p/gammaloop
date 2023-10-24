use crate::{
    cff::{generate_cff_expression, CFFExpression},
    ltd::{generate_ltd_expression, LTDExpression},
    model,
    utils::{compute_momentum, FloatLike},
};
use ahash::RandomState;
use color_eyre::{Help, Report};
use enum_dispatch::enum_dispatch;
use eyre::eyre;
use itertools::Itertools;
use lorentz_vector::LorentzVector;
use nalgebra::DMatrix;
#[allow(unused_imports)]
use num_traits::Float;
use serde::{Deserialize, Serialize};
use smartstring::{LazyCompact, SmartString};
use std::{collections::HashMap, sync::Arc};

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
                        let energy_squared =
                            emr_mom.spatial_squared() + Into::<T>::into(mass_value * mass_value);
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
    pub fn get_edge_type_list(&self) -> Vec<EdgeType> {
        self.edges.iter().map(|e| e.edge_type.clone()).collect()
    }

    pub fn generate_loop_momentum_bases(&mut self) {
        let loop_number = self.loop_momentum_basis.basis.len();
        let num_edges = self.edges.len();

        let lmbs = self
            .edges
            .iter()
            .filter(|e| e.edge_type == EdgeType::Virtual)
            .map(|e| self.get_edge_position(&e.name).unwrap())
            .combinations(loop_number)
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
                let inverted_reduced_signature_matrix =
                    reduced_signature_matrix.try_inverse().unwrap();

                let virtual_signatures_slice = self
                    .loop_momentum_basis
                    .edge_signatures
                    .iter()
                    .flat_map(|(virt, _real)| virt.iter().map(|s| *s as f64).collect_vec())
                    .collect_vec();

                let virtual_signatures_matrix =
                    DMatrix::from_row_slice(num_edges, loop_number, &virtual_signatures_slice);

                let new_signature_matrix =
                    virtual_signatures_matrix * inverted_reduced_signature_matrix.clone();

                let mut new_edge_signatures = self.loop_momentum_basis.edge_signatures.clone();
                for (edge_index, signature) in new_edge_signatures.iter_mut().enumerate() {
                    if basis.contains(&edge_index) {
                        signature.0 = inverted_reduced_signature_matrix
                            .row(basis.iter().position(|&e| e == edge_index).unwrap())
                            .iter()
                            .map(|s| *s as isize)
                            .collect_vec();
                    } else {
                        signature.0 = new_signature_matrix
                            .row(edge_index)
                            .iter()
                            .map(|s| *s as isize)
                            .collect_vec();
                    }
                }

                let sorted_basis = basis.into_iter().sorted().collect_vec();
                LoopMomentumBasis {
                    basis: sorted_basis,
                    edge_signatures: new_edge_signatures,
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
        self.derived_data.cff_expression = Some(generate_cff_expression(&self).unwrap());
    }

    #[inline]
    pub fn evaluate_ltd_expression<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
    ) -> T {
        self.derived_data
            .ltd_expression
            .as_ref()
            .unwrap()
            .evaluate(loop_moms, external_moms, self)
    }

    #[inline]
    pub fn evaluate_cff_expression<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        independent_external_momenta: &[LorentzVector<T>],
    ) -> T {
        let mut energy_cache = vec![T::zero(); self.edges.len()];

        // this can be simplified, since the flip is only needed in the external part of the enerrgy cache

        let mut flipped_externals = Vec::with_capacity(independent_external_momenta.len() + 1);
        let mut non_flipped_externals = independent_external_momenta.clone().to_vec();
        non_flipped_externals.push(LorentzVector::<T>::new());

        for (index, edge) in self
            .edges
            .iter()
            .enumerate()
            .filter(|(_index, edge)| edge.edge_type != EdgeType::Virtual)
        {
            if index < independent_external_momenta.len() - 1 {
                match edge.edge_type {
                    EdgeType::Incoming => {
                        flipped_externals.push(independent_external_momenta[index]);
                    }
                    EdgeType::Outgoing => {
                        flipped_externals.push(-independent_external_momenta[index]);
                    }
                    _ => unreachable!(),
                }
            } else {
                let dependend_momenta = -flipped_externals
                    .iter()
                    .fold(LorentzVector::<T>::new(), |acc, p| acc + p); // maybe implement std::iter::sum for LorentzVector
                flipped_externals.push(dependend_momenta);
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
                            &non_flipped_externals,
                        )
                        .spatial_squared()
                            + Into::<T>::into(mass_value * mass_value);
                        energy_squared.sqrt()
                    } else {
                        compute_momentum(
                            &self.loop_momentum_basis.edge_signatures[index],
                            loop_moms,
                            &non_flipped_externals,
                        )
                        .spatial_distance()
                    }
                }
                _ => flipped_externals[index].t,
            };
        }

        self.derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .evaluate(&energy_cache)
    }
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DerivedGraphData {
    pub loop_momentum_bases: Option<Vec<LoopMomentumBasis>>,
    pub cff_expression: Option<CFFExpression>,
    pub ltd_expression: Option<LTDExpression>,
}

impl DerivedGraphData {
    fn new_empty() -> Self {
        DerivedGraphData {
            loop_momentum_bases: None,
            cff_expression: None,
            ltd_expression: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct LoopMomentumBasis {
    pub basis: Vec<usize>,
    pub edge_signatures: Vec<(Vec<isize>, Vec<isize>)>,
}
