//! Lightweight, application-independent Feynman diagrams.
//!
//! This crate deliberately stops at the structural and symbolic diagram layer.
//! Runtime parametrizations, cuts, subtraction data, and evaluator caches belong
//! to applications consuming these diagrams.

#![forbid(unsafe_code)]

use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    fmt::Write,
};

pub use feynkit_kinematics::MomentumSignature;
use feynkit_kinematics::Signature;
use feynkit_model::{Model, ModelError, Particle};
use linnet::{
    half_edge::{
        HedgeGraph, NodeIndex,
        builder::HedgeGraphBuilder,
        involution::{Flow, HedgePair, Orientation},
    },
    parser::DotGraph,
};
use serde::{Deserialize, Serialize};
use symbolica::graph::Graph as CanonicalGraph;
use thiserror::Error;

/// Stable vertex identifier in a [`FeynmanDiagram`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct VertexId(pub usize);

/// Stable edge identifier in a [`FeynmanDiagram`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct EdgeId(pub usize);

/// A particle reference retains both common public identifiers.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct ParticleReference {
    pub name: String,
    pub pdg: i64,
}

impl ParticleReference {
    pub fn new(name: impl Into<String>, pdg: i64) -> Self {
        Self {
            name: name.into(),
            pdg,
        }
    }
}

/// Whether an external vertex is in the initial or final state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ExternalState {
    Incoming,
    Outgoing,
}

impl ExternalState {
    fn as_str(self) -> &'static str {
        match self {
            Self::Incoming => "incoming",
            Self::Outgoing => "outgoing",
        }
    }

    fn parse(value: &str) -> Result<Self, DiagramError> {
        match value {
            "incoming" | "in" | "sink" => Ok(Self::Incoming),
            "outgoing" | "out" | "source" => Ok(Self::Outgoing),
            _ => Err(DiagramError::InvalidExternalState(value.to_owned())),
        }
    }
}

/// Metadata carried by an external vertex.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct ExternalLeg {
    pub index: usize,
    pub state: ExternalState,
}

/// Structural and symbolic information associated with a vertex.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct DiagramVertex {
    pub name: String,
    pub interaction: Option<String>,
    pub external: Option<ExternalLeg>,
    pub numerator: Option<String>,
}

impl DiagramVertex {
    pub fn interaction(name: impl Into<String>, rule: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            interaction: Some(rule.into()),
            external: None,
            numerator: None,
        }
    }

    pub fn external(name: impl Into<String>, index: usize, state: ExternalState) -> Self {
        Self {
            name: name.into(),
            interaction: None,
            external: Some(ExternalLeg { index, state }),
            numerator: None,
        }
    }

    pub fn is_external(&self) -> bool {
        self.external.is_some()
    }
}

/// Structural and symbolic information associated with an edge.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct DiagramEdge {
    pub particle: ParticleReference,
    pub directed: bool,
    pub numerator: Option<String>,
}

impl DiagramEdge {
    pub fn new(particle: ParticleReference, directed: bool) -> Self {
        Self {
            particle,
            directed,
            numerator: None,
        }
    }
}

/// Source and target in the canonical orientation of an edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct EdgeEndpoints {
    pub source: VertexId,
    pub target: VertexId,
}

/// Vertex color retained by [`CanonicalDiagramKey`].
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct CanonicalDiagramVertex {
    pub interaction: Option<String>,
    pub external: Option<ExternalLeg>,
}

/// Canonically labeled edge retained by [`CanonicalDiagramKey`].
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct CanonicalDiagramEdge {
    pub endpoints: EdgeEndpoints,
    pub particle: ParticleReference,
    pub directed: bool,
}

/// Deterministic colored-topology key for diagram grouping and parity checks.
///
/// Vertex names, numerators, symmetry factors, and overall factors are payload
/// rather than graph colors and are intentionally excluded.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct CanonicalDiagramKey {
    pub vertices: Vec<CanonicalDiagramVertex>,
    pub edges: Vec<CanonicalDiagramEdge>,
}

/// A spanning-forest-induced loop momentum basis.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct LoopMomentumBasis {
    pub tree_edges: Vec<EdgeId>,
    pub loop_edges: Vec<EdgeId>,
    pub external_edges: Vec<EdgeId>,
    /// One dependent external momentum is chosen for each connected component.
    pub dependent_externals: Vec<EdgeId>,
    pub edge_signatures: BTreeMap<EdgeId, MomentumSignature>,
}

/// Errors produced while constructing or transforming diagrams.
#[derive(Debug, Error)]
pub enum DiagramError {
    #[error(transparent)]
    Model(#[from] ModelError),
    #[error("vertex {vertex} is outside a diagram with {vertices} vertices")]
    UnknownVertex { vertex: usize, vertices: usize },
    #[error("edge {edge} is outside a diagram with {edges} edges")]
    UnknownEdge { edge: usize, edges: usize },
    #[error("diagram contains duplicate external index {0}")]
    DuplicateExternalIndex(usize),
    #[error("diagram contains no external leg with index {0}")]
    UnknownExternalIndex(usize),
    #[error("invalid external state '{0}'")]
    InvalidExternalState(String),
    #[error("missing required DOT attribute '{attribute}' on {target}")]
    MissingDotAttribute {
        target: String,
        attribute: &'static str,
    },
    #[error("invalid DOT attribute '{attribute}' on {target}: {value}")]
    InvalidDotAttribute {
        target: String,
        attribute: &'static str,
        value: String,
    },
    #[error("external index {0} cannot be incremented while importing DOT")]
    DotExternalIndexOverflow(usize),
    #[error("failed to parse DOT: {0}")]
    DotParse(String),
    #[error("failed to serialize or parse diagram JSON: {0}")]
    Json(#[from] serde_json::Error),
    #[error("edge {edge} references unknown particle PDG code {pdg}")]
    UnknownParticle { edge: usize, pdg: i64 },
    #[error("edge {edge} particle name '{name}' does not match PDG code {pdg} ('{expected_name}')")]
    ParticleNameMismatch {
        edge: usize,
        name: String,
        pdg: i64,
        expected_name: String,
    },
    #[error("vertex {vertex} references unknown interaction '{interaction}'")]
    UnknownInteraction { vertex: usize, interaction: String },
    #[error(
        "vertex {vertex} interaction '{interaction}' has incident signature {actual:?}, expected {expected:?}"
    )]
    InteractionSignatureMismatch {
        vertex: usize,
        interaction: String,
        actual: Vec<(i64, Option<bool>)>,
        expected: Vec<(i64, Option<bool>)>,
    },
    #[error("external vertex {vertex} has degree {degree}, expected one")]
    InvalidExternalDegree { vertex: usize, degree: usize },
    #[error("edge {edge} connects two external vertices")]
    ExternalToExternalEdge { edge: usize },
    #[error("diagram invariant failed while {operation}: {message}")]
    Invariant {
        operation: &'static str,
        message: String,
    },
    #[error("failed to format DOT")]
    Formatting(#[from] std::fmt::Error),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DiagramSerde {
    name: String,
    vertices: Vec<DiagramVertex>,
    edges: Vec<(EdgeEndpoints, DiagramEdge)>,
    symmetry_factor: u64,
    overall_factor: String,
    numerator: Option<String>,
}

/// Model-aware Feynman-diagram IR backed by Linnet's half-edge graph.
#[derive(Debug, Clone)]
pub struct FeynmanDiagram {
    name: String,
    graph: HedgeGraph<DiagramEdge, DiagramVertex>,
    symmetry_factor: u64,
    overall_factor: String,
    numerator: Option<String>,
}

impl Serialize for FeynmanDiagram {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.serde_view().serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for FeynmanDiagram {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let data = DiagramSerde::deserialize(deserializer)?;
        Self::from_serde(data).map_err(serde::de::Error::custom)
    }
}

impl FeynmanDiagram {
    pub fn builder(name: impl Into<String>) -> FeynmanDiagramBuilder {
        FeynmanDiagramBuilder::new(name)
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn symmetry_factor(&self) -> u64 {
        self.symmetry_factor
    }

    pub fn overall_factor(&self) -> &str {
        &self.overall_factor
    }

    pub fn numerator(&self) -> Option<&str> {
        self.numerator.as_deref()
    }

    pub fn underlying(&self) -> &HedgeGraph<DiagramEdge, DiagramVertex> {
        &self.graph
    }

    pub fn vertices(&self) -> impl Iterator<Item = (VertexId, &DiagramVertex)> {
        self.graph
            .iter_nodes()
            .map(|(id, _, vertex)| (VertexId(id.0), vertex))
    }

    pub fn edges(&self) -> impl Iterator<Item = (EdgeId, EdgeEndpoints, &DiagramEdge)> {
        self.graph.iter_edges().filter_map(|(pair, id, data)| {
            self.paired_endpoints(pair, data.orientation)
                .map(|endpoints| (EdgeId(id.0), endpoints, data.data))
        })
    }

    pub fn vertex(&self, id: VertexId) -> Option<&DiagramVertex> {
        (id.0 < self.graph.n_nodes()).then(|| &self.graph[NodeIndex(id.0)])
    }

    /// Transform vertex and edge metadata without changing the diagram
    /// topology or stable identifiers.
    pub fn map_data<V, E>(&self, mut map_vertex: V, mut map_edge: E) -> Result<Self, DiagramError>
    where
        V: FnMut(VertexId, &DiagramVertex) -> DiagramVertex,
        E: FnMut(EdgeId, EdgeEndpoints, &DiagramEdge) -> DiagramEdge,
    {
        let mut builder = Self::builder(&self.name)
            .symmetry_factor(self.symmetry_factor)
            .overall_factor(&self.overall_factor)
            .numerator_option(self.numerator.clone());
        for (id, vertex) in self.vertices() {
            builder.add_vertex(map_vertex(id, vertex));
        }
        for (id, endpoints, edge) in self.edges() {
            builder.add_edge(
                endpoints.source,
                endpoints.target,
                map_edge(id, endpoints, edge),
            )?;
        }
        builder.build()
    }

    /// Build a deterministic key for the model-resolved colored topology.
    ///
    /// Canonical labeling is delegated to Symbolica's graph implementation.
    /// Particle/antiparticle edge representations are normalized so reversing
    /// a directed edge with [`Self::reverse_edge`] does not change the key.
    pub fn canonical_key(&self, model: &Model) -> Result<CanonicalDiagramKey, DiagramError> {
        self.validate(model)?;

        let mut graph = CanonicalGraph::new();
        for (_, vertex) in self.vertices() {
            graph.add_node(CanonicalDiagramVertex {
                interaction: vertex.interaction.clone(),
                external: vertex.external.clone(),
            });
        }
        for (edge_id, endpoints, edge) in self.edges() {
            let (particle, base) = Self::resolve_particle(edge_id, edge, model)?;
            let (mut source, mut target) = (endpoints.source.0, endpoints.target.0);
            if edge.directed && particle.is_antiparticle() {
                std::mem::swap(&mut source, &mut target);
            }
            graph
                .add_edge(
                    source,
                    target,
                    edge.directed,
                    ParticleReference::new(&base.name, base.pdg_code),
                )
                .map_err(|error| DiagramError::Invariant {
                    operation: "canonicalizing the diagram",
                    message: error.to_string(),
                })?;
        }

        let graph = graph.canonize().graph;
        Ok(CanonicalDiagramKey {
            vertices: graph.nodes().iter().map(|node| node.data.clone()).collect(),
            edges: graph
                .edges()
                .iter()
                .map(|edge| CanonicalDiagramEdge {
                    endpoints: EdgeEndpoints {
                        source: VertexId(edge.vertices.0),
                        target: VertexId(edge.vertices.1),
                    },
                    particle: edge.data.clone(),
                    directed: edge.directed,
                })
                .collect(),
        })
    }

    /// Reverse a directed edge while preserving its physical particle flow.
    ///
    /// Reversal swaps the endpoints and replaces the referenced particle with
    /// its antiparticle. Undirected edges are unchanged. The returned diagram
    /// is fully model-validated and retains all stable identifiers.
    pub fn reverse_edge(&self, edge_id: EdgeId, model: &Model) -> Result<Self, DiagramError> {
        self.validate(model)?;
        let Some((_, _, selected)) = self.edges().find(|(id, _, _)| *id == edge_id) else {
            return Err(DiagramError::UnknownEdge {
                edge: edge_id.0,
                edges: self.graph.n_edges(),
            });
        };
        if !selected.directed {
            return Ok(self.clone());
        }

        let (particle, _) = Self::resolve_particle(edge_id, selected, model)?;
        let antiparticle = model.antiparticle(particle)?;
        let mut builder = Self::builder(&self.name)
            .symmetry_factor(self.symmetry_factor)
            .overall_factor(&self.overall_factor)
            .numerator_option(self.numerator.clone());
        for (_, vertex) in self.vertices() {
            builder.add_vertex(vertex.clone());
        }
        for (id, endpoints, edge) in self.edges() {
            if id == edge_id {
                let mut reversed = edge.clone();
                reversed.particle =
                    ParticleReference::new(&antiparticle.name, antiparticle.pdg_code);
                builder.add_edge(endpoints.target, endpoints.source, reversed)?;
            } else {
                builder.add_edge(endpoints.source, endpoints.target, edge.clone())?;
            }
        }
        let reversed = builder.build()?;
        reversed.validate(model)?;
        Ok(reversed)
    }

    /// Atomically relabel external indices, leaving unspecified legs unchanged.
    pub fn relabel_external_legs(
        &self,
        relabeling: &BTreeMap<usize, usize>,
    ) -> Result<Self, DiagramError> {
        let external_indices: BTreeSet<_> = self
            .vertices()
            .filter_map(|(_, vertex)| vertex.external.as_ref().map(|leg| leg.index))
            .collect();
        if let Some(index) = relabeling
            .keys()
            .find(|index| !external_indices.contains(index))
        {
            return Err(DiagramError::UnknownExternalIndex(*index));
        }

        self.map_data(
            |_, vertex| {
                let mut vertex = vertex.clone();
                if let Some(external) = &mut vertex.external
                    && let Some(index) = relabeling.get(&external.index)
                {
                    external.index = *index;
                }
                vertex
            },
            |_, _, edge| edge.clone(),
        )
    }

    pub fn to_json(&self) -> Result<String, DiagramError> {
        Ok(serde_json::to_string_pretty(self)?)
    }

    pub fn from_json(input: &str) -> Result<Self, DiagramError> {
        Ok(serde_json::from_str(input)?)
    }

    /// Check model references and the structural invariants expected of a
    /// generated Feynman diagram.
    pub fn validate(&self, model: &Model) -> Result<(), DiagramError> {
        let mut degrees = vec![0_usize; self.graph.n_nodes()];
        for (edge_id, endpoints, edge) in self.edges() {
            Self::resolve_particle(edge_id, edge, model)?;
            degrees[endpoints.source.0] += 1;
            degrees[endpoints.target.0] += 1;
            if self.graph[NodeIndex(endpoints.source.0)].is_external()
                && self.graph[NodeIndex(endpoints.target.0)].is_external()
            {
                return Err(DiagramError::ExternalToExternalEdge { edge: edge_id.0 });
            }
        }
        for (vertex_id, vertex) in self.vertices() {
            if vertex.external.is_some() && degrees[vertex_id.0] != 1 {
                return Err(DiagramError::InvalidExternalDegree {
                    vertex: vertex_id.0,
                    degree: degrees[vertex_id.0],
                });
            }
            if let Some(interaction) = &vertex.interaction {
                let rule = model.vertex_rule(interaction).map_err(|_| {
                    DiagramError::UnknownInteraction {
                        vertex: vertex_id.0,
                        interaction: interaction.clone(),
                    }
                })?;
                let mut actual = Vec::new();
                for (edge_id, endpoints, edge) in self.edges() {
                    let (particle, base) = Self::resolve_particle(edge_id, edge, model)?;
                    if endpoints.source == vertex_id {
                        actual.push((
                            base.pdg_code,
                            edge.directed.then_some(!particle.is_antiparticle()),
                        ));
                    }
                    if endpoints.target == vertex_id {
                        actual.push((
                            base.pdg_code,
                            edge.directed.then_some(particle.is_antiparticle()),
                        ));
                    }
                }
                actual.sort();
                let mut expected = rule
                    .particles
                    .iter()
                    .map(|name| {
                        let particle = model.particle(name)?;
                        let base = if particle.is_antiparticle() {
                            model.antiparticle(particle)?
                        } else {
                            particle
                        };
                        Ok((
                            base.pdg_code,
                            (!particle.is_self_antiparticle())
                                .then_some(particle.is_antiparticle()),
                        ))
                    })
                    .collect::<Result<Vec<_>, feynkit_model::ModelError>>()?;
                expected.sort();
                if actual != expected {
                    return Err(DiagramError::InteractionSignatureMismatch {
                        vertex: vertex_id.0,
                        interaction: interaction.clone(),
                        actual,
                        expected,
                    });
                }
            }
        }
        Ok(())
    }

    /// Serialize to a stable DOT dialect that can be parsed by [`Self::from_dot`].
    pub fn to_dot(&self) -> Result<String, DiagramError> {
        let mut output = String::new();
        let numerator_attribute = self
            .numerator
            .as_ref()
            .map(|numerator| {
                Self::dot_string(numerator).map(|numerator| format!(", numerator={numerator}"))
            })
            .transpose()?
            .unwrap_or_default();
        writeln!(output, "digraph feynkit {{")?;
        writeln!(
            output,
            "  graph [feynkit_name={}, symmetry_factor={}, overall_factor={}{}];",
            Self::dot_string(&self.name)?,
            self.symmetry_factor,
            Self::dot_string(&self.overall_factor)?,
            numerator_attribute,
        )?;

        for (id, vertex) in self.vertices() {
            write!(
                output,
                "  v{} [id={}, feynkit_name={}",
                id.0,
                id.0,
                Self::dot_string(&vertex.name)?
            )?;
            if let Some(rule) = &vertex.interaction {
                write!(output, ", interaction={}", Self::dot_string(rule)?)?;
            }
            if let Some(external) = &vertex.external {
                write!(
                    output,
                    ", external_state={}, external_index={}",
                    Self::dot_string(external.state.as_str())?,
                    external.index
                )?;
            }
            if let Some(numerator) = &vertex.numerator {
                write!(output, ", numerator={}", Self::dot_string(numerator)?)?;
            }
            writeln!(output, "]; ")?;
        }

        for (id, endpoints, edge) in self.edges() {
            write!(
                output,
                "  v{} -> v{} [id={}, pdg={}, particle={}, directed={}",
                endpoints.source.0,
                endpoints.target.0,
                id.0,
                edge.particle.pdg,
                Self::dot_string(&edge.particle.name)?,
                edge.directed
            )?;
            if !edge.directed {
                write!(output, ", dir=none")?;
            }
            if let Some(numerator) = &edge.numerator {
                write!(output, ", numerator={}", Self::dot_string(numerator)?)?;
            }
            writeln!(output, "]; ")?;
        }
        writeln!(output, "}}")?;
        Ok(output)
    }

    /// Parse the stable FeynKit DOT dialect emitted by [`Self::to_dot`].
    pub fn from_dot(input: &str) -> Result<Self, DiagramError> {
        let parsed: DotGraph = DotGraph::from_string(input)
            .map_err(|error| DiagramError::DotParse(error.to_string()))?;
        let name = parsed
            .global_data
            .statements
            .get("feynkit_name")
            .cloned()
            .unwrap_or_else(|| parsed.global_data.name.clone());
        let symmetry_factor = parsed
            .global_data
            .statements
            .get("symmetry_factor")
            .map_or(Ok(1), |value| {
                value
                    .parse()
                    .map_err(|_| DiagramError::InvalidDotAttribute {
                        target: "graph".to_owned(),
                        attribute: "symmetry_factor",
                        value: value.clone(),
                    })
            })?;
        let overall_factor = parsed
            .global_data
            .statements
            .get("overall_factor")
            .cloned()
            .unwrap_or_else(|| "1".to_owned());
        let numerator = parsed.global_data.statements.get("numerator").cloned();

        let mut builder = Self::builder(name)
            .symmetry_factor(symmetry_factor)
            .overall_factor(overall_factor)
            .numerator_option(numerator);
        let mut node_map = BTreeMap::new();
        for (node, _, data) in parsed.iter_nodes() {
            let target = format!("vertex {}", node.0);
            let vertex_name = data
                .statements
                .get("feynkit_name")
                .cloned()
                .or_else(|| data.name.clone())
                .unwrap_or_else(|| format!("v{}", node.0));
            let interaction = data.statements.get("interaction").cloned();
            let external = data
                .statements
                .get("external_state")
                .map(|state| -> Result<ExternalLeg, DiagramError> {
                    let index = data
                        .statements
                        .get("external_index")
                        .ok_or_else(|| DiagramError::MissingDotAttribute {
                            target: target.clone(),
                            attribute: "external_index",
                        })?
                        .parse()
                        .map_err(|_| DiagramError::InvalidDotAttribute {
                            target: target.clone(),
                            attribute: "external_index",
                            value: data.statements["external_index"].clone(),
                        })?;
                    Ok(ExternalLeg {
                        index,
                        state: ExternalState::parse(state)?,
                    })
                })
                .transpose()?;
            let vertex = DiagramVertex {
                name: vertex_name,
                interaction,
                external,
                numerator: data.statements.get("numerator").cloned(),
            };
            node_map.insert(node, builder.add_vertex(vertex));
        }

        let mut next_external = 0;
        for external in builder
            .vertices
            .iter()
            .filter_map(|vertex| vertex.external.as_ref())
        {
            next_external = next_external.max(
                external
                    .index
                    .checked_add(1)
                    .ok_or(DiagramError::DotExternalIndexOverflow(external.index))?,
            );
        }
        for (pair, edge_id, data) in parsed.iter_edges() {
            let target = format!("edge {}", edge_id.0);
            let pdg: i64 = data
                .data
                .statements
                .get("pdg")
                .ok_or_else(|| DiagramError::MissingDotAttribute {
                    target: target.clone(),
                    attribute: "pdg",
                })?
                .parse()
                .map_err(|_| DiagramError::InvalidDotAttribute {
                    target: target.clone(),
                    attribute: "pdg",
                    value: data.data.statements["pdg"].clone(),
                })?;
            let particle_name = data
                .data
                .statements
                .get("particle")
                .cloned()
                .ok_or_else(|| DiagramError::MissingDotAttribute {
                    target: target.clone(),
                    attribute: "particle",
                })?;
            let mut edge = DiagramEdge {
                particle: ParticleReference::new(particle_name, pdg),
                directed: data.orientation != Orientation::Undirected,
                numerator: data.data.statements.get("numerator").cloned(),
            };
            if let Some(directed) = data.data.statements.get("directed") {
                edge.directed =
                    directed
                        .parse()
                        .map_err(|_| DiagramError::InvalidDotAttribute {
                            target: target.clone(),
                            attribute: "directed",
                            value: directed.clone(),
                        })?;
            }

            match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    let source_id = parsed.node_id(source);
                    let sink_id = parsed.node_id(sink);
                    let mut source =
                        *node_map
                            .get(&source_id)
                            .ok_or(DiagramError::UnknownVertex {
                                vertex: source_id.0,
                                vertices: node_map.len(),
                            })?;
                    let mut sink = *node_map.get(&sink_id).ok_or(DiagramError::UnknownVertex {
                        vertex: sink_id.0,
                        vertices: node_map.len(),
                    })?;
                    if data.orientation == Orientation::Reversed {
                        std::mem::swap(&mut source, &mut sink);
                    }
                    builder.add_edge(source, sink, edge)?;
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let internal_id = parsed.node_id(hedge);
                    let internal =
                        *node_map
                            .get(&internal_id)
                            .ok_or(DiagramError::UnknownVertex {
                                vertex: internal_id.0,
                                vertices: node_map.len(),
                            })?;
                    let state = match flow {
                        Flow::Sink => ExternalState::Incoming,
                        Flow::Source => ExternalState::Outgoing,
                    };
                    let external_index = next_external;
                    next_external = next_external
                        .checked_add(1)
                        .ok_or(DiagramError::DotExternalIndexOverflow(next_external))?;
                    let external = builder.add_vertex(DiagramVertex::external(
                        format!("ext{external_index}"),
                        external_index,
                        state,
                    ));
                    match state {
                        ExternalState::Incoming => builder.add_edge(external, internal, edge)?,
                        ExternalState::Outgoing => builder.add_edge(internal, external, edge)?,
                    };
                }
            }
        }
        builder.build()
    }

    pub fn loop_count(&self) -> usize {
        let topology = self.topology();
        topology.internal_edges.len() + topology.components.len() - topology.internal_vertices.len()
    }

    /// Enumerate every spanning-forest-induced loop momentum basis.
    pub fn loop_momentum_bases(&self) -> Result<Vec<LoopMomentumBasis>, DiagramError> {
        self.loop_momentum_bases_with_limit(usize::MAX)
    }

    /// Enumerate at most `limit` loop momentum bases.
    pub fn loop_momentum_bases_with_limit(
        &self,
        limit: usize,
    ) -> Result<Vec<LoopMomentumBasis>, DiagramError> {
        let topology = self.topology();
        let tree_size = topology
            .internal_vertices
            .len()
            .saturating_sub(topology.components.len());
        let mut selections = Vec::new();
        Self::spanning_forests(
            &topology,
            &topology.internal_edges,
            tree_size,
            0,
            &mut Vec::new(),
            &mut selections,
            limit,
        );

        let bases = selections
            .into_iter()
            .map(|tree_edges| topology.basis(tree_edges))
            .collect::<Result<_, _>>()?;
        Ok(bases)
    }

    fn paired_endpoints(&self, pair: HedgePair, orientation: Orientation) -> Option<EdgeEndpoints> {
        let HedgePair::Paired { source, sink } = pair else {
            return None;
        };
        let mut source = VertexId(self.graph.node_id(source).0);
        let mut target = VertexId(self.graph.node_id(sink).0);
        if orientation == Orientation::Reversed {
            std::mem::swap(&mut source, &mut target);
        }
        Some(EdgeEndpoints { source, target })
    }

    fn resolve_particle<'model>(
        edge_id: EdgeId,
        edge: &DiagramEdge,
        model: &'model Model,
    ) -> Result<(&'model Particle, &'model Particle), DiagramError> {
        let particle = model.particle_by_pdg(edge.particle.pdg).map_err(|_| {
            DiagramError::UnknownParticle {
                edge: edge_id.0,
                pdg: edge.particle.pdg,
            }
        })?;
        if particle.name != edge.particle.name {
            return Err(DiagramError::ParticleNameMismatch {
                edge: edge_id.0,
                name: edge.particle.name.clone(),
                pdg: edge.particle.pdg,
                expected_name: particle.name.clone(),
            });
        }
        let base = if particle.is_antiparticle() {
            model.antiparticle(particle)?
        } else {
            particle
        };
        Ok((particle, base))
    }

    fn serde_view(&self) -> DiagramSerde {
        DiagramSerde {
            name: self.name.clone(),
            vertices: self.vertices().map(|(_, vertex)| vertex.clone()).collect(),
            edges: self
                .edges()
                .map(|(_, endpoints, edge)| (endpoints, edge.clone()))
                .collect(),
            symmetry_factor: self.symmetry_factor,
            overall_factor: self.overall_factor.clone(),
            numerator: self.numerator.clone(),
        }
    }

    fn from_serde(data: DiagramSerde) -> Result<Self, DiagramError> {
        let mut builder = Self::builder(data.name)
            .symmetry_factor(data.symmetry_factor)
            .overall_factor(data.overall_factor)
            .numerator_option(data.numerator);
        for vertex in data.vertices {
            builder.add_vertex(vertex);
        }
        for (endpoints, edge) in data.edges {
            builder.add_edge(endpoints.source, endpoints.target, edge)?;
        }
        builder.build()
    }

    fn dot_string(value: &str) -> Result<String, DiagramError> {
        serde_json::to_string(value).map_err(DiagramError::from)
    }

    fn spanning_forests(
        topology: &Topology,
        values: &[EdgeId],
        size: usize,
        start: usize,
        current: &mut Vec<EdgeId>,
        output: &mut Vec<Vec<EdgeId>>,
        limit: usize,
    ) {
        if output.len() >= limit {
            return;
        }
        if current.len() == size {
            if topology.is_spanning_forest(current) {
                output.push(current.clone());
            }
            return;
        }
        let missing = size - current.len();
        for index in start..=values.len().saturating_sub(missing) {
            current.push(values[index]);
            Self::spanning_forests(topology, values, size, index + 1, current, output, limit);
            current.pop();
            if output.len() >= limit {
                break;
            }
        }
    }

    fn topology(&self) -> Topology {
        Topology::new(self)
    }
}

/// Incremental owner for enforcing diagram-construction invariants.
pub struct FeynmanDiagramBuilder {
    name: String,
    vertices: Vec<DiagramVertex>,
    edges: Vec<(EdgeEndpoints, DiagramEdge)>,
    symmetry_factor: u64,
    overall_factor: String,
    numerator: Option<String>,
}

impl FeynmanDiagramBuilder {
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            vertices: Vec::new(),
            edges: Vec::new(),
            symmetry_factor: 1,
            overall_factor: "1".to_owned(),
            numerator: None,
        }
    }

    pub fn symmetry_factor(mut self, symmetry_factor: u64) -> Self {
        self.symmetry_factor = symmetry_factor;
        self
    }

    pub fn overall_factor(mut self, overall_factor: impl Into<String>) -> Self {
        self.overall_factor = overall_factor.into();
        self
    }

    pub fn numerator(mut self, numerator: impl Into<String>) -> Self {
        self.numerator = Some(numerator.into());
        self
    }

    fn numerator_option(mut self, numerator: Option<String>) -> Self {
        self.numerator = numerator;
        self
    }

    pub fn add_vertex(&mut self, vertex: DiagramVertex) -> VertexId {
        let id = VertexId(self.vertices.len());
        self.vertices.push(vertex);
        id
    }

    pub fn add_edge(
        &mut self,
        source: VertexId,
        target: VertexId,
        edge: DiagramEdge,
    ) -> Result<EdgeId, DiagramError> {
        for vertex in [source, target] {
            if vertex.0 >= self.vertices.len() {
                return Err(DiagramError::UnknownVertex {
                    vertex: vertex.0,
                    vertices: self.vertices.len(),
                });
            }
        }
        let id = EdgeId(self.edges.len());
        self.edges.push((EdgeEndpoints { source, target }, edge));
        Ok(id)
    }

    pub fn build(self) -> Result<FeynmanDiagram, DiagramError> {
        let mut external_indices = BTreeSet::new();
        let mut degrees = vec![0_usize; self.vertices.len()];
        for vertex in &self.vertices {
            if let Some(external) = &vertex.external
                && !external_indices.insert(external.index)
            {
                return Err(DiagramError::DuplicateExternalIndex(external.index));
            }
        }
        for (edge, (endpoints, _)) in self.edges.iter().enumerate() {
            degrees[endpoints.source.0] += 1;
            degrees[endpoints.target.0] += 1;
            if self.vertices[endpoints.source.0].is_external()
                && self.vertices[endpoints.target.0].is_external()
            {
                return Err(DiagramError::ExternalToExternalEdge { edge });
            }
        }
        for (vertex, (data, degree)) in self.vertices.iter().zip(degrees).enumerate() {
            if data.is_external() && degree != 1 {
                return Err(DiagramError::InvalidExternalDegree { vertex, degree });
            }
        }

        let mut graph_builder = HedgeGraphBuilder::new();
        let nodes: Vec<_> = self
            .vertices
            .into_iter()
            .map(|vertex| graph_builder.add_node(vertex))
            .collect();
        for (endpoints, edge) in self.edges {
            graph_builder.add_edge(
                nodes[endpoints.source.0],
                nodes[endpoints.target.0],
                edge.clone(),
                edge.directed,
            );
        }
        Ok(FeynmanDiagram {
            name: self.name,
            graph: graph_builder.build(),
            symmetry_factor: self.symmetry_factor,
            overall_factor: self.overall_factor,
            numerator: self.numerator,
        })
    }
}

#[derive(Clone)]
struct TopologyEdge {
    source: VertexId,
    target: VertexId,
}

struct Topology {
    internal_vertices: Vec<VertexId>,
    internal_edges: Vec<EdgeId>,
    external_edges: Vec<EdgeId>,
    edges: BTreeMap<EdgeId, TopologyEdge>,
    components: Vec<Vec<VertexId>>,
    external_state: BTreeMap<EdgeId, ExternalState>,
}

impl Topology {
    fn new(diagram: &FeynmanDiagram) -> Self {
        let internal_vertices: Vec<_> = diagram
            .vertices()
            .filter_map(|(id, vertex)| (!vertex.is_external()).then_some(id))
            .collect();
        let internal_set: BTreeSet<_> = internal_vertices.iter().copied().collect();
        let mut internal_edges = Vec::new();
        let mut external_edges = Vec::new();
        let mut edges = BTreeMap::new();
        let mut external_state = BTreeMap::new();

        for (id, endpoints, _) in diagram.edges() {
            edges.insert(
                id,
                TopologyEdge {
                    source: endpoints.source,
                    target: endpoints.target,
                },
            );
            match (
                internal_set.contains(&endpoints.source),
                internal_set.contains(&endpoints.target),
            ) {
                (true, true) => internal_edges.push(id),
                (true, false) => {
                    external_edges.push(id);
                    if let Some(leg) = diagram
                        .vertex(endpoints.target)
                        .and_then(|v| v.external.as_ref())
                    {
                        external_state.insert(id, leg.state);
                    }
                }
                (false, true) => {
                    external_edges.push(id);
                    if let Some(leg) = diagram
                        .vertex(endpoints.source)
                        .and_then(|v| v.external.as_ref())
                    {
                        external_state.insert(id, leg.state);
                    }
                }
                (false, false) => {}
            }
        }

        let components = Self::components(&internal_vertices, &internal_edges, &edges);
        Self {
            internal_vertices,
            internal_edges,
            external_edges,
            edges,
            components,
            external_state,
        }
    }

    fn components(
        vertices: &[VertexId],
        edges: &[EdgeId],
        edge_data: &BTreeMap<EdgeId, TopologyEdge>,
    ) -> Vec<Vec<VertexId>> {
        let mut remaining: BTreeSet<_> = vertices.iter().copied().collect();
        let mut components = Vec::new();
        while let Some(start) = remaining.pop_first() {
            let mut component = vec![start];
            let mut queue = VecDeque::from([start]);
            while let Some(vertex) = queue.pop_front() {
                for edge_id in edges {
                    let edge = &edge_data[edge_id];
                    let next = if edge.source == vertex {
                        Some(edge.target)
                    } else if edge.target == vertex {
                        Some(edge.source)
                    } else {
                        None
                    };
                    if let Some(next) = next
                        && remaining.remove(&next)
                    {
                        component.push(next);
                        queue.push_back(next);
                    }
                }
            }
            component.sort();
            components.push(component);
        }
        components
    }

    fn is_spanning_forest(&self, selected: &[EdgeId]) -> bool {
        let components = Self::components(&self.internal_vertices, selected, &self.edges);
        components == self.components
    }

    fn basis(&self, tree_edges: Vec<EdgeId>) -> Result<LoopMomentumBasis, DiagramError> {
        let tree_set: BTreeSet<_> = tree_edges.iter().copied().collect();
        let loop_edges: Vec<_> = self
            .internal_edges
            .iter()
            .copied()
            .filter(|edge| !tree_set.contains(edge))
            .collect();
        let mut coefficients: BTreeMap<_, _> = self
            .edges
            .keys()
            .copied()
            .map(|edge| {
                (
                    edge,
                    (
                        vec![0_i8; loop_edges.len()],
                        vec![0_i8; self.external_edges.len()],
                    ),
                )
            })
            .collect();

        for (loop_index, loop_edge) in loop_edges.iter().enumerate() {
            coefficients
                .get_mut(loop_edge)
                .ok_or_else(|| DiagramError::Invariant {
                    operation: "constructing a loop-momentum basis",
                    message: format!("missing loop edge {}", loop_edge.0),
                })?
                .0[loop_index] = 1;
            let edge = self
                .edges
                .get(loop_edge)
                .ok_or_else(|| DiagramError::Invariant {
                    operation: "constructing a loop-momentum basis",
                    message: format!("missing topology edge {}", loop_edge.0),
                })?;
            for (path_edge, sign) in self.path(edge.target, edge.source, &tree_edges) {
                coefficients
                    .get_mut(&path_edge)
                    .ok_or_else(|| DiagramError::Invariant {
                        operation: "constructing a loop-momentum basis",
                        message: format!("missing path edge {}", path_edge.0),
                    })?
                    .0[loop_index] = sign;
            }
        }

        let mut dependent_externals = Vec::new();
        for component in &self.components {
            let component_set: BTreeSet<_> = component.iter().copied().collect();
            let attached: Vec<_> = self
                .external_edges
                .iter()
                .copied()
                .filter(|edge| {
                    let edge = &self.edges[edge];
                    component_set.contains(&edge.source) || component_set.contains(&edge.target)
                })
                .collect();
            let Some(root) = attached.last().copied() else {
                continue;
            };
            dependent_externals.push(root);
            let root_attachment = self.attachment(root);
            let root_flow = self.flow_into_internal(root);
            for edge_id in attached.into_iter().filter(|edge| *edge != root) {
                let external_index = self
                    .external_edges
                    .iter()
                    .position(|edge| *edge == edge_id)
                    .ok_or_else(|| DiagramError::Invariant {
                        operation: "constructing a loop-momentum basis",
                        message: format!("missing external edge {}", edge_id.0),
                    })?;
                let flow = self.flow_into_internal(edge_id);
                coefficients
                    .get_mut(&edge_id)
                    .ok_or_else(|| DiagramError::Invariant {
                        operation: "constructing a loop-momentum basis",
                        message: format!("missing external coefficient edge {}", edge_id.0),
                    })?
                    .1[external_index] = 1;
                coefficients
                    .get_mut(&root)
                    .ok_or_else(|| DiagramError::Invariant {
                        operation: "constructing a loop-momentum basis",
                        message: format!("missing dependent external edge {}", root.0),
                    })?
                    .1[external_index] = -flow * root_flow;
                for (path_edge, sign) in
                    self.path(self.attachment(edge_id), root_attachment, &tree_edges)
                {
                    coefficients
                        .get_mut(&path_edge)
                        .ok_or_else(|| DiagramError::Invariant {
                            operation: "constructing a loop-momentum basis",
                            message: format!("missing external path edge {}", path_edge.0),
                        })?
                        .1[external_index] = flow * sign;
                }
            }
        }

        let edge_signatures = coefficients
            .into_iter()
            .map(|(edge, (loops, external))| {
                Ok((
                    edge,
                    MomentumSignature::new(
                        Signature::try_from_integers(loops).map_err(|error| {
                            DiagramError::Invariant {
                                operation: "constructing a loop-momentum basis",
                                message: error.to_string(),
                            }
                        })?,
                        Signature::try_from_integers(external).map_err(|error| {
                            DiagramError::Invariant {
                                operation: "constructing a loop-momentum basis",
                                message: error.to_string(),
                            }
                        })?,
                    ),
                ))
            })
            .collect::<Result<_, DiagramError>>()?;

        Ok(LoopMomentumBasis {
            tree_edges,
            loop_edges,
            external_edges: self.external_edges.clone(),
            dependent_externals,
            edge_signatures,
        })
    }

    fn attachment(&self, edge: EdgeId) -> VertexId {
        let edge = &self.edges[&edge];
        if self.internal_vertices.contains(&edge.source) {
            edge.source
        } else {
            edge.target
        }
    }

    fn flow_into_internal(&self, edge: EdgeId) -> i8 {
        match self.external_state.get(&edge) {
            Some(ExternalState::Incoming) => 1,
            Some(ExternalState::Outgoing) => -1,
            None => 1,
        }
    }

    fn path(&self, start: VertexId, target: VertexId, tree_edges: &[EdgeId]) -> Vec<(EdgeId, i8)> {
        if start == target {
            return Vec::new();
        }
        let mut queue = VecDeque::from([start]);
        let mut previous: BTreeMap<VertexId, (VertexId, EdgeId, i8)> = BTreeMap::new();
        previous.insert(start, (start, EdgeId(usize::MAX), 0));
        while let Some(vertex) = queue.pop_front() {
            for edge_id in tree_edges {
                let edge = &self.edges[edge_id];
                let step = if edge.source == vertex {
                    Some((edge.target, 1))
                } else if edge.target == vertex {
                    Some((edge.source, -1))
                } else {
                    None
                };
                if let Some((next, sign)) = step
                    && !previous.contains_key(&next)
                {
                    previous.insert(next, (vertex, *edge_id, sign));
                    if next == target {
                        break;
                    }
                    queue.push_back(next);
                }
            }
        }

        let mut path = Vec::new();
        let mut cursor = target;
        while cursor != start {
            let Some((parent, edge, sign)) = previous.get(&cursor).copied() else {
                return Vec::new();
            };
            path.push((edge, sign));
            cursor = parent;
        }
        path.reverse();
        path
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scalar_model() -> Model {
        Model::from_json(
            r#"{
                "name":"phi3","restriction":null,"orders":[],
                "parameters":[
                    {"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal","parameter_type":"real","value":[0.0,0.0],"expression":null},
                    {"name":"M","lhablock":"MASS","lhacode":[25],"nature":"external","parameter_type":"real","value":[1.0,0.0],"expression":null}
                ],
                "particles":[{"pdg_code":25,"name":"phi","antiname":"phi","spin":1,"color":1,"mass":"M","width":"ZERO","texname":"phi","antitexname":"phi","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0}],
                "propagators":[{"name":"phi_prop","particle":"phi","numerator":"1","denominator":"P^2-M^2"}],
                "lorentz_structures":[{"name":"L1","spins":[1,1,1],"structure":"1"}],
                "couplings":[],
                "vertex_rules":[{"name":"V_1","particles":["phi","phi","phi"],"color_structures":["1"],"lorentz_structures":["L1"],"couplings":[[null]]}]
            }"#,
        )
        .unwrap()
    }

    fn one_loop() -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder("bubble");
        let incoming =
            builder.add_vertex(DiagramVertex::external("p1", 0, ExternalState::Incoming));
        let outgoing =
            builder.add_vertex(DiagramVertex::external("p2", 1, ExternalState::Outgoing));
        let left = builder.add_vertex(DiagramVertex::interaction("v0", "V_1"));
        let right = builder.add_vertex(DiagramVertex::interaction("v1", "V_1"));
        let scalar = || DiagramEdge::new(ParticleReference::new("phi", 25), false);
        builder.add_edge(incoming, left, scalar()).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        builder.add_edge(right, outgoing, scalar()).unwrap();
        builder.build().unwrap()
    }

    fn one_loop_reordered() -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder("renamed");
        let right = builder.add_vertex(DiagramVertex::interaction("right", "V_1"));
        let incoming =
            builder.add_vertex(DiagramVertex::external("in", 0, ExternalState::Incoming));
        let left = builder.add_vertex(DiagramVertex::interaction("left", "V_1"));
        let outgoing =
            builder.add_vertex(DiagramVertex::external("out", 1, ExternalState::Outgoing));
        let scalar = || DiagramEdge::new(ParticleReference::new("phi", 25), false);
        builder.add_edge(right, outgoing, scalar()).unwrap();
        builder.add_edge(right, left, scalar()).unwrap();
        builder.add_edge(incoming, left, scalar()).unwrap();
        builder.add_edge(right, left, scalar()).unwrap();
        builder.build().unwrap()
    }

    fn fermion_model() -> Model {
        Model::from_json(
            r#"{
                "name":"fermion","restriction":null,"orders":[],
                "parameters":[
                    {"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal","parameter_type":"real","value":[0.0,0.0],"expression":null},
                    {"name":"M","lhablock":"MASS","lhacode":[1],"nature":"external","parameter_type":"real","value":[1.0,0.0],"expression":null}
                ],
                "particles":[
                    {"pdg_code":1,"name":"f","antiname":"f~","spin":2,"color":1,"mass":"M","width":"ZERO","texname":"f","antitexname":"fbar","charge":0.0,"ghost_number":0,"lepton_number":1,"y_charge":0},
                    {"pdg_code":-1,"name":"f~","antiname":"f","spin":2,"color":1,"mass":"M","width":"ZERO","texname":"fbar","antitexname":"f","charge":0.0,"ghost_number":0,"lepton_number":-1,"y_charge":0}
                ],
                "propagators":[{"name":"f_prop","particle":"f","numerator":"1","denominator":"P-M"}],
                "lorentz_structures":[{"name":"L_f","spins":[2],"structure":"1"}],
                "couplings":[],
                "vertex_rules":[
                    {"name":"V_f","particles":["f"],"color_structures":["1"],"lorentz_structures":["L_f"],"couplings":[[null]]},
                    {"name":"V_af","particles":["f~"],"color_structures":["1"],"lorentz_structures":["L_f"],"couplings":[[null]]}
                ]
            }"#,
        )
        .unwrap()
    }

    fn directed_fermion_line() -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder("fermion-line");
        let anti = builder.add_vertex(DiagramVertex::interaction("anti", "V_af"));
        let particle = builder.add_vertex(DiagramVertex::interaction("particle", "V_f"));
        builder
            .add_edge(
                anti,
                particle,
                DiagramEdge::new(ParticleReference::new("f", 1), true),
            )
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn json_and_dot_round_trip() {
        let diagram = one_loop();
        let from_json = FeynmanDiagram::from_json(&diagram.to_json().unwrap()).unwrap();
        assert_eq!(from_json.loop_count(), 1);
        assert_eq!(from_json.vertices().count(), 4);

        let from_dot = FeynmanDiagram::from_dot(&diagram.to_dot().unwrap()).unwrap();
        assert_eq!(from_dot.name(), "bubble");
        assert_eq!(from_dot.loop_count(), 1);
        assert_eq!(from_dot.edges().count(), 4);
    }

    #[test]
    fn json_import_rejects_invalid_external_structure() {
        let mut missing_edge = serde_json::to_value(one_loop()).unwrap();
        missing_edge["edges"].as_array_mut().unwrap().remove(0);
        assert!(matches!(
            serde_json::from_value::<FeynmanDiagram>(missing_edge),
            Err(error) if error.to_string().contains("external vertex 0 has degree 0")
        ));

        let mut external_edge = serde_json::to_value(one_loop()).unwrap();
        external_edge["edges"][0][0]["target"] = serde_json::json!(1);
        assert!(matches!(
            serde_json::from_value::<FeynmanDiagram>(external_edge),
            Err(error) if error.to_string().contains("edge 0 connects two external vertices")
        ));
    }

    #[test]
    fn dot_import_rejects_external_index_overflow() {
        let dot = format!(
            "digraph feynkit {{ v0 [external_state=\"incoming\", external_index={}]; }}",
            usize::MAX
        );
        assert!(matches!(
            FeynmanDiagram::from_dot(&dot),
            Err(DiagramError::DotExternalIndexOverflow(index)) if index == usize::MAX
        ));
    }

    #[test]
    fn dot_import_requires_canonical_particle_attribute() {
        let dot = r#"digraph feynkit {
            v0 [interaction="V"];
            v1 [interaction="V"];
            v0 -> v1 [pdg=25, name="phi", directed=false];
        }"#;
        assert!(matches!(
            FeynmanDiagram::from_dot(dot),
            Err(DiagramError::MissingDotAttribute {
                target,
                attribute: "particle"
            }) if target == "edge 0"
        ));
    }

    #[test]
    fn validates_particle_identity_and_interaction_signature() {
        let diagram = one_loop();
        diagram.validate(&scalar_model()).unwrap();

        let wrong_name = diagram
            .map_data(
                |_, vertex| vertex.clone(),
                |id, _, edge| {
                    let mut edge = edge.clone();
                    if id == EdgeId(0) {
                        edge.particle.name = "not_phi".to_owned();
                    }
                    edge
                },
            )
            .unwrap();
        assert!(matches!(
            wrong_name.validate(&scalar_model()),
            Err(DiagramError::ParticleNameMismatch { edge: 0, .. })
        ));

        let wrong_flow = diagram
            .map_data(
                |_, vertex| vertex.clone(),
                |id, _, edge| {
                    let mut edge = edge.clone();
                    if id == EdgeId(0) {
                        edge.directed = true;
                    }
                    edge
                },
            )
            .unwrap();
        assert!(matches!(
            wrong_flow.validate(&scalar_model()),
            Err(DiagramError::InteractionSignatureMismatch { .. })
        ));
    }

    #[test]
    fn enumerates_parallel_edge_bases() {
        let diagram = one_loop();
        let bases = diagram.loop_momentum_bases().unwrap();
        assert_eq!(bases.len(), 2);
        assert!(bases.iter().all(|basis| basis.loop_edges.len() == 1));
        assert!(
            bases
                .iter()
                .all(|basis| basis.dependent_externals.len() == 1)
        );
        assert!(bases.iter().all(|basis| {
            basis
                .edge_signatures
                .values()
                .all(|signature| signature.loops.len() == 1 && signature.external.len() == 2)
        }));
    }

    #[test]
    fn basis_limit_counts_spanning_forests_instead_of_candidates() {
        let mut builder = FeynmanDiagram::builder("triangle_with_tail");
        let vertices: Vec<_> = (0..4)
            .map(|index| builder.add_vertex(DiagramVertex::interaction(format!("v{index}"), "V1")))
            .collect();
        let scalar = || DiagramEdge::new(ParticleReference::new("phi", 25), false);
        // The lexicographically first three-edge candidate is the triangle and
        // leaves vertex 3 disconnected. A limit of one must still find the
        // next candidate, which is a valid spanning tree.
        builder
            .add_edge(vertices[0], vertices[1], scalar())
            .unwrap();
        builder
            .add_edge(vertices[1], vertices[2], scalar())
            .unwrap();
        builder
            .add_edge(vertices[2], vertices[0], scalar())
            .unwrap();
        builder
            .add_edge(vertices[2], vertices[3], scalar())
            .unwrap();
        let diagram = builder.build().unwrap();

        let bases = diagram.loop_momentum_bases_with_limit(1).unwrap();
        assert_eq!(bases.len(), 1);
        assert_eq!(bases[0].tree_edges.len(), 3);
        assert_eq!(bases[0].loop_edges.len(), 1);
    }

    #[test]
    fn maps_symbolic_data_without_changing_topology() {
        let diagram = one_loop();
        let mapped = diagram
            .map_data(
                |id, vertex| {
                    let mut vertex = vertex.clone();
                    vertex.name = format!("mapped_{}_{id:?}", vertex.name);
                    vertex
                },
                |id, _, edge| {
                    let mut edge = edge.clone();
                    edge.numerator = Some(format!("N({})", id.0));
                    edge
                },
            )
            .unwrap();

        assert_eq!(mapped.loop_count(), diagram.loop_count());
        assert_eq!(mapped.vertices().count(), diagram.vertices().count());
        assert_eq!(mapped.edges().count(), diagram.edges().count());
        assert!(
            mapped
                .vertices()
                .all(|(_, vertex)| vertex.name.starts_with("mapped_"))
        );
        assert!(mapped.edges().all(|(_, _, edge)| edge.numerator.is_some()));
    }

    #[test]
    fn canonical_key_is_insertion_order_invariant_and_color_sensitive() {
        let model = scalar_model();
        let diagram = one_loop();
        let key = diagram.canonical_key(&model).unwrap();
        assert_eq!(key, one_loop_reordered().canonical_key(&model).unwrap());

        let changed_state = diagram
            .map_data(
                |_, vertex| {
                    let mut vertex = vertex.clone();
                    if let Some(external) = &mut vertex.external
                        && external.index == 0
                    {
                        external.state = ExternalState::Outgoing;
                    }
                    vertex
                },
                |_, _, edge| edge.clone(),
            )
            .unwrap();
        assert_ne!(key, changed_state.canonical_key(&model).unwrap());

        let changed_interaction = diagram
            .map_data(
                |id, vertex| {
                    let mut vertex = vertex.clone();
                    if id == VertexId(2) {
                        vertex.interaction = None;
                    }
                    vertex
                },
                |_, _, edge| edge.clone(),
            )
            .unwrap();
        assert_ne!(key, changed_interaction.canonical_key(&model).unwrap());
    }

    #[test]
    fn reverses_particle_flow_with_antiparticle_semantics() {
        let model = fermion_model();
        let diagram = directed_fermion_line();
        diagram.validate(&model).unwrap();
        let original_json = diagram.to_json().unwrap();
        let original_key = diagram.canonical_key(&model).unwrap();

        let reversed = diagram.reverse_edge(EdgeId(0), &model).unwrap();
        let (_, endpoints, edge) = reversed.edges().next().unwrap();
        assert_eq!(endpoints.source, VertexId(1));
        assert_eq!(endpoints.target, VertexId(0));
        assert_eq!(edge.particle, ParticleReference::new("f~", -1));
        assert_eq!(reversed.canonical_key(&model).unwrap(), original_key);

        let restored = reversed.reverse_edge(EdgeId(0), &model).unwrap();
        assert_eq!(restored.to_json().unwrap(), original_json);
        assert!(matches!(
            diagram.reverse_edge(EdgeId(1), &model),
            Err(DiagramError::UnknownEdge { edge: 1, edges: 1 })
        ));
    }

    #[test]
    fn relabels_external_legs_atomically_and_round_trips() {
        let diagram = one_loop();
        let swap = BTreeMap::from([(0, 1), (1, 0)]);
        let relabeled = diagram.relabel_external_legs(&swap).unwrap();
        assert_eq!(
            relabeled
                .vertex(VertexId(0))
                .unwrap()
                .external
                .as_ref()
                .unwrap()
                .index,
            1
        );
        assert_eq!(
            relabeled
                .vertex(VertexId(1))
                .unwrap()
                .external
                .as_ref()
                .unwrap()
                .index,
            0
        );
        assert_eq!(
            relabeled
                .relabel_external_legs(&swap)
                .unwrap()
                .to_json()
                .unwrap(),
            diagram.to_json().unwrap()
        );

        assert!(matches!(
            diagram.relabel_external_legs(&BTreeMap::from([(2, 3)])),
            Err(DiagramError::UnknownExternalIndex(2))
        ));
        assert!(matches!(
            diagram.relabel_external_legs(&BTreeMap::from([(0, 1)])),
            Err(DiagramError::DuplicateExternalIndex(1))
        ));
    }

    #[test]
    fn rejects_duplicate_external_indices() {
        let mut builder = FeynmanDiagram::builder("invalid");
        builder.add_vertex(DiagramVertex::external("p1", 0, ExternalState::Incoming));
        builder.add_vertex(DiagramVertex::external("p2", 0, ExternalState::Outgoing));
        assert!(matches!(
            builder.build(),
            Err(DiagramError::DuplicateExternalIndex(0))
        ));
    }
}
