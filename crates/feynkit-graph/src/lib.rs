//! Lightweight, application-independent Feynman diagrams.
//!
//! This crate owns the structural and symbolic diagram layer, including the
//! physical cut partitions selected during canonical cross-section generation.
//! Runtime cut evaluators, subtraction data, parametrizations, and caches belong
//! to applications consuming these diagrams.

#![forbid(unsafe_code)]

mod display;

use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    fmt::{self, Write},
    str::FromStr,
    sync::Arc,
};

pub use feynkit_kinematics::MomentumSignature;
use feynkit_kinematics::Signature;
use feynkit_model::{Model, ModelError, ModelFingerprint, Particle, ParticleId, VertexRuleId};
use linnet::{
    half_edge::{
        HedgeGraph, NodeIndex,
        builder::HedgeGraphBuilder,
        involution::{Flow, HedgePair, Orientation},
    },
    parser::DotGraph,
};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    graph::Graph as CanonicalGraph,
    parser::ParseSettings,
};
use thiserror::Error;

/// Stable vertex identifier in a [`FeynmanDiagram`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct VertexId(pub usize);

/// Stable edge identifier in a [`FeynmanDiagram`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct EdgeId(pub usize);

/// Stable interaction-leg position of one edge endpoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct VertexSlot(pub usize);

/// Stable content-derived identifier for a finalized diagram topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct DiagramId(pub u128);

impl DiagramId {
    fn from_key(model: ModelFingerprint, key: &CanonicalDiagramKey) -> Result<Self, DiagramError> {
        let bytes = serde_json::to_vec(&(model, key))?;
        fn fnv(bytes: &[u8], offset: u64) -> u64 {
            bytes.iter().fold(offset, |hash, byte| {
                (hash ^ u64::from(*byte)).wrapping_mul(0x0000_0100_0000_01b3)
            })
        }
        let high = fnv(&bytes, 0xcbf2_9ce4_8422_2325);
        let low = fnv(&bytes, 0x8422_2325_cbf2_9ce4);
        Ok(Self((u128::from(high) << 64) | u128::from(low)))
    }
}

impl fmt::Display for DiagramId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{:032x}", self.0)
    }
}

impl FromStr for DiagramId {
    type Err = std::num::ParseIntError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        u128::from_str_radix(value, 16).map(Self)
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
    /// Unique momentum label of this external leg.
    pub index: usize,
    pub state: ExternalState,
    /// Sewing identity. Amplitude legs use distinct values; the incoming and
    /// outgoing copies of a cross-section initial state share one value.
    pub connection: usize,
}

/// Structural and symbolic information associated with a vertex.
#[derive(Debug, Clone)]
pub struct DiagramVertex {
    pub name: String,
    pub interaction: Option<VertexRuleId>,
    pub external: Option<ExternalLeg>,
    pub numerator: Atom,
}

impl DiagramVertex {
    pub fn interaction(name: impl Into<String>, rule: VertexRuleId) -> Self {
        Self {
            name: name.into(),
            interaction: Some(rule),
            external: None,
            numerator: Atom::one(),
        }
    }

    pub fn external(name: impl Into<String>, index: usize, state: ExternalState) -> Self {
        Self::external_in_connection(name, index, state, index)
    }

    pub fn external_in_connection(
        name: impl Into<String>,
        index: usize,
        state: ExternalState,
        connection: usize,
    ) -> Self {
        Self {
            name: name.into(),
            interaction: None,
            external: Some(ExternalLeg {
                index,
                state,
                connection,
            }),
            numerator: Atom::one(),
        }
    }

    pub fn is_external(&self) -> bool {
        self.external.is_some()
    }
}

/// Structural and symbolic information associated with an edge.
#[derive(Debug, Clone)]
pub struct DiagramEdge {
    pub particle: ParticleId,
    pub directed: bool,
    pub numerator: Atom,
    source_slot: VertexSlot,
    target_slot: VertexSlot,
}

impl DiagramEdge {
    pub fn new(particle: ParticleId, directed: bool) -> Self {
        Self {
            particle,
            directed,
            numerator: Atom::one(),
            source_slot: VertexSlot(0),
            target_slot: VertexSlot(0),
        }
    }

    pub fn source_slot(&self) -> VertexSlot {
        self.source_slot
    }

    pub fn target_slot(&self) -> VertexSlot {
        self.target_slot
    }

    pub fn with_slots(mut self, source: VertexSlot, target: VertexSlot) -> Self {
        self.source_slot = source;
        self.target_slot = target;
        self
    }
}

/// Source and target in the canonical orientation of an edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct EdgeEndpoints {
    pub source: VertexId,
    pub target: VertexId,
}

/// One stable half-edge endpoint in a finalized diagram.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DiagramEndpoint {
    Source,
    Target,
}

/// A stable half-edge reference used by canonical cut metadata.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct DiagramHalfEdge {
    pub edge: EdgeId,
    pub endpoint: DiagramEndpoint,
}

/// One amplitude side retained for a physical cross-section cut.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct DiagramCutSide {
    pub half_edges: Vec<DiagramHalfEdge>,
    pub coupling_orders: BTreeMap<String, usize>,
    pub loop_count: usize,
}

/// A physical cross-section cut selected by canonical FeynKit generation.
///
/// `cut` contains the half-edge on the left side of every oriented cut edge;
/// the opposite endpoints are implied. `left` and `right` preserve the exact
/// amplitude partitions used by generation filters.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct DiagramCut {
    pub cut: Vec<DiagramHalfEdge>,
    pub left: DiagramCutSide,
    pub right: DiagramCutSide,
}

/// Side label used by the auxiliary cut-incidence graph in a canonical key.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CanonicalCutSide {
    Left,
    Right,
}

/// Vertex color retained by [`CanonicalDiagramKey`].
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CanonicalDiagramVertex {
    Diagram {
        interaction: Option<VertexRuleId>,
        external: Option<ExternalLeg>,
    },
    CutSide {
        side: CanonicalCutSide,
        coupling_orders: BTreeMap<String, usize>,
        loop_count: usize,
    },
}

/// Canonically labeled edge retained by [`CanonicalDiagramKey`].
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CanonicalDiagramEdge {
    Diagram {
        endpoints: EdgeEndpoints,
        particle: ParticleId,
        directed: bool,
        source_slot: VertexSlot,
        target_slot: VertexSlot,
    },
    CutPair {
        endpoints: EdgeEndpoints,
    },
    CutMembership {
        endpoints: EdgeEndpoints,
    },
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

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum CanonicalEdgeColor {
    Diagram {
        particle: ParticleId,
        source_slot: VertexSlot,
        target_slot: VertexSlot,
    },
    CutPair,
    CutMembership,
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

impl LoopMomentumBasis {
    /// Validate that this basis is a complete routing of `diagram`.
    ///
    /// Every edge must occur in exactly one structural role, the selected tree
    /// must be a spanning forest, and every stored signature must agree with
    /// the routing induced by that forest and the requested loop-edge order.
    pub fn validate(&self, diagram: &FeynmanDiagram) -> Result<(), DiagramError> {
        let topology = diagram.topology();
        let unique = |edges: &[EdgeId]| edges.iter().copied().collect::<BTreeSet<_>>();
        let tree = unique(&self.tree_edges);
        let loops = unique(&self.loop_edges);
        let external = unique(&self.external_edges);
        let dependent = unique(&self.dependent_externals);
        if tree.len() != self.tree_edges.len()
            || loops.len() != self.loop_edges.len()
            || external.len() != self.external_edges.len()
            || dependent.len() != self.dependent_externals.len()
        {
            return Err(DiagramError::InvalidLoopMomentumBasis(
                "edge roles contain duplicates".to_owned(),
            ));
        }

        let internal = unique(&topology.internal_edges);
        let expected_external = unique(&topology.external_edges);
        if !tree.is_subset(&internal)
            || !loops.is_subset(&internal)
            || !tree.is_disjoint(&loops)
            || tree.union(&loops).copied().collect::<BTreeSet<_>>() != internal
            || external != expected_external
            || !dependent.is_subset(&external)
            || !topology.is_spanning_forest(&self.tree_edges)
        {
            return Err(DiagramError::InvalidLoopMomentumBasis(
                "tree, loop, external, or dependent edge roles do not match the diagram topology"
                    .to_owned(),
            ));
        }

        let expected = topology
            .basis(self.tree_edges.clone())?
            .with_loop_edge_order(&self.loop_edges)?;
        if *self != expected {
            return Err(DiagramError::InvalidLoopMomentumBasis(
                "stored dependent externals or momentum signatures do not match the selected spanning forest"
                    .to_owned(),
            ));
        }
        Ok(())
    }

    /// Return the same routing with its independent loop momenta ordered by edge.
    pub fn with_loop_edge_order(mut self, requested: &[EdgeId]) -> Result<Self, DiagramError> {
        let current: BTreeSet<_> = self.loop_edges.iter().copied().collect();
        let requested_set: BTreeSet<_> = requested.iter().copied().collect();
        if requested.len() != requested_set.len() || requested_set != current {
            return Err(DiagramError::InvalidLoopMomentumEdges {
                requested: requested.to_vec(),
                available: self.loop_edges,
            });
        }

        let positions = requested
            .iter()
            .map(|edge| {
                self.loop_edges
                    .iter()
                    .position(|candidate| candidate == edge)
                    .ok_or_else(|| {
                        DiagramError::InvalidLoopMomentumBasis(format!(
                            "requested loop edge {} is not present in the stored basis",
                            edge.0
                        ))
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;
        for signature in self.edge_signatures.values_mut() {
            let loops = positions
                .iter()
                .map(|position| {
                    signature.loops.get(*position).ok_or_else(|| {
                        DiagramError::InvalidLoopMomentumBasis(format!(
                            "a momentum signature has {} loop entries, expected {}",
                            signature.loops.len(),
                            self.loop_edges.len()
                        ))
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            signature.loops = Signature::new(loops);
        }
        self.loop_edges = requested.to_vec();
        Ok(self)
    }
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
    #[error("external connection {connection} contains {legs} legs, expected one or two")]
    InvalidExternalConnectionSize { connection: usize, legs: usize },
    #[error("external connection {connection} contains two {state} legs")]
    InvalidExternalConnectionStates {
        connection: usize,
        state: &'static str,
    },
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
    #[error("failed to parse diagram {field} '{expression}': {message}")]
    SymbolicParse {
        field: &'static str,
        expression: String,
        message: String,
    },
    #[error("diagram model fingerprint {serialized} does not match supplied model {actual}")]
    ModelFingerprintMismatch {
        serialized: ModelFingerprint,
        actual: ModelFingerprint,
    },
    #[error(
        "vertex {vertex} interaction {interaction:?} has incident signature {actual:?}, expected {expected:?}"
    )]
    InteractionSignatureMismatch {
        vertex: usize,
        interaction: VertexRuleId,
        actual: Vec<(i64, Option<bool>)>,
        expected: Vec<(i64, Option<bool>)>,
    },
    #[error("external vertex {vertex} has degree {degree}, expected one")]
    InvalidExternalDegree { vertex: usize, degree: usize },
    #[error("edge {edge} connects two external vertices")]
    ExternalToExternalEdge { edge: usize },
    #[error("vertex {vertex} uses endpoint slot {slot} more than once")]
    DuplicateVertexSlot { vertex: usize, slot: usize },
    #[error("vertex {vertex} endpoint slots are {actual:?}, expected {expected:?}")]
    InvalidVertexSlots {
        vertex: usize,
        actual: Vec<usize>,
        expected: Vec<usize>,
    },
    #[error("external vertex {vertex} must have unit numerator")]
    ExternalVertexNumerator { vertex: usize },
    #[error("edge {edge} incident to an external vertex must have unit numerator")]
    ExternalEdgeNumerator { edge: usize },
    #[error("cross-section external connections require at least one finalized physical cut")]
    MissingCrossSectionCuts,
    #[error("invalid physical cut {cut}: {message}")]
    InvalidCut { cut: usize, message: String },
    #[error(
        "diagram-wide numerator differs from the product of the vertex and edge numerator fragments"
    )]
    NumeratorFragmentMismatch,
    #[error("diagram invariant failed while {operation}: {message}")]
    Invariant {
        operation: &'static str,
        message: String,
    },
    #[error("failed to format DOT")]
    Formatting(#[from] std::fmt::Error),
    #[error("the diagram has no loop-momentum basis")]
    MissingLoopMomentumBasis,
    #[error("invalid loop-momentum basis: {0}")]
    InvalidLoopMomentumBasis(String),
    #[error("serialized diagram ID {serialized} does not match its content-derived ID {actual}")]
    DiagramIdMismatch {
        serialized: DiagramId,
        actual: DiagramId,
    },
    #[error(
        "requested loop-momentum edges {requested:?} do not select the available basis {available:?}"
    )]
    InvalidLoopMomentumEdges {
        requested: Vec<EdgeId>,
        available: Vec<EdgeId>,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DiagramSerde {
    id: DiagramId,
    model: ModelFingerprint,
    name: String,
    vertices: Vec<DiagramVertexSerde>,
    edges: Vec<(EdgeEndpoints, DiagramEdgeSerde)>,
    symmetry_factor: u64,
    overall_factor: String,
    numerator: String,
    numerator_prefactor: String,
    projector: String,
    loop_momentum_basis: LoopMomentumBasis,
    cuts: Vec<DiagramCut>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DiagramVertexSerde {
    name: String,
    interaction: Option<VertexRuleId>,
    external: Option<ExternalLeg>,
    numerator: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DiagramEdgeSerde {
    particle: ParticleId,
    directed: bool,
    numerator: String,
    source_slot: VertexSlot,
    target_slot: VertexSlot,
}

/// Model-aware Feynman-diagram IR backed by Linnet's half-edge graph.
#[derive(Debug, Clone)]
pub struct FeynmanDiagram {
    model: Arc<Model>,
    id: DiagramId,
    name: String,
    graph: HedgeGraph<DiagramEdge, DiagramVertex>,
    symmetry_factor: u64,
    overall_factor: Atom,
    numerator: Atom,
    numerator_prefactor: Atom,
    projector: Atom,
    loop_momentum_basis: LoopMomentumBasis,
    cuts: Vec<DiagramCut>,
}

impl Serialize for FeynmanDiagram {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.serde_view().serialize(serializer)
    }
}

impl FeynmanDiagram {
    pub fn builder(model: impl Into<Arc<Model>>, name: impl Into<String>) -> FeynmanDiagramBuilder {
        FeynmanDiagramBuilder::new(model, name)
    }

    pub fn model(&self) -> &Model {
        &self.model
    }

    pub fn model_arc(&self) -> Arc<Model> {
        Arc::clone(&self.model)
    }

    pub fn id(&self) -> DiagramId {
        self.id
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    /// Return the raw automorphism order retained as generation provenance.
    ///
    /// This value is diagnostic only. [`Self::overall_factor`] is the
    /// authoritative symbolic weight and may also include grouping or fermion
    /// signs, so runtime consumers must not apply the symmetry factor again.
    pub fn symmetry_factor(&self) -> u64 {
        self.symmetry_factor
    }

    pub fn overall_factor(&self) -> &Atom {
        &self.overall_factor
    }

    pub fn numerator(&self) -> &Atom {
        &self.numerator
    }

    /// Global numerator multiplier supplied by the generation request.
    pub fn numerator_prefactor(&self) -> &Atom {
        &self.numerator_prefactor
    }

    /// External-state projector applied to the generated numerator.
    pub fn projector(&self) -> &Atom {
        &self.projector
    }

    /// The routing selected during diagram finalization.
    pub fn loop_momentum_basis(&self) -> &LoopMomentumBasis {
        &self.loop_momentum_basis
    }

    /// Physical cut partitions retained by cross-section generation.
    pub fn cuts(&self) -> &[DiagramCut] {
        &self.cuts
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
        let loop_edges = self.loop_momentum_basis.loop_edges.clone();
        let mut builder = Self::builder(Arc::clone(&self.model), &self.name)
            .symmetry_factor(self.symmetry_factor)
            .overall_factor(self.overall_factor.clone())
            .numerator(self.numerator.clone())
            .numerator_prefactor(self.numerator_prefactor.clone())
            .projector(self.projector.clone())
            .cuts(self.cuts.clone());
        for (id, vertex) in self.vertices() {
            builder.add_vertex(map_vertex(id, vertex));
        }
        for (id, endpoints, edge) in self.edges() {
            let edge = map_edge(id, endpoints, edge);
            builder.add_edge_with_slots(
                endpoints.source,
                endpoints.target,
                edge.clone(),
                edge.source_slot,
                edge.target_slot,
            )?;
        }
        let diagram = builder.build()?;
        if diagram.loop_momentum_basis.loop_edges == loop_edges {
            Ok(diagram)
        } else {
            diagram.with_loop_momentum_edges(&loop_edges)
        }
    }

    /// Build a deterministic key for the model-resolved colored topology.
    ///
    /// Canonical labeling is delegated to Symbolica's graph implementation.
    /// Particle/antiparticle edge representations are normalized so reversing
    /// a directed edge with [`Self::reverse_edge`] does not change the key.
    pub fn canonical_key(&self) -> Result<CanonicalDiagramKey, DiagramError> {
        self.validate()?;
        self.canonical_key_impl(true)
    }

    fn structural_key(&self) -> Result<CanonicalDiagramKey, DiagramError> {
        self.canonical_key_impl(false)
    }

    fn canonical_key_impl(
        &self,
        normalize_particle_orientation: bool,
    ) -> Result<CanonicalDiagramKey, DiagramError> {
        let mut graph = CanonicalGraph::new();
        for (_, vertex) in self.vertices() {
            graph.add_node(CanonicalDiagramVertex::Diagram {
                interaction: vertex.interaction,
                external: vertex.external.clone(),
            });
        }
        let endpoints = self
            .edges()
            .map(|(edge, endpoints, _)| (edge, endpoints))
            .collect::<BTreeMap<_, _>>();
        for cut in &self.cuts {
            let left = graph.add_node(CanonicalDiagramVertex::CutSide {
                side: CanonicalCutSide::Left,
                coupling_orders: cut.left.coupling_orders.clone(),
                loop_count: cut.left.loop_count,
            });
            let right = graph.add_node(CanonicalDiagramVertex::CutSide {
                side: CanonicalCutSide::Right,
                coupling_orders: cut.right.coupling_orders.clone(),
                loop_count: cut.right.loop_count,
            });
            graph
                .add_edge(left, right, false, CanonicalEdgeColor::CutPair)
                .map_err(|error| DiagramError::Invariant {
                    operation: "deriving the diagram ID",
                    message: error.to_string(),
                })?;
            for (side, half_edges) in [(left, &cut.left.half_edges), (right, &cut.right.half_edges)]
            {
                let vertices = half_edges
                    .iter()
                    .filter_map(|half_edge| {
                        endpoints
                            .get(&half_edge.edge)
                            .map(|endpoints| match half_edge.endpoint {
                                DiagramEndpoint::Source => endpoints.source,
                                DiagramEndpoint::Target => endpoints.target,
                            })
                    })
                    .collect::<BTreeSet<_>>();
                for vertex in vertices {
                    graph
                        .add_edge(side, vertex.0, false, CanonicalEdgeColor::CutMembership)
                        .map_err(|error| DiagramError::Invariant {
                            operation: "encoding canonical cut membership",
                            message: error.to_string(),
                        })?;
                }
            }
        }
        for (edge_id, endpoints, edge) in self.edges() {
            let mut particle = edge.particle;
            let (mut source, mut target) = (endpoints.source.0, endpoints.target.0);
            let (mut source_slot, mut target_slot) = (edge.source_slot, edge.target_slot);
            if normalize_particle_orientation {
                let (resolved, base) = self.resolve_particle(edge_id, edge)?;
                particle = base;
                if edge.directed && resolved.is_antiparticle() {
                    std::mem::swap(&mut source, &mut target);
                    std::mem::swap(&mut source_slot, &mut target_slot);
                }
            }
            graph
                .add_edge(
                    source,
                    target,
                    edge.directed,
                    CanonicalEdgeColor::Diagram {
                        particle,
                        source_slot,
                        target_slot,
                    },
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
                .map(|edge| {
                    let endpoints = EdgeEndpoints {
                        source: VertexId(edge.vertices.0),
                        target: VertexId(edge.vertices.1),
                    };
                    match edge.data {
                        CanonicalEdgeColor::Diagram {
                            particle,
                            source_slot,
                            target_slot,
                        } => CanonicalDiagramEdge::Diagram {
                            endpoints,
                            particle,
                            directed: edge.directed,
                            source_slot,
                            target_slot,
                        },
                        CanonicalEdgeColor::CutPair => CanonicalDiagramEdge::CutPair { endpoints },
                        CanonicalEdgeColor::CutMembership => {
                            CanonicalDiagramEdge::CutMembership { endpoints }
                        }
                    }
                })
                .collect(),
        })
    }

    /// Reverse a directed edge while preserving its physical particle flow.
    ///
    /// Reversal swaps the endpoints and replaces the referenced particle with
    /// its antiparticle. Undirected edges are unchanged. The returned diagram
    /// is fully model-validated and retains all stable identifiers.
    pub fn reverse_edge(&self, edge_id: EdgeId) -> Result<Self, DiagramError> {
        self.validate()?;
        let Some((_, _, selected)) = self.edges().find(|(id, _, _)| *id == edge_id) else {
            return Err(DiagramError::UnknownEdge {
                edge: edge_id.0,
                edges: self.graph.n_edges(),
            });
        };
        if !selected.directed {
            return Ok(self.clone());
        }

        let (particle, _) = self.resolve_particle(edge_id, selected)?;
        let antiparticle = self.model.antiparticle(particle)?;
        let antiparticle = self.model.particle_id(&antiparticle.name)?;
        let reverse_half_edge = |half_edge: &mut DiagramHalfEdge| {
            if half_edge.edge == edge_id {
                half_edge.endpoint = match half_edge.endpoint {
                    DiagramEndpoint::Source => DiagramEndpoint::Target,
                    DiagramEndpoint::Target => DiagramEndpoint::Source,
                };
            }
        };
        let mut cuts = self.cuts.clone();
        for cut in &mut cuts {
            for half_edge in &mut cut.cut {
                reverse_half_edge(half_edge);
            }
            for half_edge in &mut cut.left.half_edges {
                reverse_half_edge(half_edge);
            }
            for half_edge in &mut cut.right.half_edges {
                reverse_half_edge(half_edge);
            }
        }
        let mut builder = Self::builder(Arc::clone(&self.model), &self.name)
            .symmetry_factor(self.symmetry_factor)
            .overall_factor(self.overall_factor.clone())
            .numerator(self.numerator.clone())
            .numerator_prefactor(self.numerator_prefactor.clone())
            .projector(self.projector.clone())
            .cuts(cuts);
        for (_, vertex) in self.vertices() {
            builder.add_vertex(vertex.clone());
        }
        for (id, endpoints, edge) in self.edges() {
            if id == edge_id {
                let mut reversed = edge.clone();
                reversed.particle = antiparticle;
                builder.add_edge_with_slots(
                    endpoints.target,
                    endpoints.source,
                    reversed,
                    edge.target_slot,
                    edge.source_slot,
                )?;
            } else {
                builder.add_edge_with_slots(
                    endpoints.source,
                    endpoints.target,
                    edge.clone(),
                    edge.source_slot,
                    edge.target_slot,
                )?;
            }
        }
        let reversed = builder
            .build()?
            .with_loop_momentum_edges(&self.loop_momentum_basis.loop_edges)?;
        reversed.validate()?;
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

    /// Return this diagram with a finalized deterministic name.
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = name.into();
        self
    }

    /// Return this diagram with a finalized numerator-independent multiplier.
    pub fn with_overall_factor(mut self, factor: Atom) -> Self {
        self.overall_factor = factor;
        self
    }

    /// Replace the finalized physical cuts and refresh the content-derived ID.
    pub fn with_cuts(mut self, mut cuts: Vec<DiagramCut>) -> Result<Self, DiagramError> {
        for cut in &mut cuts {
            cut.cut.sort();
            cut.cut.dedup();
            cut.left.half_edges.sort();
            cut.left.half_edges.dedup();
            cut.right.half_edges.sort();
            cut.right.half_edges.dedup();
        }
        cuts.sort();
        cuts.dedup();
        self.cuts = cuts;
        self.validate()?;
        self.id = DiagramId::from_key(self.model.fingerprint(), &self.structural_key()?)?;
        Ok(self)
    }

    /// Select a spanning-forest routing by its ordered independent edges.
    pub fn with_loop_momentum_edges(mut self, requested: &[EdgeId]) -> Result<Self, DiagramError> {
        let requested_set: BTreeSet<_> = requested.iter().copied().collect();
        let basis = self
            .loop_momentum_bases()?
            .into_iter()
            .find(|basis| {
                basis.loop_edges.len() == requested.len()
                    && basis
                        .loop_edges
                        .iter()
                        .all(|edge| requested_set.contains(edge))
            })
            .ok_or_else(|| DiagramError::InvalidLoopMomentumEdges {
                requested: requested.to_vec(),
                available: self.loop_momentum_basis.loop_edges.clone(),
            })?
            .with_loop_edge_order(requested)?;
        self.loop_momentum_basis = basis;
        Ok(self)
    }

    pub fn to_json(&self) -> Result<String, DiagramError> {
        Ok(serde_json::to_string_pretty(&self.serde_view())?)
    }

    pub fn from_json(model: impl Into<Arc<Model>>, input: &str) -> Result<Self, DiagramError> {
        Self::from_serde(model.into(), serde_json::from_str(input)?)
    }

    /// Check model references and the structural invariants expected of a
    /// generated Feynman diagram.
    pub fn validate(&self) -> Result<(), DiagramError> {
        let mut degrees = vec![0_usize; self.graph.n_nodes()];
        let mut fragment_numerator = Atom::one();
        let mut external_connections = BTreeMap::<usize, Vec<ExternalState>>::new();
        let mut vertex_slots = vec![BTreeSet::new(); self.graph.n_nodes()];
        let mut vertex_half_edges = vec![Vec::new(); self.graph.n_nodes()];
        for (edge_id, endpoints, edge) in self.edges() {
            self.resolve_particle(edge_id, edge)?;
            degrees[endpoints.source.0] += 1;
            degrees[endpoints.target.0] += 1;
            let source_external = self.graph[NodeIndex(endpoints.source.0)].is_external();
            let target_external = self.graph[NodeIndex(endpoints.target.0)].is_external();
            if source_external && target_external {
                return Err(DiagramError::ExternalToExternalEdge { edge: edge_id.0 });
            }
            if (source_external || target_external) && edge.numerator != Atom::one() {
                return Err(DiagramError::ExternalEdgeNumerator { edge: edge_id.0 });
            }
            for (vertex, slot) in [
                (endpoints.source, edge.source_slot),
                (endpoints.target, edge.target_slot),
            ] {
                if !vertex_slots[vertex.0].insert(slot.0) {
                    return Err(DiagramError::DuplicateVertexSlot {
                        vertex: vertex.0,
                        slot: slot.0,
                    });
                }
            }
            vertex_half_edges[endpoints.source.0].push(DiagramHalfEdge {
                edge: edge_id,
                endpoint: DiagramEndpoint::Source,
            });
            vertex_half_edges[endpoints.target.0].push(DiagramHalfEdge {
                edge: edge_id,
                endpoint: DiagramEndpoint::Target,
            });
            fragment_numerator *= &edge.numerator;
        }
        for (vertex_id, vertex) in self.vertices() {
            if let Some(external) = &vertex.external {
                external_connections
                    .entry(external.connection)
                    .or_default()
                    .push(external.state);
            }
            if vertex.external.is_some() && degrees[vertex_id.0] != 1 {
                return Err(DiagramError::InvalidExternalDegree {
                    vertex: vertex_id.0,
                    degree: degrees[vertex_id.0],
                });
            }
            let actual = vertex_slots[vertex_id.0]
                .iter()
                .copied()
                .collect::<Vec<_>>();
            let expected = (0..degrees[vertex_id.0]).collect::<Vec<_>>();
            if actual != expected {
                return Err(DiagramError::InvalidVertexSlots {
                    vertex: vertex_id.0,
                    actual,
                    expected,
                });
            }
            if vertex.external.is_some() && vertex.numerator != Atom::one() {
                return Err(DiagramError::ExternalVertexNumerator {
                    vertex: vertex_id.0,
                });
            }
            fragment_numerator *= &vertex.numerator;
            if let Some(interaction) = vertex.interaction {
                let rule = self.model.vertex_rule_by_id(interaction)?;
                let mut actual = vec![None; degrees[vertex_id.0]];
                for (edge_id, endpoints, edge) in self.edges() {
                    let (particle, base_id) = self.resolve_particle(edge_id, edge)?;
                    let base = self.model.particle_by_id(base_id)?;
                    if endpoints.source == vertex_id {
                        actual[edge.source_slot.0] = Some((
                            base.pdg_code,
                            edge.directed.then_some(!particle.is_antiparticle()),
                        ));
                    }
                    if endpoints.target == vertex_id {
                        actual[edge.target_slot.0] = Some((
                            base.pdg_code,
                            edge.directed.then_some(particle.is_antiparticle()),
                        ));
                    }
                }
                let actual = actual.into_iter().flatten().collect::<Vec<_>>();
                let expected = rule
                    .particles
                    .iter()
                    .map(|particle_id| {
                        let particle = self.model.particle_by_id(*particle_id)?;
                        let base = if particle.is_antiparticle() {
                            self.model.antiparticle(particle)?
                        } else {
                            particle
                        };
                        Ok((
                            base.pdg_code,
                            (particle.antiparticle != *particle_id)
                                .then_some(particle.is_antiparticle()),
                        ))
                    })
                    .collect::<Result<Vec<_>, feynkit_model::ModelError>>()?;
                if actual != expected {
                    return Err(DiagramError::InteractionSignatureMismatch {
                        vertex: vertex_id.0,
                        interaction,
                        actual,
                        expected,
                    });
                }
            }
        }
        let mut has_paired_external_connection = false;
        for (connection, states) in external_connections {
            match states.as_slice() {
                [_] => {}
                [left, right] if left != right => has_paired_external_connection = true,
                [state, _] => {
                    return Err(DiagramError::InvalidExternalConnectionStates {
                        connection,
                        state: state.as_str(),
                    });
                }
                _ => {
                    return Err(DiagramError::InvalidExternalConnectionSize {
                        connection,
                        legs: states.len(),
                    });
                }
            }
        }
        if has_paired_external_connection && self.cuts.is_empty() {
            return Err(DiagramError::MissingCrossSectionCuts);
        }
        if !has_paired_external_connection && !self.cuts.is_empty() {
            return Err(DiagramError::InvalidCut {
                cut: 0,
                message: "physical cuts require paired incoming/outgoing external connections"
                    .to_owned(),
            });
        }

        let universe = self
            .edges()
            .flat_map(|(edge, _, _)| {
                [
                    DiagramHalfEdge {
                        edge,
                        endpoint: DiagramEndpoint::Source,
                    },
                    DiagramHalfEdge {
                        edge,
                        endpoint: DiagramEndpoint::Target,
                    },
                ]
            })
            .collect::<BTreeSet<_>>();
        for (cut_index, cut) in self.cuts.iter().enumerate() {
            let left = cut.left.half_edges.iter().copied().collect::<BTreeSet<_>>();
            let right = cut
                .right
                .half_edges
                .iter()
                .copied()
                .collect::<BTreeSet<_>>();
            let oriented = cut.cut.iter().copied().collect::<BTreeSet<_>>();
            if left.len() != cut.left.half_edges.len()
                || right.len() != cut.right.half_edges.len()
                || oriented.len() != cut.cut.len()
            {
                return Err(DiagramError::InvalidCut {
                    cut: cut_index,
                    message: "half-edge sets contain duplicates".to_owned(),
                });
            }
            if !left.is_disjoint(&right)
                || left.union(&right).copied().collect::<BTreeSet<_>>() != universe
            {
                return Err(DiagramError::InvalidCut {
                    cut: cut_index,
                    message:
                        "left and right half-edge partitions are not disjoint and complementary"
                            .to_owned(),
                });
            }
            let mut expected_oriented = BTreeSet::new();
            for (edge, _, _) in self.edges() {
                let source = DiagramHalfEdge {
                    edge,
                    endpoint: DiagramEndpoint::Source,
                };
                let target = DiagramHalfEdge {
                    edge,
                    endpoint: DiagramEndpoint::Target,
                };
                match (left.contains(&source), left.contains(&target)) {
                    (true, false) => {
                        expected_oriented.insert(source);
                    }
                    (false, true) => {
                        expected_oriented.insert(target);
                    }
                    _ => {}
                }
            }
            if oriented != expected_oriented {
                return Err(DiagramError::InvalidCut {
                    cut: cut_index,
                    message: "oriented cut endpoints do not equal the crossing edges".to_owned(),
                });
            }

            for (vertex, incident) in vertex_half_edges.iter().enumerate() {
                if incident.first().is_some_and(|first| {
                    incident
                        .iter()
                        .any(|edge| left.contains(edge) != left.contains(first))
                }) {
                    return Err(DiagramError::InvalidCut {
                        cut: cut_index,
                        message: format!("vertex {vertex} is split between the two sides"),
                    });
                }
            }

            let summarize = |side: &BTreeSet<DiagramHalfEdge>| {
                let vertices = vertex_half_edges
                    .iter()
                    .enumerate()
                    .filter_map(|(vertex, incident)| {
                        incident
                            .first()
                            .is_some_and(|edge| side.contains(edge))
                            .then_some(VertexId(vertex))
                    })
                    .collect::<BTreeSet<_>>();
                let internal = vertices
                    .iter()
                    .copied()
                    .filter(|vertex| {
                        self.vertex(*vertex)
                            .is_some_and(|vertex| !vertex.is_external())
                    })
                    .collect::<BTreeSet<_>>();
                let mut coupling_orders = BTreeMap::new();
                for vertex in &internal {
                    if let Some(rule) = self.vertex(*vertex).and_then(|vertex| vertex.interaction) {
                        for (name, order) in self
                            .model
                            .vertex_rule_by_id(rule)?
                            .coupling_orders(&self.model)
                        {
                            *coupling_orders.entry(name).or_insert(0) += order;
                        }
                    }
                }
                let internal_edges = self
                    .edges()
                    .filter(|(_, endpoints, _)| {
                        internal.contains(&endpoints.source) && internal.contains(&endpoints.target)
                    })
                    .map(|(edge, endpoints, _)| (edge, endpoints))
                    .collect::<Vec<_>>();
                let mut remaining = internal.clone();
                let mut components = 0;
                while let Some(start) = remaining.pop_first() {
                    components += 1;
                    let mut queue = VecDeque::from([start]);
                    while let Some(vertex) = queue.pop_front() {
                        for (_, endpoints) in &internal_edges {
                            let neighbor = if endpoints.source == vertex {
                                Some(endpoints.target)
                            } else if endpoints.target == vertex {
                                Some(endpoints.source)
                            } else {
                                None
                            };
                            if let Some(neighbor) = neighbor
                                && remaining.remove(&neighbor)
                            {
                                queue.push_back(neighbor);
                            }
                        }
                    }
                }
                Ok::<_, DiagramError>((
                    coupling_orders,
                    internal_edges.len() + components - internal.len(),
                ))
            };
            let (left_orders, left_loops) = summarize(&left)?;
            let (right_orders, right_loops) = summarize(&right)?;
            if (left_orders, left_loops) != (cut.left.coupling_orders.clone(), cut.left.loop_count)
                || (right_orders, right_loops)
                    != (cut.right.coupling_orders.clone(), cut.right.loop_count)
            {
                return Err(DiagramError::InvalidCut {
                    cut: cut_index,
                    message: "stored side coupling orders or loop counts do not match the finalized topology"
                        .to_owned(),
                });
            }
        }
        if fragment_numerator.expand() != self.numerator.expand() {
            return Err(DiagramError::NumeratorFragmentMismatch);
        }
        self.loop_momentum_basis.validate(self)?;
        Ok(())
    }

    /// Serialize to a stable DOT dialect that can be parsed by [`Self::from_dot`].
    pub fn to_dot(&self) -> Result<String, DiagramError> {
        let mut output = String::new();
        let loop_momentum_edges = self
            .loop_momentum_basis
            .loop_edges
            .iter()
            .map(|edge| edge.0.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let cuts = serde_json::to_string(&self.cuts)?;
        writeln!(output, "digraph feynkit {{")?;
        writeln!(
            output,
            "  graph [feynkit_name={}, model_fingerprint={}, symmetry_factor={}, overall_factor={}, numerator={}, numerator_prefactor={}, projector={}, loop_momentum_edges={}, cuts={}];",
            Self::dot_string(&self.name)?,
            Self::dot_string(&self.model.fingerprint().to_string())?,
            self.symmetry_factor,
            Self::dot_string(&self.overall_factor.to_canonical_string())?,
            Self::dot_string(&self.numerator.to_canonical_string())?,
            Self::dot_string(&self.numerator_prefactor.to_canonical_string())?,
            Self::dot_string(&self.projector.to_canonical_string())?,
            Self::dot_string(&loop_momentum_edges)?,
            Self::dot_string(&cuts)?,
        )?;

        for (id, vertex) in self.vertices() {
            write!(
                output,
                "  v{} [id={}, feynkit_name={}",
                id.0,
                id.0,
                Self::dot_string(&vertex.name)?
            )?;
            if let Some(rule) = vertex.interaction {
                let rule_name = &self.model.vertex_rule_by_id(rule)?.name;
                write!(
                    output,
                    ", interaction_id={}, interaction={}",
                    rule.index(),
                    Self::dot_string(rule_name)?
                )?;
            }
            if let Some(external) = &vertex.external {
                write!(
                    output,
                    ", external_state={}, external_index={}, external_connection={}",
                    Self::dot_string(external.state.as_str())?,
                    external.index,
                    external.connection,
                )?;
            }
            write!(
                output,
                ", numerator={}",
                Self::dot_string(&vertex.numerator.to_canonical_string())?
            )?;
            writeln!(output, "]; ")?;
        }

        for (id, endpoints, edge) in self.edges() {
            let particle = self.model.particle_by_id(edge.particle)?;
            write!(
                output,
                "  v{} -> v{} [id={}, particle_id={}, pdg={}, particle={}, directed={}, source_slot={}, target_slot={}",
                endpoints.source.0,
                endpoints.target.0,
                id.0,
                edge.particle.index(),
                particle.pdg_code,
                Self::dot_string(&particle.name)?,
                edge.directed,
                edge.source_slot.0,
                edge.target_slot.0,
            )?;
            if !edge.directed {
                write!(output, ", dir=none")?;
            }
            write!(
                output,
                ", numerator={}",
                Self::dot_string(&edge.numerator.to_canonical_string())?
            )?;
            writeln!(output, "]; ")?;
        }
        writeln!(output, "}}")?;
        Ok(output)
    }

    /// Parse the stable FeynKit DOT dialect emitted by [`Self::to_dot`].
    pub fn from_dot(model: impl Into<Arc<Model>>, input: &str) -> Result<Self, DiagramError> {
        let model = model.into();
        let parsed: DotGraph = DotGraph::from_string(input)
            .map_err(|error| DiagramError::DotParse(error.to_string()))?;
        let serialized_fingerprint = parsed
            .global_data
            .statements
            .get("model_fingerprint")
            .ok_or_else(|| DiagramError::MissingDotAttribute {
                target: "graph".to_owned(),
                attribute: "model_fingerprint",
            })?;
        if serialized_fingerprint != &model.fingerprint().to_string() {
            return Err(DiagramError::InvalidDotAttribute {
                target: "graph".to_owned(),
                attribute: "model_fingerprint",
                value: serialized_fingerprint.clone(),
            });
        }
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
        let numerator = parsed
            .global_data
            .statements
            .get("numerator")
            .cloned()
            .unwrap_or_else(|| "1".to_owned());
        let numerator_prefactor = parsed
            .global_data
            .statements
            .get("numerator_prefactor")
            .cloned()
            .unwrap_or_else(|| "1".to_owned());
        let projector = parsed
            .global_data
            .statements
            .get("projector")
            .cloned()
            .unwrap_or_else(|| "1".to_owned());
        let loop_momentum_edges = parsed
            .global_data
            .statements
            .get("loop_momentum_edges")
            .map(|edges| {
                edges
                    .split(',')
                    .filter(|edge| !edge.is_empty())
                    .map(|edge| {
                        edge.parse()
                            .map(EdgeId)
                            .map_err(|_| DiagramError::InvalidDotAttribute {
                                target: "graph".to_owned(),
                                attribute: "loop_momentum_edges",
                                value: edges.clone(),
                            })
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .transpose()?;
        let cuts = parsed
            .global_data
            .statements
            .get("cuts")
            .ok_or_else(|| DiagramError::MissingDotAttribute {
                target: "graph".to_owned(),
                attribute: "cuts",
            })
            .and_then(|cuts| {
                let decoded = serde_json::from_str::<String>(&format!("\"{cuts}\""));
                decoded
                    .and_then(|cuts| serde_json::from_str::<Vec<DiagramCut>>(&cuts))
                    .map_err(|_| DiagramError::InvalidDotAttribute {
                        target: "graph".to_owned(),
                        attribute: "cuts",
                        value: cuts.clone(),
                    })
            })?;

        let mut builder = Self::builder(Arc::clone(&model), name)
            .symmetry_factor(symmetry_factor)
            .overall_factor(Self::parse_expression("overall factor", overall_factor)?)
            .numerator(Self::parse_expression("numerator", numerator)?)
            .numerator_prefactor(Self::parse_expression(
                "numerator prefactor",
                numerator_prefactor,
            )?)
            .projector(Self::parse_expression("projector", projector)?)
            .cuts(cuts);
        let mut node_map = BTreeMap::new();
        for (node, _, data) in parsed.iter_nodes() {
            let target = format!("vertex {}", node.0);
            let vertex_name = data
                .statements
                .get("feynkit_name")
                .cloned()
                .or_else(|| data.name.clone())
                .unwrap_or_else(|| format!("v{}", node.0));
            let interaction = data
                .statements
                .get("interaction_id")
                .map(|id| {
                    id.parse::<usize>()
                        .map_err(|_| DiagramError::InvalidDotAttribute {
                            target: target.clone(),
                            attribute: "interaction_id",
                            value: id.clone(),
                        })
                        .and_then(|id| model.vertex_rule_id_at(id).map_err(Into::into))
                })
                .transpose()?;
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
                    let connection = data
                        .statements
                        .get("external_connection")
                        .ok_or_else(|| DiagramError::MissingDotAttribute {
                            target: target.clone(),
                            attribute: "external_connection",
                        })?
                        .parse()
                        .map_err(|_| DiagramError::InvalidDotAttribute {
                            target: target.clone(),
                            attribute: "external_connection",
                            value: data.statements["external_connection"].clone(),
                        })?;
                    Ok(ExternalLeg {
                        index,
                        state: ExternalState::parse(state)?,
                        connection,
                    })
                })
                .transpose()?;
            let vertex = DiagramVertex {
                name: vertex_name,
                interaction,
                external,
                numerator: Self::parse_expression(
                    "vertex numerator",
                    data.statements
                        .get("numerator")
                        .cloned()
                        .unwrap_or_else(|| "1".to_owned()),
                )?,
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
            let _pdg: i64 = data
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
            let _particle_name =
                data.data
                    .statements
                    .get("particle")
                    .cloned()
                    .ok_or_else(|| DiagramError::MissingDotAttribute {
                        target: target.clone(),
                        attribute: "particle",
                    })?;
            let particle_id = data
                .data
                .statements
                .get("particle_id")
                .ok_or_else(|| DiagramError::MissingDotAttribute {
                    target: target.clone(),
                    attribute: "particle_id",
                })?
                .parse::<usize>()
                .map_err(|_| DiagramError::InvalidDotAttribute {
                    target: target.clone(),
                    attribute: "particle_id",
                    value: data.data.statements["particle_id"].clone(),
                })?;
            let mut edge = DiagramEdge {
                particle: model.particle_id_at(particle_id)?,
                directed: data.orientation != Orientation::Undirected,
                numerator: Self::parse_expression(
                    "edge numerator",
                    data.data
                        .statements
                        .get("numerator")
                        .cloned()
                        .unwrap_or_else(|| "1".to_owned()),
                )?,
                source_slot: VertexSlot(0),
                target_slot: VertexSlot(0),
            };
            let source_slot = data
                .data
                .statements
                .get("source_slot")
                .ok_or_else(|| DiagramError::MissingDotAttribute {
                    target: target.clone(),
                    attribute: "source_slot",
                })?
                .parse()
                .map(VertexSlot)
                .map_err(|_| DiagramError::InvalidDotAttribute {
                    target: target.clone(),
                    attribute: "source_slot",
                    value: data.data.statements["source_slot"].clone(),
                })?;
            let target_slot = data
                .data
                .statements
                .get("target_slot")
                .ok_or_else(|| DiagramError::MissingDotAttribute {
                    target: target.clone(),
                    attribute: "target_slot",
                })?
                .parse()
                .map(VertexSlot)
                .map_err(|_| DiagramError::InvalidDotAttribute {
                    target: target.clone(),
                    attribute: "target_slot",
                    value: data.data.statements["target_slot"].clone(),
                })?;
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
                    builder.add_edge_with_slots(source, sink, edge, source_slot, target_slot)?;
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
                        ExternalState::Incoming => builder.add_edge_with_slots(
                            external,
                            internal,
                            edge,
                            source_slot,
                            target_slot,
                        )?,
                        ExternalState::Outgoing => builder.add_edge_with_slots(
                            internal,
                            external,
                            edge,
                            source_slot,
                            target_slot,
                        )?,
                    };
                }
            }
        }
        let diagram = builder.build()?;
        let diagram = match loop_momentum_edges {
            Some(edges) => diagram.with_loop_momentum_edges(&edges),
            None => Ok(diagram),
        }?;
        diagram.validate()?;
        Ok(diagram)
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

    fn resolve_particle(
        &self,
        _edge_id: EdgeId,
        edge: &DiagramEdge,
    ) -> Result<(&Particle, ParticleId), DiagramError> {
        let particle = self.model.particle_by_id(edge.particle)?;
        let base = if particle.is_antiparticle() {
            let antiparticle = self.model.antiparticle(particle)?;
            self.model.particle_id(&antiparticle.name)?
        } else {
            edge.particle
        };
        Ok((particle, base))
    }

    fn serde_view(&self) -> DiagramSerde {
        DiagramSerde {
            id: self.id,
            model: self.model.fingerprint(),
            name: self.name.clone(),
            vertices: self
                .vertices()
                .map(|(_, vertex)| DiagramVertexSerde {
                    name: vertex.name.clone(),
                    interaction: vertex.interaction,
                    external: vertex.external.clone(),
                    numerator: vertex.numerator.to_canonical_string(),
                })
                .collect(),
            edges: self
                .edges()
                .map(|(_, endpoints, edge)| {
                    (
                        endpoints,
                        DiagramEdgeSerde {
                            particle: edge.particle,
                            directed: edge.directed,
                            numerator: edge.numerator.to_canonical_string(),
                            source_slot: edge.source_slot,
                            target_slot: edge.target_slot,
                        },
                    )
                })
                .collect(),
            symmetry_factor: self.symmetry_factor,
            overall_factor: self.overall_factor.to_canonical_string(),
            numerator: self.numerator.to_canonical_string(),
            numerator_prefactor: self.numerator_prefactor.to_canonical_string(),
            projector: self.projector.to_canonical_string(),
            loop_momentum_basis: self.loop_momentum_basis.clone(),
            cuts: self.cuts.clone(),
        }
    }

    fn from_serde(model: Arc<Model>, data: DiagramSerde) -> Result<Self, DiagramError> {
        let actual_fingerprint = model.fingerprint();
        if data.model != actual_fingerprint {
            return Err(DiagramError::ModelFingerprintMismatch {
                serialized: data.model,
                actual: actual_fingerprint,
            });
        }
        let serialized_id = data.id;
        let mut builder = Self::builder(model, data.name)
            .symmetry_factor(data.symmetry_factor)
            .overall_factor(Self::parse_expression(
                "overall factor",
                data.overall_factor,
            )?)
            .numerator(Self::parse_expression("numerator", data.numerator)?)
            .numerator_prefactor(Self::parse_expression(
                "numerator prefactor",
                data.numerator_prefactor,
            )?)
            .projector(Self::parse_expression("projector", data.projector)?)
            .loop_momentum_basis(data.loop_momentum_basis)
            .cuts(data.cuts);
        for vertex in data.vertices {
            builder.add_vertex(DiagramVertex {
                name: vertex.name,
                interaction: vertex.interaction,
                external: vertex.external,
                numerator: Self::parse_expression("vertex numerator", vertex.numerator)?,
            });
        }
        for (endpoints, edge) in data.edges {
            builder.add_edge_with_slots(
                endpoints.source,
                endpoints.target,
                DiagramEdge {
                    particle: edge.particle,
                    directed: edge.directed,
                    numerator: Self::parse_expression("edge numerator", edge.numerator)?,
                    source_slot: edge.source_slot,
                    target_slot: edge.target_slot,
                },
                edge.source_slot,
                edge.target_slot,
            )?;
        }
        let diagram = builder.build()?;
        if diagram.id != serialized_id {
            return Err(DiagramError::DiagramIdMismatch {
                serialized: serialized_id,
                actual: diagram.id,
            });
        }
        diagram.validate()?;
        Ok(diagram)
    }

    fn dot_string(value: &str) -> Result<String, DiagramError> {
        serde_json::to_string(value).map_err(DiagramError::from)
    }

    fn parse_expression(field: &'static str, expression: String) -> Result<Atom, DiagramError> {
        Atom::parse(&expression, "feynkit_graph", ParseSettings::default()).map_err(|message| {
            DiagramError::SymbolicParse {
                field,
                expression,
                message,
            }
        })
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
    model: Arc<Model>,
    name: String,
    vertices: Vec<DiagramVertex>,
    edges: Vec<(EdgeEndpoints, DiagramEdge)>,
    symmetry_factor: u64,
    overall_factor: Atom,
    numerator: Atom,
    numerator_prefactor: Atom,
    projector: Atom,
    loop_momentum_basis: Option<LoopMomentumBasis>,
    cuts: Vec<DiagramCut>,
}

impl FeynmanDiagramBuilder {
    pub fn new(model: impl Into<Arc<Model>>, name: impl Into<String>) -> Self {
        Self {
            model: model.into(),
            name: name.into(),
            vertices: Vec::new(),
            edges: Vec::new(),
            symmetry_factor: 1,
            overall_factor: Atom::one(),
            numerator: Atom::one(),
            numerator_prefactor: Atom::one(),
            projector: Atom::one(),
            loop_momentum_basis: None,
            cuts: Vec::new(),
        }
    }

    /// Record the raw automorphism order as diagnostic provenance.
    ///
    /// The builder does not fold it into the authoritative `overall_factor`.
    pub fn symmetry_factor(mut self, symmetry_factor: u64) -> Self {
        self.symmetry_factor = symmetry_factor;
        self
    }

    pub fn overall_factor(mut self, overall_factor: Atom) -> Self {
        self.overall_factor = overall_factor;
        self
    }

    pub fn numerator(mut self, numerator: Atom) -> Self {
        self.numerator = numerator;
        self
    }

    pub fn numerator_prefactor(mut self, numerator_prefactor: Atom) -> Self {
        self.numerator_prefactor = numerator_prefactor;
        self
    }

    pub fn projector(mut self, projector: Atom) -> Self {
        self.projector = projector;
        self
    }

    pub fn loop_momentum_basis(mut self, basis: LoopMomentumBasis) -> Self {
        self.loop_momentum_basis = Some(basis);
        self
    }

    pub fn cuts(mut self, cuts: Vec<DiagramCut>) -> Self {
        self.cuts = cuts;
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
        mut edge: DiagramEdge,
    ) -> Result<EdgeId, DiagramError> {
        let degree = |vertex: VertexId| {
            self.edges
                .iter()
                .map(|(endpoints, _)| {
                    usize::from(endpoints.source == vertex)
                        + usize::from(endpoints.target == vertex)
                })
                .sum()
        };
        let source_slot = VertexSlot(degree(source));
        let target_slot = VertexSlot(degree(target) + usize::from(source == target));
        edge.source_slot = source_slot;
        edge.target_slot = target_slot;
        self.add_edge_with_slots(source, target, edge, source_slot, target_slot)
    }

    pub fn add_edge_with_slots(
        &mut self,
        source: VertexId,
        target: VertexId,
        mut edge: DiagramEdge,
        source_slot: VertexSlot,
        target_slot: VertexSlot,
    ) -> Result<EdgeId, DiagramError> {
        for vertex in [source, target] {
            if vertex.0 >= self.vertices.len() {
                return Err(DiagramError::UnknownVertex {
                    vertex: vertex.0,
                    vertices: self.vertices.len(),
                });
            }
        }
        edge.source_slot = source_slot;
        edge.target_slot = target_slot;
        let id = EdgeId(self.edges.len());
        self.edges.push((EdgeEndpoints { source, target }, edge));
        Ok(id)
    }

    pub fn build(mut self) -> Result<FeynmanDiagram, DiagramError> {
        for cut in &mut self.cuts {
            cut.cut.sort();
            cut.cut.dedup();
            cut.left.half_edges.sort();
            cut.left.half_edges.dedup();
            cut.right.half_edges.sort();
            cut.right.half_edges.dedup();
        }
        self.cuts.sort();
        self.cuts.dedup();
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
        let mut diagram = FeynmanDiagram {
            model: self.model,
            id: DiagramId(0),
            name: self.name,
            graph: graph_builder.build(),
            symmetry_factor: self.symmetry_factor,
            overall_factor: self.overall_factor,
            numerator: self.numerator,
            numerator_prefactor: self.numerator_prefactor,
            projector: self.projector,
            loop_momentum_basis: self.loop_momentum_basis.unwrap_or(LoopMomentumBasis {
                tree_edges: Vec::new(),
                loop_edges: Vec::new(),
                external_edges: Vec::new(),
                dependent_externals: Vec::new(),
                edge_signatures: BTreeMap::new(),
            }),
            cuts: self.cuts,
        };
        if diagram.loop_momentum_basis.edge_signatures.is_empty() && diagram.graph.n_edges() != 0 {
            diagram.loop_momentum_basis = diagram
                .loop_momentum_bases_with_limit(1)?
                .into_iter()
                .next()
                .ok_or(DiagramError::MissingLoopMomentumBasis)?;
        }
        diagram.id = DiagramId::from_key(diagram.model.fingerprint(), &diagram.structural_key()?)?;
        Ok(diagram)
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
    use std::sync::Arc;

    fn scalar_model() -> Arc<Model> {
        Arc::new(Model::from_json(
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
        .unwrap())
    }

    fn one_loop() -> FeynmanDiagram {
        let model = scalar_model();
        let rule = model.vertex_rule_id("V_1").unwrap();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "bubble");
        let incoming =
            builder.add_vertex(DiagramVertex::external("p1", 0, ExternalState::Incoming));
        let outgoing =
            builder.add_vertex(DiagramVertex::external("p2", 1, ExternalState::Outgoing));
        let left = builder.add_vertex(DiagramVertex::interaction("v0", rule));
        let right = builder.add_vertex(DiagramVertex::interaction("v1", rule));
        let scalar = || DiagramEdge::new(particle, false);
        builder.add_edge(incoming, left, scalar()).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        builder.add_edge(right, outgoing, scalar()).unwrap();
        builder.build().unwrap()
    }

    fn one_loop_reordered() -> FeynmanDiagram {
        let model = scalar_model();
        let rule = model.vertex_rule_id("V_1").unwrap();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "renamed");
        let right = builder.add_vertex(DiagramVertex::interaction("right", rule));
        let incoming =
            builder.add_vertex(DiagramVertex::external("in", 0, ExternalState::Incoming));
        let left = builder.add_vertex(DiagramVertex::interaction("left", rule));
        let outgoing =
            builder.add_vertex(DiagramVertex::external("out", 1, ExternalState::Outgoing));
        let scalar = || DiagramEdge::new(particle, false);
        builder
            .add_edge_with_slots(right, outgoing, scalar(), VertexSlot(2), VertexSlot(0))
            .unwrap();
        builder
            .add_edge_with_slots(left, right, scalar(), VertexSlot(1), VertexSlot(0))
            .unwrap();
        builder
            .add_edge_with_slots(incoming, left, scalar(), VertexSlot(0), VertexSlot(0))
            .unwrap();
        builder
            .add_edge_with_slots(left, right, scalar(), VertexSlot(2), VertexSlot(1))
            .unwrap();
        builder.build().unwrap()
    }

    fn cut_scalar_line() -> FeynmanDiagram {
        let model = scalar_model();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "cut-scalar-line");
        let incoming = builder.add_vertex(DiagramVertex::external_in_connection(
            "p-in",
            0,
            ExternalState::Incoming,
            0,
        ));
        let interaction = builder.add_vertex(DiagramVertex {
            name: "interaction".to_owned(),
            interaction: None,
            external: None,
            numerator: Atom::one(),
        });
        let outgoing = builder.add_vertex(DiagramVertex::external_in_connection(
            "p-out",
            1,
            ExternalState::Outgoing,
            0,
        ));
        let incoming_edge = builder
            .add_edge(incoming, interaction, DiagramEdge::new(particle, false))
            .unwrap();
        let outgoing_edge = builder
            .add_edge(interaction, outgoing, DiagramEdge::new(particle, false))
            .unwrap();
        let incoming_source = DiagramHalfEdge {
            edge: incoming_edge,
            endpoint: DiagramEndpoint::Source,
        };
        let incoming_target = DiagramHalfEdge {
            edge: incoming_edge,
            endpoint: DiagramEndpoint::Target,
        };
        let outgoing_source = DiagramHalfEdge {
            edge: outgoing_edge,
            endpoint: DiagramEndpoint::Source,
        };
        let outgoing_target = DiagramHalfEdge {
            edge: outgoing_edge,
            endpoint: DiagramEndpoint::Target,
        };
        let empty_side = |half_edges| DiagramCutSide {
            half_edges,
            coupling_orders: BTreeMap::new(),
            loop_count: 0,
        };
        builder
            .cuts(vec![DiagramCut {
                cut: vec![outgoing_source],
                left: empty_side(vec![incoming_source, incoming_target, outgoing_source]),
                right: empty_side(vec![outgoing_target]),
            }])
            .build()
            .unwrap()
    }

    fn fermion_model() -> Arc<Model> {
        Arc::new(Model::from_json(
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
        .unwrap())
    }

    fn directed_fermion_line() -> FeynmanDiagram {
        let model = fermion_model();
        let anti_rule = model.vertex_rule_id("V_af").unwrap();
        let particle_rule = model.vertex_rule_id("V_f").unwrap();
        let fermion = model.particle_id("f").unwrap();
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "fermion-line");
        let anti = builder.add_vertex(DiagramVertex::interaction("anti", anti_rule));
        let particle = builder.add_vertex(DiagramVertex::interaction("particle", particle_rule));
        builder
            .add_edge(anti, particle, DiagramEdge::new(fermion, true))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn json_and_dot_round_trip() {
        let diagram = one_loop();
        let from_json =
            FeynmanDiagram::from_json(diagram.model_arc(), &diagram.to_json().unwrap()).unwrap();
        assert_eq!(from_json.loop_count(), 1);
        assert_eq!(from_json.vertices().count(), 4);

        let from_dot =
            FeynmanDiagram::from_dot(diagram.model_arc(), &diagram.to_dot().unwrap()).unwrap();
        assert_eq!(from_dot.name(), "bubble");
        assert_eq!(from_dot.loop_count(), 1);
        assert_eq!(from_dot.edges().count(), 4);
    }

    #[test]
    fn validates_and_round_trips_typed_cuts() {
        let diagram = cut_scalar_line();
        diagram.validate().unwrap();

        let from_json =
            FeynmanDiagram::from_json(diagram.model_arc(), &diagram.to_json().unwrap()).unwrap();
        assert_eq!(from_json.cuts(), diagram.cuts());

        let from_dot =
            FeynmanDiagram::from_dot(diagram.model_arc(), &diagram.to_dot().unwrap()).unwrap();
        assert_eq!(from_dot.cuts(), diagram.cuts());

        let mut invalid_cut = diagram.cuts()[0].clone();
        invalid_cut.cut.clear();
        assert!(matches!(
            diagram.with_cuts(vec![invalid_cut]),
            Err(DiagramError::InvalidCut { cut: 0, message })
                if message.contains("crossing edges")
        ));
    }

    #[test]
    fn json_import_rejects_invalid_external_structure() {
        let diagram = one_loop();
        let mut missing_edge = diagram.serde_view();
        missing_edge.edges.remove(0);
        assert!(matches!(
            FeynmanDiagram::from_json(
                diagram.model_arc(),
                &serde_json::to_string(&missing_edge).unwrap()
            ),
            Err(error) if error.to_string().contains("external vertex 0 has degree 0")
        ));

        let mut external_edge = diagram.serde_view();
        external_edge.edges[0].0.target = VertexId(1);
        assert!(matches!(
            FeynmanDiagram::from_json(
                diagram.model_arc(),
                &serde_json::to_string(&external_edge).unwrap()
            ),
            Err(error) if error.to_string().contains("edge 0 connects two external vertices")
        ));
    }

    #[test]
    fn json_import_rejects_a_different_model_fingerprint() {
        let diagram = one_loop();
        let wrong_model = Arc::new(Model::empty("different-model"));

        assert!(matches!(
            FeynmanDiagram::from_json(wrong_model, &diagram.to_json().unwrap()),
            Err(DiagramError::ModelFingerprintMismatch { .. })
        ));
    }

    #[test]
    fn dot_import_rejects_external_index_overflow() {
        let diagram = one_loop();
        let dot = diagram.to_dot().unwrap().replacen(
            "external_index=0",
            &format!("external_index={}", usize::MAX),
            1,
        );
        assert!(matches!(
            FeynmanDiagram::from_dot(diagram.model_arc(), &dot),
            Err(DiagramError::DotExternalIndexOverflow(index)) if index == usize::MAX
        ));
    }

    #[test]
    fn dot_import_requires_canonical_particle_attribute() {
        let diagram = one_loop();
        let dot = diagram
            .to_dot()
            .unwrap()
            .replacen("particle=\"phi\", ", "", 1);
        assert!(matches!(
            FeynmanDiagram::from_dot(diagram.model_arc(), &dot),
            Err(DiagramError::MissingDotAttribute {
                target,
                attribute: "particle"
            }) if target == "edge 0"
        ));
    }

    #[test]
    fn validates_particle_identity_and_interaction_signature() {
        let diagram = one_loop();
        diagram.validate().unwrap();

        let unknown_particle = diagram
            .map_data(
                |_, vertex| vertex.clone(),
                |id, _, edge| {
                    let mut edge = edge.clone();
                    if id == EdgeId(0) {
                        edge.particle = ParticleId::from_index(usize::MAX);
                    }
                    edge
                },
            )
            .unwrap();
        assert!(matches!(
            unknown_particle.validate(),
            Err(DiagramError::Model(ModelError::NotFound {
                kind: feynkit_model::EntityKind::Particle,
                ..
            }))
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
            wrong_flow.validate(),
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
        let model = scalar_model();
        let rule = model.vertex_rule_id("V_1").unwrap();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "triangle_with_tail");
        let vertices: Vec<_> = (0..4)
            .map(|index| builder.add_vertex(DiagramVertex::interaction(format!("v{index}"), rule)))
            .collect();
        let scalar = || DiagramEdge::new(particle, false);
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
                    edge.numerator = Atom::parse(
                        format!("N({})", id.0),
                        "feynkit_graph_test",
                        ParseSettings::default(),
                    )
                    .unwrap();
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
        assert!(mapped.edges().all(|(_, _, edge)| !edge.numerator.is_one()));
    }

    #[test]
    fn canonical_key_is_insertion_order_invariant_and_color_sensitive() {
        let diagram = one_loop();
        let key = diagram.canonical_key().unwrap();
        assert_eq!(key, one_loop_reordered().canonical_key().unwrap());

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
        assert_ne!(key, changed_state.canonical_key().unwrap());

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
        assert_ne!(key, changed_interaction.canonical_key().unwrap());
    }

    #[test]
    fn reverses_particle_flow_with_antiparticle_semantics() {
        let diagram = directed_fermion_line();
        diagram.validate().unwrap();
        let original_json = diagram.to_json().unwrap();
        let original_key = diagram.canonical_key().unwrap();

        let reversed = diagram.reverse_edge(EdgeId(0)).unwrap();
        let (_, endpoints, edge) = reversed.edges().next().unwrap();
        assert_eq!(endpoints.source, VertexId(1));
        assert_eq!(endpoints.target, VertexId(0));
        assert_eq!(edge.particle, reversed.model().particle_id("f~").unwrap());
        assert_eq!(reversed.canonical_key().unwrap(), original_key);

        let restored = reversed.reverse_edge(EdgeId(0)).unwrap();
        assert_eq!(restored.to_json().unwrap(), original_json);
        assert!(matches!(
            diagram.reverse_edge(EdgeId(1)),
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
        let mut builder = FeynmanDiagram::builder(scalar_model(), "invalid");
        builder.add_vertex(DiagramVertex::external("p1", 0, ExternalState::Incoming));
        builder.add_vertex(DiagramVertex::external("p2", 0, ExternalState::Outgoing));
        assert!(matches!(
            builder.build(),
            Err(DiagramError::DuplicateExternalIndex(0))
        ));
    }
}
