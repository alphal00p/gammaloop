//! Lightweight, application-independent Feynman diagrams.
//!
//! This crate owns the structural and symbolic diagram layer, including the
//! physical cut partitions selected during canonical cross-section generation.
//! Runtime cut evaluators, subtraction data, parametrizations, and caches belong
//! to applications consuming these diagrams.

#![forbid(unsafe_code)]

mod compact_dot;
mod display;
pub mod expressions;
mod finalization;
mod integrals;
mod power_counting;
pub mod routing;
mod subgraph;
pub mod symbols;
pub mod thresholds;
mod uv;

pub use integrals::{IntegralFamily, IntegralFamilyError, IntegralMapping, PropagatorMapping};
pub use power_counting::DOD;

// Symbolica does not permit adding tags after a bare symbol with the same name
// has been registered. Claim momentum, index and Lorentz heads during global state
// initialization so parsing and generation always agree on their tensor types.
// Use Spenso's canonical tag names and shared printer so FeynKit interoperates
// with the Spenso instance embedded by the host.
symbolica::initialize!(|| spenso::symbolica_init::in_symbolica_initializer(|| {
    let _ = Minkowski {}.to_symbolic([Atom::Zero]);
    symbols::momentum();
    symbols::loop_momentum();
    symbols::external_momentum();
    symbols::edge_index();
    symbols::vertex_index();
    symbols::hedge_index();
    symbols::dimension();
    symbols::denominator();
    symbols::u();
    symbols::ubar();
    symbols::v();
    symbols::vbar();
    symbols::epsilon();
    symbols::epsilonbar();
}));

pub fn momentum_symbol() -> symbolica::atom::Symbol {
    symbols::momentum()
}

use std::{
    collections::{BTreeMap, BTreeSet},
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
        involution::{EdgeIndex, Flow, HedgePair, Orientation},
        subgraph::{Inclusion, SuBitGraph, SubSetLike, SubSetOps},
        tree::SimpleTraversalTree,
    },
    parser::{DotGraph, set::DotGraphSet},
};
use serde::{Deserialize, Serialize};
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::{
    atom::{Atom, AtomCore, UserData},
    graph::Graph as CanonicalGraph,
    parser::{ParseSettings, Token},
    state::{State, Workspace},
    with_default_namespace,
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

/// Whether an external momentum carrier belongs to the initial or final state.
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
}

/// Metadata carried by an external edge, including sewn initial-state connections.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct ExternalLeg {
    pub name: String,
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
    pub numerator: Atom,
}

impl DiagramVertex {
    pub fn interaction(name: impl Into<String>, rule: VertexRuleId) -> Self {
        Self {
            name: name.into(),
            interaction: Some(rule),
            numerator: Atom::one(),
        }
    }
}

/// Structural and symbolic information associated with an edge.
#[derive(Debug, Clone)]
pub struct DiagramEdge {
    pub external: Option<ExternalLeg>,
    pub is_dummy: bool,
    pub particle: ParticleId,
    pub directed: bool,
    pub numerator: Atom,
    source_slot: VertexSlot,
    target_slot: VertexSlot,
}

impl DiagramEdge {
    pub fn new(particle: ParticleId, directed: bool) -> Self {
        Self {
            external: None,
            is_dummy: false,
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
    pub source: Option<VertexId>,
    pub target: Option<VertexId>,
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

/// One topology-only threshold candidate retained during cross-section finalization.
///
/// Unlike [`DiagramCut`], this partition is not required to match the requested
/// physical final state. It records the complete, process-independent s-channel
/// cut inventory used to construct threshold counterterms downstream. Keeping
/// this metadata separate prevents topology candidates from being mistaken for
/// physical Cutkosky cuts.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct DiagramThresholdCandidate {
    /// The half-edge on the left side of every oriented crossing edge.
    pub cut: Vec<DiagramHalfEdge>,
    /// Every half-edge belonging to the source side of the partition.
    pub left: Vec<DiagramHalfEdge>,
    /// Every half-edge belonging to the target side of the partition.
    pub right: Vec<DiagramHalfEdge>,
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
    TopologyThresholdSide {
        side: CanonicalCutSide,
    },
}

/// Canonically labeled edge retained by [`CanonicalDiagramKey`].
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CanonicalDiagramEdge {
    Diagram {
        endpoints: EdgeEndpoints,
        particle: ParticleId,
        external: Option<ExternalLeg>,
        is_dummy: bool,
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
        external: Option<ExternalLeg>,
        is_dummy: bool,
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
        let expected = diagram
            .basis_from_tree(&self.tree_edges)?
            .with_loop_edge_order(&self.loop_edges)?;
        if *self != expected {
            return Err(DiagramError::InvalidLoopMomentumBasis(
                "stored signatures or edge roles do not match the shared graph routing".into(),
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
    IntegralFamily(#[from] IntegralFamilyError),
    #[error("cannot expand UV counterterm: {0}")]
    UvExpansion(String),
    #[error("cannot determine superficial UV degree: {0}")]
    UvPowerCounting(String),
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
    #[error("diagram DOT model fingerprint {serialized} does not match supplied model {actual}")]
    DotModelFingerprintMismatch {
        serialized: String,
        actual: ModelFingerprint,
    },
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
    #[error("invalid topology threshold candidate {candidate}: {message}")]
    InvalidThresholdCandidate { candidate: usize, message: String },
    #[error(
        "diagram-wide numerator differs from the product of the vertex and edge numerator fragments"
    )]
    NumeratorFragmentMismatch,
    #[error("cannot replace the numerator of a diagram without an interaction vertex")]
    MissingNumeratorAnchor,
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
    half_edge_order: Vec<DiagramHalfEdge>,
    symmetry_factor: u64,
    overall_factor: String,
    numerator: String,
    numerator_prefactor: String,
    projector: String,
    loop_momentum_basis: LoopMomentumBasis,
    cuts: Vec<DiagramCut>,
    topology_threshold_candidates: Vec<DiagramThresholdCandidate>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct DiagramVertexSerde {
    name: String,
    interaction: Option<VertexRuleId>,
    numerator: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DiagramEdgeSerde {
    external: Option<ExternalLeg>,
    is_dummy: bool,
    orientation: Orientation,
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
    topology_threshold_candidates: Vec<DiagramThresholdCandidate>,
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

    /// Product of internal quadratic propagators in GammaLoop's tagged form.
    ///
    /// Each `denom(edge, Q(edge), mass², Q(edge)² - mass²)` uses the symbolic
    /// dimension `gammalooprs::dim` and the same momenta as the numerator. Masses
    /// remain symbolic except the UFO `ZERO` parameter. External carriers,
    /// widths, and an imaginary prescription are excluded. Custom UFO
    /// propagator denominators are not instantiated here.
    pub fn denominator_expression(&self) -> Result<Atom, DiagramError> {
        self.denominator_of(&self.internal_subgraph(), &BTreeMap::new())
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

    /// Physical final-state particles in the stored order of a finalized cut.
    ///
    /// The native left side contains the incoming endpoint of each sewn initial
    /// state. A source endpoint on that side carries the stored species into
    /// the final state; a target endpoint carries its antiparticle. The same
    /// orientation determines the positive-energy cut momentum routing.
    pub fn cut_particles(&self, cut: &DiagramCut) -> Result<Vec<ParticleId>, DiagramError> {
        cut.cut
            .iter()
            .map(|half| {
                if self.half_edge_id(*half).is_none() {
                    return Err(DiagramError::Invariant {
                        operation: "reading cut particles",
                        message: format!("cut references missing half-edge {half:?}"),
                    });
                }
                let particle = self.graph[EdgeIndex(half.edge.0)].particle;
                Ok(match half.endpoint {
                    DiagramEndpoint::Source => particle,
                    DiagramEndpoint::Target => self.model.particle_by_id(particle)?.antiparticle,
                })
            })
            .collect()
    }

    /// Complete topology-only s-channel inventory used for threshold construction.
    pub fn topology_threshold_candidates(&self) -> &[DiagramThresholdCandidate] {
        &self.topology_threshold_candidates
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
        self.graph.iter_edges().map(|(pair, id, data)| {
            let endpoints = match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    EdgeEndpoints {
                        source: Some(VertexId(self.graph.node_id(source).0)),
                        target: Some(VertexId(self.graph.node_id(sink).0)),
                    }
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let vertex = Some(VertexId(self.graph.node_id(hedge).0));
                    match flow {
                        Flow::Source => EdgeEndpoints {
                            source: vertex,
                            target: None,
                        },
                        Flow::Sink => EdgeEndpoints {
                            source: None,
                            target: vertex,
                        },
                    }
                }
            };
            (EdgeId(id.0), endpoints, data.data)
        })
    }

    pub fn vertex(&self, id: impl Into<Option<VertexId>>) -> Option<&DiagramVertex> {
        id.into()
            .and_then(|id| (id.0 < self.graph.n_nodes()).then(|| &self.graph[NodeIndex(id.0)]))
    }

    /// Transform vertex and edge metadata without changing the diagram
    /// topology or stable identifiers.
    pub fn map_data<V, E>(&self, mut map_vertex: V, mut map_edge: E) -> Result<Self, DiagramError>
    where
        V: FnMut(VertexId, &DiagramVertex) -> DiagramVertex,
        E: FnMut(EdgeId, EdgeEndpoints, &DiagramEdge) -> DiagramEdge,
    {
        let mut diagram = self.clone();
        let endpoints = self
            .edges()
            .map(|(id, endpoints, _)| (id, endpoints))
            .collect::<BTreeMap<_, _>>();
        diagram.graph = self.graph.clone().map(
            |_, id, vertex| map_vertex(VertexId(id.0), &vertex),
            |_, _, _, id, edge| {
                edge.map(|edge| map_edge(EdgeId(id.0), endpoints[&EdgeId(id.0)], &edge))
            },
            |_, data| data,
        );
        diagram.id = DiagramId::from_key(diagram.model.fingerprint(), &diagram.structural_key()?)?;
        Ok(diagram)
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
                external: None,
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
                        endpoints.get(&half_edge.edge).and_then(|endpoints| {
                            match half_edge.endpoint {
                                DiagramEndpoint::Source => endpoints.source,
                                DiagramEndpoint::Target => endpoints.target,
                            }
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
        for candidate in &self.topology_threshold_candidates {
            let left = graph.add_node(CanonicalDiagramVertex::TopologyThresholdSide {
                side: CanonicalCutSide::Left,
            });
            let right = graph.add_node(CanonicalDiagramVertex::TopologyThresholdSide {
                side: CanonicalCutSide::Right,
            });
            graph
                .add_edge(left, right, false, CanonicalEdgeColor::CutPair)
                .map_err(|error| DiagramError::Invariant {
                    operation: "deriving the diagram ID",
                    message: error.to_string(),
                })?;
            for (side, half_edges) in [(left, &candidate.left), (right, &candidate.right)] {
                let vertices = half_edges
                    .iter()
                    .filter_map(|half_edge| {
                        endpoints.get(&half_edge.edge).and_then(|endpoints| {
                            match half_edge.endpoint {
                                DiagramEndpoint::Source => endpoints.source,
                                DiagramEndpoint::Target => endpoints.target,
                            }
                        })
                    })
                    .collect::<BTreeSet<_>>();
                for vertex in vertices {
                    graph
                        .add_edge(side, vertex.0, false, CanonicalEdgeColor::CutMembership)
                        .map_err(|error| DiagramError::Invariant {
                            operation: "encoding canonical topology-threshold membership",
                            message: error.to_string(),
                        })?;
                }
            }
        }
        for (edge_id, endpoints, edge) in self.edges() {
            let mut particle = edge.particle;
            let mut external = edge.external.clone();
            if let Some(leg) = &mut external {
                leg.name.clear();
            }
            let mut terminal = || {
                graph.add_node(CanonicalDiagramVertex::Diagram {
                    interaction: None,
                    external: external.clone(),
                })
            };
            let (mut source, mut target) = (
                endpoints.source.map_or_else(&mut terminal, |v| v.0),
                endpoints.target.map_or_else(&mut terminal, |v| v.0),
            );
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
                        external,
                        is_dummy: edge.is_dummy,
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
                        source: Some(VertexId(edge.vertices.0)),
                        target: Some(VertexId(edge.vertices.1)),
                    };
                    match &edge.data {
                        CanonicalEdgeColor::Diagram {
                            particle,
                            external,
                            is_dummy,
                            source_slot,
                            target_slot,
                        } => CanonicalDiagramEdge::Diagram {
                            endpoints,
                            particle: *particle,
                            external: external.clone(),
                            is_dummy: *is_dummy,
                            directed: edge.directed,
                            source_slot: *source_slot,
                            target_slot: *target_slot,
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
        let mut topology_threshold_candidates = self.topology_threshold_candidates.clone();
        for candidate in &mut topology_threshold_candidates {
            for half_edge in &mut candidate.cut {
                reverse_half_edge(half_edge);
            }
            for half_edge in &mut candidate.left {
                reverse_half_edge(half_edge);
            }
            for half_edge in &mut candidate.right {
                reverse_half_edge(half_edge);
            }
        }
        let mut builder = Self::builder(Arc::clone(&self.model), &self.name)
            .symmetry_factor(self.symmetry_factor)
            .overall_factor(self.overall_factor.clone())
            .numerator(self.numerator.clone())
            .numerator_prefactor(self.numerator_prefactor.clone())
            .projector(self.projector.clone())
            .cuts(cuts)
            .topology_threshold_candidates(topology_threshold_candidates);
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
            .edges()
            .filter_map(|(_, _, edge)| edge.external.as_ref().map(|leg| leg.index))
            .collect();
        if let Some(index) = relabeling
            .keys()
            .find(|index| !external_indices.contains(index))
        {
            return Err(DiagramError::UnknownExternalIndex(*index));
        }

        let diagram = self.map_data(
            |_, vertex| vertex.clone(),
            |_, _, edge| {
                let mut edge = edge.clone();
                if let Some(external) = &mut edge.external
                    && let Some(index) = relabeling.get(&external.index)
                {
                    external.index = *index;
                }
                edge
            },
        )?;
        diagram.validate()?;
        Ok(diagram)
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

    /// Return this diagram with a replacement scalar or tensor numerator.
    ///
    /// The replacement becomes the numerator fragment of the first
    /// interaction vertex; every other vertex and edge fragment becomes one.
    /// This keeps the aggregate numerator authoritative while preserving the
    /// interaction and particle assignments needed to interpret the topology.
    /// A diagram without an interaction vertex cannot own such a replacement.
    /// Diagram IDs identify the model-resolved topology and remain unchanged.
    pub fn with_numerator(self, numerator: Atom) -> Result<Self, DiagramError> {
        let anchor = self
            .vertices()
            .map(|(id, _)| id)
            .next()
            .ok_or(DiagramError::MissingNumeratorAnchor)?;
        let anchor_numerator = numerator.clone();
        let mut replaced = self.map_data(
            |id, vertex| {
                let mut vertex = vertex.clone();
                vertex.numerator = if id == anchor {
                    anchor_numerator.clone()
                } else {
                    Atom::one()
                };
                vertex
            },
            |_, _, edge| {
                let mut edge = edge.clone();
                edge.numerator = Atom::one();
                edge
            },
        )?;
        replaced.numerator = numerator;
        replaced.validate()?;
        Ok(replaced)
    }

    /// Return this diagram with a replacement external-state projector.
    ///
    /// The projector is numerator payload and therefore does not alter the
    /// topology-derived diagram ID.
    pub fn with_projector(mut self, projector: Atom) -> Self {
        self.projector = projector;
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
        cuts.sort_by_cached_key(|cut| self.cut_order_key(&cut.cut));
        cuts.dedup();
        self.cuts = cuts;
        self.validate()?;
        self.id = DiagramId::from_key(self.model.fingerprint(), &self.structural_key()?)?;
        Ok(self)
    }

    /// Replace the complete topology-only threshold inventory.
    pub fn with_topology_threshold_candidates(
        mut self,
        mut candidates: Vec<DiagramThresholdCandidate>,
    ) -> Result<Self, DiagramError> {
        for candidate in &mut candidates {
            candidate.cut.sort();
            candidate.cut.dedup();
            candidate.left.sort();
            candidate.left.dedup();
            candidate.right.sort();
            candidate.right.dedup();
        }
        let ignored = self.initial_state_tree().0;
        for candidate in &mut candidates {
            candidate.left.retain(|half_edge| {
                self.half_edge_id(*half_edge)
                    .is_none_or(|hedge| !ignored.includes(&hedge))
            });
            candidate.right.retain(|half_edge| {
                self.half_edge_id(*half_edge)
                    .is_none_or(|hedge| !ignored.includes(&hedge))
            });
        }
        candidates.sort_by_cached_key(|candidate| self.cut_order_key(&candidate.cut));
        candidates.dedup();
        self.topology_threshold_candidates = candidates;
        self.validate()?;
        self.id = DiagramId::from_key(self.model.fingerprint(), &self.structural_key()?)?;
        Ok(self)
    }

    /// Install topology-only threshold candidates from complete half-edge partitions.
    pub fn with_topology_threshold_partitions(
        self,
        partitions: Vec<(Vec<DiagramHalfEdge>, Vec<DiagramHalfEdge>)>,
    ) -> Result<Self, DiagramError> {
        let vertex_half_edges = self.vertex_half_edges();
        let ignored = self.initial_state_tree().0;
        let retained = |half_edge: &DiagramHalfEdge| {
            self.half_edge_id(*half_edge)
                .is_none_or(|hedge| !ignored.includes(&hedge))
        };
        let candidates = partitions
            .into_iter()
            .enumerate()
            .map(|(candidate, (mut left, mut right))| {
                left.retain(retained);
                right.retain(retained);
                self.cut_from_partitions(
                    candidate,
                    &left,
                    &right,
                    &vertex_half_edges,
                    Some(&ignored),
                )
                .map(|cut| DiagramThresholdCandidate {
                    cut: cut.cut,
                    left: cut.left.half_edges,
                    right: cut.right.half_edges,
                })
                .map_err(|error| match error {
                    DiagramError::InvalidCut { message, .. } => {
                        DiagramError::InvalidThresholdCandidate { candidate, message }
                    }
                    error => error,
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        self.with_topology_threshold_candidates(candidates)
    }

    /// Install physical cuts from authoritative left/right half-edge partitions.
    ///
    /// The oriented crossing endpoints, coupling-order summaries, and loop
    /// counts are derived from this diagram's finalized topology. This is the
    /// appropriate constructor when a valid partition is transported onto an
    /// isomorphic representative whose interaction assignment may differ.
    /// Supplied half-edge IDs must already use this diagram's edge-ID frame.
    pub fn with_cut_partitions(
        self,
        partitions: Vec<(Vec<DiagramHalfEdge>, Vec<DiagramHalfEdge>)>,
    ) -> Result<Self, DiagramError> {
        let vertex_half_edges = self.vertex_half_edges();
        let cuts = partitions
            .into_iter()
            .enumerate()
            .map(|(cut, (left, right))| {
                self.cut_from_partitions(cut, &left, &right, &vertex_half_edges, None)
            })
            .collect::<Result<Vec<_>, _>>()?;
        self.with_cuts(cuts)
    }

    /// Resolve a stable endpoint reference to its native Linnet half-edge.
    pub fn half_edge_id(
        &self,
        half_edge: DiagramHalfEdge,
    ) -> Option<linnet::half_edge::involution::Hedge> {
        if half_edge.edge.0 >= self.graph.n_edges() {
            return None;
        }
        let pair = self.graph[&EdgeIndex(half_edge.edge.0)].1;
        match (pair, half_edge.endpoint) {
            (
                HedgePair::Paired { source, .. } | HedgePair::Split { source, .. },
                DiagramEndpoint::Source,
            ) => Some(source),
            (
                HedgePair::Paired { sink, .. } | HedgePair::Split { sink, .. },
                DiagramEndpoint::Target,
            ) => Some(sink),
            (
                HedgePair::Unpaired {
                    hedge,
                    flow: Flow::Source,
                },
                DiagramEndpoint::Source,
            )
            | (
                HedgePair::Unpaired {
                    hedge,
                    flow: Flow::Sink,
                },
                DiagramEndpoint::Target,
            ) => Some(hedge),
            _ => None,
        }
    }

    /// Enumerate existing half-edges, in their native channel order.
    pub fn half_edges(&self) -> impl Iterator<Item = DiagramHalfEdge> + '_ {
        self.graph.iter_hedges().map(|(hedge, _)| DiagramHalfEdge {
            edge: EdgeId(self.graph[&hedge].0),
            endpoint: if self.graph.flow(hedge) == Flow::Source {
                DiagramEndpoint::Source
            } else {
                DiagramEndpoint::Target
            },
        })
    }

    fn vertex_half_edges(&self) -> Vec<Vec<DiagramHalfEdge>> {
        let mut result = vec![Vec::new(); self.graph.n_nodes()];
        for half_edge in self.half_edges() {
            let hedge = self.half_edge_id(half_edge).expect("existing half-edge");
            result[self.graph.node_id(hedge).0].push(half_edge);
        }
        result
    }

    fn summarize_cut_side(
        &self,
        half_edges: &BTreeSet<DiagramHalfEdge>,
        vertex_half_edges: &[Vec<DiagramHalfEdge>],
    ) -> Result<DiagramCutSide, DiagramError> {
        use linnet::half_edge::subgraph::ModifySubSet;
        let mut selected = self.graph.empty_subgraph::<SuBitGraph>();
        for half_edge in half_edges {
            if let Some(hedge) = self.half_edge_id(*half_edge) {
                selected.add(hedge);
            }
        }
        let mut coupling_orders = BTreeMap::new();
        for (vertex, incident) in vertex_half_edges.iter().enumerate() {
            if incident.iter().any(|edge| half_edges.contains(edge))
                && let Some(rule) = self.vertex(VertexId(vertex)).and_then(|v| v.interaction)
            {
                for (name, order) in self
                    .model
                    .vertex_rule_by_id(rule)?
                    .coupling_orders(&self.model)
                {
                    *coupling_orders.entry(name).or_insert(0) += order;
                }
            }
        }
        Ok(DiagramCutSide {
            half_edges: half_edges.iter().copied().collect(),
            coupling_orders,
            loop_count: self.graph.cyclotomatic_number(&selected),
        })
    }

    fn cut_from_partitions(
        &self,
        cut: usize,
        left_half_edges: &[DiagramHalfEdge],
        right_half_edges: &[DiagramHalfEdge],
        vertex_half_edges: &[Vec<DiagramHalfEdge>],
        ignored: Option<&SuBitGraph>,
    ) -> Result<DiagramCut, DiagramError> {
        let left = left_half_edges.iter().copied().collect::<BTreeSet<_>>();
        let right = right_half_edges.iter().copied().collect::<BTreeSet<_>>();
        if left.len() != left_half_edges.len() || right.len() != right_half_edges.len() {
            return Err(DiagramError::InvalidCut {
                cut,
                message: "half-edge sets contain duplicates".to_owned(),
            });
        }
        let is_ignored = |half_edge: &DiagramHalfEdge| {
            ignored.is_some_and(|ignored| {
                self.half_edge_id(*half_edge)
                    .is_some_and(|hedge| ignored.includes(&hedge))
            })
        };
        let universe = self
            .half_edges()
            .filter(|half_edge| !is_ignored(half_edge))
            .collect::<BTreeSet<_>>();
        if !left.is_disjoint(&right)
            || left.union(&right).copied().collect::<BTreeSet<_>>() != universe
        {
            return Err(DiagramError::InvalidCut {
                cut,
                message: format!(
                    "left and right half-edge partitions are not disjoint and complementary (missing {:?}, overlap {:?}, outside {:?})",
                    universe
                        .difference(&left.union(&right).copied().collect())
                        .collect::<Vec<_>>(),
                    left.intersection(&right).collect::<Vec<_>>(),
                    left.union(&right)
                        .filter(|half| !universe.contains(half))
                        .collect::<Vec<_>>()
                ),
            });
        }

        let mut oriented = BTreeSet::new();
        for (edge, endpoints, data) in self.edges() {
            if data.external.is_some() || endpoints.source.is_none() || endpoints.target.is_none() {
                continue;
            }
            let source = DiagramHalfEdge {
                edge,
                endpoint: DiagramEndpoint::Source,
            };
            let target = DiagramHalfEdge {
                edge,
                endpoint: DiagramEndpoint::Target,
            };
            // Initial-state normalization may remove one complete attachment
            // crown. A loop-carrying crossing retains its oriented boundary;
            // the loop-independent attachment tree itself contributes no cut.
            let crosses_ignored = self.loop_momentum_basis.edge_signatures[&edge]
                .loops
                .iter()
                .any(|sign| !sign.is_zero());
            let source_to_target = left.contains(&source) && right.contains(&target)
                || crosses_ignored
                    && (left.contains(&source) && is_ignored(&target)
                        || is_ignored(&source) && right.contains(&target));
            let target_to_source = left.contains(&target) && right.contains(&source)
                || crosses_ignored
                    && (left.contains(&target) && is_ignored(&source)
                        || is_ignored(&target) && right.contains(&source));
            match (source_to_target, target_to_source) {
                (true, false) => {
                    oriented.insert(source);
                }
                (false, true) => {
                    oriented.insert(target);
                }
                _ => {}
            }
        }
        for (vertex, incident) in vertex_half_edges.iter().enumerate() {
            if incident.iter().any(|edge| left.contains(edge))
                && incident.iter().any(|edge| right.contains(edge))
            {
                return Err(DiagramError::InvalidCut {
                    cut,
                    message: format!("vertex {vertex} is split between the two sides"),
                });
            }
        }

        Ok(DiagramCut {
            cut: oriented.into_iter().collect(),
            left: self.summarize_cut_side(&left, vertex_half_edges)?,
            right: self.summarize_cut_side(&right, vertex_half_edges)?,
        })
    }

    /// Select a spanning-forest routing by its ordered independent edges.
    pub fn with_loop_momentum_edges(mut self, requested: &[EdgeId]) -> Result<Self, DiagramError> {
        let requested_set: BTreeSet<_> = requested.iter().copied().collect();
        let internal = self.internal_subgraph();
        let tree_edges = self
            .graph
            .iter_edges_of(&internal)
            .filter_map(|(_, edge, _)| {
                (!requested_set.contains(&EdgeId(edge.0))).then_some(EdgeId(edge.0))
            })
            .collect::<Vec<_>>();
        self.loop_momentum_basis = self
            .basis_from_tree(&tree_edges)?
            .with_loop_edge_order(requested)?;
        Ok(self)
    }

    /// Select a spanning-forest routing by its ordered tree edges.
    ///
    /// The requested edges are validated and materialized directly. This is
    /// useful when a generator already chose a deterministic spanning forest:
    /// it does not enumerate the other spanning forests of the diagram.
    pub fn with_loop_momentum_tree_edges(
        mut self,
        requested: &[EdgeId],
    ) -> Result<Self, DiagramError> {
        self.loop_momentum_basis = self.basis_from_tree(requested)?;
        Ok(self)
    }

    /// Select the first depth-first basis for an explicit internal-edge order.
    ///
    /// Paired cross-section connections root the traversal at the attachment
    /// retained as hedge zero by legacy sewing (the outgoing side), in
    /// connection order. One-sided amplitude connections use their only
    /// attachment. No other spanning forests are enumerated. Generators can
    /// therefore reproduce a pre-existing half-edge insertion convention
    /// without moving topology or routing ownership into a downstream runtime.
    pub fn with_first_loop_momentum_basis_in_edge_order(
        self,
        ordered: &[EdgeId],
    ) -> Result<Self, DiagramError> {
        let internal = self
            .graph
            .iter_edges_of(&self.internal_subgraph())
            .map(|(_, edge, _)| EdgeId(edge.0))
            .collect::<BTreeSet<_>>();
        let ordered_set = ordered.iter().copied().collect::<BTreeSet<_>>();
        if ordered.len() != ordered_set.len() || ordered_set != internal {
            return Err(DiagramError::InvalidLoopMomentumBasis(format!(
                "ordered internal edges {ordered:?} do not match the diagram topology"
            )));
        }
        if ordered.is_empty() {
            return self.with_loop_momentum_tree_edges(&[]);
        }

        let mut builder = HedgeGraphBuilder::<EdgeId, ()>::new();
        let nodes = self
            .vertices()
            .map(|(vertex, _)| (vertex, builder.add_node(())))
            .collect::<BTreeMap<_, _>>();
        let endpoints = self
            .edges()
            .map(|(edge, endpoints, _)| (edge, endpoints))
            .collect::<BTreeMap<_, _>>();
        let mut external_attachments =
            BTreeMap::<usize, (Option<VertexId>, Option<VertexId>)>::new();
        for (_, endpoints, edge) in self.edges() {
            if let Some(external) = &edge.external {
                external_attachments
                    .insert(external.connection, (endpoints.source, endpoints.target));
            }
        }
        for edge in ordered {
            let endpoints = endpoints.get(edge).ok_or_else(|| DiagramError::Invariant {
                operation: "selecting the first ordered loop-momentum basis",
                message: format!("missing internal edge {}", edge.0),
            })?;
            builder.add_edge(
                nodes[&endpoints.source.expect("internal edge source")],
                nodes[&endpoints.target.expect("internal edge target")],
                *edge,
                Orientation::Undirected,
            );
        }
        let ordered_graph: HedgeGraph<EdgeId, ()> = builder.into();
        let full = ordered_graph.full_filter();
        let roots = external_attachments
            .values()
            .filter_map(|(incoming, outgoing)| (*outgoing).or(*incoming))
            .collect::<Vec<_>>();
        let mut remaining = full.clone();
        let mut forest: SuBitGraph = ordered_graph.empty_subgraph();
        let mut root_nodes = roots
            .into_iter()
            .filter_map(|vertex| nodes.get(&vertex).copied())
            .collect::<Vec<_>>();
        while !remaining.is_empty() {
            let root = root_nodes
                .iter()
                .position(|node| {
                    ordered_graph
                        .iter_crown(*node)
                        .any(|hedge| remaining.includes(&hedge))
                })
                .map(|position| root_nodes.remove(position))
                .unwrap_or_else(|| {
                    let hedge = remaining
                        .included_iter()
                        .next()
                        .expect("a non-empty subgraph has an included half-edge");
                    ordered_graph.node_id(hedge)
                });
            let tree =
                SimpleTraversalTree::depth_first_traverse(&ordered_graph, &full, &root, None)
                    .map_err(|error| {
                        DiagramError::InvalidLoopMomentumBasis(format!(
                            "ordered depth-first traversal failed: {error:?}"
                        ))
                    })?;
            forest.union_with(&tree.tree_subgraph);
            remaining.subtract_with(&tree.covers(&full));
        }
        let tree_edges = ordered_graph
            .iter_edges_of(&forest)
            .map(|(_, _, edge)| *edge.data)
            .collect::<Vec<_>>();
        self.with_loop_momentum_tree_edges(&tree_edges)
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
        let mut external_indices = BTreeSet::new();
        let mut external_connections = BTreeSet::new();
        let mut paired_external_edges = BTreeSet::new();
        let mut vertex_slots = vec![BTreeSet::new(); self.graph.n_nodes()];
        let vertex_half_edges = self.vertex_half_edges();
        for (edge_id, endpoints, edge) in self.edges() {
            self.resolve_particle(edge_id, edge)?;
            if let Some(external) = &edge.external {
                if !external_indices.insert(external.index) {
                    return Err(DiagramError::DuplicateExternalIndex(external.index));
                }
                if !external_connections.insert(external.connection) {
                    return Err(DiagramError::InvalidExternalConnectionSize {
                        connection: external.connection,
                        legs: 2,
                    });
                }
                if edge.numerator != Atom::one() {
                    return Err(DiagramError::ExternalEdgeNumerator { edge: edge_id.0 });
                }
                if endpoints.source.is_some() && endpoints.target.is_some() {
                    if external.state != ExternalState::Incoming {
                        return Err(DiagramError::Invariant {
                            operation: "validating sewn external flow",
                            message: format!(
                                "edge {} must carry the incoming external momentum convention",
                                edge_id.0
                            ),
                        });
                    }
                    paired_external_edges.insert(edge_id);
                } else if (external.state == ExternalState::Incoming) != endpoints.source.is_none()
                {
                    return Err(DiagramError::Invariant {
                        operation: "validating external flow",
                        message: format!("edge {} disagrees with its process state", edge_id.0),
                    });
                }
            } else if !edge.is_dummy && (endpoints.source.is_none() || endpoints.target.is_none()) {
                return Err(DiagramError::Invariant {
                    operation: "validating dangling external metadata",
                    message: format!("edge {} has no external-leg metadata", edge_id.0),
                });
            }
            for (vertex, slot) in [
                (endpoints.source, edge.source_slot),
                (endpoints.target, edge.target_slot),
            ] {
                if let Some(vertex) = vertex {
                    degrees[vertex.0] += 1;
                    if !vertex_slots[vertex.0].insert(slot.0) {
                        return Err(DiagramError::DuplicateVertexSlot {
                            vertex: vertex.0,
                            slot: slot.0,
                        });
                    }
                }
            }
            fragment_numerator *= &edge.numerator;
        }
        for (vertex_id, vertex) in self.vertices() {
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
            fragment_numerator *= &vertex.numerator;
            if let Some(interaction) = vertex.interaction {
                let rule = self.model.vertex_rule_by_id(interaction)?;
                let mut actual = vec![None; degrees[vertex_id.0]];
                for (edge_id, endpoints, edge) in self.edges() {
                    let (particle, base_id) = self.resolve_particle(edge_id, edge)?;
                    let base = self.model.particle_by_id(base_id)?;
                    if endpoints.source == Some(vertex_id) {
                        actual[edge.source_slot.0] = Some((
                            base.pdg_code,
                            edge.directed.then_some(!particle.is_antiparticle()),
                        ));
                    }
                    if endpoints.target == Some(vertex_id) {
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
        let has_paired_external_connection = !paired_external_edges.is_empty();
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
        if !has_paired_external_connection && !self.topology_threshold_candidates.is_empty() {
            return Err(DiagramError::InvalidThresholdCandidate {
                candidate: 0,
                message: "topology threshold candidates require paired incoming/outgoing external connections"
                    .to_owned(),
            });
        }

        for (cut_index, cut) in self.cuts.iter().enumerate() {
            let oriented = cut.cut.iter().copied().collect::<BTreeSet<_>>();
            if oriented.len() != cut.cut.len() {
                return Err(DiagramError::InvalidCut {
                    cut: cut_index,
                    message: "half-edge sets contain duplicates".to_owned(),
                });
            }
            let normalized = self.cut_from_partitions(
                cut_index,
                &cut.left.half_edges,
                &cut.right.half_edges,
                &vertex_half_edges,
                None,
            )?;
            if oriented != normalized.cut.iter().copied().collect::<BTreeSet<_>>() {
                return Err(DiagramError::InvalidCut {
                    cut: cut_index,
                    message: "oriented cut endpoints do not equal the crossing edges".to_owned(),
                });
            }
            if (normalized.left.coupling_orders, normalized.left.loop_count)
                != (cut.left.coupling_orders.clone(), cut.left.loop_count)
                || (
                    normalized.right.coupling_orders,
                    normalized.right.loop_count,
                ) != (cut.right.coupling_orders.clone(), cut.right.loop_count)
            {
                return Err(DiagramError::InvalidCut {
                    cut: cut_index,
                    message: "stored side coupling orders or loop counts do not match the finalized topology"
                        .to_owned(),
                });
            }
        }
        for (candidate_index, candidate) in self.topology_threshold_candidates.iter().enumerate() {
            let oriented = candidate.cut.iter().copied().collect::<BTreeSet<_>>();
            if oriented.len() != candidate.cut.len() {
                return Err(DiagramError::InvalidThresholdCandidate {
                    candidate: candidate_index,
                    message: "half-edge sets contain duplicates".to_owned(),
                });
            }
            let normalized = self
                .cut_from_partitions(
                    candidate_index,
                    &candidate.left,
                    &candidate.right,
                    &vertex_half_edges,
                    Some(&self.initial_state_tree().0),
                )
                .map_err(|error| match error {
                    DiagramError::InvalidCut { message, .. } => {
                        DiagramError::InvalidThresholdCandidate {
                            candidate: candidate_index,
                            message,
                        }
                    }
                    error => error,
                })?;
            if candidate
                .cut
                .iter()
                .filter(|half_edge| !paired_external_edges.contains(&half_edge.edge))
                .count()
                <= 1
            {
                return Err(DiagramError::InvalidThresholdCandidate {
                    candidate: candidate_index,
                    message: "threshold candidates must cross at least two non-initial-state edges"
                        .to_owned(),
                });
            }
            if oriented != normalized.cut.iter().copied().collect::<BTreeSet<_>>() {
                return Err(DiagramError::InvalidThresholdCandidate {
                    candidate: candidate_index,
                    message: "oriented cut endpoints do not equal the crossing edges".to_owned(),
                });
            }
        }
        if fragment_numerator != self.numerator
            && fragment_numerator.expand() != self.numerator.expand()
        {
            return Err(DiagramError::NumeratorFragmentMismatch);
        }
        self.loop_momentum_basis.validate(self)?;
        Ok(())
    }

    /// Serialize to a stable DOT dialect that can be parsed by [`Self::from_dot`].
    pub fn to_dot(&self) -> Result<String, DiagramError> {
        let mut output = String::new();
        let snapshot = self.serde_view();
        let loop_momentum_edges = self
            .loop_momentum_basis
            .loop_edges
            .iter()
            .map(|edge| edge.0.to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(output, "digraph feynkit {{")?;
        writeln!(
            output,
            "  graph [feynkit_name={}, model_fingerprint={}, symmetry_factor={}, overall_factor={}, numerator={}, numerator_prefactor={}, projector={}, loop_momentum_edges={}, loop_momentum_basis={}, half_edge_order={}, cuts={}, topology_threshold_candidates={}];",
            Self::dot_string(&self.name)?,
            Self::dot_string(&self.model.fingerprint().to_string())?,
            self.symmetry_factor,
            Self::dot_string(&snapshot.overall_factor)?,
            Self::dot_string(&snapshot.numerator)?,
            Self::dot_string(&snapshot.numerator_prefactor)?,
            Self::dot_string(&snapshot.projector)?,
            Self::dot_string(&loop_momentum_edges)?,
            Self::dot_string(&serde_json::to_string(&self.loop_momentum_basis)?)?,
            Self::dot_string(&serde_json::to_string(&snapshot.half_edge_order)?)?,
            Self::dot_string(&serde_json::to_string(&self.cuts)?)?,
            Self::dot_string(&serde_json::to_string(&self.topology_threshold_candidates)?)?,
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
                write!(
                    output,
                    ", interaction_id={}, interaction={}",
                    rule.index(),
                    Self::dot_string(&self.model.vertex_rule_by_id(rule)?.name)?
                )?;
            }
            writeln!(
                output,
                ", numerator={}];",
                Self::dot_string(&vertex.numerator.to_canonical_string())?
            )?;
        }
        for (id, endpoints, edge) in self.edges() {
            let source = endpoints
                .source
                .map(|vertex| format!("v{}", vertex.0))
                .unwrap_or_else(|| format!("ext{}", id.0));
            let target = endpoints
                .target
                .map(|vertex| format!("v{}", vertex.0))
                .unwrap_or_else(|| format!("ext{}", id.0));
            if endpoints.source.is_none() || endpoints.target.is_none() {
                // Linnet consumes invisible DOT endpoints as dangling half-edges,
                // never as interaction vertices in the imported graph.
                writeln!(output, "  ext{} [style=invis];", id.0)?;
            }
            let particle = self.model.particle_by_id(edge.particle)?;
            let orientation =
                linnet::half_edge::EdgeAccessors::orientation(&self.graph, EdgeIndex(id.0));
            let direction = match orientation {
                Orientation::Default => "forward",
                Orientation::Reversed => "back",
                Orientation::Undirected => "none",
            };
            writeln!(
                output,
                "  {source} -> {target} [id={}, particle_id={}, pdg={}, particle={}, directed={}, source_slot={}, target_slot={}, dir={}, external={}, is_dummy={}, numerator={}];",
                id.0,
                edge.particle.index(),
                particle.pdg_code,
                Self::dot_string(&particle.name)?,
                edge.directed,
                edge.source_slot.0,
                edge.target_slot.0,
                direction,
                Self::dot_string(&serde_json::to_string(&edge.external)?)?,
                edge.is_dummy,
                Self::dot_string(&edge.numerator.to_canonical_string())?,
            )?;
        }
        writeln!(output, "}}")?;
        Ok(output)
    }

    /// Parse compact model-aware DOT or the stable dialect emitted by [`Self::to_dot`].
    pub fn from_dot(model: impl Into<Arc<Model>>, input: &str) -> Result<Self, DiagramError> {
        let model = model.into();
        let parsed: DotGraph = DotGraph::from_string(input)
            .map_err(|error| DiagramError::DotParse(error.to_string()))?;
        if ![
            "model_fingerprint",
            "feynkit_name",
            "half_edge_order",
            "loop_momentum_basis",
            "cuts",
            "topology_threshold_candidates",
        ]
        .iter()
        .any(|key| parsed.global_data.statements.contains_key(*key))
            && !parsed
                .iter_edges()
                .any(|(_, _, edge)| edge.data.statements.contains_key("particle_id"))
        {
            return Self::from_compact_dot(model, parsed);
        }
        let serialized_fingerprint = parsed
            .global_data
            .statements
            .get("model_fingerprint")
            .ok_or_else(|| DiagramError::MissingDotAttribute {
                target: "graph".to_owned(),
                attribute: "model_fingerprint",
            })?;
        let actual_fingerprint = model.fingerprint();
        if serialized_fingerprint != &actual_fingerprint.to_string() {
            return Err(DiagramError::DotModelFingerprintMismatch {
                serialized: serialized_fingerprint.clone(),
                actual: actual_fingerprint,
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
                serde_json::from_str::<Vec<DiagramCut>>(cuts).map_err(|_| {
                    DiagramError::InvalidDotAttribute {
                        target: "graph".to_owned(),
                        attribute: "cuts",
                        value: cuts.clone(),
                    }
                })
            })?;
        let topology_threshold_candidates = parsed
            .global_data
            .statements
            .get("topology_threshold_candidates")
            .ok_or_else(|| DiagramError::MissingDotAttribute {
                target: "graph".to_owned(),
                attribute: "topology_threshold_candidates",
            })
            .and_then(|candidates| {
                serde_json::from_str::<Vec<DiagramThresholdCandidate>>(candidates).map_err(|_| {
                    DiagramError::InvalidDotAttribute {
                        target: "graph".to_owned(),
                        attribute: "topology_threshold_candidates",
                        value: candidates.clone(),
                    }
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
            .cuts(cuts)
            .topology_threshold_candidates(topology_threshold_candidates);
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
            if data.statements.contains_key("external_state") {
                return Err(DiagramError::Invariant {
                    operation: "reading finalized DOT vertices",
                    message: "external states must be represented by dangling or sewn edges, not vertices".into(),
                });
            }
            let vertex = DiagramVertex {
                name: vertex_name,
                interaction,
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
        for (_, _, data) in parsed.iter_edges() {
            if let Some(value) = data.data.statements.get("external")
                && let Some(external) = serde_json::from_str::<Option<ExternalLeg>>(value)?
            {
                next_external = next_external.max(
                    external
                        .index
                        .checked_add(1)
                        .ok_or(DiagramError::DotExternalIndexOverflow(external.index))?,
                );
            }
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
                external: data
                    .data
                    .statements
                    .get("external")
                    .map(|value| serde_json::from_str::<Option<ExternalLeg>>(value))
                    .transpose()?
                    .flatten(),
                is_dummy: data
                    .data
                    .statements
                    .get("is_dummy")
                    .map(|value| {
                        value
                            .parse::<bool>()
                            .map_err(|_| DiagramError::InvalidDotAttribute {
                                target: target.clone(),
                                attribute: "is_dummy",
                                value: value.clone(),
                            })
                    })
                    .transpose()?
                    .unwrap_or(false),
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
                    let source = *node_map
                        .get(&source_id)
                        .ok_or(DiagramError::UnknownVertex {
                            vertex: source_id.0,
                            vertices: node_map.len(),
                        })?;
                    let sink = *node_map.get(&sink_id).ok_or(DiagramError::UnknownVertex {
                        vertex: sink_id.0,
                        vertices: node_map.len(),
                    })?;
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
                    if edge.external.is_none() && !edge.is_dummy {
                        let index = next_external;
                        next_external = next_external
                            .checked_add(1)
                            .ok_or(DiagramError::DotExternalIndexOverflow(next_external))?;
                        edge.external = Some(ExternalLeg {
                            name: format!("ext{index}"),
                            index,
                            connection: index,
                            state: match flow {
                                Flow::Sink => ExternalState::Incoming,
                                Flow::Source => ExternalState::Outgoing,
                            },
                        });
                    }
                    let (source, target) = match flow {
                        Flow::Source => (Some(internal), None),
                        Flow::Sink => (None, Some(internal)),
                    };
                    builder.add_edge_with_slots(source, target, edge, source_slot, target_slot)?;
                }
            }
        }
        builder.edge_orientations = Some(
            parsed
                .iter_edges()
                .map(|(_, _, edge)| edge.orientation)
                .collect(),
        );
        builder.half_edge_order =
            Some(match parsed.global_data.statements.get("half_edge_order") {
                Some(order) => serde_json::from_str(order)?,
                None => parsed
                    .iter_hedges()
                    .map(|(hedge, _)| DiagramHalfEdge {
                        edge: EdgeId(parsed[&hedge].0),
                        endpoint: match parsed.flow(hedge) {
                            Flow::Source => DiagramEndpoint::Source,
                            Flow::Sink => DiagramEndpoint::Target,
                        },
                    })
                    .collect(),
            });
        if let Some(basis) = parsed.global_data.statements.get("loop_momentum_basis") {
            builder.loop_momentum_basis = Some(serde_json::from_str(basis)?);
        }
        let stored_basis = builder.loop_momentum_basis.is_some();
        let diagram = builder.build()?;
        let diagram = match loop_momentum_edges.filter(|_| !stored_basis) {
            Some(edges) => diagram.with_loop_momentum_edges(&edges),
            None => Ok(diagram),
        }?;
        diagram.validate()?;
        Ok(diagram)
    }

    /// Parse one or more stable FeynKit DOT diagrams from a single document.
    ///
    /// Each graph must use the dialect emitted by [`Self::to_dot`] and must
    /// carry the same model fingerprint as `model`. This is the canonical
    /// import path for a set of finalized cross-section diagrams because their
    /// typed physical cuts are retained graph by graph.
    pub fn from_dot_set(
        model: impl Into<Arc<Model>>,
        input: &str,
    ) -> Result<Vec<Self>, DiagramError> {
        let model = model.into();
        let parsed = DotGraphSet::from_string(input)
            .map_err(|error| DiagramError::DotParse(error.to_string()))?;
        parsed
            .set
            .into_iter()
            .zip(parsed.global_data)
            .map(|(graph, global_data)| {
                let graph = DotGraph { graph, global_data };
                Self::from_dot(Arc::clone(&model), &graph.debug_dot())
            })
            .collect()
    }

    pub fn loop_count(&self) -> usize {
        self.loop_momentum_basis.loop_edges.len()
    }

    /// Return the local superficial UV degree in `dimension` spacetime dimensions.
    ///
    /// This is `dimension * loops + sum(vertex degrees) + sum(edge degrees - 2)`.
    /// Degrees come from the stored local numerators; each internal propagator
    /// contributes its quadratic denominator. External legs, the projector, and
    /// diagram-wide prefactors are excluded, as in GammaLoop's local UV counting.
    /// All momenta at a vertex scale together. This is a superficial bound, not
    /// a test of cancellations in the contracted numerator or of subdivergences.
    /// Zero means logarithmic, positive means power divergent, and negative means
    /// superficially convergent.
    pub fn superficial_degree_of_divergence(&self, dimension: i32) -> Result<i32, DiagramError> {
        self.superficial_degree_of_divergence_of(&self.momentum_subgraph(), dimension)
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
        self.loop_momentum_bases_of(&self.graph.full_filter(), limit)
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
                    numerator: vertex.numerator.to_canonical_string(),
                })
                .collect(),
            edges: self
                .edges()
                .map(|(id, endpoints, edge)| {
                    (
                        endpoints,
                        DiagramEdgeSerde {
                            external: edge.external.clone(),
                            is_dummy: edge.is_dummy,
                            orientation: linnet::half_edge::EdgeAccessors::orientation(
                                &self.graph,
                                EdgeIndex(id.0),
                            ),
                            particle: edge.particle,
                            directed: edge.directed,
                            numerator: edge.numerator.to_canonical_string(),
                            source_slot: edge.source_slot,
                            target_slot: edge.target_slot,
                        },
                    )
                })
                .collect(),
            half_edge_order: self.half_edges().collect(),
            symmetry_factor: self.symmetry_factor,
            overall_factor: self.overall_factor.to_canonical_string(),
            numerator: self.numerator.to_canonical_string(),
            numerator_prefactor: self.numerator_prefactor.to_canonical_string(),
            projector: self.projector.to_canonical_string(),
            loop_momentum_basis: self.loop_momentum_basis.clone(),
            cuts: self.cuts.clone(),
            topology_threshold_candidates: self.topology_threshold_candidates.clone(),
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
            .cuts(data.cuts)
            .topology_threshold_candidates(data.topology_threshold_candidates);
        builder.half_edge_order = Some(data.half_edge_order);
        builder.edge_orientations = Some(
            data.edges
                .iter()
                .map(|(_, edge)| edge.orientation)
                .collect(),
        );
        for vertex in data.vertices {
            builder.add_vertex(DiagramVertex {
                name: vertex.name,
                interaction: vertex.interaction,
                numerator: Self::parse_expression("vertex numerator", vertex.numerator)?,
            });
        }
        for (endpoints, edge) in data.edges {
            builder.add_edge_with_slots(
                endpoints.source,
                endpoints.target,
                DiagramEdge {
                    external: edge.external,
                    is_dummy: edge.is_dummy,
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
        // Canonical strings include attributes and tags, but not Symbolica user data.
        // Reuse exact registered declarations so Spenso representation metadata survives.
        let mut symbols = State::symbol_iter()
            .filter(|(symbol, _)| !matches!(symbol.get_data(), UserData::None))
            .map(|(symbol, _)| (Atom::var(symbol).to_canonical_string().into(), symbol))
            .collect();
        let namespace = with_default_namespace!(&expression, "feynkit_graph");
        Workspace::get_local()
            .with(|workspace| {
                Token::parse_with_atom_info(
                    &expression,
                    ParseSettings::default(),
                    Some((&namespace, &mut symbols, workspace)),
                )?
                .to_atom(&namespace, &mut symbols, workspace)
            })
            .map_err(|message| DiagramError::SymbolicParse {
                field,
                expression,
                message,
            })
    }
}

/// Incremental owner for enforcing diagram-construction invariants.
pub struct FeynmanDiagramBuilder {
    model: Arc<Model>,
    name: String,
    vertices: Vec<DiagramVertex>,
    generation_externals: BTreeMap<VertexId, ExternalLeg>,
    edges: Vec<(EdgeEndpoints, DiagramEdge)>,
    half_edge_order: Option<Vec<DiagramHalfEdge>>,
    edge_orientations: Option<Vec<Orientation>>,
    symmetry_factor: u64,
    overall_factor: Atom,
    numerator: Atom,
    numerator_prefactor: Atom,
    projector: Atom,
    loop_momentum_basis: Option<LoopMomentumBasis>,
    cuts: Vec<DiagramCut>,
    topology_threshold_candidates: Vec<DiagramThresholdCandidate>,
}

impl FeynmanDiagramBuilder {
    pub fn new(model: impl Into<Arc<Model>>, name: impl Into<String>) -> Self {
        Self {
            model: model.into(),
            name: name.into(),
            vertices: Vec::new(),
            generation_externals: BTreeMap::new(),
            edges: Vec::new(),
            half_edge_order: None,
            edge_orientations: None,
            symmetry_factor: 1,
            overall_factor: Atom::one(),
            numerator: Atom::one(),
            numerator_prefactor: Atom::one(),
            projector: Atom::one(),
            loop_momentum_basis: None,
            cuts: Vec::new(),
            topology_threshold_candidates: Vec::new(),
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

    pub fn topology_threshold_candidates(
        mut self,
        candidates: Vec<DiagramThresholdCandidate>,
    ) -> Self {
        self.topology_threshold_candidates = candidates;
        self
    }

    pub fn add_vertex(&mut self, vertex: DiagramVertex) -> VertexId {
        let id = VertexId(self.vertices.len());
        self.vertices.push(vertex);
        id
    }

    /// Reserve a temporary Symbolica-generation endpoint, removed during finalization.
    /// Public diagrams carry this metadata on the resulting external edge.
    #[doc(hidden)]
    pub fn add_generation_external(
        &mut self,
        name: impl Into<String>,
        index: usize,
        state: ExternalState,
        connection: usize,
    ) -> VertexId {
        let name = name.into();
        let vertex = self.add_vertex(DiagramVertex {
            name: name.clone(),
            interaction: None,
            numerator: Atom::one(),
        });
        self.generation_externals.insert(
            vertex,
            ExternalLeg {
                name,
                index,
                state,
                connection,
            },
        );
        vertex
    }

    pub fn add_edge(
        &mut self,
        source: impl Into<Option<VertexId>>,
        target: impl Into<Option<VertexId>>,
        mut edge: DiagramEdge,
    ) -> Result<EdgeId, DiagramError> {
        let (source, target) = (source.into(), target.into());
        let degree = |vertex: Option<VertexId>| {
            self.edges
                .iter()
                .map(|(endpoints, _)| {
                    usize::from(vertex.is_some() && endpoints.source == vertex)
                        + usize::from(vertex.is_some() && endpoints.target == vertex)
                })
                .sum()
        };
        let source_slot = VertexSlot(degree(source));
        let target_slot =
            VertexSlot(degree(target) + usize::from(source.is_some() && source == target));
        edge.source_slot = source_slot;
        edge.target_slot = target_slot;
        self.add_edge_with_slots(source, target, edge, source_slot, target_slot)
    }

    pub fn add_edge_with_slots(
        &mut self,
        source: impl Into<Option<VertexId>>,
        target: impl Into<Option<VertexId>>,
        mut edge: DiagramEdge,
        source_slot: VertexSlot,
        target_slot: VertexSlot,
    ) -> Result<EdgeId, DiagramError> {
        let (source, target) = (source.into(), target.into());
        if source.is_none() && target.is_none() {
            return Err(DiagramError::Invariant {
                operation: "adding an edge",
                message: "an edge must have an attached half-edge".into(),
            });
        }
        for vertex in [source, target].into_iter().flatten() {
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

    pub fn build(self) -> Result<FeynmanDiagram, DiagramError> {
        self.finalize()
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

    #[test]
    fn compact_dot_resolves_particles_slots_numerators_ports_and_routing() {
        let model = scalar_model();
        let dot = r#"digraph triangle {
            graph [num="5"];
            ext [style=invis];
            a [num="2"];
            ext -> a:4 [id=0, pdg=25];
            b -> ext [id=1, particle="phi"];
            c -> ext [id=2, particle="phi"];
            a -> b [id=3, particle="phi", num="3", lmb_id=0];
            b -> c [id=4, particle="phi"];
            c -> a [id=5, particle="phi"];
        }"#;
        let diagram = FeynmanDiagram::from_dot(model.clone(), dot).unwrap();
        diagram.validate().unwrap();
        assert_eq!(diagram.loop_count(), 1);
        assert_eq!(diagram.numerator(), &Atom::num(6));
        assert_eq!(diagram.numerator_prefactor(), &Atom::num(5));
        assert_eq!(diagram.loop_momentum_basis().loop_edges, vec![EdgeId(3)]);
        assert!(
            diagram
                .vertices()
                .all(|(_, vertex)| vertex.interaction.is_some())
        );
        assert_eq!(
            diagram
                .half_edge_id(DiagramHalfEdge {
                    edge: EdgeId(0),
                    endpoint: DiagramEndpoint::Target
                })
                .unwrap()
                .0,
            4
        );
        assert_eq!(
            FeynmanDiagram::from_dot(model, &diagram.to_dot().unwrap())
                .unwrap()
                .to_json()
                .unwrap(),
            diagram.to_json().unwrap()
        );
    }

    #[test]
    fn compact_dot_sews_initial_states_and_selects_physical_cuts() {
        let model = scalar_model();
        let dot = r#"digraph bubble {
            graph [final_state="phi,phi"];
            ext [style=invis];
            ext -> a [particle="phi", is_cut=9];
            b -> ext [particle="phi", is_cut=9];
            a -> b [particle="phi", lmb_id=0];
            b -> a [particle="phi"];
        }"#;
        let diagram = FeynmanDiagram::from_dot(model.clone(), dot).unwrap();
        diagram.validate().unwrap();
        assert_eq!(diagram.loop_count(), 1);
        assert_eq!(diagram.cuts().len(), 1);
        assert_eq!(diagram.cuts()[0].cut.len(), 2);
        assert_eq!(diagram.cuts()[0].left.loop_count, 0);
        assert_eq!(diagram.cuts()[0].right.loop_count, 0);
        assert_eq!(diagram.topology_threshold_candidates().len(), 1);
        assert_eq!(
            diagram
                .edges()
                .filter(|(_, _, edge)| edge.external.is_some())
                .count(),
            1
        );
        assert_eq!(
            FeynmanDiagram::from_dot(model.clone(), &diagram.to_dot().unwrap())
                .unwrap()
                .to_json()
                .unwrap(),
            diagram.to_json().unwrap()
        );
        assert!(FeynmanDiagram::from_dot(model.clone(), &dot.replace("phi,phi", "phi")).is_err());
        assert!(
            FeynmanDiagram::from_dot(
                model.clone(),
                &dot.replace("graph [final_state=\"phi,phi\"];", "")
            )
            .is_err()
        );
        assert!(FeynmanDiagram::from_dot(model, &dot.replace("b -> ext", "ext -> b")).is_err());
    }

    #[test]
    fn compact_dot_rejects_inconsistent_model_data_and_loop_slots() {
        let model = scalar_model();
        for dot in [
            r#"digraph { ext [style=invis]; ext -> a [particle="missing"]; a -> ext [particle="phi"]; a -> ext [particle="phi"]; }"#,
            r#"digraph { ext [style=invis]; ext -> a [particle="phi"]; a -> ext [particle="phi"]; }"#,
            r#"digraph { ext [style=invis]; ext -> a [particle="phi", lmb_id=0]; a -> ext [particle="phi"]; a -> ext [particle="phi"]; }"#,
            r#"digraph { ext [style=invis]; ext -> a [particle="phi", pdg=999]; a -> ext [particle="phi"]; a -> ext [particle="phi"]; }"#,
        ] {
            assert!(
                FeynmanDiagram::from_dot(model.clone(), dot).is_err(),
                "{dot}"
            );
        }
    }

    #[test]
    fn canonical_dot_missing_fingerprint_does_not_become_compact() {
        let model = scalar_model();
        let diagram = FeynmanDiagram::from_dot(model.clone(), r#"digraph { ext [style=invis]; ext -> a [particle="phi"]; a -> ext [particle="phi"]; a -> ext [particle="phi"]; }"#).unwrap();
        let dot = diagram.to_dot().unwrap().replace(
            &format!("model_fingerprint=\"{}\", ", model.fingerprint()),
            "",
        );
        assert!(matches!(
            FeynmanDiagram::from_dot(model, &dot),
            Err(DiagramError::MissingDotAttribute {
                attribute: "model_fingerprint",
                ..
            })
        ));
    }

    fn external_edge(
        particle: ParticleId,
        name: &str,
        index: usize,
        state: ExternalState,
    ) -> DiagramEdge {
        let mut edge = DiagramEdge::new(particle, false);
        edge.external = Some(ExternalLeg {
            name: name.into(),
            index,
            state,
            connection: index,
        });
        edge
    }

    fn one_loop() -> FeynmanDiagram {
        let model = scalar_model();
        let rule = model.vertex_rule_id("V_1").unwrap();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "bubble");
        let left = builder.add_vertex(DiagramVertex::interaction("v0", rule));
        let right = builder.add_vertex(DiagramVertex::interaction("v1", rule));
        builder
            .add_edge(
                None,
                left,
                external_edge(particle, "p1", 0, ExternalState::Incoming),
            )
            .unwrap();
        builder
            .add_edge(left, right, DiagramEdge::new(particle, false))
            .unwrap();
        builder
            .add_edge(left, right, DiagramEdge::new(particle, false))
            .unwrap();
        builder
            .add_edge(
                right,
                None,
                external_edge(particle, "p2", 1, ExternalState::Outgoing),
            )
            .unwrap();
        builder.build().unwrap()
    }

    fn one_loop_reordered() -> FeynmanDiagram {
        let model = scalar_model();
        let rule = model.vertex_rule_id("V_1").unwrap();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "renamed");
        let right = builder.add_vertex(DiagramVertex::interaction("right", rule));
        let left = builder.add_vertex(DiagramVertex::interaction("left", rule));
        builder
            .add_edge_with_slots(
                right,
                None,
                external_edge(particle, "out", 1, ExternalState::Outgoing),
                VertexSlot(2),
                VertexSlot(0),
            )
            .unwrap();
        builder
            .add_edge_with_slots(
                left,
                right,
                DiagramEdge::new(particle, false),
                VertexSlot(1),
                VertexSlot(0),
            )
            .unwrap();
        builder
            .add_edge_with_slots(
                None,
                left,
                external_edge(particle, "in", 0, ExternalState::Incoming),
                VertexSlot(0),
                VertexSlot(0),
            )
            .unwrap();
        builder
            .add_edge_with_slots(
                left,
                right,
                DiagramEdge::new(particle, false),
                VertexSlot(2),
                VertexSlot(1),
            )
            .unwrap();
        builder.build().unwrap()
    }

    fn cut_scalar_line() -> FeynmanDiagram {
        let model = scalar_model();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "cut-scalar-line");
        let left = builder.add_vertex(DiagramVertex {
            name: "left".into(),
            interaction: None,
            numerator: Atom::one(),
        });
        let right = builder.add_vertex(DiagramVertex {
            name: "right".into(),
            interaction: None,
            numerator: Atom::one(),
        });
        let initial = builder
            .add_edge(
                left,
                right,
                external_edge(particle, "p", 0, ExternalState::Incoming),
            )
            .unwrap();
        let crossing = builder
            .add_edge(left, right, DiagramEdge::new(particle, false))
            .unwrap();
        let source = |edge| DiagramHalfEdge {
            edge,
            endpoint: DiagramEndpoint::Source,
        };
        let target = |edge| DiagramHalfEdge {
            edge,
            endpoint: DiagramEndpoint::Target,
        };
        let side = |half_edges| DiagramCutSide {
            half_edges,
            coupling_orders: BTreeMap::new(),
            loop_count: 0,
        };
        builder
            .cuts(vec![DiagramCut {
                cut: vec![target(crossing)],
                left: side(vec![target(initial), target(crossing)]),
                right: side(vec![source(initial), source(crossing)]),
            }])
            .build()
            .unwrap()
    }

    fn cut_scalar_bubble() -> FeynmanDiagram {
        let model = scalar_model();
        let rule = model.vertex_rule_id("V_1").unwrap();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "cut-scalar-bubble");
        let a = builder.add_vertex(DiagramVertex::interaction("left", rule));
        let b = builder.add_vertex(DiagramVertex::interaction("right", rule));
        let initial = builder
            .add_edge(
                a,
                b,
                external_edge(particle, "p", 0, ExternalState::Incoming),
            )
            .unwrap();
        let first = builder
            .add_edge(a, b, DiagramEdge::new(particle, false))
            .unwrap();
        let second = builder
            .add_edge(a, b, DiagramEdge::new(particle, false))
            .unwrap();
        let source = |edge| DiagramHalfEdge {
            edge,
            endpoint: DiagramEndpoint::Source,
        };
        let target = |edge| DiagramHalfEdge {
            edge,
            endpoint: DiagramEndpoint::Target,
        };
        let left = vec![target(initial), target(first), target(second)];
        let right = vec![source(initial), source(first), source(second)];
        let cut = vec![target(first), target(second)];
        let side = |half_edges| DiagramCutSide {
            half_edges,
            coupling_orders: BTreeMap::new(),
            loop_count: 0,
        };
        let diagram = builder
            .cuts(vec![DiagramCut {
                cut: cut.clone(),
                left: side(left.clone()),
                right: side(right.clone()),
            }])
            .topology_threshold_candidates(vec![DiagramThresholdCandidate { cut, left, right }])
            .build()
            .unwrap();
        diagram.validate().unwrap();
        diagram
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
    fn finalization_remaps_indices_inside_momentum_arguments() {
        let model = scalar_model();
        let particle = model.particle_id("phi").unwrap();
        let rule = model.vertex_rule_id("V_1").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "nested-momentum-indices");
        let external = builder.add_generation_external("in", 0, ExternalState::Incoming, 0);
        let vertex = builder.add_vertex(DiagramVertex::interaction("vertex", rule));
        builder
            .add_edge(vertex, vertex, DiagramEdge::new(particle, false))
            .unwrap();
        builder
            .add_edge(external, vertex, DiagramEdge::new(particle, false))
            .unwrap();
        let numerator = Atom::parse(
            "gammalooprs::Q(0,spenso::mink(4,gammalooprs::edge(0,1)))\
             *spenso::g(spenso::mink(4,gammalooprs::edge(0,1)),\
             spenso::mink(4,gammalooprs::hedge(0,1)))\
             +gammalooprs::Q(0,spenso::mink(4,gammalooprs::hedge(0,1)))",
            "feynkit_graph",
            ParseSettings::default(),
        )
        .unwrap();
        let diagram = builder.numerator(numerator).build().unwrap();
        let (edge, _, _) = diagram
            .edges()
            .find(|(_, _, edge)| edge.external.is_none())
            .unwrap();
        let hedge = diagram
            .half_edge_id(DiagramHalfEdge {
                edge,
                endpoint: DiagramEndpoint::Source,
            })
            .unwrap();
        assert_ne!(edge.0, 0);
        assert_ne!(hedge.0, 0);
        let expected = Atom::parse(
            format!(
                "gammalooprs::Q({e},spenso::mink(4,gammalooprs::edge({e},1)))\
                 *spenso::g(spenso::mink(4,gammalooprs::edge({e},1)),\
                 spenso::mink(4,gammalooprs::hedge({h},1)))\
                 +gammalooprs::Q({e},spenso::mink(4,gammalooprs::hedge({h},1)))",
                e = edge.0,
                h = hedge.0,
            ),
            "feynkit_graph",
            ParseSettings::default(),
        )
        .unwrap();
        assert_eq!(diagram.numerator(), &expected);
    }

    #[test]
    fn finalization_transports_signed_momentum_sums_once() {
        for coefficient in [1, 2] {
            let model = scalar_model();
            let particle = model.particle_id("phi").unwrap();
            let rule = model.vertex_rule_id("V_1").unwrap();
            let mut builder = FeynmanDiagram::builder(model, "signed-momentum-sum");
            let incoming = builder.add_generation_external("in", 0, ExternalState::Incoming, 0);
            let outgoing_1 = builder.add_generation_external("out1", 1, ExternalState::Outgoing, 1);
            let outgoing_2 = builder.add_generation_external("out2", 2, ExternalState::Outgoing, 2);
            let vertex = builder.add_vertex(DiagramVertex::interaction("vertex", rule));
            // All edge IDs change; the incoming edge also reverses its flow.
            builder
                .add_edge(vertex, outgoing_2, DiagramEdge::new(particle, false))
                .unwrap();
            builder
                .add_edge(vertex, incoming, DiagramEdge::new(particle, false))
                .unwrap();
            builder
                .add_edge(vertex, outgoing_1, DiagramEdge::new(particle, false))
                .unwrap();
            let numerator = Atom::parse(
                format!(
                    "-{coefficient}*gammalooprs::Q(1,spenso::mink(4,gammalooprs::hedge(2,1)))\
                     -gammalooprs::Q(0,spenso::mink(4,gammalooprs::edge(0,1)))\
                     +gammalooprs::Q(2,spenso::mink(4,gammalooprs::vertex(3,1)))"
                ),
                "feynkit_graph",
                ParseSettings::default(),
            )
            .unwrap();
            let diagram = builder.numerator(numerator).build().unwrap();
            let expected = Atom::parse(
                format!(
                    "{coefficient}*gammalooprs::Q(0,spenso::mink(4,gammalooprs::hedge(0,1)))\
                     -gammalooprs::Q(2,spenso::mink(4,gammalooprs::edge(2,1)))\
                     +gammalooprs::Q(1,spenso::mink(4,gammalooprs::vertex(0,1)))"
                ),
                "feynkit_graph",
                ParseSettings::default(),
            )
            .unwrap();
            assert_eq!(diagram.numerator(), &expected, "coefficient {coefficient}");
        }
    }

    #[test]
    fn uv_expansion_preserves_boundary_momenta() {
        let bubble = one_loop();
        let full = bubble
            .momentum_basis_of(&bubble.graph.full_filter())
            .unwrap();
        let selected = bubble
            .momentum_basis_of(&bubble.internal_subgraph())
            .unwrap();
        assert_eq!(selected.external_edges, full.external_edges);
        assert_eq!(selected.dependent_externals, full.dependent_externals);
        for edge in selected.external_edges {
            assert_eq!(selected.edge_signatures[&edge], full.edge_signatures[&edge]);
        }
    }

    #[test]
    fn uv_expansion_treats_an_omitted_chord_as_soft() {
        use linnet::half_edge::subgraph::ModifySubSet;
        let model = scalar_model();
        let particle = model.particle_id("phi").unwrap();
        let rule = model.vertex_rule_id("V_1").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "sunset");
        let left = builder.add_vertex(DiagramVertex::interaction("left", rule));
        let right = builder.add_vertex(DiagramVertex::interaction("right", rule));
        for _ in 0..3 {
            builder
                .add_edge(left, right, DiagramEdge::new(particle, false))
                .unwrap();
        }
        let diagram = builder.build().unwrap();
        let mut selected = diagram.internal_subgraph();
        selected.sub(diagram.graph[&EdgeIndex(2)].1);
        let basis = diagram.momentum_basis_of(&selected).unwrap();
        assert_eq!(basis.loop_edges.len(), 1);
        assert_eq!(basis.external_edges, vec![EdgeId(2)]);
        assert!(basis.dependent_externals.is_empty());
        let mass = Atom::var(symbolica::symbol!("feynkit_graph::test_mUV"));
        let soft = symbolica::function!(symbols::momentum(), 2, 0);
        let scalar = diagram
            .uv_expansion_of(&selected, &mass, 4, None, &BTreeMap::new())
            .unwrap();
        let tensor = diagram
            .uv_expansion_of(&selected, &mass, 4, Some(&soft), &BTreeMap::new())
            .unwrap();
        assert_eq!(tensor, &soft * scalar);
        let momenta = (0..3)
            .map(|edge| basis.route_expression(&symbols::momentum().call(edge)))
            .fold(Atom::Zero, |a, b| a + b);
        assert_eq!(momenta.expand(), Atom::Zero);
    }

    #[test]
    fn uv_expansion_retains_logarithmic_bubble_and_quadratic_mass_correction() {
        let bubble = one_loop();
        let mass = Atom::var(symbolica::symbol!("feynkit_graph::test_mUV"));
        let selected = bubble.internal_subgraph();
        let loop_edge = bubble.momentum_basis_of(&selected).unwrap().loop_edges[0];
        let q = symbols::momentum().call(loop_edge.0);
        let q2 = Minkowski {}
            .new_rep(symbols::dimension())
            .inner_product(&q, &q);
        let a = symbolica::symbol!("feynkit_graph::uv_a_");
        let b = symbolica::symbol!("feynkit_graph::uv_b_");
        let c = symbolica::symbol!("feynkit_graph::uv_c_");
        let d = symbolica::symbol!("feynkit_graph::uv_d_");
        let explicit = |atom: Atom| {
            atom.replace(symbolica::function!(symbols::denominator(), a, b, c, d))
                .with(d)
        };
        let expansion = bubble
            .uv_expansion_of(&selected, &mass, 4, None, &BTreeMap::new())
            .unwrap();
        assert_eq!(
            (explicit(expansion) - (&q2 - mass.pow(2)).pow(-2))
                .expand()
                .cancel(),
            Atom::Zero
        );
        assert_eq!(
            bubble
                .uv_expansion_of(&selected, &mass, 2, None, &BTreeMap::new())
                .unwrap(),
            Atom::Zero
        );
        assert_eq!(
            bubble
                .uv_expansion_of(&selected, &mass, 4, Some(&Atom::Zero), &BTreeMap::new())
                .unwrap(),
            Atom::Zero
        );

        // With all soft momenta set to zero, the six-dimensional bubble keeps
        // the quadratic leading term and the logarithmic physical-mass term.
        let expansion = bubble
            .uv_expansion_of(&selected, &mass, 6, None, &BTreeMap::new())
            .unwrap();
        let basis = bubble.momentum_basis_of(&selected).unwrap();
        // Before setting external momenta to zero, the subtraction must remove
        // every divergent coefficient of the original, undeformed integrand.
        let bare = Atom::one() / bubble.denominator_of(&selected, &BTreeMap::new()).unwrap();
        let scale = symbolica::symbol!("feynkit_graph::uv_test_scale"; Scalar);
        let args = symbolica::symbol!("feynkit_graph::uv_test_args___");
        let loop_momentum = symbolica::function!(symbols::loop_momentum(), args);
        let remainder = basis
            .route_expression(&explicit(bare - &expansion))
            .replace(loop_momentum.to_pattern())
            .with(&loop_momentum / scale)
            * Atom::var(scale).pow(-6);
        assert_eq!(
            remainder
                .series(scale, Atom::Zero, 0)
                .unwrap()
                .to_atom()
                .expand()
                .cancel(),
            Atom::Zero
        );
        let mut expansion = explicit(expansion);
        for edge in basis.external_edges {
            expansion = expansion
                .replace(symbolica::function!(symbols::momentum(), edge.0, args))
                .with(Atom::Zero);
        }
        let vacuum = &q2 - mass.pow(2);
        let expected = vacuum.pow(-2)
            + Atom::num(2)
                * (Atom::var(symbolica::symbol!("UFO::M")).pow(2) - mass.pow(2))
                * vacuum.pow(-3);
        assert_eq!((expansion - expected).expand().cancel(), Atom::Zero);
    }

    #[test]
    fn uv_expansion_matches_taylor_coefficients_for_signed_edge_powers() {
        let bubble = one_loop();
        let selected = bubble.internal_subgraph();
        let basis = bubble.momentum_basis_of(&selected).unwrap();
        let loop_edge = basis.loop_edges[0];
        let shifted_edge = basis.tree_edges[0];
        let external_edge = basis.external_edges[0];
        let k = symbols::momentum().call(loop_edge.0);
        let p = symbols::momentum().call(external_edge.0);
        assert_eq!(
            basis.route_expression(&symbols::momentum().call(shifted_edge.0)),
            basis.route_expression(&(&p - &k))
        );
        let metric = Minkowski {}.new_rep(symbols::dimension());
        let k2 = metric.inner_product(&k, &k);
        let kp = metric.inner_product(&k, &p);
        let p2 = metric.inner_product(&p, &p);
        let mass = Atom::var(symbolica::symbol!("feynkit_graph::test_mUV"));
        let mass_difference = Atom::var(symbolica::symbol!("UFO::M")).pow(2) - mass.pow(2);
        let vacuum = &k2 - mass.pow(2);
        let a = symbolica::symbol!("feynkit_graph::uv_a_");
        let b = symbolica::symbol!("feynkit_graph::uv_b_");
        let c = symbolica::symbol!("feynkit_graph::uv_c_");
        let d = symbolica::symbol!("feynkit_graph::uv_d_");
        for (hard_power, shifted_power, dimension) in [(1, 2, 8), (1, 0, 4), (2, -1, 4)] {
            let powers = BTreeMap::from([(loop_edge, hard_power), (shifted_edge, shifted_power)]);
            let expanded = bubble
                .uv_expansion_of(&selected, &mass, dimension, None, &powers)
                .unwrap()
                .replace(symbolica::function!(symbols::denominator(), a, b, c, d))
                .with(d);
            // Taylor-expand (V + mUV² - m²)^-a [V - 2 k.p + p² + mUV² - m²]^-b
            // through degree two. Keep the odd term before tensor integration.
            let n = (hard_power + shifted_power) as i64;
            let shifted_power = shifted_power as i64;
            let expected = vacuum.pow(-n)
                + Atom::num(2 * shifted_power) * &kp * vacuum.pow(-n - 1)
                + (Atom::num(n) * &mass_difference - Atom::num(shifted_power) * &p2)
                    * vacuum.pow(-n - 1)
                + Atom::num(2 * shifted_power * (shifted_power + 1))
                    * kp.pow(2)
                    * vacuum.pow(-n - 2);
            assert_eq!((expanded - expected).expand().cancel(), Atom::Zero);
        }
        let powers = BTreeMap::from([(shifted_edge, 2)]);
        assert_eq!(
            bubble
                .uv_expansion_of(&selected, &mass, 4, None, &powers)
                .unwrap(),
            Atom::Zero
        );
        let logarithmic = bubble
            .uv_expansion_of(&selected, &mass, 6, None, &powers)
            .unwrap();
        assert_eq!(
            (logarithmic
                .clone()
                .replace(symbolica::function!(symbols::denominator(), a, b, c, d))
                .with(d)
                - vacuum.pow(-3))
            .expand()
            .cancel(),
            Atom::Zero
        );
        // As for denominator_of, powers for non-selected or unknown edges are ignored.
        let mut unused_powers = powers.clone();
        unused_powers.insert(external_edge, -2);
        unused_powers.insert(EdgeId(usize::MAX), 3);
        assert_eq!(
            bubble.denominator_of(&selected, &unused_powers).unwrap(),
            bubble.denominator_of(&selected, &powers).unwrap()
        );
        assert_eq!(
            bubble
                .uv_expansion_of(&selected, &mass, 6, None, &unused_powers)
                .unwrap(),
            logarithmic
        );
    }

    #[test]
    fn uv_expansion_uses_only_selected_loops_and_preserves_tensor_numerators() {
        use linnet::half_edge::subgraph::ModifySubSet;
        let model = scalar_model();
        let particle = model.particle_id("phi").unwrap();
        let rule = model.vertex_rule_id("V_1").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "two-tadpoles");
        let left = builder.add_vertex(DiagramVertex::interaction("left", rule));
        let right = builder.add_vertex(DiagramVertex::interaction("right", rule));
        builder
            .add_edge(left, left, DiagramEdge::new(particle, false))
            .unwrap();
        builder
            .add_edge(left, right, DiagramEdge::new(particle, false))
            .unwrap();
        builder
            .add_edge(right, right, DiagramEdge::new(particle, false))
            .unwrap();
        let diagram = builder.build().unwrap();
        let loop_edge = diagram
            .edges()
            .find(|(_, ends, _)| ends.source == ends.target)
            .unwrap()
            .0;
        let mut selected = diagram.graph.empty_subgraph::<SuBitGraph>();
        selected.add(diagram.graph[&EdgeIndex(loop_edge.0)].1);
        assert_eq!(
            diagram
                .momentum_basis_of(&selected)
                .unwrap()
                .loop_edges
                .len(),
            1
        );
        let mass = Atom::var(symbolica::symbol!("feynkit_graph::test_mUV"));
        let mass2 = mass.pow(2);
        let vacuum = expressions::PropagatorSymbols {
            momentum: symbols::momentum(),
            denominator: symbols::denominator(),
        }
        .denominator(EdgeIndex(loop_edge.0), &mass2, symbols::dimension());
        let expected = vacuum.pow(-1)
            + (Atom::var(symbolica::symbol!("UFO::M")).pow(2) - mass2) * vacuum.pow(-2);
        assert_eq!(
            (diagram
                .uv_expansion_of(&selected, &mass, 4, None, &BTreeMap::new())
                .unwrap()
                - expected)
                .expand(),
            Atom::Zero
        );
        let index = symbolica::symbol!("feynkit_graph::uv_test_mu");
        let tensor = symbolica::function!(
            symbols::momentum(),
            loop_edge.0,
            symbolica::function!(symbolica::symbol!("spenso::mink"), 4, index)
        );
        assert_eq!(
            diagram
                .uv_expansion_of(&selected, &mass, 2, Some(&tensor), &BTreeMap::new())
                .unwrap(),
            &tensor / &vacuum
        );
        assert_eq!(
            diagram
                .uv_expansion_of(
                    &diagram.graph.empty_subgraph(),
                    &mass,
                    4,
                    None,
                    &BTreeMap::new()
                )
                .unwrap(),
            Atom::Zero
        );
        assert!(
            diagram
                .uv_expansion_of(&selected, &mass, 0, None, &BTreeMap::new())
                .is_err()
        );
    }

    #[test]
    fn superficial_uv_degree_counts_local_numerators_and_internal_edges() {
        let bubble = one_loop();
        assert_eq!(bubble.superficial_degree_of_divergence(4).unwrap(), 0);
        assert_eq!(bubble.superficial_degree_of_divergence(6).unwrap(), 2);
        assert_eq!(bubble.superficial_degree_of_divergence(2).unwrap(), -2);
        let internal_edges = bubble
            .underlying()
            .iter_edges_of(&bubble.internal_subgraph())
            .map(|(_, edge, _)| EdgeId(edge.0))
            .collect::<Vec<_>>();
        let mut with_momenta = bubble
            .map_data(
                |_, vertex| {
                    let mut vertex = vertex.clone();
                    vertex.numerator =
                        symbolica::function!(momentum_symbol(), internal_edges[0].0, 0)
                            + symbolica::symbol!("feynkit_graph::uv_test_mass");
                    vertex
                },
                |id, _, edge| {
                    let mut edge = edge.clone();
                    if internal_edges.contains(&id) {
                        edge.numerator = symbolica::function!(momentum_symbol(), id.0, 0)
                            + symbolica::symbol!("feynkit_graph::uv_test_mass");
                    }
                    edge
                },
            )
            .unwrap();
        with_momenta.numerator = with_momenta
            .vertices()
            .map(|(_, vertex)| &vertex.numerator)
            .chain(with_momenta.edges().map(|(_, _, edge)| &edge.numerator))
            .fold(Atom::num(1), |product, numerator| product * numerator);
        assert_eq!(with_momenta.superficial_degree_of_divergence(4).unwrap(), 4);
        let restored =
            FeynmanDiagram::from_json(with_momenta.model_arc(), &with_momenta.to_json().unwrap())
                .unwrap();
        assert_eq!(restored.superficial_degree_of_divergence(4).unwrap(), 4);
    }

    #[test]
    fn diagram_integral_family_matches_bubble_for_every_routing() {
        let bubble = one_loop();
        let generic = bubble
            .integral_family(&feynkit_kinematics::Kinematics::new())
            .unwrap();
        assert_eq!(generic.loop_momenta().len(), 1);
        assert_eq!(generic.external_momenta().len(), 1);
        assert_eq!(generic.denominators().len(), 2);
        let kin = generic
            .kinematics()
            .clone()
            .with_mass_squared(&generic.external_momenta()[0], symbolica::parse!("s"))
            .unwrap();
        let reference = bubble.integral_family(&kin).unwrap();
        for basis in bubble.loop_momentum_bases().unwrap() {
            let diagram = bubble
                .clone()
                .with_loop_momentum_edges(&basis.loop_edges)
                .unwrap();
            let family = diagram.integral_family(&kin).unwrap();
            assert!(family.is_complete() && family.is_independent());
            let (u, f) = family
                .symanzik(&[symbolica::parse!("x"), symbolica::parse!("y")])
                .unwrap();
            assert_eq!(u, symbolica::parse!("x+y"));
            assert!(
                (f - symbolica::parse!("UFO::M^2*(x+y)^2-s*x*y"))
                    .expand()
                    .is_zero()
            );
            assert!(family.find_mapping(&reference, 100).unwrap().is_some());
        }
        assert!(matches!(
            directed_fermion_line().integral_family(&kin),
            Err(DiagramError::IntegralFamily(IntegralFamilyError::NoLoops))
        ));
    }

    #[test]
    fn denominator_uses_only_internal_edges_and_symbolic_masses() {
        let bubble = one_loop();
        let q1 = symbolica::function!(momentum_symbol(), 1);
        let q2 = symbolica::function!(momentum_symbol(), 2);
        let mass2 = Atom::var(symbolica::symbol!("UFO::M")).pow(2);
        let metric = Minkowski {}.new_rep(symbols::dimension());
        let expected = symbolica::function!(
            symbols::denominator(),
            1,
            &q1,
            &mass2,
            metric.inner_product(&q1, &q1) - &mass2
        ) * symbolica::function!(
            symbols::denominator(),
            2,
            &q2,
            &mass2,
            metric.inner_product(&q2, &q2) - &mass2
        );
        assert_eq!(
            bubble
                .underlying()
                .iter_edges_of(&bubble.internal_subgraph())
                .map(|(_, edge, _)| EdgeId(edge.0))
                .collect::<Vec<_>>()
                .len(),
            2
        );
        assert_eq!(bubble.denominator_expression().unwrap(), expected);
        let restored =
            FeynmanDiagram::from_json(bubble.model_arc(), &bubble.to_json().unwrap()).unwrap();
        assert_eq!(restored.denominator_expression().unwrap(), expected);
        let contact = FeynmanDiagram::builder(scalar_model(), "empty")
            .build()
            .unwrap();
        assert_eq!(contact.denominator_expression().unwrap(), Atom::num(1));
    }

    #[test]
    fn json_and_dot_round_trip() {
        let diagram = one_loop();
        let from_json =
            FeynmanDiagram::from_json(diagram.model_arc(), &diagram.to_json().unwrap()).unwrap();
        assert_eq!(from_json.loop_count(), 1);
        assert_eq!(from_json.vertices().count(), 2);

        let dot = diagram.to_dot().unwrap();
        assert!(
            dot.lines().all(|line| line.trim_end() == line),
            "canonical DOT must not contain trailing whitespace"
        );
        let from_dot = FeynmanDiagram::from_dot(diagram.model_arc(), &dot).unwrap();
        assert_eq!(from_dot.name(), "bubble");
        assert_eq!(from_dot.loop_count(), 1);
        assert_eq!(from_dot.edges().count(), 4);
    }

    #[test]
    fn dot_set_round_trip() {
        let first = cut_scalar_line().with_name("first");
        let second = cut_scalar_line().with_name("second");
        let input = format!("{}\n{}", first.to_dot().unwrap(), second.to_dot().unwrap());
        let diagrams = FeynmanDiagram::from_dot_set(first.model_arc(), &input).unwrap();

        assert_eq!(
            diagrams
                .iter()
                .map(|diagram| diagram.name())
                .collect::<Vec<_>>(),
            ["first", "second"]
        );
        assert!(
            diagrams
                .iter()
                .all(|diagram| diagram.cuts() == first.cuts())
        );
    }

    #[test]
    fn validates_factored_and_expanded_numerators_and_rejects_mismatches() {
        let factored =
            Atom::parse("(x+y)*(x-y)", "validation_test", ParseSettings::default()).unwrap();
        let mut diagram = one_loop().with_numerator(factored.clone()).unwrap();
        diagram.validate().unwrap();
        diagram.numerator = factored.expand();
        diagram.validate().unwrap();
        diagram.numerator += Atom::one();
        assert!(matches!(
            diagram.validate(),
            Err(DiagramError::NumeratorFragmentMismatch)
        ));
    }

    #[test]
    fn replacing_a_numerator_preserves_identity_and_fragment_invariants() {
        let diagram = one_loop();
        let topology_id = diagram.id();
        let replaced = diagram.with_numerator(Atom::num(7)).unwrap();

        assert_eq!(replaced.id(), topology_id);
        assert_eq!(replaced.numerator(), &Atom::num(7));
        assert_eq!(
            replaced
                .vertices()
                .filter(|(_, vertex)| vertex.numerator != Atom::one())
                .count(),
            1
        );
        assert!(
            replaced
                .edges()
                .all(|(_, _, edge)| edge.numerator == Atom::one())
        );
        replaced.validate().unwrap();
        let from_json =
            FeynmanDiagram::from_json(replaced.model_arc(), &replaced.to_json().unwrap()).unwrap();
        let from_dot =
            FeynmanDiagram::from_dot(replaced.model_arc(), &replaced.to_dot().unwrap()).unwrap();
        assert_eq!(from_json.numerator(), &Atom::num(7));
        assert_eq!(from_dot.numerator(), &Atom::num(7));
    }

    #[test]
    fn replacing_a_numerator_requires_a_non_external_anchor() {
        let diagram = FeynmanDiagram::builder(scalar_model(), "empty")
            .build()
            .unwrap();
        assert!(matches!(
            diagram.with_numerator(Atom::num(7)),
            Err(DiagramError::MissingNumeratorAnchor)
        ));
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

        let source = diagram
            .to_linnest(None, &Default::default(), None, "(:)")
            .unwrap();
        assert!(source.contains("mode: \"cross-section\""));
        assert_eq!(source.matches("is_cut: 0").count(), 1);

        let mut invalid_cut = diagram.cuts()[0].clone();
        invalid_cut.cut.clear();
        assert!(matches!(
            diagram.with_cuts(vec![invalid_cut]),
            Err(DiagramError::InvalidCut { cut: 0, message })
                if message.contains("crossing edges")
        ));
    }

    #[test]
    fn validates_and_round_trips_topology_threshold_candidates() {
        let diagram = cut_scalar_bubble();
        let expected = diagram.topology_threshold_candidates().to_vec();
        assert_eq!(expected.len(), 1);
        assert_eq!(expected[0].cut.len(), 2);

        let from_json =
            FeynmanDiagram::from_json(diagram.model_arc(), &diagram.to_json().unwrap()).unwrap();
        assert_eq!(from_json.topology_threshold_candidates(), expected);

        let from_dot =
            FeynmanDiagram::from_dot(diagram.model_arc(), &diagram.to_dot().unwrap()).unwrap();
        assert_eq!(from_dot.topology_threshold_candidates(), expected);

        let rebuilt = diagram
            .clone()
            .with_topology_threshold_partitions(vec![(
                expected[0].left.clone(),
                expected[0].right.clone(),
            )])
            .unwrap();
        assert_eq!(rebuilt.topology_threshold_candidates(), expected);

        let without_candidates = diagram
            .clone()
            .with_topology_threshold_candidates(Vec::new())
            .unwrap();
        assert_ne!(diagram.id(), without_candidates.id());
        assert_ne!(
            diagram.canonical_key().unwrap(),
            without_candidates.canonical_key().unwrap()
        );
    }

    #[test]
    fn rejects_invalid_topology_threshold_candidates() {
        let diagram = cut_scalar_bubble();
        let source = diagram.topology_threshold_candidates()[0].clone();

        let mut incomplete = source.clone();
        incomplete.right.pop();
        assert!(matches!(
            diagram
                .clone()
                .with_topology_threshold_candidates(vec![incomplete]),
            Err(DiagramError::InvalidThresholdCandidate {
                candidate: 0,
                message
            }) if message.contains("disjoint and complementary")
        ));

        let mut wrong_crossing = source.clone();
        wrong_crossing.cut.pop();
        assert!(matches!(
            diagram
                .clone()
                .with_topology_threshold_candidates(vec![wrong_crossing]),
            Err(DiagramError::InvalidThresholdCandidate {
                candidate: 0,
                message
            }) if message.contains("at least two non-initial-state edges")
        ));

        let mut split = source;
        let moved = split.left.pop().unwrap();
        split.right.push(moved);
        assert!(matches!(
            diagram
                .clone()
                .with_topology_threshold_candidates(vec![split]),
            Err(DiagramError::InvalidThresholdCandidate {
                candidate: 0,
                message
            }) if message.contains("is split")
        ));

        let candidate = diagram.topology_threshold_candidates()[0].clone();
        assert!(matches!(
            one_loop().with_topology_threshold_candidates(vec![candidate]),
            Err(DiagramError::InvalidThresholdCandidate {
                candidate: 0,
                message
            }) if message.contains("paired incoming/outgoing")
        ));

        let line = cut_scalar_line();
        let source = |edge| DiagramHalfEdge {
            edge: EdgeId(edge),
            endpoint: DiagramEndpoint::Source,
        };
        let target = |edge| DiagramHalfEdge {
            edge: EdgeId(edge),
            endpoint: DiagramEndpoint::Target,
        };
        assert!(matches!(
            line.with_topology_threshold_candidates(vec![DiagramThresholdCandidate {
                cut: vec![target(1)],
                left: vec![target(0), target(1)],
                right: vec![source(0), source(1)],
            }]),
            Err(DiagramError::InvalidThresholdCandidate {
                candidate: 0,
                message
            }) if message.contains("non-initial-state")
        ));
    }

    #[test]
    fn cut_partitions_recompute_crossings_and_side_summaries() {
        let diagram = cut_scalar_line();
        let source = diagram.cuts()[0].clone();
        let rebuilt = diagram
            .clone()
            .with_cut_partitions(vec![(
                source.left.half_edges.clone(),
                source.right.half_edges.clone(),
            )])
            .unwrap();
        assert_eq!(rebuilt.cuts(), diagram.cuts());
        rebuilt.validate().unwrap();

        let mirrored = diagram
            .with_cut_partitions(vec![(
                source.right.half_edges.clone(),
                source.left.half_edges.clone(),
            )])
            .unwrap();
        assert_eq!(mirrored.cuts()[0].cut.len(), source.cut.len());
        assert_ne!(mirrored.cuts()[0].cut, source.cut);
        assert!(mirrored.cuts()[0].left.coupling_orders.is_empty());
        assert_eq!(mirrored.cuts()[0].left.loop_count, 0);
        mirrored.validate().unwrap();
    }

    #[test]
    fn cut_partitions_reject_incomplete_and_split_vertices() {
        let diagram = cut_scalar_line();
        let source = diagram.cuts()[0].clone();
        let mut incomplete_right = source.right.half_edges.clone();
        incomplete_right.pop();
        assert!(matches!(
            diagram.clone().with_cut_partitions(vec![(
                source.left.half_edges.clone(),
                incomplete_right,
            )]),
            Err(DiagramError::InvalidCut { cut: 0, message })
                if message.contains("disjoint and complementary")
        ));

        let mut split_left = source.left.half_edges.clone();
        let moved = split_left.remove(1);
        let mut split_right = source.right.half_edges.clone();
        split_right.push(moved);
        assert!(matches!(
            diagram.with_cut_partitions(vec![(split_left, split_right)]),
            Err(DiagramError::InvalidCut { cut: 0, message })
                if message.contains("is split")
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
            Err(DiagramError::UnknownEdge { .. } | DiagramError::Invariant { .. })
        ));

        let mut external_edge = diagram.serde_view();
        external_edge.edges[0].0.target = None;
        assert!(matches!(
            FeynmanDiagram::from_json(
                diagram.model_arc(),
                &serde_json::to_string(&external_edge).unwrap()
            ),
            Err(error) if error.to_string().contains("attached half-edge")
        ));
    }

    #[test]
    fn native_external_metadata_matches_dangling_and_sewn_roles() {
        let amplitude = one_loop()
            .map_data(
                |_, vertex| vertex.clone(),
                |_, _, edge| {
                    let mut edge = edge.clone();
                    edge.external = None;
                    edge
                },
            )
            .unwrap();
        assert!(matches!(
            amplitude.validate(),
            Err(DiagramError::Invariant {
                operation: "validating dangling external metadata",
                ..
            })
        ));

        let sewn = cut_scalar_line()
            .map_data(
                |_, vertex| vertex.clone(),
                |_, _, edge| {
                    let mut edge = edge.clone();
                    if let Some(external) = &mut edge.external {
                        external.state = ExternalState::Outgoing;
                    }
                    edge
                },
            )
            .unwrap();
        assert!(matches!(
            sewn.validate(),
            Err(DiagramError::Invariant {
                operation: "validating sewn external flow",
                ..
            })
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
        let metadata = diagram.edges().next().unwrap().2.external.clone();
        let mut overflow = metadata.clone();
        overflow.as_mut().unwrap().index = usize::MAX;
        let dot = diagram.to_dot().unwrap().replacen(
            &FeynmanDiagram::dot_string(&serde_json::to_string(&metadata).unwrap()).unwrap(),
            &FeynmanDiagram::dot_string(&serde_json::to_string(&overflow).unwrap()).unwrap(),
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

        let direct = diagram
            .with_loop_momentum_tree_edges(&[EdgeId(0), EdgeId(1), EdgeId(3)])
            .unwrap();
        assert_eq!(
            direct.loop_momentum_basis().tree_edges,
            vec![EdgeId(0), EdgeId(1), EdgeId(3)]
        );
        assert_eq!(direct.loop_momentum_basis().loop_edges, vec![EdgeId(2)]);
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
                |_, vertex| vertex.clone(),
                |_, _, edge| {
                    let mut edge = edge.clone();
                    if let Some(external) = &mut edge.external
                        && external.index == 0
                    {
                        external.connection = 2;
                    }
                    edge
                },
            )
            .unwrap();
        assert_ne!(key, changed_state.canonical_key().unwrap());

        let changed_interaction = diagram
            .map_data(
                |id, vertex| {
                    let mut vertex = vertex.clone();
                    if id == VertexId(0) {
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
        assert_eq!(endpoints.source, Some(VertexId(1)));
        assert_eq!(endpoints.target, Some(VertexId(0)));
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
        let indices = relabeled
            .edges()
            .filter_map(|(_, _, edge)| edge.external.as_ref().map(|leg| leg.index))
            .collect::<Vec<_>>();
        assert_eq!(indices, vec![1, 0]);
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
        let model = scalar_model();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "invalid");
        let vertex = builder.add_vertex(DiagramVertex {
            name: "v".into(),
            interaction: None,
            numerator: Atom::one(),
        });
        builder
            .add_edge(
                None,
                vertex,
                external_edge(particle, "p1", 0, ExternalState::Incoming),
            )
            .unwrap();
        builder
            .add_edge(
                vertex,
                None,
                external_edge(particle, "p2", 0, ExternalState::Outgoing),
            )
            .unwrap();
        assert!(matches!(
            builder.build().and_then(|diagram| diagram.validate()),
            Err(DiagramError::DuplicateExternalIndex(0))
        ));
    }
}
