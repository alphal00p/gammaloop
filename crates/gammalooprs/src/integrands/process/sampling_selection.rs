//! Resolution of user-facing graph sampling channel selections.
//!
//! Settings deliberately keep selectors as strings so that TOML and Python
//! have the same compact interface.  This module turns those strings into a
//! typed, deterministic selection after the graph name is known.  Numerical
//! channel construction is intentionally kept out of this layer.

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
    str::FromStr,
};

use crate::settings::runtime::{SamplingChannelDefinition, SamplingChannelSelection};
use color_eyre::eyre::{Result, eyre};
use serde::{Deserialize, Serialize};

use super::sampling_maps::combine_contracts;
use super::{
    ImplicitSurfaceRadialMap, SamplingChannelScore, SamplingMapAffine, SamplingMapComponent,
    SamplingMapContract, SamplingMapDefinition, SamplingMapEmbedding, SamplingMapEvaluation,
    SamplingMapKernel, SamplingPartition, SamplingPartitionMode, SamplingScoreFunction,
    SurfaceRadialMap,
};
use crate::momentum::sample::{LoopMomenta, MomentumSample};
use crate::settings::runtime::ParameterizationSettings;
use crate::settings::runtime::kinematic::Externals;
use crate::utils::{F, FloatLike};
use crate::{DependentMomentaConstructor, momentum::ThreeMomentum};

/// Built-in selectors understood by the channel catalogue.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SamplingChannelPreset {
    /// Every admissible loop-momentum basis channel.
    Lmb,
    /// The existing optimized LMB heuristic and its soft-coverage audit.
    OptimizedLmb,
    /// Surface-aware channels together with the ordinary coverage channels.
    Surfaces,
}

impl SamplingChannelPreset {
    pub const LMB: &'static str = "auto:lmb";
    pub const OPTIMIZED_LMB: &'static str = "auto:optimized_lmb";
    pub const SURFACES: &'static str = "auto:surfaces";

    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Lmb => Self::LMB,
            Self::OptimizedLmb => Self::OPTIMIZED_LMB,
            Self::Surfaces => Self::SURFACES,
        }
    }
}

/// One normalized selector in a graph's channel selection.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum SamplingChannelSelector {
    Preset(SamplingChannelPreset),
    Named(String),
}

impl SamplingChannelSelector {
    pub fn parse(source: &str) -> Result<Self, SamplingSelectionError> {
        let source = source.trim();
        if source.is_empty() {
            return Err(SamplingSelectionError::EmptySelector);
        }
        match source {
            SamplingChannelPreset::LMB => Ok(Self::Preset(SamplingChannelPreset::Lmb)),
            SamplingChannelPreset::OPTIMIZED_LMB => {
                Ok(Self::Preset(SamplingChannelPreset::OptimizedLmb))
            }
            SamplingChannelPreset::SURFACES => Ok(Self::Preset(SamplingChannelPreset::Surfaces)),
            value if value.starts_with("auto:") => {
                Err(SamplingSelectionError::UnknownPreset(value.to_owned()))
            }
            value => Ok(Self::Named(value.to_owned())),
        }
    }

    pub const fn preset(&self) -> Option<SamplingChannelPreset> {
        match self {
            Self::Preset(preset) => Some(*preset),
            Self::Named(_) => None,
        }
    }

    pub fn name(&self) -> Option<&str> {
        match self {
            Self::Named(name) => Some(name),
            Self::Preset(_) => None,
        }
    }
}

impl fmt::Display for SamplingChannelSelector {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Preset(preset) => formatter.write_str(preset.as_str()),
            Self::Named(name) => formatter.write_str(name),
        }
    }
}

impl FromStr for SamplingChannelSelector {
    type Err = SamplingSelectionError;

    fn from_str(source: &str) -> Result<Self, Self::Err> {
        Self::parse(source)
    }
}

/// A named channel together with its parsed Symbolica map definition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolvedNamedSamplingChannel {
    pub name: String,
    pub definition: SamplingChannelDefinition,
    pub map: SamplingMapDefinition,
}

/// The complete selection for one graph after applying the settings' fallback.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolvedSamplingChannelSelection {
    pub graph_name: String,
    pub selectors: Vec<SamplingChannelSelector>,
    pub named_channels: Vec<ResolvedNamedSamplingChannel>,
}

/// One entry in the graph-local sampling catalogue. The catalogue is the
/// migration target for grid construction and evaluation. The existing setup
/// only supplies generated entries until that runtime migration is complete;
/// it is not a second production channel universe.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SamplingCatalogueEntry {
    Lmb {
        basis_id: usize,
        edges: Vec<usize>,
        preset: SamplingChannelPreset,
    },
    /// An automatically enumerated E-surface candidate. Its geometry is
    /// prepared later from the current cut/orientation kinematics.
    Surface {
        edges: Vec<usize>,
        parent_lmb: Vec<usize>,
    },
    Named(ResolvedNamedSamplingChannel),
}

/// Graph-scoped, deterministically ordered channel catalogue.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SamplingChannelCatalogue {
    pub graph_name: String,
    pub selectors: Vec<SamplingChannelSelector>,
    pub entries: Vec<SamplingCatalogueEntry>,
}

/// Stable, side-effect-free view of the resolved catalogue for diagnostics.
///
/// This report is deliberately built from [`SamplingChannelCatalogue`] rather
/// than enumerating channels independently.  CLI and Python inspection can
/// therefore show exactly the catalogue that the production sampling path
/// will compile, without constructing kinematics or running a sample.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SamplingChannelInspection {
    pub graph_name: String,
    pub selectors: Vec<String>,
    pub entries: Vec<String>,
}

/// Kinematic data needed when compiling a graph-local surface channel.
///
/// The centre and threshold radius are intentionally supplied by the process
/// layer: they depend on the prepared external/cut kinematics and cannot be
/// inferred from a symbolic `surface(...)` expression alone.  `loop_edges`
/// are master-graph edge ids and define the local raw-frame embedding of the
/// compiled block.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingSurfaceGeometry {
    pub center: Vec<f64>,
    pub threshold_radius: Option<f64>,
    pub beta: f64,
    pub power: f64,
}

/// Input context for graph-independent compilation of catalogue entries.
///
/// A caller must identify the master graph explicitly.  This prevents a
/// surface block from being mistaken for a coordinate block in a different
/// graph when channel definitions are shared across a graph group.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingChannelCompileContext {
    pub master_graph: String,
    /// Complete ordered parent LMB in the master graph frame. Every named
    /// channel is resolved against this exact list before compilation.
    pub parent_lmb: Vec<usize>,
    pub parameterization_settings: ParameterizationSettings,
    pub e_cm: f64,
    pub n_loop_momenta: usize,
    pub surfaces: BTreeMap<Vec<usize>, SamplingSurfaceGeometry>,
    /// Host cut identity for metadata whose `on_cut` list is explicit.
    pub cut_id: Option<usize>,
    pub orientation: Option<usize>,
    pub side: Option<super::SamplingCutSide>,
    /// Exact affine maps routing a generated LMB into the master frame. The
    /// key is the generated catalogue basis id; absent entries are diagnosed
    /// when a non-master LMB is selected.
    pub lmb_frame_maps: BTreeMap<usize, SamplingMapAffine>,
    /// Explicit named `lmb(...)` definitions use their edge list rather than
    /// a generated basis id, so their prepared routing is keyed by that same
    /// master-graph edge list.
    pub lmb_frame_maps_by_edges: BTreeMap<Vec<usize>, SamplingMapAffine>,
}

impl SamplingChannelCompileContext {
    pub fn new(
        master_graph: impl Into<String>,
        parent_lmb: Vec<usize>,
        parameterization_settings: ParameterizationSettings,
        e_cm: f64,
        n_loop_momenta: usize,
    ) -> Self {
        Self {
            master_graph: master_graph.into(),
            parent_lmb,
            parameterization_settings,
            e_cm,
            n_loop_momenta,
            surfaces: BTreeMap::new(),
            lmb_frame_maps: BTreeMap::new(),
            lmb_frame_maps_by_edges: BTreeMap::new(),
            cut_id: None,
            orientation: None,
            side: None,
        }
    }
}

/// The map kernels currently compilable without a graph-specific implicit
/// solver.  Composite/cut maps stay in the typed catalogue until their
/// prepared kinematic context is supplied by the process layer.
#[derive(Clone, Debug)]
pub enum CompiledSamplingMap {
    Lmb(SamplingMapKernel),
    AffineLmb {
        lmb: SamplingMapKernel,
        frame: SamplingMapAffine,
    },
    Surface(SurfaceRadialMap),
    ImplicitSurface(ImplicitSurfaceRadialMap),
    Embedded(SamplingMapEmbedding),
}

impl CompiledSamplingMap {
    pub fn contract(&self) -> SamplingMapContract {
        match self {
            Self::Lmb(map) => map.contract(),
            Self::AffineLmb { lmb, frame } => combine_contracts(lmb.contract(), frame.contract()),
            Self::Surface(map) => map.contract(),
            Self::ImplicitSurface(map) => map.contract(),
            Self::Embedded(map) => map.contract(),
        }
    }

    pub fn dimensions(&self) -> usize {
        match self {
            Self::Lmb(map) => map.dimensions(),
            Self::AffineLmb { lmb, .. } => lmb.dimensions(),
            Self::Surface(map) => map.dimension(),
            Self::ImplicitSurface(map) => map.dimension(),
            Self::Embedded(map) => map.dimensions(),
        }
    }

    pub fn as_component(&self) -> &dyn SamplingMapComponent {
        self
    }

    pub fn forward(&self, coordinates: &[f64]) -> Result<SamplingMapEvaluation> {
        <Self as SamplingMapComponent>::forward(self, coordinates, &[])
    }

    pub fn inverse(&self, point: &[f64]) -> Result<SamplingMapEvaluation> {
        <Self as SamplingMapComponent>::inverse(self, point, &[])
    }
}

impl SamplingMapComponent for CompiledSamplingMap {
    fn dimensions(&self) -> usize {
        self.dimensions()
    }

    fn output_dimensions(&self) -> usize {
        self.dimensions()
    }

    fn contract(&self) -> SamplingMapContract {
        self.contract()
    }

    fn name(&self) -> &'static str {
        match self {
            Self::Lmb(_) => "lmb",
            Self::AffineLmb { .. } => "affine_lmb",
            Self::Surface(_) => "surface",
            Self::ImplicitSurface(_) => "implicit_surface",
            Self::Embedded(_) => "embedded",
        }
    }

    fn forward(&self, coordinates: &[f64], context: &[f64]) -> Result<SamplingMapEvaluation> {
        match self {
            Self::Lmb(map) => SamplingMapComponent::forward(map, coordinates, context),
            Self::Surface(map) => SamplingMapComponent::forward(map, coordinates, context),
            Self::ImplicitSurface(map) => SamplingMapComponent::forward(map, coordinates, context),
            Self::Embedded(map) => SamplingMapComponent::forward(map, coordinates, context),
            Self::AffineLmb { lmb, frame } => {
                let lmb_evaluation = SamplingMapComponent::forward(lmb, coordinates, context)?;
                let frame_evaluation =
                    SamplingMapComponent::forward(frame, &lmb_evaluation.point, context)?;
                Ok(SamplingMapEvaluation {
                    coordinates: lmb_evaluation.coordinates,
                    point: frame_evaluation.point,
                    jacobian: lmb_evaluation.jacobian * frame_evaluation.jacobian,
                    inverse_jacobian: lmb_evaluation.inverse_jacobian
                        * frame_evaluation.inverse_jacobian,
                    residual: lmb_evaluation.residual.max(frame_evaluation.residual),
                    support: combine_contracts(lmb.contract(), frame.contract()).support,
                    diagnostics: lmb_evaluation
                        .diagnostics
                        .into_iter()
                        .chain(frame_evaluation.diagnostics)
                        .collect(),
                })
            }
        }
    }

    fn inverse(&self, point: &[f64], context: &[f64]) -> Result<SamplingMapEvaluation> {
        match self {
            Self::Lmb(map) => SamplingMapComponent::inverse(map, point, context),
            Self::Surface(map) => SamplingMapComponent::inverse(map, point, context),
            Self::ImplicitSurface(map) => SamplingMapComponent::inverse(map, point, context),
            Self::Embedded(map) => SamplingMapComponent::inverse(map, point, context),
            Self::AffineLmb { lmb, frame } => {
                let frame_evaluation = SamplingMapComponent::inverse(frame, point, context)?;
                let lmb_evaluation =
                    SamplingMapComponent::inverse(lmb, &frame_evaluation.point, context)?;
                Ok(SamplingMapEvaluation {
                    coordinates: lmb_evaluation.coordinates,
                    point: point.to_vec(),
                    jacobian: lmb_evaluation.jacobian * frame_evaluation.jacobian,
                    inverse_jacobian: lmb_evaluation.inverse_jacobian
                        * frame_evaluation.inverse_jacobian,
                    residual: lmb_evaluation.residual.max(frame_evaluation.residual),
                    support: combine_contracts(lmb.contract(), frame.contract()).support,
                    diagnostics: lmb_evaluation
                        .diagnostics
                        .into_iter()
                        .chain(frame_evaluation.diagnostics)
                        .collect(),
                })
            }
        }
    }
}

/// A compiled graph channel and the master-graph raw-frame block it occupies.
#[derive(Clone, Debug)]
pub struct CompiledSamplingChannel {
    pub name: String,
    pub master_graph: String,
    pub basis_id: Option<usize>,
    pub definition: SamplingMapDefinition,
    /// Ordered master-graph edge ids for this raw coordinate block.
    pub embedded_edges: Vec<usize>,
    pub map: CompiledSamplingMap,
}

/// Canonical identifier for a compiled sampling channel.
///
/// The same id is used by the graph evaluator and by the sampling catalogue;
/// no second legacy channel-index domain is maintained.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Ord, PartialOrd, Serialize, Deserialize)]
pub struct SamplingChannelId(pub usize);

impl SamplingChannelId {
    pub const fn index(self) -> usize {
        self.0
    }
}

impl From<usize> for SamplingChannelId {
    fn from(index: usize) -> Self {
        Self(index)
    }
}

impl From<SamplingChannelId> for usize {
    fn from(index: SamplingChannelId) -> Self {
        index.0
    }
}

/// A bounded bridge from compiled graph channels to the existing graph
/// evaluator.  The bridge deliberately deals in the complete master raw
/// frame: it never silently embeds a lower-dimensional surface block or
/// reinterprets a channel in another graph's loop basis.
#[derive(Clone, Debug)]
pub struct SamplingChannelBridge {
    channels: Vec<CompiledSamplingChannel>,
    scores: Vec<SamplingChannelScore>,
    dimensions: usize,
}

/// One push-forward/inverse result together with the common raw-frame
/// multichannel partition.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingChannelBridgeEvaluation {
    pub channel_id: SamplingChannelId,
    pub channel_name: String,
    /// The complete master raw coordinate frame passed to graph evaluation.
    pub raw_coordinates: Vec<f64>,
    pub map: SamplingMapEvaluation,
    pub partition: SamplingPartition,
}

/// Runtime context needed to turn one bridge point into the sample consumed by
/// the graph evaluator.  The bridge itself is graph-frame aware, while this
/// small context supplies the cache and external-momentum ownership that are
/// deliberately process-local.
#[derive(Clone, Copy, Debug)]
pub struct SamplingMomentumSampleContext<'a> {
    pub loop_mom_cache_id: usize,
    pub external_moms: &'a Externals,
    pub external_mom_cache_id: usize,
    pub dependent_momenta_constructor: DependentMomentaConstructor<'a>,
    pub orientation: Option<usize>,
}

impl SamplingChannelBridgeEvaluation {
    /// Materialize this exact full-frame bridge point as a graph-evaluator
    /// momentum sample. The production sampler uses the same conversion after
    /// graph-level channel and cut preparation.
    pub fn to_momentum_sample<T: FloatLike>(
        &self,
        context: SamplingMomentumSampleContext<'_>,
    ) -> Result<MomentumSample<T>> {
        if self.raw_coordinates.len() != self.map.point.len() {
            return Err(eyre!(
                "sampling bridge point has {} coordinates but its raw frame has {}",
                self.map.point.len(),
                self.raw_coordinates.len()
            ));
        }
        if self.raw_coordinates.is_empty() || self.raw_coordinates.len() % 3 != 0 {
            return Err(eyre!(
                "sampling bridge raw frame has {} coordinates; expected a non-empty multiple of 3",
                self.raw_coordinates.len()
            ));
        }
        if self.raw_coordinates.iter().any(|value| !value.is_finite()) {
            return Err(eyre!(
                "sampling bridge raw frame contains a non-finite momentum component"
            ));
        }
        if !self.map.jacobian.is_finite() || self.map.jacobian <= 0.0 {
            return Err(eyre!(
                "sampling bridge map Jacobian must be finite and positive, got {}",
                self.map.jacobian
            ));
        }

        let loop_momenta =
            LoopMomenta::from_iter(self.raw_coordinates.chunks_exact(3).map(|components| {
                ThreeMomentum::new(
                    F::<T>::from_f64(components[0]),
                    F::<T>::from_f64(components[1]),
                    F::<T>::from_f64(components[2]),
                )
            }));
        MomentumSample::new(
            loop_momenta,
            context.loop_mom_cache_id,
            context.external_moms,
            context.external_mom_cache_id,
            F::<T>::from_f64(self.map.jacobian),
            context.dependent_momenta_constructor,
            context.orientation,
        )
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SamplingChannelBridgeError {
    Empty,
    DimensionMismatch {
        channel: String,
        dimensions: usize,
        expected: usize,
    },
    IncompatibleFrames {
        channel: String,
        edges: Vec<usize>,
        expected: Vec<usize>,
    },
    MasterGraphMismatch {
        channel: String,
        graph: String,
        expected: String,
    },
    PartialSupport {
        channel: String,
        support: super::SamplingSupport,
    },
    InexactJacobian {
        channel: String,
        jacobian: super::SamplingJacobian,
    },
    InvalidChannel {
        channel: String,
        error: String,
    },
    UnknownChannel {
        channel: SamplingChannelId,
    },
}

impl fmt::Display for SamplingChannelBridgeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => {
                formatter.write_str("sampling channel bridge requires at least one channel")
            }
            Self::DimensionMismatch {
                channel,
                dimensions,
                expected,
            } => write!(
                formatter,
                "sampling channel `{channel}` has {dimensions} raw dimensions; the bridge requires {expected}"
            ),
            Self::IncompatibleFrames {
                channel,
                edges,
                expected,
            } => write!(
                formatter,
                "sampling channel `{channel}` uses master-frame edges {edges:?}, incompatible with bridge frame {expected:?}"
            ),
            Self::MasterGraphMismatch {
                channel,
                graph,
                expected,
            } => write!(
                formatter,
                "sampling channel `{channel}` belongs to master graph `{graph}`, but the bridge frame is `{expected}`"
            ),
            Self::PartialSupport { channel, support } => write!(
                formatter,
                "sampling channel `{channel}` has {support:?} support; partial/conditional surface maps cannot be passed to the graph evaluator"
            ),
            Self::InexactJacobian { channel, jacobian } => write!(
                formatter,
                "sampling channel `{channel}` supplies {jacobian:?}; the evaluator bridge requires an exact map Jacobian"
            ),
            Self::InvalidChannel { channel, error } => {
                write!(
                    formatter,
                    "sampling channel `{channel}` is invalid: {error}"
                )
            }
            Self::UnknownChannel { channel } => {
                write!(
                    formatter,
                    "sampling channel index {} is out of range",
                    channel.0
                )
            }
        }
    }
}

impl std::error::Error for SamplingChannelBridgeError {}

impl SamplingChannelBridge {
    /// Construct an exact map-density bridge.  Every map must cover the full
    /// raw frame and have an exact forward determinant; this rejects partial
    /// surface maps before they can reach graph evaluation.
    pub fn new(
        channels: Vec<CompiledSamplingChannel>,
    ) -> std::result::Result<Self, SamplingChannelBridgeError> {
        let Some(first) = channels.first() else {
            return Err(SamplingChannelBridgeError::Empty);
        };
        let dimensions = first.dimensions();
        let frame = &first.embedded_edges;
        let master_graph = &first.master_graph;
        let mut scores = Vec::with_capacity(channels.len());
        for channel in &channels {
            if channel.master_graph != *master_graph {
                return Err(SamplingChannelBridgeError::MasterGraphMismatch {
                    channel: channel.name.clone(),
                    graph: channel.master_graph.clone(),
                    expected: master_graph.clone(),
                });
            }
            if channel.embedded_edges != *frame {
                return Err(SamplingChannelBridgeError::IncompatibleFrames {
                    channel: channel.name.clone(),
                    edges: channel.embedded_edges.clone(),
                    expected: frame.clone(),
                });
            }
            if channel.dimensions() != dimensions {
                return Err(SamplingChannelBridgeError::DimensionMismatch {
                    channel: channel.name.clone(),
                    dimensions: channel.dimensions(),
                    expected: dimensions,
                });
            }
            let contract = channel.contract();
            if contract.support != super::SamplingSupport::Full {
                return Err(SamplingChannelBridgeError::PartialSupport {
                    channel: channel.name.clone(),
                    support: contract.support,
                });
            }
            if matches!(contract.jacobian, super::SamplingJacobian::ProxyOnly) {
                return Err(SamplingChannelBridgeError::InexactJacobian {
                    channel: channel.name.clone(),
                    jacobian: contract.jacobian,
                });
            }
            let map = channel.map.clone();
            let name = channel.name.clone();
            scores.push(SamplingChannelScore::map_density(
                name,
                SamplingScoreFunction::from_positive_function(move |raw| {
                    let evaluation = map.inverse(raw)?;
                    let density = evaluation.inverse_jacobian.abs();
                    if !density.is_finite() || density <= 0.0 {
                        return Err(color_eyre::eyre::eyre!(
                            "inverse map density is not finite and positive: {density}"
                        ));
                    }
                    Ok(Some(density))
                }),
            ));
        }
        Ok(Self {
            channels,
            scores,
            dimensions,
        })
    }

    pub fn channels(&self) -> &[CompiledSamplingChannel] {
        &self.channels
    }

    pub fn dimensions(&self) -> usize {
        self.dimensions
    }

    pub fn partition(&self, raw_coordinates: &[f64]) -> Result<SamplingPartition> {
        if raw_coordinates.len() != self.dimensions {
            return Err(color_eyre::eyre::eyre!(
                "raw sampling frame has {}, expected {} dimensions",
                raw_coordinates.len(),
                self.dimensions
            ));
        }
        SamplingPartition::new(
            SamplingPartitionMode::MapDensity,
            &self.scores,
            raw_coordinates,
        )
    }

    pub fn forward(
        &self,
        channel_id: SamplingChannelId,
        coordinates: &[f64],
    ) -> Result<SamplingChannelBridgeEvaluation> {
        let channel_index = channel_id.0;
        let channel =
            self.channels
                .get(channel_index)
                .ok_or(SamplingChannelBridgeError::UnknownChannel {
                    channel: channel_id,
                })?;
        let map = channel.forward(coordinates)?;
        if map.point.len() != self.dimensions {
            return Err(SamplingChannelBridgeError::DimensionMismatch {
                channel: channel.name.clone(),
                dimensions: map.point.len(),
                expected: self.dimensions,
            }
            .into());
        }
        let partition = self.partition(&map.point)?;
        Ok(SamplingChannelBridgeEvaluation {
            channel_id,
            channel_name: channel.name.clone(),
            raw_coordinates: map.point.clone(),
            map,
            partition,
        })
    }

    pub fn inverse(
        &self,
        channel_id: SamplingChannelId,
        raw_coordinates: &[f64],
    ) -> Result<SamplingChannelBridgeEvaluation> {
        let channel_index = channel_id.0;
        let channel =
            self.channels
                .get(channel_index)
                .ok_or(SamplingChannelBridgeError::UnknownChannel {
                    channel: channel_id,
                })?;
        let map = channel.inverse(raw_coordinates)?;
        let partition = self.partition(raw_coordinates)?;
        Ok(SamplingChannelBridgeEvaluation {
            channel_id,
            channel_name: channel.name.clone(),
            raw_coordinates: raw_coordinates.to_vec(),
            map,
            partition,
        })
    }
}

impl CompiledSamplingChannel {
    pub fn contract(&self) -> SamplingMapContract {
        self.map.contract()
    }

    pub fn dimensions(&self) -> usize {
        self.map.dimensions()
    }

    pub fn forward(&self, coordinates: &[f64]) -> Result<SamplingMapEvaluation> {
        self.map.forward(coordinates)
    }

    pub fn inverse(&self, point: &[f64]) -> Result<SamplingMapEvaluation> {
        self.map.inverse(point)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SamplingChannelCompileError {
    EmptyMasterGraph,
    UnsupportedMap {
        channel: String,
        map: String,
    },
    MissingSurfaceGeometry {
        channel: String,
        edges: Vec<usize>,
    },
    MissingLmbFrameMap {
        channel: String,
        basis_id: Option<usize>,
        edges: Vec<usize>,
        parent_lmb: Vec<usize>,
    },
    InvalidLmbFrameMap {
        channel: String,
        basis_id: usize,
        error: String,
    },
    InvalidChannel {
        channel: String,
        error: String,
    },
}

impl fmt::Display for SamplingChannelCompileError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyMasterGraph => {
                formatter.write_str("sampling channel compilation requires a master graph")
            }
            Self::UnsupportedMap { channel, map } => write!(
                formatter,
                "sampling channel `{channel}` uses map `{map}`, which needs prepared graph context before it can be compiled"
            ),
            Self::MissingSurfaceGeometry { channel, edges } => write!(
                formatter,
                "sampling channel `{channel}` has no prepared surface geometry for master-graph edges {edges:?}"
            ),
            Self::MissingLmbFrameMap {
                channel,
                basis_id,
                edges,
                parent_lmb,
            } => write!(
                formatter,
                "sampling channel `{channel}` uses non-master LMB edges {edges:?} (parent LMB {parent_lmb:?}) but no affine frame map was supplied for basis {basis_id:?}"
            ),
            Self::InvalidLmbFrameMap {
                channel,
                basis_id,
                error,
            } => write!(
                formatter,
                "sampling channel `{channel}` has invalid affine frame map for basis {basis_id}: {error}"
            ),
            Self::InvalidChannel { channel, error } => write!(
                formatter,
                "sampling channel `{channel}` could not be compiled: {error}"
            ),
        }
    }
}

impl std::error::Error for SamplingChannelCompileError {}

/// Compile one radial surface in its native subspace and embed it together
/// with the ordinary map on the complementary parent-LMB edges.  The output
/// permutation is explicit, so sharing a surface across graph channels never
/// relies on an implicit edge ordering or a silent local-frame assumption.
fn compile_surface_map(
    channel: &str,
    edges: &[usize],
    context: &SamplingChannelCompileContext,
) -> Result<CompiledSamplingMap, SamplingChannelCompileError> {
    if edges.is_empty() {
        return Err(SamplingChannelCompileError::InvalidChannel {
            channel: channel.to_owned(),
            error: "surface subspace must contain at least one edge".to_owned(),
        });
    }
    let parent_positions = context
        .parent_lmb
        .iter()
        .enumerate()
        .map(|(position, edge)| (*edge, position))
        .collect::<BTreeMap<_, _>>();
    if edges
        .iter()
        .any(|edge| !parent_positions.contains_key(edge))
    {
        return Err(SamplingChannelCompileError::InvalidChannel {
            channel: channel.to_owned(),
            error: format!(
                "surface edges {edges:?} are not a subspace of parent LMB {:?}",
                context.parent_lmb
            ),
        });
    }
    let Some(geometry) = context.surfaces.get(edges) else {
        return Err(SamplingChannelCompileError::MissingSurfaceGeometry {
            channel: channel.to_owned(),
            edges: edges.to_vec(),
        });
    };
    let expected_dimension = 3 * edges.len();
    if geometry.center.len() != expected_dimension {
        return Err(SamplingChannelCompileError::InvalidChannel {
            channel: channel.to_owned(),
            error: format!(
                "surface centre has dimension {}, expected {} for edges {edges:?}",
                geometry.center.len(),
                expected_dimension
            ),
        });
    }
    let surface = SurfaceRadialMap::new(
        expected_dimension,
        geometry.center.clone(),
        geometry.threshold_radius,
        geometry.beta,
        geometry.power,
    )
    .map_err(|error| SamplingChannelCompileError::InvalidChannel {
        channel: channel.to_owned(),
        error: error.to_string(),
    })?;
    if edges.len() == context.n_loop_momenta {
        if edges == context.parent_lmb {
            return Ok(CompiledSamplingMap::Surface(surface));
        }
        let output_indices = edges
            .iter()
            .map(|edge| parent_positions[edge])
            .flat_map(|position| [3 * position, 3 * position + 1, 3 * position + 2])
            .collect();
        return SamplingMapEmbedding::product(vec![Box::new(surface)], output_indices)
            .map(CompiledSamplingMap::Embedded)
            .map_err(|error| SamplingChannelCompileError::InvalidChannel {
                channel: channel.to_owned(),
                error: error.to_string(),
            });
    }

    let complement = context
        .parent_lmb
        .iter()
        .copied()
        .filter(|edge| !edges.contains(edge))
        .collect::<Vec<_>>();
    if complement.is_empty() {
        return Err(SamplingChannelCompileError::InvalidChannel {
            channel: channel.to_owned(),
            error: format!("surface edges {edges:?} do not leave a complementary parent-LMB block"),
        });
    }
    let complement_map = SamplingMapKernel::new(
        SamplingMapDefinition::Lmb(complement.clone()),
        context.parameterization_settings.clone(),
        context.e_cm,
        complement.len(),
    )
    .map_err(|error| SamplingChannelCompileError::InvalidChannel {
        channel: channel.to_owned(),
        error: format!("complement LMB {complement:?} could not be compiled: {error}"),
    })?;

    // Product order is surface block followed by complement block.  Convert
    // each block's edge ordering to the complete parent-LMB component frame.
    let mut output_indices = Vec::with_capacity(3 * context.n_loop_momenta);
    for edge in edges.iter().chain(complement.iter()) {
        let position = parent_positions[edge];
        output_indices.extend([3 * position, 3 * position + 1, 3 * position + 2]);
    }
    SamplingMapEmbedding::product(
        vec![Box::new(surface), Box::new(complement_map)],
        output_indices,
    )
    .map(CompiledSamplingMap::Embedded)
    .map_err(|error| SamplingChannelCompileError::InvalidChannel {
        channel: channel.to_owned(),
        error: error.to_string(),
    })
}

fn compile_lmb_map(
    channel: &str,
    basis_id: Option<usize>,
    edges: &[usize],
    context: &SamplingChannelCompileContext,
) -> Result<CompiledSamplingMap, SamplingChannelCompileError> {
    let definition = SamplingMapDefinition::Lmb(edges.to_vec());
    let lmb = SamplingMapKernel::new(
        definition,
        context.parameterization_settings.clone(),
        context.e_cm,
        context.n_loop_momenta,
    )
    .map_err(|error| SamplingChannelCompileError::InvalidChannel {
        channel: channel.to_owned(),
        error: error.to_string(),
    })?;
    if edges == context.parent_lmb {
        return Ok(CompiledSamplingMap::Lmb(lmb));
    }
    let Some(basis_id) = basis_id else {
        return Err(SamplingChannelCompileError::MissingLmbFrameMap {
            channel: channel.to_owned(),
            basis_id: None,
            edges: edges.to_vec(),
            parent_lmb: context.parent_lmb.clone(),
        });
    };
    let Some(frame) = context.lmb_frame_maps.get(&basis_id) else {
        return Err(SamplingChannelCompileError::MissingLmbFrameMap {
            channel: channel.to_owned(),
            basis_id: Some(basis_id),
            edges: edges.to_vec(),
            parent_lmb: context.parent_lmb.clone(),
        });
    };
    let expected_dimension = 3 * context.n_loop_momenta;
    if frame.dimension() != expected_dimension {
        return Err(SamplingChannelCompileError::InvalidLmbFrameMap {
            channel: channel.to_owned(),
            basis_id,
            error: format!(
                "affine frame dimension {} does not match parent raw dimension {expected_dimension}",
                frame.dimension()
            ),
        });
    }
    Ok(CompiledSamplingMap::AffineLmb {
        lmb,
        frame: frame.clone(),
    })
}

impl SamplingChannelCatalogue {
    pub fn lmb_entries(&self) -> impl Iterator<Item = (usize, &[usize])> {
        self.entries.iter().filter_map(|entry| match entry {
            SamplingCatalogueEntry::Lmb {
                basis_id, edges, ..
            } => Some((*basis_id, edges.as_slice())),
            SamplingCatalogueEntry::Named(_) => None,
            SamplingCatalogueEntry::Surface { .. } => None,
        })
    }

    pub fn named_entries(&self) -> impl Iterator<Item = &ResolvedNamedSamplingChannel> {
        self.entries.iter().filter_map(|entry| match entry {
            SamplingCatalogueEntry::Named(channel) => Some(channel),
            SamplingCatalogueEntry::Lmb { .. } => None,
            SamplingCatalogueEntry::Surface { .. } => None,
        })
    }

    /// Compile the ordinary LMB and radial surface entries in this catalogue.
    ///
    /// Surface geometry is looked up by its canonical master-graph edge list.
    /// Maps such as `cut`, `intersect` and `then` are deliberately rejected
    /// here until the process has prepared their conditional kinematics; this
    /// avoids creating a numerically plausible map with an incorrect frame.
    pub fn compile(
        &self,
        context: &SamplingChannelCompileContext,
    ) -> Result<Vec<CompiledSamplingChannel>, SamplingChannelCompileError> {
        if context.master_graph.trim().is_empty() {
            return Err(SamplingChannelCompileError::EmptyMasterGraph);
        }
        let mut parent_edges = BTreeSet::new();
        if context.parent_lmb.len() != context.n_loop_momenta
            || context
                .parent_lmb
                .iter()
                .any(|edge| !parent_edges.insert(*edge))
        {
            return Err(SamplingChannelCompileError::InvalidChannel {
                channel: context.master_graph.clone(),
                error: format!(
                    "parent LMB {:?} must contain one unique edge per loop momentum (expected {})",
                    context.parent_lmb, context.n_loop_momenta
                ),
            });
        }
        let mut compiled = Vec::with_capacity(self.entries.len());
        for entry in &self.entries {
            let (name, basis_id, definition, map) = match entry {
                SamplingCatalogueEntry::Lmb {
                    basis_id, edges, ..
                } => {
                    let definition = SamplingMapDefinition::Lmb(edges.clone());
                    let map = compile_lmb_map(
                        &format!("lmb[{basis_id}]"),
                        Some(*basis_id),
                        edges,
                        context,
                    )?;
                    (format!("lmb[{basis_id}]"), Some(*basis_id), definition, map)
                }
                SamplingCatalogueEntry::Surface { edges, parent_lmb } => {
                    if parent_lmb != &context.parent_lmb {
                        return Err(SamplingChannelCompileError::InvalidChannel {
                            channel: format!("surface:{edges:?}"),
                            error: format!(
                                "parent LMB {:?} does not match master context {:?}",
                                parent_lmb, context.parent_lmb
                            ),
                        });
                    }
                    let definition = SamplingMapDefinition::Surface(edges.clone());
                    let map = compile_surface_map(&format!("surface:{edges:?}"), edges, context)?;
                    (format!("surface:{edges:?}"), None, definition, map)
                }
                SamplingCatalogueEntry::Named(channel) => {
                    if !channel.definition.on_cut.is_empty()
                        && !channel
                            .definition
                            .on_cut
                            .contains(&context.cut_id.ok_or_else(|| {
                                SamplingChannelCompileError::InvalidChannel {
                                    channel: channel.name.clone(),
                                    error: format!(
                                        "channel declares on_cut {:?}, but no prepared host cut identity was supplied",
                                        channel.definition.on_cut
                                    ),
                                }
                            })?)
                    {
                        return Err(SamplingChannelCompileError::InvalidChannel {
                            channel: channel.name.clone(),
                            error: format!(
                                "channel declares on_cut {:?}, incompatible with prepared cut {:?}",
                                channel.definition.on_cut, context.cut_id
                            ),
                        });
                    }
                    if channel.definition.parent_lmb != context.parent_lmb {
                        return Err(SamplingChannelCompileError::InvalidChannel {
                            channel: channel.name.clone(),
                            error: format!(
                                "parent LMB {:?} does not match master context {:?}",
                                channel.definition.parent_lmb, context.parent_lmb
                            ),
                        });
                    }
                    let definition = channel.map.clone();
                    match &definition {
                        SamplingMapDefinition::Left(_) | SamplingMapDefinition::Right(_) => {
                            let required_side =
                                if matches!(&definition, SamplingMapDefinition::Left(_)) {
                                    super::SamplingCutSide::Left
                                } else {
                                    super::SamplingCutSide::Right
                                };
                            if context.side != Some(required_side) {
                                return Err(SamplingChannelCompileError::InvalidChannel {
                                    channel: channel.name.clone(),
                                    error: format!(
                                        "map `{definition:?}` requires prepared side {required_side:?}, but context supplies {:?}",
                                        context.side
                                    ),
                                });
                            }
                        }
                        _ => {}
                    }
                    let map = match &definition {
                        SamplingMapDefinition::Lmb(edges) => {
                            let map = SamplingMapKernel::new(
                                definition.clone(),
                                context.parameterization_settings.clone(),
                                context.e_cm,
                                context.n_loop_momenta,
                            )
                            .map_err(|error| {
                                SamplingChannelCompileError::InvalidChannel {
                                    channel: channel.name.clone(),
                                    error: error.to_string(),
                                }
                            })?;
                            if edges == &context.parent_lmb {
                                CompiledSamplingMap::Lmb(map)
                            } else {
                                let frame = context.lmb_frame_maps_by_edges.get(edges).ok_or_else(
                                    || SamplingChannelCompileError::MissingLmbFrameMap {
                                        channel: channel.name.clone(),
                                        basis_id: None,
                                        edges: edges.clone(),
                                        parent_lmb: context.parent_lmb.clone(),
                                    },
                                )?;
                                if frame.dimension() != 3 * context.n_loop_momenta {
                                    return Err(SamplingChannelCompileError::InvalidLmbFrameMap {
                                        channel: channel.name.clone(),
                                        basis_id: usize::MAX,
                                        error: format!(
                                            "affine frame routing has dimension {}, expected {}",
                                            frame.dimension(),
                                            3 * context.n_loop_momenta
                                        ),
                                    });
                                }
                                CompiledSamplingMap::AffineLmb {
                                    lmb: map,
                                    frame: frame.clone(),
                                }
                            }
                        }
                        SamplingMapDefinition::Surface(edges) => {
                            compile_surface_map(&channel.name, edges, context)?
                        }
                        unsupported => {
                            return Err(SamplingChannelCompileError::UnsupportedMap {
                                channel: channel.name.clone(),
                                map: format!("{unsupported:?}"),
                            });
                        }
                    };
                    (channel.name.clone(), None, definition, map)
                }
            };
            let embedded_edges = match &map {
                CompiledSamplingMap::Embedded(_) | CompiledSamplingMap::AffineLmb { .. } => {
                    context.parent_lmb.clone()
                }
                _ => match &definition {
                    SamplingMapDefinition::Lmb(edges)
                    | SamplingMapDefinition::Surface(edges)
                    | SamplingMapDefinition::Complement(edges)
                    | SamplingMapDefinition::Cut(edges) => edges.clone(),
                    _ => Vec::new(),
                },
            };
            compiled.push(CompiledSamplingChannel {
                name,
                master_graph: context.master_graph.clone(),
                basis_id,
                definition,
                embedded_edges,
                map,
            });
        }
        Ok(compiled)
    }

    /// Human-readable inspection rows, stable across runs and suitable for
    /// CLI/API diagnostics.
    pub fn inspection_rows(&self) -> Vec<String> {
        self.entries
            .iter()
            .enumerate()
            .map(|(index, entry)| match entry {
                SamplingCatalogueEntry::Lmb {
                    basis_id,
                    edges,
                    preset,
                } => format!(
                    "{index}: lmb basis={basis_id} edges={edges:?} source={}",
                    preset.as_str()
                ),
                SamplingCatalogueEntry::Surface { edges, parent_lmb } => {
                    format!("{index}: surface edges={edges:?} parent_lmb={parent_lmb:?}")
                }
                SamplingCatalogueEntry::Named(channel) => format!(
                    "{index}: {} around={} parent_lmb={:?} on_cut={:?}",
                    channel.name,
                    channel.definition.around,
                    channel.definition.parent_lmb,
                    channel.definition.on_cut
                ),
            })
            .collect()
    }

    /// Return the canonical resolved selectors and catalogue rows for
    /// read-only CLI/API inspection.
    pub fn inspection(&self) -> SamplingChannelInspection {
        SamplingChannelInspection {
            graph_name: self.graph_name.clone(),
            selectors: self.selectors.iter().map(ToString::to_string).collect(),
            entries: self.inspection_rows(),
        }
    }
}

/// Expand a resolved graph selection against the generated LMB catalogue.
/// Entries are deduplicated by their resolved identity while preserving the
/// first selector's order.  A surface preset includes the optimized LMB
/// fallback so that a surface-aware selection remains defined where a surface
/// is absent; numerical use of these entries is implemented by the process
/// sampler in a later phase.
pub fn build_sampling_channel_catalogue(
    resolved: &ResolvedSamplingChannelSelection,
    all_lmbs: &[(usize, Vec<usize>)],
    optimized_lmbs: &[usize],
) -> SamplingChannelCatalogue {
    build_sampling_channel_catalogue_with_surfaces(resolved, all_lmbs, optimized_lmbs, &[], &[])
}

/// Expand a selection and append automatically enumerated surface candidates.
/// `surface_edges` and `parent_lmb` are already resolved in the master graph
/// frame; no graph mapping is inferred here.
pub fn build_sampling_channel_catalogue_with_surfaces(
    resolved: &ResolvedSamplingChannelSelection,
    all_lmbs: &[(usize, Vec<usize>)],
    optimized_lmbs: &[usize],
    surface_edges: &[Vec<usize>],
    parent_lmb: &[usize],
) -> SamplingChannelCatalogue {
    let mut entries = Vec::new();
    for selector in &resolved.selectors {
        let Some(preset) = selector.preset() else {
            if let SamplingChannelSelector::Named(name) = selector {
                if let Some(channel) = resolved.named(name) {
                    let entry = SamplingCatalogueEntry::Named(channel.clone());
                    if !entries.contains(&entry) {
                        entries.push(entry);
                    }
                }
            }
            continue;
        };
        let basis_ids: Vec<usize> = match preset {
            SamplingChannelPreset::Lmb => all_lmbs.iter().map(|(id, _)| *id).collect(),
            SamplingChannelPreset::OptimizedLmb | SamplingChannelPreset::Surfaces => {
                optimized_lmbs.to_vec()
            }
        };
        for basis_id in basis_ids {
            let Some((_, edges)) = all_lmbs.iter().find(|(id, _)| *id == basis_id) else {
                continue;
            };
            let entry = SamplingCatalogueEntry::Lmb {
                basis_id,
                edges: edges.clone(),
                preset,
            };
            let duplicate = entries.iter().any(|existing| {
                matches!(
                    existing,
                    SamplingCatalogueEntry::Lmb {
                        basis_id: existing_basis,
                        edges: existing_edges,
                        ..
                    } if *existing_basis == basis_id && existing_edges == edges
                )
            });
            if !duplicate {
                entries.push(entry);
            }
        }
        if preset == SamplingChannelPreset::Surfaces {
            for edges in surface_edges {
                if !edges.is_empty()
                    && !entries.iter().any(|entry| {
                        matches!(entry, SamplingCatalogueEntry::Surface {
                            edges: existing, parent_lmb: parent
                        } if existing == edges && parent == parent_lmb)
                    })
                {
                    entries.push(SamplingCatalogueEntry::Surface {
                        edges: edges.clone(),
                        parent_lmb: parent_lmb.to_vec(),
                    });
                }
            }
        }
    }
    SamplingChannelCatalogue {
        graph_name: resolved.graph_name.clone(),
        selectors: resolved.selectors.clone(),
        entries,
    }
}

impl ResolvedSamplingChannelSelection {
    pub fn has_preset(&self, preset: SamplingChannelPreset) -> bool {
        self.selectors
            .iter()
            .any(|selector| selector.preset() == Some(preset))
    }

    pub fn named(&self, name: &str) -> Option<&ResolvedNamedSamplingChannel> {
        self.named_channels
            .iter()
            .find(|channel| channel.name == name)
    }
}

/// Errors from selector parsing or graph-scoped resolution.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SamplingSelectionError {
    EmptyGraphName,
    EmptySelector,
    UnknownPreset(String),
    MissingChannelDefinition {
        graph: String,
        channel: String,
        available: Vec<String>,
    },
    MissingParentLmb {
        graph: String,
        channel: String,
    },
    InvalidChannelDefinition {
        graph: String,
        channel: String,
        error: String,
    },
}

impl fmt::Display for SamplingSelectionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyGraphName => {
                formatter.write_str("sampling channel selection requires a graph name")
            }
            Self::EmptySelector => {
                formatter.write_str("sampling channel selection contains an empty selector")
            }
            Self::UnknownPreset(preset) => write!(
                formatter,
                "unknown sampling channel preset `{preset}`; expected auto:lmb, auto:optimized_lmb or auto:surfaces"
            ),
            Self::MissingChannelDefinition {
                graph,
                channel,
                available,
            } => write!(
                formatter,
                "sampling channel `{graph}.{channel}` was selected but has no channel definition; available definitions: {available:?}"
            ),
            Self::MissingParentLmb { graph, channel } => write!(
                formatter,
                "sampling channel `{graph}.{channel}` must declare a non-empty parent_lmb"
            ),
            Self::InvalidChannelDefinition {
                graph,
                channel,
                error,
            } => write!(
                formatter,
                "sampling channel `{graph}.{channel}` has an invalid around expression: {error}"
            ),
        }
    }
}

impl std::error::Error for SamplingSelectionError {}

/// Resolve the selectors that apply to `graph_name`.
///
/// A graph-specific entry replaces the default list. If no entry exists for
/// this graph, the default list applies. Entries within the selected list are
/// additive and deduplicated in declaration order.
pub fn resolve_sampling_channel_selection(
    graph_name: &str,
    selection: &SamplingChannelSelection,
) -> Result<ResolvedSamplingChannelSelection, SamplingSelectionError> {
    resolve_selection(graph_name, selection, true)
}

/// Resolve with an additive default plus graph-specific list. This is useful
/// for callers that intentionally construct a shared baseline catalogue; the
/// runtime TOML contract uses replacement semantics above.
pub fn resolve_sampling_channel_selection_with_defaults(
    graph_name: &str,
    selection: &SamplingChannelSelection,
) -> Result<ResolvedSamplingChannelSelection, SamplingSelectionError> {
    resolve_selection(graph_name, selection, false)
}

/// Resolve using the explicit graph entry as a replacement for the defaults.
///
/// This is retained as a named compatibility entry point for callers that
/// want to state the replacement policy explicitly.
pub fn resolve_sampling_channel_selection_replacing_default(
    graph_name: &str,
    selection: &SamplingChannelSelection,
) -> Result<ResolvedSamplingChannelSelection, SamplingSelectionError> {
    resolve_selection(graph_name, selection, true)
}

fn resolve_selection(
    graph_name: &str,
    selection: &SamplingChannelSelection,
    replace_default: bool,
) -> Result<ResolvedSamplingChannelSelection, SamplingSelectionError> {
    let graph_name = graph_name.trim();
    if graph_name.is_empty() {
        return Err(SamplingSelectionError::EmptyGraphName);
    }

    let graph_selectors = selection.channel_selection.get(graph_name);
    let raw_selectors = if replace_default {
        graph_selectors
            .map_or(
                selection.default_channel_selection.as_slice(),
                Vec::as_slice,
            )
            .iter()
            .collect::<Vec<_>>()
    } else {
        selection
            .default_channel_selection
            .iter()
            .chain(graph_selectors.into_iter().flatten())
            .collect::<Vec<_>>()
    };
    let mut selectors = Vec::with_capacity(raw_selectors.len());
    for raw in raw_selectors {
        let selector = SamplingChannelSelector::parse(raw)?;
        if !selectors.contains(&selector) {
            selectors.push(selector);
        }
    }

    let definitions = selection
        .channel_definitions
        .get(graph_name)
        .cloned()
        .unwrap_or_default();
    let available = definitions.keys().cloned().collect::<Vec<_>>();
    let mut named_channels = Vec::new();
    for selector in &selectors {
        let SamplingChannelSelector::Named(name) = selector else {
            continue;
        };
        let Some(definition) = definitions.get(name) else {
            return Err(SamplingSelectionError::MissingChannelDefinition {
                graph: graph_name.to_owned(),
                channel: name.clone(),
                available: available.clone(),
            });
        };
        if definition.parent_lmb.is_empty() {
            return Err(SamplingSelectionError::MissingParentLmb {
                graph: graph_name.to_owned(),
                channel: name.clone(),
            });
        }
        let map = SamplingMapDefinition::parse(&definition.around).map_err(|error| {
            SamplingSelectionError::InvalidChannelDefinition {
                graph: graph_name.to_owned(),
                channel: name.clone(),
                error: error.to_string(),
            }
        })?;
        named_channels.push(ResolvedNamedSamplingChannel {
            name: name.clone(),
            definition: definition.clone(),
            map,
        });
    }

    Ok(ResolvedSamplingChannelSelection {
        graph_name: graph_name.to_owned(),
        selectors,
        named_channels,
    })
}

/// Return graph names having explicit selection entries, in deterministic order.
pub fn explicitly_selected_graphs(selection: &SamplingChannelSelection) -> Vec<&str> {
    selection
        .channel_selection
        .keys()
        .map(String::as_str)
        .collect()
}

/// Return a copy of graph definitions for diagnostics and catalogue builders.
pub fn graph_channel_definitions(
    selection: &SamplingChannelSelection,
    graph_name: &str,
) -> BTreeMap<String, SamplingChannelDefinition> {
    selection
        .channel_definitions
        .get(graph_name)
        .cloned()
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::settings::runtime::ParameterizationSettings;

    fn definition(around: &str) -> SamplingChannelDefinition {
        SamplingChannelDefinition {
            around: around.to_owned(),
            parent_lmb: vec![1, 2],
            on_cut: vec![],
        }
    }

    #[test]
    fn graph_entry_replaces_default_and_deduplicates_entries() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["auto:lmb".into(), "common".into()],
            ..Default::default()
        };
        selection.channel_selection.insert(
            "G".into(),
            vec!["auto:surfaces".into(), "named".into(), "named".into()],
        );
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("named".into(), definition("surface(1,2)"));

        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        assert_eq!(resolved.graph_name, "G");
        assert_eq!(
            resolved.selectors,
            vec![
                SamplingChannelSelector::Preset(SamplingChannelPreset::Surfaces),
                SamplingChannelSelector::Named("named".into())
            ]
        );
        assert!(resolved.named("named").is_some());
    }

    #[test]
    fn additive_resolver_keeps_default_and_graph_entries() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["auto:lmb".into(), "common".into()],
            ..Default::default()
        };
        selection.channel_selection.insert(
            "G".into(),
            vec!["auto:surfaces".into(), "common".into(), "named".into()],
        );
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("common".into(), definition("lmb(1,2)"));
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("named".into(), definition("surface(1,2)"));

        let resolved = resolve_sampling_channel_selection_with_defaults("G", &selection).unwrap();
        assert_eq!(
            resolved.selectors,
            vec![
                SamplingChannelSelector::Preset(SamplingChannelPreset::Lmb),
                SamplingChannelSelector::Named("common".into()),
                SamplingChannelSelector::Preset(SamplingChannelPreset::Surfaces),
                SamplingChannelSelector::Named("named".into())
            ]
        );
    }

    #[test]
    fn replacement_variant_ignores_defaults() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["auto:lmb".into()],
            ..Default::default()
        };
        selection
            .channel_selection
            .insert("G".into(), vec!["auto:surfaces".into()]);
        let resolved =
            resolve_sampling_channel_selection_replacing_default("G", &selection).unwrap();
        assert_eq!(
            resolved.selectors,
            vec![SamplingChannelSelector::Preset(
                SamplingChannelPreset::Surfaces
            )]
        );
    }

    #[test]
    fn default_is_used_when_graph_has_no_override() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["auto:optimized_lmb".into(), "named".into()],
            ..Default::default()
        };
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("named".into(), definition("surface(1)"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        assert!(resolved.has_preset(SamplingChannelPreset::OptimizedLmb));
        assert_eq!(resolved.named_channels.len(), 1);
    }

    #[test]
    fn missing_named_definition_reports_graph_and_available_names() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["missing".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("available".into(), definition("surface(1)"));
        let error = resolve_sampling_channel_selection("G", &selection).unwrap_err();
        assert_eq!(
            error,
            SamplingSelectionError::MissingChannelDefinition {
                graph: "G".into(),
                channel: "missing".into(),
                available: vec!["available".into()]
            }
        );
    }

    #[test]
    fn named_channel_requires_explicit_parent_lmb() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["named".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert(
                "named".into(),
                SamplingChannelDefinition {
                    around: "lmb(1,2)".into(),
                    parent_lmb: Vec::new(),
                    on_cut: Vec::new(),
                },
            );
        assert!(matches!(
            resolve_sampling_channel_selection("G", &selection),
            Err(SamplingSelectionError::MissingParentLmb { .. })
        ));
    }

    #[test]
    fn unknown_auto_preset_is_rejected() {
        assert_eq!(
            SamplingChannelSelector::parse("auto:not_a_mode").unwrap_err(),
            SamplingSelectionError::UnknownPreset("auto:not_a_mode".into())
        );
    }

    #[test]
    fn catalogue_expands_presets_and_keeps_named_diagnostics() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["auto:optimized_lmb".into(), "surface_hz".into()],
            ..Default::default()
        };
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("surface_hz".into(), definition("surface(2,4)"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue =
            build_sampling_channel_catalogue(&resolved, &[(0, vec![1, 2]), (1, vec![2, 4])], &[1]);
        assert_eq!(catalogue.entries.len(), 2);
        assert_eq!(catalogue.lmb_entries().next(), Some((1, &[2, 4][..])));
        assert_eq!(catalogue.named_entries().next().unwrap().name, "surface_hz");
        assert!(catalogue.inspection_rows()[0].contains("basis=1"));
        assert!(catalogue.inspection_rows()[1].contains("parent_lmb=[1, 2]"));
    }

    #[test]
    fn inspection_reports_resolved_selectors_and_canonical_rows() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:optimized_lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[(0, vec![1, 2])], &[0]);

        let report = catalogue.inspection();
        assert_eq!(report.graph_name, "G");
        assert_eq!(report.selectors, vec!["auto:optimized_lmb"]);
        assert_eq!(report.entries, catalogue.inspection_rows());
        assert_eq!(
            report.entries,
            vec!["0: lmb basis=0 edges=[1, 2] source=auto:optimized_lmb"]
        );
    }

    #[test]
    fn surface_preset_has_optimized_lmb_fallback() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:surfaces".into(), "auto:lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue =
            build_sampling_channel_catalogue(&resolved, &[(0, vec![1]), (1, vec![2])], &[0, 1]);
        assert_eq!(catalogue.lmb_entries().count(), 2);
    }

    #[test]
    fn surface_preset_appends_master_frame_surface_candidates() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:surfaces".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue_with_surfaces(
            &resolved,
            &[(0, vec![1, 2])],
            &[0],
            &[vec![1, 2], vec![2, 3]],
            &[1, 2],
        );
        assert!(catalogue.entries.iter().any(|entry| matches!(
            entry,
            SamplingCatalogueEntry::Surface { edges, parent_lmb }
                if edges == &vec![1, 2] && parent_lmb == &vec![1, 2]
        )));
    }

    #[test]
    fn catalogue_compiles_lmb_and_surface_with_master_embedding() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["auto:lmb".into(), "threshold".into()],
            ..Default::default()
        };
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("threshold".into(), definition("surface(1,2)"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[(0, vec![1, 2])], &[0]);
        let mut context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context.surfaces.insert(
            vec![1, 2],
            SamplingSurfaceGeometry {
                center: vec![0.0; 6],
                threshold_radius: Some(3.0),
                beta: 2.0,
                power: 1.0,
            },
        );
        let compiled = catalogue.compile(&context).unwrap();
        assert_eq!(compiled.len(), 2);
        assert!(matches!(compiled[0].map, CompiledSamplingMap::Lmb(_)));
        assert_eq!(compiled[0].embedded_edges, vec![1, 2]);
        assert!(matches!(compiled[1].map, CompiledSamplingMap::Surface(_)));
        assert_eq!(compiled[1].master_graph, "G");
        assert_eq!(compiled[1].dimensions(), 6);
    }

    #[test]
    fn catalogue_rejects_non_parent_lmb_until_affine_routing_is_compiled() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue =
            build_sampling_channel_catalogue(&resolved, &[(0, vec![1, 2]), (1, vec![2, 4])], &[0]);
        let context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        let error = catalogue.compile(&context).unwrap_err();
        assert!(matches!(
            error,
            SamplingChannelCompileError::MissingLmbFrameMap { .. }
        ));
        assert!(error.to_string().contains("no affine frame map"));
    }

    #[test]
    fn catalogue_wraps_non_parent_lmb_in_supplied_affine_frame() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[(0, vec![2, 4])], &[0]);
        let mut context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context.lmb_frame_maps.insert(
            0,
            SamplingMapAffine::new(
                vec![
                    vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                ],
                vec![0.5; 6],
            )
            .unwrap(),
        );
        let compiled = catalogue.compile(&context).unwrap();
        assert!(matches!(
            compiled[0].map,
            CompiledSamplingMap::AffineLmb { .. }
        ));
        assert_eq!(compiled[0].embedded_edges, vec![1, 2]);
        let point = compiled[0]
            .forward(&[0.31, 0.42, 0.57, 0.23, 0.68, 0.81])
            .unwrap();
        assert!(point.residual < 1.0e-10);
        assert!((point.point[0] - 0.5).abs() > 1.0e-6);
        let inverse = compiled[0].inverse(&point.point).unwrap();
        assert!(inverse.residual < 1.0e-10);
        let bridge = SamplingChannelBridge::new(compiled.clone()).unwrap();
        let bridged = bridge
            .forward(
                SamplingChannelId::from(0),
                &[0.31, 0.42, 0.57, 0.23, 0.68, 0.81],
            )
            .unwrap();
        assert_eq!(bridged.raw_coordinates, point.point);
        assert!((bridged.partition.weight_sum() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn named_non_parent_lmb_uses_edge_keyed_affine_frame() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["named".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert(
                "named".into(),
                SamplingChannelDefinition {
                    around: "lmb(2,4)".into(),
                    parent_lmb: vec![1, 2],
                    on_cut: Vec::new(),
                },
            );
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context.lmb_frame_maps_by_edges.insert(
            vec![2, 4],
            SamplingMapAffine::new(
                vec![
                    vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                ],
                vec![0.0; 6],
            )
            .unwrap(),
        );
        let compiled = catalogue.compile(&context).unwrap();
        assert!(matches!(
            compiled[0].map,
            CompiledSamplingMap::AffineLmb { .. }
        ));
        assert_eq!(compiled[0].embedded_edges, vec![1, 2]);
    }

    #[test]
    fn catalogue_embeds_partial_surface_with_complement_lmb() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["threshold".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("threshold".into(), definition("surface(2)"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context.surfaces.insert(
            vec![2],
            SamplingSurfaceGeometry {
                center: vec![0.0; 3],
                threshold_radius: Some(3.0),
                beta: 2.0,
                power: 1.0,
            },
        );
        let compiled = catalogue.compile(&context).unwrap();
        assert_eq!(compiled.len(), 1);
        assert!(matches!(compiled[0].map, CompiledSamplingMap::Embedded(_)));
        assert_eq!(compiled[0].embedded_edges, vec![1, 2]);
        let point = compiled[0]
            .forward(&[0.31, 0.42, 0.57, 0.23, 0.68, 0.81])
            .unwrap();
        let inverse = compiled[0].inverse(&point.point).unwrap();
        assert!(inverse.residual < 1.0e-10);
    }

    #[test]
    fn named_cut_metadata_requires_matching_prepared_cut_context() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["cut_channel".into()];
        let mut cut_definition = definition("lmb(1,2)");
        cut_definition.on_cut = vec![3];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("cut_channel".into(), cut_definition);
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        let error = catalogue.compile(&context).unwrap_err();
        assert!(error.to_string().contains("on_cut [3]"));

        let mut matching = context.clone();
        matching.cut_id = Some(3);
        assert!(catalogue.compile(&matching).is_ok());
        matching.cut_id = Some(4);
        assert!(catalogue.compile(&matching).is_err());
    }

    #[test]
    fn side_qualified_channel_requires_matching_prepared_side() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["right_channel".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("right_channel".into(), definition("right(lmb(1,2))"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        let error = catalogue.compile(&context).unwrap_err();
        assert!(error.to_string().contains("requires prepared side Right"));

        let mut left = context.clone();
        left.side = Some(crate::integrands::process::SamplingCutSide::Left);
        let error = catalogue.compile(&left).unwrap_err();
        assert!(error.to_string().contains("requires prepared side Right"));
    }

    #[test]
    fn surface_compilation_requires_prepared_geometry() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["threshold".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("threshold".into(), definition("surface(1,2)"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        assert!(matches!(
            catalogue.compile(&context),
            Err(SamplingChannelCompileError::MissingSurfaceGeometry { .. })
        ));
    }

    #[test]
    fn unresolved_composite_map_is_rejected_without_context() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["joint".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("joint".into(), definition("product(surface(1), lmb(1,2))"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        assert!(matches!(
            catalogue.compile(&context),
            Err(SamplingChannelCompileError::UnsupportedMap { .. })
        ));
    }

    #[test]
    fn bridge_pushes_selected_channel_and_partitions_the_raw_frame() {
        let settings = ParameterizationSettings::default();
        let map = SamplingMapKernel::new(
            SamplingMapDefinition::Lmb(vec![1]),
            settings.clone(),
            100.0,
            1,
        )
        .unwrap();
        let second =
            SamplingMapKernel::new(SamplingMapDefinition::Lmb(vec![1]), settings, 100.0, 1)
                .unwrap();
        let bridge = SamplingChannelBridge::new(vec![
            CompiledSamplingChannel {
                name: "left".into(),
                master_graph: "G".into(),
                basis_id: Some(0),
                definition: SamplingMapDefinition::Lmb(vec![1]),
                embedded_edges: vec![1],
                map: CompiledSamplingMap::Lmb(map),
            },
            CompiledSamplingChannel {
                name: "right".into(),
                master_graph: "G".into(),
                basis_id: Some(1),
                definition: SamplingMapDefinition::Lmb(vec![1]),
                embedded_edges: vec![1],
                map: CompiledSamplingMap::Lmb(second),
            },
        ])
        .unwrap();
        let evaluation = bridge
            .forward(SamplingChannelId::from(1), &[0.31, 0.42, 0.57])
            .unwrap();
        assert_eq!(evaluation.raw_coordinates, evaluation.map.point);
        assert_eq!(evaluation.partition.weights.len(), 2);
        assert!((evaluation.partition.weight_sum() - 1.0).abs() < 1.0e-12);
        assert!(
            evaluation
                .partition
                .weights
                .iter()
                .all(|weight| *weight > 0.0)
        );
        let inverse = bridge
            .inverse(SamplingChannelId::from(1), &evaluation.raw_coordinates)
            .unwrap();
        assert!(inverse.map.residual < 1.0e-10);

        let externals = Externals::default();
        let sample = evaluation
            .to_momentum_sample::<f64>(SamplingMomentumSampleContext {
                loop_mom_cache_id: 4,
                external_moms: &externals,
                external_mom_cache_id: 7,
                dependent_momenta_constructor: DependentMomentaConstructor::CrossSection,
                orientation: Some(2),
            })
            .unwrap();
        assert_eq!(sample.loop_moms().0.len(), 1);
        assert_eq!(sample.sample.orientation, Some(2));
        assert!((sample.jacobian().0 - evaluation.map.jacobian).abs() < 1.0e-12);
    }

    #[test]
    fn canonical_lmb_bridge_integrates_normalized_gaussian_with_map_partition() {
        // Build both channels through the canonical catalogue.  The selected
        // map is sampled uniformly in channel space; multiplying by
        // `n_channels * J_i * w_i` converts that sample to the exact
        // map-density mixture estimator, where w_i = rho_i / sum_j rho_j.
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue =
            build_sampling_channel_catalogue(&resolved, &[(0, vec![1]), (1, vec![1])], &[0, 1]);
        let context = SamplingChannelCompileContext::new(
            "G",
            vec![1],
            ParameterizationSettings::default(),
            2.0,
            1,
        );
        let bridge = SamplingChannelBridge::new(catalogue.compile(&context).unwrap()).unwrap();

        let bins = 24usize;
        let normalisation = (2.0 * std::f64::consts::PI).powf(-1.5);
        let channel_count = bridge.channels().len() as f64;
        let mut integral = 0.0;
        for i in 0..bins {
            for j in 0..bins {
                for k in 0..bins {
                    let coordinates = [
                        (i as f64 + 0.5) / bins as f64,
                        (j as f64 + 0.5) / bins as f64,
                        (k as f64 + 0.5) / bins as f64,
                    ];
                    let channel = SamplingChannelId::from((i + j + k) % bridge.channels().len());
                    let evaluation = bridge.forward(channel, &coordinates).unwrap();
                    let radius_squared = evaluation
                        .raw_coordinates
                        .iter()
                        .map(|component| component.powi(2))
                        .sum::<f64>();
                    let target = normalisation * (-0.5 * radius_squared).exp();
                    let partition_weight = evaluation.partition.weight(channel.index()).unwrap();
                    integral += target * evaluation.map.jacobian * channel_count * partition_weight;
                }
            }
        }
        integral /= bins.pow(3) as f64;
        assert!(
            (integral - 1.0).abs() < 2.0e-2,
            "canonical map-density mixture integral = {integral}"
        );
    }

    #[test]
    fn bridge_rejects_partial_or_mismatched_channels() {
        let settings = ParameterizationSettings::default();
        let map = SamplingMapKernel::new(SamplingMapDefinition::Lmb(vec![1]), settings, 100.0, 1)
            .unwrap();
        let channel = CompiledSamplingChannel {
            name: "one".into(),
            master_graph: "G".into(),
            basis_id: Some(0),
            definition: SamplingMapDefinition::Lmb(vec![1]),
            embedded_edges: vec![1],
            map: CompiledSamplingMap::Lmb(map),
        };
        let bridge = SamplingChannelBridge::new(vec![channel.clone()]).unwrap();
        assert!(
            bridge
                .forward(SamplingChannelId::from(1), &[0.2, 0.3, 0.4])
                .is_err()
        );
        let mut malformed = channel;
        malformed.name = "wrong".into();
        // A map with a different loop count cannot be put into the common raw frame.
        let other = SamplingMapKernel::new(
            SamplingMapDefinition::Lmb(vec![1, 2]),
            ParameterizationSettings::default(),
            100.0,
            2,
        )
        .unwrap();
        malformed.map = CompiledSamplingMap::Lmb(other);
        assert!(matches!(
            SamplingChannelBridge::new(vec![malformed, bridge.channels()[0].clone()]),
            Err(SamplingChannelBridgeError::DimensionMismatch { .. })
        ));

        let different_frame = CompiledSamplingChannel {
            name: "different-frame".into(),
            master_graph: "G".into(),
            basis_id: Some(2),
            definition: SamplingMapDefinition::Lmb(vec![2]),
            embedded_edges: vec![2],
            map: CompiledSamplingMap::Lmb(
                SamplingMapKernel::new(
                    SamplingMapDefinition::Lmb(vec![2]),
                    ParameterizationSettings::default(),
                    100.0,
                    1,
                )
                .unwrap(),
            ),
        };
        assert!(matches!(
            SamplingChannelBridge::new(vec![bridge.channels()[0].clone(), different_frame]),
            Err(SamplingChannelBridgeError::IncompatibleFrames { .. })
        ));
    }
}
