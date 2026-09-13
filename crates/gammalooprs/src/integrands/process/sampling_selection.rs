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
    ImplicitSurfaceRadialMap, PreparedCutSamplingContext, SamplingChannelScore, SamplingMapAffine,
    SamplingMapComponent, SamplingMapComposition, SamplingMapContract, SamplingMapDefinition,
    SamplingMapEmbedding, SamplingMapEvaluation, SamplingMapKernel, SamplingPartition,
    SamplingPartitionMode, SamplingScoreFunction, SurfaceRadialMap,
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
#[derive(Clone, Debug)]
pub struct SamplingChannelCompileContext {
    pub master_graph: String,
    /// Optional graph identity used to validate prepared physical-cut
    /// contexts. A cut/side map requires this to be supplied explicitly.
    pub graph_id: Option<usize>,
    /// Complete ordered parent LMB in the master graph frame. Every named
    /// channel is resolved against this exact list before compilation.
    pub parent_lmb: Vec<usize>,
    pub parameterization_settings: ParameterizationSettings,
    pub e_cm: f64,
    pub n_loop_momenta: usize,
    /// Geometry keyed by (physical energy edges, ordered active LMB edges).
    pub surfaces: BTreeMap<(Vec<usize>, Vec<usize>), SamplingSurfaceGeometry>,
    /// Optional exact direction-dependent surface maps prepared by the
    /// process layer. Their evaluators already capture the relevant external
    /// and cut data, so compilation does not infer kinematics from edge lists.
    pub implicit_surfaces: BTreeMap<(Vec<usize>, Vec<usize>), ImplicitSurfaceRadialMap>,
    /// Host cut identity for metadata whose `on_cut` list is explicit.
    pub cut_id: Option<usize>,
    pub orientation: Option<usize>,
    pub side: Option<super::SamplingCutSide>,
    /// Kinematics prepared by the physical-cut host. This remains optional
    /// for ordinary LMB/surface channels, but is mandatory for phase-space
    /// and side-qualified maps.
    pub prepared_cut_context: Option<PreparedCutSamplingContext>,
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
            graph_id: None,
            parent_lmb,
            parameterization_settings,
            e_cm,
            n_loop_momenta,
            surfaces: BTreeMap::new(),
            implicit_surfaces: BTreeMap::new(),
            lmb_frame_maps: BTreeMap::new(),
            lmb_frame_maps_by_edges: BTreeMap::new(),
            cut_id: None,
            orientation: None,
            side: None,
            prepared_cut_context: None,
        }
    }

    /// Attach a graph/cut/orientation/side context prepared from one solved
    /// Cutkosky cut. The checks are transactional: conflicting context data
    /// leave this compile context unchanged.
    pub fn with_prepared_cut_context(
        mut self,
        prepared: PreparedCutSamplingContext,
    ) -> std::result::Result<Self, SamplingChannelCompileError> {
        self.set_prepared_cut_context(prepared)?;
        Ok(self)
    }

    /// Attach prepared cut data and populate the corresponding host fields.
    /// Named phase-space/left/right maps therefore cannot accidentally reuse
    /// a stale cut, orientation, side, parent LMB, or t* value.
    pub fn set_prepared_cut_context(
        &mut self,
        prepared: PreparedCutSamplingContext,
    ) -> std::result::Result<(), SamplingChannelCompileError> {
        validate_prepared_context_identity(self, &prepared, "<prepared-context>")?;
        self.graph_id = Some(prepared.graph_id);
        self.cut_id = Some(prepared.cut_id);
        self.orientation = prepared.orientation;
        self.side = Some(prepared.side);
        self.prepared_cut_context = Some(prepared);
        Ok(())
    }

    /// Set the graph id needed when compiling a physical-cut channel.
    pub fn with_graph_id(mut self, graph_id: usize) -> Self {
        self.graph_id = Some(graph_id);
        self
    }

    pub fn insert_implicit_surface(
        &mut self,
        surface_edges: Vec<usize>,
        subspace_lmb: Vec<usize>,
        map: ImplicitSurfaceRadialMap,
    ) -> Result<()> {
        if surface_edges.is_empty() || subspace_lmb.is_empty() {
            return Err(eyre!(
                "implicit surface map needs non-empty physical and subspace edge lists"
            ));
        }
        let mut surface_seen = BTreeSet::new();
        if surface_edges.iter().any(|edge| !surface_seen.insert(*edge)) {
            return Err(eyre!(
                "implicit surface map physical edge list contains duplicates: {surface_edges:?}"
            ));
        }
        let mut subspace_seen = BTreeSet::new();
        if subspace_lmb.iter().any(|edge| !subspace_seen.insert(*edge))
            || subspace_lmb
                .iter()
                .any(|edge| !self.parent_lmb.contains(edge))
        {
            return Err(eyre!(
                "implicit surface map subspace_lmb {subspace_lmb:?} must be unique and contained in parent LMB {:?}",
                self.parent_lmb
            ));
        }
        if map.dimension() != 3 * subspace_lmb.len() {
            return Err(eyre!(
                "implicit surface map has dimension {}, expected {} for subspace {:?}",
                map.dimension(),
                3 * subspace_lmb.len(),
                subspace_lmb
            ));
        }
        self.implicit_surfaces
            .insert((surface_edges, subspace_lmb), map);
        Ok(())
    }
}

/// The map kernels currently compilable without a graph-specific implicit
/// solver.  Physical conditional/cut maps stay in the typed catalogue until
/// their prepared kinematic context is supplied by the process layer; the
/// bounded direct-product route is compiled there when its blocks are explicit.
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
        self.forward_with_context(coordinates, &[])
    }

    pub fn forward_with_context(
        &self,
        coordinates: &[f64],
        context: &[f64],
    ) -> Result<SamplingMapEvaluation> {
        <Self as SamplingMapComponent>::forward(self, coordinates, context)
    }

    pub fn inverse(&self, point: &[f64]) -> Result<SamplingMapEvaluation> {
        self.inverse_with_context(point, &[])
    }

    pub fn inverse_with_context(
        &self,
        point: &[f64],
        context: &[f64],
    ) -> Result<SamplingMapEvaluation> {
        <Self as SamplingMapComponent>::inverse(self, point, context)
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

/// Summary of a deterministic acceptance run through a complete channel
/// bridge.  The reference integrand is a normalized isotropic Gaussian in
/// the bridge's raw master frame.  Every canonical channel is sampled with
/// equal probability and the estimator includes its exact map determinant
/// and the common inverse-density partition.  This makes the report useful
/// for acceptance tests of arbitrary graph-resolved channel catalogues,
/// including products and ordered surface compositions.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingChannelBridgeAcceptanceReport {
    pub sample_count: usize,
    pub channel_count: usize,
    pub finite_sample_count: usize,
    pub normalization: f64,
    pub normalization_stderr: f64,
    pub partition_min: f64,
    pub partition_max: f64,
    pub jacobian_min: f64,
    pub jacobian_max: f64,
}

impl SamplingChannelBridgeAcceptanceReport {
    /// Integrate a normalized Gaussian over all channels in `bridge`.
    ///
    /// The returned normalization converges to one as `sample_count` grows.
    /// A sample count applies independently to each channel; no channel is
    /// silently omitted from the acceptance run.  The deterministic Halton
    /// points keep this suitable for reproducible unit and acceptance tests.
    pub fn normalized_gaussian(
        bridge: &SamplingChannelBridge,
        sample_count: usize,
        width: f64,
        center: &[f64],
    ) -> Result<Self> {
        if sample_count == 0 {
            return Err(eyre!(
                "sampling bridge acceptance harness needs at least one sample"
            ));
        }
        if !width.is_finite() || width <= 0.0 {
            return Err(eyre!(
                "sampling bridge acceptance Gaussian width must be positive and finite"
            ));
        }
        if center.len() != bridge.dimensions {
            return Err(eyre!(
                "sampling bridge acceptance Gaussian has dimension {}, expected {}",
                center.len(),
                bridge.dimensions
            ));
        }
        if center.iter().any(|component| !component.is_finite()) {
            return Err(eyre!(
                "sampling bridge acceptance Gaussian centre must be finite"
            ));
        }

        let channel_count = bridge.channels.len();
        if channel_count == 0 {
            return Err(eyre!(
                "sampling bridge acceptance harness needs at least one channel"
            ));
        }
        let dimension = bridge.dimensions as f64;
        let gaussian_normalization =
            (2.0 * std::f64::consts::PI * width * width).powf(-0.5 * dimension);
        let mut report = Self {
            sample_count,
            channel_count,
            finite_sample_count: 0,
            normalization: 0.0,
            normalization_stderr: 0.0,
            partition_min: f64::INFINITY,
            partition_max: f64::NEG_INFINITY,
            jacobian_min: f64::INFINITY,
            jacobian_max: f64::NEG_INFINITY,
        };
        let mut square_sum = 0.0;
        let total_samples = sample_count * channel_count;
        for channel_index in 0..channel_count {
            let channel_id = SamplingChannelId::from(channel_index);
            for sample in 1..=sample_count {
                let coordinates = (0..bridge.channels[channel_index].dimensions())
                    .map(|axis| bridge_halton(sample, bridge_prime(axis)))
                    .collect::<Vec<_>>();
                let evaluation = bridge.forward(channel_id, &coordinates)?;
                let partition_sum = evaluation.partition.weights.iter().sum::<f64>();
                if !partition_sum.is_finite() {
                    return Err(eyre!(
                        "sampling bridge acceptance partition is non-finite for channel {} sample {}",
                        channel_index,
                        sample
                    ));
                }
                report.partition_min = report.partition_min.min(partition_sum);
                report.partition_max = report.partition_max.max(partition_sum);
                let radius_squared = evaluation
                    .raw_coordinates
                    .iter()
                    .zip(center)
                    .map(|(point, centre)| (point - centre).powi(2))
                    .sum::<f64>();
                let weight = gaussian_normalization
                    * (-0.5 * radius_squared / width.powi(2)).exp()
                    * evaluation.map.jacobian
                    * evaluation.partition.weight(channel_index).ok_or_else(|| {
                        eyre!("sampling bridge partition has no channel {channel_index}")
                    })?
                    * channel_count as f64;
                report.jacobian_min = report.jacobian_min.min(evaluation.map.jacobian);
                report.jacobian_max = report.jacobian_max.max(evaluation.map.jacobian);
                if !weight.is_finite() {
                    continue;
                }
                report.finite_sample_count += 1;
                report.normalization += weight;
                square_sum += weight * weight;
            }
        }
        if report.finite_sample_count > 0 {
            report.normalization /= total_samples as f64;
            report.normalization_stderr =
                ((square_sum / total_samples as f64 - report.normalization.powi(2)).max(0.0)
                    / total_samples as f64)
                    .sqrt();
        }
        Ok(report)
    }
}

fn bridge_prime(index: usize) -> u64 {
    const PRIMES: [u64; 16] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
    PRIMES.get(index).copied().unwrap_or(59 + 2 * index as u64)
}

fn bridge_halton(mut index: usize, base: u64) -> f64 {
    let mut fraction = 1.0;
    let mut value = 0.0;
    while index > 0 {
        fraction /= base as f64;
        value += fraction * (index as u64 % base) as f64;
        index /= base as usize;
    }
    value
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
        self.forward_with_context(channel_id, coordinates, &[])
    }

    /// Forward a channel while supplying the output of an earlier conditional
    /// map (for example a sampled complement block) to the selected map.
    pub fn forward_with_context(
        &self,
        channel_id: SamplingChannelId,
        coordinates: &[f64],
        context: &[f64],
    ) -> Result<SamplingChannelBridgeEvaluation> {
        let channel_index = channel_id.0;
        let channel =
            self.channels
                .get(channel_index)
                .ok_or(SamplingChannelBridgeError::UnknownChannel {
                    channel: channel_id,
                })?;
        let map = channel.forward_with_context(coordinates, context)?;
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
        self.inverse_with_context(channel_id, raw_coordinates, &[])
    }

    pub fn inverse_with_context(
        &self,
        channel_id: SamplingChannelId,
        raw_coordinates: &[f64],
        context: &[f64],
    ) -> Result<SamplingChannelBridgeEvaluation> {
        let channel_index = channel_id.0;
        let channel =
            self.channels
                .get(channel_index)
                .ok_or(SamplingChannelBridgeError::UnknownChannel {
                    channel: channel_id,
                })?;
        let map = channel.inverse_with_context(raw_coordinates, context)?;
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
        self.forward_with_context(coordinates, &[])
    }

    pub fn forward_with_context(
        &self,
        coordinates: &[f64],
        context: &[f64],
    ) -> Result<SamplingMapEvaluation> {
        self.map.forward_with_context(coordinates, context)
    }

    pub fn inverse(&self, point: &[f64]) -> Result<SamplingMapEvaluation> {
        self.inverse_with_context(point, &[])
    }

    pub fn inverse_with_context(
        &self,
        point: &[f64],
        context: &[f64],
    ) -> Result<SamplingMapEvaluation> {
        self.map.inverse_with_context(point, context)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SamplingChannelCompileError {
    EmptyMasterGraph,
    MissingPreparedCutContext {
        channel: String,
        map: String,
    },
    InvalidPreparedCutContext {
        channel: String,
        error: String,
    },
    UnsupportedMap {
        channel: String,
        map: String,
    },
    /// A soft or collinear selector was parsed successfully but cannot yet be
    /// turned into a map-density channel from edge ids alone.  These targets
    /// need a graph-resolved frame and an explicit radial/angular profile;
    /// retaining this distinction avoids accidentally treating a proxy score
    /// as the selected channel's integration Jacobian.
    UnsupportedSingularPrimitive {
        channel: String,
        primitive: String,
        reason: String,
    },
    MissingSurfaceGeometry {
        channel: String,
        edges: Vec<usize>,
        subspace_lmb: Vec<usize>,
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
            Self::MissingPreparedCutContext { channel, map } => write!(
                formatter,
                "sampling channel `{channel}` uses `{map}` but no prepared cut context was supplied (graph, cut, orientation, side, parent LMB and t* are required)"
            ),
            Self::InvalidPreparedCutContext { channel, error } => write!(
                formatter,
                "sampling channel `{channel}` has invalid prepared cut context: {error}"
            ),
            Self::UnsupportedMap { channel, map } => write!(
                formatter,
                "sampling channel `{channel}` uses map `{map}`, which needs prepared graph context before it can be compiled"
            ),
            Self::UnsupportedSingularPrimitive {
                channel,
                primitive,
                reason,
            } => write!(
                formatter,
                "sampling channel `{channel}` uses singular primitive `{primitive}`, which is not compiled: {reason}"
            ),
            Self::MissingSurfaceGeometry {
                channel,
                edges,
                subspace_lmb,
            } => write!(
                formatter,
                "sampling channel `{channel}` has no prepared surface geometry for energy edges {edges:?} in subspace_lmb {subspace_lmb:?}"
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

fn validate_prepared_context_identity(
    context: &SamplingChannelCompileContext,
    prepared: &PreparedCutSamplingContext,
    channel: &str,
) -> std::result::Result<(), SamplingChannelCompileError> {
    if context.master_graph != prepared.graph_name {
        return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
            channel: channel.to_owned(),
            error: format!(
                "prepared graph `{}` does not match master graph `{}`",
                prepared.graph_name, context.master_graph
            ),
        });
    }
    if let Some(graph_id) = context.graph_id {
        if graph_id != prepared.graph_id {
            return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
                channel: channel.to_owned(),
                error: format!(
                    "prepared graph id {} does not match compile-context graph id {}",
                    prepared.graph_id, graph_id
                ),
            });
        }
    }
    if prepared.parent_lmb != context.parent_lmb {
        return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
            channel: channel.to_owned(),
            error: format!(
                "prepared parent LMB {:?} does not match compile context {:?}",
                prepared.parent_lmb, context.parent_lmb
            ),
        });
    }
    if let Some(cut_id) = context.cut_id {
        if cut_id != prepared.cut_id {
            return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
                channel: channel.to_owned(),
                error: format!(
                    "prepared cut id {} does not match compile-context cut id {}",
                    prepared.cut_id, cut_id
                ),
            });
        }
    }
    if context.orientation.is_some() && context.orientation != prepared.orientation {
        return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
            channel: channel.to_owned(),
            error: format!(
                "prepared orientation {:?} does not match compile-context orientation {:?}",
                prepared.orientation, context.orientation
            ),
        });
    }
    if context.side.is_some() && context.side != Some(prepared.side) {
        return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
            channel: channel.to_owned(),
            error: format!(
                "prepared side {:?} does not match compile-context side {:?}",
                prepared.side, context.side
            ),
        });
    }
    if !prepared.rescaling_t_star.is_finite() || prepared.rescaling_t_star <= 0.0 {
        return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
            channel: channel.to_owned(),
            error: format!(
                "cut rescaling t* must be finite and positive, got {}",
                prepared.rescaling_t_star
            ),
        });
    }
    Ok(())
}

fn prepared_map_name(definition: &SamplingMapDefinition) -> Option<&'static str> {
    match definition {
        SamplingMapDefinition::PhaseSpace(_) => Some("phase_space"),
        SamplingMapDefinition::Left(_) => Some("left"),
        SamplingMapDefinition::Right(_) => Some("right"),
        SamplingMapDefinition::Product(maps)
        | SamplingMapDefinition::Intersect(maps)
        | SamplingMapDefinition::Then(maps) => maps.iter().find_map(prepared_map_name),
        _ => None,
    }
}

/// Validate the physical-cut boundary before any graph-independent map is
/// compiled. The numerical relocation of these maps is intentionally still a
/// later runtime milestone; accepting them without this boundary would make
/// stale cut or t* data indistinguishable from valid sampling coordinates.
fn validate_prepared_map_context(
    channel: &str,
    definition: &SamplingMapDefinition,
    context: &SamplingChannelCompileContext,
) -> std::result::Result<(), SamplingChannelCompileError> {
    let Some(map_name) = prepared_map_name(definition) else {
        return Ok(());
    };
    let Some(prepared) = context.prepared_cut_context.as_ref() else {
        return Err(SamplingChannelCompileError::MissingPreparedCutContext {
            channel: channel.to_owned(),
            map: map_name.to_owned(),
        });
    };
    validate_prepared_context_identity(context, prepared, channel)?;
    if context.graph_id.is_none() {
        return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
            channel: channel.to_owned(),
            error: "compile context omitted graph id for a physical-cut map".to_owned(),
        });
    }
    if prepared.orientation.is_none() || context.orientation.is_none() {
        return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
            channel: channel.to_owned(),
            error: "physical-cut map requires an explicit orientation".to_owned(),
        });
    }
    if context.cut_id != Some(prepared.cut_id) {
        return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
            channel: channel.to_owned(),
            error: format!(
                "physical-cut map requires cut id {}, compile context supplies {:?}",
                prepared.cut_id, context.cut_id
            ),
        });
    }
    if context.side != Some(prepared.side) {
        return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
            channel: channel.to_owned(),
            error: format!(
                "physical-cut map requires prepared side {:?}, compile context supplies {:?}",
                prepared.side, context.side
            ),
        });
    }
    match definition {
        SamplingMapDefinition::PhaseSpace(inner) => {
            if !matches!(inner.as_ref(), SamplingMapDefinition::Cut(_)) {
                return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
                    channel: channel.to_owned(),
                    error: "phase_space(...) must wrap exactly one cut(...) map".to_owned(),
                });
            }
        }
        SamplingMapDefinition::Left(inner) => {
            if context.side != Some(super::SamplingCutSide::Left) {
                return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
                    channel: channel.to_owned(),
                    error: format!(
                        "left(...) requires the left prepared side, got {:?}",
                        context.side
                    ),
                });
            }
            validate_prepared_map_context(channel, inner, context)?;
        }
        SamplingMapDefinition::Right(inner) => {
            if context.side != Some(super::SamplingCutSide::Right) {
                return Err(SamplingChannelCompileError::InvalidPreparedCutContext {
                    channel: channel.to_owned(),
                    error: format!(
                        "right(...) requires the right prepared side, got {:?}",
                        context.side
                    ),
                });
            }
            validate_prepared_map_context(channel, inner, context)?;
        }
        SamplingMapDefinition::Product(maps)
        | SamplingMapDefinition::Intersect(maps)
        | SamplingMapDefinition::Then(maps) => {
            for map in maps {
                validate_prepared_map_context(channel, map, context)?;
            }
        }
        _ => {}
    }
    Ok(())
}

/// Compile one radial surface in its native subspace and optionally embed it together
/// with the ordinary map on the complementary parent-LMB edges.  The output
/// permutation is explicit, so sharing a surface across graph channels never
/// relies on an implicit edge ordering or a silent local-frame assumption.
fn compile_surface_map(
    channel: &str,
    surface_edges: &[usize],
    edges: &[usize],
    context: &SamplingChannelCompileContext,
    embed_complement: bool,
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
    let expected_dimension = 3 * edges.len();
    let key = (surface_edges.to_vec(), edges.to_vec());
    let surface = if let Some(surface) = context.implicit_surfaces.get(&key) {
        if surface.dimension() != expected_dimension {
            return Err(SamplingChannelCompileError::InvalidChannel {
                channel: channel.to_owned(),
                error: format!(
                    "prepared implicit surface has dimension {}, expected {}",
                    surface.dimension(),
                    expected_dimension
                ),
            });
        }
        CompiledSamplingMap::ImplicitSurface(surface.clone())
    } else {
        let Some(geometry) = context.surfaces.get(&key) else {
            return Err(SamplingChannelCompileError::MissingSurfaceGeometry {
                channel: channel.to_owned(),
                edges: surface_edges.to_vec(),
                subspace_lmb: edges.to_vec(),
            });
        };
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
        CompiledSamplingMap::Surface(surface)
    };
    if !embed_complement {
        return Ok(surface);
    }
    if edges.len() == context.n_loop_momenta {
        if edges == context.parent_lmb {
            return Ok(surface);
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

#[derive(Clone, Copy)]
enum PartitionedMapKind {
    Product,
    Then,
}

impl PartitionedMapKind {
    fn name(self) -> &'static str {
        match self {
            Self::Product => "product",
            Self::Then => "then",
        }
    }
}

fn compile_partitioned_children(
    kind: PartitionedMapKind,
    channel: &str,
    maps: &[SamplingMapDefinition],
    channel_subspace_lmb: &[usize],
    context: &SamplingChannelCompileContext,
) -> Result<
    (
        Vec<Box<dyn SamplingMapComponent>>,
        Vec<Vec<usize>>,
        BTreeMap<usize, usize>,
    ),
    SamplingChannelCompileError,
> {
    let parent_positions = context
        .parent_lmb
        .iter()
        .enumerate()
        .map(|(position, edge)| (*edge, position))
        .collect::<BTreeMap<_, _>>();
    let mut used_edges = BTreeMap::<usize, usize>::new();
    let mut children = Vec::<Box<dyn SamplingMapComponent>>::with_capacity(maps.len());
    let mut child_blocks = Vec::<Vec<usize>>::with_capacity(maps.len());
    let mut surface_child = None;

    let map_name = kind.name();
    for (child_index, map) in maps.iter().enumerate() {
        let (block_edges, compiled) = match map {
            SamplingMapDefinition::Lmb(edges) | SamplingMapDefinition::Complement(edges) => {
                if edges.is_empty() {
                    return Err(SamplingChannelCompileError::InvalidChannel {
                        channel: channel.to_owned(),
                        error: format!("{map_name} child {child_index} has an empty edge block"),
                    });
                }
                let compiled = SamplingMapKernel::new(
                    map.clone(),
                    context.parameterization_settings.clone(),
                    context.e_cm,
                    edges.len(),
                )
                .map_err(|error| SamplingChannelCompileError::InvalidChannel {
                    channel: channel.to_owned(),
                    error: format!("{map_name} child {child_index} {map:?} is invalid: {error}"),
                })?;
                (edges.clone(), CompiledSamplingMap::Lmb(compiled))
            }
            SamplingMapDefinition::Surface(surface_edges) => {
                if let Some(previous) = surface_child {
                    return Err(SamplingChannelCompileError::InvalidChannel {
                        channel: channel.to_owned(),
                        error: format!(
                            "{map_name} contains more than one surface child; surface children at indices {} and {}",
                            previous, child_index
                        ),
                    });
                }
                surface_child = Some(child_index);
                if channel_subspace_lmb.is_empty() {
                    return Err(SamplingChannelCompileError::InvalidChannel {
                        channel: channel.to_owned(),
                        error: format!(
                            "{map_name} surface child requires non-empty channel subspace_lmb"
                        ),
                    });
                }
                let compiled = compile_surface_map(
                    channel,
                    surface_edges,
                    channel_subspace_lmb,
                    context,
                    false,
                )?;
                (channel_subspace_lmb.to_vec(), compiled)
            }
            unsupported => {
                return Err(unsupported_map_error(
                    channel,
                    format!("{map_name} child {child_index}: {unsupported:?}"),
                    unsupported,
                ));
            }
        };

        for edge in &block_edges {
            if !parent_positions.contains_key(edge) {
                return Err(SamplingChannelCompileError::InvalidChannel {
                    channel: channel.to_owned(),
                    error: format!(
                        "{map_name} child {child_index} block {block_edges:?} contains edge {edge}, which is outside parent LMB {:?}",
                        context.parent_lmb
                    ),
                });
            }
            if let Some(previous) = used_edges.insert(*edge, child_index) {
                return Err(SamplingChannelCompileError::InvalidChannel {
                    channel: channel.to_owned(),
                    error: format!(
                        "{map_name} child {child_index} block {block_edges:?} overlaps child {previous} on parent edge {edge}; child blocks must be disjoint and cover parent LMB {:?}",
                        context.parent_lmb
                    ),
                });
            }
        }
        child_blocks.push(block_edges);
        children.push(Box::new(compiled));
    }

    let missing = context
        .parent_lmb
        .iter()
        .copied()
        .filter(|edge| !used_edges.contains_key(edge))
        .collect::<Vec<_>>();
    if !missing.is_empty() {
        return Err(SamplingChannelCompileError::InvalidChannel {
            channel: channel.to_owned(),
            error: format!(
                "{map_name} child blocks leave parent edges {missing:?} uncovered; blocks {:?} must cover parent LMB {:?}",
                child_blocks, context.parent_lmb
            ),
        });
    }

    Ok((children, child_blocks, parent_positions))
}

/// Keep parsed soft/collinear syntax from becoming a silently approximate
/// channel.  A useful singular proposal must be resolved against the actual
/// routed vectors (and, for collinearity, a relative angular frame) and carry
/// a normalized profile.  Neither can be reconstructed from the edge labels
/// in a graph-independent compiler.  Until those ingredients are supplied,
/// report a capability error rather than installing a proxy-only map in the
/// exact-density catalogue.
fn unsupported_map_error(
    channel: &str,
    map_label: String,
    map: &SamplingMapDefinition,
) -> SamplingChannelCompileError {
    match singular_primitive(map) {
        Some(SamplingMapDefinition::Soft(edge)) => {
            SamplingChannelCompileError::UnsupportedSingularPrimitive {
                channel: channel.to_owned(),
                primitive: format!("soft({edge})"),
                reason: "a graph-resolved routed three-momentum frame and a normalized radial profile are required; edge ids alone do not define an exact push-forward density".to_owned(),
            }
        }
        Some(SamplingMapDefinition::Collinear(a, b)) => {
            SamplingChannelCompileError::UnsupportedSingularPrimitive {
                channel: channel.to_owned(),
                primitive: format!("collinear({a},{b})"),
                reason: "a graph-resolved relative angular frame, branch/support rules and a normalized angular profile are required; edge ids alone do not define an exact push-forward density".to_owned(),
            }
        }
        Some(other) => SamplingChannelCompileError::UnsupportedMap {
            channel: channel.to_owned(),
            map: format!("{}: unsupported primitive {other:?}", map_label),
        },
        None => SamplingChannelCompileError::UnsupportedMap {
            channel: channel.to_owned(),
            map: map_label,
        },
    }
}

fn singular_primitive(map: &SamplingMapDefinition) -> Option<&SamplingMapDefinition> {
    match map {
        SamplingMapDefinition::Soft(_) | SamplingMapDefinition::Collinear(_, _) => Some(map),
        SamplingMapDefinition::Product(maps)
        | SamplingMapDefinition::Intersect(maps)
        | SamplingMapDefinition::Then(maps) => maps.iter().find_map(singular_primitive),
        SamplingMapDefinition::PhaseSpace(map)
        | SamplingMapDefinition::Left(map)
        | SamplingMapDefinition::Right(map) => singular_primitive(map),
        _ => None,
    }
}

/// Compile a direct-product channel whose children explicitly partition the
/// parent LMB.  The surface child is the only child whose active coordinate
/// block comes from channel metadata (`subspace_lmb`); its `surface(...)`
/// arguments identify the physical energy constraints.  Ordinary `lmb(...)`
/// and `complement(...)` children retain their own local edge ordering.
fn compile_product_map(
    channel: &str,
    maps: &[SamplingMapDefinition],
    channel_subspace_lmb: &[usize],
    context: &SamplingChannelCompileContext,
) -> Result<CompiledSamplingMap, SamplingChannelCompileError> {
    let (children, child_blocks, parent_positions) = compile_partitioned_children(
        PartitionedMapKind::Product,
        channel,
        maps,
        channel_subspace_lmb,
        context,
    )?;

    let output_indices = child_blocks
        .iter()
        .flat_map(|block| {
            block.iter().flat_map(|edge| {
                let position = parent_positions[edge];
                [3 * position, 3 * position + 1, 3 * position + 2]
            })
        })
        .collect::<Vec<_>>();
    SamplingMapEmbedding::product(children, output_indices)
        .map(CompiledSamplingMap::Embedded)
        .map_err(|error| SamplingChannelCompileError::InvalidChannel {
            channel: channel.to_owned(),
            error: format!("product embedding is invalid: {error}"),
        })
}

/// Compile an ordered conditional channel whose children partition the parent
/// LMB in the same way as [`compile_product_map`].  The only semantic
/// difference is that each later child receives all previous output blocks as
/// context.  This is deliberately kept as a bounded compiler route: ordinary
/// `lmb`/`complement` blocks and at most one active `surface` block are
/// accepted, while physical cut and multi-surface compositions remain
/// process-layer work.
fn compile_then_map(
    channel: &str,
    maps: &[SamplingMapDefinition],
    channel_subspace_lmb: &[usize],
    context: &SamplingChannelCompileContext,
) -> Result<CompiledSamplingMap, SamplingChannelCompileError> {
    let (children, child_blocks, parent_positions) = compile_partitioned_children(
        PartitionedMapKind::Then,
        channel,
        maps,
        channel_subspace_lmb,
        context,
    )?;

    let output_indices = child_blocks
        .iter()
        .flat_map(|block| {
            block.iter().flat_map(|edge| {
                let position = parent_positions[edge];
                [3 * position, 3 * position + 1, 3 * position + 2]
            })
        })
        .collect::<Vec<_>>();
    let composition = SamplingMapComposition::then(children).map_err(|error| {
        SamplingChannelCompileError::InvalidChannel {
            channel: channel.to_owned(),
            error: format!("then composition is invalid: {error}"),
        }
    })?;
    SamplingMapEmbedding::from_composition(composition, output_indices)
        .map(CompiledSamplingMap::Embedded)
        .map_err(|error| SamplingChannelCompileError::InvalidChannel {
            channel: channel.to_owned(),
            error: format!("then embedding is invalid: {error}"),
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

    /// Compile ordinary LMB/surface entries and the bounded conditional
    /// `then(lmb|complement, surface?)` route in this catalogue.
    ///
    /// Surface geometry is looked up by its canonical master-graph edge list.
    /// Physical `cut`, `intersect`, and phase/side-qualified maps remain
    /// rejected here until the process has prepared their conditional
    /// kinematics; this avoids creating a numerically plausible map with an
    /// incorrect frame.
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
                    let map = compile_surface_map(
                        &format!("surface:{edges:?}"),
                        edges,
                        edges,
                        context,
                        true,
                    )?;
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
                    validate_prepared_map_context(&channel.name, &definition, context)?;
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
                            let subspace = &channel.definition.subspace_lmb;
                            if subspace.is_empty() {
                                return Err(SamplingChannelCompileError::InvalidChannel {
                                    channel: channel.name.clone(),
                                    error: "surface channel is missing subspace_lmb metadata"
                                        .to_owned(),
                                });
                            }
                            // `edges` identifies the physical energy constraints;
                            // the radial chart itself is solved only in the
                            // explicitly supplied active loop subspace.
                            let _ = edges;
                            compile_surface_map(&channel.name, edges, subspace, context, true)?
                        }
                        SamplingMapDefinition::Product(maps) => compile_product_map(
                            &channel.name,
                            maps,
                            &channel.definition.subspace_lmb,
                            context,
                        )?,
                        SamplingMapDefinition::Then(maps) => compile_then_map(
                            &channel.name,
                            maps,
                            &channel.definition.subspace_lmb,
                            context,
                        )?,
                        unsupported => {
                            return Err(unsupported_map_error(
                                &channel.name,
                                format!("{unsupported:?}"),
                                unsupported,
                            ));
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
                    "{index}: {} around={} subspace_lmb={:?} parent_lmb={:?} on_cut={:?}",
                    channel.name,
                    channel.definition.around,
                    channel.definition.subspace_lmb,
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
    MissingSubspaceLmb {
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
            Self::MissingSubspaceLmb { graph, channel } => write!(
                formatter,
                "sampling channel `{graph}.{channel}` contains a surface/cut map and must declare a non-empty subspace_lmb"
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

fn contains_surface_map(map: &SamplingMapDefinition) -> bool {
    match map {
        SamplingMapDefinition::Surface(_) | SamplingMapDefinition::Cut(_) => true,
        SamplingMapDefinition::Product(maps)
        | SamplingMapDefinition::Intersect(maps)
        | SamplingMapDefinition::Then(maps) => maps.iter().any(contains_surface_map),
        SamplingMapDefinition::PhaseSpace(map)
        | SamplingMapDefinition::Left(map)
        | SamplingMapDefinition::Right(map) => contains_surface_map(map),
        _ => false,
    }
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
        let mut parent_seen = BTreeSet::new();
        if definition
            .parent_lmb
            .iter()
            .any(|edge| !parent_seen.insert(*edge))
        {
            return Err(SamplingSelectionError::InvalidChannelDefinition {
                graph: graph_name.to_owned(),
                channel: name.clone(),
                error: format!(
                    "parent_lmb {:?} must contain unique ordered edge ids",
                    definition.parent_lmb
                ),
            });
        }
        let map = SamplingMapDefinition::parse(&definition.around).map_err(|error| {
            SamplingSelectionError::InvalidChannelDefinition {
                graph: graph_name.to_owned(),
                channel: name.clone(),
                error: error.to_string(),
            }
        })?;
        if contains_surface_map(&map) {
            if definition.subspace_lmb.is_empty() {
                return Err(SamplingSelectionError::MissingSubspaceLmb {
                    graph: graph_name.to_owned(),
                    channel: name.clone(),
                });
            }
            let mut seen = BTreeSet::new();
            if definition
                .subspace_lmb
                .iter()
                .any(|edge| !seen.insert(*edge) || !definition.parent_lmb.contains(edge))
            {
                return Err(SamplingSelectionError::InvalidChannelDefinition {
                    graph: graph_name.to_owned(),
                    channel: name.clone(),
                    error: format!(
                        "subspace_lmb {:?} must be unique and contained in parent_lmb {:?}",
                        definition.subspace_lmb, definition.parent_lmb
                    ),
                });
            }
        }
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
    use std::sync::Arc;

    use super::*;
    use crate::integrands::process::{SamplingCutSide, SamplingSupport};
    use crate::settings::runtime::ParameterizationSettings;

    fn definition(around: &str) -> SamplingChannelDefinition {
        SamplingChannelDefinition {
            around: around.to_owned(),
            subspace_lmb: vec![1, 2],
            parent_lmb: vec![1, 2],
            on_cut: vec![],
        }
    }

    fn prepared_cut_context(
        side: SamplingCutSide,
        orientation: Option<usize>,
        parent_lmb: Vec<usize>,
    ) -> PreparedCutSamplingContext {
        PreparedCutSamplingContext::new(
            "G",
            17,
            3,
            orientation,
            side,
            parent_lmb,
            0.75,
            vec![[0.0, 0.0, 0.0]],
            vec![[100.0, 0.0, 0.0, 100.0]],
            vec![],
        )
        .unwrap()
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
                    subspace_lmb: Vec::new(),
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
    fn named_surface_requires_explicit_active_subspace() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["threshold".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert(
                "threshold".into(),
                SamplingChannelDefinition {
                    around: "surface(2,4)".into(),
                    subspace_lmb: Vec::new(),
                    parent_lmb: vec![1, 2, 4],
                    on_cut: Vec::new(),
                },
            );
        assert!(matches!(
            resolve_sampling_channel_selection("G", &selection),
            Err(SamplingSelectionError::MissingSubspaceLmb { .. })
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
            (vec![1, 2], vec![1, 2]),
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
    fn singular_primitives_fail_with_actionable_exact_density_diagnostics() {
        for (name, around, expected_primitive, expected_requirement) in [
            (
                "soft_target",
                "soft(1)",
                "soft(1)",
                "routed three-momentum frame",
            ),
            (
                "collinear_target",
                "collinear(1,2)",
                "collinear(1,2)",
                "relative angular frame",
            ),
            (
                "nested_soft_target",
                "then(lmb(1,2), soft(1))",
                "soft(1)",
                "normalized radial profile",
            ),
        ] {
            let mut selection = SamplingChannelSelection::default();
            selection.default_channel_selection = vec![name.to_owned()];
            selection
                .channel_definitions
                .entry("G".to_owned())
                .or_default()
                .insert(name.to_owned(), definition(around));
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
            match error {
                SamplingChannelCompileError::UnsupportedSingularPrimitive {
                    primitive,
                    reason,
                    ..
                } => {
                    assert_eq!(primitive, expected_primitive);
                    assert!(reason.contains(expected_requirement), "{reason}");
                    assert!(reason.contains("edge ids alone"), "{reason}");
                }
                other => panic!("expected singular primitive diagnostic, got {other:?}"),
            }
        }
    }

    #[test]
    fn catalogue_accepts_prepared_exact_implicit_surface_map() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["threshold".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("threshold".into(), definition("surface(1,2)"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context
            .insert_implicit_surface(
                vec![1, 2],
                vec![1, 2],
                ImplicitSurfaceRadialMap::new(
                    6,
                    vec![0.0; 6],
                    2.0,
                    2.0,
                    Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
                )
                .unwrap(),
            )
            .unwrap();
        let compiled = catalogue.compile(&context).unwrap();
        assert!(matches!(
            compiled[0].map,
            CompiledSamplingMap::ImplicitSurface(_)
        ));
        let bridge = SamplingChannelBridge::new(compiled).unwrap();
        let mapped = bridge.forward(
            SamplingChannelId::from(0),
            &[0.31, 0.42, 0.57, 0.23, 0.68, 0.81],
        );
        assert!(mapped.is_ok());
    }

    #[test]
    fn catalogue_compiles_then_with_context_surface_and_bridges_full_composition() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["conditional".into()];
        let mut channel = definition("then(lmb(1),surface(2,4))");
        channel.subspace_lmb = vec![2];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("conditional".into(), channel);
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context
            .insert_implicit_surface(
                vec![2, 4],
                vec![2],
                ImplicitSurfaceRadialMap::new(
                    3,
                    vec![0.0; 3],
                    1.0,
                    1.0,
                    Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
                )
                .unwrap()
                .with_context_evaluator(Arc::new(|_, radius, context| {
                    let shift = context.first().copied().unwrap_or_default().abs().min(0.2);
                    Ok((radius - (1.0 + shift), 1.0))
                })),
            )
            .unwrap();
        let compiled = catalogue.compile(&context).unwrap();
        assert_eq!(compiled.len(), 1);
        let channel = &compiled[0];
        assert_eq!(channel.embedded_edges, vec![1, 2]);
        assert_eq!(channel.map.contract().support, SamplingSupport::Full);
        let coordinates = [0.31, 0.42, 0.57, 0.23, 0.68, 0.81];
        let mapped = channel.map.forward(&coordinates).unwrap();
        assert_eq!(mapped.support, SamplingSupport::Full);
        let inverse = channel.map.inverse(&mapped.point).unwrap();
        assert!(inverse.residual < 1.0e-9, "{}", inverse.residual);
        assert!(
            inverse
                .coordinates
                .iter()
                .zip(coordinates)
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-9)
        );
        let bridge = SamplingChannelBridge::new(compiled).unwrap();
        let bridged = bridge
            .forward(SamplingChannelId::from(0), &coordinates)
            .unwrap();
        assert!((bridged.partition.weight_sum() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn distinct_physical_surfaces_share_one_active_subspace_without_overwrite() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["h".into(), "z".into()];
        let mut h = definition("surface(2,4)");
        h.subspace_lmb = vec![1, 2];
        let mut z = definition("surface(3,10)");
        z.subspace_lmb = vec![1, 2];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .extend([(String::from("h"), h), (String::from("z"), z)]);
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        for edges in [vec![2, 4], vec![3, 10]] {
            context
                .insert_implicit_surface(
                    edges,
                    vec![1, 2],
                    ImplicitSurfaceRadialMap::new(
                        6,
                        vec![0.0; 6],
                        2.0,
                        2.0,
                        Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
                    )
                    .unwrap(),
                )
                .unwrap();
        }
        let compiled = catalogue.compile(&context).unwrap();
        assert_eq!(compiled.len(), 2);
        assert_eq!(compiled[0].name, "h");
        assert_eq!(compiled[1].name, "z");
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
                    subspace_lmb: Vec::new(),
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
            .insert("threshold".into(), {
                let mut channel = definition("surface(2)");
                channel.subspace_lmb = vec![2];
                channel
            });
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context
            .insert_implicit_surface(
                vec![2],
                vec![2],
                ImplicitSurfaceRadialMap::new(
                    3,
                    vec![0.0; 3],
                    2.0,
                    1.0,
                    Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
                )
                .unwrap(),
            )
            .unwrap();
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
    fn phase_space_channel_requires_prepared_cut_context() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["cut_chart".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("cut_chart".into(), definition("phase_space(cut(4,7))"));
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
        assert!(matches!(
            error,
            SamplingChannelCompileError::MissingPreparedCutContext { .. }
        ));
        assert!(error.to_string().contains("graph, cut, orientation"));
    }

    #[test]
    fn prepared_cut_context_populates_and_validates_host_metadata() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["cut_chart".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("cut_chart".into(), definition("phase_space(cut(4,7))"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        )
        .with_prepared_cut_context(prepared_cut_context(
            SamplingCutSide::Left,
            Some(2),
            vec![1, 2],
        ))
        .unwrap();
        assert_eq!(context.graph_id, Some(17));
        assert_eq!(context.cut_id, Some(3));
        assert_eq!(context.orientation, Some(2));
        assert_eq!(context.side, Some(SamplingCutSide::Left));
        assert!(matches!(
            catalogue.compile(&context),
            Err(SamplingChannelCompileError::UnsupportedMap { .. })
        ));

        let mismatch = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        )
        .with_prepared_cut_context(prepared_cut_context(
            SamplingCutSide::Left,
            Some(2),
            vec![1, 3],
        ))
        .unwrap_err();
        assert!(mismatch.to_string().contains("parent LMB"));
    }

    #[test]
    fn side_map_rejects_missing_orientation_and_wrong_prepared_side() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["left_target".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("left_target".into(), definition("left(lmb(1,2))"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let missing_orientation = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        )
        .with_prepared_cut_context(prepared_cut_context(
            SamplingCutSide::Left,
            None,
            vec![1, 2],
        ))
        .unwrap();
        let error = catalogue.compile(&missing_orientation).unwrap_err();
        assert!(error.to_string().contains("explicit orientation"));

        let wrong_side = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        )
        .with_prepared_cut_context(prepared_cut_context(
            SamplingCutSide::Right,
            Some(2),
            vec![1, 2],
        ))
        .unwrap();
        let error = catalogue.compile(&wrong_side).unwrap_err();
        assert!(error.to_string().contains("requires prepared side Left"));
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
    fn composite_map_with_overlapping_blocks_is_rejected_with_context_diagnostic() {
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
        let mut context = context;
        context.surfaces.insert(
            (vec![1], vec![1, 2]),
            SamplingSurfaceGeometry {
                center: vec![0.0; 6],
                threshold_radius: Some(3.0),
                beta: 2.0,
                power: 1.0,
            },
        );
        let error = catalogue.compile(&context).unwrap_err().to_string();
        assert!(error.contains("overlaps child"));
        assert!(error.contains("parent LMB [1, 2]"));
    }

    #[test]
    fn product_surface_and_complement_compile_as_one_master_frame_channel() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["joint".into()];
        let mut channel = definition("product(surface(2), complement(1))");
        channel.subspace_lmb = vec![2];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("joint".into(), channel);
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context
            .insert_implicit_surface(
                vec![2],
                vec![2],
                ImplicitSurfaceRadialMap::new(
                    3,
                    vec![0.0; 3],
                    2.0,
                    1.0,
                    Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
                )
                .unwrap(),
            )
            .unwrap();
        let compiled = catalogue.compile(&context).unwrap();
        assert!(matches!(compiled[0].map, CompiledSamplingMap::Embedded(_)));
        assert_eq!(compiled[0].embedded_edges, vec![1, 2]);
        assert_eq!(
            compiled[0].map.contract().support,
            crate::integrands::process::SamplingSupport::Full
        );
        let point = compiled[0]
            .forward(&[0.31, 0.42, 0.57, 0.23, 0.68, 0.81])
            .unwrap();
        let inverse = compiled[0].inverse(&point.point).unwrap();
        assert!(point.residual < 1.0e-10);
        assert!(inverse.residual < 1.0e-10);
        let bridge = SamplingChannelBridge::new(compiled).unwrap();
        let bridged = bridge
            .forward(
                SamplingChannelId::from(0),
                &[0.31, 0.42, 0.57, 0.23, 0.68, 0.81],
            )
            .unwrap();
        assert!((bridged.partition.weight_sum() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn product_rejects_overlapping_child_blocks_with_parent_diagnostic() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["joint".into()];
        let mut channel = definition("product(surface(2), complement(2))");
        channel.subspace_lmb = vec![2];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("joint".into(), channel);
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
            (vec![2], vec![2]),
            SamplingSurfaceGeometry {
                center: vec![0.0; 3],
                threshold_radius: Some(3.0),
                beta: 2.0,
                power: 1.0,
            },
        );
        let error = catalogue.compile(&context).unwrap_err();
        let diagnostic = error.to_string();
        assert!(diagnostic.contains("overlaps child 0"));
        assert!(diagnostic.contains("parent LMB [1, 2]"));
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
    fn bridge_acceptance_report_integrates_every_canonical_channel() {
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
        let report = SamplingChannelBridgeAcceptanceReport::normalized_gaussian(
            &bridge,
            1024,
            1.0,
            &[0.0, 0.0, 0.0],
        )
        .unwrap();
        assert_eq!(report.sample_count, 1024);
        assert_eq!(report.channel_count, 2);
        assert_eq!(report.finite_sample_count, 2048);
        assert!((report.partition_min - 1.0).abs() < 1.0e-12);
        assert!((report.partition_max - 1.0).abs() < 1.0e-12);
        assert!((report.normalization - 1.0).abs() < 5.0e-2);
        assert!(report.normalization_stderr.is_finite());
    }

    #[test]
    fn bridge_acceptance_report_supports_a_single_canonical_channel() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[(0, vec![1])], &[0]);
        let context = SamplingChannelCompileContext::new(
            "G",
            vec![1],
            ParameterizationSettings::default(),
            2.0,
            1,
        );
        let bridge = SamplingChannelBridge::new(catalogue.compile(&context).unwrap()).unwrap();
        let report = SamplingChannelBridgeAcceptanceReport::normalized_gaussian(
            &bridge,
            512,
            1.0,
            &[0.0, 0.0, 0.0],
        )
        .unwrap();
        assert_eq!(report.channel_count, 1);
        assert_eq!(report.finite_sample_count, report.sample_count);
        assert!((report.normalization - 1.0).abs() < 5.0e-2);
        assert!((report.partition_min - 1.0).abs() < 1.0e-12);
        assert!((report.partition_max - 1.0).abs() < 1.0e-12);
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
