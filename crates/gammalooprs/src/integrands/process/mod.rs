use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

use crate::cff::expression::OrientationID;
use crate::graph::{FeynmanGraph, Graph, GraphGroup, GroupId, LmbIndex, LoopMomentumBasis};
use crate::integrands::evaluation::{
    EvaluationMetaData, EvaluationResult, GenericEvaluationResult, GraphEvaluationResult,
    LoopMomentaEscalationMetrics, PreciseEvaluationResult, RawBatchEvaluationResult,
    RawPreciseBatchEvaluationResult, RotatedEvaluation, StabilityFailureReason, StabilityResult,
    StabilityStatus, StatisticsCounter,
};
use crate::model::Model;
use crate::momentum::sample::{BareMomentumSample, LoopIndex, LoopMomenta, MomentumSample};
use crate::momentum::{Rotation, ThreeMomentum};
use crate::observables::{
    AdditionalWeightKey, EventProcessingRuntime, GenericEvent, HistogramProcessInfo,
    ObservableAccumulatorBundle, ObservableFileFormat, ObservableSnapshotBundle,
};
use crate::processes::{CutGroupId, GraphGroupSelectionSpec, StandaloneExportSettings};
use crate::subtraction::lu_counterterm::LUSharedOverlaps;
use crate::utils::{
    ArbPrec, F, FloatLike, RuntimeCache, SamplingFloat, SamplingPrecision, f128,
    format_for_compare_digits, get_n_dim_for_n_loop_momenta,
};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::owo_colors::OwoColorize;
use colored::Colorize;
use enum_dispatch::enum_dispatch;
use eyre::{Context, eyre};
use gammaloop_sample::{DiscreteGraphSample, GammaLoopSample, parameterize};
use itertools::Itertools;
use linnet::half_edge::involution::EdgeVec;
use linnet::half_edge::involution::Orientation;
use linnet::half_edge::subgraph::{SubSetLike, subset::SubSet};
use momtrop::SampleGenerator;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use spenso::algebra::algebraic_traits::IsZero;
use spenso::algebra::complex::Complex;
use std::sync::Once;
use std::time::{Duration, Instant};
use symbolica::numerical_integration::{ContinuousGrid, DiscreteGrid, Grid, Sample};
use tracing::{debug, warn};
use typed_index_collections::TiVec;
pub mod amplitude;
pub mod cache_debugging;
pub mod cross_section;
pub mod gammaloop_sample;
pub mod ir;
pub mod sampling_context;
pub mod sampling_joint;
pub mod sampling_maps;
pub mod sampling_partition;
pub mod sampling_reference;
pub mod sampling_selection;
use crate::{
    DependentMomentaConstructor, GammaLoopContext,
    settings::RuntimeSettings,
    settings::runtime::DiscreteGraphSamplingSettings,
    settings::runtime::DiscreteGraphSamplingType,
    settings::runtime::IntegratorSettings,
    settings::runtime::ParameterizationSettings,
    settings::runtime::Precision,
    settings::runtime::StabilityLevelSetting,
    settings::runtime::StabilitySettings,
    settings::runtime::{IntegratedPhase, SamplingChannelWeight, SamplingSettings},
};
use color_eyre::Result;
use sampling_context::PreparedLUHost;
use sampling_selection::SamplingChannelPrograms;

pub mod evaluators;
pub use evaluators::ActiveF64Backend;
pub use evaluators::{GenericEvaluator, GenericEvaluatorFloat};
pub mod sampling_evaluator;
pub use sampling_evaluator::{SamplingDualValue, SamplingExpressionEvaluator};

pub mod param_builder;
pub use param_builder::{ParamBuilder, ParamValuePairs, ThresholdParams, UpdateAndGetParams};
pub use sampling_context::{PreparedSurfaceStatus, SamplingCutSide};
pub use sampling_joint::{SharedEnergyJointGeometry, SharedEnergyJointMap};
pub use sampling_maps::{
    ImplicitSurfaceContextPreparer, ImplicitSurfaceRadialContextEvaluator,
    ImplicitSurfaceRadialEvaluator, ImplicitSurfaceRadialMap, SamplingJacobian,
    SamplingMapAcceptanceReport, SamplingMapAffine, SamplingMapComponent, SamplingMapComposition,
    SamplingMapContextTransform, SamplingMapContract, SamplingMapDefinition, SamplingMapEmbedding,
    SamplingMapEvaluation, SamplingMapKernel, SamplingMapPoint, SamplingSupport, SurfaceRadialMap,
    SurfaceRadialPoint,
};
pub use sampling_partition::{
    SamplingChannelScore, SamplingPartition, SamplingPartitionMode, SamplingScoreFunction,
};
pub use sampling_reference::{
    GaussianReferenceFunction, ReferenceMoments, ReferenceSampleEvaluation, ReferenceSamplingReport,
};
pub use sampling_selection::{
    CompiledSamplingChannel, CompiledSamplingMap, ResolvedNamedSamplingChannel,
    ResolvedSamplingBlock, ResolvedSamplingChannelSelection, SamplingCatalogueEntry,
    SamplingChannelBridge, SamplingChannelBridgeAcceptanceReport, SamplingChannelBridgeError,
    SamplingChannelBridgeEvaluation, SamplingChannelCatalogue, SamplingChannelCompileContext,
    SamplingChannelCompileError, SamplingChannelId, SamplingChannelInspection,
    SamplingChannelPreset, SamplingChannelRuntimeContexts, SamplingChannelSelector,
    SamplingCoverageReport, SamplingGeometryKey, SamplingMomentumSampleContext,
    SamplingSelectionError, build_sampling_channel_catalogue,
    build_sampling_channel_catalogue_with_surfaces,
    build_sampling_channel_catalogue_with_surfaces_and_coverage, explicitly_selected_graphs,
    graph_channel_definitions, resolve_sampling_channel_selection,
    resolve_sampling_channel_selection_replacing_default,
};

pub mod threshold_multiplier;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum OrientationProfileMode {
    #[default]
    Summed,
    PerOrientation,
}

impl OrientationProfileMode {
    pub fn profiles_per_orientation(self) -> bool {
        matches!(self, Self::PerOrientation)
    }
}

#[derive(Debug, Clone)]
pub struct MomentumSpaceEvaluationInput {
    pub loop_momenta: Vec<ThreeMomentum<F<f64>>>,
    pub integrator_weight: F<f64>,
    pub graph_id: Option<usize>,
    pub group_id: Option<GroupId>,
    pub orientation: Option<usize>,
    pub channel_id: Option<SamplingChannelId>,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
#[enum_dispatch(HasIntegrand)]
pub enum ProcessIntegrand {
    Amplitude(amplitude::AmplitudeIntegrand),
    CrossSection(cross_section::CrossSectionIntegrand),
}

fn discrete_sampling_type_name(sampling_type: &DiscreteGraphSamplingType) -> &'static str {
    match sampling_type {
        DiscreteGraphSamplingType::Default(_) => "default",
        DiscreteGraphSamplingType::MultiChanneling(_) => "multi_channeling",
        DiscreteGraphSamplingType::TropicalSampling(_) => "tropical",
        DiscreteGraphSamplingType::SamplingMultiChanneling(_) => "sampling_multi_channeling",
    }
}

pub(crate) fn discrete_sampling_depth_for_settings(
    settings: &DiscreteGraphSamplingSettings,
) -> usize {
    let orientation_depth = usize::from(settings.sample_orientations);
    match &settings.sampling_type {
        DiscreteGraphSamplingType::SamplingMultiChanneling(_) => 2 + orientation_depth,
        _ => 1 + orientation_depth,
    }
}

fn invalid_discrete_sampling_depth_error(
    settings: &DiscreteGraphSamplingSettings,
    actual_depth: usize,
) -> eyre::Report {
    let mut axes = vec!["graph group"];
    if settings.sample_orientations {
        axes.push("orientation");
    }
    if matches!(
        settings.sampling_type,
        DiscreteGraphSamplingType::SamplingMultiChanneling(_)
    ) {
        axes.push("channel");
    }

    eyre!(
        "This integrand uses discrete graph sampling (sample_orientations = {}, sampling_type = {}), so x-space evaluation requires {} discrete dimensions [{}], but got {}.",
        settings.sample_orientations,
        discrete_sampling_type_name(&settings.sampling_type),
        axes.len(),
        axes.join(", "),
        actual_depth
    )
}

pub(crate) fn resolve_discrete_selection_for_sampling(
    sampling: &SamplingSettings,
    discrete_dimensions: &[usize],
    group_count: usize,
    mut orientation_count_for_group: impl FnMut(GroupId) -> Option<usize>,
    mut channel_count_for_group: impl FnMut(GroupId) -> Result<Option<usize>>,
) -> Result<(Option<GroupId>, Option<usize>, Option<SamplingChannelId>)> {
    match sampling {
        SamplingSettings::Default(_) | SamplingSettings::MultiChanneling(_) => {
            if !discrete_dimensions.is_empty() {
                return Err(eyre!(
                    "This integrand does not use discrete graph sampling; expected no discrete dimensions, got {:?}.",
                    discrete_dimensions
                ));
            }
            Ok((None, None, None))
        }
        SamplingSettings::DiscreteGraphs(settings) => {
            let expected_depth = discrete_sampling_depth_for_settings(settings);
            if discrete_dimensions.len() != expected_depth {
                return Err(invalid_discrete_sampling_depth_error(
                    settings,
                    discrete_dimensions.len(),
                ));
            }

            let group_id = GroupId(discrete_dimensions[0]);
            if group_id.0 >= group_count {
                return Err(eyre!(
                    "Discrete graph group index {} is out of range; the integrand has {} groups.",
                    group_id.0,
                    group_count
                ));
            }

            let orientation = if settings.sample_orientations {
                let orientation = discrete_dimensions[1];
                let orientation_count = orientation_count_for_group(group_id).ok_or_else(|| {
                    eyre!(
                        "Could not determine orientation count for group {}.",
                        group_id.0
                    )
                })?;
                if orientation >= orientation_count {
                    return Err(eyre!(
                        "Orientation {} is out of range for graph group {}; the group has {} orientations.",
                        orientation,
                        group_id.0,
                        orientation_count
                    ));
                }
                Some(orientation)
            } else {
                None
            };

            let channel = match &settings.sampling_type {
                DiscreteGraphSamplingType::SamplingMultiChanneling(_) => {
                    let channel_index = *discrete_dimensions.last().expect("validated depth");
                    let channel_count = channel_count_for_group(group_id)?.ok_or_else(|| {
                        eyre!(
                            "Could not determine channel count for group {}.",
                            group_id.0
                        )
                    })?;
                    if channel_index >= channel_count {
                        return Err(eyre!(
                            "Channel {} is out of range for graph group {}; the group has {} channels.",
                            channel_index,
                            group_id.0,
                            channel_count
                        ));
                    }
                    Some(SamplingChannelId::from(channel_index))
                }
                _ => None,
            };

            Ok((Some(group_id), orientation, channel))
        }
    }
}

impl ProcessIntegrand {
    pub fn clone_with_selected_graph_groups(&self, graph_names: &[String]) -> Result<Self> {
        if graph_names.is_empty() {
            return Ok(self.clone());
        }
        if !matches!(
            self.get_settings().sampling,
            SamplingSettings::DiscreteGraphs(_)
        ) {
            return Err(eyre!(
                "Runtime graph-group selection requires graphs = 'monte_carlo'."
            ));
        }
        let mut unique_names = BTreeSet::new();
        for graph_name in graph_names {
            if !unique_names.insert(graph_name) {
                return Err(eyre!(
                    "Runtime graph-group selection contains duplicate graph name '{}'.",
                    graph_name
                ));
            }
        }

        let selection = GraphGroupSelectionSpec::from_master_graph_names(graph_names.to_vec());
        let mut selected = match self {
            Self::Amplitude(integrand) => {
                let plan = selection.plan(&integrand.data.graph_group_structure, |graph_id| {
                    integrand
                        .data
                        .graph_terms
                        .get(graph_id)
                        .map(|term| &term.graph)
                })?;
                Self::Amplitude(integrand.clone_with_graph_group_selection(&plan)?)
            }
            Self::CrossSection(integrand) => {
                let plan = selection.plan(&integrand.data.graph_group_structure, |graph_id| {
                    integrand
                        .data
                        .graph_terms
                        .get(graph_id)
                        .map(|term| &term.graph)
                })?;
                Self::CrossSection(integrand.clone_with_graph_group_selection(&plan)?)
            }
        };
        let selected_names = selected
            .graph_group_master_names()
            .into_iter()
            .map(str::to_string)
            .collect();
        let SamplingSettings::DiscreteGraphs(sampling) = &mut selected.get_mut_settings().sampling
        else {
            unreachable!("validated graph sampling before constructing the selected view")
        };
        sampling.graph_names = selected_names;
        Ok(selected)
    }

    pub fn resume_fingerprint(&self) -> Result<String> {
        let mut bytes = match self {
            Self::Amplitude(integrand) => {
                bincode::encode_to_vec(&integrand.data, bincode::config::standard())
            }
            Self::CrossSection(integrand) => {
                bincode::encode_to_vec(&integrand.data, bincode::config::standard())
            }
        }
        .map_err(|err| eyre!("Could not serialize integrand fingerprint payload: {err}"))?;

        let mut hash = 0xcbf29ce484222325u64;
        for byte in self
            .variant_tag()
            .as_bytes()
            .iter()
            .copied()
            .chain(bytes.drain(..))
        {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x100000001b3);
        }

        Ok(format!("{hash:016x}"))
    }

    fn variant_tag(&self) -> &'static str {
        match self {
            Self::Amplitude(_) => "amplitude",
            Self::CrossSection(_) => "cross_section",
        }
    }

    pub fn kind_name(&self) -> &'static str {
        self.variant_tag()
    }

    pub fn export_standalone(
        &self,
        path: impl AsRef<Path>,
        settings: &StandaloneExportSettings,
    ) -> Result<()> {
        match self {
            Self::Amplitude(a) => a.export_standalone(path, settings),
            Self::CrossSection(a) => a.export_standalone(path, settings),
        }
    }

    pub fn warm_up(&mut self, model: &Model) -> Result<()> {
        let settings = self.get_mut_settings();
        if matches!(
            &settings.sampling,
            SamplingSettings::DiscreteGraphs(settings) if settings.sample_orientations
        ) {
            if !matches!(
                settings.general.evaluator_method,
                evaluators::EvaluatorMethod::SingleParametric
            ) {
                return Err(eyre!(
                    "Monte Carlo sampling over orientations requires general.evaluator_method=SingleParametric; got {:?}.",
                    settings.general.evaluator_method
                ));
            }
            warn!("Monte Carlo sampling over orientations is using the SingleParametric evaluator");
        }

        match self {
            Self::Amplitude(a) => a.warm_up(model),
            Self::CrossSection(a) => a.warm_up(model),
        }
    }

    pub fn frozen_compilation(&self) -> &crate::settings::global::FrozenCompilationMode {
        match self {
            Self::Amplitude(a) => a.frozen_compilation(),
            Self::CrossSection(a) => a.frozen_compilation(),
        }
    }

    pub fn active_f64_backend(&self) -> ActiveF64Backend {
        match self {
            Self::Amplitude(a) => a.active_f64_backend(),
            Self::CrossSection(a) => a.active_f64_backend(),
        }
    }

    pub(crate) fn activate_runtime_backends_after_load(
        &mut self,
        allow_symjit_fallback: bool,
    ) -> Result<Option<String>> {
        match self {
            Self::Amplitude(a) => a.activate_runtime_backends_after_load(allow_symjit_fallback),
            Self::CrossSection(a) => a.activate_runtime_backends_after_load(allow_symjit_fallback),
        }
    }

    pub(crate) fn save(&self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let path = path.as_ref().join("integrand");

        let r = fs::create_dir_all(&path).with_context(|| {
            format!(
                "Trying to create directory to save amplitude {}",
                path.display()
            )
        });
        if override_existing {
            r?;
        }
        match self {
            ProcessIntegrand::Amplitude(integrand) => integrand.save(path, override_existing),
            ProcessIntegrand::CrossSection(integrand) => integrand.save(path, override_existing),
        }
    }

    pub(crate) fn compile(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        thread_pool: &rayon::ThreadPool,
    ) -> Result<Vec<(String, std::time::Duration)>> {
        let path = path.as_ref().join("integrand");

        let r = fs::create_dir_all(&path).with_context(|| {
            format!(
                "Trying to create directory to save amplitude {}",
                path.display()
            )
        });
        if override_existing {
            r?;
        }
        match self {
            ProcessIntegrand::Amplitude(integrand) => {
                integrand.compile(path, override_existing, thread_pool)
            }
            ProcessIntegrand::CrossSection(integrand) => {
                integrand.compile(path, override_existing, thread_pool)
            }
        }
    }

    pub fn get_settings(&self) -> &RuntimeSettings {
        match self {
            ProcessIntegrand::Amplitude(integrand) => &integrand.settings,
            ProcessIntegrand::CrossSection(integrand) => &integrand.settings,
        }
    }

    /// Invalidate runtime caches before exposing mutable settings. Call
    /// `warm_up` again before evaluating samples through sampling channels.
    pub fn get_mut_settings(&mut self) -> &mut RuntimeSettings {
        match self {
            ProcessIntegrand::Amplitude(integrand) => {
                integrand.invalidate_runtime_caches();
                &mut integrand.settings
            }
            ProcessIntegrand::CrossSection(integrand) => {
                integrand.invalidate_runtime_caches();
                &mut integrand.settings
            }
        }
    }

    pub fn graph_count(&self) -> usize {
        match self {
            ProcessIntegrand::Amplitude(integrand) => integrand.data.graph_terms.len(),
            ProcessIntegrand::CrossSection(integrand) => integrand.data.graph_terms.len(),
        }
    }

    pub fn graph_group_count(&self) -> usize {
        match self {
            ProcessIntegrand::Amplitude(integrand) => integrand.data.graph_group_structure.len(),
            ProcessIntegrand::CrossSection(integrand) => integrand.data.graph_group_structure.len(),
        }
    }

    pub fn graph_group_master_names(&self) -> Vec<&str> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => integrand
                .data
                .graph_group_structure
                .iter()
                .map(|group| {
                    integrand.data.graph_terms[group.master()]
                        .graph
                        .name
                        .as_str()
                })
                .collect(),
            ProcessIntegrand::CrossSection(integrand) => integrand
                .data
                .graph_group_structure
                .iter()
                .map(|group| {
                    integrand.data.graph_terms[group.master()]
                        .graph
                        .name
                        .as_str()
                })
                .collect(),
        }
    }

    pub fn find_graph_id_by_name(&self, graph_name: &str) -> Option<usize> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => integrand
                .data
                .graph_terms
                .iter()
                .position(|term| term.graph.name == graph_name),
            ProcessIntegrand::CrossSection(integrand) => integrand
                .data
                .graph_terms
                .iter()
                .position(|term| term.graph.name == graph_name),
        }
    }

    pub fn graph_name_by_id(&self, graph_id: usize) -> Option<&str> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => integrand
                .data
                .graph_terms
                .get(graph_id)
                .map(|term| term.graph.name.as_str()),
            ProcessIntegrand::CrossSection(integrand) => integrand
                .data
                .graph_terms
                .get(graph_id)
                .map(|term| term.graph.name.as_str()),
        }
    }

    pub fn graph_group_id_by_graph_id(&self, graph_id: usize) -> Option<usize> {
        self.find_group_id_containing_graph(graph_id)
            .map(usize::from)
    }

    pub fn cut_edge_ids(&self, graph_id: usize, cut_id: usize) -> Option<Vec<usize>> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => {
                (cut_id == 0 && graph_id < integrand.data.graph_terms.len()).then(Vec::new)
            }
            ProcessIntegrand::CrossSection(integrand) => {
                let graph_term = integrand.data.graph_terms.get(graph_id)?;
                let cut = graph_term.cuts.get(crate::processes::CutId::from(cut_id))?;
                Some(
                    graph_term
                        .graph
                        .underlying
                        .iter_edges_of(&cut.cut)
                        .map(|(_, edge_id, _)| edge_id.0)
                        .sorted()
                        .collect(),
                )
            }
        }
    }

    pub fn lmb_sample_id_for_channel(
        &self,
        graph_id: usize,
        sampling_channel_id: usize,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<Option<usize>> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => {
                let Some(graph_term) = integrand.data.graph_terms.get(graph_id) else {
                    return Ok(None);
                };
                Ok(graph_term
                    .multi_channeling_setup
                    .sampling_channel_lmb_id(
                        SamplingChannelId::from(sampling_channel_id),
                        &graph_term.multi_channeling_setup.graph.name,
                        parameterization_settings,
                    )?
                    .map(usize::from))
            }
            ProcessIntegrand::CrossSection(integrand) => {
                let Some(graph_term) = integrand.data.graph_terms.get(graph_id) else {
                    return Ok(None);
                };
                Ok(graph_term
                    .multi_channeling_setup
                    .sampling_channel_lmb_id(
                        SamplingChannelId::from(sampling_channel_id),
                        &graph_term.multi_channeling_setup.graph.name,
                        parameterization_settings,
                    )?
                    .map(usize::from))
            }
        }
    }

    fn find_group_id_containing_graph(&self, graph_id: usize) -> Option<GroupId> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => integrand
                .data
                .graph_group_structure
                .iter_enumerated()
                .find_map(|(group_id, group)| {
                    group.into_iter().contains(&graph_id).then_some(group_id)
                }),
            ProcessIntegrand::CrossSection(integrand) => integrand
                .data
                .graph_group_structure
                .iter_enumerated()
                .find_map(|(group_id, group)| {
                    group.into_iter().contains(&graph_id).then_some(group_id)
                }),
        }
    }

    pub fn resolve_group_id_by_master_name(&self, graph_name: &str) -> Result<GroupId> {
        match self.find_graph_id_by_name(graph_name) {
            Some(graph_id) => {
                let group_id = self
                    .find_group_id_containing_graph(graph_id)
                    .ok_or_else(|| {
                        eyre!("Could not find graph group for graph '{}'.", graph_name)
                    })?;
                let master_graph_name = match self {
                    ProcessIntegrand::Amplitude(integrand) => {
                        &integrand.data.graph_terms
                            [integrand.data.graph_group_structure[group_id].master()]
                        .graph
                        .name
                    }
                    ProcessIntegrand::CrossSection(integrand) => {
                        &integrand.data.graph_terms
                            [integrand.data.graph_group_structure[group_id].master()]
                        .graph
                        .name
                    }
                };
                if master_graph_name != graph_name {
                    return Err(eyre!(
                        "Graph '{}' is not the master graph of its group; use '{}' instead.",
                        graph_name,
                        master_graph_name
                    ));
                }
                Ok(group_id)
            }
            None => Err(eyre!(
                "Unknown graph '{}' in momentum-space evaluation.",
                graph_name
            )),
        }
    }

    pub fn resolve_discrete_selection(
        &self,
        discrete_dimensions: &[usize],
    ) -> Result<(Option<GroupId>, Option<usize>, Option<SamplingChannelId>)> {
        let group_count = match self {
            ProcessIntegrand::Amplitude(integrand) => integrand.data.graph_group_structure.len(),
            ProcessIntegrand::CrossSection(integrand) => integrand.data.graph_group_structure.len(),
        };

        resolve_discrete_selection_for_sampling(
            &self.get_settings().sampling,
            discrete_dimensions,
            group_count,
            |group_id| self.group_orientation_count(group_id),
            |group_id| Ok(self.group_channel_count(group_id)),
        )
    }

    pub fn expected_x_space_dimension(&self, discrete_dimensions: &[usize]) -> Result<usize> {
        let settings = self.get_settings();
        let (group_id, _, _) = self.resolve_discrete_selection(discrete_dimensions)?;
        if matches!(
            &settings.sampling,
            SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
                sampling_type: DiscreteGraphSamplingType::TropicalSampling(_),
                ..
            })
        ) {
            let group_id = group_id.ok_or_else(|| {
                eyre!("Tropical sampling requires a discrete graph-group selection.")
            })?;
            let (loop_number, n_edges) = match self {
                ProcessIntegrand::Amplitude(integrand) => {
                    let master_graph = &integrand.data.graph_terms
                        [integrand.data.graph_group_structure[group_id].master()];
                    (
                        master_graph.get_graph().get_loop_number(),
                        master_graph.get_graph().iter_loop_edges().count(),
                    )
                }
                ProcessIntegrand::CrossSection(integrand) => {
                    let master_graph = &integrand.data.graph_terms
                        [integrand.data.graph_group_structure[group_id].master()];
                    (
                        master_graph.get_graph().get_loop_number(),
                        master_graph.get_graph().iter_loop_edges().count(),
                    )
                }
            };
            return Ok(get_n_dim_for_n_loop_momenta(
                &settings.sampling,
                loop_number,
                Some(n_edges),
            ));
        }

        let loop_number = match self {
            ProcessIntegrand::Amplitude(integrand) => {
                integrand.data.graph_terms[0].graph.get_loop_number()
            }
            ProcessIntegrand::CrossSection(integrand) => {
                integrand.data.graph_terms[0].graph.get_loop_number()
            }
        };
        Ok(get_n_dim_for_n_loop_momenta(
            &settings.sampling,
            loop_number,
            None,
        ))
    }

    pub fn discrete_sampling_depth(&self) -> usize {
        match &self.get_settings().sampling {
            SamplingSettings::Default(_) | SamplingSettings::MultiChanneling(_) => 0,
            SamplingSettings::DiscreteGraphs(settings) => {
                discrete_sampling_depth_for_settings(settings)
            }
        }
    }

    pub fn group_orientation_count(&self, group_id: GroupId) -> Option<usize> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => Some(
                integrand.data.graph_terms[integrand.data.graph_group_structure[group_id].master()]
                    .get_num_orientations(),
            ),
            ProcessIntegrand::CrossSection(integrand) => Some(
                integrand.data.graph_terms[integrand.data.graph_group_structure[group_id].master()]
                    .get_num_orientations(),
            ),
        }
    }

    /// Return the canonical sampling-channel IDs for a graph group.
    ///
    /// The IDs come directly from the master graph's single catalogue. This
    /// helper is intended for saved-state acceptance callers constructing
    /// explicit `(group, orientation, channel)` selections; callers must pass
    /// the same parameterization settings used by the loaded state.
    pub fn group_sampling_channel_ids(
        &self,
        group_id: GroupId,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<Vec<SamplingChannelId>> {
        let group = match self {
            ProcessIntegrand::Amplitude(integrand) => {
                integrand.data.graph_group_structure.get(group_id)
            }
            ProcessIntegrand::CrossSection(integrand) => {
                integrand.data.graph_group_structure.get(group_id)
            }
        }
        .ok_or_else(|| eyre!("Unknown graph group {}.", group_id.0))?;
        let master = group.master();
        match self {
            ProcessIntegrand::Amplitude(integrand) => {
                integrand.data.graph_terms[master].sampling_channel_ids(parameterization_settings)
            }
            ProcessIntegrand::CrossSection(integrand) => {
                integrand.data.graph_terms[master].sampling_channel_ids(parameterization_settings)
            }
        }
    }

    pub fn group_channel_count(&self, group_id: GroupId) -> Option<usize> {
        let parameterization_settings = self
            .get_settings()
            .sampling
            .get_parameterization_settings()
            .unwrap_or_default();
        self.group_sampling_channel_ids(group_id, &parameterization_settings)
            .ok()
            .map(|channel_ids| channel_ids.len())
    }

    pub fn graph_orientation_count(&self, graph_id: usize) -> Option<usize> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => integrand
                .data
                .graph_terms
                .get(graph_id)
                .map(|term| term.get_num_orientations()),
            ProcessIntegrand::CrossSection(integrand) => integrand
                .data
                .graph_terms
                .get(graph_id)
                .map(|term| term.get_num_orientations()),
        }
    }

    pub fn evaluate_momentum_configuration(
        &mut self,
        model: &Model,
        input: &MomentumSpaceEvaluationInput,
        use_arb_prec: bool,
    ) -> Result<EvaluationResult> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => evaluate_momentum_configuration(
                integrand,
                model,
                input,
                input.integrator_weight,
                use_arb_prec,
                Complex::new_zero(),
            ),
            ProcessIntegrand::CrossSection(integrand) => evaluate_momentum_configuration(
                integrand,
                model,
                input,
                input.integrator_weight,
                use_arb_prec,
                Complex::new_zero(),
            ),
        }
    }

    pub fn evaluate_sample_precise(
        &mut self,
        sample: &Sample<F<f64>>,
        model: &Model,
        wgt: F<f64>,
        use_arb_prec: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<PreciseEvaluationResult> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => {
                evaluate_sample_precise(integrand, model, sample, wgt, use_arb_prec, max_eval)
            }
            ProcessIntegrand::CrossSection(integrand) => {
                evaluate_sample_precise(integrand, model, sample, wgt, use_arb_prec, max_eval)
            }
        }
    }

    /// Evaluate a normalized reference function through the real process
    /// parameterization, retaining its Jacobian and integrator weight.
    ///
    /// This is intentionally a test-only style overlay: graph routing,
    /// orientation/channel selection and all sampling maps are still exercised,
    /// while the physical graph evaluation is replaced by `reference`.
    pub fn evaluate_reference_sample(
        &mut self,
        sample: &Sample<F<f64>>,
        reference: &GaussianReferenceFunction,
    ) -> Result<EvaluationResult> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => {
                Ok(evaluate_reference_sample(integrand, sample, reference)?.evaluation)
            }
            ProcessIntegrand::CrossSection(integrand) => {
                Ok(evaluate_reference_sample(integrand, sample, reference)?.evaluation)
            }
        }
    }

    /// Evaluate a normalized reference function from unit-cube coordinates.
    ///
    /// This convenience entry point is intended for saved-state acceptance
    /// harnesses: callers can load a process normally, generate deterministic
    /// coordinates (for example Halton points), and exercise the complete
    /// process parameterization without constructing Symbolica `Sample`
    /// values themselves. The coordinate batch is interpreted as an equally
    /// weighted quadrature rule, so each sample has unit inverse-density weight
    /// and the report averages once over the `N` draws; the
    /// canonical sampling channel selected by
    /// the loaded settings remains responsible for its map and Jacobian.
    pub fn evaluate_reference_coordinates(
        &mut self,
        coordinates: &[Vec<f64>],
        reference: &GaussianReferenceFunction,
    ) -> Result<ReferenceSamplingReport> {
        if coordinates.is_empty() {
            return Err(eyre!(
                "reference acceptance coordinate batch needs at least one sample"
            ));
        }
        let samples = coordinates
            .iter()
            .enumerate()
            .map(|(sample_index, coordinate)| {
                if coordinate.iter().any(|value| !value.is_finite()) {
                    return Err(eyre!(
                        "reference acceptance sample {sample_index} contains a non-finite coordinate"
                    ));
                }
                Ok(Sample::Continuous(
                    F(1.0),
                    coordinate.iter().copied().map(F).collect(),
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        self.evaluate_reference_samples(&samples, reference)
    }

    /// Evaluate deterministic coordinates under one explicit discrete
    /// selection path. `discrete_indices` is ordered from the outermost
    /// `Sample::Discrete` selector to the innermost one, matching the
    /// canonical graph/orientation/channel order used by the process sampler.
    /// The selection itself is never re-enumerated or filtered here; callers
    /// must supply IDs obtained from the loaded process' canonical catalogue.
    /// Unit discrete weights make this a partition contribution, which need
    /// not integrate to one. A complete acceptance estimator must sum the
    /// per-draw contributions or use MC samples with reciprocal selection
    /// probabilities. Do not add report errors for channels sharing the same
    /// coordinates: their contributions must be combined before squaring.
    pub fn evaluate_reference_discrete_coordinates(
        &mut self,
        discrete_indices: &[usize],
        coordinates: &[Vec<f64>],
        reference: &GaussianReferenceFunction,
    ) -> Result<ReferenceSamplingReport> {
        if coordinates.is_empty() {
            return Err(eyre!(
                "reference acceptance coordinate batch needs at least one sample"
            ));
        }
        let samples = coordinates
            .iter()
            .enumerate()
            .map(|(sample_index, coordinate)| {
                if coordinate.iter().any(|value| !value.is_finite()) {
                    return Err(eyre!(
                        "reference acceptance sample {sample_index} contains a non-finite coordinate"
                    ));
                }
                let mut sample = Sample::Continuous(
                    F(1.0),
                    coordinate.iter().copied().map(F).collect(),
                );
                for index in discrete_indices.iter().rev() {
                    sample = Sample::Discrete(F(1.0), *index, Some(Box::new(sample)));
                }
                Ok(sample)
            })
            .collect::<Result<Vec<_>>>()?;
        self.evaluate_reference_samples(&samples, reference)
    }

    /// Evaluate a batch of samples with a normalized reference function while
    /// retaining the process maps, Jacobians and sampling-grid weights.
    pub fn evaluate_reference_samples(
        &mut self,
        samples: &[Sample<F<f64>>],
        reference: &GaussianReferenceFunction,
    ) -> Result<ReferenceSamplingReport> {
        if samples.is_empty() {
            return Err(eyre!("reference acceptance needs at least one sample"));
        }
        let evaluations = samples
            .iter()
            .map(|sample| self.evaluate_reference_sample_detailed(sample, reference))
            .collect::<Result<Vec<_>>>()?;
        let report = ReferenceSamplingReport::from_evaluations(evaluations, reference);
        if report.finite_sample_count != report.sample_count
            || !report.normalization.is_finite()
            || !report.normalization_stderr.is_finite()
            || !report.second_moment.is_finite()
            || !report.second_moment_stderr.is_finite()
        {
            return Err(eyre!(
                "reference acceptance has invalid mapped contributions or statistics: {} finite draws out of {}",
                report.finite_sample_count,
                report.sample_count
            ));
        }
        Ok(report)
    }

    /// Detailed variant used by acceptance harnesses that also check raw
    /// momentum moments.  The public scalar method above remains convenient
    /// for callers interested only in the mapped value.
    pub fn evaluate_reference_sample_detailed(
        &mut self,
        sample: &Sample<F<f64>>,
        reference: &GaussianReferenceFunction,
    ) -> Result<ReferenceSampleEvaluation> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => {
                evaluate_reference_sample(integrand, sample, reference)
            }
            ProcessIntegrand::CrossSection(integrand) => {
                evaluate_reference_sample(integrand, sample, reference)
            }
        }
    }

    pub fn evaluate_momentum_configuration_precise(
        &mut self,
        model: &Model,
        input: &MomentumSpaceEvaluationInput,
        use_arb_prec: bool,
    ) -> Result<PreciseEvaluationResult> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => evaluate_momentum_configuration_precise(
                integrand,
                model,
                input,
                input.integrator_weight,
                use_arb_prec,
                Complex::new_zero(),
            ),
            ProcessIntegrand::CrossSection(integrand) => evaluate_momentum_configuration_precise(
                integrand,
                model,
                input,
                input.integrator_weight,
                use_arb_prec,
                Complex::new_zero(),
            ),
        }
    }

    pub fn evaluate_samples_raw(
        &mut self,
        target: EvaluationTarget<'_>,
        samples: &[Sample<F<f64>>],
        _iter: usize,
        use_arb_prec: bool,
        stop_on_interrupt: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<RawBatchEvaluationResult> {
        let mut results = Vec::with_capacity(samples.len());
        for sample in samples {
            if stop_on_interrupt && crate::is_interrupted() {
                break;
            }
            macro_rules! evaluate {
                ($integrand:expr) => {
                    evaluate_from_source_precise(
                        $integrand,
                        target,
                        EvaluationSource::XSpace(sample),
                        sample.get_weight(),
                        use_arb_prec,
                        max_eval,
                    )
                };
            }
            let precise = match self {
                ProcessIntegrand::Amplitude(integrand) => evaluate!(integrand),
                ProcessIntegrand::CrossSection(integrand) => evaluate!(integrand),
            }?;
            let mut result = match target {
                EvaluationTarget::Reference(reference) => reference.integration_report(precise)?,
                EvaluationTarget::Physical(_) => precise.try_into_f64()?,
                #[cfg(test)]
                EvaluationTarget::SamplingLaw(_) => precise.try_into_f64()?,
            };

            self.process_evaluation_result(&result);
            maybe_discard_generated_events_in_result(self.get_settings(), &mut result);
            results.push(result);
            if stop_on_interrupt && crate::is_interrupted() {
                break;
            }
        }

        Ok(RawBatchEvaluationResult {
            statistics: StatisticsCounter::from_evaluation_results(&results),
            samples: results,
        })
    }

    pub fn evaluate_momentum_configurations_raw(
        &mut self,
        model: &Model,
        inputs: &[MomentumSpaceEvaluationInput],
        use_arb_prec: bool,
    ) -> Result<RawBatchEvaluationResult> {
        let mut results = Vec::with_capacity(inputs.len());
        for input in inputs {
            let mut result = match self {
                ProcessIntegrand::Amplitude(integrand) => evaluate_momentum_configuration(
                    integrand,
                    model,
                    input,
                    input.integrator_weight,
                    use_arb_prec,
                    Complex::new_zero(),
                ),
                ProcessIntegrand::CrossSection(integrand) => evaluate_momentum_configuration(
                    integrand,
                    model,
                    input,
                    input.integrator_weight,
                    use_arb_prec,
                    Complex::new_zero(),
                ),
            }?;

            self.process_evaluation_result(&result);
            maybe_discard_generated_events_in_result(self.get_settings(), &mut result);
            results.push(result);
        }

        Ok(RawBatchEvaluationResult {
            statistics: StatisticsCounter::from_evaluation_results(&results),
            samples: results,
        })
    }

    pub fn evaluate_samples_precise_raw(
        &mut self,
        model: &Model,
        samples: &[Sample<F<f64>>],
        use_arb_prec: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<RawPreciseBatchEvaluationResult> {
        let mut results = Vec::with_capacity(samples.len());
        for sample in samples {
            results.push(self.evaluate_sample_precise(
                sample,
                model,
                sample.get_weight(),
                use_arb_prec,
                max_eval,
            )?);
        }

        Ok(RawPreciseBatchEvaluationResult { samples: results })
    }

    pub fn evaluate_momentum_configurations_precise_raw(
        &mut self,
        model: &Model,
        inputs: &[MomentumSpaceEvaluationInput],
        use_arb_prec: bool,
    ) -> Result<RawPreciseBatchEvaluationResult> {
        let mut results = Vec::with_capacity(inputs.len());
        for input in inputs {
            results.push(self.evaluate_momentum_configuration_precise(
                model,
                input,
                use_arb_prec,
            )?);
        }

        Ok(RawPreciseBatchEvaluationResult { samples: results })
    }

    pub fn process_evaluation_result(&mut self, result: &EvaluationResult) {
        match self {
            ProcessIntegrand::Amplitude(integrand) => {
                process_evaluation_result_runtime(integrand, result)
            }
            ProcessIntegrand::CrossSection(integrand) => {
                process_evaluation_result_runtime(integrand, result)
            }
        }
    }

    pub fn merge_event_processing_runtime(&mut self, other: &mut Self) -> Result<()> {
        match (self, other) {
            (ProcessIntegrand::Amplitude(lhs), ProcessIntegrand::Amplitude(rhs)) => {
                merge_event_processing_runtime(lhs, rhs)
            }
            (ProcessIntegrand::CrossSection(lhs), ProcessIntegrand::CrossSection(rhs)) => {
                merge_event_processing_runtime(lhs, rhs)
            }
            _ => Err(eyre!(
                "Cannot merge event-processing runtime for incompatible process integrands."
            )),
        }
    }

    pub fn update_event_processing_runtime(&mut self, iter: usize) {
        match self {
            ProcessIntegrand::Amplitude(integrand) => {
                update_event_processing_runtime(integrand, iter)
            }
            ProcessIntegrand::CrossSection(integrand) => {
                update_event_processing_runtime(integrand, iter)
            }
        }
    }

    pub fn observable_accumulator_bundle(&self) -> Option<ObservableAccumulatorBundle> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => observable_accumulator_bundle(integrand),
            ProcessIntegrand::CrossSection(integrand) => observable_accumulator_bundle(integrand),
        }
    }

    pub fn has_observables(&self) -> bool {
        match self {
            ProcessIntegrand::Amplitude(integrand) => integrand
                .event_processing_runtime()
                .is_some_and(EventProcessingRuntime::has_observables),
            ProcessIntegrand::CrossSection(integrand) => integrand
                .event_processing_runtime()
                .is_some_and(EventProcessingRuntime::has_observables),
        }
    }

    pub fn observable_snapshot_bundle(&self) -> Option<ObservableSnapshotBundle> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => observable_snapshot_bundle(integrand),
            ProcessIntegrand::CrossSection(integrand) => observable_snapshot_bundle(integrand),
        }
    }

    pub fn build_observable_snapshots_for_result(
        &self,
        result: &EvaluationResult,
    ) -> Option<ObservableSnapshotBundle> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => {
                build_observable_snapshots_for_result(integrand, result)
            }
            ProcessIntegrand::CrossSection(integrand) => {
                build_observable_snapshots_for_result(integrand, result)
            }
        }
    }

    pub fn build_observable_snapshots_for_precise_result(
        &self,
        result: &PreciseEvaluationResult,
    ) -> Option<ObservableSnapshotBundle> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => match result {
                PreciseEvaluationResult::Double(result) => {
                    build_observable_snapshots_for_event_groups(integrand, &result.event_groups)
                }
                PreciseEvaluationResult::Quad(result) => {
                    build_observable_snapshots_for_event_groups(integrand, &result.event_groups)
                }
                PreciseEvaluationResult::Arb(result) => {
                    build_observable_snapshots_for_event_groups(integrand, &result.event_groups)
                }
            },
            ProcessIntegrand::CrossSection(integrand) => match result {
                PreciseEvaluationResult::Double(result) => {
                    build_observable_snapshots_for_event_groups(integrand, &result.event_groups)
                }
                PreciseEvaluationResult::Quad(result) => {
                    build_observable_snapshots_for_event_groups(integrand, &result.event_groups)
                }
                PreciseEvaluationResult::Arb(result) => {
                    build_observable_snapshots_for_event_groups(integrand, &result.event_groups)
                }
            },
        }
    }

    pub fn write_observable_snapshots(
        &self,
        path: impl AsRef<Path>,
        format: ObservableFileFormat,
    ) -> Result<()> {
        let Some(bundle) = self.observable_snapshot_bundle() else {
            return Ok(());
        };

        write_observable_snapshot_bundle(&bundle, path.as_ref(), format)
    }

    pub fn restore_observable_snapshot_bundle(
        &mut self,
        bundle: &ObservableSnapshotBundle,
    ) -> Result<()> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => {
                restore_observable_snapshot_bundle(integrand, bundle)
            }
            ProcessIntegrand::CrossSection(integrand) => {
                restore_observable_snapshot_bundle(integrand, bundle)
            }
        }
    }
}

fn format_orientation_label(signature: &EdgeVec<Orientation>) -> String {
    signature
        .iter()
        .map(|(_, orientation)| match *orientation {
            Orientation::Default => '+',
            Orientation::Reversed => '-',
            Orientation::Undirected => '0',
        })
        .collect()
}

pub(crate) fn resolve_visible_orientation_id(
    orientation_filter: &SubSet<OrientationID>,
    visible_orientation_id: usize,
) -> Option<OrientationID> {
    if orientation_filter.is_full() {
        Some(OrientationID::from(visible_orientation_id))
    } else {
        orientation_filter
            .included_iter()
            .nth(visible_orientation_id)
    }
}

pub(crate) fn filtered_orientation_count(
    orientation_filter: &SubSet<OrientationID>,
    orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
) -> usize {
    if orientation_filter.is_full() {
        orientations.len()
    } else {
        orientation_filter.included_iter().count()
    }
}

pub(crate) fn orientation_labels_for_graph<I: ProcessIntegrandImpl>(
    integrand: &I,
    graph_id: usize,
) -> Result<Vec<String>> {
    if graph_id >= integrand.graph_count() {
        return Err(eyre!(
            "Unknown graph '{}' while resolving orientation labels.",
            graph_id
        ));
    }
    let graph = integrand.get_graph(graph_id);
    Ok((0..graph.get_num_orientations())
        .map(|orientation_id| {
            graph
                .orientation_label(orientation_id)
                .unwrap_or_else(|| format!("#{}", orientation_id))
        })
        .collect())
}

pub(crate) fn evaluate_profile_momentum_point<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    graph_id: usize,
    orientation: Option<usize>,
    loop_momenta: Vec<ThreeMomentum<F<f64>>>,
    use_arb_prec: bool,
) -> Result<EvaluationResult> {
    let input = MomentumSpaceEvaluationInput {
        loop_momenta,
        integrator_weight: F(1.0),
        graph_id: Some(graph_id),
        group_id: None,
        orientation,
        channel_id: None,
    };
    evaluate_momentum_configuration(
        integrand,
        model,
        &input,
        F(1.0),
        use_arb_prec,
        Complex::new_re(F(100.0 * integrand.get_settings().kinematics.e_cm)),
    )
}

pub(crate) fn evaluate_profile_momentum_point_precise<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    graph_id: usize,
    orientation: Option<usize>,
    loop_momenta: Vec<ThreeMomentum<F<f64>>>,
    use_arb_prec: bool,
) -> Result<PreciseEvaluationResult> {
    let input = MomentumSpaceEvaluationInput {
        loop_momenta,
        integrator_weight: F(1.0),
        graph_id: Some(graph_id),
        group_id: None,
        orientation,
        channel_id: None,
    };
    evaluate_momentum_configuration_precise(
        integrand,
        model,
        &input,
        F(1.0),
        use_arb_prec,
        Complex::new_re(F(100.0 * integrand.get_settings().kinematics.e_cm)),
    )
}

fn format_sampling_channel_label(edge_ids: &[usize]) -> String {
    let mut sorted = edge_ids.to_vec();
    sorted.sort_unstable();
    format!(
        "({})",
        sorted
            .iter()
            .map(|edge_id| edge_id.to_string())
            .collect_vec()
            .join(",")
    )
}

pub(crate) fn histogram_process_info_for_integrand<I: ProcessIntegrandImpl>(
    integrand: &I,
) -> Result<HistogramProcessInfo> {
    let parameterization_settings = integrand
        .get_settings()
        .sampling
        .get_parameterization_settings()
        .unwrap_or_default();
    let graph_names = (0..integrand.graph_count())
        .map(|graph_id| integrand.get_graph(graph_id).name())
        .collect_vec();
    let graph_to_group_id = (0..integrand.graph_count())
        .map(|graph_id| {
            integrand
                .graph_group_id_for_graph(graph_id)
                .unwrap_or_else(|| {
                    panic!(
                        "graph {} is missing a graph-group mapping for histogram process info",
                        graph_id
                    )
                })
        })
        .collect_vec();
    let graph_group_master_names = integrand
        .get_group_structure()
        .iter_enumerated()
        .map(|(group_id, _)| integrand.get_master_graph(group_id).name())
        .collect_vec();
    let orientation_labels_by_group = integrand
        .get_group_structure()
        .iter_enumerated()
        .map(|(group_id, _)| {
            let master = integrand.get_master_graph(group_id);
            (0..master.get_num_orientations())
                .map(|orientation_id| {
                    master
                        .orientation_label(orientation_id)
                        .unwrap_or_else(|| format!("#{}", orientation_id))
                })
                .collect_vec()
        })
        .collect_vec();
    let sampling_channel_labels_by_group = integrand
        .get_group_structure()
        .iter_enumerated()
        .map(|(group_id, _)| {
            let master = integrand.get_master_graph(group_id);
            let channel_ids = master.sampling_channel_ids(&parameterization_settings)?;
            channel_ids
                .into_iter()
                .map(|channel_id| {
                    Ok(master
                        .sampling_channel_label(channel_id, &parameterization_settings)?
                        .unwrap_or_else(|| format!("#{}", channel_id.index())))
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(HistogramProcessInfo {
        graph_names,
        graph_to_group_id,
        graph_group_master_names,
        orientation_labels_by_group,
        sampling_channel_labels_by_group,
    })
}

pub(crate) fn graph_to_group_id_for_group_structure(
    group_structure: &TiVec<GroupId, GraphGroup>,
) -> Vec<usize> {
    group_structure
        .iter_enumerated()
        .flat_map(|(group_id, group)| {
            group
                .into_iter()
                .map(move |graph_id| (graph_id, group_id.0))
        })
        .sorted_by_key(|(graph_id, _)| *graph_id)
        .map(|(_, group_id)| group_id)
        .collect_vec()
}

pub(crate) struct PreparedBufferedEvent<T: FloatLike> {
    pub(crate) buffered_event: Option<GenericEvent<T>>,
    pub(crate) selectors_pass: bool,
    pub(crate) event_processing_time: Duration,
    pub(crate) generated_event_count: usize,
    pub(crate) accepted_event_count: usize,
}

impl<T: FloatLike> Default for PreparedBufferedEvent<T> {
    fn default() -> Self {
        Self {
            buffered_event: None,
            selectors_pass: true,
            event_processing_time: Duration::ZERO,
            generated_event_count: 0,
            accepted_event_count: 0,
        }
    }
}

pub(crate) fn prepare_buffered_event<T: FloatLike>(
    settings: &RuntimeSettings,
    rotation: &Rotation,
    event_processing_runtime: Option<&mut EventProcessingRuntime>,
    build_event: impl FnOnce() -> Result<GenericEvent<T>>,
) -> Result<PreparedBufferedEvent<T>> {
    let needs_selector_events = event_processing_runtime
        .as_ref()
        .is_some_and(|runtime| runtime.has_selectors());
    let should_buffer_event = rotation.is_identity() && settings.should_buffer_generated_events();
    let should_build_event = if rotation.is_identity() {
        settings.should_generate_events()
    } else {
        needs_selector_events
    };

    if !should_build_event {
        return Ok(PreparedBufferedEvent::default());
    }

    let build_start = Instant::now();
    let mut event = build_event()?;
    let mut event_processing_time = build_start.elapsed();
    let generated_event_count = usize::from(rotation.is_identity());

    let selector_start = Instant::now();
    let selectors_pass = if let Some(runtime) = event_processing_runtime {
        if rotation.is_identity() {
            runtime.process_event(&mut event)
        } else {
            runtime.process_event_for_selectors(&mut event)
        }
    } else {
        true
    };
    event_processing_time += selector_start.elapsed();

    let buffered_event = if selectors_pass && should_buffer_event {
        Some(event)
    } else {
        None
    };

    // Count selector-accepted physical events even when they are processed only
    // transiently for selectors and are not retained in the returned event buffer.
    let accepted_event_count = usize::from(rotation.is_identity() && selectors_pass);

    Ok(PreparedBufferedEvent {
        accepted_event_count,
        buffered_event,
        selectors_pass,
        event_processing_time,
        generated_event_count,
    })
}

fn process_evaluation_result_runtime<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    result: &EvaluationResult,
) {
    if let Some(runtime) = integrand.event_processing_runtime_mut()
        && runtime.has_observables()
    {
        runtime.process_event_groups(&result.event_groups);
    }
}

fn maybe_discard_generated_events_in_result(
    settings: &RuntimeSettings,
    result: &mut EvaluationResult,
) {
    if !settings.should_return_generated_events() {
        result.event_groups.clear();
    }
}

fn merge_event_processing_runtime<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    other: &mut I,
) -> Result<()> {
    match (
        integrand.event_processing_runtime_mut(),
        other.event_processing_runtime_mut(),
    ) {
        (Some(lhs), Some(rhs)) => lhs.merge_samples(rhs),
        _ => Ok(()),
    }
}

fn update_event_processing_runtime<I: ProcessIntegrandImpl>(integrand: &mut I, iter: usize) {
    if let Some(runtime) = integrand.event_processing_runtime_mut() {
        runtime.update_results(iter);
    }
}

fn observable_accumulator_bundle<I: ProcessIntegrandImpl>(
    integrand: &I,
) -> Option<ObservableAccumulatorBundle> {
    integrand
        .event_processing_runtime()
        .filter(|runtime| runtime.has_observables())
        .map(EventProcessingRuntime::accumulator_bundle)
}

fn observable_snapshot_bundle<I: ProcessIntegrandImpl>(
    integrand: &I,
) -> Option<ObservableSnapshotBundle> {
    integrand
        .event_processing_runtime()
        .filter(|runtime| runtime.has_observables())
        .map(EventProcessingRuntime::snapshot_bundle)
}

fn build_observable_snapshots_for_result<I: ProcessIntegrandImpl>(
    integrand: &I,
    result: &EvaluationResult,
) -> Option<ObservableSnapshotBundle> {
    build_observable_snapshots_for_event_groups(integrand, &result.event_groups)
}

fn build_observable_snapshots_for_event_groups<I: ProcessIntegrandImpl, T: FloatLike>(
    integrand: &I,
    event_groups: &crate::observables::GenericEventGroupList<T>,
) -> Option<ObservableSnapshotBundle> {
    let runtime = integrand.event_processing_runtime()?;
    if !runtime.has_observables() {
        return None;
    }

    let mut runtime = runtime.cleared_observable_clone();
    runtime.process_event_groups(event_groups);
    Some(runtime.snapshot_bundle())
}

fn restore_observable_snapshot_bundle<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    bundle: &ObservableSnapshotBundle,
) -> Result<()> {
    let runtime = integrand.event_processing_runtime_mut().ok_or_else(|| {
        eyre!("Cannot restore observables before the integrand has been warmed up")
    })?;
    if !runtime.has_observables() {
        return Err(eyre!(
            "Cannot restore observable snapshots for an integrand without configured observables"
        ));
    }
    runtime.restore_snapshot_bundle(bundle)
}

fn full_event_multiplicative_factor_precise<T: FloatLike>(
    parameterization_jacobian: Option<F<T>>,
    integrator_weight: F<T>,
) -> Complex<F<T>> {
    let jacobian = parameterization_jacobian.unwrap_or_else(|| integrator_weight.one());
    Complex::new_re(jacobian * integrator_weight)
}

pub(super) fn apply_full_event_multiplicative_factor_precise<T: FloatLike>(
    event_groups: &mut crate::observables::GenericEventGroupList<T>,
    full_factor: &Complex<F<T>>,
) {
    for event_group in event_groups.iter_mut() {
        for event in event_group.iter_mut() {
            event.apply_multiplicative_factor(full_factor);

            if !event.additional_weights.weights.is_empty() {
                event
                    .additional_weights
                    .weights
                    .entry(AdditionalWeightKey::FullMultiplicativeFactor)
                    .and_modify(|value| *value *= full_factor.clone())
                    .or_insert_with(|| full_factor.clone());
            }
        }
    }
}

pub(crate) fn write_observable_snapshot_bundle(
    bundle: &ObservableSnapshotBundle,
    path: &Path,
    format: ObservableFileFormat,
) -> Result<()> {
    match format {
        ObservableFileFormat::None => Ok(()),
        ObservableFileFormat::Hwu => bundle.write_hwu_file(path),
        ObservableFileFormat::Json => bundle.to_json_file(path),
    }
}

#[derive(Debug, Clone, Copy)]
pub enum IntegrandType {
    Amplitude,
    CrossSection,
}

fn create_stability_iterator(
    settings: &StabilitySettings,
    use_arb_prec: bool,
) -> Vec<StabilityLevelSetting> {
    if use_arb_prec {
        // overwrite the stability settings if use_f128 is enabled, but attempt to use user defined settings for f128
        if let Some(f128_settings_position) = settings
            .levels
            .iter()
            .position(|stability_level_setting| stability_level_setting.precision == Precision::Arb)
        {
            vec![settings.levels[f128_settings_position]]
        } else {
            vec![StabilityLevelSetting::default_arb()]
        }
    } else {
        settings.levels.clone()
    }
}

#[inline]
fn complex_from_f64<T: FloatLike>(value: &Complex<F<f64>>) -> Complex<F<T>> {
    Complex::new(F::<T>::from_ff64(value.re), F::<T>::from_ff64(value.im))
}

#[inline]
fn complex_to_f64<T: FloatLike>(value: &Complex<F<T>>) -> Complex<F<f64>> {
    Complex::new(value.re.into_ff64(), value.im.into_ff64())
}

type StabilityCheckResult<T> = (
    Complex<F<T>>,
    Option<F<T>>,
    bool,
    Option<StabilityFailureReason>,
);

#[inline]
fn stability_check<T: FloatLike>(
    ecm_scale: Option<&F<T>>,
    results: &[Complex<F<T>>],
    stability_settings: &StabilityLevelSetting,
    max_eval: Complex<F<T>>,
    wgt: F<T>,
    is_final_level: bool,
    escalate_if_exact_zero: bool,
) -> StabilityCheckResult<T> {
    stability_check_components(
        ecm_scale,
        results,
        stability_settings,
        max_eval,
        wgt,
        is_final_level,
        escalate_if_exact_zero,
        true,
        true,
    )
}

#[inline]
fn stability_check_components<T: FloatLike>(
    ecm_scale: Option<&F<T>>,
    results: &[Complex<F<T>>],
    stability_settings: &StabilityLevelSetting,
    max_eval: Complex<F<T>>,
    wgt: F<T>,
    is_final_level: bool,
    escalate_if_exact_zero: bool,
    check_real: bool,
    check_imag: bool,
) -> StabilityCheckResult<T> {
    // Nonfinite probes cannot establish stability, even at the final precision.
    if results.iter().any(|result| {
        (check_real && (result.re.is_nan() || result.re.is_infinite()))
            || (check_imag && (result.im.is_nan() || result.im.is_infinite()))
    }) {
        return (
            results[0].clone(),
            None,
            false,
            Some(StabilityFailureReason::ErrorThreshold),
        );
    }

    if results.len() == 1 {
        return (results[0].clone(), None, true, None);
    }

    let average = results
        .iter()
        .skip(1)
        .fold(results[0].clone(), |acc, x| acc + x)
        / F::<T>::from_f64(results.len() as f64);

    // Bound the reported rotation average without cancellation, applying the
    // remaining weight before any normal-range test. Native values and their
    // measured relative errors remain intact; only suppressed final components
    // may waive a relative-error failure at binary64's normal boundary.
    let weighted_absolute_average = results
        .iter()
        .fold(Complex::new_re(average.re.zero()), |sum, result| {
            sum + Complex::new((&result.re * &wgt).abs(), (&result.im * &wgt).abs())
        })
        / F::<T>::from_f64(results.len() as f64);
    let minimum_normal = F::<T>::from_f64(f64::MIN_POSITIVE);
    let real_underflow = check_real
        && weighted_absolute_average.re.0.is_finite()
        && weighted_absolute_average.re < minimum_normal;
    let imag_underflow = check_imag
        && weighted_absolute_average.im.0.is_finite()
        && weighted_absolute_average.im < minimum_normal;

    let errors = results.iter().map(|res| {
        let error_re = if !check_real {
            F::<T>::from_f64(0.0)
        } else if IsZero::is_zero(&res.re) && IsZero::is_zero(&average.re) {
            F::<T>::from_f64(0.0)
        } else {
            ((&res.re - &average.re) / &average.re).abs()
        };
        let error_im = if !check_imag {
            F::<T>::from_f64(0.0)
        } else if IsZero::is_zero(&res.im) && IsZero::is_zero(&average.im) {
            F::<T>::from_f64(0.0)
        } else {
            ((&res.im - &average.im) / &average.im).abs()
        };
        Complex::new(error_re, error_im)
    });
    let mut estimated_relative_accuracy = average.re.zero();

    let mut unstable_reason = None;
    let mut unstable_sample = None;
    for (index, error) in errors.enumerate() {
        estimated_relative_accuracy =
            estimated_relative_accuracy.max(error.re.clone().max(error.im.clone()));
        if !is_final_level
            && escalate_if_exact_zero
            && (!check_real || error.re == F::<T>::from_f64(0.0))
            && (!check_imag || error.im == F::<T>::from_f64(0.0))
        {
            unstable_reason = Some(StabilityFailureReason::ZeroError);
            unstable_sample = Some(index);
            break;
        }

        // The optional allowance is dimensionless relative to E_cm, after
        // every sample factor. Preserve the measured relative discrepancy;
        // another component's scale must not supply the allowance.
        let absolute_error = Complex::new(
            ((&results[index].re - &average.re) * &wgt).abs(),
            ((&results[index].im - &average.im) * &wgt).abs(),
        );
        let real_absolute = stability_settings.ecm_relative_tolerance_for_re > 0.0
            && wgt.0.is_finite()
            && weighted_absolute_average.re.0.is_finite()
            && absolute_error.re.0.is_finite()
            && ecm_scale.is_some_and(|scale| {
                scale.0.is_finite()
                    && scale > &scale.zero()
                    && &absolute_error.re / scale
                        <= F::<T>::from_f64(stability_settings.ecm_relative_tolerance_for_re)
            });
        let imag_absolute = stability_settings.ecm_relative_tolerance_for_im > 0.0
            && wgt.0.is_finite()
            && weighted_absolute_average.im.0.is_finite()
            && absolute_error.im.0.is_finite()
            && ecm_scale.is_some_and(|scale| {
                scale.0.is_finite()
                    && scale > &scale.zero()
                    && &absolute_error.im / scale
                        <= F::<T>::from_f64(stability_settings.ecm_relative_tolerance_for_im)
            });
        if (check_real
            && error.re > F::<T>::from_f64(stability_settings.required_precision_for_re)
            && !real_underflow
            && !real_absolute)
            || (check_imag
                && error.im > F::<T>::from_f64(stability_settings.required_precision_for_im)
                && !imag_underflow
                && !imag_absolute)
        {
            unstable_reason = Some(StabilityFailureReason::ErrorThreshold);
            unstable_sample = Some(index);
            break;
        }
    }

    if let Some(unstable_index) = unstable_sample {
        let unstable_point = &results[unstable_index];

        let ((real_formatted, rotated_real_formatted), (imag_formatted, rotated_imag_formatted)) = (
            format_for_compare_digits(
                average.re.clone().into_ff64(),
                unstable_point.re.clone().into_ff64(),
            ),
            format_for_compare_digits(
                average.im.clone().into_ff64(),
                unstable_point.im.clone().into_ff64(),
            ),
        );

        debug!("{}", "\nUnstable point detected:".red());
        debug!("\taverage result: {} + {}i", real_formatted, imag_formatted,);
        debug!(
            "\trotated result: {} + {}i",
            rotated_real_formatted, rotated_imag_formatted,
        );
    }

    let stable = unstable_sample.is_none();

    let has_active_max = (check_real && !IsZero::is_zero(&max_eval.re))
        || (check_imag && !IsZero::is_zero(&max_eval.im));
    let below_wgt_threshold =
        if stability_settings.escalate_for_large_weight_threshold > 0. && has_active_max {
            (check_real
                && average.re.abs() * wgt.clone()
                    < F::<T>::from_f64(stability_settings.escalate_for_large_weight_threshold)
                        * max_eval.re)
                || (check_imag
                    && average.im.abs() * wgt
                        < F::<T>::from_f64(stability_settings.escalate_for_large_weight_threshold)
                            * max_eval.im)
        } else {
            true
        };

    let weight_reason = if stable && !below_wgt_threshold {
        Some(StabilityFailureReason::WeightThreshold)
    } else {
        None
    };

    (
        average,
        Some(estimated_relative_accuracy),
        stable && below_wgt_threshold,
        unstable_reason.or(weight_reason),
    )
}

#[inline]
#[cfg(test)]
fn stability_check_on_norm<T: FloatLike>(
    ecm_scale: Option<&F<T>>,
    results: &[Complex<F<T>>],
    stability_settings: &StabilityLevelSetting,
    max_eval: Complex<F<T>>,
    wgt: F<T>,
    is_final_level: bool,
    escalate_if_exact_zero: bool,
) -> StabilityCheckResult<T> {
    stability_check_on_norm_components(
        ecm_scale,
        results,
        stability_settings,
        max_eval,
        wgt,
        is_final_level,
        escalate_if_exact_zero,
        true,
        true,
    )
}

#[inline]
fn stability_check_on_norm_components<T: FloatLike>(
    ecm_scale: Option<&F<T>>,
    results: &[Complex<F<T>>],
    stability_settings: &StabilityLevelSetting,
    max_eval: Complex<F<T>>,
    wgt: F<T>,
    is_final_level: bool,
    escalate_if_exact_zero: bool,
    check_real: bool,
    check_imag: bool,
) -> StabilityCheckResult<T> {
    // Nonfinite probes cannot establish stability, even at the final precision.
    if results.iter().any(|result| {
        (check_real && (result.re.is_nan() || result.re.is_infinite()))
            || (check_imag && (result.im.is_nan() || result.im.is_infinite()))
    }) {
        return (
            results[0].clone(),
            None,
            false,
            Some(StabilityFailureReason::ErrorThreshold),
        );
    }

    if results.len() == 1 {
        return (results[0].clone(), None, true, None);
    }

    let component_magnitude = |x: &Complex<F<T>>| match (check_real, check_imag) {
        (true, false) => x.re.abs(),
        (false, true) => x.im.abs(),
        (true, true) => x.norm_squared().sqrt(),
        (false, false) => x.re.zero(),
    };
    let average = results
        .iter()
        .fold(F::<T>::from_f64(0.0), |acc, x| acc + component_magnitude(x))
        / F::<T>::from_f64(results.len() as f64);

    // The norm owner returns the primary probe, so bound it as well as the
    // average. Componentwise L1 magnitudes avoid squaring tiny weighted values.
    let weighted_magnitude = |result: &Complex<F<T>>| &component_magnitude(result) * &wgt;
    let weighted_absolute_average = results.iter().fold(average.zero(), |sum, result| {
        sum + weighted_magnitude(result)
    }) / F::<T>::from_f64(results.len() as f64);
    let primary_magnitude = weighted_magnitude(&results[0]);
    let minimum_normal = F::<T>::from_f64(f64::MIN_POSITIVE);
    let underflow = weighted_absolute_average.0.is_finite()
        && primary_magnitude.0.is_finite()
        && weighted_absolute_average < minimum_normal
        && primary_magnitude < minimum_normal;

    // The E_cm-relative allowance must bound every component against the
    // primary value returned by this owner, independently of its norm.
    let absolute_agreement = (stability_settings.ecm_relative_tolerance_for_re > 0.0
        || stability_settings.ecm_relative_tolerance_for_im > 0.0)
        && wgt.0.is_finite()
        && weighted_absolute_average.0.is_finite()
        && ecm_scale.is_some_and(|scale| scale.0.is_finite() && scale > &scale.zero())
        && results.iter().all(|result| {
            let error_re = ((&result.re - &results[0].re) * &wgt).abs();
            let error_im = ((&result.im - &results[0].im) * &wgt).abs();
            (!check_real
                || (error_re.0.is_finite()
                    && error_re / ecm_scale.unwrap()
                        <= F::<T>::from_f64(stability_settings.ecm_relative_tolerance_for_re)))
                && (!check_imag
                    || (error_im.0.is_finite()
                        && error_im / ecm_scale.unwrap()
                            <= F::<T>::from_f64(stability_settings.ecm_relative_tolerance_for_im)))
        });

    let errors = results.iter().map(|res| {
        let res = component_magnitude(res);
        if IsZero::is_zero(&res) && IsZero::is_zero(&average) {
            (F::<T>::from_f64(0.0), true) // true zero is fishy -> upgrade to next precision
        } else {
            (((res - average.clone()) / average.clone()).abs(), false)
        }
    });
    let mut estimated_relative_accuracy = average.zero();

    let mut unstable_reason = None;
    let mut unstable_sample = None;
    for (index, (error, result_is_exact_zero)) in errors.enumerate() {
        estimated_relative_accuracy = estimated_relative_accuracy.max(error.clone());
        if !is_final_level
            && error == F::<T>::from_f64(0.0)
            && result_is_exact_zero
            && escalate_if_exact_zero
        {
            unstable_reason = Some(StabilityFailureReason::ZeroError);
            unstable_sample = Some(index);
            break;
        }

        if error > F::<T>::from_f64(stability_settings.required_precision_for_re)
            && !underflow
            && !absolute_agreement
        {
            unstable_reason = Some(StabilityFailureReason::ErrorThreshold);
            unstable_sample = Some(index);
            break;
        }
    }

    if let Some(unstable_index) = unstable_sample {
        let unstable_point = &results[unstable_index];

        let (real_formatted, rotated_real_formatted) = format_for_compare_digits(
            average.clone().into_ff64(),
            unstable_point.re.clone().into_ff64(),
        );

        debug!("{}", "\nUnstable point detected:".red());
        debug!("\tnormed average result: {}", real_formatted,);
        debug!("\tnormed rotated result: {}", rotated_real_formatted,);
    }

    let stable = unstable_sample.is_none();

    let max_magnitude = component_magnitude(&max_eval);
    let below_wgt_threshold = if stability_settings.escalate_for_large_weight_threshold > 0.
        && max_magnitude != max_magnitude.zero()
    {
        average.abs() * wgt
            < F::<T>::from_f64(stability_settings.escalate_for_large_weight_threshold)
                * max_magnitude
    } else {
        true
    };

    let weight_reason = if stable && !below_wgt_threshold {
        Some(StabilityFailureReason::WeightThreshold)
    } else {
        None
    };

    (
        results[0].clone(),
        Some(estimated_relative_accuracy),
        stable && below_wgt_threshold,
        unstable_reason.or(weight_reason),
    )
}

#[derive(Debug, Clone)]
struct PreciseStabilityLevelResult<T: FloatLike> {
    pub result: Complex<F<T>>,
    pub graph_result: GraphEvaluationResult<T>,
    pub stability_level_used: Precision,
    pub estimated_relative_accuracy: Option<F<T>>,
    pub sample_count: usize,
    pub total_time: Duration,
    pub parameterization_jacobian: Option<F<T>>,
    pub is_stable: bool,
    pub rotated_results: Vec<RotatedEvaluation>,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct LmbMultiChannelingSetup {
    pub lmb_basis_ids: TiVec<SamplingChannelId, LmbIndex>,
    /// Canonical group-master graph used to resolve group-level channel overrides.
    pub graph: Graph,
    pub all_bases: TiVec<LmbIndex, LoopMomentumBasis>,
    /// Current masses of the master graph; member masses never define its partition.
    pub(crate) master_edge_masses: RuntimeCache<EdgeVec<Complex<F<f64>>>>,
    pub(crate) sampling_bridge:
        RuntimeCache<Result<SamplingChannelBridge, sampling_maps::SamplingEvaluationError>>,
    pub(crate) sampling_bridge_quad:
        RuntimeCache<Result<SamplingChannelBridge<f128>, sampling_maps::SamplingEvaluationError>>,
    pub(crate) sampling_bridge_fixed256: RuntimeCache<
        Result<SamplingChannelBridge<SamplingFloat>, sampling_maps::SamplingEvaluationError>,
    >,
    pub(crate) sampling_bridge_arb: RuntimeCache<
        Result<SamplingChannelBridge<ArbPrec>, sampling_maps::SamplingEvaluationError>,
    >,
    /// One fixed source precision and strictest accuracy budget for the whole
    /// integrand epoch. RuntimeCache contributes no bytes to saved states.
    pub(crate) sampling_source: RuntimeCache<(SamplingPrecision, f64)>,
    pub(crate) sampling_catalogue: RuntimeCache<SamplingChannelCatalogue>,
    pub(crate) sampling_programs: RuntimeCache<Vec<SamplingChannelPrograms>>,
}

impl LmbMultiChannelingSetup {
    /// Refresh master-owned mass evaluators independently of member parameter slots.
    /// Keep complex values until the selected score decides whether it supports them.
    pub(crate) fn warm_up_masses(&mut self, settings: &RuntimeSettings, model: &Model) {
        self.master_edge_masses.invalidate();
        let parameters = &mut self.graph.param_builder;
        parameters.m_uv_value(Complex::new_re(F(settings.general.m_uv)));
        parameters.renormalization_localization_scale_value(Complex::new_re(F(settings
            .general
            .renormalization_localization_scale)));
        parameters.mu_r_sq_value(Complex::new_re(F(settings.general.mu_r_sq())));
        parameters.numerator_sampling_scale_value(Complex::new_re(F(settings
            .general
            .numerator_sampling_scale)));
        parameters.update_model_values(model);
        self.master_edge_masses
            .set(self.graph.new_edgevec(|edge, _, _| {
                edge.mass_value(model, &self.graph.param_builder)
                    .unwrap_or_else(|| Complex::new_re(F(0.0)))
            }));
    }

    /// Borrow the bridge compiled in the current successful warmup epoch, or
    /// replay that precision's cached numerical failure into the stability loop.
    /// Explicit constructors remain fresh and never populate this runtime cache.
    pub fn sampling_bridge<T: FloatLike>(&self) -> Result<&SamplingChannelBridge<T>> {
        T::sampling_bridge_cache(self).as_ref().ok_or_else(|| {
            eyre!(
                "sampling bridge for graph '{}' is not initialized; call warm_up after loading or changing runtime settings, model parameters, or graph routing",
                self.graph.name
            )
        })?.as_ref().map_err(|error| {
            eyre::Report::new(error.clone()).wrap_err(format!(
                "sampling binding for graph '{}' at {} precision",
                self.graph.name, T::sampling_precision()
            ))
        })
    }

    /// Invalidate all numerical bindings and their single canonical program epoch.
    pub(crate) fn invalidate_sampling(&mut self) {
        self.sampling_bridge.invalidate();
        self.sampling_bridge_quad.invalidate();
        self.sampling_bridge_fixed256.invalidate();
        self.sampling_bridge_arb.invalidate();
        self.sampling_source.invalidate();
        self.sampling_catalogue.invalidate();
        self.sampling_programs.invalidate();
    }

    /// Expand a graph-resolved selection into the canonical catalogue used by
    /// inspection and map construction. Resolve with the actual graph name
    /// before invoking this method (the setup may be shared by a graph group).
    pub fn sampling_channel_catalogue(
        &self,
        resolved: &ResolvedSamplingChannelSelection,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<SamplingChannelCatalogue> {
        self.sampling_channel_catalogue_with_surfaces(resolved, parameterization_settings, &[])
    }

    /// Resolve the canonical catalogue and return its read-only inspection
    /// view.  This is useful to CLI/Python diagnostics and intentionally does
    /// not compile maps, prepare kinematics, or enumerate a second channel
    /// domain.
    pub fn inspect_sampling_channels(
        &self,
        resolved: &ResolvedSamplingChannelSelection,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<SamplingChannelInspection> {
        Ok(self
            .sampling_channel_catalogue(resolved, parameterization_settings)?
            .inspection())
    }

    /// Report the conservative full-domain and elementary soft coverage of
    /// the canonical catalogue.  This is an inspection operation only: it
    /// does not discover physical E-surfaces or construct a second channel
    /// enumeration.
    pub fn sampling_channel_coverage_report(
        &self,
        resolved: &ResolvedSamplingChannelSelection,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<SamplingCoverageReport> {
        let catalogue = self.sampling_channel_catalogue(resolved, parameterization_settings)?;
        let all_lmbs = self
            .all_bases
            .iter_enumerated()
            .map(|(basis_id, basis)| {
                (
                    usize::from(basis_id),
                    basis.loop_edges.iter().map(|edge| edge.0).collect(),
                )
            })
            .collect::<Vec<_>>();
        let massless_edges = self
            .graph
            .underlying
            .iter_edges()
            .filter_map(|(_, edge_id, edge)| edge.data.particle.is_massless().then_some(edge_id.0))
            .collect::<Vec<_>>();
        Ok(catalogue.coverage_report(&all_lmbs, &massless_edges))
    }

    /// Resolve the same selection while supplying E-surface candidates already
    /// enumerated in the master graph frame. Surface existence and geometry are
    /// still prepared per cut/orientation before compilation.
    pub fn sampling_channel_catalogue_with_surfaces(
        &self,
        resolved: &ResolvedSamplingChannelSelection,
        parameterization_settings: &ParameterizationSettings,
        surface_edges: &[Vec<usize>],
    ) -> Result<SamplingChannelCatalogue> {
        let all_lmbs = self
            .all_bases
            .iter_enumerated()
            .map(|(basis_id, basis)| {
                (
                    usize::from(basis_id),
                    basis.loop_edges.iter().map(|edge| edge.0).collect(),
                )
            })
            .collect::<Vec<_>>();
        let optimized_lmbs = parameterization_settings
            .lmb_basis_ids
            .get(&self.graph.name)
            .cloned()
            .unwrap_or_else(|| {
                self.lmb_basis_ids
                    .iter()
                    .map(|basis| usize::from(*basis))
                    .collect()
            });
        let parent_lmb = self
            .graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .map(|edge| edge.0)
            .collect::<Vec<_>>();
        let massless_edges = self
            .graph
            .underlying
            .iter_edges()
            .filter_map(|(_, edge_id, edge)| edge.data.particle.is_massless().then_some(edge_id.0))
            .collect::<Vec<_>>();
        Ok(build_sampling_channel_catalogue_with_surfaces_and_coverage(
            resolved,
            &all_lmbs,
            &optimized_lmbs,
            surface_edges,
            &parent_lmb,
            &massless_edges,
        ))
    }

    /// Compile the selected catalogue against prepared master-graph
    /// kinematics. The context binds each surface's native equation and frame;
    /// any cut/side centers, radii and `t*` are prepared from its declared prior
    /// blocks inside the compiled map. This rejects missing or stale physical
    /// prerequisites before a mapped sample reaches the graph evaluator.
    pub fn compile_sampling_channels<T: FloatLike>(
        &self,
        catalogue: &SamplingChannelCatalogue,
        programs: &[SamplingChannelPrograms],
        context: &SamplingChannelCompileContext<T>,
    ) -> Result<Vec<CompiledSamplingChannel<T>>> {
        if context.master_graph != self.graph.name {
            return Err(eyre!(
                "sampling channel master graph `{}` does not match LMB setup graph `{}`",
                context.master_graph,
                self.graph.name
            ));
        }
        let mut context = context.clone();
        if let Some(masses) = self.master_edge_masses.as_ref() {
            context.edge_masses = masses
                .iter()
                .map(|(edge, mass)| (edge.0, mass.map_ref(|value| F::<T>::from_ff64(*value))))
                .collect();
        }
        catalogue.compile(&context, programs)
    }

    /// Compile channels after preparing the external data needed by every
    /// selected-LMB-to-parent affine routing.  This is the required entry
    /// point when a channel catalogue contains more than the parent LMB;
    /// callers must pass external data from the same solved cut/orientation
    /// context as the map geometry.
    pub fn compile_sampling_channels_with_external<T: FloatLike>(
        &self,
        catalogue: &SamplingChannelCatalogue,
        programs: &[SamplingChannelPrograms],
        context: &SamplingChannelCompileContext<T>,
        external_momenta: &[[T; 4]],
    ) -> Result<Vec<CompiledSamplingChannel<T>>> {
        let mut context = context.clone();
        for (basis_id, edges) in catalogue.lmb_basis_entries() {
            if edges != context.parent_lmb.as_slice() {
                context.lmb_frame_maps.insert(
                    basis_id,
                    self.lmb_frame_map(
                        &self.all_bases[LmbIndex::from(basis_id)],
                        external_momenta,
                    )?,
                );
            }
        }
        for channel in catalogue.named_entries() {
            let edges = match &channel.map {
                SamplingMapDefinition::Lmb(edges) => edges,
                _ => &channel.definition.parent_lmb,
            };
            if edges == context.parent_lmb.as_slice() {
                continue;
            }
            let basis = self.sampling_parent_lmb(edges)?;
            context
                .lmb_frame_maps_by_edges
                .insert(edges.clone(), self.lmb_frame_map(&basis, external_momenta)?);
        }
        self.compile_sampling_channels(catalogue, programs, &context)
    }

    /// Compile the selected channels and bind them to the raw-frame bridge.
    /// The process sampler supplies the prepared frame and external data.
    pub fn compile_sampling_channel_bridge<T: FloatLike>(
        &self,
        catalogue: &SamplingChannelCatalogue,
        programs: &[SamplingChannelPrograms],
        context: &SamplingChannelCompileContext<T>,
    ) -> Result<SamplingChannelBridge<T>> {
        let channels = self.compile_sampling_channels(catalogue, programs, context)?;
        let mode = match context.parameterization_settings.sampling_channels.weight {
            SamplingChannelWeight::MapDensity | SamplingChannelWeight::InverseJacobian => {
                SamplingPartitionMode::MapDensity
            }
            SamplingChannelWeight::Ose => SamplingPartitionMode::Ose,
            SamplingChannelWeight::SingularityProxy => SamplingPartitionMode::SingularityProxy,
        };
        SamplingChannelBridge::new_with_partition_mode(channels, mode).map_err(Into::into)
    }

    /// External-data variant of [`Self::compile_sampling_channel_bridge`].
    pub fn compile_sampling_channel_bridge_with_external<T: FloatLike>(
        &self,
        catalogue: &SamplingChannelCatalogue,
        programs: &[SamplingChannelPrograms],
        context: &SamplingChannelCompileContext<T>,
        external_momenta: &[[T; 4]],
    ) -> Result<SamplingChannelBridge<T>> {
        let channels = self.compile_sampling_channels_with_external(
            catalogue,
            programs,
            context,
            external_momenta,
        )?;
        let mode = match context.parameterization_settings.sampling_channels.weight {
            SamplingChannelWeight::MapDensity | SamplingChannelWeight::InverseJacobian => {
                SamplingPartitionMode::MapDensity
            }
            SamplingChannelWeight::Ose => SamplingPartitionMode::Ose,
            SamplingChannelWeight::SingularityProxy => SamplingPartitionMode::SingularityProxy,
        };
        SamplingChannelBridge::new_with_partition_mode(channels, mode).map_err(Into::into)
    }

    fn validate_lmb_basis_id(&self, basis_id: usize, graph_name: &str) -> Result<LmbIndex> {
        if basis_id >= self.all_bases.len() {
            return Err(eyre!(
                "Requested LMB basis id {} is out of range for graph '{}'; the graph has {} generated LMB bases.",
                basis_id,
                graph_name,
                self.all_bases.len()
            ));
        }
        Ok(LmbIndex::from(basis_id))
    }

    /// Build the one canonical catalogue used by all channel-count and LMB
    /// lookup helpers during the runtime migration.  The old generated basis
    /// list is only input data for this catalogue; it is never enumerated as a
    /// second channel universe.
    fn canonical_sampling_catalogue(
        &self,
        graph_name: &str,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<SamplingChannelCatalogue> {
        let resolved = resolve_sampling_channel_selection(
            graph_name,
            &parameterization_settings.sampling_channels,
        )?;
        self.sampling_channel_catalogue(&resolved, parameterization_settings)
    }

    /// Return the stable IDs from the one canonical catalogue.
    pub fn sampling_channel_ids(
        &self,
        graph_name: &str,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<Vec<SamplingChannelId>> {
        let catalogue = self.canonical_sampling_catalogue(graph_name, parameterization_settings)?;
        Ok((0..catalogue.entries.len())
            .map(SamplingChannelId::from)
            .collect())
    }

    pub fn sampling_channel_label(
        &self,
        channel_id: SamplingChannelId,
        graph_name: &str,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<String> {
        let catalogue = self.canonical_sampling_catalogue(graph_name, parameterization_settings)?;
        let entry = catalogue.entries.get(channel_id.index()).ok_or_else(|| {
            eyre!(
                "Requested sampling channel {} is out of range for graph '{}'",
                channel_id.index(),
                graph_name
            )
        })?;
        Ok(match entry {
            SamplingCatalogueEntry::Lmb {
                basis_id, edges, ..
            } => {
                format!("lmb[{basis_id}] {:?}", edges)
            }
            SamplingCatalogueEntry::Surface { edges, .. } => format!("surface:{edges:?}"),
            SamplingCatalogueEntry::Named(channel) => channel.name.clone(),
        })
    }

    pub fn sampling_channel_is_lmb(
        &self,
        channel_id: SamplingChannelId,
        graph_name: &str,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<bool> {
        let catalogue = self.canonical_sampling_catalogue(graph_name, parameterization_settings)?;
        let entry = catalogue.entries.get(channel_id.index()).ok_or_else(|| {
            eyre!(
                "Requested sampling channel {} is out of range for graph '{}'",
                channel_id.index(),
                graph_name
            )
        })?;
        Ok(matches!(entry, SamplingCatalogueEntry::Lmb { .. }))
    }

    /// Return the generated LMB behind a canonical channel when that channel
    /// is an LMB entry. Graph-aware channels intentionally have no LMB index;
    /// callers that only expose LMB metadata can skip those entries without
    /// creating a second channel enumeration.
    pub fn sampling_channel_lmb_id(
        &self,
        channel_id: SamplingChannelId,
        graph_name: &str,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<Option<LmbIndex>> {
        let catalogue = self.canonical_sampling_catalogue(graph_name, parameterization_settings)?;
        let entry = catalogue.entries.get(channel_id.index()).ok_or_else(|| {
            eyre!(
                "Requested sampling channel {} is out of range for graph '{}'",
                channel_id.index(),
                graph_name
            )
        })?;
        match entry {
            SamplingCatalogueEntry::Lmb { basis_id, .. } => {
                self.validate_lmb_basis_id(*basis_id, graph_name).map(Some)
            }
            SamplingCatalogueEntry::Surface { .. } | SamplingCatalogueEntry::Named(_) => Ok(None),
        }
    }

    /// Resolve a canonical sampling channel to its generated LMB basis.
    ///
    /// This strict accessor is intentionally named in terms of the canonical
    /// sampling catalogue.  It is only valid for channels whose map is an LMB;
    /// graph-aware channels must be evaluated through the sampling bridge.
    pub fn sampling_channel_lmb_basis_id(
        &self,
        channel_id: SamplingChannelId,
        graph_name: &str,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<LmbIndex> {
        self.sampling_channel_lmb_id(channel_id, graph_name, parameterization_settings)?
            .ok_or_else(|| {
                eyre!(
                    "Sampling channel {} for graph '{}' is graph-aware and has no generated LMB basis; use the canonical sampling bridge for this channel.",
                    channel_id.index(), graph_name
                )
            })
    }

    pub fn sampling_channel_edge_ids(
        &self,
        channel_id: SamplingChannelId,
        graph_name: &str,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<SmallVec<[usize; 4]>> {
        let lmb_index = self
            .sampling_channel_lmb_id(channel_id, graph_name, parameterization_settings)?
            .ok_or_else(|| {
                eyre!(
                    "Sampling channel {} for graph '{}' is graph-aware and has no generated LMB basis; use the canonical sampling bridge for this channel.",
                    channel_id.index(), graph_name
                )
            })?;
        Ok(self.all_bases[lmb_index]
            .loop_edges
            .iter()
            .map(|edge_id| edge_id.0)
            .collect())
    }

    pub fn selected_lmb_basis_id(
        &self,
        graph_name: &str,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<LmbIndex> {
        let channels = self.canonical_lmb_basis_entries(graph_name, parameterization_settings)?;
        channels.first().map(|(_, basis_id)| *basis_id).ok_or_else(|| {
            eyre!(
                "Could not select a default LMB basis for graph '{}'; the optimized LMB subset is empty.",
                graph_name
            )
        })
    }

    /// Resolve generated bases for default sampling's LMB choice. IDs retain
    /// their catalogue positions; this view never creates another channel axis.
    /// Multichannel sampling uses the compiled bridge for every channel.
    fn canonical_lmb_basis_entries(
        &self,
        graph_name: &str,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<Vec<(SamplingChannelId, LmbIndex)>> {
        let catalogue = self.canonical_sampling_catalogue(graph_name, parameterization_settings)?;
        let mut channels = Vec::with_capacity(catalogue.entries.len());
        for (index, entry) in catalogue.entries.iter().enumerate() {
            match entry {
                SamplingCatalogueEntry::Lmb { basis_id, .. } => channels.push((
                    SamplingChannelId::from(index),
                    self.validate_lmb_basis_id(*basis_id, graph_name)?,
                )),
                SamplingCatalogueEntry::Surface { .. } | SamplingCatalogueEntry::Named(_) => {
                    return Err(eyre!(
                        "LMB-only weighting was requested for graph '{graph_name}', but canonical sampling channel {} is graph-aware; use the advanced discrete sampling route for this selection",
                        index
                    ));
                }
            }
        }
        Ok(channels)
    }

    /// Resolve a complete user-ordered parent from a generated edge set without
    /// rebuilding its tree or external routing. Reorder a clone only: canonical
    /// catalogue IDs and the generated basis order are immutable.
    pub(crate) fn sampling_parent_lmb(&self, edges: &[usize]) -> Result<LoopMomentumBasis> {
        let selected = edges.iter().copied().collect::<BTreeSet<_>>();
        if edges.len() != self.graph.loop_momentum_basis.loop_edges.len()
            || selected.len() != edges.len()
        {
            return Err(eyre!(
                "sampling parent requires {} distinct loop edges, received {edges:?}",
                self.graph.loop_momentum_basis.loop_edges.len()
            ));
        }
        let candidates = self
            .all_bases
            .iter()
            .filter(|basis| {
                basis.loop_edges.len() == edges.len()
                    && basis
                        .loop_edges
                        .iter()
                        .all(|edge| selected.contains(&edge.0))
            })
            .collect_vec();
        let [basis] = candidates.as_slice() else {
            return Err(eyre!(
                "sampling parent {edges:?} matches {} generated bases for graph '{}'; expected one complete basis, available {:?}",
                candidates.len(),
                self.graph.name,
                self.all_bases
                    .iter()
                    .map(|basis| &basis.loop_edges)
                    .collect_vec()
            ));
        };
        let mut basis = (*basis).clone();
        for (index, edge) in edges.iter().enumerate() {
            let current = basis
                .loop_edges
                .iter_enumerated()
                .find(|(_, candidate)| candidate.0 == *edge)
                .unwrap()
                .0;
            basis.swap_loops(LoopIndex(index), current);
        }
        Ok(basis)
    }

    /// Build the exact affine routing from one generated LMB into this setup's
    /// parent loop frame.  The integer edge signatures provide the linear
    /// block matrix; the supplied external momenta provide its translation.
    /// Keeping this operation on the graph-aware setup ensures that a compiled
    /// channel never mistakes selected-LMB coordinates for parent-frame
    /// coordinates.
    pub fn lmb_frame_map<T: FloatLike>(
        &self,
        channel_lmb: &LoopMomentumBasis,
        external_momenta: &[[T; 4]],
    ) -> Result<SamplingMapAffine<T>> {
        let parent_loop_edges = &self.graph.loop_momentum_basis.loop_edges;
        let channel_loop_count = channel_lmb.loop_edges.len();
        if parent_loop_edges.len() != channel_loop_count {
            return Err(eyre!(
                "cannot route LMB {:?} with {} loop blocks into parent frame with {} blocks",
                &channel_lmb.loop_edges,
                channel_loop_count,
                parent_loop_edges.len()
            ));
        }
        let dimension = 3 * channel_loop_count;
        let zero = F::<T>::default().zero();
        let mut matrix = vec![vec![zero.0.clone(); dimension]; dimension];
        let mut translation = vec![zero.0.clone(); dimension];
        for (parent_block, &edge_index) in parent_loop_edges.iter().enumerate() {
            let signature = &channel_lmb.edge_signatures[edge_index];
            let internal = signature.internal.to_momtrop_format();
            if internal.len() != channel_loop_count {
                return Err(eyre!(
                    "LMB {:?} edge {} has internal signature length {}, expected {}",
                    &channel_lmb.loop_edges,
                    edge_index.0,
                    internal.len(),
                    channel_loop_count
                ));
            }
            let external = signature.external.to_momtrop_format();
            if external.len() != external_momenta.len() {
                return Err(eyre!(
                    "LMB {:?} edge {} has external signature length {}, but {} external momenta were supplied",
                    &channel_lmb.loop_edges,
                    edge_index.0,
                    external.len(),
                    external_momenta.len()
                ));
            }
            for component in 0..3 {
                let row = 3 * parent_block + component;
                for (channel_block, coefficient) in internal.iter().enumerate() {
                    matrix[row][3 * channel_block + component] =
                        zero.from_i64(*coefficient as i64).0;
                }
                translation[row] = external
                    .iter()
                    .zip(external_momenta)
                    .map(|(coefficient, momentum)| {
                        zero.from_i64(*coefficient as i64) * F(momentum[component + 1].clone())
                    })
                    .fold(zero.clone(), |sum, value| sum + value)
                    .0;
            }
        }
        SamplingMapAffine::new(matrix, translation)
    }

    fn reinterpret_loop_momenta_for_lmb_impl<T: FloatLike>(
        &self,
        lmb_index: LmbIndex,
        momentum_sample: &BareMomentumSample<T>,
        loop_mom_cache_id: usize,
    ) -> BareMomentumSample<T> {
        let channel_lmb = &self.all_bases[lmb_index];
        let new_loop_moms: LoopMomenta<F<T>> = self
            .graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .map(|&edge_index| {
                let signature_of_edge_channel_lmb = &channel_lmb.edge_signatures[edge_index];

                signature_of_edge_channel_lmb
                    .internal
                    .apply_typed(&momentum_sample.loop_moms)
                    + signature_of_edge_channel_lmb
                        .external
                        .apply(&momentum_sample.external_moms.raw)
                        .spatial
            })
            .collect();

        BareMomentumSample {
            loop_moms: new_loop_moms,
            dual_loop_moms: momentum_sample.dual_loop_moms.clone().map(|_dlm| todo!()),
            loop_mom_cache_id,
            loop_mom_base_cache_id: momentum_sample.loop_mom_base_cache_id,
            external_mom_cache_id: momentum_sample.external_mom_cache_id,
            external_mom_base_cache_id: momentum_sample.external_mom_base_cache_id,
            external_moms: momentum_sample.external_moms.clone(),
            jacobian: momentum_sample.jacobian.clone(),
            orientation: momentum_sample.orientation,
            parameterization_branch: momentum_sample.parameterization_branch,
        }
    }

    pub(crate) fn reinterpret_loop_momenta_for_lmb<T: FloatLike>(
        &self,
        lmb_index: LmbIndex,
        momentum_sample: &MomentumSample<T>,
        loop_mom_cache_id: usize,
    ) -> MomentumSample<T> {
        MomentumSample {
            sample: self.reinterpret_loop_momenta_for_lmb_impl(
                lmb_index,
                &momentum_sample.sample,
                loop_mom_cache_id,
            ),
        }
    }
}

pub trait ProcessIntegrandImpl {
    type G: GraphTerm;

    fn warm_up(&mut self, model: &Model) -> Result<()>;

    /// Natural E_cm scale of the returned physical quantity. Custom graphs
    /// declare their integrated energy dimension; raw inputs omit spatial measure.
    fn stability_reference_scale<T: FloatLike>(
        &self,
        graph_id: usize,
        sample: &MomentumSample<T>,
        missing_measure_dimension: i32,
    ) -> Result<F<T>> {
        let settings = self.get_settings();
        let dimension = settings.stability.integrated_energy_dimension.ok_or_else(|| {
            eyre!("E_cm-relative physical stability requires `stability.integrated_energy_dimension` (after flux normalization, before output-unit conversion)")
        })?;
        let exponent = dimension.checked_sub(missing_measure_dimension).ok_or_else(|| {
            eyre!("stability energy dimension overflows after subtracting the missing spatial measure")
        })?;
        let mut scale = F::<T>::from_f64(settings.kinematics.e_cm).powi(exponent);
        if matches!(
            self.get_dependent_momenta_constructor(),
            DependentMomentaConstructor::CrossSection
        ) && !settings.general.disable_flux_factor
            && sample.external_moms().len() == 2
        {
            let graph = self.get_graph(graph_id).get_graph();
            let unit = settings
                .general
                .integral_unit
                .resolve_for_cross_section(graph.initial_state_cut.iter_edges(graph).count());
            scale *= cross_section::barn_conversion_factor(unit, sample.one());
        }
        Ok(scale)
    }

    /// Choose the complete numerical proposal before any draw or physical retry.
    /// This conservative admission floor is not a pointwise map-error bound;
    /// density/root checks still fail explicitly when the fixed law is unresolved.
    fn sampling_source_policy(&self) -> Result<(SamplingPrecision, f64)> {
        let settings = self.get_settings();
        let budget = GammaLoopSample::<ArbPrec>::source_accuracy_budget(settings)?;
        let floor = F::<f128>::default().epsilon().sqrt().into_ff64().0;
        let parameterization = settings
            .sampling
            .get_parameterization_settings()
            .ok_or_else(|| eyre!("sampling source policy requires a channel parameterization"))?;
        let quad_eligible = budget >= floor
            && self.graph_count() > 0
            && parameterization.sampling_channels.weight != SamplingChannelWeight::SingularityProxy
            && (0..self.graph_count()).all(|id| {
                self.get_graph(id)
                    .sampling_setup()
                    .sampling_catalogue
                    .as_ref()
                    .is_some_and(SamplingChannelCatalogue::supports_fixed_quad_source)
            });
        let externals = settings
            .kinematics
            .externals
            .get_dependent_externals::<ArbPrec>(self.get_dependent_momenta_constructor())?;
        let mut inputs = externals
            .iter()
            .flat_map(|momentum| {
                [
                    momentum.temporal.value.clone(),
                    momentum.spatial.px.clone(),
                    momentum.spatial.py.clone(),
                    momentum.spatial.pz.clone(),
                ]
            })
            .collect_vec();
        let scale = F::<ArbPrec>::from_f64(settings.kinematics.e_cm);
        let b = F::<ArbPrec>::from_f64(parameterization.b);
        inputs.extend([
            scale.square(),
            &scale * &b,
            scale,
            b,
            F::<ArbPrec>::from_f64(parameterization.power),
            F::<ArbPrec>::from_f64(settings.lu_h_function.sigma),
        ]);
        for id in 0..self.graph_count() {
            for mass in self
                .get_graph(id)
                .get_real_mass_vector()
                .iter()
                .filter_map(|(_, mass)| mass.as_ref())
            {
                let mass = F::<ArbPrec>::from_ff64(*mass);
                inputs.extend([mass.square(), mass]);
            }
        }
        if inputs
            .iter()
            .any(|value| value.is_nan() || value.is_infinite())
        {
            return Err(eyre!(
                "sampling source has nonfinite fixed kinematic, mass or parameterization input"
            ));
        }
        let quad_representable = inputs
            .iter()
            .all(|value| F::<f128>::from_arb(&value.0).is_ok());
        let fixed256_floor = F::<SamplingFloat>::default().epsilon().sqrt().into_ff64().0;
        Ok((
            if quad_eligible && quad_representable {
                SamplingPrecision::Quad
            } else if budget >= fixed256_floor
                && inputs
                    .iter()
                    .all(|value| F::<SamplingFloat>::from_arb(&value.0).is_ok())
            {
                SamplingPrecision::Fixed256
            } else {
                SamplingPrecision::Arb
            },
            budget,
        ))
    }

    /// Compile fixed graph geometry after process warmup has prepared masses
    /// and improved externals. Publish only when the fixed canonical precision
    /// has every bridge valid. Native bindings remain available for component
    /// diagnostics; they do not define physical retry proposals. Radial roots
    /// and conditional cut contexts remain point-dependent runtime data.
    /// A failed canonical binding invalidates the epoch.
    fn warm_up_sampling(&mut self) -> Result<()> {
        for graph in self.get_terms_mut() {
            graph.sampling_setup_mut().invalidate_sampling();
        }
        if !self.get_settings().sampling.uses_sampling_channels() {
            return Ok(());
        }
        let parameterization = self
            .get_settings()
            .sampling
            .get_parameterization_settings()
            .expect("sampling channels require a parameterization");
        let prepared = (0..self.graph_count())
            .map(|id| {
                let graph = self.get_graph(id);
                let setup = graph.sampling_setup();
                let resolved = resolve_sampling_channel_selection(
                    &setup.graph.name,
                    &parameterization.sampling_channels,
                )?;
                // Inspect already resolved blocks without requiring a successful
                // binding at another precision; the original source will retain
                // one fixed canonical map for every physical rescue lane.
                let catalogue = setup.sampling_channel_catalogue(&resolved, &parameterization)?;
                let programs = catalogue.compile_programs(
                    3 * graph.get_graph().get_loop_number(),
                    &self.get_settings().lu_h_function,
                )?;
                Ok((catalogue, programs))
            })
            .collect::<Result<Vec<_>>>()?;
        for (graph, (catalogue, programs)) in self.get_terms_mut().zip(prepared) {
            let setup = graph.sampling_setup_mut();
            setup.sampling_catalogue.set(catalogue);
            setup.sampling_programs.set(programs);
        }
        let result = (|| {
            let (precision, budget) = self.sampling_source_policy()?;
            match precision {
                SamplingPrecision::Quad => self.prepare_sampling_precision::<f128>(),
                SamplingPrecision::Fixed256 => self.prepare_sampling_precision::<SamplingFloat>(),
                SamplingPrecision::Arb => self.prepare_sampling_precision::<ArbPrec>(),
                SamplingPrecision::Double => unreachable!("Double is never a canonical source"),
            }?;
            for graph in self.get_terms_mut() {
                graph
                    .sampling_setup_mut()
                    .sampling_source
                    .set((precision, budget));
            }
            crate::debug_tags!(#sampling;
                stage = "sampling_source_warmup", precision = %precision, accuracy_budget = budget,
                "selected fixed source for the complete integrand epoch"
            );
            Ok(())
        })();
        if result.is_ok() {
            // Preserve the configured native component APIs and their cached
            // numerical diagnostics, without requiring a second proposal law
            // to be usable before evaluating the canonical one.
            let precisions = self
                .get_settings()
                .stability
                .levels
                .iter()
                .map(|level| level.precision)
                .unique()
                .collect_vec();
            for precision in precisions {
                let native = match precision {
                    Precision::Double => self.prepare_sampling_precision::<f64>(),
                    Precision::Quad => self.prepare_sampling_precision::<f128>(),
                    Precision::Arb => self.prepare_sampling_precision::<ArbPrec>(),
                };
                if let Err(error) = native {
                    crate::debug_tags!(#sampling;
                        stage = "sampling_native_component_binding", precision = %precision,
                        error = %error, "native component binding unavailable; production uses the canonical draw"
                    );
                }
            }
        }
        if result.is_err() {
            // An epoch with no usable canonical binding publishes neither a
            // catalogue nor partial bindings, irrespective of native availability.
            for graph in self.get_terms_mut() {
                graph.sampling_setup_mut().invalidate_sampling();
            }
        }
        result
    }

    /// Bind a configured or explicitly requested precision once per warmup epoch,
    /// retaining typed numerical failures as well as successes. Programs and
    /// canonical IDs are reused; geometry comes from the same native
    /// improved external data consumed by physical evaluation. Worker clones own
    /// their evaluator buffers, and a precision request never reparses metadata.
    fn prepare_sampling_precision<T: FloatLike>(&mut self) -> Result<()> {
        if !self.get_settings().sampling.uses_sampling_channels() {
            return Ok(());
        }
        if (0..self.graph_count()).all(|id| {
            T::sampling_bridge_cache(self.get_graph(id).sampling_setup())
                .as_ref()
                .is_some()
        }) {
            return (0..self.graph_count()).try_for_each(|id| {
                self.get_graph(id)
                    .sampling_setup()
                    .sampling_bridge::<T>()
                    .map(|_| ())
            });
        }
        let parameterization = self
            .get_settings()
            .sampling
            .get_parameterization_settings()
            .expect("sampling channels require a parameterization");
        let externals = self
            .get_settings()
            .kinematics
            .externals
            .get_dependent_externals::<T>(self.get_dependent_momenta_constructor())?;
        let external_momenta = externals
            .iter()
            .map(|momentum| {
                [
                    momentum.temporal.value.0.clone(),
                    momentum.spatial.px.0.clone(),
                    momentum.spatial.py.0.clone(),
                    momentum.spatial.pz.0.clone(),
                ]
            })
            .collect_vec();
        // The materialized draw, map consistency and physical host adoption
        // share one accuracy budget; native physical retries need no map binding.
        let (source_precision, source_budget) = self
            .get_graph(0)
            .sampling_setup()
            .sampling_source
            .as_ref()
            .copied()
            .map_or_else(|| self.sampling_source_policy(), Ok)?;
        let density_tolerance = if source_precision == T::sampling_precision() {
            source_budget
        } else {
            GammaLoopSample::<T>::relative_accuracy_budget(self.get_settings())
        };
        let bridges = (0..self.graph_count()).map(|id| {
            let graph = self.get_graph(id);
            let setup = graph.sampling_setup();
            if T::sampling_bridge_cache(setup).as_ref().is_some() { return Ok(None); }
            let catalogue = setup.sampling_catalogue.as_ref().ok_or_else(|| eyre!(
                "sampling catalogue for graph '{}' is not initialized; call warm_up after loading or changing runtime settings, model parameters, or graph routing", graph.name()
            ))?;
            let programs = setup.sampling_programs.as_ref().ok_or_else(|| eyre!(
                "sampling programs for graph '{}' are not initialized; call warm_up", graph.name()
            ))?;
            let bridge = graph.bind_sampling_bridge(catalogue, programs, &parameterization, self.get_settings(), &external_momenta, None)
                .and_then(|bridge| bridge.with_relative_density_tolerance(density_tolerance));
            match bridge {
                Ok(bridge) => Ok(Some(Ok(bridge))),
                Err(error) => match error.downcast_ref::<sampling_maps::SamplingEvaluationError>() {
                    Some(error) => Ok(Some(Err(error.clone()))),
                    None => Err(error),
                },
            }
        }).collect::<Result<Vec<_>>>()?;
        for (graph, binding) in self.get_terms_mut().zip(bridges) {
            if let Some(binding) = binding {
                crate::debug_tags!(#sampling;
                    stage = "sampling_bridge_warmup", graph = %graph.name(),
                    channels = ?binding.as_ref().ok().map(|bridge| bridge.channels().len()),
                    precision = std::any::type_name::<T>(), valid = binding.is_ok(),
                    "prepared graph sampling binding"
                );
                T::sampling_bridge_cache_mut(graph.sampling_setup_mut()).set(binding);
            }
        }
        (0..self.graph_count()).try_for_each(|id| {
            self.get_graph(id)
                .sampling_setup()
                .sampling_bridge::<T>()
                .map(|_| ())
        })
    }

    fn get_rotations(&self) -> impl Iterator<Item = &Rotation>;

    fn increment_loop_cache_id(&mut self, val: usize);
    fn loop_cache_id(&self) -> usize;

    fn increment_external_cache_id(&mut self, val: usize);
    fn external_cache_id(&self) -> usize;

    /// Signal that external momenta configuration has actually changed
    /// This increments the external cache ID to invalidate cached computations
    ///
    /// # Usage
    /// Call this method when:
    /// - You change the external momenta values in your Monte Carlo sampling
    /// - You switch to a different kinematic configuration
    /// - You want to force recomputation of cached polarizations/external quantities
    ///
    /// # Example
    /// ```rust,ignore
    /// // When you change external momenta in your sampling
    /// integrand.signal_external_momenta_changed();
    /// let new_sample = parameterize(&sample_point, &mut integrand, &mut first_draw_metadata)?;
    ///
    /// // Rotations will automatically get new cache IDs
    /// let rotated_samples = evaluate_all_rotations(&new_sample, &mut integrand, true)?;
    ///
    /// // Next iteration - revert to base configuration to reuse cache
    /// integrand.revert_to_base_external_cache_id();
    /// let next_sample = parameterize(&next_sample_point, &mut integrand, &mut next_draw_metadata)?;
    /// ```
    fn signal_external_momenta_changed(&mut self) {
        self.increment_external_cache_id(1);
    }

    /// Get the current external cache ID for reuse (doesn't increment)
    ///
    /// This returns the current cache ID without incrementing it, allowing
    /// new samples to reuse cached computations for the same external
    /// momenta configuration.
    fn get_current_external_cache_id(&self) -> usize {
        self.external_cache_id()
    }

    /// Revert to the base external cache ID for the current configuration
    ///
    /// This allows new samples to reuse the base cache ID when they represent
    /// the same underlying external momenta configuration (e.g., after rotations)
    ///
    /// # Usage
    /// Call this method when you want to create a new sample that represents
    /// the same base external momenta configuration as before, allowing reuse
    /// of cached polarizations and other external-dependent computations.
    ///
    /// # Example Cache Flow
    /// ```text
    /// 1. Initial sample:           cache_id = 0, base_cache_id = 0
    /// 2. Rotation 1:               cache_id = 1, base_cache_id = 0
    /// 3. Rotation 2:               cache_id = 2, base_cache_id = 0
    /// 4. revert_to_base():         cache_id = 0, base_cache_id = 0  // Reuse cache!
    /// 5. signal_changed():         cache_id = 3, base_cache_id = 3  // New base config
    /// 6. Rotation of new config:   cache_id = 4, base_cache_id = 3
    /// 7. revert_to_base():         cache_id = 3, base_cache_id = 3  // Reuse new base
    /// ```
    fn revert_to_base_external_cache_id(&mut self);

    /// Check if external momenta caching is beneficial
    ///
    /// Returns true if the same external cache ID has been used multiple times,
    /// indicating that caching is providing benefits.
    fn is_external_caching_beneficial(&self) -> bool {
        // Simple heuristic: if we've created samples without incrementing
        // cache ID, then caching is being used
        self.external_cache_id() < self.loop_cache_id()
    }

    /// Force a cache consistency check with detailed reporting
    fn debug_cache_state(&self, context: &str) {
        if std::env::var("GAMMALOOP_DEBUG_CACHE").is_ok() {
            let validation = self.validate_cache_consistency();
            let stats = self.get_cache_stats();

            tracing::info!("🔍 DEBUG CACHE STATE at {}", context);
            tracing::info!("   Validation: {}", validation);
            tracing::info!("   Statistics: {}", stats);

            if !validation.is_valid {
                tracing::error!("   ❌ CACHE INCONSISTENCY DETECTED!");
                panic!(
                    "Cache corruption at {}: {}",
                    context, validation.diagnostics
                );
            }

            if validation.has_rotations {
                tracing::info!(
                    "   🔄 {} rotation variants from base cache_id {}",
                    validation.current_external_cache_id - validation.base_external_cache_id,
                    validation.base_external_cache_id
                );
            }

            if stats.efficiency_ratio < 0.3 {
                tracing::warn!(
                    "   ⚠️ Very low cache efficiency: {:.1}%",
                    stats.efficiency_ratio * 100.0
                );
            } else if stats.efficiency_ratio > 0.8 {
                tracing::info!(
                    "   ✅ Excellent cache efficiency: {:.1}%",
                    stats.efficiency_ratio * 100.0
                );
            }
        }
    }

    /// Get the base external cache ID for the current configuration
    fn get_base_external_cache_id(&self) -> usize;

    /// Validate cache ID consistency and return diagnostics
    fn validate_cache_consistency(&self) -> CacheValidationResult {
        let current_id = self.external_cache_id();
        let base_id = self.get_base_external_cache_id();
        let loop_id = self.loop_cache_id();

        let is_valid = base_id <= current_id;
        let has_rotations = current_id > base_id;
        let cache_efficiency = if loop_id > 0 {
            1.0 - (current_id as f64 / loop_id as f64)
        } else {
            0.0
        };

        CacheValidationResult {
            is_valid,
            current_external_cache_id: current_id,
            base_external_cache_id: base_id,
            loop_cache_id: loop_id,
            has_rotations,
            cache_efficiency,
            diagnostics: if is_valid {
                "Cache IDs are consistent".to_string()
            } else {
                format!(
                    "ERROR: Base cache ID ({}) > Current cache ID ({})",
                    base_id, current_id
                )
            },
        }
    }

    /// Get cache usage statistics for monitoring
    fn get_cache_stats(&self) -> CacheStats {
        let validation = self.validate_cache_consistency();
        CacheStats {
            total_external_increments: validation.current_external_cache_id,
            total_loop_increments: validation.loop_cache_id,
            base_configurations: validation.base_external_cache_id + 1,
            rotational_variants: validation.current_external_cache_id
                - validation.base_external_cache_id,
            efficiency_ratio: validation.cache_efficiency,
        }
    }

    fn get_group_masters(&self) -> impl Iterator<Item = &Self::G>;

    fn get_terms_mut(&mut self) -> impl Iterator<Item = &mut Self::G>;
    fn graph_count(&self) -> usize;
    fn get_settings(&self) -> &RuntimeSettings;
    fn get_master_graph(&self, group_id: GroupId) -> &Self::G;
    fn get_graph(&self, graph_id: usize) -> &Self::G;
    fn get_graph_mut(&mut self, graph_id: usize) -> &mut Self::G;
    fn graph_group_id_for_graph(&self, graph_id: usize) -> Option<usize>;
    fn get_group(&self, group_id: GroupId) -> &GraphGroup;
    fn get_group_structure(&self) -> &TiVec<GroupId, GraphGroup>;
    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor<'_>;
    fn take_event_processing_runtime(&mut self) -> Option<EventProcessingRuntime> {
        None
    }
    fn restore_event_processing_runtime(&mut self, _runtime: Option<EventProcessingRuntime>) {}
    fn event_processing_runtime(&self) -> Option<&EventProcessingRuntime> {
        None
    }
    fn event_processing_runtime_mut(&mut self) -> Option<&mut EventProcessingRuntime> {
        None
    }
    fn groups_default_sample_events_by_graph_group(&self) -> bool {
        false
    }

    fn uses_explicit_orientation_sum_only(&self) -> bool {
        false
    }

    // fn get_builder_cache(&self) -> &ParamBuilder<f64>;
}

pub(crate) fn validate_process_runtime_settings(
    settings: &RuntimeSettings,
    explicit_orientation_sum_only: bool,
) -> Result<()> {
    if let Some(parameterization) = settings.sampling.get_parameterization_settings() {
        let alpha = parameterization.sampling_channels.alpha;
        if !alpha.is_finite() || alpha < 0.0 {
            return Err(eyre!(
                "sampling.alpha must be finite and nonnegative; got {alpha}"
            ));
        }
    }
    for (index, level) in settings.stability.levels.iter().enumerate() {
        for (component, tolerance) in [
            ("re", level.ecm_relative_tolerance_for_re),
            ("im", level.ecm_relative_tolerance_for_im),
        ] {
            if !tolerance.is_finite() || tolerance < 0.0 {
                return Err(eyre!(
                    "`runtime.stability.levels[{index}].ecm_relative_tolerance_for_{component}` must be finite and nonnegative; got {tolerance}"
                ));
            }
        }
    }
    if settings.stability.levels.iter().any(|level| {
        level.ecm_relative_tolerance_for_re > 0.0 || level.ecm_relative_tolerance_for_im > 0.0
    }) && (!settings.kinematics.e_cm.is_finite() || settings.kinematics.e_cm <= 0.0)
    {
        return Err(eyre!(
            "E_cm-relative stability requires finite positive `runtime.kinematics.e_cm`"
        ));
    }
    if settings.general.use_ltd {
        return Err(eyre!(
            "`runtime.general.use_ltd = true` is reserved for deferred proper-LTD support; the current evaluation backend is CFF"
        ));
    }

    // The shared process parameter layout always includes M, even when unused.
    if settings.general.numerator_sampling_scale == 0.0 {
        return Err(eyre!(
            "`runtime.general.numerator_sampling_scale` must be nonzero for the auxiliary sampling scale M"
        ));
    }

    if !explicit_orientation_sum_only {
        return Ok(());
    }

    if settings.general.orientation_pat.pat.is_some() {
        return Err(eyre!(
            "`global.generation.explicit_orientation_sum_only = true` already contains the complete orientation sum; `runtime.general.orientation_pat` must be unset"
        ));
    }

    if let SamplingSettings::DiscreteGraphs(discrete_settings) = &settings.sampling
        && discrete_settings.sample_orientations
    {
        return Err(eyre!(
            "`global.generation.explicit_orientation_sum_only = true` does not support runtime individual-orientation Monte Carlo sampling; set `sampling.sample_orientations = false`"
        ));
    }

    Ok(())
}

fn validate_orientation_catalog_group<'a>(
    group_id: GroupId,
    catalogs: impl IntoIterator<Item = (String, Vec<&'a str>)>,
) -> Result<()> {
    let mut catalogs = catalogs.into_iter();
    let Some((master_name, master_catalog)) = catalogs.next() else {
        return Ok(());
    };
    for (graph_name, catalog) in catalogs {
        if catalog != master_catalog {
            let first_difference = master_catalog
                .iter()
                .zip_longest(&catalog)
                .position(|entry| match entry {
                    itertools::EitherOrBoth::Both(left, right) => left != right,
                    itertools::EitherOrBoth::Left(_) | itertools::EitherOrBoth::Right(_) => true,
                })
                .unwrap_or_default();
            return Err(eyre!(
                "Runtime orientation Monte Carlo cannot use graph group {} because graph '{}' and master '{}' have different exact residue-map catalogs ({} versus {} selected maps; first difference at channel {}). Disable `sampling.sample_orientations` so each graph explicitly sums its own complete map catalog.",
                group_id.0,
                graph_name,
                master_name,
                catalog.len(),
                master_catalog.len(),
                first_difference,
            ));
        }
    }
    Ok(())
}

pub(crate) fn validate_group_orientation_catalogs<G: GraphTerm>(
    settings: &RuntimeSettings,
    graph_terms: &[G],
    groups: &TiVec<GroupId, GraphGroup>,
) -> Result<()> {
    let SamplingSettings::DiscreteGraphs(discrete_settings) = &settings.sampling else {
        return Ok(());
    };
    if !discrete_settings.sample_orientations {
        return Ok(());
    }

    for (group_id, group) in groups.iter_enumerated() {
        let master_id = group.master();
        validate_orientation_catalog_group(
            group_id,
            std::iter::once(master_id)
                .chain(group.into_iter().filter(|graph_id| *graph_id != master_id))
                .map(|graph_id| {
                    let graph = &graph_terms[graph_id];
                    (graph.name(), graph.selected_production_orientation_keys())
                }),
        )?;
    }
    Ok(())
}

fn get_global_dimension_if_exists<I: ProcessIntegrandImpl>(integrand: &I) -> Option<usize> {
    if integrand
        .get_settings()
        .sampling
        .get_parameterization_settings()
        .is_none()
    {
        None
    } else {
        Some(
            integrand
                .get_master_graph(GroupId(0))
                .get_graph()
                .get_loop_number()
                * 3,
        )
    }
}

pub trait GraphTerm {
    /// Prepare physical center decisions at the original Arb point. Geometry
    /// is channel-independent for the same selector-accepted cut set; explicit
    /// channel selectors retain their source metadata. Amplitudes own their
    /// separate overlap preparation.
    fn prepare_physical_overlaps(
        &mut self,
        _sample: &MomentumSample<ArbPrec>,
        _context: GraphTermEvaluationContext<'_, '_, ArbPrec>,
    ) -> Result<Option<TiVec<CutGroupId, Option<LUSharedOverlaps<ArbPrec>>>>> {
        Ok(None)
    }

    fn evaluate<T: FloatLike>(
        &mut self,
        sample: &MomentumSample<T>,
        context: GraphTermEvaluationContext<'_, '_, T>,
    ) -> Result<GraphEvaluationResult<T>>;

    fn name(&self) -> String;
    fn orientation_label(&self, orientation_id: usize) -> Option<String>;
    fn sampling_channel_label(
        &self,
        channel_id: SamplingChannelId,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<Option<String>>;

    fn warm_up(&mut self, settings: &RuntimeSettings, model: &Model) -> Result<()>;
    fn get_graph(&self) -> &Graph;
    fn sampling_setup(&self) -> &LmbMultiChannelingSetup;
    fn sampling_setup_mut(&mut self) -> &mut LmbMultiChannelingSetup;
    fn get_num_orientations(&self) -> usize;
    fn production_orientation_keys(&self) -> &[String];
    fn selected_production_orientation_keys(&self) -> Vec<&str>;
    fn selected_lmb_basis_id(
        &self,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<LmbIndex>;
    fn get_tropical_sampler(&self) -> &SampleGenerator<3>;
    fn get_mut_param_builder(&mut self) -> &mut ParamBuilder<f64>;
    fn get_real_mass_vector(&self) -> EdgeVec<Option<F<f64>>>;

    /// Compile the canonical full-frame sampling bridge for this graph.
    /// Process implementations own the graph-specific external signature and
    /// LMB setup; the sampler only supplies the already-resolved kinematics.
    fn compile_sampling_bridge<T: FloatLike>(
        &self,
        parameterization_settings: &ParameterizationSettings,
        settings: &RuntimeSettings,
        external_momenta: &[[T; 4]],
        orientation: Option<usize>,
    ) -> Result<SamplingChannelBridge<T>> {
        let resolved = resolve_sampling_channel_selection(
            &self.sampling_setup().graph.name,
            &parameterization_settings.sampling_channels,
        )?;
        let catalogue = self
            .sampling_setup()
            .sampling_channel_catalogue(&resolved, parameterization_settings)?;
        let programs = catalogue.compile_programs(
            3 * self.get_graph().get_loop_number(),
            &settings.lu_h_function,
        )?;
        self.bind_sampling_bridge(
            &catalogue,
            &programs,
            parameterization_settings,
            settings,
            external_momenta,
            orientation,
        )
    }

    /// Bind the already resolved catalogue and compiled programs to native
    /// graph geometry. Canonical warmup and native component diagnostics reuse
    /// this one definition.
    fn bind_sampling_bridge<T: FloatLike>(
        &self,
        catalogue: &SamplingChannelCatalogue,
        programs: &[SamplingChannelPrograms],
        parameterization_settings: &ParameterizationSettings,
        settings: &RuntimeSettings,
        external_momenta: &[[T; 4]],
        orientation: Option<usize>,
    ) -> Result<SamplingChannelBridge<T>>;

    fn sampling_channel_is_lmb(
        &self,
        channel_id: SamplingChannelId,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<bool>;
    fn sampling_channel_ids(
        &self,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<Vec<SamplingChannelId>>;
}

#[derive(Clone, Copy)]
pub enum EvaluationTarget<'a> {
    Physical(&'a Model),
    Reference(&'a GaussianReferenceFunction),
    #[cfg(test)]
    SamplingLaw(&'a SamplingLawProbe<'a>),
}

/// A fixed normalized density and controlled body failure for exercising the
/// real original-source retry owner. It never reads the current proposal's
/// Jacobian or partition, and is absent from production builds.
#[cfg(test)]
pub struct SamplingLawProbe<'a> {
    reference: &'a LmbMultiChannelingSetup,
    channel_id: SamplingChannelId,
    radial_coordinate: usize,
    retry_lower_half: bool,
    sign_coordinate: Option<usize>,
    alternating_graph_sign: bool,
    calls: [std::sync::atomic::AtomicUsize; 3],
}

#[cfg(test)]
impl SamplingLawProbe<'_> {
    fn evaluate<T: FloatLike>(
        &self,
        sample: &MomentumSample<T>,
        rotation: &Rotation,
        canonical: Option<&DiscreteGraphSample<ArbPrec>>,
        graph_id: usize,
    ) -> Result<GraphEvaluationResult<T>> {
        let canonical =
            canonical.expect("production sampling probe receives the retained canonical row");
        let expected = canonical.sample.materialize::<T>()?.rotate(
            rotation,
            sample.sample.loop_mom_cache_id,
            sample.sample.external_mom_cache_id,
        );
        assert_eq!(
            sample.loop_moms(),
            expected.loop_moms(),
            "physical retry changed the canonical point"
        );
        assert_eq!(sample.external_moms(), expected.external_moms());
        assert_eq!(
            sample.jacobian(),
            expected.jacobian(),
            "physical retry changed the combined canonical factor"
        );
        let precision = T::sampling_precision();
        let index = match precision {
            SamplingPrecision::Double => 0,
            SamplingPrecision::Quad => 1,
            SamplingPrecision::Arb => 2,
            SamplingPrecision::Fixed256 => unreachable!("source precision is not a physical lane"),
        };
        self.calls[index].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let point = sample
            .loop_moms()
            .iter()
            .flat_map(|momentum| {
                let momentum = rotation.inverse_rotate_three(momentum);
                [momentum.px.0, momentum.py.0, momentum.pz.0]
            })
            .collect_vec();
        let mut result = GraphEvaluationResult::zero(sample.zero());
        if let Some(reference) = self.reference.sampling_bridge::<T>()?.channels()
            [self.channel_id.index()]
        .map
        .inverse(&point)?
        {
            let half = sample.one() / sample.one().from_usize(2);
            result.integrand_result.re = if self.retry_lower_half
                && precision == SamplingPrecision::Double
                && F(reference.coordinates[self.radial_coordinate].clone()) < half
            {
                sample.zero() / sample.zero()
            } else {
                F(reference.inverse_jacobian)
            };
            if let Some(axis) = self.sign_coordinate {
                if F(reference.coordinates[axis].clone()) < half {
                    result.integrand_result.re = -result.integrand_result.re;
                }
                result.integrand_result.im = -result.integrand_result.re.clone();
            }
            if self.alternating_graph_sign && graph_id % 2 == 1 {
                result.integrand_result = -result.integrand_result;
            }
        }
        Ok(result)
    }
}

struct EvaluationContext<'a, 'm> {
    target: EvaluationTarget<'a>,
    settings: &'a RuntimeSettings,
    rotation: &'a Rotation,
    evaluation_metadata: &'m mut EvaluationMetaData,
}

pub struct GraphTermEvaluationContext<'a, 'm, T: FloatLike> {
    pub model: &'a Model,
    pub settings: &'a RuntimeSettings,
    pub event_processing_runtime: Option<&'m mut EventProcessingRuntime>,
    pub rotation: &'a Rotation,
    pub evaluation_metadata: &'m mut EvaluationMetaData,
    /// The canonical channel which mapped this point into the parent frame.
    /// Its sampling partition is applied outside the physical graph evaluation.
    pub sampling_channel: Option<SamplingChannelId>,
    pub graph_id: usize,
    pub(crate) prepared_lu_hosts: &'m [PreparedLUHost<T>],
    pub(crate) canonical_sample: Option<&'a DiscreteGraphSample<ArbPrec>>,
    pub(crate) sampling_accuracy_budget: f64,
}

/// Evaluate one graph term using the canonical sampling channel contract.
fn evaluate_graph_term<T: FloatLike, I: ProcessIntegrandImpl>(
    integrand: &mut I,
    graph_id: usize,
    sample: &MomentumSample<T>,
    context: &mut EvaluationContext<'_, '_>,
    sampling_channel: Option<(SamplingChannelId, &[PreparedLUHost<T>])>,
    canonical_sample: Option<&DiscreteGraphSample<ArbPrec>>,
) -> Result<GraphEvaluationResult<T>> {
    if let Some((channel_id, hosts)) = sampling_channel {
        let graph = integrand.get_graph(graph_id);
        graph
            .sampling_setup()
            .sampling_catalogue
            .as_ref()
            .and_then(|catalogue| catalogue.entries.get(channel_id.index()))
            .ok_or_else(|| {
                eyre!(
                    "sampling channel {} is out of range for graph '{}'",
                    channel_id.index(),
                    graph.name()
                )
            })?;
        if hosts.iter().any(|host| {
            host.source.graph_id != graph_id
                || host.source.generating_channel != channel_id
                || host.source.target_channel != channel_id
        }) {
            return Err(eyre!(
                "selected LU host records do not belong to graph-group master {graph_id} and generating channel {channel_id:?}"
            ));
        }
        // Fixed canonical preparation owns the completed parent-frame map and
        // its selected host records. Physical lanes only adopt these records;
        // native map binding is not a prerequisite for evaluating the point.
    }
    // Every actual target and evaluator call counts, including failed attempts
    // and probe rotations. Physical diagnostic/event selection stays separate.
    // Host adoption is sampling work even though it runs inside this body.
    let sampling_before = context.evaluation_metadata.parameterization_time;
    let started = Instant::now();
    let result = (|| -> Result<GraphEvaluationResult<T>> {
        let model = match context.target {
            EvaluationTarget::Physical(model) => model,
            #[cfg(test)]
            EvaluationTarget::SamplingLaw(probe) => {
                return probe.evaluate(sample, context.rotation, canonical_sample, graph_id);
            }
            EvaluationTarget::Reference(reference) => {
                // A reference is defined in the original raw frame, including its
                // nonzero center. Stability rotations must not rotate the target
                // relative to the integration point.
                let raw_momenta = sample
                    .loop_moms()
                    .0
                    .iter()
                    .map(|momentum| context.rotation.inverse_rotate_three(momentum))
                    .collect::<LoopMomenta<_>>();
                let mut value = reference.evaluate(&raw_momenta)?
                    / sample.zero().from_usize(integrand.graph_count());
                if sample.sample.orientation.is_some() {
                    value /= sample
                        .zero()
                        .from_usize(integrand.get_graph(graph_id).get_num_orientations());
                }
                let radius_squared = raw_momenta.0.iter().fold(sample.zero(), |sum, momentum| {
                    sum + momentum.px.square() + momentum.py.square() + momentum.pz.square()
                });
                let mut result = GraphEvaluationResult::zero(sample.zero());
                result.reference_moments = Some(ReferenceMoments {
                    second_moment: value.clone() * radius_squared,
                    jacobian_min: 1.0,
                    jacobian_max: 1.0,
                });
                result.integrand_result = Complex::new_re(value);
                return Ok(result);
            }
        };
        let mut event_processing_runtime = integrand.take_event_processing_runtime();
        let result = {
            let graph_context = GraphTermEvaluationContext {
                model,
                settings: context.settings,
                event_processing_runtime: event_processing_runtime.as_mut(),
                rotation: context.rotation,
                evaluation_metadata: context.evaluation_metadata,
                sampling_channel: sampling_channel.map(|(channel_id, _)| channel_id),
                graph_id,
                prepared_lu_hosts: sampling_channel.map_or(&[], |(_, hosts)| hosts),
                canonical_sample,
                sampling_accuracy_budget: GammaLoopSample::<T>::relative_accuracy_budget(
                    context.settings,
                ),
            };
            integrand
                .get_graph_mut(graph_id)
                .evaluate(sample, graph_context)
        };
        integrand.restore_event_processing_runtime(event_processing_runtime);
        let mut result = result?;
        let graph_group_id = integrand.graph_group_id_for_graph(graph_id);
        for event_group in result.event_groups.iter_mut() {
            for event in event_group.iter_mut() {
                event.cut_info.graph_id = graph_id;
                event.cut_info.graph_group_id = graph_group_id;
            }
        }
        Ok(result)
    })();
    let nested_sampling = context.evaluation_metadata.parameterization_time - sampling_before;
    context.evaluation_metadata.integrand_evaluation_time +=
        started.elapsed().saturating_sub(nested_sampling);
    result
}

fn selected_lmb_basis_for_default_sampling<I: ProcessIntegrandImpl>(
    integrand: &I,
    graph_id: usize,
    use_lmb_basis: bool,
) -> Result<Option<LmbIndex>> {
    if !use_lmb_basis {
        return Ok(None);
    }

    let parameterization_settings = integrand
        .get_settings()
        .sampling
        .get_parameterization_settings()
        .expect("Default LMB-basis sampling requires a parameterization.");
    let group_id = integrand
        .graph_group_id_for_graph(graph_id)
        .map(GroupId)
        .ok_or_else(|| {
            eyre!(
                "Could not determine graph group for graph {} while selecting an LMB basis.",
                graph_id
            )
        })?;

    Ok(Some(
        integrand
            .get_master_graph(group_id)
            .selected_lmb_basis_id(&parameterization_settings)?,
    ))
}

fn evaluate_all_rotations<T: FloatLike, I: ProcessIntegrandImpl>(
    integrand: &mut I,
    target: EvaluationTarget<'_>,
    gammaloop_sample: &GammaLoopSample<T>,
    evaluation_metadata: &mut EvaluationMetaData,
    record_rotated_results: bool,
    canonical_sample: Option<&GammaLoopSample<ArbPrec>>,
) -> Result<(Vec<GraphEvaluationResult<T>>, usize, Vec<RotatedEvaluation>)> {
    let rotations = integrand.get_rotations().cloned().collect_vec();

    let cache = integrand.get_settings().general.enable_cache;

    let mut loop_mom_cache_id = integrand.loop_cache_id();
    let mut external_mom_cache_id = integrand.external_cache_id();

    // rotate the momenta for the stability tests.
    let gammaloop_samples: Vec<_> = rotations
        .iter()
        .map(|rotation| {
            if rotation.is_identity() {
                return gammaloop_sample.clone();
            }
            if cache {
                loop_mom_cache_id += gammaloop_sample.row_count();
                external_mom_cache_id += 1;
            }
            gammaloop_sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id)
        })
        .collect();

    if cache {
        integrand.increment_loop_cache_id(rotations.len() * gammaloop_sample.row_count());
        integrand.increment_external_cache_id(rotations.len());
        // External configurations are shared by subsequent draws; restore the
        // base even if a later physical rotation returns a typed retry error.
        integrand.revert_to_base_external_cache_id();
    }

    let primary_rotation_index = rotations
        .iter()
        .position(Rotation::is_identity)
        .unwrap_or(0);
    let mut evaluation_results: Vec<GraphEvaluationResult<T>> =
        Vec::with_capacity(gammaloop_samples.len());
    for (gammaloop_sample, rotation) in gammaloop_samples.iter().zip(rotations.iter()) {
        debug!("Evaluating rotation: {}", rotation.method);
        let result = evaluate_single(
            integrand,
            target,
            gammaloop_sample,
            rotation,
            evaluation_metadata,
            canonical_sample,
        )?;

        evaluation_results.push(result);
    }

    for result in &evaluation_results {
        evaluation_metadata.event_processing_time += result.event_processing_time;
    }

    let rotated_results = if record_rotated_results {
        rotations
            .iter()
            .zip(evaluation_results.iter())
            .map(|(rotation, result)| RotatedEvaluation {
                rotation: rotation.method.to_string(),
                result: complex_to_f64(&result.integrand_result),
            })
            .collect()
    } else {
        Vec::new()
    };

    Ok((evaluation_results, primary_rotation_index, rotated_results))
}

struct StabilityEvaluationContext<'a, 'm> {
    target: EvaluationTarget<'a>,
    source: &'a EvaluationSource<'a>,
    stability_level: &'a StabilityLevelSetting,
    max_eval: &'a Complex<F<f64>>,
    wgt: F<f64>,
    check_on_norm: bool,
    is_final_level: bool,
    evaluation_metadata: &'m mut EvaluationMetaData,
    record_rotated_results: bool,
    precision_label: &'static str,
    escalate_if_exact_zero: bool,
    check_real: bool,
    check_imag: bool,
}

impl StabilityEvaluationContext<'_, '_> {
    fn ecm_reference_scale<T: FloatLike, I: ProcessIntegrandImpl>(
        &self,
        integrand: &I,
        sample: &GammaLoopSample<T>,
    ) -> Result<Option<F<T>>> {
        if self.stability_level.ecm_relative_tolerance_for_re > 0.0
            || self.stability_level.ecm_relative_tolerance_for_im > 0.0
        {
            let e_cm = integrand.get_settings().kinematics.e_cm;
            if !e_cm.is_finite() || e_cm <= 0.0 {
                return Err(eyre!(
                    "E_cm-relative stability requires finite positive E_cm"
                ));
            }
            let momentum_sample = sample.get_default_sample();
            let missing_measure_dimension = if self.source.is_x_space() {
                0
            } else {
                let loops = momentum_sample.loop_moms().0.len();
                if sample
                    .groups
                    .iter()
                    .flat_map(|(_, rows)| rows)
                    .any(|row| row.sample.loop_moms().0.len() != loops)
                {
                    return Err(eyre!(
                        "E_cm-relative stability cannot compare a raw momentum sum with different spatial dimensions"
                    ));
                }
                i32::try_from(loops)?.checked_mul(3).ok_or_else(|| {
                    eyre!("missing spatial measure dimension exceeds the supported integer range")
                })?
            };
            let scale = match self.target {
                EvaluationTarget::Physical(_) => integrand.stability_reference_scale(
                    sample.groups[0].1[0].graph_id,
                    momentum_sample,
                    missing_measure_dimension,
                )?,
                EvaluationTarget::Reference(_) => {
                    F::<T>::from_f64(integrand.get_settings().kinematics.e_cm)
                        .powi(-missing_measure_dimension)
                }
                #[cfg(test)]
                EvaluationTarget::SamplingLaw(_) => momentum_sample.one(),
            };
            if !scale.0.is_finite() || scale <= momentum_sample.zero() {
                return Err(eyre!(
                    "E_cm-relative stability scale must be finite and positive"
                ));
            }
            Ok(Some(scale))
        } else {
            Ok(None)
        }
    }
}

fn evaluate_stability_level_precise<T: FloatLike, I: ProcessIntegrandImpl>(
    integrand: &mut I,
    context: &mut StabilityEvaluationContext<'_, '_>,
) -> Result<PreciseStabilityLevelResult<T>> {
    let level_start = Instant::now();
    let gammaloop_sample = context
        .source
        .build_gamma_sample::<T, I>(integrand, context.evaluation_metadata)?;
    let ecm_scale = context.ecm_reference_scale(integrand, &gammaloop_sample)?;
    debug!("{} parameterization succeeded", context.precision_label);
    debug!(
        "jacobian: {:+16e}",
        gammaloop_sample.get_default_sample().jacobian()
    );

    let (graph_results, primary_rotation_index, rotated_results) = evaluate_all_rotations(
        integrand,
        context.target,
        &gammaloop_sample,
        context.evaluation_metadata,
        context.record_rotated_results,
        context.source.canonical_sample(),
    )?;
    let threshold_counterterm_failed = context
        .evaluation_metadata
        .threshold_counterterm_error
        .is_some();
    if context.is_final_level
        && let Some(threshold_error) = &context.evaluation_metadata.threshold_counterterm_error
    {
        return Err(eyre!(
            "threshold-counterterm evaluation remained invalid after the final {} stability level: {}",
            context.stability_level.precision,
            threshold_error,
        ));
    }
    let results = graph_results
        .iter()
        .map(|result| result.integrand_result.clone())
        .collect_vec();

    let mut max_eval = complex_from_f64::<T>(context.max_eval);
    if matches!(context.target, EvaluationTarget::Reference(_)) {
        // Integration reports a second observable in Im. Its scale must not
        // alter the Gaussian value's existing independent stability check.
        max_eval.im = max_eval.im.zero();
    }
    let wgt = F::<T>::from_ff64(context.wgt);

    let (average_result, mut estimated_relative_accuracy, mut is_stable, _instability_reason) =
        if context.check_on_norm {
            stability_check_on_norm_components(
                ecm_scale.as_ref(),
                &results,
                context.stability_level,
                max_eval,
                wgt,
                context.is_final_level,
                context.escalate_if_exact_zero,
                context.check_real,
                context.check_imag,
            )
        } else {
            stability_check_components(
                ecm_scale.as_ref(),
                &results,
                context.stability_level,
                max_eval,
                wgt,
                context.is_final_level,
                context.escalate_if_exact_zero,
                context.check_real,
                context.check_imag,
            )
        };

    let mut graph_result = graph_results[primary_rotation_index].clone();
    graph_result.integrand_result = average_result.clone();
    if graph_result.absolute_integrand_result.is_some() {
        // The physical absolute integral is a separate observable: cancellation
        // between different channel points must not hide an unstable absolute
        // contribution. Reuse the configured component/norm criterion, without
        // borrowing the signed integral's maximum-weight escalation scale.
        let absolute_results = graph_results
            .iter()
            .map(|result| {
                result.absolute_integrand_result.clone().ok_or_else(|| {
                    eyre!("sampling-channel absolute contribution missing from a rotation")
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let (absolute_check_real, absolute_check_imag) = (context.check_real, context.check_imag);
        let check = if context.check_on_norm {
            stability_check_on_norm_components::<T>
        } else {
            stability_check_components::<T>
        };
        let (absolute, accuracy, stable, _) = check(
            ecm_scale.as_ref(),
            &absolute_results,
            context.stability_level,
            Complex::new_re(average_result.re.zero()),
            F::<T>::from_ff64(context.wgt),
            context.is_final_level,
            context.escalate_if_exact_zero,
            absolute_check_real,
            absolute_check_imag,
        );
        graph_result.absolute_integrand_result = Some(absolute);
        estimated_relative_accuracy = estimated_relative_accuracy
            .into_iter()
            .chain(accuracy)
            .reduce(|left, right| left.max(right));
        is_stable &= stable;
    }
    if matches!(context.target, EvaluationTarget::Reference(_)) {
        // The Gaussian value and its raw-momentum moment are independent
        // observables. Reuse the scalar stability owner for the moment too;
        // agreement of the value alone cannot certify rotation invariance.
        let moments = graph_results
            .iter()
            .map(|result| {
                result
                    .reference_moments
                    .as_ref()
                    .map(|moments| Complex::new_re(moments.second_moment.clone()))
                    .ok_or_else(|| eyre!("reference acceptance produced no mapped moment"))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut moment_level = *context.stability_level;
        moment_level.required_precision_for_re = moment_level
            .required_precision_for_re
            .min(moment_level.required_precision_for_im);
        moment_level.required_precision_for_im = moment_level.required_precision_for_re;
        // Like its relative precision, the independent moment uses the stricter
        // component allowance: a Re-only allowance does not relax this oracle.
        moment_level.ecm_relative_tolerance_for_re = moment_level
            .ecm_relative_tolerance_for_re
            .min(moment_level.ecm_relative_tolerance_for_im);
        moment_level.ecm_relative_tolerance_for_im = moment_level.ecm_relative_tolerance_for_re;
        // The raw second moment has two additional energy dimensions. Its
        // independent oracle uses E_cm^2 times the reference normalization scale.
        let moment_scale = ecm_scale.as_ref().map(|scale| {
            scale * F::<T>::from_f64(integrand.get_settings().kinematics.e_cm).square()
        });
        if moment_scale
            .as_ref()
            .is_some_and(|scale| !scale.0.is_finite() || scale <= &scale.zero())
        {
            return Err(eyre!(
                "E_cm-relative reference-moment scale must be finite and positive"
            ));
        }
        let (moment, accuracy, stable, _) = stability_check(
            moment_scale.as_ref(),
            &moments,
            &moment_level,
            Complex::new_re(average_result.re.zero()),
            F::<T>::from_ff64(context.wgt),
            context.is_final_level,
            context.escalate_if_exact_zero,
        );
        graph_result
            .reference_moments
            .as_mut()
            .unwrap()
            .second_moment = moment.re;
        estimated_relative_accuracy = estimated_relative_accuracy
            .into_iter()
            .chain(accuracy)
            .reduce(|left, right| left.max(right));
        is_stable &= stable;
    }

    Ok(PreciseStabilityLevelResult {
        result: average_result,
        graph_result,
        stability_level_used: context.stability_level.precision,
        estimated_relative_accuracy,
        sample_count: results.len(),
        total_time: level_start.elapsed(),
        parameterization_jacobian: context
            .source
            .is_x_space()
            .then(|| gammaloop_sample.get_default_sample().one()),
        is_stable: is_stable && !threshold_counterterm_failed,
        rotated_results,
    })
}

/// Result of cache validation checks
#[derive(Debug, Clone)]
pub struct CacheValidationResult {
    pub is_valid: bool,
    pub current_external_cache_id: usize,
    pub base_external_cache_id: usize,
    pub loop_cache_id: usize,
    pub has_rotations: bool,
    pub cache_efficiency: f64,
    pub diagnostics: String,
}

impl std::fmt::Display for CacheValidationResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Cache Validation: {} | Current: {} | Base: {} | Loop: {} | Rotations: {} | Efficiency: {:.1}%",
            if self.is_valid { "✓" } else { "✗" },
            self.current_external_cache_id,
            self.base_external_cache_id,
            self.loop_cache_id,
            if self.has_rotations { "Yes" } else { "No" },
            self.cache_efficiency * 100.0
        )
    }
}

/// Cache usage statistics
#[derive(Debug, Clone)]
pub struct CacheStats {
    pub total_external_increments: usize,
    pub total_loop_increments: usize,
    pub base_configurations: usize,
    pub rotational_variants: usize,
    pub efficiency_ratio: f64,
}

impl std::fmt::Display for CacheStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Cache Stats: {} base configs, {} rotations, {} loop increments, {:.1}% efficiency",
            self.base_configurations,
            self.rotational_variants,
            self.total_loop_increments,
            self.efficiency_ratio * 100.0
        )
    }
}

/// Helper macro for cache debugging
#[macro_export]
macro_rules! debug_cache {
    ($integrand:expr, $msg:expr) => {
        let validation = $integrand.validate_cache_consistency();
        tracing::debug!("{}: {}", $msg, validation);
    };
}

/// Helper macro for cache monitoring with custom conditions
#[macro_export]
macro_rules! monitor_cache {
    ($integrand:expr, $condition:expr, $msg:expr) => {
        if $condition {
            let stats = $integrand.get_cache_stats();
            tracing::info!("{}: {}", $msg, stats);
        }
    };
}

/// Helper macro for missed cache hit detection in debug mode
#[macro_export]
macro_rules! debug_cache_search {
    ($integrand:expr, $msg:expr) => {
        if std::env::var("GAMMALOOP_DEBUG_CACHE").is_ok() {
            let validation = $integrand.validate_cache_consistency();
            if !validation.is_valid {
                tracing::error!(
                    "CACHE CORRUPTION DETECTED at {}: {}",
                    $msg,
                    validation.diagnostics
                );
            } else {
                tracing::debug!("Cache search at {}: {}", $msg, validation);
            }
        }
    };
}

/// Helper macro for cache efficiency warnings
#[macro_export]
macro_rules! warn_cache_efficiency {
    ($integrand:expr, $threshold:expr, $msg:expr) => {
        let stats = $integrand.get_cache_stats();
        if stats.efficiency_ratio < $threshold {
            tracing::warn!(
                "⚠️ Low cache efficiency at {}: {:.1}% (threshold: {:.1}%)",
                $msg,
                stats.efficiency_ratio * 100.0,
                $threshold * 100.0
            );
            tracing::warn!("   Cache stats: {}", stats);
        }
    };
}

fn evaluate_single<T: FloatLike, I: ProcessIntegrandImpl>(
    integrand: &mut I,
    target: EvaluationTarget<'_>,
    gammaloop_sample: &GammaLoopSample<T>,
    rotation: &Rotation,
    evaluation_metadata: &mut EvaluationMetaData,
    canonical_sample: Option<&GammaLoopSample<ArbPrec>>,
) -> Result<GraphEvaluationResult<T>> {
    let settings = integrand.get_settings().clone();
    let zero = gammaloop_sample.get_default_sample().zero();
    let mut context = EvaluationContext {
        target,
        settings: &settings,
        rotation,
        evaluation_metadata,
    };
    let mut result = GraphEvaluationResult::zero(zero.clone());
    let mut absolute = Complex::new_re(zero.clone());
    let mut has_sampling_channels = false;
    for (group_index, (group_id, rows)) in gammaloop_sample.groups.iter().enumerate() {
        let mut grouped_events = crate::observables::GenericEventGroup::default();
        let mut channel_values = BTreeMap::new();
        for (row_index, row) in rows.iter().enumerate() {
            let canonical_row =
                canonical_sample.map(|sample| &sample.groups[group_index].1[row_index]);
            let mut graph_result = evaluate_graph_term(
                integrand,
                row.graph_id,
                &row.sample,
                &mut context,
                row.channel_id
                    .map(|id| (id, row.prepared_lu_hosts.as_slice())),
                canonical_row,
            )?;
            // Preserve tropical's integrand-only compensation. All ordinary
            // map/partition factors reach integrand, events and reference moments
            // once through the same per-row owner, including explicit sums.
            graph_result.integrand_result *= Complex::new_re(row.integrand_prefactor.clone());
            graph_result.apply_sampling_factor(row.sample.jacobian());
            has_sampling_channels |= row.channel_id.is_some();
            // The source owner prepares one master point per channel and clones
            // it across group members. Sum their complete physical bodies here;
            // only distinct channel points contribute separate absolute values.
            *channel_values
                .entry(row.channel_id)
                .or_insert_with(|| Complex::new_re(zero.clone())) +=
                graph_result.integrand_result.clone();
            if group_id.is_some() {
                for mut events in graph_result.event_groups.drain(..) {
                    grouped_events.append(&mut events);
                }
            }
            result.merge_in_place(graph_result);
        }
        if !grouped_events.is_empty() {
            result.event_groups.push(grouped_events);
        }
        for value in channel_values.into_values() {
            absolute += Complex::new(value.re.abs(), value.im.abs());
        }
    }
    result.absolute_integrand_result = Some(if has_sampling_channels {
        absolute
    } else {
        Complex::new(
            result.integrand_result.re.abs(),
            result.integrand_result.im.abs(),
        )
    });
    Ok(result)
}

fn create_grid_for_graph<G: GraphTerm>(
    graph_term: &G,
    settings: &DiscreteGraphSamplingSettings,
    integrator_settings: &IntegratorSettings,
) -> Grid<F<f64>> {
    match &settings.sampling_type {
        DiscreteGraphSamplingType::Default(_) | DiscreteGraphSamplingType::MultiChanneling(_) => {
            let continuous_grid = create_default_continous_grid(graph_term, integrator_settings);

            if settings.sample_orientations {
                let continuous_grids = (0..graph_term.get_num_orientations())
                    .map(|_| Some(continuous_grid.clone()))
                    .collect();

                Grid::Discrete(DiscreteGrid::new(
                    continuous_grids,
                    F(integrator_settings.max_prob_ratio),
                    integrator_settings.train_on_avg,
                ))
            } else {
                continuous_grid
            }
        }
        DiscreteGraphSamplingType::SamplingMultiChanneling(multichanneling_settings) => {
            let continuous_grid = create_default_continous_grid(graph_term, integrator_settings);
            let channel_count = graph_term
                .sampling_channel_ids(&multichanneling_settings.parameterization_settings)
                .map(|channel_ids| channel_ids.len())
                .unwrap_or_else(|error| panic!("cannot build the sampling channel grid: {error}"));
            let lmb_channel_grid = Grid::Discrete(DiscreteGrid::new(
                (0..channel_count)
                    .map(|_| Some(continuous_grid.clone()))
                    .collect_vec(),
                F(integrator_settings.max_prob_ratio),
                integrator_settings.train_on_avg,
            ));

            if settings.sample_orientations {
                Grid::Discrete(DiscreteGrid::new(
                    (0..graph_term.get_num_orientations())
                        .map(|_| Some(lmb_channel_grid.clone()))
                        .collect(),
                    F(integrator_settings.max_prob_ratio),
                    integrator_settings.train_on_avg,
                ))
            } else {
                lmb_channel_grid
            }
        }

        DiscreteGraphSamplingType::TropicalSampling(_) => {
            let dimension = get_n_dim_for_n_loop_momenta(
                &SamplingSettings::DiscreteGraphs(settings.clone()),
                graph_term.get_graph().get_loop_number(),
                Some(graph_term.get_graph().iter_loop_edges().count()),
            );

            let continious_grid = Grid::Continuous(ContinuousGrid::new(
                dimension,
                integrator_settings.n_bins,
                integrator_settings.min_samples_for_update,
                integrator_settings.bin_number_evolution.clone(),
                integrator_settings.train_on_avg,
            ));

            if settings.sample_orientations {
                let continuous_grids = (0..graph_term.get_num_orientations())
                    .map(|_| Some(continious_grid.clone()))
                    .collect();

                Grid::Discrete(DiscreteGrid::new(
                    continuous_grids,
                    F(integrator_settings.max_prob_ratio),
                    integrator_settings.train_on_avg,
                ))
            } else {
                continious_grid
            }
        }
    }
}

fn create_default_continous_grid<G: GraphTerm>(
    graph_term: &G,
    integrator_settings: &IntegratorSettings,
) -> Grid<F<f64>> {
    Grid::Continuous(ContinuousGrid::new(
        graph_term.get_graph().get_loop_number() * 3,
        integrator_settings.n_bins,
        integrator_settings.min_samples_for_update,
        integrator_settings.bin_number_evolution.clone(),
        integrator_settings.train_on_avg,
    ))
}

fn create_grid<I: ProcessIntegrandImpl>(integrand: &I) -> Grid<F<f64>> {
    let settings = integrand.get_settings();
    match &settings.sampling {
        SamplingSettings::Default(_) => Grid::Continuous(ContinuousGrid::new(
            get_global_dimension_if_exists(integrand).unwrap(),
            settings.integrator.n_bins,
            settings.integrator.min_samples_for_update,
            settings.integrator.bin_number_evolution.clone(),
            settings.integrator.train_on_avg,
        )),
        SamplingSettings::MultiChanneling(_) => Grid::Continuous(ContinuousGrid::new(
            get_global_dimension_if_exists(integrand).unwrap(),
            settings.integrator.n_bins,
            settings.integrator.min_samples_for_update,
            settings.integrator.bin_number_evolution.clone(),
            settings.integrator.train_on_avg,
        )),
        SamplingSettings::DiscreteGraphs(discrete_graph_sampling_settings) => {
            Grid::Discrete(DiscreteGrid::new(
                integrand
                    .get_group_masters()
                    .map(|term| {
                        Some(create_grid_for_graph(
                            term,
                            discrete_graph_sampling_settings,
                            &settings.integrator,
                        ))
                    })
                    .collect(),
                F(settings.integrator.max_prob_ratio),
                settings.integrator.train_on_avg,
            ))
        }
    }
}

#[derive(Clone, Copy)]
enum EvaluationSource<'a> {
    XSpace(&'a Sample<F<f64>>),
    Momentum(&'a MomentumSpaceEvaluationInput),
    /// Borrowed from the one original-source preparation outside the physical
    /// retry stack. Retain source kind for reporting, not for native remapping.
    Prepared {
        sample: &'a GammaLoopSample<ArbPrec>,
        original: &'a EvaluationSource<'a>,
    },
}

impl<'a> EvaluationSource<'a> {
    fn canonical_sample(&self) -> Option<&'a GammaLoopSample<ArbPrec>> {
        match self {
            Self::Prepared { sample, .. } => Some(sample),
            _ => None,
        }
    }

    fn is_x_space(&self) -> bool {
        match self {
            Self::XSpace(_) => true,
            Self::Momentum(_) => false,
            Self::Prepared { original, .. } => original.is_x_space(),
        }
    }

    /// Prepare every participating row once at the fixed proposal precision.
    /// Failure ends before physical bodies: choosing another map precision after
    /// inspecting the active point would change the generated proposal law.
    fn prepare_draw<I: ProcessIntegrandImpl>(
        &self,
        integrand: &mut I,
        metadata: &mut EvaluationMetaData,
    ) -> Result<Option<GammaLoopSample<ArbPrec>>> {
        if matches!(self, Self::Prepared { .. }) {
            return Ok(None);
        }
        let started = Instant::now();
        let history = std::mem::take(&mut metadata.radial_root_diagnostics);
        metadata.sampling_proposal_policies.begin_collection();
        let result = (|| {
            let precision = if integrand.get_settings().sampling.uses_sampling_channels()
                && !matches!(self, Self::Momentum(input) if input.channel_id.is_none())
            {
                integrand
                    .get_graph(0)
                    .sampling_setup()
                    .sampling_source
                    .as_ref()
                    .ok_or_else(|| eyre!("fixed sampling source is not initialized; call warm_up"))?
                    .0
            } else {
                SamplingPrecision::Arb
            };
            match (self, precision) {
                (Self::XSpace(sample), SamplingPrecision::Quad) => {
                    parameterize::<f128, I>(sample, integrand, metadata)?.into_canonical(
                        &integrand.get_settings().kinematics.externals,
                        integrand.get_dependent_momenta_constructor(),
                    )
                }
                (Self::Momentum(input), SamplingPrecision::Quad) => {
                    integrand.prepare_sampling_precision::<f128>()?;
                    build_direct_gamma_sample::<f128, I>(integrand, input, metadata)?
                        .into_canonical(
                            &integrand.get_settings().kinematics.externals,
                            integrand.get_dependent_momenta_constructor(),
                        )
                }
                (Self::XSpace(sample), SamplingPrecision::Fixed256) => {
                    parameterize::<SamplingFloat, I>(sample, integrand, metadata)?.into_canonical(
                        &integrand.get_settings().kinematics.externals,
                        integrand.get_dependent_momenta_constructor(),
                    )
                }
                (Self::Momentum(input), SamplingPrecision::Fixed256) => {
                    integrand.prepare_sampling_precision::<SamplingFloat>()?;
                    build_direct_gamma_sample::<SamplingFloat, I>(integrand, input, metadata)?
                        .into_canonical(
                            &integrand.get_settings().kinematics.externals,
                            integrand.get_dependent_momenta_constructor(),
                        )
                }
                (Self::XSpace(sample), SamplingPrecision::Arb) => {
                    parameterize::<ArbPrec, I>(sample, integrand, metadata)
                }
                (Self::Momentum(input), SamplingPrecision::Arb) => {
                    if input.channel_id.is_some() {
                        integrand.prepare_sampling_precision::<ArbPrec>()?;
                    }
                    build_direct_gamma_sample::<ArbPrec, I>(integrand, input, metadata)
                }
                _ => unreachable!("prepared sources and Double source policies are excluded"),
            }
        })();
        // Canonical roots and foreign inverses never consume native physical
        // retry occurrences. Preserve the discrete decisions and inclusive cost.
        metadata.radial_root_diagnostics = history;
        metadata.sampling_proposal_policies.seal();
        let elapsed = started.elapsed();
        metadata.canonical_sampling_preparation_time += elapsed;
        metadata.parameterization_time += elapsed;
        result
            .map(Some)
            .wrap_err("canonical sampling preparation failed at the fixed source precision")
    }

    fn build_gamma_sample<T: FloatLike, I: ProcessIntegrandImpl>(
        &self,
        integrand: &mut I,
        metadata: &mut EvaluationMetaData,
    ) -> Result<GammaLoopSample<T>> {
        // Raw callers outside the stability driver still use the same fixed
        // proposal preparation; the production driver retains its output once.
        let prepared = if !matches!(self, Self::Prepared { .. }) {
            self.prepare_draw(integrand, metadata)?
        } else {
            None
        };
        let started = Instant::now();
        let result = (|| {
            let mut sample = prepared
                .as_ref()
                .or(self.canonical_sample())
                .expect("every original source has a canonical row")
                .materialize::<T>(GammaLoopSample::<T>::relative_accuracy_budget(
                    integrand.get_settings(),
                ))?;
            // Distinct graph/channel rows may share one evaluator cache. Give
            // every completed point its own identity, including base identities
            // used by rotation-aware caches; retries start with a fresh range.
            let first_id = integrand.loop_cache_id();
            for (index, row) in sample
                .groups
                .iter_mut()
                .flat_map(|(_, rows)| rows)
                .enumerate()
            {
                row.sample.sample.loop_mom_cache_id = first_id + index;
                row.sample.sample.loop_mom_base_cache_id = first_id + index;
            }
            if integrand.get_settings().general.enable_cache {
                integrand.increment_loop_cache_id(sample.row_count());
            }
            Ok(sample)
        })();
        metadata.parameterization_time += started.elapsed();
        result
    }

    fn loop_norm_sum<I: ProcessIntegrandImpl>(
        &self,
        integrand: &mut I,
        metadata: &mut EvaluationMetaData,
    ) -> Result<F<f64>> {
        if let Self::Momentum(input) = self {
            return Ok(sum_loop_norms(input.loop_momenta.iter()));
        }
        // Inspect only canonical momenta: a norm prepass must not reject an
        // otherwise usable Arb draw because its Jacobian is outside f64 range.
        let prepared = if self.canonical_sample().is_none() {
            self.prepare_draw(integrand, metadata)?
        } else {
            None
        };
        let sample = self
            .canonical_sample()
            .or(prepared.as_ref())
            .expect("mapped source has prepared rows");
        let started = Instant::now();
        let zero = sample.get_default_sample().zero();
        let largest = sample
            .groups
            .iter()
            .flat_map(|(_, rows)| rows)
            .map(|row| {
                row.sample
                    .loop_moms()
                    .0
                    .iter()
                    .fold(zero.clone(), |sum, momentum| sum + momentum.norm())
            })
            .fold(zero.clone(), |a, b| if a > b { a } else { b });
        metadata.parameterization_time += started.elapsed();
        Ok(largest.into_ff64())
    }

    fn debug_sample<I: ProcessIntegrandImpl>(
        &self,
        integrand: &mut I,
        metadata: &mut EvaluationMetaData,
    ) -> Result<GammaLoopSample<f64>> {
        self.build_gamma_sample::<f64, I>(integrand, metadata)
    }
}

fn sum_loop_norms<'a>(loop_momenta: impl Iterator<Item = &'a ThreeMomentum<F<f64>>>) -> F<f64> {
    loop_momenta.fold(F(0.0), |acc, momentum| acc + momentum.norm())
}

fn build_direct_gamma_sample<T: FloatLike, I: ProcessIntegrandImpl>(
    integrand: &mut I,
    input: &MomentumSpaceEvaluationInput,
    metadata: &mut EvaluationMetaData,
) -> Result<GammaLoopSample<T>> {
    if integrand.uses_explicit_orientation_sum_only() && input.orientation.is_some() {
        return Err(eyre!(
            "`global.generation.explicit_orientation_sum_only = true` represents all orientations as a single summed contribution, so explicit orientation selection is not supported"
        ));
    }

    let expected_loop_count = if let Some(graph_id) = input.graph_id {
        let group_id = integrand
            .get_group_structure()
            .iter_enumerated()
            .find_map(|(group_id, group)| group.into_iter().contains(&graph_id).then_some(group_id))
            .ok_or_else(|| eyre!("Unknown graph '{}' in momentum-space evaluation.", graph_id))?;
        integrand
            .get_master_graph(group_id)
            .get_graph()
            .get_loop_number()
    } else if let Some(group_id) = input.group_id {
        if group_id.0 >= integrand.get_group_structure().len() {
            return Err(eyre!(
                "Unknown graph group '{}' in momentum-space evaluation.",
                group_id.0
            ));
        }
        integrand
            .get_master_graph(group_id)
            .get_graph()
            .get_loop_number()
    } else {
        integrand
            .get_group_masters()
            .next()
            .map(|graph| graph.get_graph().get_loop_number())
            .ok_or_else(|| eyre!("Cannot evaluate an integrand with no graph terms."))?
    };

    if input.loop_momenta.len() != expected_loop_count {
        return Err(eyre!(
            "Expected {} loop momenta in momentum-space evaluation, got {}.",
            expected_loop_count,
            input.loop_momenta.len()
        ));
    }

    // A direct input already denotes a binary64 point. Preserve those exact
    // tokens in every lane; the decimal conversion policy for settings is unchanged.
    let loop_momenta = input
        .loop_momenta
        .iter()
        .map(|momentum| {
            ThreeMomentum::new(
                F(T::from_f64_exact_binary(momentum.px.0)),
                F(T::from_f64_exact_binary(momentum.py.0)),
                F(T::from_f64_exact_binary(momentum.pz.0)),
            )
        })
        .collect::<LoopMomenta<F<T>>>();
    let sample = MomentumSample::new(
        loop_momenta,
        integrand.loop_cache_id(),
        &integrand.get_settings().kinematics.externals,
        integrand.get_current_external_cache_id(),
        F::<T>::from_f64(1.0),
        integrand.get_dependent_momenta_constructor(),
        input.orientation,
    )?;

    if let Some(graph_id) = input.graph_id {
        if input.group_id.is_some() || input.channel_id.is_some() {
            return Err(eyre!(
                "Explicit graph selection is mutually exclusive with discrete graph/channel selections in momentum-space evaluation."
            ));
        }
        return Ok(GammaLoopSample {
            groups: vec![(
                None,
                vec![DiscreteGraphSample {
                    graph_id,
                    channel_id: None,
                    prepared_lu_hosts: vec![],
                    physical_overlaps: None,
                    integrand_prefactor: sample.one(),
                    sample,
                }],
            )],
        });
    }

    match &integrand.get_settings().sampling {
        SamplingSettings::Default(_) | SamplingSettings::MultiChanneling(_) => {
            if input.group_id.is_some() || input.channel_id.is_some() {
                return Err(eyre!(
                    "Discrete graph/channel selections are not supported for this sampling mode."
                ));
            }

            GammaLoopSample::from_default(integrand, sample, false, None)
        }
        SamplingSettings::DiscreteGraphs(settings) => {
            let Some(group_id) = input.group_id else {
                if input.orientation.is_some() || input.channel_id.is_some() {
                    return Err(eyre!(
                        "Explicit orientation or channel selections require selecting a graph group in momentum-space evaluation."
                    ));
                }
                return GammaLoopSample::from_default(integrand, sample, false, None);
            };
            match &settings.sampling_type {
                DiscreteGraphSamplingType::Default(_) => {
                    if input.channel_id.is_some() {
                        return Err(eyre!(
                            "Channel selection is not available for this discrete-graph sampling mode."
                        ));
                    }
                    GammaLoopSample::from_default(integrand, sample, false, Some(group_id))
                }
                DiscreteGraphSamplingType::MultiChanneling(_) => {
                    if input.channel_id.is_some() {
                        return Err(eyre!(
                            "Channel selection is not available for this discrete-graph sampling mode."
                        ));
                    }
                    // Direct momentum input already contains the desired
                    // parent-frame point; only unit-cube samples replay maps.
                    GammaLoopSample::from_default(integrand, sample, false, Some(group_id))
                }
                DiscreteGraphSamplingType::TropicalSampling(_) => {
                    if input.channel_id.is_some() {
                        return Err(eyre!(
                            "Channel selection is not available for tropical discrete-graph sampling."
                        ));
                    }
                    GammaLoopSample::from_tropical(integrand, group_id, sample)
                }
                DiscreteGraphSamplingType::SamplingMultiChanneling(_) => {
                    let Some(channel_id) = input.channel_id else {
                        // An unselected raw point evaluates the complete physical
                        // graph without a proposal partition or artificial channel.
                        return GammaLoopSample::from_default(
                            integrand,
                            sample,
                            false,
                            Some(group_id),
                        );
                    };
                    let graph = integrand.get_master_graph(group_id);
                    let bridge = graph.sampling_setup().sampling_bridge::<T>()?;
                    let mut contexts = SamplingChannelRuntimeContexts::for_draw(
                        bridge.channels().len(),
                        integrand.get_group(group_id).master(),
                        channel_id,
                        metadata,
                    );
                    let mapped = bridge.inverse_with_runtime_contexts(
                        channel_id,
                        &sample
                            .loop_moms()
                            .0
                            .iter()
                            .flat_map(|momentum| {
                                [
                                    momentum.px.0.clone(),
                                    momentum.py.0.clone(),
                                    momentum.pz.0.clone(),
                                ]
                            })
                            .collect::<Vec<_>>(),
                        &mut contexts,
                    )?
                    .ok_or_else(|| eyre!(
                        "raw momentum point is outside sampling channel {}; select a full-support sibling for direct-momentum evaluation",
                        channel_id.0,
                    ))?;
                    let partition_weight =
                        mapped.partition.weight(channel_id.0).ok_or_else(|| {
                            eyre!(
                                "sampling channel partition has no weight for channel {}",
                                channel_id.0
                            )
                        })?;
                    if !partition_weight.is_finite() || partition_weight <= partition_weight.zero()
                    {
                        return Err(eyre!(
                            "sampling channel partition has invalid weight {partition_weight}"
                        ));
                    }
                    let mut sample = sample;
                    // Direct inputs have no parameterization Jacobian. Retain
                    // only the selected partition in the same per-row factor.
                    sample.sample.jacobian = F(partition_weight);
                    GammaLoopSample::from_selected(
                        integrand,
                        group_id,
                        channel_id,
                        sample,
                        mapped.prepared_lu_hosts,
                    )
                }
            }
        }
    }
}

fn log_rotated_samples<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    gammaloop_sample: &GammaLoopSample<f64>,
    rotated_results: &[RotatedEvaluation],
) {
    let mut loop_mom_cache_id = integrand.loop_cache_id();
    let mut external_mom_cache_id = integrand.external_cache_id();
    let mut shift = 0;

    let rotated_samples: Vec<_> = integrand
        .get_rotations()
        .map(|rotation| {
            if rotation.is_identity() {
                return gammaloop_sample.clone();
            }
            loop_mom_cache_id += gammaloop_sample.row_count();
            shift += 1;
            external_mom_cache_id += 1;
            gammaloop_sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id)
        })
        .collect();
    integrand.increment_external_cache_id(shift);
    integrand.increment_loop_cache_id(shift * gammaloop_sample.row_count());

    for (sample, result) in rotated_samples.iter().zip(rotated_results.iter()) {
        let default_sample = sample.get_default_sample();
        debug!(
            "loop_moms: {}, external_moms: {}",
            format!("{}", default_sample.loop_moms()).blue(),
            format!("{:?}", default_sample.external_moms()).blue()
        );

        debug!(
            "result of current level: {}",
            format!("{:16e}", result.result).blue()
        );
    }
}

fn evaluate_from_source_precise<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    target: EvaluationTarget<'_>,
    source: EvaluationSource<'_>,
    wgt: F<f64>,
    use_arb_prec: bool,
    max_eval: Complex<F<f64>>,
) -> Result<PreciseEvaluationResult> {
    let stability = &integrand.get_settings().stability;
    if matches!(target, EvaluationTarget::Physical(_))
        && stability.integrated_energy_dimension.is_none()
        && stability.levels.iter().any(|level| {
            level.ecm_relative_tolerance_for_re > 0.0 || level.ecm_relative_tolerance_for_im > 0.0
        })
    {
        return Err(eyre!(
            "E_cm-relative physical stability requires `stability.integrated_energy_dimension` (after flux normalization, before output-unit conversion)"
        ));
    }
    let start_eval = std::time::Instant::now();
    let mut escalate_if_exact_zero = integrand.get_settings().stability.escalate_if_exact_zero;
    if escalate_if_exact_zero
        && integrand
            .get_settings()
            .selectors
            .values()
            .any(|selector| selector.active)
    {
        warn_selectors_disable_zero_once();
        escalate_if_exact_zero = false;
    }
    let mut evaluation_metadata = EvaluationMetaData::new_empty();
    if matches!(target, EvaluationTarget::Reference(_))
        && matches!(&integrand.get_settings().sampling, SamplingSettings::DiscreteGraphs(settings)
            if matches!(settings.sampling_type, DiscreteGraphSamplingType::TropicalSampling(_)))
    {
        return Err(eyre!(
            "reference acceptance does not support tropical sampling"
        ));
    }
    let mut anchor = source.prepare_draw(integrand, &mut evaluation_metadata)?;
    if let (Some(sample), EvaluationTarget::Physical(model)) = (&mut anchor, target) {
        sample.prepare_physical_overlaps(integrand, model, &mut evaluation_metadata)?;
    }
    let original_source = source;
    let source = match &anchor {
        Some(sample) => EvaluationSource::Prepared {
            sample,
            original: &original_source,
        },
        None => original_source,
    };
    let (stability_iterator, loop_momenta_escalation) =
        stability_iterator_for_source(integrand, &source, use_arb_prec, &mut evaluation_metadata);

    let mut final_result = None;

    let mut stability_results = Vec::with_capacity(stability_iterator.len());
    let mut sampling_failures = Vec::new();
    let total_levels = stability_iterator.len();
    for (level_index, stability_level) in stability_iterator.into_iter().enumerate() {
        evaluation_metadata.clear_threshold_counterterm_error();
        let is_final_level = level_index + 1 == total_levels;
        let record_rotated_results = integrand
            .get_settings()
            .stability
            .recording
            .map(|recording| recording.record_rotated_results)
            .unwrap_or(false);
        let mut context = StabilityEvaluationContext {
            target,
            source: &source,
            stability_level: &stability_level,
            max_eval: &max_eval,
            wgt,
            check_on_norm: integrand.get_settings().stability.check_on_norm,
            is_final_level,
            evaluation_metadata: &mut evaluation_metadata,
            record_rotated_results,
            precision_label: match stability_level.precision {
                Precision::Double => "f64",
                Precision::Quad => "f128",
                Precision::Arb => "ArbPrec",
            },
            escalate_if_exact_zero,
            check_real: !matches!(target, EvaluationTarget::Physical(_))
                || !matches!(
                    integrand.get_settings().integrator.integrated_phase,
                    IntegratedPhase::Imag
                ),
            check_imag: !matches!(target, EvaluationTarget::Physical(_))
                || !matches!(
                    integrand.get_settings().integrator.integrated_phase,
                    IntegratedPhase::Real
                ),
        };
        let level_start = Instant::now();
        let mut is_stable = false;
        let mut rotated_results = Vec::new();
        macro_rules! evaluate_native_level {
            ($scalar:ty, $variant:ident) => {
                evaluate_stability_level_precise::<$scalar, I>(integrand, &mut context).map(
                    |mut result| {
                        is_stable = result.is_stable;
                        stability_results.push(StabilityResult {
                            precision: result.stability_level_used,
                            estimated_relative_accuracy: result
                                .estimated_relative_accuracy
                                .as_ref()
                                .map(F::into_ff64),
                            status: StabilityStatus::from_sample_count(
                                result.sample_count,
                                is_stable,
                            ),
                            total_time: result.total_time,
                        });
                        debug!(
                            "level: {}. result: {}",
                            format!("{}", result.stability_level_used).green(),
                            format!("{:16e}", result.result).blue()
                        );
                        rotated_results = std::mem::take(&mut result.rotated_results);
                        PreciseEvaluationResult::$variant(finalize_precise_evaluation_result(
                            result,
                            wgt,
                            evaluation_metadata.clone(),
                        ))
                    },
                )
            };
        }
        let result_of_level = match stability_level.precision {
            Precision::Double => evaluate_native_level!(f64, Double),
            Precision::Quad => evaluate_native_level!(f128, Quad),
            Precision::Arb => evaluate_native_level!(ArbPrec, Arb),
        };
        let result_of_level = match result_of_level {
            Ok(result) => result,
            Err(error)
                if error
                    .downcast_ref::<sampling_maps::SamplingEvaluationError>()
                    .is_some() =>
            {
                stability_results.push(StabilityResult {
                    precision: stability_level.precision,
                    estimated_relative_accuracy: None,
                    status: StabilityStatus::Unstable(0),
                    total_time: level_start.elapsed(),
                });
                crate::debug_tags!(#sampling;
                    stage = "native_sampling_retry", precision = %stability_level.precision,
                    error = %error, final_level = is_final_level,
                    "sampling reconstruction failed at this precision"
                );
                // A partial rotation/channel pass may have populated numerical
                // caches. Start the next reconstruction with fresh point IDs;
                // dropped GraphEvaluationResults contain all its uncommitted events.
                integrand.increment_loop_cache_id(integrand.get_rotations().count() + 1);
                integrand.revert_to_base_external_cache_id();
                sampling_failures.push(format!("{}: {error:#}", stability_level.precision));
                if is_final_level {
                    return Err(error.wrap_err(format!(
                        "sampling reconstruction failed; attempts [{}]; numerical errors [{}]",
                        stability_results
                            .iter()
                            .map(|result| format!("{}: {:?}", result.precision, result.status))
                            .join("; "),
                        sampling_failures.join("; ")
                    )));
                }
                continue;
            }
            Err(error) => return Err(error),
        };

        if matches!(target, EvaluationTarget::Reference(_)) && is_final_level && !is_stable {
            return Err(eyre!(
                "reference value or moment remained unstable after the final {} precision level",
                stability_level.precision
            ));
        }
        final_result = Some(result_of_level);

        if is_stable {
            break;
        } else {
            debug!("unstable at level: {}", stability_level.precision);
            if let Ok(gammaloop_sample) = source.debug_sample(integrand, &mut evaluation_metadata) {
                log_rotated_samples(integrand, &gammaloop_sample, &rotated_results);
            } else {
                debug!("failed to reconstruct sample for instability logging");
            }
        }
    }

    let mut result = final_result.ok_or_else(|| eyre!("no stability level was evaluated"))?;
    let metadata = match &mut result {
        PreciseEvaluationResult::Double(result) => &mut result.evaluation_metadata,
        PreciseEvaluationResult::Quad(result) => &mut result.evaluation_metadata,
        PreciseEvaluationResult::Arb(result) => &mut result.evaluation_metadata,
    };
    metadata.total_timing = start_eval.elapsed();
    // A last unstable-level debug replay occurs after its result was cloned.
    // Copy the sole source owner's accumulated counters, never one lane's time.
    metadata.parameterization_time = evaluation_metadata.parameterization_time;
    metadata.canonical_sampling_preparation_time =
        evaluation_metadata.canonical_sampling_preparation_time;
    metadata.integrand_evaluation_time = evaluation_metadata.integrand_evaluation_time;
    metadata.loop_momenta_escalation = loop_momenta_escalation;
    metadata.stability_results = stability_results;
    Ok(result)
}

fn evaluate_from_source<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    source: EvaluationSource<'_>,
    wgt: F<f64>,
    use_arb_prec: bool,
    max_eval: Complex<F<f64>>,
) -> Result<EvaluationResult> {
    evaluate_from_source_precise(
        integrand,
        EvaluationTarget::Physical(model),
        source,
        wgt,
        use_arb_prec,
        max_eval,
    )?
    .try_into_f64()
}

fn stability_iterator_for_source<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    source: &EvaluationSource<'_>,
    use_arb_prec: bool,
    metadata: &mut EvaluationMetaData,
) -> (
    Vec<StabilityLevelSetting>,
    Option<LoopMomentaEscalationMetrics>,
) {
    let mut stability_iterator =
        create_stability_iterator(&integrand.get_settings().stability, use_arb_prec);
    let escalation_factor = integrand
        .get_settings()
        .stability
        .loop_momenta_norm_escalation_factor;
    let record_loop_momenta_escalation = integrand
        .get_settings()
        .stability
        .recording
        .map(|recording| recording.record_loop_momenta_escalation)
        .unwrap_or(false);
    let mut loop_momenta_escalation = None;
    if escalation_factor > 0.0
        && stability_iterator.len() > 1
        && let Ok(sum_norm) = source.loop_norm_sum(integrand, metadata)
    {
        let threshold =
            F::<f64>::from_f64(escalation_factor * integrand.get_settings().kinematics.e_cm);
        if record_loop_momenta_escalation {
            loop_momenta_escalation = Some(LoopMomentaEscalationMetrics {
                sum_norm: sum_norm.0,
                threshold: threshold.0,
            });
        }
        if sum_norm > threshold {
            let escalated_level_index = stability_iterator
                .iter()
                .position(|level| level.precision == Precision::Quad)
                .or_else(|| stability_iterator.len().checked_sub(1));
            if let Some(level_index) = escalated_level_index {
                stability_iterator = stability_iterator[level_index..].to_vec();
            }
        }
    }

    (stability_iterator, loop_momenta_escalation)
}

fn finalize_precise_evaluation_result<T: FloatLike>(
    result: PreciseStabilityLevelResult<T>,
    integrator_weight: F<f64>,
    mut evaluation_metadata: EvaluationMetaData,
) -> GenericEvaluationResult<T> {
    let re_is_nan = result.result.re.is_nan()
        || result.result.re.is_infinite()
        || result
            .graph_result
            .absolute_integrand_result
            .as_ref()
            .is_some_and(|value| value.re.is_nan() || value.re.is_infinite());
    let im_is_nan = result.result.im.is_nan()
        || result.result.im.is_infinite()
        || result
            .graph_result
            .absolute_integrand_result
            .as_ref()
            .is_some_and(|value| value.im.is_nan() || value.im.is_infinite());
    if re_is_nan || im_is_nan {
        warn!(
            stage = "process_final_nonfinite_sample",
            result = %result.result,
            re_is_nan,
            im_is_nan,
            "process evaluation is nonfinite"
        );
    }
    evaluation_metadata.generated_event_count = result.graph_result.generated_event_count;
    evaluation_metadata.accepted_event_count = result.graph_result.accepted_event_count;
    evaluation_metadata.relative_instability_error = Complex::new_zero();
    evaluation_metadata.is_nan = re_is_nan || im_is_nan;
    let nanless_result = if re_is_nan && !im_is_nan {
        Complex::new(result.result.re.zero(), result.result.im)
    } else if im_is_nan && !re_is_nan {
        Complex::new(result.result.re, result.result.im.zero())
    } else if re_is_nan && im_is_nan {
        Complex::new(result.result.re.zero(), result.result.im.zero())
    } else {
        result.result
    };

    let mut event_groups = result.graph_result.event_groups;
    let integrator_weight = F::<T>::from_ff64(integrator_weight);
    let parameterization_jacobian = result.parameterization_jacobian;
    let full_factor = full_event_multiplicative_factor_precise(
        parameterization_jacobian.clone(),
        integrator_weight.clone(),
    );
    apply_full_event_multiplicative_factor_precise(&mut event_groups, &full_factor);

    GenericEvaluationResult {
        reference_moments: result.graph_result.reference_moments,
        integrand_result: nanless_result,
        absolute_integrand_result: result.graph_result.absolute_integrand_result.map(|value| {
            Complex::new(
                if re_is_nan { value.re.zero() } else { value.re },
                if im_is_nan { value.im.zero() } else { value.im },
            )
        }),
        parameterization_jacobian,
        integrator_weight,
        event_groups,
        evaluation_metadata,
    }
}

fn warn_selectors_disable_zero_once() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        warn!(
            "disabling `stability.escalate_if_exact_zero` during evaluation because selectors can legitimately zero the event weight"
        );
    });
}

fn evaluate_sample<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    sample: &Sample<F<f64>>,
    wgt: F<f64>,
    _iter: usize,
    use_arb_prec: bool,
    max_eval: Complex<F<f64>>,
) -> Result<EvaluationResult> {
    evaluate_from_source(
        integrand,
        model,
        EvaluationSource::XSpace(sample),
        wgt,
        use_arb_prec,
        max_eval,
    )
}

fn evaluate_reference_sample<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    sample: &Sample<F<f64>>,
    reference: &GaussianReferenceFunction,
) -> Result<ReferenceSampleEvaluation> {
    let result = evaluate_from_source_precise(
        integrand,
        EvaluationTarget::Reference(reference),
        EvaluationSource::XSpace(sample),
        sample.get_weight(),
        false,
        Complex::new_zero(),
    )?;
    // Select precision once in the common loop, then narrow its value and
    // moment together at the existing f64 acceptance reporting boundary.
    macro_rules! report {
        ($result:expr) => {{
            let mut result = $result;
            if result.evaluation_metadata.is_nan {
                return Err(eyre!(
                    "reference acceptance cannot report a nonfinite finalized value"
                ));
            }
            let moments = result
                .reference_moments
                .take()
                .ok_or_else(|| eyre!("reference acceptance produced no mapped moment"))?
                .try_into_f64()?;
            Ok(ReferenceSampleEvaluation {
                evaluation: result.try_into_f64()?,
                moments,
            })
        }};
    }
    match result {
        PreciseEvaluationResult::Double(result) => report!(result),
        PreciseEvaluationResult::Quad(result) => report!(result),
        PreciseEvaluationResult::Arb(result) => report!(result),
    }
}

fn evaluate_sample_precise<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    sample: &Sample<F<f64>>,
    wgt: F<f64>,
    use_arb_prec: bool,
    max_eval: Complex<F<f64>>,
) -> Result<crate::integrands::evaluation::PreciseEvaluationResult> {
    evaluate_from_source_precise(
        integrand,
        EvaluationTarget::Physical(model),
        EvaluationSource::XSpace(sample),
        wgt,
        use_arb_prec,
        max_eval,
    )
}

fn evaluate_momentum_configuration<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    input: &MomentumSpaceEvaluationInput,
    wgt: F<f64>,
    use_arb_prec: bool,
    max_eval: Complex<F<f64>>,
) -> Result<EvaluationResult> {
    evaluate_from_source(
        integrand,
        model,
        EvaluationSource::Momentum(input),
        wgt,
        use_arb_prec,
        max_eval,
    )
}

fn evaluate_momentum_configuration_precise<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    input: &MomentumSpaceEvaluationInput,
    wgt: F<f64>,
    use_arb_prec: bool,
    max_eval: Complex<F<f64>>,
) -> Result<crate::integrands::evaluation::PreciseEvaluationResult> {
    evaluate_from_source_precise(
        integrand,
        EvaluationTarget::Physical(model),
        EvaluationSource::Momentum(input),
        wgt,
        use_arb_prec,
        max_eval,
    )
}

#[cfg(test)]
pub(crate) mod tests {
    use super::{
        GraphTerm, LmbMultiChannelingSetup, RuntimeCache, SamplingChannelCompileContext,
        SamplingChannelId, filtered_orientation_count, resolve_sampling_channel_selection,
        resolve_visible_orientation_id, validate_orientation_catalog_group,
        validate_process_runtime_settings,
    };
    use crate::cff::expression::OrientationID;
    use crate::{
        dot,
        graph::{
            FeynmanGraph, Graph, GroupId, LMBext, LmbIndex, LoopMomentumBasis,
            parse::from_dot::IntoGraph,
        },
        initialisation::test_initialise,
        momentum::{
            ThreeMomentum,
            sample::{BareMomentumSample, ExternalFourMomenta, LoopMomenta, MomentumSample},
            signature::LoopExtSignature,
        },
        settings::{
            RuntimeSettings,
            global::OrientationPattern,
            runtime::{
                DiscreteGraphSamplingSettings, DiscreteGraphSamplingType, HFunctionSettings,
                MultiChannelingSettings, ParameterizationSettings, SamplingChannelDefinition,
                SamplingChannelSelection, SamplingSettings,
            },
        },
        utils::F,
    };
    use linnet::half_edge::{
        involution::{EdgeIndex, EdgeVec, Orientation},
        subgraph::{ModifySubSet, SubSetLike, subset::SubSet},
    };
    use std::sync::OnceLock;
    use typed_index_collections::TiVec;

    /// Use the existing generated scalar cut fixture for a fixed-source body
    /// retry. Its normalized q target makes a missing/doubled J*w observable.
    pub(crate) fn check_fixed_quad_source_transport(
        runtime: &mut super::ProcessIntegrand,
        model: &crate::model::Model,
        source: &symbolica::numerical_integration::Sample<F<f64>>,
    ) -> color_eyre::Result<()> {
        use super::{
            EvaluationSource, EvaluationTarget, GraphTerm, ProcessIntegrand, ProcessIntegrandImpl,
            SamplingLawProbe, evaluate_from_source_precise,
        };
        use crate::{
            integrands::evaluation::EvaluationMetaData,
            settings::runtime::{Precision, StabilityLevelSetting},
            utils::{ArbPrec, QuadFloat, SamplingFloat, SamplingPrecision},
        };
        use spenso::algebra::complex::Complex;
        use std::sync::atomic::{AtomicUsize, Ordering};
        let original_sampling = runtime.get_settings().sampling.clone();
        for fixed256 in [false, true] {
            runtime.get_mut_settings().sampling = original_sampling.clone();
            if fixed256 {
                let mut parser: crate::settings::runtime::SamplingSettingsParser =
                    toml::from_str(&toml::to_string(&original_sampling)?)?;
                parser.sampling_channel_weight =
                    crate::settings::runtime::SamplingChannelWeight::SingularityProxy;
                for definitions in parser.channel_definitions.values_mut() {
                    for definition in definitions.values_mut() {
                        definition.singularity_proxy = Some("1".into());
                    }
                }
                runtime.get_mut_settings().sampling = toml::from_str(&toml::to_string(&parser)?)?;
            }
            let expected_source = if fixed256 {
                SamplingPrecision::Fixed256
            } else {
                SamplingPrecision::Quad
            };
            let settings = runtime.get_mut_settings();
            settings.stability.levels = vec![
                StabilityLevelSetting::default_double(),
                StabilityLevelSetting::default_quad(),
                StabilityLevelSetting::default_arb(),
            ];
            settings.stability.rotation_axis.clear();
            settings.stability.escalate_if_exact_zero = false;
            settings.stability.loop_momenta_norm_escalation_factor = 0.0;
            runtime.warm_up(model)?;
            let ProcessIntegrand::CrossSection(integrand) = runtime else {
                unreachable!("generated cut fixture")
            };
            let policy = integrand
                .get_graph(0)
                .sampling_setup()
                .sampling_source
                .as_ref()
                .copied()
                .unwrap();
            assert_eq!(policy.0, expected_source);
            let saved_catalogue = integrand
                .get_graph(0)
                .sampling_setup()
                .sampling_catalogue
                .as_ref()
                .unwrap()
                .clone();
            integrand
                .get_graph_mut(0)
                .sampling_setup_mut()
                .sampling_catalogue
                .as_mut()
                .unwrap()
                .entries
                .push(super::sampling_selection::SamplingCatalogueEntry::Surface {
                    edges: vec![1, 2],
                    parent_lmb: vec![1],
                });
            assert_eq!(
                integrand.sampling_source_policy()?.0,
                SamplingPrecision::Fixed256
            );
            integrand
                .get_graph_mut(0)
                .sampling_setup_mut()
                .sampling_catalogue
                .set(saved_catalogue);
            let saved_externals = integrand.settings.kinematics.externals.clone();
            let mut external_values = saved_externals.get_dependent_externals::<ArbPrec>(
                integrand.get_dependent_momenta_constructor(),
            )?;
            external_values[crate::momentum::sample::ExternalIndex(0)]
                .spatial
                .px = F::<ArbPrec>::default().from_usize(2).powi(2000);
            let crate::settings::runtime::kinematic::Externals::Constant { arb_cache, .. } =
                &mut integrand.settings.kinematics.externals;
            arb_cache.set(external_values);
            assert_eq!(
                integrand.sampling_source_policy()?.0,
                SamplingPrecision::Fixed256
            );
            integrand.settings.kinematics.externals = saved_externals;
            assert_eq!(integrand.sampling_source_policy()?, policy);
            let mut metadata = EvaluationMetaData::new_empty();
            let expected = if fixed256 {
                super::parameterize::<SamplingFloat, _>(source, integrand, &mut metadata)?
                    .into_canonical(
                        &integrand.get_settings().kinematics.externals,
                        integrand.get_dependent_momenta_constructor(),
                    )?
            } else {
                super::parameterize::<QuadFloat, _>(source, integrand, &mut metadata)?
                    .into_canonical(
                        &integrand.get_settings().kinematics.externals,
                        integrand.get_dependent_momenta_constructor(),
                    )?
            };
            let source_kind = EvaluationSource::XSpace(source);
            let anchor = source_kind.prepare_draw(integrand, &mut metadata)?.unwrap();
            assert_eq!(format!("{anchor:?}"), format!("{expected:?}"));
            if !fixed256 {
                use super::{
                    GaussianReferenceFunction, MomentumSpaceEvaluationInput,
                    StabilityEvaluationContext,
                };
                use crate::settings::runtime::IntegralUnit;

                // Exercise the real source/target scale owners with this already
                // generated one-incoming cut and its explicitly selected channel.
                let original_settings = integrand.settings.clone();
                let mut level = StabilityLevelSetting::default_arb();
                level.ecm_relative_tolerance_for_re = 1e-4;
                level.ecm_relative_tolerance_for_im = 1e-4;
                integrand.settings.stability.levels = vec![level];
                integrand.settings.stability.integrated_energy_dimension = None;
                let error = evaluate_from_source_precise(
                    integrand,
                    EvaluationTarget::Physical(model),
                    source_kind,
                    F(1.0),
                    false,
                    Complex::new_zero(),
                )
                .unwrap_err();
                assert!(error.to_string().contains("integrated_energy_dimension"));

                let default_sample = anchor.get_default_sample();
                assert_eq!(default_sample.external_moms().len(), 1);
                let loops = default_sample.loop_moms().0.len();
                let direct = MomentumSpaceEvaluationInput {
                    loop_momenta: default_sample
                        .loop_moms()
                        .iter()
                        .map(|momentum| {
                            ThreeMomentum::new(
                                momentum.px.into_ff64(),
                                momentum.py.into_ff64(),
                                momentum.pz.into_ff64(),
                            )
                        })
                        .collect(),
                    integrator_weight: F(1.0),
                    graph_id: None,
                    group_id: Some(GroupId(0)),
                    orientation: None,
                    channel_id: Some(SamplingChannelId(0)),
                };
                integrand.settings.sampling =
                    SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
                        sampling_type: DiscreteGraphSamplingType::SamplingMultiChanneling(
                            MultiChannelingSettings {
                                parameterization_settings: original_settings
                                    .sampling
                                    .get_parameterization_settings()
                                    .unwrap(),
                            },
                        ),
                        ..Default::default()
                    });
                integrand.warm_up(model)?;
                let raw_source = EvaluationSource::Momentum(&direct);
                let raw = raw_source.prepare_draw(integrand, &mut metadata)?.unwrap();
                let prepared_raw = EvaluationSource::Prepared {
                    sample: &raw,
                    original: &raw_source,
                };
                assert!(!prepared_raw.is_x_space());
                let reference_function =
                    GaussianReferenceFunction::centered(original_settings.kinematics.e_cm, loops)?;
                let maximum = Complex::new_zero();
                for energy in [
                    original_settings.kinematics.e_cm,
                    2.0 * original_settings.kinematics.e_cm,
                ] {
                    integrand.settings.kinematics.e_cm = energy;
                    let e_cm = F::<ArbPrec>::from_f64(energy);
                    for (evaluation_source, sample, missing) in [
                        (&source_kind, &anchor, 0),
                        (&prepared_raw, &raw, 3 * loops as i32),
                    ] {
                        let mut context = StabilityEvaluationContext {
                            target: EvaluationTarget::Reference(&reference_function),
                            source: evaluation_source,
                            stability_level: &level,
                            max_eval: &maximum,
                            wgt: F(1.0),
                            check_on_norm: false,
                            is_final_level: true,
                            evaluation_metadata: &mut metadata,
                            record_rotated_results: false,
                            precision_label: "ArbPrec",
                            escalate_if_exact_zero: false,
                            check_real: true,
                            check_imag: true,
                        };
                        integrand.settings.stability.integrated_energy_dimension = None;
                        assert_eq!(
                            context.ecm_reference_scale(integrand, sample)?.unwrap(),
                            e_cm.powi(-missing)
                        );
                        context.target = EvaluationTarget::Physical(model);
                        assert!(context.ecm_reference_scale(integrand, sample).is_err());
                        for (disabled_flux, dimension) in [(false, 1), (true, 2)] {
                            integrand.settings.general.disable_flux_factor = disabled_flux;
                            integrand.settings.stability.integrated_energy_dimension =
                                Some(dimension);
                            // The actual decay branch has no barn conversion,
                            // including when an explicit barn unit was supplied.
                            for unit in [
                                IntegralUnit::None,
                                IntegralUnit::Picobarn,
                                IntegralUnit::Femtobarn,
                            ] {
                                integrand.settings.general.integral_unit = unit;
                                assert_eq!(
                                    context.ecm_reference_scale(integrand, sample)?.unwrap(),
                                    e_cm.powi(dimension - missing)
                                );
                            }
                        }
                    }
                }
                integrand.settings = original_settings.clone();
                integrand.settings.stability.levels = vec![level];
                integrand.settings.stability.integrated_energy_dimension = None;
                integrand.warm_up(model)?;
                let reference_result = evaluate_from_source_precise(
                    integrand,
                    EvaluationTarget::Reference(&reference_function),
                    source_kind,
                    F(1.0),
                    false,
                    Complex::new_zero(),
                )?;
                assert!(match &reference_result {
                    super::PreciseEvaluationResult::Double(result) =>
                        result.reference_moments.is_some(),
                    super::PreciseEvaluationResult::Quad(result) =>
                        result.reference_moments.is_some(),
                    super::PreciseEvaluationResult::Arb(result) =>
                        result.reference_moments.is_some(),
                });
                reference_result.try_into_f64()?;
                integrand.settings = original_settings;
                integrand.warm_up(model)?;
            }
            let reference = integrand.get_graph(0).sampling_setup().clone();
            let probe = SamplingLawProbe {
                reference: &reference,
                channel_id: SamplingChannelId(0),
                radial_coordinate: 0,
                retry_lower_half: true,
                sign_coordinate: None,
                alternating_graph_sign: false,
                calls: std::array::from_fn(|_| AtomicUsize::new(0)),
            };
            let value = evaluate_from_source_precise(
                integrand,
                EvaluationTarget::SamplingLaw(&probe),
                EvaluationSource::XSpace(source),
                F(1.0),
                false,
                Complex::new_zero(),
            )?
            .try_into_f64()?;
            assert_eq!(
                value.evaluation_metadata.final_precision(),
                Some(Precision::Quad)
            );
            assert_eq!(probe.calls[0].load(Ordering::Relaxed), 1);
            assert_eq!(probe.calls[1].load(Ordering::Relaxed), 1);
            assert_eq!(probe.calls[2].load(Ordering::Relaxed), 0);
            assert!((value.integrand_result.re.0 - 1.0).abs() < 1.0e-8);
            let prepared = EvaluationSource::Prepared {
                sample: &anchor,
                original: &source_kind,
            };
            let arb_sample = prepared.build_gamma_sample::<ArbPrec, _>(integrand, &mut metadata)?;
            let arb_value = super::evaluate_single(
                integrand,
                EvaluationTarget::SamplingLaw(&probe),
                &arb_sample,
                &crate::momentum::Rotation::new(crate::momentum::RotationMethod::Identity),
                &mut metadata,
                Some(&anchor),
            )?;
            assert_eq!(probe.calls[2].load(Ordering::Relaxed), 1);
            assert!((arb_value.integrand_result.re.into_ff64().0 - 1.0).abs() < 1.0e-8);
            for calls in &probe.calls {
                calls.store(0, Ordering::Relaxed);
            }
            assert!(
                integrand
                    .get_graph(0)
                    .sampling_setup()
                    .sampling_bridge::<ArbPrec>()
                    .is_ok()
            );
            // A signed normalized density q(k) has integral zero and absolute
            // integral one when its reference radial CDF selects the sign.
            // Pair u with 1-u: two distinct channel points cancel in the signed
            // estimator, while their absolute contributions must still add.
            let mut absolute_runtime = integrand.clone();
            let catalogue = absolute_runtime
                .get_graph_mut(0)
                .sampling_setup_mut()
                .sampling_catalogue
                .as_mut()
                .unwrap();
            catalogue.entries.push(catalogue.entries[0].clone());
            let mut signed_probe = SamplingLawProbe {
                reference: &reference,
                channel_id: SamplingChannelId(0),
                radial_coordinate: 0,
                retry_lower_half: false,
                sign_coordinate: Some(0),
                alternating_graph_sign: false,
                calls: std::array::from_fn(|_| AtomicUsize::new(0)),
            };
            let (selection, original_coordinates) =
                super::gammaloop_sample::unwrap_sample::<f64>(source);
            for lower in [0.125, 0.25, 0.375] {
                let mut rows = Vec::new();
                for (channel, radial) in [lower, 1.0 - lower].into_iter().enumerate() {
                    let mut coordinates = original_coordinates.clone();
                    coordinates[0] = F(radial);
                    let source = symbolica::numerical_integration::Sample::Uniform(
                        F(1.0),
                        selection.clone(),
                        coordinates,
                    );
                    let mut draw = EvaluationSource::XSpace(&source)
                        .prepare_draw(&mut absolute_runtime, &mut metadata)?
                        .unwrap();
                    let mut row = draw.groups[0].1.remove(0);
                    row.sample.sample.jacobian /= row.sample.one().from_usize(2);
                    row.channel_id = Some(SamplingChannelId(channel));
                    row.prepared_lu_hosts.clear();
                    rows.push(row);
                }
                let draw = super::GammaLoopSample {
                    groups: vec![(Some(GroupId(0)), rows)],
                };
                let identity =
                    crate::momentum::Rotation::new(crate::momentum::RotationMethod::Identity);
                let value = super::evaluate_single(
                    &mut absolute_runtime,
                    EvaluationTarget::SamplingLaw(&signed_probe),
                    &draw,
                    &identity,
                    &mut metadata,
                    Some(&draw),
                )?;
                assert!(value.integrand_result.re.into_ff64().0.abs() < 1.0e-8);
                assert!(value.integrand_result.im.into_ff64().0.abs() < 1.0e-8);
                let absolute = value.absolute_integrand_result.unwrap();
                assert!((absolute.re.into_ff64().0 - 1.0).abs() < 1.0e-8);
                assert!((absolute.im.into_ff64().0 - 1.0).abs() < 1.0e-8);
                // Selected-channel MC has probability 1/2. Its individual
                // absolute weighted samples are also one, not two or one-half.
                for row in &draw.groups[0].1 {
                    let selected = super::GammaLoopSample {
                        groups: vec![(Some(GroupId(0)), vec![row.clone()])],
                    };
                    let value = super::evaluate_single(
                        &mut absolute_runtime,
                        EvaluationTarget::SamplingLaw(&signed_probe),
                        &selected,
                        &identity,
                        &mut metadata,
                        Some(&selected),
                    )?;
                    assert!(
                        (2.0 * value.absolute_integrand_result.unwrap().re.into_ff64().0 - 1.0)
                            .abs()
                            < 1.0e-8
                    );
                }
                // Opposite graph pieces at the SAME physical channel point
                // cancel before abs, including the imaginary component.
                let mut grouped = absolute_runtime.clone();
                grouped
                    .data
                    .graph_terms
                    .push(grouped.data.graph_terms[0].clone());
                grouped.data.graph_to_group_id.push(0);
                let mut grouped_draw = draw.clone();
                for row in &draw.groups[0].1 {
                    let mut partner = row.clone();
                    partner.graph_id = 1;
                    grouped_draw.groups[0].1.push(partner);
                }
                signed_probe.alternating_graph_sign = true;
                let canceled = super::evaluate_single(
                    &mut grouped,
                    EvaluationTarget::SamplingLaw(&signed_probe),
                    &grouped_draw,
                    &identity,
                    &mut metadata,
                    Some(&grouped_draw),
                )?;
                signed_probe.alternating_graph_sign = false;
                assert_eq!(canceled.integrand_result, Complex::new_re(draw.zero()));
                assert_eq!(
                    canceled.absolute_integrand_result,
                    Some(Complex::new_re(draw.zero()))
                );
            }
            let unavailable = super::sampling_maps::SamplingEvaluationError::Unrepresentable {
                operation: "test-only fixed source poison",
                detail: "cannot redraw at Arb".into(),
            };
            let setup = integrand.get_graph_mut(0).sampling_setup_mut();
            if fixed256 {
                setup.sampling_bridge_fixed256.set(Err(unavailable));
            } else {
                setup.sampling_bridge_quad.set(Err(unavailable));
            }
            let failure = evaluate_from_source_precise(
                integrand,
                EvaluationTarget::SamplingLaw(&probe),
                EvaluationSource::XSpace(source),
                F(1.0),
                false,
                Complex::new_zero(),
            )
            .unwrap_err();
            assert!(format!("{failure:#}").contains("fixed source precision"));
            assert!(
                probe
                    .calls
                    .iter()
                    .all(|calls| calls.load(Ordering::Relaxed) == 0)
            );
            let input = super::MomentumSpaceEvaluationInput {
                loop_momenta: vec![ThreeMomentum::new(F(0.1), F(0.2), F(0.3))],
                integrator_weight: F(1.0),
                graph_id: Some(0),
                group_id: None,
                orientation: None,
                channel_id: None,
            };
            let raw = EvaluationSource::Momentum(&input)
                .prepare_draw(integrand, &mut metadata)?
                .unwrap();
            assert_eq!(
                raw.get_default_sample().loop_moms().0[0].px,
                F(0.1).to_arb_exact()?
            );
            let selected_input = super::MomentumSpaceEvaluationInput {
                channel_id: Some(SamplingChannelId(0)),
                ..input
            };
            assert!(
                EvaluationSource::Momentum(&selected_input)
                    .prepare_draw(integrand, &mut metadata)
                    .is_err()
            );
        }
        runtime.get_mut_settings().sampling = original_sampling;
        // Fixed inputs and accuracy select an epoch, never an active point.
        let floor = F::<QuadFloat>::default().epsilon().sqrt().into_ff64().0;
        let settings = runtime.get_mut_settings();
        for level in &mut settings.stability.levels {
            level.required_precision_for_re = floor * 10.0;
            level.required_precision_for_im = floor * 10.0;
        }
        runtime.warm_up(model)?;
        let ProcessIntegrand::CrossSection(integrand) = runtime else {
            unreachable!()
        };
        let (precision, budget) = integrand.sampling_source_policy()?;
        assert!(budget >= floor);
        assert_eq!(precision, SamplingPrecision::Quad);
        runtime.get_mut_settings().stability.levels[2].required_precision_for_re = 1.0e-20;
        runtime.warm_up(model)?;
        let ProcessIntegrand::CrossSection(integrand) = runtime else {
            unreachable!()
        };
        assert_eq!(
            integrand
                .get_graph(0)
                .sampling_setup()
                .sampling_source
                .as_ref()
                .unwrap()
                .0,
            SamplingPrecision::Fixed256
        );
        let fixed256_floor = F::<SamplingFloat>::default().epsilon().sqrt().into_ff64().0;
        runtime.get_mut_settings().stability.levels[2].required_precision_for_re = fixed256_floor;
        runtime.warm_up(model)?;
        let ProcessIntegrand::CrossSection(integrand) = runtime else {
            unreachable!()
        };
        assert_eq!(
            integrand.sampling_source_policy()?.0,
            SamplingPrecision::Arb
        );
        assert!(
            integrand
                .get_graph(0)
                .sampling_setup()
                .sampling_bridge_fixed256
                .as_ref()
                .is_none()
        );
        Ok(())
    }

    /// Exercise source-private lifetime boundaries from the generated host fixture.
    pub(crate) fn check_host_source_transport<I: super::ProcessIntegrandImpl>(
        integrand: &mut I,
        model: &crate::model::Model,
        source: &symbolica::numerical_integration::Sample<F<f64>>,
        graph_id: usize,
        channel_id: SamplingChannelId,
    ) -> color_eyre::Result<()> {
        use super::{EvaluationSource, gammaloop_sample::unwrap_sample};
        use crate::{
            integrands::evaluation::EvaluationMetaData,
            momentum::{Rotatable, Rotation, RotationMethod},
            utils::{QuadFloat, newton_solver::RadialRootDiagnostics},
        };
        use symbolica::numerical_integration::Sample;

        assert_eq!(
            integrand
                .get_graph(0)
                .sampling_setup()
                .sampling_source
                .as_ref()
                .unwrap()
                .0,
            super::SamplingPrecision::Fixed256,
        );
        integrand.prepare_sampling_precision::<super::ArbPrec>()?;
        let mut metadata = EvaluationMetaData::new_empty();
        // Seed an existing physical occurrence so successful and failed source
        // preparation/debug operations must preserve nonempty history as well.
        metadata
            .radial_root_diagnostics
            .solve(
                &crate::utils::newton_solver::RadialRootIdentity::new(
                    "existing physical occurrence".into(),
                ),
                &F(0.0),
                &F(2.0),
                |r| (*r - F(1.0), F(1.0)),
                &F(1.0e-8),
                32,
                8,
                &F(1.0),
            )
            .unwrap();
        let original_source = EvaluationSource::XSpace(source);
        let mut anchor = original_source
            .prepare_draw(integrand, &mut metadata)?
            .unwrap();
        let sampling_preparation_time = metadata.canonical_sampling_preparation_time;
        assert!(sampling_preparation_time > std::time::Duration::ZERO);
        assert_eq!(metadata.parameterization_time, sampling_preparation_time);
        anchor.prepare_physical_overlaps(integrand, model, &mut metadata)?;
        // Clear only the test control's retained host so the physical owner must
        // solve its original Arb cut equation independently at this same point.
        let mut independent = anchor.clone();
        for (_, rows) in &mut independent.groups {
            for row in rows {
                row.prepared_lu_hosts.clear();
                row.physical_overlaps = None;
            }
        }
        let mut independent_metadata = EvaluationMetaData::new_empty();
        independent.prepare_physical_overlaps(integrand, model, &mut independent_metadata)?;
        let identity = Rotation::new(RotationMethod::Identity);
        let adopted_arb = super::evaluate_single(
            integrand,
            super::EvaluationTarget::Physical(model),
            &anchor,
            &identity,
            &mut EvaluationMetaData::new_empty(),
            Some(&anchor),
        )?;
        let independent_arb = super::evaluate_single(
            integrand,
            super::EvaluationTarget::Physical(model),
            &independent,
            &identity,
            &mut independent_metadata,
            Some(&independent),
        )?;
        for (adopted, solved) in [
            (
                &adopted_arb.integrand_result.re,
                &independent_arb.integrand_result.re,
            ),
            (
                &adopted_arb.integrand_result.im,
                &independent_arb.integrand_result.im,
            ),
        ] {
            let scale = adopted
                .abs()
                .max(solved.abs())
                .max(adopted.one() / adopted.from_usize(10).powi(25));
            assert!((adopted - solved).abs() <= scale / adopted.from_usize(10).powi(12));
        }
        assert_eq!(metadata.generated_event_count, 0);
        assert_eq!(metadata.accepted_event_count, 0);
        let preparation_time = metadata.canonical_physical_preparation_time;
        let prepared_source = EvaluationSource::Prepared {
            sample: &anchor,
            original: &original_source,
        };
        assert!(
            prepared_source
                .prepare_draw(integrand, &mut metadata)?
                .is_none()
        );
        let original = prepared_source.build_gamma_sample::<f64, _>(integrand, &mut metadata)?;
        assert_eq!(
            metadata.canonical_sampling_preparation_time,
            sampling_preparation_time
        );
        let prepared_lu_hosts = &original.groups[0].1[0].prepared_lu_hosts;
        assert_eq!(original.groups[0].1[0].channel_id, Some(channel_id));
        assert!(!prepared_lu_hosts.is_empty());
        for host in prepared_lu_hosts {
            assert_eq!(host.source.graph_id, graph_id);
            assert_eq!(host.source.generating_channel, channel_id);
            assert_eq!(host.source.target_channel, channel_id);
        }
        let history = format!("{:?}", metadata.radial_root_diagnostics);
        assert_ne!(history, format!("{:?}", RadialRootDiagnostics::default()));
        let rotation = Rotation::new(RotationMethod::Pi2Z);
        let rotated = original.rotate(&rotation, 17, 23);
        let canonical_overlaps = anchor.groups[0].1[0].physical_overlaps.as_ref().unwrap();
        assert!(std::sync::Arc::ptr_eq(
            canonical_overlaps,
            rotated.groups[0].1[0].physical_overlaps.as_ref().unwrap(),
        ));
        let retained = &rotated.groups[0].1[0].prepared_lu_hosts;
        assert_eq!(retained.len(), prepared_lu_hosts.len());
        for (before, after) in prepared_lu_hosts.iter().zip(retained) {
            assert_eq!(before.plan, after.plan);
            assert_eq!(before.source, after.source);
            assert_eq!(before.prior, after.prior);
            assert_eq!(before.solution.solution, after.solution.solution);
            assert_eq!(
                before.solution.derivative_at_solution,
                after.solution.derivative_at_solution
            );
            assert_eq!(
                before.solution.error_of_function,
                after.solution.error_of_function
            );
            assert_eq!(
                before.solution.num_iterations_used,
                after.solution.num_iterations_used
            );
            for t in [F(0.0), F(1.0), before.solution.solution] {
                assert_eq!(before.ray.evaluate(&t), after.ray.evaluate(&t));
            }
        }
        let original_momenta = original.get_default_sample().loop_moms();
        let rotated_momenta = rotated.get_default_sample().loop_moms();
        assert!(original_momenta.0.iter().any(|p| p != &p.rotate(&rotation)));
        for (before, after) in original_momenta.0.iter().zip(&rotated_momenta.0) {
            assert_eq!(*after, before.rotate(&rotation));
        }
        assert_eq!(rotated.get_default_sample().sample.loop_mom_cache_id, 17);

        // Successful physical adoption must also respect a nonidentity probe.
        // Keep the original source/host payload and compare the complete raised
        // estimator and normal identity/probe event policy, without another
        // generation or grid.
        let mut baseline_metadata = metadata.clone();
        let baseline = super::evaluate_single(
            integrand,
            super::EvaluationTarget::Physical(model),
            &original,
            &Rotation::new(RotationMethod::Identity),
            &mut baseline_metadata,
            Some(&anchor),
        )?;
        let mut rotated_metadata = metadata.clone();
        let physical_rotated = super::evaluate_single(
            integrand,
            super::EvaluationTarget::Physical(model),
            &rotated,
            &rotation,
            &mut rotated_metadata,
            Some(&anchor),
        )?;
        assert_eq!(
            baseline_metadata.canonical_physical_preparation_time,
            preparation_time
        );
        assert_eq!(
            rotated_metadata.canonical_physical_preparation_time,
            preparation_time
        );
        // Both identity and stability-probe evaluator work contributes to E;
        // preparing or adopting the retained source does not repeat C_S.
        for physical_metadata in [&baseline_metadata, &rotated_metadata] {
            assert!(
                physical_metadata.evaluator_evaluation_time > metadata.evaluator_evaluation_time
            );
            assert_eq!(
                physical_metadata.canonical_sampling_preparation_time,
                sampling_preparation_time
            );
        }
        for result in [&baseline, &physical_rotated] {
            assert!(result.integrand_result.re.0.is_finite());
            assert!(result.integrand_result.im.0.is_finite());
            assert!(result.integrand_result.re.0.abs() + result.integrand_result.im.0.abs() > 0.0);
        }
        for (a, b) in [
            (
                baseline.integrand_result.re.0,
                physical_rotated.integrand_result.re.0,
            ),
            (
                baseline.integrand_result.im.0,
                physical_rotated.integrand_result.im.0,
            ),
        ] {
            assert!(
                (a - b).abs() <= 1.0e-8 * a.abs().max(b.abs()).max(1.0e-25),
                "rotated physical contribution {b} != {a}"
            );
        }
        let baseline_events = baseline
            .event_groups
            .iter()
            .flat_map(|group| group.iter())
            .collect::<Vec<_>>();
        assert!(!baseline_events.is_empty());
        assert!(
            baseline_events
                .iter()
                .all(|event| event.weight.re.0.is_finite() && event.weight.im.0.is_finite())
        );
        // prepare_buffered_event intentionally retains/counts identity events
        // only; a rotated probe builds transient selector data when needed.
        // The surrounding fixture compares identity selected/direct event weights.
        assert!(physical_rotated.event_groups.is_empty());
        assert_eq!(physical_rotated.generated_event_count, 0);
        assert_eq!(physical_rotated.accepted_event_count, 0);

        // Losing a canonical accepted cut is not a request to choose a new
        // overlap. The real physical boundary must reject it before CT work.
        let mut incomplete = anchor.clone();
        let overlaps = std::sync::Arc::make_mut(
            incomplete.groups[0].1[0]
                .physical_overlaps
                .as_mut()
                .unwrap(),
        );
        let accepted = overlaps.iter_mut().find(|entry| entry.is_some()).unwrap();
        *accepted = None;
        let mut incomplete_metadata = metadata.clone();
        let error = super::evaluate_single(
            integrand,
            super::EvaluationTarget::Physical(model),
            &original,
            &Rotation::new(RotationMethod::Identity),
            &mut incomplete_metadata,
            Some(&incomplete),
        )
        .unwrap_err();
        assert!(
            error
                .downcast_ref::<super::sampling_maps::SamplingEvaluationError>()
                .is_some(),
            "{error:?}"
        );
        assert!(
            error.to_string().contains("selector/cut membership"),
            "{error:?}"
        );
        assert_eq!(
            incomplete_metadata.canonical_physical_preparation_time,
            preparation_time
        );

        // A claimed selected record with another graph ID is corruption,
        // not an ordinary graph-group member for which no authority is supplied.
        let mut corrupted = original.clone();
        let claimed_hosts = &mut corrupted.groups[0].1[0].prepared_lu_hosts;
        claimed_hosts[0].source.graph_id = graph_id.wrapping_add(1);
        let reference =
            super::GaussianReferenceFunction::new(1.0, vec![0.0; 3 * original_momenta.0.len()])?;
        let error = super::evaluate_single(
            integrand,
            super::EvaluationTarget::Reference(&reference),
            &corrupted,
            &rotation,
            &mut metadata,
            None,
        )
        .unwrap_err();
        assert!(
            error.to_string().contains("graph-group master"),
            "{error:?}"
        );
        assert_eq!(history, format!("{:?}", metadata.radial_root_diagnostics));

        // A valid alternate parent still does not belong to the selected
        // immutable channel. Exercise actual physical adoption, where this
        // mismatch and a forged active representative must fail before roots/events.
        for forge_parent in [true, false] {
            let mut corrupted = original.clone();
            let claimed_hosts = &mut corrupted.groups[0].1[0].prepared_lu_hosts;
            let plan = std::sync::Arc::make_mut(&mut claimed_hosts[0].plan);
            let expected = if forge_parent {
                assert!(plan.parent_lmb.len() >= 2);
                plan.parent_lmb.swap(0, 1);
                "selected channel definition"
            } else {
                plan.representative_cut_id.0 = plan.representative_cut_id.0.wrapping_add(1);
                "incompatible active representative"
            };
            let error = super::evaluate_single(
                integrand,
                super::EvaluationTarget::Physical(model),
                &corrupted,
                &rotation,
                &mut metadata,
                None,
            )
            .unwrap_err();
            assert!(error.to_string().contains(expected), "{error:?}");
            assert_eq!(history, format!("{:?}", metadata.radial_root_diagnostics));
        }

        // Norm/debug now consume this same retained draw and cannot prepare
        // roots or advance the physical retry occurrence/history owner.
        let before_time = metadata.parameterization_time;
        prepared_source.loop_norm_sum(integrand, &mut metadata)?;
        assert_eq!(history, format!("{:?}", metadata.radial_root_diagnostics));
        prepared_source.debug_sample(integrand, &mut metadata)?;
        assert_eq!(history, format!("{:?}", metadata.radial_root_diagnostics));
        let mut large_factor = anchor.clone();
        large_factor.groups[0].1[0].sample.sample.jacobian *=
            F::<crate::utils::ArbPrec>::from_f64(2.0).powi(2048);
        let large_factor_source = EvaluationSource::Prepared {
            sample: &large_factor,
            original: &original_source,
        };
        assert_eq!(
            large_factor_source.loop_norm_sum(integrand, &mut metadata)?,
            prepared_source.loop_norm_sum(integrand, &mut metadata)?
        );
        assert!(
            large_factor_source
                .debug_sample(integrand, &mut metadata)
                .is_err()
        );
        assert_eq!(history, format!("{:?}", metadata.radial_root_diagnostics));
        let (selection, mut coordinates) = unwrap_sample::<f64>(source);
        coordinates[0] = F(f64::NAN);
        let invalid = Sample::Uniform(F(1.0), selection, coordinates);
        assert!(
            EvaluationSource::XSpace(&invalid)
                .loop_norm_sum(integrand, &mut metadata)
                .is_err()
        );
        assert_eq!(history, format!("{:?}", metadata.radial_root_diagnostics));
        assert!(
            EvaluationSource::XSpace(&invalid)
                .debug_sample(integrand, &mut metadata)
                .is_err()
        );
        assert_eq!(history, format!("{:?}", metadata.radial_root_diagnostics));
        assert!(metadata.parameterization_time >= before_time);

        // Quad materializes the original canonical record directly; no earlier
        // Double geometry is promoted and no map/root is replayed.
        let quad = prepared_source.build_gamma_sample::<QuadFloat, _>(integrand, &mut metadata)?;
        let native = &quad.groups[0].1[0].prepared_lu_hosts;
        assert_eq!(native.len(), prepared_lu_hosts.len());
        for (before, after) in prepared_lu_hosts.iter().zip(native) {
            assert_eq!(before.plan, after.plan);
            assert_eq!(before.source, after.source);
            assert_eq!(before.prior.len(), after.prior.len());
        }
        let mut direct_momenta = original_momenta.0.clone();
        direct_momenta[0].px = F(0.1);
        let direct = super::MomentumSpaceEvaluationInput {
            loop_momenta: direct_momenta,
            integrator_weight: F(1.0),
            graph_id: Some(graph_id),
            group_id: None,
            channel_id: None,
            orientation: None,
        };
        let raw_source = EvaluationSource::Momentum(&direct);
        let mut raw_anchor = raw_source.prepare_draw(integrand, &mut metadata)?.unwrap();
        raw_anchor.prepare_physical_overlaps(integrand, model, &mut metadata)?;
        assert!(
            raw_anchor
                .groups
                .iter()
                .flat_map(|(_, rows)| rows)
                .all(|row| {
                    row.channel_id.is_none()
                        && row.prepared_lu_hosts.is_empty()
                        && row.physical_overlaps.is_some()
                })
        );
        let raw_prepared = EvaluationSource::Prepared {
            sample: &raw_anchor,
            original: &raw_source,
        };
        let raw_quad = raw_prepared.build_gamma_sample::<QuadFloat, _>(integrand, &mut metadata)?;
        assert_eq!(
            raw_quad.get_default_sample().loop_moms().0[0].px.0,
            <QuadFloat as crate::utils::FloatLike>::from_f64_exact_binary(0.1)
        );
        assert!(
            raw_quad
                .groups
                .iter()
                .flat_map(|(_, rows)| rows)
                .all(|row| row.channel_id.is_none() && row.prepared_lu_hosts.is_empty())
        );
        Ok(())
    }

    /// Reuse the generated kite fixture while retaining the real source,
    /// partition and native stability owners. The final target alone is a probe.
    pub(super) fn check_joint_policy_source_replay(
        integrand: &super::ProcessIntegrand,
        model: &crate::model::Model,
        joint_id: super::SamplingChannelId,
        cube: &[f64],
    ) -> color_eyre::Result<()> {
        use super::{
            EvaluationSource, EvaluationTarget, GraphTerm, ProcessIntegrand, ProcessIntegrandImpl,
            SamplingLawProbe, evaluate_from_source_precise,
        };
        use crate::{
            integrands::evaluation::StabilityStatus,
            settings::runtime::{
                Precision, SamplingChannelWeight, SamplingSettingsParser, StabilityLevelSetting,
                SumMode,
            },
            utils::ArbPrec,
        };
        use spenso::algebra::complex::Complex;
        use std::sync::atomic::{AtomicUsize, Ordering};
        use symbolica::numerical_integration::Sample;

        let mut runtime = integrand.clone();
        let mut parser: SamplingSettingsParser =
            toml::from_str(&toml::to_string(&runtime.get_settings().sampling)?)?;
        parser.graphs = SumMode::MonteCarlo;
        parser.orientations = SumMode::Summed;
        parser.sampling_channels = SumMode::MonteCarlo;
        parser.sampling_channel_weight = SamplingChannelWeight::SingularityProxy;
        for definitions in parser.channel_definitions.values_mut() {
            for definition in definitions.values_mut() {
                definition.singularity_proxy = Some("1".into());
            }
        }
        runtime.get_mut_settings().sampling = toml::from_str(&toml::to_string(&parser)?)?;
        let stability = &mut runtime.get_mut_settings().stability;
        stability.rotation_axis.clear();
        stability.escalate_if_exact_zero = false;
        stability.loop_momenta_norm_escalation_factor = 0.0;
        stability.levels = vec![
            StabilityLevelSetting::default_double(),
            StabilityLevelSetting::default_quad(),
        ];
        runtime.warm_up(model)?;
        let ProcessIntegrand::Amplitude(amplitude) = &mut runtime else {
            unreachable!("the fixture is the generated kite amplitude")
        };
        amplitude.prepare_sampling_precision::<ArbPrec>()?;
        let reference = amplitude.get_graph(0).sampling_setup().clone();
        assert_eq!(reference.sampling_bridge::<f64>()?.channels().len(), 2);
        let ordinary_name = reference
            .sampling_bridge::<f64>()?
            .channels()
            .iter()
            .find(|channel| channel.contract().support == super::SamplingSupport::Full)
            .unwrap()
            .name
            .clone();
        let mut probe = SamplingLawProbe {
            reference: &reference,
            channel_id: joint_id,
            radial_coordinate: 3,
            retry_lower_half: true,
            sign_coordinate: None,
            alternating_graph_sign: false,
            calls: std::array::from_fn(|_| AtomicUsize::new(0)),
        };
        for graph in amplitude.get_terms_mut() {
            let unavailable = super::sampling_maps::SamplingEvaluationError::Unrepresentable {
                operation: "test-only native map poison",
                detail: "physical retry must consume the canonical draw".into(),
            };
            graph
                .sampling_setup_mut()
                .sampling_bridge
                .set(Err(unavailable.clone()));
            graph
                .sampling_setup_mut()
                .sampling_bridge_quad
                .set(Err(unavailable));
        }
        for radial in [0.25, 0.75] {
            for calls in &probe.calls {
                calls.store(0, Ordering::Relaxed);
            }
            let mut coordinates = cube.to_vec();
            coordinates[3] = radial;
            // The outer weight includes the inverse channel probability 2.
            // Both support-gated constant proxies equal one on q_rho's patch.
            let source = Sample::Discrete(
                F(2.0),
                0,
                Some(Box::new(Sample::Discrete(
                    F(1.0),
                    joint_id.index(),
                    Some(Box::new(Sample::Continuous(
                        F(1.0),
                        coordinates.into_iter().map(F).collect(),
                    ))),
                ))),
            );
            let ProcessIntegrand::Amplitude(amplitude) = &mut runtime else {
                unreachable!()
            };
            let result = evaluate_from_source_precise(
                amplitude,
                EvaluationTarget::SamplingLaw(&probe),
                EvaluationSource::XSpace(&source),
                F(2.0),
                false,
                Complex::new_zero(),
            )?
            .try_into_f64()?;
            let retried = radial < 0.5;
            assert!(
                !result.evaluation_metadata.is_nan
                    && result.integrand_result.re.0.is_finite()
                    && result.integrand_result.im.0.is_finite(),
                "nonfinite final retry result at radial {radial}: {result:?}"
            );
            assert_eq!(
                probe.calls[0].load(Ordering::Relaxed),
                1,
                "Double body call count at radial {radial}: {:?}",
                result.evaluation_metadata.stability_results
            );
            assert_eq!(
                probe.calls[1].load(Ordering::Relaxed),
                usize::from(retried),
                "Quad body call count at radial {radial}: {:?}",
                result.evaluation_metadata.stability_results
            );
            assert_eq!(
                probe.calls[2].load(Ordering::Relaxed),
                0,
                "canonical preparation must not evaluate the target"
            );
            assert_eq!(
                result.evaluation_metadata.final_precision(),
                Some(if retried {
                    Precision::Quad
                } else {
                    Precision::Double
                })
            );
            assert_eq!(
                result.evaluation_metadata.stability_results.len(),
                1 + usize::from(retried)
            );
            if retried {
                let first = &result.evaluation_metadata.stability_results[0];
                assert_eq!(first.precision, Precision::Double);
                // A single identity probe is reported as Unknown regardless
                // of its acceptance. Body calls and the two native levels
                // above establish the intended nonfinite-Double retry.
                assert_eq!(first.status, StabilityStatus::Unknown);
            }
            assert!(
                (result.integrand_result.re.0 * result.integrator_weight.0 - 1.0).abs() < 1.0e-8
            );
            let policies = result
                .evaluation_metadata
                .sampling_proposal_policies
                .clone();
            assert!(!policies.is_collecting());
            assert_eq!(policies.len(), 1);
            let mut replay_metadata = result.evaluation_metadata.clone();
            let before = replay_metadata.parameterization_time;
            EvaluationSource::XSpace(&source).debug_sample(amplitude, &mut replay_metadata)?;
            assert!(replay_metadata.parameterization_time >= before);
            assert_eq!(
                replay_metadata.integrand_evaluation_time,
                result.evaluation_metadata.integrand_evaluation_time
            );
            assert_eq!(replay_metadata.sampling_proposal_policies, policies);

            let mut quad_only = runtime.clone();
            quad_only.get_mut_settings().stability.levels =
                vec![StabilityLevelSetting::default_quad()];
            quad_only.warm_up(model)?;
            let ProcessIntegrand::Amplitude(amplitude) = &mut quad_only else {
                unreachable!()
            };
            let direct = evaluate_from_source_precise(
                amplitude,
                EvaluationTarget::SamplingLaw(&probe),
                EvaluationSource::XSpace(&source),
                F(2.0),
                false,
                Complex::new_zero(),
            )?
            .try_into_f64()?;
            assert!((result.integrand_result.re.0 - direct.integrand_result.re.0).abs() < 1.0e-10);
            assert_eq!(
                result.evaluation_metadata.sampling_proposal_policies,
                direct.evaluation_metadata.sampling_proposal_policies
            );

            // An active-point norm may choose Quad immediately, but cannot
            // bypass the canonical phase or initialize policy in that lane.
            let mut norm_escalated = runtime.clone();
            norm_escalated
                .get_mut_settings()
                .stability
                .loop_momenta_norm_escalation_factor = 1.0e-6;
            norm_escalated.warm_up(model)?;
            let ProcessIntegrand::Amplitude(amplitude) = &mut norm_escalated else {
                unreachable!()
            };
            let escalated = evaluate_from_source_precise(
                amplitude,
                EvaluationTarget::SamplingLaw(&probe),
                EvaluationSource::XSpace(&source),
                F(2.0),
                false,
                Complex::new_zero(),
            )?
            .try_into_f64()?;
            assert_eq!(
                escalated.evaluation_metadata.final_precision(),
                Some(Precision::Quad)
            );
            assert_eq!(escalated.evaluation_metadata.stability_results.len(), 1);
            assert_eq!(
                escalated.evaluation_metadata.sampling_proposal_policies,
                policies
            );
            assert_eq!(escalated.integrand_result, direct.integrand_result);
        }

        // Ordinary channels also retain one fixed canonical proposal, while
        // requiring no joint decisions even if unused joint definitions remain.
        parser.default_channel_selection = vec![ordinary_name];
        parser.channel_selection.clear();
        runtime.get_mut_settings().sampling = toml::from_str(&toml::to_string(&parser)?)?;
        runtime.warm_up(model)?;
        let ProcessIntegrand::Amplitude(amplitude) = &mut runtime else {
            unreachable!()
        };
        assert!(
            amplitude
                .get_graph(0)
                .sampling_setup()
                .sampling_bridge::<ArbPrec>()
                .is_ok()
        );
        probe.retry_lower_half = false;
        let source = Sample::Uniform(F(1.0), vec![0, 0], cube.iter().copied().map(F).collect());
        let result = evaluate_from_source_precise(
            amplitude,
            EvaluationTarget::SamplingLaw(&probe),
            EvaluationSource::XSpace(&source),
            F(1.0),
            false,
            Complex::new_zero(),
        )?
        .try_into_f64()?;
        assert!(
            result
                .evaluation_metadata
                .sampling_proposal_policies
                .is_empty()
        );
        assert!(
            amplitude
                .get_graph(0)
                .sampling_setup()
                .sampling_bridge::<ArbPrec>()
                .is_ok()
        );
        // A failed fixed proposal cannot choose a native map as a fallback,
        // regardless of which physical precision would otherwise be usable.
        for calls in &probe.calls {
            calls.store(0, Ordering::Relaxed);
        }
        amplitude
            .get_graph_mut(0)
            .sampling_setup_mut()
            .sampling_bridge_fixed256
            .set(Err(
                super::sampling_maps::SamplingEvaluationError::Unrepresentable {
                    operation: "test-only canonical map poison",
                    detail: "fixed proposal unavailable".into(),
                },
            ));
        let failure = evaluate_from_source_precise(
            amplitude,
            EvaluationTarget::SamplingLaw(&probe),
            EvaluationSource::XSpace(&source),
            F(1.0),
            false,
            Complex::new_zero(),
        )
        .unwrap_err();
        assert!(
            format!("{failure:#}").contains("fixed source precision"),
            "{failure:?}"
        );
        assert!(
            probe
                .calls
                .iter()
                .all(|calls| calls.load(Ordering::Relaxed) == 0)
        );
        Ok(())
    }

    #[test]
    fn active_dependent_retry_changes_a_normalized_radial_law() {
        // Even identical unit determinants and inverse densities do not fix a
        // law when an active-dependent retry changes a continuous chart. Circle
        // maps x and x+delta are each uniform, but switching only x<delta leaves
        // [0,delta) empty. A normalized narrow target has mean zero rather than
        // one although the same-point J*q identity holds exactly in both maps.
        for denominator in [32, 1024, 16384] {
            let delta = 1.0 / f64::from(denominator);
            let original = 0.5 * delta;
            let switched = original + delta;
            assert!(original < delta && switched >= delta);
            assert_eq!(delta * (1.0 / delta), 1.0);
            assert_eq!(delta * 0.0, 0.0);
        }

        // For f=q_rho, uniform R in (0,rho) has unit weight. If only the
        // lower-half draws replay with 2rho, they stay inside f's support but
        // receive weight two. This oracle is independent of map code.
        let retained: f64 = [0.25, 0.75].into_iter().map(|_| 1.0).sum::<f64>() / 2.0;
        let switched = [0.25, 0.75]
            .into_iter()
            .map(|u| {
                let radius = if u < 0.5 { 2.0 } else { 1.0 };
                let reference_density = if radius * u < 1.0 { 1.0 } else { 0.0 };
                radius * reference_density
            })
            .sum::<f64>()
            / 2.0;
        assert_eq!(retained, 1.0);
        assert_eq!(switched, 1.5);
        // Equal channel probabilities and a full ordinary sibling with its
        // unchanged unit conditional mean give a complete mean of 5/4.
        assert_eq!((switched + 1.0) / 2.0, 1.25);
    }

    #[test]
    fn stability_checks_reject_nonfinite_probes_at_every_level() {
        use super::{StabilityFailureReason, StabilityLevelSetting};
        use spenso::algebra::complex::Complex;

        let level = StabilityLevelSetting::default_double();
        for check_on_norm in [false, true] {
            let check = if check_on_norm {
                super::stability_check_on_norm::<f64>
            } else {
                super::stability_check::<f64>
            };
            for is_final_level in [false, true] {
                for nonfinite in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                    for invalid in [(nonfinite, 2.0), (2.0, nonfinite)] {
                        for probes in [
                            vec![invalid, (2.0, 3.0)],
                            vec![(2.0, 3.0), invalid],
                            vec![invalid],
                        ] {
                            let results = probes
                                .into_iter()
                                .map(|(re, im)| Complex::new(F(re), F(im)))
                                .collect::<Vec<_>>();
                            let (result, accuracy, stable, reason) = check(
                                None,
                                &results,
                                &level,
                                Complex::new_zero(),
                                F(1.0),
                                is_final_level,
                                false,
                            );
                            assert!(
                                !stable,
                                "nonfinite probe accepted: norm={check_on_norm}, final={is_final_level}, probes={results:?}"
                            );
                            assert!(accuracy.is_none());
                            assert_eq!(reason, Some(StabilityFailureReason::ErrorThreshold));
                            // Preserve even a nonfinite primary for the existing
                            // validity flag and later component sanitization.
                            assert_eq!(result.re.0.to_bits(), results[0].re.0.to_bits());
                            assert_eq!(result.im.0.to_bits(), results[0].im.0.to_bits());
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn real_only_stability_ignores_inactive_imaginary_probe() {
        use super::StabilityFailureReason;
        use crate::settings::runtime::StabilityLevelSetting;
        use spenso::algebra::complex::Complex;

        let level = StabilityLevelSetting::default_double();
        let results = vec![Complex::new(F(2.0), F(1.0)), Complex::new(F(2.0), F(2.0))];
        let (_, _, real_stable, _) = super::stability_check_components(
            None,
            &results,
            &level,
            Complex::new(F(10.0), F(10.0)),
            F(1.0),
            false,
            false,
            true,
            false,
        );
        let (_, _, both_stable, reason) = super::stability_check_components(
            None,
            &results,
            &level,
            Complex::new(F(10.0), F(10.0)),
            F(1.0),
            false,
            false,
            true,
            true,
        );
        assert!(real_stable);
        assert!(!both_stable);
        assert_eq!(reason, Some(StabilityFailureReason::ErrorThreshold));

        let nonfinite_im = vec![Complex::new(F(2.0), F(f64::NAN))];
        let (_, _, stable, _) = super::stability_check_components(
            None,
            &nonfinite_im,
            &level,
            Complex::new(F(10.0), F(10.0)),
            F(1.0),
            true,
            false,
            true,
            false,
        );
        assert!(stable);
    }

    #[test]
    fn stability_checks_preserve_finite_controls_at_every_level() {
        use super::{StabilityFailureReason, StabilityLevelSetting, StabilityStatus};
        use spenso::algebra::complex::Complex;

        let level = StabilityLevelSetting::default_double();
        for check_on_norm in [false, true] {
            let check = if check_on_norm {
                super::stability_check_on_norm::<f64>
            } else {
                super::stability_check::<f64>
            };
            for is_final_level in [false, true] {
                for count in [1, 2] {
                    for (re, im) in [(0.0, 0.0), (2.0, 3.0)] {
                        let results = vec![Complex::new(F(re), F(im)); count];
                        let (result, accuracy, stable, reason) = check(
                            None,
                            &results,
                            &level,
                            Complex::new_zero(),
                            F(1.0),
                            is_final_level,
                            false,
                        );
                        assert_eq!(result, results[0]);
                        assert!(stable);
                        assert_eq!(reason, None);
                        if count == 1 {
                            assert!(accuracy.is_none());
                            assert_eq!(
                                StabilityStatus::from_sample_count(count, stable),
                                StabilityStatus::Unknown
                            );
                        } else {
                            assert_eq!(accuracy, Some(F(0.0)));
                        }
                    }
                }

                let results = [Complex::new(F(2.0), F(3.0)), Complex::new(F(4.0), F(6.0))];
                let (result, accuracy, stable, reason) = check(
                    None,
                    &results,
                    &level,
                    Complex::new_zero(),
                    F(1.0),
                    is_final_level,
                    false,
                );
                assert!(!stable);
                assert!(accuracy.unwrap() > F(level.required_precision_for_re));
                assert_eq!(reason, Some(StabilityFailureReason::ErrorThreshold));
                assert_eq!(
                    result,
                    if check_on_norm {
                        results[0]
                    } else {
                        Complex::new(F(3.0), F(4.5))
                    }
                );
            }
        }
    }

    #[test]
    fn stability_checks_bound_complete_underflow_without_cancellation() {
        use crate::{settings::runtime::StabilityLevelSetting, utils::ArbPrec};
        use spenso::algebra::complex::Complex;

        let level = StabilityLevelSetting::default_arb();
        let one = F::<ArbPrec>::default().one();
        let minimum = F::<ArbPrec>::from_f64(f64::MIN_POSITIVE);
        // One probe is normal-sized, but the cancellation-free bound on the
        // returned average is 7/8 of MIN_POSITIVE, as in GL638 sample 3559.
        let probes = [
            Complex::new_re(&minimum / one.from_usize(2)),
            Complex::new_re(-&minimum * one.from_usize(5) / one.from_usize(4)),
        ];
        let average = (&probes[0] + &probes[1]) / one.from_usize(2);
        for norm in [false, true] {
            let check = if norm {
                super::stability_check_on_norm::<ArbPrec>
            } else {
                super::stability_check::<ArbPrec>
            };
            for (weight, expected) in [
                (one.clone(), true),
                (one.from_usize(2), false),
                (one.from_usize(10).powi(100), false),
            ] {
                let (result, accuracy, stable, _) = check(
                    None,
                    &probes,
                    &level,
                    Complex::new_re(one.zero()),
                    weight,
                    true,
                    false,
                );
                assert_eq!(stable, expected);
                assert!(accuracy.unwrap() > F::from_f64(level.required_precision_for_re));
                assert_eq!(
                    result,
                    if norm {
                        probes[0].clone()
                    } else {
                        average.clone()
                    }
                );
            }
            // Norm checking returns the primary rather than the average, so a
            // normal-sized primary cannot inherit the average's waiver.
            assert_eq!(
                check(
                    None,
                    &[probes[1].clone(), probes[0].clone()],
                    &level,
                    Complex::new_re(one.zero()),
                    one.clone(),
                    true,
                    false,
                )
                .2,
                !norm,
            );
            for invalid_weight in [f64::NAN, f64::INFINITY] {
                assert!(
                    !check(
                        None,
                        &probes,
                        &level,
                        Complex::new_re(one.zero()),
                        F::from_f64(invalid_weight),
                        true,
                        false,
                    )
                    .2
                );
            }
        }
        // A small signed average does not bound large cancelling probes.
        assert!(
            !super::stability_check(
                None,
                &[Complex::new_re(one.clone()), Complex::new_re(-one.clone())],
                &level,
                Complex::new_re(one.zero()),
                one,
                true,
                false,
            )
            .2
        );
    }

    #[test]
    fn stability_underflow_waivers_are_componentwise() {
        use crate::{settings::runtime::StabilityLevelSetting, utils::ArbPrec};
        use spenso::algebra::complex::Complex;

        let level = StabilityLevelSetting::default_arb();
        let one = F::<ArbPrec>::default().one();
        let tiny = F::<ArbPrec>::from_f64(f64::MIN_POSITIVE) / one.from_usize(4);
        for swap in [false, true] {
            for (second_normal, expected) in [(1, true), (2, false)] {
                let mut probes = [
                    Complex::new(tiny.clone(), one.clone()),
                    Complex::new(-tiny.clone(), one.from_usize(second_normal)),
                ];
                if swap {
                    for probe in &mut probes {
                        std::mem::swap(&mut probe.re, &mut probe.im);
                    }
                }
                assert_eq!(
                    super::stability_check(
                        None,
                        &probes,
                        &level,
                        Complex::new_re(one.zero()),
                        one.clone(),
                        true,
                        false,
                    )
                    .2,
                    expected,
                );
            }
        }
    }

    #[test]
    fn stability_ecm_tolerance_preserves_energy_and_unit_scaling() {
        use crate::settings::runtime::{IntegralUnit, StabilityLevelSetting};
        use crate::utils::ArbPrec;
        use spenso::algebra::complex::Complex;

        let one = F::<ArbPrec>::default().one();
        let mut level = StabilityLevelSetting::default_arb();
        level.ecm_relative_tolerance_for_re = 1e-5;
        for energy in [600, 1200] {
            // Cross sections, widths and a scalar two-loop amplitude have
            // different dimensions; the decision concerns a dimensionless ratio.
            for dimension in [-6, -2, 0, 1, 2] {
                for missing_measure in [0, 6] {
                    for unit in [
                        IntegralUnit::None,
                        IntegralUnit::Picobarn,
                        IntegralUnit::Femtobarn,
                    ] {
                        let scale = one.from_usize(energy).powi(dimension - missing_measure)
                            * super::cross_section::barn_conversion_factor(unit, one.clone());
                        let small = &scale / one.from_usize(1_000_000);
                        let probes = [
                            Complex::new_re(small.clone()),
                            Complex::new_re(&small * one.from_usize(3)),
                        ];
                        for norm in [false, true] {
                            let check = if norm {
                                super::stability_check_on_norm::<ArbPrec>
                            } else {
                                super::stability_check::<ArbPrec>
                            };
                            for (weight, expected) in [(1, true), (-1, true), (100, false)] {
                                let (result, accuracy, stable, _) = check(
                                    Some(&scale),
                                    &probes,
                                    &level,
                                    Complex::new_re(one.zero()),
                                    one.from_i64(weight),
                                    true,
                                    false,
                                );
                                assert_eq!(stable, expected);
                                assert!(
                                    accuracy.unwrap()
                                        > F::from_f64(level.required_precision_for_re)
                                );
                                assert_eq!(
                                    result.re,
                                    if norm {
                                        small.clone()
                                    } else {
                                        &small * one.from_usize(2)
                                    }
                                );
                            }
                        }
                    }
                }
            }
        }
        // A body below binary64 range must be multiplied in native precision:
        // its huge outer factor can put the disagreement well above the allowance.
        let tiny = one.from_usize(10).powi(-1000);
        let probes = [
            Complex::new_re(tiny.clone()),
            Complex::new_re(&tiny * one.from_usize(3)),
        ];
        level.ecm_relative_tolerance_for_re = 1e-201;
        assert!(
            !super::stability_check(
                Some(&one),
                &probes,
                &level,
                Complex::new_re(one.zero()),
                one.from_usize(10).powi(800),
                true,
                false
            )
            .2
        );
    }

    #[test]
    fn stability_ecm_tolerance_bounds_the_returned_components() {
        use crate::settings::runtime::StabilityLevelSetting;
        use spenso::algebra::complex::Complex;

        let scale = F(1.0);
        let mut level = StabilityLevelSetting::default_double();
        let probes = [Complex::new(F(1.0), F(0.0)), Complex::new(F(3.0), F(0.0))];
        for (allowance, expected_component, expected_norm) in [
            (0.0, false, false),
            (0.5, false, false),
            (1.0, true, false),
            (2.0, true, true),
        ] {
            level.ecm_relative_tolerance_for_re = allowance;
            for (norm, expected, returned) in
                [(false, expected_component, 2.0), (true, expected_norm, 1.0)]
            {
                let check = if norm {
                    super::stability_check_on_norm::<f64>
                } else {
                    super::stability_check::<f64>
                };
                let (value, accuracy, stable, _) = check(
                    Some(&scale),
                    &probes,
                    &level,
                    Complex::new_zero(),
                    F(1.0),
                    true,
                    false,
                );
                assert_eq!(stable, expected);
                assert_eq!(value.re, F(returned));
                assert!(accuracy.unwrap() > F(level.required_precision_for_re));
            }
        }
        // A real allowance cannot excuse an independently unstable imaginary
        // component; swapping Re and Im swaps the applicable setting too.
        for swap in [false, true] {
            let mut mixed = [Complex::new(F(1.0), F(4.0)), Complex::new(F(3.0), F(8.0))];
            level.ecm_relative_tolerance_for_re = 1.0;
            level.ecm_relative_tolerance_for_im = 0.0;
            if swap {
                for probe in &mut mixed {
                    std::mem::swap(&mut probe.re, &mut probe.im);
                }
                std::mem::swap(
                    &mut level.ecm_relative_tolerance_for_re,
                    &mut level.ecm_relative_tolerance_for_im,
                );
            }
            assert!(
                !super::stability_check(
                    Some(&scale),
                    &mixed,
                    &level,
                    Complex::new_zero(),
                    F(1.0),
                    true,
                    false
                )
                .2
            );
        }
        level.ecm_relative_tolerance_for_re = 0.5;
        let cancelling = [Complex::new_re(F(1.0)), Complex::new_re(F(-1.0))];
        assert!(
            !super::stability_check(
                Some(&scale),
                &cancelling,
                &level,
                Complex::new_zero(),
                F(1.0),
                true,
                false
            )
            .2
        );
        for check in [
            super::stability_check::<f64>,
            super::stability_check_on_norm::<f64>,
        ] {
            level.ecm_relative_tolerance_for_re = 10.0;
            for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                assert!(
                    !check(
                        Some(&scale),
                        &probes,
                        &level,
                        Complex::new_zero(),
                        F(invalid),
                        true,
                        false
                    )
                    .2
                );
                assert!(
                    !check(
                        Some(&F(invalid)),
                        &probes,
                        &level,
                        Complex::new_zero(),
                        F(1.0),
                        true,
                        false
                    )
                    .2
                );
            }
        }
        // Finite differences do not make overflowing complete products safe.
        level.required_precision_for_re = 0.0;
        let large = [
            Complex::new_re(F(1e200)),
            Complex::new_re(F(1e200 * (1.0 - 1e-14))),
        ];
        assert!(
            !super::stability_check(
                Some(&F(1e308)),
                &large,
                &level,
                Complex::new_zero(),
                F(1e110),
                true,
                false
            )
            .2
        );
    }

    #[test]
    fn stability_ecm_runtime_validation_precedes_orientation_checks() {
        let mut settings = RuntimeSettings::default();
        for component in ["re", "im"] {
            for invalid in [-1.0, f64::NAN, f64::INFINITY] {
                let level = &mut settings.stability.levels[0];
                level.ecm_relative_tolerance_for_re = 0.0;
                level.ecm_relative_tolerance_for_im = 0.0;
                if component == "re" {
                    level.ecm_relative_tolerance_for_re = invalid;
                } else {
                    level.ecm_relative_tolerance_for_im = invalid;
                }
                for explicit in [false, true] {
                    let error = validate_process_runtime_settings(&settings, explicit).unwrap_err();
                    assert!(
                        error
                            .to_string()
                            .contains(&format!("levels[0].ecm_relative_tolerance_for_{component}"))
                    );
                }
            }
        }
        settings.stability.levels[0].ecm_relative_tolerance_for_re = 1e-100;
        settings.stability.levels[0].ecm_relative_tolerance_for_im = 0.0;
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            settings.kinematics.e_cm = invalid;
            assert!(
                validate_process_runtime_settings(&settings, false)
                    .unwrap_err()
                    .to_string()
                    .contains("kinematics.e_cm")
            );
        }
        settings.kinematics.e_cm = 600.0;
        // The physical dimension is admitted at its target boundary; generic
        // warmup also serves reference functions with their known dimensions.
        assert!(settings.stability.integrated_energy_dimension.is_none());
        validate_process_runtime_settings(&settings, false).unwrap();
    }

    #[test]
    fn precise_event_normalization_preserves_prior_factors_and_partial_weights() {
        use crate::{
            observables::{
                AdditionalWeightKey, GenericEvent, GenericEventGroup, GenericEventGroupList,
            },
            utils::ArbPrec,
        };
        use spenso::algebra::complex::Complex;

        let one = F::<ArbPrec>::default().one();
        let original = Complex::new_re(one.from_usize(3));
        let counterterm = Complex::new_re(-one.clone());
        let prior_factor = Complex::new_re(one.from_usize(5));
        let mut event = GenericEvent::<ArbPrec> {
            weight: (&original + &counterterm) * &prior_factor,
            ..Default::default()
        };
        event.additional_weights.weights.extend([
            (AdditionalWeightKey::Original, original.clone()),
            (
                AdditionalWeightKey::ThresholdCounterterm { subset_index: 0 },
                counterterm.clone(),
            ),
            (
                AdditionalWeightKey::FullMultiplicativeFactor,
                prior_factor.clone(),
            ),
        ]);
        let mut events = GenericEventGroupList(vec![GenericEventGroup(vec![event])]);
        let final_factor = super::full_event_multiplicative_factor_precise(
            Some(one.from_usize(7)),
            one.from_usize(11),
        );
        super::apply_full_event_multiplicative_factor_precise(&mut events, &final_factor);
        let event = &events[0][0];
        let weights = &event.additional_weights.weights;
        assert_eq!(weights[&AdditionalWeightKey::Original], original);
        assert_eq!(
            weights[&AdditionalWeightKey::ThresholdCounterterm { subset_index: 0 }],
            counterterm
        );
        assert_eq!(
            weights[&AdditionalWeightKey::FullMultiplicativeFactor],
            &prior_factor * &final_factor
        );
        assert_eq!(
            event.weight,
            (&original + &counterterm) * &weights[&AdditionalWeightKey::FullMultiplicativeFactor]
        );
        assert_eq!(event.weight, Complex::new_re(one.from_usize(770)));
    }

    #[test]
    fn native_sampling_product_and_precise_reporting_preserve_full_range() {
        use crate::{
            integrands::evaluation::{
                EvaluationMetaData, GraphEvaluationResult, PreciseEvaluationResult,
            },
            integrands::process::sampling_reference::ReferenceMoments,
            observables::{GenericEvent, GenericEventGroup, GenericEventGroupList},
            settings::runtime::Precision,
            utils::ArbPrec,
        };
        use spenso::algebra::complex::Complex;
        use std::time::Duration;

        let one = F::<ArbPrec>::default().one();
        let large = one.from_usize(10).powi(400);
        let small = large.clone().inv();
        assert!(large.clone().into_ff64().0.is_infinite());
        assert_eq!(small.clone().into_ff64().0, 0.0);
        let mut graph = GraphEvaluationResult::zero(one.zero());
        graph.integrand_result = Complex::new_re(small.clone());
        graph.reference_moments = Some(ReferenceMoments {
            second_moment: small.clone(),
            jacobian_min: 1.0,
            jacobian_max: 1.0,
        });
        graph.event_groups = GenericEventGroupList(vec![GenericEventGroup(vec![GenericEvent {
            weight: Complex::new_re(small.clone()),
            ..Default::default()
        }])]);
        // This is the same single owner used by selected and summed physical
        // maps and reference moments, before any reporting conversion.
        graph.apply_sampling_factor(large.clone());
        assert!((&graph.integrand_result.re - &one).abs() < one.epsilon() * one.from_usize(4));
        assert_eq!(graph.event_groups[0][0].weight, graph.integrand_result);
        assert_eq!(
            graph.reference_moments.as_ref().unwrap().second_moment,
            graph.integrand_result.re
        );
        let level = super::PreciseStabilityLevelResult {
            result: graph.integrand_result.clone(),
            graph_result: graph,
            stability_level_used: Precision::Arb,
            estimated_relative_accuracy: None,
            sample_count: 1,
            total_time: Duration::ZERO,
            parameterization_jacobian: Some(one.clone()),
            is_stable: true,
            rotated_results: Vec::new(),
        };
        let mut metadata = EvaluationMetaData::new_empty();
        metadata.canonical_sampling_preparation_time = Duration::from_micros(7);
        metadata.parameterization_time = Duration::from_micros(7) + Duration::from_micros(11);
        metadata.integrand_evaluation_time = Duration::from_micros(19) + Duration::from_micros(23);
        let result = super::finalize_precise_evaluation_result(level.clone(), F(1.0), metadata);
        assert_eq!(
            result.evaluation_metadata.parameterization_time,
            Duration::from_micros(18)
        );
        assert_eq!(
            result
                .evaluation_metadata
                .canonical_sampling_preparation_time,
            Duration::from_micros(7)
        );
        assert_eq!(
            result.evaluation_metadata.integrand_evaluation_time,
            Duration::from_micros(42)
        );
        assert_eq!(
            result.clone().try_into_f64().unwrap().integrand_result.re,
            F(1.0)
        );
        assert_eq!(result.event_groups[0][0].weight, result.integrand_result);
        // Cancellation can leave a finite signed result while the physical
        // absolute channel sum exceeds the ordinary reporting range.
        let mut excessive_absolute = result.clone();
        excessive_absolute.absolute_integrand_result = Some(Complex::new_re(large.clone()));
        assert!(excessive_absolute.try_into_f64().is_err());
        for extreme in [large.clone(), small.clone()] {
            let mut extreme_level = level.clone();
            extreme_level.result = Complex::new_re(extreme.clone());
            extreme_level.graph_result.integrand_result = extreme_level.result.clone();
            let precise = PreciseEvaluationResult::Arb(super::finalize_precise_evaluation_result(
                extreme_level,
                F(1.0),
                EvaluationMetaData::new_empty(),
            ));
            let PreciseEvaluationResult::Arb(value) = &precise else {
                unreachable!()
            };
            assert_eq!(value.integrand_result.re, extreme);
            assert!(!value.evaluation_metadata.is_nan);
            if extreme == large {
                assert!(
                    format!("{:#}", precise.try_into_f64().unwrap_err())
                        .contains("f64 integration/reporting boundary")
                );
            } else {
                assert_eq!(precise.try_into_f64().unwrap().integrand_result.re, F(0.0));
            }
        }
        // Correct final rounding is allowed for signed and absolute values,
        // including fully weighted event totals, while precise output stays native.
        let mut underflow = result.clone();
        underflow.integrand_result = Complex::new(small.clone(), -small.clone());
        underflow.absolute_integrand_result = Some(Complex::new(small.clone(), small.clone()));
        underflow.event_groups[0][0].weight = underflow.integrand_result.clone();
        let rounded = underflow.clone().try_into_f64().unwrap();
        assert_eq!(rounded.integrand_result, Complex::new_zero());
        assert_eq!(rounded.absolute_integrand_result, Some(Complex::new_zero()));
        assert_eq!(rounded.event_groups[0][0].weight, Complex::new_zero());
        assert_eq!(underflow.integrand_result.re, small);

        let mut subnormal = result.clone();
        subnormal.integrand_result = Complex::new_re(one.from_usize(10).powi(320).inv());
        let subnormal_value = subnormal.integrand_result.re.into_ff64();
        assert!(subnormal_value.0 > 0.0 && subnormal_value.0 < f64::MIN_POSITIVE);
        assert_eq!(
            subnormal.try_into_f64().unwrap().integrand_result.re,
            subnormal_value
        );

        // An outer-grid factor can rescue an otherwise underflowed intermediate;
        // reject that loss rather than silently returning a zero observation.
        underflow.integrator_weight = one.from_usize(10).powi(100);
        let error = format!("{:#}", underflow.clone().try_into_f64().unwrap_err());
        assert!(error.contains("integrand_result.re"));
        assert!(error.contains("complete contribution"));
        underflow.integrand_result = Complex::new_re(one.zero());
        assert!(
            format!("{:#}", underflow.clone().try_into_f64().unwrap_err())
                .contains("absolute_integrand_result.re")
        );
        underflow.absolute_integrand_result = None;
        // Event totals are already weighted; do not apply the outer weight twice.
        assert_eq!(
            underflow.clone().try_into_f64().unwrap().event_groups[0][0].weight,
            Complex::new_zero()
        );
        underflow.integrator_weight = one.clone();
        underflow.parameterization_jacobian = Some(one.from_usize(10).powi(100));
        underflow.integrand_result = Complex::new_im(-small.clone());
        assert!(
            format!("{:#}", underflow.clone().try_into_f64().unwrap_err())
                .contains("integrand_result.im")
        );

        // Separately reported factors must themselves survive conversion: a
        // final underflow does not excuse an infinite serialized Jacobian.
        let mut unrepresentable_factor = result.clone();
        unrepresentable_factor.integrand_result = Complex::new_re(&small * &small);
        unrepresentable_factor.parameterization_jacobian = Some(large.clone());
        assert!(
            format!(
                "{:#}",
                unrepresentable_factor.clone().try_into_f64().unwrap_err()
            )
            .contains("parameterization_jacobian")
        );
        unrepresentable_factor.parameterization_jacobian = Some(small.clone());
        assert!(
            format!(
                "{:#}",
                unrepresentable_factor.clone().try_into_f64().unwrap_err()
            )
            .contains("parameterization_jacobian")
        );
        unrepresentable_factor.parameterization_jacobian = None;
        unrepresentable_factor.integrator_weight = large.clone();
        assert!(
            format!(
                "{:#}",
                unrepresentable_factor.clone().try_into_f64().unwrap_err()
            )
            .contains("integrator_weight")
        );
        unrepresentable_factor.integrator_weight = small.clone();
        assert!(
            format!("{:#}", unrepresentable_factor.try_into_f64().unwrap_err())
                .contains("integrator_weight")
        );

        // Auxiliary event entries deliberately remain factorized. Their raw
        // components cannot be rounded away before the stored multiplier acts.
        let mut auxiliary = result.clone();
        auxiliary.event_groups[0][0]
            .additional_weights
            .weights
            .insert(
                crate::observables::AdditionalWeightKey::Original,
                Complex::new_re(small.clone()),
            );
        assert!(
            format!("{:#}", auxiliary.try_into_f64().unwrap_err())
                .contains("additional_weights[Original].re")
        );

        // A cancelling total does not make individually unrepresentable event
        // weights representable. Native APIs retain them; ordinary output errors.
        let mut cancelling_events = result;
        cancelling_events.event_groups[0].0 = vec![
            GenericEvent {
                weight: Complex::new_re(large.clone()),
                ..Default::default()
            },
            GenericEvent {
                weight: Complex::new_re(-large),
                ..Default::default()
            },
        ];
        assert!(cancelling_events.try_into_f64().is_err());
    }

    #[test]
    fn absolute_channel_instability_cannot_hide_behind_signed_cancellation() {
        use crate::settings::runtime::StabilityLevelSetting;
        use spenso::algebra::complex::Complex;

        let mut level = StabilityLevelSetting::default_double();
        level.ecm_relative_tolerance_for_re = 0.1;
        level.ecm_relative_tolerance_for_im = 0.1;
        let signed = [Complex::new_zero(), Complex::new_zero()];
        let absolute = [Complex::new(F(2.0), F(4.0)), Complex::new(F(3.0), F(6.0))];
        for check in [
            super::stability_check::<f64>,
            super::stability_check_on_norm::<f64>,
        ] {
            assert!(
                check(
                    Some(&F(1.0)),
                    &signed,
                    &level,
                    Complex::new_zero(),
                    F(1.0),
                    false,
                    false
                )
                .2
            );
            assert!(
                !check(
                    Some(&F(1.0)),
                    &absolute,
                    &level,
                    Complex::new_zero(),
                    F(1.0),
                    false,
                    false
                )
                .2
            );
        }
    }

    #[test]
    fn explicit_orientation_sum_rejects_runtime_filters_and_ltd() {
        let mut settings = RuntimeSettings::default();
        settings.general.orientation_pat = OrientationPattern::from_user_pattern("(+)").unwrap();
        let error = validate_process_runtime_settings(&settings, true).unwrap_err();
        assert!(error.to_string().contains("orientation_pat` must be unset"));

        validate_process_runtime_settings(&settings, false).unwrap();

        settings.general.orientation_pat = OrientationPattern::default();
        settings.general.use_ltd = true;
        let error = validate_process_runtime_settings(&settings, false).unwrap_err();
        assert!(error.to_string().contains("deferred proper-LTD support"));
    }

    #[test]
    fn numerator_sampling_scale_requires_nonzero_runtime_value() {
        let mut settings = RuntimeSettings::default();
        validate_process_runtime_settings(&settings, false).unwrap();

        for scale in [0.0, -0.0] {
            settings.general.numerator_sampling_scale = scale;
            for explicit_orientation_sum_only in [false, true] {
                let error =
                    validate_process_runtime_settings(&settings, explicit_orientation_sum_only)
                        .unwrap_err();
                assert!(error.to_string().contains("sampling scale M"));
            }
        }

        settings.general.numerator_sampling_scale = -2.0;
        validate_process_runtime_settings(&settings, false).unwrap();
    }

    #[test]
    fn grouped_orientation_sampling_requires_identical_exact_map_catalogs() {
        let master = ["O[+0]|M[0]", "O[+0]|M[1]"];
        validate_orientation_catalog_group(
            GroupId(3),
            [
                ("master".to_string(), master.to_vec()),
                ("matching".to_string(), master.to_vec()),
            ],
        )
        .unwrap();

        let error = validate_orientation_catalog_group(
            GroupId(3),
            [
                ("master".to_string(), master.to_vec()),
                ("different".to_string(), vec!["O[+0]|M[0]", "O[+0]|M[2]"]),
            ],
        )
        .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("different exact residue-map catalogs")
        );
        assert!(error.to_string().contains("first difference at channel 1"));
    }

    #[test]
    fn runtime_cache_serializes_as_empty() {
        let encoded = bincode::encode_to_vec(
            RuntimeCache::<usize>::default(),
            bincode::config::standard(),
        )
        .expect("runtime cache should encode");
        assert!(encoded.is_empty());

        let (decoded, consumed): (RuntimeCache<usize>, usize) =
            bincode::decode_from_slice(&encoded, bincode::config::standard())
                .expect("runtime cache should decode");
        assert_eq!(consumed, 0);
        assert!(decoded.as_ref().is_none());
        let mut source = RuntimeCache::default();
        source.set((super::Precision::Quad, 1.0e-13));
        let encoded_source = bincode::encode_to_vec(&source, bincode::config::standard()).unwrap();
        assert!(encoded_source.is_empty());
        let (decoded_source, consumed): (RuntimeCache<(super::Precision, f64)>, usize) =
            bincode::decode_from_slice(&encoded_source, bincode::config::standard()).unwrap();
        assert_eq!(consumed, 0);
        assert!(decoded_source.as_ref().is_none());
        // Cached preparation failures are as transient as compiled bridges;
        // neither payload needs a codec or may enter a saved state.
        let mut failure: RuntimeCache<
            Result<super::SamplingChannelBridge, super::sampling_maps::SamplingEvaluationError>,
        > = RuntimeCache::default();
        failure.set(Err(
            super::sampling_maps::SamplingEvaluationError::UncertainGeometry {
                detail: "numeric binding fixture".into(),
            },
        ));
        assert!(
            bincode::encode_to_vec(&failure, bincode::config::standard())
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn filtered_orientation_helpers_map_visible_indices_into_subset_order() {
        let orientations = TiVec::<OrientationID, EdgeVec<Orientation>>::from_iter([
            EdgeVec::from_iter([Orientation::Default]),
            EdgeVec::from_iter([Orientation::Reversed]),
            EdgeVec::from_iter([Orientation::Undirected]),
            EdgeVec::from_iter([Orientation::Default]),
        ]);
        let mut filter = SubSet::empty(orientations.len());
        filter.add(OrientationID(1));
        filter.add(OrientationID(3));

        assert_eq!(filtered_orientation_count(&filter, &orientations), 2);
        assert_eq!(
            resolve_visible_orientation_id(&filter, 0),
            Some(OrientationID(1))
        );
        assert_eq!(
            resolve_visible_orientation_id(&filter, 1),
            Some(OrientationID(3))
        );
        assert_eq!(resolve_visible_orientation_id(&filter, 2), None);
    }

    #[test]
    fn discrete_acceptance_selection_decodes_canonical_channel_id() {
        let settings = SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
            sample_orientations: true,
            sampling_type: DiscreteGraphSamplingType::SamplingMultiChanneling(
                MultiChannelingSettings::default(),
            ),
            ..Default::default()
        });

        let selected = super::resolve_discrete_selection_for_sampling(
            &settings,
            &[2, 5, 7],
            3,
            |_| Some(6),
            |_| Ok::<Option<usize>, eyre::Report>(Some(8)),
        )
        .unwrap();
        assert_eq!(selected.0, Some(GroupId(2)));
        assert_eq!(selected.1, Some(5));
        assert_eq!(selected.2, Some(SamplingChannelId::from(7)));

        let error = super::resolve_discrete_selection_for_sampling(
            &settings,
            &[2, 5, 8],
            3,
            |_| Some(6),
            |_| Ok::<Option<usize>, eyre::Report>(Some(8)),
        )
        .unwrap_err();
        assert!(error.to_string().contains("Channel 8 is out of range"));
    }

    #[test]
    fn sampling_channel_basis_ids_use_graph_override_or_optimized_channels() {
        test_initialise().unwrap();
        static GRAPH: OnceLock<Graph> = OnceLock::new();
        let graph = GRAPH
            .get_or_init(|| {
                dot!(
                    digraph lmb_basis_selection {
                        edge [num=1 mass=0]
                        node [num=1]
                        A -> B [id=0]
                        A -> B [id=1]
                        A -> B [id=2]
                    }
                )
                .unwrap()
            })
            .clone();
        let lmb = |edge_id| LoopMomentumBasis {
            tree: graph.underlying.empty_subgraph(),
            loop_edges: vec![EdgeIndex::from(edge_id)].into(),
            ext_edges: Vec::new().into(),
            edge_signatures: graph.underlying.new_edgevec(|_, _, _| {
                LoopExtSignature::from((Vec::<isize>::new(), Vec::<isize>::new()))
            }),
        };
        let all_bases = vec![lmb(0), lmb(1), lmb(2)].into();
        let setup = LmbMultiChannelingSetup {
            master_edge_masses: Default::default(),
            sampling_bridge: Default::default(),
            sampling_bridge_quad: Default::default(),
            sampling_bridge_fixed256: Default::default(),
            sampling_bridge_arb: Default::default(),
            sampling_source: Default::default(),
            sampling_catalogue: Default::default(),
            sampling_programs: Default::default(),
            lmb_basis_ids: vec![LmbIndex::from(2), LmbIndex::from(0)].into(),
            graph,
            all_bases,
        };
        let default_settings = ParameterizationSettings::default();
        let override_settings = ParameterizationSettings {
            lmb_basis_ids: std::collections::BTreeMap::from([(setup.graph.name.clone(), vec![1])]),
            ..Default::default()
        };
        let out_of_range_settings = ParameterizationSettings {
            lmb_basis_ids: std::collections::BTreeMap::from([(setup.graph.name.clone(), vec![3])]),
            ..Default::default()
        };

        assert_eq!(
            setup
                .selected_lmb_basis_id(&setup.graph.name, &default_settings)
                .unwrap(),
            LmbIndex::from(2)
        );
        assert_eq!(
            setup
                .selected_lmb_basis_id(&setup.graph.name, &override_settings)
                .unwrap(),
            LmbIndex::from(1)
        );
        assert_eq!(
            setup
                .sampling_channel_ids(&setup.graph.name, &override_settings)
                .unwrap(),
            vec![
                SamplingChannelId::from(0),
                SamplingChannelId::from(1),
                SamplingChannelId::from(2)
            ]
        );
        // The canonical IDs retain their one-domain ordering while resolving
        // to the generated LMB basis shown in diagnostics.
        assert_eq!(
            [
                SamplingChannelId::from(0),
                SamplingChannelId::from(1),
                SamplingChannelId::from(2),
            ]
            .into_iter()
            .map(|channel_id| {
                setup
                    .sampling_channel_lmb_id(channel_id, &setup.graph.name, &override_settings)
                    .unwrap()
            })
            .collect::<Vec<_>>(),
            vec![
                Some(LmbIndex::from(1)),
                Some(LmbIndex::from(0)),
                Some(LmbIndex::from(2)),
            ]
        );
        assert_eq!(
            setup
                .sampling_channel_lmb_basis_id(
                    SamplingChannelId::from(0),
                    &setup.graph.name,
                    &override_settings,
                )
                .unwrap(),
            LmbIndex::from(1)
        );
        assert_eq!(
            setup
                .sampling_channel_edge_ids(
                    SamplingChannelId::from(0),
                    &setup.graph.name,
                    &override_settings,
                )
                .unwrap()
                .as_slice(),
            &[1]
        );
        assert!(
            setup
                .selected_lmb_basis_id(&setup.graph.name, &out_of_range_settings)
                .is_err()
        );

        let mut physical_settings = ParameterizationSettings::default();
        let parent_lmb = setup
            .graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .map(|edge| edge.0)
            .collect::<Vec<_>>();
        physical_settings
            .sampling_channels
            .default_channel_selection = vec!["cut".into()];
        physical_settings
            .sampling_channels
            .channel_definitions
            .entry(setup.graph.name.clone())
            .or_default()
            .insert(
                "cut".into(),
                SamplingChannelDefinition {
                    channel_weight: None,
                    around: "phase_space(cut(0))".into(),
                    subspace_lmb: parent_lmb.clone(),
                    parent_lmb,
                    on_cut: vec![],
                    singularity_proxy: None,
                    radial_profile: None,
                },
            );
        let catalogue = setup
            .canonical_sampling_catalogue(&setup.graph.name, &physical_settings)
            .unwrap();
        let dimension = 3 * setup.graph.get_loop_number();
        let programs = catalogue
            .compile_programs(dimension, &HFunctionSettings::default())
            .unwrap();
        let context = SamplingChannelCompileContext::<f64>::new(
            setup.graph.name.clone(),
            setup
                .graph
                .loop_momentum_basis
                .loop_edges
                .iter()
                .map(|edge| edge.0)
                .collect(),
            physical_settings,
            1.0,
            setup.graph.get_loop_number(),
        );
        // Physical labels alone do not provide their native cut geometry.
        assert!(
            setup
                .compile_sampling_channel_bridge(&catalogue, &programs, &context)
                .is_err()
        );

        // Graph-aware entries occupy the same canonical channel axis as LMB
        // entries.  In particular, inserting a named channel before the
        // generated LMBs must not make channel id 1 resolve as basis 1 by
        // positional filtering.
        let mut mixed_settings = ParameterizationSettings::default();
        mixed_settings.sampling_channels.default_channel_selection =
            vec!["named".into(), "auto:lmb".into()];
        mixed_settings
            .sampling_channels
            .channel_definitions
            .entry(setup.graph.name.clone())
            .or_default()
            .insert(
                "named".into(),
                SamplingChannelDefinition {
                    channel_weight: None,
                    around: "lmb(0)".into(),
                    subspace_lmb: Vec::new(),
                    parent_lmb: vec![0],
                    on_cut: Vec::new(),
                    singularity_proxy: None,
                    radial_profile: None,
                },
            );
        assert_eq!(
            setup
                .sampling_channel_ids(&setup.graph.name, &mixed_settings)
                .unwrap(),
            vec![
                SamplingChannelId::from(0),
                SamplingChannelId::from(1),
                SamplingChannelId::from(2),
                SamplingChannelId::from(3),
            ]
        );
        assert!(
            setup
                .sampling_channel_lmb_basis_id(
                    SamplingChannelId::from(0),
                    &setup.graph.name,
                    &mixed_settings,
                )
                .is_err()
        );
        assert_eq!(
            setup
                .sampling_channel_lmb_id(
                    SamplingChannelId::from(1),
                    &setup.graph.name,
                    &mixed_settings,
                )
                .unwrap(),
            Some(LmbIndex::from(0))
        );
    }

    #[test]
    fn canonical_lmb_density_partition_matches_inverse_jacobians() {
        test_initialise().unwrap();
        let mut graph: Graph = dot!(
            digraph lmb_prefactor_partition {
                edge [num=1 mass=0]
                node [num=1]
                ext [style=invis]
                ext -> A [id=0]
                A -> B [id=1]
                A -> B [id=2]
                B -> ext [id=3]
            }
        )
        .unwrap();
        let all_bases = graph.generate_loop_momentum_bases();
        assert!(all_bases.len() >= 2);
        graph.loop_momentum_basis = all_bases[LmbIndex::from(0)].clone();
        let mut setup = LmbMultiChannelingSetup {
            master_edge_masses: Default::default(),
            sampling_bridge: Default::default(),
            sampling_bridge_quad: Default::default(),
            sampling_bridge_fixed256: Default::default(),
            sampling_bridge_arb: Default::default(),
            sampling_source: Default::default(),
            sampling_catalogue: Default::default(),
            sampling_programs: Default::default(),
            lmb_basis_ids: vec![LmbIndex::from(0), LmbIndex::from(1)].into(),
            graph: graph.clone(),
            all_bases,
        };
        let loop_moms: LoopMomenta<F<f64>> = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .enumerate()
            .map(|(index, _)| {
                let offset = index as f64;
                ThreeMomentum::new(F(0.4 + offset), F(-0.3), F(0.2 - offset))
            })
            .collect();
        let external_moms: ExternalFourMomenta<F<f64>> =
            (0..graph.loop_momentum_basis.ext_edges.len())
                .map(|_| [F(0.0), F(0.2), F(-0.3), F(0.4)].into())
                .collect();
        let sample = MomentumSample {
            sample: BareMomentumSample {
                loop_moms,
                dual_loop_moms: None,
                loop_mom_cache_id: 0,
                loop_mom_base_cache_id: 0,
                external_moms,
                external_mom_cache_id: 0,
                external_mom_base_cache_id: 0,
                jacobian: F(1.0),
                orientation: None,
                parameterization_branch: None,
            },
        };
        let mut parameterization_settings = ParameterizationSettings::default();
        parameterization_settings
            .sampling_channels
            .default_channel_selection = vec!["auto:lmb".into()];
        let resolved = resolve_sampling_channel_selection(
            &graph.name,
            &parameterization_settings.sampling_channels,
        )
        .unwrap();
        let context = SamplingChannelCompileContext::new(
            graph.name.clone(),
            graph
                .loop_momentum_basis
                .loop_edges
                .iter()
                .map(|edge| edge.0)
                .collect(),
            parameterization_settings.clone(),
            1.0,
            sample.loop_moms().0.len(),
        );
        let external = sample
            .external_moms()
            .iter()
            .map(|momentum| {
                [
                    momentum.temporal.value.0,
                    momentum.spatial.px.0,
                    momentum.spatial.py.0,
                    momentum.spatial.pz.0,
                ]
            })
            .collect::<Vec<_>>();
        let bridge = setup
            .compile_sampling_channel_bridge_with_external(
                &setup
                    .sampling_channel_catalogue(&resolved, &context.parameterization_settings)
                    .unwrap(),
                &setup
                    .sampling_channel_catalogue(&resolved, &context.parameterization_settings)
                    .unwrap()
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
                &context,
                &external,
            )
            .unwrap();
        let raw = sample
            .loop_moms()
            .iter()
            .flat_map(|momentum| [momentum.px.0, momentum.py.0, momentum.pz.0])
            .collect::<Vec<_>>();
        // Check map-density weights against independent inverse parameterizations
        // in every generated LMB at the same parent momentum point. This retains
        // the former prefactor partition's pointwise normalization test without
        // keeping a separate LMB weighting implementation in production.
        let scores = setup
            .all_bases
            .iter()
            .map(|basis| {
                let momenta = basis
                    .loop_edges
                    .iter()
                    .map(|edge| {
                        let signature = &graph.loop_momentum_basis.edge_signatures[*edge];
                        signature.internal.apply_typed(sample.loop_moms())
                            + signature
                                .external
                                .apply(&sample.external_moms().raw)
                                .spatial
                    })
                    .collect::<Vec<_>>();
                crate::utils::global_inv_parameterize(&momenta, F(1.0), &parameterization_settings)
                    .1
                    .0
            })
            .collect::<Vec<_>>();
        let total = scores.iter().sum::<f64>();
        for channel_id in 0..bridge.channels().len() {
            let evaluation = bridge
                .inverse(SamplingChannelId::from(channel_id), &raw)
                .unwrap()
                .expect("full-support channel inverse");
            assert!((evaluation.partition.weights.iter().sum::<f64>() - 1.0).abs() < 1.0e-14);
            let actual = evaluation.partition.weight(channel_id).unwrap();
            let expected = scores[channel_id] / total;
            assert!(
                (actual - expected).abs() < 1.0e-14,
                "channel {channel_id}: {actual} != {expected}"
            );
        }
        // The master owns OSE masses even if a member/context supplies different ones.
        // Refresh twice to catch a stale runtime cache after a model/mass change.
        let model = crate::utils::load_generic_model("sm");
        let mut context = context;
        context.parameterization_settings.sampling_channels.weight =
            crate::settings::runtime::SamplingChannelWeight::Ose;
        context.edge_masses = graph
            .new_edgevec(|_, _, _| spenso::algebra::complex::Complex::new_re(F(1000.0)))
            .iter()
            .map(|(edge, mass)| (edge.0, *mass))
            .collect();
        for mass in [2.0, 5.0] {
            let edge = *setup.all_bases[LmbIndex::from(0)]
                .loop_edges
                .first()
                .unwrap();
            setup.graph.underlying[edge].mass = crate::graph::edge::EdgeMass::Value(
                spenso::algebra::complex::Complex::new_re(F(mass)),
            );
            setup.warm_up_masses(&RuntimeSettings::default(), &model);
            assert_eq!(setup.master_edge_masses.as_ref().unwrap()[edge].re, F(mass));
            assert_ne!(graph.get_real_mass_vector::<f64>(&model)[edge], F(mass));
            let catalogue = setup
                .sampling_channel_catalogue(&resolved, &context.parameterization_settings)
                .unwrap();
            let programs = catalogue
                .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                .unwrap();
            let bridge = setup
                .compile_sampling_channel_bridge_with_external(
                    &catalogue, &programs, &context, &external,
                )
                .unwrap();
            let energies = setup.graph.get_energy_cache(
                &model,
                sample.loop_moms(),
                sample.external_moms(),
                &setup.graph.loop_momentum_basis,
            );
            let scores = setup
                .all_bases
                .iter()
                .map(|basis| {
                    basis
                        .loop_edges
                        .iter()
                        .map(|edge| energies[*edge].0)
                        .product::<f64>()
                        .powf(-3.0)
                })
                .collect::<Vec<_>>();
            let partition = bridge.partition(&raw).unwrap();
            let total = scores.iter().sum::<f64>();
            for (actual, score) in partition.weights.iter().zip(scores) {
                assert!((actual - score / total).abs() < 1e-13);
            }
        }
    }

    #[test]
    fn process_catalogue_bridge_integrates_normalized_gaussian() {
        test_initialise().unwrap();
        let mut graph: Graph = dot!(
            digraph process_sampling_acceptance {
                edge [num=1 mass=0]
                node [num=1]
                ext [style=invis]
                ext -> A [id=0]
                A -> B [id=1]
                A -> B [id=2]
                B -> ext [id=3]
            }
        )
        .unwrap();
        let generated_bases = graph.generate_loop_momentum_bases();
        assert!(!generated_bases.is_empty());
        graph.loop_momentum_basis = generated_bases[LmbIndex::from(0)].clone();
        let all_bases = vec![graph.loop_momentum_basis.clone()].into();
        let parent_lmb = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .map(|edge| edge.0)
            .collect::<Vec<_>>();
        let setup = LmbMultiChannelingSetup {
            master_edge_masses: Default::default(),
            sampling_bridge: Default::default(),
            sampling_bridge_quad: Default::default(),
            sampling_bridge_fixed256: Default::default(),
            sampling_bridge_arb: Default::default(),
            sampling_source: Default::default(),
            sampling_catalogue: Default::default(),
            sampling_programs: Default::default(),
            lmb_basis_ids: vec![LmbIndex::from(0)].into(),
            graph: graph.clone(),
            all_bases,
        };

        let selection = SamplingChannelSelection::default();
        let resolved = resolve_sampling_channel_selection(&graph.name, &selection).unwrap();
        let mut parameterization_settings = ParameterizationSettings::default();
        parameterization_settings.sampling_channels = selection;
        let context = SamplingChannelCompileContext::new(
            graph.name.clone(),
            parent_lmb,
            parameterization_settings,
            100.0,
            graph.loop_momentum_basis.loop_edges.len(),
        );
        let bridge = setup
            .compile_sampling_channel_bridge(
                &setup
                    .sampling_channel_catalogue(&resolved, &context.parameterization_settings)
                    .unwrap(),
                &setup
                    .sampling_channel_catalogue(&resolved, &context.parameterization_settings)
                    .unwrap()
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
                &context,
            )
            .unwrap();
        assert_eq!(
            setup
                .sampling_channel_ids(&graph.name, &context.parameterization_settings)
                .unwrap()
                .len(),
            bridge.channels().len(),
            "omitted defaults must resolve identically for channel counts and compiled maps"
        );
        assert_eq!(bridge.channels().len(), 1);

        let dimensions = bridge.dimensions();
        let width = 1.5;
        let normalisation =
            (2.0 * std::f64::consts::PI * width * width).powf(-0.5 * dimensions as f64);
        let samples = 2048usize;
        let mut integral = 0.0;
        for sample in 1..=samples {
            let coordinates = (0..dimensions)
                .map(|axis| {
                    let base = [2_u64, 3, 5, 7, 11, 13, 17, 19][axis];
                    let mut index = sample;
                    let mut fraction = 1.0;
                    let mut value = 0.0;
                    while index > 0 {
                        fraction /= base as f64;
                        value += fraction * (index as u64 % base) as f64;
                        index /= base as usize;
                    }
                    value
                })
                .collect::<Vec<_>>();
            let evaluation = bridge
                .forward(SamplingChannelId::from(0), &coordinates)
                .unwrap();
            assert!((evaluation.partition.weights.iter().sum::<f64>() - 1.0).abs() < 1.0e-14);
            let radius_squared = evaluation
                .raw_coordinates
                .iter()
                .map(|component| component * component)
                .sum::<f64>();
            integral += normalisation
                * (-0.5 * radius_squared / width.powi(2)).exp()
                * evaluation.map.jacobian;
        }
        integral /= samples as f64;
        assert!(
            (integral - 1.0).abs() < 5.0e-2,
            "process sampling bridge Gaussian integral = {integral}"
        );
    }
}
