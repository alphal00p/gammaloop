use std::fs;
use std::path::Path;

use crate::graph::{FeynmanGraph, Graph, GraphGroup, GroupId, LmbIndex, LoopMomentumBasis};
use crate::integrands::evaluation::{
    EvaluationMetaData, EvaluationResult, GenericEvaluationResult, GraphEvaluationResult,
    LoopMomentaEscalationMetrics, PreciseEvaluationResult, RawBatchEvaluationResult,
    RawPreciseBatchEvaluationResult, RotatedEvaluation, StabilityFailureReason, StabilityResult,
    StabilityStatus, StatisticsCounter,
};
use crate::model::Model;
use crate::momentum::sample::{BareMomentumSample, LoopMomenta, MomentumSample};
use crate::momentum::{Rotation, ThreeMomentum};
use crate::observables::{
    AdditionalWeightKey, EventProcessingRuntime, ObservableAccumulatorBundle, ObservableFileFormat,
    ObservableSnapshotBundle,
};
use crate::processes::StandaloneExportSettings;
use crate::utils::{
    ArbPrec, F, FloatLike, f128, format_for_compare_digits, get_n_dim_for_n_loop_momenta,
};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::owo_colors::OwoColorize;
use colored::Colorize;
use derive_more::{From, Into};
use enum_dispatch::enum_dispatch;
use eyre::{Context, eyre};
use gammaloop_sample::{DiscreteGraphSample, GammaLoopSample, parameterize};
use itertools::Itertools;
use linnet::half_edge::involution::EdgeVec;
use momtrop::SampleGenerator;
use momtrop::float::MomTropFloat;
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
use crate::{
    DependentMomentaConstructor, GammaLoopContext, settings::RuntimeSettings,
    settings::runtime::DiscreteGraphSamplingSettings, settings::runtime::DiscreteGraphSamplingType,
    settings::runtime::IntegratorSettings, settings::runtime::Precision,
    settings::runtime::SamplingSettings, settings::runtime::StabilityLevelSetting,
    settings::runtime::StabilitySettings,
};
use color_eyre::Result;

pub mod evaluators;
pub use evaluators::ActiveF64Backend;
pub use evaluators::{GenericEvaluator, GenericEvaluatorFloat};

pub mod param_builder;
pub use param_builder::{ParamBuilder, ParamValuePairs, ThresholdParams, UpdateAndGetParams};

#[derive(Debug, Clone)]
pub struct MomentumSpaceEvaluationInput {
    pub loop_momenta: Vec<ThreeMomentum<F<f64>>>,
    pub integrator_weight: F<f64>,
    pub graph_id: Option<usize>,
    pub group_id: Option<GroupId>,
    pub orientation: Option<usize>,
    pub channel_id: Option<ChannelIndex>,
}

#[derive(Clone, Debug)]
pub(crate) struct RuntimeCache<T>(Option<T>);

impl<T> Default for RuntimeCache<T> {
    fn default() -> Self {
        Self(None)
    }
}

impl<T> RuntimeCache<T> {
    pub(crate) fn invalidate(&mut self) {
        self.0 = None;
    }

    pub(crate) fn set(&mut self, value: T) {
        self.0 = Some(value);
    }

    pub(crate) fn take(&mut self) -> Option<T> {
        self.0.take()
    }

    pub(crate) fn as_ref(&self) -> Option<&T> {
        self.0.as_ref()
    }

    pub(crate) fn as_mut(&mut self) -> Option<&mut T> {
        self.0.as_mut()
    }
}

impl<T> bincode::Encode for RuntimeCache<T> {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        _encoder: &mut E,
    ) -> std::result::Result<(), bincode::error::EncodeError> {
        Ok(())
    }
}

impl<C, T> bincode::Decode<C> for RuntimeCache<T> {
    fn decode<D: bincode::de::Decoder<Context = C>>(
        _decoder: &mut D,
    ) -> std::result::Result<Self, bincode::error::DecodeError> {
        Ok(Self::default())
    }
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
        DiscreteGraphSamplingType::DiscreteMultiChanneling(_) => "discrete_multi_channeling",
    }
}

pub(crate) fn discrete_sampling_depth_for_settings(
    settings: &DiscreteGraphSamplingSettings,
) -> usize {
    let orientation_depth = usize::from(settings.sample_orientations);
    match &settings.sampling_type {
        DiscreteGraphSamplingType::DiscreteMultiChanneling(_) => 2 + orientation_depth,
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
        DiscreteGraphSamplingType::DiscreteMultiChanneling(_)
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
    mut channel_count_for_group: impl FnMut(GroupId) -> Option<usize>,
) -> Result<(Option<GroupId>, Option<usize>, Option<ChannelIndex>)> {
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
                DiscreteGraphSamplingType::DiscreteMultiChanneling(_) => {
                    let channel_index = *discrete_dimensions.last().expect("validated depth");
                    let channel_count = channel_count_for_group(group_id).ok_or_else(|| {
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
                    Some(ChannelIndex::from(channel_index))
                }
                _ => None,
            };

            Ok((Some(group_id), orientation, channel))
        }
    }
}

impl ProcessIntegrand {
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

    pub fn get_mut_settings(&mut self) -> &mut RuntimeSettings {
        match self {
            ProcessIntegrand::Amplitude(integrand) => {
                integrand.invalidate_event_processing_runtime();
                &mut integrand.settings
            }
            ProcessIntegrand::CrossSection(integrand) => {
                integrand.invalidate_event_processing_runtime();
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
    ) -> Result<(Option<GroupId>, Option<usize>, Option<ChannelIndex>)> {
        let group_count = match self {
            ProcessIntegrand::Amplitude(integrand) => integrand.data.graph_group_structure.len(),
            ProcessIntegrand::CrossSection(integrand) => integrand.data.graph_group_structure.len(),
        };

        resolve_discrete_selection_for_sampling(
            &self.get_settings().sampling,
            discrete_dimensions,
            group_count,
            |group_id| self.group_orientation_count(group_id),
            |group_id| self.group_channel_count(group_id),
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

    pub fn group_channel_count(&self, group_id: GroupId) -> Option<usize> {
        match self {
            ProcessIntegrand::Amplitude(integrand) => Some(
                integrand.data.graph_terms[integrand.data.graph_group_structure[group_id].master()]
                    .get_num_channels(),
            ),
            ProcessIntegrand::CrossSection(integrand) => Some(
                integrand.data.graph_terms[integrand.data.graph_group_structure[group_id].master()]
                    .get_num_channels(),
            ),
        }
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
        model: &Model,
        samples: &[Sample<F<f64>>],
        iter: usize,
        use_arb_prec: bool,
        stop_on_interrupt: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<RawBatchEvaluationResult> {
        let mut results = Vec::with_capacity(samples.len());
        for sample in samples {
            if stop_on_interrupt && crate::is_interrupted() {
                break;
            }
            let mut result = match self {
                ProcessIntegrand::Amplitude(integrand) => evaluate_sample(
                    integrand,
                    model,
                    sample,
                    sample.get_weight(),
                    iter,
                    use_arb_prec,
                    max_eval,
                ),
                ProcessIntegrand::CrossSection(integrand) => evaluate_sample(
                    integrand,
                    model,
                    sample,
                    sample.get_weight(),
                    iter,
                    use_arb_prec,
                    max_eval,
                ),
            }?;

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

fn full_event_multiplicative_factor(
    parameterization_jacobian: Option<F<f64>>,
    integrator_weight: F<f64>,
) -> Complex<F<f64>> {
    let jacobian = parameterization_jacobian.unwrap_or(F(1.0));
    Complex::new_re(jacobian * integrator_weight)
}

fn apply_full_event_multiplicative_factor(
    event_groups: &mut crate::observables::EventGroupList,
    full_factor: &Complex<F<f64>>,
) {
    for event_group in event_groups.iter_mut() {
        for event in event_group.iter_mut() {
            event.weight *= full_factor;
            if !event.additional_weights.weights.is_empty() {
                event
                    .additional_weights
                    .weights
                    .entry(AdditionalWeightKey::FullMultiplicativeFactor)
                    .and_modify(|value| *value *= full_factor)
                    .or_insert_with(|| *full_factor);
            }
        }
    }
}

fn full_event_multiplicative_factor_precise<T: FloatLike>(
    parameterization_jacobian: Option<F<T>>,
    integrator_weight: F<T>,
) -> Complex<F<T>> {
    let jacobian = parameterization_jacobian.unwrap_or_else(|| integrator_weight.one());
    Complex::new_re(jacobian * integrator_weight)
}

fn apply_full_event_multiplicative_factor_precise<T: FloatLike>(
    event_groups: &mut crate::observables::GenericEventGroupList<T>,
    full_factor: &Complex<F<T>>,
) {
    for event_group in event_groups.iter_mut() {
        for event in event_group.iter_mut() {
            event.weight *= full_factor.clone();

            if !event.additional_weights.weights.is_empty() {
                event.additional_weights.weights.insert(
                    AdditionalWeightKey::FullMultiplicativeFactor,
                    full_factor.clone(),
                );
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
            vec![StabilityLevelSetting {
                precision: Precision::Arb,
                required_precision_for_re: 1e-5,
                required_precision_for_im: 1e-5,
                escalate_for_large_weight_threshold: -1.,
            }]
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

#[inline]
fn stability_check<T: FloatLike>(
    _settings: &RuntimeSettings,
    results: &[Complex<F<T>>],
    stability_settings: &StabilityLevelSetting,
    max_eval: Complex<F<T>>,
    wgt: F<T>,
    is_final_level: bool,
    escalate_if_exact_zero: bool,
) -> (
    Complex<F<T>>,
    Option<F<T>>,
    bool,
    Option<StabilityFailureReason>,
) {
    if results.len() == 1 {
        return (results[0].clone(), None, true, None);
    }

    let average = results
        .iter()
        .skip(1)
        .fold(results[0].clone(), |acc, x| acc + x)
        / F::<T>::from_f64(results.len() as f64);

    let errors = results.iter().map(|res| {
        let error_re = if IsZero::is_zero(&res.re) && IsZero::is_zero(&average.re) {
            F::<T>::from_f64(0.0)
        } else {
            ((&res.re - &average.re) / &average.re).abs()
        };
        let error_im = if IsZero::is_zero(&res.im) && IsZero::is_zero(&average.im) {
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
            && error.re == F::<T>::from_f64(0.0)
            && error.im == F::<T>::from_f64(0.0)
        {
            unstable_reason = Some(StabilityFailureReason::ZeroError);
            unstable_sample = Some(index);
            break;
        }

        if error.re > F::<T>::from_f64(stability_settings.required_precision_for_re)
            || error.im > F::<T>::from_f64(stability_settings.required_precision_for_im)
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

    let below_wgt_threshold = if stability_settings.escalate_for_large_weight_threshold > 0.
        && max_eval.is_non_zero()
    {
        average.re.abs() * wgt.clone()
            < F::<T>::from_f64(stability_settings.escalate_for_large_weight_threshold) * max_eval.re
            || average.im.abs() * wgt
                < F::<T>::from_f64(stability_settings.escalate_for_large_weight_threshold)
                    * max_eval.im
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
fn stability_check_on_norm<T: FloatLike>(
    _settings: &RuntimeSettings,
    results: &[Complex<F<T>>],
    stability_settings: &StabilityLevelSetting,
    max_eval: Complex<F<T>>,
    wgt: F<T>,
    is_final_level: bool,
    escalate_if_exact_zero: bool,
) -> (
    Complex<F<T>>,
    Option<F<T>>,
    bool,
    Option<StabilityFailureReason>,
) {
    if results.len() == 1 {
        return (results[0].clone(), None, true, None);
    }

    let average = results.iter().fold(F::<T>::from_f64(0.0), |acc, x| {
        acc + x.norm_squared().sqrt()
    }) / F::<T>::from_f64(results.len() as f64);

    let errors = results.iter().map(|res| {
        let res = res.norm_squared().sqrt();
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

        if error > F::<T>::from_f64(stability_settings.required_precision_for_re) {
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

    let below_wgt_threshold =
        if stability_settings.escalate_for_large_weight_threshold > 0. && max_eval.is_non_zero() {
            average.abs() * wgt
                < F::<T>::from_f64(stability_settings.escalate_for_large_weight_threshold)
                    * max_eval.norm_squared().sqrt()
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
pub struct StabilityLevelResult {
    pub result: Complex<F<f64>>,
    pub graph_result: GraphEvaluationResult<f64>,
    pub stability_level_used: Precision,
    pub estimated_relative_accuracy: Option<F<f64>>,
    pub sample_count: usize,
    pub total_time: Duration,
    pub parameterization_time: Duration,
    pub parameterization_jacobian: Option<F<f64>>,
    pub integrand_evaluation_time: Duration,
    pub evaluator_evaluation_time: Duration,
    pub is_stable: bool,
    pub instability_reason: Option<StabilityFailureReason>,
    pub rotated_results: Vec<RotatedEvaluation>,
}

#[derive(Debug, Clone)]
struct PreciseStabilityLevelResult<T: FloatLike> {
    pub result: Complex<F<T>>,
    pub graph_result: GraphEvaluationResult<T>,
    pub stability_level_used: Precision,
    pub estimated_relative_accuracy: Option<F<T>>,
    pub sample_count: usize,
    pub total_time: Duration,
    pub parameterization_time: Duration,
    pub parameterization_jacobian: Option<F<T>>,
    pub integrand_evaluation_time: Duration,
    pub evaluator_evaluation_time: Duration,
    pub is_stable: bool,
    pub instability_reason: Option<StabilityFailureReason>,
    pub rotated_results: Vec<RotatedEvaluation>,
}

impl<T: FloatLike> PreciseStabilityLevelResult<T> {
    fn into_f64(self) -> StabilityLevelResult {
        StabilityLevelResult {
            result: complex_to_f64(&self.result),
            graph_result: self.graph_result.into_f64(),
            stability_level_used: self.stability_level_used,
            estimated_relative_accuracy: self
                .estimated_relative_accuracy
                .map(|value| value.into_ff64()),
            sample_count: self.sample_count,
            total_time: self.total_time,
            parameterization_time: self.parameterization_time,
            parameterization_jacobian: self
                .parameterization_jacobian
                .map(|value| value.into_ff64()),
            integrand_evaluation_time: self.integrand_evaluation_time,
            evaluator_evaluation_time: self.evaluator_evaluation_time,
            is_stable: self.is_stable,
            instability_reason: self.instability_reason,
            rotated_results: self.rotated_results,
        }
    }
}

#[derive(
    Debug, Clone, Copy, From, Into, PartialEq, Eq, Ord, PartialOrd, Hash, Serialize, Deserialize,
)]
pub struct ChannelIndex(usize);

/// Helper struct for the LMB multi-channeling setup
#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct LmbMultiChannelingSetup {
    pub channels: TiVec<ChannelIndex, LmbIndex>,
    pub graph: Graph,
    pub all_bases: TiVec<LmbIndex, LoopMomentumBasis>,
}

impl LmbMultiChannelingSetup {
    pub(crate) fn channel_edge_ids(&self, channel_index: ChannelIndex) -> SmallVec<[usize; 4]> {
        self.all_bases[self.channels[channel_index]]
            .loop_edges
            .iter()
            .map(|edge_id| edge_id.0)
            .collect()
    }

    fn reinterpret_loop_momenta_impl<T: FloatLike>(
        &self,
        channel_index: ChannelIndex,
        momentum_sample: &BareMomentumSample<T>,
        loop_mom_cache_id: usize,
    ) -> BareMomentumSample<T> {
        let channel_lmb = &self.all_bases[self.channels[channel_index]];
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
        }
    }

    /// Note this increments the loop_mom_cache_id of all of the returned BareMomentumSample
    #[allow(dead_code)]
    pub(crate) fn reinterpret_loop_momenta_and_compute_prefactor_all_channels<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        model: &Model,
        alpha: &F<T>,
        cache: bool,
    ) -> TiVec<ChannelIndex, (MomentumSample<T>, F<T>)> {
        let mut loop_mom_cache_id = momentum_sample.sample.loop_mom_cache_id;
        self.channels
            .iter_enumerated()
            .map(|(channel_index, _)| {
                if cache {
                    loop_mom_cache_id += 1;
                }
                self.reinterpret_loop_momenta_and_compute_prefactor(
                    channel_index,
                    momentum_sample,
                    loop_mom_cache_id,
                    model,
                    alpha,
                )
            })
            .collect()
    }

    /// This function is used to do do LMB multi-channeling without fully switching to a different lmb
    /// for each channel. The momenta provided are reinterpreted as loop momenta of the lmb corresponding to the channel_index.
    /// Then we transform these loop momenta to the fixed lmb of the graph. The prefactor is immediately computed for the requested channel
    ///
    /// Note this increments the loop_mom_cache_id of the returned BareMomentumSample
    pub(crate) fn reinterpret_loop_momenta_and_compute_prefactor<T: FloatLike>(
        &self,
        channel_index: ChannelIndex,
        momentum_sample: &MomentumSample<T>,
        loop_mom_cache_id: usize,
        model: &Model,
        alpha: &F<T>,
    ) -> (MomentumSample<T>, F<T>) {
        let sample = MomentumSample {
            sample: self.reinterpret_loop_momenta_impl(
                channel_index,
                &momentum_sample.sample,
                loop_mom_cache_id,
            ), // uuid: momentum_sample.uuid,
        };

        let prefactor = self.compute_prefactor_impl(channel_index, &sample, model, alpha);

        (sample, prefactor)
    }

    /// Computes the prefactor for the given channel index and momentum sample.
    fn compute_prefactor_impl<T: FloatLike>(
        &self,
        channel_index: ChannelIndex,
        momentum_sample: &MomentumSample<T>,
        model: &Model,
        alpha: &F<T>,
    ) -> F<T> {
        let all_energies = self.graph.get_energy_cache(
            model,
            &momentum_sample.sample.loop_moms,
            &momentum_sample.sample.external_moms,
            &self.graph.loop_momentum_basis,
        );

        let mut numerator = momentum_sample.zero();

        let denominators = self
            .channels
            .iter()
            .map(|&lmb_index| {
                let channel_product = self.all_bases[lmb_index]
                    .loop_edges
                    .iter()
                    .map(|&edge_index| &all_energies[edge_index])
                    .fold(momentum_sample.one(), |product, energy| product * energy)
                    .powf(&-alpha);

                if self.channels[channel_index] == lmb_index {
                    numerator = channel_product.clone();
                }

                channel_product
            })
            .fold(momentum_sample.zero(), |sum, summand| sum + summand);

        numerator / denominators
    }
}

pub trait ProcessIntegrandImpl {
    type G: GraphTerm;

    fn warm_up(&mut self, model: &Model) -> Result<()>;
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
    /// let new_sample = parameterize(&sample_point, &mut integrand)?;
    ///
    /// // Rotations will automatically get new cache IDs
    /// let rotated_samples = evaluate_all_rotations(&new_sample, &mut integrand, true)?;
    ///
    /// // Next iteration - revert to base configuration to reuse cache
    /// integrand.revert_to_base_external_cache_id();
    /// let next_sample = parameterize(&next_sample_point, &mut integrand)?;
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
    fn get_graph_mut(&mut self, graph_id: usize) -> &mut Self::G;
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

    // fn get_builder_cache(&self) -> &ParamBuilder<f64>;
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
    fn evaluate<T: FloatLike>(
        &mut self,
        sample: &MomentumSample<T>,
        model: &Model,
        settings: &RuntimeSettings,
        event_processing_runtime: Option<&mut EventProcessingRuntime>,
        rotation: &Rotation,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
        channel_id: Option<(ChannelIndex, F<T>)>,
    ) -> Result<GraphEvaluationResult<T>>;

    fn name(&self) -> String;

    fn warm_up(&mut self, settings: &RuntimeSettings, model: &Model) -> Result<()>;
    fn get_graph(&self) -> &Graph;
    fn get_num_channels(&self) -> usize;
    fn get_num_orientations(&self) -> usize;
    fn get_tropical_sampler(&self) -> &SampleGenerator<3>;
    fn get_mut_param_builder(&mut self) -> &mut ParamBuilder<f64>;
    fn get_real_mass_vector(&self) -> EdgeVec<Option<F<f64>>>;
}

fn evaluate_graph_term<T: FloatLike, I: ProcessIntegrandImpl>(
    integrand: &mut I,
    graph_id: usize,
    sample: &MomentumSample<T>,
    model: &Model,
    settings: &RuntimeSettings,
    rotation: &Rotation,
    evaluation_metadata: &mut EvaluationMetaData,
    record_primary_timing: bool,
    channel_id: Option<(ChannelIndex, F<T>)>,
) -> Result<GraphEvaluationResult<T>> {
    let mut event_processing_runtime = integrand.take_event_processing_runtime();
    let result = integrand.get_graph_mut(graph_id).evaluate(
        sample,
        model,
        settings,
        event_processing_runtime.as_mut(),
        rotation,
        evaluation_metadata,
        record_primary_timing,
        channel_id,
    );
    integrand.restore_event_processing_runtime(event_processing_runtime);
    let mut result = result?;
    for event_group in result.event_groups.iter_mut() {
        for event in event_group.iter_mut() {
            event.cut_info.graph_id = graph_id;
        }
    }
    Ok(result)
}

fn evaluate_graph_group<T: FloatLike, I: ProcessIntegrandImpl>(
    integrand: &mut I,
    group_id: GroupId,
    sample: &DiscreteGraphSample<T>,
    model: &Model,
    settings: &RuntimeSettings,
    rotation: &Rotation,
    evaluation_metadata: &mut EvaluationMetaData,
    record_primary_timing: bool,
    zero: &F<T>,
) -> Result<GraphEvaluationResult<T>> {
    let group = integrand.get_group(group_id).into_iter().collect_vec();

    let mut result = GraphEvaluationResult::zero(zero.clone());
    let mut grouped_events = crate::observables::GenericEventGroup::default();

    for graph_id in group {
        let graph_term_result = match sample {
            DiscreteGraphSample::Default(sample) => evaluate_graph_term(
                integrand,
                graph_id,
                sample,
                model,
                settings,
                rotation,
                evaluation_metadata,
                record_primary_timing,
                None,
            ),
            DiscreteGraphSample::DiscreteMultiChanneling {
                alpha,
                channel_id,
                sample,
            } => evaluate_graph_term(
                integrand,
                graph_id,
                sample,
                model,
                settings,
                rotation,
                evaluation_metadata,
                record_primary_timing,
                Some((*channel_id, alpha.clone())),
            ),
            DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                let num_channels = integrand.get_master_graph(group_id).get_num_channels();
                (0..num_channels)
                    .map(ChannelIndex::from)
                    .map(|channel_index| {
                        evaluate_graph_term(
                            integrand,
                            graph_id,
                            sample,
                            model,
                            settings,
                            rotation,
                            evaluation_metadata,
                            record_primary_timing,
                            Some((channel_index, alpha.clone())),
                        )
                    })
                    .try_fold(
                        GraphEvaluationResult::zero(zero.clone()),
                        |mut sum, term| {
                            sum.merge_in_place(term?);
                            Ok::<GraphEvaluationResult<T>, eyre::Report>(sum)
                        },
                    )
            }
            DiscreteGraphSample::Tropical(sample) => {
                let master_graph = integrand.get_master_graph(group_id).get_graph();

                let energy_cache = master_graph.get_energy_cache(
                    model,
                    sample.loop_moms(),
                    sample.external_moms(),
                    &master_graph.loop_momentum_basis,
                );

                let prefactor = master_graph
                    .iter_loop_edges()
                    .map(|(_, edge_index, _)| edge_index)
                    .zip(
                        integrand
                            .get_master_graph(group_id)
                            .get_tropical_sampler()
                            .iter_edge_weights(),
                    )
                    .fold(sample.one(), |product, (edge_id, weight)| {
                        let energy = &energy_cache[edge_id];
                        product * energy.powf(&F::from_f64(2. * weight))
                    });

                let mut graph_result = evaluate_graph_term(
                    integrand,
                    graph_id,
                    sample,
                    model,
                    settings,
                    rotation,
                    evaluation_metadata,
                    record_primary_timing,
                    None,
                )?;
                graph_result.integrand_result *= Complex::new_re(prefactor);
                Ok(graph_result)
            }
        }?;

        let mut graph_term_result = graph_term_result;
        for mut event_group in graph_term_result.event_groups.drain(..) {
            grouped_events.append(&mut event_group);
        }
        result.merge_in_place(graph_term_result);
    }

    if !grouped_events.is_empty() {
        result.event_groups.push(grouped_events);
    }

    Ok(result)
}

fn evaluate_all_rotations<T: FloatLike, I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    gammaloop_sample: &GammaLoopSample<T>,
    evaluation_metadata: &mut EvaluationMetaData,
    is_primary_stability_level: bool,
    record_rotated_results: bool,
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
                loop_mom_cache_id += 1;
                external_mom_cache_id += 1;
            }
            gammaloop_sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id)
        })
        .collect();

    let primary_rotation_index = rotations
        .iter()
        .position(Rotation::is_identity)
        .unwrap_or(0);
    let mut original_call_timed = false;
    let mut evaluation_results: Vec<GraphEvaluationResult<T>> =
        Vec::with_capacity(gammaloop_samples.len());
    for (rotation_index, (gammaloop_sample, rotation)) in
        gammaloop_samples.iter().zip(rotations.iter()).enumerate()
    {
        debug!("Evaluating rotation: {}", rotation.method);
        let record_primary_timing = is_primary_stability_level
            && !original_call_timed
            && rotation_index == primary_rotation_index;

        let result = evaluate_single(
            integrand,
            model,
            gammaloop_sample,
            rotation,
            evaluation_metadata,
            record_primary_timing,
        )?;

        if record_primary_timing {
            original_call_timed = true;
        }
        evaluation_results.push(result);
    }

    for result in &evaluation_results {
        evaluation_metadata.event_processing_time += result.event_processing_time;
    }

    if cache {
        integrand.increment_loop_cache_id(rotations.len());
        integrand.increment_external_cache_id(rotations.len());
        // After evaluating all rotations, revert to base cache ID to enable cache reuse
        // for subsequent sample points with the same base external momenta
        integrand.revert_to_base_external_cache_id();
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

fn evaluate_stability_level_precise<T: FloatLike, I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    source: &EvaluationSource<'_>,
    stability_level: &StabilityLevelSetting,
    max_eval: &Complex<F<f64>>,
    wgt: F<f64>,
    check_on_norm: bool,
    is_final_level: bool,
    is_primary_stability_level: bool,
    evaluation_metadata: &mut EvaluationMetaData,
    record_rotated_results: bool,
    precision_label: &str,
    escalate_if_exact_zero: bool,
) -> Result<PreciseStabilityLevelResult<T>> {
    let level_start = Instant::now();
    let (gammaloop_sample, parameterization_time) = source.build_gamma_sample::<T, I>(integrand)?;
    debug!("{precision_label} parameterization succeeded");
    debug!(
        "jacobian: {:+16e}",
        gammaloop_sample.get_default_sample().jacobian()
    );

    let integrand_time_before = evaluation_metadata.integrand_evaluation_time;
    let evaluator_time_before = evaluation_metadata.evaluator_evaluation_time;
    let (graph_results, primary_rotation_index, rotated_results) = evaluate_all_rotations(
        integrand,
        model,
        &gammaloop_sample,
        evaluation_metadata,
        is_primary_stability_level,
        record_rotated_results,
    )?;
    let results = graph_results
        .iter()
        .map(|result| result.integrand_result.clone())
        .collect_vec();

    let max_eval = complex_from_f64::<T>(max_eval);
    let wgt = F::<T>::from_ff64(wgt);

    let (average_result, estimated_relative_accuracy, is_stable, instability_reason) =
        if check_on_norm {
            stability_check_on_norm(
                integrand.get_settings(),
                &results,
                stability_level,
                max_eval,
                wgt,
                is_final_level,
                escalate_if_exact_zero,
            )
        } else {
            stability_check(
                integrand.get_settings(),
                &results,
                stability_level,
                max_eval,
                wgt,
                is_final_level,
                escalate_if_exact_zero,
            )
        };

    let mut graph_result = graph_results[primary_rotation_index].clone();
    graph_result.integrand_result = average_result.clone();

    Ok(PreciseStabilityLevelResult {
        result: average_result,
        graph_result,
        stability_level_used: stability_level.precision,
        estimated_relative_accuracy,
        sample_count: results.len(),
        total_time: level_start.elapsed(),
        parameterization_time,
        parameterization_jacobian: match source {
            EvaluationSource::XSpace(_) => Some(gammaloop_sample.get_default_sample().jacobian()),
            EvaluationSource::Momentum(_) => None,
        },
        integrand_evaluation_time: evaluation_metadata
            .integrand_evaluation_time
            .saturating_sub(integrand_time_before),
        evaluator_evaluation_time: evaluation_metadata
            .evaluator_evaluation_time
            .saturating_sub(evaluator_time_before),
        is_stable,
        instability_reason,
        rotated_results,
    })
}

fn evaluate_stability_level<T: FloatLike, I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    source: &EvaluationSource<'_>,
    stability_level: &StabilityLevelSetting,
    max_eval: &Complex<F<f64>>,
    wgt: F<f64>,
    check_on_norm: bool,
    is_final_level: bool,
    is_primary_stability_level: bool,
    evaluation_metadata: &mut EvaluationMetaData,
    record_rotated_results: bool,
    precision_label: &str,
    escalate_if_exact_zero: bool,
) -> Result<StabilityLevelResult> {
    evaluate_stability_level_precise::<T, I>(
        integrand,
        model,
        source,
        stability_level,
        max_eval,
        wgt,
        check_on_norm,
        is_final_level,
        is_primary_stability_level,
        evaluation_metadata,
        record_rotated_results,
        precision_label,
        escalate_if_exact_zero,
    )
    .map(PreciseStabilityLevelResult::into_f64)
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
    model: &Model,
    gammaloop_sample: &GammaLoopSample<T>,
    rotation: &Rotation,
    evaluation_metadata: &mut EvaluationMetaData,
    record_primary_timing: bool,
) -> Result<GraphEvaluationResult<T>> {
    let settings = integrand.get_settings().clone();
    let zero = gammaloop_sample.get_default_sample().zero();
    let loop_cache_shift = 0;
    let cache = integrand.get_settings().general.enable_cache;

    let start_integrand_timing = if record_primary_timing {
        Some(std::time::Instant::now())
    } else {
        None
    };
    let result = (|| -> Result<GraphEvaluationResult<T>> {
        let result = match &gammaloop_sample {
            GammaLoopSample::Default(sample) => {
                if integrand.groups_default_sample_events_by_graph_group() {
                    integrand
                        .get_group_structure()
                        .iter_enumerated()
                        .map(|(group_id, _)| group_id)
                        .collect_vec()
                        .into_iter()
                        .try_fold(
                            GraphEvaluationResult::zero(zero.clone()),
                            |mut sum, group_id| {
                                let group_result = evaluate_graph_group(
                                    integrand,
                                    group_id,
                                    &DiscreteGraphSample::Default(sample.clone()),
                                    model,
                                    &settings,
                                    rotation,
                                    evaluation_metadata,
                                    record_primary_timing,
                                    &zero,
                                )?;
                                sum.merge_in_place(group_result);
                                Ok::<GraphEvaluationResult<T>, eyre::Report>(sum)
                            },
                        )?
                } else {
                    (0..integrand.graph_count()).try_fold(
                        GraphEvaluationResult::zero(zero.clone()),
                        |mut sum, graph_id| {
                            let graph_result = evaluate_graph_term(
                                integrand,
                                graph_id,
                                sample,
                                model,
                                &settings,
                                rotation,
                                evaluation_metadata,
                                record_primary_timing,
                                None,
                            )?;
                            sum.merge_in_place(graph_result);
                            Ok::<GraphEvaluationResult<T>, eyre::Report>(sum)
                        },
                    )?
                }
            }
            GammaLoopSample::Graph { graph_id, sample } => evaluate_graph_term(
                integrand,
                *graph_id,
                sample,
                model,
                &settings,
                rotation,
                evaluation_metadata,
                record_primary_timing,
                None,
            )?,
            GammaLoopSample::MultiChanneling { .. } => {
                unimplemented!(
                    "deprecated due to annyoing borrow issues, just set each graph to the same group"
                );
            }
            GammaLoopSample::DiscreteGraph { group_id, sample } => evaluate_graph_group(
                integrand,
                *group_id,
                sample,
                model,
                &settings,
                rotation,
                evaluation_metadata,
                record_primary_timing,
                &zero,
            )?,
        };

        if cache {
            integrand.increment_loop_cache_id(loop_cache_shift);
        }

        Ok(result)
    })();
    if record_primary_timing {
        evaluation_metadata.integrand_evaluation_time = evaluation_metadata
            .integrand_evaluation_time
            .saturating_add(
                start_integrand_timing
                    .expect("integrand timing start should exist")
                    .elapsed(),
            );
    }

    result
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
        DiscreteGraphSamplingType::DiscreteMultiChanneling(_multichanneling_settings) => {
            let continuous_grid = create_default_continous_grid(graph_term, integrator_settings);
            let lmb_channel_grid = Grid::Discrete(DiscreteGrid::new(
                (0..graph_term.get_num_channels())
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
}

impl<'a> EvaluationSource<'a> {
    fn build_gamma_sample<T: FloatLike, I: ProcessIntegrandImpl>(
        &self,
        integrand: &mut I,
    ) -> Result<(GammaLoopSample<T>, Duration)> {
        match self {
            EvaluationSource::XSpace(sample) => {
                let before_parameterization = std::time::Instant::now();
                let sample = parameterize::<T, I>(sample, integrand)?;
                Ok((sample, before_parameterization.elapsed()))
            }
            EvaluationSource::Momentum(input) => Ok((
                build_direct_gamma_sample::<T, I>(integrand, input)?,
                Duration::ZERO,
            )),
        }
    }

    fn loop_norm_sum<I: ProcessIntegrandImpl>(&self, integrand: &mut I) -> Result<F<f64>> {
        match self {
            EvaluationSource::XSpace(sample) => {
                let sample = parameterize::<f64, I>(sample, integrand)?;
                Ok(sum_loop_norms(
                    sample.get_default_sample().loop_moms().0.iter(),
                ))
            }
            EvaluationSource::Momentum(input) => Ok(sum_loop_norms(input.loop_momenta.iter())),
        }
    }

    fn debug_sample<I: ProcessIntegrandImpl>(
        &self,
        integrand: &mut I,
    ) -> Result<GammaLoopSample<f64>> {
        self.build_gamma_sample::<f64, I>(integrand)
            .map(|(sample, _)| sample)
    }
}

fn sum_loop_norms<'a>(loop_momenta: impl Iterator<Item = &'a ThreeMomentum<F<f64>>>) -> F<f64> {
    loop_momenta.fold(F(0.0), |acc, momentum| acc + momentum.norm())
}

fn build_direct_gamma_sample<T: FloatLike, I: ProcessIntegrandImpl>(
    integrand: &mut I,
    input: &MomentumSpaceEvaluationInput,
) -> Result<GammaLoopSample<T>> {
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

    let loop_momenta = input
        .loop_momenta
        .iter()
        .map(|momentum| {
            ThreeMomentum::new(
                F::<T>::from_ff64(momentum.px),
                F::<T>::from_ff64(momentum.py),
                F::<T>::from_ff64(momentum.pz),
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
        return Ok(GammaLoopSample::Graph { graph_id, sample });
    }

    match &integrand.get_settings().sampling {
        SamplingSettings::Default(_) | SamplingSettings::MultiChanneling(_) => {
            if input.group_id.is_some() || input.channel_id.is_some() {
                return Err(eyre!(
                    "Discrete graph/channel selections are not supported for this sampling mode."
                ));
            }

            Ok(GammaLoopSample::Default(sample))
        }
        SamplingSettings::DiscreteGraphs(settings) => {
            let Some(group_id) = input.group_id else {
                if input.orientation.is_some() || input.channel_id.is_some() {
                    return Err(eyre!(
                        "Explicit orientation or channel selections require selecting a graph group in momentum-space evaluation."
                    ));
                }
                return Ok(GammaLoopSample::Default(sample));
            };
            let discrete_sample = match &settings.sampling_type {
                DiscreteGraphSamplingType::Default(_) => {
                    if input.channel_id.is_some() {
                        return Err(eyre!(
                            "Channel selection is not available for this discrete-graph sampling mode."
                        ));
                    }
                    DiscreteGraphSample::Default(sample)
                }
                DiscreteGraphSamplingType::MultiChanneling(multichanneling_settings) => {
                    if input.channel_id.is_some() {
                        return Err(eyre!(
                            "Channel selection is not available for this discrete-graph sampling mode."
                        ));
                    }
                    DiscreteGraphSample::MultiChanneling {
                        alpha: F::from_f64(multichanneling_settings.alpha),
                        sample,
                    }
                }
                DiscreteGraphSamplingType::TropicalSampling(_) => {
                    if input.channel_id.is_some() {
                        return Err(eyre!(
                            "Channel selection is not available for tropical discrete-graph sampling."
                        ));
                    }
                    DiscreteGraphSample::Tropical(sample)
                }
                DiscreteGraphSamplingType::DiscreteMultiChanneling(multichanneling_settings) => {
                    let channel_id = input.channel_id.ok_or_else(|| {
                        eyre!(
                            "Momentum-space evaluation for discrete multichanneling requires selecting a channel."
                        )
                    })?;
                    DiscreteGraphSample::DiscreteMultiChanneling {
                        alpha: F::from_f64(multichanneling_settings.alpha),
                        channel_id,
                        sample,
                    }
                }
            };

            Ok(GammaLoopSample::DiscreteGraph {
                group_id,
                sample: discrete_sample,
            })
        }
    }
}

fn log_rotated_samples<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    gammaloop_sample: &GammaLoopSample<f64>,
    level_result: &StabilityLevelResult,
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
            loop_mom_cache_id += 1;
            shift += 1;
            external_mom_cache_id += 1;
            gammaloop_sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id)
        })
        .collect();
    integrand.increment_external_cache_id(shift);
    integrand.increment_loop_cache_id(shift);

    for (sample, result) in rotated_samples
        .iter()
        .zip(level_result.rotated_results.iter())
    {
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

fn evaluate_from_source<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    source: EvaluationSource<'_>,
    wgt: F<f64>,
    use_arb_prec: bool,
    max_eval: Complex<F<f64>>,
) -> Result<EvaluationResult> {
    let start_eval = std::time::Instant::now();
    let mut escalate_if_exact_zero = integrand.get_settings().stability.escalate_if_exact_zero;
    if escalate_if_exact_zero && !integrand.get_settings().selectors.is_empty() {
        warn_selectors_disable_zero_once();
        escalate_if_exact_zero = false;
    }
    let mut evaluation_metadata = EvaluationMetaData::new_empty();
    let (stability_iterator, loop_momenta_escalation) =
        stability_iterator_for_source(integrand, &source, use_arb_prec);

    let mut results_of_stability_levels = Vec::with_capacity(stability_iterator.len());

    let total_levels = stability_iterator.len();
    for (level_index, stability_level) in stability_iterator.into_iter().enumerate() {
        let is_final_level = level_index + 1 == total_levels;
        let record_rotated_results = integrand
            .get_settings()
            .stability
            .recording
            .map(|recording| recording.record_rotated_results)
            .unwrap_or(false);
        let is_primary_stability_level = level_index == 0;
        let result_of_level = match stability_level.precision {
            Precision::Double => evaluate_stability_level::<f64, I>(
                integrand,
                model,
                &source,
                &stability_level,
                &max_eval,
                wgt,
                integrand.get_settings().stability.check_on_norm,
                is_final_level,
                is_primary_stability_level,
                &mut evaluation_metadata,
                record_rotated_results,
                "f64",
                escalate_if_exact_zero,
            ),
            Precision::Quad => evaluate_stability_level::<f128, I>(
                integrand,
                model,
                &source,
                &stability_level,
                &max_eval,
                wgt,
                integrand.get_settings().stability.check_on_norm,
                is_final_level,
                is_primary_stability_level,
                &mut evaluation_metadata,
                record_rotated_results,
                "f128",
                escalate_if_exact_zero,
            ),
            Precision::Arb => evaluate_stability_level::<ArbPrec, I>(
                integrand,
                model,
                &source,
                &stability_level,
                &max_eval,
                wgt,
                integrand.get_settings().stability.check_on_norm,
                is_final_level,
                is_primary_stability_level,
                &mut evaluation_metadata,
                record_rotated_results,
                "ArbPrec",
                escalate_if_exact_zero,
            ),
        }?;

        let is_stable = result_of_level.is_stable;
        results_of_stability_levels.push(result_of_level);

        if is_stable {
            break;
        } else {
            debug!("unstable at level: {}", stability_level.precision);
            if let Ok(gammaloop_sample) = source.debug_sample(integrand) {
                let level_result = results_of_stability_levels
                    .last()
                    .expect("stability level result missing");
                log_rotated_samples(integrand, &gammaloop_sample, level_result);
            } else {
                debug!("failed to reconstruct sample for instability logging");
            }
        }
    }

    debug!("result at each level:");
    for level_result in results_of_stability_levels.iter() {
        debug!(
            "level: {}. result: {}",
            format!("{}", level_result.stability_level_used).green(),
            format!("{:16e}", level_result.result).blue()
        );
    }

    if let Some(stability_level_result) = results_of_stability_levels.last().cloned() {
        let re_is_nan = stability_level_result.result.re.is_nan()
            || stability_level_result.result.re.is_infinite();
        let im_is_nan = stability_level_result.result.im.is_nan()
            || stability_level_result.result.im.is_infinite();
        let is_nan = re_is_nan || im_is_nan;

        let stability_results = results_of_stability_levels
            .iter()
            .map(|level| StabilityResult {
                precision: level.stability_level_used,
                estimated_relative_accuracy: level.estimated_relative_accuracy,
                status: StabilityStatus::from_sample_count(level.sample_count, level.is_stable),
                total_time: level.total_time,
            })
            .collect();

        evaluation_metadata.total_timing = start_eval.elapsed();
        evaluation_metadata.parameterization_time = stability_level_result.parameterization_time;
        evaluation_metadata.generated_event_count =
            stability_level_result.graph_result.generated_event_count;
        evaluation_metadata.accepted_event_count =
            stability_level_result.graph_result.accepted_event_count;
        evaluation_metadata.relative_instability_error = Complex::new_zero();
        evaluation_metadata.is_nan = is_nan;
        evaluation_metadata.loop_momenta_escalation = loop_momenta_escalation;
        evaluation_metadata.stability_results = stability_results;

        let nanless_result = if re_is_nan && !im_is_nan {
            Complex::new(F(0.0), stability_level_result.result.im)
        } else if im_is_nan && !re_is_nan {
            Complex::new(stability_level_result.result.re, F(0.0))
        } else if im_is_nan && re_is_nan {
            Complex::new(F(0.0), F(0.0))
        } else {
            stability_level_result.result
        };
        let mut event_groups = stability_level_result.graph_result.event_groups;
        let parameterization_jacobian = stability_level_result.parameterization_jacobian;
        let full_factor = full_event_multiplicative_factor(parameterization_jacobian, wgt);
        apply_full_event_multiplicative_factor(&mut event_groups, &full_factor);
        Ok(EvaluationResult {
            integrand_result: nanless_result,
            parameterization_jacobian,
            integrator_weight: wgt,
            event_groups,
            evaluation_metadata,
        })
    } else {
        println!("Evaluation failed at all stability levels");
        Ok(EvaluationResult {
            integrand_result: Complex::new(F(0.0), F(0.0)),
            parameterization_jacobian: None,
            integrator_weight: wgt,
            event_groups: Default::default(),
            evaluation_metadata: EvaluationMetaData {
                total_timing: Duration::ZERO,
                integrand_evaluation_time: Duration::ZERO,
                evaluator_evaluation_time: Duration::ZERO,
                parameterization_time: Duration::ZERO,
                event_processing_time: Duration::ZERO,
                generated_event_count: 0,
                accepted_event_count: 0,
                relative_instability_error: Complex::new(F(0.0), F(0.0)),
                is_nan: true,
                loop_momenta_escalation: None,
                stability_results: Vec::new(),
            },
        })
    }
}

fn stability_iterator_for_source<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    source: &EvaluationSource<'_>,
    use_arb_prec: bool,
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
        && let Ok(sum_norm) = source.loop_norm_sum(integrand)
    {
        let threshold =
            F::<f64>::from_f64(escalation_factor * integrand.get_settings().kinematics.e_cm);
        if record_loop_momenta_escalation {
            loop_momenta_escalation = Some(LoopMomentaEscalationMetrics {
                sum_norm: sum_norm.0,
                threshold: threshold.0,
            });
        }
        if sum_norm > threshold
            && let Some(last) = stability_iterator.last().copied()
        {
            stability_iterator = vec![last];
        }
    }

    (stability_iterator, loop_momenta_escalation)
}

fn finalize_precise_evaluation_result<T: FloatLike>(
    result: PreciseStabilityLevelResult<T>,
    integrator_weight: F<f64>,
    evaluation_metadata: EvaluationMetaData,
) -> GenericEvaluationResult<T> {
    let re_is_nan = result.result.re.is_nan() || result.result.re.is_infinite();
    let im_is_nan = result.result.im.is_nan() || result.result.im.is_infinite();
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
        integrand_result: nanless_result,
        parameterization_jacobian,
        integrator_weight,
        event_groups,
        evaluation_metadata,
    }
}

fn evaluate_from_source_precise<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    source: EvaluationSource<'_>,
    wgt: F<f64>,
    use_arb_prec: bool,
    max_eval: Complex<F<f64>>,
) -> Result<crate::integrands::evaluation::PreciseEvaluationResult> {
    let base_result = evaluate_from_source(integrand, model, source, wgt, use_arb_prec, max_eval)?;
    let final_precision = base_result
        .evaluation_metadata
        .final_precision()
        .unwrap_or(Precision::Double);
    let base_metadata = base_result.evaluation_metadata.clone();
    let mut escalate_if_exact_zero = integrand.get_settings().stability.escalate_if_exact_zero;
    if escalate_if_exact_zero && !integrand.get_settings().selectors.is_empty() {
        escalate_if_exact_zero = false;
    }

    match final_precision {
        Precision::Double => Ok(
            crate::integrands::evaluation::PreciseEvaluationResult::Double(
                GenericEvaluationResult {
                    integrand_result: base_result.integrand_result,
                    parameterization_jacobian: base_result.parameterization_jacobian,
                    integrator_weight: base_result.integrator_weight,
                    event_groups: base_result.event_groups,
                    evaluation_metadata: base_metadata,
                },
            ),
        ),
        Precision::Quad => {
            let (stability_iterator, _) =
                stability_iterator_for_source(integrand, &source, use_arb_prec);
            let stability_level = stability_iterator
                .into_iter()
                .rev()
                .find(|level| level.precision == Precision::Quad)
                .ok_or_else(|| {
                    eyre!(
                        "Quad precision was selected for the final result, but no quad stability level is configured."
                    )
                })?;
            let mut evaluation_metadata = EvaluationMetaData::new_empty();
            let result = evaluate_stability_level_precise::<f128, I>(
                integrand,
                model,
                &source,
                &stability_level,
                &max_eval,
                wgt,
                integrand.get_settings().stability.check_on_norm,
                true,
                true,
                &mut evaluation_metadata,
                false,
                "f128",
                escalate_if_exact_zero,
            )?;
            Ok(
                crate::integrands::evaluation::PreciseEvaluationResult::Quad(
                    finalize_precise_evaluation_result(result, wgt, base_metadata),
                ),
            )
        }
        Precision::Arb => {
            let (stability_iterator, _) =
                stability_iterator_for_source(integrand, &source, use_arb_prec);
            let stability_level = stability_iterator
                .into_iter()
                .rev()
                .find(|level| level.precision == Precision::Arb)
                .ok_or_else(|| {
                    eyre!(
                        "Arbitrary precision was selected for the final result, but no Arb precision stability level is configured."
                    )
                })?;
            let mut evaluation_metadata = EvaluationMetaData::new_empty();
            let result = evaluate_stability_level_precise::<ArbPrec, I>(
                integrand,
                model,
                &source,
                &stability_level,
                &max_eval,
                wgt,
                integrand.get_settings().stability.check_on_norm,
                true,
                true,
                &mut evaluation_metadata,
                false,
                "ArbPrec",
                escalate_if_exact_zero,
            )?;
            Ok(crate::integrands::evaluation::PreciseEvaluationResult::Arb(
                finalize_precise_evaluation_result(result, wgt, base_metadata),
            ))
        }
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
        model,
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
        model,
        EvaluationSource::Momentum(input),
        wgt,
        use_arb_prec,
        max_eval,
    )
}

#[cfg(test)]
mod tests {
    use super::RuntimeCache;

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
        assert!(decoded.0.is_none());
    }
}
