#![allow(dead_code)]

//! This module steers the integration process.
//! It contains two main ways of integrating the integrand.
//! The havana_integrate function is mostly used for local runs.
//! The master node in combination with batch_integrate is for distributed runs.

mod display;
mod render_ratatui;
mod render_tabled;
mod status_update;

use bincode::Decode;
use bincode::Encode;
use color_eyre::{Report, Result};
use colored::Colorize;
use itertools::Itertools;
use itertools::izip;
use linnet::half_edge::involution::{EdgeVec, Orientation};
use rayon::ThreadPoolBuilder;
use serde::Deserialize;
use serde::Serialize;
use spenso::algebra::algebraic_traits::IsZero;
use symbolica::domains::float::Constructible;
use symbolica::numerical_integration::{
    DiscreteGrid, Grid, MonteCarloRng, Sample, StatisticsAccumulator,
};

use crate::Integrand;
use crate::graph::{GroupId, LoopMomentumBasis};
use crate::integrands::HasIntegrand;
use crate::integrands::evaluation::EvaluationResult;
use crate::integrands::evaluation::StatisticsCounter;
use crate::integrands::process::ProcessIntegrand;
use crate::model::{Model, SerializableInputParamCard};
use crate::observables::{EventGroupList, ObservableAccumulatorBundle, ObservableFileFormat};
use crate::settings::IntegratorSettings;
use crate::settings::RuntimeSettings;
use crate::settings::runtime::{
    ComponentDiscreteBreakdown, DiscreteBreakdown, DiscreteBreakdownEntry, DiscreteCoordinate,
    DiscreteGraphSamplingType, IntegralEstimate, IntegratedPhase, IntegrationResult,
    IntegrationTableComponentResult, MaxWeightInfoEntry, SamplingSettings, SlotIntegrationResult,
};
use crate::utils;
use crate::utils::F;
use crate::{
    clear_iteration_abort_request, is_interrupted, is_iteration_abort_requested, set_interrupted,
};
use rayon::prelude::*;
pub use render_ratatui::RatatuiDashboardState;
pub use render_tabled::TabledRenderOptions;
use spenso::algebra::complex::Complex;
pub use status_update::{
    ContributionSortMode, IntegrationStatusKind, IntegrationStatusPhaseDisplay,
    IntegrationStatusViewOptions, StatusUpdate,
};
use status_update::{
    StatusUpdateBuildRequest, build_saved_status_update, build_status_update,
    evaluate_target_accuracy,
};
use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::time::Duration;
use std::time::Instant;
use tabled::Tabled;
#[allow(unused_imports)]
use tracing::{debug, error, info, trace, warn};

#[derive(Clone, Copy, Debug, Default)]
pub struct WorkspaceSnapshotControl {
    pub write_iteration_archives: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct IterationBatchingSettings {
    pub batch_size: Option<usize>,
    pub batch_timing_seconds: f64,
    pub min_time_between_status_updates_seconds: f64,
    pub emit_live_status_updates: bool,
    pub emit_initial_status_update: bool,
}

impl Default for IterationBatchingSettings {
    fn default() -> Self {
        Self {
            batch_size: None,
            batch_timing_seconds: 5.0,
            min_time_between_status_updates_seconds: 0.0,
            emit_live_status_updates: true,
            emit_initial_status_update: true,
        }
    }
}

pub struct IntegrationSlot {
    pub meta: SlotMeta,
    pub settings: RuntimeSettings,
    pub model: Model,
    pub integrand: Integrand,
    pub target: Option<Complex<F<f64>>>,
}

impl IntegrationSlot {
    pub fn new(
        meta: SlotMeta,
        settings: RuntimeSettings,
        model: Model,
        integrand: Integrand,
        target: Option<Complex<F<f64>>>,
    ) -> Self {
        Self {
            meta,
            settings,
            model,
            integrand,
            target,
        }
    }
}

pub struct HavanaIntegrateRequest {
    pub slots: Vec<IntegrationSlot>,
    pub sampling_correlation_mode: SamplingCorrelationMode,
    pub n_cores: usize,
    pub state: Option<IntegrationState>,
    pub workspace: Option<PathBuf>,
    pub output_control: WorkspaceSnapshotControl,
    pub batching: IterationBatchingSettings,
    pub view_options: IntegrationStatusViewOptions,
}

// const N_INTEGRAND_ACCUMULATORS: usize = 2;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, Encode, Decode)]
pub struct SlotMeta {
    pub process_name: String,
    pub integrand_name: String,
}

impl SlotMeta {
    pub fn key(&self) -> String {
        format!("{}@{}", self.process_name, self.integrand_name)
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct IntegrationWorkspaceManifest {
    pub slots: Vec<SlotMeta>,
    pub targets: Vec<Option<Complex<F<f64>>>>,
    pub effective_model_parameters: Vec<SerializableInputParamCard<F<f64>>>,
    pub integrand_fingerprints: Vec<String>,
    pub training_slot: usize,
    pub integrator_settings_slot: usize,
    pub sampling_correlation_mode: SamplingCorrelationMode,
}

impl crate::utils::serde_utils::SmartSerde for IntegrationWorkspaceManifest {}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Encode, Decode)]
#[serde(rename_all = "snake_case")]
pub enum SamplingCorrelationMode {
    Correlated,
    Uncorrelated,
}

impl SamplingCorrelationMode {
    fn state_index(self, slot_index: usize) -> usize {
        match self {
            Self::Correlated => 0,
            Self::Uncorrelated => slot_index,
        }
    }

    fn state_count(self, n_slots: usize) -> usize {
        match self {
            Self::Correlated => 1,
            Self::Uncorrelated => n_slots,
        }
    }
}

#[derive(Serialize, Deserialize, Encode, Decode, Clone)]
struct DiscreteGridAccumulatorSummary {
    #[bincode(with_serde)]
    bins: Vec<DiscreteGridBinAccumulatorSummary>,
}

#[derive(Serialize, Deserialize, Encode, Decode, Clone)]
struct DiscreteGridBinAccumulatorSummary {
    #[bincode(with_serde)]
    accumulator: StatisticsAccumulator<F<f64>>,
    sub_summary: Option<Box<DiscreteGridAccumulatorSummary>>,
}

#[derive(Serialize, Deserialize, Encode, Decode, Clone, Default)]
struct PersistedDiscreteBreakdownMetadata {
    axis_label: String,
    #[bincode(with_serde)]
    fixed_coordinates: Vec<DiscreteCoordinate>,
    bin_labels: Vec<String>,
}

#[derive(Clone, Serialize, Deserialize, Encode, Decode)]
struct SamplingSlotState {
    #[bincode(with_serde)]
    grid: Grid<F<f64>>,
    discrete_axis_labels: Vec<String>,
}

impl SamplingSlotState {
    fn new(grid: Grid<F<f64>>, discrete_axis_labels: Vec<String>) -> Self {
        Self {
            grid,
            discrete_axis_labels,
        }
    }

    fn from_integrand(integrand: &Integrand) -> Self {
        let grid = integrand.create_grid();
        let discrete_axis_labels = match integrand {
            Integrand::ProcessIntegrand(process_integrand) => {
                discrete_axis_labels(&process_integrand.get_settings().sampling)
                    .into_iter()
                    .map(str::to_string)
                    .collect_vec()
            }
            _ => Vec::new(),
        };
        Self::new(grid, discrete_axis_labels)
    }
}

#[derive(Default)]
struct MonitoredDiscreteSetup {
    path: Option<Vec<usize>>,
    axis_label: Option<String>,
    bin_descriptions: Option<Vec<String>>,
    warning: Option<String>,
}

#[derive(Tabled)]
pub struct IntegralResult {
    id: String,
    n_samples: String,
    #[tabled(rename = "n_samples[%]")]
    n_samples_perc: String,
    #[tabled(rename = "<I>")]
    integral: String,
    #[tabled(rename = "sqrt(σ)")]
    variance: String,
    err: String,
    #[tabled(rename = "err[%]")]
    err_perc: String,
    #[tabled(rename = "PDF")]
    pdf: String,
}

#[derive(Serialize, Deserialize, Encode, Decode, Clone)]
pub struct ComplexAccumulator {
    #[bincode(with_serde)]
    pub re: StatisticsAccumulator<F<f64>>,
    #[bincode(with_serde)]
    pub im: StatisticsAccumulator<F<f64>>,
}

impl ComplexAccumulator {
    pub(crate) fn new() -> Self {
        Self {
            re: StatisticsAccumulator::new(),
            im: StatisticsAccumulator::new(),
        }
    }

    /// used to escalate the stability test on high evaluation values
    pub(crate) fn get_worst_case(&self) -> Complex<F<f64>> {
        Complex::new(
            self.re
                .max_eval_positive
                .abs()
                .max(self.re.max_eval_negative.abs()),
            self.im
                .max_eval_positive
                .abs()
                .max(self.im.max_eval_negative.abs()),
        )
    }

    /// add a sample to the accumulator
    pub(crate) fn add_sample(
        &mut self,
        result: Complex<F<f64>>,
        sample_weight: F<f64>,
        sample: Option<&Sample<F<f64>>>,
    ) {
        self.re.add_sample(result.re * sample_weight, sample);
        self.im.add_sample(result.im * sample_weight, sample);
    }

    pub(crate) fn merge(&mut self, other: &Self) {
        self.re.merge_samples_no_reset(&other.re);
        self.im.merge_samples_no_reset(&other.im);
    }

    pub(crate) fn update_iter(&mut self, use_weighted_average: bool) {
        self.re.update_iter(use_weighted_average);
        self.im.update_iter(use_weighted_average);
    }

    pub(crate) fn max_weight_rows(
        &self,
        slot_meta: &SlotMeta,
        discrete_axis_labels: &[String],
    ) -> Vec<[String; 3]> {
        let max_evals = [
            &self.re.max_eval_positive,
            &self.re.max_eval_negative,
            &self.im.max_eval_positive,
            &self.im.max_eval_negative,
        ];

        let max_eval_samples = [
            &self.re.max_eval_positive_xs,
            &self.re.max_eval_negative_xs,
            &self.im.max_eval_positive_xs,
            &self.im.max_eval_negative_xs,
        ];

        let sign_strs = ["+", "-", "+", "-"];
        let phase_strs = ["re", "re", "im", "im"];

        izip!(
            max_evals.iter(),
            max_eval_samples.iter(),
            sign_strs.iter(),
            phase_strs.iter()
        )
        .filter_map(|(max_eval, max_eval_sample, sign_str, phase_str)| {
            if max_eval.is_non_zero() {
                Some([
                    format!(
                        "{} {} [{}] ",
                        slot_label(slot_meta),
                        format!("{:<2}", phase_str).blue(),
                        format!("{:<1}", sign_str).blue()
                    ),
                    format!("{:+.16e}", max_eval),
                    if let Some(sample) = max_eval_sample {
                        format_max_eval_sample(sample, discrete_axis_labels, &[])
                    } else {
                        "N/A".to_string()
                    },
                ])
            } else {
                None
            }
        })
        .collect()
    }
}

impl DiscreteGridAccumulatorSummary {
    fn from_grid(grid: &Grid<F<f64>>) -> Option<Self> {
        match grid {
            Grid::Discrete(discrete_grid) => Some(Self {
                bins: discrete_grid
                    .bins
                    .iter()
                    .map(|bin| DiscreteGridBinAccumulatorSummary {
                        accumulator: StatisticsAccumulator::new(),
                        sub_summary: bin
                            .sub_grid
                            .as_ref()
                            .and_then(Self::from_grid)
                            .map(Box::new),
                    })
                    .collect(),
            }),
            Grid::Continuous(_) | Grid::Uniform(_, _) => None,
        }
    }

    fn merge_iteration_grid(&mut self, grid: &Grid<F<f64>>) {
        let Grid::Discrete(discrete_grid) = grid else {
            return;
        };

        for (summary_bin, grid_bin) in self.bins.iter_mut().zip(discrete_grid.bins.iter()) {
            summary_bin
                .accumulator
                .merge_samples_no_reset(&grid_bin.accumulator);
            if let (Some(sub_summary), Some(sub_grid)) =
                (summary_bin.sub_summary.as_mut(), grid_bin.sub_grid.as_ref())
            {
                sub_summary.merge_iteration_grid(sub_grid);
            }
        }
    }

    fn update_iter(&mut self) {
        for bin in &mut self.bins {
            bin.accumulator.update_iter(false);
            if let Some(sub_summary) = bin.sub_summary.as_mut() {
                sub_summary.update_iter();
            }
        }
    }

    fn first_non_trivial_breakdown(
        &self,
        metadata: &PersistedDiscreteBreakdownMetadata,
        pdfs: &[F<f64>],
    ) -> Option<DiscreteBreakdown> {
        (self.bins.len() > 1).then(|| DiscreteBreakdown {
            axis_label: metadata.axis_label.clone(),
            fixed_coordinates: metadata.fixed_coordinates.clone(),
            entries: self
                .bins
                .iter()
                .enumerate()
                .map(|(bin_index, bin)| DiscreteBreakdownEntry {
                    bin_index,
                    bin_label: metadata.bin_labels.get(bin_index).cloned(),
                    pdf: pdfs.get(bin_index).copied().unwrap_or(F(0.0)),
                    value: bin.accumulator.avg,
                    error: bin.accumulator.err,
                    chi_sq: bin.accumulator.chi_sq,
                    processed_samples: bin.accumulator.processed_samples,
                })
                .collect(),
        })
    }
}

struct IntegralResultCells {
    integrand: String,
    value: String,
    relative_error: String,
    chi_sq: String,
    delta_sigma: Option<String>,
    delta_percent: Option<String>,
    mwi: String,
}

const DEFAULT_MAX_SHARED_TABLE_WIDTH: usize = 250;

fn slot_label(slot_meta: &SlotMeta) -> String {
    format!("itg {}", slot_meta.key())
}

fn slot_key_label(slot_meta: &SlotMeta) -> String {
    slot_meta.key()
}

fn format_max_eval_coordinate(value: F<f64>) -> String {
    let formatted = format!("{:.16e}", value.0);
    let Some((mantissa, exponent)) = formatted.rsplit_once('e') else {
        return formatted;
    };
    let exponent = exponent.parse::<i32>().unwrap_or_default();
    format!("{mantissa}e{exponent:+03}")
}

fn format_max_eval_coordinates(xs: &[F<f64>]) -> String {
    if xs.len() <= 3 {
        return format!(
            "[ {} ]",
            xs.iter()
                .map(|value| format_max_eval_coordinate(*value))
                .join(" ")
        );
    }

    let rows = xs
        .chunks(3)
        .map(|chunk| {
            chunk
                .iter()
                .map(|value| format_max_eval_coordinate(*value))
                .join(" ")
        })
        .join("\n");
    format!("[\n{rows} ]")
}

pub(crate) fn discrete_axis_label(axis_labels: &[String], depth: usize) -> &str {
    axis_labels.get(depth).map(String::as_str).unwrap_or("idx")
}

fn append_max_eval_sample_parts(
    sample: &Sample<F<f64>>,
    axis_labels: &[String],
    depth: usize,
    parts: &mut Vec<String>,
) {
    match sample {
        Sample::Continuous(_, xs) => {
            parts.push(format!("xs: {}", format_max_eval_coordinates(xs)));
        }
        Sample::Discrete(_, index, Some(nested_sample)) => {
            parts.push(format!(
                "{}: {}",
                discrete_axis_label(axis_labels, depth),
                index
            ));
            append_max_eval_sample_parts(nested_sample, axis_labels, depth + 1, parts);
        }
        Sample::Discrete(_, index, None) => {
            parts.push(format!(
                "{}: {}",
                discrete_axis_label(axis_labels, depth),
                index
            ));
        }
        Sample::Uniform(_, indices, xs) => {
            for (offset, index) in indices.iter().enumerate() {
                parts.push(format!(
                    "{}: {}",
                    discrete_axis_label(axis_labels, depth + offset),
                    index
                ));
            }
            parts.push(format!("xs: {}", format_max_eval_coordinates(xs)));
        }
    }
}

fn format_max_eval_sample(
    sample: &Sample<F<f64>>,
    axis_labels: &[String],
    prefix_path: &[usize],
) -> String {
    let mut parts = prefix_path
        .iter()
        .enumerate()
        .map(|(depth, index)| format!("{}: {}", discrete_axis_label(axis_labels, depth), index))
        .collect_vec();
    append_max_eval_sample_parts(sample, axis_labels, prefix_path.len(), &mut parts);
    if parts.is_empty() {
        String::from("N/A")
    } else {
        parts.join(", ")
    }
}

fn format_iteration_points(points: usize) -> String {
    format_abbreviated_count(points)
}

fn format_total_points(points: usize) -> String {
    format_abbreviated_count(points)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
enum ComponentKind {
    Real,
    Imag,
}

impl ComponentKind {
    fn all_for_display(display: IntegrationStatusPhaseDisplay) -> Vec<Self> {
        let mut components = Vec::new();
        if display.shows_real() {
            components.push(Self::Real);
        }
        if display.shows_imag() {
            components.push(Self::Imag);
        }
        components
    }

    fn tag(self) -> &'static str {
        match self {
            Self::Real => "re",
            Self::Imag => "im",
        }
    }

    fn colorized_tag(self) -> String {
        self.tag().blue().bold().to_string()
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ContributionKind {
    All,
    Sum,
    Bin(usize),
}

impl ContributionKind {
    fn label(self, integration_state: &IntegrationState) -> String {
        match self {
            Self::All => "All".bold().green().to_string(),
            Self::Sum => "Sum".bold().green().to_string(),
            Self::Bin(bin_index) => contribution_bin_label(integration_state, bin_index),
        }
    }
}

#[derive(Clone, Debug)]
struct DiscreteLevelContext {
    path: Vec<usize>,
    pdfs: Vec<F<f64>>,
}

fn first_non_trivial_discrete_path(grid: &Grid<F<f64>>) -> Option<Vec<usize>> {
    match grid {
        Grid::Discrete(discrete_grid) => {
            if discrete_grid.bins.len() > 1 {
                Some(Vec::new())
            } else {
                discrete_grid
                    .bins
                    .first()
                    .and_then(|bin| bin.sub_grid.as_ref())
                    .and_then(first_non_trivial_discrete_path)
                    .map(|mut path| {
                        path.insert(0, 0);
                        path
                    })
            }
        }
        Grid::Continuous(_) | Grid::Uniform(_, _) => None,
    }
}

fn discrete_grid_at_path<'a>(
    grid: &'a Grid<F<f64>>,
    path: &[usize],
) -> Option<&'a DiscreteGrid<F<f64>>> {
    match (grid, path.split_first()) {
        (Grid::Discrete(discrete_grid), None) => Some(discrete_grid),
        (Grid::Discrete(discrete_grid), Some((bin_index, rest))) => discrete_grid
            .bins
            .get(*bin_index)?
            .sub_grid
            .as_ref()
            .and_then(|sub_grid| discrete_grid_at_path(sub_grid, rest)),
        _ => None,
    }
}

fn summary_at_path<'a>(
    summary: &'a DiscreteGridAccumulatorSummary,
    path: &[usize],
) -> Option<&'a DiscreteGridAccumulatorSummary> {
    match path.split_first() {
        None => Some(summary),
        Some((bin_index, rest)) => summary
            .bins
            .get(*bin_index)?
            .sub_summary
            .as_deref()
            .and_then(|sub_summary| summary_at_path(sub_summary, rest)),
    }
}

fn first_non_trivial_discrete_context(
    sampling_grid: &Grid<F<f64>>,
) -> Option<DiscreteLevelContext> {
    let path = first_non_trivial_discrete_path(sampling_grid)?;
    let discrete_grid = discrete_grid_at_path(sampling_grid, &path)?;
    Some(DiscreteLevelContext {
        path,
        pdfs: discrete_grid.bins.iter().map(|bin| bin.pdf).collect(),
    })
}

fn discrete_axis_labels(sampling: &SamplingSettings) -> Vec<&'static str> {
    match sampling {
        SamplingSettings::Default(_) | SamplingSettings::MultiChanneling(_) => Vec::new(),
        SamplingSettings::DiscreteGraphs(settings) => {
            let mut labels = vec!["graph"];
            match &settings.sampling_type {
                DiscreteGraphSamplingType::Default(_)
                | DiscreteGraphSamplingType::MultiChanneling(_)
                | DiscreteGraphSamplingType::TropicalSampling(_) => {
                    if settings.sample_orientations {
                        labels.push("orientation");
                    }
                }
                DiscreteGraphSamplingType::DiscreteMultiChanneling(_) => {
                    if settings.sample_orientations {
                        labels.push("orientation");
                    }
                    labels.push("LMB channel");
                }
            }
            labels
        }
    }
}

fn orientation_description(orientation: &EdgeVec<Orientation>) -> String {
    orientation
        .iter()
        .map(|(_, orientation)| match *orientation {
            Orientation::Default => '+',
            Orientation::Reversed => '-',
            Orientation::Undirected => '0',
        })
        .collect()
}

fn graph_group_description<I>(graph_names: I) -> String
where
    I: IntoIterator<Item = String>,
{
    let graph_names = graph_names.into_iter().collect_vec();
    if graph_names.len() <= 1 {
        return graph_names.into_iter().next().unwrap_or_default();
    }

    format!("[{}]", graph_names.join(","))
}

fn lmb_channel_description(lmb: &LoopMomentumBasis) -> String {
    format!(
        "({})",
        lmb.loop_edges
            .iter()
            .map(|edge_id| edge_id.0.to_string())
            .join(",")
    )
}

fn first_non_trivial_discrete_bin_descriptions_for_process_integrand(
    integrand: &ProcessIntegrand,
    path: &[usize],
    axis_label: &str,
) -> Option<Vec<String>> {
    match (integrand, axis_label) {
        (ProcessIntegrand::Amplitude(integrand), "graph") => {
            Some(
                integrand
                    .data
                    .graph_group_structure
                    .iter()
                    .map(|group| {
                        graph_group_description(group.into_iter().map(|graph_id| {
                            integrand.data.graph_terms[graph_id].graph.name.clone()
                        }))
                    })
                    .collect(),
            )
        }
        (ProcessIntegrand::CrossSection(integrand), "graph") => {
            Some(
                integrand
                    .data
                    .graph_group_structure
                    .iter()
                    .map(|group| {
                        graph_group_description(group.into_iter().map(|graph_id| {
                            integrand.data.graph_terms[graph_id].graph.name.clone()
                        }))
                    })
                    .collect(),
            )
        }
        (ProcessIntegrand::Amplitude(integrand), "orientation") => {
            let group_id = GroupId(*path.first()?);
            let group = integrand.data.graph_group_structure.get(group_id)?;
            let master = group.master();
            Some(
                integrand.data.graph_terms[master]
                    .orientations
                    .iter()
                    .map(orientation_description)
                    .collect(),
            )
        }
        (ProcessIntegrand::CrossSection(integrand), "orientation") => {
            let group_id = GroupId(*path.first()?);
            let group = integrand.data.graph_group_structure.get(group_id)?;
            let master = group.master();
            Some(
                integrand.data.graph_terms[master]
                    .orientations
                    .iter()
                    .map(orientation_description)
                    .collect(),
            )
        }
        (ProcessIntegrand::Amplitude(integrand), "LMB channel") => {
            let group_id = GroupId(*path.first()?);
            let group = integrand.data.graph_group_structure.get(group_id)?;
            let master = group.master();
            let graph_term = &integrand.data.graph_terms[master];
            Some(
                graph_term
                    .multi_channeling_setup
                    .channels
                    .iter()
                    .map(|&channel_lmb| {
                        lmb_channel_description(
                            &graph_term.multi_channeling_setup.all_bases[channel_lmb],
                        )
                    })
                    .collect(),
            )
        }
        (ProcessIntegrand::CrossSection(integrand), "LMB channel") => {
            let group_id = GroupId(*path.first()?);
            let group = integrand.data.graph_group_structure.get(group_id)?;
            let master = group.master();
            let graph_term = &integrand.data.graph_terms[master];
            Some(
                graph_term
                    .multi_channeling_setup
                    .channels
                    .iter()
                    .map(|&channel_lmb| {
                        lmb_channel_description(
                            &graph_term.multi_channeling_setup.all_bases[channel_lmb],
                        )
                    })
                    .collect(),
            )
        }
        _ => None,
    }
}

fn first_non_trivial_discrete_bin_descriptions_for_integrand(
    integrand: &Integrand,
    path: &[usize],
    axis_label: &str,
) -> Option<Vec<String>> {
    match integrand {
        Integrand::ProcessIntegrand(process_integrand) => {
            first_non_trivial_discrete_bin_descriptions_for_process_integrand(
                process_integrand,
                path,
                axis_label,
            )
        }
        _ => None,
    }
}

fn discrete_coordinate_label(
    integrand: &Integrand,
    path_prefix: &[usize],
    axis_label: &str,
    bin_index: usize,
) -> Option<String> {
    first_non_trivial_discrete_bin_descriptions_for_integrand(integrand, path_prefix, axis_label)?
        .get(bin_index)
        .cloned()
}

fn build_persisted_discrete_breakdown_metadata(
    integrand: &Integrand,
    path: &[usize],
    discrete_axis_labels: &[String],
) -> Option<PersistedDiscreteBreakdownMetadata> {
    let axis_label = discrete_axis_labels.get(path.len())?.clone();
    let bin_labels =
        first_non_trivial_discrete_bin_descriptions_for_integrand(integrand, path, &axis_label)?;

    let fixed_coordinates = path
        .iter()
        .enumerate()
        .map(|(axis_index, &bin_index)| {
            let fixed_axis_label = discrete_axis_labels.get(axis_index)?.clone();
            Some(DiscreteCoordinate {
                axis_label: fixed_axis_label.clone(),
                bin_index,
                bin_label: discrete_coordinate_label(
                    integrand,
                    &path[..axis_index],
                    &fixed_axis_label,
                    bin_index,
                ),
            })
        })
        .collect::<Option<Vec<_>>>()?;

    Some(PersistedDiscreteBreakdownMetadata {
        axis_label,
        fixed_coordinates,
        bin_labels,
    })
}

fn monitored_discrete_layout(
    grid: &Grid<F<f64>>,
    discrete_axis_labels: &[String],
) -> Option<(Vec<usize>, String, usize)> {
    let path = first_non_trivial_discrete_path(grid)?;
    let axis_label = discrete_axis_labels.get(path.len())?.clone();
    let discrete_grid = discrete_grid_at_path(grid, &path)?;
    Some((path, axis_label, discrete_grid.bins.len()))
}

fn coalesce_first_non_trivial_discrete_bin_descriptions(
    axis_label: &str,
    slot_descriptions: &[(String, Vec<String>)],
) -> (Option<Vec<String>>, Option<String>) {
    let Some((reference_slot, reference_descriptions)) = slot_descriptions.first() else {
        return (None, None);
    };

    if let Some((slot_key, _)) = slot_descriptions
        .iter()
        .skip(1)
        .find(|(_, descriptions)| descriptions != reference_descriptions)
    {
        return (
            None,
            Some(format!(
                "Selected integrands do not share the same semantic labels for the monitored discrete dimension '{axis_label}' ({} vs {}). Falling back to raw bin indices in integration reporting.",
                reference_slot.blue(),
                slot_key.blue()
            )),
        );
    }

    (Some(reference_descriptions.clone()), None)
}

fn resolve_first_non_trivial_discrete_bin_descriptions(
    slots: &[IntegrationSlot],
    monitored_path: &[usize],
    axis_label: &str,
) -> (Option<Vec<String>>, Option<String>) {
    let Some(slot_descriptions) = slots
        .iter()
        .map(|slot| {
            first_non_trivial_discrete_bin_descriptions_for_integrand(
                &slot.integrand,
                monitored_path,
                axis_label,
            )
            .map(|descriptions| (slot.meta.key(), descriptions))
        })
        .collect::<Option<Vec<_>>>()
    else {
        return (None, None);
    };

    coalesce_first_non_trivial_discrete_bin_descriptions(axis_label, &slot_descriptions)
}

fn resolve_monitored_discrete_setup(
    sampling_correlation_mode: SamplingCorrelationMode,
    slots: &[IntegrationSlot],
    sampling_states: &[SamplingSlotState],
) -> MonitoredDiscreteSetup {
    let Some(reference_state) = sampling_states.first() else {
        return MonitoredDiscreteSetup::default();
    };

    let Some((reference_path, reference_axis_label, reference_bin_count)) =
        monitored_discrete_layout(&reference_state.grid, &reference_state.discrete_axis_labels)
    else {
        return MonitoredDiscreteSetup::default();
    };

    if sampling_correlation_mode == SamplingCorrelationMode::Uncorrelated {
        for (slot_index, sampling_state) in sampling_states.iter().enumerate().skip(1) {
            let Some((path, axis_label, bin_count)) = monitored_discrete_layout(
                &sampling_state.grid,
                &sampling_state.discrete_axis_labels,
            ) else {
                return MonitoredDiscreteSetup {
                    warning: Some(
                        "Selected integrands do not all expose a compatible first non-trivial monitored discrete dimension. Shared discrete-bin monitoring tables will be disabled."
                            .to_string(),
                    ),
                    ..MonitoredDiscreteSetup::default()
                };
            };

            if path != reference_path
                || axis_label != reference_axis_label
                || bin_count != reference_bin_count
            {
                return MonitoredDiscreteSetup {
                    warning: Some(format!(
                        "Selected integrands do not share a compatible first non-trivial monitored discrete layout (mismatch at {}). Shared discrete-bin monitoring tables will be disabled.",
                        slots[slot_index].meta.key().blue()
                    )),
                    ..MonitoredDiscreteSetup::default()
                };
            }
        }
    }

    let (descriptions, label_warning) = resolve_first_non_trivial_discrete_bin_descriptions(
        slots,
        &reference_path,
        &reference_axis_label,
    );

    MonitoredDiscreteSetup {
        path: Some(reference_path),
        axis_label: Some(reference_axis_label),
        bin_descriptions: descriptions,
        warning: label_warning,
    }
}

fn render_orientation_description(description: &str) -> String {
    description
        .chars()
        .map(|sign| match sign {
            '+' => "+".green().bold().to_string(),
            '-' => "-".red().bold().to_string(),
            '0' => "0".dimmed().to_string(),
            other => other.to_string(),
        })
        .join("")
}

fn render_graph_description(description: &str) -> String {
    if let Some(inner) = description
        .strip_prefix('[')
        .and_then(|trimmed| trimmed.strip_suffix(']'))
    {
        let names = inner
            .split(',')
            .filter(|name| !name.is_empty())
            .collect_vec();
        if names.is_empty() {
            return description.bold().green().to_string();
        }

        let rendered = names
            .into_iter()
            .enumerate()
            .map(|(index, name)| {
                if index == 0 {
                    name.bold().green().to_string()
                } else {
                    name.bold().blue().to_string()
                }
            })
            .join(",");
        return format!("[{rendered}]");
    }

    description.bold().green().to_string()
}

fn render_bin_description(axis_label: &str, description: &str) -> String {
    match axis_label {
        "orientation" => render_orientation_description(description),
        "graph" => render_graph_description(description),
        _ => description.bold().green().to_string(),
    }
}

fn contribution_bin_label(integration_state: &IntegrationState, bin_index: usize) -> String {
    if let (Some(axis_label), Some(descriptions)) = (
        integration_state
            .first_non_trivial_discrete_label
            .as_deref(),
        integration_state
            .first_non_trivial_discrete_bin_descriptions
            .as_ref(),
    ) && let Some(description) = descriptions.get(bin_index)
    {
        return format!(
            "#{bin_index}: {}",
            render_bin_description(axis_label, description)
        );
    }

    format!("#{bin_index}")
}

fn contribution_header_label(
    integration_state: &IntegrationState,
    discrete_monitoring_enabled: bool,
) -> String {
    if discrete_monitoring_enabled
        && let Some(label) = integration_state
            .first_non_trivial_discrete_label
            .as_deref()
    {
        return format!("Contribution (idx={label})")
            .bold()
            .blue()
            .to_string();
    }

    "Contribution".bold().blue().to_string()
}

fn format_percentage_sig(value: f64, significant_digits: usize) -> String {
    if !value.is_finite() {
        return "None".red().to_string();
    }

    if value == 0.0 {
        let decimals = significant_digits.saturating_sub(1);
        return format!("{:.*}%", decimals, 0.0);
    }

    let abs_value = value.abs();
    let exponent = abs_value.log10().floor() as i32;
    if exponent < -2 || exponent >= significant_digits as i32 {
        return format!("{:.*e}%", significant_digits.saturating_sub(1), value);
    }

    let decimals = (significant_digits as i32 - exponent - 1).max(0) as usize;
    format!("{value:.decimals$}%")
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum UncertaintyNotation {
    Dynamic,
    Scientific,
}

fn format_significant_percentage(
    value: f64,
    significant_digits: usize,
    scientific_threshold: Option<(f64, f64)>,
) -> String {
    if !value.is_finite() {
        return "None".red().to_string();
    }

    if value == 0.0 {
        return format!("{:.*}%", significant_digits.saturating_sub(1), 0.0);
    }

    let abs_value = value.abs();
    if scientific_threshold.is_some_and(|(upper, lower)| abs_value >= upper || abs_value < lower) {
        return format!("{:.*e}%", significant_digits.saturating_sub(1), value);
    }

    let exponent = abs_value.log10().floor() as i32;
    let decimals = (significant_digits as i32 - exponent - 1).max(0) as usize;
    format!("{value:.decimals$}%")
}

fn should_use_scientific_uncertainty_notation(avg: F<f64>, err: F<f64>) -> bool {
    let avg_abs = avg.abs().0;
    let err_abs = err.abs().0;
    avg_abs >= 1e6
        || (avg_abs != 0.0 && avg_abs < 1e-5)
        || (avg.is_zero() && !(1e-4..1e5).contains(&err_abs))
}

fn format_uncertainty_with_notation(
    avg: F<f64>,
    err: F<f64>,
    notation: UncertaintyNotation,
) -> String {
    if !matches!(notation, UncertaintyNotation::Scientific)
        && !should_use_scientific_uncertainty_notation(avg, err)
    {
        return utils::format_uncertainty(avg, err);
    }

    if !avg.0.is_finite() || !err.0.is_finite() {
        return utils::format_uncertainty(avg, err);
    }

    let exponent = if avg.is_non_zero() {
        avg.abs().0.log10().floor() as i32
    } else if err.is_non_zero() {
        err.abs().0.log10().floor() as i32
    } else {
        0
    };
    let scale = 10_f64.powi(exponent);
    let mantissa = utils::format_uncertainty(F(avg.0 / scale), F(err.0 / scale));
    format!("{mantissa}e{exponent}")
}

fn format_signed_uncertainty(avg: F<f64>, err: F<f64>, notation: UncertaintyNotation) -> String {
    let formatted = format_uncertainty_with_notation(avg, err, notation);
    if avg.0.is_sign_negative() {
        formatted
    } else {
        format!("+{formatted}")
    }
}

fn format_abbreviated_count(value: usize) -> String {
    if value < 1_000 {
        return value.to_string();
    }

    let value = value as f64;
    if value < 1_000_000.0 {
        return format!("{:.2}K", value / 1_000.0);
    }
    if value < 1_000_000_000.0 {
        return format!("{:.2}M", value / 1_000_000.0);
    }
    if value < 1_000_000_000_000.0 {
        return format!("{:.2}B", value / 1_000_000_000.0);
    }

    format!("{:.3}T", value / 1_000_000_000_000.0)
}

fn build_integral_result_cells(
    itg: &StatisticsAccumulator<F<f64>>,
    slot_meta: &SlotMeta,
    i_iter: usize,
    tag: &str,
    trgt: Option<F<f64>>,
) -> IntegralResultCells {
    let relative_error = format_relative_error_cell(itg);
    let chi_sq = format_chi_sq_cell(itg, i_iter);
    let (delta_sigma, delta_percent) = format_delta_cells(itg, trgt);
    let mwi = format_mwi_cell(itg);

    IntegralResultCells {
        integrand: format!(
            "{} {}:",
            slot_label(slot_meta),
            format!("{:-2}", tag).blue().bold(),
        ),
        value: format_signed_uncertainty(itg.avg, itg.err, UncertaintyNotation::Scientific)
            .blue()
            .bold()
            .to_string(),
        relative_error,
        chi_sq,
        delta_sigma,
        delta_percent,
        mwi,
    }
}

fn format_relative_error_cell(itg: &StatisticsAccumulator<F<f64>>) -> String {
    format_relative_error_from_estimate(itg.avg, itg.err)
}

fn format_relative_error_from_estimate(avg: F<f64>, err: F<f64>) -> String {
    if avg.is_zero() {
        return String::new();
    }

    let formatted =
        format_significant_percentage((err / avg).abs().0 * 100.0, 3, Some((1.0e4, 1.0e-4)));
    if (err / avg).abs().0 > 0.01 {
        formatted.red().to_string()
    } else {
        formatted.green().to_string()
    }
}

fn format_chi_sq_cell(itg: &StatisticsAccumulator<F<f64>>, i_iter: usize) -> String {
    let chi_sq = format!("{:.3}", itg.chi_sq.0 / (i_iter as f64));
    if itg.chi_sq / F::<f64>::new_from_usize(i_iter) > F(5.) {
        chi_sq.red().to_string()
    } else {
        chi_sq
    }
}

fn format_delta_cells(
    itg: &StatisticsAccumulator<F<f64>>,
    trgt: Option<F<f64>>,
) -> (Option<String>, Option<String>) {
    format_delta_cells_from_estimate(itg.avg, itg.err, trgt)
}

fn format_delta_cells_from_estimate(
    avg: F<f64>,
    err: F<f64>,
    trgt: Option<F<f64>>,
) -> (Option<String>, Option<String>) {
    let Some(t) = trgt else {
        return (None, None);
    };

    let delta_in_sigmas = if err.is_zero() {
        0.0
    } else {
        (t - avg).abs().0 / err.0
    };
    let delta_in_percent = if t.abs().is_non_zero() {
        (t - avg).abs().0 / t.abs().0 * 100.
    } else {
        0.
    };
    let is_outside_target =
        delta_in_sigmas > 5. || (t.abs().is_non_zero() && ((t - avg).abs() / t.abs()).0 > 0.01);
    let sigma_text = format!("Δ = {:.3}σ", delta_in_sigmas);
    let percent_text = format!("Δ = {:.3}%", delta_in_percent);
    if is_outside_target {
        (
            Some(sigma_text.red().to_string()),
            Some(percent_text.red().to_string()),
        )
    } else {
        (
            Some(sigma_text.green().to_string()),
            Some(percent_text.green().to_string()),
        )
    }
}

fn format_mwi_cell(itg: &StatisticsAccumulator<F<f64>>) -> String {
    let mwi_value = max_weight_impact(itg);
    let formatted = format!("{:.4e}", mwi_value.0);
    if mwi_value > F(1.) {
        formatted.red().to_string()
    } else {
        formatted
    }
}

fn max_weight_impact(itg: &StatisticsAccumulator<F<f64>>) -> F<f64> {
    if itg.avg.abs().0 == 0. || itg.processed_samples == 0 {
        return F(0.0);
    }

    itg.max_eval_negative.abs().max(itg.max_eval_positive.abs())
        / (itg.avg.abs() * F::<f64>::new_from_usize(itg.processed_samples))
}

fn build_table_result_summary(
    slot_meta: &SlotMeta,
    accumulator: &ComplexAccumulator,
    iter: usize,
    target: Option<Complex<F<f64>>>,
) -> Vec<IntegrationTableComponentResult> {
    [
        ("re", &accumulator.re, target.as_ref().map(|value| value.re)),
        ("im", &accumulator.im, target.as_ref().map(|value| value.im)),
    ]
    .into_iter()
    .map(|(component, accumulator, target_component)| {
        let cells =
            build_integral_result_cells(accumulator, slot_meta, iter, component, target_component);
        IntegrationTableComponentResult {
            component: component.to_string(),
            value: accumulator.avg,
            error: accumulator.err,
            relative_error_percent: accumulator
                .avg
                .is_non_zero()
                .then(|| (accumulator.err / accumulator.avg).abs().0 * 100.0),
            chi_sq_per_dof: if iter > 0 {
                accumulator.chi_sq.0 / (iter as f64)
            } else {
                0.0
            },
            target_delta_sigma: cells.delta_sigma.as_ref().map(|_| {
                if accumulator.err.is_zero() {
                    0.0
                } else if let Some(target_value) = target_component {
                    (target_value - accumulator.avg).abs().0 / accumulator.err.0
                } else {
                    0.0
                }
            }),
            target_delta_percent: target_component.map(|target_value| {
                if target_value.is_zero() {
                    0.0
                } else {
                    (target_value - accumulator.avg).abs().0 / target_value.abs().0 * 100.0
                }
            }),
            max_weight_impact: max_weight_impact(accumulator).0,
        }
    })
    .collect()
}

fn build_max_weight_info_summary(
    discrete_axis_labels: &[String],
    accumulator: &ComplexAccumulator,
) -> Vec<MaxWeightInfoEntry> {
    [
        (
            "re",
            "+",
            accumulator.re.max_eval_positive,
            accumulator.re.max_eval_positive_xs.as_ref(),
        ),
        (
            "re",
            "-",
            accumulator.re.max_eval_negative,
            accumulator.re.max_eval_negative_xs.as_ref(),
        ),
        (
            "im",
            "+",
            accumulator.im.max_eval_positive,
            accumulator.im.max_eval_positive_xs.as_ref(),
        ),
        (
            "im",
            "-",
            accumulator.im.max_eval_negative,
            accumulator.im.max_eval_negative_xs.as_ref(),
        ),
    ]
    .into_iter()
    .filter_map(|(component, sign, max_eval, sample)| {
        if max_eval.is_zero() {
            return None;
        }

        Some(MaxWeightInfoEntry {
            component: component.to_string(),
            sign: sign.to_string(),
            max_eval,
            coordinates: sample
                .map(|sample| format_max_eval_sample(sample, discrete_axis_labels, &[])),
        })
    })
    .collect()
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct LiveIterationProgress {
    completed_points: usize,
    target_points: usize,
}

fn build_iteration_status_header_left(
    elapsed_time: Duration,
    iter: usize,
    live_progress: Option<LiveIterationProgress>,
) -> String {
    let iteration_label = if live_progress.is_some() {
        format!("Iteration #{:-4} ( running   )", iter)
            .bold()
            .green()
    } else {
        format!("Iteration #{:-4} ( completed )", iter)
            .bold()
            .green()
    };

    format!(
        "[ {} ] {}",
        format!(
            "{:^7}",
            utils::format_wdhms(elapsed_time.as_secs() as usize)
        )
        .bold(),
        iteration_label,
    )
}

fn build_iteration_status_header_middle(
    cur_points: usize,
    total_points: usize,
    live_progress: Option<LiveIterationProgress>,
) -> String {
    let per_iteration = if let Some(progress) = live_progress {
        let percentage = if progress.target_points == 0 {
            String::from("0.0%")
        } else {
            format!(
                "{:.1}%",
                (progress.completed_points as f64) / (progress.target_points as f64) * 100.0
            )
            .green()
            .to_string()
        };
        format!(
            "Iteration progress {}/{} {}",
            format_iteration_points(progress.completed_points),
            format_iteration_points(progress.target_points),
            percentage
        )
    } else if cur_points > 0 {
        format!(
            "# samples per iteration = {}",
            format_iteration_points(cur_points)
        )
        .blue()
        .bold()
        .to_string()
    } else {
        String::new()
    };

    let total = format!("# samples total = {}", format_total_points(total_points))
        .bold()
        .green()
        .to_string();
    if per_iteration.is_empty() {
        total
    } else {
        format!("{per_iteration} {total}")
    }
}

fn build_iteration_status_header_tail(
    cores: usize,
    elapsed_time: Duration,
    n_samples_evaluated: usize,
) -> String {
    let average_sample_time = if n_samples_evaluated == 0 {
        "N/A".red().to_string()
    } else {
        utils::format_evaluation_time_from_f64(
            elapsed_time.as_secs_f64() / (n_samples_evaluated as f64) * (cores as f64),
        )
        .bold()
        .blue()
        .to_string()
    };

    format!("{average_sample_time} /sample/core")
}

fn max_weight_row_descriptors(
    phase_display: IntegrationStatusPhaseDisplay,
) -> Vec<(ComponentKind, &'static str, bool)> {
    let mut rows = Vec::new();
    if phase_display.shows_real() {
        rows.push((ComponentKind::Real, "+", true));
        rows.push((ComponentKind::Real, "-", false));
    }
    if phase_display.shows_imag() {
        rows.push((ComponentKind::Imag, "+", true));
        rows.push((ComponentKind::Imag, "-", false));
    }
    rows
}

type MaxEvalEntry<'a> = Option<(F<f64>, Option<&'a Sample<F<f64>>>)>;

fn max_eval_entry(accumulator: &StatisticsAccumulator<F<f64>>, positive: bool) -> MaxEvalEntry<'_> {
    let (value, sample) = if positive {
        (
            accumulator.max_eval_positive,
            accumulator.max_eval_positive_xs.as_ref(),
        )
    } else {
        (
            accumulator.max_eval_negative,
            accumulator.max_eval_negative_xs.as_ref(),
        )
    };

    if value.is_zero() {
        None
    } else {
        Some((value, sample))
    }
}

/// struct to keep track of state, used in the havana_integrate function
/// the idea is to save this to disk after each iteration, so that the integration can be resumed
#[derive(Clone, Serialize, Deserialize, Encode, Decode)]
pub struct IntegrationState {
    pub num_points: usize,
    #[bincode(with_serde)]
    pub all_integrals: Vec<ComplexAccumulator>,
    #[bincode(with_serde)]
    slot_re_summaries: Vec<Option<DiscreteGridAccumulatorSummary>>,
    #[bincode(with_serde)]
    slot_im_summaries: Vec<Option<DiscreteGridAccumulatorSummary>>,
    pub stats: StatisticsCounter,
    pub slot_metas: Vec<SlotMeta>,
    pub sampling_correlation_mode: SamplingCorrelationMode,
    sampling_states: Vec<SamplingSlotState>,
    monitored_discrete_path: Option<Vec<usize>>,
    pub first_non_trivial_discrete_label: Option<String>,
    pub first_non_trivial_discrete_bin_descriptions: Option<Vec<String>>,
    slot_first_non_trivial_discrete_breakdown_metadata:
        Vec<Option<PersistedDiscreteBreakdownMetadata>>,
    pub iter: usize,
    pub elapsed_seconds: f64,
    pub n_cores: usize,
}

impl IntegrationState {
    fn new_from_settings(
        sampling_correlation_mode: SamplingCorrelationMode,
        sampling_states: Vec<SamplingSlotState>,
        slot_metas: Vec<SlotMeta>,
        monitored_discrete_path: Option<Vec<usize>>,
        first_non_trivial_discrete_label: Option<String>,
        first_non_trivial_discrete_bin_descriptions: Option<Vec<String>>,
        slot_first_non_trivial_discrete_breakdown_metadata: Vec<
            Option<PersistedDiscreteBreakdownMetadata>,
        >,
    ) -> Self {
        let num_points = 0;
        let iter = 0;
        let all_integrals = vec![ComplexAccumulator::new(); slot_metas.len()];
        let slot_re_summaries = (0..slot_metas.len())
            .map(|slot_index| {
                DiscreteGridAccumulatorSummary::from_grid(
                    &sampling_states[sampling_correlation_mode.state_index(slot_index)].grid,
                )
            })
            .collect();
        let slot_im_summaries = (0..slot_metas.len())
            .map(|slot_index| {
                DiscreteGridAccumulatorSummary::from_grid(
                    &sampling_states[sampling_correlation_mode.state_index(slot_index)].grid,
                )
            })
            .collect();
        let stats = StatisticsCounter::new_empty();

        Self {
            num_points,
            all_integrals,
            slot_re_summaries,
            slot_im_summaries,
            stats,
            slot_metas,
            sampling_correlation_mode,
            sampling_states,
            monitored_discrete_path,
            first_non_trivial_discrete_label,
            first_non_trivial_discrete_bin_descriptions,
            slot_first_non_trivial_discrete_breakdown_metadata,
            iter,
            elapsed_seconds: 0.0,
            n_cores: 1,
        }
    }

    fn update_iter(&mut self, use_weighted_average: bool) {
        self.all_integrals
            .iter_mut()
            .for_each(|acc| acc.update_iter(use_weighted_average));
    }

    fn sampling_state_for_slot(&self, slot_index: usize) -> &SamplingSlotState {
        &self.sampling_states[self.sampling_correlation_mode.state_index(slot_index)]
    }

    fn sampling_state_for_slot_mut(&mut self, slot_index: usize) -> &mut SamplingSlotState {
        let sampling_state_index = self.sampling_correlation_mode.state_index(slot_index);
        &mut self.sampling_states[sampling_state_index]
    }

    fn monitored_discrete_context_for_slot(
        &self,
        slot_index: usize,
    ) -> Option<DiscreteLevelContext> {
        let path = self.monitored_discrete_path.clone()?;
        let discrete_grid =
            discrete_grid_at_path(&self.sampling_state_for_slot(slot_index).grid, &path)?;
        Some(DiscreteLevelContext {
            path,
            pdfs: discrete_grid.bins.iter().map(|bin| bin.pdf).collect(),
        })
    }

    fn monitored_discrete_context(&self) -> Option<DiscreteLevelContext> {
        self.monitored_discrete_context_for_slot(0)
    }
}

struct CoreIterationState {
    slot_integrands: Vec<Integrand>,
    stats: StatisticsCounter,
    integrals: Vec<ComplexAccumulator>,
    sampling_correlation_mode: SamplingCorrelationMode,
    sampling_states: Vec<CoreSamplingSlotState>,
    slot_re_grids: Vec<Grid<F<f64>>>,
    slot_im_grids: Vec<Grid<F<f64>>>,
    remaining_points: usize,
    completed_points: usize,
}

struct CoreSamplingSlotState {
    sampling_grid: Grid<F<f64>>,
    rng: MonteCarloRng,
}

impl CoreIterationState {
    fn new(
        slot_integrands: Vec<Integrand>,
        sampling_correlation_mode: SamplingCorrelationMode,
        sampling_grid_templates: &[Grid<F<f64>>],
        seed: u64,
        sample_skip: usize,
        remaining_points: usize,
    ) -> Self {
        let n_slots = slot_integrands.len();
        let sampling_states = (0..sampling_correlation_mode.state_count(n_slots))
            .map(|sampling_state_index| {
                let slot_index = match sampling_correlation_mode {
                    SamplingCorrelationMode::Correlated => 0,
                    SamplingCorrelationMode::Uncorrelated => sampling_state_index,
                };
                let mut rng = MonteCarloRng::new(slot_seed(seed, slot_index), 0);
                let mut sampling_grid =
                    sampling_grid_templates[sampling_state_index].clone_without_samples();

                for _ in 0..sample_skip {
                    let mut sample = Sample::new();
                    sampling_grid.sample(&mut rng, &mut sample);
                }

                CoreSamplingSlotState { sampling_grid, rng }
            })
            .collect_vec();

        Self {
            slot_integrands,
            stats: StatisticsCounter::new_empty(),
            integrals: vec![ComplexAccumulator::new(); n_slots],
            sampling_correlation_mode,
            sampling_states,
            slot_re_grids: (0..n_slots)
                .map(|slot_index| {
                    sampling_grid_templates[sampling_correlation_mode.state_index(slot_index)]
                        .clone_without_samples()
                })
                .collect(),
            slot_im_grids: (0..n_slots)
                .map(|slot_index| {
                    sampling_grid_templates[sampling_correlation_mode.state_index(slot_index)]
                        .clone_without_samples()
                })
                .collect(),
            remaining_points,
            completed_points: 0,
        }
    }

    fn evaluate_chunk(
        &mut self,
        slot_settings: &[&RuntimeSettings],
        slot_models: &[&Model],
        iter: usize,
        current_max_evals: &[Complex<F<f64>>],
        chunk_size: usize,
    ) -> Result<usize> {
        let n_points = chunk_size.min(self.remaining_points);
        if n_points == 0 {
            return Ok(0);
        }

        let chunk_start = Instant::now();
        let mut batch_evaluation_time = Duration::ZERO;
        let mut batch_stats = StatisticsCounter::new_empty();
        let mut processed_points = 0;
        let mut total_sample_evaluations = 0;

        match self.sampling_correlation_mode {
            SamplingCorrelationMode::Correlated => {
                let mut samples = Vec::with_capacity(n_points);
                for _ in 0..n_points {
                    if is_interrupted() {
                        break;
                    }

                    let mut sample = Sample::new();
                    let sampling_state = &mut self.sampling_states[0];
                    sampling_state
                        .sampling_grid
                        .sample(&mut sampling_state.rng, &mut sample);
                    processed_points += 1;
                    samples.push(sample);
                }

                if processed_points == 0 {
                    return Ok(0);
                }

                let mut slot_results = Vec::with_capacity(self.slot_integrands.len());

                for (slot_index, integrand) in self.slot_integrands.iter_mut().enumerate() {
                    let evaluation_start = Instant::now();
                    let raw_batch = integrand.evaluate_samples_raw(
                        &samples,
                        slot_models[slot_index],
                        iter,
                        false,
                        true,
                        current_max_evals[slot_index],
                    )?;
                    if raw_batch.samples.len() < samples.len() {
                        return Ok(0);
                    }
                    batch_evaluation_time += evaluation_start.elapsed();
                    total_sample_evaluations += raw_batch.samples.len();
                    batch_stats = batch_stats.merged(&raw_batch.statistics);
                    slot_results.push(raw_batch.samples);
                }

                for (sample_index, sample) in samples.iter().enumerate() {
                    for (slot_index, (((core_accumulator, re_grid), im_grid), results)) in self
                        .integrals
                        .iter_mut()
                        .zip(self.slot_re_grids.iter_mut())
                        .zip(self.slot_im_grids.iter_mut())
                        .zip(slot_results.iter())
                        .enumerate()
                    {
                        let result = &results[sample_index];
                        let jacobian = result.parameterization_jacobian.unwrap_or(F(1.0));
                        let effective_integrand_result =
                            result.integrand_result * Complex::new_re(jacobian);

                        core_accumulator.add_sample(
                            effective_integrand_result,
                            sample.get_weight(),
                            Some(sample),
                        );
                        re_grid
                            .add_training_sample(sample, effective_integrand_result.re)
                            .map_err(Report::msg)?;
                        im_grid
                            .add_training_sample(sample, effective_integrand_result.im)
                            .map_err(Report::msg)?;

                        if slot_index == 0 {
                            let training_eval = match slot_settings[0].integrator.integrated_phase {
                                IntegratedPhase::Real => effective_integrand_result.re,
                                IntegratedPhase::Imag => effective_integrand_result.im,
                                IntegratedPhase::Both => unimplemented!(),
                            };

                            self.sampling_states[0]
                                .sampling_grid
                                .add_training_sample(sample, training_eval)
                                .map_err(Report::msg)?;
                        }
                    }
                }
            }
            SamplingCorrelationMode::Uncorrelated => {
                for slot_index in 0..self.slot_integrands.len() {
                    let mut samples = Vec::with_capacity(n_points);
                    for _ in 0..n_points {
                        if is_interrupted() {
                            break;
                        }

                        let mut sample = Sample::new();
                        let sampling_state = &mut self.sampling_states[slot_index];
                        sampling_state
                            .sampling_grid
                            .sample(&mut sampling_state.rng, &mut sample);
                        samples.push(sample);
                    }

                    if slot_index == 0 {
                        processed_points = samples.len();
                    } else if samples.len() != processed_points {
                        return Err(Report::msg(
                            "uncorrelated slot batch sizes diverged unexpectedly",
                        ));
                    }

                    if samples.is_empty() {
                        continue;
                    }

                    let evaluation_start = Instant::now();
                    let raw_batch = self.slot_integrands[slot_index].evaluate_samples_raw(
                        &samples,
                        slot_models[slot_index],
                        iter,
                        false,
                        true,
                        current_max_evals[slot_index],
                    )?;
                    if raw_batch.samples.len() < samples.len() {
                        return Ok(0);
                    }
                    batch_evaluation_time += evaluation_start.elapsed();
                    total_sample_evaluations += raw_batch.samples.len();
                    batch_stats = batch_stats.merged(&raw_batch.statistics);

                    for (sample, result) in samples.iter().zip(raw_batch.samples.iter()) {
                        let jacobian = result.parameterization_jacobian.unwrap_or(F(1.0));
                        let effective_integrand_result =
                            result.integrand_result * Complex::new_re(jacobian);

                        self.integrals[slot_index].add_sample(
                            effective_integrand_result,
                            sample.get_weight(),
                            Some(sample),
                        );
                        self.slot_re_grids[slot_index]
                            .add_training_sample(sample, effective_integrand_result.re)
                            .map_err(Report::msg)?;
                        self.slot_im_grids[slot_index]
                            .add_training_sample(sample, effective_integrand_result.im)
                            .map_err(Report::msg)?;

                        let training_eval =
                            match slot_settings[slot_index].integrator.integrated_phase {
                                IntegratedPhase::Real => effective_integrand_result.re,
                                IntegratedPhase::Imag => effective_integrand_result.im,
                                IntegratedPhase::Both => unimplemented!(),
                            };
                        self.sampling_states[slot_index]
                            .sampling_grid
                            .add_training_sample(sample, training_eval)
                            .map_err(Report::msg)?;
                    }
                }

                if processed_points == 0 {
                    return Ok(0);
                }
            }
        }

        self.remaining_points -= processed_points;
        self.completed_points += processed_points;
        self.stats = self.stats.merged(&batch_stats);
        self.stats.add_integrator_overhead(
            chunk_start.elapsed().saturating_sub(batch_evaluation_time),
            total_sample_evaluations,
        );

        Ok(processed_points)
    }
}

fn initial_batch_size(batching: IterationBatchingSettings) -> usize {
    batching.batch_size.unwrap_or(100).max(1)
}

fn next_batch_size(
    batching: IterationBatchingSettings,
    current_batch_size: usize,
    round_elapsed: Duration,
) -> usize {
    if let Some(batch_size) = batching.batch_size {
        return batch_size.max(1);
    }

    let round_seconds = round_elapsed.as_secs_f64();
    if round_seconds <= 0.0 {
        return current_batch_size.max(1);
    }

    let estimated = ((current_batch_size as f64) * batching.batch_timing_seconds / round_seconds)
        .round() as usize;
    estimated.max(1)
}

fn total_completed_points(core_states: &[CoreIterationState]) -> usize {
    core_states.iter().map(|state| state.completed_points).sum()
}

fn total_remaining_points(core_states: &[CoreIterationState]) -> usize {
    core_states.iter().map(|state| state.remaining_points).sum()
}

fn slot_seed(global_seed: u64, slot_index: usize) -> u64 {
    let slot = (slot_index as u64).wrapping_add(1);
    let mut mixed = global_seed.wrapping_add(0x9e37_79b9_7f4a_7c15u64.wrapping_mul(slot));
    mixed ^= mixed >> 30;
    mixed = mixed.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    mixed ^= mixed >> 27;
    mixed = mixed.wrapping_mul(0x94d0_49bb_1331_11eb);
    mixed ^ (mixed >> 31)
}

fn apply_iteration_core_states(
    integration_state: &mut IntegrationState,
    slots: &[IntegrationSlot],
    cores: usize,
    cur_points: usize,
    elapsed_seconds: f64,
    core_states: &[CoreIterationState],
) {
    let n_slots = integration_state.slot_metas.len();

    for core_state in core_states {
        integration_state.stats = integration_state.stats.merged(&core_state.stats);
        for (integral, core_integral) in integration_state
            .all_integrals
            .iter_mut()
            .zip(core_state.integrals.iter())
        {
            integral.merge(core_integral);
        }
    }

    let mut merged_sampling_grids = integration_state
        .sampling_states
        .iter()
        .map(|sampling_state| sampling_state.grid.clone_without_samples())
        .collect_vec();
    let mut merged_re_grids = (0..n_slots)
        .map(|slot_index| {
            integration_state
                .sampling_state_for_slot(slot_index)
                .grid
                .clone_without_samples()
        })
        .collect_vec();
    let mut merged_im_grids = (0..n_slots)
        .map(|slot_index| {
            integration_state
                .sampling_state_for_slot(slot_index)
                .grid
                .clone_without_samples()
        })
        .collect_vec();

    for core_state in core_states {
        for (sampling_state, merged_grid) in core_state
            .sampling_states
            .iter()
            .zip(merged_sampling_grids.iter_mut())
        {
            merged_grid
                .merge(&sampling_state.sampling_grid)
                .expect("could not merge grids");
        }
        for slot_index in 0..n_slots {
            merged_re_grids[slot_index]
                .merge(&core_state.slot_re_grids[slot_index])
                .expect("could not merge real accumulation grids");
            merged_im_grids[slot_index]
                .merge(&core_state.slot_im_grids[slot_index])
                .expect("could not merge imaginary accumulation grids");
        }
    }

    for (summary, grid) in integration_state
        .slot_re_summaries
        .iter_mut()
        .zip(merged_re_grids.iter())
    {
        if let Some(summary) = summary.as_mut() {
            summary.merge_iteration_grid(grid);
            summary.update_iter();
        }
    }
    for (summary, grid) in integration_state
        .slot_im_summaries
        .iter_mut()
        .zip(merged_im_grids.iter())
    {
        if let Some(summary) = summary.as_mut() {
            summary.merge_iteration_grid(grid);
            summary.update_iter();
        }
    }

    for (slot_index, merged_grid) in merged_sampling_grids.into_iter().enumerate() {
        let actual_slot_index = match integration_state.sampling_correlation_mode {
            SamplingCorrelationMode::Correlated => 0,
            SamplingCorrelationMode::Uncorrelated => slot_index,
        };
        let discrete_axis_labels = integration_state
            .sampling_state_for_slot(actual_slot_index)
            .discrete_axis_labels
            .clone();
        let sampling_state = integration_state.sampling_state_for_slot_mut(actual_slot_index);
        *sampling_state = SamplingSlotState {
            grid: merged_grid,
            discrete_axis_labels,
        };
        sampling_state.grid.update(
            F(slots[actual_slot_index]
                .settings
                .integrator
                .discrete_dim_learning_rate),
            F(slots[actual_slot_index]
                .settings
                .integrator
                .continuous_dim_learning_rate),
        );
    }

    integration_state.update_iter(false);
    integration_state.iter += 1;
    integration_state.elapsed_seconds = elapsed_seconds;
    integration_state.n_cores = cores;
    integration_state.num_points += cur_points;
}

fn build_preview_integration_state(
    integration_state: &IntegrationState,
    slots: &[IntegrationSlot],
    cores: usize,
    completed_points: usize,
    elapsed_seconds: f64,
    core_states: &[CoreIterationState],
) -> IntegrationState {
    let mut preview = integration_state.clone();
    apply_iteration_core_states(
        &mut preview,
        slots,
        cores,
        completed_points,
        elapsed_seconds,
        core_states,
    );
    preview
}

/// Integrate function used for local runs
pub fn havana_integrate<S>(
    request: HavanaIntegrateRequest,
    mut status_emitter: S,
) -> Result<IntegrationResult>
where
    S: FnMut(StatusUpdate) -> Result<()>,
{
    let HavanaIntegrateRequest {
        mut slots,
        sampling_correlation_mode,
        n_cores,
        state,
        workspace,
        output_control,
        batching,
        view_options,
    } = request;

    if slots.is_empty() {
        return Err(Report::msg(
            "At least one integrand must be selected for integration",
        ));
    }

    let slot_metas = slots.iter().map(|slot| slot.meta.clone()).collect_vec();
    let targets = slots.iter().map(|slot| slot.target).collect_vec();
    let sampling_states = match sampling_correlation_mode {
        SamplingCorrelationMode::Correlated => {
            vec![SamplingSlotState::from_integrand(&slots[0].integrand)]
        }
        SamplingCorrelationMode::Uncorrelated => slots
            .iter()
            .map(|slot| SamplingSlotState::from_integrand(&slot.integrand))
            .collect_vec(),
    };
    let monitored_discrete_setup =
        resolve_monitored_discrete_setup(sampling_correlation_mode, &slots, &sampling_states);
    let slot_first_non_trivial_discrete_breakdown_metadata =
        if let Some(path) = monitored_discrete_setup.path.as_deref() {
            slots
                .iter()
                .enumerate()
                .map(|(slot_index, slot)| {
                    build_persisted_discrete_breakdown_metadata(
                        &slot.integrand,
                        path,
                        &sampling_states[sampling_correlation_mode.state_index(slot_index)]
                            .discrete_axis_labels,
                    )
                })
                .collect_vec()
        } else {
            vec![None; slot_metas.len()]
        };
    if let Some(label_warning) = monitored_discrete_setup.warning.as_ref() {
        warn!("{label_warning}");
    }

    let mut integration_state = if let Some(integration_state) = state {
        integration_state
    } else {
        IntegrationState::new_from_settings(
            sampling_correlation_mode,
            sampling_states.clone(),
            slot_metas.clone(),
            monitored_discrete_setup.path.clone(),
            monitored_discrete_setup.axis_label.clone(),
            monitored_discrete_setup.bin_descriptions.clone(),
            slot_first_non_trivial_discrete_breakdown_metadata,
        )
    };
    integration_state.monitored_discrete_path = monitored_discrete_setup.path;
    integration_state.first_non_trivial_discrete_label = monitored_discrete_setup.axis_label;
    integration_state.first_non_trivial_discrete_bin_descriptions =
        monitored_discrete_setup.bin_descriptions;
    if integration_state.slot_metas != slot_metas {
        return Err(Report::msg(
            "Saved integration state slots do not match the currently selected integrands",
        ));
    }
    if integration_state.sampling_correlation_mode != sampling_correlation_mode {
        return Err(Report::msg(
            "Saved integration state sampling mode does not match the current integration mode",
        ));
    }
    if integration_state.sampling_states.len()
        != sampling_correlation_mode.state_count(slot_metas.len())
    {
        return Err(Report::msg(
            "Saved integration state sampling metadata is inconsistent with the selected slots",
        ));
    }
    if integration_state
        .slot_first_non_trivial_discrete_breakdown_metadata
        .len()
        != slot_metas.len()
    {
        return Err(Report::msg(
            "Saved integration state discrete breakdown metadata is inconsistent with the selected slots",
        ));
    }

    let primary = &slots[0];
    let sampling_str = primary.settings.sampling.describe_settings();
    let dimension = primary.integrand.get_n_dim();
    let discrete_depth = primary.settings.sampling.discrete_depth();
    let is_tropical_sampling = primary
        .settings
        .sampling
        .get_parameterization_settings()
        .is_none();
    let use_ltd = primary.settings.general.use_ltd;
    let integration_seed = primary.settings.integrator.seed;
    let integrated_phase = primary.settings.integrator.integrated_phase;
    let target_relative_accuracy = primary.settings.integrator.target_relative_accuracy;
    let target_absolute_accuracy = primary.settings.integrator.target_absolute_accuracy;
    let n_start = primary.settings.integrator.n_start;
    let n_increase = primary.settings.integrator.n_increase;
    let n_max = primary.settings.integrator.n_max;

    let cont_dim_str = if is_tropical_sampling {
        format!("a median continious dimension of {}", dimension)
    } else {
        format!("{} continuous dimensions", dimension)
    };

    let graph_string = if discrete_depth > 0 {
        let num_graphs = match &integration_state.sampling_state_for_slot(0).grid {
            Grid::Discrete(g) => g.bins.len(),
            _ => unreachable!(),
        };
        format!(
            "{discrete_depth} nested discrete grids with {} {} and ",
            num_graphs,
            if num_graphs > 1 { "graphs" } else { "graph" }
        )
    } else {
        String::new()
    };

    let correlation_str = match sampling_correlation_mode {
        SamplingCorrelationMode::Correlated => "correlated",
        SamplingCorrelationMode::Uncorrelated => "uncorrelated",
    };
    let grid_str = format!("{graph_string}{cont_dim_str} using {sampling_str} ({correlation_str})");

    let cores = n_cores.max(1);
    let n_slots = integration_state.slot_metas.len();
    let elapsed_seconds_offset = integration_state.elapsed_seconds;

    let t_start = Instant::now();

    info!(
        "Integrating using {} ltd with {} {} over {} ...",
        if use_ltd { "naive" } else { "cff" },
        cores,
        if cores > 1 { "cores" } else { "core" },
        grid_str
    );
    info!("");

    let pool = ThreadPoolBuilder::new().num_threads(cores).build().unwrap();

    let mut n_samples_evaluated = 0;
    let mut emitted_latest_observable_paths = vec![None; n_slots];
    clear_iteration_abort_request();
    'integrateLoop: while integration_state.num_points < n_max {
        // ensure we do not overshoot
        let cur_points = {
            let cur_points_not_final_iter = n_start + n_increase * integration_state.iter;
            if cur_points_not_final_iter + integration_state.num_points > n_max {
                n_max - integration_state.num_points
            } else {
                cur_points_not_final_iter
            }
        };

        // the number of points per core is the same for all cores, except for the last one
        let target_points_per_core = (cur_points - 1) / cores + 1;
        let n_points_per_core = (0..cores)
            .map(|core_id| {
                if core_id + 1 == cores {
                    cur_points - target_points_per_core * (cores - 1)
                } else {
                    target_points_per_core
                }
            })
            .collect_vec();

        let iteration_start = Instant::now();
        let mut current_batch_size = initial_batch_size(batching);
        let mut last_live_status = None;
        if batching.emit_initial_status_update && !is_interrupted() {
            status_emitter(build_status_update(StatusUpdateBuildRequest {
                kind: IntegrationStatusKind::Live,
                integration_state: &integration_state,
                cores,
                elapsed_time: t_start.elapsed(),
                iteration_elapsed_time: iteration_start.elapsed(),
                cur_points: 0,
                total_points_display: integration_state.num_points,
                n_samples_evaluated,
                targets: &targets,
                render_options: &view_options,
                live_progress: Some(status_update::LiveIterationProgress {
                    completed_points: 0,
                    target_points: cur_points,
                }),
            }))?;
        }

        let current_max_evals = integration_state
            .all_integrals
            .iter()
            .map(ComplexAccumulator::get_worst_case)
            .collect_vec();

        let mut worker_states = n_points_per_core
            .iter()
            .enumerate()
            .map(|(core_id, &n_points)| {
                CoreIterationState::new(
                    slots.iter().map(|slot| slot.integrand.clone()).collect(),
                    integration_state.sampling_correlation_mode,
                    &integration_state
                        .sampling_states
                        .iter()
                        .map(|sampling_state| sampling_state.grid.clone())
                        .collect_vec(),
                    integration_seed + integration_state.iter as u64,
                    target_points_per_core * core_id,
                    n_points,
                )
            })
            .collect_vec();
        while total_remaining_points(&worker_states) > 0 {
            let round_started_at = Instant::now();
            let processed_per_core: Vec<Result<usize>> = {
                let slot_settings = slots.iter().map(|slot| &slot.settings).collect_vec();
                let slot_models = slots.iter().map(|slot| &slot.model).collect_vec();
                pool.install(|| {
                    worker_states
                        .par_iter_mut()
                        .map(|worker_state| {
                            worker_state.evaluate_chunk(
                                &slot_settings,
                                &slot_models,
                                integration_state.iter,
                                &current_max_evals,
                                current_batch_size,
                            )
                        })
                        .collect()
                })
            };
            let processed_this_round = processed_per_core
                .into_iter()
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .sum::<usize>();

            if processed_this_round == 0 {
                break;
            }

            let completed_points = total_completed_points(&worker_states);
            if batching.emit_live_status_updates
                && completed_points < cur_points
                && !is_interrupted()
                && !is_iteration_abort_requested()
            {
                let now = Instant::now();
                let should_emit_live_status = last_live_status.is_none_or(|previous| {
                    now.duration_since(previous).as_secs_f64()
                        >= batching.min_time_between_status_updates_seconds
                });
                if should_emit_live_status {
                    let preview_state = build_preview_integration_state(
                        &integration_state,
                        &slots,
                        cores,
                        completed_points,
                        elapsed_seconds_offset + t_start.elapsed().as_secs_f64(),
                        &worker_states,
                    );
                    status_emitter(build_status_update(StatusUpdateBuildRequest {
                        kind: IntegrationStatusKind::Live,
                        integration_state: &preview_state,
                        cores,
                        elapsed_time: t_start.elapsed(),
                        iteration_elapsed_time: iteration_start.elapsed(),
                        cur_points: completed_points,
                        total_points_display: integration_state.num_points + completed_points,
                        n_samples_evaluated: n_samples_evaluated + completed_points,
                        targets: &targets,
                        render_options: &view_options,
                        live_progress: Some(status_update::LiveIterationProgress {
                            completed_points,
                            target_points: cur_points,
                        }),
                    }))?;
                    last_live_status = Some(now);
                }
            }

            if is_interrupted() {
                warn!("{}", "Integration interrupted by user".yellow());
                break 'integrateLoop;
            }

            if is_iteration_abort_requested() {
                break;
            }

            if total_remaining_points(&worker_states) > 0 {
                current_batch_size =
                    next_batch_size(batching, current_batch_size, round_started_at.elapsed());
            }
        }

        let abort_current_iteration = is_iteration_abort_requested();
        clear_iteration_abort_request();

        if is_interrupted() {
            warn!("{}", "Integration interrupted by user".yellow());
            break 'integrateLoop;
        }

        let completed_points = total_completed_points(&worker_states);
        if abort_current_iteration && completed_points < cur_points {
            if completed_points == 0 {
                warn!(
                    "{}",
                    "Current iteration abort requested before any samples completed; restarting the iteration."
                        .yellow()
                );
                continue;
            }
            warn!(
                "{}",
                format!(
                    "Current iteration aborted by user after {} completed samples.",
                    completed_points
                )
                .yellow()
            );
            continue;
        }

        n_samples_evaluated += completed_points;
        apply_iteration_core_states(
            &mut integration_state,
            &slots,
            cores,
            completed_points,
            elapsed_seconds_offset + t_start.elapsed().as_secs_f64(),
            &worker_states,
        );

        // merge runtime state per slot into the owner core
        let (owner_core_slice, other_cores) = worker_states.split_at_mut(1);
        let owner_core = &mut owner_core_slice[0].slot_integrands;
        for other_core in other_cores {
            for (owner_slot, other_slot) in owner_core
                .iter_mut()
                .zip(other_core.slot_integrands.iter_mut())
            {
                owner_slot.merge_runtime_results(other_slot)?;
            }
        }
        for (slot, integrand) in slots.iter_mut().zip(owner_core.iter()) {
            slot.integrand = integrand.clone();
        }

        // Update observable accumulators, persist authoritative resume state, then refresh
        // the latest user-facing snapshots.
        for slot in &mut slots {
            slot.integrand
                .update_runtime_results(integration_state.iter);
        }

        if let Some(ref workspace_path) = workspace {
            for (slot_index, slot) in slots.iter().enumerate() {
                write_observable_resume_state(
                    &slot.integrand,
                    Some(workspace_path.as_path()),
                    &integration_state.slot_metas[slot_index],
                    integration_state.iter,
                    output_control,
                )?;
            }
            write_integration_state_to_workspace(workspace_path, &integration_state)?;
            write_integration_result_snapshots(
                workspace_path,
                &integration_state,
                &targets,
                output_control,
            )?;
        }

        for (slot_index, slot) in slots.iter().enumerate() {
            let workspace_path = workspace
                .as_deref()
                .map(|root| slot_workspace_path(root, &integration_state.slot_metas[slot_index]));
            emitted_latest_observable_paths[slot_index] =
                write_latest_observables_output(&slot.integrand, workspace_path.as_deref())?;
            write_observable_snapshot_archive(
                &slot.integrand,
                workspace_path.as_deref(),
                integration_state.iter,
                output_control,
            )?;
        }

        status_emitter(build_status_update(StatusUpdateBuildRequest {
            kind: IntegrationStatusKind::Iteration,
            integration_state: &integration_state,
            cores,
            elapsed_time: t_start.elapsed(),
            iteration_elapsed_time: iteration_start.elapsed(),
            cur_points: completed_points,
            total_points_display: integration_state.num_points,
            n_samples_evaluated,
            targets: &targets,
            render_options: &view_options,
            live_progress: None,
        }))?;

        let target_accuracy_status = evaluate_target_accuracy(
            &integration_state,
            integration_state.num_points,
            t_start.elapsed(),
            &targets,
            match integrated_phase {
                IntegratedPhase::Real | IntegratedPhase::Both => {
                    IntegrationStatusPhaseDisplay::Real
                }
                IntegratedPhase::Imag => IntegrationStatusPhaseDisplay::Imag,
            },
            target_relative_accuracy,
            target_absolute_accuracy,
        );
        if target_accuracy_status.is_reached() {
            let reached_target = match (
                target_accuracy_status.relative_reached,
                target_accuracy_status.absolute_reached,
            ) {
                (true, true) => "relative and absolute",
                (true, false) => "relative",
                (false, true) => "absolute",
                (false, false) => unreachable!(),
            };
            info!(
                "Stopping integration after reaching the configured {reached_target} accuracy target."
            );
            break;
        }
    }
    // Reset the interrupted flag
    set_interrupted(false);
    clear_iteration_abort_request();

    if integration_state.num_points > 0 {
        let final_view_options = view_options.for_final();
        status_emitter(build_status_update(StatusUpdateBuildRequest {
            kind: IntegrationStatusKind::Final,
            integration_state: &integration_state,
            cores,
            elapsed_time: t_start.elapsed(),
            iteration_elapsed_time: Duration::ZERO,
            cur_points: 0,
            total_points_display: integration_state.num_points,
            n_samples_evaluated,
            targets: &targets,
            render_options: &final_view_options,
            live_progress: None,
        }))?;
    } else {
        info!("");
        warn!(
            "{}",
            "No final integration results to display since no iteration completed.".yellow()
        );
        info!("");
    }

    if !slots.is_empty() {
        emit_results_output_summary(
            workspace.as_deref(),
            &slots,
            &emitted_latest_observable_paths,
            output_control,
        );
    }

    Ok(build_integration_result(&integration_state, &targets))
}

/// Batch integrate function used for distributed runs, used by the worker nodes.
/// Evaluates a batch of points and returns the results in a manner specified by the user.
pub(crate) fn batch_integrate(
    integrand: &mut Integrand,
    model: &Model,
    input: BatchIntegrateInput,
) -> Result<BatchResult> {
    let samples = match input.samples {
        SampleInput::SampleList { samples } => samples,
        SampleInput::Grid {
            mut grid,
            num_points,
            seed,
            thread_id,
        } => {
            let mut rng = MonteCarloRng::new(seed, thread_id);

            (0..num_points)
                .map(|_| {
                    let mut sample = Sample::new();
                    grid.sample(&mut rng, &mut sample);
                    sample
                })
                .collect_vec()
        }
    };

    let (evaluation_results, metadata_statistics) = evaluate_sample_list(
        integrand,
        &samples,
        model,
        input.num_cores,
        input.iter,
        input.max_eval,
    )?;

    let integrand_output = generate_integrand_output(
        integrand,
        &evaluation_results,
        &samples,
        input.integrand_output_settings,
        input.settings.integrator.integrated_phase,
    );

    let event_output = generate_event_output(
        integrand,
        evaluation_results,
        input.event_output_settings,
        input.settings,
    );

    Ok(BatchResult {
        statistics: metadata_statistics,
        integrand_data: integrand_output,
        event_data: event_output,
    })
}

/// Map the evaluation result on to the right output specified by the user.
fn generate_integrand_output(
    integrand: &Integrand,
    evaluation_results: &[EvaluationResult],
    samples: &[Sample<F<f64>>],
    integrand_output_settings: IntegralOutputSettings,
    integrated_phase: IntegratedPhase,
) -> BatchIntegrateOutput {
    fn effective_integrand_result(result: &EvaluationResult) -> Complex<F<f64>> {
        let jacobian = result.parameterization_jacobian.unwrap_or(F(1.0));
        result.integrand_result * Complex::new_re(jacobian)
    }

    match integrand_output_settings {
        IntegralOutputSettings::Default => {
            let integrand_values = evaluation_results
                .iter()
                .map(effective_integrand_result)
                .collect();

            BatchIntegrateOutput::Default(integrand_values, samples.to_vec())
        }
        IntegralOutputSettings::Accumulator => {
            let mut real_accumulator = StatisticsAccumulator::new();
            let mut imag_accumulator = StatisticsAccumulator::new();
            let mut grid = integrand.create_grid();

            for (result, sample) in evaluation_results.iter().zip(samples.iter()) {
                let effective_result = effective_integrand_result(result);
                real_accumulator
                    .add_sample(effective_result.re * sample.get_weight(), Some(sample));
                imag_accumulator
                    .add_sample(effective_result.im * sample.get_weight(), Some(sample));

                match integrated_phase {
                    IntegratedPhase::Real => {
                        grid.add_training_sample(sample, effective_result.re)
                            .unwrap();
                    }
                    IntegratedPhase::Imag => {
                        grid.add_training_sample(sample, effective_result.im)
                            .unwrap();
                    }
                    IntegratedPhase::Both => {
                        unimplemented!()
                    }
                }
            }

            BatchIntegrateOutput::Accumulator(
                Box::new((real_accumulator, imag_accumulator)),
                Box::new(grid),
            )
        }
    }
}

/// Process events into histograms or event lists, as specified by the user.
fn generate_event_output(
    integrand: &Integrand,
    evaluation_results: Vec<EvaluationResult>,
    event_output_settings: EventOutputSettings,
    _settings: &RuntimeSettings,
) -> EventOutput {
    match event_output_settings {
        EventOutputSettings::None => EventOutput::None,

        EventOutputSettings::EventList => {
            let mut event_groups = EventGroupList::default();
            for mut result in evaluation_results {
                event_groups.append(&mut result.event_groups);
            }
            EventOutput::EventList { event_groups }
        }
        EventOutputSettings::Histogram => integrand
            .observable_accumulator_bundle()
            .map(|histograms| EventOutput::Histogram { histograms })
            .unwrap_or(EventOutput::None),
    }
}

pub fn slot_workspace_path(workspace: &Path, slot_meta: &SlotMeta) -> PathBuf {
    workspace.join("integrands").join(slot_meta.key())
}

pub fn workspace_manifest_path(workspace: &Path) -> PathBuf {
    workspace.join("manifest.json")
}

pub fn workspace_state_dir(workspace: &Path) -> PathBuf {
    workspace.join("state")
}

pub fn workspace_state_path(workspace: &Path) -> PathBuf {
    workspace_state_dir(workspace).join("integration_state.bin")
}

fn workspace_observable_state_dir(workspace: &Path, slot_meta: &SlotMeta) -> PathBuf {
    workspace_state_dir(workspace)
        .join("observables")
        .join(slot_meta.key())
}

pub fn latest_observable_resume_state_path(workspace: &Path, slot_meta: &SlotMeta) -> PathBuf {
    workspace_observable_state_dir(workspace, slot_meta).join("latest.json")
}

fn archived_observable_resume_state_path(
    workspace: &Path,
    slot_meta: &SlotMeta,
    iter: usize,
) -> PathBuf {
    workspace_observable_state_dir(workspace, slot_meta).join(format!("iter_{iter:04}.json"))
}

pub fn workspace_result_snapshot_path(workspace: &Path) -> PathBuf {
    workspace.join("integration_result.json")
}

fn workspace_result_archive_path(workspace: &Path, iter: usize) -> PathBuf {
    workspace
        .join("results")
        .join(format!("integration_result_iter_{iter:04}.json"))
}

fn observable_output_extension(format: ObservableFileFormat) -> Option<&'static str> {
    match format {
        ObservableFileFormat::None => None,
        ObservableFileFormat::Hwu => Some("hwu"),
        ObservableFileFormat::Json => Some("json"),
    }
}

fn latest_observable_output_path(
    workspace: &Path,
    format: ObservableFileFormat,
) -> Option<PathBuf> {
    let extension = observable_output_extension(format)?;
    Some(workspace.join(format!("observables_final.{extension}")))
}

fn archived_observable_output_path(
    workspace: &Path,
    format: ObservableFileFormat,
    iter: usize,
) -> Option<PathBuf> {
    let extension = observable_output_extension(format)?;
    Some(workspace.join(format!("observables_final_iter_{iter:04}.{extension}")))
}

fn observable_iteration_output_pattern(
    workspace: &Path,
    format: ObservableFileFormat,
) -> Option<String> {
    let extension = observable_output_extension(format)?;
    Some(
        workspace
            .join(format!("observables_final_iter_<iter>.{extension}"))
            .display()
            .to_string(),
    )
}

fn user_facing_observables_output_enabled(integrand: &Integrand) -> bool {
    let Integrand::ProcessIntegrand(process_integrand) = integrand else {
        return false;
    };
    process_integrand
        .get_settings()
        .integrator
        .observables_output
        .format
        != ObservableFileFormat::None
}

fn write_atomic_bytes(path: &Path, bytes: &[u8]) -> Result<()> {
    let parent = path
        .parent()
        .ok_or_else(|| Report::msg("Atomic write target is missing a parent directory"))?;
    fs::create_dir_all(parent)?;
    let tmp_path = path.with_extension(format!(
        "{}.tmp",
        path.extension()
            .and_then(|ext| ext.to_str())
            .unwrap_or("tmp")
    ));
    fs::write(&tmp_path, bytes)?;
    fs::rename(&tmp_path, path)?;
    Ok(())
}

fn write_atomic_json<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    write_atomic_bytes(path, &serde_json::to_vec_pretty(value)?)
}

fn write_observable_resume_state(
    integrand: &Integrand,
    workspace: Option<&Path>,
    slot_meta: &SlotMeta,
    iter: usize,
    output_control: WorkspaceSnapshotControl,
) -> Result<Option<PathBuf>> {
    let Some(bundle) = integrand.observable_snapshot_bundle() else {
        return Ok(None);
    };
    let Some(workspace) = workspace else {
        return Ok(None);
    };
    let latest_path = latest_observable_resume_state_path(workspace, slot_meta);
    write_atomic_json(&latest_path, &bundle)?;
    if output_control.write_iteration_archives {
        let archived_path = archived_observable_resume_state_path(workspace, slot_meta, iter);
        write_atomic_json(&archived_path, &bundle)?;
    }
    Ok(Some(latest_path))
}

fn write_observable_snapshot_archive(
    integrand: &Integrand,
    workspace: Option<&Path>,
    iter: usize,
    output_control: WorkspaceSnapshotControl,
) -> Result<Option<PathBuf>> {
    let Integrand::ProcessIntegrand(process_integrand) = integrand else {
        return Ok(None);
    };
    if !output_control.write_iteration_archives {
        return Ok(None);
    }
    let format = process_integrand
        .get_settings()
        .integrator
        .observables_output
        .format;
    let Some(workspace) = workspace else {
        return Ok(None);
    };
    let Some(path) = archived_observable_output_path(workspace, format, iter) else {
        return Ok(None);
    };

    integrand.write_observable_snapshots(&path, format)?;
    Ok(Some(path))
}

fn write_latest_observables_output(
    integrand: &Integrand,
    workspace: Option<&Path>,
) -> Result<Option<PathBuf>> {
    let Integrand::ProcessIntegrand(process_integrand) = integrand else {
        return Ok(None);
    };
    if !user_facing_observables_output_enabled(integrand) {
        return Ok(None);
    }
    let format = process_integrand
        .get_settings()
        .integrator
        .observables_output
        .format;
    let Some(workspace) = workspace else {
        return Ok(None);
    };
    let Some(path) = latest_observable_output_path(workspace, format) else {
        return Ok(None);
    };

    integrand.write_observable_snapshots(&path, format)?;
    Ok(Some(path))
}

fn write_integration_state_to_workspace(
    workspace_path: &Path,
    integration_state: &IntegrationState,
) -> Result<()> {
    write_atomic_bytes(
        &workspace_state_path(workspace_path),
        &bincode::encode_to_vec(integration_state, bincode::config::standard())
            .unwrap_or_else(|_| panic!("failed to serialize the integration state")),
    )?;
    Ok(())
}

fn write_integration_result_snapshots(
    workspace_path: &Path,
    integration_state: &IntegrationState,
    targets: &[Option<Complex<F<f64>>>],
    output_control: WorkspaceSnapshotControl,
) -> Result<()> {
    let result = build_integration_result(integration_state, targets);
    write_atomic_json(&workspace_result_snapshot_path(workspace_path), &result)?;
    if output_control.write_iteration_archives {
        write_atomic_json(
            &workspace_result_archive_path(workspace_path, integration_state.iter),
            &result,
        )?;
    }
    Ok(())
}

fn workspace_relative_display_path(root: &Path, path: &Path) -> String {
    path.strip_prefix(root)
        .unwrap_or(path)
        .display()
        .to_string()
}

fn emit_results_output_summary(
    workspace: Option<&Path>,
    slots: &[IntegrationSlot],
    emitted_paths: &[Option<PathBuf>],
    output_control: WorkspaceSnapshotControl,
) {
    let mut summary_rows = Vec::new();
    for (slot, emitted_path) in slots.iter().zip(emitted_paths.iter()) {
        let Some(final_path) = emitted_path else {
            continue;
        };
        let Integrand::ProcessIntegrand(process_integrand) = &slot.integrand else {
            continue;
        };
        let format = process_integrand
            .get_settings()
            .integrator
            .observables_output
            .format;
        let iteration_pattern = output_control
            .write_iteration_archives
            .then(|| {
                observable_iteration_output_pattern(
                    final_path.parent().unwrap_or_else(|| Path::new(".")),
                    format,
                )
            })
            .flatten();
        summary_rows.push((&slot.meta, final_path, iteration_pattern));
    }

    if summary_rows.is_empty() && workspace.is_none() {
        return;
    }

    info!("");
    info!("{}", "Integration results emitted:".green().bold());
    if let Some(workspace) = workspace {
        info!(
            "results -> {}",
            workspace_relative_display_path(workspace, &workspace_result_snapshot_path(workspace))
        );
        if output_control.write_iteration_archives {
            info!(
                "results iteration snapshots -> {}",
                workspace_relative_display_path(
                    workspace,
                    &workspace_result_archive_path(workspace, 1),
                )
                .replace("iter_0001", "iter_*")
            );
        }
    }
    for (slot_meta, final_path, iteration_pattern) in summary_rows {
        let display_path = workspace
            .map(|root| workspace_relative_display_path(root, final_path))
            .unwrap_or_else(|| final_path.display().to_string());
        info!("{} -> {}", slot_key_label(slot_meta), display_path);
        if let Some(iteration_pattern) = iteration_pattern {
            let display_pattern = workspace
                .map(|root| workspace_relative_display_path(root, Path::new(&iteration_pattern)))
                .unwrap_or(iteration_pattern);
            info!(
                "{} iteration snapshots -> {}",
                format!("{} ", slot_key_label(slot_meta)).dimmed(),
                display_pattern
            );
        }
    }
}

/// This function actually evaluates the list of samples in parallel.
fn evaluate_sample_list(
    integrand: &mut Integrand,
    samples: &[Sample<F<f64>>],
    model: &Model,
    num_cores: usize,
    iter: usize,
    max_eval: Complex<F<f64>>,
) -> Result<(Vec<EvaluationResult>, StatisticsCounter)> {
    let list_size = samples.len();
    let nvec_per_core = (list_size - 1) / num_cores + 1;
    let num_chunks = samples.chunks(nvec_per_core).len();

    let sample_chunks = samples.par_chunks(nvec_per_core);
    let integrands = (0..num_chunks)
        .map(|_| integrand.clone())
        .collect_vec()
        .into_par_iter();

    let evaluation_results_per_core: Vec<
        Result<(Vec<EvaluationResult>, StatisticsCounter, Integrand)>,
    > = sample_chunks
        .zip(integrands)
        .map(|(chunk, mut integrand)| {
            let raw_batch =
                integrand.evaluate_samples_raw(chunk, model, iter, false, false, max_eval)?;
            Ok((raw_batch.samples, raw_batch.statistics, integrand))
        })
        .collect();
    let mut evaluation_results_per_core: Vec<(
        Vec<EvaluationResult>,
        StatisticsCounter,
        Integrand,
    )> = evaluation_results_per_core
        .into_iter()
        .collect::<Result<_>>()?;

    for (_, _, worker_integrand) in evaluation_results_per_core.iter_mut() {
        integrand.merge_runtime_results(worker_integrand)?;
    }

    let (evaluation_results, meta_data_statistics) = evaluation_results_per_core.into_iter().fold(
        (Vec::new(), StatisticsCounter::new_empty()),
        |(mut all_results, stats), (results, worker_stats, _)| {
            all_results.extend(results);
            (all_results, stats.merged(&worker_stats))
        },
    );

    Ok((evaluation_results, meta_data_statistics))
}

/// Different ways of passing samples to the batch_integrate function
/// One is simply a list of samples, the other is a grid and the desired number of samples
/// The worker then generates the samples itself.
#[derive(Serialize, Deserialize)]
pub enum SampleInput {
    SampleList {
        samples: Vec<Sample<F<f64>>>,
    },
    Grid {
        grid: Box<Grid<F<f64>>>,
        num_points: usize,
        seed: u64,
        thread_id: usize,
    },
}

type AccumulatorPair = (StatisticsAccumulator<F<f64>>, StatisticsAccumulator<F<f64>>);

#[derive(Serialize, Deserialize)]
pub enum BatchIntegrateOutput {
    Default(Vec<Complex<F<f64>>>, Vec<Sample<F<f64>>>),
    Accumulator(Box<AccumulatorPair>, Box<Grid<F<f64>>>),
}

/// Different ways of processing events, EventList is a list of events, Histogram does accumulation of events on the worker nodes, so the
/// master node only has to merge the histograms.
#[derive(Serialize, Deserialize)]
pub enum EventOutput {
    None,
    EventList {
        event_groups: EventGroupList,
    },
    Histogram {
        histograms: ObservableAccumulatorBundle,
    },
}

/// The result of evaluating a batch of points
#[derive(Serialize, Deserialize, Encode, Decode)]
pub struct BatchResult {
    pub statistics: StatisticsCounter,
    #[bincode(with_serde)]
    pub integrand_data: BatchIntegrateOutput,
    #[bincode(with_serde)]
    pub event_data: EventOutput,
}

/// Input for the batch_integrate function, created by the master node
pub struct BatchIntegrateInput<'a> {
    // global run info:
    pub max_eval: Complex<F<f64>>,
    pub iter: usize,
    pub settings: &'a RuntimeSettings,
    // input data:
    pub samples: SampleInput,
    pub integrand_output_settings: IntegralOutputSettings,
    pub event_output_settings: EventOutputSettings,
    pub num_cores: usize,
}

/// Choose whether to output the integrand results, or a statistics accumulator
#[derive(Serialize, Deserialize)]
pub enum IntegralOutputSettings {
    Default,
    Accumulator,
}

/// Choose whether to output an eventlist, a histogram, or nothing
#[derive(Serialize, Deserialize)]
pub enum EventOutputSettings {
    None,
    EventList,
    Histogram,
}

#[derive(Serialize, Deserialize, Encode, Decode)]
pub struct SerializableBatchIntegrateInput {
    #[bincode(with_serde)]
    pub max_eval: Complex<F<f64>>,
    pub iter: usize,
    #[bincode(with_serde)]
    pub samples: SampleInput,
    #[bincode(with_serde)]
    pub integrand_output_settings: IntegralOutputSettings,
    #[bincode(with_serde)]
    pub event_output_settings: EventOutputSettings,
    pub num_cores: usize,
}

impl SerializableBatchIntegrateInput {
    pub(crate) fn into_batch_integrate_input(
        self,
        settings: &RuntimeSettings,
    ) -> BatchIntegrateInput<'_> {
        BatchIntegrateInput {
            max_eval: self.max_eval,
            iter: self.iter,
            samples: self.samples,
            integrand_output_settings: self.integrand_output_settings,
            event_output_settings: self.event_output_settings,
            num_cores: self.num_cores,
            settings,
        }
    }
}

/// Master node which accumulates data from jobs, and creates the input for jobs
#[derive(Clone)]
pub struct MasterNode {
    grid: Grid<F<f64>>,
    integrator_settings: IntegratorSettings,
    master_accumulator_re: StatisticsAccumulator<F<f64>>,
    master_accumulator_im: StatisticsAccumulator<F<f64>>,
    statistics: StatisticsCounter,
    observable_accumulators: Option<ObservableAccumulatorBundle>,
    event_groups: EventGroupList,
    current_iter: usize,
}

impl MasterNode {
    pub(crate) fn new(grid: Grid<F<f64>>, integrator_settings: IntegratorSettings) -> Self {
        MasterNode {
            grid,
            integrator_settings,
            master_accumulator_im: StatisticsAccumulator::new(),
            master_accumulator_re: StatisticsAccumulator::new(),
            statistics: StatisticsCounter::new_empty(),
            observable_accumulators: None,
            event_groups: EventGroupList::default(),
            current_iter: 0,
        }
    }

    /// Update the grid with the data from another grid.
    fn update_grid_with_grid(&mut self, other_grid: &Grid<F<f64>>) -> Result<(), String> {
        self.grid.merge(other_grid)
    }

    /// Update the grid with the data from a set of samples.
    fn update_grid_with_samples(
        &mut self,
        samples_points: &[Sample<F<f64>>],
        results: &[Complex<F<f64>>],
    ) -> Result<(), String> {
        let integrated_phase = self.integrator_settings.integrated_phase;

        for (sample_point, result) in samples_points.iter().zip(results.iter()) {
            match integrated_phase {
                IntegratedPhase::Real => self.grid.add_training_sample(sample_point, result.re)?,
                IntegratedPhase::Imag => self.grid.add_training_sample(sample_point, result.im)?,
                IntegratedPhase::Both => {
                    unimplemented!("integrated phase both not yet implemented")
                }
            }
        }

        Ok(())
    }

    /// Update the accumulators with the data from another set of accumulators.
    pub(crate) fn update_accumulators_with_accumulators(
        &mut self,
        mut real_accumulator: StatisticsAccumulator<F<f64>>,
        mut imaginary_accumulator: StatisticsAccumulator<F<f64>>,
    ) {
        self.master_accumulator_re
            .merge_samples(&mut real_accumulator);

        self.master_accumulator_im
            .merge_samples(&mut imaginary_accumulator);
    }

    /// Update the accumulators with the data from a set of samples.
    fn update_accumuators_with_samples(
        &mut self,
        sample_points: &[Sample<F<f64>>],
        results: &[Complex<F<f64>>],
    ) {
        for (sample_point, result) in sample_points.iter().zip(results.iter()) {
            self.master_accumulator_re
                .add_sample(result.re * sample_point.get_weight(), Some(sample_point));

            self.master_accumulator_im
                .add_sample(result.im * sample_point.get_weight(), Some(sample_point));
        }
    }

    /// Update the metadata statistics with the data from another set of metadata statistics.
    fn update_metadata_statistics(&mut self, statistics: StatisticsCounter) {
        self.statistics = self.statistics.merged(&statistics);
    }

    /// Finish the current iteration. This should be called after all jobs have been processed.
    pub(crate) fn update_iter(&mut self) {
        self.grid.update(
            F(self.integrator_settings.discrete_dim_learning_rate),
            F(self.integrator_settings.continuous_dim_learning_rate),
        );
        self.master_accumulator_re.update_iter(false);
        self.master_accumulator_im.update_iter(false);
        if let Some(observable_accumulators) = self.observable_accumulators.as_mut() {
            observable_accumulators.update_results();
        }

        self.current_iter += 1;
    }

    /// Write the input for a batch job to a file.
    pub(crate) fn write_batch_input(
        &mut self,
        num_cores: usize,
        num_samples: usize,
        export_grid: bool,
        output_accumulator: bool,
        workspace_path: &str,
        job_id: usize,
    ) -> Result<(), Report> {
        let max_eval = Complex::new(
            self.master_accumulator_re
                .max_eval_positive
                .max(self.master_accumulator_re.max_eval_negative),
            self.master_accumulator_im
                .max_eval_positive
                .max(self.master_accumulator_im.max_eval_negative),
        );

        let samples = if export_grid {
            SampleInput::Grid {
                grid: Box::new(self.grid.clone()),
                num_points: num_samples,
                seed: self.integrator_settings.seed,
                thread_id: job_id,
            }
        } else {
            let mut rng = rand::rng();
            let mut samples_temp = vec![Sample::new(); num_samples];
            for sample in samples_temp.iter_mut() {
                self.grid.sample(&mut rng, sample);
            }
            SampleInput::SampleList {
                samples: samples_temp,
            }
        };

        let integrand_output_settings = if output_accumulator {
            IntegralOutputSettings::Accumulator
        } else {
            IntegralOutputSettings::Default
        };

        let input = SerializableBatchIntegrateInput {
            num_cores,
            max_eval,
            iter: self.current_iter,
            samples,
            integrand_output_settings,
            event_output_settings: EventOutputSettings::None,
        };

        let input_bytes = bincode::encode_to_vec(&input, bincode::config::standard())?;
        let job_name = format!("job_{}", job_id);
        let job_path = std::path::Path::new(workspace_path).join(job_name);

        std::fs::write(job_path, input_bytes)?;

        Ok(())
    }

    /// Process the output of a batch job.
    pub(crate) fn process_batch_output(&mut self, output: BatchResult) -> Result<(), String> {
        self.update_metadata_statistics(output.statistics);

        match output.integrand_data {
            BatchIntegrateOutput::Default(results, samples) => {
                self.update_accumuators_with_samples(&samples, &results);
                self.update_grid_with_samples(&samples, &results)?;
            }
            BatchIntegrateOutput::Accumulator(accumulators, grid) => {
                let (real_accumulator, imag_accumulator) = *accumulators;
                self.update_accumulators_with_accumulators(real_accumulator, imag_accumulator);
                self.update_grid_with_grid(&grid)?;
            }
        }

        match output.event_data {
            EventOutput::None => {}
            EventOutput::EventList { mut event_groups } => {
                self.event_groups.append(&mut event_groups);
            }
            EventOutput::Histogram { mut histograms } => {
                if let Some(existing) = self.observable_accumulators.as_mut() {
                    existing
                        .merge_samples(&mut histograms)
                        .map_err(|err| err.to_string())?;
                } else {
                    self.observable_accumulators = Some(histograms);
                }
            }
        }

        Ok(())
    }

    /// Display the current status of the integration. Usually called after each iteration.
    pub(crate) fn display_status(&self) {
        let status_block = [
            render_integral_result(
                &self.master_accumulator_re,
                "itg",
                self.current_iter,
                "re",
                None,
            ),
            render_integral_result(
                &self.master_accumulator_im,
                "itg",
                self.current_iter,
                "im",
                None,
            ),
            self.statistics.render_status_table(),
        ]
        .join("\n");

        info!("\n{status_block}");
    }
}

pub fn emit_integration_status_via_tracing(
    kind: IntegrationStatusKind,
    status_block: impl AsRef<str>,
) -> Result<()> {
    let status_block = status_block.as_ref();
    match kind {
        IntegrationStatusKind::Live => {}
        IntegrationStatusKind::Iteration => {
            info!("\n{status_block}");
            info!("");
        }
        IntegrationStatusKind::Final => {
            info!("");
            info!("{}", "Final integration results:".bold().green());
            info!("");
            info!("\n{status_block}");
            info!("");
        }
    }

    Ok(())
}

pub fn render_saved_integration_summary(
    integration_state: &IntegrationState,
    targets: &[Option<Complex<F<f64>>>],
    view_options: &IntegrationStatusViewOptions,
    tabled_options: &TabledRenderOptions,
) -> String {
    render_tabled::render_status_update(
        &build_saved_status_update(integration_state, targets, view_options),
        tabled_options,
    )
}

pub fn render_status_update_tabled(
    update: &StatusUpdate,
    tabled_options: &TabledRenderOptions,
) -> String {
    render_tabled::render_status_update(update, tabled_options)
}

pub fn build_integration_result(
    integration_state: &IntegrationState,
    targets: &[Option<Complex<F<f64>>>],
) -> IntegrationResult {
    let integration_statistics = integration_state.stats.snapshot();
    let slots = integration_state
        .slot_metas
        .iter()
        .enumerate()
        .map(|(slot_index, slot_meta)| {
            let accumulator = &integration_state.all_integrals[slot_index];
            let discrete_context =
                integration_state.monitored_discrete_context_for_slot(slot_index);
            SlotIntegrationResult {
                key: slot_meta.key(),
                process: slot_meta.process_name.clone(),
                integrand: slot_meta.integrand_name.clone(),
                target: targets[slot_index],
                integral: IntegralEstimate {
                    neval: accumulator.re.processed_samples,
                    real_zero: accumulator.re.num_zero_evaluations,
                    im_zero: accumulator.im.num_zero_evaluations,
                    result: Complex::new(accumulator.re.avg, accumulator.im.avg),
                    error: Complex::new(accumulator.re.err, accumulator.im.err),
                    real_chisq: accumulator.re.chi_sq,
                    im_chisq: accumulator.im.chi_sq,
                },
                table_results: build_table_result_summary(
                    slot_meta,
                    accumulator,
                    integration_state.iter,
                    targets[slot_index],
                ),
                integration_statistics: integration_statistics.clone(),
                max_weight_info: build_max_weight_info_summary(
                    &integration_state
                        .sampling_state_for_slot(slot_index)
                        .discrete_axis_labels,
                    accumulator,
                ),
                grid_breakdown: ComponentDiscreteBreakdown {
                    re: discrete_context.as_ref().and_then(|context| {
                        integration_state.slot_re_summaries[slot_index]
                            .as_ref()
                            .zip(
                                integration_state
                                    .slot_first_non_trivial_discrete_breakdown_metadata[slot_index]
                                    .as_ref(),
                            )
                            .and_then(|(summary, metadata)| {
                                summary_at_path(summary, &context.path).and_then(|summary| {
                                    summary.first_non_trivial_breakdown(metadata, &context.pdfs)
                                })
                            })
                    }),
                    im: discrete_context.as_ref().and_then(|context| {
                        integration_state.slot_im_summaries[slot_index]
                            .as_ref()
                            .zip(
                                integration_state
                                    .slot_first_non_trivial_discrete_breakdown_metadata[slot_index]
                                    .as_ref(),
                            )
                            .and_then(|(summary, metadata)| {
                                summary_at_path(summary, &context.path).and_then(|summary| {
                                    summary.first_non_trivial_breakdown(metadata, &context.pdfs)
                                })
                            })
                    }),
                },
            }
        })
        .collect();

    IntegrationResult { slots }
}

pub fn print_integral_result(
    itg: &StatisticsAccumulator<F<f64>>,
    label: &str,
    i_iter: usize,
    tag: &str,
    trgt: Option<F<f64>>,
) {
    info!("{}", render_integral_result(itg, label, i_iter, tag, trgt));
}

#[allow(clippy::format_in_format_args)]
fn render_integral_result(
    itg: &StatisticsAccumulator<F<f64>>,
    label: &str,
    i_iter: usize,
    tag: &str,
    trgt: Option<F<f64>>,
) -> String {
    let slot_meta = SlotMeta {
        process_name: label.to_string(),
        integrand_name: String::new(),
    };
    let mut cells = build_integral_result_cells(itg, &slot_meta, i_iter, tag, trgt);
    cells.integrand = format!("{label} {}:", format!("{:-2}", tag).blue().bold());
    let delta = match (cells.delta_sigma.as_ref(), cells.delta_percent.as_ref()) {
        (Some(delta_sigma), Some(delta_percent)) => format!("{delta_sigma}, {delta_percent}"),
        _ => String::new(),
    };

    format!(
        "|  {} {} {} {} {} {}",
        cells.integrand, cells.value, cells.relative_error, cells.chi_sq, delta, cells.mwi
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use colored::control;
    use ratatui::{Terminal, backend::TestBackend};
    use symbolica::numerical_integration::ContinuousGrid;

    fn make_accumulator(
        re_avg: f64,
        re_err: f64,
        re_chi_sq: f64,
        im_avg: f64,
        im_err: f64,
        im_chi_sq: f64,
    ) -> ComplexAccumulator {
        let mut accumulator = ComplexAccumulator::new();
        accumulator.re.avg = F(re_avg);
        accumulator.re.err = F(re_err);
        accumulator.re.chi_sq = F(re_chi_sq);
        accumulator.re.processed_samples = 100_000;
        accumulator.re.max_eval_positive = F(1.0);
        accumulator.im.avg = F(im_avg);
        accumulator.im.err = F(im_err);
        accumulator.im.chi_sq = F(im_chi_sq);
        accumulator.im.processed_samples = 100_000;
        accumulator.im.max_eval_positive = F(1.0);
        accumulator
    }

    fn make_integration_state() -> IntegrationState {
        let sampling_grid = Grid::Continuous(ContinuousGrid::new(1, 64, 100, None, false));
        let mut state = IntegrationState::new_from_settings(
            SamplingCorrelationMode::Correlated,
            vec![SamplingSlotState::new(sampling_grid, Vec::new())],
            vec![
                SlotMeta {
                    process_name: "proc_a".to_string(),
                    integrand_name: "itg_a".to_string(),
                },
                SlotMeta {
                    process_name: "proc_b".to_string(),
                    integrand_name: "itg_b".to_string(),
                },
            ],
            None,
            None,
            None,
            vec![None, None],
        );
        state.iter = 1;
        state.num_points = 100_000;
        let mut evaluation = EvaluationResult::zero();
        evaluation.evaluation_metadata.integrand_evaluation_time = Duration::from_micros(414);
        evaluation.evaluation_metadata.evaluator_evaluation_time = Duration::from_micros(279);
        evaluation.evaluation_metadata.parameterization_time = Duration::from_nanos(5_800);
        evaluation.evaluation_metadata.total_timing = Duration::from_micros(462);
        state.stats = StatisticsCounter::from_evaluation_results(&[evaluation]);
        state.all_integrals = vec![
            make_accumulator(7.5e-5, 9.8e-5, 0.394, 3.2e-5, 1.5e-5, 0.378),
            make_accumulator(2.1e-5, 2.0e-6, 0.221, -1.7e-5, 3.0e-6, 0.187),
        ];
        state
    }

    fn make_discrete_integration_state() -> IntegrationState {
        let make_continuous_grid =
            || Grid::Continuous(ContinuousGrid::new(1, 64, 100, None, false));
        let sampling_grid = Grid::Discrete(DiscreteGrid::new(
            vec![Some(make_continuous_grid()), Some(make_continuous_grid())],
            F(10.0),
            false,
        ));
        let mut state = IntegrationState::new_from_settings(
            SamplingCorrelationMode::Correlated,
            vec![SamplingSlotState::new(
                sampling_grid,
                vec!["graph".to_string()],
            )],
            vec![
                SlotMeta {
                    process_name: "proc_a".to_string(),
                    integrand_name: "itg_a".to_string(),
                },
                SlotMeta {
                    process_name: "proc_b".to_string(),
                    integrand_name: "itg_b".to_string(),
                },
            ],
            Some(vec![]),
            Some("graph".to_string()),
            Some(vec!["GL0".to_string(), "GL1".to_string()]),
            vec![
                Some(PersistedDiscreteBreakdownMetadata {
                    axis_label: "graph".to_string(),
                    fixed_coordinates: Vec::new(),
                    bin_labels: vec!["GL0".to_string(), "GL1".to_string()],
                }),
                Some(PersistedDiscreteBreakdownMetadata {
                    axis_label: "graph".to_string(),
                    fixed_coordinates: Vec::new(),
                    bin_labels: vec!["GL0".to_string(), "GL1".to_string()],
                }),
            ],
        );
        state.first_non_trivial_discrete_bin_descriptions =
            Some(vec!["GL0".to_string(), "GL1".to_string()]);
        state.iter = 2;
        state.num_points = 210_000;
        let mut evaluation = EvaluationResult::zero();
        evaluation.evaluation_metadata.integrand_evaluation_time = Duration::from_micros(414);
        evaluation.evaluation_metadata.evaluator_evaluation_time = Duration::from_micros(279);
        evaluation.evaluation_metadata.parameterization_time = Duration::from_nanos(5_800);
        evaluation.evaluation_metadata.total_timing = Duration::from_micros(462);
        state.stats = StatisticsCounter::from_evaluation_results(&[evaluation]);
        state.all_integrals = vec![
            make_accumulator(7.5e-5, 9.8e-5, 0.394, 3.2e-5, 1.5e-5, 0.378),
            make_accumulator(2.1e-5, 2.0e-6, 0.221, -1.7e-5, 3.0e-6, 0.187),
        ];
        if let Grid::Discrete(discrete_grid) = &mut state.sampling_state_for_slot_mut(0).grid {
            discrete_grid.bins[0].pdf = F(0.75);
            discrete_grid.bins[1].pdf = F(0.25);
        }

        for summaries in [&mut state.slot_re_summaries, &mut state.slot_im_summaries] {
            for (slot_index, summary) in summaries.iter_mut().enumerate() {
                let summary = summary.as_mut().expect("discrete summary expected");
                summary.bins[0].accumulator.avg = F(1.0e-5 * (slot_index as f64 + 1.0));
                summary.bins[0].accumulator.err = F(2.0e-6);
                summary.bins[0].accumulator.chi_sq = F(0.2);
                summary.bins[0].accumulator.processed_samples = 150;
                summary.bins[0].accumulator.max_eval_positive = F(0.75);
                summary.bins[0].accumulator.max_eval_positive_xs =
                    Some(Sample::Continuous(F(1.0), vec![F(0.25)]));
                summary.bins[1].accumulator.avg = F(5.0e-6 * (slot_index as f64 + 1.0));
                summary.bins[1].accumulator.err = F(1.0e-6);
                summary.bins[1].accumulator.chi_sq = F(0.1);
                summary.bins[1].accumulator.processed_samples = 50;
                summary.bins[1].accumulator.max_eval_positive = F(0.25);
                summary.bins[1].accumulator.max_eval_positive_xs =
                    Some(Sample::Continuous(F(1.0), vec![F(0.75)]));
            }
        }

        state
    }

    fn default_view_options() -> IntegrationStatusViewOptions {
        IntegrationStatusViewOptions {
            phase_display: IntegrationStatusPhaseDisplay::Both,
            training_phase_display: IntegrationStatusPhaseDisplay::Real,
            training_slot: 0,
            slot_training_phase_displays: vec![IntegrationStatusPhaseDisplay::Real],
            per_slot_training_phase: false,
            target_relative_accuracy: None,
            target_absolute_accuracy: None,
            show_statistics: true,
            show_max_weight_details: true,
            show_top_discrete_grid: false,
            show_discrete_contributions_sum: false,
            contribution_sort: ContributionSortMode::Error,
            show_max_weight_info_for_discrete_bins: false,
        }
    }

    fn default_tabled_options() -> TabledRenderOptions {
        TabledRenderOptions {
            max_table_width: DEFAULT_MAX_SHARED_TABLE_WIDTH,
            show_statistics: true,
            show_max_weight_details: true,
            show_top_discrete_grid: false,
            show_discrete_contributions_sum: false,
            show_max_weight_info_for_discrete_bins: false,
        }
    }

    fn tabled_options_for_view(view_options: &IntegrationStatusViewOptions) -> TabledRenderOptions {
        TabledRenderOptions {
            show_statistics: view_options.show_statistics,
            show_max_weight_details: view_options.show_max_weight_details,
            show_top_discrete_grid: view_options.show_top_discrete_grid,
            show_discrete_contributions_sum: view_options.show_discrete_contributions_sum,
            show_max_weight_info_for_discrete_bins: view_options
                .show_max_weight_info_for_discrete_bins,
            ..default_tabled_options()
        }
    }

    fn render_update(request: StatusUpdateBuildRequest<'_>) -> String {
        let tabled_options = tabled_options_for_view(request.render_options);
        render_status_update_tabled(&build_status_update(request), &tabled_options)
    }

    fn render_ratatui_update(
        request: StatusUpdateBuildRequest<'_>,
        configure: impl FnOnce(&mut RatatuiDashboardState),
    ) -> String {
        let mut dashboard = RatatuiDashboardState::new();
        dashboard.update(build_status_update(request));
        configure(&mut dashboard);

        let backend = TestBackend::new(180, 48);
        let mut terminal = Terminal::new(backend).expect("test terminal");
        terminal
            .draw(|frame| dashboard.draw(frame))
            .expect("ratatui draw");

        let buffer = terminal.backend().buffer();
        (0..buffer.area.height)
            .map(|y| {
                let mut line = String::new();
                for x in 0..buffer.area.width {
                    line.push_str(buffer[(x, y)].symbol());
                }
                line.trim_end().to_string()
            })
            .collect::<Vec<_>>()
            .join("\n")
    }

    #[test]
    fn max_weight_details_table_renders_titled_table_without_wrapping_parentheses() {
        let mut accumulator = ComplexAccumulator::new();
        accumulator.re.max_eval_positive = F(2.3840672847728);
        accumulator.re.max_eval_positive_xs = Some(Sample::Discrete(
            F(1.0),
            0,
            Some(Box::new(Sample::Discrete(
                F(1.0),
                0,
                Some(Box::new(Sample::Continuous(
                    F(1.0),
                    vec![
                        F(0.9695746085826609),
                        F(0.5327714835649985),
                        F(0.003410284786539597),
                    ],
                ))),
            ))),
        ));

        let mut state = IntegrationState::new_from_settings(
            SamplingCorrelationMode::Correlated,
            vec![SamplingSlotState::new(
                Grid::Continuous(ContinuousGrid::new(1, 64, 100, None, false)),
                Vec::new(),
            )],
            vec![SlotMeta {
                process_name: "proc".to_string(),
                integrand_name: "itg".to_string(),
            }],
            None,
            None,
            None,
            vec![None],
        );
        state.all_integrals = vec![accumulator];
        let view_options = IntegrationStatusViewOptions {
            show_statistics: false,
            show_top_discrete_grid: false,
            show_discrete_contributions_sum: false,
            ..default_view_options()
        };
        let rendered = render_status_update_tabled(
            &build_status_update(
                StatusUpdateBuildRequest::new(
                    IntegrationStatusKind::Iteration,
                    &state,
                    &[None],
                    &view_options,
                )
                .with_timing(
                    1,
                    Duration::from_secs(0),
                    Duration::from_secs(0),
                    0,
                    0,
                    0,
                ),
            ),
            &default_tabled_options(),
        );

        assert!(rendered.contains("Maximum weight details"), "{rendered}");
        assert!(rendered.contains("Integrand"), "{rendered}");
        assert!(rendered.contains("Max eval"), "{rendered}");
        assert!(rendered.contains("Max eval coordinates"), "{rendered}");
        assert!(rendered.contains("proc@itg"), "{rendered}");
        assert!(rendered.contains("re [+]"), "{rendered}");
        assert!(rendered.contains("idx: 0, idx: 0, xs: [ "), "{rendered}");
        assert!(rendered.contains("e-01"), "{rendered}");
        assert!(rendered.contains("e-03 ]"), "{rendered}");
        assert!(!rendered.contains(", 0.532"), "{rendered}");
        assert!(!rendered.contains("( graph:"), "{rendered}");
    }

    #[test]
    fn iteration_status_block_uses_compact_header_and_optional_statistics() {
        let state = make_integration_state();
        let view_options = IntegrationStatusViewOptions {
            show_statistics: false,
            show_max_weight_details: false,
            ..default_view_options()
        };
        let rendered = render_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Iteration,
                &state,
                &[Some(Complex::new(F(1.0e-4), F(2.0e-5))), None],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(1),
                Duration::from_secs(1),
                100_000,
                100_000,
                100_000,
            ),
        );

        assert!(
            rendered.contains("Iteration #   1 ( completed )"),
            "{rendered}"
        );
        assert!(
            rendered.contains("# samples per iteration = 100.00K # samples total = 100.00K"),
            "{rendered}"
        );
        assert!(rendered.contains("/sample/core"), "{rendered}");
        assert!(rendered.contains("Contribution"), "{rendered}");
        assert!(rendered.contains("proc_a@itg_a"), "{rendered}");
        assert!(rendered.contains("proc_b@itg_b"), "{rendered}");
        assert!(rendered.contains("+7.5(9.8)e-5"), "{rendered}");
        assert!(rendered.contains("Δ = 0.255σ"), "{rendered}");
        assert!(!rendered.contains("Integration statistics"), "{rendered}");
        let lines = rendered.lines().collect::<Vec<_>>();
        assert_eq!(lines[1].matches('│').count(), 2, "{rendered}");
        assert!(lines[2].contains('┬'), "{rendered}");
        assert!(
            lines.last().is_some_and(|line| line.contains('┴')),
            "{rendered}"
        );
    }

    #[test]
    fn iteration_status_block_omits_delta_columns_when_no_target_is_provided() {
        let state = make_integration_state();
        let view_options = IntegrationStatusViewOptions {
            show_statistics: true,
            show_max_weight_details: false,
            ..default_view_options()
        };
        let rendered = render_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Iteration,
                &state,
                &[None, None],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(1),
                Duration::from_secs(1),
                100_000,
                100_000,
                100_000,
            ),
        );

        assert!(!rendered.contains("Δ [σ]"), "{rendered}");
        assert!(!rendered.contains("Δ ="), "{rendered}");
        assert!(rendered.contains("Integration statistics"), "{rendered}");
        assert!(rendered.contains("mwi"), "{rendered}");
    }

    #[test]
    fn live_iteration_status_block_uses_progress_header() {
        let mut state = make_integration_state();
        state.iter = 2;
        state.num_points = 125_000;
        let view_options = IntegrationStatusViewOptions {
            show_statistics: false,
            show_max_weight_details: false,
            ..default_view_options()
        };
        let rendered = render_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Live,
                &state,
                &[Some(Complex::new(F(1.0e-4), F(2.0e-5))), None],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(1),
                Duration::from_secs(1),
                25_000,
                125_000,
                125_000,
            )
            .with_live_progress(Some(status_update::LiveIterationProgress {
                completed_points: 25_000,
                target_points: 100_000,
            })),
        );

        assert!(
            rendered.contains("Iteration #   2 ( running   )"),
            "{rendered}"
        );
        assert!(
            rendered.contains("Iteration progress 25.00K/100.00K"),
            "{rendered}"
        );
        assert!(rendered.contains("25.0%"), "{rendered}");
        assert!(rendered.contains("# samples total = 125.00K"), "{rendered}");
    }

    #[test]
    fn iteration_status_block_shows_discrete_contributions_and_discrete_max_weights() {
        let state = make_discrete_integration_state();
        let view_options = IntegrationStatusViewOptions {
            show_statistics: false,
            show_max_weight_details: true,
            show_top_discrete_grid: true,
            show_discrete_contributions_sum: true,
            contribution_sort: ContributionSortMode::Index,
            show_max_weight_info_for_discrete_bins: true,
            ..default_view_options()
        };
        let rendered = render_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Iteration,
                &state,
                &[Some(Complex::new(F(1.0e-4), F(2.0e-5))), None],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(2),
                Duration::from_secs(2),
                110_000,
                210_000,
                210_000,
            ),
        );

        assert!(rendered.contains("Sum"), "{rendered}");
        assert!(rendered.contains("Contribution (idx=graph)"), "{rendered}");
        assert!(
            !rendered.contains("│             (idx=graph)             │"),
            "{rendered}"
        );
        assert!(rendered.contains("#0: GL0"), "{rendered}");
        assert!(rendered.contains("75.0%"), "{rendered}");
        assert!(rendered.contains("25.0%"), "{rendered}");
        assert!(
            rendered.contains("Maximum weight details by discrete bin"),
            "{rendered}"
        );
        assert!(
            rendered.contains("xs: [ 2.5000000000000000e-01 ]"),
            "{rendered}"
        );
    }

    #[test]
    fn grouped_graph_descriptions_include_all_group_members() {
        let description =
            graph_group_description(["GL0".to_string(), "GL1".to_string(), "GL2".to_string()]);

        assert_eq!(description, "[GL0,GL1,GL2]");
    }

    #[test]
    fn iteration_status_block_hides_spanned_metadata_header_separators() {
        let state = make_discrete_integration_state();
        let view_options = IntegrationStatusViewOptions {
            show_statistics: false,
            show_max_weight_details: false,
            show_top_discrete_grid: true,
            show_discrete_contributions_sum: true,
            contribution_sort: ContributionSortMode::Index,
            show_max_weight_info_for_discrete_bins: false,
            ..default_view_options()
        };
        let rendered = render_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Iteration,
                &state,
                &[Some(Complex::new(F(1.0e-4), F(2.0e-5))), None],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(2),
                Duration::from_secs(2),
                110_000,
                210_000,
                210_000,
            ),
        );

        let header_line = rendered
            .lines()
            .find(|line| line.contains("χ²/dof") && line.contains("mwi"))
            .expect("expected spanned metadata header line");

        let chi_to_mwi =
            &header_line[header_line.find("χ²/dof").unwrap()..header_line.find("mwi").unwrap()];
        assert!(!chi_to_mwi.contains('│'), "{rendered}");

        let delta_sigma_to_percent =
            &header_line[header_line.find("Δ [σ]").unwrap()..header_line.find("Δ [%]").unwrap()];
        assert!(!delta_sigma_to_percent.contains('│'), "{rendered}");
    }

    #[test]
    fn integration_result_always_contains_first_non_trivial_discrete_breakdown() {
        let state = make_discrete_integration_state();
        let result =
            build_integration_result(&state, &[Some(Complex::new(F(1.0e-4), F(2.0e-5))), None]);

        let slot = result
            .slot("proc_a@itg_a")
            .expect("discrete slot must be present");
        let breakdown = slot
            .grid_breakdown
            .re
            .as_ref()
            .expect("discrete breakdown must be persisted");

        assert_eq!(breakdown.axis_label, "graph");
        assert!(breakdown.fixed_coordinates.is_empty());
        assert_eq!(breakdown.entries.len(), 2);
        assert_eq!(breakdown.entries[0].bin_index, 0);
        assert_eq!(breakdown.entries[0].bin_label.as_deref(), Some("GL0"));
        assert_eq!(breakdown.entries[0].pdf, F(0.75));
        assert_eq!(breakdown.entries[0].processed_samples, 150);
        assert_eq!(breakdown.entries[1].bin_index, 1);
        assert_eq!(breakdown.entries[1].bin_label.as_deref(), Some("GL1"));
        assert_eq!(breakdown.entries[1].pdf, F(0.25));
        assert_eq!(breakdown.entries[1].processed_samples, 50);
    }

    #[test]
    fn orientation_descriptions_use_colored_signs() {
        let mut state = make_discrete_integration_state();
        state.first_non_trivial_discrete_label = Some("orientation".to_string());
        state.first_non_trivial_discrete_bin_descriptions =
            Some(vec!["+-0".to_string(), "0++".to_string()]);
        for metadata in &mut state.slot_first_non_trivial_discrete_breakdown_metadata {
            if let Some(metadata) = metadata.as_mut() {
                metadata.axis_label = "orientation".to_string();
                metadata.bin_labels = vec!["+-0".to_string(), "0++".to_string()];
            }
        }

        control::set_override(true);
        let expected_plus = "+".green().bold().to_string();
        let expected_minus = "-".red().bold().to_string();

        let view_options = IntegrationStatusViewOptions {
            show_statistics: false,
            show_max_weight_details: false,
            show_top_discrete_grid: true,
            show_discrete_contributions_sum: false,
            contribution_sort: ContributionSortMode::Index,
            show_max_weight_info_for_discrete_bins: false,
            ..default_view_options()
        };
        let ansi = render_tabled::render_status_update(
            &build_status_update(
                StatusUpdateBuildRequest::new(
                    IntegrationStatusKind::Iteration,
                    &state,
                    &[None, None],
                    &view_options,
                )
                .with_timing(
                    1,
                    Duration::from_secs(0),
                    Duration::from_secs(0),
                    0,
                    0,
                    0,
                ),
            ),
            &tabled_options_for_view(&view_options),
        );
        control::set_override(false);

        assert!(ansi.contains(&expected_plus), "{ansi}");
        assert!(ansi.contains(&expected_minus), "{ansi}");
        assert!(ansi.contains("0"), "{ansi}");
    }

    #[test]
    fn mismatched_bin_descriptions_fall_back_to_indices_with_warning() {
        let (descriptions, warning) = coalesce_first_non_trivial_discrete_bin_descriptions(
            "graph",
            &[
                (
                    "proc_a@itg_a".to_string(),
                    vec!["GL0".to_string(), "GL1".to_string()],
                ),
                (
                    "proc_b@itg_b".to_string(),
                    vec!["GX0".to_string(), "GX1".to_string()],
                ),
            ],
        );

        assert!(descriptions.is_none());
        let warning = warning.expect("mismatch should produce a warning");
        assert!(warning.contains("graph"), "{warning}");
        assert!(warning.contains("proc_a@itg_a"), "{warning}");
        assert!(warning.contains("proc_b@itg_b"), "{warning}");
    }

    #[test]
    fn format_max_eval_sample_keeps_full_discrete_coordinates() {
        let axis_labels = vec!["graph".to_string(), "LMB channel".to_string()];
        let full_sample = Sample::Discrete(
            F(1.0),
            0,
            Some(Box::new(Sample::Discrete(
                F(1.0),
                1,
                Some(Box::new(Sample::Continuous(F(1.0), vec![F(0.25)]))),
            ))),
        );
        let nested_sample = Sample::Discrete(
            F(1.0),
            1,
            Some(Box::new(Sample::Continuous(F(1.0), vec![F(0.75)]))),
        );

        assert_eq!(
            display::format_max_eval_sample(&full_sample, &axis_labels, &[]),
            "graph: 0, LMB channel: 1, xs: [ 2.5000000000000000e-01 ]"
        );
        assert_eq!(
            display::format_max_eval_sample(&nested_sample, &axis_labels, &[0]),
            "graph: 0, LMB channel: 1, xs: [ 7.5000000000000000e-01 ]"
        );
    }

    #[test]
    fn format_max_eval_sample_wraps_long_coordinate_lists() {
        let axis_labels = vec!["graph".to_string(), "orientation".to_string()];
        let sample = Sample::Discrete(
            F(1.0),
            0,
            Some(Box::new(Sample::Discrete(
                F(1.0),
                9,
                Some(Box::new(Sample::Continuous(
                    F(1.0),
                    vec![
                        F(0.125),
                        F(0.25),
                        F(0.375),
                        F(0.5),
                        F(0.625),
                        F(0.75),
                        F(0.875),
                    ],
                ))),
            ))),
        );

        assert_eq!(
            display::format_max_eval_sample(&sample, &axis_labels, &[]),
            "graph: 0, orientation: 9, xs: [\n1.2500000000000000e-01 2.5000000000000000e-01 3.7500000000000000e-01\n5.0000000000000000e-01 6.2500000000000000e-01 7.5000000000000000e-01\n8.7500000000000000e-01 ]"
        );
    }

    #[test]
    fn discrete_bin_prefixes_use_fixed_width_within_each_range() {
        assert_eq!(status_update::format_discrete_bin_prefix(2, 99), "#2: ");
        assert_eq!(status_update::format_discrete_bin_prefix(12, 99), "#12:");
        assert_eq!(status_update::format_discrete_bin_prefix(2, 999), "#2:  ");
        assert_eq!(status_update::format_discrete_bin_prefix(12, 999), "#12: ");
        assert_eq!(status_update::format_discrete_bin_prefix(123, 999), "#123:");
        assert_eq!(
            status_update::format_discrete_bin_prefix(2, 9_999),
            "#2:   "
        );
        assert_eq!(
            status_update::format_discrete_bin_prefix(1234, 9_999),
            "#1234:"
        );
    }

    #[test]
    fn final_summary_omits_discrete_max_weight_block_when_disabled() {
        let state = make_discrete_integration_state();
        let rendered = render_saved_integration_summary(
            &state,
            &[None, None],
            &IntegrationStatusViewOptions {
                show_statistics: true,
                show_max_weight_details: true,
                show_top_discrete_grid: false,
                show_discrete_contributions_sum: false,
                contribution_sort: ContributionSortMode::Index,
                show_max_weight_info_for_discrete_bins: false,
                ..default_view_options()
            },
            &default_tabled_options(),
        );

        assert!(rendered.contains("Maximum weight details"), "{rendered}");
        assert!(
            !rendered.contains("Maximum weight details by discrete bin"),
            "{rendered}"
        );
        assert!(
            !rendered.contains("# samples per iteration = 0"),
            "{rendered}"
        );
    }

    #[test]
    fn ratatui_overview_shows_eta_and_all_slot_metrics() {
        let state = make_integration_state();
        let view_options = IntegrationStatusViewOptions {
            show_statistics: false,
            show_max_weight_details: false,
            ..default_view_options()
        };
        let rendered = render_ratatui_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Live,
                &state,
                &[Some(Complex::new(F(1.0e-4), F(2.0e-5))), None],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(15),
                Duration::from_secs(10),
                25_000,
                125_000,
                125_000,
            )
            .with_live_progress(Some(status_update::LiveIterationProgress {
                completed_points: 25_000,
                target_points: 100_000,
            })),
            |_| {},
        );

        assert!(rendered.contains("Iteration progress"), "{rendered}");
        assert!(rendered.contains("ETA 30s"), "{rendered}");
        assert!(rendered.contains("25.00K / 100.00K (25.0%)"), "{rendered}");
        assert!(rendered.contains("# samples total 125.00K"), "{rendered}");
        assert!(rendered.contains("#samples/s"), "{rendered}");
        assert!(rendered.contains("Integrands"), "{rendered}");
        assert!(rendered.contains("Focused integrand"), "{rendered}");
        assert!(rendered.contains("Results summary"), "{rendered}");
        assert!(rendered.contains("% err"), "{rendered}");
        assert!(rendered.contains("chi^2"), "{rendered}");
        assert!(rendered.contains("m.w.i"), "{rendered}");
        assert!(rendered.contains("Timing composition"), "{rendered}");
        assert!(rendered.contains("Precision mix"), "{rendered}");
    }

    #[test]
    fn ratatui_overview_chart_phase_toggle_updates_convergence_label() {
        let state = make_integration_state();
        let view_options = IntegrationStatusViewOptions {
            show_statistics: false,
            show_max_weight_details: false,
            ..default_view_options()
        };
        let rendered = render_ratatui_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Live,
                &state,
                &[Some(Complex::new(F(1.0e-4), F(2.0e-5))), None],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(15),
                Duration::from_secs(10),
                25_000,
                125_000,
                125_000,
            )
            .with_live_progress(Some(status_update::LiveIterationProgress {
                completed_points: 25_000,
                target_points: 100_000,
            })),
            |dashboard| dashboard.toggle_chart_component(),
        );

        assert!(
            rendered.contains("Convergence : imag (not selected for training)"),
            "{rendered}"
        );
    }

    #[test]
    fn ratatui_integrand_targets_are_trimmed_for_display() {
        let state = make_integration_state();
        let view_options = IntegrationStatusViewOptions {
            show_statistics: false,
            show_max_weight_details: false,
            ..default_view_options()
        };
        let rendered = render_ratatui_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Live,
                &state,
                &[
                    Some(Complex::new(
                        F(7.600000000000001e-6),
                        F(7.430000000000004e-5),
                    )),
                    Some(Complex::new(F(6.629999999999999e-5), F(0.0))),
                ],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(15),
                Duration::from_secs(10),
                25_000,
                125_000,
                125_000,
            )
            .with_live_progress(Some(status_update::LiveIterationProgress {
                completed_points: 25_000,
                target_points: 100_000,
            })),
            |_| {},
        );

        assert!(rendered.contains("trgt"), "{rendered}");
        assert!(rendered.contains("+7.6e-6"), "{rendered}");
        assert!(rendered.contains("+7.43e-5"), "{rendered}");
        assert!(rendered.contains("+6.63e-5"), "{rendered}");
        assert!(rendered.contains("+0e0"), "{rendered}");
    }

    #[test]
    fn ratatui_overview_shows_eta_to_target_when_configured() {
        let state = make_integration_state();
        let view_options = IntegrationStatusViewOptions {
            target_relative_accuracy: Some(0.05),
            show_statistics: false,
            show_max_weight_details: false,
            ..default_view_options()
        };
        let rendered = render_ratatui_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Live,
                &state,
                &[Some(Complex::new(F(1.0e-4), F(0.0))), None],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(25),
                Duration::from_secs(10),
                25_000,
                125_000,
                125_000,
            )
            .with_live_progress(Some(status_update::LiveIterationProgress {
                completed_points: 25_000,
                target_points: 100_000,
            })),
            |_| {},
        );

        assert!(rendered.contains("ETA to target"), "{rendered}");
        assert!(rendered.contains("(% err <= 5%)"), "{rendered}");
    }

    #[test]
    fn eta_to_target_specification_formats_relative_and_absolute_targets() {
        let state = make_integration_state();
        let relative_view_options = IntegrationStatusViewOptions {
            target_relative_accuracy: Some(0.05),
            ..default_view_options()
        };
        let update = build_status_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Live,
                &state,
                &[Some(Complex::new(F(1.0e-4), F(0.0))), None],
                &relative_view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(25),
                Duration::from_secs(10),
                25_000,
                125_000,
                125_000,
            )
            .with_live_progress(Some(status_update::LiveIterationProgress {
                completed_points: 25_000,
                target_points: 100_000,
            })),
        );
        assert_eq!(
            update.meta.eta_to_target_specification(),
            Some("% err <= 5%")
        );

        let relative_scientific_view_options = IntegrationStatusViewOptions {
            target_relative_accuracy: Some(1.0e-6),
            ..default_view_options()
        };
        let update = build_status_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Live,
                &state,
                &[Some(Complex::new(F(1.0e-4), F(0.0))), None],
                &relative_scientific_view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(25),
                Duration::from_secs(10),
                25_000,
                125_000,
                125_000,
            )
            .with_live_progress(Some(status_update::LiveIterationProgress {
                completed_points: 25_000,
                target_points: 100_000,
            })),
        );
        assert_eq!(
            update.meta.eta_to_target_specification(),
            Some("% err <= 1e-4%")
        );

        let absolute_view_options = IntegrationStatusViewOptions {
            target_absolute_accuracy: Some(1.0e-6),
            ..default_view_options()
        };
        let update = build_status_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Live,
                &state,
                &[Some(Complex::new(F(1.0e-4), F(0.0))), None],
                &absolute_view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(25),
                Duration::from_secs(10),
                25_000,
                125_000,
                125_000,
            )
            .with_live_progress(Some(status_update::LiveIterationProgress {
                completed_points: 25_000,
                target_points: 100_000,
            })),
        );
        assert_eq!(
            update.meta.eta_to_target_specification(),
            Some("err <= 1e-6")
        );

        let combined_view_options = IntegrationStatusViewOptions {
            target_relative_accuracy: Some(0.05),
            target_absolute_accuracy: Some(1.0e-6),
            ..default_view_options()
        };
        let update = build_status_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Live,
                &state,
                &[Some(Complex::new(F(1.0e-4), F(0.0))), None],
                &combined_view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(25),
                Duration::from_secs(10),
                25_000,
                125_000,
                125_000,
            )
            .with_live_progress(Some(status_update::LiveIterationProgress {
                completed_points: 25_000,
                target_points: 100_000,
            })),
        );
        assert_eq!(
            update.meta.eta_to_target_specification(),
            Some("% err <= 5% or err <= 1e-6")
        );
    }

    #[test]
    fn target_accuracy_status_uses_either_absolute_or_relative_condition() {
        let state = make_integration_state();

        let reached_by_absolute = status_update::evaluate_target_accuracy(
            &state,
            100_000,
            Duration::from_secs(10),
            &[None, None],
            IntegrationStatusPhaseDisplay::Real,
            None,
            Some(1.0e-4),
        );
        assert!(reached_by_absolute.is_reached());

        let reached_by_relative = status_update::evaluate_target_accuracy(
            &state,
            100_000,
            Duration::from_secs(10),
            &[Some(Complex::new(F(1.0e-3), F(0.0))), None],
            IntegrationStatusPhaseDisplay::Real,
            Some(0.1),
            None,
        );
        assert!(reached_by_relative.is_reached());
    }

    #[test]
    fn target_accuracy_status_treats_zero_relative_reference_as_inactive() {
        let state = make_integration_state();
        let status = status_update::evaluate_target_accuracy(
            &state,
            100_000,
            Duration::from_secs(10),
            &[Some(Complex::new(F(0.0), F(0.0))), None],
            IntegrationStatusPhaseDisplay::Real,
            Some(0.05),
            None,
        );

        assert!(!status.relative_reached);
        assert!(!status.absolute_reached);
        assert_eq!(status.eta_to_target, None);
    }

    #[test]
    fn ratatui_recent_history_window_tracks_current_iteration_points() {
        let mut dashboard = RatatuiDashboardState::new();
        let mut state = make_integration_state();
        let view_options = IntegrationStatusViewOptions {
            show_statistics: false,
            show_max_weight_details: false,
            ..default_view_options()
        };

        state.iter = 1;
        dashboard.update(build_status_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Live,
                &state,
                &[None, None],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(10),
                Duration::from_secs(10),
                500_000,
                1_000_000,
                1_000_000,
            )
            .with_live_progress(Some(status_update::LiveIterationProgress {
                completed_points: 500_000,
                target_points: 1_000_000,
            })),
        ));

        state.iter = 2;
        dashboard.update(build_status_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Live,
                &state,
                &[None, None],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(20),
                Duration::from_secs(5),
                0,
                2_000_000,
                2_000_000,
            )
            .with_live_progress(Some(status_update::LiveIterationProgress {
                completed_points: 0,
                target_points: 1_000_000,
            })),
        ));
        dashboard.update(build_status_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Live,
                &state,
                &[None, None],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(25),
                Duration::from_secs(10),
                200_000,
                2_200_000,
                2_200_000,
            )
            .with_live_progress(Some(status_update::LiveIterationProgress {
                completed_points: 200_000,
                target_points: 1_000_000,
            })),
        ));

        dashboard.toggle_chart_history_window();
        for _ in 0..5 {
            dashboard.narrow_chart_history_window();
        }

        assert_eq!(
            dashboard.visible_history_sample_bounds(),
            Some((2_000_000, 2_200_000))
        );
    }

    #[test]
    fn ratatui_full_history_preserves_origin_after_compaction() {
        let mut dashboard = RatatuiDashboardState::new();
        let mut state = make_integration_state();
        let view_options = IntegrationStatusViewOptions {
            show_statistics: false,
            show_max_weight_details: false,
            ..default_view_options()
        };

        for index in 0..4_300 {
            state.iter = 1 + index / 40;
            let total_points = index * 1_000;
            dashboard.update(build_status_update(
                StatusUpdateBuildRequest::new(
                    IntegrationStatusKind::Live,
                    &state,
                    &[None, None],
                    &view_options,
                )
                .with_timing(
                    4,
                    Duration::from_secs(index as u64),
                    Duration::from_secs(10),
                    total_points % 1_000_000,
                    total_points,
                    total_points,
                )
                .with_live_progress(Some(status_update::LiveIterationProgress {
                    completed_points: total_points % 1_000_000,
                    target_points: 1_000_000,
                })),
            ));
        }

        dashboard.toggle_chart_history_window();
        dashboard.toggle_chart_history_window();

        assert_eq!(
            dashboard.visible_history_sample_bounds(),
            Some((0, 4_299_000))
        );
    }

    #[test]
    fn ratatui_convergence_title_shows_active_y_span() {
        let state = make_integration_state();
        let view_options = IntegrationStatusViewOptions {
            show_statistics: false,
            show_max_weight_details: false,
            ..default_view_options()
        };
        let rendered = render_ratatui_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Live,
                &state,
                &[None, None],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(15),
                Duration::from_secs(10),
                25_000,
                125_000,
                125_000,
            )
            .with_live_progress(Some(status_update::LiveIterationProgress {
                completed_points: 25_000,
                target_points: 100_000,
            })),
            |dashboard| dashboard.widen_chart_y_sigma_span(),
        );

        assert!(rendered.contains("y ±5σ"), "{rendered}");
    }

    #[test]
    fn ratatui_discrete_tab_renders_selected_bin_detail() {
        let state = make_discrete_integration_state();
        let view_options = IntegrationStatusViewOptions {
            show_statistics: false,
            show_max_weight_details: false,
            show_top_discrete_grid: false,
            show_discrete_contributions_sum: false,
            show_max_weight_info_for_discrete_bins: false,
            ..default_view_options()
        };
        let rendered = render_ratatui_update(
            StatusUpdateBuildRequest::new(
                IntegrationStatusKind::Iteration,
                &state,
                &[Some(Complex::new(F(1.0e-4), F(2.0e-5))), None],
                &view_options,
            )
            .with_timing(
                4,
                Duration::from_secs(12),
                Duration::from_secs(12),
                110_000,
                210_000,
                210_000,
            ),
            |dashboard| dashboard.select_tab(1),
        );

        assert!(rendered.contains("Selected bin"), "{rendered}");
        assert!(
            rendered.contains("Discrete bins for focused integrand"),
            "{rendered}"
        );
        assert!(rendered.contains("sample %"), "{rendered}");
        assert!(rendered.contains("GL0"), "{rendered}");
        assert!(rendered.contains("GL1"), "{rendered}");
        assert!(rendered.contains("pdf"), "{rendered}");
        assert!(rendered.contains("Per integrand details"), "{rendered}");
        assert!(rendered.contains("# samples"), "{rendered}");
    }
}

#[test]
fn test_threading() {
    use symbolica::numerical_integration::ContinuousGrid;

    fn test_fn(x: f64) -> f64 {
        x * x * (x * 6.).sin()
    }

    let mut acc_1 = StatisticsAccumulator::<f64>::new();
    let mut acc_2 = StatisticsAccumulator::<f64>::new();

    let samples_per_sample = 4;
    let samples_per_iter = 100000;
    let n_iter = 10;

    let mut rng = MonteCarloRng::new(42, 0);

    let mut grid = Grid::<f64>::Continuous(ContinuousGrid::new(1, 64, 100, None, false));

    for _i_iter in 0..n_iter {
        let mut multiplice_accs = vec![StatisticsAccumulator::<f64>::new(); samples_per_sample];
        for _i_sample in 0..samples_per_iter {
            let n_samples = (0..samples_per_sample)
                .map(|_| {
                    let mut sample = Sample::new();
                    grid.sample(&mut rng, &mut sample);
                    sample
                })
                .collect_vec();

            let n_evals = n_samples
                .iter()
                .map(|s| match s {
                    Sample::Continuous(_, xs) => test_fn(xs[0]),
                    _ => unreachable!(),
                })
                .collect_vec();

            for (i_eval, (sample, eval)) in n_samples.iter().zip(&n_evals).enumerate() {
                acc_1.add_sample(eval * sample.get_weight(), Some(sample));
                multiplice_accs[i_eval].add_sample(eval * sample.get_weight(), Some(sample));
            }
        }

        for acc in multiplice_accs {
            acc_2.merge_samples(&mut acc.clone());
        }
        acc_1.update_iter(false);
        acc_2.update_iter(false);
    }

    println!("acc1: {:?}", acc_1);
    println!("acc2: {:?}", acc_2);
}
