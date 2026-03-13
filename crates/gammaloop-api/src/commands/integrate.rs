use std::{
    collections::{BTreeMap, HashSet},
    fmt, fs,
    io::{self, IsTerminal, Write},
    ops::Deref,
    path::PathBuf,
};

use clap::{ArgAction, Args, ValueEnum};
use crossterm::{
    cursor::{MoveDown, MoveToColumn, MoveUp},
    queue,
    terminal::{self, Clear, ClearType},
};
use gammalooprs::utils::serde_utils::SmartSerde;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use unicode_width::UnicodeWidthChar;

use color_eyre::{eyre::eyre, Result};
use colored::Colorize;
use gammalooprs::{
    integrands::{HasIntegrand, Integrand},
    integrate::{
        build_integration_result, emit_integration_status_via_tracing, havana_integrate,
        latest_observable_resume_state_path, print_integral_result,
        render_saved_integration_summary, slot_workspace_path, workspace_manifest_path,
        workspace_state_path, ContributionSortMode, IntegrationState, IntegrationStatusKind,
        IntegrationStatusPhaseDisplay, IntegrationStatusRenderOptions,
        IntegrationWorkspaceManifest, IterationBatchingSettings, SamplingCorrelationMode, SlotMeta,
        WorkspaceSnapshotControl,
    },
    model::{Model, SerializableInputParamCard},
    observables::ObservableSnapshotBundle,
    settings::{
        runtime::{IntegratedPhase, IntegrationResult},
        RuntimeSettings,
    },
    utils::F,
};
use itertools::Itertools;
use symbolica::numerical_integration::Grid;
use tracing::{info, warn};

use crate::{
    completion::CompletionArgExt,
    state::{ProcessRef, State},
    CLISettings,
};

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, unsendable, name = "IntegrationSettings")
)]
#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Integrate {
    /// Process reference: #<id>, name:<name>, or <id>/<name>; repeat with -i to select multiple integrands
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Any)
    )]
    pub process: Vec<ProcessRef>,

    /// The integrand name to integrate; repeat with -p to select multiple integrands
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Any)
    )]
    pub integrand_name: Vec<String>,

    /// Number of cores to parallelize over
    #[arg(short = 'c', long)]
    pub n_cores: Option<usize>,

    /// The path to run the integrationg within
    #[arg(short = 'w', long, value_hint = clap::ValueHint::DirPath)]
    pub workspace_path: Option<PathBuf>,

    /// Specify a target as either `re im` for a single integrand or `process@integrand=re,im`
    #[arg(
        short = 't',
        long,
        num_args = 1..=2,
        action = ArgAction::Append,
        allow_negative_numbers = true,
        allow_hyphen_values = true,
        completion_selected_integrand_target()
    )]
    pub target: Vec<String>,

    /// Whether to restart the integration from scratch, or continue from a previous run if possible
    #[arg(short = 'r', long)]
    pub restart: bool,

    /// Use one independent sampling/training grid per selected integrand instead of a shared grid
    #[arg(long = "uncorrelated")]
    pub uncorrelated: bool,

    /// Show the overall max-weight details table
    #[arg(long = "show-max-weight-info", default_value_t = true)]
    pub show_max_weight_info: bool,

    /// Hide integration statistics during intermediate iteration updates
    #[arg(long = "no-show-integration-statistics")]
    pub no_show_integration_statistics: bool,

    /// Control which complex phase is shown in integration summary tables
    #[arg(long = "show-phase", default_value = "both")]
    pub show_phase: ShowPhaseOption,

    /// Show monitoring rows for the first non-trivial discrete dimension
    #[arg(long = "show-top-discrete-grid")]
    pub show_top_discrete_grid: bool,

    /// Show the summed discrete contributions block independently of per-bin rows
    #[arg(long = "show-discrete-contributions-sum")]
    pub show_discrete_contributions_sum: bool,

    /// Sort discrete contribution rows by index, absolute integral, or error
    #[arg(long = "sort-contributions", default_value = "error")]
    pub sort_contributions: ContributionSortOption,

    /// Show per-bin max-weight information for the first non-trivial discrete dimension
    #[arg(
        long = "show-max-weight-info-for-discrete-bins",
        default_value_t = false
    )]
    pub show_max_weight_info_for_discrete_bins: bool,

    /// Only render the saved integration summary, without performing integration work
    #[arg(long = "show-summary-only")]
    pub show_summary_only: bool,

    /// Disable streaming of end-of-iteration integration summaries in interactive terminals
    #[arg(long = "no-stream-iterations")]
    pub no_stream_iterations: bool,

    /// Disable live streaming of in-iteration integration status updates
    #[arg(long = "no-stream-updates")]
    pub no_stream_updates: bool,

    /// Evaluate each iteration in per-core batches of this size
    #[arg(long = "batch-size")]
    pub batch_size: Option<usize>,

    /// Target per-core batch wall time in seconds when batch size is not specified
    #[arg(long = "batch-timing", default_value_t = 5.0)]
    pub batch_timing: f64,

    /// Minimum time between streamed in-iteration status updates
    #[arg(long = "min-time-between-status-updates", default_value_t = 0.0)]
    pub min_time_between_status_updates: f64,

    /// Maximum width used when normalizing integration status tables
    #[arg(long = "max-table-width", default_value_t = 250)]
    pub max_table_width: usize,

    /// Keep per-iteration result and observable snapshot files in addition to the latest snapshots
    #[arg(long = "write-results-for-each-iteration")]
    pub write_results_for_each_iteration: bool,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct IntegrationOutput {
    pub result: IntegrationResult,
    pub observables: BTreeMap<String, ObservableSnapshotBundle>,
    pub workspace_path: PathBuf,
}

impl IntegrationOutput {
    pub fn slot_observables(&self, key: &str) -> Option<&ObservableSnapshotBundle> {
        self.observables.get(key)
    }

    pub fn single_slot_observables(&self) -> Option<&ObservableSnapshotBundle> {
        (self.observables.len() == 1)
            .then(|| self.observables.values().next())
            .flatten()
    }
}

impl Deref for IntegrationOutput {
    type Target = IntegrationResult;

    fn deref(&self) -> &Self::Target {
        &self.result
    }
}

impl fmt::Display for IntegrationOutput {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.result.fmt(f)
    }
}

impl Default for Integrate {
    fn default() -> Self {
        Self {
            process: Vec::new(),
            integrand_name: Vec::new(),
            n_cores: None,
            workspace_path: None,
            target: Vec::new(),
            restart: false,
            uncorrelated: false,
            show_max_weight_info: true,
            no_show_integration_statistics: false,
            show_phase: ShowPhaseOption::Both,
            show_top_discrete_grid: false,
            show_discrete_contributions_sum: false,
            sort_contributions: ContributionSortOption::Error,
            show_max_weight_info_for_discrete_bins: false,
            show_summary_only: false,
            no_stream_iterations: false,
            no_stream_updates: false,
            batch_size: None,
            batch_timing: 5.0,
            min_time_between_status_updates: 0.0,
            max_table_width: 250,
            write_results_for_each_iteration: false,
        }
    }
}

#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema, ValueEnum, Default,
)]
pub enum ShowPhaseOption {
    #[default]
    Both,
    Real,
    Imag,
    Selected,
}

impl ShowPhaseOption {
    fn resolve(self, selected_phase: IntegratedPhase) -> IntegrationStatusPhaseDisplay {
        match self {
            Self::Both => IntegrationStatusPhaseDisplay::Both,
            Self::Real => IntegrationStatusPhaseDisplay::Real,
            Self::Imag => IntegrationStatusPhaseDisplay::Imag,
            Self::Selected => match selected_phase {
                IntegratedPhase::Imag => IntegrationStatusPhaseDisplay::Imag,
                IntegratedPhase::Both => IntegrationStatusPhaseDisplay::Both,
                IntegratedPhase::Real => IntegrationStatusPhaseDisplay::Real,
            },
        }
    }
}

#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema, ValueEnum, Default,
)]
pub enum ContributionSortOption {
    Index,
    Integral,
    #[default]
    Error,
}

impl From<ContributionSortOption> for ContributionSortMode {
    fn from(value: ContributionSortOption) -> Self {
        match value {
            ContributionSortOption::Index => ContributionSortMode::Index,
            ContributionSortOption::Integral => ContributionSortMode::Integral,
            ContributionSortOption::Error => ContributionSortMode::Error,
        }
    }
}

struct StreamRenderer {
    stderr: io::Stderr,
    rendered_line_count: usize,
}

impl StreamRenderer {
    fn new() -> Self {
        Self {
            stderr: io::stderr(),
            rendered_line_count: 0,
        }
    }

    fn render(&mut self, block: &str) -> Result<()> {
        let prepared = prepare_stream_block(block, stream_terminal_width());
        if self.rendered_line_count > 0 {
            self.clear_rendered_block()?;
        }

        queue!(self.stderr, MoveToColumn(0))?;
        write!(self.stderr, "{prepared}")?;
        self.stderr.flush()?;
        self.rendered_line_count = prepared.lines().count().max(1);
        Ok(())
    }

    fn clear(&mut self) -> Result<()> {
        if self.rendered_line_count == 0 {
            return Ok(());
        }

        self.clear_rendered_block()?;
        self.rendered_line_count = 0;
        self.stderr.flush()?;
        Ok(())
    }

    fn clear_rendered_block(&mut self) -> Result<()> {
        queue!(self.stderr, MoveToColumn(0))?;
        if self.rendered_line_count > 1 {
            queue!(
                self.stderr,
                MoveUp((self.rendered_line_count.saturating_sub(1)) as u16)
            )?;
        }

        for line_index in 0..self.rendered_line_count {
            queue!(self.stderr, MoveToColumn(0), Clear(ClearType::CurrentLine))?;
            if line_index + 1 < self.rendered_line_count {
                queue!(self.stderr, MoveDown(1))?;
            }
        }
        queue!(self.stderr, MoveToColumn(0))?;
        if self.rendered_line_count > 1 {
            queue!(
                self.stderr,
                MoveUp((self.rendered_line_count.saturating_sub(1)) as u16)
            )?;
        }
        Ok(())
    }
}

impl Drop for StreamRenderer {
    fn drop(&mut self) {
        let _ = self.clear();
    }
}

fn stream_terminal_width() -> usize {
    terminal::size()
        .map(|(width, _)| width.saturating_sub(1).max(1) as usize)
        .unwrap_or(120)
}

fn prepare_stream_block(block: &str, max_width: usize) -> String {
    block
        .lines()
        .map(|line| truncate_ansi_line(line, max_width))
        .collect::<Vec<_>>()
        .join("\n")
}

fn truncate_ansi_line(line: &str, max_width: usize) -> String {
    if max_width == 0 {
        return String::new();
    }

    let mut truncated = String::new();
    let mut chars = line.chars().peekable();
    let mut visible_width = 0usize;
    let mut saw_escape = false;
    let mut was_truncated = false;

    while let Some(ch) = chars.next() {
        if ch == '\u{1b}' && chars.peek() == Some(&'[') {
            saw_escape = true;
            truncated.push(ch);
            truncated.push(chars.next().unwrap());
            for code in chars.by_ref() {
                truncated.push(code);
                if ('@'..='~').contains(&code) {
                    break;
                }
            }
            continue;
        }

        let char_width = UnicodeWidthChar::width(ch).unwrap_or(0);
        if visible_width + char_width > max_width {
            was_truncated = true;
            break;
        }

        truncated.push(ch);
        visible_width += char_width;
    }

    if was_truncated && saw_escape && !truncated.ends_with("\u{1b}[0m") {
        truncated.push_str("\u{1b}[0m");
    }

    truncated
}

#[derive(Debug, Clone)]
struct ResolvedIntegrandSlot {
    process_id: usize,
    slot_meta: SlotMeta,
}

fn slot_settings_path(workspace: &std::path::Path, slot_meta: &SlotMeta) -> PathBuf {
    slot_workspace_path(workspace, slot_meta).join("settings.toml")
}

fn read_existing_workspace_state(
    workspace_path: &std::path::Path,
) -> Result<(IntegrationWorkspaceManifest, IntegrationState)> {
    let manifest = IntegrationWorkspaceManifest::from_file(
        workspace_manifest_path(workspace_path),
        "integration manifest",
    )?;
    let state_bytes = fs::read(workspace_state_path(workspace_path)).map_err(|err| {
        eyre!(
            "Could not read integration workspace state from {}: {err}",
            workspace_path.display()
        )
    })?;
    let integration_state = bincode::decode_from_slice::<IntegrationState, _>(
        &state_bytes,
        bincode::config::standard(),
    )
    .map_err(|err| eyre!("Could not deserialize integration state: {err}"))?
    .0;
    Ok((manifest, integration_state))
}

fn grids_have_compatible_sample_shape(lhs: &Grid<F<f64>>, rhs: &Grid<F<f64>>) -> bool {
    match (lhs, rhs) {
        (Grid::Continuous(lhs), Grid::Continuous(rhs)) => {
            lhs.continuous_dimensions.len() == rhs.continuous_dimensions.len()
        }
        (Grid::Discrete(lhs), Grid::Discrete(rhs)) => {
            lhs.bins.len() == rhs.bins.len()
                && lhs
                    .bins
                    .iter()
                    .zip(rhs.bins.iter())
                    .all(
                        |(lhs_bin, rhs_bin)| match (&lhs_bin.sub_grid, &rhs_bin.sub_grid) {
                            (Some(lhs_sub_grid), Some(rhs_sub_grid)) => {
                                grids_have_compatible_sample_shape(lhs_sub_grid, rhs_sub_grid)
                            }
                            (None, None) => true,
                            _ => false,
                        },
                    )
        }
        (
            Grid::Uniform(lhs_discrete, lhs_continuous),
            Grid::Uniform(rhs_discrete, rhs_continuous),
        ) => {
            lhs_discrete == rhs_discrete
                && lhs_continuous.continuous_dimensions.len()
                    == rhs_continuous.continuous_dimensions.len()
        }
        _ => false,
    }
}

impl Integrate {
    pub fn from_slots<I, S>(slots: I) -> Self
    where
        I: IntoIterator<Item = (ProcessRef, S)>,
        S: Into<String>,
    {
        let mut integrate = Self::default();
        for (process, integrand_name) in slots {
            integrate.process.push(process);
            integrate.integrand_name.push(integrand_name.into());
        }
        integrate
    }

    pub fn with_single_target(mut self, target: Complex<F<f64>>) -> Self {
        self.target = vec![target.re.to_string(), target.im.to_string()];
        self
    }

    pub fn with_keyed_targets<I, S>(mut self, targets: I) -> Self
    where
        I: IntoIterator<Item = (S, Complex<F<f64>>)>,
        S: Into<String>,
    {
        self.target = targets
            .into_iter()
            .map(|(slot_key, target)| format!("{}={},{}", slot_key.into(), target.re, target.im))
            .collect();
        self
    }

    fn default_workspace_path(&self, global_cli_settings: &CLISettings) -> PathBuf {
        if global_cli_settings.session.read_only_state {
            let workspace_name = global_cli_settings
                .state
                .name
                .as_deref()
                .map(str::trim)
                .filter(|name| !name.is_empty())
                .map(|name| format!("{name}_integration_workspace"))
                .unwrap_or_else(|| "integration_workspace".to_string());
            PathBuf::from(".").join(workspace_name)
        } else {
            global_cli_settings
                .state
                .folder
                .join("integration_workspace")
        }
    }

    fn resolve_selected_slots(&self, state: &State) -> Result<Vec<ResolvedIntegrandSlot>> {
        let selections = if self.process.is_empty() && self.integrand_name.is_empty() {
            vec![state.find_integrand_ref(None, None)?]
        } else if self.process.len() == 1 && self.integrand_name.is_empty() {
            vec![state.find_integrand_ref(self.process.first(), None)?]
        } else if self.process.is_empty() && self.integrand_name.len() == 1 {
            vec![state.find_integrand_ref(None, self.integrand_name.first())?]
        } else if self.process.len() == self.integrand_name.len() {
            self.process
                .iter()
                .zip(self.integrand_name.iter())
                .map(|(process, integrand_name)| {
                    state.find_integrand_ref(Some(process), Some(integrand_name))
                })
                .collect::<Result<Vec<_>>>()?
        } else {
            return Err(eyre!(
                "The integrate command expects repeated `-p/-i` pairs. Received {} process selector(s) and {} integrand selector(s).",
                self.process.len(),
                self.integrand_name.len()
            ));
        };

        let mut seen = HashSet::new();
        let mut resolved = Vec::with_capacity(selections.len());
        for (process_id, integrand_name) in selections {
            let process_name = state.process_list.processes[process_id]
                .definition
                .folder_name
                .clone();
            let slot_meta = SlotMeta {
                process_name,
                integrand_name,
            };
            if !seen.insert(slot_meta.key()) {
                return Err(eyre!(
                    "The same integrand '{}' was selected more than once",
                    slot_meta.key()
                ));
            }
            resolved.push(ResolvedIntegrandSlot {
                process_id,
                slot_meta,
            });
        }

        Ok(resolved)
    }

    fn resolve_manifest_slots(
        &self,
        state: &State,
        manifest: &IntegrationWorkspaceManifest,
    ) -> Result<Vec<ResolvedIntegrandSlot>> {
        manifest
            .slots
            .iter()
            .map(|slot_meta| {
                let process_ref = ProcessRef::Name(slot_meta.process_name.clone());
                let integrand_name = slot_meta.integrand_name.clone();
                let (process_id, resolved_integrand_name) =
                    state.find_integrand_ref(Some(&process_ref), Some(&integrand_name))?;
                Ok(ResolvedIntegrandSlot {
                    process_id,
                    slot_meta: SlotMeta {
                        process_name: slot_meta.process_name.clone(),
                        integrand_name: resolved_integrand_name,
                    },
                })
            })
            .collect()
    }

    fn build_render_options(
        &self,
        settings: &RuntimeSettings,
        show_statistics: bool,
    ) -> IntegrationStatusRenderOptions {
        IntegrationStatusRenderOptions {
            phase_display: self
                .show_phase
                .resolve(settings.integrator.integrated_phase),
            show_statistics,
            show_max_weight_details: self.show_max_weight_info,
            show_top_discrete_grid: self.show_top_discrete_grid,
            show_discrete_contributions_sum: self.show_discrete_contributions_sum,
            contribution_sort: self.sort_contributions.into(),
            show_max_weight_info_for_discrete_bins: self.show_max_weight_info_for_discrete_bins,
            max_table_width: self.max_table_width,
        }
    }

    fn workspace_snapshot_control(&self) -> WorkspaceSnapshotControl {
        WorkspaceSnapshotControl {
            write_iteration_archives: self.write_results_for_each_iteration,
        }
    }

    fn sampling_correlation_mode(&self) -> SamplingCorrelationMode {
        if self.uncorrelated {
            SamplingCorrelationMode::Uncorrelated
        } else {
            SamplingCorrelationMode::Correlated
        }
    }

    fn build_batching_settings(&self, emit_live_status_updates: bool) -> IterationBatchingSettings {
        IterationBatchingSettings {
            batch_size: self.batch_size,
            batch_timing_seconds: self.batch_timing,
            min_time_between_status_updates_seconds: self.min_time_between_status_updates,
            emit_live_status_updates,
        }
    }

    fn parse_target_components(first: &str, second: Option<&str>) -> Result<Complex<F<f64>>> {
        if let Some(second) = second {
            return Ok(Complex::new(F(first.parse()?), F(second.parse()?)));
        }

        let (re, im) = first
            .split_once(',')
            .ok_or_else(|| eyre!("Targets must be provided as `re im` or `re,im`"))?;
        Ok(Complex::new(F(re.parse()?), F(im.parse()?)))
    }

    fn resolve_targets(
        &self,
        selected_slots: &[ResolvedIntegrandSlot],
    ) -> Result<Vec<Option<Complex<F<f64>>>>> {
        let mut resolved = vec![None; selected_slots.len()];
        if self.target.is_empty() {
            return Ok(resolved);
        }

        if selected_slots.len() == 1 && self.target.iter().all(|target| !target.contains('=')) {
            let target = match self.target.as_slice() {
                [single] => Self::parse_target_components(single, None)?,
                [re, im] => Self::parse_target_components(re, Some(im))?,
                _ => {
                    return Err(eyre!(
                        "A single-integrand target must be given as `--target re im` or `--target re,im`"
                    ));
                }
            };
            resolved[0] = Some(target);
            return Ok(resolved);
        }

        for target in &self.target {
            let (slot_key, values) = target.split_once('=').ok_or_else(|| {
                eyre!(
                    "Multi-integrand targets must use the keyed form `--target process@integrand=re,im`"
                )
            })?;
            let parsed = Self::parse_target_components(values, None)?;
            let slot_index = selected_slots
                .iter()
                .position(|slot| slot.slot_meta.key() == slot_key)
                .ok_or_else(|| eyre!("Unknown target slot key '{}'", slot_key))?;
            if resolved[slot_index].is_some() {
                return Err(eyre!(
                    "The target for '{}' was specified more than once",
                    slot_key
                ));
            }
            resolved[slot_index] = Some(parsed);
        }

        Ok(resolved)
    }

    fn load_or_prepare_workspace_state(
        &self,
        state: &mut State,
        selected_slots: &[ResolvedIntegrandSlot],
        current_effective_model_parameters: &[SerializableInputParamCard<F<f64>>],
        current_integrand_fingerprints: &[String],
        workspace_path: &std::path::Path,
        targets: &mut Vec<Option<Complex<F<f64>>>>,
    ) -> Result<Option<IntegrationState>> {
        let path_to_state = workspace_state_path(workspace_path);
        let manifest_path = workspace_manifest_path(workspace_path);
        match fs::read(&path_to_state) {
            Ok(state_bytes) => {
                let manifest: IntegrationWorkspaceManifest =
                    IntegrationWorkspaceManifest::from_file(
                        &manifest_path,
                        "integration manifest",
                    )?;
                if manifest.sampling_correlation_mode != self.sampling_correlation_mode() {
                    return Err(eyre!(
                        "Workspace integration sampling mode does not match the requested mode; use --restart to switch between correlated and uncorrelated integration"
                    ));
                }
                let expected_slots = selected_slots
                    .iter()
                    .map(|slot| slot.slot_meta.clone())
                    .collect_vec();
                if manifest.slots != expected_slots {
                    return Err(eyre!(
                        "Workspace integration slots do not match the currently selected integrands"
                    ));
                }
                if manifest.targets != *targets {
                    warn!("targets have changed with respect to workspace, reverting changes");
                    *targets = manifest.targets.clone();
                }
                if manifest.integrand_fingerprints.len() != selected_slots.len() {
                    return Err(eyre!(
                        "Workspace integrand fingerprint metadata is inconsistent with the selected slots"
                    ));
                }
                if manifest.effective_model_parameters.len() != selected_slots.len() {
                    return Err(eyre!(
                        "Workspace effective model parameter metadata is inconsistent with the selected slots"
                    ));
                }
                if current_effective_model_parameters.len() != selected_slots.len() {
                    return Err(eyre!(
                        "Current effective model parameter metadata is inconsistent with the selected slots"
                    ));
                }
                let mismatched_fingerprint_slots = selected_slots
                    .iter()
                    .zip(manifest.integrand_fingerprints.iter())
                    .zip(current_integrand_fingerprints.iter())
                    .filter_map(|((slot, saved), current)| {
                        (saved != current).then(|| slot.slot_meta.key())
                    })
                    .collect_vec();
                if !mismatched_fingerprint_slots.is_empty() {
                    return Err(eyre!(
                        "Workspace integrand fingerprints do not match the current generated integrands for {}. Resume requires the exact same generated integrands; use --restart or restore the previous generation.",
                        mismatched_fingerprint_slots.join(", ")
                    ));
                }
                let mismatched_slots = selected_slots
                    .iter()
                    .zip(manifest.effective_model_parameters.iter())
                    .zip(current_effective_model_parameters.iter())
                    .filter_map(|((slot, saved_card), current_card)| {
                        (saved_card != current_card).then(|| slot.slot_meta.key())
                    })
                    .collect_vec();
                if !mismatched_slots.is_empty() {
                    return Err(eyre!(
                        "Workspace effective model parameters do not match the current state for {}. Resume requires an exact match; use --restart or restore the shared/per-integrand model parameters.",
                        mismatched_slots.join(", ")
                    ));
                }

                info!(
                    "{}",
                    "Found integration state, result of previous integration:".yellow()
                );
                info!("");

                let integration_state: IntegrationState =
                    bincode::decode_from_slice::<IntegrationState, _>(
                        &state_bytes,
                        bincode::config::standard(),
                    )
                    .expect("Could not deserialize state")
                    .0;

                for (slot, target) in selected_slots.iter().zip(targets.iter()) {
                    let settings_label = format!("workspace settings for {}", slot.slot_meta.key());
                    let workspace_settings: RuntimeSettings = RuntimeSettings::from_file(
                        slot_settings_path(workspace_path, &slot.slot_meta),
                        &settings_label,
                    )?;
                    let gloop_integrand = state
                        .process_list
                        .get_integrand_mut(slot.process_id, &slot.slot_meta.integrand_name)?;
                    let current_settings = gloop_integrand.get_mut_settings();
                    if *current_settings != workspace_settings {
                        warn!(
                            "settings for {} have changed with respect to workspace, reverting non-model changes",
                            slot.slot_meta.key()
                        );
                        let preserved_model_overrides = current_settings.model.clone();
                        *current_settings = workspace_settings;
                        current_settings.model = preserved_model_overrides;
                    }

                    let label = format!("itg {}", slot.slot_meta.key());
                    let saved_slot = integration_state
                        .all_integrals
                        .get(
                            selected_slots
                                .iter()
                                .position(|candidate| candidate.slot_meta == slot.slot_meta)
                                .expect("selected slot must exist"),
                        )
                        .expect("saved slot must exist");
                    print_integral_result(
                        &saved_slot.re,
                        &label,
                        integration_state.iter,
                        "re",
                        target.as_ref().map(|value| value.re),
                    );
                    print_integral_result(
                        &saved_slot.im,
                        &label,
                        integration_state.iter,
                        "im",
                        target.as_ref().map(|value| value.im),
                    );
                }
                info!("");
                warn!(
                    "Any changes to the settings will be ignored, integrate with the {} option for changes to take effect",
                    "--restart".blue()
                );
                info!("{}", "Resuming integration".yellow());

                Ok(Some(integration_state))
            }
            Err(_) => {
                info!("No integration state found, starting new integration");
                Ok(None)
            }
        }
    }

    fn warm_and_clone_integrands(
        &self,
        state: &mut State,
        selected_slots: &[ResolvedIntegrandSlot],
        slot_models: &[Model],
    ) -> Result<Vec<gammalooprs::integrands::process::ProcessIntegrand>> {
        info!(
            "Gammaloop now integrates {}",
            selected_slots
                .iter()
                .map(|slot| slot.slot_meta.key().green().bold().to_string())
                .join(", ")
        );

        selected_slots
            .iter()
            .zip(slot_models.iter())
            .map(|slot| {
                let (slot, model) = slot;
                let gloop_integrand = state
                    .process_list
                    .get_integrand_mut(slot.process_id, &slot.slot_meta.integrand_name)?;
                gloop_integrand.warm_up(model)?;
                Ok(gloop_integrand.clone())
            })
            .collect()
    }

    fn restore_workspace_observables(
        &self,
        workspace_path: &std::path::Path,
        selected_slots: &[ResolvedIntegrandSlot],
        integration_state: Option<&IntegrationState>,
        slot_integrands: &mut [gammalooprs::integrands::process::ProcessIntegrand],
    ) -> Result<()> {
        let Some(integration_state) = integration_state else {
            return Ok(());
        };
        if integration_state.iter == 0 {
            return Ok(());
        }

        for ((slot, integrand), slot_meta) in selected_slots
            .iter()
            .zip(slot_integrands.iter_mut())
            .zip(integration_state.slot_metas.iter())
        {
            debug_assert_eq!(slot.slot_meta, *slot_meta);
            if integrand.observable_snapshot_bundle().is_none() {
                continue;
            }

            let snapshot_path =
                latest_observable_resume_state_path(workspace_path, &slot.slot_meta);
            let snapshot =
                ObservableSnapshotBundle::from_json_file(&snapshot_path).map_err(|err| {
                    eyre!(
                        "Could not restore observable checkpoint for {} from {}: {err}",
                        slot.slot_meta.key(),
                        snapshot_path.display()
                    )
                })?;
            integrand.restore_observable_snapshot_bundle(&snapshot)?;
        }

        Ok(())
    }

    fn resolve_slot_models(
        &self,
        state: &State,
        selected_slots: &[ResolvedIntegrandSlot],
    ) -> Result<Vec<Model>> {
        selected_slots
            .iter()
            .map(|slot| {
                state.resolve_model_for_integrand(slot.process_id, &slot.slot_meta.integrand_name)
            })
            .collect()
    }

    fn resolve_effective_model_parameters(
        &self,
        state: &State,
        selected_slots: &[ResolvedIntegrandSlot],
    ) -> Result<Vec<SerializableInputParamCard<F<f64>>>> {
        selected_slots
            .iter()
            .map(|slot| {
                state.resolve_effective_model_parameter_card_for_integrand(
                    slot.process_id,
                    &slot.slot_meta.integrand_name,
                )
            })
            .collect()
    }

    fn resolve_integrand_fingerprints(
        &self,
        state: &mut State,
        selected_slots: &[ResolvedIntegrandSlot],
    ) -> Result<Vec<String>> {
        selected_slots
            .iter()
            .map(|slot| {
                state
                    .process_list
                    .get_integrand_mut(slot.process_id, &slot.slot_meta.integrand_name)?
                    .resume_fingerprint()
            })
            .collect()
    }

    fn validate_slot_compatibility(
        &self,
        selected_slots: &[ResolvedIntegrandSlot],
        slot_integrands: &[gammalooprs::integrands::process::ProcessIntegrand],
    ) -> Result<()> {
        if self.uncorrelated || slot_integrands.len() <= 1 {
            return Ok(());
        }

        let reference = &slot_integrands[0];
        let reference_grid = reference.create_grid();
        let reference_sampling = &reference.get_settings().sampling;

        for (slot, integrand) in selected_slots.iter().zip(slot_integrands.iter()).skip(1) {
            if integrand.get_settings().sampling != *reference_sampling {
                return Err(eyre!(
                    "Integrand '{}' does not share the same sampling settings as the leading integrand '{}'",
                    slot.slot_meta.key(),
                    selected_slots[0].slot_meta.key()
                ));
            }
            if !grids_have_compatible_sample_shape(&reference_grid, &integrand.create_grid()) {
                return Err(eyre!(
                    "Integrand '{}' does not have a sample shape compatible with the leading integrand '{}'",
                    slot.slot_meta.key(),
                    selected_slots[0].slot_meta.key()
                ));
            }
        }

        Ok(())
    }

    fn write_workspace_manifest_and_settings(
        &self,
        slot_integrands: &[gammalooprs::integrands::process::ProcessIntegrand],
        selected_slots: &[ResolvedIntegrandSlot],
        targets: &[Option<Complex<F<f64>>>],
        effective_model_parameters: &[SerializableInputParamCard<F<f64>>],
        workspace_path: &std::path::Path,
    ) -> Result<()> {
        let integrand_fingerprints = slot_integrands
            .iter()
            .map(|integrand| integrand.resume_fingerprint())
            .collect::<Result<Vec<_>>>()?;
        let manifest = IntegrationWorkspaceManifest {
            slots: selected_slots
                .iter()
                .map(|slot| slot.slot_meta.clone())
                .collect(),
            targets: targets.to_vec(),
            effective_model_parameters: effective_model_parameters.to_vec(),
            integrand_fingerprints,
            training_slot: 0,
            integrator_settings_slot: 0,
            sampling_correlation_mode: self.sampling_correlation_mode(),
        };
        manifest.to_file(workspace_manifest_path(workspace_path), true)?;

        for (slot, integrand) in selected_slots.iter().zip(slot_integrands.iter()) {
            let slot_workspace = slot_workspace_path(workspace_path, &slot.slot_meta);
            fs::create_dir_all(&slot_workspace)?;
            integrand
                .get_settings()
                .to_file(slot_settings_path(workspace_path, &slot.slot_meta), true)?;
        }

        Ok(())
    }

    fn collect_workspace_observable_snapshots(
        &self,
        workspace_path: &std::path::Path,
        slot_metas: impl IntoIterator<Item = SlotMeta>,
    ) -> Result<BTreeMap<String, ObservableSnapshotBundle>> {
        let mut observables = BTreeMap::new();
        for slot_meta in slot_metas {
            let snapshot_path = latest_observable_resume_state_path(workspace_path, &slot_meta);
            if !snapshot_path.exists() {
                continue;
            }
            let snapshot =
                ObservableSnapshotBundle::from_json_file(&snapshot_path).map_err(|err| {
                    eyre!(
                        "Could not load observable snapshot for {} from {}: {err}",
                        slot_meta.key(),
                        snapshot_path.display()
                    )
                })?;
            observables.insert(slot_meta.key(), snapshot);
        }
        Ok(observables)
    }

    pub fn run(
        &self,
        state: &mut State,
        global_cli_settings: &CLISettings,
    ) -> Result<IntegrationOutput> {
        if self.show_summary_only && self.restart {
            return Err(eyre!(
                "The integrate options `--show-summary-only` and `--restart` cannot be used together"
            ));
        }
        if self.max_table_width == 0 {
            return Err(eyre!(
                "The integrate option `--max-table-width` must be greater than zero"
            ));
        }
        if self.batch_size == Some(0) {
            return Err(eyre!(
                "The integrate option `--batch-size` must be greater than zero"
            ));
        }
        if self.batch_timing < 0.0 {
            return Err(eyre!(
                "The integrate option `--batch-timing` must be greater than or equal to zero"
            ));
        }
        if self.min_time_between_status_updates < 0.0 {
            return Err(eyre!(
                "The integrate option `--min-time-between-status-updates` must be greater than or equal to zero"
            ));
        }

        let default_workspace_path = self.default_workspace_path(global_cli_settings);
        let workspace_path = if let Some(p) = self.workspace_path.clone() {
            p
        } else {
            default_workspace_path.clone()
        };

        if self.show_summary_only {
            let (manifest, integration_state) = read_existing_workspace_state(&workspace_path)?;
            let selected_slots = if self.process.is_empty() && self.integrand_name.is_empty() {
                self.resolve_manifest_slots(state, &manifest)?
            } else {
                self.resolve_selected_slots(state)?
            };
            let expected_slots = selected_slots
                .iter()
                .map(|slot| slot.slot_meta.clone())
                .collect_vec();
            if manifest.slots != expected_slots {
                return Err(eyre!(
                    "Workspace integration slots do not match the currently selected integrands"
                ));
            }

            let mut targets = if self.target.is_empty() {
                manifest.targets.clone()
            } else {
                self.resolve_targets(&selected_slots)?
            };
            if targets != manifest.targets {
                warn!("targets have changed with respect to workspace, reverting changes");
                targets = manifest.targets.clone();
            }

            if integration_state.num_points == 0 {
                return Err(eyre!(
                    "No completed integration iteration is available in {}",
                    workspace_path.display()
                ));
            }

            let slot0_meta = manifest.slots.first().ok_or_else(|| {
                eyre!("Integration workspace does not contain any integrand slots")
            })?;
            let slot0_settings: RuntimeSettings = RuntimeSettings::from_file(
                slot_settings_path(&workspace_path, slot0_meta),
                &format!("workspace settings for {}", slot0_meta.key()),
            )?;
            let render_options = self.build_render_options(&slot0_settings, true);
            emit_integration_status_via_tracing(
                IntegrationStatusKind::Final,
                render_saved_integration_summary(&integration_state, &targets, &render_options),
            )?;

            return Ok(IntegrationOutput {
                result: build_integration_result(&integration_state, &targets),
                observables: self.collect_workspace_observable_snapshots(
                    &workspace_path,
                    selected_slots.iter().map(|slot| slot.slot_meta.clone()),
                )?,
                workspace_path,
            });
        }

        let selected_slots = self.resolve_selected_slots(state)?;

        if self.restart && workspace_path.exists() {
            fs::remove_dir_all(&workspace_path)?;
        }

        if !workspace_path.exists() {
            fs::create_dir_all(&workspace_path)?;
            info!(
                "Created workspace directory at {}",
                workspace_path.display()
            );
        }
        let mut targets = self.resolve_targets(&selected_slots)?;
        let slot_models = self.resolve_slot_models(state, &selected_slots)?;
        let effective_model_parameters =
            self.resolve_effective_model_parameters(state, &selected_slots)?;
        let integrand_fingerprints = self.resolve_integrand_fingerprints(state, &selected_slots)?;
        let integration_state = self.load_or_prepare_workspace_state(
            state,
            &selected_slots,
            &effective_model_parameters,
            &integrand_fingerprints,
            &workspace_path,
            &mut targets,
        )?;
        let mut slot_integrands =
            self.warm_and_clone_integrands(state, &selected_slots, &slot_models)?;
        self.restore_workspace_observables(
            &workspace_path,
            &selected_slots,
            integration_state.as_ref(),
            &mut slot_integrands,
        )?;
        self.validate_slot_compatibility(&selected_slots, &slot_integrands)?;
        let slot_settings = slot_integrands
            .iter()
            .map(|integrand| integrand.get_settings().clone())
            .collect_vec();
        let render_options =
            self.build_render_options(&slot_settings[0], !self.no_show_integration_statistics);
        self.write_workspace_manifest_and_settings(
            &slot_integrands,
            &selected_slots,
            &targets,
            &effective_model_parameters,
            &workspace_path,
        )?;

        let n_cores = self
            .n_cores
            .unwrap_or(global_cli_settings.global.n_cores.integrate);

        let stderr_is_tty = io::stderr().is_terminal();
        let stream_updates = !self.no_stream_updates && stderr_is_tty;
        let stream_iterations = !self.no_stream_iterations && stderr_is_tty;
        if !stderr_is_tty && (!self.no_stream_updates || !self.no_stream_iterations) {
            info!(
                "Streaming integration updates disabled because stderr is not a TTY; live updates will be skipped and iteration summaries will be logged."
            );
        }
        let mut stream_renderer = (stream_updates || stream_iterations).then(StreamRenderer::new);

        let result = havana_integrate(
            slot_settings,
            self.sampling_correlation_mode(),
            slot_models,
            selected_slots
                .iter()
                .map(|slot| slot.slot_meta.clone())
                .collect(),
            slot_integrands
                .into_iter()
                .map(Integrand::ProcessIntegrand)
                .collect(),
            n_cores,
            targets,
            integration_state,
            Some(workspace_path.clone()),
            self.workspace_snapshot_control(),
            self.build_batching_settings(stream_updates),
            render_options,
            move |kind, status_block| {
                if let Some(renderer) = stream_renderer.as_mut() {
                    match kind {
                        IntegrationStatusKind::Live => {
                            if stream_updates {
                                renderer.render(&status_block)?;
                            }
                        }
                        IntegrationStatusKind::Iteration => {
                            if stream_iterations {
                                renderer.render(&status_block)?;
                            } else {
                                renderer.clear()?;
                                emit_integration_status_via_tracing(kind, &status_block)?;
                            }
                        }
                        IntegrationStatusKind::Final => {
                            renderer.clear()?;
                            emit_integration_status_via_tracing(kind, &status_block)?;
                        }
                    }
                } else {
                    if kind != IntegrationStatusKind::Live {
                        emit_integration_status_via_tracing(kind, &status_block)?;
                    }
                }

                Ok(())
            },
        )?;

        Ok(IntegrationOutput {
            observables: self.collect_workspace_observable_snapshots(
                &workspace_path,
                selected_slots.iter().map(|slot| slot.slot_meta.clone()),
            )?,
            result,
            workspace_path,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::Integrate;
    use super::IntegrationOutput;
    use super::ResolvedIntegrandSlot;
    use super::{ContributionSortOption, ShowPhaseOption};
    use crate::{state::ProcessRef, CLISettings, SessionSettings, StateSettings};
    use gammalooprs::{
        integrate::{ContributionSortMode, IntegrationStatusPhaseDisplay, SlotMeta},
        observables::ObservableSnapshotBundle,
        settings::{runtime::IntegratedPhase, IntegratorSettings, RuntimeSettings},
        utils::F,
    };
    use spenso::algebra::complex::Complex;
    use std::collections::BTreeMap;
    use std::path::PathBuf;
    use symbolica::numerical_integration::{ContinuousGrid, DiscreteGrid, Grid};

    fn resolved_slot(process_name: &str, integrand_name: &str) -> ResolvedIntegrandSlot {
        ResolvedIntegrandSlot {
            process_id: 0,
            slot_meta: SlotMeta {
                process_name: process_name.to_string(),
                integrand_name: integrand_name.to_string(),
            },
        }
    }

    #[test]
    fn read_only_state_uses_cwd_workspace_default() {
        let integrate = Integrate::default();
        let mut settings = CLISettings::default();
        settings.state = StateSettings {
            folder: PathBuf::from("/tmp/saved_state"),
            name: Some("gg_hhh_1l".to_string()),
        };
        settings.session = SessionSettings {
            read_only_state: true,
            ..SessionSettings::default()
        };

        assert_eq!(
            integrate.default_workspace_path(&settings),
            PathBuf::from("./gg_hhh_1l_integration_workspace")
        );

        settings.state.name = None;
        assert_eq!(
            integrate.default_workspace_path(&settings),
            PathBuf::from("./integration_workspace")
        );
    }

    #[test]
    fn integrate_defaults_to_showing_overall_max_weight_info() {
        assert!(Integrate::default().show_max_weight_info);
        assert!(!Integrate::default().show_max_weight_info_for_discrete_bins);
    }

    #[test]
    fn truncate_ansi_line_limits_visible_width() {
        let line = "\u{1b}[32mhello\u{1b}[0m world";
        let truncated = super::truncate_ansi_line(&line, 7);

        assert!(truncated.contains('\u{1b}'));
        assert!(truncated.ends_with("\u{1b}[0m"));
        assert!(!truncated.contains("world"));
    }

    #[test]
    fn prepare_stream_block_truncates_each_line_independently() {
        let block = "123456789\nabcdefghi";
        let prepared = super::prepare_stream_block(block, 5);

        assert_eq!(prepared, "12345\nabcde");
    }

    #[test]
    fn resolve_targets_supports_keyed_multi_slot_targets() {
        let integrate = Integrate {
            target: vec![
                "triangle@LO=1.0,2.0".to_string(),
                "box@scalar_box=3.0,4.0".to_string(),
            ],
            ..Integrate::default()
        };
        let slots = vec![
            resolved_slot("triangle", "LO"),
            resolved_slot("box", "scalar_box"),
        ];

        let targets = integrate.resolve_targets(&slots).unwrap();

        assert_eq!(
            targets,
            vec![
                Some(Complex::new(F(1.0), F(2.0))),
                Some(Complex::new(F(3.0), F(4.0))),
            ]
        );
    }

    #[test]
    fn resolve_targets_preserves_single_slot_legacy_target_format() {
        let integrate = Integrate {
            target: vec!["1.0".to_string(), "2.0".to_string()],
            ..Integrate::default()
        };
        let slots = vec![resolved_slot("triangle", "LO")];

        let targets = integrate.resolve_targets(&slots).unwrap();

        assert_eq!(targets, vec![Some(Complex::new(F(1.0), F(2.0)))]);
    }

    #[test]
    fn from_slots_preserves_slot_order() {
        let integrate = Integrate::from_slots([
            (ProcessRef::Id(2), "first"),
            (ProcessRef::Name("triangle".to_string()), "second"),
        ]);

        assert_eq!(
            integrate.process,
            vec![ProcessRef::Id(2), ProcessRef::Name("triangle".to_string())]
        );
        assert_eq!(
            integrate.integrand_name,
            vec!["first".to_string(), "second".to_string()]
        );
    }

    #[test]
    fn integration_output_single_slot_observables_returns_only_bundle() {
        let mut observables = BTreeMap::new();
        observables.insert(
            "box@scalar_box".to_string(),
            ObservableSnapshotBundle::default(),
        );
        let output = IntegrationOutput {
            observables,
            ..IntegrationOutput::default()
        };

        assert!(output.single_slot_observables().is_some());
    }

    #[test]
    fn sample_shape_compatibility_requires_matching_grid_topology() {
        let continuous_2d = Grid::Continuous(ContinuousGrid::new(2, 8, 0, None, false));
        let continuous_3d = Grid::Continuous(ContinuousGrid::new(3, 8, 0, None, false));
        let discrete_2 = Grid::Discrete(DiscreteGrid::new(vec![None, None], F(10.0), false));
        let discrete_nested = Grid::Discrete(DiscreteGrid::new(
            vec![
                Some(Grid::Continuous(ContinuousGrid::new(2, 8, 0, None, false))),
                Some(Grid::Continuous(ContinuousGrid::new(2, 8, 0, None, false))),
            ],
            F(10.0),
            false,
        ));

        assert!(super::grids_have_compatible_sample_shape(
            &continuous_2d,
            &Grid::Continuous(ContinuousGrid::new(2, 4, 0, None, false))
        ));
        assert!(!super::grids_have_compatible_sample_shape(
            &continuous_2d,
            &continuous_3d
        ));
        assert!(!super::grids_have_compatible_sample_shape(
            &continuous_2d,
            &discrete_2
        ));
        assert!(!super::grids_have_compatible_sample_shape(
            &discrete_2,
            &discrete_nested
        ));
    }

    #[test]
    fn selected_show_phase_resolves_from_slot_zero_phase() {
        assert_eq!(
            ShowPhaseOption::Selected.resolve(IntegratedPhase::Real),
            IntegrationStatusPhaseDisplay::Real
        );
        assert_eq!(
            ShowPhaseOption::Selected.resolve(IntegratedPhase::Imag),
            IntegrationStatusPhaseDisplay::Imag
        );
    }

    #[test]
    fn build_render_options_propagates_discrete_monitoring_flags() {
        let integrate = Integrate {
            show_phase: ShowPhaseOption::Imag,
            show_max_weight_info: true,
            show_top_discrete_grid: true,
            show_discrete_contributions_sum: true,
            sort_contributions: ContributionSortOption::Integral,
            show_max_weight_info_for_discrete_bins: true,
            ..Integrate::default()
        };
        let settings = RuntimeSettings {
            integrator: IntegratorSettings {
                integrated_phase: IntegratedPhase::Real,
                ..Default::default()
            },
            ..RuntimeSettings::default()
        };

        let render_options = integrate.build_render_options(&settings, false);

        assert_eq!(
            render_options.phase_display,
            IntegrationStatusPhaseDisplay::Imag
        );
        assert!(!render_options.show_statistics);
        assert!(render_options.show_max_weight_details);
        assert!(render_options.show_top_discrete_grid);
        assert!(render_options.show_discrete_contributions_sum);
        assert_eq!(
            render_options.contribution_sort,
            ContributionSortMode::Integral
        );
        assert!(render_options.show_max_weight_info_for_discrete_bins);
    }
}
