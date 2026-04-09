use std::time::Duration;

use itertools::Itertools;
use spenso::algebra::algebraic_traits::IsZero;
use spenso::algebra::complex::Complex;
use symbolica::numerical_integration::StatisticsAccumulator;

use crate::{
    integrands::evaluation::StatisticsCounter, settings::runtime::IntegrationStatisticsSnapshot,
    utils, utils::F,
};

use super::{
    ComplexAccumulator, DiscreteGridAccumulatorSummary, IntegrationState,
    display::{
        DisplayField, StyledText, TextStyle, UncertaintyNotation, format_abbreviated_count,
        format_iteration_points, format_max_eval_sample, format_signed_uncertainty,
        format_significant_percentage, format_total_points, styled_bin_description,
    },
    max_eval_entry, max_weight_impact, summary_at_path,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum IntegrationStatusKind {
    Live,
    Iteration,
    Final,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum IntegrationStatusPhaseDisplay {
    Both,
    Real,
    Imag,
}

impl IntegrationStatusPhaseDisplay {
    pub(crate) fn shows_real(self) -> bool {
        matches!(self, Self::Both | Self::Real)
    }

    pub(crate) fn shows_imag(self) -> bool {
        matches!(self, Self::Both | Self::Imag)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ContributionSortMode {
    Index,
    Integral,
    Error,
}

#[derive(Clone, Debug, PartialEq)]
pub struct IntegrationStatusViewOptions {
    pub phase_display: IntegrationStatusPhaseDisplay,
    pub training_phase_display: IntegrationStatusPhaseDisplay,
    pub training_slot: usize,
    pub slot_training_phase_displays: Vec<IntegrationStatusPhaseDisplay>,
    pub per_slot_training_phase: bool,
    pub target_relative_accuracy: Option<f64>,
    pub target_absolute_accuracy: Option<f64>,
    pub show_statistics: bool,
    pub show_max_weight_details: bool,
    pub show_top_discrete_grid: bool,
    pub show_discrete_contributions_sum: bool,
    pub contribution_sort: ContributionSortMode,
    pub show_max_weight_info_for_discrete_bins: bool,
}

impl IntegrationStatusViewOptions {
    pub(crate) fn for_final(self) -> Self {
        Self {
            show_statistics: true,
            ..self
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct LiveIterationProgress {
    pub(crate) completed_points: usize,
    pub(crate) target_points: usize,
}

pub(crate) struct StatusUpdateBuildRequest<'a> {
    pub(crate) kind: IntegrationStatusKind,
    pub(crate) integration_state: &'a IntegrationState,
    pub(crate) cores: usize,
    pub(crate) elapsed_time: Duration,
    pub(crate) iteration_elapsed_time: Duration,
    pub(crate) cur_points: usize,
    pub(crate) total_points_display: usize,
    pub(crate) n_samples_evaluated: usize,
    pub(crate) targets: &'a [Option<Complex<F<f64>>>],
    pub(crate) render_options: &'a IntegrationStatusViewOptions,
    pub(crate) live_progress: Option<LiveIterationProgress>,
}

impl<'a> StatusUpdateBuildRequest<'a> {
    pub(crate) fn new(
        kind: IntegrationStatusKind,
        integration_state: &'a IntegrationState,
        targets: &'a [Option<Complex<F<f64>>>],
        render_options: &'a IntegrationStatusViewOptions,
    ) -> Self {
        Self {
            kind,
            integration_state,
            cores: 0,
            elapsed_time: Duration::ZERO,
            iteration_elapsed_time: Duration::ZERO,
            cur_points: 0,
            total_points_display: 0,
            n_samples_evaluated: 0,
            targets,
            render_options,
            live_progress: None,
        }
    }

    pub(crate) fn with_timing(
        mut self,
        cores: usize,
        elapsed_time: Duration,
        iteration_elapsed_time: Duration,
        cur_points: usize,
        total_points_display: usize,
        n_samples_evaluated: usize,
    ) -> Self {
        self.cores = cores;
        self.elapsed_time = elapsed_time;
        self.iteration_elapsed_time = iteration_elapsed_time;
        self.cur_points = cur_points;
        self.total_points_display = total_points_display;
        self.n_samples_evaluated = n_samples_evaluated;
        self
    }

    pub(crate) fn with_live_progress(
        mut self,
        live_progress: Option<LiveIterationProgress>,
    ) -> Self {
        self.live_progress = live_progress;
        self
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct StatusMeta {
    pub(crate) elapsed_time: Duration,
    pub(crate) iteration_elapsed_time: Duration,
    pub(crate) iteration: usize,
    pub(crate) current_iteration_points: usize,
    pub(crate) total_points: usize,
    pub(crate) n_samples_evaluated: usize,
    pub(crate) cores: usize,
    pub(crate) training_slot: usize,
    pub(crate) training_phase_display: IntegrationStatusPhaseDisplay,
    pub(crate) live_progress: Option<LiveIterationProgress>,
    pub(crate) show_eta_to_target: bool,
    pub(crate) eta_to_target_specification: Option<String>,
    pub(crate) eta_to_target: Option<Duration>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct TargetAccuracyStatus {
    pub(crate) relative_reached: bool,
    pub(crate) absolute_reached: bool,
    pub(crate) eta_to_target: Option<Duration>,
}

impl TargetAccuracyStatus {
    pub(crate) fn is_reached(self) -> bool {
        self.relative_reached || self.absolute_reached
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub(crate) enum ComponentKind {
    Real,
    Imag,
}

impl ComponentKind {
    pub(crate) fn from_training_phase_display(
        display: IntegrationStatusPhaseDisplay,
    ) -> Option<Self> {
        match display {
            IntegrationStatusPhaseDisplay::Real => Some(Self::Real),
            IntegrationStatusPhaseDisplay::Imag => Some(Self::Imag),
            IntegrationStatusPhaseDisplay::Both => None,
        }
    }

    pub(crate) fn all_for_display(display: IntegrationStatusPhaseDisplay) -> Vec<Self> {
        let mut components = Vec::new();
        if display.shows_real() {
            components.push(Self::Real);
        }
        if display.shows_imag() {
            components.push(Self::Imag);
        }
        components
    }

    pub(crate) fn tag(self) -> &'static str {
        match self {
            Self::Real => "re",
            Self::Imag => "im",
        }
    }

    pub(crate) fn phase_name(self) -> &'static str {
        match self {
            Self::Real => "real",
            Self::Imag => "imag",
        }
    }

    pub(crate) fn text_style(self) -> TextStyle {
        match self {
            Self::Real => TextStyle::pink().bold(),
            Self::Imag => TextStyle::yellow().bold(),
        }
    }

    pub(crate) fn label_display(self) -> StyledText {
        StyledText::styled(self.tag(), self.text_style())
    }

    fn display_field(self) -> DisplayField<Self> {
        DisplayField::new(self, self.label_display())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ContributionKind {
    All,
    Sum,
    Bin(usize),
}

#[derive(Clone, Debug)]
pub struct StatusUpdate {
    pub(crate) kind: IntegrationStatusKind,
    pub(crate) meta: StatusMeta,
    pub(crate) targets: Vec<Option<Complex<F<f64>>>>,
    pub(crate) slot_training_phase_displays: Vec<IntegrationStatusPhaseDisplay>,
    pub(crate) per_slot_training_phase: bool,
    pub(crate) main_results: MainResultsSection,
    pub(crate) max_weight_details: Option<MaxWeightDetailsSection>,
    pub(crate) discrete_max_weight_details: Option<DiscreteMaxWeightDetailsSection>,
    pub(crate) statistics: Option<StatisticsSection>,
}

impl StatusUpdate {
    pub fn kind(&self) -> IntegrationStatusKind {
        self.kind
    }

    pub fn is_initial_live_status(&self) -> bool {
        self.kind == IntegrationStatusKind::Live
            && self
                .meta
                .live_progress
                .is_some_and(|progress| progress.completed_points == 0)
    }

    pub(crate) fn training_component(&self) -> Option<ComponentKind> {
        ComponentKind::from_training_phase_display(self.meta.training_phase_display)
    }

    pub(crate) fn statistics_snapshot(&self) -> Option<IntegrationStatisticsSnapshot> {
        self.statistics
            .as_ref()
            .map(StatisticsSection::global_snapshot)
    }

    pub(crate) fn training_target(&self) -> Option<F<f64>> {
        self.target_for_slot(self.meta.training_slot)
    }

    pub(crate) fn target_for_slot(&self, slot_index: usize) -> Option<F<f64>> {
        let component = self.training_component()?;
        self.target_for_slot_component(slot_index, component)
    }

    pub(crate) fn target_for_slot_component(
        &self,
        slot_index: usize,
        component: ComponentKind,
    ) -> Option<F<f64>> {
        let target = self.targets.get(slot_index)?.as_ref()?;
        Some(match component {
            ComponentKind::Real => target.re,
            ComponentKind::Imag => target.im,
        })
    }

    pub(crate) fn target_display_for_slot_component(
        &self,
        slot_index: usize,
        component: ComponentKind,
    ) -> Option<StyledText> {
        self.target_for_slot_component(slot_index, component)
            .map(|value| {
                StyledText::styled(
                    format_target_value_for_display(value.0),
                    TextStyle::blue().bold(),
                )
            })
    }

    pub(crate) fn training_component_for_slot(&self, slot_index: usize) -> Option<ComponentKind> {
        let phase = if self.per_slot_training_phase {
            self.slot_training_phase_displays
                .get(slot_index)
                .copied()
                .unwrap_or(self.meta.training_phase_display)
        } else if slot_index == self.meta.training_slot {
            self.meta.training_phase_display
        } else {
            return None;
        };
        match phase {
            IntegrationStatusPhaseDisplay::Real => Some(ComponentKind::Real),
            IntegrationStatusPhaseDisplay::Imag => Some(ComponentKind::Imag),
            IntegrationStatusPhaseDisplay::Both => None,
        }
    }

    pub(crate) fn slot_component_selected_for_training(
        &self,
        slot_index: usize,
        component: ComponentKind,
    ) -> bool {
        self.training_component_for_slot(slot_index) == Some(component)
    }

    pub(crate) fn target_deltas_for_row_slot(
        &self,
        row: &MainResultsRow,
        slot_index: usize,
    ) -> (Option<DisplayField<f64>>, Option<DisplayField<f64>>) {
        let Some(cell) = row.slot_cell(slot_index) else {
            return (None, None);
        };
        let Some(value) = cell.value.as_ref() else {
            return (None, None);
        };
        let target = self.target_for_slot_component(slot_index, row.component.raw);
        format_delta_fields_from_estimate(value.raw.0, value.raw.1, target)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum StatisticsScope {
    Global,
    Slot(usize),
}

impl StatusMeta {
    pub(crate) fn iteration_progress_ratio(&self) -> Option<f64> {
        self.live_progress.map(|progress| {
            if progress.target_points == 0 {
                0.0
            } else {
                progress.completed_points as f64 / progress.target_points as f64
            }
        })
    }

    pub(crate) fn iteration_eta(&self) -> Option<Duration> {
        let progress = self.live_progress?;
        if progress.completed_points == 0 || self.iteration_elapsed_time.is_zero() {
            return None;
        }

        let processed_per_second =
            progress.completed_points as f64 / self.iteration_elapsed_time.as_secs_f64();
        if processed_per_second <= 0.0 {
            return None;
        }

        let remaining_points = progress
            .target_points
            .saturating_sub(progress.completed_points);
        Some(utils::duration_from_secs_f64_saturating(
            remaining_points as f64 / processed_per_second,
        ))
    }

    pub(crate) fn total_sample_rate_per_second(&self) -> Option<f64> {
        if self.total_points == 0 || self.elapsed_time.is_zero() {
            return None;
        }
        Some(self.total_points as f64 / self.elapsed_time.as_secs_f64())
    }

    pub(crate) fn sample_core_time(&self) -> Option<String> {
        if self.n_samples_evaluated == 0 {
            return None;
        }

        Some(utils::format_evaluation_time_from_f64(
            self.elapsed_time.as_secs_f64() / (self.n_samples_evaluated as f64)
                * (self.cores as f64),
        ))
    }

    pub(crate) fn eta_to_target(&self) -> Option<Duration> {
        self.eta_to_target
    }

    pub(crate) fn eta_to_target_specification(&self) -> Option<&str> {
        self.eta_to_target_specification.as_deref()
    }
}

pub(crate) fn evaluate_target_accuracy(
    integration_state: &IntegrationState,
    total_points: usize,
    elapsed_time: Duration,
    targets: &[Option<Complex<F<f64>>>],
    training_phase_display: IntegrationStatusPhaseDisplay,
    target_relative_accuracy: Option<f64>,
    target_absolute_accuracy: Option<f64>,
) -> TargetAccuracyStatus {
    if target_relative_accuracy.is_none() && target_absolute_accuracy.is_none() {
        return TargetAccuracyStatus::default();
    }
    if total_points == 0 {
        return TargetAccuracyStatus::default();
    }

    let components = ComponentKind::all_for_display(training_phase_display);
    if components.is_empty() {
        return TargetAccuracyStatus::default();
    }
    if integration_state.all_integrals.is_empty() {
        return TargetAccuracyStatus::default();
    }

    let absolute_target_error = target_absolute_accuracy.filter(|target| *target >= 0.0);
    let relative_target_accuracy = target_relative_accuracy.filter(|target| *target >= 0.0);

    let mut absolute_reached = absolute_target_error.is_some();
    let mut relative_reached = relative_target_accuracy.is_some();
    let mut saw_absolute_constraint = false;
    let mut saw_relative_constraint = false;
    let mut absolute_eta_to_target = None;
    let mut relative_eta_to_target = None;

    for (slot_index, accumulator) in integration_state.all_integrals.iter().enumerate() {
        let target = targets.get(slot_index).and_then(|target| target.as_ref());
        for component in &components {
            let (avg, err) = match component {
                ComponentKind::Real => (accumulator.re.avg.0, accumulator.re.err.0),
                ComponentKind::Imag => (accumulator.im.avg.0, accumulator.im.err.0),
            };

            if let Some(target_error) = absolute_target_error {
                saw_absolute_constraint = true;
                absolute_reached &= err <= target_error;
                absolute_eta_to_target = max_duration_option(
                    absolute_eta_to_target,
                    estimate_eta_to_target(total_points, elapsed_time, err, Some(target_error)),
                );
            }

            if let Some(target_accuracy) = relative_target_accuracy {
                let reference = target
                    .map(|target| match component {
                        ComponentKind::Real => target.re.0,
                        ComponentKind::Imag => target.im.0,
                    })
                    .unwrap_or(avg)
                    .abs();
                if reference != 0.0 {
                    saw_relative_constraint = true;
                    let target_error = target_accuracy * reference;
                    relative_reached &= err <= target_error;
                    relative_eta_to_target = max_duration_option(
                        relative_eta_to_target,
                        estimate_eta_to_target(total_points, elapsed_time, err, Some(target_error)),
                    );
                }
            }
        }
    }

    if !saw_absolute_constraint {
        absolute_reached = false;
        absolute_eta_to_target = None;
    }
    if !saw_relative_constraint {
        relative_reached = false;
        relative_eta_to_target = None;
    }

    let eta_to_target = min_duration_option(absolute_eta_to_target, relative_eta_to_target);

    TargetAccuracyStatus {
        relative_reached,
        absolute_reached,
        eta_to_target,
    }
}

impl MainResultsSection {
    pub(crate) fn all_rows(&self) -> impl Iterator<Item = &MainResultsRow> {
        self.row_groups.iter().flat_map(|group| group.rows.iter())
    }

    pub(crate) fn row_groups_of_kind(
        &self,
        kind: MainResultsRowGroupKind,
    ) -> impl Iterator<Item = &MainResultsRowGroup> {
        self.row_groups
            .iter()
            .filter(move |group| group.kind == kind)
    }

    pub(crate) fn find_row(
        &self,
        contribution: ContributionKind,
        component: ComponentKind,
    ) -> Option<&MainResultsRow> {
        self.all_rows()
            .find(|row| row.contribution.raw == contribution && row.component.raw == component)
    }

    pub(crate) fn metadata_headers(&self) -> Vec<StyledText> {
        let mut headers = vec![
            styled_colored("χ²/dof", TextStyle::blue().bold()),
            styled_colored("mwi", TextStyle::blue().bold()),
        ];
        if self.has_target_columns {
            headers.push(styled_colored("Δ [σ]", TextStyle::blue().bold()));
            headers.push(styled_colored("Δ [%]", TextStyle::blue().bold()));
        }
        headers
    }

    pub(crate) fn component_header(&self) -> StyledText {
        styled_plain("")
    }

    pub(crate) fn summary_headers(&self) -> Vec<StyledText> {
        let mut headers = vec![self.contribution_header.clone(), self.component_header()];
        headers.extend(self.slot_headers.iter().cloned());
        headers
    }

    pub(crate) fn discrete_headers(&self) -> Vec<StyledText> {
        vec![
            self.contribution_header.clone(),
            self.component_header(),
            styled_plain("integral"),
            styled_plain("% err"),
            styled_plain("χ^2"),
            styled_plain("m.w.i"),
            styled_plain("sample %"),
            styled_plain("# samples"),
            styled_plain("pdf"),
        ]
    }

    pub(crate) fn selected_bin_detail_headers(&self) -> Vec<StyledText> {
        vec![
            styled_plain("Integrand"),
            styled_plain("integral"),
            styled_plain("% err"),
            styled_plain("χ^2"),
            styled_plain("m.w.i"),
            styled_plain("sample %"),
            styled_plain("# samples"),
            styled_plain("pdf"),
        ]
    }
}

impl MainResultsRow {
    pub(crate) fn slot_cell(&self, slot_index: usize) -> Option<&MainTableSlotCells> {
        self.slot_cells.get(slot_index)
    }
}

impl MaxWeightDetailsSection {
    pub(crate) fn title(&self) -> StyledText {
        styled_colored("Maximum weight details", TextStyle::green().bold())
    }

    pub(crate) fn headers(&self) -> Vec<StyledText> {
        vec![
            styled_colored("Integrand", TextStyle::blue().bold()),
            styled_plain(""),
            styled_colored("Max eval", TextStyle::blue().bold()),
            styled_colored("Max eval coordinates", TextStyle::blue().bold()),
        ]
    }
}

impl DiscreteMaxWeightDetailsSection {
    pub(crate) fn title(&self) -> StyledText {
        styled_colored(
            "Maximum weight details by discrete bin",
            TextStyle::green().bold(),
        )
    }

    pub(crate) fn summary_headers(&self) -> Vec<StyledText> {
        let mut headers = vec![self.contribution_header.clone(), styled_plain("")];
        headers.extend(self.slot_headers.iter().cloned());
        headers
    }

    pub(crate) fn coordinates_header(&self) -> StyledText {
        styled_colored("Max eval coordinates", TextStyle::blue().bold())
    }

    pub(crate) fn coordinate_headers(&self) -> Vec<StyledText> {
        vec![
            self.contribution_header.clone(),
            styled_plain(""),
            styled_plain("Integrand"),
            self.coordinates_header(),
        ]
    }
}

impl StatisticsSection {
    pub(crate) fn global_snapshot(&self) -> IntegrationStatisticsSnapshot {
        self.global.snapshot()
    }

    fn scoped_counter(&self, scope: StatisticsScope) -> Option<&StatisticsCounter> {
        match scope {
            StatisticsScope::Global => Some(&self.global),
            StatisticsScope::Slot(slot_index) => self.slot_counters.get(slot_index),
        }
    }

    fn scope_label(&self, scope: StatisticsScope) -> StyledText {
        let label = match scope {
            StatisticsScope::Global => "[global]".to_string(),
            StatisticsScope::Slot(slot_index) => self
                .slot_labels
                .get(slot_index)
                .map(|label| format!("[{label}]"))
                .unwrap_or_else(|| "[global]".to_string()),
        };
        styled_colored(format!(" {label}"), TextStyle::blue().bold())
    }

    fn scoped_title(&self, prefix: &str, scope: StatisticsScope) -> StyledText {
        let mut title = styled_colored(prefix, TextStyle::green().bold());
        title.append(self.scope_label(scope));
        title
    }

    fn uses_global_integrator_scope(scope: StatisticsScope) -> bool {
        matches!(scope, StatisticsScope::Global)
    }

    pub(crate) fn statistics_title(&self, scope: StatisticsScope) -> StyledText {
        self.scoped_title("Integration statistics", scope)
    }

    pub(crate) fn timing_title(&self, scope: StatisticsScope) -> StyledText {
        self.scoped_title("Timing composition", scope)
    }

    pub(crate) fn precision_title(&self, scope: StatisticsScope) -> StyledText {
        self.scoped_title("Precision mix", scope)
    }

    pub(crate) fn table_rows(&self, scope: StatisticsScope) -> Vec<StatisticsTableRow> {
        let snapshot = self
            .scoped_counter(scope)
            .map(StatisticsCounter::snapshot)
            .unwrap_or_else(|| self.global.snapshot());
        let integrator_value = if Self::uses_global_integrator_scope(scope) {
            styled_colored(
                utils::format_evaluation_time_from_f64(snapshot.average_integrator_time_seconds),
                TextStyle::green(),
            )
        } else {
            styled_plain("N/A")
        };
        vec![
            StatisticsTableRow {
                row_label: styled_colored("  timing", TextStyle::blue().bold()),
                entries: vec![
                    StatisticsTableEntry {
                        label: styled_plain("total"),
                        value: styled_colored(
                            utils::format_evaluation_time_from_f64(
                                snapshot.average_total_time_seconds,
                            ),
                            TextStyle::green(),
                        ),
                    },
                    StatisticsTableEntry {
                        label: styled_plain("param"),
                        value: styled_colored(
                            utils::format_evaluation_time_from_f64(
                                snapshot.average_parameterization_time_seconds,
                            ),
                            TextStyle::green(),
                        ),
                    },
                    StatisticsTableEntry {
                        label: styled_plain("itg"),
                        value: styled_colored(
                            utils::format_evaluation_time_from_f64(
                                snapshot.average_integrand_time_seconds,
                            ),
                            TextStyle::green(),
                        ),
                    },
                    StatisticsTableEntry {
                        label: styled_plain("evaluators"),
                        value: styled_colored(
                            utils::format_evaluation_time_from_f64(
                                snapshot.average_evaluator_time_seconds,
                            ),
                            TextStyle::green(),
                        ),
                    },
                ],
            },
            StatisticsTableRow {
                row_label: styled_colored("  evals", TextStyle::blue().bold()),
                entries: vec![
                    StatisticsTableEntry {
                        label: styled_plain("f64"),
                        value: styled_colored(
                            format!("{:.2}%", snapshot.f64_percentage),
                            TextStyle::green(),
                        ),
                    },
                    StatisticsTableEntry {
                        label: styled_plain("f128"),
                        value: styled_colored(
                            format!("{:.2}%", snapshot.f128_percentage),
                            TextStyle::blue(),
                        ),
                    },
                    StatisticsTableEntry {
                        label: styled_plain("arb"),
                        value: styled_plain(format!("{:.2}%", snapshot.arb_percentage)),
                    },
                    StatisticsTableEntry {
                        label: styled_plain("nans+unstable"),
                        value: styled_colored(
                            format!("{:.2}%", snapshot.nan_or_unstable_percentage),
                            if snapshot.nan_or_unstable_percentage > 0.0 {
                                TextStyle::red()
                            } else {
                                TextStyle::green()
                            },
                        ),
                    },
                ],
            },
            StatisticsTableRow {
                row_label: styled_colored("  events", TextStyle::blue().bold()),
                entries: vec![
                    StatisticsTableEntry {
                        label: styled_plain("evts #"),
                        value: styled_colored(
                            format_abbreviated_count(snapshot.generated_event_count),
                            TextStyle::green(),
                        ),
                    },
                    StatisticsTableEntry {
                        label: styled_plain("sel. %"),
                        value: snapshot
                            .selection_efficiency_percentage
                            .map(|value| styled_colored(format!("{value:.2}%"), TextStyle::green()))
                            .unwrap_or_else(|| styled_plain("N/A")),
                    },
                    StatisticsTableEntry {
                        label: styled_plain("obs"),
                        value: styled_colored(
                            utils::format_evaluation_time_from_f64(
                                snapshot.average_observable_time_seconds,
                            ),
                            TextStyle::green(),
                        ),
                    },
                    StatisticsTableEntry {
                        label: styled_plain("integrator"),
                        value: integrator_value,
                    },
                ],
            },
        ]
    }

    pub(crate) fn timing_mix_segments(&self, scope: StatisticsScope) -> Vec<StatisticsMixSegment> {
        let snapshot = self
            .scoped_counter(scope)
            .map(StatisticsCounter::snapshot)
            .unwrap_or_else(|| self.global.snapshot());
        let parameterization = snapshot.average_parameterization_time_seconds.max(0.0);
        let evaluator = snapshot.average_evaluator_time_seconds.max(0.0);
        let observable = snapshot.average_observable_time_seconds.max(0.0);
        let integrand_core =
            (snapshot.average_integrand_time_seconds.max(0.0) - observable - evaluator).max(0.0);
        let mut raw_segments = vec![
            (
                styled_colored("evaluators", TextStyle::green().bold()),
                evaluator,
            ),
            (
                styled_colored("itg_core", TextStyle::blue().bold()),
                integrand_core,
            ),
            (
                styled_colored("obs", TextStyle::yellow().bold()),
                observable,
            ),
            (
                styled_colored("param", TextStyle::red().bold()),
                parameterization,
            ),
        ];
        if Self::uses_global_integrator_scope(scope) {
            raw_segments.push((
                styled_plain("integrator"),
                snapshot.average_integrator_time_seconds.max(0.0),
            ));
        }
        let mut segments = normalize_mix_segments(raw_segments);
        segments.sort_by(|lhs, rhs| {
            rhs.percentage
                .partial_cmp(&lhs.percentage)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        segments
    }

    pub(crate) fn precision_mix_segments(
        &self,
        scope: StatisticsScope,
    ) -> Vec<StatisticsMixSegment> {
        let snapshot = self
            .scoped_counter(scope)
            .map(StatisticsCounter::snapshot)
            .unwrap_or_else(|| self.global.snapshot());
        normalize_mix_segments(vec![
            (
                styled_colored("f64", TextStyle::green()),
                snapshot.f64_percentage,
            ),
            (
                styled_colored("f128", TextStyle::blue()),
                snapshot.f128_percentage,
            ),
            (styled_plain("arb"), snapshot.arb_percentage),
            (
                styled_colored("unstbl.+nan.", TextStyle::red()),
                snapshot.nan_or_unstable_percentage,
            ),
        ])
    }

    pub(crate) fn stability_mix_segments(
        &self,
        scope: StatisticsScope,
    ) -> Vec<StatisticsMixSegment> {
        let snapshot = self
            .scoped_counter(scope)
            .map(StatisticsCounter::snapshot)
            .unwrap_or_else(|| self.global.snapshot());
        let unstable = (snapshot.nan_or_unstable_percentage - snapshot.nan_percentage).max(0.0);
        let stable = (100.0 - snapshot.nan_or_unstable_percentage).max(0.0);
        normalize_mix_segments(vec![
            (styled_colored("stable", TextStyle::green()), stable),
            (styled_colored("unstable", TextStyle::red()), unstable),
            (
                styled_colored("nan", TextStyle::red()),
                snapshot.nan_percentage,
            ),
        ])
    }
}

#[derive(Clone, Debug)]
pub(crate) struct MainResultsSection {
    pub(crate) header_left: StyledText,
    pub(crate) header_middle: StyledText,
    pub(crate) header_tail: StyledText,
    pub(crate) contribution_header: StyledText,
    pub(crate) slot_headers: Vec<StyledText>,
    pub(crate) has_discrete_columns: bool,
    pub(crate) has_target_columns: bool,
    pub(crate) row_groups: Vec<MainResultsRowGroup>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MainResultsRowGroupKind {
    All,
    Sum,
    Bins,
}

#[derive(Clone, Debug)]
pub(crate) struct MainResultsRowGroup {
    pub(crate) kind: MainResultsRowGroupKind,
    pub(crate) rows: Vec<MainResultsRow>,
}

#[derive(Clone, Debug)]
pub(crate) struct MainResultsRow {
    pub(crate) contribution: DisplayField<ContributionKind>,
    pub(crate) component: DisplayField<ComponentKind>,
    pub(crate) slot_cells: Vec<MainTableSlotCells>,
    pub(crate) chi_sq: Option<DisplayField<f64>>,
    pub(crate) delta_sigma: Option<DisplayField<f64>>,
    pub(crate) delta_percent: Option<DisplayField<f64>>,
    pub(crate) max_weight_impact: Option<DisplayField<f64>>,
}

#[derive(Clone, Debug, Default)]
pub(crate) struct MainTableSlotCells {
    pub(crate) value: Option<DisplayField<(F<f64>, F<f64>)>>,
    pub(crate) relative_error: Option<DisplayField<f64>>,
    pub(crate) chi_sq: Option<DisplayField<f64>>,
    pub(crate) max_weight_impact: Option<DisplayField<f64>>,
    pub(crate) sample_fraction: Option<DisplayField<f64>>,
    pub(crate) sample_count: Option<DisplayField<usize>>,
    pub(crate) target_pdf: Option<DisplayField<f64>>,
}

#[derive(Clone, Debug)]
pub(crate) struct MaxWeightDetailsSection {
    pub(crate) rows_by_slot: Vec<Vec<MaxWeightDetailsRow>>,
}

#[derive(Clone, Debug)]
pub(crate) struct MaxWeightDetailsRow {
    pub(crate) slot: DisplayField<String>,
    pub(crate) component_sign: DisplayField<(ComponentKind, bool)>,
    pub(crate) max_eval: DisplayField<F<f64>>,
    pub(crate) coordinates: DisplayField<String>,
}

#[derive(Clone, Debug)]
pub(crate) struct DiscreteMaxWeightDetailsSection {
    pub(crate) contribution_header: StyledText,
    pub(crate) slot_headers: Vec<StyledText>,
    pub(crate) row_groups: Vec<Vec<DiscreteMaxWeightRow>>,
}

#[derive(Clone, Debug)]
pub(crate) struct DiscreteMaxWeightRow {
    pub(crate) contribution: DisplayField<ContributionKind>,
    pub(crate) component_sign: DisplayField<(ComponentKind, bool)>,
    pub(crate) slot_values: Vec<Option<DisplayField<F<f64>>>>,
    pub(crate) slot_coordinates: Vec<SlotCoordinateEntry>,
}

#[derive(Clone, Debug)]
pub(crate) struct SlotCoordinateEntry {
    pub(crate) slot: DisplayField<String>,
    pub(crate) coordinates: DisplayField<String>,
}

#[derive(Clone, Debug)]
pub(crate) struct StatisticsSection {
    pub(crate) global: StatisticsCounter,
    pub(crate) slot_counters: Vec<StatisticsCounter>,
    pub(crate) slot_labels: Vec<String>,
}

#[derive(Clone, Debug)]
pub(crate) struct StatisticsTableRow {
    pub(crate) row_label: StyledText,
    pub(crate) entries: Vec<StatisticsTableEntry>,
}

#[derive(Clone, Debug)]
pub(crate) struct StatisticsTableEntry {
    pub(crate) label: StyledText,
    pub(crate) value: StyledText,
}

#[derive(Clone, Debug)]
pub(crate) struct StatisticsMixSegment {
    pub(crate) label: StyledText,
    pub(crate) percentage: f64,
}

fn component_accumulator(
    accumulator: &ComplexAccumulator,
    component: ComponentKind,
) -> &StatisticsAccumulator<F<f64>> {
    match component {
        ComponentKind::Real => &accumulator.re,
        ComponentKind::Imag => &accumulator.im,
    }
}

fn slot_component_summary(
    integration_state: &IntegrationState,
    slot_index: usize,
    component: ComponentKind,
) -> Option<&DiscreteGridAccumulatorSummary> {
    match component {
        ComponentKind::Real => integration_state.slot_re_summaries[slot_index].as_ref(),
        ComponentKind::Imag => integration_state.slot_im_summaries[slot_index].as_ref(),
    }
}

fn sum_estimate_error(summary: &DiscreteGridAccumulatorSummary) -> (F<f64>, F<f64>) {
    let (avg, err_sq) = summary
        .bins
        .iter()
        .fold((F(0.0), F(0.0)), |(avg, err_sq), bin| {
            (
                avg + bin.accumulator.avg,
                err_sq + bin.accumulator.err * bin.accumulator.err,
            )
        });
    (avg, F(err_sq.0.sqrt()))
}

fn total_processed_samples(summary: &DiscreteGridAccumulatorSummary) -> usize {
    summary
        .bins
        .iter()
        .map(|bin| bin.accumulator.processed_samples)
        .sum()
}

fn normalize_mix_segments(segments: Vec<(StyledText, f64)>) -> Vec<StatisticsMixSegment> {
    let total: f64 = segments.iter().map(|(_, value)| value.max(0.0)).sum();
    if total <= f64::EPSILON {
        return segments
            .into_iter()
            .map(|(label, _)| StatisticsMixSegment {
                label,
                percentage: 0.0,
            })
            .collect();
    }

    segments
        .into_iter()
        .map(|(label, value)| StatisticsMixSegment {
            label,
            percentage: value.max(0.0) / total * 100.0,
        })
        .collect()
}

fn estimate_eta_to_target(
    total_points: usize,
    elapsed_time: Duration,
    current_error: f64,
    target_error: Option<f64>,
) -> Option<Duration> {
    let target_error = target_error?;
    if target_error < 0.0 {
        return None;
    }
    if current_error <= target_error {
        return Some(Duration::ZERO);
    }
    if total_points == 0 || elapsed_time.is_zero() || current_error <= 0.0 {
        return None;
    }
    if target_error == 0.0 {
        return Some(Duration::MAX);
    }

    let rate_per_second = total_points as f64 / elapsed_time.as_secs_f64();
    if rate_per_second <= 0.0 {
        return None;
    }

    let target_points = total_points as f64 * (current_error / target_error).powi(2);
    if !target_points.is_finite() {
        return Some(Duration::MAX);
    }
    let remaining_points = (target_points - total_points as f64).max(0.0);
    Some(utils::duration_from_secs_f64_saturating(
        remaining_points / rate_per_second,
    ))
}

fn min_duration_option(lhs: Option<Duration>, rhs: Option<Duration>) -> Option<Duration> {
    match (lhs, rhs) {
        (Some(lhs), Some(rhs)) => Some(lhs.min(rhs)),
        (Some(lhs), None) => Some(lhs),
        (None, Some(rhs)) => Some(rhs),
        (None, None) => None,
    }
}

fn max_duration_option(lhs: Option<Duration>, rhs: Option<Duration>) -> Option<Duration> {
    match (lhs, rhs) {
        (Some(lhs), Some(rhs)) => Some(lhs.max(rhs)),
        (Some(lhs), None) => Some(lhs),
        (None, Some(rhs)) => Some(rhs),
        (None, None) => None,
    }
}

fn styled_plain(text: impl Into<String>) -> StyledText {
    StyledText::plain(text)
}

fn styled_colored(text: impl Into<String>, style: TextStyle) -> StyledText {
    StyledText::styled(text, style)
}

fn format_eta_to_target_specification(
    target_relative_accuracy: Option<f64>,
    target_absolute_accuracy: Option<f64>,
) -> Option<String> {
    let relative_spec = target_relative_accuracy
        .filter(|target| *target >= 0.0)
        .map(|target| {
            format!(
                "% err <= {}",
                format_percentage_target_value_for_display(target * 100.0)
            )
        });
    let absolute_spec = target_absolute_accuracy
        .filter(|target| *target >= 0.0)
        .map(|target| {
            format!(
                "err <= {}",
                format_positive_target_value_for_display(target)
            )
        });

    match (relative_spec, absolute_spec) {
        (Some(relative), Some(absolute)) => Some(format!("{relative} or {absolute}")),
        (Some(relative), None) => Some(relative),
        (None, Some(absolute)) => Some(absolute),
        (None, None) => None,
    }
}

fn format_target_value_for_display(value: f64) -> String {
    if value == 0.0 {
        return String::from("+0e0");
    }

    for precision in 0..=16 {
        let candidate = normalize_scientific_exponent(&format!("{:+.*e}", precision, value));
        let Ok(parsed) = candidate.parse::<f64>() else {
            continue;
        };
        if ulp_distance(parsed, value) <= 8 {
            return candidate;
        }
    }

    normalize_scientific_exponent(&format!("{:+.16e}", value))
}

fn normalize_scientific_exponent(formatted: &str) -> String {
    let Some((mantissa, exponent)) = formatted.rsplit_once('e') else {
        return formatted.to_string();
    };
    let exponent = exponent.parse::<i32>().unwrap_or_default();
    format!("{mantissa}e{exponent:+}")
}

fn ulp_distance(lhs: f64, rhs: f64) -> u64 {
    if lhs.is_nan() || rhs.is_nan() {
        return u64::MAX;
    }
    ordered_f64_bits(lhs).abs_diff(ordered_f64_bits(rhs))
}

fn ordered_f64_bits(value: f64) -> u64 {
    let bits = value.to_bits();
    if (bits >> 63) != 0 {
        !bits
    } else {
        bits | (1_u64 << 63)
    }
}

fn format_positive_target_value_for_display(value: f64) -> String {
    format_target_value_for_display(value)
        .trim_start_matches('+')
        .to_string()
}

fn format_percentage_target_value_for_display(value: f64) -> String {
    let formatted = if value.abs() <= 1.0e-4 {
        format_positive_target_value_for_display(value)
    } else {
        format_fixed_target_value_for_display(value)
    };
    format!("{formatted}%")
}

fn format_fixed_target_value_for_display(value: f64) -> String {
    if value == 0.0 {
        return String::from("0");
    }

    for precision in 0..=16 {
        let candidate = normalize_fixed_decimal(&format!("{value:.precision$}"));
        let Ok(parsed) = candidate.parse::<f64>() else {
            continue;
        };
        if ulp_distance(parsed, value) <= 8 {
            return candidate;
        }
    }

    format_positive_target_value_for_display(value)
}

fn normalize_fixed_decimal(formatted: &str) -> String {
    let trimmed = formatted.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() || trimmed == "-0" {
        String::from("0")
    } else {
        trimmed.to_string()
    }
}

fn format_value_field(avg: F<f64>, err: F<f64>) -> DisplayField<(F<f64>, F<f64>)> {
    let formatted = format_signed_uncertainty(avg, err, UncertaintyNotation::Scientific);
    let uncertainty_style = format_relative_error_field_from_estimate(avg, err)
        .map(|field| {
            field
                .display
                .spans
                .first()
                .map(|span| span.style)
                .unwrap_or(TextStyle::PLAIN)
                .bold()
        })
        .unwrap_or(TextStyle::PLAIN.bold());

    let display = if let Some(start) = formatted.find('(') {
        if let Some(end_offset) = formatted[start..].find(')') {
            let end = start + end_offset + 1;
            let mut styled = StyledText::new();
            styled.push_text(&formatted[..start], TextStyle::blue().bold());
            styled.push_text(&formatted[start..end], uncertainty_style);
            styled.push_text(&formatted[end..], TextStyle::blue().bold());
            styled
        } else {
            styled_colored(formatted, TextStyle::blue().bold())
        }
    } else {
        styled_colored(formatted, TextStyle::blue().bold())
    };

    DisplayField::new((avg, err), display)
}

fn format_relative_error_field_from_estimate(
    avg: F<f64>,
    err: F<f64>,
) -> Option<DisplayField<f64>> {
    if avg.is_zero() {
        return None;
    }

    let raw = (err / avg).abs().0 * 100.0;
    let formatted = format_significant_percentage(raw, 3, Some((1.0e4, 1.0e-4)));
    let style = if raw > 1.0 {
        TextStyle::red()
    } else {
        TextStyle::green()
    };
    Some(DisplayField::new(raw, styled_colored(formatted, style)))
}

fn format_relative_error_field(itg: &StatisticsAccumulator<F<f64>>) -> Option<DisplayField<f64>> {
    format_relative_error_field_from_estimate(itg.avg, itg.err)
}

fn format_chi_sq_field(
    itg: &StatisticsAccumulator<F<f64>>,
    i_iter: usize,
) -> Option<DisplayField<f64>> {
    if i_iter == 0 {
        return None;
    }
    let raw = itg.chi_sq.0 / (i_iter as f64);
    let style = if raw > 5.0 {
        TextStyle::red()
    } else {
        TextStyle::PLAIN
    };
    Some(DisplayField::new(
        raw,
        styled_colored(format!("{raw:.3}"), style),
    ))
}

fn format_delta_fields_from_estimate(
    avg: F<f64>,
    err: F<f64>,
    target: Option<F<f64>>,
) -> (Option<DisplayField<f64>>, Option<DisplayField<f64>>) {
    let Some(target) = target else {
        return (None, None);
    };

    let delta_sigma = if err.is_zero() {
        0.0
    } else {
        (target - avg).abs().0 / err.0
    };
    let delta_percent = if target.abs().is_non_zero() {
        (target - avg).abs().0 / target.abs().0 * 100.0
    } else {
        0.0
    };
    let is_outside_target = delta_sigma > 5.0
        || (target.abs().is_non_zero() && ((target - avg).abs() / target.abs()).0 > 0.01);
    let style = if is_outside_target {
        TextStyle::red()
    } else {
        TextStyle::green()
    };

    (
        Some(DisplayField::new(
            delta_sigma,
            styled_colored(format!("Δ = {:.3}σ", delta_sigma), style),
        )),
        Some(DisplayField::new(
            delta_percent,
            styled_colored(format!("Δ = {:.3}%", delta_percent), style),
        )),
    )
}

fn format_delta_fields(
    itg: &StatisticsAccumulator<F<f64>>,
    target: Option<F<f64>>,
) -> (Option<DisplayField<f64>>, Option<DisplayField<f64>>) {
    format_delta_fields_from_estimate(itg.avg, itg.err, target)
}

fn format_mwi_field(itg: &StatisticsAccumulator<F<f64>>) -> Option<DisplayField<f64>> {
    let raw = max_weight_impact(itg).0;
    let style = if raw > 1.0 {
        TextStyle::red()
    } else {
        TextStyle::PLAIN
    };
    Some(DisplayField::new(
        raw,
        styled_colored(format!("{raw:.4e}"), style),
    ))
}

fn contribution_display_field(
    contribution: ContributionKind,
    integration_state: &IntegrationState,
) -> DisplayField<ContributionKind> {
    let display = match contribution {
        ContributionKind::All => styled_colored("All", TextStyle::green().bold()),
        ContributionKind::Sum => styled_colored("Sum", TextStyle::green().bold()),
        ContributionKind::Bin(bin_index) => {
            if let (Some(axis_label), Some(descriptions)) = (
                integration_state
                    .first_non_trivial_discrete_label
                    .as_deref(),
                integration_state
                    .first_non_trivial_discrete_bin_descriptions
                    .as_ref(),
            ) && let Some(description) = descriptions.get(bin_index)
            {
                let mut text =
                    StyledText::plain(format_discrete_bin_prefix(bin_index, descriptions.len()));
                text.append(styled_bin_description(axis_label, description));
                text
            } else {
                styled_plain(format!("#{bin_index}"))
            }
        }
    };
    DisplayField::new(contribution, display)
}

pub(crate) fn format_discrete_bin_prefix(bin_index: usize, bin_count: usize) -> String {
    let width = if bin_count < 100 {
        4
    } else if bin_count < 1_000 {
        5
    } else {
        6
    };
    format!("{:<width$}", format!("#{bin_index}:"))
}

fn contribution_header(
    integration_state: &IntegrationState,
    discrete_monitoring_enabled: bool,
) -> StyledText {
    if discrete_monitoring_enabled
        && let Some(label) = integration_state
            .first_non_trivial_discrete_label
            .as_deref()
    {
        let mut header = StyledText::styled("Contribution", TextStyle::blue().bold());
        header.push_text("\n(idx=", TextStyle::PLAIN);
        header.push_text(label, TextStyle::blue().bold());
        header.push_text(")", TextStyle::PLAIN);
        return header;
    }

    styled_colored("Contribution", TextStyle::blue().bold())
}

fn header_left(
    elapsed_time: Duration,
    iter: usize,
    live_progress: Option<LiveIterationProgress>,
) -> StyledText {
    let mut text = StyledText::plain("[ ");
    text.push_text(
        format!(
            "{:^7}",
            utils::format_wdhms(elapsed_time.as_secs() as usize)
        ),
        TextStyle::PLAIN.bold(),
    );
    text.push_text(" ] ", TextStyle::PLAIN);
    let iteration_label = if live_progress.is_some() {
        format!("Iteration #{:-4} ( running   )", iter)
    } else {
        format!("Iteration #{:-4} ( completed )", iter)
    };
    text.push_text(iteration_label, TextStyle::green().bold());
    text
}

fn header_middle(
    cur_points: usize,
    total_points: usize,
    live_progress: Option<LiveIterationProgress>,
) -> StyledText {
    let mut text = StyledText::new();
    if let Some(progress) = live_progress {
        text.push_text(
            format!(
                "Iteration progress {}/{} ",
                format_iteration_points(progress.completed_points),
                format_iteration_points(progress.target_points),
            ),
            TextStyle::PLAIN,
        );
        let percentage = if progress.target_points == 0 {
            String::from("0.0%")
        } else {
            format!(
                "{:.1}%",
                (progress.completed_points as f64) / (progress.target_points as f64) * 100.0
            )
        };
        text.push_text(percentage, TextStyle::green());
        text.push_text(" ", TextStyle::PLAIN);
    } else if cur_points > 0 {
        text.push_text(
            format!(
                "# samples per iteration = {} ",
                format_iteration_points(cur_points)
            ),
            TextStyle::blue().bold(),
        );
    }
    text.push_text(
        format!("# samples total = {}", format_total_points(total_points)),
        TextStyle::green().bold(),
    );
    text
}

fn header_tail(cores: usize, elapsed_time: Duration, n_samples_evaluated: usize) -> StyledText {
    let mut text = StyledText::new();
    if n_samples_evaluated == 0 {
        text.push_text("N/A /sample/core", TextStyle::red());
    } else {
        text.push_text(
            format!(
                "{} /sample/core",
                utils::format_evaluation_time_from_f64(
                    elapsed_time.as_secs_f64() / (n_samples_evaluated as f64) * (cores as f64),
                )
            ),
            TextStyle::green().bold(),
        );
    }
    text.push_text(" ", TextStyle::PLAIN);
    text.push_text(format!("({cores} cores)"), TextStyle::blue().bold());
    text
}

fn main_results_row(
    integration_state: &IntegrationState,
    targets: &[Option<Complex<F<f64>>>],
    monitored_path: Option<&[usize]>,
    contribution: ContributionKind,
    component: ComponentKind,
    has_discrete_columns: bool,
) -> Option<MainResultsRow> {
    let slot_cells = integration_state
        .slot_metas
        .iter()
        .enumerate()
        .map(|(slot_index, _)| match contribution {
            ContributionKind::All => {
                let accumulator =
                    component_accumulator(&integration_state.all_integrals[slot_index], component);
                Some(MainTableSlotCells {
                    value: Some(format_value_field(accumulator.avg, accumulator.err)),
                    relative_error: format_relative_error_field(accumulator),
                    chi_sq: format_chi_sq_field(accumulator, integration_state.iter),
                    max_weight_impact: format_mwi_field(accumulator),
                    sample_fraction: None,
                    sample_count: None,
                    target_pdf: None,
                })
            }
            ContributionKind::Sum => {
                let summary = slot_component_summary(integration_state, slot_index, component)
                    .and_then(|summary| {
                        monitored_path.and_then(|path| summary_at_path(summary, path))
                    })?;
                let (avg, err) = sum_estimate_error(summary);
                Some(MainTableSlotCells {
                    value: Some(format_value_field(avg, err)),
                    relative_error: format_relative_error_field_from_estimate(avg, err),
                    chi_sq: None,
                    max_weight_impact: None,
                    sample_fraction: None,
                    sample_count: None,
                    target_pdf: None,
                })
            }
            ContributionKind::Bin(bin_index) => {
                let slot_context =
                    integration_state.monitored_discrete_context_for_slot(slot_index);
                let summary = slot_component_summary(integration_state, slot_index, component)
                    .and_then(|summary| {
                        monitored_path.and_then(|path| summary_at_path(summary, path))
                    })?;
                let bin = summary.bins.get(bin_index)?;
                let total_samples = total_processed_samples(summary);
                let sample_fraction = if has_discrete_columns && total_samples > 0 {
                    let raw =
                        bin.accumulator.processed_samples as f64 / total_samples as f64 * 100.0;
                    Some(DisplayField::new(
                        raw,
                        styled_colored(
                            format_significant_percentage(raw, 3, None),
                            TextStyle::blue(),
                        ),
                    ))
                } else {
                    None
                };
                let sample_count = if has_discrete_columns {
                    Some(DisplayField::new(
                        bin.accumulator.processed_samples,
                        styled_colored(
                            format_abbreviated_count(bin.accumulator.processed_samples),
                            TextStyle::blue(),
                        ),
                    ))
                } else {
                    None
                };
                let target_pdf = if has_discrete_columns {
                    slot_context
                        .as_ref()
                        .and_then(|ctx| ctx.pdfs.get(bin_index).copied())
                        .map(|pdf| {
                            DisplayField::new(
                                pdf.0 * 100.0,
                                styled_plain(format_significant_percentage(pdf.0 * 100.0, 3, None)),
                            )
                        })
                } else {
                    None
                };
                Some(MainTableSlotCells {
                    value: Some(format_value_field(bin.accumulator.avg, bin.accumulator.err)),
                    relative_error: format_relative_error_field(&bin.accumulator),
                    chi_sq: format_chi_sq_field(&bin.accumulator, integration_state.iter),
                    max_weight_impact: format_mwi_field(&bin.accumulator),
                    sample_fraction,
                    sample_count,
                    target_pdf,
                })
            }
        })
        .collect::<Option<Vec<_>>>()?;

    let slot0_target = targets
        .first()
        .and_then(|target| target.as_ref())
        .map(|target| match component {
            ComponentKind::Real => target.re,
            ComponentKind::Imag => target.im,
        });
    let (chi_sq, delta_sigma, delta_percent, max_weight_impact) = match contribution {
        ContributionKind::All => {
            let accumulator = component_accumulator(&integration_state.all_integrals[0], component);
            let (delta_sigma, delta_percent) = format_delta_fields(accumulator, slot0_target);
            (
                format_chi_sq_field(accumulator, integration_state.iter),
                delta_sigma,
                delta_percent,
                format_mwi_field(accumulator),
            )
        }
        ContributionKind::Sum => {
            let summary =
                slot_component_summary(integration_state, 0, component).and_then(|summary| {
                    monitored_path.and_then(|path| summary_at_path(summary, path))
                })?;
            let (avg, err) = sum_estimate_error(summary);
            let (delta_sigma, delta_percent) =
                format_delta_fields_from_estimate(avg, err, slot0_target);
            (None, delta_sigma, delta_percent, None)
        }
        ContributionKind::Bin(bin_index) => {
            let summary =
                slot_component_summary(integration_state, 0, component).and_then(|summary| {
                    monitored_path.and_then(|path| summary_at_path(summary, path))
                })?;
            let bin = summary.bins.get(bin_index)?;
            (
                format_chi_sq_field(&bin.accumulator, integration_state.iter),
                None,
                None,
                format_mwi_field(&bin.accumulator),
            )
        }
    };

    Some(MainResultsRow {
        contribution: contribution_display_field(contribution, integration_state),
        component: component.display_field(),
        slot_cells,
        chi_sq,
        delta_sigma,
        delta_percent,
        max_weight_impact,
    })
}

fn discrete_sort_key(
    integration_state: &IntegrationState,
    monitored_path: &[usize],
    component: ComponentKind,
    bin_index: usize,
    sort_mode: ContributionSortMode,
) -> f64 {
    let Some(summary) = slot_component_summary(integration_state, 0, component)
        .and_then(|summary| summary_at_path(summary, monitored_path))
    else {
        return 0.0;
    };
    let Some(bin) = summary.bins.get(bin_index) else {
        return 0.0;
    };
    match sort_mode {
        ContributionSortMode::Integral => bin.accumulator.avg.abs().0,
        ContributionSortMode::Error => bin.accumulator.err.0,
        ContributionSortMode::Index => bin_index as f64,
    }
}

fn build_main_results_section(request: &StatusUpdateBuildRequest<'_>) -> MainResultsSection {
    let components = ComponentKind::all_for_display(request.render_options.phase_display);
    let discrete_context = request.integration_state.monitored_discrete_context();
    let monitored_path = request.integration_state.monitored_discrete_path.as_deref();
    let has_discrete_columns = monitored_path.is_some();
    let has_target_columns = request.targets.first().is_some_and(Option::is_some);

    let mut row_groups = vec![MainResultsRowGroup {
        kind: MainResultsRowGroupKind::All,
        rows: components
            .iter()
            .filter_map(|component| {
                main_results_row(
                    request.integration_state,
                    request.targets,
                    monitored_path,
                    ContributionKind::All,
                    *component,
                    has_discrete_columns,
                )
            })
            .collect_vec(),
    }];

    if let Some(discrete_context) = discrete_context.as_ref() {
        let sum_rows = components
            .iter()
            .filter_map(|component| {
                main_results_row(
                    request.integration_state,
                    request.targets,
                    monitored_path,
                    ContributionKind::Sum,
                    *component,
                    has_discrete_columns,
                )
            })
            .collect_vec();
        if !sum_rows.is_empty() {
            row_groups.push(MainResultsRowGroup {
                kind: MainResultsRowGroupKind::Sum,
                rows: sum_rows,
            });
        }

        let bin_count = discrete_context.pdfs.len();
        match request.render_options.contribution_sort {
            ContributionSortMode::Index => {
                let mut rows = Vec::new();
                for bin_index in 0..bin_count {
                    for component in &components {
                        if let Some(row) = main_results_row(
                            request.integration_state,
                            request.targets,
                            monitored_path,
                            ContributionKind::Bin(bin_index),
                            *component,
                            has_discrete_columns,
                        ) {
                            rows.push(row);
                        }
                    }
                }
                if !rows.is_empty() {
                    row_groups.push(MainResultsRowGroup {
                        kind: MainResultsRowGroupKind::Bins,
                        rows,
                    });
                }
            }
            ContributionSortMode::Integral | ContributionSortMode::Error => {
                for component in &components {
                    let mut bin_indices = (0..bin_count).collect_vec();
                    bin_indices.sort_by(|lhs, rhs| {
                        discrete_sort_key(
                            request.integration_state,
                            &discrete_context.path,
                            *component,
                            *rhs,
                            request.render_options.contribution_sort,
                        )
                        .partial_cmp(&discrete_sort_key(
                            request.integration_state,
                            &discrete_context.path,
                            *component,
                            *lhs,
                            request.render_options.contribution_sort,
                        ))
                        .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    let rows = bin_indices
                        .into_iter()
                        .filter_map(|bin_index| {
                            main_results_row(
                                request.integration_state,
                                request.targets,
                                monitored_path,
                                ContributionKind::Bin(bin_index),
                                *component,
                                has_discrete_columns,
                            )
                        })
                        .collect_vec();
                    if !rows.is_empty() {
                        row_groups.push(MainResultsRowGroup {
                            kind: MainResultsRowGroupKind::Bins,
                            rows,
                        });
                    }
                }
            }
        }
    }

    MainResultsSection {
        header_left: header_left(
            request.elapsed_time,
            request.integration_state.iter,
            request.live_progress,
        ),
        header_middle: header_middle(
            request.cur_points,
            request.total_points_display,
            request.live_progress,
        ),
        header_tail: header_tail(
            request.cores,
            request.elapsed_time,
            request.n_samples_evaluated,
        ),
        contribution_header: contribution_header(request.integration_state, has_discrete_columns),
        slot_headers: request
            .integration_state
            .slot_metas
            .iter()
            .map(|slot_meta| styled_colored(slot_meta.key(), TextStyle::blue().bold()))
            .collect(),
        has_discrete_columns,
        has_target_columns,
        row_groups,
    }
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

fn styled_component_sign(
    component: ComponentKind,
    sign: &'static str,
    positive: bool,
) -> DisplayField<(ComponentKind, bool)> {
    let mut display = StyledText::new();
    display.append(component.label_display());
    display.push_text(" [", TextStyle::PLAIN);
    display.push_text(sign, TextStyle::blue());
    display.push_text("]", TextStyle::PLAIN);
    DisplayField::new((component, positive), display)
}

fn build_max_weight_details_section(
    integration_state: &IntegrationState,
    render_options: &IntegrationStatusViewOptions,
) -> Option<MaxWeightDetailsSection> {
    let rows_by_slot = integration_state
        .slot_metas
        .iter()
        .enumerate()
        .zip(integration_state.all_integrals.iter())
        .filter_map(|((slot_index, slot_meta), integral)| {
            let rows = max_weight_row_descriptors(render_options.phase_display)
                .into_iter()
                .filter_map(|(component, sign, positive)| {
                    let accumulator = match component {
                        ComponentKind::Real => &integral.re,
                        ComponentKind::Imag => &integral.im,
                    };
                    let (value, coordinates) = max_eval_entry(accumulator, positive)?;
                    Some(MaxWeightDetailsRow {
                        slot: DisplayField::new(
                            slot_meta.key(),
                            styled_colored(slot_meta.key(), TextStyle::green()),
                        ),
                        component_sign: styled_component_sign(component, sign, positive),
                        max_eval: DisplayField::new(
                            value,
                            styled_plain(format!("{:+.16e}", value)),
                        ),
                        coordinates: DisplayField::new(
                            coordinates
                                .map(|sample| {
                                    format_max_eval_sample(
                                        sample,
                                        &integration_state
                                            .sampling_state_for_slot(slot_index)
                                            .discrete_axis_labels,
                                        &[],
                                    )
                                })
                                .unwrap_or_else(|| "N/A".to_string()),
                            styled_plain(
                                coordinates
                                    .map(|sample| {
                                        format_max_eval_sample(
                                            sample,
                                            &integration_state
                                                .sampling_state_for_slot(slot_index)
                                                .discrete_axis_labels,
                                            &[],
                                        )
                                    })
                                    .unwrap_or_else(|| "N/A".to_string()),
                            ),
                        ),
                    })
                })
                .collect_vec();
            (!rows.is_empty()).then_some(rows)
        })
        .collect_vec();

    if rows_by_slot.is_empty() {
        None
    } else {
        Some(MaxWeightDetailsSection { rows_by_slot })
    }
}

fn build_discrete_max_weight_details_section(
    integration_state: &IntegrationState,
    render_options: &IntegrationStatusViewOptions,
) -> Option<DiscreteMaxWeightDetailsSection> {
    let discrete_context = integration_state.monitored_discrete_context()?;
    let contributions = std::iter::once(ContributionKind::All)
        .chain((0..discrete_context.pdfs.len()).map(ContributionKind::Bin))
        .collect_vec();

    let mut row_groups = Vec::new();
    for contribution in contributions {
        let mut rows = Vec::new();
        for (component, sign, positive) in max_weight_row_descriptors(render_options.phase_display)
        {
            let slot_values = integration_state
                .slot_metas
                .iter()
                .enumerate()
                .map(|(slot_index, _)| {
                    let accumulator = match contribution {
                        ContributionKind::All => {
                            let integral = &integration_state.all_integrals[slot_index];
                            match component {
                                ComponentKind::Real => &integral.re,
                                ComponentKind::Imag => &integral.im,
                            }
                        }
                        ContributionKind::Bin(bin_index) => {
                            let summary =
                                slot_component_summary(integration_state, slot_index, component)
                                    .and_then(|summary| {
                                        summary_at_path(summary, &discrete_context.path)
                                    })?;
                            &summary.bins.get(bin_index)?.accumulator
                        }
                        ContributionKind::Sum => return None,
                    };
                    Some(max_eval_entry(accumulator, positive).map(|(value, _)| {
                        DisplayField::new(value, styled_plain(format!("{:+.16e}", value)))
                    }))
                })
                .collect::<Option<Vec<_>>>()?;

            let slot_coordinates = integration_state
                .slot_metas
                .iter()
                .enumerate()
                .filter_map(|(slot_index, slot_meta)| {
                    let coordinates = match contribution {
                        ContributionKind::All => {
                            let accumulator = match component {
                                ComponentKind::Real => {
                                    &integration_state.all_integrals[slot_index].re
                                }
                                ComponentKind::Imag => {
                                    &integration_state.all_integrals[slot_index].im
                                }
                            };
                            max_eval_entry(accumulator, positive)
                                .and_then(|(_, sample)| sample)
                                .map(|sample| {
                                    format_max_eval_sample(
                                        sample,
                                        &integration_state
                                            .sampling_state_for_slot(slot_index)
                                            .discrete_axis_labels,
                                        &[],
                                    )
                                })
                        }
                        ContributionKind::Bin(bin_index) => {
                            let slot_context = integration_state
                                .monitored_discrete_context_for_slot(slot_index)?;
                            let summary =
                                slot_component_summary(integration_state, slot_index, component)
                                    .and_then(|summary| {
                                        summary_at_path(summary, &slot_context.path)
                                    })?;
                            max_eval_entry(&summary.bins.get(bin_index)?.accumulator, positive)
                                .and_then(|(_, sample)| sample)
                                .map(|sample| {
                                    format_max_eval_sample(
                                        sample,
                                        &integration_state
                                            .sampling_state_for_slot(slot_index)
                                            .discrete_axis_labels,
                                        &slot_context.path,
                                    )
                                })
                        }
                        ContributionKind::Sum => None,
                    }?;

                    Some(SlotCoordinateEntry {
                        slot: DisplayField::new(
                            slot_meta.key(),
                            styled_colored(slot_meta.key(), TextStyle::green()),
                        ),
                        coordinates: DisplayField::new(
                            coordinates.clone(),
                            styled_plain(coordinates),
                        ),
                    })
                })
                .collect_vec();

            if slot_values.iter().all(|value| value.is_none()) && slot_coordinates.is_empty() {
                continue;
            }

            rows.push(DiscreteMaxWeightRow {
                contribution: contribution_display_field(contribution, integration_state),
                component_sign: styled_component_sign(component, sign, positive),
                slot_values,
                slot_coordinates,
            });
        }
        if !rows.is_empty() {
            row_groups.push(rows);
        }
    }

    if row_groups.is_empty() {
        return None;
    }

    Some(DiscreteMaxWeightDetailsSection {
        contribution_header: contribution_header(integration_state, true),
        slot_headers: integration_state
            .slot_metas
            .iter()
            .map(|slot_meta| styled_colored(slot_meta.key(), TextStyle::blue().bold()))
            .collect(),
        row_groups,
    })
}

pub(crate) fn build_status_update(request: StatusUpdateBuildRequest<'_>) -> StatusUpdate {
    let target_accuracy_status = evaluate_target_accuracy(
        request.integration_state,
        request.total_points_display,
        request.elapsed_time,
        request.targets,
        request.render_options.training_phase_display,
        request.render_options.target_relative_accuracy,
        request.render_options.target_absolute_accuracy,
    );
    StatusUpdate {
        kind: request.kind,
        meta: StatusMeta {
            elapsed_time: request.elapsed_time,
            iteration_elapsed_time: request.iteration_elapsed_time,
            iteration: request.integration_state.iter,
            current_iteration_points: request.cur_points,
            total_points: request.total_points_display,
            n_samples_evaluated: request.n_samples_evaluated,
            cores: request.cores,
            training_slot: request.render_options.training_slot,
            training_phase_display: request.render_options.training_phase_display,
            live_progress: request.live_progress,
            show_eta_to_target: request.render_options.target_relative_accuracy.is_some()
                || request.render_options.target_absolute_accuracy.is_some(),
            eta_to_target_specification: format_eta_to_target_specification(
                request.render_options.target_relative_accuracy,
                request.render_options.target_absolute_accuracy,
            ),
            eta_to_target: target_accuracy_status.eta_to_target,
        },
        targets: request.targets.to_vec(),
        slot_training_phase_displays: request.render_options.slot_training_phase_displays.clone(),
        per_slot_training_phase: request.render_options.per_slot_training_phase,
        main_results: build_main_results_section(&request),
        max_weight_details: build_max_weight_details_section(
            request.integration_state,
            request.render_options,
        ),
        discrete_max_weight_details: build_discrete_max_weight_details_section(
            request.integration_state,
            request.render_options,
        ),
        statistics: Some(StatisticsSection {
            global: request.integration_state.stats,
            slot_counters: request.integration_state.slot_stats.clone(),
            slot_labels: request
                .integration_state
                .slot_metas
                .iter()
                .map(|slot_meta| slot_meta.key())
                .collect(),
        }),
    }
}

pub(crate) fn build_saved_status_update(
    integration_state: &IntegrationState,
    targets: &[Option<Complex<F<f64>>>],
    render_options: &IntegrationStatusViewOptions,
) -> StatusUpdate {
    let final_render_options = render_options.clone().for_final();
    build_status_update(StatusUpdateBuildRequest {
        kind: IntegrationStatusKind::Final,
        integration_state,
        cores: integration_state.n_cores.max(1),
        elapsed_time: utils::duration_from_secs_f64_saturating(integration_state.elapsed_seconds),
        iteration_elapsed_time: Duration::ZERO,
        cur_points: 0,
        total_points_display: integration_state.num_points,
        n_samples_evaluated: integration_state.num_points,
        targets,
        render_options: &final_render_options,
        live_progress: None,
    })
}
