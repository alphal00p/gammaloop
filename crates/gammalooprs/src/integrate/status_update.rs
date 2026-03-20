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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct IntegrationStatusViewOptions {
    pub phase_display: IntegrationStatusPhaseDisplay,
    pub training_phase_display: IntegrationStatusPhaseDisplay,
    pub training_slot: usize,
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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
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
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub(crate) enum ComponentKind {
    Real,
    Imag,
}

impl ComponentKind {
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

    fn display_field(self) -> DisplayField<Self> {
        DisplayField::new(
            self,
            StyledText::styled(self.tag(), TextStyle::blue().bold()),
        )
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
        match self.meta.training_phase_display {
            IntegrationStatusPhaseDisplay::Real => Some(ComponentKind::Real),
            IntegrationStatusPhaseDisplay::Imag => Some(ComponentKind::Imag),
            IntegrationStatusPhaseDisplay::Both => None,
        }
    }

    pub(crate) fn statistics_snapshot(&self) -> Option<IntegrationStatisticsSnapshot> {
        self.statistics.map(|section| section.raw.snapshot())
    }

    pub(crate) fn training_target(&self) -> Option<F<f64>> {
        let component = self.training_component()?;
        let target = self.targets.get(self.meta.training_slot)?.as_ref()?;
        Some(match component {
            ComponentKind::Real => target.re,
            ComponentKind::Imag => target.im,
        })
    }
}

impl StatusMeta {
    pub(crate) fn iteration_progress_ratio(self) -> Option<f64> {
        self.live_progress.map(|progress| {
            if progress.target_points == 0 {
                0.0
            } else {
                progress.completed_points as f64 / progress.target_points as f64
            }
        })
    }

    pub(crate) fn iteration_eta(self) -> Option<Duration> {
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
}

impl MainResultsRow {
    pub(crate) fn slot_cell(&self, slot_index: usize) -> Option<&MainTableSlotCells> {
        self.slot_cells.get(slot_index)
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

#[derive(Clone, Copy, Debug)]
pub(crate) struct StatisticsSection {
    pub(crate) raw: StatisticsCounter,
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

fn styled_plain(text: impl Into<String>) -> StyledText {
    StyledText::plain(text)
}

fn styled_colored(text: impl Into<String>, style: TextStyle) -> StyledText {
    StyledText::styled(text, style)
}

fn format_value_field(avg: F<f64>, err: F<f64>) -> DisplayField<(F<f64>, F<f64>)> {
    DisplayField::new(
        (avg, err),
        styled_colored(
            format_signed_uncertainty(avg, err, UncertaintyNotation::Scientific),
            TextStyle::blue().bold(),
        ),
    )
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
        return styled_colored(
            format!("Contribution (idx={label})"),
            TextStyle::blue().bold(),
        );
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
    if n_samples_evaluated == 0 {
        return styled_colored("N/A /sample/core", TextStyle::red());
    }

    styled_colored(
        format!(
            "{} /sample/core",
            utils::format_evaluation_time_from_f64(
                elapsed_time.as_secs_f64() / (n_samples_evaluated as f64) * (cores as f64),
            )
        ),
        TextStyle::blue().bold(),
    )
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

fn build_main_results_section(
    kind: IntegrationStatusKind,
    integration_state: &IntegrationState,
    cores: usize,
    elapsed_time: Duration,
    cur_points: usize,
    total_points_display: usize,
    n_samples_evaluated: usize,
    targets: &[Option<Complex<F<f64>>>],
    render_options: &IntegrationStatusViewOptions,
    live_progress: Option<LiveIterationProgress>,
) -> MainResultsSection {
    let components = ComponentKind::all_for_display(render_options.phase_display);
    let discrete_context = integration_state.monitored_discrete_context();
    let monitored_path = integration_state.monitored_discrete_path.as_deref();
    let has_discrete_columns = monitored_path.is_some();
    let has_target_columns = targets.first().is_some_and(Option::is_some);

    let mut row_groups = vec![MainResultsRowGroup {
        kind: MainResultsRowGroupKind::All,
        rows: components
            .iter()
            .filter_map(|component| {
                main_results_row(
                    integration_state,
                    targets,
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
                    integration_state,
                    targets,
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
        match render_options.contribution_sort {
            ContributionSortMode::Index => {
                let mut rows = Vec::new();
                for bin_index in 0..bin_count {
                    for component in &components {
                        if let Some(row) = main_results_row(
                            integration_state,
                            targets,
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
                            integration_state,
                            &discrete_context.path,
                            *component,
                            *rhs,
                            render_options.contribution_sort,
                        )
                        .partial_cmp(&discrete_sort_key(
                            integration_state,
                            &discrete_context.path,
                            *component,
                            *lhs,
                            render_options.contribution_sort,
                        ))
                        .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    let rows = bin_indices
                        .into_iter()
                        .filter_map(|bin_index| {
                            main_results_row(
                                integration_state,
                                targets,
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

    let _ = kind;
    MainResultsSection {
        header_left: header_left(elapsed_time, integration_state.iter, live_progress),
        header_middle: header_middle(cur_points, total_points_display, live_progress),
        header_tail: header_tail(cores, elapsed_time, n_samples_evaluated),
        contribution_header: contribution_header(integration_state, has_discrete_columns),
        slot_headers: integration_state
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
    display.push_text(component.tag(), TextStyle::blue().bold());
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

pub(crate) fn build_status_update(
    kind: IntegrationStatusKind,
    integration_state: &IntegrationState,
    cores: usize,
    elapsed_time: Duration,
    iteration_elapsed_time: Duration,
    cur_points: usize,
    total_points_display: usize,
    n_samples_evaluated: usize,
    targets: &[Option<Complex<F<f64>>>],
    render_options: &IntegrationStatusViewOptions,
    live_progress: Option<LiveIterationProgress>,
) -> StatusUpdate {
    StatusUpdate {
        kind,
        meta: StatusMeta {
            elapsed_time,
            iteration_elapsed_time,
            iteration: integration_state.iter,
            current_iteration_points: cur_points,
            total_points: total_points_display,
            n_samples_evaluated,
            cores,
            training_slot: render_options.training_slot,
            training_phase_display: render_options.training_phase_display,
            live_progress,
        },
        targets: targets.to_vec(),
        main_results: build_main_results_section(
            kind,
            integration_state,
            cores,
            elapsed_time,
            cur_points,
            total_points_display,
            n_samples_evaluated,
            targets,
            render_options,
            live_progress,
        ),
        max_weight_details: build_max_weight_details_section(integration_state, render_options),
        discrete_max_weight_details: build_discrete_max_weight_details_section(
            integration_state,
            render_options,
        ),
        statistics: Some(StatisticsSection {
            raw: integration_state.stats,
        }),
    }
}

pub(crate) fn build_saved_status_update(
    integration_state: &IntegrationState,
    targets: &[Option<Complex<F<f64>>>],
    render_options: &IntegrationStatusViewOptions,
) -> StatusUpdate {
    build_status_update(
        IntegrationStatusKind::Final,
        integration_state,
        integration_state.n_cores.max(1),
        utils::duration_from_secs_f64_saturating(integration_state.elapsed_seconds),
        Duration::ZERO,
        0,
        integration_state.num_points,
        integration_state.num_points,
        targets,
        &render_options.for_final(),
        None,
    )
}
