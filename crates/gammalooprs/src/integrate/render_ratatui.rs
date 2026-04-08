use std::cmp::Ordering;
use std::time::Duration;

use ratatui::{
    Frame,
    layout::{Alignment, Constraint, Direction, Layout, Rect},
    prelude::{Color, Line, Modifier, Span, Style},
    symbols::Marker,
    text::Text,
    widgets::{
        Axis, Block, Borders, Cell, Chart, Clear, Dataset, Gauge, Paragraph, Row, Table,
        TableState, Tabs, Wrap,
    },
};

use crate::utils;

use super::{
    StatusUpdate,
    display::{StyledText, TextColor, TextStyle},
    status_update::{
        ComponentKind, ContributionKind, ContributionSortMode, MainResultsRow,
        MainResultsRowGroupKind, MainTableSlotCells, StatisticsMixSegment, StatisticsScope,
    },
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum DashboardTab {
    Overview,
    Discrete,
    MaxWeight,
}

impl DashboardTab {
    fn all() -> [Self; 3] {
        [Self::Overview, Self::Discrete, Self::MaxWeight]
    }

    fn index(self) -> usize {
        match self {
            Self::Overview => 0,
            Self::Discrete => 1,
            Self::MaxWeight => 2,
        }
    }

    fn title(self) -> &'static str {
        match self {
            Self::Overview => "Overview",
            Self::Discrete => "Discrete",
            Self::MaxWeight => "Max Weight",
        }
    }

    fn from_index(index: usize) -> Self {
        Self::all().get(index).copied().unwrap_or(Self::Overview)
    }

    fn next(self) -> Self {
        Self::from_index((self.index() + 1) % Self::all().len())
    }

    fn previous(self) -> Self {
        Self::from_index((self.index() + Self::all().len() - 1) % Self::all().len())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum DensityMode {
    Compact,
    Metrics,
    Full,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
enum ChartHistoryWindow {
    #[default]
    Full,
    RecentIterations(usize),
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
enum DashboardStatisticsScope {
    #[default]
    Global,
    FocusedSlot,
}

const MAX_HISTORY_POINTS: usize = 4096;

impl ChartHistoryWindow {
    fn description(self) -> String {
        match self {
            Self::Full => "full".to_string(),
            Self::RecentIterations(count) => format!("last {count} iters"),
        }
    }

    fn toggle(self) -> Self {
        match self {
            Self::Full => Self::RecentIterations(6),
            Self::RecentIterations(_) => Self::Full,
        }
    }

    fn widen(self) -> Self {
        match self {
            Self::Full => Self::RecentIterations(6),
            Self::RecentIterations(count) => Self::RecentIterations((count + 1).min(128)),
        }
    }

    fn narrow(self) -> Self {
        match self {
            Self::Full => Self::RecentIterations(6),
            Self::RecentIterations(count) => Self::RecentIterations(count.saturating_sub(1).max(1)),
        }
    }
}

impl DensityMode {
    fn next(self) -> Self {
        match self {
            Self::Compact => Self::Metrics,
            Self::Metrics => Self::Full,
            Self::Full => Self::Compact,
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::Compact => "compact",
            Self::Metrics => "metrics",
            Self::Full => "full",
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct SlotMetricVisibility {
    relative_error: bool,
    chi_sq: bool,
    max_weight_impact: bool,
}

impl Default for SlotMetricVisibility {
    fn default() -> Self {
        Self {
            relative_error: true,
            chi_sq: true,
            max_weight_impact: true,
        }
    }
}

#[derive(Clone, Debug)]
struct HistoryPoint {
    iteration: usize,
    samples: usize,
    real_slot_values: Vec<Option<(f64, f64)>>,
    imag_slot_values: Vec<Option<(f64, f64)>>,
    completed_iteration: bool,
}

impl HistoryPoint {
    fn slot_values(&self, component: ComponentKind) -> &[Option<(f64, f64)>] {
        match component {
            ComponentKind::Real => &self.real_slot_values,
            ComponentKind::Imag => &self.imag_slot_values,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct DiscreteRowRef<'a> {
    row: &'a MainResultsRow,
}

impl<'a> DiscreteRowRef<'a> {
    fn component(self) -> ComponentKind {
        self.row.component.raw
    }

    fn contribution(self) -> ContributionKind {
        self.row.contribution.raw
    }

    fn sort_value(self, mode: ContributionSortMode, slot_index: usize) -> f64 {
        let Some(cell) = self.row.slot_cell(slot_index) else {
            return 0.0;
        };
        match mode {
            ContributionSortMode::Index => match self.contribution() {
                ContributionKind::Bin(index) => index as f64,
                ContributionKind::All => -1.0,
                ContributionKind::Sum => -0.5,
            },
            ContributionSortMode::Integral => cell
                .value
                .as_ref()
                .map(|value| value.raw.0.abs().0)
                .unwrap_or_default(),
            ContributionSortMode::Error => cell
                .value
                .as_ref()
                .map(|value| value.raw.1.abs().0)
                .unwrap_or_default(),
        }
    }
}

pub struct RatatuiDashboardState {
    latest_update: Option<StatusUpdate>,
    active_tab: DashboardTab,
    density: DensityMode,
    metric_visibility: SlotMetricVisibility,
    focused_slot: usize,
    statistics_scope: DashboardStatisticsScope,
    selected_discrete_row: usize,
    discrete_sort: ContributionSortMode,
    discrete_descending: bool,
    show_help: bool,
    history: Vec<HistoryPoint>,
    chart_history_window: ChartHistoryWindow,
    chart_y_sigma_span: usize,
    chart_component: Option<ComponentKind>,
}

impl Default for RatatuiDashboardState {
    fn default() -> Self {
        Self {
            latest_update: None,
            active_tab: DashboardTab::Overview,
            density: DensityMode::Metrics,
            metric_visibility: SlotMetricVisibility::default(),
            focused_slot: 0,
            statistics_scope: DashboardStatisticsScope::Global,
            selected_discrete_row: 0,
            discrete_sort: ContributionSortMode::Error,
            discrete_descending: true,
            show_help: false,
            history: Vec::new(),
            chart_history_window: ChartHistoryWindow::default(),
            chart_y_sigma_span: 4,
            chart_component: None,
        }
    }
}

impl RatatuiDashboardState {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn has_update(&self) -> bool {
        self.latest_update.is_some()
    }

    pub fn update(&mut self, update: StatusUpdate) {
        self.push_history_point(&update);
        self.latest_update = Some(update);
        self.clamp_state();
    }

    pub fn next_tab(&mut self) {
        self.active_tab = self.active_tab.next();
    }

    pub fn previous_tab(&mut self) {
        self.active_tab = self.active_tab.previous();
    }

    pub fn select_tab(&mut self, index: usize) {
        self.active_tab = DashboardTab::from_index(index);
    }

    pub fn cycle_density(&mut self) {
        self.density = self.density.next();
    }

    pub fn toggle_help(&mut self) {
        self.show_help = !self.show_help;
    }

    pub fn toggle_relative_error(&mut self) {
        self.metric_visibility.relative_error = !self.metric_visibility.relative_error;
    }

    pub fn toggle_chi_sq(&mut self) {
        self.metric_visibility.chi_sq = !self.metric_visibility.chi_sq;
    }

    pub fn toggle_max_weight_impact(&mut self) {
        self.metric_visibility.max_weight_impact = !self.metric_visibility.max_weight_impact;
    }

    pub fn focus_next_slot(&mut self) {
        let slot_count = self
            .latest_update
            .as_ref()
            .map(|update| update.main_results.slot_headers.len())
            .unwrap_or_default();
        if slot_count > 0 {
            self.focused_slot = (self.focused_slot + 1) % slot_count;
        }
        self.clamp_state();
    }

    pub fn focus_previous_slot(&mut self) {
        let slot_count = self
            .latest_update
            .as_ref()
            .map(|update| update.main_results.slot_headers.len())
            .unwrap_or_default();
        if slot_count > 0 {
            self.focused_slot = (self.focused_slot + slot_count - 1) % slot_count;
        }
        self.clamp_state();
    }

    pub fn toggle_statistics_scope(&mut self) {
        self.statistics_scope = match self.statistics_scope {
            DashboardStatisticsScope::Global => DashboardStatisticsScope::FocusedSlot,
            DashboardStatisticsScope::FocusedSlot => DashboardStatisticsScope::Global,
        };
    }

    pub fn select_next_discrete_row(&mut self) {
        let row_count = self.discrete_rows().len();
        if row_count > 0 {
            self.selected_discrete_row = (self.selected_discrete_row + 1) % row_count;
        }
    }

    pub fn select_previous_discrete_row(&mut self) {
        let row_count = self.discrete_rows().len();
        if row_count > 0 {
            self.selected_discrete_row = (self.selected_discrete_row + row_count - 1) % row_count;
        }
    }

    pub fn cycle_discrete_sort(&mut self) {
        self.discrete_sort = match self.discrete_sort {
            ContributionSortMode::Index => ContributionSortMode::Integral,
            ContributionSortMode::Integral => ContributionSortMode::Error,
            ContributionSortMode::Error => ContributionSortMode::Index,
        };
        self.clamp_state();
    }

    pub fn toggle_discrete_sort_direction(&mut self) {
        self.discrete_descending = !self.discrete_descending;
    }

    pub fn toggle_chart_history_window(&mut self) {
        self.chart_history_window = self.chart_history_window.toggle();
    }

    pub fn widen_chart_history_window(&mut self) {
        self.chart_history_window = self.chart_history_window.widen();
    }

    pub fn narrow_chart_history_window(&mut self) {
        self.chart_history_window = self.chart_history_window.narrow();
    }

    pub fn widen_chart_y_sigma_span(&mut self) {
        self.chart_y_sigma_span = (self.chart_y_sigma_span + 1).min(128);
    }

    pub fn narrow_chart_y_sigma_span(&mut self) {
        self.chart_y_sigma_span = self.chart_y_sigma_span.saturating_sub(1).max(1);
    }

    pub fn reset_chart_y_sigma_span(&mut self) {
        self.chart_y_sigma_span = 4;
    }

    pub fn toggle_chart_component(&mut self) {
        let Some(update) = self.latest_update.as_ref() else {
            return;
        };
        let available = self.available_chart_components(update);
        if available.len() < 2 {
            return;
        }

        let current = self.chart_component(update);
        self.chart_component = match current {
            Some(ComponentKind::Real) if available.contains(&ComponentKind::Imag) => {
                Some(ComponentKind::Imag)
            }
            Some(ComponentKind::Imag) if available.contains(&ComponentKind::Real) => {
                Some(ComponentKind::Real)
            }
            _ => available.first().copied(),
        };
    }

    pub fn draw(&self, frame: &mut Frame<'_>) {
        let area = frame.area();
        if let Some(update) = self.latest_update.as_ref() {
            self.draw_dashboard(frame, area, update);
            if self.show_help {
                self.draw_help_overlay(frame, area);
            }
        } else {
            frame.render_widget(
                Paragraph::new("Waiting for integration status updates...")
                    .block(titled_block("Integration dashboard")),
                area,
            );
        }
    }

    fn clamp_state(&mut self) {
        if let Some(update) = self.latest_update.as_ref() {
            let slot_count = update.main_results.slot_headers.len();
            if slot_count == 0 {
                self.focused_slot = 0;
            } else {
                self.focused_slot = self.focused_slot.min(slot_count - 1);
            }
        } else {
            self.focused_slot = 0;
        }

        if let Some(update) = self.latest_update.as_ref() {
            let available = self.available_chart_components(update);
            if available.is_empty() {
                self.chart_component = None;
            } else if !self
                .chart_component
                .is_some_and(|component| available.contains(&component))
            {
                self.chart_component = update
                    .training_component()
                    .filter(|component| available.contains(component))
                    .or_else(|| available.first().copied());
            }
        }

        let row_count = self.discrete_rows().len();
        if row_count == 0 {
            self.selected_discrete_row = 0;
        } else {
            self.selected_discrete_row = self.selected_discrete_row.min(row_count - 1);
        }
    }

    fn push_history_point(&mut self, update: &StatusUpdate) {
        let real_slot_values = self.collect_history_slot_values(update, ComponentKind::Real);
        let imag_slot_values = self.collect_history_slot_values(update, ComponentKind::Imag);
        if real_slot_values.iter().all(Option::is_none)
            && imag_slot_values.iter().all(Option::is_none)
        {
            return;
        }

        let point = HistoryPoint {
            iteration: update.meta.iteration,
            samples: update.meta.total_points,
            real_slot_values,
            imag_slot_values,
            completed_iteration: !matches!(update.kind(), super::IntegrationStatusKind::Live),
        };

        if self
            .history
            .last()
            .is_some_and(|previous| previous.samples == point.samples)
        {
            let preserve_completed_point = self.history.last().is_some_and(|previous| {
                previous.completed_iteration
                    && !point.completed_iteration
                    && update.is_initial_live_status()
            });
            if !preserve_completed_point {
                let _ = self.history.pop();
            }
        }
        if self
            .history
            .last()
            .is_some_and(|previous| previous.samples > point.samples)
        {
            self.history.clear();
        }
        self.history.push(point);
        while self.history.len() > MAX_HISTORY_POINTS {
            self.compact_history_preserving_span();
        }
    }

    fn compact_history_preserving_span(&mut self) {
        if self.history.len() <= MAX_HISTORY_POINTS {
            return;
        }

        let retain_recent = MAX_HISTORY_POINTS / 2;
        let split_index = self.history.len().saturating_sub(retain_recent);
        if split_index == 0 {
            return;
        }

        let older = &self.history[..split_index];
        let recent = &self.history[split_index..];
        let max_older_points = MAX_HISTORY_POINTS.saturating_sub(recent.len()).max(2);
        let step = older.len().div_ceil(max_older_points);

        let mut compacted = Vec::with_capacity(MAX_HISTORY_POINTS);
        for (index, point) in older.iter().enumerate() {
            if index == 0 || index + 1 == older.len() || index % step == 0 {
                compacted.push(point.clone());
            }
        }
        compacted.extend(recent.iter().cloned());

        if compacted.len() > MAX_HISTORY_POINTS {
            let keep_from_older = compacted.len().saturating_sub(recent.len());
            let extra = compacted.len() - MAX_HISTORY_POINTS;
            if keep_from_older > extra + 1 {
                compacted.drain(1..=extra);
            } else {
                compacted.truncate(MAX_HISTORY_POINTS);
            }
        }

        self.history = compacted;
    }

    fn collect_history_slot_values(
        &self,
        update: &StatusUpdate,
        component: ComponentKind,
    ) -> Vec<Option<(f64, f64)>> {
        let Some(row) = update
            .main_results
            .find_row(ContributionKind::All, component)
        else {
            return vec![None; update.main_results.slot_headers.len()];
        };
        (0..update.main_results.slot_headers.len())
            .map(|slot_index| {
                row.slot_cell(slot_index)
                    .and_then(|cell| cell.value.as_ref())
                    .map(|value| (value.raw.0.0, value.raw.1.0.abs()))
            })
            .collect()
    }

    fn visible_history(&self) -> Vec<&HistoryPoint> {
        match self.chart_history_window {
            ChartHistoryWindow::Full => self.history.iter().collect(),
            ChartHistoryWindow::RecentIterations(iterations) => {
                if iterations == 0 || self.history.is_empty() {
                    return self.history.iter().collect();
                }
                let latest_iteration = self
                    .history
                    .last()
                    .map(|point| point.iteration)
                    .unwrap_or_default();
                let start_iteration = latest_iteration
                    .saturating_add(1)
                    .saturating_sub(iterations);
                self.history
                    .iter()
                    .filter(|point| point.iteration >= start_iteration)
                    .collect()
            }
        }
    }

    #[cfg(test)]
    pub(crate) fn visible_history_sample_bounds(&self) -> Option<(usize, usize)> {
        let visible = self.visible_history();
        Some((visible.first()?.samples, visible.last()?.samples))
    }

    fn available_chart_components(&self, update: &StatusUpdate) -> Vec<ComponentKind> {
        [ComponentKind::Real, ComponentKind::Imag]
            .into_iter()
            .filter(|component| {
                update
                    .main_results
                    .find_row(ContributionKind::All, *component)
                    .and_then(|row| row.slot_cell(self.focused_slot))
                    .and_then(|cell| cell.value.as_ref())
                    .is_some()
            })
            .collect()
    }

    fn chart_component(&self, update: &StatusUpdate) -> Option<ComponentKind> {
        let available = self.available_chart_components(update);
        if available.is_empty() {
            return None;
        }

        self.chart_component
            .filter(|component| available.contains(component))
            .or_else(|| {
                update
                    .training_component_for_slot(self.focused_slot)
                    .filter(|component| available.contains(component))
            })
            .or_else(|| available.first().copied())
    }

    fn draw_dashboard(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let vertical = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),
                Constraint::Length(6),
                Constraint::Length(update.main_results.slot_headers.len().max(1) as u16 + 2),
                Constraint::Min(12),
                Constraint::Length(2),
            ])
            .split(area);

        self.draw_tabs(frame, vertical[0]);
        self.draw_progress(frame, vertical[1], update);
        self.draw_slot_ribbon(frame, vertical[2], update);

        match self.active_tab {
            DashboardTab::Overview => self.draw_overview_tab(frame, vertical[3], update),
            DashboardTab::Discrete => self.draw_discrete_tab(frame, vertical[3], update),
            DashboardTab::MaxWeight => self.draw_max_weight_tab(frame, vertical[3], update),
        }

        self.draw_footer(frame, vertical[4], update);
    }

    fn draw_tabs(&self, frame: &mut Frame<'_>, area: Rect) {
        let titles = DashboardTab::all()
            .into_iter()
            .map(|tab| {
                Line::from(vec![Span::styled(
                    format!(" {} ", tab.title()),
                    Style::default().add_modifier(Modifier::BOLD),
                )])
            })
            .collect::<Vec<_>>();

        let tabs = Tabs::new(titles)
            .block(titled_block("GammaLoop status"))
            .select(self.active_tab.index())
            .highlight_style(Style::default().add_modifier(Modifier::REVERSED | Modifier::BOLD));
        frame.render_widget(tabs, area);
    }

    fn draw_progress(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let block = titled_block("Iteration progress");
        let inner = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(1),
                Constraint::Length(1),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(block.inner(area));
        frame.render_widget(block, area);

        frame.render_widget(Paragraph::new(""), inner[0]);

        let mut left = StyledText::new();
        left.push_text(" ", TextStyle::PLAIN);
        left.push_text(
            format!(
                "[ {:^7} ]",
                utils::format_wdhms(update.meta.elapsed_time.as_secs() as usize)
            ),
            TextStyle::PLAIN.bold(),
        );
        left.push_text(" | ", TextStyle::PLAIN);
        left.push_text("# samples total ", TextStyle::green().bold());
        left.push_text(
            super::display::format_total_points(update.meta.total_points),
            TextStyle::PLAIN.bold(),
        );
        left.push_text(" | ", TextStyle::PLAIN);
        left.push_text("Iteration ", TextStyle::green().bold());
        left.push_text(
            format!("#{}", update.meta.iteration),
            TextStyle::PLAIN.bold(),
        );
        left.push_text(
            format!(
                " ({})",
                if matches!(update.kind(), super::IntegrationStatusKind::Live) {
                    "running"
                } else {
                    "completed"
                }
            ),
            TextStyle::green().bold(),
        );
        if matches!(update.kind(), super::IntegrationStatusKind::Live) {
            left.push_text(" | ", TextStyle::PLAIN);
            left.push_text("ETA ", TextStyle::green().bold());
            left.push_text(
                update
                    .meta
                    .iteration_eta()
                    .map(|duration| utils::format_wdhms(duration.as_secs() as usize))
                    .unwrap_or_else(|| "warming up".to_string()),
                TextStyle::PLAIN,
            );
            if update.meta.show_eta_to_target {
                left.push_text(" | ", TextStyle::PLAIN);
                left.push_text("ETA to target ", TextStyle::green().bold());
                if let Some(specification) = update.meta.eta_to_target_specification() {
                    left.push_text("(", TextStyle::green().bold());
                    left.push_text(specification, TextStyle::green().bold());
                    left.push_text(") ", TextStyle::green().bold());
                }
                left.push_text(
                    update
                        .meta
                        .eta_to_target()
                        .map(|duration| {
                            if duration == Duration::MAX {
                                "∞".to_string()
                            } else {
                                utils::format_wdhms(duration.as_secs() as usize)
                            }
                        })
                        .unwrap_or_else(|| "N/A".to_string()),
                    TextStyle::PLAIN,
                );
            }
        }

        let mut right = StyledText::new();
        right.push_text(
            update
                .meta
                .total_sample_rate_per_second()
                .map(format_samples_per_second)
                .unwrap_or_else(|| "N/A".to_string()),
            TextStyle::PLAIN.bold(),
        );
        right.push_text(" ", TextStyle::PLAIN);
        right.push_text("#samples/s", TextStyle::green().bold());
        right.push_text(" | ", TextStyle::PLAIN);
        if let Some(sample_core_time) = update.meta.sample_core_time() {
            right.push_text(sample_core_time, TextStyle::PLAIN.bold());
        } else {
            right.push_text("N/A", TextStyle::red());
        }
        right.push_text(" /sample/core", TextStyle::green().bold());
        right.push_text(" ", TextStyle::PLAIN);

        let header = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Percentage(68), Constraint::Percentage(32)])
            .split(inner[1]);
        frame.render_widget(Paragraph::new(text_from_styled_text(&left)), header[0]);
        frame.render_widget(
            Paragraph::new(text_from_styled_text(&right)).alignment(Alignment::Right),
            header[1],
        );

        let progress = update.meta.iteration_progress_ratio().unwrap_or_else(|| {
            if matches!(update.kind(), super::IntegrationStatusKind::Live) {
                0.0
            } else {
                1.0
            }
        });
        let gauge = Gauge::default()
            .gauge_style(
                Style::default()
                    .fg(Color::Green)
                    .bg(Color::Rgb(70, 70, 70))
                    .add_modifier(Modifier::BOLD),
            )
            .ratio(progress.clamp(0.0, 1.0))
            .label(Span::styled(
                format!(
                    "{} / {} ({:.1}%)",
                    abbreviate_count(
                        update
                            .meta
                            .live_progress
                            .map(|progress| progress.completed_points)
                            .unwrap_or(update.meta.current_iteration_points),
                    ),
                    abbreviate_count(
                        update
                            .meta
                            .live_progress
                            .map(|progress| progress.target_points)
                            .unwrap_or(update.meta.current_iteration_points),
                    ),
                    progress * 100.0
                ),
                Style::default()
                    .fg(Color::White)
                    .add_modifier(Modifier::BOLD),
            ));
        frame.render_widget(gauge, inner[2]);
        frame.render_widget(Paragraph::new(""), inner[3]);
    }

    fn draw_slot_ribbon(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let slot_name_width = update
            .main_results
            .slot_headers
            .iter()
            .map(|header| header.to_plain_string().chars().count())
            .max()
            .unwrap_or(10)
            .max(10) as u16;
        let value_width = (0..update.main_results.slot_headers.len())
            .flat_map(|slot_index| {
                [ComponentKind::Real, ComponentKind::Imag]
                    .into_iter()
                    .filter_map(move |component| {
                        update
                            .main_results
                            .find_row(ContributionKind::All, component)
                            .and_then(|row| row.slot_cell(slot_index))
                            .and_then(|cell| cell.value.as_ref())
                            .map(|value| value.display.to_plain_string().chars().count())
                    })
            })
            .max()
            .unwrap_or(16)
            .max(16) as u16;
        let has_targets = (0..update.main_results.slot_headers.len()).any(|slot_index| {
            update
                .target_display_for_slot_component(slot_index, ComponentKind::Real)
                .is_some()
                || update
                    .target_display_for_slot_component(slot_index, ComponentKind::Imag)
                    .is_some()
        });
        let target_value_width = (0..update.main_results.slot_headers.len())
            .flat_map(|slot_index| {
                [ComponentKind::Real, ComponentKind::Imag]
                    .into_iter()
                    .filter_map(move |component| {
                        update
                            .target_display_for_slot_component(slot_index, component)
                            .map(|value| value.to_plain_string().chars().count())
                    })
            })
            .max()
            .unwrap_or(16)
            .max(16) as u16;
        let rows = (0..update.main_results.slot_headers.len())
            .map(|slot_index| {
                let re_value = update
                    .main_results
                    .find_row(ContributionKind::All, ComponentKind::Real)
                    .and_then(|row| row.slot_cell(slot_index))
                    .and_then(|cell| cell.value.as_ref())
                    .map(|value| &value.display);
                let im_value = update
                    .main_results
                    .find_row(ContributionKind::All, ComponentKind::Imag)
                    .and_then(|row| row.slot_cell(slot_index))
                    .and_then(|cell| cell.value.as_ref())
                    .map(|value| &value.display);
                let target_re =
                    update.target_display_for_slot_component(slot_index, ComponentKind::Real);
                let target_im =
                    update.target_display_for_slot_component(slot_index, ComponentKind::Imag);
                let mut cells = vec![
                    Cell::from(Span::styled(
                        if slot_index == self.focused_slot {
                            ">"
                        } else {
                            " "
                        },
                        if slot_index == self.focused_slot {
                            Style::default()
                                .fg(Color::White)
                                .add_modifier(Modifier::BOLD)
                        } else {
                            Style::default()
                        },
                    )),
                    cell_from_styled_text(&update.main_results.slot_headers[slot_index]),
                    cell_from_styled_text(&component_label_with_colon(ComponentKind::Real)),
                    cell_from_optional_styled_text(re_value),
                    cell_from_styled_text(&component_label_with_colon(ComponentKind::Imag)),
                    cell_from_optional_styled_text(im_value),
                ];
                if has_targets {
                    cells.push(cell_from_styled_text(&StyledText::plain("|")));
                    cells.push(plain_header_cell("trgt"));
                    cells.push(cell_from_styled_text(&component_label_with_colon(
                        ComponentKind::Real,
                    )));
                    cells.push(cell_from_optional_styled_text(target_re.as_ref()));
                    cells.push(cell_from_styled_text(&component_label_with_colon(
                        ComponentKind::Imag,
                    )));
                    cells.push(cell_from_optional_styled_text(target_im.as_ref()));
                }
                Row::new(cells)
            })
            .collect::<Vec<_>>();

        let mut widths = vec![
            Constraint::Length(1),
            Constraint::Length(slot_name_width),
            Constraint::Length(3),
            Constraint::Length(value_width),
            Constraint::Length(3),
            Constraint::Length(value_width),
        ];
        if has_targets {
            widths.extend([
                Constraint::Length(1),
                Constraint::Length(5),
                Constraint::Length(3),
                Constraint::Length(target_value_width),
                Constraint::Length(3),
                Constraint::Min(target_value_width),
            ]);
        }

        frame.render_widget(
            Table::new(rows, widths)
                .block(titled_block("Integrands"))
                .column_spacing(1),
            area,
        );
    }

    fn draw_overview_tab(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let vertical = if area.width >= 120 {
            Layout::default()
                .direction(Direction::Vertical)
                .constraints([Constraint::Percentage(66), Constraint::Percentage(34)])
                .split(area)
        } else {
            Layout::default()
                .direction(Direction::Vertical)
                .constraints([
                    Constraint::Percentage(48),
                    Constraint::Length(11),
                    Constraint::Min(10),
                ])
                .split(area)
        };

        if area.width >= 120 {
            let top = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([Constraint::Percentage(62), Constraint::Percentage(38)])
                .split(vertical[0]);
            self.draw_chart(frame, top[0], update);
            self.draw_focused_slot_detail(frame, top[1], update);

            let bottom = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([Constraint::Percentage(56), Constraint::Percentage(44)])
                .split(vertical[1]);
            self.draw_results_summary_table(frame, bottom[0], update);
            self.draw_statistics_panel(frame, bottom[1], update);
        } else {
            self.draw_chart(frame, vertical[0], update);
            self.draw_focused_slot_detail(frame, vertical[1], update);
            let bottom = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([Constraint::Percentage(55), Constraint::Percentage(45)])
                .split(vertical[2]);
            self.draw_results_summary_table(frame, bottom[0], update);
            self.draw_statistics_panel(frame, bottom[1], update);
        }
    }

    fn draw_chart(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let Some(component) = self.chart_component(update) else {
            frame.render_widget(
                Paragraph::new("Training phase is not a single component.")
                    .block(titled_block("Convergence")),
                area,
            );
            return;
        };

        let focused_slot = self
            .focused_slot
            .min(update.main_results.slot_headers.len().saturating_sub(1));
        let visible_history = self.visible_history();
        if visible_history.is_empty() {
            frame.render_widget(
                Paragraph::new("Waiting for history points...").block(titled_block("Convergence")),
                area,
            );
            return;
        }

        let central_points = visible_history
            .iter()
            .filter_map(|point| {
                (*point)
                    .slot_values(component)
                    .get(focused_slot)
                    .copied()
                    .flatten()
                    .map(|(central, _)| (point.samples as f64, central))
            })
            .collect::<Vec<_>>();
        let upper_points = visible_history
            .iter()
            .filter_map(|point| {
                (*point)
                    .slot_values(component)
                    .get(focused_slot)
                    .copied()
                    .flatten()
                    .map(|(central, error)| (point.samples as f64, central + error))
            })
            .collect::<Vec<_>>();
        let lower_points = visible_history
            .iter()
            .filter_map(|point| {
                (*point)
                    .slot_values(component)
                    .get(focused_slot)
                    .copied()
                    .flatten()
                    .map(|(central, error)| (point.samples as f64, central - error))
            })
            .collect::<Vec<_>>();
        let completed_points = visible_history
            .iter()
            .filter(|point| point.completed_iteration)
            .filter_map(|point| {
                (*point)
                    .slot_values(component)
                    .get(focused_slot)
                    .copied()
                    .flatten()
                    .map(|(central, _)| (point.samples as f64, central))
            })
            .collect::<Vec<_>>();
        if central_points.is_empty() {
            frame.render_widget(
                Paragraph::new("Waiting for focused-integrand history...").block(titled_block(
                    format!("Convergence ({})", component.phase_name()),
                )),
                area,
            );
            return;
        }

        let x_min = visible_history
            .first()
            .map(|point| point.samples as f64)
            .unwrap_or(0.0);
        let x_max = visible_history
            .last()
            .map(|point| point.samples as f64)
            .unwrap_or(1.0);
        let current_value = update
            .main_results
            .find_row(ContributionKind::All, component)
            .and_then(|row| row.slot_cell(focused_slot))
            .and_then(|cell| cell.value.as_ref())
            .map(|value| (value.raw.0.0, value.raw.1.0.abs()))
            .or_else(|| {
                self.history.last().and_then(|point| {
                    point
                        .slot_values(component)
                        .get(focused_slot)
                        .copied()
                        .flatten()
                })
            })
            .unwrap_or((0.0, 1.0));
        let sigma = current_value
            .1
            .max(current_value.0.abs().max(1.0) * 1.0e-12);
        let y_span = self.chart_y_sigma_span as f64;
        let y_min = current_value.0 - y_span * sigma;
        let y_max = current_value.0 + y_span * sigma;

        let mut datasets = vec![
            Dataset::default()
                .name("central")
                .marker(Marker::Braille)
                .graph_type(ratatui::widgets::GraphType::Line)
                .style(
                    Style::default()
                        .fg(Color::Green)
                        .add_modifier(Modifier::BOLD),
                )
                .data(&central_points),
            Dataset::default()
                .name("upper")
                .graph_type(ratatui::widgets::GraphType::Line)
                .style(Style::default().fg(Color::Blue))
                .data(&upper_points),
            Dataset::default()
                .name("lower")
                .graph_type(ratatui::widgets::GraphType::Line)
                .style(Style::default().fg(Color::Blue))
                .data(&lower_points),
            Dataset::default()
                .name("iter")
                .marker(Marker::Bar)
                .graph_type(ratatui::widgets::GraphType::Scatter)
                .style(Style::default().fg(Color::LightMagenta))
                .data(&completed_points),
        ];

        let target_points = update
            .target_for_slot_component(focused_slot, component)
            .map(|target| vec![(x_min, target.0), (x_max.max(x_min + 1.0), target.0)]);
        if let Some(target_points) = target_points.as_ref() {
            datasets.push(
                Dataset::default()
                    .name("target")
                    .graph_type(ratatui::widgets::GraphType::Line)
                    .style(Style::default().fg(Color::Yellow))
                    .data(target_points),
            );
        }

        let selected_for_training =
            update.slot_component_selected_for_training(focused_slot, component);
        let mut title = StyledText::styled("Convergence : ", TextStyle::green().bold());
        title.push_text(component.phase_name(), component.text_style());
        title.push_text(
            if selected_for_training {
                " (selected for training)"
            } else {
                " (not selected for training)"
            },
            TextStyle::green().bold(),
        );
        title.push_text(" [", TextStyle::green().bold());
        title.append(update.main_results.slot_headers[focused_slot].clone());
        title.push_text(
            format!(
                " | {} | y ±{}σ]",
                self.chart_history_window.description(),
                self.chart_y_sigma_span
            ),
            TextStyle::green().bold(),
        );
        let x_bounds = [x_min, x_max.max(x_min + 1.0)];
        let chart = Chart::new(datasets)
            .block(titled_block_styled(&title))
            .x_axis(
                Axis::default()
                    .title("samples")
                    .bounds(x_bounds)
                    .labels(count_axis_labels(x_bounds[0], x_bounds[1])),
            )
            .y_axis(
                Axis::default()
                    .title("central value")
                    .bounds([y_min, y_max])
                    .labels(value_axis_labels(y_min, y_max)),
            );
        frame.render_widget(chart, area);
    }

    fn draw_focused_slot_detail(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let focused_slot = self
            .focused_slot
            .min(update.main_results.slot_headers.len().saturating_sub(1));
        let rows = self.overview_summary_rows(update);
        let show_target_columns = rows.iter().any(|row| {
            update
                .target_deltas_for_row_slot(row, focused_slot)
                .0
                .is_some()
        });
        let focused_column_count = 5 + if show_target_columns { 2 } else { 0 };

        let mut table_rows = vec![blank_table_row(focused_column_count)];
        let mut inserted_sum_gap = false;
        for row in rows {
            if matches!(row.contribution.raw, ContributionKind::Sum) && !inserted_sum_gap {
                table_rows.push(blank_table_row(focused_column_count));
                inserted_sum_gap = true;
            }
            let Some(cell) = row.slot_cell(focused_slot) else {
                continue;
            };
            let (delta_sigma, delta_percent) = update.target_deltas_for_row_slot(row, focused_slot);
            let mut label = row.contribution.display.clone();
            label.push_text(" ", TextStyle::PLAIN);
            label.append(row.component.display.clone());
            let mut cells = vec![
                cell_from_styled_text(&label),
                cell_from_optional_styled_text(cell.value.as_ref().map(|value| &value.display)),
                cell_from_optional_styled_text(
                    cell.relative_error.as_ref().map(|value| &value.display),
                ),
                cell_from_optional_styled_text(cell.chi_sq.as_ref().map(|value| &value.display)),
                cell_from_optional_styled_text(
                    cell.max_weight_impact.as_ref().map(|value| &value.display),
                ),
            ];
            if show_target_columns {
                cells.push(cell_from_optional_styled_text(
                    delta_sigma.as_ref().map(|value| &value.display),
                ));
                cells.push(cell_from_optional_styled_text(
                    delta_percent.as_ref().map(|value| &value.display),
                ));
            }
            table_rows.push(Row::new(cells));
        }

        let mut headers = vec![
            plain_header_cell(""),
            plain_header_cell("integral"),
            plain_header_cell("% err"),
            plain_header_cell("chi^2"),
            plain_header_cell("m.w.i"),
        ];
        if show_target_columns {
            headers.push(plain_header_cell("Δ [σ]"));
            headers.push(plain_header_cell("Δ [%]"));
        }

        let mut widths = vec![
            Constraint::Length(9),
            Constraint::Length(16),
            Constraint::Length(8),
            Constraint::Length(8),
            Constraint::Length(10),
        ];
        if show_target_columns {
            widths.push(Constraint::Length(8));
            widths.push(Constraint::Length(8));
        }

        let table = Table::new(table_rows, widths)
            .header(Row::new(headers))
            .block(titled_block_styled(&scoped_integrand_title(
                "Focused integrand ",
                &update.main_results.slot_headers[focused_slot],
                "",
            )))
            .column_spacing(1);
        frame.render_widget(table, area);
    }

    fn draw_results_summary_table(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let rows = self.overview_summary_rows(update);
        if rows.is_empty() {
            frame.render_widget(
                Paragraph::new("No overview summary rows are available.")
                    .block(titled_block("Results summary")),
                area,
            );
            return;
        }

        let contribution_width = update
            .main_results
            .row_groups
            .iter()
            .flat_map(|group| group.rows.iter())
            .map(|row| row.contribution.display.to_plain_string().chars().count())
            .chain(std::iter::once(
                update
                    .main_results
                    .contribution_header
                    .to_plain_string()
                    .chars()
                    .count(),
            ))
            .max()
            .unwrap_or(18)
            .max(18) as u16;
        let contribution_title = styled_text_line(&update.main_results.contribution_header, 0);
        let contribution_detail = styled_text_line(&update.main_results.contribution_header, 1);
        let mut header_top = vec![cell_from_styled_text(&contribution_title), blank_cell()];
        let mut header_bottom = vec![cell_from_styled_text(&contribution_detail), blank_cell()];
        let mut widths = vec![
            Constraint::Length(contribution_width),
            Constraint::Length(3),
        ];
        for (slot_index, slot_header) in update.main_results.slot_headers.iter().enumerate() {
            let integral_width = rows
                .iter()
                .filter_map(|row| {
                    row.slot_cell(slot_index)
                        .and_then(|cell| cell.value.as_ref())
                        .map(|value| value.display.to_plain_string().chars().count())
                })
                .chain(std::iter::once(
                    slot_header.to_plain_string().chars().count(),
                ))
                .chain(std::iter::once("integral".chars().count()))
                .max()
                .unwrap_or(16)
                .max(16) as u16;
            header_top.push(cell_from_styled_text(slot_header));
            header_bottom.push(plain_header_cell("integral"));
            widths.push(Constraint::Length(integral_width));
            if self.metric_visibility.relative_error {
                header_top.push(blank_cell());
                header_bottom.push(plain_header_cell("% err"));
                widths.push(Constraint::Length(8));
            }
            if self.metric_visibility.chi_sq {
                header_top.push(blank_cell());
                header_bottom.push(plain_header_cell("chi^2"));
                widths.push(Constraint::Length(8));
            }
            if self.metric_visibility.max_weight_impact {
                header_top.push(blank_cell());
                header_bottom.push(plain_header_cell("m.w.i"));
                widths.push(Constraint::Length(10));
            }
        }

        let mut table_rows = vec![
            blank_table_row(widths.len()),
            Row::new(header_top),
            Row::new(header_bottom),
            blank_table_row(widths.len()),
        ];
        for group in [MainResultsRowGroupKind::All, MainResultsRowGroupKind::Sum] {
            let group_rows = update
                .main_results
                .row_groups_of_kind(group)
                .flat_map(|row_group| row_group.rows.iter())
                .collect::<Vec<_>>();
            if group_rows.is_empty() {
                continue;
            }
            if table_rows.len() > 1 {
                table_rows.push(blank_table_row(widths.len()));
            }
            for row in group_rows {
                let mut cells = vec![
                    cell_from_styled_text(&row.contribution.display),
                    cell_from_styled_text(&row.component.display),
                ];
                for slot_index in 0..update.main_results.slot_headers.len() {
                    let cell = row.slot_cell(slot_index);
                    cells.push(cell_from_optional_styled_text(
                        cell.and_then(|cell| cell.value.as_ref())
                            .map(|value| &value.display),
                    ));
                    if self.metric_visibility.relative_error {
                        cells.push(cell_from_optional_styled_text(
                            cell.and_then(|cell| cell.relative_error.as_ref())
                                .map(|value| &value.display),
                        ));
                    }
                    if self.metric_visibility.chi_sq {
                        cells.push(cell_from_optional_styled_text(
                            cell.and_then(|cell| cell.chi_sq.as_ref())
                                .map(|value| &value.display),
                        ));
                    }
                    if self.metric_visibility.max_weight_impact {
                        cells.push(cell_from_optional_styled_text(
                            cell.and_then(|cell| cell.max_weight_impact.as_ref())
                                .map(|value| &value.display),
                        ));
                    }
                }
                table_rows.push(Row::new(cells));
            }
        }

        let table = Table::new(table_rows, widths)
            .block(titled_block("Results summary"))
            .column_spacing(1);
        frame.render_widget(table, area);
    }

    fn draw_statistics_panel(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let Some(statistics) = update.statistics.as_ref() else {
            frame.render_widget(
                Paragraph::new("Statistics unavailable.")
                    .block(titled_block("Integration statistics")),
                area,
            );
            return;
        };
        let statistics_scope = self.selected_statistics_scope(update);
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Percentage(38),
                Constraint::Percentage(31),
                Constraint::Percentage(31),
            ])
            .split(area);

        let mut summary_rows = vec![blank_table_row(9)];
        summary_rows.extend(
            statistics
                .table_rows(statistics_scope)
                .into_iter()
                .map(|row| {
                    let mut cells = vec![cell_from_styled_text(&row.row_label)];
                    for entry in row.entries {
                        cells.push(cell_from_styled_text(&entry.label));
                        cells.push(cell_from_styled_text(&entry.value));
                    }
                    Row::new(cells)
                })
                .collect::<Vec<_>>(),
        );
        summary_rows.push(blank_table_row(9));
        let title = statistics.statistics_title(statistics_scope);
        let summary_table = Table::new(
            summary_rows,
            [
                Constraint::Length(10),
                Constraint::Length(12),
                Constraint::Length(12),
                Constraint::Length(12),
                Constraint::Length(12),
                Constraint::Length(12),
                Constraint::Length(12),
                Constraint::Length(13),
                Constraint::Min(12),
            ],
        )
        .block(titled_block_styled(&title))
        .column_spacing(1);
        frame.render_widget(summary_table, layout[0]);

        self.draw_mix_panel(
            frame,
            layout[1],
            &statistics.timing_title(statistics_scope),
            &statistics.timing_mix_segments(statistics_scope),
        );
        self.draw_mix_panel(
            frame,
            layout[2],
            &statistics.precision_title(statistics_scope),
            &statistics.precision_mix_segments(statistics_scope),
        );
    }

    fn draw_discrete_tab(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let discrete_rows = self.discrete_rows();
        if discrete_rows.is_empty() {
            frame.render_widget(
                Paragraph::new("No monitored discrete breakdown is available for this run.")
                    .block(titled_block("Discrete breakdown"))
                    .wrap(Wrap { trim: false }),
                area,
            );
            return;
        }

        let focused_slot = self
            .focused_slot
            .min(update.main_results.slot_headers.len().saturating_sub(1));
        let horizontal = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Percentage(55), Constraint::Percentage(45)])
            .split(area);

        let contribution_header = update.main_results.contribution_header.to_plain_string();
        let contribution_width = discrete_rows
            .iter()
            .map(|entry| {
                entry
                    .row
                    .contribution
                    .display
                    .to_plain_string()
                    .chars()
                    .count()
            })
            .chain(std::iter::once(contribution_header.chars().count()))
            .max()
            .unwrap_or(24)
            .max(24) as u16;
        let integral_width = discrete_rows
            .iter()
            .filter_map(|entry| {
                entry
                    .row
                    .slot_cell(focused_slot)
                    .and_then(|cell| cell.value.as_ref())
                    .map(|value| value.display.to_plain_string().chars().count())
            })
            .chain(std::iter::once("integral".chars().count()))
            .max()
            .unwrap_or(14)
            .max(14) as u16;
        let mut rows = vec![blank_table_row(9)];
        rows.extend(
            discrete_rows
                .iter()
                .map(|entry| {
                    let cell = entry.row.slot_cell(focused_slot);
                    Row::new(vec![
                        cell_from_styled_text(&entry.row.contribution.display),
                        cell_from_styled_text(&entry.row.component.display),
                        cell_from_optional_styled_text(
                            cell.and_then(|cell| cell.value.as_ref())
                                .map(|value| &value.display),
                        ),
                        cell_from_optional_styled_text(
                            cell.and_then(|cell| cell.relative_error.as_ref())
                                .map(|value| &value.display),
                        ),
                        cell_from_optional_styled_text(
                            cell.and_then(|cell| cell.chi_sq.as_ref())
                                .map(|value| &value.display),
                        ),
                        cell_from_optional_styled_text(
                            cell.and_then(|cell| cell.max_weight_impact.as_ref())
                                .map(|value| &value.display),
                        ),
                        cell_from_optional_styled_text(
                            cell.and_then(|cell| cell.sample_fraction.as_ref())
                                .map(|value| &value.display),
                        ),
                        cell_from_optional_styled_text(
                            cell.and_then(|cell| cell.sample_count.as_ref())
                                .map(|value| &value.display),
                        ),
                        cell_from_optional_styled_text(
                            cell.and_then(|cell| cell.target_pdf.as_ref())
                                .map(|value| &value.display),
                        ),
                    ])
                })
                .collect::<Vec<_>>(),
        );

        let mut state = TableState::default().with_selected(Some(self.selected_discrete_row + 1));
        let table = Table::new(
            rows,
            [
                Constraint::Length(contribution_width),
                Constraint::Length(3),
                Constraint::Length(integral_width),
                Constraint::Length(8),
                Constraint::Length(8),
                Constraint::Length(10),
                Constraint::Length(8),
                Constraint::Length(10),
                Constraint::Min(8),
            ],
        )
        .header(Row::new(
            update
                .main_results
                .discrete_headers()
                .iter()
                .map(cell_from_styled_text)
                .collect::<Vec<_>>(),
        ))
        .row_highlight_style(Style::default().add_modifier(Modifier::REVERSED))
        .block(titled_block_styled(&scoped_integrand_title(
            "Discrete bins for focused integrand ",
            &update.main_results.slot_headers[focused_slot],
            &format!(
                " (sort: {} {})",
                discrete_sort_label(self.discrete_sort),
                if self.discrete_descending {
                    "desc"
                } else {
                    "asc"
                }
            ),
        )));
        frame.render_stateful_widget(table, horizontal[0], &mut state);

        let selected = discrete_rows
            .get(self.selected_discrete_row)
            .copied()
            .unwrap_or(discrete_rows[0]);
        self.draw_selected_discrete_detail(frame, horizontal[1], update, selected.row);
    }

    fn draw_selected_discrete_detail(
        &self,
        frame: &mut Frame<'_>,
        area: Rect,
        update: &StatusUpdate,
        row: &MainResultsRow,
    ) {
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([Constraint::Length(3), Constraint::Min(6)])
            .split(area);
        let mut selected_line = vec![Span::raw("selected: ")];
        selected_line.extend(spans_from_styled_text(&row.contribution.display));
        selected_line.push(Span::raw(" "));
        selected_line.extend(spans_from_styled_text(&row.component.display));
        frame.render_widget(
            Paragraph::new(Line::from(selected_line)).block(titled_block("Selected bin")),
            layout[0],
        );

        let slot_width = update
            .main_results
            .slot_headers
            .iter()
            .map(|slot| slot.to_plain_string().chars().count())
            .max()
            .unwrap_or(8)
            .max(8) as u16;
        let integral_width = (0..update.main_results.slot_headers.len())
            .filter_map(|slot_index| {
                row.slot_cell(slot_index)
                    .and_then(|cell| cell.value.as_ref())
                    .map(|value| value.display.to_plain_string().chars().count())
            })
            .chain(std::iter::once("integral".chars().count()))
            .max()
            .unwrap_or(14)
            .max(14) as u16;
        let mut rows = vec![blank_table_row(8)];
        rows.extend(
            (0..update.main_results.slot_headers.len())
                .filter_map(|slot_index| {
                    let cell = row.slot_cell(slot_index)?;
                    Some(Row::new(vec![
                        cell_from_styled_text(&update.main_results.slot_headers[slot_index]),
                        cell_from_optional_styled_text(
                            cell.value.as_ref().map(|value| &value.display),
                        ),
                        cell_from_optional_styled_text(
                            cell.relative_error.as_ref().map(|value| &value.display),
                        ),
                        cell_from_optional_styled_text(
                            cell.chi_sq.as_ref().map(|value| &value.display),
                        ),
                        cell_from_optional_styled_text(
                            cell.max_weight_impact.as_ref().map(|value| &value.display),
                        ),
                        cell_from_optional_styled_text(
                            cell.sample_fraction.as_ref().map(|value| &value.display),
                        ),
                        cell_from_optional_styled_text(
                            cell.sample_count.as_ref().map(|value| &value.display),
                        ),
                        cell_from_optional_styled_text(
                            cell.target_pdf.as_ref().map(|value| &value.display),
                        ),
                    ]))
                })
                .collect::<Vec<_>>(),
        );

        let table = Table::new(
            rows,
            [
                Constraint::Length(slot_width),
                Constraint::Length(integral_width),
                Constraint::Length(10),
                Constraint::Length(8),
                Constraint::Length(10),
                Constraint::Length(10),
                Constraint::Length(10),
                Constraint::Min(8),
            ],
        )
        .header(Row::new(
            update
                .main_results
                .selected_bin_detail_headers()
                .iter()
                .map(cell_from_styled_text)
                .collect::<Vec<_>>(),
        ))
        .block(titled_block("Per integrand details"))
        .column_spacing(1);
        frame.render_widget(table, layout[1]);
    }

    fn draw_max_weight_tab(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let focused_slot = self
            .focused_slot
            .min(update.main_results.slot_headers.len().saturating_sub(1));
        let focused_slot_name = update.main_results.slot_headers[focused_slot].to_plain_string();
        let top_height = if area.height > 18 {
            update
                .max_weight_details
                .as_ref()
                .map(|section| {
                    let row_count = section
                        .rows_by_slot
                        .iter()
                        .map(|rows| rows.len())
                        .sum::<usize>() as u16;
                    row_count
                        .saturating_add(6)
                        .clamp(8, area.height.saturating_sub(10))
                })
                .unwrap_or(area.height.saturating_mul(2) / 5)
        } else {
            area.height.saturating_mul(2) / 5
        };
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([Constraint::Length(top_height), Constraint::Min(8)])
            .split(area);

        if let Some(section) = update.max_weight_details.as_ref() {
            let headers = section.headers();
            let title = section.title();
            let max_eval_width = max_styled_text_width(
                std::iter::once(&headers[2]).chain(
                    section
                        .rows_by_slot
                        .iter()
                        .flat_map(|group| group.iter().map(|row| &row.max_eval.display)),
                ),
            );
            let mut rows = vec![blank_table_row(4)];
            rows.extend(
                section
                    .rows_by_slot
                    .iter()
                    .flat_map(|group| group.iter())
                    .map(|row| {
                        Row::new(vec![
                            cell_from_styled_text(&row.slot.display),
                            cell_from_styled_text(&row.component_sign.display),
                            cell_from_styled_text(&row.max_eval.display),
                            cell_from_styled_text(&row.coordinates.display),
                        ])
                        .height(styled_text_line_count(&row.coordinates.display))
                    })
                    .collect::<Vec<_>>(),
            );
            rows.push(blank_table_row(4));
            let table = Table::new(
                rows,
                [
                    Constraint::Length(24),
                    Constraint::Length(10),
                    Constraint::Length(max_eval_width),
                    Constraint::Min(24),
                ],
            )
            .header(Row::new(
                headers
                    .iter()
                    .map(cell_from_styled_text)
                    .collect::<Vec<_>>(),
            ))
            .block(titled_block_styled(&title));
            frame.render_widget(table, layout[0]);
        } else {
            frame.render_widget(
                Paragraph::new("No overall max-weight data.")
                    .block(titled_block("Overall max weight")),
                layout[0],
            );
        }

        if let Some(section) = update.discrete_max_weight_details.as_ref() {
            let sorted_rows = self.sorted_discrete_max_weight_rows(section);
            let lower = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([Constraint::Percentage(40), Constraint::Percentage(60)])
                .split(layout[1]);
            let summary_headers = section.summary_headers();
            let mut summary_rows = vec![blank_table_row(section.slot_headers.len() + 2)];
            summary_rows.extend(
                sorted_rows
                    .iter()
                    .map(|row| {
                        let mut cells = vec![
                            cell_from_styled_text(&row.contribution.display),
                            cell_from_styled_text(&row.component_sign.display),
                        ];
                        cells.extend(
                            row.slot_values
                                .iter()
                                .map(|value| {
                                    cell_from_optional_styled_text(
                                        value.as_ref().map(|field| &field.display),
                                    )
                                })
                                .collect::<Vec<_>>(),
                        );
                        Row::new(cells)
                    })
                    .collect::<Vec<_>>(),
            );
            let summary_widths = std::iter::once(Constraint::Length(24))
                .chain(std::iter::once(Constraint::Length(10)))
                .chain(
                    section
                        .slot_headers
                        .iter()
                        .enumerate()
                        .map(|(slot_index, header)| {
                            Constraint::Min(
                                max_styled_text_width(std::iter::once(header).chain(
                                    sorted_rows.iter().filter_map(|row| {
                                        row.slot_values
                                            .get(slot_index)
                                            .and_then(|value| value.as_ref())
                                            .map(|field| &field.display)
                                    }),
                                ))
                                .max(16),
                            )
                        }),
                )
                .collect::<Vec<_>>();
            let table = Table::new(summary_rows, summary_widths)
                .header(Row::new(
                    summary_headers
                        .iter()
                        .map(cell_from_styled_text)
                        .collect::<Vec<_>>(),
                ))
                .block(titled_block(format!(
                    "Per-bin max weight (sort: {} {})",
                    discrete_sort_label(self.discrete_sort),
                    if self.discrete_descending {
                        "desc"
                    } else {
                        "asc"
                    }
                )))
                .column_spacing(1);
            frame.render_widget(table, lower[0]);

            let coordinate_headers = section.coordinate_headers();
            let focused_slot_filter = focused_slot_name.clone();
            let mut coordinate_rows = vec![blank_table_row(3)];
            coordinate_rows.extend(
                sorted_rows
                    .iter()
                    .flat_map(|row| {
                        let focused_slot_filter = focused_slot_filter.clone();
                        row.slot_coordinates
                            .iter()
                            .filter(move |entry| entry.slot.raw == focused_slot_filter)
                            .map(move |entry| {
                                Row::new(vec![
                                    cell_from_styled_text(&row.contribution.display),
                                    cell_from_styled_text(&row.component_sign.display),
                                    cell_from_styled_text(&entry.coordinates.display),
                                ])
                                .height(styled_text_line_count(&entry.coordinates.display))
                            })
                    })
                    .collect::<Vec<_>>(),
            );
            let coordinates_table = Table::new(
                coordinate_rows,
                [
                    Constraint::Length(18),
                    Constraint::Length(8),
                    Constraint::Min(24),
                ],
            )
            .header(Row::new(vec![
                cell_from_styled_text(&coordinate_headers[0]),
                cell_from_styled_text(&coordinate_headers[1]),
                cell_from_styled_text(&coordinate_headers[3]),
            ]))
            .block(titled_block_styled(&scoped_integrand_title(
                "Maximum weight coordinates for integrand ",
                &update.main_results.slot_headers[focused_slot],
                "",
            )));
            frame.render_widget(coordinates_table, lower[1]);
        } else {
            frame.render_widget(
                Paragraph::new("No per-bin max-weight data.")
                    .block(titled_block("Per-bin max weight")),
                layout[1],
            );
        }
    }

    fn draw_footer(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let phase = update
            .training_component()
            .map(|component| match component {
                ComponentKind::Real => "re",
                ComponentKind::Imag => "im",
            })
            .unwrap_or("n/a");
        let footer = Paragraph::new(Line::from(vec![
            Span::styled("Tabs", label_style()),
            Span::raw(" 1/2/3 or <- ->  "),
            Span::styled("Integrand", label_style()),
            Span::raw(" [ / ]  "),
            Span::styled("Bins", label_style()),
            Span::raw(" j/k  "),
            Span::styled("Metrics", label_style()),
            Span::raw(" r c w  "),
            Span::styled("Phase", label_style()),
            Span::raw(" p  "),
            Span::styled("Sort", label_style()),
            Span::raw(" s/v  "),
            Span::styled("History", label_style()),
            Span::raw(format!(
                " g +/- ({})  ",
                self.chart_history_window.description()
            )),
            Span::styled("Y-axis", label_style()),
            Span::raw(format!(" , . 0 (±{}σ)  ", self.chart_y_sigma_span)),
            Span::styled("Help", label_style()),
            Span::raw(" ?  "),
            Span::styled("Abort", label_style()),
            Span::raw(format!(" x / Ctrl-C   training {phase}")),
        ]))
        .alignment(Alignment::Left)
        .block(Block::default().borders(Borders::TOP));
        frame.render_widget(footer, area);
    }

    fn draw_help_overlay(&self, frame: &mut Frame<'_>, area: Rect) {
        let popup = centered_rect(72, 65, area);
        frame.render_widget(Clear, popup);
        frame.render_widget(
            Paragraph::new(Text::from(vec![
                Line::from("1/2/3 or Left/Right  switch tabs"),
                Line::from("[ and ]             change focused integrand"),
                Line::from("j / k               move selected discrete row"),
                Line::from("s                   cycle discrete sort"),
                Line::from("v                   reverse discrete sort order"),
                Line::from("i                   toggle global/focused statistics"),
                Line::from("r / c / w           toggle rel err, chi^2, mwi columns"),
                Line::from("p                   toggle convergence real/imag phase"),
                Line::from("g                   toggle full/recent chart history"),
                Line::from("+ / -               widen or narrow recent-history window"),
                Line::from(", / .               tighten or widen the y-range in σ units"),
                Line::from("0                   reset the y-range to ±4σ"),
                Line::from("x                   abort the current iteration"),
                Line::from("Ctrl-C              stop the integration"),
                Line::from("?                   close help"),
            ]))
            .block(titled_block("Dashboard controls"))
            .alignment(Alignment::Left)
            .wrap(Wrap { trim: false }),
            popup,
        );
    }

    fn draw_mix_panel(
        &self,
        frame: &mut Frame<'_>,
        area: Rect,
        title: &StyledText,
        segments: &[StatisticsMixSegment],
    ) {
        let block = titled_block_styled(title);
        let inner = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(1),
                Constraint::Length(1),
                Constraint::Length(1),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(block.inner(area));
        frame.render_widget(block, area);
        frame.render_widget(Paragraph::new(""), inner[0]);
        frame.render_widget(
            Paragraph::new(mix_bar_line(segments, inner[1].width as usize)),
            inner[1],
        );
        frame.render_widget(Paragraph::new(""), inner[2]);
        if !segments.is_empty() {
            let constraints = (0..segments.len())
                .map(|_| Constraint::Percentage((100 / segments.len().max(1)) as u16))
                .collect::<Vec<_>>();
            let labels = Layout::default()
                .direction(Direction::Horizontal)
                .constraints(constraints)
                .split(inner[3]);
            for (segment, label_area) in segments.iter().zip(labels.iter().copied()) {
                frame.render_widget(
                    Paragraph::new(Line::from(vec![Span::styled(
                        format!(
                            "{} : {:.1}%",
                            segment.label.to_plain_string(),
                            segment.percentage
                        ),
                        style_from_text_style(
                            segment
                                .label
                                .spans
                                .first()
                                .map(|span| span.style)
                                .unwrap_or(TextStyle::PLAIN),
                        )
                        .add_modifier(Modifier::BOLD),
                    )]))
                    .alignment(Alignment::Center),
                    label_area,
                );
            }
        }
        frame.render_widget(Paragraph::new(""), inner[4]);
    }

    fn discrete_rows(&self) -> Vec<DiscreteRowRef<'_>> {
        let Some(update) = self.latest_update.as_ref() else {
            return Vec::new();
        };
        let focused_slot = self
            .focused_slot
            .min(update.main_results.slot_headers.len().saturating_sub(1));
        self.discrete_rows_for_slot(focused_slot)
    }

    fn selected_statistics_scope(&self, update: &StatusUpdate) -> StatisticsScope {
        match self.statistics_scope {
            DashboardStatisticsScope::Global => StatisticsScope::Global,
            DashboardStatisticsScope::FocusedSlot => StatisticsScope::Slot(
                self.focused_slot
                    .min(update.main_results.slot_headers.len().saturating_sub(1)),
            ),
        }
    }

    fn discrete_rows_for_slot(&self, slot_index: usize) -> Vec<DiscreteRowRef<'_>> {
        let Some(update) = self.latest_update.as_ref() else {
            return Vec::new();
        };
        let mut rows = update
            .main_results
            .row_groups_of_kind(MainResultsRowGroupKind::Bins)
            .flat_map(|group| group.rows.iter())
            .map(|row| DiscreteRowRef { row })
            .collect::<Vec<_>>();

        rows.sort_by(|lhs, rhs| {
            let lhs_key = lhs.sort_value(self.discrete_sort, slot_index);
            let rhs_key = rhs.sort_value(self.discrete_sort, slot_index);
            let ordering = match self.discrete_sort {
                ContributionSortMode::Index => match (lhs.contribution(), rhs.contribution()) {
                    (ContributionKind::Bin(lhs_index), ContributionKind::Bin(rhs_index)) => {
                        lhs_index
                            .cmp(&rhs_index)
                            .then(lhs.component().cmp(&rhs.component()))
                    }
                    _ => Ordering::Equal,
                },
                ContributionSortMode::Integral | ContributionSortMode::Error => lhs_key
                    .partial_cmp(&rhs_key)
                    .unwrap_or(Ordering::Equal)
                    .then(lhs.component().cmp(&rhs.component())),
            };
            if self.discrete_descending {
                ordering.reverse()
            } else {
                ordering
            }
        });
        rows
    }

    fn sorted_discrete_max_weight_rows<'a>(
        &self,
        section: &'a super::status_update::DiscreteMaxWeightDetailsSection,
    ) -> Vec<&'a super::status_update::DiscreteMaxWeightRow> {
        let discrete_order = self
            .discrete_rows_for_slot(0)
            .into_iter()
            .map(|row| (row.contribution(), row.component()))
            .collect::<Vec<_>>();
        let mut all_rows = section
            .row_groups
            .iter()
            .flat_map(|group| group.iter())
            .collect::<Vec<_>>();
        all_rows.sort_by_key(|row| {
            self.max_weight_row_sort_rank(
                &discrete_order,
                row.contribution.raw,
                row.component_sign.raw.0,
            )
        });
        all_rows
    }

    fn max_weight_row_sort_rank(
        &self,
        discrete_order: &[(ContributionKind, ComponentKind)],
        contribution: ContributionKind,
        component: ComponentKind,
    ) -> (usize, usize) {
        match contribution {
            ContributionKind::All => (0, component_rank(component)),
            ContributionKind::Bin(_) => (
                1,
                discrete_order
                    .iter()
                    .position(|entry| *entry == (contribution, component))
                    .unwrap_or(usize::MAX),
            ),
            ContributionKind::Sum => (2, component_rank(component)),
        }
    }

    fn overview_value_text(&self, update: &StatusUpdate, slot_index: usize) -> Text<'static> {
        overview_metric_text(update, slot_index, |cell| {
            cell.value.as_ref().map(|value| &value.display)
        })
    }

    fn overview_summary_rows<'a>(&self, update: &'a StatusUpdate) -> Vec<&'a MainResultsRow> {
        update
            .main_results
            .row_groups_of_kind(MainResultsRowGroupKind::All)
            .chain(
                update
                    .main_results
                    .row_groups_of_kind(MainResultsRowGroupKind::Sum),
            )
            .flat_map(|group| group.rows.iter())
            .collect()
    }

    fn results_summary_cell(&self, row: &MainResultsRow, slot_index: usize) -> Text<'static> {
        let Some(cell) = row.slot_cell(slot_index) else {
            return Text::from("N/A");
        };

        let mut lines = Vec::new();
        if let Some(value) = cell.value.as_ref() {
            lines.push(Line::from(spans_from_styled_text(&value.display)));
        } else {
            lines.push(Line::from("N/A"));
        }

        let mut metrics = Vec::new();
        if self.metric_visibility.relative_error {
            metrics.push(labeled_value_line(
                "rel err",
                cell.relative_error.as_ref().map(|value| &value.display),
            ));
        }
        if self.metric_visibility.chi_sq {
            metrics.push(labeled_value_line(
                "chi^2",
                cell.chi_sq.as_ref().map(|value| &value.display),
            ));
        }
        if self.metric_visibility.max_weight_impact {
            metrics.push(labeled_value_line(
                "mwi",
                cell.max_weight_impact.as_ref().map(|value| &value.display),
            ));
        }
        if !metrics.is_empty() && !matches!(self.density, DensityMode::Compact) {
            let mut spans = Vec::new();
            for (index, metric_line) in metrics.into_iter().enumerate() {
                if index > 0 {
                    spans.push(Span::raw("   "));
                }
                spans.extend(metric_line.spans);
            }
            lines.push(Line::from(spans));
        }

        Text::from(lines)
    }
}

fn overview_metric_text(
    update: &StatusUpdate,
    slot_index: usize,
    select: impl Fn(&MainTableSlotCells) -> Option<&StyledText>,
) -> Text<'static> {
    let select = &select;
    let lines = update
        .main_results
        .row_groups_of_kind(MainResultsRowGroupKind::All)
        .flat_map(|group| group.rows.iter())
        .filter_map(|row| {
            let text = row.slot_cell(slot_index).and_then(select)?;
            let mut spans = spans_from_styled_text(&row.component.display);
            spans.push(Span::raw(": "));
            spans.extend(spans_from_styled_text(text));
            Some(Line::from(spans))
        })
        .collect::<Vec<_>>();

    if lines.is_empty() {
        Text::from("N/A")
    } else {
        Text::from(lines)
    }
}

fn component_label_with_colon(component: ComponentKind) -> StyledText {
    let mut label = component.label_display();
    label.push_text(":", component.text_style());
    label
}

fn component_rank(component: ComponentKind) -> usize {
    match component {
        ComponentKind::Real => 0,
        ComponentKind::Imag => 1,
    }
}

fn scoped_integrand_title(prefix: &str, slot_header: &StyledText, suffix: &str) -> StyledText {
    let mut title = StyledText::styled(prefix, TextStyle::green().bold());
    title.push_text("[", TextStyle::green().bold());
    title.append(slot_header.clone());
    title.push_text("]", TextStyle::green().bold());
    title.push_text(suffix, TextStyle::green().bold());
    title
}

fn format_samples_per_second(rate: f64) -> String {
    abbreviate_count(rate.max(0.0).round() as usize)
}

fn labeled_value_line(label: &str, value: Option<&StyledText>) -> Line<'static> {
    let mut spans = vec![Span::raw(format!("{label} "))];
    if let Some(value) = value {
        spans.extend(spans_from_styled_text(value));
    } else {
        spans.push(Span::raw("N/A"));
    }
    Line::from(spans)
}

fn style_from_text_style(style: TextStyle) -> Style {
    let mut rendered = Style::default();
    if let Some(color) = color_from_text_style(style) {
        rendered = rendered.fg(color);
    }
    if style.bold {
        rendered = rendered.add_modifier(Modifier::BOLD);
    }
    if style.dimmed {
        rendered = rendered.add_modifier(Modifier::DIM);
    }
    rendered
}

fn spans_from_styled_text(text: &StyledText) -> Vec<Span<'static>> {
    text.spans
        .iter()
        .map(|span| Span::styled(span.text.clone(), style_from_text_style(span.style)))
        .collect()
}

fn line_from_styled_text(text: &StyledText) -> Line<'static> {
    Line::from(spans_from_styled_text(text))
}

fn text_from_styled_text(text: &StyledText) -> Text<'static> {
    let mut lines = vec![Line::default()];
    for span in &text.spans {
        let style = style_from_text_style(span.style);
        for (index, segment) in span.text.split('\n').enumerate() {
            if index > 0 {
                lines.push(Line::default());
            }
            if !segment.is_empty() {
                lines
                    .last_mut()
                    .expect("at least one line")
                    .spans
                    .push(Span::styled(segment.to_string(), style));
            }
        }
    }
    Text::from(lines)
}

fn styled_text_line_count(text: &StyledText) -> u16 {
    text.to_plain_string()
        .lines()
        .count()
        .max(1)
        .try_into()
        .unwrap_or(u16::MAX)
}

fn styled_text_max_line_width(text: &StyledText) -> u16 {
    text.to_plain_string()
        .lines()
        .map(|line| line.chars().count())
        .max()
        .unwrap_or(0)
        .max(1)
        .try_into()
        .unwrap_or(u16::MAX)
}

fn max_styled_text_width<'a>(texts: impl IntoIterator<Item = &'a StyledText>) -> u16 {
    texts
        .into_iter()
        .map(styled_text_max_line_width)
        .max()
        .unwrap_or(1)
}

fn styled_text_line(text: &StyledText, line_index: usize) -> StyledText {
    let mut lines = vec![StyledText::new()];
    for span in &text.spans {
        for (index, segment) in span.text.split('\n').enumerate() {
            if index > 0 {
                lines.push(StyledText::new());
            }
            if !segment.is_empty() {
                lines
                    .last_mut()
                    .expect("at least one styled line")
                    .push_text(segment, span.style);
            }
        }
    }
    lines.get(line_index).cloned().unwrap_or_default()
}

fn cell_from_styled_text(text: &StyledText) -> Cell<'static> {
    Cell::from(text_from_styled_text(text))
}

fn cell_from_optional_styled_text(text: Option<&StyledText>) -> Cell<'static> {
    text.map(cell_from_styled_text).unwrap_or_else(blank_cell)
}

fn plain_header_cell(text: impl Into<String>) -> Cell<'static> {
    Cell::from(Span::styled(text.into(), label_style()))
}

fn blank_cell() -> Cell<'static> {
    Cell::from(String::new())
}

fn blank_table_row(column_count: usize) -> Row<'static> {
    Row::new((0..column_count).map(|_| blank_cell()).collect::<Vec<_>>())
}

fn stacked_header(title: &StyledText, subtitle: &str) -> StyledText {
    let mut text = title.clone();
    text.push_text("\n", TextStyle::PLAIN);
    text.push_text(subtitle, TextStyle::PLAIN.bold());
    text
}

fn titled_block(title: impl Into<String>) -> Block<'static> {
    Block::default()
        .borders(Borders::ALL)
        .title(Line::from(vec![Span::styled(
            title.into(),
            Style::default()
                .fg(Color::Green)
                .add_modifier(Modifier::BOLD),
        )]))
        .border_style(Style::default().fg(Color::DarkGray))
}

fn titled_block_styled(title: &StyledText) -> Block<'static> {
    Block::default()
        .borders(Borders::ALL)
        .title(line_from_styled_text(title))
        .border_style(Style::default().fg(Color::DarkGray))
}

fn label_style() -> Style {
    Style::default().add_modifier(Modifier::BOLD)
}

fn abbreviate_count(value: usize) -> String {
    super::display::format_abbreviated_count(value)
}

fn mix_bar_line(segments: &[StatisticsMixSegment], width: usize) -> Line<'static> {
    let bar_width = width.max(8);
    let mut spans = Vec::new();
    let mut used = 0usize;
    for (index, segment) in segments.iter().enumerate() {
        let color = color_from_styled_text(&segment.label, Color::DarkGray);
        let width = if index + 1 == segments.len() {
            bar_width.saturating_sub(used)
        } else {
            let segment_width = ((segment.percentage / 100.0) * bar_width as f64).round() as usize;
            used += segment_width;
            segment_width
        };
        spans.push(Span::styled("█".repeat(width), Style::default().fg(color)));
    }
    Line::from(spans)
}

fn discrete_sort_label(mode: ContributionSortMode) -> &'static str {
    match mode {
        ContributionSortMode::Index => "index",
        ContributionSortMode::Integral => "integral",
        ContributionSortMode::Error => "error",
    }
}

fn color_from_text_style(style: TextStyle) -> Option<Color> {
    match style.color {
        Some(TextColor::Green) => Some(Color::Green),
        Some(TextColor::Blue) => Some(Color::Blue),
        Some(TextColor::Red) => Some(Color::Red),
        Some(TextColor::Yellow) => Some(Color::Yellow),
        Some(TextColor::Pink) => Some(Color::LightMagenta),
        None => None,
    }
}

fn color_from_styled_text(text: &StyledText, fallback: Color) -> Color {
    text.spans
        .iter()
        .find_map(|span| color_from_text_style(span.style))
        .unwrap_or(fallback)
}

fn count_axis_labels(min: f64, max: f64) -> Vec<Span<'static>> {
    evenly_spaced_values(min, max, 9)
        .into_iter()
        .map(|value| Span::raw(abbreviate_count(value.max(0.0).round() as usize)))
        .collect()
}

fn value_axis_labels(min: f64, max: f64) -> Vec<Span<'static>> {
    evenly_spaced_values(min, max, 9)
        .into_iter()
        .map(|value| Span::raw(format!("{value:.3e}")))
        .collect()
}

fn evenly_spaced_values(min: f64, max: f64, count: usize) -> Vec<f64> {
    if count <= 1 || !min.is_finite() || !max.is_finite() || (max - min).abs() < f64::EPSILON {
        return vec![min];
    }

    let step = (max - min) / (count.saturating_sub(1) as f64);
    (0..count).map(|index| min + step * index as f64).collect()
}

fn centered_rect(percent_x: u16, percent_y: u16, area: Rect) -> Rect {
    let popup_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Percentage((100 - percent_y) / 2),
            Constraint::Percentage(percent_y),
            Constraint::Percentage((100 - percent_y) / 2),
        ])
        .split(area);
    Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage((100 - percent_x) / 2),
            Constraint::Percentage(percent_x),
            Constraint::Percentage((100 - percent_x) / 2),
        ])
        .split(popup_layout[1])[1]
}
