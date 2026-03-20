use std::cmp::Ordering;

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

use crate::{settings::runtime::IntegrationStatisticsSnapshot, utils};

use super::{
    StatusUpdate,
    display::{StyledText, TextColor, TextStyle},
    status_update::{
        ComponentKind, ContributionKind, ContributionSortMode, MainResultsRow,
        MainResultsRowGroupKind, MainTableSlotCells,
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

#[derive(Clone, Copy, Debug)]
struct HistoryPoint {
    samples: usize,
    central: f64,
    error: f64,
    completed_iteration: bool,
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

    fn sort_value(self, mode: ContributionSortMode) -> f64 {
        let Some(cell) = self.row.slot_cell(0) else {
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
    selected_discrete_row: usize,
    discrete_sort: ContributionSortMode,
    discrete_descending: bool,
    show_help: bool,
    history: Vec<HistoryPoint>,
}

impl Default for RatatuiDashboardState {
    fn default() -> Self {
        Self {
            latest_update: None,
            active_tab: DashboardTab::Overview,
            density: DensityMode::Metrics,
            metric_visibility: SlotMetricVisibility::default(),
            focused_slot: 0,
            selected_discrete_row: 0,
            discrete_sort: ContributionSortMode::Error,
            discrete_descending: true,
            show_help: false,
            history: Vec::new(),
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

        let row_count = self.discrete_rows().len();
        if row_count == 0 {
            self.selected_discrete_row = 0;
        } else {
            self.selected_discrete_row = self.selected_discrete_row.min(row_count - 1);
        }
    }

    fn push_history_point(&mut self, update: &StatusUpdate) {
        let Some(component) = update.training_component() else {
            return;
        };
        let Some(row) = update
            .main_results
            .find_row(ContributionKind::All, component)
        else {
            return;
        };
        let Some(value) = row
            .slot_cell(update.meta.training_slot)
            .and_then(|cell| cell.value.as_ref())
        else {
            return;
        };

        let point = HistoryPoint {
            samples: update.meta.total_points,
            central: value.raw.0.0,
            error: value.raw.1.0.abs(),
            completed_iteration: !matches!(update.kind(), super::IntegrationStatusKind::Live),
        };

        if self
            .history
            .last()
            .is_some_and(|previous| previous.samples == point.samples)
        {
            let _ = self.history.pop();
        }
        if self
            .history
            .last()
            .is_some_and(|previous| previous.samples > point.samples)
        {
            self.history.clear();
        }
        self.history.push(point);
        if self.history.len() > 4096 {
            let drain = self.history.len() - 4096;
            self.history.drain(0..drain);
        }
    }

    fn draw_dashboard(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let vertical = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),
                Constraint::Length(4),
                Constraint::Length(if update.main_results.slot_headers.len() > 3 {
                    4
                } else {
                    3
                }),
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
            .highlight_style(
                Style::default()
                    .fg(Color::Black)
                    .bg(Color::Green)
                    .add_modifier(Modifier::BOLD),
            );
        frame.render_widget(tabs, area);
    }

    fn draw_progress(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let inner = Layout::default()
            .direction(Direction::Vertical)
            .constraints([Constraint::Length(1), Constraint::Length(2)])
            .split(area);

        let eta = update
            .meta
            .iteration_eta()
            .map(|duration| utils::format_wdhms(duration.as_secs() as usize))
            .unwrap_or_else(|| "warming up".to_string());
        let throughput = format_rate(update.meta.live_progress.map(|progress| {
            if update.meta.iteration_elapsed_time.is_zero() {
                0.0
            } else {
                progress.completed_points as f64 / update.meta.iteration_elapsed_time.as_secs_f64()
            }
        }));
        let per_sample_core = if update.meta.n_samples_evaluated == 0 {
            "N/A".to_string()
        } else {
            utils::format_evaluation_time_from_f64(
                update.meta.elapsed_time.as_secs_f64() / update.meta.n_samples_evaluated as f64
                    * update.meta.cores as f64,
            )
        };
        let header = Line::from(vec![
            Span::styled(
                format!(" Iteration #{} ", update.meta.iteration),
                Style::default()
                    .fg(Color::Green)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::raw("  "),
            Span::styled(
                format!(
                    "elapsed {}",
                    utils::format_wdhms(update.meta.elapsed_time.as_secs() as usize)
                ),
                Style::default().fg(Color::Blue),
            ),
            Span::raw("  "),
            Span::styled(
                format!("iter ETA {eta}"),
                Style::default().fg(Color::Yellow),
            ),
            Span::raw("  "),
            Span::styled(
                format!("rate {throughput}"),
                Style::default().fg(Color::Blue),
            ),
            Span::raw("  "),
            Span::styled(
                format!("{} /sample/core", per_sample_core),
                Style::default().fg(Color::Blue),
            ),
        ]);
        frame.render_widget(Paragraph::new(header), inner[0]);

        let progress = update.meta.iteration_progress_ratio().unwrap_or_else(|| {
            if matches!(update.kind(), super::IntegrationStatusKind::Live) {
                0.0
            } else {
                1.0
            }
        });
        let gauge = Gauge::default()
            .block(titled_block("Iteration progress"))
            .gauge_style(
                Style::default()
                    .fg(Color::Green)
                    .bg(Color::Black)
                    .add_modifier(Modifier::BOLD),
            )
            .ratio(progress.clamp(0.0, 1.0))
            .label(format!(
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
            ));
        frame.render_widget(gauge, inner[1]);
    }

    fn draw_slot_ribbon(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let mut lines = Vec::new();
        let focused = self.focused_slot;
        for slot_index in 0..update.main_results.slot_headers.len() {
            let title = update.main_results.slot_headers[slot_index].to_plain_string();
            let value = self
                .overview_value_text(update, slot_index)
                .lines
                .first()
                .cloned()
                .unwrap_or_else(|| Line::from("N/A"));
            let mut spans = vec![Span::styled(
                if slot_index == focused {
                    format!("> {title}")
                } else {
                    format!("  {title}")
                },
                if slot_index == focused {
                    Style::default()
                        .fg(Color::Black)
                        .bg(Color::Green)
                        .add_modifier(Modifier::BOLD)
                } else {
                    Style::default()
                        .fg(Color::Blue)
                        .add_modifier(Modifier::BOLD)
                },
            )];
            spans.push(Span::raw("  "));
            spans.extend(value.spans.iter().cloned());
            lines.push(Line::from(spans));
        }

        frame.render_widget(
            Paragraph::new(Text::from(lines))
                .block(titled_block("Slots"))
                .wrap(Wrap { trim: false }),
            area,
        );
    }

    fn draw_overview_tab(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let vertical = if area.width >= 120 {
            Layout::default()
                .direction(Direction::Vertical)
                .constraints([Constraint::Length(12), Constraint::Min(10)])
                .split(area)
        } else {
            Layout::default()
                .direction(Direction::Vertical)
                .constraints([
                    Constraint::Length(10),
                    Constraint::Length(10),
                    Constraint::Min(8),
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
            self.draw_slot_metrics_table(frame, bottom[0], update);
            self.draw_statistics_panel(frame, bottom[1], update);
        } else {
            self.draw_chart(frame, vertical[0], update);
            self.draw_focused_slot_detail(frame, vertical[1], update);
            let bottom = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([Constraint::Percentage(55), Constraint::Percentage(45)])
                .split(vertical[2]);
            self.draw_slot_metrics_table(frame, bottom[0], update);
            self.draw_statistics_panel(frame, bottom[1], update);
        }
    }

    fn draw_chart(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let Some(component) = update.training_component() else {
            frame.render_widget(
                Paragraph::new("Training phase is not a single component.")
                    .block(titled_block("Convergence")),
                area,
            );
            return;
        };

        if self.history.is_empty() {
            frame.render_widget(
                Paragraph::new("Waiting for history points...").block(titled_block("Convergence")),
                area,
            );
            return;
        }

        let central_points = self
            .history
            .iter()
            .map(|point| (point.samples as f64, point.central))
            .collect::<Vec<_>>();
        let upper_points = self
            .history
            .iter()
            .map(|point| (point.samples as f64, point.central + point.error))
            .collect::<Vec<_>>();
        let lower_points = self
            .history
            .iter()
            .map(|point| (point.samples as f64, point.central - point.error))
            .collect::<Vec<_>>();
        let completed_points = self
            .history
            .iter()
            .filter(|point| point.completed_iteration)
            .map(|point| (point.samples as f64, point.central))
            .collect::<Vec<_>>();

        let x_min = self
            .history
            .first()
            .map(|point| point.samples as f64)
            .unwrap_or(0.0);
        let x_max = self
            .history
            .last()
            .map(|point| point.samples as f64)
            .unwrap_or(1.0);
        let mut y_min = lower_points
            .iter()
            .map(|(_, y)| *y)
            .fold(f64::INFINITY, f64::min);
        let mut y_max = upper_points
            .iter()
            .map(|(_, y)| *y)
            .fold(f64::NEG_INFINITY, f64::max);
        if let Some(target) = update.training_target() {
            y_min = y_min.min(target.0);
            y_max = y_max.max(target.0);
        }
        if !y_min.is_finite() || !y_max.is_finite() || (y_max - y_min).abs() < 1.0e-14 {
            y_min -= 1.0;
            y_max += 1.0;
        } else {
            let padding = (y_max - y_min).abs() * 0.08;
            y_min -= padding;
            y_max += padding;
        }

        let mut datasets = vec![
            Dataset::default()
                .name("central")
                .marker(Marker::Braille)
                .graph_type(ratatui::widgets::GraphType::Line)
                .style(
                    Style::default()
                        .fg(Color::Blue)
                        .add_modifier(Modifier::BOLD),
                )
                .data(&central_points),
            Dataset::default()
                .name("upper")
                .graph_type(ratatui::widgets::GraphType::Line)
                .style(Style::default().fg(Color::Cyan))
                .data(&upper_points),
            Dataset::default()
                .name("lower")
                .graph_type(ratatui::widgets::GraphType::Line)
                .style(Style::default().fg(Color::Cyan))
                .data(&lower_points),
            Dataset::default()
                .name("iter")
                .marker(Marker::Dot)
                .graph_type(ratatui::widgets::GraphType::Scatter)
                .style(Style::default().fg(Color::Green))
                .data(&completed_points),
        ];

        let target_points = update
            .training_target()
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

        let title = format!(
            "Convergence ({})",
            match component {
                ComponentKind::Real => "training re",
                ComponentKind::Imag => "training im",
            }
        );
        let chart = Chart::new(datasets)
            .block(titled_block(title))
            .x_axis(
                Axis::default()
                    .title("samples")
                    .bounds([x_min, x_max.max(x_min + 1.0)])
                    .labels(vec![
                        Span::raw(abbreviate_count(x_min.max(0.0) as usize)),
                        Span::raw(abbreviate_count(x_max.max(0.0) as usize)),
                    ]),
            )
            .y_axis(
                Axis::default()
                    .title("central value")
                    .bounds([y_min, y_max])
                    .labels(vec![
                        Span::raw(format!("{y_min:.3e}")),
                        Span::raw(format!("{y_max:.3e}")),
                    ]),
            );
        frame.render_widget(chart, area);
    }

    fn draw_focused_slot_detail(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let focused_slot = self
            .focused_slot
            .min(update.main_results.slot_headers.len().saturating_sub(1));
        let slot_name = update.main_results.slot_headers[focused_slot].to_plain_string();
        let rows = update
            .main_results
            .row_groups_of_kind(MainResultsRowGroupKind::All)
            .flat_map(|group| group.rows.iter())
            .collect::<Vec<_>>();

        let mut lines = vec![Line::from(vec![Span::styled(
            slot_name,
            Style::default()
                .fg(Color::Blue)
                .add_modifier(Modifier::BOLD),
        )])];
        for row in rows {
            let component = match row.component.raw {
                ComponentKind::Real => "re",
                ComponentKind::Imag => "im",
            };
            let Some(cell) = row.slot_cell(focused_slot) else {
                continue;
            };
            let value = cell
                .value
                .as_ref()
                .map(|value| spans_from_styled_text(&value.display))
                .unwrap_or_else(|| vec![Span::raw("N/A")]);
            let rel = cell
                .relative_error
                .as_ref()
                .map(|value| value.display.to_plain_string())
                .unwrap_or_else(|| "N/A".to_string());
            let chi_sq = cell
                .chi_sq
                .as_ref()
                .map(|value| value.display.to_plain_string())
                .unwrap_or_else(|| "N/A".to_string());
            let mwi = cell
                .max_weight_impact
                .as_ref()
                .map(|value| value.display.to_plain_string())
                .unwrap_or_else(|| "N/A".to_string());

            let mut line_spans = vec![Span::styled(
                format!("{component}: "),
                Style::default()
                    .fg(Color::Green)
                    .add_modifier(Modifier::BOLD),
            )];
            line_spans.extend(value);
            lines.push(Line::from(line_spans));
            lines.push(Line::from(format!(
                "   rel err {rel}   chi^2 {chi_sq}   mwi {mwi}"
            )));
        }

        if focused_slot == 0 {
            let deltas = update
                .main_results
                .row_groups_of_kind(MainResultsRowGroupKind::All)
                .flat_map(|group| group.rows.iter())
                .filter_map(|row| {
                    row.delta_sigma
                        .as_ref()
                        .map(|delta| (row.component.raw, delta.display.to_plain_string()))
                })
                .collect::<Vec<_>>();
            if !deltas.is_empty() {
                lines.push(Line::from(""));
                for (component, delta) in deltas {
                    let label = match component {
                        ComponentKind::Real => "target re",
                        ComponentKind::Imag => "target im",
                    };
                    lines.push(Line::from(format!("{label}: {delta}")));
                }
            }
        }

        frame.render_widget(
            Paragraph::new(Text::from(lines))
                .block(titled_block("Focused slot"))
                .wrap(Wrap { trim: false }),
            area,
        );
    }

    fn draw_slot_metrics_table(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let mut header = vec![Cell::from("slot"), Cell::from("value")];
        if self.metric_visibility.relative_error {
            header.push(Cell::from("rel err"));
        }
        if self.metric_visibility.chi_sq {
            header.push(Cell::from("chi^2"));
        }
        if self.metric_visibility.max_weight_impact {
            header.push(Cell::from("mwi"));
        }

        let rows = (0..update.main_results.slot_headers.len())
            .map(|slot_index| {
                let mut cells = vec![
                    Cell::from(if slot_index == self.focused_slot {
                        format!(
                            "> {}",
                            update.main_results.slot_headers[slot_index].to_plain_string()
                        )
                    } else {
                        update.main_results.slot_headers[slot_index].to_plain_string()
                    }),
                    Cell::from(self.overview_value_text(update, slot_index)),
                ];
                if self.metric_visibility.relative_error {
                    cells.push(Cell::from(self.overview_metric_text(
                        update,
                        slot_index,
                        |cell| cell.relative_error.as_ref().map(|field| &field.display),
                    )));
                }
                if self.metric_visibility.chi_sq {
                    cells.push(Cell::from(self.overview_metric_text(
                        update,
                        slot_index,
                        |cell| cell.chi_sq.as_ref().map(|field| &field.display),
                    )));
                }
                if self.metric_visibility.max_weight_impact {
                    cells.push(Cell::from(self.overview_metric_text(
                        update,
                        slot_index,
                        |cell| cell.max_weight_impact.as_ref().map(|field| &field.display),
                    )));
                }
                Row::new(cells)
            })
            .collect::<Vec<_>>();

        let widths = match header.len() {
            2 => vec![Constraint::Length(24), Constraint::Min(20)],
            3 => vec![
                Constraint::Length(24),
                Constraint::Min(20),
                Constraint::Length(16),
            ],
            4 => vec![
                Constraint::Length(24),
                Constraint::Min(20),
                Constraint::Length(16),
                Constraint::Length(14),
            ],
            _ => vec![
                Constraint::Length(24),
                Constraint::Min(20),
                Constraint::Length(16),
                Constraint::Length(14),
                Constraint::Length(14),
            ],
        };
        let table = Table::new(rows, widths)
            .header(
                Row::new(header).style(
                    Style::default()
                        .fg(Color::Blue)
                        .add_modifier(Modifier::BOLD),
                ),
            )
            .block(titled_block(format!(
                "All-slot metrics [{} density]",
                self.density.label()
            )))
            .column_spacing(1);
        frame.render_widget(table, area);
    }

    fn draw_statistics_panel(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let Some(snapshot) = update.statistics_snapshot() else {
            frame.render_widget(
                Paragraph::new("Statistics unavailable.").block(titled_block("Statistics")),
                area,
            );
            return;
        };

        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(5),
                Constraint::Length(3),
                Constraint::Length(3),
                Constraint::Min(3),
            ])
            .split(area);

        let summary = vec![
            Line::from(format!("evals: {}", abbreviate_count(snapshot.num_evals))),
            Line::from(format!(
                "avg total {}   param {}   integ {}   eval {}",
                utils::format_evaluation_time_from_f64(snapshot.average_total_time_seconds),
                utils::format_evaluation_time_from_f64(
                    snapshot.average_parameterization_time_seconds
                ),
                utils::format_evaluation_time_from_f64(snapshot.average_integrand_time_seconds),
                utils::format_evaluation_time_from_f64(snapshot.average_evaluator_time_seconds),
            )),
            Line::from(format!(
                "obs {}   integrator {}   selection {}",
                utils::format_evaluation_time_from_f64(snapshot.average_observable_time_seconds),
                utils::format_evaluation_time_from_f64(snapshot.average_integrator_time_seconds),
                snapshot
                    .selection_efficiency_percentage
                    .map(|value| format!("{value:.2}%"))
                    .unwrap_or_else(|| "N/A".to_string()),
            )),
        ];
        frame.render_widget(
            Paragraph::new(Text::from(summary)).block(titled_block("Statistics")),
            layout[0],
        );

        frame.render_widget(
            Paragraph::new(Text::from(vec![stacked_percentage_line(
                "timings",
                &timing_percentages(&snapshot),
                area.width.saturating_sub(4) as usize,
            )]))
            .block(titled_block("Timing composition")),
            layout[1],
        );
        frame.render_widget(
            Paragraph::new(Text::from(vec![stacked_percentage_line(
                "precision",
                &precision_percentages(&snapshot),
                area.width.saturating_sub(4) as usize,
            )]))
            .block(titled_block("Precision mix")),
            layout[2],
        );
        frame.render_widget(
            Paragraph::new(Text::from(vec![stacked_percentage_line(
                "stability",
                &stability_percentages(&snapshot),
                area.width.saturating_sub(4) as usize,
            )]))
            .block(titled_block("Stability mix")),
            layout[3],
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

        let horizontal = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Percentage(55), Constraint::Percentage(45)])
            .split(area);

        let header = Row::new(["contribution", "cmp", "value", "rel err", "samples", "pdf"]).style(
            Style::default()
                .fg(Color::Blue)
                .add_modifier(Modifier::BOLD),
        );
        let rows = discrete_rows
            .iter()
            .map(|entry| {
                let cell = entry.row.slot_cell(0);
                Row::new(vec![
                    Cell::from(entry.row.contribution.display.to_plain_string()),
                    Cell::from(entry.row.component.display.to_plain_string()),
                    Cell::from(
                        cell.and_then(|cell| cell.value.as_ref())
                            .map(|value| value.display.to_plain_string())
                            .unwrap_or_default(),
                    ),
                    Cell::from(
                        cell.and_then(|cell| cell.relative_error.as_ref())
                            .map(|value| value.display.to_plain_string())
                            .unwrap_or_default(),
                    ),
                    Cell::from(
                        cell.and_then(|cell| cell.sample_count.as_ref())
                            .map(|value| value.display.to_plain_string())
                            .unwrap_or_default(),
                    ),
                    Cell::from(
                        cell.and_then(|cell| cell.target_pdf.as_ref())
                            .map(|value| value.display.to_plain_string())
                            .unwrap_or_default(),
                    ),
                ])
            })
            .collect::<Vec<_>>();

        let mut state = TableState::default().with_selected(Some(self.selected_discrete_row));
        let table = Table::new(
            rows,
            [
                Constraint::Length(24),
                Constraint::Length(6),
                Constraint::Min(18),
                Constraint::Length(12),
                Constraint::Length(10),
                Constraint::Length(9),
            ],
        )
        .header(header)
        .row_highlight_style(Style::default().bg(Color::Blue).fg(Color::Black))
        .block(titled_block(format!(
            "Discrete bins (sort: {:?} {})",
            self.discrete_sort,
            if self.discrete_descending {
                "desc"
            } else {
                "asc"
            }
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
        let mut lines = vec![Line::from(vec![
            Span::styled(
                "selected: ",
                Style::default()
                    .fg(Color::Blue)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::raw(row.contribution.display.to_plain_string()),
            Span::raw(" "),
            Span::raw(row.component.display.to_plain_string()),
        ])];

        for slot_index in 0..update.main_results.slot_headers.len() {
            let slot_name = update.main_results.slot_headers[slot_index].to_plain_string();
            let Some(cell) = row.slot_cell(slot_index) else {
                continue;
            };
            lines.push(Line::from(""));
            lines.push(Line::from(vec![Span::styled(
                slot_name,
                Style::default()
                    .fg(Color::Green)
                    .add_modifier(Modifier::BOLD),
            )]));
            if let Some(value) = cell.value.as_ref() {
                lines.push(Line::from(format!(
                    "value    {}",
                    value.display.to_plain_string()
                )));
            }
            if let Some(value) = cell.relative_error.as_ref() {
                lines.push(Line::from(format!(
                    "rel err  {}",
                    value.display.to_plain_string()
                )));
            }
            if let Some(value) = cell.chi_sq.as_ref() {
                lines.push(Line::from(format!(
                    "chi^2    {}",
                    value.display.to_plain_string()
                )));
            }
            if let Some(value) = cell.max_weight_impact.as_ref() {
                lines.push(Line::from(format!(
                    "mwi      {}",
                    value.display.to_plain_string()
                )));
            }
            if let Some(value) = cell.sample_fraction.as_ref() {
                lines.push(Line::from(format!(
                    "samples  {}",
                    value.display.to_plain_string()
                )));
            }
            if let Some(value) = cell.target_pdf.as_ref() {
                lines.push(Line::from(format!(
                    "pdf      {}",
                    value.display.to_plain_string()
                )));
            }
        }

        frame.render_widget(
            Paragraph::new(Text::from(lines))
                .block(titled_block("Selected bin detail"))
                .wrap(Wrap { trim: false }),
            area,
        );
    }

    fn draw_max_weight_tab(&self, frame: &mut Frame<'_>, area: Rect, update: &StatusUpdate) {
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([Constraint::Percentage(45), Constraint::Percentage(55)])
            .split(area);

        if let Some(section) = update.max_weight_details.as_ref() {
            let rows = section
                .rows_by_slot
                .iter()
                .flat_map(|group| group.iter())
                .map(|row| {
                    Row::new(vec![
                        Cell::from(row.slot.display.to_plain_string()),
                        Cell::from(row.component_sign.display.to_plain_string()),
                        Cell::from(row.max_eval.display.to_plain_string()),
                        Cell::from(row.coordinates.display.to_plain_string()),
                    ])
                })
                .collect::<Vec<_>>();
            let table = Table::new(
                rows,
                [
                    Constraint::Length(24),
                    Constraint::Length(10),
                    Constraint::Length(20),
                    Constraint::Min(24),
                ],
            )
            .header(
                Row::new(["slot", "cmp", "max eval", "coordinates"]).style(
                    Style::default()
                        .fg(Color::Blue)
                        .add_modifier(Modifier::BOLD),
                ),
            )
            .block(titled_block("Overall max weight"));
            frame.render_widget(table, layout[0]);
        } else {
            frame.render_widget(
                Paragraph::new("No overall max-weight data.")
                    .block(titled_block("Overall max weight")),
                layout[0],
            );
        }

        if let Some(section) = update.discrete_max_weight_details.as_ref() {
            let rows = section
                .row_groups
                .iter()
                .flat_map(|group| group.iter())
                .map(|row| {
                    let values = row
                        .slot_values
                        .iter()
                        .map(|value| {
                            value
                                .as_ref()
                                .map(|field| field.display.to_plain_string())
                                .unwrap_or_default()
                        })
                        .collect::<Vec<_>>()
                        .join(" | ");
                    let coordinates = row
                        .slot_coordinates
                        .iter()
                        .map(|entry| {
                            format!(
                                "{}: {}",
                                entry.slot.display.to_plain_string(),
                                entry.coordinates.display.to_plain_string()
                            )
                        })
                        .collect::<Vec<_>>()
                        .join("\n");
                    Row::new(vec![
                        Cell::from(row.contribution.display.to_plain_string()),
                        Cell::from(row.component_sign.display.to_plain_string()),
                        Cell::from(values),
                        Cell::from(coordinates),
                    ])
                })
                .collect::<Vec<_>>();
            let table = Table::new(
                rows,
                [
                    Constraint::Length(24),
                    Constraint::Length(10),
                    Constraint::Min(18),
                    Constraint::Min(28),
                ],
            )
            .header(
                Row::new(["contribution", "cmp", "slot values", "coordinates"]).style(
                    Style::default()
                        .fg(Color::Blue)
                        .add_modifier(Modifier::BOLD),
                ),
            )
            .block(titled_block("Per-bin max weight"));
            frame.render_widget(table, layout[1]);
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
            Span::styled(
                "Tabs",
                Style::default()
                    .fg(Color::Blue)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::raw(" 1/2/3 or <- ->  "),
            Span::styled(
                "Slot",
                Style::default()
                    .fg(Color::Blue)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::raw(" [ / ]  "),
            Span::styled(
                "Bins",
                Style::default()
                    .fg(Color::Blue)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::raw(" j/k  "),
            Span::styled(
                "Metrics",
                Style::default()
                    .fg(Color::Blue)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::raw(" r c w  "),
            Span::styled(
                "Sort",
                Style::default()
                    .fg(Color::Blue)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::raw(" s/v  "),
            Span::styled(
                "Density",
                Style::default()
                    .fg(Color::Blue)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::raw(" m  "),
            Span::styled(
                "Help",
                Style::default()
                    .fg(Color::Blue)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::raw(" ?  "),
            Span::styled(
                "Abort",
                Style::default().fg(Color::Red).add_modifier(Modifier::BOLD),
            ),
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
                Line::from("[ and ]             change focused slot"),
                Line::from("j / k               move selected discrete row"),
                Line::from("s                   cycle discrete sort"),
                Line::from("v                   reverse discrete sort order"),
                Line::from("r / c / w           toggle rel err, chi^2, mwi columns"),
                Line::from("m                   cycle density mode"),
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

    fn discrete_rows(&self) -> Vec<DiscreteRowRef<'_>> {
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
            let lhs_key = lhs.sort_value(self.discrete_sort);
            let rhs_key = rhs.sort_value(self.discrete_sort);
            let ordering = match self.discrete_sort {
                ContributionSortMode::Index => match (lhs.contribution(), rhs.contribution()) {
                    (ContributionKind::Bin(lhs_index), ContributionKind::Bin(rhs_index)) => {
                        lhs_index
                            .cmp(&rhs_index)
                            .then(lhs.component().cmp(&rhs.component()))
                    }
                    _ => Ordering::Equal,
                },
                ContributionSortMode::Integral | ContributionSortMode::Error => rhs_key
                    .partial_cmp(&lhs_key)
                    .unwrap_or(Ordering::Equal)
                    .then(lhs.component().cmp(&rhs.component())),
            };
            if self.discrete_descending {
                ordering
            } else {
                ordering.reverse()
            }
        });
        rows
    }

    fn overview_value_text(&self, update: &StatusUpdate, slot_index: usize) -> Text<'static> {
        overview_metric_text(update, slot_index, |cell| {
            cell.value.as_ref().map(|value| &value.display)
        })
    }

    fn overview_metric_text(
        &self,
        update: &StatusUpdate,
        slot_index: usize,
        select: impl Fn(&MainTableSlotCells) -> Option<&StyledText>,
    ) -> Text<'static> {
        overview_metric_text(update, slot_index, select)
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
            let component = match row.component.raw {
                ComponentKind::Real => "re",
                ComponentKind::Imag => "im",
            };
            let text = row.slot_cell(slot_index).and_then(select)?;
            let mut spans = vec![Span::styled(
                format!("{component}: "),
                Style::default()
                    .fg(Color::Green)
                    .add_modifier(Modifier::BOLD),
            )];
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

fn style_from_text_style(style: TextStyle) -> Style {
    let mut rendered = Style::default();
    rendered = match style.color {
        Some(TextColor::Green) => rendered.fg(Color::Green),
        Some(TextColor::Blue) => rendered.fg(Color::Blue),
        Some(TextColor::Red) => rendered.fg(Color::Red),
        None => rendered,
    };
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

fn titled_block(title: impl Into<String>) -> Block<'static> {
    Block::default()
        .borders(Borders::ALL)
        .title(title.into())
        .border_style(Style::default().fg(Color::Blue))
}

fn abbreviate_count(value: usize) -> String {
    super::display::format_abbreviated_count(value)
}

fn format_rate(value: Option<f64>) -> String {
    let Some(value) = value else {
        return "N/A".to_string();
    };
    if !value.is_finite() || value <= 0.0 {
        return "N/A".to_string();
    }
    format!("{}/s", abbreviate_count(value.round() as usize))
}

fn stacked_percentage_line(
    label: &str,
    segments: &[(String, f64, Color)],
    width: usize,
) -> Line<'static> {
    let prefix = format!("{label:<10}");
    let bar_width = width.saturating_sub(prefix.len() + 12).max(8);
    let mut spans = vec![Span::styled(
        prefix,
        Style::default()
            .fg(Color::Blue)
            .add_modifier(Modifier::BOLD),
    )];
    let mut used = 0usize;
    for (index, (_, percentage, color)) in segments.iter().enumerate() {
        let width = if index + 1 == segments.len() {
            bar_width.saturating_sub(used)
        } else {
            let segment = ((*percentage / 100.0) * bar_width as f64).round() as usize;
            used += segment;
            segment
        };
        spans.push(Span::styled("█".repeat(width), Style::default().fg(*color)));
    }
    spans.push(Span::raw(" "));
    spans.push(Span::raw(
        segments
            .iter()
            .map(|(name, percentage, _)| format!("{name}:{percentage:.1}%"))
            .collect::<Vec<_>>()
            .join(" "),
    ));
    Line::from(spans)
}

fn timing_percentages(snapshot: &IntegrationStatisticsSnapshot) -> Vec<(String, f64, Color)> {
    let total = snapshot.average_total_time_seconds.max(f64::EPSILON);
    vec![
        (
            "param".to_string(),
            snapshot.average_parameterization_time_seconds / total * 100.0,
            Color::Cyan,
        ),
        (
            "integrand".to_string(),
            snapshot.average_integrand_time_seconds / total * 100.0,
            Color::Green,
        ),
        (
            "evaluator".to_string(),
            snapshot.average_evaluator_time_seconds / total * 100.0,
            Color::Blue,
        ),
        (
            "obs".to_string(),
            snapshot.average_observable_time_seconds / total * 100.0,
            Color::Magenta,
        ),
        (
            "overhead".to_string(),
            snapshot.average_integrator_time_seconds / total * 100.0,
            Color::Yellow,
        ),
    ]
}

fn precision_percentages(snapshot: &IntegrationStatisticsSnapshot) -> Vec<(String, f64, Color)> {
    vec![
        ("f64".to_string(), snapshot.f64_percentage, Color::Green),
        ("f128".to_string(), snapshot.f128_percentage, Color::Blue),
        ("arb".to_string(), snapshot.arb_percentage, Color::Yellow),
    ]
}

fn stability_percentages(snapshot: &IntegrationStatisticsSnapshot) -> Vec<(String, f64, Color)> {
    let unstable = (snapshot.nan_or_unstable_percentage - snapshot.nan_percentage).max(0.0);
    let stable = (100.0 - snapshot.nan_or_unstable_percentage).max(0.0);
    vec![
        ("stable".to_string(), stable, Color::Green),
        ("unstable".to_string(), unstable, Color::Yellow),
        ("nan".to_string(), snapshot.nan_percentage, Color::Red),
    ]
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
