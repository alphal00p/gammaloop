//! Terminal presentation of the typed monitor model.

use std::io::{self, IsTerminal, Stderr, Write};
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Duration;

use crossterm::{
    cursor,
    event::{self, Event, KeyCode, KeyEventKind, KeyModifiers},
    execute,
    terminal::{
        EnterAlternateScreen, LeaveAlternateScreen, disable_raw_mode, enable_raw_mode,
        is_raw_mode_enabled,
    },
};
use ratatui::{
    Frame, Terminal,
    backend::CrosstermBackend,
    layout::{Constraint, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, BorderType, Clear, Paragraph, Wrap},
};

use super::model::{
    DfsModel, EstimateState, MonitorModel, OrderedModel, PhaseKind, PhaseRecord, StageKind, format,
};
use super::{ApiStatus, MonitorState, MonitorStatus, ScopeKind, UiAction, UiView};

// A process-local UI claim, never a reducer or observer lock.
static TERMINAL_CLAIMED: AtomicBool = AtomicBool::new(false);

struct TerminalClaim;

impl TerminalClaim {
    fn acquire() -> io::Result<Self> {
        TERMINAL_CLAIMED
            .compare_exchange(false, true, Ordering::AcqRel, Ordering::Acquire)
            .map(|_| Self)
            .map_err(|_| {
                io::Error::new(
                    io::ErrorKind::AlreadyExists,
                    "another monitor already owns this process's terminal",
                )
            })
    }
}

impl Drop for TerminalClaim {
    fn drop(&mut self) {
        TERMINAL_CLAIMED.store(false, Ordering::Release);
    }
}

/// The caller owns this guard on the terminal thread, never the native callback.
pub(super) struct TerminalGuard {
    terminal: Terminal<CrosstermBackend<Stderr>>,
    restore_raw_mode: bool,
    alternate_screen: bool,
    cursor_hidden: bool,
    // Last field: release only after our restoration and the terminal/backend Drop.
    _claim: TerminalClaim,
}

impl TerminalGuard {
    pub(super) fn enter() -> io::Result<Self> {
        let claim = TerminalClaim::acquire()?;
        if !io::stdin().is_terminal()
            || !io::stderr().is_terminal()
            || std::env::var("TERM").ok().as_deref() == Some("dumb")
        {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "the monitor requires terminal input and stderr",
            ));
        }
        #[cfg(unix)]
        {
            let attributes = rustix::termios::tcgetattr(io::stdin())?;
            let ordinary = rustix::termios::LocalModes::ICANON | rustix::termios::LocalModes::ECHO;
            if !attributes.local_modes.contains(ordinary) {
                return Err(io::Error::new(
                    io::ErrorKind::AlreadyExists,
                    "the monitor will not take over an externally raw or cbreak terminal",
                ));
            }
        }
        // Another UI may own both raw mode and the alternate screen. Refuse
        // before touching either; leaving a nested screen would disrupt it.
        if is_raw_mode_enabled()? {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists,
                "the monitor will not take over an already-raw terminal",
            ));
        }

        let terminal = Terminal::new(CrosstermBackend::new(io::stderr()))?;
        let mut guard = Self {
            terminal,
            restore_raw_mode: false,
            alternate_screen: false,
            cursor_hidden: false,
            _claim: claim,
        };
        // Mark attempted transitions first: Drop also rolls back partial entry.
        guard.restore_raw_mode = true;
        enable_raw_mode()?;
        guard.alternate_screen = true;
        execute!(guard.terminal.backend_mut(), EnterAlternateScreen)?;
        guard.cursor_hidden = true;
        execute!(guard.terminal.backend_mut(), cursor::Hide)?;
        guard.terminal.clear()?;
        Ok(guard)
    }

    pub(super) fn draw(
        &mut self,
        state: &mut MonitorState,
        now: Duration,
        no_color: bool,
    ) -> io::Result<()> {
        let model = state.model(now);
        let maximum_scroll = self.draw_model(&model, state.ui(), no_color)?;
        // Keep the next key action relative to the visible position. Recompute
        // every successful frame, including resizes and shrinking observations.
        state.clamp_scroll(maximum_scroll);
        Ok(())
    }

    /// Draw one model; return the scroll extent of the details view.
    pub(super) fn draw_model(
        &mut self,
        model: &MonitorModel,
        ui: UiView,
        no_color: bool,
    ) -> io::Result<usize> {
        let mut maximum_scroll = 0;
        self.terminal.draw(|frame| {
            maximum_scroll = draw(frame, model, ui, no_color);
        })?;
        Ok(maximum_scroll)
    }
}

impl Drop for TerminalGuard {
    fn drop(&mut self) {
        if self.cursor_hidden || self.alternate_screen {
            let _ = execute!(
                self.terminal.backend_mut(),
                crossterm::style::ResetColor,
                cursor::Show
            );
        }
        if self.alternate_screen {
            let _ = execute!(self.terminal.backend_mut(), LeaveAlternateScreen);
        }
        if self.restore_raw_mode {
            let _ = disable_raw_mode();
        }
        let _ = self.terminal.backend_mut().flush();
    }
}

/// Raw mode turns Ctrl-C into a key. After the display has restored the
/// terminal, deliver SIGINT to the foreground process group as the terminal
/// itself would have done.
pub(super) fn interrupt_foreground_group() {
    #[cfg(unix)]
    let _ = rustix::process::kill_current_process_group(rustix::process::Signal::INT);
}

/// Poll once on the terminal-owner thread; no input loop runs in the observer.
pub(super) fn poll_action(timeout: Duration) -> io::Result<Option<UiAction>> {
    if !event::poll(timeout)? {
        return Ok(None);
    }
    let Event::Key(key) = event::read()? else {
        return Ok(None);
    };
    if key.kind == KeyEventKind::Release {
        return Ok(None);
    }
    if key.modifiers.contains(KeyModifiers::CONTROL) {
        return Ok(matches!(key.code, KeyCode::Char('c' | 'C'))
            .then_some(UiAction::Interrupt)
            .filter(|_| key.kind == KeyEventKind::Press));
    }
    if key
        .modifiers
        .intersects(KeyModifiers::ALT | KeyModifiers::SUPER)
    {
        return Ok(None);
    }
    let press = key.kind == KeyEventKind::Press;
    let action = match key.code {
        KeyCode::Char('q') if press => UiAction::Detach,
        KeyCode::Char(' ') if press => UiAction::ToggleFreeze,
        KeyCode::Char('?') if press => UiAction::ToggleHelp,
        KeyCode::Char('d') if press => UiAction::ToggleDetails,
        KeyCode::Up => UiAction::ScrollUp,
        KeyCode::Down => UiAction::ScrollDown,
        KeyCode::PageUp => UiAction::PageUp,
        KeyCode::PageDown => UiAction::PageDown,
        _ => return Ok(None),
    };
    Ok(Some(action))
}

/// Semantic styles; the renderer never inspects label or value text.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Role {
    Label,
    Plain,
    Muted,
    Count,
    Time,
    Rate,
    Memory,
    Good,
    Warn,
    Bad,
    Title,
    Rule,
}

const BLUE: Color = Color::Rgb(138, 188, 255);
const CYAN: Color = Color::Rgb(115, 215, 225);
const GREEN: Color = Color::Rgb(125, 215, 165);
const AMBER: Color = Color::Rgb(240, 200, 120);
const RED: Color = Color::Rgb(240, 120, 120);
const MAGENTA: Color = Color::Rgb(210, 165, 235);

fn style(role: Role, no_color: bool) -> Style {
    let (color, modifier) = match role {
        Role::Label => (Color::Gray, Modifier::empty()),
        Role::Plain => (Color::Reset, Modifier::empty()),
        Role::Muted => (Color::DarkGray, Modifier::empty()),
        Role::Count => (Color::Reset, Modifier::BOLD),
        Role::Time => (CYAN, Modifier::BOLD),
        Role::Rate => (GREEN, Modifier::BOLD),
        Role::Memory => (BLUE, Modifier::BOLD),
        Role::Good => (GREEN, Modifier::BOLD),
        Role::Warn => (AMBER, Modifier::BOLD),
        Role::Bad => (RED, Modifier::BOLD),
        Role::Title => (MAGENTA, Modifier::BOLD),
        Role::Rule => (Color::DarkGray, Modifier::empty()),
    };
    let style = Style::default().add_modifier(modifier);
    if no_color || color == Color::Reset {
        style
    } else {
        style.fg(color)
    }
}

/// Spans of one line, built with semantic roles.
struct Text {
    spans: Vec<(String, Role)>,
}

impl Text {
    fn new() -> Self {
        Self { spans: Vec::new() }
    }

    fn push(&mut self, text: impl Into<String>, role: Role) -> &mut Self {
        self.spans.push((text.into(), role));
        self
    }

    /// A ` · ` separator unless the line is empty.
    fn sep(&mut self) -> &mut Self {
        if !self.spans.is_empty() {
            self.push(" · ", Role::Muted);
        }
        self
    }

    fn width(&self) -> usize {
        self.spans
            .iter()
            .map(|(text, _)| Span::raw(text.as_str()).width())
            .sum()
    }

    fn line(self, no_color: bool) -> Line<'static> {
        Line::from(
            self.spans
                .into_iter()
                .map(|(text, role)| Span::styled(text, style(role, no_color)))
                .collect::<Vec<_>>(),
        )
    }
}

const LABEL_WIDTH: usize = 10;

/// `Label     value…    right`, the right part flush right when it fits.
fn labelled(
    label: &str,
    value: Text,
    right: Option<Text>,
    width: u16,
    no_color: bool,
) -> Line<'static> {
    let mut text = Text::new();
    text.push(format!(" {label:<w$}", w = LABEL_WIDTH - 1), Role::Label);
    let used = LABEL_WIDTH + value.width();
    text.spans.extend(value.spans);
    if let Some(right) = right {
        // One free column at the right edge; a separator when it is tight.
        let total = used + right.width() + 1;
        if total + 2 <= usize::from(width) {
            text.push(" ".repeat(usize::from(width) - total), Role::Plain);
        } else {
            text.sep();
        }
        text.spans.extend(right.spans);
    }
    text.line(no_color)
}

fn rule(title: &str, width: u16, no_color: bool) -> Line<'static> {
    let head = format!(" ── {title} ");
    let fill = usize::from(width).saturating_sub(Span::raw(head.as_str()).width() + 1);
    Line::from(vec![
        Span::styled(head, style(Role::Rule, no_color)),
        Span::styled("─".repeat(fill), style(Role::Rule, no_color)),
    ])
}

fn engine(model: &MonitorModel) -> Option<&'static str> {
    if model.ordered.is_some() {
        Some("Ordered")
    } else if model.dfs.is_some() {
        Some("DFS")
    } else {
        None
    }
}

fn header(model: &MonitorModel, ui: UiView, width: u16, no_color: bool) -> Line<'static> {
    let mut left = Text::new();
    left.push(" vakint monitor", Role::Title);
    if let Some(engine) = engine(model) {
        left.sep().push(engine, Role::Plain);
    }
    left.sep().push(
        match model.scope {
            ScopeKind::Api => "call",
            ScopeKind::Job => "job",
        },
        Role::Muted,
    );
    if model.scope == ScopeKind::Job && model.api_status != ApiStatus::NotStarted {
        left.push(
            match model.api_status {
                ApiStatus::NotStarted => "",
                ApiStatus::Running => " (API running)",
                ApiStatus::Succeeded => " (API done)",
                ApiStatus::Failed => " (API failed)",
            },
            Role::Muted,
        );
    }
    let mut right = Text::new();
    if ui.frozen {
        right.push("display frozen", Role::Warn).sep();
    }
    let (status, role) = match model.status {
        MonitorStatus::Running => ("● RUNNING", Role::Time),
        MonitorStatus::Succeeded => ("✔ SUCCEEDED", Role::Good),
        MonitorStatus::Failed => ("✖ FAILED", Role::Bad),
        MonitorStatus::Partial => ("◐ PARTIAL", Role::Warn),
    };
    right
        .push(status, role)
        .push("  ", Role::Plain)
        .push(format::duration(model.elapsed_seconds), Role::Time)
        .push(" ", Role::Plain);
    let gap = usize::from(width)
        .saturating_sub(left.width() + right.width())
        .max(1);
    left.push(" ".repeat(gap), Role::Plain);
    left.spans.extend(right.spans);
    left.line(no_color)
}

fn run_line(model: &MonitorModel, width: u16, no_color: bool) -> Option<Line<'static>> {
    let run = model.run.as_ref()?;
    let mut value = Text::new();
    value
        .push(&run.arithmetic, Role::Plain)
        .sep()
        .push(run.n_cores.to_string(), Role::Count)
        .push(
            if run.n_cores == 1 { " core" } else { " cores" },
            Role::Plain,
        )
        .sep()
        .push("step cache ", Role::Muted)
        .push(&run.step_cache, Role::Plain);
    Some(labelled("Run", value, None, width, no_color))
}

fn stage_line(model: &MonitorModel, width: u16, no_color: bool) -> Line<'static> {
    let stage = &model.stage;
    let mut value = Text::new();
    if stage.id.is_empty() {
        value.push("starting", Role::Muted);
    } else {
        value.push(stage.kind.label(), Role::Plain);
        let detail = match stage.kind {
            StageKind::Native => stage.activity.clone(),
            _ => Some(super::model::humanize(&stage.id)),
        };
        if let Some(detail) = detail.filter(|detail| detail != stage.kind.label()) {
            value.push(" ▸ ", Role::Muted).push(detail, Role::Count);
        }
        if let Some(work) = &stage.work {
            value
                .sep()
                .push(format!("{}/{}", work.completed, work.total), Role::Count)
                .push(format!(" {}", work.label), Role::Plain);
        }
    }
    let mut right = Text::new();
    right.push(format::duration(stage.seconds), Role::Time);
    labelled("Stage", value, Some(right), width, no_color)
}

fn batch_line(model: &MonitorModel, width: u16, no_color: bool) -> Option<Line<'static>> {
    if model.batch.is_none() && model.cohort.is_none() {
        return None;
    }
    let mut value = Text::new();
    if let Some(batch) = &model.batch {
        // A batch counter counts finished batches; show the running one.
        let current = (batch.completed + 1).min(batch.total);
        value
            .push(format!("{current}/{}", batch.total), Role::Count)
            .push(format!(" {}", batch.label), Role::Plain);
    }
    if let Some(cohort) = &model.cohort {
        let current = (cohort.completed + 1).min(cohort.total);
        value
            .sep()
            .push("cohort ", Role::Muted)
            .push(format!("{current}/{}", cohort.total), Role::Count);
    }
    Some(labelled("Batch", value, None, width, no_color))
}

fn phase_name(kind: PhaseKind, support: Option<usize>, owner: Option<usize>) -> String {
    let mut name = kind.label().to_owned();
    if kind == PhaseKind::Done {
        return name;
    }
    if let Some(support) = support {
        name.push_str(&format!(" s{support}"));
    }
    if let Some(owner) = owner {
        name.push_str(&format!(" o{owner}"));
    }
    name
}

fn replay(hits: u64, misses: u64) -> Option<String> {
    format::percent(hits, hits + misses)
}

fn ordered_lines(ordered: &OrderedModel, width: u16, no_color: bool) -> Vec<Line<'static>> {
    let phase = &ordered.phase;
    let mut lines = Vec::new();

    let mut value = Text::new();
    value
        .push(
            format!("{}/{}", phase.index.min(phase.total), phase.total),
            Role::Count,
        )
        .push(" ", Role::Plain)
        .push(
            phase_name(phase.kind, phase.support, phase.owner),
            Role::Title,
        );
    if phase.kind != PhaseKind::Done {
        value
            .sep()
            .push(format::count(phase.keys), Role::Count)
            .push(" keys", Role::Plain);
        if let Some(rate) = phase
            .seconds
            .filter(|seconds| *seconds >= MIN_RATE_SECONDS)
            .map(|seconds| phase.keys as f64 / seconds)
        {
            value.sep().push(format::rate(rate), Role::Rate);
        }
    }
    let right = phase.seconds.map(|seconds| {
        let mut right = Text::new();
        right
            .push(format::duration(seconds), Role::Time)
            .push(" in phase", Role::Muted);
        right
    });
    lines.push(labelled("Phase", value, right, width, no_color));

    let mut value = Text::new();
    match ordered.recent_rate {
        Some(rate) => {
            value
                .push(format::rate(rate), Role::Rate)
                .push(" now", Role::Plain);
        }
        // Only a working core has a current rate.
        None if phase.seconds.is_some() => {
            value.push("—", Role::Muted).push(" now", Role::Plain);
        }
        None => {}
    }
    if let Some(rate) = ordered.request_rate {
        value
            .sep()
            .push(format::rate(rate), Role::Rate)
            .push(" request", Role::Plain);
    }
    value
        .sep()
        .push(format::count(ordered.processed_keys), Role::Count)
        .push(" keys done", Role::Plain);
    lines.push(labelled("Rate", value, None, width, no_color));

    let frontier = &ordered.frontier;
    let mut value = Text::new();
    value
        .push(format::count(frontier.keys), Role::Count)
        .push(" keys", Role::Plain)
        .push(
            format!(" (peak {})", format::count(frontier.peak_keys)),
            Role::Muted,
        )
        .sep()
        .push(format::count(frontier.monomials), Role::Count)
        .push(" mono", Role::Plain);
    if let Some(per_key) = frontier.monomials.checked_div(frontier.keys) {
        value.push(format!(" ({per_key}/key)"), Role::Muted);
    }
    lines.push(labelled("Frontier", value, None, width, no_color));

    let output = &ordered.output;
    let mut value = Text::new();
    value
        .push(format::count(output.keys), Role::Count)
        .push(" keys", Role::Plain)
        .sep()
        .push(format::count(output.monomials), Role::Count)
        .push(" mono", Role::Plain)
        .sep()
        .push(format::count(ordered.input_keys), Role::Count)
        .push(" input keys", Role::Muted);
    lines.push(labelled("Output", value, None, width, no_color));

    let grouping = &ordered.grouping;
    let mut value = Text::new();
    value
        .push(format::count(grouping.emitted), Role::Count)
        .push(" emitted", Role::Plain);
    if let Some(merged) = format::percent(grouping.merged, grouping.emitted) {
        value
            .sep()
            .push(merged, Role::Good)
            .push(" merged", Role::Plain);
    }
    value
        .sep()
        .push(format::count(grouping.cancelled), Role::Count)
        .push(" cancelled", Role::Plain);
    lines.push(labelled("Grouping", value, None, width, no_color));

    // The Run line names the cache mode; this line shows its activity only.
    if let Some(share) = replay(ordered.step_cache_hits, ordered.step_cache_misses) {
        let mut value = Text::new();
        value
            .push("step cache ", Role::Muted)
            .push(share, Role::Good)
            .push(" replayed", Role::Plain)
            .push(
                format!(
                    " ({} / {})",
                    format::count(ordered.step_cache_hits),
                    format::count(ordered.step_cache_misses)
                ),
                Role::Muted,
            );
        if let Some(share) = replay(phase.step_cache_hits, phase.step_cache_misses) {
            value
                .sep()
                .push(share, Role::Good)
                .push(" this phase", Role::Plain);
        }
        lines.push(labelled("Cache", value, None, width, no_color));
    }
    lines
}

fn eta_text(value: &mut Text, seconds: Option<f64>, state: EstimateState) {
    match (state, seconds) {
        (EstimateState::Live, Some(seconds)) => {
            value
                .push("ETA ", Role::Muted)
                .push(format::duration(seconds), Role::Time);
        }
        (EstimateState::Live, None) => {
            value.push("ETA —", Role::Muted);
        }
        (EstimateState::Warmup, _) => {
            value.push("ETA warming up", Role::Muted);
        }
        (EstimateState::Stale, _) => {
            value.push("ETA stale", Role::Warn);
        }
        (EstimateState::Inactive, _) => {
            value.push("ETA indeterminate", Role::Muted);
        }
    }
}

fn dfs_lines(dfs: &DfsModel, width: u16, no_color: bool) -> Vec<Line<'static>> {
    let estimates = &dfs.estimates;
    let mut lines = Vec::new();
    let count = |value: Option<u64>| value.map_or_else(|| "—".to_owned(), format::count);

    let mut value = Text::new();
    value.push(
        format!("{}/{}", count(dfs.roots_completed), count(dfs.roots_total)),
        Role::Count,
    );
    if let Some(rate) = estimates.root_rate {
        value.sep().push(format::rate(rate), Role::Rate);
    }
    value.sep();
    eta_text(&mut value, estimates.root_eta_seconds, estimates.state);
    lines.push(labelled("Roots", value, None, width, no_color));

    let mut value = Text::new();
    value
        .push(format::count(dfs.entered), Role::Count)
        .push(" entered", Role::Plain)
        .sep()
        .push(count(dfs.completed), Role::Count)
        .push(" completed", Role::Plain)
        .sep()
        .push(count(dfs.pending_unique), Role::Warn)
        .push(
            if dfs.pending_exact {
                " pending"
            } else {
                " pending (sampled)"
            },
            Role::Plain,
        );
    lines.push(labelled("Nodes", value, None, width, no_color));

    let mut value = Text::new();
    match estimates.completion_rate {
        Some(rate) => value.push(format::rate(rate), Role::Rate),
        None => value.push("—", Role::Muted),
    };
    if let Some(growth) = estimates.growth_rate {
        value
            .sep()
            .push("frontier ", Role::Muted)
            .push(format!("{growth:+.2}/s"), Role::Count);
    }
    value.sep();
    eta_text(&mut value, estimates.frontier_eta_seconds, estimates.state);
    lines.push(labelled("Rate", value, None, width, no_color));

    let mut value = Text::new();
    value
        .push("persistent ", Role::Muted)
        .push(count(dfs.persistent_hits), Role::Good)
        .push(" hits / ", Role::Plain)
        .push(count(dfs.persistent_misses), Role::Count)
        .push(" misses", Role::Plain)
        .sep()
        .push("RAM ", Role::Muted)
        .push(count(dfs.ram_hits), Role::Good)
        .push(" hits", Role::Plain);
    lines.push(labelled("Cache", value, None, width, no_color));

    if let Some(collection) = &dfs.collection {
        let mut value = Text::new();
        value
            .push(
                format!(
                    "{}/{}",
                    count(collection.completed_roots),
                    count(collection.total_roots)
                ),
                Role::Count,
            )
            .push(" roots", Role::Plain);
        if let (Some(item), Some(items)) = (collection.contribution, collection.contributions) {
            value
                .sep()
                .push("assembling ", Role::Muted)
                .push(format!("{item}/{items}"), Role::Count);
        }
        if let Some(age) = collection.age_seconds {
            value
                .sep()
                .push(format!("age {}", format::duration(age)), Role::Muted);
        }
        lines.push(labelled("Collect", value, None, width, no_color));
    }
    lines
}

fn memory_line(model: &MonitorModel, width: u16, no_color: bool) -> Line<'static> {
    let memory = &model.memory;
    let bytes = |value: Option<u64>| value.map_or_else(|| "—".to_owned(), format::bytes);
    let mut value = Text::new();
    if !memory.verified && memory.available_bytes.is_none() {
        value.push("waiting for the first sample", Role::Muted);
        return labelled("Memory", value, None, width, no_color);
    }
    if memory.verified {
        value
            .push("RSS ", Role::Muted)
            .push(bytes(memory.rss_bytes), Role::Memory)
            .push(
                format!(" (peak {})", bytes(memory.peak_rss_bytes)),
                Role::Muted,
            );
    } else {
        value.push("RSS unverified", Role::Muted);
    }
    value
        .sep()
        .push("avail ", Role::Muted)
        .push(bytes(memory.available_bytes), Role::Memory);
    if let Some(cpu) = memory.cpu_percent {
        value
            .sep()
            .push("CPU ", Role::Muted)
            .push(format!("{cpu:.0}%"), Role::Count);
    }
    if let Some(guard) = &memory.guard
        && guard.status != "active"
    {
        value
            .sep()
            .push(format!("guard {}", guard.status), Role::Warn);
    }
    labelled("Memory", value, None, width, no_color)
}

fn phase_table_header(width: u16, no_color: bool) -> Line<'static> {
    let text = format!(
        "  {:>4}  {:<14}{:>9}{:>9}{:>9}{:>10}{:>9}{:>8}",
        "#", "phase", "time", "keys", "rate", "frontier", "mono", "cache"
    );
    Line::from(Span::styled(
        truncate(&text, width),
        style(Role::Label, no_color),
    ))
}

fn truncate(text: &str, width: u16) -> String {
    text.chars().take(usize::from(width)).collect()
}

/// Shortest interval whose keys per second are shown.
const MIN_RATE_SECONDS: f64 = 0.5;

/// The state of a phase-log row.
#[derive(Clone, Copy, PartialEq, Eq)]
enum RowState {
    Finished,
    /// The phase the native core is working on.
    Running,
    /// An unfinished phase while the native core is not working.
    Idle,
}

fn phase_row(record: &PhaseRecord, state: RowState, width: u16, no_color: bool) -> Line<'static> {
    let running = state != RowState::Finished;
    let rate = (record.seconds >= MIN_RATE_SECONDS).then(|| record.keys as f64 / record.seconds);
    let cache = replay(record.step_cache_hits, record.step_cache_misses).unwrap_or_default();
    let cells = [
        (format!("  {:>4}  ", record.index), Role::Muted),
        (
            format!(
                "{:<14}",
                phase_name(record.kind, record.support, record.owner)
            ),
            Role::Plain,
        ),
        (
            format!("{:>9}", format::duration(record.seconds)),
            Role::Time,
        ),
        (format!("{:>9}", format::count(record.keys)), Role::Count),
        (
            format!("{:>9}", rate.map_or_else(|| "—".to_owned(), format::rate)),
            Role::Rate,
        ),
        if running {
            match state {
                RowState::Running => (format!("{:>10}", "running"), Role::Warn),
                _ => (format!("{:>10}", "idle"), Role::Muted),
            }
        } else {
            (
                format!("{:>10}", format::count(record.frontier_keys)),
                Role::Count,
            )
        },
        if running {
            (format!("{:>9}", ""), Role::Plain)
        } else {
            (
                format!("{:>9}", format::count(record.frontier_monomials)),
                Role::Count,
            )
        },
        (format!("{cache:>8}"), Role::Good),
    ];
    let mut spans = Vec::new();
    let mut used = 0;
    for (text, role) in cells {
        let text = truncate(&text, width.saturating_sub(used as u16));
        used += Span::raw(text.as_str()).width();
        spans.push(Span::styled(text, style(role, no_color)));
    }
    Line::from(spans)
}

/// The finished phases of the request and the running one, newest last.
fn phase_rows(ordered: &OrderedModel) -> Vec<(PhaseRecord, RowState)> {
    let mut rows = ordered
        .phases
        .iter()
        .cloned()
        .map(|record| (record, RowState::Finished))
        .collect::<Vec<_>>();
    let phase = &ordered.phase;
    if phase.kind != PhaseKind::Done {
        rows.push((
            PhaseRecord {
                index: phase.index,
                kind: phase.kind,
                support: phase.support,
                owner: phase.owner,
                seconds: phase.seconds.unwrap_or_default(),
                keys: phase.keys,
                step_cache_hits: phase.step_cache_hits,
                step_cache_misses: phase.step_cache_misses,
                frontier_keys: ordered.frontier.keys,
                frontier_monomials: ordered.frontier.monomials,
            },
            if phase.seconds.is_some() {
                RowState::Running
            } else {
                RowState::Idle
            },
        ));
    }
    rows
}

fn stages_line(model: &MonitorModel, width: u16, no_color: bool) -> Option<Line<'static>> {
    if model.stages.is_empty() {
        return None;
    }
    let mut value = Text::new();
    for total in &model.stages {
        value
            .sep()
            .push(format!("{} ", total.kind.short_label()), Role::Muted)
            .push(
                format::duration(total.seconds),
                if total.kind == model.stage.kind {
                    Role::Time
                } else {
                    Role::Plain
                },
            );
    }
    Some(labelled("Stages", value, None, width, no_color))
}

fn event_lines(model: &MonitorModel, rows: usize, no_color: bool) -> Vec<Line<'static>> {
    let start = model.events.len().saturating_sub(rows);
    model.events[start..]
        .iter()
        .map(|event| {
            Line::from(Span::styled(
                format!("  {event}"),
                style(Role::Plain, no_color),
            ))
        })
        .collect()
}

fn details_lines(model: &MonitorModel, with_events: bool, no_color: bool) -> Vec<Line<'static>> {
    let details = &model.details;
    let age = |value: Option<f64>| value.map_or_else(|| "—".to_owned(), format::duration);
    let mut rows: Vec<(&str, String)> = Vec::new();
    if let Some(ordered) = &model.ordered {
        rows.push(("native request", ordered.request.to_string()));
        rows.push(("census age", format::duration(ordered.census_age_seconds)));
        rows.push((
            "output peak",
            format!("{} mono", format::count(ordered.output.peak_monomials)),
        ));
        rows.push((
            "output grouping",
            format!(
                "{} merged · {} cancelled",
                format::count(ordered.grouping.output_merged),
                format::count(ordered.grouping.output_cancelled)
            ),
        ));
    }
    if let Some(family) = &details.family {
        rows.push(("family", family.clone()));
    }
    if let Some(key) = &details.current_key {
        rows.push(("current key", key.clone()));
    }
    if let Some([r, a, d]) = &details.current_degrees {
        rows.push(("R / A / D", format!("{r} / {a} / {d}")));
    }
    if let Some(expression) = &details.caller_expression {
        rows.push((
            "caller Atom",
            format!(
                "{} · {} summands · {}",
                format::bytes(expression.bytes as u64),
                expression
                    .summands
                    .map_or_else(|| "—".to_owned(), |summands| summands.to_string()),
                super::model::humanize(&expression.stage)
            ),
        ));
    }
    rows.push(("event age", age(details.event_age_seconds)));
    rows.push(("native age", age(details.native_observation_age_seconds)));
    rows.push(("memory age", age(model.memory.sample_age_seconds)));
    if let Some(total) = model.memory.total_bytes {
        rows.push(("host memory", format::bytes(total)));
    }
    if let Some(guard) = &model.memory.guard {
        rows.push((
            "memory guard",
            format!(
                "{} · floor {} · {}",
                guard.status,
                guard
                    .threshold_bytes
                    .map_or_else(|| "—".to_owned(), format::bytes),
                if guard.covered {
                    "covers this worker"
                } else {
                    "coverage unverified"
                }
            ),
        ));
    }
    if let Some(pid) = details.worker_pid {
        rows.push(("worker PID", pid.to_string()));
    }
    if let Some(dfs) = &model.dfs {
        rows.push((
            "DFS frames",
            format!(
                "{} active · {} pending · {} attempts",
                dfs.active,
                dfs.pending_frames,
                dfs.attempts
                    .map_or_else(|| "—".to_owned(), |v| v.to_string())
            ),
        ));
        rows.push(("since completion", age(dfs.since_completion_seconds)));
        rows.push((
            "rate window",
            format!(
                "{} · {} samples · age {}",
                age(dfs.estimates.window_seconds),
                dfs.estimates.samples,
                age(dfs.estimates.sample_age_seconds)
            ),
        ));
        if let Some(entries) = dfs.cached_entries {
            rows.push(("cached entries", entries.to_string()));
        }
    }
    let mut lines = rows
        .into_iter()
        .map(|(label, value)| {
            Line::from(vec![
                Span::styled(format!("  {label:<17}"), style(Role::Label, no_color)),
                Span::styled(value, style(Role::Plain, no_color)),
            ])
        })
        .collect::<Vec<_>>();
    if let Some(dfs) = &model.dfs
        && !dfs.path.is_empty()
    {
        lines.push(Line::from(Span::styled(
            "  dependency path",
            style(Role::Label, no_color),
        )));
        lines.extend(dfs.path.iter().map(|step| {
            Line::from(Span::styled(
                format!("    {step}"),
                style(Role::Plain, no_color),
            ))
        }));
    }
    if with_events && !model.events.is_empty() {
        lines.push(Line::from(Span::styled(
            "  events",
            style(Role::Label, no_color),
        )));
        lines.extend(event_lines(model, model.events.len(), no_color));
    }
    lines
}

fn footer(ui: UiView, no_color: bool) -> Line<'static> {
    let mut text = Text::new();
    for (key, action) in [
        ("q", "detach"),
        ("space", if ui.frozen { "unfreeze" } else { "freeze" }),
        (
            "d",
            if ui.show_details {
                "summary"
            } else {
                "details"
            },
        ),
        ("^C", "interrupt"),
        ("?", "help"),
    ] {
        text.push(" ", Role::Plain)
            .push(key, Role::Count)
            .push(format!(" {action} "), Role::Muted);
    }
    text.line(no_color)
}

/// The summary block: everything about the run that fits a few lines.
fn summary_lines(
    model: &MonitorModel,
    width: u16,
    with_stages: bool,
    no_color: bool,
) -> Vec<Line<'static>> {
    let mut lines = Vec::new();
    lines.extend(run_line(model, width, no_color));
    lines.push(stage_line(model, width, no_color));
    lines.extend(batch_line(model, width, no_color));
    if let Some(ordered) = &model.ordered {
        lines.extend(ordered_lines(ordered, width, no_color));
    } else if let Some(dfs) = &model.dfs {
        lines.extend(dfs_lines(dfs, width, no_color));
    }
    lines.push(memory_line(model, width, no_color));
    if with_stages && let Some(stages) = stages_line(model, width, no_color) {
        lines.push(Line::default());
        lines.push(stages);
    }
    lines
}

/// Width of the summary column of the wide (two-column) layout.
const WIDE_SUMMARY: u16 = 72;

/// Draw the whole monitor; return the scroll extent of the details view.
pub(super) fn draw(
    frame: &mut Frame<'_>,
    model: &MonitorModel,
    ui: UiView,
    no_color: bool,
) -> usize {
    let area = frame.area();
    if area.width == 0 || area.height == 0 {
        return 0;
    }
    let [top, body, bottom] = Layout::vertical([
        Constraint::Length(1),
        Constraint::Min(0),
        Constraint::Length(1),
    ])
    .areas(area);
    frame.render_widget(Paragraph::new(header(model, ui, area.width, no_color)), top);
    frame.render_widget(Paragraph::new(footer(ui, no_color)), bottom);
    // Wide terminals put the phase log beside the summary.
    let wide = area.width >= 2 * WIDE_SUMMARY && model.ordered.is_some();
    let (summary_area, log_area, summary) = if wide {
        let [left, right] =
            Layout::horizontal([Constraint::Length(WIDE_SUMMARY), Constraint::Min(0)]).areas(body);
        (
            left,
            right,
            summary_lines(model, left.width, true, no_color),
        )
    } else {
        let summary = summary_lines(model, body.width, false, no_color);
        let height = u16::try_from(summary.len())
            .unwrap_or(u16::MAX)
            .min(body.height);
        let [upper, lower] =
            Layout::vertical([Constraint::Length(height), Constraint::Min(0)]).areas(body);
        (upper, lower, summary)
    };
    frame.render_widget(Paragraph::new(summary), summary_area);
    let maximum_scroll = if ui.show_details {
        draw_details(frame, log_area, model, ui.scroll, no_color)
    } else {
        draw_log(frame, log_area, model, !wide, wide, no_color);
        0
    };
    if ui.show_help {
        draw_help(frame, area, no_color);
    }
    maximum_scroll
}

/// Phase log, stage totals and recent events below (or beside) the summary.
fn draw_log(
    frame: &mut Frame<'_>,
    area: Rect,
    model: &MonitorModel,
    with_stages: bool,
    with_details: bool,
    no_color: bool,
) {
    if area.height == 0 || area.width == 0 {
        return;
    }
    let width = area.width;
    let mut tail = Vec::new();
    if with_stages && let Some(stages) = stages_line(model, width, no_color) {
        tail.push(stages);
    }
    let mut lines = Vec::new();
    let height = usize::from(area.height);
    if let Some(ordered) = &model.ordered {
        let rows = phase_rows(ordered);
        // Keep at least the rule, the header and one row for the phases.
        let events = (height.saturating_sub(tail.len() + 3 + rows.len()))
            .min(model.events.len())
            .min(4);
        let room = height.saturating_sub(tail.len() + 2 + if events > 0 { events + 1 } else { 0 });
        lines.push(rule("Phases", width, no_color));
        lines.push(phase_table_header(width, no_color));
        let start = rows.len().saturating_sub(room);
        if rows.is_empty() {
            lines.push(Line::from(Span::styled(
                "  no phase finished yet",
                style(Role::Muted, no_color),
            )));
        }
        lines.extend(
            rows[start..]
                .iter()
                .map(|(record, state)| phase_row(record, *state, width, no_color)),
        );
        if events > 0 {
            lines.push(rule("Events", width, no_color));
            lines.extend(event_lines(model, events, no_color));
        }
    } else if !model.events.is_empty() {
        let events = height.saturating_sub(tail.len() + 1);
        lines.push(rule("Events", width, no_color));
        lines.extend(event_lines(model, events, no_color));
    }
    // Wide terminals use the room below the log for the details.
    let room = height.saturating_sub(lines.len() + tail.len() + 2);
    if with_details && room >= 3 {
        lines.push(Line::default());
        lines.push(rule("Details", width, no_color));
        lines.extend(details_lines(model, false, no_color).into_iter().take(room));
    }
    let free = height.saturating_sub(lines.len() + tail.len());
    lines.extend(std::iter::repeat_n(Line::default(), free));
    lines.extend(tail);
    frame.render_widget(Paragraph::new(lines), area);
}

fn draw_details(
    frame: &mut Frame<'_>,
    area: Rect,
    model: &MonitorModel,
    scroll: usize,
    no_color: bool,
) -> usize {
    if area.height < 2 || area.width == 0 {
        return 0;
    }
    let [title, rest] = Layout::vertical([Constraint::Length(1), Constraint::Min(0)]).areas(area);
    let lines = details_lines(model, true, no_color);
    let paragraph = Paragraph::new(lines).wrap(Wrap { trim: false });
    let maximum = paragraph
        .line_count(rest.width)
        .saturating_sub(usize::from(rest.height))
        // Paragraph's scroll coordinate is u16, a display encoding bound.
        .min(usize::from(u16::MAX));
    let offset = scroll.min(maximum);
    let heading = if maximum > 0 {
        format!("Details ({offset}/{maximum}, ↑↓ scroll)")
    } else {
        "Details".to_owned()
    };
    frame.render_widget(Paragraph::new(rule(&heading, area.width, no_color)), title);
    frame.render_widget(
        paragraph.scroll((u16::try_from(offset).unwrap_or(u16::MAX), 0)),
        rest,
    );
    maximum
}

fn draw_help(frame: &mut Frame<'_>, area: Rect, no_color: bool) {
    let width = area.width.min(74);
    let height = area.height.min(17);
    let popup = Rect::new(
        area.x + area.width.saturating_sub(width) / 2,
        area.y + area.height.saturating_sub(height) / 2,
        width,
        height,
    );
    frame.render_widget(Clear, popup);
    frame.render_widget(
        Paragraph::new(
            "q        detach this display; the reduction continues\n\
             space    freeze the native census; clocks and memory stay live\n\
             d        details: ages, family, current key, guard, events\n\
             ↑↓ PgUp  scroll the details\n\
             Ctrl-C   restore the terminal, then interrupt as usual\n\
             ?        close this help\n\n\
             Phase n/N is a position in a fixed schedule, not a time\n\
             fraction: phase sizes are unknown in advance, so the Ordered\n\
             engine shows rates but no ETA. Keys are Route/Apply identities;\n\
             mono = numerator+denominator monomials of the stored weights.\n\
             Grouping counts global frontier merges and exact cancellations.",
        )
        .block(
            Block::bordered()
                .border_type(BorderType::Rounded)
                .border_style(style(Role::Rule, no_color))
                .title(Span::styled(" Monitor help ", style(Role::Title, no_color))),
        )
        .wrap(Wrap { trim: false }),
        popup,
    );
}

#[cfg(test)]
pub(super) mod tests {
    use std::sync::{
        Arc, Barrier,
        atomic::{AtomicBool, AtomicU64, Ordering},
        mpsc,
    };
    use std::thread::{self, JoinHandle};
    use std::time::Instant;

    use ratatui::{backend::TestBackend, buffer::Buffer};

    use super::super::model::{
        Census, Counter, Details, DfsEstimates, DfsModel, EstimateState, Grouping, GuardSummary,
        MemoryModel, PhaseModel, RunConfig, StageModel, StageTotal,
    };
    use super::*;

    pub(in super::super) fn render(
        model: &MonitorModel,
        ui: UiView,
        width: u16,
        height: u16,
        no_color: bool,
    ) -> Buffer {
        let mut terminal = Terminal::new(TestBackend::new(width, height)).unwrap();
        terminal
            .draw(|frame| {
                draw(frame, model, ui, no_color);
            })
            .unwrap();
        terminal.backend().buffer().clone()
    }

    /// The rendered rows, right-trimmed, one per line.
    pub(in super::super) fn text(buffer: &Buffer) -> String {
        let mut rows = buffer
            .content
            .chunks(usize::from(buffer.area.width))
            .map(|cells| {
                cells
                    .iter()
                    .map(|cell| cell.symbol())
                    .collect::<String>()
                    .trim_end()
                    .to_owned()
            })
            .collect::<Vec<_>>()
            .join("\n");
        rows.push('\n');
        rows
    }

    /// Human-readable golden screens under `monitor/snapshots/`. Set
    /// `VAKINT_UPDATE_TUI_SNAPSHOTS=1` to rewrite them after a layout change.
    fn assert_screen(name: &str, rendered: &str) {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("src/rustred_evaluation/monitor/snapshots")
            .join(format!("{name}.txt"));
        if std::env::var_os("VAKINT_UPDATE_TUI_SNAPSHOTS").is_some() {
            std::fs::create_dir_all(path.parent().unwrap()).unwrap();
            std::fs::write(&path, rendered).unwrap();
            return;
        }
        let expected = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("missing golden screen {}: {error}", path.display()));
        assert_eq!(
            rendered, expected,
            "screen {name} changed; rerun with VAKINT_UPDATE_TUI_SNAPSHOTS=1 and review the diff:\n{rendered}"
        );
    }

    fn running() -> UiView {
        UiView::default()
    }

    fn phase(
        index: usize,
        kind: PhaseKind,
        support: usize,
        owner: Option<usize>,
        seconds: f64,
        keys: u64,
        frontier: (u64, u64),
    ) -> PhaseRecord {
        PhaseRecord {
            index,
            kind,
            support: Some(support),
            owner,
            seconds,
            keys,
            step_cache_hits: 0,
            step_cache_misses: 0,
            frontier_keys: frontier.0,
            frontier_monomials: frontier.1,
        }
    }

    /// The 5L zero-momentum Ordered run of 2026-10-30 about 2.5 h in: Route
    /// support 5 (phase 76 of 83) has been running since t = 993 s. Counters
    /// follow the recorded per-phase timeline (single core, exact weights).
    pub(in super::super) fn five_loop_route_five() -> MonitorModel {
        MonitorModel {
            scope: ScopeKind::Api,
            status: MonitorStatus::Running,
            api_status: ApiStatus::Running,
            elapsed_seconds: 9000.0,
            run: Some(RunConfig {
                n_cores: 1,
                arithmetic: "exact".into(),
                step_cache: "off".into(),
            }),
            stage: StageModel {
                kind: StageKind::Native,
                id: "native".into(),
                activity: Some("coefficient combination".into()),
                seconds: 8847.0,
                work: None,
            },
            stages: vec![
                StageTotal {
                    kind: StageKind::Setup,
                    seconds: 0.4,
                },
                StageTotal {
                    kind: StageKind::DataLoad,
                    seconds: 115.2,
                },
                StageTotal {
                    kind: StageKind::Ingress,
                    seconds: 37.4,
                },
                StageTotal {
                    kind: StageKind::Native,
                    seconds: 8847.0,
                },
            ],
            batch: Some(Counter {
                completed: 0,
                total: 1,
                label: "5L / Retain".into(),
            }),
            cohort: Some(Counter {
                completed: 0,
                total: 1,
                label: String::new(),
            }),
            ordered: Some(OrderedModel {
                request: 2,
                phase: PhaseModel {
                    index: 76,
                    total: 83,
                    kind: PhaseKind::Route,
                    support: Some(5),
                    owner: None,
                    seconds: Some(8007.0),
                    keys: 77_668,
                    step_cache_hits: 0,
                    step_cache_misses: 0,
                },
                input_keys: 1750,
                processed_keys: 381_512,
                recent_rate: Some(9.5),
                request_rate: Some(42.6),
                frontier: Census {
                    keys: 228_412,
                    peak_keys: 228_412,
                    monomials: 29_731_904,
                    peak_monomials: 29_731_904,
                },
                output: Census {
                    keys: 121,
                    peak_keys: 121,
                    monomials: 41_206,
                    peak_monomials: 41_206,
                },
                grouping: Grouping {
                    emitted: 14_102_331,
                    merged: 13_594_647,
                    cancelled: 61_455,
                    output_merged: 9_871,
                    output_cancelled: 12,
                },
                step_cache_hits: 0,
                step_cache_misses: 0,
                phases: vec![
                    phase(
                        64,
                        PhaseKind::Apply,
                        8,
                        Some(12),
                        4.1,
                        3_904,
                        (52_310, 0.9e6 as u64),
                    ),
                    phase(
                        66,
                        PhaseKind::Route,
                        8,
                        None,
                        6.2,
                        9_812,
                        (71_402, 1.4e6 as u64),
                    ),
                    phase(
                        67,
                        PhaseKind::Apply,
                        7,
                        Some(3),
                        3.0,
                        2_210,
                        (70_118, 1.5e6 as u64),
                    ),
                    phase(
                        70,
                        PhaseKind::Route,
                        7,
                        None,
                        18.4,
                        21_544,
                        (139_207, 3.5e6 as u64),
                    ),
                    phase(
                        71,
                        PhaseKind::Route,
                        6,
                        None,
                        350.0,
                        92_072,
                        (216_330, 7.5e6 as u64),
                    ),
                    phase(
                        72,
                        PhaseKind::Apply,
                        6,
                        Some(1),
                        371.0,
                        35_534,
                        (230_041, 16.6e6 as u64),
                    ),
                    phase(
                        73,
                        PhaseKind::Apply,
                        6,
                        Some(6),
                        22.0,
                        21_902,
                        (214_771, 16.0e6 as u64),
                    ),
                    phase(
                        74,
                        PhaseKind::Apply,
                        6,
                        Some(8),
                        24.0,
                        23_455,
                        (203_118, 15.8e6 as u64),
                    ),
                    phase(
                        75,
                        PhaseKind::Apply,
                        6,
                        Some(24),
                        24.0,
                        24_309,
                        (197_688, 15.8e6 as u64),
                    ),
                ],
                census_age_seconds: 0.6,
            }),
            dfs: None,
            memory: MemoryModel {
                verified: true,
                rss_bytes: Some(8_697_000_000),
                peak_rss_bytes: Some(8_701_000_000),
                cpu_percent: Some(99.8),
                available_bytes: Some(222_500_000_000),
                total_bytes: Some(1_081_000_000_000),
                sample_age_seconds: Some(0.4),
                guard: None,
            },
            details: Details {
                event_age_seconds: Some(0.3),
                native_observation_age_seconds: Some(0.6),
                family: Some("5L/zero-momentum/selection-v3:2f0c… (5104 chars)".into()),
                current_key: Some("Ordered [1,1,1,1,1,0,0,0,0,0,0,0,0,0,0]".into()),
                current_degrees: Some(["0".into(), "5".into(), "5".into()]),
                caller_expression: None,
                worker_pid: Some(4_123_599),
            },
            events: Vec::new(),
        }
    }

    /// A finished three-loop evaluation: four lower-loop Reduce phases.
    fn three_loop_finished() -> MonitorModel {
        MonitorModel {
            scope: ScopeKind::Api,
            status: MonitorStatus::Succeeded,
            api_status: ApiStatus::Succeeded,
            elapsed_seconds: 3.4,
            run: Some(RunConfig {
                n_cores: 10,
                arithmetic: "modular_polynomial".into(),
                step_cache: "read-write".into(),
            }),
            stage: StageModel {
                kind: StageKind::PostProcessing,
                id: "evaluation_output".into(),
                activity: None,
                seconds: 0.2,
                work: None,
            },
            stages: vec![
                StageTotal {
                    kind: StageKind::Setup,
                    seconds: 0.1,
                },
                StageTotal {
                    kind: StageKind::Ingress,
                    seconds: 0.9,
                },
                StageTotal {
                    kind: StageKind::Native,
                    seconds: 2.0,
                },
                StageTotal {
                    kind: StageKind::PostProcessing,
                    seconds: 0.4,
                },
            ],
            batch: Some(Counter {
                completed: 1,
                total: 1,
                label: "3L / Retain".into(),
            }),
            cohort: Some(Counter {
                completed: 2,
                total: 2,
                label: String::new(),
            }),
            ordered: Some(OrderedModel {
                request: 3,
                phase: PhaseModel {
                    index: 7,
                    total: 7,
                    kind: PhaseKind::Done,
                    ..Default::default()
                },
                input_keys: 2,
                processed_keys: 1_876,
                recent_rate: None,
                request_rate: Some(1_403.0),
                frontier: Census {
                    keys: 0,
                    peak_keys: 611,
                    monomials: 0,
                    peak_monomials: 48_220,
                },
                output: Census {
                    keys: 3,
                    peak_keys: 3,
                    monomials: 1_207,
                    peak_monomials: 1_207,
                },
                grouping: Grouping {
                    emitted: 21_406,
                    merged: 19_310,
                    cancelled: 214,
                    output_merged: 117,
                    output_cancelled: 0,
                },
                step_cache_hits: 1_701,
                step_cache_misses: 175,
                phases: vec![
                    PhaseRecord {
                        step_cache_hits: 812,
                        step_cache_misses: 14,
                        ..phase(1, PhaseKind::Reduce, 6, None, 0.6, 826, (402, 31_000))
                    },
                    PhaseRecord {
                        step_cache_hits: 610,
                        step_cache_misses: 92,
                        ..phase(2, PhaseKind::Reduce, 5, None, 0.7, 702, (118, 6_904))
                    },
                    PhaseRecord {
                        step_cache_hits: 279,
                        step_cache_misses: 69,
                        ..phase(3, PhaseKind::Reduce, 4, None, 0.3, 348, (0, 0))
                    },
                ],
                census_age_seconds: 0.4,
            }),
            dfs: None,
            memory: MemoryModel {
                verified: true,
                rss_bytes: Some(1_291_000_000),
                peak_rss_bytes: Some(1_420_000_000),
                cpu_percent: Some(310.0),
                available_bytes: Some(231_900_000_000),
                total_bytes: Some(1_081_000_000_000),
                sample_age_seconds: Some(0.2),
                guard: None,
            },
            details: Details::default(),
            events: vec!["request finished · caller · succeeded".into()],
        }
    }

    fn depth_first() -> MonitorModel {
        MonitorModel {
            status: MonitorStatus::Running,
            api_status: ApiStatus::Running,
            elapsed_seconds: 812.0,
            stage: StageModel {
                kind: StageKind::Native,
                id: "native".into(),
                activity: Some("transport".into()),
                seconds: 790.0,
                work: None,
            },
            stages: vec![StageTotal {
                kind: StageKind::Native,
                seconds: 790.0,
            }],
            dfs: Some(DfsModel {
                roots_completed: Some(412),
                roots_total: Some(1750),
                entered: 20_412,
                completed: Some(18_301),
                attempts: Some(18_950),
                active: 31,
                pending_frames: 2_204,
                pending_unique: Some(1_190),
                pending_exact: false,
                persistent_hits: Some(1_204),
                persistent_misses: Some(311),
                ram_hits: Some(5_633),
                cached_entries: Some(3_201),
                since_completion_seconds: Some(0.8),
                estimates: DfsEstimates {
                    state: EstimateState::Live,
                    window_seconds: Some(60.0),
                    samples: 61,
                    sample_age_seconds: Some(0.8),
                    completion_rate: Some(12.3),
                    root_rate: Some(0.52),
                    growth_rate: Some(0.5),
                    frontier_eta_seconds: Some(96.7),
                    root_eta_seconds: Some(2_573.0),
                },
                collection: None,
                path: vec!["Route [2,1,1]".into(), "Apply owner [1,1,0] [2,1,0]".into()],
            }),
            memory: MemoryModel {
                verified: true,
                rss_bytes: Some(4_400_000_000),
                peak_rss_bytes: Some(4_500_000_000),
                cpu_percent: Some(100.0),
                available_bytes: Some(212_000_000_000),
                total_bytes: Some(1_081_000_000_000),
                sample_age_seconds: Some(0.5),
                guard: Some(GuardSummary {
                    status: "active".into(),
                    threshold_bytes: Some(200_000_000_000),
                    covered: true,
                }),
            },
            ..Default::default()
        }
    }

    #[test]
    fn five_loop_route_five_at_80x24() {
        let rendered = text(&render(&five_loop_route_five(), running(), 80, 24, true));
        assert_screen("ordered_5l_route5_80x24", &rendered);
        // The essentials are on the main screen, without disclaimers.
        for expected in [
            "76/83 Route s5",
            "2:13:27 in phase",
            "9.5/s now",
            "228.4k keys",
            "96.4% merged",
            "1/1 5L / Retain",
            "RSS 8.1 GiB",
        ] {
            assert!(
                rendered.contains(expected),
                "missing {expected}:\n{rendered}"
            );
        }
        assert!(!rendered.contains("not time fraction") && !rendered.contains("chars)"));
    }

    #[test]
    fn five_loop_route_five_wide() {
        let rendered = text(&render(&five_loop_route_five(), running(), 160, 40, true));
        assert_screen("ordered_5l_route5_160x40", &rendered);
        assert!(rendered.contains("Stages"));
        assert!(rendered.contains("running"));
    }

    #[test]
    fn three_loop_finished_at_80x24() {
        let rendered = text(&render(&three_loop_finished(), running(), 80, 24, true));
        assert_screen("ordered_3l_finished_80x24", &rendered);
        assert!(rendered.contains("SUCCEEDED") && rendered.contains("7/7 done"));
        assert!(rendered.contains("90.7% replayed"));
    }

    #[test]
    fn depth_first_at_80x24() {
        let rendered = text(&render(&depth_first(), running(), 80, 24, true));
        assert_screen("dfs_80x24", &rendered);
        assert!(rendered.contains("412/1750") && rendered.contains("ETA 42:53"));
    }

    #[test]
    fn details_view_scrolls_and_clamps_at_its_content() {
        let mut model = five_loop_route_five();
        model.events = (0..30).map(|index| format!("event {index}")).collect();
        let ui = UiView {
            show_details: true,
            scroll: usize::MAX,
            ..Default::default()
        };
        let mut terminal = Terminal::new(TestBackend::new(80, 24)).unwrap();
        let mut maximum = 0;
        terminal
            .draw(|frame| maximum = draw(frame, &model, ui, true))
            .unwrap();
        let rendered = text(terminal.backend().buffer());
        assert!(maximum > 0);
        assert!(rendered.contains("event 29"), "{rendered}");
        assert!(rendered.contains(&format!("Details ({maximum}/{maximum}")));
        let top = text(&render(
            &model,
            UiView {
                show_details: true,
                ..Default::default()
            },
            100,
            40,
            true,
        ));
        for expected in ["native request", "family", "R / A / D", "worker PID"] {
            assert!(top.contains(expected), "missing {expected}:\n{top}");
        }
        assert!(top.contains("(5104 chars)"));
    }

    #[test]
    fn help_and_tiny_viewports_remain_display_only() {
        let ui = UiView {
            show_help: true,
            ..Default::default()
        };
        let rendered = text(&render(&five_loop_route_five(), ui, 100, 30, true));
        assert!(rendered.contains("reduction continues"));
        assert!(rendered.contains("no ETA"));
        for (width, height) in [(1, 1), (3, 2), (20, 5), (79, 8)] {
            let buffer = render(&five_loop_route_five(), ui, width, height, true);
            assert_eq!(
                buffer.content.len(),
                usize::from(width) * usize::from(height)
            );
        }
    }

    #[test]
    fn no_color_resets_every_cell_and_keeps_the_same_text() {
        for model in [five_loop_route_five(), three_loop_finished(), depth_first()] {
            for (width, height) in [(80, 24), (160, 40)] {
                let plain = render(&model, running(), width, height, true);
                assert!(
                    plain
                        .content
                        .iter()
                        .all(|cell| cell.fg == Color::Reset && cell.bg == Color::Reset)
                );
                let colored = render(&model, running(), width, height, false);
                assert_eq!(text(&plain), text(&colored));
                assert!(colored.content.iter().any(|cell| cell.fg != Color::Reset));
            }
        }
    }

    #[test]
    fn values_are_styled_by_their_role_not_by_their_label() {
        let buffer = render(&five_loop_route_five(), running(), 80, 24, false);
        let find = |value: &str| {
            let characters = value.chars().map(|ch| ch.to_string()).collect::<Vec<_>>();
            buffer
                .content
                .windows(characters.len())
                .find(|cells| {
                    cells
                        .iter()
                        .zip(&characters)
                        .all(|(cell, ch)| cell.symbol() == ch)
                })
                .unwrap_or_else(|| panic!("{value} not rendered"))
                .to_vec()
        };
        for (value, color) in [
            ("2:13:27", CYAN),
            ("9.5/s", GREEN),
            ("96.4%", GREEN),
            ("8.1 GiB", BLUE),
            ("Route s5", MAGENTA),
        ] {
            assert!(
                find(value).iter().all(|cell| cell.fg == color),
                "{value} is not {color:?}"
            );
        }
        let mut failed = five_loop_route_five();
        failed.status = MonitorStatus::Failed;
        let buffer = render(&failed, running(), 80, 24, false);
        let header = &buffer.content[..80];
        assert!(header.iter().any(|cell| cell.fg == RED));
    }

    #[test]
    fn frozen_display_is_announced_and_clocks_stay_in_the_header() {
        let ui = UiView {
            frozen: true,
            ..Default::default()
        };
        let rendered = text(&render(&five_loop_route_five(), ui, 80, 24, true));
        let header = rendered.lines().next().unwrap();
        assert!(header.contains("display frozen") && header.contains("2:30:00"));
        assert!(rendered.contains("unfreeze"));
    }

    #[test]
    fn terminal_claim_releases_after_error_and_unwind() {
        let failed = (|| -> io::Result<()> {
            let _claim = TerminalClaim::acquire()?;
            Err(io::Error::other("intentional UI admission failure"))
        })();
        assert!(failed.is_err());
        drop(TerminalClaim::acquire().unwrap());
        let unwound = std::panic::catch_unwind(|| {
            let _claim = TerminalClaim::acquire().unwrap();
            panic!("intentional UI claim unwind");
        });
        assert!(unwound.is_err());
        let claim = TerminalClaim::acquire().unwrap();
        assert!(
            matches!(TerminalClaim::acquire(), Err(error) if error.kind() == io::ErrorKind::AlreadyExists)
        );
        drop(claim);
        drop(TerminalClaim::acquire().unwrap());
    }

    fn pty_message(message: &str) -> io::Result<()> {
        let mut output = io::stdout().lock();
        writeln!(output, "{message}")?;
        output.flush()
    }

    fn pty_action() -> UiAction {
        let deadline = Instant::now() + Duration::from_secs(5);
        loop {
            let remaining = deadline.saturating_duration_since(Instant::now());
            assert!(
                !remaining.is_zero(),
                "PTY controller did not send an action"
            );
            if let Some(action) = poll_action(remaining).unwrap() {
                return action;
            }
        }
    }

    /// A caller-owned terminal lets the controller check nested refusal before cleanup.
    struct OuterTerminal;

    impl OuterTerminal {
        fn enter() -> io::Result<Self> {
            let guard = Self;
            enable_raw_mode()?;
            execute!(io::stderr(), EnterAlternateScreen, cursor::Hide)?;
            Ok(guard)
        }
    }

    impl Drop for OuterTerminal {
        fn drop(&mut self) {
            let _ = execute!(io::stderr(), cursor::Show, LeaveAlternateScreen);
            let _ = disable_raw_mode();
            let _ = io::stderr().flush();
        }
    }

    /// This worker counts ordinary test work; it has no solver or cancellation bridge.
    struct FakeReduction {
        stop_worker: Arc<AtomicBool>,
        counter: Arc<AtomicU64>,
        worker: Option<JoinHandle<()>>,
    }

    impl FakeReduction {
        fn start() -> Self {
            let stop_worker = Arc::new(AtomicBool::new(false));
            let counter = Arc::new(AtomicU64::new(0));
            let worker_stop = Arc::clone(&stop_worker);
            let worker_counter = Arc::clone(&counter);
            let worker = thread::spawn(move || {
                while !worker_stop.load(Ordering::Relaxed) {
                    worker_counter.fetch_add(1, Ordering::Relaxed);
                    thread::sleep(Duration::from_millis(1));
                }
            });
            Self {
                stop_worker,
                counter,
                worker: Some(worker),
            }
        }

        fn assert_advancing(&self) {
            let before = self.counter.load(Ordering::Relaxed);
            let deadline = Instant::now() + Duration::from_secs(1);
            thread::sleep(Duration::from_millis(20));
            while self.counter.load(Ordering::Relaxed) <= before {
                assert!(Instant::now() < deadline, "fake worker stopped advancing");
                thread::sleep(Duration::from_millis(1));
            }
        }
    }

    impl Drop for FakeReduction {
        fn drop(&mut self) {
            self.stop_worker.store(true, Ordering::Relaxed);
            if let Some(worker) = self.worker.take() {
                let _ = worker.join();
            }
        }
    }

    fn pty_error_case() -> io::Result<()> {
        let mut terminal = TerminalGuard::enter()?;
        terminal.draw(&mut MonitorState::default(), Duration::ZERO, true)?;
        pty_message("PTY_READY:error")?;
        Err(io::Error::other("intentional PTY component error"))
    }

    fn pty_key(expected: char) {
        let deadline = Instant::now() + Duration::from_secs(5);
        loop {
            let remaining = deadline.saturating_duration_since(Instant::now());
            assert!(
                !remaining.is_zero(),
                "PTY controller did not send its synchronization key"
            );
            assert!(
                event::poll(remaining).unwrap(),
                "PTY controller wait expired"
            );
            if let Event::Key(key) = event::read().unwrap()
                && key.kind == KeyEventKind::Press
            {
                assert_eq!(key.code, KeyCode::Char(expected));
                break;
            }
        }
    }

    fn pty_nested_case() {
        let outer = OuterTerminal::enter().unwrap();
        pty_message("PTY_OUTER_READY:nested").unwrap();
        pty_key('i');
        let error = match TerminalGuard::enter() {
            Ok(_) => panic!("nested terminal entry unexpectedly succeeded"),
            Err(error) => error,
        };
        assert_eq!(error.kind(), io::ErrorKind::AlreadyExists);
        assert!(is_raw_mode_enabled().unwrap());
        pty_message("PTY_READY:nested").unwrap();
        pty_key('c');
        assert!(is_raw_mode_enabled().unwrap());
        drop(outer);
    }

    fn pty_concurrent_case() {
        let start = Arc::new(Barrier::new(3));
        let (outcomes_tx, outcomes_rx) = mpsc::channel();
        let mut releases = Vec::new();
        let mut workers = Vec::new();
        for index in 0..2 {
            let start = Arc::clone(&start);
            let outcomes = outcomes_tx.clone();
            let (release_tx, release_rx) = mpsc::channel();
            releases.push(release_tx);
            workers.push(thread::spawn(move || {
                start.wait();
                match TerminalGuard::enter() {
                    Ok(guard) => {
                        outcomes.send((index, Ok(()))).unwrap();
                        // Keep the winner alive through both observed outcomes.
                        release_rx.recv_timeout(Duration::from_secs(8)).unwrap();
                        drop(guard);
                    }
                    Err(error) => outcomes.send((index, Err(error.kind()))).unwrap(),
                }
            }));
        }
        drop(outcomes_tx);
        start.wait();
        let outcomes = [
            outcomes_rx.recv_timeout(Duration::from_secs(5)).unwrap(),
            outcomes_rx.recv_timeout(Duration::from_secs(5)).unwrap(),
        ];
        let winners = outcomes
            .iter()
            .filter_map(|(index, result)| result.is_ok().then_some(*index))
            .collect::<Vec<_>>();
        let failures = outcomes
            .iter()
            .filter_map(|(_, result)| result.as_ref().err().copied())
            .collect::<Vec<_>>();
        assert_eq!(winners.len(), 1);
        assert_eq!(failures, vec![io::ErrorKind::AlreadyExists]);
        assert!(is_raw_mode_enabled().unwrap());
        pty_message("PTY_READY:concurrent").unwrap();
        pty_key('c');
        releases[winners[0]].send(()).unwrap();
        for worker in workers {
            worker.join().unwrap();
        }
        assert!(!is_raw_mode_enabled().unwrap());
        // This succeeds only after the winner and its backend released the claim.
        drop(TerminalGuard::enter().unwrap());
        assert!(!is_raw_mode_enabled().unwrap());
        pty_message("PTY_REENTERED:concurrent").unwrap();
    }

    fn pty_controls_case(case: &str) {
        let worker = FakeReduction::start();
        let mut terminal = Some(TerminalGuard::enter().unwrap());
        let mut state = MonitorState::default();
        let started = Instant::now();
        terminal
            .as_mut()
            .unwrap()
            .draw(&mut state, started.elapsed(), true)
            .unwrap();
        pty_message(&format!("PTY_READY:{case}")).unwrap();
        let controls = [
            UiAction::ToggleFreeze,
            UiAction::ToggleHelp,
            UiAction::ScrollDown,
            UiAction::PageDown,
            UiAction::PageUp,
            UiAction::ScrollUp,
            UiAction::ToggleHelp,
            UiAction::ToggleFreeze,
            UiAction::Detach,
        ];
        let expected = if case == "q" {
            &controls[8..]
        } else {
            &controls[..]
        };
        for (index, &expected) in expected.iter().enumerate() {
            let action = pty_action();
            assert_eq!(action, expected);
            state.apply(action, started.elapsed());
            if case == "controls" {
                let view = state.ui();
                match index {
                    0 | 7 => assert_eq!(view.frozen, index == 0),
                    1 | 6 => assert_eq!(view.show_help, index == 1),
                    // The summary view has no scroll extent: each frame
                    // clamps the offset that the next key starts from.
                    2..=5 => assert_eq!(view.scroll, [1, 10, 0, 0][index - 2]),
                    _ => {}
                }
            }
            if action == UiAction::Detach {
                assert!(state.detached());
                drop(terminal.take());
                assert!(!is_raw_mode_enabled().unwrap());
            } else {
                terminal
                    .as_mut()
                    .unwrap()
                    .draw(&mut state, started.elapsed(), true)
                    .unwrap();
            }
            if matches!(action, UiAction::ToggleFreeze | UiAction::Detach) {
                worker.assert_advancing();
            }
            if action == UiAction::Detach {
                pty_message(&format!("PTY_COUNTER_CONTINUED:{case}")).unwrap();
            }
            pty_message(&format!("PTY_ACK:{case}:{index}")).unwrap();
        }
        assert!(state.detached());
        assert!(terminal.is_none());
    }

    #[test]
    #[ignore = "requires a fresh controlling PTY and its external component controller"]
    fn real_pty_controller_case() {
        let case = std::env::var("VAKINT_MONITOR_PTY_SCENARIO")
            .expect("set VAKINT_MONITOR_PTY_SCENARIO for the PTY component controller");
        assert!(!is_raw_mode_enabled().unwrap());
        match case.as_str() {
            "normal" => {
                let mut terminal = TerminalGuard::enter().unwrap();
                terminal
                    .draw(&mut MonitorState::default(), Duration::ZERO, true)
                    .unwrap();
                pty_message("PTY_READY:normal").unwrap();
                drop(terminal);
            }
            "error" => {
                let error = pty_error_case().unwrap_err();
                assert_eq!(error.kind(), io::ErrorKind::Other);
                assert_eq!(error.to_string(), "intentional PTY component error");
            }
            "unwind" => {
                let result = std::panic::catch_unwind(|| {
                    let mut terminal = TerminalGuard::enter().unwrap();
                    terminal
                        .draw(&mut MonitorState::default(), Duration::ZERO, true)
                        .unwrap();
                    pty_message("PTY_READY:unwind").unwrap();
                    panic!("intentional PTY component unwind");
                });
                assert!(result.is_err());
                drop(TerminalClaim::acquire().unwrap());
            }
            "nested" => pty_nested_case(),
            "concurrent" => pty_concurrent_case(),
            "external" => {
                let error = match TerminalGuard::enter() {
                    Ok(_) => panic!("externally raw terminal entry unexpectedly succeeded"),
                    Err(error) => error,
                };
                assert_eq!(error.kind(), io::ErrorKind::AlreadyExists);
                drop(TerminalClaim::acquire().unwrap());
                pty_message("PTY_READY:external").unwrap();
            }
            "q" | "controls" => pty_controls_case(&case),
            _ => panic!("unknown PTY component scenario"),
        }
        assert!(!is_raw_mode_enabled().unwrap());
        pty_message(&format!("PTY_DONE:{case}")).unwrap();
    }
}
