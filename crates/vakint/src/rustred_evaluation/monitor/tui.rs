//! Terminal presentation of bounded, caller-owned monitor views.

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

use super::{MonitorState, MonitorStatus, MonitorView, UiAction};

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
        state: &MonitorState,
        now: Duration,
        no_color: bool,
    ) -> io::Result<()> {
        let view = state.view(now);
        self.draw_view(&view, no_color)
    }

    pub(super) fn draw_view(&mut self, view: &MonitorView, no_color: bool) -> io::Result<()> {
        self.terminal
            .draw(|frame| draw(frame, view, no_color))
            .map(|_| ())
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

/// Poll once on the terminal-owner thread; no input loop runs in the observer.
pub(super) fn poll_action(timeout: Duration) -> io::Result<Option<UiAction>> {
    if !event::poll(timeout)? {
        return Ok(None);
    }
    let Event::Key(key) = event::read()? else {
        return Ok(None);
    };
    if key.kind == KeyEventKind::Release
        || key
            .modifiers
            .intersects(KeyModifiers::CONTROL | KeyModifiers::ALT | KeyModifiers::SUPER)
    {
        return Ok(None);
    }
    // Raw mode captures control keys. Only these agreed display actions are
    // mapped here; cancellation and signal policy belong to the caller.
    let action = match key.code {
        KeyCode::Char('q') if key.kind == KeyEventKind::Press => UiAction::Detach,
        KeyCode::Char(' ') if key.kind == KeyEventKind::Press => UiAction::ToggleFreeze,
        KeyCode::Char('?') if key.kind == KeyEventKind::Press => UiAction::ToggleHelp,
        KeyCode::Up => UiAction::ScrollUp,
        KeyCode::Down => UiAction::ScrollDown,
        KeyCode::PageUp => UiAction::PageUp,
        KeyCode::PageDown => UiAction::PageDown,
        _ => return Ok(None),
    };
    Ok(Some(action))
}

/// Draw only supplied labels and values; estimates and metric scopes are model-owned.
pub(super) fn draw(frame: &mut Frame<'_>, view: &MonitorView, no_color: bool) {
    let area = frame.area();
    if area.width == 0 || area.height == 0 {
        return;
    }
    let [header, body, footer] = Layout::vertical([
        Constraint::Length(4),
        Constraint::Min(0),
        Constraint::Length(1),
    ])
    .areas(area);
    let (status, color) = match view.status {
        MonitorStatus::Running => ("Running", Color::Cyan),
        MonitorStatus::Succeeded => ("Succeeded", Color::Green),
        MonitorStatus::Failed => ("Failed", Color::Red),
        MonitorStatus::Partial => ("Returned partial", Color::Yellow),
        MonitorStatus::Detached => ("Detached", Color::Yellow),
    };
    let status_style = if no_color {
        Style::default().add_modifier(Modifier::BOLD)
    } else {
        Style::default().fg(color).add_modifier(Modifier::BOLD)
    };
    let mut status_line = vec![Span::styled(status, status_style)];
    if view.frozen {
        status_line.push(Span::raw("  |  display frozen; reduction continues"));
    }
    frame.render_widget(
        Paragraph::new(vec![
            Line::from(status_line),
            Line::from(view.subtitle.as_str()),
        ])
        .block(panel(view.title.clone(), no_color))
        .wrap(Wrap { trim: false }),
        header,
    );

    if body.width >= 88 && body.height >= 17 {
        draw_panels(frame, body, view, no_color);
    } else {
        draw_compact(frame, body, view, no_color);
    }
    frame.render_widget(
        Paragraph::new("q detach | space freeze | ↑↓ / PgUp PgDn scroll | ? help")
            .style(muted(no_color)),
        footer,
    );
    if view.show_help {
        draw_help(frame, area, no_color);
    }
}

fn panel(title: String, no_color: bool) -> Block<'static> {
    Block::bordered()
        .border_type(BorderType::Rounded)
        .border_style(muted(no_color))
        .title(Line::from(title).style(Style::default().add_modifier(Modifier::BOLD)))
}

fn muted(no_color: bool) -> Style {
    if no_color {
        Style::default()
    } else {
        Style::default().fg(Color::DarkGray)
    }
}

fn rows<'a>(values: &'a [(String, String)], no_color: bool) -> Vec<Line<'a>> {
    values
        .iter()
        .map(|(label, value)| {
            Line::from(vec![
                Span::styled(format!("{label}: "), muted(no_color)),
                Span::raw(value.as_str()),
            ])
        })
        .collect()
}

fn draw_rows(
    frame: &mut Frame<'_>,
    area: Rect,
    title: &str,
    values: &[(String, String)],
    scroll: usize,
    no_color: bool,
) {
    let title = if values.len() > usize::from(area.height.saturating_sub(2)) {
        format!("{title} · offset {scroll} · clipped")
    } else if scroll > 0 {
        format!("{title} · offset {scroll}")
    } else {
        title.to_owned()
    };
    let content = if values.is_empty() {
        vec![Line::from("No observation supplied").style(muted(no_color))]
    } else {
        rows(values, no_color)
    };
    frame.render_widget(
        Paragraph::new(content)
            .block(panel(title, no_color))
            .wrap(Wrap { trim: false })
            .scroll((u16::try_from(scroll).unwrap_or(u16::MAX), 0)),
        area,
    );
}

fn draw_panels(frame: &mut Frame<'_>, area: Rect, view: &MonitorView, no_color: bool) {
    let [left, right] =
        Layout::horizontal([Constraint::Percentage(52), Constraint::Percentage(48)]).areas(area);
    let [overview, frontier, estimates] = Layout::vertical([
        Constraint::Ratio(1, 3),
        Constraint::Ratio(1, 3),
        Constraint::Ratio(1, 3),
    ])
    .areas(left);
    let [cache, complexity, details] = Layout::vertical([
        Constraint::Ratio(1, 3),
        Constraint::Ratio(1, 3),
        Constraint::Ratio(1, 3),
    ])
    .areas(right);
    for (area, title, values) in [
        (overview, "Observation", &view.overview),
        (frontier, "Frontier", &view.frontier),
        (estimates, "Estimates · supplied model", &view.estimates),
        (cache, "Cache / memory · scoped", &view.cache),
        (complexity, "Complexity · scoped", &view.complexity),
    ] {
        draw_rows(frame, area, title, values, view.scroll, no_color);
    }
    let mut lines = vec![Line::from("Dependency path").style(muted(no_color))];
    lines.extend(view.path.iter().map(|line| Line::from(line.as_str())));
    lines.push(Line::from("Events").style(muted(no_color)));
    lines.extend(view.events.iter().map(|line| Line::from(line.as_str())));
    let title = format!("Path / events · line offset {}", view.scroll);
    frame.render_widget(
        Paragraph::new(lines)
            .block(panel(title, no_color))
            .wrap(Wrap { trim: false })
            .scroll((u16::try_from(view.scroll).unwrap_or(u16::MAX), 0)),
        details,
    );
}

fn draw_compact(frame: &mut Frame<'_>, area: Rect, view: &MonitorView, no_color: bool) {
    let mut lines = Vec::new();
    for (title, values) in [
        ("Observation", &view.overview),
        ("Frontier", &view.frontier),
        ("Cache / memory · scoped", &view.cache),
        ("Complexity · scoped", &view.complexity),
        ("Estimates · supplied model", &view.estimates),
    ] {
        lines.push(Line::from(title).style(Style::default().add_modifier(Modifier::BOLD)));
        lines.extend(rows(values, no_color));
        lines.push(Line::default());
    }
    lines.push(Line::from("Dependency path").style(muted(no_color)));
    lines.extend(view.path.iter().map(|line| Line::from(line.as_str())));
    lines.push(Line::from("Events").style(muted(no_color)));
    lines.extend(view.events.iter().map(|line| Line::from(line.as_str())));
    frame.render_widget(
        Paragraph::new(lines)
            .block(panel(
                format!("Observations · line offset {}", view.scroll),
                no_color,
            ))
            .wrap(Wrap { trim: false })
            .scroll((u16::try_from(view.scroll).unwrap_or(u16::MAX), 0)),
        area,
    );
}

fn draw_help(frame: &mut Frame<'_>, area: Rect, no_color: bool) {
    let width = area.width.min(72);
    let height = area.height.min(12);
    let popup = Rect::new(
        area.x + area.width.saturating_sub(width) / 2,
        area.y + area.height.saturating_sub(height) / 2,
        width,
        height,
    );
    frame.render_widget(Clear, popup);
    frame.render_widget(
        Paragraph::new(
            "q              Detach this display; reduction continues\n\
             space          Freeze / resume only the display\n\
             ↑ / ↓          Scroll all observation panels\n\
             PgUp / PgDn    Scroll by a page\n\
             ?              Close this help\n\n\
             Age measures time since a native observation. Long native\n\
             operations need not emit a heartbeat. Estimates are supplied\n\
             by the model, not a percentage of an unknown dependency DAG.",
        )
        .block(panel("Monitor help".into(), no_color))
        .wrap(Wrap { trim: false }),
        popup,
    );
}

#[cfg(test)]
mod tests {
    use std::sync::{
        Arc, Barrier,
        atomic::{AtomicBool, AtomicU64, Ordering},
        mpsc,
    };
    use std::thread::{self, JoinHandle};
    use std::time::Instant;

    use ratatui::{backend::TestBackend, buffer::Buffer};

    use super::*;

    fn view() -> MonitorView {
        MonitorView {
            title: "RustRed observation".into(),
            subtitle: "Safe-boundary sample; future work unknown".into(),
            status: MonitorStatus::Running,
            frozen: false,
            show_help: false,
            scroll: 0,
            overview: vec![("Snapshot age".into(), "3s".into())],
            frontier: vec![("Pending unique".into(), "unknown".into())],
            cache: vec![("Coefficient bytes".into(), "12 (cached scope)".into())],
            complexity: vec![("Degree".into(), "7 (active node, sampled)".into())],
            estimates: vec![
                ("Method".into(), "selected-root throughput".into()),
                ("Estimate scope".into(), "local observed work".into()),
                ("Confidence".into(), "heuristic; age 3s".into()),
            ],
            path: vec!["Route root -> Apply descendant".into()],
            events: vec!["backedge observed; no cycle proof".into()],
        }
    }

    fn render(view: &MonitorView, width: u16, height: u16, no_color: bool) -> Buffer {
        let mut terminal = Terminal::new(TestBackend::new(width, height)).unwrap();
        terminal.draw(|frame| draw(frame, view, no_color)).unwrap();
        terminal.backend().buffer().clone()
    }

    fn text(buffer: &Buffer) -> String {
        buffer
            .content
            .chunks(usize::from(buffer.area.width))
            .map(|cells| cells.iter().map(|cell| cell.symbol()).collect::<String>())
            .collect::<Vec<_>>()
            .join("\n")
    }

    #[test]
    fn memory_labels_keep_os_and_symbolic_scopes_explicit() {
        // Synthetic OS values exercise presentation only, not a reduction.
        let mut state = MonitorState::default();
        state.latest.resources = super::super::ResourceUsage::from_proc(
            "VmRSS: 17 kB\nVmHWM: 29 kB",
            "MemTotal: 1000 kB\nMemAvailable: 700 kB\nMemFree: 200 kB",
        );
        state.latest.resource_sample_at = Some(Duration::ZERO);
        let view = state.view(Duration::from_secs(3));
        let rendered = text(&render(&view, 220, 80, false));
        for label in [
            "VmRSS / VmHWM",
            "MemAvailable / MemFree",
            "usable / unused",
            "MemTotal",
            "sample age",
            "Full symbolic expression size",
            "unknown",
            "Current canonical R / A / D",
        ] {
            assert!(rendered.contains(label), "missing {label}: {rendered}");
        }
        assert!(rendered.contains("Local observations"));
    }

    #[test]
    fn wide_and_narrow_views_preserve_supplied_scopes_and_unknowns() {
        let view = view();
        for (width, height) in [(150, 40), (54, 40)] {
            let rendered = text(&render(&view, width, height, false));
            for supplied in [
                "unknown",
                "cached scope",
                "active node, sampled",
                "selected-root throughput",
                "local observed work",
                "heuristic; age 3s",
                "no cycle proof",
            ] {
                assert!(
                    rendered.contains(supplied),
                    "missing {supplied}: {rendered}"
                );
            }
            assert!(!rendered.contains('%'));
        }
    }

    #[test]
    fn no_color_resets_every_cell_foreground_and_background() {
        let view = view();
        for (width, height) in [(150, 40), (54, 40)] {
            let buffer = render(&view, width, height, true);
            assert!(
                buffer
                    .content
                    .iter()
                    .all(|cell| cell.fg == Color::Reset && cell.bg == Color::Reset)
            );
        }
    }

    #[test]
    fn scrolling_reveals_clipped_panel_rows_without_changing_the_view() {
        let mut view = view();
        view.frontier = (0..20)
            .map(|index| (format!("Row {index}"), "scoped observation".into()))
            .collect();
        assert!(!text(&render(&view, 150, 24, true)).contains("Row 19:"));
        view.scroll = 18;
        let rendered = text(&render(&view, 150, 24, true));
        assert!(rendered.contains("Row 19:"));
        assert!(rendered.contains("offset 18"));
        assert_eq!(view.frontier.len(), 20);
    }

    #[test]
    fn help_and_tiny_viewports_remain_display_only() {
        let mut view = view();
        view.show_help = true;
        let rendered = text(&render(&view, 100, 30, true));
        assert!(rendered.contains("reduction continues"));
        assert!(rendered.contains("not a percentage"));
        for (width, height) in [(1, 1), (3, 2)] {
            assert_eq!(
                render(&view, width, height, true).content.len(),
                usize::from(width) * usize::from(height)
            );
        }
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
        terminal.draw(&MonitorState::default(), Duration::ZERO, true)?;
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
            .draw(&state, started.elapsed(), true)
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
                let view = state.view(started.elapsed());
                match index {
                    0 | 7 => assert_eq!(view.frozen, index == 0),
                    1 | 6 => assert_eq!(view.show_help, index == 1),
                    2..=5 => assert_eq!(view.scroll, [1, 11, 1, 0][index - 2]),
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
                    .draw(&state, started.elapsed(), true)
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
                    .draw(&MonitorState::default(), Duration::ZERO, true)
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
                        .draw(&MonitorState::default(), Duration::ZERO, true)
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
