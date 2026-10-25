//! Supervisor-side presentation of the worker's unchanged display projection.

use std::path::Path;

/// Read-only terminal attachment. Dropping or detaching it never controls a worker.
pub struct RustRedSnapshotViewer {
    #[cfg(feature = "progress-tui")]
    active: Option<ActiveViewer>,
}

impl RustRedSnapshotViewer {
    /// A non-terminal caller remains inactive and does not read the snapshot.
    pub fn new(snapshot_path: &Path, no_color: bool) -> Self {
        #[cfg(feature = "progress-tui")]
        {
            Self {
                active: super::tui::TerminalGuard::enter().ok().map(|terminal| {
                    ActiveViewer::new(
                        terminal,
                        snapshot_path,
                        no_color || std::env::var_os("NO_COLOR").is_some(),
                    )
                }),
            }
        }
        #[cfg(not(feature = "progress-tui"))]
        {
            let _ = (snapshot_path, no_color);
            Self {}
        }
    }

    /// Bind the freshly spawned child before the first refresh of this attachment.
    pub fn bind_worker_pid(&mut self, worker_pid: u32) {
        #[cfg(feature = "progress-tui")]
        if let Some(active) = self.active.as_mut() {
            active.worker_pid = Some(worker_pid);
        }
        #[cfg(not(feature = "progress-tui"))]
        let _ = worker_pid;
    }

    /// Poll display controls and show the latest atomic snapshot. Faults detach
    /// only this attachment; the supervisor retains its ordinary worker policy.
    pub fn refresh(&mut self) {
        #[cfg(feature = "progress-tui")]
        if let Some(active) = self.active.as_mut()
            && !std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| active.refresh()))
                .is_ok_and(|result| result.unwrap_or(false))
        {
            self.active = None;
        }
    }
}

#[cfg(feature = "progress-tui")]
use std::{
    fs::File,
    io::{self, BufReader, Read},
    path::PathBuf,
    time::{Duration, SystemTime, UNIX_EPOCH},
};

#[cfg(feature = "progress-tui")]
use super::{MonitorStatus, MonitorView, UiAction, tui};

#[cfg(feature = "progress-tui")]
#[derive(serde::Deserialize)]
struct SnapshotEnvelope {
    monitor: Option<PublishedView>,
}

#[cfg(feature = "progress-tui")]
#[derive(Clone, serde::Deserialize)]
struct PublishedView {
    view: MonitorView,
    published_unix_seconds: Option<f64>,
    worker_pid: Option<u32>,
}

#[cfg(feature = "progress-tui")]
struct ActiveViewer {
    terminal: tui::TerminalGuard,
    path: PathBuf,
    no_color: bool,
    attached_unix_seconds: Option<f64>,
    worker_pid: Option<u32>,
    latest: PublishedView,
    frozen: Option<PublishedView>,
    show_help: bool,
    scroll: usize,
}

#[cfg(feature = "progress-tui")]
impl ActiveViewer {
    fn open_snapshot(path: &Path) -> io::Result<File> {
        let mut options = std::fs::OpenOptions::new();
        options.read(true);
        #[cfg(all(unix, not(target_arch = "wasm32")))]
        {
            use std::os::unix::fs::OpenOptionsExt;
            options.custom_flags(rustix::fs::OFlags::NONBLOCK.bits() as i32);
        }
        let file = options.open(path)?;
        if !file.metadata()?.is_file() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "worker snapshot must be a regular file",
            ));
        }
        Ok(file)
    }

    fn new(terminal: tui::TerminalGuard, path: &Path, no_color: bool) -> Self {
        Self {
            terminal,
            path: path.to_owned(),
            no_color,
            attached_unix_seconds: SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .ok()
                .map(|now| now.as_secs_f64()),
            worker_pid: None,
            latest: PublishedView {
                view: MonitorView {
                    title: "RustRed · worker snapshot".into(),
                    subtitle: "Waiting for a worker snapshot; no native work observed yet".into(),
                    status: MonitorStatus::Running,
                    frozen: false,
                    show_help: false,
                    scroll: 0,
                    overview: Vec::new(),
                    frontier: Vec::new(),
                    cache: Vec::new(),
                    complexity: Vec::new(),
                    estimates: Vec::new(),
                    path: Vec::new(),
                    events: Vec::new(),
                },
                published_unix_seconds: None,
                worker_pid: None,
            },
            frozen: None,
            show_help: false,
            scroll: 0,
        }
    }

    fn refresh(&mut self) -> io::Result<bool> {
        if let Some(action) = tui::poll_action(Duration::ZERO)? {
            match action {
                UiAction::Detach => return Ok(false),
                UiAction::ToggleFreeze => {
                    self.frozen = if self.frozen.is_some() {
                        None
                    } else {
                        Some(self.latest.clone())
                    };
                }
                UiAction::ToggleHelp => self.show_help = !self.show_help,
                UiAction::ScrollUp => self.scroll = self.scroll.saturating_sub(1),
                UiAction::ScrollDown => self.scroll = self.scroll.saturating_add(1),
                UiAction::PageUp => self.scroll = self.scroll.saturating_sub(10),
                UiAction::PageDown => self.scroll = self.scroll.saturating_add(10),
            }
        }
        match Self::open_snapshot(&self.path) {
            Ok(file) => {
                // This bounds display-file I/O only, never native work or output.
                const MAX_SNAPSHOT_BYTES: u64 = 16 * 1024 * 1024;
                if file.metadata()?.len() > MAX_SNAPSHOT_BYTES {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "worker snapshot exceeds the display read allowance",
                    ));
                }
                let envelope: SnapshotEnvelope =
                    serde_json::from_reader(BufReader::new(file.take(MAX_SNAPSHOT_BYTES)))?;
                if let Some(published) = envelope.monitor
                    && self.worker_pid.is_some()
                    && published.worker_pid == self.worker_pid
                    && published
                        .published_unix_seconds
                        .zip(self.attached_unix_seconds)
                        .is_some_and(|(published, attached)| {
                            published.is_finite() && published >= attached
                        })
                {
                    // Reused paths can still contain a previous worker's status.
                    // The supervisor attaches before spawning its current worker.
                    self.latest = published;
                }
            }
            Err(error) if error.kind() == io::ErrorKind::NotFound => {}
            Err(error) => return Err(error),
        }
        let published = self.frozen.as_ref().unwrap_or(&self.latest);
        let mut view = published.view.clone();
        view.frozen = self.frozen.is_some();
        view.show_help = self.show_help;
        view.scroll = self.scroll;
        if let Some(published_at) = published.published_unix_seconds
            && let Ok(now) = SystemTime::now().duration_since(UNIX_EPOCH)
        {
            let age = now.as_secs_f64() - published_at;
            if age.is_finite() && age >= 0.0 {
                view.subtitle
                    .push_str(&format!(" · snapshot published {age:.1}s ago"));
            }
        }
        self.terminal.draw_view(&view, self.no_color)?;
        Ok(true)
    }
}

#[cfg(all(test, feature = "progress-tui", unix, not(target_arch = "wasm32")))]
mod tests {
    use std::{
        io::Write,
        sync::{
            Arc,
            atomic::{AtomicU64, Ordering},
        },
        thread,
        time::Instant,
    };

    use super::*;

    #[test]
    fn fifo_snapshot_and_its_symlink_are_rejected_without_waiting_for_a_writer() {
        let directory = tempfile::tempdir().unwrap();
        let fifo = directory.path().join("snapshot.fifo");
        rustix::fs::mkfifoat(
            rustix::fs::CWD,
            &fifo,
            rustix::fs::Mode::RUSR | rustix::fs::Mode::WUSR,
        )
        .unwrap();
        let link = directory.path().join("snapshot-link");
        std::os::unix::fs::symlink(&fifo, &link).unwrap();
        for path in [&fifo, &link] {
            let started = Instant::now();
            let error = ActiveViewer::open_snapshot(path).unwrap_err();
            assert_eq!(error.kind(), io::ErrorKind::InvalidData);
            assert!(started.elapsed() < Duration::from_secs(1));
        }
    }

    #[test]
    #[ignore = "requires a fresh controlling PTY for real viewer fault isolation"]
    fn real_pty_fifo_fault_detaches_while_supervisor_work_continues() {
        let directory = tempfile::tempdir().unwrap();
        let fifo = directory.path().join("snapshot.fifo");
        rustix::fs::mkfifoat(
            rustix::fs::CWD,
            &fifo,
            rustix::fs::Mode::RUSR | rustix::fs::Mode::WUSR,
        )
        .unwrap();
        assert!(!crossterm::terminal::is_raw_mode_enabled().unwrap());
        let mut viewer = RustRedSnapshotViewer::new(&fifo, true);
        assert!(viewer.active.is_some());
        viewer.bind_worker_pid(std::process::id());
        let counter = Arc::new(AtomicU64::new(0));
        let worker_counter = Arc::clone(&counter);
        let worker = thread::spawn(move || {
            for _ in 0..20 {
                worker_counter.fetch_add(1, Ordering::Relaxed);
                thread::sleep(Duration::from_millis(1));
            }
        });
        println!("PTY_READY:fifo");
        std::io::stdout().flush().unwrap();
        let started = Instant::now();
        viewer.refresh();
        let detached = viewer.active.is_none();
        let elapsed = started.elapsed();
        worker.join().unwrap();
        assert!(detached);
        assert!(elapsed < Duration::from_secs(1));
        assert_eq!(counter.load(Ordering::Relaxed), 20);
        assert!(!crossterm::terminal::is_raw_mode_enabled().unwrap());
        println!("PTY_DONE:fifo");
        std::io::stdout().flush().unwrap();
    }
}
