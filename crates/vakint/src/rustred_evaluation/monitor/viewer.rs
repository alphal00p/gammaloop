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
        if let Some(active) = self.active.as_mut()
            && active.bind_worker(worker_pid).is_err()
        {
            self.active = None;
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
    sync::{
        Arc, Mutex,
        atomic::{AtomicBool, Ordering},
    },
    thread::{self, JoinHandle},
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};

#[cfg(feature = "progress-tui")]
use super::{
    MonitorStatus, MonitorView, ResourceUsage, UiAction,
    runtime::{ResourceBinding, ResourceTask},
    tui,
};

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
    #[serde(default)]
    worker_start_ticks: Option<u64>,
    #[serde(default)]
    phase: Option<String>,
    #[serde(default)]
    native_work_active: Option<bool>,
    #[serde(default)]
    phase_elapsed_seconds: Option<f64>,
    #[serde(default)]
    event_age_seconds: Option<f64>,
    #[serde(default)]
    native_observation_age_seconds: Option<f64>,
    #[serde(default)]
    native_quantitative_sample_age_seconds: Option<f64>,
    #[serde(default)]
    native_document_age_seconds: Option<f64>,
    #[serde(default)]
    resources: Option<ResourceUsage>,
}

#[cfg(feature = "progress-tui")]
#[derive(PartialEq, Eq)]
struct SnapshotStamp {
    modified: Option<SystemTime>,
    length: u64,
    #[cfg(unix)]
    device: u64,
    #[cfg(unix)]
    inode: u64,
}

#[cfg(feature = "progress-tui")]
impl SnapshotStamp {
    fn new(metadata: &std::fs::Metadata) -> Self {
        #[cfg(unix)]
        use std::os::unix::fs::MetadataExt;
        Self {
            modified: metadata.modified().ok(),
            length: metadata.len(),
            #[cfg(unix)]
            device: metadata.dev(),
            #[cfg(unix)]
            inode: metadata.ino(),
        }
    }
}

#[cfg(feature = "progress-tui")]
struct ReaderUpdate {
    published: Option<Arc<PublishedView>>,
    received_at: Instant,
    publication_age: Option<Duration>,
    fault: Option<io::Error>,
}

#[cfg(feature = "progress-tui")]
#[derive(Default)]
struct ReaderState {
    latest: Mutex<Option<ReaderUpdate>>,
    stopped: AtomicBool,
}

#[cfg(feature = "progress-tui")]
impl ReaderState {
    fn publish(&self, update: ReaderUpdate) -> Result<(), Box<ReaderUpdate>> {
        if self.stopped.load(Ordering::Relaxed) {
            return Ok(());
        }
        let Ok(mut slot) = self.latest.try_lock() else {
            return Err(Box::new(update));
        };
        if self.stopped.load(Ordering::Relaxed) {
            return Ok(());
        }
        let previous = slot.replace(update);
        drop(slot);
        // Releasing a large previous snapshot never holds the handoff lock.
        drop(previous);
        Ok(())
    }

    fn take(&self) -> Option<ReaderUpdate> {
        self.latest.try_lock().ok()?.take()
    }
}

/// This task owns no terminal or reducer. Drop requests stop without joining a
/// reader that may still be finishing its current bounded snapshot-file read.
#[cfg(feature = "progress-tui")]
struct SnapshotReader {
    state: Arc<ReaderState>,
    _worker: JoinHandle<()>,
}

#[cfg(feature = "progress-tui")]
impl Drop for SnapshotReader {
    fn drop(&mut self) {
        self.state.stopped.store(true, Ordering::Relaxed);
    }
}

#[cfg(feature = "progress-tui")]
impl SnapshotReader {
    fn new(
        path: PathBuf,
        worker_pid: u32,
        attached: Option<f64>,
        binding: ResourceBinding,
    ) -> io::Result<Self> {
        let mut stamp = None;
        Self::with_reader(worker_pid, attached, binding, move || {
            Self::read(&path, &mut stamp)
        })
    }

    fn with_reader(
        worker_pid: u32,
        attached: Option<f64>,
        binding: ResourceBinding,
        mut read: impl FnMut() -> io::Result<Option<PublishedView>> + Send + 'static,
    ) -> io::Result<Self> {
        let state = Arc::new(ReaderState::default());
        let writer = Arc::clone(&state);
        let worker = thread::Builder::new()
            .name("rustred-snapshot-reader".into())
            .spawn(move || {
                let mut published = None;
                let mut received_at = Instant::now();
                let mut publication_age = None;
                while !writer.stopped.load(Ordering::Relaxed) {
                    let fault = match read() {
                        Ok(Some(view))
                            if view.matches_worker(worker_pid, attached)
                                && view.bind_resources(&binding) =>
                        {
                            received_at = Instant::now();
                            publication_age = view.publication_age();
                            published = Some(Arc::new(view));
                            None
                        }
                        Ok(_) => None,
                        Err(error) => Some(error),
                    };
                    let failed = fault.is_some();
                    if writer.stopped.load(Ordering::Relaxed) {
                        break;
                    }
                    let update = ReaderUpdate {
                        published: published.clone(),
                        received_at,
                        publication_age,
                        fault,
                    };
                    if failed {
                        // A terminal fault must survive handoff contention; it
                        // is not an ordinary coalesced progress sample.
                        let mut pending = update;
                        while let Err(update) = writer.publish(pending) {
                            pending = *update;
                            thread::sleep(Duration::from_millis(10));
                        }
                        break;
                    }
                    let _ = writer.publish(update);
                    thread::sleep(Duration::from_millis(200));
                }
            })?;
        Ok(Self {
            state,
            _worker: worker,
        })
    }

    fn read(path: &Path, stamp: &mut Option<SnapshotStamp>) -> io::Result<Option<PublishedView>> {
        let file = match ActiveViewer::open_snapshot(path) {
            Ok(file) => file,
            Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(None),
            Err(error) => return Err(error),
        };
        let next = SnapshotStamp::new(&file.metadata()?);
        // A byte allowance protects presentation only, never native work.
        const MAX_SNAPSHOT_BYTES: u64 = 16 * 1024 * 1024;
        if next.length > MAX_SNAPSHOT_BYTES {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "worker snapshot exceeds the display read allowance",
            ));
        }
        if stamp.as_ref() == Some(&next) {
            return Ok(None);
        }
        let envelope: SnapshotEnvelope =
            serde_json::from_reader(BufReader::new(file.take(MAX_SNAPSHOT_BYTES)))?;
        *stamp = Some(next);
        Ok(envelope.monitor)
    }
}

#[cfg(feature = "progress-tui")]
impl PublishedView {
    fn bind_resources(&self, binding: &ResourceBinding) -> bool {
        let Some(expected) = self.worker_start_ticks else {
            // A legacy record can still be inspected, but supplies no identity.
            return true;
        };
        // The capability cannot retain the task or adopt a reused PID.
        binding.try_bind_start_ticks(expected)
    }

    fn verified_resources(&self, sample: &ResourceUsage) -> ResourceUsage {
        let mut live = sample.clone();
        if !(sample.process_identity_verified
            && self.worker_pid.is_some()
            && self.worker_start_ticks.is_some()
            && sample.process_pid == self.worker_pid
            && sample.process_start_ticks == self.worker_start_ticks)
        {
            live.process_identity_verified = false;
            live.process_pid = self.worker_pid;
            live.process_start_ticks = self.worker_start_ticks;
            live.cpu_seconds = None;
            live.cpu_percent = None;
            live.rss_bytes = None;
            live.peak_rss_bytes = None;
            live.source = live
                .source
                .as_ref()
                .map(|_| "attachment host OS sample; worker PID/start unverified or exited".into());
        }
        live
    }

    fn resource_rows(
        &self,
        live: &ResourceUsage,
        publication_age: Option<Duration>,
    ) -> Vec<(String, String)> {
        let mut rows = live.rows();
        if !live.process_identity_verified {
            let recorded = self.resources.as_ref().map(ResourceUsage::rows);
            let recorded_age = self
                .resources
                .as_ref()
                .and_then(|sample| sample.sample_age_seconds)
                .and_then(|age| Duration::try_from_secs_f64(age).ok())
                .zip(publication_age)
                .map(|(sample, publication)| sample.saturating_add(publication));
            for (label, value) in &mut rows {
                if matches!(label.as_str(), "Process CPU usage" | "Process RSS / peak") {
                    let previous = recorded
                        .as_ref()
                        .unwrap_or(&self.view.cache)
                        .iter()
                        .find(|(name, _)| name == label);
                    *value = previous.map_or_else(
                        || "unknown; live worker PID/start is unverified".into(),
                        |(_, previous)| {
                            format!(
                                "recorded/unverified · {previous} · {}",
                                recorded_age.map_or_else(
                                    || format!(
                                        "sample age unknown; publication age {}",
                                        publication_age.map_or_else(
                                            || "unknown".into(),
                                            |age| format!("{:.1}s", age.as_secs_f64())
                                        )
                                    ),
                                    |age| format!("sample age {:.1}s", age.as_secs_f64())
                                )
                            )
                        },
                    );
                }
            }
        }
        if let Some((_, detail)) = rows
            .iter_mut()
            .find(|(label, _)| label == "Memory units / raw values")
        {
            *detail = "IEC units; attachment OS samples separate from recorded worker JSON".into();
        }
        rows
    }

    fn guard_rows(
        &self,
        verified: bool,
        publication_age: Option<Duration>,
    ) -> Vec<(String, String)> {
        let mut guard = self
            .resources
            .as_ref()
            .map(|resources| resources.guard.clone())
            .unwrap_or_default();
        if let Some(age) = guard.report_age_seconds.as_mut() {
            if let Some(publication) = publication_age {
                *age += publication.as_secs_f64();
            } else {
                guard.status = "unknown".into();
                guard.current_process_covered = false;
            }
            if guard
                .sample_seconds
                .is_some_and(|interval| *age > 3.0 * interval as f64)
            {
                guard.status = "unknown".into();
                guard.current_process_covered = false;
                guard.detail = Some("stale guard report; no active coverage claim".into());
            }
        }
        if !verified {
            guard.status = "unknown".into();
            guard.current_process_covered = false;
            guard.detail =
                Some("worker PID/start is unverified; no live guard coverage claim".into());
        }
        guard.rows()
    }

    fn matches_worker(&self, worker_pid: u32, attached: Option<f64>) -> bool {
        self.worker_pid == Some(worker_pid)
            && self
                .published_unix_seconds
                .zip(attached)
                .is_some_and(|(published, attached)| published.is_finite() && published >= attached)
    }

    fn publication_age(&self) -> Option<Duration> {
        let now = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .ok()?
            .as_secs_f64();
        Duration::try_from_secs_f64(now - self.published_unix_seconds?).ok()
    }

    fn phase(&self) -> Option<&str> {
        self.phase.as_deref().or_else(|| {
            self.view
                .overview
                .iter()
                .find_map(|(label, value)| (label == "Phase").then_some(value.as_str()))
        })
    }
}

#[cfg(feature = "progress-tui")]
struct ActiveViewer {
    terminal: tui::TerminalGuard,
    path: PathBuf,
    no_color: bool,
    attached: Instant,
    attached_unix_seconds: Option<f64>,
    reader: Option<SnapshotReader>,
    resource_task: Option<ResourceTask>,
    latest: Arc<PublishedView>,
    received_at: Instant,
    publication_age: Option<Duration>,
    phase_observed_at: Instant,
    resources: Option<(ResourceUsage, Duration)>,
    frozen: Option<Arc<PublishedView>>,
    show_help: bool,
    scroll: usize,
}

#[cfg(feature = "progress-tui")]
impl Drop for ActiveViewer {
    fn drop(&mut self) {
        // Stop the reader before restoring the terminal. No join or I/O occurs.
        self.reader = None;
        self.resource_task = None;
    }
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
            attached: Instant::now(),
            attached_unix_seconds: SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .ok()
                .map(|now| now.as_secs_f64()),
            reader: None,
            resource_task: None,
            latest: Arc::new(PublishedView {
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
                worker_start_ticks: None,
                phase: None,
                native_work_active: None,
                phase_elapsed_seconds: None,
                event_age_seconds: None,
                native_observation_age_seconds: None,
                native_quantitative_sample_age_seconds: None,
                native_document_age_seconds: None,
                resources: None,
            }),
            received_at: Instant::now(),
            publication_age: None,
            phase_observed_at: Instant::now(),
            resources: None,
            frozen: None,
            show_help: false,
            scroll: 0,
        }
    }

    fn bind_worker(&mut self, worker_pid: u32) -> io::Result<()> {
        let resources = ResourceTask::new_attachment(worker_pid, self.attached)?;
        let reader = SnapshotReader::new(
            self.path.clone(),
            worker_pid,
            self.attached_unix_seconds,
            resources.binding(),
        )?;
        self.reader = Some(reader);
        self.resource_task = Some(resources);
        Ok(())
    }

    fn receive(&mut self) -> io::Result<()> {
        if let Some(sample) = self.resource_task.as_ref().and_then(ResourceTask::take) {
            self.resources = Some(sample);
        }
        if let Some(update) = self.reader.as_ref().and_then(|reader| reader.state.take()) {
            if let Some(error) = update.fault {
                return Err(error);
            }
            if let Some(published) = update.published
                && !Arc::ptr_eq(&published, &self.latest)
            {
                if published.phase() != self.latest.phase() {
                    self.phase_observed_at = Instant::now();
                }
                self.latest = published;
                self.received_at = update.received_at;
                self.publication_age = update.publication_age;
            }
        }
        Ok(())
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
        self.receive()?;
        let view = self.view(Instant::now());
        let maximum_scroll = self.terminal.draw_view(&view, self.no_color)?;
        self.scroll = self.scroll.min(maximum_scroll);
        Ok(true)
    }

    fn view(&self, now: Instant) -> MonitorView {
        let published = self.frozen.as_ref().unwrap_or(&self.latest);
        let mut view = published.view.clone();
        // Only native panels freeze. Attachment liveness, status, estimates and
        // independently sampled worker resources remain visible and current.
        view.status = self.latest.view.status;
        view.estimates = self.latest.view.estimates.clone();
        view.frozen = self.frozen.is_some();
        view.show_help = self.show_help;
        view.scroll = self.scroll;
        let since_read = now.saturating_duration_since(self.received_at);
        let publication_age = self
            .publication_age
            .map(|age| age.saturating_add(since_read));
        let event_age = self
            .latest
            .event_age_seconds
            .and_then(|age| Duration::try_from_secs_f64(age).ok())
            .zip(publication_age)
            .map(|(age, publication)| age.saturating_add(publication));
        let native_age = self
            .latest
            .native_observation_age_seconds
            .and_then(|age| Duration::try_from_secs_f64(age).ok())
            .zip(publication_age)
            .map(|(age, publication)| age.saturating_add(publication));
        let quantitative_age = self
            .latest
            .native_quantitative_sample_age_seconds
            .and_then(|age| Duration::try_from_secs_f64(age).ok())
            .zip(publication_age)
            .map(|(age, publication)| age.saturating_add(publication));
        let projection_age = self
            .latest
            .native_document_age_seconds
            .and_then(|age| Duration::try_from_secs_f64(age).ok())
            .zip(publication_age)
            .map(|(age, publication)| age.saturating_add(publication));
        let phase_elapsed = self
            .latest
            .phase_elapsed_seconds
            .and_then(|age| Duration::try_from_secs_f64(age).ok())
            .zip(publication_age)
            .map(|(age, publication)| age.saturating_add(publication));
        let seconds = |age: Option<Duration>| {
            age.map_or_else(
                || "unknown".into(),
                |age| format!("{:.1}s", age.as_secs_f64()),
            )
        };
        let phase_clock = phase_elapsed.map_or_else(
            || {
                format!(
                    "observed {:.1}s (attachment only)",
                    now.saturating_duration_since(self.phase_observed_at)
                        .as_secs_f64()
                )
            },
            |elapsed| seconds(Some(elapsed)),
        );
        view.subtitle = format!(
            "{} · active {} · observer age {} · Attachment live {:.1}s",
            self.latest
                .phase()
                .unwrap_or("waiting for worker")
                .replace('_', " "),
            phase_clock,
            seconds(event_age),
            now.saturating_duration_since(self.attached).as_secs_f64()
        );
        let live_rows = [
            (
                "Phase".into(),
                self.latest.phase().unwrap_or("unobserved").into(),
            ),
            ("Worker publication age".into(), seconds(publication_age)),
            ("Active phase elapsed".into(), phase_clock),
            ("Worker observer event age".into(), seconds(event_age)),
            ("Snapshot age".into(), seconds(event_age)),
            ("Last native observation age".into(), seconds(native_age)),
            (
                "Native quantitative sample age".into(),
                seconds(quantitative_age),
            ),
            ("Native projection age".into(), seconds(projection_age)),
            (
                "Progress state".into(),
                if self.latest.published_unix_seconds.is_none() {
                    "awaiting worker observation; indeterminate"
                } else if event_age.is_none_or(|age| age >= Duration::from_secs(1))
                    && view.status == MonitorStatus::Running
                {
                    "no newer worker observer event; indeterminate"
                } else {
                    "worker observations; no global completion certificate"
                }
                .into(),
            ),
        ];
        view.overview
            .retain(|(label, _)| !live_rows.iter().any(|(live, _)| label == live));
        view.overview.extend(live_rows);
        for (label, value) in &self.latest.view.overview {
            if matches!(
                label.as_str(),
                "Caller phase progress"
                    | "Request finished"
                    | "Monitor scope"
                    | "Reduction API status"
                    | "Whole job finished"
            ) && let Some(row) = view.overview.iter_mut().find(|(name, _)| name == label)
            {
                row.1 = value.clone();
            }
        }
        if let Some((_, age)) = view
            .estimates
            .iter_mut()
            .find(|(label, _)| label == "Estimate sample age")
        {
            *age = seconds(quantitative_age);
        }
        if (self.latest.native_work_active != Some(true)
            || quantitative_age.is_none_or(|age| age >= Duration::from_secs(1)))
            && view.status == MonitorStatus::Running
        {
            for (label, estimate) in &mut view.estimates {
                if label.ends_with("ETA") {
                    *estimate = "indeterminate · no newer quantitative sample".into();
                }
            }
        }
        let sample = if let Some((resources, sampled_at)) = &self.resources {
            let mut resources = resources.clone();
            resources.sample_age_seconds = resources.source.as_ref().map(|_| {
                now.saturating_duration_since(self.attached)
                    .saturating_sub(*sampled_at)
                    .as_secs_f64()
            });
            resources
        } else {
            ResourceUsage::default()
        };
        let live = self.latest.verified_resources(&sample);
        view.subtitle.push('\n');
        view.subtitle.push_str(&live.headline());
        let mut rows = self.latest.resource_rows(&live, publication_age);
        rows.extend(
            self.latest
                .guard_rows(live.process_identity_verified, publication_age),
        );
        view.cache
            .retain(|(label, _)| !rows.iter().any(|(live, _)| label == live));
        rows.extend(view.cache);
        view.cache = rows;
        view
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

    struct PausedChild(std::process::Child);

    impl Drop for PausedChild {
        fn drop(&mut self) {
            if self.0.try_wait().ok().flatten().is_none() {
                let _ = self.0.kill();
                let _ = self.0.wait();
            }
        }
    }

    fn paused_child() -> (PausedChild, std::io::BufReader<std::process::ChildStdout>) {
        use std::io::BufRead;
        use std::process::{Command, Stdio};

        let test = format!(
            "{}::real_pty_silent_worker_keeps_live_resources_and_detaches_without_control",
            module_path!().split_once("::").unwrap().1
        );
        let mut child = PausedChild(
            Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    &test,
                    "--ignored",
                    "--nocapture",
                    "--format",
                    "terse",
                ])
                .env("VAKINT_MONITOR_PTY_SCENARIO", "viewer-silence-child")
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::null())
                .spawn()
                .unwrap(),
        );
        let mut output = std::io::BufReader::new(child.0.stdout.take().unwrap());
        let mut line = String::new();
        loop {
            line.clear();
            assert!(output.read_line(&mut line).unwrap() > 0);
            if line.contains("PTY_CHILD_READY:viewer-silence") {
                return (child, output);
            }
        }
    }

    fn snapshot(pid: u32, published: f64) -> serde_json::Value {
        let mut view = super::super::MonitorState::default().view(Duration::ZERO);
        view.frontier = vec![("Synthetic unchanged identities".into(), "7".into())];
        serde_json::json!({"monitor": {
            "view":view, "worker_pid":pid, "published_unix_seconds":published,
            "worker_start_ticks":super::super::runtime::ResourceSampler::current_start_ticks(pid),
            "phase":"synthetic_paused_producer", "phase_elapsed_seconds":0.0,
            "event_age_seconds":0.0, "native_document_age_seconds":0.0,
            "native_observation_age_seconds":0.0, "native_quantitative_sample_age_seconds":0.0
        }})
    }

    fn update() -> ReaderUpdate {
        ReaderUpdate {
            published: None,
            received_at: Instant::now(),
            publication_age: None,
            fault: None,
        }
    }

    #[test]
    fn snapshot_identity_and_attachment_freshness_remain_required() {
        let envelope: SnapshotEnvelope = serde_json::from_value(snapshot(42, 10.0)).unwrap();
        let mut view = envelope.monitor.unwrap();
        assert!(view.matches_worker(42, Some(10.0)));
        assert!(!view.matches_worker(43, Some(10.0)));
        assert!(!view.matches_worker(42, Some(11.0)));
        assert!(!view.matches_worker(42, None));
        view.published_unix_seconds = Some(f64::NAN);
        assert!(!view.matches_worker(42, Some(0.0)));
    }

    #[test]
    fn legacy_process_values_stay_recorded_and_unverified_while_host_memory_refreshes() {
        let envelope: SnapshotEnvelope = serde_json::from_value(snapshot(42, 10.0)).unwrap();
        let mut view = envelope.monitor.unwrap();
        view.worker_start_ticks = None;
        view.resources = Some(ResourceUsage {
            cpu_percent: Some(12.5),
            rss_bytes: Some(44),
            sample_age_seconds: Some(2.0),
            ..Default::default()
        });
        let sample = ResourceUsage {
            process_pid: Some(42),
            process_start_ticks: Some(100),
            process_identity_verified: true,
            cpu_percent: Some(99.0),
            rss_bytes: Some(991),
            mem_available_bytes: Some(2048),
            ..Default::default()
        };
        let live = view.verified_resources(&sample);
        assert!(!live.process_identity_verified);
        assert!(live.cpu_percent.is_none() && live.rss_bytes.is_none());
        assert_eq!(live.mem_available_bytes, Some(2048));
        let rows = view.resource_rows(&live, Some(Duration::from_secs(35)));
        for label in ["Process CPU usage", "Process RSS / peak"] {
            let value = &rows.iter().find(|(key, _)| key == label).unwrap().1;
            assert!(value.contains("recorded/unverified") && value.contains("37.0s"));
            assert!(!value.contains("99.0%") && !value.contains("991 B"));
        }
        let recorded = view.resources.as_mut().unwrap();
        recorded.guard.status = "armed".into();
        recorded.guard.current_process_covered = true;
        recorded.guard.report_age_seconds = Some(0.0);
        recorded.guard.sample_seconds = Some(10);
        assert!(
            view.guard_rows(live.process_identity_verified, Some(Duration::ZERO))
                .iter()
                .any(
                    |(label, value)| label == "Local guard state / current-PID coverage"
                        && value == "unknown / unverified"
                )
        );
        view.resources = None;
        view.view.cache = vec![(
            "Process RSS / peak".into(),
            "original legacy memory row".into(),
        )];
        assert!(
            view.resource_rows(&live, Some(Duration::from_secs(35)))
                .iter()
                .any(|(key, value)| key == "Process RSS / peak"
                    && value.contains("original legacy memory row")
                    && value.contains("unverified"))
        );
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn admitted_start_identity_never_rebinds_to_a_reused_pid() {
        let resources = ResourceTask::new_attachment(u32::MAX, Instant::now()).unwrap();
        let binding = resources.binding();
        let envelope: SnapshotEnvelope = serde_json::from_value(snapshot(u32::MAX, 10.0)).unwrap();
        let mut view = envelope.monitor.unwrap();
        assert!(view.bind_resources(&binding));
        assert!(
            binding.bound_start_ticks().is_none(),
            "legacy snapshot bound a PID"
        );
        view.worker_start_ticks = Some(100);
        assert!(view.bind_resources(&binding));
        view.worker_start_ticks = Some(101);
        assert!(!view.bind_resources(&binding));
        assert_eq!(binding.bound_start_ticks(), Some(100));
        let live = view.verified_resources(&ResourceUsage {
            process_pid: Some(u32::MAX),
            process_start_ticks: Some(100),
            process_identity_verified: true,
            cpu_percent: Some(99.0),
            rss_bytes: Some(991),
            ..Default::default()
        });
        assert!(!live.process_identity_verified);
        assert!(live.cpu_percent.is_none() && live.rss_bytes.is_none());
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn exited_worker_does_not_leave_live_cpu_or_rss() {
        use std::io::Read;

        let (mut child, mut output) = paused_child();
        let mut sampler = super::super::runtime::ResourceSampler::new(child.0.id());
        let resources = ResourceTask::new_attachment(child.0.id(), Instant::now()).unwrap();
        let binding = resources.binding();
        let start =
            super::super::runtime::ResourceSampler::current_start_ticks(child.0.id()).unwrap();
        let envelope: SnapshotEnvelope =
            serde_json::from_value(snapshot(child.0.id(), 10.0)).unwrap();
        let mut view = envelope.monitor.unwrap();
        assert!(sampler.bound_start_ticks().is_none());
        assert!(!sampler.sample().process_identity_verified);
        assert!(view.bind_resources(&binding));
        sampler.bind_start_ticks(binding.bound_start_ticks().unwrap());
        let recorded = sampler.sample();
        assert!(view.verified_resources(&recorded).process_identity_verified);
        view.resources = Some(recorded);
        child
            .0
            .stdin
            .as_mut()
            .unwrap()
            .write_all(b"continue\n")
            .unwrap();
        let mut complete = String::new();
        output.read_to_string(&mut complete).unwrap();
        assert!(complete.contains("PTY_CHILD_DONE:viewer-silence:value=17"));
        let deadline = Instant::now() + Duration::from_secs(1);
        loop {
            let stat = std::fs::read_to_string(format!("/proc/{}/stat", child.0.id())).unwrap();
            if stat
                .rsplit_once(')')
                .and_then(|(_, fields)| fields.split_whitespace().next())
                .is_some_and(|state| matches!(state, "Z" | "X"))
            {
                break;
            }
            assert!(Instant::now() < deadline, "tiny child did not exit");
            thread::sleep(Duration::from_millis(1));
        }
        // The unreaped zombie still has a PID/stat file, but is not a live target.
        let zombie = sampler.sample();
        assert!(!zombie.process_identity_verified);
        assert!(zombie.cpu_percent.is_none() && zombie.rss_bytes.is_none());
        let unverified = view.verified_resources(&zombie);
        assert!(!unverified.process_identity_verified);
        assert!(
            view.resource_rows(&unverified, Some(Duration::from_secs(1)))
                .iter()
                .any(|(label, value)| label == "Process RSS / peak"
                    && value.contains("recorded/unverified"))
        );
        assert!(child.0.wait().unwrap().success());
        let sample = sampler.sample();
        assert!(!sample.process_identity_verified);
        assert!(sample.cpu_percent.is_none() && sample.rss_bytes.is_none());
        assert!(sample.mem_available_bytes.is_some());
        assert_eq!(sampler.bound_start_ticks(), Some(start));
    }

    #[test]
    fn unchanged_atomic_snapshot_is_not_deserialized_again() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("latest.json");
        let sink =
            super::super::LatestJsonFile::from_options(&super::super::RustRedMonitorOptions {
                progress_snapshot: Some(path.clone()),
                ..Default::default()
            });
        sink.write(&snapshot(42, 10.0)).unwrap();
        let mut stamp = None;
        assert!(SnapshotReader::read(&path, &mut stamp).unwrap().is_some());
        assert!(SnapshotReader::read(&path, &mut stamp).unwrap().is_none());
        sink.write(&snapshot(42, 11.0)).unwrap();
        assert_eq!(
            SnapshotReader::read(&path, &mut stamp)
                .unwrap()
                .unwrap()
                .published_unix_seconds,
            Some(11.0)
        );
    }

    #[test]
    fn bounded_handoff_coalesces_and_never_waits_on_the_reader() {
        let state = ReaderState::default();
        let _ = state.publish(update());
        let mut replacement = update();
        replacement.publication_age = Some(Duration::from_secs(123));
        let _ = state.publish(replacement);
        assert_eq!(
            state.take().unwrap().publication_age,
            Some(Duration::from_secs(123))
        );
        assert!(state.take().is_none());
        let locked = state.latest.lock().unwrap();
        let started = Instant::now();
        assert!(state.take().is_none());
        let held = state.publish(update());
        assert!(held.is_err());
        let mut fault = update();
        fault.fault = Some(io::Error::other("contended reader fault"));
        let fault = state.publish(fault).err().unwrap();
        assert!(started.elapsed() < Duration::from_secs(1));
        drop(locked);
        assert!(state.take().is_none());
        assert!(state.publish(*fault).is_ok());
        assert_eq!(
            state.take().unwrap().fault.unwrap().to_string(),
            "contended reader fault"
        );
        state.stopped.store(true, Ordering::Relaxed);
        let _ = state.publish(update());
        assert!(state.take().is_none());
    }

    #[test]
    fn dropping_attachment_does_not_join_an_inflight_reader() {
        let state = Arc::new(ReaderState::default());
        let writer = Arc::clone(&state);
        let (release, paused) = std::sync::mpsc::channel();
        let (done, finished) = std::sync::mpsc::channel();
        let worker = thread::spawn(move || {
            paused.recv().unwrap();
            let _ = writer.publish(update());
            done.send(()).unwrap();
        });
        let reader = SnapshotReader {
            state: Arc::clone(&state),
            _worker: worker,
        };
        let started = Instant::now();
        drop(reader);
        assert!(started.elapsed() < Duration::from_secs(1));
        assert!(state.stopped.load(Ordering::Relaxed));
        release.send(()).unwrap();
        finished.recv_timeout(Duration::from_secs(1)).unwrap();
        assert!(state.take().is_none());
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn blocked_snapshot_read_does_not_block_os_samples_or_detach() {
        let started = Instant::now();
        let pid = std::process::id();
        let published = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_secs_f64();
        let envelope: SnapshotEnvelope = serde_json::from_value(snapshot(pid, published)).unwrap();
        let mut first = envelope.monitor;
        let resources = ResourceTask::new_attachment(pid, started).unwrap();
        let binding = resources.binding();
        let reader_binding = binding.clone();
        let (blocked, in_read) = std::sync::mpsc::channel();
        let (release, wait) = std::sync::mpsc::channel();
        let (returned, read_finished) = std::sync::mpsc::channel();
        let reader = SnapshotReader::with_reader(pid, Some(published), reader_binding, move || {
            if let Some(first) = first.take() {
                return Ok(Some(first));
            }
            blocked.send(()).unwrap();
            wait.recv().unwrap();
            returned.send(()).unwrap();
            Ok(None)
        })
        .unwrap();
        in_read.recv_timeout(Duration::from_secs(2)).unwrap();
        let state = Arc::clone(&reader.state);
        let initial = state.take().unwrap();
        assert_eq!(initial.published.unwrap().worker_pid, Some(pid));
        let deadline = Instant::now() + Duration::from_secs(8);
        let mut samples = Vec::new();
        while samples.len() < 3 {
            if let Some((sample, at)) = resources.take()
                && sample.process_identity_verified
            {
                assert!(sample.rss_bytes.is_some() && sample.mem_available_bytes.is_some());
                assert!(started.elapsed().saturating_sub(at) <= Duration::from_millis(2500));
                assert!(samples.last().is_none_or(|previous| at > *previous));
                samples.push(at);
            }
            assert!(
                Instant::now() < deadline,
                "independent OS sampling stopped during blocked read"
            );
            assert!(state.take().is_none(), "blocked reader fabricated progress");
            thread::sleep(Duration::from_millis(5));
        }
        assert!(samples.last().unwrap().saturating_sub(samples[0]) >= Duration::from_millis(1500));
        let detached_at = Instant::now();
        drop(reader);
        drop(resources);
        assert!(detached_at.elapsed() < Duration::from_secs(1));
        assert!(
            !binding.try_bind_start_ticks(binding.bound_start_ticks().unwrap()),
            "reader retained the sampler owner after detach"
        );
        assert!(state.stopped.load(Ordering::Relaxed));
        release.send(()).unwrap();
        read_finished.recv_timeout(Duration::from_secs(1)).unwrap();
        assert!(
            state.take().is_none(),
            "stopped reader published after its I/O returned"
        );
    }

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
        let deadline = started + Duration::from_secs(1);
        while viewer.active.is_some() && Instant::now() < deadline {
            viewer.refresh();
            thread::sleep(Duration::from_millis(1));
        }
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

    #[test]
    #[ignore = "requires a fresh controlling PTY; synthetic producer is silent for over 35s"]
    fn real_pty_silent_worker_keeps_live_resources_and_detaches_without_control() {
        use std::io::BufRead;

        const SCENARIO: &str = "VAKINT_MONITOR_PTY_SCENARIO";
        if std::env::var(SCENARIO).as_deref() == Ok("viewer-silence-child") {
            println!("PTY_CHILD_READY:viewer-silence");
            std::io::stdout().flush().unwrap();
            let mut command = String::new();
            std::io::stdin().read_line(&mut command).unwrap();
            assert_eq!(command.trim(), "continue");
            println!("PTY_CHILD_DONE:viewer-silence:value=17");
            return;
        }
        assert_eq!(std::env::var(SCENARIO).as_deref(), Ok("viewer-silence"));
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("latest.json");
        let mut viewer = RustRedSnapshotViewer::new(&path, true);
        assert!(viewer.active.is_some());
        let (mut child, child_output) = paused_child();
        viewer.bind_worker_pid(child.0.id());
        let published = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_secs_f64();
        super::super::LatestJsonFile::from_options(&super::super::RustRedMonitorOptions {
            progress_snapshot: Some(path),
            ..Default::default()
        })
        .write(&snapshot(child.0.id(), published))
        .unwrap();
        println!("PTY_READY:viewer-silence");
        std::io::stdout().flush().unwrap();
        let started = Instant::now();
        let mut observed_silence = false;
        let mut saw_help = false;
        let mut saw_freeze = false;
        while viewer.active.is_some() && started.elapsed() < Duration::from_secs(45) {
            viewer.refresh();
            if let Some(active) = viewer.active.as_ref()
                && active.latest.published_unix_seconds.is_some()
            {
                let view = active.view(Instant::now());
                saw_help |= view.show_help;
                saw_freeze |= view.frozen;
                assert_eq!(
                    view.frontier,
                    vec![("Synthetic unchanged identities".into(), "7".into())]
                );
                if started.elapsed() > Duration::from_secs(35) {
                    assert!(view.subtitle.contains("Attachment live"));
                    assert!(
                        view.overview
                            .iter()
                            .any(|(key, value)| key == "Progress state"
                                && value.contains("indeterminate"))
                    );
                    for label in [
                        "Active phase elapsed",
                        "Worker observer event age",
                        "Last native observation age",
                        "Native quantitative sample age",
                    ] {
                        let age = view
                            .overview
                            .iter()
                            .find(|(key, _)| key == label)
                            .unwrap()
                            .1
                            .strip_suffix('s')
                            .unwrap()
                            .parse::<f64>()
                            .unwrap();
                        assert!(age >= 35.0, "{label} stopped during native silence");
                    }
                    #[cfg(target_os = "linux")]
                    {
                        let (resources, at) = active.resources.as_ref().unwrap();
                        assert!(resources.rss_bytes.is_some());
                        assert!(resources.mem_available_bytes.is_some());
                        assert!(
                            active.attached.elapsed().saturating_sub(*at)
                                <= Duration::from_millis(2500)
                        );
                    }
                    observed_silence = true;
                }
            }
            thread::sleep(Duration::from_millis(20));
        }
        let detached = viewer.active.is_none();
        drop(viewer);
        assert!(detached, "PTY driver did not detach the viewer");
        assert!(
            observed_silence,
            "detach happened before the >35s quiet interval"
        );
        assert!(
            saw_help && saw_freeze,
            "PTY driver did not exercise help and freeze"
        );
        assert!(
            child.0.try_wait().unwrap().is_none(),
            "detach controlled the paused child"
        );
        assert!(!crossterm::terminal::is_raw_mode_enabled().unwrap());
        println!("PTY_DETACHED:viewer-silence");
        std::io::stdout().flush().unwrap();
        let mut command = String::new();
        std::io::stdin().read_line(&mut command).unwrap();
        assert_eq!(command.trim(), "continue");
        child
            .0
            .stdin
            .as_mut()
            .unwrap()
            .write_all(b"continue\n")
            .unwrap();
        let mut finished = false;
        for line in child_output.lines() {
            finished |= line
                .unwrap()
                .contains("PTY_CHILD_DONE:viewer-silence:value=17");
        }
        assert!(child.0.wait().unwrap().success() && finished);
        println!("PTY_DONE:viewer-silence");
        std::io::stdout().flush().unwrap();
    }
}
