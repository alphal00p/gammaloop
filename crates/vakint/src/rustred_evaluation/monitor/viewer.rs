//! Supervisor-side presentation of a worker's published monitor model.

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
        if let Some(active) = self.active.as_mut() {
            let outcome =
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| active.refresh()))
                    .ok()
                    .and_then(Result::ok)
                    .unwrap_or(Refresh::Detach);
            if outcome != Refresh::Continue {
                // Restore the terminal before any signal is delivered.
                self.active = None;
                if outcome == Refresh::Interrupt {
                    super::tui::interrupt_foreground_group();
                }
            }
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
    MonitorStatus, ResourceUsage, UiAction, UiView,
    model::{MODEL_SCHEMA, MemoryModel, MonitorModel},
    runtime::{ResourceBinding, ResourceTask},
    tui,
};

#[cfg(feature = "progress-tui")]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Refresh {
    Continue,
    Detach,
    Interrupt,
}

#[cfg(feature = "progress-tui")]
#[derive(serde::Deserialize)]
struct SnapshotEnvelope {
    monitor: Option<PublishedView>,
}

/// The parts of a worker's `monitor` object the viewer reads.
#[cfg(feature = "progress-tui")]
#[derive(Clone, serde::Deserialize)]
struct PublishedView {
    #[serde(default)]
    schema: Option<String>,
    #[serde(default)]
    model: Option<MonitorModel>,
    published_unix_seconds: Option<f64>,
    worker_pid: Option<u32>,
    #[serde(default)]
    worker_start_ticks: Option<u64>,
    #[serde(default)]
    phase: Option<String>,
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

    /// The attachment's OS sample, with process values only for the
    /// published live worker identity.
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
        }
        live
    }

    /// Memory of the attachment: live host values, live process values of
    /// a verified worker, and the worker's own guard report only then.
    fn memory(&self, sample: &ResourceUsage) -> MemoryModel {
        let live = self.verified_resources(sample);
        let mut memory = live.memory();
        memory.guard = live
            .process_identity_verified
            .then(|| {
                self.model
                    .as_ref()
                    .and_then(|model| model.memory.guard.clone())
            })
            .flatten();
        memory
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

    /// The published model; an older worker without one shows its phase.
    fn model(&self) -> MonitorModel {
        match &self.model {
            Some(model) if self.schema.as_deref() == Some(MODEL_SCHEMA) => model.clone(),
            _ => {
                let id = self.phase.clone().unwrap_or_default();
                MonitorModel {
                    stage: super::model::StageModel {
                        kind: super::model::classify_stage(&id),
                        id,
                        ..Default::default()
                    },
                    events: vec![format!(
                        "worker snapshot schema {} has no {MODEL_SCHEMA} model",
                        self.schema.as_deref().unwrap_or("v1")
                    )],
                    ..Default::default()
                }
            }
        }
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
    latest: Option<Arc<PublishedView>>,
    received_at: Instant,
    publication_age: Option<Duration>,
    resources: Option<(ResourceUsage, Duration)>,
    frozen: Option<MonitorModel>,
    ui: UiView,
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
            latest: None,
            received_at: Instant::now(),
            publication_age: None,
            resources: None,
            frozen: None,
            ui: UiView::default(),
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
                && !self
                    .latest
                    .as_ref()
                    .is_some_and(|latest| Arc::ptr_eq(&published, latest))
            {
                self.latest = Some(published);
                self.received_at = update.received_at;
                self.publication_age = update.publication_age;
            }
        }
        Ok(())
    }

    fn refresh(&mut self) -> io::Result<Refresh> {
        if let Some(action) = tui::poll_action(Duration::ZERO)? {
            match action {
                UiAction::Detach => return Ok(Refresh::Detach),
                UiAction::Interrupt => return Ok(Refresh::Interrupt),
                UiAction::ToggleFreeze => {
                    self.frozen = if self.frozen.is_some() {
                        None
                    } else {
                        Some(self.model(Instant::now()))
                    };
                }
                UiAction::ToggleHelp => self.ui.show_help = !self.ui.show_help,
                UiAction::ToggleDetails => {
                    self.ui.show_details = !self.ui.show_details;
                    self.ui.scroll = 0;
                }
                UiAction::ScrollUp => self.ui.scroll = self.ui.scroll.saturating_sub(1),
                UiAction::ScrollDown => self.ui.scroll = self.ui.scroll.saturating_add(1),
                UiAction::PageUp => self.ui.scroll = self.ui.scroll.saturating_sub(10),
                UiAction::PageDown => self.ui.scroll = self.ui.scroll.saturating_add(10),
            }
        }
        self.receive()?;
        let model = self.model(Instant::now());
        let ui = UiView {
            frozen: self.frozen.is_some(),
            ..self.ui
        };
        let maximum_scroll = self.terminal.draw_model(&model, ui, self.no_color)?;
        self.ui.scroll = self.ui.scroll.min(maximum_scroll);
        Ok(Refresh::Continue)
    }

    /// The worker's model with its clocks advanced to `now` and this
    /// attachment's own memory sample; native panels stay frozen on request.
    fn model(&self, now: Instant) -> MonitorModel {
        let Some(published) = &self.latest else {
            return MonitorModel {
                elapsed_seconds: now.saturating_duration_since(self.attached).as_secs_f64(),
                events: vec!["waiting for the worker's first snapshot".into()],
                ..Default::default()
            };
        };
        let mut model = published.model();
        if let Some(age) = self.publication_age {
            model.advance(
                age.saturating_add(now.saturating_duration_since(self.received_at))
                    .as_secs_f64(),
            );
        }
        let sample = self.resources.as_ref().map_or_else(
            ResourceUsage::default,
            |(resources, sampled_at)| {
                let mut resources = resources.clone();
                resources.sample_age_seconds = resources.source.as_ref().map(|_| {
                    now.saturating_duration_since(self.attached)
                        .saturating_sub(*sampled_at)
                        .as_secs_f64()
                });
                resources
            },
        );
        model.memory = published.memory(&sample);
        // A worker that stopped publishing is not reported as running forever.
        if model.status == MonitorStatus::Running
            && self.publication_age.is_some_and(|age| {
                age.saturating_add(now.saturating_duration_since(self.received_at))
                    >= Duration::from_secs(10)
            })
        {
            model.events.push("no worker snapshot for over 10s".into());
        }
        match &self.frozen {
            Some(frozen) => model.with_frozen_native(frozen),
            None => model,
        }
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

    use super::super::model::MODEL_SCHEMA;
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
        let mut model = super::super::MonitorState::default().live_model(Duration::ZERO);
        model.stage.id = "synthetic_paused_producer".into();
        model.details.event_age_seconds = Some(0.0);
        model.ordered = Some(super::super::model::OrderedModel {
            processed_keys: 7,
            ..Default::default()
        });
        serde_json::json!({"monitor": {
            "schema": MODEL_SCHEMA, "model": model,
            "worker_pid":pid, "published_unix_seconds":published,
            "worker_start_ticks":super::super::runtime::ResourceSampler::current_start_ticks(pid),
            "phase":"synthetic_paused_producer"
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
    fn published_clocks_advance_by_publication_age_without_changing_counters() {
        let mut document = snapshot(42, 10.0);
        document["monitor"]["model"]["stage"]["seconds"] = 2.0.into();
        document["monitor"]["model"]["ordered"]["census_age_seconds"] = 2.0.into();
        document["monitor"]["model"]["dfs"] = serde_json::json!({
            "collection": {"completed_roots": 177, "total_roots": 1750, "age_seconds": 2.0}
        });
        let envelope: SnapshotEnvelope = serde_json::from_value(document).unwrap();
        let published = envelope.monitor.unwrap();
        let mut model = published.model();
        model.advance(36.0);
        assert_eq!(model.stage.seconds, 38.0);
        let ordered = model.ordered.as_ref().unwrap();
        assert_eq!(ordered.census_age_seconds, 38.0);
        assert_eq!(ordered.processed_keys, 7);
        let collection = model.dfs.as_ref().unwrap().collection.as_ref().unwrap();
        assert_eq!(collection.age_seconds, Some(38.0));
        assert_eq!(collection.completed_roots, Some(177));
        // A finished worker's clocks stop; only ages keep growing.
        let mut finished = published.model();
        finished.status = MonitorStatus::Succeeded;
        finished.advance(36.0);
        assert_eq!(finished.stage.seconds, 2.0);
        assert_eq!(finished.ordered.unwrap().census_age_seconds, 38.0);
    }

    #[test]
    fn legacy_snapshot_without_a_model_shows_its_phase_and_schema() {
        let document = serde_json::json!({"monitor": {
            "view": {"title": "RustRed · request"}, "worker_pid": 42,
            "published_unix_seconds": 10.0, "phase": "prepared_input_loading"
        }});
        let envelope: SnapshotEnvelope = serde_json::from_value(document).unwrap();
        let model = envelope.monitor.unwrap().model();
        assert_eq!(model.stage.id, "prepared_input_loading");
        assert_eq!(model.stage.kind, super::super::model::StageKind::DataLoad);
        assert!(model.events[0].contains(MODEL_SCHEMA));
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
    fn unverified_process_values_are_hidden_while_host_memory_refreshes() {
        let envelope: SnapshotEnvelope = serde_json::from_value(snapshot(42, 10.0)).unwrap();
        let mut view = envelope.monitor.unwrap();
        view.worker_start_ticks = None;
        let mut model = view.model.take().unwrap();
        model.memory.guard = Some(super::super::model::GuardSummary {
            status: "active".into(),
            threshold_bytes: Some(1),
            covered: true,
        });
        view.model = Some(model);
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
        let memory = view.memory(&sample);
        assert!(!memory.verified && memory.rss_bytes.is_none() && memory.cpu_percent.is_none());
        assert_eq!(memory.available_bytes, Some(2048));
        // No guard coverage is claimed for an unverified worker.
        assert!(memory.guard.is_none());
        view.worker_start_ticks = Some(100);
        let memory = view.memory(&sample);
        assert!(memory.verified && memory.rss_bytes == Some(991));
        assert_eq!(memory.guard.unwrap().status, "active");
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
        let view = envelope.monitor.unwrap();
        assert!(sampler.bound_start_ticks().is_none());
        assert!(!sampler.sample().process_identity_verified);
        assert!(view.bind_resources(&binding));
        sampler.bind_start_ticks(binding.bound_start_ticks().unwrap());
        let recorded = sampler.sample();
        assert!(view.verified_resources(&recorded).process_identity_verified);
        assert!(view.memory(&recorded).rss_bytes.is_some());
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
        assert!(view.memory(&zombie).rss_bytes.is_none());
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
                && active.latest.is_some()
            {
                let model = active.model(Instant::now());
                saw_help |= active.ui.show_help;
                saw_freeze |= active.frozen.is_some();
                assert_eq!(model.ordered.as_ref().unwrap().processed_keys, 7);
                if started.elapsed() > Duration::from_secs(35) {
                    assert_eq!(model.stage.id, "synthetic_paused_producer");
                    assert!(model.stage.seconds >= 35.0, "stage clock stopped");
                    assert!(
                        model.details.event_age_seconds.unwrap() >= 35.0,
                        "event age stopped during native silence"
                    );
                    assert!(
                        model
                            .events
                            .iter()
                            .any(|event| event.contains("no worker snapshot"))
                    );
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
