//! Per-call presentation of native progress; never a reduction or cache owner.

use std::{
    collections::VecDeque,
    io::{self, Write},
    path::PathBuf,
    time::Duration,
};

mod runtime;
#[cfg(test)]
mod tests;
#[cfg(feature = "progress-tui")]
mod tui;
mod viewer;
pub(crate) use runtime::MonitorContext;
pub use viewer::RustRedSnapshotViewer;

/// Optional diagnostics for a single call. Omission leaves observation disabled.
#[derive(Clone, Debug)]
pub struct RustRedMonitorOptions {
    pub display: RustRedMonitorDisplay,
    /// Atomic latest native JSON, separate from cache/checkpoint reports.
    /// A path enables JSON even when terminal display is off.
    pub progress_snapshot: Option<PathBuf>,
    /// Minimum native safe-boundary interval, not a CAS heartbeat deadline.
    pub interval: Duration,
    pub no_color: bool,
    /// Explicit lossless benchmark metadata collection. Logger verbosity never
    /// enables this diagnostic. Presentation and snapshots must remain off.
    #[cfg(feature = "experimental-rustred")]
    pub cache_hit_diagnostic: bool,
}

impl Default for RustRedMonitorOptions {
    fn default() -> Self {
        Self {
            display: RustRedMonitorDisplay::Off,
            progress_snapshot: None,
            interval: Duration::from_secs(1),
            no_color: false,
            #[cfg(feature = "experimental-rustred")]
            cache_hit_diagnostic: false,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RustRedMonitorDisplay {
    #[default]
    Off,
    Tui,
    Json,
}

impl RustRedMonitorOptions {
    pub fn enabled(&self) -> bool {
        self.display != RustRedMonitorDisplay::Off || self.progress_snapshot.is_some()
    }

    pub fn validate(&self) -> io::Result<()> {
        #[cfg(feature = "experimental-rustred")]
        if self.cache_hit_diagnostic && self.enabled() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "lossless cache-hit diagnostics cannot be combined with presentation or snapshots",
            ));
        }
        if self.display == RustRedMonitorDisplay::Tui && !cfg!(feature = "progress-tui") {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "terminal monitoring requires the progress-tui feature",
            ));
        }
        if self.display == RustRedMonitorDisplay::Json && self.progress_snapshot.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "JSON monitoring requires a progress_snapshot destination",
            ));
        }
        Ok(())
    }
}

/// JSONL stream adapter for the native projection, without a consumer schema.
/// This is not the atomic latest-JSON file sink selected by `progress_snapshot`.
pub(super) struct SnapshotWriter<W>(W);

impl<W: Write> SnapshotWriter<W> {
    pub(super) fn new(writer: W) -> Self {
        Self(writer)
    }

    pub(super) fn write(&mut self, native: &serde_json::Value) -> io::Result<()> {
        serde_json::to_writer(&mut self.0, native)?;
        self.0.write_all(b"\n")?;
        self.0.flush()
    }
}

/// Atomically replace the latest native JSON, without a consumer-owned schema.
pub(super) struct LatestJsonFile {
    destination: Option<PathBuf>,
}

impl LatestJsonFile {
    pub(super) fn from_options(options: &RustRedMonitorOptions) -> Self {
        Self {
            destination: options.progress_snapshot.clone(),
        }
    }

    pub(super) fn write(&self, native: &serde_json::Value) -> io::Result<()> {
        self.write_with(|file| {
            serde_json::to_writer(&mut *file, native)?;
            file.write_all(b"\n")
        })
    }

    fn write_with(
        &self,
        write: impl FnOnce(&mut std::fs::File) -> io::Result<()>,
    ) -> io::Result<()> {
        let Some(destination) = self.destination.as_ref() else {
            return Ok(());
        };
        #[cfg(all(unix, not(target_arch = "wasm32")))]
        {
            let parent = destination
                .parent()
                .filter(|path| !path.as_os_str().is_empty())
                .unwrap_or_else(|| std::path::Path::new("."));
            let mut temporary = tempfile::NamedTempFile::new_in(parent)?;
            write(temporary.as_file_mut())?;
            temporary.as_file_mut().flush()?;
            temporary.as_file().sync_all()?;
            // Sibling persist replaces atomically. On failure its owned tempfile
            // is dropped and removed; no fallible step follows publication.
            temporary
                .persist(destination)
                .map(|_| ())
                .map_err(|error| error.error)
        }
        #[cfg(not(all(unix, not(target_arch = "wasm32"))))]
        {
            let _ = (destination, write);
            Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "atomic monitor snapshots are unavailable on this target",
            ))
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub(super) enum MonitorStatus {
    #[default]
    Running,
    Succeeded,
    Failed,
    Partial,
    Detached,
}

#[cfg(any(feature = "progress-tui", test))]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum UiAction {
    Detach,
    ToggleFreeze,
    ToggleHelp,
    ScrollUp,
    ScrollDown,
    PageUp,
    PageDown,
}

/// A display projection, not an observer event or serialized native contract.
/// The native adapter copies these bounded, owned values inside its callback.
#[derive(Clone, Debug, Default)]
pub(super) struct DisplaySnapshot {
    pub native_epoch: u64,
    /// The current native envelope's key and degrees, never maxima or history.
    pub current_key: Option<String>,
    pub current_degrees: Option<(u128, u128, i128)>,
    /// Identity/time of the last coherent counters + bounded frontier sample.
    /// Activity notices and forwarded summaries never advance these fields.
    pub quantitative_sample_at: Option<Duration>,
    pub quantitative_sample_sequence: u64,
    pub sequence: u64,
    pub elapsed: Duration,
    pub phase: String,
    pub status: MonitorStatus,
    pub request_finished: bool,
    pub unique_entered: u64,
    pub completed_unique: Option<u64>,
    pub completion_attempts: Option<u64>,
    pub active_nodes: u64,
    pub pending_frames: u64,
    pub pending_unique: Option<FrontierCount>,
    pub completed_roots: Option<u64>,
    pub total_roots: Option<u64>,
    pub root_output_terms: Option<u64>,
    pub cache: Vec<(String, String)>,
    pub complexity: Vec<(String, String)>,
    pub path: Vec<String>,
    pub events: Vec<String>,
    pub rss_bytes: Option<u64>,
    pub native_observed_at: Option<Duration>,
}

#[derive(Clone, Debug)]
pub(super) struct FrontierCount {
    pub count: u64,
    pub exact: bool,
    /// Only comparable Route/Apply identities may be divided by node throughput.
    pub operational_scope: bool,
    pub scope: String,
}

#[derive(Clone, Copy, Debug)]
struct RateSample {
    sequence: u64,
    elapsed: Duration,
    entered: u64,
    completed: Option<u64>,
    roots: Option<u64>,
}

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub(super) struct MonitorView {
    pub title: String,
    pub subtitle: String,
    pub status: MonitorStatus,
    pub frozen: bool,
    pub show_help: bool,
    pub scroll: usize,
    pub overview: Vec<(String, String)>,
    pub frontier: Vec<(String, String)>,
    pub cache: Vec<(String, String)>,
    pub complexity: Vec<(String, String)>,
    pub estimates: Vec<(String, String)>,
    pub path: Vec<String>,
    pub events: Vec<String>,
}

/// No native callback, solver, registry, coefficients or application is retained.
#[derive(Default)]
pub(super) struct MonitorState {
    latest: DisplaySnapshot,
    received_at: Duration,
    samples: VecDeque<RateSample>,
    last_completion: Option<Duration>,
    frozen: Option<MonitorView>,
    show_help: bool,
    scroll: usize,
    detached: bool,
}

impl MonitorState {
    pub(super) fn observe(&mut self, snapshot: DisplaySnapshot, now: Duration) {
        let resets = snapshot.native_epoch != self.latest.native_epoch
            || snapshot.elapsed < self.latest.elapsed
            || snapshot.sequence < self.latest.sequence
            || snapshot.quantitative_sample_sequence < self.latest.quantitative_sample_sequence
            || snapshot.unique_entered < self.latest.unique_entered
            || matches!((snapshot.completed_unique, self.latest.completed_unique),
                (Some(next), Some(previous)) if next < previous)
            || matches!((snapshot.completed_roots, self.latest.completed_roots),
                (Some(next), Some(previous)) if next < previous);
        if resets {
            self.samples.clear();
            self.last_completion = None;
        }
        if let Some(sample_at) = snapshot.quantitative_sample_at
            && self
                .samples
                .back()
                .is_none_or(|previous| previous.sequence != snapshot.quantitative_sample_sequence)
        {
            if !resets
                && matches!((snapshot.completed_unique, self.samples.back().and_then(|sample| sample.completed)),
                    (Some(next), Some(previous)) if next > previous)
            {
                self.last_completion = Some(sample_at);
            }
            self.samples.push_back(RateSample {
                sequence: snapshot.quantitative_sample_sequence,
                elapsed: sample_at,
                entered: snapshot.unique_entered,
                completed: snapshot.completed_unique,
                roots: snapshot.completed_roots,
            });
            // This bounds only display history, never reduction work.
            if self.samples.len() > 64 {
                self.samples.pop_front();
            }
        }
        self.latest = snapshot;
        self.received_at = now;
    }

    #[cfg(any(feature = "progress-tui", test))]
    pub(super) fn apply(&mut self, action: UiAction, now: Duration) {
        match action {
            UiAction::Detach => self.detached = true,
            UiAction::ToggleFreeze => {
                self.frozen = if self.frozen.is_some() {
                    None
                } else {
                    Some(self.live_view(now))
                };
            }
            UiAction::ToggleHelp => self.show_help = !self.show_help,
            UiAction::ScrollUp => self.scroll = self.scroll.saturating_sub(1),
            UiAction::ScrollDown => self.scroll = self.scroll.saturating_add(1),
            UiAction::PageUp => self.scroll = self.scroll.saturating_sub(10),
            UiAction::PageDown => self.scroll = self.scroll.saturating_add(10),
        }
    }

    #[cfg(any(feature = "progress-tui", test))]
    pub(super) fn detached(&self) -> bool {
        self.detached
    }

    pub(super) fn view(&self, now: Duration) -> MonitorView {
        let mut view = self.frozen.clone().unwrap_or_else(|| self.live_view(now));
        view.frozen = self.frozen.is_some();
        view.show_help = self.show_help;
        view.scroll = self.scroll;
        if self.detached {
            view.status = MonitorStatus::Detached;
        }
        view
    }

    fn live_view(&self, now: Duration) -> MonitorView {
        let snapshot = &self.latest;
        let age = now.saturating_sub(self.received_at);
        let count = |value: Option<u64>| value.map_or_else(|| "unknown".into(), |n| n.to_string());
        let native_count = |value: u64| {
            if snapshot.native_observed_at.is_some() || snapshot.completed_unique.is_some() {
                value.to_string()
            } else {
                "unobserved".into()
            }
        };
        let pending = snapshot.pending_unique.as_ref().map_or_else(
            || "unknown".into(),
            |value| {
                format!(
                    "{} ({}, {})",
                    value.count,
                    value.scope,
                    if value.exact { "exact" } else { "sampled" }
                )
            },
        );
        let current = snapshot.current_key.as_ref().zip(snapshot.current_degrees);
        let mut complexity = vec![
            (
                "Current canonical native key".into(),
                current.map_or_else(
                    || "unknown / no current native entry".into(),
                    |(key, _)| key.clone(),
                ),
            ),
            (
                "Current canonical R / A / D".into(),
                current.map_or_else(
                    || "unknown".into(),
                    |(_, (r, a, d))| format!("{r} / {a} / {d}"),
                ),
            ),
        ];
        complexity.extend(snapshot.complexity.iter().cloned());
        MonitorView {
            title: "RustRed · request".into(),
            subtitle:
                "Local observations; future discoveries and other backend phases are unobserved"
                    .into(),
            // A native core can finish before saved terminal substitution.
            status: if snapshot.status == MonitorStatus::Succeeded && !snapshot.request_finished {
                MonitorStatus::Running
            } else {
                snapshot.status
            },
            frozen: false,
            show_help: false,
            scroll: 0,
            overview: vec![
                ("Phase".into(), snapshot.phase.clone()),
                ("Elapsed at observation".into(), seconds(snapshot.elapsed)),
                ("Snapshot age".into(), seconds(age)),
                (
                    "Last native observation age".into(),
                    snapshot.native_observed_at.map_or_else(
                        || "unknown / no native dispatch".into(),
                        |at| seconds(snapshot.elapsed.saturating_sub(at).saturating_add(age)),
                    ),
                ),
                (
                    "Request finished".into(),
                    snapshot.request_finished.to_string(),
                ),
                (
                    "Roots collected / selected".into(),
                    format!(
                        "{} / {}",
                        count(snapshot.completed_roots),
                        count(snapshot.total_roots)
                    ),
                ),
                (
                    "Current root output terms".into(),
                    count(snapshot.root_output_terms),
                ),
                (
                    "Process RSS bytes (not coefficient payload)".into(),
                    count(snapshot.rss_bytes),
                ),
            ],
            frontier: vec![
                (
                    "Entered unique Route/Apply identities".into(),
                    native_count(snapshot.unique_entered),
                ),
                (
                    "Completed unique identities (not closure)".into(),
                    count(snapshot.completed_unique),
                ),
                (
                    "Completion attempts (includes recomputation)".into(),
                    count(snapshot.completion_attempts),
                ),
                (
                    "Active unique identities".into(),
                    native_count(snapshot.active_nodes),
                ),
                (
                    "Pending operational frames".into(),
                    native_count(snapshot.pending_frames),
                ),
                ("Pending unique identities".into(), pending),
                (
                    "Since last observed new completion".into(),
                    self.last_completion.map_or_else(
                        || "unknown".into(),
                        |last| seconds(snapshot.elapsed.saturating_sub(last).saturating_add(age)),
                    ),
                ),
            ],
            cache: snapshot.cache.clone(),
            complexity,
            estimates: self.estimates(age),
            path: snapshot.path.clone(),
            events: snapshot.events.clone(),
        }
    }

    fn estimates(&self, age: Duration) -> Vec<(String, String)> {
        let mut rows = vec![
            (
                "Estimate scope".into(),
                "Observed local work only; future discovery excluded".into(),
            ),
            (
                "Estimate sample age".into(),
                self.samples.back().map_or_else(
                    || "unknown".into(),
                    |sample| {
                        seconds(
                            self.latest
                                .elapsed
                                .saturating_sub(sample.elapsed)
                                .saturating_add(age),
                        )
                    },
                ),
            ),
            ("Discovered-frontier ETA".into(), "warmup / unknown".into()),
            ("Root-rate ETA".into(), "warmup / unknown".into()),
            ("Net frontier growth".into(), "unknown".into()),
            (
                "Confidence".into(),
                "Low confidence: heuristic; no remaining-work certificate or global bound".into(),
            ),
        ];
        let (Some(first), Some(last)) = (self.samples.front(), self.samples.back()) else {
            return rows;
        };
        let observed_window = last.elapsed.saturating_sub(first.elapsed);
        let window = observed_window.as_secs_f64();
        rows.push((
            "Observation window".into(),
            format!("{window:.1}s / {} snapshots", self.samples.len()),
        ));
        if window <= 0.0 {
            return rows;
        }
        let sample_age = self
            .latest
            .elapsed
            .saturating_sub(last.elapsed)
            .saturating_add(age);
        if sample_age > observed_window.max(Duration::from_secs(1)) {
            rows[2].1 = "stale / unknown".into();
            rows[3].1 = "stale / unknown".into();
            rows[4].1 = "stale / unknown".into();
            return rows;
        }
        if let (Some(start), Some(end)) = (first.completed, last.completed)
            && end >= start
        {
            let completed = end - start;
            let growth =
                (last.entered.saturating_sub(first.entered) as f64 - completed as f64) / window;
            rows[4].1 = format!(
                "{growth:+.2} identities/s · {}",
                if growth > 0.0 {
                    "growing"
                } else {
                    "not growing in sample"
                }
            );
            if completed > 0 {
                let rate = completed as f64 / window;
                if let Some(frontier) = &self.latest.pending_unique
                    && frontier.operational_scope
                    && let Some(eta) = duration(frontier.count as f64 / rate)
                {
                    rows[2].1 = format!(
                        "{} · pending identities / {rate:.2} completions/s · {} frontier",
                        seconds(eta),
                        if frontier.exact { "exact" } else { "sampled" }
                    );
                }
            }
        }
        if let (Some(start), Some(end), Some(total)) =
            (first.roots, last.roots, self.latest.total_roots)
            && end > start
            && total >= end
        {
            let rate = (end - start) as f64 / window;
            if let Some(eta) = duration((total - end) as f64 / rate) {
                rows[3].1 = format!(
                    "{} · known selected roots / {rate:.2} collected roots/s; root costs vary",
                    seconds(eta)
                );
            }
        }
        rows
    }
}

fn duration(seconds: f64) -> Option<Duration> {
    Duration::try_from_secs_f64(seconds).ok()
}

fn seconds(value: Duration) -> String {
    format!("{:.1}s", value.as_secs_f64())
}
