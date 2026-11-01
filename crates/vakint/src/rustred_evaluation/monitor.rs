//! Presentation of native and caller work; never a reduction or cache owner.

use rustred::solver::{CandidateOrderedReductionProgress, CandidateRoutedStoredState};
use std::{
    collections::VecDeque,
    io::{self, Write},
    path::PathBuf,
    time::Duration,
};

use model::{
    Collection, Counter, Details, DfsEstimates, DfsModel, EstimateState, ExpressionSize,
    GuardSummary, MemoryModel, MonitorModel, NATIVE_STAGE, OrderedModel, PhaseKind, PhaseModel,
    PhaseRecord, Rate, RunConfig, StageKind, StageModel, StageTotal,
};

mod model;
mod publication;
mod runtime;
#[cfg(test)]
mod tests;
#[cfg(feature = "progress-tui")]
mod tui;
mod viewer;
pub(crate) use runtime::{CallerProgress, MonitorContext};
pub use runtime::{RustRedMonitorScope, RustRedMonitorWork};
pub use viewer::RustRedSnapshotViewer;

/// Optional diagnostics for a single call or an explicitly scoped caller job.
/// Omission leaves observation disabled.
#[derive(Clone, Debug)]
pub struct RustRedMonitorOptions {
    pub display: RustRedMonitorDisplay,
    /// Atomic latest native JSON, separate from cache/checkpoint reports.
    /// A path enables JSON even when terminal display is off.
    pub progress_snapshot: Option<PathBuf>,
    /// Minimum native safe-boundary interval, not a CAS heartbeat deadline.
    pub interval: Duration,
    pub no_color: bool,
}

impl Default for RustRedMonitorOptions {
    fn default() -> Self {
        Self {
            display: RustRedMonitorDisplay::Off,
            progress_snapshot: None,
            interval: Duration::from_secs(1),
            no_color: false,
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
    lease: Option<publication::FileLease>,
}

impl LatestJsonFile {
    #[cfg(test)]
    pub(super) fn from_options(options: &RustRedMonitorOptions) -> Self {
        Self::from_owner(options, publication::FileLease::next_owner())
    }

    fn from_owner(options: &RustRedMonitorOptions, owner: u64) -> Self {
        Self {
            destination: options.progress_snapshot.clone(),
            lease: options
                .progress_snapshot
                .as_deref()
                .map(|path| publication::FileLease::with_owner(path, owner)),
        }
    }

    #[cfg(test)]
    pub(super) fn write(&self, native: &serde_json::Value) -> io::Result<()> {
        self.write_if(native, &std::sync::atomic::AtomicBool::new(true))
    }

    fn write_if(
        &self,
        native: &serde_json::Value,
        active: &std::sync::atomic::AtomicBool,
    ) -> io::Result<()> {
        self.write_with_check(
            |file| {
                serde_json::to_writer(&mut *file, native)?;
                file.write_all(b"\n")
            },
            active,
        )
    }

    #[cfg(test)]
    fn write_with(
        &self,
        write: impl FnOnce(&mut std::fs::File) -> io::Result<()>,
    ) -> io::Result<()> {
        self.write_with_check(write, &std::sync::atomic::AtomicBool::new(true))
    }

    fn write_with_check(
        &self,
        write: impl FnOnce(&mut std::fs::File) -> io::Result<()>,
        active: &std::sync::atomic::AtomicBool,
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
            // Atomic replacement is enough for a progress file: readers see
            // the old or the new document, never a partial one. No fsync.
            temporary.as_file_mut().flush()?;
            let lease = self
                .lease
                .as_ref()
                .expect("destination owns a publication lease");
            let _commit = lease.commit();
            // Recheck AFTER write/flush. Superseded or detached publishers
            // must never commit a delayed old snapshot over a newer owner.
            if !active.load(std::sync::atomic::Ordering::Acquire) || !lease.is_current() {
                return Ok(());
            }
            // Sibling persist replaces atomically. On failure its owned tempfile
            // is dropped and removed; no fallible step follows publication.
            temporary
                .persist(destination)
                .map(|_| ())
                .map_err(|error| error.error)
        }
        #[cfg(not(all(unix, not(target_arch = "wasm32"))))]
        {
            let _ = (destination, write, active);
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
}

#[cfg(any(feature = "progress-tui", test))]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum UiAction {
    Detach,
    /// Ctrl-C: restore the terminal, then deliver SIGINT as a terminal would.
    Interrupt,
    ToggleFreeze,
    ToggleHelp,
    ToggleDetails,
    ScrollUp,
    ScrollDown,
    PageUp,
    PageDown,
}

/// Display-only state of one terminal attachment.
#[cfg(any(feature = "progress-tui", test))]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct UiView {
    pub frozen: bool,
    pub show_help: bool,
    pub show_details: bool,
    pub scroll: usize,
}

/// A display projection, not an observer event or serialized native contract.
/// The native adapter copies these bounded, owned values inside its callback.
#[derive(Clone, Debug, Default)]
pub(super) struct DisplaySnapshot {
    pub scope_kind: ScopeKind,
    pub api_status: ApiStatus,
    pub native_epoch: u64,
    pub run: Option<RunConfig>,
    /// The current native envelope's key and degrees, never maxima or history.
    pub current_key: Option<String>,
    pub current_degrees: Option<(u128, u128, i128)>,
    /// Identity/time of the last coherent counters + bounded frontier sample.
    /// Activity notices and forwarded summaries never advance these fields.
    pub quantitative_sample_at: Option<Duration>,
    pub quantitative_sample_sequence: u64,
    pub native_collection: Option<NativeCollectionSample>,
    pub native_collection_sample_at: Option<Duration>,
    pub ordered: Option<OrderedReductionSample>,
    pub sequence: u64,
    pub elapsed: Duration,
    /// The caller phase, or `native <Activity>` for native envelopes.
    pub phase: String,
    pub phase_started_at: Duration,
    /// Cumulative wall time of every finished phase label, in first-seen
    /// order. The running phase is added on presentation; labels are bounded.
    pub phase_ledger: Vec<(String, Duration)>,
    /// The caller phase identifier or [`NATIVE_STAGE`]: unlike `phase`, all
    /// native activities share one stage and one clock.
    pub stage: String,
    pub stage_started_at: Duration,
    /// Cumulative wall time of every finished stage, per stage kind.
    pub stage_ledger: Vec<(StageKind, Duration)>,
    pub native_work_active: bool,
    pub caller_work: Option<CallerWork>,
    pub ordered_batches: Option<CallerWork>,
    pub reduction_cohorts: Option<CallerWork>,
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
    pub dfs_cache: Option<DfsCacheCounts>,
    pub path: Vec<String>,
    pub events: Vec<String>,
    pub resources: ResourceUsage,
    pub resource_sample_at: Option<Duration>,
    pub native_observed_at: Option<Duration>,
    pub stored_state: Option<StoredStateSample>,
    pub caller_expression: Option<CallerExpressionSample>,
}

/// Depth-first cache counters of the last coherent native sample.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct DfsCacheCounts {
    pub persistent_hits: u64,
    pub persistent_misses: u64,
    pub ram_hits: u64,
    pub cached_entries: u64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub(super) enum ScopeKind {
    #[default]
    Api,
    Job,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub(super) enum ApiStatus {
    #[default]
    NotStarted,
    Running,
    Succeeded,
    Failed,
}

/// Cheap actual collection activity, independent of the frontier/rate census.
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub(super) struct NativeCollectionSample {
    pub completed_roots: Option<u64>,
    pub total_roots: Option<u64>,
    /// One-based item being assembled, never a completed contribution count.
    pub contribution_index: Option<usize>,
    pub contribution_total: Option<usize>,
    pub source: String,
    pub sequence: u64,
    pub sample_age_seconds: Option<f64>,
}

/// Caller loop units are distinct from native integrals and frontier metrics.
#[derive(Clone, Debug)]
pub(super) struct CallerWork {
    pub phase: String,
    pub completed: usize,
    pub total: usize,
    pub unit: String,
    pub sampled_at: Duration,
}

/// Separate owner scopes, never a globally combined algebraic expression.
#[derive(Clone, Debug)]
pub(super) struct StoredStateSample {
    pub counts: CandidateRoutedStoredState,
    pub family_fingerprint: String,
    pub sampled_at: Duration,
}

#[derive(Clone, Debug)]
pub(super) struct CallerExpressionSample {
    pub encoded_bytes: usize,
    pub stored_summands: Option<usize>,
    pub phase: String,
    pub sampled_at: Duration,
}

/// Counters of one native request when its current phase started.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct PhaseStart {
    pub at: Duration,
    pub processed: usize,
    pub hits: usize,
    pub misses: usize,
}

impl PhaseStart {
    fn of(counts: &CandidateOrderedReductionProgress, at: Duration) -> Self {
        Self {
            at,
            processed: counts.processed_keys,
            hits: counts.step_cache_hits,
            misses: counts.step_cache_misses,
        }
    }
}

/// Finished phases retained per request.
const PHASE_LOG: usize = 16;
/// Span of the recent-rate window and its minimum spacing between entries.
const RATE_WINDOW: Duration = Duration::from_secs(30);
const RATE_SPACING: Duration = Duration::from_millis(500);

/// One exact native frontier census, separate from caller expressions and from
/// previous requests in a monitor spanning a complete job. The per-phase
/// bookkeeping relies on RustRed delivering one census at every phase end.
#[derive(Clone, Debug)]
pub(super) struct OrderedReductionSample {
    pub counts: CandidateOrderedReductionProgress,
    pub native_epoch: u64,
    pub family_fingerprint: String,
    pub sampled_at: Duration,
    pub phase_start: PhaseStart,
    /// The most recent finished phases with work, oldest first.
    pub phases: Vec<PhaseRecord>,
    /// First census of the request: (time, processed keys).
    pub first: (Duration, usize),
    /// Recent censuses (time, processed keys), spaced and bounded.
    pub window: VecDeque<(Duration, usize)>,
}

impl OrderedReductionSample {
    /// The next census of a request, carrying the phase log, the phase start
    /// and the rate window over from the previous census of the same request.
    pub(super) fn next(
        previous: Option<&Self>,
        counts: CandidateOrderedReductionProgress,
        native_epoch: u64,
        family_fingerprint: &str,
        sampled_at: Duration,
    ) -> Self {
        let mut sample = Self {
            counts,
            native_epoch,
            family_fingerprint: family_fingerprint.to_owned(),
            sampled_at,
            phase_start: PhaseStart::of(&counts, sampled_at),
            phases: Vec::new(),
            first: (sampled_at, counts.processed_keys),
            window: VecDeque::from([(sampled_at, counts.processed_keys)]),
        };
        let Some(previous) = previous.filter(|previous| {
            previous.native_epoch == native_epoch
                && previous.counts.completed_phases <= counts.completed_phases
                && previous.counts.processed_keys <= counts.processed_keys
        }) else {
            return sample;
        };
        sample.phases.clone_from(&previous.phases);
        sample.first = previous.first;
        sample.window.clone_from(&previous.window);
        if sample
            .window
            .back()
            .is_none_or(|(at, _)| sampled_at.saturating_sub(*at) >= RATE_SPACING)
        {
            sample.window.push_back((sampled_at, counts.processed_keys));
        }
        while sample.window.len() > 2
            && sample
                .window
                .get(1)
                .is_some_and(|(at, _)| sampled_at.saturating_sub(*at) >= RATE_WINDOW)
        {
            sample.window.pop_front();
        }
        if previous.counts.completed_phases == counts.completed_phases {
            sample.phase_start = previous.phase_start;
            return sample;
        }
        // RustRed reports every phase end; this census closes the phase the
        // previous census observed (intervening phases, if any, were empty).
        let closed = previous.counts;
        let start = previous.phase_start;
        let record = PhaseRecord {
            index: closed.completed_phases,
            kind: PhaseKind::of(
                closed.phase_kind,
                closed.owner_index,
                closed.completed_phases,
                closed.total_phases,
            ),
            support: closed.active_support,
            owner: closed.owner_index,
            seconds: sampled_at.saturating_sub(start.at).as_secs_f64(),
            keys: counts.processed_keys.saturating_sub(start.processed) as u64,
            step_cache_hits: counts.step_cache_hits.saturating_sub(start.hits) as u64,
            step_cache_misses: counts.step_cache_misses.saturating_sub(start.misses) as u64,
            frontier_keys: counts.frontier_keys as u64,
            frontier_monomials: counts.frontier_coefficient_monomials as u64,
        };
        if record.keys > 0 {
            sample.phases.push(record);
            let excess = sample.phases.len().saturating_sub(PHASE_LOG);
            sample.phases.drain(..excess);
        }
        sample
    }

    /// Step-cache hits and misses of the current phase.
    pub(super) fn phase_step_cache(&self) -> (usize, usize) {
        (
            self.counts
                .step_cache_hits
                .saturating_sub(self.phase_start.hits),
            self.counts
                .step_cache_misses
                .saturating_sub(self.phase_start.misses),
        )
    }

    /// Keys over the recent window and over the whole request, `silence`
    /// seconds after this census (a silent core lowers the rates).
    fn rates(&self, silence: f64) -> (Option<Rate>, Option<Rate>) {
        let rate = |(at, processed): (Duration, usize)| {
            Rate::over(
                self.counts.processed_keys.saturating_sub(processed) as u64,
                self.sampled_at.saturating_sub(at).as_secs_f64() + silence,
            )
        };
        (
            self.window.front().copied().and_then(rate),
            rate(self.first),
        )
    }

    /// The model `census_age` seconds after this census. `phase_seconds` is
    /// the time in the unfinished phase, `active` whether the core works on it.
    fn model(&self, phase_seconds: Option<f64>, active: bool, census_age: f64) -> OrderedModel {
        let c = &self.counts;
        let (phase_hits, phase_misses) = self.phase_step_cache();
        // Only a working core accumulates silence; a stopped one keeps the
        // rates of its last census.
        let silence = if active { census_age } else { 0.0 };
        let (recent_rate, request_rate) = self.rates(silence);
        // A current rate only exists while the native core is working.
        let recent_rate = recent_rate.filter(|_| active);
        OrderedModel {
            request: self.native_epoch,
            phase: PhaseModel {
                index: c.completed_phases,
                total: c.total_phases,
                kind: PhaseKind::of(
                    c.phase_kind,
                    c.owner_index,
                    c.completed_phases,
                    c.total_phases,
                ),
                support: c.active_support,
                owner: c.owner_index,
                seconds: phase_seconds,
                active,
                keys: c.processed_keys.saturating_sub(self.phase_start.processed) as u64,
                step_cache_hits: phase_hits as u64,
                step_cache_misses: phase_misses as u64,
            },
            input_keys: c.input_keys as u64,
            processed_keys: c.processed_keys as u64,
            recent_rate,
            request_rate,
            frontier: model::Census {
                keys: c.frontier_keys as u64,
                peak_keys: c.peak_frontier_keys as u64,
                monomials: c.frontier_coefficient_monomials as u64,
                peak_monomials: c.peak_frontier_coefficient_monomials as u64,
            },
            output: model::Census {
                keys: c.terminal_keys as u64,
                peak_keys: c.terminal_keys as u64,
                monomials: c.terminal_coefficient_monomials as u64,
                peak_monomials: c.peak_terminal_coefficient_monomials as u64,
            },
            grouping: model::Grouping {
                emitted: c.emitted_terms as u64,
                merged: c.merged_terms as u64,
                cancelled: c.cancelled_terms as u64,
                output_merged: c.terminal_merges as u64,
                output_cancelled: c.terminal_cancellations as u64,
            },
            step_cache_hits: c.step_cache_hits as u64,
            step_cache_misses: c.step_cache_misses as u64,
            phases: self.phases.clone(),
            census_age_seconds: census_age,
        }
    }
}

/// OS observations, independent of native progress or coefficient payloads.
#[derive(Clone, Debug, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub(super) struct ResourceUsage {
    pub process_pid: Option<u32>,
    pub process_start_ticks: Option<u64>,
    #[serde(default)]
    pub process_identity_verified: bool,
    pub cpu_seconds: Option<f64>,
    pub cpu_percent: Option<f64>,
    pub rss_bytes: Option<u64>,
    pub peak_rss_bytes: Option<u64>,
    pub mem_total_bytes: Option<u64>,
    pub mem_available_bytes: Option<u64>,
    pub mem_free_bytes: Option<u64>,
    pub source: Option<String>,
    pub sample_age_seconds: Option<f64>,
    #[serde(default)]
    pub guard: ResourceGuardView,
}

/// A read-only report of a separately owned, explicitly bound local guard.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub(super) struct ResourceGuardView {
    pub status: String,
    pub threshold_bytes: Option<u64>,
    pub sample_seconds: Option<u64>,
    pub report_age_seconds: Option<f64>,
    pub reported_latched_abort: Option<bool>,
    pub guard_pid: Option<u64>,
    pub current_process_covered: bool,
    pub registered_worker: Option<serde_json::Value>,
    pub registered_slots: Vec<String>,
    pub source: Option<String>,
    pub detail: Option<String>,
}

impl Default for ResourceGuardView {
    fn default() -> Self {
        Self {
            status: model::GUARD_UNKNOWN.into(),
            threshold_bytes: None,
            sample_seconds: None,
            report_age_seconds: None,
            reported_latched_abort: None,
            guard_pid: None,
            current_process_covered: false,
            registered_worker: None,
            registered_slots: Vec::new(),
            source: None,
            detail: Some("no explicit task guard binding".into()),
        }
    }
}

impl ResourceUsage {
    /// The memory line of the model. Process values require a verified live
    /// worker identity; host values are always current.
    pub(super) fn memory(&self) -> MemoryModel {
        let verified = self.process_identity_verified;
        let guard = &self.guard;
        MemoryModel {
            verified,
            rss_bytes: self.rss_bytes.filter(|_| verified),
            peak_rss_bytes: self.peak_rss_bytes.filter(|_| verified),
            cpu_percent: self.cpu_percent.filter(|_| verified),
            available_bytes: self.mem_available_bytes,
            total_bytes: self.mem_total_bytes,
            sample_age_seconds: self.sample_age_seconds,
            guard: (guard.status != model::GUARD_UNKNOWN || guard.threshold_bytes.is_some()).then(
                || GuardSummary {
                    status: guard.status.clone(),
                    threshold_bytes: guard.threshold_bytes,
                    covered: guard.current_process_covered,
                    report_age_seconds: guard.report_age_seconds,
                    interval_seconds: guard.sample_seconds,
                },
            ),
        }
    }

    #[cfg(any(target_os = "linux", test))]
    fn from_proc(status: &str, meminfo: &str) -> Self {
        let mut sample = Self {
            process_pid: None,
            process_start_ticks: None,
            process_identity_verified: false,
            cpu_seconds: None,
            cpu_percent: None,
            rss_bytes: Self::proc_bytes(status, "VmRSS"),
            peak_rss_bytes: Self::proc_bytes(status, "VmHWM"),
            mem_total_bytes: Self::proc_bytes(meminfo, "MemTotal"),
            mem_available_bytes: Self::proc_bytes(meminfo, "MemAvailable"),
            mem_free_bytes: Self::proc_bytes(meminfo, "MemFree"),
            source: Some("Linux /proc/self/status + /proc/meminfo".into()),
            sample_age_seconds: None,
            guard: Default::default(),
        };
        if [
            sample.rss_bytes,
            sample.peak_rss_bytes,
            sample.mem_total_bytes,
            sample.mem_available_bytes,
            sample.mem_free_bytes,
        ]
        .iter()
        .all(Option::is_none)
        {
            sample.source = None;
        }
        sample
    }

    #[cfg(any(target_os = "linux", test))]
    fn proc_bytes(document: &str, key: &str) -> Option<u64> {
        let mut fields = document.lines().filter_map(|line| {
            let (name, value) = line.split_once(':')?;
            (name == key).then_some(value)
        });
        let value = fields.next()?;
        if fields.next().is_some() {
            return None;
        }
        let mut parts = value.split_whitespace();
        let bytes = parts.next()?.parse::<u64>().ok()?.checked_mul(1024)?;
        (parts.next()? == "kB" && parts.next().is_none()).then_some(bytes)
    }
}

#[derive(Clone, Debug)]
pub(super) struct FrontierCount {
    pub count: u64,
    pub exact: bool,
    /// Only comparable Route/Apply identities may be divided by node throughput.
    pub operational_scope: bool,
}

#[derive(Clone, Copy, Debug)]
struct RateSample {
    sequence: u64,
    elapsed: Duration,
    entered: u64,
    completed: Option<u64>,
    roots: Option<u64>,
}

/// No native callback, solver, registry, coefficients or application is retained.
#[derive(Default)]
pub(super) struct MonitorState {
    latest: DisplaySnapshot,
    received_at: Duration,
    samples: VecDeque<RateSample>,
    last_completion: Option<Duration>,
    #[cfg(any(feature = "progress-tui", test))]
    terminal: TerminalDisplay,
}

/// Display state of the in-process terminal, independent of the model.
#[cfg(any(feature = "progress-tui", test))]
#[derive(Default)]
struct TerminalDisplay {
    frozen: Option<MonitorModel>,
    ui: UiView,
    detached: bool,
    interrupted: bool,
}

impl MonitorState {
    /// Machine-readable state of the monitor. The typed model is always live;
    /// freezing is a property of one terminal display only.
    pub(super) fn presentation(&self, now: Duration) -> serde_json::Value {
        let ordered_reduction = self.latest.ordered.as_ref().map(|sample| {
            let c = sample.counts;
            serde_json::json!({
                "scope":"one native request; Route/Apply identities are distinct",
                "native_request":sample.native_epoch,
                "family_fingerprint":sample.family_fingerprint,
                "sample_age_seconds":self.observation_age(sample.sampled_at, now).as_secs_f64(),
                "completed_phases":c.completed_phases, "total_phases":c.total_phases,
                "phase_kind":format!("{:?}", c.phase_kind),
                "active_support":c.active_support, "owner_index":c.owner_index,
                "input_keys":c.input_keys, "frontier_keys":c.frontier_keys,
                "peak_frontier_keys":c.peak_frontier_keys, "terminal_keys":c.terminal_keys,
                "frontier_coefficient_monomials":c.frontier_coefficient_monomials,
                "peak_frontier_coefficient_monomials":c.peak_frontier_coefficient_monomials,
                "terminal_coefficient_monomials":c.terminal_coefficient_monomials,
                "peak_terminal_coefficient_monomials":c.peak_terminal_coefficient_monomials,
                "coefficient_monomials_scope":"materialized numerator+denominator monomials; excludes guards, scratch, keys, Atom encoding and bytes/RSS",
                "emitted_terms":c.emitted_terms, "merged_terms":c.merged_terms,
                "cancelled_terms":c.cancelled_terms, "terminal_merges":c.terminal_merges,
                "terminal_cancellations":c.terminal_cancellations,
                "processed_keys":c.processed_keys,
                "grouping_scope":"global frontier/output; excludes input and intra-rule coalescing",
                "phase_fraction_is_time_fraction":false,
                "step_cache":{
                    "scope":"Ordered step cache of this native request; misses run the checked step and are recorded when the cache populates",
                    "hits":c.step_cache_hits, "misses":c.step_cache_misses,
                    "phase_hits":sample.phase_step_cache().0,
                    "phase_misses":sample.phase_step_cache().1,
                    // Phases with step-cache activity; `phase_kind` is the
                    // native kind (as above), `kind` the displayed one.
                    "completed_phases":sample.phases.iter()
                        .filter(|phase| phase.step_cache_hits + phase.step_cache_misses > 0)
                        .map(|phase| serde_json::json!({
                            "phase":phase.index,
                            "phase_kind":phase.kind.native().map(|kind| format!("{kind:?}")),
                            "kind":phase.kind,
                            "active_support":phase.support, "owner_index":phase.owner,
                            "hits":phase.step_cache_hits, "misses":phase.step_cache_misses,
                        })).collect::<Vec<_>>(),
                },
            })
        });
        let mut document = serde_json::json!({
            "schema": model::MODEL_SCHEMA,
            "scope_kind": self.latest.scope_kind,
            "api_status": self.latest.api_status,
            "job_finished": (self.latest.scope_kind == ScopeKind::Job).then_some(self.latest.request_finished),
            "phase": self.latest.phase,
            "native_work_active": self.latest.native_work_active,
            "phase_elapsed_seconds": (!self.latest.phase.is_empty()).then(|| self.observation_age(self.latest.phase_started_at, now).as_secs_f64()),
            "phase_seconds": self.phase_seconds(now),
            "event_age_seconds": (!self.latest.phase.is_empty()).then(|| now.saturating_sub(self.received_at).as_secs_f64()),
            "native_observation_age_seconds": self.latest.native_observed_at.map(|at| self.observation_age(at, now).as_secs_f64()),
            "native_quantitative_sample_age_seconds": self.latest.quantitative_sample_at.map(|at| self.observation_age(at, now).as_secs_f64()),
            "native_collection": self.live_collection(now),
            "caller_work": self.latest.caller_work.as_ref().map(|work| serde_json::json!({
                "scope": "current caller loop units; not native integral/frontier counts",
                "phase": work.phase, "completed": work.completed, "total": work.total, "unit": work.unit,
                "sample_age_seconds": now.saturating_sub(work.sampled_at).as_secs_f64(),
            })),
            "status": self.status(),
            "resources": self.live_resources(now),
            "resource_sample_elapsed_seconds": self.latest.resource_sample_at.map(|at| at.as_secs_f64()),
            "stored_symbolic_state": self.latest.stored_state.as_ref().map(|sample| serde_json::json!({
                "scope": "native stored partitions; overlapping occurrences, not a global expression",
                "family_fingerprint": sample.family_fingerprint,
                "sample_age_seconds": self.observation_age(sample.sampled_at, now).as_secs_f64(),
                "unstarted_roots": sample.counts.unstarted_roots,
                "deferred_rule_maps": sample.counts.deferred_rule_maps,
                "deferred_child_integral_occurrences": sample.counts.deferred_child_integral_occurrences,
                "collected_output_integrals": sample.counts.collected_output_integrals,
            })),
            "caller_expression": self.latest.caller_expression.as_ref().map(|sample| serde_json::json!({
                "scope": "materialized caller Atom; encoded bytes exclude symbol metadata",
                "phase": sample.phase,
                "encoded_bytes": sample.encoded_bytes,
                "stored_summands": sample.stored_summands,
                "sample_age_seconds": self.observation_age(sample.sampled_at, now).as_secs_f64(),
            })),
            "model": self.live_model(now),
        });
        document["ordered_reduction"] = ordered_reduction.unwrap_or(serde_json::Value::Null);
        document["ordered_batches"] = self.latest.ordered_batches.as_ref().map_or(
            serde_json::Value::Null,
            |work| serde_json::json!({
                "scope":"ordered loop/mode batches in this evaluation; excludes other backends and finalization",
                "label":work.phase, "completed":work.completed, "total":work.total,
                "sample_age_seconds":now.saturating_sub(work.sampled_at).as_secs_f64(),
                "fraction_is_time_fraction":false,
            }),
        );
        document["reduction_cohorts"] = self.latest.reduction_cohorts.as_ref().map_or(
            serde_json::Value::Null,
            |work| serde_json::json!({
                "scope":"compatible family, mass and spectator carriers within the current collected ordered batch",
                "completed":work.completed, "total":work.total,
                "sample_age_seconds":now.saturating_sub(work.sampled_at).as_secs_f64(),
                "fraction_is_time_fraction":false,
            }),
        );
        document
    }

    /// A native core can finish before saved terminal substitution and caller
    /// post-processing: only the outer request finishing ends the run.
    fn status(&self) -> MonitorStatus {
        if self.latest.status == MonitorStatus::Succeeded && !self.latest.request_finished {
            MonitorStatus::Running
        } else {
            self.latest.status
        }
    }

    pub(super) fn live_resources(&self, now: Duration) -> ResourceUsage {
        let mut resources = self.latest.resources.clone();
        let age = self
            .latest
            .resource_sample_at
            .map(|at| now.saturating_sub(at).as_secs_f64());
        resources.sample_age_seconds = resources.source.as_ref().and(age);
        if let Some(age) = resources.guard.report_age_seconds.as_mut() {
            *age += self
                .latest
                .resource_sample_at
                .map_or(0.0, |at| now.saturating_sub(at).as_secs_f64());
            if model::guard_report_stale(*age, resources.guard.sample_seconds) {
                resources.guard.status = model::GUARD_UNKNOWN.into();
                resources.guard.current_process_covered = false;
                resources.guard.detail =
                    Some("stale guard report; no active coverage claim".into());
            }
        }
        resources
    }

    fn live_collection(&self, now: Duration) -> Option<NativeCollectionSample> {
        self.latest.native_collection.as_ref().map(|sample| {
            let mut sample = sample.clone();
            sample.sample_age_seconds = self
                .latest
                .native_collection_sample_at
                .map(|at| self.observation_age(at, now).as_secs_f64());
            sample
        })
    }

    /// Stage breakdown of this scope: cumulative seconds per phase label,
    /// including the running phase up to `now`.
    fn phase_seconds(&self, now: Duration) -> serde_json::Value {
        let mut ledger = self.latest.phase_ledger.clone();
        if !self.latest.phase.is_empty() {
            let running = self.observation_age(self.latest.phase_started_at, now);
            match ledger
                .iter_mut()
                .find(|(label, _)| *label == self.latest.phase)
            {
                Some((_, total)) => *total += running,
                None => ledger.push((self.latest.phase.clone(), running)),
            }
        }
        serde_json::json!({
            "scope": "cumulative wall seconds per observed phase label (caller phases and native activities are disjoint), first-seen order",
            "phases": ledger.iter().map(|(phase, seconds)| serde_json::json!({
                "phase": phase, "seconds": seconds.as_secs_f64(),
            })).collect::<Vec<_>>(),
        })
    }

    fn observation_age(&self, sampled_at: Duration, now: Duration) -> Duration {
        // Native/caller samples use the observer clock; reception and rendering
        // use the presentation clock. Their Instant origins need not coincide.
        self.latest
            .elapsed
            .saturating_sub(sampled_at)
            .saturating_add(now.saturating_sub(self.received_at))
    }

    /// Time since `started_at` on the observer clock; a finished run stops
    /// its clocks at the final observation.
    fn clock(&self, started_at: Duration, now: Duration) -> Duration {
        if self.status() == MonitorStatus::Running {
            self.observation_age(started_at, now)
        } else {
            self.latest.elapsed.saturating_sub(started_at)
        }
    }

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
        let frozen = matches!(action, UiAction::ToggleFreeze)
            .then(|| self.live_model(now))
            .filter(|_| self.terminal.frozen.is_none());
        let display = &mut self.terminal;
        let ui = &mut display.ui;
        match action {
            UiAction::Detach => display.detached = true,
            UiAction::Interrupt => {
                display.detached = true;
                display.interrupted = true;
            }
            UiAction::ToggleFreeze => display.frozen = frozen,
            UiAction::ToggleHelp => ui.show_help = !ui.show_help,
            UiAction::ToggleDetails => {
                ui.show_details = !ui.show_details;
                ui.scroll = 0;
            }
            UiAction::ScrollUp => ui.scroll = ui.scroll.saturating_sub(1),
            UiAction::ScrollDown => ui.scroll = ui.scroll.saturating_add(1),
            UiAction::PageUp => ui.scroll = ui.scroll.saturating_sub(10),
            UiAction::PageDown => ui.scroll = ui.scroll.saturating_add(10),
        }
        ui.frozen = display.frozen.is_some();
    }

    #[cfg(any(feature = "progress-tui", test))]
    pub(super) fn detached(&self) -> bool {
        self.terminal.detached
    }

    #[cfg(any(feature = "progress-tui", test))]
    pub(super) fn interrupted(&self) -> bool {
        self.terminal.interrupted
    }

    /// Display-only state of the terminal attachment.
    #[cfg(any(feature = "progress-tui", test))]
    pub(super) fn ui(&self) -> UiView {
        self.terminal.ui
    }

    #[cfg(any(feature = "progress-tui", test))]
    pub(super) fn clamp_scroll(&mut self, maximum: usize) {
        self.terminal.ui.scroll = self.terminal.ui.scroll.min(maximum);
    }

    /// The model the terminal shows: live, with frozen native panels while
    /// the display is frozen.
    #[cfg(any(feature = "progress-tui", test))]
    pub(super) fn model(&self, now: Duration) -> MonitorModel {
        let live = self.live_model(now);
        match &self.terminal.frozen {
            Some(frozen) => live.with_frozen_native(frozen),
            None => live,
        }
    }

    pub(super) fn live_model(&self, now: Duration) -> MonitorModel {
        let snapshot = &self.latest;
        let seconds = |at: Duration| self.observation_age(at, now).as_secs_f64();
        let stage_kind = model::classify_stage(&snapshot.stage);
        let stage_seconds = self.clock(snapshot.stage_started_at, now);
        let stages = StageKind::ALL
            .iter()
            .filter_map(|&kind| {
                let total = snapshot
                    .stage_ledger
                    .iter()
                    .filter(|(finished, _)| *finished == kind)
                    .map(|(_, spent)| *spent)
                    .sum::<Duration>()
                    + if kind == stage_kind && !snapshot.stage.is_empty() {
                        stage_seconds
                    } else {
                        Duration::ZERO
                    };
                (!total.is_zero()).then_some(StageTotal {
                    kind,
                    seconds: total.as_secs_f64(),
                })
            })
            .collect::<Vec<_>>();
        let counter = |work: &CallerWork, label: &str| Counter {
            completed: work.completed,
            total: work.total,
            label: label.to_owned(),
        };
        let fingerprint = snapshot
            .ordered
            .as_ref()
            .map(|ordered| ordered.family_fingerprint.as_str())
            .or(snapshot
                .stored_state
                .as_ref()
                .map(|state| state.family_fingerprint.as_str()))
            .unwrap_or_default();
        let resources = self.live_resources(now);
        MonitorModel {
            scope: snapshot.scope_kind,
            status: self.status(),
            api_status: snapshot.api_status,
            elapsed_seconds: self.clock(Duration::ZERO, now).as_secs_f64(),
            // The worker's own model is current; only a supervisor ages it.
            snapshot_age_seconds: None,
            run: snapshot.run.clone(),
            stage: StageModel {
                kind: stage_kind,
                id: snapshot.stage.clone(),
                activity: snapshot
                    .phase
                    .strip_prefix("native ")
                    .filter(|_| snapshot.stage == NATIVE_STAGE)
                    .map(model::humanize),
                seconds: stage_seconds.as_secs_f64(),
                work: snapshot
                    .caller_work
                    .as_ref()
                    .map(|work| counter(work, &work.unit)),
            },
            stages,
            batch: snapshot
                .ordered_batches
                .as_ref()
                .map(|work| counter(work, &work.phase)),
            cohort: snapshot
                .reduction_cohorts
                .as_ref()
                .map(|work| counter(work, "")),
            ordered: snapshot.ordered.as_ref().map(|ordered| {
                // The phase clock runs from the census that opened the phase
                // while the core works on it; a phase the core left unfinished
                // (failure, end of the core, caller work) keeps the time up to
                // its last census.
                let open = ordered.counts.completed_phases < ordered.counts.total_phases;
                let active = open && snapshot.native_work_active;
                let phase_seconds = open.then(|| {
                    if active {
                        self.clock(ordered.phase_start.at, now)
                    } else {
                        ordered.sampled_at.saturating_sub(ordered.phase_start.at)
                    }
                    .as_secs_f64()
                });
                ordered.model(phase_seconds, active, seconds(ordered.sampled_at))
            }),
            dfs: (snapshot.ordered.is_none()
                && (snapshot.quantitative_sample_at.is_some()
                    || snapshot.native_collection.is_some()))
            .then(|| self.dfs_model(now)),
            memory: resources.memory(),
            details: Details {
                event_age_seconds: (!snapshot.phase.is_empty())
                    .then(|| now.saturating_sub(self.received_at).as_secs_f64()),
                native_observation_age_seconds: snapshot.native_observed_at.map(seconds),
                family: model::short_fingerprint(fingerprint),
                current_key: snapshot.current_key.clone(),
                current_degrees: snapshot
                    .current_degrees
                    .map(|(r, a, d)| [r.to_string(), a.to_string(), d.to_string()]),
                caller_expression: snapshot.caller_expression.as_ref().map(|sample| {
                    ExpressionSize {
                        bytes: sample.encoded_bytes,
                        summands: sample.stored_summands,
                        stage: sample.phase.clone(),
                    }
                }),
                worker_pid: resources.process_pid,
            },
            events: snapshot.events.clone(),
        }
    }

    fn dfs_model(&self, now: Duration) -> DfsModel {
        let snapshot = &self.latest;
        let age = now.saturating_sub(self.received_at);
        let cache = snapshot.dfs_cache;
        DfsModel {
            roots_completed: snapshot.completed_roots,
            roots_total: snapshot.total_roots,
            entered: snapshot.unique_entered,
            completed: snapshot.completed_unique,
            attempts: snapshot.completion_attempts,
            active: snapshot.active_nodes,
            pending_frames: snapshot.pending_frames,
            pending_unique: snapshot
                .pending_unique
                .as_ref()
                .map(|pending| pending.count),
            pending_exact: snapshot
                .pending_unique
                .as_ref()
                .is_some_and(|pending| pending.exact),
            persistent_hits: cache.map(|cache| cache.persistent_hits),
            persistent_misses: cache.map(|cache| cache.persistent_misses),
            ram_hits: cache.map(|cache| cache.ram_hits),
            cached_entries: cache.map(|cache| cache.cached_entries),
            since_completion_seconds: self.last_completion.map(|last| {
                snapshot
                    .elapsed
                    .saturating_sub(last)
                    .saturating_add(age)
                    .as_secs_f64()
            }),
            estimates: self.dfs_estimates(age),
            collection: self.live_collection(now).map(|sample| Collection {
                completed_roots: sample.completed_roots,
                total_roots: sample.total_roots,
                contribution: sample.contribution_index,
                contributions: sample.contribution_total,
                age_seconds: sample.sample_age_seconds,
            }),
            path: snapshot.path.clone(),
        }
    }

    /// Throughput estimates of native depth-first work. They extrapolate the
    /// observed window; future discovery and caller work are excluded.
    fn dfs_estimates(&self, age: Duration) -> DfsEstimates {
        let mut estimates = DfsEstimates {
            samples: self.samples.len(),
            sample_age_seconds: self.samples.back().map(|sample| {
                self.latest
                    .elapsed
                    .saturating_sub(sample.elapsed)
                    .saturating_add(age)
                    .as_secs_f64()
            }),
            ..Default::default()
        };
        if !self.latest.native_work_active
            || self.latest.request_finished
            || self.latest.status != MonitorStatus::Running
        {
            estimates.state = EstimateState::Inactive;
            return estimates;
        }
        let (Some(first), Some(last)) = (self.samples.front(), self.samples.back()) else {
            return estimates;
        };
        let observed_window = last.elapsed.saturating_sub(first.elapsed);
        let window = observed_window.as_secs_f64();
        estimates.window_seconds = Some(window);
        if window <= 0.0 {
            return estimates;
        }
        let sample_age = self
            .latest
            .elapsed
            .saturating_sub(last.elapsed)
            .saturating_add(age);
        if sample_age > observed_window.max(Duration::from_secs(1)) {
            estimates.state = EstimateState::Stale;
            return estimates;
        }
        estimates.state = EstimateState::Live;
        if let (Some(start), Some(end)) = (first.completed, last.completed)
            && end >= start
        {
            let completed = end - start;
            estimates.growth_rate = Some(
                (last.entered.saturating_sub(first.entered) as f64 - completed as f64) / window,
            );
            if completed > 0 {
                let rate = completed as f64 / window;
                estimates.completion_rate = Some(rate);
                if let Some(frontier) = &self.latest.pending_unique
                    && frontier.operational_scope
                {
                    estimates.frontier_eta_seconds = finite(frontier.count as f64 / rate);
                }
            }
        }
        if let (Some(start), Some(end), Some(total)) =
            (first.roots, last.roots, self.latest.total_roots)
            && end > start
            && total >= end
        {
            let rate = (end - start) as f64 / window;
            estimates.root_rate = Some(rate);
            estimates.root_eta_seconds = finite((total - end) as f64 / rate);
        }
        estimates
    }
}

fn finite(seconds: f64) -> Option<f64> {
    (seconds.is_finite() && seconds >= 0.0).then_some(seconds)
}
