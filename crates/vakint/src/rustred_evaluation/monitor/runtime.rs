//! One per-call native observer feeds optional presentation and benchmark metadata.

use std::{
    collections::VecDeque,
    io,
    sync::{
        Arc, Mutex,
        atomic::{AtomicBool, Ordering},
        mpsc::{self, SyncSender},
    },
    thread::{self, JoinHandle},
    time::{Duration, Instant},
};

use rustred::solver::{
    CandidateRoutedActivity as Activity, CandidateRoutedApplicationEnvelope as Envelope,
    CandidateRoutedApplicationEvent as Event, CandidateRoutedApplicationObserver as Observer,
    CandidateRoutedApplicationStatus as Status, CandidateRoutedObservationOptions,
    CandidateRoutedObserverControl as Control, CandidateRoutedObserverStatus,
    CandidateRoutedWorkRef,
};
use symbolica::atom::AtomView;

use super::{
    CallerExpressionSample, DisplaySnapshot, FrontierCount, LatestJsonFile, MonitorState,
    MonitorStatus, RustRedMonitorDisplay, RustRedMonitorOptions, SnapshotWriter, StoredStateSample,
};
use crate::VakintError;

#[derive(Default)]
struct Pending {
    display: DisplaySnapshot,
    document: Option<serde_json::Value>,
    dirty: bool,
    caller_phase: String,
    last_document_at: Option<Duration>,
    events: VecDeque<String>,
    fault: Option<String>,
    last_cycle: Option<serde_json::Value>,
    last_node_completion: Option<serde_json::Value>,
    force_publish: bool,
}

struct Shared {
    pending: Mutex<Pending>,
    active: AtomicBool,
    stopping: AtomicBool,
    started: Instant,
}

impl Shared {
    fn pending(&self) -> std::sync::MutexGuard<'_, Pending> {
        self.pending
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    fn fault(&self, message: String) {
        self.pending().fault = Some(message);
        self.active.store(false, Ordering::Relaxed);
    }
}

/// Owned queue capability, not a second observer or a stored callback.
#[derive(Clone)]
struct MonitorHandle {
    shared: Arc<Shared>,
    wake: SyncSender<()>,
    interval: Duration,
}

impl MonitorHandle {
    fn observe(&self, envelope: &Envelope<'_>) -> Result<Control, ()> {
        if !self.shared.active.load(Ordering::Relaxed) {
            return if self.shared.pending().fault.is_some() {
                Err(())
            } else {
                Ok(Control::Detach)
            };
        }
        let mut pending = self.shared.pending();
        if envelope.activity == Activity::Caller
            && let Some(phase) = envelope.caller_phase
        {
            pending.caller_phase = phase.to_owned();
        }
        let initial = MonitorContext::is_initial_snapshot(envelope);
        let admission_start = MonitorContext::is_saved_admission_start(envelope);
        let forced = initial || admission_start
            // A fast native batch still has one coherent final quantitative
            // sample, even when its regular sampling interval has not elapsed.
            || (matches!(envelope.event, Event::CoreFinished)
                && envelope.activity != Activity::Caller
                && envelope.caller_phase.is_none())
            || matches!(
                envelope.event,
                Event::RequestFinished | Event::CycleDetected
            );
        pending.force_publish |= initial
            || admission_start
            || matches!(envelope.event, Event::CycleDetected)
            || (envelope.activity == Activity::Caller && envelope.request_finished);
        let due = forced
            || pending
                .last_document_at
                .is_none_or(|last| envelope.elapsed.saturating_sub(last) >= self.interval);
        if let Event::NodeCompleted(source) = envelope.event
            && let Some(work) = envelope.current
        {
            let key = work.target();
            let (kind, owner) = match work {
                CandidateRoutedWorkRef::Route(_) => ("route", None),
                CandidateRoutedWorkRef::Apply { owner_sector, .. } => (
                    "apply",
                    Some(owner_sector.iter().copied().take(15).collect::<Vec<_>>()),
                ),
            };
            pending.last_node_completion = Some(serde_json::json!({
                "scope": "last observed native completion; not complete event history",
                "sequence": envelope.sequence,
                "elapsed_seconds": envelope.elapsed.as_secs_f64(),
                "family_fingerprint": envelope.family_fingerprint,
                "kind": kind, "owner_sector": owner,
                "powers": key.powers().iter().copied().take(15).collect::<Vec<_>>(),
                "arity": key.powers().len(), "truncated": key.powers().len() > 15,
                "degrees": envelope.degrees.map(|d| serde_json::json!({
                    "R": d.numerator_rank.to_string(), "A": d.positive_power_sum.to_string(), "D": d.power_difference.to_string()
                })),
                "source": format!("{source:?}"), "activity": format!("{:?}", envelope.activity),
            }));
        }
        update_display(&mut pending, envelope, due, initial);
        // Forced activity notices can be very frequent. Their cheap display
        // metadata is coalesced; native JSON is projected only at this cadence.
        if due {
            let document = rustred_app::saved_application_observation_document(envelope);
            if matches!(envelope.event, Event::CycleDetected) {
                pending.last_cycle = Some(document.clone());
            }
            pending.document = Some(document);
            pending.last_document_at = Some(envelope.elapsed);
        }
        pending.dirty = true;
        drop(pending);
        let _ = self.wake.try_send(());
        Ok(Control::Continue)
    }
}

pub(crate) struct MonitorContext {
    handle: Option<MonitorHandle>,
    worker: Option<JoinHandle<()>>,
}

impl MonitorContext {
    pub(crate) fn run<T>(
        options: &RustRedMonitorOptions,
        operation: impl FnOnce(Option<&mut Observer<'_>>) -> Result<T, VakintError>,
    ) -> Result<T, VakintError> {
        // Collection is an explicit per-call option. The private log target
        // controls delivery only; general log verbosity cannot enable an observer.
        #[cfg(feature = "experimental-rustred")]
        let benchmark = options.cache_hit_diagnostic;
        #[cfg(not(feature = "experimental-rustred"))]
        let benchmark = false;
        if benchmark {
            options.validate().map_err(|error| {
                VakintError::EvaluationError(format!("invalid RustRed monitor options: {error}"))
            })?;
        }
        if !options.enabled() && !benchmark {
            return operation(None);
        }
        // Classification-only collection has no presenter, worker or JSON sink.
        let mut context = options.enabled().then(|| Self::new(options));
        let handle = context.as_ref().and_then(|context| context.handle.clone());
        if options.enabled() && handle.is_none() {
            return operation(None);
        }
        let phase_handle = handle.clone();
        let mut callback = move |envelope: &Envelope<'_>| {
            #[cfg(feature = "experimental-rustred")]
            if benchmark {
                Self::log_benchmark_event(envelope);
            }
            match &handle {
                Some(handle) => handle.observe(envelope),
                None => Ok(Control::Continue),
            }
        };
        let (result, observer_status) = {
            let observation = if benchmark {
                CandidateRoutedObservationOptions {
                    interval: Duration::ZERO,
                    max_frontier_frames: 0,
                    max_path_nodes: 0,
                }
            } else {
                CandidateRoutedObservationOptions {
                    interval: options.interval,
                    ..Default::default()
                }
            };
            let mut observer = Observer::new(observation, &mut callback);
            Self::native_phase(Some(&mut observer), "request_start", None);
            let result = operation(Some(&mut observer));
            let phase = phase_handle
                .as_ref()
                .map(|handle| handle.shared.pending().caller_phase.clone())
                .unwrap_or_default();
            let mut final_event = Envelope::new(Event::RequestFinished, "");
            final_event.activity = Activity::Caller;
            final_event.caller_phase = Some(&phase);
            final_event.request_finished = true;
            final_event.status = if result.is_ok() {
                Status::Succeeded
            } else {
                Status::Failed
            };
            final_event.snapshot = observer.last_snapshot();
            final_event.metrics = observer.last_metrics();
            final_event.stored_state = observer.last_stored_state();
            observer.emit(final_event);
            (result, observer.status())
        };
        if let CandidateRoutedObserverStatus::Faulted(fault) = observer_status
            && let Some(handle) = &phase_handle
        {
            handle
                .shared
                .fault(format!("native observer disabled: {fault:?}"));
        }
        #[cfg(feature = "experimental-rustred")]
        if benchmark {
            // A faulted/detached observer suppresses callbacks. Report its
            // actual final state after dropping it, without replacing the math
            // result if a custom logger panics during this last attempt.
            let operation_succeeded = result.is_ok();
            let _ = std::panic::catch_unwind(|| {
                let document =
                    Self::benchmark_finish_document(observer_status, operation_succeeded);
                log::debug!(target: "vakint::benchmark_cache_hits", "{document}");
            });
        }
        if let Some(context) = &mut context {
            context.finish();
        }
        result
    }

    fn is_initial_snapshot(envelope: &Envelope<'_>) -> bool {
        envelope.activity != Activity::Caller
            && envelope.caller_phase.is_none()
            && matches!(envelope.event, Event::Snapshot)
            && envelope
                .snapshot
                .is_some_and(|snapshot| snapshot.sequence == 1)
    }

    fn is_saved_admission_start(envelope: &Envelope<'_>) -> bool {
        // PreparedSavedApplication clears the observer's native summary before
        // this boundary. Preflight can fail without producing a new core sample.
        envelope.activity == Activity::Admission
            && envelope.event == Event::ActivityStarted
            && envelope.caller_phase == Some("entry_admission")
            && envelope.snapshot.is_none()
            && envelope.stored_state.is_none()
    }

    #[cfg(feature = "experimental-rustred")]
    fn benchmark_work_document(work: CandidateRoutedWorkRef<'_>) -> serde_json::Value {
        let (stage, owner) = match work {
            CandidateRoutedWorkRef::Route(_) => ("route", None),
            CandidateRoutedWorkRef::Apply { owner_sector, .. } => ("apply", Some(owner_sector)),
        };
        serde_json::json!({
            "stage": stage,
            "powers": work.target().powers(),
            "owner_sector": owner,
        })
    }

    #[cfg(feature = "experimental-rustred")]
    fn benchmark_event_document(envelope: &Envelope<'_>) -> Option<serde_json::Value> {
        if envelope.activity == Activity::Caller || envelope.caller_phase.is_some() {
            return None;
        }
        let event = match envelope.event {
            Event::Snapshot if Self::is_initial_snapshot(envelope) => "dispatch_start",
            Event::NodeCompleted(rustred::solver::CandidateRoutedCompletionSource::BackingHit) => {
                "backing_hit"
            }
            Event::RootCollected => "root_collected",
            Event::CoreFinished => "dispatch_finish",
            _ => return None,
        };
        let persistent = envelope.snapshot.map(|snapshot| {
            let counters = snapshot.statistics.persistent_cache;
            serde_json::json!({
                "lookups": counters.lookups,
                "route_hits": counters.route_hits,
                "apply_hits": counters.apply_hits,
                "requested_root_hits": counters.requested_root_hits,
                "intermediate_hits": counters.intermediate_hits,
                "misses": counters.misses,
            })
        });
        let root = if matches!(envelope.event, Event::RootCollected) {
            envelope.root.map(Self::benchmark_work_document)
        } else {
            None
        };
        Some(serde_json::json!({
            "schema": "vakint.cache-hit-events.v1",
            "event": event,
            "family_fingerprint": envelope.family_fingerprint,
            "sequence": envelope.sequence,
            "work": envelope.current.map(Self::benchmark_work_document),
            "rootwork": root,
            "requested_roots": envelope.metrics.requested_roots,
            "completed_roots": envelope.metrics.completed_roots,
            "persistent": persistent,
            "status": match envelope.status {
                Status::Running => "running",
                Status::Succeeded => "succeeded",
                Status::Failed => "failed",
            },
        }))
    }

    #[cfg(feature = "experimental-rustred")]
    fn log_benchmark_event(envelope: &Envelope<'_>) {
        if let Some(document) = Self::benchmark_event_document(envelope) {
            log::debug!(target: "vakint::benchmark_cache_hits", "{document}");
        }
    }

    #[cfg(feature = "experimental-rustred")]
    fn benchmark_finish_document(
        observer_status: CandidateRoutedObserverStatus,
        operation_succeeded: bool,
    ) -> serde_json::Value {
        serde_json::json!({
            "schema": "vakint.cache-hit-events.v1",
            "event": "observer_finish",
            "sequence": null,
            "observer_status": format!("{observer_status:?}"),
            "operation_succeeded": operation_succeeded,
            "interval_ns": 0,
        })
    }

    pub(crate) fn native_phase(
        observer: Option<&mut Observer<'_>>,
        phase: &str,
        expression: Option<AtomView<'_>>,
    ) {
        if let Some(observer) = observer
            && observer.status() == CandidateRoutedObserverStatus::Active
        {
            let mut event = Envelope::new(Event::CallerProgress, "");
            event.activity = Activity::Caller;
            event.caller_phase = Some(phase);
            if let Some(expression) = expression {
                // Read existing encoded headers, never format or construct an
                // expression for telemetry. These are caller Atom units, not
                // Native coefficient clone bounds or distinct integral keys.
                event.complexity.terms = Some(expression.nterms());
                event.complexity.encoded_expression_bytes = Some(expression.get_byte_size());
            }
            observer.emit(event);
        }
    }

    fn new(options: &RustRedMonitorOptions) -> Self {
        let shared = Arc::new(Shared {
            pending: Mutex::new(Pending::default()),
            active: AtomicBool::new(true),
            stopping: AtomicBool::new(false),
            started: Instant::now(),
        });
        let (wake, receiver) = mpsc::sync_channel(1);
        let handle = MonitorHandle {
            shared: Arc::clone(&shared),
            wake,
            interval: options.interval,
        };
        let options = options.clone();
        let worker_shared = Arc::clone(&shared);
        match thread::Builder::new()
            .name("vakint-monitor".into())
            .spawn(move || {
                let worker_state = Arc::clone(&worker_shared);
                if std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    worker(worker_state, receiver, options)
                }))
                .is_err()
                {
                    worker_shared.fault("presentation worker panicked".into());
                }
            }) {
            Ok(worker) => Self {
                handle: Some(handle),
                worker: Some(worker),
            },
            Err(error) => {
                log::warn!("RustRed monitor could not start: {error}");
                Self {
                    handle: None,
                    worker: None,
                }
            }
        }
    }

    fn finish(&mut self) {
        if let Some(handle) = &self.handle {
            handle.shared.stopping.store(true, Ordering::Relaxed);
            let _ = handle.wake.try_send(());
        }
        if let Some(worker) = self.worker.take()
            && worker.join().is_err()
        {
            log::warn!("RustRed monitor worker unwound; scientific result is unchanged");
        }
        if let Some(handle) = &self.handle
            && let Some(fault) = &handle.shared.pending().fault
        {
            log::warn!("RustRed monitor telemetry stopped: {fault}");
        }
        self.handle = None;
    }
}

impl Drop for MonitorContext {
    fn drop(&mut self) {
        self.finish();
    }
}

fn worker(shared: Arc<Shared>, receiver: mpsc::Receiver<()>, options: RustRedMonitorOptions) {
    let mut state = MonitorState::default();
    let file = LatestJsonFile::from_options(&options);
    let file_enabled = options.progress_snapshot.is_some();
    let mut stderr_json = false;
    #[cfg(feature = "progress-tui")]
    let mut terminal = if options.display == RustRedMonitorDisplay::Tui {
        match super::tui::TerminalGuard::enter() {
            Ok(terminal) => Some(terminal),
            Err(error) if error.kind() == io::ErrorKind::AlreadyExists => {
                // Another UI owns this terminal; never write JSON into it.
                if !file_enabled {
                    shared.active.store(false, Ordering::Relaxed);
                }
                None
            }
            Err(_) => {
                stderr_json = !file_enabled;
                None
            }
        }
    } else {
        None
    };
    #[cfg(not(feature = "progress-tui"))]
    if options.display == RustRedMonitorDisplay::Tui {
        stderr_json = !file_enabled;
    }
    if options.display == RustRedMonitorDisplay::Json && !file_enabled {
        stderr_json = true;
    }
    let mut stream = SnapshotWriter::new(io::stderr());
    let no_color = options.no_color || std::env::var_os("NO_COLOR").is_some();
    let mut last_document = None;
    let mut last_cycle = None;
    let mut last_node_completion = None;
    let mut document_received = Duration::ZERO;
    let mut last_publish = None::<Instant>;
    let mut last_resource_sample = None::<Instant>;
    // Task-scoped presentation input only, captured once. This never arms a guard.
    let guard_snapshot =
        std::env::var_os("VAKINT_RESOURCE_GUARD_SNAPSHOT").map(std::path::PathBuf::from);
    loop {
        if !shared.stopping.load(Ordering::Relaxed) {
            let _ = receiver.recv_timeout(Duration::from_millis(50));
        }
        let update = {
            let mut pending = shared.pending();
            if pending.dirty {
                pending.dirty = false;
                let force_publish = pending.force_publish;
                pending.force_publish = false;
                Some((
                    pending.display.clone(),
                    pending.document.take(),
                    pending.last_cycle.clone(),
                    pending.last_node_completion.clone(),
                    force_publish,
                ))
            } else {
                None
            }
        };
        let mut force_publish = false;
        if let Some((mut display, document, cycle, completion, force)) = update {
            force_publish = force;
            // A native notice does not reset the independent OS sample or its age.
            display.resources = state.latest.resources.clone();
            display.resource_sample_at = state.latest.resource_sample_at;
            state.observe(display, shared.started.elapsed());
            if let Some(document) = document {
                last_document = Some(document);
                document_received = shared.started.elapsed();
            }
            last_cycle = cycle;
            last_node_completion = completion;
        }
        if force_publish
            || shared.stopping.load(Ordering::Relaxed)
            || last_resource_sample.is_none_or(|last| {
                last.elapsed() >= options.interval.max(Duration::from_millis(100))
            })
        {
            // Presentation-thread I/O only; no native progress/rate sample is added.
            state.latest.resources = resource_usage(guard_snapshot.as_deref());
            state.latest.resource_sample_at =
                cfg!(target_os = "linux").then(|| shared.started.elapsed());
            last_resource_sample = Some(Instant::now());
        }
        #[cfg(feature = "progress-tui")]
        if terminal.is_some() {
            let action = super::tui::poll_action(Duration::ZERO);
            match action {
                Ok(Some(action)) => {
                    state.apply(action, shared.started.elapsed());
                    if state.detached() {
                        terminal = None;
                        if !file_enabled {
                            shared.active.store(false, Ordering::Relaxed);
                        }
                    }
                }
                Ok(None) => {}
                Err(error) => {
                    shared.fault(error.to_string());
                    break;
                }
            }
            if let Some(ui) = terminal.as_mut()
                && let Err(error) = ui.draw(&state, shared.started.elapsed(), no_color)
            {
                shared.fault(error.to_string());
                break;
            }
        }
        #[cfg(not(feature = "progress-tui"))]
        let _ = no_color;
        let refresh_due = last_publish
            .is_none_or(|last| last.elapsed() >= options.interval.max(Duration::from_millis(100)));
        if (force_publish || refresh_due || shared.stopping.load(Ordering::Relaxed))
            && let Some(native) = &last_document
        {
            let mut document = native.clone();
            // Only the nested visual view may be frozen, never status/ETA/resources.
            let mut monitor = state.presentation(shared.started.elapsed());
            monitor["scope"] =
                "observational consumer presentation; native fields unchanged".into();
            monitor["worker_pid"] = std::process::id().into();
            monitor["published_unix_seconds"] = serde_json::json!(
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .ok()
                    .map(|value| value.as_secs_f64())
            );
            monitor["native_document_age_seconds"] = shared
                .started
                .elapsed()
                .saturating_sub(document_received)
                .as_secs_f64()
                .into();
            monitor["last_cycle"] = serde_json::json!(last_cycle);
            monitor["last_node_completion"] = serde_json::json!(last_node_completion);
            monitor["global_remaining_work_certificate"] = false.into();
            document["monitor"] = monitor;
            if let Err(error) = file.write(&document).and_then(|()| {
                if stderr_json {
                    stream.write(&document)
                } else {
                    Ok(())
                }
            }) {
                shared.fault(error.to_string());
                break;
            }
            last_publish = Some(Instant::now());
        }
        if !shared.active.load(Ordering::Relaxed) {
            break;
        }
        // The final callback can arrive while a previous frame is being
        // written. Finish only after its queued final document is drained.
        if shared.stopping.load(Ordering::Relaxed) {
            if !shared.pending().dirty {
                break;
            }
            continue;
        }
        // Bounded coalescing is presentation-only. No renderer/file I/O runs on
        // the native callback thread; input is polled at least every 50ms.
        if let Some(published) = last_publish {
            thread::sleep(Duration::from_millis(50).saturating_sub(published.elapsed()));
        }
    }
}

fn work_label(work: CandidateRoutedWorkRef<'_>) -> String {
    let key = work.target();
    let mut powers = key
        .powers()
        .iter()
        .take(15)
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(",");
    if key.powers().len() > 15 {
        powers.push_str(",… (display truncated)");
    }
    match work {
        CandidateRoutedWorkRef::Route(_) => format!("Route [{powers}]"),
        CandidateRoutedWorkRef::Apply { owner_sector, .. } => format!(
            "Apply owner {:?} [{powers}]",
            &owner_sector[..owner_sector.len().min(15)]
        ),
    }
}

fn update_display(pending: &mut Pending, envelope: &Envelope<'_>, detail: bool, initial: bool) {
    let display = &mut pending.display;
    display.sequence = envelope.sequence;
    display.elapsed = envelope.elapsed;
    // Cheap activity notices identify long native operations even when no new
    // frontier census is sampled. Caller/saved phases have no current native key.
    display.current_key = None;
    display.current_degrees = None;
    if envelope.activity != Activity::Caller
        && envelope.caller_phase.is_none()
        && let Some((work, degrees)) = envelope.current.zip(envelope.degrees)
    {
        display.current_key = Some(work_label(work));
        display.current_degrees = Some((
            degrees.numerator_rank,
            degrees.positive_power_sum,
            degrees.power_difference,
        ));
    }
    if let Some(age) = envelope.native_summary_age {
        display.native_observed_at = Some(envelope.elapsed.saturating_sub(age));
    }
    display.phase = envelope
        .caller_phase
        .map_or_else(|| format!("native {:?}", envelope.activity), str::to_owned);
    let outer_finished = envelope.activity == Activity::Caller && envelope.request_finished;
    display.request_finished = outer_finished;
    display.status = if outer_finished {
        match envelope.status {
            Status::Failed => MonitorStatus::Failed,
            Status::Succeeded
                if display.phase.contains("cancelled") || display.phase.contains("incomplete") =>
            {
                MonitorStatus::Partial
            }
            Status::Succeeded => MonitorStatus::Succeeded,
            Status::Running => MonitorStatus::Running,
        }
    } else {
        MonitorStatus::Running
    };
    let admission_start = MonitorContext::is_saved_admission_start(envelope);
    if initial || admission_start {
        display.native_epoch = display.native_epoch.saturating_add(1);
        display.quantitative_sample_at = None;
        display.quantitative_sample_sequence = 0;
        display.pending_unique = None;
        display.stored_state = None;
    }
    if admission_start {
        display.native_observed_at = None;
        display.unique_entered = 0;
        display.completed_unique = None;
        display.completion_attempts = None;
        display.active_nodes = 0;
        display.pending_frames = 0;
        display.completed_roots = None;
        display.total_roots = None;
        display.root_output_terms = None;
        display.cache.clear();
        display.complexity.clear();
        display.path.clear();
    }
    if envelope.activity == Activity::Caller
        && let Some(phase) = envelope.caller_phase
        && let Some(encoded_bytes) = envelope.complexity.encoded_expression_bytes
    {
        display.caller_expression = Some(CallerExpressionSample {
            encoded_bytes,
            stored_summands: envelope.complexity.terms,
            phase: phase.to_owned(),
            sampled_at: envelope.elapsed,
        });
    }
    if envelope.activity != Activity::Caller
        && envelope.caller_phase.is_none()
        && let Some(counts) = envelope.stored_state
    {
        // These O(1) structural counters are current at their own native safe
        // boundary. A forwarded saved/caller summary never refreshes this age.
        display.stored_state = Some(StoredStateSample {
            counts,
            family_fingerprint: envelope.family_fingerprint.to_owned(),
            sampled_at: envelope.elapsed,
        });
    }
    if detail
        && envelope.activity != Activity::Caller
        && envelope.caller_phase.is_none()
        && let Some(snapshot) = envelope.snapshot
        && let Some(frontier) = envelope.frontier
    {
        // Keep the complete quantitative tuple coherent even if a cheap notice
        // is coalesced before the worker receives this genuine frontier sample.
        display.quantitative_sample_at = Some(envelope.elapsed);
        display.quantitative_sample_sequence = envelope.sequence;
        display.unique_entered = snapshot.unique_nodes as u64;
        display.active_nodes = snapshot.active_nodes as u64;
        display.pending_frames = snapshot.pending_frames as u64;
        display.completed_unique = Some(envelope.metrics.completed_unique_operational_nodes as u64);
        display.completion_attempts = Some(envelope.metrics.completion_attempts as u64);
        display.completed_roots = Some(envelope.metrics.completed_roots as u64);
        display.total_roots = envelope.metrics.requested_roots.map(|v| v as u64);
        display.root_output_terms = envelope.metrics.output_terms.map(|v| v as u64);
        display.pending_unique = Some(FrontierCount {
            count: frontier.unique_nodes as u64,
            exact: frontier.exact,
            operational_scope: true,
            scope: format!(
                "Route/Apply identities; {}/{} frames sampled",
                frontier.frames_examined, frontier.total_pending_frames
            ),
        });
        display.cache = vec![
            (
                "Persistent lookups".into(),
                snapshot.statistics.persistent_cache.lookups.to_string(),
            ),
            (
                "Persistent requested-root hits".into(),
                snapshot
                    .statistics
                    .persistent_cache
                    .requested_root_hits
                    .to_string(),
            ),
            (
                "Persistent intermediate hits".into(),
                snapshot
                    .statistics
                    .persistent_cache
                    .intermediate_hits
                    .to_string(),
            ),
            (
                "Persistent Route hits".into(),
                snapshot.statistics.persistent_cache.route_hits.to_string(),
            ),
            (
                "Persistent Apply hits".into(),
                snapshot.statistics.persistent_cache.apply_hits.to_string(),
            ),
            (
                "Persistent misses".into(),
                snapshot.statistics.persistent_cache.misses.to_string(),
            ),
            (
                "Request RAM/work hits".into(),
                snapshot.statistics.work.cache_hits().to_string(),
            ),
            (
                "Cached native coefficient entries".into(),
                snapshot.cached.entries.to_string(),
            ),
            (
                "Cached coefficient terms (guards included)".into(),
                snapshot.cached.coefficient_terms.to_string(),
            ),
            (
                "Cached coefficient bytes (not RSS/encoding)".into(),
                snapshot.cached.coefficient_bytes.to_string(),
            ),
            (
                "Pinned entries".into(),
                snapshot.pinned.map_or_else(
                    || "unknown (policy does not count pins)".into(),
                    |v| v.entries.to_string(),
                ),
            ),
            (
                "Eviction-eligible entries".into(),
                snapshot
                    .eligible
                    .map_or_else(|| "unknown".into(), |v| v.entries.to_string()),
            ),
            ("Nested factor phases".into(), "unobserved".into()),
        ];
        let degrees = envelope.metrics.max_seen_degrees;
        display.complexity = vec![
            (
                "Max seen canonical A".into(),
                degrees.positive_power_sum.to_string(),
            ),
            (
                "Max seen canonical R".into(),
                degrees.numerator_rank.to_string(),
            ),
            (
                "Max seen canonical D=A-R".into(),
                degrees.power_difference.to_string(),
            ),
            (
                "Physical input rank".into(),
                "caller expression; not inferred from canonical R".into(),
            ),
            (
                "Native output terms (sampled)".into(),
                envelope
                    .complexity
                    .terms
                    .map_or_else(|| "unknown".into(), |v| v.to_string()),
            ),
            (
                "Native coefficient terms (sampled)".into(),
                envelope.complexity.coefficient_terms.map_or_else(
                    || "unknown / not materialized for telemetry".into(),
                    |v| v.to_string(),
                ),
            ),
            (
                "Native coefficient bytes (sampled)".into(),
                envelope
                    .complexity
                    .coefficient_bytes
                    .map_or_else(|| "unknown".into(), |v| v.to_string()),
            ),
            (
                "Complexity scope".into(),
                if envelope.complexity.factorized {
                    "factorized native storage"
                } else {
                    "current native operation"
                }
                .into(),
            ),
        ];
    }
    if !envelope.active_path.is_empty() || matches!(envelope.event, Event::CycleDetected) {
        display.path = envelope
            .active_path
            .iter()
            .take(64)
            .copied()
            .map(work_label)
            .collect();
        display.path.push(format!(
            "{} / {} ancestor nodes retained",
            envelope.active_path.len(),
            envelope.path_total_nodes
        ));
        if let Some(backedge) = envelope.cycle_backedge {
            display
                .path
                .push(format!("Native cycle backedge to ancestor {backedge}"));
        }
    }
    if matches!(
        envelope.event,
        Event::CycleDetected | Event::RequestFinished
    ) {
        pending.events.push_back(format!(
            "{:?} / {:?} / {:?}",
            envelope.event, envelope.activity, envelope.status
        ));
        if pending.events.len() > 64 {
            pending.events.pop_front();
        }
    }
    display.events = pending.events.iter().cloned().collect();
}

fn resource_usage(guard_snapshot: Option<&std::path::Path>) -> super::ResourceUsage {
    #[cfg(target_os = "linux")]
    {
        let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
        let meminfo = std::fs::read_to_string("/proc/meminfo").unwrap_or_default();
        let mut sample = super::ResourceUsage::from_proc(&status, &meminfo);
        if let Some(path) = guard_snapshot {
            sample.guard = GuardBinding { path }.sample();
        }
        sample
    }
    #[cfg(not(target_os = "linux"))]
    {
        let _ = guard_snapshot;
        super::ResourceUsage::default()
    }
}

#[cfg(target_os = "linux")]
#[derive(Clone, Copy)]
struct ProcessIdentity {
    pid: u64,
    start_ticks: u64,
    uid: u32,
}

#[cfg(target_os = "linux")]
impl ProcessIdentity {
    fn read(pid: u64) -> io::Result<Self> {
        use std::os::unix::fs::MetadataExt;
        let path = std::path::PathBuf::from(format!("/proc/{pid}"));
        let stat = std::fs::read_to_string(path.join("stat"))?;
        let fields = stat
            .rsplit_once(')')
            .ok_or_else(|| io::Error::other("invalid process stat"))?
            .1
            .split_whitespace()
            .collect::<Vec<_>>();
        if matches!(fields.first(), Some(&"Z" | &"X")) {
            return Err(io::Error::other("registered process is terminal"));
        }
        let start_ticks = fields
            .get(19)
            .and_then(|value| value.parse().ok())
            .ok_or_else(|| io::Error::other("missing process start ticks"))?;
        Ok(Self {
            pid,
            start_ticks,
            uid: std::fs::metadata(path)?.uid(),
        })
    }
}

/// Read-only attachment to the existing external guard, never a policy executor.
#[cfg(target_os = "linux")]
struct GuardBinding<'a> {
    path: &'a std::path::Path,
}

#[cfg(target_os = "linux")]
impl GuardBinding<'_> {
    fn document(
        path: &std::path::Path,
    ) -> io::Result<(serde_json::Value, Vec<u8>, std::time::SystemTime)> {
        use std::{io::Read, os::unix::fs::OpenOptionsExt};
        let file = std::fs::OpenOptions::new()
            .read(true)
            .custom_flags(rustix::fs::OFlags::NONBLOCK.bits() as i32)
            .open(path)?;
        let metadata = file.metadata()?;
        if !metadata.is_file() || metadata.len() > 65536 {
            return Err(io::Error::other(
                "guard input must be a bounded regular file",
            ));
        }
        let modified = metadata.modified()?;
        let mut bytes = Vec::new();
        file.take(65537).read_to_end(&mut bytes)?;
        if bytes.len() > 65536 {
            return Err(io::Error::other("guard input grew beyond bound"));
        }
        let value = serde_json::from_slice(&bytes).map_err(io::Error::other)?;
        Ok((value, bytes, modified))
    }

    fn sample(&self) -> super::ResourceGuardView {
        let result = (|| {
            use sha2::{Digest, Sha256};
            if !self.path.is_absolute() {
                return Err(io::Error::other("guard path must be absolute"));
            }
            let parent = self
                .path
                .parent()
                .ok_or_else(|| io::Error::other("missing guard directory"))?;
            let (heartbeat, _, modified) = Self::document(self.path)?;
            let (registration, _, _) = Self::document(&parent.join("registration.json"))?;
            let (config, bytes, _) = Self::document(&parent.join("guard-config.json"))?;
            let pid = registration["guard"]["pid"]
                .as_u64()
                .ok_or_else(|| io::Error::other("missing guard PID"))?;
            let guard_proc = std::path::PathBuf::from(format!("/proc/{pid}"));
            let cwd = std::fs::read_link(guard_proc.join("cwd"))?;
            let cgroup = std::fs::read_to_string(guard_proc.join("cgroup"))?;
            if config["output_directory"]
                .as_str()
                .map(std::path::Path::new)
                != Some(parent)
                || cwd != parent
                || registration["guard"]["cwd"]
                    .as_str()
                    .map(std::path::Path::new)
                    != Some(cwd.as_path())
                || registration["guard"]["cgroup"].as_str() != Some(cgroup.trim())
            {
                return Err(io::Error::other(
                    "guard workspace or cgroup binding differs",
                ));
            }
            let identities = (
                ProcessIdentity::read(pid)?,
                ProcessIdentity::read(u64::from(std::process::id()))?,
            );
            let age = modified.elapsed().map_err(io::Error::other)?;
            self.project(
                (&config, &format!("{:x}", Sha256::digest(bytes))),
                &registration,
                &heartbeat,
                identities,
                age,
            )
        })();
        result.unwrap_or_else(|error: io::Error| super::ResourceGuardView {
            source: Some(self.path.to_string_lossy().chars().take(512).collect()),
            detail: Some(error.to_string().chars().take(256).collect()),
            ..Default::default()
        })
    }

    fn project(
        &self,
        config: (&serde_json::Value, &str),
        registration: &serde_json::Value,
        heartbeat: &serde_json::Value,
        identities: (ProcessIdentity, ProcessIdentity),
        age: Duration,
    ) -> io::Result<super::ResourceGuardView> {
        let (config, digest) = config;
        let (guard, current) = identities;
        let interval = config["sample_seconds"]
            .as_u64()
            .filter(|value| *value > 0)
            .ok_or_else(|| io::Error::other("invalid guard interval"))?;
        let threshold = config["threshold_bytes"]
            .as_u64()
            .filter(|value| *value > 0)
            .ok_or_else(|| io::Error::other("invalid guard floor"))?;
        if config["schema"] != "vakint.physical-memory-floor-config.v1"
            || config["activation_approved"] != true
            || registration["schema"] != "vakint.physical-memory-floor-registration.v1"
            || heartbeat["schema"] != "vakint.physical-memavailable-guard.v1"
            || registration["config_sha256"] != digest
            || registration["guard"]["pid"] != guard.pid
            || heartbeat["guard_pid"] != guard.pid
            || registration["guard"]["start_ticks"] != guard.start_ticks
            || config["uid"] != guard.uid
            || guard.uid != current.uid
            || registration["guard"]["uids"]
                != serde_json::json!([guard.uid, guard.uid, guard.uid, guard.uid])
            || registration["threshold_bytes"] != threshold
            || heartbeat["threshold_bytes"] != threshold
            || registration["sample_seconds"] != interval
            || heartbeat["sample_seconds"] != interval
            || heartbeat["latched_abort"].as_bool().is_none()
            || heartbeat.get("sample_error").is_none()
        {
            return Err(io::Error::other(
                "guard config, registration, process or snapshot binding differs",
            ));
        }
        let slots = config["slots"]
            .as_array()
            .filter(|slots| slots.len() <= 16)
            .ok_or_else(|| io::Error::other("missing or excessive registered guard slots"))?;
        let mut report = super::ResourceGuardView {
            threshold_bytes: Some(threshold),
            sample_seconds: Some(interval),
            report_age_seconds: Some(age.as_secs_f64()),
            reported_latched_abort: heartbeat["latched_abort"].as_bool(),
            guard_pid: Some(guard.pid),
            source: Some(self.path.to_string_lossy().chars().take(512).collect()),
            detail: Some("current process has no verified registered coverage".into()),
            ..Default::default()
        };
        for slot in slots {
            let name = slot["name"]
                .as_str()
                .ok_or_else(|| io::Error::other("missing guard slot name"))?;
            let observed = &heartbeat["slots"][name];
            let state = observed["state"]
                .as_str()
                .ok_or_else(|| io::Error::other("missing guard slot state"))?;
            let target = observed
                .get("target")
                .or_else(|| observed["event"].get("target"));
            report.registered_slots.push(format!(
                "{}: {}{}",
                name.chars().take(64).collect::<String>(),
                state.chars().take(64).collect::<String>(),
                target.map_or_else(String::new, |target| format!(
                    " PID {}/{}",
                    target["pid"], target["start_ticks"]
                ))
            ));
            if registration["slots"].get(name).is_some()
                && slot["uid"] == current.uid
                && target.is_some_and(|target| {
                    target["pid"] == current.pid && target["start_ticks"] == current.start_ticks
                })
                && (slot["target"].is_null() || target == Some(&slot["target"]))
            {
                report.registered_worker = Some(
                    serde_json::json!({"slot":name,"pid":current.pid,"start_ticks":current.start_ticks}),
                );
                report.current_process_covered = state == "ARMED_EXACT_PIDFD";
                report.status = if report.current_process_covered {
                    "active"
                } else {
                    "fault"
                }
                .into();
                report.detail = Some(format!("registered slot reports {state}"));
            }
        }
        if heartbeat["sample_error"] != serde_json::Value::Null
            || !heartbeat["MemAvailable_bytes"].is_u64()
        {
            report.status = "fault".into();
            report.current_process_covered = false;
            report.detail = Some("guard reports a missing or failed OS sample".into());
        } else if report.reported_latched_abort == Some(true) {
            report.status = "latched_abort".into();
            report.current_process_covered = false;
        }
        if age.as_secs_f64() > 3.0 * interval as f64 {
            report.status = "unknown".into();
            report.current_process_covered = false;
            report.detail = Some("stale guard report; no active coverage claim".into());
        }
        Ok(report)
    }
}

#[cfg(all(test, unix, not(target_arch = "wasm32")))]
mod tests {
    use super::*;

    #[test]
    fn caller_atom_sizes_read_existing_representation_and_detachment_stops_delivery() {
        let expression = crate::vakint_parse!("(size_x+size_y)^9+size_z").unwrap();
        let mut records = Vec::new();
        {
            let mut callback = |envelope: &Envelope<'_>| {
                records.push((
                    envelope.activity,
                    envelope.caller_phase.unwrap().to_owned(),
                    envelope.complexity,
                ));
                Ok(Control::Detach)
            };
            let mut observer = Observer::new(Default::default(), &mut callback);
            MonitorContext::native_phase(None, "disabled", Some(expression.as_view()));
            MonitorContext::native_phase(
                Some(&mut observer),
                "factored_input",
                Some(expression.as_view()),
            );
            assert_eq!(observer.status(), CandidateRoutedObserverStatus::Detached);
            MonitorContext::native_phase(
                Some(&mut observer),
                "detached",
                Some(expression.as_view()),
            );
        }
        assert_eq!(records.len(), 1);
        let (activity, phase, size) = &records[0];
        assert_eq!(*activity, Activity::Caller);
        assert_eq!(phase, "factored_input");
        assert_eq!(size.terms, Some(expression.nterms()));
        assert_eq!(
            size.encoded_expression_bytes,
            Some(expression.as_view().get_byte_size())
        );
        assert_eq!(size.coefficient_bytes, None);
        assert_eq!(size.coefficient_terms, None);
    }

    #[test]
    fn stored_state_age_is_not_refreshed_by_forwarding_or_caller_atom_samples() {
        use rustred::solver::CandidateRoutedStoredState;
        let mut pending = Pending::default();
        let mut envelope = Envelope::new(Event::ActivityStarted, "native-family");
        envelope.sequence = 2;
        envelope.elapsed = Duration::from_secs(3);
        envelope.stored_state = Some(CandidateRoutedStoredState {
            unstarted_roots: Some(2),
            deferred_rule_maps: Some(3),
            deferred_child_integral_occurrences: Some(7),
            collected_output_integrals: Some(1),
        });
        update_display(&mut pending, &envelope, false, false);
        assert_eq!(
            pending.display.stored_state.as_ref().unwrap().sampled_at,
            Duration::from_secs(3)
        );
        assert!(pending.display.quantitative_sample_at.is_none());

        envelope.activity = Activity::TerminalSubstitution;
        envelope.caller_phase = Some("terminal_substitution");
        envelope.elapsed = Duration::from_secs(8);
        envelope.native_summary_age = Some(Duration::from_secs(5));
        update_display(&mut pending, &envelope, true, false);
        assert_eq!(
            pending.display.stored_state.as_ref().unwrap().sampled_at,
            Duration::from_secs(3)
        );

        envelope.activity = Activity::Caller;
        envelope.caller_phase = Some("retained_output");
        envelope.elapsed = Duration::from_secs(11);
        envelope.complexity.terms = Some(4);
        envelope.complexity.encoded_expression_bytes = Some(217);
        update_display(&mut pending, &envelope, true, false);
        assert_eq!(
            pending
                .display
                .caller_expression
                .as_ref()
                .unwrap()
                .sampled_at,
            Duration::from_secs(11)
        );
        assert_eq!(
            pending.display.stored_state.as_ref().unwrap().sampled_at,
            Duration::from_secs(3)
        );

        // The saved application clears its native summary before preflight.
        // A failure here never produces a first core snapshot.
        envelope = Envelope::new(Event::ActivityStarted, "next-family");
        envelope.caller_phase = Some("entry_admission");
        envelope.elapsed = Duration::from_secs(12);
        update_display(&mut pending, &envelope, true, false);
        assert!(pending.display.stored_state.is_none());
        assert!(pending.display.native_observed_at.is_none());
        assert!(pending.display.completed_roots.is_none());
        assert!(pending.display.quantitative_sample_at.is_none());
        envelope.event = Event::RequestFinished;
        envelope.activity = Activity::Caller;
        envelope.status = Status::Failed;
        envelope.request_finished = true;
        envelope.elapsed = Duration::from_secs(13);
        update_display(&mut pending, &envelope, true, false);
        assert!(pending.display.stored_state.is_none());
        assert_eq!(pending.display.status, MonitorStatus::Failed);
        assert_eq!(
            pending
                .display
                .caller_expression
                .as_ref()
                .unwrap()
                .encoded_bytes,
            217
        );
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn guard_attachment_requires_bound_live_identity_freshness_and_actual_slot_coverage() {
        let binding = GuardBinding {
            path: std::path::Path::new("/workspace/guard-fixture/heartbeat.json"),
        };
        let config = serde_json::json!({"schema":"vakint.physical-memory-floor-config.v1",
            "activation_approved":true,"uid":7,"sample_seconds":1,"threshold_bytes":50000000000_u64,
            "slots":[{"name":"q0","uid":7,"target":{"pid":33,"start_ticks":44}},
                {"name":"future","uid":7,"target":null}]});
        let registration = serde_json::json!({"schema":"vakint.physical-memory-floor-registration.v1",
            "config_sha256":"fixture-digest","threshold_bytes":50000000000_u64,"sample_seconds":1,
            "guard":{"pid":11,"start_ticks":22,"uids":[7,7,7,7]},"slots":{"q0":{},"future":{}}});
        let heartbeat = serde_json::json!({"schema":"vakint.physical-memavailable-guard.v1",
            "guard_pid":11,"threshold_bytes":50000000000_u64,"sample_seconds":1,"latched_abort":false,
            "sample_error":null,"MemAvailable_bytes":60000000000_u64,"slots":{
                "q0":{"state":"ARMED_EXACT_PIDFD","target":{"pid":33,"start_ticks":44}},
                "future":{"state":"WAITING_AUTHORIZED_FUTURE_LAUNCH"}}});
        let identities = (
            ProcessIdentity {
                pid: 11,
                start_ticks: 22,
                uid: 7,
            },
            ProcessIdentity {
                pid: 33,
                start_ticks: 44,
                uid: 7,
            },
        );
        let observe = |config: &serde_json::Value,
                       registration: &serde_json::Value,
                       heartbeat: &serde_json::Value,
                       identities,
                       age| {
            binding.project(
                (config, "fixture-digest"),
                registration,
                heartbeat,
                identities,
                age,
            )
        };
        let report = observe(
            &config,
            &registration,
            &heartbeat,
            identities,
            Duration::ZERO,
        )
        .unwrap();
        assert_eq!(report.status, "active");
        assert!(report.current_process_covered);
        assert_eq!(report.threshold_bytes, Some(50000000000));
        assert_eq!(report.registered_worker.as_ref().unwrap()["slot"], "q0");
        let other_worker = (
            identities.0,
            ProcessIdentity {
                pid: 34,
                ..identities.1
            },
        );
        let unregistered = observe(
            &config,
            &registration,
            &heartbeat,
            other_worker,
            Duration::ZERO,
        )
        .unwrap();
        assert_eq!(unregistered.status, "unknown");
        assert!(!unregistered.current_process_covered);
        let stale = observe(
            &config,
            &registration,
            &heartbeat,
            identities,
            Duration::from_secs(4),
        )
        .unwrap();
        assert_eq!(stale.status, "unknown");
        assert!(!stale.current_process_covered);
        for (field, value) in [
            ("start_ticks", serde_json::json!(23)),
            ("uids", serde_json::json!([8, 8, 8, 8])),
        ] {
            let mut wrong = registration.clone();
            wrong["guard"][field] = value;
            assert!(observe(&config, &wrong, &heartbeat, identities, Duration::ZERO).is_err());
        }
        let mut wrong_digest = registration.clone();
        wrong_digest["config_sha256"] = "other".into();
        assert!(
            observe(
                &config,
                &wrong_digest,
                &heartbeat,
                identities,
                Duration::ZERO
            )
            .is_err()
        );
        let mut latched = heartbeat.clone();
        latched["latched_abort"] = true.into();
        let report = observe(&config, &registration, &latched, identities, Duration::ZERO).unwrap();
        assert_eq!(report.status, "latched_abort");
        assert!(!report.current_process_covered);
        let mut fault = heartbeat.clone();
        fault["sample_error"] = "synthetic OS error".into();
        assert_eq!(
            observe(&config, &registration, &fault, identities, Duration::ZERO)
                .unwrap()
                .status,
            "fault"
        );
        assert_eq!(super::super::ResourceGuardView::default().status, "unknown");
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn guard_attachment_rejects_fifo_and_oversized_inputs_without_a_writer() {
        let directory = tempfile::tempdir().unwrap();
        let fifo = directory.path().join("heartbeat.json");
        rustix::fs::mkfifoat(
            rustix::fs::CWD,
            &fifo,
            rustix::fs::Mode::RUSR | rustix::fs::Mode::WUSR,
        )
        .unwrap();
        let started = Instant::now();
        assert!(GuardBinding::document(&fifo).is_err());
        assert!(started.elapsed() < Duration::from_secs(1));
        let large = directory.path().join("large.json");
        std::fs::write(&large, vec![b' '; 65537]).unwrap();
        assert!(GuardBinding::document(&large).is_err());
        let report = GuardBinding { path: &fifo }.sample();
        assert_eq!(report.status, "unknown");
        assert!(!report.current_process_covered);
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn real_proc_memory_refreshes_during_synthetic_native_silence_and_keeps_cycle_failure() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("resources.json");
        let options = RustRedMonitorOptions {
            progress_snapshot: Some(path.clone()),
            interval: Duration::from_millis(100),
            ..Default::default()
        };
        let failure = MonitorContext::run::<()>(&options, |mut observer| {
            // Real OS sampling/file publication; synthetic control events only.
            // No saved application, coefficient, cache or reduction is loaded.
            MonitorContext::native_phase(observer.as_deref_mut(), "synthetic_busy_phase", None);
            let deadline = Instant::now() + Duration::from_secs(5);
            let first = loop {
                if let Ok(bytes) = std::fs::read(&path) {
                    let document: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
                    if document["monitor"]["resources"]["rss_bytes"].is_u64() {
                        break document;
                    }
                }
                assert!(Instant::now() < deadline, "resource publication timed out");
                thread::sleep(Duration::from_millis(5));
            };
            let refreshed = loop {
                let document: serde_json::Value =
                    serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
                if document["monitor"]["published_unix_seconds"]
                    .as_f64()
                    .unwrap()
                    > first["monitor"]["published_unix_seconds"].as_f64().unwrap()
                    && document["monitor"]["native_document_age_seconds"]
                        .as_f64()
                        .unwrap()
                        >= 0.1
                {
                    break document;
                }
                assert!(
                    Instant::now() < deadline,
                    "independent OS refresh timed out"
                );
                thread::sleep(Duration::from_millis(5));
            };
            assert!(
                refreshed["monitor"]["resources"]["sample_age_seconds"]
                    .as_f64()
                    .unwrap()
                    < 0.1
            );
            assert!(refreshed["monitor"]["last_cycle"].is_null());
            let mut cycle = Envelope::new(Event::CycleDetected, "synthetic-resource-fixture");
            cycle.activity = Activity::CoefficientCombination;
            cycle.status = Status::Failed;
            observer.as_mut().unwrap().emit(cycle);
            Err(VakintError::EvaluationError(
                "synthetic cycle fixture".into(),
            ))
        })
        .unwrap_err();
        assert!(matches!(failure, VakintError::EvaluationError(message)
            if message == "synthetic cycle fixture"));
        let document: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
        let resources = &document["monitor"]["resources"];
        for field in [
            "rss_bytes",
            "peak_rss_bytes",
            "mem_total_bytes",
            "mem_available_bytes",
            "mem_free_bytes",
        ] {
            assert!(resources[field].is_u64(), "{field}");
        }
        assert!(
            resources["peak_rss_bytes"].as_u64().unwrap()
                >= resources["rss_bytes"].as_u64().unwrap()
        );
        assert!(
            resources["source"]
                .as_str()
                .unwrap()
                .contains("/proc/meminfo")
        );
        assert_eq!(document["monitor"]["status"], "Failed");
        assert_eq!(
            document["monitor"]["global_remaining_work_certificate"],
            false
        );
        assert!(document["monitor"]["last_cycle"].is_object());
        let view: super::super::MonitorView =
            serde_json::from_value(document["monitor"]["view"].clone()).unwrap();
        assert!(
            view.estimates
                .iter()
                .any(|(label, value)| label == "Discovered-frontier ETA"
                    && value.contains("unknown"))
        );
        assert!(
            view.events
                .iter()
                .any(|event| event.contains("CycleDetected"))
        );
        #[cfg(feature = "progress-tui")]
        {
            let backend = ratatui::backend::TestBackend::new(220, 80);
            let mut terminal = ratatui::Terminal::new(backend).unwrap();
            terminal
                .draw(|frame| super::super::tui::draw(frame, &view, false))
                .unwrap();
            let rendered = terminal
                .backend()
                .buffer()
                .content
                .iter()
                .map(|cell| cell.symbol())
                .collect::<String>();
            for label in [
                "VmRSS / VmHWM",
                "MemAvailable / MemFree",
                "usable / unused",
                "MemTotal",
                "sample age",
                "Caller Atom bytes / stored summands",
                "CycleDetected",
                "Failed",
            ] {
                assert!(rendered.contains(label), "missing {label}");
            }
        }
    }

    #[test]
    fn dirty_activity_respects_file_cadence_and_final_document_is_drained() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("latest.json");
        let options = RustRedMonitorOptions {
            progress_snapshot: Some(path.clone()),
            interval: Duration::from_secs(60),
            ..Default::default()
        };
        let result = MonitorContext::run(&options, |mut observer| {
            let deadline = Instant::now() + Duration::from_secs(5);
            while !path.exists() {
                assert!(
                    Instant::now() < deadline,
                    "first snapshot was not published"
                );
                thread::sleep(Duration::from_millis(1));
            }
            let first = std::fs::read(&path).unwrap();
            for _ in 0..100 {
                MonitorContext::native_phase(observer.as_deref_mut(), "busy_phase", None);
            }
            thread::sleep(Duration::from_millis(100));
            assert_eq!(
                std::fs::read(&path).unwrap(),
                first,
                "ordinary dirty notices bypassed the interval"
            );
            Ok(())
        });
        result.unwrap();
        let final_document: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
        assert_eq!(final_document["status"], "succeeded");
        assert_eq!(final_document["request_finished"], true);
        assert_eq!(final_document["caller_phase"], "busy_phase");
        assert_eq!(
            final_document["monitor"]["view"]["phase"],
            serde_json::Value::Null
        );
        assert!(final_document["monitor"]["estimates"].is_array());
    }

    #[test]
    fn cheap_and_forwarded_notices_retain_the_coherent_quantitative_sample() {
        use rustred::solver::{
            CandidateRoutedApplicationPhase, CandidateRoutedApplicationSnapshot,
            CandidateRoutedFrontierSample,
        };
        let mut pending = Pending::default();
        let mut envelope = Envelope::new(Event::Snapshot, "fake-only");
        envelope.sequence = 1;
        envelope.elapsed = Duration::from_secs(10);
        envelope.snapshot = Some(CandidateRoutedApplicationSnapshot {
            sequence: 1,
            elapsed: Duration::from_secs(10),
            phase: CandidateRoutedApplicationPhase::Expansion,
            status: Status::Running,
            statistics: Default::default(),
            unique_nodes: 20,
            active_nodes: 1,
            pending_frames: 12,
            cached: Default::default(),
            pinned: None,
            eligible: None,
        });
        envelope.frontier = Some(CandidateRoutedFrontierSample {
            unique_nodes: 12,
            ..Default::default()
        });
        envelope.metrics.completed_unique_operational_nodes = 5;
        update_display(&mut pending, &envelope, true, true);
        let epoch = pending.display.native_epoch;
        envelope.event = Event::ActivityStarted;
        envelope.sequence = 2;
        envelope.elapsed = Duration::from_secs(30);
        envelope.native_summary_age = Some(Duration::ZERO);
        envelope.frontier = None;
        envelope.metrics.completed_unique_operational_nodes = 999;
        update_display(&mut pending, &envelope, true, false);
        assert_eq!(pending.display.completed_unique, Some(5));
        assert_eq!(pending.display.pending_unique.as_ref().unwrap().count, 12);
        assert_eq!(
            pending.display.quantitative_sample_at,
            Some(Duration::from_secs(10))
        );
        assert_eq!(pending.display.quantitative_sample_sequence, 1);
        // Saved forwarding can be non-Caller but always has a caller phase.
        envelope.activity = Activity::TerminalSubstitution;
        envelope.caller_phase = Some("saved_terminal_forward");
        envelope.frontier = Some(CandidateRoutedFrontierSample {
            unique_nodes: 999,
            ..Default::default()
        });
        envelope.sequence = 3;
        envelope.elapsed = Duration::from_secs(40);
        update_display(&mut pending, &envelope, true, false);
        assert_eq!(pending.display.native_epoch, epoch);
        assert_eq!(pending.display.completed_unique, Some(5));
        assert_eq!(pending.display.pending_unique.as_ref().unwrap().count, 12);
        assert_eq!(
            pending.display.quantitative_sample_at,
            Some(Duration::from_secs(10))
        );
        assert_eq!(pending.display.quantitative_sample_sequence, 1);
    }

    #[test]
    fn fast_core_finish_retains_final_counts_before_outer_caller_finish() {
        use rustred::solver::{
            CandidateRoutedApplicationPhase, CandidateRoutedApplicationSnapshot,
            CandidateRoutedFrontierSample,
        };
        // Only fake envelopes and the bounded in-memory pending queue; no
        // solver, worker thread, terminal, file, sleep or reduction is started.
        let shared = Arc::new(Shared {
            pending: Mutex::new(Pending::default()),
            active: AtomicBool::new(true),
            stopping: AtomicBool::new(false),
            started: Instant::now(),
        });
        let (wake, _receiver) = mpsc::sync_channel(1);
        let handle = MonitorHandle {
            shared: Arc::clone(&shared),
            wake,
            interval: Duration::from_secs(1),
        };
        let mut snapshot = CandidateRoutedApplicationSnapshot {
            sequence: 1,
            elapsed: Duration::ZERO,
            phase: CandidateRoutedApplicationPhase::Admission,
            status: Status::Running,
            statistics: Default::default(),
            unique_nodes: 0,
            active_nodes: 0,
            pending_frames: 0,
            cached: Default::default(),
            pinned: None,
            eligible: None,
        };
        let mut envelope = Envelope::new(Event::Snapshot, "fake-fast-finish");
        envelope.sequence = 1;
        envelope.snapshot = Some(snapshot);
        envelope.frontier = Some(CandidateRoutedFrontierSample::default());
        envelope.metrics.requested_roots = Some(1);
        assert!(matches!(handle.observe(&envelope), Ok(Control::Continue)));
        {
            // Simulate the worker draining only the initial publication flag.
            let mut pending = shared.pending();
            pending.force_publish = false;
            pending.document = None;
        }
        snapshot.sequence = 2;
        snapshot.elapsed = Duration::from_millis(100);
        snapshot.phase = CandidateRoutedApplicationPhase::RootCollection;
        snapshot.status = Status::Succeeded;
        snapshot.unique_nodes = 3;
        snapshot.cached.entries = 3;
        snapshot.statistics.persistent_cache.requested_root_hits = 1;
        envelope.event = Event::CoreFinished;
        envelope.activity = Activity::RootCollection;
        envelope.sequence = 2;
        envelope.elapsed = Duration::from_millis(100);
        envelope.status = Status::Succeeded;
        envelope.snapshot = Some(snapshot);
        envelope.metrics.completed_unique_operational_nodes = 3;
        envelope.metrics.completion_attempts = 3;
        envelope.metrics.completed_roots = 1;
        envelope.frontier = Some(CandidateRoutedFrontierSample {
            exact: true,
            ..Default::default()
        });
        assert!(matches!(handle.observe(&envelope), Ok(Control::Continue)));
        {
            let pending = shared.pending();
            assert_eq!(pending.last_document_at, Some(Duration::from_millis(100)));
            assert!(pending.document.is_some());
            assert!(
                !pending.force_publish,
                "core completion must not force outer success publication"
            );
            assert_eq!(pending.display.status, MonitorStatus::Running);
            assert!(!pending.display.request_finished);
            assert_eq!(pending.display.completed_unique, Some(3));
        }
        envelope.event = Event::RequestFinished;
        envelope.activity = Activity::Caller;
        envelope.caller_phase = Some("request_done");
        envelope.sequence = 3;
        envelope.elapsed = Duration::from_millis(200);
        envelope.request_finished = true;
        envelope.frontier = None;
        assert!(matches!(handle.observe(&envelope), Ok(Control::Continue)));
        let pending = shared.pending();
        assert_eq!(pending.display.native_epoch, 1);
        assert_eq!(pending.display.completed_unique, Some(3));
        assert_eq!(pending.display.pending_unique.as_ref().unwrap().count, 0);
        assert_eq!(
            pending.display.quantitative_sample_at,
            Some(Duration::from_millis(100))
        );
        assert_eq!(pending.display.quantitative_sample_sequence, 2);
        assert_eq!(
            pending
                .display
                .cache
                .iter()
                .find(|(label, _)| label == "Persistent requested-root hits")
                .unwrap()
                .1,
            "1"
        );
        assert_eq!(
            pending
                .display
                .cache
                .iter()
                .find(|(label, _)| label == "Cached native coefficient entries")
                .unwrap()
                .1,
            "3"
        );
    }
    #[test]
    fn activity_refreshes_current_native_key_without_fabricating_a_census() {
        use rustred::{family::IntegralKey, solver::CandidateRoutedKeyDegrees};
        let first = IntegralKey::try_new([-2, 4, 1]).unwrap();
        let second = IntegralKey::try_new([-1, 2, 8]).unwrap();
        let mut pending = Pending::default();
        let mut event = Envelope::new(Event::ActivityStarted, "fake-current-only");
        event.current = Some(CandidateRoutedWorkRef::Route(&first));
        event.degrees = Some(CandidateRoutedKeyDegrees::from_key(&first));
        event.metrics.max_seen_degrees.numerator_rank = 99;
        update_display(&mut pending, &event, false, false);
        assert_eq!(pending.display.current_degrees, Some((2, 5, 3)));
        assert!(
            pending
                .display
                .current_key
                .as_ref()
                .unwrap()
                .contains("-2,4,1")
        );
        assert!(pending.display.quantitative_sample_at.is_none());
        event.current = Some(CandidateRoutedWorkRef::Route(&second));
        event.degrees = Some(CandidateRoutedKeyDegrees::from_key(&second));
        update_display(&mut pending, &event, false, false);
        assert_eq!(pending.display.current_degrees, Some((1, 10, 9)));
        assert!(
            pending
                .display
                .current_key
                .as_ref()
                .unwrap()
                .contains("-1,2,8")
        );
        event.activity = Activity::Caller;
        event.caller_phase = Some("numerical_materialization");
        update_display(&mut pending, &event, false, false);
        assert!(pending.display.current_key.is_none());
        assert!(pending.display.current_degrees.is_none());
        assert!(pending.display.quantitative_sample_at.is_none());
    }
}

#[cfg(all(test, not(feature = "experimental-rustred")))]
mod disabled_benchmark_tests {
    use super::*;

    #[test]
    fn default_off_invokes_the_unobserved_operation_once() {
        let mut calls = 0;
        let value = MonitorContext::run(&RustRedMonitorOptions::default(), |observation| {
            calls += 1;
            assert!(observation.is_none());
            Ok(17)
        })
        .unwrap();
        assert_eq!(value, 17);
        assert_eq!(calls, 1);
    }
}

#[cfg(all(test, feature = "experimental-rustred"))]
mod benchmark_tests {
    use super::*;
    use rustred::{
        family::IntegralKey,
        solver::{
            CandidateRoutedApplicationPhase, CandidateRoutedApplicationSnapshot,
            CandidateRoutedCompletionSource, CandidateRoutedObserverFault,
        },
    };

    #[test]
    fn diagnostic_observation_requires_explicit_opt_in_with_broad_debug_logger() {
        const CHILD: &str = "VAKINT_TEST_EXPLICIT_CACHE_HIT_OBSERVATION";
        if std::env::var_os(CHILD).is_none() {
            // log::set_logger is process-global and set-once. Isolate this
            // regression from other library tests that install test_log first.
            let name = concat!(
                module_path!(),
                "::diagnostic_observation_requires_explicit_opt_in_with_broad_debug_logger"
            );
            let name = name.split_once("::").unwrap().1;
            let output = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", name, "--nocapture", "--test-threads=1"])
                .env(CHILD, "1")
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "isolated logger regression failed: {}{}",
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
            assert!(String::from_utf8_lossy(&output.stdout).contains("1 passed; 0 failed"));
            return;
        }
        struct BroadDebug;
        impl log::Log for BroadDebug {
            fn enabled(&self, _: &log::Metadata<'_>) -> bool {
                true
            }
            fn log(&self, _: &log::Record<'_>) {}
            fn flush(&self) {}
        }
        static LOGGER: BroadDebug = BroadDebug;
        log::set_logger(&LOGGER).unwrap();
        log::set_max_level(log::LevelFilter::Debug);
        assert!(log::log_enabled!(target: "vakint::benchmark_cache_hits", log::Level::Debug));
        let normal = RustRedMonitorOptions::default();
        assert!(!normal.enabled());
        assert!(!normal.cache_hit_diagnostic);
        let mut calls = 0;
        let value = MonitorContext::run(&normal, |observer| {
            calls += 1;
            assert!(observer.is_none());
            Ok(17)
        })
        .unwrap();
        assert_eq!((value, calls), (17, 1));
        let diagnostic = RustRedMonitorOptions {
            cache_hit_diagnostic: true,
            ..Default::default()
        };
        assert!(!diagnostic.enabled());
        assert!(diagnostic.progress_snapshot.is_none());
        let mut calls = 0;
        let value = MonitorContext::run(&diagnostic, |observer| {
            calls += 1;
            assert!(observer.is_some());
            Ok(23)
        })
        .unwrap();
        assert_eq!((value, calls), (23, 1));
        let failure = MonitorContext::run::<()>(&diagnostic, |observer| {
            assert!(observer.is_some());
            Err(VakintError::EvaluationError(
                "metadata-only sentinel".into(),
            ))
        })
        .unwrap_err();
        assert!(
            matches!(failure, VakintError::EvaluationError(message) if message == "metadata-only sentinel")
        );
        let invalid = RustRedMonitorOptions {
            cache_hit_diagnostic: true,
            display: RustRedMonitorDisplay::Json,
            progress_snapshot: Some("must-not-be-created.json".into()),
            ..Default::default()
        };
        assert_eq!(
            invalid.validate().unwrap_err().kind(),
            std::io::ErrorKind::InvalidInput
        );
        let error = MonitorContext::run::<()>(&invalid, |_| {
            panic!("invalid diagnostic must not execute operation")
        })
        .unwrap_err();
        assert!(
            matches!(error, VakintError::EvaluationError(message) if message.contains("cannot be combined"))
        );
    }

    fn snapshot(sequence: u64) -> CandidateRoutedApplicationSnapshot {
        let statistics = rustred::solver::CandidateRoutedApplicationStatistics {
            persistent_cache: rustred::solver::CandidateRoutedPersistentCacheStatistics {
                lookups: 7,
                route_hits: 2,
                apply_hits: 3,
                requested_root_hits: 1,
                intermediate_hits: 4,
                misses: 2,
                ..Default::default()
            },
            ..Default::default()
        };
        CandidateRoutedApplicationSnapshot {
            sequence,
            elapsed: Duration::ZERO,
            phase: CandidateRoutedApplicationPhase::Expansion,
            status: Status::Running,
            statistics,
            unique_nodes: 0,
            active_nodes: 0,
            pending_frames: 0,
            cached: Default::default(),
            pinned: None,
            eligible: None,
        }
    }

    #[test]
    fn benchmark_events_keep_full_work_identity_and_native_persistent_counters() {
        // Deliberately exceed the UI preview's fifteen slots: the proof stream
        // must retain the complete physical key and Apply owner, not a card.
        let powers = (-2..16).collect::<Vec<i64>>();
        let key = IntegralKey::try_new(powers.clone()).unwrap();
        let owner = powers.iter().map(|power| *power > 0).collect::<Vec<_>>();
        let mut event = Envelope::new(
            Event::NodeCompleted(CandidateRoutedCompletionSource::BackingHit),
            "fake-metadata-only",
        );
        event.sequence = 23;
        event.snapshot = Some(snapshot(5));
        event.current = Some(CandidateRoutedWorkRef::Apply {
            owner_sector: &owner,
            target: &key,
        });
        event.root = Some(CandidateRoutedWorkRef::Route(&key));
        event.metrics.requested_roots = Some(2);
        event.metrics.completed_roots = 1;
        let document = MonitorContext::benchmark_event_document(&event).unwrap();
        assert_eq!(document["schema"], "vakint.cache-hit-events.v1");
        assert_eq!(document["event"], "backing_hit");
        assert_eq!(document["family_fingerprint"], "fake-metadata-only");
        assert_eq!(document["sequence"], 23);
        assert_eq!(document["work"]["stage"], "apply");
        assert_eq!(document["work"]["powers"], serde_json::json!(powers));
        assert_eq!(document["work"]["owner_sector"], serde_json::json!(owner));
        assert!(document["rootwork"].is_null());
        assert_eq!(document["requested_roots"], 2);
        assert_eq!(document["completed_roots"], 1);
        assert_eq!(document["persistent"]["lookups"], 7);
        assert_eq!(document["persistent"]["route_hits"], 2);
        assert_eq!(document["persistent"]["apply_hits"], 3);
        assert_eq!(document["persistent"]["requested_root_hits"], 1);
        assert_eq!(document["persistent"]["intermediate_hits"], 4);
        assert_eq!(document["persistent"]["misses"], 2);
        assert_eq!(document["status"], "running");

        event.event = Event::RootCollected;
        event.current = Some(CandidateRoutedWorkRef::Route(&key));
        let document = MonitorContext::benchmark_event_document(&event).unwrap();
        assert_eq!(document["event"], "root_collected");
        assert_eq!(document["rootwork"]["stage"], "route");
        assert_eq!(document["rootwork"]["powers"], serde_json::json!(powers));
        assert!(document["rootwork"]["owner_sector"].is_null());
    }

    #[test]
    fn benchmark_filter_uses_the_native_initial_snapshot_and_only_backing_hits() {
        let mut event = Envelope::new(Event::Snapshot, "fake-metadata-only");
        event.snapshot = Some(snapshot(1));
        assert_eq!(
            MonitorContext::benchmark_event_document(&event).unwrap()["event"],
            "dispatch_start"
        );
        event.snapshot = Some(snapshot(2));
        assert!(MonitorContext::benchmark_event_document(&event).is_none());
        for ignored in [
            Event::ActivityStarted,
            Event::CallerProgress,
            Event::RequestFinished,
            Event::NodeCompleted(CandidateRoutedCompletionSource::Computed),
            Event::NodeCompleted(CandidateRoutedCompletionSource::RamHit),
        ] {
            event.event = ignored;
            assert!(MonitorContext::benchmark_event_document(&event).is_none());
        }
        event.event = Event::CoreFinished;
        event.status = Status::Failed;
        assert_eq!(
            MonitorContext::benchmark_event_document(&event).unwrap()["status"],
            "failed"
        );
        event.caller_phase = Some("forwarded_summary");
        assert!(MonitorContext::benchmark_event_document(&event).is_none());
        event.caller_phase = None;
        event.activity = Activity::Caller;
        assert!(MonitorContext::benchmark_event_document(&event).is_none());
    }

    #[test]
    fn missing_proof_metadata_is_explicit_null_instead_of_a_dropped_event() {
        let mut event = Envelope::new(
            Event::NodeCompleted(CandidateRoutedCompletionSource::BackingHit),
            "fake-metadata-only",
        );
        let document = MonitorContext::benchmark_event_document(&event).unwrap();
        assert_eq!(document["event"], "backing_hit");
        assert!(document["work"].is_null());
        assert!(document["rootwork"].is_null());
        assert!(document["persistent"].is_null());
        assert!(document["requested_roots"].is_null());
        event.event = Event::RootCollected;
        let document = MonitorContext::benchmark_event_document(&event).unwrap();
        assert_eq!(document["event"], "root_collected");
        assert!(document["rootwork"].is_null());
    }

    #[test]
    fn finish_record_preserves_actual_observer_fault_detach_and_operation_failure() {
        for (action, expected) in [
            (Ok(Control::Continue), CandidateRoutedObserverStatus::Active),
            (Ok(Control::Detach), CandidateRoutedObserverStatus::Detached),
            (
                Err(()),
                CandidateRoutedObserverStatus::Faulted(CandidateRoutedObserverFault::CallbackError),
            ),
        ] {
            let status = {
                let mut callback = |_: &Envelope<'_>| action;
                let mut observer = Observer::new(Default::default(), &mut callback);
                observer.emit(Envelope::new(Event::CallerProgress, "fake-metadata-only"));
                observer.status()
            };
            assert_eq!(status, expected);
            let document = MonitorContext::benchmark_finish_document(status, false);
            assert_eq!(document["event"], "observer_finish");
            assert!(document["sequence"].is_null());
            assert_eq!(document["observer_status"], format!("{expected:?}"));
            assert_eq!(document["operation_succeeded"], false);
            assert_eq!(document["interval_ns"], 0);
        }
    }
}
