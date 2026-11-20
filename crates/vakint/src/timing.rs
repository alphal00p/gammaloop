//! Always-on, thread-local, cumulative wall-clock counters that split the
//! time spent inside vakint into the integral REDUCTION and the LOWERING.
//!
//! * `reduce`: every path into an integral reduction (the RustRed p library
//!   and ordered vacuum reductions, the MATAD / FMFT / FORCER / alphaLoop /
//!   pySecDec backends), including the master-coefficient finalisation that
//!   follows (master substitution, Laurent expansion, numerical resolution).
//! * `lower`: the rest of vakint's work on the way there: turning numerators
//!   into integral keys (topology matching, tensor reduction, direct lowering,
//!   block reads, folding).
//!
//! Guards mark scopes of one kind; time is attributed to the INNERMOST active
//! guard, so a reduction scope inside a lowering entry point counts as
//! reduction, and a lowering function called from inside a coarse reduction
//! scope counts as lowering. A guard of the kind already active is a no-op
//! (re-entrancy: nested entry points are not double counted), so for every
//! interval `reduce + lower <= wall`. `reduce_calls` counts the entries into
//! a reduction scope from outside any reduction scope.
//!
//! SUB-STEPS, a second axis with the same guard rules: [`contract`] marks the
//! numerator processing vakint does (the idenso / spenso contraction of the
//! factorised entries and of the hand-over path, a record's Taylor step and q
//! averaging), [`finalize`] the master substitution, ε expansion and assembly
//! after a reduction. A sub-step guard never changes the reduce / lower
//! attribution: `reduce_seconds` and `lower_seconds` mean exactly what they
//! meant before the sub-steps existed (`reduce_seconds` stays inclusive of the
//! finalisation). The sub-step time is also recorded per enclosing kind, so
//! [`Delta::stage`] gives FOUR EXCLUSIVE stages: CONTRACTION, LOWERING (P12's
//! lower without the sub-steps inside it), REDUCTION (P12's reduce without the
//! sub-steps inside it) and FINALISATION; their sum is `<= wall`.
//!
//! RESOURCES (opt-in: [`set_resources`] or `VAKINT_TIMING_RESOURCES=1`, read
//! once; off by default): at every transition also the calling thread's CPU
//! time (the kernel's per-thread on-CPU nanoseconds, `/proc/thread-self/
//! schedstat`), attributed like the wall time, and the process's resident set
//! size (`VmRSS` of `/proc/self/status`), sampled at the stage
//! boundaries: per stage the largest RSS sampled at the end of one of its
//! intervals ([`StageTime::peak_rss_bytes`], cumulative on the thread since
//! the resources were enabled: with one diagram per thread, that diagram's
//! peak at the stage boundaries; the RSS is the PROCESS's, so with several
//! threads it includes the others' memory). Off: no system call, one relaxed
//! atomic load per transition.
//!
//! The counters only advance on the thread that holds the guards: callers read
//! a [`Snapshot`] before and after a unit on the same thread. Callers measure
//! wall time, so work a reduction hands to helper threads is included in the
//! caller's reduction seconds (it waits for them); its CPU time is not (it is
//! the helper threads'). Cost: two clock reads per scope transition, no
//! allocation.

use std::cell::Cell;
use std::marker::PhantomData;
use std::sync::atomic::{AtomicU8, Ordering};
use std::time::Instant;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    None,
    Lower,
    Reduce,
}

impl Kind {
    fn slot(self) -> usize {
        match self {
            Kind::None => 0,
            Kind::Lower => 1,
            Kind::Reduce => 2,
        }
    }
}

/// The sub-step axis.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Step {
    None,
    Contract,
    Finalize,
}

impl Step {
    fn row(self) -> Option<usize> {
        match self {
            Step::None => None,
            Step::Contract => Some(0),
            Step::Finalize => Some(1),
        }
    }
}

/// Per (sub-step, enclosing kind) nanoseconds: `[contract | finalize][none |
/// lower | reduce]`.
type Matrix = [[u64; 3]; 2];

/// One of the four exclusive stages of [`Delta::stage`].
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Stage {
    /// Numerator processing: idenso / spenso contraction called by vakint.
    Contraction,
    /// Lowering onto RustRed's slots (matching, tensor reduction, direct
    /// lowering, folding), without the sub-steps inside it.
    Lowering,
    /// The reduction (RustRed and the other backends), without the
    /// finalisation inside it.
    Reduction,
    /// Master substitution, ε expansion, assembly.
    Finalization,
}

impl Stage {
    pub const ALL: [Stage; 4] = [
        Stage::Contraction,
        Stage::Lowering,
        Stage::Reduction,
        Stage::Finalization,
    ];

    fn index(self) -> usize {
        self as usize
    }

    /// The stage the active (kind, sub-step) belongs to; `None` outside every
    /// scope.
    fn of(kind: Kind, step: Step) -> Option<Stage> {
        match (step, kind) {
            (Step::Contract, _) => Some(Stage::Contraction),
            (Step::Finalize, _) => Some(Stage::Finalization),
            (Step::None, Kind::Lower) => Some(Stage::Lowering),
            (Step::None, Kind::Reduce) => Some(Stage::Reduction),
            (Step::None, Kind::None) => None,
        }
    }
}

struct State {
    reduce_ns: Cell<u64>,
    lower_ns: Cell<u64>,
    calls: Cell<u64>,
    mode: Cell<Kind>,
    since: Cell<Option<Instant>>,
    reduce_depth: Cell<u32>,
    step: Cell<Step>,
    step_ns: Cell<Matrix>,
    // resources (advanced only while enabled)
    cpu_since: Cell<Option<u64>>,
    stage_cpu_ns: Cell<[u64; 4]>,
    cpu_ns: Cell<u64>,
    stage_rss: Cell<[u64; 4]>,
}

thread_local! {
    static STATE: State = const {
        State {
            reduce_ns: Cell::new(0),
            lower_ns: Cell::new(0),
            calls: Cell::new(0),
            mode: Cell::new(Kind::None),
            since: Cell::new(None),
            reduce_depth: Cell::new(0),
            step: Cell::new(Step::None),
            step_ns: Cell::new([[0; 3]; 2]),
            cpu_since: Cell::new(None),
            stage_cpu_ns: Cell::new([0; 4]),
            cpu_ns: Cell::new(0),
            stage_rss: Cell::new([0; 4]),
        }
    };
}

/// 0 = not read yet, 1 = off, 2 = on.
static RESOURCES: AtomicU8 = AtomicU8::new(0);

/// Whether thread CPU time and the RSS samples are recorded.
pub fn resources_enabled() -> bool {
    match RESOURCES.load(Ordering::Relaxed) {
        1 => false,
        2 => true,
        _ => {
            let on = matches!(
                std::env::var("VAKINT_TIMING_RESOURCES")
                    .as_deref()
                    .map(str::trim),
                Ok("1")
            );
            RESOURCES.store(if on { 2 } else { 1 }, Ordering::Relaxed);
            on
        }
    }
}

/// Switch the resource counters on or off for the whole process (overrides
/// `VAKINT_TIMING_RESOURCES`). Counters already accumulated are kept.
pub fn set_resources(on: bool) {
    RESOURCES.store(if on { 2 } else { 1 }, Ordering::Relaxed);
}

/// CPU nanoseconds of the calling thread: the kernel's on-CPU time of the
/// thread (`/proc/thread-self/schedstat`, first field, nanoseconds; the
/// scheduler clock behind `CLOCK_THREAD_CPUTIME_ID`). `None` where it is not
/// readable.
pub fn thread_cpu_ns() -> Option<u64> {
    let text = std::fs::read_to_string("/proc/thread-self/schedstat").ok()?;
    text.split_whitespace().next()?.parse().ok()
}

/// A `kB` field of `/proc/self/status`, in bytes.
fn status_bytes(field: &str) -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find(|l| l.starts_with(field))?;
    let kb: u64 = line[field.len()..]
        .trim()
        .trim_end_matches("kB")
        .trim()
        .parse()
        .ok()?;
    Some(kb.saturating_mul(1024))
}

/// The process's resident set size in bytes (`VmRSS`).
pub fn process_rss_bytes() -> Option<u64> {
    status_bytes("VmRSS:")
}

/// The process's resident high-water mark in bytes (`VmHWM`), process-wide
/// and monotone.
pub fn process_peak_rss_bytes() -> Option<u64> {
    status_bytes("VmHWM:")
}

fn nanos(from: Option<Instant>, to: Instant) -> u64 {
    from.map(|from| to.saturating_duration_since(from).as_nanos() as u64)
        .unwrap_or(0)
}

impl State {
    /// Add the time since the last transition to the active kind's counter
    /// and the active sub-step's; with resources, the thread CPU to the
    /// active stage and the RSS sample to its maximum.
    fn flush(&self, now: Instant) {
        let delta = nanos(self.since.get(), now);
        let mode = self.mode.get();
        match mode {
            Kind::Reduce => self
                .reduce_ns
                .set(self.reduce_ns.get().saturating_add(delta)),
            Kind::Lower => self.lower_ns.set(self.lower_ns.get().saturating_add(delta)),
            Kind::None => {}
        }
        if let Some(row) = self.step.get().row() {
            let mut m = self.step_ns.get();
            m[row][mode.slot()] = m[row][mode.slot()].saturating_add(delta);
            self.step_ns.set(m);
        }
        self.since.set(Some(now));
        if resources_enabled() {
            let cpu = thread_cpu_ns();
            let stage = Stage::of(mode, self.step.get());
            if let (Some(c0), Some(c1)) = (self.cpu_since.get(), cpu) {
                let dc = c1.saturating_sub(c0);
                self.cpu_ns.set(self.cpu_ns.get().saturating_add(dc));
                if let Some(stage) = stage {
                    let mut a = self.stage_cpu_ns.get();
                    a[stage.index()] = a[stage.index()].saturating_add(dc);
                    self.stage_cpu_ns.set(a);
                }
            }
            self.cpu_since.set(cpu);
            if let (Some(stage), Some(rss)) = (stage, process_rss_bytes()) {
                let mut a = self.stage_rss.get();
                a[stage.index()] = a[stage.index()].max(rss);
                self.stage_rss.set(a);
            }
        } else {
            self.cpu_since.set(None);
        }
    }
}

/// A scope of one kind (see the module documentation). Not `Send`: it must
/// be dropped on the thread that created it, in reverse creation order.
#[must_use = "a timing guard measures until it is dropped"]
pub struct Guard {
    /// The kind active before this guard; `None` for a no-op guard.
    previous: Option<Kind>,
    kind: Kind,
    _not_send: PhantomData<*const ()>,
}

impl Guard {
    fn enter(kind: Kind) -> Self {
        let previous = STATE
            .try_with(|state| {
                let previous = state.mode.get();
                if previous == kind {
                    return None;
                }
                state.flush(Instant::now());
                if kind == Kind::Reduce {
                    if state.reduce_depth.get() == 0 {
                        state.calls.set(state.calls.get().saturating_add(1));
                    }
                    state
                        .reduce_depth
                        .set(state.reduce_depth.get().saturating_add(1));
                }
                state.mode.set(kind);
                Some(previous)
            })
            .ok()
            .flatten();
        Self {
            previous,
            kind,
            _not_send: PhantomData,
        }
    }
}

impl Drop for Guard {
    fn drop(&mut self) {
        let Some(previous) = self.previous else {
            return;
        };
        let kind = self.kind;
        let _ = STATE.try_with(|state| {
            state.flush(Instant::now());
            if kind == Kind::Reduce {
                state
                    .reduce_depth
                    .set(state.reduce_depth.get().saturating_sub(1));
            }
            state.mode.set(previous);
        });
    }
}

/// A sub-step scope ([`contract`] or [`finalize`]). Not `Send`; dropped in
/// reverse creation order.
#[must_use = "a timing guard measures until it is dropped"]
pub struct StepGuard {
    previous: Option<Step>,
    _not_send: PhantomData<*const ()>,
}

impl StepGuard {
    fn enter(step: Step) -> Self {
        let previous = STATE
            .try_with(|state| {
                let previous = state.step.get();
                if previous == step {
                    return None;
                }
                state.flush(Instant::now());
                state.step.set(step);
                Some(previous)
            })
            .ok()
            .flatten();
        Self {
            previous,
            _not_send: PhantomData,
        }
    }
}

impl Drop for StepGuard {
    fn drop(&mut self) {
        let Some(previous) = self.previous else {
            return;
        };
        let _ = STATE.try_with(|state| {
            state.flush(Instant::now());
            state.step.set(previous);
        });
    }
}

/// Enter a reduction scope (reduction proper and its finalisation).
pub fn reduce() -> Guard {
    Guard::enter(Kind::Reduce)
}

/// Enter a lowering scope (vakint's numerator-to-key work).
pub fn lower() -> Guard {
    Guard::enter(Kind::Lower)
}

/// Enter a contraction sub-step (numerator processing: idenso / spenso
/// called by vakint, a record's Taylor step and q averaging).
pub fn contract() -> StepGuard {
    StepGuard::enter(Step::Contract)
}

/// Enter a finalisation sub-step (master substitution, ε expansion,
/// assembly after a reduction).
pub fn finalize() -> StepGuard {
    StepGuard::enter(Step::Finalize)
}

/// The cumulative counters of this thread, including the running scope.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Snapshot {
    reduce_ns: u64,
    lower_ns: u64,
    calls: u64,
    step_ns: Matrix,
    cpu: Option<(u64, [u64; 4])>,
    stage_rss: Option<[u64; 4]>,
}

impl Snapshot {
    pub fn now() -> Self {
        STATE
            .try_with(|state| {
                let now = Instant::now();
                let resources = resources_enabled();
                if resources {
                    // fold the running interval into the resource counters
                    state.flush(now);
                }
                let mut snapshot = Snapshot {
                    reduce_ns: state.reduce_ns.get(),
                    lower_ns: state.lower_ns.get(),
                    calls: state.calls.get(),
                    step_ns: state.step_ns.get(),
                    cpu: resources.then(|| (state.cpu_ns.get(), state.stage_cpu_ns.get())),
                    stage_rss: resources.then(|| state.stage_rss.get()),
                };
                if !resources {
                    // the running interval, without touching the state
                    let running = nanos(state.since.get(), now);
                    let mode = state.mode.get();
                    match mode {
                        Kind::Reduce => snapshot.reduce_ns += running,
                        Kind::Lower => snapshot.lower_ns += running,
                        Kind::None => {}
                    }
                    if let Some(row) = state.step.get().row() {
                        snapshot.step_ns[row][mode.slot()] += running;
                    }
                }
                snapshot
            })
            .unwrap_or_default()
    }

    /// The counters accumulated on this thread since `self`.
    pub fn elapsed(&self) -> Delta {
        Snapshot::now().since(self)
    }

    /// `self - earlier`.
    pub fn since(&self, earlier: &Snapshot) -> Delta {
        let mut step_ns = [[0u64; 3]; 2];
        for (row, out) in step_ns.iter_mut().enumerate() {
            for (col, value) in out.iter_mut().enumerate() {
                *value = self.step_ns[row][col].saturating_sub(earlier.step_ns[row][col]);
            }
        }
        let cpu = match (self.cpu, earlier.cpu) {
            (Some((total, stages)), Some((total0, stages0))) => Some((
                total.saturating_sub(total0),
                std::array::from_fn(|i| stages[i].saturating_sub(stages0[i])),
            )),
            _ => None,
        };
        Delta {
            reduce_seconds: self.reduce_ns.saturating_sub(earlier.reduce_ns) as f64 * 1e-9,
            lower_seconds: self.lower_ns.saturating_sub(earlier.lower_ns) as f64 * 1e-9,
            reduce_calls: self.calls.saturating_sub(earlier.calls),
            reduce_ns: self.reduce_ns.saturating_sub(earlier.reduce_ns),
            lower_ns: self.lower_ns.saturating_sub(earlier.lower_ns),
            step_ns,
            cpu,
            stage_rss: self.stage_rss,
        }
    }

    pub fn reduce_seconds(&self) -> f64 {
        self.reduce_ns as f64 * 1e-9
    }

    pub fn lower_seconds(&self) -> f64 {
        self.lower_ns as f64 * 1e-9
    }

    pub fn reduce_calls(&self) -> u64 {
        self.calls
    }

    /// Cumulative contraction (numerator processing) seconds of this thread.
    pub fn contract_seconds(&self) -> f64 {
        self.step_ns[0].iter().sum::<u64>() as f64 * 1e-9
    }

    /// Cumulative finalisation seconds of this thread (part of
    /// [`Self::reduce_seconds`] where they run inside a reduction scope).
    pub fn finalize_seconds(&self) -> f64 {
        self.step_ns[1].iter().sum::<u64>() as f64 * 1e-9
    }
}

/// One stage of an interval ([`Delta::stage`]).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StageTime {
    /// Wall nanoseconds, exclusive.
    pub wall_ns: u64,
    /// Thread CPU nanoseconds, exclusive; `None` when the resources were not
    /// enabled for the whole interval.
    pub cpu_ns: Option<u64>,
    /// The largest process RSS sampled at the end of one of this stage's
    /// intervals on this thread since the resources were enabled (module
    /// documentation); `None` when not sampled.
    pub peak_rss_bytes: Option<u64>,
}

/// The counters of one interval.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Delta {
    pub reduce_seconds: f64,
    pub lower_seconds: f64,
    pub reduce_calls: u64,
    reduce_ns: u64,
    lower_ns: u64,
    step_ns: Matrix,
    cpu: Option<(u64, [u64; 4])>,
    stage_rss: Option<[u64; 4]>,
}

impl Delta {
    /// ` reduce_s=<s> lower_s=<s> contract_s=<s> finalize_s=<s>` (leading
    /// space), and ` cpu_s=<s>` (the thread's CPU over the interval) when the
    /// resources are enabled, for the log lines. The first two fields keep
    /// their text and order.
    pub fn fields(&self) -> String {
        let mut out = format!(
            " reduce_s={:.3} lower_s={:.3} contract_s={:.3} finalize_s={:.3}",
            self.reduce_seconds,
            self.lower_seconds,
            self.contract_seconds(),
            self.finalize_seconds()
        );
        if let Some((cpu, _)) = self.cpu {
            out.push_str(&format!(" cpu_s={:.3}", cpu as f64 * 1e-9));
        }
        out
    }

    /// Contraction (numerator processing) wall seconds, inclusive of where it
    /// ran (= [`Stage::Contraction`]).
    pub fn contract_seconds(&self) -> f64 {
        self.step_ns[0].iter().sum::<u64>() as f64 * 1e-9
    }

    /// Finalisation wall seconds (= [`Stage::Finalization`]; inside
    /// `reduce_seconds` where it runs in a reduction scope).
    pub fn finalize_seconds(&self) -> f64 {
        self.step_ns[1].iter().sum::<u64>() as f64 * 1e-9
    }

    /// The calling thread's CPU seconds over the interval (`None` unless the
    /// resources were enabled for the whole interval).
    pub fn cpu_seconds(&self) -> Option<f64> {
        self.cpu.map(|(cpu, _)| cpu as f64 * 1e-9)
    }

    /// One of the four EXCLUSIVE stages (module documentation): wall and
    /// thread-CPU nanoseconds, the sampled peak RSS.
    pub fn stage(&self, stage: Stage) -> StageTime {
        let m = &self.step_ns;
        let wall_ns = match stage {
            Stage::Contraction => m[0].iter().sum(),
            Stage::Finalization => m[1].iter().sum(),
            Stage::Lowering => self
                .lower_ns
                .saturating_sub(m[0][Kind::Lower.slot()] + m[1][Kind::Lower.slot()]),
            Stage::Reduction => self
                .reduce_ns
                .saturating_sub(m[0][Kind::Reduce.slot()] + m[1][Kind::Reduce.slot()]),
        };
        StageTime {
            wall_ns,
            cpu_ns: self.cpu.map(|(_, stages)| stages[stage.index()]),
            peak_rss_bytes: self
                .stage_rss
                .map(|rss| rss[stage.index()])
                .filter(|bytes| *bytes > 0),
        }
    }

    /// ` contraction_s=.. lowering_s=.. reduction_s=.. finalization_s=..`
    /// (the four exclusive wall stages; leading space).
    pub fn stage_fields(&self) -> String {
        let s = |stage| self.stage(stage).wall_ns as f64 * 1e-9;
        format!(
            " contraction_s={:.3} lowering_s={:.3} reduction_s={:.3} finalization_s={:.3}",
            s(Stage::Contraction),
            s(Stage::Lowering),
            s(Stage::Reduction),
            s(Stage::Finalization)
        )
    }
}

impl std::ops::AddAssign for Delta {
    fn add_assign(&mut self, other: Delta) {
        self.reduce_seconds += other.reduce_seconds;
        self.lower_seconds += other.lower_seconds;
        self.reduce_calls += other.reduce_calls;
        self.reduce_ns += other.reduce_ns;
        self.lower_ns += other.lower_ns;
        for row in 0..2 {
            for col in 0..3 {
                self.step_ns[row][col] += other.step_ns[row][col];
            }
        }
        // CPU: summed when both carry it, else whichever does (an interval
        // partly without resources reports what it has)
        self.cpu = match (self.cpu, other.cpu) {
            (Some((a, sa)), Some((b, sb))) => Some((a + b, std::array::from_fn(|i| sa[i] + sb[i]))),
            (a, b) => a.or(b),
        };
        self.stage_rss = match (self.stage_rss, other.stage_rss) {
            (Some(a), Some(b)) => Some(std::array::from_fn(|i| a[i].max(b[i]))),
            (a, b) => a.or(b),
        };
    }
}

/// Cumulative reduction seconds of this thread.
pub fn reduce_seconds() -> f64 {
    Snapshot::now().reduce_seconds()
}

/// Cumulative lowering seconds of this thread.
pub fn lower_seconds() -> f64 {
    Snapshot::now().lower_seconds()
}

/// Reduction scopes entered on this thread (outermost entries only).
pub fn reduce_calls() -> u64 {
    Snapshot::now().reduce_calls()
}

/// Cumulative contraction seconds of this thread.
pub fn contract_seconds() -> f64 {
    Snapshot::now().contract_seconds()
}

/// Cumulative finalisation seconds of this thread.
pub fn finalize_seconds() -> f64 {
    Snapshot::now().finalize_seconds()
}

/// Whether the factorised entries print one `VAKINT_RECORD` line per record
/// on stderr (`VAKINT_RECORD_LOG=1`, read once; off by default: no output).
pub fn record_log_enabled() -> bool {
    static ON: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| {
        matches!(
            std::env::var("VAKINT_RECORD_LOG").as_deref().map(str::trim),
            Ok("1")
        )
    });
    *ON
}

/// The per-record line the beta's runs read (their WH-RECORD-ANSWERS Q6):
/// `VAKINT_RECORD kind= label="" loops= depth= patterns= contract_s=
/// lower_s= reduce_s= finalize_s= cpu_s= wall_s= rss_gib=
/// contracted_bytes=`. `lower_s` / `reduce_s` are P12's inclusive counters,
/// `contract_s` / `finalize_s` the sub-steps; `cpu_s` and `rss_gib` (the
/// largest process RSS sampled at a stage boundary) are `-` unless the
/// resources are enabled; `patterns` is `-` where the entry has none.
#[allow(clippy::too_many_arguments)]
pub fn record_line(
    kind: &str,
    label: &str,
    loops: usize,
    depth: Option<i64>,
    patterns: Option<usize>,
    delta: &Delta,
    wall_seconds: f64,
    contracted_bytes: usize,
) -> String {
    let rss = Stage::ALL
        .iter()
        .filter_map(|stage| delta.stage(*stage).peak_rss_bytes)
        .max();
    let dash = |v: Option<String>| v.unwrap_or_else(|| "-".into());
    format!(
        "VAKINT_RECORD kind={kind} label={label:?} loops={loops} depth={} patterns={} contract_s={:.3} lower_s={:.3} reduce_s={:.3} finalize_s={:.3} cpu_s={} wall_s={wall_seconds:.3} rss_gib={} contracted_bytes={contracted_bytes}",
        dash(depth.map(|d| d.to_string())),
        dash(patterns.map(|p| p.to_string())),
        delta.contract_seconds(),
        delta.lower_seconds,
        delta.reduce_seconds,
        delta.finalize_seconds(),
        dash(delta.cpu_seconds().map(|c| format!("{c:.3}"))),
        dash(rss.map(|b| format!("{:.2}", b as f64 / (1u64 << 30) as f64))),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spin(ms: u64) {
        let until = Instant::now() + std::time::Duration::from_millis(ms);
        while Instant::now() < until {}
    }

    #[test]
    fn innermost_scope_wins_and_reentry_is_not_double_counted() {
        let start = Snapshot::now();
        {
            let _lower = lower();
            spin(20);
            {
                let _reduce = reduce();
                spin(30);
                {
                    let _again = reduce();
                    spin(10);
                }
                {
                    let _inner_lower = lower();
                    spin(10);
                }
            }
            spin(5);
        }
        spin(10);
        let delta = start.elapsed();
        assert!((delta.reduce_seconds - 0.040).abs() < 0.008, "{delta:?}");
        assert!((delta.lower_seconds - 0.035).abs() < 0.008, "{delta:?}");
        assert_eq!(delta.reduce_calls, 1);
    }

    /// The sub-step guards leave reduce / lower as P12 attributes them, and
    /// give the four exclusive stages.
    #[test]
    fn sub_steps_split_without_changing_reduce_and_lower() {
        let start = Snapshot::now();
        let wall = Instant::now();
        {
            let _outside = contract();
            spin(10); // contraction outside any kind
        }
        {
            let _lower = lower();
            spin(10); // lowering
            {
                let _contract = contract();
                spin(15); // contraction inside lowering
                {
                    let _again = contract();
                    spin(5); // same-kind nesting: still contraction
                }
            }
            {
                let _reduce = reduce();
                spin(20); // reduction
                {
                    let _finalize = finalize();
                    spin(25); // finalisation inside reduction
                    {
                        let _inner = contract();
                        spin(5); // innermost sub-step wins
                    }
                }
            }
        }
        let wall = wall.elapsed().as_secs_f64();
        let delta = start.elapsed();
        assert!((delta.lower_seconds - 0.030).abs() < 0.008, "{delta:?}");
        assert!((delta.reduce_seconds - 0.050).abs() < 0.008, "{delta:?}");
        let s = |stage| delta.stage(stage).wall_ns as f64 * 1e-9;
        assert!((s(Stage::Contraction) - 0.035).abs() < 0.008, "{delta:?}");
        assert!((s(Stage::Lowering) - 0.010).abs() < 0.006, "{delta:?}");
        assert!((s(Stage::Reduction) - 0.020).abs() < 0.006, "{delta:?}");
        assert!((s(Stage::Finalization) - 0.025).abs() < 0.006, "{delta:?}");
        let sum: f64 = Stage::ALL.iter().map(|stage| s(*stage)).sum();
        assert!(sum <= wall + 1e-6, "{sum} > wall {wall}");
        // exclusive stages add back up to P12's counters plus the
        // contraction outside every kind
        assert!(
            (sum - (delta.lower_seconds + delta.reduce_seconds + 0.010)).abs() < 0.004,
            "{delta:?}"
        );
        assert!((delta.contract_seconds() - s(Stage::Contraction)).abs() < 1e-9);
        assert!((delta.finalize_seconds() - s(Stage::Finalization)).abs() < 1e-9);
        // the log fields: P12's two first, in their order
        let fields = delta.fields();
        assert!(fields.starts_with(" reduce_s="), "{fields}");
        assert!(fields.contains(" lower_s=") && fields.contains(" contract_s="));
        assert!(fields.contains(" finalize_s="), "{fields}");
    }

    #[test]
    fn the_record_line_carries_every_field() {
        let start = Snapshot::now();
        {
            let _lower = lower();
            let _contract = contract();
            spin(1);
        }
        let line = record_line("direct", "x", 2, Some(3), None, &start.elapsed(), 0.5, 1234);
        for field in [
            "VAKINT_RECORD kind=direct label=\"x\" loops=2 depth=3 patterns=- contract_s=",
            " lower_s=",
            " reduce_s=",
            " finalize_s=",
            " cpu_s=",
            " wall_s=0.500",
            " rss_gib=",
            " contracted_bytes=1234",
        ] {
            assert!(line.contains(field), "{line}");
        }
    }

    #[test]
    fn counters_are_monotone() {
        let mut previous = Snapshot::now();
        for _ in 0..5 {
            {
                let _lower = lower();
                let _contract = contract();
                spin(1);
            }
            let now = Snapshot::now();
            let delta = now.since(&previous);
            assert!(delta.stage(Stage::Contraction).wall_ns > 0);
            assert!(now.lower_seconds() >= previous.lower_seconds());
            assert!(now.contract_seconds() >= previous.contract_seconds());
            previous = now;
        }
    }

    /// Resources off: no CPU, no RSS, no cpu_s field. On (own thread): the
    /// CPU of a busy single-threaded stage is close to its wall time, and the
    /// RSS is sampled.
    #[test]
    fn resources_follow_the_wall_time_on_one_thread() {
        std::thread::spawn(|| {
            set_resources(false);
            let start = Snapshot::now();
            {
                let _reduce = reduce();
                spin(5);
            }
            let off = start.elapsed();
            assert_eq!(off.cpu_seconds(), None);
            assert_eq!(off.stage(Stage::Reduction).cpu_ns, None);
            assert!(!off.fields().contains("cpu_s="));
            set_resources(true);
            let start = Snapshot::now();
            {
                let _reduce = reduce();
                let _finalize = finalize();
                spin(40);
            }
            let delta = start.elapsed();
            let stage = delta.stage(Stage::Finalization);
            let cpu = stage.cpu_ns.unwrap() as f64 * 1e-9;
            assert!(cpu > 0.025, "{stage:?}");
            assert!(cpu <= stage.wall_ns as f64 * 1e-9 + 0.005, "{stage:?}");
            assert!(stage.peak_rss_bytes.unwrap() > 0, "{stage:?}");
            assert!(delta.cpu_seconds().unwrap() >= cpu);
            assert!(delta.fields().contains(" cpu_s="));
            assert!(process_peak_rss_bytes().unwrap() > 0);
            set_resources(false);
        })
        .join()
        .unwrap();
    }
}
