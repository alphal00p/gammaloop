//! The per-diagram accounting of `GAMMALOOP_UNCONTRACTED_NUMERATOR` (U7a
//! delta: Astra's U7a verdict finding 4, helper VH's round 2 R2-I1, R2-M2
//! and R2-M3; the beta project's results/consolidation/U7-PLAN.md
//! requirements 16 and 17).
//!
//! # The ledger of a diagram
//!
//! A diagram's records are counted in ONE ledger ([`DiagramLedger`]), keyed
//! by the graph name and the diagram key: the beta opens one per diagram job
//! with the diagram's integral index (the `NNNNNN` of `integral-NNNNNN.json`);
//! a root evaluation outside such a ledger opens its own
//! (`diagram=root#<n>`, a process-wide ordinal), so that two root
//! evaluations of one graph name never share or split a line (R2-M2). The
//! ledger is a thread-local stack of frames: an evaluation and the R* of its
//! UV children run on one thread (a plan computes its children in a serial
//! loop, a Δ-cache fill runs on the claiming thread), so every record of the
//! diagram lands in its frame, whichever of its root evaluations (insertion
//! surrogates included) made it. The frames are restored on unwinding (a
//! guard), so a panicking diagram never leaves a frame to the next diagram
//! of its worker. A record outside every frame of its thread prints an
//! `UNCONTRACTED_ORPHAN` line (never silently uncounted).
//!
//! Every EXPECTED record (a remainder of an enabled kind at the hook) is
//! exactly one disposition:
//!
//! * `compared` (check: our value compared with today's), `taken` (on:
//!   vakint's value, or a zero Δ, replaced today's);
//! * `declined:<code>` (a designed refusal; today's path took it);
//! * `error:<code>` (not designed: the diagram fails);
//! * `unchecked` (reached the hook but never compared: check mode whose
//!   today's path failed after our attempt, or a count the ledger misses).
//!
//! A K3 term that builds no record by design (`zero_loops`: Δ(γ) inserted,
//! no vakint call today; `order_zero`: K8) is not expected; the lines list
//! it as `no_record=<code>:<n>`.
//!
//! The close prints ONE `UNCONTRACTED_DISPOSITION` line per diagram and
//! enabled kind (also when nothing was expected, so a missing line is a
//! missing diagram) and one `UNCONTRACTED_SUMMARY` line per diagram. On and
//! check are strict (the only behaviour, R2-I1): an error outcome fails the
//! diagram.

use std::{
    cell::RefCell,
    collections::BTreeMap,
    sync::atomic::{AtomicUsize, Ordering},
};

use eyre::{Result, eyre};
use linnet::half_edge::subgraph::SuBitGraph;
use symbolica::atom::Atom;

use crate::{graph::Graph, model::Model, uv::counterterm_cache::CountertermCache};

use super::uncontracted::{self, Attempt, Check, Kind, Mode, Outcome, Switch};

// ---------------------------------------------------------------------------
// The switches of the ledger: strict (a no-op)
// ---------------------------------------------------------------------------

/// `GAMMALOOP_UNCONTRACTED_STRICT` (R2-I1): on and check always fail a
/// diagram on an error outcome, so `1` / `on` (and unset, empty) are
/// accepted as a no-op and `0` / `off` are refused with the reason.
pub(crate) fn parse_strict(value: Option<&str>) -> std::result::Result<(), String> {
    match value.map(str::trim) {
        None | Some("" | "1" | "on") => Ok(()),
        Some("0" | "off") => Err(format!(
            "{}={} is refused: the non-strict mode is removed (helper VH round 2 R2-I1); with {} on or check an error outcome always fails its diagram. Unset {} (1 is accepted as a no-op)",
            uncontracted::STRICT_VARIABLE,
            value.unwrap_or_default().trim(),
            uncontracted::VARIABLE,
            uncontracted::STRICT_VARIABLE
        )),
        Some(other) => Err(format!(
            "{} must be unset or 1 (a no-op; strict is the only behaviour), not `{other}`",
            uncontracted::STRICT_VARIABLE
        )),
    }
}

/// What the switch's initialisation checks besides the switch itself
/// ([`uncontracted::validate`]): `GAMMALOOP_UNCONTRACTED_STRICT` (the
/// environment's value), refused when `0` in every mode.
pub(crate) fn initialise(_switch: Switch, strict: Option<&str>) -> std::result::Result<(), String> {
    parse_strict(strict)
}

// ---------------------------------------------------------------------------
// Counts and frames
// ---------------------------------------------------------------------------

/// The records of one kind in a frame: the dispositions (exclusive; module
/// doc) and what the attempts gave.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct Counts {
    /// Records that reached the hook.
    pub expected: usize,
    /// Records built and handed to vakint, or refused before (this frame).
    pub attempted: usize,
    /// Attempts vakint evaluated (outcome uncontracted) and zero Δs, in
    /// either mode (check mode compares them).
    pub uncontracted: usize,
    pub zero: usize,
    /// Comparisons that agreed.
    pub equal: usize,
    /// Disposition -> records (`compared`, `taken`, `declined:<code>`,
    /// `error:<code>`, `unchecked`).
    pub classes: BTreeMap<String, usize>,
    /// K3 terms of the kind that build no record by design, by code
    /// (`zero_loops`: Δ(γ) inserted, no vakint call today; `order_zero`: K8):
    /// not expected and never compared, today's path takes them.
    pub no_record: BTreeMap<String, usize>,
}

impl Counts {
    fn add(&mut self, class: impl Into<String>, n: usize) {
        if n > 0 {
            *self.classes.entry(class.into()).or_default() += n;
        }
    }

    fn merge(&mut self, other: &Counts) {
        self.expected += other.expected;
        self.attempted += other.attempted;
        self.uncontracted += other.uncontracted;
        self.zero += other.zero;
        self.equal += other.equal;
        for (class, n) in &other.classes {
            self.add(class.clone(), *n);
        }
        for (code, n) in &other.no_record {
            *self.no_record.entry(code.clone()).or_default() += n;
        }
    }

    fn class(&self, name: &str) -> usize {
        self.classes.get(name).copied().unwrap_or(0)
    }

    /// The classes `<prefix><rest>` as `rest -> n`.
    fn with_prefix(&self, prefix: &str) -> BTreeMap<String, usize> {
        self.classes
            .iter()
            .filter_map(|(class, n)| Some((class.strip_prefix(prefix)?.to_owned(), *n)))
            .collect()
    }

    fn accounted(&self) -> usize {
        self.classes.values().sum()
    }

    /// `unchecked` records and any record no disposition accounts for.
    fn unchecked(&self) -> usize {
        self.class("unchecked") + self.expected.saturating_sub(self.accounted())
    }

    fn compared(&self) -> usize {
        self.class("compared")
    }
}

/// A frame of this thread's ledger stack.
enum Scope {
    /// A diagram (closed by its [`DiagramLedger`]).
    Diagram { graph: String, diagram: String },
}

struct Frame {
    id: usize,
    scope: Scope,
    counts: BTreeMap<Kind, Counts>,
}

thread_local! {
    static FRAMES: RefCell<Vec<Frame>> = const { RefCell::new(Vec::new()) };
    /// The evaluation path (module doc), one entry per R* child fill.
    static PATH: RefCell<Vec<String>> = const { RefCell::new(Vec::new()) };
}

static FRAME_IDS: AtomicUsize = AtomicUsize::new(1);

/// Whether the ledger runs: the switch on or check.
fn active() -> bool {
    uncontracted::switch().mode != Mode::Off
}

fn push(scope: Scope) -> usize {
    let id = FRAME_IDS.fetch_add(1, Ordering::Relaxed);
    FRAMES.with(|frames| {
        frames.borrow_mut().push(Frame {
            id,
            scope,
            counts: BTreeMap::new(),
        })
    });
    id
}

/// Remove the frame `id` (and any frame left above it); never panics (it
/// runs on unwinding).
fn pop(id: usize) -> Option<Frame> {
    FRAMES
        .try_with(|frames| {
            let mut frames = frames.try_borrow_mut().ok()?;
            let at = frames.iter().rposition(|frame| frame.id == id)?;
            let mut popped = frames.split_off(at);
            Some(popped.swap_remove(0))
        })
        .ok()
        .flatten()
}

/// `f` on the innermost frame (nothing without one: a record outside any
/// diagram or fill, e.g. a test calling the emitter directly).
fn with_top(f: impl FnOnce(&mut Frame)) {
    let _ = FRAMES.try_with(|frames| {
        if let Ok(mut frames) = frames.try_borrow_mut()
            && let Some(frame) = frames.last_mut()
        {
            f(frame);
        }
    });
}

fn with_counts(kind: Kind, f: impl FnOnce(&mut Counts)) {
    let mut counted = false;
    with_top(|frame| {
        counted = true;
        f(frame.counts.entry(kind).or_default())
    });
    if !counted && uncontracted::switch().mode != Mode::Off {
        // a record outside every diagram ledger and fill of this thread (a
        // test calling the emitter directly; in a run, a thread the ledger
        // does not reach): loud, never silently uncounted
        eprintln!(
            "UNCONTRACTED_ORPHAN kind={} thread={}",
            kind.as_str(),
            std::thread::current().name().unwrap_or("?")
        );
    }
}

fn in_diagram() -> bool {
    FRAMES.with(|frames| {
        frames
            .borrow()
            .iter()
            .any(|frame| matches!(frame.scope, Scope::Diagram { .. }))
    })
}

// ---------------------------------------------------------------------------
// The records
// ---------------------------------------------------------------------------

/// A remainder of the enabled `kind` reached the hook.
pub(super) fn expect(kind: Kind) {
    with_counts(kind, |counts| counts.expected += 1);
}

/// The disposition of one reported attempt ([`uncontracted::report`]).
pub(super) fn record(attempt: &Attempt, check: Option<&Check>, mode: Mode) {
    with_counts(attempt.kind, |counts| {
        counts.attempted += 1;
        let class = match &attempt.outcome {
            Outcome::Value(_) => {
                counts.uncontracted += 1;
                None
            }
            Outcome::Zero => {
                counts.zero += 1;
                None
            }
            Outcome::Declined(decline) => Some(format!("declined:{}", decline.code)),
            Outcome::Error(decline) => Some(format!("error:{}", decline.code)),
        };
        let class = class.unwrap_or_else(|| match (mode, check) {
            (Mode::Check, Some(_)) => "compared".to_owned(),
            // check mode without a comparison: today's path failed after
            // our attempt
            (Mode::Check, None) => "unchecked".to_owned(),
            _ => "taken".to_owned(),
        });
        counts.equal += usize::from(check.is_some_and(|check| check.equal));
        counts.add(class, 1);
    });
}

/// vakint's unknown-master refusal of an attempt of `kind` (a hard error).
pub(super) fn record_unknown_masters(kind: Kind) {
    with_counts(kind, |counts| {
        counts.attempted += 1;
        counts.add("error:unknown_masters", 1);
    });
}

/// A K3 term of the enabled `kind` that builds no record by design (not
/// expected: listed as `no_record=`).
pub(super) fn record_no_record(kind: Kind, code: &str) {
    with_counts(kind, |counts| {
        *counts.no_record.entry(code.to_owned()).or_default() += 1
    });
}

/// `cache.counterterm_with` for an R* child Δ (`RStarPlan::rstar_child_counterterm`):
/// with the switch on or check, the evaluation path names the child inside
/// the fill (its records' labels). Off: exactly `cache.counterterm_with`.
pub(super) fn through_cache(
    cache: &CountertermCache,
    graph: &Graph,
    support: &SuBitGraph,
    model: &Model,
    settings_key: String,
    compute: impl FnOnce() -> Result<Atom>,
) -> Result<Atom> {
    if uncontracted::switch().mode == Mode::Off {
        return cache.counterterm_with(graph, support, model, settings_key, compute);
    }
    cache.counterterm_with(graph, support, model, settings_key, || {
        let _path = ChildPath::enter(graph, support);
        compute()
    })
}

/// The evaluation path inside an R* child fill.
struct ChildPath;

impl ChildPath {
    fn enter(graph: &Graph, support: &SuBitGraph) -> Self {
        let segment = format!(
            "{}/child{{{}}}",
            current_path(),
            super::plan::join_edge_ids(&super::plan::uv_child_edge_ids(graph, support))
        );
        PATH.with(|path| path.borrow_mut().push(segment));
        ChildPath
    }
}

impl Drop for ChildPath {
    fn drop(&mut self) {
        let _ = PATH.try_with(|path| {
            if let Ok(mut path) = path.try_borrow_mut() {
                path.pop();
            }
        });
    }
}

/// The evaluation path of this thread (`root` outside every child fill).
pub(super) fn current_path() -> String {
    PATH.try_with(|path| path.try_borrow().ok().and_then(|path| path.last().cloned()))
        .ok()
        .flatten()
        .unwrap_or_else(|| "root".to_owned())
}

// ---------------------------------------------------------------------------
// The diagram ledger
// ---------------------------------------------------------------------------

/// One diagram's ledger (module doc): opened by the beta per diagram job, or
/// by a root evaluation outside one. [`Self::close`] prints its lines;
/// dropped without a close (an error
/// return or an unwind) it prints them with `outcome=error`.
pub struct DiagramLedger {
    id: usize,
    closed: bool,
}

impl DiagramLedger {
    /// A ledger for the diagram `diagram` (the beta: its integral index) of
    /// `graph`; `None` with the switch off.
    pub fn open(graph: &str, diagram: &str) -> Option<Self> {
        active().then(|| Self {
            id: push(Scope::Diagram {
                graph: graph.to_owned(),
                diagram: diagram.to_owned(),
            }),
            closed: false,
        })
    }

    /// The ledger of a root evaluation of `graph` that runs outside every
    /// diagram ledger (`diagram=root#<n>`); `None` inside one.
    pub(super) fn implicit(graph: &str) -> Option<Self> {
        static ROOTS: AtomicUsize = AtomicUsize::new(1);
        if !active() || in_diagram() {
            return None;
        }
        Self::open(
            graph,
            &format!("root#{}", ROOTS.fetch_add(1, Ordering::Relaxed)),
        )
    }

    /// Close the ledger of a diagram whose evaluation succeeded (`ok`) or
    /// failed: one `UNCONTRACTED_DISPOSITION` line per enabled kind, one
    /// `UNCONTRACTED_SUMMARY` line.
    pub fn close(mut self, ok: bool) -> Result<()> {
        self.closed = true;
        let Some(frame) = pop(self.id) else {
            return Err(eyre!(
                "UNCONTRACTED: the diagram ledger is not on this thread's stack"
            ));
        };
        finish(frame, ok)
    }
}

impl Drop for DiagramLedger {
    fn drop(&mut self) {
        if !self.closed
            && let Some(frame) = pop(self.id)
        {
            let _ = finish(frame, false);
        }
    }
}

/// The verdict on one diagram-kind.
struct Verdict {
    outcome: &'static str,
}

fn list(map: &BTreeMap<String, usize>) -> String {
    if map.is_empty() {
        "none".to_owned()
    } else {
        map.iter()
            .map(|(code, n)| format!("{code}:{n}"))
            .collect::<Vec<_>>()
            .join(",")
    }
}

fn verdict(counts: &Counts, ok: bool) -> Verdict {
    let outcome = if !ok || !counts.with_prefix("error:").is_empty() {
        // the diagram fails with its own error
        "error"
    } else if counts.unchecked() > 0 {
        "unchecked"
    } else {
        "ok"
    };
    Verdict { outcome }
}

/// The `UNCONTRACTED_DISPOSITION` line of one diagram-kind.
fn disposition_line(
    graph: &str,
    diagram: &str,
    mode: Mode,
    kind: Kind,
    counts: &Counts,
    verdict: &Verdict,
) -> String {
    format!(
        "UNCONTRACTED_DISPOSITION graph={graph} diagram={diagram} mode={} kind={} expected={} attempted={} uncontracted={} zero={} compared={} equal={} declined={} no_record={} unchecked={} error={} outcome={}",
        mode.as_str(),
        kind.as_str(),
        counts.expected,
        counts.attempted,
        counts.uncontracted,
        counts.zero,
        counts.compared(),
        counts.equal,
        list(&counts.with_prefix("declined:")),
        list(&counts.no_record),
        counts.unchecked(),
        list(&counts.with_prefix("error:")),
        verdict.outcome
    )
}

/// The `UNCONTRACTED_SUMMARY` line of a diagram (its kinds summed).
fn summary_line(
    graph: &str,
    diagram: &str,
    mode: Mode,
    kinds: &[Kind],
    total: &Counts,
    outcome: &str,
) -> String {
    format!(
        "UNCONTRACTED_SUMMARY graph={graph} diagram={diagram} mode={} kinds={} attempted={} uncontracted={} zero={} compared={} equal={} declined={} error={} expected={} no_record={} unchecked={} outcome={outcome}",
        mode.as_str(),
        kinds
            .iter()
            .map(|kind| kind.as_str())
            .collect::<Vec<_>>()
            .join(","),
        total.attempted,
        total.uncontracted,
        total.zero,
        total.compared(),
        total.equal,
        list(&total.with_prefix("declined:")),
        list(&total.with_prefix("error:")),
        total.expected,
        list(&total.no_record),
        total.unchecked(),
    )
}

/// The diagram outcome of its kinds' outcomes (the worst).
fn worst<'o>(outcomes: impl IntoIterator<Item = &'o str>) -> &'static str {
    let rank = |outcome: &str| match outcome {
        "error" => 2,
        "unchecked" => 1,
        _ => 0,
    };
    match outcomes.into_iter().map(rank).max().unwrap_or(0) {
        2 => "error",
        1 => "unchecked",
        _ => "ok",
    }
}

fn print(line: String) {
    eprintln!("{line}");
    #[cfg(test)]
    {
        if line.starts_with("UNCONTRACTED_SUMMARY ") {
            SUMMARIES.with(|lines| lines.borrow_mut().push(line.clone()));
        }
        LINES.with(|lines| lines.borrow_mut().push(line));
    }
}

#[cfg(test)]
thread_local! {
    static LINES: RefCell<Vec<String>> = const { RefCell::new(Vec::new()) };
    static SUMMARIES: RefCell<Vec<String>> = const { RefCell::new(Vec::new()) };
}

/// Tests: the `UNCONTRACTED_SUMMARY` lines printed on this thread since the
/// last call (kept apart from [`take_lines`]).
#[cfg(test)]
pub(crate) fn take_summaries() -> Vec<String> {
    SUMMARIES.with(|lines| std::mem::take(&mut *lines.borrow_mut()))
}

/// Tests: the ledger's lines (dispositions, summaries, resumed lines)
/// printed on this thread since the last call.
#[cfg(test)]
pub(crate) fn take_lines() -> Vec<String> {
    LINES.with(|lines| std::mem::take(&mut *lines.borrow_mut()))
}

/// The lines of a closed diagram frame.
fn finish(frame: Frame, ok: bool) -> Result<()> {
    let switch = uncontracted::switch();
    if switch.mode == Mode::Off {
        return Ok(());
    }
    let Scope::Diagram { graph, diagram } = frame.scope;
    let kinds = switch.enabled_kinds();
    let mut total = Counts::default();
    let mut outcomes = Vec::new();
    for kind in &kinds {
        let counts = frame.counts.get(kind).cloned().unwrap_or_default();
        let verdict = verdict(&counts, ok);
        print(disposition_line(
            &graph,
            &diagram,
            switch.mode,
            *kind,
            &counts,
            &verdict,
        ));
        outcomes.push(verdict.outcome);
        total.merge(&counts);
    }
    print(summary_line(
        &graph,
        &diagram,
        switch.mode,
        &kinds,
        &total,
        worst(outcomes),
    ));
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests: an injected outcome at the hook
// ---------------------------------------------------------------------------

#[cfg(test)]
thread_local! {
    static FAULT: RefCell<Option<(Kind, String, bool)>> = const { RefCell::new(None) };
}

/// Tests: run `f` with every attempt of `kind` refused by `code` at the
/// hook (`designed`: a designed decline, else an error outcome), as if the
/// emitter or vakint's entry had refused it.
#[cfg(test)]
pub(crate) fn with_fault<T>(kind: Kind, code: &str, designed: bool, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<(Kind, String, bool)>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0.take();
            FAULT.with(|slot| *slot.borrow_mut() = previous);
        }
    }
    let _restore =
        Restore(FAULT.with(|slot| slot.replace(Some((kind, code.to_owned(), designed)))));
    f()
}

/// Tests: the attempt with the injected outcome of [`with_fault`], if any.
#[cfg(test)]
pub(super) fn inject(mut attempt: Attempt) -> Attempt {
    FAULT.with(|slot| {
        if let Some((kind, code, designed)) = slot.borrow().as_ref()
            && *kind == attempt.kind
        {
            let decline = uncontracted::Decline {
                code: code.clone(),
                detail: "injected by the test".to_owned(),
                designed: *designed,
            };
            attempt.text = None;
            attempt.outcome = if *designed {
                Outcome::Declined(decline)
            } else {
                Outcome::Error(decline)
            };
        }
    });
    attempt
}
