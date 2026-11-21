//! `GAMMALOOP_STAGE_TIME=1`: one machine-readable `STAGE_SPAN` line per timed
//! unit (a *span*: a diagram, an R* evaluation, a UV child, a remainder term,
//! a record, an IR hard child, a quotient, a boundary pass, a surrogate).
//! LOGGING ONLY: no value, cache or control flow depends on anything here.
//!
//! Each span reports only its SELF time, split into stage leaves of the
//! three fixed categories (`N` numerator processing, `L` lowering, `R`
//! reduction; timing-breakdown DESIGN §2), so that for every span
//!
//! ```text
//! wall_ns = Σ self leaves + Σ children wall_ns          (exactly, in ns)
//! ```
//!
//! and the leaves of a diagram sum to its total. The leaves:
//!
//! * vakint's own exclusive stages ([`vakint::timing::Delta::stage`]) of the
//!   span minus those of its children: `N.contract.vakint` (Contraction),
//!   `L.vakint` (Lowering), `R.ibp` (Reduction), `R.finalise`
//!   (Finalisation);
//! * gammaloop laps ([`lap`], [`Mark`]): wall time of a named scope NET of the
//!   vakint stages and the child spans inside it (`N.taylor`, `N.project`,
//!   `N.contract.gl`, `N.convert`, `N.record.build`, ...);
//! * `N.other`: the explicit residual `wall − Σ children − Σ leaves`. A
//!   negative residual (clock reads of nested scopes) is emitted as 0 with
//!   `overlap_ns=`.
//!
//! The counters are thread-local, like vakint's: a span measures the work of
//! the thread that opened it (all production and gate runs use one thread per
//! diagram, `--tensor-parallelism serial`).
//!
//! Occurrence keys: `<fnv64 hex>#<ordinal>` over (parent key, site, graph,
//! label), the ordinal counting equal digests within the diagram (reset by
//! [`diagram`]); a remainder-term span takes the uncontracted ledger's
//! occurrence when the ledger has one (`GAMMALOOP_UNCONTRACTED_TIMING=1` or
//! the switch on), so OFF and ON runs pair by key.

use std::cell::RefCell;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::marker::PhantomData;
use std::sync::LazyLock;
use std::time::Instant;

use vakint::timing::{Snapshot, Stage};

pub(crate) const VARIABLE: &str = "GAMMALOOP_STAGE_TIME";

static ENABLED: LazyLock<bool> = LazyLock::new(|| {
    matches!(
        std::env::var(VARIABLE).as_deref().map(str::trim),
        Ok("1" | "on" | "true")
    )
});

/// Whether `GAMMALOOP_STAGE_TIME=1` (read once).
pub fn enabled() -> bool {
    *ENABLED
}

/// Display order of the self leaves (taxonomy order, DESIGN §2); unknown lap
/// names follow the `N` ones.
const ORDER: &[&str] = &[
    "N.forest",
    "N.taylor",
    "N.project",
    "N.contract.gl",
    "N.contract.vakint",
    "N.convert",
    "N.boundary.build",
    "N.record.build",
    "N.record.finish",
    "N.backend",
    "N.finish",
    "N.other",
    "L.vakint",
    "R.ibp",
    "R.finalise",
];

fn rank(leaf: &str) -> usize {
    ORDER
        .iter()
        .position(|known| *known == leaf)
        .unwrap_or(ORDER.len() - 4)
}

fn vakint_ns(delta: &vakint::timing::Delta) -> [u64; 4] {
    [
        delta.stage(Stage::Contraction).wall_ns,
        delta.stage(Stage::Lowering).wall_ns,
        delta.stage(Stage::Reduction).wall_ns,
        delta.stage(Stage::Finalization).wall_ns,
    ]
}

fn vakint_total(snapshot: &Snapshot) -> (u64, u64) {
    let delta = snapshot.since(&Snapshot::default());
    let ns = vakint_ns(&delta);
    (ns.iter().sum(), delta.reduce_calls)
}

struct Frame {
    seq: u64,
    site: &'static str,
    graph: String,
    occurrence: String,
    parent: String,
    depth: usize,
    fields: Vec<(&'static str, String)>,
    started: Instant,
    snapshot: Snapshot,
    child_wall: u64,
    child_vakint: [u64; 4],
    child_calls: u64,
    child_spans: u64,
    laps: Vec<(&'static str, u64)>,
    lap_open: bool,
}

#[derive(Default)]
struct ThreadState {
    stack: Vec<Frame>,
    diagram: String,
    ordinals: HashMap<u64, u32>,
    seq: u64,
}

thread_local! {
    static STATE: RefCell<ThreadState> = RefCell::new(ThreadState::default());
}

/// An open span; it closes (and prints its line) when dropped. Not `Send`.
#[must_use = "a span measures until it is dropped"]
pub struct Span {
    seq: Option<u64>,
    _not_send: PhantomData<*const ()>,
}

impl Span {
    const INERT: Span = Span {
        seq: None,
        _not_send: PhantomData,
    };

    /// Open a span of `site` on `graph`; `label` (the unit's identity within
    /// its parent: edges, loops, ids) is only built when on.
    pub fn open(site: &'static str, graph: &str, label: impl FnOnce() -> String) -> Span {
        if !enabled() {
            return Span::INERT;
        }
        let label = label();
        let started = Instant::now();
        let snapshot = Snapshot::now();
        STATE
            .try_with(|state| {
                let Ok(mut state) = state.try_borrow_mut() else {
                    return Span::INERT;
                };
                let parent = state
                    .stack
                    .last()
                    .map_or_else(|| "-".to_owned(), |frame| frame.occurrence.clone());
                let digest = super::uncontracted_ledger::fnv64([
                    parent.as_bytes(),
                    site.as_bytes(),
                    graph.as_bytes(),
                    label.as_bytes(),
                ]);
                let ordinal = {
                    let next = state.ordinals.entry(digest).or_insert(0);
                    *next += 1;
                    *next
                };
                state.seq += 1;
                let seq = state.seq;
                let depth = state.stack.len();
                state.stack.push(Frame {
                    seq,
                    site,
                    graph: graph.to_owned(),
                    occurrence: format!("{digest:016x}#{ordinal}"),
                    parent,
                    depth,
                    fields: Vec::new(),
                    started,
                    snapshot,
                    child_wall: 0,
                    child_vakint: [0; 4],
                    child_calls: 0,
                    child_spans: 0,
                    laps: Vec::new(),
                    lap_open: false,
                });
                Span {
                    seq: Some(seq),
                    _not_send: PhantomData,
                }
            })
            .unwrap_or(Span::INERT)
    }

    fn with_frame(&self, f: impl FnOnce(&mut Frame)) {
        let Some(seq) = self.seq else {
            return;
        };
        let _ = STATE.try_with(|state| {
            if let Ok(mut state) = state.try_borrow_mut()
                && let Some(frame) = state.stack.iter_mut().rev().find(|frame| frame.seq == seq)
            {
                f(frame);
            }
        });
    }

    /// A `key=value` field of the line (printed in insertion order; a later
    /// value of the same key replaces the earlier one).
    pub fn set(&self, key: &'static str, value: impl std::fmt::Display) {
        if self.seq.is_none() {
            return;
        }
        let value = value.to_string();
        self.with_frame(|frame| {
            if let Some(slot) = frame.fields.iter_mut().find(|(k, _)| *k == key) {
                slot.1 = value;
            } else {
                frame.fields.push((key, value));
            }
        });
    }

    /// Use the uncontracted ledger's occurrence key (OFF / ON pairing).
    pub fn set_occurrence(&self, occurrence: &str) {
        if self.seq.is_none() || occurrence.is_empty() {
            return;
        }
        let occurrence = occurrence.to_owned();
        self.with_frame(|frame| frame.occurrence = occurrence);
    }
}

impl Drop for Span {
    fn drop(&mut self) {
        let Some(seq) = self.seq else {
            return;
        };
        let ended = Instant::now();
        let snapshot = Snapshot::now();
        let _ = STATE.try_with(|state| {
            let Ok(mut state) = state.try_borrow_mut() else {
                return;
            };
            if !state.stack.iter().any(|frame| frame.seq == seq) {
                return;
            }
            // Close every frame above this one first (a span dropped out of
            // order): each is flagged.
            while let Some(frame) = state.stack.pop() {
                let own = frame.seq == seq;
                let (line, wall, vakint, calls, spans) = close(frame, ended, &snapshot, !own, &state.diagram);
                eprintln!("{line}");
                if let Some(parent) = state.stack.last_mut() {
                    parent.child_wall = parent.child_wall.saturating_add(wall);
                    for (slot, ns) in parent.child_vakint.iter_mut().zip(vakint) {
                        *slot = slot.saturating_add(ns);
                    }
                    parent.child_calls = parent.child_calls.saturating_add(calls);
                    parent.child_spans = parent.child_spans.saturating_add(spans + 1);
                }
                if own {
                    break;
                }
            }
            if state.stack.is_empty() {
                // the diagram (or a span opened outside one) is over
                state.diagram.clear();
            }
        });
    }
}

/// The line of a closing frame; also its wall, inclusive vakint stage ns,
/// reduce calls and descendant span count for the parent.
fn close(
    frame: Frame,
    ended: Instant,
    snapshot: &Snapshot,
    out_of_order: bool,
    diagram: &str,
) -> (String, u64, [u64; 4], u64, u64) {
    let wall = ended.saturating_duration_since(frame.started).as_nanos() as u64;
    let delta = snapshot.since(&frame.snapshot);
    let vakint = vakint_ns(&delta);
    let mut leaves: Vec<(&'static str, u64)> = Vec::new();
    fn add(leaves: &mut Vec<(&'static str, u64)>, name: &'static str, ns: u64) {
        if ns == 0 {
            return;
        }
        if let Some(slot) = leaves.iter_mut().find(|(leaf, _)| *leaf == name) {
            slot.1 = slot.1.saturating_add(ns);
        } else {
            leaves.push((name, ns));
        }
    }
    let own_vakint: [u64; 4] =
        std::array::from_fn(|i| vakint[i].saturating_sub(frame.child_vakint[i]));
    for (name, ns) in &frame.laps {
        add(&mut leaves, name, *ns);
    }
    for (name, ns) in ["N.contract.vakint", "L.vakint", "R.ibp", "R.finalise"]
        .into_iter()
        .zip(own_vakint)
    {
        add(&mut leaves, name, ns);
    }
    let claimed: u64 = leaves.iter().map(|(_, ns)| *ns).sum();
    let accounted = claimed.saturating_add(frame.child_wall);
    let overlap = accounted.saturating_sub(wall);
    add(&mut leaves, "N.other", wall.saturating_sub(accounted));
    leaves.sort_by_key(|(leaf, _)| rank(leaf));
    let calls = delta.reduce_calls;
    let mut line = String::with_capacity(256);
    let _ = write!(
        line,
        "STAGE_SPAN v=1 graph={} diagram={} site={} occurrence={} parent={} depth={}",
        frame.graph,
        if diagram.is_empty() { "-" } else { diagram },
        frame.site,
        frame.occurrence,
        frame.parent,
        frame.depth
    );
    for (key, value) in &frame.fields {
        let _ = write!(line, " {key}={value}");
    }
    let _ = write!(line, " wall_ns={wall} self=");
    if leaves.is_empty() {
        line.push('-');
    }
    for (n, (leaf, ns)) in leaves.iter().enumerate() {
        if n > 0 {
            line.push(',');
        }
        let _ = write!(line, "{leaf}:{ns}");
    }
    let _ = write!(
        line,
        " calls=reduce:{calls},self_reduce:{},spans:{}",
        calls.saturating_sub(frame.child_calls),
        frame.child_spans
    );
    if overlap > 0 {
        let _ = write!(line, " overlap_ns={overlap}");
    }
    if out_of_order {
        line.push_str(" flag=out_of_order");
    }
    (line, wall, vakint, calls, frame.child_spans)
}

/// A field of the innermost open span of this thread (e.g. `cache=fill` from
/// inside a cache fill closure).
pub fn set_top(key: &'static str, value: &str) {
    if !enabled() {
        return;
    }
    let _ = STATE.try_with(|state| {
        if let Ok(mut state) = state.try_borrow_mut()
            && let Some(frame) = state.stack.last_mut()
        {
            if let Some(slot) = frame.fields.iter_mut().find(|(k, _)| *k == key) {
                slot.1 = value.to_owned();
            } else {
                frame.fields.push((key, value.to_owned()));
            }
        }
    });
}

/// The diagram span (the host binary's whole diagram): resets the ordinals
/// and names the diagram on every span of this thread until it closes.
pub fn diagram(graph: &str, diagram: &str) -> Span {
    if !enabled() {
        return Span::INERT;
    }
    let _ = STATE.try_with(|state| {
        if let Ok(mut state) = state.try_borrow_mut()
            && state.stack.is_empty()
        {
            state.ordinals.clear();
            state.diagram = diagram.to_owned();
        }
    });
    let span = Span::open("diagram", graph, || diagram.to_owned());
    span.set("index", diagram);
    span
}

/// The current top frame's accumulators, for a lap.
#[derive(Clone, Copy)]
struct Point {
    seq: u64,
    vakint: u64,
    child_wall: u64,
    child_vakint: u64,
}

fn point() -> Option<Point> {
    let (vakint, _) = vakint_total(&Snapshot::now());
    STATE
        .try_with(|state| {
            let state = state.try_borrow().ok()?;
            let frame = state.stack.last()?;
            Some(Point {
                seq: frame.seq,
                vakint,
                child_wall: frame.child_wall,
                child_vakint: frame.child_vakint.iter().sum(),
            })
        })
        .ok()
        .flatten()
}

/// Add `wall` ns of `leaf`, net of the vakint stages and the child spans since
/// `from`, to the top frame when it is still `from`'s.
fn charge(from: Point, leaf: &'static str, wall: u64) {
    let (vakint, _) = vakint_total(&Snapshot::now());
    let _ = STATE.try_with(|state| {
        let Ok(mut state) = state.try_borrow_mut() else {
            return;
        };
        let Some(frame) = state.stack.last_mut() else {
            return;
        };
        if frame.seq != from.seq {
            return;
        }
        let child_wall = frame.child_wall.saturating_sub(from.child_wall);
        let child_vakint = frame
            .child_vakint
            .iter()
            .sum::<u64>()
            .saturating_sub(from.child_vakint);
        let own_vakint = vakint
            .saturating_sub(from.vakint)
            .saturating_sub(child_vakint);
        let net = wall
            .saturating_sub(child_wall)
            .saturating_sub(own_vakint);
        if net == 0 {
            return;
        }
        if let Some(slot) = frame.laps.iter_mut().find(|(name, _)| *name == leaf) {
            slot.1 = slot.1.saturating_add(net);
        } else {
            frame.laps.push((leaf, net));
        }
    });
}

/// A named gammaloop scope of the current span (see [`lap`]).
#[must_use = "a lap measures until it is dropped"]
pub struct Lap {
    from: Option<(Point, Instant)>,
    leaf: &'static str,
    _not_send: PhantomData<*const ()>,
}

/// Charge the wall time until the guard drops to `leaf` of the current span,
/// net of vakint and of child spans. A lap inside a lap of the same span is
/// inert (the outer one keeps the time).
pub fn lap(leaf: &'static str) -> Lap {
    let inert = Lap {
        from: None,
        leaf,
        _not_send: PhantomData,
    };
    if !enabled() {
        return inert;
    }
    let Some(from) = point() else {
        return inert;
    };
    let opened = STATE
        .try_with(|state| {
            let Ok(mut state) = state.try_borrow_mut() else {
                return false;
            };
            match state.stack.last_mut() {
                Some(frame) if !frame.lap_open => {
                    frame.lap_open = true;
                    true
                }
                _ => false,
            }
        })
        .unwrap_or(false);
    if !opened {
        return inert;
    }
    Lap {
        from: Some((from, Instant::now())),
        leaf,
        _not_send: PhantomData,
    }
}

impl Drop for Lap {
    fn drop(&mut self) {
        let Some((from, started)) = self.from else {
            return;
        };
        let wall = started.elapsed().as_nanos() as u64;
        charge(from, self.leaf, wall);
        let _ = STATE.try_with(|state| {
            if let Ok(mut state) = state.try_borrow_mut()
                && let Some(frame) = state.stack.iter_mut().find(|frame| frame.seq == from.seq)
            {
                frame.lap_open = false;
            }
        });
    }
}

/// A running mark for code that already times consecutive phases
/// ([`super::expression::IntegrationTimings`]): [`Mark::phase`] charges the
/// interval since the previous mark to a leaf.
pub(crate) struct Mark {
    from: Option<Point>,
}

impl Mark {
    pub(crate) fn now() -> Mark {
        Mark {
            from: if enabled() { point() } else { None },
        }
    }

    /// Charge `wall` (the caller's own interval since the previous mark) to
    /// `leaf`, net of vakint and child spans; the mark moves to now.
    pub(crate) fn phase(&mut self, leaf: &'static str, wall: std::time::Duration) {
        let Some(from) = self.from else {
            return;
        };
        let open = STATE
            .try_with(|state| {
                state
                    .try_borrow()
                    .ok()
                    .and_then(|state| state.stack.last().map(|frame| frame.lap_open))
            })
            .ok()
            .flatten()
            .unwrap_or(true);
        if !open {
            charge(from, leaf, wall.as_nanos() as u64);
        }
        self.from = point();
    }
}

/// The leaf of an outer-integration phase of `IntegrationTimings`.
pub(crate) fn outer_phase_leaf(phase: &str) -> &'static str {
    match phase {
        "taylor" | "materialization" => "N.taylor",
        "tensor_projection" => "N.project",
        "algebra" => "N.contract.gl",
        "routing" | "conversion" => "N.convert",
        "backend" => "N.backend",
        "finish" => "N.finish",
        _ => "N.other",
    }
}
