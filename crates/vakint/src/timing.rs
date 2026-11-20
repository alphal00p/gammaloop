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
//! The counters only advance on the thread that holds the guards: callers read
//! a [`Snapshot`] before and after a unit on the same thread. Callers measure
//! wall time, so work a reduction hands to helper threads is included in the
//! caller's reduction seconds (it waits for them). Cost: two clock reads per
//! scope transition, no allocation.

use std::cell::Cell;
use std::marker::PhantomData;
use std::time::Instant;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    None,
    Lower,
    Reduce,
}

struct State {
    reduce_ns: Cell<u64>,
    lower_ns: Cell<u64>,
    calls: Cell<u64>,
    mode: Cell<Kind>,
    since: Cell<Option<Instant>>,
    reduce_depth: Cell<u32>,
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
        }
    };
}

fn nanos(from: Option<Instant>, to: Instant) -> u64 {
    from.map(|from| to.saturating_duration_since(from).as_nanos() as u64)
        .unwrap_or(0)
}

impl State {
    /// Add the time since the last transition to the active kind's counter.
    fn flush(&self, now: Instant) {
        let delta = nanos(self.since.get(), now);
        match self.mode.get() {
            Kind::Reduce => self.reduce_ns.set(self.reduce_ns.get().saturating_add(delta)),
            Kind::Lower => self.lower_ns.set(self.lower_ns.get().saturating_add(delta)),
            Kind::None => {}
        }
        self.since.set(Some(now));
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
                    state.reduce_depth.set(state.reduce_depth.get().saturating_add(1));
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
                state.reduce_depth.set(state.reduce_depth.get().saturating_sub(1));
            }
            state.mode.set(previous);
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

/// The cumulative counters of this thread, including the running scope.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Snapshot {
    reduce_ns: u64,
    lower_ns: u64,
    calls: u64,
}

impl Snapshot {
    pub fn now() -> Self {
        STATE
            .try_with(|state| {
                let mut snapshot = Snapshot {
                    reduce_ns: state.reduce_ns.get(),
                    lower_ns: state.lower_ns.get(),
                    calls: state.calls.get(),
                };
                let running = nanos(state.since.get(), Instant::now());
                match state.mode.get() {
                    Kind::Reduce => snapshot.reduce_ns += running,
                    Kind::Lower => snapshot.lower_ns += running,
                    Kind::None => {}
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
        Delta {
            reduce_seconds: self.reduce_ns.saturating_sub(earlier.reduce_ns) as f64 * 1e-9,
            lower_seconds: self.lower_ns.saturating_sub(earlier.lower_ns) as f64 * 1e-9,
            reduce_calls: self.calls.saturating_sub(earlier.calls),
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
}

/// The counters of one interval.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Delta {
    pub reduce_seconds: f64,
    pub lower_seconds: f64,
    pub reduce_calls: u64,
}

impl Delta {
    /// ` reduce_s=<s> lower_s=<s>` (leading space), for the log lines.
    pub fn fields(&self) -> String {
        format!(
            " reduce_s={:.3} lower_s={:.3}",
            self.reduce_seconds, self.lower_seconds
        )
    }
}

impl std::ops::AddAssign for Delta {
    fn add_assign(&mut self, other: Delta) {
        self.reduce_seconds += other.reduce_seconds;
        self.lower_seconds += other.lower_seconds;
        self.reduce_calls += other.reduce_calls;
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
}
