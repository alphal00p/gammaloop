//! Wall-time split of one evaluation for profiling (`VAKINT_P_PHASES=1`):
//! each timed phase prints `vakint-phase <name> <seconds>` on stderr. Nothing
//! is measured unless the variable is set; values never depend on it.

use std::sync::OnceLock;
use std::time::{Duration, Instant};

/// Whether `VAKINT_P_PHASES` is set (read once per process).
pub(crate) fn enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var_os("VAKINT_P_PHASES").is_some_and(|value| value != "0"))
}

/// Whether `VAKINT_P_PHASES=2`: one line per p reduction request too.
pub(crate) fn per_request() -> bool {
    static PER_REQUEST: OnceLock<bool> = OnceLock::new();
    *PER_REQUEST.get_or_init(|| std::env::var("VAKINT_P_PHASES").is_ok_and(|value| value == "2"))
}

/// The start of a timed phase; `None` when profiling is off.
pub(crate) fn start() -> Option<Instant> {
    enabled().then(Instant::now)
}

/// Print the wall time of `name` since `started`.
pub(crate) fn end(name: &str, started: Option<Instant>) {
    if let Some(started) = started {
        report(name, started.elapsed());
    }
}

/// Print an accumulated duration (and an optional count).
pub(crate) fn report(name: &str, elapsed: Duration) {
    if enabled() {
        eprintln!("vakint-phase {name} {:.6}", elapsed.as_secs_f64());
    }
}

/// Time `work` as `name`.
pub(crate) fn time<T>(name: &str, work: impl FnOnce() -> T) -> T {
    let started = start();
    let value = work();
    end(name, started);
    value
}

/// An accumulator for a phase repeated many times (per channel).
#[derive(Default)]
pub(crate) struct Accumulator {
    pub total: Duration,
    pub count: usize,
}

impl Accumulator {
    pub(crate) fn time<T>(&mut self, work: impl FnOnce() -> T) -> T {
        if !enabled() {
            return work();
        }
        let started = Instant::now();
        let value = work();
        self.total += started.elapsed();
        self.count += 1;
        value
    }

    pub(crate) fn report(&self, name: &str) {
        if enabled() {
            eprintln!(
                "vakint-phase {name} {:.6} count={}",
                self.total.as_secs_f64(),
                self.count
            );
        }
    }
}
