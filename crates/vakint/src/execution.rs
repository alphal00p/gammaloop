//! Worker cores, Ordered weight arithmetic and cancellation of Vakint calls.
//!
//! [`VakintSettings::n_cores`](crate::VakintSettings::n_cores) is the single
//! source of every Vakint and RustRed worker count. One core runs every
//! computation on the calling thread: no thread pool is built and no thread is
//! spawned. [`OrderedWeightArithmetic`] selects how the Ordered scheduler
//! combines frontier weights; every arithmetic returns the same exact result.
//! [`Cancellation`] stops a running request at its next checkpoint.
use std::{
    cell::RefCell,
    fmt,
    num::NonZeroUsize,
    str::FromStr,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};

use symbolica::license::LicenseManager;
use symbolica::utils::resources::{self, ResourceBudget, ResourceLimitError, ResourceLimits};

use crate::{RustRedOrderedWeightArithmetic, VakintError};

/// Message of the checkpoint error a cancelled request fails with.
pub const CANCELLED_REQUEST_MESSAGE: &str = "request cancelled";

type Poll = Box<dyn Fn() -> bool + Send + Sync>;

struct Request {
    cancelled: AtomicBool,
    poll: Option<Poll>,
}

impl Request {
    /// Cancelled, or the poll asks to cancel now (then it stays cancelled).
    #[inline]
    fn tripped(&self) -> bool {
        if self.cancelled.load(Ordering::Relaxed) {
            return true;
        }
        let asked = self.poll.as_ref().is_some_and(|poll| poll());
        if asked {
            self.cancelled.store(true, Ordering::Release);
        }
        asked
    }
}

/// One entered request scope and the scopes it is nested in.
struct Scope {
    request: Cancellation,
    outer: Option<Arc<Scope>>,
}

impl Scope {
    /// This request or an enclosing one is cancelled (polls included).
    fn tripped(&self) -> bool {
        let mut scope = Some(self);
        while let Some(current) = scope {
            if current.request.0.tripped() {
                return true;
            }
            scope = current.outer.as_deref();
        }
        false
    }

    fn is_cancelled(&self) -> bool {
        let mut scope = Some(self);
        while let Some(current) = scope {
            if current.request.is_cancelled() {
                return true;
            }
            scope = current.outer.as_deref();
        }
        false
    }
}

thread_local! {
    /// The request scopes of Vakint's coarse checkpoints on this thread,
    /// innermost first.
    static REQUEST: RefCell<Option<Arc<Scope>>> = const { RefCell::new(None) };
}

/// Install `scopes` as this thread's request scopes for `operation`.
fn with_scopes<R>(scopes: Option<Arc<Scope>>, operation: impl FnOnce() -> R) -> R {
    struct Restore(Option<Arc<Scope>>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0.take();
            REQUEST.with(|request| *request.borrow_mut() = previous);
        }
    }
    let previous = REQUEST.with(|request| request.replace(scopes));
    let _restore = Restore(previous);
    operation()
}

/// The request scopes of a thread, carried to the workers it dispatches.
#[derive(Clone)]
pub(crate) struct RequestScopes(Arc<Scope>);

impl RequestScopes {
    /// Run `operation` on this (worker) thread inside the carried scopes.
    pub(crate) fn run<R>(&self, operation: impl FnOnce() -> R) -> R {
        with_scopes(Some(Arc::clone(&self.0)), operation)
    }

    /// Whether the request or an enclosing one was cancelled.
    pub(crate) fn is_cancelled(&self) -> bool {
        self.0.is_cancelled()
    }
}

/// Cooperative cancellation of one request, shared by its clones.
///
/// Two scopes, with different costs:
/// - [`Self::run`] observes the flag at EVERY resource checkpoint: it also
///   installs it as a Symbolica resource budget, so Vakint's, RustRed's and
///   Symbolica's own checkpoints fail once [`Self::cancel`] was called, on
///   the calling thread and on every worker that captures the caller's budgets
///   (Vakint's frontend workers, RustRed's Ordered scheduler and routed
///   engine). Symbolica then checks the budget at each of its checkpoints,
///   which are dense in some paths: about 2.4 s more on a ten-second
///   four-loop data load, a few percent on reductions.
/// - [`Self::run_coarse`] observes it only at Vakint's own checkpoints (per
///   data file, pack, program and expression term) and in Vakint's frontend
///   workers, at no cost elsewhere. A long native call between two of them
///   (e.g. loading the five-loop application) finishes first.
///
/// A computation without checkpoints (one CAS call, the modular lift)
/// finishes before the request stops. A cancelled request returns an error,
/// whose type depends on where the checkpoint failed (often
/// [`VakintError::ResourceLimit`] with [`CANCELLED_REQUEST_MESSAGE`]). Once
/// [`Self::is_cancelled`] holds, treat any error of the operation as the
/// cancellation. Caches never record a step that did not finish, so the
/// engine stays usable and later requests return the same results as on a
/// fresh engine. (RustRed keeps a failed one-time preparation of a
/// five-loop owner; a stopped request that left one makes the engine reload
/// the five-loop application before its next use.)
///
/// ```
/// use vakint::Cancellation;
/// use symbolica::utils::resources::checkpoint;
///
/// let cancellation = Cancellation::default();
/// assert!(cancellation.run(checkpoint).is_ok());
/// cancellation.cancel();
/// assert!(cancellation.run(checkpoint).is_err());
/// // Symbolica's own checkpoints do not see a coarse scope.
/// assert!(cancellation.run_coarse(checkpoint).is_ok());
/// ```
#[derive(Clone)]
pub struct Cancellation(Arc<Request>);

impl Default for Cancellation {
    fn default() -> Self {
        Self(Arc::new(Request {
            cancelled: AtomicBool::new(false),
            poll: None,
        }))
    }
}

impl fmt::Debug for Cancellation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Cancellation")
            .field("cancelled", &self.is_cancelled())
            .field("polled", &self.0.poll.is_some())
            .finish()
    }
}

impl Cancellation {
    /// A cancellation that also asks `poll` at each observing checkpoint
    /// whether to cancel (`true` cancels, as [`Self::cancel`] would).
    ///
    /// `poll` runs on the checkpoint's thread, also on dispatched workers,
    /// up to millions of times per second under [`Self::run`]: keep it to a
    /// few atomic loads in the common case. Like every Symbolica resource
    /// monitor it must not enter a budget scope or call checkpoints.
    pub fn with_poll(poll: impl Fn() -> bool + Send + Sync + 'static) -> Self {
        Self(Arc::new(Request {
            cancelled: AtomicBool::new(false),
            poll: Some(Box::new(poll)),
        }))
    }

    /// Request cancellation; operations running under [`Self::run`] or
    /// [`Self::run_coarse`] fail at their next observing checkpoint.
    pub fn cancel(&self) {
        self.0.cancelled.store(true, Ordering::Release);
    }

    /// Whether the request was cancelled (by [`Self::cancel`] or the poll).
    pub fn is_cancelled(&self) -> bool {
        self.0.cancelled.load(Ordering::Acquire)
    }

    /// Run `operation` with every resource checkpoint observing this flag.
    /// Budgets installed inside it (e.g. the engine's
    /// [`ResourceLimits`](crate::ResourceLimits)) keep their meaning.
    pub fn run<R>(&self, operation: impl FnOnce() -> R) -> R {
        let request = Arc::clone(&self.0);
        let budget = ResourceBudget::with_memory_monitor(
            // Only the monitor's verdict matters: it reports zero bytes while
            // the request runs, which no allocation check can exceed.
            ResourceLimits {
                max_process_rss_bytes: Some(u64::MAX),
                max_elapsed: None,
            },
            move || {
                if request.tripped() {
                    Err(CANCELLED_REQUEST_MESSAGE.to_owned())
                } else {
                    Ok(0)
                }
            },
        );
        self.run_coarse(|| budget.run(operation))
    }

    /// Run `operation` with Vakint's own (coarse) checkpoints observing this
    /// flag, at no cost to Symbolica's and RustRed's checkpoints.
    ///
    /// Scopes nest like Symbolica budgets: an inner request never hides an
    /// outer one, so cancelling either stops the inner operation.
    pub fn run_coarse<R>(&self, operation: impl FnOnce() -> R) -> R {
        let outer = REQUEST.with(|request| request.borrow().clone());
        let scope = Arc::new(Scope {
            request: self.clone(),
            outer,
        });
        with_scopes(Some(scope), operation)
    }

    /// The request scopes of this thread, for workers it dispatches.
    pub(crate) fn current() -> Option<RequestScopes> {
        REQUEST.with(|request| request.borrow().clone().map(RequestScopes))
    }
}

/// Vakint's own checkpoint: every Symbolica budget of this thread, then the
/// coarse request scope ([`Cancellation::run_coarse`]).
pub(crate) fn checkpoint() -> Result<(), ResourceLimitError> {
    resources::checkpoint()?;
    let cancelled = REQUEST.with(|request| {
        request
            .borrow()
            .as_ref()
            .is_some_and(|scope| scope.tripped())
    });
    if cancelled {
        Err(ResourceLimitError::MonitorUnavailable {
            resource: "request",
            message: CANCELLED_REQUEST_MESSAGE.to_owned(),
        })
    } else {
        Ok(())
    }
}

/// Upper bound of [`default_n_cores`].
pub const MAX_DEFAULT_N_CORES: usize = 10;

/// Ordered weight arithmetic selected by [`crate::VakintSettings::default`]
/// and the Python bindings, in the text form of [`OrderedWeightArithmetic`].
/// Changing this constant changes the default of every Vakint entry point.
/// `modular_polynomial` is the fastest arithmetic up to four loops; five-loop
/// requests (vacuum and p) default to `modular_sampled`, see
/// `five_loop_ordered_weight_arithmetic`; every arithmetic returns the same
/// exact result.
pub const DEFAULT_ORDERED_WEIGHT_ARITHMETIC: &str = "modular_polynomial";

/// Default [`crate::VakintSettings::n_cores`]: the cores available to this
/// process (affinity masks and cgroup quotas included), at most
/// [`MAX_DEFAULT_N_CORES`]. Unlicensed Symbolica is restricted to one thread,
/// so the default is one core without a Symbolica license.
pub fn default_n_cores() -> usize {
    let cores = std::thread::available_parallelism()
        .map_or(1, NonZeroUsize::get)
        .min(MAX_DEFAULT_N_CORES);
    // The license is only consulted for a multicore budget.
    if cores > 1 {
        LicenseManager::max_threads(cores).max(1)
    } else {
        1
    }
}

/// Check one worker count against the cores available to this process and
/// the Symbolica license, as RustRed's own executors do before building a pool.
pub(crate) fn validate_n_cores(n_cores: usize) -> Result<(), VakintError> {
    rustred::campaign::ParallelExecution::preflight_requested_core_budget(n_cores).map_err(
        |error| VakintError::InvalidNCores {
            n_cores,
            reason: error.to_string(),
        },
    )
}

/// How Vakint's Ordered RustRed requests combine frontier weights.
///
/// One of RustRed's three modes ([`RustRedOrderedWeightArithmetic`]):
/// `exact` keeps normalized rational functions in `d` on the frontier;
/// `modular_sampled` records the weight-independent trace once, replays it at
/// sample points over 64-bit prime fields and interpolates every terminal
/// coefficient per prime; `modular_polynomial` runs scheduler passes with
/// polynomial weights over prime fields. Both modular modes lift their prime
/// images over the rationals (CRT and rational reconstruction) and verify them
/// on fresh primes, so every mode returns the same exact result.
///
/// The text form is RustRed's canonical syntax
/// `exact | modular_sampled[:key=value,...] | modular_polynomial[:key=value,...]`
/// (`modular` is an alias of `modular_sampled`); the keys are the fields of
/// RustRed's `ModularCommon` and of the mode's options. A modular mode runs on
/// [`crate::VakintSettings::n_cores`] workers unless the selection fixes
/// `threads`.
///
/// Equality compares the selection (mode, options, fixed workers); whether it
/// is the default selection ([`Self::is_default_selection`]) is not part of it.
#[derive(Clone, Copy, Debug)]
pub struct OrderedWeightArithmetic {
    /// Mode and options with the worker count normalized to one.
    mode: RustRedOrderedWeightArithmetic,
    /// Worker count fixed by the selection instead of `n_cores`.
    threads: Option<usize>,
    /// Made by [`Default`] (no explicit choice of the mode).
    by_default: bool,
}

impl PartialEq for OrderedWeightArithmetic {
    fn eq(&self, other: &Self) -> bool {
        self.mode == other.mode && self.threads == other.threads
    }
}

impl Eq for OrderedWeightArithmetic {}

impl OrderedWeightArithmetic {
    /// Select `arithmetic`; a modular mode runs on `n_cores` workers, whatever
    /// the worker count of `arithmetic`.
    pub fn new(arithmetic: RustRedOrderedWeightArithmetic) -> Self {
        Self {
            mode: arithmetic.with_threads(1),
            threads: None,
            by_default: false,
        }
    }

    /// Exact normalized rational-function weights.
    pub fn exact() -> Self {
        Self::new(RustRedOrderedWeightArithmetic::Exact)
    }

    /// RustRed's sampled modular arithmetic with default options.
    pub fn modular_sampled() -> Self {
        Self::new(RustRedOrderedWeightArithmetic::modular_sampled())
    }

    /// RustRed's polynomial modular arithmetic with default options.
    pub fn modular_polynomial() -> Self {
        Self::new(RustRedOrderedWeightArithmetic::modular_polynomial())
    }

    /// Fix the worker count of a modular mode instead of following `n_cores`;
    /// `exact` has no workers and is unchanged.
    pub fn with_threads(self, threads: usize) -> Self {
        Self {
            threads: self.mode.common().map(|_| threads),
            ..self
        }
    }

    /// Mode name: `exact`, `modular_sampled` or `modular_polynomial`.
    pub fn name(&self) -> &'static str {
        self.mode.name()
    }

    /// Worker count fixed by the selection, if any.
    pub fn threads(&self) -> Option<usize> {
        self.threads
    }

    /// Whether the mode is the default one ([`Default`], also with workers
    /// fixed by [`Self::with_threads`]) rather than an explicit choice: an
    /// explicitly selected `modular_polynomial` is not the default selection,
    /// although it equals it. Five-loop requests use `modular_sampled` for the
    /// default selection.
    pub fn is_default_selection(&self) -> bool {
        self.by_default
    }

    /// RustRed's arithmetic of a call with `n_cores` workers.
    pub fn resolve(&self, n_cores: usize) -> RustRedOrderedWeightArithmetic {
        self.mode.with_threads(self.threads.unwrap_or(n_cores))
    }
}

impl Default for OrderedWeightArithmetic {
    fn default() -> Self {
        let selection: Self = DEFAULT_ORDERED_WEIGHT_ARITHMETIC
            .parse()
            .expect("the default ordered weight arithmetic is valid");
        Self {
            by_default: true,
            ..selection
        }
    }
}

impl From<RustRedOrderedWeightArithmetic> for OrderedWeightArithmetic {
    fn from(arithmetic: RustRedOrderedWeightArithmetic) -> Self {
        Self::new(arithmetic)
    }
}

impl FromStr for OrderedWeightArithmetic {
    type Err = VakintError;

    /// Parse RustRed's canonical text; a `threads` option fixes the worker
    /// count, otherwise the selection follows `n_cores`.
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let arithmetic = text
            .parse::<RustRedOrderedWeightArithmetic>()
            .map_err(|error| VakintError::InvalidOrderedWeightArithmetic(error.to_string()))?;
        // RustRed accepted the syntax, so its options are `key=value` items
        // with unique keys.
        let fixes_threads = text.split_once(':').is_some_and(|(_, options)| {
            options.split(',').any(|item| {
                item.split_once('=')
                    .is_some_and(|(key, _)| key.trim().eq_ignore_ascii_case("threads"))
            })
        });
        let selection = Self::new(arithmetic);
        Ok(if fixes_threads {
            selection.with_threads(arithmetic.threads())
        } else {
            selection
        })
    }
}

impl fmt::Display for OrderedWeightArithmetic {
    /// RustRed's canonical text of the options that differ from their
    /// defaults, led by `threads` when the selection fixes it.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let text = self.mode.to_string();
        match (self.threads, text.split_once(':')) {
            (None, _) => f.write_str(&text),
            (Some(threads), Some((mode, options))) => {
                write!(f, "{mode}:threads={threads},{options}")
            }
            (Some(threads), None) => write!(f, "{text}:threads={threads}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{RustRedModularCommon, RustRedSampledOptions};
    use std::time::{Duration, Instant};
    use symbolica::utils::resources::{self, ResourceLimitError};

    /// Checkpoint until one fails; the time from `cancelled_at` to that failure.
    fn checkpoint_until_failure(cancelled_at: &std::sync::OnceLock<Instant>) -> Duration {
        loop {
            if let Err(error) = resources::checkpoint() {
                assert!(
                    matches!(
                        &error,
                        ResourceLimitError::MonitorUnavailable { message, .. }
                            if message == CANCELLED_REQUEST_MESSAGE
                    ),
                    "{error}"
                );
                let cancelled = cancelled_at.get().expect("failure only after cancel()");
                return cancelled.elapsed();
            }
            std::hint::spin_loop();
        }
    }

    #[test]
    fn cancellation_fails_next_checkpoint_on_caller_and_captured_workers() {
        let cancellation = Cancellation::default();
        let cancelled_at = std::sync::OnceLock::new();
        let latencies = std::thread::scope(|scope| {
            let canceller = cancellation.clone();
            let stamp = &cancelled_at;
            scope.spawn(move || {
                std::thread::sleep(Duration::from_millis(50));
                stamp.set(Instant::now()).unwrap();
                canceller.cancel();
            });
            cancellation.run(|| {
                // Four frontend workers dispatched from a thread that entered
                // the caller's captured budgets (as every worker pool does),
                // plus the caller itself.
                let captured = ResourceBudget::capture_current();
                let workers = std::thread::scope(|inner| {
                    let handle = inner.spawn(|| {
                        captured.run_captured(|| {
                            crate::rustred_evaluation::parallel::map_ordered(4, &[(); 4], |_, _| {
                                checkpoint_until_failure(&cancelled_at)
                            })
                        })
                    });
                    let caller = checkpoint_until_failure(&cancelled_at);
                    let mut latencies = handle.join().unwrap();
                    latencies.push(caller);
                    latencies
                });
                assert!(cancellation.is_cancelled());
                workers
            })
        });
        assert_eq!(latencies.len(), 5);
        for latency in latencies {
            assert!(latency < Duration::from_millis(50), "{latency:?}");
        }
        // Outside the request nothing observes the flag.
        resources::checkpoint().unwrap();
    }

    #[test]
    fn polling_cancels_at_the_first_true_and_stays_cancelled() {
        let calls = Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let counted = Arc::clone(&calls);
        let cancellation =
            Cancellation::with_poll(move || counted.fetch_add(1, Ordering::Relaxed) == 2);
        cancellation.run(|| {
            resources::checkpoint().unwrap();
            resources::checkpoint().unwrap();
            assert!(!cancellation.is_cancelled());
            assert!(resources::checkpoint().is_err());
            assert!(cancellation.is_cancelled());
            // Cancelled: later checkpoints fail without asking again.
            assert!(resources::checkpoint().is_err());
            assert!(checkpoint().is_err());
        });
        assert_eq!(calls.load(Ordering::Relaxed), 3);
    }

    #[test]
    fn coarse_scopes_nest_and_are_seen_by_vakint_checkpoints_and_frontend_workers_only() {
        let cancellation = Cancellation::default();
        cancellation.run_coarse(|| {
            checkpoint().unwrap();
            assert!(Cancellation::current().is_some());
            cancellation.cancel();
            // Symbolica's own checkpoints carry no budget for it.
            resources::checkpoint().unwrap();
            assert!(matches!(
                checkpoint(),
                Err(ResourceLimitError::MonitorUnavailable {
                    resource: "request",
                    ..
                })
            ));
            // Frontend workers inherit the scope.
            let seen = crate::rustred_evaluation::parallel::map_ordered(3, &[(); 3], |_, _| {
                checkpoint().is_err()
            });
            assert_eq!(seen, [true, true, true]);
            // A nested scope never hides the cancelled outer one, and the
            // outer one is restored after it.
            Cancellation::default().run_coarse(|| {
                assert!(checkpoint().is_err());
                assert!(Cancellation::current().unwrap().is_cancelled());
            });
            assert!(checkpoint().is_err());
        });
        assert!(Cancellation::current().is_none());
        checkpoint().unwrap();
        // Cancelling only the inner scope stops the inner operation only.
        let outer = Cancellation::default();
        outer.run_coarse(|| {
            let inner = Cancellation::default();
            inner.run_coarse(|| {
                checkpoint().unwrap();
                inner.cancel();
                assert!(checkpoint().is_err());
                // Workers carry the whole chain.
                let seen = crate::rustred_evaluation::parallel::map_ordered(2, &[(); 2], |_, _| {
                    checkpoint().is_err()
                });
                assert_eq!(seen, [true, true]);
            });
            checkpoint().unwrap();
            // An outer cancellation reaches a scope entered before it.
            Cancellation::default().run_coarse(|| {
                checkpoint().unwrap();
                outer.cancel();
                assert!(checkpoint().is_err());
            });
        });
    }

    #[test]
    fn cancellation_nests_with_user_limits() {
        let cancellation = Cancellation::default();
        let memory = |limit: u64, observed: u64| {
            ResourceBudget::with_memory_monitor(
                ResourceLimits {
                    max_process_rss_bytes: Some(limit),
                    ..Default::default()
                },
                move || Ok(observed),
            )
        };
        cancellation.run(|| {
            // Not cancelled: no headroom check can trip the cancellation budget.
            resources::check_allocation(u64::MAX - 1).unwrap();
            resources::check_allocation(u64::MAX).unwrap();
            // User budgets inside the request report their own errors.
            assert!(matches!(
                memory(5, 6).run(resources::checkpoint),
                Err(ResourceLimitError::Memory { .. })
            ));
            assert!(matches!(
                ResourceBudget::new(ResourceLimits {
                    max_elapsed: Some(Duration::ZERO),
                    ..Default::default()
                })
                .run(resources::checkpoint),
                Err(ResourceLimitError::Elapsed { .. })
            ));
            memory(100, 1).run(resources::checkpoint).unwrap();
            // Once cancelled, the (outer) cancellation fails first.
            cancellation.cancel();
            let error = memory(5, 6).run(resources::checkpoint).unwrap_err();
            assert_eq!(
                error,
                ResourceLimitError::MonitorUnavailable {
                    resource: "process RSS",
                    message: CANCELLED_REQUEST_MESSAGE.into(),
                }
            );
            assert!(matches!(
                VakintError::from(error),
                VakintError::ResourceLimit(_)
            ));
        });
        // Clones share the flag; a fresh cancellation does not.
        assert!(cancellation.clone().is_cancelled());
        assert!(!Cancellation::default().is_cancelled());
        Cancellation::default().run(resources::checkpoint).unwrap();
    }

    #[test]
    fn default_n_cores_is_bounded_by_the_available_cores() {
        let available = std::thread::available_parallelism().map_or(1, NonZeroUsize::get);
        let n_cores = default_n_cores();
        assert!((1..=available.min(MAX_DEFAULT_N_CORES)).contains(&n_cores));
        validate_n_cores(n_cores).unwrap();
    }

    #[test]
    fn invalid_core_counts_are_rejected_with_their_reason() {
        let zero = validate_n_cores(0).unwrap_err().to_string();
        assert!(
            zero.contains("n_cores = 0") && zero.contains("positive"),
            "{zero}"
        );
        let available = std::thread::available_parallelism().map_or(1, NonZeroUsize::get);
        let excessive = validate_n_cores(available + 1).unwrap_err().to_string();
        assert!(excessive.contains("exceeds"), "{excessive}");
        validate_n_cores(1).unwrap();
    }

    #[test]
    fn the_default_arithmetic_is_the_default_constant() {
        assert_eq!(
            OrderedWeightArithmetic::default(),
            DEFAULT_ORDERED_WEIGHT_ARITHMETIC.parse().unwrap()
        );
        assert_eq!(OrderedWeightArithmetic::default().threads(), None);
        assert_eq!(
            OrderedWeightArithmetic::default(),
            OrderedWeightArithmetic::modular_polynomial()
        );
        // Equal, but only the default is the default selection.
        assert!(OrderedWeightArithmetic::default().is_default_selection());
        assert!(
            OrderedWeightArithmetic::default()
                .with_threads(2)
                .is_default_selection()
        );
        assert!(!OrderedWeightArithmetic::modular_polynomial().is_default_selection());
        let parsed: OrderedWeightArithmetic = DEFAULT_ORDERED_WEIGHT_ARITHMETIC.parse().unwrap();
        assert!(!parsed.is_default_selection());
    }

    #[test]
    fn modular_workers_follow_n_cores_unless_fixed() {
        let sampled: OrderedWeightArithmetic = "modular_sampled".parse().unwrap();
        assert_eq!(sampled, OrderedWeightArithmetic::modular_sampled());
        assert_eq!(sampled.threads(), None);
        assert_eq!(sampled.resolve(1).threads(), 1);
        assert_eq!(sampled.resolve(6).threads(), 6);
        let fixed: OrderedWeightArithmetic = "modular:batch=64,threads=1".parse().unwrap();
        assert_eq!(fixed.threads(), Some(1));
        assert_eq!(
            fixed.resolve(6),
            RustRedOrderedWeightArithmetic::ModularSampled(
                RustRedModularCommon {
                    threads: 1,
                    ..Default::default()
                },
                RustRedSampledOptions {
                    batch_points: 64,
                    ..Default::default()
                }
            )
        );
        let exact = OrderedWeightArithmetic::exact().with_threads(4);
        assert_eq!(exact.threads(), None);
        assert_eq!(exact.resolve(8), RustRedOrderedWeightArithmetic::Exact);
        // Every worker count resolves to the same mode and options.
        let polynomial = OrderedWeightArithmetic::modular_polynomial();
        assert_eq!(
            polynomial.resolve(3).with_threads(1),
            RustRedOrderedWeightArithmetic::modular_polynomial()
        );
    }

    #[test]
    fn text_form_round_trips_and_keeps_the_worker_choice() {
        for (text, canonical) in [
            ("exact", "exact"),
            (" Exact ", "exact"),
            ("modular", "modular_sampled"),
            ("modular_sampled:threads=4", "modular_sampled:threads=4"),
            ("modular_sampled:threads=1", "modular_sampled:threads=1"),
            (
                "modular:cross_check=1,threads=16,batch=64",
                "modular_sampled:threads=16,cross_check_exact=true,batch_points=64",
            ),
            (
                "modular_polynomial:denominators=dense",
                "modular_polynomial:denominators=dense",
            ),
            (
                "modular_polynomial:lanes=2,THREADS=3",
                "modular_polynomial:threads=3,lanes=2",
            ),
        ] {
            let selection: OrderedWeightArithmetic = text.parse().unwrap();
            assert_eq!(selection.to_string(), canonical, "{text}");
            assert_eq!(
                canonical.parse::<OrderedWeightArithmetic>().unwrap(),
                selection
            );
        }
    }

    #[test]
    fn invalid_text_is_rejected() {
        for text in [
            "",
            "dfs",
            "exact:threads=2",
            "modular:threads",
            "modular:threads=0",
            "modular:threads=1,threads=2",
            "modular_sampled:lanes=2",
            "modular_polynomial:batch_points=8",
        ] {
            let error = text.parse::<OrderedWeightArithmetic>().unwrap_err();
            assert!(
                matches!(error, VakintError::InvalidOrderedWeightArithmetic(_)),
                "{text}: {error}"
            );
        }
    }
}
