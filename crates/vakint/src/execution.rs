//! Worker cores and Ordered weight arithmetic of Vakint calls.
//!
//! [`VakintSettings::n_cores`](crate::VakintSettings::n_cores) is the single
//! source of every Vakint and RustRed worker count. One core runs every
//! computation on the calling thread: no thread pool is built and no thread is
//! spawned. [`OrderedWeightArithmetic`] selects how the Ordered scheduler
//! combines frontier weights; every arithmetic returns the same exact result.
use std::{fmt, num::NonZeroUsize, str::FromStr};

use symbolica::license::LicenseManager;

use crate::{RustRedOrderedWeightArithmetic, VakintError};

/// Upper bound of [`default_n_cores`].
pub const MAX_DEFAULT_N_CORES: usize = 10;

/// Ordered weight arithmetic selected by [`crate::VakintSettings::default`]
/// and the Python bindings, in the text form of [`OrderedWeightArithmetic`].
/// Changing this constant changes the default of every Vakint entry point.
/// `modular_polynomial` is the fastest arithmetic (at five loops about 28 times
/// less replay CPU than `modular_sampled`); every arithmetic returns the same
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
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OrderedWeightArithmetic {
    /// Mode and options with the worker count normalized to one.
    mode: RustRedOrderedWeightArithmetic,
    /// Worker count fixed by the selection instead of `n_cores`.
    threads: Option<usize>,
}

impl OrderedWeightArithmetic {
    /// Select `arithmetic`; a modular mode runs on `n_cores` workers, whatever
    /// the worker count of `arithmetic`.
    pub fn new(arithmetic: RustRedOrderedWeightArithmetic) -> Self {
        Self {
            mode: arithmetic.with_threads(1),
            threads: None,
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

    /// RustRed's arithmetic of a call with `n_cores` workers.
    pub fn resolve(&self, n_cores: usize) -> RustRedOrderedWeightArithmetic {
        self.mode.with_threads(self.threads.unwrap_or(n_cores))
    }
}

impl Default for OrderedWeightArithmetic {
    fn default() -> Self {
        DEFAULT_ORDERED_WEIGHT_ARITHMETIC
            .parse()
            .expect("the default ordered weight arithmetic is valid")
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
