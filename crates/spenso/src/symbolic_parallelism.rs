use std::sync::{
    OnceLock,
    atomic::{AtomicU8, Ordering},
};

use symbolica::LicenseManager;

/// Policy used to configure Rayon for operations that manipulate Symbolica atoms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SymbolicParallelism {
    /// Permit Rayon when licensed and use workload heuristics where available.
    Auto,
    /// Keep symbolic operations on the calling thread.
    Serial,
    /// Force Rayon without `Auto`'s Symbolica license safety check.
    Parallel,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
enum ResolvedSymbolicParallelism {
    Serial,
    Adaptive,
    Parallel,
}

impl ResolvedSymbolicParallelism {
    fn from_raw(raw: u8) -> Self {
        match raw {
            raw if raw == Self::Serial as u8 => Self::Serial,
            raw if raw == Self::Adaptive as u8 => Self::Adaptive,
            raw if raw == Self::Parallel as u8 => Self::Parallel,
            _ => unreachable!("invalid cached symbolic parallelism state"),
        }
    }

    fn rayon_enabled(self) -> bool {
        !matches!(self, Self::Serial)
    }

    fn adapts_to_workload(self) -> bool {
        matches!(self, Self::Adaptive)
    }

    fn rayon_enabled_for(self, is_profitable: impl FnOnce() -> bool) -> bool {
        match self {
            Self::Serial => false,
            Self::Adaptive => is_profitable(),
            Self::Parallel => true,
        }
    }
}

static SYMBOLICA_RAYON_POLICY: OnceLock<AtomicU8> = OnceLock::new();

fn cached_setting() -> &'static AtomicU8 {
    SYMBOLICA_RAYON_POLICY.get_or_init(|| {
        AtomicU8::new(SymbolicParallelism::Auto.resolve_with(LicenseManager::is_licensed) as u8)
    })
}

fn resolved_setting() -> ResolvedSymbolicParallelism {
    ResolvedSymbolicParallelism::from_raw(cached_setting().load(Ordering::Relaxed))
}

impl SymbolicParallelism {
    fn resolve_with(self, is_licensed: impl FnOnce() -> bool) -> ResolvedSymbolicParallelism {
        match self {
            Self::Auto if is_licensed() => ResolvedSymbolicParallelism::Adaptive,
            Self::Auto | Self::Serial => ResolvedSymbolicParallelism::Serial,
            Self::Parallel => ResolvedSymbolicParallelism::Parallel,
        }
    }

    /// Decide whether one symbolic operation should use Rayon.
    ///
    /// The profitability candidate is evaluated only by licensed [`Self::Auto`].
    /// Explicit serial and parallel policies override it without evaluating it.
    pub(crate) fn rayon_enabled_for(is_profitable: impl FnOnce() -> bool) -> bool {
        resolved_setting().rayon_enabled_for(is_profitable)
    }

    /// Return whether the current policy needs a profitability estimate.
    pub(crate) fn adapts_to_workload() -> bool {
        resolved_setting().adapts_to_workload()
    }
}

/// Configure whether operations involving Symbolica atoms may use Rayon.
///
/// [`SymbolicParallelism::Auto`] queries the Symbolica license exactly once,
/// when this function is called. A licensed configuration permits Rayon;
/// operations with a workload model may still choose serial execution, while
/// other operations retain their existing parallel behavior. An unlicensed
/// configuration remains serial. Later operations only read the cached state.
/// Configure this before starting tensor operations; changing it concurrently
/// with an active operation is not supported.
pub fn set_symbolica_rayon_enabled(policy: SymbolicParallelism) {
    set_symbolica_rayon_enabled_with(policy, LicenseManager::is_licensed);
}

fn set_symbolica_rayon_enabled_with(
    policy: SymbolicParallelism,
    is_licensed: impl FnOnce() -> bool,
) {
    let resolved = policy.resolve_with(is_licensed);
    SYMBOLICA_RAYON_POLICY
        .get_or_init(|| AtomicU8::new(resolved as u8))
        .store(resolved as u8, Ordering::Relaxed);
}

#[cfg(test)]
pub(crate) struct SymbolicParallelismTestGuard {
    previous: ResolvedSymbolicParallelism,
    _lock: std::sync::MutexGuard<'static, ()>,
}

#[cfg(test)]
static SYMBOLIC_PARALLELISM_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

#[cfg(test)]
impl Drop for SymbolicParallelismTestGuard {
    fn drop(&mut self) {
        cached_setting().store(self.previous as u8, Ordering::Relaxed);
    }
}

#[cfg(test)]
pub(crate) fn scoped_symbolica_rayon_setting_for_test(
    policy: SymbolicParallelism,
    is_licensed: impl FnOnce() -> bool,
) -> SymbolicParallelismTestGuard {
    let lock = SYMBOLIC_PARALLELISM_TEST_LOCK
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let previous = resolved_setting();
    set_symbolica_rayon_enabled_with(policy, is_licensed);
    SymbolicParallelismTestGuard {
        previous,
        _lock: lock,
    }
}

/// Return whether the resolved setting permits symbolic Rayon work.
///
/// If no policy has been configured yet, this performs the default
/// [`SymbolicParallelism::Auto`] initialization and queries the license once.
/// Subsequent reads only load the cached state. For licensed `Auto`, `true`
/// means that operations may use Rayon; operations with a workload model may
/// still choose serial execution.
pub fn symbolica_rayon_enabled() -> bool {
    resolved_setting().rayon_enabled()
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;

    use super::*;

    #[test]
    fn auto_checks_the_license_once_when_resolved() {
        let calls = Cell::new(0);
        let resolved = SymbolicParallelism::Auto.resolve_with(|| {
            calls.set(calls.get() + 1);
            true
        });

        assert_eq!(resolved, ResolvedSymbolicParallelism::Adaptive);
        assert_eq!(calls.get(), 1);
    }

    #[test]
    fn unlicensed_auto_resolves_to_serial() {
        assert_eq!(
            SymbolicParallelism::Auto.resolve_with(|| false),
            ResolvedSymbolicParallelism::Serial
        );
    }

    #[test]
    fn explicit_policies_do_not_check_the_license() {
        for (policy, expected) in [
            (
                SymbolicParallelism::Serial,
                ResolvedSymbolicParallelism::Serial,
            ),
            (
                SymbolicParallelism::Parallel,
                ResolvedSymbolicParallelism::Parallel,
            ),
        ] {
            let calls = Cell::new(0);
            let resolved = policy.resolve_with(|| {
                calls.set(calls.get() + 1);
                false
            });

            assert_eq!(resolved, expected);
            assert_eq!(calls.get(), 0);
        }
    }

    #[test]
    fn resolved_policy_applies_the_profitability_candidate_only_in_adaptive_mode() {
        for (resolved, profitable, expected) in [
            (ResolvedSymbolicParallelism::Serial, false, false),
            (ResolvedSymbolicParallelism::Serial, true, false),
            (ResolvedSymbolicParallelism::Adaptive, false, false),
            (ResolvedSymbolicParallelism::Adaptive, true, true),
            (ResolvedSymbolicParallelism::Parallel, false, true),
            (ResolvedSymbolicParallelism::Parallel, true, true),
        ] {
            assert_eq!(resolved.rayon_enabled_for(|| profitable), expected);
        }
    }

    #[test]
    fn forced_policies_do_not_evaluate_the_profitability_candidate() {
        for (resolved, expected) in [
            (ResolvedSymbolicParallelism::Serial, false),
            (ResolvedSymbolicParallelism::Parallel, true),
        ] {
            let calls = Cell::new(0);
            assert_eq!(
                resolved.rayon_enabled_for(|| {
                    calls.set(calls.get() + 1);
                    false
                }),
                expected
            );

            assert_eq!(calls.get(), 0);
        }
    }

    #[test]
    fn setter_caches_safe_serial_policy() {
        let _guard = scoped_symbolica_rayon_setting_for_test(SymbolicParallelism::Serial, || false);
        assert!(!symbolica_rayon_enabled());

        set_symbolica_rayon_enabled_with(SymbolicParallelism::Auto, || false);
        assert!(!symbolica_rayon_enabled());
    }

    #[test]
    fn resolved_modes_report_rayon_capability() {
        assert!(!ResolvedSymbolicParallelism::Serial.rayon_enabled());
        assert!(ResolvedSymbolicParallelism::Adaptive.rayon_enabled());
        assert!(ResolvedSymbolicParallelism::Parallel.rayon_enabled());
    }

    #[test]
    fn only_licensed_auto_requests_workload_estimates() {
        for (policy, licensed, expected) in [
            (SymbolicParallelism::Serial, false, false),
            (SymbolicParallelism::Auto, false, false),
            (SymbolicParallelism::Auto, true, true),
            (SymbolicParallelism::Parallel, false, false),
        ] {
            let resolved = policy.resolve_with(|| licensed);
            assert_eq!(resolved.adapts_to_workload(), expected);
        }
    }
}
