use std::sync::{
    OnceLock,
    atomic::{AtomicBool, Ordering},
};

use symbolica::LicenseManager;

/// Policy used to configure Rayon for operations that manipulate Symbolica atoms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SymbolicParallelism {
    /// Enable Rayon only when Symbolica is licensed at configuration time.
    Auto,
    /// Keep symbolic operations on the calling thread.
    Serial,
    /// Allow symbolic operations to use Rayon workers.
    Parallel,
}

static SYMBOLICA_RAYON_ENABLED: OnceLock<AtomicBool> = OnceLock::new();

fn cached_setting() -> &'static AtomicBool {
    SYMBOLICA_RAYON_ENABLED.get_or_init(|| {
        AtomicBool::new(SymbolicParallelism::Auto.resolve_with(LicenseManager::is_licensed))
    })
}

impl SymbolicParallelism {
    fn resolve_with(self, is_licensed: impl FnOnce() -> bool) -> bool {
        match self {
            Self::Auto => is_licensed(),
            Self::Serial => false,
            Self::Parallel => true,
        }
    }
}

/// Configure whether operations involving Symbolica atoms may use Rayon.
///
/// [`SymbolicParallelism::Auto`] queries the Symbolica license exactly once,
/// when this function is called. Later operations only read the cached boolean.
/// Configure this before starting tensor operations; changing it concurrently
/// with an active operation is not supported.
pub fn set_symbolica_rayon_enabled(policy: SymbolicParallelism) {
    set_symbolica_rayon_enabled_with(policy, LicenseManager::is_licensed);
}

fn set_symbolica_rayon_enabled_with(
    policy: SymbolicParallelism,
    is_licensed: impl FnOnce() -> bool,
) {
    let enabled = policy.resolve_with(is_licensed);
    SYMBOLICA_RAYON_ENABLED
        .get_or_init(|| AtomicBool::new(enabled))
        .store(enabled, Ordering::Relaxed);
}

#[cfg(test)]
pub(crate) struct SymbolicParallelismTestGuard {
    previous: bool,
    _lock: std::sync::MutexGuard<'static, ()>,
}

#[cfg(test)]
static SYMBOLIC_PARALLELISM_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

#[cfg(test)]
impl Drop for SymbolicParallelismTestGuard {
    fn drop(&mut self) {
        cached_setting().store(self.previous, Ordering::Relaxed);
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
    let previous = symbolica_rayon_enabled();
    set_symbolica_rayon_enabled_with(policy, is_licensed);
    SymbolicParallelismTestGuard {
        previous,
        _lock: lock,
    }
}

/// Return the resolved symbolic Rayon setting.
///
/// If no policy has been configured yet, this performs the default
/// [`SymbolicParallelism::Auto`] initialization and queries the license once.
/// Subsequent reads only load the cached boolean.
pub fn symbolica_rayon_enabled() -> bool {
    cached_setting().load(Ordering::Relaxed)
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;

    use super::*;

    #[test]
    fn auto_checks_the_license_once_when_resolved() {
        let calls = Cell::new(0);
        let enabled = SymbolicParallelism::Auto.resolve_with(|| {
            calls.set(calls.get() + 1);
            false
        });

        assert!(!enabled);
        assert_eq!(calls.get(), 1);
    }

    #[test]
    fn explicit_policies_do_not_check_the_license() {
        for (policy, expected) in [
            (SymbolicParallelism::Serial, false),
            (SymbolicParallelism::Parallel, true),
        ] {
            let calls = Cell::new(0);
            let enabled = policy.resolve_with(|| {
                calls.set(calls.get() + 1);
                !expected
            });

            assert_eq!(enabled, expected);
            assert_eq!(calls.get(), 0);
        }
    }

    #[test]
    fn setter_caches_the_resolved_boolean() {
        let _guard = scoped_symbolica_rayon_setting_for_test(SymbolicParallelism::Serial, || false);
        assert!(!symbolica_rayon_enabled());

        set_symbolica_rayon_enabled(SymbolicParallelism::Parallel);
        assert!(symbolica_rayon_enabled());
    }
}
