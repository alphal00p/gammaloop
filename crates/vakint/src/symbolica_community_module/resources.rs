//! Python boundary for the shared native RSS and elapsed-time budget.

use std::time::Duration;

use pyo3::PyResult;
use pyo3::exceptions::PyValueError;
use symbolica::utils::resources as native_resources;
use symbolica::utils::resources::ResourceBudget;

use super::vakint_to_python_error;
use crate::{ResourceLimits, VakintError};

pub(super) fn limits(
    max_process_rss_bytes: Option<u64>,
    max_elapsed_seconds: Option<f64>,
) -> PyResult<ResourceLimits> {
    let max_elapsed = max_elapsed_seconds
        .map(Duration::try_from_secs_f64)
        .transpose()
        .map_err(|_| PyValueError::new_err("max_elapsed_seconds must be finite and nonnegative"))?;
    Ok(ResourceLimits {
        max_process_rss_bytes,
        max_elapsed,
    })
}

fn checkpoint() -> PyResult<()> {
    native_resources::checkpoint()
        .map_err(VakintError::from)
        .map_err(vakint_to_python_error)
}

pub(super) fn run<T>(
    limits: ResourceLimits,
    operation: impl FnOnce() -> PyResult<T>,
) -> PyResult<T> {
    ResourceBudget::new(limits).run(|| {
        checkpoint()?;
        let result = operation();
        let checked = checkpoint();
        match result {
            Ok(value) => {
                checked?;
                Ok(value)
            }
            // Python callback errors and signals retain their original type,
            // identity and durable cache_report even if the request also timed out.
            Err(error) => Err(error),
        }
    })
}

/// Native commit-aware operations check immediately before publication. A
/// later facade check could report failure after an attachment or file changed.
pub(super) fn run_committing<T>(
    limits: ResourceLimits,
    operation: impl FnOnce() -> PyResult<T>,
) -> PyResult<T> {
    ResourceBudget::new(limits).run(|| {
        checkpoint()?;
        operation()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::Python;
    use pyo3::types::PyAnyMethods;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    #[test]
    fn resource_settings_are_unbounded_by_default_and_validate_time() {
        assert_eq!(limits(None, None).unwrap(), ResourceLimits::default());
        assert_eq!(
            limits(Some(1), Some(0.0)).unwrap().max_elapsed,
            Some(Duration::ZERO)
        );
        for seconds in [f64::NAN, f64::INFINITY, -1.0] {
            assert!(limits(None, Some(seconds)).is_err());
        }
    }

    #[test]
    fn callback_exception_keeps_identity_when_the_postcheck_also_fails() {
        Python::initialize();
        Python::attach(|py| {
            let original = PyValueError::new_err("original callback failure");
            let expected = original.clone_ref(py);
            let checks = Arc::new(AtomicUsize::new(0));
            let monitor = Arc::clone(&checks);
            let outer = ResourceBudget::with_memory_monitor(
                ResourceLimits {
                    max_process_rss_bytes: Some(10),
                    ..ResourceLimits::default()
                },
                move || {
                    Ok(if monitor.fetch_add(1, Ordering::SeqCst) == 0 {
                        0
                    } else {
                        11
                    })
                },
            );
            let result: PyResult<()> =
                outer.run(|| run(ResourceLimits::default(), || Err(original)));
            let error = result.unwrap_err();
            assert!(error.value(py).is(expected.value(py)));
            assert!(checks.load(Ordering::SeqCst) >= 2);
        });
    }
}
