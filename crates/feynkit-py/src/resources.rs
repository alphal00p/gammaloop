//! Conversion of Python resource settings to the shared native controls.

use std::time::Duration;

use pyo3::exceptions::PyValueError;
use pyo3::{PyErr, PyResult};
use symbolica::utils::resources::{
    self as native_resources, ResourceBudget, ResourceLimitError, ResourceLimits,
};

pub(crate) fn limits(
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

pub(crate) fn checkpoint() -> PyResult<()> {
    native_resources::checkpoint().map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Scope native work and Rust-side result conversion for a Python request.
/// The caller preserves its domain-specific exception type for resource failures.
/// PyO3 input extraction and final Python marshalling remain boundary operations.
pub(crate) fn run<T>(
    limits: ResourceLimits,
    resource_error: impl Fn(ResourceLimitError) -> PyErr,
    operation: impl FnOnce() -> PyResult<T>,
) -> PyResult<T> {
    ResourceBudget::new(limits).run(|| {
        native_resources::checkpoint().map_err(&resource_error)?;
        let result = operation()?;
        native_resources::checkpoint().map_err(resource_error)?;
        Ok(result)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resource_cause_remains_visible_through_native_error_sources() {
        use std::error::Error;

        use feynkit_cff::CffError;
        use feynkit_graph::IntegralFamilyError;
        use feynkit_tensor::{TensorReductionError, WeingartenError};

        let cause = ResourceLimitError::Memory {
            rss_bytes: 2,
            additional_bytes: 0,
            limit_bytes: 1,
        };
        let errors: Vec<Box<dyn Error>> = vec![
            Box::new(IntegralFamilyError::Resource(cause.clone())),
            Box::new(TensorReductionError::Resource(cause.clone())),
            Box::new(WeingartenError::Resource(cause.clone())),
            Box::new(TensorReductionError::from(WeingartenError::Resource(
                cause.clone(),
            ))),
            Box::new(CffError::Resource(cause.clone())),
        ];
        for error in errors {
            assert_eq!(error.to_string(), cause.to_string());
            assert_eq!(
                error.source().and_then(|source| source.downcast_ref()),
                Some(&cause),
            );
        }
    }

    #[test]
    fn defaults_disable_resource_limits_and_zero_elapsed_is_explicit() {
        assert_eq!(limits(None, None).unwrap(), ResourceLimits::default());
        let configured = limits(Some(1024), Some(0.0)).unwrap();
        assert_eq!(configured.max_process_rss_bytes, Some(1024));
        assert_eq!(configured.max_elapsed, Some(Duration::ZERO));
    }

    #[test]
    fn invalid_elapsed_settings_are_errors() {
        for value in [f64::NAN, f64::INFINITY, -1.0] {
            assert!(limits(None, Some(value)).is_err());
        }
    }
}
