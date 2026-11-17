//! Preserve native failure chains at the shared Python presentation boundary.
use std::error::Error;

use pyo3::exceptions::PyRuntimeError;
use pyo3::{PyErr, Python};

/// Attach a Rust failure chain without changing the outer Python exception type.
///
/// ```
/// use std::io::{Error, ErrorKind};
/// use linnet_render_py::PyErrExt;
/// use pyo3::exceptions::PyValueError;
/// use pyo3::Python;
/// Python::initialize();
/// Python::attach(|py| {
///     let source = Error::new(ErrorKind::InvalidData, "invalid drawing");
///     let error = PyValueError::new_err("cannot render").with_rust_cause(&source);
///     assert!(error.is_instance_of::<PyValueError>(py));
///     assert!(error.cause(py).is_some());
/// });
/// ```
pub trait PyErrExt {
    /// Retain every source as a Python cause. Existing Python exceptions retain
    /// their type, traceback and cause chain. Like other Python operations,
    /// this requires an initialized interpreter (the host initializes it).
    fn with_rust_cause(self, source: &(dyn Error + 'static)) -> Self;
}

impl PyErrExt for PyErr {
    fn with_rust_cause(self, source: &(dyn Error + 'static)) -> Self {
        Python::attach(|py| {
            let cause = if let Some(error) = source.downcast_ref::<PyErr>() {
                error.clone_ref(py)
            } else {
                let error = PyRuntimeError::new_err(source.to_string());
                if let Some(source) = source.source() {
                    error.with_rust_cause(source)
                } else {
                    error
                }
            };
            self.set_cause(py, Some(cause));
            self
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::exceptions::{PyTypeError, PyValueError};
    use thiserror::Error;

    #[derive(Debug, Error)]
    #[error("invalid record")]
    struct RecordError(#[source] PyErr);

    #[test]
    fn native_sources_preserve_python_exception_types_and_nested_causes() {
        Python::initialize();
        Python::attach(|py| {
            let source = RecordError(PyTypeError::new_err("bad coordinate"));
            let error = PyValueError::new_err("cannot render").with_rust_cause(&source);
            assert!(error.is_instance_of::<PyValueError>(py));
            let record = error.cause(py).unwrap();
            assert!(record.to_string().contains("invalid record"));
            let coordinate = record.cause(py).unwrap();
            assert!(coordinate.is_instance_of::<PyTypeError>(py));
            assert!(coordinate.to_string().contains("bad coordinate"));
        });
    }
}
