#[cfg(not(feature = "python_stubgen"))]
use pyo3::create_exception;
use pyo3::{
    Borrowed, Bound, FromPyObject, PyAny, PyErr, PyResult,
    exceptions::{PyTypeError, PyValueError},
    types::{PyModule, PyModuleMethods},
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::create_exception;

create_exception!(
    symbolica.community.tensor,
    CanonicalizationError,
    PyValueError,
    python_doc!("CanonicalizationError")
);
create_exception!(
    symbolica.community.tensor,
    CookingError,
    PyTypeError,
    python_doc!("CookingError")
);
create_exception!(
    symbolica.community.tensor,
    DiracAdjointError,
    PyValueError,
    python_doc!("DiracAdjointError")
);
create_exception!(
    symbolica.community.tensor,
    NetworkToolingError,
    PyValueError,
    python_doc!("NetworkToolingError")
);

/// Index payload encoding accepted by the Python tensor entry points.
#[derive(Clone, Copy)]
pub(crate) enum Intern {
    Indices,
    Flattened,
}

impl Intern {
    pub(crate) fn rust(self) -> idenso::CookSettings {
        let settings = idenso::CookSettings::indices();
        match self {
            Self::Indices => settings.with_mode(idenso::CookMode::ReversibleEncoding),
            Self::Flattened => settings,
        }
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for Intern {
    type Error = PyErr;

    fn extract(value: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        match value.extract::<String>()?.as_str() {
            "indices" => Ok(Self::Indices),
            "flattened" => Ok(Self::Flattened),
            _ => Err(PyValueError::new_err(
                "intern must be None, 'indices', or 'flattened'",
            )),
        }
    }
}

#[cfg(feature = "python_stubgen")]
impl pyo3_stub_gen::PyStubType for Intern {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        pyo3_stub_gen::TypeInfo {
            name: "typing.Literal['indices', 'flattened']".into(),
            import: ["typing".into()].into(),
        }
    }
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "CanonicalizationError",
        module.py().get_type::<CanonicalizationError>(),
    )?;
    module.add("CookingError", module.py().get_type::<CookingError>())?;
    module.add(
        "DiracAdjointError",
        module.py().get_type::<DiracAdjointError>(),
    )?;
    module.add(
        "NetworkToolingError",
        module.py().get_type::<NetworkToolingError>(),
    )?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::Python;

    #[test]
    fn network_tooling_exception_is_a_value_error() {
        Python::initialize();
        Python::attach(|py| {
            assert!(
                CanonicalizationError::new_err("canonicalization failure")
                    .is_instance_of::<PyValueError>(py)
            );
            assert!(
                DiracAdjointError::new_err("adjoint failure").is_instance_of::<PyValueError>(py)
            );
            assert!(
                NetworkToolingError::new_err("network failure").is_instance_of::<PyValueError>(py)
            );
        });
    }
}
