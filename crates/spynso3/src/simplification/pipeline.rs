//! Composition of existing algebra passes; each pass retains its native dimension rules.
use super::{PyColorSimplifySettings, PyGammaSimplifySettings};
use pyo3::{exceptions::PyValueError, prelude::*};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// Immutable choice of algebra passes for TensorExpression.simplify().
/// Dimensions come from tensor slots. No dimensional substitution or gamma5
/// prescription is chosen by this settings object. Individual algebra identities
/// can introduce sums; expand controls additional full polynomial expansion.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    get_all,
    name = "SimplifySettings",
    module = "symbolica.community.spenso"
)]
#[derive(Clone, Copy)]
pub(crate) struct PySimplifySettings {
    pub(crate) metrics: bool,
    pub(crate) gamma: Option<PyGammaSimplifySettings>,
    pub(crate) color: Option<PyColorSimplifySettings>,
    pub(crate) epsilon: bool,
    pub(crate) expand: bool,
    pub(crate) max_passes: usize,
}

impl Default for PySimplifySettings {
    fn default() -> Self {
        Self {
            metrics: true,
            gamma: None,
            color: None,
            epsilon: false,
            expand: false,
            max_passes: 16,
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySimplifySettings {
    /// Select passes explicitly. Defaults contract metrics without expanding scalar algebra.
    #[new]
    #[pyo3(signature=(*, metrics=true, gamma=None, color=None, epsilon=false, expand=false, max_passes=16))]
    fn new(
        metrics: bool,
        gamma: Option<PyGammaSimplifySettings>,
        color: Option<PyColorSimplifySettings>,
        epsilon: bool,
        expand: bool,
        max_passes: usize,
    ) -> PyResult<Self> {
        if max_passes == 0 {
            return Err(PyValueError::new_err("max_passes must be positive"));
        }
        Ok(Self {
            metrics,
            gamma,
            color,
            epsilon,
            expand,
            max_passes,
        })
    }

    /// Enable metric, gamma, color, and epsilon passes with their native defaults.
    /// Full polynomial expansion and the optional three-gamma epsilon identity stay off.
    #[staticmethod]
    fn hep() -> Self {
        Self {
            gamma: Some(PyGammaSimplifySettings::repeated_pairs()),
            color: Some(PyColorSimplifySettings::new(true, true, false)),
            epsilon: true,
            ..Self::default()
        }
    }

    fn __repr__(self_: PyRef<'_, Self>) -> PyResult<String> {
        let py = self_.py();
        let object = self_.into_pyobject(py)?;
        crate::display::constructor_repr(
            object.as_any(),
            &[
                ("metrics", "metrics"),
                ("gamma", "gamma"),
                ("color", "color"),
                ("epsilon", "epsilon"),
                ("expand", "expand"),
                ("max_passes", "max_passes"),
            ],
        )
    }
}
