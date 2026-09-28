//! Python conversion for the shared tensor-algebra scheduler.
use super::{PyColorSimplifySettings, PyGammaSimplifySettings};
use idenso::tensor::simplification::SimplifySettings;
use pyo3::prelude::*;
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// Immutable selection of tensor identities. Dimensions come from tensor slots.
/// Materialize scalar polynomials explicitly with the result's `expand()` method.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "SimplifySettings",
    module = "symbolica.community.spenso"
)]
#[derive(Clone, Copy, Default)]
pub(crate) struct PySimplifySettings {
    inner: SimplifySettings,
}

impl PySimplifySettings {
    pub(crate) fn rust(&self) -> &SimplifySettings {
        &self.inner
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySimplifySettings {
    /// Select passes explicitly. Defaults contract metrics without expanding scalar algebra.
    #[new]
    #[pyo3(signature=(*, metrics=true, gamma=None, color=None, epsilon=false, max_passes=16))]
    pub(crate) fn new(
        metrics: bool,
        gamma: Option<PyGammaSimplifySettings>,
        color: Option<PyColorSimplifySettings>,
        epsilon: bool,
        max_passes: usize,
    ) -> PyResult<Self> {
        let inner = SimplifySettings {
            metrics,
            gamma: gamma.map(|settings| settings.rust()),
            color: color.map(|settings| settings.rust()),
            epsilon,
            max_passes,
        };
        inner
            .validate()
            .map_err(crate::expression::TensorExpression::inference_error)?;
        Ok(Self { inner })
    }

    /// Enable metric, gamma, color, and epsilon identities with their native defaults.
    #[staticmethod]
    fn hep() -> Self {
        Self {
            inner: SimplifySettings::hep(),
        }
    }

    #[getter]
    fn metrics(&self) -> bool {
        self.inner.metrics
    }
    #[getter]
    fn gamma(&self) -> Option<PyGammaSimplifySettings> {
        self.inner.gamma.map(PyGammaSimplifySettings::from_rust)
    }
    #[getter]
    fn color(&self) -> Option<PyColorSimplifySettings> {
        self.inner.color.map(PyColorSimplifySettings::from_rust)
    }
    #[getter]
    fn epsilon(&self) -> bool {
        self.inner.epsilon
    }
    #[getter]
    fn max_passes(&self) -> usize {
        self.inner.max_passes
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
                ("max_passes", "max_passes"),
            ],
        )
    }
}
