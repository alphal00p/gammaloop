use pyo3::{
    Bound, PyResult,
    exceptions::PyValueError,
    pyclass, pyfunction, pymethods,
    types::{PyModule, PyModuleMethods},
    wrap_pyfunction,
};
use spenso::structure::abstract_index::AbstractIndex;
use symbolica::{api::python::PythonExpression, atom::Symbol};

use crate::{
    color::{ColorCasimirSettings, ColorSimplifier, ColorSimplifySettings},
    dirac::{GammaChainOrdering, GammaSimplifier, GammaSimplifySettings},
    epsilon::EpsilonSimplifier,
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{
    gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pyfunction, gen_stub_pymethods,
};

/// Controls how open gamma chains are reordered during simplification.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    from_py_object,
    eq,
    eq_int,
    name = "GammaChainOrdering",
    module = "symbolica.community.idenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum PyGammaChainOrdering {
    /// Only move repeated gamma matrices toward their matching partner.
    RepeatedPairs,
    /// Canonically order gamma matrices with Clifford swaps.
    Canonical,
}

impl From<PyGammaChainOrdering> for GammaChainOrdering {
    fn from(value: PyGammaChainOrdering) -> Self {
        match value {
            PyGammaChainOrdering::RepeatedPairs => Self::RepeatedPairs,
            PyGammaChainOrdering::Canonical => Self::Canonical,
        }
    }
}

impl From<GammaChainOrdering> for PyGammaChainOrdering {
    fn from(value: GammaChainOrdering) -> Self {
        match value {
            GammaChainOrdering::RepeatedPairs => Self::RepeatedPairs,
            GammaChainOrdering::Canonical => Self::Canonical,
        }
    }
}

/// Configuration for gamma-chain simplification.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    name = "GammaSimplifySettings",
    module = "symbolica.community.idenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct PyGammaSimplifySettings {
    inner: GammaSimplifySettings,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyGammaSimplifySettings {
    #[new]
    #[pyo3(signature = (
        *,
        chain_ordering = PyGammaChainOrdering::RepeatedPairs,
        evaluate_traces = true,
        expand_three_gamma_epsilon = false
    ))]
    fn new(
        chain_ordering: PyGammaChainOrdering,
        evaluate_traces: bool,
        expand_three_gamma_epsilon: bool,
    ) -> Self {
        Self {
            inner: GammaSimplifySettings {
                chain_ordering: chain_ordering.into(),
                evaluate_traces,
                expand_three_gamma_epsilon,
            },
        }
    }

    /// FORM-like repeated-pair ordering with trace evaluation enabled.
    #[staticmethod]
    fn repeated_pairs() -> Self {
        Self {
            inner: GammaSimplifySettings::repeated_pairs(),
        }
    }

    /// Canonical gamma ordering with trace evaluation enabled.
    #[staticmethod]
    fn canonical() -> Self {
        Self {
            inner: GammaSimplifySettings::canonical(),
        }
    }

    /// Return a copy that leaves collected trace nodes unevaluated.
    fn without_trace_evaluation(&self) -> Self {
        Self {
            inner: self.inner.without_trace_evaluation(),
        }
    }

    /// Return a copy that expands three 4D gammas into gamma5 and epsilon.
    fn with_gamma5_epsilon_expansion(&self) -> Self {
        Self {
            inner: self.inner.with_gamma5_epsilon_expansion(),
        }
    }

    #[getter]
    fn chain_ordering(&self) -> PyGammaChainOrdering {
        self.inner.chain_ordering.into()
    }

    #[getter]
    fn evaluate_traces(&self) -> bool {
        self.inner.evaluate_traces
    }

    #[getter]
    fn expand_three_gamma_epsilon(&self) -> bool {
        self.inner.expand_three_gamma_epsilon
    }
}

impl From<PyGammaSimplifySettings> for GammaSimplifySettings {
    fn from(value: PyGammaSimplifySettings) -> Self {
        value.inner
    }
}

/// Configuration for SU(N) color simplification.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    name = "ColorSimplifySettings",
    module = "symbolica.community.idenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct PyColorSimplifySettings {
    inner: ColorSimplifySettings,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyColorSimplifySettings {
    #[new]
    #[pyo3(signature = (
        *,
        evaluate_traces = true,
        expand_cross_chain_fierz = true,
        substitute_cof_dimension_invariants = false
    ))]
    fn new(
        evaluate_traces: bool,
        expand_cross_chain_fierz: bool,
        substitute_cof_dimension_invariants: bool,
    ) -> Self {
        Self {
            inner: ColorSimplifySettings {
                evaluate_traces,
                expand_cross_chain_fierz,
                substitute_cof_dimension_invariants,
            },
        }
    }

    /// Return a copy that leaves collected color traces unevaluated.
    fn without_trace_evaluation(&self) -> Self {
        Self {
            inner: self.inner.without_trace_evaluation(),
        }
    }

    /// Return a copy that keeps separate open chains instead of applying Fierz.
    fn without_cross_chain_fierz_expansion(&self) -> Self {
        Self {
            inner: self.inner.without_cross_chain_fierz_expansion(),
        }
    }

    /// Return a copy that substitutes cof(N) dimension invariants.
    fn with_cof_dimension_invariants(&self) -> Self {
        Self {
            inner: self.inner.with_cof_dimension_invariants(),
        }
    }

    #[getter]
    fn evaluate_traces(&self) -> bool {
        self.inner.evaluate_traces
    }

    #[getter]
    fn expand_cross_chain_fierz(&self) -> bool {
        self.inner.expand_cross_chain_fierz
    }

    #[getter]
    fn substitute_cof_dimension_invariants(&self) -> bool {
        self.inner.substitute_cof_dimension_invariants
    }
}

impl From<PyColorSimplifySettings> for ColorSimplifySettings {
    fn from(value: PyColorSimplifySettings) -> Self {
        value.inner
    }
}

/// Configuration for rewriting color invariants into the Casimir basis.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    name = "ColorCasimirSettings",
    module = "symbolica.community.idenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct PyColorCasimirSettings {
    inner: ColorCasimirSettings,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyColorCasimirSettings {
    #[new]
    #[pyo3(signature = (*, rewrite_nc = true, rewrite_na = true, substitute_tr = false))]
    fn new(rewrite_nc: bool, rewrite_na: bool, substitute_tr: bool) -> Self {
        Self {
            inner: ColorCasimirSettings {
                rewrite_nc,
                rewrite_na,
                substitute_tr,
            },
        }
    }

    /// Return a copy that substitutes the conventional TR = 1/2 normalization.
    fn with_trace_normalization(&self) -> Self {
        Self {
            inner: self.inner.with_trace_normalization(),
        }
    }

    #[getter]
    fn rewrite_nc(&self) -> bool {
        self.inner.rewrite_nc
    }

    #[getter]
    fn rewrite_na(&self) -> bool {
        self.inner.rewrite_na
    }

    #[getter]
    fn substitute_tr(&self) -> bool {
        self.inner.substitute_tr
    }
}

impl From<PyColorCasimirSettings> for ColorCasimirSettings {
    fn from(value: PyColorCasimirSettings) -> Self {
        value.inner
    }
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Collect products of gamma matrices into typed open chains or closed traces.
fn collect_gamma_chains(expression: &PythonExpression) -> PythonExpression {
    expression.expr.collect_gamma_chains().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Simplify products involving the distinguished gamma0 matrix.
fn simplify_gamma0(expression: &PythonExpression) -> PythonExpression {
    expression.expr.simplify_gamma0().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Rewrite conjugated gamma matrices using gamma0 insertions.
fn simplify_gamma_conjugate(expression: &PythonExpression) -> PyResult<PythonExpression> {
    expression
        .expr
        .simplify_gamma_conj::<AbstractIndex>()
        .map(Into::into)
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Simplify epsilon-metric contractions and products of epsilon tensors.
fn simplify_epsilon(expression: &PythonExpression) -> PythonExpression {
    expression.expr.simplify_epsilon().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (expression, settings = None))]
/// Rewrite color scalar factors into the CA/CF Casimir basis.
fn to_color_casimir(
    expression: &PythonExpression,
    settings: Option<PyColorCasimirSettings>,
) -> PythonExpression {
    expression
        .expr
        .to_color_casimir_with(settings.map(Into::into).unwrap_or_default())
        .into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Rewrite supported cof(N) invariants as explicit functions of N.
fn to_cof_dimension_invariants(expression: &PythonExpression) -> PythonExpression {
    expression.expr.to_cof_dimension_invariants().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Collect terms by their indexed color structure.
fn collect_color(expression: &PythonExpression) -> PythonExpression {
    expression.expr.collect_color().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Collect terms by their color-scalar invariant factors.
fn collect_color_constants(expression: &PythonExpression) -> PythonExpression {
    expression.expr.collect_color_constants().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Move each color-dependent coefficient into a wrapper function.
fn wrap_color(expression: &PythonExpression, symbol: Symbol) -> PythonExpression {
    expression.expr.wrap_color(symbol).into()
}

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyGammaChainOrdering>()?;
    m.add_class::<PyGammaSimplifySettings>()?;
    m.add_class::<PyColorSimplifySettings>()?;
    m.add_class::<PyColorCasimirSettings>()?;

    m.add_function(wrap_pyfunction!(collect_gamma_chains, m)?)?;
    m.add_function(wrap_pyfunction!(simplify_gamma0, m)?)?;
    m.add_function(wrap_pyfunction!(simplify_gamma_conjugate, m)?)?;
    m.add_function(wrap_pyfunction!(simplify_epsilon, m)?)?;
    m.add_function(wrap_pyfunction!(to_color_casimir, m)?)?;
    m.add_function(wrap_pyfunction!(to_cof_dimension_invariants, m)?)?;
    m.add_function(wrap_pyfunction!(collect_color, m)?)?;
    m.add_function(wrap_pyfunction!(collect_color_constants, m)?)?;
    m.add_function(wrap_pyfunction!(wrap_color, m)?)?;

    Ok(())
}
