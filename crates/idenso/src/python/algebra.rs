#[cfg(not(feature = "python_stubgen"))]
use pyo3::create_exception;
use pyo3::{
    Bound, PyResult,
    exceptions::PyValueError,
    pyclass, pyfunction, pymethods,
    types::{PyModule, PyModuleMethods},
    wrap_pyfunction,
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::create_exception;
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{
    gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pyfunction, gen_stub_pymethods,
};
use spenso::structure::abstract_index::AbstractIndex;
use symbolica::{api::python::PythonExpression, atom::Symbol};

use crate::{
    color::{ColorCasimirSettings, ColorSimplifier, ColorSimplifySettings},
    dirac::{GammaChainOrdering, GammaSimplifier, GammaSimplifySettings},
    epsilon::EpsilonSimplifier,
    python::tooling::RegisteredRepresentation,
};

create_exception!(
    symbolica.community.idenso,
    GammaConjugationError,
    PyValueError,
    "Raised when conjugated gamma matrices cannot be rewritten consistently."
);

/// Controls how open gamma chains are reordered during simplification.
///
/// Available values are `RepeatedPairs`, which only moves matching matrices together, and
/// `Canonical`, which canonically orders the complete open chain.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    frozen,
    from_py_object,
    eq,
    eq_int,
    name = "GammaChainOrdering",
    module = "symbolica.community.idenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PyGammaChainOrdering {
    /// Move repeated gamma matrices toward each other without reordering unrelated factors.
    RepeatedPairs,
    /// Canonically order open chains using adjacent Clifford-algebra swaps.
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

/// Immutable configuration for gamma-chain simplification.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "GammaSimplifySettings",
    module = "symbolica.community.idenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PyGammaSimplifySettings {
    inner: GammaSimplifySettings,
}

impl PyGammaSimplifySettings {
    /// Construct settings from explicit Rust-side values.
    pub fn new(
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

    pub fn rust(&self) -> GammaSimplifySettings {
        self.inner
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyGammaSimplifySettings {
    /// Configure gamma-chain ordering, trace evaluation, and the optional 4D three-gamma identity.
    ///
    /// `chain_ordering=None` selects `GammaChainOrdering.RepeatedPairs`.
    #[new]
    #[pyo3(
        signature = (
            *,
            chain_ordering = None,
            evaluate_traces = true,
            expand_three_gamma_epsilon = false
        ),
        text_signature = "(*, chain_ordering=None, evaluate_traces=True, expand_three_gamma_epsilon=False)"
    )]
    pub fn py_new(
        chain_ordering: Option<PyGammaChainOrdering>,
        evaluate_traces: bool,
        expand_three_gamma_epsilon: bool,
    ) -> Self {
        Self::new(
            chain_ordering.unwrap_or(PyGammaChainOrdering::RepeatedPairs),
            evaluate_traces,
            expand_three_gamma_epsilon,
        )
    }

    /// Use FORM-like repeated-pair ordering and evaluate closed traces.
    #[staticmethod]
    pub fn repeated_pairs() -> Self {
        Self {
            inner: GammaSimplifySettings::repeated_pairs(),
        }
    }

    /// Canonically order open gamma chains and evaluate closed traces.
    #[staticmethod]
    pub fn canonical() -> Self {
        Self {
            inner: GammaSimplifySettings::canonical(),
        }
    }

    /// Ordering strategy used for open gamma chains.
    #[getter]
    pub fn chain_ordering(&self) -> PyGammaChainOrdering {
        self.inner.chain_ordering.into()
    }

    /// Whether closed gamma chains are evaluated as traces.
    #[getter]
    pub fn evaluate_traces(&self) -> bool {
        self.inner.evaluate_traces
    }

    /// Whether three four-dimensional gammas expand into a gamma5-epsilon basis.
    #[getter]
    pub fn expand_three_gamma_epsilon(&self) -> bool {
        self.inner.expand_three_gamma_epsilon
    }
}

/// Immutable configuration for color simplification.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "ColorSimplifySettings",
    module = "symbolica.community.idenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PyColorSimplifySettings {
    inner: ColorSimplifySettings,
}

impl PyColorSimplifySettings {
    pub fn rust(&self) -> ColorSimplifySettings {
        self.inner
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyColorSimplifySettings {
    /// Configure color-trace evaluation, cross-chain Fierz expansion, and invariant substitution.
    #[new]
    #[pyo3(signature = (
        *,
        evaluate_traces = true,
        expand_cross_chain_fierz = true,
        substitute_cof_dimension_invariants = false
    ))]
    pub fn new(
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

    /// Whether closed color chains are evaluated as traces.
    #[getter]
    pub fn evaluate_traces(&self) -> bool {
        self.inner.evaluate_traces
    }

    /// Whether generators on different open chains are expanded with the Fierz identity.
    #[getter]
    pub fn expand_cross_chain_fierz(&self) -> bool {
        self.inner.expand_cross_chain_fierz
    }

    /// Whether supported `cof(N)` invariants are replaced by explicit dimension formulas.
    #[getter]
    pub fn substitute_cof_dimension_invariants(&self) -> bool {
        self.inner.substitute_cof_dimension_invariants
    }
}

/// Immutable configuration for rewriting color invariants into a Casimir basis.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "ColorCasimirSettings",
    module = "symbolica.community.idenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PyColorCasimirSettings {
    inner: ColorCasimirSettings,
}

impl PyColorCasimirSettings {
    pub fn rust(&self) -> ColorCasimirSettings {
        self.inner
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyColorCasimirSettings {
    /// Configure the SU(N) dimension and fundamental-index normalizations used by Casimir rewriting.
    #[new]
    #[pyo3(signature = (
        *,
        rewrite_fundamental_dimension = true,
        substitute_fundamental_index = false
    ))]
    pub fn new(rewrite_fundamental_dimension: bool, substitute_fundamental_index: bool) -> Self {
        Self {
            inner: ColorCasimirSettings {
                rewrite_fundamental_dimension,
                substitute_fundamental_index,
            },
        }
    }

    /// Whether the fundamental dimension is rewritten with the SU(N) relation `d_F = C_A`.
    #[getter]
    pub fn rewrite_fundamental_dimension(&self) -> bool {
        self.inner.rewrite_fundamental_dimension
    }

    /// Whether the fundamental Dynkin index is replaced by `T_F = 1/2`.
    #[getter]
    pub fn substitute_fundamental_index(&self) -> bool {
        self.inner.substitute_fundamental_index
    }
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Convert bispinor tensors into chain/trace shorthands and join adjacent gamma chains.
pub fn collect_gamma_chains(expression: &PythonExpression) -> PythonExpression {
    expression.expr.collect_gamma_chains().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Simplify products and linear combinations involving the time-like gamma matrix `gamma0`.
pub fn simplify_gamma0(expression: &PythonExpression) -> PythonExpression {
    expression.expr.simplify_gamma0().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Rewrite conjugated gamma matrices as gamma0-sandwiched matrices with fresh spinor indices.
///
/// Raises `GammaConjugationError` when the expression cannot be rewritten consistently.
pub fn simplify_gamma_conjugate(expression: &PythonExpression) -> PyResult<PythonExpression> {
    expression
        .expr
        .simplify_gamma_conj::<AbstractIndex>()
        .map(Into::into)
        .map_err(|error| GammaConjugationError::new_err(error.to_string()))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Simplify Levi-Civita/metric contractions and pairs of Levi-Civita tensors to a fixed point.
pub fn simplify_epsilon(expression: &PythonExpression) -> PythonExpression {
    expression.expr.simplify_epsilon().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (expression, *, fundamental, adjoint, settings = None))]
/// Rewrite the supplied color-representation dimensions and invariants into a Casimir basis.
///
/// `fundamental` and `adjoint` accept Spynso `Representation` objects or their symbolic
/// expressions. Only scalar coefficient positions are rewritten.
pub fn to_color_casimir(
    expression: &PythonExpression,
    fundamental: RegisteredRepresentation,
    adjoint: RegisteredRepresentation,
    settings: Option<&PyColorCasimirSettings>,
) -> PythonExpression {
    let fundamental = fundamental.to_atom();
    let adjoint = adjoint.to_atom();
    expression
        .expr
        .to_color_casimir_with(
            fundamental.as_view(),
            adjoint.as_view(),
            settings.map_or_else(ColorCasimirSettings::default, |s| s.rust()),
        )
        .into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Replace supported `cof(N)` Casimir, Dynkin-index, and Gram invariants by dimension formulas.
pub fn to_cof_dimension_invariants(expression: &PythonExpression) -> PythonExpression {
    expression.expr.to_cof_dimension_invariants().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Factor an expression around tensors carrying fundamental, antifundamental, or adjoint color.
pub fn collect_color(expression: &PythonExpression) -> PythonExpression {
    expression.expr.collect_color().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Factor an expression around recognized scalar color invariants such as Casimirs and indices.
pub fn collect_color_constants(expression: &PythonExpression) -> PythonExpression {
    expression.expr.collect_color_constants().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Expand around color structures and wrap each resulting scalar coefficient with `symbol`.
pub fn wrap_color(expression: &PythonExpression, symbol: Symbol) -> PythonExpression {
    expression.expr.wrap_color(symbol).into()
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "GammaConjugationError",
        module.py().get_type::<GammaConjugationError>(),
    )?;
    module.add_class::<PyGammaChainOrdering>()?;
    module.add_class::<PyGammaSimplifySettings>()?;
    module.add_class::<PyColorSimplifySettings>()?;
    module.add_class::<PyColorCasimirSettings>()?;

    module.add_function(wrap_pyfunction!(collect_gamma_chains, module)?)?;
    module.add_function(wrap_pyfunction!(simplify_gamma0, module)?)?;
    module.add_function(wrap_pyfunction!(simplify_gamma_conjugate, module)?)?;
    module.add_function(wrap_pyfunction!(simplify_epsilon, module)?)?;
    module.add_function(wrap_pyfunction!(to_color_casimir, module)?)?;
    module.add_function(wrap_pyfunction!(to_cof_dimension_invariants, module)?)?;
    module.add_function(wrap_pyfunction!(collect_color, module)?)?;
    module.add_function(wrap_pyfunction!(collect_color_constants, module)?)?;
    module.add_function(wrap_pyfunction!(wrap_color, module)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_settings_match_rust_defaults() {
        assert_eq!(
            PyGammaSimplifySettings::new(PyGammaChainOrdering::RepeatedPairs, true, false).rust(),
            GammaSimplifySettings::default()
        );
        assert_eq!(
            PyColorSimplifySettings::new(true, true, false).rust(),
            ColorSimplifySettings::default()
        );
        assert_eq!(
            PyColorCasimirSettings::new(true, false).rust(),
            ColorCasimirSettings::default()
        );
    }
}
