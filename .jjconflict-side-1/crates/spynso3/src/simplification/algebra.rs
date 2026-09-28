use idenso::{
    color::{ColorCasimirSettings, ColorSimplifySettings},
    dirac::{GammaChainOrdering, GammaOutput, GammaSimplifySettings},
};
use pyo3::{
    Bound, PyResult,
    exceptions::PyValueError,
    pyclass, pymethods,
    types::{PyModule, PyModuleMethods},
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

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
    module = "symbolica.community.spenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyGammaChainOrdering {
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
    module = "symbolica.community.spenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PyGammaSimplifySettings {
    inner: GammaSimplifySettings,
}

impl PyGammaSimplifySettings {
    pub(crate) fn from_rust(inner: GammaSimplifySettings) -> Self {
        Self { inner }
    }

    /// Construct settings from explicit Rust-side values.
    pub(crate) fn new(
        chain_ordering: PyGammaChainOrdering,
        evaluate_traces: bool,
        gamma0: bool,
        conjugate: bool,
        expand_three_gamma_epsilon: bool,
        output: GammaOutput,
    ) -> Self {
        Self {
            inner: GammaSimplifySettings {
                chain_ordering: chain_ordering.into(),
                evaluate_traces,
                gamma0,
                conjugate,
                expand_three_gamma_epsilon,
                output,
            },
        }
    }

    pub(crate) fn rust(&self) -> GammaSimplifySettings {
        self.inner
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyGammaSimplifySettings {
    fn __repr__(self_: pyo3::PyRef<'_, Self>) -> PyResult<String> {
        let py = self_.py();
        let object = pyo3::IntoPyObject::into_pyobject(self_, py)?;
        crate::display::constructor_repr(
            object.as_any(),
            &[
                ("output", "output"),
                ("chain_ordering", "chain_ordering"),
                ("evaluate_traces", "evaluate_traces"),
                ("gamma0", "gamma0"),
                ("conjugate", "conjugate"),
                ("expand_three_gamma_epsilon", "expand_three_gamma_epsilon"),
            ],
        )
    }

    /// Configure gamma identities, trace evaluation, and the optional 4D epsilon identity.
    /// Conjugation runs before gamma0 reduction and ordinary chain identities.
    /// Expand the returned alias carrier explicitly when a polynomial is needed.
    ///
    /// `output="chains"` only joins gamma chains, without Clifford reduction or trace evaluation.
    /// `output="reduced"` applies the configured gamma identities.
    /// `chain_ordering=None` selects `GammaChainOrdering.RepeatedPairs`.
    #[new]
    #[pyo3(
        signature = (
            *,
            output = "reduced",
            chain_ordering = None,
            evaluate_traces = true,
            gamma0 = false,
            conjugate = false,
            expand_three_gamma_epsilon = false
        ),
        text_signature = "(*, output='reduced', chain_ordering=None, evaluate_traces=True, gamma0=False, conjugate=False, expand_three_gamma_epsilon=False)"
    )]
    pub(crate) fn py_new(
        #[gen_stub(override_type(type_repr="typing.Literal['reduced', 'chains']", imports=("typing")))]
        output: &str,
        chain_ordering: Option<PyGammaChainOrdering>,
        evaluate_traces: bool,
        gamma0: bool,
        conjugate: bool,
        expand_three_gamma_epsilon: bool,
    ) -> PyResult<Self> {
        let output = match output {
            "reduced" => GammaOutput::Reduced,
            "chains" => GammaOutput::Chains,
            _ => {
                return Err(PyValueError::new_err(
                    "output must be 'reduced' or 'chains'",
                ));
            }
        };
        Ok(Self::new(
            chain_ordering.unwrap_or(PyGammaChainOrdering::RepeatedPairs),
            evaluate_traces,
            gamma0,
            conjugate,
            expand_three_gamma_epsilon,
            output,
        ))
    }

    /// Requested reduced algebra or chain-only output.
    #[getter]
    #[gen_stub(override_return_type(type_repr="typing.Literal['reduced', 'chains']", imports=("typing")))]
    pub(crate) fn output(&self) -> &'static str {
        match self.inner.output {
            GammaOutput::Reduced => "reduced",
            GammaOutput::Chains => "chains",
        }
    }

    /// Use FORM-like repeated-pair ordering and evaluate closed traces.
    #[staticmethod]
    pub(crate) fn repeated_pairs() -> Self {
        Self {
            inner: GammaSimplifySettings::repeated_pairs(),
        }
    }

    /// Canonically order open gamma chains and evaluate closed traces.
    #[staticmethod]
    pub(crate) fn canonical() -> Self {
        Self {
            inner: GammaSimplifySettings::canonical(),
        }
    }

    /// Ordering strategy used for open gamma chains.
    #[getter]
    pub(crate) fn chain_ordering(&self) -> PyGammaChainOrdering {
        self.inner.chain_ordering.into()
    }

    /// Whether closed gamma chains are evaluated as traces.
    #[getter]
    pub(crate) fn evaluate_traces(&self) -> bool {
        self.inner.evaluate_traces
    }

    /// Whether to reduce gamma0 sandwiches before ordinary gamma identities.
    #[getter]
    pub(crate) fn gamma0(&self) -> bool {
        self.inner.gamma0
    }

    /// Whether to rewrite complex-conjugated gamma matrices before gamma0 reduction.
    #[getter]
    pub(crate) fn conjugate(&self) -> bool {
        self.inner.conjugate
    }

    /// Whether three four-dimensional gammas expand into a gamma5-epsilon basis.
    #[getter]
    pub(crate) fn expand_three_gamma_epsilon(&self) -> bool {
        self.inner.expand_three_gamma_epsilon
    }
}

/// Immutable configuration for color simplification.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "ColorSimplifySettings",
    module = "symbolica.community.spenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PyColorSimplifySettings {
    inner: ColorSimplifySettings,
}

impl PyColorSimplifySettings {
    pub(crate) fn from_rust(inner: ColorSimplifySettings) -> Self {
        Self { inner }
    }

    pub(crate) fn rust(&self) -> ColorSimplifySettings {
        self.inner
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyColorSimplifySettings {
    fn __repr__(self_: pyo3::PyRef<'_, Self>) -> PyResult<String> {
        let py = self_.py();
        let object = pyo3::IntoPyObject::into_pyobject(self_, py)?;
        crate::display::constructor_repr(
            object.as_any(),
            &[
                ("evaluate_traces", "evaluate_traces"),
                ("expand_cross_chain_fierz", "expand_cross_chain_fierz"),
                (
                    "substitute_cof_dimension_invariants",
                    "substitute_cof_dimension_invariants",
                ),
            ],
        )
    }

    /// Configure color-trace evaluation, cross-chain Fierz expansion, and invariant substitution.
    #[new]
    #[pyo3(signature = (
        *,
        evaluate_traces = true,
        expand_cross_chain_fierz = true,
        substitute_cof_dimension_invariants = false
    ))]
    pub(crate) fn new(
        evaluate_traces: bool,
        expand_cross_chain_fierz: bool,
        substitute_cof_dimension_invariants: bool,
    ) -> Self {
        Self {
            inner: ColorSimplifySettings {
                evaluate_traces,
                expand_cross_chain_fierz,
                substitute_cof_dimension_invariants,
                ..ColorSimplifySettings::default()
            },
        }
    }

    /// Whether closed color chains are evaluated as traces.
    #[getter]
    pub(crate) fn evaluate_traces(&self) -> bool {
        self.inner.evaluate_traces
    }

    /// Whether generators on different open chains or traces are expanded with the Fierz identity.
    #[getter]
    pub(crate) fn expand_cross_chain_fierz(&self) -> bool {
        self.inner.expand_cross_chain_fierz
    }

    /// Whether supported `cof(N)` invariants are replaced by explicit dimension formulas.
    #[getter]
    pub(crate) fn substitute_cof_dimension_invariants(&self) -> bool {
        self.inner.substitute_cof_dimension_invariants
    }
}

/// Immutable configuration for rewriting color invariants into a Casimir basis.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "ColorCasimirSettings",
    module = "symbolica.community.spenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PyColorCasimirSettings {
    inner: ColorCasimirSettings,
}

impl PyColorCasimirSettings {
    pub(crate) fn rust(&self) -> ColorCasimirSettings {
        self.inner
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyColorCasimirSettings {
    fn __repr__(self_: pyo3::PyRef<'_, Self>) -> PyResult<String> {
        let py = self_.py();
        let object = pyo3::IntoPyObject::into_pyobject(self_, py)?;
        crate::display::constructor_repr(
            object.as_any(),
            &[
                (
                    "rewrite_fundamental_dimension",
                    "rewrite_fundamental_dimension",
                ),
                (
                    "substitute_fundamental_index",
                    "substitute_fundamental_index",
                ),
            ],
        )
    }

    /// Configure the SU(N) dimension and fundamental-index normalizations used by Casimir rewriting.
    #[new]
    #[pyo3(signature = (
        *,
        rewrite_fundamental_dimension = true,
        substitute_fundamental_index = false
    ))]
    pub(crate) fn new(
        rewrite_fundamental_dimension: bool,
        substitute_fundamental_index: bool,
    ) -> Self {
        Self {
            inner: ColorCasimirSettings {
                rewrite_fundamental_dimension,
                substitute_fundamental_index,
            },
        }
    }

    /// Whether the fundamental dimension is rewritten with the SU(N) relation `d_F = C_A`.
    #[getter]
    pub(crate) fn rewrite_fundamental_dimension(&self) -> bool {
        self.inner.rewrite_fundamental_dimension
    }

    /// Whether the fundamental Dynkin index is replaced by `T_F = 1/2`.
    #[getter]
    pub(crate) fn substitute_fundamental_index(&self) -> bool {
        self.inner.substitute_fundamental_index
    }
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyGammaChainOrdering>()?;
    module.add_class::<PyGammaSimplifySettings>()?;
    module.add_class::<PyColorSimplifySettings>()?;
    module.add_class::<PyColorCasimirSettings>()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_settings_match_rust_defaults() {
        assert_eq!(
            PyGammaSimplifySettings::new(
                PyGammaChainOrdering::RepeatedPairs,
                true,
                false,
                false,
                false,
                GammaOutput::Reduced,
            )
            .rust(),
            GammaSimplifySettings::default()
        );
        let conjugated = PyGammaSimplifySettings::new(
            PyGammaChainOrdering::RepeatedPairs,
            false,
            true,
            true,
            false,
            GammaOutput::Reduced,
        );
        assert!(conjugated.gamma0());
        assert!(conjugated.conjugate());
        assert!(!conjugated.evaluate_traces());
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
