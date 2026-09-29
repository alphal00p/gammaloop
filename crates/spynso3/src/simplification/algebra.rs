use idenso::{
    color::ColorSimplifySettings,
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

/// Ordering strategy for open Dirac-matrix chains.
///
/// ``RepeatedPairs`` moves matching matrices together while preserving unrelated
/// order where possible. ``Canonical`` orders the complete chain using Clifford
/// anticommutation, which may introduce additional metric terms.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import GammaChainOrdering, GammaSimplifySettings
/// >>> settings = GammaSimplifySettings(chain_ordering=GammaChainOrdering.Canonical)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    frozen,
    from_py_object,
    eq,
    eq_int,
    name = "GammaChainOrdering",
    module = "symbolica.community.tensor"
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

/// Options for Dirac Clifford-algebra simplification.
///
/// Instances are immutable. Pass a configured instance to the corresponding
/// tensor operation; its properties expose the selected values.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import GammaSimplifySettings
/// >>> settings = GammaSimplifySettings()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "GammaSimplifySettings",
    module = "symbolica.community.tensor"
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
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyGammaSimplifySettings {
    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.GammaSimplifySettings()
    /// >>> text = repr(settings)
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

    /// Configure options for Dirac Clifford-algebra simplification.
    ///
    /// Parameters
    /// ----------
    /// output : str, default "reduced"
    ///     "reduced" applies identities. "chains" only joins chains, without
    ///     Clifford reduction or trace evaluation.
    /// chain_ordering : GammaChainOrdering or None, default None
    ///     Ordering of open Dirac chains. None selects RepeatedPairs.
    /// evaluate_traces : bool, default True
    ///     Evaluate closed Dirac traces when output is "reduced".
    /// gamma0 : bool, default False
    ///     Reduce four-dimensional gamma^0 sandwiches before ordinary identities.
    /// conjugate : bool, default False
    ///     Rewrite complex-conjugated Dirac matrices before gamma^0 reduction.
    /// expand_three_gamma_epsilon : bool, default False
    ///     Expand three four-dimensional gamma matrices into a gamma^5/epsilon basis.
    ///
    /// Returns
    /// -------
    /// GammaSimplifySettings
    ///     Immutable settings with the requested options.
    ///
    /// Notes
    /// -----
    /// Four-dimensional identities require explicit four-dimensional index spaces;
    /// these settings do not define a gamma^5 scheme in noninteger dimension.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import GammaSimplifySettings
    /// >>> settings = GammaSimplifySettings()
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

    /// "reduced" applies identities. "chains" only joins chains, without
    ///
    /// Returns
    /// -------
    /// str
    ///     "reduced" applies identities. "chains" only joins chains, without
    ///     Clifford reduction or trace evaluation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import GammaSimplifySettings
    /// >>> settings = GammaSimplifySettings()
    /// >>> value = settings.output
    #[getter]
    #[gen_stub(override_return_type(type_repr="typing.Literal['reduced', 'chains']", imports=("typing")))]
    pub(crate) fn output(&self) -> &'static str {
        match self.inner.output {
            GammaOutput::Reduced => "reduced",
            GammaOutput::Chains => "chains",
        }
    }

    /// Select repeated-pair reduction and evaluate traces.
    ///
    /// Returns
    /// -------
    /// GammaSimplifySettings
    ///     Reduced output with the selected ordering; optional gamma^0,
    ///     conjugation, and three-gamma epsilon passes remain disabled.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import GammaSimplifySettings
    /// >>> settings = GammaSimplifySettings.repeated_pairs()
    /// >>> settings.evaluate_traces
    /// True
    #[staticmethod]
    pub(crate) fn repeated_pairs() -> Self {
        Self {
            inner: GammaSimplifySettings::repeated_pairs(),
        }
    }

    /// Select canonical open-chain ordering and evaluate traces.
    ///
    /// Returns
    /// -------
    /// GammaSimplifySettings
    ///     Reduced output with the selected ordering; optional gamma^0,
    ///     conjugation, and three-gamma epsilon passes remain disabled.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import GammaSimplifySettings
    /// >>> settings = GammaSimplifySettings.canonical()
    /// >>> settings.evaluate_traces
    /// True
    #[staticmethod]
    pub(crate) fn canonical() -> Self {
        Self {
            inner: GammaSimplifySettings::canonical(),
        }
    }

    /// Ordering of open Dirac chains. None selects RepeatedPairs.
    ///
    /// Returns
    /// -------
    /// GammaChainOrdering
    ///     Ordering of open Dirac chains. None selects RepeatedPairs.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import GammaSimplifySettings
    /// >>> settings = GammaSimplifySettings()
    /// >>> value = settings.chain_ordering
    #[getter]
    pub(crate) fn chain_ordering(&self) -> PyGammaChainOrdering {
        self.inner.chain_ordering.into()
    }

    /// Evaluate closed Dirac traces when output is "reduced".
    ///
    /// Returns
    /// -------
    /// bool
    ///     Evaluate closed Dirac traces when output is "reduced".
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import GammaSimplifySettings
    /// >>> settings = GammaSimplifySettings()
    /// >>> value = settings.evaluate_traces
    #[getter]
    pub(crate) fn evaluate_traces(&self) -> bool {
        self.inner.evaluate_traces
    }

    /// Reduce four-dimensional gamma^0 sandwiches before ordinary identities.
    ///
    /// Returns
    /// -------
    /// bool
    ///     Reduce four-dimensional gamma^0 sandwiches before ordinary identities.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import GammaSimplifySettings
    /// >>> settings = GammaSimplifySettings()
    /// >>> value = settings.gamma0
    #[getter]
    pub(crate) fn gamma0(&self) -> bool {
        self.inner.gamma0
    }

    /// Rewrite complex-conjugated Dirac matrices before gamma^0 reduction.
    ///
    /// Returns
    /// -------
    /// bool
    ///     Rewrite complex-conjugated Dirac matrices before gamma^0 reduction.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import GammaSimplifySettings
    /// >>> settings = GammaSimplifySettings()
    /// >>> value = settings.conjugate
    #[getter]
    pub(crate) fn conjugate(&self) -> bool {
        self.inner.conjugate
    }

    /// Expand three four-dimensional gamma matrices into a gamma^5/epsilon basis.
    ///
    /// Returns
    /// -------
    /// bool
    ///     Expand three four-dimensional gamma matrices into a gamma^5/epsilon basis.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import GammaSimplifySettings
    /// >>> settings = GammaSimplifySettings()
    /// >>> value = settings.expand_three_gamma_epsilon
    #[getter]
    pub(crate) fn expand_three_gamma_epsilon(&self) -> bool {
        self.inner.expand_three_gamma_epsilon
    }
}

/// Options for SU(N) color-algebra simplification.
///
/// Instances are immutable. Pass a configured instance to the corresponding
/// tensor operation; its properties expose the selected values.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import ColorSimplifySettings
/// >>> settings = ColorSimplifySettings()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "ColorSimplifySettings",
    module = "symbolica.community.tensor"
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
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyColorSimplifySettings {
    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.ColorSimplifySettings()
    /// >>> text = repr(settings)
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

    /// Configure options for SU(N) color-algebra simplification.
    ///
    /// Parameters
    /// ----------
    /// evaluate_traces : bool, default True
    ///     Evaluate supported closed color traces.
    /// expand_cross_chain_fierz : bool, default True
    ///     Apply Fierz identities to contracted generators on separate open
    ///     chains or traces; this may introduce sums.
    /// substitute_cof_dimension_invariants : bool, default False
    ///     Replace supported fundamental SU(N) representation invariants by
    ///     explicit dimension formulas with T_R = 1/2.
    ///
    /// Returns
    /// -------
    /// ColorSimplifySettings
    ///     Immutable settings with the requested options.
    ///
    /// Notes
    /// -----
    /// With invariant substitution disabled, Casimir and Dynkin invariants remain
    /// symbolic wherever possible. Dimensions are read from the tensor slots.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import ColorSimplifySettings
    /// >>> settings = ColorSimplifySettings()
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
            },
        }
    }

    /// Evaluate supported closed color traces.
    ///
    /// Returns
    /// -------
    /// bool
    ///     Evaluate supported closed color traces.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import ColorSimplifySettings
    /// >>> settings = ColorSimplifySettings()
    /// >>> value = settings.evaluate_traces
    #[getter]
    pub(crate) fn evaluate_traces(&self) -> bool {
        self.inner.evaluate_traces
    }

    /// Apply Fierz identities to contracted generators on separate open
    ///
    /// Returns
    /// -------
    /// bool
    ///     Apply Fierz identities to contracted generators on separate open
    ///     chains or traces; this may introduce sums.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import ColorSimplifySettings
    /// >>> settings = ColorSimplifySettings()
    /// >>> value = settings.expand_cross_chain_fierz
    #[getter]
    pub(crate) fn expand_cross_chain_fierz(&self) -> bool {
        self.inner.expand_cross_chain_fierz
    }

    /// Replace supported fundamental SU(N) representation invariants by
    ///
    /// Returns
    /// -------
    /// bool
    ///     Replace supported fundamental SU(N) representation invariants by
    ///     explicit dimension formulas with T_R = 1/2.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import ColorSimplifySettings
    /// >>> settings = ColorSimplifySettings()
    /// >>> value = settings.substitute_cof_dimension_invariants
    #[getter]
    pub(crate) fn substitute_cof_dimension_invariants(&self) -> bool {
        self.inner.substitute_cof_dimension_invariants
    }
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyGammaChainOrdering>()?;
    module.add_class::<PyGammaSimplifySettings>()?;
    module.add_class::<PyColorSimplifySettings>()?;

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
    }
}
