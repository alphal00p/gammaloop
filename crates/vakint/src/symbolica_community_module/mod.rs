use symbolica::api::python::{Citation, SymbolicaCommunityModule};
mod cache;
mod diagram_integral;
mod inspection;
mod numerical;
mod presentation;
mod resources;

use std::collections::{BTreeMap, HashMap};
use std::env;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use pyo3::types::{PyAnyMethods, PyDict, PyMapping, PyModule, PyModuleMethods, PyType};
use pyo3::{
    Bound, FromPyObject, PyAny, PyClass, PyErr, PyRef, Python, exceptions, pyclass, pymethods,
};
use pyo3::{Py, PyResult};
use symbolica::api::python::PythonExpression;
use symbolica::atom::{Atom, AtomView, Symbol};
use symbolica::domains::float::{Complex, Float, RealLike};

use crate::presentation::Present;
use crate::symbols::S;
use crate::{
    AlphaLoopOptions, EvaluationMethod, EvaluationOrder, FMFTOptions, LoopNormalizationFactor,
    MATADOptions, NumericalEvaluationResult, PySecDecOptions, RustRedEvaluationError,
    RustRedEvaluationOptions, TensorReductionMethod, Vakint, VakintError, VakintExpression,
    VakintSettings, vakint_symbol,
};
use inspection::{
    ConstantSymbolWrapper, MasterSymbolWrapper, NumericalCoefficientWrapper,
    PrecisionWarningWrapper, VakintConstantWrapper, VakintMasterWrapper,
};

use pyo3::types::PyListMethods;
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::*;

impl SymbolicaCommunityModule for VakintWrapper {
    fn get_citations() -> Vec<Citation> {
        crate::citations::USED_CITATIONS.get_citations()
    }

    fn get_name() -> String {
        "hepkit_vakint".to_string()
    }

    fn register_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
        initialize_vakint(m)
    }

    fn initialize(_py: Python) -> PyResult<()> {
        let _ = S.dot;
        Ok(())
    }
}

trait ModuleInit: PyClass {
    fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<Self>()
    }
}
impl ModuleInit for VakintWrapper {}
impl ModuleInit for NumericalEvaluationResultWrapper {}
impl ModuleInit for VakintExpressionWrapper {}
impl ModuleInit for VakintEvaluationMethodWrapper {}

macro_rules! define_vakint_python_surface {
    ($($class:ty),+ $(,)?) => {
        pub(crate) fn initialize_vakint(m: &Bound<'_, PyModule>) -> PyResult<()> {
            $(<$class as ModuleInit>::init(m)?;)+
            m.add_function(pyo3::wrap_pyfunction!(diagram_integral::integral_from_diagram, m)?)?;
            Ok(())
        }

        /// Public names registered on `symbolica.community.hepkit.vakint`.
        #[cfg(feature = "python_stubgen")]
        pub const PYTHON_STUB_SURFACE: &[&str] = &[
            $(<$class as PyClass>::NAME,)+
            "integral_from_diagram",
        ];
    };
}

define_vakint_python_surface! {
    VakintWrapper,
    cache::RunningCacheControlWrapper,
    cache::CacheFileOptionsWrapper,
    NumericalEvaluationResultWrapper,
    VakintExpressionWrapper,
    VakintEvaluationMethodWrapper,
    inspection::AvailablePrecisionWrapper,
    inspection::ErrorQualityWrapper,
    inspection::ComponentPrecisionWrapper,
    inspection::CoefficientPrecisionWrapper,
    inspection::PrecisionSourceWrapper,
    inspection::PrecisionWarningWrapper,
    inspection::DataProvenanceWrapper,
    inspection::NumericalCoefficientWrapper,
    inspection::ConstantSymbolWrapper,
    inspection::MasterSymbolWrapper,
    inspection::MasterSymbolicExpressionWrapper,
    inspection::VakintConstantWrapper,
    inspection::MasterCoefficientWrapper,
    inspection::VakintMasterWrapper,
}

fn vakint_to_python_error(vakint_error: VakintError) -> PyErr {
    exceptions::PyValueError::new_err(format!("Vakint error | {vakint_error}"))
}

fn rustred_cache_policy(
    value: &str,
) -> Result<rustred::solver::CandidateRoutedCachePolicy, String> {
    use rustred::solver::CandidateRoutedCachePolicy;
    match value {
        "evict-lru" => Ok(CandidateRoutedCachePolicy::EvictLeastRecentlyUsed),
        "evict-aged-lfu" => Ok(CandidateRoutedCachePolicy::EvictAgedLeastFrequentlyUsed),
        "fail-on-limit" => Ok(CandidateRoutedCachePolicy::FailOnLimit),
        _ => Err(format!(
            "unknown rustred_cache_policy {value:?}; expected 'evict-lru', 'evict-aged-lfu' or 'fail-on-limit'"
        )),
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "Vakint", module = "symbolica.community.hepkit.vakint")]
/// Vakint engine and settings used for matching, reduction, and evaluation.
///
/// Construct one instance and reuse it. Initialization prepares lower-loop data;
/// the five-loop catalogue is attached only after successful artifact loading.
pub struct VakintWrapper {
    pub vakint: Vakint,
    pub settings: VakintSettings,
}

impl VakintWrapper {
    /// Stage publication separately from the current wrapper's mutable application slot.
    fn staged_engine(&self) -> PyResult<Vakint> {
        let loaded = self
            .vakint
            .rustred_application
            .lock()
            .map_err(|_| vakint_to_python_error(RustRedEvaluationError::CacheLockPoisoned.into()))?
            .clone();
        let mut staged = self.vakint.clone();
        staged.rustred_application = Arc::new(Mutex::new(loaded));
        Ok(staged)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "VakintNumericalResult",
    module = "symbolica.community.hepkit.vakint",
    frozen
)]
/// Numerical Laurent series in the dimensional-regularization parameter epsilon.
pub struct NumericalEvaluationResultWrapper {
    pub value: NumericalEvaluationResult,
}

impl<'a, 'py> FromPyObject<'a, 'py> for NumericalEvaluationResultWrapper {
    type Error = PyErr;

    fn extract(ob: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(a) = ob.cast::<NumericalEvaluationResultWrapper>() {
            Ok(NumericalEvaluationResultWrapper {
                value: a.borrow().value.clone(),
            })
        } else {
            Err(exceptions::PyValueError::new_err(
                "Not a valid Vakint numerical result (Laurent series in epsilon)",
            ))
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(record_usage)]
#[pymethods]
impl NumericalEvaluationResultWrapper {
    /// Return a Symbolica expression with independently achieved component precision.
    ///
    /// The default variable is `vakint::ε`. Use `Vakint.numerical_result_to_expression`
    /// when the originating engine has a custom epsilon symbol.
    /// Raises ValueError for uncertain zeros or components without reliably encodable digits.
    #[pyo3(signature = (epsilon_symbol = None))]
    pub fn to_expression(
        &self,
        epsilon_symbol: Option<PythonExpression>,
    ) -> PyResult<PythonExpression> {
        let epsilon = match epsilon_symbol.as_ref().map(|e| e.expr.as_view()) {
            None => vakint_symbol!("ε"),
            Some(AtomView::Var(v)) => v.get_symbol(),
            Some(_) => {
                return Err(exceptions::PyValueError::new_err(
                    "epsilon_symbol must be a Symbolica variable",
                ));
            }
        };
        self.value
            .to_atom(epsilon)
            .map(Into::into)
            .map_err(vakint_to_python_error)
    }

    /// Render raw central values, achieved component precision, and warnings.
    pub fn to_html(&self, py: Python<'_>) -> PyResult<String> {
        presentation::render(py, &self.value.presentation(), None)
    }

    /// Display remains available even when checked expression conversion fails.
    pub fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        self.to_html(py)
    }

    /// Read-only records keyed by epsilon power; absent orders are not fabricated.
    #[getter]
    pub fn coefficients(&self) -> BTreeMap<i64, NumericalCoefficientWrapper> {
        self.value
            .coefficients()
            .iter()
            .map(|(power, coefficient)| (*power, coefficient.clone().into()))
            .collect()
    }

    #[getter]
    pub fn requested_decimal_precision(&self) -> u32 {
        self.value.requested_decimal_precision()
    }

    #[getter]
    pub fn warnings(&self) -> Vec<PrecisionWarningWrapper> {
        self.value
            .warnings()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    pub fn __repr__(&self) -> String {
        self.value.to_string()
    }

    /// Return Symbolica's `FormattedOutput`; keyword options are forwarded unchanged.
    #[gen_stub(override_return_type(
        type_repr = "symbolica.FormattedOutput",
        imports = ("symbolica")
    ))]
    #[pyo3(signature = (**options))]
    pub fn formatted(
        &self,
        py: Python<'_>,
        options: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Py<PyAny>> {
        Ok(Py::new(py, self.to_expression(None)?)?
            .bind(py)
            .call_method("formatted", (), options)?
            .unbind())
    }

    /// String representation of the numerical result.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.hepkit.vakint import VakintNumericalResult
    /// >>> result = VakintNumericalResult([
    /// ...     (-3, (0.0, -11440.53140354612)),
    /// ...     (-2, (0.0, 57169.95521898031)),
    /// ...     (-1, (0.0, -178748.9838377694)),
    /// ...     (0, (0.0, 321554.1122184795)),
    /// ... ])
    /// >>> "Vakint numerical result" in str(result)
    /// True
    /// ```
    pub fn __str__(&self) -> PyResult<String> {
        Ok(format!("{}", self.value))
    }

    /// Lossy conversion to Python float tuples; use coefficients to retain precision metadata.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.hepkit.vakint import VakintNumericalResult
    /// >>> result = VakintNumericalResult([
    /// ...     (-3, (0.0, -11440.53140354612)),
    /// ...     (-2, (0.0, 57169.95521898031)),
    /// ...     (-1, (0.0, -178748.9838377694)),
    /// ...     (0, (0.0, 321554.1122184795)),
    /// ... ])
    ///
    /// >>> values = result.to_list()
    /// >>> values[0]
    /// (-3, (0.0, -11440.53140354612))
    /// ```
    pub fn to_list(&self) -> PyResult<Vec<(i64, (f64, f64))>> {
        Ok(self
            .value
            .coefficients()
            .iter()
            .map(|(exp, coeff)| (*exp, (coeff.value().re.to_f64(), coeff.value().im.to_f64())))
            .collect())
    }
    #[new]
    /// Create a numerical Laurent series from `(epsilon exponent, (real, imaginary))` tuples.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.hepkit.vakint import VakintNumericalResult
    /// >>> result = VakintNumericalResult([
    /// ...     (-3, (0.0, -11440.53140354612)),
    /// ...     (-2, (0.0, 57169.95521898031)),
    /// ...     (-1, (0.0, -178748.9838377694)),
    /// ...     (0, (0.0, 321554.1122184795)),
    /// ... ])
    /// >>> len(result.to_list())
    /// 4
    /// ```
    ///
    /// Parameters
    /// ----------
    ///
    /// values : List[Tuple[int, Tuple[float, float]]]
    ///    A list of tuples, each containing an integer exponent of epsilon and a tuple of two floats
    ///    representing the real and imaginary parts of the coefficient.
    pub fn __new__(values: Vec<(i64, (f64, f64))>) -> PyResult<NumericalEvaluationResultWrapper> {
        let r = NumericalEvaluationResult::from_values(
            values
                .iter()
                .map(|(eps_pwr, (re, im))| {
                    (
                        *eps_pwr,
                        Complex::new(Float::with_val(53, re), Float::with_val(53, im)),
                    )
                })
                .collect::<Vec<_>>(),
            16,
        );
        Ok(NumericalEvaluationResultWrapper { value: r })
    }

    #[pyo3(signature = (other, relative_threshold, error = None, max_pull = None))]
    /// Compare this numerical result to another, returning a tuple of (bool, str) where the bool indicates whether the results match within the specified thresholds,
    /// and the str provides details of the comparison.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.hepkit.vakint import VakintNumericalResult
    /// >>> result1 = VakintNumericalResult([
    /// ...     (-3, (0.0, -11440.53140354612)),
    /// ... ])
    /// >>> result2 = VakintNumericalResult([
    /// ...     (-3, (0.0, -11440.53140354612)),
    /// ...     (-2, (0.0, 2.0)),
    /// ... ])
    /// >>> matches, details = result1.compare_to(result2, relative_threshold=1e-5)
    /// >>> matches
    /// False
    /// >>> "ε^-2" in details
    /// True
    /// ```
    ///
    /// Parameters
    /// ----------
    ///
    /// other : VakintNumericalResult
    ///    The other numerical result to compare to.
    /// relative_threshold : float
    ///    The relative threshold for comparison.
    /// error : Optional[VakintNumericalResult]
    ///    An optional numerical result representing the error in the evaluation.
    /// max_pull : Optional[float]
    ///    The maximum pull for comparison. Default is 3.0.
    pub fn compare_to(
        &self,
        other: PyRef<NumericalEvaluationResultWrapper>,
        relative_threshold: f64,
        error: Option<PyRef<NumericalEvaluationResultWrapper>>,
        max_pull: Option<f64>,
    ) -> PyResult<(bool, String)> {
        Ok(self.value.does_approx_match(
            &other.value,
            error.map(|e| e.value.clone()).as_ref(),
            relative_threshold,
            max_pull.unwrap_or(3.0),
        ))
    }
}

/// A Vakint integral split into its numerator and normalized topology structure.
///
/// Construct this wrapper from a Symbolica expression before applying Vakint operations.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "VakintExpression",
    module = "symbolica.community.hepkit.vakint"
)]
pub struct VakintExpressionWrapper {
    pub value: VakintExpression,
}

impl<'a, 'py> FromPyObject<'a, 'py> for VakintExpressionWrapper {
    type Error = PyErr;

    fn extract(ob: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(a) = ob.cast::<VakintExpressionWrapper>() {
            Ok(VakintExpressionWrapper {
                value: a.borrow().value.clone(),
            })
        } else {
            Err(exceptions::PyValueError::new_err(
                "Not a valid vakint expression",
            ))
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(record_usage)]
#[pymethods]
impl VakintExpressionWrapper {
    /// Display through Symbolica's native rich expression printer.
    pub fn _repr_html_(&self) -> PyResult<String> {
        self.to_expression()?._repr_html_()
    }

    /// Return Symbolica's `FormattedOutput`; keyword options are forwarded unchanged.
    #[gen_stub(override_return_type(
        type_repr = "symbolica.FormattedOutput",
        imports = ("symbolica")
    ))]
    #[pyo3(signature = (**options))]
    pub fn formatted(
        &self,
        py: Python<'_>,
        options: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Py<PyAny>> {
        Ok(Py::new(py, self.to_expression()?)?
            .bind(py)
            .call_method("formatted", (), options)?
            .unbind())
    }

    /// String representation of the VakintExpression.
    ///
    pub fn __str__(&self) -> PyResult<String> {
        Ok(format!("{}", self.value))
    }

    /// Convert the VakintExpression back to a Symbolica Expression.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.hepkit.vakint import VakintExpression
    /// >>> integral = VakintExpression(E('''
    /// ...     k(1,11)*k(1,11)
    /// ...     *topo(prop(1,edge(1,1),k(1),muvsq,1))
    /// ... ''', default_namespace="vakint"))
    /// >>> "topo(" in str(integral.to_expression())
    /// True
    /// ```
    pub fn to_expression(&self) -> PyResult<PythonExpression> {
        let a: Atom = self.value.clone().into();
        Ok(a.into())
    }

    #[new]
    /// Split a Symbolica expression into Vakint numerator and topology components.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.hepkit.vakint import VakintExpression
    /// >>> integral = E('''
    /// ...     (
    /// ...         k(1,11)*k(2,11)*k(1,22)*k(2,22)
    /// ...       + p(1,11)*k(3,11)*k(3,22)*p(2,22)
    /// ...       + p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)
    /// ...     )*topo(
    /// ...          prop(1,edge(1,2),k(1),muvsq,1)
    /// ...         * prop(2,edge(2,3),k(2),muvsq,1)
    /// ...         * prop(3,edge(3,1),k(3),muvsq,1)
    /// ...         * prop(4,edge(1,4),k(3)-k(1),muvsq,1)
    /// ...         * prop(5,edge(2,4),k(1)-k(2),muvsq,1)
    /// ...         * prop(6,edge(3,4),k(2)-k(3),muvsq,1)
    /// ...     )
    /// ... ''', default_namespace="vakint")
    /// >>> wrapped = VakintExpression(integral)
    /// >>> "topo(" in str(wrapped)
    /// True
    /// ```
    ///
    /// Parameters
    /// ----------
    ///
    /// atom : Expression
    ///   A Symbolica Expression containing a vakint integral, i.e. a sum of terms, each a product of a numerator and a `vakint::topo(...)` structure.
    pub fn new(atom: Py<PyAny>) -> PyResult<VakintExpressionWrapper> {
        Python::attach(|py| {
            let a = atom.extract::<PythonExpression>(py)?;
            Ok(VakintExpressionWrapper {
                value: VakintExpression::try_from(a.expr).map_err(vakint_to_python_error)?,
            })
        })
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "VakintEvaluationMethod",
    module = "symbolica.community.hepkit.vakint"
)]
/// One configured backend in a `Vakint` instance's evaluation order.
pub struct VakintEvaluationMethodWrapper {
    pub method: EvaluationMethod,
}

impl<'a, 'py> FromPyObject<'a, 'py> for VakintEvaluationMethodWrapper {
    type Error = PyErr;

    fn extract(ob: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(a) = ob.cast::<VakintEvaluationMethodWrapper>() {
            Ok(VakintEvaluationMethodWrapper {
                method: a.borrow().method.clone(),
            })
        } else {
            Err(exceptions::PyValueError::new_err(
                "Not a valid vakint evaluation method",
            ))
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(record_usage)]
#[pymethods]
impl VakintEvaluationMethodWrapper {
    /// String representation of the evaluation method.
    pub fn __str__(&self) -> PyResult<String> {
        Ok(format!("{}", self.method))
    }

    #[classmethod]
    /// Create a new VakintEvaluationMethod instance representing the AlphaLoop method.
    /// This method does not take any parameters.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.hepkit.vakint import VakintEvaluationMethod
    /// >>> alphaloop_method = VakintEvaluationMethod.new_alphaloop_method()
    /// >>> "AlphaLoop" in str(alphaloop_method)
    /// True
    /// ```
    pub fn new_alphaloop_method(
        _cls: &Bound<'_, PyType>,
    ) -> PyResult<VakintEvaluationMethodWrapper> {
        Ok(VakintEvaluationMethodWrapper {
            method: EvaluationMethod::AlphaLoop(AlphaLoopOptions::default()),
        })
    }

    #[pyo3(signature = (master_mode = "numerical"))]
    #[classmethod]
    /// Create the opt-in FORM-independent RustRed scalar-integral evaluation method.
    ///
    /// This method leaves Vakint's tensor-reduction selection unchanged. It uses
    /// shipped closing artifacts for supported one- through three-loop equal-mass
    /// vacuum families and guarded, finite-tested candidate programs at four loops.
    /// Four-loop success is pointwise reduction, not arbitrary-index closure.
    /// The default FeynKit tensor prepass is also FORM-independent; explicitly
    /// choosing the historical AlphaLoop tensor prepass still requires FORM.
    ///
    /// ## Examples
    /// ```python
    /// rustred_method = VakintEvaluationMethod.new_rustred_method(
    ///     master_mode="numerical"
    /// )
    /// ```
    ///
    /// Parameters
    /// ----------
    ///
    /// master_mode : str
    ///    "retain", "symbolic", or "numerical" (default).
    pub fn new_rustred_method(
        _cls: &Bound<'_, PyType>,
        master_mode: &str,
    ) -> PyResult<VakintEvaluationMethodWrapper> {
        Ok(VakintEvaluationMethodWrapper {
            method: EvaluationMethod::RustRed(RustRedEvaluationOptions {
                master_mode: match master_mode {
                    "retain" => crate::RustRedMasterMode::Retain,
                    "symbolic" => crate::RustRedMasterMode::Symbolic,
                    "numerical" => crate::RustRedMasterMode::Numerical,
                    _ => {
                        return Err(pyo3::exceptions::PyValueError::new_err(
                            "master_mode must be retain, symbolic, or numerical",
                        ));
                    }
                },
            }),
        })
    }

    #[pyo3(signature = (expand_masters = None, susbstitute_masters = None, substitute_hpls = None, direct_numerical_substition = None))]
    #[classmethod]
    /// Create a new VakintEvaluationMethod instance representing the MATAD method.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.hepkit.vakint import VakintEvaluationMethod
    /// >>> matad_method = VakintEvaluationMethod.new_matad_method(
    /// ...     expand_masters=True,
    /// ...     susbstitute_masters=True,
    /// ...     substitute_hpls=True,
    /// ...     direct_numerical_substition=True,
    /// ... )
    /// >>> "MATAD" in str(matad_method)
    /// True
    /// ```
    ///
    /// Parameters
    /// ----------
    ///
    /// expand_masters : Optional[bool]
    ///    Whether to expand master integrals. Default is True.
    /// susbstitute_masters : Optional[bool]
    ///    Whether to substitute master integrals. Default is True.
    /// substitute_hpls : Optional[bool]
    ///    Whether to substitute harmonic polylogarithms. Default is True.
    /// direct_numerical_substition : Optional[bool]
    ///    Whether to perform direct numerical substitution. Default is True.
    ///
    /// Notes
    /// -----
    /// `susbstitute_masters` and `direct_numerical_substition` retain their historical
    /// misspellings for API compatibility. They mean `substitute_masters` and
    /// `direct_numerical_substitution`, respectively.
    pub fn new_matad_method(
        _cls: &Bound<'_, PyType>,
        expand_masters: Option<bool>,
        susbstitute_masters: Option<bool>,
        substitute_hpls: Option<bool>,
        direct_numerical_substition: Option<bool>,
    ) -> PyResult<VakintEvaluationMethodWrapper> {
        Ok(VakintEvaluationMethodWrapper {
            method: EvaluationMethod::MATAD(MATADOptions {
                expand_masters: expand_masters.unwrap_or(true),
                susbstitute_masters: susbstitute_masters.unwrap_or(true),
                substitute_hpls: substitute_hpls.unwrap_or(true),
                direct_numerical_substition: direct_numerical_substition.unwrap_or(true),
            }),
        })
    }

    #[pyo3(signature = (expand_masters = None, susbstitute_masters = None))]
    #[classmethod]
    /// Create a new VakintEvaluationMethod instance representing the FMFT method.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.hepkit.vakint import VakintEvaluationMethod
    /// >>> fmft_method = VakintEvaluationMethod.new_fmft_method(
    /// ...     expand_masters=True,
    /// ...     susbstitute_masters=True,
    /// ... )
    /// >>> "FMFT" in str(fmft_method)
    /// True
    /// ```
    ///
    /// Parameters
    /// ----------
    ///
    /// expand_masters : Optional[bool]
    ///   Whether to expand master integrals. Default is True.
    /// susbstitute_masters : Optional[bool]
    ///   Whether to substitute master integrals. Default is True.
    ///
    /// Notes
    /// -----
    /// `susbstitute_masters` retains its historical misspelling for API compatibility;
    /// it means `substitute_masters`.
    pub fn new_fmft_method(
        _cls: &Bound<'_, PyType>,
        expand_masters: Option<bool>,
        susbstitute_masters: Option<bool>,
    ) -> PyResult<VakintEvaluationMethodWrapper> {
        Ok(VakintEvaluationMethodWrapper {
            method: EvaluationMethod::FMFT(FMFTOptions {
                expand_masters: expand_masters.unwrap_or(true),
                susbstitute_masters: susbstitute_masters.unwrap_or(true),
            }),
        })
    }

    #[pyo3(signature = (quiet = None, relative_precision = None, min_n_evals = None, max_n_evals = None, reuse_existing_output = None, numerical_masses = None, numerical_external_momenta = None))]
    #[allow(clippy::too_many_arguments)]
    #[classmethod]
    /// Create a new VakintEvaluationMethod instance representing the numerical pySecDec method.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.hepkit.vakint import VakintEvaluationMethod
    /// >>> pysecdec_method = VakintEvaluationMethod.new_pysecdec_method(
    /// ...     quiet=True,
    /// ...     relative_precision=1e-7,
    /// ...     min_n_evals=10_000,
    /// ...     max_n_evals=1_000_000_000_000,
    /// ...     reuse_existing_output=None,
    /// ...     numerical_masses={"muvsq": 1.0},
    /// ...     numerical_external_momenta={
    /// ...         1: (1.0, 0.0, 0.0, 0.0),
    /// ...         2: (0.0, 1.0, 0.0, 0.0),
    /// ...     },
    /// ... )
    /// >>> "PySecDec" in str(pysecdec_method)
    /// True
    /// ```
    ///
    /// pySecDec performs numerical evaluation, so every required mass and external momentum
    /// must have a numerical value. Constructing this method does not run pySecDec; evaluation
    /// requires a working Python/pySecDec installation.
    ///
    /// Parameters
    /// ----------
    ///
    /// quiet : Optional[bool]
    ///    Whether to suppress output from pySecDec. Default is True.
    /// relative_precision : Optional[float]
    ///    The relative precision to be achieved in the numerical integration. Default is 1e-7.
    /// min_n_evals : Optional[int]
    ///    The minimum number of evaluations to be performed in the numerical integration. Default is 10,000.
    /// max_n_evals : Optional[int]
    ///    The maximum number of evaluations to be performed in the numerical integration. Default is 1,000,000,000,000.
    /// reuse_existing_output : Optional[str]
    ///    Path to existing pySecDec output to reuse. Default is None.
    /// numerical_masses : Optional[Dict[str, float]]
    ///    A dictionary mapping mass parameter names to their numerical values. Default is an empty dictionary.
    /// numerical_external_momenta : Optional[Dict[int, Tuple[float, float, float, float]]]
    ///    A dictionary mapping external momentum indices to their numerical 4-vector values. Default is an empty dictionary.
    pub fn new_pysecdec_method(
        _cls: &Bound<'_, PyType>,
        quiet: Option<bool>,
        relative_precision: Option<f64>,
        min_n_evals: Option<u64>,
        max_n_evals: Option<u64>,
        reuse_existing_output: Option<String>,
        numerical_masses: Option<HashMap<String, f64>>,
        numerical_external_momenta: Option<HashMap<usize, (f64, f64, f64, f64)>>,
    ) -> PyResult<VakintEvaluationMethodWrapper> {
        let ext_mom = if let Some(em) = numerical_external_momenta {
            HashMap::from_iter(em.iter().map(|(i_ext, ps)| (format!("p{}", i_ext), *ps)))
        } else {
            HashMap::new()
        };
        Ok(VakintEvaluationMethodWrapper {
            method: EvaluationMethod::PySecDec(PySecDecOptions {
                quiet: quiet.unwrap_or(true),
                reuse_existing_output,
                relative_precision: relative_precision.unwrap_or(1e-7),
                min_n_evals: min_n_evals.unwrap_or(10_000),
                max_n_evals: max_n_evals.unwrap_or(1_000_000_000_000),
                numerical_masses: numerical_masses.unwrap_or_default(),
                numerical_external_momenta: ext_mom,
            }),
        })
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(record_usage)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl VakintWrapper {
    /// Load a local external-data checkout without network access.
    #[pyo3(signature = (path, max_loop_count))]
    pub fn load_external_data(
        &mut self,
        path: std::path::PathBuf,
        max_loop_count: usize,
    ) -> PyResult<()> {
        let staged = resources::run(self.settings.rustred_resources, || {
            let mut staged = self.staged_engine()?;
            staged
                .load_external_data(path, max_loop_count)
                .map_err(vakint_to_python_error)?;
            Ok(staged)
        })?;
        self.vakint = staged;
        Ok(())
    }

    /// Load a local loop pack, or a sealed lower-loop native rules file.
    pub fn load_rustred_artifacts(
        &mut self,
        loop_count: usize,
        path: std::path::PathBuf,
    ) -> PyResult<()> {
        let staged = resources::run(self.settings.rustred_resources, || {
            let mut staged = self.staged_engine()?;
            staged
                .load_rustred_artifacts(loop_count, path)
                .map_err(vakint_to_python_error)?;
            Ok(staged)
        })?;
        self.vakint = staged;
        Ok(())
    }

    /// Install all official loop packs through the requested maximum, or load a
    /// supplied checkout directly. Existing validated installations work offline.
    /// Prepared native caches load by default. Set `load_prepared_caches=False`
    /// to activate the same validated loop packs without importing their caches.
    #[cfg(all(unix, not(target_arch = "wasm32")))]
    #[pyo3(signature = (max_loop_count, version="latest", local_artifacts_path=None, refresh=false, *, load_prepared_caches=true))]
    pub fn load_official_external_data<'py>(
        &mut self,
        py: Python<'py>,
        max_loop_count: usize,
        version: &str,
        local_artifacts_path: Option<std::path::PathBuf>,
        refresh: bool,
        load_prepared_caches: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        let (staged, result) = resources::run(self.settings.rustred_resources, || {
            let mut staged = self.staged_engine()?;
            let report = staged
                .load_official_external_data(
                    max_loop_count,
                    &crate::OfficialExternalDataOptions {
                        version: version.to_owned(),
                        local_artifacts_path,
                        refresh,
                        load_prepared_caches,
                        ..Default::default()
                    },
                )
                .map_err(vakint_to_python_error)?;
            let result = PyDict::new(py);
            result.set_item("version", report.version)?;
            result.set_item("path", report.path)?;
            result.set_item("load_prepared_caches", report.load_prepared_caches)?;
            let loops = pyo3::types::PyList::empty(py);
            for item in report.loops {
                let record = PyDict::new(py);
                record.set_item("loop_count", item.loop_count)?;
                record.set_item("path", item.path)?;
                record.set_item("cache_available", item.cache_available)?;
                record.set_item("cache_loaded", item.cache_loaded)?;
                record.set_item("cache_path", item.cache_path)?;
                record.set_item("cache_format", item.cache_format)?;
                record.set_item("cache_source", item.cache_source)?;
                record.set_item("source", item.source)?;
                loops.append(record)?;
            }
            result.set_item("loops", loops)?;
            Ok((staged, result))
        })?;
        self.vakint = staged;
        Ok(result)
    }

    #[pyo3(signature = (run_time_decimal_precision = None, evaluation_order = None, tensor_reduction_method = None, epsilon_symbol = None, mu_r_sq_symbol = None, form_exe_path = None, python_exe_path = None, verify_numerator_identification = None, integral_normalization_factor = None, allow_unknown_integrals = None, clean_tmp_dir = None, number_of_terms_in_epsilon_expansion = None, use_dot_product_notation = None, temporary_directory = None, rustred_artifact_path = None, rustred_cache_bytes = None, rustred_cache_terms = None, rustred_product_cache_bytes = None, rustred_product_cache_terms = None, rustred_cache_policy = "evict-lru", rustred_scalar_max_work = None, max_process_rss_bytes = None, max_elapsed_seconds = None))]
    #[allow(clippy::too_many_arguments)]
    #[new]
    /// Create a new Vakint instance, specifying details of the evaluation stack. Note that the same instance can be recycled across multiple evaluations.
    /// Note that the creation of a Vakint instance involves the processing and creation of the library of all known topologies, which can be time consuming.
    /// External executables are validated when an operation needs them, so the FeynKit tensor backend can be used on systems without FORM.
    ///
    /// ## Examples
    /// ```python
    /// vakint = Vakint(
    ///     integral_normalization_factor="MSbar",
    ///     mu_r_sq_symbol=S("mursq"),
    ///     # If you select 5 terms, then MATAD will be used, but for 4 and fewer, alphaLoop is will be used as
    ///     # it is first in the evaluation_order supplied.
    ///     number_of_terms_in_epsilon_expansion=4,
    ///     evaluation_order=[
    ///         VakintEvaluationMethod.new_alphaloop_method(),
    ///         VakintEvaluationMethod.new_matad_method(),
    ///         VakintEvaluationMethod.new_fmft_method(),
    ///         VakintEvaluationMethod.new_pysecdec_method(
    ///             min_n_evals=10_000,
    ///             max_n_evals=1000_000,
    ///             numerical_masses=masses,
    ///             numerical_external_momenta=external_momenta
    ///         ),
    ///     ],
    ///     tensor_reduction_method="feynkit",
    ///     form_exe_path="form",
    ///     python_exe_path="python3",
    /// )
    /// ```
    ///
    /// An empty evaluation order is appropriate for matching, canonicalization, and tensor
    /// reduction. Add explicit `VakintEvaluationMethod` entries before evaluating an integral;
    /// construction validates the executables required by those entries.
    ///
    /// Parameters
    /// ----------
    ///
    /// run_time_decimal_precision : Optional[int]
    ///     The decimal precision to be used during the evaluation. Default is 17.
    /// evaluation_order : Optional[Sequence[VakintEvaluationMethod]]
    ///     A list of `VakintEvaluationMethod` instances specifying the order in which evaluation methods are to be applied. Default is all available methods in a sensible order.
    /// tensor_reduction_method : Optional[str]
    ///     Numerator tensor-reduction backend: "feynkit" is the default, native backend and does not require FORM; "alphaloop" explicitly selects the historical FORM projector.
    /// epsilon_symbol : Optional[Expression]
    ///     The symbol to be used for the dimensional regularisation parameter epsilon. Default is "ε".
    /// mu_r_sq_symbol : Optional[Expression]
    ///     The symbol to be used for the renormalisation scale squared. Default is "mursq".
    /// form_exe_path : Optional[str]
    ///     The path to the FORM executable. Default is "form".
    /// python_exe_path : Optional[str]
    ///     The path to the Python executable. Default is "python3".
    /// verify_numerator_identification : Optional[bool]
    ///     Whether to verify the identification of numerator structures. Default is True.
    /// integral_normalization_factor : Optional[str]
    ///     The normalization factor to be used for integrals. Can be "MSbar", "pySecDec", "FMFTandMATAD" or a custom string. Default is "MSbar".
    /// allow_unknown_integrals : Optional[bool]
    ///     Whether to allow unknown integrals to be processed. Default is True.
    /// clean_tmp_dir : Optional[bool]
    ///     Whether to clean the temporary directory after evaluation. Default is True, unless the environment variable VAKINT_NO_CLEAN_TMP_DIR is set.
    /// number_of_terms_in_epsilon_expansion : Optional[int]
    ///     The number of terms in the epsilon expansion to be computed. Default is 4.
    /// use_dot_product_notation : Optional[bool]
    ///     Whether to use dot product notation for scalar products. Default is False.
    /// temporary_directory : Optional[str]
    ///     The path to the temporary directory to be used. Default is None, in which case a system temporary directory will be used.
    /// rustred_artifact_path : Optional[str]
    ///     Saved application package for exact symbolic five-loop RustRed reduction.
    ///     Use master_mode="retain" for reduction without coefficient evaluation.
    /// rustred_cache_bytes : Optional[int]
    ///     Saved-application coefficient-cache payload-byte allowance. None
    ///     imposes no payload-byte ceiling and does not enable count diagnostics.
    /// rustred_cache_terms : Optional[int]
    ///     Optional diagnostic aggregate coefficient-term ceiling. Supplying a
    ///     count ceiling explicitly enables the corresponding native diagnostic profile.
    /// rustred_product_cache_bytes : Optional[int]
    ///     Separate lower-factor cache budget; not an aggregate with the saved
    ///     application cache. None imposes no payload-byte ceiling.
    /// rustred_product_cache_terms : Optional[int]
    ///     Optional lower-factor aggregate-term diagnostic. None leaves count
    ///     diagnostics disabled; a byte override alone does not enable them.
    /// rustred_cache_policy : str
    ///     "evict-lru" (default) evicts inactive completed saved-rule results;
    ///     "fail-on-limit" preserves them and reports a cache-limit error.
    ///     A protected working set that cannot fit still errors under LRU.
    /// rustred_scalar_max_work : Optional[int]
    ///     Optional scalar-work diagnostic per five-loop caller-numerator lowering.
    ///     None disables count diagnostics. This does not change campaign domains.
    /// max_process_rss_bytes : Optional[int]
    ///     Cooperative whole-process resident-memory ceiling for each request.
    ///     None imposes no RSS ceiling; cache payload budgets remain separate.
    /// max_elapsed_seconds : Optional[float]
    ///     Cooperative elapsed-time ceiling per public operation, including data
    ///     loading and lookup. None disables it. Native kernels are checked at
    ///     boundaries, so an external supervisor is needed for a hard deadline.
    pub fn new(
        run_time_decimal_precision: Option<u32>,
        evaluation_order: Option<Vec<PyRef<VakintEvaluationMethodWrapper>>>,
        tensor_reduction_method: Option<String>,
        epsilon_symbol: Option<Symbol>,
        mu_r_sq_symbol: Option<Symbol>,
        form_exe_path: Option<String>,
        python_exe_path: Option<String>,
        verify_numerator_identification: Option<bool>,
        integral_normalization_factor: Option<String>,
        allow_unknown_integrals: Option<bool>,
        clean_tmp_dir: Option<bool>,
        number_of_terms_in_epsilon_expansion: Option<i64>,
        use_dot_product_notation: Option<bool>,
        temporary_directory: Option<String>,
        rustred_artifact_path: Option<String>,
        rustred_cache_bytes: Option<usize>,
        rustred_cache_terms: Option<usize>,
        rustred_product_cache_bytes: Option<usize>,
        rustred_product_cache_terms: Option<usize>,
        rustred_cache_policy: &str,
        rustred_scalar_max_work: Option<usize>,
        max_process_rss_bytes: Option<u64>,
        max_elapsed_seconds: Option<f64>,
    ) -> PyResult<VakintWrapper> {
        let rustred_resources = resources::limits(max_process_rss_bytes, max_elapsed_seconds)?;
        resources::run(rustred_resources, || {
            let mut rustred_application_options = rustred_app::SavedApplicationOptions {
                cache_policy: self::rustred_cache_policy(rustred_cache_policy)
                    .map_err(exceptions::PyValueError::new_err)?,
                ..Default::default()
            };
            let mut rustred_scalar_numerator_limits =
                crate::RustRedScalarNumeratorLimits::default();
            if let Some(work) = rustred_scalar_max_work {
                rustred_scalar_numerator_limits =
                    rustred_scalar_numerator_limits.with_diagnostics();
                rustred_scalar_numerator_limits
                    .exact_algebra
                    .max_term_operations = work;
            }
            if let Some(bytes) = rustred_cache_bytes {
                rustred_application_options
                    .reduction_limits
                    .max_cached_coefficient_bytes = bytes;
            }
            if let Some(terms) = rustred_cache_terms {
                rustred_application_options.reduction_limits = rustred_application_options
                    .reduction_limits
                    .with_diagnostics();
                rustred_application_options
                    .reduction_limits
                    .max_cached_coefficient_terms = terms;
            }
            let mut rustred_product_limits =
                rustred::reduction::terminal_normalization::VacuumProductReductionLimits::default();
            if let Some(bytes) = rustred_product_cache_bytes {
                rustred_product_limits
                    .reduction
                    .max_cached_coefficient_bytes = bytes;
            }
            if let Some(terms) = rustred_product_cache_terms {
                rustred_product_limits = rustred_product_limits.with_diagnostics();
                rustred_product_limits
                    .reduction
                    .max_cached_coefficient_terms = terms;
            }
            let norm_factor = if let Some(nf) = integral_normalization_factor {
                match nf.as_str() {
                    "MSbar" => LoopNormalizationFactor::MSbar,
                    "pySecDec" => LoopNormalizationFactor::pySecDec,
                    "FMFTandMATAD" => LoopNormalizationFactor::FMFTandMATAD,
                    s => LoopNormalizationFactor::Custom(s.into()),
                }
            } else {
                LoopNormalizationFactor::MSbar
            };
            // let eval_order = if let Some(eo) = evaluation_order {
            //     let mut eval_stack = vec![];
            //     for m in eo.iter() {
            //         let method = Python::with_gil(|py| m.extract::<VakintEvaluationMethodWrapper>(py))?;
            //         eval_stack.push(method.method);
            //     }
            //     EvaluationOrder(eval_stack)
            // } else {
            //     EvaluationOrder::all()
            // };
            let eval_order = if let Some(eo) = evaluation_order {
                EvaluationOrder(
                    eo.iter()
                        .map(|m| m.method.clone())
                        .collect::<Vec<EvaluationMethod>>(),
                )
            } else {
                EvaluationOrder::all()
            };
            let eps_symbol = if let Some(es) = epsilon_symbol {
                es.to_string()
            } else {
                "ε".into()
            };
            let mu_r_sq_sym = if let Some(ms) = mu_r_sq_symbol {
                ms.to_string()
            } else {
                "mursq".into()
            };
            #[allow(clippy::needless_update)]
            let settings = VakintSettings {
                run_time_decimal_precision: run_time_decimal_precision.unwrap_or(17),
                epsilon_symbol: eps_symbol,
                mu_r_sq_symbol: mu_r_sq_sym,
                form_exe_path: form_exe_path.unwrap_or("form".into()),
                python_exe_path: python_exe_path.unwrap_or("python3".into()),
                verify_numerator_identification: verify_numerator_identification.unwrap_or(true),
                integral_normalization_factor: norm_factor,
                allow_unknown_integrals: allow_unknown_integrals.unwrap_or(true),
                clean_tmp_dir: clean_tmp_dir
                    .unwrap_or(env::var("VAKINT_NO_CLEAN_TMP_DIR").is_err()),
                temporary_directory,
                evaluation_order: eval_order,
                tensor_reduction_method: tensor_reduction_method
                    .as_deref()
                    .map(str::parse::<TensorReductionMethod>)
                    .transpose()
                    .map_err(vakint_to_python_error)?
                    .unwrap_or_default(),
                // This quantity is typically set equal to *one plus the maximum loop count* of the UV regularisation problem considered.
                // For example when considering a 2-loop problem, then:
                //   a) for the nested one-loop integrals appearing, the single pole, finite term *and* order-epsilon term will need to be considered.
                //   b) for the two-loop integrals, the double pole, single pole and finite terms will be needed, so again three terms
                number_of_terms_in_epsilon_expansion: number_of_terms_in_epsilon_expansion
                    .unwrap_or(4),
                use_dot_product_notation: use_dot_product_notation.unwrap_or(false),
                rustred_artifact_path: rustred_artifact_path.map(Into::into),
                rustred_application_options,
                rustred_scalar_numerator_limits,
                rustred_product_limits,
                rustred_resources,
                ..VakintSettings::default()
            };
            let vakint = Vakint::new().map_err(vakint_to_python_error)?;
            let wrapper = VakintWrapper { vakint, settings };
            Ok(wrapper)
        })
    }

    /// Explicitly load a local five-loop application package once.
    /// Also selects this package in the engine settings.
    pub fn load_rustred_application(&mut self, path: String) -> PyResult<()> {
        let (staged, settings) = resources::run(self.settings.rustred_resources, || {
            let staged = self.staged_engine()?;
            let mut settings = self.settings.clone();
            let staged = staged
                .with_rustred_application(&mut settings, path)
                .map_err(vakint_to_python_error)?;
            Ok((staged, settings))
        })?;
        self.vakint = staged;
        self.settings = settings;
        Ok(())
    }

    /// Saved-rule cache occupancy and cumulative work counters. Byte counts
    /// describe charged coefficient payloads, not total process RSS.
    /// `recomputed_nodes` counts repeated cache misses within a single request;
    /// no unbounded cross-request history is retained.
    pub fn rustred_application_cache_statistics(&self) -> PyResult<HashMap<String, usize>> {
        resources::run(self.settings.rustred_resources, || {
            let cache = self
                .vakint
                .rustred_application_cache_statistics()
                .map_err(vakint_to_python_error)?;
            let work = self
                .vakint
                .rustred_application_cache_work_statistics()
                .map_err(vakint_to_python_error)?;
            Ok([
                ("cached_integrals", cache.cached_integrals()),
                ("cached_coefficient_terms", cache.cached_coefficient_terms()),
                ("cached_coefficient_bytes", cache.cached_coefficient_bytes()),
                ("rule_applications", cache.rule_applications()),
                ("cache_hits", cache.cache_hits()),
                ("evicted_entries", work.evicted_entries),
                ("evicted_coefficient_terms", work.evicted_coefficient_terms),
                ("evicted_coefficient_bytes", work.evicted_coefficient_bytes),
                ("recomputed_nodes", work.recomputed_nodes),
            ]
            .into_iter()
            .map(|(key, value)| (key.to_owned(), value))
            .collect())
        })
    }

    /// Clear saved-rule and lower-factor memoization without changing the
    /// admitted scope. Cumulative work/eviction counters are not reset.
    pub fn clear_rustred_application_cache(&self) -> PyResult<()> {
        resources::run_committing(self.settings.rustred_resources, || {
            self.vakint
                .clear_rustred_application_cache()
                .map_err(vakint_to_python_error)
        })
    }

    /// Explicitly prepare a finite native cache for one through five loops.
    /// No cache attaches automatically and no numerical masters are computed.
    /// Inspect the returned status: resource/time/cancellation stops are partial.
    /// `workload` is a sequence of `(namespace, complete_integer_powers)` pairs;
    /// otherwise the native all-pinch relative-dot/rank primer is used for 1L–4L.
    /// Five-loop preparation requires an explicit selected workload and writes
    /// one compressed cache file. `checkpoint_path` selects its resumable work
    /// directory; it is only valid for 5L and never changes the final file format.
    ///
    /// The GIL is released during native work. `progress(report)` runs between
    /// roots (throttled) and may return False to cancel; `cancel()` is a polled
    /// predicate. Python signals and callback errors stop cooperatively after
    /// the current native operation, checkpoint completed roots, then re-raise
    /// with `cache_report` attached when possible. No hard deadline is promised.
    /// Count overrides (`cache_terms`, `cache_integrals`, `max_rule_applications`
    /// and `max_roots`) explicitly enable their native diagnostic profiles;
    /// omission leaves count stops disabled. `batch_records` selects the flush
    /// layout, while `extra_dots` and `numerator_rank` select a workload.
    /// `monitor` selects optional per-call native diagnostics; detaching its TUI
    /// never requests cancellation. `progress_snapshot` writes latest native JSON.
    #[pyo3(signature = (n_loops, output_path, *, ram_bytes=None, disk_bytes=None, cache_policy="evict-lru", cache_terms=None, cache_integrals=None, max_rule_applications=None, time_limit_seconds=None, checkpoint_interval_seconds=None, max_roots=None, batch_records=None, extra_dots=None, numerator_rank=None, workload=None, checkpoint_path=None, resume=true, progress=None, cancel=None, progress_interval_seconds=1.0, monitor="off", progress_snapshot=None))]
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn prepare_cache(
        &self,
        py: Python<'_>,
        n_loops: usize,
        output_path: String,
        ram_bytes: Option<u64>,
        disk_bytes: Option<u64>,
        cache_policy: &str,
        cache_terms: Option<usize>,
        cache_integrals: Option<usize>,
        max_rule_applications: Option<usize>,
        time_limit_seconds: Option<u64>,
        checkpoint_interval_seconds: Option<u64>,
        max_roots: Option<usize>,
        batch_records: Option<usize>,
        extra_dots: Option<u32>,
        numerator_rank: Option<u32>,
        workload: Option<Vec<(usize, Vec<i64>)>>,
        checkpoint_path: Option<PathBuf>,
        resume: bool,
        progress: Option<Py<PyAny>>,
        cancel: Option<Py<PyAny>>,
        progress_interval_seconds: f64,
        monitor: &str,
        progress_snapshot: Option<PathBuf>,
    ) -> PyResult<Py<PyDict>> {
        resources::run_committing(self.settings.rustred_resources, || {
            cache::check_loops(n_loops)?;
            let monitor =
                cache::monitor_options(monitor, progress_snapshot, progress_interval_seconds)?;
            let mut options = crate::RustRedCachePreparationOptions::default();
            options.storage.ram_bytes = ram_bytes;
            options.storage.disk_bytes = disk_bytes;
            options.storage.cache_policy =
                rustred_cache_policy(cache_policy).map_err(exceptions::PyValueError::new_err)?;
            options.storage.resume = resume;
            options.checkpoint_path = checkpoint_path;
            if let Some(value) = time_limit_seconds {
                options.storage.time_limit_seconds = value;
            }
            if let Some(value) = checkpoint_interval_seconds {
                options.storage.checkpoint_interval_seconds = value;
            }
            if let Some(value) = max_roots {
                options.storage = options.storage.with_diagnostics();
                options.storage.max_roots = value;
            }
            if let Some(value) = batch_records {
                options.storage.batch_records = value;
            }
            cache::apply_native_overrides(
                &mut options.native,
                None,
                cache_terms,
                cache_integrals,
                max_rule_applications,
            )?;
            options.workload = cache::workload(extra_dots, numerator_rank, workload)?;
            cache::prepare(
                py,
                self.vakint.clone(),
                n_loops,
                output_path.into(),
                options,
                progress,
                cancel,
                monitor,
            )
        })
    }

    /// Attach an exact read-only native cache bundle. Omitted settings preserve
    /// the active RAM, work, eviction and lookup policies. A new attachment
    /// uses native defaults, not preparation's larger RAM default.
    #[pyo3(signature = (n_loops, path, *, ram_bytes=None, disk_bytes=None, cache_policy=None, lookup_budget_policy=None, cache_terms=None, cache_integrals=None, max_rule_applications=None, file_options=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn load_cache(
        &mut self,
        py: Python<'_>,
        n_loops: usize,
        path: String,
        ram_bytes: Option<u64>,
        disk_bytes: Option<u64>,
        cache_policy: Option<&str>,
        lookup_budget_policy: Option<&str>,
        cache_terms: Option<usize>,
        cache_integrals: Option<usize>,
        max_rule_applications: Option<usize>,
        file_options: Option<PyRef<'_, cache::CacheFileOptionsWrapper>>,
    ) -> PyResult<()> {
        resources::run_committing(self.settings.rustred_resources, || {
            cache::check_loops(n_loops)?;
            let mut options = crate::RustRedCacheLoadOptions {
                file_options: file_options.map(|options| options.native.clone()),
                ..Default::default()
            };
            if ram_bytes.is_some()
                || cache_terms.is_some()
                || cache_integrals.is_some()
                || max_rule_applications.is_some()
            {
                let mut native = self
                    .vakint
                    .configured_cache_limits(n_loops)
                    .map_err(vakint_to_python_error)?
                    .unwrap_or_default();
                cache::apply_native_overrides(
                    &mut native,
                    ram_bytes,
                    cache_terms,
                    cache_integrals,
                    max_rule_applications,
                )?;
                options.native = Some(native);
            }
            options.policy = cache_policy
                .map(rustred_cache_policy)
                .transpose()
                .map_err(exceptions::PyValueError::new_err)?;
            options.lookup_budget_policy =
                lookup_budget_policy.map(cache::lookup_policy).transpose()?;
            if let Some(value) = disk_bytes {
                options.storage.max_disk_bytes = value;
            }
            let engine = self.vakint.clone();
            py.detach(move || engine.load_cache(n_loops, path, &options))
                .map_err(vakint_to_python_error)
        })
    }

    /// Attach incremental native capture, creating its delta if necessary.
    /// Immutable bases remain separate. `read_backing=False` disables all
    /// persistent reads while retaining ordinary RAM memoization and capture.
    /// `read_delta=False` reads only the pinned bases and snapshots. The returned
    /// control can report or checkpoint while evaluation holds the engine.
    /// Queue and checkpoint options use the native defaults when omitted.
    #[pyo3(signature = (n_loops, path, *, read_backing=true, read_delta=true, base_paths=None, snapshots=None, disk_bytes=None, max_pending_records=None, max_pending_bytes=None, batch_records=None, batch_bytes=None, checkpoint_interval_seconds=None, file_options=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn attach_running_cache(
        &mut self,
        py: Python<'_>,
        n_loops: usize,
        path: PathBuf,
        read_backing: bool,
        read_delta: bool,
        base_paths: Option<Vec<PathBuf>>,
        snapshots: Option<Vec<Py<PyDict>>>,
        disk_bytes: Option<u64>,
        max_pending_records: Option<usize>,
        max_pending_bytes: Option<usize>,
        batch_records: Option<usize>,
        batch_bytes: Option<usize>,
        checkpoint_interval_seconds: Option<f64>,
        file_options: Option<PyRef<'_, cache::CacheFileOptionsWrapper>>,
    ) -> PyResult<cache::RunningCacheControlWrapper> {
        resources::run_committing(self.settings.rustred_resources, || {
            cache::check_loops(n_loops)?;
            let mut options = cache::running_options(
                py,
                read_backing,
                read_delta,
                base_paths,
                snapshots,
                disk_bytes,
                max_pending_records,
                max_pending_bytes,
                batch_records,
                batch_bytes,
                checkpoint_interval_seconds,
            )?;
            options.file_options = file_options.map(|options| options.native.clone());
            let engine = self.vakint.clone();
            let control = py
                .detach(move || engine.attach_running_cache(n_loops, path, &options))
                .map_err(vakint_to_python_error)?;
            Ok(cache::RunningCacheControlWrapper {
                control,
                resources: self.settings.rustred_resources,
            })
        })
    }

    /// Resume an existing incremental native cache; never create a new delta.
    /// This has the same read policy and queue options as `attach_running_cache`.
    #[pyo3(signature = (n_loops, path, *, read_backing=true, read_delta=true, base_paths=None, snapshots=None, disk_bytes=None, max_pending_records=None, max_pending_bytes=None, batch_records=None, batch_bytes=None, checkpoint_interval_seconds=None, file_options=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn load_running_cache(
        &mut self,
        py: Python<'_>,
        n_loops: usize,
        path: PathBuf,
        read_backing: bool,
        read_delta: bool,
        base_paths: Option<Vec<PathBuf>>,
        snapshots: Option<Vec<Py<PyDict>>>,
        disk_bytes: Option<u64>,
        max_pending_records: Option<usize>,
        max_pending_bytes: Option<usize>,
        batch_records: Option<usize>,
        batch_bytes: Option<usize>,
        checkpoint_interval_seconds: Option<f64>,
        file_options: Option<PyRef<'_, cache::CacheFileOptionsWrapper>>,
    ) -> PyResult<cache::RunningCacheControlWrapper> {
        resources::run_committing(self.settings.rustred_resources, || {
            cache::check_loops(n_loops)?;
            let mut options = cache::running_options(
                py,
                read_backing,
                read_delta,
                base_paths,
                snapshots,
                disk_bytes,
                max_pending_records,
                max_pending_bytes,
                batch_records,
                batch_bytes,
                checkpoint_interval_seconds,
            )?;
            options.file_options = file_options.map(|options| options.native.clone());
            let engine = self.vakint.clone();
            let control = py
                .detach(move || engine.load_running_cache(n_loops, path, &options))
                .map_err(vakint_to_python_error)?;
            Ok(cache::RunningCacheControlWrapper {
                control,
                resources: self.settings.rustred_resources,
            })
        })
    }

    /// Attach immutable bases and frozen native descriptors without a writer.
    /// Frozen generations do not refresh, and evaluation captures no records.
    #[pyo3(signature = (n_loops, *, base_paths=None, snapshots=None, disk_bytes=None, file_options=None))]
    pub fn load_cache_snapshot(
        &mut self,
        py: Python<'_>,
        n_loops: usize,
        base_paths: Option<Vec<PathBuf>>,
        snapshots: Option<Vec<Py<PyDict>>>,
        disk_bytes: Option<u64>,
        file_options: Option<PyRef<'_, cache::CacheFileOptionsWrapper>>,
    ) -> PyResult<()> {
        resources::run_committing(self.settings.rustred_resources, || {
            cache::check_loops(n_loops)?;
            let mut options = cache::running_options(
                py, true, false, base_paths, snapshots, disk_bytes, None, None, None, None, None,
            )?;
            options.file_options = file_options.map(|options| options.native.clone());
            let engine = self.vakint.clone();
            py.detach(move || engine.load_cache_snapshot(n_loops, &options))
                .map_err(vakint_to_python_error)
        })
    }

    /// Flush accepted native records and return the acknowledged durable report.
    /// A flush failure is an error; pending records must not be reported saved.
    pub fn save_running_cache(&self, py: Python<'_>, n_loops: usize) -> PyResult<Py<PyDict>> {
        resources::run_committing(self.settings.rustred_resources, || {
            cache::check_loops(n_loops)?;
            let engine = self.vakint.clone();
            let report = py
                .detach(move || engine.save_running_cache(n_loops))
                .map_err(vakint_to_python_error)?;
            cache::running_report(py, &report)
        })
    }

    /// Flush and freeze a small native descriptor for this delta generation.
    /// Immutable base payloads are referenced separately and are never copied.
    pub fn freeze_running_cache(
        &self,
        py: Python<'_>,
        n_loops: usize,
        sidecar_path: PathBuf,
    ) -> PyResult<Py<PyDict>> {
        resources::run_committing(self.settings.rustred_resources, || {
            cache::check_loops(n_loops)?;
            let engine = self.vakint.clone();
            let snapshot = py
                .detach(move || engine.freeze_running_cache(n_loops, sidecar_path))
                .map_err(vakint_to_python_error)?;
            cache::json_report(py, &snapshot.to_document().to_string())
        })
    }

    /// Explicitly save `resident-roots` or a portable `full-snapshot`.
    /// A full snapshot includes captured bases and completed resident nodes.
    /// Five loops stream one shared-table file; a running writer must already
    /// be closed. Lower loops retain native byte-copy/checkpoint reporting.
    /// Every full snapshot requires a new destination.
    #[pyo3(signature = (n_loops, output_path, *, mode, disk_bytes=None, resume=false, file_options=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn save_cache(
        &self,
        py: Python<'_>,
        n_loops: usize,
        output_path: String,
        mode: &str,
        disk_bytes: Option<u64>,
        resume: bool,
        file_options: Option<PyRef<'_, cache::CacheFileOptionsWrapper>>,
    ) -> PyResult<Py<PyDict>> {
        resources::run_committing(self.settings.rustred_resources, || {
            cache::check_loops(n_loops)?;
            let mode = cache::save_mode(mode)?;
            let mut options = crate::RustRedCacheSaveOptions::default();
            options.disk_bytes = disk_bytes;
            options.resume = resume;
            if let Some(file_options) = file_options {
                options = options.with_file_options(file_options.native.clone());
            }
            let engine = self.vakint.clone();
            let report = py
                .detach(move || engine.save_cache(n_loops, output_path, mode, &options))
                .map_err(vakint_to_python_error)?;
            cache::json_report(py, &report.to_document().to_string())
        })
    }

    /// Import a trusted SHA-256 checked archive or native multipart manifest.
    /// The imported snapshot is not attached to evaluation automatically.
    #[cfg(all(unix, not(target_arch = "wasm32")))]
    #[pyo3(signature = (n_loops, archive_path, output_path, *, sha256, format="archive", max_manifest_bytes=None, max_archive_bytes=None, max_extracted_bytes=None, max_entries=None, max_temporary_bytes=None, disk_bytes=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn import_cache(
        &self,
        py: Python<'_>,
        n_loops: usize,
        archive_path: String,
        output_path: String,
        sha256: &str,
        format: &str,
        max_manifest_bytes: Option<u64>,
        max_archive_bytes: Option<u64>,
        max_extracted_bytes: Option<u64>,
        max_entries: Option<usize>,
        max_temporary_bytes: Option<u64>,
        disk_bytes: Option<u64>,
    ) -> PyResult<Py<PyDict>> {
        resources::run_committing(self.settings.rustred_resources, || {
            cache::check_loops(n_loops)?;
            let options = cache::transfer_options(
                sha256,
                format,
                max_manifest_bytes,
                max_archive_bytes,
                max_extracted_bytes,
                max_entries,
                max_temporary_bytes,
                disk_bytes,
                None,
            )?;
            let engine = self.vakint.clone();
            let report = py
                .detach(move || engine.import_cache(n_loops, archive_path, output_path, &options))
                .map_err(vakint_to_python_error)?;
            cache::transfer_report(py, report, format)
        })
    }

    /// Explicit download with a required trusted checksum. Archive downloads
    /// refuse redirects; HTTPS multipart downloads follow at most eight redirects.
    /// Load the published snapshot separately.
    #[cfg(all(unix, not(target_arch = "wasm32")))]
    #[pyo3(signature = (n_loops, url, output_path, *, sha256, format="archive", max_manifest_bytes=None, max_archive_bytes=None, max_extracted_bytes=None, max_entries=None, max_temporary_bytes=None, disk_bytes=None, timeout_seconds=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn download_cache(
        &self,
        py: Python<'_>,
        n_loops: usize,
        url: String,
        output_path: String,
        sha256: &str,
        format: &str,
        max_manifest_bytes: Option<u64>,
        max_archive_bytes: Option<u64>,
        max_extracted_bytes: Option<u64>,
        max_entries: Option<usize>,
        max_temporary_bytes: Option<u64>,
        disk_bytes: Option<u64>,
        timeout_seconds: Option<f64>,
    ) -> PyResult<Py<PyDict>> {
        resources::run_committing(self.settings.rustred_resources, || {
            cache::check_loops(n_loops)?;
            let options = cache::transfer_options(
                sha256,
                format,
                max_manifest_bytes,
                max_archive_bytes,
                max_extracted_bytes,
                max_entries,
                max_temporary_bytes,
                disk_bytes,
                timeout_seconds,
            )?;
            let engine = self.vakint.clone();
            let report = py
                .detach(move || engine.download_cache(n_loops, &url, output_path, &options))
                .map_err(vakint_to_python_error)?;
            cache::transfer_report(py, report, format)
        })
    }

    /// Clear only the attached native RAM cache; preserve its read-only files
    /// and cumulative counters. No program is loaded for an empty attachment.
    pub fn clear_cache(&self, py: Python<'_>, n_loops: usize) -> PyResult<()> {
        resources::run_committing(self.settings.rustred_resources, || {
            cache::check_loops(n_loops)?;
            let engine = self.vakint.clone();
            py.detach(move || engine.clear_cache(n_loops))
                .map_err(vakint_to_python_error)
        })
    }

    /// Detach persistent backing and clear resident reduction records.
    /// The loaded reduction campaign remains available for uncached evaluation.
    pub fn detach_cache(&self, py: Python<'_>, n_loops: usize) -> PyResult<()> {
        resources::run_committing(self.settings.rustred_resources, || {
            cache::check_loops(n_loops)?;
            let engine = self.vakint.clone();
            py.detach(move || engine.detach_cache(n_loops))
                .map_err(vakint_to_python_error)
        })
    }

    /// Return aggregate and exact-namespace cache telemetry. Coefficient payload
    /// and I/O buffers are separate; neither is a measurement of process RSS.
    pub fn cache_statistics(&self, py: Python<'_>, n_loops: usize) -> PyResult<Py<PyDict>> {
        resources::run(self.settings.rustred_resources, || {
            cache::check_loops(n_loops)?;
            let engine = self.vakint.clone();
            let stats = py
                .detach(move || engine.cache_statistics(n_loops))
                .map_err(vakint_to_python_error)?;
            cache::statistics(py, &stats)
        })
    }

    /// Apply the selected five-loop package and retain its symbolic conditions.
    ///
    /// Returns `(value, nonzero_conditions, families)`, where `families` contains
    /// `(master_label, full_family_fingerprint)` pairs. Conditions are exact
    /// expressions required to be nonzero at generic dimension. No numerical
    /// master values are substituted. The ordinary native tensor prepass runs
    /// before scalar application, with the same settings as `evaluate`.
    /// Optional `monitor` diagnostics are local to this call. A snapshot path
    /// enables JSON output independently of the terminal display.
    #[pyo3(signature = (integral_expression, *, monitor="off", progress_snapshot=None, progress_interval_seconds=1.0))]
    #[allow(clippy::type_complexity)]
    pub fn evaluate_rustred_symbolic(
        &self,
        py: Python<'_>,
        integral_expression: PythonExpression,
        monitor: &str,
        progress_snapshot: Option<PathBuf>,
        progress_interval_seconds: f64,
    ) -> PyResult<(
        PythonExpression,
        Vec<PythonExpression>,
        Vec<(PythonExpression, String)>,
    )> {
        resources::run(self.settings.rustred_resources, || {
            let monitor =
                cache::monitor_options(monitor, progress_snapshot, progress_interval_seconds)?;
            let result = if monitor.enabled() {
                let engine = self.vakint.clone();
                let settings = self.settings.clone();
                let input = integral_expression.expr.clone();
                py.detach(move || {
                    engine.evaluate_rustred_symbolic_with_monitor(
                        &settings,
                        input.as_view(),
                        &monitor,
                    )
                })
            } else {
                self.vakint
                    .evaluate_rustred_symbolic(&self.settings, integral_expression.expr.as_view())
            }
            .map_err(vakint_to_python_error)?;
            Ok((
                result.value.into(),
                result
                    .nonzero_conditions
                    .into_iter()
                    .map(Into::into)
                    .collect(),
                result
                    .families
                    .into_iter()
                    .map(|(label, fingerprint)| (label.into(), fingerprint))
                    .collect(),
            ))
        })
    }

    /// Inspect a registered RustRed constant in the active data snapshot.
    ///
    /// Precision shortfalls are returned as structured warnings. No data is downloaded.
    #[pyo3(signature = (symbol, prec=16))]
    pub fn get_constant(
        &self,
        symbol: &ConstantSymbolWrapper,
        prec: u32,
    ) -> PyResult<VakintConstantWrapper> {
        resources::run(self.settings.rustred_resources, || {
            self.vakint
                .get_constant(&symbol.value, prec)
                .map(Into::into)
                .map_err(vakint_to_python_error)
        })
    }

    /// Inspect every available Laurent coefficient and the authenticated master graph.
    ///
    /// The returned object keeps its data snapshot even after this engine loads another release.
    #[pyo3(signature = (symbol, prec=16))]
    pub fn get_master(
        &self,
        symbol: &MasterSymbolWrapper,
        prec: u32,
    ) -> PyResult<VakintMasterWrapper> {
        resources::run(self.settings.rustred_resources, || {
            self.vakint
                .get_master(&symbol.value, prec)
                .map(Into::into)
                .map_err(vakint_to_python_error)
        })
    }

    /// Discover all registered constants in the loaded or vendored RustRed data.
    pub fn constant_symbols(&self) -> PyResult<Vec<ConstantSymbolWrapper>> {
        resources::run(self.settings.rustred_resources, || {
            self.vakint
                .constant_symbols()
                .map(|symbols| symbols.into_iter().map(Into::into).collect())
                .map_err(vakint_to_python_error)
        })
    }

    /// Discover qualified master identities, including declared uncomputed terminals.
    pub fn master_symbols(&self) -> PyResult<Vec<MasterSymbolWrapper>> {
        resources::run(self.settings.rustred_resources, || {
            self.vakint
                .master_symbols()
                .map(|symbols| symbols.into_iter().map(Into::into).collect())
                .map_err(vakint_to_python_error)
        })
    }

    #[pyo3(signature = (expr))]
    /// Interpret a Symbolica expression as a numerical Laurent series in epsilon.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.hepkit.vakint import Vakint
    /// >>> vakint = Vakint(evaluation_order=[])
    /// >>> result = vakint.numerical_result_from_expression(
    /// ...     E("vakint::ε^-2 + 1 + 0.12*vakint::ε^-1")
    /// ... )
    /// >>> sorted(exponent for exponent, _ in result.to_list())
    /// [-2, -1, 0]
    /// ```
    ///
    /// Parameters
    /// ----------
    ///
    /// expr : Expression
    ///   A Symbolica expression representing a Laurent series in the dimensional regularisation parameter epsilon specified in the vakint engine.
    pub fn numerical_result_from_expression(
        &self,
        expr: PythonExpression,
    ) -> PyResult<NumericalEvaluationResultWrapper> {
        resources::run(self.settings.rustred_resources, || {
            let value = NumericalEvaluationResult::from_atom(
                expr.expr.as_view(),
                vakint_symbol!(self.settings.epsilon_symbol.clone()),
                &self.settings,
            )
            .map_err(vakint_to_python_error)?;
            Ok(NumericalEvaluationResultWrapper { value })
        })
    }

    #[pyo3(signature = (evaluated_integral, params, externals = None))]
    /// Substitute numerical parameters into an integral already evaluated parametrically by Vakint.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.hepkit.vakint import Vakint
    /// >>> vakint = Vakint(evaluation_order=[])
    /// >>> evaluated = E(
    /// ...     "muvsq*vakint::ε^-1 + mursq",
    /// ...     default_namespace="vakint",
    /// ... )
    /// >>> result, error = vakint.numerical_evaluation(
    /// ...     evaluated,
    /// ...     {"muvsq": 2.0, "mursq": 3.0},
    /// ... )
    /// >>> sorted(exponent for exponent, _ in result.to_list())
    /// [-1, 0]
    /// >>> error is None
    /// True
    /// ```
    ///
    /// Parameters
    /// ----------
    ///
    /// evaluated_integral : Expression
    ///   A Symbolica expression representing an integral that has been evaluated parametrically by Vakint.
    /// params : Mapping[str, Float | ComplexFloat | Expression | int | float | complex]
    ///   Native floats retain precision; Python floats/complex values carry 53-bit input precision.
    ///   Integers and exact expressions are substituted symbolically; inexact expression leaves
    ///   retain named source errors instead of becoming exact values.
    /// externals : Optional[Mapping[int, Sequence[Float | ComplexFloat | Expression | int | float | complex]]]
    ///   Four components in (time, x, y, z) order. Contractions use tracked arithmetic.
    pub fn numerical_evaluation(
        &self,
        evaluated_integral: PythonExpression,
        #[gen_stub(override_type(type_repr="typing.Mapping[str, symbolica.Float | symbolica.ComplexFloat | symbolica.Expression | int | float | complex]", imports=("typing", "symbolica")))]
        params: &Bound<'_, PyMapping>,
        #[gen_stub(override_type(type_repr="typing.Mapping[int, typing.Sequence[symbolica.Float | symbolica.ComplexFloat | symbolica.Expression | int | float | complex]] | None", imports=("typing", "symbolica")))]
        externals: Option<&Bound<'_, PyMapping>>,
    ) -> PyResult<(
        NumericalEvaluationResultWrapper,
        Option<NumericalEvaluationResultWrapper>,
    )> {
        resources::run(self.settings.rustred_resources, || {
            let (numerical_result, numerical_error) = numerical::evaluate(
                &self.vakint,
                &self.settings,
                &evaluated_integral.expr,
                params,
                externals,
            )?;
            Ok((
                NumericalEvaluationResultWrapper {
                    value: numerical_result,
                },
                numerical_error.map(|value| NumericalEvaluationResultWrapper { value }),
            ))
        })
    }

    #[pyo3(signature = (result))]
    /// Convert a Vakint numerical result to a Symbolica Laurent-series expression.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.hepkit.vakint import Vakint, VakintNumericalResult
    /// >>> vakint = Vakint(evaluation_order=[])
    /// >>> result = VakintNumericalResult([
    /// ...     (-1, (2.0, 0.0)),
    /// ...     (0, (3.0, 0.0)),
    /// ... ])
    /// >>> expression = vakint.numerical_result_to_expression(result)
    /// >>> "ε" in str(expression)
    /// True
    /// ```
    pub fn numerical_result_to_expression(
        &self,
        result: PyRef<NumericalEvaluationResultWrapper>,
    ) -> PyResult<PythonExpression> {
        resources::run(self.settings.rustred_resources, || {
            let value = result
                .value
                .to_atom(vakint_symbol!(self.settings.epsilon_symbol.clone()))
                .map_err(vakint_to_python_error)?;

            Ok(value.into())
        })
    }

    #[pyo3(signature = (integral_expression, short_form = None))]
    /// Convert a Vakint expression to canonical momentum routing and topology numbering.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.hepkit.vakint import Vakint
    /// >>> vakint = Vakint(evaluation_order=[])
    /// >>> integral = E(
    /// ...     "topo(prop(18,edge(7,7),k(99),muvsq,1))",
    /// ...     default_namespace="vakint",
    /// ... )
    /// >>> canonical = vakint.to_canonical(integral, short_form=True)
    /// >>> "I1L" in str(canonical)
    /// True
    /// ```
    ///
    /// Parameters
    /// ----------
    ///
    /// integral_expression : Expression
    ///   A Symbolica expression representing a vakint integral.
    /// short_form : Optional[bool]
    ///   Whether to use the short form for the topology representation. Default is False.
    pub fn to_canonical(
        &self,
        integral_expression: PythonExpression,
        short_form: Option<bool>,
    ) -> PyResult<PythonExpression> {
        resources::run(self.settings.rustred_resources, || {
            let result = self
                .vakint
                .to_canonical(
                    &self.settings,
                    integral_expression.expr.as_view(),
                    short_form.unwrap_or(false),
                )
                .map_err(vakint_to_python_error)?;
            Ok(result.into())
        })
    }

    #[pyo3(signature = (integral_expression))]
    /// Convert a vakint expression to a form where tensor integrals are reduced to scalar integrals.
    /// The backend is selected by `tensor_reduction_method` when constructing `Vakint`.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.hepkit.vakint import Vakint
    /// >>> vakint = Vakint(evaluation_order=[])
    /// >>> integral = E(
    /// ...     "k(1,101)*k(1,102)*topo(prop(1,edge(1,1),k(1),muvsq,1))",
    /// ...     default_namespace="vakint",
    /// ... )
    /// >>> reduced = vakint.tensor_reduce(integral)
    /// >>> "g(101,102)" in str(reduced)
    /// True
    /// ```
    ///
    /// Parameters
    /// ----------
    /// integral_expression : Expression
    ///    A Symbolica expression representing a vakint integral.
    pub fn tensor_reduce(
        &self,
        integral_expression: PythonExpression,
    ) -> PyResult<PythonExpression> {
        resources::run(self.settings.rustred_resources, || {
            let result = self
                .vakint
                .tensor_reduce(&self.settings, integral_expression.expr.as_view())
                .map_err(vakint_to_python_error)?;
            Ok(result.into())
        })
    }

    #[pyo3(signature = (integral_expression))]
    /// Perform the parametric evaluation of *only the integral* appearing in the Symbolica expression given in input representing a vakint integral.
    /// The numerator is left unchanged.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.hepkit.vakint import Vakint, VakintEvaluationMethod
    /// >>> vakint = Vakint(
    /// ...     evaluation_order=[VakintEvaluationMethod.new_alphaloop_method()]
    /// ... )
    /// >>> integral = E(
    /// ...     "topo(prop(1,edge(1,1),k(1),muvsq,1))",
    /// ...     default_namespace="vakint",
    /// ... )
    /// >>> evaluated = vakint.evaluate_integral(integral)
    /// >>> "ε" in str(evaluated)
    /// True
    /// ```
    ///
    /// The AlphaLoop evaluation method used here invokes FORM. Configure `form_exe_path` if
    /// FORM is not available as `form` on `PATH`.
    ///
    /// Parameters
    /// ----------
    ///
    /// integral_expression : Expression
    ///   A Symbolica expression representing a vakint integral.
    pub fn evaluate_integral(
        &self,
        integral_expression: PythonExpression,
    ) -> PyResult<PythonExpression> {
        resources::run(self.settings.rustred_resources, || {
            let result = self
                .vakint
                .evaluate_integral(&self.settings, integral_expression.expr.as_view())
                .map_err(vakint_to_python_error)?;
            Ok(result.into())
        })
    }

    #[pyo3(signature = (integral_expression, *, monitor="off", progress_snapshot=None, progress_interval_seconds=1.0))]
    /// Perform the complete parametric evaluation of the Vakint integral represented by the Symbolica expression given in input.
    /// Note that the tensor reduction will be automatically performed on the input given.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.hepkit.vakint import Vakint, VakintEvaluationMethod
    /// >>> vakint = Vakint(
    /// ...     evaluation_order=[VakintEvaluationMethod.new_alphaloop_method()]
    /// ... )
    /// >>> integral = E(
    /// ...     "k(1,101)*k(1,102)*topo(prop(1,edge(1,1),k(1),muvsq,1))",
    /// ...     default_namespace="vakint",
    /// ... )
    /// >>> evaluated = vakint.evaluate(integral)
    /// >>> "g(101,102)" in str(evaluated)
    /// True
    /// ```
    ///
    /// This complete path performs tensor reduction before integral evaluation and therefore
    /// has the same FORM requirement as `evaluate_integral` for the AlphaLoop method.
    ///
    /// Parameters
    /// ----------
    ///
    /// integral_expression : Expression
    ///   A Symbolica expression representing a vakint integral.
    /// Optional native diagnostics use `monitor="tui"` or `monitor="json"`.
    /// A `progress_snapshot` path also enables JSON with terminal monitoring off.
    pub fn evaluate(
        &self,
        py: Python<'_>,
        integral_expression: PythonExpression,
        monitor: &str,
        progress_snapshot: Option<PathBuf>,
        progress_interval_seconds: f64,
    ) -> PyResult<PythonExpression> {
        resources::run(self.settings.rustred_resources, || {
            let monitor =
                cache::monitor_options(monitor, progress_snapshot, progress_interval_seconds)?;
            let result = if monitor.enabled() {
                let engine = self.vakint.clone();
                let settings = self.settings.clone();
                let input = integral_expression.expr.clone();
                py.detach(move || {
                    engine.evaluate_with_monitor(&settings, input.as_view(), &monitor)
                })
            } else {
                self.vakint
                    .evaluate(&self.settings, integral_expression.expr.as_view())
            }
            .map_err(vakint_to_python_error)?;
            Ok(result.into())
        })
    }
}

fn record_usage() {
    crate::citations::USED_CITATIONS.record(crate::citations::CitationSource::Vakint);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modern_module_identity_preserves_default_backend_order() {
        assert_eq!(VakintWrapper::get_name(), "hepkit_vakint");
        let defaults = VakintSettings::default();
        assert!(matches!(
            defaults.evaluation_order.0.as_slice(),
            [
                EvaluationMethod::AlphaLoop(_),
                EvaluationMethod::MATAD(_),
                EvaluationMethod::FMFT(_),
                EvaluationMethod::PySecDec(_)
            ]
        ));
        assert!(matches!(
            defaults.tensor_reduction_method,
            TensorReductionMethod::FeynKit
        ));
    }

    #[test]
    fn native_payloads_preserve_numerical_precision_without_python() {
        // Python rendering/extraction is tested in an installed-host subprocess,
        // never by importing a second Symbolica kernel into this Rust binary.
        let value = Float::parse("1.2345678901234567890123456789012345", Some(160)).unwrap();
        let coefficient = Complex::new(value.clone(), Float::with_val(160, 0));
        let result = NumericalEvaluationResult::from_values(vec![(-2, coefficient.clone())], 32);
        let epsilon = vakint_symbol!("ε");
        // Checked storage may round to the requested 32 digits, but must not
        // pass through binary64 on the way to a Python numerical expression.
        assert!(result.coefficients()[&-2].value().re.prec() >= 100);
        let rounded = NumericalEvaluationResult::from_values(
            vec![(
                -2,
                Complex::new(
                    Float::with_val(160, value.to_f64()),
                    Float::with_val(160, 0),
                ),
            )],
            32,
        );
        assert_ne!(
            result.to_atom(epsilon).unwrap(),
            rounded.to_atom(epsilon).unwrap()
        );

        let topology = symbolica::function!(
            S.topo,
            symbolica::function!(
                S.prop,
                1,
                symbolica::function!(S.edge, 1, 1),
                symbolica::function!(S.k, 1),
                Atom::num(1),
                1
            )
        );
        let split = VakintExpression::try_from(topology.clone()).unwrap();
        assert_eq!(Atom::from(split.clone()), topology);
        assert_eq!(Atom::from(split), topology);
    }

    #[test]
    fn python_rustred_method_constructor_preserves_its_option() {
        Python::initialize();
        Python::attach(|py| {
            let class = py.get_type::<VakintEvaluationMethodWrapper>();

            let default_method =
                VakintEvaluationMethodWrapper::new_rustred_method(&class, "numerical").unwrap();
            let EvaluationMethod::RustRed(default_options) = default_method.method else {
                panic!("constructor returned a non-RustRed evaluation method");
            };
            assert!(default_options.master_mode == crate::RustRedMasterMode::Numerical);

            let unsubstituted_method =
                VakintEvaluationMethodWrapper::new_rustred_method(&class, "retain").unwrap();
            let EvaluationMethod::RustRed(unsubstituted_options) = unsubstituted_method.method
            else {
                panic!("constructor returned a non-RustRed evaluation method");
            };
            assert!(unsubstituted_options.master_mode == crate::RustRedMasterMode::Retain);
        });
    }

    #[test]
    fn python_saved_cache_policy_values_are_explicit_and_validated() {
        use rustred::solver::CandidateRoutedCachePolicy;
        assert_eq!(
            rustred_cache_policy("evict-lru").unwrap(),
            CandidateRoutedCachePolicy::EvictLeastRecentlyUsed
        );
        assert_eq!(
            rustred_cache_policy("fail-on-limit").unwrap(),
            CandidateRoutedCachePolicy::FailOnLimit
        );
        assert_eq!(
            rustred_cache_policy("evict-aged-lfu").unwrap(),
            CandidateRoutedCachePolicy::EvictAgedLeastFrequentlyUsed
        );
        assert!(rustred_cache_policy("unbounded").is_err());
    }
}
