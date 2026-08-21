use symbolica::api::python::SymbolicaCommunityModule;

use std::collections::HashMap;
use std::env;

use pyo3::types::{PyModule, PyModuleMethods, PyType};
use pyo3::{
    Bound, FromPyObject, PyAny, PyClass, PyErr, PyRef, Python, exceptions, pyclass, pymethods,
};
use pyo3::{Py, PyResult};
use symbolica::api::python::PythonExpression;
use symbolica::atom::{Atom, Symbol};
use symbolica::domains::float::{Complex, Float, RealLike};

use crate::symbols::S;
use crate::{
    AlphaLoopOptions, EvaluationMethod, EvaluationOrder, FMFTOptions, LoopNormalizationFactor,
    MATADOptions, NumericalEvaluationResult, PySecDecOptions, Vakint, VakintError,
    VakintExpression, VakintSettings, vakint_symbol,
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::*;

impl SymbolicaCommunityModule for VakintWrapper {
    fn get_name() -> String {
        "vakint".to_string()
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
            Ok(())
        }

        /// The classes registered on `symbolica.community.vakint`.
        #[cfg(feature = "python_stubgen")]
        pub const PYTHON_STUB_SURFACE: &[&str] = &[$(<$class as PyClass>::NAME,)+];
    };
}

define_vakint_python_surface! {
    VakintWrapper,
    NumericalEvaluationResultWrapper,
    VakintExpressionWrapper,
    VakintEvaluationMethodWrapper,
}

fn vakint_to_python_error(vakint_error: VakintError) -> PyErr {
    exceptions::PyValueError::new_err(format!("Vakint error | {vakint_error}"))
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "Vakint", module = "symbolica.community.vakint")]
/// Vakint engine and settings used for matching, reduction, and evaluation.
///
/// Construct one instance and reuse it: initialization processes the complete topology library.
pub struct VakintWrapper {
    pub vakint: Vakint,
    pub settings: VakintSettings,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "VakintNumericalResult", module = "symbolica.community.vakint")]
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
#[pymethods]
impl NumericalEvaluationResultWrapper {
    /// String representation of the numerical result.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.vakint import VakintNumericalResult
    /// >>> result = VakintNumericalResult([
    /// ...     (-3, (0.0, -11440.53140354612)),
    /// ...     (-2, (0.0, 57169.95521898031)),
    /// ...     (-1, (0.0, -178748.9838377694)),
    /// ...     (0, (0.0, 321554.1122184795)),
    /// ... ])
    /// >>> str(result)
    /// ε^-3 : (0+-11440.5314035461i)
    /// ε^-2 : (0+57169.9552189803i)
    /// ε^-1 : (0+-178748.983837769i)
    /// ε^ 0 : (0+321554.112218480i)
    /// ```
    pub fn __str__(&self) -> PyResult<String> {
        Ok(format!("{}", self.value))
    }

    /// Convert the numerical result to a native Python list of (epsilon exponent, (real, imag)) tuples.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.vakint import VakintNumericalResult
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
            .0
            .iter()
            .map(|(exp, coeff)| (*exp, (coeff.re.to_f64(), coeff.im.to_f64())))
            .collect())
    }
    #[new]
    /// Create a numerical Laurent series from `(epsilon exponent, (real, imaginary))` tuples.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.vakint import VakintNumericalResult
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
        let r = NumericalEvaluationResult(
            values
                .iter()
                .map(|(eps_pwr, (re, im))| {
                    (
                        *eps_pwr,
                        Complex::new(Float::with_val(53, re), Float::with_val(53, im)),
                    )
                })
                .collect::<Vec<_>>(),
        );
        Ok(NumericalEvaluationResultWrapper { value: r })
    }

    #[pyo3(signature = (other, relative_threshold, error = None, max_pull = None))]
    /// Compare this numerical result to another, returning a tuple of (bool, str) where the bool indicates whether the results match within the specified thresholds,
    /// and the str provides details of the comparison.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.vakint import VakintNumericalResult
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
#[pyclass(name = "VakintExpression", module = "symbolica.community.vakint")]
pub struct VakintExpressionWrapper {
    pub value: VakintExpression,
}

impl<'a, 'py> FromPyObject<'a, 'py> for VakintExpressionWrapper {
    type Error = PyErr;

    fn extract(ob: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(a) = ob.extract::<VakintExpressionWrapper>() {
            Ok(VakintExpressionWrapper { value: a.value })
        } else {
            Err(exceptions::PyValueError::new_err(
                "Not a valid vakint expression",
            ))
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl VakintExpressionWrapper {
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
    /// >>> from symbolica.community.vakint import VakintExpression
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
    /// >>> from symbolica.community.vakint import VakintExpression
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
#[pyclass(name = "VakintEvaluationMethod", module = "symbolica.community.vakint")]
/// One configured backend in a `Vakint` instance's evaluation order.
pub struct VakintEvaluationMethodWrapper {
    pub method: EvaluationMethod,
}

impl<'a, 'py> FromPyObject<'a, 'py> for VakintEvaluationMethodWrapper {
    type Error = PyErr;

    fn extract(ob: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(a) = ob.extract::<VakintEvaluationMethodWrapper>() {
            Ok(VakintEvaluationMethodWrapper { method: a.method })
        } else {
            Err(exceptions::PyValueError::new_err(
                "Not a valid vakint evaluation method",
            ))
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
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
    /// >>> from symbolica.community.vakint import VakintEvaluationMethod
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

    #[pyo3(signature = (expand_masters = None, susbstitute_masters = None, substitute_hpls = None, direct_numerical_substition = None))]
    #[classmethod]
    /// Create a new VakintEvaluationMethod instance representing the MATAD method.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.vakint import VakintEvaluationMethod
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
    /// >>> from symbolica.community.vakint import VakintEvaluationMethod
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
    /// >>> from symbolica.community.vakint import VakintEvaluationMethod
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
#[pymethods]
impl VakintWrapper {
    #[pyo3(signature = (run_time_decimal_precision = None, evaluation_order = None, epsilon_symbol = None, mu_r_sq_symbol = None, form_exe_path = None, python_exe_path = None, verify_numerator_identification = None, integral_normalization_factor = None, allow_unknown_integrals = None, clean_tmp_dir = None, number_of_terms_in_epsilon_expansion = None, use_dot_product_notation = None, temporary_directory = None))]
    #[allow(clippy::too_many_arguments)]
    #[new]
    /// Create a new Vakint instance, specifying details of the evaluation stack. Note that the same instance can be recycled across multiple evaluations.
    /// Note that the creation of a Vakint instance involves the processing and creation of the library of all known topologies, which can be time consuming.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.vakint import Vakint
    /// >>> vakint = Vakint(evaluation_order=[])
    /// >>> vakint is not None
    /// True
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
    pub fn new(
        run_time_decimal_precision: Option<u32>,
        evaluation_order: Option<Vec<PyRef<VakintEvaluationMethodWrapper>>>,
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
    ) -> PyResult<VakintWrapper> {
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
            clean_tmp_dir: clean_tmp_dir.unwrap_or(env::var("VAKINT_NO_CLEAN_TMP_DIR").is_err()),
            temporary_directory,
            evaluation_order: eval_order,
            // This quantity is typically set equal to *one plus the maximum loop count* of the UV regularisation problem considered.
            // For example when considering a 2-loop problem, then:
            //   a) for the nested one-loop integrals appearing, the single pole, finite term *and* order-epsilon term will need to be considered.
            //   b) for the two-loop integrals, the double pole, single pole and finite terms will be needed, so again three terms
            number_of_terms_in_epsilon_expansion: number_of_terms_in_epsilon_expansion.unwrap_or(4),
            use_dot_product_notation: use_dot_product_notation.unwrap_or(false),
            ..VakintSettings::default()
        };
        let vakint = Vakint::new().map_err(vakint_to_python_error)?;
        vakint
            .validate_settings(&settings)
            .map_err(vakint_to_python_error)?;
        let wrapper = VakintWrapper { vakint, settings };
        Ok(wrapper)
    }

    #[pyo3(signature = (expr))]
    /// Interpret a Symbolica expression as a numerical Laurent series in epsilon.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.vakint import Vakint
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
        let value = NumericalEvaluationResult::from_atom(
            expr.expr.as_view(),
            vakint_symbol!(self.settings.epsilon_symbol.clone()),
            &self.settings,
        )
        .map_err(vakint_to_python_error)?;
        Ok(NumericalEvaluationResultWrapper { value })
    }

    #[pyo3(signature = (evaluated_integral, params, externals = None))]
    /// Substitute numerical parameters into an integral already evaluated parametrically by Vakint.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.vakint import Vakint
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
    /// params : Dict[str, float]
    ///   A dictionary mapping parameter names to their numerical values.
    /// externals : Optional[Dict[int, Tuple[float, float, float, float]]]
    ///   An optional dictionary mapping external momentum indices to their numerical 4-vector values.
    pub fn numerical_evaluation(
        &self,
        evaluated_integral: Py<PyAny>,
        params: HashMap<String, f64>,
        externals: Option<HashMap<usize, (f64, f64, f64, f64)>>,
    ) -> PyResult<(
        NumericalEvaluationResultWrapper,
        Option<NumericalEvaluationResultWrapper>,
    )> {
        Python::attach(|py| {
            let evaluated_integral_atom = evaluated_integral.extract::<PythonExpression>(py)?;
            let p = self.vakint.params_from_f64(&self.settings, &params);
            let e = externals.map(|ext| self.vakint.externals_from_f64(&self.settings, &ext));
            let (numerical_result, numerical_error) = self
                .vakint
                .numerical_evaluation(
                    &self.settings,
                    evaluated_integral_atom.expr.as_view(),
                    &p,
                    &HashMap::default(),
                    e.as_ref(),
                )
                .map_err(vakint_to_python_error)?;
            Ok((
                NumericalEvaluationResultWrapper {
                    value: numerical_result,
                },
                numerical_error.map(|e| NumericalEvaluationResultWrapper { value: e }),
            ))
        })
    }

    #[pyo3(signature = (result))]
    /// Convert a Vakint numerical result to a Symbolica Laurent-series expression.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica.community.vakint import Vakint, VakintNumericalResult
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
        let value = result
            .value
            .to_atom(vakint_symbol!(self.settings.epsilon_symbol.clone()));

        Ok(value.into())
    }

    #[pyo3(signature = (integral_expression, short_form = None))]
    /// Convert a Vakint expression to canonical momentum routing and topology numbering.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.vakint import Vakint
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
        let result = self
            .vakint
            .to_canonical(
                &self.settings,
                integral_expression.expr.as_view(),
                short_form.unwrap_or(false),
            )
            .map_err(vakint_to_python_error)?;
        Ok(result.into())
    }

    #[pyo3(signature = (integral_expression))]
    /// Reduce the tensor integrals in a Vakint expression to scalar integrals.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.vakint import Vakint
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
        let result = self
            .vakint
            .tensor_reduce(&self.settings, integral_expression.expr.as_view())
            .map_err(vakint_to_python_error)?;
        Ok(result.into())
    }

    #[pyo3(signature = (integral_expression))]
    /// Perform the parametric evaluation of *only the integral* appearing in the Symbolica expression given in input representing a vakint integral.
    /// The numerator is left unchanged.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.vakint import Vakint, VakintEvaluationMethod
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
        let result = self
            .vakint
            .evaluate_integral(&self.settings, integral_expression.expr.as_view())
            .map_err(vakint_to_python_error)?;
        Ok(result.into())
    }

    #[pyo3(signature = (integral_expression))]
    /// Perform the complete parametric evaluation of the Vakint integral represented by the Symbolica expression given in input.
    /// Note that the tensor reduction will be automatically performed on the input given.
    ///
    /// ## Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.vakint import Vakint, VakintEvaluationMethod
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
    pub fn evaluate(&self, integral_expression: PythonExpression) -> PyResult<PythonExpression> {
        let result = self
            .vakint
            .evaluate(&self.settings, integral_expression.expr.as_view())
            .map_err(vakint_to_python_error)?;
        Ok(result.into())
    }
}
