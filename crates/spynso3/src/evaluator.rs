//! Tensor-shaped results backed by Symbolica's public Python evaluator.
use std::collections::HashMap;

use pyo3::{
    exceptions::PyValueError,
    prelude::*,
    types::{PyList, PyTuple},
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::*;
#[cfg(not(feature = "python_stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;
use symbolica::api::python::{
    PythonExpression, PythonExpressionEvaluator, PythonFunctionDefinition,
};
use symbolica::atom::{Atom, Symbol};

use crate::{
    AbstractIndex, Complex, DataTensor, DenseTensor, HasStructure, MixedTensor, PartialStructure,
    RealOrComplexTensor, ShadowedStructure, Spensor, SymbolicTensor, tensor_data_layout,
    tensor_data_layout_error,
};

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl Spensor {
    /// Optimise all components together using Symbolica's evaluator.
    ///
    /// Accepts the same parameters, FunctionDefinitions, optimisation controls,
    /// and JIT settings as Expression.evaluator. Fixed values can be substituted
    /// before construction. Evaluation returns one Tensor per input row, retaining
    /// the logical axes and data identity. JIT compilation occurs on first use.
    ///
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> values.evaluator([x]).evaluate([[2.0]])[0][1]
    /// 4.0
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (params, functions=Vec::default(), iterations=1, cpe_iterations=None,
        n_cores=4, verbose=false, jit_compile=true, direct_translation=true,
        jit_direct_translation=false, jit_optimization_level=3, jit_options=HashMap::default(),
        max_horner_scheme_variables=500, max_common_pair_cache_entries=1_000_000,
        max_common_pair_distance=100))]
    pub fn evaluator(
        &self,
        params: Vec<PythonExpression>,
        functions: Vec<PythonFunctionDefinition>,
        iterations: usize,
        cpe_iterations: Option<usize>,
        n_cores: usize,
        verbose: bool,
        jit_compile: bool,
        direct_translation: bool,
        jit_direct_translation: bool,
        jit_optimization_level: u8,
        jit_options: HashMap<String, String>,
        max_horner_scheme_variables: usize,
        max_common_pair_cache_entries: usize,
        max_common_pair_distance: usize,
        py: Python<'_>,
    ) -> PyResult<SpensoExpressionEvaluator> {
        // Public component order may differ from canonical storage order.
        let expressions = (0..self.__len__())
            .map(|i| {
                self.__getitem__(crate::SliceOrIntOrExpanded::Int(i))?
                    .bind(py)
                    .extract::<symbolica::api::python::ConvertibleToExpression>()
                    .map(|value| value.to_expression())
            })
            .collect::<PyResult<Vec<_>>>()?;
        let evaluator = PythonExpression::evaluator_multiple(
            &py.get_type::<PythonExpression>(),
            expressions,
            params.clone(),
            functions,
            iterations,
            cpe_iterations,
            n_cores,
            verbose,
            jit_compile,
            direct_translation,
            jit_direct_translation,
            jit_optimization_level,
            jit_options.into_iter().collect(),
            max_horner_scheme_variables,
            max_common_pair_cache_entries,
            max_common_pair_distance,
            py,
        )?;
        Ok(SpensoExpressionEvaluator {
            evaluator: Py::new(py, evaluator)?,
            layout: TensorEvaluationLayout {
                storage: self.tensor.structure().clone(),
                descriptor: self.descriptor.clone(),
                name: self.descriptor_name,
                args: self.descriptor_args.clone(),
                parameters: params,
            },
        })
    }
}

#[derive(Clone)]
struct TensorEvaluationLayout {
    storage: ShadowedStructure<AbstractIndex>,
    descriptor: SymbolicTensor<PartialStructure>,
    name: Option<Symbol>,
    args: Vec<Atom>,
    parameters: Vec<PythonExpression>,
}

impl TensorEvaluationLayout {
    fn wrap(&self, rows: &Bound<'_, PyAny>, complex: bool) -> PyResult<Vec<Spensor>> {
        let layout = tensor_data_layout(self.descriptor.structure())?;
        // Symbolica returns nested lists when NumPy is unavailable.
        let rows = if rows.is_instance_of::<PyList>() {
            rows.clone()
        } else {
            rows.call_method0("tolist")?
        };
        let tensors: Vec<MixedTensor<f64, ShadowedStructure<AbstractIndex>>> = if complex {
            rows.extract::<Vec<Vec<Complex<f64>>>>()?
                .into_iter()
                .map(|row| {
                    let data = layout
                        .reorder_to_storage(row)
                        .map_err(tensor_data_layout_error)?;
                    let dense = DenseTensor::from_storage_data(data, self.storage.clone())
                        .map_err(|e| PyValueError::new_err(e.to_string()))?;
                    Ok(MixedTensor::Concrete(RealOrComplexTensor::Complex(
                        DataTensor::Dense(dense),
                    )))
                })
                .collect::<PyResult<_>>()?
        } else {
            rows.extract::<Vec<Vec<f64>>>()?
                .into_iter()
                .map(|row| {
                    let data = layout
                        .reorder_to_storage(row)
                        .map_err(tensor_data_layout_error)?;
                    let dense = DenseTensor::from_storage_data(data, self.storage.clone())
                        .map_err(|e| PyValueError::new_err(e.to_string()))?;
                    Ok(dense.into())
                })
                .collect::<PyResult<_>>()?
        };
        Ok(tensors
            .into_iter()
            .map(|tensor| {
                Spensor::from_storage_with_descriptor(
                    tensor,
                    self.descriptor.clone(),
                    self.name,
                    self.args.clone(),
                )
            })
            .collect())
    }

    fn output_shape<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(
            py,
            tensor_data_layout(self.descriptor.structure())?.logical_shape(),
        )
    }
}

/// Symbolica evaluator returning component tensors in their original logical layout.
///
/// Construct with Tensor.evaluator. The scalar_evaluator property exposes the
/// underlying Symbolica Evaluator for precision evaluation, export and inspection.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "TensorEvaluator", module = "symbolica.community.tensor")]
pub struct SpensoExpressionEvaluator {
    evaluator: Py<PythonExpressionEvaluator>,
    layout: TensorEvaluationLayout,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl SpensoExpressionEvaluator {
    /// Symbolic inputs in evaluation order.
    #[getter]
    fn parameters(&self) -> Vec<PythonExpression> {
        self.layout.parameters.clone()
    }

    /// Number of entries in an input row.
    #[getter]
    fn input_size(&self) -> usize {
        self.layout.parameters.len()
    }

    /// Logical shape of each output tensor, excluding the batch dimension.
    #[getter]
    fn output_shape<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        self.layout.output_shape(py)
    }

    /// Whether all exact coefficients are real.
    #[getter]
    fn supports_real(&self, py: Python<'_>) -> bool {
        self.evaluator
            .borrow(py)
            .rational_constants
            .iter()
            .all(|c| c.im.is_zero())
    }

    /// Underlying Symbolica Evaluator, with outputs in logical component order.
    ///
    /// Use for arbitrary-precision evaluation, export, or instruction inspection.
    /// These operations return scalar components; they do not wrap tensor metadata.
    #[getter]
    fn scalar_evaluator(&self, py: Python<'_>) -> Py<PythonExpressionEvaluator> {
        self.evaluator.clone_ref(py)
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "TensorEvaluator({}, real={}, complex=True)",
            crate::display::format_structured(&self.layout.descriptor, false),
            if self.supports_real(py) {
                "True"
            } else {
                "False"
            }
        )
    }

    /// Evaluate real inputs using Symbolica's ArrayLike input conventions.
    ///
    /// Returns one Tensor per row. A flat input array is also accepted, with
    /// the same reshaping and validation as Evaluator.evaluate.
    fn evaluate(
        &self,
        #[gen_stub(override_type(type_repr = "numpy.typing.ArrayLike", imports = ("numpy.typing",)))]
        inputs: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<Spensor>> {
        self.layout.wrap(
            &self
                .evaluator
                .bind(inputs.py())
                .call_method1("evaluate", (inputs,))?,
            false,
        )
    }

    /// Evaluate complex inputs; returns one Tensor per input row.
    fn evaluate_complex(
        &self,
        #[gen_stub(override_type(type_repr = "numpy.typing.ArrayLike", imports = ("numpy.typing",)))]
        inputs: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<Spensor>> {
        self.layout.wrap(
            &self
                .evaluator
                .bind(inputs.py())
                .call_method1("evaluate_complex", (inputs,))?,
            true,
        )
    }

    /// Enable or disable JIT compilation with Symbolica's settings.
    #[pyo3(signature = (jit_compile, direct_translation=None, optimization_level=None, options=None))]
    fn jit_compile(
        &self,
        py: Python<'_>,
        jit_compile: bool,
        direct_translation: Option<bool>,
        optimization_level: Option<u8>,
        options: Option<HashMap<String, String>>,
    ) -> PyResult<()> {
        self.evaluator.bind(py).call_method1(
            "jit_compile",
            (jit_compile, direct_translation, optimization_level, options),
        )?;
        Ok(())
    }

    /// Mark real parameters using the same assumptions as Evaluator.set_real_params.
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (real_params, sqrt_real=false, log_real=false, powf_real=false, real_if_args_real=false, verbose=false))]
    fn set_real_params(
        &self,
        py: Python<'_>,
        real_params: Vec<usize>,
        sqrt_real: bool,
        log_real: bool,
        powf_real: bool,
        real_if_args_real: bool,
        verbose: bool,
    ) -> PyResult<()> {
        self.evaluator.bind(py).call_method1(
            "set_real_params",
            (
                real_params,
                sqrt_real,
                log_real,
                powf_real,
                real_if_args_real,
                verbose,
            ),
        )?;
        Ok(())
    }
}

#[cfg(feature = "native")]
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl SpensoExpressionEvaluator {
    /// Compile a native evaluator with the same options as Evaluator.compile.
    ///
    /// number_type selects real, complex, real_4x, complex_4x, cuda_real or
    /// cuda_complex. The returned object's evaluate method uses that number type
    /// and returns tensors. Source and library files are written to the given paths.
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (function_name, filename, library_name, number_type,
        inline_asm="default", optimization_level=3, native=true, compiler_path=None,
        compiler_flags=None, custom_header=None, cuda_number_of_evaluations=1, cuda_block_size=512))]
    fn compile(
        &self,
        py: Python<'_>,
        function_name: &str,
        filename: &str,
        library_name: &str,
        number_type: &str,
        inline_asm: &str,
        optimization_level: u8,
        native: bool,
        compiler_path: Option<&str>,
        compiler_flags: Option<Vec<String>>,
        custom_header: Option<String>,
        cuda_number_of_evaluations: usize,
        cuda_block_size: usize,
    ) -> PyResult<SpensoCompiledExpressionEvaluator> {
        let evaluator = self
            .evaluator
            .bind(py)
            .call_method1(
                "compile",
                (
                    function_name,
                    filename,
                    library_name,
                    number_type,
                    inline_asm,
                    optimization_level,
                    native,
                    compiler_path,
                    compiler_flags,
                    custom_header,
                    cuda_number_of_evaluations,
                    cuda_block_size,
                ),
            )?
            .unbind();
        Ok(SpensoCompiledExpressionEvaluator {
            evaluator,
            layout: self.layout.clone(),
            complex: number_type.contains("complex"),
        })
    }
}

/// Native Symbolica evaluator returning tensors. Construct with TensorEvaluator.compile.
///
/// evaluate uses the number_type selected at compilation, as in Symbolica's compiled
/// evaluators. scalar_evaluator exposes the underlying native evaluator.
#[cfg(feature = "native")]
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CompiledTensorEvaluator",
    module = "symbolica.community.tensor"
)]
pub struct SpensoCompiledExpressionEvaluator {
    evaluator: Py<PyAny>,
    layout: TensorEvaluationLayout,
    complex: bool,
}

#[cfg(feature = "native")]
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl SpensoCompiledExpressionEvaluator {
    /// Symbolic inputs in evaluation order.
    #[getter]
    fn parameters(&self) -> Vec<PythonExpression> {
        self.layout.parameters.clone()
    }

    /// Number of entries in each input row.
    #[getter]
    fn input_size(&self) -> usize {
        self.layout.parameters.len()
    }

    /// Logical shape of each result tensor.
    #[getter]
    fn output_shape<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        self.layout.output_shape(py)
    }

    /// Whether compilation selected a real number type.
    #[getter]
    fn supports_real(&self) -> bool {
        !self.complex
    }

    /// Underlying Symbolica compiled evaluator.
    #[getter]
    fn scalar_evaluator(&self, py: Python<'_>) -> Py<PyAny> {
        self.evaluator.clone_ref(py)
    }

    /// Evaluate real or complex inputs, according to the compiled number type.
    fn evaluate(
        &self,
        #[gen_stub(override_type(type_repr = "numpy.typing.ArrayLike", imports = ("numpy.typing",)))]
        inputs: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<Spensor>> {
        self.layout.wrap(
            &self
                .evaluator
                .bind(inputs.py())
                .call_method1("evaluate", (inputs,))?,
            self.complex,
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "CompiledTensorEvaluator({}, number_type={})",
            crate::display::format_structured(&self.layout.descriptor, false),
            if self.complex { "complex" } else { "real" }
        )
    }
}
