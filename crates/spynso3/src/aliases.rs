//! Python conversion and result wrapping for Idenso's typed alias registry.

use std::collections::HashMap;

use idenso::tensor::{SymbolicTensor, aliases::AliasInterfaces, inference::TensorInferenceError};
use pyo3::{exceptions::PyValueError, prelude::*};
use spenso::structure::partial::PartialStructure;
use symbolica::{
    api::python::{PythonExpression, PythonExpressionEvaluator},
    atom::{AliasedAtom, Atom, Symbol},
    domains::float::Complex,
    evaluate::OptimizationSettings,
};

use crate::{ModuleInit, expression::TensorExpression};

type Descriptor = (Option<Symbol>, Vec<Atom>);

/// A typed root and literal tensor definitions, materialized only on request.
/// Scalar results use this same class. `evaluator` consumes the DAG directly.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(
    frozen,
    name = "AliasedTensorExpression",
    module = "symbolica.community.spenso"
)]
pub struct AliasedTensorExpression {
    pub(crate) value: SymbolicTensor<AliasInterfaces, AliasedAtom>,
    descriptor: Descriptor,
    descriptors: HashMap<Atom, (Descriptor, Descriptor)>,
}

impl ModuleInit for AliasedTensorExpression {}

impl AliasedTensorExpression {
    fn wrap(
        py: Python<'_>,
        tensor: SymbolicTensor<PartialStructure>,
        descriptor: &Descriptor,
    ) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_parts_unchecked(
            py,
            tensor.expression,
            tensor.structure,
            descriptor.0,
            descriptor.1.clone(),
        )
    }

    fn descriptor(value: &PyRef<'_, TensorExpression>) -> Descriptor {
        (value.name, value.name_args.clone())
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl AliasedTensorExpression {
    #[new]
    #[pyo3(signature = (root, aliases=Vec::new()))]
    fn new(
        py: Python<'_>,
        root: PyRef<'_, TensorExpression>,
        aliases: Vec<(Py<TensorExpression>, Py<TensorExpression>)>,
    ) -> PyResult<Self> {
        let mut descriptors = HashMap::new();
        let aliases = aliases
            .into_iter()
            .map(|(handle, body)| {
                let handle = handle.borrow(py);
                let body = body.borrow(py);
                descriptors.insert(
                    handle.as_super().expr.clone(),
                    (Self::descriptor(&handle), Self::descriptor(&body)),
                );
                (
                    TensorExpression::structured(&handle),
                    TensorExpression::structured(&body),
                )
            })
            .collect::<Vec<_>>();
        Ok(Self {
            value: TensorExpression::structured(&root)
                .with_aliases(aliases)
                .map_err(TensorExpression::inference_error)?,
            descriptor: Self::descriptor(&root),
            descriptors,
        })
    }

    /// Give a tensor one fresh opaque handle, retaining its original definition.
    #[staticmethod]
    fn from_expression(expression: PyRef<'_, TensorExpression>) -> PyResult<Self> {
        let body = TensorExpression::structured(&expression);
        let handle = body
            .alias_handle()
            .map_err(TensorExpression::inference_error)?;
        let descriptor = Self::descriptor(&expression);
        let descriptors = HashMap::from([(
            handle.expression.clone(),
            ((None, Vec::new()), descriptor.clone()),
        )]);
        Ok(Self {
            value: handle
                .clone()
                .with_aliases([(handle, body)])
                .map_err(TensorExpression::inference_error)?,
            descriptor,
            descriptors,
        })
    }

    #[getter]
    fn root(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        Self::wrap(py, self.value.root(), &self.descriptor)
    }

    #[getter]
    fn aliases(
        &self,
        py: Python<'_>,
    ) -> PyResult<Vec<(Py<TensorExpression>, Py<TensorExpression>)>> {
        self.value
            .aliases()
            .map_err(TensorExpression::inference_error)?
            .into_iter()
            .map(|(handle, body)| {
                let empty = ((None, Vec::new()), (None, Vec::new()));
                let (handle_descriptor, body_descriptor) =
                    self.descriptors.get(&handle.expression).unwrap_or(&empty);
                Ok((
                    Self::wrap(py, handle, handle_descriptor)?,
                    Self::wrap(py, body, body_descriptor)?,
                ))
            })
            .collect()
    }

    fn to_expression(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        Self::wrap(
            py,
            self.value
                .resolved()
                .map_err(TensorExpression::inference_error)?,
            &self.descriptor,
        )
    }

    fn expand(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        Self::wrap(
            py,
            self.value
                .expanded()
                .map_err(TensorExpression::inference_error)?,
            &self.descriptor,
        )
    }

    fn get_byte_size(&self) -> usize {
        self.value.expression.get_byte_size()
    }

    /// Apply a Python callback once per definition, retaining its typed interface.
    fn map_aliases(&self, py: Python<'_>, callback: Py<PyAny>) -> PyResult<Self> {
        let mut callback_error = None;
        let mut descriptors = self.descriptors.clone();
        let value = self.value.map_aliases(|handle, body| {
            let result = (|| {
                let empty = ((None, Vec::new()), (None, Vec::new()));
                let descriptor = self.descriptors.get(&handle.expression).unwrap_or(&empty);
                let argument = Self::wrap(py, body, &descriptor.1)?;
                let output = callback.call1(py, (argument,))?;
                let tensor = output.extract::<PyRef<'_, TensorExpression>>(py)?;
                descriptors.insert(
                    handle.expression.clone(),
                    (descriptor.0.clone(), Self::descriptor(&tensor)),
                );
                Ok::<_, PyErr>(TensorExpression::structured(&tensor))
            })();
            result.map_err(|error| {
                callback_error = Some(error);
                TensorInferenceError::Invalid("alias callback failed".into())
            })
        });
        if let Some(error) = callback_error {
            return Err(error);
        }
        Ok(Self {
            value: value.map_err(TensorExpression::inference_error)?,
            descriptor: self.descriptor.clone(),
            descriptors,
        })
    }

    /// Build Symbolica's evaluator directly from the root and its definitions.
    #[pyo3(signature = (params, *, iterations=1, n_cores=1))]
    fn evaluator(
        &self,
        py: Python<'_>,
        params: Vec<PythonExpression>,
        iterations: usize,
        n_cores: usize,
    ) -> PyResult<PythonExpressionEvaluator> {
        let params = params
            .into_iter()
            .map(|value| value.expr)
            .collect::<Vec<_>>();
        let settings = OptimizationSettings::new()
            .horner_iterations(iterations)
            .cores(n_cores);
        let evaluator = py
            .detach(|| {
                self.value
                    .evaluator(&params)
                    .map_err(|error| error.to_string())?
                    .optimization_settings(settings)
                    .build()
                    .map_err(|error| error.to_string())
            })
            .map_err(PyValueError::new_err)?;
        Ok(PythonExpressionEvaluator {
            rational_constants: evaluator.get_constants().to_vec(),
            eval_complex: evaluator.map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64())),
            eval_real: None,
            #[cfg(feature = "native")]
            jit_real: None,
            #[cfg(feature = "native")]
            jit_complex: None,
            eval_double_float: None,
            eval_double_float_complex: None,
            eval_arb_prec: None,
            eval_arb_prec_complex: None,
            jit_compile: false,
            #[cfg(feature = "native")]
            jit_settings: symbolica::evaluate::JITCompilationSettings::new(),
        })
    }
}
