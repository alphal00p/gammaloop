//! Python conversion and result wrapping for Idenso's typed alias registry.

use std::{collections::HashMap, sync::Arc};

use idenso::tensor::{
    ContractionSettings, SymbolicTensor, aliases::AliasInterfaces, inference::TensorInferenceError,
};
use pyo3::prelude::*;
use spenso::structure::partial::PartialStructure;
use symbolica::{
    api::python::{PythonExpression, PythonExpressionEvaluator},
    atom::{AliasedAtom, Atom, Symbol},
};

use crate::{
    ModuleInit,
    expression::TensorExpression,
    simplification::{PyColorSimplifySettings, PyGammaSimplifySettings, PySimplifySettings},
};

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
    pub(crate) value: Arc<SymbolicTensor<AliasInterfaces, AliasedAtom>>,
    descriptor: Descriptor,
    descriptors: Arc<HashMap<Atom, (Descriptor, Descriptor)>>,
}

impl ModuleInit for AliasedTensorExpression {}

impl AliasedTensorExpression {
    pub(crate) fn from_parts(
        value: Arc<SymbolicTensor<AliasInterfaces, AliasedAtom>>,
        name: Option<Symbol>,
        arguments: Vec<Atom>,
    ) -> Self {
        Self {
            value,
            descriptor: (name, arguments),
            descriptors: Arc::new(HashMap::new()),
        }
    }

    fn with_value(&self, value: Arc<SymbolicTensor<AliasInterfaces, AliasedAtom>>) -> Self {
        Self {
            value,
            descriptor: self.descriptor.clone(),
            descriptors: Arc::clone(&self.descriptors),
        }
    }

    fn wrap(
        py: Python<'_>,
        tensor: SymbolicTensor<PartialStructure>,
        descriptor: &Descriptor,
    ) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_shared(py, Arc::new(tensor), descriptor.0, descriptor.1.clone())
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
                    handle.atom().clone(),
                    (Self::descriptor(&handle), Self::descriptor(&body)),
                );
                (
                    TensorExpression::structured(&handle).clone(),
                    TensorExpression::structured(&body).clone(),
                )
            })
            .collect::<Vec<_>>();
        Ok(Self {
            value: Arc::new(
                TensorExpression::structured(&root)
                    .clone()
                    .with_aliases(aliases)
                    .map_err(TensorExpression::inference_error)?,
            ),
            descriptor: Self::descriptor(&root),
            descriptors: Arc::new(descriptors),
        })
    }

    /// Give a tensor one fresh opaque handle, retaining its original definition.
    #[staticmethod]
    fn from_expression(expression: PyRef<'_, TensorExpression>) -> PyResult<Self> {
        let body = TensorExpression::structured(&expression).clone();
        let handle = body
            .alias_handle()
            .map_err(TensorExpression::inference_error)?;
        let descriptor = Self::descriptor(&expression);
        let descriptors = HashMap::from([(
            handle.expression().clone(),
            ((None, Vec::new()), descriptor.clone()),
        )]);
        Ok(Self {
            value: Arc::new(
                handle
                    .clone()
                    .with_aliases([(handle, body)])
                    .map_err(TensorExpression::inference_error)?,
            ),
            descriptor,
            descriptors: Arc::new(descriptors),
        })
    }

    #[getter]
    fn root(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        Self::wrap(py, self.value.root(), &self.descriptor)
    }

    /// Whether the metric/vector contractor certified completion.
    /// False also covers an uncontracted value or an exact retained frontier
    /// stopped by its budget; explicit materialization is a separate operation.
    #[getter]
    fn contraction_complete(&self) -> bool {
        self.value.contraction_complete()
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
                    self.descriptors.get(handle.expression()).unwrap_or(&empty);
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
        self.value.expression().get_byte_size()
    }

    fn __repr__(&self) -> String {
        format!(
            "AliasedTensorExpression(root={}, definitions={}, bytes={})",
            crate::display::format_structured(&self.value.root(), false),
            self.value.expression().get_aliases().len(),
            self.get_byte_size(),
        )
    }

    /// Apply a Python callback once per definition, retaining its typed interface.
    fn map_aliases(&self, py: Python<'_>, callback: Py<PyAny>) -> PyResult<Self> {
        let mut callback_error = None;
        let mut descriptors = (*self.descriptors).clone();
        let value = self.value.map_aliases(|handle, body| {
            let result = (|| {
                let empty = ((None, Vec::new()), (None, Vec::new()));
                let descriptor = self.descriptors.get(handle.expression()).unwrap_or(&empty);
                let argument = Self::wrap(py, body, &descriptor.1)?;
                let output = callback.call1(py, (argument,))?;
                let tensor = output.extract::<PyRef<'_, TensorExpression>>(py)?;
                descriptors.insert(
                    handle.expression().clone(),
                    (descriptor.0.clone(), Self::descriptor(&tensor)),
                );
                Ok::<_, PyErr>(TensorExpression::structured(&tensor).clone())
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
            value: Arc::new(value.map_err(TensorExpression::inference_error)?),
            descriptor: self.descriptor.clone(),
            descriptors: Arc::new(descriptors),
        })
    }

    /// Apply ordered tensor rules to the root and reachable definitions without resolving them.
    /// The shared owner preserves first-match priority and per-call callback caches.
    #[pyo3(signature=(rules))]
    fn replace(&self, rules: &Bound<'_, PyAny>) -> PyResult<Self> {
        let rules = crate::tensor_rule::PyTensorRule::extract_rules(rules)?;
        let value = self
            .value
            .replace_rules(&rules)
            .map_err(TensorExpression::inference_error)?;
        Ok(self.with_value(value))
    }

    /// Contract the root and reachable definitions, reusing certified completed results.
    #[pyo3(signature=(order=None, *, rank_one=true))]
    fn contract(&self, order: Option<Vec<usize>>, rank_one: bool) -> PyResult<Self> {
        let mut settings = ContractionSettings::default();
        if let Some(order) = &order {
            settings = settings.with_order(order);
        }
        if !rank_one {
            settings = settings.without_rank_one_tensors();
        }
        let value = self
            .value
            .contract(settings)
            .map_err(TensorExpression::inference_error)?;
        Ok(self.with_value(value))
    }

    /// Render compact metric products in the root and reachable definitions.
    fn to_dots(&self) -> PyResult<Self> {
        self.value
            .to_dots()
            .map(|value| self.with_value(value))
            .map_err(TensorExpression::inference_error)
    }

    /// Open dots into indexed contractions without resolving the alias registry.
    fn undo_dots(&self) -> PyResult<Self> {
        self.value
            .undo_dots()
            .map(|value| self.with_value(value))
            .map_err(TensorExpression::inference_error)
    }

    /// Apply shared tensor identities to the root and reachable definitions.
    #[pyo3(signature=(settings=None))]
    fn simplify(&self, settings: Option<&PySimplifySettings>) -> PyResult<Self> {
        let settings = settings.copied().unwrap_or_default();
        let value = self
            .value
            .simplify(settings.rust())
            .map_err(TensorExpression::inference_error)?;
        Ok(self.with_value(value))
    }

    /// Apply Dirac identities while retaining tensor-valued trace definitions.
    #[pyo3(signature=(settings=None))]
    fn simplify_gamma(&self, settings: Option<&PyGammaSimplifySettings>) -> PyResult<Self> {
        let settings = settings.map_or_else(Default::default, PyGammaSimplifySettings::rust);
        let value = self
            .value
            .simplify_gamma(settings)
            .map_err(TensorExpression::inference_error)?;
        Ok(self.with_value(value))
    }

    /// Apply color identities to the root and reachable definitions.
    #[pyo3(signature=(settings=None))]
    fn simplify_color(&self, settings: Option<&PyColorSimplifySettings>) -> PyResult<Self> {
        let settings = settings.map_or_else(Default::default, PyColorSimplifySettings::rust);
        let value = self
            .value
            .simplify_color(settings)
            .map_err(TensorExpression::inference_error)?;
        Ok(self.with_value(value))
    }

    /// Apply epsilon identities to the root and reachable definitions.
    fn simplify_epsilon(&self) -> PyResult<Self> {
        let value = self
            .value
            .simplify_epsilon()
            .map_err(TensorExpression::inference_error)?;
        Ok(self.with_value(value))
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
        let builder = self
            .value
            .evaluator(&params)
            .map_err(TensorExpression::inference_error)?;
        TensorExpression::build_evaluator(py, builder, iterations, n_cores)
    }
}
