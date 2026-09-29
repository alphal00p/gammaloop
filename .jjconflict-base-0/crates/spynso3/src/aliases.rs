//! Python conversion and result wrapping for Idenso's typed alias registry.

use std::{collections::HashMap, sync::Arc};

use idenso::tensor::{SymbolicTensor, aliases::AliasInterfaces, inference::TensorInferenceError};
use pyo3::prelude::*;
use spenso::structure::partial::PartialStructure;
use symbolica::{
    api::python::{PythonExpression, PythonExpressionEvaluator},
    atom::{AliasedAtom, Atom, Symbol},
};

use crate::{ModuleInit, expression::TensorExpression, simplification::PyAlgebraSettings};

type Descriptor = (Option<Symbol>, Vec<Atom>);

/// A typed root and literal tensor definitions, materialized only on request.
/// Scalar results use this same class. `evaluator` consumes the DAG directly.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> from symbolica import E, S
/// >>> x = S("docs::x")
/// >>> handle = sp.TensorExpression(S("docs::handle"))
/// >>> body = sp.TensorExpression((x + 1)**2)
/// >>> value = sp.AliasedTensorExpression(handle, [(handle, body)])
/// >>> resolved = value.to_expression()
///
/// Parameters
/// ----------
/// root : TensorExpression
///     Tensor expression in which the alias handles occur.
/// aliases : Sequence[tuple[TensorExpression, TensorExpression]]
///     Sequence of (handle, body) TensorExpression pairs. Each body must have the same external interface as its handle.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(
    frozen,
    name = "AliasedTensorExpression",
    module = "symbolica.community.tensor"
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
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl AliasedTensorExpression {
    /// Build a factored expression from a root and named subexpressions.
    /// Each `(handle, body)` pair defines one reusable subexpression with the same tensor interface.
    /// Definitions may refer to other handles, but cyclic definitions are rejected.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> handle = sp.TensorExpression(S("docs::handle"))
    /// >>> body = sp.TensorExpression((x + 1)**2)
    /// >>> value = sp.AliasedTensorExpression(handle, [(handle, body)])
    /// >>> resolved = value.to_expression()
    ///
    /// Parameters
    /// ----------
    /// root : TensorExpression
    ///     Tensor expression in which the alias handles occur.
    /// aliases : Sequence[tuple[TensorExpression, TensorExpression]]
    ///     Sequence of (handle, body) TensorExpression pairs. Each body must have the same external interface as its handle.
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
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> body = sp.TensorExpression((x + 1)**2)
    /// >>> value = sp.AliasedTensorExpression.from_expression(body)
    /// >>> assert value.to_expression() == body
    ///
    /// Parameters
    /// ----------
    /// expression : TensorExpression
    ///     TensorExpression to store under one fresh alias handle.
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

    /// Tensor expression in which the alias handles occur.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> value = sp.AliasedTensorExpression.from_expression(sp.TensorExpression((x + 1)**2))
    /// >>> result = value.root
    #[getter]
    fn root(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        Self::wrap(py, self.value.root(), &self.descriptor)
    }

    /// Whether the metric/vector contractor certified completion.
    /// False also covers an uncontracted value or an exact retained frontier
    /// stopped by its budget; explicit materialization is a separate operation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> value = sp.AliasedTensorExpression.from_expression(sp.TensorExpression((x + 1)**2))
    /// >>> result = value.contraction_complete
    #[getter]
    fn contraction_complete(&self) -> bool {
        self.value.contraction_complete()
    }

    /// Sequence of (handle, body) TensorExpression pairs. Each body must have the same external interface as its handle.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> value = sp.AliasedTensorExpression.from_expression(sp.TensorExpression((x + 1)**2))
    /// >>> result = value.aliases
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

    /// Replace every alias by its definition without expanding scalar products or powers.
    /// The result is a TensorExpression with the same external indices.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> value = sp.AliasedTensorExpression.from_expression(sp.TensorExpression((x + 1)**2))
    /// >>> resolved = value.to_expression()
    /// >>> assert resolved.to_expression() == (x + 1)**2
    fn to_expression(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        Self::wrap(
            py,
            self.value
                .resolved()
                .map_err(TensorExpression::inference_error)?,
            &self.descriptor,
        )
    }

    /// Resolve aliases and expand the resulting scalar algebra into explicit terms.
    /// Use `.to_expression()` instead when you want to keep products and powers factored.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> value = sp.AliasedTensorExpression.from_expression(sp.TensorExpression((x + 1)**2))
    /// >>> expanded = value.expand()
    /// >>> assert expanded.to_expression() == ((x + 1)**2).expand()
    fn expand(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        Self::wrap(
            py,
            self.value
                .expanded()
                .map_err(TensorExpression::inference_error)?,
            &self.descriptor,
        )
    }

    /// Report the storage size in bytes of the factored symbolic expression.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> value = sp.AliasedTensorExpression.from_expression(sp.TensorExpression((x + 1)**2))
    /// >>> bytes_used = value.get_byte_size()
    /// >>> assert bytes_used > 0
    fn get_byte_size(&self) -> usize {
        self.value.expression().get_byte_size()
    }

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> value = sp.AliasedTensorExpression.from_expression(sp.TensorExpression((x + 1)**2))
    /// >>> text = repr(value)
    fn __repr__(&self) -> String {
        format!(
            "AliasedTensorExpression(root={}, definitions={}, bytes={})",
            crate::display::format_structured(&self.value.root(), false),
            self.value.expression().get_aliases().len(),
            self.get_byte_size(),
        )
    }

    /// Apply a Python callback once per definition, retaining its typed interface.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> value = sp.AliasedTensorExpression.from_expression(sp.TensorExpression((x + 1)**2))
    /// >>> mapped = value.map_aliases(lambda body: body.expand())
    ///
    /// Parameters
    /// ----------
    /// callback : Any
    ///     Function receiving one TensorExpression definition and returning its replacement. Keep its external interface unchanged; Python exceptions propagate to the caller.
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

    /// Apply ordered tensor rules to a factored expression and its definitions.
    /// Aliases remain factored; call `.to_expression()` or `.expand()` to materialize them.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> value = sp.AliasedTensorExpression.from_expression(sp.TensorExpression((x + 1)**2))
    /// >>> changed = value.replace(sp.TensorRule(x, E("2")))
    ///
    /// Parameters
    /// ----------
    /// rules : Any
    ///     One TensorRule or an ordered list of rules. The first matching rule wins; the replacement is checked against the matched tensor interface.
    #[pyo3(signature=(rules))]
    fn replace(&self, rules: &Bound<'_, PyAny>) -> PyResult<Self> {
        let rules = crate::tensor_rule::PyTensorRule::extract_rules(rules)?;
        let value = self
            .value
            .replace_rules(&rules)
            .map_err(TensorExpression::inference_error)?;
        Ok(self.with_value(value))
    }

    /// Whether selected reduction work is complete, deferred, or budget-capped.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> value = sp.AliasedTensorExpression.from_expression(sp.TensorExpression((x + 1)**2))
    /// >>> result = value.reduction_status
    #[getter]
    fn reduction_status(&self) -> crate::simplification::PyReductionStatus {
        self.value.reduction_status().into()
    }

    /// Contract compatible indices, retaining shared subexpressions in the result.
    /// Use `ContractSettings` to select representations or control chain and trace collection.
    /// Call `.expand()` to materialize the factored result.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> metric = sp.TensorExpression.g(r)("i", "j")
    /// >>> vector = sp.TensorName.vector("docs::p")(r("j"))
    /// >>> value = sp.AliasedTensorExpression.from_expression(metric * vector)
    /// >>> contracted = value.contract()
    /// >>> explicit = contracted.expand()
    ///
    /// Parameters
    /// ----------
    /// settings : ContractSettings, optional
    ///     ContractSettings selecting contractions and collection. None enables the default contraction policy.
    #[pyo3(signature=(settings=None))]
    fn contract(
        &self,
        settings: Option<&crate::simplification::PyContractSettings>,
    ) -> PyResult<Self> {
        let settings = settings.cloned().unwrap_or_default();
        let value = self
            .value
            .contract(settings.rust())
            .map_err(TensorExpression::inference_error)?;
        Ok(self.with_value(value))
    }

    /// Normalize scalar-product notation in the root and reachable definitions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> p = sp.TensorName.vector("docs::p")(r)
    /// >>> q = sp.TensorName.vector("docs::q")(r)
    /// >>> value = sp.dot(p, q)
    /// >>> value = sp.AliasedTensorExpression.from_expression(value)
    /// >>> rewritten = value.to_dots()
    fn to_dots(&self) -> PyResult<Self> {
        self.value
            .to_dots()
            .map(|value| self.with_value(value))
            .map_err(TensorExpression::inference_error)
    }

    /// Open dots into indexed contractions without resolving the alias registry.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> p = sp.TensorName.vector("docs::p")(r)
    /// >>> q = sp.TensorName.vector("docs::q")(r)
    /// >>> value = sp.dot(p, q)
    /// >>> value = sp.AliasedTensorExpression.from_expression(value)
    /// >>> rewritten = value.undo_dots()
    fn undo_dots(&self) -> PyResult<Self> {
        self.value
            .undo_dots()
            .map(|value| self.with_value(value))
            .map_err(TensorExpression::inference_error)
    }

    /// Expose collected chain factors with fresh compatible dummy indices.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> value = sp.chain(r("i"), r("j"), A, A)
    /// >>> value = sp.AliasedTensorExpression.from_expression(value)
    /// >>> opened = value.undo_chain()
    fn undo_chain(&self) -> PyResult<Self> {
        self.value
            .undo_chain()
            .map(|value| self.with_value(value))
            .map_err(TensorExpression::inference_error)
    }

    /// Open traces as closed chains without evaluating them.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> value = sp.trace(r, A, A)
    /// >>> value = sp.AliasedTensorExpression.from_expression(value)
    /// >>> opened = value.undo_trace()
    fn undo_trace(&self) -> PyResult<Self> {
        self.value
            .undo_trace()
            .map(|value| self.with_value(value))
            .map_err(TensorExpression::inference_error)
    }

    /// Apply the selected tensor identities while keeping intermediate expressions factored.
    /// Choose `AlgebraSettings.hep()` for gamma, color, and epsilon identities.
    /// The result retains shared subexpressions; call `.expand()` when an explicit result is needed.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> gamma = sp.TensorExpression.gamma(4)
    /// >>> trace = gamma("a", "b", "mu") * gamma("b", "a", "nu")
    /// >>> value = sp.AliasedTensorExpression.from_expression(trace)
    /// >>> reduced = value.reduce_algebra(sp.AlgebraSettings.hep())
    /// >>> expanded = reduced.expand()
    ///
    /// Parameters
    /// ----------
    /// settings : AlgebraSettings, optional
    ///     AlgebraSettings selecting identity families. None uses the default settings with no optional gamma, color, or epsilon family enabled.
    #[pyo3(signature=(settings=None))]
    fn reduce_algebra(&self, settings: Option<&PyAlgebraSettings>) -> PyResult<Self> {
        let settings = settings.copied().unwrap_or_default();
        let value = self
            .value
            .reduce_algebra(settings.rust())
            .map_err(TensorExpression::inference_error)?;
        Ok(self.with_value(value))
    }

    /// Build a scalar numerical evaluator directly from the factored expression.
    /// Shared definitions remain available to the evaluator without first expanding
    /// the expression. Input rows follow the order in `params`.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> value = sp.AliasedTensorExpression.from_expression(sp.TensorExpression((x + 1)**2))
    /// >>> evaluator = value.evaluator([x])
    /// >>> results = evaluator.evaluate([[2.0], [3.0]])
    ///
    /// Parameters
    /// ----------
    /// params : Sequence[Expression]
    ///     Symbolica variables in the exact order expected in each numerical input row. Include every free numerical parameter that remains after substitutions.
    /// iterations : int
    ///     Number of optimization iterations when constructing the evaluator. More iterations can improve repeated evaluations at a higher setup cost.
    /// n_cores : int
    ///     Worker count used while optimizing the evaluator. Use 1 for serial setup.
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
