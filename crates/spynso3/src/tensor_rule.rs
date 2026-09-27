use idenso::tensor::TensorRule;
use pyo3::prelude::*;
use symbolica::api::python::{ConvertibleToPattern, ConvertibleToReplaceWith};

use crate::{ModuleInit, expression::TensorExpression, network::ReplacementCondition};

/// A reusable whole-tensor replacement rule.
///
/// Construction checks wildcard closure and retains proofs for literal right-hand
/// sides. Binding-dependent interfaces are checked when the rule is applied.
/// Conditions run for each candidate match. RHS callbacks are cached only within
/// one application; use `rhs_cache_size=0` for callbacks with side effects.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(frozen, name = "TensorRule", module = "symbolica.community.spenso")]
pub struct PyTensorRule {
    pub(crate) rule: TensorRule<'static>,
}

impl ModuleInit for PyTensorRule {}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyTensorRule {
    #[new]
    #[pyo3(signature = (pattern, rhs, *, cond=None, rhs_cache_size=100))]
    pub(crate) fn new(
        pattern: ConvertibleToPattern,
        rhs: ConvertibleToReplaceWith,
        cond: Option<ReplacementCondition>,
        rhs_cache_size: usize,
    ) -> PyResult<Self> {
        Ok(Self {
            rule: TensorRule::new(
                pattern.to_pattern()?.expr,
                rhs.to_replace_with()?,
                cond.map(|condition| condition.0),
                rhs_cache_size,
            )
            .map_err(TensorExpression::inference_error)?,
        })
    }
}
