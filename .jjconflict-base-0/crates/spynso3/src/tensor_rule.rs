use idenso::tensor::TensorRule;
use pyo3::{prelude::*, types::PyCFunction};
use symbolica::{
    api::python::{ConvertibleToPattern, ConvertibleToReplaceWith, PythonExpression},
    atom::AtomCore,
};

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

impl PyTensorRule {
    pub(crate) fn extract_rules(value: &Bound<'_, PyAny>) -> PyResult<Vec<TensorRule<'static>>> {
        if let Ok(rule) = value.extract::<PyRef<'_, Self>>() {
            return Ok(vec![rule.rule.clone()]);
        }
        value.extract::<Vec<Py<Self>>>().map(|rules| {
            rules
                .iter()
                .map(|rule| rule.borrow(value.py()).rule.clone())
                .collect()
        })
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyTensorRule {
    fn __repr__(&self) -> String {
        format!("{:?}", self.rule)
    }

    #[new]
    #[pyo3(signature = (pattern, rhs, *, cond=None, rhs_cache_size=100))]
    pub(crate) fn new(
        py: Python<'_>,
        pattern: &Bound<'_, PyAny>,
        rhs: &Bound<'_, PyAny>,
        cond: Option<ReplacementCondition>,
        rhs_cache_size: usize,
    ) -> PyResult<Self> {
        let pattern = if let Ok(value) = pattern.extract::<PyRef<'_, TensorExpression>>() {
            value.atom().to_pattern()
        } else {
            pattern
                .extract::<ConvertibleToPattern>()?
                .to_pattern()?
                .expr
        };
        let rhs = if let Ok(value) = rhs.extract::<PyRef<'_, TensorExpression>>() {
            value.atom().clone().into()
        } else {
            match rhs.extract::<ConvertibleToReplaceWith>()? {
                ConvertibleToReplaceWith::Pattern(pattern) => {
                    ConvertibleToReplaceWith::Pattern(pattern).to_replace_with()?
                }
                ConvertibleToReplaceWith::Map(callback) => {
                    // Keep Symbolica's match-map and error policy; adapt only
                    // the typed return value at this explicit binding boundary.
                    let callback = PyCFunction::new_closure(
                        py,
                        None,
                        None,
                        move |args, _| -> PyResult<Py<PyAny>> {
                            let py = args.py();
                            let result = callback.call1(py, args)?;
                            if let Ok(value) = result.extract::<PyRef<'_, TensorExpression>>(py) {
                                Ok(PythonExpression {
                                    expr: value.atom().clone(),
                                }
                                .into_pyobject(py)?
                                .unbind()
                                .into_any())
                            } else {
                                Ok(result)
                            }
                        },
                    )?;
                    ConvertibleToReplaceWith::Map(callback.unbind().into_any()).to_replace_with()?
                }
            }
        };
        Ok(Self {
            rule: TensorRule::new(
                pattern,
                rhs,
                cond.map(|condition| condition.0),
                rhs_cache_size,
            )
            .map_err(TensorExpression::inference_error)?,
        })
    }
}
