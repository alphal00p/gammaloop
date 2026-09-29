use idenso::tensor::TensorRule;
use pyo3::{prelude::*, types::PyCFunction};
use symbolica::{
    api::python::{ConvertibleToPattern, ConvertibleToReplaceWith, PythonExpression},
    atom::AtomCore,
};

use crate::{ModuleInit, expression::TensorExpression, network::ReplacementCondition};

/// Create a reusable rule for replacing tensors or scalar subexpressions.
/// Index tensor ports before using them in a rule. Apply the rule with `.replace()`;
/// its replacement must keep the matched tensor interface, while a scalar zero
/// can inherit that interface. Callbacks receive a dictionary of wildcard bindings.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> r = sp.Representation.euc(2)
/// >>> A = sp.TensorName("docs::A")(r, r)
/// >>> A = A("i", "j")
/// >>> rule = sp.TensorRule(A, 2 * A)
/// >>> result = A.replace(rule)
///
/// Parameters
/// ----------
/// pattern : TensorExpression | Expression or scalar | HeldExpression
///     Expression or held expression to match, including Symbolica wildcards where needed.
/// rhs : callable
///     Replacement expression or callback receiving a dict of wildcard bindings. A callback may return a TensorExpression or a Symbolica scalar expression.
/// cond : PatternRestriction | Condition, optional
///     Optional Symbolica condition or pattern restriction. Only candidate matches satisfying it are replaced.
/// rhs_cache_size : int
///     Maximum cached right-hand-side callback results within one application. Use 0 for callbacks with side effects.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(frozen, name = "TensorRule", module = "symbolica.community.tensor")]
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
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyTensorRule {
    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> rule = sp.TensorRule(A("i", "j"), 2 * A("i", "j"))
    /// >>> text = repr(rule)
    fn __repr__(&self) -> String {
        format!("{:?}", self.rule)
    }

    /// Create a reusable rule for replacing tensors or scalar subexpressions.
    /// Index tensor ports before using them in a rule. Apply the rule with `.replace()`;
    /// its replacement must keep the matched tensor interface, while a scalar zero
    /// can inherit that interface. Callbacks receive a dictionary of wildcard bindings.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> A = A("i", "j")
    /// >>> rule = sp.TensorRule(A, 2 * A)
    /// >>> result = A.replace(rule)
    ///
    /// Parameters
    /// ----------
    /// pattern : TensorExpression | Expression or scalar | HeldExpression
    ///     Expression or held expression to match, including Symbolica wildcards where needed.
    /// rhs : callable
    ///     Replacement expression or callback receiving a dict of wildcard bindings. A callback may return a TensorExpression or a Symbolica scalar expression.
    /// cond : PatternRestriction | Condition, optional
    ///     Optional Symbolica condition or pattern restriction. Only candidate matches satisfying it are replaced.
    /// rhs_cache_size : int
    ///     Maximum cached right-hand-side callback results within one application. Use 0 for callbacks with side effects.
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
