use idenso::tensor::TensorRule;
use pyo3::{prelude::*, types::PyCFunction};
use symbolica::{
    api::python::{ConvertibleToPattern, ConvertibleToReplaceWith, PythonExpression},
    atom::AtomCore,
};

use crate::{ModuleInit, expression::TensorExpression, network::ReplacementCondition};

/// A reusable replacement rule that preserves tensor interfaces.
///
/// Use with TensorExpression.replace.
/// A replacement must preserve the matched tensor's external indices, including
/// their representations. Literal replacements are checked where possible at
/// construction; wildcard-dependent replacements are checked when matched.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
/// >>> space = Representation.euc(2)
/// >>> A = TensorName("M")(space, space)
/// >>> from symbolica.community.tensor import TensorRule
/// >>> indexed = A("i", "j")
/// >>> rule = TensorRule(indexed, 2 * indexed)
/// >>> indexed.replace(rule).rank
/// 2
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
#[spenso_macros::track_usage(crate::record_usage, on_success)]
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

    /// Construct a checked whole-tensor replacement rule.
    ///
    /// Parameters
    /// ----------
    /// pattern : TensorExpression, scalar expression, or HeldExpression
    ///     Tensor to match. TensorPattern can describe wildcard tensor syntax.
    /// rhs : TensorExpression, scalar expression, HeldExpression, or callable
    ///     Replacement with the same external interface. A callable receives
    ///     a dict from wildcard Expressions to matched Expressions and returns
    ///     the replacement. Zero retains the matched tensor's interface.
    /// cond : PatternRestriction or Condition, optional
    ///     Condition evaluated for each candidate match.
    /// rhs_cache_size : int, default 100
    ///     Maximum number of cached replacement results per application. Use
    ///     zero for a callback whose side effects must run at every match.
    ///
    /// Returns
    /// -------
    /// TensorRule
    ///     Reusable replacement, without applying it yet.
    ///
    /// Notes
    /// -----
    /// Right-hand-side wildcards must be bound by the pattern. Index multiplicity
    /// and external representations must remain valid after substitution.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import TensorRule
    /// >>> indexed = A("i", "j")
    /// >>> rule = TensorRule(indexed, 2 * indexed)
    /// >>> result = indexed.replace(rule)
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
