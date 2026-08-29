use feynkit_tensor::TensorReducer;
use pyo3::{
    prelude::*,
    types::{PyAny, PyModule},
};
use symbolica::api::python::PythonExpression;

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen_derive::{gen_stub_pyclass, gen_stub_pymethods};

use crate::error;

/// Project Lorentz-tensor integrands onto Spenso invariants.
///
/// The reducer implements the symmetry-orbit form of the orthogonal
/// Weingarten projector. Repeated loop and projector momenta are kept in
/// compact contraction classes, making the common rank-20 vacuum projections
/// practical without constructing the full ``19!!`` pairing matrix.
/// Fully contracted projectors yield scalar ``spenso::dot`` invariants. If
/// projector indices remain free, the returned expression retains them as
/// explicit ``spenso::g`` tensors. Keep mixed high-rank reductions symbolic
/// in ``D`` until afterward: at fixed positive integer dimension below half
/// the rank, dimension-specific identities make the universal metric basis
/// singular. The all-equal isotropic fast path remains well defined.
///
/// Examples
/// --------
/// Select every native FeynKit momentum in a pure vacuum numerator:
///
/// >>> from symbolica import E
/// >>> reducer = fk.TensorReducer.feynkit(E("4"))
/// >>> scalar_numerator = reducer.reduce(vacuum_numerator)
/// >>> scalar_numerator
///
/// Parameters
/// ----------
/// dimension : Expression
///     Lorentz-space dimension, commonly ``D`` or ``4 - 2*eps``.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "TensorReducer",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyTensorReducer {
    pub(crate) inner: TensorReducer,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyTensorReducer {
    /// Construct a reducer without selecting an integrated momentum.
    ///
    /// Add selectors with :meth:`with_integrated_head` or
    /// :meth:`with_integrated_vector`. The object is immutable; each selector
    /// returns a new reducer.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import E
    /// >>> reducer = fk.TensorReducer(E("4 - 2*eps")).with_integrated_head("Loop::k")
    ///
    /// Parameters
    /// ----------
    /// dimension : Expression
    ///     Lorentz-space dimension used by every ``spenso::mink`` slot.
    #[new]
    fn new(dimension: &PythonExpression) -> Self {
        Self {
            inner: TensorReducer::new(dimension.expr.clone()),
        }
    }

    /// Construct a reducer that selects every ``FeynKit::Momentum`` tensor.
    ///
    /// This is the convenient constructor for vacuum numerators produced by
    /// FeynKit's native Feynman-rule generator. It selects the entire
    /// ``FeynKit::Momentum`` head and is therefore intended for pure vacuum
    /// numerators, where every such momentum is integrated. If a graph still
    /// contains external ``FeynKit::Momentum`` tensors, construct a reducer
    /// with exact :meth:`with_integrated_vector` selectors for its internal
    /// momenta instead.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import E
    /// >>> reducer = fk.TensorReducer.feynkit(E("4"))
    /// >>> reduced = reducer.reduce(vacuum_diagram.numerator_expression())
    ///
    /// Parameters
    /// ----------
    /// dimension : Expression
    ///     Lorentz-space dimension used by every ``spenso::mink`` slot. It
    ///     must match the dimension carried by the input expression.
    #[staticmethod]
    fn feynkit(dimension: &PythonExpression) -> Self {
        Self {
            inner: TensorReducer::feynkit(dimension.expr.clone()),
        }
    }

    /// Return the configured Lorentz-space dimension.
    #[getter]
    fn dimension(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.dimension().clone(),
        }
    }

    /// Select all rank-one tensors with the qualified Symbolica head name.
    ///
    /// Examples
    /// --------
    /// >>> reducer = fk.TensorReducer(D).with_integrated_head("Loop::k")
    /// >>> scalar = reducer.reduce(numerator)
    ///
    /// Parameters
    /// ----------
    /// head : str
    ///     Qualified function name of the integrated vector.
    fn with_integrated_head(&self, head: &str) -> Self {
        Self {
            inner: self.inner.clone().with_integrated_head_name(head),
        }
    }

    /// Select one exact compact Spenso vector.
    ///
    /// A compact vector has a representation but no explicit index, for
    /// example ``K(1, spenso::mink(D))``. This distinguishes different loop
    /// momenta that share one tensor head.
    ///
    /// Examples
    /// --------
    /// >>> k1 = K(1, mink(D))
    /// >>> reducer = fk.TensorReducer(D).with_integrated_vector(k1)
    /// >>> scalar = reducer.reduce(numerator)
    ///
    /// Parameters
    /// ----------
    /// vector : Expression
    ///     Exact indexed-free vector to integrate.
    fn with_integrated_vector(&self, vector: &PythonExpression) -> Self {
        Self {
            inner: self
                .inner
                .clone()
                .with_integrated_vector(vector.expr.clone()),
        }
    }

    /// Set the labeled-pairing budget for unsymmetrized or free-index output.
    ///
    /// Symmetric contraction-orbit paths do not consume this budget.
    ///
    /// Examples
    /// --------
    /// >>> reducer = reducer.with_pairing_limit(200_000)
    ///
    /// Parameters
    /// ----------
    /// limit : int
    ///     Maximum number of labeled perfect matchings to enumerate.
    fn with_pairing_limit(&self, limit: usize) -> Self {
        Self {
            inner: self.inner.clone().with_pairing_limit(limit as u128),
        }
    }

    /// Set the relative-pairing budget for residual free-index output.
    ///
    /// Examples
    /// --------
    /// >>> reducer = reducer.with_pairing_product_limit(150_000_000)
    ///
    /// Parameters
    /// ----------
    /// limit : int
    ///     Maximum Cartesian product of internal and projector matchings.
    fn with_pairing_product_limit(&self, limit: usize) -> Self {
        Self {
            inner: self.inner.clone().with_pairing_product_limit(limit as u128),
        }
    }

    /// Set the maximum number of distinct invariant output terms.
    ///
    /// Examples
    /// --------
    /// >>> reducer = reducer.with_output_term_limit(20_000)
    ///
    /// Parameters
    /// ----------
    /// limit : int
    ///     Maximum compact contraction classes to materialize.
    fn with_output_term_limit(&self, limit: usize) -> Self {
        Self {
            inner: self.inner.clone().with_output_term_limit(limit),
        }
    }

    /// Reduce a Spenso tensor expression to one Symbolica expression.
    ///
    /// Rank-one tensors must carry a final ``spenso::mink(D,index)`` argument.
    /// Fully contracted projectors are returned as scalar ``spenso::dot``
    /// invariants. Residual free projector pairs remain explicit
    /// ``spenso::g`` tensors; this method intentionally does not reject
    /// tensor-valued output. Odd-rank vacuum tensors vanish.
    ///
    /// Examples
    /// --------
    /// A rank-two vacuum projection becomes a product of dot products divided
    /// by the dimension:
    ///
    /// >>> reduced = reducer.reduce(k(mu) * k(nu) * p(mu) * p(nu))
    /// >>> reduced
    ///
    /// Parameters
    /// ----------
    /// expression : Expression
    ///     Tensor numerator or projected tensor numerator to reduce.
    fn reduce(&self, expression: &PythonExpression) -> PyResult<PythonExpression> {
        self.inner
            .reduce(expression.expr.as_view())
            .map(|reduction| PythonExpression {
                expr: reduction.into_expression(),
            })
            .map_err(error::tensor)
    }

    /// Return a concise description of the reducer configuration.
    ///
    /// Examples
    /// --------
    /// >>> print(reducer)
    ///
    fn __repr__(&self) -> String {
        format!("TensorReducer(dimension={})", self.inner.dimension())
    }

    /// Write the reducer summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// IPython calls this method when formatting a reducer for text display.
    ///
    /// Parameters
    /// ----------
    /// pretty : Any
    ///     The IPython pretty-printer object.
    /// cycle : bool
    ///     Whether this object is part of a recursive formatting cycle.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1(
            "text",
            (if cycle {
                "...".to_owned()
            } else {
                self.__repr__()
            },),
        )?;
        Ok(())
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyTensorReducer>()?;
    Ok(())
}
