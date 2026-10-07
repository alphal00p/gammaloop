use feynkit_tensor::TensorReducer;
use pyo3::{
    prelude::*,
    types::{PyAny, PyModule},
};
use spynso3::structure::SpensoName;
use symbolica::{api::python::PythonExpression, atom::AtomView};

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
/// For denominators depending on external momenta, add an independent basis
/// with the ``external`` constructor keyword. Only the transverse components are then
/// rotationally averaged; odd total ranks need not vanish.
///
/// Examples
/// --------
/// Average a rank-two vacuum numerator over the loop direction. The result
/// is proportional to the metric tensor times ``k.k / D``.
///
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hepkit as hep
/// >>> from symbolica.community.tensor import TensorName, Representation, PortPattern, TensorPattern
/// >>> D, mu, nu = S("D", "mu", "nu")
/// >>> k, p = (TensorName.vector("hep_reducer_docs::" + name).to_expression() for name in ("k", "p"))
/// >>> lorentz = Representation.mink(D)
/// >>> reducer = hep.TensorReducer(D, integrated=[k(PortPattern.exact(lorentz))])
/// >>> numerator = k(PortPattern.exact(lorentz, mu)) * k(PortPattern.exact(lorentz, nu))
/// >>> projected = reducer.reduce(numerator)
///
/// Parameters
/// ----------
/// dimension : Expression
///     Lorentz-space dimension, commonly ``D`` or ``4 - 2*eps``.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "TensorReducer",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyTensorReducer {
    // CPython's wasm allocator guarantees only 8-byte alignment, while the
    // reducer's u128 pairing limits require 16. Keep that payload on the Rust heap.
    pub(crate) inner: Box<TensorReducer>,
}

#[cfg(target_arch = "wasm32")]
const _: () = {
    assert!(std::mem::align_of::<PyTensorReducer>() <= 8);
    assert!(
        std::mem::align_of::<<PyTensorReducer as pyo3::impl_::pyclass::PyClassImpl>::Layout>() <= 8
    );
};

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_tensor_usage, on_success)]
#[pymethods]
impl PyTensorReducer {
    /// Configure the integrated momenta and independent external basis.
    ///
    /// Bare symbols and ``TensorName`` entries in ``integrated`` select every
    /// vector with that head. Compact vector expressions select one exact
    /// momentum, including its scalar arguments. External entries are compact
    /// vectors and take precedence over integrated-head selections.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``TensorReducer`` class example:
    ///
    /// >>> reducer = hep.TensorReducer(D, integrated=[k])
    /// >>> reducer = hep.TensorReducer(D, integrated=[k(PortPattern.exact(lorentz))], external=[p(PortPattern.exact(lorentz))])
    ///
    /// Parameters
    /// ----------
    /// dimension : Expression
    ///     Lorentz-space dimension used by every ``spenso::mink`` slot.
    /// integrated : sequence of TensorName or Expression or TensorExpression, optional
    ///     Whole vector heads or exact compact vectors to integrate.
    /// external : sequence of Expression or TensorExpression, optional
    ///     Independent compact vectors in the denominator's external basis.
    #[new]
    #[pyo3(signature = (dimension, *, integrated=Vec::new(), external=Vec::new()))]
    fn new(
        dimension: &PythonExpression,
        #[gen_stub(override_type(type_repr="typing.Sequence[symbolica.community.tensor.TensorName | symbolica.Expression | symbolica.community.tensor.TensorExpression]", imports=("typing", "symbolica", "symbolica.community.tensor")))]
        integrated: Vec<Bound<'_, PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Sequence[symbolica.Expression | symbolica.community.tensor.TensorExpression]", imports=("typing", "symbolica", "symbolica.community.tensor")))]
        external: Vec<Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let mut inner = TensorReducer::new(dimension.expr.clone());
        for selector in integrated {
            if let Ok(name) = selector.extract::<SpensoName>() {
                inner = inner.with_integrated_head(name.name);
                continue;
            }
            // TensorExpression inherits Expression and shares its immutable atom.
            let expression = selector.extract::<PyRef<'_, PythonExpression>>()?;
            inner = match expression.expr.as_view() {
                AtomView::Var(head) => inner.with_integrated_head(head.get_symbol()),
                AtomView::Fun(_) => inner.with_integrated_vector(expression.expr.clone()),
                _ => {
                    return Err(pyo3::exceptions::PyValueError::new_err(
                        "integrated selectors must be bare symbols or compact vector expressions",
                    ));
                }
            };
        }
        for vector in external {
            let expression = vector.extract::<PyRef<'_, PythonExpression>>()?;
            inner = inner.with_external_vector(expression.expr.clone());
        }
        Ok(Self {
            inner: Box::new(inner),
        })
    }

    /// Construct a reducer that selects every ``gammalooprs::Q`` tensor.
    ///
    /// This is the convenient constructor for vacuum numerators produced by
    /// HEP native Feynman-rule generator. It selects the entire
    /// ``gammalooprs::Q`` head and is therefore intended for pure vacuum
    /// numerators, where every such momentum is integrated. If a graph still
    /// contains external ``gammalooprs::Q`` tensors, construct a reducer
    /// with exact compact vectors in ``integrated`` for its internal
    /// momenta instead.
    ///
    /// Examples
    /// --------
    /// The equivalent explicit selector shows which momentum head is integrated:
    ///
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hepkit as hep
    /// >>> model = hep.Model.phi4()
    /// >>> vacuum_diagram = model.process([], []).generate_diagrams(loops=2, factorized_loop_topologies_count_range=None).diagrams[0]
    /// >>> reducer = hep.TensorReducer(E("4"), integrated=[E("gammalooprs::Q")])
    /// >>> reduced = reducer.reduce(vacuum_diagram.numerator_expression().to_expression())
    ///
    /// Parameters
    /// ----------
    /// dimension : Expression
    ///     Lorentz-space dimension used by every ``spenso::mink`` slot. It
    ///     must match the dimension carried by the input expression.
    #[staticmethod]
    fn feynkit(dimension: &PythonExpression) -> Self {
        Self {
            inner: Box::new(TensorReducer::feynkit(dimension.expr.clone())),
        }
    }

    /// Return the configured Lorentz-space dimension.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``TensorReducer`` class example:
    ///
    /// >>> assert reducer.dimension == D
    #[getter]
    fn dimension(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.dimension().clone(),
        }
    }

    /// Set the labeled-pairing budget for unsymmetrized or free-index output.
    ///
    /// Symmetric contraction-orbit paths do not consume this budget.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``TensorReducer`` class example:
    ///
    /// >>> reducer = reducer.with_pairing_limit(200_000)
    ///
    /// Parameters
    /// ----------
    /// limit : int
    ///     Maximum number of labeled perfect matchings to enumerate.
    fn with_pairing_limit(&self, limit: usize) -> Self {
        Self {
            inner: Box::new(
                self.inner
                    .as_ref()
                    .clone()
                    .with_pairing_limit(limit as u128),
            ),
        }
    }

    /// Set the relative-pairing budget for residual free-index output.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``TensorReducer`` class example:
    ///
    /// >>> reducer = reducer.with_pairing_product_limit(150_000_000)
    ///
    /// Parameters
    /// ----------
    /// limit : int
    ///     Maximum Cartesian product of internal and projector matchings.
    fn with_pairing_product_limit(&self, limit: usize) -> Self {
        Self {
            inner: Box::new(
                self.inner
                    .as_ref()
                    .clone()
                    .with_pairing_product_limit(limit as u128),
            ),
        }
    }

    /// Set the maximum number of distinct invariant output terms.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``TensorReducer`` class example:
    ///
    /// >>> reducer = reducer.with_output_term_limit(20_000)
    ///
    /// Parameters
    /// ----------
    /// limit : int
    ///     Maximum compact contraction classes to materialize.
    fn with_output_term_limit(&self, limit: usize) -> Self {
        Self {
            inner: Box::new(self.inner.as_ref().clone().with_output_term_limit(limit)),
        }
    }

    /// Reduce a Spenso tensor expression to one Symbolica expression.
    ///
    /// Rank-one tensors must carry a final ``spenso::mink(D,index)`` argument.
    /// Fully contracted projectors are returned as scalar ``spenso::dot``
    /// invariants. Residual free projector pairs remain explicit
    /// ``spenso::g`` tensors; this method intentionally does not reject
    /// tensor-valued output. Compact dots between an integrated vector and a
    /// spectator are reduced directly, including nonnegative integer powers.
    /// Dots between integrated vectors or with declared external basis vectors
    /// remain scalar invariants. Negative or noninteger powers of spectator
    /// contractions are rejected. Odd-rank vacuum tensors vanish.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``TensorReducer`` class example:
    ///
    /// A rank-two vacuum projection becomes a product of dot products divided
    /// by the dimension:
    ///
    /// >>> from symbolica import S
    /// >>> from symbolica.community import hepkit as hep
    /// >>> D, mu, nu = S("hep_docs::D", "hep_docs::mu", "hep_docs::nu")
    /// >>> from symbolica.community.tensor import TensorName, Representation, PortPattern, TensorPattern
    /// >>> k, p = (TensorName.vector("hep_reducer_docs::" + name).to_expression() for name in ("k", "p"))
    /// >>> lorentz, dot = Representation.mink(D), TensorPattern.dot
    /// >>> k_compact = k(PortPattern.exact(lorentz))
    /// >>> p_compact = p(PortPattern.exact(lorentz))
    /// >>> numerator = (
    /// ...     k(PortPattern.exact(lorentz, mu)) * k(PortPattern.exact(lorentz, nu))
    /// ...     * p(PortPattern.exact(lorentz, mu)) * p(PortPattern.exact(lorentz, nu))
    /// ... )
    /// >>> reducer = hep.TensorReducer(D, integrated=[k_compact])
    /// >>> reduced = reducer.reduce(numerator)
    /// >>> expected = dot(k_compact, k_compact) * dot(p_compact, p_compact) / D
    /// >>> assert reduced == expected
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
    /// Using the setup in the ``TensorReducer`` class example:
    ///
    /// >>> print(reducer)
    fn __repr__(&self) -> String {
        format!("TensorReducer(dimension={})", self.inner.dimension())
    }

    /// Write the reducer summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``TensorReducer`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(reducer)
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

#[cfg(test)]
mod tests {
    use super::PyTensorReducer;

    #[test]
    fn python_wrapper_fits_wasm_allocator_alignment() {
        assert!(std::mem::align_of::<PyTensorReducer>() <= 8);
        assert!(
            std::mem::align_of::<<PyTensorReducer as pyo3::impl_::pyclass::PyClassImpl>::Layout>()
                <= 8
        );
    }
}
