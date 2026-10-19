use feynkit_graph::IntegralFamily;
use pyo3::prelude::*;
use symbolica::{api::python::PythonExpression, atom::Atom};

use crate::{error, graph::PyFeynmanDiagram, symbolic_kinematics::PyKinematics};

/// Ordered inverse propagators sharing loop momenta and external kinematics.
///
/// Denominators are expressions such as ``k.k - m2``, not reciprocals.
/// Integral powers follow the denominator order: zero omits a denominator,
/// while a negative power puts it in the numerator. Completing a family appends
/// auxiliary denominators without reordering its physical propagators.
///
/// Examples
/// --------
/// >>> from symbolica import S
/// >>> from symbolica.community import feynkit as fk
/// >>> k, m2 = S("k", "m2")
/// >>> kin = fk.Kinematics(momenta=[k])
/// >>> family = fk.IntegralFamily([k], [], [kin.scalar_product(k, k) - m2], kinematics=kin)
/// >>> assert family.is_complete and family.denominators[0] == kin.scalar_product(k, k) - m2
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(
    name = "IntegralFamily",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyIntegralFamily {
    pub(crate) inner: IntegralFamily,
}

impl PyIntegralFamily {
    /// Borrow the native family without serializing its symbolic expressions.
    pub fn as_family(&self) -> &IntegralFamily {
        &self.inner
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyIntegralFamily {
    /// Construct a family and determine the exact rank of its denominator basis.
    ///
    /// Examples
    /// --------
    /// Using the class example:
    ///
    /// >>> assert family.rank == 1
    ///
    /// Parameters
    /// ----------
    /// loop_momenta : list[Expression]
    ///     Distinct integrated momentum names, in their intended order.
    /// external_momenta : list[Expression]
    ///     Independent external momentum names.
    /// denominators : list[Expression]
    ///     Ordered inverse propagators.
    /// kinematics : Kinematics or None
    ///     Scalar-product assumptions and dimension; defaults to unconstrained 4D.
    #[new]
    #[pyo3(signature = (loop_momenta, external_momenta, denominators, *, kinematics=None))]
    fn new(
        loop_momenta: Vec<PythonExpression>,
        external_momenta: Vec<PythonExpression>,
        denominators: Vec<PythonExpression>,
        kinematics: Option<&PyKinematics>,
    ) -> PyResult<Self> {
        let default = feynkit_kinematics::Kinematics::new();
        Ok(Self {
            inner: IntegralFamily::new(
                loop_momenta.into_iter().map(|p| p.expr).collect(),
                external_momenta.into_iter().map(|p| p.expr).collect(),
                denominators.into_iter().map(|d| d.expr).collect(),
                kinematics.map(|context| &context.inner).unwrap_or(&default),
            )
            .map_err(error::integral_family)?,
        })
    }

    /// Complete the diagram's routed physical propagators with auxiliary slots.
    ///
    /// Physical denominators retain edge order and model masses. Preferred
    /// scalar products are appended if independent, then remaining directions
    /// are filled automatically. Tree graphs and dependent propagators fail.
    ///
    /// Examples
    /// --------
    /// >>> family = fk.IntegralFamily.from_diagram(diagram)
    /// >>> assert family.is_complete
    ///
    /// Parameters
    /// ----------
    /// diagram : FeynmanDiagram
    ///     Graph with its stored loop-momentum routing.
    /// independent_dot_products : list[Expression] or None
    ///     Preferred auxiliary denominators, tried in the supplied order.
    /// kinematics : Kinematics or None
    ///     Assumptions and dimension; defaults to the graph's symbolic dimension.
    #[staticmethod]
    #[pyo3(signature = (diagram, independent_dot_products=None, *, kinematics=None))]
    pub(crate) fn from_diagram(
        diagram: &PyFeynmanDiagram,
        independent_dot_products: Option<Vec<PythonExpression>>,
        kinematics: Option<&PyKinematics>,
    ) -> PyResult<Self> {
        let default = feynkit_kinematics::Kinematics::in_dimension(&Atom::var(
            feynkit_graph::symbols::dimension(),
        ))
        .expect("the shared dimension is a Symbolica symbol");
        let candidates = independent_dot_products
            .unwrap_or_default()
            .into_iter()
            .map(|product| product.expr)
            .collect::<Vec<_>>();
        Ok(Self {
            inner: IntegralFamily::from_diagram(
                diagram.as_diagram(),
                kinematics.map(|context| &context.inner).unwrap_or(&default),
                &candidates,
            )
            .map_err(error::diagram)?,
        })
    }

    /// Integrated momentum names in their original order.
    #[getter]
    fn loop_momenta(&self) -> Vec<PythonExpression> {
        self.inner
            .loop_momenta()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Independent external momentum names in their original order.
    #[getter]
    fn external_momenta(&self) -> Vec<PythonExpression> {
        self.inner
            .external_momenta()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Independent scalar-product coordinates involving at least one loop.
    #[getter]
    fn scalar_products(&self) -> Vec<PythonExpression> {
        self.inner
            .scalar_products()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Inverse propagators in the order used by integral powers.
    #[getter]
    fn denominators(&self) -> Vec<PythonExpression> {
        self.inner
            .denominators()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Immutable scalar-product assumptions and Lorentz dimension.
    #[getter]
    fn kinematics(&self) -> PyKinematics {
        PyKinematics {
            inner: self.inner.kinematics().clone(),
        }
    }

    /// Rank of the denominators over external invariants.
    #[getter]
    fn rank(&self) -> usize {
        self.inner.rank()
    }

    /// Whether the denominator span includes every loop scalar product.
    #[getter]
    fn is_complete(&self) -> bool {
        self.inner.is_complete()
    }

    /// Whether the denominators are linearly independent.
    #[getter]
    fn is_independent(&self) -> bool {
        self.inner.is_independent()
    }

    /// Append independent auxiliary denominators without changing the original slots.
    ///
    /// Examples
    /// --------
    /// >>> completed = family.complete()
    /// >>> assert completed.is_complete
    ///
    /// Parameters
    /// ----------
    /// candidates : list[Expression] or None
    ///     Preferred scalar products; automatic products fill the remaining span.
    #[pyo3(signature = (candidates=None))]
    fn complete(&self, candidates: Option<Vec<PythonExpression>>) -> PyResult<Self> {
        Ok(Self {
            inner: self
                .inner
                .complete(
                    &candidates
                        .unwrap_or_default()
                        .into_iter()
                        .map(|d| d.expr)
                        .collect::<Vec<_>>(),
                )
                .map_err(error::integral_family)?,
        })
    }

    /// Express a scalar numerator using labels for this complete denominator basis.
    ///
    /// Examples
    /// --------
    /// Using the class example:
    ///
    /// >>> d = S("d")
    /// >>> assert family.rewrite_numerator(kin.scalar_product(k, k), [d]) == d + m2
    ///
    /// Parameters
    /// ----------
    /// numerator : Expression
    ///     Scalar expression to rewrite with the family's assumptions.
    /// labels : list[Expression]
    ///     One distinct symbol or labeled call per denominator.
    fn rewrite_numerator(
        &self,
        numerator: &PythonExpression,
        labels: Vec<PythonExpression>,
    ) -> PyResult<PythonExpression> {
        self.inner
            .rewrite_numerator(
                &numerator.expr,
                &labels
                    .into_iter()
                    .map(|label| label.expr)
                    .collect::<Vec<_>>(),
            )
            .map(Into::into)
            .map_err(error::integral_family)
    }
}
