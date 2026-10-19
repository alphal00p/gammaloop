use pyo3::prelude::*;
use symbolica::api::python::PythonExpression;

use crate::error;

/// Immutable symbolic scalar-product assumptions in a declared Lorentz dimension.
///
/// Momentum inputs are unindexed names. Declare them before using linear
/// combinations so that other symbols are interpreted as scalar coefficients.
///
/// Examples
/// --------
/// >>> from symbolica import S
/// >>> from symbolica.community import feynkit as fk
/// >>> k, p, s = S("k", "p", "s")
/// >>> kin = fk.Kinematics(momenta=[k, p]).with_scalar_product(p, p, s)
/// >>> assert kin.scalar_product(2*p, p) == 2*s
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(
    name = "Kinematics",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyKinematics {
    pub(crate) inner: feynkit_kinematics::Kinematics,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyKinematics {
    /// Begin without scalar-product assumptions in the requested dimension.
    ///
    /// Examples
    /// --------
    /// >>> kin = fk.Kinematics(S("D"), momenta=[k, p])
    ///
    /// Parameters
    /// ----------
    /// dimension : Expression or None
    ///     Positive integer or symbolic dimension; defaults to four.
    /// momenta : list[Expression] or None
    ///     Momentum names used in scalar products and linear combinations.
    #[new]
    #[pyo3(signature = (dimension=None, *, momenta=None))]
    fn new(
        dimension: Option<&PythonExpression>,
        momenta: Option<Vec<PythonExpression>>,
    ) -> PyResult<Self> {
        let inner = match dimension {
            None => feynkit_kinematics::Kinematics::new(),
            Some(dimension) => feynkit_kinematics::Kinematics::in_dimension(&dimension.expr)
                .map_err(error::kinematics)?,
        };
        Ok(Self {
            inner: inner
                .with_momenta(
                    momenta
                        .unwrap_or_default()
                        .into_iter()
                        .map(|momentum| momentum.expr),
                )
                .map_err(error::kinematics)?,
        })
    }

    /// Lorentz dimension as an exact integer or symbol.
    #[getter]
    fn dimension(&self) -> PythonExpression {
        self.inner.dimension().to_symbolic().into()
    }

    /// Return a context with one additional scalar-product assumption.
    ///
    /// Examples
    /// --------
    /// >>> massive = kin.with_scalar_product(p, p, s)
    /// >>> assert massive.scalar_product(p, p) == s
    ///
    /// Parameters
    /// ----------
    /// left : Expression
    ///     First unindexed momentum name.
    /// right : Expression
    ///     Second unindexed momentum name.
    /// value : Expression
    ///     Exact scalar product to assign.
    fn with_scalar_product(
        &self,
        left: &PythonExpression,
        right: &PythonExpression,
        value: &PythonExpression,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: self
                .inner
                .clone()
                .with_scalar_product(&left.expr, &right.expr, value.expr.clone())
                .map_err(error::kinematics)?,
        })
    }

    /// Expand a bilinear scalar product and apply this context's assumptions.
    ///
    /// Nonlinear momentum expressions fail with KinematicsError.
    ///
    /// Examples
    /// --------
    /// >>> assert kin.scalar_product(2*p, p) == 2*s
    ///
    /// Parameters
    /// ----------
    /// left : Expression
    ///     First momentum or linear combination of declared momenta.
    /// right : Expression
    ///     Second momentum or linear combination of declared momenta.
    fn scalar_product(
        &self,
        left: &PythonExpression,
        right: &PythonExpression,
    ) -> PyResult<PythonExpression> {
        self.inner
            .scalar_product(&left.expr, &right.expr)
            .map(Into::into)
            .map_err(error::kinematics)
    }

    /// Apply the declared scalar-product assumptions to an expression.
    ///
    /// Examples
    /// --------
    /// >>> assert kin.apply(kin.scalar_product(p, p)) == s
    ///
    /// Parameters
    /// ----------
    /// expression : Expression
    ///     Scalar expression in native FeynKit scalar-product notation.
    fn apply(&self, expression: &PythonExpression) -> PythonExpression {
        self.inner.apply(&expression.expr).into()
    }
}
