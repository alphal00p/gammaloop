use feynkit_kinematics::{
    Axis, Boost, ClusteringResult, FourMomentum, Helicity, Jet, JetAlgorithm, JetDefinition,
    Rotation, ThreeMomentum,
};
use pyo3::{
    exceptions::PyIndexError,
    prelude::*,
    types::{PyAny, PyList, PyModule},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    derive::{
        gen_methods_from_python, gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods,
    },
    inventory::submit,
};

use crate::error;
use spynso3::expression::TensorExpression;
use symbolica::api::python::{ConvertibleToExpression, PythonExpression};

/// Scoped symbolic scalar products and two-to-two Mandelstam kinematics.
///
/// The immutable object uses Spenso's dot products and metric
/// shorthand. Apply it after contracting tensors with Idenso. Momentum inputs
/// are unindexed names, and mass inputs are squared masses.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> from symbolica import S, E
/// >>> p1, p2, p3, p4, s, t, u = S("p1", "p2", "p3", "p4", "s", "t", "u")
/// >>> kin = hep.Kinematics.mandelstam([p1, p2, p3, p4], [E("0")]*4, [s, t, u])
/// >>> assert kin.scalar_product(p1, p2) == s/2
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
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

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyKinematics {
    /// Start with no scalar-product assumptions in the chosen dimension.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Kinematics`` class example:
    ///
    /// >>> kin = hep.Kinematics()
    /// >>> dimensional = hep.Kinematics(S("D"))
    ///
    /// Parameters
    /// ----------
    /// dimension : Expression | None
    ///     Integer or symbolic Lorentz dimension; defaults to four.
    /// momenta : list[Expression] | None
    ///     Momentum names used in linear combinations with scalar coefficients.
    #[new]
    #[pyo3(signature = (dimension=None, *, momenta=None))]
    fn new(
        dimension: Option<&PythonExpression>,
        momenta: Option<Vec<PythonExpression>>,
    ) -> PyResult<Self> {
        let inner = match dimension {
            None => feynkit_kinematics::Kinematics::new(),
            Some(dimension) => feynkit_kinematics::Kinematics::in_dimension(&dimension.expr)
                .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?,
        };
        let inner = inner
            .with_momenta(momenta.unwrap_or_default().into_iter().map(|p| p.expr))
            .map_err(error::kinematics)?;
        Ok(Self { inner })
    }

    /// Lorentz dimension as a Symbolica integer or symbol.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Kinematics`` class example:
    ///
    /// >>> assert hep.Kinematics(S("D")).dimension == S("D")
    #[getter]
    fn dimension(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.dimension().to_symbolic(),
        }
    }

    /// Set the invariants for ``p1 + p2 -> p3 + p4``.
    ///
    /// The convention is ``s=(p1+p2)^2``, ``t=(p1-p3)^2``, and
    /// ``u=(p1-p4)^2``. These obey ``s+t+u=sum(mass_squared)``; use Symbolica
    /// substitution when you want to eliminate one invariant.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Kinematics`` class example:
    ///
    /// >>> kin = hep.Kinematics.mandelstam([p1, p2, p3, p4], [E("0")]*4, [s, t, u])
    ///
    /// Parameters
    /// ----------
    /// momenta : list[Expression]
    ///     Four unindexed momenta, with the incoming pair first.
    /// mass_squared : list[Expression]
    ///     Four squared masses in the same order.
    /// invariants : list[Expression]
    ///     The three symbols or expressions ``s, t, u``.
    #[staticmethod]
    fn mandelstam(
        momenta: [PythonExpression; 4],
        mass_squared: [PythonExpression; 4],
        invariants: [PythonExpression; 3],
    ) -> PyResult<Self> {
        Ok(Self {
            inner: feynkit_kinematics::Kinematics::mandelstam(
                momenta.each_ref().map(|p| &p.expr),
                mass_squared.map(|m| m.expr),
                invariants.map(|s| s.expr),
            )
            .map_err(error::kinematics)?,
        })
    }

    /// Return a new context with one scalar product set.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Kinematics`` class example:
    ///
    /// >>> p, m = S("p", "m")
    /// >>> kin = hep.Kinematics(momenta=[p]).with_scalar_product(p, p, m**2)
    /// >>> assert kin.scalar_product(p, p) == m**2
    ///
    /// Parameters
    /// ----------
    /// left : Expression
    ///     First unindexed momentum.
    /// right : Expression
    ///     Second unindexed momentum.
    /// value : Expression
    ///     Assumed scalar product.
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

    /// Expand a bilinear scalar product and apply known assumptions.
    ///
    /// For linear combinations, declare the momentum names in the constructor
    /// or by setting scalar products. Other symbols are scalar coefficients.
    /// Nonlinear momentum expressions raise ``KinematicsError``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Kinematics`` class example:
    ///
    /// >>> assert kin.scalar_product(p1, p2) == s/2
    /// >>> assert kin.scalar_product(p1 + p2, p1 + p2) == s
    ///
    /// Parameters
    /// ----------
    /// left : Expression
    ///     First unindexed momentum or linear combination.
    /// right : Expression
    ///     Second unindexed momentum or linear combination.
    fn scalar_product(
        &self,
        left: &PythonExpression,
        right: &PythonExpression,
    ) -> PyResult<PythonExpression> {
        Ok(PythonExpression {
            expr: self
                .inner
                .scalar_product(&left.expr, &right.expr)
                .map_err(error::kinematics)?,
        })
    }

    /// Return the initial-state denominator for a cross section or decay rate.
    ///
    /// Two momenta give ``4*sqrt((p1.p2)**2-p1**2*p2**2)``. One momentum
    /// gives ``2*sqrt(p**2)`` for a decay in its rest frame. Divide the squared
    /// matrix element times phase space by this value. Inputs must be physical,
    /// future-directed on-shell momenta. Symbolica retains square-root branches;
    /// declare positive invariants with ``S("s", is_positive=True)`` when known.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Kinematics`` class example:
    ///
    /// >>> p1, p2, s = S("p1", "p2", "s")
    /// >>> kin = hep.Kinematics(momenta=[p1, p2])
    /// >>> kin = kin.with_scalar_product(p1, p1, E("0")).with_scalar_product(p2, p2, E("0"))
    /// >>> kin = kin.with_scalar_product(p1, p2, s/2)
    /// >>> flux = kin.flux(p1, p2)
    ///
    /// Parameters
    /// ----------
    /// first : Expression
    ///     Incoming unindexed momentum or declared linear combination.
    /// second : Expression | None
    ///     Other incoming momentum; None selects a rest-frame decay.
    #[pyo3(signature = (first, second=None))]
    fn flux(
        &self,
        first: &PythonExpression,
        second: Option<&PythonExpression>,
    ) -> PyResult<PythonExpression> {
        Ok(PythonExpression {
            expr: self
                .inner
                .flux(&first.expr, second.map(|p| &p.expr))
                .map_err(error::kinematics)?,
        })
    }

    /// Return four-dimensional two-body phase space per unit solid angle.
    ///
    /// This is ``dPhi_2/dOmega`` in the final pair's rest frame, with
    /// ``(2*pi)**4*delta**4(P-p1-p2)`` and ``d**3p/((2*pi)**3*2E)`` for each
    /// final particle. Use physical on-shell momenta above threshold. Flux,
    /// spin/color averages and identical-particle factors are separate. For an
    /// angle-independent amplitude, integrating this measure gives ``4*pi``
    /// times the returned expression. Non-four-dimensional contexts are rejected.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Kinematics`` class example:
    ///
    /// >>> p1, p2, p3, p4, s, t, u = S("p1", "p2", "p3", "p4", "s", "t", "u")
    /// >>> kin = hep.Kinematics.mandelstam([p1, p2, p3, p4], [E("0")]*4, [s, t, u])
    /// >>> density = kin.two_body_phase_space(p3, p4)
    ///
    /// Parameters
    /// ----------
    /// first : Expression
    ///     First outgoing unindexed momentum or declared linear combination.
    /// second : Expression
    ///     Second outgoing unindexed momentum or declared linear combination.
    fn two_body_phase_space(
        &self,
        first: &PythonExpression,
        second: &PythonExpression,
    ) -> PyResult<PythonExpression> {
        Ok(PythonExpression {
            expr: self
                .inner
                .two_body_phase_space(&first.expr, &second.expr)
                .map_err(error::kinematics)?,
        })
    }

    /// Return four-dimensional three-body phase space per two Dalitz invariants.
    ///
    /// This is ``dPhi_3/(ds12*ds23) = 1/(128*pi**3*P**2)``, where
    /// ``P=first+second+third`` and ``sij=(pi+pj)**2``. The overall spatial
    /// orientation is integrated. Use an orientation-independent or
    /// orientation-averaged squared amplitude and physical on-shell momenta.
    /// The measure includes ``(2*pi)**4*delta**4`` and one
    /// ``d**3p/((2*pi)**3*2E)`` for each final particle.
    ///
    /// Masses constrain the allowed Dalitz region; its boundaries are not
    /// imposed here. Flux, spin/color averages and identical-particle factors
    /// remain separate. Non-four-dimensional contexts are rejected.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Kinematics`` class example:
    ///
    /// >>> k1, k2, k3 = S("k1", "k2", "k3")
    /// >>> kin = hep.Kinematics(momenta=[k1, k2, k3])
    /// >>> density = kin.three_body_phase_space(k1, k2, k3)
    ///
    /// Parameters
    /// ----------
    /// first, second, third : Expression
    ///     Final-state on-shell momenta. Their scalar products determine the
    ///     total invariant mass squared through this kinematic context.
    fn three_body_phase_space(
        &self,
        first: &PythonExpression,
        second: &PythonExpression,
        third: &PythonExpression,
    ) -> PyResult<PythonExpression> {
        Ok(PythonExpression {
            expr: self
                .inner
                .three_body_phase_space(&first.expr, &second.expr, &third.expr)
                .map_err(error::kinematics)?,
        })
    }

    /// Substitute scalar products without mutating global assumptions.
    ///
    /// Accepts a scalar Expression or a tensor expression. Tensor results retain
    /// the ordered open slots, including when the result is zero.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Kinematics`` class example:
    ///
    /// >>> p1, p2, p3, p4, s, t, u = S("p1", "p2", "p3", "p4", "s", "t", "u")
    /// >>> free = hep.Kinematics(momenta=[p1, p2, p3, p4])
    /// >>> kin = hep.Kinematics.mandelstam([p1, p2, p3, p4], [E("0")]*4, [s, t, u])
    /// >>> contracted_expression = free.scalar_product(p1, p2)
    /// >>> assert kin.apply(contracted_expression) == s/2
    ///
    /// Parameters
    /// ----------
    /// expression : Expression
    ///     Expression with compact scalar products.
    #[gen_stub(skip)]
    fn apply(&self, py: Python<'_>, expression: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        let atom = expression
            .extract::<ConvertibleToExpression>()?
            .to_expression()
            .expr;
        let result = self.inner.apply(&atom);
        if expression.is_instance_of::<TensorExpression>() {
            let tensor = expression.extract::<PyRef<'_, TensorExpression>>()?;
            TensorExpression::preserving_interface(&tensor, py, result).map(Py::into_any)
        } else {
            Py::new(py, PythonExpression { expr: result }).map(Py::into_any)
        }
    }
}

#[cfg(feature = "python_stubgen")]
submit! {
    gen_methods_from_python! {
        r#"
        import typing

        class PyKinematics:
            @typing.overload
            def apply(
                self,
                expression: pyo3_stub_gen.RustType["TensorExpression"],
            ) -> pyo3_stub_gen.RustType["TensorExpression"]:
                """
                Substitute scalar products without mutating global assumptions.

                Accepts a scalar Expression or a tensor expression. Tensor results retain
                the ordered open slots, including when the result is zero.

                Examples
                --------
                Using the setup in the ``Kinematics`` class example:

                >>> p1, p2, p3, p4, s, t, u = S("p1", "p2", "p3", "p4", "s", "t", "u")
                >>> free = hep.Kinematics(momenta=[p1, p2, p3, p4])
                >>> kin = hep.Kinematics.mandelstam([p1, p2, p3, p4], [E("0")]*4, [s, t, u])
                >>> contracted_expression = free.scalar_product(p1, p2)
                >>> assert kin.apply(contracted_expression) == s/2

                Parameters
                ----------
                expression : Expression
                    Expression with compact scalar products.
                """

            @typing.overload
            def apply(
                self,
                expression: pyo3_stub_gen.RustType["ConvertibleToExpression"],
            ) -> pyo3_stub_gen.RustType["PythonExpression"]:
                """
                Substitute scalar products without mutating global assumptions.

                Accepts a scalar Expression or a tensor expression. Tensor results retain
                the ordered open slots, including when the result is zero.

                Examples
                --------
                Using the setup in the ``Kinematics`` class example:

                >>> p1, p2, p3, p4, s, t, u = S("p1", "p2", "p3", "p4", "s", "t", "u")
                >>> free = hep.Kinematics(momenta=[p1, p2, p3, p4])
                >>> kin = hep.Kinematics.mandelstam([p1, p2, p3, p4], [E("0")]*4, [s, t, u])
                >>> contracted_expression = free.scalar_product(p1, p2)
                >>> assert kin.apply(contracted_expression) == s/2

                Parameters
                ----------
                expression : Expression
                    Expression with compact scalar products.
                """
        "#
    }
}

/// A spin projection along a particle's direction of motion.
///
/// HEP represents the physical minus, longitudinal, and plus helicity
/// states by the integers -1, 0, and 1.
///
/// Examples
/// --------
/// >>> from symbolica.community import hep
/// >>> incoming_fermion_helicity = hep.Helicity(-1)
/// >>> incoming_fermion_helicity == hep.Helicity.MINUS
/// True
///
/// Parameters
/// ----------
/// value : int
///     Helicity value; must be -1, 0, or 1.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Helicity",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Copy)]
pub struct PyHelicity {
    inner: Helicity,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyHelicity {
    /// Construct a physical helicity from -1, 0, or 1.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Helicity`` class example:
    ///
    /// >>> helicity = hep.Helicity(-1)
    ///
    /// Parameters
    /// ----------
    /// value : int
    ///     Integer helicity value.
    #[new]
    fn new(value: i8) -> PyResult<Self> {
        Helicity::try_from(value)
            .map(|inner| Self { inner })
            .map_err(error::kinematics)
    }

    /// Parse a named, signed, or integer helicity value.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Helicity`` class example:
    ///
    /// >>> hep.Helicity.parse("+") == hep.Helicity.PLUS
    /// True
    ///
    /// Parameters
    /// ----------
    /// value : str
    ///     Helicity spelling such as ``"-"``, ``"0"``, or ``"+"``.
    #[staticmethod]
    fn parse(value: &str) -> PyResult<Self> {
        value
            .parse()
            .map(|inner| Self { inner })
            .map_err(error::kinematics)
    }

    /// Return the negative-helicity singleton.
    ///
    /// Examples
    /// --------
    /// >>> int(Helicity.MINUS)
    /// -1
    ///
    #[classattr]
    #[pyo3(name = "MINUS")]
    fn minus() -> PyHelicity {
        Self {
            inner: Helicity::MINUS,
        }
    }

    /// Return the zero-helicity singleton.
    ///
    /// Examples
    /// --------
    /// >>> int(Helicity.ZERO)
    /// 0
    ///
    #[classattr]
    #[pyo3(name = "ZERO")]
    fn zero() -> PyHelicity {
        Self {
            inner: Helicity::ZERO,
        }
    }

    /// Return the positive-helicity singleton.
    ///
    /// Examples
    /// --------
    /// >>> int(Helicity.PLUS)
    /// 1
    ///
    #[classattr]
    #[pyo3(name = "PLUS")]
    fn plus() -> PyHelicity {
        Self {
            inner: Helicity::PLUS,
        }
    }

    /// Return the helicity as -1, 0, or 1.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Helicity`` class example:
    ///
    /// >>> hep.Helicity.PLUS.value
    /// 1
    #[getter]
    fn value(&self) -> i8 {
        self.inner.integer()
    }

    /// Convert this helicity to its integer value.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Helicity`` class example:
    ///
    /// >>> int(hep.Helicity.MINUS)
    /// -1
    fn __int__(&self) -> i8 {
        self.inner.integer()
    }

    /// Compare with another helicity value.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Helicity`` class example:
    ///
    /// >>> hep.Helicity(0) == hep.Helicity.ZERO
    /// True
    ///
    /// Parameters
    /// ----------
    /// other : object
    ///     Object to compare with this helicity value.
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        other
            .cast::<Self>()
            .is_ok_and(|other| self.inner == other.get().inner)
    }

    /// Return an evaluable-style representation of this helicity.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Helicity`` class example:
    ///
    /// >>> print(f"selected external helicity: {hep.Helicity.PLUS!r}")
    fn __repr__(&self) -> String {
        format!("Helicity({})", self.inner.integer())
    }
}

/// A Cartesian axis used to specify spatial rotations.
///
/// Examples
/// --------
/// >>> from symbolica.community import hep
/// >>> beam_axis = hep.Axis.Z
/// >>> rotation = hep.Rotation.quarter_turn(beam_axis)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    name = "Axis",
    module = "symbolica.community.feynkit",
    frozen,
    eq,
    eq_int,
    from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PyAxis {
    X,
    Y,
    Z,
}

impl From<PyAxis> for Axis {
    fn from(value: PyAxis) -> Self {
        match value {
            PyAxis::X => Self::X,
            PyAxis::Y => Self::Y,
            PyAxis::Z => Self::Z,
        }
    }
}

/// A sequential-recombination algorithm for collider jet clustering.
///
/// The available choices are kT, Cambridge--Aachen, and anti-kT.
///
/// Examples
/// --------
/// >>> from symbolica.community import hep
/// >>> algorithm = hep.JetAlgorithm.AntiKt
/// >>> definition = hep.JetDefinition(algorithm, radius=0.4)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    name = "JetAlgorithm",
    module = "symbolica.community.feynkit",
    frozen,
    eq,
    eq_int,
    from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PyJetAlgorithm {
    Kt,
    CambridgeAachen,
    AntiKt,
}

impl From<PyJetAlgorithm> for JetAlgorithm {
    fn from(value: PyJetAlgorithm) -> Self {
        match value {
            PyJetAlgorithm::Kt => Self::Kt,
            PyJetAlgorithm::CambridgeAachen => Self::CambridgeAachen,
            PyJetAlgorithm::AntiKt => Self::AntiKt,
        }
    }
}

impl From<JetAlgorithm> for PyJetAlgorithm {
    fn from(value: JetAlgorithm) -> Self {
        match value {
            JetAlgorithm::Kt => Self::Kt,
            JetAlgorithm::CambridgeAachen => Self::CambridgeAachen,
            JetAlgorithm::AntiKt => Self::AntiKt,
        }
    }
}

/// A Cartesian spatial momentum ``(px, py, pz)``.
///
/// Three-momenta are used for boost velocities, spatial rotations, and the
/// three-vector part of a relativistic four-momentum.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> momentum = hep.ThreeMomentum(3.0, 4.0, 0.0)
/// >>> p = momentum
/// >>> first = hep.ThreeMomentum(0.0, 1.0, 0.0)
/// >>> second = hep.ThreeMomentum(1.0, 0.0, 0.0)
/// >>> assert momentum.pt == 5.0
///
/// Parameters
/// ----------
/// px : float
///     Momentum component along the x axis.
/// py : float
///     Momentum component along the y axis.
/// pz : float
///     Momentum component along the z axis.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ThreeMomentum",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyThreeMomentum {
    pub(crate) inner: ThreeMomentum<f64>,
}

impl From<ThreeMomentum<f64>> for PyThreeMomentum {
    fn from(inner: ThreeMomentum<f64>) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyThreeMomentum {
    /// Construct a Cartesian three-momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> momentum = hep.ThreeMomentum(3.0, 4.0, 0.0)
    ///
    /// Parameters
    /// ----------
    /// px : float
    ///     Momentum along the x axis.
    /// py : float
    ///     Momentum along the y axis.
    /// pz : float
    ///     Momentum along the z axis.
    #[new]
    fn new(px: f64, py: f64, pz: f64) -> Self {
        ThreeMomentum::new(px, py, pz).into()
    }

    /// Return the x component.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> assert momentum.px == 3.0
    #[getter]
    fn px(&self) -> f64 {
        self.inner.px
    }
    /// Return the y component.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> assert momentum.py == 4.0
    #[getter]
    fn py(&self) -> f64 {
        self.inner.py
    }
    /// Return the z component.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> assert momentum.pz == 0.0
    #[getter]
    fn pz(&self) -> f64 {
        self.inner.pz
    }
    /// Return the squared Euclidean norm.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> hep.ThreeMomentum(3.0, 4.0, 0.0).norm_squared
    /// 25.0
    #[getter]
    fn norm_squared(&self) -> f64 {
        self.inner.norm_squared()
    }
    /// Return the Euclidean norm.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> hep.ThreeMomentum(3.0, 4.0, 0.0).norm
    /// 5.0
    #[getter]
    fn norm(&self) -> f64 {
        self.inner.norm()
    }
    /// Return the transverse-momentum magnitude.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> hep.ThreeMomentum(3.0, 4.0, 12.0).pt
    /// 5.0
    #[getter]
    fn pt(&self) -> f64 {
        self.inner.pt()
    }
    /// Return the azimuthal angle in radians.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> hep.ThreeMomentum(1.0, 0.0, 0.0).phi
    /// 0.0
    #[getter]
    fn phi(&self) -> f64 {
        self.inner.phi()
    }
    /// Return the pseudorapidity derived from the momentum direction.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> abs(hep.ThreeMomentum(1.0, 0.0, 0.0).pseudorapidity) < 1e-12
    /// True
    #[getter]
    fn pseudorapidity(&self) -> f64 {
        self.inner.pseudorapidity()
    }
    /// Return the Euclidean dot product with another three-momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> hep.ThreeMomentum(1.0, 0.0, 0.0).dot(hep.ThreeMomentum(2.0, 0.0, 0.0))
    /// 2.0
    ///
    /// Parameters
    /// ----------
    /// other : ThreeMomentum
    ///     Momentum to contract with this one.
    fn dot(&self, other: &Self) -> f64 {
        self.inner.dot(&other.inner)
    }
    /// Return the vector cross product with another three-momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> hep.ThreeMomentum(1.0, 0.0, 0.0).cross(hep.ThreeMomentum(0.0, 1.0, 0.0))
    /// ThreeMomentum(0, 0, 1)
    ///
    /// Parameters
    /// ----------
    /// other : ThreeMomentum
    ///     Right-hand operand of the cross product.
    fn cross(&self, other: &Self) -> Self {
        self.inner.cross(&other.inner).into()
    }
    /// Return the wrapped azimuthal separation from another momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> first.delta_phi(second)
    /// 1.5707963267948966
    ///
    /// Parameters
    /// ----------
    /// other : ThreeMomentum
    ///     Momentum whose azimuth is compared.
    fn delta_phi(&self, other: &Self) -> f64 {
        self.inner.delta_phi(&other.inner)
    }
    /// Return the angular distance in pseudorapidity-azimuth space.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> separation = first.delta_r(second)
    /// >>> passes_isolation = separation > 0.4
    ///
    /// Parameters
    /// ----------
    /// other : ThreeMomentum
    ///     Momentum to compare with this one.
    fn delta_r(&self, other: &Self) -> f64 {
        self.inner.delta_r(&other.inner)
    }

    /// Add two spatial momenta component by component.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> total = first + second
    ///
    /// Parameters
    /// ----------
    /// other : ThreeMomentum
    ///     Momentum to add.
    fn __add__(&self, other: &Self) -> Self {
        (self.inner + other.inner).into()
    }

    /// Subtract another spatial momentum component by component.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> incoming = hep.ThreeMomentum(3.0, 4.0, 0.0)
    /// >>> outgoing = hep.ThreeMomentum(1.0, 0.0, 0.0)
    /// >>> transfer = incoming - outgoing
    /// >>> assert transfer.px == 2.0
    ///
    /// Parameters
    /// ----------
    /// other : ThreeMomentum
    ///     Momentum to subtract.
    fn __sub__(&self, other: &Self) -> Self {
        (self.inner - other.inner).into()
    }

    /// Reverse every spatial momentum component.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> incoming_convention = hep.ThreeMomentum(3.0, 4.0, 0.0)
    /// >>> outgoing_convention = -incoming_convention
    /// >>> assert outgoing_convention.px == -3.0
    fn __neg__(&self) -> Self {
        (-self.inner).into()
    }

    /// Scale every spatial momentum component.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> half_momentum = momentum * 0.5
    ///
    /// Parameters
    /// ----------
    /// scalar : float
    ///     Multiplicative scale factor.
    fn __mul__(&self, scalar: f64) -> Self {
        (self.inner * scalar).into()
    }

    /// Scale every spatial momentum component from the left.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> half_momentum = 0.5 * momentum
    ///
    /// Parameters
    /// ----------
    /// scalar : float
    ///     Multiplicative scale factor.
    fn __rmul__(&self, scalar: f64) -> Self {
        (self.inner * scalar).into()
    }

    /// Lift this spatial momentum to an on-shell four-momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> hep.ThreeMomentum(3.0, 4.0, 0.0).on_shell().energy
    /// 5.0
    ///
    /// Parameters
    /// ----------
    /// mass : float, optional
    ///     On-shell mass; omitted values describe a massless momentum.
    #[pyo3(signature = (mass=None))]
    fn on_shell(&self, mass: Option<f64>) -> PyFourMomentum {
        self.inner.on_shell(mass.as_ref()).into()
    }

    /// Return a constructor-style representation of the components.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> track_momentum = hep.ThreeMomentum(1.0, 2.0, 3.0)
    /// >>> print(f"track momentum: {track_momentum!r}")
    fn __repr__(&self) -> String {
        format!(
            "ThreeMomentum({}, {}, {})",
            self.inner.px, self.inner.py, self.inner.pz
        )
    }

    /// Render the momentum as a mathematical three-vector.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> latex = momentum._repr_latex_()
    fn _repr_latex_(&self) -> String {
        format!(
            r"$\vec{{p}}=\left({},{},{}\right)$",
            self.inner.px, self.inner.py, self.inner.pz
        )
    }

    /// Write the constructor-style form to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ThreeMomentum`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(momentum)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
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

/// A relativistic four-momentum in ``(energy, px, py, pz)`` order.
///
/// The class provides collider observables, invariant products, rotations, and
/// boosts while keeping the component convention explicit.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> momentum = hep.FourMomentum(5.0, 3.0, 4.0, 0.0)
/// >>> p = momentum
/// >>> first = hep.FourMomentum(5.0, 0.0, 5.0, 0.0)
/// >>> second = hep.FourMomentum(5.0, 5.0, 0.0, 0.0)
/// >>> assert momentum.mass_squared == 0.0
///
/// Parameters
/// ----------
/// energy : float
///     Energy component.
/// px : float
///     Momentum component along the x axis.
/// py : float
///     Momentum component along the y axis.
/// pz : float
///     Momentum component along the z axis.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "FourMomentum",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyFourMomentum {
    pub(crate) inner: FourMomentum<f64>,
}

impl From<FourMomentum<f64>> for PyFourMomentum {
    fn from(inner: FourMomentum<f64>) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyFourMomentum {
    /// Construct a four-momentum from energy and Cartesian spatial components.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> momentum = hep.FourMomentum(5.0, 3.0, 4.0, 0.0)
    ///
    /// Parameters
    /// ----------
    /// energy : float
    ///     Temporal component.
    /// px : float
    ///     Momentum along the x axis.
    /// py : float
    ///     Momentum along the y axis.
    /// pz : float
    ///     Momentum along the z axis.
    #[new]
    fn new(energy: f64, px: f64, py: f64, pz: f64) -> Self {
        FourMomentum::from_args(energy, px, py, pz).into()
    }

    /// Return the energy component.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> assert momentum.energy == 5.0
    #[getter]
    fn energy(&self) -> f64 {
        self.inner.temporal.value
    }
    /// Return the x component of spatial momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> assert momentum.px == 3.0
    #[getter]
    fn px(&self) -> f64 {
        self.inner.spatial.px
    }
    /// Return the y component of spatial momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> assert momentum.py == 4.0
    #[getter]
    fn py(&self) -> f64 {
        self.inner.spatial.py
    }
    /// Return the z component of spatial momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> assert momentum.pz == 0.0
    #[getter]
    fn pz(&self) -> f64 {
        self.inner.spatial.pz
    }
    /// Return the spatial three-momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> hep.FourMomentum(5.0, 3.0, 4.0, 0.0).spatial.pt
    /// 5.0
    #[getter]
    fn spatial(&self) -> PyThreeMomentum {
        self.inner.spatial.into()
    }
    /// Return ``(energy, px, py, pz)``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> hep.FourMomentum(5.0, 3.0, 4.0, 0.0).components()
    /// (5.0, 3.0, 4.0, 0.0)
    fn components(&self) -> (f64, f64, f64, f64) {
        self.inner.into()
    }
    /// Return the Minkowski dot product with another four-momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> momentum.dot(momentum)
    /// 0.0
    ///
    /// Parameters
    /// ----------
    /// other : FourMomentum
    ///     Momentum to contract with this one.
    fn dot(&self, other: &Self) -> f64 {
        self.inner.dot(&other.inner)
    }
    /// Return the invariant mass squared.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> hep.FourMomentum(5.0, 3.0, 4.0, 0.0).mass_squared
    /// 0.0
    #[getter]
    fn mass_squared(&self) -> f64 {
        self.inner.mass_squared()
    }

    /// Return the decay denominator ``2E`` or the invariant two-particle flux.
    ///
    /// Inputs must be physical, future-directed on-shell momenta. Decay rates
    /// refer to this momentum's frame; the rest-frame result is ``2M``.
    /// No unit conversion, spin/color average or symmetry factor is included.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> p = hep.FourMomentum(5.0, 0.0, 0.0, 5.0)
    /// >>> q = hep.FourMomentum(5.0, 0.0, 0.0, -5.0)
    /// >>> p.flux(q)
    /// 200.0
    ///
    /// Parameters
    /// ----------
    /// other : FourMomentum | None
    ///     Other incoming momentum; None selects a decay in this frame.
    #[pyo3(signature = (other=None))]
    fn flux(&self, other: Option<&PyFourMomentum>) -> f64 {
        self.inner.flux(other.map(|p| &p.inner))
    }

    /// Return the invariant mass.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> hep.FourMomentum(5.0, 0.0, 0.0, 0.0).mass
    /// 5.0
    #[getter]
    fn mass(&self) -> f64 {
        self.inner.mass()
    }
    /// Return the transverse-momentum magnitude.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> hep.FourMomentum(13.0, 3.0, 4.0, 12.0).pt
    /// 5.0
    #[getter]
    fn pt(&self) -> f64 {
        self.inner.pt()
    }
    /// Return the azimuthal angle in radians.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> hep.FourMomentum(1.0, 1.0, 0.0, 0.0).phi
    /// 0.0
    #[getter]
    fn phi(&self) -> f64 {
        self.inner.phi()
    }
    /// Return the pseudorapidity of the spatial momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> abs(hep.FourMomentum(1.0, 1.0, 0.0, 0.0).pseudorapidity) < 1e-12
    /// True
    #[getter]
    fn pseudorapidity(&self) -> f64 {
        self.inner.pseudorapidity()
    }
    /// Return the longitudinal rapidity.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> hep.FourMomentum(1.0, 1.0, 0.0, 0.0).rapidity
    /// 0.0
    #[getter]
    fn rapidity(&self) -> f64 {
        self.inner.rapidity()
    }
    /// Return the wrapped azimuthal separation from another momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> first.delta_phi(second)
    /// 1.5707963267948966
    ///
    /// Parameters
    /// ----------
    /// other : FourMomentum
    ///     Momentum whose azimuth is compared.
    fn delta_phi(&self, other: &Self) -> f64 {
        self.inner.delta_phi(&other.inner)
    }
    /// Return the distance in rapidity-azimuth space.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> separation = first.delta_r(second)
    /// >>> same_jet = separation < 0.4
    ///
    /// Parameters
    /// ----------
    /// other : FourMomentum
    ///     Momentum to compare with this one.
    fn delta_r(&self, other: &Self) -> f64 {
        self.inner.delta_r(&other.inner)
    }
    /// Add two four-momenta component by component.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> (first + second).energy == first.energy + second.energy
    /// True
    ///
    /// Parameters
    /// ----------
    /// other : FourMomentum
    ///     Momentum to add.
    fn __add__(&self, other: &Self) -> Self {
        (self.inner + other.inner).into()
    }
    /// Subtract another four-momentum component by component.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> (first - second).energy == first.energy - second.energy
    /// True
    ///
    /// Parameters
    /// ----------
    /// other : FourMomentum
    ///     Momentum to subtract.
    fn __sub__(&self, other: &Self) -> Self {
        (self.inner - other.inner).into()
    }

    /// Reverse the four-momentum flow convention.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> incoming_convention = hep.FourMomentum(5.0, 3.0, 4.0, 0.0)
    /// >>> outgoing_convention = -incoming_convention
    /// >>> assert outgoing_convention.energy == -5.0
    fn __neg__(&self) -> Self {
        (-self.inner).into()
    }

    /// Scale every four-momentum component.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> half_momentum = momentum * 0.5
    ///
    /// Parameters
    /// ----------
    /// scalar : float
    ///     Multiplicative scale factor.
    fn __mul__(&self, scalar: f64) -> Self {
        (self.inner * scalar).into()
    }

    /// Scale every four-momentum component from the left.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> half_momentum = 0.5 * momentum
    ///
    /// Parameters
    /// ----------
    /// scalar : float
    ///     Multiplicative scale factor.
    fn __rmul__(&self, scalar: f64) -> Self {
        (self.inner * scalar).into()
    }

    /// Return a constructor-style representation of the components.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> muon_momentum = hep.FourMomentum(5.0, 3.0, 4.0, 0.0)
    /// >>> print(f"muon four-momentum: {muon_momentum!r}")
    fn __repr__(&self) -> String {
        let (energy, px, py, pz): (f64, f64, f64, f64) = self.inner.into();
        format!("FourMomentum({energy}, {px}, {py}, {pz})")
    }

    /// Render the momentum as a contravariant four-vector.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> latex = momentum._repr_latex_()
    fn _repr_latex_(&self) -> String {
        let (energy, px, py, pz): (f64, f64, f64, f64) = self.inner.into();
        format!(r"$p^\mu=\left({energy},{px},{py},{pz}\right)$")
    }

    /// Write the constructor-style form to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FourMomentum`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(momentum)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
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

/// A spatial rotation acting on three- and four-momenta.
///
/// Rotations may be built from Euler angles, an identity transformation, or a
/// quarter turn around a Cartesian axis.
///
/// Examples
/// --------
/// >>> from symbolica.community import hep
/// >>> rotation = hep.Rotation.quarter_turn(hep.Axis.Z)
/// >>> rotated = rotation.apply_three(hep.ThreeMomentum(1.0, 0.0, 0.0))
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Rotation",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyRotation {
    inner: Rotation<f64>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyRotation {
    /// Construct the identity rotation.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Rotation`` class example:
    ///
    /// >>> rotation = hep.Rotation.identity()
    #[staticmethod]
    fn identity() -> Self {
        Self {
            inner: Rotation::Identity,
        }
    }

    /// Construct a rotation from three Euler angles in radians.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Rotation`` class example:
    ///
    /// >>> rotation = hep.Rotation.euler(0.1, 0.2, 0.3)
    ///
    /// Parameters
    /// ----------
    /// alpha : float
    ///     First Euler angle.
    /// beta : float
    ///     Second Euler angle.
    /// gamma : float
    ///     Third Euler angle.
    #[staticmethod]
    fn euler(alpha: f64, beta: f64, gamma: f64) -> Self {
        Self {
            inner: Rotation::euler(alpha, beta, gamma),
        }
    }

    /// Construct a positive quarter turn around a Cartesian axis.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Rotation`` class example:
    ///
    /// >>> rotation = hep.Rotation.quarter_turn(hep.Axis.Z)
    ///
    /// Parameters
    /// ----------
    /// axis : Axis
    ///     Axis about which to rotate by pi/2.
    #[staticmethod]
    fn quarter_turn(axis: PyAxis) -> Self {
        Self {
            inner: Rotation::quarter_turn(axis.into()),
        }
    }

    /// Apply this rotation to a three-momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Rotation`` class example:
    ///
    /// >>> momentum = hep.ThreeMomentum(1.0, 0.0, 0.0)
    /// >>> rotated = rotation.apply_three(momentum)
    ///
    /// Parameters
    /// ----------
    /// momentum : ThreeMomentum
    ///     Spatial momentum to rotate.
    fn apply_three(&self, momentum: &PyThreeMomentum) -> PyThreeMomentum {
        self.inner.rotate_three(&momentum.inner).into()
    }

    /// Apply this spatial rotation to a four-momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Rotation`` class example:
    ///
    /// >>> momentum = hep.FourMomentum(2.0, 1.0, 0.0, 0.0)
    /// >>> rotated = rotation.apply_four(momentum)
    ///
    /// Parameters
    /// ----------
    /// momentum : FourMomentum
    ///     Four-momentum whose spatial components are rotated.
    fn apply_four(&self, momentum: &PyFourMomentum) -> PyFourMomentum {
        self.inner.rotate_four(&momentum.inner).into()
    }

    /// Apply the inverse rotation to a three-momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Rotation`` class example:
    ///
    /// >>> momentum = hep.ThreeMomentum(1.0, 0.0, 0.0)
    /// >>> original = rotation.apply_inverse_three(rotation.apply_three(momentum))
    ///
    /// Parameters
    /// ----------
    /// momentum : ThreeMomentum
    ///     Spatial momentum to inverse-rotate.
    fn apply_inverse_three(&self, momentum: &PyThreeMomentum) -> PyThreeMomentum {
        self.inner.inverse_rotate_three(&momentum.inner).into()
    }

    /// Apply the inverse spatial rotation to a four-momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Rotation`` class example:
    ///
    /// >>> momentum = hep.FourMomentum(2.0, 1.0, 0.0, 0.0)
    /// >>> original = rotation.apply_inverse_four(rotation.apply_four(momentum))
    ///
    /// Parameters
    /// ----------
    /// momentum : FourMomentum
    ///     Four-momentum whose spatial components are inverse-rotated.
    fn apply_inverse_four(&self, momentum: &PyFourMomentum) -> PyFourMomentum {
        self.inner.inverse_rotate_four(&momentum.inner).into()
    }
}

/// A proper Lorentz boost specified by a three-velocity ``beta``.
///
/// Use boosts to move four-momenta between the laboratory frame and a useful
/// rest frame, with units chosen so that the speed of light is one.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> momentum = hep.FourMomentum(10.0, 0.0, 0.0, 0.0)
/// >>> boost = hep.Boost(hep.ThreeMomentum(0.0, 0.0, 0.5))
/// >>> boosted = boost.apply(momentum)
/// >>> assert abs(boosted.mass_squared - momentum.mass_squared) < 1e-10
///
/// Parameters
/// ----------
/// beta : ThreeMomentum
///     Dimensionless boost velocity, whose magnitude must be smaller than one.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Boost",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyBoost {
    inner: Boost<f64>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyBoost {
    /// Construct a Lorentz boost from a dimensionless three-velocity.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Boost`` class example:
    ///
    /// >>> boost = hep.Boost(hep.ThreeMomentum(0.0, 0.0, 0.5))
    ///
    /// Parameters
    /// ----------
    /// beta : ThreeMomentum
    ///     Boost velocity in units where the speed of light is one.
    #[new]
    fn new(beta: &PyThreeMomentum) -> PyResult<Self> {
        Boost::new(beta.inner)
            .map(|inner| Self { inner })
            .map_err(error::kinematics)
    }

    /// Return the dimensionless boost velocity.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Boost`` class example:
    ///
    /// >>> assert boost.beta.pz == 0.5
    #[getter]
    fn beta(&self) -> PyThreeMomentum {
        self.inner.beta().to_owned().into()
    }
    /// Apply this Lorentz boost to a four-momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Boost`` class example:
    ///
    /// >>> boosted = boost.apply(momentum)
    ///
    /// Parameters
    /// ----------
    /// momentum : FourMomentum
    ///     Four-momentum to boost.
    fn apply(&self, momentum: &PyFourMomentum) -> PyFourMomentum {
        self.inner.apply(&momentum.inner).into()
    }
    /// Apply the inverse Lorentz boost to a four-momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Boost`` class example:
    ///
    /// >>> original = boost.apply_inverse(boost.apply(momentum))
    ///
    /// Parameters
    /// ----------
    /// momentum : FourMomentum
    ///     Four-momentum to inverse-boost.
    fn apply_inverse(&self, momentum: &PyFourMomentum) -> PyFourMomentum {
        self.inner.apply_inverse(&momentum.inner).into()
    }
}

/// A reconstructed collider jet and its input-particle constituents.
///
/// Jets are returned by ``JetDefinition.cluster`` and expose the recombined
/// four-momentum together with the indices of their original inputs.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> particles = [hep.FourMomentum(50.0, 30.0, 40.0, 0.0),
/// ...              hep.FourMomentum(25.0, -15.0, -20.0, 0.0)]
/// >>> definition = hep.JetDefinition.anti_kt(radius=0.4, minimum_pt=20.0)
/// >>> clustering = definition.cluster(particles)
/// >>> jets = clustering.jets
/// >>> leading_jet = jets[0]
/// >>> assert leading_jet.pt == 50.0
/// >>> inputs = [particles[i] for i in leading_jet.constituent_indices]
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Jet",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyJet {
    inner: Jet<f64>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyJet {
    /// Return the recombined four-momentum of this jet.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Jet`` class example:
    ///
    /// >>> assert leading_jet.momentum.pt == leading_jet.pt
    #[getter]
    fn momentum(&self) -> PyFourMomentum {
        self.inner.momentum.into()
    }
    /// Return sorted positions of the input momenta assigned to this jet.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Jet`` class example:
    ///
    /// >>> inputs = [particles[i] for i in leading_jet.constituent_indices]
    #[getter]
    fn constituent_indices(&self) -> Vec<usize> {
        self.inner.constituent_indices().to_vec()
    }
    /// Return the jet transverse momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Jet`` class example:
    ///
    /// >>> assert leading_jet.pt == 50.0
    #[getter]
    fn pt(&self) -> f64 {
        self.inner.pt()
    }
    /// Return the jet rapidity.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Jet`` class example:
    ///
    /// >>> assert leading_jet.rapidity == 0.0
    #[getter]
    fn rapidity(&self) -> f64 {
        self.inner.rapidity()
    }
    /// Return the jet azimuthal angle in radians.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Jet`` class example:
    ///
    /// >>> azimuth = leading_jet.phi
    #[getter]
    fn phi(&self) -> f64 {
        self.inner.phi()
    }

    /// Return a concise jet summary with its kinematics and constituent count.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Jet`` class example:
    ///
    /// >>> summary = repr(leading_jet)
    fn __repr__(&self) -> String {
        format!(
            "Jet(pt={}, rapidity={}, phi={}, constituents={})",
            self.inner.pt(),
            self.inner.rapidity(),
            self.inner.phi(),
            self.inner.constituent_indices().len(),
        )
    }

    /// Render the jet kinematics as a compact notebook table.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Jet`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(leading_jet)
    fn _repr_html_(&self) -> String {
        format!(
            "<table class=\"feynkit-jet\" style=\"border-collapse:collapse\"><thead><tr>\
             <th style=\"padding:.2rem .55rem;text-align:right\">p<sub>T</sub></th>\
             <th style=\"padding:.2rem .55rem;text-align:right\">y</th><th style=\"\
             padding:.2rem .55rem;text-align:right\">&phi;</th><th style=\"padding:.2rem \
             .55rem;text-align:left\">constituents</th></tr></thead><tbody><tr><td style=\"\
             padding:.2rem .55rem;text-align:right\">{:.6}</td><td style=\"padding:.2rem \
             .55rem;text-align:right\">{:.6}</td><td style=\"padding:.2rem .55rem;\
             text-align:right\">{:.6}</td><td style=\"padding:.2rem .55rem\"><code>{:?}\
             </code></td></tr></tbody></table>",
            self.inner.pt(),
            self.inner.rapidity(),
            self.inner.phi(),
            self.inner.constituent_indices(),
        )
    }

    /// Write the concise jet summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Jet`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(leading_jet)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
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

/// The selected jets from one clustering operation.
///
/// Jets are ordered by decreasing transverse momentum, and each jet's
/// ``constituent_indices`` map back to positions in the supplied momentum list.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> particles = [hep.FourMomentum(50.0, 30.0, 40.0, 0.0),
/// ...              hep.FourMomentum(25.0, -15.0, -20.0, 0.0)]
/// >>> definition = hep.JetDefinition.anti_kt(radius=0.4, minimum_pt=20.0)
/// >>> clustering = definition.cluster(particles)
/// >>> jets = clustering.jets
/// >>> assert len(jets) == 2
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ClusteringResult",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyClusteringResult {
    inner: ClusteringResult<f64>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyClusteringResult {
    /// Return inclusive jets ordered by decreasing transverse momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ClusteringResult`` class example:
    ///
    /// >>> jets = clustering.jets
    /// >>> all(left.pt >= right.pt for left, right in zip(jets, jets[1:]))
    /// True
    #[getter]
    fn jets(&self) -> Vec<PyJet> {
        self.inner
            .jets
            .iter()
            .cloned()
            .map(|inner| PyJet { inner })
            .collect()
    }
    /// Return the number of clustered jets.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ClusteringResult`` class example:
    ///
    /// >>> jet_multiplicity = len(clustering)
    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Return one jet by decreasing-transverse-momentum index.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ClusteringResult`` class example:
    ///
    /// >>> leading_jet = clustering[0]
    ///
    /// Parameters
    /// ----------
    /// index : int
    ///     Zero-based index; negative indices count from the end.
    fn __getitem__(&self, index: isize) -> PyResult<PyJet> {
        let length = self.inner.jets.len() as isize;
        let index = if index < 0 { length + index } else { index };
        if !(0..length).contains(&index) {
            return Err(PyIndexError::new_err("jet index out of range"));
        }
        Ok(PyJet {
            inner: self.inner.jets[index as usize].clone(),
        })
    }

    /// Iterate over jets from highest to lowest transverse momentum.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ClusteringResult`` class example:
    ///
    /// >>> transverse_momenta = [jet.pt for jet in clustering]
    #[gen_stub(override_return_type(
        type_repr = "collections.abc.Iterator[Jet]",
        imports = ("collections.abc")
    ))]
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyList::new(py, self.jets())?.call_method0("__iter__")
    }

    /// Return a concise summary of the clustered jet collection.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ClusteringResult`` class example:
    ///
    /// >>> print(clustering)
    fn __repr__(&self) -> String {
        format!("ClusteringResult(jets={})", self.inner.len())
    }

    /// Render a bounded table of clustered jet kinematics.
    ///
    /// At most 20 jets are included so notebook display remains responsive for
    /// unusually large clustering results.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ClusteringResult`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(clustering)
    fn _repr_html_(&self) -> String {
        const DISPLAY_LIMIT: usize = 20;
        let rows = self
            .inner
            .jets
            .iter()
            .take(DISPLAY_LIMIT)
            .enumerate()
            .map(|(index, jet)| {
                format!(
                    "<tr><td style=\"padding:.2rem .55rem;text-align:right\">{index}</td>\
                     <td style=\"padding:.2rem .55rem;text-align:right\">{:.6}</td><td \
                     style=\"padding:.2rem .55rem;text-align:right\">{:.6}</td><td style=\"\
                     padding:.2rem .55rem;text-align:right\">{:.6}</td><td style=\"padding:\
                     .2rem .55rem\"><code>{:?}</code></td></tr>",
                    jet.pt(),
                    jet.rapidity(),
                    jet.phi(),
                    jet.constituent_indices(),
                )
            })
            .collect::<String>();
        let omitted = self.inner.len().saturating_sub(DISPLAY_LIMIT);
        let note = if omitted > 0 {
            format!("<div style=\"opacity:.7\">{omitted} additional jets omitted</div>")
        } else {
            String::new()
        };
        format!(
            "<div class=\"feynkit-clustering-result\" style=\"display:inline-block;max-width:\
             100%;overflow-x:auto\"><strong>Clustered jets ({})</strong><table style=\"\
             border-collapse:collapse;margin-top:.3rem\"><thead><tr><th style=\"padding:.2rem \
             .55rem;text-align:right\">#</th><th style=\"padding:.2rem .55rem;text-align:\
             right\">p<sub>T</sub></th><th style=\"padding:.2rem .55rem;text-align:right\">\
             y</th><th style=\"padding:.2rem .55rem;text-align:right\">&phi;</th><th style=\"\
             padding:.2rem .55rem;text-align:left\">constituents</th></tr></thead><tbody>{rows}\
             </tbody></table>{note}</div>",
            self.inner.len(),
        )
    }

    /// Write the concise collection summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ClusteringResult`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(clustering)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
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

/// A generalized-kT jet definition for sequential recombination.
///
/// The definition selects the distance measure, jet radius, and transverse-
/// momentum threshold used to reconstruct jets from final-state momenta.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> particles = [hep.FourMomentum(50.0, 30.0, 40.0, 0.0),
/// ...              hep.FourMomentum(25.0, -15.0, -20.0, 0.0)]
/// >>> definition = hep.JetDefinition.anti_kt(radius=0.4, minimum_pt=20.0)
/// >>> clustering = definition.cluster(particles)
/// >>> assert definition.radius == 0.4
/// >>> assert len(clustering.jets) == 2
///
/// Parameters
/// ----------
/// algorithm : JetAlgorithm
///     Generalized-kT algorithm to use.
/// radius : float
///     Jet-radius parameter.
/// minimum_pt : float, optional
///     Minimum transverse momentum of returned jets.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "JetDefinition",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyJetDefinition {
    inner: JetDefinition<f64>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyJetDefinition {
    /// Construct a generalized-kT jet definition.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``JetDefinition`` class example:
    ///
    /// >>> definition = hep.JetDefinition(hep.JetAlgorithm.AntiKt, 0.4, 20.0)
    ///
    /// Parameters
    /// ----------
    /// algorithm : JetAlgorithm
    ///     Generalized-kT algorithm.
    /// radius : float
    ///     Jet-radius parameter.
    /// minimum_pt : float, optional
    ///     Minimum transverse momentum for returned jets.
    #[new]
    #[pyo3(signature = (algorithm, radius, minimum_pt=0.0))]
    fn new(algorithm: PyJetAlgorithm, radius: f64, minimum_pt: f64) -> Self {
        Self {
            inner: JetDefinition::new(algorithm.into(), radius).with_minimum_pt(minimum_pt),
        }
    }

    /// Construct a kT jet definition.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``JetDefinition`` class example:
    ///
    /// >>> definition = hep.JetDefinition.kt(0.4)
    ///
    /// Parameters
    /// ----------
    /// radius : float
    ///     Jet-radius parameter.
    /// minimum_pt : float, optional
    ///     Minimum transverse momentum for returned jets.
    #[staticmethod]
    #[pyo3(signature = (radius, minimum_pt=0.0))]
    fn kt(radius: f64, minimum_pt: f64) -> Self {
        Self {
            inner: JetDefinition::kt(radius).with_minimum_pt(minimum_pt),
        }
    }

    /// Construct a Cambridge-Aachen jet definition.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``JetDefinition`` class example:
    ///
    /// >>> definition = hep.JetDefinition.cambridge_aachen(0.4)
    ///
    /// Parameters
    /// ----------
    /// radius : float
    ///     Jet-radius parameter.
    /// minimum_pt : float, optional
    ///     Minimum transverse momentum for returned jets.
    #[staticmethod]
    #[pyo3(signature = (radius, minimum_pt=0.0))]
    fn cambridge_aachen(radius: f64, minimum_pt: f64) -> Self {
        Self {
            inner: JetDefinition::cambridge_aachen(radius).with_minimum_pt(minimum_pt),
        }
    }

    /// Construct an anti-kT jet definition.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``JetDefinition`` class example:
    ///
    /// >>> definition = hep.JetDefinition.anti_kt(0.4, minimum_pt=20.0)
    ///
    /// Parameters
    /// ----------
    /// radius : float
    ///     Jet-radius parameter.
    /// minimum_pt : float, optional
    ///     Minimum transverse momentum for returned jets.
    #[staticmethod]
    #[pyo3(signature = (radius, minimum_pt=0.0))]
    fn anti_kt(radius: f64, minimum_pt: f64) -> Self {
        Self {
            inner: JetDefinition::anti_kt(radius).with_minimum_pt(minimum_pt),
        }
    }

    /// Return the selected clustering algorithm.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``JetDefinition`` class example:
    ///
    /// >>> assert definition.algorithm == hep.JetAlgorithm.AntiKt
    #[getter]
    fn algorithm(&self) -> PyJetAlgorithm {
        self.inner.algorithm().into()
    }
    /// Return the jet-radius parameter.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``JetDefinition`` class example:
    ///
    /// >>> assert definition.radius == 0.4
    #[getter]
    fn radius(&self) -> f64 {
        *self.inner.radius()
    }
    /// Return the minimum transverse momentum for retained jets.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``JetDefinition`` class example:
    ///
    /// >>> assert definition.minimum_pt == 20.0
    #[getter]
    fn minimum_pt(&self) -> f64 {
        *self.inner.minimum_pt()
    }

    /// Cluster four-momenta with this sequential-recombination definition.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``JetDefinition`` class example:
    ///
    /// >>> result = definition.cluster(particles)
    /// >>> assert len(result.jets) == 2
    ///
    /// Parameters
    /// ----------
    /// momenta : sequence of FourMomentum
    ///     Input four-momenta to cluster.
    fn cluster(
        &self,
        py: Python<'_>,
        momenta: Vec<PyFourMomentum>,
    ) -> PyResult<PyClusteringResult> {
        let definition = self.inner.clone();
        let momenta = momenta
            .into_iter()
            .map(|momentum| momentum.inner)
            .collect::<Vec<_>>();
        py.detach(move || definition.cluster(&momenta))
            .map(|inner| PyClusteringResult { inner })
            .map_err(error::kinematics)
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyKinematics>()?;
    module.add_class::<PyHelicity>()?;
    module.add_class::<PyAxis>()?;
    module.add_class::<PyJetAlgorithm>()?;
    module.add_class::<PyThreeMomentum>()?;
    module.add_class::<PyFourMomentum>()?;
    module.add_class::<PyRotation>()?;
    module.add_class::<PyBoost>()?;
    module.add_class::<PyJet>()?;
    module.add_class::<PyClusteringResult>()?;
    module.add_class::<PyJetDefinition>()?;
    Ok(())
}
