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
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

use crate::error;

/// A spin projection along a particle's direction of motion.
///
/// FeynKit represents the physical minus, longitudinal, and plus helicity
/// states by the integers -1, 0, and 1.
///
/// Examples
/// --------
/// >>> import symbolica.community.feynkit as fk
/// >>> incoming_fermion_helicity = fk.Helicity(-1)
/// >>> incoming_fermion_helicity == fk.Helicity.MINUS
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
    /// >>> helicity = Helicity(-1)
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
    /// >>> Helicity.parse("+") == Helicity.PLUS
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
    /// >>> Helicity.PLUS.value
    /// 1
    ///
    #[getter]
    fn value(&self) -> i8 {
        self.inner.integer()
    }

    /// Convert this helicity to its integer value.
    ///
    /// Examples
    /// --------
    /// >>> int(Helicity.MINUS)
    /// -1
    ///
    fn __int__(&self) -> i8 {
        self.inner.integer()
    }

    /// Compare two helicity values.
    ///
    /// Examples
    /// --------
    /// >>> Helicity(0) == Helicity.ZERO
    /// True
    ///
    /// Parameters
    /// ----------
    /// other : Helicity
    ///     Helicity to compare with this value.
    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    /// Return an evaluable-style representation of this helicity.
    ///
    /// Examples
    /// --------
    /// >>> print(f"selected external helicity: {Helicity.PLUS!r}")
    ///
    fn __repr__(&self) -> String {
        format!("Helicity({})", self.inner.integer())
    }
}

/// A Cartesian axis used to specify spatial rotations.
///
/// Examples
/// --------
/// >>> import symbolica.community.feynkit as fk
/// >>> beam_axis = fk.Axis.Z
/// >>> rotation = fk.Rotation.quarter_turn(beam_axis)
///
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
/// >>> import symbolica.community.feynkit as fk
/// >>> algorithm = fk.JetAlgorithm.AntiKt
/// >>> definition = fk.JetDefinition(algorithm, radius=0.4)
///
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
/// >>> import symbolica.community.feynkit as fk
/// >>> p = fk.ThreeMomentum(30.0, 40.0, 10.0)
/// >>> p.pt
/// 50.0
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
    /// >>> momentum = ThreeMomentum(3.0, 4.0, 0.0)
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
    #[getter]
    fn px(&self) -> f64 {
        self.inner.px
    }
    /// Return the y component.
    #[getter]
    fn py(&self) -> f64 {
        self.inner.py
    }
    /// Return the z component.
    #[getter]
    fn pz(&self) -> f64 {
        self.inner.pz
    }
    /// Return the squared Euclidean norm.
    ///
    /// Examples
    /// --------
    /// >>> ThreeMomentum(3.0, 4.0, 0.0).norm_squared
    /// 25.0
    ///
    #[getter]
    fn norm_squared(&self) -> f64 {
        self.inner.norm_squared()
    }
    /// Return the Euclidean norm.
    ///
    /// Examples
    /// --------
    /// >>> ThreeMomentum(3.0, 4.0, 0.0).norm
    /// 5.0
    ///
    #[getter]
    fn norm(&self) -> f64 {
        self.inner.norm()
    }
    /// Return the transverse-momentum magnitude.
    ///
    /// Examples
    /// --------
    /// >>> ThreeMomentum(3.0, 4.0, 12.0).pt
    /// 5.0
    ///
    #[getter]
    fn pt(&self) -> f64 {
        self.inner.pt()
    }
    /// Return the azimuthal angle in radians.
    ///
    /// Examples
    /// --------
    /// >>> ThreeMomentum(1.0, 0.0, 0.0).phi
    /// 0.0
    ///
    #[getter]
    fn phi(&self) -> f64 {
        self.inner.phi()
    }
    /// Return the pseudorapidity derived from the momentum direction.
    ///
    /// Examples
    /// --------
    /// >>> ThreeMomentum(1.0, 0.0, 0.0).pseudorapidity
    /// 0.0
    ///
    #[getter]
    fn pseudorapidity(&self) -> f64 {
        self.inner.pseudorapidity()
    }
    /// Return the Euclidean dot product with another three-momentum.
    ///
    /// Examples
    /// --------
    /// >>> ThreeMomentum(1.0, 0.0, 0.0).dot(ThreeMomentum(2.0, 0.0, 0.0))
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
    /// >>> ThreeMomentum(1.0, 0.0, 0.0).cross(ThreeMomentum(0.0, 1.0, 0.0))
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
    /// >>> transfer = incoming - outgoing
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
    /// >>> outgoing_convention = -incoming_convention
    ///
    fn __neg__(&self) -> Self {
        (-self.inner).into()
    }

    /// Scale every spatial momentum component.
    ///
    /// Examples
    /// --------
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
    /// >>> ThreeMomentum(3.0, 4.0, 0.0).on_shell().energy
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
    /// >>> track_momentum = ThreeMomentum(1.0, 2.0, 3.0)
    /// >>> print(f"track momentum: {track_momentum!r}")
    ///
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
    /// Leave ``momentum`` as the final expression in a notebook cell to render
    /// its Cartesian components.
    ///
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
    /// IPython invokes this method when only a text representation is supported.
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
/// >>> import symbolica.community.feynkit as fk
/// >>> p = fk.FourMomentum(50.0, 30.0, 40.0, 0.0)
/// >>> p.mass_squared
/// 0.0
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
    /// >>> momentum = FourMomentum(5.0, 3.0, 4.0, 0.0)
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
    #[getter]
    fn energy(&self) -> f64 {
        self.inner.temporal.value
    }
    /// Return the x component of spatial momentum.
    #[getter]
    fn px(&self) -> f64 {
        self.inner.spatial.px
    }
    /// Return the y component of spatial momentum.
    #[getter]
    fn py(&self) -> f64 {
        self.inner.spatial.py
    }
    /// Return the z component of spatial momentum.
    #[getter]
    fn pz(&self) -> f64 {
        self.inner.spatial.pz
    }
    /// Return the spatial three-momentum.
    ///
    /// Examples
    /// --------
    /// >>> FourMomentum(5.0, 3.0, 4.0, 0.0).spatial.pt
    /// 5.0
    ///
    #[getter]
    fn spatial(&self) -> PyThreeMomentum {
        self.inner.spatial.into()
    }
    /// Return ``(energy, px, py, pz)``.
    ///
    /// Examples
    /// --------
    /// >>> FourMomentum(5.0, 3.0, 4.0, 0.0).components()
    /// (5.0, 3.0, 4.0, 0.0)
    ///
    fn components(&self) -> (f64, f64, f64, f64) {
        self.inner.into()
    }
    /// Return the Minkowski dot product with another four-momentum.
    ///
    /// Examples
    /// --------
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
    /// >>> FourMomentum(5.0, 3.0, 4.0, 0.0).mass_squared
    /// 0.0
    ///
    #[getter]
    fn mass_squared(&self) -> f64 {
        self.inner.mass_squared()
    }
    /// Return the invariant mass.
    ///
    /// Examples
    /// --------
    /// >>> FourMomentum(5.0, 0.0, 0.0, 0.0).mass
    /// 5.0
    ///
    #[getter]
    fn mass(&self) -> f64 {
        self.inner.mass()
    }
    /// Return the transverse-momentum magnitude.
    ///
    /// Examples
    /// --------
    /// >>> FourMomentum(13.0, 3.0, 4.0, 12.0).pt
    /// 5.0
    ///
    #[getter]
    fn pt(&self) -> f64 {
        self.inner.pt()
    }
    /// Return the azimuthal angle in radians.
    ///
    /// Examples
    /// --------
    /// >>> FourMomentum(1.0, 1.0, 0.0, 0.0).phi
    /// 0.0
    ///
    #[getter]
    fn phi(&self) -> f64 {
        self.inner.phi()
    }
    /// Return the pseudorapidity of the spatial momentum.
    ///
    /// Examples
    /// --------
    /// >>> FourMomentum(1.0, 1.0, 0.0, 0.0).pseudorapidity
    /// 0.0
    ///
    #[getter]
    fn pseudorapidity(&self) -> f64 {
        self.inner.pseudorapidity()
    }
    /// Return the longitudinal rapidity.
    ///
    /// Examples
    /// --------
    /// >>> FourMomentum(1.0, 1.0, 0.0, 0.0).rapidity
    /// 0.0
    ///
    #[getter]
    fn rapidity(&self) -> f64 {
        self.inner.rapidity()
    }
    /// Return the wrapped azimuthal separation from another momentum.
    ///
    /// Examples
    /// --------
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
    /// >>> outgoing_convention = -incoming_convention
    ///
    fn __neg__(&self) -> Self {
        (-self.inner).into()
    }

    /// Scale every four-momentum component.
    ///
    /// Examples
    /// --------
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
    /// >>> muon_momentum = FourMomentum(5.0, 3.0, 4.0, 0.0)
    /// >>> print(f"muon four-momentum: {muon_momentum!r}")
    ///
    fn __repr__(&self) -> String {
        let (energy, px, py, pz): (f64, f64, f64, f64) = self.inner.into();
        format!("FourMomentum({energy}, {px}, {py}, {pz})")
    }

    /// Render the momentum as a contravariant four-vector.
    ///
    /// Examples
    /// --------
    /// Leave ``momentum`` as the final expression in a notebook cell to render
    /// its energy and Cartesian momentum components.
    ///
    fn _repr_latex_(&self) -> String {
        let (energy, px, py, pz): (f64, f64, f64, f64) = self.inner.into();
        format!(r"$p^\mu=\left({energy},{px},{py},{pz}\right)$")
    }

    /// Write the constructor-style form to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// IPython invokes this method when only a text representation is supported.
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
/// >>> import symbolica.community.feynkit as fk
/// >>> rotation = fk.Rotation.quarter_turn(fk.Axis.Z)
/// >>> rotated = rotation.apply_three(fk.ThreeMomentum(1.0, 0.0, 0.0))
///
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
    /// >>> rotation = Rotation.identity()
    ///
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
    /// >>> rotation = Rotation.euler(0.1, 0.2, 0.3)
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
    /// >>> rotation = Rotation.quarter_turn(Axis.Z)
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
    /// >>> rotated = Rotation.quarter_turn(Axis.Z).apply_three(momentum)
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
    /// >>> rotated = Rotation.quarter_turn(Axis.Z).apply_four(momentum)
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
/// >>> import symbolica.community.feynkit as fk
/// >>> boost = fk.Boost(fk.ThreeMomentum(0.0, 0.0, 0.5))
/// >>> boosted = boost.apply(fk.FourMomentum(10.0, 0.0, 0.0, 0.0))
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
    /// >>> boost = Boost(ThreeMomentum(0.0, 0.0, 0.5))
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
    #[getter]
    fn beta(&self) -> PyThreeMomentum {
        self.inner.beta().to_owned().into()
    }
    /// Apply this Lorentz boost to a four-momentum.
    ///
    /// Examples
    /// --------
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
/// >>> jets = fk.JetDefinition.anti_kt(0.4).cluster(particles).jets
/// >>> leading_jet = jets[0]
/// >>> leading_jet.momentum
///
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
    #[getter]
    fn momentum(&self) -> PyFourMomentum {
        self.inner.momentum.into()
    }
    /// Return sorted positions of the input momenta assigned to this jet.
    #[getter]
    fn constituent_indices(&self) -> Vec<usize> {
        self.inner.constituent_indices().to_vec()
    }
    /// Return the jet transverse momentum.
    #[getter]
    fn pt(&self) -> f64 {
        self.inner.pt()
    }
    /// Return the jet rapidity.
    #[getter]
    fn rapidity(&self) -> f64 {
        self.inner.rapidity()
    }
    /// Return the jet azimuthal angle in radians.
    #[getter]
    fn phi(&self) -> f64 {
        self.inner.phi()
    }

    /// Return a concise jet summary with its kinematics and constituent count.
    ///
    /// Examples
    /// --------
    /// >>> print(jet)
    ///
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
    /// Leave ``jet`` as the final expression in a notebook cell to display its
    /// transverse momentum, rapidity, azimuth, and constituents.
    ///
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
    /// IPython invokes this method when only a text representation is supported.
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
/// >>> definition = fk.JetDefinition.anti_kt(0.4, minimum_pt=20.0)
/// >>> clustering = definition.cluster(particles)
/// >>> jets = clustering.jets
///
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
#[pymethods]
impl PyClusteringResult {
    /// Return inclusive jets ordered by decreasing transverse momentum.
    ///
    /// Examples
    /// --------
    /// >>> jets = result.jets
    /// >>> all(left.pt >= right.pt for left, right in zip(jets, jets[1:]))
    /// True
    ///
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
    /// >>> jet_multiplicity = len(result)
    ///
    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Return one jet by decreasing-transverse-momentum index.
    ///
    /// Examples
    /// --------
    /// >>> leading_jet = result[0]
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
    /// >>> transverse_momenta = [jet.pt for jet in result]
    ///
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
    /// >>> print(result)
    ///
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
    /// Leave ``result`` as the final expression in a notebook cell to display
    /// the jet collection.
    ///
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
    /// IPython invokes this method when only a text representation is supported.
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
/// >>> import symbolica.community.feynkit as fk
/// >>> definition = fk.JetDefinition.anti_kt(radius=0.4, minimum_pt=20.0)
/// >>> result = definition.cluster(particles)
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
    /// >>> definition = JetDefinition(JetAlgorithm.AntiKt, 0.4, 20.0)
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
    /// >>> definition = JetDefinition.kt(0.4)
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
    /// >>> definition = JetDefinition.cambridge_aachen(0.4)
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
    /// >>> definition = JetDefinition.anti_kt(0.4, minimum_pt=20.0)
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
    #[getter]
    fn algorithm(&self) -> PyJetAlgorithm {
        self.inner.algorithm().into()
    }
    /// Return the jet-radius parameter.
    #[getter]
    fn radius(&self) -> f64 {
        *self.inner.radius()
    }
    /// Return the minimum transverse momentum for retained jets.
    #[getter]
    fn minimum_pt(&self) -> f64 {
        *self.inner.minimum_pt()
    }

    /// Cluster four-momenta with this sequential-recombination definition.
    ///
    /// Examples
    /// --------
    /// >>> result = JetDefinition.anti_kt(0.4).cluster(momenta)
    /// >>> leading_jet = result.jets[0]
    /// >>> leading_jet.momentum
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
