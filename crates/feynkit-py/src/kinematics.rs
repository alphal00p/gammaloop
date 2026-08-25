use feynkit_kinematics::{
    Axis, Boost, ClusteringResult, FourMomentum, Helicity, Jet, JetAlgorithm, JetDefinition,
    Rotation, ThreeMomentum,
};
use pyo3::{prelude::*, types::PyModule};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

use crate::error;

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
    #[new]
    fn new(value: i8) -> PyResult<Self> {
        Helicity::try_from(value)
            .map(|inner| Self { inner })
            .map_err(error::kinematics)
    }

    #[staticmethod]
    fn parse(value: &str) -> PyResult<Self> {
        value
            .parse()
            .map(|inner| Self { inner })
            .map_err(error::kinematics)
    }

    #[classattr]
    #[pyo3(name = "MINUS")]
    fn minus() -> PyHelicity {
        Self {
            inner: Helicity::MINUS,
        }
    }

    #[classattr]
    #[pyo3(name = "ZERO")]
    fn zero() -> PyHelicity {
        Self {
            inner: Helicity::ZERO,
        }
    }

    #[classattr]
    #[pyo3(name = "PLUS")]
    fn plus() -> PyHelicity {
        Self {
            inner: Helicity::PLUS,
        }
    }

    #[getter]
    fn value(&self) -> i8 {
        self.inner.integer()
    }

    fn __int__(&self) -> i8 {
        self.inner.integer()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    fn __repr__(&self) -> String {
        format!("Helicity({})", self.inner.integer())
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    name = "Axis",
    module = "symbolica.community.feynkit",
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

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    name = "JetAlgorithm",
    module = "symbolica.community.feynkit",
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
    #[new]
    fn new(px: f64, py: f64, pz: f64) -> Self {
        ThreeMomentum::new(px, py, pz).into()
    }

    #[getter]
    fn px(&self) -> f64 {
        self.inner.px
    }
    #[getter]
    fn py(&self) -> f64 {
        self.inner.py
    }
    #[getter]
    fn pz(&self) -> f64 {
        self.inner.pz
    }
    fn norm_squared(&self) -> f64 {
        self.inner.norm_squared()
    }
    fn norm(&self) -> f64 {
        self.inner.norm()
    }
    fn pt(&self) -> f64 {
        self.inner.pt()
    }
    fn phi(&self) -> f64 {
        self.inner.phi()
    }
    fn pseudorapidity(&self) -> f64 {
        self.inner.pseudorapidity()
    }
    fn dot(&self, other: &Self) -> f64 {
        self.inner.dot(&other.inner)
    }
    fn cross(&self, other: &Self) -> Self {
        self.inner.cross(&other.inner).into()
    }
    fn delta_phi(&self, other: &Self) -> f64 {
        self.inner.delta_phi(&other.inner)
    }
    fn delta_r(&self, other: &Self) -> f64 {
        self.inner.delta_r(&other.inner)
    }

    #[pyo3(signature = (mass=None))]
    fn on_shell(&self, mass: Option<f64>) -> PyFourMomentum {
        self.inner.on_shell(mass.as_ref()).into()
    }

    fn __repr__(&self) -> String {
        format!(
            "ThreeMomentum({}, {}, {})",
            self.inner.px, self.inner.py, self.inner.pz
        )
    }
}

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
    #[new]
    fn new(energy: f64, px: f64, py: f64, pz: f64) -> Self {
        FourMomentum::from_args(energy, px, py, pz).into()
    }

    #[getter]
    fn energy(&self) -> f64 {
        self.inner.temporal.value
    }
    #[getter]
    fn px(&self) -> f64 {
        self.inner.spatial.px
    }
    #[getter]
    fn py(&self) -> f64 {
        self.inner.spatial.py
    }
    #[getter]
    fn pz(&self) -> f64 {
        self.inner.spatial.pz
    }
    #[getter]
    fn spatial(&self) -> PyThreeMomentum {
        self.inner.spatial.into()
    }
    fn components(&self) -> (f64, f64, f64, f64) {
        self.inner.into()
    }
    fn dot(&self, other: &Self) -> f64 {
        self.inner.dot(&other.inner)
    }
    fn mass_squared(&self) -> f64 {
        self.inner.mass_squared()
    }
    fn mass(&self) -> f64 {
        self.inner.mass()
    }
    fn pt(&self) -> f64 {
        self.inner.pt()
    }
    fn phi(&self) -> f64 {
        self.inner.phi()
    }
    fn pseudorapidity(&self) -> f64 {
        self.inner.pseudorapidity()
    }
    fn rapidity(&self) -> f64 {
        self.inner.rapidity()
    }
    fn delta_phi(&self, other: &Self) -> f64 {
        self.inner.delta_phi(&other.inner)
    }
    fn delta_r(&self, other: &Self) -> f64 {
        self.inner.delta_r(&other.inner)
    }
    fn __add__(&self, other: &Self) -> Self {
        (self.inner + other.inner).into()
    }
    fn __sub__(&self, other: &Self) -> Self {
        (self.inner - other.inner).into()
    }

    fn __repr__(&self) -> String {
        let (energy, px, py, pz): (f64, f64, f64, f64) = self.inner.into();
        format!("FourMomentum({energy}, {px}, {py}, {pz})")
    }
}

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
    #[staticmethod]
    fn identity() -> Self {
        Self {
            inner: Rotation::Identity,
        }
    }

    #[staticmethod]
    fn euler(alpha: f64, beta: f64, gamma: f64) -> Self {
        Self {
            inner: Rotation::euler(alpha, beta, gamma),
        }
    }

    #[staticmethod]
    fn quarter_turn(axis: PyAxis) -> Self {
        Self {
            inner: Rotation::quarter_turn(axis.into()),
        }
    }

    fn apply_three(&self, momentum: &PyThreeMomentum) -> PyThreeMomentum {
        self.inner.rotate_three(&momentum.inner).into()
    }

    fn apply_four(&self, momentum: &PyFourMomentum) -> PyFourMomentum {
        self.inner.rotate_four(&momentum.inner).into()
    }

    fn apply_inverse_three(&self, momentum: &PyThreeMomentum) -> PyThreeMomentum {
        self.inner.inverse_rotate_three(&momentum.inner).into()
    }

    fn apply_inverse_four(&self, momentum: &PyFourMomentum) -> PyFourMomentum {
        self.inner.inverse_rotate_four(&momentum.inner).into()
    }
}

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
    #[new]
    fn new(beta: &PyThreeMomentum) -> PyResult<Self> {
        Boost::new(beta.inner)
            .map(|inner| Self { inner })
            .map_err(error::kinematics)
    }

    #[getter]
    fn beta(&self) -> PyThreeMomentum {
        self.inner.beta().to_owned().into()
    }
    fn apply(&self, momentum: &PyFourMomentum) -> PyFourMomentum {
        self.inner.apply(&momentum.inner).into()
    }
    fn apply_inverse(&self, momentum: &PyFourMomentum) -> PyFourMomentum {
        self.inner.apply_inverse(&momentum.inner).into()
    }
}

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
    #[getter]
    fn momentum(&self) -> PyFourMomentum {
        self.inner.momentum.into()
    }
    #[getter]
    fn constituent_indices(&self) -> Vec<usize> {
        self.inner.constituent_indices().to_vec()
    }
    #[getter]
    fn pt(&self) -> f64 {
        self.inner.pt()
    }
    #[getter]
    fn rapidity(&self) -> f64 {
        self.inner.rapidity()
    }
    #[getter]
    fn phi(&self) -> f64 {
        self.inner.phi()
    }
}

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
    #[getter]
    fn jets(&self) -> Vec<PyJet> {
        self.inner
            .jets
            .iter()
            .cloned()
            .map(|inner| PyJet { inner })
            .collect()
    }
    fn __len__(&self) -> usize {
        self.inner.len()
    }
}

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
    #[new]
    #[pyo3(signature = (algorithm, radius, minimum_pt=0.0))]
    fn new(algorithm: PyJetAlgorithm, radius: f64, minimum_pt: f64) -> Self {
        Self {
            inner: JetDefinition::new(algorithm.into(), radius).with_minimum_pt(minimum_pt),
        }
    }

    #[staticmethod]
    #[pyo3(signature = (radius, minimum_pt=0.0))]
    fn kt(radius: f64, minimum_pt: f64) -> Self {
        Self {
            inner: JetDefinition::kt(radius).with_minimum_pt(minimum_pt),
        }
    }

    #[staticmethod]
    #[pyo3(signature = (radius, minimum_pt=0.0))]
    fn cambridge_aachen(radius: f64, minimum_pt: f64) -> Self {
        Self {
            inner: JetDefinition::cambridge_aachen(radius).with_minimum_pt(minimum_pt),
        }
    }

    #[staticmethod]
    #[pyo3(signature = (radius, minimum_pt=0.0))]
    fn anti_kt(radius: f64, minimum_pt: f64) -> Self {
        Self {
            inner: JetDefinition::anti_kt(radius).with_minimum_pt(minimum_pt),
        }
    }

    #[getter]
    fn algorithm(&self) -> PyJetAlgorithm {
        self.inner.algorithm().into()
    }
    #[getter]
    fn radius(&self) -> f64 {
        *self.inner.radius()
    }
    #[getter]
    fn minimum_pt(&self) -> f64 {
        *self.inner.minimum_pt()
    }

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
