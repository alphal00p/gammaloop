use std::{
    collections::{BTreeMap, VecDeque},
    sync::{Arc, Mutex},
    time::{Duration, Instant},
};

use feynkit_generator::{
    CancellationToken, DiagramGroup, EdgeColor, FilterScope, GenerationControl, GenerationFilter,
    GenerationOptions, GenerationProgress, GenerationReport, GenerationResult, GenerationType,
    Generator, GraphGroupingOptions, GroupMember, NodeColor, NumeratorGrouping, ParticleSelector,
    Process, SelfEnergyFilterOptions, SewnFilterOptions, SnailFilterOptions, TadpoleFilterOptions,
    VertexSelector,
};
use feynkit_graph::{DiagramId, EdgeId};
use feynkit_model::Model;
use pyo3::{
    FromPyObject, IntoPyObjectExt,
    exceptions::{PyIndexError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyAny, PyBool, PyEllipsis, PyList, PyModule},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    PyStubType, TypeInfo,
    derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods},
};

use crate::{
    display::render_diagram_html,
    error,
    graph::PyFeynmanDiagram,
    model::{PyModel, PyParticle, PyVertexRule},
};
use symbolica::{
    api::python::{PythonExpression, PythonGraph},
    graph::Graph,
};

/// Select whether diagrams describe an amplitude or a squared cross section.
///
/// The generation type determines which graph construction and external-state
/// conventions are applied to a particle-physics process.
///
/// Examples
/// --------
/// >>> import symbolica.community.feynkit as fk
/// >>> kind = fk.GenerationType.AMPLITUDE
/// >>> kind.value
/// 'amplitude'
///
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    name = "GenerationType",
    module = "symbolica.community.feynkit",
    rename_all = "SCREAMING_SNAKE_CASE",
    frozen,
    from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PyGenerationType {
    Amplitude,
    CrossSection,
}

impl From<GenerationType> for PyGenerationType {
    fn from(value: GenerationType) -> Self {
        match value {
            GenerationType::Amplitude => Self::Amplitude,
            GenerationType::CrossSection => Self::CrossSection,
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyGenerationType {
    /// Return the stable lowercase value used to identify the generation type.
    ///
    /// Examples
    /// --------
    /// >>> fk.GenerationType.AMPLITUDE.value
    /// 'amplitude'
    ///
    #[getter]
    fn value(&self) -> &'static str {
        match self {
            Self::Amplitude => "amplitude",
            Self::CrossSection => "cross_section",
        }
    }

    /// Format the generation type as its stable lowercase value.
    ///
    /// Examples
    /// --------
    /// >>> str(fk.GenerationType.CROSS_SECTION)
    /// 'cross_section'
    ///
    fn __str__(&self) -> &'static str {
        self.value()
    }

    /// Compare with another generation type or its lowercase string value.
    ///
    /// Examples
    /// --------
    /// >>> fk.GenerationType.AMPLITUDE == "amplitude"
    /// True
    ///
    /// Parameters
    /// ----------
    /// other : object
    ///     Generation type or string to compare with.
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        if let Ok(other) = other.cast::<Self>() {
            self == other.get()
        } else {
            other
                .extract::<String>()
                .is_ok_and(|other| self.value() == other)
        }
    }
}

/// A model-independent way to identify an external particle.
///
/// Selectors may use a UFO particle name or a signed PDG code, making process
/// definitions convenient while deferring validation to a concrete model.
///
/// Examples
/// --------
/// >>> import symbolica.community.feynkit as fk
/// >>> electron = fk.ParticleSelector.by_pdg(11)
/// >>> positron = fk.ParticleSelector.by_name("e+")
///
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ParticleSelector",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PyParticleSelector {
    inner: ParticleSelector,
}

impl From<ParticleSelector> for PyParticleSelector {
    fn from(inner: ParticleSelector) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyParticleSelector {
    /// Select a particle by its model name.
    ///
    /// Examples
    /// --------
    /// >>> fk.ParticleSelector.by_name("e-").name
    /// 'e-'
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Particle name as it appears in the model.
    #[staticmethod]
    fn by_name(name: String) -> Self {
        Self {
            inner: ParticleSelector::Name(name),
        }
    }

    /// Select a particle by its signed PDG code.
    ///
    /// Examples
    /// --------
    /// >>> fk.ParticleSelector.by_pdg(11).pdg
    /// 11
    ///
    /// Parameters
    /// ----------
    /// pdg : int
    ///     Signed PDG code of the particle.
    #[staticmethod]
    fn by_pdg(pdg: i64) -> Self {
        Self {
            inner: ParticleSelector::Pdg(pdg),
        }
    }

    /// Return the selected particle name.
    ///
    /// Raises :class:`TypeError` for a PDG selector. Use :attr:`is_name` to
    /// distinguish the selector variants before accessing their payloads.
    #[getter]
    fn name(&self) -> PyResult<&str> {
        match &self.inner {
            ParticleSelector::Id { .. } => Err(PyTypeError::new_err(
                "an ID particle selector does not carry a name",
            )),
            ParticleSelector::Name(name) => Ok(name),
            ParticleSelector::Pdg(_) => Err(PyTypeError::new_err(
                "a PDG particle selector does not carry a name",
            )),
        }
    }

    /// Return the selected PDG code.
    ///
    /// Raises :class:`TypeError` for a name selector. Use :attr:`is_pdg` to
    /// distinguish the selector variants before accessing their payloads.
    #[getter]
    fn pdg(&self) -> PyResult<i64> {
        match &self.inner {
            ParticleSelector::Id { .. } | ParticleSelector::Name(_) => Err(PyTypeError::new_err(
                "this particle selector does not carry a PDG code",
            )),
            ParticleSelector::Pdg(pdg) => Ok(*pdg),
        }
    }

    /// Report whether this selector identifies a particle by name.
    #[getter]
    fn is_name(&self) -> bool {
        matches!(&self.inner, ParticleSelector::Name(_))
    }

    /// Report whether this selector identifies a particle by PDG code.
    #[getter]
    fn is_pdg(&self) -> bool {
        matches!(&self.inner, ParticleSelector::Pdg(_))
    }

    /// Format the selector as its particle name or signed PDG code.
    ///
    /// Examples
    /// --------
    /// >>> str(fk.ParticleSelector.by_pdg(11))
    /// '11'
    ///
    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    /// Return a constructor-style representation of the selector.
    ///
    /// Examples
    /// --------
    /// >>> print(fk.ParticleSelector.by_pdg(11))  # electron selector in a process
    ///
    fn __repr__(&self) -> String {
        match &self.inner {
            ParticleSelector::Id { particle, model } => {
                format!("ParticleSelector(<particle#{}@{model}>)", particle.index())
            }
            ParticleSelector::Name(name) => format!("ParticleSelector.by_name('{name}')"),
            ParticleSelector::Pdg(pdg) => format!("ParticleSelector.by_pdg({pdg})"),
        }
    }

    /// Compare with another selector, a particle name, or a PDG code.
    ///
    /// Examples
    /// --------
    /// >>> fk.ParticleSelector.by_pdg(11) == 11
    /// True
    ///
    /// Parameters
    /// ----------
    /// other : object
    ///     Selector, particle-name string, or integer PDG code to compare with.
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        if let Ok(other) = other.cast::<Self>() {
            self == other.get()
        } else if let Ok(name) = other.extract::<String>() {
            self.inner == ParticleSelector::Name(name)
        } else {
            other
                .extract::<i64>()
                .is_ok_and(|pdg| self.inner == ParticleSelector::Pdg(pdg))
        }
    }
}

#[derive(Clone)]
pub(crate) struct ParticleInput(ParticleSelector);

impl<'a, 'py> FromPyObject<'a, 'py> for ParticleInput {
    type Error = PyErr;

    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(particle) = value.cast::<PyParticle>() {
            Ok(Self(ParticleSelector::Pdg(particle.get().signed_pdg())))
        } else if let Ok(pdg) = value.extract::<i64>() {
            Ok(Self(ParticleSelector::Pdg(pdg)))
        } else {
            value
                .extract::<String>()
                .map(ParticleSelector::Name)
                .map(Self)
        }
    }
}

#[derive(FromPyObject)]
pub(crate) enum SelectorInput {
    Particle(ParticleInput),
    Selector(PyParticleSelector),
    Name(String),
    Pdg(i64),
}

#[derive(FromPyObject)]
pub(crate) enum VertexInput {
    Vertex(PyVertexRule),
    Name(String),
}

// PyO3 can extract the concrete pyclass here, but does not implement
// `FromPyObject` for `Box<T>`, so boxing the large variant would break Python input.
#[allow(clippy::large_enum_variant)]
#[derive(FromPyObject)]
pub(crate) enum DiagramSelectionInput {
    Diagram(PyFeynmanDiagram),
    Text(String),
}

impl DiagramSelectionInput {
    fn split(self, ids: &mut Vec<DiagramId>, names: &mut Vec<String>) {
        match self {
            Self::Diagram(diagram) => ids.push(diagram.inner.id()),
            Self::Text(text) => match text.parse() {
                Ok(id) => ids.push(id),
                Err(_) => names.push(text),
            },
        }
    }
}

pub(crate) struct OrderRangeInput<const UNBOUNDED: bool = false> {
    pub(crate) minimum: usize,
    pub(crate) maximum: Option<usize>,
}

impl<const UNBOUNDED: bool> Default for OrderRangeInput<UNBOUNDED> {
    fn default() -> Self {
        Self {
            minimum: 0,
            maximum: Some(0),
        }
    }
}

impl<'py, const UNBOUNDED: bool> IntoPyObject<'py> for OrderRangeInput<UNBOUNDED> {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        if Some(self.minimum) == self.maximum {
            self.minimum.into_bound_py_any(py)
        } else {
            (self.minimum, self.maximum).into_bound_py_any(py)
        }
    }
}

impl<'a, 'py, const UNBOUNDED: bool> FromPyObject<'a, 'py> for OrderRangeInput<UNBOUNDED> {
    type Error = PyErr;

    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if value.is_instance_of::<PyBool>() {
            return Err(PyTypeError::new_err(
                "order must be a non-negative integer or a (minimum, maximum) pair",
            ));
        }
        let (minimum, maximum) = if let Ok(exact) = value.extract::<usize>() {
            (exact, Some(exact))
        } else if let Ok(range) = value.extract::<(usize, Option<usize>)>() {
            range
        } else {
            return Err(PyTypeError::new_err(
                "order must be a non-negative integer or a (minimum, maximum) pair",
            ));
        };
        if !UNBOUNDED && maximum.is_none() {
            return Err(PyValueError::new_err(
                "this order requires a finite maximum",
            ));
        }
        if maximum.is_some_and(|maximum| maximum < minimum) {
            return Err(PyValueError::new_err(
                "order bounds must satisfy minimum <= maximum",
            ));
        }
        Ok(Self { minimum, maximum })
    }
}

impl From<SelectorInput> for ParticleSelector {
    fn from(value: SelectorInput) -> Self {
        match value {
            SelectorInput::Particle(particle) => particle.0,
            SelectorInput::Selector(selector) => selector.inner,
            SelectorInput::Name(name) => Self::Name(name),
            SelectorInput::Pdg(pdg) => Self::Pdg(pdg),
        }
    }
}

impl From<VertexInput> for VertexSelector {
    fn from(value: VertexInput) -> Self {
        match value {
            VertexInput::Vertex(vertex) => Self::Name(vertex.vertex_rule_name().to_owned()),
            VertexInput::Name(name) => Self::Name(name),
        }
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for SelectorInput {
    fn type_input() -> TypeInfo {
        PyParticle::type_input()
            | PyParticleSelector::type_input()
            | String::type_input()
            | i64::type_input()
    }

    fn type_output() -> TypeInfo {
        PyParticle::type_output()
            | PyParticleSelector::type_output()
            | String::type_output()
            | i64::type_output()
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for ParticleInput {
    fn type_input() -> TypeInfo {
        PyParticle::type_input() | String::type_input() | i64::type_input()
    }

    fn type_output() -> TypeInfo {
        PyParticle::type_output() | String::type_output() | i64::type_output()
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for VertexInput {
    fn type_input() -> TypeInfo {
        PyVertexRule::type_input() | String::type_input()
    }

    fn type_output() -> TypeInfo {
        PyVertexRule::type_output() | String::type_output()
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for DiagramSelectionInput {
    fn type_input() -> TypeInfo {
        PyFeynmanDiagram::type_input() | String::type_input()
    }

    fn type_output() -> TypeInfo {
        PyFeynmanDiagram::type_output() | String::type_output()
    }
}

#[cfg(feature = "python_stubgen")]
impl<const UNBOUNDED: bool> PyStubType for OrderRangeInput<UNBOUNDED> {
    fn type_input() -> TypeInfo {
        usize::type_input()
            | if UNBOUNDED {
                <(usize, Option<usize>)>::type_input()
            } else {
                <(usize, usize)>::type_input()
            }
    }

    fn type_output() -> TypeInfo {
        Self::type_input()
    }
}

/// A scattering or decay process to pass to the diagram generator.
///
/// A process records its incoming and outgoing particles, loop-order range,
/// and optional external-state symmetrizations. External states accept loaded
/// :class:`Particle` objects as well as names, signed PDG codes, and explicit
/// :class:`ParticleSelector` objects.
///
/// Examples
/// --------
/// >>> import symbolica.community.feynkit as fk
/// >>> process = fk.Process.amplitude(["e-", "e+"], ["mu-", "mu+"])
/// >>> one_loop = process.with_loop_count(1, 1)
///
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Process",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyProcess {
    pub(crate) inner: Process,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyProcess {
    /// Define an amplitude with ordered incoming and outgoing particles.
    ///
    /// Examples
    /// --------
    /// >>> electron = model.particle_by_pdg(11)
    /// >>> positron = model.particle_by_pdg(-11)
    /// >>> process = fk.Process.amplitude([electron, positron], [22])
    ///
    /// Parameters
    /// ----------
    /// incoming : sequence[Particle | ParticleSelector | str | int]
    ///     Incoming particles in external-leg order.
    /// outgoing : sequence[Particle | ParticleSelector | str | int]
    ///     Outgoing particles in external-leg order.
    #[staticmethod]
    fn amplitude(incoming: Vec<SelectorInput>, outgoing: Vec<SelectorInput>) -> Self {
        Self {
            inner: Process::amplitude(
                incoming.into_iter().map(ParticleSelector::from),
                outgoing.into_iter().map(ParticleSelector::from),
            ),
        }
    }

    /// Define a cross section with ordered incoming and outgoing particles.
    ///
    /// Examples
    /// --------
    /// >>> process = fk.Process.cross_section([11, -11], [13, -13])
    ///
    /// Parameters
    /// ----------
    /// incoming : sequence[Particle | ParticleSelector | str | int]
    ///     Incoming particles in external-leg order.
    /// outgoing : sequence[Particle | ParticleSelector | str | int]
    ///     Outgoing particles in external-leg order.
    #[staticmethod]
    fn cross_section(incoming: Vec<SelectorInput>, outgoing: Vec<SelectorInput>) -> Self {
        Self {
            inner: Process::cross_section(
                incoming.into_iter().map(ParticleSelector::from),
                outgoing.into_iter().map(ParticleSelector::from),
            )
            .symmetrize_final(true),
        }
    }

    /// Return a process restricted to an inclusive loop-count range.
    ///
    /// Examples
    /// --------
    /// >>> loop_process = process.with_loop_count(1, 2)
    ///
    /// Parameters
    /// ----------
    /// minimum : int
    ///     Minimum number of loops to generate.
    /// maximum : int
    ///     Maximum number of loops to generate, inclusive.
    fn with_loop_count(&self, minimum: usize, maximum: usize) -> PyResult<Self> {
        self.inner
            .clone()
            .with_loop_count(minimum, maximum)
            .map(|inner| Self { inner })
            .map_err(error::process)
    }

    /// Return a process accepting any of the supplied final states.
    /// Amplitudes accept exactly one alternative; cross sections may include an empty
    /// alternative for a vacuum final state.
    ///
    /// Examples
    /// --------
    /// >>> process = fk.Process.cross_section([11, -11], [22, 22])
    /// >>> inclusive = process.with_final_state_alternatives([[22, 22], [13, -13]])
    ///
    /// Parameters
    /// ----------
    /// alternatives : sequence[sequence[Particle | ParticleSelector | str | int]]
    ///     Allowed outgoing particle lists.
    fn with_final_state_alternatives(
        &self,
        alternatives: Vec<Vec<SelectorInput>>,
    ) -> PyResult<Self> {
        self.inner
            .clone()
            .with_final_state_alternatives(
                alternatives
                    .into_iter()
                    .map(|state| state.into_iter().map(ParticleSelector::from)),
            )
            .map(|inner| Self { inner })
            .map_err(error::process)
    }

    /// Return a process configured with the selected graph symmetries.
    ///
    /// Examples
    /// --------
    /// >>> process = process.with_symmetrization(initial=True, final_state=True)
    ///
    /// Parameters
    /// ----------
    /// initial : bool, optional
    ///     Identify graphs related by permutations of initial-state particles.
    /// final_state : bool, optional
    ///     Identify graphs related by permutations of final-state particles.
    /// left_right : bool, optional
    ///     Identify cross-section graphs related by exchanging amplitude sides.
    /// external_fermions : bool, optional
    ///     Include amplitude fermions in enabled external-state symmetry classes.
    #[pyo3(signature = (*, initial=false, final_state=false, left_right=false, external_fermions=false))]
    fn with_symmetrization(
        &self,
        initial: bool,
        final_state: bool,
        left_right: bool,
        external_fermions: bool,
    ) -> Self {
        Self {
            inner: self
                .inner
                .clone()
                .symmetrize_initial(initial)
                .symmetrize_final(final_state)
                .symmetrize_left_right(left_right)
                .symmetrize_external_fermions(external_fermions),
        }
    }

    /// Return whether this process generates amplitudes or cross sections.
    ///
    /// Examples
    /// --------
    /// >>> process.generation_type == fk.GenerationType.AMPLITUDE
    /// True
    ///
    #[getter]
    fn generation_type(&self) -> PyGenerationType {
        self.inner.generation_type().into()
    }

    /// Return the ordered incoming-particle selectors.
    ///
    /// Examples
    /// --------
    /// >>> [selector.pdg for selector in process.incoming]
    /// [11, -11]
    ///
    #[getter]
    fn incoming(&self) -> Vec<PyParticleSelector> {
        self.inner
            .incoming()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Return every allowed ordered final-state alternative.
    #[getter]
    fn outgoing_alternatives(&self) -> Vec<Vec<PyParticleSelector>> {
        self.inner
            .outgoing_alternatives()
            .iter()
            .map(|state| state.iter().cloned().map(Into::into).collect())
            .collect()
    }

    /// Return the inclusive minimum and maximum loop counts.
    ///
    /// Examples
    /// --------
    /// >>> process.with_loop_count(1, 2).loop_count
    /// (1, 2)
    ///
    #[getter]
    fn loop_count(&self) -> (usize, usize) {
        let range = self.inner.loop_count();
        (*range.start(), *range.end())
    }

    /// Report whether initial-state permutations are identified.
    ///
    /// Examples
    /// --------
    /// >>> process.with_symmetrization(initial=True).symmetrizes_initial
    /// True
    ///
    #[getter]
    fn symmetrizes_initial(&self) -> bool {
        self.inner.symmetrizes_initial()
    }

    /// Report whether final-state permutations are identified.
    ///
    /// Examples
    /// --------
    /// >>> process.with_symmetrization(final_state=True).symmetrizes_final
    /// True
    ///
    #[getter]
    fn symmetrizes_final(&self) -> bool {
        self.inner.symmetrizes_final()
    }

    /// Report whether exchanging the two cross-section sides is identified.
    ///
    /// Examples
    /// --------
    /// >>> process.with_symmetrization(left_right=True).symmetrizes_left_right
    /// True
    ///
    #[getter]
    fn symmetrizes_left_right(&self) -> bool {
        self.inner.symmetrizes_left_right()
    }

    /// Report whether amplitude fermions participate in enabled state symmetries.
    ///
    /// Examples
    /// --------
    /// >>> process.with_symmetrization(external_fermions=True).symmetrizes_external_fermions
    /// True
    ///
    #[getter]
    fn symmetrizes_external_fermions(&self) -> bool {
        self.inner.symmetrizes_external_fermions()
    }
}

/// A thread-safe signal for cancelling a long diagram-generation job.
///
/// Pass one token as ``cancellation_token`` and call ``cancel`` from a
/// controlling thread when a large topology search should stop early.
///
/// Examples
/// --------
/// >>> import symbolica.community.feynkit as fk
/// >>> token = fk.CancellationToken()
/// >>> token.is_cancelled
/// False
///
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CancellationToken",
    module = "symbolica.community.feynkit",
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyCancellationToken {
    inner: CancellationToken,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCancellationToken {
    /// Create an independent token that can cancel a generation request.
    ///
    /// Examples
    /// --------
    /// >>> token = fk.CancellationToken()
    ///
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// Mark this token as cancelled for all generation requests using it.
    ///
    /// Examples
    /// --------
    /// >>> token.cancel()
    ///
    fn cancel(&self) {
        self.inner.cancel();
    }

    /// Report whether cancellation has been requested.
    ///
    /// Examples
    /// --------
    /// >>> token = fk.CancellationToken()
    /// >>> token.cancel()
    /// >>> token.is_cancelled
    /// True
    ///
    #[getter]
    fn is_cancelled(&self) -> bool {
        self.inner.is_cancelled()
    }
}

/// Configure rejection of self-energy subgraphs by mass category.
///
/// Examples
/// --------
/// >>> fk.SelfEnergyFilterOptions(veto_massive=True, veto_massless=True)
///
/// Parameters
/// ----------
/// veto_massive : bool, optional
///     Reject self energies carried by massive particles.
/// veto_massless : bool, optional
///     Reject self energies carried by massless particles.
/// only_scaleless : bool, optional
///     Currently unsupported; ``True`` makes generation return an error.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "SelfEnergyFilterOptions",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PySelfEnergyFilterOptions {
    inner: SelfEnergyFilterOptions,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySelfEnergyFilterOptions {
    /// Configure rejection of self-energy subgraphs by mass category.
    ///
    /// Examples
    /// --------
    /// >>> fk.SelfEnergyFilterOptions(veto_massive=True, veto_massless=True)
    ///
    /// Parameters
    /// ----------
    /// veto_massive : bool, optional
    ///     Reject self energies carried by massive particles.
    /// veto_massless : bool, optional
    ///     Reject self energies carried by massless particles.
    /// only_scaleless : bool, optional
    ///     Currently unsupported; ``True`` makes generation return an error.
    #[new]
    #[pyo3(signature = (*, veto_massive=true, veto_massless=true, only_scaleless=false))]
    fn new(veto_massive: bool, veto_massless: bool, only_scaleless: bool) -> Self {
        Self {
            inner: SelfEnergyFilterOptions {
                veto_massive,
                veto_massless,
                only_scaleless,
            },
        }
    }
}

/// Configure rejection of tadpoles by the mass of their attachment.
///
/// Examples
/// --------
/// >>> fk.TadpoleFilterOptions(veto_attached_to_massless=True)
///
/// Parameters
/// ----------
/// veto_attached_to_massive : bool, optional
///     Reject tadpoles attached through a massive particle.
/// veto_attached_to_massless : bool, optional
///     Reject tadpoles attached through a massless particle.
/// only_scaleless : bool, optional
///     Currently unsupported; ``True`` makes generation return an error.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "TadpoleFilterOptions",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyTadpoleFilterOptions {
    inner: TadpoleFilterOptions,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyTadpoleFilterOptions {
    /// Configure rejection of tadpoles by the mass of their attachment.
    ///
    /// Examples
    /// --------
    /// >>> fk.TadpoleFilterOptions(veto_attached_to_massless=True)
    ///
    /// Parameters
    /// ----------
    /// veto_attached_to_massive : bool, optional
    ///     Reject tadpoles attached through a massive particle.
    /// veto_attached_to_massless : bool, optional
    ///     Reject tadpoles attached through a massless particle.
    /// only_scaleless : bool, optional
    ///     Currently unsupported; ``True`` makes generation return an error.
    #[new]
    #[pyo3(signature = (*, veto_attached_to_massive=true, veto_attached_to_massless=true, only_scaleless=false))]
    fn new(
        veto_attached_to_massive: bool,
        veto_attached_to_massless: bool,
        only_scaleless: bool,
    ) -> Self {
        Self {
            inner: TadpoleFilterOptions {
                veto_attached_to_massive,
                veto_attached_to_massless,
                only_scaleless,
            },
        }
    }
}

/// Configure rejection of zero-momentum snail subgraphs.
///
/// Examples
/// --------
/// >>> fk.SnailFilterOptions(veto_attached_to_massless=True)
///
/// Parameters
/// ----------
/// veto_attached_to_massive : bool, optional
///     Reject zero-momentum snails attached through a massive particle.
/// veto_attached_to_massless : bool, optional
///     Reject zero-momentum snails attached through a massless particle.
/// only_scaleless : bool, optional
///     Currently unsupported; ``True`` makes generation return an error.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "SnailFilterOptions",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PySnailFilterOptions {
    inner: SnailFilterOptions,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySnailFilterOptions {
    /// Configure rejection of zero-momentum snail subgraphs.
    ///
    /// Examples
    /// --------
    /// >>> fk.SnailFilterOptions(veto_attached_to_massless=True)
    ///
    /// Parameters
    /// ----------
    /// veto_attached_to_massive : bool, optional
    ///     Reject zero-momentum snails attached through a massive particle.
    /// veto_attached_to_massless : bool, optional
    ///     Reject zero-momentum snails attached through a massless particle.
    /// only_scaleless : bool, optional
    ///     Currently unsupported; ``True`` makes generation return an error.
    #[new]
    #[pyo3(signature = (*, veto_attached_to_massive=false, veto_attached_to_massless=true, only_scaleless=false))]
    fn new(
        veto_attached_to_massive: bool,
        veto_attached_to_massless: bool,
        only_scaleless: bool,
    ) -> Self {
        Self {
            inner: SnailFilterOptions {
                veto_attached_to_massive,
                veto_attached_to_massless,
                only_scaleless,
            },
        }
    }
}

/// Choose numerator zero detection and cross-diagram grouping.
///
/// Examples
/// --------
/// >>> fk.NumeratorGrouping("identical", number_of_numerical_samples=7)
///
/// Parameters
/// ----------
/// mode : {"none", "zeroes", "identical", "up_to_sign", "up_to_scalar"}
///     Disable parsing/grouping, detect only zeroes, or compare numerators
///     exactly, up to a sign, or up to a scalar factor.
/// numerical_sample_seed : int, optional
///     Deterministic seed used to choose numerical substitution values.
/// number_of_numerical_samples : int, optional
///     Number of independent substitutions used to compare numerators.
/// differentiate_particle_masses_only : bool, optional
///     Treat internal species with equal mass and spin as interchangeable.
/// fully_numerical_substitution : bool, optional
///     Substitute scalar parameters as well as nonscalar indeterminates.
/// check_canonical_numerator : bool, optional
///     Try an exact canonical comparison before numerical sampling.
/// symmetric_polarizations : bool, optional
///     Reuse wavefunction samples across the two sides of a sewn external state.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "NumeratorGrouping",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyNumeratorGrouping {
    inner: NumeratorGrouping,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyNumeratorGrouping {
    /// Choose numerator zero detection and cross-diagram grouping.
    ///
    /// Examples
    /// --------
    /// >>> fk.NumeratorGrouping("identical", number_of_numerical_samples=7)
    ///
    /// Parameters
    /// ----------
    /// mode : {"none", "zeroes", "identical", "up_to_sign", "up_to_scalar"}
    ///     Disable parsing/grouping, detect only zeroes, or compare numerators
    ///     exactly, up to a sign, or up to a scalar factor.
    /// numerical_sample_seed : int, optional
    ///     Deterministic seed used to choose numerical substitution values.
    /// number_of_numerical_samples : int, optional
    ///     Number of independent substitutions used to compare numerators.
    /// differentiate_particle_masses_only : bool, optional
    ///     Treat internal species with equal mass and spin as interchangeable.
    /// fully_numerical_substitution : bool, optional
    ///     Substitute scalar parameters as well as nonscalar indeterminates.
    /// check_canonical_numerator : bool, optional
    ///     Try an exact canonical comparison before numerical sampling.
    /// symmetric_polarizations : bool, optional
    ///     Reuse wavefunction samples across the two sides of a sewn external state.
    #[new]
    #[pyo3(signature = (mode, *, numerical_sample_seed=3, number_of_numerical_samples=5, differentiate_particle_masses_only=true, fully_numerical_substitution=false, check_canonical_numerator=false, symmetric_polarizations=false))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        mode: &str,
        numerical_sample_seed: u16,
        number_of_numerical_samples: usize,
        differentiate_particle_masses_only: bool,
        fully_numerical_substitution: bool,
        check_canonical_numerator: bool,
        symmetric_polarizations: bool,
    ) -> PyResult<Self> {
        let options = GraphGroupingOptions {
            numerical_sample_seed,
            number_of_numerical_samples,
            differentiate_particle_masses_only,
            fully_numerical_substitution,
            check_canonical_numerator,
            symmetric_polarizations,
        };
        let inner = match mode {
            "none" => NumeratorGrouping::None,
            "zeroes" => NumeratorGrouping::OnlyDetectZeroes,
            "identical" => NumeratorGrouping::Identical(options),
            "up_to_sign" => NumeratorGrouping::UpToSign(options),
            "up_to_scalar" => NumeratorGrouping::UpToScalar(options),
            _ => {
                return Err(PyValueError::new_err(
                    "grouping mode must be 'none', 'zeroes', 'identical', 'up_to_sign', or 'up_to_scalar'",
                ));
            }
        };
        Ok(Self { inner })
    }
}

/// Configuration for Feynman-diagram generation and filtering.
///
/// The Python entry points construct the same Rust options from keyword arguments;
/// resource limits, filters, grouping, and cancellation have one implementation.
pub(crate) struct GenerationSettings {
    pub(crate) inner: GenerationOptions,
}

impl GenerationSettings {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        py: Python<'_>,
        process: &Process,
        threads: Option<usize>,
        max_vertices: Option<usize>,
        allow_self_loops: bool,
        allow_zero_flow_edges: bool,
        graph_prefix: Option<String>,
        particle_veto: Option<Vec<ParticleInput>>,
        vertex_allow: Option<Vec<VertexInput>>,
        vertex_veto: Option<Vec<VertexInput>>,
        maximum_bridges: Option<usize>,
        self_energy: Option<Py<PyAny>>,
        tadpoles: Option<Py<PyAny>>,
        zero_snails: Option<Py<PyAny>>,
        coupling_orders: Option<BTreeMap<String, OrderRangeInput<true>>>,
        fermion_loop_count_range: Option<(usize, usize)>,
        factorized_loop_topologies_count_range: Option<Py<PyAny>>,
        blob_range: Option<Py<PyAny>>,
        spectator_range: Option<Py<PyAny>>,
        perturbative_orders: Option<BTreeMap<String, usize>>,
        sewn_tadpoles: Option<bool>,
        cut_amplitude_coupling_orders: Option<BTreeMap<String, OrderRangeInput<true>>>,
        cut_amplitude_loop_count_range: Option<(usize, usize)>,
        select_diagrams: Option<Vec<DiagramSelectionInput>>,
        veto_diagrams: Option<Vec<DiagramSelectionInput>>,
        loop_momentum_bases: Option<Vec<(DiagramSelectionInput, Vec<usize>)>>,
        numerator_prefactor: Option<PythonExpression>,
        projector: Option<PythonExpression>,
        numerator_grouping: Option<Py<PyAny>>,
        cancellation_token: Option<PyCancellationToken>,
    ) -> PyResult<Self> {
        // Ported from GammaLoop's CLI policy. Explicit None disables a default;
        // Ellipsis selects the process-dependent default without a mutable preset.
        let cross_section = process.generation_type() == GenerationType::CrossSection;
        let vacuum = process.incoming().is_empty()
            && (cross_section || process.outgoing_alternatives().iter().all(Vec::is_empty));
        let self_energy = match self_energy {
            Some(value) if value.bind(py).is_instance_of::<PyEllipsis>() => {
                (!vacuum).then(SelfEnergyFilterOptions::default)
            }
            value => value
                .map(|value| {
                    value
                        .extract::<PySelfEnergyFilterOptions>(py)
                        .map_err(PyErr::from)
                        .map(|value| value.inner)
                })
                .transpose()?,
        };
        let tadpoles = match tadpoles {
            Some(value) if value.bind(py).is_instance_of::<PyEllipsis>() => {
                (!vacuum).then(TadpoleFilterOptions::default)
            }
            value => value
                .map(|value| {
                    value
                        .extract::<PyTadpoleFilterOptions>(py)
                        .map_err(PyErr::from)
                        .map(|value| value.inner)
                })
                .transpose()?,
        };
        let zero_snails = match zero_snails {
            Some(value) if value.bind(py).is_instance_of::<PyEllipsis>() => {
                (!vacuum).then(SnailFilterOptions::default)
            }
            value => value
                .map(|value| {
                    value
                        .extract::<PySnailFilterOptions>(py)
                        .map_err(PyErr::from)
                        .map(|value| value.inner)
                })
                .transpose()?,
        };
        let numerator_grouping = match numerator_grouping {
            Some(value) if value.bind(py).is_instance_of::<PyEllipsis>() => Some(
                NumeratorGrouping::UpToScalar(GraphGroupingOptions::default()),
            ),
            value => value
                .map(|value| {
                    value
                        .extract::<PyNumeratorGrouping>(py)
                        .map_err(PyErr::from)
                        .map(|value| value.inner)
                })
                .transpose()?,
        };
        let factorized_loop_topologies_count_range = match factorized_loop_topologies_count_range {
            Some(value) if value.bind(py).is_instance_of::<PyEllipsis>() => {
                vacuum.then_some((1, 1))
            }
            value => value
                .map(|value| value.extract::<(usize, usize)>(py))
                .transpose()?,
        };
        let blob_range = match blob_range {
            Some(value) if value.bind(py).is_instance_of::<PyEllipsis>() => {
                cross_section.then_some((1, 1))
            }
            value => value
                .map(|value| value.extract::<(usize, usize)>(py))
                .transpose()?,
        };
        let spectator_range = match spectator_range {
            Some(value) if value.bind(py).is_instance_of::<PyEllipsis>() => {
                cross_section.then_some((0, 0))
            }
            value => value
                .map(|value| value.extract::<(usize, usize)>(py))
                .transpose()?,
        };
        let mut inner = GenerationOptions::default()
            .allow_self_loops(allow_self_loops)
            .allow_zero_flow_edges(allow_zero_flow_edges);
        if let Some(value) = threads {
            inner = inner.threads(value);
        }
        if let Some(value) = max_vertices {
            inner = inner.max_vertices(value);
        }
        if let Some(value) = graph_prefix {
            inner = inner.graph_prefix(value);
        }
        if let Some(value) = particle_veto {
            inner = inner.with_graph_filter(GenerationFilter::ParticleVeto(
                value.into_iter().map(|particle| particle.0).collect(),
            ));
        }
        if let Some(value) = vertex_allow {
            inner = inner.with_graph_filter(GenerationFilter::VertexAllow(
                value.into_iter().map(Into::into).collect(),
            ));
        }
        if let Some(value) = vertex_veto {
            inner = inner.with_graph_filter(GenerationFilter::VertexVeto(
                value.into_iter().map(Into::into).collect(),
            ));
        }
        if let Some(value) = maximum_bridges {
            inner = inner.with_graph_filter(GenerationFilter::MaxNumberOfBridges(value));
        }
        if let Some(value) = self_energy {
            inner = inner.with_graph_filter(GenerationFilter::SelfEnergy(value));
        }
        if let Some(value) = tadpoles {
            inner = inner.with_graph_filter(GenerationFilter::Tadpoles(value));
        }
        if let Some(value) = zero_snails {
            inner = inner.with_graph_filter(GenerationFilter::ZeroSnails(value));
        }
        for (scope, orders) in [
            (FilterScope::Graph, coupling_orders),
            (FilterScope::CutAmplitude, cut_amplitude_coupling_orders),
        ] {
            if let Some(orders) = orders {
                let orders = orders
                    .into_iter()
                    .map(|(name, range)| (name, (range.minimum, range.maximum)))
                    .collect();
                inner = inner.with_scoped_filter(scope, GenerationFilter::CouplingOrders(orders));
            }
        }
        if let Some(value) = fermion_loop_count_range {
            inner = inner.with_graph_filter(GenerationFilter::FermionLoopCountRange(value));
        }
        if let Some(value) = factorized_loop_topologies_count_range {
            inner = inner
                .with_graph_filter(GenerationFilter::FactorizedLoopTopologiesCountRange(value));
        }
        if let Some(value) = blob_range {
            inner = inner.with_graph_filter(GenerationFilter::BlobRange(value.0..=value.1));
        }
        if let Some(value) = spectator_range {
            inner = inner.with_graph_filter(GenerationFilter::SpectatorRange(value.0..=value.1));
        }
        if let Some(value) = perturbative_orders {
            inner = inner.with_graph_filter(GenerationFilter::PerturbativeOrders(value));
        }
        if let Some(value) = sewn_tadpoles {
            inner = inner.with_graph_filter(GenerationFilter::Sewn(SewnFilterOptions {
                filter_tadpoles: value,
            }));
        }
        if let Some(value) = cut_amplitude_loop_count_range {
            inner = inner.with_cut_amplitude_filter(GenerationFilter::LoopCountRange(value));
        }
        if let Some(diagrams) = select_diagrams {
            let (mut ids, mut names) = (Vec::new(), Vec::new());
            for diagram in diagrams {
                diagram.split(&mut ids, &mut names);
            }
            if !ids.is_empty() {
                inner = inner.select_diagram_ids(ids);
            }
            if !names.is_empty() {
                inner = inner.select_diagrams(names);
            }
        }
        if let Some(diagrams) = veto_diagrams {
            let (mut ids, mut names) = (Vec::new(), Vec::new());
            for diagram in diagrams {
                diagram.split(&mut ids, &mut names);
            }
            inner = inner.veto_diagram_ids(ids).veto_diagrams(names);
        }
        for (diagram, edges) in loop_momentum_bases.unwrap_or_default() {
            let (mut ids, mut names) = (Vec::new(), Vec::new());
            diagram.split(&mut ids, &mut names);
            let edges = edges.into_iter().map(EdgeId).collect::<Vec<_>>();
            inner = if let Some(id) = ids.into_iter().next() {
                inner.with_loop_momentum_basis(id, edges)
            } else {
                inner.with_named_loop_momentum_basis(names.pop().expect("one selector"), edges)
            };
        }
        if let Some(value) = numerator_prefactor {
            inner = inner.numerator_prefactor(value.expr);
        }
        if let Some(value) = projector {
            inner = inner.projector(value.expr);
        }
        if let Some(value) = numerator_grouping {
            inner = inner.numerator_grouping(value);
        }
        if let Some(value) = cancellation_token {
            inner = inner.cancellation_token(value.inner);
        }
        Ok(Self { inner })
    }
}

/// A progress snapshot delivered on the Python thread running generation.
///
/// Counts restart at each stage and measure processed work, not retained diagrams.
/// ``total`` is None when the amount of work is not yet known.
///
/// Examples
/// --------
/// >>> def report(progress):
/// ...     print(progress.stage, progress.completed, progress.total)
/// >>> result = generator.generate(process, progress=report)
///
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "GenerationProgress",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyGenerationProgress {
    inner: GenerationProgress,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyGenerationProgress {
    /// Pipeline stage: topologies, topology_filters, interactions,
    /// interaction_filters, numerators, selection, grouping, complete or cancelled.
    #[getter]
    fn stage(&self) -> &'static str {
        self.inner.stage
    }

    /// Work items processed within this stage.
    #[getter]
    fn completed(&self) -> usize {
        self.inner.completed
    }

    /// Stage total, or None while the amount of work is unknown.
    #[getter]
    fn total(&self) -> Option<usize> {
        self.inner.total
    }
}

/// Counts and completion status from a diagram-generation run.
///
/// The report distinguishes explored topologies and interaction assignments
/// from the physical diagrams retained after numerator checks.
///
/// Examples
/// --------
/// >>> report = result.report
/// >>> print(report)
///
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "GenerationReport",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyGenerationReport {
    inner: GenerationReport,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyGenerationReport {
    /// Return the number of distinct topologies considered during generation.
    #[getter]
    fn topology_count(&self) -> usize {
        self.inner.topology_count
    }
    /// Return the number of interaction assignments examined.
    #[getter]
    fn interaction_assignment_count(&self) -> usize {
        self.inner.interaction_assignment_count
    }
    /// Return the number of diagrams retained after removing zero numerators.
    #[getter]
    fn retained_count(&self) -> usize {
        self.inner.retained_count
    }
    /// Return the number of diagrams removed for having an exact zero numerator.
    #[getter]
    fn zero_numerator_count(&self) -> usize {
        self.inner.zero_numerator_count
    }
    /// Report whether generation finished without cancellation.
    #[getter]
    fn completed(&self) -> bool {
        self.inner.completed
    }

    /// Return a concise constructor-style generation summary.
    ///
    /// Examples
    /// --------
    /// >>> print(result.report)
    ///
    fn __repr__(&self) -> String {
        format!(
            "GenerationReport(topology_count={}, interaction_assignment_count={}, retained_count={}, zero_numerator_count={}, completed={})",
            self.inner.topology_count,
            self.inner.interaction_assignment_count,
            self.inner.retained_count,
            self.inner.zero_numerator_count,
            self.inner.completed,
        )
    }

    /// Render generation statistics as a compact HTML table.
    ///
    /// Examples
    /// --------
    /// Leave ``result.report`` as the final expression in a notebook cell.
    ///
    fn _repr_html_(&self) -> String {
        let status = if self.inner.completed {
            "completed"
        } else {
            "cancelled"
        };
        format!(
            "<div class=\"feynkit-generation-report\" style=\"display:inline-block;max-width:100%;overflow-x:auto\">\
             <strong>Generation report</strong>\
             <table style=\"border-collapse:collapse;margin-top:.25rem\"><tbody>\
             <tr><th style=\"padding:.2rem .65rem;text-align:left\">status</th><td style=\"padding:.2rem .65rem\">{status}</td></tr>\
             <tr><th style=\"padding:.2rem .65rem;text-align:left\">topologies</th><td style=\"padding:.2rem .65rem;text-align:right\">{}</td></tr>\
             <tr><th style=\"padding:.2rem .65rem;text-align:left\">interaction assignments</th><td style=\"padding:.2rem .65rem;text-align:right\">{}</td></tr>\
             <tr><th style=\"padding:.2rem .65rem;text-align:left\">retained diagrams</th><td style=\"padding:.2rem .65rem;text-align:right\">{}</td></tr>\
             <tr><th style=\"padding:.2rem .65rem;text-align:left\">zero numerators removed</th><td style=\"padding:.2rem .65rem;text-align:right\">{}</td></tr>\
             </tbody></table></div>",
            self.inner.topology_count,
            self.inner.interaction_assignment_count,
            self.inner.retained_count,
            self.inner.zero_numerator_count,
        )
    }

    /// Write a concise generation summary to an IPython pretty printer.
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

/// One diagram and its numerator ratio inside a grouped result.
///
/// A group member points into ``GenerationResult.diagrams`` and expresses its
/// numerator relative to the group's master diagram.
///
/// Examples
/// --------
/// >>> member = next(iter(next(iter(result.groups)).members))
/// >>> diagram = result.diagrams[member.diagram]
///
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "GroupMember",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyGroupMember {
    inner: GroupMember,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyGroupMember {
    /// Return the generated-order index from before zero-numerator removal.
    #[getter]
    fn source_diagram(&self) -> usize {
        self.inner.source_diagram
    }
    /// Return the source diagram's stable content-derived ID.
    #[getter]
    fn source_id(&self) -> String {
        self.inner.source_id.to_string()
    }
    /// Return the finalized display name assigned to the source diagram.
    #[getter]
    fn source_name(&self) -> &str {
        &self.inner.source_name
    }
    /// Return the collapsed master index in ``GenerationResult.diagrams``.
    ///
    /// Examples
    /// --------
    /// >>> result.diagrams[result.groups[0].members[0].diagram]
    ///
    #[getter]
    fn diagram(&self) -> usize {
        self.inner.diagram
    }
    /// Return the member numerator divided by the group master numerator.
    #[getter]
    fn ratio(&self) -> String {
        self.inner.ratio.to_plain_string()
    }
    /// Parse the numerator ratio as a native Symbolica expression.
    ///
    /// Examples
    /// --------
    /// >>> ratio = result.groups[0].members[0].ratio_expression()
    ///
    fn ratio_expression(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.ratio.clone(),
        }
    }
    /// Parse the source diagram's numerator-independent factor.
    ///
    /// Examples
    /// --------
    /// >>> group = result.groups[0]
    /// >>> member = group.members[0]
    /// >>> master = result.diagrams[group.master]
    /// >>> reconstructed = (member.overall_factor_expression()
    /// ...                  * member.ratio_expression()
    /// ...                  * master.numerator_expression())
    /// >>> reconstructed
    ///
    fn overall_factor_expression(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.overall_factor.clone(),
        }
    }
}

/// Diagrams whose numerators are related by known scalar factors.
///
/// Grouping lets amplitude calculations evaluate one master numerator and
/// reconstruct related diagrams from each member's symbolic ratio.
///
/// Examples
/// --------
/// >>> group = next(iter(result.groups))
/// >>> master_diagram = result.diagrams[group.master]
///
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "DiagramGroup",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyDiagramGroup {
    inner: DiagramGroup,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyDiagramGroup {
    /// Return the retained-diagram index used as the numerator reference.
    ///
    /// Examples
    /// --------
    /// >>> result.diagrams[result.groups[0].master]
    ///
    #[getter]
    fn master(&self) -> usize {
        self.inner.master
    }
    /// Return the deterministically ordered group members, including the master.
    #[getter]
    fn members(&self) -> Vec<PyGroupMember> {
        self.inner
            .members
            .iter()
            .cloned()
            .map(|inner| PyGroupMember { inner })
            .collect()
    }
}

/// Feynman diagrams and diagnostics produced for one process.
///
/// The result retains generated diagrams in deterministic order and may also
/// contain numerator-equivalence groups for efficient downstream evaluation.
///
/// Examples
/// --------
/// >>> import symbolica.community.feynkit as fk
/// >>> result = fk.Generator(model).generate(process)
/// >>> diagrams = result.diagrams
///
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "GenerationResult",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyGenerationResult {
    inner: GenerationResult,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyGenerationResult {
    /// Return every retained diagram in generated order.
    #[getter]
    fn diagrams(&self) -> Vec<PyFeynmanDiagram> {
        self.inner
            .diagrams
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Return the numerator groups in deterministic master-index order.
    #[getter]
    fn groups(&self) -> Vec<PyDiagramGroup> {
        self.inner
            .groups
            .iter()
            .cloned()
            .map(|inner| PyDiagramGroup { inner })
            .collect()
    }

    /// Return generation counts and completion status.
    #[getter]
    fn report(&self) -> PyGenerationReport {
        PyGenerationReport {
            inner: self.inner.report.clone(),
        }
    }

    /// Return the number of retained diagrams.
    ///
    /// Examples
    /// --------
    /// >>> number_of_diagrams = len(result)
    ///
    fn __len__(&self) -> usize {
        self.inner.diagrams.len()
    }

    /// Return one retained diagram by generated-order index.
    ///
    /// Examples
    /// --------
    /// >>> first_diagram = result[0]
    ///
    /// Parameters
    /// ----------
    /// index : int
    ///     Zero-based index; negative indices count from the end.
    fn __getitem__(&self, index: isize) -> PyResult<PyFeynmanDiagram> {
        let length = self.inner.diagrams.len() as isize;
        let index = if index < 0 { length + index } else { index };
        if !(0..length).contains(&index) {
            return Err(PyIndexError::new_err("diagram index out of range"));
        }
        Ok(self.inner.diagrams[index as usize].clone().into())
    }

    /// Iterate over retained diagrams in deterministic generated order.
    ///
    /// Examples
    /// --------
    /// >>> one_loop = [diagram for diagram in result if diagram.loop_count == 1]
    ///
    #[gen_stub(override_return_type(
        type_repr = "collections.abc.Iterator[FeynmanDiagram]",
        imports = ("collections.abc")
    ))]
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyList::new(py, self.diagrams())?.call_method0("__iter__")
    }

    /// Return a concise summary of the retained diagrams and groups.
    ///
    /// Examples
    /// --------
    /// >>> print(result)
    ///
    fn __repr__(&self) -> String {
        format!(
            "GenerationResult(diagrams={}, groups={}, completed={})",
            self.inner.diagrams.len(),
            self.inner.groups.len(),
            self.inner.report.completed,
        )
    }

    /// Render generation statistics and a bounded diagram gallery as HTML.
    ///
    /// At most six diagrams are rendered so that displaying a large generation
    /// result remains responsive. Access ``result.diagrams`` to inspect the
    /// complete collection.
    ///
    /// Examples
    /// --------
    /// Leave ``result`` as the final expression in a notebook cell.
    ///
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        const PREVIEW_LIMIT: usize = 6;

        let report = PyGenerationReport {
            inner: self.inner.report.clone(),
        }
        ._repr_html_();
        let diagrams = self
            .inner
            .diagrams
            .iter()
            .take(PREVIEW_LIMIT)
            .map(|diagram| {
                render_diagram_html(py, diagram).map(|html| {
                    format!("<div style=\"min-width:0;overflow-x:auto\">{}</div>", html)
                })
            })
            .collect::<PyResult<String>>()?;
        let gallery = if diagrams.is_empty() {
            "<p style=\"margin:.5rem 0;opacity:.75\">No diagrams retained.</p>".to_owned()
        } else {
            format!(
                "<div style=\"display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,20rem),1fr));gap:.75rem;margin-top:.5rem\">{diagrams}</div>"
            )
        };
        let remainder = self.inner.diagrams.len().saturating_sub(PREVIEW_LIMIT);
        let omitted = if remainder == 0 {
            String::new()
        } else {
            format!(
                "<p style=\"margin:.5rem 0 0;opacity:.75\">{remainder} additional diagram{} not shown.</p>",
                if remainder == 1 { "" } else { "s" },
            )
        };

        Ok(format!(
            "<section class=\"feynkit-generation-result\" style=\"max-width:100%\">\
             <h3 style=\"margin:.25rem 0\">Generation result</h3>{report}{gallery}{omitted}</section>"
        ))
    }

    /// Write a concise result summary to an IPython pretty printer.
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

/// Generate Feynman diagrams from a particle model and a process definition.
///
/// The generator combines model interactions with graph topologies at the
/// requested loop orders. It instantiates the selected Feynman rules as
/// numerator annotations on interaction vertices and propagator edges, and
/// stores their product as the diagram-wide numerator. The resulting typed
/// diagrams are ready for symbolic manipulation and CFF operations.
///
/// Examples
/// --------
/// >>> import symbolica.community.feynkit as fk
/// >>> generator = fk.Generator(model)
/// >>> result = generator.generate(process)
///
/// Parameters
/// ----------
/// model : Model
///     Particle model supplying fields, propagators, and interaction vertices.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Generator",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyGenerator {
    inner: Generator,
    model: Arc<Model>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyGenerator {
    /// Create a diagram generator backed by a loaded particle model.
    ///
    /// Examples
    /// --------
    /// >>> generator = fk.Generator(model)
    ///
    /// Parameters
    /// ----------
    /// model : Model
    ///     Particle model supplying particles, interactions, and parameters.
    #[new]
    fn new(model: &PyModel) -> Self {
        Self {
            inner: Generator::new(model.inner.clone()),
            model: model.inner.clone(),
        }
    }

    /// Return the particle model used by this generator.
    #[getter]
    fn model(&self) -> PyModel {
        self.model.clone().into()
    }

    /// Generate and optionally group all diagrams matching a process.
    ///
    /// Model Feynman rules are instantiated into each returned diagram: inspect
    /// ``vertex.numerator_expression()`` and ``edge.numerator_expression()`` for
    /// individual factors, or ``diagram.numerator_expression()`` for their
    /// combined numerator.
    ///
    /// Examples
    /// --------
    /// >>> result = generator.generate(process, max_vertices=6)
    /// >>> combined_numerator = result.diagrams[0].numerator_expression()
    ///
    /// Parameters
    /// ----------
    /// process : Process
    ///     Scattering process and loop range to generate.
    /// threads : int or None, optional
    ///     Number of worker threads; None uses the generator default.
    /// max_vertices : int or None, optional
    ///     Maximum interaction vertices; None applies no override.
    /// allow_self_loops : bool, optional
    ///     Permit propagators that start and end on the same vertex; defaults to True.
    /// allow_zero_flow_edges : bool, optional
    ///     Permit internal edges with identically zero momentum flow.
    /// graph_prefix : str or None, optional
    ///     Prefix assigned to generated diagram names.
    /// particle_veto : sequence[Particle | str | int] or None, optional
    ///     Reject graphs containing these particles, model names, or signed PDG codes.
    /// vertex_allow : sequence[VertexRule | str] or None, optional
    ///     Keep only graphs whose vertices use these model rules or names.
    /// vertex_veto : sequence[VertexRule | str] or None, optional
    ///     Reject graphs containing these interaction vertices.
    /// maximum_bridges : int or None, optional
    ///     Largest allowed number of internal graph bridges; defaults to 0.
    ///     Pass None to allow unrestricted bridges, including exchange-channel trees.
    /// self_energy : SelfEnergyFilterOptions or None, optional
    ///     Reject self-energy subgraphs. Omission enables the default filter for
    ///     non-vacuum processes; explicit None disables it. Ellipsis selects automatic defaults.
    /// tadpoles : TadpoleFilterOptions or None, optional
    ///     Reject tadpole subgraphs. Omission enables the default filter for
    ///     non-vacuum processes; explicit None disables it. Ellipsis selects automatic defaults.
    /// zero_snails : SnailFilterOptions or None, optional
    ///     Reject zero-snail subgraphs. Omission enables the default filter for
    ///     non-vacuum processes; explicit None disables it. Ellipsis selects automatic defaults.
    /// coupling_orders : dict[str, int | tuple[int, int or None]] or None, optional
    ///     Exact coupling powers or inclusive ranges; an upper None is unbounded.
    /// fermion_loop_count_range : tuple[int, int] or None, optional
    ///     Inclusive range of closed fermion loops.
    /// factorized_loop_topologies_count_range : tuple[int, int] or None, optional
    ///     Inclusive range of factorized loop-topology components. Defaults to
    ///     ``(1, 1)`` for vacuum processes; ``None`` disables the restriction.
    /// blob_range : tuple[int, int] or None, optional
    ///     Inclusive cross-section blob-count range. Defaults to ``(1, 1)`` for
    ///     cross sections; ``None`` disables the restriction.
    /// spectator_range : tuple[int, int] or None, optional
    ///     Inclusive cross-section spectator-count range. Defaults to ``(0, 0)`` for
    ///     cross sections; ``None`` disables the restriction.
    /// perturbative_orders : dict[str, int] or None, optional
    ///     Exact perturbative powers required for cross-section graphs.
    /// sewn_tadpoles : bool or None, optional
    ///     Reject tadpoles revealed by sewing cross-section sides; None applies no filter.
    /// cut_amplitude_coupling_orders : dict[str, int | tuple[int, int or None]] or None, optional
    ///     Coupling-order bounds applied independently within every cut amplitude.
    /// cut_amplitude_loop_count_range : tuple[int, int] or None, optional
    ///     Inclusive combined loop count across both sides of every cut.
    /// select_diagrams : sequence[FeynmanDiagram | str] or None, optional
    ///     Retain only these diagram objects, content-derived IDs, or finalized names.
    /// veto_diagrams : sequence[FeynmanDiagram | str] or None, optional
    ///     Remove these diagram objects, content-derived IDs, or finalized names.
    /// loop_momentum_bases : sequence[tuple[FeynmanDiagram | str, sequence[int]]] or None, optional
    ///     Diagram selectors paired with ordered stable edge IDs for independent loop momenta.
    /// numerator_prefactor : Expression or None, optional
    ///     Scalar multiplier retained on every finalized diagram numerator.
    /// projector : Expression or None, optional
    ///     Override external-state contraction; S("1") disables external wavefunctions.
    /// numerator_grouping : NumeratorGrouping or None, optional
    ///     Omission groups up to scalar rescaling, matching the GammaLoop CLI.
    ///     Explicit None disables comparison, but diagrams still contain numerators.
    /// progress : Callable[[GenerationProgress], None] or None, optional
    ///     Observe stage changes and coalesced counts on the calling Python thread.
    ///     Callback exceptions propagate and stop generation.
    /// filter : Callable[[symbolica.core.Graph, int], bool] or None, optional
    ///     Prune partial topologies during enumeration. The first N vertices are
    ///     complete. False rejects only this search branch. Edge data is the base
    ///     particle PDG code; node data is 0 internally, -(index+1) for incoming
    ///     legs and +(index+1) for outgoing legs. Mutating the snapshot does not
    ///     change enumeration. Keep callbacks cheap: each snapshot is constructed
    ///     using Symbolica's Python Graph API.
    /// cancellation_token : CancellationToken or None, optional
    ///     Shared token for cancelling a running generation task. Token cancellation
    ///     returns an incomplete result; Python signal-handler exceptions, including
    ///     KeyboardInterrupt, stop generation and propagate to the caller.
    #[pyo3(signature = (process, *, threads=None, max_vertices=None, allow_self_loops=true, allow_zero_flow_edges=false, graph_prefix=None, particle_veto=None, vertex_allow=None, vertex_veto=None, maximum_bridges=0, self_energy=Some(Python::attach(|py| py.Ellipsis())), tadpoles=Some(Python::attach(|py| py.Ellipsis())), zero_snails=Some(Python::attach(|py| py.Ellipsis())), coupling_orders=None, fermion_loop_count_range=None, factorized_loop_topologies_count_range=Some(Python::attach(|py| py.Ellipsis())), blob_range=Some(Python::attach(|py| py.Ellipsis())), spectator_range=Some(Python::attach(|py| py.Ellipsis())), perturbative_orders=None, sewn_tadpoles=None, cut_amplitude_coupling_orders=None, cut_amplitude_loop_count_range=None, select_diagrams=None, veto_diagrams=None, loop_momentum_bases=None, numerator_prefactor=None, projector=None, numerator_grouping=Some(Python::attach(|py| py.Ellipsis())), cancellation_token=None, progress=None, filter=None))]
    #[allow(clippy::too_many_arguments)]
    fn generate(
        &self,
        py: Python<'_>,
        process: &PyProcess,
        threads: Option<usize>,
        max_vertices: Option<usize>,
        allow_self_loops: bool,
        allow_zero_flow_edges: bool,
        graph_prefix: Option<String>,
        particle_veto: Option<Vec<ParticleInput>>,
        vertex_allow: Option<Vec<VertexInput>>,
        vertex_veto: Option<Vec<VertexInput>>,
        #[gen_stub(override_type(type_repr = "int | None"))] maximum_bridges: Option<usize>,
        #[gen_stub(override_type(type_repr = "SelfEnergyFilterOptions | types.EllipsisType | None", imports = ("types")))]
        self_energy: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "TadpoleFilterOptions | types.EllipsisType | None", imports = ("types")))]
        tadpoles: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "SnailFilterOptions | types.EllipsisType | None", imports = ("types")))]
        zero_snails: Option<Py<PyAny>>,
        coupling_orders: Option<BTreeMap<String, OrderRangeInput<true>>>,
        fermion_loop_count_range: Option<(usize, usize)>,
        #[gen_stub(override_type(type_repr = "tuple[int, int] | types.EllipsisType | None", imports = ("types")))]
        factorized_loop_topologies_count_range: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "tuple[int, int] | types.EllipsisType | None", imports = ("types")))]
        blob_range: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "tuple[int, int] | types.EllipsisType | None", imports = ("types")))]
        spectator_range: Option<Py<PyAny>>,
        perturbative_orders: Option<BTreeMap<String, usize>>,
        sewn_tadpoles: Option<bool>,
        cut_amplitude_coupling_orders: Option<BTreeMap<String, OrderRangeInput<true>>>,
        cut_amplitude_loop_count_range: Option<(usize, usize)>,
        select_diagrams: Option<Vec<DiagramSelectionInput>>,
        veto_diagrams: Option<Vec<DiagramSelectionInput>>,
        loop_momentum_bases: Option<Vec<(DiagramSelectionInput, Vec<usize>)>>,
        numerator_prefactor: Option<PythonExpression>,
        projector: Option<PythonExpression>,
        #[gen_stub(override_type(type_repr = "NumeratorGrouping | types.EllipsisType | None", imports = ("types")))]
        numerator_grouping: Option<Py<PyAny>>,
        cancellation_token: Option<PyCancellationToken>,
        #[gen_stub(override_type(type_repr = "collections.abc.Callable[[GenerationProgress], None] | None", imports = ("collections.abc")))]
        progress: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "collections.abc.Callable[[symbolica.core.Graph, int], bool] | None", imports = ("collections.abc", "symbolica.core")))]
        filter: Option<Py<PyAny>>,
    ) -> PyResult<PyGenerationResult> {
        let generator = self.inner.clone();
        let process = process.inner.clone();
        let options = GenerationSettings::new(
            py,
            &process,
            threads,
            max_vertices,
            allow_self_loops,
            allow_zero_flow_edges,
            graph_prefix,
            particle_veto,
            vertex_allow,
            vertex_veto,
            maximum_bridges,
            self_energy,
            tadpoles,
            zero_snails,
            coupling_orders,
            fermion_loop_count_range,
            factorized_loop_topologies_count_range,
            blob_range,
            spectator_range,
            perturbative_orders,
            sewn_tadpoles,
            cut_amplitude_coupling_orders,
            cut_amplitude_loop_count_range,
            select_diagrams,
            veto_diagrams,
            loop_momentum_bases,
            numerator_prefactor,
            projector,
            numerator_grouping,
            cancellation_token,
        )?
        .inner;
        generate_diagrams(py, generator, process, options, progress, filter)
    }
}

/// Run both entry points with callbacks and signals on the calling Python thread.
pub(crate) fn generate_diagrams(
    py: Python<'_>,
    generator: Generator,
    process: Process,
    mut options: GenerationOptions,
    progress: Option<Py<PyAny>>,
    filter: Option<Py<PyAny>>,
) -> PyResult<PyGenerationResult> {
    py.check_signals()?;
    for (name, callback) in [("progress", &progress), ("filter", &filter)] {
        if callback
            .as_ref()
            .is_some_and(|callback| !callback.bind(py).is_callable())
        {
            return Err(PyTypeError::new_err(format!("{name} must be callable")));
        }
    }
    let interruption = Arc::new(Mutex::new(None));
    let snapshots = Arc::new(Mutex::new(VecDeque::<GenerationProgress>::new()));
    let delivery = Mutex::new((Instant::now(), None));
    let pending = snapshots.clone();
    let has_progress = progress.is_some();
    let deliver_progress = Arc::new(move |force: bool| -> PyResult<()> {
        let mut pending = pending.lock().unwrap();
        let mut delivery = delivery.lock().unwrap();
        let changed = pending.front().is_some_and(|p| Some(p.stage) != delivery.1);
        let finished = pending.back().is_some_and(|p| p.total == Some(p.completed));
        if !force
            && !changed
            && !finished
            && pending.len() < 2
            && delivery.0.elapsed() < Duration::from_millis(150)
        {
            return Ok(());
        }
        let updates = std::mem::take(&mut *pending);
        if let Some(last) = updates.back() {
            *delivery = (Instant::now(), Some(last.stage));
        }
        drop(pending);
        drop(delivery);
        if let Some(callback) = &progress {
            Python::attach(|py| -> PyResult<()> {
                for inner in updates {
                    callback.call1(py, (PyGenerationProgress { inner },))?;
                }
                Ok(())
            })?;
        }
        Ok(())
    });
    if has_progress {
        #[cfg(target_arch = "wasm32")]
        let (deliver, error) = (deliver_progress.clone(), interruption.clone());
        options = options.progress(move |snapshot| {
            let mut pending = snapshots.lock().unwrap();
            if pending
                .back()
                .is_some_and(|last| last.stage == snapshot.stage)
            {
                *pending.back_mut().unwrap() = snapshot;
            } else {
                pending.push_back(snapshot);
            }
            drop(pending);
            #[cfg(target_arch = "wasm32")]
            if error.lock().unwrap().is_none() {
                let result = deliver(false);
                if let Err(exception) = result {
                    *error.lock().unwrap() = Some(exception);
                    return GenerationControl::Cancel;
                }
            }
            GenerationControl::Continue
        });
    }
    let has_filter = filter.is_some();
    let pdgs: Vec<_> = generator
        .model()
        .particles()
        .iter()
        .map(|p| p.pdg_code)
        .collect();
    let apply_filter =
        move |graph: &Graph<NodeColor, EdgeColor>, completed_vertices| -> PyResult<bool> {
            Python::attach(|py| {
                // Replace these Python API calls with PythonGraph::from once Symbolica
                // exposes its constructor. Keep the callback's Graph type unchanged.
                let snapshot = py.get_type::<PythonGraph>().call0()?;
                for node in graph.nodes() {
                    let label = node.data.external.as_ref().map_or(0, |external| {
                        let index = external.index as i64 + 1;
                        match external.state {
                            feynkit_graph::ExternalState::Incoming => -index,
                            feynkit_graph::ExternalState::Outgoing => index,
                        }
                    });
                    snapshot.call_method1("add_node", (label,))?;
                }
                for edge in graph.edges() {
                    snapshot.call_method1(
                        "add_edge",
                        (
                            edge.vertices.0,
                            edge.vertices.1,
                            edge.directed,
                            pdgs[edge.data.particle.index()],
                        ),
                    )?;
                }
                filter
                    .as_ref()
                    .unwrap()
                    .call1(py, (snapshot, completed_vertices))?
                    .extract(py)
            })
        };
    #[cfg(not(target_arch = "wasm32"))]
    let result = {
        let cancellation = CancellationToken::new();
        let check = cancellation.clone();
        options = options.cancellation_check(move || check.is_cancelled());
        let (requests, receiver) = std::sync::mpsc::sync_channel(1);
        if has_filter {
            options = options.filter(move |graph, completed_vertices| {
                let (reply, response) = std::sync::mpsc::sync_channel(1);
                requests
                    .send((graph.clone(), completed_vertices, reply))
                    .is_ok()
                    && response.recv().unwrap_or(false)
            });
        }
        let interruption = interruption.clone();
        let deliver_progress = deliver_progress.clone();
        py.detach(move || {
            std::thread::scope(|scope| {
                let worker = scope.spawn(|| generator.generate(&process, &options));
                // Python handles signals only on its main thread, including while
                // generation is inside Rayon's parallel interaction assignment.
                // Filter requests wake this thread immediately; never impose the
                // progress refresh interval on each enumeration decision.
                while !worker.is_finished() {
                    if interruption.lock().unwrap().is_none() {
                        let result = Python::attach(|py| py.check_signals())
                            .and_then(|()| deliver_progress(false));
                        if let Err(error) = result {
                            *interruption.lock().unwrap() = Some(error);
                            cancellation.cancel();
                        }
                    }
                    if let Ok((graph, completed_vertices, reply)) =
                        receiver.recv_timeout(Duration::from_millis(25))
                    {
                        let accepted = if interruption.lock().unwrap().is_some() {
                            false
                        } else {
                            match apply_filter(&graph, completed_vertices) {
                                Ok(accepted) => accepted,
                                Err(error) => {
                                    *interruption.lock().unwrap() = Some(error);
                                    cancellation.cancel();
                                    false
                                }
                            }
                        };
                        let _ = reply.send(accepted);
                    }
                }
                worker
                    .join()
                    .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
            })
        })
    };
    #[cfg(target_arch = "wasm32")]
    let result = {
        // Browser kernels execute generation on the calling thread. Poll the
        // existing cancellation hook there instead of creating a native thread.
        let signal_error = interruption.clone();
        options = options.cancellation_check(move || {
            if signal_error.lock().unwrap().is_none() {
                let error = Python::attach(|py| py.check_signals()).err();
                *signal_error.lock().unwrap() = error;
            }
            signal_error.lock().unwrap().is_some()
        });
        if has_filter {
            let error = interruption.clone();
            options = options.filter(move |graph, completed_vertices| {
                match apply_filter(graph, completed_vertices) {
                    Ok(accepted) => accepted,
                    Err(exception) => {
                        *error.lock().unwrap() = Some(exception);
                        false
                    }
                }
            });
        }
        py.detach(|| generator.generate(&process, &options))
    };
    if let Some(error) = interruption.lock().unwrap().take() {
        return Err(error);
    }
    deliver_progress(true)?;
    py.check_signals()?;
    result
        .map(|inner| PyGenerationResult { inner })
        .map_err(error::generation)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyGenerationType>()?;
    module.add_class::<PyParticleSelector>()?;
    module.add_class::<PyProcess>()?;
    module.add_class::<PyCancellationToken>()?;
    module.add_class::<PySelfEnergyFilterOptions>()?;
    module.add_class::<PyTadpoleFilterOptions>()?;
    module.add_class::<PySnailFilterOptions>()?;
    module.add_class::<PyNumeratorGrouping>()?;
    module.add_class::<PyGenerationProgress>()?;
    module.add_class::<PyGenerationReport>()?;
    module.add_class::<PyGroupMember>()?;
    module.add_class::<PyDiagramGroup>()?;
    module.add_class::<PyGenerationResult>()?;
    module.add_class::<PyGenerator>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use pyo3::types::PyDict;

    use super::*;

    #[test]
    fn generation_defaults_match_the_ported_cli_policy() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "symbolica.community.feynkit").unwrap();
            crate::initialize_feynkit(&module).unwrap();
            let locals = PyDict::new(py);
            locals.set_item("fk", &module).unwrap();
            locals
                .set_item(
                    "MODEL_JSON",
                    include_str!("../tests/fixtures/scalars_2p_3p.json"),
                )
                .unwrap();
            let code = CString::new(r#"
import inspect
model = fk.Model.from_json(MODEL_JSON)
generator = fk.Generator(model)
settings = dict(max_vertices=3, vertex_allow=["V_3_SCALAR_000"])
explicit = dict(
    maximum_bridges=0, allow_self_loops=True,
    self_energy=fk.SelfEnergyFilterOptions(), tadpoles=fk.TadpoleFilterOptions(),
    zero_snails=fk.SnailFilterOptions(), numerator_grouping=fk.NumeratorGrouping("up_to_scalar"),
)
for via_model in (True, False):
    def generate(incoming, outgoing, loops=0, **kwargs):
        if via_model:
            return model.generate_diagrams(incoming, outgoing, loops=loops, **settings, **kwargs)
        process = fk.Process.amplitude(incoming, outgoing).with_loop_count(loops, loops)
        return generator.generate(process, **settings, **kwargs)

    implicit = generate([1000], [1000, 1000], loops=1)
    configured = generate([1000], [1000, 1000], loops=1, **explicit)
    assert implicit.report.completed and len(implicit) > 0
    assert [d.to_json() for d in implicit.diagrams] == [d.to_json() for d in configured.diagrams]

    # A tree with two cubic vertices needs an internal bridge. None explicitly
    # restores unrestricted generation; omitting maximum_bridges rejects it.
    assert len(generate([1000, 1000], [1000, 1000], numerator_grouping=None)) == 0
    unrestricted = generate([1000, 1000], [1000, 1000], maximum_bridges=None, numerator_grouping=None)
    assert len(unrestricted) > 0
    assert len(unrestricted.groups) == len(unrestricted)
    assert all(len(group.members) == 1 for group in unrestricted.groups)

    vacuum = generate([], [], loops=2, numerator_grouping=None)
    configured_vacuum = generate([], [], loops=2, numerator_grouping=None,
        maximum_bridges=0, allow_self_loops=True, self_energy=None,
        tadpoles=None, zero_snails=None, factorized_loop_topologies_count_range=(1, 1))
    assert [d.to_json() for d in vacuum.diagrams] == [d.to_json() for d in configured_vacuum.diagrams]

assert fk.Process.cross_section([1000], [1000]).symmetrizes_final
for generate in (model.generate_diagrams, generator.generate):
    parameters = inspect.signature(generate).parameters
    assert parameters["maximum_bridges"].default == 0
    assert parameters["allow_self_loops"].default is True
    assert parameters["numerator_grouping"].default is Ellipsis
"#).unwrap();
            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
    }

    #[test]
    fn progress_and_partial_filters_use_the_calling_python_thread() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "symbolica.community.feynkit").unwrap();
            crate::initialize_feynkit(&module).unwrap();
            let locals = PyDict::new(py);
            locals.set_item("fk", &module).unwrap();
            locals
                .set_item("Graph", py.get_type::<PythonGraph>())
                .unwrap();
            locals
                .set_item(
                    "MODEL_JSON",
                    include_str!("../tests/fixtures/scalars_2p_3p.json"),
                )
                .unwrap();
            let code = CString::new(r#"
import threading

model = fk.Model.from_json(MODEL_JSON)
generator = fk.Generator(model)
process = fk.Process.amplitude([1000], [1000, 1000]).with_loop_count(1, 1)
settings = dict(max_vertices=3, threads=2, vertex_allow=["V_3_SCALAR_000"])
caller = threading.get_ident()

for via_model in (True, False):
    def generate(**callbacks):
        if via_model:
            return model.generate_diagrams([1000], [1000, 1000], loops=1, **settings, **callbacks)
        return generator.generate(process, **settings, **callbacks)

    baseline = generate()
    assert baseline.report.completed and len(baseline) > 0
    updates = []
    topologies = []

    def report(p):
        assert threading.get_ident() == caller
        assert isinstance(p, fk.GenerationProgress)
        assert p.completed >= 0
        assert p.total is None or p.completed <= p.total
        updates.append(p)

    def keep(graph, completed_vertices):
        assert threading.get_ident() == caller
        assert isinstance(graph, Graph)
        assert 0 <= completed_vertices <= graph.num_nodes()
        assert all(data == 1000 for source, target, directed, data in graph.edges())
        topologies.append((graph, completed_vertices))
        return True

    result = generate(progress=report, filter=keep)
    assert result.report.completed and len(result) == len(baseline)
    assert topologies and any(n < g.num_nodes() for g, n in topologies)
    assert updates[0].stage == "topologies"
    assert updates[-1].stage == "complete"
    assert updates[-1].completed == len(result)
    stages = ["topologies", "topology_filters", "interactions", "interaction_filters", "numerators", "selection", "grouping", "complete"]
    assert list(dict.fromkeys(p.stage for p in updates)) == stages
    for before, after in zip(updates, updates[1:]):
        if before.stage == after.stage:
            assert before.completed <= after.completed
    # Snapshots outlive enumeration and never mutate its internal graph.
    topologies[0][0].add_node(42)
    assert len(generate()) == len(baseline)

    updates.clear()
    empty = generate(progress=report, filter=lambda graph, n: False)
    assert empty.report.completed and len(empty) == 0
    assert updates[-1].stage == "complete" and updates[-1].completed == 0

    for keyword in ("progress", "filter"):
        expected = RuntimeError("callback failed")
        def fail(*args):
            raise expected
        try:
            generate(**{keyword: fail})
        except RuntimeError as error:
            assert error is expected
        else:
            raise AssertionError("callback exception was lost")
        try:
            generate(**{keyword: 1})
        except TypeError:
            pass
        else:
            raise AssertionError("non-callable callback was accepted")
    assert len(generate()) == len(baseline)
"#).unwrap();
            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
    }

    #[cfg(unix)]
    #[test]
    fn python_signals_interrupt_both_generation_entry_points() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "symbolica.community.feynkit").unwrap();
            crate::initialize_feynkit(&module).unwrap();
            let locals = PyDict::new(py);
            locals.set_item("fk", &module).unwrap();
            locals
                .set_item(
                    "MODEL_JSON",
                    include_str!("../tests/fixtures/scalars_2p_3p.json"),
                )
                .unwrap();
            let code = CString::new(
                r#"
import signal

model = fk.Model.from_json(MODEL_JSON)
generator = fk.Generator(model)
process = fk.Process.amplitude([1000], [1000, 1000]).with_loop_count(6, 6)
settings = dict(max_vertices=13, threads=2, vertex_allow=["V_3_SCALAR_000"])

def interrupt(signum, frame):
    raise expected_error("generation interrupted")

previous_handler = signal.signal(signal.SIGALRM, interrupt)
try:
    for expected_error in (KeyboardInterrupt, RuntimeError):
        for via_model in (True, False):
            signal.setitimer(signal.ITIMER_REAL, 0.05)
            try:
                if via_model:
                    model.generate_diagrams([1000], [1000, 1000], loops=6, **settings)
                else:
                    generator.generate(process, **settings)
            except expected_error as error:
                assert str(error) == "generation interrupted"
            else:
                raise AssertionError("generation ignored its Python signal handler")
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
finally:
    signal.signal(signal.SIGALRM, previous_handler)

# Interruption belongs to one call; subsequent generation must still complete.
result = model.generate_diagrams([1000], [1000, 1000], loops=0, **settings)
assert result.report.completed
assert len(result) == 1
"#,
            )
            .unwrap();
            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
    }

    #[test]
    fn process_metadata_and_selectors_are_typed_and_round_trip() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "symbolica.community.feynkit").unwrap();
            crate::initialize_feynkit(&module).unwrap();
            let locals = PyDict::new(py);
            locals.set_item("fk", &module).unwrap();
            locals
                .set_item(
                    "MODEL_JSON",
                    include_str!("../tests/fixtures/scalars_2p_3p.json"),
                )
                .unwrap();
            locals
                .set_item(
                    "SM_MODEL_JSON",
                    include_str!("../../feynkit-model/tests/fixtures/sm.json"),
                )
                .unwrap();
            let code = CString::new(
                r#"
by_name = fk.ParticleSelector.by_name("1")
by_pdg = fk.ParticleSelector.by_pdg(1)
model = fk.Model.from_json(MODEL_JSON)
particle = model.particle("scalar_0")

assert by_name.name == "1"
try:
    by_name.pdg
except TypeError:
    pass
else:
    raise AssertionError("named selectors must reject PDG access")
assert by_name.is_name and not by_name.is_pdg
try:
    by_pdg.name
except TypeError:
    pass
else:
    raise AssertionError("PDG selectors must reject name access")
assert by_pdg.pdg == 1
assert by_pdg.is_pdg and not by_pdg.is_name
assert by_name == "1"
assert by_pdg == 1
assert by_name != by_pdg
try:
    by_pdg.pdg = 2
except AttributeError:
    pass
else:
    raise AssertionError("particle selectors must be immutable")

mixed_process = fk.Process.amplitude(
    [particle, by_name, "scalar_2", 1001],
    [particle.antiparticle],
)
assert [
    selector.name if selector.is_name else selector.pdg
    for selector in mixed_process.incoming
] == [1000, "1", "scalar_2", 1001]
assert mixed_process.outgoing_alternatives[0][0].pdg == 1000

sm_model = fk.Model.from_json(SM_MODEL_JSON)
gluon = sm_model.particle_by_pdg(21)
gluon_process = fk.Process.amplitude([gluon], [gluon.antiparticle])
assert gluon_process.incoming[0].pdg == 21
assert gluon_process.outgoing_alternatives[0][0].pdg == 21

bottom = sm_model.particle_by_pdg(5)
bottom_process = fk.Process.amplitude([bottom], [bottom.antiparticle])
assert bottom_process.incoming[0].pdg == 5
assert bottom_process.outgoing_alternatives[0][0].pdg == -5

process = fk.Process.amplitude(
    [by_name, by_pdg],
    [fk.ParticleSelector.by_name("out"), fk.ParticleSelector.by_pdg(-1)],
).with_loop_count(1, 2)
assert process.generation_type == fk.GenerationType.AMPLITUDE
assert process.generation_type == "amplitude"
assert isinstance(process.generation_type, fk.GenerationType)
assert all(isinstance(selector, fk.ParticleSelector) for selector in process.incoming)
assert process.incoming == [by_name, by_pdg]
assert process.outgoing_alternatives[0][0].name == "out"
assert process.outgoing_alternatives[0][1].pdg == -1

round_tripped = fk.Process.amplitude(
    process.incoming,
    process.outgoing_alternatives[0],
).with_loop_count(*process.loop_count)
assert round_tripped.generation_type == process.generation_type
assert round_tripped.incoming == process.incoming
assert round_tripped.outgoing_alternatives == process.outgoing_alternatives
assert round_tripped.loop_count == process.loop_count

cross_section = fk.Process.cross_section([particle, by_pdg], [by_name])
cross_section = cross_section.with_final_state_alternatives(
    [[particle.antiparticle], [by_name], [by_pdg], ["scalar_2"], [1002]],
)
cross_section = cross_section.with_symmetrization(
    initial=True,
    final_state=True,
    left_right=True,
    external_fermions=True,
)
assert cross_section.generation_type == fk.GenerationType.CROSS_SECTION
assert cross_section.generation_type.value == "cross_section"
assert cross_section.incoming[0].pdg == 1000
assert cross_section.outgoing_alternatives[0][0].pdg == 1000
assert cross_section.outgoing_alternatives[1:] == [
    [by_name],
    [by_pdg],
    [fk.ParticleSelector.by_name("scalar_2")],
    [fk.ParticleSelector.by_pdg(1002)],
]
assert cross_section.symmetrizes_initial
assert cross_section.symmetrizes_final
assert cross_section.symmetrizes_left_right
assert cross_section.symmetrizes_external_fermions

for particle_veto in ([particle, 1001], ["scalar_0"]):
    assert len(model.generate_diagrams([particle], [particle, particle], particle_veto=particle_veto)) == 0

foreign_model = fk.Model.from_json(MODEL_JSON.replace("scalar_0", "foreign_scalar_0"))
foreign_particle = foreign_model.particle("foreign_scalar_0")
foreign_process = fk.Process.amplitude([foreign_particle], [foreign_particle.antiparticle])
assert foreign_process.incoming[0].pdg == 1000
assert foreign_process.outgoing_alternatives[0][0].pdg == 1000
"#,
            )
            .unwrap();
            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
    }
}
