use std::{collections::BTreeMap, sync::Arc};

use feynkit_generator::{
    CancellationToken, DiagramGroup, GenerationFilter, GenerationOptions, GenerationReport,
    GenerationResult, GenerationType, Generator, GraphGroupingOptions, GroupMember,
    NumeratorGrouping, ParticleSelector, Process, SelfEnergyFilterOptions, SewnFilterOptions,
    SnailFilterOptions, TadpoleFilterOptions,
};
use feynkit_model::Model;
use pyo3::{
    FromPyObject, IntoPyObjectExt,
    exceptions::{PyIndexError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyAny, PyBool, PyList, PyModule},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    PyStubType, TypeInfo,
    derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods},
};

use crate::{
    display::render_diagram_html,
    error,
    graph::{PyFeynmanDiagram, parse_symbolic_annotation},
    model::PyModel,
};
use symbolica::api::python::PythonExpression;

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
            ParticleSelector::Name(_) => Err(PyTypeError::new_err(
                "a named particle selector does not carry a PDG code",
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

#[derive(FromPyObject)]
pub(crate) enum SelectorInput {
    Selector(PyParticleSelector),
    Name(String),
    Pdg(i64),
}

#[derive(Default)]
pub(crate) struct LoopOrderInput {
    minimum: usize,
    maximum: usize,
}

impl<'py> IntoPyObject<'py> for LoopOrderInput {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        if self.minimum == self.maximum {
            self.minimum.into_bound_py_any(py)
        } else {
            (self.minimum, self.maximum).into_bound_py_any(py)
        }
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for LoopOrderInput {
    type Error = PyErr;

    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if value.is_instance_of::<PyBool>() {
            return Err(PyTypeError::new_err(
                "loops must be a non-negative integer or a (minimum, maximum) pair",
            ));
        }
        let (minimum, maximum) = if let Ok(exact) = value.extract::<usize>() {
            (exact, exact)
        } else if let Ok(range) = value.extract::<(usize, usize)>() {
            range
        } else {
            return Err(PyTypeError::new_err(
                "loops must be a non-negative integer or a (minimum, maximum) pair",
            ));
        };
        if maximum < minimum {
            return Err(PyValueError::new_err(
                "loop bounds must satisfy minimum <= maximum",
            ));
        }
        Ok(Self { minimum, maximum })
    }
}

impl From<SelectorInput> for ParticleSelector {
    fn from(value: SelectorInput) -> Self {
        match value {
            SelectorInput::Selector(selector) => selector.inner,
            SelectorInput::Name(name) => Self::Name(name),
            SelectorInput::Pdg(pdg) => Self::Pdg(pdg),
        }
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for SelectorInput {
    fn type_input() -> TypeInfo {
        PyParticleSelector::type_input() | String::type_input() | i64::type_input()
    }

    fn type_output() -> TypeInfo {
        PyParticleSelector::type_output() | String::type_output() | i64::type_output()
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for LoopOrderInput {
    fn type_input() -> TypeInfo {
        usize::type_input() | <(usize, usize)>::type_input()
    }

    fn type_output() -> TypeInfo {
        usize::type_output() | <(usize, usize)>::type_output()
    }
}

/// A scattering or decay process to pass to the diagram generator.
///
/// A process records its incoming and outgoing particles, loop-order range,
/// and optional external-state symmetrizations.
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
    /// >>> process = fk.Process.amplitude([11, -11], [22])
    ///
    /// Parameters
    /// ----------
    /// incoming : sequence[ParticleSelector | str | int]
    ///     Incoming particles in external-leg order.
    /// outgoing : sequence[ParticleSelector | str | int]
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
    /// incoming : sequence[ParticleSelector | str | int]
    ///     Incoming particles in external-leg order.
    /// outgoing : sequence[ParticleSelector | str | int]
    ///     Outgoing particles in external-leg order.
    #[staticmethod]
    fn cross_section(incoming: Vec<SelectorInput>, outgoing: Vec<SelectorInput>) -> Self {
        Self {
            inner: Process::cross_section(
                incoming.into_iter().map(ParticleSelector::from),
                outgoing.into_iter().map(ParticleSelector::from),
            ),
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
    /// Amplitudes accept exactly one alternative; cross-section alternatives cannot be empty.
    ///
    /// Examples
    /// --------
    /// >>> process = process.with_final_state_alternatives([[22, 22], [23]])
    ///
    /// Parameters
    /// ----------
    /// alternatives : sequence[sequence[ParticleSelector | str | int]]
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
/// Share one token through ``GenerationOptions`` and call ``cancel`` from a
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

/// Configuration for Feynman-diagram generation and filtering.
///
/// Options control parallelism, topology limits, graph filters, numerator
/// grouping, and cancellation without changing the physical process itself.
///
/// Examples
/// --------
/// >>> import symbolica.community.feynkit as fk
/// >>> options = fk.GenerationOptions(threads=4, max_vertices=8)
///
/// Parameters
/// ----------
/// threads : int, optional
///     Number of worker threads used during generation.
/// max_vertices : int, optional
///     Maximum number of interaction vertices in a generated topology.
/// allow_self_loops : bool, optional
///     Permit propagators that start and end on the same vertex.
/// allow_zero_flow_edges : bool, optional
///     Permit internal edges with identically zero momentum flow.
/// graph_prefix : str, optional
///     Prefix assigned to generated diagram names.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "GenerationOptions",
    module = "symbolica.community.feynkit",
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyGenerationOptions {
    pub(crate) inner: GenerationOptions,
}

impl PyGenerationOptions {
    fn add_graph_filter(&mut self, filter: GenerationFilter) {
        self.inner = self.inner.clone().with_graph_filter(filter);
    }

    fn add_cut_amplitude_filter(&mut self, filter: GenerationFilter) {
        self.inner = self.inner.clone().with_cut_amplitude_filter(filter);
    }

    fn grouping_options(
        numerical_sample_seed: u16,
        number_of_numerical_samples: usize,
        differentiate_particle_masses_only: bool,
        fully_numerical_substitution: bool,
        check_canonical_numerator: bool,
    ) -> GraphGroupingOptions {
        GraphGroupingOptions {
            numerical_sample_seed,
            number_of_numerical_samples,
            differentiate_particle_masses_only,
            fully_numerical_substitution,
            check_canonical_numerator,
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyGenerationOptions {
    /// Create generation options with resource limits and topology allowances.
    ///
    /// Examples
    /// --------
    /// >>> options = fk.GenerationOptions(threads=4, max_vertices=8)
    ///
    /// Parameters
    /// ----------
    /// threads : int or None, optional
    ///     Number of worker threads; ``None`` uses the generator default.
    /// max_vertices : int or None, optional
    ///     Maximum number of interaction vertices; ``None`` applies no override.
    /// allow_self_loops : bool, optional
    ///     Allow edges whose two endpoints are the same vertex.
    /// allow_zero_flow_edges : bool, optional
    ///     Allow edges carrying zero momentum flow.
    /// graph_prefix : str or None, optional
    ///     Prefix assigned to generated graph names.
    #[new]
    #[pyo3(signature = (*, threads=None, max_vertices=None, allow_self_loops=false, allow_zero_flow_edges=false, graph_prefix=None))]
    fn new(
        threads: Option<usize>,
        max_vertices: Option<usize>,
        allow_self_loops: bool,
        allow_zero_flow_edges: bool,
        graph_prefix: Option<String>,
    ) -> Self {
        let mut inner = GenerationOptions::default()
            .allow_self_loops(allow_self_loops)
            .allow_zero_flow_edges(allow_zero_flow_edges);
        if let Some(threads) = threads {
            inner = inner.threads(threads);
        }
        if let Some(max_vertices) = max_vertices {
            inner = inner.max_vertices(max_vertices);
        }
        if let Some(prefix) = graph_prefix {
            inner = inner.graph_prefix(prefix);
        }
        Self { inner }
    }

    /// Use a shared token to make generation cancellable.
    ///
    /// Examples
    /// --------
    /// >>> options.set_cancellation_token(token)
    ///
    /// Parameters
    /// ----------
    /// token : CancellationToken
    ///     Token whose cancellation state is checked during generation.
    fn set_cancellation_token(&mut self, token: &PyCancellationToken) {
        self.inner = self.inner.clone().cancellation_token(token.inner.clone());
    }

    /// Reject every graph containing any listed particle PDG code.
    ///
    /// Examples
    /// --------
    /// >>> options.add_particle_veto([6, -6])
    ///
    /// Parameters
    /// ----------
    /// pdg_codes : sequence[int]
    ///     Signed PDG codes forbidden on graph edges.
    fn add_particle_veto(&mut self, pdg_codes: Vec<i64>) {
        self.add_graph_filter(GenerationFilter::ParticleVeto(pdg_codes));
    }

    /// Keep only graphs whose interaction vertices use allowed vertex names.
    ///
    /// Examples
    /// --------
    /// >>> options.add_vertex_allow(["QED_vertex"])
    ///
    /// Parameters
    /// ----------
    /// vertices : sequence[str]
    ///     Model vertex names allowed in generated graphs.
    fn add_vertex_allow(&mut self, vertices: Vec<String>) {
        self.add_graph_filter(GenerationFilter::VertexAllow(vertices));
    }

    /// Reject graphs containing any listed interaction vertex.
    ///
    /// Examples
    /// --------
    /// >>> options.add_vertex_veto(["effective_vertex"])
    ///
    /// Parameters
    /// ----------
    /// vertices : sequence[str]
    ///     Model vertex names forbidden in generated graphs.
    fn add_vertex_veto(&mut self, vertices: Vec<String>) {
        self.add_graph_filter(GenerationFilter::VertexVeto(vertices));
    }

    /// Reject graphs with more than the specified number of bridge edges.
    ///
    /// Examples
    /// --------
    /// >>> options.set_maximum_bridges(2)
    ///
    /// Parameters
    /// ----------
    /// maximum : int
    ///     Largest allowed number of graph bridges.
    fn set_maximum_bridges(&mut self, maximum: usize) {
        self.add_graph_filter(GenerationFilter::MaxNumberOfBridges(maximum));
    }

    /// Configure rejection of self-energy subgraphs by mass category.
    ///
    /// Examples
    /// --------
    /// >>> options.set_self_energy_filter(veto_massive=True, veto_massless=True)
    ///
    /// Parameters
    /// ----------
    /// veto_massive : bool, optional
    ///     Reject self energies carried by massive particles.
    /// veto_massless : bool, optional
    ///     Reject self energies carried by massless particles.
    /// only_scaleless : bool, optional
    ///     Currently unsupported; ``True`` makes generation return an error.
    #[pyo3(signature = (*, veto_massive=true, veto_massless=true, only_scaleless=false))]
    fn set_self_energy_filter(
        &mut self,
        veto_massive: bool,
        veto_massless: bool,
        only_scaleless: bool,
    ) {
        self.add_graph_filter(GenerationFilter::SelfEnergy(SelfEnergyFilterOptions {
            veto_massive,
            veto_massless,
            only_scaleless,
        }));
    }

    /// Configure rejection of tadpoles by the mass of their attachment.
    ///
    /// Examples
    /// --------
    /// >>> options.set_tadpole_filter(veto_attached_to_massless=True)
    ///
    /// Parameters
    /// ----------
    /// veto_attached_to_massive : bool, optional
    ///     Reject tadpoles attached through a massive particle.
    /// veto_attached_to_massless : bool, optional
    ///     Reject tadpoles attached through a massless particle.
    /// only_scaleless : bool, optional
    ///     Currently unsupported; ``True`` makes generation return an error.
    #[pyo3(signature = (*, veto_attached_to_massive=true, veto_attached_to_massless=true, only_scaleless=false))]
    fn set_tadpole_filter(
        &mut self,
        veto_attached_to_massive: bool,
        veto_attached_to_massless: bool,
        only_scaleless: bool,
    ) {
        self.add_graph_filter(GenerationFilter::Tadpoles(TadpoleFilterOptions {
            veto_attached_to_massive,
            veto_attached_to_massless,
            only_scaleless,
        }));
    }

    /// Configure rejection of zero-momentum snail subgraphs.
    ///
    /// Examples
    /// --------
    /// >>> options.set_zero_snail_filter(veto_attached_to_massless=True)
    ///
    /// Parameters
    /// ----------
    /// veto_attached_to_massive : bool, optional
    ///     Reject zero-momentum snails attached through a massive particle.
    /// veto_attached_to_massless : bool, optional
    ///     Reject zero-momentum snails attached through a massless particle.
    /// only_scaleless : bool, optional
    ///     Currently unsupported; ``True`` makes generation return an error.
    #[pyo3(signature = (*, veto_attached_to_massive=false, veto_attached_to_massless=true, only_scaleless=false))]
    fn set_zero_snail_filter(
        &mut self,
        veto_attached_to_massive: bool,
        veto_attached_to_massless: bool,
        only_scaleless: bool,
    ) {
        self.add_graph_filter(GenerationFilter::ZeroSnails(SnailFilterOptions {
            veto_attached_to_massive,
            veto_attached_to_massless,
            only_scaleless,
        }));
    }

    /// Restrict total coupling-order powers to inclusive ranges.
    ///
    /// Examples
    /// --------
    /// >>> options.set_coupling_orders({"QED": (2, 4), "QCD": (0, None)})
    ///
    /// Parameters
    /// ----------
    /// orders : dict[str, tuple[int, int or None]]
    ///     Coupling name mapped to its minimum and optional inclusive maximum power.
    fn set_coupling_orders(&mut self, orders: BTreeMap<String, (usize, Option<usize>)>) {
        self.add_graph_filter(GenerationFilter::CouplingOrders(orders));
    }

    /// Restrict generated graphs to an inclusive loop-count range.
    ///
    /// Examples
    /// --------
    /// >>> options.set_loop_count_range(1, 2)
    ///
    /// Parameters
    /// ----------
    /// minimum : int
    ///     Minimum number of loops.
    /// maximum : int
    ///     Maximum number of loops, inclusive.
    fn set_loop_count_range(&mut self, minimum: usize, maximum: usize) {
        self.add_graph_filter(GenerationFilter::LoopCountRange((minimum, maximum)));
    }

    /// Restrict generated graphs to an inclusive fermion-loop-count range.
    ///
    /// Examples
    /// --------
    /// >>> options.set_fermion_loop_count_range(0, 1)
    ///
    /// Parameters
    /// ----------
    /// minimum : int
    ///     Minimum number of closed fermion loops.
    /// maximum : int
    ///     Maximum number of closed fermion loops, inclusive.
    fn set_fermion_loop_count_range(&mut self, minimum: usize, maximum: usize) {
        self.add_graph_filter(GenerationFilter::FermionLoopCountRange((minimum, maximum)));
    }

    /// Restrict the number of factorized loop-topology components.
    ///
    /// Examples
    /// --------
    /// >>> options.set_factorized_loop_topologies_count_range(1, 2)
    ///
    /// Parameters
    /// ----------
    /// minimum : int
    ///     Minimum number of factorized loop-topology components.
    /// maximum : int
    ///     Maximum number of components, inclusive.
    fn set_factorized_loop_topologies_count_range(&mut self, minimum: usize, maximum: usize) {
        self.add_graph_filter(GenerationFilter::FactorizedLoopTopologiesCountRange((
            minimum, maximum,
        )));
    }

    /// Restrict cross-section graphs to an inclusive blob-count range.
    ///
    /// Examples
    /// --------
    /// >>> options.set_blob_range(1, 2)
    ///
    /// Parameters
    /// ----------
    /// minimum : int
    ///     Minimum number of blobs.
    /// maximum : int
    ///     Maximum number of blobs, inclusive.
    fn set_blob_range(&mut self, minimum: usize, maximum: usize) {
        self.add_graph_filter(GenerationFilter::BlobRange(minimum..=maximum));
    }

    /// Restrict cross-section graphs to an inclusive spectator-count range.
    ///
    /// Examples
    /// --------
    /// >>> options.set_spectator_range(0, 1)
    ///
    /// Parameters
    /// ----------
    /// minimum : int
    ///     Minimum number of spectator lines.
    /// maximum : int
    ///     Maximum number of spectator lines, inclusive.
    fn set_spectator_range(&mut self, minimum: usize, maximum: usize) {
        self.add_graph_filter(GenerationFilter::SpectatorRange(minimum..=maximum));
    }

    /// Require exact perturbative orders for cross-section graphs.
    ///
    /// Examples
    /// --------
    /// >>> options.set_perturbative_orders({"QED": 2, "QCD": 1})
    ///
    /// Parameters
    /// ----------
    /// orders : dict[str, int]
    ///     Coupling name mapped to its required perturbative power.
    fn set_perturbative_orders(&mut self, orders: BTreeMap<String, usize>) {
        self.add_graph_filter(GenerationFilter::PerturbativeOrders(orders));
    }

    /// Configure rejection of tadpole topologies revealed by sewing cross-section sides.
    ///
    /// Examples
    /// --------
    /// >>> options.set_sewn_filter(filter_tadpoles=True)
    ///
    /// Parameters
    /// ----------
    /// filter_tadpoles : bool, optional
    ///     Reject sewn tadpole topologies; ``False`` disables this check.
    #[pyo3(signature = (*, filter_tadpoles=true))]
    fn set_sewn_filter(&mut self, filter_tadpoles: bool) {
        self.add_graph_filter(GenerationFilter::Sewn(SewnFilterOptions {
            filter_tadpoles,
        }));
    }

    /// Restrict coupling orders independently within every cut amplitude.
    ///
    /// Examples
    /// --------
    /// >>> options.set_cut_amplitude_coupling_orders({"QED": (1, 2)})
    ///
    /// Parameters
    /// ----------
    /// orders : dict[str, tuple[int, int or None]]
    ///     Coupling name mapped to its minimum and optional inclusive maximum power.
    fn set_cut_amplitude_coupling_orders(
        &mut self,
        orders: BTreeMap<String, (usize, Option<usize>)>,
    ) {
        self.add_cut_amplitude_filter(GenerationFilter::CouplingOrders(orders));
    }

    /// Restrict the summed loop count across both sides of every cut.
    ///
    /// Examples
    /// --------
    /// >>> options.set_cut_amplitude_loop_count_range(0, 1)
    ///
    /// Parameters
    /// ----------
    /// minimum : int
    ///     Minimum combined loop count of the two cut amplitudes.
    /// maximum : int
    ///     Maximum combined loop count, inclusive.
    fn set_cut_amplitude_loop_count_range(&mut self, minimum: usize, maximum: usize) {
        self.add_cut_amplitude_filter(GenerationFilter::LoopCountRange((minimum, maximum)));
    }

    /// Disable numerator parsing, zero detection, and cross-diagram grouping.
    ///
    /// Examples
    /// --------
    /// >>> options.disable_numerator_grouping()
    ///
    fn disable_numerator_grouping(&mut self) {
        self.inner = self
            .inner
            .clone()
            .numerator_grouping(NumeratorGrouping::None);
    }

    /// Detect zero numerators without grouping the remaining diagrams.
    ///
    /// Examples
    /// --------
    /// >>> options.detect_zero_numerators()
    ///
    fn detect_zero_numerators(&mut self) {
        self.inner = self
            .inner
            .clone()
            .numerator_grouping(NumeratorGrouping::OnlyDetectZeroes);
    }

    /// Group diagrams only when their numerators are identical.
    ///
    /// Examples
    /// --------
    /// >>> options.group_identical_numerators(number_of_numerical_samples=7)
    ///
    /// Parameters
    /// ----------
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
    #[pyo3(signature = (*, numerical_sample_seed=3, number_of_numerical_samples=5, differentiate_particle_masses_only=true, fully_numerical_substitution=false, check_canonical_numerator=false))]
    fn group_identical_numerators(
        &mut self,
        numerical_sample_seed: u16,
        number_of_numerical_samples: usize,
        differentiate_particle_masses_only: bool,
        fully_numerical_substitution: bool,
        check_canonical_numerator: bool,
    ) {
        let options = Self::grouping_options(
            numerical_sample_seed,
            number_of_numerical_samples,
            differentiate_particle_masses_only,
            fully_numerical_substitution,
            check_canonical_numerator,
        );
        self.inner = self
            .inner
            .clone()
            .numerator_grouping(NumeratorGrouping::Identical(options));
    }

    /// Group diagrams whose numerators differ only by an overall sign.
    ///
    /// Examples
    /// --------
    /// >>> options.group_numerators_up_to_sign(check_canonical_numerator=True)
    ///
    /// Parameters
    /// ----------
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
    #[pyo3(signature = (*, numerical_sample_seed=3, number_of_numerical_samples=5, differentiate_particle_masses_only=true, fully_numerical_substitution=false, check_canonical_numerator=false))]
    fn group_numerators_up_to_sign(
        &mut self,
        numerical_sample_seed: u16,
        number_of_numerical_samples: usize,
        differentiate_particle_masses_only: bool,
        fully_numerical_substitution: bool,
        check_canonical_numerator: bool,
    ) {
        let options = Self::grouping_options(
            numerical_sample_seed,
            number_of_numerical_samples,
            differentiate_particle_masses_only,
            fully_numerical_substitution,
            check_canonical_numerator,
        );
        self.inner = self
            .inner
            .clone()
            .numerator_grouping(NumeratorGrouping::UpToSign(options));
    }

    /// Group diagrams whose numerators differ by a scalar factor.
    ///
    /// Examples
    /// --------
    /// >>> options.group_numerators_up_to_scalar(fully_numerical_substitution=True)
    ///
    /// Parameters
    /// ----------
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
    #[pyo3(signature = (*, numerical_sample_seed=3, number_of_numerical_samples=5, differentiate_particle_masses_only=true, fully_numerical_substitution=false, check_canonical_numerator=false))]
    fn group_numerators_up_to_scalar(
        &mut self,
        numerical_sample_seed: u16,
        number_of_numerical_samples: usize,
        differentiate_particle_masses_only: bool,
        fully_numerical_substitution: bool,
        check_canonical_numerator: bool,
    ) {
        let options = Self::grouping_options(
            numerical_sample_seed,
            number_of_numerical_samples,
            differentiate_particle_masses_only,
            fully_numerical_substitution,
            check_canonical_numerator,
        );
        self.inner = self
            .inner
            .clone()
            .numerator_grouping(NumeratorGrouping::UpToScalar(options));
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
/// >>> member = result.groups[0].members[0]
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
    /// Return the member's index in the retained ``GenerationResult.diagrams``.
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
    fn ratio(&self) -> &str {
        &self.inner.ratio
    }
    /// Parse the numerator ratio as a native Symbolica expression.
    ///
    /// Examples
    /// --------
    /// >>> ratio = result.groups[0].members[0].ratio_expression()
    ///
    fn ratio_expression(&self) -> PyResult<PythonExpression> {
        parse_symbolic_annotation(&self.inner.ratio).map_err(error::GenerationError::new_err)
    }
}

/// Diagrams whose numerators are related by known scalar factors.
///
/// Grouping lets amplitude calculations evaluate one master numerator and
/// reconstruct related diagrams from each member's symbolic ratio.
///
/// Examples
/// --------
/// >>> group = result.groups[0]
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
    /// >>> result = generator.generate(process, options)
    /// >>> combined_numerator = result.diagrams[0].numerator_expression()
    ///
    /// Parameters
    /// ----------
    /// process : Process
    ///     Scattering process and loop range to generate.
    /// options : GenerationOptions or None, optional
    ///     Filters, limits, grouping mode, and cancellation state to apply.
    #[pyo3(signature = (process, options=None))]
    fn generate(
        &self,
        py: Python<'_>,
        process: &PyProcess,
        options: Option<&PyGenerationOptions>,
    ) -> PyResult<PyGenerationResult> {
        let generator = self.inner.clone();
        let process = process.inner.clone();
        let options = options.cloned().unwrap_or_default().inner;
        py.detach(move || generator.generate(&process, &options))
            .map(|inner| PyGenerationResult { inner })
            .map_err(error::generation)
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn generate_diagrams_for_model(
    py: Python<'_>,
    model: &PyModel,
    incoming: Vec<SelectorInput>,
    outgoing: Vec<SelectorInput>,
    kind: &str,
    loops: LoopOrderInput,
    options: Option<&PyGenerationOptions>,
    final_state_alternatives: Option<Vec<Vec<SelectorInput>>>,
) -> PyResult<PyGenerationResult> {
    let LoopOrderInput { minimum, maximum } = loops;
    let incoming = incoming
        .into_iter()
        .map(ParticleSelector::from)
        .collect::<Vec<_>>();
    let outgoing = outgoing
        .into_iter()
        .map(ParticleSelector::from)
        .collect::<Vec<_>>();
    let alternatives = final_state_alternatives
        .unwrap_or_default()
        .into_iter()
        .map(|alternative| {
            alternative
                .into_iter()
                .map(ParticleSelector::from)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let process = match kind {
        "amplitude" => {
            if !alternatives.is_empty() {
                return Err(PyValueError::new_err(
                    "final_state_alternatives is only supported for cross sections",
                ));
            }
            Process::amplitude(incoming, outgoing)
                .with_loop_count(minimum, maximum)
                .map_err(error::process)?
        }
        "cross_section" => {
            let mut process = Process::cross_section(incoming, outgoing.clone())
                .with_loop_count(minimum, maximum)
                .map_err(error::process)?;
            if !alternatives.is_empty() {
                let mut final_states = Vec::with_capacity(alternatives.len() + 1);
                final_states.push(outgoing);
                final_states.extend(alternatives);
                process = process
                    .with_final_state_alternatives(final_states)
                    .map_err(error::process)?;
            }
            process
        }
        _ => {
            return Err(PyValueError::new_err(
                "kind must be 'amplitude' or 'cross_section'",
            ));
        }
    };

    let generator = Generator::new(model.inner.clone());
    let options = options.cloned().unwrap_or_default().inner;
    py.detach(move || generator.generate(&process, &options))
        .map(|inner| PyGenerationResult { inner })
        .map_err(error::generation)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyGenerationType>()?;
    module.add_class::<PyParticleSelector>()?;
    module.add_class::<PyProcess>()?;
    module.add_class::<PyCancellationToken>()?;
    module.add_class::<PyGenerationOptions>()?;
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
    fn process_metadata_and_selectors_are_typed_and_round_trip() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "symbolica.community.feynkit").unwrap();
            crate::initialize_feynkit(&module).unwrap();
            let locals = PyDict::new(py);
            locals.set_item("fk", &module).unwrap();
            let code = CString::new(
                r#"
by_name = fk.ParticleSelector.by_name("1")
by_pdg = fk.ParticleSelector.by_pdg(1)

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

cross_section = fk.Process.cross_section([by_pdg], [by_name])
cross_section = cross_section.with_final_state_alternatives([[by_name], [by_pdg]])
cross_section = cross_section.with_symmetrization(
    initial=True,
    final_state=True,
    left_right=True,
    external_fermions=True,
)
assert cross_section.generation_type == fk.GenerationType.CROSS_SECTION
assert cross_section.generation_type.value == "cross_section"
assert cross_section.outgoing_alternatives == [[by_name], [by_pdg]]
assert cross_section.symmetrizes_initial
assert cross_section.symmetrizes_final
assert cross_section.symmetrizes_left_right
assert cross_section.symmetrizes_external_fermions
"#,
            )
            .unwrap();
            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
    }
}
