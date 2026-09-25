use std::{
    collections::{BTreeMap, VecDeque},
    sync::{Arc, Mutex},
    time::{Duration, Instant},
};

use feynkit_generator::{
    CancellationToken, DiagramGroup, EdgeColor, FilterScope, GenerationControl, GenerationFilter,
    GenerationOptions, GenerationProgress, GenerationReport, GenerationResult, GenerationType,
    GraphGroupingOptions, GroupMember, NodeColor, NumeratorGrouping, ParticleSelector, Process,
    SelfEnergyFilterOptions, SewnFilterOptions, SnailFilterOptions, TadpoleFilterOptions,
    VertexSelector,
};
use feynkit_graph::{DiagramId, EdgeId};
use feynkit_model::Model;
use pyo3::{
    FromPyObject, IntoPyObjectExt,
    exceptions::{PyIndexError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyAny, PyBool, PyDict, PyEllipsis, PyList, PyModule, PyString},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    PyStubType, TypeInfo,
    derive::{gen_stub_pyclass, gen_stub_pymethods},
};

use crate::{
    amplitude::PyAmplitude,
    error,
    graph::PyFeynmanDiagram,
    model::{PyModel, PyParticle, PyVertexRule},
};
use symbolica::{
    api::python::{ConvertibleToExpression, PythonExpression, PythonGraph},
    graph::Graph,
};

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
/// >>> process = fk.Process(model, ["e-", "e+"], ["mu-", "mu+"])
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
    model: Arc<Model>,
    final_state_symmetry: Option<bool>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyProcess {
    /// Define a model-bound process independently of what will be generated.
    ///
    /// Examples
    /// --------
    /// >>> process = fk.Process(model, ["e-", "e+"], ["a", "a"])
    /// >>> amplitude = process.generate_amplitude()
    ///
    /// Parameters
    /// ----------
    /// model : Model
    ///     Particle model supplying the Feynman rules.
    /// incoming : sequence[Particle | ParticleSelector | str | int]
    ///     Ordered incoming external states.
    /// outgoing : sequence[Particle | ParticleSelector | str | int]
    ///     Ordered outgoing external states.
    /// loops : int or tuple[int, int], optional
    ///     Default exact loop order or inclusive range, initially zero.
    #[new]
    #[pyo3(signature = (model, incoming, outgoing, *, loops=OrderRangeInput::default()))]
    fn new(
        model: &PyModel,
        incoming: Vec<SelectorInput>,
        outgoing: Vec<SelectorInput>,
        loops: OrderRangeInput,
    ) -> PyResult<Self> {
        let inner = Process::new(
            incoming.into_iter().map(ParticleSelector::from),
            outgoing.into_iter().map(ParticleSelector::from),
        )
        .with_loop_count(loops.minimum, loops.maximum.unwrap())
        .map_err(error::process)?;
        Ok(Self {
            inner,
            model: model.inner.clone(),
            final_state_symmetry: None,
        })
    }

    /// The particle model supplying this process's Feynman rules.
    #[getter]
    fn model(&self) -> PyModel {
        self.model.clone().into()
    }

    /// Format the external states, model, and default loop range.
    ///
    /// Examples
    /// --------
    /// >>> repr(process)
    /// 'Process("sm": [e-, e+] -> [a, a]; loops=(0, 0))'
    fn __repr__(&self) -> String {
        let names = |state: &[ParticleSelector]| {
            state
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        };
        format!(
            "Process({:?}: [{}] -> {}; loops={:?})",
            self.model.name(),
            names(self.inner.incoming()),
            self.inner
                .outgoing_alternatives()
                .iter()
                .map(|state| format!("[{}]", names(state)))
                .collect::<Vec<_>>()
                .join(" | "),
            self.loop_count()
        )
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
            .map(|inner| Self {
                inner,
                ..self.clone()
            })
            .map_err(error::process)
    }

    /// Return a process accepting any of the supplied final states.
    /// Generating amplitudes requires exactly one alternative; cross sections may include an empty
    /// alternative for a vacuum final state.
    ///
    /// Examples
    /// --------
    /// >>> process = fk.Process(model, [11, -11], [22, 22])
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
            .map(|inner| Self {
                inner,
                ..self.clone()
            })
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
            final_state_symmetry: Some(final_state),
            ..self.clone()
        }
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
    /// >>> result = process.generate_diagrams(max_vertices=6)
    /// >>> combined_numerator = result.diagrams[0].numerator_expression()
    ///
    /// Parameters
    /// ----------
    /// loops : int or tuple[int, int] or None, optional
    ///     Override the process loop range for this call. Cross sections count
    ///     loops in the sewn forward graph; two-particle tree cuts need loops=1.
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
    /// maximum_bridges : int, None, or Ellipsis, optional
    ///     Omission or Ellipsis requires one-particle irreducibility only for diagrams
    ///     with loops; tree exchanges are allowed, including in mixed loop ranges.
    ///     An integer limits internal bridges at every loop order; None disables it.
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
    ///     Defaults to None: no numerator comparison or grouping. Diagrams still
    ///     contain numerators. Pass NumeratorGrouping to enable zero detection or grouping.
    /// progress : {"auto"}, Callable[[GenerationProgress], None] or None, optional
    ///     Defaults to "auto": show progress when marimo.running_in_notebook()
    ///     is true, with stage, counts, and elapsed time. None disables progress.
    ///     Known totals use a progress bar; unknown totals use a spinner.
    ///     The display closes on completion, cancellation, or error.
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
    #[pyo3(signature = (*, loops=None, threads=None, max_vertices=None, allow_self_loops=true, allow_zero_flow_edges=false, graph_prefix=None, particle_veto=None, vertex_allow=None, vertex_veto=None, maximum_bridges=Some(Python::attach(|py| py.Ellipsis())), self_energy=Some(Python::attach(|py| py.Ellipsis())), tadpoles=Some(Python::attach(|py| py.Ellipsis())), zero_snails=Some(Python::attach(|py| py.Ellipsis())), coupling_orders=None, fermion_loop_count_range=None, factorized_loop_topologies_count_range=Some(Python::attach(|py| py.Ellipsis())), blob_range=Some(Python::attach(|py| py.Ellipsis())), spectator_range=Some(Python::attach(|py| py.Ellipsis())), perturbative_orders=None, sewn_tadpoles=None, cut_amplitude_coupling_orders=None, cut_amplitude_loop_count_range=None, select_diagrams=None, veto_diagrams=None, loop_momentum_bases=None, numerator_prefactor=None, projector=None, numerator_grouping=None, cancellation_token=None, progress=Some(Python::attach(|py| PyString::new(py, "auto").into_any().unbind())), filter=None))]
    #[pyo3(
        text_signature = "($self, *, loops=None, threads=None, max_vertices=None, allow_self_loops=True, allow_zero_flow_edges=False, graph_prefix=None, particle_veto=None, vertex_allow=None, vertex_veto=None, maximum_bridges=..., self_energy=..., tadpoles=..., zero_snails=..., coupling_orders=None, fermion_loop_count_range=None, factorized_loop_topologies_count_range=..., blob_range=..., spectator_range=..., perturbative_orders=None, sewn_tadpoles=None, cut_amplitude_coupling_orders=None, cut_amplitude_loop_count_range=None, select_diagrams=None, veto_diagrams=None, loop_momentum_bases=None, numerator_prefactor=None, projector=None, numerator_grouping=None, cancellation_token=None, progress='auto', filter=None)"
    )]
    #[allow(clippy::too_many_arguments)]
    fn generate_diagrams(
        &self,
        py: Python<'_>,
        loops: Option<OrderRangeInput>,
        threads: Option<usize>,
        max_vertices: Option<usize>,
        allow_self_loops: bool,
        allow_zero_flow_edges: bool,
        graph_prefix: Option<String>,
        particle_veto: Option<Vec<ParticleInput>>,
        vertex_allow: Option<Vec<VertexInput>>,
        vertex_veto: Option<Vec<VertexInput>>,
        #[gen_stub(override_type(type_repr = "int | None | types.EllipsisType"))]
        maximum_bridges: Option<Py<PyAny>>,
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
        numerator_grouping: Option<PyNumeratorGrouping>,
        cancellation_token: Option<PyCancellationToken>,
        #[gen_stub(override_type(type_repr = "typing.Literal['auto'] | collections.abc.Callable[[GenerationProgress], None] | None", imports = ("collections.abc", "typing")))]
        progress: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "collections.abc.Callable[[symbolica.core.Graph, int], bool] | None", imports = ("collections.abc", "symbolica.core")))]
        filter: Option<Py<PyAny>>,
    ) -> PyResult<PyGenerationResult> {
        let generation_type = GenerationType::Amplitude;
        let mut process = self.inner.clone().symmetrize_final(
            self.final_state_symmetry
                .unwrap_or(generation_type == GenerationType::CrossSection),
        );
        if let Some(loops) = loops {
            process = process
                .with_loop_count(loops.minimum, loops.maximum.unwrap())
                .map_err(error::process)?;
        }
        let options = GenerationSettings::new(
            py,
            &process,
            generation_type,
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
        run_generation(
            py,
            self.model.clone(),
            process,
            generation_type,
            options,
            progress,
            filter,
        )
    }
    /// Generate a coherent symbolic amplitude for this process.
    /// Cancelled or empty generation cannot produce an amplitude.
    ///
    /// Model Feynman rules are instantiated into each returned diagram: inspect
    /// ``vertex.numerator_expression()`` and ``edge.numerator_expression()`` for
    /// individual factors, or ``diagram.numerator_expression()`` for their
    /// combined numerator.
    ///
    /// Examples
    /// --------
    /// >>> result = process.generate_amplitude(max_vertices=6)
    /// >>> operator = result.expression()
    ///
    /// Parameters
    /// ----------
    /// dimension : int or Expression, optional
    ///     Lorentz dimension of the amplitude; default four.
    /// real : list[Expression] or None, optional
    ///     Additional scalars assumed real under conjugation.
    /// loops : int or tuple[int, int] or None, optional
    ///     Override the process loop range for this call. Cross sections count
    ///     loops in the sewn forward graph; two-particle tree cuts need loops=1.
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
    /// maximum_bridges : int, None, or Ellipsis, optional
    ///     Omission or Ellipsis requires one-particle irreducibility only for diagrams
    ///     with loops; tree exchanges are allowed, including in mixed loop ranges.
    ///     An integer limits internal bridges at every loop order; None disables it.
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
    ///     Defaults to None: no numerator comparison or grouping. Diagrams still
    ///     contain numerators. Pass NumeratorGrouping to enable zero detection or grouping.
    /// progress : {"auto"}, Callable[[GenerationProgress], None] or None, optional
    ///     Defaults to "auto": show progress when marimo.running_in_notebook()
    ///     is true, with stage, counts, and elapsed time. None disables progress.
    ///     Known totals use a progress bar; unknown totals use a spinner.
    ///     The display closes on completion, cancellation, or error.
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
    #[pyo3(signature = (*, dimension=None, real=None, loops=None, threads=None, max_vertices=None, allow_self_loops=true, allow_zero_flow_edges=false, graph_prefix=None, particle_veto=None, vertex_allow=None, vertex_veto=None, maximum_bridges=Some(Python::attach(|py| py.Ellipsis())), self_energy=Some(Python::attach(|py| py.Ellipsis())), tadpoles=Some(Python::attach(|py| py.Ellipsis())), zero_snails=Some(Python::attach(|py| py.Ellipsis())), coupling_orders=None, fermion_loop_count_range=None, factorized_loop_topologies_count_range=Some(Python::attach(|py| py.Ellipsis())), blob_range=Some(Python::attach(|py| py.Ellipsis())), spectator_range=Some(Python::attach(|py| py.Ellipsis())), perturbative_orders=None, sewn_tadpoles=None, cut_amplitude_coupling_orders=None, cut_amplitude_loop_count_range=None, select_diagrams=None, veto_diagrams=None, loop_momentum_bases=None, numerator_prefactor=None, projector=None, numerator_grouping=None, cancellation_token=None, progress=Some(Python::attach(|py| PyString::new(py, "auto").into_any().unbind())), filter=None))]
    #[pyo3(
        text_signature = "($self, *, dimension=None, real=None, loops=None, threads=None, max_vertices=None, allow_self_loops=True, allow_zero_flow_edges=False, graph_prefix=None, particle_veto=None, vertex_allow=None, vertex_veto=None, maximum_bridges=..., self_energy=..., tadpoles=..., zero_snails=..., coupling_orders=None, fermion_loop_count_range=None, factorized_loop_topologies_count_range=..., blob_range=..., spectator_range=..., perturbative_orders=None, sewn_tadpoles=None, cut_amplitude_coupling_orders=None, cut_amplitude_loop_count_range=None, select_diagrams=None, veto_diagrams=None, loop_momentum_bases=None, numerator_prefactor=None, projector=None, numerator_grouping=None, cancellation_token=None, progress='auto', filter=None)"
    )]
    #[allow(clippy::too_many_arguments)]
    fn generate_amplitude(
        &self,
        py: Python<'_>,
        dimension: Option<ConvertibleToExpression>,
        real: Option<Vec<PythonExpression>>,
        loops: Option<OrderRangeInput>,
        threads: Option<usize>,
        max_vertices: Option<usize>,
        allow_self_loops: bool,
        allow_zero_flow_edges: bool,
        graph_prefix: Option<String>,
        particle_veto: Option<Vec<ParticleInput>>,
        vertex_allow: Option<Vec<VertexInput>>,
        vertex_veto: Option<Vec<VertexInput>>,
        #[gen_stub(override_type(type_repr = "int | None | types.EllipsisType"))]
        maximum_bridges: Option<Py<PyAny>>,
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
        numerator_grouping: Option<PyNumeratorGrouping>,
        cancellation_token: Option<PyCancellationToken>,
        #[gen_stub(override_type(type_repr = "typing.Literal['auto'] | collections.abc.Callable[[GenerationProgress], None] | None", imports = ("collections.abc", "typing")))]
        progress: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "collections.abc.Callable[[symbolica.core.Graph, int], bool] | None", imports = ("collections.abc", "symbolica.core")))]
        filter: Option<Py<PyAny>>,
    ) -> PyResult<PyAmplitude> {
        let generation_type = GenerationType::Amplitude;
        let mut process = self.inner.clone().symmetrize_final(
            self.final_state_symmetry
                .unwrap_or(generation_type == GenerationType::CrossSection),
        );
        if let Some(loops) = loops {
            process = process
                .with_loop_count(loops.minimum, loops.maximum.unwrap())
                .map_err(error::process)?;
        }
        let options = GenerationSettings::new(
            py,
            &process,
            generation_type,
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
        let result = run_generation(
            py,
            self.model.clone(),
            process,
            generation_type,
            options,
            progress,
            filter,
        )?;
        if !result.inner.report.completed {
            return Err(error::generation(
                feynkit_generator::GenerationError::IncompleteAmplitude,
            ));
        }
        PyAmplitude::new(result.diagrams(), dimension, real)
    }
    /// Generate sewn forward diagrams and their physical final-state cuts.
    /// The result contains diagrams and cut metadata, before phase-space integration.
    ///
    /// Model Feynman rules are instantiated into each returned diagram: inspect
    /// ``vertex.numerator_expression()`` and ``edge.numerator_expression()`` for
    /// individual factors, or ``diagram.numerator_expression()`` for their
    /// combined numerator.
    ///
    /// Examples
    /// --------
    /// >>> result = process.generate_cross_section(max_vertices=6)
    /// >>> combined_numerator = result.diagrams[0].numerator_expression()
    ///
    /// Parameters
    /// ----------
    /// loops : int or tuple[int, int] or None, optional
    ///     Override the process loop range for this call. Cross sections count
    ///     loops in the sewn forward graph; two-particle tree cuts need loops=1.
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
    /// maximum_bridges : int, None, or Ellipsis, optional
    ///     Omission or Ellipsis requires one-particle irreducibility only for diagrams
    ///     with loops; tree exchanges are allowed, including in mixed loop ranges.
    ///     An integer limits internal bridges at every loop order; None disables it.
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
    ///     Defaults to None: no numerator comparison or grouping. Diagrams still
    ///     contain numerators. Pass NumeratorGrouping to enable zero detection or grouping.
    /// progress : {"auto"}, Callable[[GenerationProgress], None] or None, optional
    ///     Defaults to "auto": show progress when marimo.running_in_notebook()
    ///     is true, with stage, counts, and elapsed time. None disables progress.
    ///     Known totals use a progress bar; unknown totals use a spinner.
    ///     The display closes on completion, cancellation, or error.
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
    #[pyo3(signature = (*, loops=None, threads=None, max_vertices=None, allow_self_loops=true, allow_zero_flow_edges=false, graph_prefix=None, particle_veto=None, vertex_allow=None, vertex_veto=None, maximum_bridges=Some(Python::attach(|py| py.Ellipsis())), self_energy=Some(Python::attach(|py| py.Ellipsis())), tadpoles=Some(Python::attach(|py| py.Ellipsis())), zero_snails=Some(Python::attach(|py| py.Ellipsis())), coupling_orders=None, fermion_loop_count_range=None, factorized_loop_topologies_count_range=Some(Python::attach(|py| py.Ellipsis())), blob_range=Some(Python::attach(|py| py.Ellipsis())), spectator_range=Some(Python::attach(|py| py.Ellipsis())), perturbative_orders=None, sewn_tadpoles=None, cut_amplitude_coupling_orders=None, cut_amplitude_loop_count_range=None, select_diagrams=None, veto_diagrams=None, loop_momentum_bases=None, numerator_prefactor=None, projector=None, numerator_grouping=None, cancellation_token=None, progress=Some(Python::attach(|py| PyString::new(py, "auto").into_any().unbind())), filter=None))]
    #[pyo3(
        text_signature = "($self, *, loops=None, threads=None, max_vertices=None, allow_self_loops=True, allow_zero_flow_edges=False, graph_prefix=None, particle_veto=None, vertex_allow=None, vertex_veto=None, maximum_bridges=..., self_energy=..., tadpoles=..., zero_snails=..., coupling_orders=None, fermion_loop_count_range=None, factorized_loop_topologies_count_range=..., blob_range=..., spectator_range=..., perturbative_orders=None, sewn_tadpoles=None, cut_amplitude_coupling_orders=None, cut_amplitude_loop_count_range=None, select_diagrams=None, veto_diagrams=None, loop_momentum_bases=None, numerator_prefactor=None, projector=None, numerator_grouping=None, cancellation_token=None, progress='auto', filter=None)"
    )]
    #[allow(clippy::too_many_arguments)]
    fn generate_cross_section(
        &self,
        py: Python<'_>,
        loops: Option<OrderRangeInput>,
        threads: Option<usize>,
        max_vertices: Option<usize>,
        allow_self_loops: bool,
        allow_zero_flow_edges: bool,
        graph_prefix: Option<String>,
        particle_veto: Option<Vec<ParticleInput>>,
        vertex_allow: Option<Vec<VertexInput>>,
        vertex_veto: Option<Vec<VertexInput>>,
        #[gen_stub(override_type(type_repr = "int | None | types.EllipsisType"))]
        maximum_bridges: Option<Py<PyAny>>,
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
        numerator_grouping: Option<PyNumeratorGrouping>,
        cancellation_token: Option<PyCancellationToken>,
        #[gen_stub(override_type(type_repr = "typing.Literal['auto'] | collections.abc.Callable[[GenerationProgress], None] | None", imports = ("collections.abc", "typing")))]
        progress: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "collections.abc.Callable[[symbolica.core.Graph, int], bool] | None", imports = ("collections.abc", "symbolica.core")))]
        filter: Option<Py<PyAny>>,
    ) -> PyResult<PyGenerationResult> {
        let generation_type = GenerationType::CrossSection;
        let mut process = self.inner.clone().symmetrize_final(
            self.final_state_symmetry
                .unwrap_or(generation_type == GenerationType::CrossSection),
        );
        if let Some(loops) = loops {
            process = process
                .with_loop_count(loops.minimum, loops.maximum.unwrap())
                .map_err(error::process)?;
        }
        let options = GenerationSettings::new(
            py,
            &process,
            generation_type,
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
        run_generation(
            py,
            self.model.clone(),
            process,
            generation_type,
            options,
            progress,
            filter,
        )
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

    /// Final-state symmetry override; None uses False for amplitudes and True for cross sections.
    ///
    /// Examples
    /// --------
    /// >>> process.with_symmetrization(final_state=True).symmetrizes_final
    /// True
    ///
    #[getter]
    fn symmetrizes_final(&self) -> Option<bool> {
        self.final_state_symmetry
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
struct GenerationSettings {
    inner: GenerationOptions,
}

impl GenerationSettings {
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        process: &Process,
        generation_type: GenerationType,
        threads: Option<usize>,
        max_vertices: Option<usize>,
        allow_self_loops: bool,
        allow_zero_flow_edges: bool,
        graph_prefix: Option<String>,
        particle_veto: Option<Vec<ParticleInput>>,
        vertex_allow: Option<Vec<VertexInput>>,
        vertex_veto: Option<Vec<VertexInput>>,
        maximum_bridges: Option<Py<PyAny>>,
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
        numerator_grouping: Option<PyNumeratorGrouping>,
        cancellation_token: Option<PyCancellationToken>,
    ) -> PyResult<Self> {
        // Ported from GammaLoop's CLI policy. Explicit None disables a default;
        // Ellipsis selects the process-dependent default without a mutable preset.
        let cross_section = generation_type == GenerationType::CrossSection;
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
            let filter = if value.bind(py).is_instance_of::<PyEllipsis>() {
                GenerationFilter::LoopOneParticleIrreducible
            } else {
                GenerationFilter::MaxNumberOfBridges(value.extract::<usize>(py)?)
            };
            inner = inner.with_graph_filter(filter);
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
            inner = inner.numerator_grouping(value.inner);
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
/// >>> result = process.generate_diagrams(progress=report)
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
    /// interaction_filters, numerators, selection, grouping_preparation,
    /// grouping_samples, grouping_comparison, grouping, complete or cancelled.
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
/// >>> result = process.generate_diagrams()
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
                PyFeynmanDiagram::from(diagram.clone())
                    ._repr_html_(py)
                    .map(|html| {
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

/// Run generation with callbacks and signals on the calling Python thread.
#[allow(clippy::too_many_arguments)]
fn run_generation(
    py: Python<'_>,
    model: Arc<Model>,
    process: Process,
    generation_type: GenerationType,
    mut options: GenerationOptions,
    progress: Option<Py<PyAny>>,
    filter: Option<Py<PyAny>>,
) -> PyResult<PyGenerationResult> {
    py.check_signals()?;
    let automatic = progress.as_ref().is_some_and(|value| {
        value
            .bind(py)
            .extract::<String>()
            .is_ok_and(|s| s == "auto")
    });
    let progress = if automatic { None } else { progress };
    for (name, callback) in [("progress", &progress), ("filter", &filter)] {
        if callback
            .as_ref()
            .is_some_and(|callback| !callback.bind(py).is_callable())
        {
            return Err(PyTypeError::new_err(if name == "progress" {
                "progress must be 'auto', callable, or None".to_owned()
            } else {
                "filter must be callable".to_owned()
            }));
        }
    }
    // A running Marimo notebook has already imported marimo. Avoid importing an
    // optional notebook dependency (and its side effects) in ordinary scripts.
    let marimo = py.import("sys")?.getattr("modules")?;
    let marimo = marimo.cast::<PyDict>()?.get_item("marimo")?;
    let marimo_progress = if automatic
        && let Some(marimo) = marimo.filter(|module| !module.is_none())
        && marimo.call_method0("running_in_notebook")?.is_truthy()?
    {
        let flush = marimo
            .getattr("output")?
            .getattr("_output")?
            .getattr("flush")?
            .unbind();
        let status = marimo.getattr("status")?;
        let context =
            status.call_method1("spinner", ("Generating diagrams", "Starting generation"))?;
        let indicator = context.call_method0("__enter__")?.unbind();
        let state = Arc::new(Mutex::new((
            context.unbind(),
            indicator,
            None::<usize>,
            0_usize,
        )));
        Some((status.unbind(), state, flush))
    } else {
        None
    };
    let display = marimo_progress
        .as_ref()
        .map(|(status, state, flush)| (status.clone_ref(py), state.clone(), flush.clone_ref(py)));
    let started = Instant::now();
    let interruption = Arc::new(Mutex::new(None));
    let snapshots = Arc::new(Mutex::new(VecDeque::<GenerationProgress>::new()));
    let delivery = Mutex::new((Instant::now(), None));
    let pending = snapshots.clone();
    let has_progress = progress.is_some() || display.is_some();
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
        if has_progress {
            Python::attach(|py| -> PyResult<()> {
                for inner in updates {
                    if let Some((status, state, flush)) = &display {
                        let title = match inner.stage {
                            "topologies" => "Enumerating topologies",
                            "topology_filters" => "Filtering topologies",
                            "interactions" => "Assigning interactions",
                            "interaction_filters" => "Filtering interactions",
                            "numerators" => "Constructing numerators",
                            "selection" => "Selecting diagrams",
                            "grouping_preparation" => "Preparing numerators for grouping",
                            "grouping_samples" => "Sampling numerators for grouping",
                            "grouping_comparison" => "Comparing numerators",
                            "grouping" => "Finalizing diagram groups",
                            "complete" => "Diagram generation complete",
                            "cancelled" => "Diagram generation cancelled",
                            stage => stage,
                        };
                        let counts = match (inner.stage, inner.total) {
                            ("complete" | "cancelled", _) => {
                                format!("{} diagrams retained", inner.completed)
                            }
                            (_, Some(total)) => format!("{} / {total} processed", inner.completed),
                            (_, None) => format!("{} processed", inner.completed),
                        };
                        let subtitle =
                            format!("{counts} · {:.1}s elapsed", started.elapsed().as_secs_f64());
                        let mut state = state.lock().unwrap();
                        let (context, indicator, total, completed) = &mut *state;
                        // Unknown or empty work uses a spinner. Reuse a bar when
                        // totals agree; a new stage can reset its count to zero.
                        let next_total = inner.total.filter(|total| *total > 0);
                        if next_total != *total {
                            context.call_method1(
                                py,
                                "__exit__",
                                (py.None(), py.None(), py.None()),
                            )?;
                            let kwargs = PyDict::new(py);
                            kwargs.set_item("title", title)?;
                            kwargs.set_item("subtitle", &subtitle)?;
                            kwargs.set_item("remove_on_exit", true)?;
                            let widget = if let Some(total) = next_total {
                                kwargs.set_item("total", total)?;
                                kwargs.set_item("show_rate", false)?;
                                kwargs.set_item("show_eta", false)?;
                                "progress_bar"
                            } else {
                                "spinner"
                            };
                            *context = status.call_method(py, widget, (), Some(&kwargs))?;
                            *indicator = context.call_method0(py, "__enter__")?;
                            *total = next_total;
                            *completed = 0;
                        }
                        if total.is_some() {
                            let increment = inner.completed as i128 - *completed as i128;
                            indicator.call_method1(py, "update", (increment, title, subtitle))?;
                        } else {
                            indicator.call_method1(py, "update", (title, subtitle))?;
                        }
                        *completed = inner.completed;
                        // Delivery is already coalesced above. Marimo's own
                        // throttle can otherwise drop a stage's only update.
                        flush.call0(py)?;
                    }
                    if let Some(callback) = &progress {
                        callback.call1(py, (PyGenerationProgress { inner },))?;
                    }
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
    let pdgs: Vec<_> = model.particles().iter().map(|p| p.pdg_code).collect();
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
                let worker = scope.spawn(|| match generation_type {
                    GenerationType::Amplitude => process.generate_diagrams(model.clone(), &options),
                    GenerationType::CrossSection => {
                        process.generate_cross_section(model.clone(), &options)
                    }
                });
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
        py.detach(|| match generation_type {
            GenerationType::Amplitude => process.generate_diagrams(model.clone(), &options),
            GenerationType::CrossSection => process.generate_cross_section(model.clone(), &options),
        })
    };
    let result = (|| {
        if let Some(error) = interruption.lock().unwrap().take() {
            return Err(error);
        }
        deliver_progress(true)?;
        py.check_signals()?;
        result
            .map(|inner| PyGenerationResult { inner })
            .map_err(error::generation)
    })();
    if let Some((_, state, _)) = marimo_progress {
        let state = state.lock().unwrap();
        let (context, indicator, total, _) = &*state;
        // Always close the context, including when a callback or Python signal
        // interrupted generation. Cleanup must not replace the original error.
        let cleanup = match &result {
            Ok(_) => context.call_method1(py, "__exit__", (py.None(), py.None(), py.None())),
            Err(error) => {
                let title = "Diagram generation stopped";
                let subtitle = error.to_string();
                let _ = if total.is_some() {
                    indicator.call_method1(py, "update", (0, title, subtitle))
                } else {
                    indicator.call_method1(py, "update", (title, subtitle))
                };
                context.call_method1(
                    py,
                    "__exit__",
                    (error.get_type(py), error.value(py), error.traceback(py)),
                )
            }
        };
        if result.is_ok() {
            cleanup?;
        }
    }
    result
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
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
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use pyo3::types::PyDict;

    use super::*;

    #[test]
    fn automatic_progress_tracks_marimo_lifecycle() {
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
import inspect
import sys
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

model = fk.Model.from_json(MODEL_JSON)
process = fk.Process(model, [1000], [1000, 1000]).with_loop_count(1, 1)
settings = dict(max_vertices=3, threads=2, vertex_allow=["V_3_SCALAR_000"])
context = MagicMock()
indicator = context.__enter__.return_value
marimo = SimpleNamespace(
    running_in_notebook=MagicMock(return_value=True),
    status=SimpleNamespace(
        spinner=MagicMock(return_value=context),
        progress_bar=MagicMock(return_value=context),
    ),
    output=SimpleNamespace(_output=SimpleNamespace(flush=MagicMock())),
)

for via_model in (True, False):
    def generate(**kwargs):
        if via_model:
            return fk.Process(model, [1000], [1000, 1000]).generate_diagrams(loops=1, **settings, **kwargs)
        return process.generate_diagrams(**settings, **kwargs)

    method = process.generate_diagrams
    assert inspect.signature(method).parameters["progress"].default == "auto"
    with patch.dict(sys.modules, {"marimo": None}):
        baseline = generate()
    with patch.dict(sys.modules, {"marimo": marimo}):
        context.reset_mock()
        marimo.status.spinner.reset_mock()
        marimo.status.progress_bar.reset_mock()
        marimo.running_in_notebook.return_value = False
        assert len(generate()) == len(baseline)
        marimo.status.spinner.assert_not_called()
        marimo.status.progress_bar.assert_not_called()
        marimo.running_in_notebook.return_value = True
        assert len(generate(progress=None)) == len(baseline)
        updates = []
        generate(progress=updates.append)
        assert updates[-1].stage == "complete"
        marimo.status.spinner.assert_not_called()
        marimo.status.progress_bar.assert_not_called()

        visible = []
        marimo.output._output.flush.side_effect = lambda: visible.append(indicator.update.call_args.args[-2])
        result = generate()
        marimo.status.spinner.assert_called()
        marimo.status.progress_bar.assert_called()
        for call in marimo.status.progress_bar.call_args_list:
            assert call.kwargs["total"] > 0
            assert call.kwargs["remove_on_exit"] is True
            assert call.kwargs["show_eta"] is False
            assert call.kwargs["show_rate"] is False
        context.__enter__.assert_called()
        context.__exit__.assert_called_with(None, None, None)
        assert context.__exit__.call_count == context.__enter__.call_count
        titles = list(dict.fromkeys(call.args[-2] for call in indicator.update.call_args_list))
        assert titles == ["Enumerating topologies", "Filtering topologies",
            "Assigning interactions", "Filtering interactions", "Constructing numerators",
            "Selecting diagrams", "Preparing numerators for grouping",
            "Sampling numerators for grouping", "Comparing numerators",
            "Finalizing diagram groups", "Diagram generation complete"]
        # Every update must be flushed, even rapid transitions that Marimo's
        # own 150 ms throttle would suppress until another update arrives.
        assert visible == [call.args[-2] for call in indicator.update.call_args_list]
        assert indicator.update.call_args.args[-1].startswith(f"{len(result)} diagrams retained")
        assert all("elapsed" in call.args[-1] for call in indicator.update.call_args_list)
        assert any(" / " in call.args[-1] for call in indicator.update.call_args_list)

        context.reset_mock()
        empty = generate(progress="auto", filter=lambda graph, n: False)
        assert len(empty) == 0
        assert indicator.update.call_args.args[-1].startswith("0 diagrams retained")
        context.__exit__.assert_called_with(None, None, None)
        assert context.__exit__.call_count == context.__enter__.call_count

        context.reset_mock()
        token = fk.CancellationToken()
        token.cancel()
        result = generate(cancellation_token=token)
        assert not result.report.completed
        assert indicator.update.call_args.args[-2] == "Diagram generation cancelled"
        context.__exit__.assert_called_with(None, None, None)
        assert context.__exit__.call_count == context.__enter__.call_count

        for exception in (RuntimeError("filter failed"), KeyboardInterrupt()):
            context.reset_mock()
            def fail(*args):
                raise exception
            try:
                generate(filter=fail)
            except BaseException as caught:
                assert caught is exception
            else:
                raise AssertionError("generation exception was lost")
            assert context.__exit__.call_count == context.__enter__.call_count
            assert context.__exit__.call_args.args[:2] == (type(exception), exception)
            assert indicator.update.call_args.args[-2] == "Diagram generation stopped"

        context.reset_mock()
        exception = RuntimeError("spinner failed")
        indicator.update.side_effect = exception
        try:
            generate()
        except RuntimeError as caught:
            assert caught is exception
        else:
            raise AssertionError("spinner exception was lost")
        assert context.__exit__.call_count == context.__enter__.call_count
        indicator.update.side_effect = None

        for invalid in ("invalid", True, 1):
            context.reset_mock()
            try:
                generate(progress=invalid)
            except TypeError:
                pass
            else:
                raise AssertionError("invalid progress was accepted")
            context.__enter__.assert_not_called()
"#,
            )
            .unwrap();
            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
    }

    #[test]
    fn generation_defaults_leave_numerator_grouping_opt_in() {
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
settings = dict(max_vertices=3, vertex_allow=["V_3_SCALAR_000"])
explicit = dict(
    maximum_bridges=0, allow_self_loops=True,
    self_energy=fk.SelfEnergyFilterOptions(), tadpoles=fk.TadpoleFilterOptions(),
    zero_snails=fk.SnailFilterOptions(), numerator_grouping=None,
)
for via_model in (True, False):
    def generate(incoming, outgoing, loops=0, **kwargs):
        if via_model:
            return fk.Process(model, incoming, outgoing).generate_diagrams(loops=loops, **settings, **kwargs)
        process = fk.Process(model, incoming, outgoing).with_loop_count(loops, loops)
        return process.generate_diagrams(**settings, **kwargs)

    stages = []
    implicit = generate([1000], [1000, 1000], loops=1, progress=lambda p: stages.append(p.stage))
    configured = generate([1000], [1000, 1000], loops=1, **explicit)
    assert implicit.report.completed and len(implicit) > 0
    assert all(len(group.members) == 1 for group in implicit.groups)
    assert not any(stage.startswith("grouping_") for stage in stages)
    stages.clear()
    grouped = generate([1000], [1000, 1000], loops=1,
        numerator_grouping=fk.NumeratorGrouping("up_to_scalar"),
        progress=lambda p: stages.append(p.stage))
    assert grouped.report.completed and len(grouped) > 0
    assert "grouping_preparation" in stages
    assert [d.to_json() for d in implicit.diagrams] == [d.to_json() for d in configured.diagrams]

    # Automatic 1PI filtering applies only to loops. Explicit bridge limits still
    # apply to trees, while None removes the restriction at every loop order.
    tree = generate([1000, 1000], [1000, 1000], numerator_grouping=None)
    assert len(tree) > 0
    assert len(generate([1000, 1000], [1000, 1000], maximum_bridges=0)) == 0
    unrestricted = generate([1000, 1000], [1000, 1000], maximum_bridges=None, numerator_grouping=None)
    assert len(unrestricted) > 0
    assert len(unrestricted.groups) == len(unrestricted)
    assert all(len(group.members) == 1 for group in unrestricted.groups)

    vacuum = generate([], [], loops=2, numerator_grouping=None)
    configured_vacuum = generate([], [], loops=2, numerator_grouping=None,
        maximum_bridges=0, allow_self_loops=True, self_energy=None,
        tadpoles=None, zero_snails=None, factorized_loop_topologies_count_range=(1, 1))
    assert [d.to_json() for d in vacuum.diagrams] == [d.to_json() for d in configured_vacuum.diagrams]

assert fk.Process(model, [1000], [1000]).symmetrizes_final is None
for generate in (fk.Process.generate_diagrams, fk.Process.generate_amplitude, fk.Process.generate_cross_section):
    parameters = inspect.signature(generate).parameters
    assert parameters["maximum_bridges"].default is Ellipsis
    assert parameters["allow_self_loops"].default is True
    assert parameters["numerator_grouping"].default is None
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
process = fk.Process(model, [1000], [1000, 1000]).with_loop_count(1, 1)
settings = dict(max_vertices=3, threads=2, vertex_allow=["V_3_SCALAR_000"])
caller = threading.get_ident()

for via_model in (True, False):
    def generate(**callbacks):
        if via_model:
            return fk.Process(model, [1000], [1000, 1000]).generate_diagrams(loops=1, **settings, **callbacks)
        return process.generate_diagrams(**settings, **callbacks)

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
    stages = ["topologies", "topology_filters", "interactions", "interaction_filters", "numerators", "selection", "grouping_preparation", "grouping_samples", "grouping_comparison", "grouping", "complete"]
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
process = fk.Process(model, [1000], [1000, 1000]).with_loop_count(6, 6)
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
                    fk.Process(model, [1000], [1000, 1000]).generate_diagrams(loops=6, **settings)
                else:
                    process.generate_diagrams(**settings)
            except expected_error as error:
                assert str(error) == "generation interrupted"
            else:
                raise AssertionError("generation ignored its Python signal handler")
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
finally:
    signal.signal(signal.SIGALRM, previous_handler)

# Interruption belongs to one call; subsequent generation must still complete.
result = fk.Process(model, [1000], [1000, 1000]).generate_diagrams(loops=0, **settings)
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

mixed_process = fk.Process(model,
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
gluon_process = fk.Process(sm_model, [gluon], [gluon.antiparticle])
assert gluon_process.incoming[0].pdg == 21
assert gluon_process.outgoing_alternatives[0][0].pdg == 21

bottom = sm_model.particle_by_pdg(5)
bottom_process = fk.Process(sm_model, [bottom], [bottom.antiparticle])
assert bottom_process.incoming[0].pdg == 5
assert bottom_process.outgoing_alternatives[0][0].pdg == -5

process = fk.Process(model,
    [by_name, by_pdg],
    [fk.ParticleSelector.by_name("out"), fk.ParticleSelector.by_pdg(-1)],
).with_loop_count(1, 2)
assert all(isinstance(selector, fk.ParticleSelector) for selector in process.incoming)
assert process.incoming == [by_name, by_pdg]
assert process.outgoing_alternatives[0][0].name == "out"
assert process.outgoing_alternatives[0][1].pdg == -1

round_tripped = fk.Process(model, process.incoming, process.outgoing_alternatives[0]).with_loop_count(*process.loop_count)
assert round_tripped.incoming == process.incoming
assert round_tripped.outgoing_alternatives == process.outgoing_alternatives
assert round_tripped.loop_count == process.loop_count

cross_section = fk.Process(model, [particle, by_pdg], [by_name])
cross_section = cross_section.with_final_state_alternatives(
    [[particle.antiparticle], [by_name], [by_pdg], ["scalar_2"], [1002]],
)
cross_section = cross_section.with_symmetrization(
    initial=True,
    final_state=True,
    left_right=True,
    external_fermions=True,
)
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
    assert len(fk.Process(model, [particle], [particle, particle]).generate_diagrams(particle_veto=particle_veto)) == 0

foreign_model = fk.Model.from_json(MODEL_JSON.replace("scalar_0", "foreign_scalar_0"))
foreign_particle = foreign_model.particle("foreign_scalar_0")
foreign_process = fk.Process(model, [foreign_particle], [foreign_particle.antiparticle])
assert foreign_process.incoming[0].pdg == 1000
assert foreign_process.outgoing_alternatives[0][0].pdg == 1000
"#,
            )
            .unwrap();
            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
    }
}
