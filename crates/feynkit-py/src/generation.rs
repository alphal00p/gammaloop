//! Staged physics requests, compact results, and explicit numerator construction.
mod execution;
mod groups;
mod reports;
mod results;
mod specs;

pub use groups::{PyDiagramGroup, PyGroupMember};
pub use reports::{
    PyCounts, PyGenerationProgress, PyGenerationReport, PyLoopCounts, PyPlanReport, PyStage,
};
pub use results::{PyBuiltDiagramSet, PyDiagramSet, PyPreparedDiagram, PySkeleton};
pub use specs::*;

use crate::{
    error,
    graph::PyFeynmanDiagram,
    model::{PyModel, PyParticle, PyVertexRule},
    render_settings::PyDiagramStyle,
};
use feynkit_generator::{
    CancellationToken, DiagramGroup, GroupMember, ModelProcessExt, ParticleSelector, Process,
    Request, Run, VertexSelector,
};
use feynkit_graph::DiagramId;
use feynkit_model::Model;
use linnet::half_edge::generate::batch::{Hook, HookMask};
use linnet_py::{PyDiagramRender, PyRenderSettings};
use pyo3::{
    FromPyObject, IntoPyObjectExt, PyTraverseError, PyVisit,
    exceptions::{PyTypeError, PyValueError},
    prelude::*,
    types::{PyAny, PyBool, PyDict, PyModule},
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    PyStubType, TypeInfo,
    derive::{gen_stub_pyclass, gen_stub_pymethods},
};
use std::sync::Arc;
use symbolica::api::python::PythonExpression;

/// Select one model particle by name, signed PDG code, or stable model ID.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> selector = hep.ParticleSelector.by_name("e-")
/// >>> assert selector.name == "e-"
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ParticleSelector",
    module = "symbolica.community.hepkit",
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
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyParticleSelector {
    /// Select a particle by its model name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParticleSelector`` class example:
    ///
    /// >>> hep.ParticleSelector.by_name("e-").name
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
    /// Using the setup in the ``ParticleSelector`` class example:
    ///
    /// >>> hep.ParticleSelector.by_pdg(11).pdg
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
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParticleSelector`` class example:
    ///
    /// >>> assert hep.ParticleSelector.by_name("e-").name == "e-"
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
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParticleSelector`` class example:
    ///
    /// >>> assert hep.ParticleSelector.by_pdg(11).pdg == 11
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
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParticleSelector`` class example:
    ///
    /// >>> assert hep.ParticleSelector.by_name("e-").is_name
    #[getter]
    fn is_name(&self) -> bool {
        matches!(&self.inner, ParticleSelector::Name(_))
    }

    /// Report whether this selector identifies a particle by PDG code.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParticleSelector`` class example:
    ///
    /// >>> assert hep.ParticleSelector.by_pdg(11).is_pdg
    #[getter]
    fn is_pdg(&self) -> bool {
        matches!(&self.inner, ParticleSelector::Pdg(_))
    }
    /// Whether this selector retains an exact model fingerprint and particle ID.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParticleSelector`` class example:
    ///
    /// >>> model = hep.Model.standard_model()
    /// >>> selector = model.process([model.particle("e-")], []).incoming[0]
    /// >>> assert selector.is_id
    #[getter]
    fn is_id(&self) -> bool {
        matches!(&self.inner, ParticleSelector::Id { .. })
    }
    /// Resolve this selector in a model, validating any stored model identity.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParticleSelector`` class example:
    ///
    /// >>> model = hep.Model.standard_model()
    /// >>> assert electron.resolve(model).pdg_code == 11
    ///
    /// Parameters
    /// ----------
    /// model : Model
    ///     Model in which to resolve and validate the selector.
    fn resolve(&self, model: &PyModel) -> PyResult<PyParticle> {
        let id = self
            .inner
            .resolve(&model.inner)
            .map_err(|e| error::generation(e.into()))?;
        Ok(PyParticle::new(id, model.inner.clone()))
    }

    /// Format the selector as its particle name or signed PDG code.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParticleSelector`` class example:
    ///
    /// >>> str(hep.ParticleSelector.by_pdg(11))
    /// '11'
    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    /// Return a constructor-style representation of the selector.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParticleSelector`` class example:
    ///
    /// >>> print(hep.ParticleSelector.by_pdg(11))  # electron selector in a process
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
    /// Using the setup in the ``ParticleSelector`` class example:
    ///
    /// >>> hep.ParticleSelector.by_pdg(11) == 11
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
            Ok(Self(particle.get().selector()))
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
    fn split(self, ids: &mut Vec<DiagramId>, names: &mut Vec<String>) -> PyResult<()> {
        match self {
            Self::Diagram(diagram) => ids.push(diagram.inner.id().map_err(error::diagram)?),
            // IDs use their displayed width; short hexadecimal names such as D0 stay names.
            Self::Text(text)
                if text.len() == 32 && text.bytes().all(|byte| byte.is_ascii_hexdigit()) =>
            {
                ids.push(text.parse().expect("32 hexadecimal digits"));
            }
            Self::Text(text) => names.push(text),
        }
        Ok(())
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
            VertexInput::Vertex(vertex) => vertex.selector(),
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

/// Validated external states. Physics restrictions belong to a Request.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> process = model.process(["e-", "e+"], ["a", "a"])
/// >>> request = process.amplitude(loops=hep.Loops.total(0))
/// >>> assert request.count().complete
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Process",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyProcess {
    pub(crate) inner: Process,
    model: Arc<Model>,
}

impl PyProcess {
    pub(crate) fn from_model(
        model: &PyModel,
        incoming: Vec<SelectorInput>,
        outgoing: Vec<SelectorInput>,
    ) -> PyResult<Self> {
        let inner = model
            .inner
            .process(
                incoming.into_iter().map(ParticleSelector::from),
                outgoing.into_iter().map(ParticleSelector::from),
            )
            .map_err(error::generation)?;
        inner.validate_in(&model.inner).map_err(error::generation)?;
        Ok(Self {
            inner,
            model: model.inner.clone(),
        })
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyProcess {
    /// Model owning these external states.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Process`` class example:
    ///
    /// >>> assert process.model.name == model.name
    #[getter]
    fn model(&self) -> PyModel {
        self.model.clone().into()
    }
    /// Incoming particle selectors.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Process`` class example:
    ///
    /// >>> particles = [selector.resolve(model) for selector in process.incoming]
    #[getter]
    fn incoming(&self) -> Vec<PyParticleSelector> {
        self.inner
            .incoming()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }
    /// Accepted outgoing particle lists.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Process`` class example:
    ///
    /// >>> states = [[selector.resolve(model).name for selector in state] for state in process.outgoing_alternatives]
    #[getter]
    fn outgoing_alternatives(&self) -> Vec<Vec<PyParticleSelector>> {
        self.inner
            .outgoing_alternatives()
            .iter()
            .map(|s| s.iter().cloned().map(Into::into).collect())
            .collect()
    }
    /// Describe the current value.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Process`` class example:
    ///
    /// >>> print(process)
    fn __repr__(&self) -> String {
        let names = |s: &[ParticleSelector]| {
            s.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        };
        format!(
            "Process({:?}: [{}] -> {})",
            self.model.name(),
            names(self.inner.incoming()),
            self.inner
                .outgoing_alternatives()
                .iter()
                .map(|s| format!("[{}]", names(s)))
                .collect::<Vec<_>>()
                .join(" | ")
        )
    }
    /// Create a displayable blob with the process's physical incoming and outgoing particles.
    /// Alternative final states are displayed as separate schematics.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Process`` class example:
    ///
    /// >>> drawing = process.render()
    /// >>> drawing
    /// >>> svg = drawing.to_svg()
    ///
    /// Parameters
    /// ----------
    /// style : DiagramStyle or None, optional
    ///     Particle, momentum, and physical index presentation.
    /// config : RenderSettings or None, optional
    ///     Particle-label, layout and drawing overrides shared with Feynman diagrams.
    #[pyo3(signature = (*, config=None, style=None))]
    fn render(
        &self,
        config: Option<&PyRenderSettings>,
        style: Option<&PyDiagramStyle>,
    ) -> PyResult<PyDiagramRender> {
        let resolve = |state: &[ParticleSelector]| {
            state
                .iter()
                .map(|p| {
                    p.resolve(&self.model)
                        .map_err(|e| error::generation(e.into()))
                })
                .collect::<PyResult<Vec<_>>>()
        };
        let incoming = resolve(self.inner.incoming())?;
        let outgoing = self
            .inner
            .outgoing_alternatives()
            .iter()
            .map(|state| resolve(state))
            .collect::<PyResult<Vec<_>>>()?;
        let svg = crate::display::process_svg(&self.model, &incoming, &outgoing, config, style)?;
        let html = format!(
            "<figure class=\"feynkit-process\" style=\"max-width:100%;margin:.5rem 0\"><div style=\"width:360px;max-width:100%;overflow-x:auto\">{svg}</div><figcaption style=\"font-size:.85em;opacity:.75\">{}</figcaption></figure>",
            crate::display::escape_html(&self.__repr__())
        );
        Ok(PyDiagramRender::new(svg, html))
    }

    /// Display the process schematic in SVG-aware frontends.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Process`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(process)
    fn _repr_svg_(&self) -> PyResult<String> {
        self.render(None, None)?.to_svg()
    }

    /// Display the process blob, model and active restrictions in notebooks.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Process`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(process)
    fn _repr_html_(&self) -> PyResult<String> {
        self.render(None, None)?.to_html()
    }

    /// Return a process accepting any of the supplied final states.
    /// Generating amplitudes requires exactly one alternative; cross sections may include an empty
    /// alternative for a vacuum final state.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Process`` class example:
    ///
    /// >>> process = model.process([11, -11], [22, 22])
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
        let inner = self
            .inner
            .clone()
            .with_final_state_alternatives(
                alternatives
                    .into_iter()
                    .map(|state| state.into_iter().map(ParticleSelector::from)),
            )
            .map_err(error::process)?;
        inner.validate_in(&self.model).map_err(error::generation)?;
        Ok(Self {
            inner,
            ..self.clone()
        })
    }
    /// Create an amplitude (including vacuum) discovery request with explicit total loops.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Process`` class example:
    ///
    /// >>> request = process.amplitude(loops=hep.Loops.total(0))
    /// >>> compact = request.generate()
    ///
    /// Parameters
    /// ----------
    /// loops : Loops
    ///     Required finite total loop range for amplitudes or forward loop range for cross sections.
    /// orders : Orders or None, optional
    ///     Plain, squared, and applicable perturbative order constraints.
    /// restrict : Restrict or None, optional
    ///     Model-validated particle and vertex restrictions.
    /// topology : TopologyRules or None, optional
    ///     Explicit topology criteria; None keeps minimal unrestricted defaults.
    /// symmetrize : Symmetrize or None, optional
    ///     External-label permutation symmetries, disabled by default.
    /// prune_skeletons : callable or None, optional
    ///     Hereditary bool predicate receiving a scoped graph.SkeletonView.
    /// prune_dressings : callable or None, optional
    ///     Hereditary bool predicate receiving a scoped graph.DressingView.
    /// select : callable or None, optional
    ///     Arbitrary bool predicate receiving a complete graph.DiagramView.
    #[pyo3(signature=(*, loops, orders=None, restrict=None, topology=None, symmetrize=None, prune_skeletons=None, prune_dressings=None, select=None))]
    #[allow(clippy::too_many_arguments)]
    fn amplitude(
        &self,
        loops: &PyLoops,
        orders: Option<&PyOrders>,
        restrict: Option<&PyRestrict>,
        topology: Option<&PyTopologyRules>,
        symmetrize: Option<&PySymmetrize>,
        prune_skeletons: Option<Py<PyAny>>,
        prune_dressings: Option<Py<PyAny>>,
        select: Option<Py<PyAny>>,
    ) -> PyResult<PyRequest> {
        PyRequest::new(
            self,
            Request::amplitude(self.inner.clone(), loops.inner.clone()),
            orders,
            restrict,
            topology,
            None,
            symmetrize,
            prune_skeletons,
            prune_dressings,
            select,
        )
    }
    /// Create a cross-section discovery request with explicit forward loops.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Process`` class example:
    ///
    /// >>> forward = process.cross_section(loops=hep.Loops.forward(1),
    /// ...     cuts=hep.Cuts(blobs=hep.Range.exact(1), spectators=hep.Range.exact(0)))
    /// >>> compact = forward.generate()
    ///
    /// Parameters
    /// ----------
    /// loops : Loops
    ///     Required finite total loop range for amplitudes or forward loop range for cross sections.
    /// orders : Orders or None, optional
    ///     Plain, squared, and applicable perturbative order constraints.
    /// restrict : Restrict or None, optional
    ///     Model-validated particle and vertex restrictions.
    /// topology : TopologyRules or None, optional
    ///     Explicit topology criteria; None keeps minimal unrestricted defaults.
    /// cuts : Cuts or None, optional
    ///     Physical cut content and ranges; blobs and spectators default unrestricted.
    /// symmetrize : Symmetrize or None, optional
    ///     External-label permutation symmetries, disabled by default.
    /// prune_skeletons : callable or None, optional
    ///     Hereditary bool predicate receiving a scoped graph.SkeletonView.
    /// prune_dressings : callable or None, optional
    ///     Hereditary bool predicate receiving a scoped graph.DressingView.
    /// select : callable or None, optional
    ///     Arbitrary bool predicate receiving a complete graph.DiagramView.
    #[pyo3(signature=(*, loops, orders=None, restrict=None, topology=None, cuts=None, symmetrize=None, prune_skeletons=None, prune_dressings=None, select=None))]
    #[allow(clippy::too_many_arguments)]
    fn cross_section(
        &self,
        loops: &PyLoops,
        orders: Option<&PyOrders>,
        restrict: Option<&PyRestrict>,
        topology: Option<&PyTopologyRules>,
        cuts: Option<&PyCuts>,
        symmetrize: Option<&PySymmetrize>,
        prune_skeletons: Option<Py<PyAny>>,
        prune_dressings: Option<Py<PyAny>>,
        select: Option<Py<PyAny>>,
    ) -> PyResult<PyRequest> {
        PyRequest::new(
            self,
            Request::cross_section(self.inner.clone(), loops.inner.clone()),
            orders,
            restrict,
            topology,
            cuts,
            symmetrize,
            prune_skeletons,
            prune_dressings,
            select,
        )
    }
}

/// Immutable physics discovery criteria with replaceable scoped Python predicates.
/// Register hooks as decorators or constructor arguments; registering a hook again
/// replaces its previous callable and worker permission. The callable is returned unchanged.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> process = model.process(["e-", "e+"], ["a", "a"])
/// >>> request = process.amplitude(loops=hep.Loops.total(0))
/// >>> discovered = request.generate(run=hep.Run(threads=1))
/// >>> assert discovered.report.completed and len(discovered) > 0
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "Request", module = "symbolica.community.hepkit")]
pub struct PyRequest {
    inner: Request,
    model: Arc<Model>,
    prune_skeletons: Option<Py<PyAny>>,
    prune_dressings: Option<Py<PyAny>>,
    select: Option<Py<PyAny>>,
    worker_hooks: HookMask,
    check_heredity: bool,
}
impl PyRequest {
    #[allow(clippy::too_many_arguments)]
    fn new(
        process: &PyProcess,
        mut inner: Request,
        orders: Option<&PyOrders>,
        restrict: Option<&PyRestrict>,
        topology: Option<&PyTopologyRules>,
        cuts: Option<&PyCuts>,
        symmetrize: Option<&PySymmetrize>,
        prune_skeletons: Option<Py<PyAny>>,
        prune_dressings: Option<Py<PyAny>>,
        select: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        if let Some(v) = orders {
            inner = inner.with_orders(v.inner.clone());
        }
        if let Some(v) = restrict {
            inner = inner.with_restrict(v.inner.clone());
        }
        if let Some(v) = topology {
            inner = inner.with_topology(v.inner.clone());
        }
        if let Some(v) = cuts {
            inner = inner.with_cuts(v.inner.clone());
        }
        if let Some(v) = symmetrize {
            inner = inner.with_symmetrize(v.inner);
        }
        Python::attach(|py| -> PyResult<()> {
            for callback in [&prune_skeletons, &prune_dressings, &select]
                .into_iter()
                .flatten()
            {
                Self::validate_hook(callback.bind(py))?;
            }
            Ok(())
        })?;
        inner
            .validate_in(&process.model)
            .map_err(error::generation)?;
        Ok(Self {
            inner,
            model: process.model.clone(),
            prune_skeletons,
            prune_dressings,
            select,
            worker_hooks: HookMask::default(),
            check_heredity: false,
        })
    }

    fn validate_hook(callback: &Bound<'_, PyAny>) -> PyResult<()> {
        if callback.is_callable() {
            Ok(())
        } else {
            Err(PyTypeError::new_err(
                "generation predicates must be callable",
            ))
        }
    }

    fn register_hook(
        slf: Py<PyRequest>,
        py: Python<'_>,
        hook: Hook,
        callback: Option<Py<PyAny>>,
        thread_safe: bool,
    ) -> PyResult<Py<PyAny>> {
        if let Some(callback) = callback {
            Self::validate_hook(callback.bind(py))?;
            let mut request = slf.try_borrow_mut(py)?;
            match hook {
                Hook::Skeleton => {
                    request.prune_skeletons = Some(callback.clone_ref(py));
                    request.worker_hooks.skeletons = thread_safe;
                }
                Hook::Dressing => {
                    request.prune_dressings = Some(callback.clone_ref(py));
                    request.worker_hooks.dressings = thread_safe;
                }
                Hook::Select => {
                    request.select = Some(callback.clone_ref(py));
                    request.worker_hooks.select = thread_safe;
                }
            }
            Ok(callback)
        } else {
            // functools owns parameterized decorator binding and participates in Python GC.
            let name = match hook {
                Hook::Skeleton => "prune_skeletons",
                Hook::Dressing => "prune_dressings",
                Hook::Select => "select",
            };
            let kwargs = PyDict::new(py);
            kwargs.set_item("thread_safe", thread_safe)?;
            py.import("functools")?
                .getattr("partial")?
                .call((slf.bind(py).getattr(name)?,), Some(&kwargs))
                .map(Bound::unbind)
        }
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyRequest {
    /// Register a hereditary topology predicate and return the original callable.
    ///
    /// Use ``@request.prune_skeletons`` or ``@request.prune_skeletons(thread_safe=True)``.
    /// Registration replaces the previous topology hook. Worker permission applies only
    /// to this hook; unmarked hooks and progress stay on the invoking thread.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> @request.prune_skeletons
    /// ... def no_self_edges(view):
    /// ...     return all(view.self_loops(v) == 0 for v in range(view.n_vertices))
    /// >>> assert request.count().complete
    ///
    /// Parameters
    /// ----------
    /// callback : callable or None
    ///     Receives a scoped SkeletonView. None returns a parameterized decorator.
    /// thread_safe : bool
    ///     Permit the existing predicate worker; mutable graph snapshots require False.
    #[pyo3(signature=(callback=None, *, thread_safe=false))]
    fn prune_skeletons(
        slf: Py<PyRequest>,
        py: Python<'_>,
        callback: Option<Py<PyAny>>,
        thread_safe: bool,
    ) -> PyResult<Py<PyAny>> {
        Self::register_hook(slf, py, Hook::Skeleton, callback, thread_safe)
    }

    /// Register a hereditary field-assignment predicate and return the original callable.
    ///
    /// Use ``@request.prune_dressings`` or ``@request.prune_dressings(thread_safe=True)``.
    /// Registration replaces the previous dressing hook. The scoped DressingView expires
    /// after each call; retain an explicit graph snapshot only on the invoking thread.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> @request.prune_dressings(thread_safe=True)
    /// ... def no_gluon_ghosts(view):
    /// ...     return all(abs(stub.data.pdg_code) != 9000005
    /// ...                for _, _, stub in view.edges() if stub is not None)
    /// >>> assert request.count(run=hep.Run(threads=1)).complete
    ///
    /// Parameters
    /// ----------
    /// callback : callable or None
    ///     Receives a scoped DressingView. None returns a parameterized decorator.
    /// thread_safe : bool
    ///     Permit this hook to run on the existing predicate worker.
    #[pyo3(signature=(callback=None, *, thread_safe=false))]
    fn prune_dressings(
        slf: Py<PyRequest>,
        py: Python<'_>,
        callback: Option<Py<PyAny>>,
        thread_safe: bool,
    ) -> PyResult<Py<PyAny>> {
        Self::register_hook(slf, py, Hook::Dressing, callback, thread_safe)
    }

    /// Register a complete-assignment predicate and return the original callable.
    ///
    /// Use ``@request.select`` or ``@request.select(thread_safe=True)``. Registration
    /// replaces the previous selection hook. This terminal predicate need not be hereditary.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> @request.select
    /// ... def two_vertices(view):
    /// ...     return view.n_vertices == 2
    /// >>> selected = request.generate()
    /// >>> assert all(len(diagram.rules) == 2 for diagram in selected)
    ///
    /// Parameters
    /// ----------
    /// callback : callable or None
    ///     Receives a scoped DiagramView. None returns a parameterized decorator.
    /// thread_safe : bool
    ///     Permit this hook to run on the existing predicate worker.
    #[pyo3(signature=(callback=None, *, thread_safe=false))]
    fn select(
        slf: Py<PyRequest>,
        py: Python<'_>,
        callback: Option<Py<PyAny>>,
        thread_safe: bool,
    ) -> PyResult<Py<PyAny>> {
        Self::register_hook(slf, py, Hook::Select, callback, thread_safe)
    }

    /// Enable native hereditary-pruning checks for this Python request.
    /// Either this preference or ``Run(check_heredity=True)`` enables checking; False
    /// never disables an explicitly enabled Run. This preference is not persisted.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> request.check_heredity = True
    /// >>> assert request.check_heredity and request.count().complete
    #[getter]
    fn check_heredity(&self) -> bool {
        self.check_heredity
    }

    /// Set the request-local preference for native hereditary-pruning checks.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> request.check_heredity = True
    /// >>> checked = request.count(run=hep.Run(threads=1))
    /// >>> assert checked.complete
    ///
    /// Parameters
    /// ----------
    /// value : bool
    ///     Enable checking for this request; False does not override an enabled Run.
    #[setter]
    fn set_check_heredity(&mut self, value: bool) {
        self.check_heredity = value;
    }

    /// Validated external-state definition.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> assert request.process.incoming == process.incoming
    #[getter]
    fn process(&self) -> PyProcess {
        PyProcess {
            inner: self.inner.process().clone(),
            model: self.model.clone(),
        }
    }
    /// Required loop domain and inclusive limits.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> bounds = request.loops
    #[getter]
    fn loops(&self) -> PyLoops {
        PyLoops {
            inner: self.inner.loops().clone(),
        }
    }
    /// Coupling-order constraints.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> criteria = request.orders
    #[getter]
    fn orders(&self) -> PyOrders {
        PyOrders {
            inner: self.inner.orders().clone(),
        }
    }
    /// Particle and interaction restrictions.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> species = request.restrict
    #[getter]
    fn restrict(&self) -> PyRestrict {
        PyRestrict {
            inner: self.inner.restrictions().clone(),
        }
    }
    /// Topology criteria.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> assert request.topology.self_loops
    #[getter]
    fn topology(&self) -> PyTopologyRules {
        PyTopologyRules {
            inner: self.inner.topology().clone(),
        }
    }
    /// Physical cut criteria.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> assert request.cuts.blobs.maximum is None
    #[getter]
    fn cuts(&self) -> PyCuts {
        PyCuts {
            inner: self.inner.cuts().clone(),
        }
    }
    /// External-state symmetry choices.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> symmetry = request.symmetrize
    #[getter]
    fn symmetrize(&self) -> PySymmetrize {
        PySymmetrize {
            inner: self.inner.symmetrization(),
        }
    }
    /// Discover compact physical diagrams without numerators.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> compact = request.generate(run=hep.Run(threads=1))
    /// >>> assert len(compact) == request.count().diagrams
    ///
    /// Parameters
    /// ----------
    /// run : Run or None, optional
    ///     Execution controls; None uses Run defaults.
    #[pyo3(signature=(*,run=None))]
    fn generate(&self, py: Python<'_>, run: Option<&PyRun>) -> PyResult<PyDiagramSet> {
        let model = self.model.clone();
        let inner = execution::execute(py, model.clone(), Some(self), run, move |request, run| {
            request.expect("request execution").generate(model, &run)
        })?;
        Ok(PyDiagramSet {
            inner: Arc::new(inner),
            model: self.model.clone(),
        })
    }
    /// Count physical assignments and exact weights without constructing diagrams.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> counts = request.count(run=hep.Run(threads=1))
    /// >>> exact_weight = counts.sum_inv_aut
    ///
    /// Parameters
    /// ----------
    /// run : Run or None, optional
    ///     Execution controls; no FeynmanDiagram or numerator is constructed.
    #[pyo3(signature=(*,run=None))]
    fn count(&self, py: Python<'_>, run: Option<&PyRun>) -> PyResult<PyCounts> {
        let model = self.model.clone();
        execution::execute(py, model.clone(), Some(self), run, move |request, run| {
            request.expect("request execution").count(model, &run)
        })
        .map(PyCounts)
    }
    /// Explain active criteria without searching.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> report = request.explain()
    /// >>> criteria = report.criteria
    fn explain(&self, py: Python<'_>) -> PyResult<PyPlanReport> {
        let model = self.model.clone();
        execution::execute(py, model.clone(), Some(self), None, move |request, _| {
            request.expect("request execution").explain(&model)
        })
        .map(Into::into)
    }
    /// Describe the current value.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Request`` class example:
    ///
    /// >>> print(request)
    fn __repr__(&self) -> PyResult<String> {
        // The native serialization owns all physical defaults and omits callable handles.
        // Add the actual Python execution preferences, without reporting dormant native hooks.
        let mut effective = serde_json::to_value(&self.inner)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        effective["prune_skeletons"] = self.prune_skeletons.is_some().into();
        effective["prune_dressings"] = self.prune_dressings.is_some().into();
        effective["select"] = self.select.is_some().into();
        effective["thread_safe_hooks"] = serde_json::json!({
            "prune_skeletons": self.worker_hooks.skeletons,
            "prune_dressings": self.worker_hooks.dressings,
            "select": self.worker_hooks.select,
        });
        effective["check_heredity"] = self.check_heredity.into();
        Ok(format!("Request({effective})"))
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.prune_skeletons)?;
        visit.call(&self.prune_dressings)?;
        visit.call(&self.select)
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.prune_skeletons = None;
        self.prune_dressings = None;
        self.select = None;
        self.worker_hooks = HookMask::default();
    }
}

/// Shared cancellation flag, independent of request criteria.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> token = hep.CancellationToken()
/// >>> run = hep.Run(cancellation=token)
/// >>> token.cancel()
/// >>> assert run.cancellation.is_cancelled
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CancellationToken",
    module = "symbolica.community.hepkit",
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyCancellationToken {
    inner: CancellationToken,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCancellationToken {
    /// Create an uncancelled shared token.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CancellationToken`` class example:
    ///
    /// >>> token = hep.CancellationToken()
    #[new]
    fn new() -> Self {
        Self::default()
    }
    /// Request cooperative cancellation.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CancellationToken`` class example:
    ///
    /// >>> token.cancel()
    /// >>> assert token.is_cancelled
    fn cancel(&self) {
        self.inner.cancel();
    }
    /// Whether cancellation has been requested.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CancellationToken`` class example:
    ///
    /// >>> token.cancel()
    /// >>> assert token.is_cancelled
    #[getter]
    fn is_cancelled(&self) -> bool {
        self.inner.is_cancelled()
    }
}
/// Execution controls and an optional invoking-thread progress callback.
/// ``thread_safe=True`` permits typed predicates to run on one reused callback worker.
/// Worker predicates may inspect typed fields, but ``graph()`` and ``complete_graph()``
/// raise ``RuntimeError`` because mutable Graph snapshots require the calling thread.
/// Progress remains on the invoking thread; WebAssembly always runs callbacks inline.
/// This permission makes no parallel acceleration guarantee. Progress accepts None or
/// True to continue and False to cancel; other return values raise TypeError.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> run = hep.Run(threads=1)
/// >>> process = hep.Model.standard_model().process(["e-", "e+"], ["a", "a"])
/// >>> counts = process.amplitude(loops=hep.Loops.total(0)).count(run=run)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "Run", module = "symbolica.community.hepkit")]
pub struct PyRun {
    inner: Run,
    progress: Option<Py<PyAny>>,
    thread_safe: bool,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyRun {
    /// Construct execution controls.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Run`` class example:
    ///
    /// >>> run = hep.Run(threads=1, progress=lambda update: None)
    ///
    /// Parameters
    /// ----------
    /// threads : int or None, optional
    ///     Positive native worker limit; None uses the native default.
    /// parallel : bool, optional
    ///     Permit native parallel search, default True; WebAssembly runs inline.
    /// check_heredity : bool, optional
    ///     Check that pruning callbacks preserve hereditary rejection.
    /// cancellation : CancellationToken or None, optional
    ///     Shared cancellation flag; None creates an execution control token.
    /// progress : callable, str or None, optional
    ///     Calling-thread progress observer, "auto" for notebook display, or None.
    /// thread_safe : bool, optional
    ///     Opt into one reused predicate worker; Graph snapshots require False.
    #[new]
    #[pyo3(signature=(*,threads=None,parallel=true,check_heredity=false,cancellation=None,progress=None,thread_safe=false))]
    fn new(
        py: Python<'_>,
        threads: Option<usize>,
        parallel: bool,
        check_heredity: bool,
        cancellation: Option<&PyCancellationToken>,
        progress: Option<Py<PyAny>>,
        thread_safe: bool,
    ) -> PyResult<Self> {
        if threads == Some(0) {
            return Err(PyValueError::new_err("threads must be positive"));
        }
        if let Some(v) = &progress {
            if !v.bind(py).is_callable()
                && !v.bind(py).extract::<String>().is_ok_and(|s| s == "auto")
            {
                return Err(PyTypeError::new_err(
                    "progress must be callable, 'auto', or None",
                ));
            }
        }
        let mut inner = Run::default()
            .parallel(parallel)
            .check_heredity(check_heredity);
        if let Some(v) = threads {
            inner = inner.threads(v);
        }
        if let Some(v) = cancellation {
            inner = inner.cancellation_token(v.inner.clone());
        }
        Ok(Self {
            inner,
            progress,
            thread_safe,
        })
    }
    /// Shared cancellation token used by these controls.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Run`` class example:
    ///
    /// >>> run.cancellation.cancel()
    /// >>> assert run.cancellation.is_cancelled
    #[getter]
    fn cancellation(&self) -> PyCancellationToken {
        PyCancellationToken {
            inner: self.inner.cancellation().clone(),
        }
    }
    /// Describe the current value.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Run`` class example:
    ///
    /// >>> print(run)
    fn __repr__(&self, py: Python<'_>) -> String {
        let progress = self.progress.as_ref().map_or("None", |value| {
            if value
                .bind(py)
                .extract::<String>()
                .is_ok_and(|s| s == "auto")
            {
                "auto"
            } else {
                "callable"
            }
        });
        format!(
            "Run(controls={:?}, progress={progress}, thread_safe={})",
            self.inner, self.thread_safe
        )
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.progress)
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.progress = None;
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    specs::register(module.py(), module)?;
    reports::register(module)?;
    results::register(module)?;
    module.add_class::<PyParticleSelector>()?;
    module.add_class::<PyProcess>()?;
    module.add_class::<PyRequest>()?;
    module.add_class::<PyRun>()?;
    module.add_class::<PyCancellationToken>()?;
    module.add_class::<PyGroupMember>()?;
    module.add_class::<PyDiagramGroup>()?;
    Ok(())
}

#[cfg(test)]
mod tests;
