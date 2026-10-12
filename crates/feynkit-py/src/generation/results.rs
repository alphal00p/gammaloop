use super::*;
use crate::amplitude::PyAmplitude;
use feynkit_generator::{
    BuiltDiagramSet, DiagramSet, ExternalLabel, ExternalNode, PreparedDiagram, PreparedEdge,
    PreparedEndpoint,
};
use feynkit_model::ParticleId;
use linnet::half_edge::generate::SkeletonSnapshot;
use linnet_py::PyGraph;
use pyo3::{
    exceptions::{PyIndexError, PyOverflowError},
    types::PyTuple,
};

/// An external state from a prepared diagram. Obtain records from PreparedDiagram.legs or external edge endpoints.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> discovered = request.generate()
/// >>> diagram = discovered[0]
/// >>> attachment, leg, edge_index = diagram.legs[0]
/// >>> assert leg.index >= 0 and leg.particle.name
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "PreparedLeg", module = "symbolica.community.hepkit", frozen)]
pub struct PyPreparedLeg {
    inner: ExternalNode,
    model: Arc<Model>,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyPreparedLeg {
    /// External-state index retained by the physical process.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert leg.index >= 0
    #[getter]
    fn index(&self) -> usize {
        self.inner.index
    }
    /// External-state role, such as incoming or outgoing.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert leg.state in {"incoming", "outgoing"}
    #[getter]
    fn state(&self) -> &'static str {
        self.inner.state.as_str()
    }
    /// Particle record from the originating interaction model.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert abs(leg.particle.pdg_code) in {11, 22}
    #[getter]
    fn particle(&self) -> PyParticle {
        PyParticle::new(self.inner.particle, self.model.clone())
    }
}

/// A compact edge from PreparedDiagram.edges. Endpoints are vertex indices or PreparedLeg records.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> discovered = request.generate()
/// >>> diagram = discovered[0]
/// >>> edge = diagram.edges[0]
/// >>> assert abs(edge.particle.pdg_code) in {11, 22}
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "PreparedEdge", module = "symbolica.community.hepkit", frozen)]
pub struct PyPreparedEdge {
    inner: PreparedEdge,
    model: Arc<Model>,
}
impl PyPreparedEdge {
    fn endpoint<'py>(
        &self,
        py: Python<'py>,
        endpoint: PreparedEndpoint,
    ) -> PyResult<Bound<'py, PyAny>> {
        match endpoint {
            PreparedEndpoint::Vertex(index) => index.into_bound_py_any(py),
            PreparedEndpoint::External(inner) => PyPreparedLeg {
                inner,
                model: self.model.clone(),
            }
            .into_bound_py_any(py),
        }
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyPreparedEdge {
    /// The source endpoint: an internal vertex index or an external leg.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> endpoints = [(edge.source, edge.sink) for edge in diagram.edges]
    /// >>> assert len(endpoints) == len(diagram.edges)
    #[getter]
    #[gen_stub(override_return_type(type_repr = "int | PreparedLeg"))]
    fn source<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.endpoint(py, self.inner.source)
    }
    /// The sink endpoint: an internal vertex index or an external leg.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> endpoints = [(edge.source, edge.sink) for edge in diagram.edges]
    /// >>> assert len(endpoints) == len(diagram.edges)
    #[getter]
    #[gen_stub(override_return_type(type_repr = "int | PreparedLeg"))]
    fn sink<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.endpoint(py, self.inner.sink)
    }
    /// Particle species assigned to this edge.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert abs(edge.particle.pdg_code) in {11, 22}
    #[getter]
    fn particle(&self) -> PyParticle {
        PyParticle::new(self.inner.particle, self.model.clone())
    }
    /// Whether the particle assignment orients the edge.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> directed_particles = [edge.particle.pdg_code for edge in diagram.edges if edge.directed]
    /// >>> assert all(abs(pdg) == 11 for pdg in directed_particles)
    #[getter]
    fn directed(&self) -> bool {
        self.inner.directed
    }
}

/// A shared skeleton snapshot independent of an active search. Obtain it from a prepared diagram or DiagramSet.skeletons.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> discovered = request.generate()
/// >>> diagram = discovered[0]
/// >>> skeleton = diagram.skeleton
/// >>> assert skeleton.loops == 0
/// >>> identities = {skeleton.id}
/// >>> assert skeleton.id in identities
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "Skeleton", module = "symbolica.community.hepkit", frozen)]
pub struct PySkeleton {
    inner: SkeletonSnapshot<ParticleId, ExternalLabel>,
    model: Arc<Model>,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySkeleton {
    /// Hashable native canonical key, including model identity and complete key tables. This is a skeleton key, not a built Feynman-diagram identifier.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> identities = {skeleton.id}
    /// >>> assert diagram.skeleton.id in identities
    #[getter]
    fn id(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        let form = self.inner.id();
        let words = PyTuple::new(py, form.words())?;
        let nodes = PyTuple::new(py, form.node_keys())?;
        let edges = PyTuple::new(
            py,
            form.edge_keys()
                .iter()
                .map(|field| field.map(|id| id.index())),
        )?;
        let labels = PyTuple::new(
            py,
            form.hedge_keys().iter().map(|label| {
                label.map(|label| match label {
                    ExternalLabel::Exact { state, index } => (state.as_str(), Some(index)),
                    ExternalLabel::WithinState(state) => (state.as_str(), None),
                })
            }),
        )?;
        (
            self.model.fingerprint().to_string(),
            words,
            nodes,
            edges,
            labels,
        )
            .into_pyobject(py)
            .map(Bound::unbind)
            .map_err(Into::into)
    }
    /// Loop order of this skeleton.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert skeleton.loops == 0
    #[getter]
    fn loops(&self) -> usize {
        self.inner.view().loops()
    }
    /// Exact skeleton automorphism-group order. Raises OverflowError if the native order exceeds u128.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert skeleton.automorphisms >= 1
    #[getter]
    fn automorphisms(&self) -> PyResult<u128> {
        self.inner
            .automorphisms()
            .ok_or_else(|| PyOverflowError::new_err("skeleton automorphism order exceeds u128"))
    }
    /// Materialize an independent graph. Nodes carry optional engine signature ranks; assigned edges carry Particle objects and external labels retain state and index.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> graph = skeleton.graph
    /// >>> assert graph.n_edges > 0
    #[getter]
    fn graph(&self, py: Python<'_>) -> PyResult<PyGraph> {
        PyGraph::from_hedge_graph(
            py,
            self.inner.graph(),
            |rank| rank.into_bound_py_any(py).map(Bound::unbind),
            |field| {
                field
                    .map(|id| PyParticle::new(id, self.model.clone()))
                    .into_bound_py_any(py)
                    .map(Bound::unbind)
            },
            |label| {
                label
                    .map(|label| match label {
                        ExternalLabel::Exact { state, index } => (state.as_str(), Some(index)),
                        ExternalLabel::WithinState(state) => (state.as_str(), None),
                    })
                    .into_bound_py_any(py)
                    .map(Bound::unbind)
            },
        )
    }
    /// Describe the loop and vertex counts.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert "loops=0" in repr(skeleton)
    fn __repr__(&self) -> String {
        format!(
            "Skeleton(loops={}, vertices={})",
            self.inner.view().loops(),
            self.inner.view().n_vertices()
        )
    }
}

/// A physically prepared assignment without numerators, names, or canonical diagram IDs. Obtained by indexing or iterating DiagramSet.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> discovered = request.generate()
/// >>> diagram = discovered[0]
/// >>> assert diagram.loops == 0 and diagram.aut >= 1
/// >>> assert len(diagram.legs) == 4
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "PreparedDiagram",
    module = "symbolica.community.hepkit",
    frozen
)]
pub struct PyPreparedDiagram {
    set: Arc<DiagramSet>,
    index: usize,
    model: Arc<Model>,
}
impl PyPreparedDiagram {
    fn inner(&self) -> &PreparedDiagram {
        self.set.get(self.index).expect("retained compact index")
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyPreparedDiagram {
    /// Loop order of the prepared diagram.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert diagram.loops == 0
    #[getter]
    fn loops(&self) -> usize {
        self.inner().loops()
    }
    /// Automorphism-group order of this field and interaction assignment.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> from fractions import Fraction
    /// >>> assert Fraction(1, diagram.aut) > 0
    #[getter]
    fn aut(&self) -> u128 {
        self.inner().aut()
    }
    /// Assigned model vertex rules in native vertex order.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert len(diagram.rules) == diagram.graph.n_nodes
    #[getter]
    fn rules(&self) -> Vec<PyVertexRule> {
        self.inner()
            .rules()
            .map(|id| PyVertexRule::new(id, self.model.clone()))
            .collect()
    }
    /// Compact edge records in native edge order.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert len(diagram.edges) == diagram.graph.n_edges
    #[getter]
    fn edges(&self) -> Vec<PyPreparedEdge> {
        self.inner()
            .edges()
            .map(|inner| PyPreparedEdge {
                inner,
                model: self.model.clone(),
            })
            .collect()
    }
    /// Attachment vertex, physical external-state record, and edge index tuples in native order.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert len(diagram.legs) == 4
    /// >>> assert all(0 <= edge < len(diagram.edges) for vertex, leg, edge in diagram.legs)
    #[getter]
    fn legs(&self) -> Vec<(usize, PyPreparedLeg, usize)> {
        self.inner()
            .legs()
            .map(|(vertex, inner, edge)| {
                (
                    vertex,
                    PyPreparedLeg {
                        inner,
                        model: self.model.clone(),
                    },
                    edge,
                )
            })
            .collect()
    }
    /// Coupling powers of the assigned interaction rules.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert diagram.orders["QED"] == 2
    #[getter]
    fn orders(&self) -> std::collections::BTreeMap<String, usize> {
        self.inner().orders().clone()
    }
    /// Shared skeleton snapshot backing this assignment.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert diagram.skeleton.loops == diagram.loops
    #[getter]
    fn skeleton(&self) -> PySkeleton {
        PySkeleton {
            inner: self.inner().skeleton().clone(),
            model: self.model.clone(),
        }
    }
    /// Materialize an independent Graph with actual model rules and particles. Mutating that graph cannot alter this compact result.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> graph = diagram.graph
    /// >>> assert graph.n_edges == len(diagram.edges)
    #[getter]
    fn graph(&self, py: Python<'_>) -> PyResult<PyGraph> {
        PyGraph::from_hedge_graph(
            py,
            self.inner().graph(),
            |&id| {
                PyVertexRule::new(id, self.model.clone())
                    .into_bound_py_any(py)
                    .map(Bound::unbind)
            },
            |&id| {
                PyParticle::new(id, self.model.clone())
                    .into_bound_py_any(py)
                    .map(Bound::unbind)
            },
            |external| {
                external
                    .map(|inner| PyPreparedLeg {
                        inner,
                        model: self.model.clone(),
                    })
                    .into_bound_py_any(py)
                    .map(Bound::unbind)
            },
        )
    }
    /// Describe the loop count and automorphism order.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert "loops=0" in repr(diagram)
    fn __repr__(&self) -> String {
        format!(
            "PreparedDiagram(loops={}, aut={})",
            self.inner().loops(),
            self.inner().aut()
        )
    }
}

/// Compact physical discoveries returned by Request.generate(). Iteration does not construct Feynman diagrams; build() explicitly constructs numerators.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> discovered = request.generate()
/// >>> assert discovered.report.completed
/// >>> assert all(diagram.loops == 0 for diagram in discovered)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "DiagramSet", module = "symbolica.community.hepkit", frozen)]
pub struct PyDiagramSet {
    pub(super) inner: Arc<DiagramSet>,
    pub(super) model: Arc<Model>,
}
impl PyDiagramSet {
    fn at(&self, index: usize) -> PyPreparedDiagram {
        PyPreparedDiagram {
            set: self.inner.clone(),
            index,
            model: self.model.clone(),
        }
    }
    fn index(index: isize, len: usize) -> PyResult<usize> {
        let index = if index < 0 {
            len.checked_sub(index.unsigned_abs())
        } else {
            Some(index as usize)
        };
        index
            .filter(|&i| i < len)
            .ok_or_else(|| PyIndexError::new_err("diagram index out of range"))
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyDiagramSet {
    /// Number of retained items in this shared sequence.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert len(discovered) > 0
    fn __len__(&self) -> usize {
        self.inner.len()
    }
    /// Access one retained item by index. Negative indices count from the end; invalid indices raise IndexError.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert discovered[0].loops == 0
    ///
    /// Parameters
    /// ----------
    /// index : int
    ///     Sequence index, with Python negative-index semantics.
    fn __getitem__(&self, index: isize) -> PyResult<PyPreparedDiagram> {
        Ok(self.at(Self::index(index, self.inner.len())?))
    }
    /// Iterate the retained sequence without an eager Python list.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert sum(1 for item in discovered) == len(discovered)
    fn __iter__(&self) -> PyDiagramIterator {
        PyDiagramIterator {
            inner: self.inner.clone(),
            model: self.model.clone(),
            index: 0,
        }
    }
    /// Owned discovery report, including completion and retained counts.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert discovered.report.retained_count == len(discovered)
    #[getter]
    fn report(&self) -> PyGenerationReport {
        PyGenerationReport {
            inner: self.inner.report().clone(),
        }
    }
    /// Shared sequence of retained skeleton snapshots.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert all(skeleton.loops == 0 for skeleton in discovered.skeletons)
    #[getter]
    fn skeletons(&self) -> PySkeletons {
        PySkeletons {
            set: self.inner.clone(),
            model: self.model.clone(),
        }
    }
    /// Construct Feynman diagrams and numerators from these compact discoveries.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> built = discovered.build(hep.NumeratorOptions(naming=hep.Naming.named("sample")))
    /// >>> assert len(built) == len(discovered)
    ///
    /// Parameters
    /// ----------
    /// options : NumeratorOptions | None
    ///     Numerator, naming, zero-color filtering, and loop-momentum-basis settings.
    /// run : Run | None
    ///     Execution settings, cancellation, and progress callback for this build.
    #[pyo3(signature=(options=None,*,run=None))]
    fn build(
        &self,
        py: Python<'_>,
        options: Option<&PyNumeratorOptions>,
        run: Option<&PyRun>,
    ) -> PyResult<PyBuiltDiagramSet> {
        let options = options.map_or_else(Default::default, |v| v.inner.clone());
        let set = self.inner.clone();
        let inner = execution::execute(py, self.model.clone(), None, run, move |_, run| {
            set.build(&options, &run)
        })?;
        Ok(PyBuiltDiagramSet {
            inner: Arc::new(inner),
            model: self.model.clone(),
        })
    }
    /// Render discovery counts and completion without constructing diagrams or numerators.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> html = discovered._repr_html_()
    /// >>> assert "retained diagrams" in html
    fn _repr_html_(&self) -> String {
        self.report()._repr_html_()
    }
    /// Describe the compact result count and completion state.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert "completed=true" in repr(discovered)
    fn __repr__(&self) -> String {
        format!(
            "DiagramSet(diagrams={}, completed={})",
            self.inner.len(),
            self.inner.report().completed
        )
    }
}
/// An owning iterator returned by iter(discovered); its backing result remains alive while the iterator is retained.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> discovered = request.generate()
/// >>> iterator = iter(discovered)
/// >>> assert next(iterator).loops == 0
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "DiagramIterator", module = "symbolica.community.hepkit")]
pub struct PyDiagramIterator {
    inner: Arc<DiagramSet>,
    model: Arc<Model>,
    index: usize,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyDiagramIterator {
    /// Return this iterator without resetting its position.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert iter(iterator) is iterator
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    /// Return the next retained item, or raise StopIteration at exhaustion.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> iterator = iter(discovered)
    /// >>> assert next(iterator).loops == 0
    fn __next__(&mut self) -> Option<PyPreparedDiagram> {
        if self.index >= self.inner.len() {
            return None;
        }
        let result = PyPreparedDiagram {
            set: self.inner.clone(),
            index: self.index,
            model: self.model.clone(),
        };
        self.index += 1;
        Some(result)
    }
}
/// A shared sequence of skeleton snapshots obtained from DiagramSet.skeletons.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> discovered = request.generate()
/// >>> skeletons = discovered.skeletons
/// >>> assert len(skeletons) > 0
/// >>> assert all(skeleton.loops == 0 for skeleton in skeletons)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "Skeletons", module = "symbolica.community.hepkit", frozen)]
pub struct PySkeletons {
    set: Arc<DiagramSet>,
    model: Arc<Model>,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySkeletons {
    /// Number of retained items in this shared sequence.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert len(skeletons) > 0
    fn __len__(&self) -> usize {
        self.set.skeletons().len()
    }
    /// Access one retained item by index. Negative indices count from the end; invalid indices raise IndexError.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert skeletons[0].loops == 0
    ///
    /// Parameters
    /// ----------
    /// index : int
    ///     Sequence index, with Python negative-index semantics.
    fn __getitem__(&self, index: isize) -> PyResult<PySkeleton> {
        let i = PyDiagramSet::index(index, self.set.skeletons().len())?;
        Ok(PySkeleton {
            inner: self.set.skeletons()[i].clone(),
            model: self.model.clone(),
        })
    }
    /// Iterate the retained sequence without an eager Python list.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert sum(1 for item in skeletons) == len(skeletons)
    fn __iter__(&self) -> PySkeletonIterator {
        PySkeletonIterator {
            set: self.set.clone(),
            model: self.model.clone(),
            index: 0,
        }
    }
}
/// An owning iterator returned by iter(skeletons); its backing result remains alive while the iterator is retained.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> discovered = request.generate()
/// >>> skeletons = discovered.skeletons
/// >>> iterator = iter(skeletons)
/// >>> assert next(iterator).loops == 0
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "SkeletonIterator", module = "symbolica.community.hepkit")]
pub struct PySkeletonIterator {
    set: Arc<DiagramSet>,
    model: Arc<Model>,
    index: usize,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySkeletonIterator {
    /// Return this iterator without resetting its position.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert iter(iterator) is iterator
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    /// Return the next retained item, or raise StopIteration at exhaustion.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> iterator = iter(skeletons)
    /// >>> assert next(iterator).loops == 0
    fn __next__(&mut self) -> Option<PySkeleton> {
        let inner = self.set.skeletons().get(self.index)?.clone();
        self.index += 1;
        Some(PySkeleton {
            inner,
            model: self.model.clone(),
        })
    }
}

/// Explicitly constructed Feynman diagrams with provenance and retained squared-order interference policy. Obtain it through DiagramSet.build().
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> discovered = request.generate()
/// >>> built = discovered.build()
/// >>> assert len(built) == len(discovered)
/// >>> built.validate_groups()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "BuiltDiagramSet",
    module = "symbolica.community.hepkit",
    frozen
)]
pub struct PyBuiltDiagramSet {
    inner: Arc<BuiltDiagramSet>,
    model: Arc<Model>,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyBuiltDiagramSet {
    /// Number of retained items in this shared sequence.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert len(built) > 0
    fn __len__(&self) -> usize {
        self.inner.len()
    }
    /// Access one retained item by index. Negative indices count from the end; invalid indices raise IndexError.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> first = built[0]
    /// >>> last = built[-1]
    /// >>> assert len(built) > 0
    ///
    /// Parameters
    /// ----------
    /// index : int
    ///     Sequence index, with Python negative-index semantics.
    fn __getitem__(&self, index: isize) -> PyResult<PyFeynmanDiagram> {
        let i = PyDiagramSet::index(index, self.inner.len())?;
        Ok(self.inner.diagrams()[i].clone().into())
    }
    /// Iterate the retained sequence without an eager Python list.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert sum(1 for item in built) == len(built)
    fn __iter__(&self) -> PyBuiltDiagramIterator {
        PyBuiltDiagramIterator {
            inner: self.inner.clone(),
            index: 0,
        }
    }
    /// Shared sequence view of constructed diagrams, without an eager Python list.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert len(built.diagrams) == len(built)
    #[getter]
    fn diagrams(&self) -> PyBuiltDiagramSet {
        Self {
            inner: self.inner.clone(),
            model: self.model.clone(),
        }
    }
    /// Numerator-equivalence groups and their member weights.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> grouped = built.group(hep.Grouping("identical"))
    /// >>> grouped.validate_groups()
    /// >>> groups = grouped.groups
    #[getter]
    fn groups(&self) -> Vec<PyDiagramGroup> {
        self.inner
            .groups()
            .iter()
            .cloned()
            .map(|inner| PyDiagramGroup { inner })
            .collect()
    }
    /// Owned report from construction, selection, or grouping.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert built.report.completed
    #[getter]
    fn report(&self) -> PyGenerationReport {
        PyGenerationReport {
            inner: self.inner.report().clone(),
        }
    }
    /// Return the subset selected by IDs, names, vetoes, and forced cuts.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> empty = built.select(hep.Selection(select=[]))
    /// >>> assert len(empty) == 0
    ///
    /// Parameters
    /// ----------
    /// selection : Selection
    ///     Selection policy; explicit empty inclusion retains no diagrams.
    /// run : Run | None
    ///     Execution settings, cancellation, and progress callback for selection.
    #[pyo3(signature=(selection,*,run=None))]
    fn select(
        &self,
        py: Python<'_>,
        selection: &PySelection,
        run: Option<&PyRun>,
    ) -> PyResult<Self> {
        let set = self.inner.clone();
        let selection = selection.inner.clone();
        execution::execute(py, self.model.clone(), None, run, move |_, run| {
            set.select(&selection, &run)
        })
        .map(|inner| Self {
            inner: Arc::new(inner),
            model: self.model.clone(),
        })
    }
    /// Return numerator-equivalence groups under the requested comparison relation.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> grouped = built.group(hep.Grouping("identical"))
    /// >>> grouped.validate_groups()
    ///
    /// Parameters
    /// ----------
    /// grouping : Grouping
    ///     Numerator equivalence relation and comparison settings.
    /// run : Run | None
    ///     Execution settings, cancellation, and progress callback for grouping.
    #[pyo3(signature=(grouping,*,run=None))]
    fn group(&self, py: Python<'_>, grouping: &PyGrouping, run: Option<&PyRun>) -> PyResult<Self> {
        let set = self.inner.clone();
        let grouping = grouping.inner.clone();
        execution::execute(py, self.model.clone(), None, run, move |_, run| {
            set.group(&grouping, &run)
        })
        .map(|inner| Self {
            inner: Arc::new(inner),
            model: self.model.clone(),
        })
    }
    /// Check that group masters, members, and weights obey native invariants; raise on inconsistency.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> grouped = built.group(hep.Grouping("up_to_sign"))
    /// >>> grouped.validate_groups()
    fn validate_groups(&self) -> PyResult<()> {
        self.inner.validate_groups().map_err(error::generation)
    }
    /// Construct an amplitude while preserving resolved allowed interference pairs.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> amplitude = built.into_amplitude(hep.AmplitudeOptions(dimension=4))
    /// >>> squared = amplitude.squared()
    ///
    /// Parameters
    /// ----------
    /// options : AmplitudeOptions | None
    ///     Lorentz dimension and additional real expressions.
    #[pyo3(signature=(options=None))]
    fn into_amplitude(
        &self,
        py: Python<'_>,
        options: Option<&PyAmplitudeOptions>,
    ) -> PyResult<PyAmplitude> {
        let set = self.inner.clone();
        let options = options.map_or_else(Default::default, |v| v.inner.clone());
        py.detach(move || (*set).clone().into_amplitude(options))
            .map(|inner| PyAmplitude { inner })
            .map_err(|err| match err {
                feynkit_generator::GenerationError::Amplitude(err) => error::amplitude(*err),
                err => error::generation(err),
            })
    }
    /// Describe diagram count, group count, and completion state.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert "BuiltDiagramSet" in repr(built)
    fn __repr__(&self) -> String {
        format!(
            "BuiltDiagramSet(diagrams={}, groups={}, completed={})",
            self.inner.len(),
            self.inner.groups().len(),
            self.inner.report().completed
        )
    }
    /// Render a notebook preview of the constructed diagram collection.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert "Built diagrams" in built._repr_html_()
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        let subtitle = format!(
            "{} diagrams · {} groups · {}",
            self.inner.len(),
            self.inner.groups().len(),
            if self.inner.report().completed {
                "complete"
            } else {
                "incomplete"
            }
        );
        crate::display::collection_html(
            py,
            "Built diagrams",
            &subtitle,
            self.inner
                .diagrams()
                .iter()
                .cloned()
                .map(PyFeynmanDiagram::from),
            None,
            None,
            None,
            crate::display::PREVIEW_LIMIT,
        )
        .map(|(html, _)| html)
    }
}
/// An owning iterator returned by iter(built); its backing result remains alive while the iterator is retained.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> discovered = request.generate()
/// >>> built = discovered.build()
/// >>> iterator = iter(built)
/// >>> first = next(iterator)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "BuiltDiagramIterator", module = "symbolica.community.hepkit")]
pub struct PyBuiltDiagramIterator {
    inner: Arc<BuiltDiagramSet>,
    index: usize,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyBuiltDiagramIterator {
    /// Return this iterator without resetting its position.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert iter(iterator) is iterator
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    /// Return the next retained item, or raise StopIteration at exhaustion.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> iterator = iter(built)
    /// >>> first = next(iterator)
    fn __next__(&mut self) -> Option<PyFeynmanDiagram> {
        let result = self.inner.diagrams().get(self.index)?.clone().into();
        self.index += 1;
        Some(result)
    }
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyPreparedLeg>()?;
    module.add_class::<PyPreparedEdge>()?;
    module.add_class::<PyPreparedDiagram>()?;
    module.add_class::<PySkeleton>()?;
    module.add_class::<PyDiagramSet>()?;
    module.add_class::<PyDiagramIterator>()?;
    module.add_class::<PySkeletons>()?;
    module.add_class::<PySkeletonIterator>()?;
    module.add_class::<PyBuiltDiagramSet>()?;
    module.add_class::<PyBuiltDiagramIterator>()?;
    Ok(())
}
