//! Graph generation with Linnet's native generator: fields and leg labels are Python keys,
//! interned and ranked like canonization keys, and the generated graphs stream to Python.
use std::cell::Cell;
use std::ops::Bound as Limit;
use std::ptr::NonNull;
use std::sync::atomic::{AtomicBool, AtomicU32};
use std::sync::{Arc, Mutex};

use linnet::half_edge::generate::{self, GraphGenerator, GraphView, Stage, Stub, VertexSet};
use pyo3::types::{PyBool, PyRange, PyRangeMethods};

use super::algorithms::automorphism_order;
use super::keys::PayloadKeys;
use super::*;

/// A typed slot of a vertex signature or an external leg: a field and its direction relative to
/// its vertex, outgoing (True), incoming (False) or undirected (None). Two slots pair into an
/// edge when their fields are equal and both are undirected, or one is outgoing and the other
/// incoming. Fields are hashable keys, compared and ordered like `Graph.canonize` keys.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "EdgeSignature", frozen)]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PyEdgeSignature {
    data: Py<PyAny>,
    flow: Option<Flow>,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyEdgeSignature {
    /// A field `data` with an outgoing (True), incoming (False) or undirected (None) direction.
    #[new]
    #[pyo3(signature=(data, direction=None))]
    fn new(
        #[gen_stub(override_type(type_repr = "_Key"))] data: Py<PyAny>,
        direction: Option<bool>,
    ) -> Self {
        Self {
            data,
            flow: direction.map(|out| if out { Flow::Source } else { Flow::Sink }),
        }
    }

    /// The same field with the opposite direction.
    fn flip(&self, py: Python<'_>) -> Self {
        Self {
            data: self.data.clone_ref(py),
            flow: self.flow.map(|flow| -flow),
        }
    }

    /// The field, as given.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "_Key"))]
    fn data(&self, py: Python<'_>) -> Py<PyAny> {
        self.data.clone_ref(py)
    }

    #[getter]
    fn direction(&self) -> Option<bool> {
        self.flow.map(|flow| flow == Flow::Source)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let direction = match self.direction() {
            None => "None",
            Some(true) => "True",
            Some(false) => "False",
        };
        Ok(format!(
            "EdgeSignature({}, direction={direction})",
            self.data.bind(py).repr()?
        ))
    }
}

/// How far a `GenerationView` is built.
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyclass_enum
)]
#[pyclass(from_py_object, frozen, eq, eq_int, hash, name = "GenerationStage")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum PyGenerationStage {
    /// A connected untyped topology, partial or complete: vertices with their degrees, attached
    /// legs, edges without fields, and open slots that later vertices or legs take. Three or more
    /// legs with distinct labels and a common slot are open slots and unplaced legs until every
    /// vertex and edge is placed, then labelled one at a time.
    Topology,
    /// A complete topology with fields on some of its edges.
    Fields,
    /// A complete graph: every edge has a field and every vertex a signature.
    Complete,
}

impl From<Stage> for PyGenerationStage {
    fn from(stage: Stage) -> Self {
        match stage {
            Stage::Topology => Self::Topology,
            Stage::Fields => Self::Fields,
            Stage::Complete => Self::Complete,
        }
    }
}

/// The Python fields and leg labels, indexed by the ranks the generator sees.
struct Objects(Vec<Option<Py<PyAny>>>);

impl Objects {
    fn get(&self, py: Python<'_>, rank: u32) -> Py<PyAny> {
        self.0[rank as usize]
            .as_ref()
            .expect("a rank of a field or label")
            .clone_ref(py)
    }

    fn signature(&self, py: Python<'_>, stub: Stub<&u32>) -> PyEdgeSignature {
        PyEdgeSignature {
            data: self.get(py, *stub.field),
            flow: stub.flow,
        }
    }

    /// A generated or partial graph with signature indices on its nodes, the fields of its edges
    /// as Python objects (or None) and the labels on its legs.
    fn graph<E, V>(
        &self,
        py: Python<'_>,
        graph: &HedgeGraph<E, V, Option<u32>>,
        field: impl Fn(&E) -> Option<u32>,
    ) -> PyResult<PyGraph>
    where
        V: Clone + for<'py> IntoPyObject<'py>,
    {
        PyGraph::from_hedge_graph(
            py,
            graph,
            |node| node.clone().into_py_any(py),
            |edge| field(edge).map(|f| self.get(py, f)).into_py_any(py),
            |label| label.map(|l| self.get(py, l)).into_py_any(py),
        )
    }
}

/// A Python callback of the search, held where the garbage collector can visit and clear it.
#[derive(Default)]
struct Hook(Mutex<Option<Py<PyAny>>>);

impl Hook {
    fn get(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.0
            .lock()
            .unwrap()
            .as_ref()
            .map(|hook| hook.clone_ref(py))
    }

    /// Only a search step locks a hook, and it holds `&mut` on the iterator, so PyO3 skips the
    /// traversal of a locked hook already.
    fn traverse(&self, visit: &PyVisit<'_>) -> Result<(), PyTraverseError> {
        match self.0.try_lock() {
            Ok(hook) => visit.call(&*hook),
            Err(_) => Ok(()),
        }
    }

    fn clear(&self) {
        self.0.lock().unwrap().take();
    }
}

/// The state of the Python callbacks shared with the search: the callbacks, the first exception
/// they raised, which stops the generation, and a counter that spaces out signal checks.
#[derive(Default)]
struct Callbacks {
    filter: Hook,
    progress: Hook,
    error: Mutex<Option<PyErr>>,
    failed: AtomicBool,
    polls: AtomicU32,
}

/// Abort checks between two checks for Python signals.
const SIGNAL_POLLS: u32 = 1024;

impl Callbacks {
    fn failed(&self) -> bool {
        self.failed.load(Ordering::Relaxed)
    }

    /// The answer of a Python callback; an exception is kept for the caller, answers `false`
    /// and stops the generation, and no callback runs after one.
    fn ask(&self, callback: impl FnOnce(Python<'_>) -> PyResult<bool>) -> bool {
        if self.failed() {
            return false;
        }
        Python::attach(|py| {
            callback(py).unwrap_or_else(|error| {
                self.error.lock().unwrap().get_or_insert(error);
                self.failed.store(true, Ordering::Relaxed);
                false
            })
        })
    }

    /// Whether to stop: after an exception, or on a pending signal such as KeyboardInterrupt.
    fn interrupted(&self) -> bool {
        let poll = self.polls.fetch_add(1, Ordering::Relaxed);
        self.failed()
            || poll.is_multiple_of(SIGNAL_POLLS)
                && !self.ask(|py| py.check_signals().map(|()| true))
    }

    fn raise(&self) -> PyResult<()> {
        self.error.lock().unwrap().take().map_or(Ok(()), Err)
    }
}

/// A search state of `Graph.generate`, as its `filter_fn` sees it. Valid only during the call
/// that received it.
///
/// Vertices are numbered `0..n_vertices` and keep their numbers in every completion of the state:
/// the topology stage adds vertices with higher numbers, the later stages fields and signatures.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(unsendable, name = "GenerationView")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PyGenerationView {
    /// The viewed state; cleared before the filter call that borrows it returns.
    view: Cell<Option<NonNull<GraphView<'static, u32, u32>>>>,
    objects: Arc<Objects>,
}

impl PyGenerationView {
    /// Call `filter` on `view` as a Python object that is invalidated when the call returns.
    fn scoped(
        py: Python<'_>,
        view: &GraphView<'_, u32, u32>,
        objects: &Arc<Objects>,
        filter: &Py<PyAny>,
    ) -> PyResult<bool> {
        let view = NonNull::from(view).cast::<GraphView<'static, u32, u32>>();
        let object = Bound::new(
            py,
            Self {
                view: Cell::new(Some(view)),
                objects: objects.clone(),
            },
        )?;
        /// Clears the view when the call returns or unwinds.
        struct Expire<'a, 'py>(&'a Bound<'py, PyGenerationView>);
        impl Drop for Expire<'_, '_> {
            fn drop(&mut self) {
                self.0.borrow().view.set(None);
            }
        }
        let object = Expire(&object);
        filter.bind(py).call1((object.0,))?.extract()
    }

    /// `f` of the viewed state.
    fn with<R>(
        &self,
        f: impl for<'v> FnOnce(&'v GraphView<'v, u32, u32>) -> PyResult<R>,
    ) -> PyResult<R> {
        let view = self.view.get().ok_or_else(|| {
            PyReferenceError::new_err(
                "a generation view is only valid during the filter call that received it",
            )
        })?;
        // SAFETY: `scoped` clears the pointer before the state it borrows can change or end,
        // the class is unsendable, and `f` cannot keep a borrow beyond this call.
        f(unsafe { view.as_ref() })
    }

    /// `f` of the viewed state and a checked vertex of it.
    fn at<R>(
        &self,
        vertex: usize,
        f: impl for<'v> FnOnce(&'v GraphView<'v, u32, u32>, usize) -> R,
    ) -> PyResult<R> {
        self.with(|view| {
            Self::check_vertex(view, vertex)?;
            Ok(f(view, vertex))
        })
    }

    fn check_vertex(view: &GraphView<'_, u32, u32>, vertex: usize) -> PyResult<()> {
        if vertex < view.n_vertices() {
            Ok(())
        } else {
            Err(PyIndexError::new_err(format!(
                "vertex {vertex} is out of range for {} vertices",
                view.n_vertices()
            )))
        }
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyGenerationView {
    #[getter]
    fn stage(&self) -> PyResult<PyGenerationStage> {
        self.with(|view| Ok(view.stage().into()))
    }

    /// The loop order being generated.
    #[getter]
    fn loop_order(&self) -> PyResult<usize> {
        self.with(|view| Ok(view.loop_order()))
    }

    #[getter]
    fn n_vertices(&self) -> PyResult<usize> {
        self.with(|view| Ok(view.n_vertices()))
    }

    /// Independent cycles of the edges placed so far.
    #[getter]
    fn loops(&self) -> PyResult<usize> {
        self.with(|view| Ok(view.loops()))
    }

    /// Legs not attached to a vertex yet; zero from `GenerationStage.Fields` on.
    #[getter]
    fn unplaced_legs(&self) -> PyResult<usize> {
        self.with(|view| Ok(view.unplaced_legs()))
    }

    /// Every vertex, edge and leg is placed; always true from `GenerationStage.Fields` on.
    #[getter]
    fn is_topology_complete(&self) -> PyResult<bool> {
        self.with(|view| Ok(view.is_topology_complete()))
    }

    fn degree(&self, vertex: usize) -> PyResult<usize> {
        self.at(vertex, |view, v| view.degree(v))
    }

    /// Slots of `vertex` not taken by an edge or a leg yet; zero from `GenerationStage.Fields`
    /// on.
    fn open_stubs(&self, vertex: usize) -> PyResult<usize> {
        self.at(vertex, |view, v| view.open_stubs(v))
    }

    fn self_loops(&self, vertex: usize) -> PyResult<usize> {
        self.at(vertex, |view, v| view.self_loops(v))
    }

    /// Number of edges between the distinct vertices `u` and `v`.
    fn multiplicity(&self, u: usize, v: usize) -> PyResult<usize> {
        self.with(|view| {
            Self::check_vertex(view, u)?;
            Self::check_vertex(view, v)?;
            if u == v {
                return Err(PyValueError::new_err(
                    "multiplicity needs distinct vertices; use self_loops",
                ));
            }
            Ok(view.multiplicity(u, v))
        })
    }

    /// The labels of the legs attached to `vertex`, with their slots at it.
    #[gen_stub(override_return_type(type_repr = "list[tuple[_Key, EdgeSignature]]"))]
    fn legs(&self, py: Python<'_>, vertex: usize) -> PyResult<Vec<(Py<PyAny>, PyEdgeSignature)>> {
        self.at(vertex, |view, vertex| {
            view.legs(vertex)
                .map(|(&label, stub)| {
                    (
                        self.objects.get(py, label),
                        self.objects.signature(py, stub),
                    )
                })
                .collect()
        })
    }

    /// Every edge between vertices as `(u, v, slot)`, equal ends for a self-loop, with the field
    /// and direction at `u` once assigned: an outgoing slot points the edge from `u` to `v`. From
    /// `GenerationStage.Fields` on, the edges with a field come first.
    fn edges(&self, py: Python<'_>) -> PyResult<Vec<(usize, usize, Option<PyEdgeSignature>)>> {
        self.with(|view| {
            Ok(view
                .edges()
                .map(|edge| {
                    let [u, v] = edge.vertices;
                    (u, v, edge.stub.map(|stub| self.objects.signature(py, stub)))
                })
                .collect())
        })
    }

    /// The indices of the signatures that fit `vertex` once every slot of it has a field (from
    /// `GenerationStage.Fields` on), and in a complete graph the signature of `vertex`.
    fn rules(&self, vertex: usize) -> PyResult<Option<Vec<u32>>> {
        self.at(vertex, |view, vertex| {
            view.rules(vertex).map(<[u32]>::to_vec)
        })
    }

    /// The bridges with a final side, as `(outside, inside, side)`: `side` holds `inside` and has
    /// no open slots, so no later vertex or edge attaches to it and the bridge stays a bridge. A
    /// bridge of a complete topology is listed once per side.
    fn bridges(&self) -> PyResult<Vec<(usize, usize, Vec<usize>)>> {
        self.with(|view| {
            Ok(view
                .bridges()
                .into_iter()
                .map(|bridge| {
                    let [outside, inside] = bridge.ends;
                    (outside, inside, bridge.side.iter().collect())
                })
                .collect())
        })
    }

    /// Independent cycles of the subgraph on the non-empty connected vertex set `side`, such as
    /// the side of a bridge; raises ValueError for another set.
    fn loops_in(&self, side: Vec<usize>) -> PyResult<usize> {
        self.with(|view| {
            for &v in &side {
                Self::check_vertex(view, v)?;
            }
            let side = side.into_iter().collect::<VertexSet>();
            if !view.is_connected(side) {
                return Err(PyValueError::new_err(
                    "loops_in needs a non-empty set of vertices connected by placed edges",
                ));
            }
            Ok(view.loops_in(side))
        })
    }

    /// The state as a graph: the edges placed so far with their field once assigned (None
    /// before), the legs labelled as in generated graphs, and on each node the index of its
    /// signature once it is the only one that fits (None before). Open slots are left out.
    fn graph(&self, py: Python<'_>) -> PyResult<PyGraph> {
        self.with(|view| self.objects.graph(py, &view.graph(), |&field| field))
    }

    /// The complete graph as `Graph.generate` yields it; None before
    /// `GenerationStage.Complete`.
    fn complete_graph(&self, py: Python<'_>) -> PyResult<Option<PyGraph>> {
        self.with(|view| {
            view.complete_graph()
                .map(|graph| self.objects.graph(py, graph, |&field| Some(field)))
                .transpose()
        })
    }
}

self_cell::self_cell!(
    /// A generator with its running search.
    struct Search {
        owner: GraphGenerator<u32, u32>,
        #[not_covariant]
        dependent: Diagrams,
    }
);

type Diagrams<'g> = generate::Diagrams<'g, u32, u32>;

/// The graphs of `Graph.generate` with their automorphism-group orders, in depth-first search
/// order; each step searches up to the next graph.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "GeneratedGraphs")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PyGeneratedGraphs {
    search: Search,
    objects: Arc<Objects>,
    callbacks: Arc<Callbacks>,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyGeneratedGraphs {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        self.callbacks.filter.traverse(&visit)?;
        self.callbacks.progress.traverse(&visit)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.callbacks.filter.clear();
        self.callbacks.progress.clear();
    }

    #[gen_stub(override_return_type(type_repr = "tuple[Graph, int]"))]
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<(PyGraph, u128)>> {
        if self.callbacks.failed() {
            return Ok(None);
        }
        let search = &mut self.search;
        let diagram = py.detach(|| search.with_dependent_mut(|_, diagrams| diagrams.next()));
        self.callbacks.raise()?;
        let Some(diagram) = diagram else {
            return Ok(None);
        };
        let count = automorphism_order(diagram.automorphism_count)?;
        let graph = self
            .objects
            .graph(py, diagram.graph(), |&field| Some(field))?;
        Ok(Some((graph, count)))
    }
}

/// The `loops` of `Graph.generate`: a loop order or a `range` of them.
struct LoopOrders(Limit<usize>, Limit<usize>);

impl<'a, 'py> FromPyObject<'a, 'py> for LoopOrders {
    type Error = PyErr;

    fn extract(loops: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        if loops.is_instance_of::<PyBool>() {
            return Err(PyTypeError::new_err(
                "loops must be an int or a range, not a bool",
            ));
        }
        let Ok(range) = loops.cast::<PyRange>() else {
            let order = loops.extract()?;
            return Ok(Self(Limit::Included(order), Limit::Included(order)));
        };
        if range.step()? != 1 {
            return Err(PyValueError::new_err("a range of loop orders needs step 1"));
        }
        let bound = |value: isize| {
            usize::try_from(value)
                .map_err(|_| PyValueError::new_err("loop orders are non-negative"))
        };
        Ok(Self(
            Limit::Included(bound(range.start()?)?),
            Limit::Excluded(bound(range.stop()?)?),
        ))
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyGraph {
    /// Stream the connected graphs built from `vertex_signatures` and the external legs, one per
    /// isomorphism class, with their automorphism-group orders (the symmetry factors).
    ///
    /// `external_edges` pairs each leg label with its slot at the vertex the leg attaches to;
    /// legs with equal labels and slots are interchangeable. Fields and labels are hashable keys,
    /// compared and ordered like `canonize` keys. A generated graph has a node per vertex with the
    /// index of its signature as data, an edge per propagator with its field as data, directed
    /// from the outgoing to the incoming slot, and per leg a dangling half-edge with the leg
    /// label as data and the direction of the leg's slot; paired half-edges carry None.
    /// Signatures with identical slots are distinct vertex rules and give distinct graphs.
    ///
    /// `loops` is a loop order or a `range` of them. `loops` or `max_vertices` (legs are not
    /// vertices) must bound the search, and `max_vertices` is needed with signatures of two or
    /// fewer slots. `max_bridges=0` keeps one-particle-irreducible graphs; legs are never
    /// bridges. `allow_zero_flow_edges` allows bridges with no leg on one side.
    ///
    /// `filter_fn` prunes the search. It sees every search state as a `GenerationView`: partial
    /// untyped topologies (`GenerationStage.Topology`), then for each accepted complete topology
    /// partial field assignments (`Fields`) and complete graphs (`Complete`). A rejected state is
    /// dropped with everything that completes it, so the filter must be hereditary within each
    /// stage: when it rejects a partial topology, it rejects every complete topology completing
    /// it, and when it rejects a partial field assignment, every complete graph completing it.
    /// The output is then exactly the graphs whose complete topology and complete graph it
    /// accepts. Each state costs a call: test a condition at the stage that decides it and
    /// accept the earlier states, and keep verdicts on complete graphs independent of the vertex
    /// numbering.
    ///
    /// `progress_fn(loops, topologies, graphs)` runs after every complete topology and every
    /// yielded graph; returning False stops the generation. The search runs without the GIL,
    /// which the callbacks take while they run. An exception in a callback, or a signal such as
    /// KeyboardInterrupt, stops the generation and is raised by the step that saw it.
    #[staticmethod]
    #[pyo3(signature=(external_edges, vertex_signatures, *, loops=None, max_vertices=None, max_bridges=None, allow_self_loops=false, allow_zero_flow_edges=false, filter_fn=None, progress_fn=None))]
    #[allow(clippy::too_many_arguments)]
    fn generate<'py>(
        py: Python<'py>,
        #[gen_stub(override_type(type_repr = "typing.Sequence[tuple[_Key, EdgeSignature]]", imports = ("typing")))]
        external_edges: Vec<(Bound<'py, PyAny>, Bound<'py, PyEdgeSignature>)>,
        vertex_signatures: Vec<Vec<Bound<'py, PyEdgeSignature>>>,
        #[gen_stub(override_type(type_repr = "int | range | None"))] loops: Option<LoopOrders>,
        max_vertices: Option<usize>,
        max_bridges: Option<usize>,
        allow_self_loops: bool,
        allow_zero_flow_edges: bool,
        #[gen_stub(override_type(type_repr = "typing.Callable[[GenerationView], bool] | None", imports = ("typing")))]
        filter_fn: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "typing.Callable[[int, int, int], bool] | None", imports = ("typing")))]
        progress_fn: Option<Py<PyAny>>,
    ) -> PyResult<PyGeneratedGraphs> {
        let mut keys = PayloadKeys::new(py, None, None, None);
        let mut firsts = vec![];
        let mut key = |element, value: &Bound<'py, PyAny>| -> PyResult<u32> {
            let id = keys.term(element, value)?;
            firsts.push((id, value.clone().unbind()));
            Ok(id)
        };
        let labels = external_edges
            .iter()
            .map(|(label, _)| key("label", label))
            .collect::<PyResult<Vec<_>>>()?;
        let mut stub = |signature: &Bound<'py, PyEdgeSignature>| -> PyResult<Stub<u32>> {
            let signature = signature.get();
            Ok(Stub {
                field: key("field", signature.data.bind(py))?,
                flow: signature.flow,
            })
        };
        let signatures = vertex_signatures
            .iter()
            .map(|signature| signature.iter().map(&mut stub).collect())
            .collect::<PyResult<Vec<Vec<_>>>>()?;
        let externals = labels
            .into_iter()
            .zip(&external_edges)
            .map(|(label, (_, signature))| Ok((label, stub(signature)?)))
            .collect::<PyResult<Vec<_>>>()?;
        let rank = keys.ranks()?;
        let mut objects = Objects(rank.iter().map(|_| None).collect());
        for (id, object) in firsts {
            objects.0[rank[id as usize] as usize].get_or_insert(object);
        }
        let objects = Arc::new(objects);
        let ranked = |stub: Stub<u32>| Stub {
            field: rank[stub.field as usize],
            flow: stub.flow,
        };

        let mut generator = GraphGenerator::new(
            externals
                .into_iter()
                .map(|(label, stub)| (rank[label as usize], ranked(stub))),
            signatures
                .into_iter()
                .map(|signature| signature.into_iter().map(ranked).collect()),
        )
        .allow_self_loops(allow_self_loops)
        .allow_zero_flow_edges(allow_zero_flow_edges);
        if let Some(LoopOrders(min, max)) = loops {
            generator = generator.loops((min, max));
        }
        if let Some(max_vertices) = max_vertices {
            generator = generator.max_vertices(max_vertices);
        }
        if let Some(max_bridges) = max_bridges {
            generator = generator.max_bridges(max_bridges);
        }
        let (filter, progress) = (filter_fn.is_some(), progress_fn.is_some());
        let callbacks = Arc::new(Callbacks {
            filter: Hook(Mutex::new(filter_fn)),
            progress: Hook(Mutex::new(progress_fn)),
            ..Callbacks::default()
        });
        // a callback cleared by the garbage collector accepts: its iterator is being destroyed
        if filter {
            let (callbacks, objects) = (callbacks.clone(), objects.clone());
            generator = generator.filter(move |view| {
                callbacks.ask(|py| match callbacks.filter.get(py) {
                    Some(filter) => PyGenerationView::scoped(py, view, &objects, &filter),
                    None => Ok(true),
                })
            });
        }
        if progress {
            let callbacks = callbacks.clone();
            generator = generator.progress(move |p| {
                callbacks.ask(|py| match callbacks.progress.get(py) {
                    Some(progress) => progress
                        .call1(py, (p.loops, p.topologies, p.diagrams))?
                        .extract(py),
                    None => Ok(true),
                })
            });
        }
        let interrupts = callbacks.clone();
        generator = generator.abort_check(move || interrupts.interrupted());
        let search = Search::try_new(generator, |generator| generator.iter())
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok(PyGeneratedGraphs {
            search,
            objects,
            callbacks,
        })
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyEdgeSignature>()?;
    module.add_class::<PyGenerationStage>()?;
    module.add_class::<PyGenerationView>()?;
    module.add_class::<PyGeneratedGraphs>()?;
    Ok(())
}
