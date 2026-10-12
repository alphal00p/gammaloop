//! Staged generation. Python keys are interned once per theory; native search and compact
//! outputs carry only integer ranks. Python tuples own the original keys and participate in GC.

use std::ops::Bound as Limit;
use std::sync::atomic::{AtomicBool, AtomicU64};
use std::sync::Arc;

use linnet::half_edge::generate::{self, Stub};
use pyo3::exceptions::PyOverflowError;
use pyo3::types::{PyBool, PyRange, PyRangeMethods, PyTuple};

use super::algorithms::automorphism_order;
use super::keys::PayloadKeys;
use super::*;

mod callbacks;
pub use callbacks::{CancelCheck, GenerationCallbacks};
mod reports;
#[doc(hidden)]
pub use reports::PyPlanReport;
mod views;

/// A field and its direction at a vertex: outgoing True, incoming False, undirected None.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "EdgeSignature")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PyEdgeSignature {
    data: Option<Py<PyAny>>,
    flow: Option<Flow>,
}

impl PyEdgeSignature {
    fn field(&self) -> PyResult<&Py<PyAny>> {
        self.data.as_ref().ok_or_else(cleared)
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyEdgeSignature {
    #[new]
    #[pyo3(signature=(data, direction=None))]
    fn new(
        #[gen_stub(override_type(type_repr = "_Key"))] data: Py<PyAny>,
        direction: Option<bool>,
    ) -> Self {
        Self {
            data: Some(data),
            flow: direction.map(|out| if out { Flow::Source } else { Flow::Sink }),
        }
    }
    fn flip(&self, py: Python<'_>) -> PyResult<Self> {
        Ok(Self {
            data: Some(self.field()?.clone_ref(py)),
            flow: self.flow.map(|flow| -flow),
        })
    }
    #[getter]
    #[gen_stub(override_return_type(type_repr = "_Key"))]
    fn data(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(self.field()?.clone_ref(py))
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
            self.field()?.bind(py).repr()?
        ))
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.data)
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.data = None;
    }
}

fn cleared() -> PyErr {
    PyReferenceError::new_err("generation object was cleared by garbage collection")
}

/// Immutable tuple ownership is visible to Python's collector; each wrapper owns its own
/// references to the tuples, instead of traversing one hidden Arc-held reference repeatedly.
struct Objects {
    fields: Py<PyTuple>,
    labels: Py<PyTuple>,
    graph_thread: Option<std::thread::ThreadId>,
}
impl Objects {
    fn clone_ref(&self, py: Python<'_>) -> Self {
        Self {
            fields: self.fields.clone_ref(py),
            labels: self.labels.clone_ref(py),
            graph_thread: self.graph_thread,
        }
    }
    fn check_graph_snapshots(&self) -> PyResult<()> {
        if self
            .graph_thread
            .is_none_or(|thread| thread == std::thread::current().id())
        {
            Ok(())
        } else {
            Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Graph snapshots require the calling thread; use Run(thread_safe=False)",
            ))
        }
    }
    fn field(&self, py: Python<'_>, rank: u32) -> PyResult<Py<PyAny>> {
        Ok(self.fields.bind(py).get_item(rank as usize)?.unbind())
    }
    fn label(&self, py: Python<'_>, rank: u32) -> PyResult<Py<PyAny>> {
        Ok(self.labels.bind(py).get_item(rank as usize)?.unbind())
    }
    fn signature(&self, py: Python<'_>, stub: Stub<&u32>) -> PyResult<PyEdgeSignature> {
        Ok(PyEdgeSignature {
            data: Some(self.field(py, *stub.field)?),
            flow: stub.flow,
        })
    }
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
            |edge| {
                field(edge)
                    .map(|f| self.field(py, f))
                    .transpose()?
                    .into_py_any(py)
            },
            |label| {
                label
                    .map(|l| self.label(py, l))
                    .transpose()?
                    .into_py_any(py)
            },
        )
    }
    fn skeleton_form(
        &self,
        py: Python<'_>,
        form: &linnet::half_edge::canon::CanonicalForm<Option<usize>, Option<u32>, Option<u32>>,
    ) -> PyResult<Py<PyTuple>> {
        let words = PyTuple::new(py, form.words())?;
        let nodes = PyTuple::new(py, form.node_keys())?;
        let edges = PyTuple::new(
            py,
            form.edge_keys()
                .iter()
                .map(|f| match f {
                    Some(f) => PyTuple::new(py, [true.into_py_any(py)?, self.field(py, *f)?]),
                    None => PyTuple::new(py, [false.into_py_any(py)?]),
                })
                .collect::<PyResult<Vec<_>>>()?,
        )?;
        let hedges = PyTuple::new(
            py,
            form.hedge_keys()
                .iter()
                .map(|l| match l {
                    Some(l) => PyTuple::new(py, [true.into_py_any(py)?, self.label(py, *l)?]),
                    None => PyTuple::new(py, [false.into_py_any(py)?]),
                })
                .collect::<PyResult<Vec<_>>>()?,
        )?;
        Ok(PyTuple::new(py, [words, nodes, edges, hedges])?.unbind())
    }
    fn traverse(&self, visit: &PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.fields)?;
        visit.call(&self.labels)
    }
}

/// Canonical ordering of the supplied top-level keys. Nested tuple terms used by PayloadKeys
/// are not themselves fields or labels and therefore do not occupy output-table entries.
fn intern(
    py: Python<'_>,
    kind: &'static str,
    values: &[Py<PyAny>],
) -> PyResult<(Vec<u32>, Py<PyTuple>)> {
    let mut keys = PayloadKeys::new(py, None, None, None);
    let ids = values
        .iter()
        .map(|value| keys.term(kind, value.bind(py)))
        .collect::<PyResult<Vec<_>>>()?;
    let ranks = keys.ranks()?;
    let mut ordered: Vec<_> = ids.iter().copied().enumerate().collect();
    ordered.sort_unstable_by_key(|&(position, id)| (ranks[id as usize], position));
    ordered.dedup_by_key(|entry| entry.1);
    let mut dense = BTreeMap::new();
    for (rank, &(_, id)) in ordered.iter().enumerate() {
        dense.insert(id, rank as u32);
    }
    let table = PyTuple::new(py, ordered.iter().map(|&(i, _)| values[i].clone_ref(py)))?.unbind();
    Ok((ids.into_iter().map(|id| dense[&id]).collect(), table))
}

struct Slot {
    label: Py<PyAny>,
    field: Py<PyAny>,
    flow: Option<Flow>,
}
impl Slot {
    fn new(py: Python<'_>, label: Py<PyAny>, stub: &PyEdgeSignature) -> PyResult<Self> {
        Ok(Self {
            label,
            field: stub.field()?.clone_ref(py),
            flow: stub.flow,
        })
    }
    fn clone_ref(&self, py: Python<'_>) -> Self {
        Self {
            label: self.label.clone_ref(py),
            field: self.field.clone_ref(py),
            flow: self.flow,
        }
    }
    fn traverse(&self, visit: &PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.label)?;
        visit.call(&self.field)
    }
}

/// Ordered external slots. Equal labels and signatures describe interchangeable legs.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Legs")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PyLegs {
    slots: Vec<Slot>,
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyLegs {
    #[new]
    #[pyo3(signature=(slots=None))]
    fn new(
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "typing.Sequence[tuple[_Key, EdgeSignature]] | None", imports=("typing")))]
        slots: Option<Vec<(Py<PyAny>, Py<PyEdgeSignature>)>>,
    ) -> PyResult<Self> {
        Ok(Self {
            slots: slots
                .unwrap_or_default()
                .into_iter()
                .map(|(label, stub)| Slot::new(py, label, &stub.borrow(py)))
                .collect::<PyResult<_>>()?,
        })
    }
    fn leg<'py>(
        mut slf: PyRefMut<'py, Self>,
        #[gen_stub(override_type(type_repr = "_Key"))] label: Py<PyAny>,
        stub: &PyEdgeSignature,
    ) -> PyResult<PyRefMut<'py, Self>> {
        let slot = Slot::new(slf.py(), label, stub)?;
        slf.slots.push(slot);
        Ok(slf)
    }
    fn identical<'py>(
        mut slf: PyRefMut<'py, Self>,
        #[gen_stub(override_type(type_repr = "_Key"))] label: Py<PyAny>,
        stub: &PyEdgeSignature,
        multiplicity: usize,
    ) -> PyResult<PyRefMut<'py, Self>> {
        let slot = Slot::new(slf.py(), label, stub)?;
        for _ in 0..multiplicity {
            let copy = slot.clone_ref(slf.py());
            slf.slots.push(copy);
        }
        Ok(slf)
    }
    fn __len__(&self) -> usize {
        self.slots.len()
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        for slot in &self.slots {
            slot.traverse(&visit)?;
        }
        Ok(())
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.slots.clear();
    }
}

/// Bridge policy; dangling external legs never count as bridges.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Bridges", frozen, from_py_object)]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
#[derive(Clone)]
#[doc(hidden)]
pub struct PyBridges(generate::Bridges);

impl From<generate::Bridges> for PyBridges {
    fn from(value: generate::Bridges) -> Self {
        Self(value)
    }
}

impl From<PyBridges> for generate::Bridges {
    fn from(value: PyBridges) -> Self {
        value.0
    }
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyBridges {
    #[staticmethod]
    fn any() -> Self {
        Self(generate::Bridges::Any)
    }
    #[staticmethod]
    fn at_most(count: usize) -> Self {
        Self(generate::Bridges::AtMost(count))
    }
    #[staticmethod]
    fn one_pi() -> Self {
        Self(generate::Bridges::OnePI)
    }
    #[staticmethod]
    fn one_pi_in_loops() -> Self {
        Self(generate::Bridges::OnePIInLoops)
    }
    fn __repr__(&self) -> String {
        format!("Bridges({:?})", self.0)
    }
}

/// Restrictions decided by topology alone. Defaults are the native Structure defaults.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Structure", frozen, from_py_object)]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
#[derive(Clone)]
pub(crate) struct PyStructure(generate::Structure);
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyStructure {
    #[new]
    #[pyo3(signature=(*, bridges=None, self_loops=false, zero_flow_edges=false, max_vertices=None))]
    fn new(
        bridges: Option<&PyBridges>,
        self_loops: bool,
        zero_flow_edges: bool,
        max_vertices: Option<usize>,
    ) -> Self {
        Self(generate::Structure {
            bridges: bridges.map_or(generate::Bridges::Any, |b| b.0),
            self_loops,
            zero_flow_edges,
            max_vertices,
        })
    }
    #[staticmethod]
    fn one_pi() -> Self {
        Self(generate::Structure::one_pi())
    }
    #[getter]
    fn bridges(&self) -> PyBridges {
        PyBridges(self.0.bridges)
    }
    #[getter]
    fn self_loops(&self) -> bool {
        self.0.self_loops
    }
    #[getter]
    fn zero_flow_edges(&self) -> bool {
        self.0.zero_flow_edges
    }
    #[getter]
    fn max_vertices(&self) -> Option<usize> {
        self.0.max_vertices
    }
}

enum CountKind {
    Vertices(Vec<usize>),
    Propagators,
    Chords,
    ClosedLines,
}
/// A native rule/propagator/chord/closed-line count selector.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Countable")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PyCountable {
    kind: CountKind,
    fields: Vec<Py<PyAny>>,
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyCountable {
    #[staticmethod]
    fn vertices(rules: Vec<usize>) -> Self {
        Self {
            kind: CountKind::Vertices(rules),
            fields: Vec::new(),
        }
    }
    #[staticmethod]
    fn propagators(
        #[gen_stub(override_type(type_repr="typing.Sequence[_Key]", imports=("typing")))]
        fields: Vec<Py<PyAny>>,
    ) -> Self {
        Self {
            kind: CountKind::Propagators,
            fields,
        }
    }
    #[staticmethod]
    fn chords(
        #[gen_stub(override_type(type_repr="typing.Sequence[_Key]", imports=("typing")))]
        fields: Vec<Py<PyAny>>,
    ) -> Self {
        Self {
            kind: CountKind::Chords,
            fields,
        }
    }
    #[staticmethod]
    fn closed_lines(
        #[gen_stub(override_type(type_repr="typing.Sequence[_Key]", imports=("typing")))]
        fields: Vec<Py<PyAny>>,
    ) -> Self {
        Self {
            kind: CountKind::ClosedLines,
            fields,
        }
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        for field in &self.fields {
            visit.call(field)?;
        }
        Ok(())
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.fields.clear();
    }
}

#[derive(Clone)]
struct LoopOrders(Limit<usize>, Limit<usize>);
impl<'a, 'py> FromPyObject<'a, 'py> for LoopOrders {
    type Error = PyErr;
    fn extract(value: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        if value.is_instance_of::<PyBool>() {
            return Err(PyTypeError::new_err(
                "loops must be an int or a range, not a bool",
            ));
        }
        let Ok(range) = value.cast::<PyRange>() else {
            let n = value.extract()?;
            return Ok(Self(Limit::Included(n), Limit::Included(n)));
        };
        if range.step()? != 1 {
            return Err(PyValueError::new_err("a range of loop orders needs step 1"));
        }
        let bound = |n| {
            usize::try_from(n).map_err(|_| PyValueError::new_err("loop orders are non-negative"))
        };
        Ok(Self(
            Limit::Included(bound(range.start()?)?),
            Limit::Excluded(bound(range.stop()?)?),
        ))
    }
}
type WeightLimits = Vec<(String, (Option<i32>, Option<i32>))>;
type CountLimits = Vec<(Py<PyCountable>, (Option<usize>, Option<usize>))>;

/// A request snapshot. Later edits to the supplied Legs do not change it.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Query")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PyQuery {
    slots: Vec<Slot>,
    loops: LoopOrders,
    structure: generate::Structure,
    weights: WeightLimits,
    counts: CountLimits,
    hooks: Option<Py<callbacks::Hooks>>,
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyQuery {
    #[new]
    #[pyo3(signature=(legs, *, loops=None, structure=None, weights=None, counts=None, prune_skeletons=None, prune_dressings=None, select=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        legs: &PyLegs,
        #[gen_stub(override_type(type_repr = "int | range | None"))] loops: Option<LoopOrders>,
        structure: Option<&PyStructure>,
        #[gen_stub(override_type(type_repr = "dict[str, tuple[int | None, int | None]] | None"))]
        weights: Option<&Bound<'_, PyDict>>,
        #[gen_stub(override_type(type_repr="typing.Sequence[tuple[Countable, tuple[int | None, int | None]]] | None", imports=("typing")))]
        counts: Option<CountLimits>,
        #[gen_stub(override_type(type_repr="typing.Callable[[SkeletonView], bool] | None", imports=("typing")))]
        prune_skeletons: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[DressingView], bool] | None", imports=("typing")))]
        prune_dressings: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[DiagramView], bool] | None", imports=("typing")))]
        select: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        Ok(Self {
            slots: legs.slots.iter().map(|s| s.clone_ref(py)).collect(),
            loops: loops.unwrap_or(LoopOrders(Limit::Unbounded, Limit::Unbounded)),
            structure: structure.map_or_else(generate::Structure::default, |s| s.0.clone()),
            weights: weights
                .map(|d| {
                    d.iter()
                        .map(|(k, v)| Ok((k.extract()?, v.extract()?)))
                        .collect::<PyResult<_>>()
                })
                .transpose()?
                .unwrap_or_default(),
            counts: counts.unwrap_or_default(),
            hooks: callbacks::Hooks::new(py, prune_skeletons, prune_dressings, select)?,
        })
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        for slot in &self.slots {
            slot.traverse(&visit)?;
        }
        for (count, _) in &self.counts {
            visit.call(count)?;
        }
        visit.call(&self.hooks)
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.slots.clear();
        self.counts.clear();
        self.hooks = None;
    }
}

/// Ordered rules and reusable structural automata. Equal signatures remain distinct rules.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Theory")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PyTheory {
    inner: Arc<generate::Theory<u32>>,
    fields: Option<Py<PyTuple>>,
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyTheory {
    #[new]
    #[pyo3(signature=(vertices, *, weights=None))]
    fn new(
        py: Python<'_>,
        vertices: Vec<Vec<Py<PyEdgeSignature>>>,
        #[gen_stub(override_type(type_repr="dict[str, typing.Sequence[int]] | None", imports=("typing")))]
        weights: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Self> {
        let columns: Vec<(String, Vec<i32>)> = weights
            .map(|d| {
                d.iter()
                    .map(|(k, v)| Ok((k.extract()?, v.extract()?)))
                    .collect::<PyResult<_>>()
            })
            .transpose()?
            .unwrap_or_default();
        for (name, values) in &columns {
            if values.len() != vertices.len() {
                return Err(PyValueError::new_err(format!(
                    "weight {name:?} needs {} values, found {}",
                    vertices.len(),
                    values.len()
                )));
            }
        }
        let values = vertices
            .iter()
            .flatten()
            .map(|s| s.borrow(py).field().map(|f| f.clone_ref(py)))
            .collect::<PyResult<Vec<_>>>()?;
        let (ranks, fields) = intern(py, "field", &values)?;
        let mut ranks = ranks.into_iter();
        let signatures: Vec<Vec<_>> = vertices
            .iter()
            .map(|v| {
                v.iter()
                    .map(|s| Stub {
                        field: ranks.next().expect("one rank per stub"),
                        flow: s.borrow(py).flow,
                    })
                    .collect()
            })
            .collect();
        let theory = py
            .detach(move || {
                let mut builder = generate::Theory::builder()
                    .weights(columns.iter().map(|(name, _)| name.clone()));
                for (i, signature) in signatures.into_iter().enumerate() {
                    builder = builder.vertex(signature, columns.iter().map(|(_, v)| v[i]));
                }
                builder.build()
            })
            .map_err(|error| generation_error(py, error))?;
        Ok(Self {
            inner: Arc::new(theory),
            fields: Some(fields),
        })
    }
    #[getter]
    fn n_rules(&self) -> usize {
        self.inner.n_rules()
    }
    #[getter]
    fn weight_names(&self) -> Vec<String> {
        self.inner.weight_names().map(str::to_owned).collect()
    }
    #[pyo3(signature=(legs, *, loops=None, structure=None, weights=None, counts=None, prune_skeletons=None, prune_dressings=None, select=None))]
    #[allow(clippy::too_many_arguments)]
    fn query(
        &self,
        py: Python<'_>,
        legs: &PyLegs,
        #[gen_stub(override_type(type_repr = "int | range | None"))] loops: Option<LoopOrders>,
        structure: Option<&PyStructure>,
        #[gen_stub(override_type(type_repr = "dict[str, tuple[int | None, int | None]] | None"))]
        weights: Option<&Bound<'_, PyDict>>,
        #[gen_stub(override_type(type_repr="typing.Sequence[tuple[Countable, tuple[int | None, int | None]]] | None", imports=("typing")))]
        counts: Option<CountLimits>,
        #[gen_stub(override_type(type_repr="typing.Callable[[SkeletonView], bool] | None", imports=("typing")))]
        prune_skeletons: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[DressingView], bool] | None", imports=("typing")))]
        prune_dressings: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[DiagramView], bool] | None", imports=("typing")))]
        select: Option<Py<PyAny>>,
    ) -> PyResult<PyPlan> {
        PyPlan::new(
            py,
            self,
            &PyQuery::new(
                py,
                legs,
                loops,
                structure,
                weights,
                counts,
                prune_skeletons,
                prune_dressings,
                select,
            )?,
        )
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.fields)
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.fields = None;
    }
}

fn generation_error(py: Python<'_>, error: generate::GenerationError) -> PyErr {
    if matches!(error, generate::GenerationError::AutomorphismOverflow) {
        return PyOverflowError::new_err(error.to_string());
    }
    let exception = PyValueError::new_err(error.to_string());
    if let generate::GenerationError::Infeasible(report) = error {
        if let Err(error) = exception
            .value(py)
            .setattr("report", reports::PyPlanReport(report))
        {
            return error;
        }
    }
    exception
}

/// Transient execution controls. Cancellation and dispatch statistics are shared by executions
/// using this Run; callbacks are invoked on the Python thread advancing the execution.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Run")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PyRun {
    parallel: bool,
    check_heredity: bool,
    cancelled: Arc<AtomicBool>,
    progress: Option<Py<PyAny>>,
    round_trips: Arc<AtomicU64>,
}
impl PyRun {
    fn native(&self) -> generate::Run<'static> {
        let flag = Arc::clone(&self.cancelled);
        generate::Run::default()
            .parallel(self.parallel)
            .check_heredity(self.check_heredity)
            .deferred_progress(self.progress.is_some())
            .abort_check(move || flag.load(Ordering::Relaxed))
    }
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyRun {
    #[new]
    #[pyo3(signature=(*, parallel=false, check_heredity=false, progress=None))]
    fn new(
        py: Python<'_>,
        parallel: bool,
        check_heredity: bool,
        #[gen_stub(override_type(type_repr="typing.Callable[[Progress], bool] | None", imports=("typing")))]
        progress: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        if progress.as_ref().is_some_and(|p| !p.bind(py).is_callable()) {
            return Err(PyTypeError::new_err("progress must be callable"));
        }
        Ok(Self {
            parallel,
            check_heredity,
            progress,
            cancelled: Arc::new(AtomicBool::new(false)),
            round_trips: Arc::new(AtomicU64::new(0)),
        })
    }
    fn cancel(&self) {
        self.cancelled.store(true, Ordering::Relaxed);
    }
    #[getter]
    fn cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed)
    }
    /// Number of outer predicate-batch dispatches across executions using this Run.
    #[getter]
    fn callback_round_trips(&self) -> u64 {
        self.round_trips.load(Ordering::Relaxed)
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.progress)
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.cancel();
        self.progress = None;
    }
}

/// A compiled request, reusable across independent executions.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Plan")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PyPlan {
    inner: Arc<generate::Plan<u32, u32>>,
    objects: Option<Objects>,
    hooks: Option<Py<callbacks::Hooks>>,
}
impl PyPlan {
    fn objects(&self) -> PyResult<&Objects> {
        self.objects.as_ref().ok_or_else(cleared)
    }
    fn execution(
        &self,
        py: Python<'_>,
        run: Option<&PyRun>,
        parallel: Option<bool>,
        check_heredity: Option<bool>,
    ) -> PyResult<(Arc<Execution>, Controls)> {
        if run.is_some() && (parallel.is_some() || check_heredity.is_some()) {
            return Err(PyValueError::new_err(
                "pass either Run or execution keyword controls",
            ));
        }
        let default;
        let run = match run {
            Some(run) => run,
            None => {
                default = PyRun::new(
                    py,
                    parallel.unwrap_or(false),
                    check_heredity.unwrap_or(false),
                    None,
                )?;
                &default
            }
        };
        let controls = Controls {
            hooks: self.hooks.as_ref().map(|h| h.clone_ref(py)),
            progress: run.progress.as_ref().map(|p| p.clone_ref(py)),
            round_trips: Arc::clone(&run.round_trips),
            cancelled: {
                let token = Arc::clone(&run.cancelled);
                Arc::new(move || token.load(Ordering::Relaxed))
            },
        };
        Ok((
            Arc::new(Execution {
                plan: Arc::clone(&self.inner),
                run: run.native(),
            }),
            controls,
        ))
    }
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyPlan {
    #[new]
    fn new(py: Python<'_>, theory: &PyTheory, query: &PyQuery) -> PyResult<Self> {
        let fields = theory.fields.as_ref().ok_or_else(cleared)?;
        let mut keys = PayloadKeys::new(py, None, None, None);
        let mut field_ids = BTreeMap::new();
        for (rank, field) in fields.bind(py).iter().enumerate() {
            field_ids.insert(keys.term("field", &field)?, rank as u32);
        }
        let mut rank = |value: &Py<PyAny>| -> PyResult<u32> {
            let id = keys.term("field", value.bind(py))?;
            let next = field_ids.len() as u32;
            Ok(*field_ids.entry(id).or_insert(next))
        };
        let labels: Vec<_> = query.slots.iter().map(|s| s.label.clone_ref(py)).collect();
        let (label_ids, labels) = intern(py, "label", &labels)?;
        let mut legs = generate::Legs::new();
        for (slot, label) in query.slots.iter().zip(label_ids) {
            legs = legs.leg(
                label,
                Stub {
                    field: rank(&slot.field)?,
                    flow: slot.flow,
                },
            );
        }
        let mut native = generate::Query::new(&legs)
            .loops((query.loops.0, query.loops.1))
            .structure(query.structure.clone());
        for (name, (low, high)) in &query.weights {
            native = native.weight(
                name,
                (
                    low.map_or(Limit::Unbounded, Limit::Included),
                    high.map_or(Limit::Unbounded, Limit::Included),
                ),
            );
        }
        for (selector, (low, high)) in &query.counts {
            let selector = selector.borrow(py);
            let selected = selector
                .fields
                .iter()
                .map(&mut rank)
                .collect::<PyResult<Vec<_>>>()?;
            let count = match &selector.kind {
                CountKind::Vertices(rules) => generate::Countable::Vertices(rules.clone()),
                CountKind::Propagators => generate::Countable::Propagators(selected),
                CountKind::Chords => generate::Countable::Chords(selected),
                CountKind::ClosedLines => generate::Countable::ClosedLines(selected),
            };
            native = native.count(
                count,
                (
                    low.map_or(Limit::Unbounded, Limit::Included),
                    high.map_or(Limit::Unbounded, Limit::Included),
                ),
            );
        }
        if let Some(hooks) = &query.hooks {
            native = native.deferred_hooks(hooks.borrow(py).mask());
        }
        keys.ranks()?; // Preserve canonical custom-key validation for later query fields too.
        let inner = py
            .detach(|| generate::Plan::new(&theory.inner, &native))
            .map_err(|e| generation_error(py, e))?;
        Ok(Self {
            inner: Arc::new(inner),
            hooks: query.hooks.as_ref().map(|h| h.clone_ref(py)),
            objects: Some(Objects {
                fields: fields.clone_ref(py),
                labels,
                graph_thread: None,
            }),
        })
    }
    fn explain(&self) -> reports::PyPlanReport {
        reports::PyPlanReport(self.inner.explain().clone())
    }
    #[pyo3(signature=(run=None, *, parallel=None, check_heredity=None))]
    fn count(
        &self,
        py: Python<'_>,
        run: Option<&PyRun>,
        parallel: Option<bool>,
        check_heredity: Option<bool>,
    ) -> PyResult<reports::PyCounts> {
        let objects = self.objects()?;
        let (execution, controls) = self.execution(py, run, parallel, check_heredity)?;
        let mut search = DriverSearch::try_new(execution, |e| {
            e.plan.batched(&e.run, generate::batch::Output::Count)
        })
        .map_err(|e| generation_error(py, e))?;
        match search.advance(py, objects, &controls)? {
            generate::batch::Event::Done => search
                .with_dependent(|_, d| d.counts().cloned())
                .map(reports::PyCounts)
                .map_err(|error| generation_error(py, error)),
            _ => Err(pyo3::exceptions::PyRuntimeError::new_err(
                "count execution emitted an object",
            )),
        }
    }
    #[pyo3(signature=(run=None, *, parallel=None, check_heredity=None))]
    fn diagrams(
        &self,
        py: Python<'_>,
        run: Option<&PyRun>,
        parallel: Option<bool>,
        check_heredity: Option<bool>,
    ) -> PyResult<PyDiagrams> {
        let objects = self.objects()?.clone_ref(py);
        let (execution, controls) = self.execution(py, run, parallel, check_heredity)?;
        let search = DriverSearch::try_new(execution, |e| {
            e.plan.batched(&e.run, generate::batch::Output::Diagrams)
        })
        .map_err(|e| generation_error(py, e))?;
        Ok(PyDiagrams {
            search: Some(search),
            objects: Some(objects),
            controls,
        })
    }
    #[pyo3(signature=(run=None, *, parallel=None, check_heredity=None))]
    fn skeletons(
        &self,
        py: Python<'_>,
        run: Option<&PyRun>,
        parallel: Option<bool>,
        check_heredity: Option<bool>,
    ) -> PyResult<PySkeletons> {
        let objects = self.objects()?.clone_ref(py);
        let (execution, controls) = self.execution(py, run, parallel, check_heredity)?;
        let search = DriverSearch::try_new(execution, |e| {
            e.plan.batched(&e.run, generate::batch::Output::Skeletons)
        })
        .map_err(|e| generation_error(py, e))?;
        Ok(PySkeletons {
            search: Some(search),
            objects: Some(objects),
            controls,
        })
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        if let Some(o) = &self.objects {
            o.traverse(&visit)?;
        }
        visit.call(&self.hooks)
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.objects = None;
        self.hooks = None;
    }
}

struct Execution {
    plan: Arc<generate::Plan<u32, u32>>,
    run: generate::Run<'static>,
}
struct Controls {
    cancelled: CancelCheck,
    hooks: Option<Py<callbacks::Hooks>>,
    progress: Option<Py<PyAny>>,
    round_trips: Arc<AtomicU64>,
}
impl Controls {
    fn clone_ref(&self, py: Python<'_>) -> Self {
        Self {
            hooks: self.hooks.as_ref().map(|h| h.clone_ref(py)),
            progress: self.progress.as_ref().map(|p| p.clone_ref(py)),
            round_trips: Arc::clone(&self.round_trips),
            cancelled: Arc::clone(&self.cancelled),
        }
    }
    fn traverse(&self, visit: &PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.hooks)?;
        visit.call(&self.progress)
    }
    fn clear(&mut self) {
        self.hooks = None;
        self.progress = None;
    }
}
type Driver<'a> = generate::batch::Driver<'a, u32, u32>;
self_cell::self_cell!(
    struct DriverSearch {
        owner: Arc<Execution>,
        #[not_covariant]
        dependent: Driver,
    }
);
impl DriverSearch {
    fn abort(&mut self) {
        self.with_dependent_mut(|_, d| d.abort());
    }
    fn native<T: Send>(
        &mut self,
        py: Python<'_>,
        operation: impl FnOnce(&mut Driver<'_>) -> Result<T, generate::GenerationError> + Send,
    ) -> PyResult<T> {
        py.detach(|| {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                self.with_dependent_mut(|_, driver| operation(driver))
            }))
        })
        .map_err(|panic| {
            let message = panic
                .downcast_ref::<String>()
                .map(String::as_str)
                .or_else(|| panic.downcast_ref::<&str>().copied())
                .unwrap_or("generation panicked");
            pyo3::exceptions::PyRuntimeError::new_err(message.to_owned())
        })?
        .map_err(|error| generation_error(py, error))
    }
    fn advance(
        &mut self,
        py: Python<'_>,
        objects: &Objects,
        controls: &Controls,
    ) -> PyResult<generate::batch::Event<u32, u32>> {
        let result = (|| loop {
            py.check_signals()?;
            let event = self.native(py, |driver| driver.poll())?;
            match event {
                generate::batch::Event::Batch(requests) => {
                    let hooks = controls.hooks.as_ref().ok_or_else(cleared)?;
                    controls.round_trips.fetch_add(1, Ordering::Relaxed);
                    let Some(replies) =
                        callbacks::dispatch(py, requests, objects, hooks, &controls.cancelled)?
                    else {
                        self.abort();
                        return Ok(generate::batch::Event::Done);
                    };
                    self.native(py, |driver| driver.resolve(&replies))?;
                }
                generate::batch::Event::Progress(progress) => {
                    let callback = controls.progress.as_ref().ok_or_else(cleared)?;
                    if !callback
                        .call1(py, (reports::PyProgress(progress),))?
                        .extract::<bool>(py)?
                    {
                        self.with_dependent_mut(|_, driver| driver.stop_after_progress());
                    }
                }
                generate::batch::Event::Yield => {}
                event => return Ok(event),
            }
        })();
        if result.is_err() {
            self.abort();
        }
        result
    }
}

/// Complete untyped skeletons, streamed without decorating them.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Skeletons")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PySkeletons {
    search: Option<DriverSearch>,
    controls: Controls,
    objects: Option<Objects>,
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PySkeletons {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<PySkeleton>> {
        let Some(search) = &mut self.search else {
            return Ok(None);
        };
        match search.advance(
            py,
            self.objects.as_ref().ok_or_else(cleared)?,
            &self.controls,
        )? {
            generate::batch::Event::Skeleton(snapshot) => Ok(Some(PySkeleton {
                snapshot,
                execution: Some(Arc::clone(search.borrow_owner())),
                objects: self.objects.as_ref().map(|o| o.clone_ref(py)),
                controls: self.controls.clone_ref(py),
            })),
            generate::batch::Event::Done => Ok(None),
            _ => Err(pyo3::exceptions::PyRuntimeError::new_err(
                "skeleton execution emitted a diagram",
            )),
        }
    }
    #[getter]
    fn aborted(&self) -> bool {
        self.search
            .as_ref()
            .is_some_and(|s| s.with_dependent(|_, i| i.aborted()))
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        if let Some(o) = &self.objects {
            o.traverse(&visit)?;
        }
        self.controls.traverse(&visit)
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.search = None;
        self.objects = None;
        self.controls.clear();
    }
}

/// An immutable compact skeleton and its originating plan, independent of its source iterator.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Skeleton")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PySkeleton {
    snapshot: generate::SkeletonSnapshot<u32, u32>,
    execution: Option<Arc<Execution>>,
    objects: Option<Objects>,
    controls: Controls,
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PySkeleton {
    #[getter]
    fn loops(&self) -> usize {
        self.snapshot.view().loop_order()
    }
    fn automorphisms(&self) -> PyResult<u128> {
        automorphism_order(self.snapshot.automorphisms())
    }
    fn graph(&self, py: Python<'_>) -> PyResult<PyGraph> {
        self.objects
            .as_ref()
            .ok_or_else(cleared)?
            .graph(py, self.snapshot.graph(), |&field| field)
    }
    #[gen_stub(override_return_type(type_repr="tuple[typing.Any, ...]", imports=("typing")))]
    fn id(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        self.objects
            .as_ref()
            .ok_or_else(cleared)?
            .skeleton_form(py, self.snapshot.id())
    }
    fn view(&self, py: Python<'_>) -> PyResult<views::PySkeletonView> {
        Ok(views::PySkeletonView::from_snapshot(
            self.snapshot.clone(),
            self.objects.as_ref().ok_or_else(cleared)?.clone_ref(py),
        ))
    }
    /// Start fresh decoration without repeating topology search, using this execution's hooks.
    fn diagrams(&self, py: Python<'_>) -> PyResult<PyDiagrams> {
        let objects = self.objects.as_ref().ok_or_else(cleared)?.clone_ref(py);
        let execution = Arc::clone(self.execution.as_ref().ok_or_else(cleared)?);
        let snapshot = self.snapshot.clone();
        let search =
            DriverSearch::try_new(execution, |e| e.plan.batched_dressings(&snapshot, &e.run))
                .map_err(|e| generation_error(py, e))?;
        Ok(PyDiagrams {
            search: Some(search),
            objects: Some(objects),
            controls: self.controls.clone_ref(py),
        })
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        if let Some(o) = &self.objects {
            o.traverse(&visit)?;
        }
        self.controls.traverse(&visit)
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.execution = None;
        self.objects = None;
        self.controls.clear();
    }
}

/// A lazy stream of compact diagrams; no Python Graph or canonical form is created by next().
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Diagrams")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PyDiagrams {
    search: Option<DriverSearch>,
    controls: Controls,
    objects: Option<Objects>,
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyDiagrams {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<PyDiagram>> {
        let Some(search) = &mut self.search else {
            return Ok(None);
        };
        match search.advance(
            py,
            self.objects.as_ref().ok_or_else(cleared)?,
            &self.controls,
        )? {
            generate::batch::Event::Diagram(inner) => Ok(Some(PyDiagram {
                inner: Arc::new(inner),
                objects: self.objects.as_ref().map(|o| o.clone_ref(py)),
            })),
            generate::batch::Event::Done => Ok(None),
            _ => Err(pyo3::exceptions::PyRuntimeError::new_err(
                "diagram execution emitted a skeleton",
            )),
        }
    }
    #[getter]
    fn aborted(&self) -> bool {
        self.search
            .as_ref()
            .is_some_and(|s| s.with_dependent(|_, d| d.aborted()))
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        if let Some(o) = &self.objects {
            o.traverse(&visit)?;
        }
        self.controls.traverse(&visit)
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.search = None;
        self.objects = None;
        self.controls.clear();
    }
}

/// Owned compact output. Graph and canonical identity are computed only by their accessors.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Diagram")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(crate) struct PyDiagram {
    inner: Arc<generate::Diagram<u32, u32>>,
    objects: Option<Objects>,
}
impl PyDiagram {
    fn objects(&self) -> PyResult<&Objects> {
        self.objects.as_ref().ok_or_else(cleared)
    }
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyDiagram {
    #[getter]
    fn loops(&self) -> usize {
        self.inner.loops
    }
    fn symmetry_factor(&self) -> PyResult<u128> {
        automorphism_order(self.inner.automorphism_count)
    }
    fn rules(&self) -> Vec<u32> {
        self.inner.rules().to_vec()
    }
    fn edges(&self, py: Python<'_>) -> PyResult<Vec<(usize, usize, PyEdgeSignature)>> {
        let objects = self.objects()?;
        self.inner
            .edges()
            .map(|e| {
                Ok((
                    e.vertices[0],
                    e.vertices[1],
                    objects.signature(py, e.stub.expect("complete edge"))?,
                ))
            })
            .collect()
    }
    #[gen_stub(override_return_type(type_repr = "list[tuple[int, _Key, EdgeSignature]]"))]
    fn legs(&self, py: Python<'_>) -> PyResult<Vec<(usize, Py<PyAny>, PyEdgeSignature)>> {
        let objects = self.objects()?;
        self.inner
            .legs()
            .map(|(v, l, s)| Ok((v, objects.label(py, *l)?, objects.signature(py, s)?)))
            .collect()
    }
    fn skeleton(&self, py: Python<'_>) -> PyResult<views::PySkeletonView> {
        Ok(views::PySkeletonView::from_diagram(
            Arc::clone(&self.inner),
            self.objects()?.clone_ref(py),
        ))
    }
    fn graph(&self, py: Python<'_>) -> PyResult<PyGraph> {
        self.objects()?.graph(py, self.inner.graph(), |&f| Some(f))
    }
    /// Hashable native canonical words and key tables; no graph materialization is needed.
    #[gen_stub(override_return_type(type_repr="tuple[tuple[int, ...], tuple[int, ...], tuple[_Key, ...], tuple[tuple[typing.Any, ...], ...]]", imports=("typing")))]
    fn canonical_form(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        let o = self.objects()?;
        let form = self.inner.canonical_form();
        let words = PyTuple::new(py, form.words())?;
        let nodes = PyTuple::new(py, form.node_keys())?;
        let edges = PyTuple::new(
            py,
            form.edge_keys()
                .iter()
                .map(|&f| o.field(py, f))
                .collect::<PyResult<Vec<_>>>()?,
        )?;
        // Tag absent half-edge labels, preserving the distinction from a leg whose key is None.
        let hedges = PyTuple::new(
            py,
            form.hedge_keys()
                .iter()
                .map(|l| match l {
                    Some(l) => PyTuple::new(py, [true.into_py_any(py)?, o.label(py, *l)?]),
                    None => PyTuple::new(py, [false.into_py_any(py)?]),
                })
                .collect::<PyResult<Vec<_>>>()?,
        )?;
        Ok(PyTuple::new(py, [words, nodes, edges, hedges])?.unbind())
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        if let Some(o) = &self.objects {
            o.traverse(&visit)?;
        }
        Ok(())
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.objects = None;
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyEdgeSignature>()?;
    module.add_class::<PyLegs>()?;
    module.add_class::<PyBridges>()?;
    module.add_class::<PyStructure>()?;
    module.add_class::<PyCountable>()?;
    module.add_class::<PyQuery>()?;
    module.add_class::<PyTheory>()?;
    module.add_class::<PyRun>()?;
    module.add_class::<PyPlan>()?;
    module.add_class::<PyDiagrams>()?;
    module.add_class::<PySkeletons>()?;
    module.add_class::<PySkeleton>()?;
    module.add_class::<PyDiagram>()?;
    views::register(module)?;
    reports::register(module)
}
