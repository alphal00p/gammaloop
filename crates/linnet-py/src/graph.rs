use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use linnet::half_edge::builder::{HedgeData, HedgeGraphBuilder};
use linnet::half_edge::involution::{EdgeIndex, Flow, Hedge, HedgePair, Orientation};
use linnet::half_edge::swap::Swap;
use linnet::half_edge::NodeIndex;
use linnet::permutation::Permutation;
use pyo3::class::gc::{PyTraverseError, PyVisit};
use pyo3::exceptions::{PyIndexError, PyKeyError, PyReferenceError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyTuple, PyType};
use pyo3_stub_gen::derive::{
    gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pyfunction, gen_stub_pymethods,
};

use crate::dot::{self, PyDotCodec, PyGlobalData};
use crate::drawing::{
    copy_drawing, drawing_dict, PyEdgeDrawing, PyHalfEdgeDrawing, PyNodeDrawing, EDGE_FIELDS,
    HEDGE_FIELDS, NODE_FIELDS,
};
use crate::native_graph::{PyHedgeGraph, PyNodeStore};

static SPEC_ID: AtomicU64 = AtomicU64::new(0);

/// The directed role of a half-edge within its edge.
#[gen_stub_pyclass_enum]
#[pyclass(from_py_object, eq, eq_int, name = "Flow")]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PyFlow {
    Source,
    Sink,
}

impl From<PyFlow> for Flow {
    fn from(value: PyFlow) -> Self {
        match value {
            PyFlow::Source => Self::Source,
            PyFlow::Sink => Self::Sink,
        }
    }
}

impl From<Flow> for PyFlow {
    fn from(value: Flow) -> Self {
        match value {
            Flow::Source => Self::Source,
            Flow::Sink => Self::Sink,
        }
    }
}

/// How an edge's declared endpoints determine its logical direction.
#[gen_stub_pyclass_enum]
#[pyclass(from_py_object, eq, eq_int, name = "Orientation")]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum PyOrientation {
    #[default]
    Default,
    Reversed,
    Undirected,
}

impl From<PyOrientation> for Orientation {
    fn from(value: PyOrientation) -> Self {
        match value {
            PyOrientation::Default => Self::Default,
            PyOrientation::Reversed => Self::Reversed,
            PyOrientation::Undirected => Self::Undirected,
        }
    }
}

impl From<Orientation> for PyOrientation {
    fn from(value: Orientation) -> Self {
        match value {
            Orientation::Default => Self::Default,
            Orientation::Reversed => Self::Reversed,
            Orientation::Undirected => Self::Undirected,
        }
    }
}

/// A declarative node description consumed by `build()`.
#[gen_stub_pyclass]
#[pyclass(unsendable, name = "NodeSpec")]
pub struct PyNodeSpec {
    token: u64,
    name: Option<String>,
    data: Option<Py<PyAny>>,
    drawing: Option<Py<PyDict>>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyNodeSpec {
    #[getter]
    fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    #[getter]
    fn data(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.data
            .as_ref()
            .map(|data| data.clone_ref(py))
            .ok_or_else(|| PyReferenceError::new_err("node specification has been cleared"))
    }

    #[setter]
    fn set_data(&mut self, value: Py<PyAny>) {
        self.data = Some(value);
    }

    #[getter]
    fn drawing(&self, py: Python<'_>) -> PyResult<PyNodeDrawing> {
        self.drawing
            .as_ref()
            .map(|drawing| PyNodeDrawing::detached(py, drawing))
            .ok_or_else(|| PyReferenceError::new_err("node specification has been cleared"))
    }

    fn __repr__(&self) -> String {
        match &self.name {
            Some(name) => format!("node({name:?})"),
            None => "node()".to_string(),
        }
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.data)?;
        visit.call(&self.drawing)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.data = None;
        self.drawing = None;
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum EndpointRole {
    Source,
    Sink,
}

impl EndpointRole {
    fn flow(self) -> PyFlow {
        match self {
            Self::Source => PyFlow::Source,
            Self::Sink => PyFlow::Sink,
        }
    }
}

/// A declarative edge endpoint produced by `source()` or `sink()`.
#[gen_stub_pyclass]
#[pyclass(unsendable, name = "HalfEdgeSpec")]
pub struct PyHalfEdgeSpec {
    node: Option<Py<PyNodeSpec>>,
    role: EndpointRole,
    data: Option<Py<PyAny>>,
    drawing: Option<Py<PyDict>>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyHalfEdgeSpec {
    #[getter]
    fn flow(&self) -> PyFlow {
        self.role.flow()
    }

    #[getter]
    fn data(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.data
            .as_ref()
            .map(|data| data.clone_ref(py))
            .ok_or_else(|| PyReferenceError::new_err("half-edge specification has been cleared"))
    }

    #[setter]
    fn set_data(&mut self, value: Py<PyAny>) {
        self.data = Some(value);
    }

    #[getter]
    fn drawing(&self, py: Python<'_>) -> PyResult<PyHalfEdgeDrawing> {
        self.drawing
            .as_ref()
            .map(|drawing| PyHalfEdgeDrawing::detached(py, drawing))
            .ok_or_else(|| PyReferenceError::new_err("half-edge specification has been cleared"))
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.node)?;
        visit.call(&self.data)?;
        visit.call(&self.drawing)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.node = None;
        self.data = None;
        self.drawing = None;
    }
}

/// A declarative internal or external edge description consumed by `build()`.
#[gen_stub_pyclass]
#[pyclass(unsendable, name = "EdgeSpec")]
pub struct PyEdgeSpec {
    name: Option<String>,
    first: Option<Py<PyHalfEdgeSpec>>,
    second: Option<Py<PyHalfEdgeSpec>>,
    data: Option<Py<PyAny>>,
    drawing: Option<Py<PyDict>>,
    orientation: PyOrientation,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyEdgeSpec {
    #[getter]
    fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    #[getter]
    fn data(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.data
            .as_ref()
            .map(|data| data.clone_ref(py))
            .ok_or_else(|| PyReferenceError::new_err("edge specification has been cleared"))
    }

    #[setter]
    fn set_data(&mut self, value: Py<PyAny>) {
        self.data = Some(value);
    }

    #[getter]
    fn drawing(&self, py: Python<'_>) -> PyResult<PyEdgeDrawing> {
        self.drawing
            .as_ref()
            .map(|drawing| PyEdgeDrawing::detached(py, drawing))
            .ok_or_else(|| PyReferenceError::new_err("edge specification has been cleared"))
    }

    #[getter]
    fn orientation(&self) -> PyOrientation {
        self.orientation
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.first)?;
        visit.call(&self.second)?;
        visit.call(&self.data)?;
        visit.call(&self.drawing)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.first = None;
        self.second = None;
        self.data = None;
        self.drawing = None;
    }
}

/// Describe a node while preserving its arbitrary Python data by identity.
#[gen_stub_pyfunction(python = r#"
    import typing

    def node(name: _OptionalString = None, *, data: typing.Any = None, label: _OptionalStaticContent = ..., placement: _PlacementValue = ..., shift: _DrawingPoint = ..., rank: _OptionalInteger = ..., minimum_size: _OptionalNumber = ..., maximum_size: _OptionalNumber = ..., style: _OptionalStyle = ..., label_style: _OptionalStyle = ..., extensions: _NativeDict = ...) -> NodeSpec:
        """Describe a node while preserving its arbitrary Python data by identity."""
        ...
"#)]
#[pyfunction(signature = (name=None, *, data=None, **drawing))]
pub fn node(
    py: Python<'_>,
    name: Option<String>,
    data: Option<Py<PyAny>>,
    drawing: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyNodeSpec> {
    Ok(PyNodeSpec {
        token: SPEC_ID.fetch_add(1, Ordering::Relaxed),
        name,
        data: Some(data.unwrap_or_else(|| py.None())),
        drawing: Some(drawing_dict(py, NODE_FIELDS, drawing)?),
    })
}

fn endpoint(
    py: Python<'_>,
    node: Py<PyNodeSpec>,
    role: EndpointRole,
    data: Option<Py<PyAny>>,
    drawing: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyHalfEdgeSpec> {
    Ok(PyHalfEdgeSpec {
        node: Some(node),
        role,
        data: Some(data.unwrap_or_else(|| py.None())),
        drawing: Some(drawing_dict(py, HEDGE_FIELDS, drawing)?),
    })
}

/// Attach a source endpoint and its data and drawing metadata to a node spec.
#[gen_stub_pyfunction(python = r#"
    import typing

    def source(node: NodeSpec, *, data: typing.Any = None, label: _OptionalStaticContent = ..., statement: _DrawingString = ..., port_label: _DrawingString = ..., compass: _CompassValue = ..., anchor: _AnchorValue = ..., routing: _RoutingValue = ..., style: _OptionalStyleLayers = ..., extensions: _NativeDict = ...) -> HalfEdgeSpec:
        """Attach a source endpoint and its metadata to a node spec."""
        ...
"#)]
#[pyfunction(signature = (node, *, data=None, **drawing))]
pub fn source(
    py: Python<'_>,
    node: Py<PyNodeSpec>,
    data: Option<Py<PyAny>>,
    drawing: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyHalfEdgeSpec> {
    endpoint(py, node, EndpointRole::Source, data, drawing)
}

/// Attach a sink endpoint and its data and drawing metadata to a node spec.
#[gen_stub_pyfunction(python = r#"
    import typing

    def sink(node: NodeSpec, *, data: typing.Any = None, label: _OptionalStaticContent = ..., statement: _DrawingString = ..., port_label: _DrawingString = ..., compass: _CompassValue = ..., anchor: _AnchorValue = ..., routing: _RoutingValue = ..., style: _OptionalStyleLayers = ..., extensions: _NativeDict = ...) -> HalfEdgeSpec:
        """Attach a sink endpoint and its metadata to a node spec."""
        ...
"#)]
#[pyfunction(signature = (node, *, data=None, **drawing))]
pub fn sink(
    py: Python<'_>,
    node: Py<PyNodeSpec>,
    data: Option<Py<PyAny>>,
    drawing: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyHalfEdgeSpec> {
    endpoint(py, node, EndpointRole::Sink, data, drawing)
}

/// Describe an edge from one or two endpoint specs.
#[gen_stub_pyfunction(python = r#"
    import typing

    def edge(first: HalfEdgeSpec, name: _OptionalString = None, second: _OptionalHalfEdgeSpec = None, *, data: typing.Any = None, orientation: Orientation = Orientation.Default, label: _OptionalStaticContent = ..., particle: _DrawingString = ..., momentum: _DrawingScalar = ..., cut_id: _DrawingScalar = ..., placement: _PlacementValue = ..., label_position: _DrawingPoint = ..., label_offset: _OptionalNumber = ..., label_angle: _DrawingAngle = ..., bend: _DrawingAngle = ..., routing: _RoutingValue = ..., minimum_length: _OptionalInteger = ..., same_rank: _OptionalBoolean = ..., style: _OptionalStyleLayers = ..., label_style: _OptionalStyle = ..., momentum_style: _OptionalStyleLayers = ..., decoration: _DrawingDecoration = ..., extensions: _NativeDict = ...) -> EdgeSpec:
        """Describe an edge from one or two endpoint specs."""
        ...
"#)]
#[pyfunction(signature = (first, name=None, second=None, *, data=None, orientation=PyOrientation::Default, **drawing))]
pub fn edge(
    py: Python<'_>,
    first: Py<PyHalfEdgeSpec>,
    name: Option<String>,
    second: Option<Py<PyHalfEdgeSpec>>,
    data: Option<Py<PyAny>>,
    orientation: PyOrientation,
    drawing: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyEdgeSpec> {
    if let Some(second) = &second {
        let first_role = first.borrow(py).role;
        let second_role = second.borrow(py).role;
        if first_role == second_role {
            return Err(PyValueError::new_err(
                "an internal edge needs one source and one sink",
            ));
        }
    }
    Ok(PyEdgeSpec {
        name,
        first: Some(first),
        second,
        data: Some(data.unwrap_or_else(|| py.None())),
        drawing: Some(drawing_dict(py, EDGE_FIELDS, drawing)?),
        orientation,
    })
}

pub(crate) struct NodeRecord {
    pub(crate) name: Option<String>,
    pub(crate) data: Py<PyAny>,
    pub(crate) drawing: Py<PyDict>,
}

impl NodeRecord {
    pub(crate) fn copy(&self, py: Python<'_>) -> PyResult<Self> {
        Ok(Self {
            name: self.name.clone(),
            data: self.data.clone_ref(py),
            drawing: copy_drawing(py, &self.drawing)?,
        })
    }
}

pub(crate) struct HalfEdgeRecord {
    pub(crate) data: Py<PyAny>,
    pub(crate) drawing: Py<PyDict>,
}

impl HalfEdgeRecord {
    pub(crate) fn copy(&self, py: Python<'_>) -> PyResult<Self> {
        Ok(Self {
            data: self.data.clone_ref(py),
            drawing: copy_drawing(py, &self.drawing)?,
        })
    }
}

pub(crate) struct EdgeRecord {
    pub(crate) name: Option<String>,
    pub(crate) data: Py<PyAny>,
    pub(crate) drawing: Py<PyDict>,
}

impl EdgeRecord {
    pub(crate) fn copy(&self, py: Python<'_>) -> PyResult<Self> {
        Ok(Self {
            name: self.name.clone(),
            data: self.data.clone_ref(py),
            drawing: copy_drawing(py, &self.drawing)?,
        })
    }
}

pub(crate) struct GraphState {
    pub(crate) graph: PyHedgeGraph,
    pub(crate) name: Option<String>,
    pub(crate) global_data: PyGlobalData,
    pub(crate) codec: Option<Py<PyDotCodec>>,
    pub(crate) render_config: Py<PyAny>,
    pub(crate) revision: u64,
}

impl GraphState {
    pub(crate) fn copy_graph(&self, py: Python<'_>) -> PyResult<PyHedgeGraph> {
        self.graph.copy(py)
    }

    pub(crate) fn copy(&self, py: Python<'_>) -> PyResult<Self> {
        let revision = self.revision;
        let mut state = self.derived(py, self.copy_graph(py)?)?;
        state.revision = revision;
        Ok(state)
    }

    pub(crate) fn derived(&self, py: Python<'_>, graph: PyHedgeGraph) -> PyResult<Self> {
        Ok(Self {
            graph,
            name: self.name.clone(),
            global_data: self.global_data.clone(),
            codec: self.codec.as_ref().map(|codec| codec.clone_ref(py)),
            render_config: crate::typst::render_config_copy(py, &self.render_config)?,
            revision: 0,
        })
    }
}

/// An owned topology with arbitrary element data and typed rendering configuration.
#[gen_stub_pyclass]
#[pyclass(unsendable, name = "Graph")]
pub struct PyGraph {
    pub(crate) state: RefCell<Option<GraphState>>,
}

impl PyGraph {
    pub(crate) fn from_state(state: GraphState) -> Self {
        Self {
            state: RefCell::new(Some(state)),
        }
    }

    pub(crate) fn check_revision(&self, revision: u64) -> PyResult<()> {
        let state = self.state.borrow();
        let state = state
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        if state.revision != revision {
            return Err(PyReferenceError::new_err(
                "this graph view is stale after a topology mutation",
            ));
        }
        Ok(())
    }

    pub(crate) fn revision(&self) -> PyResult<u64> {
        self.state
            .borrow()
            .as_ref()
            .map(|state| state.revision)
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))
    }

    pub(crate) fn resolve_node(&self, key: &Bound<'_, PyAny>) -> PyResult<usize> {
        let state = self.state.borrow();
        let state = state
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        if let Ok(index) = key.extract::<usize>() {
            return (index < state.graph.n_nodes())
                .then_some(index)
                .ok_or_else(|| PyIndexError::new_err("node index out of range"));
        }
        if let Ok(name) = key.extract::<String>() {
            return state
                .graph
                .iter_nodes()
                .find_map(|(index, _, node)| {
                    (node.name.as_deref() == Some(name.as_str())).then_some(index.0)
                })
                .ok_or_else(|| PyKeyError::new_err(format!("unknown node {name:?}")));
        }
        Err(PyTypeError::new_err("node key must be an integer or name"))
    }

    pub(crate) fn resolve_edge(&self, key: &Bound<'_, PyAny>) -> PyResult<usize> {
        let state = self.state.borrow();
        let state = state
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        if let Ok(index) = key.extract::<usize>() {
            return (index < state.graph.n_edges())
                .then_some(index)
                .ok_or_else(|| PyIndexError::new_err("edge index out of range"));
        }
        if let Ok(name) = key.extract::<String>() {
            return state
                .graph
                .iter_edges()
                .find_map(|(_, index, edge)| {
                    (edge.data.name.as_deref() == Some(name.as_str())).then_some(index.0)
                })
                .ok_or_else(|| PyKeyError::new_err(format!("unknown edge {name:?}")));
        }
        Err(PyTypeError::new_err("edge key must be an integer or name"))
    }
}

#[derive(Clone, Copy)]
enum ViewIndex {
    Node(usize),
    Edge(usize),
    Hedge(usize),
}

macro_rules! graph_view {
    ($doc:literal, $rust:ident, $python:literal, $index_variant:ident) => {
        #[doc = $doc]
        #[gen_stub_pyclass]
        #[pyclass(unsendable, name = $python)]
        pub struct $rust {
            graph: Option<Py<PyGraph>>,
            index: usize,
            revision: u64,
        }

        impl $rust {
            pub(crate) fn new(graph: Py<PyGraph>, index: usize, revision: u64) -> Self {
                Self {
                    graph: Some(graph),
                    index,
                    revision,
                }
            }

            fn owner<'py>(&self, py: Python<'py>) -> PyResult<&Bound<'py, PyGraph>> {
                let graph = self
                    .graph
                    .as_ref()
                    .ok_or_else(|| PyReferenceError::new_err("graph view has been cleared"))?;
                graph.borrow(py).check_revision(self.revision)?;
                Ok(graph.bind(py))
            }

            fn clone_graph(&self, py: Python<'_>) -> PyResult<Py<PyGraph>> {
                self.owner(py)?;
                Ok(self.graph.as_ref().expect("checked").clone_ref(py))
            }
        }
    };
}

graph_view!(
    "A live node view whose data and drawing metadata update its graph.",
    PyNode,
    "Node",
    Node
);
graph_view!(
    "A live edge view whose data and drawing metadata update its graph.",
    PyEdge,
    "Edge",
    Edge
);
graph_view!(
    "A live half-edge view whose data and drawing metadata update its graph.",
    PyHalfEdge,
    "HalfEdge",
    Hedge
);

impl PyEdge {
    fn endpoint_view(
        &self,
        py: Python<'_>,
        role: EndpointRole,
    ) -> PyResult<Option<Py<PyHalfEdge>>> {
        let graph_obj = self.clone_graph(py)?;
        let graph = graph_obj.borrow(py);
        let state = graph.state.borrow();
        let state = state.as_ref().expect("checked");
        let pair = state.graph[&EdgeIndex(self.index)].1;
        let hedge = match (pair, role) {
            (HedgePair::Paired { source, .. }, EndpointRole::Source)
            | (HedgePair::Split { source, .. }, EndpointRole::Source) => Some(source),
            (HedgePair::Paired { sink, .. }, EndpointRole::Sink)
            | (HedgePair::Split { sink, .. }, EndpointRole::Sink) => Some(sink),
            (
                HedgePair::Unpaired {
                    hedge,
                    flow: Flow::Source,
                },
                EndpointRole::Source,
            ) => Some(hedge),
            (
                HedgePair::Unpaired {
                    hedge,
                    flow: Flow::Sink,
                },
                EndpointRole::Sink,
            ) => Some(hedge),
            _ => None,
        };
        hedge
            .map(|index| {
                Py::new(
                    py,
                    PyHalfEdge::new(graph_obj.clone_ref(py), index.0, self.revision),
                )
            })
            .transpose()
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl PyNode {
    #[getter]
    fn index(&self, py: Python<'_>) -> PyResult<usize> {
        self.owner(py)?;
        Ok(self.index)
    }

    #[getter]
    fn name(&self, py: Python<'_>) -> PyResult<Option<String>> {
        let graph = self.owner(py)?.borrow();
        let state = graph.state.borrow();
        Ok(
            state.as_ref().expect("checked").graph[NodeIndex(self.index)]
                .name
                .clone(),
        )
    }

    #[getter]
    fn data(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let graph = self.owner(py)?.borrow();
        let state = graph.state.borrow();
        Ok(
            state.as_ref().expect("checked").graph[NodeIndex(self.index)]
                .data
                .clone_ref(py),
        )
    }

    #[setter]
    fn set_data(slf: &Bound<'_, Self>, value: Py<PyAny>) -> PyResult<()> {
        let py = slf.py();
        let view = slf.borrow();
        let graph = view.owner(py)?.borrow();
        graph.state.borrow_mut().as_mut().expect("checked").graph[NodeIndex(view.index)].data =
            value;
        Ok(())
    }

    #[getter]
    fn drawing(&self, py: Python<'_>) -> PyResult<PyNodeDrawing> {
        let graph_obj = self.clone_graph(py)?;
        let graph = graph_obj.borrow(py);
        let state = graph.state.borrow();
        Ok(PyNodeDrawing::live(
            py,
            &state.as_ref().expect("checked").graph[NodeIndex(self.index)].drawing,
            graph_obj.clone_ref(py),
            self.revision,
        ))
    }

    #[getter]
    fn incidence(&self, py: Python<'_>) -> PyResult<Vec<Py<PyHalfEdge>>> {
        let graph_obj = self.clone_graph(py)?;
        let graph = graph_obj.borrow(py);
        let state = graph.state.borrow();
        let state = state.as_ref().expect("checked");
        state
            .graph
            .iter_crown(NodeIndex(self.index))
            .map(|hedge| {
                Py::new(
                    py,
                    PyHalfEdge::new(graph_obj.clone_ref(py), hedge.0, self.revision),
                )
            })
            .collect()
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(match self.name(py)? {
            Some(name) => format!("Node({}: {name:?})", self.index),
            None => format!("Node({})", self.index),
        })
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.graph)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.graph = None;
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl PyEdge {
    #[getter]
    fn index(&self, py: Python<'_>) -> PyResult<usize> {
        self.owner(py)?;
        Ok(self.index)
    }

    #[getter]
    fn name(&self, py: Python<'_>) -> PyResult<Option<String>> {
        let graph = self.owner(py)?.borrow();
        let state = graph.state.borrow();
        Ok(
            state.as_ref().expect("checked").graph[EdgeIndex(self.index)]
                .name
                .clone(),
        )
    }

    #[getter]
    fn data(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let graph = self.owner(py)?.borrow();
        let state = graph.state.borrow();
        Ok(
            state.as_ref().expect("checked").graph[EdgeIndex(self.index)]
                .data
                .clone_ref(py),
        )
    }

    #[setter]
    fn set_data(slf: &Bound<'_, Self>, value: Py<PyAny>) -> PyResult<()> {
        let py = slf.py();
        let view = slf.borrow();
        let graph = view.owner(py)?.borrow();
        graph.state.borrow_mut().as_mut().expect("checked").graph[EdgeIndex(view.index)].data =
            value;
        Ok(())
    }

    #[getter]
    fn drawing(&self, py: Python<'_>) -> PyResult<PyEdgeDrawing> {
        let graph_obj = self.clone_graph(py)?;
        let graph = graph_obj.borrow(py);
        let state = graph.state.borrow();
        Ok(PyEdgeDrawing::live(
            py,
            &state.as_ref().expect("checked").graph[EdgeIndex(self.index)].drawing,
            graph_obj.clone_ref(py),
            self.revision,
        ))
    }

    #[getter]
    fn orientation(&self, py: Python<'_>) -> PyResult<PyOrientation> {
        let graph = self.owner(py)?.borrow();
        let state = graph.state.borrow();
        Ok(state
            .as_ref()
            .expect("checked")
            .graph
            .orientation(EdgeIndex(self.index))
            .into())
    }

    #[getter]
    fn source(&self, py: Python<'_>) -> PyResult<Option<Py<PyHalfEdge>>> {
        self.endpoint_view(py, EndpointRole::Source)
    }

    #[getter]
    fn sink(&self, py: Python<'_>) -> PyResult<Option<Py<PyHalfEdge>>> {
        self.endpoint_view(py, EndpointRole::Sink)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(match self.name(py)? {
            Some(name) => format!("Edge({}: {name:?})", self.index),
            None => format!("Edge({})", self.index),
        })
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.graph)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.graph = None;
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl PyHalfEdge {
    #[getter]
    fn index(&self, py: Python<'_>) -> PyResult<usize> {
        self.owner(py)?;
        Ok(self.index)
    }

    #[getter]
    fn flow(&self, py: Python<'_>) -> PyResult<PyFlow> {
        let graph = self.owner(py)?.borrow();
        let state = graph.state.borrow();
        Ok(state
            .as_ref()
            .expect("checked")
            .graph
            .flow(Hedge(self.index))
            .into())
    }

    #[getter]
    fn data(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let graph = self.owner(py)?.borrow();
        let state = graph.state.borrow();
        Ok(state.as_ref().expect("checked").graph[Hedge(self.index)]
            .data
            .clone_ref(py))
    }

    #[setter]
    fn set_data(slf: &Bound<'_, Self>, value: Py<PyAny>) -> PyResult<()> {
        let py = slf.py();
        let view = slf.borrow();
        let graph = view.owner(py)?.borrow();
        graph.state.borrow_mut().as_mut().expect("checked").graph[Hedge(view.index)].data = value;
        Ok(())
    }

    #[getter]
    fn drawing(&self, py: Python<'_>) -> PyResult<PyHalfEdgeDrawing> {
        let graph_obj = self.clone_graph(py)?;
        let graph = graph_obj.borrow(py);
        let state = graph.state.borrow();
        Ok(PyHalfEdgeDrawing::live(
            py,
            &state.as_ref().expect("checked").graph[Hedge(self.index)].drawing,
            graph_obj.clone_ref(py),
            self.revision,
        ))
    }

    #[getter]
    fn node(&self, py: Python<'_>) -> PyResult<Py<PyNode>> {
        let graph_obj = self.clone_graph(py)?;
        let graph = graph_obj.borrow(py);
        let state = graph.state.borrow();
        let node = state
            .as_ref()
            .expect("checked")
            .graph
            .node_id(Hedge(self.index));
        Py::new(
            py,
            PyNode::new(graph_obj.clone_ref(py), node.0, self.revision),
        )
    }

    #[getter]
    fn edge(&self, py: Python<'_>) -> PyResult<Py<PyEdge>> {
        let graph_obj = self.clone_graph(py)?;
        let graph = graph_obj.borrow(py);
        let state = graph.state.borrow();
        let edge = state.as_ref().expect("checked").graph[&Hedge(self.index)];
        Py::new(
            py,
            PyEdge::new(graph_obj.clone_ref(py), edge.0, self.revision),
        )
    }

    #[getter]
    fn pair(&self, py: Python<'_>) -> PyResult<Option<Py<PyHalfEdge>>> {
        let graph_obj = self.clone_graph(py)?;
        let graph = graph_obj.borrow(py);
        let state = graph.state.borrow();
        let hedge = Hedge(self.index);
        let pair = state.as_ref().expect("checked").graph.inv(hedge);
        (pair != hedge)
            .then_some(pair)
            .map(|pair| {
                Py::new(
                    py,
                    PyHalfEdge::new(graph_obj.clone_ref(py), pair.0, self.revision),
                )
            })
            .transpose()
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.owner(py)?;
        Ok(format!("HalfEdge({})", self.index))
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.graph)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.graph = None;
    }
}

fn checked_order(order: Vec<usize>, len: usize, kind: &str) -> PyResult<Vec<usize>> {
    if order.len() != len {
        return Err(PyValueError::new_err(format!(
            "{kind} order has length {}, expected {len}",
            order.len()
        )));
    }
    let values = order.iter().copied().collect::<BTreeSet<_>>();
    if values.len() != len || values.iter().copied().ne(0..len) {
        return Err(PyValueError::new_err(format!(
            "{kind} order must be a permutation of 0..{len}"
        )));
    }
    Ok(order)
}

fn view_for(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    revision: u64,
    index: ViewIndex,
) -> PyResult<Py<PyAny>> {
    Ok(match index {
        ViewIndex::Node(index) => {
            Py::new(py, PyNode::new(graph.clone_ref(py), index, revision))?.into_any()
        }
        ViewIndex::Edge(index) => {
            Py::new(py, PyEdge::new(graph.clone_ref(py), index, revision))?.into_any()
        }
        ViewIndex::Hedge(index) => {
            Py::new(py, PyHalfEdge::new(graph.clone_ref(py), index, revision))?.into_any()
        }
    })
}

#[gen_stub_pymethods]
#[pymethods]
impl PyGraph {
    #[new]
    #[pyo3(signature = (*, name=None, global_data=None, codec=None, render_config=None, node_store=PyNodeStore::Vec))]
    fn new(
        py: Python<'_>,
        name: Option<String>,
        global_data: Option<PyGlobalData>,
        codec: Option<Py<PyDotCodec>>,
        #[gen_stub(override_type(type_repr = "RenderConfig | None"))] render_config: Option<
            Py<PyAny>,
        >,
        node_store: PyNodeStore,
    ) -> PyResult<Self> {
        Ok(Self::from_state(GraphState {
            graph: PyHedgeGraph::empty(node_store),
            name,
            global_data: global_data.unwrap_or_default(),
            codec,
            render_config: match render_config {
                Some(config) => {
                    crate::typst::validate_render_config(py, &config)?;
                    config
                }
                None => crate::typst::default_render_config(py)?,
            },
            revision: 0,
        }))
    }

    #[getter]
    fn name(&self) -> PyResult<Option<String>> {
        Ok(self
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .name
            .clone())
    }

    #[getter]
    fn global_data(&self) -> PyResult<PyGlobalData> {
        Ok(self
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .global_data
            .clone())
    }

    #[setter]
    fn set_global_data(&self, value: PyGlobalData) -> PyResult<()> {
        self.state
            .borrow_mut()
            .as_mut()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .global_data = value;
        Ok(())
    }

    #[getter]
    #[gen_stub(override_return_type(type_repr = "RenderConfig"))]
    fn render_config(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(self
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .render_config
            .clone_ref(py))
    }

    #[setter]
    fn set_render_config(
        slf: &Bound<'_, Self>,
        #[gen_stub(override_type(type_repr = "RenderConfig"))] value: Py<PyAny>,
    ) -> PyResult<()> {
        let py = slf.py();
        crate::typst::validate_render_config(py, &value)?;
        slf.borrow()
            .state
            .borrow_mut()
            .as_mut()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .render_config = value;
        Ok(())
    }

    fn __len__(&self) -> PyResult<usize> {
        Ok(self
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .graph
            .n_nodes())
    }

    #[getter]
    fn n_nodes(&self) -> PyResult<usize> {
        self.__len__()
    }

    #[getter]
    fn n_edges(&self) -> PyResult<usize> {
        Ok(self
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .graph
            .n_edges())
    }

    #[getter]
    fn n_half_edges(&self) -> PyResult<usize> {
        Ok(self
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .graph
            .n_hedges())
    }

    /// Return the node-storage strategy used by this graph.
    #[getter]
    fn node_store(&self) -> PyResult<PyNodeStore> {
        Ok(self
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .graph
            .node_store())
    }

    /// Copy this graph into another node-storage strategy.
    #[pyo3(signature = (node_store))]
    fn to_node_store(&self, py: Python<'_>, node_store: PyNodeStore) -> PyResult<Self> {
        let state = self.state.borrow();
        let state = state
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        let graph = state
            .copy_graph(py)?
            .into_node_store(node_store)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok(Self::from_state(state.derived(py, graph)?))
    }

    #[pyo3(signature = (key))]
    fn node(
        slf: Py<PyGraph>,
        #[gen_stub(override_type(type_repr="builtins.int | builtins.str", imports=("builtins")))]
        key: &Bound<'_, PyAny>,
        py: Python<'_>,
    ) -> PyResult<Py<PyNode>> {
        let index = slf.borrow(py).resolve_node(key)?;
        let revision = slf.borrow(py).revision()?;
        Py::new(py, PyNode::new(slf.clone_ref(py), index, revision))
    }

    #[pyo3(signature = (key))]
    fn edge(
        slf: Py<PyGraph>,
        #[gen_stub(override_type(type_repr="builtins.int | builtins.str", imports=("builtins")))]
        key: &Bound<'_, PyAny>,
        py: Python<'_>,
    ) -> PyResult<Py<PyEdge>> {
        let index = slf.borrow(py).resolve_edge(key)?;
        let revision = slf.borrow(py).revision()?;
        Py::new(py, PyEdge::new(slf.clone_ref(py), index, revision))
    }

    #[pyo3(signature = (index))]
    fn half_edge(slf: Py<PyGraph>, index: usize, py: Python<'_>) -> PyResult<Py<PyHalfEdge>> {
        let graph = slf.borrow(py);
        let revision = graph.revision()?;
        let len = graph
            .state
            .borrow()
            .as_ref()
            .expect("checked")
            .graph
            .n_hedges();
        if index >= len {
            return Err(PyIndexError::new_err("half-edge index out of range"));
        }
        drop(graph);
        Py::new(py, PyHalfEdge::new(slf.clone_ref(py), index, revision))
    }

    #[pyo3(signature = ())]
    fn nodes(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<Vec<Py<PyNode>>> {
        let graph = slf.borrow(py);
        let revision = graph.revision()?;
        let len = graph
            .state
            .borrow()
            .as_ref()
            .expect("checked")
            .graph
            .n_nodes();
        drop(graph);
        (0..len)
            .map(|index| Py::new(py, PyNode::new(slf.clone_ref(py), index, revision)))
            .collect()
    }

    #[pyo3(signature = ())]
    fn edges(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<Vec<Py<PyEdge>>> {
        let graph = slf.borrow(py);
        let revision = graph.revision()?;
        let len = graph
            .state
            .borrow()
            .as_ref()
            .expect("checked")
            .graph
            .n_edges();
        drop(graph);
        (0..len)
            .map(|index| Py::new(py, PyEdge::new(slf.clone_ref(py), index, revision)))
            .collect()
    }

    #[pyo3(signature = ())]
    fn half_edges(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<Vec<Py<PyHalfEdge>>> {
        let graph = slf.borrow(py);
        let revision = graph.revision()?;
        let len = graph
            .state
            .borrow()
            .as_ref()
            .expect("checked")
            .graph
            .n_hedges();
        drop(graph);
        (0..len)
            .map(|index| Py::new(py, PyHalfEdge::new(slf.clone_ref(py), index, revision)))
            .collect()
    }

    fn reorder_nodes(&self, order: Vec<usize>) -> PyResult<()> {
        let mut state = self.state.borrow_mut();
        let state = state
            .as_mut()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        let order = checked_order(order, state.graph.n_nodes(), "node")?;
        <PyHedgeGraph as Swap<NodeIndex>>::permute(&mut state.graph, &Permutation::from_inv(order));
        state.revision += 1;
        Ok(())
    }

    fn reorder_edges(&self, order: Vec<usize>) -> PyResult<()> {
        let mut state = self.state.borrow_mut();
        let state = state
            .as_mut()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        let order = checked_order(order, state.graph.n_edges(), "edge")?;
        let hedge_order = order
            .iter()
            .flat_map(|index| match state.graph[&EdgeIndex(*index)].1 {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    [Some(source.0), Some(sink.0)]
                }
                HedgePair::Unpaired { hedge, .. } => [Some(hedge.0), None],
            })
            .flatten()
            .collect();
        <PyHedgeGraph as Swap<Hedge>>::permute(
            &mut state.graph,
            &Permutation::from_inv(hedge_order),
        );
        <PyHedgeGraph as Swap<EdgeIndex>>::permute(&mut state.graph, &Permutation::from_inv(order));
        state.revision += 1;
        Ok(())
    }

    fn reverse_edge(
        &self,
        #[gen_stub(override_type(type_repr="builtins.int | builtins.str", imports=("builtins")))]
        key: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let index = self.resolve_edge(key)?;
        let mut state = self.state.borrow_mut();
        let state = state.as_mut().expect("checked");
        let edge = EdgeIndex(index);
        let pair = state.graph[&edge].1;
        let orientation = state.graph.orientation(edge);
        match pair {
            HedgePair::Paired { source, .. } | HedgePair::Split { source, .. } => {
                state.graph.set_flow(source, Flow::Sink)
            }
            HedgePair::Unpaired { hedge, flow } => state.graph.set_flow(hedge, -flow),
        }
        state.graph.set_orientation(edge, orientation);
        state.revision += 1;
        Ok(())
    }

    #[pyo3(signature = (*, node=None, edge=None, source=None, sink=None))]
    fn map(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Callable[[Node], typing.Any] | None", imports=("typing")))]
        node: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[Edge], typing.Any] | None", imports=("typing")))]
        edge: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdge], typing.Any] | None", imports=("typing")))]
        source: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdge], typing.Any] | None", imports=("typing")))]
        sink: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        let (revision, mut mapped) = {
            let graph = slf.borrow(py);
            let state_ref = graph.state.borrow();
            let state = state_ref
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
            (state.revision, state.copy(py)?)
        };

        for index in 0..mapped.graph.n_nodes() {
            if let Some(callback) = &node {
                mapped.graph[NodeIndex(index)].data =
                    callback.call1(py, (view_for(py, &slf, revision, ViewIndex::Node(index))?,))?;
            }
        }
        let source_callback = source.as_ref();
        let sink_callback = sink.as_ref();
        for index in 0..mapped.graph.n_edges() {
            if let Some(callback) = &edge {
                mapped.graph[EdgeIndex(index)].data =
                    callback.call1(py, (view_for(py, &slf, revision, ViewIndex::Edge(index))?,))?;
            }
            let pair = mapped.graph[&EdgeIndex(index)].1;
            for (hedge, callback) in match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    [(Some(source), source_callback), (Some(sink), sink_callback)]
                }
                HedgePair::Unpaired { hedge, flow } => match flow {
                    Flow::Source => [(Some(hedge), source_callback), (None, None)],
                    Flow::Sink => [(None, None), (Some(hedge), sink_callback)],
                },
            } {
                if let (Some(hedge), Some(callback)) = (hedge, callback) {
                    mapped.graph[hedge].data = callback.call1(
                        py,
                        (view_for(py, &slf, revision, ViewIndex::Hedge(hedge.0))?,),
                    )?;
                }
            }
        }
        slf.borrow(py).check_revision(revision)?;
        mapped.revision = 0;
        Ok(Self::from_state(mapped))
    }

    #[classmethod]
    #[pyo3(signature = (source, codec, *, node_store=PyNodeStore::Vec))]
    fn from_dot(
        _cls: &Bound<'_, PyType>,
        py: Python<'_>,
        source: &str,
        codec: Py<PyDotCodec>,
        node_store: PyNodeStore,
    ) -> PyResult<Self> {
        dot::decode_dot(py, source, codec, node_store)
    }

    #[classmethod]
    #[pyo3(signature = (source, codec, *, node_store=PyNodeStore::Vec))]
    fn from_dot_set(
        _cls: &Bound<'_, PyType>,
        py: Python<'_>,
        source: &str,
        codec: Py<PyDotCodec>,
        node_store: PyNodeStore,
    ) -> PyResult<Vec<Self>> {
        dot::decode_dot_set(py, source, codec, node_store)
    }

    #[classmethod]
    #[pyo3(signature = (path, codec, *, node_store=PyNodeStore::Vec))]
    fn from_dot_file(
        _cls: &Bound<'_, PyType>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="builtins.str | os.PathLike[builtins.str]", imports=("builtins", "os")))]
        path: PathBuf,
        codec: Py<PyDotCodec>,
        node_store: PyNodeStore,
    ) -> PyResult<Self> {
        let source = std::fs::read_to_string(&path).map_err(|error| {
            PyValueError::new_err(format!("failed to read {}: {error}", path.display()))
        })?;
        dot::decode_dot(py, &source, codec, node_store)
    }

    #[pyo3(signature = (codec=None))]
    fn to_dot(&self, py: Python<'_>, codec: Option<Py<PyDotCodec>>) -> PyResult<String> {
        dot::encode_graph(py, self, codec)
    }

    #[pyo3(signature = (output, *, config=None))]
    fn render(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="builtins.str | os.PathLike[builtins.str]", imports=("builtins", "os")))]
        output: PathBuf,
        #[gen_stub(override_type(type_repr = "RenderConfig | None"))] config: Option<
            &Bound<'_, PyAny>,
        >,
    ) -> PyResult<PathBuf> {
        crate::render::render_graph(py, &slf, output, config)
    }

    #[pyo3(signature = (*, config=None))]
    fn to_svg(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "RenderConfig | None"))] config: Option<
            &Bound<'_, PyAny>,
        >,
    ) -> PyResult<String> {
        crate::render::graph_to_svg(py, &slf, config)
    }

    #[pyo3(signature = ())]
    fn _repr_svg_(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<String> {
        crate::render::graph_to_svg(py, &slf, None)
    }

    fn __repr__(&self) -> PyResult<String> {
        let state = self.state.borrow();
        let state = state
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        Ok(format!(
            "Graph(nodes={}, edges={}, half_edges={}, node_store={:?})",
            state.graph.n_nodes(),
            state.graph.n_edges(),
            state.graph.n_hedges(),
            state.graph.node_store(),
        ))
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        let Ok(state) = self.state.try_borrow() else {
            return Ok(());
        };
        let Some(state) = state.as_ref() else {
            return Ok(());
        };
        for (_, _, node) in state.graph.iter_nodes() {
            visit.call(&node.data)?;
            visit.call(&node.drawing)?;
        }
        for (_, _, edge) in state.graph.iter_edges() {
            visit.call(&edge.data.data)?;
            visit.call(&edge.data.drawing)?;
        }
        for (_, half_edge) in state.graph.iter_hedges() {
            visit.call(&half_edge.data)?;
            visit.call(&half_edge.drawing)?;
        }
        visit.call(&state.codec)?;
        visit.call(&state.render_config)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        *self.state.get_mut() = None;
    }
}

fn endpoint_record(
    py: Python<'_>,
    endpoint: &Py<PyHalfEdgeSpec>,
    node_indices: &BTreeMap<u64, NodeIndex>,
) -> PyResult<(EndpointRole, HedgeData<HalfEdgeRecord>)> {
    let endpoint = endpoint.borrow(py);
    let node = endpoint
        .node
        .as_ref()
        .ok_or_else(|| PyReferenceError::new_err("half-edge specification has been cleared"))?
        .borrow(py);
    let index = node_indices.get(&node.token).copied().ok_or_else(|| {
        PyValueError::new_err("an edge endpoint refers to a node that was not passed to build()")
    })?;
    Ok((
        endpoint.role,
        HedgeData {
            node: index,
            is_in_subgraph: false,
            data: HalfEdgeRecord {
                data: endpoint
                    .data
                    .as_ref()
                    .ok_or_else(|| {
                        PyReferenceError::new_err("half-edge specification has been cleared")
                    })?
                    .clone_ref(py),
                drawing: copy_drawing(
                    py,
                    endpoint.drawing.as_ref().ok_or_else(|| {
                        PyReferenceError::new_err("half-edge specification has been cleared")
                    })?,
                )?,
            },
        },
    ))
}

/// Build a graph from declarative node and edge specs.
#[gen_stub_pyfunction(python = r#"
    import typing

    def build(*items: _GraphItem, name: _OptionalString = None, global_data: _OptionalGlobalData = None, codec: _OptionalDotCodec = None, render_config: _OptionalRenderConfig = None, node_store: NodeStore = NodeStore.Vec) -> Graph:
        """Build a graph from declarative node and edge specs."""
        ...
"#)]
#[pyfunction(signature = (*items, name=None, global_data=None, codec=None, render_config=None, node_store=PyNodeStore::Vec))]
pub fn build(
    py: Python<'_>,
    items: &Bound<'_, PyTuple>,
    name: Option<String>,
    global_data: Option<PyGlobalData>,
    codec: Option<Py<PyDotCodec>>,
    render_config: Option<Py<PyAny>>,
    node_store: PyNodeStore,
) -> PyResult<PyGraph> {
    let mut node_specs = Vec::<Py<PyNodeSpec>>::new();
    let mut edge_specs = Vec::<Py<PyEdgeSpec>>::new();
    for item in items.iter() {
        if item.is_instance_of::<PyNodeSpec>() {
            node_specs.push(item.extract()?);
        } else if item.is_instance_of::<PyEdgeSpec>() {
            edge_specs.push(item.extract()?);
        } else {
            return Err(PyTypeError::new_err(
                "build() accepts only values returned by node() and edge()",
            ));
        }
    }

    let mut names = BTreeSet::new();
    let mut node_indices = BTreeMap::new();
    let mut builder = HedgeGraphBuilder::<EdgeRecord, NodeRecord, HalfEdgeRecord>::new();
    for spec in &node_specs {
        let spec = spec.borrow(py);
        if let Some(name) = &spec.name {
            if !names.insert(name.clone()) {
                return Err(PyValueError::new_err(format!(
                    "duplicate node name {name:?}"
                )));
            }
        }
        let index = builder.add_node(NodeRecord {
            name: spec.name.clone(),
            data: spec
                .data
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("node specification has been cleared"))?
                .clone_ref(py),
            drawing: copy_drawing(
                py,
                spec.drawing.as_ref().ok_or_else(|| {
                    PyReferenceError::new_err("node specification has been cleared")
                })?,
            )?,
        });
        if node_indices.insert(spec.token, index).is_some() {
            return Err(PyValueError::new_err(
                "the same node specification was passed to build() more than once",
            ));
        }
    }

    names.clear();
    for spec in edge_specs {
        let spec = spec.borrow(py);
        if let Some(name) = &spec.name {
            if !names.insert(name.clone()) {
                return Err(PyValueError::new_err(format!(
                    "duplicate edge name {name:?}"
                )));
            }
        }
        let first = endpoint_record(
            py,
            spec.first
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("edge specification has been cleared"))?,
            &node_indices,
        )?;
        let second = spec
            .second
            .as_ref()
            .map(|endpoint| endpoint_record(py, endpoint, &node_indices))
            .transpose()?;
        let (source, sink) = match (first, second) {
            ((EndpointRole::Source, source), Some((EndpointRole::Sink, sink)))
            | ((EndpointRole::Sink, sink), Some((EndpointRole::Source, source))) => {
                (Some(source), Some(sink))
            }
            ((EndpointRole::Source, source), None) => (Some(source), None),
            ((EndpointRole::Sink, sink), None) => (None, Some(sink)),
            _ => {
                return Err(PyValueError::new_err(
                    "an internal edge needs one source and one sink",
                ));
            }
        };
        let record = EdgeRecord {
            name: spec.name.clone(),
            data: spec
                .data
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("edge specification has been cleared"))?
                .clone_ref(py),
            drawing: copy_drawing(
                py,
                spec.drawing.as_ref().ok_or_else(|| {
                    PyReferenceError::new_err("edge specification has been cleared")
                })?,
            )?,
        };
        match (source, sink) {
            (Some(source), Some(sink)) => builder.add_edge(source, sink, record, spec.orientation),
            (Some(source), None) => {
                builder.add_external_edge(source, record, spec.orientation, Flow::Source)
            }
            (None, Some(sink)) => {
                builder.add_external_edge(sink, record, spec.orientation, Flow::Sink)
            }
            (None, None) => unreachable!("validated edge has an endpoint"),
        }
    }

    Ok(PyGraph::from_state(GraphState {
        graph: PyHedgeGraph::build(builder, node_store),
        name,
        global_data: global_data.unwrap_or_default(),
        codec,
        render_config: match render_config {
            Some(config) => {
                crate::typst::validate_render_config(py, &config)?;
                config
            }
            None => crate::typst::default_render_config(py)?,
        },
        revision: 0,
    }))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyFlow>()?;
    module.add_class::<PyOrientation>()?;
    module.add_class::<PyNodeStore>()?;
    module.add_class::<PyNodeSpec>()?;
    module.add_class::<PyHalfEdgeSpec>()?;
    module.add_class::<PyEdgeSpec>()?;
    module.add_class::<PyGraph>()?;
    module.add_class::<PyNode>()?;
    module.add_class::<PyEdge>()?;
    module.add_class::<PyHalfEdge>()?;
    module.add_function(wrap_pyfunction!(node, module)?)?;
    module.add_function(wrap_pyfunction!(source, module)?)?;
    module.add_function(wrap_pyfunction!(sink, module)?)?;
    module.add_function(wrap_pyfunction!(edge, module)?)?;
    module.add_function(wrap_pyfunction!(build, module)?)?;
    Ok(())
}
