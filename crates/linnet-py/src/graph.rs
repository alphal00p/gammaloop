mod algorithms;
mod generate;
pub use generate::{CancelCheck, GenerationCallbacks, PyBridges, PyPlanReport};
mod keys;
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use linnet::half_edge::builder::{HedgeData, HedgeGraphBuilder};
use linnet::half_edge::involution::{EdgeIndex, Flow, Hedge, HedgePair, Orientation};
use linnet::half_edge::remap::{HedgeRemap, NodeRemap};
use linnet::half_edge::swap::Swap;
use linnet::half_edge::{EdgeAccessors, HedgeGraph, NodeIndex};
use linnet::permutation::Permutation;
use pyo3::class::gc::{PyTraverseError, PyVisit};
use pyo3::exceptions::{PyIndexError, PyKeyError, PyReferenceError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyType};
use pyo3::IntoPyObjectExt;

use crate::dot::{self, PyDotCodec, PyEdgeValue, PyGlobalData, PyNodeValue};
use crate::drawing::{
    copy_drawing, drawing_dict, PyEdgeDrawing, PyHalfEdgeDrawing, PyNodeDrawing, EDGE_FIELDS,
    HEDGE_FIELDS, NODE_FIELDS,
};
use crate::topology::{PyCutPartition, PyCycle, PyDirectionBasis, PySubgraph, PyTraversalTree};

/// Python graphs compact every identification and never open topology
/// checkpoints, so renumbering them cannot fail.
pub(crate) const SETTLED: &str = "Python graphs keep no pending identifications";

static SPEC_ID: AtomicU64 = AtomicU64::new(0);

/// The directed role of a half-edge within its edge.
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyclass_enum
)]
#[pyclass(from_py_object, eq, eq_int, name = "Flow")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
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
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyclass_enum
)]
#[pyclass(from_py_object, eq, eq_int, name = "Orientation")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
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

/// A reusable node description accepted by `Graph` and `Graph.add_node()`.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(unsendable, name = "Node")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub struct PyNodeSpec {
    token: u64,
    name: Option<String>,
    data: Option<Py<PyAny>>,
    drawing: Option<Py<PyDict>>,
}

impl PyNodeSpec {
    pub(crate) fn snapshot(&self, py: Python<'_>) -> PyResult<NodeRecord> {
        Ok(NodeRecord {
            name: self.name.clone(),
            data: self
                .data
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("node specification has been cleared"))?
                .clone_ref(py),
            drawing: copy_drawing(
                py,
                self.drawing.as_ref().ok_or_else(|| {
                    PyReferenceError::new_err("node specification has been cleared")
                })?,
            )?,
        })
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyNodeSpec {
    #[new]
    #[pyo3(signature = (name=None, *, data=None, drawing=None, **options))]
    #[gen_stub(skip)]
    fn new(
        py: Python<'_>,
        name: Option<String>,
        data: Option<Py<PyAny>>,
        drawing: Option<PyRef<'_, PyNodeDrawing>>,
        options: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Self> {
        let options = drawing_dict(py, NODE_FIELDS, options)?;
        let values = if let Some(drawing) = drawing {
            let values = drawing.clone_values(py)?;
            values.bind(py).update(options.bind(py).as_mapping())?;
            values
        } else {
            options
        };
        Ok(Self {
            token: SPEC_ID.fetch_add(1, Ordering::Relaxed),
            name,
            data: Some(data.unwrap_or_else(|| py.None())),
            drawing: Some(values),
        })
    }

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
    fn set_data(slf: &Bound<'_, Self>, value: Py<PyAny>) {
        let old = slf.borrow_mut().data.replace(value);
        drop(old);
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
            Some(name) => format!("Node({name:?})"),
            None => "Node()".to_string(),
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
pub(crate) enum EndpointRole {
    Source,
    Sink,
}

/// A reusable incidence description whose role is chosen by its edge endpoint.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(unsendable, name = "HalfEdge")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub struct PyHalfEdgeSpec {
    node: Option<Py<PyAny>>,
    data: Option<Py<PyAny>>,
    drawing: Option<Py<PyDict>>,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyHalfEdgeSpec {
    #[new]
    #[pyo3(signature = (node, *, data=None, drawing=None, **options))]
    #[gen_stub(skip)]
    fn new(
        py: Python<'_>,
        node: Py<PyAny>,
        data: Option<Py<PyAny>>,
        drawing: Option<PyRef<'_, PyHalfEdgeDrawing>>,
        options: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Self> {
        let target = node.bind(py);
        if !target.is_instance_of::<PyNodeSpec>()
            && !target.is_instance_of::<PyNode>()
            && target.extract::<usize>().is_err()
            && target.extract::<String>().is_err()
        {
            return Err(PyTypeError::new_err(
                "an endpoint must be a Node, NodeView, integer index, or node name",
            ));
        }
        let options = drawing_dict(py, HEDGE_FIELDS, options)?;
        let values = if let Some(drawing) = drawing {
            let values = drawing.clone_values(py)?;
            values.bind(py).update(options.bind(py).as_mapping())?;
            values
        } else {
            options
        };
        Ok(Self {
            node: Some(node),
            data: Some(data.unwrap_or_else(|| py.None())),
            drawing: Some(values),
        })
    }

    #[getter]
    #[gen_stub(override_return_type(type_repr = "_EndpointTarget"))]
    fn node(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.node
            .as_ref()
            .map(|node| node.clone_ref(py))
            .ok_or_else(|| PyReferenceError::new_err("half-edge specification has been cleared"))
    }

    #[getter]
    fn data(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.data
            .as_ref()
            .map(|data| data.clone_ref(py))
            .ok_or_else(|| PyReferenceError::new_err("half-edge specification has been cleared"))
    }

    #[setter]
    fn set_data(slf: &Bound<'_, Self>, value: Py<PyAny>) {
        let old = slf.borrow_mut().data.replace(value);
        drop(old);
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

/// A reusable edge description accepted by `Graph` and `Graph.add_edge()`.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(unsendable, name = "Edge")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub struct PyEdgeSpec {
    name: Option<String>,
    source: Option<Py<PyHalfEdgeSpec>>,
    sink: Option<Py<PyHalfEdgeSpec>>,
    data: Option<Py<PyAny>>,
    drawing: Option<Py<PyDict>>,
    orientation: PyOrientation,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyEdgeSpec {
    #[new]
    #[pyo3(signature = (source=None, sink=None, *, name=None, data=None, orientation=PyOrientation::Default, drawing=None, **options))]
    #[gen_stub(skip)]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        source: Option<Py<PyAny>>,
        sink: Option<Py<PyAny>>,
        name: Option<String>,
        data: Option<Py<PyAny>>,
        orientation: PyOrientation,
        drawing: Option<PyRef<'_, PyEdgeDrawing>>,
        options: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Self> {
        if source.is_none() && sink.is_none() {
            return Err(PyValueError::new_err("an edge needs at least one endpoint"));
        }
        let source = source.map(|value| Self::endpoint(py, value)).transpose()?;
        let sink = sink.map(|value| Self::endpoint(py, value)).transpose()?;
        let options = drawing_dict(py, EDGE_FIELDS, options)?;
        let values = if let Some(drawing) = drawing {
            let values = drawing.clone_values(py)?;
            values.bind(py).update(options.bind(py).as_mapping())?;
            values
        } else {
            options
        };
        Ok(Self {
            name,
            source,
            sink,
            data: Some(data.unwrap_or_else(|| py.None())),
            drawing: Some(values),
            orientation,
        })
    }

    #[getter]
    fn source(&self, py: Python<'_>) -> Option<Py<PyHalfEdgeSpec>> {
        self.source.as_ref().map(|value| value.clone_ref(py))
    }

    #[getter]
    fn sink(&self, py: Python<'_>) -> Option<Py<PyHalfEdgeSpec>> {
        self.sink.as_ref().map(|value| value.clone_ref(py))
    }

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
    fn set_data(slf: &Bound<'_, Self>, value: Py<PyAny>) {
        let old = slf.borrow_mut().data.replace(value);
        drop(old);
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
        visit.call(&self.source)?;
        visit.call(&self.sink)?;
        visit.call(&self.data)?;
        visit.call(&self.drawing)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.source = None;
        self.sink = None;
        self.data = None;
        self.drawing = None;
    }
}

impl PyEdgeSpec {
    fn endpoint(py: Python<'_>, value: Py<PyAny>) -> PyResult<Py<PyHalfEdgeSpec>> {
        if value.bind(py).is_instance_of::<PyHalfEdgeSpec>() {
            value.extract(py).map_err(Into::into)
        } else {
            Py::new(py, PyHalfEdgeSpec::new(py, value, None, None, None)?)
        }
    }
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

/// The Linnet graph behind a Python graph.
pub(crate) type PyHedgeGraph = HedgeGraph<EdgeRecord, NodeRecord, HalfEdgeRecord>;

pub(crate) struct IncrementalEdgeRecords {
    pub(crate) source: Option<HedgeData<HalfEdgeRecord>>,
    pub(crate) sink: Option<HedgeData<HalfEdgeRecord>>,
    pub(crate) edge: EdgeRecord,
    pub(crate) orientation: PyOrientation,
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
    pub(crate) element_epoch: u64,
}

impl GraphState {
    pub(crate) fn copy_graph(&self, py: Python<'_>) -> PyResult<PyHedgeGraph> {
        self.graph.map_data_ref_result(
            |_, _, node| node.copy(py),
            |_, _, _, edge| edge.map_result(|edge| edge.copy(py)),
            |(_, half_edge)| half_edge.copy(py),
        )
    }

    pub(crate) fn copy(&self, py: Python<'_>) -> PyResult<Self> {
        let revision = self.revision;
        let mut state = self.derived(py, self.copy_graph(py)?)?;
        state.revision = revision;
        state.element_epoch = self.element_epoch;
        Ok(state)
    }

    /// Move node `order[i]` to position `i`, keeping every payload with its element.
    pub(crate) fn reorder_nodes(&mut self, order: Vec<usize>) -> PyResult<()> {
        let order = checked_order(order, self.graph.node_bound(), "node")?;
        (self.graph)
            .permute_nodes(&NodeRemap::from_permutation(&Permutation::from_inv(order)))
            .expect(SETTLED);
        Ok(())
    }

    /// Move edge `order[i]` to position `i`, and its half-edges along, source first.
    pub(crate) fn reorder_edges(&mut self, order: Vec<usize>) -> PyResult<()> {
        let order = checked_order(order, self.graph.n_edges(), "edge")?;
        let hedge_order = order
            .iter()
            .flat_map(|index| match self.graph[&EdgeIndex(*index)].1 {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    [Some(source.0), Some(sink.0)]
                }
                HedgePair::Unpaired { hedge, .. } => [Some(hedge.0), None],
            })
            .flatten()
            .collect();
        (self.graph)
            .permute_hedges(&HedgeRemap::from_permutation(&Permutation::from_inv(
                hedge_order,
            )))
            .expect(SETTLED);
        <PyHedgeGraph as Swap<EdgeIndex>>::permute(&mut self.graph, &Permutation::from_inv(order));
        Ok(())
    }

    pub(crate) fn derived(&self, py: Python<'_>, graph: PyHedgeGraph) -> PyResult<Self> {
        Ok(Self {
            graph,
            name: self.name.clone(),
            global_data: self.global_data.clone(),
            codec: self.codec.as_ref().map(|codec| codec.clone_ref(py)),
            render_config: crate::typst::render_config_copy(py, &self.render_config)?,
            revision: 0,
            element_epoch: 0,
        })
    }
}

/// An owned topology with arbitrary element data and typed rendering configuration.
///
/// Examples
/// --------
/// ```python
/// from symbolica.community import graph as linnet
/// graph = linnet.Graph(nodes=[linnet.Node("source"), linnet.Node("sink")])
/// ```
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(unsendable, name = "Graph")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub struct PyGraph {
    pub(crate) state: RefCell<Option<GraphState>>,
}

impl PyGraph {
    // PyO3 calls this only for supplied arguments, preserving explicit Python None.
    fn supplied_data(value: &Bound<'_, PyAny>) -> PyResult<Option<Py<PyAny>>> {
        Ok(Some(value.clone().unbind()))
    }

    pub(crate) fn from_state(state: GraphState) -> Self {
        Self {
            state: RefCell::new(Some(state)),
        }
    }

    /// Construct the common Python graph from a native Linnet graph, mapping every node, edge
    /// and half-edge payload to its Python data. Structure and orientations are kept.
    pub fn from_hedge_graph<E, V, H>(
        py: Python<'_>,
        graph: &HedgeGraph<E, V, H>,
        mut node: impl FnMut(&V) -> PyResult<Py<PyAny>>,
        mut edge: impl FnMut(&E) -> PyResult<Py<PyAny>>,
        mut hedge: impl FnMut(&H) -> PyResult<Py<PyAny>>,
    ) -> PyResult<Self> {
        let graph = graph.map_data_ref_result(
            |_, _, data| {
                Ok::<_, PyErr>(NodeRecord {
                    name: None,
                    data: node(data)?,
                    drawing: drawing_dict(py, NODE_FIELDS, None)?,
                })
            },
            |_, _, _, data| {
                data.map_result(|data| {
                    Ok(EdgeRecord {
                        name: None,
                        data: edge(data)?,
                        drawing: drawing_dict(py, EDGE_FIELDS, None)?,
                    })
                })
            },
            |(_, data)| {
                Ok(HalfEdgeRecord {
                    data: hedge(data)?,
                    drawing: drawing_dict(py, HEDGE_FIELDS, None)?,
                })
            },
        )?;
        Ok(Self::from_state(GraphState {
            graph,
            name: None,
            global_data: PyGlobalData::default(),
            codec: None,
            render_config: crate::typst::default_render_config(py)?,
            revision: 0,
            element_epoch: 0,
        }))
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

    pub(crate) fn element_epoch(&self) -> PyResult<u64> {
        self.state
            .borrow()
            .as_ref()
            .map(|state| state.element_epoch)
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))
    }

    pub(crate) fn check_element_epoch(&self, epoch: u64) -> PyResult<()> {
        if self.element_epoch()? != epoch {
            return Err(PyReferenceError::new_err(
                "this graph element view is stale after a topology transformation",
            ));
        }
        Ok(())
    }

    pub(crate) fn resolve_node(&self, key: &Bound<'_, PyAny>) -> PyResult<usize> {
        // Integer conversion can invoke Python __index__ and mutate this graph.
        let index = key.extract::<usize>().ok();
        let name = if index.is_none() {
            key.extract::<String>().ok()
        } else {
            None
        };
        let state = self.state.borrow();
        let state = state
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        if let Some(index) = index {
            return (index < state.graph.node_bound())
                .then_some(index)
                .ok_or_else(|| PyIndexError::new_err("node index out of range"));
        }
        if let Some(name) = name {
            return state
                .graph
                .iter_node_ids()
                .find_map(|index| {
                    (state.graph[index].name.as_deref() == Some(name.as_str())).then_some(index.0)
                })
                .ok_or_else(|| PyKeyError::new_err(format!("unknown node {name:?}")));
        }
        Err(PyTypeError::new_err("node key must be an integer or name"))
    }

    pub(crate) fn resolve_edge(&self, key: &Bound<'_, PyAny>) -> PyResult<usize> {
        // Integer conversion can invoke Python __index__ and mutate this graph.
        let index = key.extract::<usize>().ok();
        let name = if index.is_none() {
            key.extract::<String>().ok()
        } else {
            None
        };
        let state = self.state.borrow();
        let state = state
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        if let Some(index) = index {
            return (index < state.graph.n_edges())
                .then_some(index)
                .ok_or_else(|| PyIndexError::new_err("edge index out of range"));
        }
        if let Some(name) = name {
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

    pub(crate) fn resolve_half_edge(
        &self,
        py: Python<'_>,
        owner: &Py<PyGraph>,
        key: &Bound<'_, PyAny>,
    ) -> PyResult<usize> {
        let epoch = self.element_epoch()?;
        if let Ok(view) = key.extract::<PyRef<'_, PyHalfEdge>>() {
            return view.index_for_graph(py, owner, epoch);
        }
        let index = key.extract::<usize>().map_err(|_| {
            PyTypeError::new_err("half-edge key must be an integer or HalfEdgeView")
        })?;
        let len = self
            .state
            .borrow()
            .as_ref()
            .expect("checked")
            .graph
            .n_hedges();
        (index < len)
            .then_some(index)
            .ok_or_else(|| PyIndexError::new_err("half-edge index out of range"))
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
        #[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
        #[pyclass(unsendable, name = $python)]
        #[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
        #[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
        pub struct $rust {
            graph: Option<Py<PyGraph>>,
            index: usize,
            epoch: u64,
        }

        impl $rust {
            pub(crate) fn new(graph: Py<PyGraph>, index: usize, epoch: u64) -> Self {
                Self {
                    graph: Some(graph),
                    index,
                    epoch,
                }
            }

            fn owner<'py>(&self, py: Python<'py>) -> PyResult<&Bound<'py, PyGraph>> {
                let graph = self
                    .graph
                    .as_ref()
                    .ok_or_else(|| PyReferenceError::new_err("graph view has been cleared"))?;
                graph.borrow(py).check_element_epoch(self.epoch)?;
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
    "NodeView",
    Node
);
graph_view!(
    "A live edge view whose data and drawing metadata update its graph.",
    PyEdge,
    "EdgeView",
    Edge
);
graph_view!(
    "A live half-edge view whose data and drawing metadata update its graph.",
    PyHalfEdge,
    "HalfEdgeView",
    Hedge
);

impl PyNode {
    pub(crate) fn index_for_graph(
        &self,
        py: Python<'_>,
        graph: &Py<PyGraph>,
        epoch: u64,
    ) -> PyResult<usize> {
        let owner = self.owner(py)?;
        if self.epoch != epoch || !owner.is(graph.bind(py)) {
            return Err(PyValueError::new_err(
                "node endpoint belongs to a different graph",
            ));
        }
        Ok(self.index)
    }
}

impl PyHalfEdge {
    pub(crate) fn index_for_graph(
        &self,
        py: Python<'_>,
        graph: &Py<PyGraph>,
        epoch: u64,
    ) -> PyResult<usize> {
        let owner = self.owner(py)?;
        if self.epoch != epoch || !owner.is(graph.bind(py)) {
            return Err(PyValueError::new_err(
                "half-edge belongs to a different graph",
            ));
        }
        Ok(self.index)
    }
}

impl PyEdge {
    pub(crate) fn index_for_graph(
        &self,
        py: Python<'_>,
        graph: &Py<PyGraph>,
        epoch: u64,
    ) -> PyResult<usize> {
        let owner = self.owner(py)?;
        if self.epoch != epoch || !owner.is(graph.bind(py)) {
            return Err(PyValueError::new_err("edge belongs to a different graph"));
        }
        Ok(self.index)
    }

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
                    PyHalfEdge::new(graph_obj.clone_ref(py), index.0, self.epoch),
                )
            })
            .transpose()
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
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
        let old = {
            let view = slf.borrow();
            let graph = view.owner(slf.py())?.borrow();
            let mut state = graph.state.borrow_mut();
            std::mem::replace(
                &mut state.as_mut().expect("checked").graph[NodeIndex(view.index)].data,
                value,
            )
        };
        drop(old);
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
            self.epoch,
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
                    PyHalfEdge::new(graph_obj.clone_ref(py), hedge.0, self.epoch),
                )
            })
            .collect()
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(match self.name(py)? {
            Some(name) => format!("NodeView({}: {name:?})", self.index),
            None => format!("NodeView({})", self.index),
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

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
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
        let old = {
            let view = slf.borrow();
            let graph = view.owner(slf.py())?.borrow();
            let mut state = graph.state.borrow_mut();
            std::mem::replace(
                &mut state.as_mut().expect("checked").graph[EdgeIndex(view.index)].data,
                value,
            )
        };
        drop(old);
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
            self.epoch,
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
            Some(name) => format!("EdgeView({}: {name:?})", self.index),
            None => format!("EdgeView({})", self.index),
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

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
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
        let old = {
            let view = slf.borrow();
            let graph = view.owner(slf.py())?.borrow();
            let mut state = graph.state.borrow_mut();
            std::mem::replace(
                &mut state.as_mut().expect("checked").graph[Hedge(view.index)].data,
                value,
            )
        };
        drop(old);
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
            self.epoch,
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
        Py::new(py, PyNode::new(graph_obj.clone_ref(py), node.0, self.epoch))
    }

    #[getter]
    fn edge(&self, py: Python<'_>) -> PyResult<Py<PyEdge>> {
        let graph_obj = self.clone_graph(py)?;
        let graph = graph_obj.borrow(py);
        let state = graph.state.borrow();
        let edge = state.as_ref().expect("checked").graph[&Hedge(self.index)];
        Py::new(py, PyEdge::new(graph_obj.clone_ref(py), edge.0, self.epoch))
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
                    PyHalfEdge::new(graph_obj.clone_ref(py), pair.0, self.epoch),
                )
            })
            .transpose()
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.owner(py)?;
        Ok(format!("HalfEdgeView({})", self.index))
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
    epoch: u64,
    index: ViewIndex,
) -> PyResult<Py<PyAny>> {
    Ok(match index {
        ViewIndex::Node(index) => {
            Py::new(py, PyNode::new(graph.clone_ref(py), index, epoch))?.into_any()
        }
        ViewIndex::Edge(index) => {
            Py::new(py, PyEdge::new(graph.clone_ref(py), index, epoch))?.into_any()
        }
        ViewIndex::Hedge(index) => {
            Py::new(py, PyHalfEdge::new(graph.clone_ref(py), index, epoch))?.into_any()
        }
    })
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyGraph {
    /// Build a graph once from reusable descriptions. Referenced Node descriptions
    /// are included automatically; named and indexed endpoints must resolve.
    #[new]
    #[pyo3(signature = (*, nodes=None, edges=None, name=None, global_data=None, codec=None, render_config=None))]
    fn new(
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "typing.Iterable[Node] | None", imports=("typing")))]
        nodes: Option<&Bound<'_, PyAny>>,
        #[gen_stub(override_type(type_repr = "typing.Iterable[Edge] | None", imports=("typing")))]
        edges: Option<&Bound<'_, PyAny>>,
        name: Option<String>,
        global_data: Option<PyGlobalData>,
        codec: Option<Py<PyDotCodec>>,
        #[gen_stub(override_type(type_repr = "RenderSettings | None"))] render_config: Option<
            Py<PyAny>,
        >,
    ) -> PyResult<Self> {
        let nodes = nodes
            .map(|items| {
                items
                    .try_iter()?
                    .map(|item| item?.extract().map_err(Into::into))
                    .collect::<PyResult<Vec<_>>>()
            })
            .transpose()?
            .unwrap_or_default();
        let edges = edges
            .map(|items| {
                items
                    .try_iter()?
                    .map(|item| item?.extract().map_err(Into::into))
                    .collect::<PyResult<Vec<_>>>()
            })
            .transpose()?
            .unwrap_or_default();
        Self::from_specs(py, nodes, edges, name, global_data, codec, render_config)
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
    #[gen_stub(override_return_type(type_repr = "RenderSettings"))]
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
        #[gen_stub(override_type(type_repr = "RenderSettings"))] value: Py<PyAny>,
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
            .n_live())
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

    #[pyo3(signature = (key))]
    fn node(
        slf: Py<PyGraph>,
        #[gen_stub(override_type(type_repr="builtins.int | builtins.str", imports=("builtins")))]
        key: &Bound<'_, PyAny>,
        py: Python<'_>,
    ) -> PyResult<Py<PyNode>> {
        let index = slf.borrow(py).resolve_node(key)?;
        let epoch = slf.borrow(py).element_epoch()?;
        Py::new(py, PyNode::new(slf.clone_ref(py), index, epoch))
    }

    #[pyo3(signature = (key))]
    fn edge(
        slf: Py<PyGraph>,
        #[gen_stub(override_type(type_repr="builtins.int | builtins.str", imports=("builtins")))]
        key: &Bound<'_, PyAny>,
        py: Python<'_>,
    ) -> PyResult<Py<PyEdge>> {
        let index = slf.borrow(py).resolve_edge(key)?;
        let epoch = slf.borrow(py).element_epoch()?;
        Py::new(py, PyEdge::new(slf.clone_ref(py), index, epoch))
    }

    #[pyo3(signature = (index))]
    fn half_edge(slf: Py<PyGraph>, index: usize, py: Python<'_>) -> PyResult<Py<PyHalfEdge>> {
        let graph = slf.borrow(py);
        let epoch = graph.element_epoch()?;
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
        Py::new(py, PyHalfEdge::new(slf.clone_ref(py), index, epoch))
    }

    #[pyo3(signature = ())]
    fn nodes(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<Vec<Py<PyNode>>> {
        let graph = slf.borrow(py);
        let epoch = graph.element_epoch()?;
        let len = graph
            .state
            .borrow()
            .as_ref()
            .expect("checked")
            .graph
            .node_bound();
        drop(graph);
        (0..len)
            .map(|index| Py::new(py, PyNode::new(slf.clone_ref(py), index, epoch)))
            .collect()
    }

    #[pyo3(signature = ())]
    fn edges(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<Vec<Py<PyEdge>>> {
        let graph = slf.borrow(py);
        let epoch = graph.element_epoch()?;
        let len = graph
            .state
            .borrow()
            .as_ref()
            .expect("checked")
            .graph
            .n_edges();
        drop(graph);
        (0..len)
            .map(|index| Py::new(py, PyEdge::new(slf.clone_ref(py), index, epoch)))
            .collect()
    }

    #[pyo3(signature = ())]
    fn half_edges(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<Vec<Py<PyHalfEdge>>> {
        let graph = slf.borrow(py);
        let epoch = graph.element_epoch()?;
        let len = graph
            .state
            .borrow()
            .as_ref()
            .expect("checked")
            .graph
            .n_hedges();
        drop(graph);
        (0..len)
            .map(|index| Py::new(py, PyHalfEdge::new(slf.clone_ref(py), index, epoch)))
            .collect()
    }

    fn reorder_nodes(&self, order: Vec<usize>) -> PyResult<()> {
        let mut state = self.state.borrow_mut();
        let state = state
            .as_mut()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        state.reorder_nodes(order)?;
        state.revision += 1;
        state.element_epoch += 1;
        Ok(())
    }

    fn reorder_edges(&self, order: Vec<usize>) -> PyResult<()> {
        let mut state = self.state.borrow_mut();
        let state = state
            .as_mut()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        state.reorder_edges(order)?;
        state.revision += 1;
        state.element_epoch += 1;
        Ok(())
    }

    /// Change an edge's direction and invalidate existing graph-element views.
    fn set_orientation(
        &self,
        #[gen_stub(override_type(type_repr="builtins.int | builtins.str", imports=("builtins")))]
        key: &Bound<'_, PyAny>,
        orientation: PyOrientation,
    ) -> PyResult<()> {
        let index = self.resolve_edge(key)?;
        let mut state = self.state.borrow_mut();
        let state = state.as_mut().expect("checked");
        state
            .graph
            .set_orientation(EdgeIndex(index), orientation.into());
        state.revision += 1;
        state.element_epoch += 1;
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
        state.element_epoch += 1;
        Ok(())
    }

    #[pyo3(signature = (*, node=None, edge=None, source=None, sink=None))]
    fn map(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Callable[[NodeView], typing.Any] | None", imports=("typing")))]
        node: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[EdgeView], typing.Any] | None", imports=("typing")))]
        edge: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdgeView], typing.Any] | None", imports=("typing")))]
        source: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdgeView], typing.Any] | None", imports=("typing")))]
        sink: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        let (revision, epoch, mut mapped) = {
            let graph = slf.borrow(py);
            let state_ref = graph.state.borrow();
            let state = state_ref
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
            (state.revision, state.element_epoch, state.copy(py)?)
        };

        for index in 0..mapped.graph.node_bound() {
            if let Some(callback) = &node {
                mapped.graph[NodeIndex(index)].data =
                    callback.call1(py, (view_for(py, &slf, epoch, ViewIndex::Node(index))?,))?;
            }
        }
        let source_callback = source.as_ref();
        let sink_callback = sink.as_ref();
        for index in 0..mapped.graph.n_edges() {
            if let Some(callback) = &edge {
                mapped.graph[EdgeIndex(index)].data =
                    callback.call1(py, (view_for(py, &slf, epoch, ViewIndex::Edge(index))?,))?;
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
                    mapped.graph[hedge].data = callback
                        .call1(py, (view_for(py, &slf, epoch, ViewIndex::Hedge(hedge.0))?,))?;
                }
            }
        }
        slf.borrow(py).check_revision(revision)?;
        mapped.revision = 0;
        mapped.element_epoch = 0;
        Ok(Self::from_state(mapped))
    }

    #[classmethod]
    fn from_dot(
        _cls: &Bound<'_, PyType>,
        py: Python<'_>,
        source: &str,
        codec: Py<PyDotCodec>,
    ) -> PyResult<Self> {
        dot::decode_dot(py, source, codec)
    }

    #[classmethod]
    fn from_dot_set(
        _cls: &Bound<'_, PyType>,
        py: Python<'_>,
        source: &str,
        codec: Py<PyDotCodec>,
    ) -> PyResult<Vec<Self>> {
        dot::decode_dot_set(py, source, codec)
    }

    #[classmethod]
    fn from_dot_file(
        _cls: &Bound<'_, PyType>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="builtins.str | os.PathLike[builtins.str]", imports=("builtins", "os")))]
        path: PathBuf,
        codec: Py<PyDotCodec>,
    ) -> PyResult<Self> {
        let source = std::fs::read_to_string(&path).map_err(|error| {
            PyValueError::new_err(format!("failed to read {}: {error}", path.display()))
        })?;
        dot::decode_dot(py, &source, codec)
    }

    #[pyo3(signature = (codec=None))]
    fn to_dot(&self, py: Python<'_>, codec: Option<Py<PyDotCodec>>) -> PyResult<String> {
        dot::encode_graph(py, self, codec)
    }

    /// Snapshot a drawing, using native SVG geometry for ordinary graph settings.
    #[pyo3(signature = (*, config=None))]
    fn render(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "RenderSettings | None"))] config: Option<
            &Bound<'_, PyAny>,
        >,
    ) -> PyResult<crate::PyDiagramRender> {
        crate::render::render_graph(py, &slf, config, None)
    }

    /// Render one static SVG page; ``interactive=True`` adds the hover,
    /// selection and keyboard targets of notebook display.
    #[pyo3(signature = (*, config=None, interactive=false))]
    fn to_svg(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "RenderSettings | None"))] config: Option<
            &Bound<'_, PyAny>,
        >,
        interactive: bool,
    ) -> PyResult<String> {
        crate::render::graph_to_svg(py, &slf, config, interactive)
    }

    #[pyo3(signature = ())]
    fn _repr_svg_(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<String> {
        crate::render::graph_to_svg(py, &slf, None, false)
    }

    /// Display the SVG with hover details and local interactive selection.
    #[pyo3(signature = ())]
    fn _repr_html_(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<String> {
        crate::render::graph_to_svg(py, &slf, None, true)
    }

    #[pyo3(name = "full_subgraph", signature = ())]
    fn py_full_subgraph(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<PySubgraph> {
        Self::full_subgraph(slf, py)
    }

    #[pyo3(name = "empty_subgraph", signature = ())]
    fn py_empty_subgraph(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<PySubgraph> {
        Self::empty_subgraph(slf, py)
    }

    #[pyo3(
        name = "subgraph",
        signature = (*, nodes=Vec::new(), edges=Vec::new(), half_edges=Vec::new())
    )]
    fn py_subgraph(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Sequence[builtins.int | builtins.str]", imports=("builtins", "typing")))]
        nodes: Vec<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Sequence[builtins.int | builtins.str]", imports=("builtins", "typing")))]
        edges: Vec<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Sequence[builtins.int]", imports=("builtins", "typing")))]
        half_edges: Vec<usize>,
    ) -> PyResult<PySubgraph> {
        Self::subgraph(slf, py, nodes, edges, half_edges)
    }

    #[pyo3(
        name = "filter",
        signature = (*, node=None, edge=None, half_edge=None)
    )]
    fn py_filter(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Callable[[NodeView], builtins.bool] | None", imports=("builtins", "typing")))]
        node: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[EdgeView], builtins.bool] | None", imports=("builtins", "typing")))]
        edge: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdgeView], builtins.bool] | None", imports=("builtins", "typing")))]
        half_edge: Option<Py<PyAny>>,
    ) -> PyResult<PySubgraph> {
        Self::filter(slf, py, node, edge, half_edge)
    }

    #[pyo3(name = "nodes_of", signature = (subgraph))]
    fn py_nodes_of(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
    ) -> PyResult<Vec<Py<PyNode>>> {
        Self::nodes_of(slf, py, subgraph)
    }

    #[pyo3(name = "edges_of", signature = (subgraph))]
    fn py_edges_of(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
    ) -> PyResult<Vec<Py<PyEdge>>> {
        Self::edges_of(slf, py, subgraph)
    }

    #[pyo3(name = "half_edges_of", signature = (subgraph))]
    fn py_half_edges_of(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
    ) -> PyResult<Vec<Py<PyHalfEdge>>> {
        Self::half_edges_of(slf, py, subgraph)
    }

    /// Return selected half-edges that lie on the selection's internal boundary.
    #[pyo3(name = "internal_boundary", signature = (subgraph))]
    fn py_internal_boundary(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
    ) -> PyResult<PySubgraph> {
        Self::internal_boundary(slf, py, subgraph)
    }

    /// Return all boundary half-edges incident to nodes touched by the selection.
    #[pyo3(name = "boundary", signature = (subgraph))]
    fn py_boundary(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
    ) -> PyResult<PySubgraph> {
        Self::boundary(slf, py, subgraph)
    }

    /// Return dangling/external half-edges as a composable structural selection.
    #[pyo3(name = "external_half_edges", signature = ())]
    fn py_external_half_edges(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<PySubgraph> {
        Self::external_half_edges(slf, py)
    }

    #[pyo3(name = "connected_components", signature = (subgraph=None))]
    fn py_connected_components(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<Vec<PySubgraph>> {
        Self::connected_components(slf, py, subgraph)
    }

    #[pyo3(name = "count_connected_components", signature = (subgraph=None))]
    fn py_count_connected_components(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<usize> {
        Self::count_connected_components(slf, py, subgraph)
    }

    #[pyo3(name = "is_connected", signature = (subgraph=None))]
    fn py_is_connected(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<bool> {
        Self::is_connected(slf, py, subgraph)
    }

    #[pyo3(name = "cyclomatic_number", signature = (subgraph=None))]
    fn py_cyclomatic_number(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<usize> {
        Self::cyclomatic_number(slf, py, subgraph)
    }

    /// Test directed reachability within an optional structural selection.
    ///
    /// Both endpoints must belong to the selected subgraph; otherwise this raises `ValueError`.
    #[pyo3(
        name = "is_reachable",
        signature = (source, target, *, subgraph=None, direction=PyDirectionBasis::Underlying)
    )]
    fn py_is_reachable(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="NodeView | builtins.int | builtins.str", imports=("builtins")))]
        source: &Bound<'_, PyAny>,
        #[gen_stub(override_type(type_repr="NodeView | builtins.int | builtins.str", imports=("builtins")))]
        target: &Bound<'_, PyAny>,
        subgraph: Option<&PySubgraph>,
        direction: PyDirectionBasis,
    ) -> PyResult<bool> {
        Self::is_reachable(slf, py, source, target, subgraph, direction)
    }

    /// Return a deterministic topological order for the selected directed graph.
    ///
    /// Explicitly selected zero-crown nodes appear first in ascending graph order.
    #[pyo3(
        name = "topological_order",
        signature = (*, subgraph=None, direction=PyDirectionBasis::Underlying)
    )]
    fn py_topological_order(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
        direction: PyDirectionBasis,
    ) -> PyResult<Vec<Py<PyNode>>> {
        Self::topological_order(slf, py, subgraph, direction)
    }

    /// Return a new graph with directionally redundant DAG edges removed.
    #[pyo3(
        name = "transitive_reduction",
        signature = (*, direction=PyDirectionBasis::Underlying)
    )]
    fn py_transitive_reduction(
        slf: Py<PyGraph>,
        py: Python<'_>,
        direction: PyDirectionBasis,
    ) -> PyResult<PyGraph> {
        Self::transitive_reduction(slf, py, direction)
    }

    #[pyo3(
        name = "depth_first_traverse",
        signature = (root, *, subgraph=None, include=None)
    )]
    fn py_depth_first_traverse(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="builtins.int | builtins.str", imports=("builtins")))]
        root: &Bound<'_, PyAny>,
        subgraph: Option<&PySubgraph>,
        include: Option<usize>,
    ) -> PyResult<PyTraversalTree> {
        Self::depth_first_traverse(slf, py, root, subgraph, include)
    }

    #[pyo3(
        name = "breadth_first_traverse",
        signature = (root, *, subgraph=None, include=None)
    )]
    fn py_breadth_first_traverse(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="builtins.int | builtins.str", imports=("builtins")))]
        root: &Bound<'_, PyAny>,
        subgraph: Option<&PySubgraph>,
        include: Option<usize>,
    ) -> PyResult<PyTraversalTree> {
        Self::breadth_first_traverse(slf, py, root, subgraph, include)
    }

    #[pyo3(name = "bridges", signature = (subgraph=None))]
    fn py_bridges(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<PySubgraph> {
        Self::bridges(slf, py, subgraph)
    }

    #[pyo3(name = "cycle_basis", signature = (subgraph=None))]
    fn py_cycle_basis(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<(Vec<PyCycle>, PySubgraph)> {
        Self::cycle_basis(slf, py, subgraph)
    }

    #[pyo3(name = "all_spanning_forests", signature = (subgraph=None))]
    fn py_all_spanning_forests(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<Vec<PySubgraph>> {
        Self::all_spanning_forests(slf, py, subgraph)
    }

    /// Enumerate inclusion-minimal cutsets with an inclusive boundary-size range.
    ///
    /// Bounds default to `[1, unbounded]`; `min_size=0` is rejected because a bond is non-empty.
    #[pyo3(
        name = "all_bonds",
        signature = (*, subgraph=None, min_size=None, max_size=None)
    )]
    fn py_all_bonds(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
        min_size: Option<usize>,
        max_size: Option<usize>,
    ) -> PyResult<Vec<PySubgraph>> {
        Self::all_bonds(slf, py, subgraph, min_size, max_size)
    }

    /// Find one native inclusion-minimal cutset within an inclusive size range.
    #[pyo3(
        name = "find_bond",
        signature = (*, subgraph=None, min_size=None, max_size=None)
    )]
    fn py_find_bond(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
        min_size: Option<usize>,
        max_size: Option<usize>,
    ) -> PyResult<Option<PySubgraph>> {
        Self::find_bond(slf, py, subgraph, min_size, max_size)
    }

    /// Enumerate circuits only when the complete cycle-space candidate count fits the bound.
    ///
    /// `max_results` bounds all non-empty combinations of a cycle basis, not only
    /// the combinations that ultimately form circuits.
    #[pyo3(
        name = "all_cycles",
        signature = (*, max_results, subgraph=None)
    )]
    fn py_all_cycles(
        slf: Py<PyGraph>,
        py: Python<'_>,
        max_results: usize,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<Vec<PyCycle>> {
        Self::all_cycles(slf, py, max_results, subgraph)
    }

    /// Enumerate tadpole components after identifying a non-empty terminal node set.
    #[pyo3(name = "tadpoles", signature = (externals))]
    fn py_tadpoles(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Sequence[NodeView | builtins.int | builtins.str]", imports=("builtins", "typing")))]
        externals: Vec<Py<PyAny>>,
    ) -> PyResult<Vec<PySubgraph>> {
        Self::tadpoles(slf, py, externals)
    }

    /// Enumerate all native separating partitions between two disjoint, non-empty node groups.
    #[pyo3(name = "all_cuts", signature = (source, target))]
    fn py_all_cuts(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Sequence[builtins.int | builtins.str]", imports=("builtins", "typing")))]
        source: Vec<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Sequence[builtins.int | builtins.str]", imports=("builtins", "typing")))]
        target: Vec<Py<PyAny>>,
    ) -> PyResult<Vec<PyCutPartition>> {
        Self::all_cuts(slf, py, source, target)
    }

    /// Append a node description or construct one from a name, payload, and drawing.
    /// Existing element views remain valid; topology selections must be recomputed.
    #[pyo3(name = "add_node", signature = (node=None, *, name=None, data=None, drawing=None, **options))]
    #[gen_stub(skip)]
    fn py_add_node(
        slf: Py<PyGraph>,
        py: Python<'_>,
        node: Option<Py<PyAny>>,
        name: Option<String>,
        #[pyo3(from_py_with = PyGraph::supplied_data)] data: Option<Py<PyAny>>,
        drawing: Option<PyRef<'_, PyNodeDrawing>>,
        options: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Py<PyNode>> {
        let spec = if let Some(node) = node {
            if node.bind(py).is_instance_of::<PyNodeSpec>() {
                if name.is_some()
                    || data.is_some()
                    || drawing.is_some()
                    || options.is_some_and(|options| !options.is_empty())
                {
                    return Err(PyTypeError::new_err(
                        "a Node description cannot be combined with node options",
                    ));
                }
                node.extract::<Py<PyNodeSpec>>(py)?
            } else {
                if name.is_some() {
                    return Err(PyTypeError::new_err("node name was supplied twice"));
                }
                let name = node.extract::<String>(py).map_err(|_| {
                    PyTypeError::new_err("add_node() expects a Node description or a node name")
                })?;
                Py::new(py, PyNodeSpec::new(py, Some(name), data, drawing, options)?)?
            }
        } else {
            Py::new(py, PyNodeSpec::new(py, name, data, drawing, options)?)?
        };
        Self::add_node(slf, py, spec.borrow(py))
    }

    /// Append an edge description or construct an internal or dangling edge.
    /// Node endpoints create fresh incidences; HalfEdge endpoints supply richer data.
    #[pyo3(name = "add_edge", signature = (source=None, sink=None, *, name=None, data=None, orientation=None, drawing=None, **options))]
    #[gen_stub(skip)]
    #[allow(clippy::too_many_arguments)]
    fn py_add_edge(
        slf: Py<PyGraph>,
        py: Python<'_>,
        source: Option<Py<PyAny>>,
        sink: Option<Py<PyAny>>,
        name: Option<String>,
        #[pyo3(from_py_with = PyGraph::supplied_data)] data: Option<Py<PyAny>>,
        orientation: Option<PyOrientation>,
        drawing: Option<PyRef<'_, PyEdgeDrawing>>,
        options: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Py<PyEdge>> {
        let spec = if source
            .as_ref()
            .is_some_and(|value| value.bind(py).is_instance_of::<PyEdgeSpec>())
        {
            if sink.is_some()
                || name.is_some()
                || data.is_some()
                || orientation.is_some()
                || drawing.is_some()
                || options.is_some_and(|options| !options.is_empty())
            {
                return Err(PyTypeError::new_err(
                    "an Edge description cannot be combined with edge options",
                ));
            }
            source.expect("checked").extract::<Py<PyEdgeSpec>>(py)?
        } else {
            Py::new(
                py,
                PyEdgeSpec::new(
                    py,
                    source,
                    sink,
                    name,
                    data,
                    orientation.unwrap_or(PyOrientation::Default),
                    drawing,
                    options,
                )?,
            )?
        };
        Self::add_edge(slf, py, spec.borrow(py))
    }

    /// Split a paired edge at one half-edge, preserving the opposite edge value.
    #[pyo3(
        name = "split_edge",
        signature = (at, replacement, *, name=None, orientation=PyOrientation::Default)
    )]
    fn py_split_edge(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "_HalfEdgeTarget"))] at: &Bound<'_, PyAny>,
        replacement: PyRef<'_, PyEdgeValue>,
        name: Option<String>,
        orientation: PyOrientation,
    ) -> PyResult<(Py<PyEdge>, Py<PyEdge>)> {
        Self::split_edge(slf, py, at, replacement, name, orientation)
    }

    /// Connect two dangling half-edges; the first becomes the source endpoint.
    #[pyo3(
        name = "connect",
        signature = (source, sink, replacement, *, name=None, orientation=PyOrientation::Default)
    )]
    fn py_connect(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "_HalfEdgeTarget"))] source: &Bound<'_, PyAny>,
        #[gen_stub(override_type(type_repr = "_HalfEdgeTarget"))] sink: &Bound<'_, PyAny>,
        replacement: PyRef<'_, PyEdgeValue>,
        name: Option<String>,
        orientation: PyOrientation,
    ) -> PyResult<Py<PyEdge>> {
        Self::connect(slf, py, source, sink, replacement, name, orientation)
    }

    /// Copy a structural subgraph into an independent graph.
    #[pyo3(name = "concretize", signature = (subgraph))]
    fn py_concretize(slf: Py<PyGraph>, py: Python<'_>, subgraph: &PySubgraph) -> PyResult<PyGraph> {
        Self::concretize(slf, py, subgraph)
    }

    /// Remove a structural subgraph and return it as an independent graph.
    #[pyo3(name = "extract", signature = (subgraph))]
    fn py_extract(slf: Py<PyGraph>, py: Python<'_>, subgraph: &PySubgraph) -> PyResult<PyGraph> {
        Self::extract(slf, py, subgraph)
    }

    /// Delete a structural subgraph.
    #[pyo3(name = "delete", signature = (subgraph))]
    fn py_delete(slf: Py<PyGraph>, py: Python<'_>, subgraph: &PySubgraph) -> PyResult<()> {
        Self::delete(slf, py, subgraph)
    }

    /// Contract a non-empty subgraph into one replacement node.
    #[pyo3(
        name = "contract",
        signature = (subgraph, replacement, *, name=None)
    )]
    fn py_contract(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
        replacement: PyRef<'_, PyNodeValue>,
        name: Option<String>,
    ) -> PyResult<()> {
        Self::contract(slf, py, subgraph, replacement, name)
    }

    /// Append another graph without matching dangling half-edges.
    #[pyo3(name = "append", signature = (other))]
    fn py_append(slf: Py<PyGraph>, py: Python<'_>, other: Py<PyGraph>) -> PyResult<PyGraph> {
        Self::append(slf, py, other)
    }

    /// Append another graph in place without matching dangling half-edges.
    #[pyo3(name = "append_mut", signature = (other))]
    fn py_append_mut(slf: Py<PyGraph>, py: Python<'_>, other: Py<PyGraph>) -> PyResult<()> {
        Self::append_mut(slf, py, other)
    }

    /// Join dangling half-edges selected by Python callbacks.
    #[pyo3(name = "join", signature = (other, *, matching, merge))]
    fn py_join(
        slf: Py<PyGraph>,
        py: Python<'_>,
        other: Py<PyGraph>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdgeView, HalfEdgeView], builtins.bool]", imports=("builtins", "typing")))]
        matching: Py<PyAny>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdgeView, HalfEdgeView], tuple[Flow, Orientation, builtins.str | None, EdgeValue]]", imports=("builtins", "typing")))]
        merge: Py<PyAny>,
    ) -> PyResult<PyGraph> {
        Self::join(slf, py, other, matching, merge)
    }

    /// Join another graph in place while leaving the other graph unchanged.
    #[pyo3(name = "join_mut", signature = (other, *, matching, merge))]
    fn py_join_mut(
        slf: Py<PyGraph>,
        py: Python<'_>,
        other: Py<PyGraph>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdgeView, HalfEdgeView], builtins.bool]", imports=("builtins", "typing")))]
        matching: Py<PyAny>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdgeView, HalfEdgeView], tuple[Flow, Orientation, builtins.str | None, EdgeValue]]", imports=("builtins", "typing")))]
        merge: Py<PyAny>,
    ) -> PyResult<()> {
        Self::join_mut(slf, py, other, matching, merge)
    }

    fn __repr__(&self) -> PyResult<String> {
        let state = self.state.borrow();
        let state = state
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        Ok(format!(
            "Graph(nodes={}, edges={}, half_edges={})",
            state.graph.n_live(),
            state.graph.n_edges(),
            state.graph.n_hedges(),
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

impl PyHalfEdgeSpec {
    fn snapshot(&self, py: Python<'_>, node: NodeIndex) -> PyResult<HedgeData<HalfEdgeRecord>> {
        Ok(HedgeData {
            node,
            is_in_subgraph: false,
            data: HalfEdgeRecord {
                data: self
                    .data
                    .as_ref()
                    .ok_or_else(|| {
                        PyReferenceError::new_err("half-edge specification has been cleared")
                    })?
                    .clone_ref(py),
                drawing: copy_drawing(
                    py,
                    self.drawing.as_ref().ok_or_else(|| {
                        PyReferenceError::new_err("half-edge specification has been cleared")
                    })?,
                )?,
            },
        })
    }

    fn record_for_build(
        &self,
        py: Python<'_>,
        node_indices: &BTreeMap<u64, NodeIndex>,
        node_names: &BTreeMap<String, NodeIndex>,
        node_count: usize,
    ) -> PyResult<HedgeData<HalfEdgeRecord>> {
        let target = self.node(py)?;
        let target = target.bind(py);
        let node = if let Ok(spec) = target.extract::<PyRef<'_, PyNodeSpec>>() {
            *node_indices
                .get(&spec.token)
                .expect("referenced specifications were collected")
        } else if target.is_instance_of::<PyNode>() {
            return Err(PyValueError::new_err(
                "Graph() cannot use a NodeView from another graph; use its name or index explicitly",
            ));
        } else if let Ok(index) = target.extract::<usize>() {
            (index < node_count)
                .then_some(NodeIndex(index))
                .ok_or_else(|| PyIndexError::new_err("node endpoint index out of range"))?
        } else if let Ok(name) = target.extract::<String>() {
            node_names
                .get(&name)
                .copied()
                .ok_or_else(|| PyKeyError::new_err(format!("unknown node {name:?}")))?
        } else {
            return Err(PyTypeError::new_err(
                "an endpoint must be a Node, integer index, or node name",
            ));
        };
        self.snapshot(py, node)
    }

    fn record_for_graph(
        &self,
        py: Python<'_>,
        graph: &Py<PyGraph>,
        epoch: u64,
    ) -> PyResult<HedgeData<HalfEdgeRecord>> {
        let target = self.node(py)?;
        let target = target.bind(py);
        let node = if let Ok(view) = target.extract::<PyRef<'_, PyNode>>() {
            view.index_for_graph(py, graph, epoch)?
        } else if let Ok(spec) = target.extract::<PyRef<'_, PyNodeSpec>>() {
            let name = spec.name.as_deref().ok_or_else(|| {
                PyTypeError::new_err(
                    "Graph.add_edge() cannot resolve an unnamed Node; use the NodeView returned by add_node()",
                )
            })?;
            graph
                .borrow(py)
                .resolve_node(name.into_pyobject(py)?.as_any())?
        } else {
            graph.borrow(py).resolve_node(target)?
        };
        self.snapshot(py, NodeIndex(node))
    }
}

impl PyEdgeSpec {
    fn snapshot(&self, py: Python<'_>) -> PyResult<EdgeRecord> {
        Ok(EdgeRecord {
            name: self.name.clone(),
            data: self
                .data
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("edge specification has been cleared"))?
                .clone_ref(py),
            drawing: copy_drawing(
                py,
                self.drawing.as_ref().ok_or_else(|| {
                    PyReferenceError::new_err("edge specification has been cleared")
                })?,
            )?,
        })
    }

    pub(crate) fn incremental_records(
        &self,
        py: Python<'_>,
        graph: &Py<PyGraph>,
        epoch: u64,
    ) -> PyResult<IncrementalEdgeRecords> {
        let edge = self.snapshot(py)?;
        let source = self
            .source
            .as_ref()
            .map(|endpoint| endpoint.borrow(py).record_for_graph(py, graph, epoch))
            .transpose()?;
        let sink = self
            .sink
            .as_ref()
            .map(|endpoint| endpoint.borrow(py).record_for_graph(py, graph, epoch))
            .transpose()?;
        Ok(IncrementalEdgeRecords {
            source,
            sink,
            edge,
            orientation: self.orientation,
        })
    }
}

impl PyGraph {
    fn from_specs(
        py: Python<'_>,
        mut nodes: Vec<Py<PyNodeSpec>>,
        edges: Vec<Py<PyEdgeSpec>>,
        name: Option<String>,
        global_data: Option<PyGlobalData>,
        codec: Option<Py<PyDotCodec>>,
        render_config: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        let mut tokens = BTreeSet::new();
        for node in &nodes {
            if !tokens.insert(node.borrow(py).token) {
                return Err(PyValueError::new_err(
                    "the same Node was supplied more than once",
                ));
            }
        }
        // A referenced Node is a reusable description, never adopted by the graph.
        // Explicit nodes retain their order; inferred nodes follow first occurrence.
        for edge in &edges {
            let edge = edge.borrow(py);
            for endpoint in [&edge.source, &edge.sink].into_iter().flatten() {
                let target = endpoint.borrow(py).node(py)?;
                if let Ok(node) = target.extract::<Py<PyNodeSpec>>(py) {
                    if tokens.insert(node.borrow(py).token) {
                        nodes.push(node);
                    }
                }
            }
        }

        let mut names = BTreeSet::new();
        let mut node_indices = BTreeMap::new();
        let mut node_names = BTreeMap::new();
        let mut builder = HedgeGraphBuilder::<EdgeRecord, NodeRecord, HalfEdgeRecord>::new();
        for spec in &nodes {
            let spec = spec.borrow(py);
            if let Some(name) = &spec.name {
                if !names.insert(name.clone()) {
                    return Err(PyValueError::new_err(format!(
                        "duplicate node name {name:?}"
                    )));
                }
            }
            let index = builder.add_node(spec.snapshot(py)?);
            if let Some(name) = &spec.name {
                node_names.insert(name.clone(), index);
            }
            node_indices.insert(spec.token, index);
        }

        names.clear();
        for spec in edges {
            let spec = spec.borrow(py);
            if let Some(name) = &spec.name {
                if !names.insert(name.clone()) {
                    return Err(PyValueError::new_err(format!(
                        "duplicate edge name {name:?}"
                    )));
                }
            }
            let source = spec
                .source
                .as_ref()
                .map(|endpoint| {
                    endpoint.borrow(py).record_for_build(
                        py,
                        &node_indices,
                        &node_names,
                        nodes.len(),
                    )
                })
                .transpose()?;
            let sink = spec
                .sink
                .as_ref()
                .map(|endpoint| {
                    endpoint.borrow(py).record_for_build(
                        py,
                        &node_indices,
                        &node_names,
                        nodes.len(),
                    )
                })
                .transpose()?;
            let record = spec.snapshot(py)?;
            match (source, sink) {
                (Some(source), Some(sink)) => {
                    builder.add_edge(source, sink, record, spec.orientation)
                }
                (Some(source), None) => {
                    builder.add_external_edge(source, record, spec.orientation, Flow::Source)
                }
                (None, Some(sink)) => {
                    builder.add_external_edge(sink, record, spec.orientation, Flow::Sink)
                }
                (None, None) => {
                    return Err(PyValueError::new_err("an edge needs at least one endpoint"))
                }
            }
        }

        Ok(Self::from_state(GraphState {
            graph: builder.build(),
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
            element_epoch: 0,
        }))
    }
}

#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        import typing

        class PyNodeSpec:
            def __new__(cls, name: _OptionalString = None, *, data: typing.Any = None, drawing: _OptionalNodeDrawing = None, label: _OptionalStaticContent = ..., placement: _PlacementValue = ..., shift: _DrawingPoint = ..., rank: _OptionalInteger = ..., minimum_size: _OptionalNumber = ..., maximum_size: _OptionalNumber = ..., style: _OptionalStyle = ..., label_style: _OptionalStyle = ..., extensions: _NativeDict = ...) -> Node: ...
    "# }
}

#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        import typing

        class PyHalfEdgeSpec:
            def __new__(cls, node: _EndpointTarget, *, data: typing.Any = None, drawing: _OptionalHalfEdgeDrawing = None, label: _OptionalStaticContent = ..., statement: _DrawingString = ..., port_label: _DrawingString = ..., compass: _CompassValue = ..., route_points: _DrawingPoints = ..., anchor: _AnchorValue = ..., routing: _RoutingValue = ..., style: _OptionalStyleLayers = ..., extensions: _NativeDict = ...) -> HalfEdge: ...
    "# }
}

#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        import typing

        class PyEdgeSpec:
            def __new__(cls, source: _EdgeEndpoint = None, sink: _EdgeEndpoint = None, *, name: _OptionalString = None, data: typing.Any = None, orientation: Orientation = Orientation.Default, drawing: _OptionalEdgeDrawing = None, label: _OptionalStaticContent = ..., placement: _PlacementValue = ..., label_position: _DrawingPoint = ..., label_offset: _OptionalNumber = ..., label_angle: _DrawingAngle = ..., bend: _DrawingAngle = ..., routing: _RoutingValue = ..., minimum_length: _OptionalInteger = ..., same_rank: _OptionalBoolean = ..., style: _OptionalStyleLayers = ..., label_style: _OptionalStyle = ..., decoration: _DrawingDecoration = ..., extensions: _NativeDict = ...) -> Edge: ...
    "# }
}

#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        import typing

        class PyGraph:
            def add_node(self, node: _NodeInput = None, *, name: _OptionalString = None, data: typing.Any = None, drawing: _OptionalNodeDrawing = None, label: _OptionalStaticContent = ..., placement: _PlacementValue = ..., shift: _DrawingPoint = ..., rank: _OptionalInteger = ..., minimum_size: _OptionalNumber = ..., maximum_size: _OptionalNumber = ..., style: _OptionalStyle = ..., label_style: _OptionalStyle = ..., extensions: _NativeDict = ...) -> NodeView: ...
            def add_edge(self, source: _EdgeInput = None, sink: _EdgeEndpoint = None, *, name: _OptionalString = None, data: typing.Any = None, orientation: _OptionalOrientation = None, drawing: _OptionalEdgeDrawing = None, label: _OptionalStaticContent = ..., placement: _PlacementValue = ..., label_position: _DrawingPoint = ..., label_offset: _OptionalNumber = ..., label_angle: _DrawingAngle = ..., bend: _DrawingAngle = ..., routing: _RoutingValue = ..., minimum_length: _OptionalInteger = ..., same_rank: _OptionalBoolean = ..., style: _OptionalStyleLayers = ..., label_style: _OptionalStyle = ..., decoration: _DrawingDecoration = ..., extensions: _NativeDict = ...) -> EdgeView: ...
    "# }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyFlow>()?;
    module.add_class::<PyOrientation>()?;
    module.add_class::<PyNodeSpec>()?;
    module.add_class::<PyHalfEdgeSpec>()?;
    module.add_class::<PyEdgeSpec>()?;
    module.add_class::<PyGraph>()?;
    generate::register(module)?;
    module.add_class::<PyNode>()?;
    module.add_class::<PyEdge>()?;
    module.add_class::<PyHalfEdge>()?;
    Ok(())
}
