use std::cell::RefCell;
use std::collections::BTreeMap;

use linnet::half_edge::builder::{HedgeData, HedgeGraphBuilder};
use linnet::half_edge::involution::{EdgeData, EdgeIndex, Flow, Hedge, HedgePair, Orientation};
use linnet::half_edge::nodestore::DefaultNodeStore;
use linnet::half_edge::subgraph::{
    cut::OrientedCut, cycle::Cycle, internal::InternalSubGraph, node::HedgeNode, Inclusion,
    ModifySubSet, SuBitGraph, SubSetLike, SubSetOps,
};
use linnet::half_edge::tree::SimpleTraversalTree;
use linnet::half_edge::{HedgeGraphError, NodeIndex};
use linnet::parser::{
    set::DotGraphSet, DotEdgeData, DotGraph, DotHedgeData, DotVertexData, GlobalData,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyList, PyString, PyTuple, PyType};
use pyo3_stub_gen::{
    define_stub_info_gatherer,
    derive::{gen_stub_pyclass, gen_stub_pymethods},
    inventory::submit,
    type_info::{
        MethodInfo, MethodType, ParameterDefault, ParameterInfo, ParameterKind, PyMethodsInfo,
    },
    PyStubType,
};

/// Half-edge identifier.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "Hedge")]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct PyHedge {
    hedge: Hedge,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyHedge {
    /// Create a hedge from a zero-based index.
    #[new]
    fn new(value: usize) -> Self {
        Self {
            hedge: Hedge(value),
        }
    }

    /// Numeric value of this hedge.
    #[getter]
    fn value(&self) -> usize {
        self.hedge.0
    }

    /// Convert to int.
    fn __int__(&self) -> usize {
        self.hedge.0
    }

    /// Debug-style representation.
    fn __repr__(&self) -> String {
        format!("Hedge({})", self.hedge.0)
    }
}

/// Node index identifier.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "NodeIndex")]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct PyNodeIndex {
    node: NodeIndex,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyNodeIndex {
    /// Create a node index from a zero-based index.
    #[new]
    fn new(value: usize) -> Self {
        Self {
            node: NodeIndex(value),
        }
    }

    /// Numeric value of this node index.
    #[getter]
    fn value(&self) -> usize {
        self.node.0
    }

    /// Convert to int.
    fn __int__(&self) -> usize {
        self.node.0
    }

    /// Debug-style representation.
    fn __repr__(&self) -> String {
        format!("NodeIndex({})", self.node.0)
    }
}

/// Edge index identifier.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "EdgeIndex")]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct PyEdgeIndex {
    edge: EdgeIndex,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyEdgeIndex {
    /// Create an edge index from a zero-based index.
    #[new]
    fn new(value: usize) -> Self {
        Self {
            edge: EdgeIndex::from(value),
        }
    }

    /// Numeric value of this edge index.
    #[getter]
    fn value(&self) -> usize {
        self.edge.0
    }

    /// Convert to int.
    fn __int__(&self) -> usize {
        self.edge.0
    }

    /// Debug-style representation.
    fn __repr__(&self) -> String {
        format!("EdgeIndex({})", self.edge.0)
    }
}

/// Flow direction (Source/Sink).
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "Flow")]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct PyFlow {
    flow: Flow,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyFlow {
    /// Flow::Source.
    #[staticmethod]
    fn source() -> Self {
        Self { flow: Flow::Source }
    }

    /// Flow::Sink.
    #[staticmethod]
    fn sink() -> Self {
        Self { flow: Flow::Sink }
    }

    /// Debug-style representation.
    fn __repr__(&self) -> String {
        match self.flow {
            Flow::Source => "Flow.Source".to_string(),
            Flow::Sink => "Flow.Sink".to_string(),
        }
    }
}

/// Edge orientation (Default/Reversed/Undirected).
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "Orientation")]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct PyOrientation {
    orientation: Orientation,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyOrientation {
    /// Orientation::Default.
    #[staticmethod]
    fn default() -> Self {
        Self {
            orientation: Orientation::Default,
        }
    }

    /// Orientation::Reversed.
    #[staticmethod]
    fn reversed() -> Self {
        Self {
            orientation: Orientation::Reversed,
        }
    }

    /// Orientation::Undirected.
    #[staticmethod]
    fn undirected() -> Self {
        Self {
            orientation: Orientation::Undirected,
        }
    }

    /// Debug-style representation.
    fn __repr__(&self) -> String {
        match self.orientation {
            Orientation::Default => "Orientation.Default".to_string(),
            Orientation::Reversed => "Orientation.Reversed".to_string(),
            Orientation::Undirected => "Orientation.Undirected".to_string(),
        }
    }
}

#[derive(Debug)]
enum DotVertexDataSource {
    Owned(DotVertexData),
    Graph { graph: Py<PyAny>, node: NodeIndex },
}

/// DOT vertex data (name, optional index, and attribute statements).
#[gen_stub_pyclass]
#[pyclass(name = "DotVertexData")]
#[derive(Debug)]
struct PyDotVertexData {
    data: DotVertexDataSource,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyDotVertexData {
    /// Create vertex data from fields.
    #[new]
    #[pyo3(signature = (name=None, index=None, statements=None))]
    fn new(
        name: Option<String>,
        index: Option<usize>,
        statements: Option<BTreeMap<String, String>>,
    ) -> Self {
        Self {
            data: DotVertexDataSource::Owned(DotVertexData {
                name,
                index: index.map(NodeIndex),
                statements: statements.unwrap_or_default(),
            }),
        }
    }

    /// Optional vertex name.
    #[getter]
    fn name(&self, py: Python<'_>) -> Option<String> {
        self.read_data(py, |d| d.name.clone())
    }

    /// Optional node index.
    #[getter]
    fn index(&self, py: Python<'_>) -> Option<PyNodeIndex> {
        self.read_data(py, |d| d.index.map(|i| PyNodeIndex { node: i }))
    }

    /// Attribute statements as a dict-like proxy.
    #[getter]
    fn statements(slf: PyRef<'_, Self>) -> PyResult<PyNodeStatements> {
        let py = slf.py();
        let data = slf.into_pyobject(py)?.unbind().into();
        Ok(PyNodeStatements { data })
    }

    /// Insert or update an attribute statement.
    fn add_statement(&mut self, key: String, value: String) {
        let _ = Python::try_attach(|py| self.insert_statement(key, value, py));
    }

    /// Debug-style representation.
    fn __repr__(&self, py: Python<'_>) -> String {
        let (name, index, count) = self.read_data(py, |d| {
            (d.name.clone(), d.index.map(|i| i.0), d.statements.len())
        });
        format!(
            "DotVertexData(name={:?}, index={:?}, statements={})",
            name, index, count
        )
    }
}

/// Dict-like proxy for node attribute statements.
#[gen_stub_pyclass]
#[pyclass(name = "NodeStatements")]
#[derive(Debug)]
struct PyNodeStatements {
    data: Py<PyAny>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyNodeStatements {
    fn __getitem__(&self, key: String, py: Python<'_>) -> PyResult<String> {
        let value = self.get_value(py, &key);
        value.ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(key))
    }

    fn __setitem__(&mut self, key: String, value: String, py: Python<'_>) -> PyResult<()> {
        let mut data = self.data.bind(py).extract::<PyRefMut<PyDotVertexData>>()?;
        data.insert_statement(key, value, py);
        Ok(())
    }

    fn __delitem__(&mut self, key: String, py: Python<'_>) -> PyResult<()> {
        let mut data = self.data.bind(py).extract::<PyRefMut<PyDotVertexData>>()?;
        if data.remove_statement(&key, py) {
            Ok(())
        } else {
            Err(pyo3::exceptions::PyKeyError::new_err(key))
        }
    }

    fn __contains__(&self, key: String, py: Python<'_>) -> bool {
        self.get_value(py, &key).is_some()
    }

    fn __len__(&self, py: Python<'_>) -> usize {
        self.snapshot(py).len()
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let keys = self.keys(py)?;
        keys.call_method0("__iter__")
    }

    fn get<'py>(
        &self,
        key: String,
        default: Option<Bound<'py, PyAny>>,
        py: Python<'py>,
    ) -> Bound<'py, PyAny> {
        if let Some(value) = self.get_value(py, &key) {
            PyString::new(py, &value).into_any()
        } else {
            default.unwrap_or_else(|| py.None().into_bound(py))
        }
    }

    fn keys<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let keys = self
            .snapshot(py)
            .into_iter()
            .map(|(k, _)| k)
            .collect::<Vec<_>>();
        Ok(PyList::new(py, keys)?.into_any())
    }

    fn values<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let values = self
            .snapshot(py)
            .into_iter()
            .map(|(_, v)| v)
            .collect::<Vec<_>>();
        Ok(PyList::new(py, values)?.into_any())
    }

    fn items<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let items = self.snapshot(py);
        Ok(PyList::new(py, items)?.into_any())
    }

    fn clear(&mut self, py: Python<'_>) {
        if let Ok(mut data) = self.data.bind(py).extract::<PyRefMut<PyDotVertexData>>() {
            data.clear_statements(py);
        }
    }

    fn update(&mut self, other: &Bound<'_, PyAny>, py: Python<'_>) -> PyResult<()> {
        let map = other.extract::<BTreeMap<String, String>>()?;
        let mut data = self.data.bind(py).extract::<PyRefMut<PyDotVertexData>>()?;
        for (k, v) in map {
            data.insert_statement(k, v, py);
        }
        Ok(())
    }

    #[pyo3(signature = (key, default=None))]
    fn pop<'py>(
        &mut self,
        key: String,
        default: Option<Bound<'py, PyAny>>,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let value = self.get_value(py, &key);
        if let Some(value) = value {
            let mut data = self.data.bind(py).extract::<PyRefMut<PyDotVertexData>>()?;
            data.remove_statement(&key, py);
            Ok(PyString::new(py, &value).into_any())
        } else if let Some(default) = default {
            Ok(default)
        } else {
            Err(pyo3::exceptions::PyKeyError::new_err(key))
        }
    }

    fn popitem<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let items = self.snapshot(py);
        if let Some((key, value)) = items.into_iter().next() {
            let mut data = self.data.bind(py).extract::<PyRefMut<PyDotVertexData>>()?;
            data.remove_statement(&key, py);
            Ok(PyTuple::new(py, [key, value])?.into_any())
        } else {
            Err(pyo3::exceptions::PyKeyError::new_err("popitem(): empty"))
        }
    }

    #[pyo3(signature = (key, default=None))]
    fn setdefault<'py>(
        &mut self,
        key: String,
        default: Option<String>,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if let Some(value) = self.get_value(py, &key) {
            Ok(PyString::new(py, &value).into_any())
        } else {
            let value = default.unwrap_or_default();
            let mut data = self.data.bind(py).extract::<PyRefMut<PyDotVertexData>>()?;
            data.insert_statement(key, value.clone(), py);
            Ok(PyString::new(py, &value).into_any())
        }
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!("NodeStatements(len={})", self.snapshot(py).len())
    }
}

impl PyDotVertexData {
    fn read_data<R>(&self, py: Python<'_>, f: impl FnOnce(&DotVertexData) -> R) -> R {
        match &self.data {
            DotVertexDataSource::Owned(data) => f(data),
            DotVertexDataSource::Graph { graph, node } => {
                let graph_ref = graph
                    .bind(py)
                    .extract::<PyRef<PyDotGraph>>()
                    .expect("DotVertexData: invalid graph reference");
                f(&graph_ref.graph.graph[*node])
            }
        }
    }

    fn write_data(&mut self, py: Python<'_>, f: impl FnOnce(&mut DotVertexData)) {
        match &mut self.data {
            DotVertexDataSource::Owned(data) => f(data),
            DotVertexDataSource::Graph { graph, node } => {
                let mut graph_ref = graph
                    .bind(py)
                    .extract::<PyRefMut<PyDotGraph>>()
                    .expect("DotVertexData: invalid graph reference");
                f(&mut graph_ref.graph.graph[*node])
            }
        }
    }

    fn snapshot(&self, py: Python<'_>) -> DotVertexData {
        self.read_data(py, |d| d.clone())
    }

    fn get_statement(&self, key: &str, py: Python<'_>) -> Option<String> {
        self.read_data(py, |d| d.statements.get(key).cloned())
    }

    fn snapshot_statements(&self, py: Python<'_>) -> Vec<(String, String)> {
        self.read_data(py, |d| {
            d.statements
                .iter()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect()
        })
    }

    fn insert_statement(&mut self, key: String, value: String, py: Python<'_>) {
        self.write_data(py, |d| {
            d.statements.insert(key, value);
        });
    }

    fn remove_statement(&mut self, key: &str, py: Python<'_>) -> bool {
        let mut removed = false;
        self.write_data(py, |d| {
            removed = d.statements.remove(key).is_some();
        });
        removed
    }

    fn clear_statements(&mut self, py: Python<'_>) {
        self.write_data(py, |d| d.statements.clear());
    }
}

impl PyNodeStatements {
    fn get_value(&self, py: Python<'_>, key: &str) -> Option<String> {
        let data = self
            .data
            .bind(py)
            .extract::<PyRef<PyDotVertexData>>()
            .expect("NodeStatements: invalid data reference");
        data.get_statement(key, py)
    }

    fn snapshot(&self, py: Python<'_>) -> Vec<(String, String)> {
        let data = self
            .data
            .bind(py)
            .extract::<PyRef<PyDotVertexData>>()
            .expect("NodeStatements: invalid data reference");
        data.snapshot_statements(py)
    }
}

/// DOT hedge (half-edge) data.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "DotHedgeData")]
#[derive(Clone, Debug, PartialEq, Eq)]
struct PyDotHedgeData {
    data: DotHedgeData,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyDotHedgeData {
    /// Create hedge data from fields.
    #[new]
    #[pyo3(signature = (statement=None, id=None, port_label=None, compasspt=None))]
    fn new(
        statement: Option<String>,
        id: Option<usize>,
        port_label: Option<String>,
        compasspt: Option<String>,
    ) -> Self {
        let compasspt = compasspt.and_then(|value| parse_compasspt(&value));
        Self {
            data: DotHedgeData {
                statement,
                id: id.map(Hedge),
                port_label,
                compasspt,
            },
        }
    }

    /// Optional statement string.
    #[getter]
    fn statement(&self) -> Option<String> {
        self.data.statement.clone()
    }

    /// Optional hedge id.
    #[getter]
    fn id(&self) -> Option<PyHedge> {
        self.data.id.map(|h| PyHedge { hedge: h })
    }

    /// Optional port label.
    #[getter]
    fn port_label(&self) -> Option<String> {
        self.data.port_label.clone()
    }

    /// Optional compass point as a string.
    #[getter]
    fn compasspt(&self) -> Option<String> {
        self.data.compasspt.as_ref().map(|c| format!("{c:?}"))
    }

    /// Debug-style representation.
    fn __repr__(&self) -> String {
        format!(
            "DotHedgeData(statement={:?}, id={:?}, port_label={:?}, compasspt={:?})",
            self.data.statement,
            self.data.id.map(|h| h.0),
            self.data.port_label,
            self.data.compasspt.as_ref().map(|c| format!("{c:?}"))
        )
    }
}

#[derive(Debug)]
enum DotEdgeDataSource {
    Owned(DotEdgeData),
    Graph { graph: Py<PyAny>, edge: EdgeIndex },
}

/// DOT edge data (attributes + optional edge id).
#[gen_stub_pyclass]
#[pyclass(name = "DotEdgeData")]
#[derive(Debug)]
struct PyDotEdgeData {
    data: DotEdgeDataSource,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyDotEdgeData {
    /// Create edge data from fields.
    #[new]
    #[pyo3(signature = (statements=None, local_statements=None, edge_id=None))]
    fn new(
        statements: Option<BTreeMap<String, String>>,
        local_statements: Option<BTreeMap<String, String>>,
        edge_id: Option<usize>,
    ) -> Self {
        Self {
            data: DotEdgeDataSource::Owned(DotEdgeData {
                statements: statements.unwrap_or_default(),
                local_statements: local_statements.unwrap_or_default(),
                edge_id: edge_id.map(EdgeIndex::from),
            }),
        }
    }

    /// Edge attributes as a dict-like proxy.
    #[getter]
    fn statements(slf: PyRef<'_, Self>) -> PyResult<PyEdgeStatements> {
        let py = slf.py();
        let data = slf.into_pyobject(py)?.unbind().into();
        Ok(PyEdgeStatements { data, local: false })
    }

    /// Local edge attributes as a dict-like proxy.
    #[getter]
    fn local_statements(slf: PyRef<'_, Self>) -> PyResult<PyEdgeStatements> {
        let py = slf.py();
        let data = slf.into_pyobject(py)?.unbind().into();
        Ok(PyEdgeStatements { data, local: true })
    }

    /// Optional edge id.
    #[getter]
    fn edge_id(&self, py: Python<'_>) -> Option<PyEdgeIndex> {
        self.read_data(py, |d| d.edge_id.map(|i| PyEdgeIndex { edge: i }))
    }

    /// Insert or update an attribute statement.
    fn add_statement(&mut self, key: String, value: String) {
        let _ = Python::try_attach(|py| self.insert_statement(false, key, value, py));
    }

    /// Debug-style representation.
    fn __repr__(&self, py: Python<'_>) -> String {
        let (statements, local, edge_id) = self.read_data(py, |d| {
            (
                d.statements.len(),
                d.local_statements.len(),
                d.edge_id.map(|i| i.0),
            )
        });
        format!(
            "DotEdgeData(statements={}, local_statements={}, edge_id={:?})",
            statements, local, edge_id
        )
    }
}

/// Dict-like proxy for edge attribute statements.
#[gen_stub_pyclass]
#[pyclass(name = "EdgeStatements")]
#[derive(Debug)]
struct PyEdgeStatements {
    data: Py<PyAny>,
    local: bool,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyEdgeStatements {
    fn __getitem__(&self, key: String, py: Python<'_>) -> PyResult<String> {
        let value = self.get_value(py, &key);
        value.ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(key))
    }

    fn __setitem__(&mut self, key: String, value: String, py: Python<'_>) -> PyResult<()> {
        let mut data = self.data.bind(py).extract::<PyRefMut<PyDotEdgeData>>()?;
        data.insert_statement(self.local, key, value, py);
        Ok(())
    }

    fn __delitem__(&mut self, key: String, py: Python<'_>) -> PyResult<()> {
        let mut data = self.data.bind(py).extract::<PyRefMut<PyDotEdgeData>>()?;
        if data.remove_statement(self.local, &key, py) {
            Ok(())
        } else {
            Err(pyo3::exceptions::PyKeyError::new_err(key))
        }
    }

    fn __contains__(&self, key: String, py: Python<'_>) -> bool {
        self.get_value(py, &key).is_some()
    }

    fn __len__(&self, py: Python<'_>) -> usize {
        self.snapshot(py).len()
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let keys = self.keys(py)?;
        keys.call_method0("__iter__")
    }

    fn get<'py>(
        &self,
        key: String,
        default: Option<Bound<'py, PyAny>>,
        py: Python<'py>,
    ) -> Bound<'py, PyAny> {
        if let Some(value) = self.get_value(py, &key) {
            PyString::new(py, &value).into_any()
        } else {
            default.unwrap_or_else(|| py.None().into_bound(py))
        }
    }

    fn keys<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let keys = self
            .snapshot(py)
            .into_iter()
            .map(|(k, _)| k)
            .collect::<Vec<_>>();
        Ok(PyList::new(py, keys)?.into_any())
    }

    fn values<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let values = self
            .snapshot(py)
            .into_iter()
            .map(|(_, v)| v)
            .collect::<Vec<_>>();
        Ok(PyList::new(py, values)?.into_any())
    }

    fn items<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let items = self.snapshot(py);
        Ok(PyList::new(py, items)?.into_any())
    }

    fn clear(&mut self, py: Python<'_>) {
        if let Ok(mut data) = self.data.bind(py).extract::<PyRefMut<PyDotEdgeData>>() {
            data.clear_statements(self.local, py);
        }
    }

    fn update(&mut self, other: &Bound<'_, PyAny>, py: Python<'_>) -> PyResult<()> {
        let map = other.extract::<BTreeMap<String, String>>()?;
        let mut data = self.data.bind(py).extract::<PyRefMut<PyDotEdgeData>>()?;
        for (k, v) in map {
            data.insert_statement(self.local, k, v, py);
        }
        Ok(())
    }

    #[pyo3(signature = (key, default=None))]
    fn pop<'py>(
        &mut self,
        key: String,
        default: Option<Bound<'py, PyAny>>,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let value = self.get_value(py, &key);
        if let Some(value) = value {
            let mut data = self.data.bind(py).extract::<PyRefMut<PyDotEdgeData>>()?;
            data.remove_statement(self.local, &key, py);
            Ok(PyString::new(py, &value).into_any())
        } else if let Some(default) = default {
            Ok(default)
        } else {
            Err(pyo3::exceptions::PyKeyError::new_err(key))
        }
    }

    fn popitem<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let items = self.snapshot(py);
        if let Some((key, value)) = items.into_iter().next() {
            let mut data = self.data.bind(py).extract::<PyRefMut<PyDotEdgeData>>()?;
            data.remove_statement(self.local, &key, py);
            Ok(PyTuple::new(py, [key, value])?.into_any())
        } else {
            Err(pyo3::exceptions::PyKeyError::new_err("popitem(): empty"))
        }
    }

    #[pyo3(signature = (key, default=None))]
    fn setdefault<'py>(
        &mut self,
        key: String,
        default: Option<String>,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if let Some(value) = self.get_value(py, &key) {
            Ok(PyString::new(py, &value).into_any())
        } else {
            let value = default.unwrap_or_default();
            let mut data = self.data.bind(py).extract::<PyRefMut<PyDotEdgeData>>()?;
            data.insert_statement(self.local, key, value.clone(), py);
            Ok(PyString::new(py, &value).into_any())
        }
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        let label = if self.local {
            "EdgeStatements(local)"
        } else {
            "EdgeStatements"
        };
        format!("{}(len={})", label, self.snapshot(py).len())
    }
}

impl PyDotEdgeData {
    fn read_data<R>(&self, py: Python<'_>, f: impl FnOnce(&DotEdgeData) -> R) -> R {
        match &self.data {
            DotEdgeDataSource::Owned(data) => f(data),
            DotEdgeDataSource::Graph { graph, edge } => {
                let graph_ref = graph
                    .bind(py)
                    .extract::<PyRef<PyDotGraph>>()
                    .expect("DotEdgeData: invalid graph reference");
                f(&graph_ref.graph.graph[*edge])
            }
        }
    }

    fn write_data(&mut self, py: Python<'_>, f: impl FnOnce(&mut DotEdgeData)) {
        match &mut self.data {
            DotEdgeDataSource::Owned(data) => f(data),
            DotEdgeDataSource::Graph { graph, edge } => {
                let mut graph_ref = graph
                    .bind(py)
                    .extract::<PyRefMut<PyDotGraph>>()
                    .expect("DotEdgeData: invalid graph reference");
                f(&mut graph_ref.graph.graph[*edge])
            }
        }
    }

    fn snapshot(&self, py: Python<'_>) -> DotEdgeData {
        self.read_data(py, |d| d.clone())
    }

    fn get_statement(&self, local: bool, key: &str, py: Python<'_>) -> Option<String> {
        self.read_data(py, |d| {
            let map = if local {
                &d.local_statements
            } else {
                &d.statements
            };
            map.get(key).cloned()
        })
    }

    fn snapshot_statements(&self, local: bool, py: Python<'_>) -> Vec<(String, String)> {
        self.read_data(py, |d| {
            let map = if local {
                &d.local_statements
            } else {
                &d.statements
            };
            map.iter().map(|(k, v)| (k.clone(), v.clone())).collect()
        })
    }

    fn insert_statement(&mut self, local: bool, key: String, value: String, py: Python<'_>) {
        self.write_data(py, |d| {
            let map = if local {
                &mut d.local_statements
            } else {
                &mut d.statements
            };
            map.insert(key, value);
        });
    }

    fn remove_statement(&mut self, local: bool, key: &str, py: Python<'_>) -> bool {
        let mut removed = false;
        self.write_data(py, |d| {
            let map = if local {
                &mut d.local_statements
            } else {
                &mut d.statements
            };
            removed = map.remove(key).is_some();
        });
        removed
    }

    fn clear_statements(&mut self, local: bool, py: Python<'_>) {
        self.write_data(py, |d| {
            let map = if local {
                &mut d.local_statements
            } else {
                &mut d.statements
            };
            map.clear();
        });
    }
}

impl PyEdgeStatements {
    fn get_value(&self, py: Python<'_>, key: &str) -> Option<String> {
        let data = self
            .data
            .bind(py)
            .extract::<PyRef<PyDotEdgeData>>()
            .expect("EdgeStatements: invalid data reference");
        data.get_statement(self.local, key, py)
    }

    fn snapshot(&self, py: Python<'_>) -> Vec<(String, String)> {
        let data = self
            .data
            .bind(py)
            .extract::<PyRef<PyDotEdgeData>>()
            .expect("EdgeStatements: invalid data reference");
        data.snapshot_statements(self.local, py)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum StatementKind {
    Graph,
    Edge,
    Node,
}

/// Dict-like proxy for global statements.
#[gen_stub_pyclass]
#[pyclass(name = "GlobalStatements")]
#[derive(Debug)]
struct PyGlobalStatements {
    data: Py<PyAny>,
    kind: StatementKind,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyGlobalStatements {
    fn __getitem__(&self, key: String, py: Python<'_>) -> PyResult<String> {
        let value = self.get_value(py, &key);
        value.ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(key))
    }

    fn __setitem__(&mut self, key: String, value: String, py: Python<'_>) -> PyResult<()> {
        let mut data = self.data.bind(py).extract::<PyRefMut<PyGlobalData>>()?;
        data.insert_statement(self.kind, key, value, py);
        Ok(())
    }

    fn __delitem__(&mut self, key: String, py: Python<'_>) -> PyResult<()> {
        let mut data = self.data.bind(py).extract::<PyRefMut<PyGlobalData>>()?;
        if data.remove_statement(self.kind, &key, py) {
            Ok(())
        } else {
            Err(pyo3::exceptions::PyKeyError::new_err(key))
        }
    }

    fn __contains__(&self, key: String, py: Python<'_>) -> bool {
        self.get_value(py, &key).is_some()
    }

    fn __len__(&self, py: Python<'_>) -> usize {
        self.snapshot(py).len()
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let keys = self.keys(py)?;
        keys.call_method0("__iter__")
    }

    fn get<'py>(
        &self,
        key: String,
        default: Option<Bound<'py, PyAny>>,
        py: Python<'py>,
    ) -> Bound<'py, PyAny> {
        if let Some(value) = self.get_value(py, &key) {
            PyString::new(py, &value).into_any()
        } else {
            default.unwrap_or_else(|| py.None().into_bound(py))
        }
    }

    fn keys<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let keys = self
            .snapshot(py)
            .into_iter()
            .map(|(k, _)| k)
            .collect::<Vec<_>>();
        Ok(PyList::new(py, keys)?.into_any())
    }

    fn values<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let values = self
            .snapshot(py)
            .into_iter()
            .map(|(_, v)| v)
            .collect::<Vec<_>>();
        Ok(PyList::new(py, values)?.into_any())
    }

    fn items<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let items = self.snapshot(py);
        Ok(PyList::new(py, items)?.into_any())
    }

    fn clear(&mut self, py: Python<'_>) {
        if let Ok(mut data) = self.data.bind(py).extract::<PyRefMut<PyGlobalData>>() {
            data.clear_statements(self.kind, py);
        }
    }

    fn update(&mut self, other: &Bound<'_, PyAny>, py: Python<'_>) -> PyResult<()> {
        let map = other.extract::<BTreeMap<String, String>>()?;
        let mut data = self.data.bind(py).extract::<PyRefMut<PyGlobalData>>()?;
        for (k, v) in map {
            data.insert_statement(self.kind, k, v, py);
        }
        Ok(())
    }

    #[pyo3(signature = (key, default=None))]
    fn pop<'py>(
        &mut self,
        key: String,
        default: Option<Bound<'py, PyAny>>,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let value = self.get_value(py, &key);
        if let Some(value) = value {
            let mut data = self.data.bind(py).extract::<PyRefMut<PyGlobalData>>()?;
            data.remove_statement(self.kind, &key, py);
            Ok(PyString::new(py, &value).into_any())
        } else if let Some(default) = default {
            Ok(default)
        } else {
            Err(pyo3::exceptions::PyKeyError::new_err(key))
        }
    }

    fn popitem<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let items = self.snapshot(py);
        if let Some((key, value)) = items.into_iter().next() {
            let mut data = self.data.bind(py).extract::<PyRefMut<PyGlobalData>>()?;
            data.remove_statement(self.kind, &key, py);
            Ok(PyTuple::new(py, [key, value])?.into_any())
        } else {
            Err(pyo3::exceptions::PyKeyError::new_err("popitem(): empty"))
        }
    }

    #[pyo3(signature = (key, default=None))]
    fn setdefault<'py>(
        &mut self,
        key: String,
        default: Option<String>,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if let Some(value) = self.get_value(py, &key) {
            Ok(PyString::new(py, &value).into_any())
        } else {
            let value = default.unwrap_or_default();
            let mut data = self.data.bind(py).extract::<PyRefMut<PyGlobalData>>()?;
            data.insert_statement(self.kind, key, value.clone(), py);
            Ok(PyString::new(py, &value).into_any())
        }
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!("GlobalStatements(len={})", self.snapshot(py).len())
    }
}

/// Global graph attributes (graph/node/edge statements).
#[gen_stub_pyclass]
#[pyclass(name = "GlobalData")]
#[derive(Debug)]
struct PyGlobalData {
    data: GlobalData,
    #[gen_stub(skip)]
    graph: Option<Py<PyAny>>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyGlobalData {
    /// Create global data from fields.
    #[new]
    #[pyo3(signature = (name=None, statements=None, edge_statements=None, node_statements=None))]
    fn new(
        name: Option<String>,
        statements: Option<BTreeMap<String, String>>,
        edge_statements: Option<BTreeMap<String, String>>,
        node_statements: Option<BTreeMap<String, String>>,
    ) -> Self {
        Self {
            data: GlobalData {
                name: name.unwrap_or_default(),
                statements: statements.unwrap_or_default(),
                edge_statements: edge_statements.unwrap_or_default(),
                node_statements: node_statements.unwrap_or_default(),
            },
            graph: None,
        }
    }

    /// Graph name.
    #[getter]
    fn name(&self, py: Python<'_>) -> String {
        self.read_global(|g| g.name.clone(), py)
    }

    #[setter]
    fn set_name(&mut self, name: String) {
        let _ = Python::try_attach(|py| self.write_global(|g| g.name = name, py));
    }

    /// Graph-level statements.
    #[getter]
    fn statements(slf: PyRef<'_, Self>) -> PyResult<PyGlobalStatements> {
        let py = slf.py();
        let data = slf.into_pyobject(py)?.unbind().into();
        Ok(PyGlobalStatements {
            data,
            kind: StatementKind::Graph,
        })
    }

    #[setter]
    fn set_statements(&mut self, statements: BTreeMap<String, String>) {
        let _ = Python::try_attach(|py| self.write_global(|g| g.statements = statements, py));
    }

    /// Default edge statements.
    #[getter]
    fn edge_statements(slf: PyRef<'_, Self>) -> PyResult<PyGlobalStatements> {
        let py = slf.py();
        let data = slf.into_pyobject(py)?.unbind().into();
        Ok(PyGlobalStatements {
            data,
            kind: StatementKind::Edge,
        })
    }

    #[setter]
    fn set_edge_statements(&mut self, statements: BTreeMap<String, String>) {
        let _ = Python::try_attach(|py| self.write_global(|g| g.edge_statements = statements, py));
    }

    /// Default node statements.
    #[getter]
    fn node_statements(slf: PyRef<'_, Self>) -> PyResult<PyGlobalStatements> {
        let py = slf.py();
        let data = slf.into_pyobject(py)?.unbind().into();
        Ok(PyGlobalStatements {
            data,
            kind: StatementKind::Node,
        })
    }

    #[setter]
    fn set_node_statements(&mut self, statements: BTreeMap<String, String>) {
        let _ = Python::try_attach(|py| self.write_global(|g| g.node_statements = statements, py));
    }

    /// Insert or update a graph-level statement.
    fn add_statement(&mut self, key: String, value: String) {
        let _ =
            Python::try_attach(|py| self.insert_statement(StatementKind::Graph, key, value, py));
    }

    /// Insert or update a default edge statement.
    fn add_edge_statement(&mut self, key: String, value: String) {
        let _ = Python::try_attach(|py| self.insert_statement(StatementKind::Edge, key, value, py));
    }

    /// Insert or update a default node statement.
    fn add_node_statement(&mut self, key: String, value: String) {
        let _ = Python::try_attach(|py| self.insert_statement(StatementKind::Node, key, value, py));
    }

    /// Debug-style representation.
    fn __repr__(&self, py: Python<'_>) -> String {
        self.read_global(
            |g| {
                format!(
                    "GlobalData(name=\"{}\", statements={}, edge_statements={}, node_statements={})",
                    g.name,
                    g.statements.len(),
                    g.edge_statements.len(),
                    g.node_statements.len()
                )
            },
            py,
        )
    }
}

impl PyGlobalData {
    fn snapshot(&self, py: Python<'_>) -> GlobalData {
        self.read_global(|g| g.clone(), py)
    }

    fn read_global<R>(&self, f: impl FnOnce(&GlobalData) -> R, py: Python<'_>) -> R {
        if let Some(graph) = &self.graph {
            let graph_ref = graph
                .bind(py)
                .extract::<PyRef<PyDotGraph>>()
                .expect("GlobalData: invalid graph reference");
            f(&graph_ref.graph.global_data)
        } else {
            f(&self.data)
        }
    }

    fn write_global(&mut self, f: impl FnOnce(&mut GlobalData), py: Python<'_>) {
        if let Some(graph) = &self.graph {
            let mut graph_ref = graph
                .bind(py)
                .extract::<PyRefMut<PyDotGraph>>()
                .expect("GlobalData: invalid graph reference");
            f(&mut graph_ref.graph.global_data)
        } else {
            f(&mut self.data)
        }
    }

    fn get_statement(&self, kind: StatementKind, key: &str, py: Python<'_>) -> Option<String> {
        self.read_global(
            |g| match kind {
                StatementKind::Graph => g.statements.get(key).cloned(),
                StatementKind::Edge => g.edge_statements.get(key).cloned(),
                StatementKind::Node => g.node_statements.get(key).cloned(),
            },
            py,
        )
    }

    fn snapshot_statements(&self, kind: StatementKind, py: Python<'_>) -> Vec<(String, String)> {
        self.read_global(
            |g| {
                let map = match kind {
                    StatementKind::Graph => &g.statements,
                    StatementKind::Edge => &g.edge_statements,
                    StatementKind::Node => &g.node_statements,
                };
                map.iter().map(|(k, v)| (k.clone(), v.clone())).collect()
            },
            py,
        )
    }

    fn insert_statement(
        &mut self,
        kind: StatementKind,
        key: String,
        value: String,
        py: Python<'_>,
    ) {
        self.write_global(
            |g| match kind {
                StatementKind::Graph => {
                    g.statements.insert(key, value);
                }
                StatementKind::Edge => {
                    g.edge_statements.insert(key, value);
                }
                StatementKind::Node => {
                    g.node_statements.insert(key, value);
                }
            },
            py,
        );
    }

    fn remove_statement(&mut self, kind: StatementKind, key: &str, py: Python<'_>) -> bool {
        let mut removed = false;
        self.write_global(
            |g| {
                removed = match kind {
                    StatementKind::Graph => g.statements.remove(key).is_some(),
                    StatementKind::Edge => g.edge_statements.remove(key).is_some(),
                    StatementKind::Node => g.node_statements.remove(key).is_some(),
                };
            },
            py,
        );
        removed
    }

    fn clear_statements(&mut self, kind: StatementKind, py: Python<'_>) {
        self.write_global(
            |g| match kind {
                StatementKind::Graph => g.statements.clear(),
                StatementKind::Edge => g.edge_statements.clear(),
                StatementKind::Node => g.node_statements.clear(),
            },
            py,
        );
    }
}

impl PyGlobalStatements {
    fn get_value(&self, py: Python<'_>, key: &str) -> Option<String> {
        let data = self
            .data
            .bind(py)
            .extract::<PyRef<PyGlobalData>>()
            .expect("GlobalStatements: invalid data reference");
        data.get_statement(self.kind, key, py)
    }

    fn snapshot(&self, py: Python<'_>) -> Vec<(String, String)> {
        let data = self
            .data
            .bind(py)
            .extract::<PyRef<PyGlobalData>>()
            .expect("GlobalStatements: invalid data reference");
        data.snapshot_statements(self.kind, py)
    }
}

/// Edge data with orientation.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "EdgeData")]
#[derive(Clone, Debug, PartialEq, Eq)]
struct PyEdgeData {
    data: EdgeData<DotEdgeData>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyEdgeData {
    /// Create edge data from orientation and DotEdgeData.
    #[new]
    fn new(orientation: &Bound<'_, PyAny>, data: &Bound<'_, PyAny>) -> PyResult<Self> {
        let orientation = extract_orientation(orientation)?;
        let data = extract_dot_edge_data(data)?;
        Ok(Self {
            data: EdgeData { orientation, data },
        })
    }

    /// Orientation of this edge.
    #[getter]
    fn orientation(&self) -> PyOrientation {
        PyOrientation {
            orientation: self.data.orientation,
        }
    }

    /// DOT edge data.
    #[getter]
    fn data(&self) -> PyDotEdgeData {
        PyDotEdgeData {
            data: DotEdgeDataSource::Owned(self.data.data.clone()),
        }
    }

    /// Debug-style representation.
    fn __repr__(&self) -> String {
        format!("EdgeData(orientation={:?})", self.data.orientation)
    }
}

/// Hedge pairing (paired, unpaired, split).
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "HedgePair")]
#[derive(Clone, Debug, PartialEq, Eq)]
struct PyHedgePair {
    pair: HedgePair,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyHedgePair {
    /// Kind of hedge pair: \"paired\", \"unpaired\", or \"split\".
    #[getter]
    fn kind(&self) -> String {
        match self.pair {
            HedgePair::Paired { .. } => "paired",
            HedgePair::Unpaired { .. } => "unpaired",
            HedgePair::Split { .. } => "split",
        }
        .to_string()
    }

    /// Source hedge (paired/split).
    #[getter]
    fn source(&self) -> Option<PyHedge> {
        match self.pair {
            HedgePair::Paired { source, .. } => Some(PyHedge { hedge: source }),
            HedgePair::Split { source, .. } => Some(PyHedge { hedge: source }),
            HedgePair::Unpaired { .. } => None,
        }
    }

    /// Sink hedge (paired/split).
    #[getter]
    fn sink(&self) -> Option<PyHedge> {
        match self.pair {
            HedgePair::Paired { sink, .. } => Some(PyHedge { hedge: sink }),
            HedgePair::Split { sink, .. } => Some(PyHedge { hedge: sink }),
            HedgePair::Unpaired { .. } => None,
        }
    }

    /// Hedge (unpaired only).
    #[getter]
    fn hedge(&self) -> Option<PyHedge> {
        match self.pair {
            HedgePair::Unpaired { hedge, .. } => Some(PyHedge { hedge }),
            _ => None,
        }
    }

    /// Flow for unpaired hedge.
    #[getter]
    fn flow(&self) -> Option<PyFlow> {
        match self.pair {
            HedgePair::Unpaired { flow, .. } => Some(PyFlow { flow }),
            _ => None,
        }
    }

    /// Which side is split for split edges.
    #[getter]
    fn split(&self) -> Option<PyFlow> {
        match self.pair {
            HedgePair::Split { split, .. } => Some(PyFlow { flow: split }),
            _ => None,
        }
    }

    /// Debug-style representation.
    fn __repr__(&self) -> String {
        format!("HedgePair({})", self.kind())
    }
}

/// Subgraph represented as a bitset of hedges.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "Subgraph")]
#[derive(Clone, Debug, PartialEq, Eq)]
struct PySubgraph {
    subgraph: SuBitGraph,
}

#[gen_stub_pymethods]
#[pymethods]
impl PySubgraph {
    /// Create an empty subgraph of the given size (number of hedges).
    #[classmethod]
    fn empty(_cls: &Bound<'_, PyType>, size: usize) -> Self {
        Self {
            subgraph: SuBitGraph::empty(size),
        }
    }

    /// Create a full subgraph including all hedges.
    #[classmethod]
    fn full(_cls: &Bound<'_, PyType>, size: usize) -> Self {
        Self {
            subgraph: SuBitGraph::full(size),
        }
    }

    /// Create a subgraph from a list of hedges.
    #[classmethod]
    fn from_hedges(_cls: &Bound<'_, PyType>, size: usize, hedges: Vec<PyHedge>) -> PyResult<Self> {
        let mut subgraph = SuBitGraph::empty(size);
        for h in hedges {
            subgraph.add(h.hedge);
        }
        Ok(Self { subgraph })
    }

    /// List included hedges.
    fn to_hedges(&self) -> Vec<PyHedge> {
        self.subgraph
            .included_iter()
            .map(|h| PyHedge { hedge: h })
            .collect()
    }

    /// Total number of hedges in the parent graph.
    fn size(&self) -> usize {
        self.subgraph.size()
    }

    /// Number of included hedges.
    fn n_included(&self) -> usize {
        self.subgraph.n_included()
    }

    /// Whether a hedge is included.
    fn includes(&self, hedge: &Bound<'_, PyAny>) -> PyResult<bool> {
        let hedge = extract_hedge(hedge)?;
        Ok(self.subgraph.includes(&hedge))
    }

    /// Union with another subgraph.
    fn union(&self, other: &PySubgraph) -> Self {
        Self {
            subgraph: self.subgraph.union(&other.subgraph),
        }
    }

    /// Intersection with another subgraph.
    fn intersection(&self, other: &PySubgraph) -> Self {
        Self {
            subgraph: self.subgraph.intersection(&other.subgraph),
        }
    }

    /// Symmetric difference with another subgraph.
    fn sym_diff(&self, other: &PySubgraph) -> Self {
        Self {
            subgraph: self.subgraph.sym_diff(&other.subgraph),
        }
    }

    /// Subtract another subgraph.
    fn subtract(&self, other: &PySubgraph) -> Self {
        Self {
            subgraph: self.subgraph.subtract(&other.subgraph),
        }
    }

    /// Debug-style representation.
    fn __repr__(&self) -> String {
        format!(
            "Subgraph(size={}, included={})",
            self.size(),
            self.n_included()
        )
    }
}

/// Cycle represented as a subgraph and optional loop count.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "Cycle")]
#[derive(Clone, Debug, PartialEq, Eq)]
struct PyCycle {
    cycle: Cycle,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyCycle {
    /// Subgraph filter for this cycle.
    #[getter]
    fn filter(&self) -> PySubgraph {
        PySubgraph {
            subgraph: self.cycle.filter.clone(),
        }
    }

    /// Optional loop count for this cycle.
    #[getter]
    fn loop_count(&self) -> Option<usize> {
        self.cycle.loop_count
    }

    /// Debug-style representation.
    fn __repr__(&self) -> String {
        let count = match self.cycle.loop_count {
            Some(v) => v.to_string(),
            None => "None".to_string(),
        };
        format!(
            "Cycle(included={}, loop_count={})",
            self.cycle.filter.n_included(),
            count
        )
    }
}

/// Oriented cut represented by left/right subgraphs.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "OrientedCut")]
#[derive(Clone, Debug, PartialEq, Eq)]
struct PyOrientedCut {
    cut: OrientedCut,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyOrientedCut {
    /// Left side of the cut.
    #[getter]
    fn left(&self) -> PySubgraph {
        PySubgraph {
            subgraph: self.cut.left.clone(),
        }
    }

    /// Right side of the cut.
    #[getter]
    fn right(&self) -> PySubgraph {
        PySubgraph {
            subgraph: self.cut.right.clone(),
        }
    }

    /// Debug-style representation.
    fn __repr__(&self) -> String {
        format!("OrientedCut({})", self.cut)
    }
}

/// Hedge node with internal graph and hairs.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "HedgeNode")]
#[derive(Clone, Debug, PartialEq, Eq)]
struct PyHedgeNode {
    node: HedgeNode,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyHedgeNode {
    /// Create a hedge node from internal graph and hairs subgraphs.
    #[new]
    fn new(internal_graph: &PySubgraph, hairs: &PySubgraph) -> Self {
        let internal = InternalSubGraph {
            filter: internal_graph.subgraph.clone(),
            loopcount: None,
        };
        Self {
            node: HedgeNode {
                internal_graph: internal,
                hairs: hairs.subgraph.clone(),
            },
        }
    }

    /// Internal subgraph.
    #[getter]
    fn internal_graph(&self) -> PySubgraph {
        PySubgraph {
            subgraph: self.node.internal_graph.filter.clone(),
        }
    }

    /// Hair subgraph.
    #[getter]
    fn hairs(&self) -> PySubgraph {
        PySubgraph {
            subgraph: self.node.hairs.clone(),
        }
    }

    /// Debug-style representation.
    fn __repr__(&self) -> String {
        format!(
            "HedgeNode(internal={}, hairs={})",
            self.node.internal_graph.filter.n_included(),
            self.node.hairs.n_included()
        )
    }
}

/// Traversal tree produced by DFS/BFS.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "TraversalTree")]
#[derive(Clone, Debug)]
struct PyTraversalTree {
    tree: SimpleTraversalTree,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyTraversalTree {
    /// Subgraph corresponding to the tree edges.
    fn tree_subgraph(&self) -> PySubgraph {
        PySubgraph {
            subgraph: self.tree.tree_subgraph.clone(),
        }
    }

    /// Node order from traversal.
    fn node_order(&self) -> Vec<PyNodeIndex> {
        self.tree
            .node_order()
            .into_iter()
            .map(|n| PyNodeIndex { node: n })
            .collect()
    }

    /// Covers a subgraph with the traversal tree.
    fn covers(&self, subgraph: &PySubgraph) -> PySubgraph {
        PySubgraph {
            subgraph: self.tree.covers(&subgraph.subgraph),
        }
    }

    /// Iterate hedges as (hedge, kind, root_hedge).
    fn iter_hedges(&self) -> Vec<(PyHedge, String, Option<PyHedge>)> {
        self.tree
            .iter_hedges()
            .map(|(h, root, root_hedge)| {
                let kind = match root {
                    linnet::half_edge::tree::TTRoot::Root => "root",
                    linnet::half_edge::tree::TTRoot::Child(_) => "child",
                    linnet::half_edge::tree::TTRoot::None => "none",
                }
                .to_string();
                (
                    PyHedge { hedge: h },
                    kind,
                    root_hedge.map(|rh| PyHedge { hedge: rh }),
                )
            })
            .collect()
    }
}

/// DOT-backed hedge graph.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "DotGraph")]
#[derive(Clone, Debug)]
struct PyDotGraph {
    graph: DotGraph,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyDotGraph {
    /// Parse a DOT string into a graph.
    #[classmethod]
    fn from_string(_cls: &Bound<'_, PyType>, s: &str) -> PyResult<Self> {
        let graph = DotGraph::from_string(s).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { graph })
    }

    /// Parse a DOT string into multiple graphs.
    #[classmethod]
    fn from_string_set(_cls: &Bound<'_, PyType>, s: &str) -> PyResult<Vec<Self>> {
        let set = DotGraphSet::from_string(s).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(set.into_iter().map(|graph| Self { graph }).collect())
    }

    /// Parse a DOT file into a graph.
    #[classmethod]
    fn from_file(_cls: &Bound<'_, PyType>, path: &str) -> PyResult<Self> {
        let graph = DotGraph::from_file(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { graph })
    }

    /// Serialize graph to DOT for debugging.
    fn debug_dot(&self) -> String {
        self.graph.debug_dot()
    }

    /// Global graph attributes.
    #[getter]
    fn global_data(slf: PyRef<'_, Self>) -> PyResult<PyGlobalData> {
        let py = slf.py();
        let data = slf.graph.global_data.clone();
        let graph_obj = slf.into_pyobject(py)?.unbind().into();
        Ok(PyGlobalData {
            data,
            graph: Some(graph_obj),
        })
    }

    #[setter]
    fn set_global_data(&mut self, data: &PyGlobalData) {
        let _ = Python::try_attach(|py| {
            self.graph.global_data = data.snapshot(py);
        });
    }

    /// Serialize the full graph to DOT.
    fn dot(&self) -> String {
        self.graph.dot_of(&self.graph.full_filter())
    }

    /// Serialize a subgraph to DOT.
    fn dot_of(&self, subgraph: &PySubgraph) -> String {
        self.graph.dot_of(&subgraph.subgraph)
    }

    /// Access hedge/node/edge data via indexing.
    #[gen_stub(skip)]
    fn __getitem__(
        slf: Py<PyDotGraph>,
        key: &Bound<'_, PyAny>,
        py: Python<'_>,
    ) -> PyResult<Py<PyAny>> {
        if let Ok(h) = key.extract::<PyRef<PyHedge>>() {
            let graph_ref = slf.borrow(py);
            let data = graph_ref.graph.graph[h.hedge].clone();
            let obj = Py::new(py, PyDotHedgeData { data })?;
            return Ok(obj.into_any());
        }
        if let Ok(n) = key.extract::<PyRef<PyNodeIndex>>() {
            let graph_obj: Py<PyAny> = slf.into();
            let obj = Py::new(
                py,
                PyDotVertexData {
                    data: DotVertexDataSource::Graph {
                        graph: graph_obj,
                        node: n.node,
                    },
                },
            )?;
            return Ok(obj.into_any());
        }
        if let Ok(e) = key.extract::<PyRef<PyEdgeIndex>>() {
            let graph_obj: Py<PyAny> = slf.into();
            let obj = Py::new(
                py,
                PyDotEdgeData {
                    data: DotEdgeDataSource::Graph {
                        graph: graph_obj,
                        edge: e.edge,
                    },
                },
            )?;
            return Ok(obj.into_any());
        }
        Err(PyValueError::new_err(
            "expected Hedge, NodeIndex, or EdgeIndex",
        ))
    }

    /// Number of nodes.
    fn n_nodes(&self) -> usize {
        self.graph.n_nodes()
    }

    /// Number of edges.
    fn n_edges(&self) -> usize {
        self.graph.n_edges()
    }

    /// Number of hedges.
    fn n_hedges(&self) -> usize {
        self.graph.n_hedges()
    }

    /// Number of external hedges.
    fn n_externals(&self) -> usize {
        self.graph.n_externals()
    }

    /// Number of internal hedges.
    fn n_internals(&self) -> usize {
        self.graph.n_internals()
    }

    /// Subgraph including all hedges.
    fn full_filter(&self) -> PySubgraph {
        PySubgraph {
            subgraph: self.graph.full_filter(),
        }
    }

    /// Empty subgraph of this graph.
    fn empty_subgraph(&self) -> PySubgraph {
        PySubgraph {
            subgraph: self.graph.empty_subgraph(),
        }
    }

    /// Iterate edges in the full graph.
    fn iter_edges(&self) -> Vec<(PyHedgePair, PyEdgeIndex, PyEdgeData)> {
        self.graph
            .iter_edges()
            .map(|(pair, eid, data)| {
                let owned = EdgeData {
                    orientation: data.orientation,
                    data: data.data.clone(),
                };
                (
                    PyHedgePair { pair },
                    PyEdgeIndex { edge: eid },
                    PyEdgeData { data: owned },
                )
            })
            .collect()
    }

    /// Iterate edges within a subgraph.
    fn iter_edges_of(&self, subgraph: &PySubgraph) -> Vec<(PyHedgePair, PyEdgeIndex, PyEdgeData)> {
        self.graph
            .iter_edges_of(&subgraph.subgraph)
            .map(|(pair, eid, data)| {
                let owned = EdgeData {
                    orientation: data.orientation,
                    data: data.data.clone(),
                };
                (
                    PyHedgePair { pair },
                    PyEdgeIndex { edge: eid },
                    PyEdgeData { data: owned },
                )
            })
            .collect()
    }

    /// Iterate nodes in the full graph.
    fn iter_nodes(&self) -> Vec<(PyNodeIndex, Vec<PyHedge>, PyDotVertexData)> {
        self.graph
            .iter_nodes()
            .map(|(node, neighbors, data)| {
                let hedges = neighbors.map(|h| PyHedge { hedge: h }).collect();
                (
                    PyNodeIndex { node },
                    hedges,
                    PyDotVertexData {
                        data: DotVertexDataSource::Owned(data.clone()),
                    },
                )
            })
            .collect()
    }

    /// Iterate nodes within a subgraph.
    fn iter_nodes_of(
        &self,
        subgraph: &PySubgraph,
    ) -> Vec<(PyNodeIndex, Vec<PyHedge>, PyDotVertexData)> {
        self.graph
            .iter_nodes_of(&subgraph.subgraph)
            .map(|(node, neighbors, data)| {
                let hedges = neighbors.map(|h| PyHedge { hedge: h }).collect();
                (
                    PyNodeIndex { node },
                    hedges,
                    PyDotVertexData {
                        data: DotVertexDataSource::Owned(data.clone()),
                    },
                )
            })
            .collect()
    }

    /// Connected components of a subgraph.
    fn connected_components(&self, subgraph: &PySubgraph) -> Vec<PySubgraph> {
        self.graph
            .connected_components(&subgraph.subgraph)
            .into_iter()
            .map(|sg| PySubgraph { subgraph: sg })
            .collect()
    }

    /// Count connected components.
    fn count_connected_components(&self, subgraph: &PySubgraph) -> usize {
        self.graph.count_connected_components(&subgraph.subgraph)
    }

    /// Whether a subgraph is connected.
    fn is_connected(&self, subgraph: &PySubgraph) -> bool {
        self.graph.is_connected(&subgraph.subgraph)
    }

    /// Depth-first traversal from a root node.
    fn depth_first_traverse(
        &self,
        subgraph: &PySubgraph,
        root_node: &PyNodeIndex,
        include_hedge: Option<&PyHedge>,
    ) -> PyResult<PyTraversalTree> {
        let root = root_node.node;
        let include = include_hedge.map(|h| h.hedge);
        let tree = SimpleTraversalTree::depth_first_traverse(
            &self.graph.graph,
            &subgraph.subgraph,
            &root,
            include,
        )
        .map_err(to_py_err)?;
        Ok(PyTraversalTree { tree })
    }

    /// Breadth-first traversal from a root node.
    fn breadth_first_traverse(
        &self,
        subgraph: &PySubgraph,
        root_node: &PyNodeIndex,
        include_hedge: Option<&PyHedge>,
    ) -> PyResult<PyTraversalTree> {
        let root = root_node.node;
        let include = include_hedge.map(|h| h.hedge);
        let tree = SimpleTraversalTree::breadth_first_traverse(
            &self.graph.graph,
            &subgraph.subgraph,
            &root,
            include,
        )
        .map_err(to_py_err)?;
        Ok(PyTraversalTree { tree })
    }

    /// Bridges in the full graph.
    fn bridges(&self) -> PySubgraph {
        PySubgraph {
            subgraph: self.graph.bridges(),
        }
    }

    /// Bridges within a subgraph.
    fn bridges_of(&self, subgraph: &PySubgraph) -> PySubgraph {
        PySubgraph {
            subgraph: self.graph.bridges_of(&subgraph.subgraph),
        }
    }

    /// Cycle basis of the full graph.
    fn cycle_basis(&self) -> (Vec<PyCycle>, PyTraversalTree) {
        let (cycles, tree) = self.graph.cycle_basis();
        let cycles = cycles.into_iter().map(|c| PyCycle { cycle: c }).collect();
        (cycles, PyTraversalTree { tree })
    }

    /// Cycle basis of a subgraph.
    fn cycle_basis_of(&self, subgraph: &PySubgraph) -> (Vec<PyCycle>, PyTraversalTree) {
        let (cycles, tree) = self.graph.cycle_basis_of(&subgraph.subgraph);
        let cycles = cycles.into_iter().map(|c| PyCycle { cycle: c }).collect();
        (cycles, PyTraversalTree { tree })
    }

    /// All spanning forests of the full graph.
    fn all_spanning_forests(&self) -> Vec<PySubgraph> {
        let full = self.graph.full_filter();
        self.graph
            .all_spanning_forests_of(&full)
            .into_iter()
            .map(|sg| PySubgraph { subgraph: sg })
            .collect()
    }

    /// All spanning forests of a subgraph.
    fn all_spanning_forests_of(&self, subgraph: &PySubgraph) -> Vec<PySubgraph> {
        self.graph
            .all_spanning_forests_of(&subgraph.subgraph)
            .into_iter()
            .map(|sg| PySubgraph { subgraph: sg })
            .collect()
    }

    /// Combine nodes into a single hedge node.
    fn combine_to_single_hedgenode(&self, nodes: Vec<PyNodeIndex>) -> PyResult<PyHedgeNode> {
        let ids = nodes.into_iter().map(|n| n.node).collect::<Vec<_>>();
        let node = self.graph.combine_to_single_hedgenode(&ids);
        Ok(PyHedgeNode { node })
    }

    /// All cuts between two hedge nodes.
    fn all_cuts(
        &self,
        source: &PyHedgeNode,
        target: &PyHedgeNode,
    ) -> Vec<(PySubgraph, PyOrientedCut, PySubgraph)> {
        self.graph
            .all_cuts(source.node.clone(), target.node.clone())
            .into_iter()
            .map(|(s_side, cut, t_side)| {
                (
                    PySubgraph { subgraph: s_side },
                    PyOrientedCut { cut },
                    PySubgraph { subgraph: t_side },
                )
            })
            .collect()
    }

    /// All cuts between two sets of node indices.
    fn all_cuts_from_ids(
        &self,
        source: Vec<PyNodeIndex>,
        target: Vec<PyNodeIndex>,
    ) -> PyResult<Vec<(PySubgraph, PyOrientedCut, PySubgraph)>> {
        let source_ids = source.into_iter().map(|n| n.node).collect::<Vec<_>>();
        let target_ids = target.into_iter().map(|n| n.node).collect::<Vec<_>>();
        let cuts = self.graph.all_cuts_from_ids(&source_ids, &target_ids);
        Ok(cuts
            .into_iter()
            .map(|(s_side, cut, t_side)| {
                (
                    PySubgraph { subgraph: s_side },
                    PyOrientedCut { cut },
                    PySubgraph { subgraph: t_side },
                )
            })
            .collect())
    }

    /// Contract a subgraph into a single node, deleting its edges.
    #[pyo3(signature = (subgraph, node_data_merge=None))]
    fn contract_subgraph(
        &mut self,
        subgraph: &PySubgraph,
        node_data_merge: Option<&PyDotVertexData>,
    ) {
        let node_data = match node_data_merge {
            Some(obj) => {
                Python::try_attach(|py| obj.snapshot(py)).unwrap_or_else(DotVertexData::empty)
            }
            None => DotVertexData::empty(),
        };
        self.graph
            .graph
            .contract_subgraph(&subgraph.subgraph, node_data);
    }

    /// Join two graphs, matching dangling edges via a Python callback.
    fn join(
        &self,
        other: &PyDotGraph,
        matching_fn: Py<PyAny>,
        merge_fn: Py<PyAny>,
    ) -> PyResult<Self> {
        let error: RefCell<Option<PyErr>> = RefCell::new(None);
        let left = self.graph.graph.clone();
        let right = other.graph.graph.clone();
        let matching = matching_fn;
        let merging = merge_fn;

        let result = left.join(
            right,
            |f1, d1, f2, d2| {
                if error.borrow().is_some() {
                    return false;
                }
                Python::attach(|py| {
                    let py_f1 = PyFlow { flow: f1 };
                    let py_f2 = PyFlow { flow: f2 };
                    let py_d1 = PyEdgeData {
                        data: EdgeData {
                            orientation: d1.orientation,
                            data: d1.data.clone(),
                        },
                    };
                    let py_d2 = PyEdgeData {
                        data: EdgeData {
                            orientation: d2.orientation,
                            data: d2.data.clone(),
                        },
                    };

                    match matching.bind(py).call1((py_f1, py_d1, py_f2, py_d2)) {
                        Ok(val) => match val.extract::<bool>() {
                            Ok(b) => b,
                            Err(e) => {
                                *error.borrow_mut() = Some(e);
                                false
                            }
                        },
                        Err(e) => {
                            *error.borrow_mut() = Some(e);
                            false
                        }
                    }
                })
            },
            |f1, d1, f2, d2| {
                if error.borrow().is_some() {
                    return (
                        Flow::Source,
                        EdgeData {
                            orientation: Orientation::Undirected,
                            data: DotEdgeData::empty(),
                        },
                    );
                }
                Python::attach(|py| {
                    let py_f1 = PyFlow { flow: f1 };
                    let py_f2 = PyFlow { flow: f2 };
                    let py_d1 = PyEdgeData { data: d1 };
                    let py_d2 = PyEdgeData { data: d2 };

                    match merging.bind(py).call1((py_f1, py_d1, py_f2, py_d2)) {
                        Ok(val) => match extract_merge_result(&val) {
                            Ok(out) => out,
                            Err(e) => {
                                *error.borrow_mut() = Some(e);
                                (
                                    Flow::Source,
                                    EdgeData {
                                        orientation: Orientation::Undirected,
                                        data: DotEdgeData::empty(),
                                    },
                                )
                            }
                        },
                        Err(e) => {
                            *error.borrow_mut() = Some(e);
                            (
                                Flow::Source,
                                EdgeData {
                                    orientation: Orientation::Undirected,
                                    data: DotEdgeData::empty(),
                                },
                            )
                        }
                    }
                })
            },
        );

        if let Some(err) = error.into_inner() {
            return Err(err);
        }

        let graph = result.map_err(to_py_err)?;
        Ok(Self {
            graph: DotGraph {
                graph,
                global_data: self.graph.global_data.clone(),
            },
        })
    }

    /// In-place join, matching dangling edges via a Python callback.
    fn join_mut(
        &mut self,
        other: &PyDotGraph,
        matching_fn: Py<PyAny>,
        merge_fn: Py<PyAny>,
    ) -> PyResult<()> {
        let error: RefCell<Option<PyErr>> = RefCell::new(None);
        let other_graph = other.graph.graph.clone();
        let matching = matching_fn;
        let merging = merge_fn;

        let result = self.graph.graph.join_mut(
            other_graph,
            |f1, d1, f2, d2| {
                if error.borrow().is_some() {
                    return false;
                }
                Python::attach(|py| {
                    let py_f1 = PyFlow { flow: f1 };
                    let py_f2 = PyFlow { flow: f2 };
                    let py_d1 = PyEdgeData {
                        data: EdgeData {
                            orientation: d1.orientation,
                            data: d1.data.clone(),
                        },
                    };
                    let py_d2 = PyEdgeData {
                        data: EdgeData {
                            orientation: d2.orientation,
                            data: d2.data.clone(),
                        },
                    };

                    match matching.bind(py).call1((py_f1, py_d1, py_f2, py_d2)) {
                        Ok(val) => match val.extract::<bool>() {
                            Ok(b) => b,
                            Err(e) => {
                                *error.borrow_mut() = Some(e);
                                false
                            }
                        },
                        Err(e) => {
                            *error.borrow_mut() = Some(e);
                            false
                        }
                    }
                })
            },
            |f1, d1, f2, d2| {
                if error.borrow().is_some() {
                    return (
                        Flow::Source,
                        EdgeData {
                            orientation: Orientation::Undirected,
                            data: DotEdgeData::empty(),
                        },
                    );
                }
                Python::attach(|py| {
                    let py_f1 = PyFlow { flow: f1 };
                    let py_f2 = PyFlow { flow: f2 };
                    let py_d1 = PyEdgeData { data: d1 };
                    let py_d2 = PyEdgeData { data: d2 };

                    match merging.bind(py).call1((py_f1, py_d1, py_f2, py_d2)) {
                        Ok(val) => match extract_merge_result(&val) {
                            Ok(out) => out,
                            Err(e) => {
                                *error.borrow_mut() = Some(e);
                                (
                                    Flow::Source,
                                    EdgeData {
                                        orientation: Orientation::Undirected,
                                        data: DotEdgeData::empty(),
                                    },
                                )
                            }
                        },
                        Err(e) => {
                            *error.borrow_mut() = Some(e);
                            (
                                Flow::Source,
                                EdgeData {
                                    orientation: Orientation::Undirected,
                                    data: DotEdgeData::empty(),
                                },
                            )
                        }
                    }
                })
            },
        );

        if let Some(err) = error.into_inner() {
            return Err(err);
        }

        result.map_err(to_py_err)
    }

    /// Extract a subgraph with Python callbacks to transform edge/node data.
    fn extract(
        &mut self,
        subgraph: &PySubgraph,
        split_edge_fn: Py<PyAny>,
        internal_data: Py<PyAny>,
        split_node: Py<PyAny>,
        owned_node: Py<PyAny>,
    ) -> PyResult<Self> {
        let error: RefCell<Option<PyErr>> = RefCell::new(None);
        let split_edge = split_edge_fn;
        let internal = internal_data;
        let split_node_fn = split_node;
        let owned_node_fn = owned_node;

        let extracted = self.graph.graph.extract(
            &subgraph.subgraph,
            |edge_ref| {
                if error.borrow().is_some() {
                    return EdgeData {
                        orientation: Orientation::Undirected,
                        data: DotEdgeData::empty(),
                    };
                }
                Python::attach(|py| {
                    let py_edge = PyEdgeData {
                        data: EdgeData {
                            orientation: edge_ref.orientation,
                            data: edge_ref.data.clone(),
                        },
                    };
                    match split_edge.bind(py).call1((py_edge,)) {
                        Ok(val) => match extract_edge_data(&val) {
                            Ok(out) => out,
                            Err(e) => {
                                *error.borrow_mut() = Some(e);
                                EdgeData {
                                    orientation: Orientation::Undirected,
                                    data: DotEdgeData::empty(),
                                }
                            }
                        },
                        Err(e) => {
                            *error.borrow_mut() = Some(e);
                            EdgeData {
                                orientation: Orientation::Undirected,
                                data: DotEdgeData::empty(),
                            }
                        }
                    }
                })
            },
            |edge_owned| {
                if error.borrow().is_some() {
                    return EdgeData {
                        orientation: Orientation::Undirected,
                        data: DotEdgeData::empty(),
                    };
                }
                Python::attach(|py| {
                    let py_edge = PyEdgeData { data: edge_owned };
                    match internal.bind(py).call1((py_edge,)) {
                        Ok(val) => match extract_edge_data(&val) {
                            Ok(out) => out,
                            Err(e) => {
                                *error.borrow_mut() = Some(e);
                                EdgeData {
                                    orientation: Orientation::Undirected,
                                    data: DotEdgeData::empty(),
                                }
                            }
                        },
                        Err(e) => {
                            *error.borrow_mut() = Some(e);
                            EdgeData {
                                orientation: Orientation::Undirected,
                                data: DotEdgeData::empty(),
                            }
                        }
                    }
                })
            },
            |node_ref| {
                if error.borrow().is_some() {
                    return DotVertexData::empty();
                }
                Python::attach(|py| {
                    let py_node = PyDotVertexData {
                        data: DotVertexDataSource::Owned(node_ref.clone()),
                    };
                    match split_node_fn.bind(py).call1((py_node,)) {
                        Ok(val) => match extract_dot_vertex_data(&val) {
                            Ok(out) => out,
                            Err(e) => {
                                *error.borrow_mut() = Some(e);
                                DotVertexData::empty()
                            }
                        },
                        Err(e) => {
                            *error.borrow_mut() = Some(e);
                            DotVertexData::empty()
                        }
                    }
                })
            },
            |node_owned| {
                if error.borrow().is_some() {
                    return DotVertexData::empty();
                }
                Python::attach(|py| {
                    let py_node = PyDotVertexData {
                        data: DotVertexDataSource::Owned(node_owned),
                    };
                    match owned_node_fn.bind(py).call1((py_node,)) {
                        Ok(val) => match extract_dot_vertex_data(&val) {
                            Ok(out) => out,
                            Err(e) => {
                                *error.borrow_mut() = Some(e);
                                DotVertexData::empty()
                            }
                        },
                        Err(e) => {
                            *error.borrow_mut() = Some(e);
                            DotVertexData::empty()
                        }
                    }
                })
            },
        );

        if let Some(err) = error.into_inner() {
            return Err(err);
        }

        Ok(Self {
            graph: DotGraph {
                graph: extracted,
                global_data: self.graph.global_data.clone(),
            },
        })
    }
}

/// Builder for constructing DOT graphs programmatically.
#[gen_stub_pyclass]
#[pyclass(name = "DotGraphBuilder")]
struct PyDotGraphBuilder {
    builder: HedgeGraphBuilder<DotEdgeData, DotVertexData, DotHedgeData>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyDotGraphBuilder {
    /// Create a new, empty graph builder.
    #[new]
    fn new() -> Self {
        Self {
            builder: HedgeGraphBuilder::new(),
        }
    }

    /// Add a node and return its index.
    #[pyo3(signature = (data=None))]
    fn add_node(&mut self, data: Option<&PyDotVertexData>) -> PyResult<PyNodeIndex> {
        let data = match data {
            Some(obj) => {
                Python::try_attach(|py| obj.snapshot(py)).unwrap_or_else(DotVertexData::empty)
            }
            None => DotVertexData::empty(),
        };
        let node = self.builder.add_node(data);
        Ok(PyNodeIndex { node })
    }

    /// Add an edge between two nodes.
    #[pyo3(signature = (source, sink, data=None, orientation=None, source_hedge=None, sink_hedge=None))]
    fn add_edge(
        &mut self,
        source: &PyNodeIndex,
        sink: &PyNodeIndex,
        data: Option<&PyDotEdgeData>,
        orientation: Option<&PyOrientation>,
        source_hedge: Option<&PyDotHedgeData>,
        sink_hedge: Option<&PyDotHedgeData>,
    ) -> PyResult<()> {
        let data = match data {
            Some(obj) => {
                Python::try_attach(|py| obj.snapshot(py)).unwrap_or_else(DotEdgeData::empty)
            }
            None => DotEdgeData::empty(),
        };
        let orientation = match orientation {
            Some(obj) => obj.orientation,
            None => Orientation::Default,
        };
        let source_hedge = match source_hedge {
            Some(obj) => obj.data.clone(),
            None => DotHedgeData::default(),
        };
        let sink_hedge = match sink_hedge {
            Some(obj) => obj.data.clone(),
            None => DotHedgeData::default(),
        };

        self.builder.add_edge(
            HedgeData {
                data: source_hedge,
                is_in_subgraph: false,
                node: source.node,
            },
            HedgeData {
                data: sink_hedge,
                is_in_subgraph: false,
                node: sink.node,
            },
            data,
            orientation,
        );
        Ok(())
    }

    /// Add a dangling (external) edge incident to a node.
    #[pyo3(signature = (source, data=None, orientation=None, flow=None, hedge=None))]
    fn add_external_edge(
        &mut self,
        source: &PyNodeIndex,
        data: Option<&PyDotEdgeData>,
        orientation: Option<&PyOrientation>,
        flow: Option<&PyFlow>,
        hedge: Option<&PyDotHedgeData>,
    ) -> PyResult<()> {
        let data = match data {
            Some(obj) => {
                Python::try_attach(|py| obj.snapshot(py)).unwrap_or_else(DotEdgeData::empty)
            }
            None => DotEdgeData::empty(),
        };
        let orientation = match orientation {
            Some(obj) => obj.orientation,
            None => Orientation::Default,
        };
        let flow = match flow {
            Some(obj) => obj.flow,
            None => Flow::Source,
        };
        let hedge = match hedge {
            Some(obj) => obj.data.clone(),
            None => DotHedgeData::default(),
        };

        self.builder.add_external_edge(
            HedgeData {
                data: hedge,
                is_in_subgraph: false,
                node: source.node,
            },
            data,
            orientation,
            flow,
        );
        Ok(())
    }

    /// Build the graph and reset the builder.
    fn build(&mut self) -> PyDotGraph {
        let builder = std::mem::take(&mut self.builder);
        let graph = builder.build::<DefaultNodeStore<DotVertexData>>();
        PyDotGraph {
            graph: DotGraph {
                graph,
                global_data: GlobalData::from(()),
            },
        }
    }
}

#[pymodule]
fn linnet_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyHedge>()?;
    m.add_class::<PyNodeIndex>()?;
    m.add_class::<PyEdgeIndex>()?;
    m.add_class::<PyFlow>()?;
    m.add_class::<PyOrientation>()?;
    m.add_class::<PyDotVertexData>()?;
    m.add_class::<PyNodeStatements>()?;
    m.add_class::<PyDotHedgeData>()?;
    m.add_class::<PyDotEdgeData>()?;
    m.add_class::<PyEdgeStatements>()?;
    m.add_class::<PyGlobalData>()?;
    m.add_class::<PyGlobalStatements>()?;
    m.add_class::<PyEdgeData>()?;
    m.add_class::<PyHedgePair>()?;
    m.add_class::<PySubgraph>()?;
    m.add_class::<PyCycle>()?;
    m.add_class::<PyOrientedCut>()?;
    m.add_class::<PyHedgeNode>()?;
    m.add_class::<PyTraversalTree>()?;
    m.add_class::<PyDotGraph>()?;
    m.add_class::<PyDotGraphBuilder>()?;
    Ok(())
}

submit! {
    PyMethodsInfo {
        struct_id: std::any::TypeId::of::<PyDotGraph>,
        attrs: &[],
        getters: &[],
        setters: &[],
        methods: &[
            MethodInfo {
                name: "__getitem__",
                parameters: &[
                    ParameterInfo {
                        name: "key",
                        kind: ParameterKind::PositionalOrKeyword,
                        type_info: <PyHedge as PyStubType>::type_input,
                        default: ParameterDefault::None,
                    }
                ],
                r#return: <PyDotHedgeData as PyStubType>::type_output,
                doc: "",
                r#type: MethodType::Instance,
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
            MethodInfo {
                name: "__getitem__",
                parameters: &[
                    ParameterInfo {
                        name: "key",
                        kind: ParameterKind::PositionalOrKeyword,
                        type_info: <PyNodeIndex as PyStubType>::type_input,
                        default: ParameterDefault::None,
                    }
                ],
                r#return: <PyDotVertexData as PyStubType>::type_output,
                doc: "",
                r#type: MethodType::Instance,
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
            MethodInfo {
                name: "__getitem__",
                parameters: &[
                    ParameterInfo {
                        name: "key",
                        kind: ParameterKind::PositionalOrKeyword,
                        type_info: <PyEdgeIndex as PyStubType>::type_input,
                        default: ParameterDefault::None,
                    }
                ],
                r#return: <PyDotEdgeData as PyStubType>::type_output,
                doc: "",
                r#type: MethodType::Instance,
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
        ],
        file: file!(),
        line: line!(),
        column: column!(),
    }
}

define_stub_info_gatherer!(stub_info);

fn extract_hedge(obj: &Bound<'_, PyAny>) -> PyResult<Hedge> {
    if let Ok(h) = obj.extract::<PyRef<PyHedge>>() {
        Ok(h.hedge)
    } else {
        Ok(Hedge(obj.extract::<usize>()?))
    }
}

fn extract_flow(obj: &Bound<'_, PyAny>) -> PyResult<Flow> {
    if let Ok(f) = obj.extract::<PyRef<PyFlow>>() {
        Ok(f.flow)
    } else if let Ok(s) = obj.extract::<&str>() {
        match s.to_ascii_lowercase().as_str() {
            "source" => Ok(Flow::Source),
            "sink" => Ok(Flow::Sink),
            _ => Err(PyValueError::new_err("invalid flow")),
        }
    } else {
        Err(PyValueError::new_err("invalid flow"))
    }
}

fn extract_orientation(obj: &Bound<'_, PyAny>) -> PyResult<Orientation> {
    if let Ok(o) = obj.extract::<PyRef<PyOrientation>>() {
        Ok(o.orientation)
    } else if let Ok(s) = obj.extract::<&str>() {
        match s.to_ascii_lowercase().as_str() {
            "default" => Ok(Orientation::Default),
            "reversed" => Ok(Orientation::Reversed),
            "undirected" => Ok(Orientation::Undirected),
            _ => Err(PyValueError::new_err("invalid orientation")),
        }
    } else {
        Err(PyValueError::new_err("invalid orientation"))
    }
}

fn extract_dot_edge_data(obj: &Bound<'_, PyAny>) -> PyResult<DotEdgeData> {
    let data = obj.extract::<PyRef<PyDotEdgeData>>()?;
    Ok(data.snapshot(obj.py()))
}

fn extract_dot_vertex_data(obj: &Bound<'_, PyAny>) -> PyResult<DotVertexData> {
    let data = obj.extract::<PyRef<PyDotVertexData>>()?;
    Ok(data.snapshot(obj.py()))
}

fn extract_edge_data(obj: &Bound<'_, PyAny>) -> PyResult<EdgeData<DotEdgeData>> {
    if let Ok(ed) = obj.extract::<PyRef<PyEdgeData>>() {
        Ok(ed.data.clone())
    } else if let Ok(tuple) = obj.cast::<PyTuple>() {
        if tuple.len() != 2 {
            return Err(PyValueError::new_err("expected (orientation, data)"));
        }
        let item0 = tuple.get_item(0)?;
        let item1 = tuple.get_item(1)?;
        let orientation = extract_orientation(&item0)?;
        let data = extract_dot_edge_data(&item1)?;
        Ok(EdgeData { orientation, data })
    } else {
        Err(PyValueError::new_err("invalid edge data"))
    }
}

fn extract_merge_result(obj: &Bound<'_, PyAny>) -> PyResult<(Flow, EdgeData<DotEdgeData>)> {
    if let Ok(tuple) = obj.cast::<PyTuple>() {
        if tuple.len() != 2 {
            return Err(PyValueError::new_err("expected (flow, edge_data)"));
        }
        let item0 = tuple.get_item(0)?;
        let item1 = tuple.get_item(1)?;
        let flow = extract_flow(&item0)?;
        let data = extract_edge_data(&item1)?;
        Ok((flow, data))
    } else {
        Err(PyValueError::new_err("expected (flow, edge_data)"))
    }
}

fn to_py_err(err: HedgeGraphError) -> PyErr {
    PyValueError::new_err(err.to_string())
}

fn parse_compasspt(value: &str) -> Option<dot_parser::ast::CompassPt> {
    use dot_parser::ast::CompassPt;
    match value.to_ascii_lowercase().as_str() {
        "n" => Some(CompassPt::N),
        "ne" => Some(CompassPt::NE),
        "e" => Some(CompassPt::E),
        "se" => Some(CompassPt::SE),
        "s" => Some(CompassPt::S),
        "sw" => Some(CompassPt::SW),
        "w" => Some(CompassPt::W),
        "nw" => Some(CompassPt::NW),
        "c" => Some(CompassPt::C),
        "_" => Some(CompassPt::Underscore),
        _ => None,
    }
}
