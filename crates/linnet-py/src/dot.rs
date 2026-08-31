use std::collections::{BTreeMap, BTreeSet};

use dot_parser::ast::CompassPt;
use linnet::half_edge::involution::{EdgeIndex, Flow, Hedge, HedgePair};
use linnet::half_edge::nodestore::NodeStorageVec;
use linnet::half_edge::{HedgeGraph, NodeIndex};
use linnet::parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData, GlobalData};
use pyo3::class::gc::{PyTraverseError, PyVisit};
use pyo3::exceptions::{PyReferenceError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyDictMethods, PyType};
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use serde::{Deserialize, Serialize};

use crate::drawing::{
    copy_drawing, drawing_dict, PyEdgeDrawing, PyHalfEdgeDrawing, PyNodeDrawing, EDGE_FIELDS,
    HEDGE_FIELDS, NODE_FIELDS,
};
use crate::graph::{EdgeRecord, GraphState, HalfEdgeRecord, NodeRecord, PyGraph};
use crate::native_graph::{PyHedgeGraph, PyNodeStore};

/// DOT graph metadata kept separate from arbitrary Python element data.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "GlobalData")]
#[derive(Clone, Debug)]
pub struct PyGlobalData {
    pub(crate) inner: GlobalData,
}

impl Default for PyGlobalData {
    fn default() -> Self {
        Self {
            inner: GlobalData::from(()),
        }
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl PyGlobalData {
    #[new]
    #[pyo3(signature = (name="", *, payload=None, statements=None, edge_statements=None, node_statements=None))]
    fn new(
        name: &str,
        payload: Option<Vec<u8>>,
        statements: Option<BTreeMap<String, String>>,
        edge_statements: Option<BTreeMap<String, String>>,
        node_statements: Option<BTreeMap<String, String>>,
    ) -> Self {
        Self {
            inner: GlobalData {
                name: name.to_string(),
                payload,
                statements: statements.unwrap_or_default(),
                edge_statements: edge_statements.unwrap_or_default(),
                node_statements: node_statements.unwrap_or_default(),
            },
        }
    }

    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    #[setter]
    fn set_name(&mut self, name: String) {
        self.inner.name = name;
    }

    #[getter]
    fn payload(&self) -> Option<Vec<u8>> {
        self.inner.payload.clone()
    }

    #[setter]
    fn set_payload(&mut self, payload: Option<Vec<u8>>) {
        self.inner.payload = payload;
    }

    #[getter]
    fn statements(&self) -> BTreeMap<String, String> {
        self.inner.statements.clone()
    }

    #[setter]
    fn set_statements(&mut self, value: BTreeMap<String, String>) {
        self.inner.statements = value;
    }

    #[getter]
    fn edge_statements(&self) -> BTreeMap<String, String> {
        self.inner.edge_statements.clone()
    }

    #[setter]
    fn set_edge_statements(&mut self, value: BTreeMap<String, String>) {
        self.inner.edge_statements = value;
    }

    #[getter]
    fn node_statements(&self) -> BTreeMap<String, String> {
        self.inner.node_statements.clone()
    }

    #[setter]
    fn set_node_statements(&mut self, value: BTreeMap<String, String>) {
        self.inner.node_statements = value;
    }

    fn __repr__(&self) -> String {
        format!("GlobalData(name={:?})", self.inner.name)
    }
}

/// The DOT-representable record exchanged for one node by a `DotCodec`.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "DotVertexData")]
#[derive(Clone, Debug)]
pub struct PyDotVertexData {
    pub(crate) inner: DotVertexData,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyDotVertexData {
    #[new]
    #[pyo3(signature = (*, name=None, index=None, payload=None, statements=None))]
    fn new(
        name: Option<String>,
        index: Option<usize>,
        payload: Option<Vec<u8>>,
        statements: Option<BTreeMap<String, String>>,
    ) -> Self {
        Self {
            inner: DotVertexData {
                name,
                index: index.map(NodeIndex),
                payload,
                statements: statements.unwrap_or_default(),
            },
        }
    }

    #[getter]
    fn name(&self) -> Option<&str> {
        self.inner.name.as_deref()
    }

    #[getter]
    fn index(&self) -> Option<usize> {
        self.inner.index.map(|index| index.0)
    }

    #[getter]
    fn payload(&self) -> Option<Vec<u8>> {
        self.inner.payload.clone()
    }

    #[getter]
    fn statements(&self) -> BTreeMap<String, String> {
        self.inner.statements.clone()
    }
}

/// The DOT-representable record exchanged for one edge by a `DotCodec`.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "DotEdgeData")]
#[derive(Clone, Debug)]
pub struct PyDotEdgeData {
    pub(crate) inner: DotEdgeData,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyDotEdgeData {
    #[new]
    #[pyo3(signature = (*, edge_id=None, payload=None, statements=None, local_statements=None))]
    fn new(
        edge_id: Option<usize>,
        payload: Option<Vec<u8>>,
        statements: Option<BTreeMap<String, String>>,
        local_statements: Option<BTreeMap<String, String>>,
    ) -> Self {
        Self {
            inner: DotEdgeData {
                payload,
                statements: statements.unwrap_or_default(),
                local_statements: local_statements.unwrap_or_default(),
                edge_id: edge_id.map(EdgeIndex::from),
            },
        }
    }

    #[getter]
    fn edge_id(&self) -> Option<usize> {
        self.inner.edge_id.map(|index| index.0)
    }

    #[getter]
    fn payload(&self) -> Option<Vec<u8>> {
        self.inner.payload.clone()
    }

    #[getter]
    fn statements(&self) -> BTreeMap<String, String> {
        self.inner.statements.clone()
    }

    #[getter]
    fn local_statements(&self) -> BTreeMap<String, String> {
        self.inner.local_statements.clone()
    }
}

fn parse_compass(value: Option<&str>) -> PyResult<Option<CompassPt>> {
    value
        .map(|value| {
            Ok(match value.to_ascii_lowercase().as_str() {
                "n" => CompassPt::N,
                "ne" => CompassPt::NE,
                "e" => CompassPt::E,
                "se" => CompassPt::SE,
                "s" => CompassPt::S,
                "sw" => CompassPt::SW,
                "w" => CompassPt::W,
                "nw" => CompassPt::NW,
                "c" => CompassPt::C,
                "_" => CompassPt::Underscore,
                _ => return Err(PyValueError::new_err("invalid DOT compass point")),
            })
        })
        .transpose()
}

fn compass_name(value: CompassPt) -> &'static str {
    match value {
        CompassPt::N => "n",
        CompassPt::NE => "ne",
        CompassPt::E => "e",
        CompassPt::SE => "se",
        CompassPt::S => "s",
        CompassPt::SW => "sw",
        CompassPt::W => "w",
        CompassPt::NW => "nw",
        CompassPt::C => "c",
        CompassPt::Underscore => "_",
    }
}

/// The DOT-representable record exchanged for one half-edge by a `DotCodec`.
#[gen_stub_pyclass]
#[pyclass(from_py_object, name = "DotHalfEdgeData")]
#[derive(Clone, Debug)]
pub struct PyDotHalfEdgeData {
    pub(crate) inner: DotHedgeData,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyDotHalfEdgeData {
    #[new]
    #[pyo3(signature = (*, statement=None, index=None, payload=None, port_label=None, compass=None))]
    fn new(
        statement: Option<String>,
        index: Option<usize>,
        payload: Option<Vec<u8>>,
        port_label: Option<String>,
        compass: Option<&str>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: DotHedgeData {
                statement,
                id: index.map(Hedge),
                payload,
                port_label,
                compasspt: parse_compass(compass)?,
            },
        })
    }

    #[getter]
    fn statement(&self) -> Option<&str> {
        self.inner.statement.as_deref()
    }

    #[getter]
    fn index(&self) -> Option<usize> {
        self.inner.id.map(|index| index.0)
    }

    #[getter]
    fn payload(&self) -> Option<Vec<u8>> {
        self.inner.payload.clone()
    }

    #[getter]
    fn port_label(&self) -> Option<&str> {
        self.inner.port_label.as_deref()
    }

    #[getter]
    fn compass(&self) -> Option<&'static str> {
        self.inner.compasspt.map(compass_name)
    }
}

macro_rules! value_class {
    ($doc:literal, $rust:ident, $python:literal, $drawing:ident, $fields:ident) => {
        #[doc = $doc]
        #[gen_stub_pyclass]
        #[pyclass(unsendable, name = $python)]
        pub struct $rust {
            pub(crate) data: Option<Py<PyAny>>,
            pub(crate) drawing: Option<Py<PyDict>>,
        }

        impl $rust {
            pub(crate) fn snapshot(
                py: Python<'_>,
                data: &Py<PyAny>,
                drawing: &Py<PyDict>,
            ) -> PyResult<Self> {
                Ok(Self {
                    data: Some(data.clone_ref(py)),
                    drawing: Some(copy_drawing(py, drawing)?),
                })
            }
        }

        #[gen_stub_pymethods]
        #[pymethods]
        impl $rust {
            #[new]
            #[pyo3(signature = (*, data=None, drawing=None))]
            fn new(
                py: Python<'_>,
                #[gen_stub(override_type(type_repr="typing.Any", imports=("typing")))] data: Option<
                    Py<PyAny>,
                >,
                drawing: Option<PyRef<'_, $drawing>>,
            ) -> PyResult<Self> {
                Ok(Self {
                    data: Some(data.unwrap_or_else(|| py.None())),
                    drawing: Some(match drawing {
                        Some(drawing) => drawing.clone_values(py)?,
                        None => drawing_dict(py, $fields, None)?,
                    }),
                })
            }

            #[getter]
            fn data(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                self.data
                    .as_ref()
                    .map(|data| data.clone_ref(py))
                    .ok_or_else(|| PyReferenceError::new_err("detached value has been cleared"))
            }

            #[setter]
            fn set_data(&mut self, value: Py<PyAny>) {
                self.data = Some(value);
            }

            #[getter]
            fn drawing(&self, py: Python<'_>) -> PyResult<$drawing> {
                self.drawing
                    .as_ref()
                    .map(|drawing| $drawing::detached(py, drawing))
                    .ok_or_else(|| PyReferenceError::new_err("detached value has been cleared"))
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
    };
}

value_class!(
    "Detached node data and drawing metadata passed across a `DotCodec` boundary.",
    PyNodeValue,
    "NodeValue",
    PyNodeDrawing,
    NODE_FIELDS
);
value_class!(
    "Detached edge data and drawing metadata passed across a `DotCodec` boundary.",
    PyEdgeValue,
    "EdgeValue",
    PyEdgeDrawing,
    EDGE_FIELDS
);
value_class!(
    "Detached half-edge data and drawing metadata passed across a `DotCodec` boundary.",
    PyHalfEdgeValue,
    "HalfEdgeValue",
    PyHalfEdgeDrawing,
    HEDGE_FIELDS
);

enum CodecKind {
    Custom {
        encode_node: Py<PyAny>,
        decode_node: Py<PyAny>,
        encode_edge: Py<PyAny>,
        decode_edge: Py<PyAny>,
        encode_half_edge: Py<PyAny>,
        decode_half_edge: Py<PyAny>,
    },
    Linnest,
    Topology,
}

/// An explicit mapping between graph values and DOT-representable records.
#[gen_stub_pyclass]
#[pyclass(unsendable, name = "DotCodec")]
pub struct PyDotCodec {
    kind: Option<CodecKind>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyDotCodec {
    #[new]
    #[pyo3(signature = (*, encode_node, decode_node, encode_edge, decode_edge, encode_half_edge, decode_half_edge))]
    fn new(
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Callable[[NodeValue], DotVertexData]", imports=("typing")))]
        encode_node: Py<PyAny>,
        #[gen_stub(override_type(type_repr="typing.Callable[[DotVertexData], NodeValue]", imports=("typing")))]
        decode_node: Py<PyAny>,
        #[gen_stub(override_type(type_repr="typing.Callable[[EdgeValue], DotEdgeData]", imports=("typing")))]
        encode_edge: Py<PyAny>,
        #[gen_stub(override_type(type_repr="typing.Callable[[DotEdgeData], EdgeValue]", imports=("typing")))]
        decode_edge: Py<PyAny>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdgeValue], DotHalfEdgeData]", imports=("typing")))]
        encode_half_edge: Py<PyAny>,
        #[gen_stub(override_type(type_repr="typing.Callable[[DotHalfEdgeData], HalfEdgeValue]", imports=("typing")))]
        decode_half_edge: Py<PyAny>,
    ) -> PyResult<Self> {
        for (name, callback) in [
            ("encode_node", &encode_node),
            ("decode_node", &decode_node),
            ("encode_edge", &encode_edge),
            ("decode_edge", &decode_edge),
            ("encode_half_edge", &encode_half_edge),
            ("decode_half_edge", &decode_half_edge),
        ] {
            if !callback.bind(py).is_callable() {
                return Err(PyTypeError::new_err(format!(
                    "{name} must be a Python callable"
                )));
            }
        }
        Ok(Self {
            kind: Some(CodecKind::Custom {
                encode_node,
                decode_node,
                encode_edge,
                decode_edge,
                encode_half_edge,
                decode_half_edge,
            }),
        })
    }

    #[classmethod]
    fn linnest(_cls: &Bound<'_, PyType>) -> Self {
        Self {
            kind: Some(CodecKind::Linnest),
        }
    }

    #[classmethod]
    fn topology(_cls: &Bound<'_, PyType>) -> Self {
        Self {
            kind: Some(CodecKind::Topology),
        }
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        if let Some(CodecKind::Custom {
            encode_node,
            decode_node,
            encode_edge,
            decode_edge,
            encode_half_edge,
            decode_half_edge,
        }) = &self.kind
        {
            visit.call(encode_node)?;
            visit.call(decode_node)?;
            visit.call(encode_edge)?;
            visit.call(decode_edge)?;
            visit.call(encode_half_edge)?;
            visit.call(decode_half_edge)?;
        }
        Ok(())
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.kind = None;
    }
}

const NODE_PREFIX: &str = "linnet_node_";
const EDGE_PREFIX: &str = "linnet_edge_";
const GLOBAL_METADATA: &str = "linnet_codec_global_metadata";
const NODE_METADATA: &str = "linnet_codec_node_metadata";
const EDGE_METADATA: &str = "linnet_codec_edge_metadata";
const ELEMENT_NAME: &str = "linnet_name";
const SOURCE_HEDGE_DATA: &str = "linnet_codec_source_half_edge";
const SINK_HEDGE_DATA: &str = "linnet_codec_sink_half_edge";

#[derive(Serialize, Deserialize)]
struct GlobalEnvelope {
    graph_name: Option<String>,
    dot_name: String,
    payload: Option<Vec<u8>>,
}

#[derive(Serialize, Deserialize)]
struct NodeEnvelope {
    codec_name: Option<String>,
    payload: Option<Vec<u8>>,
}

#[derive(Serialize, Deserialize)]
struct EdgeEnvelope {
    payload: Option<Vec<u8>>,
    local_statements: BTreeMap<String, String>,
}

#[derive(Serialize, Deserialize)]
struct HalfEdgeEnvelope {
    statement: Option<String>,
    payload: Option<Vec<u8>>,
    port_label: Option<String>,
    compass: Option<String>,
}

fn encode_hex(bytes: &[u8]) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut output = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        output.push(DIGITS[(byte >> 4) as usize] as char);
        output.push(DIGITS[(byte & 0xf) as usize] as char);
    }
    output
}

fn decode_hex(value: &str) -> PyResult<Vec<u8>> {
    if !value.len().is_multiple_of(2) {
        return Err(PyValueError::new_err("invalid hex-encoded DOT metadata"));
    }
    value
        .as_bytes()
        .chunks_exact(2)
        .map(|digits| {
            let digit = |value: u8| match value {
                b'0'..=b'9' => Some(value - b'0'),
                b'a'..=b'f' => Some(value - b'a' + 10),
                b'A'..=b'F' => Some(value - b'A' + 10),
                _ => None,
            };
            Ok((digit(digits[0])
                .ok_or_else(|| PyValueError::new_err("invalid hex-encoded DOT metadata"))?
                << 4)
                | digit(digits[1])
                    .ok_or_else(|| PyValueError::new_err("invalid hex-encoded DOT metadata"))?)
        })
        .collect()
}

fn decode_text(value: &str, kind: &str) -> PyResult<String> {
    String::from_utf8(decode_hex(value)?)
        .map_err(|_| PyValueError::new_err(format!("invalid UTF-8 {kind} metadata")))
}

fn encode_half_edge_envelope(value: &DotHedgeData) -> PyResult<String> {
    let envelope = HalfEdgeEnvelope {
        statement: value.statement.clone(),
        payload: value.payload.clone(),
        port_label: value.port_label.clone(),
        compass: value.compasspt.map(compass_name).map(str::to_owned),
    };
    serde_json::to_vec(&envelope)
        .map(|bytes| encode_hex(&bytes))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

fn restore_half_edge_envelope(value: &mut DotHedgeData, encoded: &str) -> PyResult<()> {
    let envelope: HalfEdgeEnvelope = serde_json::from_slice(&decode_hex(encoded)?)
        .map_err(|error| PyValueError::new_err(format!("invalid half-edge metadata: {error}")))?;
    value.statement = envelope.statement;
    value.payload = envelope.payload;
    value.port_label = envelope.port_label;
    value.compasspt = parse_compass(envelope.compass.as_deref())?;
    Ok(())
}

fn encode_envelope(value: &impl Serialize) -> PyResult<String> {
    serde_json::to_vec(value)
        .map(|bytes| encode_hex(&bytes))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

fn decode_envelope<T: for<'de> Deserialize<'de>>(value: &str, kind: &str) -> PyResult<T> {
    serde_json::from_slice(&decode_hex(value)?)
        .map_err(|error| PyValueError::new_err(format!("invalid {kind} metadata: {error}")))
}

fn encode_drawing(
    py: Python<'_>,
    drawing: &Py<PyDict>,
    prefix: &str,
    statements: &mut BTreeMap<String, String>,
) -> PyResult<()> {
    for (key, value) in drawing.bind(py).iter() {
        let key = key.extract::<String>()?;
        statements.insert(
            format!("{prefix}{key}"),
            encode_hex(crate::typst::encode_dot_native(py, &value)?.as_bytes()),
        );
    }
    Ok(())
}

fn decode_drawing(
    py: Python<'_>,
    statements: &BTreeMap<String, String>,
    prefix: &str,
    fields: &[&str],
) -> PyResult<Py<PyDict>> {
    let values = PyDict::new(py);
    for (key, value) in statements {
        if let Some(key) = key.strip_prefix(prefix) {
            let encoded = decode_hex(value)?;
            let encoded = std::str::from_utf8(&encoded)
                .map_err(|_| PyValueError::new_err("invalid Linnest drawing metadata"))?;
            values.set_item(key, crate::typst::decode_dot_native(py, encoded)?)?;
        }
    }
    drawing_dict(py, fields, Some(&values))
}

impl PyDotCodec {
    fn kind(&self) -> PyResult<&CodecKind> {
        self.kind
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("DOT codec has been cleared"))
    }

    fn encode_node(
        &self,
        py: Python<'_>,
        value: PyNodeValue,
        index: usize,
        name: Option<String>,
    ) -> PyResult<DotVertexData> {
        match self.kind()? {
            CodecKind::Custom { encode_node, .. } => {
                let value = Py::new(py, value)?;
                let result = encode_node.call1(py, (value,))?;
                let result = result.extract::<PyRef<'_, PyDotVertexData>>(py)?;
                Ok(result.inner.clone())
            }
            CodecKind::Linnest => {
                let mut output = DotVertexData {
                    name,
                    index: Some(NodeIndex(index)),
                    payload: None,
                    statements: BTreeMap::new(),
                };
                encode_drawing(
                    py,
                    value.drawing.as_ref().expect("snapshot"),
                    NODE_PREFIX,
                    &mut output.statements,
                )?;
                Ok(output)
            }
            CodecKind::Topology => Ok(DotVertexData {
                name,
                index: Some(NodeIndex(index)),
                payload: None,
                statements: BTreeMap::new(),
            }),
        }
    }

    fn encode_edge(
        &self,
        py: Python<'_>,
        value: PyEdgeValue,
        index: usize,
    ) -> PyResult<DotEdgeData> {
        match self.kind()? {
            CodecKind::Custom { encode_edge, .. } => {
                let result = encode_edge.call1(py, (Py::new(py, value)?,))?;
                let result = result.extract::<PyRef<'_, PyDotEdgeData>>(py)?;
                Ok(result.inner.clone())
            }
            CodecKind::Linnest => {
                let mut output = DotEdgeData {
                    payload: None,
                    statements: BTreeMap::new(),
                    local_statements: BTreeMap::new(),
                    edge_id: Some(EdgeIndex::from(index)),
                };
                encode_drawing(
                    py,
                    value.drawing.as_ref().expect("snapshot"),
                    EDGE_PREFIX,
                    &mut output.statements,
                )?;
                Ok(output)
            }
            CodecKind::Topology => Ok(DotEdgeData {
                payload: None,
                statements: BTreeMap::new(),
                local_statements: BTreeMap::new(),
                edge_id: Some(EdgeIndex::from(index)),
            }),
        }
    }

    fn encode_half_edge(
        &self,
        py: Python<'_>,
        value: PyHalfEdgeValue,
        index: usize,
    ) -> PyResult<DotHedgeData> {
        match self.kind()? {
            CodecKind::Custom {
                encode_half_edge, ..
            } => {
                let result = encode_half_edge.call1(py, (Py::new(py, value)?,))?;
                let result = result.extract::<PyRef<'_, PyDotHalfEdgeData>>(py)?;
                Ok(result.inner.clone())
            }
            CodecKind::Linnest => {
                let drawing = value.drawing.as_ref().expect("snapshot");
                // Mirror literal strings into DOT's structural half-edge
                // fields; the native payload below preserves richer values.
                let statement = drawing
                    .bind(py)
                    .get_item("statement")?
                    .and_then(|value| value.extract::<String>().ok());
                let port_label = drawing
                    .bind(py)
                    .get_item("port_label")?
                    .and_then(|value| value.extract::<String>().ok());
                let compass = drawing
                    .bind(py)
                    .get_item("compass")?
                    .map(|value| crate::typst::compass_keyword(&value))
                    .transpose()?
                    .flatten();
                Ok(DotHedgeData {
                    statement,
                    id: Some(Hedge(index)),
                    payload: Some(
                        crate::typst::encode_dot_native(py, drawing.bind(py).as_any())?
                            .into_bytes(),
                    ),
                    port_label,
                    compasspt: parse_compass(compass.as_deref())?,
                })
            }
            CodecKind::Topology => Ok(DotHedgeData {
                id: Some(Hedge(index)),
                ..Default::default()
            }),
        }
    }

    fn decode_node(&self, py: Python<'_>, value: &DotVertexData) -> PyResult<PyNodeValue> {
        match self.kind()? {
            CodecKind::Custom { decode_node, .. } => {
                let raw = Py::new(
                    py,
                    PyDotVertexData {
                        inner: value.clone(),
                    },
                )?;
                let result = decode_node.call1(py, (raw,))?;
                snapshot_value::<PyNodeValue>(py, result.bind(py))
            }
            CodecKind::Linnest => Ok(PyNodeValue {
                data: Some(py.None()),
                drawing: Some(decode_drawing(
                    py,
                    &value.statements,
                    NODE_PREFIX,
                    NODE_FIELDS,
                )?),
            }),
            CodecKind::Topology => Ok(PyNodeValue {
                data: Some(py.None()),
                drawing: Some(drawing_dict(py, NODE_FIELDS, None)?),
            }),
        }
    }

    fn decode_edge(&self, py: Python<'_>, value: &DotEdgeData) -> PyResult<PyEdgeValue> {
        match self.kind()? {
            CodecKind::Custom { decode_edge, .. } => {
                let raw = Py::new(
                    py,
                    PyDotEdgeData {
                        inner: value.clone(),
                    },
                )?;
                let result = decode_edge.call1(py, (raw,))?;
                snapshot_value::<PyEdgeValue>(py, result.bind(py))
            }
            CodecKind::Linnest => Ok(PyEdgeValue {
                data: Some(py.None()),
                drawing: Some(decode_drawing(
                    py,
                    &value.statements,
                    EDGE_PREFIX,
                    EDGE_FIELDS,
                )?),
            }),
            CodecKind::Topology => Ok(PyEdgeValue {
                data: Some(py.None()),
                drawing: Some(drawing_dict(py, EDGE_FIELDS, None)?),
            }),
        }
    }

    fn decode_half_edge(&self, py: Python<'_>, value: &DotHedgeData) -> PyResult<PyHalfEdgeValue> {
        match self.kind()? {
            CodecKind::Custom {
                decode_half_edge, ..
            } => {
                let raw = Py::new(
                    py,
                    PyDotHalfEdgeData {
                        inner: value.clone(),
                    },
                )?;
                let result = decode_half_edge.call1(py, (raw,))?;
                snapshot_value::<PyHalfEdgeValue>(py, result.bind(py))
            }
            CodecKind::Linnest => {
                if let Some(payload) = &value.payload {
                    let encoded = std::str::from_utf8(payload).map_err(|_| {
                        PyValueError::new_err("invalid Linnest half-edge drawing payload")
                    })?;
                    let decoded = crate::typst::decode_dot_native(py, encoded)?;
                    let decoded = decoded.bind(py).cast::<PyDict>().map_err(|_| {
                        PyValueError::new_err(
                            "Linnest half-edge drawing payload is not a dictionary",
                        )
                    })?;
                    return Ok(PyHalfEdgeValue {
                        data: Some(py.None()),
                        drawing: Some(drawing_dict(py, HEDGE_FIELDS, Some(decoded))?),
                    });
                }
                let values = PyDict::new(py);
                if let Some(statement) = &value.statement {
                    values.set_item("statement", statement)?;
                }
                if let Some(port_label) = &value.port_label {
                    values.set_item("port_label", port_label)?;
                }
                if let Some(compass) = value.compasspt {
                    values.set_item(
                        "compass",
                        crate::typst::compass_value(py, compass_name(compass))?,
                    )?;
                }
                Ok(PyHalfEdgeValue {
                    data: Some(py.None()),
                    drawing: Some(drawing_dict(py, HEDGE_FIELDS, Some(&values))?),
                })
            }
            CodecKind::Topology => Ok(PyHalfEdgeValue {
                data: Some(py.None()),
                drawing: Some(drawing_dict(py, HEDGE_FIELDS, None)?),
            }),
        }
    }
}

fn reject_reserved(
    statements: &BTreeMap<String, String>,
    reserved: &[&str],
    kind: &str,
) -> PyResult<()> {
    if let Some(key) = reserved.iter().find(|key| statements.contains_key(**key)) {
        return Err(PyValueError::new_err(format!(
            "DOT {kind} statement name {key:?} is reserved by the codec transport"
        )));
    }
    if let Some((key, _)) = statements
        .iter()
        .find(|(key, value)| key.ends_with('\\') || value.ends_with('\\'))
    {
        return Err(PyValueError::new_err(format!(
            "DOT {kind} statement {key:?} cannot end in a backslash with the supported DOT parser"
        )));
    }
    Ok(())
}

trait DetachedValue: Sized {
    fn clone_detached(py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<Self>;
}

macro_rules! detached_impl {
    ($rust:ident) => {
        impl DetachedValue for $rust {
            fn clone_detached(py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<Self> {
                let value = value.extract::<PyRef<'_, $rust>>().map_err(|_| {
                    PyTypeError::new_err(concat!("codec callback must return ", stringify!($rust)))
                })?;
                Ok(Self {
                    data: Some(
                        value
                            .data
                            .as_ref()
                            .ok_or_else(|| PyReferenceError::new_err("detached value was cleared"))?
                            .clone_ref(py),
                    ),
                    drawing: Some(copy_drawing(
                        py,
                        value.drawing.as_ref().ok_or_else(|| {
                            PyReferenceError::new_err("detached value was cleared")
                        })?,
                    )?),
                })
            }
        }
    };
}

detached_impl!(PyNodeValue);
detached_impl!(PyEdgeValue);
detached_impl!(PyHalfEdgeValue);

fn snapshot_value<T: DetachedValue>(py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<T> {
    T::clone_detached(py, value)
}

type RawGraph = HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData, NodeStorageVec<DotVertexData>>;
type RawDotGraph = DotGraph<NodeStorageVec<DotVertexData>>;
type RawDotGraphSet = linnet::parser::set::GraphSet<
    DotEdgeData,
    DotVertexData,
    DotHedgeData,
    GlobalData,
    NodeStorageVec<DotVertexData>,
>;

fn decode_parsed(
    py: Python<'_>,
    mut parsed: RawDotGraph,
    codec: Py<PyDotCodec>,
    node_store: PyNodeStore,
) -> PyResult<PyGraph> {
    let parsed_name =
        (!parsed.global_data.name.is_empty()).then(|| parsed.global_data.name.clone());
    let graph_name = if let Some(encoded) = parsed.global_data.statements.remove(GLOBAL_METADATA) {
        let envelope: GlobalEnvelope = decode_envelope(&encoded, "global")?;
        parsed.global_data.name = envelope.dot_name;
        parsed.global_data.payload = envelope.payload;
        envelope.graph_name
    } else {
        parsed_name
    };
    let codec_ref = codec.borrow(py);
    let mut node_names = BTreeSet::new();
    let mut nodes = Vec::with_capacity(parsed.n_nodes());
    for (_, _, raw) in parsed.graph.iter_nodes() {
        let mut raw = raw.clone();
        let name = raw
            .statements
            .remove(ELEMENT_NAME)
            .map(|value| decode_text(&value, "node name"))
            .transpose()?
            .or_else(|| raw.name.clone());
        if let Some(encoded) = raw.statements.remove(NODE_METADATA) {
            let envelope: NodeEnvelope = decode_envelope(&encoded, "node")?;
            raw.name = envelope.codec_name;
            raw.payload = envelope.payload;
        }
        if let Some(name) = &name {
            if !node_names.insert(name.clone()) {
                return Err(PyValueError::new_err(format!(
                    "duplicate imported node name {name:?}"
                )));
            }
        }
        let value = codec_ref.decode_node(py, &raw)?;
        nodes.push(Some(NodeRecord {
            name,
            data: value.data.expect("decoded snapshot"),
            drawing: value.drawing.expect("decoded snapshot"),
        }));
    }
    let mut edge_names = BTreeSet::new();
    let mut edges = Vec::with_capacity(parsed.n_edges());
    for (_, _, raw_edge) in parsed.graph.iter_edges() {
        let mut clean_edge = raw_edge.data.clone();
        clean_edge.statements.remove(SOURCE_HEDGE_DATA);
        clean_edge.statements.remove(SINK_HEDGE_DATA);
        let name = clean_edge
            .statements
            .remove(ELEMENT_NAME)
            .map(|value| decode_text(&value, "edge name"))
            .transpose()?;
        if let Some(encoded) = clean_edge.statements.remove(EDGE_METADATA) {
            let envelope: EdgeEnvelope = decode_envelope(&encoded, "edge")?;
            clean_edge.payload = envelope.payload;
            clean_edge.local_statements = envelope.local_statements;
        }
        if let Some(name) = &name {
            if !edge_names.insert(name.clone()) {
                return Err(PyValueError::new_err(format!(
                    "duplicate imported edge name {name:?}"
                )));
            }
        }
        let value = codec_ref.decode_edge(py, &clean_edge)?;
        edges.push(Some(EdgeRecord {
            name,
            data: value.data.expect("decoded snapshot"),
            drawing: value.drawing.expect("decoded snapshot"),
        }));
    }
    let mut half_edges = Vec::with_capacity(parsed.n_hedges());
    for (hedge, raw_half_edge) in parsed.graph.iter_hedges() {
        let mut raw_half_edge = raw_half_edge.clone();
        let edge = &parsed.graph[parsed.graph[&hedge]];
        let metadata = match parsed.graph.flow(hedge) {
            Flow::Source => SOURCE_HEDGE_DATA,
            Flow::Sink => SINK_HEDGE_DATA,
        };
        if let Some(encoded) = edge.statements.get(metadata) {
            restore_half_edge_envelope(&mut raw_half_edge, encoded)?;
        }
        let value = codec_ref.decode_half_edge(py, &raw_half_edge)?;
        half_edges.push(Some(HalfEdgeRecord {
            data: value.data.expect("decoded snapshot"),
            drawing: value.drawing.expect("decoded snapshot"),
        }));
    }
    let graph = parsed.graph.map(
        |_, index, _| nodes[index.0].take().expect("decoded node"),
        |_, _, _, index, edge| edge.map(|_| edges[index.0].take().expect("decoded edge")),
        |hedge, _| half_edges[hedge.0].take().expect("decoded half-edge"),
    );
    let graph = PyHedgeGraph::Vec(graph)
        .into_node_store(node_store)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    let global_data = PyGlobalData {
        inner: parsed.global_data,
    };
    drop(codec_ref);
    Ok(PyGraph::from_state(GraphState {
        graph,
        name: graph_name,
        global_data,
        codec: Some(codec),
        render_config: crate::typst::default_render_config(py)?,
        revision: 0,
    }))
}

pub(crate) fn decode_dot(
    py: Python<'_>,
    source: &str,
    codec: Py<PyDotCodec>,
    node_store: PyNodeStore,
) -> PyResult<PyGraph> {
    let parsed = RawDotGraph::from_string(source)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    decode_parsed(py, parsed, codec, node_store)
}

pub(crate) fn decode_dot_set(
    py: Python<'_>,
    source: &str,
    codec: Py<PyDotCodec>,
    node_store: PyNodeStore,
) -> PyResult<Vec<PyGraph>> {
    let parsed = RawDotGraphSet::from_string(source)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    parsed
        .into_iter()
        .map(|graph| decode_parsed(py, graph, codec.clone_ref(py), node_store))
        .collect()
}

pub(crate) fn encode_graph(
    py: Python<'_>,
    graph: &PyGraph,
    codec: Option<Py<PyDotCodec>>,
) -> PyResult<String> {
    let (codec, topology, graph_name, global_data) = {
        let state = graph.state.borrow();
        let state = state
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        let codec = codec
            .or_else(|| state.codec.as_ref().map(|codec| codec.clone_ref(py)))
            .ok_or_else(|| {
                PyValueError::new_err("Graph.to_dot() requires an explicit or stored DotCodec")
            })?;
        (
            codec,
            state.copy_graph(py)?,
            state.name.clone(),
            state.global_data.clone(),
        )
    };
    let topology = topology
        .into_node_store(PyNodeStore::Vec)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    let PyHedgeGraph::Vec(topology) = topology else {
        unreachable!("requested Vec node storage")
    };
    let codec = codec.borrow(py);
    let mut raw_nodes = Vec::with_capacity(topology.n_nodes());
    for (index, _, node) in topology.iter_nodes() {
        let value = PyNodeValue::snapshot(py, &node.data, &node.drawing)?;
        let mut raw = codec.encode_node(py, value, index.0, node.name.clone())?;
        reject_reserved(
            &raw.statements,
            &["id", ELEMENT_NAME, NODE_METADATA],
            "node",
        )?;
        if let Some(name) = &node.name {
            raw.statements
                .insert(ELEMENT_NAME.to_owned(), encode_hex(name.as_bytes()));
        } else {
            raw.statements.remove(ELEMENT_NAME);
        }
        if raw.name.is_some() || raw.payload.is_some() {
            raw.statements.insert(
                NODE_METADATA.to_owned(),
                encode_envelope(&NodeEnvelope {
                    codec_name: raw.name.clone(),
                    payload: raw.payload.clone(),
                })?,
            );
        }
        raw.name = None;
        raw_nodes.push(Some(raw));
    }
    let mut raw_edges = Vec::with_capacity(topology.n_edges());
    for (_, index, edge) in topology.iter_edges() {
        let value = PyEdgeValue::snapshot(py, &edge.data.data, &edge.data.drawing)?;
        let mut raw = codec.encode_edge(py, value, index.0)?;
        reject_reserved(
            &raw.statements,
            &[
                "id",
                ELEMENT_NAME,
                EDGE_METADATA,
                SOURCE_HEDGE_DATA,
                SINK_HEDGE_DATA,
            ],
            "edge",
        )?;
        if let Some(name) = &edge.data.name {
            raw.statements
                .insert(ELEMENT_NAME.to_owned(), encode_hex(name.as_bytes()));
        } else {
            raw.statements.remove(ELEMENT_NAME);
        }
        if raw.payload.is_some() || !raw.local_statements.is_empty() {
            raw.statements.insert(
                EDGE_METADATA.to_owned(),
                encode_envelope(&EdgeEnvelope {
                    payload: raw.payload.clone(),
                    local_statements: raw.local_statements.clone(),
                })?,
            );
        }
        raw_edges.push(Some(raw));
    }
    let mut raw_half_edges = Vec::with_capacity(topology.n_hedges());
    for (hedge, half_edge) in topology.iter_hedges() {
        let value = PyHalfEdgeValue::snapshot(py, &half_edge.data, &half_edge.drawing)?;
        raw_half_edges.push(Some(codec.encode_half_edge(py, value, hedge.0)?));
    }
    let mut raw: RawGraph = topology.map(
        |_, index, _| raw_nodes[index.0].take().expect("encoded node"),
        |_, _, _, index, edge| edge.map(|_| raw_edges[index.0].take().expect("encoded edge")),
        |hedge, _| raw_half_edges[hedge.0].take().expect("encoded half-edge"),
    );
    let edges = raw
        .iter_edges()
        .map(|(pair, index, _)| (pair, index))
        .collect::<Vec<_>>();
    for (pair, edge) in edges {
        let endpoints = match pair {
            HedgePair::Paired { source, sink } => [
                (Some(source), SOURCE_HEDGE_DATA),
                (Some(sink), SINK_HEDGE_DATA),
            ],
            HedgePair::Unpaired { hedge, flow } => match flow {
                Flow::Source => [(Some(hedge), SOURCE_HEDGE_DATA), (None, SINK_HEDGE_DATA)],
                Flow::Sink => [(None, SOURCE_HEDGE_DATA), (Some(hedge), SINK_HEDGE_DATA)],
            },
            HedgePair::Split { .. } => {
                return Err(PyValueError::new_err(
                    "a full graph cannot serialize a split edge to DOT",
                ));
            }
        };
        for (hedge, metadata) in endpoints {
            let Some(hedge) = hedge else {
                continue;
            };
            let envelope = encode_half_edge_envelope(&raw[hedge])?;
            raw[edge].statements.insert(metadata.to_owned(), envelope);
            raw[hedge] = DotHedgeData {
                id: raw[hedge].id,
                ..Default::default()
            };
        }
    }
    let mut global = global_data.inner.clone();
    reject_reserved(&global.statements, &[GLOBAL_METADATA], "global")?;
    reject_reserved(&global.edge_statements, &[], "global edge default")?;
    reject_reserved(&global.node_statements, &[], "global node default")?;
    global.statements.insert(
        GLOBAL_METADATA.to_owned(),
        encode_envelope(&GlobalEnvelope {
            graph_name,
            dot_name: global_data.inner.name,
            payload: global_data.inner.payload,
        })?,
    );
    global.name = "linnet".to_owned();
    let mut graph = DotGraph {
        global_data: global,
        graph: raw,
    };
    graph
        .apply_explicit_id_ordering()
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    Ok(graph.debug_dot())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyGlobalData>()?;
    module.add_class::<PyDotVertexData>()?;
    module.add_class::<PyDotEdgeData>()?;
    module.add_class::<PyDotHalfEdgeData>()?;
    module.add_class::<PyNodeValue>()?;
    module.add_class::<PyEdgeValue>()?;
    module.add_class::<PyHalfEdgeValue>()?;
    module.add_class::<PyDotCodec>()?;
    Ok(())
}
