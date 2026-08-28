use std::collections::BTreeSet;

use pyo3::class::gc::{PyTraverseError, PyVisit};
use pyo3::exceptions::{PyReferenceError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyDictMethods};
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use crate::graph::PyGraph;

pub(crate) const NODE_FIELDS: &[&str] = &[
    "label",
    "placement",
    "shift",
    "rank",
    "minimum_size",
    "maximum_size",
    "style",
    "label_style",
];

pub(crate) const EDGE_FIELDS: &[&str] = &[
    "label",
    "placement",
    "label_position",
    "label_offset",
    "label_angle",
    "bend",
    "routing",
    "minimum_length",
    "same_rank",
    "style",
    "label_style",
    "decoration",
];

pub(crate) const HEDGE_FIELDS: &[&str] = &[
    "label",
    "statement",
    "port_label",
    "compass",
    "anchor",
    "routing",
    "style",
];

#[derive(Clone, Copy)]
pub(crate) enum DrawingKind {
    Node,
    Edge,
    HalfEdge,
}

impl DrawingKind {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Node => "node",
            Self::Edge => "edge",
            Self::HalfEdge => "half-edge",
        }
    }

    pub(crate) fn transport_key(self, key: &str) -> &str {
        match (self, key) {
            (Self::Node | Self::Edge, "placement") => "pos",
            (Self::Node, "minimum_size") => "minimum-size",
            (Self::Node, "maximum_size") => "maximum-size",
            (Self::Node, "style") => "node-style",
            (Self::Node, "label_style") => "node-label-style",
            (Self::Edge, "label_position") => "label-pos",
            (Self::Edge, "label_offset") => "label-offset",
            (Self::Edge, "label_angle") => "label-angle",
            (Self::Edge, "minimum_length") => "minimum-length",
            (Self::Edge, "same_rank") => "same-rank",
            (Self::Edge, "style") => "edge-style",
            (Self::Edge, "label_style") => "edge-label-style",
            (Self::HalfEdge, "port_label") => "port-label",
            _ => key,
        }
    }
}

fn drawing_kind(fields: &[&str]) -> DrawingKind {
    if fields == NODE_FIELDS {
        DrawingKind::Node
    } else if fields == EDGE_FIELDS {
        DrawingKind::Edge
    } else {
        DrawingKind::HalfEdge
    }
}

pub(crate) struct DrawingGuard {
    graph: Py<PyGraph>,
    revision: u64,
}

impl DrawingGuard {
    pub(crate) fn new(graph: Py<PyGraph>, revision: u64) -> Self {
        Self { graph, revision }
    }

    fn check(&self, py: Python<'_>) -> PyResult<()> {
        self.graph.borrow(py).check_revision(self.revision)
    }
}

fn validate_extensions(extensions: &Bound<'_, PyDict>, fields: &[&str]) -> PyResult<()> {
    let kind = drawing_kind(fields);
    let mut reserved = fields.iter().copied().collect::<BTreeSet<_>>();
    reserved.extend(fields.iter().map(|field| kind.transport_key(field)));
    match kind {
        DrawingKind::Node => {
            reserved.insert("statements");
        }
        DrawingKind::Edge => {
            reserved.extend(["shift", "statements"]);
        }
        DrawingKind::HalfEdge => {}
    }
    for (key, value) in extensions.iter() {
        let key = key
            .extract::<String>()
            .map_err(|_| PyValueError::new_err("drawing extension keys must be strings"))?;
        if reserved.contains(key.as_str()) || key == "extensions" {
            return Err(PyValueError::new_err(format!(
                "drawing extension name {key:?} is reserved"
            )));
        }
        crate::typst::validate_native(&value)?;
    }
    Ok(())
}

fn copy_extensions(py: Python<'_>, extensions: &Bound<'_, PyDict>) -> PyResult<Py<PyDict>> {
    let copied = PyDict::new(py);
    for (key, value) in extensions.iter() {
        copied.set_item(key, crate::typst::copy_native(py, &value.unbind())?)?;
    }
    Ok(copied.unbind())
}

pub(crate) fn drawing_dict(
    py: Python<'_>,
    fields: &[&str],
    values: Option<&Bound<'_, PyDict>>,
) -> PyResult<Py<PyDict>> {
    let output = PyDict::new(py);
    let Some(values) = values else {
        return Ok(output.unbind());
    };
    let allowed = fields.iter().copied().collect::<BTreeSet<_>>();
    for (key, value) in values.iter() {
        let key = key
            .extract::<String>()
            .map_err(|_| PyValueError::new_err("drawing field names must be strings"))?;
        if key == "extensions" {
            let extensions = value.cast::<PyDict>().map_err(|_| {
                PyValueError::new_err("extensions must be a string-keyed dictionary")
            })?;
            validate_extensions(extensions, fields)?;
            output.set_item("extensions", copy_extensions(py, extensions)?)?;
        } else if allowed.contains(key.as_str()) {
            crate::typst::validate_drawing_field(drawing_kind(fields).name(), &key, &value)?;
            output.set_item(key, crate::typst::copy_native(py, &value.unbind())?)?;
        } else {
            return Err(PyValueError::new_err(format!(
                "unknown drawing field {key:?}; template-specific fields belong in extensions"
            )));
        }
    }
    Ok(output.unbind())
}

pub(crate) fn copy_drawing(py: Python<'_>, values: &Py<PyDict>) -> PyResult<Py<PyDict>> {
    let copied = PyDict::new(py);
    for (key, value) in values.bind(py).iter() {
        copied.set_item(key, crate::typst::copy_native(py, &value.unbind())?)?;
    }
    Ok(copied.unbind())
}

macro_rules! drawing_class {
    (
        $doc:literal,
        $rust:ident,
        $python:literal,
        $fields:ident,
        [$($getter:ident, $setter:ident, $property:ident, $name:literal, $stub_type:literal);+ $(;)?]
    ) => {
        #[doc = $doc]
        #[gen_stub_pyclass]
        #[pyclass(unsendable, name = $python)]
        pub struct $rust {
            values: Option<Py<PyDict>>,
            guard: Option<DrawingGuard>,
        }

        impl $rust {
            pub(crate) fn live(
                py: Python<'_>,
                values: &Py<PyDict>,
                graph: Py<PyGraph>,
                revision: u64,
            ) -> Self {
                Self {
                    values: Some(values.clone_ref(py)),
                    guard: Some(DrawingGuard::new(graph, revision)),
                }
            }

            pub(crate) fn detached(py: Python<'_>, values: &Py<PyDict>) -> Self {
                Self {
                    values: Some(values.clone_ref(py)),
                    guard: None,
                }
            }

            pub(crate) fn values<'py>(&self, py: Python<'py>) -> PyResult<&Bound<'py, PyDict>> {
                self.check(py)?;
                self.values
                    .as_ref()
                    .map(|values| values.bind(py))
                    .ok_or_else(|| PyReferenceError::new_err("drawing view has been cleared"))
            }

            fn check(&self, py: Python<'_>) -> PyResult<()> {
                if let Some(guard) = &self.guard {
                    guard.check(py)?;
                }
                if self.values.is_none() {
                    return Err(PyReferenceError::new_err("drawing view has been cleared"));
                }
                Ok(())
            }

            fn value(&self, py: Python<'_>, name: &str) -> PyResult<Py<PyAny>> {
                self.values(py)?
                    .get_item(name)?
                    .map(|value| crate::typst::copy_native(py, &value.unbind()))
                    .transpose()
                    .map(|value| value.unwrap_or_else(|| crate::typst::inherit(py)))
            }

            fn set_value(
                &self,
                py: Python<'_>,
                name: &str,
                value: &Bound<'_, PyAny>,
            ) -> PyResult<()> {
                crate::typst::validate_drawing_field(drawing_kind($fields).name(), name, value)?;
                self.values(py)?
                    .set_item(name, crate::typst::copy_native(py, &value.clone().unbind())?)
            }

            pub(crate) fn clone_values(&self, py: Python<'_>) -> PyResult<Py<PyDict>> {
                copy_drawing(py, self.values(py)?.as_unbound())
            }
        }

        #[gen_stub_pymethods]
        #[pymethods]
        impl $rust {
            #[new]
            #[pyo3(signature = (**values))]
            #[gen_stub(skip)]
            fn new(py: Python<'_>, values: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
                Ok(Self {
                    values: Some(drawing_dict(py, $fields, values)?),
                    guard: None,
                })
            }

            $(
                #[getter($property)]
                #[gen_stub(override_return_type(type_repr=$stub_type))]
                fn $getter(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                    self.value(py, $name)
                }

                #[setter($property)]
                fn $setter(
                    slf: &Bound<'_, Self>,
                    #[gen_stub(override_type(type_repr=$stub_type))]
                    value: &Bound<'_, PyAny>,
                ) -> PyResult<()> {
                    slf.borrow().set_value(slf.py(), $name, value)
                }
            )+

            #[getter]
            #[gen_stub(override_return_type(type_repr="builtins.dict[builtins.str, _NativeValue]", imports=("builtins")))]
            fn extensions(&self, py: Python<'_>) -> PyResult<Py<PyDict>> {
                let output = self
                    .values(py)?
                    .get_item("extensions")?
                    .map(|value| value.cast_into::<PyDict>())
                    .transpose()?
                    .map(|value| copy_extensions(py, &value))
                    .transpose()?
                    .unwrap_or_else(|| PyDict::new(py).unbind());
                Ok(output)
            }

            #[setter]
            fn set_extensions(
                slf: &Bound<'_, Self>,
                #[gen_stub(override_type(type_repr="builtins.dict[builtins.str, _NativeValue]", imports=("builtins")))]
                extensions: &Bound<'_, PyDict>,
            ) -> PyResult<()> {
                validate_extensions(extensions, $fields)?;
                slf.borrow()
                    .values(slf.py())?
                    .set_item("extensions", copy_extensions(slf.py(), extensions)?)
            }

            fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
                Ok(format!("{}({})", $python, self.values(py)?.repr()?))
            }

            #[gen_stub(skip)]
            fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
                visit.call(&self.values)?;
                if let Some(guard) = &self.guard {
                    visit.call(&guard.graph)?;
                }
                Ok(())
            }

            #[gen_stub(skip)]
            fn __clear__(&mut self) {
                self.values = None;
                self.guard = None;
            }
        }
    };
}

drawing_class!(
    "Mutable typed drawing metadata for a node.",
    PyNodeDrawing,
    "NodeDrawing",
    NODE_FIELDS,
    [
        label, set_label, label, "label", "_OptionalStaticContent";
        placement, set_placement, placement, "placement", "_PlacementValue";
        shift, set_shift, shift, "shift", "_DrawingPoint";
        rank, set_rank, rank, "rank", "_OptionalInteger";
        minimum_size, set_minimum_size, minimum_size, "minimum_size", "_OptionalNumber";
        maximum_size, set_maximum_size, maximum_size, "maximum_size", "_OptionalNumber";
        style, set_style, style, "style", "_OptionalStyle";
        label_style, set_label_style, label_style, "label_style", "_OptionalStyle";
    ]
);

drawing_class!(
    "Mutable typed drawing metadata for an edge.",
    PyEdgeDrawing,
    "EdgeDrawing",
    EDGE_FIELDS,
    [
        label, set_label, label, "label", "_OptionalStaticContent";
        placement, set_placement, placement, "placement", "_PlacementValue";
        label_position, set_label_position, label_position, "label_position", "_DrawingPoint";
        label_offset, set_label_offset, label_offset, "label_offset", "_OptionalNumber";
        label_angle, set_label_angle, label_angle, "label_angle", "_DrawingAngle";
        bend, set_bend, bend, "bend", "_DrawingAngle";
        routing, set_routing, routing, "routing", "_RoutingValue";
        minimum_length, set_minimum_length, minimum_length, "minimum_length", "_OptionalInteger";
        same_rank, set_same_rank, same_rank, "same_rank", "_OptionalBoolean";
        style, set_style, style, "style", "_OptionalStyleLayers";
        label_style, set_label_style, label_style, "label_style", "_OptionalStyle";
        decoration, set_decoration, decoration, "decoration", "_DrawingDecoration";
    ]
);

drawing_class!(
    "Mutable typed drawing metadata for one endpoint of an edge.",
    PyHalfEdgeDrawing,
    "HalfEdgeDrawing",
    HEDGE_FIELDS,
    [
        label, set_label, label, "label", "_OptionalStaticContent";
        statement, set_statement, statement, "statement", "_DrawingString";
        port_label, set_port_label, port_label, "port_label", "_DrawingString";
        compass, set_compass, compass, "compass", "_CompassValue";
        get_anchor, set_anchor, anchor, "anchor", "_AnchorValue";
        routing, set_routing, routing, "routing", "_RoutingValue";
        style, set_style, style, "style", "_OptionalStyleLayers";
    ]
);

pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        import typing

        class PyNodeDrawing:
            def __new__(cls, *, label: _OptionalStaticContent = ..., placement: _PlacementValue = ..., shift: _DrawingPoint = ..., rank: _OptionalInteger = ..., minimum_size: _OptionalNumber = ..., maximum_size: _OptionalNumber = ..., style: _OptionalStyle = ..., label_style: _OptionalStyle = ..., extensions: _NativeDict = ...) -> NodeDrawing: ...
    "# }
}

pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        import typing

        class PyEdgeDrawing:
            def __new__(cls, *, label: _OptionalStaticContent = ..., placement: _PlacementValue = ..., label_position: _DrawingPoint = ..., label_offset: _OptionalNumber = ..., label_angle: _DrawingAngle = ..., bend: _DrawingAngle = ..., routing: _RoutingValue = ..., minimum_length: _OptionalInteger = ..., same_rank: _OptionalBoolean = ..., style: _OptionalStyleLayers = ..., label_style: _OptionalStyle = ..., decoration: _DrawingDecoration = ..., extensions: _NativeDict = ...) -> EdgeDrawing: ...
    "# }
}

pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        import typing

        class PyHalfEdgeDrawing:
            def __new__(cls, *, label: _OptionalStaticContent = ..., statement: _DrawingString = ..., port_label: _DrawingString = ..., compass: _CompassValue = ..., anchor: _AnchorValue = ..., routing: _RoutingValue = ..., style: _OptionalStyleLayers = ..., extensions: _NativeDict = ...) -> HalfEdgeDrawing: ...
    "# }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyNodeDrawing>()?;
    module.add_class::<PyEdgeDrawing>()?;
    module.add_class::<PyHalfEdgeDrawing>()?;
    Ok(())
}
