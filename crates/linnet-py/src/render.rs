use std::path::PathBuf;

use clinnet::{TypstModule, TypstRenderRequest, TypstRenderer};
use pyo3::exceptions::{PyReferenceError, PyRuntimeError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyDictMethods, PyList, PyListMethods};

use crate::drawing::DrawingKind;
use crate::graph::{PyEdge, PyGraph, PyHalfEdge, PyNode};
use crate::typst::{
    evaluate_selector, render_config_transport, RenderConfigTransport, SelectorCallbacks,
    TypstImport, TypstModuleSource,
};

fn element_drawing<'py>(
    py: Python<'py>,
    source: &Bound<'py, PyDict>,
    kind: DrawingKind,
) -> PyResult<Bound<'py, PyDict>> {
    let output = PyDict::new(py);
    for (key, value) in source.iter() {
        let key = key.extract::<String>()?;
        if key == "extensions" {
            let extensions = value.cast::<PyDict>()?;
            for (extension, value) in extensions.iter() {
                output.set_item(extension, crate::typst::copy_native(py, &value.unbind())?)?;
            }
        } else {
            let value = if key == "placement" {
                crate::typst::normalize_placement(py, &value)?
            } else {
                crate::typst::copy_native(py, &value.unbind())?
            };
            output.set_item(kind.transport_key(&key), value)?;
        }
    }
    Ok(output)
}

fn base_elements(py: Python<'_>, graph: &PyGraph) -> PyResult<Py<PyDict>> {
    let state = graph.state.borrow();
    let state = state
        .as_ref()
        .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
    let nodes = PyList::empty(py);
    for node in &state.nodes {
        nodes.append(element_drawing(
            py,
            node.drawing.bind(py),
            DrawingKind::Node,
        )?)?;
    }
    let edges = PyList::empty(py);
    let half_edges = PyList::empty(py);
    for edge in &state.edges {
        edges.append(element_drawing(
            py,
            edge.drawing.bind(py),
            DrawingKind::Edge,
        )?)?;
        for endpoint in [edge.source.as_ref(), edge.sink.as_ref()]
            .into_iter()
            .flatten()
        {
            half_edges.append(element_drawing(
                py,
                endpoint.drawing.bind(py),
                DrawingKind::HalfEdge,
            )?)?;
        }
    }
    let elements = PyDict::new(py);
    elements.set_item("graph", PyDict::new(py))?;
    elements.set_item("nodes", nodes)?;
    elements.set_item("edges", edges)?;
    elements.set_item("hedges", half_edges)?;
    Ok(elements.unbind())
}

fn apply_selector(
    py: Python<'_>,
    selector: Option<&Py<PyAny>>,
    elements: &Bound<'_, PyList>,
    views: impl IntoIterator<Item = PyResult<Py<PyAny>>>,
    kind: DrawingKind,
) -> PyResult<()> {
    let Some(selector) = selector else {
        return Ok(());
    };
    for (index, view) in views.into_iter().enumerate() {
        let value = evaluate_selector(py, selector, view?.bind(py), kind)?;
        let element = elements.get_item(index)?.cast_into::<PyDict>()?;
        if element.get_item("selector-style")?.is_none() {
            element.set_item("selector-style", value)?;
        }
    }
    Ok(())
}

fn apply_selectors(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    elements: &Bound<'_, PyDict>,
    selectors: &SelectorCallbacks,
) -> PyResult<()> {
    let graph_ref = graph.borrow(py);
    let (revision, node_count, edge_count, half_edge_count, half_edge_roles) = {
        let state = graph_ref.state.borrow();
        let state = state
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        (
            state.revision,
            state.nodes.len(),
            state.edges.len(),
            state.n_hedges(),
            state
                .edges
                .iter()
                .flat_map(|edge| {
                    [
                        edge.source.as_ref().map(|_| true),
                        edge.sink.as_ref().map(|_| false),
                    ]
                })
                .flatten()
                .collect::<Vec<_>>(),
        )
    };
    drop(graph_ref);

    let nodes = elements
        .get_item("nodes")?
        .expect("base elements")
        .cast_into::<PyList>()?;
    apply_selector(
        py,
        selectors.node.as_ref(),
        &nodes,
        (0..node_count).map(|index| {
            Ok(Py::new(py, PyNode::new(graph.clone_ref(py), index, revision))?.into_any())
        }),
        DrawingKind::Node,
    )?;

    let edges = elements
        .get_item("edges")?
        .expect("base elements")
        .cast_into::<PyList>()?;
    apply_selector(
        py,
        selectors.edge.as_ref(),
        &edges,
        (0..edge_count).map(|index| {
            Ok(Py::new(py, PyEdge::new(graph.clone_ref(py), index, revision))?.into_any())
        }),
        DrawingKind::Edge,
    )?;

    let half_edges = elements
        .get_item("hedges")?
        .expect("base elements")
        .cast_into::<PyList>()?;
    for (index, is_source) in half_edge_roles
        .into_iter()
        .enumerate()
        .take(half_edge_count)
    {
        let view = Py::new(py, PyHalfEdge::new(graph.clone_ref(py), index, revision))?.into_any();
        let selector = if is_source {
            selectors.source.as_ref()
        } else {
            selectors.sink.as_ref()
        };
        if let Some(selector) = selector {
            let value = evaluate_selector(py, selector, view.bind(py), DrawingKind::HalfEdge)?;
            let element = half_edges.get_item(index)?.cast_into::<PyDict>()?;
            if element.get_item("selector-style")?.is_none() {
                element.set_item("selector-style", value)?;
            }
        }
    }
    graph.borrow(py).check_revision(revision).map_err(|_| {
        PyReferenceError::new_err("drawing selectors must not mutate graph topology")
    })?;
    Ok(())
}

fn clinnet_module(import: TypstImport) -> PyResult<TypstModule> {
    let result = match import.source {
        TypstModuleSource::File(path) => TypstModule::new(import.alias, path),
        TypstModuleSource::Package(package) => TypstModule::package(import.alias, package),
    };
    result.map_err(|error| PyRuntimeError::new_err(error.to_string()))
}

fn request(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    overlay: Option<&Bound<'_, PyAny>>,
) -> PyResult<(TypstRenderer, TypstRenderRequest, tempfile::TempDir)> {
    let graph_ref = graph.borrow(py);
    let elements = base_elements(py, &graph_ref)?;
    let base = graph_ref
        .state
        .borrow()
        .as_ref()
        .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
        .render_config
        .clone_ref(py);
    let initial = render_config_transport(py, &base, overlay, elements.bind(py))?;
    drop(graph_ref);
    apply_selectors(py, graph, elements.bind(py), &initial.selectors)?;
    let transport = render_config_transport(py, &base, overlay, elements.bind(py))?;
    build_request(py, graph, transport)
}

fn build_request(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    transport: RenderConfigTransport,
) -> PyResult<(TypstRenderer, TypstRenderRequest, tempfile::TempDir)> {
    let dot = crate::dot::topology_dot(py, &graph.borrow(py))?;
    let build_dir =
        tempfile::tempdir().map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
    let renderer =
        TypstRenderer::new(build_dir.path()).typst_executable(transport.typst_executable);
    let mut request = TypstRenderRequest::new(dot, transport.config);
    if let Some(template) = transport.template {
        request = request.template(template);
    }
    for import in transport.imports {
        request = request.module(clinnet_module(import)?);
    }
    Ok((renderer, request, build_dir))
}

pub(crate) fn render_graph(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    output: PathBuf,
    config: Option<&Bound<'_, PyAny>>,
) -> PyResult<PathBuf> {
    let (renderer, request, build_dir) = request(py, graph, config)?;
    let rendered = output.clone();
    py.detach(move || {
        let _build_dir = build_dir;
        renderer.render(&request, &rendered)
    })
    .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
    Ok(output)
}

pub(crate) fn graph_to_svg(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    config: Option<&Bound<'_, PyAny>>,
) -> PyResult<String> {
    let (renderer, request, build_dir) = request(py, graph, config)?;
    py.detach(move || {
        let _build_dir = build_dir;
        renderer.to_svg(&request)
    })
    .map_err(|error| PyRuntimeError::new_err(error.to_string()))
}
