use std::env;
use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};

use clinnet::{PreparedTypstRender, TypstModule, TypstRenderRequest, TypstRenderer};
use linnet::half_edge::involution::{Flow, Hedge};
use pyo3::exceptions::{PyReferenceError, PyRuntimeError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyBytes, PyDict, PyDictMethods, PyList, PyListMethods, PyModule};

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
    for (_, _, node) in state.graph.iter_nodes() {
        nodes.append(element_drawing(
            py,
            node.drawing.bind(py),
            DrawingKind::Node,
        )?)?;
    }
    let edges = PyList::empty(py);
    let half_edges = PyList::empty(py);
    for (_, _, edge) in state.graph.iter_edges() {
        edges.append(element_drawing(
            py,
            edge.data.drawing.bind(py),
            DrawingKind::Edge,
        )?)?;
    }
    for index in 0..state.graph.n_hedges() {
        half_edges.append(element_drawing(
            py,
            state.graph[Hedge(index)].drawing.bind(py),
            DrawingKind::HalfEdge,
        )?)?;
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
        let Some(value) = evaluate_selector(py, selector, view?.bind(py), kind)? else {
            continue;
        };
        let patch = element_drawing(py, value.bind(py), kind)?;
        let element = elements.get_item(index)?.cast_into::<PyDict>()?;
        for (key, value) in patch.iter() {
            if element.get_item(&key)?.is_none() {
                element.set_item(key, value)?;
            }
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
            state.graph.n_nodes(),
            state.graph.n_edges(),
            state.graph.n_hedges(),
            (0..state.graph.n_hedges())
                .map(|index| state.graph.flow(Hedge(index)) == Flow::Source)
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
            let Some(value) =
                evaluate_selector(py, selector, view.bind(py), DrawingKind::HalfEdge)?
            else {
                continue;
            };
            let patch = element_drawing(py, value.bind(py), DrawingKind::HalfEdge)?;
            let element = half_edges.get_item(index)?.cast_into::<PyDict>()?;
            for (key, value) in patch.iter() {
                if element.get_item(&key)?.is_none() {
                    element.set_item(key, value)?;
                }
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
    let renderer = TypstRenderer::new(build_dir.path());
    let mut request = TypstRenderRequest::new(dot, transport.config);
    if let Some(template) = transport.template {
        request = request.template(template);
    }
    for import in transport.imports {
        request = request.module(clinnet_module(import)?);
    }
    Ok((renderer, request, build_dir))
}

fn output_format(output: &Path) -> PyResult<&'static str> {
    match output
        .extension()
        .and_then(OsStr::to_str)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        Some("pdf") => Ok("pdf"),
        Some("svg") => Ok("svg"),
        Some("png") => Ok("png"),
        _ => Err(PyRuntimeError::new_err(format!(
            "unsupported Typst output {}; expected a .pdf, .svg, or .png suffix",
            output.display()
        ))),
    }
}

fn compile_typst<'py>(
    py: Python<'py>,
    prepared: &PreparedTypstRender,
    output: Option<&Path>,
    format: &str,
) -> PyResult<Bound<'py, PyAny>> {
    let kwargs = PyDict::new(py);
    kwargs.set_item("input", prepared.entrypoint().to_path_buf())?;
    kwargs.set_item("root", prepared.root().to_path_buf())?;
    kwargs.set_item("format", format)?;
    if let Some(output) = output {
        kwargs.set_item("output", output.to_path_buf())?;
    }
    if let Some(path) = env::var_os("TYPST_PACKAGE_PATH") {
        kwargs.set_item("package_path", PathBuf::from(path))?;
    }
    if let Some(path) = env::var_os("TYPST_PACKAGE_CACHE_PATH") {
        kwargs.set_item("package_cache_path", PathBuf::from(path))?;
    }
    if let Some(paths) = env::var_os("TYPST_FONT_PATHS") {
        kwargs.set_item("font_paths", env::split_paths(&paths).collect::<Vec<_>>())?;
    }
    PyModule::import(py, "typst")?
        .getattr("compile")?
        .call((), Some(&kwargs))
}

pub(crate) fn render_graph(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    output: PathBuf,
    config: Option<&Bound<'_, PyAny>>,
) -> PyResult<PathBuf> {
    let (renderer, request, build_dir) = request(py, graph, config)?;
    let format = output_format(&output)?;
    if let Some(parent) = output
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent).map_err(|error| {
            PyRuntimeError::new_err(format!(
                "failed to create output directory {}: {error}",
                parent.display()
            ))
        })?;
    }
    let prepared = renderer
        .prepare(&request)
        .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
    let _build_dir = build_dir;
    compile_typst(py, &prepared, Some(&output), format)?;
    Ok(output)
}

pub(crate) fn graph_to_svg(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    config: Option<&Bound<'_, PyAny>>,
) -> PyResult<String> {
    let (renderer, request, build_dir) = request(py, graph, config)?;
    let prepared = renderer
        .prepare(&request)
        .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
    let _build_dir = build_dir;
    let rendered = compile_typst(py, &prepared, None, "svg")?;
    let bytes = if let Ok(bytes) = rendered.cast::<PyBytes>() {
        bytes.as_bytes().to_vec()
    } else if let Ok(pages) = rendered.cast::<PyList>() {
        if pages.len() != 1 {
            return Err(PyRuntimeError::new_err(format!(
                "Typst SVG render produced {} pages; expected exactly one",
                pages.len()
            )));
        }
        pages
            .get_item(0)?
            .cast_into::<PyBytes>()?
            .as_bytes()
            .to_vec()
    } else {
        return Err(PyRuntimeError::new_err(
            "typst.compile(format='svg') did not return SVG bytes",
        ));
    };
    String::from_utf8(bytes).map_err(|error| {
        PyRuntimeError::new_err(format!("Typst returned invalid UTF-8 SVG: {error}"))
    })
}
