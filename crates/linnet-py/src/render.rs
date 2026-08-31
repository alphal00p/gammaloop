use std::collections::BTreeMap;
use std::env;
use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};

use linnest::{
    encode_graph_spec_bytes, TypstEdgeSpec, TypstEndpointSpec, TypstGraphSpec, TypstNodeSpec,
};
use linnet::half_edge::involution::{Flow, Hedge, HedgePair, Orientation};
use pyo3::exceptions::{PyReferenceError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyBytes, PyDict, PyDictMethods, PyList, PyListMethods, PyModule};
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use rust_embed::RustEmbed;
use walkdir::WalkDir;

use crate::drawing::DrawingKind;
use crate::graph::{PyEdge, PyGraph, PyHalfEdge, PyNode};
use crate::native_graph::PyHedgeGraph;
use crate::typst::{
    evaluate_selector, render_config_transport, typst_string, RenderConfigTransport,
    SelectorCallbacks, TypstModuleSource,
};

const LINNEST_PACKAGE_DIR: &str = "crates/linnest/typst";
const KURVST_PACKAGE_DIR: &str = "crates/kurvst/typst";
const TYPST_PACKAGES_DIR: &str = "typst-packages";
const USER_SOURCES_DIR: &str = "user-sources";
const DEFAULT_TEMPLATE: &str = "crates/linnest/typst/src/render/figure.typ";

#[derive(RustEmbed)]
#[folder = "$CARGO_MANIFEST_DIR/../linnest/typst"]
#[include = "src/*.typ"]
#[include = "src/**/*.typ"]
#[include = "typst.toml"]
#[include = "linnest.wasm"]
#[include = "LICENSE"]
struct EmbeddedLinnestPackage;

#[derive(RustEmbed)]
#[folder = "$CARGO_MANIFEST_DIR/../kurvst/typst"]
#[include = "src/*.typ"]
#[include = "src/**/*.typ"]
#[include = "typst.toml"]
#[include = "kurvst.wasm"]
#[include = "LICENSE"]
struct EmbeddedKurvstPackage;

#[derive(RustEmbed)]
#[folder = "$CARGO_MANIFEST_DIR/vendor/typst-packages"]
struct EmbeddedTypstPackages;

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
            if element
                .get_item(&key)?
                .as_ref()
                .is_none_or(|value| crate::typst::is_inherit(value))
            {
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
                if element
                    .get_item(&key)?
                    .as_ref()
                    .is_none_or(|value| crate::typst::is_inherit(value))
                {
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

fn request(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    overlay: Option<&Bound<'_, PyAny>>,
) -> PyResult<(Vec<u8>, RenderConfigTransport)> {
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
    let topology = topology_spec(&graph.borrow(py))?;
    Ok((topology, transport))
}

fn topology_endpoint(graph: &PyHedgeGraph, hedge: Hedge) -> TypstEndpointSpec {
    TypstEndpointSpec {
        node: graph.node_id(hedge).0,
        statement: None,
        id: Some(hedge.0),
        data: None,
        port_label: None,
        compass: None,
        in_subgraph: false,
    }
}

fn topology_spec(graph: &PyGraph) -> PyResult<Vec<u8>> {
    let state = graph.state.borrow();
    let state = state
        .as_ref()
        .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
    let nodes = state
        .graph
        .iter_nodes()
        .map(|(index, _, node)| TypstNodeSpec {
            name: node.name.clone(),
            index: Some(index.0),
            data: None,
            pos: None,
            statements: BTreeMap::new(),
        })
        .collect();
    let mut edges = Vec::with_capacity(state.graph.n_edges());
    for (pair, index, edge) in state.graph.iter_edges() {
        let (source, sink, flow) = match pair {
            HedgePair::Paired { source, sink } => (
                Some(topology_endpoint(&state.graph, source)),
                Some(topology_endpoint(&state.graph, sink)),
                None,
            ),
            HedgePair::Unpaired { hedge, flow } => match flow {
                Flow::Source => (
                    Some(topology_endpoint(&state.graph, hedge)),
                    None,
                    Some("source".to_owned()),
                ),
                Flow::Sink => (
                    None,
                    Some(topology_endpoint(&state.graph, hedge)),
                    Some("sink".to_owned()),
                ),
            },
            HedgePair::Split { .. } => {
                return Err(PyValueError::new_err(
                    "a full graph cannot render a split edge",
                ));
            }
        };
        let orientation = match edge.orientation {
            Orientation::Default => "default",
            Orientation::Reversed => "reversed",
            Orientation::Undirected => "undirected",
        };
        edges.push(TypstEdgeSpec {
            name: edge.data.name.clone(),
            source,
            sink,
            data: None,
            orientation: Some(orientation.to_owned()),
            flow,
            id: Some(index.0),
            pos: None,
            statements: BTreeMap::new(),
        });
    }
    // GlobalData is a DOT-codec concern. Rendering stages only topology;
    // typed drawing state travels separately in the V1 configuration.
    let spec = TypstGraphSpec {
        name: Some(state.name.clone().unwrap_or_else(|| "linnet".to_owned())),
        data: None,
        statements: BTreeMap::new(),
        default_edge_statements: BTreeMap::new(),
        default_node_statements: BTreeMap::new(),
        nodes,
        edges,
    };
    encode_graph_spec_bytes(&spec).map_err(PyRuntimeError::new_err)
}

fn prepare(topology: Vec<u8>, transport: RenderConfigTransport) -> PyResult<PreparedRender> {
    let build_dir =
        tempfile::tempdir().map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
    stage_default_assets(build_dir.path())?;
    let build_root = canonicalize(build_dir.path(), "render directory")?;
    let package_store = canonicalize(
        &build_root.join(TYPST_PACKAGES_DIR),
        "bundled Typst package store",
    )?;
    for (variable, description) in [
        ("TYPST_PACKAGE_CACHE_PATH", "Typst package cache"),
        ("TYPST_PACKAGE_PATH", "Typst package path"),
    ] {
        if let Some(path) = env::var_os(variable) {
            let path = Path::new(&path);
            if variable == "TYPST_PACKAGE_CACHE_PATH" && !path.exists() {
                continue;
            }
            let path = canonicalize(path, description)?;
            copy_directory(&path, &package_store, description)?;
        }
    }
    let topology_path = build_root.join("diagram.cbor");
    fs::write(&topology_path, topology).map_err(|error| {
        PyRuntimeError::new_err(format!(
            "failed to stage topology {}: {error}",
            topology_path.display()
        ))
    })?;

    let mut source_paths = transport.template.iter().cloned().collect::<Vec<_>>();
    source_paths.extend(transport.imports.iter().filter_map(|import| {
        if let TypstModuleSource::File(path) = &import.source {
            Some(path.clone())
        } else {
            None
        }
    }));
    let mut staged_sources =
        stage_user_sources(&build_root, &source_paths, transport.source_root.as_deref())?
            .into_iter();
    let template = if transport.template.is_some() {
        staged_sources
            .next()
            .ok_or_else(|| PyRuntimeError::new_err("failed to stage Typst template"))?
    } else {
        canonicalize(
            &build_root.join(DEFAULT_TEMPLATE),
            "default Linnest template",
        )?
    };
    let module_files = transport
        .imports
        .iter()
        .map(|import| match import.source {
            TypstModuleSource::File(_) => staged_sources
                .next()
                .map(Some)
                .ok_or_else(|| PyRuntimeError::new_err("failed to stage Typst module")),
            TypstModuleSource::Package(_) => Ok(None),
        })
        .collect::<PyResult<Vec<_>>>()?;
    let entrypoint = build_root.join("render.typ");
    let source = entrypoint_source(
        &transport,
        &template,
        &module_files,
        &topology_path,
        &build_root,
    )?;
    fs::write(&entrypoint, source).map_err(|error| {
        PyRuntimeError::new_err(format!(
            "failed to stage Typst entrypoint {}: {error}",
            entrypoint.display()
        ))
    })?;
    Ok(PreparedRender {
        _build_dir: build_dir,
        entrypoint,
        root: build_root,
        package_store,
    })
}

fn stage_default_assets(root: &Path) -> PyResult<()> {
    write_embedded_assets::<EmbeddedLinnestPackage>(&root.join(LINNEST_PACKAGE_DIR))?;
    write_embedded_assets::<EmbeddedKurvstPackage>(&root.join(KURVST_PACKAGE_DIR))?;
    write_embedded_assets::<EmbeddedTypstPackages>(&root.join(TYPST_PACKAGES_DIR))
}

fn write_embedded_assets<E: RustEmbed>(root: &Path) -> PyResult<()> {
    for path in E::iter() {
        let target = root.join(path.as_ref());
        let contents = E::get(path.as_ref()).ok_or_else(|| {
            PyRuntimeError::new_err(format!("embedded render asset {path} is missing"))
        })?;
        if let Some(parent) = target.parent() {
            fs::create_dir_all(parent).map_err(|error| {
                PyRuntimeError::new_err(format!(
                    "failed to create render asset directory {}: {error}",
                    parent.display()
                ))
            })?;
        }
        fs::write(&target, contents.data.as_ref()).map_err(|error| {
            PyRuntimeError::new_err(format!(
                "failed to stage render asset {}: {error}",
                target.display()
            ))
        })?;
    }
    Ok(())
}

fn canonicalize(path: &Path, description: &str) -> PyResult<PathBuf> {
    fs::canonicalize(path).map_err(|error| {
        PyRuntimeError::new_err(format!(
            "failed to resolve {description} {}: {error}",
            path.display()
        ))
    })
}

fn copy_directory(source: &Path, target: &Path, description: &str) -> PyResult<()> {
    fs::create_dir_all(target).map_err(|error| {
        PyRuntimeError::new_err(format!(
            "failed to create staged {description} {}: {error}",
            target.display()
        ))
    })?;
    for entry in WalkDir::new(source)
        .follow_links(true)
        .into_iter()
        .filter_entry(|entry| !entry.path().starts_with(target))
    {
        let entry = entry.map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
        let relative = entry.path().strip_prefix(source).map_err(|error| {
            PyRuntimeError::new_err(format!(
                "failed to stage {description} {}: {error}",
                entry.path().display()
            ))
        })?;
        let destination = target.join(relative);
        if entry.file_type().is_dir() {
            fs::create_dir_all(&destination).map_err(|error| {
                PyRuntimeError::new_err(format!(
                    "failed to create staged {description} {}: {error}",
                    destination.display()
                ))
            })?;
        } else if entry.file_type().is_file() {
            fs::copy(entry.path(), &destination).map_err(|error| {
                PyRuntimeError::new_err(format!(
                    "failed to stage {description} {}: {error}",
                    entry.path().display()
                ))
            })?;
        }
    }
    Ok(())
}

fn stage_user_sources(
    build_root: &Path,
    paths: &[PathBuf],
    configured_root: Option<&Path>,
) -> PyResult<Vec<PathBuf>> {
    let mut roots = Vec::<(PathBuf, PathBuf)>::new();
    let configured_root = configured_root
        .map(|root| canonicalize(root, "Typst source root"))
        .transpose()?;
    if let Some(root) = configured_root.as_ref().filter(|root| !root.is_dir()) {
        return Err(PyRuntimeError::new_err(format!(
            "Typst source root {} is not a directory",
            root.display()
        )));
    }
    paths
        .iter()
        .map(|path| {
            let source = canonicalize(path, "Typst source")?;
            let parent = source.parent().ok_or_else(|| {
                PyRuntimeError::new_err(format!(
                    "Typst source {} has no parent directory",
                    source.display()
                ))
            })?;
            let source_root = configured_root.clone().unwrap_or_else(|| {
                roots
                    .iter()
                    .find(|(root, _)| source.starts_with(root))
                    .map(|(root, _)| root.clone())
                    .unwrap_or_else(|| {
                        parent
                            .ancestors()
                            .find(|ancestor| ancestor.join("typst.toml").is_file())
                            .unwrap_or(parent)
                            .to_path_buf()
                    })
            });
            let relative = source.strip_prefix(&source_root).map_err(|_| {
                PyRuntimeError::new_err(format!(
                    "Typst source {} is outside source root {}",
                    source.display(),
                    source_root.display()
                ))
            })?;
            let target_root = if let Some((_, target)) =
                roots.iter().find(|(existing, _)| existing == &source_root)
            {
                target.clone()
            } else {
                let target = build_root
                    .join(USER_SOURCES_DIR)
                    .join(roots.len().to_string());
                copy_directory(&source_root, &target, "Typst source tree")?;
                roots.push((source_root.clone(), target.clone()));
                target
            };
            Ok(target_root.join(relative))
        })
        .collect()
}

fn typst_root_path(path: &Path, root: &Path) -> PyResult<String> {
    let path = canonicalize(path, "Typst input")?;
    let relative = path.strip_prefix(root).map_err(|_| {
        PyRuntimeError::new_err(format!(
            "Typst input {} is outside project root {}",
            path.display(),
            root.display()
        ))
    })?;
    Ok(typst_string(&format!(
        "/{}",
        relative.to_string_lossy().replace('\\', "/")
    )))
}

fn entrypoint_source(
    transport: &RenderConfigTransport,
    template: &Path,
    module_files: &[Option<PathBuf>],
    topology_path: &Path,
    root: &Path,
) -> PyResult<String> {
    let template = typst_root_path(template, root)?;
    let topology_path = typst_root_path(topology_path, root)?;
    let mut source = format!("#import {template} as _linnet_template\n");
    for (module, file) in transport.imports.iter().zip(module_files) {
        let module_source = match (&module.source, file) {
            (TypstModuleSource::File(_), Some(path)) => typst_root_path(path, root)?,
            (TypstModuleSource::Package(package), None) => typst_string(package),
            _ => {
                return Err(PyRuntimeError::new_err(
                    "internal Typst module source mismatch",
                ));
            }
        };
        source.push_str(&format!("#import {module_source} as {}\n", module.alias));
    }
    source.push_str("\n#let _linnet_config = {\n  let value = (");
    source.push_str(&transport.config_source);
    source.push_str(
        ")\n  if type(value) != dictionary {\n    panic(\"Linnet render config must be a dictionary\")\n  }\n  value + (graph-spec-path: ",
    );
    source.push_str(&topology_path);
    source.push_str(",)\n}\n\n#_linnet_template.render(_linnet_config)\n");
    Ok(source)
}

/// One Typst render whose generated entrypoint, topology, and source modules share a lifetime.
#[gen_stub_pyclass]
#[pyclass(module = "linnet_py", frozen)]
pub(crate) struct PreparedRender {
    _build_dir: tempfile::TempDir,
    entrypoint: PathBuf,
    root: PathBuf,
    package_store: PathBuf,
}

impl PreparedRender {
    fn typst_source_value(&self) -> PyResult<String> {
        fs::read_to_string(&self.entrypoint)
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))
    }

    fn render_to(&self, py: Python<'_>, output: PathBuf) -> PyResult<PathBuf> {
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
        compile_typst(py, self, Some(&output), format)?;
        Ok(output)
    }

    fn svg(&self, py: Python<'_>) -> PyResult<String> {
        let rendered = compile_typst(py, self, None, "svg")?;
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
}

#[gen_stub_pymethods]
#[pymethods]
impl PreparedRender {
    /// Return the exact generated Typst entrypoint for this preparation.
    #[getter]
    fn typst_source(&self) -> PyResult<String> {
        self.typst_source_value()
    }

    /// Compile this preparation to a PDF, SVG, or PNG selected by the suffix.
    fn render(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="builtins.str | os.PathLike[builtins.str]", imports=("builtins", "os")))]
        output: PathBuf,
    ) -> PyResult<PathBuf> {
        self.render_to(py, output)
    }

    /// Compile this preparation and return its one-page SVG document.
    fn to_svg(&self, py: Python<'_>) -> PyResult<String> {
        self.svg(py)
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PreparedRender>()?;
    Ok(())
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
    prepared: &PreparedRender,
    output: Option<&Path>,
    format: &str,
) -> PyResult<Bound<'py, PyAny>> {
    let kwargs = PyDict::new(py);
    kwargs.set_item("input", &prepared.entrypoint)?;
    kwargs.set_item("root", &prepared.root)?;
    kwargs.set_item("format", format)?;
    if let Some(output) = output {
        kwargs.set_item("output", output.to_path_buf())?;
    }
    kwargs.set_item("package_path", &prepared.package_store)?;
    kwargs.set_item("package_cache_path", &prepared.package_store)?;
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
    prepare_graph(py, graph, config)?.render_to(py, output)
}

pub(crate) fn graph_to_svg(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    config: Option<&Bound<'_, PyAny>>,
) -> PyResult<String> {
    prepare_graph(py, graph, config)?.svg(py)
}

pub(crate) fn prepare_graph(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    config: Option<&Bound<'_, PyAny>>,
) -> PyResult<PreparedRender> {
    let (topology, transport) = request(py, graph, config)?;
    prepare(topology, transport)
}
