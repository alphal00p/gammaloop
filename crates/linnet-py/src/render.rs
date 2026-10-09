use std::collections::BTreeMap;
use std::env;
use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};

use linnest::{
    TypstEdgeSpec, TypstEndpointSpec, TypstGraphSpec, TypstNodeSpec, encode_graph_spec_bytes,
};
use linnet::half_edge::involution::{Flow, Hedge, HedgePair, Orientation};
use linnet::half_edge::subgraph::{Inclusion, SubSetLike};
use pyo3::exceptions::{PyReferenceError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyDictMethods, PyList, PyListMethods, PyModule};
use rust_embed::RustEmbed;
use walkdir::WalkDir;

use crate::drawing::DrawingKind;
use crate::graph::{PyEdge, PyGraph, PyHalfEdge, PyNode};
use crate::native_graph::PyHedgeGraph;
use crate::topology::PySubgraph;
use crate::typst::{
    RenderConfigTransport, SelectorCallbacks, TypstModuleSource, default_render_config,
    evaluate_selector, render_config_transport, typst_string,
};

const LINNEST_PACKAGE_DIR: &str = "crates/linnest/typst";
const KURVST_PACKAGE_DIR: &str = "crates/kurvst/typst";
const TYPST_PACKAGES_DIR: &str = "typst-packages";
const USER_SOURCES_DIR: &str = "user-sources";
const DEFAULT_TEMPLATE: &str = "crates/linnest/typst/src/render/figure.typ";
const ENTRYPOINT: &str = "main.typ";
const TOPOLOGY: &str = "diagram.cbor";
const SUBGRAPH_STYLE: &str = "subgraph.typ";

#[derive(RustEmbed)]
#[folder = "$CARGO_MANIFEST_DIR/../linnest/typst"]
#[include = "src/*.typ"]
#[include = "src/**/*.typ"]
#[include = "typst.toml"]
#[include = "linnest.wasm"]
#[include = "ec-layout.wasm"]
#[include = "LICENSE"]
#[include = "LICENSE.ec-layout"]
#[include = "LICENSE.clarabel"]
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
#[include = "preview/**"]
struct EmbeddedTypstPackages;

impl EmbeddedTypstPackages {
    const CETZ_ARCHIVE: &'static [u8] =
        include_bytes!("../vendor/typst-packages/archives/cetz-0.5.2.tar.gz");
    const CETZ_PACKAGE: &'static str = "preview/cetz/0.5.2";

    /// Remove an owned staging entry without following a symlink to its target.
    fn remove_staged_path(path: &Path) -> std::io::Result<()> {
        match fs::symlink_metadata(path) {
            Ok(metadata) if metadata.is_dir() => fs::remove_dir_all(path),
            Ok(_) => fs::remove_file(path),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
            Err(error) => Err(error),
        }
    }

    fn stage_cetz(package_store: &Path) -> std::io::Result<()> {
        use std::io::{Error, ErrorKind};
        use std::path::Component;

        // Only this compile-time public release archive enters this traversal.
        // Validate every entry before changing the private staged package.
        let decoder = flate2::read::GzDecoder::new(Self::CETZ_ARCHIVE);
        let mut archive = tar::Archive::new(decoder);
        let mut assets = BTreeMap::new();
        for entry in archive.entries()? {
            let mut entry = entry?;
            let path = entry.path()?.into_owned();
            if path
                .components()
                .any(|part| !matches!(part, Component::Normal(_) | Component::CurDir))
            {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    "unsafe CeTZ archive path",
                ));
            }
            let kind = entry.header().entry_type();
            if kind.is_dir() {
                continue;
            }
            if !kind.is_file() || path == Path::new(".") {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    "unsupported CeTZ archive entry",
                ));
            }
            let mut contents = Vec::new();
            std::io::Read::read_to_end(&mut entry, &mut contents)?;
            if assets.insert(path, contents).is_some() {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    "duplicate CeTZ archive path",
                ));
            }
        }
        let root = package_store.join(Self::CETZ_PACKAGE);
        // Remove copied files (including read-only files) rather than overwriting
        // them. The package store is an owned temporary directory, not a cache.
        Self::remove_staged_path(&root)?;
        fs::create_dir_all(&root)?;
        for (path, contents) in assets {
            let target = root.join(path);
            if let Some(parent) = target.parent() {
                fs::create_dir_all(parent)?;
            }
            fs::write(target, contents)?;
        }
        Ok(())
    }
}

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
        route_points: Vec::new(),
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
    // GlobalData is a DOT-codec concern. Rendering transports only topology;
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

fn prepare(
    topology: Vec<u8>,
    transport: RenderConfigTransport,
    selection: Option<(Vec<bool>, Vec<usize>)>,
) -> PyResult<PreparedRender> {
    let mut prepared = PreparedRender::from_sources(BTreeMap::new())?;
    let build_root = &prepared.root;
    write_project_asset(build_root, TOPOLOGY, &topology)?;

    let mut source =
        prepared.configuration_source(&transport, Some(DEFAULT_TEMPLATE), Some(TOPOLOGY))?;
    let files = &mut prepared.files;
    if let Some((hedges, nodes)) = selection {
        files.insert(
            SUBGRAPH_STYLE.to_owned(),
            include_bytes!("../typst/subgraph.typ").to_vec(),
        );
        let hedges = hedges
            .iter()
            .map(|value| format!("{value},"))
            .collect::<String>();
        let nodes = nodes
            .iter()
            .map(|value| format!("{value},"))
            .collect::<String>();
        source.push_str(&format!(
            "\n#import \"/{SUBGRAPH_STYLE}\" as _linnet_subgraph\n\
             #set page(fill: none)\n\
             #_linnet_template.render(_linnet_subgraph.focus(_linnet_config, ({hedges}), ({nodes})))\n"
        ));
    } else {
        source.push_str("\n#_linnet_template.render(_linnet_config)\n");
    }
    files.insert(ENTRYPOINT.to_owned(), source.into_bytes());
    Ok(prepared)
}

fn insert_embedded_assets<E: RustEmbed>(
    files: &mut BTreeMap<String, Vec<u8>>,
    build_root: &Path,
    root: &str,
) -> PyResult<()> {
    for path in E::iter() {
        let contents = E::get(path.as_ref()).ok_or_else(|| {
            PyRuntimeError::new_err(format!("embedded render asset {path} is missing"))
        })?;
        insert_project_asset(
            files,
            build_root,
            &format!("{root}/{}", path.replace('\\', "/")),
            contents.data.as_ref(),
        )?;
    }
    Ok(())
}

fn insert_project_asset(
    files: &mut BTreeMap<String, Vec<u8>>,
    build_root: &Path,
    path: &str,
    contents: &[u8],
) -> PyResult<()> {
    // typst-py's multi-file input decodes every value as UTF-8. Keep binary
    // assets in the project filesystem (MEMFS in Pyodide) at the same paths.
    if std::str::from_utf8(contents).is_ok() {
        files.insert(path.to_owned(), contents.to_vec());
        Ok(())
    } else {
        write_project_asset(build_root, path, contents)
    }
}

fn write_project_asset(build_root: &Path, path: &str, contents: &[u8]) -> PyResult<()> {
    let target = build_root.join(path);
    if let Some(parent) = target.parent() {
        fs::create_dir_all(parent).map_err(|error| {
            PyRuntimeError::new_err(format!(
                "failed to create render asset directory {}: {error}",
                parent.display()
            ))
        })?;
    }
    fs::write(&target, contents).map_err(|error| {
        PyRuntimeError::new_err(format!(
            "failed to stage render asset {}: {error}",
            target.display()
        ))
    })
}

fn write_embedded_assets<E: RustEmbed>(root: &Path) -> PyResult<()> {
    let mut staged = std::collections::BTreeSet::new();
    for path in E::iter() {
        // Bundled identities replace complete packages, not a mixture of
        // embedded files and potentially incompatible external remnants.
        let package: PathBuf = Path::new(path.as_ref()).components().take(3).collect();
        if staged.insert(package.clone()) {
            let destination = root.join(package);
            EmbeddedTypstPackages::remove_staged_path(&destination).map_err(|error| {
                PyRuntimeError::new_err(format!(
                    "failed to replace bundled package {}: {error}",
                    destination.display()
                ))
            })?;
        }
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
        // Overlay authority applies to entry types too. Replace only the
        // private staged destination; neither external source store is changed.
        if destination.exists()
            && (entry.file_type().is_file()
                || (entry.file_type().is_dir() && !destination.is_dir()))
        {
            EmbeddedTypstPackages::remove_staged_path(&destination).map_err(|error| {
                PyRuntimeError::new_err(format!(
                    "failed to replace staged {description} {}: {error}",
                    destination.display()
                ))
            })?;
        }
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

fn collect_user_sources(
    files: &mut BTreeMap<String, Vec<u8>>,
    build_root: &Path,
    paths: &[PathBuf],
    configured_root: Option<&Path>,
) -> PyResult<Vec<String>> {
    let mut roots = Vec::<(PathBuf, String)>::new();
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
                let target = format!("{USER_SOURCES_DIR}/{}", roots.len());
                collect_directory(
                    files,
                    build_root,
                    &source_root,
                    &target,
                    "Typst source tree",
                )?;
                roots.push((source_root.clone(), target.clone()));
                target
            };
            Ok(format!(
                "{target_root}/{}",
                relative.to_string_lossy().replace('\\', "/")
            ))
        })
        .collect()
}

fn collect_directory(
    files: &mut BTreeMap<String, Vec<u8>>,
    build_root: &Path,
    source: &Path,
    target: &str,
    description: &str,
) -> PyResult<()> {
    for entry in WalkDir::new(source).follow_links(true) {
        let entry = entry.map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
        if !entry.file_type().is_file() {
            continue;
        }
        let relative = entry.path().strip_prefix(source).map_err(|error| {
            PyRuntimeError::new_err(format!(
                "failed to collect {description} {}: {error}",
                entry.path().display()
            ))
        })?;
        let path = format!("{target}/{}", relative.to_string_lossy().replace('\\', "/"));
        let contents = fs::read(entry.path()).map_err(|error| {
            PyRuntimeError::new_err(format!(
                "failed to collect {description} {}: {error}",
                entry.path().display()
            ))
        })?;
        insert_project_asset(files, build_root, &path, &contents)?;
    }
    Ok(())
}

fn typst_project_path(path: &str) -> String {
    typst_string(&format!("/{path}"))
}

fn entrypoint_source(
    transport: &RenderConfigTransport,
    template: Option<&str>,
    module_files: &[Option<String>],
    topology_path: Option<&str>,
) -> PyResult<String> {
    let mut source = template.map_or_else(String::new, |path| {
        format!("#import {} as _linnet_template\n", typst_project_path(path))
    });
    for (module, file) in transport.imports.iter().zip(module_files) {
        let module_source = match (&module.source, file) {
            (TypstModuleSource::File(_), Some(path)) => typst_project_path(path),
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
    source.push_str(")\n  if type(value) != dictionary {\n    panic(\"Linnet render config must be a dictionary\")\n  }\n  value");
    if let Some(path) = topology_path {
        source.push_str(&format!(
            " + (graph-spec-path: {},)",
            typst_project_path(path)
        ));
    }
    source.push_str("\n}\n");
    Ok(source)
}

/// One Typst render whose virtual project and generated entrypoint share a lifetime.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(module = "linnet", frozen)]
pub struct PreparedRender {
    _build_dir: tempfile::TempDir,
    files: BTreeMap<String, Vec<u8>>,
    root: PathBuf,
    package_store: PathBuf,
}

impl PreparedRender {
    /// Prepare a virtual Typst project with bundled Linnest, Kurvst, and packages.
    /// Sources must include `main.typ` before compilation; binary assets are staged on disk.
    pub fn from_sources(sources: BTreeMap<String, Vec<u8>>) -> PyResult<Self> {
        let build_dir =
            tempfile::tempdir().map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
        let build_root = canonicalize(build_dir.path(), "render directory")?;
        let package_store = build_root.join(TYPST_PACKAGES_DIR);
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
        // External stores add packages; bundled source and Wasm stay authoritative.
        write_embedded_assets::<EmbeddedTypstPackages>(&package_store)?;
        EmbeddedTypstPackages::stage_cetz(&package_store).map_err(|error| {
            PyRuntimeError::new_err(format!("failed to stage stock CeTZ: {error}"))
        })?;

        let mut files = BTreeMap::new();
        insert_embedded_assets::<EmbeddedLinnestPackage>(
            &mut files,
            &build_root,
            LINNEST_PACKAGE_DIR,
        )?;
        insert_embedded_assets::<EmbeddedKurvstPackage>(
            &mut files,
            &build_root,
            KURVST_PACKAGE_DIR,
        )?;
        for (path, contents) in sources {
            insert_project_asset(&mut files, &build_root, &path, &contents)?;
        }
        Ok(Self {
            _build_dir: build_dir,
            files,
            root: build_root,
            package_store,
        })
    }

    fn configuration_source(
        &mut self,
        transport: &RenderConfigTransport,
        default_template: Option<&str>,
        topology_path: Option<&str>,
    ) -> PyResult<String> {
        let files = &mut self.files;
        let build_root = &self.root;
        let mut source_paths = transport.template.iter().cloned().collect::<Vec<_>>();
        source_paths.extend(transport.imports.iter().filter_map(|import| {
            if let TypstModuleSource::File(path) = &import.source {
                Some(path.clone())
            } else {
                None
            }
        }));
        let mut staged_sources = collect_user_sources(
            files,
            build_root,
            &source_paths,
            transport.source_root.as_deref(),
        )?
        .into_iter();
        let template = if transport.template.is_some() {
            Some(
                staged_sources
                    .next()
                    .ok_or_else(|| PyRuntimeError::new_err("failed to collect Typst template"))?,
            )
        } else {
            default_template.map(str::to_owned)
        };
        let module_files = transport
            .imports
            .iter()
            .map(|import| match import.source {
                TypstModuleSource::File(_) => staged_sources
                    .next()
                    .map(Some)
                    .ok_or_else(|| PyRuntimeError::new_err("failed to collect Typst module")),
                TypstModuleSource::Package(_) => Ok(None),
            })
            .collect::<PyResult<Vec<_>>>()?;
        entrypoint_source(transport, template.as_deref(), &module_files, topology_path)
    }

    fn typst_source_value(&self) -> PyResult<String> {
        String::from_utf8(
            self.files
                .get(ENTRYPOINT)
                .expect("prepared render has an entrypoint")
                .clone(),
        )
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
        if format == "svg" {
            fs::write(&output, self.svg(py)?)
                .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
        } else {
            let pages = self.compile(format)?;
            let [page] = pages.as_slice() else {
                return Err(PyValueError::new_err(
                    "file output requires exactly one page",
                ));
            };
            fs::write(&output, page).map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
        }
        Ok(output)
    }

    /// Compile the prepared project to a single SVG page.
    pub fn svg(&self, py: Python<'_>) -> PyResult<String> {
        let pages = self.svg_pages(py)?;
        let [svg] = pages.as_slice() else {
            return Err(PyRuntimeError::new_err(format!(
                "Typst SVG render produced {} pages; expected exactly one",
                pages.len()
            )));
        };
        Self::interactive_svg(svg)
    }

    /// Compile the prepared project to one SVG document per page.
    pub fn svg_pages(&self, _py: Python<'_>) -> PyResult<Vec<String>> {
        self.compile("svg")?
            .into_iter()
            .map(|bytes| {
                String::from_utf8(bytes).map_err(|error| PyRuntimeError::new_err(error.to_string()))
            })
            .collect()
    }

    /// Compile with the embedded Rust compiler and offline package store.
    pub fn compile(&self, format: &str) -> PyResult<Vec<Vec<u8>>> {
        typst_renderer::Document::new(&self.root, &self.package_store, &self.files)
            .compile(format)
            .map_err(PyRuntimeError::new_err)
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PreparedRender {
    /// Prepare an authored main.typ document with the shared renderer assets.
    ///
    /// The document can read ``_linnet_config`` for the typed layout, drawing,
    /// style and template options. Referenced Typst modules are snapshotted.
    /// A template or selectors require Graph.prepare_render instead.
    #[staticmethod]
    #[pyo3(name = "from_sources", signature = (sources, *, config=None))]
    pub fn from_source_files(
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "builtins.dict[builtins.str, builtins.bytes]", imports=("builtins")))]
        sources: BTreeMap<String, Vec<u8>>,
        #[gen_stub(override_type(type_repr = "RenderConfig | None"))] config: Option<
            &Bound<'_, PyAny>,
        >,
    ) -> PyResult<Self> {
        let base = default_render_config(py)?;
        let config = config
            .map(|config| Bound::new(py, crate::RenderConfig::from_authored_config(config)?))
            .transpose()?;
        let transport = render_config_transport(
            py,
            &base,
            config.as_ref().map(Bound::as_any),
            &PyDict::new(py),
        )?;
        if transport.template.is_some()
            || transport.selectors.node.is_some()
            || transport.selectors.edge.is_some()
            || transport.selectors.source.is_some()
            || transport.selectors.sink.is_some()
        {
            return Err(PyValueError::new_err(
                "authored rendering accepts layout, drawing, style and template options; use Graph.prepare_render for templates or selectors",
            ));
        }
        let mut prepared = Self::from_sources(sources)?;
        let main = prepared
            .files
            .remove(ENTRYPOINT)
            .ok_or_else(|| PyValueError::new_err("render sources must include main.typ"))?;
        let mut source = prepared
            .configuration_source(&transport, None, None)?
            .into_bytes();
        source.extend(main);
        prepared.files.insert(ENTRYPOINT.to_owned(), source);
        Ok(prepared)
    }

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

pub(crate) fn render_graph(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    output: PathBuf,
    config: Option<&Bound<'_, PyAny>>,
) -> PyResult<PathBuf> {
    prepare_graph(py, graph, config, None)?.render_to(py, output)
}

pub(crate) fn graph_to_svg(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    config: Option<&Bound<'_, PyAny>>,
) -> PyResult<String> {
    prepare_graph(py, graph, config, None)?.svg(py)
}

pub(crate) fn prepare_graph(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    config: Option<&Bound<'_, PyAny>>,
    subgraph: Option<&PySubgraph>,
) -> PyResult<PreparedRender> {
    // Validate before callbacks and again after preparation of the graph records.
    if let Some(subgraph) = subgraph {
        subgraph.selection_for(py, graph, graph.borrow(py).revision()?)?;
    }
    let (topology, transport) = request(py, graph, config)?;
    let selection = if let Some(subgraph) = subgraph {
        let owner = graph.borrow(py);
        let (selected, isolated) = subgraph.selection_for(py, graph, owner.revision()?)?;
        let state = owner.state.borrow();
        let state = state.as_ref().expect("checked selection owner");
        let hedges = (0..selected.size())
            .map(|index| selected.includes(&Hedge(index)))
            .collect();
        let nodes = state
            .graph
            .iter_nodes()
            .filter_map(|(node, mut crown, _)| {
                (isolated.contains(&node.0) || crown.any(|hedge| selected.includes(&hedge)))
                    .then_some(node.0)
            })
            .collect();
        Some((hedges, nodes))
    } else {
        None
    };
    prepare(topology, transport, selection)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn package_path_replaces_read_only_cached_copies_without_changing_sources() {
        let directory = tempfile::tempdir().unwrap();
        let store = directory.path().join("staged");
        let relative = "preview/shared/1.0.0/lib.typ";
        for (name, contents) in [("cache", b"cache".as_slice()), ("path", b"path")] {
            let source = directory.path().join(name);
            let asset = source.join(relative);
            fs::create_dir_all(asset.parent().unwrap()).unwrap();
            fs::write(&asset, contents).unwrap();
            let mut permissions = fs::metadata(&asset).unwrap().permissions();
            permissions.set_readonly(true);
            fs::set_permissions(&asset, permissions).unwrap();
            copy_directory(&source, &store, "test package store").unwrap();
            assert_eq!(fs::read(&asset).unwrap(), contents);
            assert!(fs::metadata(asset).unwrap().permissions().readonly());
        }
        assert_eq!(fs::read(store.join(relative)).unwrap(), b"path");
    }

    #[test]
    fn package_overlays_replace_conflicting_entry_types_without_changing_sources() {
        for cache_is_directory in [false, true] {
            let directory = tempfile::tempdir().unwrap();
            let relative = "preview/cetz/0.5.2/src/lib.typ";
            let store = directory.path().join("staged");
            for (name, is_directory) in
                [("cache", cache_is_directory), ("path", !cache_is_directory)]
            {
                let source = directory.path().join(name);
                let asset = source.join(relative);
                let payload = if is_directory {
                    fs::create_dir_all(&asset).unwrap();
                    asset.join("contents")
                } else {
                    fs::create_dir_all(asset.parent().unwrap()).unwrap();
                    asset.clone()
                };
                fs::write(&payload, name).unwrap();
                copy_directory(&source, &store, "test package store").unwrap();
                assert_eq!(fs::read(payload).unwrap(), name.as_bytes());
                assert_eq!(asset.is_dir(), is_directory);
            }
            let staged = store.join(relative);
            assert_eq!(staged.is_dir(), !cache_is_directory);
            assert_eq!(
                fs::read(if staged.is_dir() {
                    staged.join("contents")
                } else {
                    staged
                })
                .unwrap(),
                b"path"
            );
            EmbeddedTypstPackages::stage_cetz(&store).unwrap();
            assert!(store.join(relative).is_file());
        }
    }

    #[test]
    fn bundled_packages_replace_invalid_roots_and_stale_external_entries() {
        let directory = tempfile::tempdir().unwrap();
        let stock = directory.path().join(EmbeddedTypstPackages::CETZ_PACKAGE);
        fs::create_dir_all(stock.parent().unwrap()).unwrap();
        fs::write(&stock, "external package root is a file").unwrap();
        EmbeddedTypstPackages::stage_cetz(directory.path()).unwrap();
        assert!(stock.join("typst.toml").is_file());

        let oxifmt = directory.path().join("preview/oxifmt/1.0.0");
        fs::create_dir_all(oxifmt.join("typst.toml")).unwrap();
        fs::write(oxifmt.join("external-only"), "stale").unwrap();
        write_embedded_assets::<EmbeddedTypstPackages>(directory.path()).unwrap();
        assert!(oxifmt.join("typst.toml").is_file());
        assert!(!oxifmt.join("external-only").exists());
    }

    #[test]
    fn stock_cetz_archive_replaces_read_only_copies_with_exact_release_bytes() {
        let directory = tempfile::tempdir().unwrap();
        let store = directory.path();
        let package = store.join(EmbeddedTypstPackages::CETZ_PACKAGE);
        for asset in ["src/lib.typ", "cetz-core/cetz_core.wasm"] {
            let target = package.join(asset);
            fs::create_dir_all(target.parent().unwrap()).unwrap();
            fs::write(&target, b"external override").unwrap();
            let mut permissions = fs::metadata(&target).unwrap().permissions();
            permissions.set_readonly(true);
            fs::set_permissions(target, permissions).unwrap();
        }
        fs::write(package.join("external-only.typ"), b"stale fork asset").unwrap();
        let unrelated = store.join("preview/unrelated/1.0.0/lib.typ");
        fs::create_dir_all(unrelated.parent().unwrap()).unwrap();
        fs::write(&unrelated, b"external package").unwrap();
        EmbeddedTypstPackages::stage_cetz(store).unwrap();
        EmbeddedTypstPackages::stage_cetz(store).unwrap();
        assert!(!package.join("external-only.typ").exists());
        assert_eq!(fs::read(unrelated).unwrap(), b"external package");

        let decoder = flate2::read::GzDecoder::new(EmbeddedTypstPackages::CETZ_ARCHIVE);
        let mut archive = tar::Archive::new(decoder);
        let mut file_count = 0;
        for entry in archive.entries().unwrap() {
            let mut entry = entry.unwrap();
            if !entry.header().entry_type().is_file() {
                continue;
            }
            let path = entry.path().unwrap().into_owned();
            let mut expected = Vec::new();
            std::io::Read::read_to_end(&mut entry, &mut expected).unwrap();
            assert_eq!(fs::read(package.join(&path)).unwrap(), expected, "{path:?}");
            if path == Path::new("LICENSE") {
                assert_eq!(
                    expected,
                    include_bytes!("../vendor/typst-packages/archives/LICENSE.cetz")
                );
            }
            file_count += 1;
        }
        assert!(file_count > 40);
        assert_eq!(
            WalkDir::new(&package)
                .into_iter()
                .map(Result::unwrap)
                .filter(|entry| entry.file_type().is_file())
                .count(),
            file_count
        );
    }

    #[test]
    fn embedded_package_tree_excludes_distribution_notices_and_archives() {
        let directory = tempfile::tempdir().unwrap();
        write_embedded_assets::<EmbeddedTypstPackages>(directory.path()).unwrap();
        for package in ["preview/oxifmt/1.0.0", "preview/mitex/0.2.6"] {
            assert!(directory.path().join(package).join("typst.toml").is_file());
        }
        assert!(!directory.path().join("PROVENANCE.typ").exists());
        assert!(!directory.path().join("archives").exists());
        assert!(!directory.path().join("preview/cetz/0.5.1").exists());
    }

    #[test]
    fn compiled_graph_inspection_preserves_native_drawing() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let linnet = PyModule::new(py, "linnet")?;
            crate::linnet_py(&linnet)?;
            py.import("sys")?
                .getattr("modules")?
                .set_item("linnet", &linnet)?;
            let source =
                std::ffi::CString::new(include_str!("../tests/test_svg_interaction.py")).unwrap();
            let fixture = PyModule::from_code(
                py,
                &source,
                c"test_svg_interaction.py",
                c"test_svg_interaction",
            )?;
            let case = fixture.getattr("SvgInteractionTests")?.call0()?;
            case.call_method0("setUp")?;
            let prepared = case.getattr("graph")?.call_method0("prepare_render")?;
            let prepared = prepared.extract::<PyRef<'_, PreparedRender>>()?;
            let pages = prepared.svg_pages(py)?;
            assert_eq!(pages.len(), 1);
            case.call_method1(
                "assert_native_drawing_unchanged",
                (&pages[0], prepared.svg(py)?),
            )?;
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn compiles_bundled_mitex_without_python_or_downloads() {
        let sources = BTreeMap::from([(
            "main.typ".to_owned(),
            br##"#set page(width: auto, height: auto)
#import "@preview/mitex:0.2.6": mi
#mi("\\frac{x^2}{1+y}")"##
                .to_vec(),
        )]);
        let prepared = PreparedRender::from_sources(sources).unwrap();
        let pages = prepared.compile("svg").unwrap();
        assert_eq!(pages.len(), 1);
        assert!(std::str::from_utf8(&pages[0]).unwrap().contains("<svg"));
        assert!(prepared.compile("pdf").unwrap()[0].starts_with(b"%PDF"));
        assert!(prepared.compile("png").unwrap()[0].starts_with(b"\x89PNG"));
    }

    #[test]
    fn extracts_typst_packages_for_offline_rendering() {
        let prepared = PreparedRender::from_sources(BTreeMap::new()).unwrap();
        assert!(
            prepared
                .package_store
                .join("preview/cetz/0.5.2/typst.toml")
                .is_file()
        );
    }

    #[test]
    fn embeds_linnest_and_kurvst_with_their_wasm_modules() {
        let prepared = PreparedRender::from_sources(BTreeMap::new()).unwrap();
        assert!(
            prepared
                .files
                .contains_key("crates/linnest/typst/src/graph.typ")
        );
        assert!(
            fs::metadata(prepared.root.join("crates/linnest/typst/linnest.wasm"))
                .unwrap()
                .len()
                > 0
        );
        assert!(
            fs::metadata(prepared.root.join("crates/linnest/typst/ec-layout.wasm"))
                .unwrap()
                .len()
                > 0
        );
        assert!(
            prepared
                .files
                .contains_key("crates/kurvst/typst/src/lib.typ")
        );
        assert!(
            fs::metadata(prepared.root.join("crates/kurvst/typst/kurvst.wasm"))
                .unwrap()
                .len()
                > 0
        );
    }
}
