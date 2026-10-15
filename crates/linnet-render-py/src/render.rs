use std::collections::BTreeMap;
use std::env;
use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};

use crate::typst::{
    default_render_config, render_config_transport, typst_string, RenderSettingsTransport,
    TypstModuleSource,
};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use rust_embed::RustEmbed;
use walkdir::WalkDir;
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
struct EmbeddedTypstPackages;

pub fn prepare(
    topology: Vec<u8>,
    transport: RenderSettingsTransport,
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
        if target.exists() {
            fs::remove_file(&target).map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
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
    transport: &RenderSettingsTransport,
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
        fs::create_dir_all(build_root.join(TYPST_PACKAGES_DIR))
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
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

        // Bundled packages include required patches; user caches may add packages,
        // but cannot replace these compiler inputs with an incompatible release.
        write_embedded_assets::<EmbeddedTypstPackages>(&package_store)?;
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
        transport: &RenderSettingsTransport,
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

    pub fn typst_source_value(&self) -> PyResult<String> {
        String::from_utf8(
            self.files
                .get(ENTRYPOINT)
                .expect("prepared render has an entrypoint")
                .clone(),
        )
        .map_err(|error| PyRuntimeError::new_err(error.to_string()))
    }

    pub fn render_to(&self, py: Python<'_>, output: PathBuf) -> PyResult<PathBuf> {
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
            .map_err(PyRuntimeError::new_err)?
            .compile(format)
            .map_err(PyRuntimeError::new_err)
    }
}

impl PreparedRender {
    /// Prepare an authored main.typ document with the shared renderer assets.
    ///
    /// The document can read ``_linnet_config`` for the typed layout, drawing,
    /// style and template options. Referenced Typst modules are snapshotted.
    /// A template or selectors require Graph.render instead.
    pub fn from_source_files(
        py: Python<'_>,

        sources: BTreeMap<String, Vec<u8>>,
        config: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let base = default_render_config(py)?;
        let transport = render_config_transport(py, &base, config, &PyDict::new(py))?;
        if transport.template.is_some()
            || transport.selectors.node.is_some()
            || transport.selectors.edge.is_some()
            || transport.selectors.source.is_some()
            || transport.selectors.sink.is_some()
        {
            return Err(PyValueError::new_err(
                "authored rendering accepts layout, drawing, style and template options; use Graph.render for templates or selectors",
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

#[cfg(test)]
mod tests {
    use super::*;
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
        assert!(prepared
            .package_store
            .join("preview/cetz/0.5.1/typst.toml")
            .is_file());
    }

    #[test]
    fn embeds_linnest_and_kurvst_with_their_wasm_modules() {
        let prepared = PreparedRender::from_sources(BTreeMap::new()).unwrap();
        assert!(prepared
            .files
            .contains_key("crates/linnest/typst/src/graph.typ"));
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
        assert!(prepared
            .files
            .contains_key("crates/kurvst/typst/src/lib.typ"));
        assert!(
            fs::metadata(prepared.root.join("crates/kurvst/typst/kurvst.wasm"))
                .unwrap()
                .len()
                > 0
        );
    }
}
