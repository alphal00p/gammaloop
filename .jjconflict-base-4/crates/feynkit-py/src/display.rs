use std::{fs, path::Path, sync::OnceLock};

use feynkit_graph::FeynmanDiagram;
use include_dir::{Dir, include_dir};
use pyo3::{
    exceptions::{PyImportError, PyRuntimeError, PyValueError},
    prelude::*,
    types::{PyBytes, PyDict, PyModule},
};

const LINNEST_PACKAGE_DIR: &str = "crates/linnest/typst";
const KURVST_PACKAGE_DIR: &str = "crates/kurvst/typst";

static LINNEST_SOURCE: Dir<'_> = include_dir!("$CARGO_MANIFEST_DIR/../linnest/typst/src");
static KURVST_SOURCE: Dir<'_> = include_dir!("$CARGO_MANIFEST_DIR/../kurvst/typst/src");
static LINNEST_WASM: &[u8] = include_bytes!("../../linnest/typst/linnest.wasm");
static KURVST_WASM: &[u8] = include_bytes!("../../kurvst/typst/kurvst.wasm");
static TYPST_ASSETS: OnceLock<Result<tempfile::TempDir, String>> = OnceLock::new();

pub(crate) fn escape_html(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            '\'' => escaped.push_str("&#39;"),
            _ => escaped.push(character),
        }
    }
    escaped
}

fn extract_typst_assets() -> Result<tempfile::TempDir, String> {
    let root = tempfile::Builder::new()
        .prefix("feynkit-typst-")
        .tempdir()
        .map_err(|error| format!("could not create the FeynKit Typst asset directory: {error}"))?;
    let linnest = root.path().join(LINNEST_PACKAGE_DIR);
    let kurvst = root.path().join(KURVST_PACKAGE_DIR);
    fs::create_dir_all(linnest.join("src"))
        .map_err(|error| format!("could not create the Linnest package directory: {error}"))?;
    fs::create_dir_all(kurvst.join("src"))
        .map_err(|error| format!("could not create the Kurvst package directory: {error}"))?;
    LINNEST_SOURCE
        .extract(linnest.join("src"))
        .map_err(|error| format!("could not extract the embedded Linnest package: {error}"))?;
    KURVST_SOURCE
        .extract(kurvst.join("src"))
        .map_err(|error| format!("could not extract the embedded Kurvst package: {error}"))?;
    fs::write(linnest.join("linnest.wasm"), LINNEST_WASM)
        .map_err(|error| format!("could not extract Linnest's Wasm module: {error}"))?;
    fs::write(kurvst.join("kurvst.wasm"), KURVST_WASM)
        .map_err(|error| format!("could not extract Kurvst's Wasm module: {error}"))?;
    Ok(root)
}

fn typst_asset_root() -> PyResult<&'static Path> {
    match TYPST_ASSETS.get_or_init(extract_typst_assets) {
        Ok(root) => Ok(root.path()),
        Err(message) => Err(PyRuntimeError::new_err(message.clone())),
    }
}

/// Compile a diagram's Linnest source with typst-py and return its SVG page.
pub(crate) fn render_diagram_svg(py: Python<'_>, diagram: &FeynmanDiagram) -> PyResult<String> {
    let typst = PyModule::import(py, "typst").map_err(|error| {
        if error.is_instance_of::<PyImportError>(py) {
            PyImportError::new_err(
                "FeynKit diagram rendering requires typst-py; install the Python package `typst>=0.15,<0.16`",
            )
        } else {
            error
        }
    })?;
    let source = diagram.to_linnest();
    let kwargs = PyDict::new(py);
    kwargs.set_item("root", typst_asset_root()?.to_string_lossy().as_ref())?;
    kwargs.set_item("format", "svg")?;
    kwargs.set_item("ignore_system_fonts", true)?;
    kwargs.set_item("pretty", false)?;

    let output = typst
        .getattr("compile")?
        .call((PyBytes::new(py, source.as_bytes()),), Some(&kwargs))?;
    let bytes = output.cast::<PyBytes>().map_err(|_| {
        PyRuntimeError::new_err(
            "typst-py returned multiple pages while rendering a single FeynKit diagram",
        )
    })?;
    String::from_utf8(bytes.as_bytes().to_vec())
        .map_err(|error| PyValueError::new_err(format!("Typst returned non-UTF-8 SVG: {error}")))
}

pub(crate) fn render_diagram_html(py: Python<'_>, diagram: &FeynmanDiagram) -> PyResult<String> {
    let svg = render_diagram_svg(py, diagram)?;
    Ok(format!(
        "<figure class=\"feynkit-diagram\" style=\"max-width:100%;margin:.5rem 0\">\
         <div style=\"max-width:100%;overflow-x:auto\">{svg}</div>\
         <figcaption style=\"font-size:.85em;opacity:.75;text-align:center\">Feynman diagram {} \
         ({} loop{})</figcaption></figure>",
        escape_html(diagram.name()),
        diagram.loop_count(),
        if diagram.loop_count() == 1 { "" } else { "s" },
    ))
}

#[cfg(test)]
mod tests {
    use super::{KURVST_SOURCE, KURVST_WASM, LINNEST_SOURCE, LINNEST_WASM, escape_html};

    #[test]
    fn escapes_html_metadata() {
        assert_eq!(
            escape_html("<script data-x='a&b'>\"x\"</script>"),
            "&lt;script data-x=&#39;a&amp;b&#39;&gt;&quot;x&quot;&lt;/script&gt;"
        );
    }

    #[test]
    fn embeds_linnest_and_kurvst_with_their_wasm_modules() {
        assert!(LINNEST_SOURCE.get_file("graph.typ").is_some());
        assert!(!LINNEST_WASM.is_empty());
        assert!(KURVST_SOURCE.get_file("lib.typ").is_some());
        assert!(!KURVST_WASM.is_empty());
    }
}
