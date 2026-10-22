//! Optional Python Typst compilation and source-only figure exports.
use std::{
    collections::BTreeMap,
    path::{Component, Path},
};

use pyo3::{
    exceptions::{PyImportError, PyRuntimeError, PyValueError},
    prelude::*,
    types::{PyBytes, PyDict, PyList},
};

/// Compile named source files, importing the optional Python compiler on demand.
/// Files are staged byte-for-byte in a scoped project rooted at `main.typ`.
/// Byte output is one document; list output retains every SVG or PNG page.
pub fn compile_sources(
    py: Python<'_>,
    files: &BTreeMap<String, Vec<u8>>,
    format: &str,
) -> PyResult<Vec<Vec<u8>>> {
    let typst = py.import("typst").map_err(|cause| {
        if !cause.is_instance_of::<PyImportError>(py) {
            return cause;
        }
        let error = PyImportError::new_err(
            "Typst exports require the optional Python package 'typst'; install it with `python -m pip install typst`",
        );
        error.set_cause(py, Some(cause));
        error
    })?;
    compile(&typst, files, format)
}

fn compile(
    typst: &Bound<'_, PyModule>,
    files: &BTreeMap<String, Vec<u8>>,
    format: &str,
) -> PyResult<Vec<Vec<u8>>> {
    let py = typst.py();
    if !files.contains_key("main.typ") {
        return Err(PyValueError::new_err(
            "Typst export requires a 'main.typ' project file",
        ));
    }
    for name in files.keys() {
        if name.is_empty()
            || Path::new(name)
                .components()
                .any(|component| !matches!(component, Component::Normal(_)))
        {
            return Err(PyValueError::new_err(format!(
                "Typst project filename must be relative without '..' components: {name:?}"
            )));
        }
    }
    // typst-py 0.15 parses every dictionary value as UTF-8, including binary
    // assets. A file-backed project keeps CBOR and other assets byte-exact.
    let project = tempfile::tempdir()?;
    let root = project.path().canonicalize()?;
    for (name, bytes) in files {
        let path = root.join(name);
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        std::fs::write(path, bytes)?;
    }
    let main = root.join("main.typ");
    let kwargs = PyDict::new(py);
    kwargs.set_item("root", root.as_os_str())?;
    kwargs.set_item("format", format)?;
    let output = typst
        .call_method("compile", (main.as_os_str(),), Some(&kwargs))
        .map_err(|cause| {
            let error = PyRuntimeError::new_err(format!("Typst export failed: {cause}"));
            error.set_cause(py, Some(cause));
            error
        })?;
    if let Ok(bytes) = output.cast::<PyBytes>() {
        return Ok(vec![bytes.as_bytes().to_vec()]);
    }
    output
        .cast::<PyList>()?
        .iter()
        .map(|page| Ok(page.cast::<PyBytes>()?.as_bytes().to_vec()))
        .collect()
}

/// A self-contained Typst export of an already drawn SVG figure.
pub fn svg_source(svg: &str) -> String {
    format!("#set page(width: auto, height: auto, margin: 0pt)\n#image(bytes({svg:?}))")
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::exceptions::{PyRuntimeError, PyTypeError};
    use std::path::PathBuf;

    #[test]
    fn svg_source_preserves_escaped_image_content() {
        assert_eq!(
            svg_source("<svg data-name=\"α\\β\">\n</svg>"),
            "#set page(width: auto, height: auto, margin: 0pt)\n#image(bytes(\"<svg data-name=\\\"α\\\\β\\\">\\n</svg>\"))"
        );
        assert_eq!(
            svg_source(""),
            "#set page(width: auto, height: auto, margin: 0pt)\n#image(bytes(\"\"))"
        );
    }

    #[test]
    fn optional_compiler_preserves_binary_files_pages_and_cleans_up_projects() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let typst = PyModule::from_code(
                py,
                cr#"
from pathlib import Path

roots = []

def compile(main, *, root, format):
    main, root = Path(main), Path(root)
    assert main.is_absolute() and root.is_absolute()
    assert main == root / 'main.typ'
    assert {path.relative_to(root).as_posix(): path.read_bytes()
            for path in root.rglob('*') if path.is_file()} == {
        'main.typ': b'$ x $',
        'nested/tree.cbor': bytes([0, 255, 128]),
        'unused.bin': bytes([255]),
    }
    roots.append(root)
    if format == 'svg':
        return [b'<svg>first</svg>', b'<svg>second</svg>']
    if format == 'png':
        return [bytes([137, 80, 78, 71, 0]), bytes([137, 80, 78, 71, 255])]
    if format == 'pdf':
        return b'%PDF-1.7' + bytes([0, 255])
    if format == 'html':
        return b'<html><body>one document</body></html>'
    if format == 'invalid':
        return [b'page', 'not bytes']
    raise RuntimeError('compiler diagnostic')
"#,
                c"test_typst_export.py",
                c"test_typst_export",
            )?;
            let files = BTreeMap::from([
                ("main.typ".into(), b"$ x $".to_vec()),
                ("nested/tree.cbor".into(), vec![0, 255, 128]),
                ("unused.bin".into(), vec![255]),
            ]);
            let assert_cleaned = || -> PyResult<()> {
                let root = typst
                    .getattr("roots")?
                    .call_method0("pop")?
                    .extract::<PathBuf>()?;
                assert!(!root.exists(), "temporary project retained: {root:?}");
                Ok(())
            };
            for (format, expected) in [
                (
                    "svg",
                    vec![b"<svg>first</svg>".to_vec(), b"<svg>second</svg>".to_vec()],
                ),
                (
                    "png",
                    vec![vec![137, 80, 78, 71, 0], vec![137, 80, 78, 71, 255]],
                ),
                ("pdf", vec![b"%PDF-1.7\0\xff".to_vec()]),
                (
                    "html",
                    vec![b"<html><body>one document</body></html>".to_vec()],
                ),
            ] {
                assert_eq!(compile(&typst, &files, format)?, expected);
                assert_cleaned()?;
            }
            let error = compile(&typst, &files, "invalid").unwrap_err();
            assert_cleaned()?;
            assert!(error.is_instance_of::<PyTypeError>(py));
            let error = compile(&typst, &files, "failure").unwrap_err();
            assert_cleaned()?;
            assert!(error.is_instance_of::<PyRuntimeError>(py));
            assert!(error
                .value(py)
                .str()?
                .to_str()?
                .contains("compiler diagnostic"));
            assert_eq!(
                error.cause(py).unwrap().value(py).str()?.to_str()?,
                "compiler diagnostic"
            );
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn optional_compiler_requires_main_and_relative_project_paths() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let typst = PyModule::from_code(
                py,
                c"def compile(*args, **kwargs):\n    raise AssertionError('invalid project reached compiler')\n",
                c"test_typst_export_paths.py",
                c"test_typst_export_paths",
            )?;
            let error = compile(&typst, &BTreeMap::new(), "svg").unwrap_err();
            assert!(error.is_instance_of::<PyValueError>(py));
            assert!(error.to_string().contains("main.typ"));
            for name in ["", ".", "..", "../escape.bin", "nested/../escape.bin", "/escape.bin"] {
                let files = BTreeMap::from([
                    ("main.typ".into(), b"$ x $".to_vec()),
                    (name.into(), vec![255]),
                ]);
                let error = compile(&typst, &files, "svg").unwrap_err();
                assert!(error.is_instance_of::<PyValueError>(py), "{name:?}: {error}");
                assert!(error.to_string().contains("relative"));
            }
            Ok(())
        })
        .unwrap();
    }
}
