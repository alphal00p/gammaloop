//! Immutable render snapshots, shared by generic graphs and domain renderers.
use std::fs;
use std::path::PathBuf;
use std::sync::Arc;

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::{pyclass, pymethods, PyResult};

use crate::typst_export;

/// An immutable native SVG rendering with file exports and notebook display.
/// Geometry and typeset labels are captured when the result is created.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(
    module = "symbolica.community.graph",
    name = "DiagramRender",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub struct PyDiagramRender {
    pub(crate) pages: Arc<Vec<String>>,
    pub(crate) html: Option<String>,
    /// Configured child snapshots for a rendered collection, in display order.
    #[pyo3(get)]
    diagrams: Vec<PyDiagramRender>,
}
impl PyDiagramRender {
    pub fn new(svg: String, html: String) -> Self {
        Self {
            pages: Arc::new(vec![svg]),
            html: Some(Self::without_font(&html)),
            diagrams: Vec::new(),
        }
    }
    pub fn from_diagrams(diagrams: Vec<Self>, html: String) -> PyResult<Self> {
        let pages: Vec<String> = diagrams
            .iter()
            .map(Self::to_svg_pages)
            .collect::<PyResult<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect();
        Ok(Self {
            pages: Arc::new(pages),
            html: Some(Self::without_font(&html)),
            diagrams,
        })
    }
    // Only capture's exact font block is a resource, never arbitrary label text.
    pub(crate) fn without_font(html: &str) -> String {
        let mut output = String::with_capacity(html.len());
        let mut rest = html;
        while let Some(start) = rest.find("<style data-linnet-font=\"true\">") {
            let Some(end) = rest[start..].find("</style>") else {
                break;
            };
            output.push_str(&rest[..start]);
            rest = &rest[start + end + "</style>".len()..];
        }
        output.push_str(rest);
        output
    }

    pub(crate) fn pages(&self) -> &[String] {
        &self.pages
    }
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyDiagramRender {
    /// Return the one-page SVG. Use to_svg_pages for multipage documents.
    pub fn to_svg(&self) -> PyResult<String> {
        let pages = self.pages();
        match pages {
            [svg] => Ok(svg.clone()),
            _ => Err(PyValueError::new_err(format!(
                "render has {} pages; use to_svg_pages()",
                pages.len()
            ))),
        }
    }
    pub fn to_svg_pages(&self) -> PyResult<Vec<String>> {
        Ok(self.pages().to_vec())
    }
    pub fn to_html(&self) -> PyResult<String> {
        let html = self.html.as_ref().map_or_else(
            || {
                format!(
                    "<figure class=\"linnet-graph\">{}</figure>",
                    Self::without_font(&self.pages().join("\n"))
                )
            },
            String::clone,
        );
        let font = self.pages().iter().find_map(|page| {
            let (_, rest) = page.split_once("<style data-linnet-font=\"true\"><![CDATA[")?;
            rest.split_once("]]></style>").map(|(css, _)| css)
        });
        Ok(match font {
            Some(css) => format!("<style data-linnet-font=\"true\">{css}</style>{html}"),
            None => html,
        })
    }

    fn to_linnest(&self) -> PyResult<String> {
        if self
            .pages()
            .iter()
            .any(|svg| svg.contains("class=\"linnet-math-label\""))
        {
            return Err(PyValueError::new_err(
                "Typst SVG imports require outlined labels; call view.outlined() before exporting a Typst document",
            ));
        }
        Ok(self
            .pages()
            .iter()
            .map(|svg| typst_export::svg_source(svg))
            .collect::<Vec<_>>()
            .join("\n#pagebreak()\n"))
    }
    /// Export SVG, HTML, Typst, PDF, or PNG, selected by the filename suffix.
    fn save(&self, output: PathBuf) -> PyResult<PathBuf> {
        let extension = output.extension().and_then(|v| v.to_str()).unwrap_or("");
        let text = match extension {
            "svg" => Some(self.to_svg()?),
            "html" => Some(self.to_html()?),
            "typ" => Some(self.to_linnest()?),
            "pdf" | "png" => None,
            _ => {
                return Err(PyValueError::new_err(
                    "expected .svg, .html, .typ, .pdf, or .png",
                ));
            }
        };
        if let Some(parent) = output.parent().filter(|path| !path.as_os_str().is_empty()) {
            fs::create_dir_all(parent).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        }
        if let Some(text) = text {
            fs::write(&output, text).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
            Ok(output)
        } else {
            let bytes = self.export_pages(extension)?;
            match bytes.as_slice() {
                [page] => {
                    fs::write(&output, page).map_err(|e| PyRuntimeError::new_err(e.to_string()))?
                }
                _ => {
                    return Err(PyValueError::new_err(
                        "use save_pages() for multiple PNG pages",
                    ));
                }
            }
            Ok(output)
        }
    }

    /// Export every SVG or PNG page to a directory and return the written paths.
    #[pyo3(signature = (directory, *, format="svg"))]
    fn save_pages(&self, directory: PathBuf, format: &str) -> PyResult<Vec<PathBuf>> {
        let pages = match format {
            "svg" => self
                .pages()
                .iter()
                .map(|page| page.as_bytes().to_vec())
                .collect(),
            "png" => self.export_pages("png")?,
            _ => {
                return Err(PyValueError::new_err(
                    "page format must be svg or png; save a PDF for a multipage document",
                ));
            }
        };
        fs::create_dir_all(&directory).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        pages
            .into_iter()
            .enumerate()
            .map(|(index, page)| {
                let path = directory.join(format!("page-{}.{format}", index + 1));
                fs::write(&path, page).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
                Ok(path)
            })
            .collect()
    }

    fn _repr_html_(&self) -> PyResult<String> {
        self.to_html()
    }
    fn _repr_svg_(&self) -> PyResult<Option<String>> {
        Ok(match self.pages() {
            [svg] => Some(svg.clone()),
            _ => None,
        })
    }
    fn _mime_(&self) -> PyResult<(&str, String)> {
        Ok(("text/html", self.to_html()?))
    }
    fn __repr__(&self) -> &'static str {
        "DiagramRender()"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn html_shares_one_font_while_svg_pages_remain_standalone() {
        let font =
            "<style data-linnet-font=\"true\"><![CDATA[@font-face{font-family:test}]]></style>";
        let svg = format!("<svg>{font}<path/></svg>");
        let child = PyDiagramRender::new(svg.clone(), format!("<figure>{svg}</figure>"));
        assert!(!child.html.as_ref().unwrap().contains("@font-face"));
        assert_eq!(child.to_html().unwrap().matches("@font-face").count(), 1);
        let collection = PyDiagramRender::from_diagrams(
            vec![child.clone(), child.clone(), child],
            format!("{svg}{svg}{svg}"),
        )
        .unwrap();
        assert_eq!(
            collection.to_html().unwrap().matches("@font-face").count(),
            1
        );
        assert!(collection
            .pages()
            .iter()
            .all(|page| page.matches("@font-face").count() == 1));
    }
}
