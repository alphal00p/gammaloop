//! Immutable render snapshots, shared by generic graphs and domain renderers.
use crate::render::PreparedRender;
use pyo3::{
    exceptions::{PyRuntimeError, PyValueError},
    prelude::*,
};
use std::{
    collections::BTreeMap,
    path::PathBuf,
    sync::{Arc, OnceLock},
};

/// An immutable graph rendering with source inspection, file exports, and notebook display.
/// Render inputs and referenced assets are captured when the result is created.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(
    module = "symbolica.community.graph",
    name = "DiagramRender",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub struct PyDiagramRender {
    project: Option<Arc<PreparedRender>>,
    pub(crate) pages: Arc<OnceLock<Vec<String>>>,
    pub(crate) html: Option<String>,
    /// Configured child snapshots for a rendered collection, in display order.
    #[pyo3(get)]
    diagrams: Vec<PyDiagramRender>,
}
impl PyDiagramRender {
    pub fn new(svg: String, html: String) -> Self {
        Self {
            project: None,
            pages: Arc::new(OnceLock::from(vec![svg])),
            html: Some(html),
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
            project: None,
            pages: Arc::new(OnceLock::from(pages)),
            html: Some(html),
            diagrams,
        })
    }
    pub(crate) fn prepared(project: PreparedRender) -> Self {
        Self {
            project: Some(Arc::new(project)),
            pages: Arc::new(OnceLock::new()),
            html: None,
            diagrams: Vec::new(),
        }
    }
    pub(crate) fn pages(&self) -> PyResult<&Vec<String>> {
        if self.pages.get().is_none() {
            let pages = Python::attach(|py| self.project.as_ref().unwrap().svg_pages(py))?;
            let pages = pages
                .iter()
                .map(|svg| {
                    linnest::svg::Scene::interactive_svg(svg).map_err(PyRuntimeError::new_err)
                })
                .collect::<PyResult<Vec<_>>>()?;
            let _ = self.pages.set(pages);
        }
        Ok(self.pages.get().unwrap())
    }
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyDiagramRender {
    /// Snapshot an authored Typst project, including its referenced assets.
    #[staticmethod]
    #[pyo3(signature = (sources, *, config=None))]
    fn from_sources(
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "dict[str, bytes]"))] sources: BTreeMap<
            String,
            Vec<u8>,
        >,
        config: Option<&crate::PyRenderSettings>,
    ) -> PyResult<Self> {
        let config = config.map(|c| Py::new(py, c.clone())).transpose()?;
        Ok(Self::prepared(PreparedRender::from_source_files(
            py,
            sources,
            config.as_ref().map(|c| c.bind(py).as_any()),
        )?))
    }
    /// Return the one-page SVG. Use to_svg_pages for multipage documents.
    pub fn to_svg(&self) -> PyResult<String> {
        let pages = self.pages()?;
        match pages.as_slice() {
            [svg] => Ok(svg.clone()),
            _ => Err(PyValueError::new_err(format!(
                "render has {} pages; use to_svg_pages()",
                pages.len()
            ))),
        }
    }
    pub fn to_svg_pages(&self) -> PyResult<Vec<String>> {
        Ok(self.pages()?.clone())
    }
    pub fn to_html(&self) -> PyResult<String> {
        if let Some(html) = &self.html {
            return Ok(html.clone());
        }
        Ok(format!(
            "<figure class=\"linnet-graph\">{}</figure>",
            self.pages()?.join("\n")
        ))
    }
    #[getter]
    fn typst_source(&self) -> PyResult<String> {
        match &self.project {
            Some(project) => project.typst_source_value(),
            None => self.to_linnest(),
        }
    }
    fn to_linnest(&self) -> PyResult<String> {
        Ok(self
            .pages()?
            .iter()
            .map(|svg| typst_renderer::Document::svg_source(svg))
            .collect::<Vec<_>>()
            .join("\n#pagebreak()\n"))
    }
    /// Export SVG, HTML, Typst, PDF, or PNG, selected by the filename suffix.
    fn save(&self, py: Python<'_>, output: PathBuf) -> PyResult<PathBuf> {
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
            std::fs::create_dir_all(parent).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        }
        if let Some(text) = text {
            std::fs::write(&output, text).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
            Ok(output)
        } else if let Some(project) = &self.project {
            project.render_to(py, output)
        } else {
            PreparedRender::from_sources(BTreeMap::from([(
                "main.typ".into(),
                self.to_linnest()?.into_bytes(),
            )]))?
            .render_to(py, output)
        }
    }
    /// Export every SVG or PNG page to a directory and return the written paths.
    #[pyo3(signature = (directory, *, format="svg"))]
    fn save_pages(&self, directory: PathBuf, format: &str) -> PyResult<Vec<PathBuf>> {
        let pages = match format {
            "svg" => self
                .pages()?
                .iter()
                .map(|page| page.as_bytes().to_vec())
                .collect(),
            "png" => match &self.project {
                Some(project) => project.compile("png")?,
                None => typst_renderer::Document::compile_sources(
                    &BTreeMap::from([("main.typ".into(), self.to_linnest()?.into_bytes())]),
                    "png",
                )
                .map_err(PyRuntimeError::new_err)?,
            },
            _ => {
                return Err(PyValueError::new_err(
                    "page format must be svg or png; save a PDF for a multipage document",
                ));
            }
        };
        std::fs::create_dir_all(&directory).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        pages
            .into_iter()
            .enumerate()
            .map(|(index, page)| {
                let path = directory.join(format!("page-{}.{format}", index + 1));
                std::fs::write(&path, page).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
                Ok(path)
            })
            .collect()
    }

    fn _repr_html_(&self) -> PyResult<String> {
        self.to_html()
    }
    fn _repr_svg_(&self) -> PyResult<Option<String>> {
        Ok(match self.pages()?.as_slice() {
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
