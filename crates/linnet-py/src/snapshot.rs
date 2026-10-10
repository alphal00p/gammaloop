//! Immutable render snapshots, shared by generic graphs and domain renderers.
//!
//! A snapshot draws two outputs on demand. Static pages hold the drawing alone,
//! for files, documents and exports. Interactive pages add the decorated
//! inspection targets of notebook and HTML display: hover details, selection
//! and keyboard focus. Static output needs no link per target, so its files stay
//! small and its PDFs carry no link annotations.
use crate::render::PreparedRender;
use linnest::svg::{Figure, Scene};
use pyo3::{
    exceptions::{PyRuntimeError, PyValueError},
    prelude::*,
};
use std::{
    collections::BTreeMap,
    path::PathBuf,
    sync::{Arc, OnceLock},
};

/// The source of a snapshot's pages.
#[derive(Clone)]
enum Drawing {
    /// A prepared Typst project, compiled once with inspection targets for both outputs.
    Project(Arc<PreparedRender>),
    /// A natively drawn figure, serialized once per output.
    Native(Arc<Figure>),
    /// The pages of the child snapshots in `diagrams`, in display order.
    Pages,
    /// Independently drawn channels side by side on one page.
    Row(Vec<PyDiagramRender>),
}

/// An immutable graph rendering with source inspection, file exports, and notebook display.
/// Render inputs and referenced assets are captured when the result is created.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "DiagramRender", frozen, skip_from_py_object)]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
#[derive(Clone)]
pub struct PyDiagramRender {
    drawing: Drawing,
    /// Original configuration source when geometry was rendered natively.
    pub(crate) source: Option<String>,
    /// Rewrites applied in order to every page of both outputs.
    themes: Vec<fn(String) -> String>,
    /// Static pages, drawn on first use.
    pages: Arc<OnceLock<Vec<String>>>,
    /// Interactive pages, drawn on first use.
    interactive_pages: Arc<OnceLock<Vec<String>>>,
    html: Option<String>,
    /// Configured child snapshots for a rendered collection, in display order.
    #[pyo3(get)]
    diagrams: Vec<PyDiagramRender>,
}
impl PyDiagramRender {
    fn new(drawing: Drawing) -> Self {
        Self {
            drawing,
            source: None,
            themes: Vec::new(),
            pages: Arc::default(),
            interactive_pages: Arc::default(),
            html: None,
            diagrams: Vec::new(),
        }
    }
    pub(crate) fn prepared(project: PreparedRender) -> Self {
        Self::new(Drawing::Project(Arc::new(project)))
    }
    pub fn from_figure(figure: Figure) -> Self {
        Self::new(Drawing::Native(Arc::new(figure)))
    }
    pub fn from_diagrams(diagrams: Vec<Self>, html: String) -> Self {
        Self {
            diagrams,
            ..Self::new(Drawing::Pages)
        }
        .with_html(html)
    }
    /// Arrange single-page snapshots side by side on one page.
    pub fn row(channels: Vec<Self>) -> Self {
        Self::new(Drawing::Row(channels))
    }
    pub fn with_html(mut self, html: String) -> Self {
        self.html = Some(html);
        self
    }
    /// Rewrite every page of both outputs, for example to apply a theme.
    pub fn map_svg(mut self, theme: fn(String) -> String) -> Self {
        self.themes.push(theme);
        self.pages = Arc::default();
        self.interactive_pages = Arc::default();
        self
    }
    /// Compile the snapshot's own static pages: its project, or its SVG pages,
    /// which need no project, package or plugin.
    fn compile(&self, format: &str) -> PyResult<Vec<Vec<u8>>> {
        match &self.drawing {
            Drawing::Project(project) => project.compile(format),
            _ => typst_renderer::Document::compile_sources(
                &BTreeMap::from([("main.typ".into(), self.to_linnest()?.into_bytes())]),
                format,
            )
            .map_err(PyRuntimeError::new_err),
        }
    }
    /// The static or interactive pages, drawn and themed on first use.
    pub(crate) fn pages(&self, interactive: bool) -> PyResult<&[String]> {
        let cell = if interactive {
            &self.interactive_pages
        } else {
            &self.pages
        };
        if let Some(pages) = cell.get() {
            return Ok(pages);
        }
        let pages = match &self.drawing {
            Drawing::Project(project) => project.svg_pages(interactive)?,
            Drawing::Native(figure) if interactive => {
                vec![figure.interactive_svg().map_err(PyRuntimeError::new_err)?]
            }
            Drawing::Native(figure) => vec![figure.svg()],
            Drawing::Pages => self
                .diagrams
                .iter()
                .map(|diagram| diagram.pages(interactive).map(<[String]>::to_vec))
                .collect::<PyResult<Vec<_>>>()?
                .concat(),
            Drawing::Row(channels) => vec![Scene::combine_svgs(
                &channels
                    .iter()
                    .map(|channel| channel.to_svg(interactive))
                    .collect::<PyResult<Vec<_>>>()?,
            )
            .map_err(PyRuntimeError::new_err)?],
        };
        let pages = pages
            .into_iter()
            .map(|page| self.themes.iter().fold(page, |page, theme| theme(page)))
            .collect();
        Ok(cell.get_or_init(|| pages))
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
    ///
    /// The SVG is static by default: the drawing alone, for files and documents.
    /// ``interactive=True`` adds the hover, selection and keyboard targets of
    /// notebook display, with their script and style.
    #[pyo3(signature = (*, interactive=false))]
    pub fn to_svg(&self, interactive: bool) -> PyResult<String> {
        match self.pages(interactive)? {
            [svg] => Ok(svg.clone()),
            pages => Err(PyValueError::new_err(format!(
                "render has {} pages; use to_svg_pages()",
                pages.len()
            ))),
        }
    }
    /// Return one SVG per page, static or ``interactive`` as in ``to_svg``.
    #[pyo3(signature = (*, interactive=false))]
    pub fn to_svg_pages(&self, interactive: bool) -> PyResult<Vec<String>> {
        Ok(self.pages(interactive)?.to_vec())
    }
    /// Return the interactive HTML display.
    pub fn to_html(&self) -> PyResult<String> {
        if let Some(html) = &self.html {
            return Ok(html.clone());
        }
        Ok(format!(
            "<figure class=\"linnet-graph\">{}</figure>",
            self.pages(true)?.join("\n")
        ))
    }
    #[getter]
    fn typst_source(&self) -> PyResult<String> {
        if let Some(source) = &self.source {
            return Ok(source.clone());
        }
        match &self.drawing {
            Drawing::Project(project) => project.typst_source_value(),
            _ => self.to_linnest(),
        }
    }
    fn to_linnest(&self) -> PyResult<String> {
        Ok(self
            .pages(false)?
            .iter()
            .map(|svg| typst_renderer::Document::svg_source(svg))
            .collect::<Vec<_>>()
            .join("\n#pagebreak()\n"))
    }
    /// Export SVG, HTML, Typst, PDF, or PNG, selected by the filename suffix.
    /// HTML is the interactive display; the other formats are static.
    fn save(&self, output: PathBuf) -> PyResult<PathBuf> {
        let extension = output.extension().and_then(|v| v.to_str()).unwrap_or("");
        let text = match extension {
            "svg" => Some(self.to_svg(false)?),
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
        let contents = match text {
            Some(text) => text.into_bytes(),
            None => match self.compile(extension)?.as_slice() {
                [page] => page.clone(),
                _ => {
                    return Err(PyValueError::new_err(
                        "file output requires exactly one page",
                    ));
                }
            },
        };
        std::fs::write(&output, contents).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(output)
    }
    /// Export every static SVG or PNG page to a directory and return the written paths.
    #[pyo3(signature = (directory, *, format="svg"))]
    fn save_pages(&self, directory: PathBuf, format: &str) -> PyResult<Vec<PathBuf>> {
        let pages = match format {
            "svg" => self
                .pages(false)?
                .iter()
                .map(|page| page.as_bytes().to_vec())
                .collect(),
            "png" => self.compile("png")?,
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
    /// The static SVG, for frontends that show or export images.
    fn _repr_svg_(&self) -> PyResult<Option<String>> {
        Ok(match self.pages(false)? {
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
