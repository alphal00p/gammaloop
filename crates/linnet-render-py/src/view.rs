//! Prepared native scenes become immutable snapshots after real label measurement.
use std::collections::BTreeMap;
use std::ffi::CString;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, OnceLock};

use linnest::svg::{LabelMetrics, LabelSource, Scene};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyAny, PyAnyMethods, PyBytes, PyDict, PyDictMethods, PyModule};
use pyo3::{pyclass, pymethods, Bound, Py, PyResult, Python};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
#[cfg(not(feature = "python_stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;
use serde::Deserialize;
use serde_json::json;

use crate::data_assets::{STIX_TWO_MATH, STIX_TWO_MATH_LICENSE};
use crate::{typst_export, PyDiagramRender};

type SvgMap = Arc<dyn Fn(String) -> String + Send + Sync>;

enum Source {
    Scenes(Vec<Scene>),
    Views(Vec<PyDiagramView>),
    Captured(PyDiagramRender),
}

struct Completed {
    metrics: Vec<LabelMetrics>,
    snapshot: PyDiagramRender,
}

struct FontResource {
    css: String,
    key: String,
}

/// A resolved drawing awaiting its notebook browser's label measurements.
///
/// Clones retain the same scenes, receipt identity, and completion cells. A
/// widget and its caller therefore observe the same completed snapshot. Each
/// freshly constructed view has a new identity; receipts belong to that view,
/// and a completed view rejects different measurements. SVG/HTML exports are
/// immutable snapshots and never trigger implicit remeasurement.
///
/// ```no_run
/// use linnest::svg::Scene;
/// use linnet_render_py::PyDiagramView;
/// use pyo3::{PyResult, Python};
/// # fn example(scene: Scene) -> PyResult<()> {
/// Python::initialize();
/// Python::attach(|py| -> PyResult<()> {
///     let view = PyDiagramView::from_scene(py, scene, None)?;
///     let retained = view.clone(); // Shares completion with the notebook owner.
///     let request = view.measurement_request(py)?;
///     // Send request to the browser; measured() accepts its complete receipt.
///     // Scripts can explicitly request Chromium instead (requires Playwright).
///     let snapshot = view.measured(py, None, Some("chromium"))?;
///     assert_eq!(retained.to_svg_pages(py)?, snapshot.to_svg_pages()?);
///     let collection = PyDiagramView::from_views(vec![view], None)?;
///     assert_eq!(collection.diagrams().len(), 1);
///     let svg = snapshot.to_svg()?;
///     // Optional vector compilation is separate from the measured snapshot.
///     let outlined = collection.outlined(py)?;
///     assert_eq!(snapshot.to_svg()?, svg);
///     # Ok(())
/// })
/// # }
/// ```
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    module = "symbolica.community.graph",
    name = "DiagramView",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub struct PyDiagramView {
    id: String,
    source: Arc<Source>,
    html_template: Option<String>,
    maps: Vec<SvgMap>,
    completed: Arc<OnceLock<Completed>>,
    vector: Arc<OnceLock<PyDiagramRender>>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Receipt {
    id: String,
    metrics: Vec<MeasuredLabel>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct MeasuredLabel {
    id: usize,
    width: f64,
    height: f64,
    baseline: f64,
}

impl PyDiagramView {
    fn new(source: Source, html_template: Option<String>) -> Self {
        static NEXT_ID: AtomicU64 = AtomicU64::new(1);
        Self {
            id: format!("linnet-view-{}", NEXT_ID.fetch_add(1, Ordering::Relaxed)),
            source: Arc::new(source),
            html_template,
            maps: Vec::new(),
            completed: Arc::new(OnceLock::new()),
            vector: Arc::new(OnceLock::new()),
        }
    }

    pub(crate) fn prepared(scenes: Vec<Scene>) -> Self {
        Self::new(Source::Scenes(scenes), None)
    }

    /// Collections retain child owners; their template slots contain child SVGs.
    ///
    /// See [`Self`] for an example. The collection gets its own receipt identity;
    /// its ordered probes recursively include all children. Completing it also
    /// completes those retained children. Construction itself is infallible;
    /// the result type is retained as part of the external rendering API.
    pub fn from_views(views: Vec<Self>, html_template: Option<String>) -> PyResult<Self> {
        Ok(Self::new(Source::Views(views), html_template))
    }

    /// Include a captured render without reconstructing or measuring its geometry.
    ///
    /// ```no_run
    /// use linnet_render_py::{PyDiagramRender, PyDiagramView};
    /// # fn captured(snapshot: PyDiagramRender) {
    /// let view = PyDiagramView::from_render(snapshot);
    /// assert!(view.diagrams().is_empty());
    /// # }
    /// ```
    #[must_use]
    pub fn from_render(snapshot: PyDiagramRender) -> Self {
        Self::new(Source::Captured(snapshot), None)
    }

    /// Slots are `<!--linnet-svg-0-->`, `<!--linnet-svg-1-->`, and so on.
    /// Existing snapshots are reframed without rerunning layout or measurement.
    /// Other clones retain their original wrappers.
    ///
    /// ```no_run
    /// use linnet_render_py::PyDiagramView;
    /// # fn wrap(view: PyDiagramView) {
    /// let view = view.with_html("<figure><!--linnet-svg-0--></figure>".into());
    /// # }
    /// ```
    ///
    /// # Panics
    /// Panics if a completed collection lacks its completed child snapshots;
    /// that would violate the collection completion invariant.
    #[must_use]
    pub fn with_html(mut self, template: String) -> Self {
        self.html_template = Some(template);
        let previous = self.completed.get().map(|done| Completed {
            metrics: done.metrics.clone(),
            snapshot: self.reframe(done.snapshot.clone(), false),
        });
        let vector = self
            .vector
            .get()
            .map(|snapshot| self.reframe(snapshot.clone(), true));
        self.completed = Arc::new(OnceLock::new());
        self.vector = Arc::new(OnceLock::new());
        if let Some(done) = previous {
            let _ = self.completed.set(done);
        }
        if let Some(snapshot) = vector {
            let _ = self.vector.set(snapshot);
        }
        self
    }

    /// Configure an owned SVG transformation before either completion route runs.
    ///
    /// ```no_run
    /// use linnet_render_py::PyDiagramView;
    /// use pyo3::PyResult;
    /// # fn style(view: PyDiagramView) -> PyResult<PyDiagramView> {
    /// view.map_svg(|svg| svg.replacen("<svg ", "<svg class=\"highlighted\" ", 1))
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns `ValueError` if this view is already measured or outlined. The
    /// transformation must preserve SVG validity and any geometry it needs.
    pub fn map_svg(
        mut self,
        f: impl Fn(String) -> String + Send + Sync + 'static,
    ) -> PyResult<Self> {
        if self.completed.get().is_some() || self.vector.get().is_some() {
            return Err(PyValueError::new_err(
                "configure SVG styling before measurement or outlined export",
            ));
        }
        self.maps.push(Arc::new(f));
        self.completed = Arc::new(OnceLock::new());
        self.vector = Arc::new(OnceLock::new());
        Ok(self)
    }

    // Changing an HTML wrapper never repeats native layout or label compilation.
    fn reframe(&self, snapshot: PyDiagramRender, vector: bool) -> PyDiagramRender {
        if let Source::Views(views) = self.source.as_ref() {
            let diagrams = views
                .iter()
                .map(|view| {
                    if vector {
                        view.vector
                            .get()
                            .expect("an outlined collection has outlined every child")
                            .clone()
                    } else {
                        view.completed
                            .get()
                            .expect("a measured collection has completed every child")
                            .snapshot
                            .clone()
                    }
                })
                .collect();
            return self
                .present(diagrams, true)
                .expect("reframing captured pages does not compile or measure them");
        }
        let html = self
            .html_template
            .as_ref()
            .expect("with_html installs a template before reframing")
            .replace("<!--linnet-svg-0-->", &snapshot.pages().join("\n"));
        snapshot.with_html(html)
    }

    fn labels<'a>(&'a self, labels: &mut Vec<&'a LabelSource>) {
        match self.source.as_ref() {
            Source::Scenes(scenes) => {
                for scene in scenes {
                    labels.extend(scene.title.iter().chain(&scene.pages));
                }
            }
            Source::Views(views) => {
                for view in views {
                    view.labels(labels);
                }
            }
            Source::Captured(_) => {}
        }
    }

    fn label_count(&self) -> usize {
        match self.source.as_ref() {
            Source::Scenes(scenes) => scenes
                .iter()
                .map(|scene| scene.pages.len() + usize::from(scene.title.is_some()))
                .sum(),
            Source::Views(views) => views.iter().map(Self::label_count).sum(),
            Source::Captured(_) => 0,
        }
    }

    fn font_resource(py: Python<'_>) -> PyResult<&'static FontResource> {
        static FONT: OnceLock<FontResource> = OnceLock::new();
        if let Some(font) = FONT.get() {
            return Ok(font);
        }
        let bytes = STIX_TWO_MATH
            .bytes()
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
        let encoded: String = py
            .import("base64")?
            .call_method1("b64encode", (PyBytes::new(py, &bytes),))?
            .call_method1("decode", ("ascii",))?
            .extract()?;
        let css = format!(
            "/* {} */\n@font-face{{font-family:\"Linnet Offline Math\";src:url(\"data:font/woff2;base64,{encoded}\") format(\"woff2\");font-display:block}}\n.linnet-math-label,.linnet-math-label math{{font-family:\"Linnet Offline Math\";font-weight:400;font-style:normal;letter-spacing:normal}}",
            STIX_TWO_MATH_LICENSE,
        );
        // Include the rules as well as the font bytes in browser cache identity.
        let key = py
            .import("hashlib")?
            .call_method1("sha256", (PyBytes::new(py, css.as_bytes()),))?
            .call_method0("hexdigest")?
            .extract()?;
        let _ = FONT.set(FontResource { css, key });
        Ok(FONT.get().expect("font resource was initialized"))
    }

    fn python<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyModule>> {
        static MODULE: PyOnceLock<Py<PyModule>> = PyOnceLock::new();
        MODULE
            .get_or_try_init(py, || {
                let source =
                    CString::new(include_str!("view.py")).expect("Python source has no NUL");
                PyModule::from_code(
                    py,
                    &source,
                    c"linnet_diagram_view.py",
                    c"_linnet_diagram_view",
                )
                .map(Bound::unbind)
            })
            .map(|module| module.bind(py).clone())
    }

    fn receipt(&self, py: Python<'_>, value: &Bound<'_, PyDict>) -> PyResult<Vec<LabelMetrics>> {
        let json: String = py
            .import("json")?
            .call_method1("dumps", (value,))?
            .extract()?;
        let receipt: Receipt = serde_json::from_str(&json).map_err(|error| {
            PyValueError::new_err(format!("invalid label measurements: {error}"))
        })?;
        if receipt.id != self.id {
            return Err(PyValueError::new_err(
                "stale label measurements: request id does not match this view",
            ));
        }
        if receipt.metrics.len() != self.label_count() {
            return Err(PyValueError::new_err(
                "label measurements must cover every requested label exactly once",
            ));
        }
        receipt
            .metrics
            .into_iter()
            .enumerate()
            .map(|(index, label)| {
                if label.id != index {
                    return Err(PyValueError::new_err(
                        "label measurements must contain every label id in request order",
                    ));
                }
                let metrics = LabelMetrics {
                    width: label.width,
                    height: label.height,
                    baseline: label.baseline,
                };
                metrics
                    .validate()
                    .map_err(|error| PyValueError::new_err(error.to_string()))?;
                Ok(metrics)
            })
            .collect()
    }

    fn complete(
        &self,
        py: Python<'_>,
        metrics: &[LabelMetrics],
        css: &str,
    ) -> PyResult<PyDiagramRender> {
        if self.completed.get().is_none() {
            let snapshot = self.capture(py, Some(metrics), css)?;
            // Callbacks can reenter completion. The first published snapshot
            // wins, and every caller must observe that same immutable result.
            let _ = self.completed.set(Completed {
                metrics: metrics.to_vec(),
                snapshot,
            });
        }
        let done = self
            .completed
            .get()
            .expect("this call or a reentrant/concurrent caller published completion");
        if done.metrics.as_slice() != metrics {
            return Err(PyValueError::new_err(
                "view already has different label measurements; create a new view to measure again",
            ));
        }
        Ok(done.snapshot.clone())
    }

    fn capture(
        &self,
        py: Python<'_>,
        metrics: Option<&[LabelMetrics]>,
        css: &str,
    ) -> PyResult<PyDiagramRender> {
        let mut offset = 0;
        let (diagrams, collection) = match self.source.as_ref() {
            Source::Scenes(scenes) => {
                let mut figures = Vec::with_capacity(scenes.len());
                for scene in scenes {
                    let count = scene.pages.len() + usize::from(scene.title.is_some());
                    let typeset = if let Some(metrics) = metrics {
                        scene
                            .typeset_mathml(&metrics[offset..offset + count])
                            .map_err(|error| PyValueError::new_err(error.to_string()))?
                    } else if count == 0 {
                        scene
                            .typeset(&[])
                            .map_err(|error| PyValueError::new_err(error.to_string()))?
                    } else {
                        let source = scene
                            .label_document()
                            .map_err(|error| PyValueError::new_err(error.to_string()))?;
                        let pages = typst_export::compile_sources(
                            py,
                            &BTreeMap::from([("main.typ".into(), source.into_bytes())]),
                            "svg",
                        )?
                        .into_iter()
                        .map(|bytes| {
                            String::from_utf8(bytes)
                                .map_err(|error| PyRuntimeError::new_err(error.to_string()))
                        })
                        .collect::<PyResult<Vec<_>>>()?;
                        scene
                            .typeset(&pages)
                            .map_err(|error| PyValueError::new_err(error.to_string()))?
                    };
                    offset += count;
                    let svg = scene.render(&typeset).map_err(PyRuntimeError::new_err)?;
                    // Identity, focus and selection belong to each channel's SVG.
                    // The combined outer viewport must not own their targets.
                    figures.push(Scene::interactive_svg(&svg).map_err(PyRuntimeError::new_err)?);
                }
                let mut svg = Scene::combine_svgs(&figures).map_err(PyRuntimeError::new_err)?;
                if scenes.len() > 1 {
                    svg = svg.replacen(
                        "<svg ",
                        "<svg style=\"display:block;max-width:100%;height:auto\" ",
                        1,
                    );
                }
                if !css.is_empty() {
                    svg = svg.replacen(
                        '>',
                        &format!("><style data-linnet-font=\"true\"><![CDATA[{css}]]></style>"),
                        1,
                    );
                }
                (
                    vec![PyDiagramRender::new(
                        svg.clone(),
                        format!("<figure>{svg}</figure>"),
                    )],
                    false,
                )
            }
            Source::Views(views) => {
                let mut diagrams = Vec::with_capacity(views.len());
                for view in views {
                    let count = view.label_count();
                    diagrams.push(if let Some(metrics) = metrics {
                        view.complete(py, &metrics[offset..offset + count], css)?
                    } else {
                        view.outlined(py)?
                    });
                    offset += count;
                }
                (diagrams, true)
            }
            Source::Captured(snapshot) => {
                if metrics.is_none()
                    && snapshot
                        .pages()
                        .iter()
                        .any(|page| page.contains("linnet-math-label"))
                {
                    return Err(PyValueError::new_err("captured MathML pages have no retained scene for outlining; call outlined() on their original DiagramView"));
                }
                (vec![snapshot.clone()], false)
            }
        };
        self.present(diagrams, collection)
    }

    // Estimates never publish completion, including for retained child views.
    fn estimated(&self, py: Python<'_>) -> PyResult<PyDiagramRender> {
        if let Some(done) = self.completed.get() {
            return Ok(done.snapshot.clone());
        }
        if let Source::Views(views) = self.source.as_ref() {
            let diagrams = views
                .iter()
                .map(|view| view.estimated(py))
                .collect::<PyResult<Vec<_>>>()?;
            return self.present(diagrams, true);
        }
        let mut labels = Vec::new();
        self.labels(&mut labels);
        let metrics = labels
            .iter()
            .map(|label| {
                label
                    .estimated_metrics()
                    .map_err(|error| PyValueError::new_err(error.to_string()))
            })
            .collect::<PyResult<Vec<_>>>()?;
        let css = if metrics.is_empty() {
            ""
        } else {
            &Self::font_resource(py)?.css
        };
        self.capture(py, Some(&metrics), css)
    }

    fn present(
        &self,
        diagrams: Vec<PyDiagramRender>,
        collection: bool,
    ) -> PyResult<PyDiagramRender> {
        let diagrams = diagrams
            .into_iter()
            .map(|mut diagram| {
                for map in &self.maps {
                    diagram = diagram.map_svg(|svg| map(svg))?;
                }
                Ok(diagram)
            })
            .collect::<PyResult<Vec<_>>>()?;
        let html = if let Some(template) = &self.html_template {
            let mut html = template.clone();
            for (index, diagram) in diagrams.iter().enumerate() {
                html = html.replace(
                    &format!("<!--linnet-svg-{index}-->"),
                    &diagram.pages().join("\n"),
                );
            }
            html
        } else {
            diagrams
                .iter()
                .map(PyDiagramRender::to_html)
                .collect::<PyResult<Vec<_>>>()?
                .join("\n")
        };
        if collection {
            PyDiagramRender::from_diagrams(diagrams, html)
        } else {
            Ok(diagrams
                .into_iter()
                .next()
                .expect("one captured drawing")
                .with_html(html))
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl PyDiagramView {
    /// Configured child views for a collection, retaining their original scenes.
    /// See [`Self`] for an example. Returned clones share child completion state.
    #[getter]
    pub fn diagrams(&self) -> Vec<PyDiagramView> {
        match self.source.as_ref() {
            Source::Views(views) => views.clone(),
            _ => Vec::new(),
        }
    }

    /// The exact ordered MathML probes, fixed font, and receipt identity.
    ///
    /// See [`Self`] for the measurement lifecycle. A receipt must echo this
    /// request's `id` and contain one finite, nonnegative rectangle per label in
    /// request order. All dimensions and the baseline are SVG user units (points).
    ///
    /// # Errors
    /// Returns a Python error if bundled font loading, hashing, or request
    /// serialization fails. No browser or notebook dependency is imported here.
    pub fn measurement_request(&self, py: Python<'_>) -> PyResult<Py<PyDict>> {
        let mut labels = Vec::new();
        self.labels(&mut labels);
        let labels = labels
            .iter()
            .enumerate()
            .map(|(id, source)| json!({"id": id, "mathml": source.mathml}))
            .collect::<Vec<_>>();
        let font = if labels.is_empty() {
            None
        } else {
            Some(Self::font_resource(py)?)
        };
        let request = json!({
            "id": self.id,
            "labels": labels,
            "font_css": font.map_or("", |font| font.css.as_str()),
            "font_key": font.map_or("", |font| font.key.as_str()),
            "font_style": {"fontFamily": "Linnet Offline Math", "fontSize": "9px", "lineHeight": "1", "whiteSpace": "nowrap", "fontWeight": "400", "fontStyle": "normal", "letterSpacing": "normal"},
            "coordinate_system": "svg-user-units",
        });
        Ok(py
            .import("json")?
            .call_method1("loads", (request.to_string(),))?
            .extract()?)
    }

    /// Obtain a completed immutable snapshot, or explicitly measure using Chromium.
    ///
    /// See [`Self`] for an example. Calling without arguments reuses completion;
    /// pending labels need either an external receipt or explicit browser mode.
    /// The same receipt is idempotent, including through a retained clone.
    ///
    /// # Errors
    /// Returns `ValueError` for conflicting arguments, an unsupported browser,
    /// stale/malformed receipts, changed measurements, or invalid label sources.
    /// Returns `RuntimeError` while measurement is pending or native rendering
    /// fails. Optional browser setup and execution errors propagate as Python
    /// errors. Failed capture does not publish this view's completed snapshot;
    /// child views already completed remain independently usable.
    #[pyo3(signature = (measurements=None, *, browser=None))]
    pub fn measured(
        &self,
        py: Python<'_>,
        measurements: Option<&Bound<'_, PyDict>>,
        browser: Option<&str>,
    ) -> PyResult<PyDiagramRender> {
        if measurements.is_some() && browser.is_some() {
            return Err(PyValueError::new_err(
                "pass measurements or browser, not both",
            ));
        }
        if let Some(browser) = browser {
            if browser != "chromium" {
                return Err(PyValueError::new_err(
                    "the headless measurement browser must be 'chromium'",
                ));
            }
        }
        let metrics = if let Some(measurements) = measurements {
            self.receipt(py, measurements)?
        } else if let Some(done) = self.completed.get() {
            return Ok(done.snapshot.clone());
        } else if self.label_count() == 0 {
            Vec::new()
        } else if browser.is_some() {
            let request = self.measurement_request(py)?;
            let receipt =
                Self::python(py)?.call_method1("headless", (request, include_str!("view.js")))?;
            self.receipt(py, receipt.cast::<PyDict>()?)?
        } else {
            return Err(PyRuntimeError::new_err("diagram measurement is pending; display this view in a live notebook, pass measurement_request() results to measured(measurements=...), or explicitly use measured(browser='chromium')"));
        };
        let css = if metrics.is_empty() {
            ""
        } else {
            Self::font_resource(py)?.css.as_str()
        };
        self.complete(py, &metrics, css)
    }

    /// Outline labels with the optional Python Typst compiler; no browser is used.
    ///
    /// See [`Self`] for an example. The vector snapshot has its own completion
    /// cell, shared by clones, and never changes an existing measured export.
    ///
    /// # Errors
    /// Returns `ValueError` if labels lack Typst sources or captured MathML has
    /// no retained scene, `ImportError` if the optional `typst` package is absent,
    /// or a Python error on compilation/rendering failure. Compiler exceptions
    /// retain their Python cause via [`typst_export::compile_sources`].
    pub fn outlined(&self, py: Python<'_>) -> PyResult<PyDiagramRender> {
        if self.vector.get().is_none() {
            let snapshot = self.capture(py, None, "")?;
            let _ = self.vector.set(snapshot);
        }
        Ok(self
            .vector
            .get()
            .expect("this call or a reentrant/concurrent caller published outlined completion")
            .clone())
    }

    /// Static export owner: completed measurement, optional Typst, then estimates.
    ///
    /// Estimates use [`LabelSource::estimated_metrics`] at the request's fixed
    /// 9px STIX style. They never complete measurement; a later receipt replaces
    /// them. See [`Self::to_svg`] for an export example.
    ///
    /// # Errors
    /// Propagates label validation, optional compiler, and rendering errors.
    #[pyo3(signature = (*, browser=None))]
    pub fn snapshot(&self, py: Python<'_>, browser: Option<&str>) -> PyResult<PyDiagramRender> {
        if browser.is_some() || self.completed.get().is_some() {
            return self.measured(py, None, browser);
        }
        if Self::python(py)?
            .call_method1("available", ("typst",))?
            .extract::<bool>()?
        {
            self.outlined(py)
        } else {
            self.estimated(py)
        }
    }

    /// Return the completed single-page SVG; see [`Self`] for an example.
    ///
    /// # Errors
    /// Propagates [`Self::snapshot`] errors and rejects multipage snapshots.
    pub fn to_svg(&self, py: Python<'_>) -> PyResult<String> {
        self.snapshot(py, None)?.to_svg()
    }
    /// Return independent strings for all completed pages; see [`Self`].
    ///
    /// # Errors
    /// Propagates [`Self::snapshot`] errors if completion is unavailable.
    pub fn to_svg_pages(&self, py: Python<'_>) -> PyResult<Vec<String>> {
        self.snapshot(py, None)?.to_svg_pages()
    }
    /// Return completed HTML with its font resource; see [`Self`].
    ///
    /// # Errors
    /// Propagates [`Self::snapshot`] errors if completion is unavailable.
    pub fn to_html(&self, py: Python<'_>) -> PyResult<String> {
        self.snapshot(py, None)?.to_html()
    }
    /// Save the completed snapshot, selecting the format by suffix.
    ///
    /// ```no_run
    /// use linnet_render_py::PyDiagramView;
    /// use pyo3::{PyResult, Python};
    /// # fn save(py: Python<'_>, completed: PyDiagramView) -> PyResult<()> {
    /// completed.save(py, "diagram.svg".into())?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Propagates [`Self::snapshot`] errors, filesystem/export errors, and
    /// invalid format/page-count errors. PDF/PNG need an outlined snapshot.
    pub fn save(&self, py: Python<'_>, output: PathBuf) -> PyResult<PathBuf> {
        Py::new(py, self.snapshot(py, None)?)?
            .call_method1(py, "save", (output,))?
            .extract(py)
    }

    /// Optional notebook dependencies are imported only when displaying a view.
    /// The widget retains a clone of the backend owner until it is closed.
    ///
    /// ```no_run
    /// use linnet_render_py::PyDiagramView;
    /// use pyo3::{PyResult, Python};
    /// # fn notebook(py: Python<'_>, view: PyDiagramView) -> PyResult<()> {
    /// let widget = view.widget(py)?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Propagates font/request errors and Python import or widget setup errors.
    pub fn widget(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let completed = Some(self.to_html(py)?);
        let request = self.measurement_request(py)?;
        // Captured children can already contain MathML without pending probes.
        if completed
            .as_ref()
            .is_some_and(|html| html.contains("<style data-linnet-font=\"true\">"))
        {
            let font = Self::font_resource(py)?;
            request.bind(py).set_item("font_css", &font.css)?;
            request.bind(py).set_item("font_key", &font.key)?;
        }
        Self::python(py)?
            .call_method1(
                "widget",
                (
                    Py::new(py, self.clone())?,
                    request,
                    include_str!("view.js"),
                    completed,
                ),
            )
            .map(Bound::unbind)
    }

    /// Jupyter MIME bundle hook; see [`Self::widget`] for notebook construction.
    ///
    /// # Errors
    /// Propagates widget setup and Python MIME protocol errors.
    #[pyo3(signature = (include=None, exclude=None))]
    pub fn _repr_mimebundle_(
        &self,
        py: Python<'_>,
        include: Option<Vec<String>>,
        exclude: Option<Vec<String>>,
    ) -> PyResult<Py<PyAny>> {
        let wanted = |mime: &str| {
            include
                .as_ref()
                .is_none_or(|values| values.iter().any(|value| value == mime))
                && !exclude
                    .as_ref()
                    .is_some_and(|values| values.iter().any(|value| value == mime))
        };
        let bundle = PyDict::new(py);
        let snapshot = self.snapshot(py, None)?;
        if wanted("text/plain") {
            bundle.set_item("text/plain", self.__repr__())?;
        }
        if wanted("text/html") {
            bundle.set_item("text/html", snapshot.to_html()?)?;
        }
        if wanted("image/svg+xml") && snapshot.pages().len() == 1 {
            bundle.set_item("image/svg+xml", snapshot.to_svg()?)?;
        }
        if wanted("application/vnd.jupyter.widget-view+json")
            && Self::python(py)?
                .call_method1("available", ("anywidget",))?
                .extract::<bool>()?
        {
            let live = self.widget(py)?.call_method0(py, "_repr_mimebundle_")?;
            let live = live.bind(py).cast::<PyDict>()?;
            if let Some(value) = live.get_item("application/vnd.jupyter.widget-view+json")? {
                bundle.set_item("application/vnd.jupyter.widget-view+json", value)?;
            }
        }
        Ok(bundle.into_any().unbind())
    }

    /// Marimo adapter for the retained widget; see [`Self::widget`].
    ///
    /// # Errors
    /// Propagates widget errors and missing/failed Marimo imports or adapters.
    pub fn _display_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let marimo = py.import("marimo")?;
        if Self::python(py)?
            .call_method1("available", ("anywidget",))?
            .extract::<bool>()?
        {
            marimo
                .getattr("ui")?
                .call_method1("anywidget", (self.widget(py)?,))
                .map(Bound::unbind)
        } else {
            marimo
                .call_method1("Html", (self.to_html(py)?,))
                .map(Bound::unbind)
        }
    }
    /// Marimo MIME hook for [`Self::_display_`]; see [`Self::widget`].
    ///
    /// # Errors
    /// Propagates display adapter and Python MIME protocol errors.
    pub fn _mime_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self._display_(py)?.call_method0(py, "_mime_")
    }
    fn __repr__(&self) -> &'static str {
        if self.completed.get().is_some() {
            "DiagramView(measured)"
        } else {
            "DiagramView(awaiting browser measurement)"
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use linnest::svg::{Dash, Details, EdgeDrawing, NodeDrawing, Stroke};
    use roxmltree::Document;
    use std::sync::Weak;

    #[test]
    fn reentrant_completion_returns_the_first_published_snapshot() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            for mode in ["measured", "outlined"] {
                let owner = Arc::new(OnceLock::<Weak<PyDiagramView>>::new());
                let callback_owner = Arc::clone(&owner);
                let calls = Arc::new(AtomicU64::new(0));
                let callback_calls = Arc::clone(&calls);
                let svg = r#"<svg xmlns="http://www.w3.org/2000/svg"/>"#;
                let snapshot = PyDiagramRender::new(svg.into(), format!("<figure>{svg}</figure>"));
                let view = Arc::new(PyDiagramView::from_render(snapshot).map_svg(move |svg| {
                    let call = callback_calls.fetch_add(1, Ordering::Relaxed);
                    if call == 0 {
                        let view = callback_owner.get().unwrap().upgrade().unwrap();
                        Python::attach(|py| {
                            let snapshot = if mode == "measured" {
                                view.measured(py, None, None)
                            } else {
                                view.outlined(py)
                            }
                            .unwrap();
                            assert!(snapshot.to_svg().unwrap().contains("data-pass=\"1\""));
                        });
                    }
                    svg.replacen("<svg ", &format!("<svg data-pass=\"{call}\" "), 1)
                })?);
                assert!(owner.set(Arc::downgrade(&view)).is_ok());
                for _ in 0..2 {
                    let snapshot = if mode == "measured" {
                        view.measured(py, None, None)?
                    } else {
                        view.outlined(py)?
                    };
                    assert!(snapshot.to_svg()?.contains("data-pass=\"1\""));
                    assert!(!snapshot.to_svg()?.contains("data-pass=\"0\""));
                }
                assert_eq!(calls.load(Ordering::Relaxed), 2);
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn reentrant_measurement_rejects_the_losing_receipt() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let prepared = channel_view(py, 1)?;
            let Source::Scenes(mut scenes) = Arc::try_unwrap(prepared.source).ok().unwrap() else {
                panic!("fixture must own prepared scenes");
            };
            scenes[0].title = Some(LabelSource::text("receipt"));
            // Explicit external-metrics fixtures exercise receipt conflict,
            // independently of browser/font setup.
            let first = LabelMetrics {
                width: 12.0,
                height: 9.0,
                baseline: 7.0,
            };
            let winner = LabelMetrics {
                width: 13.0,
                ..first
            };
            let owner = Arc::new(OnceLock::<Weak<PyDiagramView>>::new());
            let callback_owner = Arc::clone(&owner);
            let calls = AtomicU64::new(0);
            let view = Arc::new(PyDiagramView::prepared(scenes).map_svg(move |svg| {
                if calls.fetch_add(1, Ordering::Relaxed) == 0 {
                    let view = callback_owner.get().unwrap().upgrade().unwrap();
                    Python::attach(|py| view.complete(py, &[winner], "")).unwrap();
                }
                svg
            })?);
            assert!(owner.set(Arc::downgrade(&view)).is_ok());
            let Err(error) = view.complete(py, &[first], "") else {
                panic!("a different winning receipt must reject the losing call");
            };
            assert!(error.is_instance_of::<PyValueError>(py));
            assert!(error.to_string().contains("different label measurements"));
            assert_eq!(view.completed.get().unwrap().metrics, [winner]);
            let saved = view.measured(py, None, None)?.to_svg()?;
            assert_eq!(view.complete(py, &[winner], "")?.to_svg()?, saved);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn measurement_request_uses_the_saved_svgs_local_font_and_coordinate_system() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            // The protocol is available without loading fonts or notebook tools.
            let view = PyDiagramView::prepared(Vec::new());
            let locals = PyDict::new(py);
            locals.set_item("request", view.measurement_request(py)?)?;
            py.run(
                cr#"
assert request['font_style']['fontSize'] == '9px'
assert request['font_style']['lineHeight'] == '1'
assert request['coordinate_system'] == 'svg-user-units'
assert 'points_per_css_pixel' not in request
assert request['labels'] == []
assert request['font_css'] == request['font_key'] == ''
"#,
                Some(&locals),
                Some(&locals),
            )?;
            Ok(())
        })
        .unwrap();
    }

    #[test]
    #[ignore = "requires the optional AnyWidget package and bundled font"]
    fn notebook_transport_reuses_its_class_and_preserves_complete_exports() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let original = channel_view(py, 1)?;
            let Source::Scenes(mut scenes) = Arc::try_unwrap(original.source).ok().unwrap() else {
                unreachable!()
            };
            scenes[0].title = Some(LabelSource::text("transport <font>"));
            let view = PyDiagramView::prepared(scenes);
            let request = view.measurement_request(py)?;
            let receipt = py.import("json")?.call_method1(
                "loads",
                (format!(
                    r#"{{"id":"{}","metrics":[{{"id":0,"width":80,"height":10,"baseline":8}}]}}"#,
                    view.id
                ),),
            )?;
            let snapshot = view.measured(py, Some(receipt.cast()?), None)?;
            let locals = PyDict::new(py);
            locals.set_item("request", request)?;
            locals.set_item("view", Py::new(py, view.clone())?)?;
            locals.set_item("widget", view.widget(py)?)?;
            locals.set_item("other", view.widget(py)?)?;
            locals.set_item("svg", snapshot.to_svg()?)?;
            locals.set_item("html", snapshot.to_html()?)?;
            locals.set_item("receipt", receipt)?;
            py.run(
                cr#"
import base64, hashlib, json
try:
    # AnyWidget adds an instance-specific trait subclass; its class owner is shared.
    assert type(widget).__base__ is type(other).__base__
    assert 'data:text/javascript;base64,' in widget._esm
    assert len(json.dumps(widget.request)) < 4096
    assert 'font_css' not in widget.request
    assert request['font_key'] == hashlib.sha256(request['font_css'].encode()).hexdigest()
    assert widget.font == other.font == {}
    widget._receive(widget, {'action': 'font', 'key': 'wrong'}, [])
    assert widget.font == {}
    widget._receive(widget, {'action': 'font', 'key': request['font_key']}, [])
    assert widget.font == {'key': request['font_key'], 'css': request['font_css']}
    block = '<style data-linnet-font="true"><![CDATA[' + request['font_css'] + ']]></style>'
    assert block in svg and block in html
    assert widget.rendered['html'] == html.replace(block, '')
    assert 'transport &lt;font&gt;' in widget.rendered['html']
    widget._receive(widget, {'action':'measure', 'measurements':receipt}, [])
    assert widget.rendered['html'] == html.replace(block, '')
    assert view.to_svg() == svg and view.to_html() == html
    # Exact owned-block removal preserves text and arbitrary other styles.
    literal = '<mtext>data-linnet-font="true"</mtext><style>mtext{color:red}</style>'
    assert widget._notebook_html(literal + block) == literal
    owner, comm = widget._view, widget.comm
    widget._receive(widget, {'action':'attach','view':'first'}, [])
    widget._receive(widget, {'action':'dispose','view':'first'}, [])
    assert widget._view is owner and widget.comm is comm
    assert widget._receive in widget._msg_callbacks.callbacks
finally:
    widget.close()
    other.close()
assert widget._view is None and widget._font_css == ''
"#,
                Some(&locals),
                Some(&locals),
            )?;
            Ok(())
        })
        .unwrap();
    }

    fn channel_view(py: Python<'_>, count: usize) -> PyResult<PyDiagramView> {
        let stroke = Stroke {
            paint: "#333333".into(),
            width: 0.5,
            dash: Dash::Solid,
            round_cap: false,
        };
        let node = NodeDrawing {
            radius: 0.2,
            label_padding: None,
            label: None,
            rectangular: false,
            fill: "none".into(),
            stroke: stroke.clone(),
            details: Details::default(),
        };
        let scenes = (0..count).map(|_| Scene {
                graph: serde_json::from_value(json!({
                    "nodes": [{"index": 0}, {"index": 1}],
                    "edges": [{"id": 0, "source": {"node": 0, "id": 0}, "sink": {"node": 1, "id": 1}}],
                })).unwrap(),
                nodes: vec![node.clone(), node.clone()],
                edges: vec![EdgeDrawing {
                    stroke: stroke.clone(),
                    pattern: None,
                    flow: None,
                    momentum: false,
                    momentum_offset: EdgeDrawing::DEFAULT_MOMENTUM_OFFSET,
                    label: None,
                    details: Details::default(),
                }],
                preamble: String::new(),
                title: None,
                pages: Vec::new(),
                layout: serde_json::from_value(json!({"layout-algo": "dot"})).unwrap(),
                layout_passes: Vec::new(),
                layout_edges: None,
                label_feedback: false,
            }).collect();
        PyDiagramView::from_scenes(py, scenes, None)
    }

    #[test]
    fn combined_channels_keep_separate_edge_identity_and_keyboard_focus() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let view = channel_view(py, 2)?;
            let svg = view.measured(py, None, None)?.to_svg()?;
            let document = Document::parse(&svg).unwrap();
            assert!(!document
                .root_element()
                .has_attribute("data-linnet-interactive"));
            assert!(document
                .root_element()
                .attribute("style")
                .unwrap()
                .contains("max-width:100%"));
            let owners = document
                .descendants()
                .filter(|node| node.has_attribute("data-linnet-interactive"))
                .collect::<Vec<_>>();
            assert_eq!(owners.len(), 2);
            for owner in owners {
                let targets = owner
                    .descendants()
                    .filter(|node| {
                        node.attribute("data-linnet-kind") == Some("edge")
                            && node.attribute("data-linnet-id") == Some("0")
                    })
                    .collect::<Vec<_>>();
                assert!(!targets.is_empty());
                assert_eq!(
                    targets
                        .iter()
                        .filter(|node| node.attribute("tabindex") == Some("0"))
                        .count(),
                    1
                );
                for target in targets {
                    assert_eq!(
                        target
                            .ancestors()
                            .find(|node| node.has_attribute("data-linnet-interactive")),
                        Some(owner)
                    );
                }
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    #[ignore = "requires the optional AnyWidget package"]
    fn retained_widget_keeps_its_backend_owner_and_comm_after_detach() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let view = channel_view(py, 1)?;
            let locals = PyDict::new(py);
            locals.set_item("widget", view.widget(py)?)?;
            py.run(
                cr#"
comm = widget.comm
assert comm is not None
handler = widget._receive
assert handler in widget._msg_callbacks.callbacks

def deliver(message):
    widget._msg_callbacks(widget, message, [])

try:
    deliver({'action': 'attach', 'view': 'first'})
    receipt = {'id': widget.request['id'], 'metrics': []}
    deliver({'action': 'measure', 'measurements': receipt})
    html = widget.rendered['html']
    deliver({'action': 'dispose', 'view': 'first'})
    assert widget._views == set()
    assert widget.comm is comm
    assert widget._view is not None
    assert handler in widget._msg_callbacks.callbacks

    deliver({'action': 'attach', 'view': 'second'})
    assert widget._views == {'second'}
    deliver({'action': 'measure', 'measurements': receipt})
    assert widget.rendered['html'] == html
    assert widget._view.measured().to_html() == html
    deliver({'action': 'dispose', 'view': 'second'})
    assert widget.comm is comm
finally:
    widget.close()
assert widget.comm is None
assert widget._view is None
assert handler not in widget._msg_callbacks.callbacks
"#,
                Some(&locals),
                Some(&locals),
            )?;
            Ok(())
        })
        .unwrap();
    }

    #[test]
    #[ignore = "requires Playwright and installed Chromium"]
    fn channel_bounds_and_widget_frame_height_survive_zoom_and_reset() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let documents = (1..=2).map(|count| {
                channel_view(py, count)?
                    .with_html("<figure style=\"margin:17px 0 29px\"><!--linnet-svg-0--></figure>".into())
                    .measured(py, None, None)?.to_html()
            }).collect::<PyResult<Vec<_>>>()?;
            let locals = PyDict::new(py);
            locals.set_item("documents", documents)?;
            locals.set_item("javascript", include_str!("view.js"))?;
            py.run(cr#"
import os
from playwright.sync_api import sync_playwright

with sync_playwright() as driver:
    browser = driver.chromium.launch(executable_path=os.environ.get('LINNET_TEST_CHROMIUM'))
    try:
        page = browser.new_page()
        for html in documents:
            page.evaluate('''async ({javascript, html}) => {
                const url = URL.createObjectURL(new Blob([javascript], {type: 'text/javascript'}));
                const el = document.createElement('div');
                el.style.width = '350px';
                document.body.append(el);
                let cleanup, timer;
                const assert = (value, message) => { if (!value) throw new Error(message); };
                const attributes = svg => ['x', 'y', 'width', 'height', 'viewBox', 'style']
                    .map(name => svg.getAttribute(name));
                try {
                    const module = await import(url);
                    const original = [...new DOMParser().parseFromString(html, 'text/html')
                        .querySelectorAll('svg[data-linnet-interactive]')].map(attributes);
                    const state = {request: {id: 'channels', labels: []}, rendered: {id: 'channels', html}};
                    const ready = new Promise((resolve, reject) => {
                        timer = setTimeout(() => reject(new Error('render event timed out')), 10000);
                        el.addEventListener('linnet:rendered', resolve, {once: true});
                    });
                    const model = {get: name => state[name], on() {}, off() {}, send() {}};
                    cleanup = await module.default.render({model, el});
                    const doc = (await ready).detail.document;
                    const settle = async () => {
                        for (let i = 0; i < 4; i++) await new Promise(resolve => doc.defaultView.requestAnimationFrame(resolve));
                    };
                    await settle();
                    assert(doc.body.hasAttribute('data-linnet-frame-owner'), 'widget must own frame sizing before controls initialize');
                    const owners = [...doc.querySelectorAll('svg[data-linnet-interactive]')];
                    assert(owners.length === original.length, 'lost channel controllers');
                    if (owners.length > 1) {
                        assert(JSON.stringify(owners.map(attributes)) === JSON.stringify(original), 'controller resized a composite channel');
                        const first = owners[0].getBoundingClientRect(), second = owners[1].getBoundingClientRect();
                        assert(first.right <= second.left, 'composite channels overlap');
                        const outer = owners[0].parentElement.getBoundingClientRect();
                        assert(second.right <= outer.right + .01 && second.bottom <= outer.bottom + .01, 'composite channel is clipped');
                    }
                    const before = owners.map(owner => owner.querySelector('.linnet-viewport').getAttribute('viewBox'));
                    for (const owner of owners) {
                        const viewport = owner.querySelector('.linnet-viewport');
                        viewport.dispatchEvent(new doc.defaultView.KeyboardEvent('keydown', {key: '+', bubbles: true}));
                        viewport.dispatchEvent(new doc.defaultView.KeyboardEvent('keydown', {key: '0', bubbles: true}));
                    }
                    await settle();
                    assert(JSON.stringify(owners.map(owner => owner.querySelector('.linnet-viewport').getAttribute('viewBox'))) === JSON.stringify(before), 'reset changed the fitted channel cameras');
                    if (owners.length > 1) {
                        assert(JSON.stringify(owners.map(attributes)) === JSON.stringify(original), 'reset changed composite geometry');
                        owners[0].querySelector('[data-linnet-kind="edge"][data-linnet-id="0"][tabindex="0"]')
                            .dispatchEvent(new doc.defaultView.MouseEvent('click', {bubbles: true, shiftKey: true}));
                        assert(owners[0].linnetSelection.edges.includes(0), 'first channel selection is inert');
                        assert(owners[1].linnetSelection.edges.length === 0, 'channel selection leaked');
                    }
                    await settle();
                    const figure = doc.querySelector('figure');
                    const bottom = figure.getBoundingClientRect().bottom + parseFloat(doc.defaultView.getComputedStyle(figure).marginBottom);
                    const height = el.querySelector('iframe').getBoundingClientRect().height;
                    assert(height >= bottom, 'a controller overwrote frame height and clipped the figure margin');
                } finally {
                    clearTimeout(timer);
                    cleanup?.();
                    el.remove();
                    URL.revokeObjectURL(url);
                }
            }''', {'javascript': javascript, 'html': html})
    finally:
        browser.close()
"#, Some(&locals), Some(&locals))?;
            Ok(())
        }).unwrap();
    }

    #[test]
    fn explicit_headless_measurement_joins_workers_with_an_active_asyncio_loop() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = CString::new(include_str!("view.py")).unwrap();
            let adapter = PyModule::from_code(py, &source, c"view.py", c"_view_worker_fixture")?;
            let locals = PyDict::new(py);
            locals.set_item("adapter", adapter)?;
            py.run(
                cr#"
import asyncio
import threading

main = threading.current_thread()
workers = []

def measure_fixture(request, javascript):
    assert javascript == 'source'
    current = threading.current_thread()
    if request == 'direct':
        assert current is main
    else:
        assert current is not main
        workers.append(current)
        try:
            asyncio.get_running_loop()
        except RuntimeError:
            pass
        else:
            raise AssertionError('worker inherited the caller event loop')
    if request == 'failure':
        raise RuntimeError('measurement failed')
    return request

adapter._headless = measure_fixture
assert adapter.headless('direct', 'source') == 'direct'

async def notebook():
    assert adapter.headless('success', 'source') == 'success'
    assert not workers[-1].is_alive()
    try:
        adapter.headless('failure', 'source')
    except RuntimeError as error:
        assert str(error) == 'measurement failed'
    else:
        raise AssertionError('measurement failure was swallowed')
    assert not workers[-1].is_alive()

asyncio.run(notebook())
assert len(workers) == 2
assert all(not worker.is_alive() for worker in workers)
"#,
                Some(&locals),
                Some(&locals),
            )?;
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn headless_browser_closes_after_success_and_evaluation_failure() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = CString::new(include_str!("view.py")).unwrap();
            let adapter =
                PyModule::from_code(py, &source, c"view.py", c"_view_browser_cleanup_fixture")?;
            let locals = PyDict::new(py);
            locals.set_item("adapter", adapter)?;
            py.run(
                cr#"
import sys
import types
from unittest.mock import patch

browsers = []
exits = []

class Browser:
    def __init__(self):
        self.closed = False
        browsers.append(self)
    def new_page(self):
        return self
    def route(self, pattern, handler):
        assert pattern == '**/*'
    def evaluate(self, javascript, payload):
        if payload['request'] == 'failure':
            raise RuntimeError('evaluation failed')
        return {'id': payload['request'], 'metrics': []}
    def close(self):
        self.closed = True

class Driver:
    chromium = types.SimpleNamespace(launch=Browser)
    def __enter__(self):
        return self
    def __exit__(self, *error):
        exits.append(error)

playwright = types.ModuleType('playwright')
sync_api = types.ModuleType('playwright.sync_api')
sync_api.sync_playwright = Driver
with patch.dict(sys.modules, {'playwright': playwright, 'playwright.sync_api': sync_api}):
    assert adapter.headless('success', 'source') == {'id': 'success', 'metrics': []}
    try:
        adapter.headless('failure', 'source')
    except RuntimeError as error:
        assert str(error) == 'evaluation failed'
    else:
        raise AssertionError('evaluation failure was swallowed')
assert len(browsers) == len(exits) == 2
assert all(browser.closed for browser in browsers)
"#,
                Some(&locals),
                Some(&locals),
            )?;
            Ok(())
        })
        .unwrap();
    }

    #[test]
    #[ignore = "requires Playwright and installed Chromium"]
    fn notebook_iframe_height_includes_outer_figure_margins() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let locals = PyDict::new(py);
            locals.set_item("javascript", include_str!("view.js"))?;
            py.run(cr#"
import os
from playwright.sync_api import sync_playwright

with sync_playwright() as driver:
    browser = driver.chromium.launch(executable_path=os.environ.get('LINNET_TEST_CHROMIUM'))
    try:
        page = browser.new_page()
        dimensions = page.evaluate('''async javascript => {
            const url = URL.createObjectURL(new Blob([javascript], {type: 'text/javascript'}));
            const el = document.createElement('div');
            document.body.append(el);
            let cleanup;
            let timer;
            try {
                const module = await import(url);
                const state = {
                    request: {id: 'margins', labels: []},
                    rendered: {id: 'margins', html: '<figure style="margin:17px 0 29px"><svg style="display:block" width="80" height="40"></svg></figure>'},
                };
                const ready = new Promise((resolve, reject) => {
                    timer = setTimeout(() => reject(new Error('render event timed out')), 10000);
                    el.addEventListener('linnet:rendered', resolve, {once: true});
                });
                const model = {get: name => state[name], on() {}, off() {}, send() {}};
                cleanup = await module.default.render({model, el});
                const event = await ready;
                const doc = event.detail.document;
                const figure = doc.querySelector('figure').getBoundingClientRect();
                return {
                    frame: el.querySelector('iframe').getBoundingClientRect().height,
                    body: doc.body.getBoundingClientRect().height,
                    bottom: figure.bottom + parseFloat(doc.defaultView.getComputedStyle(doc.querySelector('figure')).marginBottom),
                };
            } finally {
                clearTimeout(timer);
                cleanup?.();
                el.remove();
                URL.revokeObjectURL(url);
            }
        }''', javascript)
        assert dimensions['body'] == 86, dimensions
        assert dimensions['frame'] >= dimensions['bottom'], dimensions
    finally:
        browser.close()
"#, Some(&locals), Some(&locals))?;
            Ok(())
        }).unwrap();
    }
}
