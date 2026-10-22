//! Prepared native scenes become immutable snapshots after real label measurement.
use std::{
    collections::BTreeMap,
    ffi::CString,
    path::PathBuf,
    sync::{
        atomic::{AtomicU64, Ordering},
        Arc, OnceLock,
    },
};

use linnest::svg::{LabelMetrics, LabelSource, Scene};
use pyo3::{
    exceptions::{PyRuntimeError, PyValueError},
    prelude::*,
    sync::PyOnceLock,
    types::{PyBytes, PyDict},
};
use serde::Deserialize;

use crate::PyDiagramRender;

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

/// A resolved drawing awaiting its notebook browser's label measurements.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
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
    pub fn from_views(views: Vec<Self>, html_template: Option<String>) -> PyResult<Self> {
        Ok(Self::new(Source::Views(views), html_template))
    }

    /// Include a captured render without reconstructing or measuring its geometry.
    pub fn from_render(snapshot: PyDiagramRender) -> Self {
        Self::new(Source::Captured(snapshot), None)
    }

    /// Slots are `<!--linnet-svg-0-->`, `<!--linnet-svg-1-->`, and so on.
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
                        view.vector.get().expect("completed child vector").clone()
                    } else {
                        view.completed
                            .get()
                            .expect("completed child measurement")
                            .snapshot
                            .clone()
                    }
                })
                .collect();
            return self.present(diagrams, true).expect("captured presentation");
        }
        let html = self
            .html_template
            .as_ref()
            .expect("HTML template")
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

    fn font_css(py: Python<'_>) -> PyResult<String> {
        static CSS: OnceLock<String> = OnceLock::new();
        if let Some(css) = CSS.get() {
            return Ok(css.clone());
        }
        let bytes = crate::data_assets::STIX_TWO_MATH
            .bytes()
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
        let encoded: String = py
            .import("base64")?
            .call_method1("b64encode", (PyBytes::new(py, &bytes),))?
            .call_method1("decode", ("ascii",))?
            .extract()?;
        let css = format!(
            "/* {} */\n@font-face{{font-family:\"Linnet Offline Math\";src:url(\"data:font/woff2;base64,{encoded}\") format(\"woff2\");font-display:block}}\n.linnet-math-label,.linnet-math-label math{{font-family:\"Linnet Offline Math\";font-weight:400;font-style:normal;letter-spacing:normal}}",
            crate::data_assets::STIX_TWO_MATH_LICENSE,
        );
        let _ = CSS.set(css.clone());
        Ok(css)
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
                metrics.validate().map_err(PyValueError::new_err)?;
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
        if let Some(done) = self.completed.get() {
            if done.metrics.as_slice() != metrics {
                return Err(PyValueError::new_err("view already has different label measurements; create a new view to measure again"));
            }
            return Ok(done.snapshot.clone());
        }
        let snapshot = self.capture(py, Some(metrics), css)?;
        let _ = self.completed.set(Completed {
            metrics: metrics.to_vec(),
            snapshot: snapshot.clone(),
        });
        Ok(snapshot)
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
                            .map_err(PyValueError::new_err)?
                    } else if count == 0 {
                        scene.typeset(&[]).map_err(PyValueError::new_err)?
                    } else {
                        let source = scene.label_document().map_err(PyValueError::new_err)?;
                        let pages = crate::typst_export::compile_sources(
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
                        scene.typeset(&pages).map_err(PyValueError::new_err)?
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
                    &diagram.to_svg_pages()?.join("\n"),
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

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyDiagramView {
    /// Configured child views for a collection, retaining their original scenes.
    #[getter]
    pub fn diagrams(&self) -> Vec<PyDiagramView> {
        match self.source.as_ref() {
            Source::Views(views) => views.clone(),
            _ => Vec::new(),
        }
    }

    /// The exact ordered MathML probes, fixed font, and receipt identity.
    pub fn measurement_request(&self, py: Python<'_>) -> PyResult<Py<PyDict>> {
        let mut labels = Vec::new();
        self.labels(&mut labels);
        let labels = labels
            .iter()
            .enumerate()
            .map(|(id, source)| serde_json::json!({"id": id, "mathml": source.mathml}))
            .collect::<Vec<_>>();
        let css = if labels.is_empty() {
            String::new()
        } else {
            Self::font_css(py)?
        };
        let request = serde_json::json!({
            "id": self.id,
            "labels": labels,
            "font_css": css,
            "font_style": {"fontFamily": "Linnet Offline Math", "fontSize": "9px", "lineHeight": "1", "whiteSpace": "nowrap", "fontWeight": "400", "fontStyle": "normal", "letterSpacing": "normal"},
            "coordinate_system": "svg-user-units",
        });
        Ok(py
            .import("json")?
            .call_method1("loads", (request.to_string(),))?
            .extract()?)
    }

    /// Obtain a completed immutable snapshot, or explicitly measure using Chromium.
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
            String::new()
        } else {
            Self::font_css(py)?
        };
        self.complete(py, &metrics, &css)
    }

    /// Outline labels with the optional Python Typst compiler; no browser is used.
    pub fn outlined(&self, py: Python<'_>) -> PyResult<PyDiagramRender> {
        if let Some(snapshot) = self.vector.get() {
            return Ok(snapshot.clone());
        }
        let snapshot = self.capture(py, None, "")?;
        let _ = self.vector.set(snapshot.clone());
        Ok(snapshot)
    }

    pub fn to_svg(&self, py: Python<'_>) -> PyResult<String> {
        self.measured(py, None, None)?.to_svg()
    }
    pub fn to_svg_pages(&self, py: Python<'_>) -> PyResult<Vec<String>> {
        self.measured(py, None, None)?.to_svg_pages()
    }
    pub fn to_html(&self, py: Python<'_>) -> PyResult<String> {
        self.measured(py, None, None)?.to_html()
    }
    pub fn save(&self, py: Python<'_>, output: PathBuf) -> PyResult<PathBuf> {
        Py::new(py, self.measured(py, None, None)?)?
            .call_method1(py, "save", (output,))?
            .extract(py)
    }

    /// Optional notebook dependencies are imported only when displaying a view.
    pub fn widget(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let completed = self
            .completed
            .get()
            .map(|done| done.snapshot.to_html())
            .transpose()?;
        Self::python(py)?
            .call_method1(
                "widget",
                (
                    Py::new(py, self.clone())?,
                    self.measurement_request(py)?,
                    include_str!("view.js"),
                    completed,
                ),
            )
            .map(Bound::unbind)
    }

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
        if !wanted("application/vnd.jupyter.widget-view+json") {
            let bundle = PyDict::new(py);
            if wanted("text/plain") {
                bundle.set_item("text/plain", self.__repr__())?;
            }
            if wanted("text/html") {
                if let Some(done) = self.completed.get() {
                    bundle.set_item("text/html", done.snapshot.to_html()?)?;
                }
            }
            return Ok(bundle.into_any().unbind());
        }
        let kwargs = PyDict::new(py);
        kwargs.set_item("include", include)?;
        kwargs.set_item("exclude", exclude)?;
        self.widget(py)?
            .call_method(py, "_repr_mimebundle_", (), Some(&kwargs))
    }

    pub fn _display_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        py.import("marimo")?
            .getattr("ui")?
            .call_method1("anywidget", (self.widget(py)?,))
            .map(Bound::unbind)
    }
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
                graph: serde_json::from_value(serde_json::json!({
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
                layout: serde_json::from_value(serde_json::json!({"layout-algo": "dot"})).unwrap(),
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
            let document = roxmltree::Document::parse(&svg).unwrap();
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
