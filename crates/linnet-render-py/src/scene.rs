//! Native graph drawing with Typst used only to typeset labels and titles.
use crate::{PyDiagramRender, PyRenderSettings};
use linnest::svg::Scene;
use pyo3::{
    exceptions::PyRuntimeError,
    prelude::*,
    types::{PyDict, PyDictMethods, PyList},
};
use std::collections::BTreeMap;

impl PyDiagramRender {
    pub fn with_html(mut self, html: String) -> Self {
        self.html = Some(html);
        self
    }
    pub fn map_svg(mut self, f: impl Fn(String) -> String) -> PyResult<Self> {
        let pages: Vec<String> = self.pages().iter().cloned().map(f).collect();
        self.pages = std::sync::Arc::new(pages);
        Ok(self)
    }
    pub fn from_scene(
        py: Python<'_>,
        mut scene: Scene,
        settings: Option<&PyRenderSettings>,
    ) -> PyResult<Self> {
        let settings = settings.cloned().unwrap_or_default();
        if settings.has_selectors() {
            let selector = SCENE_SELECTORS.get().ok_or_else(|| {
                PyRuntimeError::new_err(
                    "drawing selectors require the optional Linnet Python graph bindings",
                )
            })?;
            let elements = PyDict::new(py);
            for (name, details) in [
                (
                    "nodes",
                    scene
                        .nodes
                        .iter()
                        .map(|drawing| &drawing.details)
                        .collect::<Vec<_>>(),
                ),
                (
                    "edges",
                    scene
                        .edges
                        .iter()
                        .map(|drawing| &drawing.details)
                        .collect::<Vec<_>>(),
                ),
            ] {
                let values = PyList::empty(py);
                for details in details {
                    let drawing = PyDict::new(py);
                    let extensions = PyDict::new(py);
                    let data = py.import("json")?.call_method1(
                        "loads",
                        (serde_json::to_string(
                            &details
                                .fields()
                                .iter()
                                .cloned()
                                .collect::<serde_json::Map<_, _>>(),
                        )
                        .map_err(|error| PyRuntimeError::new_err(error.to_string()))?,),
                    )?;
                    extensions.set_item("inspection", data)?;
                    drawing.set_item("extensions", extensions)?;
                    values.append(drawing)?;
                }
                elements.set_item(name, values)?;
            }
            let transport = selector(py, &scene.graph, &elements, &settings)?;
            scene = transport.native_scene(scene.graph.clone(), None, Some(scene))?;
        } else {
            settings
                .config()?
                .apply(&mut scene)
                .map_err(PyRuntimeError::new_err)?;
        }
        let pages = if scene.pages.is_empty() && scene.title.is_none() {
            vec![]
        } else {
            typst_renderer::Document::compile_sources(
                &BTreeMap::from([("main.typ".into(), scene.label_document().into_bytes())]),
                "svg",
            )
            .map_err(PyRuntimeError::new_err)?
            .into_iter()
            .map(|b| String::from_utf8(b).map_err(|e| PyRuntimeError::new_err(e.to_string())))
            .collect::<PyResult<Vec<_>>>()?
        };
        let typeset = scene.typeset(&pages).map_err(PyRuntimeError::new_err)?;
        let svg = Scene::interactive_svg(&scene.render(&typeset).map_err(PyRuntimeError::new_err)?)
            .map_err(PyRuntimeError::new_err)?;
        Ok(Self::new(svg.clone(), format!("<figure>{svg}</figure>")))
    }
}

/// Optional graph bindings provide live views only when selectors need them.
pub type SceneSelectors = fn(
    Python<'_>,
    &linnest::TypstGraphSpec,
    &Bound<'_, PyDict>,
    &PyRenderSettings,
) -> PyResult<crate::typst::RenderSettingsTransport>;
static SCENE_SELECTORS: std::sync::OnceLock<SceneSelectors> = std::sync::OnceLock::new();
pub fn register_scene_selectors(selectors: SceneSelectors) {
    let _ = SCENE_SELECTORS.set(selectors);
}
