//! Resolve drawing settings once, before browser measurement or vector export.
use crate::{PyDiagramRender, PyRenderSettings};
use linnest::svg::Scene;
use pyo3::{
    exceptions::PyRuntimeError,
    prelude::*,
    types::{PyDict, PyDictMethods, PyList},
};

impl PyDiagramRender {
    pub fn with_html(mut self, html: String) -> Self {
        self.html = Some(html);
        self
    }
    pub fn map_svg(mut self, f: impl Fn(String) -> String) -> PyResult<Self> {
        let pages: Vec<String> = self.pages().iter().cloned().map(f).collect();
        if let Some(html) = &mut self.html {
            for (old, new) in self.pages.iter().zip(&pages) {
                *html = html.replace(old, new);
            }
        }
        self.pages = std::sync::Arc::new(pages);
        Ok(self)
    }
}

impl crate::PyDiagramView {
    pub fn from_scene(
        py: Python<'_>,
        scene: Scene,
        settings: Option<&PyRenderSettings>,
    ) -> PyResult<Self> {
        Self::from_scenes(py, vec![scene], settings)
    }

    /// Resolve each channel once; measurement later combines their native SVGs.
    pub fn from_scenes(
        py: Python<'_>,
        scenes: Vec<Scene>,
        settings: Option<&PyRenderSettings>,
    ) -> PyResult<Self> {
        let scenes = scenes
            .into_iter()
            .map(|scene| Self::resolve_scene(py, scene, settings))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self::prepared(scenes))
    }

    fn resolve_scene(
        py: Python<'_>,
        mut scene: Scene,
        settings: Option<&PyRenderSettings>,
    ) -> PyResult<Scene> {
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
        Ok(scene)
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
