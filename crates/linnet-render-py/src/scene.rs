//! Resolve drawing settings once, before browser measurement or vector export.
use std::sync::{Arc, OnceLock};

use linnest::svg::Scene;
use linnest::TypstGraphSpec;
use pyo3::exceptions::PyRuntimeError;
use pyo3::types::{PyAnyMethods, PyDict, PyDictMethods, PyList, PyListMethods};
use pyo3::{Bound, PyResult, Python};
use serde_json::Map;

use crate::typst::RenderSettingsTransport;
use crate::{PyDiagramRender, PyDiagramView, PyRenderSettings};

impl PyDiagramRender {
    /// Replace the snapshot's wrapper without changing its captured SVG pages.
    ///
    /// See [`PyDiagramView::with_html`] for the equivalent preparation step.
    #[must_use]
    pub fn with_html(mut self, html: String) -> Self {
        self.html = Some(html);
        self
    }
    /// Transform independent page strings and update their occurrences in HTML.
    /// Existing clones retain their original captured pages.
    ///
    /// See [`PyDiagramView::map_svg`] for an SVG transformation example. This
    /// operation does not validate the resulting SVG; the result type is part
    /// of the external rendering API.
    pub fn map_svg(mut self, f: impl Fn(String) -> String) -> PyResult<Self> {
        let pages: Vec<String> = self.pages().iter().cloned().map(f).collect();
        if let Some(html) = &mut self.html {
            for (old, new) in self.pages.iter().zip(&pages) {
                *html = html.replace(old, new);
            }
        }
        self.pages = Arc::new(pages);
        Ok(self)
    }
}

impl PyDiagramView {
    /// Resolve one scene's settings and retain it for browser/vector completion.
    /// See [`Self`] for a complete measurement example.
    ///
    /// # Errors
    /// Returns a Python error for unsupported settings or failing selectors.
    /// Selectors require the optional Linnet graph bindings to be registered.
    pub fn from_scene(
        py: Python<'_>,
        scene: Scene,
        settings: Option<&PyRenderSettings>,
    ) -> PyResult<Self> {
        Self::from_scenes(py, vec![scene], settings)
    }

    /// Resolve each channel once; measurement later combines their native SVGs.
    ///
    /// See [`Self`] for the lifecycle. Multiple channels share one presentation
    /// while retaining separate graph identity and interaction targets.
    ///
    /// # Errors
    /// Returns the same configuration/selector errors as [`Self::from_scene`].
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
        let defaults = PyRenderSettings::default();
        let settings = settings.unwrap_or(&defaults);
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
                            &details.fields().iter().cloned().collect::<Map<_, _>>(),
                        )
                        .map_err(|error| PyRuntimeError::new_err(error.to_string()))?,),
                    )?;
                    extensions.set_item("inspection", data)?;
                    drawing.set_item("extensions", extensions)?;
                    values.append(drawing)?;
                }
                elements.set_item(name, values)?;
            }
            let transport = selector(py, &scene.graph, &elements, settings)?;
            // The transport needs its own editable topology while the consumed
            // scene supplies the original drawing defaults.
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
    &TypstGraphSpec,
    &Bound<'_, PyDict>,
    &PyRenderSettings,
) -> PyResult<RenderSettingsTransport>;
static SCENE_SELECTORS: OnceLock<SceneSelectors> = OnceLock::new();
pub fn register_scene_selectors(selectors: SceneSelectors) {
    let _ = SCENE_SELECTORS.set(selectors);
}
