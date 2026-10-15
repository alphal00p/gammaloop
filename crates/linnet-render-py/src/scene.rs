//! Domain scenes use the native renderer when possible and the common Typst pipeline otherwise.
use crate::{PyDiagramRender, PyRenderSettings};
use linnest::svg::{Dash, Pattern, Scene, Stroke};
use pyo3::{
    exceptions::PyRuntimeError,
    prelude::*,
    types::{PyDict, PyDictMethods, PyList},
};
use std::collections::BTreeMap;

fn stroke<'py>(py: Python<'py>, value: &Stroke) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("paint", crate::typst::scene_color(py, &value.paint)?)?;
    result.set_item("thickness", crate::typst::point_length(py, value.width)?)?;
    result.set_item(
        "dash",
        match value.dash {
            Dash::Solid => "solid",
            Dash::Dotted => "dotted",
            Dash::Dashed(..) => "dashed",
        },
    )?;
    if value.round_cap {
        result.set_item("cap", "round")?;
    }
    Ok(result)
}
impl PyDiagramRender {
    pub fn with_html(mut self, html: String) -> Self {
        self.html = Some(html);
        self
    }
    pub fn map_svg(mut self, f: impl Fn(String) -> String) -> PyResult<Self> {
        let pages: Vec<String> = self.pages()?.iter().cloned().map(f).collect();
        self.pages = std::sync::Arc::new(std::sync::OnceLock::from(pages));
        Ok(self)
    }
    pub fn from_scene(
        py: Python<'_>,
        mut scene: Scene,
        settings: Option<&PyRenderSettings>,
    ) -> PyResult<Self> {
        let settings = settings.cloned().unwrap_or_default();
        // Probe on a copy: a rejected late option must not partially alter fallback input.
        let probe = Scene {
            graph: scene.graph.clone(),
            nodes: scene.nodes.clone(),
            edges: scene.edges.clone(),
            preamble: scene.preamble.clone(),
            title: scene.title.clone(),
            pages: scene.pages.clone(),
            layout: scene.layout.clone(),
            layout_edges: scene.layout_edges.clone(),
            label_feedback: scene.label_feedback,
        };
        let native = settings.config().and_then(|config| {
            let mut probe = probe;
            config.apply(&mut probe).map_err(PyRuntimeError::new_err)?;
            if !config.template_options.is_empty() {
                return Err(PyRuntimeError::new_err("template options require Typst"));
            }
            Ok(probe)
        });
        if let Ok(native) = native {
            scene = native;
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
            let svg =
                Scene::interactive_svg(&scene.render(&typeset).map_err(PyRuntimeError::new_err)?)
                    .map_err(PyRuntimeError::new_err)?;
            return Ok(Self::new(svg.clone(), format!("<figure>{svg}</figure>")));
        }
        let dir = tempfile::tempdir().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        let labels = dir.path().join("labels.typ");
        let settings = settings.with_scene_defaults(py, &scene, labels.clone())?;
        let mut source = include_str!("../typst/scene.typ").to_owned();
        source.push_str(&scene.preamble);
        for (i, label) in scene.pages.iter().enumerate() {
            source.push_str(&format!(
                "\n#let label{i} = [\n{}\n#{label}\n]\n",
                scene.preamble
            ));
        }
        std::fs::write(&labels, source).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        let label = |drawing: &Bound<'_, PyDict>, index: Option<usize>| -> PyResult<()> {
            drawing.set_item(
                "label",
                match index {
                    Some(i) => {
                        crate::typst::label_reference(py, labels.clone(), &format!("label{i}"))?
                    }
                    None => py.None(),
                },
            )
        };
        let nodes = PyList::empty(py);
        let edges = PyList::empty(py);
        for node in &scene.nodes {
            let drawing = PyDict::new(py);
            label(&drawing, node.label)?;
            let style = PyDict::new(py);
            style.set_item(
                "fill",
                if node.fill == "url(#linnest-process-hatch)" {
                    crate::typst::scene_hatch(py, labels.clone())?
                } else {
                    crate::typst::scene_color(py, &node.fill)?
                },
            )?;
            let extensions = PyDict::new(py);
            extensions.set_item("scene-rectangular", node.rectangular)?;
            drawing.set_item("extensions", &extensions)?;
            style.set_item("radius", node.radius)?;
            style.set_item("stroke", stroke(py, &node.stroke)?)?;
            drawing.set_item(
                "style",
                settings.scene_style(py, labels.clone(), "node", &style)?,
            )?;
            drawing.set_item("minimum_size", 2.0 * node.radius)?;
            let details = py.import("json")?.call_method1(
                "loads",
                (serde_json::to_string(
                    &node
                        .details
                        .fields()
                        .iter()
                        .cloned()
                        .collect::<serde_json::Map<_, _>>(),
                )
                .map_err(|e| PyRuntimeError::new_err(e.to_string()))?,),
            )?;
            extensions.set_item("inspection", &details)?;
            nodes.append(drawing)?;
        }
        for edge in &scene.edges {
            let drawing = PyDict::new(py);
            label(&drawing, edge.label)?;
            let extensions = PyDict::new(py);
            extensions.set_item("scene-momentum", edge.momentum)?;
            extensions.set_item("scene-momentum-offset", edge.momentum_offset)?;
            drawing.set_item("extensions", &extensions)?;
            let style = PyDict::new(py);
            style.set_item("stroke", stroke(py, &edge.stroke)?)?;
            if let Some(pattern) = edge.pattern {
                let (name, amplitude, wavelength, scale) = match pattern {
                    Pattern::Coil {
                        amplitude,
                        wavelength,
                        longitudinal_scale,
                    } => ("coil", amplitude, wavelength, Some(longitudinal_scale)),
                    Pattern::Wave {
                        amplitude,
                        wavelength,
                    } => ("wave", amplitude, wavelength, None),
                    Pattern::Zigzag {
                        amplitude,
                        wavelength,
                    } => ("zigzag", amplitude, wavelength, None),
                };
                style.set_item("pattern", name)?;
                style.set_item("pattern-amplitude", amplitude)?;
                style.set_item("pattern-wavelength", wavelength)?;
                if let Some(scale) = scale {
                    style.set_item("pattern-coil-longitudinal-scale", scale)?;
                }
            }
            if let Some(forward) = edge.flow {
                let mark = PyDict::new(py);
                mark.set_item(if forward { "end" } else { "start" }, ">")?;
                style.set_item("mark", mark)?;
            }
            drawing.set_item(
                "style",
                settings.scene_style(py, labels.clone(), "edge", &style)?,
            )?;
            let details = py.import("json")?.call_method1(
                "loads",
                (serde_json::to_string(
                    &edge
                        .details
                        .fields()
                        .iter()
                        .cloned()
                        .collect::<serde_json::Map<_, _>>(),
                )
                .map_err(|e| PyRuntimeError::new_err(e.to_string()))?,),
            )?;
            extensions.set_item("inspection", &details)?;
            edges.append(drawing)?;
        }
        let elements = PyDict::new(py);
        elements.set_item("graph", PyDict::new(py))?;
        elements.set_item("nodes", &nodes)?;
        elements.set_item("edges", &edges)?;
        let hedges = PyList::empty(py);
        for spec in &scene.graph.edges {
            for _ in spec.source.iter().chain(&spec.sink) {
                hedges.append(PyDict::new(py))?;
            }
        }
        elements.set_item("hedges", &hedges)?;
        let transport = if settings.has_selectors(py) {
            let selector = SCENE_SELECTORS.get().ok_or_else(|| {
                PyRuntimeError::new_err(
                    "drawing selectors require the optional Linnet Python graph bindings",
                )
            })?;
            selector(py, &scene.graph, &elements, &settings)?
        } else {
            let normalized = PyDict::new(py);
            normalized.set_item("graph", PyDict::new(py))?;
            for (name, kind) in [
                ("nodes", crate::drawing::DrawingKind::Node),
                ("edges", crate::drawing::DrawingKind::Edge),
                ("hedges", crate::drawing::DrawingKind::HalfEdge),
            ] {
                let values = PyList::empty(py);
                for value in elements.get_item(name)?.unwrap().cast::<PyList>()? {
                    values.append(crate::drawing::element_drawing(
                        py,
                        value.cast::<PyDict>()?,
                        kind,
                    )?)?;
                }
                normalized.set_item(name, values)?;
            }
            let base = Py::new(py, settings)?.into_any();
            crate::typst::render_config_transport(py, &base, None, &normalized)?
        };
        // Preparation snapshots every referenced label file before the temporary module is dropped.
        Ok(Self::prepared(crate::render::prepare(
            linnest::encode_graph_spec_bytes(&scene.graph).map_err(PyRuntimeError::new_err)?,
            transport,
            None,
        )?))
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
