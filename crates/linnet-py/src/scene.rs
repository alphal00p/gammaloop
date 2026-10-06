//! Domain scenes use the native renderer when possible and the common Typst pipeline otherwise.
use crate::{
    dot::PyGlobalData,
    graph::{EdgeRecord, GraphState, HalfEdgeRecord, NodeRecord, PyGraph},
    native_graph::{PyHedgeGraph, PyNodeStore},
    PyDiagramRender, PyRenderSettings,
};
use linnest::svg::{Dash, Pattern, Scene, Stroke};
use linnet::half_edge::{
    builder::{HedgeData, HedgeGraphBuilder},
    involution::{Flow, Orientation},
};
use pyo3::{
    exceptions::PyRuntimeError,
    prelude::*,
    types::{PyDict, PyDictMethods},
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
        let mut builder = HedgeGraphBuilder::new();
        let mut nodes = BTreeMap::new();
        for (index, (spec, node)) in scene.graph.nodes.iter().zip(&scene.nodes).enumerate() {
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
            nodes.insert(
                spec.index.unwrap_or(index),
                builder.add_node(NodeRecord {
                    name: spec.name.clone(),
                    data: details.unbind(),
                    drawing: drawing.unbind(),
                }),
            );
        }
        for (spec, edge) in scene.graph.edges.iter().zip(&scene.edges) {
            let drawing = PyDict::new(py);
            label(&drawing, edge.label)?;
            let extensions = PyDict::new(py);
            extensions.set_item("scene-momentum", edge.momentum)?;
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
            let record = EdgeRecord {
                name: spec.name.clone(),
                data: details.unbind(),
                drawing: drawing.unbind(),
            };
            let endpoint = |end: &linnest::TypstEndpointSpec| HedgeData {
                node: nodes[&end.node],
                data: HalfEdgeRecord {
                    data: py.None(),
                    drawing: PyDict::new(py).unbind(),
                },
                is_in_subgraph: end.in_subgraph,
            };
            let orientation = match spec.orientation.as_deref() {
                Some("undirected") => Orientation::Undirected,
                Some("reversed") => Orientation::Reversed,
                _ => Orientation::Default,
            };
            match (&spec.source, &spec.sink) {
                (Some(a), Some(b)) => {
                    builder.add_edge(endpoint(a), endpoint(b), record, orientation)
                }
                (Some(a), None) => {
                    builder.add_external_edge(endpoint(a), record, orientation, Flow::Source);
                }
                (None, Some(b)) => {
                    builder.add_external_edge(endpoint(b), record, orientation, Flow::Sink);
                }
                _ => return Err(PyRuntimeError::new_err("scene edge has no endpoints")),
            }
        }
        let graph = Py::new(
            py,
            PyGraph::from_state(GraphState {
                graph: PyHedgeGraph::build(builder, PyNodeStore::Vec),
                name: scene.graph.name.clone(),
                global_data: PyGlobalData::default(),
                codec: None,
                render_config: Py::new(py, settings)?.into_any(),
                revision: 0,
            }),
        )?;
        // Preparation snapshots every referenced label file before the temporary module is dropped.
        Ok(Self::prepared(crate::render::prepare_scene(
            py,
            &graph,
            &scene.graph,
        )?))
    }
}
