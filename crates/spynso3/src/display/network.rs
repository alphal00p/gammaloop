//! Direct SVG network drawing; Typst typesets only node and port labels.
use super::{DisplaySettings, IndexAliases, TensorDisplayMode, format_atom_with_mode};
use crate::{
    metadata::SpensoRepresentationName,
    network::{SpensoNet, execution::ExecutionStatus},
};
use linnest::{
    TypstEdgeSpec, TypstEndpointSpec, TypstGraphSpec, TypstNodeSpec,
    svg::{Config, Dash, Details, EdgeDrawing, NodeDrawing, Scene, Stroke},
};
use linnet::half_edge::involution::{Flow, HedgePair, Orientation};
use pyo3::{
    prelude::*,
    types::{PyDict, PyList},
};
use serde_json::Value;
use spenso::{
    network::{
        graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp},
        parsing::ShadowedStructure,
        store::TensorScalarStore,
    },
    structure::{
        HasName,
        abstract_index::AbstractIndex,
        partial::PartialIndex,
        slot::{IsAbstractSlot, ParseableAind},
    },
    tensors::{complex::RealOrComplexTensor, data::DataTensor, parametric::MixedTensor},
};
use std::collections::BTreeMap;
use symbolica::atom::{Atom, FunctionBuilder, Symbol};

/// Domain details feed Linnet's common hover and click inspector.
struct NetworkLabel {
    kind: &'static str,
    title: String,
    value: String,
    typst: Option<String>,
    properties: Vec<(String, String)>,
}

impl NetworkLabel {
    fn new(kind: &'static str, title: impl Into<String>, value: impl Into<String>) -> Self {
        Self {
            kind,
            title: title.into(),
            value: value.into(),
            typst: None,
            properties: Vec::new(),
        }
    }

    fn property(&mut self, label: &str, value: impl ToString) {
        self.properties.push((label.into(), value.to_string()));
    }

    fn atom(&mut self, atom: &Atom) {
        self.value = format_atom_with_mode(atom, TensorDisplayMode::Plain, false);
        self.typst = Some(format_atom_with_mode(atom, TensorDisplayMode::Typst, false));
    }

    fn named(&mut self, named: &impl HasName<Name = Symbol, Args = Vec<Atom>>) {
        if let Some(name) = named.name() {
            let args = named.args().unwrap_or_default();
            let atom = if args.is_empty() {
                Atom::var(name)
            } else {
                FunctionBuilder::new(name).add_args(&args).finish()
            };
            self.atom(&atom);
            self.property("Tensor name", name.get_name());
            if !args.is_empty() {
                self.property(
                    "Arguments",
                    args.iter()
                        .map(Atom::to_plain_string)
                        .collect::<Vec<_>>()
                        .join(", "),
                );
            }
        }
    }

    fn storage<T, S>(&mut self, data: &DataTensor<T, S>) {
        let (storage, count) = match data {
            DataTensor::Dense(data) => ("Dense", data.data.len()),
            DataTensor::Sparse(data) => ("Sparse", data.elements.len()),
        };
        self.property("Storage", storage);
        self.property("Stored components", count);
    }

    fn tensor(
        &mut self,
        tensor: &MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
        index: usize,
    ) {
        self.named(tensor);
        self.property("Store index", index);
        match tensor {
            MixedTensor::Param(data) => {
                self.property("Components", "Symbolic");
                self.storage(&data.tensor);
            }
            MixedTensor::Concrete(RealOrComplexTensor::Real(data)) => {
                self.property("Components", "Real");
                self.storage(data);
            }
            MixedTensor::Concrete(RealOrComplexTensor::Complex(data)) => {
                self.property("Components", "Complex");
                self.storage(data);
            }
        }
    }

    fn inspection<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let details = PyDict::new(py);
        details.set_item("title", &self.title)?;
        let mut summary = format!("{} · {}", self.title, self.value);
        for (label, value) in &self.properties {
            if matches!(
                label.as_str(),
                "Rank"
                    | "Port dimensions"
                    | "Components"
                    | "Inputs"
                    | "Representation"
                    | "Dimension"
                    | "Index"
            ) {
                summary.push_str(&format!("\n{label}: {value}"));
            }
        }
        details.set_item("summary", summary)?;
        details.set_item("properties", &self.properties)?;
        Ok(details)
    }
}

impl SpensoNet {
    /// Copy native topology and leaf metadata into Linnet's typed configuration.
    /// No textual graph format is involved; original node/edge/half-edge IDs survive.
    fn render_snapshot<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let graph = &self.network.graph.graph;
        let store = &self.network.store;
        let settings = DisplaySettings::default();
        // Resolve compound indices together so different edges cannot receive
        // the same alphabet alias just because each slot was printed alone.
        let mut ports = FunctionBuilder::new(spenso::structure::abstract_index::AIND_SYMBOLS.aind);
        for (_, _, edge) in graph.iter_edges() {
            if let NetworkEdge::Slot(slot) = edge.data {
                ports = ports.add_arg(slot.to_atom());
            }
        }
        let aliases = IndexAliases::for_atom(&ports.finish(), &settings.index_style);
        let nodes = PyList::empty(py);
        for (id, _, node) in graph.iter_nodes() {
            let mut label = match node {
                NetworkNode::Op(op) => {
                    let title = match op {
                        NetworkOp::Product => "Product",
                        NetworkOp::Sum => "Sum",
                        NetworkOp::Neg => "Negation",
                        NetworkOp::Power(_) => "Power",
                        NetworkOp::Function(_) => "Function",
                    };
                    let mut label =
                        NetworkLabel::new("operator", title, op.display_with(ToString::to_string));
                    let inputs = graph
                        .iter_crown(id)
                        .filter(|&h| {
                            graph.get_edge_data(h).is_head() && graph.flow(h) == Flow::Sink
                        })
                        .count();
                    label.property("Inputs", inputs);
                    match op {
                        NetworkOp::Power(power) => label.property("Exponent", power),
                        NetworkOp::Function(name) => label.property("Function", name.get_name()),
                        _ => {}
                    }
                    label
                }
                NetworkNode::Leaf(leaf) => match leaf {
                    NetworkLeaf::LibraryKey { key, .. } => {
                        let mut label = NetworkLabel::new("library", "Library tensor", "Tensor");
                        label.named(key.canonical());
                        label.property("Resolution", "Library reference");
                        label
                    }
                    NetworkLeaf::LocalTensor(index) => {
                        let mut label = NetworkLabel::new("tensor", "Stored tensor", "Tensor");
                        label.tensor(store.get_tensor(*index), *index);
                        label
                    }
                    NetworkLeaf::TensorSum(indices) => {
                        let mut label = NetworkLabel::new(
                            "tensor",
                            "Tensor sum",
                            format!("Sum of {} tensors", indices.len()),
                        );
                        label.property("Terms", indices.len());
                        label.property("Store indices", format!("{indices:?}"));
                        label
                    }
                    NetworkLeaf::ScaledTensor(term) => {
                        let mut label = NetworkLabel::new("tensor", "Scaled tensor", "Tensor");
                        label.tensor(store.get_tensor(term.tensor), term.tensor);
                        if let Some(scale) = term.scale {
                            label.property(
                                "Coefficient",
                                store.get_scalar_ref(scale).to_plain_string(),
                            );
                        }
                        label
                    }
                    NetworkLeaf::ScaledTensorSum(terms) => {
                        let mut label = NetworkLabel::new(
                            "tensor",
                            "Scaled tensor sum",
                            format!("Sum of {} terms", terms.len()),
                        );
                        label.property("Terms", terms.len());
                        label.property(
                            "Store indices",
                            terms
                                .iter()
                                .map(|term| term.tensor.to_string())
                                .collect::<Vec<_>>()
                                .join(", "),
                        );
                        label
                    }
                    NetworkLeaf::Scalar(scalar) => {
                        let mut label = NetworkLabel::new("scalar", "Scalar", "");
                        let scalar = store.get_scalar_ref(*scalar);
                        label.atom(scalar);
                        label.property("Expression", scalar.to_plain_string());
                        label
                    }
                },
            };
            let slots = graph
                .iter_crown(id)
                .filter_map(|h| match graph.get_edge_data(h) {
                    NetworkEdge::Slot(slot) => Some(*slot),
                    _ => None,
                })
                .collect::<Vec<_>>();
            if !slots.is_empty() || label.kind == "tensor" || label.kind == "library" {
                label.property("Rank", slots.len());
                label.property(
                    "Port dimensions",
                    slots
                        .iter()
                        .map(|slot| slot.rep().dim.to_string())
                        .collect::<Vec<_>>()
                        .join(" × "),
                );
                label.property(
                    "Slots",
                    slots
                        .iter()
                        .map(|slot| slot.to_atom().to_plain_string())
                        .collect::<Vec<_>>()
                        .join(", "),
                );
            }
            let record = PyDict::new(py);
            record.set_item("id", id.0)?;
            record.set_item("kind", label.kind)?;
            record.set_item("value", &label.value)?;
            record.set_item("label-typst", &label.typst)?;
            record.set_item("inspection", label.inspection(py)?)?;
            nodes.append(record)?;
        }
        let edges = PyList::empty(py);
        for (pair, id, edge) in graph.iter_edges() {
            let (source, sink) = match pair {
                HedgePair::Paired { source, sink } => (Some(source), Some(sink)),
                HedgePair::Unpaired {
                    hedge,
                    flow: Flow::Source,
                } => (Some(hedge), None),
                HedgePair::Unpaired {
                    hedge,
                    flow: Flow::Sink,
                } => (None, Some(hedge)),
                HedgePair::Split { .. } => unreachable!("whole network graphs have no split edges"),
            };
            let record = PyDict::new(py);
            record.set_item("id", id.0)?;
            record.set_item("source", source.map(|h| (graph.node_id(h).0, h.0)))?;
            record.set_item("sink", sink.map(|h| (graph.node_id(h).0, h.0)))?;
            record.set_item(
                "orientation",
                match edge.orientation {
                    Orientation::Default => "default",
                    Orientation::Reversed => "reversed",
                    Orientation::Undirected => "undirected",
                },
            )?;
            record.set_item("tree", edge.data.is_head())?;
            record.set_item(
                "slot",
                match edge.data {
                    NetworkEdge::Head => String::new(),
                    NetworkEdge::Slot(slot) | NetworkEdge::BoundPort { slot, .. } => {
                        slot.to_string()
                    }
                },
            )?;
            record.set_item(
                "label-typst",
                match edge.data {
                    NetworkEdge::Head => None,
                    NetworkEdge::Slot(slot) | NetworkEdge::BoundPort { slot, .. } => {
                        Some(aliases.port_typst(
                            slot.rep().slot(PartialIndex::Explicit(slot.aind())),
                            &settings,
                        ))
                    }
                },
            )?;
            let mut detail = match edge.data {
                NetworkEdge::Head => NetworkLabel::new(
                    "dependency",
                    if source.is_none() || sink.is_none() {
                        "Network output"
                    } else {
                        "Expression dependency"
                    },
                    "Carries an expression result",
                ),
                NetworkEdge::BoundPort { slot, value } => {
                    let mut detail =
                        NetworkLabel::new("binding", "Bound tensor port", slot.to_string());
                    detail.atom(store.get_scalar_ref(*value));
                    detail.property("Slot", slot.to_string());
                    detail
                }
                NetworkEdge::Slot(slot) => {
                    let rep = slot.rep();
                    let mut detail = NetworkLabel::new(
                        "slot",
                        if source.is_none() || sink.is_none() {
                            "Free tensor slot"
                        } else {
                            "Tensor contraction"
                        },
                        slot.to_string(),
                    );
                    let name = SpensoRepresentationName { rep: rep.rep };
                    detail.property("Representation", rep.rep.symbol().get_name());
                    detail.property("Dimension", rep.dim);
                    detail.property("Index", slot.aind().to_atom().to_plain_string());
                    detail.property("Duality", name.duality());
                    detail.property("Metric", super::metadata::metric(rep));
                    detail
                }
            };
            if let Some(h) = source {
                detail.property("Source half-edge", h.0);
            }
            if let Some(h) = sink {
                detail.property("Sink half-edge", h.0);
            }
            record.set_item("inspection", detail.inspection(py)?)?;
            edges.append(record)?;
        }
        let snapshot = PyDict::new(py);
        snapshot.set_item("nodes", nodes)?;
        snapshot.set_item("edges", edges)?;
        snapshot.set_item("root", graph.node_id(self.network.graph.head()).0)?;
        Ok(snapshot)
    }

    pub(crate) fn render_graph(
        &self,
        py: Python<'_>,
        config: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<String> {
        let json = py.import("json")?;
        let options = match config {
            Some(config) => json.call_method1("dumps", (config,))?.extract::<String>()?,
            None => "{}".into(),
        };
        let config =
            Config::from_json(&options).map_err(pyo3::exceptions::PyValueError::new_err)?;
        if !config.template_options.is_empty() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "tensor networks have no physics template options",
            ));
        }
        let snapshot: Value = serde_json::from_str(
            &json
                .call_method1("dumps", (self.render_snapshot(py)?,))?
                .extract::<String>()?,
        )
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        let runtime = pyo3::exceptions::PyRuntimeError::new_err;
        let mut scene = Self::network_scene(&snapshot).map_err(runtime)?;
        config
            .apply(&mut scene)
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        let files = BTreeMap::from([("main.typ".into(), scene.label_document().into_bytes())]);
        let pages = typst_renderer::Document::compile_sources(&files, "svg")
            .map_err(runtime)?
            .into_iter()
            .map(String::from_utf8)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| runtime(e.to_string()))?;
        let typeset = scene.typeset(&pages).map_err(runtime)?;
        let svg = scene.render(&typeset).map_err(runtime)?;
        Scene::interactive_svg(&svg).map_err(runtime)
    }

    fn network_scene(snapshot: &Value) -> Result<Scene, String> {
        let nodes = snapshot["nodes"]
            .as_array()
            .ok_or("missing network nodes")?;
        let edges = snapshot["edges"]
            .as_array()
            .ok_or("missing network edges")?;
        let mut scene = Scene {
            graph: TypstGraphSpec {
                name: None,
                data: None,
                statements: BTreeMap::new(),
                default_edge_statements: BTreeMap::new(),
                default_node_statements: BTreeMap::new(),
                nodes: vec![],
                edges: vec![],
            },
            nodes: vec![],
            edges: vec![],
            preamble: "#set text(size: 9pt, fill: rgb(\"#000000\"))".into(),
            title: None,
            pages: vec![],
            layout: serde_json::from_value(serde_json::json!({
                "layout-algo": "dot",
                "tree-dx": 0.35,
                "tree-dy": 2.2,
                "route-label-width-cap": 0.0,
                "label-steps": 40,
                "internal-label-length-scale": 0.35,
                "external-label-length-scale": 0.35,
            }))
            .unwrap(),
            layout_edges: Some(Vec::new()),
            label_feedback: false,
        };
        let mut indices = BTreeMap::new();
        for (index, record) in nodes.iter().enumerate() {
            let original = record["id"].as_u64().ok_or("invalid network node ID")? as usize;
            indices.insert(original, index);
            let page = scene.pages.len();
            scene.pages.push(match record["label-typst"].as_str() {
                Some(math) => format!("[$ {math} $]"),
                None => format!("[#{:?}]", record["value"].as_str().unwrap_or("")),
            });
            let mut details: Details = record["inspection"]
                .as_object()
                .ok_or("missing node inspection")?
                .iter()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect();
            details.insert("node", original);
            if let Some(label) = record["label-typst"].as_str() {
                details.insert("label-typst", label);
            }
            let operator = record["kind"] == "operator";
            scene.graph.nodes.push(TypstNodeSpec {
                name: Some(format!("n{original}")),
                index: Some(index),
                data: None,
                pos: None,
                statements: BTreeMap::new(),
            });
            scene.nodes.push(NodeDrawing {
                radius: if operator { 0.6 } else { 0.35 },
                label: Some(page),
                rectangular: !operator,
                fill: if operator { "#f5f5f5" } else { "#ffffff" }.into(),
                stroke: Stroke {
                    paint: if operator { "#666666" } else { "#aeb4bd" }.into(),
                    width: 0.3,
                    dash: Dash::Solid,
                    round_cap: false,
                },
                details,
            });
        }
        let root = snapshot["root"].as_u64().ok_or("missing network root")? as usize;
        let root = *indices.get(&root).ok_or("unknown network root")?;
        scene
            .layout
            .insert("layout-roots".into(), serde_json::json!([root]));
        for (index, record) in edges.iter().enumerate() {
            let tree = record["tree"] == true;
            if tree {
                scene.layout_edges.as_mut().unwrap().push(index);
            }
            let endpoint =
                |value: &Value, source: bool| -> Result<Option<TypstEndpointSpec>, String> {
                    if value.is_null() {
                        return Ok(None);
                    }
                    let node = value[0].as_u64().ok_or("invalid network endpoint")? as usize;
                    let node = *indices.get(&node).ok_or("unknown network node")?;
                    Ok(Some(TypstEndpointSpec {
                        node,
                        id: value[1].as_u64().map(|h| h as usize),
                        statement: None,
                        data: None,
                        port_label: None,
                        compass: if tree {
                            Some(if source { "n" } else { "s" }.into())
                        } else {
                            scene.nodes[node].rectangular.then(|| "s".into())
                        },
                        in_subgraph: false,
                        route_points: vec![],
                    }))
                };
            let source = endpoint(&record["source"], true)?;
            let sink = endpoint(&record["sink"], false)?;
            let mut details: Details = record["inspection"]
                .as_object()
                .ok_or("missing edge inspection")?
                .iter()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect();
            details.insert("edge", record["id"].clone());
            for endpoint in ["source", "sink"] {
                details.insert(
                    endpoint,
                    record[endpoint].get(0).cloned().unwrap_or(Value::Null),
                );
            }
            if let Some(label) = record["label-typst"].as_str() {
                details.insert("label-typst", label);
            }
            details.insert(
                "source-hedge",
                record["source"].get(1).cloned().unwrap_or(Value::Null),
            );
            details.insert(
                "sink-hedge",
                record["sink"].get(1).cloned().unwrap_or(Value::Null),
            );
            let label = record["label-typst"].as_str().map(|math| {
                let page = scene.pages.len();
                scene.pages.push(format!("[$ {math} $]"));
                page
            });
            let orientation = record["orientation"].as_str().unwrap_or("default");
            scene.graph.edges.push(TypstEdgeSpec {
                name: None,
                source,
                sink,
                data: None,
                orientation: Some(orientation.into()),
                flow: None,
                id: Some(index),
                pos: None,
                statements: BTreeMap::from([(
                    "route".into(),
                    if tree { "direct" } else { "hobby-through" }.into(),
                )]),
            });
            scene.edges.push(EdgeDrawing {
                stroke: Stroke {
                    paint: if tree { "#555555" } else { "#7f95b8" }.into(),
                    width: if tree { 0.65 } else { 1.0 },
                    dash: Dash::Solid,
                    round_cap: true,
                },
                pattern: None,
                flow: tree.then_some(orientation != "reversed"),
                momentum: false,
                label,
                details,
            });
        }
        Ok(scene)
    }
}

pub(crate) fn html(svg: &str, status: &ExecutionStatus) -> String {
    let state = if status.complete {
        "Graph reduced"
    } else {
        "Pending"
    };
    let counts = [
        (status.nodes, "node"),
        (status.operations, "operation"),
        (status.contractions, "contraction"),
    ]
    .map(|(count, noun)| format!("{count} {noun}{}", if count == 1 { "" } else { "s" }))
    .join(" · ");
    let ready = if status.complete {
        String::new()
    } else {
        format!(" · Ready operations: {}", status.ready.len())
    };
    format!(
        "<figure class=\"spenso-network\" style=\"max-width:100%;margin:.5rem 0\">\
         <figcaption style=\"font:11px/1.45 ui-monospace,monospace;opacity:.7;margin-bottom:10px\">TensorNetwork</figcaption>\
         <div class=\"spenso-network-status\" data-complete=\"{}\" role=\"status\" \
         style=\"display:flex;flex-wrap:wrap;align-items:baseline;gap:6px 12px;font:12px/1.5 system-ui,sans-serif;margin-bottom:6px\" \
         title=\"Current graph snapshot; counts are not time or cost estimates. Graph reduced means no operations or contractions remain; library and deferred results may still need materialization.\">\
         <strong style=\"font-weight:600\">{state}</strong><span>{counts}{ready}</span></div>\
         <div style=\"max-width:100%;overflow:auto\">{svg}</div></figure>",
        status.complete
    )
}

pub(crate) fn svg_theme(svg: &str) -> String {
    // Keep the example's light palette as the exported SVG fallback. Linnet's
    // existing theme bridge sets data-theme on interactive notebook SVGs.
    let mut css = String::from(
        ".spenso-network-svg{color-scheme:light dark;max-width:100%;height:auto}\
         .spenso-network-svg[data-theme=light]{color-scheme:light}\
         .spenso-network-svg[data-theme=dark]{color-scheme:dark}",
    );
    for (light, dark) in [
        ("#000000", "#e6ebf1"),
        ("#ffffff", "#1c2025"),
        ("#ffffff80", "#1c202580"),
        ("#666666", "#a6b3c5"),
        ("#555555", "#b4bdc9"),
        ("#aeb4bd", "#566477"),
        ("#7f95b8", "#a4cdf7"),
        ("#d8dce3", "#596578"),
        ("#737985", "#b2bdce"),
        ("#fde8e8", "#482e35"),
        ("#e7f6e9", "#263e30"),
        ("#e7f0ff", "#293d55"),
    ] {
        for paint in ["fill", "stroke"] {
            css.push_str(&format!(
                ".spenso-network-svg [{paint}=\"{light}\"]{{{paint}:light-dark({light},{dark})}}"
            ));
        }
    }
    svg.replacen("<svg ", "<svg class=\"spenso-network-svg\" ", 1)
        .replacen('>', &format!("><style>{css}</style>"), 1)
}
