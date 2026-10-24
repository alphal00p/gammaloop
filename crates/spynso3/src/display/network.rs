//! Native SVG network drawing with direct MathML node and port labels.
use crate::display;
use crate::display::{DisplaySettings, IndexAliases, TensorDisplayMode, metadata};
use crate::metadata::SpensoRepresentationName;
use crate::network::SpensoNet;
use crate::network::execution::ExecutionStatus;
use linnest::svg::{Dash, Details, EdgeDrawing, LabelSource, NodeDrawing, Scene, Stroke};
use linnest::{TypstEdgeSpec, TypstEndpointSpec, TypstGraphSpec, TypstNodeSpec};
use linnet::half_edge::involution::{Flow, HedgePair, Orientation};
use linnet_render_py::{PyDiagramView, PyRenderSettings};
use pyo3::exceptions::PyRuntimeError;
use pyo3::types::{PyAnyMethods, PyDictMethods, PyListMethods};
use pyo3::types::{PyDict, PyList};
use pyo3::{Bound, PyErr, PyResult, Python};
use serde_json::{Value, json};
use spenso::network::graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp};
use spenso::network::parsing::ShadowedStructure;
use spenso::network::store::TensorScalarStore;
use spenso::network::tags::mathml;
use spenso::structure::HasName;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::PartialIndex;
use spenso::structure::slot::{IsAbstractSlot, ParseableAind};
use spenso::tensors::complex::RealOrComplexTensor;
use spenso::tensors::data::DataTensor;
use spenso::tensors::parametric::MixedTensor;
use std::cell::Cell;
use std::collections::BTreeMap;
use symbolica::atom::{Atom, FunctionBuilder, Symbol};
use symbolica::printer::PrintState;

/// Domain details feed Linnet's common hover and click inspector.
struct NetworkLabel {
    kind: &'static str,
    title: String,
    value: String,
    source: LabelSource,
    properties: Vec<(String, String)>,
}

impl NetworkLabel {
    fn new(kind: &'static str, title: impl Into<String>, value: impl Into<String>) -> Self {
        let value = value.into();
        Self {
            kind,
            title: title.into(),
            source: LabelSource::text(&value),
            value,
            properties: Vec::new(),
        }
    }

    fn property(&mut self, label: &str, value: impl ToString) {
        self.properties.push((label.into(), value.to_string()));
    }

    fn atom(&mut self, atom: &Atom, settings: &DisplaySettings, aliases: &IndexAliases) {
        self.value = display::format_atom_with_mode(atom, TensorDisplayMode::Plain, false);
        let options = settings.mathml_options();
        let invalid_label = Cell::new(false);
        let mathml = mathml::expression(
            atom.as_view(),
            &options,
            &PrintState::new(),
            &|representation, index| {
                let key = (representation, index.to_owned());
                if let Some((_, display)) = aliases.entries.get(&key) {
                    let rendered = display.to_mathml();
                    invalid_label.set(invalid_label.get() || rendered.is_none());
                    return rendered;
                }
                let source = aliases.source_indices.get(&key)?;
                let rendered =
                    mathml::expression(source.as_view(), &options, &PrintState::new(), &|_, _| {
                        None
                    });
                invalid_label.set(invalid_label.get() || rendered.is_err());
                rendered.ok()
            },
        );
        self.source = LabelSource {
            // Graph labels are bounded offers: a rejected export keeps literal text.
            mathml: mathml
                .ok()
                .filter(|_| !invalid_label.get())
                .unwrap_or_else(|| LabelSource::text(&self.value).mathml),
            typst: Some(format!(
                "${}$",
                display::format_atom_with_mode(atom, TensorDisplayMode::Typst, false)
            )),
        };
    }

    fn named(
        &mut self,
        named: &impl HasName<Name = Symbol, Args = Vec<Atom>>,
        settings: &DisplaySettings,
        aliases: &IndexAliases,
    ) {
        if let Some(name) = named.name() {
            let args = named.args().unwrap_or_default();
            let atom = if args.is_empty() {
                Atom::var(name)
            } else {
                FunctionBuilder::new(name).add_args(&args).finish()
            };
            self.atom(&atom, settings, aliases);
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
        settings: &DisplaySettings,
        aliases: &IndexAliases,
    ) {
        self.named(tensor, settings, aliases);
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
                        NetworkOp::Power(power) => {
                            // The parenthesized dot denotes the incoming expression,
                            // without inventing a variable name for that subnetwork.
                            label.value = format!("(·)^{power}");
                            label.source = LabelSource {
                                mathml: format!(
                                    "<msup><mrow><mo>(</mo><mo>·</mo><mo>)</mo></mrow><mn>{}</mn></msup>",
                                    display::escape_html(&power.to_string())
                                ),
                                typst: Some(format!("$(dot.c)^({power})$")),
                            };
                            label.property("Exponent", power);
                        }
                        NetworkOp::Function(name) => label.property("Function", name.get_name()),
                        _ => {}
                    }
                    label
                }
                NetworkNode::Leaf(leaf) => match leaf {
                    NetworkLeaf::LibraryKey { key, .. } => {
                        let mut label = NetworkLabel::new("library", "Library tensor", "Tensor");
                        label.named(key.canonical(), &settings, &aliases);
                        label.property("Resolution", "Library reference");
                        label
                    }
                    NetworkLeaf::LocalTensor(index) => {
                        let mut label = NetworkLabel::new("tensor", "Stored tensor", "Tensor");
                        label.tensor(store.get_tensor(*index), *index, &settings, &aliases);
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
                        label.tensor(
                            store.get_tensor(term.tensor),
                            term.tensor,
                            &settings,
                            &aliases,
                        );
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
                        label.atom(scalar, &settings, &aliases);
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
            // Preserve any pending Python exception at the snapshot boundary.
            // Registered Python printer failures are instead logged by Symbolica
            // and decline the override without leaving a pending exception.
            if let Some(error) = PyErr::take(py) {
                return Err(error);
            }
            let record = PyDict::new(py);
            record.set_item("id", id.0)?;
            record.set_item("kind", label.kind)?;
            record.set_item("value", &label.value)?;
            record.set_item("label-mathml", &label.source.mathml)?;
            record.set_item("label-typst", &label.source.typst)?;
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
            let port_source = match edge.data {
                NetworkEdge::Head => None,
                NetworkEdge::Slot(slot) | NetworkEdge::BoundPort { slot, .. } => {
                    let slot = slot.rep().slot(PartialIndex::Explicit(slot.aind()));
                    let (index, display) = aliases.port_display(slot);
                    let mut label = NetworkLabel::new("slot", "Tensor port", "");
                    label.atom(&index, &settings, &aliases);
                    if let Some(display) = display {
                        label.source.mathml = display.to_mathml().unwrap_or_else(|| {
                            LabelSource::text(&display.to_native_string()).mathml
                        });
                    }
                    label.source.typst = Some(format!("${}$", aliases.port_typst(slot, &settings)));
                    Some(label.source)
                }
            };
            if let Some(error) = PyErr::take(py) {
                return Err(error);
            }
            record.set_item(
                "label-mathml",
                port_source.as_ref().map(|source| &source.mathml),
            )?;
            record.set_item(
                "label-typst",
                port_source
                    .as_ref()
                    .and_then(|source| source.typst.as_ref()),
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
                    detail.atom(store.get_scalar_ref(*value), &settings, &aliases);
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
                    detail.property("Metric", metadata::metric(rep));
                    detail
                }
            };
            if let Some(h) = source {
                detail.property("Source half-edge", h.0);
            }
            if let Some(h) = sink {
                detail.property("Sink half-edge", h.0);
            }
            if let Some(error) = PyErr::take(py) {
                return Err(error);
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
        config: Option<&PyRenderSettings>,
    ) -> PyResult<PyDiagramView> {
        let json = py.import("json")?;
        let snapshot: Value = serde_json::from_str(
            &json
                .call_method1("dumps", (self.render_snapshot(py)?,))?
                .extract::<String>()?,
        )
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        let runtime = PyRuntimeError::new_err;
        let scene = Self::network_scene(&snapshot).map_err(runtime)?;
        PyDiagramView::from_scene(py, scene, config)
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
            layout: serde_json::from_value(json!({
                "layout-algo": "dot",
                "tree-dx": 0.35,
                "tree-dy": 2.2,
                "route-label-width-cap": 0.0,
                "label-steps": 40,
                "internal-label-length-scale": 0.35,
                "external-label-length-scale": 0.35,
            }))
            .expect("the built-in network layout settings are valid"),
            layout_passes: Vec::new(),
            layout_edges: Some(Vec::new()),
            label_feedback: false,
        };
        let mut indices = BTreeMap::new();
        for (index, record) in nodes.iter().enumerate() {
            let original = record["id"].as_u64().ok_or("invalid network node ID")? as usize;
            indices.insert(original, index);
            let page = scene.pages.len();
            scene.pages.push(LabelSource {
                mathml: record["label-mathml"]
                    .as_str()
                    .ok_or("missing node MathML label")?
                    .into(),
                typst: record["label-typst"].as_str().map(str::to_owned),
            });
            let mut details: Details = record["inspection"]
                .as_object()
                .ok_or("missing node inspection")?
                .iter()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect();
            details.insert("node", original);
            for field in ["label-mathml", "label-typst"] {
                if let Some(label) = record[field].as_str() {
                    details.insert(field, label);
                }
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
                label_padding: Some(0.25),
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
        scene.layout.insert("layout-roots".into(), json!([root]));
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
            for field in ["label-mathml", "label-typst"] {
                if let Some(label) = record[field].as_str() {
                    details.insert(field, label);
                }
            }
            details.insert(
                "source-hedge",
                record["source"].get(1).cloned().unwrap_or(Value::Null),
            );
            details.insert(
                "sink-hedge",
                record["sink"].get(1).cloned().unwrap_or(Value::Null),
            );
            let label = record["label-mathml"].as_str().map(|mathml| {
                let page = scene.pages.len();
                scene.pages.push(LabelSource {
                    mathml: mathml.into(),
                    typst: record["label-typst"].as_str().map(str::to_owned),
                });
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
                momentum_offset: EdgeDrawing::DEFAULT_MOMENTUM_OFFSET,
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
        ("#3d2645", "#c8b6dc"),
        ("#f5f5f5", "#2b3139"),
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

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::network::tags::SPENSO_TAG;
    use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
    use symbolica::atom::{AtomCore, NamespacedSymbol, SymbolBuilder};
    use symbolica::{symbol, wrap_symbol};

    #[test]
    fn network_labels_keep_literal_and_complete_expression_sources() {
        let text = "<&>\" #import \"evil.typ\"";
        let label = NetworkLabel::new("tensor", "Tensor", text);
        assert_eq!(label.source, LabelSource::text(text));
        let atom = Atom::num(2) + Atom::var(symbol!("network_mathml_tests::x"));
        let mut label = NetworkLabel::new("scalar", "Scalar", "");
        label.atom(&atom, &DisplaySettings::default(), &IndexAliases::default());
        assert!(
            label.source.mathml.contains("<mo>+</mo>"),
            "{}",
            label.source.mathml
        );
        assert!(label.source.mathml.contains("<mn>2</mn>"));
        assert!(label.source.typst.as_deref().unwrap().starts_with('$'));
        assert!(!label.source.mathml.contains("<math"));

        let large = Atom::var(
            SymbolBuilder::new(NamespacedSymbol::parse(&format!(
                "network_mathml_tests::{}",
                "x".repeat(5000)
            )))
            .build()
            .unwrap(),
        );
        label.atom(
            &large,
            &DisplaySettings::default(),
            &IndexAliases::default(),
        );
        assert_eq!(label.source.mathml, LabelSource::text(&label.value).mathml);
        assert!(label.source.typst.is_some());
    }

    #[test]
    fn network_aliases_are_shared_without_rewriting_source_atoms() {
        let head = SymbolBuilder::new(wrap_symbol!("network_mathml_tests::hedge"))
            .with_tags([SPENSO_TAG.index.clone()])
            .build()
            .unwrap();
        let index = |edge| {
            FunctionBuilder::new(head)
                .add_arg(Atom::num(edge))
                .add_arg(Atom::num(1))
                .finish()
        };
        let rep = LibraryRep::from(Minkowski {});
        let tensor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("network_mathml_tests::T"))
            .add_arg(rep.to_symbolic([Atom::num(4), index(0)]))
            .add_arg(rep.to_symbolic([Atom::num(4), index(1)]))
            .finish();
        let original = tensor.to_canonical_string();
        let aliases = IndexAliases::for_atom(&tensor, "alphabet");
        assert_eq!(aliases.entries.len(), 2);
        let mut label = NetworkLabel::new("tensor", "Tensor", "");
        label.atom(&tensor, &DisplaySettings::default(), &aliases);
        for (_, display) in aliases.entries.values() {
            assert!(
                label.source.mathml.contains(&display.to_mathml().unwrap()),
                "{}",
                label.source.mathml
            );
        }
        assert!(!label.source.mathml.contains("hedge"));
        assert_eq!(tensor.to_canonical_string(), original);

        // A decoded index which cannot be exported must decline the entire
        // label offer, rather than silently displaying its encoded identity.
        let mut aliases = aliases;
        let key = aliases.entries.keys().next().unwrap().clone();
        aliases.entries.remove(&key);
        let source = Atom::var(
            SymbolBuilder::new(NamespacedSymbol::parse(&format!(
                "network_mathml_tests::{}",
                "s".repeat(5000)
            )))
            .build()
            .unwrap(),
        );
        aliases.source_indices.insert(key, source);
        label.atom(&tensor, &DisplaySettings::default(), &aliases);
        assert_eq!(label.source.mathml, LabelSource::text(&label.value).mathml);
        assert_eq!(tensor.to_canonical_string(), original);
    }

    #[test]
    fn mathml_network_scene_preserves_topology_styles_and_inspection() {
        let source = LabelSource::text("A <&>");
        let snapshot = json!({
            "root": 19,
            "nodes": [
                {"id": 19, "kind": "operator", "label-mathml": "<mo>∏</mo>", "label-typst": "$product$", "inspection": {"title": "Product", "summary": "2 inputs"}},
                {"id": 7, "kind": "tensor", "label-mathml": source.mathml, "label-typst": source.typst, "inspection": {"title": "Tensor", "properties": [["Rank", "1"]]}},
            ],
            "edges": [
                {"id": 23, "source": [7, 31], "sink": [19, 47], "tree": true, "orientation": "default", "label-mathml": null, "label-typst": null, "inspection": {"title": "Dependency"}},
                {"id": 41, "source": [7, 53], "sink": null, "tree": false, "orientation": "reversed", "label-mathml": "<mi>μ</mi>", "label-typst": "$mu$", "inspection": {"title": "Free tensor slot"}},
            ],
        });
        let before = snapshot.clone();
        let scene = SpensoNet::network_scene(&snapshot).unwrap();
        assert_eq!(snapshot, before);
        assert_eq!(
            scene
                .graph
                .nodes
                .iter()
                .map(|node| node.name.as_deref())
                .collect::<Vec<_>>(),
            vec![Some("n19"), Some("n7")]
        );
        let edge = &scene.graph.edges[0];
        assert_eq!(
            (
                edge.source.as_ref().unwrap().node,
                edge.source.as_ref().unwrap().id
            ),
            (1, Some(31))
        );
        assert_eq!(
            (
                edge.sink.as_ref().unwrap().node,
                edge.sink.as_ref().unwrap().id
            ),
            (0, Some(47))
        );
        assert_eq!(scene.graph.edges[1].source.as_ref().unwrap().id, Some(53));
        assert!(scene.graph.edges[1].sink.is_none());
        assert_eq!(scene.layout_edges, Some(vec![0]));
        assert_eq!(scene.layout["layout-roots"], json!([0]));
        assert_eq!(scene.pages[1], source);
        assert_eq!(scene.pages[2].mathml, "<mi>μ</mi>");
        assert_eq!(scene.pages[2].typst.as_deref(), Some("$mu$"));
        assert_eq!(scene.nodes[0].fill, "#f5f5f5");
        assert_eq!(scene.nodes[1].fill, "#ffffff");
        assert_eq!(scene.edges[0].stroke.paint, "#555555");
        assert_eq!(scene.edges[1].stroke.paint, "#7f95b8");
        assert_eq!(scene.nodes[0].details.get("title"), Some(&json!("Product")));
        assert_eq!(
            scene.nodes[1].details.get("properties"),
            Some(&json!([["Rank", "1"]]))
        );
        assert_eq!(scene.edges[1].details.get("edge"), Some(&json!(41)));
        assert_eq!(scene.edges[1].details.get("source-hedge"), Some(&json!(53)));
    }
}
