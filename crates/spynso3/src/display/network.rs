//! Network graph rendering through Linnet's shared asset and SVG pipeline.
use super::{DisplaySettings, IndexAliases, TensorDisplayMode, format_atom_with_mode};
use crate::{metadata::SpensoRepresentationName, network::SpensoNet};
use linnet::half_edge::involution::{Flow, HedgePair, Orientation};
use pyo3::{
    prelude::*,
    types::{PyBytes, PyDict, PyList},
};
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
                    NetworkEdge::Slot(slot) => slot.to_string(),
                },
            )?;
            record.set_item(
                "label-typst",
                match edge.data {
                    NetworkEdge::Head => None,
                    NetworkEdge::Slot(slot) => Some(aliases.port_typst(
                        slot.rep().slot(PartialIndex::Explicit(slot.aind())),
                        &settings,
                    )),
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

    pub(crate) fn prepare_render<'py>(
        &self,
        py: Python<'py>,
        config: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let linnet = py.import("linnet")?;
        let options = PyDict::new(py);
        options.set_item("network", self.render_snapshot(py)?)?;
        let kwargs = PyDict::new(py);
        kwargs.set_item("template_options", options)?;
        let mut effective = linnet.getattr("RenderConfig")?.call((), Some(&kwargs))?;
        if let Some(config) = config {
            effective = config.call_method1("overlay", (effective,))?;
        }
        let sources = PyDict::new(py);
        sources.set_item(
            "main.typ",
            PyBytes::new(
                py,
                concat!(
                    "#import \"crates/linnest/typst/src/render/network.typ\": render-network\n",
                    "#render-network(_linnet_config.options.at(\"network\"), config: _linnet_config)\n",
                )
                .as_bytes(),
            ),
        )?;
        let kwargs = PyDict::new(py);
        kwargs.set_item("config", effective)?;
        linnet
            .getattr("PreparedRender")?
            .call_method("from_sources", (sources,), Some(&kwargs))
    }
}

pub(crate) fn html(svg: &str) -> String {
    format!(
        "<figure class=\"spenso-network\" style=\"max-width:100%;margin:.5rem 0\">\
         <figcaption style=\"font:11px/1.45 ui-monospace,monospace;opacity:.7;margin-bottom:10px\">TensorNetwork</figcaption>\
         <div style=\"max-width:100%;overflow:auto\">{svg}</div></figure>"
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
