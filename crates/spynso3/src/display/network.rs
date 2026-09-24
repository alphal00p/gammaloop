//! Network graph rendering through Linnet's shared asset and SVG pipeline.
use linnet::half_edge::involution::{Flow, HedgePair, Orientation};
use pyo3::{
    prelude::*,
    types::{PyBytes, PyDict, PyList},
};
use spenso::network::{
    StructureLessDisplay,
    graph::{NetworkEdge, NetworkLeaf, NetworkNode},
    store::TensorScalarStore,
};

use crate::network::SpensoNet;

impl SpensoNet {
    /// Copy native topology and leaf metadata into Linnet's typed configuration.
    /// No textual graph format is involved; original node/edge/half-edge IDs survive.
    fn render_snapshot<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let graph = &self.network.graph.graph;
        let store = &self.network.store;
        let nodes = PyList::empty(py);
        for (id, _, node) in graph.iter_nodes() {
            let (kind, value) = match node {
                NetworkNode::Op(op) => ("operator", op.display_with(ToString::to_string)),
                NetworkNode::Leaf(leaf) => match leaf {
                    NetworkLeaf::LibraryKey { key, .. } => ("library", key.canonical().display()),
                    NetworkLeaf::LocalTensor(index) => {
                        ("tensor", store.get_tensor(*index).display())
                    }
                    NetworkLeaf::TensorSum(indices) => {
                        ("tensor", format!("Sum of {} tensors", indices.len()))
                    }
                    NetworkLeaf::ScaledTensor(term) => (
                        "tensor",
                        match term.scale {
                            Some(scale) => format!(
                                "{} · {}",
                                store.get_scalar_ref(scale),
                                store.get_tensor(term.tensor).display()
                            ),
                            None => store.get_tensor(term.tensor).display(),
                        },
                    ),
                    NetworkLeaf::ScaledTensorSum(terms) => {
                        ("tensor", format!("Sum of {} scaled tensors", terms.len()))
                    }
                    NetworkLeaf::Scalar(scalar) => {
                        ("scalar", store.get_scalar_ref(*scalar).to_string())
                    }
                },
            };
            let record = PyDict::new(py);
            record.set_item("id", id.0)?;
            record.set_item("kind", kind)?;
            record.set_item("value", value)?;
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
