//! Network graph rendering through Linnet's shared asset and SVG pipeline.
use pyo3::{
    prelude::*,
    types::{PyBytes, PyDict},
};

use crate::network::SpensoNet;

impl SpensoNet {
    pub(crate) fn prepare_render<'py>(
        &self,
        py: Python<'py>,
        config: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let linnet = py.import("linnet")?;
        let options = PyDict::new(py);
        options.set_item("network-dot", self.to_dot())?;
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
                    "#import \"crates/linnest/typst/src/render/network.typ\": render\n",
                    "#render(_linnet_config.options.at(\"network-dot\"), config: _linnet_config)\n",
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
