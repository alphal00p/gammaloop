use feynkit_graph::FeynmanDiagram;
use feynkit_model::{Model, ParticleId};
use pyo3::{
    prelude::*,
    types::{PyBytes, PyDict},
};
use std::{collections::BTreeSet, fmt::Write};

pub(crate) fn escape_html(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            '\'' => escaped.push_str("&#39;"),
            _ => escaped.push(character),
        }
    }
    escaped
}

/// Compile a diagram's Linnest source with typst-py and return its SVG page.
pub(crate) fn render_diagram_svg(py: Python<'_>, prepared: &Bound<'_, PyAny>) -> PyResult<String> {
    py.import("typst").map_err(|error| {
        if error.is_instance_of::<pyo3::exceptions::PyImportError>(py) {
            pyo3::exceptions::PyImportError::new_err(format!(
                "diagram rendering requires typst-py: {error}"
            ))
        } else {
            error
        }
    })?;
    let svg: String = prepared.call_method0("to_svg")?.extract()?;
    // Use the website SVG palette from docs/assets/typst/theme.typ, including
    // the 45% lightened sink strokes from the shared physics style. Keep SVG
    // paint attributes as light-mode fallbacks for viewers without CSS support.
    let mut styles = String::from("<style>");
    let mut light = String::new();
    let mut dark = String::new();
    for (name, light_color, dark_color) in [
        ("ink", "#3d2645", "#f8effa"),
        ("accent", "#6f4d85", "#d8b9e3"),
        ("ink-sink", "#948899", "#fbf6fc"),
        ("accent-sink", "#b09dbc", "#ead9f0"),
        // Selection strokes stay visible over dark notebook backgrounds.
        ("highlight", "#ffd166", "#f0bb4f"),
        ("outside", "#77777773", "#aaaaaa73"),
    ] {
        light.push_str(&format!("--feynkit-{name}:{light_color};"));
        dark.push_str(&format!("--feynkit-{name}:{dark_color};"));
        for paint in ["fill", "stroke"] {
            styles.push_str(&format!(
                ".feynkit-diagram-svg [{paint}=\"{light_color}\"]{{{paint}:var(--feynkit-{name},{light_color});}}"
            ));
        }
    }
    // Explicit notebook themes override the browser preference for saved SVGs.
    // The interactive renderer mirrors an accessible iframe parent's theme onto
    // the SVG itself, so the same palette also applies inside notebook frames.
    styles.push_str(&format!(
        r#".feynkit-diagram-svg{{{light}}}
@media(prefers-color-scheme:dark){{.feynkit-diagram-svg{{{dark}}}}}
.feynkit-diagram-svg[data-theme="dark"],:is(.dark,[data-theme="dark"],[data-jp-theme-light="false"]) .feynkit-diagram-svg{{{dark}}}
.feynkit-diagram-svg[data-theme="light"],:is(.light,[data-theme="light"],[data-jp-theme-light="true"]) .feynkit-diagram-svg{{{light}}}
</style>"#
    ));
    Ok(svg
        .replacen("<svg ", "<svg class=\"feynkit-diagram-svg\" ", 1)
        .replacen('>', &format!(">{styles}"), 1))
}

pub(crate) fn render_diagram_html(diagram: &FeynmanDiagram, svg: &str) -> String {
    format!(
        "<figure class=\"feynkit-diagram\" style=\"max-width:100%;margin:.5rem 0\">\
         <div style=\"max-width:100%;overflow-x:auto\">{svg}</div>\
         <figcaption style=\"font-size:.85em;opacity:.75;text-align:center\">Feynman diagram {} \
         ({} loop{})</figcaption></figure>",
        escape_html(diagram.name()),
        diagram.loop_count(),
        if diagram.loop_count() == 1 { "" } else { "s" },
    )
}

#[cfg(test)]
mod tests {
    use super::escape_html;

    #[test]
    fn escapes_html_metadata() {
        assert_eq!(
            escape_html("<script data-x='a&b'>\"x\"</script>"),
            "&lt;script data-x=&#39;a&amp;b&#39;&gt;&quot;x&quot;&lt;/script&gt;"
        );
    }
}

/// Package the shared physics renderer once for diagram and process displays.
pub(crate) fn prepare_physics_render<'py>(
    py: Python<'py>,
    source: &str,
    config: Option<&Bound<'_, PyAny>>,
) -> PyResult<Bound<'py, PyAny>> {
    let sources = PyDict::new(py);
    sources.set_item("main.typ", PyBytes::new(py, source.as_bytes()))?;
    for (path, source) in [
        (
            "assets/embedded/drawing/templates/layout-core.typ",
            include_bytes!("../../../assets/embedded/drawing/templates/layout-core.typ").as_slice(),
        ),
        (
            "assets/embedded/drawing/templates/physics-edge-style.typ",
            include_bytes!("../../../assets/embedded/drawing/templates/physics-edge-style.typ")
                .as_slice(),
        ),
        (
            "assets/embedded/drawing/templates/impl/physics-edge-style.typ",
            include_bytes!(
                "../../../assets/embedded/drawing/templates/impl/physics-edge-style.typ"
            )
            .as_slice(),
        ),
    ] {
        sources.set_item(path, PyBytes::new(py, source))?;
    }
    let kwargs = PyDict::new(py);
    kwargs.set_item("config", config)?;
    py.import("linnet")?.getattr("PreparedRender")?.call_method(
        "from_sources",
        (sources,),
        Some(&kwargs),
    )
}

/// A process schematic uses the same particle labels, arrows and lines as diagrams.
/// It is a native Linnest star graph, not a partially initialized FeynmanDiagram.
pub(crate) fn process_svg(
    py: Python<'_>,
    model: &Model,
    incoming: &[ParticleId],
    outgoing: &[Vec<ParticleId>],
    config: Option<&Bound<'_, PyAny>>,
) -> PyResult<String> {
    let mut source = String::from(
        r##"#set page(width: auto, height: auto, margin: 2mm, fill: none)
#set text(size: 10pt)
#import "crates/linnest/typst/src/graph.typ" as graph
#import "crates/linnest/typst/src/render/layout.typ" as renderer
#import "assets/embedded/drawing/templates/layout-core.typ" as physics-layout
#import "assets/embedded/drawing/templates/physics-edge-style.typ" as physics
#import physics: mi, palette, massive, massless, dashed, dotted, source-stroke, sink-stroke, fermion-flow, wave, coil, zigzag
#import graph: build, edge, node, sink, source
#set text(fill: palette.ink)
#let particle-map = (
"##,
    );
    let particles: BTreeSet<_> = incoming
        .iter()
        .chain(outgoing.iter().flatten())
        .copied()
        .collect();
    if particles.is_empty() {
        source.push(':');
    }
    for id in particles {
        let p = model
            .particle_by_id(id)
            .expect("validated process particle");
        writeln!(
            source,
            "{:?}: {},",
            p.name,
            p.generate_edge_typst_dict(model)
        )
        .unwrap();
    }
    source.push_str(")\n#context { stack(dir: ttb, spacing: 10pt,\n");
    for state in outgoing {
        let names = |state: &[ParticleId]| {
            state
                .iter()
                .map(|id| model.particle_by_id(*id).unwrap().name.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        };
        let initial = names(incoming);
        let final_state = names(state);
        let summary = format!("{initial} → {final_state}");
        writeln!(source, "{{ let raw = build({{\n node(<process>, inspection: (title: \"Process\", summary: {summary:?}, properties: ((\"Model\", {:?}), (\"Incoming\", {initial:?}), (\"Outgoing\", {final_state:?}))))", model.name()).unwrap();
        for (index, (id, is_incoming)) in incoming
            .iter()
            .map(|id| (id, true))
            .chain(state.iter().map(|id| (id, false)))
            .enumerate()
        {
            let p = model
                .particle_by_id(*id)
                .expect("validated process particle");
            let endpoints = if is_incoming {
                format!("<e{index}>, sink(<process>)")
            } else {
                format!("source(<process>), <e{index}>")
            };
            let orientation = if p.antiparticle == *id {
                "undirected"
            } else if p.is_antiparticle() {
                "reversed"
            } else {
                "default"
            };
            writeln!(
                source,
                "edge({endpoints}, particle: {:?}, orientation: {orientation:?}, inspection: (title: {:?}, summary: {:?}, properties: ((\"PDG\", {:?}), (\"Mass\", {:?}))))",
                p.name, p.name, format!("{} leg {}: {}", if is_incoming { "Incoming" } else { "Outgoing" }, index + 1, p.name), p.pdg_code.to_string(), model.parameter_by_id(p.mass).unwrap().name
            )
            .unwrap();
        }
        source.push_str(r#"})
 let config = _linnet_config
 let title = config.at("title", default: auto)
 let config = config + (title: if title == auto { none } else { title })
 let options = (mode: "amplitude", label-fill: palette.ink) + config.at("options", default: (:))
 let style = (node-label: none, node-style: (radius: 3.0, fill: none, stroke: (paint: palette.ink, thickness: 0.7pt, dash: "dashed"))) + config.at("style", default: (:))
 physics-layout.render-layout(config + (options: options, style: style), input: raw, graph: graph, renderer: renderer, physics: physics, edge-style: (map: particle-map, default-edge: physics.default-edge))
},
"#);
    }
    source.push_str(") }");
    let prepared = prepare_physics_render(py, &source, config)?;
    render_diagram_svg(py, &prepared)
}

/// Values are escaped text or fragments from Symbolica's native expression printer.
pub(crate) fn model_record_html(kind: &str, name: &str, rows: &[(&str, String)]) -> String {
    let mut html = format!(
        r#"<style>.feynkit-record [data-spenso-math]{{display:inline-block;vertical-align:middle;padding:0;}}</style><section class="feynkit-record" style="max-width:100%;overflow-x:auto;padding:.6rem 0;color:inherit"><div style="margin-bottom:.55rem"><strong>{}</strong><span style="margin-left:.6rem;opacity:.65">{}</span></div><table style="border-collapse:collapse;width:auto;max-width:100%"><tbody>"#,
        escape_html(name),
        escape_html(kind)
    );
    for (label, value) in rows {
        write!(html, r#"<tr><th scope="row" style="text-align:left;vertical-align:top;font-weight:500;opacity:.7;padding:.25rem 1rem .25rem 0;white-space:nowrap">{}</th><td style="text-align:left;padding:.25rem 0;overflow-wrap:anywhere">{value}</td></tr>"#, escape_html(label)).unwrap();
    }
    html.push_str("</tbody></table></section>");
    html
}

/// Model formulas use Spenso's existing MathML/Typst printer without recooking indices.
pub(crate) fn model_expression_html(
    py: Python<'_>,
    expression: symbolica::api::python::PythonExpression,
) -> PyResult<String> {
    py.import("symbolica.community.spenso")?
        .getattr("to_html")?
        .call1((expression,))?
        .extract()
}
