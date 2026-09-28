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

pub(crate) const PREVIEW_LIMIT: usize = 6;

/// Compact collections share the graph renderer and Spenso's expression printer.
/// SVG templates remain inert until a row or thumbnail is selected.
pub(crate) fn collection_html(
    py: Python<'_>,
    title: &str,
    subtitle: &str,
    diagrams: impl ExactSizeIterator<Item = crate::graph::PyFeynmanDiagram>,
    terms: Option<&[String]>,
) -> PyResult<String> {
    let count = diagrams.len();
    let mut html = format!(
        "<style>{}</style><section class=\"feynkit-collection\"><header><strong>{}</strong><small>{}</small></header>",
        include_str!("collection.css"),
        escape_html(title),
        escape_html(subtitle),
    );
    let mut templates = String::new();
    if terms.is_none() {
        html.push_str("<div class=\"fk-strip\" role=\"group\" aria-label=\"Choose diagram\">");
    }
    for (index, diagram) in diagrams.take(PREVIEW_LIMIT).enumerate() {
        let svg = diagram.render(py, None, false, None, None)?;
        write!(templates, "<template>{svg}</template>").unwrap();
        let label = format!("Diagram {}", index + 1);
        let mut caption = format!(
            "{} · {} loop{}",
            diagram.inner.name(),
            diagram.inner.loop_count(),
            if diagram.inner.loop_count() == 1 {
                ""
            } else {
                "s"
            },
        );
        let cuts = diagram.inner.cuts().len();
        if cuts > 0 {
            write!(caption, " · {cuts} cut{}", if cuts == 1 { "" } else { "s" }).unwrap();
        }
        let caption = escape_html(&caption);
        let thumbnail = format!("<img class=\"fk-thumbnail\" alt=\"{label}\">");
        if let Some(terms) = terms {
            write!(html, "<details class=\"fk-row\"><summary>{thumbnail}<div class=\"fk-term\"><small>{label} · Contribution</small>{}</div></summary><div class=\"fk-stage\"></div><div class=\"fk-caption\">{caption}</div></details>", terms[index]).unwrap();
        } else {
            write!(html, "<button type=\"button\" aria-pressed=\"false\" data-caption=\"{caption}\">{thumbnail}{label}</button>").unwrap();
        }
    }
    if terms.is_none() {
        html.push_str("</div><div class=\"fk-stage\"></div><div class=\"fk-caption\" aria-live=\"polite\"></div>");
    }
    if count == 0 {
        html.push_str("<p>No diagrams retained.</p>");
    }
    if count > PREVIEW_LIMIT {
        write!(html, "<small>Showing {PREVIEW_LIMIT} of {count} diagrams. Access .diagrams to inspect the complete collection.</small>").unwrap();
    }
    write!(
        html,
        "{templates}</section><script>{}</script>",
        include_str!("collection.js")
    )
    .unwrap();
    Ok(html)
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
#import "crates/linnest/typst/src/draw.typ" as drawing-api
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
        writeln!(
            source,
            "}})\n let process-mode = {:?}",
            if incoming.is_empty() || state.is_empty() {
                "generic"
            } else {
                "amplitude"
            }
        )
        .unwrap();
        source.push_str(r#" let config = _linnet_config
 let title = config.at("title", default: auto)
 let config = config + (title: if title == auto { none } else { title })
 let options = (mode: process-mode, label-fill: palette.ink) + config.at("options", default: (:))
 let drawing = config.at("draw", default: (:))
 let blob-fill = tiling(size: (5pt, 5pt), {
   place(line(start: (0pt, 5pt), end: (5pt, 0pt), stroke: palette.ink + 0.35pt))
 })
 let style = (node-label: none, node-style: (radius: drawing.at("node-radius", default: 3.0), fill: blob-fill, stroke: (paint: palette.ink, thickness: 0.7pt, dash: "dashed"))) + config.at("style", default: (:))
 // Size the process springs from the final blob, including element styles.
 let effective = renderer._effective-options((unit: 1.5, node-label-style: (padding: 0.08)), style, drawing)
 let callbacks = renderer._callbacks(effective)
 let measured = graph.style(renderer.attach-elements(raw, config.at("elements", default: (:))), ..(effective + callbacks))
 let node = graph.nodes(measured).first()
 let bounds = node.statements
 let node-data = graph._impl._node-style-data(node, effective)
 let fitted-radius = drawing-api._impl._node-radius(
   (length: graph._impl._canvas-length(effective.unit).to-absolute()),
   (callbacks.node-label)(node-data), (callbacks.node-style)(node-data),
   drawing.at("node-min-radius", default: 0.16), drawing.at("node-label-padding", default: 0.08),
 )
 let outset = drawing.at("node-outset", default: auto)
 let radius = calc.max(1.0, drawing-api._impl._radius-outset(fitted-radius), float(bounds.at("layout-width")) / 2, float(bounds.at("layout-height")) / 2, if outset == auto { 0 } else { outset })
 // At the default 10 x 10 viewport, dangling spring rest length is 20 times length-scale.
 // Repulsion leaves a visible leg beyond the blob; mixed states keep left/right groups.
 let defaults = (length-scale: radius / 20, external-centroid-distance: 1.5)
 let layouts = renderer._layout-passes(config.at("layouts", default: ((:),))).map(pass => defaults + pass)
 physics-layout.render-layout(config + (options: options, style: style, layouts: layouts), input: raw, graph: graph, renderer: renderer, physics: physics, edge-style: (map: particle-map, default-edge: physics.default-edge))
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
