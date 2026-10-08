use feynkit_graph::FeynmanDiagram;
use feynkit_model::{Model, ParticleId};
use linnest::svg::Scene;
use linnet_py::PyRenderSettings;

use crate::render_settings::PyDiagramStyle;
use pyo3::prelude::*;
use std::fmt::Write;

use linnet_py::PyDiagramRender;

/// Energy representations reuse Linnet's layout with their own orientation marks.
pub(crate) fn energy_graph_svg(
    py: Python<'_>,
    diagram: &FeynmanDiagram,
    style: &feynkit_graph::SceneOptions,
) -> PyResult<String> {
    let options = feynkit_graph::SceneOptions {
        show_particle: false,
        show_edge_index: true,
        show_node_index: true,
        split_initial_state: false,
        show_momentum: Some(false),
        ..style.clone()
    };
    let mut scene = diagram
        .to_scene(None, &Default::default(), None, &options)
        .map_err(crate::error::diagram)?;
    scene.title = None;
    // Particle-flow arrows would compete with the energy-routing annotations.
    for edge in &mut scene.edges {
        edge.flow = None;
        edge.pattern = None;
    }
    PyDiagramRender::from_scene(py, scene, None)?.to_svg()
}

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

/// Use the website SVG palette from docs/assets/typst/theme.typ, including
/// the 45% lightened sink strokes from the shared physics style. Keep SVG
/// paint attributes as light-mode fallbacks for viewers without CSS support.
pub(crate) fn themed_svg(svg: String) -> String {
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
    svg.replacen("<svg ", "<svg class=\"feynkit-diagram-svg\" ", 1)
        .replacen('>', &format!(">{styles}"), 1)
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
#[allow(clippy::too_many_arguments)]
pub(crate) fn collection_html(
    py: Python<'_>,
    title: &str,
    subtitle: &str,
    diagrams: impl ExactSizeIterator<Item = crate::graph::PyFeynmanDiagram>,
    terms: Option<&[String]>,
    config: Option<&PyRenderSettings>,
    style: Option<&PyDiagramStyle>,
    limit: usize,
) -> PyResult<(String, Vec<PyDiagramRender>)> {
    let count = diagrams.len();
    let mut html = format!(
        "<style>{}</style><section class=\"feynkit-collection\" data-feynkit-notebook data-linnet-frame-owner><header><strong>{}</strong><small>{}</small></header>",
        include_str!("collection.css"),
        escape_html(title),
        escape_html(subtitle),
    );
    let mut templates = String::new();
    if terms.is_none() {
        html.push_str("<div class=\"fk-strip\" role=\"group\" aria-label=\"Choose diagram\">");
    }
    let mut drawings = Vec::new();
    for (index, diagram) in diagrams.take(limit).enumerate() {
        let drawing = diagram.render(py, config, style, false, None, None)?;
        let svg = drawing.to_svg_pages()?.join("\n");
        write!(templates, "<template>{svg}</template>").unwrap();
        drawings.push(drawing);
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
    if count > limit {
        write!(html, "<small>Showing {limit} of {count} diagrams. The source collection retains all diagrams.</small>").unwrap();
    }
    write!(
        html,
        "{templates}</section><script>{}</script><script>{}</script>",
        include_str!("collection.js"),
        include_str!("notebook.js")
    )
    .unwrap();
    Ok((html, drawings))
}

/// Draw each process channel as a native star graph with a hatched interaction blob.
pub(crate) fn process_svg(
    model: &Model,
    incoming: &[ParticleId],
    outgoing: &[Vec<ParticleId>],
    config: Option<&PyRenderSettings>,
    style: Option<&PyDiagramStyle>,
) -> PyResult<String> {
    let options = PyDiagramStyle::resolve(style, false);
    let config = config.cloned().unwrap_or_default().config()?;
    let mut figures = Vec::new();
    for state in outgoing {
        let scene = options
            .process_scene(model, incoming, state, &config)
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        figures.push(Python::attach(|py| {
            PyDiagramRender::from_scene(py, scene, None)?
                .map_svg(themed_svg)?
                .to_svg()
        })?);
    }
    Scene::combine_svgs(&figures).map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

/// Values are escaped text or fragments from Symbolica's native expression printer.
pub(crate) fn record_html(kind: &str, name: &str, rows: &[(&str, String)]) -> String {
    let mut html = format!(
        r#"<style>.feynkit-record [data-spenso-math]{{display:inline-block;vertical-align:middle;padding:0;overflow:visible;}}</style><section class="feynkit-record" style="max-width:100%;overflow-x:auto;padding:.6rem 0;color:inherit"><div style="margin-bottom:.55rem"><strong>{}</strong><span style="margin-left:.6rem;opacity:.65">{}</span></div><table style="border-collapse:collapse;width:auto;max-width:100%"><tbody>"#,
        escape_html(name),
        escape_html(kind)
    );
    for (label, value) in rows {
        write!(html, r#"<tr><th scope="row" style="text-align:left;vertical-align:top;font-weight:500;opacity:.7;padding:.25rem 1rem .25rem 0;white-space:nowrap">{}</th><td style="text-align:left;padding:.25rem 0;overflow-wrap:anywhere">{value}</td></tr>"#, escape_html(label)).unwrap();
    }
    html.push_str("</tbody></table></section>");
    html
}

/// Domain formulas use Spenso's existing MathML/Typst printer without recooking indices.
pub(crate) fn expression_html(
    py: Python<'_>,
    expression: symbolica::api::python::PythonExpression,
) -> PyResult<String> {
    py.import("symbolica.community.tensor")?
        .getattr("to_html")?
        .call1((expression,))?
        .extract()
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
