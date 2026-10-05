use feynkit_graph::FeynmanDiagram;
use feynkit_model::{Model, ParticleId};
use linnest::svg::Scene;

use crate::render_settings::PyRenderSettings;
use pyo3::prelude::*;
use std::{collections::BTreeMap, fmt::Write};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// A rendered diagram snapshot with its configured notebook display.
/// Displaying or exporting the snapshot reuses the rendered SVG and labels.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> process = hep.Model.phi3().process(["phi"], ["phi", "phi"])
/// >>> drawing = process.render(config=hep.RenderSettings(node_radius=5))
/// >>> drawing
/// >>> svg = drawing.to_svg()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "DiagramRender", module = "symbolica.community.hepkit", frozen)]
pub struct PyDiagramRender {
    svg: String,
    html: String,
}

impl PyDiagramRender {
    pub(crate) fn new(svg: String, html: String) -> Self {
        Self { svg, html }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyDiagramRender {
    /// Export the configured SVG, including typeset labels and hover information.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``DiagramRender`` class example:
    ///
    /// >>> from pathlib import Path
    /// >>> Path("diagram.svg").write_text(drawing.to_svg())
    pub(crate) fn to_svg(&self) -> &str {
        &self.svg
    }

    /// Export the notebook HTML figure with its caption and interactive SVG.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``DiagramRender`` class example:
    ///
    /// >>> html = drawing.to_html()
    pub(crate) fn to_html(&self) -> &str {
        &self.html
    }

    /// Export a self-contained Typst document embedding this rendered SVG.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``DiagramRender`` class example:
    ///
    /// >>> from pathlib import Path
    /// >>> Path("diagram.typ").write_text(drawing.to_linnest())
    fn to_linnest(&self) -> String {
        typst_renderer::Document::svg_source(&self.svg)
    }

    /// Display the configured figure in IPython and Jupyter.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``DiagramRender`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(drawing)
    fn _repr_html_(&self) -> &str {
        self.to_html()
    }

    /// Provide the configured SVG to SVG-aware notebook frontends.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``DiagramRender`` class example:
    ///
    /// >>> svg = drawing._repr_svg_()
    fn _repr_svg_(&self) -> &str {
        self.to_svg()
    }

    /// Display the configured HTML figure in Marimo.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``DiagramRender`` class example:
    ///
    /// >>> import marimo as mo
    /// >>> mo.as_html(drawing)
    fn _mime_(&self) -> (&str, &str) {
        ("text/html", self.to_html())
    }

    /// Summarize the rendered snapshot in text-only frontends.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``DiagramRender`` class example:
    ///
    /// >>> repr(drawing)
    /// 'DiagramRender()'
    fn __repr__(&self) -> &str {
        "DiagramRender()"
    }
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

/// Typst handles only label pages; Rust owns all graph layout and geometry.
pub(crate) fn render_scene(scene: &Scene) -> PyResult<String> {
    let runtime = pyo3::exceptions::PyRuntimeError::new_err;
    let files = BTreeMap::from([("main.typ".into(), scene.label_document().into_bytes())]);
    if scene.pages.is_empty() && scene.title.is_none() {
        let svg = scene.render(&Default::default()).map_err(runtime)?;
        return Ok(themed_svg(Scene::interactive_svg(&svg).map_err(runtime)?));
    }
    let pages = typst_renderer::Document::compile_sources(&files, "svg")
        .map_err(runtime)?
        .into_iter()
        .map(String::from_utf8)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| runtime(e.to_string()))?;
    let typeset = scene.typeset(&pages).map_err(runtime)?;
    let svg = scene.render(&typeset).map_err(runtime)?;
    Ok(themed_svg(Scene::interactive_svg(&svg).map_err(runtime)?))
}

/// Use the website SVG palette from docs/assets/typst/theme.typ, including
/// the 45% lightened sink strokes from the shared physics style. Keep SVG
/// paint attributes as light-mode fallbacks for viewers without CSS support.
fn themed_svg(svg: String) -> String {
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
        let drawing = diagram.render(py, None, false, None, None)?;
        let svg = drawing.to_svg();
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

/// Draw each process channel as a native star graph with a hatched interaction blob.
pub(crate) fn process_svg(
    model: &Model,
    incoming: &[ParticleId],
    outgoing: &[Vec<ParticleId>],
    config: Option<&PyRenderSettings>,
) -> PyResult<String> {
    let (config, options) = PyRenderSettings::resolve(config, false);
    let mut figures = Vec::new();
    for state in outgoing {
        let scene = options
            .process_scene(model, incoming, state, &config)
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        figures.push(render_scene(&scene)?);
    }
    Scene::combine_svgs(&figures).map_err(pyo3::exceptions::PyRuntimeError::new_err)
}

/// Values are escaped text or fragments from Symbolica's native expression printer.
pub(crate) fn model_record_html(kind: &str, name: &str, rows: &[(&str, String)]) -> String {
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

/// Model formulas use Spenso's existing MathML/Typst printer without recooking indices.
pub(crate) fn model_expression_html(
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
