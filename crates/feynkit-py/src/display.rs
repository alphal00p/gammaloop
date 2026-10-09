use feynkit_graph::{FeynmanDiagram, SceneOptions};
use feynkit_model::{Model, ParticleId};
use linnest::svg::{Config, Scene};

use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyFloat, PyInt, PyList, PyString, PyTuple};
use std::{collections::BTreeMap, fmt::Write};

/// Identify data that the linked adapter can read without extension type identity.
fn is_json_data(value: &Bound<'_, PyAny>, depth: usize) -> PyResult<bool> {
    if depth > 64 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "native mark data exceeds maximum nesting depth 64",
        ));
    }
    if value.is_none()
        || value.is_instance_of::<PyBool>()
        || value.is_instance_of::<PyInt>()
        || value.is_instance_of::<PyFloat>()
        || value.is_instance_of::<PyString>()
    {
        return Ok(true);
    }
    if let Ok(dict) = value.cast::<PyDict>() {
        for (key, value) in dict.iter() {
            if !key.is_instance_of::<PyString>() || !is_json_data(&value, depth + 1)? {
                return Ok(false);
            }
        }
        return Ok(true);
    }
    if let Ok(list) = value.cast::<PyList>() {
        for value in list.iter() {
            if !is_json_data(&value, depth + 1)? {
                return Ok(false);
            }
        }
        return Ok(true);
    }
    if let Ok(tuple) = value.cast::<PyTuple>() {
        for value in tuple.iter() {
            if !is_json_data(&value, depth + 1)? {
                return Ok(false);
            }
        }
        return Ok(true);
    }
    Ok(false)
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

/// Parse the native, serializable SVG configuration.
pub(crate) fn render_config(py: Python<'_>, config: Option<&Bound<'_, PyAny>>) -> PyResult<Config> {
    let mut config = config.cloned();
    if let Some(value) = &config
        && let Ok(root) = value.cast::<PyDict>()
        && let Some(style) = root.get_item("style")?
        && let Ok(style) = style.cast::<PyDict>()
        && let Some(edges) = style.get_item("edge-style")?
        && let Ok(edges) = edges.cast::<PyDict>()
    {
        // Copy only the affected levels; rendering must not rewrite caller data.
        let root = root.copy()?;
        let style = style.copy()?;
        let edges = edges.copy()?;
        for key in ["flow-arrow", "momentum-arrow"] {
            if let Some(mark) = edges.get_item(key)? {
                let payload = if is_json_data(&mark, 0)? {
                    let native = linnet_py::native_mark_json(&mark)?;
                    py.import("json")?
                        .call_method1("loads", (native.to_string(),))?
                } else {
                    // The installed extension owns its type objects. Its public
                    // codec removes typed values before crossing into this host.
                    let public_mark = py.import("linnet")?.getattr("Mark")?;
                    let mark = if mark.is_instance(&public_mark)? {
                        mark
                    } else if mark.is_instance_of::<PyDict>() {
                        public_mark.call_method1("from_dict", (&mark,))?
                    } else {
                        return Err(pyo3::exceptions::PyTypeError::new_err(
                            "native marks must be linnet.Mark or public mark dictionaries",
                        ));
                    };
                    mark.call_method0("to_native")?
                };
                let payload = payload.cast::<PyDict>()?;
                edges.set_item(key, payload.get_item("mark")?.unwrap())?;
                edges.set_item(
                    format!("{key}-paints"),
                    payload.get_item("paints")?.unwrap(),
                )?;
            }
        }
        style.set_item("edge-style", edges)?;
        root.set_item("style", style)?;
        config = Some(root.into_any());
    }
    let source = match &config {
        Some(config) => py
            .import("json")?
            .call_method1("dumps", (config,))?
            .extract::<String>()?,
        None => "{}".into(),
    };
    Config::from_json(&source).map_err(pyo3::exceptions::PyValueError::new_err)
}

pub(crate) fn scene_options(config: &Config, momenta: bool) -> PyResult<SceneOptions> {
    let mut options = SceneOptions {
        momentum_arrows: momenta,
        ..Default::default()
    };
    for (key, value) in &config.template_options {
        if key == "mode" && value.as_str() == Some("auto") {
            continue;
        }
        let flag = value.as_bool().ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!("{key} must be boolean"))
        })?;
        match key.as_str() {
            "momentum-arrows" => options.momentum_arrows = flag,
            "show-momentum" => options.show_momentum = Some(flag),
            "show-particle" => options.show_particle = flag,
            "show-edge-index" => options.show_edge_index = flag,
            "show-node-index" => options.show_node_index = flag,
            "split-initial-state" => options.split_initial_state = flag,
            "debug" => {
                options.show_node_index = flag;
                options.show_edge_index = flag;
            }
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "unsupported SVG physics option {key:?}"
                )));
            }
        }
    }
    Ok(options)
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

/// Draw each process channel as a native star graph with a hatched interaction blob.
pub(crate) fn process_svg(
    py: Python<'_>,
    model: &Model,
    incoming: &[ParticleId],
    outgoing: &[Vec<ParticleId>],
    config: Option<&Bound<'_, PyAny>>,
) -> PyResult<String> {
    let config = render_config(py, config)?;
    let options = scene_options(&config, false)?;
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
    use super::{escape_html, render_config};
    use pyo3::{prelude::*, types::PyDict};

    #[test]
    fn escapes_html_metadata() {
        assert_eq!(
            escape_html("<script data-x='a&b'>\"x\"</script>"),
            "&lt;script data-x=&#39;a&amp;b&#39;&gt;&quot;x&quot;&lt;/script&gt;"
        );
    }

    #[test]
    fn render_config_normalizes_marks_without_mutating_caller_styles() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let input = py.eval(
                c"{'style': {'edge-style': {'flow-arrow': {'kind': 'kurvst-mark', 'shape': 'triangle', 'length': '0.21cm', 'fill': 'red'}, 'momentum-arrow': {'kind': 'kurvst-mark', 'shape': 'straight', 'width': '450%', 'stroke': 'blue'}}}}",
                None,
                None,
            )?;
            let before = input.repr()?.to_string();
            let config = render_config(py, Some(&input))?;
            let edges = &config.style["edge-style"];
            assert_eq!(edges["flow-arrow"]["shape"], "triangle");
            assert_eq!(
                edges["flow-arrow"]["length"]["points"],
                serde_json::json!(0.21 * (72.0 / 2.54))
            );
            assert_eq!(edges["flow-arrow-paints"][0]["fill"], "#ff4136");
            assert_eq!(
                edges["momentum-arrow"]["width"],
                serde_json::json!({"points": 0.0, "ratio": 4.5})
            );
            assert_eq!(edges["momentum-arrow-paints"][0]["stroke"], "#0074d9");
            assert_eq!(input.repr()?.to_string(), before);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn render_config_rejects_legacy_marks_and_preserves_unrelated_options() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let legacy = py.eval(
                c"{'style': {'edge-style': {'flow-arrow': {'end': 'straight', 'scale': 0.8}}}}",
                None,
                None,
            )?;
            let error = render_config(py, Some(&legacy)).err().unwrap();
            assert!(error.to_string().contains("CeTZ"));

            let input = PyDict::new(py);
            input.set_item("title", "Configured diagram")?;
            let config = render_config(py, Some(input.as_any()))?;
            assert_eq!(config.title.as_deref(), Some("Configured diagram"));
            assert!(config.style.is_empty());
            Ok(())
        })
        .unwrap();
    }

    #[test]
    #[ignore = "requires the freshly built public linnet extension on PYTHONPATH"]
    fn render_config_accepts_marks_from_the_installed_linnet_extension() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let locals = PyDict::new(py);
            locals.set_item("lp", py.import("linnet")?)?;
            let mark = py.eval(
                c"lp.Mark.combine(lp.Mark.combine(lp.Mark('bar', stroke=lp.Color('blue')), lp.Length.pt(2), lp.Mark('bar')), lp.Ratio.from_fraction(0.5), lp.Mark('stealth', length=lp.RelativeLength(lp.Ratio.from_fraction(4.5), lp.Length.pt(3)), inset=lp.Ratio.from_fraction(0.4), fill=lp.Color.rgb(1, 2, 3, 4)), fit='bend', shorten='80%')",
                None,
                Some(&locals),
            )?;
            let data = mark.call_method0("to_dict")?;
            let edges = PyDict::new(py);
            edges.set_item("flow-arrow", &mark)?;
            edges.set_item("momentum-arrow", &data)?;
            let style = PyDict::new(py);
            style.set_item("edge-style", &edges)?;
            let input = PyDict::new(py);
            input.set_item("style", &style)?;
            input.set_item("title", "Portable mark data")?;
            let before = input.repr()?.to_string();
            let config = render_config(py, Some(input.as_any()))?;
            let normalized = &config.style["edge-style"];
            assert_eq!(normalized["flow-arrow"], normalized["momentum-arrow"]);
            assert_eq!(
                normalized["flow-arrow"],
                serde_json::json!({
                    "shape": "combine", "fit": "bend", "shorten": 0.8,
                    "parts": [
                        {"shape": "combine", "parts": [
                            {"shape": "bar", "stroke": true},
                            {"gap": {"points": 2.0, "ratio": 0.0}},
                            {"shape": "bar"}
                        ]},
                        {"gap": {"points": 0.0, "ratio": 0.5}},
                        {"shape": "stealth", "length": {"points": 3.0, "ratio": 4.5},
                         "inset": 0.4, "fill": true}
                    ]
                })
            );
            for key in ["flow-arrow-paints", "momentum-arrow-paints"] {
                assert_eq!(
                    normalized[key],
                    serde_json::json!([{"stroke": "#0074d9"}, {}, {"fill": "#01020304"}])
                );
            }
            assert_eq!(config.title.as_deref(), Some("Portable mark data"));
            assert_eq!(input.repr()?.to_string(), before);

            // Neither objects advertising a serializer nor callable payloads
            // are part of the closed public data contract.
            py.run(
                c"class Pretender:\n    def to_native(self):\n        raise AssertionError('must not call arbitrary serializer')",
                Some(&locals),
                Some(&locals),
            )?;
            for source in [
                c"Pretender()",
                c"lambda: lp.Mark('triangle')",
                c"{'kind': 'kurvst-mark', 'shape': 'triangle', 'fill': Pretender()}",
                c"{'kind': 'kurvst-mark', 'shape': 'triangle', 'fill': lambda: 'red'}",
            ] {
                let invalid = py.eval(source, Some(&locals), Some(&locals))?;
                edges.set_item("flow-arrow", invalid)?;
                let before = input.repr()?.to_string();
                let error = render_config(py, Some(input.as_any())).err().unwrap();
                assert!(error.is_instance_of::<pyo3::exceptions::PyTypeError>(py));
                assert_eq!(input.repr()?.to_string(), before);
            }
            edges.set_item(
                "flow-arrow",
                py.eval(c"{'end': 'straight', 'scale': 0.8}", None, None)?,
            )?;
            let before = input.repr()?.to_string();
            let error = render_config(py, Some(input.as_any())).err().unwrap();
            assert!(error.to_string().contains("CeTZ"));
            assert_eq!(input.repr()?.to_string(), before);
            Ok(())
        })
        .unwrap();
    }
}
