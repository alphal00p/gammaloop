use feynkit_graph::FeynmanDiagram;
use linnet_py::PreparedRender;
use pyo3::prelude::*;

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
pub(crate) fn render_diagram_svg(py: Python<'_>, diagram: &FeynmanDiagram) -> PyResult<String> {
    let svg = PreparedRender::from_sources(
        [
            ("main.typ", diagram.to_linnest()),
            (
                "assets/embedded/drawing/templates/layout-core.typ",
                include_str!("../../../assets/embedded/drawing/templates/layout-core.typ")
                    .to_owned(),
            ),
            (
                "assets/embedded/drawing/templates/physics-edge-style.typ",
                include_str!("../../../assets/embedded/drawing/templates/physics-edge-style.typ")
                    .to_owned(),
            ),
            (
                "assets/embedded/drawing/templates/impl/physics-edge-style.typ",
                include_str!(
                    "../../../assets/embedded/drawing/templates/impl/physics-edge-style.typ"
                )
                .to_owned(),
            ),
        ]
        .into_iter()
        .map(|(path, source)| (path.to_owned(), source.into_bytes()))
        .collect(),
    )?
    .svg(py)?;
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
    styles.push_str(&format!(
        r#".feynkit-diagram-svg{{{light}}}
@media(prefers-color-scheme:dark){{.feynkit-diagram-svg{{{dark}}}}}
:is(.dark,[data-theme="dark"],[data-jp-theme-light="false"]) .feynkit-diagram-svg{{{dark}}}
:is(.light,[data-theme="light"],[data-jp-theme-light="true"]) .feynkit-diagram-svg{{{light}}}
</style>"#
    ));
    Ok(svg
        .replacen("<svg ", "<svg class=\"feynkit-diagram-svg\" ", 1)
        .replacen('>', &format!(">{styles}"), 1))
}

pub(crate) fn render_diagram_html(py: Python<'_>, diagram: &FeynmanDiagram) -> PyResult<String> {
    let svg = render_diagram_svg(py, diagram)?;
    Ok(format!(
        "<figure class=\"feynkit-diagram\" style=\"max-width:100%;margin:.5rem 0\">\
         <div style=\"max-width:100%;overflow-x:auto\">{svg}</div>\
         <figcaption style=\"font-size:.85em;opacity:.75;text-align:center\">Feynman diagram {} \
         ({} loop{})</figcaption></figure>",
        escape_html(diagram.name()),
        diagram.loop_count(),
        if diagram.loop_count() == 1 { "" } else { "s" },
    ))
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
