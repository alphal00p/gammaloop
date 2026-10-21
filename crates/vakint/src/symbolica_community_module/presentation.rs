//! Notebook rendering of the same records used by Rust's terminal display.

use feynkit_graph::FeynmanDiagram;
use feynkit_py::PyFeynmanDiagram;
use pyo3::{exceptions::PyImportError, prelude::*};
use symbolica::api::python::PythonExpression;

use crate::presentation::{Content, Presentation};

fn escape(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}

fn content(value: &Content) -> PyResult<String> {
    let plain = value.plain();
    let full = match value {
        Content::Text(_) => format!("<span class=\"vk-text\">{}</span>", escape(&plain)),
        Content::Expression(expr) => format!(
            "<div class=\"vk-math\">{}</div>",
            PythonExpression {
                expr: Content::display_expression(expr)
            }
            ._repr_html_()?
        ),
    };
    if plain.chars().count() > 160 {
        Ok(format!(
            "<details class=\"vk-detail\"><summary>{}</summary><div class=\"vk-formula\">{full}</div></details>",
            escape(&Presentation::shortened(&plain, 110))
        ))
    } else {
        Ok(format!("<div class=\"vk-formula\">{full}</div>"))
    }
}

pub(super) fn render(
    py: Python<'_>,
    card: &Presentation,
    graph: Option<&FeynmanDiagram>,
) -> PyResult<String> {
    let mut html = String::from(
        r#"<section class="vakint-card"><style>
.vakint-card{--vk-ink:#18364d;--vk-line:#d5e1e9;--vk-muted:#526777;font-family:system-ui,sans-serif;color:var(--vk-ink);border:1px solid var(--vk-line);border-radius:14px;overflow:hidden;max-width:100%;background:#fff;box-shadow:0 4px 18px #1230440d;margin:.6rem 0}
.vakint-card h3{margin:0;padding:1rem 1.2rem;background:linear-gradient(110deg,#16384e,#28627a);color:#fff;font-size:1.15rem;letter-spacing:.01em}
.vakint-card .vk-fields{display:grid;grid-template-columns:max-content minmax(0,1fr);gap:.55rem 1rem;padding:1rem 1.2rem;margin:0}
.vakint-card dt{color:var(--vk-muted);font-size:.88rem}.vakint-card dd{margin:0;min-width:0}.vakint-card .vk-table-wrap{overflow:auto;max-height:32rem;margin:0 1.2rem 1rem;border:1px solid var(--vk-line);border-radius:8px}
.vakint-card table{border-collapse:collapse;width:100%;font-size:.9rem}.vakint-card th{position:sticky;top:0;background:#edf4f8;color:var(--vk-ink);text-align:left;padding:.65rem;white-space:nowrap}.vakint-card td{border-top:1px solid var(--vk-line);padding:.65rem;vertical-align:top}.vakint-card tbody tr:nth-child(even){background:#f8fbfd}
.vakint-card .vk-formula{max-width:100%;overflow-x:auto}.vakint-card .vk-math{font-family:"STIX Two Math","Cambria Math","DejaVu Math TeX Gyre","Noto Sans Math",serif}.vakint-card .vk-text{white-space:pre-wrap;overflow-wrap:anywhere}.vakint-card .vk-detail summary{cursor:pointer;font-family:ui-monospace,monospace;font-size:.8rem;color:var(--vk-muted);overflow-wrap:anywhere}.vakint-card .vk-detail[open] summary{margin-bottom:.5rem}
.vakint-card .vk-warning{margin:.7rem 1.2rem;padding:.75rem 1rem;border-left:4px solid #bd7414;border-radius:5px;background:#fff5df;color:#74470b;overflow-wrap:anywhere}.vakint-card .vk-graph{padding:.4rem 1.2rem 1rem;overflow:auto}.vakint-card .vk-badge{display:inline-block;border:1px solid #bbd9cd;border-radius:100px;padding:.15rem .55rem;background:#edf8f2;color:#226147;font-size:.85rem}
.vakint-card .vk-graph .feynkit-diagram>div{display:flex;justify-content:center}.vakint-card .vk-graph svg{max-width:100%;max-height:22rem;height:auto}
@media(prefers-color-scheme:dark){.vakint-card{--vk-ink:#e1edf4;--vk-muted:#aac0ce;--vk-line:#395361;background:#172b36}.vakint-card th{background:#263e4b}.vakint-card tbody tr:nth-child(even){background:#1d3340}.vakint-card .vk-warning{background:#41341c;color:#f4d295}.vakint-card .vk-badge{background:#234137;color:#caefdf;border-color:#487563}}
</style>"#,
    );
    html.push_str(&format!(
        "<h3>{}</h3><dl class=\"vk-fields\">",
        escape(&card.title)
    ));
    for (label, value) in &card.fields {
        let rendered = if label == "Available precision" || label == "Status" {
            format!("<span class=\"vk-badge\">{}</span>", escape(&value.plain()))
        } else {
            content(value)?
        };
        html.push_str(&format!("<dt>{}</dt><dd>{rendered}</dd>", escape(label)));
    }
    html.push_str("</dl>");
    if !card.columns.is_empty() {
        html.push_str("<div class=\"vk-table-wrap\"><table><thead><tr>");
        for column in &card.columns {
            html.push_str(&format!("<th>{}</th>", escape(column)));
        }
        html.push_str("</tr></thead><tbody>");
        for row in &card.rows {
            html.push_str("<tr>");
            for cell in row {
                html.push_str(&format!("<td>{}</td>", content(cell)?));
            }
            html.push_str("</tr>");
        }
        html.push_str("</tbody></table></div>");
    }
    for warning in &card.warnings {
        html.push_str(&format!(
            "<div class=\"vk-warning\"><strong>Precision warning</strong><br>{}</div>",
            escape(warning)
        ));
    }
    if let Some(graph) = graph {
        let object = Py::new(py, PyFeynmanDiagram::from(graph.clone()))?;
        match object.bind(py).call_method0("to_html") {
            Ok(value) => html.push_str(&format!(
                "<div class=\"vk-graph\">{}</div>",
                value.extract::<String>()?
            )),
            Err(error) if error.is_instance_of::<PyImportError>(py) => html.push_str(&format!(
                "<div class=\"vk-warning\">Graph preview unavailable: {}</div>",
                escape(&error.to_string())
            )),
            Err(error) => return Err(error),
        }
    }
    html.push_str("</section>");
    Ok(html)
}

#[cfg(test)]
mod tests {
    use super::escape;

    #[test]
    fn escapes_untrusted_names_and_provenance() {
        assert_eq!(escape("<x a='&'>\""), "&lt;x a=&#39;&amp;&#39;&gt;&quot;");
    }
}
