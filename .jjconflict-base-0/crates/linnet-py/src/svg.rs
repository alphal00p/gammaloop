use std::collections::BTreeSet;

use pyo3::exceptions::PyRuntimeError;
use pyo3::PyResult;

use crate::PreparedRender;

impl PreparedRender {
    /// Decorate only the identity links emitted by Linnest, leaving its drawing intact.
    pub(crate) fn interactive_svg(svg: &str) -> PyResult<String> {
        let document = roxmltree::Document::parse(svg)
            .map_err(|error| PyRuntimeError::new_err(format!("invalid Typst SVG: {error}")))?;
        let root = document.root_element();
        let mut edits = Vec::new();
        let mut seen = BTreeSet::new();
        for node in root.descendants().filter(|node| node.has_tag_name("a")) {
            let Some(href) = node
                .attribute("href")
                .or_else(|| node.attribute(("http://www.w3.org/1999/xlink", "href")))
            else {
                continue;
            };
            let Some(identity) = href.strip_prefix("#linnet-") else {
                continue;
            };
            let (identity, detail) = identity.split_once('?').unwrap_or((identity, "{}"));
            let Some((kind, id)) = identity.split_once('-') else {
                continue;
            };
            if !matches!(kind, "node" | "edge" | "halfedge")
                || id.is_empty()
                || !id.bytes().all(|byte| byte.is_ascii_digit())
            {
                continue;
            }
            let label = match kind {
                "node" => "Node",
                "halfedge" => "Half-edge",
                _ => "Edge",
            };
            let detail = serde_json::from_str::<serde_json::Value>(detail)
                .ok()
                .filter(serde_json::Value::is_object)
                .unwrap_or_else(|| serde_json::json!({}));
            let title = match detail.get("summary").and_then(serde_json::Value::as_str) {
                Some(summary) => format!("{label} {id} · {summary}"),
                None => format!("{label} {id}"),
            };
            let title = title
                .replace('&', "&amp;")
                .replace('"', "&quot;")
                .replace('<', "&lt;")
                .replace('>', "&gt;");
            // Attribute text is encoded, not interpreted as markup or executable source.
            let detail = detail
                .to_string()
                .replace('&', "&amp;")
                .replace('"', "&quot;")
                .replace('<', "&lt;")
                .replace('>', "&gt;");
            let tabindex = if seen.insert(identity.to_owned()) {
                0
            } else {
                -1
            };
            let range = node.range();
            let open_end = range.start + svg[range.clone()].find('>').expect("SVG anchor") + 1;
            // Keep the original transform and link geometry byte-for-byte. Removing the
            // destination prevents fragment navigation when scripts are unavailable.
            let mut opening = svg[range.start..open_end].to_owned();
            let mut destinations = node
                .attributes()
                .filter(|attribute| attribute.name() == "href")
                .map(|attribute| attribute.range())
                .collect::<Vec<_>>();
            destinations.sort_by_key(|range| range.start);
            for attribute in destinations.into_iter().rev() {
                opening.replace_range(
                    attribute.start - range.start..attribute.end - range.start,
                    "",
                );
            }
            opening.insert_str(opening.len() - 1, &format!(
                " role=\"button\" tabindex=\"{tabindex}\" aria-label=\"{title}\" aria-pressed=\"false\" data-linnet-kind=\"{kind}\" data-linnet-id=\"{id}\" data-linnet-detail=\"{detail}\""
            ));
            opening.push_str(&format!("<title>{title}</title>"));
            edits.push((range.start..open_end, opening));
        }
        if edits.is_empty() {
            return Ok(svg.to_owned());
        }
        let start = root.range().start + "<svg".len();
        edits.push((start..start, " data-linnet-interactive=\"true\"".to_owned()));
        let end = root.range().end - "</svg>".len();
        edits.push((
            end..end,
            format!(
                "<style>{}</style><script><![CDATA[{}]]></script>",
                include_str!("svg.css"),
                include_str!("svg.js")
            ),
        ));
        edits.sort_by_key(|(range, _)| range.start);
        let mut output = svg.to_owned();
        for (range, replacement) in edits.into_iter().rev() {
            output.replace_range(range, &replacement);
        }
        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    use crate::PreparedRender;

    #[test]
    fn svg_inspection_preserves_drawing_and_scopes_identity_links() {
        let drawing = r##"<path d="M 1 2 L 4 5" stroke="#123456"/>"##;
        let svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10">{drawing}<a href="#linnet-node-7" transform="translate(1 2)"><rect width="4" height="4" fill="transparent"/></a><a href="#linnet-node-7"><rect width="2" height="2" fill="transparent"/></a><a href="https://example.org">reference</a><a href="#linnet-node-nope">untouched</a></svg>"##
        );
        let output = PreparedRender::interactive_svg(&svg).unwrap();
        assert!(output.contains(drawing));
        let document = roxmltree::Document::parse(&output).unwrap();
        let targets = document
            .descendants()
            .filter(|node| node.attribute("data-linnet-id") == Some("7"))
            .collect::<Vec<_>>();
        assert_eq!(targets.len(), 2);
        assert_eq!(targets[0].attribute("tabindex"), Some("0"));
        assert_eq!(targets[1].attribute("tabindex"), Some("-1"));
        assert_eq!(targets[0].attribute("transform"), Some("translate(1 2)"));
        assert_eq!(targets[0].attribute("href"), None);
        assert!(output.contains("href=\"https://example.org\""));
        assert!(output.contains("href=\"#linnet-node-nope\""));
    }

    #[test]
    fn svg_hover_summaries_are_escaped_and_preserve_domain_details() {
        let summary = "Stored tensor · A<\"x\">\nRank: 2";
        let details = serde_json::json!({"summary": summary, "properties": [["Rank", "2"]]});
        let encoded = details
            .to_string()
            .replace('&', "&amp;")
            .replace('"', "&quot;")
            .replace('<', "&lt;")
            .replace('>', "&gt;");
        let svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg"><a href="#linnet-node-0?{encoded}"><rect width="3" height="3"/></a></svg>"##
        );
        let output = PreparedRender::interactive_svg(&svg).unwrap();
        let document = roxmltree::Document::parse(&output).unwrap();
        let title = document
            .descendants()
            .find(|node| node.has_tag_name("title"))
            .unwrap();
        assert_eq!(title.text(), Some(format!("Node 0 · {summary}").as_str()));
        let target = title.parent().unwrap();
        let restored: serde_json::Value =
            serde_json::from_str(target.attribute("data-linnet-detail").unwrap()).unwrap();
        assert_eq!(restored, details);
    }

    #[test]
    fn svg_without_graph_targets_is_unchanged() {
        let svg = "<svg xmlns=\"http://www.w3.org/2000/svg\"><path d=\"M0 0\"/></svg>";
        assert_eq!(PreparedRender::interactive_svg(svg).unwrap(), svg);
    }

    #[test]
    fn svg_identity_links_remove_both_destination_attributes() {
        let svg = r##"<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"><a href="#linnet-edge-0" xlink:href="#linnet-edge-0"><rect width="3" height="3"/></a></svg>"##;
        let output = PreparedRender::interactive_svg(svg).unwrap();
        let document = roxmltree::Document::parse(&output).unwrap();
        let target = document
            .descendants()
            .find(|node| node.has_tag_name("a"))
            .unwrap();
        assert!(target
            .attributes()
            .all(|attribute| attribute.name() != "href"));
        assert_eq!(target.attribute("data-linnet-kind"), Some("edge"));
    }
}
