use std::collections::BTreeSet;

impl super::Scene {
    /// Decorate only the identity links emitted by Linnest, leaving its drawing intact.
    pub fn interactive_svg(svg: &str) -> Result<String, String> {
        let document = roxmltree::Document::parse(svg)
            .map_err(|error| format!("invalid Typst SVG: {error}"))?;
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
                include_str!("interactive.css"),
                include_str!("interactive.js")
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
