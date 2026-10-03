use std::collections::{BTreeMap, BTreeSet};

impl super::Scene {
    /// Decorate only the identity links emitted by Linnest, leaving its drawing intact.
    pub fn interactive_svg(svg: &str) -> Result<String, String> {
        let document = roxmltree::Document::parse(svg)
            .map_err(|error| format!("invalid Typst SVG: {error}"))?;
        let root = document.root_element();
        let mut edits = Vec::new();
        let mut seen = BTreeSet::new();
        let mut details = BTreeMap::new();
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
            // An edge's many hit targets carry the same inspection payload.
            // Cache by the full href: one identity may have distinct details.
            let (title, detail) = details.entry(href).or_insert_with(|| {
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
                (title, detail)
            });
            let tabindex = if seen.insert(identity) { 0 } else { -1 };
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
        let capacity = edits.iter().fold(svg.len(), |size, (range, replacement)| {
            size + replacement.len() - range.len()
        });
        let mut output = String::with_capacity(capacity);
        let mut cursor = 0;
        for (range, replacement) in edits {
            output.push_str(&svg[cursor..range.start]);
            output.push_str(&replacement);
            cursor = range.end;
        }
        output.push_str(&svg[cursor..]);
        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    use super::super::Scene;

    #[test]
    fn repeated_targets_retain_exact_markup_and_identity_focus_order() {
        // These openings are the previous decorator's byte-for-byte output,
        // including whitespace left behind when href attributes are removed.
        let links = [
            (
                r##"<a href="#linnet-edge-2?{&quot;summary&quot;:&quot;&lt;&amp;&gt;&quot;}" transform="translate(1 2)">"##,
                r#"<a  transform="translate(1 2)" role="button" tabindex="0" aria-label="Edge 2 · &lt;&amp;&gt;" aria-pressed="false" data-linnet-kind="edge" data-linnet-id="2" data-linnet-detail="{&quot;summary&quot;:&quot;&lt;&amp;&gt;&quot;}"><title>Edge 2 · &lt;&amp;&gt;</title>"#,
            ),
            (
                r##"<a href="#linnet-edge-2?{&quot;summary&quot;:&quot;&lt;&amp;&gt;&quot;}">"##,
                r#"<a  role="button" tabindex="-1" aria-label="Edge 2 · &lt;&amp;&gt;" aria-pressed="false" data-linnet-kind="edge" data-linnet-id="2" data-linnet-detail="{&quot;summary&quot;:&quot;&lt;&amp;&gt;&quot;}"><title>Edge 2 · &lt;&amp;&gt;</title>"#,
            ),
            (
                r##"<a href="#linnet-edge-2?{&quot;summary&quot;:&quot;different&quot;}">"##,
                r#"<a  role="button" tabindex="-1" aria-label="Edge 2 · different" aria-pressed="false" data-linnet-kind="edge" data-linnet-id="2" data-linnet-detail="{&quot;summary&quot;:&quot;different&quot;}"><title>Edge 2 · different</title>"#,
            ),
            (
                r##"<a xlink:href="#linnet-halfedge-03?[]" class="target">"##,
                r#"<a  class="target" role="button" tabindex="0" aria-label="Half-edge 03" aria-pressed="false" data-linnet-kind="halfedge" data-linnet-id="03" data-linnet-detail="{}"><title>Half-edge 03</title>"#,
            ),
            (
                r##"<a href="#linnet-node-7?invalid-json">"##,
                r#"<a  role="button" tabindex="0" aria-label="Node 7" aria-pressed="false" data-linnet-kind="node" data-linnet-id="7" data-linnet-detail="{}"><title>Node 7</title>"#,
            ),
        ];
        let mut source = String::from(r#"<svg xmlns:xlink="http://www.w3.org/1999/xlink">"#);
        let mut expected = String::from(
            r#"<svg data-linnet-interactive="true" xmlns:xlink="http://www.w3.org/1999/xlink">"#,
        );
        for (original, decorated) in links {
            source.push_str(original);
            expected.push_str(decorated);
            source.push_str("<rect width=\"8\" height=\"8\"/></a>");
            expected.push_str("<rect width=\"8\" height=\"8\"/></a>");
        }
        source.push_str("</svg>");
        expected.push_str(&format!(
            "<style>{}</style><script><![CDATA[{}]]></script></svg>",
            include_str!("interactive.css"),
            include_str!("interactive.js"),
        ));
        assert_eq!(Scene::interactive_svg(&source).unwrap(), expected);
    }

    #[test]
    fn unrelated_links_and_xml_errors_are_unchanged() {
        let source = r##"<svg><a href="https://example.com"><path/></a><a href="#linnet-other-2"/><a href="#linnet-edge-x"/><a href="#linnet-node-"/></svg>"##;
        assert_eq!(Scene::interactive_svg(source).unwrap(), source);
        let invalid = "<svg><a></svg>";
        let error = roxmltree::Document::parse(invalid).unwrap_err();
        assert_eq!(
            Scene::interactive_svg(invalid).unwrap_err(),
            format!("invalid Typst SVG: {error}"),
        );
    }
}
