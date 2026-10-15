use crate::PyDiagramRender;
use krilla::{geom::Size, page::PageSettings, Document};
use krilla_svg::{SurfaceExt, SvgSettings};
use pyo3::{exceptions::PyRuntimeError, PyResult};
use resvg::{tiny_skia, usvg};

impl PyDiagramRender {
    /// Convert the captured SVG geometry directly; labels are already glyph paths.
    pub(crate) fn export_pages(&self, format: &str) -> PyResult<Vec<Vec<u8>>> {
        let options = usvg::Options {
            // Native SVG uses points; retain the same physical page dimensions.
            dpi: 72.0,
            ..Default::default()
        };
        let trees = self
            .pages()
            .iter()
            .map(|svg| {
                usvg::Tree::from_str(svg, &options)
                    .map_err(|error| PyRuntimeError::new_err(error.to_string()))
            })
            .collect::<PyResult<Vec<_>>>()?;
        if format == "pdf" {
            let mut document = Document::new();
            for tree in &trees {
                let size = Size::from_wh(tree.size().width(), tree.size().height())
                    .ok_or_else(|| PyRuntimeError::new_err("invalid SVG page size"))?;
                let mut page = document.start_page_with(PageSettings::new(size));
                let mut surface = page.surface();
                surface.draw_svg(tree, size, SvgSettings::default());
                surface.finish();
                page.finish();
            }
            return document
                .finish()
                .map(|pdf| vec![pdf])
                .map_err(|error| PyRuntimeError::new_err(error.to_string()));
        }
        trees
            .iter()
            .map(|tree| {
                let size = tree.size().to_int_size();
                let mut pixels = tiny_skia::Pixmap::new(size.width(), size.height())
                    .ok_or_else(|| PyRuntimeError::new_err("SVG page is too large to rasterize"))?;
                resvg::render(tree, tiny_skia::Transform::identity(), &mut pixels.as_mut());
                pixels
                    .encode_png()
                    .map_err(|error| PyRuntimeError::new_err(error.to_string()))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use linnest::svg::Scene;

    #[test]
    fn svg_inspection_preserves_drawing_and_scopes_identity_links() {
        let drawing = r##"<path d="M 1 2 L 4 5" stroke="#123456"/>"##;
        let svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10">{drawing}<a href="#linnet-node-7" transform="translate(1 2)"><rect width="4" height="4" fill="transparent"/></a><a href="#linnet-node-7"><rect width="2" height="2" fill="transparent"/></a><a href="https://example.org">reference</a><a href="#linnet-node-nope">untouched</a></svg>"##
        );
        let output = Scene::interactive_svg(&svg).unwrap();
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
        let output = Scene::interactive_svg(&svg).unwrap();
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
        assert_eq!(Scene::interactive_svg(svg).unwrap(), svg);
    }

    #[test]
    fn svg_identity_links_remove_both_destination_attributes() {
        let svg = r##"<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"><a href="#linnet-edge-0" xlink:href="#linnet-edge-0"><rect width="3" height="3"/></a></svg>"##;
        let output = Scene::interactive_svg(svg).unwrap();
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
