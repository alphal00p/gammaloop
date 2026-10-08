//! CFF presentation owns only inspection state; native trees and graph IDs own the data.
use std::collections::BTreeMap;

use feynkit_graph::SceneOptions;
use linnet_py::PyDiagramRender;
use pyo3::prelude::*;
use serde_json::{Value, json};

use super::{PyCrossFreeFamily, SurfaceId};
use crate::{display::escape_html, error};

impl PyCrossFreeFamily {
    fn explorer_data(&self) -> Value {
        // Energy and H indices occupy separate arenas. Give the browser one dense
        // key, while keeping the native category and index for the displayed label.
        let surfaces = self.surfaces();
        let keys: BTreeMap<_, _> = surfaces
            .iter()
            .enumerate()
            .map(|(key, s)| (s.id, key))
            .collect();
        let orientations: Vec<_> = self.orientations().into_iter().map(|orientation| {
            let terms: Vec<Vec<_>> = orientation.inner.denominator_products().into_iter()
                .filter(|term| !term.contains(&SurfaceId::Infinite))
                .map(|term| term.into_iter().filter(|id| *id != SurfaceId::Unit).map(|id| keys[&id]).collect())
                .collect();
            json!({"id": orientation.id(), "directions": orientation.edge_orientations().into_iter().collect::<BTreeMap<_, _>>(), "terms": terms})
        }).collect();
        json!({
            "name": self.diagram.name(),
            "loops": self.diagram.loop_count(),
            "orientations": orientations,
            "surfaces": surfaces.iter().map(|s| json!({
                "index": s.index().expect("surface arenas contain no sentinels"),
                "kind": s.kind(), "e": s.positive_energies(), "negative": s.negative_energies(),
                "q": s.external_shift(), "v": s.vertices(),
            })).collect::<Vec<_>>(),
            "edges": self.diagram.edges().map(|(id, ends, edge)| json!({
                "id": id.0, "source": ends.source.map(|v| v.0), "target": ends.target.map(|v| v.0),
                "external": edge.external.is_some(),
            })).collect::<Vec<_>>(),
        })
    }

    pub(super) fn explorer_html(&self, py: Python<'_>) -> PyResult<String> {
        if self.drawing.get().is_none() {
            let options = SceneOptions {
                show_particle: false,
                show_edge_index: true,
                show_node_index: true,
                split_initial_state: false,
                ..SceneOptions::default()
            };
            let mut scene = self
                .diagram
                .to_scene(None, &Default::default(), None, &options)
                .map_err(error::diagram)?;
            scene.title = None;
            // This view shows energy flow. Particle-flow arrows would compete with
            // the editable orientation arrows; keep Linnet's geometry and labels.
            for edge in &mut scene.edges {
                edge.flow = None;
                edge.pattern = None;
            }
            let svg = PyDiagramRender::from_scene(py, scene, None)?.to_svg()?;
            let _ = self.drawing.set(svg);
        }
        // JSON is inert data, but a literal closing script tag still ends an HTML
        // raw-text element. Escape '<' even in graph names supplied by the user.
        let data = self.explorer_data().to_string().replace('<', "\\u003c");
        Ok(format!(
            "<style>{}</style><section class=\"feynkit-cff-result\" data-feynkit-notebook data-linnet-frame-owner><header class=\"hs-meta\"><code title=\"{}\">CrossFreeFamily</code><span data-metadata></span></header>{}<template data-drawing>{}</template><script type=\"application/json\" data-cff>{data}</script><noscript>Enable notebook JavaScript to explore orientations and families. Use to_expression() for the symbolic expression.</noscript></section><script>{}</script><script>{}</script>",
            include_str!("display.css"),
            escape_html(self.diagram.name()),
            include_str!("display.html"),
            self.drawing.get().expect("rendered above"),
            include_str!("display.js"),
            include_str!("../notebook.js"),
        ))
    }
}
