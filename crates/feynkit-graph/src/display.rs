use std::{
    collections::BTreeMap,
    fmt::Write,
    hash::{DefaultHasher, Hash, Hasher},
};

use linnet::half_edge::{
    NodeIndex,
    involution::EdgeIndex,
    layout::layered::{LayeredConfig, LayeredGeometry, LayeredProfile, LayeredRouteExit},
    subgraph::SuBitGraph,
};

use crate::{ExternalState, FeynmanDiagram};

const WIDTH: f64 = 640.0;
const PADDING: f64 = 58.0;

fn escape_xml(value: &str) -> String {
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

fn canvas_point(point: (f64, f64), bounds: (f64, f64, f64, f64), height: f64) -> (f64, f64) {
    let (min_x, max_x, min_y, max_y) = bounds;
    let range_x = (max_x - min_x).max(1.0);
    let range_y = (max_y - min_y).max(1.0);
    let scale = ((WIDTH - 2.0 * PADDING) / range_x).min((height - 2.0 * PADDING) / range_y);
    let offset_x = (WIDTH - range_x * scale) / 2.0 - min_x * scale;
    let offset_y = (height - range_y * scale) / 2.0 - min_y * scale;
    (point.0 * scale + offset_x, point.1 * scale + offset_y)
}

impl FeynmanDiagram {
    /// Render this diagram as a self-contained, responsive SVG.
    ///
    /// Linnet's deterministic layered layout supplies the node and edge
    /// positions. The SVG renderer has no browser-side or system dependency.
    pub fn to_svg(&self) -> String {
        let graph = self.underlying();
        let incoming: Vec<_> = self
            .vertices()
            .filter_map(|(id, vertex)| {
                vertex
                    .external
                    .as_ref()
                    .filter(|leg| leg.state == ExternalState::Incoming)
                    .map(|_| NodeIndex(id.0))
            })
            .collect();
        let outgoing: Vec<_> = self
            .vertices()
            .filter_map(|(id, vertex)| {
                vertex
                    .external
                    .as_ref()
                    .filter(|leg| leg.state == ExternalState::Outgoing)
                    .map(|_| NodeIndex(id.0))
            })
            .collect();
        let rank_same = [incoming.clone(), outgoing.clone()]
            .into_iter()
            .filter(|rank| rank.len() > 1)
            .collect();
        let config = LayeredConfig {
            profile: LayeredProfile::Dot,
            layer_gap: 2.0,
            node_gap: 1.2,
            edge_gap: 0.45,
            roots: incoming,
            rank_same,
            ..LayeredConfig::default()
        };
        let geometry = LayeredGeometry {
            node_ranks: graph.new_nodevec(|_, _, _| None),
            node_widths: graph.new_nodevec(|_, _, _| 0.4),
            node_heights: graph.new_nodevec(|_, _, _| 0.4),
            edge_label_widths: graph
                .new_edgevec(|edge, _, _| 0.35 + edge.particle.name.chars().count() as f64 * 0.11),
            edge_label_heights: graph.new_edgevec(|_, _, _| 0.25),
            edge_minlens: graph.new_edgevec(|_, _, _| 1),
            edge_weights: graph.new_edgevec(|_, _, _| 1.0),
            edge_source_exits: graph.new_edgevec(|_, _, _| LayeredRouteExit::Auto),
            edge_sink_exits: graph.new_edgevec(|_, _, _| LayeredRouteExit::Auto),
            edge_constrained: graph.new_edgevec(|_, _, _| true),
        };
        let layout = graph.layered_layout(&SuBitGraph::full(graph.n_hedges()), &config, &geometry);

        // Linnet lays ranks down the page. Swap its axes so scattering
        // processes flow from incoming states on the left to outgoing states
        // on the right.
        let raw_positions: BTreeMap<_, _> = layout
            .node_positions
            .iter()
            .filter_map(|(id, point)| point.map(|point| (id.0, (point.y, point.x))))
            .collect();
        let (min_x, max_x, min_y, max_y) = raw_positions.values().fold(
            (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ),
            |(min_x, max_x, min_y, max_y), (x, y)| {
                (min_x.min(*x), max_x.max(*x), min_y.min(*y), max_y.max(*y))
            },
        );
        let bounds = if raw_positions.is_empty() {
            (0.0, 1.0, 0.0, 1.0)
        } else {
            (min_x, max_x, min_y, max_y)
        };
        let height = (self.vertices().count() as f64 * 54.0).clamp(280.0, 720.0);
        let positions: BTreeMap<_, _> = raw_positions
            .iter()
            .map(|(id, point)| (*id, canvas_point(*point, bounds, height)))
            .collect();

        let mut hasher = DefaultHasher::new();
        self.name().hash(&mut hasher);
        for (id, endpoints, edge) in self.edges() {
            (
                id.0,
                endpoints.source.0,
                endpoints.target.0,
                &edge.particle.name,
            )
                .hash(&mut hasher);
        }
        let identifier = format!("feynkit-{:x}", hasher.finish());
        let marker = format!("{identifier}-arrow");
        let mut svg = format!(
            r#"<svg id="{identifier}" class="feynkit-diagram-svg" viewBox="0 0 {WIDTH:.0} {height:.0}" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Feynman diagram {}" style="display:block;width:100%;height:auto;max-width:{WIDTH:.0}px;color:inherit;color-scheme:light dark"><title>Feynman diagram {}</title><defs><marker id="{marker}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="currentColor"/></marker></defs><style>.fk-edge{{fill:none;stroke:currentColor;stroke-width:2;opacity:.82;vector-effect:non-scaling-stroke}}.fk-node{{stroke:currentColor;stroke-width:2;vector-effect:non-scaling-stroke}}.fk-internal{{fill:currentColor}}.fk-external{{fill:Canvas}}.fk-label,.fk-external-label{{fill:CanvasText;stroke:Canvas;stroke-width:4px;paint-order:stroke;stroke-linejoin:round;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px}}.fk-external-label{{font-size:11px}}</style>"#,
            escape_xml(self.name()),
            escape_xml(self.name()),
        );

        let mut parallel_edges: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
        for (id, endpoints, _) in self.edges() {
            let key = if endpoints.source <= endpoints.target {
                (endpoints.source.0, endpoints.target.0)
            } else {
                (endpoints.target.0, endpoints.source.0)
            };
            parallel_edges.entry(key).or_default().push(id.0);
        }

        for (id, endpoints, edge) in self.edges() {
            let Some(&(x1, y1)) = positions.get(&endpoints.source.0) else {
                continue;
            };
            let Some(&(x2, y2)) = positions.get(&endpoints.target.0) else {
                continue;
            };
            let key = if endpoints.source <= endpoints.target {
                (endpoints.source.0, endpoints.target.0)
            } else {
                (endpoints.target.0, endpoints.source.0)
            };
            let siblings = &parallel_edges[&key];
            let sibling_index = siblings
                .iter()
                .position(|edge_id| *edge_id == id.0)
                .unwrap_or(0);
            let parallel_offset =
                (sibling_index as f64 - (siblings.len() as f64 - 1.0) / 2.0) * 26.0;
            let (path, label_x, label_y) = if endpoints.source == endpoints.target {
                let radius = 38.0 + sibling_index as f64 * 18.0;
                (
                    format!(
                        "M {x1:.2} {:.2} C {:.2} {:.2}, {:.2} {:.2}, {x1:.2} {:.2}",
                        y1 - 7.0,
                        x1 + radius,
                        y1 - radius,
                        x1 - radius,
                        y1 - radius,
                        y1 - 7.0,
                    ),
                    x1,
                    y1 - radius - 5.0,
                )
            } else {
                let dx = x2 - x1;
                let dy = y2 - y1;
                let length = (dx * dx + dy * dy).sqrt().max(1.0);
                let raw_control = layout.edge_positions[EdgeIndex(id.0)]
                    .map(|point| canvas_point((point.y, point.x), bounds, height))
                    .unwrap_or(((x1 + x2) / 2.0, (y1 + y2) / 2.0));
                let control_x = raw_control.0 - dy / length * parallel_offset;
                let control_y = raw_control.1 + dx / length * parallel_offset;
                (
                    format!("M {x1:.2} {y1:.2} Q {control_x:.2} {control_y:.2} {x2:.2} {y2:.2}"),
                    (x1 + 2.0 * control_x + x2) / 4.0,
                    (y1 + 2.0 * control_y + y2) / 4.0 - 7.0,
                )
            };
            let arrow = if edge.directed {
                format!(r##" marker-end="url(#{marker})""##)
            } else {
                String::new()
            };
            write!(
                svg,
                r#"<path class="fk-edge" d="{path}"{arrow}><title>edge {}: {} ({} → {})</title></path><text class="fk-label" x="{label_x:.2}" y="{label_y:.2}" text-anchor="middle">{}</text>"#,
                id.0,
                escape_xml(&edge.particle.name),
                endpoints.source.0,
                endpoints.target.0,
                escape_xml(&edge.particle.name),
            )
            .expect("writing to a string cannot fail");
        }

        for (id, vertex) in self.vertices() {
            let Some(&(x, y)) = positions.get(&id.0) else {
                continue;
            };
            let (class, radius) = if vertex.is_external() {
                ("fk-node fk-external", 5)
            } else {
                ("fk-node fk-internal", 7)
            };
            let details = vertex.interaction.as_ref().map_or_else(
                || vertex.name.clone(),
                |interaction| format!("{}; {interaction}", vertex.name),
            );
            write!(
                svg,
                r#"<circle class="{class}" cx="{x:.2}" cy="{y:.2}" r="{radius}"><title>vertex {}: {}</title></circle>"#,
                id.0,
                escape_xml(&details),
            )
            .expect("writing to a string cannot fail");
            if let Some(external) = &vertex.external {
                let (state, label_x, anchor) = match external.state {
                    ExternalState::Incoming => ("in", x - 11.0, "end"),
                    ExternalState::Outgoing => ("out", x + 11.0, "start"),
                };
                write!(
                    svg,
                    r#"<text class="fk-external-label" x="{label_x:.2}" y="{:.2}" text-anchor="{anchor}">{state} {}</text>"#,
                    y + 4.0,
                    external.index,
                )
                .expect("writing to a string cannot fail");
            }
        }
        svg.push_str("</svg>");
        svg
    }

    /// Render this diagram as a self-contained HTML figure.
    pub fn to_html(&self) -> String {
        let loop_description = match self.loop_count() {
            0 => "tree level".to_owned(),
            1 => "1 loop".to_owned(),
            loops => format!("{loops} loops"),
        };
        format!(
            r#"<figure class="feynkit-diagram" style="margin:.35rem 0;max-width:{WIDTH:.0}px">{}<figcaption style="margin-top:.25rem;font-size:.88em;opacity:.78"><strong>{}</strong> · {} · symmetry factor {}</figcaption></figure>"#,
            self.to_svg(),
            escape_xml(self.name()),
            loop_description,
            self.symmetry_factor(),
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::{DiagramEdge, DiagramVertex, ExternalState, FeynmanDiagram, ParticleReference};

    fn diagram(particle_name: &str) -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder("one-to-two");
        let incoming =
            builder.add_vertex(DiagramVertex::external("in", 0, ExternalState::Incoming));
        let outgoing_a =
            builder.add_vertex(DiagramVertex::external("out1", 1, ExternalState::Outgoing));
        let outgoing_b =
            builder.add_vertex(DiagramVertex::external("out2", 2, ExternalState::Outgoing));
        let interaction = builder.add_vertex(DiagramVertex::interaction("v", "V_3"));
        let edge = || DiagramEdge::new(ParticleReference::new(particle_name, 25), false);
        builder.add_edge(incoming, interaction, edge()).unwrap();
        builder.add_edge(interaction, outgoing_a, edge()).unwrap();
        builder.add_edge(interaction, outgoing_b, edge()).unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn svg_is_self_contained_and_labels_the_physics_graph() {
        let svg = diagram("phi").to_svg();
        assert!(svg.starts_with("<svg"));
        assert!(svg.ends_with("</svg>"));
        assert_eq!(svg.matches("class=\"fk-edge\"").count(), 3);
        assert!(svg.contains("phi"));
        assert!(svg.contains("in 0"));
        assert!(svg.contains("out 2"));
        assert!(!svg.contains("<script"));
    }

    #[test]
    fn svg_escapes_model_metadata() {
        let svg = diagram("<script>alert('x')</script>").to_svg();
        assert!(svg.contains("&lt;script&gt;"));
        assert!(!svg.contains("<script>"));
    }
}
