//! Native SVG drawing of prepared graphs, the host-side counterpart of
//! `draw.typ`: ImPrEd layout, Kurvst geometry, the shared annotation search
//! and SVG output. Typst only typesets the label pages a [`Scene`] lists.
mod curves;
mod labels;
mod marks;
mod output;

use kurbo::{BezPath, Point};
use linnet::half_edge::layout::impred::{EdgeLabel, ImpredConfig};
use serde_json::Value;

use crate::{ImpredRun, TypstDotEdge, TypstDotEndpoint, TypstDotNode, TypstGraph, TypstGraphSpec};
use labels::{Anchor, Bounds, Placement, Stroke as CollisionStroke};
pub use labels::{LabelPage, Typeset};
pub use output::Details;

/// Drawing unit in points: Linnest's default 1.5em at 9pt text.
const UNIT: f64 = 13.5;
/// Vertex radius and edge trim, in drawing units.
const NODE_RADIUS: f64 = 0.16;

/// A prepared graph with its paint and label sources, drawn by [`Scene::render`].
pub struct Scene {
    /// Layout graph: sizes and pins as statements, in graph-spec form.
    pub graph: TypstGraphSpec,
    /// Drawing per node and edge, in spec order.
    pub nodes: Vec<NodeDrawing>,
    pub edges: Vec<EdgeDrawing>,
    /// Typst source preceding the label pages: imports, text style, definitions.
    pub preamble: String,
    /// Typst content of the title row above the drawing.
    pub title: Option<String>,
    /// Typst content of each label page.
    pub pages: Vec<String>,
    /// Explicit `impred-*` options; omitted ones keep ImPrEd's defaults.
    pub layout: serde_json::Map<String, Value>,
    /// Refine the layout around drawn labels in short warm ImPrEd passes.
    pub label_feedback: bool,
}

pub struct NodeDrawing {
    pub stroke: Stroke,
    /// Inspection fields after the node's incident edges.
    pub details: Details,
}

pub struct EdgeDrawing {
    pub stroke: Stroke,
    pub pattern: Option<Pattern>,
    /// Particle-flow arrowhead along (`true`) or against the edge orientation.
    pub flow: Option<bool>,
    /// Momentum arrow beside the edge, carrying the label.
    pub momentum: bool,
    pub label: Option<usize>,
    /// Inspection fields after the edge's own identity.
    pub details: Details,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Stroke {
    /// SVG paint, `#rrggbb` or `#rrggbbaa`.
    pub paint: String,
    /// Width in points.
    pub width: f64,
    pub dash: Dash,
    pub round_cap: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Dash {
    Solid,
    /// Dash and gap lengths in points.
    Dashed(f64, f64),
    /// Typst's `"dotted"`: dots one stroke width long, 2pt apart.
    Dotted,
}

/// Line decorations in drawing units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Pattern {
    /// A coil fitted to the complete edge.
    Coil {
        amplitude: f64,
        wavelength: f64,
        longitudinal_scale: f64,
    },
    Wave {
        amplitude: f64,
        wavelength: f64,
    },
    Zigzag {
        amplitude: f64,
        wavelength: f64,
    },
}

impl Pattern {
    fn amplitude(&self) -> f64 {
        match *self {
            Self::Coil { amplitude, .. }
            | Self::Wave { amplitude, .. }
            | Self::Zigzag { amplitude, .. } => amplitude.abs(),
        }
    }
}

/// Warm label feedback: passes, cooling start and steps of each pass.
const FEEDBACK_PASSES: usize = 3;
const FEEDBACK_WARM: f64 = 0.8;
const FEEDBACK_STEPS: usize = 100;

impl Scene {
    /// The Typst document typesetting every label page, in order.
    pub fn label_document(&self) -> String {
        let pages: Vec<String> = self
            .title
            .iter()
            .map(|title| format!("#box[{title}]"))
            .chain(self.pages.iter().map(|page| LabelPage::measured(page)))
            .collect();
        format!(
            "#set page(width: auto, height: auto, margin: 0pt, fill: none)\n{}\n{}",
            self.preamble,
            pages.join("\n#pagebreak()\n")
        )
    }

    /// Read the SVG pages compiled from [`Scene::label_document`].
    pub fn typeset(&self, svgs: &[String]) -> Result<Typeset, String> {
        let typeset = Typeset::read(svgs, self.title.is_some())?;
        if typeset.pages.len() != self.pages.len() {
            return Err("typeset pages do not match the scene's label pages".to_owned());
        }
        Ok(typeset)
    }

    /// Lay out and draw the scene with its typeset pages; one SVG document.
    pub fn render(&self, typeset: &Typeset) -> Result<String, String> {
        let graph = TypstGraph::from_spec(self.graph.clone())?;
        let mut run = ImpredRun::seeded(graph, &Value::Object(self.layout.clone()))?;
        let mut drawing = Drawing::new(self, typeset, &run)?;
        if !self.label_feedback {
            return Ok(drawing.svg);
        }
        // Keep the least-overlapping drawing; later passes may not improve it.
        let (mut best, mut svg) = (drawing.score(), std::mem::take(&mut drawing.svg));
        for _ in 0..FEEDBACK_PASSES {
            if drawing.overlaps.0 == 0 {
                break;
            }
            // Pin each drawn movable internal label beside its carrier and let a
            // short warm pass make room around it.
            run.layout.labels = drawing.pinned_labels(&run);
            run.layout.solve(ImpredConfig {
                labels: true,
                warm_start: FEEDBACK_WARM,
                steps: FEEDBACK_STEPS,
                ..run.config
            })?;
            drawing = Drawing::new(self, typeset, &run)?;
            if drawing.score() < best {
                best = drawing.score();
                svg = std::mem::take(&mut drawing.svg);
            }
        }
        Ok(svg)
    }
}

/// The laid-out graph: node positions, edge records and complete carriers.
struct Laid {
    nodes: Vec<TypstDotNode>,
    edges: Vec<TypstDotEdge>,
    carriers: Vec<Vec<[f64; 2]>>,
}

impl Laid {
    fn new(run: &ImpredRun) -> Result<Self, String> {
        let graph = run
            .session
            .apply_impred_layout(&run.layout)
            .map_err(|error| error.to_string())?;
        Ok(Self {
            nodes: graph.node_records()?,
            edges: graph.edge_records()?,
            carriers: run
                .layout
                .routes
                .iter()
                .map(|route| {
                    route
                        .iter()
                        .map(|&point| run.layout.positions[point])
                        .collect()
                })
                .collect(),
        })
    }

    fn node(&self, index: usize) -> [f64; 2] {
        self.nodes[index]
            .pos
            .as_ref()
            .map_or([0.0; 2], |pos| [pos.x, pos.y])
    }
}

/// Painted elements in drawing order.
enum Element {
    Path {
        path: BezPath,
        stroke: Stroke,
    },
    Chevron([Point; 3]),
    Triangle([Point; 3]),
    Label {
        page: usize,
        bounds: Bounds,
        href: String,
    },
    Node {
        at: [f64; 2],
        stroke: Stroke,
    },
}

/// A transparent hover target.
struct Target {
    at: [f64; 2],
    size: f64,
    href: String,
}

/// The momentum arrow paint: 1pt round ink.
fn arrow_stroke() -> Stroke {
    Stroke {
        paint: output::INK.to_owned(),
        width: 1.0,
        dash: Dash::Solid,
        round_cap: true,
    }
}

struct Drawing {
    svg: String,
    placements: Vec<Placement>,
    chosen: Vec<Bounds>,
    /// Overlapping labels, line hits, label pairs.
    overlaps: (usize, usize, usize),
}

impl Drawing {
    /// Overlapping labels first, then label pairs, then line hits.
    fn score(&self) -> (usize, usize, usize) {
        (self.overlaps.0, self.overlaps.2, self.overlaps.1)
    }

    /// Each drawn movable internal label as a box pinned beside its carrier.
    fn pinned_labels(&self, run: &ImpredRun) -> Vec<Option<EdgeLabel>> {
        let mut labels = vec![None; run.layout.routes.len()];
        for (placement, chosen) in self.placements.iter().zip(&self.chosen) {
            if placement.candidates.len() > 1 && run.layout.external[placement.edge].is_none() {
                labels[placement.edge] = Some(EdgeLabel {
                    extents: [
                        (chosen.right - chosen.left) / 2.0,
                        (chosen.top - chosen.bottom) / 2.0,
                    ],
                    center: Some([
                        (chosen.left + chosen.right) / 2.0,
                        (chosen.bottom + chosen.top) / 2.0,
                    ]),
                    pinned: true,
                    pin: None,
                });
            }
        }
        labels
    }

    fn new(scene: &Scene, typeset: &Typeset, run: &ImpredRun) -> Result<Self, String> {
        let laid = Laid::new(run)?;
        let mut layers = Vec::new();
        let mut targets = Vec::new();
        let mut strokes = Vec::new();
        let mut placements = Vec::new();
        let mut fixed_arrows = Vec::new();
        let mut region_hrefs = Vec::new();
        for (edge, carrier) in laid.edges.iter().zip(&laid.carriers) {
            let drawing = scene
                .edges
                .get(edge.edge)
                .ok_or("scene edges do not match the layout graph")?;
            let hrefs = output::edge_hrefs(edge, &drawing.details);
            let anchor = edge.pos.as_ref().map_or([0.0; 2], |p| [p.x, p.y]);
            let (visible, parts) = Self::edge_paths(&laid, edge, carrier, anchor)?;
            let part_refs: Vec<&BezPath> = parts.iter().collect();
            for (region, points) in curves::region_samples(&part_refs, UNIT)? {
                targets.extend(points.into_iter().map(|at| Target {
                    at,
                    size: 8.0,
                    href: hrefs.regions[region].clone(),
                }));
            }
            Self::paint_edge(drawing, &visible, &parts, &mut layers, &mut strokes)?;

            // The bend relative to the chord sets the preferred momentum side.
            let side = match (&edge.source, &edge.sink) {
                (Some(source), Some(sink)) => {
                    let (a, b) = (laid.node(source.node), laid.node(sink.node));
                    let cross =
                        (b[0] - a[0]) * (anchor[1] - a[1]) - (b[1] - a[1]) * (anchor[0] - a[0]);
                    if cross.abs() > 1e-9 && cross < 0.0 {
                        -1.0
                    } else {
                        1.0
                    }
                }
                _ => 1.0,
            };
            let label = drawing
                .label
                .map(|page| {
                    typeset.pages[page]
                        .metrics
                        .map(|metrics| (page, metrics))
                        .ok_or("label page has no measurements")
                })
                .transpose()?;
            let paired = edge.source.is_some() && edge.sink.is_some();
            if drawing.momentum && !(paired && label.is_some()) {
                // External arrows are fixed decorations, as are internal ones
                // without a label to carry them.
                let arrow = curves::layer(
                    &visible,
                    labels::ARROW_OFFSET * side,
                    Some(labels::ARROW_WINDOW),
                    0.0,
                )?;
                fixed_arrows.extend(labels::arrow_footprint(&arrow));
                Self::paint_arrow(&arrow, &mut layers);
            }
            let Some((page, metrics)) = label else {
                continue;
            };
            let placement = |candidates, arrow_carriers| Placement {
                edge: edge.edge,
                page,
                candidates,
                href: hrefs.label.clone(),
                arrow_carriers,
            };
            match (&edge.source, &edge.sink) {
                (Some(_), Some(_)) if drawing.momentum => {
                    let (candidates, carriers) =
                        labels::momentum_candidates(&visible, metrics, side)?;
                    placements.push(placement(candidates, Some(carriers)));
                    region_hrefs.push((edge.edge, hrefs.regions.clone()));
                }
                (Some(_), Some(_)) => {
                    // Clear the painted band and at least the label gap.
                    let band = 0.06 + drawing.pattern.as_ref().map_or(0.0, Pattern::amplitude);
                    let candidates = labels::carrier_candidates(&visible, metrics, band.max(0.15))?;
                    placements.push(placement(candidates, None));
                }
                (source, sink) => {
                    let owner = source
                        .as_ref()
                        .or(sink.as_ref())
                        .ok_or("edge without endpoints")?;
                    let node = laid.node(owner.node);
                    let segments = curves::cubics(&visible);
                    let (dx, dy) = (anchor[0] - node[0], anchor[1] - node[1]);
                    let mut direction = match (source.is_some(), segments.first(), segments.last())
                    {
                        (true, _, Some(last)) => curves::cubic_tangent(last, 1.0),
                        (false, Some(first), _) => curves::cubic_tangent(first, 0.0).map(|v| -v),
                        _ => [dx, dy],
                    };
                    if direction[0].hypot(direction[1]) <= 1e-9 {
                        direction = [dx, dy];
                    }
                    if direction[0].hypot(direction[1]) <= 1e-9 {
                        direction = [1.0, 0.0];
                    }
                    let label_pos = edge
                        .label_pos
                        .as_ref()
                        .or(edge.pos.as_ref())
                        .map_or(anchor, |p| [p.x, p.y]);
                    let candidates = labels::endpoint_candidates(
                        label_pos,
                        direction,
                        Anchor::outward(dx, dy),
                        metrics,
                    )?;
                    placements.push(placement(candidates, None));
                }
            }
        }

        let obstacles: Vec<Bounds> = laid
            .nodes
            .iter()
            .map(|node| {
                let [x, y] = laid.node(node.node);
                let size = |key: &str| {
                    node.statements
                        .get(key)
                        .and_then(|value| value.parse::<f64>().ok())
                        .unwrap_or(2.0 * NODE_RADIUS)
                };
                let (w, h) = (size("layout-width") / 2.0, size("layout-height") / 2.0);
                Bounds {
                    left: x - w,
                    right: x + w,
                    bottom: y - h,
                    top: y + h,
                }
            })
            .collect();
        let searched = labels::search(&placements, &obstacles, &strokes, &fixed_arrows)?;
        let mut chosen = Vec::new();
        for (placement, &choice) in placements.iter().zip(&searched.choices) {
            let candidate = &placement.candidates[choice];
            if let Some(carriers) = &placement.arrow_carriers {
                let arrow = curves::layer(
                    &carriers[candidate.path_index],
                    0.0,
                    Some(labels::ARROW_WINDOW),
                    candidate.path_shift,
                )?;
                Self::paint_arrow(&arrow, &mut layers);
                let regions = &region_hrefs
                    .iter()
                    .find(|(edge, _)| *edge == placement.edge)
                    .ok_or("momentum arrow without hrefs")?
                    .1;
                for (region, points) in curves::region_samples(&[&arrow], UNIT)? {
                    targets.extend(points.into_iter().map(|at| Target {
                        at,
                        size: 8.0,
                        href: regions[region].clone(),
                    }));
                }
            }
            layers.push(Element::Label {
                page: placement.page,
                bounds: candidate.bounds,
                href: placement.href.clone(),
            });
            chosen.push(candidate.bounds);
        }
        for (node, drawing) in laid.nodes.iter().zip(&scene.nodes) {
            let at = laid.node(node.node);
            layers.push(Element::Node {
                at,
                stroke: drawing.stroke.clone(),
            });
            targets.push(Target {
                at,
                size: 10.0,
                href: output::node_href(node, &laid.edges, &drawing.details),
            });
        }
        let overlaps = labels::overlaps(&placements, &chosen, &searched.lines);
        Ok(Self {
            svg: output::svg(typeset, &layers, &targets),
            placements,
            chosen,
            overlaps,
        })
    }

    /// The visible (node-trimmed) edge path, and its source and sink parts.
    fn edge_paths(
        laid: &Laid,
        edge: &TypstDotEdge,
        carrier: &[[f64; 2]],
        anchor: [f64; 2],
    ) -> Result<(BezPath, Vec<BezPath>), String> {
        let interior = carrier
            .get(1..carrier.len().saturating_sub(1))
            .unwrap_or(&[]);
        let dangling = |points: Vec<[f64; 2]>, start: f64, end: f64| -> Result<_, String> {
            let path = curves::trim_routed(&curves::routed_curve(&points)?, start, end)?;
            Ok((path.clone(), vec![path]))
        };
        match (&edge.source, &edge.sink) {
            (Some(_), Some(_)) => {
                let curve = curves::routed_curve(carrier)?;
                let half = curves::length(&curve) / 2.0;
                let source =
                    curves::trim_routed(&curves::trim(&curve, 0.0, half)?, NODE_RADIUS, 0.0)?;
                let sink =
                    curves::trim_routed(&curves::trim(&curve, half, 0.0)?, 0.0, NODE_RADIUS)?;
                // `layer(curve, outsets)`: one trimmed window, reassembled from its cubics.
                let visible = curves::from_cubics(&curves::cubics(
                    &curves::windows(&curve, &[(NODE_RADIUS, NODE_RADIUS)])?.remove(0),
                ));
                Ok((visible, vec![source, sink]))
            }
            (Some(TypstDotEndpoint { node, .. }), None) => {
                let mut points = vec![laid.node(*node)];
                points.extend_from_slice(interior);
                points.push(anchor);
                dangling(points, NODE_RADIUS, 0.0)
            }
            (None, Some(TypstDotEndpoint { node, .. })) => {
                let mut points = vec![anchor];
                points.extend_from_slice(interior);
                points.push(laid.node(*node));
                dangling(points, 0.0, NODE_RADIUS)
            }
            (None, None) => Err("edge without endpoints".to_owned()),
        }
    }

    /// The base layer: a decorated edge, or each cubic as its own stroke, then
    /// the particle-flow arrowhead at the source/sink split or the center.
    fn paint_edge(
        drawing: &EdgeDrawing,
        visible: &BezPath,
        parts: &[BezPath],
        layers: &mut Vec<Element>,
        strokes: &mut Vec<CollisionStroke>,
    ) -> Result<(), String> {
        let radius = drawing.stroke.width / UNIT / 2.0;
        if let Some(pattern) = &drawing.pattern {
            let path = curves::pattern(visible, pattern)?;
            strokes.push(CollisionStroke {
                radius,
                path: path.clone(),
            });
            layers.push(Element::Path {
                path,
                stroke: drawing.stroke.clone(),
            });
        } else {
            for segment in curves::cubics(visible) {
                let path = curves::from_cubics(&[segment]);
                strokes.push(CollisionStroke {
                    radius,
                    path: path.clone(),
                });
                layers.push(Element::Path {
                    path,
                    stroke: drawing.stroke.clone(),
                });
            }
        }
        if let Some(forward) = drawing.flow {
            let ratio = match parts {
                [source, _] => {
                    let total = curves::length(visible);
                    if total <= curves::ACCURACY {
                        0.5
                    } else {
                        curves::length(source) / total
                    }
                }
                _ => 0.5,
            };
            if let Some(points) = marks::TRIANGLE.centered(visible, ratio, forward) {
                strokes.push(CollisionStroke {
                    radius: output::MARK_STROKE / UNIT / 2.0,
                    path: BezPath::from_vec(vec![
                        kurbo::PathEl::MoveTo(points[0]),
                        kurbo::PathEl::LineTo(points[1]),
                        kurbo::PathEl::LineTo(points[2]),
                        kurbo::PathEl::ClosePath,
                    ]),
                });
                layers.push(Element::Triangle(points));
            }
        }
        Ok(())
    }

    fn paint_arrow(arrow: &BezPath, layers: &mut Vec<Element>) {
        layers.push(Element::Path {
            path: arrow.clone(),
            stroke: arrow_stroke(),
        });
        if let Some(points) = marks::CHEVRON.at_end(arrow) {
            layers.push(Element::Chevron(points));
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;
    use crate::{TypstEdgeSpec, TypstEndpointSpec, TypstNodeSpec};

    fn ink(width: f64) -> Stroke {
        Stroke {
            paint: output::INK.to_owned(),
            width,
            dash: Dash::Solid,
            round_cap: true,
        }
    }

    /// A gluon and a fermion between two vertices, with a wave and a zigzag leg.
    fn bubble(momentum: bool) -> Scene {
        let end = |node| {
            Some(TypstEndpointSpec {
                node,
                statement: None,
                id: None,
                data: None,
                port_label: None,
                compass: None,
                in_subgraph: false,
                route_points: Vec::new(),
            })
        };
        let leg = |side: &str, x: i32| -> BTreeMap<String, String> {
            [
                ("group-start-x", "true".to_owned()),
                ("pin", format!("x:@{side}")),
                ("pos", format!("{x},0")),
                ("pos-x-set", "true".to_owned()),
                ("pos-y-set", "true".to_owned()),
                ("pos-z", "0".to_owned()),
                ("pos-z-mode", "pin".to_owned()),
            ]
            .map(|(key, value)| (key.to_owned(), value))
            .into()
        };
        let edge = |id, source, sink, statements| TypstEdgeSpec {
            name: None,
            source,
            sink,
            data: None,
            orientation: Some("default".to_owned()),
            flow: None,
            id: Some(id),
            pos: None,
            statements,
        };
        let node = |index: usize| TypstNodeSpec {
            name: Some(format!("v{index}")),
            index: Some(index),
            data: None,
            pos: None,
            statements: [("layout-width", "0.16"), ("layout-height", "0.16")]
                .map(|(key, value)| (key.to_owned(), value.to_owned()))
                .into(),
        };
        let drawing = |pattern, flow| EdgeDrawing {
            stroke: ink(1.0),
            pattern,
            flow,
            momentum,
            label: Some(0),
            details: [("particle", "x")].into_iter().collect(),
        };
        let coil = Pattern::Coil {
            amplitude: 0.15,
            wavelength: 0.45,
            longitudinal_scale: 1.4,
        };
        Scene {
            graph: TypstGraphSpec {
                name: Some("bubble".to_owned()),
                data: None,
                statements: BTreeMap::new(),
                default_edge_statements: BTreeMap::new(),
                default_node_statements: BTreeMap::new(),
                nodes: vec![node(0), node(1)],
                edges: vec![
                    edge(0, None, end(0), leg("-left", -10)),
                    edge(1, end(0), end(1), BTreeMap::new()),
                    edge(2, end(0), end(1), BTreeMap::new()),
                    edge(3, end(1), None, leg("+right", 10)),
                ],
            },
            nodes: (0..2)
                .map(|_| NodeDrawing {
                    stroke: ink(1.45),
                    details: Details::default(),
                })
                .collect(),
            edges: vec![
                drawing(
                    Some(Pattern::Wave {
                        amplitude: 0.14,
                        wavelength: 0.55,
                    }),
                    None,
                ),
                drawing(Some(coil), None),
                drawing(None, Some(true)),
                drawing(
                    Some(Pattern::Zigzag {
                        amplitude: 0.14,
                        wavelength: 0.55,
                    }),
                    None,
                ),
            ],
            preamble: String::new(),
            title: Some("#\"bubble\"".to_owned()),
            pages: vec!["[x]".to_owned()],
            layout: serde_json::Map::new(),
            label_feedback: false,
        }
    }

    fn page(body: &str) -> String {
        format!(
            r#"<svg class="typst-doc" viewBox="0 0 10 12" width="10pt" height="12pt" xmlns="http://www.w3.org/2000/svg">{body}<defs id="glyph"><symbol id="g1" overflow="visible"><path d="M0 0"/></symbol></defs></svg>"#
        )
    }

    fn typeset(scene: &Scene) -> Typeset {
        let label = page(r#"<g/><a href="linnest-metrics:10,8,8,12"><rect/></a>"#);
        scene.typeset(&[page("<g/>"), label]).unwrap()
    }

    #[test]
    fn reads_title_and_measured_label_pages() {
        let scene = bubble(false);
        let typeset = typeset(&scene);
        assert!(
            typeset
                .title
                .as_ref()
                .is_some_and(|title| title.metrics.is_none())
        );
        assert_eq!(typeset.pages[0].metrics, Some([10.0, 8.0, 8.0, 12.0]));
        assert_eq!(typeset.pages[0].body, "<g/>");
        assert_eq!(typeset.defs.matches("<symbol").count(), 1);
        assert!(scene.typeset(&[page("<g/>")]).is_err());
        let document = scene.label_document();
        assert!(document.starts_with("#set page(width: auto, height: auto, margin: 0pt"));
        assert!(document.contains("#box[#\"bubble\"]\n#pagebreak()\n#context"));
    }

    #[test]
    fn draws_decorations_flow_marks_labels_and_links() {
        let scene = bubble(false);
        let svg = scene.render(&typeset(&scene)).unwrap();
        assert!(svg.starts_with("<svg viewBox=\"0 0 "));
        assert_eq!(svg.matches("<circle ").count(), 2);
        // One fermion arrowhead; every decoration is a single path.
        assert_eq!(svg.matches(r##"<path fill="#3d2645""##).count(), 1);
        // Every edge is labelled and linked, with its half-edges.
        for edge in 0..4 {
            assert!(svg.contains(&format!("#linnet-edge-{edge}?")));
        }
        assert!(svg.contains("#linnet-halfedge-1?"));
        assert!(svg.contains("#linnet-node-1?"));
        assert!(svg.contains("&quot;particle&quot;:&quot;x&quot;"));
    }

    #[test]
    fn rides_momentum_labels_on_their_arrows() {
        let scene = bubble(true);
        let svg = scene.render(&typeset(&scene)).unwrap();
        // One chevron per edge: searched internal arrows and fixed external ones.
        let chevron = r##"fill="none" stroke="#3d2645" stroke-width="1" stroke-linecap="round" stroke-linejoin="miter""##;
        assert_eq!(svg.matches(chevron).count(), 4);
        // Edges without labels keep their arrows as fixed decorations.
        let mut unlabelled = bubble(true);
        for edge in &mut unlabelled.edges {
            edge.label = None;
        }
        let svg = unlabelled.render(&typeset(&unlabelled)).unwrap();
        assert_eq!(svg.matches(chevron).count(), 4);
    }
}
