//! Native SVG drawing of prepared graphs, the host-side counterpart of
//! `draw.typ`: native layout, Kurvst geometry, the shared annotation search
//! and SVG output. Typst only typesets the label pages a [`Scene`] lists.
mod config;
mod curves;
pub use config::Config;
mod interactive;
mod labels;
mod marks;
mod output;

use kurbo::{BezPath, Point};
use linnet::half_edge::{
    involution::EdgeIndex,
    layout::impred::{EdgeLabel, ImpredConfig},
    subgraph::{ModifySubSet, SuBitGraph},
};
use serde_json::Value;
use std::sync::Arc;

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
    /// Layout options; `layout-algo` selects `dot` or `impred` (the default).
    pub layout: serde_json::Map<String, Value>,
    /// Edges defining the layered hierarchy, in spec order. Other edges are
    /// routed around it without changing the ranks. `None` uses all edges.
    pub layout_edges: Option<Vec<usize>>,
    /// Refine the layout around drawn labels in short warm ImPrEd passes.
    pub label_feedback: bool,
}

#[derive(Clone)]
pub struct NodeDrawing {
    /// Radius/minimum half-size in drawing units.
    pub radius: f64,
    /// Fit labels with this padding, or keep an explicit fixed radius.
    pub label_padding: Option<f64>,
    pub label: Option<usize>,
    pub rectangular: bool,
    pub fill: String,
    pub stroke: Stroke,
    /// Inspection fields after the node's incident edges.
    pub details: Details,
}

#[derive(Clone)]
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
        let mut spec = self.graph.clone();
        for (node, drawing) in spec.nodes.iter_mut().zip(&self.nodes) {
            let (w, h) = drawing.size(typeset);
            node.statements
                .insert("layout-width".into(), (2.0 * w).to_string());
            node.statements
                .insert("layout-height".into(), (2.0 * h).to_string());
        }
        let graph = TypstGraph::from_spec(spec)?;
        match self.layout.get("layout-algo").and_then(Value::as_str) {
            Some("dot") => {
                let laid = Laid::layered(graph, self)?;
                return Ok(Drawing::new(self, typeset, &laid)?.svg(typeset));
            }
            Some("impred") | None => {}
            Some(algo) => return Err(format!("unsupported SVG layout algorithm {algo:?}")),
        }
        let mut run = ImpredRun::seeded(graph, &Value::Object(self.layout.clone()))?;
        let drawing = Drawing::new(self, typeset, &Laid::new(&run)?)?;
        if !self.label_feedback {
            return Ok(drawing.svg(typeset));
        }
        // Keep the least-overlapping drawing; later passes may not improve it.
        // A rejected pass still supplies the next pass's pinned labels.
        let mut best = drawing;
        let mut rejected = None;
        for _ in 0..FEEDBACK_PASSES {
            let drawing = rejected.as_ref().unwrap_or(&best);
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
            let drawing = Drawing::new(self, typeset, &Laid::new(&run)?)?;
            if drawing.score() < best.score() {
                best = drawing;
                rejected = None;
            } else {
                rejected = Some(drawing);
            }
        }
        Ok(best.svg(typeset))
    }
}

/// The laid-out graph: node positions, edge records and complete carriers.
struct Laid {
    nodes: Vec<TypstDotNode>,
    edges: Vec<TypstDotEdge>,
    carriers: Vec<Vec<[f64; 2]>>,
}

impl Laid {
    fn layered(mut graph: TypstGraph, scene: &Scene) -> Result<Self, String> {
        graph.layout_config = serde_json::from_value(Value::Object(scene.layout.clone()))
            .map_err(|error| format!("invalid layered layout options: {error}"))?;
        let selected = scene
            .layout_edges
            .as_ref()
            .map(|edges| {
                let mut selected = graph.empty_subgraph::<SuBitGraph>();
                for &edge in edges {
                    if edge >= graph.n_edges() {
                        return Err("layout edge is outside the scene graph".to_owned());
                    }
                    selected.add(graph.graph[&EdgeIndex(edge)].1);
                }
                Ok(selected)
            })
            .transpose()?;
        graph.layout_with_subgraph(selected.as_ref())?;
        let nodes = graph.node_records()?;
        let edges = graph.edge_records()?;
        let point = |p: &crate::TypstPoint| [p.x, p.y];
        let carriers = edges
            .iter()
            .map(|edge| {
                let mut route = Vec::new();
                if let Some(source) = &edge.source {
                    route.push(point(nodes[source.node].pos.as_ref().unwrap()));
                    route.extend(source.route_points.iter().map(point));
                }
                route.push(point(edge.pos.as_ref().unwrap()));
                if let Some(sink) = &edge.sink {
                    route.extend(sink.route_points.iter().rev().map(point));
                    route.push(point(nodes[sink.node].pos.as_ref().unwrap()));
                }
                route
            })
            .collect();
        Ok(Self {
            nodes,
            edges,
            carriers,
        })
    }

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
    /// Invisible structural paths for overlays, independent of particle decoration.
    Carrier {
        path: BezPath,
        href: Arc<str>,
    },
    Path {
        path: BezPath,
        stroke: Stroke,
    },
    Chevron([Point; 3]),
    Triangle([Point; 3]),
    Label {
        page: usize,
        bounds: Bounds,
        href: Option<String>,
    },
    Node {
        id: usize,
        at: [f64; 2],
        size: (f64, f64),
        rectangular: bool,
        fill: String,
        stroke: Stroke,
    },
}

/// A transparent hover target.
struct Target {
    at: [f64; 2],
    size: f64,
    href: Arc<str>,
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
    layers: Vec<Element>,
    targets: Vec<Target>,
    placements: Vec<Placement>,
    chosen: Vec<Bounds>,
    /// Overlapping labels, line hits, label pairs.
    overlaps: (usize, usize, usize),
}

impl Drawing {
    /// Serialize only the selected feedback pass, preserving painting order.
    fn svg(&self, typeset: &Typeset) -> String {
        output::svg(typeset, &self.layers, &self.targets)
    }

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

    fn new(scene: &Scene, typeset: &Typeset, laid: &Laid) -> Result<Self, String> {
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
            let (visible, parts) = Self::edge_paths(scene, typeset, laid, edge, carrier, anchor)?;
            let part_refs: Vec<&BezPath> = parts.iter().collect();
            for (region, points) in curves::region_samples(&part_refs, UNIT)? {
                targets.extend(points.into_iter().map(|at| Target {
                    at,
                    size: 8.0,
                    href: hrefs.regions[region].clone(),
                }));
            }
            Self::paint_edge(drawing, &visible, &parts, &mut layers, &mut strokes)?;
            for (index, path) in parts.iter().enumerate() {
                layers.push(Element::Carrier {
                    path: path.clone(),
                    href: hrefs.regions[if index == 0 { 0 } else { 3 }].clone(),
                });
            }

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
                href: Some(placement.href.clone()),
            });
            chosen.push(candidate.bounds);
        }
        for (node, drawing) in laid.nodes.iter().zip(&scene.nodes) {
            let at = laid.node(node.node);
            let href = output::node_href(node, &laid.edges, &drawing.details);
            layers.push(Element::Node {
                id: drawing
                    .details
                    .get("node")
                    .and_then(Value::as_u64)
                    .map_or(node.node, |id| id as usize),
                at,
                size: drawing.size(typeset),
                rectangular: drawing.rectangular,
                fill: drawing.fill.clone(),
                stroke: drawing.stroke.clone(),
            });
            if let Some(page) = drawing.label {
                let label = &typeset.pages[page];
                layers.push(Element::Label {
                    page,
                    bounds: Bounds {
                        left: at[0] - label.width / UNIT / 2.0,
                        right: at[0] + label.width / UNIT / 2.0,
                        bottom: at[1] - label.height / UNIT / 2.0,
                        top: at[1] + label.height / UNIT / 2.0,
                    },
                    href: None,
                });
            }
            // One target covers both the node and its label. Separate
            // overlapping anchors steal clicks from the keyboard-focusable one.
            let (width, height) = drawing.size(typeset);
            let label_size = drawing.label.map_or(0.0, |page| {
                typeset.pages[page].width.max(typeset.pages[page].height)
            });
            targets.push(Target {
                at,
                size: (width.max(height) * UNIT * 2.0).max(label_size).max(10.0),
                href: href.into(),
            });
        }
        let overlaps = labels::overlaps(&placements, &chosen, &searched.lines);
        Ok(Self {
            layers,
            targets,
            placements,
            chosen,
            overlaps,
        })
    }

    /// The visible (node-trimmed) edge path, and its source and sink parts.
    fn edge_paths(
        scene: &Scene,
        typeset: &Typeset,
        laid: &Laid,
        edge: &TypstDotEdge,
        carrier: &[[f64; 2]],
        anchor: [f64; 2],
    ) -> Result<(BezPath, Vec<BezPath>), String> {
        // Compass ports attach to the measured boundary. In particular, a
        // dependency can enter/leave vertically without following the chord
        // between its parent and child, as in the Typst anchored-edge drawing.
        let port = |end: &TypstDotEndpoint| {
            let direction = match end.compass.as_deref() {
                Some("n") => [0.0, 1.0],
                Some("s") => [0.0, -1.0],
                Some("e") => [1.0, 0.0],
                Some("w") => [-1.0, 0.0],
                _ => return None,
            };
            let center = laid.node(end.node);
            let (w, h) = scene.nodes[end.node].size(typeset);
            Some((
                [center[0] + direction[0] * w, center[1] + direction[1] * h],
                direction,
            ))
        };
        let outset = |end: &TypstDotEndpoint| {
            if port(end).is_some() {
                return 0.0;
            }
            let (w, h) = scene.nodes[end.node].size(typeset);
            if !scene.nodes[end.node].rectangular {
                return w;
            }
            let center = laid.node(end.node);
            let neighbor = if edge
                .source
                .as_ref()
                .is_some_and(|source| source.hedge == end.hedge)
            {
                carrier.iter().find(|point| **point != center)
            } else {
                carrier.iter().rev().find(|point| **point != center)
            };
            let Some(neighbor) = neighbor else { return 0.0 };
            let dx = (neighbor[0] - center[0]).abs();
            let dy = (neighbor[1] - center[1]).abs();
            (w / dx).min(h / dy) * dx.hypot(dy)
        };
        let anchored = edge
            .source
            .iter()
            .chain(&edge.sink)
            .any(|end| port(end).is_some());
        if anchored {
            let endpoint = |end: &Option<TypstDotEndpoint>| {
                end.as_ref().map_or((anchor, None), |end| {
                    port(end).map_or((laid.node(end.node), None), |(point, dir)| {
                        (point, Some(dir))
                    })
                })
            };
            let (start, source_dir) = endpoint(&edge.source);
            let (end, sink_dir) = endpoint(&edge.sink);
            let paired = edge.source.is_some() && edge.sink.is_some();
            let direct = edge
                .statements
                .get("route")
                .is_some_and(|route| route == "direct");
            let amount = if direct {
                (0.18 * (end[0] - start[0]).abs() + 0.3 * (end[1] - start[1]).abs())
                    .clamp(0.45, 4.0)
            } else {
                0.01
            };
            let guide = |point: [f64; 2], dir: Option<[f64; 2]>| {
                Point::new(
                    point[0] + dir.map_or((anchor[0] - point[0]) / 3.0, |d| d[0] * amount),
                    point[1] + dir.map_or((anchor[1] - point[1]) / 3.0, |d| d[1] * amount),
                )
            };
            let (a, b, middle) = (
                Point::new(start[0], start[1]),
                Point::new(end[0], end[1]),
                Point::new(anchor[0], anchor[1]),
            );
            let curve = if !paired {
                let mut points = vec![start];
                points.extend_from_slice(
                    carrier
                        .get(1..carrier.len().saturating_sub(1))
                        .unwrap_or(&[]),
                );
                points.push(end);
                curves::routed_curve(&points)?
            } else if direct {
                let chord = b - a;
                let handle = if chord.hypot() == 0.0 {
                    chord
                } else {
                    chord / chord.hypot()
                        * amount.min(a.distance(middle).min(b.distance(middle)) / 3.0)
                };
                curves::from_cubics(&[
                    kurbo::CubicBez::new(a, guide(start, source_dir), middle - handle, middle),
                    kurbo::CubicBez::new(middle, middle + handle, guide(end, sink_dir), b),
                ])
            } else {
                let mut points = vec![start];
                if source_dir.is_some() {
                    let p = guide(start, source_dir);
                    points.push([p.x, p.y]);
                }
                points.extend_from_slice(
                    carrier
                        .get(1..carrier.len().saturating_sub(1))
                        .unwrap_or(&[]),
                );
                if sink_dir.is_some() {
                    let p = guide(end, sink_dir);
                    points.push([p.x, p.y]);
                }
                points.push(end);
                curves::routed_curve(&points)?
            };
            let start_trim = edge.source.as_ref().map_or(0.0, outset);
            let end_trim = edge.sink.as_ref().map_or(0.0, outset);
            return Self::split_edge(&curve, start_trim, end_trim, paired);
        }
        let interior = carrier
            .get(1..carrier.len().saturating_sub(1))
            .unwrap_or(&[]);
        let dangling = |points: Vec<[f64; 2]>, start: f64, end: f64| -> Result<_, String> {
            let path = curves::trim_routed(&curves::routed_curve(&points)?, start, end)?;
            Ok((path.clone(), vec![path]))
        };
        match (&edge.source, &edge.sink) {
            (Some(source_node), Some(sink_node)) => {
                let source_radius = outset(source_node);
                let sink_radius = outset(sink_node);
                let curve = curves::routed_curve(carrier)?;
                Self::split_edge(&curve, source_radius, sink_radius, true)
            }
            (Some(TypstDotEndpoint { node, .. }), None) => {
                let mut points = vec![laid.node(*node)];
                points.extend_from_slice(interior);
                points.push(anchor);
                dangling(points, outset(edge.source.as_ref().unwrap()), 0.0)
            }
            (None, Some(TypstDotEndpoint { node, .. })) => {
                let mut points = vec![anchor];
                points.extend_from_slice(interior);
                points.push(laid.node(*node));
                dangling(points, 0.0, outset(edge.sink.as_ref().unwrap()))
            }
            (None, None) => Err("edge without endpoints".to_owned()),
        }
    }

    fn split_edge(
        curve: &BezPath,
        start: f64,
        end: f64,
        paired: bool,
    ) -> Result<(BezPath, Vec<BezPath>), String> {
        if !paired {
            let path = curves::trim_routed(curve, start, end)?;
            return Ok((path.clone(), vec![path]));
        }
        let half = curves::length(curve) / 2.0;
        let source = curves::trim_routed(&curves::trim(curve, 0.0, half)?, start, 0.0)?;
        let sink = curves::trim_routed(&curves::trim(curve, half, 0.0)?, 0.0, end)?;
        let visible = curves::from_cubics(&curves::cubics(
            &curves::windows(curve, &[(start, end)])?.remove(0),
        ));
        Ok((visible, vec![source, sink]))
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

impl NodeDrawing {
    fn size(&self, typeset: &Typeset) -> (f64, f64) {
        let Some(padding) = self.label_padding else {
            return (self.radius, self.radius);
        };
        let Some(label) = self.label.map(|i| &typeset.pages[i]) else {
            return (self.radius, self.radius);
        };
        let (w, h) = (
            (label.width / UNIT / 2.0 + padding).max(self.radius),
            (label.height / UNIT / 2.0 + padding).max(self.radius),
        );
        if self.rectangular {
            (w, h)
        } else {
            (w.hypot(h), w.hypot(h))
        }
    }
}

impl Scene {
    /// Arrange independently rendered channels in one valid SVG document.
    pub fn combine_svgs(figures: &[String]) -> Result<String, String> {
        if let [figure] = figures {
            return Ok(figure.clone());
        }
        let mut body = String::new();
        let (mut width, mut height) = (0.0_f64, 0.0_f64);
        for figure in figures {
            let document = roxmltree::Document::parse(figure).map_err(|e| e.to_string())?;
            let size = document
                .root_element()
                .attribute("viewBox")
                .ok_or("SVG has no viewBox")?
                .split_whitespace()
                .map(str::parse::<f64>)
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| e.to_string())?;
            if size.len() != 4 {
                return Err("invalid SVG viewBox".into());
            }
            // Nested viewports use the outer SVG's graph units, not CSS points.
            let end = figure.find('>').ok_or("SVG has no opening tag")?;
            let mut opening = figure[..end].to_owned();
            let mut dimensions = document
                .root_element()
                .attributes()
                .filter(|attribute| matches!(attribute.name(), "width" | "height"))
                .map(|attribute| attribute.range())
                .collect::<Vec<_>>();
            dimensions.sort_by_key(|range| range.start);
            for range in dimensions.into_iter().rev() {
                opening.replace_range(range, "");
            }
            body.push_str(&format!(
                "{opening} x=\"{width}\" y=\"0\" width=\"{}\" height=\"{}\"{}",
                size[2],
                size[3],
                &figure[end..]
            ));
            width += size[2] + 10.0;
            height = height.max(size[3]);
        }
        Ok(format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width} {height}\" width=\"{width}pt\" height=\"{height}pt\">{body}</svg>"
        ))
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
                    radius: NODE_RADIUS,
                    label_padding: Some(0.25),
                    label: None,
                    rectangular: false,
                    fill: "none".into(),
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
            layout_edges: None,
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
        assert!(typeset
            .title
            .as_ref()
            .is_some_and(|title| title.metrics.is_none()));
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
    fn node_labels_share_one_click_target_with_their_node() {
        let mut scene = bubble(false);
        scene.nodes[0].label = Some(0);
        let svg = scene.render(&typeset(&scene)).unwrap();
        assert_eq!(svg.matches("#linnet-node-0?").count(), 1);
        assert_eq!(svg.matches("#linnet-node-1?").count(), 1);
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

    #[test]
    fn layered_scene_keeps_dependency_ranks_and_routes_contractions_below_leaves() {
        let mut scene = bubble(false);
        scene.title = None;
        scene.pages.clear();
        let mut node = scene.graph.nodes[0].clone();
        node.statements.clear();
        scene.graph.nodes = (0..4)
            .map(|index| TypstNodeSpec {
                index: Some(index),
                name: Some(format!("n{index}")),
                ..node.clone()
            })
            .collect();
        let mut drawing = scene.nodes[0].clone();
        drawing.rectangular = true;
        drawing.radius = 0.5;
        scene.nodes = vec![drawing; 4];
        // The root deliberately is not node zero. A contraction between the
        // leaves must not put either leaf below its sibling.
        scene.layout = serde_json::from_value(serde_json::json!({
            "layout-algo": "dot", "layout-roots": [3], "tree-dx": 0.35, "tree-dy": 2.2,
        }))
        .unwrap();
        scene.layout_edges = Some(vec![0, 1, 2, 3]);
        let endpoint = |node, compass: &str| {
            Some(TypstEndpointSpec {
                node,
                compass: Some(compass.into()),
                id: None,
                statement: None,
                data: None,
                port_label: None,
                in_subgraph: false,
                route_points: vec![],
            })
        };
        let template = scene.graph.edges[1].clone();
        scene.graph.edges = [
            (endpoint(0, "n"), endpoint(3, "s")),
            (endpoint(1, "n"), endpoint(0, "s")),
            (endpoint(2, "n"), endpoint(0, "s")),
            (endpoint(3, "n"), None),
            (endpoint(1, "s"), endpoint(2, "s")),
            (endpoint(2, "s"), None),
        ]
        .into_iter()
        .enumerate()
        .map(|(index, (source, sink))| TypstEdgeSpec {
            id: Some(index),
            source,
            sink,
            statements: [(
                "route".into(),
                if index < 4 { "direct" } else { "hobby-through" }.into(),
            )]
            .into(),
            ..template.clone()
        })
        .collect();
        let mut edge = scene.edges[0].clone();
        edge.pattern = None;
        edge.label = None;
        scene.edges = vec![edge; 6];
        let typeset = scene.typeset(&[]).unwrap();
        let mut spec = scene.graph.clone();
        for node in &mut spec.nodes {
            node.statements.extend([
                ("layout-width".into(), "1".into()),
                ("layout-height".into(), "1".into()),
            ]);
        }
        let laid = Laid::layered(TypstGraph::from_spec(spec).unwrap(), &scene).unwrap();
        assert!(laid.node(3)[1] > laid.node(0)[1]);
        assert!(laid.node(0)[1] > laid.node(1)[1]);
        assert_eq!(laid.node(1)[1], laid.node(2)[1]);
        assert_ne!(laid.node(1)[0], laid.node(2)[0]);
        for (index, record) in laid.edges.iter().enumerate() {
            let pos = record.pos.as_ref().unwrap();
            let (path, _) = Drawing::edge_paths(
                &scene,
                &typeset,
                &laid,
                record,
                &laid.carriers[index],
                [pos.x, pos.y],
            )
            .unwrap();
            let cubics = curves::cubics(&path);
            assert!(!cubics.is_empty());
            let source = record.source.as_ref().unwrap();
            let center = laid.node(source.node);
            assert!((cubics[0].p0.x - center[0]).abs() < 1e-9);
            let sign = if index < 4 { 1.0 } else { -1.0 };
            assert!((cubics[0].p0.y - center[1] - sign * 0.5).abs() < 1e-9);
            if index < 3 {
                assert!((cubics[0].p1.x - cubics[0].p0.x).abs() < 1e-9);
                assert!(cubics[0].p1.y > cubics[0].p0.y);
            } else if index == 4 {
                assert!(cubics.iter().any(|c| c.p3.y < center[1] - 0.5));
            } else if index == 5 {
                // A free tensor slot must leave below its box, not double back
                // through it or turn a collinear Hobby guide into a huge loop.
                assert!(cubics.last().unwrap().p3.y < cubics[0].p0.y);
                assert!(curves::length(&path) < (laid.node(3)[1] - center[1]).abs());
            }
        }
        // Enabling force-label feedback must not relax away the selected ranks.
        let expected = scene.render(&typeset).unwrap();
        scene.label_feedback = true;
        assert_eq!(scene.render(&typeset).unwrap(), expected);
        assert!(expected.contains("#linnet-node-3?"));
    }

    #[test]
    fn feedback_emits_the_same_svg_as_serializing_every_pass() {
        // Keep an eager reference to check both painting order and selection:
        // later rejected/tied passes must still drive the next warm layout.
        fn eager(scene: &Scene, typeset: &Typeset) -> (String, Vec<(usize, usize, usize)>) {
            let mut spec = scene.graph.clone();
            for (node, drawing) in spec.nodes.iter_mut().zip(&scene.nodes) {
                let (w, h) = drawing.size(typeset);
                node.statements
                    .insert("layout-width".into(), (2.0 * w).to_string());
                node.statements
                    .insert("layout-height".into(), (2.0 * h).to_string());
            }
            let graph = TypstGraph::from_spec(spec).unwrap();
            let mut run = ImpredRun::seeded(graph, &Value::Object(scene.layout.clone())).unwrap();
            let mut drawing = Drawing::new(scene, typeset, &Laid::new(&run).unwrap()).unwrap();
            let mut scores = vec![drawing.score()];
            let (mut best, mut svg) = (drawing.score(), drawing.svg(typeset));
            for _ in 0..FEEDBACK_PASSES {
                if drawing.overlaps.0 == 0 {
                    break;
                }
                run.layout.labels = drawing.pinned_labels(&run);
                run.layout
                    .solve(ImpredConfig {
                        labels: true,
                        warm_start: FEEDBACK_WARM,
                        steps: FEEDBACK_STEPS,
                        ..run.config
                    })
                    .unwrap();
                drawing = Drawing::new(scene, typeset, &Laid::new(&run).unwrap()).unwrap();
                let rendered = drawing.svg(typeset);
                scores.push(drawing.score());
                if drawing.score() < best {
                    best = drawing.score();
                    svg = rendered;
                }
            }
            (svg, scores)
        }

        let mut continued_after_rejection = false;
        for momentum in [false, true] {
            for width in [10.0, 80.0] {
                let mut scene = bubble(momentum);
                scene.label_feedback = true;
                let mut typeset = typeset(&scene);
                typeset.pages[0].width = width;
                typeset.pages[0].metrics = Some([width, 8.0, 8.0, 12.0]);
                let (expected, scores) = eager(&scene, &typeset);
                continued_after_rejection |= (1..scores.len().saturating_sub(1))
                    .any(|i| scores[i] >= *scores[..i].iter().min().unwrap());
                assert_eq!(scene.render(&typeset).unwrap(), expected, "{scores:?}");
            }
        }
        assert!(
            continued_after_rejection,
            "exercise another feedback pass after rejecting a drawing"
        );
    }
}
