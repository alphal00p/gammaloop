//! Native SVG drawing of prepared graphs, the host-side counterpart of
//! `draw.typ`: ImPrEd layout, Kurvst geometry, the shared annotation search
//! and SVG output. Typst only typesets the label pages a [`Scene`] lists.
mod config;
mod curves;
pub use config::Config;
mod interactive;
mod labels;
mod marks;
mod output;

use kurbo::BezPath;
use linnet::half_edge::layout::impred::{EdgeLabel, ImpredConfig};
pub use marks::MarkPaint;
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

#[derive(Clone)]
pub struct NodeDrawing {
    /// Radius/minimum half-size in drawing units.
    pub radius: f64,
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
    pub flow_arrow: kurvst::marks::MarkSpec,
    pub momentum_arrow: kurvst::marks::MarkSpec,
    pub flow_arrow_paints: Vec<MarkPaint>,
    pub momentum_arrow_paints: Vec<MarkPaint>,
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
    Mark {
        drawable: kurvst::marks::MarkDrawable,
        paint: MarkPaint,
        width: f64,
    },
    HitBox {
        bounds: Bounds,
        href: String,
    },
    Label {
        page: usize,
        bounds: Bounds,
        href: String,
    },
    Node {
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

/// Stable placement indices shared by the two drawing-wide geometry batches.
struct EdgeGeometry {
    edge: usize,
    visible: BezPath,
    parts: Vec<BezPath>,
    hrefs: output::EdgeHrefs,
    flow: Option<usize>,
    fixed_arrow: Option<usize>,
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
        let mut batch = marks::Batch::new();
        let mut edge_geometry = Vec::new();
        let mut momentum = Vec::new();
        for (edge, carrier) in laid.edges.iter().zip(&laid.carriers) {
            let drawing = scene
                .edges
                .get(edge.edge)
                .ok_or("scene edges do not match the layout graph")?;
            let hrefs = output::edge_hrefs(edge, &drawing.details);
            let anchor = edge.pos.as_ref().map_or([0.0; 2], |p| [p.x, p.y]);
            let (visible, parts) = Self::edge_paths(scene, typeset, &laid, edge, carrier, anchor)?;
            let flow = if let Some(forward) = drawing.flow {
                let template = batch.template(&drawing.flow_arrow, &drawing.stroke)?;
                let ratio = match parts.as_slice() {
                    [source, _] => {
                        let total = curves::length(&visible);
                        if total <= curves::ACCURACY {
                            0.5
                        } else {
                            curves::length(source) / total
                        }
                    }
                    _ => 0.5,
                };
                Some(batch.centered(template, &visible, ratio, forward))
            } else {
                None
            };

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
            let mut fixed_arrow = None;
            if drawing.momentum && !(paired && label.is_some()) {
                // External arrows are fixed decorations, as are internal ones
                // without a label to carry them.
                let arrow = curves::layer(
                    &visible,
                    labels::ARROW_OFFSET * side,
                    Some(labels::ARROW_WINDOW),
                    0.0,
                )?;
                let template = batch.template(&drawing.momentum_arrow, &arrow_stroke())?;
                fixed_arrow =
                    Some(batch.push(template, &arrow, kurvst::marks::MarkStation::End, true));
            }
            edge_geometry.push(EdgeGeometry {
                edge: edge.edge,
                visible: visible.clone(),
                parts,
                hrefs: output::edge_hrefs(edge, &drawing.details),
                flow,
                fixed_arrow,
            });
            let Some((page, metrics)) = label else {
                continue;
            };
            let placement = |candidates| Placement {
                edge: edge.edge,
                page,
                candidates,
                href: hrefs.label.clone(),
            };
            match (&edge.source, &edge.sink) {
                (Some(_), Some(_)) if drawing.momentum => {
                    let template = batch.template(&drawing.momentum_arrow, &arrow_stroke())?;
                    let proposals = labels::MomentumCandidates::new(
                        &visible, metrics, side, template, &mut batch,
                    )?;
                    momentum.push((placements.len(), proposals));
                    placements.push(placement(Vec::new()));
                    region_hrefs.push((edge.edge, hrefs.regions.clone()));
                }
                (Some(_), Some(_)) => {
                    // Clear the painted band and at least the label gap.
                    let band = 0.06 + drawing.pattern.as_ref().map_or(0.0, Pattern::amplitude);
                    let candidates = labels::carrier_candidates(&visible, metrics, band.max(0.15))?;
                    placements.push(placement(candidates));
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
                    placements.push(placement(candidates));
                }
            }
        }

        // Every head, fixed decoration and searched proposal shares this one
        // candidate preparation pass. Do not call the engine inside edge loops.
        let candidates = batch.geometry(kurvst::marks::MarkGeometryMode::Candidates)?;
        for (index, proposals) in momentum {
            placements[index].candidates = proposals.resolve(&candidates.marks)?;
        }
        for edge in &edge_geometry {
            let drawing = &scene.edges[edge.edge];
            let mark = edge.flow.map(|index| &candidates.marks[index]);
            let shaft = mark.map_or(&edge.visible, |mark| &mark.shaft.path);
            Self::paint_edge(
                drawing,
                shaft,
                mark,
                &edge.hrefs.label,
                &mut Vec::new(),
                &mut strokes,
            )?;
            if let Some(index) = edge.fixed_arrow {
                fixed_arrows.push(Bounds::painted(&candidates.marks[index].footprint.path));
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
        let mut selected: Vec<_> = edge_geometry
            .iter()
            .flat_map(|edge| [edge.flow, edge.fixed_arrow].into_iter().flatten())
            .collect();
        selected.extend(
            placements
                .iter()
                .zip(&searched.choices)
                .filter_map(|(placement, &choice)| placement.candidates[choice].arrow_mark),
        );
        let geometry = batch.selected(&selected)?;
        let heads: std::collections::BTreeMap<_, _> =
            selected.into_iter().zip(geometry.marks).collect();
        let shafts: std::collections::BTreeMap<_, _> = geometry
            .shafts
            .into_iter()
            .map(|shaft| (shaft.carrier, shaft))
            .collect();
        for edge in &edge_geometry {
            let drawing = &scene.edges[edge.edge];
            let mark = edge.flow.map(|index| &heads[&index]);
            let shaft = mark.map_or(&edge.visible, |mark| &shafts[&mark.carrier].shaft.path);
            Self::paint_edge(
                drawing,
                shaft,
                mark,
                &edge.hrefs.label,
                &mut layers,
                &mut Vec::new(),
            )?;
            let part_refs: Vec<_> = edge.parts.iter().collect();
            for (region, points) in curves::region_samples(&part_refs, UNIT)? {
                targets.extend(points.into_iter().map(|at| Target {
                    at,
                    size: 8.0,
                    href: edge.hrefs.regions[region].clone(),
                }));
            }
            if let Some(index) = edge.fixed_arrow {
                let mark = &heads[&index];
                let shaft = &shafts[&mark.carrier].shaft.path;
                Self::paint_arrow(drawing, mark, shaft, &edge.hrefs.label, &mut layers)?;
                for (region, points) in curves::region_samples(&[shaft], UNIT)? {
                    targets.extend(points.into_iter().map(|at| Target {
                        at,
                        size: 8.0,
                        href: edge.hrefs.regions[region].clone(),
                    }));
                }
            }
        }
        let mut chosen = Vec::new();
        for (placement, &choice) in placements.iter().zip(&searched.choices) {
            let candidate = &placement.candidates[choice];
            if let Some(index) = candidate.arrow_mark {
                let mark = &heads[&index];
                let arrow = &shafts[&mark.carrier].shaft.path;
                Self::paint_arrow(
                    &scene.edges[placement.edge],
                    mark,
                    arrow,
                    &placement.href,
                    &mut layers,
                )?;
                let regions = &region_hrefs
                    .iter()
                    .find(|(edge, _)| *edge == placement.edge)
                    .ok_or("momentum arrow without hrefs")?
                    .1;
                for (region, points) in curves::region_samples(&[arrow], UNIT)? {
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
                    href: output::node_href(node, &laid.edges, &drawing.details),
                });
            }
            targets.push(Target {
                at,
                size: (drawing.radius * UNIT * 2.0).max(10.0),
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
        scene: &Scene,
        typeset: &Typeset,
        laid: &Laid,
        edge: &TypstDotEdge,
        carrier: &[[f64; 2]],
        anchor: [f64; 2],
    ) -> Result<(BezPath, Vec<BezPath>), String> {
        let outset = |end: &TypstDotEndpoint| {
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
                let half = curves::length(&curve) / 2.0;
                let source =
                    curves::trim_routed(&curves::trim(&curve, 0.0, half)?, source_radius, 0.0)?;
                let sink =
                    curves::trim_routed(&curves::trim(&curve, half, 0.0)?, 0.0, sink_radius)?;
                // `layer(curve, outsets)`: one trimmed window, reassembled from its cubics.
                let visible = curves::from_cubics(&curves::cubics(
                    &curves::windows(&curve, &[(source_radius, sink_radius)])?.remove(0),
                ));
                Ok((visible, vec![source, sink]))
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

    /// The base layer: one authoritative (possibly decorated) shaft, then
    /// the particle-flow arrowhead at the source/sink split or the center.
    fn paint_edge(
        drawing: &EdgeDrawing,
        shaft: &BezPath,
        mark: Option<&kurvst::marks::MarkGeometryOutput>,
        href: &str,
        layers: &mut Vec<Element>,
        strokes: &mut Vec<CollisionStroke>,
    ) -> Result<(), String> {
        let mut heads = Vec::new();
        if let Some(mark) = mark {
            for drawable in &mark.paths {
                if !drawable.outline.path.is_empty() {
                    layers.push(Element::HitBox {
                        bounds: Bounds::painted(&drawable.outline.path),
                        href: href.to_owned(),
                    });
                }
                strokes.push(CollisionStroke {
                    radius: 0.0,
                    path: drawable.outline.path.clone(),
                });
            }
            Self::paint_head(
                mark.clone(),
                &drawing.flow_arrow_paints,
                drawing.stroke.width,
                &mut heads,
            )?;
        }
        let radius = drawing.stroke.width / UNIT / 2.0;
        if let Some(pattern) = &drawing.pattern {
            let path = curves::pattern(shaft, pattern)?;
            strokes.push(CollisionStroke {
                radius,
                path: path.clone(),
            });
            layers.push(Element::Path {
                path,
                stroke: drawing.stroke.clone(),
            });
        } else {
            strokes.push(CollisionStroke {
                radius,
                path: shaft.clone(),
            });
            layers.push(Element::Path {
                path: shaft.clone(),
                stroke: drawing.stroke.clone(),
            });
        }
        layers.extend(heads);
        Ok(())
    }

    fn paint_head(
        mark: kurvst::marks::MarkGeometryOutput,
        paints: &[MarkPaint],
        width: f64,
        layers: &mut Vec<Element>,
    ) -> Result<(), String> {
        if !paints.is_empty() && paints.len() != mark.paths.len() {
            return Err("mark paints must match the flattened mark path count".into());
        }
        for (index, drawable) in mark.paths.into_iter().enumerate() {
            layers.push(Element::Mark {
                drawable,
                paint: paints.get(index).cloned().unwrap_or_default(),
                width,
            });
        }
        Ok(())
    }

    fn paint_arrow(
        drawing: &EdgeDrawing,
        mark: &kurvst::marks::MarkGeometryOutput,
        shaft: &BezPath,
        href: &str,
        layers: &mut Vec<Element>,
    ) -> Result<(), String> {
        let stroke = arrow_stroke();
        layers.push(Element::Path {
            path: shaft.clone(),
            stroke,
        });
        if !mark.footprint.path.is_empty() {
            layers.push(Element::HitBox {
                bounds: Bounds::painted(&mark.footprint.path),
                href: href.to_owned(),
            });
        }
        Self::paint_head(mark.clone(), &drawing.momentum_arrow_paints, 1.0, layers)?;
        Ok(())
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
            flow_arrow: EdgeDrawing::default_flow_arrow(),
            momentum_arrow: EdgeDrawing::default_momentum_arrow(),
            flow_arrow_paints: Vec::new(),
            momentum_arrow_paints: Vec::new(),
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
    fn four_loop_scalar_cube_with_default_momentum_marks() {
        let mut scene = bubble(true);
        let node = scene.graph.nodes[0].clone();
        scene.graph.nodes = (0..8)
            .map(|index| {
                let mut node = node.clone();
                node.index = Some(index);
                node.name = Some(format!("v{index}"));
                node
            })
            .collect();
        scene.nodes = vec![scene.nodes[0].clone(); 8];
        let endpoint = scene.graph.edges[1].source.clone().unwrap();
        let edge = scene.graph.edges[1].clone();
        // ext->a; b->ext; b->c; c->d; d->a; e->f; f->g;
        // g->h; h->e; a->e; b->f; c->g; d->h.
        scene.graph.edges = [
            (None, Some(0)),
            (Some(1), None),
            (Some(1), Some(2)),
            (Some(2), Some(3)),
            (Some(3), Some(0)),
            (Some(4), Some(5)),
            (Some(5), Some(6)),
            (Some(6), Some(7)),
            (Some(7), Some(4)),
            (Some(0), Some(4)),
            (Some(1), Some(5)),
            (Some(2), Some(6)),
            (Some(3), Some(7)),
        ]
        .into_iter()
        .enumerate()
        .map(|(id, (source, sink))| {
            let mut edge = edge.clone();
            edge.id = Some(id);
            let end = |node| {
                let mut end = endpoint.clone();
                end.node = node;
                end
            };
            edge.source = source.map(end);
            edge.sink = sink.map(end);
            edge
        })
        .collect();
        let mut drawing = scene.edges[2].clone();
        drawing.flow = None;
        drawing.pattern = None;
        drawing.label = None;
        scene.edges = vec![drawing; 13];
        scene.layout.insert("steps".into(), 100.into());
        for momentum in [false, true] {
            for edge in &mut scene.edges {
                edge.momentum = momentum;
            }
            scene.render(&typeset(&scene)).unwrap();
        }
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
    fn one_candidate_and_one_selected_geometry_call_cover_the_entire_drawing() {
        marks::Batch::take_calls();
        let scene = bubble(true);
        scene.render(&typeset(&scene)).unwrap();
        assert_eq!(
            marks::Batch::take_calls(),
            [
                kurvst::marks::MarkGeometryMode::Candidates,
                kurvst::marks::MarkGeometryMode::Selected,
            ]
        );
    }

    #[test]
    fn normalized_catalogue_config_and_geometry_options_render_natively() {
        for shape in [
            "triangle", "straight", "stealth", "round", "tikz", "barb", "hooks", "bar", "bracket",
            "circle", "square", "diamond", "rays",
        ] {
            let mut scene = bubble(true);
            let mut mark = serde_json::json!({
                "shape": shape, "stroke": true, "fit": "bend", "shorten": 0.6
            });
            let fields = match shape {
                "triangle" | "stealth" | "round" => "length width rev fill",
                "straight" | "bracket" => "length width rev",
                "tikz" => "width",
                "barb" | "hooks" => "width rev",
                "bar" => "width align",
                "circle" | "square" | "diamond" => "length width align fill",
                "rays" => "length align",
                _ => unreachable!(),
            };
            for field in fields.split_whitespace() {
                mark[field] = match field {
                    "length" => serde_json::json!({"points":3,"ratio":2}),
                    "width" => serde_json::json!({"points":2,"ratio":1}),
                    "align" => "end".into(),
                    _ => true.into(),
                };
            }
            match shape {
                "triangle" | "stealth" | "round" => mark["inset"] = 0.25.into(),
                "barb" | "hooks" => mark["arc"] = 2.5.into(),
                "rays" => {
                    mark["n"] = 5.into();
                    mark["phase"] = 0.1.into();
                }
                _ => {}
            }
            let config = Config::from_json(
                &serde_json::json!({
                    "style": {"edge-style": {"flow-arrow": mark, "momentum-arrow": mark}}
                })
                .to_string(),
            )
            .unwrap();
            config.apply(&mut scene).unwrap();
            let svg = scene.render(&typeset(&scene)).unwrap();
            assert!(!svg.contains("NaN") && !svg.contains("inf"), "{shape}");
            assert_eq!(svg.matches("stroke-miterlimit=").count(), 5, "{shape}");
        }
    }

    #[test]
    fn nested_combine_preserves_independent_leaf_paints_and_hit_areas() {
        let mut scene = bubble(true);
        let config = Config::from_json(&serde_json::json!({
            "style": {"edge-style": {
                "momentum-arrow": {"shape":"combine","fit":"bend","parts":[
                    {"shape":"circle","fill":true,"stroke":true},
                    {"gap":{"points":2,"ratio":0}},
                    {"shape":"combine","parts":[{"shape":"bar"},{"shape":"triangle","fill":true,"stroke":true}]}
                ]},
                "momentum-arrow-paints":[
                    {"fill":"red","stroke":"blue"}, {"stroke":"green"},
                    {"fill":"yellow","stroke":"purple"}
                ]
            }}
        }).to_string()).unwrap();
        config.apply(&mut scene).unwrap();
        let svg = scene.render(&typeset(&scene)).unwrap();
        assert_eq!(svg.matches(r#"fill="red" stroke="blue""#).count(), 4);
        assert_eq!(svg.matches(r#"fill="yellow" stroke="purple""#).count(), 4);
        let interactive = Scene::interactive_svg(&svg).unwrap();
        assert!(interactive.contains("data-linnet-kind=\"edge\""));
        assert_eq!(interactive.matches("stroke-miterlimit=").count(), 13);
    }

    #[test]
    fn obsolete_and_unsupported_mark_fields_fail_with_replacement() {
        for mark in [
            serde_json::json!({"kind":"mark","symbol":">"}),
            serde_json::json!({"shape":"triangle","scale":1}),
            serde_json::json!({"shape":"triangle","anchor":"tip"}),
            serde_json::json!({"shape":"triangle","shorten_to":1}),
            serde_json::json!({"shape":"triangle","arc":1}),
            serde_json::json!({"shape":"bar","length":"2pt"}),
            serde_json::json!({"shape":"rays","n":2.5}),
            serde_json::json!({"shape":"combine","parts":[{"shape":"triangle","fit":"chord"}]}),
        ] {
            let mut scene = bubble(true);
            let config = Config::from_json(
                &serde_json::json!({
                    "style":{"edge-style":{"flow-arrow":mark}}
                })
                .to_string(),
            )
            .unwrap();
            let error = config.apply(&mut scene).unwrap_err();
            assert!(error.contains("linnet.Mark"), "{error}");
        }
        let mut scene = bubble(true);
        let config =
            Config::from_json(r#"{"style":{"edge-style":{"momentum-arrow-paints":[{},{}]}}}"#)
                .unwrap();
        assert!(config.apply(&mut scene).unwrap_err().contains("flattened"));
    }

    #[test]
    fn label_candidate_bounds_are_from_selected_painted_geometry() {
        let mut path = BezPath::new();
        path.move_to((0.0, 0.0));
        path.curve_to((0.0, 1.0), (3.0, 2.0), (3.0, 0.0));
        let mut drawing = bubble(true).edges.remove(0);
        drawing.momentum_arrow = serde_json::from_value(serde_json::json!({
            "shape":"combine","fit":"bend","shorten":0.5,
            "parts":[{"shape":"circle"},{"gap":{"points":2,"ratio":0}},{"shape":"stealth"}]
        }))
        .unwrap();
        let mut batch = marks::Batch::new();
        let template = batch
            .template(&drawing.momentum_arrow, &arrow_stroke())
            .unwrap();
        let proposals = labels::MomentumCandidates::new(
            &path,
            [10.0, 8.0, 8.0, 12.0],
            1.0,
            template,
            &mut batch,
        )
        .unwrap();
        let geometry = batch
            .geometry(kurvst::marks::MarkGeometryMode::Candidates)
            .unwrap();
        let candidates = proposals.resolve(&geometry.marks).unwrap();
        for candidate in candidates {
            let geometry = batch.selected(&[candidate.arrow_mark.unwrap()]).unwrap();
            let mark = &geometry.marks[0];
            let shaft = &geometry.shafts[0];
            let mut layers = Vec::new();
            Drawing::paint_arrow(
                &drawing,
                mark,
                &shaft.shaft.path,
                "#linnet-edge-0",
                &mut layers,
            )
            .unwrap();
            let bounds = [
                Bounds::painted(&mark.footprint.path),
                Bounds::painted(&shaft.footprint.path),
            ];
            let candidate = candidate.arrow_bounds[0];
            for bound in bounds {
                assert!(candidate.left <= bound.left && candidate.right >= bound.right);
                assert!(candidate.bottom <= bound.bottom && candidate.top >= bound.top);
            }
            assert_eq!(
                layers
                    .iter()
                    .filter(|layer| matches!(layer, Element::Path { .. }))
                    .count(),
                1
            );
        }
    }
}

impl NodeDrawing {
    fn size(&self, typeset: &Typeset) -> (f64, f64) {
        let Some(label) = self.label.map(|i| &typeset.pages[i]) else {
            return (self.radius, self.radius);
        };
        let (w, h) = (
            (label.width / UNIT / 2.0 + 0.25).max(self.radius),
            (label.height / UNIT / 2.0 + 0.25).max(self.radius),
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
