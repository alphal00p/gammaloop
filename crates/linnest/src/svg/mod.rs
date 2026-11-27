//! Native SVG drawing of prepared graphs, the host-side counterpart of
//! `draw.typ`: native layout, Kurvst geometry, the shared annotation search
//! and SVG output. The host measures the MathML label pages a [`Scene`] lists.
//!
//! See [`Scene::typeset_mathml`] for the shared label/measurement contract.
mod config;
mod curves;
pub use config::Config;
pub use highlight::{
    CountertermScheme, HighlightDetails, HighlightDrawing, HighlightOverlap, HighlightPalette,
    HighlightRegion, HighlightStroke, HighlightStyle, HullOptions, HullRegion, LaneRegion,
    SlotPaint, UnknownCountertermScheme, UnknownHighlightOverlap, UnknownHighlightStyle,
};
mod highlight;
mod hulls;
mod interactive;
mod labels;
mod lanes;
mod marks;
mod output;

use kurbo::{Arc as CurveArc, BezPath, Point, Rect, Shape, Vec2};
use linnet::half_edge::involution::EdgeIndex;
use linnet::half_edge::layout::impred::{EdgeLabel, ImpredConfig};
use linnet::half_edge::subgraph::{ModifySubSet, SuBitGraph};
use roxmltree::Document;
use serde_json::{Map, Value};
use std::collections::BTreeMap;
use std::f64::consts::{FRAC_PI_2, TAU};
use std::sync::Arc;

use crate::{ImpredRun, TypstDotEdge, TypstDotEndpoint, TypstDotNode, TypstGraph, TypstGraphSpec};
use labels::{Anchor, Bounds, Placement, Stroke as CollisionStroke};
pub use labels::{LabelError, LabelMetrics, LabelPage, LabelSource, Typeset};
use lanes::{Carrier as LaneCarrier, Station as LaneStation};
pub use output::Details;
use output::INK;

/// Drawing unit in points: Linnest's default 1.5em at 9pt text.
const UNIT: f64 = 13.5;
/// Vertex radius and edge trim, in drawing units.
const NODE_RADIUS: f64 = 0.16;

/// A prepared graph with its paint and label sources, drawn by [`Scene::render`].
/// See [`Self::typeset_mathml`] for a complete measured-label example.
pub struct Scene {
    /// Highlight drawing state. Initialize with `HighlightDrawing::default()`.
    pub highlights: HighlightDrawing,
    /// Layout graph: sizes and pins as statements, in graph-spec form.
    pub graph: TypstGraphSpec,
    /// Drawing per node and edge, in spec order.
    pub nodes: Vec<NodeDrawing>,
    pub edges: Vec<EdgeDrawing>,
    /// Typst source preceding the label pages: imports, text style, definitions.
    pub preamble: String,
    /// Label source of the title row above the drawing.
    pub title: Option<LabelSource>,
    /// Label source of each page, in measurement order after the title.
    pub pages: Vec<LabelSource>,
    /// Layout options; `layout-algo` selects `dot` or `impred` (the default).
    pub layout: Map<String, Value>,
    /// Further passes, applied in order to the previous pass's geometry.
    pub layout_passes: Vec<Map<String, Value>>,
    /// Edges defining the layered hierarchy, in spec order. Other edges are
    /// routed around it without changing the ranks. `None` uses all edges.
    pub layout_edges: Option<Vec<usize>>,
    /// Refine the layout around drawn labels in short warm ImPrEd passes.
    pub label_feedback: bool,
}

#[derive(Clone)]
pub struct NodeDrawing {
    /// Ordered paints for a segmented highlight ring.
    pub highlight_ring: Vec<String>,
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
    /// Resolved strokes in global layer order, outermost first.
    pub highlight_strokes: Vec<HighlightStroke>,
    pub pattern: Option<Pattern>,
    /// Particle-flow arrowhead along (`true`) or against the edge orientation.
    pub flow: Option<bool>,
    /// Momentum arrow beside the edge, carrying the label.
    pub momentum: bool,
    /// Distance from the edge to its parallel momentum arrow, in drawing units.
    pub momentum_offset: f64,
    pub label: Option<usize>,
    /// Inspection fields after the edge's own identity.
    pub details: Details,
}

impl EdgeDrawing {
    fn has_bands(&self) -> bool {
        self.highlight_strokes
            .first()
            .is_some_and(|stroke| stroke.offset.is_some())
    }

    /// Default clearance used by the shared physics renderer.
    pub const DEFAULT_MOMENTUM_OFFSET: f64 = 0.35;
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
    /// The optional vector-export document, with measured title first.
    ///
    /// See [`Self::typeset_mathml`] for constructing a scene, and [`Self::typeset`]
    /// for importing the pages produced by an external Typst compiler.
    ///
    /// # Errors
    /// Returns [`LabelError::MissingTypstSource`] if any label has only MathML.
    pub fn label_document(&self) -> Result<String, LabelError> {
        let pages: Vec<String> = self
            .title
            .iter()
            .chain(&self.pages)
            .enumerate()
            .map(|(page, source)| {
                source
                    .typst
                    .as_deref()
                    .map(LabelPage::measured)
                    .ok_or(LabelError::MissingTypstSource { page })
            })
            .collect::<Result<_, _>>()?;
        Ok(format!(
            "#set page(width: auto, height: auto, margin: 0pt, fill: none)\n{}\n{}",
            self.preamble,
            pages.join("\n#pagebreak()\n")
        ))
    }

    /// Import host-measured MathML boxes, title first, in points.
    ///
    /// # Errors
    /// Returns [`LabelError`] when the number of measurements is wrong, a
    /// rectangle is invalid, or a label is malformed MathML or contains foreign
    /// elements. MathML does not need a Typst source.
    ///
    /// ```
    /// use linnest::svg::{LabelMetrics, LabelSource, Scene};
    /// use serde_json;
    /// use std::error::Error;
    /// let scene = Scene {
    ///     highlights: Default::default(),
    ///     graph: serde_json::from_str(r#"{"nodes":[],"edges":[]}"#)?,
    ///     nodes: Vec::new(), edges: Vec::new(), preamble: String::new(),
    ///     title: Some(LabelSource::text("Measured title")), pages: Vec::new(),
    ///     layout: Default::default(), layout_passes: Vec::new(),
    ///     layout_edges: None, label_feedback: false,
    /// };
    /// // Supply the actual browser's measurements, in title/label order.
    /// let metrics = LabelMetrics { width: 65.0, height: 9.0, baseline: 7.0 };
    /// let labels = scene.typeset_mathml(&[metrics])?;
    /// assert_eq!(labels.title.as_ref().map(|title| title.metrics), Some(metrics));
    /// assert!(scene.label_document()?.contains("Measured title"));
    /// # Ok::<(), Box<dyn Error>>(())
    /// ```
    pub fn typeset_mathml(&self, metrics: &[LabelMetrics]) -> Result<Typeset, LabelError> {
        let expected = self.pages.len() + usize::from(self.title.is_some());
        if metrics.len() != expected {
            return Err(LabelError::MeasurementCount {
                expected,
                actual: metrics.len(),
            });
        }
        let mut pages = self
            .title
            .iter()
            .chain(&self.pages)
            .zip(metrics)
            .map(|(source, metrics)| source.typeset_mathml(*metrics));
        let title = if self.title.is_some() {
            Some(pages.next().ok_or(LabelError::MissingTitle)??)
        } else {
            None
        };
        Ok(Typeset {
            title,
            pages: pages.collect::<Result<_, _>>()?,
            defs: String::new(),
        })
    }

    /// Read the SVG pages compiled from [`Scene::label_document`].
    ///
    /// See [`Typeset::read`] for an example of the measured SVG page protocol.
    ///
    /// # Errors
    /// Returns [`LabelError`] for malformed or unmeasured SVG pages, a missing
    /// title, or a page count that differs from this scene's labels.
    pub fn typeset(&self, svgs: &[String]) -> Result<Typeset, LabelError> {
        let typeset = Typeset::read(svgs, self.title.is_some())?;
        if typeset.pages.len() != self.pages.len() {
            return Err(LabelError::PageCount {
                expected: self.pages.len(),
                actual: typeset.pages.len(),
            });
        }
        Ok(typeset)
    }

    /// Embed this scene's computed SVG paths in a self-contained Typst document.
    /// Supply outlined label pages via [`Self::typeset`] for vector text.
    ///
    /// # Errors
    /// Returns a layout or geometry failure from [`Self::render`].
    pub fn render_typst(&self, typeset: &Typeset) -> Result<String, String> {
        self.render(typeset).map(|svg| Self::svg_typst_source(&svg))
    }

    /// Wrap a static SVG as a Typst image. See [`Self::render_typst`].
    pub fn svg_typst_source(svg: &str) -> String {
        format!("#set page(width: auto, height: auto, margin: 0pt)\n#image(bytes({svg:?}))")
    }

    /// Lay out and draw the scene with its typeset pages; one SVG document.
    pub fn render(&self, typeset: &Typeset) -> Result<String, String> {
        self.validate_hull_radii()?;
        let mut spec = self.graph.clone();
        for (node, drawing) in spec.nodes.iter_mut().zip(&self.nodes) {
            let (w, h) = drawing.size(typeset);
            node.statements
                .insert("layout-width".into(), (2.0 * w).to_string());
            node.statements
                .insert("layout-height".into(), (2.0 * h).to_string());
        }
        let mut graph = TypstGraph::from_spec(spec)?;
        let mut run: Option<ImpredRun> = None;
        let mut laid = None;
        for options in std::iter::once(&self.layout).chain(&self.layout_passes) {
            match options
                .get("layout-algo")
                .and_then(Value::as_str)
                .unwrap_or("impred")
            {
                "impred" => {
                    let options = Value::Object(options.clone());
                    if let Some(run) = &mut run {
                        run.refine(&options)?;
                    } else {
                        run = Some(ImpredRun::seeded(graph, &options)?);
                    }
                    let current = run.as_ref().unwrap();
                    graph = current
                        .session
                        .apply_impred_layout(&current.layout)
                        .map_err(|error| error.to_string())?;
                    laid = Some(Laid::new(current)?);
                }
                "dot" | "force" | "anneal" | "tree" | "stable-layered" | "railroad" => {
                    laid = Some(Laid::layered(
                        &mut graph,
                        options,
                        self.layout_edges.as_ref(),
                    )?);
                    run = None;
                }
                algo => return Err(format!("unsupported SVG layout algorithm {algo:?}")),
            }
        }
        let Some(mut run) = run else {
            return Drawing::new(self, typeset, &laid.unwrap())?.svg(typeset);
        };
        let drawing = Drawing::new(self, typeset, &Laid::new(&run)?)?;
        if !self.label_feedback {
            return drawing.svg(typeset);
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
        best.svg(typeset)
    }
}

/// The laid-out graph: node positions, edge records and complete carriers.
struct Laid {
    nodes: Vec<TypstDotNode>,
    edges: Vec<TypstDotEdge>,
    carriers: Vec<Vec<[f64; 2]>>,
}

impl Laid {
    fn layered(
        graph: &mut TypstGraph,
        options: &Map<String, Value>,
        edges: Option<&Vec<usize>>,
    ) -> Result<Self, String> {
        graph.layout_config = serde_json::from_value(Value::Object(options.clone()))
            .map_err(|error| format!("invalid layered layout options: {error}"))?;
        let selected = edges
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
        highlight: bool,
    },
    /// Momentum decorations retain their physical edge identity for overlays.
    Momentum {
        href: String,
        path: BezPath,
    },
    Triangle([Point; 3]),
    BeginGhost(f64),
    EndGhost,
    Hull {
        path: BezPath,
        paint: String,
        details: HighlightDetails,
    },
    Lane {
        slot: usize,
        path: BezPath,
        paint: String,
        details: HighlightDetails,
        width: f64,
    },
    LaneCrossings(usize),
    LaneHalo {
        center: Point,
        radius: f64,
    },
    Counterterm {
        center: Point,
        paint: String,
        details: HighlightDetails,
    },
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
    recolour_stamp: bool,
    layers: Vec<Element>,
    targets: Vec<Target>,
    placements: Vec<Placement>,
    chosen: Vec<Bounds>,
    /// Overlapping labels, line hits, label pairs.
    overlaps: (usize, usize, usize),
}

impl Drawing {
    /// Serialize only the selected feedback pass, preserving painting order.
    fn svg(&self, typeset: &Typeset) -> Result<String, String> {
        output::svg(typeset, &self.layers, &self.targets, self.recolour_stamp)
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
        let mut highlight_layers = BTreeMap::new();
        let mut overlay_marks = Vec::new();
        let mut targets = Vec::new();
        let mut strokes = Vec::new();
        let mut placements = Vec::new();
        let mut fixed_arrows = Vec::new();
        let mut region_hrefs = Vec::new();
        let mut hull_paths = BTreeMap::new();
        let mut drawn_lines = Vec::new();
        let mut lane_paths = BTreeMap::new();
        for (edge, carrier) in laid.edges.iter().zip(&laid.carriers) {
            let drawing = scene
                .edges
                .get(edge.edge)
                .ok_or("scene edges do not match the layout graph")?;
            let hrefs = output::edge_hrefs(edge, &drawing.details);
            let anchor = edge.pos.as_ref().map_or([0.0; 2], |p| [p.x, p.y]);
            let (visible, parts) = Self::edge_paths(scene, typeset, laid, edge, carrier, anchor)?;
            let member = scene
                .highlights
                .hulls
                .iter()
                .any(|region| region.edges.contains(&edge.edge));
            let contracted = scene
                .highlights
                .hulls
                .iter()
                .any(|region| region.details.contracted && region.edges.contains(&edge.edge));
            if member {
                hull_paths.insert(edge.edge, visible.clone());
            }
            drawn_lines.push(visible.clone());
            if let (Some(source), Some(sink)) = (&edge.source, &edge.sink) {
                if scene
                    .highlights
                    .lanes
                    .iter()
                    .any(|region| region.edges.contains(&edge.edge))
                {
                    lane_paths.insert(
                        edge.edge,
                        LaneCarrier {
                            nodes: [source.node, sink.node],
                            path: visible.clone(),
                        },
                    );
                }
            }
            let part_refs: Vec<&BezPath> = parts.iter().collect();
            for (region, points) in curves::region_samples(&part_refs, UNIT)? {
                targets.extend(points.into_iter().map(|at| Target {
                    at,
                    size: 8.0,
                    href: hrefs.regions[region].clone(),
                }));
            }
            if contracted {
                let opacity = if scene.highlights.recolour_stamp {
                    0.62
                } else {
                    0.35
                };
                layers.push(Element::BeginGhost(opacity));
            }
            Self::paint_edge(
                drawing,
                &visible,
                &parts,
                &mut layers,
                &mut highlight_layers,
                &mut overlay_marks,
                &mut strokes,
            )?;
            if contracted {
                layers.push(Element::EndGhost);
            }
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
                .map(|page| (page, typeset.pages[page].metrics));
            let paired = edge.source.is_some() && edge.sink.is_some();
            if drawing.momentum && !(paired && label.is_some()) {
                // External arrows are fixed decorations, as are internal ones
                // without a label to carry them.
                let arrow = curves::layer(
                    &visible,
                    drawing.momentum_offset * side,
                    Some(labels::ARROW_WINDOW),
                    0.0,
                )?;
                fixed_arrows.extend(labels::arrow_footprint(&arrow));
                let paint = if drawing.highlight_strokes.is_empty() {
                    &mut layers
                } else {
                    &mut overlay_marks
                };
                paint.push(Element::Momentum {
                    href: hrefs.label.clone(),
                    path: arrow,
                });
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
                    let (candidates, carriers) = labels::momentum_candidates(
                        &visible,
                        metrics,
                        side,
                        drawing.momentum_offset,
                    )?;
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

        let highlight_elements = highlight_layers.into_values().flatten();
        if scene.edges.iter().any(EdgeDrawing::has_bands) {
            layers.splice(..0, highlight_elements);
        } else {
            layers.extend(highlight_elements);
        }
        layers.extend(overlay_marks);
        if !scene.highlights.lanes.is_empty() {
            let stations: Vec<_> = scene
                .nodes
                .iter()
                .enumerate()
                .map(|(node, drawing)| {
                    let [x, y] = laid.node(node);
                    let (width, height) = drawing.size(typeset);
                    LaneStation {
                        center: Point::new(x, y),
                        radius: width.max(height),
                    }
                })
                .collect();
            layers.splice(
                ..0,
                lanes::draw(
                    &scene.highlights.lanes,
                    &lane_paths,
                    &stations,
                    scene.edges.len(),
                )?,
            );
        }
        let mut built = Vec::new();
        for region in &scene.highlights.hulls {
            let paths: Vec<_> = region
                .edges
                .iter()
                .filter_map(|index| hull_paths.get(index))
                .filter(|path| !path.elements().is_empty())
                .collect();
            let vertices: Vec<_> = laid
                .nodes
                .iter()
                .filter(|node| region.nodes.contains(&node.node))
                .map(|node| {
                    let [x, y] = laid.node(node.node);
                    Point::new(x, y)
                })
                .collect();
            let path = hulls::union(&paths, &vertices, region.radius / UNIT)?;
            built.push((region, path, vertices));
        }
        // A counterterm glyph keeps clear of every drawn line, including the
        // outlines of all hulls (holes and nested levels), not only edges.
        let mut obstacle_lines: Vec<&BezPath> = drawn_lines.iter().collect();
        obstacle_lines.extend(built.iter().map(|(_, path, _)| path));
        let all_vertices: Vec<_> = laid
            .nodes
            .iter()
            .map(|node| {
                let [x, y] = laid.node(node.node);
                Point::new(x, y)
            })
            .collect();
        // The glyph's drawn radius: miniatures scale it with their canvas.
        let content = obstacle_lines
            .iter()
            .map(|path| path.bounding_box())
            .chain(
                all_vertices
                    .iter()
                    .map(|point| Rect::from_points(*point, *point)),
            )
            .reduce(|bounds, rect| bounds.union(rect))
            .unwrap_or(Rect::ZERO);
        let glyph_radius = output::glyph_radius(content.width(), scene.highlights.recolour_stamp);
        let mut regions = Vec::new();
        let mut counterterms = Vec::new();
        for (region, path, vertices) in &built {
            if region.details.contracted {
                let center = if vertices.is_empty() {
                    // For a fully trimmed or empty selection, retain the glyph
                    // at the selected edge anchors, or at the layout centre.
                    if path.elements().is_empty() {
                        let anchors: Vec<_> = laid
                            .edges
                            .iter()
                            .filter(|edge| region.edges.contains(&edge.edge))
                            .filter_map(|edge| edge.pos.as_ref().map(|p| Point::new(p.x, p.y)))
                            .collect();
                        let positions = if anchors.is_empty() {
                            laid.nodes
                                .iter()
                                .filter_map(|node| node.pos.as_ref().map(|p| Point::new(p.x, p.y)))
                                .collect()
                        } else {
                            anchors
                        };
                        positions
                            .split_first()
                            .map_or(Point::ORIGIN, |(&first, others)| {
                                others
                                    .iter()
                                    .fold(Rect::from_points(first, first), |bounds, &p| {
                                        bounds.union_pt(p)
                                    })
                                    .center()
                            })
                    } else {
                        path.bounding_box().center()
                    }
                } else {
                    let sum = vertices
                        .iter()
                        .fold(Vec2::ZERO, |sum, point| sum + point.to_vec2());
                    let centroid = (sum / vertices.len() as f64).to_point();
                    hulls::glyph_point(
                        path,
                        &obstacle_lines,
                        &all_vertices,
                        NODE_RADIUS,
                        glyph_radius,
                        centroid,
                    )
                };
                counterterms.push(Element::Counterterm {
                    center,
                    paint: region.paint.clone(),
                    details: region.details.clone(),
                });
            }
            if !path.elements().is_empty() {
                regions.push(Element::Hull {
                    path: path.clone(),
                    paint: region.paint.clone(),
                    details: region.details.clone(),
                });
            }
        }
        layers.splice(..0, regions);

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
                layers.push(Element::Momentum {
                    href: placement.href.clone(),
                    path: arrow,
                });
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
            let is_station = !drawing.rectangular
                && lane_paths
                    .values()
                    .any(|carrier| carrier.nodes.contains(&node.node));
            if is_station {
                layers.push(Element::LaneHalo {
                    center: Point::new(at[0], at[1]),
                    radius: drawing.size(typeset).0,
                });
            }
            layers.push(Element::Node {
                id: drawing
                    .details
                    .get("node")
                    .and_then(Value::as_u64)
                    .map_or(node.node, |id| id as usize),
                at,
                size: drawing.size(typeset),
                rectangular: drawing.rectangular,
                fill: if is_station {
                    INK.to_owned()
                } else {
                    drawing.fill.clone()
                },
                stroke: drawing.stroke.clone(),
            });
            let (width, height) = drawing.size(typeset);
            let radius = width.max(height) + 1.15 / UNIT;
            let count = drawing.highlight_ring.len();
            for (index, paint) in drawing.highlight_ring.iter().enumerate() {
                let step = TAU / count as f64;
                let gap = if count == 1 {
                    0.0
                } else {
                    0.035_f64.min(step / 4.0)
                };
                let arc = CurveArc {
                    center: Point::new(at[0], at[1]),
                    radii: Vec2::new(radius, radius),
                    start_angle: FRAC_PI_2 - index as f64 * step - gap,
                    sweep_angle: -step + 2.0 * gap,
                    x_rotation: 0.0,
                };
                let start = Point::new(
                    at[0] + radius * arc.start_angle.cos(),
                    at[1] + radius * arc.start_angle.sin(),
                );
                let mut path = BezPath::new();
                path.move_to(start);
                path.extend(arc.append_iter(curves::ACCURACY));
                layers.push(Element::Path {
                    path,
                    stroke: Stroke {
                        paint: paint.clone(),
                        width: 2.5,
                        dash: Dash::Solid,
                        round_cap: false,
                    },
                    highlight: true,
                });
            }
            if let Some(page) = drawing.label {
                let label = &typeset.pages[page];
                layers.push(Element::Label {
                    page,
                    bounds: Bounds {
                        left: at[0] - label.metrics.width / UNIT / 2.0,
                        right: at[0] + label.metrics.width / UNIT / 2.0,
                        bottom: at[1] - label.metrics.height / UNIT / 2.0,
                        top: at[1] + label.metrics.height / UNIT / 2.0,
                    },
                    href: None,
                });
            }
            // One target covers both the node and its label. Separate
            // overlapping anchors steal clicks from the keyboard-focusable one.
            let (width, height) = drawing.size(typeset);
            let label_size = drawing.label.map_or(0.0, |page| {
                let metrics = typeset.pages[page].metrics;
                metrics.width.max(metrics.height)
            });
            targets.push(Target {
                at,
                size: (width.max(height) * UNIT * 2.0).max(label_size).max(10.0),
                href: href.into(),
            });
        }
        layers.extend(counterterms);
        let overlaps = labels::overlaps(&placements, &chosen, &searched.lines);
        Ok(Self {
            recolour_stamp: scene.highlights.recolour_stamp,
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
        highlight_layers: &mut BTreeMap<usize, Vec<Element>>,
        overlay_marks: &mut Vec<Element>,
        strokes: &mut Vec<CollisionStroke>,
    ) -> Result<(), String> {
        let bands = drawing.has_bands();
        if bands && !visible.elements().is_empty() {
            for band in &drawing.highlight_strokes {
                let offset = band.offset.unwrap_or_default() / UNIT;
                let path = if offset == 0.0 {
                    visible.clone()
                } else {
                    curves::parallel(visible, offset)?
                };
                strokes.push(CollisionStroke {
                    radius: band.stroke.width / UNIT / 2.0,
                    path: path.clone(),
                });
                highlight_layers.entry(0).or_default().push(Element::Path {
                    path,
                    stroke: band.stroke.clone(),
                    highlight: true,
                });
            }
        }
        let outermost = drawing
            .highlight_strokes
            .first()
            .filter(|_| !bands)
            .map_or(&drawing.stroke, |highlight| &highlight.stroke);
        let radius = outermost.width / UNIT / 2.0;
        let paths = if let Some(pattern) = &drawing.pattern {
            vec![curves::pattern(visible, pattern)?]
        } else {
            curves::cubics(visible)
                .into_iter()
                .map(|segment| curves::from_cubics(&[segment]))
                .collect()
        };
        for path in &paths {
            strokes.push(CollisionStroke {
                radius,
                path: path.clone(),
            });
        }
        let base_stroke = (bands
            || drawing.highlight_strokes.is_empty()
            || matches!(layers.last(), Some(Element::BeginGhost(_))))
        .then_some((None, &drawing.stroke));
        let highlights = drawing
            .highlight_strokes
            .iter()
            .filter(|_| !bands)
            .map(|highlight| (Some(highlight.layer), &highlight.stroke));
        for (index, stroke) in base_stroke.into_iter().chain(highlights) {
            let paint = match index {
                Some(index) => highlight_layers.entry(index).or_default(),
                None => &mut *layers,
            };
            paint.extend(paths.iter().map(|path| Element::Path {
                path: path.clone(),
                stroke: stroke.clone(),
                highlight: index.is_some(),
            }));
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
                let paint = if drawing.highlight_strokes.is_empty() {
                    layers
                } else {
                    overlay_marks
                };
                paint.push(Element::Triangle(points));
            }
        }
        Ok(())
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
            (label.metrics.width / UNIT / 2.0 + padding).max(self.radius),
            (label.metrics.height / UNIT / 2.0 + padding).max(self.radius),
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
            let document = Document::parse(figure).map_err(|e| e.to_string())?;
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
    use std::collections::{BTreeMap, BTreeSet};
    use std::slice;

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
            highlight_strokes: Vec::new(),
            stroke: ink(1.0),
            pattern,
            flow,
            momentum,
            momentum_offset: EdgeDrawing::DEFAULT_MOMENTUM_OFFSET,
            label: Some(0),
            details: [("particle", "x")].into_iter().collect(),
        };
        let coil = Pattern::Coil {
            amplitude: 0.15,
            wavelength: 0.45,
            longitudinal_scale: 1.4,
        };
        Scene {
            highlights: Default::default(),
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
                    highlight_ring: Vec::new(),
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
            title: Some(LabelSource::text("bubble")),
            pages: vec![LabelSource::text("x")],
            layout: Map::new(),
            layout_passes: Vec::new(),
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
        scene
            .typeset(&[
                page(r#"<g/><a href="linnest-metrics:10,8,8,12"><rect/></a>"#),
                label,
            ])
            .unwrap()
    }

    #[test]
    fn hulls_cover_parallel_carriers_and_coil_with_escaped_counterterm_details() {
        let mut scene = bubble(false);
        let mut inner = HighlightRegion::new([1].into(), [0, 1].into(), HighlightStyle::Treated);
        inner.details = HighlightDetails {
            label: "γ₂ & <local>".into(),
            description: "local \"subtraction\"".into(),
            contracted: true,
            scheme: CountertermScheme::OnShell,
        };
        let outer = HighlightRegion::new([1, 2].into(), [0, 1].into(), HighlightStyle::ActedOn);
        scene.apply_highlights_with_overlap(
            &[inner, outer],
            &HighlightPalette::new("#ffd166"),
            true,
            HighlightOverlap::Hulls,
        );
        assert!(scene.highlights.hulls[0].radius > scene.highlights.hulls[1].radius);
        let svg = scene.render(&typeset(&scene)).unwrap();
        let doc = Document::parse(&svg).unwrap();
        let regions: Vec<_> = doc
            .descendants()
            .filter(|node| node.attribute("data-linnet-hull") == Some("true"))
            .collect();
        assert_eq!(regions.len(), 2);
        assert_eq!(
            regions[0].last_element_child().unwrap().attribute("fill"),
            Some("#3b82f6")
        );
        assert_eq!(
            regions[1]
                .last_element_child()
                .unwrap()
                .attribute("fill-opacity"),
            Some("0.22")
        );
        assert_ne!(
            regions[0].last_element_child().unwrap().attribute("d"),
            regions[1].last_element_child().unwrap().attribute("d")
        );
        let marker = doc
            .descendants()
            .find(|node| node.attribute("data-linnet-counterterm") == Some("true"))
            .unwrap();
        assert_eq!(marker.attribute("data-scheme"), Some("on-shell"));
        let disc = marker
            .children()
            .find(|node| node.has_tag_name("circle"))
            .unwrap();
        assert_eq!(disc.attribute("fill"), Some("#e4504f"));
        assert_eq!(
            marker.last_element_child().unwrap().attribute("stroke"),
            Some("#ffffff")
        );
        assert_eq!(marker.attribute("data-label"), Some("γ₂ & <local>"));
        assert_eq!(
            marker.first_element_child().unwrap().text(),
            Some("γ₂ & <local> — local \"subtraction\" — on-shell")
        );
        assert!(doc
            .descendants()
            .any(|node| node.attribute("opacity") == Some("0.35")));
    }

    #[test]
    fn recolour_stamp_bounds_enclose_scaled_circles_and_rectangles() {
        let mut scene = bubble(false);
        scene.nodes[1].rectangular = true;
        let mut labels = typeset(&scene);
        labels.title.as_mut().unwrap().metrics.width = 3200.0;
        let region = HighlightRegion::new([].into(), [0].into(), HighlightStyle::Treated);
        scene.apply_highlights_with_overlap(
            &[region],
            &HighlightPalette::new("#ffd166"),
            true,
            HighlightOverlap::Recolour,
        );
        let svg = scene.render(&labels).unwrap();
        let doc = Document::parse(&svg).unwrap();
        let canvas: Vec<f64> = doc
            .root_element()
            .attribute("viewBox")
            .unwrap()
            .split_whitespace()
            .map(|value| value.parse().unwrap())
            .collect();
        let nodes: Vec<_> = doc
            .descendants()
            .filter(|node| node.has_attribute("data-linnet-node"))
            .collect();
        assert_eq!(nodes.len(), 2);
        for node in nodes {
            let value = |key| node.attribute(key).unwrap().parse::<f64>().unwrap();
            let half_stroke = value("stroke-width") / 2.0;
            let (left, top, right, bottom) = if node.has_tag_name("circle") {
                let (x, y, r) = (value("cx"), value("cy"), value("r"));
                (x - r, y - r, x + r, y + r)
            } else {
                let (x, y) = (value("x"), value("y"));
                (x, y, x + value("width"), y + value("height"))
            };
            assert!(left >= half_stroke && top >= half_stroke);
            assert!(right + half_stroke <= canvas[2] && bottom + half_stroke <= canvas[3]);
        }
    }

    #[test]
    fn highlight_defaults_and_reset_restore_the_unstyled_scene() {
        let mut scene = bubble(false);
        let labels = typeset(&scene);
        let original = scene.render(&labels).unwrap();
        let palette = HighlightPalette::new("#ffd166");
        let region = HighlightRegion::new([1].into(), [0, 1].into(), HighlightStyle::Treated);
        scene.apply_highlights(slice::from_ref(&region), &palette, true);
        let single = scene.render(&labels).unwrap();
        assert!(!single.contains("data-linnet-hull"));
        assert!(single.contains("#e4504f") && single.contains("#77777773"));
        scene.apply_highlights_with_overlap(
            slice::from_ref(&region),
            &palette,
            true,
            HighlightOverlap::Recolour,
        );
        assert!(scene.highlights.recolour_stamp);
        scene.apply_highlights_with_overlap(
            slice::from_ref(&region),
            &palette,
            true,
            HighlightOverlap::Lanes,
        );
        assert!(!scene.highlights.recolour_stamp);
        assert!(!scene.highlights.lanes.is_empty());
        scene.apply_highlights(&[], &palette, true);
        assert_eq!(scene.render(&labels).unwrap(), original);

        let mut empty = HighlightRegion::new(
            Default::default(),
            Default::default(),
            HighlightStyle::Treated,
        );
        empty.details.contracted = true;
        scene.apply_highlights(&[empty], &palette, true);
        let svg = scene.render(&labels).unwrap();
        assert!(!svg.contains("data-linnet-hull"));
        assert!(svg.contains("data-linnet-counterterm"));
    }

    #[test]
    fn bands_are_offset_from_the_carrier_and_keep_node_arcs_separate() {
        let mut scene = bubble(false);
        let regions: Vec<_> = (1..=3)
            .map(|slot| {
                HighlightRegion::new([0, 1].into(), [0, 1].into(), HighlightStyle::Slot(slot))
            })
            .collect();
        let palette = HighlightPalette::new("#ffd166");
        let original = scene.edges[1].stroke.clone();
        scene.apply_highlights_with_overlap(
            &regions,
            &palette,
            false,
            HighlightOverlap::BandsMiniature,
        );
        assert_eq!(scene.edges[1].stroke, original);
        assert_eq!(
            scene.edges[1]
                .highlight_strokes
                .iter()
                .map(|stroke| stroke.offset)
                .collect::<Vec<_>>(),
            [Some(-4.0), Some(0.0), Some(4.0)]
        );
        let svg = scene.render(&typeset(&scene)).unwrap();
        let document = Document::parse(&svg).unwrap();
        let paths = |width: &str| {
            document
                .descendants()
                .filter(|element| {
                    element.has_tag_name("path") && element.attribute("stroke-width") == Some(width)
                })
                .filter_map(|element| element.attribute("d"))
                .collect::<Vec<_>>()
        };
        let strokes: Vec<_> = document
            .descendants()
            .filter(|element| element.has_tag_name("path"))
            .collect();
        let last_band = strokes
            .iter()
            .rposition(|element| element.attribute("stroke-width") == Some("2.4"))
            .unwrap();
        let first_particle = strokes
            .iter()
            .position(|element| element.attribute("stroke") == Some(original.paint.as_str()))
            .unwrap();
        assert!(
            last_band < first_particle,
            "all particle strokes must remain above the bands"
        );
        let bands = paths("2.4");
        assert_eq!(bands.len(), 6);
        for edge in bands.as_chunks::<3>().0 {
            assert_ne!(edge[0], edge[1]);
            assert_ne!(edge[1], edge[2]);
        }
        let rings = paths("2.5");
        assert_eq!(rings.len(), 6);
        assert!(rings
            .as_chunks::<3>()
            .0
            .iter()
            .all(|node| node[0] != node[1] && node[1] != node[2]));
        let palette = (1..=16).map(SlotPaint::new).collect::<Vec<_>>();
        assert_eq!(palette[0].light, "#0077bb");
        assert_eq!(
            palette
                .iter()
                .map(|paint| &paint.light)
                .collect::<BTreeSet<_>>()
                .len(),
            palette.len()
        );
        assert!(palette.iter().all(|paint| paint.light != paint.dark));
    }

    #[test]
    fn casings_keep_the_decorated_path_and_paint_outer_strokes_first() {
        let mut scene = bubble(false);
        let regions = [
            HighlightRegion {
                edges: [1, 2].into(),
                nodes: [0, 1].into(),
                style: HighlightStyle::ActedOn,
                details: HighlightDetails::default(),
            },
            HighlightRegion {
                edges: [1].into(),
                nodes: [0, 1].into(),
                style: HighlightStyle::Middle,
                details: HighlightDetails::default(),
            },
            HighlightRegion {
                edges: [1].into(),
                nodes: [0, 1].into(),
                style: HighlightStyle::Treated,
                details: HighlightDetails::default(),
            },
        ];
        let palette = HighlightPalette::new("#ffd166");
        scene.apply_highlights_with_overlap(&regions, &palette, false, HighlightOverlap::Casings);
        let expected = [
            HighlightStyle::ActedOn,
            HighlightStyle::Middle,
            HighlightStyle::Treated,
        ]
        .map(|style| style.light_paint().unwrap());
        let svg = scene.render(&typeset(&scene)).unwrap();
        let document = Document::parse(&svg).unwrap();
        let paths: Vec<_> = document
            .descendants()
            .filter(|element| element.has_tag_name("path"))
            .collect();
        let layered_colors: Vec<_> = paths
            .iter()
            .filter_map(|path| {
                path.attribute("stroke")
                    .filter(|paint| expected.contains(paint))
            })
            .collect();
        assert!(
            layered_colors.is_sorted_by_key(|paint| expected.iter().position(|role| role == paint))
        );
        let mut count = 0;
        for top in paths
            .iter()
            .filter(|path| path.attribute("stroke") == Some(expected[2]))
        {
            let casings: Vec<_> = paths
                .iter()
                .filter(|path| {
                    path.attribute("d") == top.attribute("d")
                        && path
                            .attribute("stroke")
                            .is_some_and(|paint| expected.contains(&paint))
                })
                .copied()
                .collect();
            let casings: [_; 3] = casings.try_into().unwrap();
            assert_eq!(
                casings.map(|path| path.attribute("stroke")),
                expected.map(Some)
            );
            assert!(casings
                .iter()
                .all(|path| path.attribute("d") == top.attribute("d")));
            let widths = casings.map(|path| {
                path.attribute("stroke-width")
                    .unwrap()
                    .parse::<f64>()
                    .unwrap()
            });
            assert!(widths.windows(2).all(|pair| pair[0] > pair[1]));
            count += 1;
        }
        assert_eq!(count, 1, "the decorated edge has one common path");
        assert_eq!(scene.edges[0].stroke, ink(1.0), "outside stays unmuted");
        assert_eq!(
            scene.nodes[0].stroke,
            ink(1.45),
            "named styles preserve nodes"
        );
    }

    #[test]
    fn layout_passes_continue_from_the_previous_geometry() {
        let mut scene = bubble(true);
        scene.layout.insert("impred-steps".into(), Value::from(3));
        let labels = typeset(&scene);
        let first = scene.render(&labels).unwrap();

        // A zero-step second pass must retain the solved geometry, not reseed it.
        scene.layout_passes.push(
            serde_json::from_value(serde_json::json!({
                "impred-steps": 0
            }))
            .unwrap(),
        );
        assert_eq!(scene.render(&labels).unwrap(), first);

        scene.layout_passes[0].insert("impred-steps".into(), Value::from(3));
        let refined = scene.render(&labels).unwrap();
        assert_ne!(refined, first);
        assert_eq!(refined, scene.render(&labels).unwrap());
    }

    #[test]
    fn a_later_layout_pass_reports_invalid_options() {
        let mut scene = bubble(false);
        scene.layout.insert("impred-steps".into(), Value::from(0));
        scene.layout_passes.push(
            serde_json::from_value(serde_json::json!({
                "impred-step-scale": 0
            }))
            .unwrap(),
        );
        assert!(scene
            .render(&typeset(&scene))
            .unwrap_err()
            .contains("step scale"));
    }

    #[test]
    fn an_impred_pass_preserves_a_preceding_layered_layout() {
        let mut scene = bubble(false);
        // A single carrier is valid input for either solver. Layered parallel
        // carriers can overlap, which ImPrEd correctly refuses to refine.
        let mut edge = scene.graph.edges.remove(1);
        edge.id = Some(0);
        scene.graph.edges = vec![edge];
        scene.edges = vec![scene.edges.remove(1)];
        scene
            .layout
            .insert("layout-algo".into(), Value::from("dot"));
        let labels = typeset(&scene);
        let layered = scene.render(&labels).unwrap();
        scene.layout_passes.push(
            serde_json::from_value(serde_json::json!({
                "layout-algo": "impred", "impred-steps": 0
            }))
            .unwrap(),
        );
        assert_eq!(scene.render(&labels).unwrap(), layered);
    }

    #[test]
    fn reads_title_and_measured_label_pages() {
        let scene = bubble(false);
        let typeset = typeset(&scene);
        assert_eq!(
            typeset.title.as_ref().unwrap().metrics,
            LabelMetrics {
                width: 10.0,
                height: 12.0,
                baseline: 8.0
            }
        );
        assert_eq!(
            typeset.pages[0].metrics,
            LabelMetrics {
                width: 10.0,
                height: 12.0,
                baseline: 8.0
            }
        );
        assert_eq!(typeset.pages[0].body, "<g/>");
        assert_eq!(typeset.defs.matches("<symbol").count(), 1);
        assert!(scene.typeset(&[page("<g/>")]).is_err());
        let document = scene.label_document().unwrap();
        assert!(document.starts_with("#set page(width: auto, height: auto, margin: 0pt"));
        assert!(document.contains("let body = \"bubble\";"));
        assert_eq!(document.matches("#context").count(), 2);
    }

    #[test]
    fn mathml_sources_and_measurements_roundtrip_without_vector_sources() {
        let mut scene = bubble(false);
        scene.title = Some(LabelSource::text("a<&\"b"));
        scene.pages[0] = LabelSource {
            mathml: "<mfrac><mi>α</mi><msub><mi>p</mi><mn>2</mn></msub></mfrac>".into(),
            typst: None,
        };
        let measurements = [
            LabelMetrics {
                width: 31.25,
                height: 9.5,
                baseline: 7.75,
            },
            LabelMetrics {
                width: 11.5,
                height: 18.25,
                baseline: 12.0,
            },
        ];
        let encoded = serde_json::to_string(&measurements).unwrap();
        let restored: Vec<LabelMetrics> = serde_json::from_str(&encoded).unwrap();
        assert_eq!(restored, measurements);
        let typeset = scene.typeset_mathml(&restored).unwrap();
        assert_eq!(typeset.title.as_ref().unwrap().metrics, measurements[0]);
        assert_eq!(typeset.pages[0].metrics, measurements[1]);
        assert!(typeset.defs.is_empty());
        assert!(matches!(
            scene.label_document(),
            Err(LabelError::MissingTypstSource { page: 1 })
        ));

        let svg = scene.render(&typeset).unwrap();
        let document = Document::parse(&svg).unwrap();
        let objects: Vec<_> = document
            .descendants()
            .filter(|node| node.has_tag_name(("http://www.w3.org/2000/svg", "foreignObject")))
            .collect();
        // The title and the four edge labels all embed the measured source.
        assert_eq!(objects.len(), 5);
        for (index, object) in objects.iter().enumerate() {
            let metrics = measurements[usize::from(index != 0)];
            assert_eq!(object.attribute("x"), Some("0"));
            assert_eq!(object.attribute("y"), Some("0"));
            assert_eq!(
                object.attribute("width").unwrap().parse::<f64>().unwrap(),
                metrics.width
            );
            assert_eq!(
                object.attribute("height").unwrap().parse::<f64>().unwrap(),
                metrics.height
            );
            let div = object.first_element_child().unwrap();
            assert!(div.has_tag_name(("http://www.w3.org/1999/xhtml", "div")));
            assert_eq!(
                div.attribute("style"),
                Some("font-size:9px;line-height:1;width:max-content")
            );
            let math = div.first_element_child().unwrap();
            assert!(math.has_tag_name(("http://www.w3.org/1998/Math/MathML", "math")));
            assert_eq!(
                math.first_element_child().unwrap().tag_name().name(),
                if index == 0 { "mtext" } else { "mfrac" }
            );
            if index == 0 {
                assert_eq!(math.first_element_child().unwrap().text(), Some("a<&\"b"));
            }
        }
        assert!(svg.contains("#linnet-node-0?"));
        assert!(svg.contains("#linnet-halfedge-1?"));
        assert!(svg.contains("&quot;particle&quot;:&quot;x&quot;"));
    }

    #[test]
    fn mathml_requires_exact_measurements_and_valid_fragments() {
        let mut scene = bubble(false);
        let metrics = LabelMetrics {
            width: 10.0,
            height: 12.0,
            baseline: 8.0,
        };
        for count in [0, 1, 3] {
            assert!(scene.typeset_mathml(&vec![metrics; count]).is_err());
        }
        for invalid in [
            LabelMetrics {
                width: f64::NAN,
                ..metrics
            },
            LabelMetrics {
                width: f64::INFINITY,
                ..metrics
            },
            LabelMetrics {
                width: -1.0,
                ..metrics
            },
            LabelMetrics {
                height: f64::NEG_INFINITY,
                ..metrics
            },
            LabelMetrics {
                height: -1.0,
                ..metrics
            },
            LabelMetrics {
                baseline: f64::NAN,
                ..metrics
            },
            LabelMetrics {
                baseline: -1.0,
                ..metrics
            },
            LabelMetrics {
                baseline: 13.0,
                ..metrics
            },
        ] {
            // Both the title and content boundary validate their measurements.
            assert!(scene.typeset_mathml(&[invalid, metrics]).is_err());
            assert!(scene.typeset_mathml(&[metrics, invalid]).is_err());
        }
        for source in [
            "<mfrac><mi>x</mi></mrow>",
            "<math><mi>x</mi></math>",
            "<div xmlns=\"http://www.w3.org/1999/xhtml\">x</div>",
        ] {
            scene.pages[0].mathml = source.into();
            assert!(scene.typeset_mathml(&[metrics, metrics]).is_err());
        }
        scene.pages[0] = LabelSource::text("");
        scene.title = None;
        let zero = LabelMetrics {
            width: 0.0,
            height: 0.0,
            baseline: 0.0,
        };
        assert_eq!(
            scene.typeset_mathml(&[zero]).unwrap().pages[0].metrics,
            zero
        );
        assert!(scene.typeset_mathml(&[zero, zero]).is_err());
        scene.pages.clear();
        let empty = scene.typeset_mathml(&[]).unwrap();
        assert!(empty.title.is_none() && empty.pages.is_empty());
        assert!(serde_json::from_str::<LabelMetrics>(r#"{"width":10,"height":12}"#).is_err());
    }

    #[test]
    fn native_mathml_and_vector_metrics_produce_the_same_geometry() {
        for momentum in [false, true] {
            let mut scene = bubble(momentum);
            scene.nodes[0].label = Some(0);
            let vector = typeset(&scene);
            let mathml = scene
                .typeset_mathml(&[
                    vector.title.as_ref().unwrap().metrics,
                    vector.pages[0].metrics,
                ])
                .unwrap();
            // Equal measured rectangles must feed the same node sizing,
            // collision search, arrow placement, targets, and page bounds.
            let mut geometry = mathml.clone();
            geometry.title.as_mut().unwrap().body = vector.title.as_ref().unwrap().body.clone();
            geometry.pages[0].body = vector.pages[0].body.clone();
            geometry.defs.clone_from(&vector.defs);
            assert_eq!(
                scene.render(&geometry).unwrap(),
                scene.render(&vector).unwrap()
            );
            assert_eq!(scene.nodes[0].size(&mathml), scene.nodes[0].size(&vector));
        }
    }

    #[test]
    fn vector_import_rejects_missing_and_invalid_measurements() {
        let scene = bubble(false);
        for metrics in [
            "NaN,8,8,12",
            "10,-8,8,12",
            "10,8,8,4",
            "10,8,8,inf",
            "10,8,8",
            "10,8,8,12,0",
        ] {
            let invalid = page(&format!(
                r#"<a href="linnest-metrics:{metrics}"><rect/></a>"#
            ));
            assert!(
                scene.typeset(&[invalid.clone(), invalid]).is_err(),
                "{metrics}"
            );
        }
        assert!(scene.typeset(&[page("<g/>"), page("<g/>")]).is_err());
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
        assert_eq!(svg.matches("data-linnet-momentum").count(), 4);
        for edge in 0..4 {
            assert!(svg.contains(&format!("<a href=\"#linnet-edge-{edge}?")));
        }
        // Edges without labels keep their arrows as fixed decorations.
        let mut unlabelled = bubble(true);
        for edge in &mut unlabelled.edges {
            edge.label = None;
        }
        let svg = unlabelled.render(&typeset(&unlabelled)).unwrap();
        assert_eq!(svg.matches(chevron).count(), 4);
    }

    #[test]
    fn momentum_clearance_moves_fixed_and_labelled_arrows_without_moving_edges() {
        let mut scene = bubble(true);
        let typeset = typeset(&scene);
        let graph = TypstGraph::from_spec(scene.graph.clone()).unwrap();
        let run = ImpredRun::seeded(graph, &Value::Object(scene.layout.clone())).unwrap();
        let laid = Laid::new(&run).unwrap();
        let before = Drawing::new(&scene, &typeset, &laid).unwrap();
        for edge in &mut scene.edges {
            edge.momentum_offset = 0.7;
        }
        let after = Drawing::new(&scene, &typeset, &laid).unwrap();
        let mut arrows = 0;
        for (before, after) in before.layers.iter().zip(&after.layers) {
            match (before, after) {
                (Element::Momentum { path: a, .. }, Element::Momentum { path: b, .. }) => {
                    assert_ne!(a, b);
                    arrows += 1;
                }
                (Element::Carrier { path: a, .. }, Element::Carrier { path: b, .. }) => {
                    assert_eq!(a, b);
                }
                _ => {}
            }
        }
        assert_eq!(arrows, 4);
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
        let laid = Laid::layered(
            &mut TypstGraph::from_spec(spec).unwrap(),
            &scene.layout,
            scene.layout_edges.as_ref(),
        )
        .unwrap();
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
            let (mut best, mut svg) = (drawing.score(), drawing.svg(typeset).unwrap());
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
                let rendered = drawing.svg(typeset).unwrap();
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
                typeset.pages[0].metrics.width = width;
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
