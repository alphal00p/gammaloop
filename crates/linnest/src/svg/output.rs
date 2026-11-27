//! SVG output: painted layers, typeset labels, and the `#linnet-*` inspection
//! links the interactive viewer decorates, as `draw.typ` emits them.
use std::fmt::Write as _;
use std::sync::Arc;

use kurbo::{BezPath, PathEl, Point, Rect, Shape};
use serde_json::Value;

use super::hulls::RESOLUTION_ERROR;
use super::{arrow_stroke, labels::Bounds, marks, Dash, Element, Stroke, Target, Typeset, UNIT};
use super::{CountertermScheme, HighlightDetails};
use crate::{TypstDotEdge, TypstDotEndpoint, TypstDotNode};

/// Canvas padding in drawing units; margin (2mm) and title gutter (1em) in points.
const PAD: f64 = 0.4;
const MARGIN: f64 = 5.669291339;
const GUTTER: f64 = 9.0;
/// Arrowhead paint and the particle-flow triangle's outline width in points.
pub(super) const INK: &str = "#3d2645";
pub(super) const MARK_STROKE: f64 = 0.3;
const HULL_STROKE: f64 = 1.1;

/// The ink weight of a drawing whose nominal canvas is `width` points wide:
/// recoloured stamps are drawn at reading size and scaled down, so their ink
/// grows with the canvas. The nominal canvas sets it before any bounds expansion.
fn ink_scale(width: f64, recolour_stamp: bool) -> f64 {
    if recolour_stamp {
        (width / 80.0).max(1.0)
    } else {
        1.0
    }
}

/// Outer radius of a counterterm glyph, in drawing units, on a drawing whose
/// nominal content is `content_width` drawing units wide.
pub(super) fn glyph_radius(content_width: f64, recolour_stamp: bool) -> f64 {
    let width = (content_width + 2.0 * PAD) * UNIT + 2.0 * MARGIN;
    GLYPH_EXTENT * ink_scale(width, recolour_stamp) / UNIT
}

/// Disc radius plus half its stroke, in points.
const GLYPH_EXTENT: f64 = 5.925;

// SVG's default miter limit is four times the stroke half-width.
const MITER_EXTENT: f64 = 2.0;

/// Ordered inspection fields shown by the interactive viewer.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Details(Vec<(String, Value)>);

impl Details {
    /// Ordered inspection properties, shared with graph binding snapshots.
    pub fn fields(&self) -> &[(String, Value)] {
        &self.0
    }

    /// Set a field, keeping its first position.
    pub fn insert(&mut self, key: &str, value: impl Into<Value>) {
        let value = value.into();
        match self.0.iter_mut().find(|(existing, _)| existing == key) {
            Some(slot) => slot.1 = value,
            None => self.0.push((key.to_owned(), value)),
        }
    }

    pub fn get(&self, key: &str) -> Option<&Value> {
        self.0
            .iter()
            .find(|(existing, _)| existing == key)
            .map(|(_, value)| value)
    }

    fn extended(&self, other: &Self) -> Self {
        let mut merged = self.clone();
        for (key, value) in &other.0 {
            merged.insert(key, value.clone());
        }
        merged
    }

    fn with(&self, key: &str, value: impl Into<Value>) -> Self {
        let mut copy = self.clone();
        copy.insert(key, value);
        copy
    }

    /// `#linnet-{kind}-{id}?{json}`, escaped for an attribute.
    fn href(&self, kind: &str, id: impl std::fmt::Display) -> String {
        let mut json = String::from("{");
        for (index, (key, value)) in self.0.iter().enumerate() {
            if index > 0 {
                json.push(',');
            }
            json.push_str(&Value::from(key.as_str()).to_string());
            json.push(':');
            json.push_str(&value.to_string());
        }
        json.push('}');
        xml_escape(&format!("#linnet-{kind}-{id}?{json}"))
    }
}

impl<K: Into<String>, V: Into<Value>> FromIterator<(K, V)> for Details {
    fn from_iter<I: IntoIterator<Item = (K, V)>>(iter: I) -> Self {
        let mut details = Self::default();
        for (key, value) in iter {
            details.insert(&key.into(), value);
        }
        details
    }
}

pub(super) fn xml_escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('"', "&quot;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

fn hover_details(details: &HighlightDetails) -> (String, String, String, String) {
    let escape = |text: &str| {
        let clean: String = text
            .chars()
            .map(|ch| match ch {
                '\u{0}'..='\u{8}'
                | '\u{b}'..='\u{c}'
                | '\u{e}'..='\u{1f}'
                | '\u{fffe}'
                | '\u{ffff}' => '\u{fffd}',
                _ => ch,
            })
            .collect();
        xml_escape(&clean)
            .replace('\t', "&#9;")
            .replace('\r', "&#13;")
            .replace('\n', "&#10;")
    };
    let label = escape(&details.label);
    let description = escape(&details.description);
    let mut title = match (label.is_empty(), description.is_empty()) {
        (true, _) => description.clone(),
        (_, true) => label.clone(),
        _ => format!("{label} — {description}"),
    };
    let scheme_attribute = if let Some(scheme) = details.scheme.name() {
        if !title.is_empty() {
            title.push_str(" — ");
        }
        title.push_str(scheme);
        format!(r#" data-scheme="{scheme}""#)
    } else {
        String::new()
    };
    (label, description, title, scheme_attribute)
}

/// Links of an edge: its four hover regions (source half-edge, source edge,
/// sink edge, sink half-edge) and its label.
pub(super) struct EdgeHrefs {
    pub regions: [Arc<str>; 4],
    pub label: String,
}

pub(super) fn edge_hrefs(edge: &TypstDotEdge, fields: &Details) -> EdgeHrefs {
    let node =
        |end: &Option<TypstDotEndpoint>| end.as_ref().map_or(Value::Null, |end| end.node.into());
    let details = [
        ("edge", Value::from(edge.edge)),
        ("source", node(&edge.source)),
        ("sink", node(&edge.sink)),
        ("orientation", edge.orientation.clone().into()),
    ]
    .into_iter()
    .collect::<Details>()
    .extended(fields);
    let mut halves = Vec::new();
    let mut edges = Vec::new();
    for (side, half, pair) in [
        ("source", &edge.source, &edge.sink),
        ("sink", &edge.sink, &edge.source),
    ] {
        let Some(owner) = half.as_ref().or(pair.as_ref()) else {
            continue;
        };
        let flow = match (half, side) {
            (Some(_), side) => side,
            (None, "source") => "sink",
            (None, _) => "source",
        };
        let other = if flow == "source" { "sink" } else { "source" };
        let hedge = details
            .get(&format!("{flow}-hedge"))
            .cloned()
            .unwrap_or_else(|| owner.hedge.into());
        let half_details = details.with("half-edge", hedge.clone());
        edges.push(
            half_details.href(
                "edge",
                details
                    .get("edge")
                    .cloned()
                    .unwrap_or_else(|| edge.edge.into()),
            ),
        );
        let pair_hedge = details
            .get(&format!("{other}-hedge"))
            .cloned()
            .unwrap_or_else(|| match (half, pair) {
                (Some(_), Some(pair)) => pair.hedge.into(),
                _ => Value::Null,
            });
        let half_details = half_details
            .with(
                "node",
                details
                    .get(flow)
                    .cloned()
                    .unwrap_or_else(|| owner.node.into()),
            )
            .with("flow", flow)
            .with("pair", pair_hedge);
        halves.push(half_details.href("halfedge", hedge));
    }
    let [source_half, sink_half] = <[String; 2]>::try_from(halves).unwrap_or_default();
    let [source_edge, sink_edge] = <[String; 2]>::try_from(edges).unwrap_or_default();
    EdgeHrefs {
        regions: [source_half, source_edge, sink_edge, sink_half].map(Into::into),
        label: details.href(
            "edge",
            details
                .get("edge")
                .cloned()
                .unwrap_or_else(|| edge.edge.into()),
        ),
    }
}

pub(super) fn node_href(node: &TypstDotNode, edges: &[TypstDotEdge], fields: &Details) -> String {
    let incident: Vec<Value> = edges
        .iter()
        .filter(|edge| {
            [&edge.source, &edge.sink]
                .into_iter()
                .flatten()
                .any(|end| end.node == node.node)
        })
        .map(|edge| edge.edge.into())
        .collect();
    std::iter::once(("edges", Value::from(incident)))
        .collect::<Details>()
        .extended(fields)
        .href(
            "node",
            fields
                .get("node")
                .cloned()
                .unwrap_or_else(|| node.node.into()),
        )
}

const SVG_DECIMALS: usize = 4;

fn number(value: f64) -> String {
    let text = format!("{value:.SVG_DECIMALS$}");
    let text = text.trim_end_matches('0').trim_end_matches('.');
    if text == "-0" {
        "0".into()
    } else {
        text.into()
    }
}

/// Validate the area of each closed contour in the coordinates actually
/// written to SVG. Subtract the first point to avoid cancellation on translation.
fn hull_path_data(
    path: &BezPath,
    coordinates: impl Fn(Point) -> (String, String),
) -> Result<String, String> {
    let mut data = String::new();
    let ticks = |text: &str| -> Result<i128, &'static str> {
        let value: f64 = text.parse().map_err(|_| RESOLUTION_ERROR)?;
        let value = (value * 10_i64.pow(SVG_DECIMALS as u32) as f64).round();
        if value.abs() < (1_i64 << 49) as f64 {
            Ok(value as i128)
        } else {
            Err(RESOLUTION_ERROR)
        }
    };
    let mut ring = Vec::new();
    for element in path.iter() {
        match element {
            PathEl::MoveTo(point) | PathEl::LineTo(point) => {
                let (x, y) = coordinates(point);
                let rounded = (ticks(&x)?, ticks(&y)?);
                ring.push(HullVertex { rounded, x, y });
            }
            PathEl::ClosePath => {
                write_hull_ring(&mut data, &ring)?;
                ring.clear();
            }
            _ => return Err("hull contour must be polygonal".into()),
        }
    }
    Ok(data)
}

struct HullVertex {
    rounded: (i128, i128),
    x: String,
    y: String,
}

fn write_hull_ring(data: &mut String, ring: &[HullVertex]) -> Result<(), &'static str> {
    let first = ring.first().ok_or(RESOLUTION_ERROR)?.rounded;
    let mut area = 0_i128;
    let mut command = 'M';
    for (index, vertex) in ring.iter().enumerate() {
        let count = ring.len();
        let previous = ring[(index + count - 1) % count].rounded;
        let current = vertex.rounded;
        let next = ring[(index + 1) % count].rounded;
        let incoming = (current.0 - previous.0, current.1 - previous.1);
        let outgoing = (next.0 - current.0, next.1 - current.1);
        area +=
            (current.0 - first.0) * (next.1 - first.1) - (current.1 - first.1) * (next.0 - first.0);
        // Only omit a point that lies on its neighbors' segment in emitted SVG
        // coordinates. A zero turn at a reversal is part of the contour.
        let cross = incoming.0 * outgoing.1 - incoming.1 * outgoing.0;
        let dot = incoming.0 * outgoing.0 + incoming.1 * outgoing.1;
        if cross != 0 || dot <= 0 {
            let _ = write!(data, "{command}{} {}", vertex.x, vertex.y);
            command = 'L';
        }
    }
    if area == 0 {
        return Err(RESOLUTION_ERROR);
    }
    data.push('Z');
    Ok(())
}

/// The page: an optional centered title above the padded drawing canvas.
pub(super) fn svg(
    typeset: &Typeset,
    layers: &[Element],
    targets: &[Target],
    recolour_stamp: bool,
) -> Result<String, String> {
    let mut bounds = Rect::new(
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    let mut include = |rect: Rect| bounds = bounds.union(rect);
    for element in layers {
        match element {
            Element::Carrier { .. } | Element::BeginGhost(_) | Element::EndGhost => {}
            Element::LaneCrossings(_) => {}
            Element::Hull { path, .. } | Element::Lane { path, .. } => include(path.bounding_box()),
            Element::LaneHalo { center, radius } => {
                let radius = radius + 1.6 / UNIT;
                include(Rect::from_points(*center, *center).inflate(radius, radius));
            }
            Element::Counterterm { center, .. } => {
                include(Rect::from_points(*center, *center).inflate(5.2 / UNIT, 5.2 / UNIT));
            }
            Element::Path {
                path,
                stroke,
                highlight,
            } => {
                let bounds = path.bounding_box();
                if *highlight {
                    let overflow = (stroke.width / (2.0 * UNIT) - PAD).max(0.0);
                    include(bounds.inflate(overflow, overflow));
                } else {
                    include(bounds);
                }
            }
            Element::Momentum { path, .. } => {
                include(path.bounding_box());
                if let Some(points) = marks::CHEVRON.at_end(path) {
                    for point in points {
                        include(Rect::from_points(point, point));
                    }
                }
            }
            Element::Triangle(points) => {
                for point in points {
                    include(Rect::from_points(*point, *point));
                }
            }
            Element::Label { bounds, .. } => include(Rect::from(*bounds)),
            Element::Node {
                at: [x, y],
                size: (w, h),
                ..
            } => include(Rect::new(x - w, y - h, x + w, y + h)),
        }
    }
    if !bounds.is_finite() {
        bounds = Rect::ZERO;
    }
    let title = typeset.title.as_ref();
    let title_width = title.map_or(0.0, |title| title.metrics.width);
    let title_height = title.map_or(0.0, |title| title.metrics.height + GUTTER);
    let width = ((bounds.width() + 2.0 * PAD) * UNIT).max(title_width) + 2.0 * MARGIN;
    let ink_scale = ink_scale(width, recolour_stamp);
    if recolour_stamp {
        let inflate = |rect: Rect, width: f64| rect.inflate(width / UNIT, width / UNIT);
        for element in layers {
            let painted = match element {
                Element::Node {
                    at: [x, y],
                    size: (w, h),
                    rectangular,
                    stroke,
                    ..
                } => {
                    let scale = if *rectangular { 1.0 } else { ink_scale };
                    let rect =
                        Rect::new(x - w * scale, y - h * scale, x + w * scale, y + h * scale);
                    inflate(rect, stroke.width * ink_scale / 2.0)
                }
                Element::Counterterm { center, .. } => inflate(
                    Rect::from_points(*center, *center),
                    GLYPH_EXTENT * ink_scale,
                ),
                Element::Path { path, stroke, .. } => {
                    inflate(path.bounding_box(), MITER_EXTENT * stroke.width * ink_scale)
                }
                Element::Momentum { path, .. } => {
                    let mut rect = inflate(
                        path.bounding_box(),
                        MITER_EXTENT * arrow_stroke().width * ink_scale,
                    );
                    if let Some(points) = marks::CHEVRON.at_end(path) {
                        let tip = points
                            .into_iter()
                            .fold(Rect::from_points(points[0], points[0]), |rect, point| {
                                rect.union_pt(point)
                            });
                        rect = rect.union(inflate(tip, MITER_EXTENT));
                    }
                    rect
                }
                Element::Triangle(points) => {
                    let triangle = points
                        .iter()
                        .copied()
                        .fold(Rect::from_points(points[0], points[0]), |rect, point| {
                            rect.union_pt(point)
                        });
                    inflate(triangle, MITER_EXTENT * MARK_STROKE)
                }
                Element::Hull { path, .. } => inflate(path.bounding_box(), HULL_STROKE / 2.0),
                Element::Lane { path, width, .. } => inflate(path.bounding_box(), width / 2.0),
                // The nominal bounds already include the unscaled halo, labels and nodes.
                _ => continue,
            };
            bounds = bounds.union(painted);
        }
    }
    let canvas_width = (bounds.width() + 2.0 * PAD) * UNIT;
    let canvas_height = (bounds.height() + 2.0 * PAD) * UNIT;
    let column = canvas_width.max(title_width);
    let (width, height) = (
        column + 2.0 * MARGIN,
        title_height + canvas_height + 2.0 * MARGIN,
    );
    let ox = MARGIN + (column - canvas_width) / 2.0;
    let oy = MARGIN + title_height;
    let x = |u: f64| number(ox + (u - bounds.x0 + PAD) * UNIT);
    let y = |v: f64| number(oy + (bounds.y1 + PAD - v) * UNIT);
    let mut svg = String::with_capacity(256 * 1024);
    let _ = write!(
        svg,
        r#"<svg viewBox="0 0 {w} {h}" data-linnet-renderer="native" width="{w}pt" height="{h}pt" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">"#,
        w = number(width),
        h = number(height)
    );
    if let Some(title) = title {
        let _ = write!(
            svg,
            r#"<g transform="translate({} {})">{}</g>"#,
            number(MARGIN + (column - title.metrics.width) / 2.0),
            number(MARGIN),
            title.body
        );
    }
    let path_data = |path: &BezPath| {
        let mut data = String::new();
        for element in path.elements() {
            let _ = match *element {
                PathEl::MoveTo(p) => write!(data, "M{} {}", x(p.x), y(p.y)),
                PathEl::LineTo(p) => write!(data, "L{} {}", x(p.x), y(p.y)),
                PathEl::QuadTo(c, p) => {
                    write!(data, "Q{} {} {} {}", x(c.x), y(c.y), x(p.x), y(p.y))
                }
                PathEl::CurveTo(c1, c2, p) => write!(
                    data,
                    "C{} {} {} {} {} {}",
                    x(c1.x),
                    y(c1.y),
                    x(c2.x),
                    y(c2.y),
                    x(p.x),
                    y(p.y)
                ),
                PathEl::ClosePath => write!(data, "Z"),
            };
        }
        data
    };
    let stroke_attributes = |stroke: &Stroke| {
        let mut attributes = format!(
            r#"stroke="{}" stroke-width="{}""#,
            xml_escape(&stroke.paint),
            number(stroke.width * ink_scale)
        );
        if stroke.round_cap {
            attributes.push_str(r#" stroke-linecap="round""#);
        }
        match stroke.dash {
            Dash::Solid => {}
            Dash::Dashed(on, off) => {
                let _ = write!(
                    attributes,
                    r#" stroke-dasharray="{} {}""#,
                    number(on),
                    number(off)
                );
            }
            Dash::Dotted => {
                let _ = write!(
                    attributes,
                    r#" stroke-dasharray="{} 2""#,
                    number(stroke.width)
                );
            }
        }
        attributes
    };
    let polyline = |points: &[kurbo::Point; 3]| {
        format!(
            "M{} {}L{} {}L{} {}",
            x(points[0].x),
            y(points[0].y),
            x(points[1].x),
            y(points[1].y),
            x(points[2].x),
            y(points[2].y)
        )
    };
    let has_lanes = layers
        .iter()
        .any(|element| matches!(element, Element::Lane { .. }));
    let write_targets = |svg: &mut String| {
        for Target { at, size, href } in targets {
            let _ = write!(
                svg,
                r#"<a href="{href}" transform="translate({} {})"><rect width="{s}" height="{s}" fill="transparent" stroke="none"/></a>"#,
                number(ox + (at[0] - bounds.x0 + PAD) * UNIT - size / 2.0),
                number(oy + (bounds.y1 + PAD - at[1]) * UNIT - size / 2.0),
                s = number(*size)
            );
        }
    };
    if has_lanes {
        write_targets(&mut svg);
    }
    for element in layers {
        match element {
            Element::Carrier { path, href } => {
                let _ = write!(
                    svg,
                    r#"<a href="{href}" data-linnet-carrier="true" style="display:none"><path d="{}"/></a>"#,
                    path_data(path)
                );
            }
            Element::Path { path, stroke, .. } => {
                let _ = write!(
                    svg,
                    r#"<path fill="none" {} d="{}"/>"#,
                    stroke_attributes(stroke),
                    path_data(path)
                );
            }
            Element::Momentum { path, href } => {
                let _ = write!(
                    svg,
                    r#"<a href="{href}" data-linnet-momentum="true"><path fill="none" {} d="{}"/>"#,
                    stroke_attributes(&arrow_stroke()),
                    path_data(path)
                );
                if let Some(points) = marks::CHEVRON.at_end(path) {
                    let _ = write!(
                        svg,
                        r#"<path fill="none" stroke="{INK}" stroke-width="1" stroke-linecap="round" stroke-linejoin="miter" d="{}"/>"#,
                        polyline(&points)
                    );
                }
                svg.push_str("</a>");
            }
            Element::BeginGhost(opacity) => {
                let _ = write!(svg, r#"<g data-linnet-ghost="true" opacity="{opacity}">"#);
            }
            Element::EndGhost => svg.push_str("</g>"),
            Element::Hull {
                path,
                paint,
                details,
            } => {
                let opacity = if details.contracted { 0.22 } else { 0.16 };
                let paint = xml_escape(paint);
                let d = hull_path_data(path, |point| (x(point.x), y(point.y)))?;
                let (label, description, title, scheme_attribute) = hover_details(details);
                let _ = write!(
                    svg,
                    r#"<g data-linnet-hull="true" data-label="{label}" data-description="{description}" data-contracted="{}"{scheme_attribute}><title>{title}</title><path d="{d}" fill="{paint}" fill-rule="evenodd" fill-opacity="{opacity}" stroke="{paint}" stroke-width="{HULL_STROKE}" stroke-linejoin="round"/></g>"#,
                    details.contracted
                );
            }
            Element::Lane {
                slot,
                path,
                paint,
                details,
                width,
            } => {
                let paint = xml_escape(paint);
                let (label, description, title, _) = hover_details(details);
                let _ = write!(
                    svg,
                    r#"<g data-linnet-lane="{slot}" data-label="{label}" data-description="{description}"><title>{title}</title><path d="{}" fill="none" stroke="{paint}" stroke-width="{width}" stroke-linecap="round" stroke-linejoin="round"/></g>"#,
                    path_data(path)
                );
            }
            Element::LaneCrossings(crossings) => {
                let _ = write!(
                    svg,
                    r#"<metadata data-linnet-lane-crossings="{crossings}"/>"#
                );
            }
            Element::LaneHalo { center, radius } => {
                let _ = write!(
                    svg,
                    r##"<circle cx="{}" cy="{}" r="{}" fill="none" stroke="#ffffff" stroke-width="1.6" data-linnet-paper="true" pointer-events="none"/>"##,
                    x(center.x),
                    y(center.y),
                    number(radius * UNIT + 0.8)
                );
            }
            Element::Counterterm {
                center,
                paint,
                details,
            } => {
                let paint = xml_escape(paint);
                let (label, description, title, scheme_attribute) = hover_details(details);
                let cx = x(center.x);
                let cy = y(center.y);
                let stamp_scale = if recolour_stamp {
                    format!(" scale({})", number(ink_scale))
                } else {
                    String::new()
                };
                // Filled on shell, open for MS-bar, and a dashed ring while the
                // scheme is unspecified.
                let (disc, cross, paper, ring) = match details.scheme {
                    CountertermScheme::OnShell => (
                        paint.as_str(),
                        "#ffffff",
                        r#" data-linnet-paper="true""#,
                        "",
                    ),
                    CountertermScheme::MsBar => ("#ffffff", paint.as_str(), "", ""),
                    CountertermScheme::Unspecified => (
                        "#ffffff",
                        paint.as_str(),
                        "",
                        r#" stroke-dasharray="2.2 1.6""#,
                    ),
                };
                let _ = write!(
                    svg,
                    r##"<g data-linnet-counterterm="true" data-label="{label}" data-description="{description}"{scheme_attribute} transform="translate({cx} {cy}){stamp_scale}"><title>{title}</title><circle r="5.2" fill="{disc}" stroke="{paint}" stroke-width="1.45"{ring}/><path d="M-3.65 -3.65L3.65 3.65M-3.65 3.65L3.65 -3.65" fill="none" stroke="{cross}" stroke-width="1.45"{paper}/></g>"##
                );
            }
            Element::Triangle(points) => {
                let _ = write!(
                    svg,
                    r#"<path fill="{INK}" stroke="{INK}" stroke-width="{MARK_STROKE}" stroke-linejoin="miter" d="{}Z"/>"#,
                    polyline(points)
                );
            }
            Element::Label {
                page,
                bounds: b,
                href,
            } => {
                let page = &typeset.pages[*page];
                let (lx, ly) = (x(b.left), y(b.top));
                let _ = write!(
                    svg,
                    r#"<g transform="translate({lx} {ly})">{body}</g>"#,
                    body = page.body,
                );
                if let Some(href) = href {
                    let _ = write!(
                        svg,
                        r#"<a href="{href}" transform="translate({lx} {ly})"><rect width="{w}" height="{h}" fill="transparent" stroke="none"/></a>"#,
                        w = number(page.metrics.width),
                        h = number(page.metrics.height)
                    );
                }
            }
            Element::Node {
                id,
                at: [nx, ny],
                size: (w, h),
                rectangular,
                fill,
                stroke,
            } => {
                if *rectangular {
                    let _ = write!(
                        svg,
                        r#"<rect data-linnet-node="{id}" x="{}" y="{}" width="{}" height="{}" rx="2" fill="{}" {}/>"#,
                        x(nx - w),
                        y(ny + h),
                        number(2.0 * w * UNIT),
                        number(2.0 * h * UNIT),
                        xml_escape(fill),
                        stroke_attributes(stroke)
                    );
                } else {
                    let _ = write!(
                        svg,
                        r#"<circle data-linnet-node="{id}" cx="{}" cy="{}" r="{}" fill="{}" {}/>"#,
                        x(*nx),
                        y(*ny),
                        number(w * UNIT * ink_scale),
                        xml_escape(fill),
                        stroke_attributes(stroke)
                    );
                }
            }
        }
    }
    if !has_lanes {
        write_targets(&mut svg);
    }
    let _ = write!(
        svg,
        r##"<defs><pattern id="linnest-process-hatch" width="5" height="5" patternUnits="userSpaceOnUse"><path d="M0 5L5 0" stroke="#3d2645" stroke-width="0.35"/></pattern>{}</defs></svg>"##,
        typeset.defs
    );
    Ok(svg)
}

impl From<Bounds> for Rect {
    fn from(bounds: Bounds) -> Self {
        Rect::new(bounds.left, bounds.bottom, bounds.right, bounds.top)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hull_serialization_only_omits_points_collinear_after_rounding() {
        let mut path = BezPath::new();
        path.move_to((0.0, 0.0));
        path.line_to((0.1, 0.10004));
        path.line_to((0.2, 0.20008));
        path.line_to((0.3, 0.30012));
        path.line_to((0.4, 0.30012));
        path.line_to((0.5, 0.30012));
        path.close_path();
        assert_eq!(
            hull_path_data(&path, |p| (number(p.x), number(p.y))).unwrap(),
            "M0 0L0.1 0.1L0.2 0.2001L0.3 0.3001L0.5 0.3001Z"
        );
    }

    #[test]
    fn hull_serialization_rejects_a_contour_that_rounds_to_zero_area() {
        let mut path = BezPath::new();
        path.move_to((44.4829, 58.8441));
        path.line_to((44.482901, 58.8441));
        path.line_to((44.4829, 58.844101));
        path.close_path();
        assert_eq!(
            hull_path_data(&path, |p| (number(p.x), number(p.y))).unwrap_err(),
            RESOLUTION_ERROR
        );
    }
}
