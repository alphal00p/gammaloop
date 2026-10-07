//! SVG output: painted layers, typeset labels, and the `#linnet-*` inspection
//! links the interactive viewer decorates, as `draw.typ` emits them.
use std::fmt::Write as _;
use std::sync::Arc;

use kurbo::{BezPath, PathEl, Rect, Shape};
use serde_json::Value;

use super::{labels::Bounds, Dash, Element, Stroke, Target, Typeset, UNIT};
use crate::{TypstDotEdge, TypstDotEndpoint, TypstDotNode};

/// Canvas padding in drawing units; margin (2mm) and title gutter (1em) in points.
const PAD: f64 = 0.4;
const MARGIN: f64 = 5.669291339;
const GUTTER: f64 = 9.0;
/// Arrowhead paint and the particle-flow triangle's outline width in points.
pub(super) const INK: &str = "#3d2645";
pub(super) const MARK_STROKE: f64 = 0.3;

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

fn xml_escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('"', "&quot;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
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

fn number(value: f64) -> String {
    let text = format!("{value:.4}");
    let text = text.trim_end_matches('0').trim_end_matches('.');
    if text == "-0" {
        "0".into()
    } else {
        text.into()
    }
}

/// The page: an optional centered title above the padded drawing canvas.
pub(super) fn svg(typeset: &Typeset, layers: &[Element], targets: &[Target]) -> String {
    let mut bounds = Rect::new(
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    let mut include = |rect: Rect| bounds = bounds.union(rect);
    for element in layers {
        match element {
            Element::Path { path, .. } => include(path.bounding_box()),
            Element::Chevron(points) | Element::Triangle(points) => {
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
    let canvas_width = (bounds.width() + 2.0 * PAD) * UNIT;
    let canvas_height = (bounds.height() + 2.0 * PAD) * UNIT;
    let title = typeset.title.as_ref();
    let column = title.map_or(canvas_width, |title| canvas_width.max(title.width));
    let title_height = title.map_or(0.0, |title| title.height + GUTTER);
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
            number(MARGIN + (column - title.width) / 2.0),
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
            number(stroke.width)
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
    for element in layers {
        match element {
            Element::Path { path, stroke } => {
                let _ = write!(
                    svg,
                    r#"<path fill="none" {} d="{}"/>"#,
                    stroke_attributes(stroke),
                    path_data(path)
                );
            }
            Element::Chevron(points) => {
                let _ = write!(
                    svg,
                    r#"<path fill="none" stroke="{INK}" stroke-width="1" stroke-linecap="round" stroke-linejoin="miter" d="{}"/>"#,
                    polyline(points)
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
                        w = number(page.width),
                        h = number(page.height)
                    );
                }
            }
            Element::Node {
                at: [nx, ny],
                size: (w, h),
                rectangular,
                fill,
                stroke,
            } => {
                if *rectangular {
                    let _ = write!(
                        svg,
                        r#"<rect x="{}" y="{}" width="{}" height="{}" rx="2" fill="{}" {}/>"#,
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
                        r#"<circle cx="{}" cy="{}" r="{}" fill="{}" {}/>"#,
                        x(*nx),
                        y(*ny),
                        number(w * UNIT),
                        xml_escape(fill),
                        stroke_attributes(stroke)
                    );
                }
            }
        }
    }
    for Target { at, size, href } in targets {
        let _ = write!(
            svg,
            r#"<a href="{href}" transform="translate({} {})"><rect width="{s}" height="{s}" fill="transparent" stroke="none"/></a>"#,
            number(ox + (at[0] - bounds.x0 + PAD) * UNIT - size / 2.0),
            number(oy + (bounds.y1 + PAD - at[1]) * UNIT - size / 2.0),
            s = number(*size)
        );
    }
    let _ = write!(
        svg,
        r##"<defs><pattern id="linnest-process-hatch" width="5" height="5" patternUnits="userSpaceOnUse"><path d="M0 5L5 0" stroke="#3d2645" stroke-width="0.35"/></pattern>{}</defs></svg>"##,
        typeset.defs
    );
    svg
}

impl From<Bounds> for Rect {
    fn from(bounds: Bounds) -> Self {
        Rect::new(bounds.left, bounds.bottom, bounds.right, bounds.top)
    }
}
