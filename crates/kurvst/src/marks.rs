//! Data-only mark geometry. Every head is built in a local frame whose origin
//! is the tip and whose positive x axis points along the line, which arrives
//! from negative x. A shape's paint, stroked or not, keeps its nominal length
//! and width unless the head is too small to hold its stroke, when it grows
//! just enough. The shaft stops where the drawn path meets the line at the
//! rear. Coordinates are drawing units; angles are radians. Paint stays with
//! the caller. A combined mark is prepared locally and fitted once.

use self::HeadShape::{
    Bar, Barb, Bracket, Circle, Diamond, Hooks, Rays, Round, Square, Stealth, Straight, Tikz,
    Triangle,
};
use crate::curve_api::{
    CurvePathOutput, CurvePoint, PathTrimmer, deserialize_f64, parse_stroke_cap, parse_stroke_join,
};
use kurbo::{
    Affine, Arc, BezPath, Cap, Ellipse, Join, ParamCurve, PathEl, Point, Shape, Stroke, StrokeOpts,
    Vec2,
};
use serde::{Deserialize, Serialize};
use std::f64::consts::{FRAC_PI_2, PI, SQRT_2, TAU};
use std::iter;

/// Mixed size components: physical points plus a dimensionless stroke ratio.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MarkSize {
    #[serde(deserialize_with = "deserialize_f64")]
    pub points: f64,
    #[serde(deserialize_with = "deserialize_f64")]
    pub ratio: f64,
}

impl MarkSize {
    pub fn resolve(self, context: MarkContext) -> Result<f64, String> {
        context.validate()?;
        if !self.points.is_finite() || !self.ratio.is_finite() {
            return Err("mark size components must be finite".into());
        }
        let value = self.points * context.units_per_pt + self.ratio * context.line_thickness;
        if !value.is_finite() {
            return Err("resolved mark size must be finite".into());
        }
        Ok(value)
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct MarkContext {
    #[serde(deserialize_with = "deserialize_f64")]
    pub units_per_pt: f64,
    #[serde(deserialize_with = "deserialize_f64")]
    pub line_thickness: f64,
    #[serde(default)]
    pub shaft_stroke: MarkStrokeStyle,
}

/// Normalized, solid shaft geometry. Paint, dash patterns, and callbacks are
/// deliberately absent. A solid outline is conservative for a dashed shaft.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(default, rename_all = "kebab-case", deny_unknown_fields)]
pub struct MarkStrokeStyle {
    #[serde(
        serialize_with = "MarkStrokeStyle::serialize_cap",
        deserialize_with = "MarkStrokeStyle::deserialize_cap"
    )]
    pub cap: Cap,
    #[serde(
        serialize_with = "MarkStrokeStyle::serialize_join",
        deserialize_with = "MarkStrokeStyle::deserialize_join"
    )]
    pub join: Join,
    #[serde(deserialize_with = "deserialize_f64")]
    pub miter_limit: f64,
}

impl Default for MarkStrokeStyle {
    fn default() -> Self {
        Self {
            cap: Cap::Butt,
            join: Join::Miter,
            miter_limit: 4.0,
        }
    }
}

impl MarkStrokeStyle {
    fn deserialize_cap<'de, D: serde::Deserializer<'de>>(deserializer: D) -> Result<Cap, D::Error> {
        let value = String::deserialize(deserializer)?;
        parse_stroke_cap(&value).map_err(serde::de::Error::custom)
    }

    fn deserialize_join<'de, D: serde::Deserializer<'de>>(
        deserializer: D,
    ) -> Result<Join, D::Error> {
        let value = String::deserialize(deserializer)?;
        parse_stroke_join(&value).map_err(serde::de::Error::custom)
    }

    fn serialize_cap<S: serde::Serializer>(cap: &Cap, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(match cap {
            Cap::Butt => "butt",
            Cap::Round => "round",
            Cap::Square => "square",
        })
    }

    fn serialize_join<S: serde::Serializer>(join: &Join, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(match join {
            Join::Miter => "miter",
            Join::Round => "round",
            Join::Bevel => "bevel",
        })
    }

    fn stroke(self, width: f64) -> Stroke {
        Stroke::new(width)
            .with_caps(self.cap)
            .with_join(self.join)
            .with_miter_limit(self.miter_limit)
    }
}

impl MarkContext {
    fn validate(self) -> Result<(), String> {
        if !self.units_per_pt.is_finite()
            || self.units_per_pt <= 0.0
            || !self.line_thickness.is_finite()
            || self.line_thickness < 0.0
        {
            return Err("context requires positive finite units-per-pt and nonnegative finite line-thickness".into());
        }
        if !self.shaft_stroke.miter_limit.is_finite() || self.shaft_stroke.miter_limit <= 0.0 {
            return Err("shaft stroke miter-limit must be finite and positive".into());
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum MarkFit {
    #[default]
    Chord,
    Bend,
}

#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum MarkAlign {
    #[default]
    Center,
    End,
}

#[derive(Debug, Clone, Serialize)]
#[serde(untagged)]
pub enum MarkPart {
    Mark(Box<MarkSpec>),
    Gap(MarkGap),
}

impl<'de> Deserialize<'de> for MarkPart {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        // Dispatch by the data discriminator rather than trying alternatives:
        // untagged trial deserialization erases useful nested shape errors.
        let value = ciborium::value::Value::deserialize(deserializer)?;
        let has_shape = value
            .as_map()
            .is_some_and(|fields| fields.iter().any(|(key, _)| key.as_text() == Some("shape")));
        if has_shape {
            value
                .deserialized::<MarkSpec>()
                .map(|mark| Self::Mark(Box::new(mark)))
        } else {
            value.deserialized::<MarkGap>().map(Self::Gap)
        }
        .map_err(serde::de::Error::custom)
    }
}

/// A signed size gap between successive composite children.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MarkGap {
    pub gap: MarkSize,
}

// Deriving serde against Self keeps one field inventory while the trait
// implementation below additionally enforces shape-specific data rules.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(remote = "Self", deny_unknown_fields)]
pub struct MarkSpec {
    pub shape: String,
    /// Omission resolves to chord. Composite children must omit this field,
    /// even if an explicit value would equal the default.
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "MarkSpec::present"
    )]
    pub fit: Option<MarkFit>,
    /// Omission resolves to 1. Composite shortening belongs to the root only.
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "MarkSpec::present_number"
    )]
    pub shorten: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "MarkSpec::present"
    )]
    pub length: Option<MarkSize>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "MarkSpec::present"
    )]
    pub width: Option<MarkSize>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "MarkSpec::present_number"
    )]
    pub inset: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "MarkSpec::present_number"
    )]
    pub arc: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "MarkSpec::present_number"
    )]
    pub phase: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "MarkSpec::present"
    )]
    pub n: Option<usize>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "MarkSpec::present"
    )]
    pub rev: Option<bool>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "MarkSpec::present"
    )]
    pub align: Option<MarkAlign>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "MarkSpec::present"
    )]
    pub fill: Option<bool>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "MarkSpec::present"
    )]
    pub stroke: Option<bool>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "MarkSpec::present"
    )]
    pub parts: Option<Vec<MarkPart>>,
}

impl Serialize for MarkSpec {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        Self::serialize(self, serializer)
    }
}

impl<'de> Deserialize<'de> for MarkSpec {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let spec = Self::deserialize(deserializer)?;
        spec.validate_fields().map_err(serde::de::Error::custom)?;
        Ok(spec)
    }
}

/// Stations refer to the caller's visible carrier, before mark fitting.
/// Ratios are dimensionless; distances and signed placement shifts are in
/// drawing units. Numeric stations locate the painted longitudinal center and
/// always use chord fit, even at an endpoint. Start and End locate the template
/// origin/contact and retain the requested endpoint fit.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "lowercase", deny_unknown_fields)]
pub enum MarkStation {
    Start,
    End,
    Ratio {
        #[serde(deserialize_with = "deserialize_f64")]
        value: f64,
    },
    Distance {
        #[serde(deserialize_with = "deserialize_f64")]
        value: f64,
    },
}

#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum MarkDirection {
    #[default]
    Forward,
    Backward,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MarkPlacement {
    pub template: usize,
    pub carrier: usize,
    pub station: MarkStation,
    #[serde(default, deserialize_with = "deserialize_f64")]
    pub shift: f64,
    #[serde(default)]
    pub direction: MarkDirection,
}

/// Each template owns the context that resolves its dimensions. Different
/// carrier stroke widths can therefore share a single bounded geometry batch.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MarkTemplateSpec {
    pub mark: MarkSpec,
    pub context: MarkContext,
}

/// Carriers must have one subpath. Endpoint bend on a closed carrier is
/// explicitly unsupported; numeric interior stations still use chord fit.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MarkGeometrySpec {
    #[serde(default)]
    pub mode: MarkGeometryMode,
    pub templates: Vec<MarkTemplateSpec>,
    pub carriers: Vec<CurvePathOutput>,
    pub placements: Vec<MarkPlacement>,
}

#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum MarkGeometryMode {
    #[default]
    Candidates,
    Selected,
}

/// Candidates have independent per-mark painted shafts and no grouped shafts.
/// Selected marks are head-only; paint their paths and grouped shafts once.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MarkGeometryBatchOutput {
    pub marks: Vec<MarkGeometryOutput>,
    pub shafts: Vec<MarkShaftOutput>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub struct MarkShaftOutput {
    pub carrier: usize,
    pub shaft: CurvePathOutput,
    pub shaft_outline: CurvePathOutput,
    pub shaft_style: MarkStrokeStyle,
    pub footprint: CurvePathOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub struct MarkDrawable {
    pub path: CurvePathOutput,
    pub outline: CurvePathOutput,
    pub outline_bounds: Vec<MarkBounds>,
    pub fill: bool,
    pub stroke: bool,
    pub join: String,
    pub cap: String,
    pub miter_limit: f64,
}

/// `back` is where the full-size mark stops the shaft, `end` behind the tip on
/// the line through the head. At shorten=1 it coincides with `shaft-contact`
/// on an ordinary carrier. If a carrier is too short to fit the requested head,
/// only the shaft contact is clamped inward;
/// the mark is not silently rescaled. Empty carriers use their move point (or
/// the origin) and a deterministic horizontal orientation.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub struct MarkGeometryOutput {
    pub template: usize,
    pub carrier: usize,
    pub end: f64,
    pub paths: Vec<MarkDrawable>,
    pub tip: CurvePoint,
    pub back: CurvePoint,
    pub shaft_contact: CurvePoint,
    pub shaft: CurvePathOutput,
    pub shaft_outline: CurvePathOutput,
    pub shaft_style: MarkStrokeStyle,
    pub footprint: CurvePathOutput,
    pub footprint_bounds: Vec<MarkBounds>,
}

/// Ordered conservative control hulls, one per nonempty subpath.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct MarkBounds {
    pub left: f64,
    pub right: f64,
    pub bottom: f64,
    pub top: f64,
}

impl CurvePathOutput {
    /// Match the wire element traversal: moves flush a region; endpoints precede
    /// controls. Strict comparisons retain the first equal coordinate (including
    /// signed zero), as Typst's min/max do. Close elements have no coordinates.
    fn control_hull_bounds(&self) -> Vec<MarkBounds> {
        let mut regions = Vec::new();
        let mut current: Option<MarkBounds> = None;
        for element in self.path.elements() {
            if matches!(element, PathEl::MoveTo(_))
                && let Some(bounds) = current.take()
            {
                regions.push(bounds);
            }
            let points: &[Point] = match element {
                PathEl::MoveTo(point) | PathEl::LineTo(point) => std::slice::from_ref(point),
                PathEl::QuadTo(control, end) => &[*end, *control],
                PathEl::CurveTo(first, second, end) => &[*end, *first, *second],
                PathEl::ClosePath => &[],
            };
            for point in points {
                if let Some(bounds) = &mut current {
                    if point.x < bounds.left {
                        bounds.left = point.x;
                    }
                    if point.x > bounds.right {
                        bounds.right = point.x;
                    }
                    if point.y < bounds.bottom {
                        bounds.bottom = point.y;
                    }
                    if point.y > bounds.top {
                        bounds.top = point.y;
                    }
                } else {
                    current = Some(MarkBounds {
                        left: point.x,
                        right: point.x,
                        bottom: point.y,
                        top: point.y,
                    });
                }
            }
        }
        regions.extend(current);
        regions
    }
}

pub struct PreparedMark {
    paths: Vec<MarkDrawable>,
    longitudinal_center: f64,
    end: f64,
    fit: MarkFit,
    shorten: f64,
    line_thickness: f64,
    shaft_style: MarkStrokeStyle,
}

struct PreparedCarrier {
    path: BezPath,
    trimmer: PathTrimmer,
    length: f64,
}

#[derive(Default)]
struct ShaftWindow {
    start_outset: f64,
    end_outset: f64,
    start_bend: Option<Point>,
    end_bend: Option<Point>,
}

enum ShaftEdit {
    None,
    TrimStart(f64),
    TrimEnd(f64),
    BendStart(Point),
    BendEnd(Point),
}

impl ShaftWindow {
    fn apply(&mut self, edit: &ShaftEdit) {
        match edit {
            ShaftEdit::None => {}
            ShaftEdit::TrimStart(value) => self.start_outset = *value,
            ShaftEdit::TrimEnd(value) => self.end_outset = *value,
            ShaftEdit::BendStart(point) => self.start_bend = Some(*point),
            ShaftEdit::BendEnd(point) => self.end_bend = Some(*point),
        }
    }
}

struct FittedMark {
    geometry: MarkGeometryOutput,
    shaft_edit: ShaftEdit,
}

struct SelectedShaft {
    template: usize,
    window: ShaftWindow,
    start_selected: bool,
    end_selected: bool,
}

impl PreparedCarrier {
    fn shaft(&self, window: &ShaftWindow) -> Result<BezPath, String> {
        // Fit all stations on the untouched carrier. Apply its combined chord
        // window once, then translate only bent endpoints; retained controls
        // are not changed by bending. Overlapping chord windows leave no shaft.
        if self.length > 0.0 && window.start_outset + window.end_outset >= self.length {
            return Ok(BezPath::new());
        }
        let path = if window.start_outset == 0.0 && window.end_outset == 0.0 {
            self.path.clone()
        } else {
            BezPath::from_path_segments(
                self.trimmer
                    .trim(window.start_outset, window.end_outset)?
                    .into_iter(),
            )
        };
        let mut elements = path.elements().to_vec();
        if let Some(point) = window.start_bend
            && let Some(PathEl::MoveTo(start)) = elements.first_mut()
        {
            *start = point;
        }
        if let Some(point) = window.end_bend
            && let Some(PathEl::LineTo(end) | PathEl::QuadTo(_, end) | PathEl::CurveTo(_, _, end)) =
                elements.last_mut()
        {
            *end = point;
        }
        Ok(BezPath::from_vec(elements))
    }

    fn new(path: BezPath) -> Result<Self, String> {
        if !path.is_finite() {
            return Err("carrier coordinates must be finite".into());
        }
        if path.subpaths().count() > 1 {
            let mut subpaths = path.subpaths();
            let first = BezPath::from_vec(subpaths.next().unwrap().to_vec());
            let start = match subpaths.next().unwrap()[0] {
                PathEl::MoveTo(point) => point,
                _ => unreachable!("a subpath starts with a move"),
            };
            let end = first.segments().last().map(|segment| segment.end());
            return Err(format!(
                "mark carriers must have one connected subpath; next start {start:?}, previous end {end:?}"
            ));
        }
        let trimmer = PathTrimmer::new(path.segments(), 1e-6);
        let length = trimmer.length();
        if !length.is_finite() {
            return Err("carrier length must be finite".into());
        }
        Ok(Self {
            path,
            trimmer,
            length,
        })
    }

    fn point(&self, station: f64) -> Result<Point, String> {
        Ok(self
            .trimmer
            .frame(station)?
            .map(|frame| frame.point.into())
            .unwrap_or_else(|| {
                self.path
                    .elements()
                    .iter()
                    .find_map(|el| match el {
                        PathEl::MoveTo(point) => Some(*point),
                        _ => None,
                    })
                    .unwrap_or(Point::ZERO)
            }))
    }

    /// Reuse the carrier owner's isolated and refined geometric contact.
    /// When a carrier is too short, retain full-sized head geometry, clamp its
    /// shaft contact inward, and use the available chord for orientation.
    fn chord_contact(&self, station: f64, end: f64, sign: f64) -> Result<(Point, f64), String> {
        let travel_sign = if end < 0.0 { sign } else { -sign };
        if let Some((point, at)) =
            self.trimmer
                .chord_contact(station, end.abs(), travel_sign > 0.0)?
        {
            return Ok((point.into(), at));
        }
        let at = if travel_sign > 0.0 { self.length } else { 0.0 };
        Ok((self.point(at)?, at))
    }
}

impl MarkSpec {
    fn present_number<'de, D: serde::Deserializer<'de>>(
        deserializer: D,
    ) -> Result<Option<f64>, D::Error> {
        deserialize_f64(deserializer).map(Some)
    }

    // Unlike Option's normal wire behavior, explicit null is not omission.
    fn present<'de, D: serde::Deserializer<'de>, T: Deserialize<'de>>(
        deserializer: D,
    ) -> Result<Option<T>, D::Error> {
        T::deserialize(deserializer).map(Some)
    }

    /// Check the fields against the shape, which is returned parsed.
    fn validate_fields(&self) -> Result<MarkShape, String> {
        let shape = MarkShape::parse(&self.shape)?;
        for (name, present) in [
            ("length", self.length.is_some()),
            ("width", self.width.is_some()),
            ("inset", self.inset.is_some()),
            ("arc", self.arc.is_some()),
            ("phase", self.phase.is_some()),
            ("n", self.n.is_some()),
            ("rev", self.rev.is_some()),
            ("align", self.align.is_some()),
            ("fill", self.fill.is_some()),
            ("stroke", self.stroke.is_some()),
            ("parts", self.parts.is_some()),
        ] {
            if present && !shape.options().contains(&name) {
                return Err(format!("{}.{name} is not supported", self.shape));
            }
        }
        if self
            .shorten
            .is_some_and(|value| !value.is_finite() || value < 0.0)
        {
            return Err("shorten must be finite and nonnegative".into());
        }
        for size in [self.length, self.width].into_iter().flatten() {
            if !size.points.is_finite() || !size.ratio.is_finite() {
                return Err("mark size components must be finite".into());
            }
        }
        if self
            .inset
            .is_some_and(|value| !value.is_finite() || !(0.0..=1.0).contains(&value))
        {
            return Err("inset must be finite in [0,1]".into());
        }
        if self
            .arc
            .is_some_and(|value| !value.is_finite() || value.abs() > TAU)
        {
            return Err("arc must be finite with magnitude at most 2π".into());
        }
        if self.phase.is_some_and(|value| !value.is_finite()) {
            return Err("phase must be finite".into());
        }
        if self.n.is_some_and(|value| value == 0 || value > 4096) {
            return Err("rays.n must be in 1..=4096".into());
        }
        if shape == MarkShape::Combine && self.parts.as_deref().unwrap_or_default().is_empty() {
            return Err("combine.parts must be nonempty".into());
        }
        if self.parts.as_deref().unwrap_or_default().iter().any(|part| {
            matches!(part,MarkPart::Mark(child) if child.fit.is_some() || child.shorten.is_some())
        }) {
            return Err("composite children must omit fit and shorten, including their explicit defaults; set them on the composite".into());
        }
        Ok(shape)
    }

    /// Resolve dimensions and build the painted local geometry once.
    pub fn prepare(&self, context: MarkContext) -> Result<PreparedMark, String> {
        context.validate()?;
        self.prepare_nested(context, 0)
    }

    fn prepare_nested(&self, context: MarkContext, depth: usize) -> Result<PreparedMark, String> {
        if depth > 64 {
            return Err("mark nesting exceeds 64 levels".into());
        }
        let shape = self.validate_fields()?;
        let (paths, end) = match shape {
            MarkShape::Combine => {
                // Children follow one another down the shaft: each starts
                // where the previous one's shaft stop (plus any signed gap) lies.
                let mut paths = Vec::new();
                let mut cursor = 0.0;
                for part in self.parts.as_deref().unwrap_or_default() {
                    match part {
                        MarkPart::Gap(gap) => cursor += gap.gap.resolve(context)?,
                        MarkPart::Mark(spec) => {
                            let child = spec.prepare_nested(context, depth + 1)?;
                            let shift = Affine::translate((-cursor, 0.0));
                            for mut drawable in child.paths {
                                drawable.path.path.apply_affine(shift);
                                drawable.outline.path.apply_affine(shift);
                                paths.push(drawable);
                            }
                            cursor += child.end;
                        }
                    }
                    if !cursor.is_finite() {
                        return Err("composite contact must be finite".into());
                    }
                }
                (paths, cursor)
            }
            MarkShape::Head(shape) => {
                let (drawable, end) = self.leaf(shape, context)?;
                (vec![drawable], end)
            }
        };
        Ok(PreparedMark {
            longitudinal_center: PreparedMark::longitudinal_center(&paths),
            paths,
            end,
            fit: self.fit.unwrap_or_default(),
            shorten: self.shorten.unwrap_or(1.0),
            line_thickness: context.line_thickness,
            shaft_style: context.shaft_stroke,
        })
    }

    /// One leaf head: its drawn path in the local frame, turned around for
    /// `rev`, its paint coverage, and the shaft stop distance behind the tip.
    fn leaf(&self, shape: HeadShape, context: MarkContext) -> Result<(MarkDrawable, f64), String> {
        let closed = shape.closed();
        let fill = self.fill.unwrap_or(closed);
        // A triangle or stealth named without paint flags is a plain solid;
        // naming either flag selects the general defaults, which also stroke.
        let solid = matches!(shape, HeadShape::Triangle | HeadShape::Stealth)
            && self.fill.is_none()
            && self.stroke.is_none();
        let stroke = self.stroke.unwrap_or(!solid);
        let thickness = if stroke { context.line_thickness } else { 0.0 };
        let (length, width) = self.size(shape, context)?;
        let mut head = self.head(shape, length, width, thickness / 2.0)?;
        if self.rev.unwrap_or(false) {
            head.turn_around();
        }
        if !head.path.is_finite() || !head.rear.is_finite() {
            return Err("mark dimensions produced nonfinite geometry".into());
        }
        let (join, cap) = shape.joint();
        let miter_limit = head.miter * MITER_HEADROOM;
        let mut outline = BezPath::new();
        if fill {
            // Orient the filled region positively: a mirrored head must not
            // cancel stroke coverage where the two paints overlap.
            let region = if head.path.area() < 0.0 {
                head.path.reverse_subpaths()
            } else {
                head.path.clone()
            };
            outline.extend(region.elements().iter().copied());
        }
        if thickness > 0.0 {
            let pen = Stroke::new(thickness)
                .with_join(join)
                .with_caps(cap)
                .with_miter_limit(miter_limit);
            let painted = kurbo::stroke(
                head.path.elements().iter().copied(),
                &pen,
                &StrokeOpts::default(),
                1e-6,
            );
            outline.extend(painted.elements().iter().copied());
        }
        let outline = CurvePathOutput { path: outline };
        let drawable = MarkDrawable {
            path: CurvePathOutput { path: head.path },
            outline_bounds: outline.control_hull_bounds(),
            outline,
            fill,
            stroke,
            join: if join == Join::Round {
                "round"
            } else {
                "miter"
            }
            .into(),
            cap: if cap == Cap::Round { "round" } else { "butt" }.into(),
            miter_limit,
        };
        // Every head lies behind its tip, so the shaft stops |rear| behind it.
        debug_assert!(head.rear <= 0.0, "{} reaches past its tip", self.shape);
        Ok((drawable, head.rear.abs()))
    }

    /// Painted length and width in drawing units, applying shape defaults.
    fn size(&self, shape: HeadShape, context: MarkContext) -> Result<(f64, f64), String> {
        // Pointed heads grow with the line and stay legible on hairlines;
        // cross strokes are a little shorter. Centred shapes scale with the line.
        let pointed = MarkSize {
            points: 3.0,
            ratio: 4.5,
        };
        let cross = MarkSize {
            points: 2.4,
            ratio: 3.6,
        };
        let lines = |ratio| MarkSize { points: 0.0, ratio };
        let given = |size: Option<MarkSize>| size.map(|size| size.resolve(context)).transpose();
        let (length, width) = match shape {
            Triangle | Stealth | Round | Straight => {
                let length = given(self.length)?.unwrap_or(pointed.resolve(context)?);
                (length, given(self.width)?.unwrap_or(0.8 * length))
            }
            Tikz | Barb | Hooks => (0.0, given(self.width)?.unwrap_or(pointed.resolve(context)?)),
            Bar => (0.0, given(self.width)?.unwrap_or(cross.resolve(context)?)),
            Bracket => {
                let width = given(self.width)?.unwrap_or(cross.resolve(context)?);
                (given(self.length)?.unwrap_or(0.3 * width), width)
            }
            Circle | Square | Diamond => {
                // The default diamond is the default square turned by 45°.
                let side = if shape == Diamond { 4.0 * SQRT_2 } else { 4.0 };
                let length = given(self.length)?.unwrap_or(lines(side).resolve(context)?);
                (length, given(self.width)?.unwrap_or(length))
            }
            Rays => (
                given(self.length)?.unwrap_or(lines(2.8).resolve(context)?),
                0.0,
            ),
        };
        if length < 0.0 || width < 0.0 {
            return Err("mark dimensions must be nonnegative".into());
        }
        Ok((length, width))
    }

    /// Drawn geometry in the local frame: the tip is the origin and the line
    /// arrives from negative x. The paint of a head fills its nominal
    /// `length` x `width` box: strokes of half width `half_stroke` run inside
    /// the box, so a stroked head is as large as an unstroked one. A head too
    /// small to hold its stroke that way grows just enough: a polygon about
    /// the origin, keeping its shape, and any other shape in each dimension
    /// that is too small.
    fn head(
        &self,
        shape: HeadShape,
        length: f64,
        width: f64,
        half_stroke: f64,
    ) -> Result<LocalHead, String> {
        let (l, half, s) = (length, width / 2.0, half_stroke);
        // Rounded and straight parts sit half a stroke inside the box; a
        // dimension too small for that grows to keep them `STROKE_ROOM`
        // clear of collapsing onto the line or a point. `long` and `wide` are
        // half the box's length and width after that growth.
        let least = STROKE_ROOM * s;
        let (long, wide) = ((l / 2.0).max(least), half.max(least));
        // Centred shapes straddle the tip. `align: end` puts the front of a
        // bar, circle, square or diamond on the tip, and the hub of the rays
        // one ray length behind it.
        let aligned = self.align == Some(MarkAlign::End);
        let mut path = BezPath::new();
        Ok(match shape {
            Triangle | Stealth => {
                let inset = self
                    .inset
                    .unwrap_or(if shape == Triangle { 0.0 } else { 0.4 });
                return Polygon::dart(Point::ORIGIN, l, width, inset).stroked_head(s, l);
            }
            Square | Diamond => {
                let scale = l.max(width);
                let scale = if scale == 0.0 { 1.0 } else { scale };
                let long = 0.5 * (l / scale);
                let half = 0.5 * (width / scale);
                let center = if aligned { -long } else { 0.0 };
                let (front, back) = (center + long, center - long);
                let corners = if shape == Square {
                    vec![
                        Point::new(back, -half),
                        Point::new(front, -half),
                        Point::new(front, half),
                        Point::new(back, half),
                    ]
                } else {
                    vec![
                        Point::new(front, 0.0),
                        Point::new(center, half),
                        Point::new(back, 0.0),
                        Point::new(center, -half),
                    ]
                };
                return Polygon {
                    corners,
                    scale,
                    positive_area: l > 0.0 && width > 0.0,
                }
                .stroked_head(s, 0.0);
            }
            Round => {
                // Round joins add half a stroke around every corner, so the
                // core dart sits that far inside each side of the box and the
                // rounded paint still touches the tip and the box.
                let core = 2.0 * (long - s);
                let inset = self.inset.unwrap_or(0.4);
                Polygon::dart(Point::new(-s, 0.0), core, 2.0 * (wide - s), inset)
                    .into_head(2.0 * long)
            }
            Straight => LocalHead::chevron(2.0 * long, wide, s),
            Tikz => {
                // Two quarter circles leave the tip along the line in a cusp
                // and curl outward to the full width.
                let radius = wide - s;
                let radii = (radius, radius);
                let upper = Arc::new((-s, radius), radii, PI, FRAC_PI_2, 0.0);
                let lower = Arc::new((-s, -radius), radii, FRAC_PI_2, FRAC_PI_2, 0.0);
                let tolerance = LocalHead::arc_tolerance(radius);
                path.move_to((-s - radius, radius));
                path.extend(upper.append_iter(tolerance));
                path.extend(lower.append_iter(tolerance));
                LocalHead::open(path, -s, 0.0)
            }
            // The arcs of barbs and hooks keep a radius of at least
            // `STROKE_ROOM` stroke radii, so that their strokes do not wrap
            // past the arcs' centres and out of the box.
            Barb => {
                // An arc of the full width bulges toward the tip, swept
                // symmetrically about the line; its box is half as long as wide.
                let radius = (wide - s).max(least);
                let center = -radius - s;
                let sweep = self.arc.unwrap_or(PI);
                let arc = Arc::new((center, 0.0), (radius, radius), -sweep / 2.0, sweep, 0.0);
                path.extend(arc.path_elements(LocalHead::arc_tolerance(radius)));
                LocalHead::arc_contacts(path, center, radius, -sweep / 2.0, sweep, -center)
            }
            Hooks => {
                // Two hooks, each half the width, leave one root on the line
                // toward the tip and curl away from the line by `arc`.
                let radius = ((wide - s) / 2.0).max(least);
                let root = -radius - s;
                let sweep = self.arc.unwrap_or(PI);
                let radii = (radius, radius);
                let tolerance = LocalHead::arc_tolerance(radius);
                let upper = Arc::new((root, radius), radii, -FRAC_PI_2, sweep, 0.0);
                let lower = Arc::new((root, -radius), radii, FRAC_PI_2, -sweep, 0.0);
                path.extend(upper.path_elements(tolerance));
                path.extend(lower.path_elements(tolerance));
                LocalHead::open(path, root, -root)
            }
            Bar => {
                let x = if aligned { -s } else { 0.0 };
                path.move_to((x, -half));
                path.line_to((x, half));
                LocalHead::open(path, x, 0.0)
            }
            Bracket => {
                // A spine across the line at the front, with arms running back
                // from its two right-angled corners.
                let (end, arm) = (2.0 * long, wide - s);
                path.move_to((-end, -arm));
                path.line_to((-s, -arm));
                path.line_to((-s, arm));
                path.line_to((-end, arm));
                LocalHead {
                    miter: SQRT_2,
                    ..LocalHead::open(path, -s, end)
                }
            }
            Circle => {
                let center = if aligned { -long } else { 0.0 };
                let radii = Vec2::new(long - s, wide - s);
                let ellipse = Ellipse::new((center, 0.0), radii, 0.0);
                LocalHead {
                    path: ellipse.to_path(LocalHead::arc_tolerance(radii.x.max(radii.y))),
                    rear: center - radii.x,
                    front: center + radii.x,
                    span: 0.0,
                    miter: 1.0,
                }
            }
            Rays => {
                // Evenly spaced spokes from a hub on the line, none along the
                // line, where the shaft would hide it. An even count starts
                // at π/n and is symmetric about the line; an odd count cannot
                // be without a spoke on it, so it starts perpendicular to the
                // line and is symmetric about that perpendicular.
                let count = self.n.unwrap_or(4);
                let phase = self.phase.unwrap_or(if count.is_multiple_of(2) {
                    PI / count as f64
                } else {
                    -FRAC_PI_2
                });
                let hub = Point::new(if aligned { -l } else { 0.0 }, 0.0);
                for index in 0..count {
                    let angle = phase + TAU * index as f64 / count as f64;
                    path.move_to(hub);
                    path.line_to(hub + Vec2::from_angle(angle) * l);
                }
                LocalHead::open(path, hub.x, 0.0)
            }
        })
    }
}

/// A head grown to hold its stroke keeps this much room, so its drawn outline
/// stays clear of collapsing onto a point or the line.
const STROKE_ROOM: f64 = 1.01;

/// Renderers bevel a corner whose miter ratio exceeds the stroke's miter
/// limit. Each head's limit sits this factor above its sharpest designed
/// corner, leaving room for single-precision renderers, so every corner keeps
/// its point however slender the head.
const MITER_HEADROOM: f64 = 1.001;

/// A mark is one head or a composite of marks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MarkShape {
    Head(HeadShape),
    Combine,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HeadShape {
    Triangle,
    Stealth,
    Round,
    Straight,
    Tikz,
    Barb,
    Hooks,
    Bar,
    Bracket,
    Circle,
    Square,
    Diamond,
    Rays,
}

impl MarkShape {
    fn parse(name: &str) -> Result<Self, String> {
        Ok(Self::Head(match name {
            "triangle" => Triangle,
            "stealth" => Stealth,
            "round" => Round,
            "straight" => Straight,
            "tikz" => Tikz,
            "barb" => Barb,
            "hooks" => Hooks,
            "bar" => Bar,
            "bracket" => Bracket,
            "circle" => Circle,
            "square" => Square,
            "diamond" => Diamond,
            "rays" => Rays,
            "combine" => return Ok(Self::Combine),
            _ => return Err(format!("unknown mark shape: {name}")),
        }))
    }

    /// Options read by the shape; every shape also takes fit and shorten.
    fn options(self) -> &'static [&'static str] {
        let Self::Head(shape) = self else {
            return &["parts"];
        };
        match shape {
            Triangle | Stealth | Round => &["length", "width", "inset", "rev", "fill", "stroke"],
            Straight | Bracket => &["length", "width", "rev", "stroke"],
            Tikz => &["width", "stroke"],
            Barb | Hooks => &["width", "arc", "rev", "stroke"],
            Bar => &["width", "align", "stroke"],
            Circle | Square | Diamond => &["length", "width", "align", "fill", "stroke"],
            Rays => &["length", "n", "phase", "align", "stroke"],
        }
    }
}

impl HeadShape {
    /// Closed outlines enclose a region, which is filled by default.
    fn closed(self) -> bool {
        matches!(self, Triangle | Stealth | Round | Circle | Square | Diamond)
    }

    /// Stroke join and cap. Shapes drawn as rounded or smooth lines get round
    /// joins or caps; the others keep mitred corners and butt ends.
    fn joint(self) -> (Join, Cap) {
        match self {
            Self::Round => (Join::Round, Cap::Butt),
            Self::Tikz => (Join::Round, Cap::Round),
            Self::Straight => (Join::Miter, Cap::Round),
            _ => (Join::Miter, Cap::Butt),
        }
    }
}

/// A leaf head's drawn path in the local frame, with the first and last
/// points (`rear` <= `front`) where it meets the line, the longitudinal span
/// that `rev` turns it around in, and the largest miter ratio among its
/// corners (1 when it has none), which round joins ignore.
#[derive(Debug)]
struct LocalHead {
    path: BezPath,
    rear: f64,
    front: f64,
    span: f64,
    miter: f64,
}

impl LocalHead {
    /// Open paths built here meet the line in a single point.
    fn open(path: BezPath, axis: f64, span: f64) -> Self {
        Self {
            path,
            rear: axis,
            front: axis,
            span,
            miter: 1.0,
        }
    }

    /// Axis crossings of a circular arc, including its endpoint crossings.
    /// Enumerate multiples of π in the accepted angular interval instead of
    /// assuming that every sweep meets the line only at the front. This also
    /// retains the open gap immediately before a complete positive/negative turn.
    fn arc_contacts(
        path: BezPath,
        center: f64,
        radius: f64,
        start: f64,
        sweep: f64,
        span: f64,
    ) -> Self {
        let finish = start + sweep;
        let first = (start.min(finish) / PI).ceil() as i32;
        let last = (start.max(finish) / PI).floor() as i32;
        let (mut rear, mut front) = (f64::INFINITY, f64::NEG_INFINITY);
        for crossing in first..=last {
            let x = center + if crossing % 2 == 0 { radius } else { -radius };
            rear = rear.min(x);
            front = front.max(x);
        }
        Self {
            path,
            rear,
            front,
            span,
            miter: 1.0,
        }
    }

    /// An open chevron stroked `half_stroke` to each side whose arms end in
    /// round caps touching the back corners of its `length` x 2 `half` box.
    /// Each arm's outer edge is the tangent from the tip to its cap, so the
    /// mitred point lands on the tip: an arm leaves the line at the angle of
    /// its end seen from the tip, opened by the tangent's offset, and the arms'
    /// centrelines meet `half_stroke / sin(angle)` behind the tip. The caller
    /// keeps the box large enough for the caps (`length > 2 half_stroke`,
    /// `half > half_stroke`).
    fn chevron(length: f64, half: f64, half_stroke: f64) -> Self {
        let s = half_stroke;
        let arm = Point::new(s - length, half - s);
        let apex = if s == 0.0 {
            Point::ORIGIN
        } else {
            let angle = arm.y.atan2(-arm.x) + (s / arm.to_vec2().hypot()).asin();
            Point::new(-s / angle.sin(), 0.0)
        };
        let mut path = BezPath::new();
        path.move_to(arm);
        path.line_to(apex);
        path.line_to((arm.x, -arm.y));
        Self {
            miter: if arm.y > 0.0 {
                (arm - apex).hypot() / arm.y
            } else {
                1.0
            },
            ..Self::open(path, apex.x, length)
        }
    }

    /// Flatten arcs to within 0.1% of their radius: one cubic per quarter turn.
    fn arc_tolerance(radius: f64) -> f64 {
        (radius * 1e-3).max(1e-9)
    }

    /// Point the head the other way within its span, `x -> -span - x`.
    fn turn_around(&mut self) {
        self.path
            .apply_affine(Affine::new([-1.0, 0.0, 0.0, 1.0, -self.span, 0.0]));
        (self.rear, self.front) = (-self.span - self.front, -self.span - self.rear);
    }
}

/// Closed polygon with counterclockwise vertices.
#[derive(Debug, Clone)]
struct Polygon {
    // Coordinates stay normalized until growth and inset have been resolved.
    // In particular, a tiny positive dart must not lose its notch to underflow.
    corners: Vec<Point>,
    scale: f64,
    // Parameter-derived, without multiplying tiny dimensions or relying on
    // coordinates that may have underflowed during normalization.
    positive_area: bool,
}

impl Polygon {
    /// Tip, upper wing, notch and lower wing. Form the notch using its fraction
    /// in normalized coordinates, rather than subtracting two physical lengths.
    fn dart(tip: Point, length: f64, width: f64, inset: f64) -> Self {
        let positive_area = length > 0.0 && width > 0.0 && inset < 1.0;
        let scale = length.max(width).max(tip.x.abs()).max(tip.y.abs());
        let scale = if scale == 0.0 { 1.0 } else { scale };
        let tip = Point::new(tip.x / scale, tip.y / scale);
        let length = length / scale;
        let half_width = 0.5 * (width / scale);
        Self {
            corners: vec![
                tip,
                tip + Vec2::new(-length, half_width),
                tip + Vec2::new(length * (inset - 1.0), 0.0),
                tip + Vec2::new(-length, -half_width),
            ],
            scale,
            positive_area,
        }
    }

    /// Grow before offsetting a positive-area polygon. Every moved edge is
    /// linear in the inset, so its first collapse gives an analytic limit.
    /// This has no absolute search resolution and never uses failed numerical
    /// offsetting as permission to stroke the nominal positive-area outline.
    fn stroked_head(mut self, half_stroke: f64, mut span: f64) -> Result<LocalHead, String> {
        if half_stroke == 0.0 || self.degenerate() {
            return Ok(self.into_head(span));
        }
        let count = self.corners.len();
        let slides = (0..count)
            .map(|i| self.corner_slide(i))
            .collect::<Option<Vec<_>>>()
            .ok_or("positive-area mark has an unrepresentable offset corner")?;
        let mut limit = f64::INFINITY;
        for i in 0..count {
            let next = (i + 1) % count;
            let edge = self.corners[next] - self.corners[i];
            let length = edge.hypot();
            let direction = Vec2::new(edge.x / length, edge.y / length);
            let contraction = slides[i].dot(direction) - slides[next].dot(direction);
            if contraction > 0.0 {
                limit = limit.min(length / contraction);
            }
        }
        if !limit.is_finite() || limit <= 0.0 {
            return Err("positive-area mark has an unrepresentable inset limit".into());
        }
        let scale = self.scale.max(STROKE_ROOM * half_stroke / limit);
        span = (span / self.scale) * scale;
        if !scale.is_finite() || !span.is_finite() {
            return Err("grown mark dimensions exceed finite drawing coordinates".into());
        }
        let distance = half_stroke / scale;
        for (corner, slide) in self.corners.iter_mut().zip(&slides) {
            *corner += *slide * distance;
        }
        self.scale = scale;
        let miter = slides.into_iter().map(Vec2::hypot).fold(1.0, f64::max);
        let mut head = self.into_head(span);
        head.miter = miter;
        Ok(head)
    }

    /// Scale each vector separately before orientation tests, so neither tiny
    /// dimensions nor very different edge lengths square into underflow.
    fn scaled(vector: Vec2) -> Option<Vec2> {
        let largest = vector.x.abs().max(vector.y.abs());
        if largest == 0.0 {
            return None;
        }
        // An arbitrary divisor can round two nearly opposite edges into a
        // different turn. Binary scaling preserves every representable bit.
        let bits = largest.to_bits();
        let exponent = (bits >> 52) & 0x7ff;
        let scale = if exponent == 0 {
            f64::from_bits(1_u64 << (63 - bits.leading_zeros()))
        } else {
            f64::from_bits(exponent << 52)
        };
        Some(Vec2::new(vector.x / scale, vector.y / scale))
    }

    /// Compensated 2x2 determinant. Near a sharp corner, the two products
    /// almost cancel; retain the product rounding residual instead of deriving
    /// the small turn from 1 + cos(turn).
    fn determinant(a: Vec2, b: Vec2) -> f64 {
        let product = a.y * b.x;
        a.x.mul_add(b.y, -product) - a.y.mul_add(b.x, -product)
    }

    fn degenerate(&self) -> bool {
        !self.positive_area
    }

    /// Offset-line intersection, using the stable half-angle identity on
    /// either side of a right angle. None means an exact doubled-back/zero
    /// edge or an unrepresentable result; the caller distinguishes true zero
    /// area before permitting unstroked-as-given geometry.
    fn corner_slide(&self, index: usize) -> Option<Vec2> {
        let count = self.corners.len();
        let here = self.corners[index];
        let before = Self::scaled(here - self.corners[(index + count - 1) % count])?;
        let after = Self::scaled(self.corners[(index + 1) % count] - here)?;
        let norm = before.hypot() * after.hypot();
        let sine = Self::determinant(before, after) / norm;
        let cosine = before.dot(after) / norm;
        let tangent = if cosine < 0.0 {
            if sine == 0.0 {
                return None;
            }
            (cosine - 1.0) / sine
        } else {
            -sine / (1.0 + cosine)
        };
        let direction = before / before.hypot();
        let slide = Vec2::new(-direction.y, direction.x) + direction * tangent;
        slide.is_finite().then_some(slide)
    }

    fn miter(&self) -> f64 {
        (0..self.corners.len())
            .filter_map(|i| self.corner_slide(i))
            .map(Vec2::hypot)
            .fold(1.0, f64::max)
    }

    /// The closed outline as path elements.
    fn elements(&self) -> impl Iterator<Item = PathEl> + '_ {
        let corners = self.corners.iter().enumerate().map(|(index, &corner)| {
            if index == 0 {
                PathEl::MoveTo((corner.to_vec2() * self.scale).to_point())
            } else {
                PathEl::LineTo((corner.to_vec2() * self.scale).to_point())
            }
        });
        corners.chain(iter::once(PathEl::ClosePath))
    }

    /// The polygon as a drawn path, meeting the line between its outermost
    /// crossings of `y = 0`, turned around by `rev` within `span`. The
    /// polygons built here are symmetric about the line, so they cross it.
    fn into_head(self, span: f64) -> LocalHead {
        let count = self.corners.len();
        let (mut rear, mut front) = (f64::INFINITY, f64::NEG_INFINITY);
        for index in 0..count {
            let (a, b) = (self.corners[index], self.corners[(index + 1) % count]);
            let crossing = if a.y == 0.0 {
                Some(a.x)
            } else if b.y != 0.0 && (a.y < 0.0) != (b.y < 0.0) {
                Some(a.x + (b.x - a.x) * a.y / (a.y - b.y))
            } else {
                None
            };
            if let Some(x) = crossing {
                rear = rear.min(x);
                front = front.max(x);
            }
        }
        LocalHead {
            path: self.elements().collect(),
            rear: rear * self.scale,
            front: front * self.scale,
            span,
            miter: self.miter(),
        }
    }
}

impl PreparedMark {
    fn longitudinal_center(paths: &[MarkDrawable]) -> f64 {
        paths
            .iter()
            .filter(|drawable| !drawable.outline.path.is_empty())
            .map(|drawable| drawable.outline.path.bounding_box())
            .reduce(|a, b| a.union(b))
            .map_or(0.0, |bounds| (bounds.x0 + bounds.x1) / 2.0)
    }
}

impl MarkGeometrySpec {
    pub fn geometry(self) -> Result<MarkGeometryBatchOutput, String> {
        let templates = self
            .templates
            .iter()
            .map(|spec| spec.mark.prepare(spec.context))
            .collect::<Result<Vec<_>, _>>()?;
        let carriers = self
            .carriers
            .into_iter()
            .enumerate()
            .map(|(index, carrier)| {
                PreparedCarrier::new(carrier.path)
                    .map_err(|error| format!("mark carrier {index}: {error}"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut results = Vec::with_capacity(self.placements.len());
        let mut groups = std::collections::BTreeMap::<usize, SelectedShaft>::new();
        for placement in self.placements {
            let template = templates
                .get(placement.template)
                .ok_or("invalid template index")?;
            let carrier = carriers
                .get(placement.carrier)
                .ok_or("invalid carrier index")?;
            let mut fitted = template.fit_prepared(carrier, &placement, self.mode)?;
            if self.mode == MarkGeometryMode::Selected {
                let group = groups
                    .entry(placement.carrier)
                    .or_insert_with(|| SelectedShaft {
                        template: placement.template,
                        window: ShaftWindow::default(),
                        start_selected: false,
                        end_selected: false,
                    });
                let first = &templates[group.template];
                if first.line_thickness != template.line_thickness
                    || first.shaft_style != template.shaft_style
                {
                    return Err(format!(
                        "selected marks on carrier {} have conflicting shaft thickness or styles",
                        placement.carrier
                    ));
                }
                let endpoint = match placement.station {
                    MarkStation::Start => Some(&mut group.start_selected),
                    MarkStation::End => Some(&mut group.end_selected),
                    _ => None,
                };
                if let Some(selected) = endpoint {
                    if *selected {
                        return Err(format!(
                            "multiple selected choices for one endpoint of carrier {}",
                            placement.carrier
                        ));
                    }
                    *selected = true;
                }
                group.window.apply(&fitted.shaft_edit);
                fitted.geometry.shaft.path = BezPath::new();
            }
            results.push(fitted.geometry);
        }
        let shafts = groups
            .into_iter()
            .map(|(carrier, group)| {
                let template = &templates[group.template];
                let shaft = carriers[carrier].shaft(&group.window)?;
                let outline = template.shaft_outline(&shaft)?;
                Ok(MarkShaftOutput {
                    carrier,
                    shaft: CurvePathOutput { path: shaft },
                    footprint: CurvePathOutput {
                        path: outline.clone(),
                    },
                    shaft_outline: CurvePathOutput { path: outline },
                    shaft_style: template.shaft_style,
                })
            })
            .collect::<Result<Vec<_>, String>>()?;
        Ok(MarkGeometryBatchOutput {
            marks: results,
            shafts,
        })
    }
}

impl PreparedMark {
    fn shaft_outline(&self, shaft: &BezPath) -> Result<BezPath, String> {
        if !shaft.is_finite() {
            return Err("shaft geometry must be finite".into());
        }
        let outline = if self.line_thickness > 0.0 {
            kurbo::stroke(
                shaft.elements().iter().copied(),
                &self.shaft_style.stroke(self.line_thickness),
                &StrokeOpts::default(),
                1e-6,
            )
        } else {
            BezPath::new()
        };
        if !outline.is_finite() {
            return Err("painted shaft geometry must be finite".into());
        }
        Ok(outline)
    }

    pub fn paths(&self) -> &[MarkDrawable] {
        &self.paths
    }

    pub fn end(&self) -> f64 {
        self.end
    }

    /// Scalar placement shares the batched fitting implementation. The batch
    /// additionally reuses one arc table for every placement on each carrier.
    pub fn place(
        &self,
        carrier: &CurvePathOutput,
        placement: &MarkPlacement,
    ) -> Result<MarkGeometryOutput, String> {
        self.place_prepared(&PreparedCarrier::new(carrier.path.clone())?, placement)
    }

    fn place_prepared(
        &self,
        carrier: &PreparedCarrier,
        placement: &MarkPlacement,
    ) -> Result<MarkGeometryOutput, String> {
        Ok(self
            .fit_prepared(carrier, placement, MarkGeometryMode::Candidates)?
            .geometry)
    }

    fn fit_prepared(
        &self,
        carrier: &PreparedCarrier,
        placement: &MarkPlacement,
        mode: MarkGeometryMode,
    ) -> Result<FittedMark, String> {
        let template = self;
        let length = carrier.length;
        let station = match placement.station {
            MarkStation::Start => 0.0,
            MarkStation::End => length,
            MarkStation::Ratio { value } => value * length,
            MarkStation::Distance { value } => value,
        };
        if !station.is_finite() || !placement.shift.is_finite() {
            return Err("placement station and shift must be finite".into());
        }
        let sign = if placement.direction == MarkDirection::Backward {
            -1.0
        } else {
            1.0
        };
        let endpoint = matches!(placement.station, MarkStation::Start | MarkStation::End);
        // Convert a painted-center station to the local origin before inward
        // clamping and fitting. User shifts remain signed carrier distances,
        // independent of the mark's direction.
        let center_offset = if endpoint {
            0.0
        } else {
            -sign * template.longitudinal_center
        };
        let station = (station + center_offset + placement.shift).clamp(0.0, length);
        let origin = carrier.point(station)?;
        let (contact, back_at) = carrier.chord_contact(station, template.end, sign)?;
        let tangent = carrier
            .trimmer
            .frame(station)?
            .map(|frame| Vec2::new(frame.tangent.x, frame.tangent.y) * sign)
            .unwrap_or(Vec2::new(sign, 0.0));
        let bend = endpoint && template.fit == MarkFit::Bend;
        if bend
            && carrier
                .path
                .elements()
                .iter()
                .any(|el| matches!(el, PathEl::ClosePath))
        {
            return Err("endpoint bend on a closed carrier is unsupported; use chord fit or an interior station".into());
        }
        let vector = if bend {
            tangent
        } else {
            (origin - contact) * if template.end < 0.0 { -1.0 } else { 1.0 }
        };
        let direction = if vector.hypot() > 1e-12 {
            vector / vector.hypot()
        } else if tangent.hypot() > 1e-12 {
            tangent / tangent.hypot()
        } else {
            Vec2::new(sign, 0.0)
        };
        let transform = Affine::new([
            direction.x,
            direction.y,
            -direction.y,
            direction.x,
            origin.x,
            origin.y,
        ]);
        let mut paths = template.paths.clone();
        for path in &mut paths {
            path.path.path.apply_affine(transform);
            path.outline.path.apply_affine(transform);
            path.outline_bounds = path.outline.control_hull_bounds();
        }
        let retraction = template.end * template.shorten;
        if !retraction.is_finite() {
            return Err("mark retraction must be finite".into());
        }
        // The full-size mark is never silently scaled to fit. On too-short
        // carriers its geometric back can be off-carrier; shaft_contact
        // separately describes the inward-clamped painted shaft endpoint.
        let back = origin - direction * template.end;
        let shaft_edit = if !endpoint || retraction == 0.0 {
            ShaftEdit::None
        } else if bend {
            // Clamp endpoint retraction to the carrier's available inward
            // arc length. Controls, including a short final segment's
            // controls, remain exactly as supplied.
            let inward = match placement.station {
                MarkStation::Start => -1.0,
                _ => 1.0,
            };
            let available = if matches!(placement.station, MarkStation::Start) {
                length - station
            } else {
                station
            };
            let moved = origin
                - direction * (retraction * sign * inward).clamp(0.0, available) * sign * inward;
            if matches!(placement.station, MarkStation::Start) {
                ShaftEdit::BendStart(moved)
            } else {
                ShaftEdit::BendEnd(moved)
            }
        } else {
            let painted_at = (station + (back_at - station) * template.shorten).clamp(0.0, length);
            if matches!(placement.station, MarkStation::Start) {
                ShaftEdit::TrimStart(painted_at)
            } else {
                ShaftEdit::TrimEnd(length - painted_at)
            }
        };
        let mut window = ShaftWindow::default();
        window.apply(&shaft_edit);
        let shaft = carrier.shaft(&window)?;
        let mut footprint = BezPath::new();
        for drawable in &paths {
            footprint.extend(drawable.outline.path.elements().iter().copied());
        }
        let shaft_outline = if mode == MarkGeometryMode::Candidates {
            self.shaft_outline(&shaft)?
        } else {
            BezPath::new()
        };
        footprint.extend(shaft_outline.elements().iter().copied());
        if !shaft.is_finite()
            || !footprint.is_finite()
            || paths
                .iter()
                .any(|path| !path.path.path.is_finite() || !path.outline.path.is_finite())
        {
            return Err("mark placement produced nonfinite geometry".into());
        }
        let shaft_contact = if !endpoint {
            contact
        } else if matches!(placement.station, MarkStation::Start) {
            shaft
                .segments()
                .next()
                .map(|segment| segment.start())
                .unwrap_or(carrier.point(window.start_outset)?)
        } else {
            shaft
                .segments()
                .last()
                .map(|segment| segment.end())
                .unwrap_or(carrier.point(length - window.end_outset)?)
        };
        let footprint = CurvePathOutput { path: footprint };
        Ok(FittedMark {
            shaft_edit,
            geometry: MarkGeometryOutput {
                template: placement.template,
                carrier: placement.carrier,
                end: template.end,
                paths,
                tip: origin.into(),
                back: back.into(),
                shaft_contact: shaft_contact.into(),
                shaft: CurvePathOutput { path: shaft },
                shaft_outline: CurvePathOutput {
                    path: shaft_outline,
                },
                shaft_style: self.shaft_style,
                footprint_bounds: footprint.control_hull_bounds(),
                footprint,
            },
        })
    }
}

/// Borrowed transport view: only path-valued fields become native-path CBOR
/// byte strings. Metadata keeps the native DTO's names, order, and values.
struct MarkPacket<'a, T>(&'a T);

struct NativePathPacket<'a>(&'a CurvePathOutput);

impl Serialize for NativePathPacket<'_> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut bytes = Vec::new();
        ciborium::ser::into_writer(self.0, &mut bytes).map_err(serde::ser::Error::custom)?;
        serializer.serialize_bytes(&bytes)
    }
}

// Keep the owning native path serializer authoritative, without constructing
// and rewriting a second coordinate tree for the batch.
macro_rules! packet_fields {
    ($serializer:ident, $($key:literal => $value:expr),+ $(,)?) => {{
        use serde::ser::SerializeStruct;
        let mut output = $serializer.serialize_struct("MarkPacket", [$($key),+].len())?;
        $(output.serialize_field($key, &$value)?;)+
        output.end()
    }};
}

impl Serialize for MarkPacket<'_, MarkGeometryBatchOutput> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let marks: Vec<_> = self.0.marks.iter().map(MarkPacket).collect();
        let shafts: Vec<_> = self.0.shafts.iter().map(MarkPacket).collect();
        packet_fields!(serializer, "marks" => marks, "shafts" => shafts)
    }
}

impl Serialize for MarkPacket<'_, MarkGeometryOutput> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mark = self.0;
        let paths: Vec<_> = mark.paths.iter().map(MarkPacket).collect();
        packet_fields!(serializer,
            "template" => mark.template, "carrier" => mark.carrier, "end" => mark.end,
            "paths" => paths, "tip" => mark.tip, "back" => mark.back,
            "shaft-contact" => mark.shaft_contact,
            "shaft" => NativePathPacket(&mark.shaft),
            "shaft-outline" => NativePathPacket(&mark.shaft_outline),
            "shaft-style" => mark.shaft_style,
            "footprint" => NativePathPacket(&mark.footprint),
            "footprint-bounds" => mark.footprint_bounds,
        )
    }
}

impl Serialize for MarkPacket<'_, MarkShaftOutput> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let shaft = self.0;
        packet_fields!(serializer,
            "carrier" => shaft.carrier, "shaft" => NativePathPacket(&shaft.shaft),
            "shaft-outline" => NativePathPacket(&shaft.shaft_outline),
            "shaft-style" => shaft.shaft_style,
            "footprint" => NativePathPacket(&shaft.footprint),
        )
    }
}

impl Serialize for MarkPacket<'_, MarkDrawable> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let path = self.0;
        packet_fields!(serializer,
            "path" => NativePathPacket(&path.path), "outline" => NativePathPacket(&path.outline),
            "outline-bounds" => path.outline_bounds,
            "fill" => path.fill, "stroke" => path.stroke,
            "join" => path.join, "cap" => path.cap, "miter-limit" => path.miter_limit,
        )
    }
}

/// Packed companion to `mark_geometry_bytes`; runs the same fitter once.
pub fn mark_geometry_packed_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: MarkGeometrySpec =
        ciborium::de::from_reader(arg).map_err(|e| format!("invalid mark geometry CBOR: {e}"))?;
    let output = spec.geometry()?;
    let mut bytes = Vec::new();
    ciborium::ser::into_writer(&MarkPacket(&output), &mut bytes)
        .map_err(|e| format!("mark geometry serialization failed: {e}"))?;
    Ok(bytes)
}

pub fn mark_geometry_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: MarkGeometrySpec =
        ciborium::de::from_reader(arg).map_err(|e| format!("invalid mark geometry CBOR: {e}"))?;
    let output = spec.geometry()?;
    let mut bytes = Vec::new();
    ciborium::ser::into_writer(&output, &mut bytes)
        .map_err(|e| format!("mark geometry serialization failed: {e}"))?;
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ciborium::value::Value;
    use kurbo::PathSeg;

    const SHAPES: [&str; 13] = [
        "triangle", "straight", "stealth", "round", "tikz", "barb", "hooks", "bar", "bracket",
        "circle", "square", "diamond", "rays",
    ];
    const CONTEXT: MarkContext = MarkContext {
        units_per_pt: 2.0,
        line_thickness: 2.0,
        shaft_stroke: MarkStrokeStyle {
            cap: Cap::Butt,
            join: Join::Miter,
            miter_limit: 4.0,
        },
    };

    fn mark(shape: &str) -> MarkSpec {
        if shape == "combine" {
            let mut spec = mark("bar");
            spec.shape = shape.into();
            spec.parts = Some(vec![MarkPart::Mark(Box::new(mark("bar")))]);
            return spec;
        }
        let mut bytes = Vec::new();
        ciborium::ser::into_writer(
            &Value::Map(vec![(
                Value::Text("shape".into()),
                Value::Text(shape.into()),
            )]),
            &mut bytes,
        )
        .unwrap();
        ciborium::de::from_reader(bytes.as_slice()).unwrap()
    }

    fn line(start: Point, end: Point) -> BezPath {
        let mut path = BezPath::new();
        path.move_to(start);
        path.line_to(end);
        path
    }

    fn curve() -> BezPath {
        let mut path = BezPath::new();
        path.move_to((0.0, 0.0));
        path.curve_to((0.0, 40.0), (60.0, 40.0), (60.0, 0.0));
        path
    }

    fn batch(spec: MarkSpec, path: BezPath) -> MarkGeometrySpec {
        MarkGeometrySpec {
            mode: MarkGeometryMode::Candidates,
            templates: vec![MarkTemplateSpec {
                mark: spec,
                context: CONTEXT,
            }],
            carriers: vec![CurvePathOutput { path }],
            placements: vec![MarkPlacement {
                template: 0,
                carrier: 0,
                station: MarkStation::End,
                shift: 0.0,
                direction: MarkDirection::Forward,
            }],
        }
    }

    fn geometry(spec: MarkSpec, path: BezPath) -> MarkGeometryOutput {
        batch(spec, path).geometry().unwrap().marks.remove(0)
    }

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < 2e-5, "{a} != {b}");
    }

    fn point_close(a: Point, b: Point) {
        close(a.x, b.x);
        close(a.y, b.y);
    }

    fn wire(value: &impl Serialize) -> Vec<u8> {
        let mut bytes = Vec::new();
        ciborium::ser::into_writer(value, &mut bytes).unwrap();
        bytes
    }

    // Independent oracle for the former Typst traversal, using serialized wire
    // elements rather than the owning method's PathEl representation.
    fn wire_bounds(path: &CurvePathOutput) -> Vec<MarkBounds> {
        let value: Value = ciborium::de::from_reader(wire(path).as_slice()).unwrap();
        let Value::Map(path) = value else { panic!() };
        let Value::Map(path) = &path[0].1 else {
            panic!()
        };
        let Value::Array(elements) = &path[0].1 else {
            panic!()
        };
        let mut regions = vec![Vec::<Point>::new()];
        for element in elements {
            let Value::Map(fields) = element else {
                panic!()
            };
            let field = |name: &str| {
                fields
                    .iter()
                    .find_map(|(key, value)| (key == &Value::Text(name.into())).then_some(value))
            };
            if field("kind") == Some(&Value::Text("move".into())) {
                regions.push(Vec::new());
            }
            for key in ["start", "end", "control", "control-start", "control-end"] {
                if let Some(Value::Array(point)) = field(key) {
                    let coordinate = |value: &Value| match value {
                        Value::Float(value) => *value,
                        Value::Integer(value) => i128::from(*value) as f64,
                        _ => panic!(),
                    };
                    regions
                        .last_mut()
                        .unwrap()
                        .push(Point::new(coordinate(&point[0]), coordinate(&point[1])));
                }
            }
        }
        regions
            .into_iter()
            .filter(|points| !points.is_empty())
            .map(|points| {
                let extremum = |x: bool, minimum: bool| {
                    points
                        .iter()
                        .map(|point| if x { point.x } else { point.y })
                        .reduce(|a, b| {
                            if if minimum { b < a } else { b > a } {
                                b
                            } else {
                                a
                            }
                        })
                        .unwrap()
                };
                MarkBounds {
                    left: extremum(true, true),
                    right: extremum(true, false),
                    bottom: extremum(false, true),
                    top: extremum(false, false),
                }
            })
            .collect()
    }

    #[test]
    fn packed_transport_preserves_native_paths_and_metadata_exactly() {
        fn unpack(value: &mut Value, packets: &mut usize) {
            match value {
                Value::Map(fields) => {
                    for (key, value) in fields {
                        if matches!(
                            key.as_text(),
                            Some("path" | "outline" | "shaft" | "shaft-outline" | "footprint")
                        ) {
                            let Value::Bytes(bytes) = value else {
                                panic!("path-valued field is not an opaque byte string: {key:?}");
                            };
                            let path: CurvePathOutput =
                                ciborium::de::from_reader(bytes.as_slice()).unwrap();
                            assert_eq!(wire(&path), *bytes);
                            *value = ciborium::de::from_reader(bytes.as_slice()).unwrap();
                            *packets += 1;
                        } else {
                            unpack(value, packets);
                        }
                    }
                }
                Value::Array(values) => {
                    for value in values {
                        unpack(value, packets);
                    }
                }
                _ => {}
            }
        }
        let mut composite = mark("combine");
        composite.parts = Some(vec![
            MarkPart::Mark(Box::new(mark("triangle"))),
            MarkPart::Gap(MarkGap {
                gap: MarkSize {
                    points: -2.0,
                    ratio: 0.25,
                },
            }),
            MarkPart::Mark(Box::new(mark("circle"))),
        ]);
        for spec in SHAPES.into_iter().map(mark).chain([composite]) {
            for hidden in [false, true] {
                let mut spec = spec.clone();
                if hidden && spec.shape != "combine" {
                    if matches!(
                        spec.shape.as_str(),
                        "triangle" | "stealth" | "round" | "circle" | "square" | "diamond"
                    ) {
                        spec.fill = Some(false);
                    }
                    spec.stroke = Some(false);
                }
                for carrier in [
                    curve(),
                    line(Point::new(100.0, -20.0), Point::new(-30.0, 40.0)),
                    line(Point::new(-0.0, 0.0), Point::new(0.0, -0.0)),
                    BezPath::new(),
                ] {
                    for mode in [MarkGeometryMode::Candidates, MarkGeometryMode::Selected] {
                        let mut request = batch(spec.clone(), carrier.clone());
                        request.mode = mode;
                        request.placements[0].direction = MarkDirection::Backward;
                        let native = mark_geometry_bytes(&wire(&request)).unwrap();
                        let packed = mark_geometry_packed_bytes(&wire(&request)).unwrap();
                        let mut unpacked: Value =
                            ciborium::de::from_reader(packed.as_slice()).unwrap();
                        let mut packets = 0;
                        unpack(&mut unpacked, &mut packets);
                        assert!(packets >= 3);
                        assert_eq!(wire(&unpacked), native);
                    }
                }
            }
        }
    }

    #[test]
    fn control_hulls_match_wire_traversal_for_catalogue_transforms_and_modes() {
        let mut nested = mark("combine");
        nested.parts = Some(vec![
            MarkPart::Mark(Box::new(mark("circle"))),
            MarkPart::Gap(MarkGap {
                gap: MarkSize {
                    points: -3.0,
                    ratio: 0.25,
                },
            }),
            MarkPart::Mark(Box::new(mark("combine"))),
        ]);
        let mut reversed = mark("triangle");
        reversed.rev = Some(true);
        reversed.stroke = Some(true);
        let specs = SHAPES.into_iter().map(mark).chain([nested, reversed]);
        for spec in specs {
            for carrier in [
                curve(),
                line(Point::new(100.0, -20.0), Point::new(-30.0, 40.0)),
                line(Point::new(-0.0, 0.0), Point::new(0.0, -0.0)),
                BezPath::new(),
            ] {
                for direction in [MarkDirection::Forward, MarkDirection::Backward] {
                    for mode in [MarkGeometryMode::Candidates, MarkGeometryMode::Selected] {
                        let mut request = batch(spec.clone(), carrier.clone());
                        request.mode = mode;
                        request.placements[0].direction = direction;
                        request.placements[0].shift = -3.0;
                        let output = request.geometry().unwrap();
                        for mark in output.marks {
                            assert_eq!(
                                wire(&mark.footprint_bounds),
                                wire(&wire_bounds(&mark.footprint))
                            );
                            for drawable in mark.paths {
                                assert_eq!(
                                    wire(&drawable.outline_bounds),
                                    wire(&wire_bounds(&drawable.outline))
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn control_hulls_preserve_region_order_controls_degeneracy_and_signed_zeros() {
        let path = CurvePathOutput {
            path: BezPath::from_vec(vec![
                PathEl::MoveTo(Point::new(-0.0, 0.0)),
                PathEl::LineTo(Point::new(0.0, -0.0)),
                PathEl::ClosePath,
                PathEl::MoveTo(Point::new(10.0, 20.0)),
                PathEl::QuadTo(Point::new(-40.0, 80.0), Point::new(30.0, -60.0)),
                PathEl::CurveTo(
                    Point::new(-100.0, 200.0),
                    Point::new(300.0, -400.0),
                    Point::ZERO,
                ),
                PathEl::MoveTo(Point::new(7.0, 8.0)),
            ]),
        };
        let bounds = path.control_hull_bounds();
        assert_eq!(wire(&bounds), wire(&wire_bounds(&path)));
        assert_eq!(bounds.len(), 3);
        assert_eq!(bounds[0].left.to_bits(), (-0.0_f64).to_bits());
        assert_eq!(bounds[0].right.to_bits(), (-0.0_f64).to_bits());
        assert_eq!(bounds[0].bottom.to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            bounds[1],
            MarkBounds {
                left: -100.0,
                right: 300.0,
                bottom: -400.0,
                top: 200.0
            }
        );
        assert_eq!(
            bounds[2],
            MarkBounds {
                left: 7.0,
                right: 7.0,
                bottom: 8.0,
                top: 8.0
            }
        );
        assert!(
            CurvePathOutput {
                path: BezPath::new()
            }
            .control_hull_bounds()
            .is_empty()
        );
    }

    #[test]
    fn every_shape_default_stops_the_shaft_where_its_drawn_path_meets_the_line() {
        // Defaults at line thickness 2 (stroke radius 1): pointed length 15,
        // centred side 8, diamond diagonal 8√2, cross width 12.
        let l = 15.0;
        let straight = {
            // The chevron apex lies on the line, one stroke radius from the
            // origin measured across the arm.
            let arm = Point::new(1.0 - l, 0.4 * l - 1.0);
            let apex = mark("straight").prepare(CONTEXT).unwrap().end;
            let direction = (arm - Point::new(-apex, 0.0)).normalize();
            close(direction.cross(Vec2::new(apex, 0.0)).abs(), 1.0);
            apex
        };
        let contacts = [
            l,
            straight,
            0.6 * l,
            1.0 + 0.6 * (l - 2.0),
            1.0,
            1.0,
            (l / 2.0 - 1.0) / 2.0 + 1.0,
            0.0,
            1.0,
            3.0,
            3.0,
            3.0 * SQRT_2,
            0.0,
        ];
        for (shape, end) in SHAPES.into_iter().zip(contacts) {
            let prepared = mark(shape).prepare(CONTEXT).unwrap();
            close(prepared.end, end);
            assert_eq!(prepared.fit, MarkFit::Chord);
            close(prepared.shorten, 1.0);
            assert!(
                prepared
                    .paths
                    .iter()
                    .all(|p| p.path.path.is_finite() && p.outline.path.is_finite()),
                "{shape}"
            );
            assert!(!prepared.paths[0].path.path.is_empty(), "{shape}");
        }
    }

    #[test]
    fn stroked_heads_paint_their_nominal_boxes_in_both_directions() {
        // Line thickness 2: each stroke reaches 1 beyond its drawn path.
        let r = 4.0 * SQRT_2;
        for (shape, bounds) in [
            ("triangle", [-15.0, 0.0, 6.0]),
            ("stealth", [-15.0, 0.0, 6.0]),
            ("round", [-15.0, 0.0, 6.0]),
            ("straight", [-15.0, 0.0, 6.0]),
            ("tikz", [-8.5, 0.0, 7.5]),
            ("barb", [-7.5, 0.0, 7.5]),
            ("hooks", [-4.25, 0.0, 7.5]),
            ("bar", [-1.0, 1.0, 6.0]),
            ("bracket", [-3.6, 0.0, 6.0]),
            ("circle", [-4.0, 4.0, 4.0]),
            ("square", [-4.0, 4.0, 4.0]),
            ("diamond", [-r, r, r]),
        ] {
            let mut spec = mark(shape);
            spec.stroke = Some(true);
            let reversible = MarkShape::parse(shape).unwrap().options().contains(&"rev");
            for rev in [false, true].into_iter().filter(|rev| !rev || reversible) {
                spec.rev = reversible.then_some(rev);
                let prepared = spec.prepare(CONTEXT).unwrap();
                assert_bounds(&prepared, bounds, &format!("{shape} rev={rev}"));
            }
        }
        // Mitred right angles keep their points: the paint fills the outer
        // corners of a bracket's spine, either way round, and of a square.
        for (shape, rev, corners) in [
            ("bracket", false, vec![(0.0, 6.0), (0.0, -6.0)]),
            ("bracket", true, vec![(-3.6, 6.0), (-3.6, -6.0)]),
            (
                "square",
                false,
                vec![(4.0, 4.0), (-4.0, 4.0), (-4.0, -4.0), (4.0, -4.0)],
            ),
        ] {
            let mut spec = mark(shape);
            spec.stroke = Some(true);
            spec.rev = (shape == "bracket").then_some(rev);
            let outline = &spec.prepare(CONTEXT).unwrap().paths[0].outline.path;
            for (x, y) in corners {
                // A point just inside the corner, outside any bevel.
                let point = Point::new(x - 0.05 * f64::signum(x), y - 0.05 * f64::signum(y));
                assert_ne!(outline.winding(point), 0, "{shape} rev={rev} at {point:?}");
            }
        }
    }

    /// A stroked `shape` sized in points, where the shape takes the option.
    fn sized(shape: &str, length: f64, width: f64) -> MarkSpec {
        let mut spec = mark(shape);
        let options = MarkShape::parse(shape).unwrap().options();
        let size = |points| Some(MarkSize { points, ratio: 0.0 });
        if options.contains(&"length") {
            spec.length = size(length);
        }
        if options.contains(&"width") {
            spec.width = size(width);
        }
        spec.stroke = Some(true);
        spec
    }

    fn painted(spec: &MarkSpec, line_thickness: f64) -> PreparedMark {
        spec.prepare(MarkContext {
            line_thickness,
            ..CONTEXT
        })
        .unwrap()
    }

    /// How far the stroker's cubic round joins and caps may overshoot; the
    /// error grows with the stroke radius.
    fn round_tolerance(line_thickness: f64) -> f64 {
        5e-4 * line_thickness.max(2.0) / 2.0
    }

    fn assert_bounds(prepared: &PreparedMark, [left, right, half]: [f64; 3], label: &str) {
        let bounds = prepared.paths[0].outline.path.bounding_box();
        let tolerance = round_tolerance(prepared.line_thickness);
        for (actual, expected) in [
            (bounds.x0, left),
            (bounds.x1, right),
            (bounds.y0, -half),
            (bounds.y1, half),
        ] {
            assert!(
                (actual - expected).abs() < tolerance,
                "{label}: {actual} != {expected}"
            );
        }
    }

    #[test]
    fn slender_stroked_heads_keep_their_mitred_points() {
        // 7.5pt heads 2pt or 1.5pt wide on a 0.5pt line: their tips and wings
        // are sharper than a fixed miter limit of 7 would keep.
        for (shape, width) in [
            ("triangle", 2.0),
            ("triangle", 1.5),
            ("straight", 2.0),
            ("straight", 1.5),
            ("stealth", 2.0),
        ] {
            let prepared = painted(&sized(shape, 7.5, width), 1.0);
            assert_bounds(&prepared, [-15.0, 0.0, width], &format!("{shape} {width}"));
            assert!(prepared.paths[0].miter_limit > 7.0, "{shape} {width}");
        }
    }

    #[test]
    fn heads_too_small_for_their_stroke_grow_just_enough() {
        // A 3pt x 2.4pt triangle has an inradius of about 0.812pt: it holds
        // strokes up to 1.6pt (with 1% room) and grows about its tip beyond
        // that, keeping its shape and growing continuously with the stroke.
        let inradius = 2.0 * 14.4 / (4.8 + 2.0 * 41.76_f64.sqrt());
        let triangle = sized("triangle", 3.0, 2.4);
        assert_bounds(&painted(&triangle, 3.2), [-6.0, 0.0, 2.4], "holding");
        let mut previous = 6.0;
        for line in [3.26, 4.0, 6.0] {
            let bounds = painted(&triangle, line).paths[0]
                .outline
                .path
                .bounding_box();
            assert!(bounds.x1.abs() < 1e-9, "{line}: tip at {}", bounds.x1);
            close(bounds.y1 / -bounds.x0, 0.4);
            assert!(-bounds.x0 > previous, "{line}");
            previous = -bounds.x0;
        }
        let grown = painted(&triangle, 3.26).paths[0]
            .outline
            .path
            .bounding_box();
        close(-grown.x0 / 6.0, 1.01 * 1.63 / inradius);
        // Turned around, a grown triangle or stealth keeps its grown box.
        for shape in ["triangle", "stealth"] {
            let mut spec = sized(shape, 3.0, 2.4);
            let forward = painted(&spec, 6.0).paths[0].outline.path.bounding_box();
            spec.rev = Some(true);
            let bounds = [forward.x0, 0.0, forward.y1];
            assert_bounds(&painted(&spec, 6.0), bounds, &format!("{shape} rev"));
            assert!(forward.x1.abs() < 1e-9 && forward.x0 < -6.0, "{shape}");
        }
        // Other shapes grow only the dimensions too small to hold their cores
        // a stroke radius inside the box, and barb and hook arcs keep a radius
        // of at least a stroke radius, all with the same room.
        let room = 1.01 * 3.0;
        for (shape, bounds) in [
            ("straight", [-2.0 * room, 0.0, room]),
            ("round", [-2.0 * room, 0.0, room]),
            ("tikz", [-room - 3.0, 0.0, room]),
            ("barb", [-room - 3.0, 0.0, room + 3.0]),
            ("hooks", [-room - 3.0, 0.0, 2.0 * room + 3.0]),
            ("bracket", [-2.0 * room, 0.0, room]),
            ("circle", [-room, room, room]),
        ] {
            assert_bounds(&painted(&sized(shape, 3.0, 2.4), 6.0), bounds, shape);
        }
    }

    #[test]
    fn positive_area_offsets_are_scale_and_angle_aware() {
        let context = MarkContext {
            units_per_pt: 1.0,
            line_thickness: 1.0,
            shaft_stroke: MarkStrokeStyle::default(),
        };
        for (shape, length, width, inset) in [
            ("triangle", 1e-200, 1e-200, 0.0),
            ("triangle", f64::from_bits(1), f64::from_bits(1), 0.0),
            ("triangle", 1e-20, 1e-20, 0.0),
            ("triangle", 10.0, 1e-5, 0.0),
            ("stealth", 10.0, 8.0, 0.999999),
            ("stealth", 10.0, 8.0, 0.9999999),
            ("stealth", 10.0, 8.0, f64::from_bits(1.0_f64.to_bits() - 1)),
            ("stealth", 1e-200, 8e-201, 0.999999),
            ("square", 1e-20, 1e-20, 0.0),
            ("diamond", 1e-20, 1e-20, 0.0),
        ] {
            let mut spec = sized(shape, length, width);
            spec.fill = Some(false);
            if matches!(shape, "triangle" | "stealth") {
                spec.inset = Some(inset);
            }
            let prepared = spec.prepare(context).unwrap();
            let head = &prepared.paths[0];
            assert!(
                head.path.path.is_finite() && head.outline.path.is_finite(),
                "{shape}"
            );
            let bounds = head.outline.path.bounding_box();
            let centered = matches!(shape, "square" | "diamond");
            if !centered {
                assert!(
                    bounds.x1.abs() < 1e-8,
                    "{shape} {length}x{width} {inset}: {bounds:?}"
                );
                assert!(bounds.x0 < -length, "{shape} must grow: {bounds:?}");
                assert!(prepared.end > 0.0);
            }
            assert!(
                head.path.path.area() > 0.0,
                "positive polygon collapsed: {shape}"
            );
            // Uniform growth retains the intended polygon, rather than simply
            // adding a thick stroke to its nominal (possibly microscopic) edges.
            if shape == "triangle" && width < length {
                assert!(bounds.y1 > 0.5);
                assert!(head.miter_limit > 1e6);
                assert!((bounds.y1 / -bounds.x0 - width / (2.0 * length)).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn painted_miters_preserve_complete_grown_polygon_bounds() {
        // Independent lobe inradius: area / semiperimeter. No offset-normal
        // or stroker calculation is used to predict these outer vertices.
        for scale in [1e-6, 1.0, 1e6] {
            for (shape, width, inset) in [
                ("triangle", 1e-5, 0.0),
                ("triangle", 1e-7, 0.0),
                ("triangle", 1e-9, 0.0),
                ("stealth", 8.0, 0.999999),
                ("stealth", 8.0, 0.99999999),
            ] {
                let length: f64 = 10.0;
                let half = width / 2.0;
                let radius = half * length * (1.0 - inset)
                    / (length.hypot(half) + (length * inset).hypot(half));
                let grow = (1.01 * 0.5 / radius).max(1.0);
                let expected = [-grow * length * scale, 0.0, grow * half * scale];
                let mut spec = sized(shape, length * scale, width * scale);
                spec.inset = (shape == "stealth").then_some(inset);
                spec.fill = Some(false);
                let prepared = spec
                    .prepare(MarkContext {
                        units_per_pt: 1.0,
                        line_thickness: scale,
                        shaft_stroke: MarkStrokeStyle::default(),
                    })
                    .unwrap();
                let paint = &prepared.paths[0].outline.path;
                let bounds = paint.bounding_box();
                let label = format!("{shape} width={width} inset={inset} scale={scale}");
                // Subtracting the grown tip's coordinates has absolute error
                // proportional to its extent. The wing expectation additionally
                // accounts for edge-coordinate rounding amplified by the
                // reciprocal lobe angle near a 100% notch. This is separate
                // from the exact-emitted-edge high-precision oracle.
                let tip_error = 128.0 * f64::EPSILON * expected[0].abs();
                let wing_error = 32.0 * f64::EPSILON / (1.0 - inset);
                assert!(bounds.x1.abs() <= tip_error, "{label}: {bounds:?}");
                for (actual, wanted) in [
                    (bounds.x0, expected[0]),
                    (bounds.y1, expected[2]),
                    (bounds.y0, -expected[2]),
                ] {
                    assert!(
                        (actual - wanted).abs() <= wing_error * wanted.abs() + tip_error,
                        "{label}: {actual} versus {wanted}"
                    );
                }
                // A real interior point near the grown wing must be painted;
                // matching a metadata rectangle alone cannot satisfy this.
                let near_wing = Point::new(
                    expected[0] * (0.99 + 0.005 * (1.0 - inset)),
                    expected[2] * 0.99,
                );
                assert_ne!(paint.winding(near_wing), 0, "{label}: missing wing");
                assert!(prepared.end().is_finite() && prepared.end() > 0.0);
            }
        }
    }

    #[test]
    fn actual_zero_area_polygons_keep_the_as_given_exception() {
        for (length, width, inset) in [(0.0, 8.0, 0.4), (10.0, 0.0, 0.4), (10.0, 8.0, 1.0)] {
            let polygon = Polygon::dart(Point::ORIGIN, length, width, inset);
            assert!(polygon.degenerate());
            let expected = polygon.elements().collect::<BezPath>();
            let actual = polygon.stroked_head(0.5, length).unwrap();
            assert_eq!(actual.path.elements(), expected.elements());
            assert!(actual.path.is_finite());
        }
        // A finite positive input whose aspect ratio requires an infinite
        // miter must report numerical inability, never the zero-area exception.
        let impossible = Polygon::dart(Point::ORIGIN, 1.0, f64::from_bits(1), 0.0);
        assert!(!impossible.degenerate());
        assert!(
            impossible
                .stroked_head(0.5, 1.0)
                .unwrap_err()
                .contains("unrepresentable")
        );
        // Zero stroke leaves even an extremely small positive polygon alone.
        let polygon = Polygon::dart(Point::ORIGIN, 1e-200, 1e-200, 0.0);
        assert!(!polygon.degenerate());
        let expected = polygon.elements().collect::<BezPath>();
        assert_eq!(polygon.stroked_head(0.0, 1e-200).unwrap().path, expected);
    }

    #[test]
    fn barb_contacts_follow_all_axis_crossings_and_signed_sweep_boundaries() {
        let almost_turn = f64::from_bits(TAU.to_bits() - 1);
        let context = MarkContext {
            units_per_pt: 1.0,
            line_thickness: 1.0,
            shaft_stroke: MarkStrokeStyle::default(),
        };
        for sweep in [
            -TAU,
            -almost_turn,
            -1.5 * PI,
            -PI,
            -1e-10,
            0.0,
            1e-10,
            PI,
            1.5 * PI,
            almost_turn,
            TAU,
        ] {
            for rev in [false, true] {
                let mut spec = sized("barb", 0.0, 10.0);
                spec.arc = Some(sweep);
                spec.rev = Some(rev);
                let full = sweep.abs() == TAU;
                let expected = if rev {
                    4.5
                } else if full {
                    9.5
                } else {
                    0.5
                };
                let prepared = spec.prepare(context).unwrap();
                close(prepared.end(), expected);
                let carrier = CurvePathOutput {
                    path: line(Point::ZERO, Point::new(100.0, 0.0)),
                };
                let placement = MarkPlacement {
                    template: 0,
                    carrier: 0,
                    station: MarkStation::End,
                    direction: MarkDirection::Forward,
                    shift: 0.0,
                };
                let placed = prepared.place(&carrier, &placement).unwrap();
                point_close(placed.back.into(), Point::new(100.0 - expected, 0.0));
                point_close(
                    placed.shaft_contact.into(),
                    Point::new(100.0 - expected, 0.0),
                );
                if full {
                    // The full circular drawn path's leftmost axis crossing is
                    // its rear. This also checks the existing reversal frame.
                    close(prepared.paths[0].path.path.bounding_box().x0, -expected);
                }
                if full {
                    let curved = CurvePathOutput { path: curve() };
                    let on_curve = prepared.place(&curved, &placement).unwrap();
                    let tip: Point = on_curve.tip.into();
                    let back: Point = on_curve.back.into();
                    let shaft: Point = on_curve.shaft_contact.into();
                    close(tip.distance(back), expected);
                    point_close(back, shaft);
                }
                let composite = MarkSpec {
                    shape: "combine".into(),
                    parts: Some(vec![
                        MarkPart::Mark(Box::new(spec)),
                        MarkPart::Mark(Box::new(mark("bar"))),
                    ]),
                    ..mark("combine")
                }
                .prepare(context)
                .unwrap();
                close(composite.end(), expected);
                let PathEl::MoveTo(bar) = composite.paths[1].path.path.elements()[0] else {
                    panic!("bar starts with move")
                };
                close(bar.x, -expected);
            }
        }
    }

    #[test]
    fn deep_notches_grow_rather_than_turn_inside_out() {
        // The tip, notch and upper wing corner of a 10 x 8 stealth notched 90%
        // deep meet on the line at an inset of about 0.194. Beyond it the
        // moved edges reverse, so a thicker stroke grows the head instead.
        let mut stealth = sized("stealth", 5.0, 4.0);
        stealth.inset = Some(0.9);
        assert_bounds(&painted(&stealth, 0.36), [-10.0, 0.0, 4.0], "holding");
        let prepared = painted(&stealth, 0.4);
        let drawn = &prepared.paths[0].path.path;
        let PathEl::MoveTo(tip) = drawn.elements()[0] else {
            panic!("outline starts at the tip")
        };
        assert!(drawn.area() > 0.0, "drawn outline turned inside out");
        assert!(
            -prepared.end < tip.x,
            "notch {} ahead of tip {tip:?}",
            -prepared.end
        );
        let bounds = prepared.paths[0].outline.path.bounding_box();
        assert!(bounds.x1.abs() < 1e-9 && bounds.x0 < -10.0, "{bounds:?}");
        close(bounds.y1 / -bounds.x0, 0.4);
    }

    #[test]
    fn stroked_heads_of_any_size_stay_finite_and_behind_the_tip() {
        // Only a polygon without area, here one of zero length, cannot hold a
        // stroke and is stroked as given, across the tip. A head turned around
        // within its span, grown or not, stays behind the tip as well.
        for shape in SHAPES {
            let reversible = MarkShape::parse(shape).unwrap().options().contains(&"rev");
            let centred = matches!(shape, "bar" | "circle" | "square" | "diamond" | "rays");
            for (length, width) in [
                (3.0, 2.4),
                (3.0, 0.6),
                (0.5, 0.5),
                (0.0, 1.0),
                (1.0, 0.0),
                (7.5, 1.5),
            ] {
                let flat = matches!(shape, "triangle" | "stealth") && length == 0.0;
                for line in [0.0, 1.0, 3.26, 6.0, 12.0] {
                    for rev in [false, true].into_iter().filter(|rev| !rev || reversible) {
                        let label = format!("{shape} {length}x{width} at {line} rev={rev}");
                        let mut spec = sized(shape, length, width);
                        spec.rev = reversible.then_some(rev);
                        let prepared = painted(&spec, line);
                        let outline = &prepared.paths[0].outline.path;
                        assert!(outline.is_finite() && prepared.end.is_finite(), "{label}");
                        // A mitred tip lands exactly on the line end; a turned
                        // head's rear there may be a rounded cap or join.
                        let tolerance = if rev { round_tolerance(line) } else { 1e-6 };
                        if !centred && !flat && !outline.is_empty() {
                            let front = outline.bounding_box().x1;
                            assert!(front < tolerance, "{label}: front at {front}");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn mixed_sizes_resolve_in_drawing_units_and_width_defaults_follow_length() {
        close(
            MarkSize {
                points: 3.0,
                ratio: 4.5,
            }
            .resolve(CONTEXT)
            .unwrap(),
            15.0,
        );
        for shape in SHAPES {
            let mut spec = mark(shape);
            if matches!(shape, "tikz" | "barb" | "hooks" | "bar") {
                spec.width = Some(MarkSize {
                    points: 3.0,
                    ratio: 2.0,
                });
            } else {
                spec.length = Some(MarkSize {
                    points: 3.0,
                    ratio: 2.0,
                });
                if shape != "rays" {
                    spec.width = Some(MarkSize {
                        points: 2.0,
                        ratio: 1.0,
                    });
                }
            }
            let prepared = spec.prepare(CONTEXT).unwrap();
            assert!(prepared.end.is_finite(), "{shape}");
            assert!(prepared.paths[0].path.path.is_finite(), "{shape}");
        }
        let mut triangle = mark("triangle");
        triangle.length = Some(MarkSize {
            points: 3.0,
            ratio: 2.0,
        });
        let prepared = triangle.prepare(CONTEXT).unwrap();
        close(prepared.end, 10.0);
        assert_eq!(
            prepared.paths[0].path.path.elements()[1],
            PathEl::LineTo(Point::new(-10.0, 4.0))
        );
    }

    #[test]
    fn triangle_and_stealth_use_actual_notch_not_half_length() {
        for (shape, expected) in [("triangle", 15.0), ("stealth", 9.0)] {
            let prepared = mark(shape).prepare(CONTEXT).unwrap();
            close(prepared.end, expected);
            assert!(!prepared.paths[0].stroke);
            assert!(prepared.paths[0].fill);
            let mut reverse = mark(shape);
            reverse.rev = Some(true);
            close(reverse.prepare(CONTEXT).unwrap().end, 15.0);
        }
    }

    #[test]
    fn explicit_triangle_stroke_runs_inside_the_outline_and_stops_the_shaft_beneath_it() {
        let mut spec = mark("triangle");
        spec.stroke = Some(true);
        let normal = spec.prepare(CONTEXT).unwrap();
        // The drawn outline is the nominal one moved one stroke radius inward:
        // the tip by 1/sin of its half angle, the edges by 1.
        let q = 1.0 / 0.4_f64.atan().sin();
        let elements = normal.paths[0].path.path.elements();
        let PathEl::MoveTo(tip) = elements[0] else {
            panic!("outline starts at the tip")
        };
        point_close(tip, Point::new(-q, 0.0));
        let PathEl::LineTo(wing) = elements[1] else {
            panic!("outline continues to the upper wing")
        };
        close(wing.x, -14.0);
        close(Vec2::new(-15.0, 6.0).normalize().cross(wing.to_vec2()), 1.0);
        assert_eq!(elements[2], PathEl::LineTo(Point::new(-14.0, 0.0)));
        // The shaft stops on the drawn back edge, under the stroke.
        close(normal.end, 14.0);
        assert_eq!(normal.paths[0].join, "miter");
        spec.rev = Some(true);
        close(spec.prepare(CONTEXT).unwrap().end, 15.0 - q);
    }

    #[test]
    fn tikz_arms_are_quarter_circles_meeting_in_a_cusp_and_arc_is_rejected() {
        let prepared = mark("tikz").prepare(CONTEXT).unwrap();
        // Width 15 and stroke radius 1: arms of radius 6.5 about (-1, ±6.5),
        // one cubic per quarter turn, tangent to the line at the cusp.
        let drawable = &prepared.paths[0];
        let segments: Vec<_> = drawable.path.path.segments().collect();
        assert_eq!(segments.len(), 2);
        for (segment, center) in segments
            .iter()
            .zip([Point::new(-1.0, 6.5), Point::new(-1.0, -6.5)])
        {
            for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
                assert!((segment.eval(t).distance(center) - 6.5).abs() < 2e-3);
            }
        }
        let PathSeg::Cubic(upper) = segments[0] else {
            panic!("upper arm is a cubic")
        };
        let PathSeg::Cubic(lower) = segments[1] else {
            panic!("lower arm is a cubic")
        };
        point_close(upper.p0, Point::new(-7.5, 6.5));
        point_close(upper.p3, Point::new(-1.0, 0.0));
        point_close(lower.p3, Point::new(-7.5, -6.5));
        close(upper.p2.y, 0.0);
        close(lower.p1.y, 0.0);
        assert!(upper.p2.x < upper.p3.x && lower.p1.x < lower.p0.x);
        close(prepared.end, 1.0);
        assert_eq!(
            (drawable.join.as_str(), drawable.cap.as_str()),
            ("round", "round")
        );
        let mut spec = mark("tikz");
        spec.arc = Some(PI);
        let bytes = wire(&spec);
        let err = ciborium::de::from_reader::<MarkSpec, _>(bytes.as_slice())
            .unwrap_err()
            .to_string();
        assert!(err.contains("tikz.arc"), "{err}");
        let mut combine = mark("combine");
        combine.parts = Some(vec![MarkPart::Mark(Box::new(spec))]);
        let err = ciborium::de::from_reader::<MarkSpec, _>(wire(&combine).as_slice())
            .unwrap_err()
            .to_string();
        assert!(err.contains("tikz.arc"), "{err}");
    }

    #[test]
    fn shape_inapplicable_and_unknown_fields_are_rejected_during_deserialization() {
        for (shape, field, value) in [
            ("tikz", "rev", Value::Bool(false)),
            ("straight", "fill", Value::Bool(true)),
            ("bar", "length", Value::Map(vec![])),
            ("rays", "width", Value::Map(vec![])),
            ("circle", "rev", Value::Bool(false)),
            ("hooks", "align", Value::Text("end".into())),
            ("triangle", "extra", Value::Bool(false)),
            ("triangle", "arc", Value::Null),
        ] {
            let value = Value::Map(vec![
                (Value::Text("shape".into()), Value::Text(shape.into())),
                (Value::Text(field.into()), value),
            ]);
            assert!(
                ciborium::de::from_reader::<MarkSpec, _>(wire(&value).as_slice()).is_err(),
                "{shape}.{field}"
            );
        }
        let value = Value::Map(vec![(
            Value::Text("shape".into()),
            Value::Text("unknown".into()),
        )]);
        assert!(ciborium::de::from_reader::<MarkSpec, _>(wire(&value).as_slice()).is_err());
    }

    #[test]
    fn finite_validation_covers_context_sizes_angles_shorten_and_stations() {
        let path = line(Point::ZERO, Point::new(100.0, 0.0));
        for context in [
            MarkContext {
                units_per_pt: 0.0,
                ..CONTEXT
            },
            MarkContext {
                units_per_pt: f64::NAN,
                ..CONTEXT
            },
            MarkContext {
                line_thickness: -1.0,
                ..CONTEXT
            },
            MarkContext {
                line_thickness: f64::INFINITY,
                ..CONTEXT
            },
        ] {
            let mut request = batch(mark("triangle"), path.clone());
            request.templates[0].context = context;
            assert!(request.geometry().is_err());
        }
        for value in [f64::NAN, f64::INFINITY, -1.0] {
            let mut spec = mark("triangle");
            spec.shorten = Some(value);
            assert!(batch(spec, path.clone()).geometry().is_err());
        }
        let mut spec = mark("triangle");
        spec.length = Some(MarkSize {
            points: f64::NAN,
            ratio: 0.0,
        });
        assert!(batch(spec, path.clone()).geometry().is_err());
        for shape in ["barb", "hooks"] {
            let mut spec = mark(shape);
            spec.arc = Some(f64::NAN);
            assert!(batch(spec, path.clone()).geometry().is_err());
        }
        let mut request = batch(mark("rays"), path);
        request.placements[0].station = MarkStation::Ratio { value: f64::NAN };
        assert!(request.geometry().is_err());
    }

    #[test]
    fn round_full_inset_is_meaningful_and_out_of_range_is_rejected() {
        let mut spec = mark("round");
        spec.inset = Some(1.0);
        close(spec.prepare(CONTEXT).unwrap().end, 1.0);
        for inset in [-0.1, 1.1, f64::NAN] {
            spec.inset = Some(inset);
            assert!(spec.prepare(CONTEXT).is_err());
        }
    }

    #[test]
    fn align_and_rays_phase_defaults_match_local_geometry() {
        // `align: end` puts the front of a bar, circle, square, or diamond's
        // paint on the tip, with the shaft stopping at its drawn rear, and the
        // rays' hub one ray length (5.6) behind the tip.
        for (shape, end) in [
            ("bar", 1.0),
            ("circle", 7.0),
            ("square", 7.0),
            ("diamond", 7.0 * SQRT_2),
            ("rays", 5.6),
        ] {
            let mut spec = mark(shape);
            spec.align = Some(MarkAlign::End);
            let aligned = spec.prepare(CONTEXT).unwrap();
            close(aligned.end, end);
            if shape == "rays" {
                for element in aligned.paths[0].path.path.elements() {
                    if let PathEl::MoveTo(hub) = element {
                        point_close(*hub, Point::new(-5.6, 0.0));
                    }
                }
            } else {
                let front = aligned.paths[0].outline.path.bounding_box().x1;
                assert!(front.abs() < 5e-4, "{shape}: front at {front}");
            }
        }
        // No default spoke runs along the line, where the shaft would hide
        // it: even counts start at π/n and mirror about the line, odd counts
        // start perpendicular to it and mirror about that perpendicular.
        for count in 1..=8 {
            let mut spec = mark("rays");
            spec.n = Some(count);
            let prepared = spec.prepare(CONTEXT).unwrap();
            let ends: Vec<_> = prepared.paths[0]
                .path
                .path
                .elements()
                .iter()
                .filter_map(|element| match element {
                    PathEl::LineTo(point) => Some(*point),
                    _ => None,
                })
                .collect();
            assert_eq!(ends.len(), count);
            let first = if count.is_multiple_of(2) {
                PI / count as f64
            } else {
                -FRAC_PI_2
            };
            point_close(ends[0], Point::new(5.6 * first.cos(), 5.6 * first.sin()));
            let mirror = |end: Point| {
                if count.is_multiple_of(2) {
                    Point::new(end.x, -end.y)
                } else {
                    Point::new(-end.x, end.y)
                }
            };
            for end in &ends {
                close(end.to_vec2().hypot(), 5.6);
                assert!(end.y.abs() > 1e-9, "{count}: {end:?}");
                assert!(ends.iter().any(|other| other.distance(mirror(*end)) < 1e-9));
            }
        }
    }

    #[test]
    fn reverse_offsets_are_shape_specific() {
        for shape in ["barb", "hooks", "bracket", "straight", "round"] {
            let normal = mark(shape).prepare(CONTEXT).unwrap();
            let mut spec = mark(shape);
            spec.rev = Some(true);
            let reversed = spec.prepare(CONTEXT).unwrap();
            let offset = match shape {
                "barb" => 7.5,
                "hooks" => 4.25,
                "bracket" => 0.3 * (2.4 * 2.0 + 3.6 * 2.0),
                _ => 15.0,
            };
            let mut expected = normal.paths[0].path.path.clone();
            expected.apply_affine(Affine::new([-1.0, 0.0, 0.0, 1.0, -offset, 0.0]));
            assert_eq!(expected, reversed.paths[0].path.path, "{shape}");
        }
    }

    #[test]
    fn combined_marks_accumulate_contacts_and_signed_mixed_gaps_once() {
        let mut nested = mark("combine");
        let child = mark("triangle");
        nested.parts = Some(vec![
            MarkPart::Mark(Box::new(child)),
            MarkPart::Gap(MarkGap {
                gap: MarkSize {
                    points: -1.0,
                    ratio: 0.5,
                },
            }),
            MarkPart::Mark(Box::new(mark("stealth"))),
        ]);
        let mut spec = mark("combine");
        spec.parts = Some(vec![
            MarkPart::Mark(Box::new(nested)),
            MarkPart::Mark(Box::new(mark("bar"))),
        ]);
        let prepared = spec.prepare(CONTEXT).unwrap();
        close(prepared.end, 23.0);
        assert_eq!(prepared.paths.len(), 3);
        assert_eq!(
            prepared.paths[1].path.path.elements()[0],
            PathEl::MoveTo(Point::new(-14.0, 0.0))
        );
        let result = geometry(spec, line(Point::ZERO, Point::new(100.0, 0.0)));
        point_close(result.back.into(), Point::new(77.0, 0.0));
        assert_eq!(
            result.shaft.path.elements()[1],
            PathEl::LineTo(Point::new(77.0, 0.0))
        );
    }

    #[test]
    fn straight_chord_shortening_zero_and_one_reuse_visible_carrier() {
        let path = line(Point::new(10.0, 0.0), Point::new(110.0, 0.0));
        let mut spec = mark("triangle");
        let result = geometry(spec.clone(), path.clone());
        point_close(result.tip.into(), Point::new(110.0, 0.0));
        point_close(result.back.into(), Point::new(95.0, 0.0));
        assert_eq!(
            result.shaft.path.elements()[1],
            PathEl::LineTo(Point::new(95.0, 0.0))
        );
        spec.shorten = Some(0.0);
        let unshortened = geometry(spec, path.clone());
        assert_eq!(unshortened.shaft.path, path);
        assert_eq!(wire(&unshortened.paths), wire(&result.paths));
    }

    #[test]
    fn curved_chord_contacts_are_on_curve_and_separated_by_geometric_end() {
        let path = curve();
        let result = geometry(mark("triangle"), path.clone());
        let carrier = PreparedCarrier::new(path).unwrap();
        let (contact, at) = carrier.chord_contact(carrier.length, 15.0, 1.0).unwrap();
        point_close(result.back.into(), contact);
        point_close(carrier.point(at).unwrap(), contact);
        close(Point::from(result.tip).distance(contact), 15.0);
        assert!(carrier.length - at > 15.0);
        let last = result.shaft.path.segments().last().unwrap();
        point_close(last.end(), contact);
    }

    #[test]
    fn bend_retains_original_tip_tangent_and_controls_even_with_short_last_segment() {
        let mut path = curve();
        path.curve_to((60.0, -0.1), (60.0, -0.2), (60.0, -0.3));
        let original = path.elements().to_vec();
        let mut spec = mark("triangle");
        spec.fit = Some(MarkFit::Bend);
        let result = geometry(spec, path);
        point_close(result.tip.into(), Point::new(60.0, -0.3));
        let PathEl::CurveTo(c1, c2, end) = result.shaft.path.elements().last().unwrap() else {
            panic!("last cubic");
        };
        assert_eq!(*c1, Point::new(60.0, -0.1));
        assert_eq!(*c2, Point::new(60.0, -0.2));
        point_close(*end, Point::new(60.0, 14.7));
        assert_eq!(
            &result.shaft.path.elements()[..original.len() - 1],
            &original[..original.len() - 1]
        );
    }

    #[test]
    fn start_backward_retracts_inward_and_interior_bend_always_fits_chord() {
        for fit in [MarkFit::Chord, MarkFit::Bend] {
            let mut spec = mark("triangle");
            spec.fit = Some(fit);
            let mut request = batch(spec, line(Point::ZERO, Point::new(100.0, 0.0)));
            request.placements[0].station = MarkStation::Start;
            request.placements[0].direction = MarkDirection::Backward;
            let result = request.geometry().unwrap().marks.remove(0);
            point_close(result.back.into(), Point::new(15.0, 0.0));
            assert_eq!(
                result.shaft.path.elements()[0],
                PathEl::MoveTo(Point::new(15.0, 0.0))
            );
        }
        let mut spec = mark("triangle");
        spec.fit = Some(MarkFit::Bend);
        let mut request = batch(spec, curve());
        request.placements[0].station = MarkStation::Ratio { value: 0.5 };
        let bend = request.clone().geometry().unwrap();
        request.templates[0].mark.fit = Some(MarkFit::Chord);
        let chord = request.geometry().unwrap();
        assert_eq!(wire(&bend), wire(&chord));
    }

    #[test]
    fn ratio_distance_shift_and_clamping_use_original_visible_arc_stations() {
        let mut request = batch(
            mark("triangle"),
            line(Point::new(10.0, 0.0), Point::new(110.0, 0.0)),
        );
        request.placements[0].station = MarkStation::Ratio { value: 0.5 };
        request.placements[0].shift = -5.0;
        let ratio = request.clone().geometry().unwrap();
        let head = &ratio.marks[0];
        close((head.tip.x + head.back.x) / 2.0, 55.0);
        close(head.tip.x - head.back.x, head.end);
        request.placements[0].station = MarkStation::Distance { value: 50.0 };
        assert_eq!(wire(&ratio), wire(&request.clone().geometry().unwrap()));
        request.placements[0].shift = 1000.0;
        close(request.geometry().unwrap().marks[0].tip.x, 110.0);
    }

    #[test]
    fn numeric_stations_center_painted_catalogue_and_composites_in_both_directions() {
        let mut specs: Vec<_> = SHAPES.into_iter().map(mark).collect();
        let mut composite = mark("combine");
        composite.parts = Some(vec![
            MarkPart::Mark(Box::new(mark("circle"))),
            MarkPart::Gap(MarkGap {
                gap: MarkSize {
                    points: 2.0,
                    ratio: 0.5,
                },
            }),
            MarkPart::Mark(Box::new(mark("triangle"))),
        ]);
        specs.push(composite);
        for mut spec in specs {
            let supports_rev = matches!(
                spec.shape.as_str(),
                "triangle" | "straight" | "stealth" | "round" | "barb" | "hooks" | "bracket"
            );
            for rev in [false, true].into_iter().filter(|rev| !rev || supports_rev) {
                if supports_rev {
                    spec.rev = Some(rev);
                }
                for direction in [MarkDirection::Forward, MarkDirection::Backward] {
                    for reversed_carrier in [false, true] {
                        let (start, end) = if reversed_carrier {
                            (110.0, 10.0)
                        } else {
                            (10.0, 110.0)
                        };
                        let mut request = batch(
                            spec.clone(),
                            line(Point::new(start, 0.0), Point::new(end, 0.0)),
                        );
                        request.mode = MarkGeometryMode::Selected;
                        request.templates[0].context.line_thickness = 3.0;
                        let placement = &mut request.placements[0];
                        placement.station = MarkStation::Ratio { value: 0.5 };
                        placement.direction = direction;
                        placement.shift = -4.0;
                        let ratio = request.clone().geometry().unwrap();
                        let head = &ratio.marks[0];
                        let bounds = head.footprint.path.bounding_box();
                        let expected_center = if reversed_carrier { 64.0 } else { 56.0 };
                        assert!(
                            ((bounds.x0 + bounds.x1) / 2.0 - expected_center).abs() < 1e-12,
                            "{} rev={rev} direction={direction:?}",
                            spec.shape
                        );
                        close((head.tip.x - head.back.x).abs(), head.end.abs());
                        assert_eq!(ratio.shafts[0].shaft.path, request.carriers[0].path);
                        request.placements[0].station = MarkStation::Distance { value: 50.0 };
                        assert_eq!(wire(&ratio), wire(&request.geometry().unwrap()));
                    }
                }
            }
        }
    }

    #[test]
    fn cubicized_straight_carriers_preserve_numeric_painted_centers() {
        for shape in SHAPES.into_iter().chain(["combine"]) {
            for direction in [MarkDirection::Forward, MarkDirection::Backward] {
                for (start, end, ratio, shift, expected) in
                    [(0.0, 3.0, 0.5, -0.07, 1.43), (5.0, 10.0, 0.2, 0.0, 6.0)]
                {
                    let mut path = BezPath::new();
                    path.move_to((start, 0.0));
                    // Typst's line carrier is represented by a cubic with
                    // nonuniform parameter speed, not a native Line segment.
                    path.curve_to((start, 0.0), (end, 0.0), (end, 0.0));
                    let mut request = batch(mark(shape), path.clone());
                    request.mode = MarkGeometryMode::Selected;
                    request.templates[0].context = MarkContext {
                        units_per_pt: 2.54 / 72.0,
                        line_thickness: 0.8 * 2.54 / 72.0,
                        ..CONTEXT
                    };
                    request.placements[0].station = MarkStation::Ratio { value: ratio };
                    request.placements[0].shift = shift;
                    request.placements[0].direction = direction;
                    let result = request.geometry().unwrap();
                    let bounds = result.marks[0].footprint.path.bounding_box();
                    assert!(
                        ((bounds.x0 + bounds.x1) / 2.0 - expected).abs() < 1e-12,
                        "{shape} {direction:?}"
                    );
                    assert_eq!(result.shafts[0].shaft.path, path);
                }
            }
        }
    }

    #[test]
    fn numeric_endpoint_stations_keep_chord_fit_and_clamp_origin_without_resizing() {
        for value in [0.0, 1.0] {
            let mut request = batch(mark("triangle"), curve());
            request.placements[0].station = MarkStation::Ratio { value };
            request.templates[0].mark.fit = Some(MarkFit::Bend);
            let bend = request.clone().geometry().unwrap();
            request.templates[0].mark.fit = Some(MarkFit::Chord);
            assert_eq!(wire(&bend), wire(&request.geometry().unwrap()));
        }
        for direction in [MarkDirection::Forward, MarkDirection::Backward] {
            let mut request = batch(mark("triangle"), line(Point::ZERO, Point::new(2.0, 0.0)));
            request.placements[0].station = MarkStation::Ratio { value: 0.5 };
            request.placements[0].direction = direction;
            let head = request.geometry().unwrap().marks.remove(0);
            let (tip, back) = if direction == MarkDirection::Forward {
                (2.0, -13.0)
            } else {
                (0.0, 15.0)
            };
            close(head.tip.x, tip);
            close(head.back.x, back);
            close(head.end, 15.0);
            assert!(head.shaft_contact.x >= 0.0 && head.shaft_contact.x <= 2.0);
        }
    }

    #[test]
    fn short_zero_empty_and_reversed_carriers_have_finite_painted_bounds() {
        for path in [
            line(Point::ZERO, Point::new(0.01, 0.0)),
            line(Point::new(7.0, 2.0), Point::new(7.0, 2.0)),
            BezPath::new(),
            line(Point::new(100.0, 0.0), Point::ZERO),
        ] {
            for shape in SHAPES {
                for fit in [MarkFit::Chord, MarkFit::Bend] {
                    let mut spec = mark(shape);
                    spec.fit = Some(fit);
                    let result = geometry(spec, path.clone());
                    assert!(result.shaft.path.is_finite(), "{shape}");
                    assert!(result.footprint.path.is_finite(), "{shape}");
                    assert!(result.footprint.path.bounding_box().is_finite(), "{shape}");
                    assert!(result.tip.x.is_finite() && result.back.x.is_finite());
                }
            }
        }
    }

    #[test]
    fn footprints_include_exact_painted_head_and_shaft_geometry() {
        for shape in SHAPES {
            let result = geometry(mark(shape), curve());
            let bounds = result.footprint.path.bounding_box();
            for drawable in &result.paths {
                let painted = drawable.outline.path.bounding_box();
                assert!(bounds.x0 <= painted.x0 && bounds.x1 >= painted.x1);
                assert!(bounds.y0 <= painted.y0 && bounds.y1 >= painted.y1);
            }
            let shaft = kurbo::stroke(
                result.shaft.path.elements().iter().copied(),
                &result.shaft_style.stroke(CONTEXT.line_thickness),
                &StrokeOpts::default(),
                1e-6,
            )
            .bounding_box();
            assert!(bounds.x0 <= shaft.x0 && bounds.x1 >= shaft.x1);
            assert!(bounds.y0 <= shaft.y0 && bounds.y1 >= shaft.y1);
        }
    }

    #[test]
    fn cbor_roundtrip_matches_native_and_batches_match_scalar_engine() {
        let mut request = batch(mark("triangle"), curve());
        request.templates = SHAPES
            .into_iter()
            .map(|shape| MarkTemplateSpec {
                mark: mark(shape),
                context: CONTEXT,
            })
            .collect();
        request.placements = (0..SHAPES.len())
            .map(|template| MarkPlacement {
                template,
                carrier: 0,
                station: MarkStation::End,
                shift: 0.0,
                direction: MarkDirection::Forward,
            })
            .collect();
        let native = request.clone().geometry().unwrap();
        let encoded = mark_geometry_bytes(&wire(&request)).unwrap();
        let decoded: MarkGeometryBatchOutput =
            ciborium::de::from_reader(encoded.as_slice()).unwrap();
        assert_eq!(encoded, wire(&native));
        assert_eq!(encoded, wire(&decoded));
        for (placement, expected) in request.placements.iter().zip(native.marks) {
            let template = &request.templates[placement.template];
            let prepared = template.mark.prepare(template.context).unwrap();
            close(prepared.end(), expected.end);
            assert_eq!(prepared.paths().len(), expected.paths.len());
            let scalar = prepared
                .place(&request.carriers[placement.carrier], placement)
                .unwrap();
            assert_eq!(wire(&scalar), wire(&expected));
        }
    }

    #[test]
    fn zero_line_thickness_and_invalid_batch_indices_are_handled() {
        for shape in SHAPES {
            let mut request = batch(mark(shape), curve());
            request.templates[0].context.line_thickness = 0.0;
            assert!(
                request.geometry().unwrap().marks[0]
                    .footprint
                    .path
                    .is_finite()
            );
        }
        let mut request = batch(mark("triangle"), curve());
        request.placements[0].template = 1;
        assert!(request.geometry().is_err());
        let mut request = batch(mark("triangle"), curve());
        request.placements[0].carrier = 1;
        assert!(request.geometry().is_err());
    }

    #[test]
    fn short_carriers_preserve_full_size_and_separate_geometric_and_shaft_contacts() {
        let path = line(Point::ZERO, Point::new(2.0, 0.0));
        for fit in [MarkFit::Chord, MarkFit::Bend] {
            let mut spec = mark("triangle");
            spec.fit = Some(fit);
            let result = geometry(spec, path.clone());
            close(result.end, 15.0);
            point_close(result.tip.into(), Point::new(2.0, 0.0));
            point_close(result.back.into(), Point::new(-13.0, 0.0));
            close(result.shaft_contact.x, 0.0);
            let bounds = result.footprint.path.bounding_box();
            assert!(bounds.x0 <= -13.0 && bounds.x1 >= 2.0);
            assert_eq!(
                result.paths[0].path.path.elements()[1],
                PathEl::LineTo(Point::new(-13.0, 6.0))
            );
        }
        let normal = geometry(mark("triangle"), line(Point::ZERO, Point::new(100.0, 0.0)));
        point_close(normal.back.into(), normal.shaft_contact.into());
    }

    #[test]
    fn one_batch_supports_distinct_line_contexts_and_candidate_footprints() {
        let mut request = batch(mark("triangle"), curve());
        request.templates.push(MarkTemplateSpec {
            mark: mark("triangle"),
            context: MarkContext {
                line_thickness: 4.0,
                ..CONTEXT
            },
        });
        request.placements.push(MarkPlacement {
            template: 1,
            carrier: 0,
            station: MarkStation::Ratio { value: 0.7 },
            shift: 0.0,
            direction: MarkDirection::Backward,
        });
        request.placements.push(request.placements[0].clone());
        let result = request.clone().geometry().unwrap().marks;
        close(result[0].end, 15.0);
        close(result[1].end, 24.0);
        assert_eq!(wire(&result[0]), wire(&result[2]));
        for (placement, expected) in request.placements.iter().zip(result) {
            let template = &request.templates[placement.template];
            let scalar = template
                .mark
                .prepare(template.context)
                .unwrap()
                .place(&request.carriers[placement.carrier], placement)
                .unwrap();
            assert_eq!(wire(&scalar.footprint), wire(&expected.footprint));
        }
    }

    #[test]
    fn gaps_can_make_signed_contacts_but_must_remain_finite() {
        let mut spec = mark("combine");
        spec.parts = Some(vec![
            MarkPart::Mark(Box::new(mark("triangle"))),
            MarkPart::Gap(MarkGap {
                gap: MarkSize {
                    points: -10.0,
                    ratio: 0.0,
                },
            }),
        ]);
        let mut request = batch(spec.clone(), line(Point::ZERO, Point::new(100.0, 0.0)));
        request.placements[0].station = MarkStation::Ratio { value: 0.5 };
        let result = request.geometry().unwrap().marks;
        close(result[0].end, -5.0);
        close(result[0].back.x - result[0].tip.x, 5.0);
        let bounds = result[0].paths[0].outline.path.bounding_box();
        close((bounds.x0 + bounds.x1) / 2.0, 50.0);
        spec.parts.as_mut().unwrap().push(MarkPart::Gap(MarkGap {
            gap: MarkSize {
                points: f64::INFINITY,
                ratio: 0.0,
            },
        }));
        assert!(spec.prepare(CONTEXT).is_err());
    }

    #[test]
    fn integer_cbor_numeric_fields_and_recursive_parts_roundtrip() {
        let size = Value::Map(vec![
            (Value::Text("points".into()), Value::Integer(3.into())),
            (Value::Text("ratio".into()), Value::Integer(4.into())),
        ]);
        let value = Value::Map(vec![
            (Value::Text("shape".into()), Value::Text("triangle".into())),
            (Value::Text("length".into()), size),
            (Value::Text("shorten".into()), Value::Integer(1.into())),
        ]);
        let mut parsed: MarkSpec = ciborium::de::from_reader(wire(&value).as_slice()).unwrap();
        close(parsed.prepare(CONTEXT).unwrap().end, 14.0);
        parsed.shorten = None;
        let mut combine = mark("combine");
        combine.parts = Some(vec![
            MarkPart::Mark(Box::new(parsed)),
            MarkPart::Gap(MarkGap {
                gap: MarkSize {
                    points: -1.0,
                    ratio: 0.5,
                },
            }),
        ]);
        let request = batch(combine, line(Point::ZERO, Point::new(100.0, 0.0)));
        assert_eq!(
            mark_geometry_bytes(&wire(&request)).unwrap(),
            wire(&request.geometry().unwrap())
        );
    }

    #[test]
    fn shaft_outline_and_footprint_honor_butt_round_and_square_caps() {
        // Hide the head to isolate a very short shaft's cap geometry.
        let mut spec = mark("triangle");
        spec.fill = Some(false);
        spec.stroke = Some(false);
        spec.shorten = Some(0.0);
        let mut request = batch(spec, line(Point::ZERO, Point::new(0.5, 0.0)));
        for (cap, x0, x1) in [
            (Cap::Butt, 0.0, 0.5),
            (Cap::Round, -1.0, 1.5),
            (Cap::Square, -1.0, 1.5),
        ] {
            request.templates[0].context.shaft_stroke.cap = cap;
            let batch_result = request.clone().geometry().unwrap();
            let result = batch_result.marks[0].clone();
            assert_eq!(result.shaft_style.cap, cap);
            let bounds = result.shaft_outline.path.bounding_box();
            close(bounds.x0, x0);
            close(bounds.x1, x1);
            assert_eq!(result.footprint.path, result.shaft_outline.path);
            let authoritative = kurbo::stroke(
                result.shaft.path.elements().iter().copied(),
                &result.shaft_style.stroke(CONTEXT.line_thickness),
                &StrokeOpts::default(),
                1e-6,
            );
            assert_eq!(result.shaft_outline.path, authoritative);
            assert_eq!(
                mark_geometry_bytes(&wire(&request)).unwrap(),
                wire(&batch_result)
            );
        }
    }

    #[test]
    fn shaft_style_defaults_and_strict_wire_parameters_are_normalized() {
        let value = Value::Map(vec![
            (Value::Text("units-per-pt".into()), Value::Integer(2.into())),
            (
                Value::Text("line-thickness".into()),
                Value::Integer(2.into()),
            ),
        ]);
        let context: MarkContext = ciborium::de::from_reader(wire(&value).as_slice()).unwrap();
        assert_eq!(context.shaft_stroke, MarkStrokeStyle::default());
        for (field, value) in [
            ("cap", Value::Text("invalid".into())),
            ("join", Value::Text("invalid".into())),
            ("paint", Value::Bool(true)),
        ] {
            let style = Value::Map(vec![(Value::Text(field.into()), value)]);
            assert!(
                ciborium::de::from_reader::<MarkStrokeStyle, _>(wire(&style).as_slice()).is_err()
            );
        }
        let mut request = batch(mark("triangle"), curve());
        request.templates[0].context.shaft_stroke.miter_limit = f64::INFINITY;
        assert!(request.geometry().is_err());
    }

    #[test]
    fn mirrored_fill_and_stroke_regions_cannot_cancel_painted_coverage() {
        for stroke in [false, true] {
            let mut spec = mark("triangle");
            spec.rev = Some(true);
            spec.shorten = Some(0.0);
            spec.stroke = Some(stroke);
            let result = geometry(spec, line(Point::ZERO, Point::new(100.0, 0.0)));
            assert!(result.footprint.path.winding(Point::new(95.0, 0.5)) > 0);
            let head = &result.paths[0];
            let head_stroke = kurbo::stroke(
                head.path.path.elements().iter().copied(),
                &Stroke::new(CONTEXT.line_thickness)
                    .with_join(Join::Miter)
                    .with_miter_limit(head.miter_limit),
                &StrokeOpts::default(),
                1e-6,
            );
            for x in 0..38 {
                for y in 0..32 {
                    let point =
                        Point::new(83.0 + x as f64 * 0.5 + 0.13, -8.0 + y as f64 * 0.5 + 0.17);
                    let painted = (head.fill && head.path.path.winding(point) != 0)
                        || (head.stroke && head_stroke.winding(point) != 0)
                        || result.shaft_outline.path.winding(point) != 0;
                    assert_eq!(
                        result.footprint.path.winding(point) != 0,
                        painted,
                        "{point:?}, stroke={stroke}"
                    );
                }
            }
        }
        // Normalizing fill contours must not turn a stroke's inner hole into
        // another filled region. The shaft crossing its hole still paints.
        let mut spec = mark("circle");
        spec.fill = Some(false);
        spec.stroke = Some(true);
        let mut request = batch(spec, line(Point::ZERO, Point::new(100.0, 0.0)));
        request.placements[0].station = MarkStation::Ratio { value: 0.5 };
        let result = request.geometry().unwrap().marks.remove(0);
        assert_eq!(result.footprint.path.winding(Point::new(50.0, 1.5)), 0);
        assert!(result.footprint.path.winding(Point::new(50.0, 0.5)) > 0);
        assert!(result.footprint.path.winding(Point::new(50.0, 3.0)) > 0);
    }

    #[test]
    fn chord_contact_refines_the_near_tangent_polyline_and_circle() {
        let mut path = BezPath::new();
        path.move_to((-19.999, 10.0));
        path.line_to((-19.999, 0.0));
        path.line_to((0.0, 0.0));
        let mut spec = mark("triangle");
        spec.length = Some(MarkSize {
            points: 10.0,
            ratio: 0.0,
        });
        let result = geometry(spec.clone(), path);
        let expected = Point::new(-19.999, (400.0 - 19.999_f64 * 19.999).sqrt());
        point_close(result.back.into(), expected);
        point_close(result.shaft_contact.into(), expected);
        assert!((Point::from(result.tip).distance(result.back.into()) - 20.0).abs() < 1e-6);
        let circle = kurbo::Circle::new(Point::ZERO, 10.0).to_path(1e-6);
        let carrier = PreparedCarrier::new(circle).unwrap();
        let origin = carrier.point(carrier.length).unwrap();
        let (contact, at) = carrier.chord_contact(carrier.length, 20.0, 1.0).unwrap();
        assert!((origin.distance(contact) - 20.0).abs() < 1e-6);
        point_close(contact, Point::new(-10.0, 0.0));
        assert!(at > 0.0 && at < carrier.length);
    }

    #[test]
    fn disconnected_carriers_and_closed_endpoint_bend_are_rejected_explicitly() {
        let mut path = line(Point::ZERO, Point::new(10.0, 0.0));
        path.move_to((100.0, 0.0));
        path.line_to((110.0, 0.0));
        let request = batch(mark("triangle"), path);
        assert!(
            request
                .clone()
                .geometry()
                .unwrap_err()
                .contains("connected")
        );
        assert!(
            mark_geometry_bytes(&wire(&request))
                .unwrap_err()
                .contains("connected")
        );
        let mut closed = line(Point::ZERO, Point::new(30.0, 0.0));
        closed.line_to((30.0, 30.0));
        closed.close_path();
        let mut spec = mark("triangle");
        spec.fit = Some(MarkFit::Bend);
        let mut request = batch(spec, closed);
        for station in [MarkStation::Start, MarkStation::End] {
            request.placements[0].station = station;
            assert!(request.clone().geometry().unwrap_err().contains("closed"));
        }
        request.placements[0].station = MarkStation::Ratio { value: 0.5 };
        assert!(request.geometry().is_ok());
    }

    #[test]
    fn shifted_endpoint_bend_clamps_to_available_inward_station() {
        let mut spec = mark("triangle");
        spec.fit = Some(MarkFit::Bend);
        let mut request = batch(spec, line(Point::ZERO, Point::new(100.0, 0.0)));
        request.placements[0].shift = -95.0;
        let end = request.clone().geometry().unwrap().marks.remove(0);
        close(end.tip.x, 5.0);
        close(end.back.x, -10.0);
        close(end.shaft_contact.x, 0.0);
        request.placements[0].station = MarkStation::Start;
        request.placements[0].direction = MarkDirection::Backward;
        request.placements[0].shift = 95.0;
        let start = request.geometry().unwrap().marks.remove(0);
        close(start.tip.x, 95.0);
        close(start.back.x, 110.0);
        close(start.shaft_contact.x, 100.0);
    }

    #[test]
    fn composites_require_nonempty_parts_at_the_shared_boundary() {
        let mut spec = mark("bar");
        spec.shape = "combine".into();
        assert!(spec.prepare(CONTEXT).err().unwrap().contains("nonempty"));
        spec.parts = Some(Vec::new());
        assert!(spec.prepare(CONTEXT).err().unwrap().contains("nonempty"));
    }

    #[test]
    fn composite_children_reject_explicit_fit_and_shorten_even_defaults() {
        for (fit, shorten) in [
            (Some(MarkFit::Chord), None),
            (Some(MarkFit::Bend), None),
            (None, Some(1.0)),
            (None, Some(0.0)),
        ] {
            let mut child = mark("triangle");
            child.fit = fit;
            child.shorten = shorten;
            let mut spec = mark("combine");
            spec.parts = Some(vec![MarkPart::Mark(Box::new(child))]);
            assert!(spec.prepare(CONTEXT).is_err());
            let err = ciborium::de::from_reader::<MarkSpec, _>(wire(&spec).as_slice())
                .unwrap_err()
                .to_string();
            assert!(err.contains("children"), "{err}");
        }
        let mut spec = mark("combine");
        spec.parts = Some(vec![
            MarkPart::Mark(Box::new(mark("triangle"))),
            MarkPart::Mark(Box::new(mark("stealth"))),
        ]);
        spec.fit = Some(MarkFit::Bend);
        let result = geometry(spec, line(Point::ZERO, Point::new(100.0, 0.0)));
        close(result.end, 24.0);
        close(result.shaft_contact.x, 76.0);
    }

    #[test]
    fn selected_heads_share_one_shaft_while_candidates_remain_independent() {
        for fit in [MarkFit::Chord, MarkFit::Bend] {
            let mut spec = mark("triangle");
            spec.fit = Some(fit);
            let mut request = batch(spec, line(Point::ZERO, Point::new(100.0, 0.0)));
            let mut start = request.placements[0].clone();
            start.station = MarkStation::Start;
            start.direction = MarkDirection::Backward;
            request.placements.push(start);
            let candidates = request.clone().geometry().unwrap();
            assert!(candidates.shafts.is_empty());
            assert_eq!(
                candidates.marks[0].shaft.path.elements()[1],
                PathEl::LineTo(Point::new(85.0, 0.0))
            );
            assert_eq!(
                candidates.marks[1].shaft.path.elements()[0],
                PathEl::MoveTo(Point::new(15.0, 0.0))
            );
            request.mode = MarkGeometryMode::Selected;
            let selected = request.clone().geometry().unwrap();
            assert_eq!(selected.shafts.len(), 1);
            assert_eq!(
                selected.shafts[0].shaft.path.elements(),
                &[
                    PathEl::MoveTo(Point::new(15.0, 0.0)),
                    PathEl::LineTo(Point::new(85.0, 0.0)),
                ]
            );
            for head in &selected.marks {
                assert!(head.shaft.path.is_empty() && head.shaft_outline.path.is_empty());
                assert_eq!(
                    wire(&head.paths),
                    wire(
                        &candidates
                            .marks
                            .iter()
                            .find(|mark| mark.tip == head.tip)
                            .unwrap()
                            .paths
                    )
                );
            }
            assert_eq!(
                selected.shafts[0].footprint,
                selected.shafts[0].shaft_outline
            );
            assert!(
                selected.shafts[0]
                    .footprint
                    .path
                    .winding(Point::new(50.0, 0.5))
                    > 0
            );
            assert_eq!(
                selected.shafts[0]
                    .footprint
                    .path
                    .winding(Point::new(5.0, 0.5)),
                0
            );
            assert_eq!(
                mark_geometry_bytes(&wire(&request)).unwrap(),
                wire(&selected)
            );
            request.placements.reverse();
            assert_eq!(
                wire(&request.geometry().unwrap().shafts),
                wire(&selected.shafts)
            );
        }
    }

    #[test]
    fn selected_endpoint_conflicts_and_mismatched_shaft_styles_are_rejected() {
        let mut request = batch(mark("triangle"), line(Point::ZERO, Point::new(100.0, 0.0)));
        request.mode = MarkGeometryMode::Selected;
        request.placements.push(request.placements[0].clone());
        assert!(
            request
                .clone()
                .geometry()
                .unwrap_err()
                .contains("multiple selected")
        );
        request.placements[1].station = MarkStation::Start;
        request.placements[1].direction = MarkDirection::Backward;
        request.templates.push(request.templates[0].clone());
        request.placements[1].template = 1;
        request.templates[1].context.shaft_stroke.cap = Cap::Round;
        assert!(
            request
                .clone()
                .geometry()
                .unwrap_err()
                .contains("conflicting shaft")
        );
        request.templates[1].context.shaft_stroke = MarkStrokeStyle::default();
        request.templates[1].context.line_thickness = 4.0;
        assert!(
            request
                .geometry()
                .unwrap_err()
                .contains("conflicting shaft")
        );
    }
}
