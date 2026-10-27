//! Data-only mark geometry, freshly implemented from tiptoe 0.4's MIT geometry.
//! Attribution: Mc-Zen, copyright 2024–2025, MIT license.
//! All geometry coordinates are drawing units; angles are radians. Paint stays
//! with the caller. A combined mark is prepared locally and fitted only once.
//!
//! The reproduced geometry is covered by this permission notice:
//! Permission is hereby granted, free of charge, to any person obtaining a copy
//! of this software and associated documentation files (the "Software"), to deal
//! in the Software without restriction, including without limitation the rights
//! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//! copies of the Software, and to permit persons to whom the Software is
//! furnished to do so, subject to the following conditions:
//! The above copyright notice and this permission notice shall be included in
//! all copies or substantial portions of the Software.
//! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//! THE SOFTWARE.

use crate::curve_api::{
    CurvePathOutput, CurvePoint, PathTrimmer, deserialize_f64, parse_stroke_cap, parse_stroke_join,
};
use kurbo::{
    Affine, BezPath, Cap, Join, ParamCurve, PathEl, Point, Shape, Stroke, StrokeOpts, Vec2,
};
use serde::{Deserialize, Serialize};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

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

/// `back` is the full-size mark's geometric rear/notch contact. At shorten=1
/// it coincides with `shaft-contact` on an ordinary carrier. If a carrier is too
/// short to fit the requested head, only the shaft contact is clamped inward;
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

    fn validate_fields(&self) -> Result<(), String> {
        let shape = self.shape.as_str();
        let allowed = match shape {
            "triangle" | "stealth" | "round" => "length width inset rev fill stroke",
            "straight" => "length width rev stroke",
            "tikz" => "width stroke",
            "barb" | "hooks" => "width arc rev stroke",
            "bar" => "width align stroke",
            "bracket" => "length width rev stroke",
            "circle" | "square" | "diamond" => "length width align fill stroke",
            "rays" => "length n phase align stroke",
            "combine" => "parts",
            _ => return Err(format!("unknown mark shape: {shape}")),
        };
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
            if present && !allowed.split_whitespace().any(|field| field == name) {
                return Err(if shape == "tikz" && name == "arc" {
                    "tikz.arc is unused in tiptoe 0.4 and is not supported".into()
                } else {
                    format!("{shape} does not support {name}")
                });
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
        if shape == "combine" && self.parts.as_deref().unwrap_or_default().is_empty() {
            return Err("combine.parts must be nonempty".into());
        }
        if self.parts.as_deref().unwrap_or_default().iter().any(|part| {
            matches!(part,MarkPart::Mark(child) if child.fit.is_some() || child.shorten.is_some())
        }) {
            return Err("composite children must omit fit and shorten, including their explicit defaults; set them on the composite".into());
        }
        Ok(())
    }

    /// Prepare dimensions and painted local outlines once. `inset` is restricted
    /// to [0,1]; round's full-inset case retains tiptoe's L-thickness rule.
    pub fn prepare(&self, context: MarkContext) -> Result<PreparedMark, String> {
        context.validate()?;
        self.prepare_nested(context, 0)
    }

    fn prepare_nested(&self, context: MarkContext, depth: usize) -> Result<PreparedMark, String> {
        if depth > 64 {
            return Err("mark nesting exceeds 64 levels".into());
        }
        self.validate_fields()?;
        let shape = self.shape.as_str();
        if shape == "combine" {
            let mut paths = Vec::new();
            let mut cursor = 0.0;
            for part in self.parts.as_deref().unwrap_or_default() {
                match part {
                    MarkPart::Gap(gap) => cursor += gap.gap.resolve(context)?,
                    MarkPart::Mark(spec) => {
                        let child = spec.prepare_nested(context, depth + 1)?;
                        for mut drawable in child.paths {
                            drawable
                                .path
                                .path
                                .apply_affine(Affine::translate((-cursor, 0.0)));
                            drawable
                                .outline
                                .path
                                .apply_affine(Affine::translate((-cursor, 0.0)));
                            paths.push(drawable);
                        }
                        cursor += child.end;
                    }
                }
                if !cursor.is_finite() {
                    return Err("composite contact must be finite".into());
                }
            }
            return Ok(PreparedMark {
                longitudinal_center: PreparedMark::longitudinal_center(&paths),
                paths,
                end: cursor,
                fit: self.fit.unwrap_or_default(),
                shorten: self.shorten.unwrap_or(1.0),
                line_thickness: context.line_thickness,
                shaft_style: context.shaft_stroke,
            });
        }
        let closed = matches!(
            shape,
            "triangle" | "stealth" | "round" | "circle" | "square" | "diamond"
        );
        let fill = self.fill.unwrap_or(closed);
        // The common filled triangle/stealth branch has no boundary stroke.
        // Explicit normalized flags select the general stroked geometry.
        let optimized =
            matches!(shape, "triangle" | "stealth") && self.fill.is_none() && self.stroke.is_none();
        let stroke = self.stroke.unwrap_or(!optimized);
        let t = if stroke { context.line_thickness } else { 0.0 };
        let s = t / 2.0;
        let default_length = match shape {
            "circle" | "square" => MarkSize {
                points: 0.0,
                ratio: 4.0,
            },
            "diamond" => MarkSize {
                points: 0.0,
                ratio: 5.6569,
            },
            "rays" => MarkSize {
                points: 0.0,
                ratio: 2.8,
            },
            _ => MarkSize {
                points: 3.0,
                ratio: 4.5,
            },
        };
        let mut l = self.length.unwrap_or(default_length).resolve(context)?;
        let default_width = match shape {
            "bar" | "bracket" => MarkSize {
                points: 2.4,
                ratio: 3.6,
            }
            .resolve(context)?,
            "tikz" | "barb" | "hooks" => default_length.resolve(context)?,
            "circle" | "square" | "diamond" => l,
            _ => 0.8 * l,
        };
        let w = self
            .width
            .map(|size| size.resolve(context))
            .transpose()?
            .unwrap_or(default_width);
        if shape == "bracket" && self.length.is_none() {
            l = 0.3 * w;
        }
        if l < 0.0 || w < 0.0 {
            return Err("mark dimensions must be nonnegative".into());
        }
        let rev = self.rev.unwrap_or(false);
        let aligned = self.align.unwrap_or_default() == MarkAlign::End;
        let join = if matches!(shape, "round" | "tikz") {
            Join::Round
        } else {
            Join::Miter
        };
        let cap = if matches!(shape, "straight" | "tikz") {
            Cap::Round
        } else {
            Cap::Butt
        };
        let mut path = BezPath::new();
        let mut end;
        let mut reverse_offset = l;
        let mut polygon = |points: &[(f64, f64)], close: bool| {
            if let Some(first) = points.first() {
                path.move_to(*first);
                for point in &points[1..] {
                    path.line_to(*point);
                }
                if close {
                    path.close_path();
                }
            }
        };
        match shape {
            "triangle" | "stealth" | "round" => {
                let inset = self
                    .inset
                    .unwrap_or(if shape == "triangle" { 0.0 } else { 0.4 });
                let d = if inset == 1.0 && shape == "round" {
                    (l - t).max(0.0)
                } else {
                    l * inset
                };
                if shape == "round" {
                    polygon(
                        &[
                            (-s, 0.0),
                            (-l + s, (w / 2.0 - s).max(0.0)),
                            (-l + d + s, 0.0),
                            (-l + s, -(w / 2.0 - s).max(0.0)),
                        ],
                        true,
                    );
                    end = if rev { l - s } else { l - d - s };
                } else if s > 0.0 && l > 0.0 && w > 0.0 && inset < 1.0 {
                    let tan_a = w / 2.0 / l;
                    let x3 = s / tan_a.atan().sin();
                    let (x, x4, y) = if inset == 0.0 {
                        let x = l - s;
                        (x, x, tan_a * (x - x3))
                    } else {
                        let tan_b = w / 2.0 / d;
                        let x1 = s / tan_b.atan().sin();
                        let x2 = l - d - x1 - x3;
                        let x = -x2 * tan_a / (tan_b - tan_a);
                        (l - d - x1 - x, l - d - x1, tan_b * x)
                    };
                    polygon(&[(-x3, 0.0), (-x, y), (-x4, 0.0), (-x, -y)], true);
                    end = if rev { l - x3 } else { l - d };
                } else {
                    polygon(
                        &[(0.0, 0.0), (-l, w / 2.0), (-l + d, 0.0), (-l, -w / 2.0)],
                        true,
                    );
                    end = if rev { l } else { l - d };
                }
            }
            "straight" => {
                let q = if l > 0.0 && w > 0.0 {
                    s / (w / 2.0).atan2(l).sin()
                } else {
                    0.0
                };
                polygon(
                    &[
                        (-l + s, (w / 2.0 - s).max(0.0)),
                        (-q, 0.0),
                        (-l + s, -(w / 2.0 - s).max(0.0)),
                    ],
                    false,
                );
                end = if rev { l - q } else { q };
            }
            "tikz" => {
                path.move_to((-0.42 * w, -w / 2.0 + s));
                path.curve_to((-0.32 * w, -0.1 * w + s), (-s, 0.0), (-s, 0.0));
                path.curve_to(
                    (-s, 0.0),
                    (-0.32 * w, 0.1 * w - s),
                    (-0.42 * w, w / 2.0 - s),
                );
                end = t;
            }
            "barb" | "hooks" => {
                let arc = self.arc.unwrap_or(PI);
                let radius = (w / 2.0 - s).max(0.0);
                if shape == "barb" {
                    PreparedMark::arc(
                        &mut path,
                        Point::new(-w / 2.0, 0.0),
                        radius,
                        -arc / 2.0,
                        arc,
                    );
                    end = if rev { radius } else { s };
                    reverse_offset = w / 2.0;
                } else {
                    let r = radius / 2.0;
                    PreparedMark::arc(&mut path, Point::new(-r - s, r), r, -FRAC_PI_2, arc);
                    PreparedMark::arc(&mut path, Point::new(-r - s, -r), r, FRAC_PI_2, -arc);
                    end = if rev { s } else { r };
                    reverse_offset = r + s;
                }
            }
            "bar" => {
                let x = if aligned { -s } else { 0.0 };
                polygon(&[(x, -w / 2.0), (x, w / 2.0)], false);
                end = if aligned { s } else { 0.0 };
            }
            "bracket" => {
                polygon(
                    &[
                        (-l, -(w / 2.0 - s).max(0.0)),
                        (-s, -(w / 2.0 - s).max(0.0)),
                        (-s, (w / 2.0 - s).max(0.0)),
                        (-l, (w / 2.0 - s).max(0.0)),
                    ],
                    false,
                );
                end = if rev { l } else { s };
            }
            "circle" | "square" => {
                let offset = if aligned { 0.0 } else { l / 2.0 };
                let left = -l + s + offset;
                let right = -s + offset;
                let half = (w / 2.0 - s).max(0.0);
                if shape == "circle" {
                    path = kurbo::Ellipse::new(
                        Point::new((left + right) / 2.0, 0.0),
                        ((right - left).max(0.0) / 2.0, half),
                        0.0,
                    )
                    .to_path(1e-6);
                } else {
                    polygon(
                        &[(left, -half), (right, -half), (right, half), (left, half)],
                        true,
                    );
                }
                end = l - s - offset;
            }
            "diamond" => {
                let angle = w.atan2(l);
                let q = if angle.sin().abs() > 1e-12 {
                    s / angle.sin()
                } else {
                    0.0
                };
                let u = if angle.cos().abs() > 1e-12 {
                    s / angle.cos()
                } else {
                    0.0
                };
                polygon(
                    &[
                        (-l / 2.0, w / 2.0 - u),
                        (-q, 0.0),
                        (-l / 2.0, -w / 2.0 + u),
                        (-l + q, 0.0),
                    ],
                    true,
                );
                let offset = if aligned { 0.0 } else { l / 2.0 };
                path.apply_affine(Affine::translate((offset, 0.0)));
                end = l - q - offset;
            }
            "rays" => {
                let n = self.n.unwrap_or(4);
                let phase = self
                    .phase
                    .unwrap_or(if n == 4 { PI / 4.0 } else { -FRAC_PI_2 });
                for i in 0..n {
                    let a = phase + TAU * i as f64 / n as f64;
                    let x = if aligned { -l } else { 0.0 };
                    path.move_to((x, 0.0));
                    path.line_to((x + l * a.cos(), l * a.sin()));
                }
                end = if aligned { l } else { 0.0 };
            }
            _ => unreachable!(),
        }
        if rev {
            path.apply_affine(Affine::new([-1.0, 0.0, 0.0, 1.0, -reverse_offset, 0.0]));
        }
        if !path.is_finite() || !end.is_finite() {
            return Err("mark dimensions produced nonfinite geometry".into());
        }
        // Stroke-dominated degenerate dimensions cannot require a negative
        // shaft retraction. Signed composite gaps remain signed independently.
        end = end.max(0.0);
        let mut outline = BezPath::new();
        if fill {
            // Every filled primitive is one region. Keep its original drawable
            // commands, but orient its coverage positively so a mirrored head
            // cannot cancel a positive stroke region where the paints overlap.
            let region = if path.area() < 0.0 {
                path.reverse_subpaths()
            } else {
                path.clone()
            };
            outline.extend(region.elements().iter().copied());
        }
        if stroke && t > 0.0 {
            outline.extend(
                kurbo::stroke(
                    path.elements().iter().copied(),
                    &Stroke::new(t)
                        .with_join(join)
                        .with_caps(cap)
                        .with_miter_limit(7.0),
                    &StrokeOpts::default(),
                    1e-6,
                )
                .elements()
                .iter()
                .copied(),
            );
        }
        let outline = CurvePathOutput { path: outline };
        let paths = vec![MarkDrawable {
            path: CurvePathOutput { path },
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
            miter_limit: 7.0,
        }];
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

    fn arc(path: &mut BezPath, center: Point, radius: f64, start: f64, sweep: f64) {
        let count = (sweep.abs() / FRAC_PI_2).ceil().max(1.0) as usize;
        let delta = sweep / count as f64;
        let at = |a: f64| center + Vec2::new(a.cos(), a.sin()) * radius;
        path.move_to(at(start));
        for i in 0..count {
            let a = start + i as f64 * delta;
            let b = a + delta;
            let k = 4.0 / 3.0 * (delta / 4.0).tan() * radius;
            path.curve_to(
                at(a) + Vec2::new(-a.sin(), a.cos()) * k,
                at(b) - Vec2::new(-b.sin(), b.cos()) * k,
                at(b),
            );
        }
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
    fn every_shape_default_matches_contact_and_geometry_rules() {
        let l = 15.0;
        let q = 1.0 / (0.4_f64.atan().sin());
        let contacts = [
            l,
            q,
            0.6 * l,
            0.6 * l - 1.0,
            2.0,
            1.0,
            (l / 2.0 - 1.0) / 2.0,
            0.0,
            1.0,
            3.0,
            3.0,
            5.6569 - 2.0_f64.sqrt(),
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
    fn explicit_triangle_stroke_accounts_for_tip_miter_and_reverse() {
        let mut spec = mark("triangle");
        spec.stroke = Some(true);
        let normal = spec.prepare(CONTEXT).unwrap();
        let q = 1.0 / 0.4_f64.atan().sin();
        assert_eq!(
            normal.paths[0].path.path.elements()[0],
            PathEl::MoveTo(Point::new(-q, 0.0))
        );
        close(normal.end, 15.0);
        assert_eq!(normal.paths[0].join, "miter");
        spec.rev = Some(true);
        close(spec.prepare(CONTEXT).unwrap().end, 15.0 - q);
    }

    #[test]
    fn tikz_has_degenerate_tip_handles_and_rejects_unused_arc() {
        let prepared = mark("tikz").prepare(CONTEXT).unwrap();
        let tip = Point::new(-1.0, 0.0);
        assert_eq!(
            prepared.paths[0].path.path.elements()[1],
            PathEl::CurveTo(Point::new(-4.8, -0.5), tip, tip)
        );
        assert_eq!(
            prepared.paths[0].path.path.elements()[2],
            PathEl::CurveTo(tip, Point::new(-4.8, 0.5), Point::new(-6.3, 6.5))
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
        for shape in ["bar", "circle", "square", "diamond", "rays"] {
            let normal = mark(shape).prepare(CONTEXT).unwrap();
            let mut spec = mark(shape);
            spec.align = Some(MarkAlign::End);
            let aligned = spec.prepare(CONTEXT).unwrap();
            assert!(aligned.end > normal.end, "{shape}");
        }
        let prepared = mark("rays").prepare(CONTEXT).unwrap();
        if let PathEl::LineTo(point) = prepared.paths[0].path.path.elements()[1] {
            close(point.x, 5.6 / 2.0_f64.sqrt());
            close(point.y, 5.6 / 2.0_f64.sqrt());
        } else {
            panic!("ray must be a line");
        }
        let mut spec = mark("rays");
        spec.n = Some(3);
        let prepared = spec.prepare(CONTEXT).unwrap();
        if let PathEl::LineTo(point) = prepared.paths[0].path.path.elements()[1] {
            close(point.y, -5.6);
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
                    .with_miter_limit(7.0),
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
