use kurbo::{
    BezPath, CubicBez, Line, ParamCurve, ParamCurveArclen, ParamCurveDeriv, PathEl, PathSeg, Point,
    Vec2, fit_to_bezpath, fit_to_bezpath_opt, offset::CubicOffset,
};
use serde::{
    Deserialize, Deserializer, Serialize, Serializer,
    de::{self, Visitor},
};
use std::fmt;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct CurvePoint {
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct CubicBezierSpec {
    pub start: CurvePoint,
    pub end: CurvePoint,
    pub ctrl_a: CurvePoint,
    pub ctrl_b: CurvePoint,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct CubicPathSpec {
    #[serde(flatten)]
    pub curve: CubicBezierSpec,
    #[serde(
        default = "default_arclen_accuracy",
        deserialize_with = "deserialize_f64"
    )]
    pub accuracy: f64,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct TrimPathSpec {
    #[serde(deserialize_with = "deserialize_bez_path")]
    pub path: BezPath,
    #[serde(default, deserialize_with = "deserialize_f64")]
    pub start_outset: f64,
    #[serde(default, deserialize_with = "deserialize_f64")]
    pub end_outset: f64,
    #[serde(
        default = "default_arclen_accuracy",
        deserialize_with = "deserialize_f64"
    )]
    pub accuracy: f64,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct PatternPathSpec {
    #[serde(deserialize_with = "deserialize_bez_path")]
    pub path: BezPath,
    #[serde(default = "default_pattern")]
    pub pattern: PatternInput,
    #[serde(
        default = "default_pattern_amplitude",
        deserialize_with = "deserialize_f64"
    )]
    pub amplitude: f64,
    #[serde(
        default = "default_pattern_wavelength",
        deserialize_with = "deserialize_f64"
    )]
    pub wavelength: f64,
    #[serde(default, deserialize_with = "deserialize_f64")]
    pub phase: f64,
    #[serde(default = "default_samples_per_period")]
    pub samples_per_period: usize,
    #[serde(
        default = "default_coil_longitudinal_scale",
        deserialize_with = "deserialize_f64"
    )]
    pub coil_longitudinal_scale: f64,
    #[serde(default = "default_anchor_endpoint")]
    pub anchor_start: bool,
    #[serde(default = "default_anchor_endpoint")]
    pub anchor_end: bool,
    #[serde(
        default = "default_arclen_accuracy",
        deserialize_with = "deserialize_f64"
    )]
    pub accuracy: f64,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct ParallelPathSpec {
    #[serde(deserialize_with = "deserialize_bez_path")]
    pub path: BezPath,
    #[serde(default, deserialize_with = "deserialize_f64")]
    pub distance: f64,
    #[serde(default, deserialize_with = "deserialize_f64")]
    pub start_outset: f64,
    #[serde(default, deserialize_with = "deserialize_f64")]
    pub end_outset: f64,
    #[serde(
        default = "default_arclen_accuracy",
        deserialize_with = "deserialize_f64"
    )]
    pub accuracy: f64,
    #[serde(default = "default_parallel_optimize")]
    pub optimize: bool,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
#[serde(untagged)]
pub enum PatternInput {
    Name(String),
    Points(PointPatternInput),
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct PointPatternInput {
    #[serde(default = "default_points_pattern_kind")]
    pub kind: String,
    #[serde(default)]
    pub name: Option<String>,
    pub points: Vec<PatternPointInput>,
    #[serde(default = "default_pattern_interpolation")]
    pub interpolation: String,
    #[serde(default)]
    pub endpoint_ramp: bool,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct PatternPointInput {
    #[serde(
        default,
        alias = "t",
        alias = "phase",
        deserialize_with = "deserialize_optional_f64"
    )]
    pub at: Option<f64>,
    #[serde(default, deserialize_with = "deserialize_f64")]
    pub x: f64,
    #[serde(default, deserialize_with = "deserialize_f64")]
    pub y: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct HobbyThroughSpec {
    pub start: CurvePoint,
    pub through: CurvePoint,
    pub end: CurvePoint,
    #[serde(default = "default_hobby_omega")]
    pub omega: f64,
    #[serde(
        default = "default_arclen_accuracy",
        deserialize_with = "deserialize_f64"
    )]
    pub accuracy: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HobbySplineSpec {
    pub points: Vec<CurvePoint>,
    #[serde(default = "default_hobby_omega")]
    pub omega: f64,
    #[serde(
        default = "default_arclen_accuracy",
        deserialize_with = "deserialize_f64"
    )]
    pub accuracy: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct LineSegmentSpec {
    pub start: CurvePoint,
    pub end: CurvePoint,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PatternPathOutput {
    #[serde(
        serialize_with = "serialize_bez_path",
        deserialize_with = "deserialize_bez_path"
    )]
    pub path: BezPath,
    pub pattern: String,
    pub length: f64,
    pub points: Vec<CurvePoint>,
    pub curves: Vec<CubicBezierSpec>,
    pub segments: Vec<LineSegmentSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CurvePathOutput {
    #[serde(
        serialize_with = "serialize_bez_path",
        deserialize_with = "deserialize_bez_path"
    )]
    pub path: BezPath,
    pub length: f64,
    pub points: Vec<CurvePoint>,
    pub curves: Vec<CubicBezierSpec>,
    pub segments: Vec<LineSegmentSpec>,
}

fn default_hobby_omega() -> f64 {
    1.0
}

fn default_arclen_accuracy() -> f64 {
    1e-3
}

fn default_pattern() -> PatternInput {
    PatternInput::Name("wave".to_string())
}

fn default_points_pattern_kind() -> String {
    "points".to_string()
}

fn default_pattern_interpolation() -> String {
    "smooth".to_string()
}

fn default_pattern_amplitude() -> f64 {
    0.1
}

fn default_pattern_wavelength() -> f64 {
    1.0
}

fn default_samples_per_period() -> usize {
    16
}

fn default_coil_longitudinal_scale() -> f64 {
    1.25
}

fn default_anchor_endpoint() -> bool {
    true
}

fn default_parallel_optimize() -> bool {
    true
}

fn encode_cbor<T: Serialize>(value: &T) -> Result<Vec<u8>, String> {
    let mut bytes = Vec::new();
    ciborium::ser::into_writer(value, &mut bytes)
        .map_err(|err| format!("Failed to serialize CBOR value: {err}"))?;
    Ok(bytes)
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
struct WirePath {
    elements: Vec<WirePathElement>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "kebab-case")]
enum WirePathElement {
    MoveTo {
        point: CurvePoint,
    },
    LineTo {
        point: CurvePoint,
    },
    QuadTo {
        ctrl: CurvePoint,
        end: CurvePoint,
    },
    CurveTo {
        ctrl_a: CurvePoint,
        ctrl_b: CurvePoint,
        end: CurvePoint,
    },
    ClosePath,
}

impl From<&BezPath> for WirePath {
    fn from(path: &BezPath) -> Self {
        let elements = path
            .elements()
            .iter()
            .map(|element| match *element {
                PathEl::MoveTo(point) => WirePathElement::MoveTo {
                    point: point.into(),
                },
                PathEl::LineTo(point) => WirePathElement::LineTo {
                    point: point.into(),
                },
                PathEl::QuadTo(ctrl, end) => WirePathElement::QuadTo {
                    ctrl: ctrl.into(),
                    end: end.into(),
                },
                PathEl::CurveTo(ctrl_a, ctrl_b, end) => WirePathElement::CurveTo {
                    ctrl_a: ctrl_a.into(),
                    ctrl_b: ctrl_b.into(),
                    end: end.into(),
                },
                PathEl::ClosePath => WirePathElement::ClosePath,
            })
            .collect();
        Self { elements }
    }
}

impl From<WirePath> for BezPath {
    fn from(value: WirePath) -> Self {
        let mut path = BezPath::new();
        for element in value.elements {
            match element {
                WirePathElement::MoveTo { point } => path.move_to(Point::from(point)),
                WirePathElement::LineTo { point } => path.line_to(Point::from(point)),
                WirePathElement::QuadTo { ctrl, end } => {
                    path.quad_to(Point::from(ctrl), Point::from(end));
                }
                WirePathElement::CurveTo {
                    ctrl_a,
                    ctrl_b,
                    end,
                } => path.curve_to(Point::from(ctrl_a), Point::from(ctrl_b), Point::from(end)),
                WirePathElement::ClosePath => path.close_path(),
            }
        }
        path
    }
}

fn serialize_bez_path<S>(path: &BezPath, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    WirePath::from(path).serialize(serializer)
}

fn deserialize_bez_path<'de, D>(deserializer: D) -> Result<BezPath, D::Error>
where
    D: Deserializer<'de>,
{
    WirePath::deserialize(deserializer).map(Into::into)
}

fn deserialize_f64<'de, D>(deserializer: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    struct F64Visitor;

    impl Visitor<'_> for F64Visitor {
        type Value = f64;

        fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            formatter.write_str("an integer or floating-point number")
        }

        fn visit_f64<E>(self, value: f64) -> Result<Self::Value, E> {
            Ok(value)
        }

        fn visit_i64<E>(self, value: i64) -> Result<Self::Value, E>
        where
            E: de::Error,
        {
            Ok(value as f64)
        }

        fn visit_u64<E>(self, value: u64) -> Result<Self::Value, E>
        where
            E: de::Error,
        {
            Ok(value as f64)
        }
    }

    deserializer.deserialize_any(F64Visitor)
}

fn deserialize_optional_f64<'de, D>(deserializer: D) -> Result<Option<f64>, D::Error>
where
    D: Deserializer<'de>,
{
    Option::<F64OrInt>::deserialize(deserializer).map(|value| value.map(|value| value.0))
}

#[derive(Debug, Clone, Copy)]
struct F64OrInt(f64);

impl<'de> Deserialize<'de> for F64OrInt {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserialize_f64(deserializer).map(Self)
    }
}

impl From<CurvePoint> for Point {
    fn from(value: CurvePoint) -> Self {
        Point::new(value.x, value.y)
    }
}

impl From<Point> for CurvePoint {
    fn from(value: Point) -> Self {
        Self {
            x: value.x,
            y: value.y,
        }
    }
}

impl From<CubicBezierSpec> for CubicBez {
    fn from(value: CubicBezierSpec) -> Self {
        CubicBez::new(value.start, value.ctrl_a, value.ctrl_b, value.end)
    }
}

impl From<CubicBez> for CubicBezierSpec {
    fn from(value: CubicBez) -> Self {
        Self {
            start: value.p0.into(),
            ctrl_a: value.p1.into(),
            ctrl_b: value.p2.into(),
            end: value.p3.into(),
        }
    }
}

impl From<CubicBezierSpec> for PathSeg {
    fn from(value: CubicBezierSpec) -> Self {
        PathSeg::Cubic(value.into())
    }
}

impl From<LineSegmentSpec> for Line {
    fn from(value: LineSegmentSpec) -> Self {
        Line::new(value.start, value.end)
    }
}

impl From<LineSegmentSpec> for PathSeg {
    fn from(value: LineSegmentSpec) -> Self {
        PathSeg::Line(value.into())
    }
}

pub fn curve_cubic_path_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: CubicPathSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize cubic path spec: {err}"))?;
    let output = cubic_path(spec)?;
    encode_cbor(&output)
}

pub fn curve_trim_path_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: TrimPathSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize path trim spec: {err}"))?;
    let output = trim_path(spec)?;
    encode_cbor(&output)
}

pub fn curve_hobby_through_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: HobbyThroughSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize Hobby curve spec: {err}"))?;
    let output = hobby_through(spec)?;
    encode_cbor(&output)
}

pub fn curve_hobby_spline_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: HobbySplineSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize Hobby spline spec: {err}"))?;
    let output = hobby_spline(spec)?;
    encode_cbor(&output)
}

pub fn curve_pattern_path_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: PatternPathSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize path pattern spec: {err}"))?;
    let output = pattern_path(spec)?;
    encode_cbor(&output)
}

pub fn curve_parallel_path_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: ParallelPathSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize parallel path spec: {err}"))?;
    let output = parallel_path(spec)?;
    encode_cbor(&output)
}

fn cubic_path(spec: CubicPathSpec) -> Result<CurvePathOutput, String> {
    let accuracy = validate_positive_accuracy(spec.accuracy)?;
    curve_path_from_segments(std::iter::once(PathSeg::Cubic(spec.curve.into())), accuracy)
}

fn trim_path(spec: TrimPathSpec) -> Result<CurvePathOutput, String> {
    let accuracy = validate_positive_accuracy(spec.accuracy)?;
    let start_outset = validate_outset(spec.start_outset, "start")?;
    let end_outset = validate_outset(spec.end_outset, "end")?;
    let segments = trim_path_segments(spec.path.segments(), start_outset, end_outset, accuracy)?;
    curve_path_from_segments(segments, accuracy)
}

fn validate_positive_accuracy(value: f64) -> Result<f64, String> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err("curve arc-length accuracy must be finite and positive".to_string())
    }
}

fn validate_outset(value: f64, end: &str) -> Result<f64, String> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(format!(
            "{end} curve outset must be finite and non-negative"
        ))
    }
}

fn pattern_path(spec: PatternPathSpec) -> Result<PatternPathOutput, String> {
    let amplitude = validate_finite(spec.amplitude, "pattern amplitude")?;
    let wavelength = validate_positive(spec.wavelength, "pattern wavelength")?;
    let phase = validate_finite(spec.phase, "pattern phase")?;
    validate_finite(spec.coil_longitudinal_scale, "coil longitudinal scale")?;
    let accuracy = validate_positive_accuracy(spec.accuracy)?;
    if spec.samples_per_period == 0 {
        return Err("pattern samples-per-period must be positive".to_string());
    }
    let pattern = PointPattern::from_input(spec.pattern)?;

    let path_segments = spec.path.segments().collect::<Vec<_>>();
    let segment_lengths = path_segments
        .iter()
        .map(|segment| segment.arclen(accuracy))
        .collect::<Vec<_>>();
    let length: f64 = segment_lengths.iter().sum();
    if length <= f64::EPSILON {
        let point = path_segments
            .first()
            .map(PathSeg::start)
            .unwrap_or(Point::new(0.0, 0.0));
        return Ok(PatternPathOutput {
            path: spec.path,
            pattern: pattern.name.clone(),
            length,
            points: vec![point.into()],
            curves: Vec::new(),
            segments: Vec::new(),
        });
    }

    let distances = pattern.distances(length, wavelength, spec.samples_per_period);
    let mut points = Vec::with_capacity(distances.len());
    for distance in distances {
        let (segment, segment_distance) =
            path_segment_at_distance(&path_segments, &segment_lengths, distance);
        let t = segment.inv_arclen(segment_distance, accuracy);
        let base = segment.eval(t);
        let tangent = normalized_tangent_between(
            path_seg_tangent(segment, t),
            segment.start(),
            segment.end(),
        );
        let normal = Vec2::new(-tangent.y, tangent.x);
        let pattern_point = pattern.evaluate(distance / wavelength + phase / std::f64::consts::TAU);
        let tapered_amplitude = amplitude
            * pattern_endpoint_envelope(
                &pattern,
                distance,
                length,
                wavelength,
                spec.anchor_start,
                spec.anchor_end,
            );
        let longitudinal = tapered_amplitude * pattern_point.x;
        let lateral = tapered_amplitude * pattern_point.y;
        points.push((base + tangent * longitudinal + normal * lateral).into());
    }
    if !pattern.endpoint_ramp && spec.anchor_start {
        points[0] = path_segments
            .first()
            .map(PathSeg::start)
            .unwrap_or(Point::new(0.0, 0.0))
            .into();
    }
    if !pattern.endpoint_ramp && spec.anchor_end {
        let last_index = points.len() - 1;
        points[last_index] = path_segments
            .last()
            .map(PathSeg::end)
            .unwrap_or(Point::new(0.0, 0.0))
            .into();
    }

    let segments = line_segments(&points);
    let curves = match pattern.interpolation {
        PatternInterpolation::Linear => Vec::new(),
        PatternInterpolation::Smooth => cubic_spline_through_points(&points),
    };
    let path = if curves.is_empty() {
        BezPath::from_path_segments(segments.iter().copied().map(Into::into))
    } else {
        BezPath::from_path_segments(curves.iter().copied().map(Into::into))
    };

    Ok(PatternPathOutput {
        path,
        pattern: pattern.name,
        length,
        curves,
        segments,
        points,
    })
}

fn parallel_path(spec: ParallelPathSpec) -> Result<CurvePathOutput, String> {
    let distance = validate_finite(spec.distance, "parallel path distance")?;
    let accuracy = validate_positive_accuracy(spec.accuracy)?;
    let start_outset = validate_outset(spec.start_outset, "start")?;
    let end_outset = validate_outset(spec.end_outset, "end")?;
    let source_length: f64 = spec
        .path
        .segments()
        .map(|segment| segment.arclen(accuracy))
        .sum();
    if source_length <= f64::EPSILON {
        let point = spec
            .path
            .segments()
            .next()
            .map(|segment| segment.start().into())
            .unwrap_or(CurvePoint { x: 0.0, y: 0.0 });
        return Ok(CurvePathOutput {
            path: spec.path,
            length: 0.0,
            points: vec![point],
            curves: Vec::new(),
            segments: Vec::new(),
        });
    }

    let offset_segments = spec
        .path
        .segments()
        .flat_map(|segment| {
            let offset = CubicOffset::new_regularized(segment.to_cubic(), distance, accuracy);
            let path = if spec.optimize {
                fit_to_bezpath_opt(&offset, accuracy)
            } else {
                fit_to_bezpath(&offset, accuracy)
            };
            path.segments().collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let segments = trim_path_segments(offset_segments, start_outset, end_outset, accuracy)?;
    curve_path_from_segments(segments, accuracy)
}

#[derive(Debug, Clone)]
struct PointPattern {
    name: String,
    points: Vec<ResolvedPatternPoint>,
    interpolation: PatternInterpolation,
    endpoint_ramp: bool,
}

#[derive(Debug, Clone, Copy)]
struct ResolvedPatternPoint {
    at: f64,
    x: f64,
    y: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PatternInterpolation {
    Linear,
    Smooth,
}

impl PointPattern {
    fn from_input(input: PatternInput) -> Result<Self, String> {
        match input {
            PatternInput::Name(name) => Err(format!(
                "Unsupported unresolved path pattern `{}`; use a point pattern object",
                name
            )),
            PatternInput::Points(input) => Self::from_points(input),
        }
    }

    fn from_points(input: PointPatternInput) -> Result<Self, String> {
        if !input.kind.trim().eq_ignore_ascii_case("points") {
            return Err(format!("Unsupported point pattern kind: {}", input.kind));
        }
        if input.points.len() < 2 {
            return Err("point pattern requires at least two points".to_string());
        }

        let interpolation = PatternInterpolation::parse(&input.interpolation)?;
        let denom = (input.points.len() - 1) as f64;
        let mut points = input
            .points
            .into_iter()
            .enumerate()
            .map(|(index, point)| {
                let at = point.at.unwrap_or(index as f64 / denom);
                validate_pattern_point(at, point.x, point.y)?;
                Ok(ResolvedPatternPoint {
                    at,
                    x: point.x,
                    y: point.y,
                })
            })
            .collect::<Result<Vec<_>, String>>()?;
        points.sort_by(|a, b| a.at.partial_cmp(&b.at).unwrap());
        if points.first().is_some_and(|point| point.at > 0.0) {
            let first = *points.first().unwrap();
            points.insert(
                0,
                ResolvedPatternPoint {
                    at: 0.0,
                    x: first.x,
                    y: first.y,
                },
            );
        }
        if points.last().is_some_and(|point| point.at < 1.0) {
            let last = *points.last().unwrap();
            points.push(ResolvedPatternPoint {
                at: 1.0,
                x: last.x,
                y: last.y,
            });
        }

        Ok(Self {
            name: input.name.unwrap_or_else(|| "points".to_string()),
            points,
            interpolation,
            endpoint_ramp: input.endpoint_ramp,
        })
    }

    fn distances(&self, length: f64, wavelength: f64, samples_per_period: usize) -> Vec<f64> {
        match self.interpolation {
            PatternInterpolation::Smooth => {
                sampled_distances(length, wavelength, samples_per_period)
            }
            PatternInterpolation::Linear => self.linear_distances(length, wavelength),
        }
    }

    fn linear_distances(&self, length: f64, wavelength: f64) -> Vec<f64> {
        let periods = (length / wavelength).ceil().max(1.0) as usize;
        let mut distances = Vec::new();
        for period in 0..=periods {
            let offset = period as f64 * wavelength;
            for point in &self.points {
                let distance = offset + point.at * wavelength;
                if distance <= length {
                    distances.push(distance);
                }
            }
        }
        distances.push(length);
        distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
        distances.dedup_by(|a, b| (*a - *b).abs() <= f64::EPSILON);
        distances
    }

    fn evaluate(&self, period_position: f64) -> ResolvedPatternPoint {
        let at = period_position.rem_euclid(1.0);
        let next_index = self
            .points
            .iter()
            .position(|point| point.at >= at)
            .unwrap_or(0);
        let (start, end) = if next_index == 0 {
            (*self.points.last().unwrap(), self.points[0])
        } else {
            (self.points[next_index - 1], self.points[next_index])
        };
        let span = if end.at >= start.at {
            end.at - start.at
        } else {
            end.at + 1.0 - start.at
        };
        if span <= f64::EPSILON {
            return end;
        }

        let local_at = if at >= start.at {
            at - start.at
        } else {
            at + 1.0 - start.at
        };
        let t = local_at / span;
        ResolvedPatternPoint {
            at,
            x: start.x + (end.x - start.x) * t,
            y: start.y + (end.y - start.y) * t,
        }
    }
}

fn pattern_endpoint_envelope(
    pattern: &PointPattern,
    distance: f64,
    length: f64,
    wavelength: f64,
    anchor_start: bool,
    anchor_end: bool,
) -> f64 {
    if !pattern.endpoint_ramp {
        return 1.0;
    }

    let ramp = (wavelength * 0.5).min(length * 0.5);
    if ramp <= f64::EPSILON {
        return 1.0;
    }

    let mut envelope: f64 = 1.0;
    if anchor_start {
        envelope = envelope.min(smoothstep((distance / ramp).clamp(0.0, 1.0)));
    }
    if anchor_end {
        envelope = envelope.min(smoothstep(((length - distance) / ramp).clamp(0.0, 1.0)));
    }
    envelope
}

impl PatternInterpolation {
    fn parse(value: &str) -> Result<Self, String> {
        match value.trim().to_ascii_lowercase().as_str() {
            "linear" | "line" => Ok(Self::Linear),
            "smooth" | "spline" | "cubic" => Ok(Self::Smooth),
            other => Err(format!("Unsupported point pattern interpolation: {other}")),
        }
    }
}

fn sampled_distances(length: f64, wavelength: f64, samples_per_period: usize) -> Vec<f64> {
    let step = wavelength / samples_per_period as f64;
    let count = (length / step).ceil().max(1.0) as usize;
    let mut distances = Vec::with_capacity(count + 1);
    for i in 0..=count {
        distances.push((i as f64 * step).min(length));
    }
    if distances
        .last()
        .is_none_or(|last| (length - last).abs() > f64::EPSILON)
    {
        distances.push(length);
    }
    distances.dedup_by(|a, b| (*a - *b).abs() <= f64::EPSILON);
    distances
}

fn validate_pattern_point(at: f64, x: f64, y: f64) -> Result<(), String> {
    if !(at.is_finite() && (0.0..=1.0).contains(&at)) {
        return Err("point pattern `at` values must be finite and between 0 and 1".to_string());
    }
    validate_finite(x, "point pattern x")?;
    validate_finite(y, "point pattern y")?;
    Ok(())
}

fn smoothstep(t: f64) -> f64 {
    t * t * (3.0 - 2.0 * t)
}

fn normalized_tangent_between(tangent: Vec2, start: Point, end: Point) -> Vec2 {
    if tangent.hypot2() > f64::EPSILON {
        tangent.normalize()
    } else {
        let chord = end - start;
        if chord.hypot2() > f64::EPSILON {
            chord.normalize()
        } else {
            Vec2::new(1.0, 0.0)
        }
    }
}

fn path_segment_at_distance<'a>(
    segments: &'a [PathSeg],
    lengths: &[f64],
    distance: f64,
) -> (&'a PathSeg, f64) {
    let mut cursor = 0.0;
    for (index, (segment, segment_length)) in segments.iter().zip(lengths).enumerate() {
        if distance <= cursor + segment_length || index + 1 == segments.len() {
            return (segment, (distance - cursor).clamp(0.0, *segment_length));
        }
        cursor += segment_length;
    }
    let last = segments
        .last()
        .expect("path_segment_at_distance requires a non-empty segment list");
    (last, 0.0)
}

fn path_seg_tangent(segment: &PathSeg, t: f64) -> Vec2 {
    match *segment {
        PathSeg::Line(line) => line.p1 - line.p0,
        PathSeg::Quad(quad) => quad.deriv().eval(t).to_vec2(),
        PathSeg::Cubic(cubic) => cubic.deriv().eval(t).to_vec2(),
    }
}

fn line_segments(points: &[CurvePoint]) -> Vec<LineSegmentSpec> {
    points
        .windows(2)
        .map(|window| LineSegmentSpec {
            start: window[0],
            end: window[1],
        })
        .collect()
}

fn curve_path_from_segments(
    segments: impl IntoIterator<Item = PathSeg>,
    accuracy: f64,
) -> Result<CurvePathOutput, String> {
    let path = BezPath::from_path_segments(segments.into_iter());
    let mut length = 0.0;
    let mut points = Vec::new();
    let mut curves = Vec::new();
    let mut line_segments = Vec::new();

    for segment in path.segments() {
        length += segment.arclen(accuracy);
        if points.is_empty() {
            points.push(segment.start().into());
        }
        points.push(segment.end().into());
        match segment {
            PathSeg::Line(line) => line_segments.push(LineSegmentSpec {
                start: line.p0.into(),
                end: line.p1.into(),
            }),
            PathSeg::Quad(_) | PathSeg::Cubic(_) => curves.push(segment.to_cubic().into()),
        }
    }

    if points.is_empty() {
        return Err("path fitting produced no visible path".to_string());
    }

    Ok(CurvePathOutput {
        path,
        length,
        points,
        curves,
        segments: line_segments,
    })
}

fn trim_path_segments(
    segments: impl IntoIterator<Item = PathSeg>,
    start_outset: f64,
    end_outset: f64,
    accuracy: f64,
) -> Result<Vec<PathSeg>, String> {
    let segments: Vec<_> = segments.into_iter().collect();
    if start_outset == 0.0 && end_outset == 0.0 {
        return Ok(segments);
    }

    let lengths: Vec<_> = segments
        .iter()
        .map(|segment| segment.arclen(accuracy))
        .collect();
    let length: f64 = lengths.iter().sum();
    if length <= f64::EPSILON {
        return Ok(segments);
    }

    let (start_outset, end_outset) = fit_outsets_to_length(start_outset, end_outset, length);
    let visible_start = start_outset;
    let visible_end = length - end_outset;
    if visible_end <= visible_start {
        return Err("path trim distances leave no visible path".to_string());
    }

    let mut cursor = 0.0;
    let mut trimmed = Vec::new();
    for (segment, segment_length) in segments.into_iter().zip(lengths) {
        let segment_start = cursor;
        let segment_end = cursor + segment_length;
        cursor = segment_end;

        let keep_start = visible_start.max(segment_start);
        let keep_end = visible_end.min(segment_end);
        if keep_end <= keep_start {
            continue;
        }

        let local_start = keep_start - segment_start;
        let local_end = keep_end - segment_start;
        let t0 = if local_start <= f64::EPSILON {
            0.0
        } else {
            segment.inv_arclen(local_start, accuracy)
        };
        let t1 = if segment_length - local_end <= f64::EPSILON {
            1.0
        } else {
            segment.inv_arclen(local_end, accuracy)
        };
        if t1 > t0 {
            trimmed.push(segment.subsegment(t0..t1));
        }
    }

    if trimmed.is_empty() {
        return Err("path trimming produced no visible path".to_string());
    }
    Ok(trimmed)
}

fn fit_outsets_to_length(start_outset: f64, end_outset: f64, length: f64) -> (f64, f64) {
    let total = start_outset + end_outset;
    let max_total = length * (1.0 - 1e-9);
    if total > max_total && total > 0.0 {
        let scale = max_total / total;
        (start_outset * scale, end_outset * scale)
    } else {
        (start_outset, end_outset)
    }
}

fn cubic_spline_through_points(points: &[CurvePoint]) -> Vec<CubicBezierSpec> {
    if points.len() < 2 {
        return Vec::new();
    }

    let points: Vec<Point> = points.iter().copied().map(Into::into).collect();
    let tangents: Vec<Vec2> = (0..points.len())
        .map(|i| {
            if i == 0 {
                points[1] - points[0]
            } else if i == points.len() - 1 {
                points[i] - points[i - 1]
            } else {
                (points[i + 1] - points[i - 1]) * 0.5
            }
        })
        .collect();

    points
        .windows(2)
        .zip(tangents.windows(2))
        .map(|(point_pair, tangent_pair)| {
            CubicBez::new(
                point_pair[0],
                point_pair[0] + tangent_pair[0] / 3.0,
                point_pair[1] - tangent_pair[1] / 3.0,
                point_pair[1],
            )
            .into()
        })
        .collect()
}

fn validate_finite(value: f64, name: &str) -> Result<f64, String> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(format!("{name} must be finite"))
    }
}

fn validate_positive(value: f64, name: &str) -> Result<f64, String> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(format!("{name} must be finite and positive"))
    }
}

fn hobby_through(spec: HobbyThroughSpec) -> Result<CurvePathOutput, String> {
    hobby_to_cubic_open(
        &[spec.start.into(), spec.through.into(), spec.end.into()],
        spec.omega,
    )
    .and_then(|segments| cubic_spline_output(segments, spec.accuracy))
}

fn hobby_spline(spec: HobbySplineSpec) -> Result<CurvePathOutput, String> {
    let accuracy = spec.accuracy;
    let points = spec.points.into_iter().map(Into::into).collect::<Vec<_>>();
    hobby_to_cubic_open(&points, spec.omega)
        .and_then(|segments| cubic_spline_output(segments, accuracy))
}

fn cubic_spline_output(segments: Vec<CubicBez>, accuracy: f64) -> Result<CurvePathOutput, String> {
    let accuracy = validate_positive_accuracy(accuracy)?;
    curve_path_from_segments(segments.into_iter().map(PathSeg::Cubic), accuracy)
}

fn hobby_to_cubic_open(points: &[Point], omega: f64) -> Result<Vec<CubicBez>, String> {
    if points.len() < 2 {
        return Err("Hobby curve requires at least two points".to_string());
    }
    if !omega.is_finite() || omega < 0.0 {
        return Err("Hobby omega must be finite and non-negative".to_string());
    }
    if points.len() == 2 {
        let start = points[0];
        let end = points[1];
        return Ok(vec![CubicBez::new(start, start, end, end)]);
    }

    let n = points.len() - 1;
    let mut chords = Vec::with_capacity(n);
    let mut distances = Vec::with_capacity(n);
    for segment in points.windows(2) {
        let chord = segment[1] - segment[0];
        let distance = chord.hypot();
        if distance <= f64::EPSILON {
            return Err("Hobby curve points must be distinct".to_string());
        }
        chords.push(chord);
        distances.push(distance);
    }

    let mut gamma = vec![0.0; n + 1];
    for i in 1..n {
        gamma[i] = signed_angle(chords[i - 1], chords[i]);
    }

    // This is CeTZ's open Hobby spline system with default unit tensions.
    let mut lower = vec![0.0; n + 1];
    let mut diagonal = vec![0.0; n + 1];
    let mut upper = vec![0.0; n + 1];
    let mut rhs = vec![0.0; n + 1];

    let c0 = omega + 2.0;
    let d0 = 2.0 * omega + 1.0;
    diagonal[0] = c0;
    upper[0] = d0;
    rhs[0] = -d0 * gamma[1];

    for i in 1..n {
        let prev = distances[i - 1];
        let next = distances[i];
        lower[i] = 1.0 / prev;
        let b = 2.0 / prev;
        let c = 2.0 / next;
        upper[i] = 1.0 / next;
        diagonal[i] = b + c;
        rhs[i] = -b * gamma[i] - upper[i] * gamma[i + 1];
    }

    lower[n] = 2.0 * omega + 1.0;
    diagonal[n] = omega + 2.0;

    let alpha = solve_tridiagonal(&lower, &diagonal, &upper, &rhs)?;
    let beta: Vec<_> = (0..n).map(|i| -(alpha[i + 1] + gamma[i + 1])).collect();

    let mut cubics = Vec::with_capacity(n);
    for i in 0..n {
        let start = points[i];
        let end = points[i + 1];
        let chord = chords[i];
        let ctrl_a = start + rotate(chord, alpha[i]) * (hobby_rho(alpha[i], beta[i]) / 3.0);
        let ctrl_b = end - rotate(chord, -beta[i]) * (hobby_rho(beta[i], alpha[i]) / 3.0);
        cubics.push(CubicBez::new(start, ctrl_a, ctrl_b, end));
    }

    Ok(cubics)
}

fn signed_angle(a: kurbo::Vec2, b: kurbo::Vec2) -> f64 {
    (a.x * b.y - a.y * b.x).atan2(a.x * b.x + a.y * b.y)
}

fn rotate(v: kurbo::Vec2, angle: f64) -> kurbo::Vec2 {
    let (sin, cos) = angle.sin_cos();
    kurbo::Vec2::new(v.x * cos - v.y * sin, v.x * sin + v.y * cos)
}

fn hobby_rho(alpha: f64, beta: f64) -> f64 {
    let sqrt_2 = 2.0_f64.sqrt();
    let sqrt_5 = 5.0_f64.sqrt();
    let numerator = 2.0
        + sqrt_2
            * (alpha.sin() - beta.sin() / 16.0)
            * (beta.sin() - alpha.sin() / 16.0)
            * (alpha.cos() - beta.cos());
    let denominator = 1.0 + alpha.cos() * (sqrt_5 - 1.0) / 2.0 + beta.cos() * (3.0 - sqrt_5) / 2.0;
    numerator / denominator
}

fn solve_tridiagonal(
    lower: &[f64],
    diagonal: &[f64],
    upper: &[f64],
    rhs: &[f64],
) -> Result<Vec<f64>, String> {
    let n = diagonal.len();
    if n == 0 || lower.len() != n || upper.len() != n || rhs.len() != n {
        return Err("invalid tridiagonal system dimensions".to_string());
    }

    let mut diagonal = diagonal.to_vec();
    let mut rhs = rhs.to_vec();
    for i in 1..n {
        if diagonal[i - 1].abs() <= f64::EPSILON {
            return Err("singular Hobby spline system".to_string());
        }
        let w = lower[i] / diagonal[i - 1];
        diagonal[i] -= w * upper[i - 1];
        rhs[i] -= w * rhs[i - 1];
    }

    if diagonal[n - 1].abs() <= f64::EPSILON {
        return Err("singular Hobby spline system".to_string());
    }

    let mut solution = vec![0.0; n];
    solution[n - 1] = rhs[n - 1] / diagonal[n - 1];
    for i in (0..n - 1).rev() {
        if diagonal[i].abs() <= f64::EPSILON {
            return Err("singular Hobby spline system".to_string());
        }
        solution[i] = (rhs[i] - upper[i] * solution[i + 1]) / diagonal[i];
    }
    Ok(solution)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point(x: f64, y: f64) -> CurvePoint {
        CurvePoint { x, y }
    }

    #[test]
    fn cubic_path_returns_path_dictionary() {
        let output = cubic_path(CubicPathSpec {
            curve: straight_curve(3.0),
            accuracy: 1e-6,
        })
        .unwrap();

        assert!((output.length - 3.0).abs() < 1e-6);
        assert_eq!(output.points.first().copied(), Some(point(0.0, 0.0)));
        assert_eq!(output.points.last().copied(), Some(point(3.0, 0.0)));
        assert_eq!(output.curves.len(), 1);
    }

    #[test]
    fn trim_path_trims_by_curve_arclength() {
        let output = trim_path(TrimPathSpec {
            path: straight_path(3.0),
            start_outset: 0.25,
            end_outset: 0.5,
            accuracy: 1e-6,
        })
        .unwrap();

        assert!((output.points.first().unwrap().x - 0.25).abs() < 1e-6);
        assert_eq!(output.points.first().unwrap().y, 0.0);
        assert!((output.points.last().unwrap().x - 2.5).abs() < 1e-6);
        assert_eq!(output.points.last().unwrap().y, 0.0);
    }

    #[test]
    fn parallel_path_offsets_straight_curve_to_left_normal() {
        let output = parallel_path(ParallelPathSpec {
            path: straight_path(3.0),
            distance: 0.5,
            start_outset: 0.0,
            end_outset: 0.0,
            accuracy: 1e-6,
            optimize: true,
        })
        .unwrap();

        assert!(output.points.len() >= 2);
        assert!((output.points.first().unwrap().x - 0.0).abs() < 1e-6);
        assert!((output.points.first().unwrap().y - 0.5).abs() < 1e-6);
        assert!((output.points.last().unwrap().x - 3.0).abs() < 1e-6);
        assert!((output.points.last().unwrap().y - 0.5).abs() < 1e-6);
        assert!(!output.curves.is_empty() || !output.segments.is_empty());
    }

    #[test]
    fn parallel_path_accepts_negative_distance() {
        let output = parallel_path(ParallelPathSpec {
            path: straight_path(3.0),
            distance: -0.25,
            start_outset: 0.0,
            end_outset: 0.0,
            accuracy: 1e-6,
            optimize: false,
        })
        .unwrap();

        assert!((output.points.first().unwrap().y + 0.25).abs() < 1e-6);
        assert!((output.points.last().unwrap().y + 0.25).abs() < 1e-6);
    }

    #[test]
    fn parallel_path_trims_offset_path_by_arclength() {
        let output = parallel_path(ParallelPathSpec {
            path: straight_path(4.0),
            distance: 0.25,
            start_outset: 1.0,
            end_outset: 1.0,
            accuracy: 1e-6,
            optimize: false,
        })
        .unwrap();

        assert!((output.points.first().unwrap().x - 1.0).abs() < 1e-6);
        assert!((output.points.first().unwrap().y - 0.25).abs() < 1e-6);
        assert!((output.points.last().unwrap().x - 3.0).abs() < 1e-6);
        assert!((output.points.last().unwrap().y - 0.25).abs() < 1e-6);
        assert!((output.length - 2.0).abs() < 1e-6);
    }

    fn straight_curve(length: f64) -> CubicBezierSpec {
        CubicBezierSpec {
            start: point(0.0, 0.0),
            ctrl_a: point(length / 3.0, 0.0),
            ctrl_b: point(2.0 * length / 3.0, 0.0),
            end: point(length, 0.0),
        }
    }

    fn straight_path(length: f64) -> BezPath {
        BezPath::from_path_segments(std::iter::once(PathSeg::Cubic(
            straight_curve(length).into(),
        )))
    }

    fn nearest_x(points: &[CurvePoint], x: f64) -> CurvePoint {
        *points
            .iter()
            .min_by(|a, b| (a.x - x).abs().partial_cmp(&(b.x - x).abs()).unwrap())
            .unwrap()
    }

    fn point_pattern(
        name: &str,
        interpolation: &str,
        points: Vec<PatternPointInput>,
        endpoint_ramp: bool,
    ) -> PatternInput {
        PatternInput::Points(PointPatternInput {
            kind: "points".to_string(),
            name: Some(name.to_string()),
            points,
            interpolation: interpolation.to_string(),
            endpoint_ramp,
        })
    }

    fn sampled_pattern(
        name: &str,
        samples_per_period: usize,
        mut offset: impl FnMut(f64) -> (f64, f64),
        endpoint_ramp: bool,
    ) -> PatternInput {
        let points = (0..=samples_per_period)
            .map(|index| {
                let at = index as f64 / samples_per_period as f64;
                let (x, y) = offset(std::f64::consts::TAU * at);
                PatternPointInput { at: Some(at), x, y }
            })
            .collect();
        point_pattern(name, "smooth", points, endpoint_ramp)
    }

    fn wave_pattern(samples_per_period: usize) -> PatternInput {
        sampled_pattern(
            "wave",
            samples_per_period,
            |theta| (0.0, theta.sin()),
            false,
        )
    }

    fn coil_pattern(samples_per_period: usize, longitudinal_scale: f64) -> PatternInput {
        sampled_pattern(
            "coil",
            samples_per_period,
            |theta| (longitudinal_scale * theta.cos(), theta.sin()),
            true,
        )
    }

    fn zigzag_pattern() -> PatternInput {
        point_pattern(
            "zigzag",
            "linear",
            vec![
                PatternPointInput {
                    at: Some(0.0),
                    x: 0.0,
                    y: 0.0,
                },
                PatternPointInput {
                    at: Some(0.25),
                    x: 0.0,
                    y: 1.0,
                },
                PatternPointInput {
                    at: Some(0.75),
                    x: 0.0,
                    y: -1.0,
                },
                PatternPointInput {
                    at: Some(1.0),
                    x: 0.0,
                    y: 0.0,
                },
            ],
            false,
        )
    }

    #[test]
    fn wave_pattern_offsets_along_curve_normal() {
        let output = pattern_path(PatternPathSpec {
            path: straight_path(4.0),
            pattern: wave_pattern(4),
            amplitude: 0.5,
            wavelength: 4.0,
            phase: 0.0,
            samples_per_period: 4,
            coil_longitudinal_scale: default_coil_longitudinal_scale(),
            anchor_start: default_anchor_endpoint(),
            anchor_end: default_anchor_endpoint(),
            accuracy: 1e-6,
        })
        .unwrap();

        assert_eq!(output.pattern, "wave");
        assert!(output.points.len() >= 5);
        assert!((output.points[0].y - 0.0).abs() < 1e-9);
        let peak = nearest_x(&output.points, 1.0);
        assert!((peak.x - 1.0).abs() < 1e-6);
        assert!((peak.y - 0.5).abs() < 1e-6);
        assert_eq!(output.segments.len(), output.points.len() - 1);
        assert_eq!(output.curves.len(), output.points.len() - 1);
    }

    #[test]
    fn zigzag_pattern_samples_corners() {
        let output = pattern_path(PatternPathSpec {
            path: straight_path(4.0),
            pattern: zigzag_pattern(),
            amplitude: 1.0,
            wavelength: 4.0,
            phase: 0.0,
            samples_per_period: 2,
            coil_longitudinal_scale: default_coil_longitudinal_scale(),
            anchor_start: default_anchor_endpoint(),
            anchor_end: default_anchor_endpoint(),
            accuracy: 1e-6,
        })
        .unwrap();

        assert_eq!(output.pattern, "zigzag");
        assert!(output.points.len() >= 5);
        assert!((nearest_x(&output.points, 1.0).y - 1.0).abs() < 1e-9);
        assert!((nearest_x(&output.points, 3.0).y + 1.0).abs() < 1e-9);
        assert!(output.curves.is_empty());
    }

    #[test]
    fn point_pattern_interpolates_declared_points() {
        let output = pattern_path(PatternPathSpec {
            path: straight_path(4.0),
            pattern: PatternInput::Points(PointPatternInput {
                kind: "points".to_string(),
                name: None,
                interpolation: "linear".to_string(),
                endpoint_ramp: false,
                points: vec![
                    PatternPointInput {
                        at: Some(0.0),
                        x: 0.0,
                        y: 0.0,
                    },
                    PatternPointInput {
                        at: Some(0.5),
                        x: 0.0,
                        y: 1.0,
                    },
                    PatternPointInput {
                        at: Some(1.0),
                        x: 0.0,
                        y: 0.0,
                    },
                ],
            }),
            amplitude: 0.5,
            wavelength: 4.0,
            phase: 0.0,
            samples_per_period: 4,
            coil_longitudinal_scale: default_coil_longitudinal_scale(),
            anchor_start: true,
            anchor_end: true,
            accuracy: 1e-6,
        })
        .unwrap();

        assert_eq!(output.pattern, "points");
        assert!((nearest_x(&output.points, 2.0).y - 0.5).abs() < 1e-9);
        assert!(output.curves.is_empty());
    }

    #[test]
    fn coil_pattern_adds_longitudinal_and_lateral_offsets() {
        let output = pattern_path(PatternPathSpec {
            path: straight_path(4.0),
            pattern: coil_pattern(4, 0.5),
            amplitude: 0.5,
            wavelength: 4.0,
            phase: 0.0,
            samples_per_period: 4,
            coil_longitudinal_scale: 0.5,
            anchor_start: false,
            anchor_end: false,
            accuracy: 1e-6,
        })
        .unwrap();

        assert_eq!(output.pattern, "coil");
        assert!((output.points[0].x - 0.25).abs() < 1e-9);
        assert!((output.points[0].y - 0.0).abs() < 1e-9);
        assert!((output.points[1].x - 1.0).abs() < 1e-6);
        assert!((output.points[1].y - 0.5).abs() < 1e-6);
        assert_eq!(output.curves.len(), output.points.len() - 1);
    }

    #[test]
    fn coil_default_turns_back_over_the_baseline() {
        let output = pattern_path(PatternPathSpec {
            path: straight_path(2.0),
            pattern: coil_pattern(16, default_coil_longitudinal_scale()),
            amplitude: 0.08,
            wavelength: 0.55,
            phase: 0.0,
            samples_per_period: 16,
            coil_longitudinal_scale: default_coil_longitudinal_scale(),
            anchor_start: default_anchor_endpoint(),
            anchor_end: default_anchor_endpoint(),
            accuracy: 1e-6,
        })
        .unwrap();

        assert!(
            output
                .points
                .windows(2)
                .any(|window| window[1].x < window[0].x)
        );
    }

    #[test]
    fn coil_pattern_anchors_requested_endpoints() {
        let curve = straight_curve(2.0);
        let output = pattern_path(PatternPathSpec {
            path: straight_path(2.0),
            pattern: coil_pattern(16, default_coil_longitudinal_scale()),
            amplitude: 0.08,
            wavelength: 0.55,
            phase: 0.0,
            samples_per_period: 16,
            coil_longitudinal_scale: default_coil_longitudinal_scale(),
            anchor_start: true,
            anchor_end: true,
            accuracy: 1e-6,
        })
        .unwrap();

        assert_eq!(output.points[0], curve.start);
        assert_eq!(*output.points.last().unwrap(), curve.end);
        assert!(output.points[1].y.abs() < 0.08);
    }

    #[test]
    fn hobby_through_returns_two_smooth_cubic_segments() {
        let output = hobby_through(HobbyThroughSpec {
            start: point(0.0, 0.0),
            through: point(1.0, 1.0),
            end: point(2.0, 0.0),
            omega: 1.0,
            accuracy: 1e-6,
        })
        .unwrap();

        assert_eq!(output.curves.len(), 2);
        assert_eq!(output.curves[0].start, point(0.0, 0.0));
        assert_eq!(output.curves[0].end, point(1.0, 1.0));
        assert_eq!(output.curves[1].start, point(1.0, 1.0));
        assert_eq!(output.curves[1].end, point(2.0, 0.0));

        let incoming = Point::from(output.curves[0].end) - Point::from(output.curves[0].ctrl_b);
        let outgoing = Point::from(output.curves[1].ctrl_a) - Point::from(output.curves[1].start);
        assert!(signed_angle(incoming, outgoing).abs() < 1e-12);
    }

    #[test]
    fn hobby_spline_returns_smooth_cubic_segments_through_all_points() {
        let points = vec![
            point(0.0, 0.0),
            point(1.0, 1.0),
            point(2.0, -0.4),
            point(3.0, 0.2),
        ];
        let output = hobby_spline(HobbySplineSpec {
            points: points.clone(),
            omega: 1.0,
            accuracy: 1e-3,
        })
        .unwrap();

        assert_eq!(output.curves.len(), points.len() - 1);
        for (index, segment) in output.curves.iter().enumerate() {
            assert_eq!(segment.start, points[index]);
            assert_eq!(segment.end, points[index + 1]);
        }

        for pair in output.curves.windows(2) {
            let incoming = Point::from(pair[0].end) - Point::from(pair[0].ctrl_b);
            let outgoing = Point::from(pair[1].ctrl_a) - Point::from(pair[1].start);
            assert!(signed_angle(incoming, outgoing).abs() < 1e-12);
        }
    }

    #[test]
    fn cbor_api_returns_explicit_wire_path() {
        let input = HobbyThroughSpec {
            start: point(0.0, 0.0),
            through: point(1.0, 1.0),
            end: point(2.0, 0.0),
            omega: 1.0,
            accuracy: 1e-6,
        };
        let bytes = encode_cbor(&input).unwrap();
        let output: CurvePathOutput =
            ciborium::de::from_reader(&curve_hobby_through_bytes(&bytes).unwrap()[..]).unwrap();

        assert_eq!(output.curves[0].start, point(0.0, 0.0));
        assert_eq!(output.curves[0].end, point(1.0, 1.0));
        assert_eq!(output.curves[1].start, point(1.0, 1.0));
        assert_eq!(output.curves[1].end, point(2.0, 0.0));
        assert_eq!(output.path.segments().count(), 2);
    }
}
