use kurbo::{CubicBez, ParamCurve, ParamCurveArclen, ParamCurveDeriv, Point, QuadBez, Vec2};
use serde::{
    Deserialize, Deserializer, Serialize,
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
pub struct SplitCubicSpec {
    #[serde(flatten)]
    pub curve: CubicBezierSpec,
    #[serde(default = "default_split_t")]
    pub t: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct TrimCubicSpec {
    #[serde(flatten)]
    pub curve: CubicBezierSpec,
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
pub struct PatternCubicSpec {
    #[serde(flatten)]
    pub curve: CubicBezierSpec,
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
pub struct QuadraticThroughSpec {
    pub start: CurvePoint,
    pub through: CurvePoint,
    pub end: CurvePoint,
    #[serde(default = "default_split_t")]
    pub t: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct HobbyThroughSpec {
    pub start: CurvePoint,
    pub through: CurvePoint,
    pub end: CurvePoint,
    #[serde(default = "default_hobby_omega")]
    pub omega: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SplitCurveOutput {
    pub curve: CubicBezierSpec,
    pub segments: [CubicBezierSpec; 2],
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CubicSplineOutput {
    pub segments: Vec<CubicBezierSpec>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct LineSegmentSpec {
    pub start: CurvePoint,
    pub end: CurvePoint,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PatternPathOutput {
    pub pattern: String,
    pub length: f64,
    pub points: Vec<CurvePoint>,
    pub curves: Vec<CubicBezierSpec>,
    pub segments: Vec<LineSegmentSpec>,
}

fn default_split_t() -> f64 {
    0.5
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

fn encode_cbor<T: Serialize>(value: &T) -> Result<Vec<u8>, String> {
    let mut bytes = Vec::new();
    ciborium::ser::into_writer(value, &mut bytes)
        .map_err(|err| format!("Failed to serialize CBOR value: {err}"))?;
    Ok(bytes)
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

pub fn curve_split_cubic_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: SplitCubicSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize cubic curve spec: {err}"))?;
    let output = split_cubic(spec)?;
    encode_cbor(&output)
}

pub fn curve_split_quadratic_through_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: QuadraticThroughSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize quadratic-through curve spec: {err}"))?;
    let curve = quadratic_through_to_cubic(spec)?;
    split_cubic(SplitCubicSpec { curve, t: spec.t }).and_then(|output| encode_cbor(&output))
}

pub fn curve_trim_cubic_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: TrimCubicSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize cubic curve trim spec: {err}"))?;
    let output = trim_cubic(spec)?;
    encode_cbor(&output)
}

pub fn curve_hobby_through_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: HobbyThroughSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize Hobby curve spec: {err}"))?;
    let output = hobby_through(spec)?;
    encode_cbor(&output)
}

pub fn curve_pattern_cubic_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: PatternCubicSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize cubic curve pattern spec: {err}"))?;
    let output = pattern_cubic(spec)?;
    encode_cbor(&output)
}

fn split_cubic(spec: SplitCubicSpec) -> Result<SplitCurveOutput, String> {
    let t = validate_t(spec.t)?;
    let curve = CubicBez::from(spec.curve);
    let first = curve.subsegment(0.0..t);
    let second = curve.subsegment(t..1.0);

    Ok(SplitCurveOutput {
        curve: spec.curve,
        segments: [first.into(), second.into()],
    })
}

fn trim_cubic(spec: TrimCubicSpec) -> Result<CubicBezierSpec, String> {
    let curve = CubicBez::from(spec.curve);
    let accuracy = validate_positive_accuracy(spec.accuracy)?;
    let length = curve.arclen(accuracy);
    if length <= f64::EPSILON {
        return Ok(spec.curve);
    }

    let start_outset = validate_outset(spec.start_outset, "start")?;
    let end_outset = validate_outset(spec.end_outset, "end")?;
    let max_outset = length * 0.45;
    let start_outset = start_outset.min(max_outset);
    let end_outset = end_outset.min(max_outset);

    let t0 = curve.inv_arclen(start_outset, accuracy);
    let t1 = curve.inv_arclen(length - end_outset, accuracy);
    if t1 <= t0 {
        return Err("cubic trim distances leave no visible curve segment".to_string());
    }

    Ok(curve.subsegment(t0..t1).into())
}

fn quadratic_through_to_cubic(spec: QuadraticThroughSpec) -> Result<CubicBezierSpec, String> {
    let t = validate_t(spec.t)?;
    let denom = 2.0 * t * (1.0 - t);
    if denom <= f64::EPSILON {
        return Err(
            "quadratic-through split parameter must be strictly between 0 and 1".to_string(),
        );
    }

    let start = Point::from(spec.start);
    let through = Point::from(spec.through);
    let end = Point::from(spec.end);
    let ctrl =
        (through.to_vec2() - start.to_vec2() * (1.0 - t).powi(2) - end.to_vec2() * t.powi(2))
            / denom;

    Ok(QuadBez::new(start, ctrl.to_point(), end).raise().into())
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

fn pattern_cubic(spec: PatternCubicSpec) -> Result<PatternPathOutput, String> {
    let amplitude = validate_finite(spec.amplitude, "pattern amplitude")?;
    let wavelength = validate_positive(spec.wavelength, "pattern wavelength")?;
    let phase = validate_finite(spec.phase, "pattern phase")?;
    validate_finite(spec.coil_longitudinal_scale, "coil longitudinal scale")?;
    let accuracy = validate_positive_accuracy(spec.accuracy)?;
    if spec.samples_per_period == 0 {
        return Err("pattern samples-per-period must be positive".to_string());
    }
    let pattern = PointPattern::from_input(spec.pattern)?;

    let curve = CubicBez::from(spec.curve);
    let length = curve.arclen(accuracy);
    if length <= f64::EPSILON {
        return Ok(PatternPathOutput {
            pattern: pattern.name.clone(),
            length,
            points: vec![spec.curve.start],
            curves: Vec::new(),
            segments: Vec::new(),
        });
    }

    let distances = pattern.distances(length, wavelength, spec.samples_per_period);
    let deriv = curve.deriv();
    let mut points = Vec::with_capacity(distances.len());
    for distance in distances {
        let t = curve.inv_arclen(distance, accuracy);
        let base = curve.eval(t);
        let tangent = normalized_tangent(deriv.eval(t).to_vec2(), spec.curve);
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
        points[0] = spec.curve.start;
    }
    if !pattern.endpoint_ramp && spec.anchor_end {
        let last_index = points.len() - 1;
        points[last_index] = spec.curve.end;
    }

    let segments = line_segments(&points);
    let curves = match pattern.interpolation {
        PatternInterpolation::Linear => Vec::new(),
        PatternInterpolation::Smooth => cubic_spline_through_points(&points),
    };

    Ok(PatternPathOutput {
        pattern: pattern.name,
        length,
        curves,
        segments,
        points,
    })
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

fn normalized_tangent(tangent: Vec2, curve: CubicBezierSpec) -> Vec2 {
    if tangent.hypot2() > f64::EPSILON {
        tangent.normalize()
    } else {
        let chord = Point::from(curve.end) - Point::from(curve.start);
        if chord.hypot2() > f64::EPSILON {
            chord.normalize()
        } else {
            Vec2::new(1.0, 0.0)
        }
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

fn hobby_through(spec: HobbyThroughSpec) -> Result<CubicSplineOutput, String> {
    hobby_to_cubic_open(
        &[spec.start.into(), spec.through.into(), spec.end.into()],
        spec.omega,
    )
    .map(|segments| CubicSplineOutput {
        segments: segments.into_iter().map(Into::into).collect(),
    })
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

fn validate_t(t: f64) -> Result<f64, String> {
    if !t.is_finite() {
        return Err("curve split parameter must be finite".to_string());
    }
    if !(0.0..=1.0).contains(&t) {
        return Err(format!("curve split parameter must be in [0, 1], got {t}"));
    }
    Ok(t)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point(x: f64, y: f64) -> CurvePoint {
        CurvePoint { x, y }
    }

    #[test]
    fn splits_cubic_at_requested_parameter() {
        let output = split_cubic(SplitCubicSpec {
            curve: CubicBezierSpec {
                start: point(0.0, 0.0),
                ctrl_a: point(1.0, 0.0),
                ctrl_b: point(2.0, 0.0),
                end: point(3.0, 0.0),
            },
            t: 0.5,
        })
        .unwrap();

        assert_eq!(output.segments[0].start, point(0.0, 0.0));
        assert_eq!(output.segments[0].end, point(1.5, 0.0));
        assert_eq!(output.segments[1].start, point(1.5, 0.0));
        assert_eq!(output.segments[1].end, point(3.0, 0.0));
    }

    #[test]
    fn quadratic_through_curve_passes_through_split_point() {
        let curve = quadratic_through_to_cubic(QuadraticThroughSpec {
            start: point(0.0, 0.0),
            through: point(1.0, 1.0),
            end: point(2.0, 0.0),
            t: 0.5,
        })
        .unwrap();
        let output = split_cubic(SplitCubicSpec { curve, t: 0.5 }).unwrap();

        assert_eq!(output.segments[0].end, point(1.0, 1.0));
        assert_eq!(output.segments[1].start, point(1.0, 1.0));
    }

    #[test]
    fn trims_cubic_by_curve_arclength() {
        let output = trim_cubic(TrimCubicSpec {
            curve: CubicBezierSpec {
                start: point(0.0, 0.0),
                ctrl_a: point(1.0, 0.0),
                ctrl_b: point(2.0, 0.0),
                end: point(3.0, 0.0),
            },
            start_outset: 0.25,
            end_outset: 0.5,
            accuracy: 1e-6,
        })
        .unwrap();

        assert!((output.start.x - 0.25).abs() < 1e-6);
        assert_eq!(output.start.y, 0.0);
        assert!((output.end.x - 2.5).abs() < 1e-6);
        assert_eq!(output.end.y, 0.0);
    }

    fn straight_curve(length: f64) -> CubicBezierSpec {
        CubicBezierSpec {
            start: point(0.0, 0.0),
            ctrl_a: point(length / 3.0, 0.0),
            ctrl_b: point(2.0 * length / 3.0, 0.0),
            end: point(length, 0.0),
        }
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
        let output = pattern_cubic(PatternCubicSpec {
            curve: straight_curve(4.0),
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
        let output = pattern_cubic(PatternCubicSpec {
            curve: straight_curve(4.0),
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
        let output = pattern_cubic(PatternCubicSpec {
            curve: straight_curve(4.0),
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
        let output = pattern_cubic(PatternCubicSpec {
            curve: straight_curve(4.0),
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
        let output = pattern_cubic(PatternCubicSpec {
            curve: straight_curve(2.0),
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
        let output = pattern_cubic(PatternCubicSpec {
            curve,
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
        })
        .unwrap();

        assert_eq!(output.segments.len(), 2);
        assert_eq!(output.segments[0].start, point(0.0, 0.0));
        assert_eq!(output.segments[0].end, point(1.0, 1.0));
        assert_eq!(output.segments[1].start, point(1.0, 1.0));
        assert_eq!(output.segments[1].end, point(2.0, 0.0));

        let incoming = Point::from(output.segments[0].end) - Point::from(output.segments[0].ctrl_b);
        let outgoing =
            Point::from(output.segments[1].ctrl_a) - Point::from(output.segments[1].start);
        assert!(signed_angle(incoming, outgoing).abs() < 1e-12);
    }

    #[test]
    fn cbor_api_returns_smooth_halves() {
        let input = QuadraticThroughSpec {
            start: point(0.0, 0.0),
            through: point(1.0, 1.0),
            end: point(2.0, 0.0),
            t: 0.5,
        };
        let bytes = encode_cbor(&input).unwrap();
        let output: SplitCurveOutput =
            ciborium::de::from_reader(&curve_split_quadratic_through_bytes(&bytes).unwrap()[..])
                .unwrap();

        assert_eq!(output.segments[0].start, point(0.0, 0.0));
        assert_eq!(output.segments[0].end, point(1.0, 1.0));
        assert_eq!(output.segments[1].start, point(1.0, 1.0));
        assert_eq!(output.segments[1].end, point(2.0, 0.0));
    }
}
