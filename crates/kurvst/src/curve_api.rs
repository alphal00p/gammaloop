use kurbo::{
    BezPath, Cap, CubicBez, Join, Line, ParamCurve, ParamCurveArclen, ParamCurveDeriv,
    ParamCurveExtrema, PathEl, PathSeg, Point, Stroke, StrokeOpts, Vec2, offset::offset_cubic,
};
use serde::{
    Deserialize, Deserializer, Serialize, Serializer,
    de::{self, SeqAccess, Visitor},
    ser::SerializeTuple,
};
use std::{cell::OnceCell, fmt};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurvePoint {
    pub x: f64,
    pub y: f64,
}

impl Serialize for CurvePoint {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut tuple = serializer.serialize_tuple(2)?;
        tuple.serialize_element(&self.x)?;
        tuple.serialize_element(&self.y)?;
        tuple.end()
    }
}

impl<'de> Deserialize<'de> for CurvePoint {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_any(CurvePointVisitor)
    }
}

struct CurvePointVisitor;

struct F64Value(f64);

impl<'de> Deserialize<'de> for F64Value {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserialize_f64(deserializer).map(Self)
    }
}

impl<'de> Visitor<'de> for CurvePointVisitor {
    type Value = CurvePoint;

    fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("a two-item point tuple")
    }

    fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
    where
        A: SeqAccess<'de>,
    {
        let x = seq
            .next_element::<F64Value>()?
            .ok_or_else(|| de::Error::invalid_length(0, &self))?
            .0;
        let y = seq
            .next_element::<F64Value>()?
            .ok_or_else(|| de::Error::invalid_length(1, &self))?
            .0;
        if seq.next_element::<de::IgnoredAny>()?.is_some() {
            return Err(de::Error::invalid_length(3, &self));
        }
        Ok(CurvePoint { x, y })
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct CubicBezierSpec {
    pub start: CurvePoint,
    pub control_start: CurvePoint,
    pub control_end: CurvePoint,
    pub end: CurvePoint,
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

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct TrimPathsSpec {
    #[serde(deserialize_with = "deserialize_bez_path")]
    path: BezPath,
    outsets: Vec<PathOutsets>,
    #[serde(
        default = "default_arclen_accuracy",
        deserialize_with = "deserialize_f64"
    )]
    accuracy: f64,
    #[serde(default)]
    format: LayerFormat,
    #[serde(default = "default_unit", deserialize_with = "deserialize_f64")]
    unit: f64,
}

#[derive(Debug, Clone, Copy, Default, Deserialize)]
#[serde(rename_all = "kebab-case")]
enum LayerFormat {
    #[default]
    Array,
    Cbor,
}

fn default_unit() -> f64 {
    1.0
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct PathOutsets {
    #[serde(deserialize_with = "deserialize_f64")]
    start_outset: f64,
    #[serde(deserialize_with = "deserialize_f64")]
    end_outset: f64,
}

#[derive(Debug, Serialize)]
struct PathLayerOutput {
    path: CurvePathOutput,
    segments: Vec<CubicBezierSpec>,
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
    #[serde(default, deserialize_with = "deserialize_f64")]
    pub endpoint_slope: f64,
    /// Arc distances along the base path at which to cut the already-fitted pattern.
    #[serde(default, deserialize_with = "deserialize_f64_vec")]
    pub split_at: Vec<f64>,
    #[serde(
        default = "default_arclen_accuracy",
        deserialize_with = "deserialize_f64"
    )]
    pub accuracy: f64,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct StrokeOutlineSpec {
    #[serde(deserialize_with = "deserialize_bez_path")]
    pub path: BezPath,
    #[serde(deserialize_with = "deserialize_f64")]
    pub width: f64,
    #[serde(default = "default_stroke_join")]
    pub join: String,
    #[serde(
        default = "default_stroke_miter_limit",
        deserialize_with = "deserialize_f64"
    )]
    pub miter_limit: f64,
    #[serde(default = "default_stroke_cap")]
    pub start_cap: String,
    #[serde(default = "default_stroke_cap")]
    pub end_cap: String,
    #[serde(
        default = "default_arclen_accuracy",
        deserialize_with = "deserialize_f64"
    )]
    pub accuracy: f64,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
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
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct PathLengthSpec {
    #[serde(deserialize_with = "deserialize_bez_path")]
    pub path: BezPath,
    #[serde(
        default = "default_arclen_accuracy",
        deserialize_with = "deserialize_f64"
    )]
    pub accuracy: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct PathIntersectionsSpec {
    #[serde(
        serialize_with = "serialize_bez_path",
        deserialize_with = "deserialize_bez_path"
    )]
    pub a: BezPath,
    #[serde(
        serialize_with = "serialize_bez_path",
        deserialize_with = "deserialize_bez_path"
    )]
    pub b: BezPath,
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
    FittedCoil(FittedCoilInput),
}

impl PatternInput {
    /// Kurvst's `wave()`: a smooth sine.
    pub fn wave(samples_per_period: usize) -> Self {
        Self::sampled(
            "wave",
            samples_per_period,
            |theta| (0.0, theta.sin()),
            false,
        )
    }

    /// Kurvst's `coil()`: smooth loops, `longitudinal_scale` times as long as wide.
    pub fn coil(samples_per_period: usize, longitudinal_scale: f64) -> Self {
        Self::sampled(
            "coil",
            samples_per_period,
            |theta| (longitudinal_scale * theta.cos(), theta.sin()),
            true,
        )
    }

    /// Kurvst's `zigzag()`: straight segments through the extremes.
    pub fn zigzag() -> Self {
        Self::points(
            "zigzag",
            "linear",
            [(0.0, 0.0), (0.25, 1.0), (0.75, -1.0), (1.0, 0.0)]
                .map(|(at, y)| PatternPointInput {
                    at: Some(at),
                    x: 0.0,
                    y,
                })
                .into(),
            false,
        )
    }

    fn points(
        name: &str,
        interpolation: &str,
        points: Vec<PatternPointInput>,
        endpoint_ramp: bool,
    ) -> Self {
        Self::Points(PointPatternInput {
            kind: "points".to_string(),
            name: Some(name.to_string()),
            points,
            interpolation: interpolation.to_string(),
            endpoint_ramp,
        })
    }

    /// One period sampled at `samples_per_period` equal phase steps.
    fn sampled(
        name: &str,
        samples_per_period: usize,
        mut offset: impl FnMut(f64) -> (f64, f64),
        endpoint_ramp: bool,
    ) -> Self {
        let samples = samples_per_period.max(1);
        let points = (0..=samples)
            .map(|index| {
                let at = index as f64 / samples as f64;
                let (x, y) = offset(std::f64::consts::TAU * at);
                PatternPointInput { at: Some(at), x, y }
            })
            .collect();
        Self::points(name, "smooth", points, endpoint_ramp)
    }
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
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

/// A complete natural coil, fitted before any painting splits are applied.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct FittedCoilInput {
    pub kind: String,
    #[serde(deserialize_with = "deserialize_f64")]
    pub fit_length: f64,
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
    #[serde(
        default = "default_coil_longitudinal_scale",
        deserialize_with = "deserialize_f64"
    )]
    pub longitudinal_scale: f64,
    #[serde(default = "FittedCoilInput::default_samples_per_period")]
    pub samples_per_period: i64,
}

impl FittedCoilInput {
    fn default_samples_per_period() -> i64 {
        default_samples_per_period() as i64
    }

    fn points(&self) -> Result<Vec<PatternPointInput>, String> {
        if self.kind != "fitted-coil" {
            return Err(format!("Unsupported fitted coil kind: {}", self.kind));
        }
        let length = validate_positive(self.fit_length, "coil fit-length")?;
        let amplitude = validate_finite(self.amplitude, "coil amplitude")?;
        let wavelength = validate_positive(self.wavelength, "coil wavelength")?;
        let longitudinal_scale =
            validate_finite(self.longitudinal_scale, "coil longitudinal-scale")?;
        if longitudinal_scale < 0.0 {
            return Err("coil: fitted longitudinal-scale must be non-negative".to_string());
        }
        let span = length + amplitude.abs() * longitudinal_scale * 2.0;
        validate_finite(span, "coil fitted span")?;
        // Removing half a turn preserves the visible-loop count of integer-fitted tapered coils.
        let periods = (length / wavelength).round().max(1.0) - 0.5;
        let samples = (periods * self.samples_per_period.max(1) as f64)
            .ceil()
            .max(2.0);
        if !samples.is_finite() || samples >= usize::MAX as f64 {
            return Err("coil: fitted sample count is too large".to_string());
        }
        let samples = samples as usize;
        let scale = longitudinal_scale * (length / span) * if amplitude < 0.0 { -1.0 } else { 1.0 };
        let mut points = Vec::new();
        points
            .try_reserve_exact(samples + 1)
            .map_err(|_| "coil: fitted sample count is too large")?;
        for index in 0..=samples {
            let at = index as f64 / samples as f64;
            let theta = std::f64::consts::TAU * at;
            // Keep the public Typst callback's operation order, including its recovered parameter.
            let position = theta / std::f64::consts::TAU;
            let (x, y) = if position == 0.0 || position == 1.0 {
                (0.0, 0.0)
            } else {
                let phase = std::f64::consts::PI + periods * theta;
                // Normalize offsets so the ordinary pattern API still applies amplitude once.
                (
                    scale * (1.0 + libm::cos(phase) - 2.0 * position),
                    libm::sin(phase),
                )
            };
            points.push(PatternPointInput { at: Some(at), x, y });
        }
        Ok(points)
    }
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PatternPathOutput {
    #[serde(
        serialize_with = "serialize_bez_path",
        deserialize_with = "deserialize_bez_path"
    )]
    pub path: BezPath,
    pub pattern: String,
    /// `path` cut at sorted, clamped carrier distances, retaining empty intervals.
    pub parts: Vec<CurvePathOutput>,
}

impl PatternPathOutput {
    fn split_at_carrier_distances(
        &mut self,
        distances: &[f64],
        mut split_at: Vec<f64>,
    ) -> Result<(), String> {
        let length = distances.last().copied().unwrap_or(0.0);
        for distance in &mut split_at {
            validate_finite(*distance, "pattern split-at distance")?;
            *distance = distance.clamp(0.0, length);
        }
        split_at.sort_by(f64::total_cmp);
        split_at.insert(0, 0.0);
        split_at.push(length);

        let segments = self.path.segments().collect::<Vec<_>>();
        self.parts = split_at
            .windows(2)
            .map(|interval| {
                let path = BezPath::from_path_segments(
                    segments
                        .iter()
                        .zip(distances.windows(2))
                        .filter_map(|(segment, sample)| {
                            let start = interval[0].max(sample[0]);
                            let end = interval[1].min(sample[1]);
                            if end <= start {
                                return None;
                            }
                            if start == sample[0] && end == sample[1] {
                                return Some(*segment);
                            }
                            // The full pattern has already been sampled and fitted. Map
                            // carrier distance onto that segment's parameter without
                            // adding knots or measuring the longer decorated curve.
                            let span = sample[1] - sample[0];
                            let t0 = (start - sample[0]) / span;
                            let t1 = (end - sample[0]) / span;
                            Some(segment.subsegment(t0..t1))
                        }),
                );
                CurvePathOutput { path }
            })
            .collect();
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CurvePathOutput {
    #[serde(
        serialize_with = "serialize_bez_path",
        deserialize_with = "deserialize_bez_path"
    )]
    pub path: BezPath,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub struct PathIntersection {
    pub point: CurvePoint,
    pub distance_a: f64,
    pub distance_b: f64,
    pub segment_a: usize,
    pub segment_b: usize,
    pub t_a: f64,
    pub t_b: f64,
}

#[derive(Debug, Clone, Copy)]
struct SegmentInfo {
    segment: PathSeg,
    offset: f64,
    length: f64,
}

#[derive(Debug, Clone, Copy)]
struct FlatPiece {
    line: Line,
    t0: f64,
    t1: f64,
}

const MAX_INTERSECTION_FLATTEN_DEPTH: usize = 20;
const TANGENT_SINE_TOLERANCE: f64 = 1e-7;

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

fn default_stroke_join() -> String {
    "miter".to_string()
}

fn default_stroke_cap() -> String {
    "butt".to_string()
}

fn default_stroke_miter_limit() -> f64 {
    4.0
}

fn encode_cbor<T: Serialize>(value: &T) -> Result<Vec<u8>, String> {
    let mut bytes = Vec::new();
    ciborium::ser::into_writer(value, &mut bytes)
        .map_err(|err| format!("Failed to serialize CBOR value: {err}"))?;
    Ok(bytes)
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct WirePath {
    elements: Vec<WirePathElement>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "kebab-case", deny_unknown_fields)]
enum WirePathElement {
    Move {
        start: CurvePoint,
    },
    Line {
        end: CurvePoint,
    },
    Quad {
        control: CurvePoint,
        end: CurvePoint,
    },
    Cubic {
        #[serde(rename = "control-start")]
        control_start: CurvePoint,
        #[serde(rename = "control-end")]
        control_end: CurvePoint,
        end: CurvePoint,
    },
    // A struct variant makes Serde reject unsupported fields on close elements.
    Close {},
}

impl From<&BezPath> for WirePath {
    fn from(path: &BezPath) -> Self {
        let elements = path
            .elements()
            .iter()
            .map(|element| match *element {
                PathEl::MoveTo(point) => WirePathElement::Move {
                    start: point.into(),
                },
                PathEl::LineTo(point) => WirePathElement::Line { end: point.into() },
                PathEl::QuadTo(control, end) => WirePathElement::Quad {
                    control: control.into(),
                    end: end.into(),
                },
                PathEl::CurveTo(control_start, control_end, end) => WirePathElement::Cubic {
                    control_start: control_start.into(),
                    control_end: control_end.into(),
                    end: end.into(),
                },
                PathEl::ClosePath => WirePathElement::Close {},
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
                WirePathElement::Move { start } => path.move_to(Point::from(start)),
                WirePathElement::Line { end } => path.line_to(Point::from(end)),
                WirePathElement::Quad { control, end } => {
                    path.quad_to(Point::from(control), Point::from(end));
                }
                WirePathElement::Cubic {
                    control_start,
                    control_end,
                    end,
                } => {
                    path.curve_to(
                        Point::from(control_start),
                        Point::from(control_end),
                        Point::from(end),
                    );
                }
                WirePathElement::Close {} => path.close_path(),
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

pub(crate) fn deserialize_f64<'de, D>(deserializer: D) -> Result<f64, D::Error>
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

fn deserialize_f64_vec<'de, D>(deserializer: D) -> Result<Vec<f64>, D::Error>
where
    D: Deserializer<'de>,
{
    Vec::<F64Value>::deserialize(deserializer)
        .map(|values| values.into_iter().map(|value| value.0).collect())
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
        CubicBez::new(
            value.start,
            value.control_start,
            value.control_end,
            value.end,
        )
    }
}

impl From<CubicBez> for CubicBezierSpec {
    fn from(value: CubicBez) -> Self {
        Self {
            start: value.p0.into(),
            control_start: value.p1.into(),
            control_end: value.p2.into(),
            end: value.p3.into(),
        }
    }
}

impl From<CubicBezierSpec> for PathSeg {
    fn from(value: CubicBezierSpec) -> Self {
        PathSeg::Cubic(value.into())
    }
}

pub fn curve_trim_path_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: TrimPathSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize path trim spec: {err}"))?;
    let output = spec.trimmed()?;
    encode_cbor(&output)
}

pub fn curve_trim_paths_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: TrimPathsSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize path trim batch: {err}"))?;
    let (format, unit) = (spec.format, spec.unit);
    let layers = spec.trimmed()?;
    match format {
        LayerFormat::Array => encode_cbor(&layers),
        LayerFormat::Cbor => encode_cbor(&PackedPathLayers::new(&layers, unit)?),
    }
}

pub fn curve_hobby_through_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: HobbyThroughSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize Hobby curve spec: {err}"))?;
    let output = spec.curve()?;
    encode_cbor(&output)
}

pub fn curve_hobby_spline_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: HobbySplineSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize Hobby spline spec: {err}"))?;
    let output = spec.curve()?;
    encode_cbor(&output)
}

pub fn curve_fitted_coil_points_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: FittedCoilInput = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize fitted coil spec: {err}"))?;
    encode_cbor(&spec.points()?)
}

pub fn curve_pattern_path_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: PatternPathSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize path pattern spec: {err}"))?;
    let output = spec.patterned()?;
    encode_cbor(&output)
}

#[derive(Deserialize)]
struct PatternCetzSpec {
    pattern: PatternPathSpec,
    #[serde(deserialize_with = "deserialize_f64")]
    unit: f64,
}

#[derive(Serialize)]
struct CetzSubpath<P = [f64; 3]>(P, bool, Vec<CetzSegment<P>>);

#[derive(Serialize)]
#[serde(untagged)]
enum CetzSegment<P> {
    Line((&'static str, P)),
    Cubic((&'static str, P, P, P)),
}

#[derive(Serialize)]
#[serde(untagged)]
enum PatternCetzOutput {
    Data(Vec<CetzSubpath>),
    Original(PatternPathOutput),
}

impl<P> CetzSubpath<P> {
    // Match Typst's _to-cetz-data numeric representation, including subpath
    // boundaries, quadratic promotion, and multiplication before float coercion.
    fn from_path(path: &BezPath, point: impl Fn(Point) -> P) -> Option<Vec<Self>> {
        let mut subpaths = Vec::new();
        let mut origin = None;
        let mut current = Point::ZERO;
        let mut segments = Vec::new();
        for &element in path.elements() {
            if let PathEl::MoveTo(start) = element {
                if let Some(previous) = origin {
                    subpaths.push(CetzSubpath(
                        point(previous),
                        false,
                        std::mem::take(&mut segments),
                    ));
                }
                origin = Some(start);
                current = start;
                continue;
            }
            // Typst synthesizes an integer (0, 0) origin when absent. Retain
            // those number types and signed-zero behavior for coordinate hooks.
            origin?;
            match element {
                PathEl::MoveTo(_) => unreachable!(),
                PathEl::LineTo(end) => {
                    segments.push(CetzSegment::Line(("l", point(end))));
                    current = end;
                }
                PathEl::QuadTo(control, end) => {
                    // Keep the scalar operation order of _quad-cubic-segment;
                    // Kurbo's general quadratic conversion uses another order.
                    let first = Point::new(
                        current.x + (control.x - current.x) * 2.0 / 3.0,
                        current.y + (control.y - current.y) * 2.0 / 3.0,
                    );
                    let second = Point::new(
                        end.x + (control.x - end.x) * 2.0 / 3.0,
                        end.y + (control.y - end.y) * 2.0 / 3.0,
                    );
                    segments.push(CetzSegment::Cubic((
                        "c",
                        point(first),
                        point(second),
                        point(end),
                    )));
                    current = end;
                }
                PathEl::CurveTo(first, second, end) => {
                    segments.push(CetzSegment::Cubic((
                        "c",
                        point(first),
                        point(second),
                        point(end),
                    )));
                    current = end;
                }
                PathEl::ClosePath => {
                    let start = origin.take().unwrap();
                    subpaths.push(CetzSubpath(
                        point(start),
                        true,
                        std::mem::take(&mut segments),
                    ));
                    current = start;
                }
            }
        }
        if let Some(start) = origin {
            subpaths.push(CetzSubpath(point(start), false, segments));
        }
        Some(subpaths)
    }
}

impl PatternCetzOutput {
    fn from_pattern(output: PatternPathOutput, unit: f64) -> Self {
        match CetzSubpath::from_path(&output.path, |p| [p.x * unit, p.y * unit, 0.0]) {
            Some(subpaths) => Self::Data(subpaths),
            None => Self::Original(output),
        }
    }
}

#[derive(Serialize)]
struct FootprintPath {
    segments: Vec<CetzSubpath<[f64; 2]>>,
    single: bool,
}

#[derive(Serialize)]
#[serde(rename_all = "kebab-case")]
struct PackedPathLayers {
    #[serde(serialize_with = "PackedPathLayers::serialize_bytes")]
    layers: Vec<u8>,
    #[serde(serialize_with = "PackedPathLayers::serialize_bytes")]
    footprints: Vec<u8>,
    count: usize,
    nonzero: bool,
    supported: bool,
    all_single: bool,
}

impl PackedPathLayers {
    fn serialize_bytes<S: Serializer>(bytes: &[u8], serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_bytes(bytes)
    }

    fn new(layers: &[PathLayerOutput], unit: f64) -> Result<Self, String> {
        let mut supported = true;
        let mut nonzero = false;
        let mut all_single = !layers.is_empty();
        let mut footprints = Vec::with_capacity(layers.len());
        for layer in layers {
            nonzero |= layer.segments.iter().any(|segment| {
                segment.start != segment.control_start
                    || segment.start != segment.control_end
                    || segment.start != segment.end
            });
            let single = layer.segments.len() == 1;
            all_single &= single;
            // The drawing owner elevates one segment to a Bezier before unit
            // conversion; longer carriers retain their original path commands.
            let segments = if single {
                let cubic = layer.segments[0];
                let mut path = BezPath::new();
                path.move_to(cubic.start);
                path.curve_to(cubic.control_start, cubic.control_end, cubic.end);
                CetzSubpath::from_path(&path, |p| [p.x * 1.0, p.y * 1.0])
            } else if unit.is_finite() && unit != 0.0 {
                CetzSubpath::from_path(&layer.path.path, |p| [p.x * unit, p.y * unit])
            } else {
                None
            };
            let segments = match segments {
                Some(segments) => segments,
                None => {
                    supported = false;
                    Vec::new()
                }
            };
            // Move-only subpaths require Kurvst's primitive drawing semantics.
            supported &= !segments.is_empty() && segments.iter().all(|path| !path.2.is_empty());
            footprints.push(FootprintPath { segments, single });
        }
        Ok(Self {
            layers: encode_cbor(&layers)?,
            footprints: encode_cbor(&footprints)?,
            count: layers.len(),
            nonzero,
            supported,
            all_single,
        })
    }
}

pub fn curve_pattern_cetz_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: PatternCetzSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize CeTZ pattern spec: {err}"))?;
    let output = spec.pattern.patterned()?;
    encode_cbor(&PatternCetzOutput::from_pattern(output, spec.unit))
}

pub fn curve_parallel_path_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: ParallelPathSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize parallel path spec: {err}"))?;
    let output = spec.parallel()?;
    encode_cbor(&output)
}

pub fn curve_stroke_outline_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: StrokeOutlineSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize stroke outline spec: {err}"))?;
    let output = stroke_outline(spec)?;
    encode_cbor(&output)
}

pub fn curve_path_length_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: PathLengthSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize path length spec: {err}"))?;
    let output = spec.length()?;
    encode_cbor(&output)
}

#[derive(Deserialize)]
pub struct PathFramesSpec {
    #[serde(deserialize_with = "deserialize_bez_path")]
    pub path: BezPath,
    #[serde(deserialize_with = "deserialize_f64_vec")]
    pub distances: Vec<f64>,
    #[serde(
        default = "default_arclen_accuracy",
        deserialize_with = "deserialize_f64"
    )]
    pub accuracy: f64,
}

#[derive(Debug, Serialize, PartialEq)]
pub struct PathFrame {
    pub point: CurvePoint,
    pub tangent: CurvePoint,
}

impl CubicBezierSpec {
    // Match Typst's segments() conversions and de Casteljau arithmetic exactly:
    // equivalent rearrangements can perturb the label optimizer's tie-breaking.
    fn drawable(segment: PathSeg) -> Self {
        let lerp_third = |a: Point, b: Point, twice: bool| CurvePoint {
            x: a.x
                + if twice {
                    (b.x - a.x) * 2.0 / 3.0
                } else {
                    (b.x - a.x) / 3.0
                },
            y: a.y
                + if twice {
                    (b.y - a.y) * 2.0 / 3.0
                } else {
                    (b.y - a.y) / 3.0
                },
        };
        match segment {
            PathSeg::Cubic(cubic) => cubic.into(),
            PathSeg::Line(line) => Self {
                start: line.p0.into(),
                end: line.p1.into(),
                control_start: lerp_third(line.p0, line.p1, false),
                control_end: lerp_third(line.p0, line.p1, true),
            },
            PathSeg::Quad(quad) => Self {
                start: quad.p0.into(),
                end: quad.p2.into(),
                control_start: lerp_third(quad.p0, quad.p1, true),
                control_end: lerp_third(quad.p2, quad.p1, true),
            },
        }
    }

    fn endpoint_frame(self, at_end: bool) -> PathFrame {
        let t = if at_end { 1.0 } else { 0.0 };
        let lerp = |a: CurvePoint, b: CurvePoint| CurvePoint {
            x: a.x + (b.x - a.x) * t,
            y: a.y + (b.y - a.y) * t,
        };
        let ab = lerp(self.start, self.control_start);
        let bc = lerp(self.control_start, self.control_end);
        let cd = lerp(self.control_end, self.end);
        let abc = lerp(ab, bc);
        let bcd = lerp(bc, cd);
        PathFrame {
            point: if at_end { self.end } else { self.start },
            tangent: CurvePoint {
                x: bcd.x - abc.x,
                y: bcd.y - abc.y,
            },
        }
    }
}

impl PathFramesSpec {
    /// Point and tangent at each arc distance, clamped onto the path.
    pub fn frames(self) -> Result<Vec<Option<PathFrame>>, String> {
        let accuracy = validate_positive_accuracy(self.accuracy)?;
        let trimmer = PathTrimmer::new(self.path.segments(), accuracy);
        self.distances
            .into_iter()
            .map(|distance| trimmer.frame(distance))
            .collect()
    }
}

pub fn curve_path_frames_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: PathFramesSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize path frames spec: {err}"))?;
    encode_cbor(&spec.frames()?)
}

pub fn curve_path_intersections_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: PathIntersectionsSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize path intersections spec: {err}"))?;
    let output = path_intersections(spec)?;
    encode_cbor(&output)
}

impl TrimPathSpec {
    /// Trim the path by arc-length outsets at both ends.
    pub fn trimmed(self) -> Result<CurvePathOutput, String> {
        let accuracy = validate_positive_accuracy(self.accuracy)?;
        let start_outset = validate_outset(self.start_outset, "start")?;
        let end_outset = validate_outset(self.end_outset, "end")?;
        let segments =
            PathTrimmer::new(self.path.segments(), accuracy).trim(start_outset, end_outset)?;
        curve_path_from_segments(segments)
    }
}

impl TrimPathsSpec {
    /// Trim one measured path to every outset window.
    pub fn trimmed(self) -> Result<Vec<PathLayerOutput>, String> {
        let accuracy = validate_positive_accuracy(self.accuracy)?;
        let trimmer = PathTrimmer::new(self.path.segments(), accuracy);
        self.outsets
            .into_iter()
            .map(|outsets| {
                let start = validate_outset(outsets.start_outset, "start")?;
                let end = validate_outset(outsets.end_outset, "end")?;
                let path = curve_path_from_segments(trimmer.trim(start, end)?)?;
                let segments = path
                    .path
                    .segments()
                    .map(CubicBezierSpec::drawable)
                    .collect();
                Ok(PathLayerOutput { path, segments })
            })
            .collect()
    }
}

impl PathLengthSpec {
    /// Arc length of the path.
    pub fn length(self) -> Result<f64, String> {
        let accuracy = validate_positive_accuracy(self.accuracy)?;
        Ok(self
            .path
            .segments()
            .map(|segment| segment.arclen(accuracy))
            .sum())
    }
}

fn path_intersections(spec: PathIntersectionsSpec) -> Result<Vec<PathIntersection>, String> {
    let accuracy = validate_positive_accuracy(spec.accuracy)?;
    let a = intersection_segment_info(&spec.a, accuracy)?;
    let b = intersection_segment_info(&spec.b, accuracy)?;
    let total_a = a.last().map_or(0.0, |item| item.offset + item.length);
    let total_b = b.last().map_or(0.0, |item| item.offset + item.length);
    let b_pieces = b
        .iter()
        .map(|info| {
            let mut pieces = Vec::new();
            flatten_intersection_segment(info.segment, 0.0, 1.0, accuracy * 0.25, 0, &mut pieces)?;
            Ok(pieces)
        })
        .collect::<Result<Vec<_>, String>>()?;
    let mut result = Vec::new();

    for (segment_a, a_info) in a.iter().enumerate() {
        for (segment_b, (b_info, pieces)) in b.iter().zip(&b_pieces).enumerate() {
            for &piece in pieces {
                if (piece.line.p1 - piece.line.p0).hypot() <= f64::EPSILON {
                    continue;
                }
                for hit in a_info.segment.intersect_line(piece.line) {
                    let t_a = hit.segment_t.clamp(0.0, 1.0);
                    let t_b = piece.t0 + hit.line_t * (piece.t1 - piece.t0);
                    let Some((t_a, t_b, point)) =
                        refine_crossing(a_info.segment, b_info.segment, t_a, t_b, accuracy)
                    else {
                        continue;
                    };
                    let distance_a =
                        a_info.offset + a_info.segment.subsegment(0.0..t_a).arclen(accuracy);
                    let distance_b =
                        b_info.offset + b_info.segment.subsegment(0.0..t_b).arclen(accuracy);
                    if path_endpoint_contact(distance_a, total_a, distance_b, total_b, accuracy) {
                        continue;
                    }
                    result.push(PathIntersection {
                        point: point.into(),
                        distance_a,
                        distance_b,
                        segment_a,
                        segment_b,
                        t_a,
                        t_b,
                    });
                }
            }
        }
    }

    result.sort_by(|left, right| {
        left.distance_a
            .total_cmp(&right.distance_a)
            .then(left.distance_b.total_cmp(&right.distance_b))
    });
    result.dedup_by(|left, right| same_crossing(*left, *right, accuracy));
    Ok(result)
}

fn intersection_segment_info(path: &BezPath, accuracy: f64) -> Result<Vec<SegmentInfo>, String> {
    let mut offset = 0.0;
    path.segments()
        .map(|segment| {
            if !segment.is_finite() {
                return Err("curve intersection paths must contain only finite points".to_string());
            }
            let length = segment.arclen(accuracy);
            if !length.is_finite() {
                return Err("curve intersection path length must be finite".to_string());
            }
            let info = SegmentInfo {
                segment,
                offset,
                length,
            };
            offset += length;
            if !offset.is_finite() {
                return Err("curve intersection cumulative path length must be finite".to_string());
            }
            Ok(info)
        })
        .collect()
}

fn flatten_intersection_segment(
    segment: PathSeg,
    t0: f64,
    t1: f64,
    tolerance: f64,
    depth: usize,
    output: &mut Vec<FlatPiece>,
) -> Result<(), String> {
    if matches!(segment, PathSeg::Line(_)) || intersection_flatness(segment) <= tolerance {
        output.push(FlatPiece {
            line: Line::new(segment.start(), segment.end()),
            t0,
            t1,
        });
        return Ok(());
    }
    if depth == MAX_INTERSECTION_FLATTEN_DEPTH {
        return Err("curve intersection accuracy is too small for the path scale".to_string());
    }

    let middle = (t0 + t1) * 0.5;
    flatten_intersection_segment(
        segment.subsegment(0.0..0.5),
        t0,
        middle,
        tolerance,
        depth + 1,
        output,
    )?;
    flatten_intersection_segment(
        segment.subsegment(0.5..1.0),
        middle,
        t1,
        tolerance,
        depth + 1,
        output,
    )
}

fn intersection_flatness(segment: PathSeg) -> f64 {
    match segment {
        PathSeg::Line(_) => 0.0,
        PathSeg::Quad(quad) => point_chord_distance(quad.p1, quad.p0, quad.p2)
            .max(control_polygon_excess(&[quad.p0, quad.p1, quad.p2])),
        PathSeg::Cubic(cubic) => point_chord_distance(cubic.p1, cubic.p0, cubic.p3)
            .max(point_chord_distance(cubic.p2, cubic.p0, cubic.p3))
            .max(control_polygon_excess(&[
                cubic.p0, cubic.p1, cubic.p2, cubic.p3,
            ])),
    }
}

fn point_chord_distance(point: Point, start: Point, end: Point) -> f64 {
    let chord = end - start;
    let length = chord.hypot();
    if length <= f64::EPSILON {
        (point - start).hypot()
    } else {
        point_line_distance(point, start, chord, length)
    }
}

fn control_polygon_excess(points: &[Point]) -> f64 {
    let polygon: f64 = points
        .windows(2)
        .map(|pair| (pair[1] - pair[0]).hypot())
        .sum();
    (polygon - (points[points.len() - 1] - points[0]).hypot()).max(0.0)
}

fn refine_crossing(
    a: PathSeg,
    b: PathSeg,
    mut t_a: f64,
    mut t_b: f64,
    accuracy: f64,
) -> Option<(f64, f64, Point)> {
    for _ in 0..16 {
        let point_a = a.eval(t_a);
        let point_b = b.eval(t_b);
        let residual = point_a - point_b;
        if residual.hypot() <= accuracy * 0.05 {
            break;
        }
        let tangent_a = path_seg_tangent(&a, t_a);
        let tangent_b = path_seg_tangent(&b, t_b);
        let determinant = tangent_a.cross(tangent_b);
        if determinant.abs() <= f64::EPSILON * tangent_a.hypot() * tangent_b.hypot() {
            return None;
        }
        t_a -= residual.cross(tangent_b) / determinant;
        t_b += tangent_a.cross(residual) / determinant;
        if !(-1e-7..=1.0 + 1e-7).contains(&t_a) || !(-1e-7..=1.0 + 1e-7).contains(&t_b) {
            return None;
        }
        t_a = t_a.clamp(0.0, 1.0);
        t_b = t_b.clamp(0.0, 1.0);
    }

    let point_a = a.eval(t_a);
    let point_b = b.eval(t_b);
    if (point_a - point_b).hypot() > accuracy {
        return None;
    }
    let tangent_a = path_seg_tangent(&a, t_a);
    let tangent_b = path_seg_tangent(&b, t_b);
    let tangent_scale = tangent_a.hypot() * tangent_b.hypot();
    if tangent_scale <= f64::EPSILON
        || tangent_a.cross(tangent_b).abs() / tangent_scale <= TANGENT_SINE_TOLERANCE
    {
        return None;
    }
    Some((
        t_a,
        t_b,
        Point::new((point_a.x + point_b.x) * 0.5, (point_a.y + point_b.y) * 0.5),
    ))
}

fn path_endpoint_contact(
    distance_a: f64,
    total_a: f64,
    distance_b: f64,
    total_b: f64,
    accuracy: f64,
) -> bool {
    let endpoint = |distance: f64, total: f64| distance <= accuracy || total - distance <= accuracy;
    endpoint(distance_a, total_a) || endpoint(distance_b, total_b)
}

fn same_crossing(a: PathIntersection, b: PathIntersection, accuracy: f64) -> bool {
    let a_point = Point::from(a.point);
    let b_point = Point::from(b.point);
    (a_point - b_point).hypot() <= accuracy
        && (a.distance_a - b.distance_a).abs() <= accuracy * 2.0
        && (a.distance_b - b.distance_b).abs() <= accuracy * 2.0
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

impl PatternPathSpec {
    /// Decorate the path with its pattern, split at carrier distances.
    pub fn patterned(self) -> Result<PatternPathOutput, String> {
        let amplitude = validate_finite(self.amplitude, "pattern amplitude")?;
        let wavelength = validate_positive(self.wavelength, "pattern wavelength")?;
        let phase = validate_finite(self.phase, "pattern phase")?;
        validate_finite(self.coil_longitudinal_scale, "coil longitudinal scale")?;
        let endpoint_slope = validate_finite(self.endpoint_slope, "pattern endpoint-slope")?;
        if !(0.0..=3.0).contains(&endpoint_slope) {
            return Err("pattern endpoint-slope must be between 0 and 3".to_string());
        }
        let accuracy = validate_positive_accuracy(self.accuracy)?;
        if self.samples_per_period == 0 {
            return Err("pattern samples-per-period must be positive".to_string());
        }
        let fitted_coil = matches!(&self.pattern, PatternInput::FittedCoil(_));
        let pattern = PointPattern::from_input(self.pattern)?;
        let samples_per_period = if fitted_coil {
            pattern.points.len() - 1
        } else {
            self.samples_per_period
        };

        let path_segments = self.path.segments().collect::<Vec<_>>();
        let segment_lengths = path_segments
            .iter()
            .map(|segment| segment.arclen(accuracy))
            .collect::<Vec<_>>();
        let length: f64 = segment_lengths.iter().sum();
        if length <= f64::EPSILON {
            let mut output = PatternPathOutput {
                path: self.path,
                pattern: pattern.name.clone(),
                parts: Vec::new(),
            };
            output.split_at_carrier_distances(&[0.0], self.split_at)?;
            return Ok(output);
        }

        let distances = pattern.distances(length, wavelength, samples_per_period);
        let mut points = Vec::with_capacity(distances.len());
        for &distance in &distances {
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
            let pattern_point =
                pattern.evaluate(distance / wavelength + phase / std::f64::consts::TAU);
            let envelope = pattern_endpoint_envelope(
                &pattern,
                distance,
                length,
                wavelength,
                self.anchor_start,
                self.anchor_end,
                endpoint_slope,
            );
            // Delay the doubling-back motion until the coil has opened sideways.
            let longitudinal = amplitude * envelope * envelope * pattern_point.x;
            let lateral = amplitude * envelope * pattern_point.y;
            points.push((base + tangent * longitudinal + normal * lateral).into());
        }
        if !pattern.endpoint_ramp && self.anchor_start {
            points[0] = path_segments
                .first()
                .map(PathSeg::start)
                .unwrap_or(Point::new(0.0, 0.0))
                .into();
        }
        if !pattern.endpoint_ramp && self.anchor_end {
            let last_index = points.len() - 1;
            points[last_index] = path_segments
                .last()
                .map(PathSeg::end)
                .unwrap_or(Point::new(0.0, 0.0))
                .into();
        }

        let path = match pattern.interpolation {
            PatternInterpolation::Linear => BezPath::from_path_segments(
                points
                    .windows(2)
                    .map(|window| PathSeg::Line(Line::new(window[0], window[1]))),
            ),
            PatternInterpolation::Smooth => BezPath::from_path_segments(
                cubic_spline_through_points(&points)
                    .into_iter()
                    .map(PathSeg::from),
            ),
        };

        let mut output = PatternPathOutput {
            path,
            pattern: pattern.name,
            parts: Vec::new(),
        };
        output.split_at_carrier_distances(&distances, self.split_at)?;
        Ok(output)
    }
}

impl ParallelPathSpec {
    /// Offset each source subpath (positive to the left), then trim by arc length.
    /// Differing fitted endpoints are joined by straight bevels, never snapped.
    /// Explicit source moves remain boundaries, even at identical coordinates.
    pub fn parallel(self) -> Result<CurvePathOutput, String> {
        let distance = validate_finite(self.distance, "parallel path distance")?;
        let accuracy = validate_positive_accuracy(self.accuracy)?;
        let start_outset = validate_outset(self.start_outset, "start")?;
        let end_outset = validate_outset(self.end_outset, "end")?;
        let source_length: f64 = self
            .path
            .segments()
            .map(|segment| segment.arclen(accuracy))
            .sum();
        if source_length <= f64::EPSILON {
            return Ok(CurvePathOutput { path: self.path });
        }

        let mut subpaths = Vec::new();
        for elements in self.path.subpaths() {
            let source = BezPath::from_vec(elements.to_vec());
            let mut segments: Vec<PathSeg> = Vec::new();
            for segment in source.segments() {
                // A stationary segment has no offset normal and must not
                // introduce a spurious detour between meaningful pieces.
                if segment.arclen(accuracy) == 0.0 {
                    continue;
                }
                for fitted in offset_path_segment(segment, distance, accuracy) {
                    if let Some(previous) = segments.last()
                        && previous.end() != fitted.start()
                    {
                        segments.push(PathSeg::Line(Line::new(previous.end(), fitted.start())));
                    }
                    segments.push(fitted);
                }
            }
            let closed = elements.last() == Some(&PathEl::ClosePath);
            if closed
                && let (Some(first), Some(last)) = (segments.first(), segments.last())
                && last.end() != first.start()
            {
                segments.push(PathSeg::Line(Line::new(last.end(), first.start())));
            }
            subpaths.push((source, closed, PathTrimmer::new(segments, accuracy)));
        }

        // Trim the joined offset, not the original carrier or independently
        // fitted pieces. Moves consume no arc length and are never bridged.
        let length = subpaths.iter().map(|(_, _, trim)| trim.length()).sum();
        let (start, end) = fit_outsets_to_length(start_outset, end_outset, length);
        let mut cursor = 0.0;
        let mut path = BezPath::new();
        for (source, closed, trimmer) in subpaths {
            let subpath_length = trimmer.length();
            let local_start = (start - cursor).max(0.0);
            let local_end = (cursor + subpath_length - (length - end)).max(0.0);
            cursor += subpath_length;
            if subpath_length == 0.0 {
                if start == 0.0 && end == 0.0 {
                    path.extend(source.elements().iter().copied());
                }
                continue;
            }
            if local_start + local_end >= subpath_length {
                continue;
            }
            let segments = trimmer.trim(local_start, local_end)?;
            let fitted = BezPath::from_path_segments(segments.into_iter());
            path.extend(fitted.elements().iter().copied());
            if closed && local_start == 0.0 && local_end == 0.0 {
                path.close_path();
            }
        }
        if path.segments().next().is_none() {
            return Err("path fitting produced no visible path".to_string());
        }
        Ok(CurvePathOutput { path })
    }
}

/// Expand a stroked path into a closed fill outline with joins and caps.
///
/// Open subpaths become one closed contour; closed subpaths become an outer and
/// an inner contour of opposite winding, so fill the result with the non-zero rule.
fn stroke_outline(spec: StrokeOutlineSpec) -> Result<CurvePathOutput, String> {
    let width = validate_positive(spec.width, "stroke outline width")?;
    let accuracy = validate_positive_accuracy(spec.accuracy)?;
    let miter_limit = validate_positive(spec.miter_limit, "stroke outline miter limit")?;
    let style = Stroke {
        join: parse_stroke_join(&spec.join)?,
        miter_limit,
        start_cap: parse_stroke_cap(&spec.start_cap)?,
        end_cap: parse_stroke_cap(&spec.end_cap)?,
        ..Stroke::new(width)
    };
    let path = kurbo::stroke(spec.path, &style, &StrokeOpts::default(), accuracy);
    Ok(CurvePathOutput { path })
}

pub(crate) fn parse_stroke_join(value: &str) -> Result<Join, String> {
    match value.trim().to_ascii_lowercase().as_str() {
        "miter" => Ok(Join::Miter),
        "round" => Ok(Join::Round),
        "bevel" => Ok(Join::Bevel),
        other => Err(format!("Unsupported stroke join: {other}")),
    }
}

pub(crate) fn parse_stroke_cap(value: &str) -> Result<Cap, String> {
    match value.trim().to_ascii_lowercase().as_str() {
        "butt" => Ok(Cap::Butt),
        "square" => Ok(Cap::Square),
        "round" => Ok(Cap::Round),
        other => Err(format!("Unsupported stroke cap: {other}")),
    }
}

fn offset_path_segment(segment: PathSeg, distance: f64, accuracy: f64) -> Vec<PathSeg> {
    match segment {
        PathSeg::Line(line) => vec![offset_line(line.p0, line.p1, distance)],
        PathSeg::Quad(_) | PathSeg::Cubic(_) => {
            let cubic = segment.to_cubic();
            if cubic_is_collinear(cubic, accuracy) {
                vec![offset_line(cubic.p0, cubic.p3, distance)]
            } else {
                let mut path = BezPath::new();
                offset_cubic(cubic, distance, accuracy, &mut path);
                path.segments().collect()
            }
        }
    }
}

fn offset_line(start: Point, end: Point, distance: f64) -> PathSeg {
    let tangent = normalized_tangent_between(end - start, start, end);
    let normal = tangent.turn_90() * distance;
    PathSeg::Line(Line::new(start + normal, end + normal))
}

fn cubic_is_collinear(cubic: CubicBez, accuracy: f64) -> bool {
    let chord = cubic.p3 - cubic.p0;
    let chord_len = chord.hypot();
    if chord_len <= f64::EPSILON {
        return true;
    }
    let tolerance = accuracy.max(f64::EPSILON) * chord_len;
    point_line_distance(cubic.p1, cubic.p0, chord, chord_len) <= tolerance
        && point_line_distance(cubic.p2, cubic.p0, chord, chord_len) <= tolerance
}

fn point_line_distance(point: Point, origin: Point, direction: Vec2, direction_len: f64) -> f64 {
    let delta = point - origin;
    (direction.x * delta.y - direction.y * delta.x).abs() / direction_len
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
            PatternInput::FittedCoil(input) => Self::from_points(PointPatternInput {
                kind: "points".to_string(),
                name: Some("coil".to_string()),
                points: input.points()?,
                interpolation: "smooth".to_string(),
                endpoint_ramp: false,
            }),
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
    endpoint_slope: f64,
) -> f64 {
    if !pattern.endpoint_ramp {
        return 1.0;
    }

    // Ease in over three quarters of a turn; the longitudinal offset ramps more slowly.
    let ramp = (wavelength * 0.75).min(length * 0.5);
    if ramp <= f64::EPSILON {
        return 1.0;
    }

    let mut envelope: f64 = 1.0;
    if anchor_start {
        let t = (distance / ramp).clamp(0.0, 1.0);
        envelope = envelope.min(smoothstep(t) + endpoint_slope * t * (1.0 - t).powi(2));
    }
    if anchor_end {
        let t = ((length - distance) / ramp).clamp(0.0, 1.0);
        envelope = envelope.min(smoothstep(t) + endpoint_slope * t * (1.0 - t).powi(2));
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
    // Coalesce terminal roundoff at the path's scale, keeping the exact endpoint below.
    let end_tolerance = 4.0 * f64::EPSILON * length;
    for i in 0..count {
        let distance = i as f64 * step;
        if length - distance <= end_tolerance {
            break;
        }
        distances.push(distance);
    }
    distances.push(length);
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

fn curve_path_from_segments(
    segments: impl IntoIterator<Item = PathSeg>,
) -> Result<CurvePathOutput, String> {
    let path = BezPath::from_path_segments(segments.into_iter());
    if path.segments().next().is_none() {
        return Err("path fitting produced no visible path".to_string());
    }

    Ok(CurvePathOutput { path })
}

// One measured carrier can serve many trim windows. Scalar trims use this
// same owner; untouched paths still avoid arc-length measurement entirely.
pub struct PathTrimmer {
    segments: Vec<PathSeg>,
    accuracy: f64,
    lengths: OnceCell<(Vec<f64>, f64)>,
}

impl PathTrimmer {
    pub fn new(segments: impl IntoIterator<Item = PathSeg>, accuracy: f64) -> Self {
        Self {
            segments: segments.into_iter().collect(),
            accuracy,
            lengths: OnceCell::new(),
        }
    }

    fn measured_lengths(&self) -> &(Vec<f64>, f64) {
        self.lengths.get_or_init(|| {
            let lengths: Vec<_> = self
                .segments
                .iter()
                .map(|segment| segment.arclen(self.accuracy))
                .collect();
            let length = lengths.iter().sum();
            (lengths, length)
        })
    }

    /// Reuse the same lazy arc table for all placements and trim windows.
    pub fn length(&self) -> f64 {
        self.measured_lengths().1
    }

    /// A monotone, exactly collinear Bezier is a line geometrically, even when
    /// its parameter speed is nonuniform. Solve its scalar polynomial rather
    /// than introducing arc-inversion tolerance into straight-carrier points.
    /// Keep the original segment so trims and frame tangents retain its controls.
    fn parameter_at_distance(&self, segment: PathSeg, distance: f64) -> f64 {
        let chord = segment.end() - segment.start();
        let length = chord.hypot();
        let controls = match segment {
            PathSeg::Line(_) => return segment.inv_arclen(distance, self.accuracy),
            PathSeg::Quad(quad) => [quad.p1, quad.p2],
            PathSeg::Cubic(cubic) => [cubic.p1, cubic.p2],
        };
        let axis = if chord.x.abs() >= chord.y.abs() {
            chord.x
        } else {
            chord.y
        };
        let coordinate = |point: Point| {
            let offset = point - segment.start();
            if chord.x.abs() >= chord.y.abs() {
                offset.x / chord.x
            } else {
                offset.y / chord.y
            }
        };
        // Exact collinearity is intentional: near-straight curves and
        // collinear backtracking controls must retain ordinary arc inversion.
        if length == 0.0
            || !length.is_finite()
            || !axis.is_normal()
            || controls.iter().any(|point| {
                let offset = *point - segment.start();
                let left = offset.x * chord.y;
                let right = offset.y * chord.x;
                // Compare both products and their exact rounding residuals.
                // A rounded-zero determinant alone can accept a narrow curve.
                !left.is_finite()
                    || !right.is_finite()
                    || left.is_subnormal()
                    || right.is_subnormal()
                    || (left == 0.0 && offset.x != 0.0 && chord.y != 0.0)
                    || (right == 0.0 && offset.y != 0.0 && chord.x != 0.0)
                    || left != right
                    || offset.x.mul_add(chord.y, -left) != offset.y.mul_add(chord.x, -right)
            })
        {
            return segment.inv_arclen(distance, self.accuracy);
        }
        let [a, b] = controls.map(coordinate);
        if !(0.0 <= a && a <= b && b <= 1.0)
            || a.is_subnormal()
            || b.is_subnormal()
            || (a == 0.0 && controls[0] != segment.start())
            || (b == 0.0 && controls[1] != segment.start())
        {
            return segment.inv_arclen(distance, self.accuracy);
        }
        let target = distance / length;
        if !target.is_finite() {
            return segment.inv_arclen(distance, self.accuracy);
        }
        let (c1, c2, c3) = match segment {
            PathSeg::Quad(_) => (2.0 * a, 1.0 - 2.0 * a, 0.0),
            PathSeg::Cubic(_) => (3.0 * a, 3.0 * (b - 2.0 * a), 1.0 + 3.0 * (a - b)),
            PathSeg::Line(_) => unreachable!(),
        };
        let roots = kurbo::common::solve_cubic(-target, c1, c2, c3);
        if roots.iter().any(|t| !t.is_finite()) {
            return segment.inv_arclen(distance, self.accuracy);
        }
        let mut roots = roots.into_iter().filter(|t| (0.0..=1.0).contains(t));
        if let Some(t) = roots.next()
            && roots.next().is_none()
            && (((c3 * t + c2) * t + c1) * t - target).abs() <= 32.0 * f64::EPSILON
        {
            return t;
        }
        segment.inv_arclen(distance, self.accuracy)
    }

    /// Only subdivide boundary pieces; untouched endpoints stay exact so
    /// rounding at 0/1 cannot introduce new subpaths or change endpoint frames.
    fn retained_segment(segment: PathSeg, start: f64, end: f64) -> PathSeg {
        if start == 0.0 && end == 1.0 {
            return segment;
        }
        let mut retained = segment.subsegment(start..end);
        let (first, last) = match &mut retained {
            PathSeg::Line(line) => (&mut line.p0, &mut line.p1),
            PathSeg::Quad(quad) => (&mut quad.p0, &mut quad.p2),
            PathSeg::Cubic(cubic) => (&mut cubic.p0, &mut cubic.p3),
        };
        if start == 0.0 {
            *first = segment.start();
        }
        if end == 1.0 {
            *last = segment.end();
        }
        retained
    }

    /// Point and tangent at an arc distance, with the same endpoint arithmetic
    /// as batched frames. Empty carriers have no frame.
    pub fn frame(&self, distance: f64) -> Result<Option<PathFrame>, String> {
        let Some(first) = self.segments.first() else {
            return Ok(None);
        };
        if distance.is_nan() {
            return Err("frame distance must not be NaN".to_owned());
        }
        let total = self.length();
        let at = distance.clamp(0.0, total);
        let (segment, at_end) = if at == 0.0 {
            (*first, false)
        } else if at == total {
            (*self.segments.last().unwrap(), true)
        } else if total <= self.accuracy {
            (*first, false)
        } else {
            // Same retained prefix endpoint as trim(), without allocating all
            // preceding segments for every candidate frame.
            let (_, end_outset) = fit_outsets_to_length(0.0, total - at, total);
            let visible_end = total - end_outset;
            let (lengths, _) = self.measured_lengths();
            let mut cursor = 0.0;
            let mut last = None;
            for (segment, length) in self.segments.iter().copied().zip(lengths.iter().copied()) {
                if length == 0.0 {
                    continue;
                }
                let local_end = (visible_end - cursor).min(length);
                cursor += length;
                if local_end <= 0.0 {
                    break;
                }
                let t = if visible_end >= cursor || length - local_end <= f64::EPSILON {
                    1.0
                } else {
                    self.parameter_at_distance(segment, local_end)
                };
                if t > 0.0 {
                    last = Some(Self::retained_segment(segment, 0.0, t));
                }
                if cursor >= visible_end {
                    break;
                }
            }
            (last.unwrap_or(*first), last.is_some())
        };
        Ok(Some(
            CubicBezierSpec::drawable(segment).endpoint_frame(at_end),
        ))
    }

    /// First contact, in the requested traversal direction, with a circle
    /// centered at the station. Bezier convex-hull bounds isolate roots,
    /// including tangencies, before geometric refinement. None means that the
    /// available carrier cannot fit this full-size chord.
    pub fn chord_contact(
        &self,
        station: f64,
        radius: f64,
        forward: bool,
    ) -> Result<Option<(CurvePoint, f64)>, String> {
        let Some(origin) = self.frame(station)? else {
            return Ok(None);
        };
        let origin: Point = origin.point.into();
        if !radius.is_finite() || radius < 0.0 {
            return Err("chord radius must be finite and nonnegative".into());
        }
        if radius == 0.0 {
            return Ok(Some((origin.into(), station.clamp(0.0, self.length()))));
        }
        let station = station.clamp(0.0, self.length());
        let (lengths, _) = self.measured_lengths();
        let mut cursor = if forward { 0.0 } else { self.length() };
        for ordinal in 0..self.segments.len() {
            let index = if forward {
                ordinal
            } else {
                self.segments.len() - 1 - ordinal
            };
            let length = lengths[index];
            let start = if forward { cursor } else { cursor - length };
            cursor = if forward { cursor + length } else { start };
            let end = start + length;
            if length == 0.0 || (forward && end < station) || (!forward && start > station) {
                continue;
            }
            let segment = self.segments[index];
            let at = if station <= start {
                0.0
            } else if station >= end {
                1.0
            } else {
                self.parameter_at_distance(segment, station - start)
            };
            let (a, b) = if forward { (at, 1.0) } else { (at, 0.0) };
            if let Some(t) =
                Self::circle_contact_parameter(segment, origin, radius, a, b, self.accuracy, 0)?
            {
                let point = segment.eval(t);
                if (point.distance(origin) - radius).abs() > self.accuracy {
                    return Err("chord contact refinement failed its distance check".into());
                }
                let distance = start
                    + if t == 0.0 {
                        0.0
                    } else if t == 1.0 {
                        length
                    } else {
                        segment.subsegment(0.0..t).arclen(self.accuracy)
                    };
                return Ok(Some((point.into(), distance)));
            }
        }
        Ok(None)
    }

    fn circle_contact_parameter(
        segment: PathSeg,
        origin: Point,
        radius: f64,
        a: f64,
        b: f64,
        accuracy: f64,
        depth: usize,
    ) -> Result<Option<f64>, String> {
        if let PathSeg::Line(line) = segment {
            if line.p0.distance(origin).max(line.p1.distance(origin)) < radius {
                return Ok(None);
            }
            let offset = line.p0 - origin;
            let delta = line.p1 - line.p0;
            // Collinear stations are the common straight-carrier case. Avoid
            // quadratic cancellation and keep the original trim arithmetic.
            if offset.cross(delta) == 0.0 {
                let center = -offset.dot(delta) / delta.hypot2();
                let step = radius / delta.hypot();
                let root = [center - step, center + step]
                    .into_iter()
                    .filter(|t| *t >= a.min(b) && *t <= a.max(b))
                    .min_by(|x, y| {
                        if a <= b {
                            x.total_cmp(y)
                        } else {
                            y.total_cmp(x)
                        }
                    });
                return Ok(root);
            }
            let roots = kurbo::common::solve_quadratic(
                offset.hypot2() - radius * radius,
                2.0 * offset.dot(delta),
                delta.hypot2(),
            );
            let root = roots
                .into_iter()
                .filter(|t| *t >= a.min(b) && *t <= a.max(b))
                .min_by(|x, y| {
                    if a <= b {
                        x.total_cmp(y)
                    } else {
                        y.total_cmp(x)
                    }
                });
            return Ok(root);
        }
        let piece = segment.subsegment(a..b).to_cubic();
        let controls = [piece.p0, piece.p1, piece.p2, piece.p3];
        let upper = controls
            .into_iter()
            .map(|point| point.distance(origin))
            .fold(0.0, f64::max);
        // The distance norm is convex; the control hull bounds all points on
        // this interval. Tolerance here admits floating-point tangent contacts.
        if upper < radius - radius * f64::EPSILON * 4.0 {
            return Ok(None);
        }
        let bounds = PathSeg::Cubic(piece).bounding_box();
        let nearest = Point::new(
            origin.x.clamp(bounds.x0, bounds.x1),
            origin.y.clamp(bounds.y0, bounds.y1),
        );
        if nearest.distance(origin) > radius + accuracy {
            return Ok(None);
        }
        if bounds.width().hypot(bounds.height()) <= accuracy {
            let values = [a, (a + b) / 2.0, b];
            let best = values
                .into_iter()
                .min_by(|x, y| {
                    (segment.eval(*x).distance(origin) - radius)
                        .abs()
                        .total_cmp(&(segment.eval(*y).distance(origin) - radius).abs())
                })
                .unwrap();
            if (segment.eval(best).distance(origin) - radius).abs() <= accuracy {
                return Ok(Some(best));
            }
            return Ok(None);
        }
        if depth >= 64 {
            return Err("chord contact refinement did not converge".into());
        }
        let middle = (a + b) / 2.0;
        if let Some(t) =
            Self::circle_contact_parameter(segment, origin, radius, a, middle, accuracy, depth + 1)?
        {
            return Ok(Some(t));
        }
        Self::circle_contact_parameter(segment, origin, radius, middle, b, accuracy, depth + 1)
    }

    pub fn trim(&self, start_outset: f64, end_outset: f64) -> Result<Vec<PathSeg>, String> {
        if start_outset == 0.0 && end_outset == 0.0 {
            return Ok(self.segments.clone());
        }
        let (lengths, length) = self.measured_lengths();
        let length = *length;
        if length <= f64::EPSILON {
            return Ok(self.segments.clone());
        }
        let (start_outset, end_outset) = fit_outsets_to_length(start_outset, end_outset, length);
        let visible_start = start_outset;
        let visible_end = length - end_outset;
        if visible_end <= visible_start {
            return Err("path trim distances leave no visible path".to_string());
        }

        let mut cursor = 0.0;
        let mut trimmed = Vec::new();
        for (segment, segment_length) in self.segments.iter().copied().zip(lengths.iter().copied())
        {
            let segment_start = cursor;
            let segment_end = cursor + segment_length;
            cursor = segment_end;

            // A positive source span can be smaller than the cumulative
            // cursor's resolution. Keep interior connectors verbatim: dropping
            // one would disconnect its neighbours despite a connected input.
            // Actual zero-length segments and trim-boundary connectors retain
            // the ordinary empty-intersection behavior.
            if segment_length > 0.0
                && segment_end == segment_start
                && visible_start < segment_start
                && segment_end < visible_end
            {
                trimmed.push(segment);
                continue;
            }

            let keep_start = visible_start.max(segment_start);
            let keep_end = visible_end.min(segment_end);
            if keep_end <= keep_start {
                continue;
            }

            let local_start = keep_start - segment_start;
            let local_end = keep_end - segment_start;
            let t0 = if keep_start == segment_start || local_start <= f64::EPSILON {
                0.0
            } else {
                self.parameter_at_distance(segment, local_start)
            };
            let t1 = if keep_end == segment_end || segment_length - local_end <= f64::EPSILON {
                1.0
            } else {
                self.parameter_at_distance(segment, local_end)
            };
            if t1 > t0 {
                trimmed.push(Self::retained_segment(segment, t0, t1));
            }
        }

        if trimmed.is_empty() {
            return Err("path trimming produced no visible path".to_string());
        }
        Ok(trimmed)
    }
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

impl HobbyThroughSpec {
    /// The Hobby curve from start through one point to end.
    pub fn curve(self) -> Result<CurvePathOutput, String> {
        hobby_to_cubic_open(
            &[self.start.into(), self.through.into(), self.end.into()],
            self.omega,
        )
        .and_then(|segments| cubic_spline_output(segments, self.accuracy))
    }
}

impl HobbySplineSpec {
    /// The open Hobby spline through every point.
    pub fn curve(self) -> Result<CurvePathOutput, String> {
        let accuracy = self.accuracy;
        let points = self.points.into_iter().map(Into::into).collect::<Vec<_>>();
        hobby_to_cubic_open(&points, self.omega)
            .and_then(|segments| cubic_spline_output(segments, accuracy))
    }
}

fn cubic_spline_output(segments: Vec<CubicBez>, accuracy: f64) -> Result<CurvePathOutput, String> {
    validate_positive_accuracy(accuracy)?;
    curve_path_from_segments(segments.into_iter().map(PathSeg::Cubic))
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

#[derive(Deserialize)]
pub struct RegionPart {
    pub segments: Vec<CubicBezierSpec>,
    pub visible: bool,
}

#[derive(Deserialize)]
pub struct RegionSamplesSpec {
    pub parts: Vec<RegionPart>,
    pub regions: usize,
    #[serde(deserialize_with = "deserialize_f64")]
    pub unit: f64,
    #[serde(deserialize_with = "deserialize_f64")]
    pub step: f64,
    #[serde(
        default = "default_arclen_accuracy",
        deserialize_with = "deserialize_f64"
    )]
    pub accuracy: f64,
}

impl RegionSamplesSpec {
    /// Sample points along each visible part, grouped by arc-length region.
    pub fn samples(self) -> Result<Vec<(usize, Vec<CurvePoint>)>, String> {
        let accuracy = validate_positive_accuracy(self.accuracy)?;
        if self.regions == 0
            || !self.unit.is_finite()
            || !(self.step.is_finite() && self.step > 0.0)
        {
            return Err("Region sampling needs positive regions/step and a finite unit".into());
        }
        let paths: Vec<_> = self
            .parts
            .iter()
            .map(|part| {
                let mut path = BezPath::new();
                let mut previous = None;
                for &segment in &part.segments {
                    if previous != Some(segment.start) {
                        path.move_to(segment.start);
                    }
                    path.curve_to(segment.control_start, segment.control_end, segment.end);
                    previous = Some(segment.end);
                }
                path
            })
            .collect();
        let lengths: Vec<f64> = paths
            .iter()
            .map(|p| p.segments().map(|s| s.arclen(accuracy)).sum())
            .collect();
        let total: f64 = lengths.iter().sum();
        let mut offset = 0.0;
        let mut samples = Vec::new();
        for ((part, path), length) in self.parts.iter().zip(paths).zip(lengths) {
            if part.visible && length > 0.0 {
                for region in 0..self.regions {
                    let start = 0.0f64.max(total * region as f64 / self.regions as f64 - offset);
                    let end =
                        length.min(total * (region + 1) as f64 / self.regions as f64 - offset);
                    if start >= end {
                        continue;
                    }
                    let trimmed = TrimPathSpec {
                        path: path.clone(),
                        start_outset: start,
                        end_outset: length - end,
                        accuracy,
                    }
                    .trimmed()?;
                    let mut points = Vec::new();
                    for segment in trimmed.path.segments().map(CubicBezierSpec::drawable) {
                        let distance = |a: CurvePoint, b: CurvePoint| {
                            let dx = a.x - b.x;
                            let dy = a.y - b.y;
                            (dx * dx + dy * dy).sqrt()
                        };
                        let speed = 3.0
                            * distance(segment.start, segment.control_start)
                                .max(distance(segment.control_start, segment.control_end))
                                .max(distance(segment.control_end, segment.end));
                        let steps = (speed * self.unit.abs() / self.step).ceil().max(1.0);
                        if !(steps.is_finite() && steps < usize::MAX as f64) {
                            return Err("Region sampling count is too large".into());
                        }
                        let steps = steps as usize;
                        for step in 0..=steps {
                            let t = step as f64 / steps as f64;
                            let lerp = |a: CurvePoint, b: CurvePoint| CurvePoint {
                                x: a.x + (b.x - a.x) * t,
                                y: a.y + (b.y - a.y) * t,
                            };
                            let ab = lerp(segment.start, segment.control_start);
                            let bc = lerp(segment.control_start, segment.control_end);
                            let cd = lerp(segment.control_end, segment.end);
                            points.push(lerp(lerp(ab, bc), lerp(bc, cd)));
                        }
                    }
                    samples.push((region, points));
                }
            }
            offset += length;
        }
        Ok(samples)
    }
}

pub fn curve_region_samples_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: RegionSamplesSpec = ciborium::de::from_reader(arg)
        .map_err(|err| format!("Failed to deserialize region sample spec: {err}"))?;
    encode_cbor(&spec.samples()?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cetz_pattern_promotes_quadratics_before_scaling() {
        let mut path = BezPath::new();
        path.move_to((-2.5, 4.0));
        path.quad_to((0.3, -1.2), (7.1, 0.8));
        let output = PatternPathOutput {
            path,
            pattern: "wave".into(),
            parts: Vec::new(),
        };
        let PatternCetzOutput::Data(subpaths) = PatternCetzOutput::from_pattern(output, -1.75)
        else {
            panic!("expected compact output")
        };
        assert_eq!(subpaths.len(), 1);
        let CetzSegment::Cubic((kind, first, second, end)) = &subpaths[0].2[0] else {
            panic!("expected cubic promotion")
        };
        assert_eq!(*kind, "c");
        let expected_first: [f64; 3] = [
            (-2.5 + (0.3 + 2.5) * 2.0 / 3.0) * -1.75,
            (4.0 + (-1.2 - 4.0) * 2.0 / 3.0) * -1.75,
            0.0,
        ];
        let expected_second: [f64; 3] = [
            (7.1 + (0.3 - 7.1) * 2.0 / 3.0) * -1.75,
            (0.8 + (-1.2 - 0.8) * 2.0 / 3.0) * -1.75,
            0.0,
        ];
        assert_eq!(first.map(f64::to_bits), expected_first.map(f64::to_bits));
        assert_eq!(second.map(f64::to_bits), expected_second.map(f64::to_bits));
        assert_eq!(*end, [7.1 * -1.75, 0.8 * -1.75, 0.0]);
    }

    #[test]
    fn cetz_pattern_keeps_subpath_boundaries_and_signed_zeros() {
        let mut path = BezPath::new();
        path.move_to((-0.0, 0.0));
        path.line_to((1.0, 2.0));
        path.close_path();
        path.move_to((3.0, 4.0));
        let output = PatternPathOutput {
            path,
            pattern: "zigzag".into(),
            parts: Vec::new(),
        };
        let PatternCetzOutput::Data(subpaths) = PatternCetzOutput::from_pattern(output, -1.0)
        else {
            panic!("expected compact output")
        };
        assert_eq!(subpaths.len(), 2);
        assert_eq!(
            subpaths[0].0.map(f64::to_bits),
            [0.0f64, -0.0, 0.0].map(f64::to_bits)
        );
        assert!(subpaths[0].1);
        assert_eq!(subpaths[0].2.len(), 1);
        assert!(!subpaths[1].1);
        assert!(subpaths[1].2.is_empty());
    }

    #[test]
    fn cetz_pattern_retains_original_commands_for_implicit_integer_origins() {
        let path = BezPath::from_vec(vec![
            PathEl::MoveTo(Point::ZERO),
            PathEl::ClosePath,
            PathEl::LineTo(Point::ZERO),
        ]);
        let output = PatternPathOutput {
            path,
            pattern: "coil".into(),
            parts: Vec::new(),
        };
        let expected = encode_cbor(&output).unwrap();
        for unit in [-1.0, -0.0, 0.0, 1.0] {
            let converted = PatternCetzOutput::from_pattern(output.clone(), unit);
            assert!(matches!(converted, PatternCetzOutput::Original(_)));
            assert_eq!(encode_cbor(&converted).unwrap(), expected);
        }
    }

    #[test]
    fn batched_frames_clamp_and_handle_empty_or_short_paths() {
        let samples = |path, distances: &[f64], accuracy| {
            PathFramesSpec {
                path,
                distances: distances.to_vec(),
                accuracy,
            }
            .frames()
            .unwrap()
        };
        let line = BezPath::from_path_segments(
            [PathSeg::Line(Line::new((0.0, 0.0), (3.0, 0.0)))].into_iter(),
        );
        let frames = samples(line, &[-1.0, 0.0, 1.5, 3.0, 4.0], 1e-6);
        assert_eq!(frames[0], frames[1]);
        assert_eq!(frames[3], frames[4]);
        for (frame, x) in frames.iter().zip([0.0, 0.0, 1.5, 3.0, 3.0]) {
            assert_point_close(frame.as_ref().unwrap().point, point(x, 0.0));
        }
        assert_eq!(samples(BezPath::new(), &[0.0, 1.0], 1e-6), [None, None]);
        let short = BezPath::from_path_segments(
            [PathSeg::Line(Line::new((0.0, 0.0), (1e-8, 0.0)))].into_iter(),
        );
        let frames = samples(short, &[5e-9, 1e-8], 1e-6);
        assert_eq!(frames[0].as_ref().unwrap().point, point(0.0, 0.0));
        assert_eq!(frames[1].as_ref().unwrap().point, point(1e-8, 0.0));
    }

    fn point(x: f64, y: f64) -> CurvePoint {
        CurvePoint { x, y }
    }

    #[test]
    fn trims_preserve_positive_source_connectors_below_cumulative_resolution() {
        let a = Point::new(1.0, 0.0);
        let b = Point::new(f64::from_bits(1.0f64.to_bits() - 1), 0.0);
        let first = PathSeg::Line(Line::new((-1.0, 0.0), a));
        let connector = PathSeg::Cubic(CubicBez::new(a, a, b, b));
        let zero = PathSeg::Cubic(CubicBez::new(b, b, b, b));
        let last = PathSeg::Line(Line::new(b, (3.0, 0.0)));
        let source = [first, connector, zero, last];
        let trimmer = PathTrimmer::new(source, 1e-6);
        let (lengths, total) = trimmer.measured_lengths();
        assert!(lengths[1] > 0.0);
        assert_eq!(lengths[0] + lengths[1], lengths[0]);
        assert_eq!(lengths[2], 0.0);
        assert_eq!(*total, 4.0);
        assert_eq!(trimmer.trim(0.0, 0.0).unwrap(), source);

        let interior = trimmer.trim(0.5, 0.5).unwrap();
        assert_eq!(interior.len(), 3);
        assert_eq!(interior[1], connector);
        assert_eq!(interior[0].end(), interior[1].start());
        assert_eq!(interior[1].end(), interior[2].start());
        // Neither positive collapsed spans nor actual zeros belong to a
        // window whose boundary lands at their cumulative station.
        assert_eq!(trimmer.trim(0.0, 2.0).unwrap(), [first]);
        assert_eq!(trimmer.trim(2.0, 0.0).unwrap(), [last]);
        for distance in [1.5, 2.0, 2.5, 3.5] {
            let prefix = trimmer.trim(0.0, trimmer.length() - distance).unwrap();
            assert_eq!(
                trimmer.frame(distance).unwrap().unwrap(),
                CubicBezierSpec::drawable(*prefix.last().unwrap()).endpoint_frame(true),
            );
        }

        // Without a source connector the exact endpoint gap is intentional,
        // even if it is only one ULP. Trimming must not invent a join.
        let disjoint = PathTrimmer::new([first, last], 1e-6);
        let trimmed = disjoint.trim(0.5, 0.5).unwrap();
        assert_eq!(trimmed.len(), 2);
        assert_eq!(trimmed[0].end(), a);
        assert_eq!(trimmed[1].start(), b);
        let path = BezPath::from_path_segments(trimmed.into_iter());
        assert_eq!(
            path.elements()
                .iter()
                .filter(|el| matches!(el, PathEl::MoveTo(_)))
                .count(),
            2,
        );
    }

    #[test]
    fn straight_bezier_frames_and_trims_use_exact_distance_without_rewriting_controls() {
        let segments = [
            PathSeg::Cubic(CubicBez::new(
                (0.0, 0.0),
                (0.0, 0.0),
                (3.0, 0.0),
                (3.0, 0.0),
            )),
            PathSeg::Cubic(CubicBez::new(
                (0.0, 0.0),
                (0.1, 0.0),
                (1.8, 0.0),
                (3.0, 0.0),
            )),
            PathSeg::Quad(kurbo::QuadBez::new((0.0, 0.0), (0.2, 0.0), (3.0, 0.0))),
        ];
        for segment in segments {
            let trimmer = PathTrimmer::new([segment], 1e-6);
            assert_eq!(trimmer.segments, [segment]);
            for distance in [0.21, 1.43, 1.57, 2.73] {
                let frame = trimmer.frame(distance).unwrap().unwrap();
                assert!((frame.point.x - distance).abs() < 1e-12);
                assert_eq!(frame.point.y, 0.0);
                let prefix = trimmer.trim(0.0, trimmer.length() - distance).unwrap();
                let retained = *prefix.last().unwrap();
                assert_eq!(
                    frame,
                    CubicBezierSpec::drawable(retained).endpoint_frame(true)
                );
                assert!(matches!(
                    (segment, retained),
                    (PathSeg::Cubic(_), PathSeg::Cubic(_)) | (PathSeg::Quad(_), PathSeg::Quad(_))
                ));
                let t = trimmer.parameter_at_distance(segment, distance);
                let expected =
                    CubicBezierSpec::drawable(segment.subsegment(0.0..t)).endpoint_frame(true);
                assert_eq!(frame.tangent, expected.tangent);
                assert!(frame.tangent.x > 0.0 && frame.tangent.y == 0.0);
            }
            for forward in [false, true] {
                let (contact, distance) =
                    trimmer.chord_contact(1.43, 0.2, forward).unwrap().unwrap();
                assert!((contact.x - distance).abs() < 1e-12);
                let frame = trimmer.frame(distance).unwrap().unwrap();
                assert!((frame.point.x - contact.x).abs() < 1e-12);
                assert!(((contact.x - 1.43).abs() - 0.2).abs() <= 1e-6);
            }
        }
    }

    #[test]
    fn straight_parameter_fast_path_rejects_backtracking_collapsed_and_narrow_curves() {
        let segments = [
            CubicBez::new((0.0, 0.0), (4.0, 0.0), (-1.0, 0.0), (3.0, 0.0)),
            CubicBez::new((0.0, 0.0), (1.0, 0.0), (-1.0, 0.0), (0.0, 0.0)),
            CubicBez::new((0.0, 0.0), (0.1, 1e-12), (1.8, 0.0), (3.0, 0.0)),
            // Its ordinary floating-point determinant rounds to zero, but
            // the exact products differ: it must not be treated as a line.
            CubicBez::new(
                (0.0, 0.0),
                (0.1, 0.018181818181818184),
                (1.1, 0.2),
                (1.1, 0.2),
            ),
            // Unsafe determinant products (underflow and overflow).
            CubicBez::new(
                (0.0, 0.0),
                (1e-200, 1e-200),
                (3e-200, 3e-200),
                (3e-200, 3e-200),
            ),
            CubicBez::new((0.0, 0.0), (1e200, 1e200), (3e200, 3e200), (3e200, 3e200)),
        ];
        for cubic in segments {
            let segment = PathSeg::Cubic(cubic);
            let trimmer = PathTrimmer::new([segment], 1e-6);
            let distance = trimmer.length() * 0.37;
            assert_eq!(
                trimmer.parameter_at_distance(segment, distance),
                segment.inv_arclen(distance, 1e-6),
                "{cubic:?}"
            );
        }
    }

    #[test]
    fn reusable_frames_preserve_prefix_trim_arithmetic_without_allocating_prefixes() {
        let mut path = BezPath::new();
        path.move_to((0.0, 0.0));
        path.line_to((0.0, 0.0));
        path.curve_to((0.0, 8.0), (12.0, 8.0), (12.0, 4.0));
        path.line_to((12.0, 4.0));
        path.quad_to((16.0, -8.0), (20.0, 0.0));
        path.line_to((20.0, 0.0));
        for scale in [1.0, 1e-12] {
            let mut scaled = path.clone();
            scaled.apply_affine(kurbo::Affine::scale(scale));
            let segments: Vec<_> = scaled.segments().collect();
            let trimmer = PathTrimmer::new(segments.iter().copied(), 1e-6);
            let total = trimmer.length();
            for fraction in [0.0, 1e-11, 0.23, 0.5, 0.9, 1.0] {
                let at = total * fraction;
                let (segment, at_end) = if at == 0.0 {
                    (segments[0], false)
                } else if at == total {
                    (*segments.last().unwrap(), true)
                } else if total <= 1e-6 {
                    (segments[0], false)
                } else {
                    let prefix = trimmer.trim(0.0, total - at).unwrap();
                    (*prefix.last().unwrap_or(&segments[0]), !prefix.is_empty())
                };
                let expected = CubicBezierSpec::drawable(segment).endpoint_frame(at_end);
                assert_eq!(trimmer.frame(at).unwrap().unwrap(), expected);
            }
        }
    }

    #[test]
    fn chord_refinement_cap_returns_an_error_not_an_unchecked_contact() {
        let segment = PathSeg::Cubic(CubicBez::new(
            (0.0, 0.0),
            (0.0, 8.0),
            (12.0, 8.0),
            (12.0, 4.0),
        ));
        let result =
            PathTrimmer::circle_contact_parameter(segment, Point::ZERO, 5.0, 0.0, 1.0, 1e-6, 64);
        assert!(result.unwrap_err().contains("did not converge"));
    }

    fn assert_point_close(actual: CurvePoint, expected: CurvePoint) {
        assert!((actual.x - expected.x).abs() < 1e-6);
        assert!((actual.y - expected.y).abs() < 1e-6);
    }

    fn path_points(path: &BezPath) -> Vec<CurvePoint> {
        let mut points = Vec::new();
        for segment in path.segments() {
            if points.is_empty() {
                points.push(segment.start().into());
            }
            points.push(segment.end().into());
        }
        points
    }

    fn path_cubics(path: &BezPath) -> Vec<CubicBezierSpec> {
        path.segments()
            .map(|segment| segment.to_cubic().into())
            .collect()
    }

    fn path_length_value(path: &BezPath, accuracy: f64) -> f64 {
        path.segments()
            .map(|segment| segment.arclen(accuracy))
            .sum()
    }

    fn cbor_map_keys(value: &ciborium::Value) -> Vec<&str> {
        match value {
            ciborium::Value::Map(entries) => entries
                .iter()
                .filter_map(|(key, _)| match key {
                    ciborium::Value::Text(key) => Some(key.as_str()),
                    _ => None,
                })
                .collect(),
            _ => panic!("expected CBOR map"),
        }
    }

    fn cbor_map_get<'a>(value: &'a ciborium::Value, field: &str) -> &'a ciborium::Value {
        match value {
            ciborium::Value::Map(entries) => entries
                .iter()
                .find_map(|(key, value)| match key {
                    ciborium::Value::Text(key) if key == field => Some(value),
                    _ => None,
                })
                .unwrap_or_else(|| panic!("missing CBOR field `{field}`")),
            _ => panic!("expected CBOR map"),
        }
    }

    fn assert_cbor_point_tuple(value: &ciborium::Value) {
        match value {
            ciborium::Value::Array(values) => assert_eq!(values.len(), 2),
            _ => panic!("expected CBOR point tuple"),
        }
    }

    #[test]
    fn wire_path_uses_curve_command_elements() {
        let path = straight_path(3.0);
        let wire = WirePath::from(&path);

        assert_eq!(wire.elements.len(), 2);
        match &wire.elements[0] {
            WirePathElement::Move { start } => assert_point_close(*start, point(0.0, 0.0)),
            _ => panic!("expected move element"),
        }
        match &wire.elements[1] {
            WirePathElement::Cubic {
                control_start,
                control_end,
                end,
            } => {
                assert_point_close(*control_start, point(1.0, 0.0));
                assert_point_close(*control_end, point(2.0, 0.0));
                assert_point_close(*end, point(3.0, 0.0));
            }
            _ => panic!("expected cubic element"),
        }

        let roundtrip = BezPath::from(wire);
        assert_eq!(roundtrip.elements(), path.elements());
    }

    #[test]
    fn wire_path_preserves_quadratic_curve_commands() {
        let mut path = BezPath::new();
        path.move_to(Point::new(0.0, 0.0));
        path.quad_to(Point::new(1.0, 3.0), Point::new(3.0, 0.0));

        let wire = WirePath::from(&path);
        assert_eq!(wire.elements.len(), 2);
        match &wire.elements[1] {
            WirePathElement::Quad { control, end } => {
                assert_point_close(*control, point(1.0, 3.0));
                assert_point_close(*end, point(3.0, 0.0));
            }
            _ => panic!("expected quadratic element"),
        }

        let roundtrip = BezPath::from(wire);
        assert!(matches!(roundtrip.elements()[1], PathEl::QuadTo(_, _)));
    }

    #[test]
    fn wire_close_is_straight_and_rejects_modes() {
        let bytes = encode_cbor(&std::collections::BTreeMap::from([("kind", "close")])).unwrap();
        let close: WirePathElement = ciborium::de::from_reader(&bytes[..]).unwrap();
        assert_eq!(close, WirePathElement::Close {});

        let mut path = BezPath::new();
        path.move_to((0.0, 0.0));
        path.line_to((3.0, 0.0));
        path.line_to((3.0, 4.0));
        path.close_path();
        let wire = WirePath::from(&path);
        assert_eq!(wire.elements.last(), Some(&WirePathElement::Close {}));
        assert_eq!(BezPath::from(wire).elements(), path.elements());
        assert!((path_length_value(&path, 1e-6) - 12.0).abs() < 1e-6);

        for mode in ["straight", "smooth"] {
            let bytes = encode_cbor(&std::collections::BTreeMap::from([
                ("kind", "close"),
                ("mode", mode),
            ]))
            .unwrap();
            let result: Result<WirePathElement, _> = ciborium::de::from_reader(&bytes[..]);
            assert!(result.is_err(), "close must reject mode {mode}");
        }
    }

    #[test]
    fn curve_point_rejects_dictionary_shape() {
        let input = std::collections::BTreeMap::from([("x", 0.0), ("y", 1.0)]);
        let bytes = encode_cbor(&input).unwrap();
        let result: Result<CurvePoint, _> = ciborium::de::from_reader(&bytes[..]);

        assert!(result.is_err());
    }

    fn segment_path(segment: impl Into<PathSeg>) -> BezPath {
        BezPath::from_path_segments(std::iter::once(segment.into()))
    }

    fn intersection_spec(a: BezPath, b: BezPath) -> PathIntersectionsSpec {
        PathIntersectionsSpec {
            a,
            b,
            accuracy: 1e-6,
        }
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!((actual - expected).abs() < 2e-5, "{actual} != {expected}");
    }

    #[test]
    fn line_cubic_intersection_reports_point_and_arc_distances() {
        let a = segment_path(Line::new((-2.0, 0.0), (2.0, 0.0)));
        let b = segment_path(CubicBez::new(
            (0.0, -2.0),
            (0.0, -1.0),
            (0.0, 1.0),
            (0.0, 2.0),
        ));

        let hits = path_intersections(intersection_spec(a, b)).unwrap();
        assert_eq!(hits.len(), 1);
        assert_point_close(hits[0].point, point(0.0, 0.0));
        assert_close(hits[0].distance_a, 2.0);
        assert_close(hits[0].distance_b, 2.0);
    }

    #[test]
    fn cubic_cubic_intersection_finds_all_crossings() {
        let a = segment_path(CubicBez::new(
            (-2.0, 0.0),
            (-1.0, 0.0),
            (1.0, 0.0),
            (2.0, 0.0),
        ));
        let b = segment_path(CubicBez::new(
            (-2.0, -1.0),
            (-2.0, 4.0),
            (2.0, -4.0),
            (2.0, 1.0),
        ));

        let hits = path_intersections(intersection_spec(a, b)).unwrap();
        assert_eq!(hits.len(), 3, "{hits:#?}");
        assert!(
            hits.windows(2)
                .all(|pair| pair[0].distance_a < pair[1].distance_a)
        );
        assert_point_close(hits[1].point, point(0.0, 0.0));
    }

    #[test]
    fn collinear_backtracking_cubic_keeps_distinct_parameter_crossings() {
        let a = segment_path(Line::new((0.0, -1.0), (0.0, 1.0)));
        let b = segment_path(CubicBez::new(
            (-1.0, 0.0),
            (4.0, 0.0),
            (-4.0, 0.0),
            (1.0, 0.0),
        ));

        let hits = path_intersections(intersection_spec(a, b)).unwrap();
        assert_eq!(hits.len(), 3, "{hits:#?}");
        assert!(
            hits.windows(2)
                .all(|pair| pair[0].distance_b < pair[1].distance_b)
        );
    }

    #[test]
    fn intersection_coalesces_segment_boundary_duplicates() {
        let mut a = BezPath::new();
        a.move_to((-3.0, 0.0));
        a.line_to((0.0, 0.0));
        a.line_to((3.0, 0.0));
        let mut b = BezPath::new();
        for x in [-2.0, 0.0, 2.0] {
            b.move_to((x, -1.0));
            b.line_to((x, 1.0));
        }

        let hits = path_intersections(intersection_spec(a, b)).unwrap();
        assert_eq!(hits.len(), 3, "{hits:#?}");
        assert_close(hits[0].point.x, -2.0);
        assert_close(hits[1].point.x, 0.0);
        assert_close(hits[2].point.x, 2.0);
    }

    #[test]
    fn intersection_filters_endpoint_contacts_tangencies_and_misses() {
        let horizontal = segment_path(Line::new((-2.0, 0.0), (2.0, 0.0)));
        let shared = segment_path(Line::new((-2.0, 0.0), (-2.0, 1.0)));
        let t_junction = segment_path(Line::new((0.0, 0.0), (0.0, 1.0)));
        let tangent = segment_path(kurbo::QuadBez::new((-1.0, 1.0), (0.0, -1.0), (1.0, 1.0)));
        let miss = segment_path(Line::new((-2.0, 1.0), (2.0, 1.0)));

        assert!(
            path_intersections(intersection_spec(horizontal.clone(), shared))
                .unwrap()
                .is_empty()
        );
        assert!(
            path_intersections(intersection_spec(horizontal.clone(), t_junction))
                .unwrap()
                .is_empty()
        );
        assert!(
            path_intersections(intersection_spec(horizontal.clone(), tangent))
                .unwrap()
                .is_empty()
        );
        assert!(
            path_intersections(intersection_spec(horizontal, miss))
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn reversing_a_path_reverses_its_intersection_distance() {
        let a = segment_path(Line::new((-2.0, 0.0), (2.0, 0.0)));
        let forward = segment_path(Line::new((-1.0, -3.0), (-1.0, 1.0)));
        let reversed = segment_path(Line::new((-1.0, 1.0), (-1.0, -3.0)));

        let forward_hit = path_intersections(intersection_spec(a.clone(), forward)).unwrap()[0];
        let reverse_hit = path_intersections(intersection_spec(a, reversed)).unwrap()[0];
        assert_close(forward_hit.distance_b, 3.0);
        assert_close(reverse_hit.distance_b, 1.0);
    }

    #[test]
    fn intersection_rejects_invalid_accuracy_and_geometry() {
        let finite = segment_path(Line::new((0.0, 0.0), (1.0, 0.0)));
        for accuracy in [0.0, f64::NAN] {
            let mut spec = intersection_spec(finite.clone(), finite.clone());
            spec.accuracy = accuracy;
            assert!(path_intersections(spec).is_err());
        }
        let invalid = segment_path(Line::new((f64::INFINITY, 0.0), (1.0, 0.0)));
        assert!(path_intersections(intersection_spec(invalid, finite.clone())).is_err());
        let overflowing = segment_path(Line::new((-f64::MAX, 0.0), (f64::MAX, 0.0)));
        assert!(path_intersections(intersection_spec(overflowing, finite)).is_err());

        let mut overflowing_total = BezPath::new();
        for y in [0.0, 1.0] {
            overflowing_total.move_to((0.0, y));
            overflowing_total.line_to((f64::MAX * 0.75, y));
        }
        let finite = segment_path(Line::new((0.0, 0.0), (1.0, 0.0)));
        assert!(path_intersections(intersection_spec(overflowing_total, finite)).is_err());
    }

    #[test]
    fn intersection_cbor_api_uses_typed_records() {
        let spec = intersection_spec(
            segment_path(Line::new((-1.0, 0.0), (1.0, 0.0))),
            segment_path(Line::new((0.0, -1.0), (0.0, 1.0))),
        );
        let bytes = encode_cbor(&spec).unwrap();
        let output_bytes = curve_path_intersections_bytes(&bytes).unwrap();
        let output: Vec<PathIntersection> = ciborium::de::from_reader(&output_bytes[..]).unwrap();
        let value: ciborium::Value = ciborium::de::from_reader(&output_bytes[..]).unwrap();

        assert_eq!(output.len(), 1);
        assert_point_close(output[0].point, point(0.0, 0.0));
        let record = match &value {
            ciborium::Value::Array(records) => &records[0],
            _ => panic!("expected CBOR record array"),
        };
        assert_eq!(
            cbor_map_keys(record),
            vec![
                "point",
                "distance-a",
                "distance-b",
                "segment-a",
                "segment-b",
                "t-a",
                "t-b",
            ]
        );
        assert_cbor_point_tuple(cbor_map_get(record, "point"));
    }

    #[test]
    fn trim_path_trims_by_curve_arclength() {
        let output = TrimPathSpec {
            path: straight_path(3.0),
            start_outset: 0.25,
            end_outset: 0.5,
            accuracy: 1e-6,
        }
        .trimmed()
        .unwrap();
        let points = path_points(&output.path);

        assert_point_close(points[0], point(0.25, 0.0));
        assert_point_close(*points.last().unwrap(), point(2.5, 0.0));
    }

    #[test]
    fn shifted_trim_batch_reuses_measurements_without_changing_windows() {
        let mut disconnected = straight_path(3.0);
        disconnected.move_to((5.0, 1.0));
        disconnected.quad_to((6.0, 3.0), (8.0, 1.0));
        disconnected.close_path();
        let mut curved = BezPath::new();
        curved.move_to((0.0, 0.0));
        curved.curve_to((0.5, 2.0), (2.5, -1.0), (3.0, 0.0));
        for path in [straight_path(3.0), curved, disconnected, straight_path(0.0)] {
            let outsets = [
                (0.0, -0.0),
                (0.25, 0.5),
                (1.25, 0.25),
                (20.0, 30.0),
                (0.0, 0.0),
            ];
            let batch = TrimPathsSpec {
                path: path.clone(),
                outsets: outsets
                    .iter()
                    .map(|&(start_outset, end_outset)| PathOutsets {
                        start_outset,
                        end_outset,
                    })
                    .collect(),
                accuracy: 1e-6,
                format: LayerFormat::Array,
                unit: 1.0,
            }
            .trimmed()
            .unwrap();
            for (actual, (start_outset, end_outset)) in batch.iter().zip(outsets) {
                let expected = TrimPathSpec {
                    path: path.clone(),
                    start_outset,
                    end_outset,
                    accuracy: 1e-6,
                }
                .trimmed()
                .unwrap();
                assert_eq!(
                    encode_cbor(&actual.path).unwrap(),
                    encode_cbor(&expected).unwrap()
                );
                assert_eq!(
                    actual.segments,
                    expected
                        .path
                        .segments()
                        .map(CubicBezierSpec::drawable)
                        .collect::<Vec<_>>()
                );
            }
        }

        let path = straight_path(3.0);
        let trimmer = PathTrimmer::new(path.segments(), 1e-6);
        trimmer.trim(0.0, 0.0).unwrap();
        assert!(trimmer.lengths.get().is_none());
        let first = trimmer.trim(0.25, 0.5).unwrap();
        assert!(trimmer.lengths.get().is_some());
        assert_point_close(first[0].start().into(), point(0.25, 0.0));
        assert_point_close(first.last().unwrap().end().into(), point(2.5, 0.0));
        let second = trimmer.trim(1.0, 0.25).unwrap();
        assert_point_close(second[0].start().into(), point(1.0, 0.0));
        assert_point_close(second.last().unwrap().end().into(), point(2.75, 0.0));
    }

    #[test]
    fn shifted_trim_batch_retains_validation_and_empty_batches() {
        let spec = |outsets, accuracy| TrimPathsSpec {
            path: straight_path(3.0),
            outsets,
            accuracy,
            format: LayerFormat::Array,
            unit: 1.0,
        };
        assert!(spec(vec![], 1e-6).trimmed().unwrap().is_empty());
        assert!(spec(vec![], 0.0).trimmed().is_err());
        for start_outset in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                spec(
                    vec![PathOutsets {
                        start_outset,
                        end_outset: 0.0,
                    }],
                    1e-6,
                )
                .trimmed()
                .is_err()
            );
        }
    }

    #[test]
    fn packed_layers_keep_native_records_and_share_cetz_conversion() {
        let mut path = straight_path(3.0);
        path.line_to((4.0, 1.0));
        let layers = TrimPathsSpec {
            path,
            outsets: vec![
                PathOutsets {
                    start_outset: 0.0,
                    end_outset: 0.0,
                },
                PathOutsets {
                    start_outset: 0.0,
                    end_outset: 2.0,
                },
            ],
            accuracy: 1e-6,
            format: LayerFormat::Cbor,
            unit: -2.0,
        }
        .trimmed()
        .unwrap();
        let packet = PackedPathLayers::new(&layers, -2.0).unwrap();
        assert_eq!(packet.layers, encode_cbor(&layers).unwrap());
        assert_eq!(packet.count, 2);
        assert!(packet.supported && packet.nonzero && !packet.all_single);
        let value: ciborium::Value = ciborium::de::from_reader(&packet.footprints[..]).unwrap();
        let ciborium::Value::Array(paths) = value else {
            panic!("expected footprint paths")
        };
        let expected =
            CetzSubpath::from_path(&layers[0].path.path, |p| [p.x * -2.0, p.y * -2.0]).unwrap();
        assert_eq!(
            encode_cbor(cbor_map_get(&paths[0], "segments")).unwrap(),
            encode_cbor(&expected).unwrap()
        );
        assert_eq!(
            cbor_map_get(&paths[0], "single"),
            &ciborium::Value::Bool(false)
        );
        assert_eq!(
            cbor_map_get(&paths[1], "single"),
            &ciborium::Value::Bool(true)
        );
        let encoded = encode_cbor(&packet).unwrap();
        let value: ciborium::Value = ciborium::de::from_reader(&encoded[..]).unwrap();
        assert_eq!(
            cbor_map_get(&value, "layers"),
            &ciborium::Value::Bytes(packet.layers)
        );
        assert_eq!(
            cbor_map_get(&value, "footprints"),
            &ciborium::Value::Bytes(packet.footprints)
        );
    }

    #[test]
    fn parallel_path_offsets_straight_curve_to_left_normal() {
        let output = ParallelPathSpec {
            path: straight_path(3.0),
            distance: 0.5,
            start_outset: 0.0,
            end_outset: 0.0,
            accuracy: 1e-6,
        }
        .parallel()
        .unwrap();
        let points = path_points(&output.path);

        assert!(points.len() >= 2);
        assert_point_close(points[0], point(0.0, 0.5));
        assert_point_close(*points.last().unwrap(), point(3.0, 0.5));
        assert!(output.path.segments().next().is_some());
    }

    #[test]
    fn parallel_path_accepts_negative_distance() {
        let output = ParallelPathSpec {
            path: straight_path(3.0),
            distance: -0.25,
            start_outset: 0.0,
            end_outset: 0.0,
            accuracy: 1e-6,
        }
        .parallel()
        .unwrap();
        let points = path_points(&output.path);

        assert!((points.first().unwrap().y + 0.25).abs() < 1e-6);
        assert!((points.last().unwrap().y + 0.25).abs() < 1e-6);
    }

    #[test]
    fn parallel_path_trims_offset_path_by_arclength() {
        let output = ParallelPathSpec {
            path: straight_path(4.0),
            distance: 0.25,
            start_outset: 1.0,
            end_outset: 1.0,
            accuracy: 1e-6,
        }
        .parallel()
        .unwrap();
        let points = path_points(&output.path);

        assert_point_close(points[0], point(1.0, 0.25));
        assert_point_close(*points.last().unwrap(), point(3.0, 0.25));
        assert!((path_length_value(&output.path, 1e-6) - 2.0).abs() < 1e-6);
    }

    #[test]
    fn parallel_path_bevels_near_collinear_fitting_gaps_without_snapping() {
        let mut path = BezPath::new();
        path.move_to((0.0, 0.0));
        path.curve_to((0.3, 0.0), (0.7, 0.0), (1.0, 0.0));
        path.curve_to((1.3, 0.0), (1.7, 0.0001), (2.0, 0.0001));
        let fitted: Vec<_> = path
            .segments()
            .flat_map(|segment| offset_path_segment(segment, -0.35, 0.001))
            .collect();
        let gap = fitted[0].end().distance(fitted[1].start());
        assert!(gap > 1e-9 && gap < 0.001, "{gap}");
        assert_eq!(
            BezPath::from_path_segments(fitted.iter().copied())
                .subpaths()
                .count(),
            2
        );
        for trim in [0.0, 0.1] {
            let output = ParallelPathSpec {
                path: path.clone(),
                distance: -0.35,
                start_outset: trim,
                end_outset: trim,
                accuracy: 0.001,
            }
            .parallel()
            .unwrap();
            assert_eq!(output.path.subpaths().count(), 1);
            let segments: Vec<_> = output.path.segments().collect();
            assert_eq!(segments.len(), 3);
            assert_eq!(
                segments[1],
                PathSeg::Line(Line::new(fitted[0].end(), fitted[1].start()))
            );
            if trim == 0.0 {
                assert_eq!(segments[0], fitted[0]);
                assert_eq!(segments[2], fitted[1]);
            }
        }
    }

    #[test]
    fn parallel_path_keeps_fitted_cubic_controls_at_bevel_join() {
        let mut path = BezPath::new();
        path.move_to((0.0, -1.0));
        path.curve_to((0.3, -1.0), (0.7, 0.0), (1.0, 0.0));
        path.curve_to((1.3, 0.0), (1.7, 0.0001), (2.0, 0.0001));
        let fitted: Vec<_> = path
            .segments()
            .map(|segment| offset_path_segment(segment, -0.35, 0.001))
            .collect();
        assert!(
            fitted[0]
                .iter()
                .any(|segment| matches!(segment, PathSeg::Cubic(_)))
        );
        let gap = fitted[0]
            .last()
            .unwrap()
            .end()
            .distance(fitted[1][0].start());
        assert!(gap > 1e-9 && gap < 0.001, "{gap}");
        for trim in [0.0, 0.1] {
            let output = ParallelPathSpec {
                path: path.clone(),
                distance: -0.35,
                start_outset: trim,
                end_outset: trim,
                accuracy: 0.001,
            }
            .parallel()
            .unwrap();
            assert_eq!(output.path.subpaths().count(), 1);
            if trim == 0.0 {
                let segments: Vec<_> = output.path.segments().collect();
                assert_eq!(&segments[..fitted[0].len()], fitted[0]);
                assert_eq!(&segments[fitted[0].len() + 1..], fitted[1]);
            }
        }
    }

    #[test]
    fn parallel_path_preserves_explicit_moves_even_at_touching_endpoints() {
        for next_start in [1.0, 5.0] {
            for trim in [0.0, 0.25] {
                let mut path = BezPath::new();
                path.move_to((0.0, 0.0));
                path.line_to((1.0, 0.0));
                path.move_to((next_start, 0.0));
                path.line_to((next_start + 1.0, 0.0));
                let output = ParallelPathSpec {
                    path,
                    distance: 0.5,
                    start_outset: trim,
                    end_outset: trim,
                    accuracy: 1e-6,
                }
                .parallel()
                .unwrap();
                assert_eq!(output.path.subpaths().count(), 2);
                assert_eq!(output.path.segments().count(), 2);
                assert!((path_length_value(&output.path, 1e-6) - (2.0 - 2.0 * trim)).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn parallel_path_bevels_sharp_corners_and_trims_joined_length() {
        for distance in [-0.5, 0.5] {
            let mut path = BezPath::new();
            path.move_to((0.0, 0.0));
            path.line_to((2.0, 0.0));
            path.line_to((2.0, 2.0));
            for trim in [0.0, 0.25, 2.1] {
                let output = ParallelPathSpec {
                    path: path.clone(),
                    distance,
                    start_outset: trim,
                    end_outset: trim,
                    accuracy: 1e-6,
                }
                .parallel()
                .unwrap();
                assert_eq!(output.path.subpaths().count(), 1);
                let length = 4.0 + distance.abs() * 2.0_f64.sqrt();
                assert!(
                    (path_length_value(&output.path, 1e-6) - (length - 2.0 * trim)).abs() < 1e-6
                );
                if trim == 0.0 {
                    let segments: Vec<_> = output.path.segments().collect();
                    assert_eq!(
                        segments[1],
                        PathSeg::Line(Line::new((2.0, distance), (2.0 - distance, 0.0)))
                    );
                }
            }
        }
    }

    #[test]
    fn parallel_path_retains_meaningful_closure_only_without_trimming() {
        let mut path = BezPath::new();
        path.move_to((0.0, 0.0));
        path.line_to((2.0, 0.0));
        path.line_to((2.0, 2.0));
        path.close_path();
        for trim in [0.0, 0.1] {
            let output = ParallelPathSpec {
                path: path.clone(),
                distance: -0.25,
                start_outset: trim,
                end_outset: 0.0,
                accuracy: 1e-6,
            }
            .parallel()
            .unwrap();
            assert_eq!(
                output.path.subpaths().count(),
                1,
                "trim={trim}: {:?}",
                output.path
            );
            assert_eq!(
                output.path.elements().last() == Some(&PathEl::ClosePath),
                trim == 0.0
            );
            if trim == 0.0 {
                let segments: Vec<_> = output.path.segments().collect();
                assert_eq!(
                    segments.first().unwrap().start(),
                    segments.last().unwrap().end()
                );
            }
        }
    }

    #[test]
    fn parallel_path_floating_joins_keep_moves_closure_and_exact_prefix_frames() {
        let mut path = BezPath::new();
        let origin = Point::new(-0.21169587343232443, 0.39154203702426227);
        path.move_to(origin);
        path.line_to((1.823223304703363, 2.176776695296637));
        path.line_to((-0.17677669529663687, 0.17677669529663687));
        path.close_path();
        // Even a move to the closed subpath's exact start is a boundary.
        path.move_to(origin);
        path.line_to((1.823223304703363, 2.176776695296637));
        for trim in [0.0, 0.1] {
            let output = ParallelPathSpec {
                path: path.clone(),
                distance: -0.35,
                start_outset: trim,
                end_outset: trim,
                accuracy: 0.001,
            }
            .parallel()
            .unwrap();
            assert_eq!(output.path.subpaths().count(), 2);
            assert_eq!(
                output
                    .path
                    .elements()
                    .iter()
                    .filter(|el| matches!(el, PathEl::MoveTo(_)))
                    .count(),
                2
            );
            assert_eq!(
                output
                    .path
                    .elements()
                    .iter()
                    .filter(|el| matches!(el, PathEl::ClosePath))
                    .count(),
                usize::from(trim == 0.0)
            );
            for elements in output.path.subpaths() {
                let subpath = BezPath::from_vec(elements.to_vec());
                let trimmer = PathTrimmer::new(subpath.segments(), 0.001);
                let mut cursor = 0.0;
                for segment in subpath.segments() {
                    let length = segment.arclen(0.001);
                    for station in [cursor + length * 0.5, cursor + length] {
                        let prefix = trimmer.trim(0.0, trimmer.length() - station).unwrap();
                        assert_eq!(
                            trimmer.frame(station).unwrap().unwrap(),
                            CubicBezierSpec::drawable(*prefix.last().unwrap()).endpoint_frame(true)
                        );
                    }
                    cursor += length;
                }
            }
        }
    }

    #[test]
    fn parallel_path_handles_empty_and_stationary_paths() {
        let empty = BezPath::new();
        let mut stationary = BezPath::new();
        stationary.move_to((1.0, 2.0));
        stationary.line_to((1.0, 2.0));
        for path in [empty, stationary] {
            let output = ParallelPathSpec {
                path: path.clone(),
                distance: -0.35,
                start_outset: 0.1,
                end_outset: 0.1,
                accuracy: 0.001,
            }
            .parallel()
            .unwrap();
            assert_eq!(output.path, path);
        }
    }

    #[test]
    fn parallel_path_rejects_removed_options() {
        let path_bytes = encode_cbor(&WirePath::from(&straight_path(3.0))).unwrap();
        let path: ciborium::Value = ciborium::de::from_reader(&path_bytes[..]).unwrap();
        let mut input = std::collections::BTreeMap::from([("path", path)]);
        assert!(curve_parallel_path_bytes(&encode_cbor(&input).unwrap()).is_ok());
        input.insert("optimize", ciborium::Value::Bool(true));
        let error = curve_parallel_path_bytes(&encode_cbor(&input).unwrap()).unwrap_err();
        assert!(error.contains("unknown field `optimize`"), "{error}");
    }

    #[test]
    fn pattern_split_at_cuts_one_patterned_path_into_joined_parts() {
        let spec = |amplitude: f64, split_at: Vec<f64>| PatternPathSpec {
            path: straight_path(3.0),
            pattern: PatternInput::wave(16),
            amplitude,
            wavelength: 0.7,
            phase: 0.0,
            samples_per_period: 16,
            coil_longitudinal_scale: default_coil_longitudinal_scale(),
            anchor_start: true,
            anchor_end: true,
            endpoint_slope: 0.0,
            accuracy: 1e-9,
            split_at,
        };

        let whole = spec(0.2, Vec::new()).patterned().unwrap();
        assert_eq!(whole.parts.len(), 1);
        assert_eq!(whole.parts[0].path, whole.path);

        // Cuts are sorted and clamped; repeated/end cuts retain empty part slots.
        let split = spec(0.2, vec![2.1, 1.3, 1.3, 0.0, 3.0, 7.0])
            .patterned()
            .unwrap();
        assert_eq!(split.parts.len(), 7);
        assert_eq!(split.path, whole.path);
        for index in [0, 2, 5, 6] {
            assert!(split.parts[index].path.elements().is_empty());
        }
        let joined: Vec<_> = split
            .parts
            .iter()
            .flat_map(|part| part.path.segments())
            .collect();
        for pair in joined.windows(2) {
            assert_eq!(pair[0].end(), pair[1].start());
        }
        // Subdividing fitted segments changes commands, not the painted curve.
        let joined = BezPath::from_path_segments(joined.into_iter());
        let length = path_length_value(&whole.path, 1e-9);
        assert!((path_length_value(&joined, 1e-9) - length).abs() < 1e-7);
        let frames = [whole.path.clone(), joined].map(|path| {
            PathFramesSpec {
                path,
                distances: (0..=32).map(|i| length * i as f64 / 32.0).collect(),
                accuracy: 1e-9,
            }
            .frames()
            .unwrap()
        });
        for (original, joined) in frames[0].iter().zip(&frames[1]) {
            assert_point_close(
                original.as_ref().unwrap().point,
                joined.as_ref().unwrap().point,
            );
        }

        // With zero amplitude the cut points sit exactly at the base distances.
        let flat = spec(0.0, vec![1.3, 2.1]).patterned().unwrap();
        let ends: Vec<_> = flat
            .parts
            .iter()
            .map(|part| part.path.segments().last().unwrap().end())
            .collect();
        assert!((ends[0].x - 1.3).abs() < 1e-9 && ends[0].y.abs() < 1e-9);
        assert!((ends[1].x - 2.1).abs() < 1e-9 && ends[1].y.abs() < 1e-9);
    }

    fn outline_spec(path: BezPath, width: f64, cap: &str) -> StrokeOutlineSpec {
        StrokeOutlineSpec {
            path,
            width,
            join: "miter".to_string(),
            miter_limit: 4.0,
            start_cap: cap.to_string(),
            end_cap: cap.to_string(),
            accuracy: 1e-6,
        }
    }

    #[test]
    fn stroke_outline_of_open_line_is_closed_rectangle() {
        use kurbo::Shape;

        let butt = stroke_outline(outline_spec(straight_path(3.0), 0.5, "butt")).unwrap();
        let bbox = butt.path.bounding_box();
        assert_close(bbox.x0, 0.0);
        assert_close(bbox.x1, 3.0);
        assert_close(bbox.y0, -0.25);
        assert_close(bbox.y1, 0.25);
        assert_close(butt.path.area().abs(), 1.5);
        assert_eq!(butt.path.elements().last(), Some(&PathEl::ClosePath));

        let square = stroke_outline(outline_spec(straight_path(3.0), 0.5, "square")).unwrap();
        let bbox = square.path.bounding_box();
        assert_close(bbox.x0, -0.25);
        assert_close(bbox.x1, 3.25);
    }

    #[test]
    fn stroke_outline_miters_sharp_corners() {
        use kurbo::Shape;

        let mut corner = BezPath::new();
        corner.move_to((0.0, 2.0));
        corner.line_to((0.0, 0.0));
        corner.line_to((2.0, 0.0));
        let output = stroke_outline(outline_spec(corner, 0.5, "butt")).unwrap();
        let bbox = output.path.bounding_box();
        assert_close(bbox.x0, -0.25);
        assert_close(bbox.y0, -0.25);
        assert_ne!(output.path.winding(Point::new(-0.2, -0.2)), 0);
        assert_ne!(output.path.winding(Point::new(0.1, 1.0)), 0);
        assert_eq!(output.path.winding(Point::new(1.0, 1.0)), 0);
    }

    #[test]
    fn stroke_outline_rejects_unknown_styles() {
        let mut spec = outline_spec(straight_path(1.0), 0.5, "butt");
        spec.join = "sharp".to_string();
        assert!(stroke_outline(spec).is_err());
        assert!(stroke_outline(outline_spec(straight_path(1.0), 0.5, "flat")).is_err());
        assert!(stroke_outline(outline_spec(straight_path(1.0), 0.0, "butt")).is_err());
    }

    fn straight_curve(length: f64) -> CubicBezierSpec {
        CubicBezierSpec {
            start: point(0.0, 0.0),
            control_start: point(length / 3.0, 0.0),
            control_end: point(2.0 * length / 3.0, 0.0),
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

    #[test]
    fn pattern_splits_preserve_full_fitted_geometry() {
        let carrier = segment_path(CubicBez::new(
            (0.0, 0.0),
            (0.0, 2.0),
            (2.0, 2.0),
            (2.0, 0.0),
        ));
        for pattern in [
            PatternInput::coil(16, 1.4),
            PatternInput::sampled(
                "untapered-coil",
                16,
                |theta| (1.4 * theta.cos(), theta.sin()),
                false,
            ),
            PatternInput::wave(16),
            PatternInput::zigzag(),
        ] {
            let mut spec = PatternPathSpec {
                path: carrier.clone(),
                pattern,
                amplitude: 0.15,
                wavelength: 0.53,
                phase: 0.37,
                samples_per_period: 16,
                coil_longitudinal_scale: 1.4,
                anchor_start: true,
                anchor_end: true,
                endpoint_slope: 1.0,
                split_at: Vec::new(),
                accuracy: 1e-9,
            };
            let whole = spec.clone().patterned().unwrap();
            assert_eq!(whole.parts[0].path, whole.path);
            let length = path_length_value(&carrier, spec.accuracy);
            let distances = PointPattern::from_input(spec.pattern.clone())
                .unwrap()
                .distances(length, spec.wavelength, spec.samples_per_period);
            let cut = distances[4] + 0.37 * (distances[5] - distances[4]);
            let gap_end = distances[10] + 0.61 * (distances[11] - distances[10]);
            spec.split_at = vec![cut, gap_end];
            let split = spec.patterned().unwrap();
            assert_eq!(split.path, whole.path);
            assert_eq!(split.parts.len(), 3);
            let original = whole.path.segments().collect::<Vec<_>>();
            let first = split.parts[0].path.segments().collect::<Vec<_>>();
            let middle = split.parts[1].path.segments().collect::<Vec<_>>();
            let last = split.parts[2].path.segments().collect::<Vec<_>>();
            assert_eq!(first[..4], original[..4]);
            assert_eq!(middle[1..6], original[5..10]);
            assert_eq!(last[1..], original[11..]);
            for (part, full, start, end) in [
                (first[4], original[4], 0.0, 0.37),
                (middle[0], original[4], 0.37, 1.0),
                (middle[6], original[10], 0.0, 0.61),
                (last[0], original[10], 0.61, 1.0),
            ] {
                for t in [0.0, 0.2, 0.5, 0.9, 1.0] {
                    let error = (part.eval(t) - full.eval(start + (end - start) * t)).hypot();
                    assert!(error < 1e-12, "{} split changed by {error}", split.pattern);
                }
            }
        }
    }

    #[test]
    fn pattern_splits_sort_clamp_and_retain_empty_intervals() {
        let mut spec = PatternPathSpec {
            path: straight_path(2.0),
            pattern: PatternInput::coil(16, 1.4),
            amplitude: 0.15,
            wavelength: 0.5,
            phase: 0.0,
            samples_per_period: 16,
            coil_longitudinal_scale: 1.4,
            anchor_start: true,
            anchor_end: true,
            endpoint_slope: 0.0,
            split_at: vec![3.0, 1.0, -1.0, 1.0],
            accuracy: 1e-9,
        };
        let output = spec.clone().patterned().unwrap();
        assert_eq!(output.parts.len(), 5);
        for index in [0, 2, 4] {
            assert!(output.parts[index].path.is_empty());
        }
        let pieces = output.parts[1]
            .path
            .segments()
            .chain(output.parts[3].path.segments());
        assert_eq!(BezPath::from_path_segments(pieces), output.path);

        for path in [BezPath::new(), straight_path(0.0)] {
            spec.path = path;
            let output = spec.clone().patterned().unwrap();
            assert_eq!(output.parts.len(), 5);
            assert!(output.parts.iter().all(|part| part.path.is_empty()));
        }
        for distance in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            spec.split_at = vec![distance];
            assert!(spec.clone().patterned().unwrap_err().contains("split-at"));
        }
    }

    #[test]
    fn pattern_split_distances_accept_integer_and_float_cbor() {
        let bytes = encode_cbor(&CurvePathOutput {
            path: straight_path(2.0),
        })
        .unwrap();
        let default: PatternPathSpec = ciborium::de::from_reader(&bytes[..]).unwrap();
        assert!(default.split_at.is_empty());
        let mut input: std::collections::BTreeMap<String, ciborium::Value> =
            ciborium::de::from_reader(&bytes[..]).unwrap();
        input.insert(
            "split-at".to_string(),
            ciborium::Value::Array(vec![
                ciborium::Value::Integer(1.into()),
                ciborium::Value::Float(1.5),
            ]),
        );
        let bytes = encode_cbor(&input).unwrap();
        let spec: PatternPathSpec = ciborium::de::from_reader(&bytes[..]).unwrap();
        assert_eq!(spec.split_at, [1.0, 1.5]);
    }

    #[test]
    fn wave_pattern_offsets_along_curve_normal() {
        let output = PatternPathSpec {
            path: straight_path(4.0),
            pattern: PatternInput::wave(4),
            amplitude: 0.5,
            wavelength: 4.0,
            phase: 0.0,
            samples_per_period: 4,
            coil_longitudinal_scale: default_coil_longitudinal_scale(),
            anchor_start: default_anchor_endpoint(),
            anchor_end: default_anchor_endpoint(),
            endpoint_slope: 0.0,
            accuracy: 1e-6,
            split_at: Vec::new(),
        }
        .patterned()
        .unwrap();
        let points = path_points(&output.path);
        let segments = output.path.segments().collect::<Vec<_>>();

        assert_eq!(output.pattern, "wave");
        assert!(points.len() >= 5);
        assert!((points[0].y - 0.0).abs() < 1e-9);
        let peak = nearest_x(&points, 1.0);
        assert!((peak.x - 1.0).abs() < 1e-6);
        assert!((peak.y - 0.5).abs() < 1e-6);
        assert_eq!(segments.len(), points.len() - 1);
        assert!(
            segments
                .iter()
                .all(|segment| matches!(segment, PathSeg::Cubic(_)))
        );
    }

    #[test]
    fn zigzag_pattern_samples_corners() {
        let output = PatternPathSpec {
            path: straight_path(4.0),
            pattern: PatternInput::zigzag(),
            amplitude: 1.0,
            wavelength: 4.0,
            phase: 0.0,
            samples_per_period: 2,
            coil_longitudinal_scale: default_coil_longitudinal_scale(),
            anchor_start: default_anchor_endpoint(),
            anchor_end: default_anchor_endpoint(),
            endpoint_slope: 0.0,
            accuracy: 1e-6,
            split_at: Vec::new(),
        }
        .patterned()
        .unwrap();
        let points = path_points(&output.path);

        assert_eq!(output.pattern, "zigzag");
        assert!(points.len() >= 5);
        assert!((nearest_x(&points, 1.0).y - 1.0).abs() < 1e-9);
        assert!((nearest_x(&points, 3.0).y + 1.0).abs() < 1e-9);
        assert!(
            output
                .path
                .segments()
                .all(|segment| matches!(segment, PathSeg::Line(_)))
        );
    }

    #[test]
    fn point_pattern_interpolates_declared_points() {
        let output = PatternPathSpec {
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
            endpoint_slope: 0.0,
            accuracy: 1e-6,
            split_at: Vec::new(),
        }
        .patterned()
        .unwrap();
        let points = path_points(&output.path);

        assert_eq!(output.pattern, "points");
        assert!((nearest_x(&points, 2.0).y - 0.5).abs() < 1e-9);
        assert!(
            output
                .path
                .segments()
                .all(|segment| matches!(segment, PathSeg::Line(_)))
        );
    }

    #[test]
    fn coil_pattern_adds_longitudinal_and_lateral_offsets() {
        let output = PatternPathSpec {
            path: straight_path(4.0),
            pattern: PatternInput::coil(4, 0.5),
            amplitude: 0.5,
            wavelength: 4.0,
            phase: 0.0,
            samples_per_period: 4,
            coil_longitudinal_scale: 0.5,
            anchor_start: false,
            anchor_end: false,
            endpoint_slope: 0.0,
            accuracy: 1e-6,
            split_at: Vec::new(),
        }
        .patterned()
        .unwrap();
        let points = path_points(&output.path);

        assert_eq!(output.pattern, "coil");
        assert!((points[0].x - 0.25).abs() < 1e-9);
        assert!((points[0].y - 0.0).abs() < 1e-9);
        assert!((points[1].x - 1.0).abs() < 1e-6);
        assert!((points[1].y - 0.5).abs() < 1e-6);
        assert_eq!(output.path.segments().count(), points.len() - 1);
    }

    #[test]
    fn coil_default_turns_back_over_the_baseline() {
        let output = PatternPathSpec {
            path: straight_path(2.0),
            pattern: PatternInput::coil(16, default_coil_longitudinal_scale()),
            amplitude: 0.08,
            wavelength: 0.55,
            phase: 0.0,
            samples_per_period: 16,
            coil_longitudinal_scale: default_coil_longitudinal_scale(),
            anchor_start: default_anchor_endpoint(),
            anchor_end: default_anchor_endpoint(),
            endpoint_slope: 0.0,
            accuracy: 1e-6,
            split_at: Vec::new(),
        }
        .patterned()
        .unwrap();
        let points = path_points(&output.path);

        assert!(points.windows(2).any(|window| window[1].x < window[0].x));
    }

    #[test]
    fn coil_endpoint_taper_spans_three_quarters_of_a_turn() {
        let pattern = PointPattern::from_input(PatternInput::coil(16, 1.4)).unwrap();
        for distance in [0.0, 0.0625, 0.125, 0.25, 0.375, 0.5, 1.0] {
            let expected = smoothstep((distance / 0.375_f64).min(1.0));
            assert_eq!(
                pattern_endpoint_envelope(&pattern, distance, 2.0, 0.5, true, true, 0.0),
                expected,
            );
            assert_eq!(
                pattern_endpoint_envelope(&pattern, 2.0 - distance, 2.0, 0.5, true, true, 0.0),
                expected,
            );
        }
        assert_eq!(
            pattern_endpoint_envelope(&pattern, 0.0, 2.0, 0.5, false, true, 0.0),
            1.0,
        );
        assert_eq!(
            pattern_endpoint_envelope(&pattern, 2.0, 2.0, 0.5, true, false, 0.0),
            1.0,
        );
        assert_eq!(
            pattern_endpoint_envelope(&pattern, 0.125, 0.25, 0.5, true, true, 0.0),
            1.0,
        );
    }

    #[test]
    fn coil_endpoint_taper_delays_longitudinal_motion() {
        let output = PatternPathSpec {
            path: straight_path(2.0),
            pattern: PatternInput::coil(16, 1.4),
            amplitude: 0.15,
            wavelength: 0.5,
            phase: 0.0,
            samples_per_period: 16,
            coil_longitudinal_scale: 1.4,
            anchor_start: true,
            anchor_end: true,
            endpoint_slope: 0.0,
            accuracy: 1e-6,
            split_at: Vec::new(),
        }
        .patterned()
        .unwrap();
        let points = path_points(&output.path);
        let offset = 0.15 * 1.4 * smoothstep(1.0 / 6.0).powi(2) * std::f64::consts::FRAC_1_SQRT_2;
        for (index, distance) in [(2, 0.0625), (62, 1.9375)] {
            assert!((points[index].x - (distance + offset)).abs() < 1e-6);
        }
        assert!((points[4].y - 0.15 * smoothstep(1.0 / 3.0)).abs() < 1e-6);
        assert!((points[32].x - (1.0 + 0.15 * 1.4)).abs() < 1e-6);
    }

    #[test]
    fn coil_pattern_anchors_requested_endpoints() {
        let curve = straight_curve(2.0);
        let output = PatternPathSpec {
            path: straight_path(2.0),
            pattern: PatternInput::coil(16, default_coil_longitudinal_scale()),
            amplitude: 0.08,
            wavelength: 0.55,
            phase: 0.0,
            samples_per_period: 16,
            coil_longitudinal_scale: default_coil_longitudinal_scale(),
            anchor_start: true,
            anchor_end: true,
            endpoint_slope: 0.0,
            accuracy: 1e-6,
            split_at: Vec::new(),
        }
        .patterned()
        .unwrap();
        let points = path_points(&output.path);

        assert_eq!(points[0], curve.start);
        assert_eq!(*points.last().unwrap(), curve.end);
        assert!(points[1].y.abs() < 0.08);
    }

    #[test]
    fn sampled_distances_coalesce_terminal_roundoff_but_keep_partial_steps() {
        for grid_end in [0.5_f64, 2.0, 2048.0] {
            let step = grid_end / 4.0;
            for length in [grid_end.next_down(), grid_end, grid_end.next_up()] {
                assert_eq!(
                    sampled_distances(length, grid_end, 4),
                    vec![0.0, step, 2.0 * step, 3.0 * step, length],
                );
            }
            let length = grid_end + step * 0.25;
            assert_eq!(
                sampled_distances(length, grid_end, 4),
                vec![0.0, step, 2.0 * step, 3.0 * step, grid_end, length],
            );
        }
    }

    #[test]
    fn pattern_endpoint_slopes_preserve_anchors_and_flatten_into_the_interior() {
        let pattern = PointPattern::from_input(PatternInput::coil(16, 1.4)).unwrap();
        for slope in [0.0, 1.0, 2.0, 3.0] {
            for length in [2.0, 0.25] {
                let ramp = 0.375_f64.min(length * 0.5);
                for t in [-0.25_f64, 0.0, 0.25, 0.5, 1.0, 1.25] {
                    let distance = t * ramp;
                    let t = t.clamp(0.0, 1.0);
                    let expected = smoothstep(t) + slope * t * (1.0 - t).powi(2);
                    for (distance, start) in [(distance, true), (length - distance, false)] {
                        let envelope = pattern_endpoint_envelope(
                            &pattern, distance, length, 0.5, start, !start, slope,
                        );
                        assert_close(envelope, expected);
                    }
                }
                let step = ramp * 1e-6;
                for distance in [step, length - step] {
                    let envelope = pattern_endpoint_envelope(
                        &pattern, distance, length, 0.5, true, true, slope,
                    );
                    assert_close(envelope / 1e-6, slope);
                }
                for distance in [ramp - step, length - ramp + step] {
                    let envelope = pattern_endpoint_envelope(
                        &pattern, distance, length, 0.5, true, true, slope,
                    );
                    assert_close((1.0 - envelope) / 1e-6, 0.0);
                }
            }
            assert_eq!(
                pattern_endpoint_envelope(&pattern, 0.0, 0.0, 0.5, true, true, slope),
                1.0,
            );
        }
    }

    #[test]
    fn coil_endpoint_slopes_allow_nonzero_transverse_approach() {
        for endpoint_slope in [0.0, 1.0, 2.0, 3.0] {
            let output = PatternPathSpec {
                path: straight_path(2.0),
                pattern: PatternInput::coil(16, 1.4),
                amplitude: 0.15,
                wavelength: 0.5,
                phase: std::f64::consts::FRAC_PI_2,
                samples_per_period: 4096,
                coil_longitudinal_scale: 1.4,
                anchor_start: true,
                anchor_end: true,
                endpoint_slope,
                accuracy: 1e-6,
                split_at: Vec::new(),
            }
            .patterned()
            .unwrap();
            let curves = path_cubics(&output.path);
            let first = curves.first().unwrap();
            let last = curves.last().unwrap();
            assert_eq!(first.start, point(0.0, 0.0));
            assert_eq!(last.end, point(2.0, 0.0));
            let outgoing = Point::from(first.control_start) - Point::from(first.start);
            let incoming = Point::from(last.end) - Point::from(last.control_end);
            let expected = 0.15 * endpoint_slope / 0.375;
            // Endpoint spline tangents use sampled chords, so allow discretization error.
            assert!((outgoing.y / outgoing.x - expected).abs() < 1e-3);
            assert!((incoming.y / incoming.x + expected).abs() < 1e-3);
        }
    }

    #[test]
    fn endpoint_slope_leaves_unanchored_and_unramped_paths_unchanged() {
        for (endpoint_ramp, anchor_start, anchor_end) in [
            (true, false, false),
            (false, false, false),
            (false, true, false),
            (false, false, true),
            (false, true, true),
        ] {
            let mut spec = PatternPathSpec {
                path: straight_path(2.0),
                pattern: PatternInput::sampled(
                    "coil",
                    16,
                    |theta| (1.4 * theta.cos(), theta.sin()),
                    endpoint_ramp,
                ),
                amplitude: 0.15,
                wavelength: 0.5,
                phase: std::f64::consts::FRAC_PI_2,
                samples_per_period: 16,
                coil_longitudinal_scale: 1.4,
                anchor_start,
                anchor_end,
                endpoint_slope: 0.0,
                accuracy: 1e-6,
                split_at: Vec::new(),
            };
            let expected = spec.clone().patterned().unwrap();
            for endpoint_slope in [1.0, 2.0, 3.0] {
                spec.endpoint_slope = endpoint_slope;
                assert_eq!(spec.clone().patterned().unwrap(), expected);
            }
        }
    }

    #[test]
    fn pattern_path_rejects_invalid_endpoint_slopes() {
        for endpoint_slope in [-1.0, 3.01, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let error = PatternPathSpec {
                path: straight_path(2.0),
                pattern: PatternInput::coil(16, 1.4),
                amplitude: 0.15,
                wavelength: 0.5,
                phase: 0.0,
                samples_per_period: 16,
                coil_longitudinal_scale: 1.4,
                anchor_start: true,
                anchor_end: true,
                endpoint_slope,
                accuracy: 1e-6,
                split_at: Vec::new(),
            }
            .patterned()
            .unwrap_err();
            assert!(
                error.contains("endpoint-slope"),
                "{endpoint_slope}: {error}"
            );
        }
    }

    #[test]
    fn pattern_endpoint_slope_deserializes_default_and_kebab_case_numbers() {
        let bytes = encode_cbor(&CurvePathOutput {
            path: straight_path(2.0),
        })
        .unwrap();
        let default: PatternPathSpec = ciborium::de::from_reader(&bytes[..]).unwrap();
        assert_eq!(default.endpoint_slope, 0.0);
        let mut input: std::collections::BTreeMap<String, ciborium::Value> =
            ciborium::de::from_reader(&bytes[..]).unwrap();
        for (value, endpoint_slope) in [
            (ciborium::Value::Float(1.5), 1.5),
            (ciborium::Value::Integer(2.into()), 2.0),
        ] {
            input.insert("endpoint-slope".to_string(), value);
            let bytes = encode_cbor(&input).unwrap();
            let spec: PatternPathSpec = ciborium::de::from_reader(&bytes[..]).unwrap();
            assert_eq!(
                spec,
                PatternPathSpec {
                    endpoint_slope,
                    ..default.clone()
                }
            );
        }
    }

    fn fitted_coil_input() -> FittedCoilInput {
        FittedCoilInput {
            kind: "fitted-coil".to_string(),
            fit_length: 3.0,
            amplitude: 0.12,
            wavelength: 0.55,
            longitudinal_scale: 1.25,
            samples_per_period: 16,
        }
    }

    #[test]
    fn fitted_coil_descriptor_matches_public_points_on_curved_split_paths() {
        let base = segment_path(CubicBez::new(
            (0.0, 0.0),
            (0.0, 2.0),
            (2.0, -1.0),
            (3.0, 0.0),
        ));
        let length = path_length_value(&base, 1e-9);
        for amplitude in [-0.12, -0.0, 0.0, 0.12] {
            for samples_per_period in [-2, 0, 1, 3, 16, 31] {
                let coil = FittedCoilInput {
                    fit_length: length,
                    amplitude,
                    samples_per_period,
                    ..fitted_coil_input()
                };
                let public_bytes =
                    curve_fitted_coil_points_bytes(&encode_cbor(&coil).unwrap()).unwrap();
                let points: Vec<PatternPointInput> =
                    ciborium::de::from_reader(public_bytes.as_slice()).unwrap();
                let sample_count = points.len() - 1;
                assert_eq!(points.first().unwrap().x, 0.0);
                assert_eq!(points.last().unwrap().y, 0.0);
                let descriptor = PatternPathSpec {
                    path: base.clone(),
                    pattern: PatternInput::FittedCoil(coil),
                    amplitude,
                    wavelength: length,
                    phase: 0.0,
                    samples_per_period: 16,
                    coil_longitudinal_scale: 1.25,
                    anchor_start: true,
                    anchor_end: true,
                    endpoint_slope: 0.0,
                    split_at: vec![length * 0.413, length * 0.627],
                    accuracy: 1e-9,
                };
                let explicit = PatternPathSpec {
                    pattern: PatternInput::Points(PointPatternInput {
                        kind: "points".to_string(),
                        name: Some("coil".to_string()),
                        points,
                        interpolation: "smooth".to_string(),
                        endpoint_ramp: false,
                    }),
                    samples_per_period: sample_count,
                    ..descriptor.clone()
                };
                assert_eq!(
                    descriptor.patterned().unwrap(),
                    explicit.patterned().unwrap()
                );
            }
        }
    }

    #[test]
    fn fitted_coil_rejects_invalid_geometry_and_unrepresentable_sample_counts() {
        let good = fitted_coil_input();
        for bad in [
            FittedCoilInput {
                fit_length: 0.0,
                ..good.clone()
            },
            FittedCoilInput {
                fit_length: f64::NAN,
                ..good.clone()
            },
            FittedCoilInput {
                wavelength: -1.0,
                ..good.clone()
            },
            FittedCoilInput {
                wavelength: f64::MIN_POSITIVE,
                ..good.clone()
            },
            FittedCoilInput {
                amplitude: f64::INFINITY,
                ..good.clone()
            },
            FittedCoilInput {
                amplitude: f64::MAX,
                ..good.clone()
            },
            FittedCoilInput {
                longitudinal_scale: -1.0,
                ..good.clone()
            },
            FittedCoilInput {
                longitudinal_scale: f64::NAN,
                ..good.clone()
            },
            FittedCoilInput {
                samples_per_period: i64::MAX,
                ..good.clone()
            },
            FittedCoilInput {
                kind: "unknown".to_string(),
                ..good
            },
        ] {
            assert!(
                bad.points().is_err(),
                "accepted invalid fitted coil: {bad:?}"
            );
        }
    }

    #[test]
    fn fitted_point_pattern_matches_straight_formula_and_endpoints() {
        for (length, wavelength) in [(2.0_f64, 0.4), (2.0, 0.43), (0.1, 1.0)] {
            let periods = (length / wavelength).round().max(1.0) - 0.5;
            let samples_per_period = (periods * 16.0) as usize;
            for amplitude in [-0.15_f64, 0.0, 0.15] {
                let pattern = PatternInput::sampled(
                    "fitted",
                    samples_per_period,
                    |position| {
                        let u = position / std::f64::consts::TAU;
                        let theta = std::f64::consts::PI + periods * position;
                        if u == 0.0 || u == 1.0 {
                            (0.0, 0.0)
                        } else {
                            (
                                amplitude.signum() * 1.4 * length
                                    / (length + 2.0 * amplitude.abs() * 1.4)
                                    * (1.0 + theta.cos() - 2.0 * u),
                                theta.sin(),
                            )
                        }
                    },
                    false,
                );
                let output = PatternPathSpec {
                    path: straight_path(length),
                    pattern,
                    amplitude,
                    wavelength: length,
                    phase: 0.0,
                    samples_per_period,
                    coil_longitudinal_scale: default_coil_longitudinal_scale(),
                    anchor_start: true,
                    anchor_end: true,
                    endpoint_slope: 3.0,
                    accuracy: 1e-9,
                    split_at: Vec::new(),
                }
                .patterned()
                .unwrap();
                let points = path_points(&output.path);
                assert_eq!(output.pattern, "fitted");
                assert_eq!(points.len(), samples_per_period + 1);
                assert_eq!(points[0], point(0.0, 0.0));
                assert_eq!(*points.last().unwrap(), point(length, 0.0));
                for (index, actual) in points.iter().enumerate() {
                    let distance = length * index as f64 / samples_per_period as f64;
                    let theta = std::f64::consts::PI + std::f64::consts::TAU * index as f64 / 16.0;
                    let x = (distance + amplitude.abs() * 1.4 * (theta.cos() + 1.0)) * length
                        / (length + 2.0 * amplitude.abs() * 1.4);
                    assert_point_close(*actual, point(x, amplitude * theta.sin()));
                }
                let curves = path_cubics(&output.path);
                let first = curves.first().unwrap().control_start;
                let last = curves.last().unwrap().control_end;
                assert!(first.x > 0.0 && last.x < length);
                if amplitude != 0.0 {
                    assert!(amplitude * first.y < 0.0 && amplitude * last.y < 0.0);
                }
            }
        }
    }

    #[test]
    fn fitted_point_pattern_uses_local_tangents_and_normals_on_curves() {
        let curve = CubicBez::new((0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0));
        let length = curve.arclen(1e-9);
        for amplitude in [-0.12_f64, 0.0, 0.12] {
            let scale =
                amplitude.signum() * 1.25 * length / (length + 2.0 * amplitude.abs() * 1.25);
            let spec = PatternPathSpec {
                path: segment_path(curve),
                pattern: PatternInput::sampled(
                    "fitted",
                    144,
                    |position| {
                        let u = position / std::f64::consts::TAU;
                        let theta = std::f64::consts::PI + 4.5 * position;
                        if u == 0.0 || u == 1.0 {
                            (0.0, 0.0)
                        } else {
                            (scale * (1.0 + theta.cos() - 2.0 * u), theta.sin())
                        }
                    },
                    false,
                ),
                amplitude,
                wavelength: length,
                phase: 0.0,
                samples_per_period: 144,
                coil_longitudinal_scale: default_coil_longitudinal_scale(),
                anchor_start: true,
                anchor_end: true,
                endpoint_slope: 3.0,
                accuracy: 1e-9,
                split_at: Vec::new(),
            };
            for (anchor_start, anchor_end) in
                [(true, true), (false, false), (false, true), (true, false)]
            {
                let output = PatternPathSpec {
                    anchor_start,
                    anchor_end,
                    ..spec.clone()
                }
                .patterned()
                .unwrap();
                let points = path_points(&output.path);
                assert_eq!(output.pattern, "fitted");
                assert_eq!(points.len(), 145);
                assert_eq!(points[0], point(0.0, 0.0));
                assert_eq!(*points.last().unwrap(), point(2.0, 0.0));
                for (index, actual) in points.iter().enumerate() {
                    let u = index as f64 / 144.0;
                    let theta = std::f64::consts::PI + std::f64::consts::TAU * index as f64 / 32.0;
                    let t = curve.inv_arclen(length * u, 1e-9);
                    let tangent = curve.deriv().eval(t).to_vec2().normalize();
                    let expected = curve.eval(t)
                        + tangent * (amplitude * scale * (1.0 + theta.cos() - 2.0 * u))
                        + Vec2::new(-tangent.y, tangent.x) * (amplitude * theta.sin());
                    assert_point_close(*actual, expected.into());
                }
                let segments = output.path.segments().collect::<Vec<_>>();
                for segment in &segments {
                    let midpoint = segment.eval(0.5);
                    assert!(midpoint.x.is_finite() && midpoint.y.is_finite());
                }
                for pair in segments.windows(2) {
                    assert_eq!(pair[0].end(), pair[1].start());
                    assert!(
                        (path_seg_tangent(&pair[0], 1.0) - path_seg_tangent(&pair[1], 0.0)).hypot()
                            < 1e-8
                    );
                }
            }
        }
    }

    #[test]
    fn point_patterns_handle_short_and_degenerate_paths() {
        let mut move_only = BezPath::new();
        move_only.move_to((3.0, 4.0));
        for path in [
            BezPath::new(),
            move_only,
            straight_path(0.0),
            segment_path(Line::new((1.0, 2.0), (1.0, 2.0))),
            straight_path(f64::EPSILON * 0.5),
            straight_path(1e-8),
        ] {
            let length = path_length_value(&path, 1e-9);
            for amplitude in [-0.12, 0.0, 0.12] {
                let output = PatternPathSpec {
                    path: path.clone(),
                    pattern: PatternInput::sampled(
                        "custom",
                        16,
                        |theta| (theta.cos(), theta.sin()),
                        false,
                    ),
                    amplitude,
                    wavelength: length.max(f64::EPSILON),
                    phase: 0.37,
                    samples_per_period: 16,
                    coil_longitudinal_scale: default_coil_longitudinal_scale(),
                    anchor_start: true,
                    anchor_end: true,
                    endpoint_slope: 3.0,
                    accuracy: 1e-9,
                    split_at: Vec::new(),
                }
                .patterned()
                .unwrap();
                assert_eq!(output.pattern, "custom");
                if length <= f64::EPSILON {
                    assert_eq!(output.path, path);
                } else {
                    let points = path_points(&output.path);
                    assert_eq!(points.len(), 17);
                    assert_eq!(points[0], point(0.0, 0.0));
                    assert_eq!(*points.last().unwrap(), point(length, 0.0));
                }
                for curve in path_cubics(&output.path) {
                    for point in [
                        curve.start,
                        curve.control_start,
                        curve.control_end,
                        curve.end,
                    ] {
                        assert!(point.x.is_finite() && point.y.is_finite());
                    }
                }
            }
        }
    }

    #[test]
    fn hobby_through_returns_two_smooth_cubic_segments() {
        let output = HobbyThroughSpec {
            start: point(0.0, 0.0),
            through: point(1.0, 1.0),
            end: point(2.0, 0.0),
            omega: 1.0,
            accuracy: 1e-6,
        }
        .curve()
        .unwrap();
        let curves = path_cubics(&output.path);

        assert_eq!(curves.len(), 2);
        assert_eq!(curves[0].start, point(0.0, 0.0));
        assert_eq!(curves[0].end, point(1.0, 1.0));
        assert_eq!(curves[1].start, point(1.0, 1.0));
        assert_eq!(curves[1].end, point(2.0, 0.0));

        let incoming = Point::from(curves[0].end) - Point::from(curves[0].control_end);
        let outgoing = Point::from(curves[1].control_start) - Point::from(curves[1].start);
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
        let output = HobbySplineSpec {
            points: points.clone(),
            omega: 1.0,
            accuracy: 1e-3,
        }
        .curve()
        .unwrap();
        let curves = path_cubics(&output.path);

        assert_eq!(curves.len(), points.len() - 1);
        for (index, segment) in curves.iter().enumerate() {
            assert_eq!(segment.start, points[index]);
            assert_eq!(segment.end, points[index + 1]);
        }

        for pair in curves.windows(2) {
            let incoming = Point::from(pair[0].end) - Point::from(pair[0].control_end);
            let outgoing = Point::from(pair[1].control_start) - Point::from(pair[1].start);
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
        let output_bytes = curve_hobby_through_bytes(&bytes).unwrap();
        let output_value: ciborium::Value = ciborium::de::from_reader(&output_bytes[..]).unwrap();
        let output: CurvePathOutput = ciborium::de::from_reader(&output_bytes[..]).unwrap();
        let curves = path_cubics(&output.path);
        let path_value = cbor_map_get(&output_value, "path");
        let elements = match cbor_map_get(path_value, "elements") {
            ciborium::Value::Array(elements) => elements,
            _ => panic!("expected path element array"),
        };
        let first_start = cbor_map_get(&elements[0], "start");

        assert_eq!(cbor_map_keys(&output_value), vec!["path"]);
        assert_cbor_point_tuple(first_start);
        assert_eq!(curves[0].start, point(0.0, 0.0));
        assert_eq!(curves[0].end, point(1.0, 1.0));
        assert_eq!(curves[1].start, point(1.0, 1.0));
        assert_eq!(curves[1].end, point(2.0, 0.0));
        assert_eq!(output.path.segments().count(), 2);
    }
}

#[cfg(test)]
mod region_sampling_tests {
    use super::*;

    fn part(start: f64, end: f64, visible: bool) -> RegionPart {
        RegionPart {
            visible,
            segments: vec![CubicBezierSpec {
                start: CurvePoint { x: start, y: 0.0 },
                end: CurvePoint { x: end, y: 0.0 },
                control_start: CurvePoint {
                    x: start + (end - start) / 3.0,
                    y: 0.0,
                },
                control_end: CurvePoint {
                    x: start + (end - start) * 2.0 / 3.0,
                    y: 0.0,
                },
            }],
        }
    }

    fn request(parts: Vec<RegionPart>) -> RegionSamplesSpec {
        RegionSamplesSpec {
            parts,
            regions: 4,
            unit: 0.0,
            step: 4.0,
            accuracy: 0.001,
        }
    }

    #[test]
    fn hidden_parts_retain_their_arc_length_regions() {
        let samples = request(vec![part(0.0, 4.0, false), part(4.0, 8.0, true)])
            .samples()
            .unwrap();
        assert_eq!(
            samples.iter().map(|group| group.0).collect::<Vec<_>>(),
            vec![2, 3]
        );
        assert_eq!(samples[0].1.first().unwrap().x, 4.0);
        assert!((samples[1].1.last().unwrap().x - 8.0).abs() < 0.001);
    }

    #[test]
    fn both_ends_of_each_region_are_sampled() {
        let samples = request(vec![part(0.0, 4.0, true)]).samples().unwrap();
        assert!(samples.iter().all(|(_, points)| points.len() == 2));
        for adjacent in samples.windows(2) {
            assert!(
                (adjacent[0].1.last().unwrap().x - adjacent[1].1.first().unwrap().x).abs() < 0.001
            );
        }
    }

    #[test]
    fn sampling_density_uses_scale_magnitude() {
        let samples = [10.0, -10.0].map(|unit| {
            let mut input = request(vec![part(0.0, 3.0, true)]);
            input.regions = 1;
            input.unit = unit;
            input.step = 1.0;
            input.samples().unwrap()
        });
        assert_eq!(samples[0], samples[1]);
        assert_eq!(samples[0][0].1.len(), 31);
        for adjacent in samples[0][0].1.windows(2) {
            assert!((adjacent[1].x - adjacent[0].x).abs() * 10.0 <= 1.0 + 1e-12);
        }
    }

    #[test]
    fn empty_and_zero_length_parts_produce_no_samples() {
        assert!(request(vec![]).samples().unwrap().is_empty());
        assert!(
            request(vec![part(0.0, 0.0, true)])
                .samples()
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn invalid_sampling_parameters_are_rejected() {
        let mut zero_regions = request(vec![]);
        zero_regions.regions = 0;
        assert!(zero_regions.samples().is_err());
        let mut zero_step = request(vec![]);
        zero_step.step = 0.0;
        assert!(zero_step.samples().is_err());
        let mut infinite_scale = request(vec![]);
        infinite_scale.unit = f64::INFINITY;
        assert!(infinite_scale.samples().is_err());
    }
}
