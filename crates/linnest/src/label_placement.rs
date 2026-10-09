//! Numeric geometry for the Typst label-placement search.
//!
//! Keep the operation order in sync with CeTZ's matrix/vector arithmetic. Text
//! measurement and candidate enumeration remain in Typst; this batches the
//! per-candidate transforms, collision costs and search without changing positions,
//! costs, or tie-breaking. Typst certifies every acceptance decision by the search.
use serde::{Deserialize, Deserializer, Serialize, de::Error};
use std::{cell::Cell, collections::BTreeMap};

#[derive(Deserialize)]
struct Frame {
    point: [f64; 2],
    tangent: [f64; 2],
}
impl Frame {
    fn deserialize_batch<'de, D: Deserializer<'de>>(de: D) -> Result<Vec<Self>, D::Error> {
        match ciborium::Value::deserialize(de)? {
            ciborium::Value::Bytes(bytes) => {
                ciborium::de::from_reader(bytes.as_slice()).map_err(D::Error::custom)
            }
            value => value.deserialized().map_err(D::Error::custom),
        }
    }
}

#[derive(Deserialize)]
#[serde(transparent)]
struct Transform([[f64; 4]; 4]);

impl Transform {
    fn point(&self, [x, y]: [f64; 2]) -> [f64; 3] {
        if self.0
            == [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        {
            return [x, y, 0.0];
        }
        std::array::from_fn(|i| {
            let [a, b, c, d] = self.0[i];
            a * x + b * y + c * 0.0 + d * 1.0
        })
    }
}

#[derive(Clone, Copy, Deserialize)]
struct Point(
    #[serde(deserialize_with = "crate::deserialize_f64")] f64,
    #[serde(deserialize_with = "crate::deserialize_f64")] f64,
);

impl Point {
    fn lerp(self, other: Self, at: f64) -> Self {
        Self(
            self.0 + (other.0 - self.0) * at,
            self.1 + (other.1 - self.1) * at,
        )
    }

    fn distance(self, other: Self) -> f64 {
        let x = self.0 - other.0;
        let y = self.1 - other.1;
        (x * x + y * y).sqrt()
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct CubicSegment {
    start: Point,
    control_start: Point,
    control_end: Point,
    end: Point,
}

impl CubicSegment {
    fn point(&self, at: f64) -> [f64; 2] {
        // Match Kurvst's Typst cubic-point de Casteljau evaluation exactly.
        let ab = self.start.lerp(self.control_start, at);
        let bc = self.control_start.lerp(self.control_end, at);
        let cd = self.control_end.lerp(self.end, at);
        let abc = ab.lerp(bc, at);
        let bcd = bc.lerp(cd, at);
        let Point(x, y) = abc.lerp(bcd, at);
        [x, y]
    }

    fn lines(&self, transform: &Transform, accuracy: f64) -> Vec<[[f64; 2]; 2]> {
        // Transform before testing flatness, retaining the Typst de Casteljau
        // operation order and depth cap. Push right first to emit left to right.
        let points =
            [self.start, self.control_start, self.control_end, self.end].map(|Point(x, y)| {
                let [x, y, _] = transform.point([x, y]);
                Point(x, y)
            });
        let mut pending = vec![(points, 0)];
        let mut lines = Vec::new();
        while let Some(([a, b, c, d], depth)) = pending.pop() {
            let chord = Point(d.0 - a.0, d.1 - a.1);
            let squared = chord.0 * chord.0 + chord.1 * chord.1;
            let distance = |point: Point| {
                let relative = Point(point.0 - a.0, point.1 - a.1);
                let at = if squared <= 1e-18 {
                    0.0
                } else {
                    ((relative.0 * chord.0 + relative.1 * chord.1) / squared).clamp(0.0, 1.0)
                };
                point.distance(a.lerp(d, at))
            };
            if distance(b).max(distance(c)) <= accuracy || depth >= 10 {
                lines.push([[a.0, a.1], [d.0, d.1]]);
            } else {
                let ab = a.lerp(b, 0.5);
                let bc = b.lerp(c, 0.5);
                let cd = c.lerp(d, 0.5);
                let abc = ab.lerp(bc, 0.5);
                let bcd = bc.lerp(cd, 0.5);
                let mid = abc.lerp(bcd, 0.5);
                pending.push(([mid, bcd, cd, d], depth + 1));
                pending.push(([a, ab, abc, mid], depth + 1));
            }
        }
        lines
    }
}

#[derive(Deserialize)]
struct PathLinesSpec {
    segments: Vec<CubicSegment>,
    transform: Transform,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    accuracy: f64,
}

impl PathLinesSpec {
    fn lines(&self) -> Vec<Vec<[[f64; 2]; 2]>> {
        self.segments
            .iter()
            .map(|segment| segment.lines(&self.transform, self.accuracy))
            .collect()
    }

    fn flattened(
        segments: &[CubicSegment],
        transform: &Transform,
        accuracy: f64,
    ) -> Vec<[Point; 2]> {
        segments
            .iter()
            .flat_map(|segment| segment.lines(transform, accuracy))
            .map(|line| line.map(|[x, y]| Point(x, y)))
            .collect()
    }
}

impl Point {
    // Keep returned records in their original dimensions; only collision math
    // projects canvas corners into the 2D plane.
    fn deserialize_quad<'de, D: Deserializer<'de>>(de: D) -> Result<Option<[Self; 4]>, D::Error> {
        use serde::de::{IgnoredAny, MapAccess, SeqAccess, Visitor};

        #[derive(Deserialize)]
        struct Coordinate(#[serde(deserialize_with = "crate::deserialize_f64")] f64);
        #[derive(Deserialize)]
        #[serde(field_identifier, rename_all = "lowercase")]
        enum Field {
            X,
            Y,
            #[serde(other)]
            Other,
        }
        struct Corner(Point);
        struct CornerVisitor;
        impl<'de> Visitor<'de> for CornerVisitor {
            type Value = Corner;

            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                formatter.write_str("a corner array with at least two coordinates or an x/y map")
            }

            fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Corner, A::Error> {
                let x = seq.next_element::<Coordinate>()?.ok_or_else(|| {
                    A::Error::custom("label corners require at least two coordinates")
                })?;
                let y = seq.next_element::<Coordinate>()?.ok_or_else(|| {
                    A::Error::custom("label corners require at least two coordinates")
                })?;
                while seq.next_element::<IgnoredAny>()?.is_some() {}
                Ok(Corner(Point(x.0, y.0)))
            }

            fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Corner, A::Error> {
                let (mut x, mut y) = (None, None);
                while let Some(field) = map.next_key::<Field>()? {
                    match field {
                        Field::X => x = Some(map.next_value::<Coordinate>()?.0),
                        Field::Y => y = Some(map.next_value::<Coordinate>()?.0),
                        Field::Other => {
                            map.next_value::<IgnoredAny>()?;
                        }
                    }
                }
                Ok(Corner(Point(
                    x.ok_or_else(|| A::Error::missing_field("x"))?,
                    y.ok_or_else(|| A::Error::missing_field("y"))?,
                )))
            }
        }
        impl<'de> Deserialize<'de> for Corner {
            fn deserialize<D: Deserializer<'de>>(de: D) -> Result<Self, D::Error> {
                de.deserialize_any(CornerVisitor)
            }
        }
        Ok(Option::<[Corner; 4]>::deserialize(de)?.map(|corners| corners.map(|corner| corner.0)))
    }

    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0, self.1 - other.1)
    }

    fn dot(self, other: Self) -> f64 {
        self.0 * other.0 + self.1 * other.1
    }

    fn cross(self, other: Self) -> f64 {
        self.0 * other.1 - self.1 * other.0
    }

    fn scaled(self, factor: f64) -> Self {
        Self(self.0 * factor, self.1 * factor)
    }

    fn length(self) -> f64 {
        (self.0 * self.0 + self.1 * self.1).sqrt()
    }
}

#[derive(Deserialize)]
struct AttachmentSpec {
    initial: ciborium::Value,
    lines: Vec<[Point; 2]>,
    corners: [Point; 4],
    frame: Point,
    outward: Point,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    gap: f64,
}

impl AttachmentSpec {
    // Intersect the normal ray with the same finite capsule as the Typst
    // clearance calculation: its endpoint discs and straight rectangle.
    fn capsule(a: Point, b: Point, radius: f64) -> Option<(f64, f64)> {
        if a.0.min(b.0) > radius || a.0.max(b.0) < -radius {
            return None;
        }
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for point in [a, b] {
            if point.0.abs() <= radius {
                let reach = (radius * radius - point.0 * point.0).max(0.0).sqrt();
                low = low.min(point.1 - reach);
                high = high.max(point.1 + reach);
            }
        }
        let delta = b.sub(a);
        let length = delta.length();
        if length > 1e-12 {
            let along = delta.scaled(1.0 / length);
            let across = Point(-along.1, along.0);
            let mut bottom = f64::NEG_INFINITY;
            let mut top = f64::INFINITY;
            for (axis, minimum, maximum) in [(along, 0.0, length), (across, -radius, radius)] {
                let intercept = a.dot(axis);
                let slope = axis.1;
                if slope.abs() <= 1e-12 {
                    if -intercept < minimum || -intercept > maximum {
                        bottom = f64::INFINITY;
                        top = f64::NEG_INFINITY;
                        break;
                    }
                } else {
                    let first = (minimum + intercept) / slope;
                    let last = (maximum + intercept) / slope;
                    bottom = bottom.max(first.min(last));
                    top = top.min(first.max(last));
                }
            }
            if bottom <= top {
                low = low.min(bottom);
                high = high.max(top);
            }
        }
        (low <= high).then_some((low, high))
    }

    fn offset(&self) -> Option<f64> {
        if self.lines.is_empty() {
            return None;
        }
        let scale = self.outward.length();
        if scale <= 1e-9 {
            return None;
        }
        let normal = self.outward.scaled(1.0 / scale);
        let tangent = Point(-normal.1, normal.0);
        let project = |point: Point| Point(point.dot(tangent), point.dot(normal));
        let quad = [
            self.corners[0],
            self.corners[1],
            self.corners[3],
            self.corners[2],
        ]
        .map(project);
        let orientation = quad[1].sub(quad[0]).cross(quad[2].sub(quad[1]));
        let left = quad
            .iter()
            .map(|point| point.0)
            .fold(f64::INFINITY, f64::min);
        let right = quad
            .iter()
            .map(|point| point.0)
            .fold(f64::NEG_INFINITY, f64::max);
        let radius = self.gap * scale;
        let mut intervals = Vec::with_capacity(self.lines.len());
        for &[start, end] in &self.lines {
            let start = project(start.sub(self.frame));
            let end = project(end.sub(self.frame));
            if start.0.min(end.0) > right + radius || start.0.max(end.0) < left - radius {
                continue;
            }
            // Sweep the measured quad along this finite segment, preserving
            // boundary enumeration and the degenerate-quad cases.
            let first = quad.map(|corner| start.sub(corner));
            let last = quad.map(|corner| end.sub(corner));
            let delta = end.sub(start);
            let at_end: [bool; 4] = std::array::from_fn(|i| {
                first[(i + 1) % 4].sub(first[i]).cross(delta) * orientation < 0.0
            });
            let mut low = f64::INFINITY;
            let mut high = f64::NEG_INFINITY;
            for i in 0..4 {
                let next = (i + 1) % 4;
                let points = if at_end[i] { last } else { first };
                for boundary in [
                    Some((points[i], points[next])),
                    (orientation.abs() <= 1e-18).then_some((last[i], last[next])),
                    (at_end[i] != at_end[(i + 3) % 4] || orientation.abs() <= 1e-18)
                        .then_some((first[i], last[i])),
                ]
                .into_iter()
                .flatten()
                {
                    if let Some((bottom, top)) = Self::capsule(boundary.0, boundary.1, radius) {
                        low = low.min(bottom);
                        high = high.max(top);
                    }
                }
            }
            if low <= high {
                intervals.push((low, high));
            }
        }
        // Keep the connected collision interval touching the attachment frame.
        // A remote bend must not pull the label across an already clear gap.
        intervals.sort_by(|first, second| first.0.partial_cmp(&second.0).unwrap());
        let mut reach: f64 = 0.0;
        let mut attached = false;
        for (low, high) in intervals {
            if high < 0.0 {
                continue;
            }
            if low > reach + 1e-10 {
                break;
            }
            attached = true;
            reach = reach.max(high);
        }
        // Independent manual arrow/label shifts may put the finite arrow wholly
        // away from this ray; preserve the deliberate tangential placement.
        if attached { Some(reach / scale) } else { None }
    }
}

pub fn attachment_offset_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: AttachmentSpec = ciborium::de::from_reader(arg)
        .map_err(|error| format!("Invalid finite label attachment geometry: {error}"))?;
    let mut bytes = Vec::new();
    let offset = spec
        .offset()
        .map_or_else(|| spec.initial.clone(), ciborium::Value::Float);
    ciborium::ser::into_writer(&offset, &mut bytes).map_err(|error| error.to_string())?;
    Ok(bytes)
}

#[cfg(test)]
mod attachment_tests {
    use super::*;

    #[test]
    fn finite_attachment_clears_only_the_connected_carrier() {
        let mut spec = AttachmentSpec {
            initial: ciborium::Value::Integer(4.into()),
            lines: vec![[Point(-10.0, 0.0), Point(10.0, 0.0)]],
            corners: [
                Point(-1.0, 0.5),
                Point(1.0, 0.5),
                Point(-1.0, -0.5),
                Point(1.0, -0.5),
            ],
            frame: Point(0.0, 0.0),
            outward: Point(0.0, 1.0),
            gap: 0.2,
        };
        assert_eq!(spec.offset(), Some(0.7));
        spec.lines.push([Point(-10.0, 3.0), Point(10.0, 3.0)]);
        assert_eq!(spec.offset(), Some(0.7));
        spec.lines.remove(0);
        assert_eq!(spec.offset(), None);
        spec.lines = vec![[Point(4.0, 0.0), Point(5.0, 0.0)]];
        assert_eq!(spec.offset(), None);
    }

    #[test]
    fn finite_attachment_preserves_degenerate_and_scaled_clearance() {
        let mut spec = AttachmentSpec {
            initial: ciborium::Value::Integer(4.into()),
            lines: vec![[Point(0.0, 0.0), Point(0.0, 0.0)]],
            corners: [Point(0.0, 0.0); 4],
            frame: Point(0.0, 0.0),
            outward: Point(0.0, 2.0),
            gap: 0.2,
        };
        assert_eq!(spec.offset(), Some(0.2));
        spec.outward = Point(0.0, 0.0);
        assert_eq!(spec.offset(), None);
        spec.lines.clear();
        assert_eq!(spec.offset(), None);
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct ObstacleSpec {
    segments: Vec<CubicSegment>,
    transform: Transform,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    pad_x: f64,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    pad_y: f64,
}

impl ObstacleSpec {
    fn boxes(&self) -> Vec<Bounds> {
        let mut boxes = Vec::with_capacity(self.segments.len() * 12);
        for segment in &self.segments {
            let mut start = self.transform.point(segment.point(0.0));
            for step in 1..=12 {
                let end = self.transform.point(segment.point(f64::from(step) / 12.0));
                boxes.push(Bounds {
                    left: start[0].min(end[0]) - self.pad_x,
                    right: start[0].max(end[0]) + self.pad_x,
                    bottom: start[1].min(end[1]) - self.pad_y,
                    top: start[1].max(end[1]) + self.pad_y,
                });
                start = end;
            }
        }
        boxes
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct CandidateSpec {
    #[serde(default)]
    attachment: Option<CandidateAttachment>,
    #[serde(default)]
    endpoint_direction: Option<Point>,
    #[serde(deserialize_with = "Frame::deserialize_batch")]
    frames: Vec<Frame>,
    positions: Vec<f64>,
    transform: Transform,
    origin: [f64; 3],
    corners: [[f64; 3]; 4],
    total: f64,
    preferred: f64,
    accuracy: f64,
    side: f64,
    preferred_side: f64,
    path_index: usize,
    path_shift: f64,
    gap: f64,
    clear_box: bool,
    fixed: bool,
}

#[derive(Deserialize)]
struct CandidateAttachment {
    carrier: Vec<CubicSegment>,
    #[serde(deserialize_with = "CandidateAttachment::deserialize_paths")]
    paths: Vec<Vec<CubicSegment>>,
}
impl CandidateAttachment {
    fn deserialize_paths<'de, D: Deserializer<'de>>(
        de: D,
    ) -> Result<Vec<Vec<CubicSegment>>, D::Error> {
        match ciborium::Value::deserialize(de)? {
            ciborium::Value::Bytes(bytes) => {
                #[derive(Deserialize)]
                struct Layer {
                    segments: Vec<CubicSegment>,
                }
                let layers: Vec<Layer> =
                    ciborium::de::from_reader(bytes.as_slice()).map_err(D::Error::custom)?;
                Ok(layers.into_iter().map(|layer| layer.segments).collect())
            }
            value => value.deserialized().map_err(D::Error::custom),
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
struct Bounds {
    #[serde(deserialize_with = "crate::deserialize_f64")]
    left: f64,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    right: f64,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    bottom: f64,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    top: f64,
}

impl Bounds {
    fn expanded(self, padding: f64) -> Self {
        Self {
            left: self.left - padding,
            right: self.right + padding,
            bottom: self.bottom - padding,
            top: self.top + padding,
        }
    }

    fn overlap(self, other: Self, padding: f64) -> f64 {
        let width = (self.right.min(other.right) - self.left.max(other.left) + padding).max(0.0);
        if width == 0.0 {
            return 0.0;
        }
        width * (self.top.min(other.top) - self.bottom.max(other.bottom) + padding).max(0.0)
    }

    fn union(self, other: Self) -> Self {
        Self {
            left: self.left.min(other.left),
            right: self.right.max(other.right),
            bottom: self.bottom.min(other.bottom),
            top: self.top.max(other.top),
        }
    }
}

#[derive(Deserialize)]
struct CollisionCandidate {
    bounds: Bounds,
    cost: ciborium::Value,
    #[serde(default, deserialize_with = "Point::deserialize_quad")]
    corners: Option<[Point; 4]>,
    #[serde(default, rename = "arrow-bounds")]
    arrow_bounds: Vec<Bounds>,
}

#[derive(Deserialize)]
struct StrokeLine {
    start: Point,
    end: Point,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    radius: f64,
}

impl StrokeLine {
    /// The stroke envelope tested by `edge_intersection`, with its exact roundings.
    fn extent(&self) -> Bounds {
        Bounds {
            left: self.start.0.min(self.end.0) - self.radius,
            right: self.start.0.max(self.end.0) + self.radius,
            bottom: self.start.1.min(self.end.1) - self.radius,
            top: self.start.1.max(self.end.1) + self.radius,
        }
    }
}

impl Bounds {
    fn is_finite(self) -> bool {
        [self.left, self.right, self.bottom, self.top]
            .iter()
            .all(|value| value.is_finite())
    }
}

/// Uniform grid over stroke extents. A line that passes the bounding test in
/// `edge_intersection` shares a cell with the query, because the cell index is
/// monotone in each coordinate. Ascending indices keep every sum unchanged.
struct LineGrid {
    origin: Point,
    cell: Point,
    side: usize,
    cells: Vec<Vec<usize>>,
    /// Lines with non-finite extents are visited by every query.
    always: Vec<usize>,
    count: usize,
}

impl LineGrid {
    fn new(lines: &[StrokeLine]) -> Self {
        let side = (lines.len() as f64).sqrt().ceil().max(1.0) as usize;
        let total = lines
            .iter()
            .map(StrokeLine::extent)
            .filter(|extent| extent.is_finite())
            .reduce(Bounds::union);
        let (origin, cell) = total.map_or((Point(0.0, 0.0), Point(1.0, 1.0)), |total| {
            let size = |span: f64| (span / side as f64).max(f64::MIN_POSITIVE);
            (
                Point(total.left, total.bottom),
                Point(
                    size(total.right - total.left),
                    size(total.top - total.bottom),
                ),
            )
        });
        let mut grid = Self {
            origin,
            cell,
            side,
            cells: vec![Vec::new(); side * side],
            always: Vec::new(),
            count: lines.len(),
        };
        for (index, line) in lines.iter().enumerate() {
            let extent = line.extent();
            if !extent.is_finite() {
                grid.always.push(index);
                continue;
            }
            let ([left, right], [bottom, top]) = grid.span(extent);
            for row in bottom..=top {
                for column in left..=right {
                    grid.cells[row * side + column].push(index);
                }
            }
        }
        grid
    }

    fn span(&self, bounds: Bounds) -> ([usize; 2], [usize; 2]) {
        let index = |value: f64, origin: f64, size: f64| {
            (((value - origin) / size).floor().max(0.0) as usize).min(self.side - 1)
        };
        (
            [bounds.left, bounds.right].map(|x| index(x, self.origin.0, self.cell.0)),
            [bounds.bottom, bounds.top].map(|y| index(y, self.origin.1, self.cell.1)),
        )
    }

    /// Ascending indices of every line whose extent may meet `bounds`.
    fn near(&self, bounds: Bounds) -> Vec<usize> {
        if !bounds.is_finite() {
            return (0..self.count).collect();
        }
        let ([left, right], [bottom, top]) = self.span(bounds);
        let mut indices = self.always.clone();
        for row in bottom..=top {
            for column in left..=right {
                indices.extend_from_slice(&self.cells[row * self.side + column]);
            }
        }
        indices.sort_unstable();
        indices.dedup();
        indices
    }
}

impl CollisionCandidate {
    fn edge_intersection(&self, lines: &[StrokeLine], grid: &LineGrid) -> Cost {
        let Some(corners) = self.corners.filter(|_| !lines.is_empty()) else {
            return Cost::Integer(0);
        };
        let origin = corners[0];
        let x = corners[1].sub(origin);
        let y = corners[2].sub(origin);
        let determinant = x.cross(y);
        if determinant.abs() < 1e-12 {
            return Cost::Integer(0);
        }
        let width = x.length();
        let height = y.length();
        let local = |point: Point| {
            let delta = point.sub(origin);
            [delta.cross(y) / determinant, x.cross(delta) / determinant]
        };
        let mut length = 0.0;
        for line in grid
            .near(self.bounds)
            .into_iter()
            .map(|index| &lines[index])
        {
            let bounds = self.bounds;
            if line.start.0.max(line.end.0) + line.radius < bounds.left
                || line.start.0.min(line.end.0) - line.radius > bounds.right
                || line.start.1.max(line.end.1) + line.radius < bounds.bottom
                || line.start.1.min(line.end.1) - line.radius > bounds.top
            {
                continue;
            }
            let a = local(line.start);
            let b = local(line.end);
            let mut low: f64 = 0.0;
            let mut high: f64 = 1.0;
            for axis in 0..2 {
                // The inverse transform's row norm also handles sheared canvases.
                let padding =
                    line.radius * (if axis == 0 { height } else { width }) / determinant.abs();
                let delta = b[axis] - a[axis];
                if delta.abs() < 1e-12 {
                    if a[axis] < -padding || a[axis] > 1.0 + padding {
                        high = -1.0;
                        break;
                    }
                } else {
                    let from = (-padding - a[axis]) / delta;
                    let to = (1.0 + padding - a[axis]) / delta;
                    low = low.max(from.min(to));
                    high = high.min(from.max(to));
                }
            }
            length += (high - low).max(0.0) * line.start.distance(line.end);
        }
        Cost::Float(6.0 * length / width.min(height).max(1e-9))
    }
}

/// Native batches retain typed geometry until selection. Expanded custom
/// records remain opaque so their extra fields and corner representations survive.
enum CandidateRecords {
    Expanded(Vec<ciborium::Value>),
    Packed(Vec<PackedCandidate>),
}

struct PackedCandidate {
    geometry: Candidate,
    at: ciborium::Value,
    side: ciborium::Value,
    path_index: usize,
    footprint: ciborium::Value,
    arrow_bounds: Vec<Bounds>,
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct CandidateBatch {
    candidates: ciborium::Value,
    positions: Vec<ciborium::Value>,
    side: ciborium::Value,
    path_index: usize,
}
#[derive(Deserialize)]
struct PackedCandidates {
    batches: Vec<CandidateBatch>,
    #[serde(default)]
    footprints: Option<ciborium::Value>,
    #[serde(default)]
    interleave: bool,
}
impl CandidateRecords {
    fn from_value(value: ciborium::Value) -> Result<Self, String> {
        use ciborium::Value;
        if let Value::Array(records) = value {
            if records
                .iter()
                .any(|record| !matches!(record, Value::Map(_)))
            {
                return Err("label candidates must be dictionaries".to_owned());
            }
            return Ok(Self::Expanded(records));
        }
        let packed: PackedCandidates = value.deserialized().map_err(|e| e.to_string())?;
        let footprints = match packed.footprints {
            Some(Value::Bytes(bytes)) => Some(
                ciborium::de::from_reader::<Value, _>(bytes.as_slice())
                    .map_err(|e| format!("Invalid packed label footprints: {e}"))?,
            ),
            value => value,
        };
        let footprints = footprints
            .map(|value| match value {
                Value::Array(values) => Ok(values),
                _ => Err("label footprints must be an array".to_owned()),
            })
            .transpose()?;
        let expected: usize = packed
            .batches
            .iter()
            .map(|batch| batch.positions.len())
            .sum();
        if footprints
            .as_ref()
            .is_some_and(|values| values.len() != expected)
        {
            return Err("label footprints and candidates must have equal lengths".to_owned());
        }
        let mut footprints = footprints.unwrap_or_default().into_iter();
        let mut batches = Vec::with_capacity(packed.batches.len());
        for batch in packed.batches {
            let candidates: Vec<Candidate> = match batch.candidates {
                Value::Bytes(bytes) => ciborium::de::from_reader(bytes.as_slice())
                    .map_err(|e| format!("Invalid packed label candidates: {e}"))?,
                value => value.deserialized().map_err(|e| e.to_string())?,
            };
            if candidates.len() != batch.positions.len() {
                return Err(
                    "label candidate batch and positions must have equal lengths".to_owned(),
                );
            }
            Cost::from_value(&batch.side)?;
            let records = candidates
                .into_iter()
                .zip(batch.positions)
                .map(|(candidate, at)| {
                    Cost::from_value(&at)?;
                    let footprint = footprints
                        .next()
                        .unwrap_or_else(|| Value::Array(Vec::new()));
                    // Validate while retaining the original footprint numeric types.
                    let arrow_bounds = footprint
                        .deserialized::<Vec<Bounds>>()
                        .map_err(|e| e.to_string())?;
                    Ok(PackedCandidate {
                        geometry: candidate,
                        at,
                        side: batch.side.clone(),
                        path_index: batch.path_index,
                        footprint,
                        arrow_bounds,
                    })
                })
                .collect::<Result<Vec<_>, String>>()?;
            batches.push(records);
        }
        let records = if packed.interleave {
            if batches.len() != 2 || batches[0].len() != batches[1].len() {
                return Err("interleaving requires two equally sized candidate batches".to_owned());
            }
            let right = batches.pop().unwrap();
            let left = batches.pop().unwrap();
            left.into_iter()
                .zip(right)
                .flat_map(|(a, b)| [a, b])
                .collect()
        } else {
            batches.into_iter().flatten().collect()
        };
        Ok(Self::Packed(records))
    }

    fn len(&self) -> usize {
        match self {
            Self::Expanded(records) => records.len(),
            Self::Packed(records) => records.len(),
        }
    }

    fn collision_candidates(&self) -> Result<Vec<CollisionCandidate>, String> {
        match self {
            Self::Expanded(records) => records
                .iter()
                .map(|record| {
                    record
                        .deserialized()
                        .map_err(|error| format!("Invalid label collision candidate: {error}"))
                })
                .collect(),
            Self::Packed(records) => Ok(records
                .iter()
                .map(|record| CollisionCandidate {
                    bounds: record.geometry.bounds,
                    cost: ciborium::Value::Float(record.geometry.cost),
                    corners: Some(record.geometry.corners.map(|[x, y, _]| Point(x, y))),
                    arrow_bounds: record.arrow_bounds.clone(),
                })
                .collect()),
        }
    }

    fn selected(
        &self,
        index: usize,
        cost: Option<&ciborium::Value>,
    ) -> Result<ciborium::Value, String> {
        let mut record = match self {
            Self::Expanded(records) => records[index].clone(),
            Self::Packed(records) => {
                let record = &records[index];
                record.geometry.record(
                    record.at.clone(),
                    record.side.clone(),
                    record.path_index,
                    record.footprint.clone(),
                )?
            }
        };
        if let Some(cost) = cost {
            // A real search has already required the numeric cost field.
            record
                .as_map_mut()
                .unwrap()
                .iter_mut()
                .find(|(key, _)| key.as_text() == Some("cost"))
                .unwrap()
                .1 = cost.clone();
        }
        Ok(record)
    }
}
impl<'de> Deserialize<'de> for CandidateRecords {
    fn deserialize<D: Deserializer<'de>>(de: D) -> Result<Self, D::Error> {
        Self::from_value(ciborium::Value::deserialize(de)?).map_err(D::Error::custom)
    }
}

struct CollisionLabel {
    candidates: Vec<CollisionCandidate>,
    records: CandidateRecords,
    edge: Option<ciborium::Value>,
}
impl<'de> Deserialize<'de> for CollisionLabel {
    fn deserialize<D: Deserializer<'de>>(de: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        struct Input {
            candidates: CandidateRecords,
            edge: Option<ciborium::Value>,
        }
        let Input {
            candidates: records,
            edge,
        } = Input::deserialize(de)?;
        // The all-singleton helper returns records unchanged, even when no
        // collision fields exist. Build the numeric view only for a real search.
        let candidates = Vec::new();
        Ok(Self {
            candidates,
            records,
            edge,
        })
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct Obstacle {
    #[serde(flatten)]
    bounds: Bounds,
    edge: Option<ciborium::Value>,
    #[serde(default)]
    self_loop: bool,
}

impl Obstacle {
    fn belongs_to(&self, label: &CollisionLabel) -> bool {
        use ciborium::Value::{Float, Integer};
        match (&self.edge, &label.edge) {
            (_, None) => false,
            // Typst compares integer and float identifiers numerically.
            (Some(Integer(left)), Some(Float(right)))
            | (Some(Float(right)), Some(Integer(left))) => i128::from(*left) as f64 == *right,
            (left, right) => left == right,
        }
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct CollisionSpec {
    #[serde(default)]
    fixed: bool,
    #[serde(default)]
    edge_lines: Vec<StrokeLine>,
    placements: Vec<CollisionLabel>,
    obstacles: Vec<Obstacle>,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    label_padding: f64,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    obstacle_padding: f64,
    #[serde(default)]
    fixed_arrows: Vec<Bounds>,
    #[serde(default)]
    coordinated: bool,
}

struct CollisionGeometry {
    pair_padding: f64,
    pair_envelopes: Vec<Vec<Option<Bounds>>>,
    pair_offsets: Vec<usize>,
    pair_cost_cache: Vec<Cell<(usize, usize, usize, f64)>>,
    boxes: Vec<Vec<Bounds>>,
    areas: Vec<Vec<f64>>,
    extents: Vec<Bounds>,
    costs: Vec<Vec<ciborium::Value>>,
    arrows: Vec<Vec<Vec<Bounds>>>,
    coordinated: bool,
}

#[derive(Clone, Copy)]
enum Cost {
    Integer(i64),
    Float(f64),
}

impl Cost {
    fn from_value(value: &ciborium::Value) -> Result<Self, String> {
        match value {
            ciborium::Value::Integer(value) => i64::try_from(*value)
                .map(Self::Integer)
                .map_err(|_| "label costs exceed the Typst integer range".to_owned()),
            ciborium::Value::Float(value) => Ok(Self::Float(*value)),
            _ => Err("label candidate costs must be numbers".to_owned()),
        }
    }

    fn value(self) -> ciborium::Value {
        match self {
            Self::Integer(value) => ciborium::Value::Integer(value.into()),
            Self::Float(value) => ciborium::Value::Float(value),
        }
    }

    fn float(self) -> f64 {
        match self {
            Self::Integer(value) => value as f64,
            Self::Float(value) => value,
        }
    }

    fn add(self, other: Self) -> Result<Self, String> {
        match (self, other) {
            (Self::Integer(left), Self::Integer(right)) => left
                .checked_add(right)
                .map(Self::Integer)
                .ok_or_else(|| "label search integer overflow".to_owned()),
            _ => Ok(Self::Float(self.float() + other.float())),
        }
    }

    fn subtract(self, other: Self) -> Result<Self, String> {
        match (self, other) {
            (Self::Integer(left), Self::Integer(right)) => left
                .checked_sub(right)
                .map(Self::Integer)
                .ok_or_else(|| "label search integer overflow".to_owned()),
            _ => Ok(Self::Float(self.float() - other.float())),
        }
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct SearchMath {
    temperatures: Vec<f64>,
    exponentials: Vec<(f64, f64)>,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    pair_padding: f64,
}

#[derive(Serialize)]
struct SearchResult {
    choices: Vec<usize>,
    costs: Vec<ciborium::Value>,
    missing: Vec<ExponentialCheck>,
    selected: Vec<ciborium::Value>,
}

#[derive(Serialize)]
struct ExponentialCheck {
    argument: f64,
    lower: Option<f64>,
    upper: Option<f64>,
    #[serde(skip)]
    probability: f64,
}

impl ExponentialCheck {
    fn new(argument: f64) -> Self {
        Self {
            argument,
            lower: None,
            upper: None,
            probability: argument.exp(),
        }
    }

    fn accept(&mut self, threshold: f64) -> bool {
        let accept = threshold < self.probability;
        if accept {
            self.lower = Some(self.lower.map_or(threshold, |lower| lower.max(threshold)));
        } else {
            self.upper = Some(self.upper.map_or(threshold, |upper| upper.min(threshold)));
        }
        accept
    }
}

impl CollisionGeometry {
    fn pair_cost(&self, i: usize, choice: usize, j: usize, other: usize, padding: f64) -> f64 {
        // Keep the most recent partner choice for each (candidate, label).
        // Costs are directional: reversing the operands changes addition order.
        // Colliding slots only evict work; all three indices certify a hit.
        if padding.to_bits() != self.pair_padding.to_bits() {
            return self.compute_pair_cost(i, choice, j, other, padding);
        }
        let candidate = self.pair_offsets[i] + choice;
        let index = candidate.wrapping_mul(self.boxes.len()).wrapping_add(j)
            & (self.pair_cost_cache.len() - 1);
        let entry = &self.pair_cost_cache[index];
        let cached = entry.get();
        if (cached.0, cached.1, cached.2) == (candidate, j, other) {
            return cached.3;
        }
        // Finite text/arrow envelopes can prove every overlap term is zero.
        // local_cost still adds that zero, preserving numeric types and -0.0.
        let separated = match (
            self.pair_envelopes[i][choice],
            self.pair_envelopes[j][other],
        ) {
            (Some(left), Some(right)) => left.overlap(right, padding) == 0.0,
            _ => false,
        };
        let cost = if separated {
            0.0
        } else {
            self.compute_pair_cost(i, choice, j, other, padding)
        };
        entry.set((candidate, j, other, cost));
        cost
    }

    fn compute_pair_cost(
        &self,
        i: usize,
        choice: usize,
        j: usize,
        other: usize,
        padding: f64,
    ) -> f64 {
        let bounds = self.boxes[i][choice];
        let other_bounds = self.boxes[j][other];
        let area = self.areas[i][choice];
        let other_area = self.areas[j][other];
        let mut cost = 4.0 * bounds.overlap(other_bounds, padding) / area.min(other_area).max(1e-9);
        for &arrow in &self.arrows[i][choice] {
            let arrow_area = (arrow.right - arrow.left) * (arrow.top - arrow.bottom);
            cost +=
                2.0 * arrow.overlap(other_bounds, padding) / arrow_area.min(other_area).max(1e-9);
            for &other_arrow in &self.arrows[j][other] {
                let other_arrow_area =
                    (other_arrow.right - other_arrow.left) * (other_arrow.top - other_arrow.bottom);
                cost += arrow.overlap(other_arrow, padding)
                    / arrow_area.min(other_arrow_area).max(1e-9);
            }
        }
        for &arrow in &self.arrows[j][other] {
            let arrow_area = (arrow.right - arrow.left) * (arrow.top - arrow.bottom);
            cost += 2.0 * bounds.overlap(arrow, padding) / area.min(arrow_area).max(1e-9);
        }
        cost
    }

    fn local_cost(
        &self,
        costs: &[Vec<Cost>],
        i: usize,
        choice: usize,
        choices: &[usize],
        excluding: Option<usize>,
        padding: f64,
    ) -> Cost {
        let mut cost = costs[i][choice];
        for (j, &other) in choices.iter().enumerate() {
            if j == i || Some(j) == excluding {
                continue;
            }
            // Preserve the original label-only search's zero-term elision.
            // Arrows may extend beyond the text's candidate envelope.
            if !self.coordinated && self.extents[i].overlap(self.boxes[j][other], padding) == 0.0 {
                continue;
            }
            cost = Cost::Float(cost.float() + self.pair_cost(i, choice, j, other, padding));
        }
        cost
    }

    fn search(&self, math: &SearchMath) -> Result<SearchResult, String> {
        let count = self.boxes.len();
        if math.temperatures.len() != 84
            || self.areas.len() != count
            || self.extents.len() != count
            || self.costs.len() != count
            || self.arrows.len() != count
            || (0..count).any(|i| {
                self.boxes[i].is_empty()
                    || self.areas[i].len() != self.boxes[i].len()
                    || self.costs[i].len() != self.boxes[i].len()
                    || self.arrows[i].len() != self.boxes[i].len()
            })
        {
            return Err("invalid label search geometry or temperature schedule".to_owned());
        }
        let costs: Vec<Vec<_>> = self
            .costs
            .iter()
            .map(|costs| costs.iter().map(Cost::from_value).collect::<Result<_, _>>())
            .collect::<Result<_, _>>()?;
        let exponentials: BTreeMap<_, _> = math
            .exponentials
            .iter()
            .map(|&(argument, value)| (argument.to_bits(), value))
            .collect();
        let mut missing = Vec::new();
        let mut pending = BTreeMap::new();
        let mut accept = |delta: Cost, temperature: f64, anneal: bool, random: u32| {
            if delta.float() < -1e-12 {
                return true;
            }
            if !anneal {
                return false;
            }
            // Typst's max(0, delta) keeps the integer zero on ties.
            let argument = if delta.float() <= 0.0 {
                0.0
            } else {
                -delta.float() / temperature
            };
            let threshold = f64::from(random) / 4_294_967_296.0;
            match exponentials.get(&argument.to_bits()) {
                Some(&value) => threshold < value,
                None => {
                    let index = *pending.entry(argument.to_bits()).or_insert_with(|| {
                        missing.push(ExponentialCheck::new(argument));
                        missing.len() - 1
                    });
                    // Exp affects only acceptance comparisons. The host verifies
                    // their bounds or supplies exact values and replays.
                    missing[index].accept(threshold)
                }
            }
        };
        let score = |i, choice, choices: &[usize], excluding| {
            self.local_cost(&costs, i, choice, choices, excluding, math.pair_padding)
        };
        let mut choices = vec![0; count];
        let mut best_choices = choices.clone();
        let mut energy = Cost::Integer(0);
        let mut best_energy = energy;
        let mut random = 42_u32;
        for (sweep, &temperature) in math.temperatures.iter().enumerate() {
            if sweep == 80 {
                choices.clone_from(&best_choices);
                energy = best_energy;
            }
            for i in 0..count {
                let candidates = self.boxes[i].len();
                if candidates <= 1 {
                    continue;
                }
                random = random.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                let first = if sweep < 80 {
                    (random >> 16) as usize % candidates
                } else {
                    0
                };
                let last = if sweep < 80 { first + 1 } else { candidates };
                let mut current_score = score(i, choices[i], &choices, None);
                for proposal in first..last {
                    if proposal == choices[i] {
                        continue;
                    }
                    let proposal_score = score(i, proposal, &choices, None);
                    let delta = proposal_score.subtract(current_score)?;
                    random = random.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    if accept(delta, temperature, sweep < 80, random) {
                        choices[i] = proposal;
                        current_score = proposal_score;
                        energy = energy.add(delta)?;
                        if energy.float() < best_energy.float() - 1e-12 {
                            best_energy = energy;
                            best_choices.clone_from(&choices);
                        }
                    }
                }
                if self.coordinated && sweep < 80 && sweep % 4 == 0 && count > 1 {
                    random = random.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    let j = (i + 1 + (random >> 16) as usize % (count - 1)) % count;
                    if self.boxes[j].len() <= 1 {
                        continue;
                    }
                    let mut proposed = choices.clone();
                    for index in [i, j] {
                        random = random.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                        proposed[index] = (random >> 16) as usize % self.boxes[index].len();
                    }
                    let before = score(i, choices[i], &choices, None).add(score(
                        j,
                        choices[j],
                        &choices,
                        Some(i),
                    ))?;
                    let after = score(i, proposed[i], &proposed, None).add(score(
                        j,
                        proposed[j],
                        &proposed,
                        Some(i),
                    ))?;
                    let delta = after.subtract(before)?;
                    random = random.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    if accept(delta, temperature, true, random) {
                        choices = proposed;
                        energy = energy.add(delta)?;
                        if energy.float() < best_energy.float() - 1e-12 {
                            best_energy = energy;
                            best_choices.clone_from(&choices);
                        }
                    }
                }
            }
        }
        if self.coordinated {
            choices.clone_from(&best_choices);
            for _ in 0..2 {
                for i in 0..count {
                    if self.boxes[i].len() <= 1 {
                        continue;
                    }
                    for j in i + 1..count {
                        if self.boxes[j].len() <= 1
                            || self.pair_cost(i, choices[i], j, choices[j], math.pair_padding)
                                <= 1e-12
                        {
                            continue;
                        }
                        // Rank without the partner so both annotations can
                        // leave a crowded basin together. Bound the joint search.
                        let ranked = [i, j].map(|index| {
                            let other = if index == i { j } else { i };
                            let mut ranked: Vec<_> = (0..self.boxes[index].len())
                                .map(|choice| {
                                    (choice, score(index, choice, &choices, Some(other)).float())
                                })
                                .collect();
                            ranked.sort_by(|left, right| {
                                left.1
                                    .partial_cmp(&right.1)
                                    .unwrap_or(std::cmp::Ordering::Equal)
                            });
                            ranked.truncate(6);
                            ranked.push((
                                choices[index],
                                score(index, choices[index], &choices, Some(other)).float(),
                            ));
                            ranked
                        });
                        let mut best = (choices[i], choices[j]);
                        let mut cost = score(i, best.0, &choices, None)
                            .add(score(j, best.1, &choices, Some(i)))?
                            .float();
                        for &(left, left_cost) in &ranked[0] {
                            for &(right, right_cost) in &ranked[1] {
                                let proposed = left_cost
                                    + right_cost
                                    + self.pair_cost(i, left, j, right, math.pair_padding);
                                if proposed < cost - 1e-12 {
                                    cost = proposed;
                                    best = (left, right);
                                }
                            }
                        }
                        choices[i] = best.0;
                        choices[j] = best.1;
                    }
                }
            }
            best_choices = choices;
        }
        Ok(SearchResult {
            costs: best_choices
                .iter()
                .enumerate()
                .map(|(i, &choice)| self.costs[i][choice].clone())
                .collect(),
            choices: best_choices,
            missing,
            selected: Vec::new(),
        })
    }
}

impl CollisionSpec {
    fn geometry(&self, pair_padding: f64) -> Result<CollisionGeometry, String> {
        let mut geometry = CollisionGeometry {
            pair_padding,
            pair_envelopes: Vec::new(),
            pair_offsets: Vec::new(),
            pair_cost_cache: Vec::new(),
            boxes: Vec::with_capacity(self.placements.len()),
            areas: Vec::with_capacity(self.placements.len()),
            extents: Vec::with_capacity(self.placements.len()),
            costs: Vec::with_capacity(self.placements.len()),
            arrows: Vec::with_capacity(self.placements.len()),
            coordinated: self.coordinated,
        };
        let grid = LineGrid::new(&self.edge_lines);
        for label in &self.placements {
            let boxes: Vec<_> = label
                .candidates
                .iter()
                .map(|candidate| candidate.bounds.expanded(self.label_padding))
                .collect();
            let extent = boxes
                .iter()
                .copied()
                .reduce(Bounds::union)
                .ok_or("labels must have at least one candidate")?;
            let relevant: Vec<_> = self
                .obstacles
                .iter()
                .filter(|obstacle| !obstacle.belongs_to(label) || obstacle.self_loop)
                // An obstacle outside the union cannot overlap any candidate.
                // Retain its zero term so the additions below stay unchanged.
                .map(|obstacle| {
                    (self.coordinated
                        || extent.overlap(obstacle.bounds, self.obstacle_padding) != 0.0)
                        .then_some(obstacle)
                })
                .collect();
            let mut areas = Vec::with_capacity(boxes.len());
            let mut costs = Vec::with_capacity(boxes.len());
            for (candidate, bounds) in label.candidates.iter().zip(&boxes) {
                let raw = candidate.bounds;
                let area = (raw.right - raw.left) * (raw.top - raw.bottom);
                areas.push(area);
                let initial = Cost::from_value(&candidate.cost)?
                    .add(candidate.edge_intersection(&self.edge_lines, &grid))?;
                let mut cost = initial.float();
                // Keep zero terms too: adding them converts integer costs to
                // floats. Preserve both obstacle order and numeric types.
                for obstacle in &relevant {
                    let overlap = obstacle.map_or(0.0, |obstacle| {
                        bounds.overlap(obstacle.bounds, self.obstacle_padding)
                    });
                    cost += overlap / area.max(1e-9);
                    for &arrow in &candidate.arrow_bounds {
                        let arrow_area = (arrow.right - arrow.left) * (arrow.top - arrow.bottom);
                        cost += obstacle.map_or(0.0, |obstacle| {
                            arrow.overlap(obstacle.bounds, self.obstacle_padding)
                        }) / arrow_area.max(1e-9);
                    }
                }
                for &fixed in &self.fixed_arrows {
                    let fixed_area = (fixed.right - fixed.left) * (fixed.top - fixed.bottom);
                    cost +=
                        2.0 * bounds.overlap(fixed, pair_padding) / area.min(fixed_area).max(1e-9);
                    for &arrow in &candidate.arrow_bounds {
                        let arrow_area = (arrow.right - arrow.left) * (arrow.top - arrow.bottom);
                        cost += arrow.overlap(fixed, pair_padding)
                            / arrow_area.min(fixed_area).max(1e-9);
                    }
                }
                costs.push(if relevant.is_empty() && self.fixed_arrows.is_empty() {
                    initial.value()
                } else {
                    ciborium::Value::Float(cost)
                });
            }
            geometry.arrows.push(
                label
                    .candidates
                    .iter()
                    .map(|candidate| candidate.arrow_bounds.clone())
                    .collect(),
            );
            geometry.boxes.push(boxes);
            geometry.areas.push(areas);
            geometry.extents.push(extent);
            geometry.costs.push(costs);
        }
        // Build envelopes only where the original area/division arithmetic is
        // finite. Nonfinite geometry keeps the original per-piece evaluation.
        geometry.pair_envelopes = geometry
            .boxes
            .iter()
            .enumerate()
            .map(|(i, boxes)| {
                boxes
                    .iter()
                    .enumerate()
                    .map(|(c, bounds)| {
                        let finite = |b: &Bounds| {
                            [b.left, b.right, b.bottom, b.top]
                                .iter()
                                .all(|v| v.is_finite())
                                && ((b.right - b.left) * (b.top - b.bottom)).is_finite()
                        };
                        (finite(bounds)
                            && geometry.areas[i][c].is_finite()
                            && geometry.arrows[i][c].iter().all(finite))
                        .then(|| {
                            geometry.arrows[i][c]
                                .iter()
                                .fold(*bounds, |all, &arrow| all.union(arrow))
                        })
                    })
                    .collect()
            })
            .collect();
        let mut candidates = 0;
        geometry.pair_offsets = geometry
            .boxes
            .iter()
            .map(|boxes| {
                let offset = candidates;
                candidates += boxes.len();
                offset
            })
            .collect();
        // Bound cache memory independently of graph/candidate counts. The
        // wrapping slot calculation remains safe because hits compare full keys.
        let capacity = candidates
            .saturating_mul(geometry.boxes.len())
            .clamp(1, 65_536)
            .next_power_of_two();
        geometry.pair_cost_cache = (0..capacity)
            .map(|_| Cell::new((usize::MAX, usize::MAX, usize::MAX, 0.0)))
            .collect();
        Ok(geometry)
    }
}

#[derive(Deserialize, Serialize)]
#[serde(rename_all = "kebab-case")]
struct Candidate {
    position: [f64; 2],
    normal: [f64; 2],
    outward: [f64; 3],
    nearest: f64,
    corners: [[f64; 3]; 4],
    at: f64,
    side: f64,
    path_index: usize,
    path_shift: f64,
    bounds: Bounds,
    cost: f64,
}

impl Candidate {
    fn record(
        &self,
        at: ciborium::Value,
        side: ciborium::Value,
        path_index: usize,
        arrow_bounds: ciborium::Value,
    ) -> Result<ciborium::Value, String> {
        // This is the public annotation record, in the same field order as its
        // Typst constructor. Candidate-only clearance intermediates stay native.
        #[derive(Serialize)]
        #[serde(rename_all = "kebab-case")]
        struct Selection {
            position: [f64; 2],
            at: ciborium::Value,
            side: ciborium::Value,
            path_index: usize,
            path_shift: f64,
            corners: [[f64; 3]; 4],
            bounds: Bounds,
            cost: f64,
            arrow_bounds: ciborium::Value,
        }
        ciborium::Value::serialized(&Selection {
            position: self.position,
            at,
            side,
            path_index,
            path_shift: self.path_shift,
            corners: self.corners,
            bounds: self.bounds,
            cost: self.cost,
            arrow_bounds,
        })
        .map_err(|e| e.to_string())
    }
}

impl CandidateSpec {
    fn candidates(&self) -> Result<Vec<Candidate>, String> {
        if self.frames.len() != self.positions.len() {
            return Err("label frames and positions must have equal lengths".to_owned());
        }
        if self
            .attachment
            .as_ref()
            .is_some_and(|attachment| attachment.paths.len() != self.frames.len())
        {
            return Err("label attachment paths and frames must have equal lengths".to_owned());
        }
        let attachment = self
            .attachment
            .as_ref()
            .filter(|_| self.clear_box && !self.fixed);
        let carrier_lines = attachment.map_or_else(Vec::new, |attachment| {
            PathLinesSpec::flattened(&attachment.carrier, &self.transform, self.accuracy)
        });
        Ok(self
            .positions
            .iter()
            .zip(&self.frames)
            .enumerate()
            .map(|(index, (&at, frame))| {
                let normal = if let Some(direction) = self.endpoint_direction {
                    let Point(x, y) = direction.scaled(1.0 / direction.length());
                    [x, y]
                } else {
                    let [tx, ty] = frame.tangent;
                    let length = (tx * tx + ty * ty).sqrt();
                    if length <= 1e-9 {
                        [0.0, self.side]
                    } else {
                        let scale = self.side / length;
                        [-ty * scale, tx * scale]
                    }
                };
                let transformed = self.transform.point(normal);
                let outward: [f64; 3] = std::array::from_fn(|i| transformed[i] - self.origin[i]);
                let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
                let normal_squared = dot(outward, outward);
                let nearest = if (self.clear_box || self.endpoint_direction.is_some())
                    && normal_squared > 1e-18
                {
                    self.corners
                        .iter()
                        .map(|&corner| dot(corner, outward))
                        .fold(f64::INFINITY, f64::min)
                        / normal_squared
                } else {
                    0.0
                };
                let initial = self.gap - nearest;
                let offset = attachment.map_or(initial, |attachment| {
                    // The finite arrow precedes its physical carrier, exactly as
                    // in the scalar attachment boundary. Flatten the carrier once
                    // per label batch instead of encoding its lines per candidate.
                    let mut lines = PathLinesSpec::flattened(
                        &attachment.paths[index],
                        &self.transform,
                        self.accuracy,
                    );
                    lines.extend_from_slice(&carrier_lines);
                    let [x, y, _] = self.transform.point(frame.point);
                    AttachmentSpec {
                        initial: ciborium::Value::Float(initial),
                        lines,
                        corners: self.corners.map(|[x, y, _]| Point(x, y)),
                        frame: Point(x, y),
                        outward: Point(outward[0], outward[1]),
                        gap: self.gap,
                    }
                    .offset()
                    .unwrap_or(initial)
                });
                let position = if self.fixed && self.endpoint_direction.is_none() {
                    frame.point
                } else {
                    std::array::from_fn(|i| frame.point[i] + normal[i] * offset)
                };
                let center = self.transform.point(position);
                let bounds = self
                    .corners
                    .map(|corner| std::array::from_fn(|i| center[i] + corner[i]));
                let relative = (at - self.preferred) / self.total;
                Candidate {
                    position,
                    normal,
                    outward,
                    nearest,
                    corners: bounds,
                    at,
                    side: self.side,
                    path_index: self.path_index,
                    path_shift: at - self.total / 2.0 + self.path_shift,
                    bounds: Bounds {
                        left: bounds.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min),
                        right: bounds
                            .iter()
                            .map(|p| p[0])
                            .fold(f64::NEG_INFINITY, f64::max),
                        bottom: bounds.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min),
                        top: bounds
                            .iter()
                            .map(|p| p[1])
                            .fold(f64::NEG_INFINITY, f64::max),
                    },
                    cost: (if self.total <= self.accuracy {
                        0.0
                    } else {
                        0.002 * relative.powi(2)
                    }) + if self.side == self.preferred_side {
                        0.0
                    } else {
                        0.00002
                    },
                }
            })
            .collect())
    }
}

pub fn candidates_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: CandidateSpec = ciborium::de::from_reader(arg)
        .map_err(|error| format!("Invalid label candidate geometry: {error}"))?;
    crate::graph_api::encode_cbor(&spec.candidates()?)
}

pub fn first_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let records: CandidateRecords = ciborium::de::from_reader(arg)
        .map_err(|error| format!("Invalid label candidates: {error}"))?;
    if records.len() == 0 {
        return Err("labels must have at least one candidate".to_owned());
    }
    crate::graph_api::encode_cbor(&records.selected(0, None)?)
}

pub fn search_bytes(arg: &[u8], math: &[u8]) -> Result<Vec<u8>, String> {
    let mut spec: CollisionSpec = ciborium::de::from_reader(arg)
        .map_err(|error| format!("Invalid label collision geometry: {error}"))?;
    let math: SearchMath = ciborium::de::from_reader(math)
        .map_err(|error| format!("Invalid label search math: {error}"))?;
    let result = if spec.fixed {
        if spec.placements.iter().any(|label| label.records.len() != 1) {
            return Err(
                "fixed label selection requires exactly one candidate per label".to_owned(),
            );
        }
        let selected: Vec<_> = spec
            .placements
            .iter()
            .map(|label| label.records.selected(0, None))
            .collect::<Result<_, _>>()?;
        SearchResult {
            choices: vec![0; spec.placements.len()],
            costs: selected
                .iter()
                .map(|record| {
                    record
                        .as_map()
                        .unwrap()
                        .iter()
                        .find(|(key, _)| key.as_text() == Some("cost"))
                        .map(|(_, value)| value.clone())
                        .unwrap_or_else(|| ciborium::Value::Integer(0.into()))
                })
                .collect(),
            missing: Vec::new(),
            selected,
        }
    } else {
        for label in &mut spec.placements {
            label.candidates = label.records.collision_candidates()?;
        }
        let mut result = spec.geometry(math.pair_padding)?.search(&math)?;
        result.selected = spec
            .placements
            .iter()
            .zip(&result.choices)
            .zip(&result.costs)
            .map(|((label, &choice), cost)| label.records.selected(choice, Some(cost)))
            .collect::<Result<_, _>>()?;
        result
    };
    crate::graph_api::encode_cbor(&result)
}

pub fn obstacle_boxes_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: ObstacleSpec = ciborium::de::from_reader(arg)
        .map_err(|error| format!("Invalid label obstacle geometry: {error}"))?;
    crate::graph_api::encode_cbor(&spec.boxes())
}

pub fn path_lines_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: PathLinesSpec = ciborium::de::from_reader(arg)
        .map_err(|error| format!("Invalid label path geometry: {error}"))?;
    crate::graph_api::encode_cbor(&spec.lines())
}

/// Painted CeTZ strokes, flattened to the collision lines of `CollisionSpec`.
#[derive(Deserialize)]
struct StrokeLinesSpec {
    strokes: Vec<PaintedStroke>,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    accuracy: f64,
}

/// One drawable: CeTZ subpaths `(origin, closed, commands)`, whose commands
/// are `("l", ..ends)` or `("c", control-start, control-end, end)`.
#[derive(Deserialize)]
struct PaintedStroke {
    #[serde(deserialize_with = "crate::deserialize_f64")]
    radius: f64,
    segments: Vec<(ciborium::Value, bool, Vec<Vec<ciborium::Value>>)>,
}

#[derive(Serialize)]
struct CollisionLine {
    start: [f64; 2],
    end: [f64; 2],
    radius: f64,
}

impl StrokeLinesSpec {
    fn vertex(value: &ciborium::Value) -> Result<Vec<f64>, String> {
        let coordinates = value
            .as_array()
            .filter(|coordinates| coordinates.len() >= 2)
            .ok_or("painted vertices need at least two coordinates")?;
        coordinates
            .iter()
            .map(|coordinate| match coordinate {
                ciborium::Value::Float(value) => Ok(*value),
                ciborium::Value::Integer(value) => Ok(i128::from(*value) as f64),
                _ => Err("painted coordinates must be numbers".to_owned()),
            })
            .collect()
    }

    /// Lines in drawable, subpath and command order, as the former Typst
    /// traversal emitted them; cubics use the same identity-transform flattening.
    fn lines(&self) -> Result<Vec<CollisionLine>, String> {
        if !self.accuracy.is_finite() || self.accuracy <= 0.0 {
            return Err("painted stroke accuracy must be finite and positive".to_owned());
        }
        let identity = Transform(std::array::from_fn(|row| {
            std::array::from_fn(|column| if row == column { 1.0 } else { 0.0 })
        }));
        let xy = |vertex: &[f64]| Point(vertex[0], vertex[1]);
        let mut lines = Vec::new();
        for stroke in &self.strokes {
            if !stroke.radius.is_finite() || stroke.radius < 0.0 {
                return Err("painted stroke radius must be finite and nonnegative".to_owned());
            }
            let mut push = |start: Point, end: Point| {
                lines.push(CollisionLine {
                    start: [start.0, start.1],
                    end: [end.0, end.1],
                    radius: stroke.radius,
                })
            };
            for (origin, closed, commands) in &stroke.segments {
                let origin = Self::vertex(origin)?;
                let mut start = origin.clone();
                for command in commands {
                    let kind = command.first().and_then(ciborium::Value::as_text);
                    if kind == Some("l") {
                        for end in &command[1..] {
                            let end = Self::vertex(end)?;
                            push(xy(&start), xy(&end));
                            start = end;
                        }
                    } else if kind == Some("c") {
                        if command.len() != 4 {
                            return Err("painted cubics need two controls and an end".to_owned());
                        }
                        let segment = CubicSegment {
                            start: xy(&start),
                            control_start: xy(&Self::vertex(&command[1])?),
                            control_end: xy(&Self::vertex(&command[2])?),
                            end: xy(&Self::vertex(&command[3])?),
                        };
                        for [a, b] in segment.lines(&identity, self.accuracy) {
                            push(Point(a[0], a[1]), Point(b[0], b[1]));
                        }
                        start = Self::vertex(command.last().expect("checked length"))?;
                    } else {
                        return Err(format!("unsupported painted stroke command: {kind:?}"));
                    }
                }
                if *closed && start != origin {
                    push(xy(&start), xy(&origin));
                }
            }
        }
        Ok(lines)
    }
}

pub fn stroke_lines_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: StrokeLinesSpec = ciborium::de::from_reader(arg)
        .map_err(|error| format!("Invalid painted stroke geometry: {error}"))?;
    crate::graph_api::encode_cbor(&spec.lines()?)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wire_value(value: serde_json::Value) -> ciborium::Value {
        ciborium::Value::serialized(&value).unwrap()
    }
    fn wire_field<'a>(value: &'a mut ciborium::Value, name: &str) -> &'a mut ciborium::Value {
        &mut value
            .as_map_mut()
            .unwrap()
            .iter_mut()
            .find(|(key, _)| key.as_text() == Some(name))
            .unwrap()
            .1
    }
    fn candidate_wire() -> ciborium::Value {
        wire_value(serde_json::json!({
            "frames": [{"point":[0.0,0.0],"tangent":[1.0,0.0]}, {"point":[1.0,0.0],"tangent":[1.0,0.0]}],
            "positions": [0.0,1.0],
            "transform": [[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]],
            "origin": [0.0,0.0,0.0],
            "corners": [[-0.2,0.1,0.7],[0.2,0.1,0.7],[-0.2,-0.1,0.7],[0.2,-0.1,0.7]],
            "total":2.0,"preferred":1.0,"accuracy":0.01,"side":1.0,"preferred-side":1.0,
            "path-index":0,"path-shift":0.0,"gap":0.2,"clear-box":true,"fixed":false,
            "attachment":{"carrier":[],"paths":[[],[]]}
        }))
    }
    fn packed_candidate_wire() -> ciborium::Value {
        use ciborium::Value;
        let bytes =
            candidates_bytes(&crate::graph_api::encode_cbor(&candidate_wire()).unwrap()).unwrap();
        Value::Map(vec![
            (
                "batches".into(),
                Value::Array(
                    (0..2)
                        .map(|index| {
                            Value::Map(vec![
                                ("candidates".into(), Value::Bytes(bytes.clone())),
                                (
                                    "positions".into(),
                                    Value::Array(vec![
                                        Value::Integer(0.into()),
                                        Value::Float(-0.0),
                                    ]),
                                ),
                                (
                                    "side".into(),
                                    if index == 0 {
                                        Value::Integer(1.into())
                                    } else {
                                        Value::Float(-1.0)
                                    },
                                ),
                                ("path-index".into(), Value::Integer((index + 7).into())),
                            ])
                        })
                        .collect(),
                ),
            ),
            ("footprints".into(), Value::Null),
            ("interleave".into(), Value::Bool(true)),
        ])
    }
    fn candidate_records(value: ciborium::Value) -> Vec<ciborium::Value> {
        let records = CandidateRecords::from_value(value).unwrap();
        (0..records.len())
            .map(|index| records.selected(index, None).unwrap())
            .collect()
    }
    fn collision_wire(candidates: ciborium::Value) -> ciborium::Value {
        let mut value = wire_value(serde_json::json!({
            "placements":[{"candidates":[],"edge":null}], "obstacles":[],
            "label-padding":0.0,"obstacle-padding":0.0,"fixed-arrows":[],
            "coordinated":true,"edge-lines":[],"fixed":false
        }));
        let placement = &mut wire_field(&mut value, "placements").as_array_mut().unwrap()[0];
        *wire_field(placement, "candidates") = candidates;
        value
    }
    #[test]
    fn candidate_frames_and_attachment_layers_accept_opaque_packets() {
        use ciborium::Value;
        let mut structured = candidate_wire();
        let segment = wire_value(serde_json::json!({
            "start":[-1,0],"control-start":[0,0],"control-end":[1,0],"end":[2,0]
        }));
        let attachment = wire_field(&mut structured, "attachment");
        *wire_field(attachment, "paths") =
            Value::Array(vec![Value::Array(vec![segment]), Value::Array(Vec::new())]);
        let encode = |value: &Value| crate::graph_api::encode_cbor(value).unwrap();
        let expected = candidates_bytes(&encode(&structured)).unwrap();
        let mut packed = structured;
        let frames = wire_field(&mut packed, "frames");
        *frames = Value::Bytes(encode(frames));
        let paths = wire_field(wire_field(&mut packed, "attachment"), "paths");
        let layers = Value::Array(
            paths
                .as_array()
                .unwrap()
                .iter()
                .cloned()
                .map(|segments| {
                    Value::Map(vec![
                        ("path".into(), Value::Null),
                        ("segments".into(), segments),
                    ])
                })
                .collect(),
        );
        *paths = Value::Bytes(encode(&layers));
        assert_eq!(candidates_bytes(&encode(&packed)).unwrap(), expected);
    }
    #[test]
    fn packed_records_preserve_numeric_types_footprints_and_interleaving() {
        use ciborium::Value;
        let mut packed = packed_candidate_wire();
        let footprints = wire_value(serde_json::json!(
            (0..4)
                .map(|i| vec![serde_json::json!({
                    "left":i,"right":i+1,"bottom":0,"top":1
                })])
                .collect::<Vec<_>>()
        ));
        *wire_field(&mut packed, "footprints") = footprints.clone();
        let mut records = candidate_records(packed.clone());
        let expected_sides = [
            Value::Integer(1.into()),
            Value::Float(-1.0),
            Value::Integer(1.into()),
            Value::Float(-1.0),
        ];
        for (i, record) in records.iter_mut().enumerate() {
            assert_eq!(wire_field(record, "side"), &expected_sides[i]);
            assert_eq!(
                wire_field(record, "path-index"),
                &Value::Integer((7 + i % 2).into())
            );
            assert_eq!(
                wire_field(record, "arrow-bounds"),
                &footprints.as_array().unwrap()[[0, 2, 1, 3][i]]
            );
            let at = wire_field(record, "at");
            if i < 2 {
                assert_eq!(at, &Value::Integer(0.into()));
            } else {
                assert_eq!(at.as_float().unwrap().to_bits(), (-0.0_f64).to_bits());
            }
            assert!(
                wire_field(record, "corners")
                    .as_array()
                    .unwrap()
                    .iter()
                    .all(|p| p.as_array().unwrap().len() == 3)
            );
        }
        let expected = crate::graph_api::encode_cbor(&records[0]).unwrap();
        assert_eq!(
            first_bytes(&crate::graph_api::encode_cbor(&packed).unwrap()).unwrap(),
            expected
        );
        *wire_field(&mut packed, "footprints") =
            Value::Bytes(crate::graph_api::encode_cbor(&footprints).unwrap());
        assert_eq!(
            crate::graph_api::encode_cbor(&candidate_records(packed)).unwrap(),
            crate::graph_api::encode_cbor(&records).unwrap()
        );
    }
    #[test]
    fn packed_and_structured_search_keep_identical_selected_records() {
        use ciborium::Value;
        let packed = packed_candidate_wire();
        let records = candidate_records(packed.clone());
        let math = crate::graph_api::encode_cbor(&serde_json::json!({
            "temperatures":(0..84).map(|s|0.15*0.9_f64.powi(s)).collect::<Vec<_>>(),
            "exponentials":[],"pair-padding":0.35
        }))
        .unwrap();
        let request =
            |candidates| crate::graph_api::encode_cbor(&collision_wire(candidates)).unwrap();
        let expected = search_bytes(&request(Value::Array(records.clone())), &math).unwrap();
        assert_eq!(search_bytes(&request(packed), &math).unwrap(), expected);
        let mut custom = records;
        for record in &mut custom {
            record.as_map_mut().unwrap().push((
                "custom-key".into(),
                Value::Integer(9007199254740993_i64.into()),
            ));
        }
        let reply = search_bytes(&request(Value::Array(custom)), &math).unwrap();
        let mut reply: Value = ciborium::de::from_reader(reply.as_slice()).unwrap();
        let selected = &mut wire_field(&mut reply, "selected").as_array_mut().unwrap()[0];
        assert_eq!(
            wire_field(selected, "custom-key"),
            &Value::Integer(9007199254740993_i64.into())
        );
        assert!(
            wire_field(selected, "corners")
                .as_array()
                .unwrap()
                .iter()
                .all(|p| p.as_array().unwrap()[2] == Value::Float(0.7))
        );
    }
    #[test]
    fn ordinary_search_projects_corner_maps_and_extra_coordinates_without_rewriting_them() {
        use ciborium::Value;
        let mut reference = candidate_records(packed_candidate_wire());
        for record in &mut reference {
            for corner in wire_field(record, "corners").as_array_mut().unwrap() {
                corner.as_array_mut().unwrap().truncate(2);
            }
        }
        let math = crate::graph_api::encode_cbor(&serde_json::json!({
            "temperatures":(0..84).map(|s|0.15*0.9_f64.powi(s)).collect::<Vec<_>>(),
            "exponentials":[],"pair-padding":0.35
        }))
        .unwrap();
        let request = |records| {
            let mut geometry = collision_wire(Value::Array(records));
            *wire_field(&mut geometry, "edge-lines") = wire_value(serde_json::json!([
                {"start":[-1,0.2],"end":[2,0.2],"radius":0.015}
            ]));
            crate::graph_api::encode_cbor(&geometry).unwrap()
        };
        let expected = search_bytes(&request(reference.clone()), &math).unwrap();
        let mut expected_value: Value = ciborium::de::from_reader(expected.as_slice()).unwrap();
        let expected_corners = wire_field(
            &mut wire_field(&mut expected_value, "selected")
                .as_array_mut()
                .unwrap()[0],
            "corners",
        )
        .clone();
        for map_corners in [false, true] {
            let mut records = reference.clone();
            for record in &mut records {
                for corner in wire_field(record, "corners").as_array_mut().unwrap() {
                    let coordinates = corner.as_array().unwrap();
                    // Typst's _point only reads x/y, even when other components
                    // are nonnumeric. Preserve those components in the result.
                    *corner = if map_corners {
                        Value::Map(vec![
                            ("x".into(), coordinates[0].clone()),
                            ("y".into(), coordinates[1].clone()),
                            ("ignored".into(), Value::Text("kept".into())),
                        ])
                    } else {
                        Value::Array(vec![
                            coordinates[0].clone(),
                            coordinates[1].clone(),
                            Value::Float(0.7),
                            Value::Text("kept".into()),
                        ])
                    };
                }
            }
            let actual = search_bytes(&request(records.clone()), &math).unwrap();
            let mut actual: Value = ciborium::de::from_reader(actual.as_slice()).unwrap();
            let chosen = usize::try_from(
                wire_field(&mut actual, "choices").as_array().unwrap()[0]
                    .as_integer()
                    .unwrap(),
            )
            .unwrap();
            let selected = &mut wire_field(&mut actual, "selected").as_array_mut().unwrap()[0];
            assert_eq!(
                wire_field(selected, "corners"),
                wire_field(&mut records[chosen], "corners")
            );
            *wire_field(selected, "corners") = expected_corners.clone();
            assert_eq!(crate::graph_api::encode_cbor(&actual).unwrap(), expected);
        }
    }
    #[test]
    fn singleton_selection_preserves_arbitrary_records_without_scoring() {
        use ciborium::Value;
        for record in [
            Value::Map(vec![("custom".into(), Value::Integer(7.into()))]),
            Value::Map(vec![
                ("cost".into(), Value::Float(-0.0)),
                ("custom".into(), Value::Text("kept".into())),
            ]),
            Value::Map(vec![("cost".into(), Value::Float(f64::INFINITY))]),
        ] {
            let mut geometry = collision_wire(Value::Array(vec![record.clone()]));
            *wire_field(&mut geometry, "fixed") = Value::Bool(true);
            let geometry = crate::graph_api::encode_cbor(&geometry).unwrap();
            let math = crate::graph_api::encode_cbor(
                &serde_json::json!({"temperatures":[],"exponentials":[],"pair-padding":0.35}),
            )
            .unwrap();
            let output = search_bytes(&geometry, &math).unwrap();
            let mut output: Value = ciborium::de::from_reader(output.as_slice()).unwrap();
            assert_eq!(
                crate::graph_api::encode_cbor(wire_field(&mut output, "selected")).unwrap(),
                crate::graph_api::encode_cbor(&Value::Array(vec![record])).unwrap()
            );
            assert_eq!(
                wire_field(&mut output, "missing"),
                &Value::Array(Vec::new())
            );
        }
    }
    #[test]
    fn malformed_packed_candidate_counts_are_rejected() {
        use ciborium::Value;
        let mut packed = packed_candidate_wire();
        *wire_field(&mut packed, "footprints") = Value::Array(Vec::new());
        assert!(CandidateRecords::from_value(packed).is_err());
        let mut packed = packed_candidate_wire();
        wire_field(&mut packed, "batches")
            .as_array_mut()
            .unwrap()
            .pop();
        assert!(CandidateRecords::from_value(packed).is_err());
        let mut packed = packed_candidate_wire();
        let batch = &mut wire_field(&mut packed, "batches").as_array_mut().unwrap()[0];
        wire_field(batch, "positions").as_array_mut().unwrap().pop();
        assert!(CandidateRecords::from_value(packed).is_err());
    }

    #[test]
    fn cached_pair_costs_preserve_ordered_bits_and_eviction() {
        let box_at = |x: f64, y: f64, width: f64, height: f64| Bounds {
            left: x,
            right: x + width,
            bottom: y,
            top: y + height,
        };
        for padding in [0.0, -0.0, 0.2, -0.1, f64::INFINITY, f64::NAN] {
            let spec = CollisionSpec {
                fixed: false,
                edge_lines: Vec::new(),
                obstacles: Vec::new(),
                fixed_arrows: Vec::new(),
                label_padding: 0.0,
                obstacle_padding: 0.0,
                coordinated: true,
                placements: (0..4)
                    .map(|label| CollisionLabel {
                        records: CandidateRecords::Expanded(Vec::new()),
                        edge: None,
                        candidates: (0..5)
                            .map(|choice| {
                                let x = label as f64 * 0.3 + choice as f64 * 0.15;
                                let y = label as f64 * 0.2 - choice as f64 * 0.25;
                                CollisionCandidate {
                                    bounds: if choice == 4 {
                                        box_at(x, y, f64::INFINITY, f64::INFINITY)
                                    } else {
                                        box_at(x, y, 0.7, 0.3)
                                    },
                                    cost: ciborium::Value::Integer(0.into()),
                                    corners: None,
                                    arrow_bounds: vec![
                                        box_at(x - 0.2, y, 0.8, 0.2),
                                        box_at(x + 0.1, y - 0.1, 0.3, 0.8),
                                    ],
                                }
                            })
                            .collect(),
                    })
                    .collect(),
            };
            let mut geometry = spec.geometry(padding).unwrap();
            for capacity in [geometry.pair_cost_cache.len(), 1] {
                geometry.pair_cost_cache.truncate(capacity);
                for i in 0..4 {
                    for choice in 0..5 {
                        for j in (0..4).rev() {
                            for other in (0..5).rev() {
                                let expected =
                                    geometry.compute_pair_cost(i, choice, j, other, padding);
                                for _ in 0..2 {
                                    assert_eq!(
                                        geometry.pair_cost(i, choice, j, other, padding).to_bits(),
                                        expected.to_bits()
                                    );
                                }
                            }
                        }
                    }
                }
            }
            // A different padding must not reuse an entry from the prepared geometry.
            assert_eq!(
                geometry.pair_cost(0, 0, 1, 0, 0.7).to_bits(),
                geometry.compute_pair_cost(0, 0, 1, 0, 0.7).to_bits()
            );
        }
    }

    /// Deterministic coordinates in [-span, span) for the grid invariant.
    fn sampler(seed: u64) -> impl FnMut(f64) -> f64 {
        let mut state = seed;
        move |span| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            span * ((state >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0)
        }
    }

    fn sample_box(next: &mut impl FnMut(f64) -> f64, span: f64) -> Bounds {
        let (x, y) = (next(span), next(span));
        let (width, height) = (next(1.0).abs(), next(1.0).abs());
        Bounds {
            left: x,
            right: x + width,
            bottom: y,
            top: y + height,
        }
    }

    #[test]
    fn line_grid_returns_every_line_passing_the_bounding_test() {
        let mut next = sampler(11);
        let lines: Vec<_> = (0..400)
            .map(|index| StrokeLine {
                start: Point(next(5.0), next(5.0)),
                end: if index % 7 == 0 {
                    Point(f64::NAN, 0.0)
                } else {
                    Point(next(5.0), next(5.0))
                },
                radius: next(0.2).abs(),
            })
            .collect();
        let grid = LineGrid::new(&lines);
        for _ in 0..500 {
            let bounds = sample_box(&mut next, 5.0);
            let near = grid.near(bounds);
            assert!(near.windows(2).all(|pair| pair[0] < pair[1]));
            for (index, line) in lines.iter().enumerate() {
                let extent = line.extent();
                let rejected = extent.right < bounds.left
                    || extent.left > bounds.right
                    || extent.top < bounds.bottom
                    || extent.bottom > bounds.top;
                assert!(
                    rejected || near.binary_search(&index).is_ok(),
                    "line {index}"
                );
            }
        }
    }

    #[test]
    fn painted_strokes_byte_batch_preserves_gaps_widths_and_rejects_unknown_commands() {
        let input = serde_json::json!({
            "accuracy": 0.005,
            "strokes": [
                {"radius": 0.25, "segments": [
                    [[0, 0, 0], false, [["l", [1, 0, 0]]]],
                    [[2, 0, 0], false, [["l", [3, 0, 0]]]]
                ]},
                {"radius": 0.5, "segments": [
                    [[4, 1], true, [["l", [5, 1]]]]
                ]}
            ]
        });
        let encode = |value: &serde_json::Value| crate::graph_api::encode_cbor(value).unwrap();
        let bytes = stroke_lines_bytes(&encode(&input)).unwrap();
        let lines: serde_json::Value = ciborium::de::from_reader(bytes.as_slice()).unwrap();
        assert_eq!(
            lines,
            serde_json::json!([
                {"start":[0.0,0.0],"end":[1.0,0.0],"radius":0.25},
                {"start":[2.0,0.0],"end":[3.0,0.0],"radius":0.25},
                {"start":[4.0,1.0],"end":[5.0,1.0],"radius":0.5},
                {"start":[5.0,1.0],"end":[4.0,1.0],"radius":0.5}
            ])
        );
        let mut invalid = input.clone();
        invalid["strokes"][0]["segments"][0][2][0][0] = "x".into();
        assert!(
            stroke_lines_bytes(&encode(&invalid))
                .unwrap_err()
                .contains("unsupported")
        );
        invalid = input.clone();
        invalid["accuracy"] = 0.into();
        assert!(stroke_lines_bytes(&encode(&invalid)).is_err());
        invalid = input;
        invalid["strokes"][0]["radius"] = (-1).into();
        assert!(stroke_lines_bytes(&encode(&invalid)).is_err());
    }

    #[test]
    fn painted_strokes_flatten_in_drawable_order_with_closures() {
        use ciborium::Value;
        let vertex = |x: f64, y: f64| {
            Value::Array(vec![Value::Float(x), Value::Float(y), Value::Float(0.0)])
        };
        let command = |kind: &str, points: Vec<Value>| {
            std::iter::once(Value::Text(kind.into()))
                .chain(points)
                .collect::<Vec<_>>()
        };
        let spec = StrokeLinesSpec {
            strokes: vec![PaintedStroke {
                radius: 0.25,
                segments: vec![(
                    vertex(0.0, 0.0),
                    true,
                    vec![
                        command("l", vec![vertex(1.0, 0.0), vertex(1.0, 1.0)]),
                        command(
                            "c",
                            vec![vertex(1.0, 2.0), vertex(0.0, 2.0), vertex(0.0, 1.0)],
                        ),
                    ],
                )],
            }],
            accuracy: 0.005,
        };
        let lines = spec.lines().unwrap();
        let cubic = CubicSegment {
            start: Point(1.0, 1.0),
            control_start: Point(1.0, 2.0),
            control_end: Point(0.0, 2.0),
            end: Point(0.0, 1.0),
        };
        let identity = Transform(std::array::from_fn(|row| {
            std::array::from_fn(|column| if row == column { 1.0 } else { 0.0 })
        }));
        let mut expected = vec![[[0.0, 0.0], [1.0, 0.0]], [[1.0, 0.0], [1.0, 1.0]]];
        expected.extend(cubic.lines(&identity, 0.005));
        expected.push([[0.0, 1.0], [0.0, 0.0]]);
        let actual: Vec<_> = lines.iter().map(|line| [line.start, line.end]).collect();
        assert_eq!(actual, expected);
        assert!(lines.iter().all(|line| line.radius == 0.25));
    }

    #[test]
    fn finite_edge_intersection_preserves_cost_types_and_search_addition() {
        let mut candidate = CollisionCandidate {
            bounds: Bounds {
                left: -1.0,
                right: 1.0,
                bottom: -0.5,
                top: 0.5,
            },
            corners: Some([
                Point(-1.0, 0.5),
                Point(1.0, 0.5),
                Point(-1.0, -0.5),
                Point(1.0, -0.5),
            ]),
            cost: ciborium::Value::Integer(2.into()),
            arrow_bounds: Vec::new(),
        };
        assert!(matches!(
            candidate.edge_intersection(&[], &LineGrid::new(&[])),
            Cost::Integer(0)
        ));
        let lines = vec![StrokeLine {
            start: Point(-2.0, 0.0),
            end: Point(2.0, 0.0),
            radius: 0.0,
        }];
        let grid = LineGrid::new(&lines);
        assert_eq!(candidate.edge_intersection(&lines, &grid).float(), 12.0);
        candidate.corners = Some([Point(0.0, 0.0); 4]);
        assert!(matches!(
            candidate.edge_intersection(&lines, &grid),
            Cost::Integer(0)
        ));
        candidate.corners = None;
        let mut spec = CollisionSpec {
            fixed: false,
            edge_lines: lines,
            placements: vec![CollisionLabel {
                records: CandidateRecords::Expanded(Vec::new()),
                candidates: vec![candidate],
                edge: None,
            }],
            obstacles: Vec::new(),
            label_padding: 0.0,
            obstacle_padding: 0.0,
            fixed_arrows: Vec::new(),
            coordinated: false,
        };
        assert_eq!(
            spec.geometry(0.0).unwrap().costs[0][0],
            ciborium::Value::Integer(2.into())
        );
        spec.placements[0].candidates[0].corners = Some([
            Point(-1.0, 0.5),
            Point(1.0, 0.5),
            Point(-1.0, -0.5),
            Point(1.0, -0.5),
        ]);
        assert_eq!(
            spec.geometry(0.0).unwrap().costs[0][0],
            ciborium::Value::Float(14.0)
        );
        spec.edge_lines[0].start.1 = 2.0;
        spec.edge_lines[0].end.1 = 2.0;
        assert_eq!(
            spec.geometry(0.0).unwrap().costs[0][0],
            ciborium::Value::Float(2.0)
        );
    }

    #[test]
    fn path_lines_keep_empty_segments_and_cap_subdivision_depth() {
        let mut spec = PathLinesSpec {
            segments: Vec::new(),
            transform: Transform([
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]),
            accuracy: 0.0,
        };
        assert!(spec.lines().is_empty());
        spec.segments.push(CubicSegment {
            start: Point(2.0, 3.0),
            control_start: Point(2.0, 3.0),
            control_end: Point(2.0, 3.0),
            end: Point(2.0, 3.0),
        });
        spec.segments.push(CubicSegment {
            start: Point(0.0, 0.0),
            control_start: Point(0.0, 1.0),
            control_end: Point(1.0, 1.0),
            end: Point(1.0, 0.0),
        });
        let lines = spec.lines();
        assert_eq!(lines[0], [[[2.0, 3.0], [2.0, 3.0]]]);
        assert_eq!(lines[1].len(), 1 << 10);
        assert_eq!(lines[1][0][0], [0.0, 0.0]);
        assert_eq!(lines[1].last().unwrap()[1], [1.0, 0.0]);
        for pair in lines[1].windows(2) {
            assert_eq!(pair[0][1], pair[1][0]);
        }
        spec.transform = Transform([
            [-2.0, 0.5, 0.0, 1.0],
            [0.0, 3.0, 0.0, -4.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);
        assert_eq!(spec.lines()[0], [[[-1.5, 5.0], [-1.5, 5.0]]]);
    }

    #[test]
    fn exponential_checks_keep_strict_acceptance_and_inclusive_rejection_bounds() {
        let mut check = ExponentialCheck::new(0.0);
        assert!(check.accept(0.75));
        assert!(!check.accept(1.0));
        assert!(check.accept(0.0));
        assert!(!check.accept(2.0));
        assert_eq!((check.lower, check.upper), (Some(0.75), Some(1.0)));
        for argument in [f64::NEG_INFINITY, f64::NAN] {
            let mut check = ExponentialCheck::new(argument);
            assert!(!check.accept(0.0));
            assert!(!check.accept(0.5));
            assert_eq!((check.lower, check.upper), (None, Some(0.0)));
        }
    }

    #[test]
    fn host_exponentials_override_speculation_and_finish_a_verified_search() {
        let candidate = |left, cost| CollisionCandidate {
            corners: None,
            bounds: Bounds {
                left,
                right: left + 1.0,
                bottom: 0.0,
                top: 1.0,
            },
            cost: ciborium::Value::Float(cost),
            arrow_bounds: Vec::new(),
        };
        // Either single move raises the energy, but taking both lowers it.
        let geometry = CollisionSpec {
            fixed: false,
            edge_lines: Vec::new(),
            placements: vec![
                CollisionLabel {
                    records: CandidateRecords::Expanded(Vec::new()),
                    candidates: vec![candidate(0.0, 0.001), candidate(10.9995, 0.0)],
                    edge: None,
                },
                CollisionLabel {
                    records: CandidateRecords::Expanded(Vec::new()),
                    candidates: vec![candidate(10.0, 0.001), candidate(0.9995, 0.0)],
                    edge: None,
                },
            ],
            obstacles: vec![],
            label_padding: 0.0,
            obstacle_padding: 0.0,
            fixed_arrows: Vec::new(),
            coordinated: false,
        }
        .geometry(0.0)
        .unwrap();
        let mut math = SearchMath {
            temperatures: (0..84).map(|sweep| 0.15 * 0.9_f64.powi(sweep)).collect(),
            exponentials: vec![],
            pair_padding: 0.0,
        };
        let predicted = geometry.search(&math).unwrap();
        assert_eq!(predicted.choices, [1, 1]);
        assert!(!predicted.missing.is_empty());
        let mut result = predicted;
        // Deliberately supply host answers that reject every uphill proposal.
        // This forces a different path, exercising correction after a prediction
        // differs, rather than only the usual two identical executions.
        for _ in 0..=160 {
            if result.missing.is_empty() {
                break;
            }
            math.exponentials
                .extend(result.missing.iter().map(|check| (check.argument, 0.0)));
            result = geometry.search(&math).unwrap();
        }
        assert!(result.missing.is_empty());
        assert_eq!(result.choices, [0, 0]);
        assert_eq!(
            result.costs,
            [ciborium::Value::Float(0.001), ciborium::Value::Float(0.001)]
        );
        math.temperatures.pop();
        assert!(geometry.search(&math).is_err());
        assert_eq!(
            Cost::Integer(9_007_199_254_740_993)
                .subtract(Cost::Integer(9_007_199_254_740_992))
                .unwrap()
                .float(),
            1.0
        );
    }

    #[test]
    fn paired_annotations_include_arrows_and_fixed_decorations_in_collisions() {
        let bounds = |left| Bounds {
            left,
            right: left + 1.0,
            bottom: 0.0,
            top: 1.0,
        };
        let spec = CollisionSpec {
            fixed: false,
            edge_lines: Vec::new(),
            placements: [(0.0, 4.0), (4.0, 0.0)]
                .into_iter()
                .map(|(text, arrow)| CollisionLabel {
                    records: CandidateRecords::Expanded(Vec::new()),
                    candidates: vec![CollisionCandidate {
                        corners: None,
                        bounds: bounds(text),
                        cost: ciborium::Value::Integer(0.into()),
                        arrow_bounds: vec![bounds(arrow)],
                    }],
                    edge: None,
                })
                .collect(),
            obstacles: vec![Obstacle {
                bounds: bounds(4.0),
                edge: None,
                self_loop: false,
            }],
            label_padding: 0.25,
            obstacle_padding: 0.0,
            fixed_arrows: vec![bounds(4.0)],
            coordinated: true,
        };
        let geometry = spec.geometry(0.0).unwrap();
        // The two texts are disjoint, but each arrow intersects the other text.
        assert_eq!(geometry.pair_cost(0, 0, 1, 0, 0.0), 4.0);
        assert_eq!(geometry.pair_cost(1, 0, 0, 0, 0.0), 4.0);
        // An arrow outside the text envelope still collides with fixed geometry.
        assert_eq!(geometry.costs[0][0], ciborium::Value::Float(2.0));
        assert_eq!(geometry.costs[1][0], ciborium::Value::Float(3.0));
        assert_eq!(geometry.arrows[0][0], [bounds(4.0)]);
        assert_eq!(geometry.boxes[0][0], bounds(0.0).expanded(0.25));
    }

    #[test]
    fn coordinated_proposals_cross_a_barrier_without_single_uphill_moves() {
        let candidate = |left, cost| CollisionCandidate {
            corners: None,
            bounds: Bounds {
                left,
                right: left + 1.0,
                bottom: 0.0,
                top: 1.0,
            },
            cost: ciborium::Value::Float(cost),
            arrow_bounds: Vec::new(),
        };
        let spec = CollisionSpec {
            fixed: false,
            edge_lines: Vec::new(),
            placements: vec![
                CollisionLabel {
                    records: CandidateRecords::Expanded(Vec::new()),
                    candidates: vec![candidate(0.0, 1.0), candidate(10.0, 0.0)],
                    edge: None,
                },
                CollisionLabel {
                    records: CandidateRecords::Expanded(Vec::new()),
                    candidates: vec![candidate(10.0, 1.0), candidate(0.0, 0.0)],
                    edge: None,
                },
            ],
            obstacles: Vec::new(),
            label_padding: 0.0,
            obstacle_padding: 0.0,
            fixed_arrows: Vec::new(),
            coordinated: true,
        };
        let geometry = spec.geometry(0.0).unwrap();
        let mut math = SearchMath {
            temperatures: (0..84).map(|sweep| 0.15 * 0.9_f64.powi(sweep)).collect(),
            exponentials: Vec::new(),
            pair_padding: 0.0,
        };
        let mut result = geometry.search(&math).unwrap();
        // Certify a host policy rejecting all uphill/tied moves. A simultaneous
        // exchange is the only descending route out of the initial arrangement.
        for _ in 0..1000 {
            if result.missing.is_empty() {
                break;
            }
            math.exponentials
                .extend(result.missing.iter().map(|check| (check.argument, 0.0)));
            result = geometry.search(&math).unwrap();
        }
        assert!(result.missing.is_empty());
        assert_eq!(result.choices, [1, 1]);
        assert_eq!(
            result.costs,
            [ciborium::Value::Float(0.0), ciborium::Value::Float(0.0)]
        );
    }

    #[test]
    fn obstacle_boxes_keep_degenerate_segments_and_handle_empty_paths() {
        let mut spec = ObstacleSpec {
            segments: vec![CubicSegment {
                start: Point(2.0, 3.0),
                control_start: Point(2.0, 3.0),
                control_end: Point(2.0, 3.0),
                end: Point(2.0, 3.0),
            }],
            transform: Transform([
                [2.0, 0.5, 0.0, 1.0],
                [-0.25, 3.0, 0.0, -4.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]),
            pad_x: 0.25,
            pad_y: 0.5,
        };
        assert_eq!(
            spec.boxes(),
            vec![
                Bounds {
                    left: 6.25,
                    right: 6.75,
                    bottom: 4.0,
                    top: 5.0,
                };
                12
            ]
        );
        spec.segments.clear();
        assert!(spec.boxes().is_empty());
    }

    #[test]
    fn collision_geometry_preserves_edge_ownership_and_self_loops() {
        let bounds = Bounds {
            left: 0.0,
            right: 1.0,
            bottom: 0.0,
            top: 1.0,
        };
        let edge = Some(ciborium::Value::Integer(7.into()));
        let mut spec = CollisionSpec {
            fixed: false,
            edge_lines: Vec::new(),
            placements: vec![CollisionLabel {
                records: CandidateRecords::Expanded(Vec::new()),
                candidates: vec![CollisionCandidate {
                    corners: None,
                    bounds,
                    cost: ciborium::Value::Integer(0.into()),
                    arrow_bounds: Vec::new(),
                }],
                edge: edge.clone(),
            }],
            obstacles: vec![Obstacle {
                bounds,
                edge,
                self_loop: false,
            }],
            label_padding: 0.25,
            obstacle_padding: 0.08,
            fixed_arrows: Vec::new(),
            coordinated: false,
        };
        let own_edge = spec.geometry(0.0).unwrap();
        assert_eq!(own_edge.costs[0][0], ciborium::Value::Integer(0.into()));
        assert_eq!(own_edge.boxes[0][0], bounds.expanded(0.25));
        assert_eq!(own_edge.areas[0], [1.0]);
        spec.obstacles[0].edge = Some(ciborium::Value::Float(7.0));
        assert_eq!(spec.geometry(0.0).unwrap().costs, own_edge.costs);
        spec.obstacles[0].self_loop = true;
        let expected = ciborium::Value::Float(1.08 * 1.08);
        assert_eq!(spec.geometry(0.0).unwrap().costs[0][0], expected);
        spec.obstacles[0].self_loop = false;
        spec.obstacles[0].edge = None;
        assert_eq!(spec.geometry(0.0).unwrap().costs[0][0], expected);
        spec.obstacles[0].bounds.left = 10.0;
        spec.obstacles[0].bounds.right = 11.0;
        assert_eq!(
            spec.geometry(0.0).unwrap().costs[0][0],
            ciborium::Value::Float(0.0)
        );
        spec.placements[0].candidates[0].cost = ciborium::Value::Float(-0.0);
        let ciborium::Value::Float(cost) = spec.geometry(0.0).unwrap().costs[0][0] else {
            panic!("an unowned obstacle must convert the cost to a float");
        };
        assert_eq!(cost.to_bits(), 0.0_f64.to_bits());
        spec.placements[0].candidates[0].cost = ciborium::Value::Integer(0.into());
        spec.obstacles.clear();
        assert_eq!(
            spec.geometry(0.0).unwrap().costs[0][0],
            ciborium::Value::Integer(0.into())
        );
    }

    #[test]
    fn collision_cost_additions_retain_the_original_order() {
        let bounds = Bounds {
            left: 0.0,
            right: 1.0,
            bottom: 0.0,
            top: 1.0,
        };
        let cost = 9_007_199_254_740_992.0;
        let spec = CollisionSpec {
            fixed: false,
            edge_lines: Vec::new(),
            placements: vec![CollisionLabel {
                records: CandidateRecords::Expanded(Vec::new()),
                candidates: vec![CollisionCandidate {
                    corners: None,
                    bounds,
                    cost: ciborium::Value::Float(cost),
                    arrow_bounds: Vec::new(),
                }],
                edge: None,
            }],
            obstacles: (0..3)
                .map(|_| Obstacle {
                    bounds,
                    edge: None,
                    self_loop: false,
                })
                .collect(),
            label_padding: 0.0,
            obstacle_padding: 0.0,
            fixed_arrows: Vec::new(),
            coordinated: false,
        };
        // Each individual addition rounds back to the initial value. Summing
        // obstacle penalties before adding the initial cost would change it.
        assert_eq!(
            spec.geometry(0.0).unwrap().costs[0][0],
            ciborium::Value::Float(cost)
        );
        assert_ne!(cost + 3.0, cost);
    }

    #[test]
    fn batched_finite_candidates_match_separate_geometry_calls() {
        let input = serde_json::json!({
            "frames": [
                {"point": [0.0, 0.0], "tangent": [1.0, 1.0]},
                {"point": [1.0, 0.0], "tangent": [0.0, -0.0]},
                {"point": [2.0, 0.0], "tangent": [-2.0, 1.0]}
            ],
            "positions": [0.5, 1.0, 2.0],
            "transform": [[2.0, -0.25, 0.0, 10.0], [0.5, -1.5, 0.0, -4.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
            "origin": [10.0, -4.0, 0.0],
            "corners": [[-2.0, 0.3, 0.0], [2.0, 0.3, 0.0], [-2.0, -0.3, 0.0], [2.0, -0.3, 0.0]],
            "total": 3.0, "preferred": 1.5, "accuracy": 0.001, "side": 1.0,
            "preferred-side": -1.0, "path-index": 2, "path-shift": 0.25,
            "gap": 0.2, "clear-box": true, "fixed": false,
            "attachment": {
                "carrier": [{"start": [-2, -2], "control-start": [-1, -1], "control-end": [1, 1], "end": [2, 2]}],
                "paths": [
                    [{"start": [-0.2, -0.2], "control-start": [-0.1, -0.1], "control-end": [0.1, 0.1], "end": [0.2, 0.2]}],
                    [{"start": [0, 0], "control-start": [2, -1], "control-end": [1, 3], "end": [2, 1]}],
                    []
                ]
            }
        });
        let encode = |value: &serde_json::Value| crate::graph_api::encode_cbor(value).unwrap();
        let flatten = |segments: &serde_json::Value| {
            let bytes = path_lines_bytes(&encode(&serde_json::json!({
                "segments": segments, "transform": input["transform"], "accuracy": input["accuracy"],
            }))).unwrap();
            let parts: Vec<Vec<[Point; 2]>> = ciborium::de::from_reader(bytes.as_slice()).unwrap();
            parts.into_iter().flatten().collect::<Vec<_>>()
        };
        let carrier = flatten(&input["attachment"]["carrier"]);
        let spec: CandidateSpec = serde_json::from_value(input.clone()).unwrap();
        let mut base: CandidateSpec = serde_json::from_value(input.clone()).unwrap();
        base.attachment = None;
        let mut expected = base.candidates().unwrap();
        let mut moved = false;
        for (index, candidate) in expected.iter_mut().enumerate() {
            let mut lines = flatten(&input["attachment"]["paths"][index]);
            lines.extend_from_slice(&carrier);
            let [x, y, _] = spec.transform.point(spec.frames[index].point);
            let request = serde_json::json!({
                "initial": spec.gap - candidate.nearest,
                "lines": lines.iter().map(|line| line.map(|Point(x,y)| [x,y])).collect::<Vec<_>>(),
                "corners": spec.corners.map(|[x,y,_]| [x,y]), "frame": [x,y],
                "outward": [candidate.outward[0],candidate.outward[1]], "gap": spec.gap,
            });
            let bytes = attachment_offset_bytes(&encode(&request)).unwrap();
            let offset: f64 = ciborium::de::from_reader(bytes.as_slice()).unwrap();
            let position =
                std::array::from_fn(|i| spec.frames[index].point[i] + candidate.normal[i] * offset);
            moved |= position != candidate.position;
            candidate.position = position;
            let center = spec.transform.point(position);
            candidate.corners = spec
                .corners
                .map(|corner| std::array::from_fn(|i| center[i] + corner[i]));
            candidate.bounds = Bounds {
                left: candidate
                    .corners
                    .iter()
                    .map(|p| p[0])
                    .fold(f64::INFINITY, f64::min),
                right: candidate
                    .corners
                    .iter()
                    .map(|p| p[0])
                    .fold(f64::NEG_INFINITY, f64::max),
                bottom: candidate
                    .corners
                    .iter()
                    .map(|p| p[1])
                    .fold(f64::INFINITY, f64::min),
                top: candidate
                    .corners
                    .iter()
                    .map(|p| p[1])
                    .fold(f64::NEG_INFINITY, f64::max),
            };
        }
        assert!(moved, "the fixture must exercise a finite-path correction");
        assert_eq!(
            candidates_bytes(&encode(&input)).unwrap(),
            crate::graph_api::encode_cbor(&expected).unwrap()
        );
        for (clear_box, fixed) in [(false, false), (true, true)] {
            let mut with_attachment: CandidateSpec = serde_json::from_value(input.clone()).unwrap();
            with_attachment.clear_box = clear_box;
            with_attachment.fixed = fixed;
            let attached =
                crate::graph_api::encode_cbor(&with_attachment.candidates().unwrap()).unwrap();
            with_attachment.attachment = None;
            assert_eq!(
                attached,
                crate::graph_api::encode_cbor(&with_attachment.candidates().unwrap()).unwrap()
            );
        }
        let mut invalid: CandidateSpec = serde_json::from_value(input).unwrap();
        invalid.attachment.as_mut().unwrap().paths.pop();
        assert!(
            matches!(invalid.candidates(), Err(error) if error.contains("attachment paths and frames"))
        );
    }

    #[test]
    fn transformed_label_clears_its_box_and_preserves_fixed_position() {
        let mut spec = CandidateSpec {
            attachment: None,
            endpoint_direction: None,
            frames: vec![Frame {
                point: [1.0, 2.0],
                tangent: [3.0, 0.0],
            }],
            positions: vec![5.0],
            transform: Transform([
                [2.0, 0.0, 0.0, 10.0],
                [0.0, 3.0, 0.0, -4.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]),
            origin: [10.0, -4.0, 0.0],
            corners: [
                [-2.0, -3.0, 0.0],
                [2.0, -3.0, 0.0],
                [2.0, 3.0, 0.0],
                [-2.0, 3.0, 0.0],
            ],
            total: 10.0,
            preferred: 5.0,
            accuracy: 1e-6,
            side: 1.0,
            preferred_side: 1.0,
            path_index: 2,
            path_shift: 0.5,
            gap: 0.25,
            clear_box: true,
            fixed: false,
        };
        let candidate = spec.candidates().unwrap().remove(0);
        assert_eq!(candidate.position, [1.0, 3.25]);
        assert_eq!(candidate.normal, [0.0, 1.0]);
        assert_eq!(candidate.outward, [0.0, 3.0, 0.0]);
        assert_eq!(candidate.nearest, -1.0);
        assert_eq!(
            candidate.corners,
            [
                [10.0, 2.75, 0.0],
                [14.0, 2.75, 0.0],
                [14.0, 8.75, 0.0],
                [10.0, 8.75, 0.0],
            ]
        );
        assert_eq!(
            (candidate.bounds.left, candidate.bounds.right),
            (10.0, 14.0)
        );
        assert_eq!(
            (candidate.bounds.bottom, candidate.bounds.top),
            (2.75, 8.75)
        );
        assert_eq!(candidate.path_shift, 0.5);
        assert_eq!(candidate.cost, 0.0);
        spec.fixed = true;
        assert_eq!(spec.candidates().unwrap()[0].position, [1.0, 2.0]);
        // A dangling endpoint is a frame for an outward label, even when pinned.
        spec.endpoint_direction = Some(Point(3.0, 0.0));
        spec.clear_box = false;
        let endpoint = spec.candidates().unwrap().remove(0);
        assert_eq!(endpoint.normal, [1.0, 0.0]);
        assert_eq!(endpoint.position, [2.25, 2.0]);
        assert_eq!(endpoint.nearest, -1.0);
        spec.positions.clear();
        assert!(matches!(spec.candidates(), Err(error) if error.contains("equal lengths")));
    }
}
