//! Sampled unions of round tubes and vertex discs.

use std::collections::{BTreeMap, BTreeSet};

use kurbo::{flatten, BezPath, Line, ParamCurveNearest, PathEl, Point, Rect, Shape, Vec2};

type Sample = (i64, i64);

pub(super) const RESOLUTION_ERROR: &str =
    "hull is too thin for its drawing extent; increase hull_radius or select a smaller region";

// Every marching-square crossing is indexed by the doubled integer position
// of its grid-edge midpoint. Adjacent cells therefore join exactly.
pub(super) fn union(
    paths: &[&BezPath],
    vertices: &[Point],
    radius: f64,
) -> Result<BezPath, String> {
    const MAX_CELLS_PER_SIDE: f64 = 768.0;
    let mut bounds = vertices
        .iter()
        .copied()
        .map(|p| Rect::from_points(p, p))
        .chain(paths.iter().map(|path| path.bounding_box()));
    let Some(first) = bounds.next() else {
        return Ok(BezPath::new());
    };
    let bounds = bounds
        .fold(first, |bounds, rect| bounds.union(rect))
        .inflate(radius, radius);
    // Resolve at least three samples per radius. Keep the grid bounded by
    // the selected region's own extent, including the tube and disc radii.
    let min_step = bounds.width().max(bounds.height()) / MAX_CELLS_PER_SIDE;
    let max_step = radius / 3.0;
    if min_step > max_step {
        return Err(RESOLUTION_ERROR.into());
    }
    let step = 0.018_f64.min(max_step).max(min_step);
    let max_coordinate = [bounds.x0, bounds.x1, bounds.y0, bounds.y1]
        .into_iter()
        .map(f64::abs)
        .fold(0.0, f64::max);
    if !bounds.is_finite() || !step.is_finite() || max_coordinate / step >= (1_i64 << 49) as f64 {
        return Err(RESOLUTION_ERROR.into());
    }
    let left = (bounds.x0 / step).floor() as i64 - 1;
    let top = (bounds.y0 / step).floor() as i64 - 1;
    let right = (bounds.x1 / step).ceil() as i64 + 1;
    let bottom = (bounds.y1 / step).ceil() as i64 + 1;
    let stride = (right - left + 1) as usize;
    let index = |x: i64, y: i64| (y - top) as usize * stride + (x - left) as usize;
    let mut inside = vec![false; stride * (bottom - top + 1) as usize];
    let mut add_capsule = |a: Point, b: Point| {
        let line = Line::new(a, b);
        let rect = line.bounding_box().inflate(radius, radius);
        for y in (rect.y0 / step).ceil() as i64..=(rect.y1 / step).floor() as i64 {
            for x in (rect.x0 / step).ceil() as i64..=(rect.x1 / step).floor() as i64 {
                let cell = &mut inside[index(x, y)];
                if *cell {
                    continue;
                }
                let p = Point::new(x as f64 * step, y as f64 * step);
                let distance = if a == b {
                    p.distance_squared(a)
                } else {
                    line.nearest(p, step).distance_sq
                };
                *cell = distance <= radius * radius;
            }
        }
    };
    for path in paths {
        let mut previous = Point::ORIGIN;
        let mut start = Point::ORIGIN;
        flatten(path.iter(), step / 2.0, |element| match element {
            PathEl::MoveTo(p) => {
                start = p;
                previous = p;
            }
            PathEl::ClosePath => add_capsule(previous, start),
            PathEl::LineTo(p) => {
                // Split long diagonals to bound each capsule's scan area.
                let count = (previous.distance(p) / (radius * 2.0).max(step))
                    .ceil()
                    .max(1.0) as usize;
                for i in 0..count {
                    add_capsule(
                        previous.lerp(p, i as f64 / count as f64),
                        previous.lerp(p, (i + 1) as f64 / count as f64),
                    );
                }
                previous = p;
            }
            _ => {}
        });
    }
    for &vertex in vertices {
        add_capsule(vertex, vertex);
    }
    let mut neighbors: BTreeMap<Sample, BTreeSet<Sample>> = BTreeMap::new();
    let cells = (top..bottom).flat_map(|y| (left..right).map(move |x| (x, y)));
    for (x, y) in cells {
        let corners =
            [(x, y), (x + 1, y), (x + 1, y + 1), (x, y + 1)].map(|(x, y)| inside[index(x, y)]);
        if corners.iter().all(|&inside| inside == corners[0]) {
            continue;
        }
        let mid = [
            (2 * x + 1, 2 * y),
            (2 * x + 2, 2 * y + 1),
            (2 * x + 1, 2 * y + 2),
            (2 * x, 2 * y + 1),
        ];
        let cuts: Vec<_> = (0..4)
            .filter(|&i| corners[i] != corners[(i + 1) % 4])
            .collect();
        let mut join = |a: usize, b: usize| {
            neighbors.entry(mid[a]).or_default().insert(mid[b]);
            neighbors.entry(mid[b]).or_default().insert(mid[a]);
        };
        match *cuts.as_slice() {
            [a, b] => join(a, b),
            // Separate diagonal islands rather than creating a self-intersection.
            [0, 1, 2, 3] if corners[0] => {
                join(3, 0);
                join(1, 2);
            }
            [0, 1, 2, 3] => {
                join(0, 1);
                join(2, 3);
            }
            _ => {}
        }
    }
    let mut result = BezPath::new();
    while let Some((&start, _)) = neighbors.first_key_value() {
        let mut samples = Vec::new();
        let mut at = start;
        while let Some(next) = neighbors.get(&at).and_then(|next| next.first()).copied() {
            samples.push(Point::new(
                at.0 as f64 * step / 2.0,
                at.1 as f64 * step / 2.0,
            ));
            for (a, b) in [(at, next), (next, at)] {
                if let Some(entry) = neighbors.get_mut(&a) {
                    entry.remove(&b);
                    if entry.is_empty() {
                        neighbors.remove(&a);
                    }
                }
            }
            at = next;
            if at == start {
                break;
            }
        }
        if at == start && samples.len() >= 3 {
            // Soften the grid staircase. The SVG writer can omit collinear
            // vertices once it knows the rounded canvas coordinates.
            let count = samples.len();
            for (index, point) in samples.iter().enumerate() {
                let softened = point.midpoint(samples[(index + 1) % count]);
                if index == 0 {
                    result.move_to(softened);
                } else {
                    result.line_to(softened);
                }
            }
            result.close_path();
        }
    }
    if result.elements().is_empty() {
        return Err(RESOLUTION_ERROR.into());
    }
    Ok(result)
}

/// Whether `point` lies in the painted (even-odd) region of `hull`, as the
/// SVG draws it: inside the outer contour and outside every hole.
#[cfg(test)]
fn hull_contains(hull: &BezPath, point: Point) -> bool {
    hull.winding(point).rem_euclid(2) == 1
}

/// The lines and vertices that cut the drawing into faces: the region's own
/// and the rest of the graph's.
pub(super) struct FaceWalls<'a> {
    pub own: &'a [&'a BezPath],
    pub foreign: &'a [&'a BezPath],
    pub own_vertices: &'a [Point],
    pub foreign_vertices: &'a [Point],
}

/// Wall cells of a face raster: free, a line or vertex of the region itself,
/// or one of the rest of the drawing.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Cell {
    Free,
    Own,
    Foreign,
}

/// Where a counterterm glyph sits in a face of the drawing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum FacePlacement {
    /// In the largest face whose ring fits clear of lines and labels.
    Fits(Point),
    /// At the centre of the largest face whose ring there covers no label,
    /// although it may cross lines.
    Centre(Point),
}

impl FacePlacement {
    pub(super) fn point(self) -> Point {
        match self {
            Self::Fits(point) | Self::Centre(point) => point,
        }
    }
}

/// Where to put a counterterm glyph in a face of the drawing enclosed only by
/// the region's own lines and vertices (`own`, `own_vertices`).
///
/// First the largest such face in which the glyph's ring fits, at its point
/// farthest from every wall and label: lines are drawn `line_half_width` wide
/// on each side and the ring clears that ink too, and label boxes in `labels`
/// neither split nor reject a face but take clearance from the ring.
/// Otherwise the faces from largest to smallest, each at its centre (the
/// point farthest from its lines alone), taking the first whose ring there
/// covers no label. `None` when neither holds, which leaves the band rule.
///
/// Faces are the bounded regions the drawn lines cut the plane into; a face
/// touching any line or vertex of `foreign` or `foreign_vertices` (another
/// line of the graph, or a vertex outside the region) does not count. Lines of
/// nested spinneys belong to the region, so they count as its own.
pub(super) fn face_point(
    walls: &FaceWalls<'_>,
    labels: &[Rect],
    vertex_radius: f64,
    line_half_width: f64,
    glyph_radius: f64,
) -> Option<FacePlacement> {
    let search = |step_floor: f64, fit: bool| {
        face_search(
            walls,
            labels,
            [vertex_radius, line_half_width, glyph_radius, step_floor],
            fit,
        )
    };
    // Cells of an eighth of the glyph radius resolve faces large enough for
    // the ring; loops smaller than the glyph, as in miniatures, need the
    // drawing's own resolution to find their centres.
    search(glyph_radius / 8.0, true).or_else(|| search(0.0, false))
}

/// [`face_point`] on one raster whose cells are at least `step_floor` wide;
/// `fit` also tries faces that hold the ring before taking face centres.
fn face_search(
    walls: &FaceWalls<'_>,
    labels: &[Rect],
    [vertex_radius, line_half_width, glyph_radius, step_floor]: [f64; 4],
    fit: bool,
) -> Option<FacePlacement> {
    /// Cells along the larger side of the raster, at most.
    const CELLS: f64 = 240.0;
    let FaceWalls {
        own,
        foreign,
        own_vertices,
        foreign_vertices,
    } = *walls;
    let bounds = own
        .iter()
        .chain(foreign)
        .map(|path| path.bounding_box())
        .chain(
            own_vertices
                .iter()
                .chain(foreign_vertices)
                .map(|point| Rect::from_points(*point, *point)),
        )
        .reduce(|bounds, rect| bounds.union(rect))?;
    // One free cell of margin keeps the outer face connected around the drawing.
    let step = (bounds.width().max(bounds.height()) / CELLS).max(step_floor);
    if !(step.is_finite() && step > 0.0) {
        return None;
    }
    let bounds = bounds.inflate(2.0 * step, 2.0 * step);
    let columns = (bounds.width() / step).ceil() as usize + 1;
    let rows = (bounds.height() / step).ceil() as usize + 1;
    let at =
        |i: usize, j: usize| Point::new(bounds.x0 + step * i as f64, bounds.y0 + step * j as f64);
    let mut cells = vec![Cell::Free; columns * rows];
    // A line blocks every cell it passes through, so a 4-connected path of free
    // cells cannot cross it.
    let reach = step * std::f64::consts::FRAC_1_SQRT_2 * 1.05;
    let mut paint = |segment: Line, radius: f64, cell: Cell| {
        let area = segment.bounding_box().inflate(radius, radius);
        let (i0, j0) = (
            ((area.x0 - bounds.x0) / step).floor().max(0.0) as usize,
            ((area.y0 - bounds.y0) / step).floor().max(0.0) as usize,
        );
        let (i1, j1) = (
            (((area.x1 - bounds.x0) / step).ceil() as usize).min(columns - 1),
            (((area.y1 - bounds.y0) / step).ceil() as usize).min(rows - 1),
        );
        for j in j0..=j1 {
            for i in i0..=i1 {
                let index = j * columns + i;
                if cells[index] != Cell::Foreign
                    && segment.nearest(at(i, j), 1e-9).distance_sq <= radius * radius
                {
                    cells[index] = cell;
                }
            }
        }
    };
    let segments = |path: &BezPath| {
        let mut segments = Vec::new();
        let (mut start, mut last) = (Point::ORIGIN, Point::ORIGIN);
        flatten(path, step / 4.0, |element| match element {
            PathEl::MoveTo(point) => (start, last) = (point, point),
            PathEl::LineTo(point) => {
                segments.push(Line::new(last, point));
                last = point;
            }
            PathEl::ClosePath => {
                segments.push(Line::new(last, start));
                last = start;
            }
            _ => {}
        });
        segments
    };
    for path in own {
        for segment in segments(path) {
            paint(segment, reach, Cell::Own);
        }
    }
    for path in foreign {
        for segment in segments(path) {
            paint(segment, reach, Cell::Foreign);
        }
    }
    let dot = |point: Point| Line::new(point, point);
    for point in own_vertices {
        paint(dot(*point), vertex_radius.max(reach), Cell::Own);
    }
    for point in foreign_vertices {
        paint(dot(*point), vertex_radius.max(reach), Cell::Foreign);
    }
    // Labels are free space for faces but no room for the ring.
    let labelled: Vec<bool> = (0..cells.len())
        .map(|index| {
            let point = at(index % columns, index / columns);
            labels
                .iter()
                .any(|label| label.inflate(reach, reach).contains(point))
        })
        .collect();
    let neighbours = |index: usize| {
        let (i, j) = (index % columns, index / columns);
        [
            (i > 0).then(|| index - 1),
            (i + 1 < columns).then(|| index + 1),
            (j > 0).then(|| index - columns),
            (j + 1 < rows).then(|| index + columns),
        ]
        .into_iter()
        .flatten()
    };
    let around = |index: usize| {
        let (i, j) = ((index % columns) as isize, (index / columns) as isize);
        (-1..=1_isize)
            .flat_map(move |dj| (-1..=1_isize).map(move |di| (i + di, j + dj)))
            .filter(move |&(x, y)| {
                (x, y) != (i, j)
                    && x >= 0
                    && y >= 0
                    && (x as usize) < columns
                    && (y as usize) < rows
            })
            .map(|(x, y)| y as usize * columns + x as usize)
    };
    // Distance of every cell to the nearest seed cell, in cells: chessboard
    // steps, which never exceed the Euclidean distance.
    let distances = |seed: &dyn Fn(usize) -> bool| {
        let mut distance = vec![usize::MAX; cells.len()];
        let mut queue = std::collections::VecDeque::new();
        for (index, value) in distance.iter_mut().enumerate() {
            if seed(index) {
                *value = 0;
                queue.push_back(index);
            }
        }
        while let Some(index) = queue.pop_front() {
            for next in around(index) {
                if distance[next] == usize::MAX {
                    distance[next] = distance[index] + 1;
                    queue.push_back(next);
                }
            }
        }
        distance
    };
    // Clearance from lines and labels decides whether the ring fits; the
    // centre of a face is its point farthest from its lines alone.
    let to_walls_and_labels = distances(&|index| cells[index] != Cell::Free || labelled[index]);
    let to_walls = distances(&|index| cells[index] != Cell::Free);
    // Faces: 4-connected free components that stay off the border and touch
    // only the region's own walls.
    let mut face = vec![usize::MAX; cells.len()];
    let mut faces: Vec<Vec<usize>> = Vec::new();
    for seed in 0..cells.len() {
        if cells[seed] != Cell::Free || face[seed] != usize::MAX {
            continue;
        }
        let mut members = vec![seed];
        face[seed] = seed;
        let (mut bounded, mut enclosed) = (true, true);
        let mut cursor = 0;
        while cursor < members.len() {
            let index = members[cursor];
            cursor += 1;
            let (i, j) = (index % columns, index / columns);
            if i == 0 || j == 0 || i + 1 == columns || j + 1 == rows {
                bounded = false;
            }
            for next in neighbours(index) {
                match cells[next] {
                    Cell::Free if face[next] == usize::MAX => {
                        face[next] = seed;
                        members.push(next);
                    }
                    Cell::Foreign => enclosed = false,
                    _ => {}
                }
            }
        }
        if bounded && enclosed {
            faces.push(members);
        }
    }
    // Largest first; a stable sort keeps raster order between equal sizes.
    faces.sort_by_key(|members| std::cmp::Reverse(members.len()));
    // The clearest cell of a face by `distance`; on the plateau chessboard
    // distances leave, the one nearest the face's centroid.
    let clearest = |members: &[usize], distance: &[usize]| {
        let centroid = members.iter().fold(Vec2::ZERO, |sum, index| {
            sum + at(index % columns, index / columns).to_vec2()
        }) / members.len() as f64;
        let off_centre =
            |index: usize| at(index % columns, index / columns).distance(centroid.to_point());
        members
            .iter()
            .copied()
            .max_by(|a, b| {
                distance[*a]
                    .cmp(&distance[*b])
                    .then_with(|| off_centre(*b).total_cmp(&off_centre(*a)))
            })
            .expect("a face has a cell")
    };
    // First the largest face whose ring fits clear of lines and labels; a wall
    // cell reaches up to `reach` beyond its centre.
    for members in faces.iter().filter(|_| fit) {
        let centre = clearest(members, &to_walls_and_labels);
        let clearance = to_walls_and_labels[centre] as f64 * step - reach - line_half_width;
        if clearance >= glyph_radius {
            return Some(FacePlacement::Fits(at(centre % columns, centre / columns)));
        }
    }
    // Otherwise the centre of the largest face whose ring there covers no
    // label, even if the ring crosses lines.
    faces.iter().find_map(|members| {
        let index = clearest(members, &to_walls);
        let centre = at(index % columns, index / columns);
        labels
            .iter()
            .all(|label| {
                let nearest = Point::new(
                    centre.x.clamp(label.x0, label.x1),
                    centre.y.clamp(label.y0, label.y1),
                );
                nearest.distance(centre) >= glyph_radius
            })
            .then_some(FacePlacement::Centre(centre))
    })
}

/// Where to draw a counterterm glyph of radius `glyph_radius` for a contracted
/// region with hull `hull`. The centre lies in the hull's painted region, never
/// in one of its holes; among such sampled points it prefers those whose glyph
/// covers no label box in `labels`, then keeps the largest share of the
/// glyph's ring inside that region, then stays farthest from every drawn line
/// in `lines` and every vertex, then nearest `fallback`. The glyph may
/// cover its own region's lines. Returns `fallback` when no sample lies in the
/// painted region.
pub(super) fn glyph_point(
    hull: &BezPath,
    lines: &[&BezPath],
    vertices: &[Point],
    labels: &[Rect],
    vertex_radius: f64,
    glyph_radius: f64,
    fallback: Point,
) -> Point {
    const RING: usize = 24;
    /// Candidates whose clearance is compared, at most.
    const CLEARANCE_CANDIDATES: usize = 64;
    let bounds = hull.bounding_box();
    if !(bounds.width() > 0.0 && bounds.height() > 0.0) {
        return fallback;
    }
    let polyline = |path: &BezPath| {
        let mut segments = Vec::new();
        let (mut start, mut last) = (Point::ORIGIN, Point::ORIGIN);
        flatten(path, 1e-3, |element| match element {
            PathEl::MoveTo(point) => (start, last) = (point, point),
            PathEl::LineTo(point) => {
                segments.push(Line::new(last, point));
                last = point;
            }
            PathEl::ClosePath => {
                segments.push(Line::new(last, start));
                last = start;
            }
            _ => {}
        });
        segments
    };
    let edges = polyline(hull);
    let obstacles: Vec<Line> = lines.iter().flat_map(|path| polyline(path)).collect();
    let clearance = |point: Point| {
        let line = obstacles
            .iter()
            .map(|segment| segment.nearest(point, 1e-6).distance_sq.sqrt())
            .fold(f64::INFINITY, f64::min);
        let vertex = vertices
            .iter()
            .map(|vertex| vertex.distance(point) - vertex_radius)
            .fold(f64::INFINITY, f64::min);
        line.min(vertex)
    };
    // Rasterize the painted region once; a finer raster only when the coarse
    // one misses a thin region.
    for cells in [120_usize, 480] {
        let step = bounds.width().max(bounds.height()) / cells as f64;
        let (columns, rows) = (
            (bounds.width() / step).ceil() as usize + 1,
            (bounds.height() / step).ceil() as usize + 1,
        );
        let at = |i: usize, j: usize| {
            Point::new(bounds.x0 + step * i as f64, bounds.y0 + step * j as f64)
        };
        // Even-odd scanlines: a cell is painted when an odd number of hull
        // edges cross its row to the left of it.
        let mut painted = vec![false; rows * columns];
        let mut crossings = Vec::new();
        for j in 0..rows {
            let y = bounds.y0 + step * j as f64;
            crossings.clear();
            for edge in &edges {
                let (a, b) = (edge.p0, edge.p1);
                if (a.y > y) != (b.y > y) {
                    crossings.push(a.x + (y - a.y) * (b.x - a.x) / (b.y - a.y));
                }
            }
            crossings.sort_by(f64::total_cmp);
            for [enter, leave] in crossings.as_chunks::<2>().0 {
                let first = ((enter - bounds.x0) / step).ceil().max(0.0);
                let last = ((leave - bounds.x0) / step)
                    .floor()
                    .min((columns - 1) as f64);
                if first <= last {
                    for i in first as usize..=last as usize {
                        painted[j * columns + i] = true;
                    }
                }
            }
        }
        let painted_at = |point: Point| {
            let (u, v) = (
                ((point.x - bounds.x0) / step).round(),
                ((point.y - bounds.y0) / step).round(),
            );
            u >= 0.0
                && v >= 0.0
                && (u as usize) < columns
                && (v as usize) < rows
                && painted[v as usize * columns + u as usize]
        };
        let coverage = |center: Point| {
            (0..RING)
                .filter(|index| {
                    let angle = std::f64::consts::TAU * *index as f64 / RING as f64;
                    painted_at(center + glyph_radius * Vec2::new(angle.cos(), angle.sin()))
                })
                .count()
        };
        // Whether the glyph's disc would cover part of a label.
        let clear_of_labels = |center: Point| {
            labels.iter().all(|label| {
                let nearest = Point::new(
                    center.x.clamp(label.x0, label.x1),
                    center.y.clamp(label.y0, label.y1),
                );
                nearest.distance(center) >= glyph_radius
            })
        };
        // First clear of labels where any candidate is, then the share of the
        // ring inside the painted region.
        let mut best = (false, 0);
        let mut candidates = Vec::new();
        for j in 0..rows {
            for i in 0..columns {
                if !painted[j * columns + i] {
                    continue;
                }
                let point = at(i, j);
                let key = (clear_of_labels(point), coverage(point));
                if key > best {
                    best = key;
                    candidates.clear();
                }
                if key == best {
                    candidates.push(point);
                }
            }
        }
        if candidates.is_empty() {
            continue;
        }
        // Among the best covered, prefer clear of lines, then near the fallback.
        let stride = candidates.len().div_ceil(CLEARANCE_CANDIDATES);
        return candidates
            .into_iter()
            .step_by(stride)
            .map(|point| (clearance(point) - 0.05 * point.distance(fallback), point))
            .max_by(|(a, _), (b, _)| a.total_cmp(b))
            .map_or(fallback, |(_, point)| point);
    }
    fallback
}

#[cfg(test)]
mod tests {
    use kurbo::{Circle, ParamCurveNearest, PathEl, PathSeg, Vec2};

    /// Two arcs from (-2, 0) to (2, 0) bulging by `bulge`: a bubble's lines.
    fn bubble(bulge: f64) -> [BezPath; 2] {
        [bulge, -bulge].map(|height| {
            let mut arc = BezPath::new();
            arc.move_to((-2.0, 0.0));
            arc.quad_to((0.0, 2.0 * height), (2.0, 0.0));
            arc
        })
    }

    #[test]
    fn counterterm_glyphs_centre_in_the_face_their_own_lines_enclose() {
        let [upper, lower] = bubble(1.0);
        let ends = [Point::new(-2.0, 0.0), Point::new(2.0, 0.0)];
        let legs = [
            Line::new((-2.0, 0.0), (-4.0, 1.0)).into(),
            Line::new((2.0, 0.0), (4.0, 1.0)).into(),
        ];
        let legs: Vec<BezPath> = legs.iter().map(|leg: &PathSeg| leg.to_path(1e-3)).collect();
        let point = super::face_point(
            &FaceWalls {
                own: &[&upper, &lower],
                foreign: &legs.iter().collect::<Vec<_>>(),
                own_vertices: &ends,
                foreign_vertices: &[],
            },
            &[],
            0.1,
            0.0,
            0.4,
        )
        .expect("the loop face holds the ring");
        assert!(matches!(point, FacePlacement::Fits(_)), "{point:?}");
        assert!(point.point().distance(Point::ORIGIN) < 0.1, "{point:?}");

        // A foreign line inside the face, or a foreign vertex in it, rejects it.
        let chord: BezPath = PathSeg::from(Line::new((0.0, -0.3), (0.0, 0.3))).to_path(1e-3);
        assert_eq!(
            super::face_point(
                &FaceWalls {
                    own: &[&upper, &lower],
                    foreign: &[&chord],
                    own_vertices: &ends,
                    foreign_vertices: &[]
                },
                &[],
                0.1,
                0.0,
                0.2
            ),
            None
        );
        assert_eq!(
            super::face_point(
                &FaceWalls {
                    own: &[&upper, &lower],
                    foreign: &[],
                    own_vertices: &ends,
                    foreign_vertices: &[Point::new(0.5, 0.0)]
                },
                &[],
                0.1,
                0.0,
                0.2
            ),
            None
        );
        // A label in the face takes room from the ring without rejecting the
        // face: the cross moves off it, still inside the loop.
        let label = Rect::new(-0.3, -0.15, 0.3, 0.15);
        let point = super::face_point(
            &FaceWalls {
                own: &[&upper, &lower],
                foreign: &[],
                own_vertices: &ends,
                foreign_vertices: &[],
            },
            &[label],
            0.1,
            0.0,
            0.4,
        )
        .expect("the loop face still holds the ring beside the label");
        assert!(matches!(point, FacePlacement::Fits(_)), "{point:?}");
        let point = point.point();
        let nearest = Point::new(
            point.x.clamp(label.x0, label.x1),
            point.y.clamp(label.y0, label.y1),
        );
        assert!(nearest.distance(point) >= 0.4, "{point:?}");
        assert!(point.y.abs() < 0.5 && point.x.abs() < 2.0, "{point:?}");
        // A loop too small for the ring, as in a miniature, still takes the
        // cross at its centre, the ring crossing the loop's lines.
        let [upper, lower] = bubble(0.15);
        let thin = FaceWalls {
            own: &[&upper, &lower],
            foreign: &[],
            own_vertices: &ends,
            foreign_vertices: &[],
        };
        let Some(FacePlacement::Centre(point)) = super::face_point(&thin, &[], 0.1, 0.0, 0.4)
        else {
            panic!("a thin loop centres the cross");
        };
        assert!(point.distance(Point::ORIGIN) < 0.1, "{point:?}");
        // A label at that centre leaves the band rule to place the cross.
        let centred = Rect::new(-0.2, -0.05, 0.2, 0.05);
        assert_eq!(super::face_point(&thin, &[centred], 0.1, 0.0, 0.4), None);
    }

    #[test]
    fn a_tadpole_too_small_for_the_ring_centres_the_cross_in_its_loop() {
        // A tadpole: one loop through its vertex, with a leg leaving it.
        let vertex = Point::new(0.0, -0.3);
        let loop_line = Circle::new((0.0, 0.0), 0.3).to_path(1e-4);
        let leg: BezPath = PathSeg::from(Line::new((0.0, -0.3), (0.0, -2.0))).to_path(1e-3);
        let walls = FaceWalls {
            own: &[&loop_line],
            foreign: &[&leg],
            own_vertices: &[vertex],
            foreign_vertices: &[],
        };
        let Some(FacePlacement::Centre(point)) = super::face_point(&walls, &[], 0.05, 0.0, 0.5)
        else {
            panic!("the tadpole's loop holds the cross centre");
        };
        assert!(point.distance(Point::ORIGIN) < 0.1, "{point:?}");
    }

    #[test]
    fn counterterm_glyphs_sit_in_the_painted_band_of_a_ring_hull() {
        // A ring hull: the centroid of its vertices falls in the hole.
        let mut ring = BezPath::new();
        ring.extend(Circle::new((0.0, 0.0), 3.0).path_elements(1e-3));
        ring.extend(Circle::new((0.0, 0.0), 2.0).path_elements(1e-3));
        let centroid = Point::ORIGIN;
        assert!(!super::hull_contains(&ring, centroid));
        let point = super::glyph_point(&ring, &[], &[], &[], 0.0, 0.4, centroid);
        assert!(super::hull_contains(&ring, point), "{point:?}");
        // The glyph sits in the middle of the band, wholly inside it.
        let distance = point.distance(centroid);
        assert!((2.4..=2.6).contains(&distance), "{distance}");
    }

    use super::super::UNIT;

    use super::*;

    #[test]
    fn self_loop_tube_keeps_the_interior_hole() {
        let curve = Circle::new((0.0, 0.0), 2.0).to_path(0.001);
        let region = union(&[&curve], &[Point::new(2.0, 0.0)], 0.4).unwrap();
        assert_eq!(
            region
                .iter()
                .filter(|element| matches!(element, PathEl::ClosePath))
                .count(),
            2
        );
        assert!(region.bounding_box().contains(Point::new(2.35, 0.0)));
    }
    #[test]
    fn hull_island_inside_a_tube_preserves_the_hole() {
        let center = Point::new(0.003, 0.005);
        let tube = Circle::new(center, 2.0).to_path(0.0001);
        let radius = 4.2 / UNIT;
        let angle = 39_f64.to_radians();
        let island = center + Vec2::new(angle.cos(), angle.sin()) * (2.0 - 2.0 * radius - 0.018);
        let hull = union(&[&tube], &[island], radius).unwrap();
        let mut rings = Vec::new();
        let mut ring = BezPath::new();
        for element in hull {
            ring.push(element);
            if element == PathEl::ClosePath {
                rings.push(ring);
                ring = BezPath::new();
            }
        }
        assert_eq!(rings.len(), 3);
        for (index, a) in rings.iter().enumerate() {
            for b in &rings[..index] {
                for segment in a.segments() {
                    let PathSeg::Line(line) = segment else {
                        panic!("hulls use line segments");
                    };
                    assert!(b
                        .segments()
                        .all(|segment| segment.intersect_line(line).is_empty()));
                }
            }
        }
        let witness = Point::new(16.90222058823529 / UNIT, 14.943926470588238 / UNIT);
        assert_eq!(
            rings.iter().filter(|ring| ring.contains(witness)).count(),
            2
        );
    }

    #[test]
    fn hull_outline_preserves_clearance_from_nonmember_ink() {
        let tube = Circle::new((0.003, 0.005), 2.0).to_path(0.0001);
        let hull = union(&[&tube], &[Point::new(4.0, 1.0)], 4.2 / UNIT).unwrap();
        let center = Point::new(-18.67131753365224 / UNIT, 10.65687952542826 / UNIT);
        let distance = hull
            .segments()
            .map(|segment| segment.nearest(center, 1e-9).distance_sq)
            .fold(f64::INFINITY, f64::min)
            .sqrt()
            * UNIT;
        assert!((distance - 1.2796154531954695).abs() < 1e-9, "{distance}");
        // Even the ends of the 0.1 pt non-member segment clear both strokes.
        assert!(distance - 0.05 > (1.1 + 1.0) / 2.0);
    }

    #[test]
    fn hull_diagonal_tubes_stay_connected_or_report_the_resolution_limit() {
        for (length, radius, resolvable) in [
            (1.0, 0.02, true),
            (120.0, 0.1, false),
            (1.0, 1e-6, false),
            (0.0001, 1e-6, true),
        ] {
            let mut curve = BezPath::new();
            curve.move_to((0.023, 0.031));
            curve.line_to((0.023 + length, 0.031 + length));
            let result = union(&[&curve], &[], radius);
            if resolvable {
                assert_eq!(
                    result
                        .unwrap()
                        .iter()
                        .filter(|el| matches!(el, PathEl::ClosePath))
                        .count(),
                    1
                );
            } else {
                assert!(result
                    .unwrap_err()
                    .contains("increase hull_radius or select a smaller region"));
            }
        }
    }
}
