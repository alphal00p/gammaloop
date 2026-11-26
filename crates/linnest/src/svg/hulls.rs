//! Sampled unions of round tubes and vertex discs.

use std::collections::{BTreeMap, BTreeSet};

use kurbo::{flatten, BezPath, Line, ParamCurveNearest, PathEl, Point, Rect, Shape};

type Sample = (i64, i64);

// Every marching-square crossing is indexed by the doubled integer position
// of its grid-edge midpoint. Adjacent cells therefore join exactly.
pub(super) fn union(
    paths: &[&BezPath],
    vertices: &[Point],
    radius: f64,
) -> Result<BezPath, String> {
    const MAX_CELLS_PER_SIDE: f64 = 768.0;
    const RESOLUTION_ERROR: &str =
        "hull is too thin for its drawing extent; increase hull_radius or select a smaller region";
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
    if !bounds.is_finite() || !step.is_finite() || max_coordinate / step >= (i64::MAX / 4) as f64 {
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
            // Soften the grid staircase before simplifying the closed path.
            let count = samples.len();
            let softened: Vec<_> = (0..count)
                .map(|i| samples[i].midpoint(samples[(i + 1) % count]))
                .collect();
            let simplified = simplify(&softened, (radius / 12.0).min(0.03));
            let contour = if simplified.len() < 3 {
                &softened
            } else {
                &simplified
            };
            result.move_to(contour[0]);
            for &point in &contour[1..] {
                result.line_to(point);
            }
            result.close_path();
        }
    }
    if result.elements().is_empty() {
        return Err(RESOLUTION_ERROR.into());
    }
    Ok(result)
}

/// Ramer–Douglas–Peucker on both halves of a closed contour. The tolerance
/// removes collinear samples as well as deviations below 0.4 pt.
fn simplify(points: &[Point], tolerance: f64) -> Vec<Point> {
    let mut keep = vec![false; points.len() + 1];
    let middle = points.len() / 2;
    let mut spans = vec![(0, middle), (middle, points.len())];
    keep[0] = true;
    keep[middle] = true;
    let limit = tolerance * tolerance;
    while let Some((start, end)) = spans.pop() {
        let a = points[start];
        let b = points[end % points.len()];
        let chord = Line::new(a, b);
        let mut furthest = None;
        let mut max_distance = limit;
        for (index, &point) in points.iter().enumerate().take(end).skip(start + 1) {
            let distance = if a == b {
                a.distance_squared(point)
            } else {
                chord.nearest(point, tolerance).distance_sq
            };
            if distance > max_distance {
                max_distance = distance;
                furthest = Some(index);
            }
        }
        if let Some(index) = furthest {
            keep[index] = true;
            spans.extend([(start, index), (index, end)]);
        }
    }
    points
        .iter()
        .zip(keep)
        .filter_map(|(&point, keep)| keep.then_some(point))
        .collect()
}

#[cfg(test)]
mod tests {
    use kurbo::{Circle, PathEl};

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
