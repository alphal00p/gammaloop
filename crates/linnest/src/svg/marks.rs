//! Arrowheads in drawing units, sized like CeTZ marks on a 13.5pt canvas.
use kurbo::{BezPath, Point, Vec2};

use super::{curves, labels::Bounds};

/// A mark's length along its path and half its width across it.
#[derive(Clone, Copy)]
pub(super) struct Mark {
    pub length: f64,
    pub half_width: f64,
}

/// CeTZ's open "straight" momentum chevron at scale 0.8 (0.16cm by 0.12cm).
pub(super) const CHEVRON: Mark = Mark {
    length: 0.3359580052493438,
    half_width: 0.12598425196850394,
};

/// CeTZ's filled ">" particle-flow triangle at scale 1.05 (0.21cm by 0.1575cm).
pub(super) const TRIANGLE: Mark = Mark {
    length: 0.4409448818897638,
    half_width: 0.16535433070866143,
};

impl Mark {
    /// Wing, tip, wing for a mark whose tip and back lie at `tip` and `back`.
    fn outline(self, tip: Point, back: Point) -> [Point; 3] {
        let along = (back - tip).normalize();
        let across = Vec2::new(-along.y, along.x) * self.half_width;
        let base = tip + along * self.length;
        [base + across, tip, base - across]
    }

    /// The mark at a path's end, flexed so its back also lies on the path.
    pub(super) fn at_end(self, path: &BezPath) -> Option<[Point; 3]> {
        let mut points = Vec::new();
        kurbo::flatten(path, 1e-4, |element| match element {
            kurbo::PathEl::MoveTo(p) | kurbo::PathEl::LineTo(p) => points.push(p),
            _ => {}
        });
        let tip = *points.last()?;
        let mut back = *points.first()?;
        for window in points.windows(2).rev() {
            let (a, b) = (window[0], window[1]);
            if a.distance(tip) >= self.length {
                // Solve |b + (a - b) t - tip| = length on the crossing segment.
                let (d, f) = (a - b, b - tip);
                let (qa, qb, qc) = (d.dot(d), 2.0 * d.dot(f), f.dot(f) - self.length.powi(2));
                let t = (-qb + (qb * qb - 4.0 * qa * qc).max(0.0).sqrt()) / (2.0 * qa);
                back = b + d * t.clamp(0.0, 1.0);
                break;
            }
        }
        (back != tip).then(|| self.outline(tip, back))
    }

    /// The mark centered at an arc-length ratio, pointing along or against the path.
    pub(super) fn centered(self, path: &BezPath, ratio: f64, forward: bool) -> Option<[Point; 3]> {
        let total = curves::length(path);
        let at = total * ratio.clamp(0.0, 1.0);
        let (center, _) = curves::point_at(path, at)?;
        let (behind, _) = curves::point_at(path, (at - self.length / 2.0).max(0.0))?;
        let (ahead, _) = curves::point_at(path, (at + self.length / 2.0).min(total))?;
        let along = if forward {
            ahead - behind
        } else {
            behind - ahead
        };
        if along.hypot() <= 1e-12 {
            return None;
        }
        let half = along.normalize() * (self.length / 2.0);
        Some(self.outline(center + half, center - half))
    }
}

/// Bounds of points, padded by a stroke radius.
pub(super) fn footprint(points: &[Point], radius: f64) -> Bounds {
    let mut bounds = Bounds::EMPTY;
    for point in points {
        bounds = bounds.union(Bounds::point(*point));
    }
    bounds.padded(radius)
}
