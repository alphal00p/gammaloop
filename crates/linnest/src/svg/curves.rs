//! Kurvst geometry with the arguments `draw.typ` passes.
use kurbo::{BezPath, CubicBez, ParamCurveArclen, Point};
use kurvst::{
    CubicBezierSpec, CurvePoint, FittedCoilInput, HobbySplineSpec, ParallelPathSpec,
    PathFramesSpec, PathTrimmer, PatternInput, PatternPathSpec, RegionPart, RegionSamplesSpec,
};

/// Kurvst's arc-length accuracy for drawing.
pub(super) const ACCURACY: f64 = 0.001;

/// A point and unnormalized tangent on a path.
#[derive(Clone, Copy)]
pub(super) struct Frame {
    pub point: [f64; 2],
    pub tangent: [f64; 2],
}

pub(super) fn length(path: &BezPath) -> f64 {
    path.segments()
        .map(|segment| segment.arclen(ACCURACY))
        .sum()
}

pub(super) fn cubics(path: &BezPath) -> Vec<CubicBez> {
    path.segments().map(|segment| segment.to_cubic()).collect()
}

/// `curve-api.path(..segments.map(from-cubic))`: contiguous cubics share a move.
pub(super) fn from_cubics(segments: &[CubicBez]) -> BezPath {
    let mut path = BezPath::new();
    let mut previous = None;
    for segment in segments {
        // Offset cubics can differ at shared knots by floating-point roundoff.
        // Keep those knots connected without bridging genuinely separate parts.
        if previous.is_none_or(|end: Point| end.distance(segment.p0) > 1e-9) {
            path.move_to(segment.p0);
        }
        path.curve_to(segment.p1, segment.p2, segment.p3);
        previous = Some(segment.p3);
    }
    path
}

/// The arc-length window `[start, length - end]`; empty when nothing remains.
pub(super) fn trim(path: &BezPath, start: f64, end: f64) -> Result<BezPath, String> {
    PathTrimmer::new(path.segments(), ACCURACY)
        .trim(start, end)
        .map(|segments| BezPath::from_path_segments(segments.into_iter()))
}

/// One measured path trimmed to every `(start, end)` outset window.
pub(super) fn windows(path: &BezPath, outsets: &[(f64, f64)]) -> Result<Vec<BezPath>, String> {
    let trimmer = PathTrimmer::new(path.segments(), ACCURACY);
    outsets
        .iter()
        .map(|&(start, end)| {
            trimmer
                .trim(start, end)
                .map(|segments| BezPath::from_path_segments(segments.into_iter()))
        })
        .collect()
}

/// `_trim-routed-path`: empty when the outsets consume the path.
pub(super) fn trim_routed(path: &BezPath, start: f64, end: f64) -> Result<BezPath, String> {
    if path.elements().is_empty() || (start == 0.0 && end == 0.0) {
        return Ok(path.clone());
    }
    if start + end >= length(path) {
        return Ok(BezPath::new());
    }
    trim(path, start, end)
}

/// `_routed-split-through(..).curve`: drop solver-collapsed knots, then Hobby.
pub(super) fn routed_curve(points: &[[f64; 2]]) -> Result<BezPath, String> {
    let mut knots: Vec<CurvePoint> = Vec::new();
    for &[x, y] in points {
        if knots
            .last()
            .is_none_or(|last| Point::new(last.x, last.y).distance(Point::new(x, y)) > f64::EPSILON)
        {
            knots.push(CurvePoint { x, y });
        }
    }
    if knots.len() < 2 {
        return Ok(BezPath::new());
    }
    HobbySplineSpec {
        points: knots,
        omega: 1.0,
        accuracy: ACCURACY,
    }
    .curve()
    .map(|output| output.path)
}

/// `kurvst.parallel`: positive distances lie left of the path direction.
pub(super) fn parallel(path: &BezPath, distance: f64) -> Result<BezPath, String> {
    ParallelPathSpec {
        path: path.clone(),
        distance,
        start_outset: 0.0,
        end_outset: 0.0,
        accuracy: ACCURACY,
    }
    .parallel()
    .map(|output| from_cubics(&cubics(&output.path)))
}

/// Kurvst's Typst `layers`: offset, then shifted windows of `length`/`ratio`
/// resolved by "min", from one trim batch. `None` keeps the full path.
pub(super) fn layers(
    path: &BezPath,
    shifts: &[f64],
    offset: f64,
    window: Option<(f64, f64)>,
) -> Result<Vec<BezPath>, String> {
    let path = if offset != 0.0 {
        parallel(path, offset)?
    } else {
        path.clone()
    };
    if path.segments().next().is_none() {
        return Ok(shifts.iter().map(|_| path.clone()).collect());
    }
    let center_trim = window.map_or(0.0, |(fixed, ratio)| {
        let total = length(&path);
        ((total - fixed.min(total * ratio)) / 2.0).max(0.0)
    });
    let outsets: Vec<_> = shifts
        .iter()
        .map(|&shift| {
            let shift = (-center_trim).max(center_trim.min(shift));
            (center_trim + shift, center_trim - shift)
        })
        .collect();
    windows(&path, &outsets)
}

pub(super) fn layer(
    path: &BezPath,
    offset: f64,
    window: Option<(f64, f64)>,
    shift: f64,
) -> Result<BezPath, String> {
    Ok(layers(path, &[shift], offset, window)?.remove(0))
}

/// Points and tangents at arc distances, clamped onto the path.
pub(super) fn frames(path: &BezPath, distances: Vec<f64>) -> Result<Vec<Frame>, String> {
    let frames = PathFramesSpec {
        path: path.clone(),
        distances,
        accuracy: ACCURACY,
    }
    .frames()?;
    Ok(frames
        .into_iter()
        .map(|frame| {
            frame.map_or(
                Frame {
                    point: [0.0; 2],
                    tangent: [0.0; 2],
                },
                |frame| Frame {
                    point: [frame.point.x, frame.point.y],
                    tangent: [frame.tangent.x, frame.tangent.y],
                },
            )
        })
        .collect())
}

/// A physics line decoration over a whole edge.
pub(super) fn pattern(path: &BezPath, decoration: &super::Pattern) -> Result<BezPath, String> {
    use super::Pattern;
    let length = length(path);
    if length <= f64::EPSILON {
        return Ok(path.clone());
    }
    let spec = |pattern, amplitude, wavelength, longitudinal_scale| PatternPathSpec {
        path: path.clone(),
        pattern,
        amplitude,
        wavelength,
        phase: 0.0,
        samples_per_period: 16,
        coil_longitudinal_scale: longitudinal_scale,
        anchor_start: true,
        anchor_end: true,
        endpoint_slope: 0.0,
        split_at: Vec::new(),
        accuracy: ACCURACY,
    };
    // Coils fit their natural endpoints to the complete edge before painting.
    let spec = match *decoration {
        Pattern::Coil {
            amplitude,
            wavelength,
            longitudinal_scale,
        } => spec(
            PatternInput::FittedCoil(FittedCoilInput {
                kind: "fitted-coil".to_owned(),
                fit_length: length,
                amplitude,
                wavelength,
                longitudinal_scale,
                samples_per_period: 16,
            }),
            amplitude,
            length,
            longitudinal_scale,
        ),
        Pattern::Wave {
            amplitude,
            wavelength,
        } => spec(PatternInput::wave(16), amplitude, wavelength, 1.25),
        Pattern::Zigzag {
            amplitude,
            wavelength,
        } => spec(PatternInput::zigzag(), amplitude, wavelength, 1.25),
    };
    spec.patterned().map(|output| output.path)
}

/// Hover sample points per arc-length region (0..4) of the parts together.
pub(super) type RegionSamples = Vec<(usize, Vec<[f64; 2]>)>;

/// Hover sample points along each visible part, in four arc-length regions.
pub(super) fn region_samples(parts: &[&BezPath], unit: f64) -> Result<RegionSamples, String> {
    RegionSamplesSpec {
        parts: parts
            .iter()
            .map(|path| RegionPart {
                segments: cubics(path)
                    .into_iter()
                    .map(CubicBezierSpec::from)
                    .collect(),
                visible: true,
            })
            .collect(),
        regions: 4,
        unit,
        step: 4.0,
        accuracy: ACCURACY,
    }
    .samples()
    .map(|samples| {
        samples
            .into_iter()
            .map(|(region, points)| (region, points.into_iter().map(|p| [p.x, p.y]).collect()))
            .collect()
    })
}

/// Kurvst's Typst `cubic-tangent`, operation for operation.
pub(super) fn cubic_tangent(segment: &CubicBez, t: f64) -> [f64; 2] {
    let lerp = |a: Point, b: Point| a.lerp(b, t);
    let ab = lerp(segment.p0, segment.p1);
    let bc = lerp(segment.p1, segment.p2);
    let cd = lerp(segment.p2, segment.p3);
    let (abc, bcd) = (lerp(ab, bc), lerp(bc, cd));
    [bcd.x - abc.x, bcd.y - abc.y]
}

#[cfg(test)]
mod tests {
    use super::*;
    use kurbo::PathEl;

    #[test]
    fn windows_retain_connected_offset_join_below_cumulative_length_resolution() {
        // Actual offset carrier from the nested-combine native SVG scene.
        // Its middle cubic connects distinct endpoints, but its length is too
        // small to advance the cumulative arc cursor after the first cubic.
        let mut path = BezPath::new();
        path.move_to((-2.3093354422559713, 0.45125171352077403));
        path.curve_to(
            (-1.7364145273475438, 1.0067548034067653),
            (-0.9681573656345025, 1.3196773655325762),
            (-0.16598397554102728, 1.319677365532576),
        );
        path.curve_to(
            (-0.16598397554102728, 1.319677365532576),
            (-0.16598397554102734, 1.319677365532576),
            (-0.16598397554102734, 1.319677365532576),
        );
        path.curve_to(
            (0.6361894145524469, 1.319677365532576),
            (1.4044465762654879, 1.0067548034067655),
            (1.9773674911739154, 0.4512517135207751),
        );
        let outset = 1.6707269964337421;
        let window = windows(&path, &[(outset, outset)]).unwrap().remove(0);
        assert_eq!(window.segments().count(), 3);
        assert_eq!(
            window
                .elements()
                .iter()
                .filter(|el| matches!(el, PathEl::MoveTo(_)))
                .count(),
            1,
        );
        assert_eq!(window.elements()[2], path.elements()[2]);
    }
}
