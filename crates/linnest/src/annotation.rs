//! Label candidates built in Rust as `draw.typ` builds them: the native drawing
//! places its labels with them, and the label room check predicts a layout's
//! placement with them. Candidates exchange the same records as the Typst
//! renderer through the shared annotation search.
use kurbo::{BezPath, Shape};

use crate::curves;
pub(crate) use crate::label_placement::Bounds;
use crate::label_placement::{
    CandidateAttachment, CandidateSpec, CubicSegment, Frame, Point as PlanePoint, Transform,
};

/// Clearance between an external label and its free endpoint, as
/// `draw(external-label-gap:)` sets it.
pub(crate) const EXTERNAL_GAP: f64 = 0.25;
/// The smallest vertex radius, where edges stop, as `draw(node-min-radius:)`
/// sets it.
pub(crate) const NODE_RADIUS: f64 = 0.16;

/// A momentum arrow beside its carrier and the label riding it: the arrow's
/// offset, the visible length and ratio cap of its window, and the label's gap.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct ArrowStyle {
    pub offset: f64,
    pub window: (f64, f64),
    pub gap: f64,
}

fn wire(path: &BezPath) -> Vec<CubicSegment> {
    curves::cubics(path)
        .iter()
        .map(|cubic| {
            CubicSegment::new([cubic.p0, cubic.p1, cubic.p2, cubic.p3].map(|p| [p.x, p.y]))
        })
        .collect()
}

/// A batch of candidates at `positions` along a carrier with `frames`, as
/// `draw.typ` requests them, on an untransformed canvas.
fn spec(frames: Vec<curves::Frame>, positions: Vec<f64>, corners: [[f64; 3]; 4]) -> CandidateSpec {
    CandidateSpec {
        frames: frames
            .into_iter()
            .map(|frame| Frame {
                point: frame.point,
                tangent: frame.tangent,
            })
            .collect(),
        positions,
        transform: Transform::IDENTITY,
        origin: [0.0; 3],
        corners,
        total: 0.0,
        preferred: 0.0,
        accuracy: curves::ACCURACY,
        side: 1.0,
        preferred_side: 1.0,
        path_index: 0,
        path_shift: 0.0,
        gap: 0.0,
        clear_box: true,
        fixed: false,
        endpoint_direction: None,
        attachment: None,
    }
}

/// A placed candidate: the shared placement's record and, for a momentum
/// label, its arrow. Label room reads the text and arrow boxes and the cost;
/// the native drawing reads every field.
#[derive(Clone)]
#[cfg_attr(not(feature = "svg"), allow(dead_code))]
pub(crate) struct Candidate {
    pub bounds: Bounds,
    pub cost: f64,
    /// The distance along the carrier from its preferred place, relative to
    /// its length, that the cost weighs.
    pub relative: f64,
    pub corners: [[f64; 3]; 4],
    /// The frame point on the carrier: where a leader meets it, leaving
    /// along the normal toward the candidate's side.
    pub anchor: [f64; 2],
    pub normal: Option<[f64; 2]>,
    pub side: f64,
    pub arrow_bounds: Vec<Bounds>,
    pub arrow_mark: Option<usize>,
}

fn candidates(spec: &CandidateSpec) -> Result<Vec<Candidate>, String> {
    Ok(spec
        .candidates()?
        .into_iter()
        .map(|candidate| Candidate {
            bounds: candidate.bounds,
            cost: candidate.cost,
            relative: candidate.relative,
            corners: candidate.corners,
            anchor: candidate.anchor,
            normal: candidate.normal().map(|PlanePoint(x, y)| [x, y]),
            side: candidate.side,
            arrow_bounds: Vec::new(),
            arrow_mark: None,
        })
        .collect())
}

#[derive(Clone, Copy)]
pub(crate) enum Anchor {
    Center,
    East,
    West,
    North,
    South,
}

impl Anchor {
    /// The side of an external label facing away from its free endpoint.
    pub(crate) fn outward(dx: f64, dy: f64) -> Self {
        if dx < -1e-9 {
            Self::East
        } else if dx > 1e-9 {
            Self::West
        } else if dy > 0.0 {
            Self::South
        } else if dy < 0.0 {
            Self::North
        } else {
            Self::Center
        }
    }

    /// CeTZ `content` corners (NW, NE, SW, SE) of a box with half extents
    /// `[w, h]`, relative to its placement point at this anchor.
    pub(crate) fn corners(self, [w, h]: [f64; 2]) -> [[f64; 3]; 4] {
        let offset = match self {
            Self::East => [-w, 0.0],
            Self::West => [w, 0.0],
            Self::North => [0.0, -h],
            Self::South => [0.0, h],
            Self::Center => [0.0, 0.0],
        };
        [[-w, h], [w, h], [-w, -h], [w, -h]].map(|[x, y]| [x + offset[0], y + offset[1], 0.0])
    }
}

/// The candidates a label takes: every place `draw.typ` tries, or only its
/// preferred place on its preferred side, the first of every.
#[derive(Clone, Copy, PartialEq)]
pub(crate) enum Places {
    Every,
    Preferred,
}

impl Places {
    /// The sides tried, the preferred `side` first.
    fn sides(self, side: f64) -> &'static [f64] {
        match (self, side < 0.0) {
            (Self::Every, false) => &[1.0, -1.0],
            (Self::Every, true) => &[-1.0, 1.0],
            (Self::Preferred, false) => &[1.0],
            (Self::Preferred, true) => &[-1.0],
        }
    }
}

/// The preferred position and the 3/32..29/32 arc positions draw.typ tries,
/// clamped into `[low, high]`, deduplicated and away from the preferred one.
fn positions(total: f64, preferred: f64, low: f64, high: f64, places: Places) -> Vec<f64> {
    let mut positions = vec![preferred];
    if places == Places::Preferred {
        return positions;
    }
    let mut seen = Vec::new();
    for i in 3..30 {
        let at = (total * f64::from(i) / 32.0).clamp(low, high);
        if seen.contains(&at) {
            continue;
        }
        seen.push(at);
        if (at - preferred).abs() > curves::ACCURACY {
            positions.push(at);
        }
    }
    positions
}

/// A particle label sliding along its carrier at `gap`, on both sides,
/// interleaved by position, the preferred `side` first.
pub(crate) fn carrier_candidates(
    path: &BezPath,
    corners: [[f64; 3]; 4],
    gap: f64,
    side: f64,
    places: Places,
) -> Result<Vec<Candidate>, String> {
    let total = curves::length(path);
    let preferred = (total / 2.0).clamp(0.0, total);
    let positions = positions(total, preferred, 0.0, total, places);
    let frames = curves::frames(path, positions.clone())?;
    let sides = places
        .sides(side)
        .iter()
        .map(|&candidate_side| {
            candidates(&CandidateSpec {
                total,
                preferred,
                side: candidate_side,
                preferred_side: side,
                path_index: usize::from(candidate_side != side),
                gap,
                attachment: Some(CandidateAttachment {
                    carrier: wire(path),
                    paths: positions.iter().map(|_| Vec::new()).collect(),
                }),
                ..spec(frames.clone(), positions.clone(), corners)
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok((0..positions.len())
        .flat_map(|position| sides.iter().map(move |side| side[position].clone()))
        .collect())
}

/// A momentum label riding its arrow: per side, an offset carrier and an
/// arrow window at each position. Sides are not interleaved. The windows'
/// painted footprints resolve the candidates: the label clears them, and the
/// search weighs them against other annotations.
pub(crate) struct MomentumCandidates {
    specs: Vec<CandidateSpec>,
    windows: Vec<BezPath>,
}

impl MomentumCandidates {
    pub(crate) fn new(
        path: &BezPath,
        corners: [[f64; 3]; 4],
        side: f64,
        arrow: ArrowStyle,
        places: Places,
    ) -> Result<Self, String> {
        let mut specs = Vec::new();
        let mut windows = Vec::new();
        for (index, &carrier_side) in places.sides(side).iter().enumerate() {
            let carrier = curves::layer(path, arrow.offset * carrier_side, None, 0.0)?;
            let total = curves::length(&carrier);
            let centered = curves::layer(&carrier, 0.0, Some(arrow.window), 0.0)?;
            let half = curves::length(&centered) / 2.0;
            let low = half.clamp(0.0, total);
            let high = (total - half).clamp(low, total);
            let preferred = (total / 2.0).clamp(low, high);
            let positions = positions(total, preferred, low, high, places);
            let frames = curves::frames(&carrier, positions.clone())?;
            let shifts: Vec<f64> = positions.iter().map(|at| at - total / 2.0).collect();
            windows.extend(curves::layers(&carrier, &shifts, 0.0, Some(arrow.window))?);
            specs.push(CandidateSpec {
                total,
                preferred,
                side: carrier_side,
                preferred_side: side,
                path_index: index,
                gap: arrow.gap,
                attachment: Some(CandidateAttachment {
                    carrier: wire(path),
                    paths: Vec::new(),
                }),
                ..spec(frames, positions, corners)
            });
        }
        Ok(Self { specs, windows })
    }

    /// The arrow window of every candidate, in candidate order.
    pub(crate) fn windows(&self) -> &[BezPath] {
        &self.windows
    }

    /// Candidates clearing each window's painted arrow, in window order: per
    /// window, the path the label clears and the arrow's box.
    pub(crate) fn resolve(self, arrows: &[(BezPath, Bounds)]) -> Result<Vec<Candidate>, String> {
        if arrows.len() != self.windows.len() {
            return Err("momentum arrows must match their windows".to_owned());
        }
        let mut all = Vec::new();
        let mut arrows = arrows.iter();
        for mut spec in self.specs {
            let painted: Vec<_> = arrows.by_ref().take(spec.frames.len()).collect();
            spec.attachment.as_mut().expect("momentum attachment").paths =
                painted.iter().map(|(path, _)| wire(path)).collect();
            let mut resolved = candidates(&spec)?;
            for (candidate, (_, bounds)) in resolved.iter_mut().zip(painted) {
                candidate.arrow_bounds = vec![*bounds];
            }
            all.extend(resolved);
        }
        Ok(all)
    }
}

/// A fixed external label, cleared outward from its free endpoint.
pub(crate) fn endpoint_candidates(
    at: [f64; 2],
    direction: [f64; 2],
    corners: [[f64; 3]; 4],
) -> Result<Vec<Candidate>, String> {
    candidates(&CandidateSpec {
        gap: EXTERNAL_GAP,
        clear_box: false,
        fixed: true,
        endpoint_direction: Some(PlanePoint(direction[0], direction[1])),
        ..spec(
            vec![curves::Frame {
                point: at,
                tangent: [1.0, 0.0],
            }],
            vec![0.0],
            corners,
        )
    })
}

impl Bounds {
    /// The axis-aligned bounds of a painted path.
    pub(crate) fn painted(path: &BezPath) -> Self {
        let rect = path.bounding_box();
        Self {
            left: rect.x0,
            right: rect.x1,
            bottom: rect.y0,
            top: rect.y1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn straight(from: (f64, f64), to: (f64, f64)) -> BezPath {
        let mut path = BezPath::new();
        path.move_to(from);
        path.line_to(to);
        path
    }

    #[test]
    fn anchors_place_text_corners_beside_their_point() {
        let half = [2.0, 1.0];
        let corners = |anchor: Anchor| anchor.corners(half).map(|[x, y, _]| [x, y]);
        assert_eq!(
            corners(Anchor::Center),
            [[-2.0, 1.0], [2.0, 1.0], [-2.0, -1.0], [2.0, -1.0]]
        );
        // An external label faces away from its free end: text beside it.
        assert_eq!(corners(Anchor::East)[1], [0.0, 1.0]);
        assert_eq!(corners(Anchor::West)[0], [0.0, 1.0]);
        assert_eq!(corners(Anchor::North)[0], [-2.0, 0.0]);
        assert_eq!(corners(Anchor::South)[2], [-2.0, 0.0]);
        assert!(matches!(Anchor::outward(-1.0, 0.0), Anchor::East));
        assert!(matches!(Anchor::outward(0.0, 1.0), Anchor::South));
        assert!(matches!(Anchor::outward(0.0, 0.0), Anchor::Center));
    }

    #[test]
    fn places_start_at_the_preferred_one_and_stay_in_range() {
        let every = positions(32.0, 16.0, 4.0, 28.0, Places::Every);
        assert_eq!(every[0], 16.0);
        assert!(every.iter().all(|&at| (4.0..=28.0).contains(&at)));
        // Clamped stations collapse, and none repeats the preferred one.
        let mut sorted = every.clone();
        sorted.sort_by(f64::total_cmp);
        sorted.dedup();
        assert_eq!(sorted.len(), every.len());
        assert_eq!(positions(32.0, 16.0, 4.0, 28.0, Places::Preferred), [16.0]);
        assert_eq!(Places::Every.sides(-1.0), [-1.0, 1.0]);
        assert_eq!(Places::Preferred.sides(1.0), [1.0]);
    }

    #[test]
    fn carrier_candidates_alternate_sides_from_the_preferred_one() {
        let path = straight((0.0, 0.0), (4.0, 0.0));
        let corners = Anchor::Center.corners([0.5, 0.2]);
        let candidates = carrier_candidates(&path, corners, 0.15, -1.0, Places::Every).unwrap();
        // Both sides at every station, the preferred side first; the first
        // candidate is the preferred place.
        assert_eq!(candidates.len() % 2, 0);
        assert!(candidates[0].side < 0.0 && candidates[1].side > 0.0);
        assert_eq!(candidates[0].relative, 0.0);
        assert_eq!(candidates[0].anchor, [2.0, 0.0]);
        for pair in candidates.chunks(2) {
            assert_eq!(pair[0].anchor, pair[1].anchor);
            // Text clears the carrier by the gap on its side.
            assert!(
                (pair[0].bounds.top + 0.15).abs() < 1e-9,
                "{:?}",
                pair[0].bounds
            );
            assert!(
                (pair[1].bounds.bottom - 0.15).abs() < 1e-9,
                "{:?}",
                pair[1].bounds
            );
            assert_eq!(pair[0].normal, Some([0.0, -1.0]));
        }
        // Costs grow away from the preferred place, as its relative distance.
        assert!(candidates
            .iter()
            .all(|candidate| candidate.cost >= candidates[0].cost));
        assert!(candidates
            .iter()
            .any(|candidate| candidate.relative.abs() > 0.2));
        let preferred = carrier_candidates(&path, corners, 0.15, 1.0, Places::Preferred).unwrap();
        assert_eq!(preferred.len(), 1);
        assert_eq!(preferred[0].bounds, candidates[1].bounds);
    }

    #[test]
    fn momentum_candidates_ride_one_arrow_window_each() {
        let path = straight((0.0, 0.0), (6.0, 0.0));
        let style = ArrowStyle {
            offset: 0.35,
            window: (1.4, 0.5),
            gap: 0.2,
        };
        let corners = Anchor::Center.corners([0.5, 0.2]);
        let proposals = MomentumCandidates::new(&path, corners, 1.0, style, Places::Every).unwrap();
        let windows = proposals.windows().to_vec();
        assert!(proposals.resolve(&[]).is_err());
        let proposals = MomentumCandidates::new(&path, corners, 1.0, style, Places::Every).unwrap();
        let arrows: Vec<_> = windows
            .iter()
            .map(|window| (window.clone(), Bounds::painted(window)))
            .collect();
        let candidates = proposals.resolve(&arrows).unwrap();
        assert_eq!(candidates.len(), windows.len());
        for (candidate, (window, bounds)) in candidates.iter().zip(&arrows) {
            assert_eq!(candidate.arrow_bounds, [*bounds]);
            // Each label clears its own arrow by the gap and anchors on it.
            let y = window.elements()[0].end_point().unwrap().y;
            assert!((candidate.anchor[1] - y).abs() < 1e-9);
            assert!(
                candidate.bounds.bottom >= y + 0.2 - 1e-9 || candidate.bounds.top <= y - 0.2 + 1e-9
            );
        }
        // The preferred side rides the arrow above the carrier first.
        assert!(candidates[0].anchor[1] > 0.0 && candidates[0].side > 0.0);
    }
}
