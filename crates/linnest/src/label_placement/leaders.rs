//! Leaders: a label whose chosen place still meets a painted line, another
//! label, a vertex or an arrow sits in free space nearby instead, with a
//! hairline back to its carrier. Layout room clears what it can first; a
//! leader is the drawing's last resort, and a label without a qualifying
//! leader keeps its place, where the label report records what it meets.
use serde::{Deserialize, Serialize};

use super::{
    same_edge, segment_gap, segments_meet, Bounds, CollisionSpec, LabelOverlaps, LineGrid, Point,
    StrokeLine, TextFrame, Transform,
};

/// What a leadered label keeps clear, in canvas units: its text from painted
/// lines and from arrows; its text, leader and anchor from vertices; and its
/// whole leader from other edges' lines, other arrows and other text, and,
/// once risen off it, from its own carrier. Other text keeps one line height
/// from leadered text.
const TEXT_LINE_CLEARANCE: f64 = 0.15;
const TEXT_ARROW_CLEARANCE: f64 = 0.1;
pub(super) const VERTEX_CLEARANCE: f64 = 0.3;
const LEADER_CLEARANCE: f64 = 0.1;
/// A leader leaves its anchor, on the carrier or the middle of the momentum
/// arrow, steeply toward the label's side, so it never runs alongside an
/// edge: the least sine of its angle to the carrier, and how far it rises
/// off the carrier, past the band or arrow it starts on, before it must
/// clear its own carrier's painted lines.
pub(super) const LEADER_DEPARTURE: f64 = 0.6;
const LEADER_RISE: f64 = 0.3;
/// The gap between a leader and its text, and the shortest and longest leader.
const LEADER_GAP: f64 = 0.08;
const LEADER_SHORTEST: f64 = 0.3;
const LEADER_LONGEST: f64 = 6.0;
/// The polar search around each anchor: ring step, directions, the spacing
/// of anchors tried along a carrier, and how far text may leave the drawing.
const RING_STEP: f64 = 0.2;
const DIRECTIONS: u32 = 36;
const ANCHOR_SPACING: f64 = 0.4;
const OUTSIDE_REACH: f64 = 3.0;
/// Cost beyond the leader's length: per unit of anchor distance from the
/// carrier's preferred place, relative to the carrier; per unit of painted
/// line within one unit of the text, preferring open space where the label
/// reads on its own; and per unit the text extends the drawing.
const ANCHOR_WEIGHT: f64 = 3.0;
const CONGESTION_WEIGHT: f64 = 0.25;
const OUTSIDE_WEIGHT: f64 = 1.5;

/// Leaders as the host draws them.
#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
pub(super) struct LeaderSpec {
    /// One line of text: the clearance between leadered and other text.
    #[serde(deserialize_with = "crate::deserialize_f64")]
    line_height: f64,
    /// The canvas transform of local drawing coordinates; identity if absent.
    #[serde(default)]
    transform: Option<Transform>,
    /// The edge of each painted line: a label's own carrier is the lines of
    /// its edge.
    #[serde(default)]
    line_edges: Vec<ciborium::Value>,
}

/// A label set off its carrier, beside free space, with a leader back to the
/// anchor of one of its candidates: a frame point on its carrier or momentum
/// arrow, whose arrow is drawn with it.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(rename_all = "kebab-case")]
pub(crate) struct Leader {
    pub(crate) choice: usize,
    /// The text center, where the label is drawn, in local coordinates.
    pub(crate) position: Point,
    /// The leader from the anchor to the text, in local coordinates.
    pub(crate) from: Point,
    pub(crate) to: Point,
    /// The moved text box and, in canvas coordinates, the leader from its
    /// anchor, past where it has risen off its carrier, to its end: what
    /// must clear painted lines and other text.
    pub(crate) bounds: Bounds,
    #[serde(deserialize_with = "Leader::deserialize_corners")]
    pub(crate) corners: [Point; 4],
    pub(crate) line: [Point; 3],
}

impl Leader {
    fn deserialize_corners<'de, D: serde::Deserializer<'de>>(
        de: D,
    ) -> Result<[Point; 4], D::Error> {
        use serde::de::Error;
        Point::deserialize_quad(de)?.ok_or_else(|| D::Error::custom("leader corners are required"))
    }
}

impl Bounds {
    fn around(point: Point, reach: f64) -> Self {
        Self {
            left: point.0 - reach,
            right: point.0 + reach,
            bottom: point.1 - reach,
            top: point.1 + reach,
        }
    }

    fn shifted(self, delta: Point) -> Self {
        Self {
            left: self.left + delta.0,
            right: self.right + delta.0,
            bottom: self.bottom + delta.1,
            top: self.top + delta.1,
        }
    }

    fn contains(self, other: Self) -> bool {
        self.left <= other.left
            && other.right <= self.right
            && self.bottom <= other.bottom
            && other.top <= self.top
    }

    /// How far `self` reaches beyond `extent`, summed over its four sides.
    fn excess(self, extent: Self) -> f64 {
        (extent.left - self.left).max(0.0)
            + (self.right - extent.right).max(0.0)
            + (extent.bottom - self.bottom).max(0.0)
            + (self.top - extent.top).max(0.0)
    }
}

impl Transform {
    /// The local point drawn at canvas `point`. The z column and row play no
    /// part for drawing-plane points.
    fn local(&self, Point(x, y): Point) -> Point {
        let [[a, b, _, e], [c, d, _, f], _, _] = self.0;
        let determinant = a * d - b * c;
        if determinant.abs() <= 1e-18 {
            return Point(x, y);
        }
        let (x, y) = (x - e, y - f);
        Point((d * x - b * y) / determinant, (a * y - c * x) / determinant)
    }
}

/// The cell of the blocking raster, in canvas units.
const CELL: f64 = 0.1;

/// An occupancy raster over a region, with prefix sums that count the marked
/// cells under an axis-aligned box in constant time.
struct Occupancy {
    region: Bounds,
    columns: usize,
    rows: usize,
    cells: Vec<u8>,
    sums: Vec<u32>,
}

impl Occupancy {
    fn new(region: Bounds) -> Self {
        let count = |span: f64| ((span / CELL).ceil() as usize).max(1);
        let (columns, rows) = (
            count(region.right - region.left),
            count(region.top - region.bottom),
        );
        Self {
            region,
            columns,
            rows,
            cells: vec![0; columns * rows],
            sums: Vec::new(),
        }
    }

    /// The inclusive cells a box reaches or, with `inside`, lies over
    /// entirely; `None` when there are none.
    fn span(&self, bounds: Bounds, inside: bool) -> Option<([usize; 2], [usize; 2])> {
        if !bounds.is_finite() {
            return None;
        }
        // Cell coordinates clamped into the raster, where truncation floors.
        let range = |low: f64, high: f64, origin: f64, count: usize| {
            let cells = count as f64;
            let (low, high) = ((low - origin) / CELL, (high - origin) / CELL);
            if high < 0.0 || low >= cells {
                return None;
            }
            let (low, high) = (low.clamp(0.0, cells), high.clamp(0.0, cells));
            let (first, last) = (low as usize, high as usize);
            let (first, last) = if inside {
                // Cells wholly inside: from the ceiling of `low` below `high`'s floor.
                (
                    first + usize::from((first as f64) < low),
                    last.checked_sub(1)?,
                )
            } else {
                (first, last)
            };
            (first <= last && first < count).then(|| [first, last.min(count - 1)])
        };
        Some((
            range(bounds.left, bounds.right, self.region.left, self.columns)?,
            range(bounds.bottom, bounds.top, self.region.bottom, self.rows)?,
        ))
    }

    fn mark(&mut self, bounds: Bounds, inside: bool) {
        if let Some(([left, right], [bottom, top])) = self.span(bounds, inside) {
            for row in bottom..=top {
                self.cells[row * self.columns + left..=row * self.columns + right].fill(1);
            }
        }
    }

    /// Samples of a stroke, half a cell apart: every point of the stroke is
    /// within a quarter cell of one.
    fn samples(line: &StrokeLine) -> impl Iterator<Item = Point> + '_ {
        let samples = (line.start.distance(line.end) / (CELL / 2.0))
            .ceil()
            .clamp(1.0, 1e6) as usize;
        (0..=samples).map(move |sample| line.start.lerp(line.end, sample as f64 / samples as f64))
    }

    fn summed(mut self) -> Self {
        let width = self.columns + 1;
        self.sums = vec![0; width * (self.rows + 1)];
        for row in 0..self.rows {
            let mut line = 0;
            for column in 0..self.columns {
                line += u32::from(self.cells[row * self.columns + column]);
                self.sums[(row + 1) * width + column + 1] =
                    self.sums[row * width + column + 1] + line;
            }
        }
        self
    }

    /// Marked cells a box reaches.
    fn count(&self, bounds: Bounds) -> u32 {
        let Some(([left, right], [bottom, top])) = self.span(bounds, false) else {
            return 0;
        };
        let width = self.columns + 1;
        let at = |row: usize, column: usize| self.sums[row * width + column];
        at(top + 1, right + 1) + at(bottom, left) - at(bottom, right + 1) - at(top + 1, left)
    }
}

/// Where blocking geometry lies: `reach` marks every cell it may reach and
/// `covers` only cells it surely covers. A box over no `reach` cell is clear,
/// one over a `covers` cell is blocked; only boxes in between need exact tests.
struct Blocking {
    reach: Occupancy,
    covers: Occupancy,
}

#[derive(Clone, Copy, PartialEq)]
enum Clearance {
    Clear,
    Blocked,
    Unsure,
}

impl Blocking {
    fn new(region: Bounds) -> Self {
        Self {
            reach: Occupancy::new(region),
            covers: Occupancy::new(region),
        }
    }

    fn mark(&mut self, bounds: Bounds) {
        self.reach.mark(bounds, false);
        self.covers.mark(bounds, true);
    }

    /// Mark the cells within `clearance` of a stroke: those it may reach, and
    /// those wholly within reach of one of its samples.
    fn mark_line(&mut self, line: &StrokeLine, clearance: f64) {
        let reach = line.radius + clearance;
        for at in Occupancy::samples(line) {
            self.reach
                .mark(Bounds::around(at, reach + CELL / 4.0), false);
            self.covers
                .mark(Bounds::around(at, reach / std::f64::consts::SQRT_2), true);
        }
    }

    fn summed(self) -> Self {
        Self {
            reach: self.reach.summed(),
            covers: self.covers.summed(),
        }
    }

    fn clearance(&self, bounds: Bounds) -> Clearance {
        if self.reach.count(bounds) == 0 {
            Clearance::Clear
        } else if self.covers.count(bounds) > 0 {
            Clearance::Blocked
        } else {
            Clearance::Unsure
        }
    }
}

/// A placed annotation's text: its box and exact corners.
#[derive(Clone, Copy)]
struct Text {
    bounds: Bounds,
    corners: [Point; 4],
}

impl Text {
    fn shifted(self, delta: Point) -> Self {
        Self {
            bounds: self.bounds.shifted(delta),
            corners: self.corners.map(|corner| corner.add(delta)),
        }
    }

    fn center(self) -> Point {
        self.corners[0].lerp(self.corners[3], 0.5)
    }
}

impl LeaderSpec {
    /// Leaders for the two-sided labels whose chosen place meets something,
    /// as `found` lists it. In label order, each such label takes the
    /// cheapest spot around the anchor of any of its candidates whose arrow,
    /// if any, is clear, keeping every clearance above; the cost is the
    /// leader's length plus the anchor's distance from the carrier's
    /// preferred place, painted lines near the text, and growth of the
    /// drawing. Labels without such a spot keep their place.
    pub(super) fn place(
        &self,
        spec: &CollisionSpec,
        choices: &[usize],
        found: &[LabelOverlaps],
    ) -> Vec<Option<Leader>> {
        let mut leaders = vec![None; choices.len()];
        let meeting: Vec<usize> = found
            .iter()
            .enumerate()
            .filter(|&(index, found)| {
                found.contacts() > 0
                    && spec.placements[index].candidates.len() > 1
                    && spec.placements[index].two_sided()
            })
            .map(|(index, _)| index)
            .collect();
        if meeting.is_empty() {
            return leaders;
        }
        let lines = &spec.edge_lines;
        let nodes = spec.node_boxes();
        let mut texts: Vec<Text> = spec
            .placements
            .iter()
            .zip(choices)
            .map(|(label, &choice)| {
                let candidate = &label.candidates[choice];
                Text {
                    bounds: candidate.bounds,
                    corners: candidate.corners.unwrap_or(candidate.bounds.corners()),
                }
            })
            .collect();
        let mut arrows: Vec<Vec<Bounds>> = spec
            .placements
            .iter()
            .zip(choices)
            .map(|(label, &choice)| label.candidates[choice].arrow_bounds.clone())
            .collect();
        let Some(extent) = lines
            .iter()
            .map(StrokeLine::extent)
            .chain(nodes.iter().copied())
            .chain(texts.iter().map(|text| text.bounds))
            .chain(spec.fixed_arrows.iter().copied())
            .chain(arrows.iter().flatten().copied())
            .filter(|bounds| bounds.is_finite())
            .reduce(Bounds::union)
        else {
            return leaders;
        };
        let region = extent.expanded(OUTSIDE_REACH);
        let grid = LineGrid::new(lines);
        // What no leader moves, rasterized with its clearance from moved
        // text: painted lines, vertices, fixed arrows, and the labels that
        // keep their place with their arrows. Labels that may move are
        // tested exactly at every spot.
        let mut moving = vec![false; choices.len()];
        for &index in &meeting {
            moving[index] = true;
        }
        let mut fixed = Blocking::new(region);
        for line in lines {
            fixed.mark_line(line, TEXT_LINE_CLEARANCE);
        }
        for node in &nodes {
            fixed.mark(node.expanded(VERTEX_CLEARANCE));
        }
        for arrow in &spec.fixed_arrows {
            fixed.mark(arrow.expanded(TEXT_ARROW_CLEARANCE));
        }
        for (index, text) in texts.iter().enumerate() {
            if !moving[index] {
                fixed.mark(text.bounds.expanded(self.line_height));
                for arrow in &arrows[index] {
                    fixed.mark(arrow.expanded(TEXT_ARROW_CLEARANCE));
                }
            }
        }
        let fixed = fixed.summed();
        let directions: Vec<Point> = (0..DIRECTIONS)
            .map(|turn| {
                let angle = std::f64::consts::TAU * f64::from(turn) / f64::from(DIRECTIONS);
                Point(angle.cos(), angle.sin())
            })
            .collect();
        let local = |point: Point| {
            self.transform
                .as_ref()
                .map_or(point, |transform| transform.local(point))
        };
        // Leaders drawn so far, from anchor to end.
        let mut drawn: Vec<[Point; 2]> = Vec::new();
        for index in meeting {
            let label = &spec.placements[index];
            let text = texts[index];
            // The label's own carrier: its leader rises off it.
            let carrier = |line: usize| {
                label
                    .edge
                    .as_ref()
                    .zip(self.line_edges.get(line))
                    .is_some_and(|(edge, owner)| same_edge(edge, owner))
            };
            // Whether a segment keeps clear of every painted line, its own
            // carrier's once risen from `risen`.
            let lines_clear = |anchor: Point, risen: Point, end: Point| {
                grid.near(Bounds::spanning(anchor, end).expanded(1.0))
                    .into_iter()
                    .all(|line| {
                        let from = if carrier(line) { risen } else { anchor };
                        let line = &lines[line];
                        segment_gap(from, end, line.start, line.end)
                            >= line.radius + LEADER_CLEARANCE
                    })
            };
            let others: Vec<Text> = texts
                .iter()
                .enumerate()
                .filter_map(|(other, &text)| (other != index).then_some(text))
                .collect();
            let other_arrows: Vec<Bounds> = spec
                .fixed_arrows
                .iter()
                .chain(
                    arrows
                        .iter()
                        .enumerate()
                        .filter(|&(other, _)| other != index)
                        .flat_map(|(_, arrows)| arrows.iter()),
                )
                .copied()
                .collect();
            // Labels that may still move: their current text and arrows.
            let movable: Vec<(Bounds, f64)> = (0..choices.len())
                .filter(|&other| other != index && moving[other])
                .flat_map(|other| {
                    std::iter::once((texts[other].bounds, self.line_height)).chain(
                        arrows[other]
                            .iter()
                            .map(|&arrow| (arrow, TEXT_ARROW_CLEARANCE)),
                    )
                })
                .collect();
            let center = text.center();
            let reach = text.corners[0].distance(text.corners[3]) / 2.0 + LEADER_GAP;
            // The anchors of candidates whose arrows are clear, at least the
            // anchor spacing apart on each side, chosen candidate first, with
            // the normal toward their side. Another candidate's arrow is drawn instead
            // of the chosen one's, so it must also clear painted lines other
            // than its carrier and the leaders drawn so far.
            let mut anchors: Vec<(usize, Point, Point, f64)> = Vec::new();
            for choice in std::iter::once(choices[index])
                .chain((0..label.candidates.len()).filter(|&choice| choice != choices[index]))
            {
                let candidate = &label.candidates[choice];
                let (Some(anchor), Some(normal)) = (candidate.anchor, candidate.normal) else {
                    continue;
                };
                let moved_arrow_clear = |arrow: Bounds| {
                    choice == choices[index]
                        || (grid.near(arrow).into_iter().all(|line| {
                            carrier(line)
                                || arrow.segment_distance(lines[line].start, lines[line].end)
                                    > lines[line].radius
                        }) && drawn
                            .iter()
                            .all(|&[a, b]| arrow.segment_distance(a, b) >= LEADER_CLEARANCE))
                };
                let clear = candidate.arrow_bounds.iter().all(|&arrow| {
                    nodes
                        .iter()
                        .chain(&other_arrows)
                        .chain(others.iter().map(|text| &text.bounds))
                        .all(|&other| arrow.overlap(other, TEXT_ARROW_CLEARANCE) == 0.0)
                        && moved_arrow_clear(arrow)
                }) && nodes
                    .iter()
                    .all(|node| node.segment_distance(anchor, anchor) >= VERTEX_CLEARANCE);
                if clear
                    && anchors.iter().all(|taken| {
                        taken.1.distance(anchor) >= ANCHOR_SPACING || taken.2.dot(normal) <= 0.0
                    })
                {
                    let relative = candidate.relative.map_or(0.0, f64::abs);
                    anchors.push((choice, anchor, normal, ANCHOR_WEIGHT * relative));
                }
            }
            let least = anchors
                .iter()
                .map(|anchor| anchor.3)
                .fold(f64::INFINITY, f64::min);
            let mut best: Option<(f64, usize, Point, Point, Point, Point)> = None;
            let rings = ((LEADER_LONGEST + reach) / RING_STEP).ceil() as u32;
            // Grow rings around every anchor together, so the first clear
            // spots bound the rest: a leader is at least its ring minus the
            // text's reach long.
            for ring in 1..=rings {
                let distance = f64::from(ring) * RING_STEP;
                if best.is_some_and(|(cost, ..)| distance - reach + least > cost) {
                    break;
                }
                for &(choice, anchor, normal, preference) in &anchors {
                    if best.is_some_and(|(cost, ..)| distance - reach + preference > cost) {
                        continue;
                    }
                    let own_arrows = &label.candidates[choice].arrow_bounds;
                    for &direction in &directions {
                        let departure = direction.dot(normal);
                        if departure < LEADER_DEPARTURE {
                            continue;
                        }
                        let spot = anchor.add(direction.scaled(distance));
                        let moved = text.shifted(spot.sub(center));
                        if !region.contains(moved.bounds) {
                            continue;
                        }
                        let clearance = fixed.clearance(moved.bounds);
                        if clearance == Clearance::Blocked {
                            continue;
                        }
                        let Some(frame) = TextFrame::new(moved.bounds, moved.corners) else {
                            continue;
                        };
                        // The leader ends a gap before the text, on the ray to its center.
                        let inside = frame
                            .clip(&StrokeLine {
                                start: anchor,
                                end: spot,
                                radius: 0.0,
                            })
                            .unwrap_or(0.0);
                        let length = distance - inside - LEADER_GAP;
                        let bound =
                            length + preference + OUTSIDE_WEIGHT * moved.bounds.excess(extent);
                        if !(LEADER_SHORTEST..=LEADER_LONGEST).contains(&length)
                            || best.is_some_and(|(best, ..)| bound >= best - 1e-12)
                        {
                            continue;
                        }
                        let end = anchor.add(direction.scaled(length));
                        let risen =
                            anchor.add(direction.scaled((LEADER_RISE / departure).min(length)));
                        // Exact tests: what the raster may have missed, the
                        // leaders drawn so far, then the leader itself.
                        let text_clear =
                            own_arrows.iter().all(|&arrow| {
                                moved.bounds.overlap(arrow, TEXT_ARROW_CLEARANCE) == 0.0
                            }) && movable.iter().all(|&(other, clearance)| {
                                moved.bounds.overlap(other, clearance) == 0.0
                            }) && drawn.iter().all(|&[a, b]| {
                                moved.bounds.segment_distance(a, b) >= LEADER_CLEARANCE
                            }) && (clearance == Clearance::Clear
                                || (grid
                                    .near(moved.bounds.expanded(TEXT_LINE_CLEARANCE))
                                    .into_iter()
                                    .all(|line| {
                                        frame
                                            .clip(&StrokeLine {
                                                radius: lines[line].radius + TEXT_LINE_CLEARANCE,
                                                ..lines[line]
                                            })
                                            .is_none()
                                    })
                                    && nodes.iter().all(|&node| {
                                        moved.bounds.overlap(node, VERTEX_CLEARANCE) == 0.0
                                    })
                                    && other_arrows.iter().all(|&arrow| {
                                        moved.bounds.overlap(arrow, TEXT_ARROW_CLEARANCE) == 0.0
                                    })
                                    && others.iter().all(|other| {
                                        moved.bounds.overlap(other.bounds, self.line_height) == 0.0
                                    })));
                        if !text_clear {
                            continue;
                        }
                        let leader_clear = nodes
                            .iter()
                            .all(|node| node.segment_distance(anchor, end) >= VERTEX_CLEARANCE)
                            && other_arrows.iter().all(|arrow| {
                                arrow.segment_distance(anchor, end) >= LEADER_CLEARANCE
                            })
                            && others.iter().all(|other| {
                                other.bounds.segment_distance(anchor, end) >= LEADER_CLEARANCE
                            })
                            && drawn
                                .iter()
                                .all(|&[a, b]| !segments_meet(anchor, end, a, b))
                            && lines_clear(anchor, risen, end);
                        if !leader_clear {
                            continue;
                        }
                        let near = moved.bounds.expanded(1.0);
                        let congestion: f64 =
                            TextFrame::new(near, near.corners()).map_or(0.0, |frame| {
                                grid.near(near)
                                    .into_iter()
                                    .filter_map(|line| frame.clip(&lines[line]))
                                    .sum()
                            });
                        let cost = bound + CONGESTION_WEIGHT * congestion;
                        if best.is_none_or(|(best, ..)| cost < best - 1e-12) {
                            best = Some((cost, choice, anchor, risen, spot, end));
                        }
                    }
                }
            }
            let Some((_, choice, anchor, risen, spot, end)) = best else {
                continue;
            };
            let moved = text.shifted(spot.sub(center));
            texts[index] = moved;
            arrows[index].clone_from(&label.candidates[choice].arrow_bounds);
            drawn.push([anchor, end]);
            leaders[index] = Some(Leader {
                choice,
                position: local(spot),
                from: local(anchor),
                to: local(end),
                bounds: moved.bounds,
                corners: moved.corners,
                line: [anchor, risen, end],
            });
        }
        leaders
    }
}

impl LabelOverlaps {
    /// What a leader `[anchor, risen, end]` crosses: the painted lines it
    /// comes within their stroke of, its own `carrier`'s only once risen off
    /// it, and the labels whose text it passes through.
    pub(super) fn leader_contacts(
        [anchor, risen, end]: [Point; 3],
        carrier: impl Fn(usize) -> bool,
        lines: &[StrokeLine],
        grid: &LineGrid,
        texts: &[(usize, Bounds)],
    ) -> (Vec<usize>, Vec<usize>) {
        (
            grid.near(Bounds::spanning(anchor, end).expanded(1.0))
                .into_iter()
                .filter(|&index| {
                    let from = if carrier(index) { risen } else { anchor };
                    let line = &lines[index];
                    segment_gap(from, end, line.start, line.end) <= line.radius
                })
                .collect(),
            texts
                .iter()
                .filter(|(_, text)| text.segment_distance(anchor, end) <= 0.0)
                .map(|&(other, _)| other)
                .collect(),
        )
    }
}
