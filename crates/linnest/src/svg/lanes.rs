//! B4l lane topology and continuous carrier geometry.

use std::collections::{BTreeMap, VecDeque};
use std::f64::consts::TAU;

use kurbo::{BezPath, Point, Vec2};

use super::curves::{Frame, ACCURACY};
use super::{curves, Element, LaneRegion, UNIT};

const PITCH: f64 = 3.0 / UNIT;
const WIDTH: f64 = 2.1;

pub(super) struct Carrier {
    pub(super) nodes: [usize; 2],
    pub(super) path: BezPath,
}

pub(super) struct Station {
    pub(super) center: Point,
    pub(super) radius: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct End {
    edge: usize,
    side: usize,
}

struct Turn {
    slot: usize,
    a: End,
    b: End,
}

#[derive(Clone, Copy)]
struct Constraint {
    a: usize,
    b: usize,
    different: bool,
}

fn turns(ends: &mut [(End, f64)], slot: usize) -> impl Iterator<Item = Turn> + '_ {
    ends.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    let n = ends.len();
    let gap = |i: usize| {
        let delta = (ends[(i + 1) % n].1 - ends[i].1).rem_euclid(TAU);
        if delta == 0.0 {
            TAU
        } else {
            delta
        }
    };
    let skip = (0..n).reduce(|a, b| if gap(a) < gap(b) { b } else { a });
    (0..n)
        .filter(move |i| n >= 2 && Some(*i) != skip)
        .map(move |i| Turn {
            slot,
            a: ends[i].0,
            b: ends[(i + 1) % n].0,
        })
}

fn constraints(turns: &[Turn]) -> Vec<Constraint> {
    // Compare orientations in the outward frame at this vertex. Two source
    // ends (or two sink ends) require different flips. A source and sink
    // already point head-to-tail and therefore require equal flips.
    turns
        .iter()
        .filter(|turn| turn.a.edge != turn.b.edge)
        .map(|turn| Constraint {
            a: turn.a.edge,
            b: turn.b.edge,
            different: turn.a.side == turn.b.side,
        })
        .collect()
}

fn crossings(flips: &[bool], constraints: &[Constraint]) -> usize {
    constraints
        .iter()
        .filter(|c| (flips[c.a] ^ flips[c.b]) != c.different)
        .count()
}

fn orient(size: usize, constraints: &[Constraint]) -> Vec<bool> {
    let mut adjacent = vec![Vec::new(); size];
    for &c in constraints {
        adjacent[c.a].push((c.b, c.different));
        adjacent[c.b].push((c.a, c.different));
    }
    let mut colors = vec![None; size];
    let mut flips = vec![false; size];
    for root in 0..size {
        if colors[root].is_some() || adjacent[root].is_empty() {
            continue;
        }
        colors[root] = Some(false);
        let mut queue = VecDeque::from([root]);
        let mut component = Vec::new();
        let mut conflict = false;
        while let Some(a) = queue.pop_front() {
            component.push(a);
            let color = colors[a].unwrap_or_default();
            for &(b, different) in &adjacent[a] {
                let required = color ^ different;
                match colors[b] {
                    Some(color) => conflict |= color != required,
                    None => {
                        colors[b] = Some(required);
                        queue.push_back(b);
                    }
                }
            }
        }
        let invert = component
            .iter()
            .filter(|&&i| colors[i] == Some(true))
            .count()
            > component.len() / 2;
        for &i in &component {
            flips[i] = colors[i].unwrap_or_default() ^ invert;
        }
        if conflict && component.len() <= 16 {
            // Exact weighted max-cut for the small contradictory components
            // seen in Feynman diagrams. Count every turn, including repetitions.
            component.sort_unstable();
            let edges: Vec<_> = constraints
                .iter()
                .filter(|c| component.binary_search(&c.a).is_ok())
                .collect();
            let ranked: Vec<_> = edges
                .iter()
                .map(|c| {
                    (
                        component.binary_search(&c.a).unwrap(),
                        component.binary_search(&c.b).unwrap(),
                        c.different,
                    )
                })
                .collect();
            let best = (0_u32..1 << component.len())
                .min_by_key(|&mask| {
                    let violated = ranked
                        .iter()
                        .filter(|&&(a, b, different)| {
                            (((mask >> a) ^ (mask >> b)) & 1 != 0) != different
                        })
                        .count();
                    (violated, mask.count_ones(), mask)
                })
                .unwrap_or_default();
            for (bit, &i) in component.iter().enumerate() {
                flips[i] = best & (1 << bit) != 0;
            }
        }
    }
    flips
}

#[derive(Clone, Copy)]
struct Endpoint {
    point: Point,
    toward: Vec2,
}

fn point(frame: Frame, side: usize, offset: f64) -> Endpoint {
    let tangent = Vec2::new(frame.tangent[0], frame.tangent[1]);
    let norm = tangent.hypot();
    let tangent = if norm > f64::EPSILON {
        tangent / norm
    } else {
        Vec2::ZERO
    };
    let position = Point::new(frame.point[0], frame.point[1]);
    Endpoint {
        point: position + Vec2::new(-tangent.y, tangent.x) * offset,
        toward: if side == 0 { -tangent } else { tangent },
    }
}

pub(super) fn draw(
    regions: &[LaneRegion],
    carriers: &BTreeMap<usize, Carrier>,
    stations: &[Station],
    edge_count: usize,
) -> Result<Vec<Element>, String> {
    let mut joins = Vec::new();
    let mut lengths = BTreeMap::new();
    let mut angles = BTreeMap::new();
    for (&edge, carrier) in carriers {
        let length = curves::length(&carrier.path);
        if length <= ACCURACY {
            continue;
        }
        lengths.insert(edge, length);
        let sample = (8.0 / UNIT).min(length / 3.0);
        let frames = curves::frames(&carrier.path, vec![sample, length - sample])?;
        for (side, frame) in frames.into_iter().enumerate() {
            let delta = point(frame, side, 0.0).point - stations[carrier.nodes[side]].center;
            angles.insert(End { edge, side }, delta.y.atan2(delta.x));
        }
    }
    for (slot, region) in regions.iter().enumerate() {
        let mut by_node = BTreeMap::<_, Vec<_>>::new();
        for &edge in &region.edges {
            if let Some(carrier) = carriers.get(&edge) {
                for (side, &node) in carrier.nodes.iter().enumerate() {
                    let end = End { edge, side };
                    if let Some(&angle) = angles.get(&end) {
                        by_node.entry(node).or_default().push((end, angle));
                    }
                }
            }
        }
        for ends in by_node.values_mut() {
            joins.extend(turns(ends, slot));
        }
    }
    let constraints = constraints(&joins);
    let flips = orient(edge_count, &constraints);
    let mut result = vec![Element::LaneCrossings(crossings(&flips, &constraints))];
    let cut = |node: usize| stations[node].radius + 2.0 / UNIT + regions.len() as f64 * PITCH / 2.0;
    let draw_path = |slot: usize, path| Element::Lane {
        slot,
        path,
        paint: regions[slot].paint.clone(),
        details: regions[slot].details.clone(),
        width: WIDTH,
    };
    let mut endpoints = BTreeMap::new();
    for (&edge, carrier) in carriers {
        let Some(&length) = lengths.get(&edge) else {
            continue;
        };
        let outsets = carrier
            .nodes
            .map(|node| (cut(node) - stations[node].radius).min(length / 2.0));
        let shared = length - outsets[0] - outsets[1] <= ACCURACY;
        let window = if shared {
            None
        } else {
            Some(curves::trim(&carrier.path, outsets[0], outsets[1])?)
        };
        for (slot, region) in regions.iter().enumerate() {
            if !region.edges.contains(&edge) {
                continue;
            }
            let offset = (slot as f64 - (regions.len() - 1) as f64 / 2.0)
                * PITCH
                * if flips[edge] { -1.0 } else { 1.0 };
            if let Some(window) = &window {
                let lane = if offset == 0.0 {
                    window.clone()
                } else {
                    curves::parallel(window, offset)?
                };
                let frames = curves::frames(&lane, vec![0.0, curves::length(&lane)])?;
                for (side, frame) in frames.into_iter().enumerate() {
                    endpoints.insert((slot, End { edge, side }), point(frame, side, 0.0));
                }
                result.push(draw_path(slot, lane));
            } else {
                let frames = curves::frames(&carrier.path, vec![length / 2.0])?;
                for (side, frame) in (0..2).zip(frames.into_iter().cycle()) {
                    endpoints.insert((slot, End { edge, side }), point(frame, side, offset));
                }
            }
        }
    }
    for turn in joins {
        if let (Some(a), Some(b)) = (
            endpoints.get(&(turn.slot, turn.a)),
            endpoints.get(&(turn.slot, turn.b)),
        ) {
            let handle = cut(carriers[&turn.a.edge].nodes[turn.a.side]) * 0.62;
            let mut path = BezPath::new();
            path.move_to(a.point);
            path.curve_to(
                a.point + a.toward * handle,
                b.point + b.toward * handle,
                b.point,
            );
            result.push(draw_path(turn.slot, path));
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use kurbo::ParamCurve;
    use serde::Deserialize;

    use super::*;

    #[derive(Deserialize)]
    struct Reference {
        case: String,
        edges: usize,
        crossings: usize,
        stations: Vec<Vec<(usize, usize, f64)>>,
    }

    #[test]
    fn lanes_match_round_three_crossing_counts() {
        // Native y-up angles sampled from the archived B4l prototype's carrier
        // curves. Each row is the ends belonging to one layer at one vertex.
        let cases: Vec<Reference> =
            serde_json::from_str(include_str!("fixtures/round3-lane-stations.json")).unwrap();
        for case in cases {
            let mut joins = Vec::new();
            for (slot, station) in case.stations.into_iter().enumerate() {
                let mut ends: Vec<_> = station
                    .into_iter()
                    .map(|(edge, side, angle)| (End { edge, side }, angle))
                    .collect();
                joins.extend(turns(&mut ends, slot));
            }
            let constraints = constraints(&joins);
            let flips = orient(case.edges, &constraints);
            assert_eq!(
                crossings(&flips, &constraints),
                case.crossings,
                "{}",
                case.case
            );
        }
    }

    #[test]
    fn lanes_join_short_parallel_edges_and_self_loops() {
        let mut paths = BTreeMap::new();
        for (edge, bend) in [(0, -0.2), (1, 0.2)] {
            let mut path = BezPath::new();
            path.move_to((0.1, 0.0));
            path.curve_to((0.12, bend), (0.28, bend), (0.3, 0.0));
            paths.insert(
                edge,
                Carrier {
                    nodes: [0, 1],
                    path,
                },
            );
        }
        let mut path = BezPath::new();
        path.move_to((-0.1, 0.0));
        path.curve_to((-0.9, -1.0), (-0.9, 1.0), (-0.1, 0.0));
        paths.insert(
            2,
            Carrier {
                nodes: [0, 0],
                path,
            },
        );
        let stations = [
            Station {
                center: Point::ORIGIN,
                radius: 0.1,
            },
            Station {
                center: Point::new(0.4, 0.0),
                radius: 0.1,
            },
        ];
        let regions: Vec<_> = [vec![0, 1], vec![2]]
            .into_iter()
            .map(|edges| LaneRegion {
                edges: edges.into_iter().collect(),
                paint: String::new(),
                details: Default::default(),
            })
            .collect();
        let layers = draw(&regions, &paths, &stations, 3).unwrap();
        assert!(matches!(layers[0], Element::LaneCrossings(0)));
        for slot in 0..2 {
            let paths: Vec<_> = layers
                .iter()
                .filter_map(|element| match element {
                    Element::Lane {
                        slot: owner, path, ..
                    } if *owner == slot => Some(path),
                    _ => None,
                })
                .collect();
            assert!(paths.len() >= 2);
            let endpoints: Vec<_> = paths
                .iter()
                .flat_map(|path| {
                    let mut segments = path.segments();
                    let start = segments.next().unwrap();
                    let end = segments.last().unwrap_or(start);
                    [start.start(), end.end()]
                })
                .collect();
            // Every branch end meets another path, including the self-loop.
            for (i, p) in endpoints.iter().enumerate() {
                assert!(endpoints
                    .iter()
                    .enumerate()
                    .any(|(j, q)| j != i && p.distance(*q) < 1e-6));
            }
        }
    }
}
