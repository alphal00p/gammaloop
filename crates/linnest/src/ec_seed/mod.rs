//! Embedding-constrained (EC) seed for ImPrEd layouts: a planarized straight-line
//! drawing of a diagram's physical incidence. External legs are grouped by flow
//! on two exterior sides, and crossings of nonplanar diagrams are explicit
//! route points of both crossing edges.
mod coordinates;
mod embedding;
mod spqr;
#[cfg(test)]
mod tests;

use coordinates::{cross, distance, max, min, point_segment, subtract, Point, Positions};
use embedding::{
    as_set, replace, require, ConstraintTree, Ends, Expansion, Graph, Id, Ids, InsertionRecord,
    OrderedMap, Planarization, Result, Set,
};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use spqr::{SplitDecomposer, SpqrDecomposer};
use std::collections::{BTreeMap, BTreeSet};

const EXTERIOR: &str = "__ec_exterior__";
const INCOMING: &str = "incoming";
const OUTGOING: &str = "outgoing";
const RADIAL: &str = "radial";

/// Request of the EC seed plugin protocol.
#[derive(Clone, Debug, Deserialize)]
pub struct SeedRequest {
    pub diagram: Diagram,
    /// Median physical route length of the seed.
    #[serde(default = "SeedRequest::default_scale")]
    pub scale: f64,
    /// Place incoming and outgoing external legs on opposite exterior sides.
    #[serde(default = "SeedRequest::default_external_sides")]
    pub external_sides: bool,
}

/// Physical incidence: vertices are numbered, edges have at most one missing
/// endpoint (external legs) and an optional `incoming`/`outgoing` state.
#[derive(Clone, Debug, Deserialize)]
#[serde(try_from = "Value")]
pub struct Diagram {
    vertices: usize,
    edges: Vec<DiagramEdge>,
}

#[derive(Clone, Debug)]
struct DiagramEdge {
    id: Id,
    source: Value,
    target: Value,
    state: Value,
    /// External legs sharing a row are drawn at one height on opposite sides.
    row: Option<Id>,
}

impl TryFrom<Value> for Diagram {
    type Error = String;

    fn try_from(diagram: Value) -> Result<Self> {
        let vertices = diagram
            .get("vertices")
            .and_then(Value::as_array)
            .ok_or("Diagram vertices must be an array")?
            .len();
        let edges = diagram
            .get("edges")
            .and_then(Value::as_array)
            .ok_or("Diagram edges must be an array")?
            .iter()
            .enumerate()
            .map(|(index, record)| DiagramEdge::from_record(index, record))
            .collect::<Result<_>>()?;
        Ok(Self { vertices, edges })
    }
}

/// Identity of a JSON key or edge ID: strings verbatim, other values as JSON.
fn json_key(value: &Value) -> Id {
    match value {
        Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

impl DiagramEdge {
    /// Read an edge object, or a `[edge, data]` pair whose ID is its position
    /// and whose state comes from `data.external.state`.
    fn from_record(index: usize, record: &Value) -> Result<Self> {
        let mut e = record.clone();
        if let Value::Array(pair) = record {
            let [edge, data] = &pair[..] else {
                return Err("Edge pair must contain edge and data".into());
            };
            e = if edge.is_null() {
                Value::Object(Default::default())
            } else {
                edge.clone()
            };
            let fields = e.as_object_mut().ok_or("Edge must be an object")?;
            fields.insert("id".into(), index.into());
            let data = data.as_object().ok_or("Edge data must be an object")?;
            match data.get("external") {
                None => {
                    fields.insert("state".into(), Value::Null);
                }
                Some(Value::Object(external)) => {
                    fields.insert(
                        "state".into(),
                        external.get("state").cloned().unwrap_or(Value::Null),
                    );
                }
                Some(_) => {}
            }
        }
        let endpoint = |key: &str| {
            e.get(key)
                .cloned()
                .ok_or_else(|| format!("Edge needs a {key}"))
        };
        Ok(Self {
            id: e.get("id").map_or_else(|| index.to_string(), json_key),
            source: endpoint("source")?,
            target: endpoint("target")?,
            state: e.get("state").cloned().unwrap_or(Value::Null),
            row: e.get("row").filter(|row| !row.is_null()).map(json_key),
        })
    }
}

/// External legs of one flow at one vertex, drawn as a fan.
struct Fan {
    flow: &'static str,
    owner: Id,
    spoke: Id,
    endpoints: Ids,
    rows: Vec<Option<Id>>,
}

/// Top-to-bottom order of paired rows, with the spokes it visits down the
/// incoming and the outgoing side.
struct RowOrder {
    rows: Vec<Id>,
    sides: [Ids; 2],
}

impl RowOrder {
    /// Rows in increasing height: the order runs the way its spokes already do
    /// on a side with several, so fans spread consistently on both sides.
    fn rising(self, fans: &[Fan], drawing: &Drawing) -> Result<Vec<Id>> {
        let height = |spoke: &Id| -> Result<f64> {
            let fan = fans
                .iter()
                .find(|f| f.spoke == *spoke)
                .ok_or("Unknown exterior spoke")?;
            Ok(drawing.at(&fan.endpoints[0])?[1])
        };
        let mut rows = self.rows;
        if let Some([first, .., last]) = self
            .sides
            .iter()
            .find(|side| side.len() > 1)
            .map(Vec::as_slice)
        {
            if height(first)? > height(last)? {
                rows.reverse();
            }
        }
        Ok(rows)
    }
}

impl Fan {
    /// Incoming and outgoing spokes of each row, when every leg has a row and
    /// rows pair exactly one incoming with one outgoing leg.
    fn row_spokes(fans: &[Self]) -> Option<BTreeMap<&Id, [&Id; 2]>> {
        let mut rows: BTreeMap<&Id, [Option<&Id>; 2]> = BTreeMap::new();
        for fan in fans {
            for row in &fan.rows {
                let slot =
                    &mut rows.entry(row.as_ref()?).or_default()[usize::from(fan.flow == OUTGOING)];
                if slot.replace(&fan.spoke).is_some() {
                    return None;
                }
            }
        }
        let rows = rows
            .into_iter()
            .map(|(row, [incoming, outgoing])| Some((row, [incoming?, outgoing?])))
            .collect::<Option<BTreeMap<_, _>>>()?;
        (rows.len() >= 2).then_some(rows)
    }

    /// Candidate top-to-bottom row orders, up to the mirror image, that keep the
    /// rows of each fan together on both sides. Larger row sets keep their given
    /// order.
    fn row_orders(spokes: &BTreeMap<&Id, [&Id; 2]>) -> Vec<RowOrder> {
        fn permutations(rest: Vec<usize>, prefix: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
            if rest.is_empty() {
                out.push(prefix.clone());
            }
            for (i, &next) in rest.iter().enumerate() {
                let mut remaining = rest.clone();
                remaining.remove(i);
                prefix.push(next);
                permutations(remaining, prefix, out);
                prefix.pop();
            }
        }
        let rows: Vec<&Id> = spokes.keys().copied().collect();
        let mut orders = Vec::new();
        if rows.len() <= 5 {
            permutations((0..rows.len()).collect(), &mut Vec::new(), &mut orders);
            orders.retain(|order| order[0] < order[order.len() - 1]);
        } else {
            orders.push((0..rows.len()).collect());
        }
        orders
            .into_iter()
            .filter_map(|order| {
                let mut sides: [Ids; 2] = Default::default();
                for (side, visited) in sides.iter_mut().enumerate() {
                    for &i in &order {
                        let spoke = spokes[rows[i]][side];
                        if visited.last() != Some(spoke) {
                            if visited.contains(spoke) {
                                return None;
                            }
                            visited.push(spoke.clone());
                        }
                    }
                }
                Some(RowOrder {
                    rows: order.into_iter().map(|i| rows[i].clone()).collect(),
                    sides,
                })
            })
            .collect()
    }
}

impl ConstraintTree {
    /// Cyclic order of the exterior spokes: flows on two consecutive sides, or
    /// all spokes in one free group. Paired rows run down one side and back up
    /// the other, `i₁…iₖ oₖ…o₁`, so each pair can share a height.
    fn exterior(
        spokes: &BTreeMap<&'static str, Ids>,
        side_mode: bool,
        rows: Option<&RowOrder>,
    ) -> Self {
        let leaves = |ids: &Ids| ids.iter().cloned().map(Self::Leaf).collect::<Vec<_>>();
        match (side_mode, rows) {
            (
                true,
                Some(RowOrder {
                    sides: [incoming, outgoing],
                    ..
                }),
            ) => Self::group(vec![
                Self::oriented(leaves(incoming)),
                Self::oriented(outgoing.iter().rev().cloned().map(Self::Leaf).collect()),
            ]),
            (true, None) => Self::group(
                spokes
                    .values()
                    .filter(|s| !s.is_empty())
                    .map(|s| Self::group(leaves(s)))
                    .collect(),
            ),
            (false, _) => Self::group(spokes.values().flat_map(leaves).collect()),
        }
    }
}

#[derive(Clone, Debug, Serialize)]
struct EdgeEndpoints {
    source: Value,
    target: Value,
}

/// Artificial crossing opened into one route point pair per crossing edge.
#[derive(Clone, Debug, Serialize)]
struct Crossing {
    id: Id,
    edges: Ends,
    point: Point,
    route_points: Vec<Ids>,
    clearance_radius: f64,
}

/// Route segment `(edge, index)`.
type SegmentRef = (Id, usize);

#[derive(Clone, Debug, Serialize)]
struct GeometryAudit {
    valid: bool,
    crossings: Vec<[SegmentRef; 2]>,
    contacts: Vec<[SegmentRef; 2]>,
    overlaps: Vec<[SegmentRef; 2]>,
    degenerate: Vec<SegmentRef>,
    segments: usize,
}

#[derive(Clone, Debug, Serialize)]
struct SideConstraints {
    requested: bool,
    realized: bool,
}

#[derive(Clone, Debug, Serialize)]
struct Report {
    status: &'static str,
    method: &'static str,
    constraint_tree: Option<ConstraintTree>,
    constraint_scope: &'static str,
    external_routes_straight: bool,
    external_side_constraints: SideConstraints,
    insertions: Vec<InsertionRecord>,
    geometry: GeometryAudit,
}

/// Seed positions of vertices and route points, with the physical routes and
/// identities consumed by the ImPrEd projection.
#[derive(Clone, Debug, Serialize)]
pub struct Seed {
    positions: Positions,
    routes: OrderedMap<Ids>,
    node_ids: OrderedMap<Id>,
    external_ids: OrderedMap<Id>,
    incoming_ids: Ids,
    outgoing_ids: Ids,
    edge_endpoints: OrderedMap<EdgeEndpoints>,
    crossings: Vec<Crossing>,
    report: Report,
}

/// Seed drawing: positions of all route points and the physical routes.
struct Drawing {
    positions: Positions,
    routes: OrderedMap<Ids>,
}

impl Drawing {
    fn at(&self, n: &str) -> Result<Point> {
        self.positions.at(n).copied()
    }

    /// Pairwise classification of route segments: proper crossings, contacts
    /// of a segment with a non-adjacent one, collinear overlaps and
    /// zero-length segments.
    fn audit(&self, tolerance: f64) -> Result<GeometryAudit> {
        let segments: Vec<(&Id, usize, &Id, &Id)> = self
            .routes
            .iter()
            .flat_map(|(e, path)| {
                path.windows(2)
                    .enumerate()
                    .map(move |(i, w)| (e, i, &w[0], &w[1]))
            })
            .collect();
        let (mut crossings, mut contacts, mut overlaps, mut degenerate) =
            (vec![], vec![], vec![], vec![]);
        let orient = |a: Point, b: Point, c: Point| cross(subtract(b, a), subtract(c, a));
        let on = |a: Point, b: Point, p: Point| {
            orient(a, b, p).abs() <= tolerance
                && p[0] >= min(a[0], b[0]) - tolerance
                && p[0] <= max(a[0], b[0]) + tolerance
                && p[1] >= min(a[1], b[1]) - tolerance
                && p[1] <= max(a[1], b[1]) + tolerance
        };
        for (i, first) in segments.iter().enumerate() {
            let (a, b) = (self.at(first.2)?, self.at(first.3)?);
            if distance(a, b) <= tolerance {
                degenerate.push((first.0.clone(), first.1));
            }
            for second in &segments[i + 1..] {
                let (c, d) = (self.at(second.2)?, self.at(second.3)?);
                let pair = [(first.0.clone(), first.1), (second.0.clone(), second.1)];
                let values = [
                    orient(a, b, c),
                    orient(a, b, d),
                    orient(c, d, a),
                    orient(c, d, b),
                ];
                if values.iter().all(|v| v.abs() <= tolerance) {
                    let axis = usize::from((a[0] - b[0]).abs() < (a[1] - b[1]).abs());
                    let extent = min(max(a[axis], b[axis]), max(c[axis], d[axis]))
                        - max(min(a[axis], b[axis]), min(c[axis], d[axis]));
                    if extent > tolerance {
                        overlaps.push(pair);
                        continue;
                    }
                }
                let signs = values.map(|v| {
                    if v.abs() <= tolerance {
                        0
                    } else if v > 0.0 {
                        1
                    } else {
                        -1
                    }
                });
                if signs[0] * signs[1] < 0 && signs[2] * signs[3] < 0 {
                    crossings.push(pair);
                } else if first.2 != second.2
                    && first.2 != second.3
                    && first.3 != second.2
                    && first.3 != second.3
                    && (on(a, b, c) || on(a, b, d) || on(c, d, a) || on(c, d, b))
                {
                    contacts.push(pair);
                }
            }
        }
        Ok(GeometryAudit {
            valid: crossings.is_empty()
                && contacts.is_empty()
                && overlaps.is_empty()
                && degenerate.is_empty(),
            crossings,
            contacts,
            overlaps,
            degenerate,
            segments: segments.len(),
        })
    }

    /// Replace every crossing vertex by a pair of route points on each of its
    /// two route spans, inside a clearance disk around the former vertex.
    fn open_crossings(&mut self, nodes: &Set) -> Result<Vec<Crossing>> {
        let original = Self {
            positions: self.positions.clone(),
            routes: self.routes.clone(),
        };
        let segments: BTreeSet<(&Id, &Id)> = original
            .routes
            .values()
            .flat_map(|r| r.windows(2).map(|w| (&w[0], &w[1])))
            .collect();
        let mut replacements: BTreeMap<(&Id, usize), Ids> = BTreeMap::new();
        let mut records = Vec::new();
        for n in nodes {
            let occurrences: Vec<(&Id, usize)> = original
                .routes
                .iter()
                .flat_map(|(e, r)| {
                    r.iter()
                        .enumerate()
                        .filter(|(_, x)| *x == n)
                        .map(move |(i, _)| (e, i))
                })
                .collect();
            require(
                occurrences.len() == 2,
                "Crossing vertex needs two route spans",
            )?;
            let mut neighbors = Ids::new();
            for &(e, i) in &occurrences {
                let r = original.routes.at(e)?;
                require(i > 0 && i + 1 < r.len(), "Crossing at route endpoint")?;
                neighbors.extend([r[i - 1].clone(), r[i + 1].clone()]);
            }
            require(
                as_set(&neighbors).len() == 4,
                "Crossing needs four neighbors",
            )?;
            let center = original.at(n)?;
            let mut radius = f64::INFINITY;
            for (key, &p) in original.positions.iter() {
                if key != n {
                    radius = min(radius, distance(center, p) / 8.0);
                }
            }
            for &(a, b) in &segments {
                if a != n && b != n {
                    radius = min(
                        radius,
                        point_segment(center, original.at(a)?, original.at(b)?)? / 4.0,
                    );
                }
            }
            require(
                radius.is_finite() && radius > 0.0,
                "Crossing lacks clearance",
            )?;
            let mut knots: Vec<Ids> = Vec::new();
            for (occurrence, &(e, i)) in occurrences.iter().enumerate() {
                let r = original.routes.at(e)?;
                let mut pair = Ids::new();
                for (end, neighbor) in [&r[i - 1], &r[i + 1]].into_iter().enumerate() {
                    let direction = subtract(original.at(neighbor)?, center);
                    let length = direction[0].hypot(direction[1]);
                    let key = format!("b:{e}:cross:{n}:{occurrence}:{end}");
                    require(!self.positions.contains(&key), "Crossing knot collision")?;
                    self.positions.insert(
                        key.clone(),
                        [
                            center[0] + radius * direction[0] / length,
                            center[1] + radius * direction[1] / length,
                        ],
                    );
                    pair.push(key);
                }
                replacements.insert((e, i), pair.clone());
                knots.push(pair);
            }
            let (a, b) = (self.at(&knots[0][0])?, self.at(&knots[0][1])?);
            let (c, d) = (self.at(&knots[1][0])?, self.at(&knots[1][1])?);
            let (u, v, offset) = (subtract(b, a), subtract(d, c), subtract(c, a));
            let determinant = cross(u, v);
            require(determinant != 0.0, "Crossing chords are parallel")?;
            let t = cross(offset, v) / determinant;
            let s = cross(offset, u) / determinant;
            require(
                t > 0.0 && t < 1.0 && s > 0.0 && s < 1.0,
                "Routes do not alternate at crossing",
            )?;
            records.push(Crossing {
                id: n.clone(),
                edges: [occurrences[0].0.clone(), occurrences[1].0.clone()],
                point: [a[0] + t * u[0], a[1] + t * u[1]],
                route_points: knots,
                clearance_radius: radius,
            });
        }
        for (e, r) in self.routes.iter_mut() {
            let mut result = Ids::new();
            for (i, x) in r.iter().enumerate() {
                match replacements.get(&(e, i)) {
                    None => result.push(x.clone()),
                    Some(pair) => result.extend(pair.iter().cloned()),
                }
            }
            *r = result;
        }
        for n in nodes {
            self.positions.erase(n);
        }
        Ok(records)
    }

    /// Spread fans of several external legs vertically around their first
    /// endpoint, within a radius clear of every other vertex and segment, with
    /// paired legs in the order of the `rising` rows.
    fn spread_fans(&mut self, groups: &[Fan], rising: &[Id]) -> Result<()> {
        if groups.iter().all(|f| f.endpoints.len() <= 1) {
            return Ok(());
        }
        let segments: Vec<(&Id, &Id)> = self
            .routes
            .values()
            .filter(|r| r.iter().all(|n| self.positions.contains(n)))
            .flat_map(|r| r.windows(2).map(|w| (&w[0], &w[1])))
            .collect();
        let mut radius = f64::INFINITY;
        for &(a, b) in &segments {
            radius = min(radius, distance(self.at(a)?, self.at(b)?) / 4.0);
        }
        for f in groups {
            if f.endpoints.len() == 1 {
                continue;
            }
            let endpoint = &f.endpoints[0];
            let (a, b) = (self.at(&f.owner)?, self.at(endpoint)?);
            let u = subtract(b, a);
            for (key, &p) in self.positions.iter() {
                if *key != f.owner && key != endpoint {
                    radius = min(radius, point_segment(p, a, b)? / 4.0);
                }
            }
            for &(first, second) in &segments {
                let spoke = (&f.owner, endpoint);
                if (first, second) == spoke || (second, first) == spoke {
                    continue;
                }
                let (c, d) = (self.at(first)?, self.at(second)?);
                if f.owner == *first || f.owner == *second {
                    let w = subtract(if f.owner == *first { d } else { c }, a);
                    let det = cross(u, w);
                    let (ul, wl) = (u[0].hypot(u[1]), w[0].hypot(w[1]));
                    if det.abs() <= 1e-12 * ul * wl {
                        require(u[0] * w[0] + u[1] * w[1] < 0.0, "Auxiliary spokes overlap")?;
                        radius = min(min(radius, ul / 4.0), wl / 4.0);
                    } else {
                        radius = min(radius, det.abs() / (4.0 * (u[0].abs() + w[0].abs())));
                    }
                } else {
                    let gaps = [
                        point_segment(a, c, d)?,
                        point_segment(b, c, d)?,
                        point_segment(c, a, b)?,
                        point_segment(d, a, b)?,
                    ];
                    radius = min(
                        radius,
                        gaps[1..].iter().fold(gaps[0], |m, &g| min(m, g)) / 4.0,
                    );
                }
            }
        }
        require(
            radius.is_finite() && radius > 0.0,
            "External fan lacks clearance",
        )?;
        for f in groups.iter().filter(|f| f.endpoints.len() > 1) {
            let p = self.at(&f.endpoints[0])?;
            let last = (f.endpoints.len() - 1) as f64;
            let mut order: Vec<usize> = (0..f.endpoints.len()).collect();
            order.sort_by_key(|&i| {
                let row = f.rows[i].as_ref();
                row.and_then(|row| rising.iter().position(|r| r == row))
            });
            for (i, endpoint) in order.into_iter().map(|k| &f.endpoints[k]).enumerate() {
                self.positions.insert(
                    endpoint.clone(),
                    [p[0], p[1] + radius * (2.0 * i as f64 / last - 1.0)],
                );
            }
        }
        Ok(())
    }
}

impl Positions {
    /// Coordinate-wise mean, correctly rounded. The reference accumulated in
    /// `long double`; its binary128 sum of seed coordinates is exact.
    fn center(&self) -> Point {
        [0, 1].map(|axis| Magnitude::mean(self.values().map(|p| p[axis]), self.len()))
    }
}

/// Unsigned integer in little-endian 64-bit words.
#[derive(Default)]
struct Magnitude(Vec<u64>);

impl Magnitude {
    /// Correctly rounded mean of `count` finite values (zero when empty).
    fn mean(values: impl Iterator<Item = f64>, count: usize) -> f64 {
        // A finite double is m·2^e with an integer m and e >= -1074: sum the
        // magnitudes of each sign exactly, scaled by 2^1074.
        let (mut positive, mut negative) = (Self::default(), Self::default());
        for v in values {
            let bits = v.to_bits();
            let exponent = ((bits >> 52) & 0x7ff) as usize;
            let fraction = bits & ((1 << 52) - 1);
            let (m, shift) = match exponent {
                0 => (fraction, 0),
                _ => (fraction | 1 << 52, exponent - 1),
            };
            let sum = if bits >> 63 == 0 {
                &mut positive
            } else {
                &mut negative
            };
            sum.add_shifted(m, shift);
        }
        let (sign, mut sum) = if positive.is_less(&negative) {
            (-1.0, negative.minus(&positive))
        } else {
            (1.0, positive.minus(&negative))
        };
        if count == 0 || sum.bits() == 0 {
            return 0.0;
        }
        // Quotient with at least 66 significant bits; the remainder is sticky.
        let n = count as u64;
        let k = (66 + 64 - n.leading_zeros() as usize).saturating_sub(sum.bits());
        sum = sum.shifted_left(k);
        let remainder = sum.divide(n);
        let drop = sum.bits() - 53;
        let mut mantissa = (0..53).fold(0u64, |m, i| m | u64::from(sum.bit(drop + i)) << i);
        let sticky = remainder != 0 || (0..drop - 1).any(|i| sum.bit(i));
        if sum.bit(drop - 1) && (sticky || mantissa & 1 == 1) {
            mantissa += 1;
        }
        // mean = mantissa · 2^(drop - k - 1074), scaled in two exact steps.
        let exponent = drop as i32 - k as i32 - 1074;
        sign * (mantissa as f64 * 2f64.powi(exponent / 2) * 2f64.powi(exponent - exponent / 2))
    }

    fn word(&self, i: usize) -> u64 {
        self.0.get(i).copied().unwrap_or_default()
    }

    fn is_less(&self, other: &Self) -> bool {
        let words = self.0.len().max(other.0.len());
        (0..words)
            .rev()
            .map(|i| (self.word(i), other.word(i)))
            .find(|(a, b)| a != b)
            .is_some_and(|(a, b)| a < b)
    }

    fn bit(&self, i: usize) -> bool {
        (self.word(i / 64) >> (i % 64)) & 1 == 1
    }

    fn bits(&self) -> usize {
        self.0
            .iter()
            .rposition(|&w| w != 0)
            .map_or(0, |i| 64 * i + 64 - self.0[i].leading_zeros() as usize)
    }

    /// Add `m · 2^shift`.
    fn add_shifted(&mut self, m: u64, shift: usize) {
        let wide = u128::from(m) << (shift % 64);
        let mut carry = 0;
        for (i, part) in [wide as u64, (wide >> 64) as u64, 0]
            .into_iter()
            .enumerate()
        {
            let at = shift / 64 + i;
            if self.0.len() <= at {
                self.0.resize(at + 1, 0);
            }
            let total = u128::from(self.0[at]) + u128::from(part) + carry;
            self.0[at] = total as u64;
            carry = total >> 64;
        }
        let mut at = shift / 64 + 3;
        while carry != 0 {
            if self.0.len() <= at {
                self.0.push(0);
            }
            let total = u128::from(self.0[at]) + carry;
            self.0[at] = total as u64;
            carry = total >> 64;
            at += 1;
        }
    }

    /// `self - other` for `self >= other`.
    fn minus(&self, other: &Self) -> Self {
        let mut borrow = false;
        Self(
            self.0
                .iter()
                .enumerate()
                .map(|(i, &x)| {
                    let (d, b1) = x.overflowing_sub(other.word(i));
                    let (d, b2) = d.overflowing_sub(u64::from(borrow));
                    borrow = b1 || b2;
                    d
                })
                .collect(),
        )
    }

    fn shifted_left(&self, k: usize) -> Self {
        let (words, bits) = (k / 64, (k % 64) as u32);
        let mut out = vec![0; words];
        let mut carry = 0;
        for &w in &self.0 {
            out.push((w << bits) | carry);
            carry = if bits == 0 { 0 } else { w >> (64 - bits) };
        }
        out.push(carry);
        Self(out)
    }

    /// Divide by `n` in place and return the remainder.
    fn divide(&mut self, n: u64) -> u64 {
        let mut remainder = 0u128;
        for word in self.0.iter_mut().rev() {
            let current = (remainder << 64) | u128::from(*word);
            *word = (current / u128::from(n)) as u64;
            remainder = current % u128::from(n);
        }
        remainder as u64
    }
}

impl Graph {
    /// Contract the wheel an ordered exterior side expanded to into its hub,
    /// whose rotation keeps the side's order.
    fn contract_wheel(&mut self, expansion: &Expansion, port: &Id) -> Result<Id> {
        let members = expansion
            .members
            .get(EXTERIOR)
            .ok_or("Exterior was not expanded")?;
        let mut wheel = Set::from([port.clone()]);
        let mut pending = vec![port.clone()];
        while let Some(n) = pending.pop() {
            for e in self.rotation_at(&n)?.clone() {
                let other = self.other(&e, &n)?.clone();
                if !expansion.tree_edges.contains(&e)
                    && members.contains(&other)
                    && wheel.insert(other.clone())
                {
                    pending.push(other);
                }
            }
        }
        let hub = wheel
            .iter()
            .find(|n| self.wheel_hubs.contains(*n))
            .cloned()
            .ok_or("Ordered exterior side has no wheel hub")?;
        self.contract(&wheel, &hub)?;
        Ok(hub)
    }

    /// Replace the exterior group center by one protected bridge between the
    /// incoming and outgoing side hubs; a side with a single leg gets a metric
    /// hub on its spoke, and an ordered side its contracted wheel. Returns the
    /// hub of each side.
    fn insert_side_hubs(
        &mut self,
        expansion: &Expansion,
        sides: &BTreeMap<&'static str, Ids>,
    ) -> Result<BTreeMap<&'static str, Id>> {
        let ports = expansion
            .leaf_ports
            .get(EXTERIOR)
            .ok_or("Exterior was not expanded")?;
        let port = |spoke: &Id| {
            ports
                .get(spoke)
                .ok_or_else(|| format!("Unknown exterior spoke {spoke}"))
        };
        let mut hubs = BTreeMap::new();
        for side in [INCOMING, OUTGOING] {
            let spokes = &sides[side];
            let mut node = port(&spokes[0])?.clone();
            if let [spoke] = &spokes[..] {
                let owner = self.other(spoke, &node)?.clone();
                let hub = format!("__metric_{side}_hub__");
                let reference = format!("__metric_{side}_reference__");
                require(
                    !self.nodes.contains(&hub) && !self.edges.contains(&reference),
                    "Metric hub identity collision",
                )?;
                self.nodes.insert(hub.clone());
                self.edges.insert(spoke.clone(), [owner, hub.clone()]);
                self.edges
                    .insert(reference.clone(), [hub.clone(), node.clone()]);
                replace(
                    self.rotation.entry(node).or_default(),
                    spoke,
                    std::slice::from_ref(&reference),
                )?;
                self.rotation
                    .insert(hub.clone(), vec![spoke.clone(), reference.clone()]);
                self.protected_edges.insert(reference);
                node = hub;
            } else {
                if spokes.iter().any(|s| port(s).is_ok_and(|p| *p != node)) {
                    node = self.contract_wheel(expansion, &node)?;
                }
                for s in spokes {
                    require(
                        self.edges.at(s)?.contains(&node),
                        "Side did not expand to one group",
                    )?;
                }
            }
            hubs.insert(side, node);
        }
        let hub_ids: Set = hubs.values().cloned().collect();
        let members = expansion
            .members
            .get(EXTERIOR)
            .ok_or("Exterior was not expanded")?;
        let centers: Set = members
            .iter()
            .filter(|n| self.nodes.contains(*n) && !hub_ids.contains(*n))
            .cloned()
            .collect();
        let center = match &centers.into_iter().collect::<Vec<_>>()[..] {
            [center] => center.clone(),
            _ => return Err("Exterior has no unique center".into()),
        };
        let refs = self.rotation_at(&center)?.clone();
        require(
            refs.len() == 2 && refs.iter().all(|e| self.protected_edges.contains(e)),
            "Exterior center is not protected degree two",
        )?;
        let bridge = "__metric_exterior_bridge__".to_owned();
        for e in &refs {
            let other = self.other(e, &center)?.clone();
            replace(
                self.rotation.entry(other).or_default(),
                e,
                std::slice::from_ref(&bridge),
            )?;
            self.edges.erase(e);
            self.protected_edges.remove(e);
        }
        self.rotation.remove(&center);
        self.nodes.remove(&center);
        self.edges.insert(
            bridge.clone(),
            [hubs[INCOMING].clone(), hubs[OUTGOING].clone()],
        );
        self.protected_edges.insert(bridge);
        self.validate()?;
        Ok(hubs)
    }
}

impl SeedRequest {
    fn default_scale() -> f64 {
        2.4
    }

    fn default_external_sides() -> bool {
        true
    }

    /// Embedding-constrained seed: Gutwenger–Klein–Mutzel EC expansion of the
    /// external grouping, optimal individual insertion of edges that cannot be
    /// embedded, Chrobak–Payne coordinates, and straight external legs fanned
    /// out on their side rails. The drawing is centered and scaled to the
    /// requested median route length.
    pub fn initialize(&self) -> std::result::Result<Seed, String> {
        self.initialize_with(&SplitDecomposer)
    }

    fn initialize_with(&self, spqr: &dyn SpqrDecomposer) -> Result<Seed> {
        let Self {
            diagram,
            scale,
            external_sides,
        } = self;
        require(
            scale.is_finite() && *scale > 0.0,
            "Spacing must be finite and positive",
        )?;
        let mut nodes = Set::new();
        let mut node_ids = OrderedMap::default();
        for i in 0..diagram.vertices {
            let n = format!("v:{i}");
            node_ids.insert(i.to_string(), n.clone());
            nodes.insert(n);
        }
        let mut edges = OrderedMap::default();
        let mut routes: OrderedMap<Ids> = OrderedMap::default();
        let mut segments: OrderedMap<Ids> = OrderedMap::default();
        let mut external_ids = OrderedMap::default();
        let mut endpoints = OrderedMap::default();
        let (mut incoming, mut outgoing) = (Ids::new(), Ids::new());
        let mut spokes: BTreeMap<&'static str, Ids> =
            BTreeMap::from([(INCOMING, Ids::new()), (OUTGOING, Ids::new())]);
        let mut fans: Vec<Fan> = Vec::new();
        let mut unique = Set::new();
        for e in &diagram.edges {
            let id = &e.id;
            require(unique.insert(id.clone()), "Duplicate physical edge ID")?;
            let (source, target) = (&e.source, &e.target);
            require(
                !source.is_null() || !target.is_null(),
                "Edge has no endpoint",
            )?;
            let vertex = |v: &Value| {
                if v.is_null() {
                    return Ok(format!("x:{id}"));
                }
                node_ids
                    .get(&json_key(v))
                    .cloned()
                    .ok_or_else(|| "Unknown edge endpoint".to_owned())
            };
            let (a, b) = (vertex(source)?, vertex(target)?);
            endpoints.insert(
                id.clone(),
                EdgeEndpoints {
                    source: source.clone(),
                    target: target.clone(),
                },
            );
            if source.is_null() || target.is_null() {
                let (endpoint, owner) = if source.is_null() {
                    (a.clone(), b.clone())
                } else {
                    (b.clone(), a.clone())
                };
                let flow = match &e.state {
                    Value::Null if source.is_null() => INCOMING,
                    Value::Null => OUTGOING,
                    Value::String(s) if s == INCOMING => INCOMING,
                    Value::String(s) if s == OUTGOING => OUTGOING,
                    _ => return Err("Unknown external flow".into()),
                };
                if flow == INCOMING {
                    &mut incoming
                } else {
                    &mut outgoing
                }
                .push(endpoint.clone());
                external_ids.insert(id.clone(), endpoint.clone());
                routes.insert(id.clone(), vec![a, b]);
                match fans.iter_mut().find(|f| f.flow == flow && f.owner == owner) {
                    Some(f) => {
                        f.endpoints.push(endpoint);
                        f.rows.push(e.row.clone());
                    }
                    None => {
                        let spoke = format!("spoke:{flow}:{owner}");
                        spokes.entry(flow).or_default().push(spoke.clone());
                        edges.insert(spoke.clone(), [owner.clone(), EXTERIOR.to_owned()]);
                        fans.push(Fan {
                            flow,
                            owner,
                            spoke,
                            endpoints: vec![endpoint],
                            rows: vec![e.row.clone()],
                        });
                    }
                }
            } else {
                let bends = if a == b { 2 } else { 1 };
                let mut route = vec![a];
                route.extend((0..bends).map(|i| format!("b:{id}:{i}")));
                route.push(b);
                nodes.extend(route.iter().cloned());
                let mut parts = Ids::new();
                for (i, w) in route.windows(2).enumerate() {
                    let s = format!("route:{id}:{i}");
                    edges.insert(s.clone(), [w[0].clone(), w[1].clone()]);
                    parts.push(s);
                }
                segments.insert(id.clone(), parts);
                routes.insert(id.clone(), route);
            }
        }
        let side_mode = *external_sides && !incoming.is_empty() && !outgoing.is_empty();
        if !fans.is_empty() {
            nodes.insert(EXTERIOR.to_owned());
        }
        let mut initial = Graph::new(nodes, edges)?;
        initial
            .protected_edges
            .extend(fans.iter().map(|f| f.spoke.clone()));
        // Rows pairing the two sides fix one order down both of them; choose the
        // order whose planarization needs the fewest crossings.
        let mut orders: Vec<Option<RowOrder>> = match Fan::row_spokes(&fans) {
            Some(rows) if side_mode => Fan::row_orders(&rows).into_iter().map(Some).collect(),
            _ => Vec::new(),
        };
        if orders.is_empty() {
            orders.push(None);
        }
        let mut best: Option<(
            Expansion,
            Planarization,
            OrderedMap<ConstraintTree>,
            Option<RowOrder>,
        )> = None;
        for rows in orders {
            let mut constraints = OrderedMap::default();
            if !fans.is_empty() {
                constraints.insert(
                    EXTERIOR.to_owned(),
                    ConstraintTree::exterior(&spokes, side_mode, rows.as_ref()),
                );
            }
            let expansion = Expansion::new(&initial, constraints.clone())?;
            let planarization = expansion.graph.planarize(spqr)?;
            let crossings = planarization.crossing_nodes.len();
            if best
                .as_ref()
                .is_none_or(|(_, chosen, ..)| crossings < chosen.crossing_nodes.len())
            {
                best = Some((expansion, planarization, constraints, rows));
                if crossings == 0 {
                    break;
                }
            }
        }
        let (
            expansion,
            Planarization {
                embedding: mut embedded,
                chains,
                crossing_nodes,
                insertions,
            },
            constraints,
            rows,
        ) = best.ok_or("No exterior order")?;
        let mut hubs = BTreeMap::new();
        let mut outer = None;
        if side_mode {
            hubs = embedded.insert_side_hubs(&expansion, &spokes)?;
            outer = Some(vec![hubs[INCOMING].clone(), hubs[OUTGOING].clone()]);
        } else if !fans.is_empty() {
            let hub = expansion
                .members
                .get(EXTERIOR)
                .and_then(|m| m.first())
                .ok_or("Exterior was not expanded")?
                .clone();
            hubs.insert(RADIAL, hub.clone());
            outer = Some(vec![hub]);
        }
        let positions = embedded.draw(outer.as_deref())?;
        for (id, parts) in segments.iter() {
            let route = routes.at_mut(id)?;
            let mut current = route[0].clone();
            let mut planar = vec![current.clone()];
            for segment in parts {
                for part in chains.at(segment)? {
                    current = embedded.other(part, &current)?.clone();
                    planar.push(current.clone());
                }
            }
            require(
                Some(&current) == route.last(),
                "Planarization changed endpoint",
            )?;
            *route = planar;
        }
        let mut drawing = Drawing { positions, routes };
        let mut hub_positions = BTreeMap::new();
        for (&side, hub) in &hubs {
            hub_positions.insert(side, drawing.at(hub)?);
            drawing.positions.erase(hub);
        }
        let mut rails = BTreeMap::new();
        if side_mode {
            let (left, right) = (hub_positions[INCOMING], hub_positions[OUTGOING]);
            let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
            for p in drawing.positions.values() {
                low = min(low, p[0]);
                high = max(high, p[0]);
            }
            require(
                left[0] < low && low <= high && high < right[0],
                "Exterior hub edge did not bound drawing",
            )?;
            rails.insert(INCOMING, (left[0] + low) / 2.0);
            rails.insert(OUTGOING, (right[0] + high) / 2.0);
        }
        for f in &fans {
            let p = drawing.at(&f.owner)?;
            let hub = hub_positions[if side_mode { f.flow } else { RADIAL }];
            let fraction = if side_mode {
                (rails[f.flow] - hub[0]) / (p[0] - hub[0])
            } else {
                0.15
            };
            drawing.positions.insert(
                f.endpoints[0].clone(),
                [
                    hub[0] + fraction * (p[0] - hub[0]),
                    hub[1] + fraction * (p[1] - hub[1]),
                ],
            );
        }
        let rising = match rows {
            Some(rows) => rows.rising(&fans, &drawing)?,
            None => Vec::new(),
        };
        drawing.spread_fans(&fans, &rising)?;
        let mut crossings = drawing.open_crossings(&crossing_nodes.key_set())?;
        let center = drawing.positions.center();
        let mut lengths = drawing
            .routes
            .values()
            .map(|r| {
                let mut length = 0.0;
                for w in r.windows(2) {
                    length += distance(drawing.at(&w[0])?, drawing.at(&w[1])?);
                }
                Ok(length)
            })
            .collect::<Result<Vec<f64>>>()?;
        lengths.sort_unstable_by(f64::total_cmp);
        let mut factor = 1.0;
        if !lengths.is_empty() {
            let m = lengths.len() / 2;
            let median = if lengths.len() % 2 == 1 {
                lengths[m]
            } else {
                (lengths[m - 1] + lengths[m]) / 2.0
            };
            require(median > 0.0, "Zero median route length")?;
            factor = scale / median;
        }
        let rescale = |p: &mut Point| {
            for axis in 0..2 {
                p[axis] = factor * (p[axis] - center[axis]);
            }
        };
        for (_, p) in drawing.positions.iter_mut() {
            rescale(p);
        }
        for c in &mut crossings {
            rescale(&mut c.point);
            c.clearance_radius *= factor;
        }
        let audit = drawing.audit(1e-8)?;
        require(
            audit.contacts.is_empty() && audit.overlaps.is_empty() && audit.degenerate.is_empty(),
            "EC coordinates have geometric degeneracy",
        )?;
        let multiset = |pairs: Vec<Ends>| {
            let mut counts: BTreeMap<Ends, usize> = BTreeMap::new();
            for mut pair in pairs {
                pair.sort_unstable();
                *counts.entry(pair).or_default() += 1;
            }
            counts
        };
        let actual = multiset(
            audit
                .crossings
                .iter()
                .map(|[a, b]| [a.0.clone(), b.0.clone()])
                .collect(),
        );
        let expected = multiset(crossings.iter().map(|c| c.edges.clone()).collect());
        require(
            actual == expected,
            "Geometric crossings disagree with provenance",
        )?;
        Ok(Seed {
            positions: drawing.positions,
            routes: drawing.routes,
            node_ids,
            external_ids,
            incoming_ids: incoming,
            outgoing_ids: outgoing,
            edge_endpoints: endpoints,
            report: Report {
                status: if crossings.is_empty() {
                    "ec-planar"
                } else {
                    "ec-planarized"
                },
                method:
                    "Gutwenger–Klein–Mutzel EC expansion and optimal individual edge insertion; \
                         Chrobak–Payne coordinates",
                constraint_tree: constraints.get(EXTERIOR).cloned(),
                constraint_scope: if side_mode {
                    "consecutive exterior incoming/outgoing groups; metric rails constructed separately"
                } else {
                    "unconstrained radial/cofacial external order"
                },
                external_routes_straight: true,
                external_side_constraints: SideConstraints {
                    requested: *external_sides,
                    realized: side_mode,
                },
                insertions,
                geometry: audit,
            },
            crossings,
        })
    }
}
