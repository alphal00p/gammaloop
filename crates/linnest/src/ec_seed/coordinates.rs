//! Planar straight-line coordinates for an embedded drawing graph.
use super::embedding::{require, Graph, Id, OrderedMap, Result};
use indexmap::IndexMap;
use std::collections::{btree_map::Entry, BTreeMap, BTreeSet};

pub(crate) type Point = [f64; 2];
pub(crate) type Positions = OrderedMap<Point>;

pub(crate) fn distance(a: Point, b: Point) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1])
}

pub(crate) fn cross(a: Point, b: Point) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

pub(crate) fn subtract(a: Point, b: Point) -> Point {
    [a[0] - b[0], a[1] - b[1]]
}

/// `std::min`: the first argument unless the second is strictly smaller.
pub(crate) fn min(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}

/// `std::max`: the first argument unless the second is strictly larger.
pub(crate) fn max(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

pub(crate) fn point_segment(p: Point, a: Point, b: Point) -> Result<f64> {
    let d = subtract(b, a);
    let squared = d[0] * d[0] + d[1] * d[1];
    require(squared > 0.0, "Zero-length planar segment")?;
    let f = (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / squared).clamp(0.0, 1.0);
    Ok(distance(p, [a[0] + f * d[0], a[1] + f * d[1]]))
}

const EMPTY: i32 = -1;
const DUMMY: i32 = -2;

/// Integer-set iteration follows CPython's stable positive-integer hash table.
/// The reference's canonical ordering selects with set.pop(); reproducing its
/// tie order keeps the initial geometry unchanged across the native transition.
struct IntSet {
    table: Vec<i32>,
    used: usize,
    fill: usize,
    finger: usize,
}

impl IntSet {
    fn new(values: &[i32]) -> Self {
        let mut out = Self {
            table: vec![EMPTY; 8],
            used: 0,
            fill: 0,
            finger: 0,
        };
        for &v in values {
            out.add(v);
        }
        out
    }

    fn slot(&self, value: i32, insertion: bool) -> usize {
        let mask = self.table.len() - 1;
        let mut i = value as usize & mask;
        let mut perturb = value as usize;
        let mut free = self.table.len();
        loop {
            let mut j = i;
            let mut probes = if i + 9 <= mask { 9 } else { 0 };
            loop {
                let entry = self.table[j];
                if entry == value {
                    return j;
                }
                if entry == EMPTY {
                    return if insertion && free < self.table.len() {
                        free
                    } else {
                        j
                    };
                }
                if entry == DUMMY {
                    free = j;
                }
                j += 1;
                if probes == 0 {
                    break;
                }
                probes -= 1;
            }
            perturb >>= 5;
            i = (i * 5 + 1 + perturb) & mask;
        }
    }

    fn resize(&mut self) {
        let mut size = 8;
        while size <= self.used * 4 {
            size <<= 1;
        }
        let old = std::mem::replace(&mut self.table, vec![EMPTY; size]);
        self.fill = self.used;
        for value in old.into_iter().filter(|&v| v >= 0) {
            let i = self.slot(value, true);
            self.table[i] = value;
        }
    }

    fn contains(&self, v: i32) -> bool {
        self.table[self.slot(v, false)] == v
    }

    fn add(&mut self, v: i32) {
        let i = self.slot(v, true);
        if self.table[i] == v {
            return;
        }
        let empty = self.table[i] == EMPTY;
        self.table[i] = v;
        self.used += 1;
        if empty {
            self.fill += 1;
            if self.fill * 5 >= (self.table.len() - 1) * 3 {
                self.resize();
            }
        }
    }

    fn erase(&mut self, v: i32) {
        let i = self.slot(v, false);
        if self.table[i] == v {
            self.table[i] = DUMMY;
            self.used -= 1;
        }
    }

    fn pop(&mut self) -> Result<i32> {
        require(self.used > 0, "No canonical-ordering candidate")?;
        let mask = self.table.len() - 1;
        let mut i = self.finger & mask;
        while self.table[i] < 0 {
            i = (i + 1) & mask;
        }
        let v = self.table[i];
        self.table[i] = DUMMY;
        self.used -= 1;
        self.finger = i + 1;
        Ok(v)
    }

    fn iteration(&self) -> Vec<i32> {
        self.table.iter().copied().filter(|&v| v >= 0).collect()
    }
}

/// Chrobak–Payne triangulation/canonical ordering translated from NetworkX 3.7.
/// Copyright NetworkX Developers, BSD-3-Clause; see
/// `crates/linnest/layout/native/LICENSE.networkx`.
#[derive(Default)]
struct Embedding {
    rotation: BTreeMap<i32, Vec<i32>>,
    nodes: Vec<i32>,
}

/// Contour bookkeeping of the canonical ordering.
#[derive(Default)]
struct Contour {
    ccw: BTreeMap<i32, i32>,
    cw: BTreeMap<i32, i32>,
    marked: BTreeSet<i32>,
    v1: i32,
}

impl Contour {
    fn on(&self, x: i32) -> bool {
        !self.marked.contains(&x) && (self.ccw.contains_key(&x) || x == self.v1)
    }

    fn cw_at(&self, x: i32) -> Result<i32> {
        self.cw
            .get(&x)
            .copied()
            .ok_or_else(|| format!("Missing clockwise contour neighbor of {x}"))
    }

    fn adjacent(&self, x: i32, y: i32) -> Result<bool> {
        Ok(match (self.ccw.get(&x), self.cw.get(&x)) {
            (None, _) => self.cw_at(x)? == y,
            (Some(&ccw), None) => ccw == y,
            (Some(&ccw), Some(&cw)) => ccw == y || cw == y,
        })
    }
}

impl Embedding {
    fn node(&mut self, n: i32) {
        if let Entry::Vacant(slot) = self.rotation.entry(n) {
            slot.insert(Vec::new());
            self.nodes.push(n);
        }
    }

    fn rotation_at(&self, u: i32) -> Result<&Vec<i32>> {
        self.rotation
            .get(&u)
            .ok_or_else(|| format!("Unknown planar vertex {u}"))
    }

    fn add(&mut self, u: i32, v: i32, reference: Option<i32>, before: bool) -> Result<()> {
        self.node(u);
        self.node(v);
        let r = self.rotation.entry(u).or_default();
        if r.is_empty() {
            require(reference.is_none(), "Reference on empty rotation")?;
            r.push(v);
            return Ok(());
        }
        let reference = reference.ok_or("Missing cyclic insertion reference")?;
        let i = r
            .iter()
            .position(|&x| x == reference)
            .ok_or("Missing cyclic neighbor")?;
        r.insert(if before { i } else { i + 1 }, v);
        Ok(())
    }

    fn neighbor(&self, u: i32, v: i32, delta: isize) -> Result<i32> {
        let r = self.rotation_at(u)?;
        let i = r.iter().position(|&x| x == v).ok_or("Missing halfedge")?;
        Ok(r[(i as isize + r.len() as isize + delta) as usize % r.len()])
    }

    fn edge(&self, u: i32, v: i32) -> Result<bool> {
        Ok(self.rotation_at(u)?.contains(&v))
    }

    fn face(&self, mut u: i32, mut v: i32) -> Result<Vec<i32>> {
        let (start, next) = (u, v);
        let mut result = Vec::new();
        loop {
            result.push(u);
            let w = self.neighbor(v, u, -1)?;
            u = v;
            v = w;
            require(
                result.len() <= self.rotation.len() * self.rotation.len() + 2,
                "Face walk did not close",
            )?;
            if u == start && v == next {
                return Ok(result);
            }
        }
    }

    fn from(data: &BTreeMap<i32, Vec<i32>>) -> Result<Self> {
        let mut out = Self::default();
        for (&u, r) in data {
            let mut reference = None;
            for &v in r.iter().rev() {
                out.add(u, v, reference, true)?;
                reference = Some(v);
            }
        }
        for &u in data.keys() {
            out.node(u);
        }
        Ok(out)
    }

    fn biconnect(
        &mut self,
        start: i32,
        outgoing: i32,
        counted: &mut BTreeSet<(i32, i32)>,
    ) -> Result<Vec<i32>> {
        if !counted.insert((start, outgoing)) {
            return Ok(Vec::new());
        }
        let (mut a, mut b) = (start, outgoing);
        let mut c = self.neighbor(b, a, -1)?;
        let mut face = vec![start];
        let mut seen = BTreeSet::from([start]);
        while b != start || c != outgoing {
            require(a != b, "Invalid halfedge")?;
            if seen.contains(&b) {
                self.add(a, c, Some(b), false)?;
                self.add(c, a, Some(b), true)?;
                counted.insert((b, c));
                counted.insert((c, a));
                b = a;
            } else {
                seen.insert(b);
                face.push(b);
            }
            a = b;
            b = c;
            c = self.neighbor(b, a, -1)?;
            counted.insert((a, b));
        }
        Ok(face)
    }

    fn triangulate_face(&mut self, mut a: i32, mut b: i32) -> Result<()> {
        let mut c = self.neighbor(b, a, -1)?;
        let mut d = self.neighbor(c, b, -1)?;
        if a == b || a == c {
            return Ok(());
        }
        while a != d {
            if self.edge(a, c)? {
                a = b;
                b = c;
                c = d;
            } else {
                self.add(a, c, Some(b), false)?;
                self.add(c, a, Some(b), true)?;
                b = c;
                c = d;
            }
            d = self.neighbor(c, b, -1)?;
        }
        Ok(())
    }

    fn triangulate(&mut self, full: bool) -> Result<Vec<i32>> {
        let mut counted = BTreeSet::new();
        let mut faces: Vec<Vec<i32>> = Vec::new();
        let mut outer = 0;
        for u in self.nodes.clone() {
            let Some(&first) = self.rotation_at(u)?.first() else {
                continue;
            };
            let mut v = first;
            loop {
                let face = self.biconnect(u, v, &mut counted)?;
                if !face.is_empty() {
                    faces.push(face);
                    if faces[faces.len() - 1].len() > faces[outer].len() {
                        outer = faces.len() - 1;
                    }
                }
                v = self.neighbor(u, v, 1)?;
                if v == first {
                    break;
                }
            }
        }
        require(!faces.is_empty(), "No planar face")?;
        for (i, face) in faces.iter().enumerate() {
            if i != outer || full {
                self.triangulate_face(face[0], face[1])?;
            }
        }
        let mut result = faces.swap_remove(outer);
        if full {
            result = vec![
                result[0],
                result[1],
                self.neighbor(result[1], result[0], -1)?,
            ];
        }
        Ok(result)
    }

    fn canonical(&self, outer: &[i32]) -> Result<Vec<(i32, Vec<i32>)>> {
        let (v1, v2) = (outer[0], outer[1]);
        let mut chords: BTreeMap<i32, i32> = BTreeMap::new();
        let mut contour = Contour {
            v1,
            ..Contour::default()
        };
        let mut ready = IntSet::new(outer);
        let mut prev = v2;
        for &n in &outer[2..] {
            contour.ccw.insert(prev, n);
            prev = n;
        }
        contour.ccw.insert(prev, v1);
        prev = v1;
        for &n in outer[1..].iter().rev() {
            contour.cw.insert(prev, n);
            prev = n;
        }
        for &v in outer {
            for &w in self.rotation_at(v)? {
                if contour.on(w) && !contour.adjacent(v, w)? {
                    *chords.entry(v).or_default() += 1;
                    ready.erase(v);
                }
            }
        }
        let mut order = vec![(0, Vec::new()); self.nodes.len()];
        order[0] = (v1, Vec::new());
        order[1] = (v2, Vec::new());
        ready.erase(v1);
        ready.erase(v2);
        for k in (2..self.nodes.len()).rev() {
            let v = ready.pop()?;
            contour.marked.insert(v);
            let (mut wp, mut wq) = (-1, -1);
            for &n in self.rotation_at(v)? {
                if contour.marked.contains(&n) {
                    continue;
                }
                if contour.on(n) {
                    if n == v1 {
                        wp = n;
                    } else if n == v2 {
                        wq = n;
                    } else if contour.cw_at(n)? == v {
                        wp = n;
                    } else {
                        wq = n;
                    }
                }
                if wp >= 0 && wq >= 0 {
                    break;
                }
            }
            require(wp >= 0 && wq >= 0, "Canonical boundary not found")?;
            let mut path = vec![wp];
            let mut n = wp;
            while n != wq {
                let next = self.neighbor(v, n, -1)?;
                path.push(next);
                contour.cw.insert(n, next);
                contour.ccw.insert(next, n);
                n = next;
                require(
                    path.len() <= self.nodes.len(),
                    "Canonical contour did not close",
                )?;
            }
            if path.len() == 2 {
                for w in [wp, wq] {
                    let count = chords.entry(w).or_default();
                    *count -= 1;
                    if *count == 0 {
                        ready.add(w);
                    }
                }
            } else {
                let newly = IntSet::new(&path[1..path.len() - 1]);
                for w in newly.iteration() {
                    ready.add(w);
                    for &other in self.rotation_at(w)? {
                        if contour.on(other) && !contour.adjacent(w, other)? {
                            *chords.entry(w).or_default() += 1;
                            ready.erase(w);
                            if !newly.contains(other) {
                                *chords.entry(other).or_default() += 1;
                                ready.erase(other);
                            }
                        }
                    }
                }
            }
            order[k] = (v, path);
        }
        Ok(order)
    }

    fn coordinates(
        &mut self,
        outer_edge: Option<[i32; 2]>,
        outer_vertex: Option<i32>,
    ) -> Result<IndexMap<i32, Point>> {
        let mut out = IndexMap::new();
        if self.nodes.len() < 4 {
            const TRIANGLE: [Point; 3] = [[0.0, 0.0], [2.0, 0.0], [1.0, 1.0]];
            let mut sequence = self.nodes.clone();
            if let (3, Some([a, b])) = (self.nodes.len(), outer_edge) {
                sequence = vec![a, b];
                sequence.extend(self.nodes.iter().copied().filter(|&n| n != a && n != b));
            }
            for (&n, &p) in sequence.iter().zip(&TRIANGLE) {
                out.insert(n, p);
            }
            return Ok(out);
        }
        let mut outer = self.triangulate(outer_edge.is_some() || outer_vertex.is_some())?;
        if let Some([a, b]) = outer_edge {
            outer = self.face(a, b)?;
        } else if let Some(v) = outer_vertex {
            let first = *self.rotation_at(v)?.first().ok_or("Missing halfedge")?;
            outer = self.face(v, first)?;
        }
        let order = self.canonical(&outer)?;
        let mut left: BTreeMap<i32, Option<i32>> = BTreeMap::new();
        let mut right: BTreeMap<i32, Option<i32>> = BTreeMap::new();
        let mut dx: BTreeMap<i32, i64> = BTreeMap::new();
        let mut y: BTreeMap<i32, i64> = BTreeMap::new();
        let (a, b, c) = (order[0].0, order[1].0, order[2].0);
        dx.insert(a, 0);
        y.insert(a, 0);
        right.insert(a, Some(c));
        left.insert(a, None);
        dx.insert(b, 1);
        y.insert(b, 0);
        right.insert(b, None);
        left.insert(b, None);
        dx.insert(c, 1);
        y.insert(c, 1);
        right.insert(c, Some(b));
        left.insert(c, None);
        for (v, r) in &order[3..] {
            let v = *v;
            let (wp, wp1, wq, wq1) = (r[0], r[1], r[r.len() - 1], r[r.len() - 2]);
            *dx.entry(wp1).or_default() += 1;
            *dx.entry(wq).or_default() += 1;
            let mut gap = 0;
            for n in &r[1..] {
                gap += *dx.entry(*n).or_default();
            }
            let (ywp, ywq) = (*y.entry(wp).or_default(), *y.entry(wq).or_default());
            let dxv = (-ywp + gap + ywq) / 2;
            dx.insert(v, dxv);
            y.insert(v, (ywp + gap + ywq) / 2);
            dx.insert(wq, gap - dxv);
            if r.len() > 2 {
                *dx.entry(wp1).or_default() -= dxv;
            }
            right.insert(wp, Some(v));
            right.insert(v, Some(wq));
            if r.len() > 2 {
                left.insert(v, Some(wp1));
                right.insert(wq1, None);
            } else {
                left.insert(v, None);
            }
        }
        out.insert(a, [0.0, *y.entry(a).or_default() as f64]);
        let mut pending = vec![a];
        while let Some(parent) = pending.pop() {
            let children = [
                *left.entry(parent).or_default(),
                *right.entry(parent).or_default(),
            ];
            for child in children.into_iter().flatten() {
                let x = out.get(&parent).ok_or("Missing parent coordinate")?[0]
                    + *dx.entry(child).or_default() as f64;
                out.insert(child, [x, *y.entry(child).or_default() as f64]);
                pending.push(child);
            }
        }
        Ok(out)
    }
}

impl Graph {
    /// Chrobak–Payne straight-line drawing of every connected component; later
    /// components are scaled into the width of the drawing and stacked above it.
    /// The component of `outer` comes first, with that edge (or vertex) outside.
    pub(crate) fn draw(&self, outer: Option<&[Id]>) -> Result<Positions> {
        let names: Vec<&Id> = self.nodes.iter().collect();
        let ids: BTreeMap<&str, i32> = names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), i as i32))
            .collect();
        let id = |n: &str| {
            ids.get(n)
                .copied()
                .ok_or_else(|| format!("Unknown drawing vertex {n}"))
        };
        let mut rotations: BTreeMap<i32, Vec<i32>> = BTreeMap::new();
        for n in &names {
            let mut r = self.rotation_at(n)?.clone();
            // std::min_element: the first smallest identity starts the rotation.
            if let Some((first, _)) = r.iter().enumerate().min_by_key(|(_, e)| *e) {
                r.rotate_left(first);
            }
            let neighbors = r
                .iter()
                .map(|e| id(self.other(e, n)?))
                .collect::<Result<Vec<i32>>>()?;
            require(
                neighbors.iter().collect::<BTreeSet<_>>().len() == neighbors.len(),
                "Subdivide parallel drawing edges",
            )?;
            rotations.insert(id(n)?, neighbors);
        }
        let full = Embedding::from(&rotations)?;
        let mut seen: BTreeSet<i32> = BTreeSet::new();
        let mut components: Vec<BTreeSet<i32>> = Vec::new();
        for root in 0..names.len() as i32 {
            if seen.contains(&root) {
                continue;
            }
            let mut component = BTreeSet::new();
            let mut pending = vec![root];
            while let Some(n) = pending.pop() {
                if !component.insert(n) {
                    continue;
                }
                pending.extend(&rotations[&n]);
            }
            seen.extend(&component);
            components.push(component);
        }
        let outer_front = outer
            .and_then(|outer| outer.first())
            .map(|front| id(front))
            .transpose()?;
        components.sort_unstable_by_key(|c| {
            (
                outer_front.is_some_and(|front| !c.contains(&front)),
                c.first().copied(),
            )
        });
        let mut positions = Positions::default();
        let mut top = 0.0;
        for component in &components {
            let data = component
                .iter()
                .map(|&n| Ok((n, full.rotation_at(n)?.clone())))
                .collect::<Result<BTreeMap<_, _>>>()?;
            let mut local = Embedding::from(&data)?;
            let mut edge = None;
            let mut vertex = None;
            if let (Some(outer), Some(front)) = (outer, outer_front) {
                if component.contains(&front) {
                    if outer.len() == 2 {
                        edge = Some([id(&outer[0])?, id(&outer[1])?]);
                    } else {
                        vertex = Some(front);
                    }
                }
            }
            let mut coordinates = local.coordinates(edge, vertex)?;
            if !positions.is_empty() {
                let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
                let (mut minimum, mut maximum) = (low, high);
                for p in positions.values() {
                    low = min(low, p[0]);
                    high = max(high, p[0]);
                }
                for p in coordinates.values() {
                    minimum = min(minimum, p[0]);
                    maximum = max(maximum, p[0]);
                }
                let factor = (high - low) / (2.0 * max(1.0, maximum - minimum));
                let center = (minimum + maximum) / 2.0;
                for p in coordinates.values_mut() {
                    *p = [
                        (low + high) / 2.0 + (p[0] - center) * factor,
                        top + 4.0 + p[1] * factor,
                    ];
                }
            }
            for (n, p) in coordinates {
                positions.insert(names[n as usize].clone(), p);
            }
            for p in positions.values() {
                top = max(top, p[1]);
            }
        }
        Ok(positions)
    }
}
