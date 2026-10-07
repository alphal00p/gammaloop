//! Stage 1: untyped topologies by McKay's canonical augmentation, one vertex at a time.
//!
//! An object is a connected partial graph: vertices of the signature degrees, each carrying a
//! multiset of external legs, self-loops, edges to other vertices and open stubs still to be
//! paired. An augmentation adds one vertex (degree, legs, self-loops) and sews some of its stubs to
//! open stubs of existing vertices, given as a multiplicity per existing vertex; its inverse
//! deletes a vertex whose removal keeps the rest connected.
//!
//! Rule 1 keeps one augmentation per orbit of the parent's node automorphism group: the
//! multiplicity vector must be lexicographically minimal among its images. Rule 2 accepts a child
//! only if the added vertex is in the automorphism orbit of the child's canonical deletion: the
//! deletable vertex with the largest cheap invariant (open stubs first), ties broken by the largest
//! canonical position. The invariants of the parent's vertices determine most of the child's, so
//! rule 2 mostly rejects before the child is built. Every prune is a necessary condition for
//! completing to a valid output, so every ancestor of an output along its canonical deletion path
//! survives and every output is generated exactly once.
//!
//! The search is a depth-first traversal with an explicit stack of levels, so it can stop at every
//! complete topology and resume later: generation streams without buffering subtrees.

use std::sync::atomic::{AtomicBool, Ordering};

use ahash::AHashSet;

use super::{model::Model, AbortCheck};
use crate::half_edge::canon::{mix, NodeCanon, Orbits};

/// Largest number of vertices: neighbour sets are `u64` bit sets.
pub(crate) const MAX_VERTICES: usize = 64;

/// Bounds for one loop order.
#[derive(Debug)]
pub(crate) struct Limits {
    max_vertices: usize,
    max_bridges: Option<usize>,
    zero_flow_edges: bool,
    self_loops: bool,
    /// `sum (degree - 2)` over all vertices.
    budget: i32,
    /// `most[k * width + (r - low)]`: the most vertices, at most `k`, whose `sum (degree - 2)` can
    /// be `r` (`-1` if none).
    most: Vec<i16>,
    low: i32,
    width: usize,
}

impl Limits {
    pub(crate) fn new(
        model: &Model,
        budget: i32,
        max_vertices: usize,
        max_bridges: Option<usize>,
        zero_flow_edges: bool,
        self_loops: bool,
    ) -> Self {
        let excess = |d: &u8| *d as i32 - 2;
        let lo = model.degrees.iter().map(excess).min().unwrap_or(0).min(0);
        let hi = model.degrees.iter().map(excess).max().unwrap_or(0).max(0);
        let low = lo * max_vertices as i32;
        let width = ((hi - lo) as usize) * max_vertices + 1;
        let mut reach = vec![false; width];
        reach[(-low) as usize] = true;
        let mut most = vec![-1i16; (max_vertices + 1) * width];
        for m in 0..=max_vertices {
            for i in 0..width {
                let previous = if m > 0 { most[(m - 1) * width + i] } else { -1 };
                most[m * width + i] = if reach[i] { m as i16 } else { previous };
            }
            let mut next = vec![false; width];
            for i in (0..width).filter(|&i| reach[i]) {
                for d in &model.degrees {
                    let j = i as i32 + excess(d);
                    if (0..width as i32).contains(&j) {
                        next[j as usize] = true;
                    }
                }
            }
            reach = next;
        }
        Limits {
            max_vertices,
            max_bridges,
            zero_flow_edges,
            self_loops,
            budget,
            most,
            low,
            width,
        }
    }

    /// Some further vertices can bring `n` vertices with `open` open stubs, `legs` unplaced legs
    /// and remaining budget `r` to a connected complete graph (necessary conditions). `m` new
    /// vertices bring `r + 2m` stubs, which must absorb the open stubs and the legs and pair up
    /// among themselves, so the largest feasible `m` decides. Whatever `m`, `r + open - legs` is
    /// twice the loops still to close: the partial graph cannot have more loops than its
    /// connected completions.
    fn completable(&self, n: usize, r: i32, open: u32, legs: u32) -> bool {
        if n > self.max_vertices {
            return false;
        }
        if open == 0 {
            return r == 0 && legs == 0;
        }
        let i = r - self.low;
        let unclosed = r + open as i32 - legs as i32;
        if !(0..self.width as i32).contains(&i) || unclosed < 0 || unclosed % 2 != 0 {
            return false;
        }
        let k = self.max_vertices.saturating_sub(n);
        let m = self.most[k * self.width + i as usize] as i32;
        m >= 1 && r + 2 * m >= (open + legs) as i32
    }
}

/// Packed vertex invariant: open stubs, degree, self-loops, legs, distinct neighbours (fewer is
/// larger) and a hash of the neighbour degrees and multiplicities.
fn invariant(open: u8, degree: u8, loops: u8, legs: u8, neighbours: u32, weight: u64) -> u64 {
    (open as u64) << 56
        | (degree as u64) << 48
        | (loops as u64) << 40
        | (legs as u64) << 32
        | ((64 - neighbours) as u64) << 24
        | (weight & WEIGHT)
}

const WEIGHT: u64 = 0xFF_FFFF;

/// Contribution of `multiplicity` edges to a neighbour of `degree` to the invariant weight.
fn edge_weight(multiplicity: u8, degree: u8) -> u64 {
    (multiplicity as u64) << (4 * (degree as u64 % 6))
}

/// A connected partial graph.
#[derive(Clone, Debug)]
pub(crate) struct Topology {
    pub(crate) n: usize,
    cap: usize,
    pub(crate) vertices: Vec<Vertex>,
    /// Multiplicities between distinct vertices, `cap x cap`.
    adj: Vec<u8>,
    unplaced: u64,
    unplaced_count: u32,
    /// Remaining `sum (degree - 2)`.
    budget: i32,
    open_total: u32,
}

/// One vertex of a [`Topology`].
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Vertex {
    pub(crate) degree: u8,
    pub(crate) open: u8,
    pub(crate) loops: u8,
    n_legs: u8,
    /// Packed leg multiset.
    pub(crate) legs: u64,
    /// Sum of the edge weights of the neighbours, for the invariant.
    weight: u64,
    nbr: u64,
}

/// The vertex an augmentation adds.
#[derive(Clone, Copy, Debug, Default)]
struct Corolla {
    degree: u8,
    legs: u64,
    n_legs: u32,
    loops: u8,
    /// Stubs not used by legs or self-loops.
    free: u8,
    /// Remaining `sum (degree - 2)` and unplaced legs once it is added.
    budget: i32,
    unplaced: u32,
}

/// Neighbours of `u` with a larger index.
fn above(u: usize) -> u64 {
    u64::MAX.checked_shl(u as u32 + 1).unwrap_or(0)
}

impl Topology {
    pub(crate) fn empty(model: &Model, limits: &Limits) -> Self {
        let cap = limits.max_vertices.max(1);
        Topology {
            n: 0,
            cap,
            vertices: vec![Vertex::default(); cap],
            adj: vec![0; cap * cap],
            unplaced: model.legs.all,
            unplaced_count: model.legs.total,
            budget: limits.budget,
            open_total: 0,
        }
    }

    pub(crate) fn multiplicity(&self, u: usize, v: usize) -> u8 {
        self.adj[u * self.cap + v]
    }

    pub(crate) fn is_complete(&self) -> bool {
        self.n > 0 && self.budget == 0 && self.open_total == 0 && self.unplaced_count == 0
    }

    /// Add `corolla` as vertex `n`, sewing `counts[u]` of its stubs to vertex `u`.
    fn push(&mut self, corolla: &Corolla, counts: &[u8]) {
        let w = self.n;
        let sewn: u8 = counts[..w].iter().sum();
        self.vertices[w].degree = corolla.degree;
        self.vertices[w].loops = corolla.loops;
        self.vertices[w].legs = corolla.legs;
        self.vertices[w].n_legs = corolla.n_legs as u8;
        self.vertices[w].open = corolla.free - sewn;
        self.vertices[w].nbr = 0;
        self.vertices[w].weight = 0;
        for (u, &c) in counts[..w].iter().enumerate() {
            if c > 0 {
                self.adj[u * self.cap + w] = c;
                self.adj[w * self.cap + u] = c;
                self.vertices[u].open -= c;
                self.vertices[u].weight += edge_weight(c, corolla.degree);
                self.vertices[w].weight += edge_weight(c, self.vertices[u].degree);
                self.vertices[u].nbr |= 1 << w;
                self.vertices[w].nbr |= 1 << u;
            }
        }
        self.n += 1;
        self.unplaced -= corolla.legs;
        self.unplaced_count -= corolla.n_legs;
        self.budget -= corolla.degree as i32 - 2;
        self.open_total = self.open_total + self.vertices[w].open as u32 - sewn as u32;
    }

    fn pop(&mut self, corolla: &Corolla, counts: &[u8]) {
        self.n -= 1;
        let w = self.n;
        let sewn: u8 = counts[..w].iter().sum();
        self.open_total = self.open_total + sewn as u32 - self.vertices[w].open as u32;
        for (u, &c) in counts[..w].iter().enumerate() {
            if c > 0 {
                self.adj[u * self.cap + w] = 0;
                self.adj[w * self.cap + u] = 0;
                self.vertices[u].open += c;
                self.vertices[u].weight -= edge_weight(c, corolla.degree);
                self.vertices[u].nbr &= !(1 << w);
            }
        }
        self.unplaced += corolla.legs;
        self.unplaced_count += corolla.n_legs;
        self.budget += corolla.degree as i32 - 2;
    }

    /// The vertices other than `x` are connected.
    fn connected_without(&self, x: usize) -> bool {
        let all = (u64::MAX >> (64 - self.n)) & !(1 << x);
        if all == 0 {
            return true;
        }
        let mut seen = 1u64 << all.trailing_zeros();
        let mut frontier = seen;
        while frontier != 0 {
            let mut next = 0;
            let mut f = frontier;
            while f != 0 {
                next |= self.vertices[f.trailing_zeros() as usize].nbr;
                f &= f - 1;
            }
            next &= all & !seen;
            seen |= next;
            frontier = next;
        }
        seen == all
    }

    /// Cheap isomorphism invariant of a vertex, with its leg multiset as the tie breaker; the
    /// canonical deletion maximizes it.
    fn invariant(&self, v: usize) -> (u64, u64) {
        let packed = invariant(
            self.vertices[v].open,
            self.vertices[v].degree,
            self.vertices[v].loops,
            self.vertices[v].n_legs,
            self.vertices[v].nbr.count_ones(),
            self.vertices[v].weight,
        );
        (packed, self.vertices[v].legs)
    }

    /// Second-order invariant: a commutative hash of the neighbours' invariants and multiplicities.
    fn refined(&self, v: usize, invariants: &[(u64, u64)]) -> u64 {
        let mut sum = 0u64;
        let mut f = self.vertices[v].nbr;
        while f != 0 {
            let u = f.trailing_zeros() as usize;
            f &= f - 1;
            let (packed, legs) = invariants[u];
            sum = sum.wrapping_add(mix(packed ^ mix(legs ^ self.multiplicity(u, v) as u64)));
        }
        sum
    }

    /// Whether the bridges allow a valid completion. A bridge with a side without open stubs is
    /// permanent: no later edge can close a cycle through it. Permanent bridges count against the
    /// bridge limit, and one is a zero-flow edge if a side ends up without external legs.
    fn bridges_allow(&self, limits: &Limits, scratch: &mut BridgeScratch) -> bool {
        scratch.reset(self.n);
        let mut permanent = 0;
        let total_open = self.open_total;
        let total_legs: u32 = self.vertices[..self.n]
            .iter()
            .map(|v| v.n_legs as u32)
            .sum();
        let unplaced = self.unplaced_count;
        let mut ok = true;
        self.lowlink(0, usize::MAX, scratch, &mut |sub_open, sub_legs| {
            let (rest_open, rest_legs) = (total_open - sub_open, total_legs - sub_legs);
            if sub_open > 0 && rest_open > 0 {
                return;
            }
            permanent += 1;
            if limits.max_bridges.is_some_and(|max| permanent > max) {
                ok = false;
            }
            if !limits.zero_flow_edges {
                let zero_flow = |closed_legs: u32, other_legs: u32| {
                    closed_legs == 0 || other_legs + unplaced == 0
                };
                if (sub_open == 0 && zero_flow(sub_legs, rest_legs))
                    || (rest_open == 0 && zero_flow(rest_legs, sub_legs))
                {
                    ok = false;
                }
            }
        });
        ok
    }

    /// Depth-first low-link pass from `v`; calls `bridge(open, legs)` with the totals of the far
    /// side of every bridge. Returns the subtree totals.
    fn lowlink(
        &self,
        v: usize,
        parent: usize,
        s: &mut BridgeScratch,
        bridge: &mut dyn FnMut(u32, u32),
    ) -> (u32, u32) {
        s.time += 1;
        s.disc[v] = s.time;
        s.low[v] = s.time;
        let mut totals = (self.vertices[v].open as u32, self.vertices[v].n_legs as u32);
        let mut f = self.vertices[v].nbr;
        while f != 0 {
            let u = f.trailing_zeros() as usize;
            f &= f - 1;
            if s.disc[u] == 0 {
                let sub = self.lowlink(u, v, s, bridge);
                totals.0 += sub.0;
                totals.1 += sub.1;
                s.low[v] = s.low[v].min(s.low[u]);
                // a parallel edge to the parent already lowered `low[u]`
                if s.low[u] > s.disc[v] {
                    bridge(sub.0, sub.1);
                }
            } else if u != parent || self.multiplicity(u, v) > 1 {
                s.low[v] = s.low[v].min(s.disc[u]);
            }
        }
        totals
    }
}

#[derive(Debug, Default)]
struct BridgeScratch {
    disc: Vec<u32>,
    low: Vec<u32>,
    time: u32,
}

impl BridgeScratch {
    fn reset(&mut self, n: usize) {
        self.disc.clear();
        self.disc.resize(n, 0);
        self.low.clear();
        self.low.resize(n, 0);
        self.time = 0;
    }
}

/// An explicit permutation group on the vertices, identity first.
#[derive(Debug, Default)]
pub(crate) struct Group {
    n: usize,
    elements: Vec<u8>,
    seen: AHashSet<Vec<u8>>,
}

impl Group {
    pub(crate) fn len(&self) -> usize {
        self.elements.len().checked_div(self.n).unwrap_or(1)
    }

    pub(crate) fn elements(&self) -> impl Iterator<Item = &[u8]> {
        self.elements.chunks_exact(self.n.max(1))
    }

    /// Close the generators (images of `0..n`, back to back) under composition.
    pub(crate) fn close(&mut self, n: usize, generators: &[u32]) {
        self.n = n;
        self.elements.clear();
        self.elements.extend(0..n as u8);
        if generators.is_empty() || n == 0 {
            return;
        }
        self.seen.clear();
        self.seen.insert(self.elements.clone());
        let mut i = 0;
        let mut image = vec![0u8; n];
        while i * n < self.elements.len() {
            for g in generators.chunks_exact(n) {
                for x in 0..n {
                    image[x] = g[self.elements[i * n + x] as usize] as u8;
                }
                if !self.seen.contains(&image) {
                    self.seen.insert(image.clone());
                    self.elements.extend_from_slice(&image);
                }
            }
            i += 1;
        }
    }
}

/// Per-depth state of the search: the object at this depth and the cursor over its children.
#[derive(Debug, Default)]
struct Level {
    /// Node automorphism generators of the object at this depth (stride `n`), when known.
    generators: Vec<u32>,
    known: bool,
    /// Vertex invariants of the object at this depth, while it is expanded.
    invariants: Vec<(u64, u64)>,
    /// Its deletable vertices by decreasing invariant.
    deletable: Vec<usize>,
    subsets: Vec<(u64, u32)>,
    open_vertices: OpenStubs,
    /// The corollas that can be added, in search order, and the position of the next one.
    corollas: Vec<Corolla>,
    next: usize,
    /// The current augmentation: `corolla` sewn with multiplicities `counts`, `sewn` in total.
    corolla: Corolla,
    counts: Vec<u8>,
    sewn: u8,
    /// The orbits of the sewings of size `sewn` met so far.
    sewings: SewingOrbits,
    /// Whether `sewn` and `counts` enumerate the sewings of `corolla`.
    sewing: bool,
    /// Whether the current augmentation is applied to the topology.
    pushed: bool,
}

/// The open vertices of an object with their numbers of open stubs. Multiplicity vectors with a
/// fixed sum are enumerated over them in lexicographic order.
#[derive(Debug, Default)]
struct OpenStubs(Vec<(usize, u8)>);

impl OpenStubs {
    /// Write the smallest vector placing `rest` stubs on the vertices from position `from` on;
    /// `false` if they have fewer open stubs.
    fn fill(&self, counts: &mut [u8], from: usize, mut rest: u32) -> bool {
        let mut capacity: u32 = self.0[from..].iter().map(|&(_, b)| b as u32).sum();
        if rest > capacity {
            return false;
        }
        for &(u, bound) in &self.0[from..] {
            capacity -= bound as u32;
            let c = rest.saturating_sub(capacity);
            counts[u] = c as u8;
            rest -= c;
        }
        true
    }

    /// Advance to the next vector with the same sum; `false` after the last.
    fn next(&self, counts: &mut [u8]) -> bool {
        let Some(&(last, _)) = self.0.last() else {
            return false;
        };
        let mut tail = counts[last] as u32;
        for i in (0..self.0.len() - 1).rev() {
            let (u, bound) = self.0[i];
            if counts[u] < bound && tail > 0 {
                counts[u] += 1;
                return self.fill(counts, i + 1, tail - 1);
            }
            tail += counts[u] as u32;
        }
        false
    }
}

/// Rule 1 without listing the parent's automorphism group, whose order grows factorially with
/// interchangeable vertices. The sewings of one size are met in lexicographic order and every
/// check before rule 1 is invariant under the parent's automorphisms, so the first sewing met of
/// an orbit is its smallest: marking its orbit, closed under the generators, rejects the others.
#[derive(Debug, Default)]
struct SewingOrbits {
    /// The sewings met so far, each packed as the sorted multiset of the vertices its stubs go to.
    marked: AHashSet<u128>,
    stack: Vec<u128>,
}

// a corolla sews at most `MAX_DEGREE` stubs, six bits per vertex
const _: () = assert!(6 * super::model::MAX_DEGREE <= 128 && MAX_VERTICES <= 1 << 6);

impl SewingOrbits {
    fn clear(&mut self) {
        self.marked.clear();
    }

    fn pack(stubs: &mut [u8]) -> u128 {
        stubs.sort_unstable();
        stubs.iter().fold(0, |key, &v| key << 6 | v as u128)
    }

    /// Whether `counts` is the first sewing met of its orbit under the group generated by
    /// `generators` (stride `counts.len()`); marks that orbit.
    fn first(&mut self, counts: &[u8], generators: &[u32]) -> bool {
        let n = counts.len();
        // a fixed sewing is alone in its orbit
        if generators
            .chunks_exact(n)
            .all(|g| (0..n).all(|u| counts[g[u] as usize] == counts[u]))
        {
            return true;
        }
        let mut stubs = [0u8; super::model::MAX_DEGREE];
        let mut sewn = 0;
        for (u, &c) in counts.iter().enumerate() {
            stubs[sewn..sewn + c as usize].fill(u as u8);
            sewn += c as usize;
        }
        let key = Self::pack(&mut stubs[..sewn]);
        if !self.marked.insert(key) {
            return false;
        }
        self.stack.push(key);
        while let Some(key) = self.stack.pop() {
            for g in generators.chunks_exact(n) {
                for (i, stub) in stubs[..sewn].iter_mut().enumerate() {
                    *stub = g[(key >> (6 * i) & 63) as usize] as u8;
                }
                let image = Self::pack(&mut stubs[..sewn]);
                if self.marked.insert(image) {
                    self.stack.push(image);
                }
            }
        }
        true
    }
}

/// Reusable buffers of the search.
#[derive(Debug, Default)]
pub(crate) struct Buffers {
    levels: Vec<Level>,
    canon: NodeCanon,
    colours: Vec<u32>,
    /// Invariants of the current object.
    invariants: Vec<(u64, u64)>,
    keys: Vec<((u64, u64, u64), u32)>,
    edges: Vec<[u32; 3]>,
    orbits: Orbits,
    bridges: BridgeScratch,
}

impl Buffers {
    fn level(&mut self, depth: usize) -> &mut Level {
        if self.levels.len() <= depth {
            self.levels.resize_with(depth + 1, Level::default);
        }
        &mut self.levels[depth]
    }
}

/// What a search step needs besides its own state.
pub(crate) struct Context<'a> {
    pub(crate) model: &'a Model,
    pub(crate) limits: &'a Limits,
    pub(crate) stop: &'a AtomicBool,
    pub(crate) abort_check: Option<&'a AbortCheck>,
}

impl Context<'_> {
    /// Whether the generation stops; polls the abort check.
    fn stopped(&self) -> bool {
        if self.stop.load(Ordering::Relaxed) {
            return true;
        }
        if self.abort_check.is_some_and(|check| check()) {
            self.stop.store(true, Ordering::Relaxed);
            return true;
        }
        false
    }
}

/// Where a search stands.
#[derive(Clone, Copy, Debug)]
enum State {
    /// The root object is not examined yet.
    Start,
    /// Enumerating the children of the object at this depth.
    At(usize),
    Done,
}

/// A resumable search of the subtree under one object.
#[derive(Debug)]
pub(crate) struct Search {
    topology: Topology,
    buf: Buffers,
    state: State,
    /// Depth of the object last returned by [`Self::next`].
    found: usize,
    /// Whether complete topologies need their automorphism generators.
    need_generators: bool,
    /// Return the children of the root instead of expanding them.
    children_only: bool,
}

impl Search {
    /// Search the subtree under `topology`, reusing `buf`.
    pub(crate) fn new(model: &Model, topology: Topology, buf: Buffers) -> Self {
        Search {
            topology,
            buf,
            state: State::Start,
            found: 0,
            need_generators: !model.single_diagram(),
            children_only: false,
        }
    }

    /// Enumerate the accepted children of `topology` (a complete one is its own only child).
    pub(crate) fn children(model: &Model, topology: Topology, buf: Buffers) -> Self {
        Search {
            children_only: true,
            ..Search::new(model, topology, buf)
        }
    }

    pub(crate) fn into_buffers(self) -> Buffers {
        self.buf
    }

    /// The object last returned by [`Self::next`].
    pub(crate) fn topology(&self) -> &Topology {
        &self.topology
    }

    /// Node automorphism generators of the complete topology last returned by [`Self::next`]
    /// (stride `n`), unless the model allows a single diagram per topology.
    pub(crate) fn generators(&self) -> &[u32] {
        if self.need_generators {
            &self.buf.levels[self.found].generators
        } else {
            &[]
        }
    }

    /// Advance to the next complete topology, or to the next child of the root with
    /// [`Self::children`]; `false` when the subtree is exhausted or the generation stops.
    pub(crate) fn next(&mut self, cx: &Context<'_>) -> bool {
        loop {
            match self.state {
                State::Done => return false,
                State::Start => {
                    self.buf.level(0).known = false;
                    if self.topology.is_complete() {
                        self.state = State::Done;
                        self.found(0);
                        return true;
                    }
                    if cx.stopped() {
                        self.state = State::Done;
                        return false;
                    }
                    self.prepare(cx, 0);
                    self.state = State::At(0);
                }
                State::At(depth) => {
                    let level = &mut self.buf.levels[depth];
                    if level.pushed {
                        level.pushed = false;
                        self.topology.pop(&level.corolla, &level.counts);
                    }
                    if cx.stop.load(Ordering::Relaxed) {
                        self.state = State::Done;
                        return false;
                    }
                    if !self.next_augmentation(cx, depth) {
                        self.state = match depth {
                            0 => State::Done,
                            _ => State::At(depth - 1),
                        };
                        continue;
                    }
                    let level = &mut self.buf.levels[depth];
                    self.topology.push(&level.corolla, &level.counts);
                    level.pushed = true;
                    let child = depth + 1;
                    self.buf.level(child).known = false;
                    if !self.accepted(cx.limits, depth) {
                        continue;
                    }
                    if self.topology.is_complete() || self.children_only {
                        self.found(child);
                        return true;
                    }
                    if cx.stopped() {
                        self.state = State::Done;
                        return false;
                    }
                    self.prepare(cx, child);
                    self.state = State::At(child);
                }
            }
        }
    }

    /// Record the object at `depth` as found, with its generators if it is a complete topology.
    fn found(&mut self, depth: usize) {
        self.found = depth;
        if self.need_generators && !self.children_only && !self.buf.levels[depth].known {
            self.canonize(depth);
        }
    }

    /// Canonize the object (node colours from the vertex invariants, edges coloured by their
    /// multiplicity) and store its node automorphism generators at `depth`. Distinct invariants
    /// already make the object rigid.
    fn canonize(&mut self, depth: usize) {
        let t = &self.topology;
        let n = t.n;
        let buf = &mut self.buf;
        buf.invariants.clear();
        buf.invariants.extend((0..n).map(|v| t.invariant(v)));
        buf.keys.clear();
        buf.keys.extend((0..n).map(|v| {
            let (packed, legs) = buf.invariants[v];
            ((packed, legs, t.refined(v, &buf.invariants)), v as u32)
        }));
        buf.keys.sort_unstable();
        buf.colours.clear();
        buf.colours.resize(n, 0);
        let mut rank = 0;
        for i in 0..n {
            if i > 0 && buf.keys[i].0 != buf.keys[i - 1].0 {
                rank += 1;
            }
            buf.colours[buf.keys[i].1 as usize] = rank;
        }
        let rigid = rank as usize + 1 == n;
        if !rigid {
            buf.edges.clear();
            for u in 0..n {
                let mut f = t.vertices[u].nbr & above(u);
                while f != 0 {
                    let v = f.trailing_zeros() as usize;
                    f &= f - 1;
                    buf.edges
                        .push([u as u32, v as u32, t.multiplicity(u, v) as u32]);
                }
            }
            buf.canon.run(&buf.colours, &buf.edges);
        }
        let Buffers { canon, levels, .. } = buf;
        let level = &mut levels[depth];
        level.generators.clear();
        if !rigid {
            level.generators.extend(canon.generators().flatten());
        }
        level.known = true;
    }

    /// Set up the enumeration of the children of the object at `depth`: the open vertices, the
    /// invariants for rule 2 and the corollas, in search order.
    fn prepare(&mut self, cx: &Context<'_>, depth: usize) {
        let (model, limits) = (cx.model, cx.limits);
        let t = &self.topology;
        let n = t.n;
        let level = self.buf.level(depth);
        level.counts.clear();
        level.counts.resize(n + 1, 0);
        level.open_vertices.0.clear();
        level.open_vertices.0.extend(
            (0..n)
                .filter(|&u| t.vertices[u].open > 0)
                .map(|u| (u, t.vertices[u].open)),
        );
        level.invariants.clear();
        level.invariants.extend((0..n).map(|v| t.invariant(v)));
        level.deletable.clear();
        level
            .deletable
            .extend((0..n).filter(|&x| n > 1 && t.connected_without(x)));
        let invariants = &level.invariants;
        level
            .deletable
            .sort_unstable_by(|&a, &b| invariants[b].cmp(&invariants[a]));
        // a vertex holds at most its degree in legs
        let max_degree = model.degrees.last().copied().unwrap_or(0);
        model
            .legs
            .subsets(t.unplaced, max_degree as u32, &mut level.subsets);
        level.corollas.clear();
        for &degree in &model.degrees {
            let budget = t.budget - (degree as i32 - 2);
            for &(legs, n_legs) in &level.subsets {
                if n_legs > degree as u32 || !model.legs_fit(degree, legs) {
                    continue;
                }
                let max_loops = if limits.self_loops {
                    (degree as u32 - n_legs) / 2
                } else {
                    0
                };
                let unplaced = t.unplaced_count - n_legs;
                for loops in 0..=max_loops {
                    let free = (degree as u32 - n_legs - 2 * loops) as u8;
                    if n > 0 && free == 0
                        || n == 0 && !limits.completable(1, budget, free as u32, unplaced)
                    {
                        continue;
                    }
                    level.corollas.push(Corolla {
                        degree,
                        legs,
                        n_legs,
                        loops: loops as u8,
                        free,
                        budget,
                        unplaced,
                    });
                }
            }
        }
        level.next = 0;
        level.sewing = false;
        level.pushed = false;
    }

    /// Move the cursor of the level at `depth` to its next augmentation that passes rule 1 and
    /// the checks before the push; `false` after the last.
    fn next_augmentation(&mut self, cx: &Context<'_>, depth: usize) -> bool {
        loop {
            if self.buf.levels[depth].sewing {
                if self.next_sewing(cx, depth) {
                    return true;
                }
                self.buf.levels[depth].sewing = false;
            }
            let level = &mut self.buf.levels[depth];
            let Some(&corolla) = level.corollas.get(level.next) else {
                return false;
            };
            level.next += 1;
            level.corolla = corolla;
            if self.topology.n == 0 {
                // the first vertex: nothing to sew, completability checked in `prepare`
                return true;
            }
            level.sewing = true;
            level.sewn = 0;
        }
    }

    /// Advance to the next multiplicity vector sewing the level's corolla to the open stubs that
    /// passes the rule-2 precheck and rule 1, by increasing number of sewn stubs.
    fn next_sewing(&mut self, cx: &Context<'_>, depth: usize) -> bool {
        let n = self.topology.n;
        let open_total = self.topology.open_total;
        let level = &mut self.buf.levels[depth];
        let corolla = level.corolla;
        let max_sewn = corolla.free.min(open_total.min(u8::MAX as u32) as u8);
        let opens = std::mem::take(&mut level.open_vertices);
        let mut counts = std::mem::take(&mut level.counts);
        let mut sewn = level.sewn;
        let found = loop {
            if sewn == 0 || !opens.next(&mut counts) {
                // the first vector of the next sum that can still be completed: completability
                // depends only on the number of sewn stubs
                let filled = loop {
                    if sewn == max_sewn {
                        break false;
                    }
                    sewn += 1;
                    let open_after = open_total - sewn as u32 + (corolla.free - sewn) as u32;
                    if cx
                        .limits
                        .completable(n + 1, corolla.budget, open_after, corolla.unplaced)
                        && opens.fill(&mut counts, 0, sewn as u32)
                    {
                        break true;
                    }
                };
                if !filled {
                    break false;
                }
                self.buf.levels[depth].sewings.clear();
            }
            let new_open = corolla.free - sewn;
            let forced = corolla.budget == 0
                && open_total - sewn as u32 + new_open as u32 == 0
                && corolla.unplaced == 0;
            if !self.beaten(depth, &corolla, &counts, &opens, new_open)
                && (forced || self.minimal(depth, &counts[..n]))
            {
                break true;
            }
        };
        let level = &mut self.buf.levels[depth];
        level.sewn = sewn;
        level.counts = counts;
        level.open_vertices = opens;
        found
    }

    /// Rule 2 before building the child: a vertex of the parent that is certainly deletable in
    /// the child has a larger invariant there than the new vertex.
    fn beaten(
        &self,
        depth: usize,
        corolla: &Corolla,
        counts: &[u8],
        opens: &OpenStubs,
        new_open: u8,
    ) -> bool {
        let t = &self.topology;
        let level = &self.buf.levels[depth];
        let mut support = 0;
        let mut weight = 0;
        for &(u, _) in &opens.0 {
            if counts[u] > 0 {
                support += 1;
                weight += edge_weight(counts[u], t.vertices[u].degree);
            }
        }
        let key = (
            invariant(
                new_open,
                corolla.degree,
                corolla.loops,
                corolla.n_legs as u8,
                support,
                weight,
            ),
            corolla.legs,
        );
        // invariants only decrease from parent to child, so the scan stops at the first vertex
        // that cannot beat the new one even unchanged
        for &x in &level.deletable {
            let (packed, legs) = level.invariants[x];
            if (packed, legs) <= key {
                return false;
            }
            let c = counts[x];
            if c == 0 {
                // untouched and still deletable
                return true;
            }
            if support > 1 {
                // still deletable: the new vertex also hangs on another vertex
                let weight = (packed + edge_weight(c, corolla.degree)) & WEIGHT;
                let child = ((packed & !WEIGHT) - ((c as u64) << 56) - (1 << 24)) | weight;
                if (child, legs) > key {
                    return true;
                }
            }
        }
        false
    }

    /// Rule 1: `counts` is lexicographically minimal among its images under the parent's group,
    /// as the first sewing met of its orbit ([`SewingOrbits`]).
    fn minimal(&mut self, depth: usize, counts: &[u8]) -> bool {
        if !self.buf.levels[depth].known {
            self.canonize(depth);
        }
        let level = &mut self.buf.levels[depth];
        level.sewings.first(counts, &level.generators)
    }

    /// Bridge prunes and rule 2 for the child just pushed onto the object at `parent`.
    fn accepted(&mut self, limits: &Limits, parent: usize) -> bool {
        let counts = std::mem::take(&mut self.buf.levels[parent].counts);
        let accepted = self.accepted_by(limits, parent, &counts);
        self.buf.levels[parent].counts = counts;
        accepted
    }

    fn accepted_by(&mut self, limits: &Limits, parent: usize, counts: &[u8]) -> bool {
        let t = &self.topology;
        let n = t.n;
        if (limits.max_bridges.is_some() || !limits.zero_flow_edges)
            && !t.bridges_allow(limits, &mut self.buf.bridges)
        {
            return false;
        }
        if n == 1 {
            return true;
        }
        let w = n - 1;
        let key = t.invariant(w);
        let mut ties = 0u64;
        for (x, &c) in counts[..w].iter().enumerate() {
            let k = if c == 0 {
                self.buf.levels[parent].invariants[x]
            } else {
                t.invariant(x)
            };
            if k < key || !t.connected_without(x) {
                continue;
            }
            if k > key {
                return false;
            }
            ties |= 1 << x;
        }
        if ties == 0 {
            return true;
        }
        // second-order invariants of the tied vertices
        let buf = &mut self.buf;
        buf.invariants.clear();
        buf.invariants.extend((0..n).map(|x| {
            if x < w && counts[x] == 0 {
                buf.levels[parent].invariants[x]
            } else {
                t.invariant(x)
            }
        }));
        let refined = t.refined(w, &buf.invariants);
        let mut f = ties;
        while f != 0 {
            let x = f.trailing_zeros() as usize;
            f &= f - 1;
            match t.refined(x, &buf.invariants).cmp(&refined) {
                std::cmp::Ordering::Greater => return false,
                std::cmp::Ordering::Less => ties &= !(1 << x),
                std::cmp::Ordering::Equal => {}
            }
        }
        if ties == 0 {
            return true;
        }
        self.canonize(parent + 1);
        let lab = self.buf.canon.lab();
        let star = (0..n)
            .rev()
            .map(|q| lab[q] as usize)
            .find(|&v| v == w || ties >> v & 1 == 1)
            .expect("the new vertex is a candidate");
        if star == w {
            return true;
        }
        let Buffers { orbits, levels, .. } = &mut self.buf;
        orbits.reset(n);
        for g in levels[parent + 1].generators.chunks_exact(n) {
            orbits.join(g);
        }
        orbits.find(star as u32) == orbits.find(w as u32)
    }
}

#[cfg(test)]
mod tests {
    use super::OpenStubs;

    #[test]
    fn compositions_with_a_fixed_sum_are_enumerated_once() {
        let opens = OpenStubs(vec![(0, 2), (2, 1), (3, 3)]);
        for sum in 0..=6 {
            let mut counts = [0u8; 4];
            let mut seen = vec![];
            if opens.fill(&mut counts, 0, sum) {
                loop {
                    seen.push(counts);
                    if !opens.next(&mut counts) {
                        break;
                    }
                }
            }
            let mut expected = vec![];
            for a in 0..=2u8 {
                for b in 0..=1u8 {
                    for c in 0..=3u8 {
                        if (a + b + c) as u32 == sum {
                            expected.push([a, 0, b, c]);
                        }
                    }
                }
            }
            seen.sort();
            assert_eq!(seen, expected, "sum {sum}");
        }
    }
}
