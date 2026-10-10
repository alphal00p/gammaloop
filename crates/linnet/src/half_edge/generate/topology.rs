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
//! deletable vertex with the largest cheap invariant (legs first, then open stubs), ties broken by
//! the largest canonical position. The invariants of the parent's vertices determine most of the
//! child's, so rule 2 mostly rejects before the child is built. Every prune is a necessary
//! condition for completing to a valid output, so every ancestor of an output along its canonical
//! deletion path survives and every output is generated exactly once. The generator's filter is
//! called on every accepted child and prunes it too: it is hereditary, so a rejected object has no
//! accepted output below it.
//!
//! Two kinds of steps keep processes with many external legs from searching dead ends:
//!
//! - Pendants: once the counts force every remaining vertex to be a pendant (one sewn stub, no
//!   self-loop) on an open stub, the pendants are hung all at once, one assignment of pendant types
//!   per orbit of the object's group (rule 1), and kept if canonical deletions of the result remove
//!   exactly the pendants (rule 2). Added one at a time, every subset of them would be a canonical
//!   object, most without a canonical completion.
//! - Labels: three or more legs with distinct labels and a common stub are placed as anonymous
//!   legs of one slot class, so the search sees few, symmetric objects. Once the topology is
//!   shaped (every vertex, edge and leg placed), the labels are assigned one at a time in a fixed
//!   order, each to one vertex per orbit of the group fixing the vertices of the previous labels;
//!   the last of these stabilizers is the node group of the labelled topology.
//!
//! Anonymous legs make many twins: vertices with equal corollas, open stubs and multiplicities to
//! every other vertex, such as leaves of one type on one neighbour. Exchanging two twins is an
//! automorphism, so a rule-2 tie among twins needs no canonization, and neither does the group of
//! an object whose colour classes are classes of twins: their transpositions generate it.
//!
//! The search is a depth-first traversal with an explicit stack of levels, so it can stop at every
//! complete topology and resume later: generation streams without buffering subtrees.

use std::sync::atomic::{AtomicBool, Ordering};

use ahash::AHashSet;

use super::{
    model::{Legs, Model},
    view::{Prune, State as View, VertexSet},
    AbortCheck,
};
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

    /// The most vertices that can be added to `n` vertices with remaining budget `r`; negative if
    /// none.
    fn most_new(&self, n: usize, r: i32) -> i32 {
        let i = r - self.low;
        if n > self.max_vertices || !(0..self.width as i32).contains(&i) {
            return -1;
        }
        self.most[(self.max_vertices - n) * self.width + i as usize] as i32
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

/// Packed vertex invariant: legs, open stubs (four bits each: a vertex has at most 15 stubs),
/// degree, self-loops, distinct neighbours (fewer is larger) and a hash of the neighbour degrees
/// and multiplicities.
///
/// Legs come first so that canonical deletions take the vertices with the most legs: the search
/// adds them last, and those that complete the topology as pendants all at once.
fn invariant(open: u8, degree: u8, loops: u8, legs: u8, neighbours: u32, weight: u64) -> u64 {
    (legs as u64) << 60
        | (open as u64) << 56
        | (degree as u64) << 48
        | (loops as u64) << 40
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
    /// Slot legs not attached to a vertex yet, packed by slot class, and their number.
    unplaced: u64,
    unplaced_count: u32,
    /// Labels of anonymous slots not assigned yet.
    unlabelled: u32,
    /// Remaining `sum (degree - 2)`.
    budget: i32,
    open_total: u32,
    /// While checking the filter's heredity: the filter rejected this partial topology or one
    /// on the way to it.
    pub(crate) vetoed: bool,
}

/// One vertex of a [`Topology`].
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Vertex {
    pub(crate) degree: u8,
    pub(crate) open: u8,
    pub(crate) loops: u8,
    pub(crate) n_legs: u8,
    /// Leg multiset, packed by slot class.
    pub(crate) legs: u64,
    /// Its legs with a label, packed by label class, and the anonymous ones among them, packed by
    /// slot class.
    pub(crate) labelled: u64,
    labelled_slots: u64,
    /// Sum of the edge weights of the neighbours, for the invariant.
    weight: u64,
    nbr: u64,
}

impl Vertex {
    /// Its anonymous legs without a label.
    pub(crate) fn unlabelled(&self, model: &Model) -> u32 {
        model.slots.size(model.anonymous(self.legs)) - model.slots.size(self.labelled_slots)
    }
}

/// The vertex an augmentation adds.
#[derive(Clone, Copy, Debug, Default)]
struct Corolla {
    degree: u8,
    legs: u64,
    n_legs: u32,
    /// The legs of its labelled slot classes, packed by label class.
    labelled: u64,
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
            unplaced: model.slots.all,
            unplaced_count: model.slots.total,
            unlabelled: model.labelling.len() as u32,
            budget: limits.budget,
            open_total: 0,
            vetoed: false,
        }
    }

    pub(crate) fn multiplicity(&self, u: usize, v: usize) -> u8 {
        self.adj[u * self.cap + v]
    }

    /// Every vertex, edge and leg is placed and every leg labelled.
    pub(crate) fn is_complete(&self) -> bool {
        self.is_shaped() && self.unlabelled == 0
    }

    /// Every vertex, edge and leg is placed; anonymous legs may still lack their labels.
    fn is_shaped(&self) -> bool {
        self.n > 0 && self.budget == 0 && self.open_total == 0 && self.unplaced_count == 0
    }

    /// Stubs of `v` not taken by an edge or a labelled leg yet: open stubs and anonymous legs
    /// without a label.
    pub(crate) fn free_stubs(&self, model: &Model, v: usize) -> u32 {
        self.vertices[v].open as u32 + self.vertices[v].unlabelled(model)
    }

    /// Index of the next label to assign in [`Model::labelling`].
    fn next_label(&self, model: &Model) -> usize {
        model.labelling.len() - self.unlabelled as usize
    }

    /// Give the anonymous leg of slot class `slot` at `v` the label of class `label`.
    fn assign(&mut self, model: &Model, v: usize, (slot, label): (usize, usize)) {
        let vertex = &mut self.vertices[v];
        vertex.labelled_slots += model.slots.one(slot);
        vertex.labelled += model.legs.one(label);
        self.unlabelled -= 1;
    }

    fn unassign(&mut self, model: &Model, v: usize, (slot, label): (usize, usize)) {
        let vertex = &mut self.vertices[v];
        vertex.labelled_slots -= model.slots.one(slot);
        vertex.labelled -= model.legs.one(label);
        self.unlabelled += 1;
    }

    /// Add `corolla` as vertex `n`, sewing `counts[u]` of its stubs to vertex `u`.
    #[inline(always)]
    fn push(&mut self, corolla: &Corolla, counts: &[u8]) {
        let w = self.n;
        let sewn: u8 = counts[..w].iter().sum();
        self.vertices[w].degree = corolla.degree;
        self.vertices[w].loops = corolla.loops;
        self.vertices[w].legs = corolla.legs;
        self.vertices[w].n_legs = corolla.n_legs as u8;
        self.vertices[w].labelled = corolla.labelled;
        self.vertices[w].labelled_slots = 0;
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

    #[inline(always)]
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

    /// Exchange the numbers of vertices `a` and `b`.
    fn swap(&mut self, a: usize, b: usize) {
        let cap = self.cap;
        self.vertices.swap(a, b);
        for v in &mut self.vertices[..self.n] {
            if (v.nbr >> a ^ v.nbr >> b) & 1 == 1 {
                v.nbr ^= 1 << a | 1 << b;
            }
        }
        for v in 0..self.n {
            self.adj.swap(a * cap + v, b * cap + v);
        }
        for v in 0..self.n {
            self.adj.swap(v * cap + a, v * cap + b);
        }
    }

    /// `a` and `b` are twins: equal corollas, open stubs and labels, and equal multiplicities to
    /// every other vertex, so exchanging them is an automorphism. Twins of one vertex are twins of
    /// each other.
    fn twins(&self, a: usize, b: usize) -> bool {
        let (x, y) = (&self.vertices[a], &self.vertices[b]);
        let others = !(1u64 << a | 1u64 << b);
        (x.degree, x.open, x.loops, x.legs, x.labelled)
            == (y.degree, y.open, y.loops, y.legs, y.labelled)
            && x.nbr & others == y.nbr & others
            && VertexSet(x.nbr & others)
                .iter()
                .all(|v| self.multiplicity(a, v) == self.multiplicity(b, v))
    }

    /// The vertices other than `x` are connected.
    fn connected_without(&self, x: usize) -> bool {
        let all = (u64::MAX >> (64 - self.n)) & !(1 << x);
        all == 0 || self.reach(all.trailing_zeros() as usize, all) == all
    }

    /// The vertices that `from` reaches within the set `within`, which contains it.
    fn reach(&self, from: usize, within: u64) -> u64 {
        let mut seen = 1u64 << from;
        let mut frontier = seen;
        while frontier != 0 {
            let mut next = 0;
            let mut f = frontier;
            while f != 0 {
                next |= self.vertices[f.trailing_zeros() as usize].nbr;
                f &= f - 1;
            }
            next &= within & !seen;
            seen |= next;
            frontier = next;
        }
        seen
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
        for u in VertexSet(self.vertices[v].nbr).iter() {
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
        let stubs = |v: &Vertex| (v.open as u32, v.n_legs as u32);
        self.lowlink(
            0,
            usize::MAX,
            scratch,
            &stubs,
            &mut |_, sub_open, sub_legs| {
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
            },
        );
        ok
    }

    /// The bridges whose `side` (the vertices beyond `ends[1]`) has no free stubs: no later edge
    /// or label reaches that side, so the bridge and the subgraph on that side are final. A
    /// complete topology reports every bridge once per side.
    pub(crate) fn closed_sides(&self, model: &Model, mut side: impl FnMut([usize; 2], u64)) {
        if self.n == 0 {
            return;
        }
        let mut scratch = BridgeScratch::default();
        scratch.reset(self.n);
        let all = u64::MAX >> (64 - self.n);
        let total_open = (0..self.n).map(|v| self.free_stubs(model, v)).sum::<u32>();
        let stubs = |v: &Vertex| (v.open as u32 + v.unlabelled(model), 0);
        self.lowlink(
            0,
            usize::MAX,
            &mut scratch,
            &stubs,
            &mut |[v, u], open, _| {
                if open != 0 && open != total_open {
                    return;
                }
                // without the bridge, `u` reaches its side only
                let far = self.reach(u, all & !(1 << v));
                if open == 0 {
                    side([v, u], far);
                }
                if open == total_open {
                    side([u, v], all & !far);
                }
            },
        );
    }

    /// Calls `bridge(u, v)` for every bridge `(u, v)` of a connected topology.
    pub(crate) fn bridges(&self, mut bridge: impl FnMut(usize, usize)) {
        let mut scratch = BridgeScratch::default();
        scratch.reset(self.n);
        self.lowlink(
            0,
            usize::MAX,
            &mut scratch,
            &|_| (0, 0),
            &mut |[v, u], _, _| bridge(v, u),
        );
    }

    /// Depth-first low-link pass from `v`; calls `bridge([parent, child], open, legs)` with the
    /// totals of `stubs` (open stubs and legs per vertex) over the child's side of every bridge.
    /// Returns the subtree totals.
    fn lowlink(
        &self,
        v: usize,
        parent: usize,
        s: &mut BridgeScratch,
        stubs: &impl Fn(&Vertex) -> (u32, u32),
        bridge: &mut impl FnMut([usize; 2], u32, u32),
    ) -> (u32, u32) {
        s.time += 1;
        s.disc[v] = s.time;
        s.low[v] = s.time;
        let mut totals = stubs(&self.vertices[v]);
        let mut f = self.vertices[v].nbr;
        while f != 0 {
            let u = f.trailing_zeros() as usize;
            f &= f - 1;
            if s.disc[u] == 0 {
                let sub = self.lowlink(u, v, s, stubs, bridge);
                totals.0 += sub.0;
                totals.1 += sub.1;
                s.low[v] = s.low[v].min(s.low[u]);
                // a parallel edge to the parent already lowered `low[u]`
                if s.low[u] > s.disc[v] {
                    bridge([v, u], sub.0, sub.1);
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
    /// How the children of the object at this depth are made.
    mode: Mode,
    /// [`Mode::Label`]: the vertices that can take the label, the current one at `next - 1`.
    targets: Vec<usize>,
    /// [`Mode::Pendants`]: the current assignment of pendants to open stubs.
    pendants: Pendants,
    /// Whether the current augmentation is applied to the topology.
    pushed: bool,
    /// The children found so far, and whether the rest was split off (see [`Search::split`]).
    children: usize,
    done: bool,
    /// While checking the filter's heredity: the filter rejected the object at this depth or one
    /// on the way to it. Kept only with a filter.
    vetoed: bool,
}

impl Level {
    /// Remove the child applied at this level from `t`.
    fn unapply(&mut self, t: &mut Topology, model: &Model) {
        match self.mode {
            Mode::Vertex => t.pop(&self.corolla, &self.counts),
            Mode::Pendants => self.pop_pendants(t),
            Mode::Label => {
                let label = model.labelling[t.next_label(model) - 1];
                t.unassign(model, self.targets[self.next - 1], label);
            }
        }
    }

    /// Hang the pendants of the current assignment, stub by stub.
    fn push_pendants(&mut self, t: &mut Topology) {
        let p = &self.pendants;
        for (&u, &c) in p.owners.iter().zip(&p.choice) {
            self.counts[u] = 1;
            t.push(&p.types[c], &self.counts);
            self.counts[u] = 0;
        }
    }

    fn pop_pendants(&mut self, t: &mut Topology) {
        for i in (0..self.pendants.owners.len()).rev() {
            self.pop_pendant(t, i);
        }
    }

    /// Remove the last vertex of `t`, the pendant of stub `i`.
    fn pop_pendant(&mut self, t: &mut Topology, i: usize) {
        let p = &self.pendants;
        let u = p.owners[i];
        self.counts[u] = 1;
        t.pop(&p.types[p.choice[i]], &self.counts);
        self.counts[u] = 0;
    }
}

/// How the children of an object are made.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum Mode {
    /// Add one vertex, sewn to open stubs.
    #[default]
    Vertex,
    /// Hang a pendant vertex on every open stub at once: every completion of the object consists
    /// of pendants.
    Pendants,
    /// Label an anonymous leg of a shaped object.
    Label,
}

/// Pendant vertices (one sewn stub, no self-loop) hung on the open stubs of an object, one per
/// stub. Assignments of pendant types to stubs are enumerated in lexicographic order,
/// non-decreasing at each vertex, using exactly the unplaced legs.
#[derive(Debug, Default)]
struct Pendants {
    /// The pendant types.
    types: Vec<Corolla>,
    /// The vertex of each open stub, ascending.
    owners: Vec<usize>,
    /// The type at each stub.
    choice: Vec<usize>,
    /// `used[i]`: the legs of the pendants at the stubs before `i`, packed by slot class.
    used: Vec<u64>,
    /// Per vertex of the object: its first stub in `owners`.
    start: Vec<usize>,
    /// The stub of each pendant vertex while rule 2 removes them.
    order: Vec<usize>,
    /// The assignments met so far, as for [`SewingOrbits`].
    marked: AHashSet<Vec<usize>>,
    stack: Vec<Vec<usize>>,
}

impl Pendants {
    /// Advance to the next assignment (the first if `first`) that places exactly `legs`; `false`
    /// after the last.
    fn next(&mut self, slots: &Legs, legs: u64, first: bool) -> bool {
        let k = self.owners.len();
        let mut i = if first {
            self.choice.clear();
            self.choice.resize(k, 0);
            self.marked.clear();
            self.used.clear();
            self.used.resize(k + 1, 0);
            0
        } else {
            self.choice[k - 1] += 1;
            k - 1
        };
        loop {
            if self.choice[i] == self.types.len() {
                if i == 0 {
                    return false;
                }
                i -= 1;
                self.choice[i] += 1;
                continue;
            }
            let add = self.types[self.choice[i]].legs;
            if !slots.fits(self.used[i], add, legs) {
                self.choice[i] += 1;
                continue;
            }
            self.used[i + 1] = self.used[i] + add;
            if i + 1 == k {
                if self.used[k] == legs {
                    return true;
                }
                self.choice[i] += 1;
                continue;
            }
            i += 1;
            // non-decreasing at each vertex: pendants on one vertex are interchangeable
            self.choice[i] = if self.owners[i] == self.owners[i - 1] {
                self.choice[i - 1]
            } else {
                0
            };
        }
    }

    /// The stubs of vertex `u` in `owners`.
    fn block(&self, u: usize) -> std::ops::Range<usize> {
        self.start[u]..self.start.get(u + 1).copied().unwrap_or(self.owners.len())
    }

    /// Rule 1, as for sewings ([`SewingOrbits`]): assignments are met in lexicographic order and
    /// only their legs are checked before, so the first assignment met of an orbit under the group
    /// generated by `generators` (stride: the object's vertices) is its smallest. Marks that orbit.
    fn first(&mut self, generators: &[u32]) -> bool {
        let n = self.start.len();
        let fixes = |g: &[u32]| {
            (0..n).all(|u| self.choice[self.block(g[u] as usize)] == self.choice[self.block(u)])
        };
        if generators.chunks_exact(n).all(fixes) {
            // a fixed assignment is alone in its orbit
            return true;
        }
        if !self.marked.insert(self.choice.clone()) {
            return false;
        }
        let mut stack = std::mem::take(&mut self.stack);
        stack.push(self.choice.clone());
        while let Some(choice) = stack.pop() {
            for g in generators.chunks_exact(n) {
                let mut image = choice.clone();
                for u in 0..n {
                    image[self.block(g[u] as usize)].copy_from_slice(&choice[self.block(u)]);
                }
                if !self.marked.contains(&image) {
                    self.marked.insert(image.clone());
                    stack.push(image);
                }
            }
        }
        self.stack = stack;
        true
    }
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
    keys: Vec<((u64, u64, u64, u64), u32)>,
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
    pub(crate) filter: Option<&'a dyn Prune>,
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

    /// The accepted children of `topology` after the first `skip` (a complete one is its own only
    /// child), searched with `buf`.
    pub(crate) fn children(
        cx: &Context<'_>,
        topology: Topology,
        skip: usize,
        buf: &mut Buffers,
    ) -> Vec<Topology> {
        let mut search = Search {
            children_only: true,
            ..Search::new(cx.model, topology, std::mem::take(buf))
        };
        let mut children = vec![];
        let mut found = 0;
        while search.next(cx) {
            if found >= skip {
                children.push(search.topology().clone());
            }
            found += 1;
        }
        *buf = search.buf;
        children
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
                    let root = self.buf.level(0);
                    root.known = false;
                    root.vetoed = self.topology.vetoed;
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
                    self.undo(cx.model, depth);
                    if cx.stop.load(Ordering::Relaxed) {
                        self.state = State::Done;
                        return false;
                    }
                    let child = depth + 1;
                    let Level { mode, done, .. } = self.buf.levels[depth];
                    let advanced = !done
                        && match mode {
                            Mode::Vertex => self.next_augmentation(cx, depth),
                            Mode::Pendants => self.next_pendants(cx, depth),
                            Mode::Label => self.next_label(cx.model, depth),
                        };
                    if !advanced {
                        self.state = match depth {
                            0 => State::Done,
                            _ => State::At(depth - 1),
                        };
                        continue;
                    }
                    if mode == Mode::Vertex {
                        let level = &mut self.buf.levels[depth];
                        self.topology.push(&level.corolla, &level.counts);
                        level.pushed = true;
                        self.buf.level(child).known = false;
                        if !self.accepted(cx.limits, depth) {
                            continue;
                        }
                    }
                    if let Some(filter) = cx.filter {
                        self.buf.levels[child].vetoed = self.buf.levels[depth].vetoed;
                        if !self.filter_keeps(filter, child) {
                            continue;
                        }
                    }
                    self.buf.levels[depth].children += 1;
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

    /// Whether `filter` keeps the object at `depth`. While checking heredity, a rejected partial
    /// topology is marked instead of pruned, and a complete topology that the filter keeps must
    /// have no rejected part.
    fn filter_keeps(&mut self, filter: &dyn Prune, depth: usize) -> bool {
        let keeps = filter.keeps(View::Topology(&self.topology));
        if !filter.checks_heredity() {
            return keeps;
        }
        let vetoed = &mut self.buf.levels[depth].vetoed;
        if !self.topology.is_complete() {
            *vetoed |= !keeps;
            return true;
        }
        assert!(
            !(keeps && *vetoed),
            "the filter is not hereditary: it accepts a complete topology but rejected a partial \
             topology it completes"
        );
        keeps
    }

    /// Record the object at `depth` as found, with its generators if it is a complete topology.
    fn found(&mut self, depth: usize) {
        self.found = depth;
        self.topology.vetoed = self.buf.levels[depth].vetoed;
        if self.need_generators && !self.children_only && !self.buf.levels[depth].known {
            self.canonize(depth, false);
        }
    }

    /// Store the node automorphism generators of the object at `depth`, and with `labelling` its
    /// canonical labelling (`canon.lab()`). Node colours come from the vertex invariants and the
    /// labelled legs, edges are coloured by their multiplicity. The canonizer runs only if the
    /// labelling is asked for or some colour class is not a class of twins: distinct colours
    /// make the object rigid, and since automorphisms preserve colours, the transpositions of
    /// twins then generate the group.
    fn canonize(&mut self, depth: usize, labelling: bool) {
        let t = &self.topology;
        let n = t.n;
        let Buffers {
            canon,
            levels,
            invariants,
            keys,
            colours,
            edges,
            ..
        } = &mut self.buf;
        invariants.clear();
        invariants.extend((0..n).map(|v| t.invariant(v)));
        keys.clear();
        keys.extend((0..n).map(|v| {
            let (packed, legs) = invariants[v];
            let refined = t.refined(v, invariants);
            ((packed, legs, t.vertices[v].labelled, refined), v as u32)
        }));
        keys.sort_unstable();
        let level = &mut levels[depth];
        level.generators.clear();
        level.known = true;
        colours.clear();
        colours.resize(n, 0);
        let (mut rank, mut first, mut twins) = (0, 0, true);
        for i in 0..n {
            if i > 0 && keys[i].0 != keys[i - 1].0 {
                rank += 1;
                first = i;
            }
            let (a, b) = (keys[first].1, keys[i].1);
            colours[b as usize] = rank;
            if i > first && twins {
                twins = t.twins(a as usize, b as usize);
                level.generators.extend((0..n as u32).map(|x| match x {
                    _ if x == a => b,
                    _ if x == b => a,
                    _ => x,
                }));
            }
        }
        // distinct colours leave no tie to break, so `twins` holds then too
        if twins && !labelling {
            return;
        }
        level.generators.clear();
        edges.clear();
        for u in 0..n {
            for v in VertexSet(t.vertices[u].nbr & above(u)).iter() {
                edges.push([u as u32, v as u32, t.multiplicity(u, v) as u32]);
            }
        }
        canon.run(colours, edges);
        level.generators.extend(canon.generators().flatten());
    }

    /// Set up the enumeration of the children of the object at `depth`: labels for a shaped
    /// object, augmentations otherwise.
    fn prepare(&mut self, cx: &Context<'_>, depth: usize) {
        if self.topology.is_shaped() {
            self.prepare_labels(cx.model, depth);
        } else {
            self.prepare_augmentations(cx, depth);
        }
    }

    /// The vertices that can take the next label of the shaped object at `depth`: one per orbit
    /// of the object's node automorphism group, the group fixing the vertices of the labels
    /// assigned so far. Labels are assigned in a fixed order, so choosing each one up to the
    /// stabilizer of the previous ones yields every labelling once up to isomorphism.
    fn prepare_labels(&mut self, model: &Model, depth: usize) {
        if !self.buf.levels[depth].known {
            self.canonize(depth, false);
        }
        let t = &self.topology;
        let n = t.n;
        let (slot, _) = model.labelling[t.next_label(model)];
        let Buffers { orbits, levels, .. } = &mut self.buf;
        let level = &mut levels[depth];
        orbits.reset(n);
        for g in level.generators.chunks_exact(n) {
            orbits.join(g);
        }
        level.targets.clear();
        level.targets.extend((0..n).filter(|&v| {
            let vertex = &t.vertices[v];
            model.slots.count(vertex.legs, slot) > model.slots.count(vertex.labelled_slots, slot)
                && orbits.find(v as u32) == v as u32
        }));
        level.mode = Mode::Label;
        level.next = 0;
        level.pushed = false;
        level.children = 0;
        level.done = false;
    }

    /// Give the next label of the shaped object at `depth` to its next target vertex; `false`
    /// after the last. The child is rigid if the object is.
    fn next_label(&mut self, model: &Model, depth: usize) -> bool {
        let level = &mut self.buf.levels[depth];
        let Some(&v) = level.targets.get(level.next) else {
            return false;
        };
        level.next += 1;
        level.pushed = true;
        let rigid = level.generators.is_empty();
        let label = model.labelling[self.topology.next_label(model)];
        self.topology.assign(model, v, label);
        let child = self.buf.level(depth + 1);
        child.known = rigid;
        if rigid {
            child.generators.clear();
        }
        true
    }

    /// Undo the augmentation or label applied at `depth`, if any.
    fn undo(&mut self, model: &Model, depth: usize) {
        let level = &mut self.buf.levels[depth];
        if level.pushed {
            level.pushed = false;
            level.unapply(&mut self.topology, model);
        }
    }

    /// Hand the children of the object at the shallowest depth that has some left to the caller,
    /// as objects to search separately: in the search order they come after everything else this
    /// search still finds, which no longer includes them. Their search uses `buf`.
    pub(crate) fn split(&mut self, cx: &Context<'_>, buf: &mut Buffers) -> Vec<Topology> {
        let State::At(deepest) = self.state else {
            return vec![];
        };
        for depth in 0..=deepest {
            if self.buf.levels[depth].done {
                continue;
            }
            // the object at `depth`, and the children it has after those already found
            let mut object = self.topology.clone();
            for level in self.buf.levels[depth..=deepest].iter_mut().rev() {
                if level.pushed {
                    level.unapply(&mut object, cx.model);
                }
            }
            object.vetoed = self.buf.levels[depth].vetoed;
            let rest = Search::children(cx, object, self.buf.levels[depth].children, buf);
            self.buf.levels[depth].done = true;
            if !rest.is_empty() {
                return rest;
            }
        }
        vec![]
    }

    /// Hang the next assignment of pendants on the open stubs of the object at `depth` that is
    /// minimal under the object's group (rule 1) and whose pendants are the last vertices on the
    /// canonical deletion path of the result (rule 2); `false` after the last.
    fn next_pendants(&mut self, cx: &Context<'_>, depth: usize) -> bool {
        loop {
            let level = &mut self.buf.levels[depth];
            let first = level.next == 0;
            level.next = 1;
            if !level
                .pendants
                .next(&cx.model.slots, self.topology.unplaced, first)
            {
                return false;
            }
            if !self.buf.levels[depth].known {
                self.canonize(depth, false);
            }
            let level = &mut self.buf.levels[depth];
            if !level.pendants.first(&level.generators) {
                continue;
            }
            let level = &mut self.buf.levels[depth];
            level.push_pendants(&mut self.topology);
            level.pushed = true;
            self.buf.level(depth + 1).known = false;
            let limits = cx.limits;
            if (limits.max_bridges.is_none() && limits.zero_flow_edges
                || self.topology.bridges_allow(limits, &mut self.buf.bridges))
                && self.pendants_canonical(depth)
            {
                return true;
            }
            let level = &mut self.buf.levels[depth];
            level.pop_pendants(&mut self.topology);
            level.pushed = false;
        }
    }

    /// Rule 2 for the pendants just hung on the object at `depth`: canonical deletions of the
    /// result remove the pendants one by one, so the object is an ancestor of the result on its
    /// canonical construction path. Leaves the result as it found it.
    fn pendants_canonical(&mut self, depth: usize) -> bool {
        let k = self.buf.levels[depth].pendants.owners.len();
        let base = self.topology.n - k;
        let mut order = std::mem::take(&mut self.buf.levels[depth].pendants.order);
        order.clear();
        order.extend(0..k);
        let mut left = k;
        while left > 0 {
            let n = base + left;
            let pendants = (u64::MAX >> (64 - n)) & !((1u64 << base) - 1);
            let Some(p) = self.canonical_deletion(depth + 1, pendants) else {
                break;
            };
            if p != n - 1 {
                self.topology.swap(p, n - 1);
                order.swap(p - base, n - 1 - base);
            }
            self.buf.levels[depth].pop_pendant(&mut self.topology, order[left - 1]);
            left -= 1;
        }
        let canonical = left == 0;
        let level = &mut self.buf.levels[depth];
        for &i in order[..left].iter().rev() {
            level.pop_pendant(&mut self.topology, i);
        }
        level.pendants.order = order;
        level.push_pendants(&mut self.topology);
        self.buf.levels[depth + 1].known = false;
        canonical
    }

    /// The canonical deletion of the object, as rule 2 chooses it (the deletable vertex with the
    /// largest invariant, then second-order invariant, then canonical position), if a vertex in
    /// its orbit is in `candidates`: that vertex. The automorphism generators are stored at
    /// `scratch`.
    fn canonical_deletion(&mut self, scratch: usize, candidates: u64) -> Option<usize> {
        let t = &self.topology;
        let n = t.n;
        let buf = &mut self.buf;
        buf.invariants.clear();
        buf.invariants.extend((0..n).map(|v| t.invariant(v)));
        let mut best = None;
        let mut ties = 0u64;
        for x in (0..n).filter(|&x| t.connected_without(x)) {
            let key = (buf.invariants[x], t.refined(x, &buf.invariants));
            match best.map_or(std::cmp::Ordering::Greater, |b| key.cmp(&b)) {
                std::cmp::Ordering::Greater => {
                    best = Some(key);
                    ties = 1 << x;
                }
                std::cmp::Ordering::Equal => ties |= 1 << x,
                std::cmp::Ordering::Less => {}
            }
        }
        // a vertex with a unique key is alone in its orbit, and twins share one
        let first = ties.trailing_zeros() as usize;
        if VertexSet(ties).iter().all(|x| t.twins(first, x)) {
            return VertexSet(ties & candidates).iter().next();
        }
        self.canonize(scratch, true);
        let lab = self.buf.canon.lab();
        let star = (0..n)
            .rev()
            .map(|q| lab[q] as usize)
            .find(|&v| ties >> v & 1 == 1)
            .expect("a deletable vertex");
        let Buffers { orbits, levels, .. } = &mut self.buf;
        orbits.reset(n);
        for g in levels[scratch].generators.chunks_exact(n) {
            orbits.join(g);
        }
        let orbit = orbits.find(star as u32);
        (0..n).find(|&x| candidates >> x & 1 == 1 && orbits.find(x as u32) == orbit)
    }

    /// Set up the augmentations of the object at `depth`: the open vertices, the invariants for
    /// rule 2 and the corollas, in search order.
    fn prepare_augmentations(&mut self, cx: &Context<'_>, depth: usize) {
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
            .slots
            .subsets(t.unplaced, max_degree as u32, &mut level.subsets);
        level.corollas.clear();
        for &degree in &model.degrees {
            let budget = t.budget - (degree as i32 - 2);
            for &(legs, n_legs) in &level.subsets {
                if n_legs > degree as u32 || model.leg_state(degree, legs).is_none() {
                    continue;
                }
                let max_loops = if limits.self_loops {
                    (degree as u32 - n_legs) / 2
                } else {
                    0
                };
                let unplaced = t.unplaced_count - n_legs;
                let labelled = model.labelled(legs);
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
                        labelled,
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
        level.children = 0;
        level.done = false;
        level.mode = Mode::Vertex;
        // every completion hangs one pendant on each open stub: no further cycle (open stubs and
        // budget take exactly the unplaced legs) and no more new vertices than open stubs
        let open = t.open_total as i32;
        if n > 0
            && open > 0
            && open + t.budget == t.unplaced_count as i32
            && limits.most_new(n, t.budget) == open
        {
            level.mode = Mode::Pendants;
            let pendants = &mut level.pendants;
            pendants.types.clear();
            pendants.types.extend(
                level
                    .corollas
                    .iter()
                    .filter(|c| c.free == 1 && c.loops == 0),
            );
            pendants.owners.clear();
            pendants.start.clear();
            for u in 0..n {
                pendants.start.push(pendants.owners.len());
                let open = t.vertices[u].open as usize;
                pendants.owners.extend(std::iter::repeat_n(u, open));
            }
            level.counts.resize(t.cap, 0);
        }
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
            self.canonize(depth, false);
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
        for x in VertexSet(ties).iter() {
            match t.refined(x, &buf.invariants).cmp(&refined) {
                std::cmp::Ordering::Greater => return false,
                std::cmp::Ordering::Less => ties &= !(1 << x),
                std::cmp::Ordering::Equal => {}
            }
        }
        // twins of the new vertex are in its orbit
        if ties == 0 || VertexSet(ties).iter().all(|x| t.twins(w, x)) {
            return true;
        }
        self.canonize(parent + 1, true);
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
    use super::{Corolla, Legs, Limits, Model, OpenStubs, Pendants, Topology};
    use crate::half_edge::generate::budget::Request;

    #[test]
    fn twins_have_equal_corollas_and_neighbours() {
        // phi^3 with anonymous legs
        let model = Model::new(1, &[vec![0; 3]], &[(0, 1); 5], &Request::default()).unwrap();
        let limits = Limits::new(&model, 3, 3, None, true, true);
        let corolla = |n_legs: u32| Corolla {
            degree: 3,
            legs: n_legs as u64 * model.slots.one(0),
            n_legs,
            free: 3 - n_legs as u8,
            ..Corolla::default()
        };
        // adjacent, with one leg and one open stub each
        let mut t = Topology::empty(&model, &limits);
        t.push(&corolla(1), &[]);
        t.push(&corolla(1), &[1]);
        assert!(t.twins(0, 1));
        // the path 0 - 1 - 2 with two legs at its ends and one in the middle
        let mut t = Topology::empty(&model, &limits);
        t.push(&corolla(2), &[]);
        t.push(&corolla(1), &[1]);
        t.push(&corolla(2), &[0, 1]);
        assert!(t.twins(0, 2) && t.twins(2, 0));
        assert!(!t.twins(0, 1) && !t.twins(1, 2));
    }

    #[test]
    fn twins_have_equal_multiplicities() {
        // phi^4: vertices 2 and 3 share the neighbours 0 and 1, which look alike, but with swapped
        // multiplicities, so exchanging 2 and 3 alone is not an automorphism
        let model = Model::new(1, &[vec![0; 4]], &[], &Request::default()).unwrap();
        let limits = Limits::new(&model, 8, 4, None, true, true);
        let phi4 = Corolla {
            degree: 4,
            free: 4,
            ..Corolla::default()
        };
        let mut t = Topology::empty(&model, &limits);
        t.push(&phi4, &[]);
        t.push(&phi4, &[0]);
        t.push(&phi4, &[2, 1]);
        t.push(&phi4, &[1, 2, 0]);
        assert_eq!(t.invariant(0), t.invariant(1));
        assert_eq!(t.invariant(2), t.invariant(3));
        assert!(!t.twins(2, 3) && !t.twins(0, 1));
    }

    #[test]
    fn pendant_assignments_are_kept_once_per_orbit() {
        // two vertices with two open stubs each, exchanged by the group; pendants with one or two
        // legs, six legs in all
        let legs = Legs::new(&[(0, 6)]).unwrap();
        let mut pendants = Pendants {
            types: (1..=2)
                .map(|k| Corolla {
                    legs: k,
                    ..Corolla::default()
                })
                .collect(),
            owners: vec![0, 0, 1, 1],
            start: vec![0, 2],
            ..Pendants::default()
        };
        for (generators, expected) in [
            (&[][..], &[[0, 0, 1, 1], [0, 1, 0, 1], [1, 1, 0, 0]][..]),
            (&[1, 0], &[[0, 0, 1, 1], [0, 1, 0, 1]]),
        ] {
            let mut kept = vec![];
            let mut first = true;
            while pendants.next(&legs, legs.all, first) {
                first = false;
                if pendants.first(generators) {
                    kept.push(pendants.choice.clone());
                }
            }
            assert_eq!(kept, expected);
        }
    }

    #[test]
    fn pendant_assignments_use_the_legs_exactly_once_per_multiset() {
        // one leg class with five legs; pendants carrying one, two or three legs; stubs 0 and 1
        // on vertex 0, stub 2 on vertex 1
        let legs = Legs::new(&[(0, 5)]).unwrap();
        let mut pendants = Pendants {
            types: (1..=3)
                .map(|k| Corolla {
                    legs: k,
                    ..Corolla::default()
                })
                .collect(),
            owners: vec![0, 0, 1],
            ..Pendants::default()
        };
        let mut seen = vec![];
        let mut first = true;
        while pendants.next(&legs, legs.all, first) {
            first = false;
            seen.push(pendants.choice.clone());
        }
        // legs per type are 1, 2, 3; the two pendants on vertex 0 are a multiset
        let mut expected = vec![];
        for a in 0..3 {
            for b in a..3 {
                for c in 0..3 {
                    if a + b + c + 3 == 5 {
                        expected.push(vec![a, b, c]);
                    }
                }
            }
        }
        assert_eq!(seen, expected);
    }

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
