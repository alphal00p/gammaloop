//! Stage 2: fields, orientations and vertex rules on a complete topology, one diagram per orbit of
//! the topology's node automorphism group `H`.
//!
//! Edges are grouped into bundles (all edges between one vertex pair, or all self-loops at one
//! vertex) and labelled edge by edge with non-decreasing labels inside a bundle, so permutations
//! of parallel edges and self-loop flips never produce duplicates. Bundles are ordered
//! breadth-first from the legs so that vertices complete early; the order is fixed, since the
//! symmetry breaking below compares labellings in it.
//!
//! A vertex is a state of the model's vertex automaton, and an edge tries only the labels that fit
//! the states of both of its ends (an intersection of their label masks). Labelling an edge moves
//! its ends to new states; the states before are kept on a trail for backtracking. Forward
//! checking then tests the later edges at these two vertices: each needs a label that fits its
//! ends and is not below the labels before it in its bundle, else the partial labelling has no
//! completion and is dropped at once. No other edge can lose a label, as its ends keep their
//! states. Before the first label, every edge needs a label that fits its ends with just their
//! legs. The edges of a bundle share their ends and lower bound, so one of them is tested for all.
//!
//! A complete labelling is kept if it is lexicographically minimal among its images under `H`; the
//! vertex rules are then chosen as minimal tuples under the stabilizer of the labelling. The
//! generator's filter sees every partial labelling that passes forward checking and prunes the
//! labellings it rejects with all their extensions.

use std::cmp::Ordering;

use super::{
    model::{Fitting, Model, VertexState},
    topology::{Group, Topology},
    view::{Prune, State},
};

/// One edge of the topology: `u <= v`, a self-loop when equal.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Edge {
    pub(crate) u: u8,
    pub(crate) v: u8,
    bundle: u16,
}

/// A diagram on a complete topology: labels on a prefix of its edges, and vertex rules once all
/// edges are labelled and the rules are chosen.
#[derive(Clone, Copy)]
pub(crate) struct Decoration<'d> {
    pub(crate) topology: &'d Topology,
    pub(crate) edges: &'d [Edge],
    /// Per labelled edge: a [`super::model::Label`] index, or a [`super::model::LoopLabel`] index
    /// for self-loops.
    pub(crate) labels: &'d [u8],
    /// Per vertex: the automaton state of its legs and labelled edges.
    pub(crate) states: &'d [VertexState],
    /// Per vertex: the signature index; empty until the rules are chosen.
    pub(crate) rules: &'d [u32],
    /// While checking the filter's heredity: the filter rejected a partial labelling of this
    /// diagram.
    pub(crate) vetoed: bool,
}

#[derive(Debug, Default)]
pub(crate) struct Decorator {
    edges: Vec<Edge>,
    /// Edge ranges of the bundles.
    bundles: Vec<(usize, usize)>,
    /// The first edge of each bundle.
    heads: Vec<Edge>,
    bundle_at: Vec<u16>,
    rank: Vec<usize>,
    /// Bundle sort keys: (later rank, earlier rank, u, v).
    pairs: Vec<(usize, usize, usize, usize)>,
    queue: Vec<usize>,
    states: Vec<VertexState>,
    /// Per labelled edge: the states of its ends before its label.
    trail: Vec<[VertexState; 2]>,
    /// The first edges of the later bundles sharing a vertex with bundle `b` are
    /// `later[later_at[b]..later_at[b + 1]]`: forward checking tests them once an edge of `b` is
    /// labelled.
    later: Vec<Edge>,
    later_at: Vec<usize>,
    labels: Vec<u8>,
    group: Group,
    /// Per non-identity group element and bundle: the preimage bundle and whether the pair flips.
    preimage: Vec<u16>,
    flips: Vec<bool>,
    stabilizer: Vec<usize>,
    image: Vec<u8>,
    choices: Vec<usize>,
    rules: Vec<u32>,
    rule_image: Vec<u32>,
    phase: Phase,
    /// Number of edges labelled so far, and whether the labelling extends (or backtracks next).
    assigned: usize,
    descending: bool,
    /// While checking the filter's heredity: the first edge whose partial labelling the filter
    /// rejected.
    vetoed: Option<usize>,
}

/// Where the decoration of a topology stands.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum Phase {
    #[default]
    Done,
    /// The one diagram of a single-label model is ready.
    Single,
    /// Searching labellings.
    Labels,
    /// Enumerating the rule tuples of a minimal labelling; `true` before the first.
    Rules(bool),
}

impl Decorator {
    /// Start decorating `topology`, whose node automorphism generators (stride `n`) are given
    /// unless the model allows a single diagram per topology. [`Self::next`] yields one diagram
    /// per isomorphism class.
    pub(crate) fn start(&mut self, model: &Model, topology: &Topology, generators: &[u32]) {
        self.phase = Phase::Done;
        self.vetoed = None;
        let n = topology.n;
        let single = model.single_diagram();
        if single {
            self.list_edges(topology);
        } else {
            self.order_bundles(topology);
        }
        self.states.clear();
        for vertex in &topology.vertices[..n] {
            let Some(state) = model.leg_state(vertex.degree, vertex.legs) else {
                return;
            };
            self.states.push(state);
        }
        self.labels.clear();
        self.labels.resize(self.edges.len(), 0);
        if single {
            if self.single(model) {
                self.phase = Phase::Single;
            }
            return;
        }
        // forward checking from the legs on: every bundle needs a label that fits its ends
        if self
            .heads
            .iter()
            .any(|&f| self.fitting(model, f, 0).next().is_none())
        {
            return;
        }
        self.list_later();
        self.trail.resize(self.edges.len(), [0; 2]);
        self.group.close(n, generators);
        self.lift();
        self.assigned = 0;
        self.descending = true;
        self.phase = Phase::Labels;
    }

    /// The next diagram on the topology given to [`Self::start`], with its partial labellings
    /// pruned by `filter`.
    pub(crate) fn next<'d>(
        &'d mut self,
        model: &Model,
        t: &'d Topology,
        filter: Option<&dyn Prune>,
    ) -> Option<Decoration<'d>> {
        loop {
            match self.phase {
                Phase::Done => return None,
                Phase::Single => {
                    self.phase = Phase::Done;
                    break;
                }
                Phase::Rules(first) => {
                    if self.next_rules(model, first) {
                        self.phase = Phase::Rules(false);
                        break;
                    }
                    self.phase = Phase::Labels;
                }
                Phase::Labels => {
                    // without a filter, every labelling is kept at no cost
                    let labelled = match filter {
                        None => self.next_labelling(model, &mut |_, _| true),
                        Some(filter) => {
                            self.next_labelling(model, &mut |d: &mut Self, e| d.keeps(t, filter, e))
                        }
                    };
                    if !labelled {
                        self.phase = Phase::Done;
                    } else if self.labelling_minimal(model) {
                        self.choices.clear();
                        self.choices.resize(t.n, 0);
                        self.phase = Phase::Rules(true);
                    }
                }
            }
        }
        Some(Decoration {
            rules: &self.rules,
            vetoed: self.vetoed.is_some(),
            ..self.labelled(t, self.edges.len())
        })
    }

    /// The diagram with its first `n` edges labelled and no rules.
    fn labelled<'d>(&'d self, t: &'d Topology, n: usize) -> Decoration<'d> {
        Decoration {
            topology: t,
            edges: &self.edges,
            labels: &self.labels[..n],
            states: &self.states,
            rules: &[],
            vetoed: false,
        }
    }

    /// Label every edge with the only label and every vertex with its only rule; `false` if the
    /// label does not fit an edge or a vertex has no rule.
    fn single(&mut self, model: &Model) -> bool {
        for e in 0..self.edges.len() {
            let edge = self.edges[e];
            let Some(l) = self.fitting(model, edge, 0).next() else {
                return false;
            };
            self.add(model, edge, l, self.ends(edge));
        }
        self.rules.clear();
        for &state in &self.states {
            let Some(&rule) = model.rules(state).first() else {
                return false;
            };
            self.rules.push(rule);
        }
        true
    }

    /// Edges by vertex, without bundle order.
    fn list_edges(&mut self, t: &Topology) {
        self.edges.clear();
        for u in 0..t.n {
            for v in u..t.n {
                let m = if u == v {
                    t.vertices[u].loops
                } else {
                    t.multiplicity(u, v)
                };
                for _ in 0..m {
                    self.edges.push(Edge {
                        u: u as u8,
                        v: v as u8,
                        bundle: 0,
                    });
                }
            }
        }
    }

    /// Bundles in breadth-first order from the vertices with legs (or vertex 0), each bundle's
    /// edges consecutive.
    fn order_bundles(&mut self, t: &Topology) {
        let n = t.n;
        self.rank.clear();
        self.rank.resize(n, usize::MAX);
        self.queue.clear();
        self.queue
            .extend((0..n).filter(|&v| t.vertices[v].legs != 0));
        if self.queue.is_empty() {
            self.queue.push(0);
        }
        for (r, &v) in self.queue.iter().enumerate() {
            self.rank[v] = r;
        }
        let mut i = 0;
        while i < self.queue.len() {
            let x = self.queue[i];
            i += 1;
            for y in 0..n {
                if y != x && t.multiplicity(x, y) > 0 && self.rank[y] == usize::MAX {
                    self.rank[y] = self.queue.len();
                    self.queue.push(y);
                }
            }
        }
        let pairs = &mut self.pairs;
        pairs.clear();
        for u in 0..n {
            if t.vertices[u].loops > 0 {
                pairs.push((self.rank[u], self.rank[u], u, u));
            }
            for v in u + 1..n {
                if t.multiplicity(u, v) > 0 {
                    let (a, b) = (self.rank[u], self.rank[v]);
                    pairs.push((a.max(b), a.min(b), u, v));
                }
            }
        }
        pairs.sort_unstable();
        self.edges.clear();
        self.bundles.clear();
        self.heads.clear();
        self.bundle_at.clear();
        self.bundle_at.resize(n * n, u16::MAX);
        for (b, &(_, _, u, v)) in self.pairs.iter().enumerate() {
            let m = if u == v {
                t.vertices[u].loops
            } else {
                t.multiplicity(u, v)
            };
            let head = Edge {
                u: u as u8,
                v: v as u8,
                bundle: b as u16,
            };
            let start = self.edges.len();
            self.edges.extend(std::iter::repeat_n(head, m as usize));
            self.bundles.push((start, self.edges.len()));
            self.heads.push(head);
            self.bundle_at[u * n + v] = b as u16;
        }
    }

    /// Per bundle, the first edges of the later bundles that share a vertex with it.
    fn list_later(&mut self) {
        let heads = &self.heads;
        self.later.clear();
        self.later_at.clear();
        self.later_at.push(0);
        for (b, a) in heads.iter().enumerate() {
            let shares = |f: &&Edge| a.u == f.u || a.u == f.v || a.v == f.u || a.v == f.v;
            self.later.extend(heads[b + 1..].iter().filter(shares));
            self.later_at.push(self.later.len());
        }
    }

    /// Bundle permutations of the non-identity group elements.
    fn lift(&mut self) {
        let nb = self.bundles.len();
        let n = self.rank.len();
        self.preimage.clear();
        self.flips.clear();
        for g in self.group.elements().skip(1) {
            let start = self.preimage.len();
            self.preimage.resize(start + nb, 0);
            self.flips.resize(start + nb, false);
            for (i, &Edge { u, v, .. }) in self.heads.iter().enumerate() {
                let (a, b) = (g[u as usize] as usize, g[v as usize] as usize);
                let j = self.bundle_at[a.min(b) * n + a.max(b)] as usize;
                self.preimage[start + j] = i as u16;
                self.flips[start + j] = a > b;
            }
        }
    }

    /// Advance to the next complete labelling whose vertices all fit a signature and whose
    /// prefixes are kept (`keeps(self, e)` judges the labelling up to edge `e`), depth-first over
    /// the edges; `false` after the last.
    fn next_labelling(
        &mut self,
        model: &Model,
        keeps: &mut impl FnMut(&mut Self, usize) -> bool,
    ) -> bool {
        loop {
            if self.descending {
                let e = self.assigned;
                if e == self.edges.len() {
                    self.descending = false;
                    return true;
                }
                let bundle = self.edges[e].bundle;
                let first = if e > 0 && self.edges[e - 1].bundle == bundle {
                    self.labels[e - 1]
                } else {
                    0
                };
                self.descending = self.place(model, e, first as usize, keeps);
            } else {
                // backtrack: the last labelled edge takes its next label
                let Some(e) = self.assigned.checked_sub(1) else {
                    return false;
                };
                self.assigned = e;
                self.set_ends(self.edges[e], self.trail[e]);
                self.descending = self.place(model, e, self.labels[e] as usize + 1, keeps);
            }
        }
    }

    /// Give edge `e` the first label from `first` on that fits its ends, passes forward checking
    /// and that `keeps`.
    fn place(
        &mut self,
        model: &Model,
        e: usize,
        first: usize,
        keeps: &mut impl FnMut(&mut Self, usize) -> bool,
    ) -> bool {
        let edge = self.edges[e];
        let before = self.ends(edge);
        self.trail[e] = before;
        for l in self.fitting(model, edge, first) {
            self.add(model, edge, l, before);
            self.labels[e] = l as u8;
            if self.forward_check(model, e) && keeps(self, e) {
                self.assigned = e + 1;
                return true;
            }
        }
        self.set_ends(edge, before);
        false
    }

    /// The states of the ends of `edge`.
    #[inline(always)]
    fn ends(&self, edge: Edge) -> [VertexState; 2] {
        [self.states[edge.u as usize], self.states[edge.v as usize]]
    }

    /// Put the ends of `edge` into these states (equal for a self-loop).
    #[inline(always)]
    fn set_ends(&mut self, edge: Edge, [a, b]: [VertexState; 2]) {
        self.states[edge.u as usize] = a;
        self.states[edge.v as usize] = b;
    }

    /// The labels from `from` on that fit `edge` between the current states of its ends.
    #[inline(always)]
    fn fitting<'m>(&self, model: &'m Model, edge: Edge, from: usize) -> Fitting<'m> {
        let [a, b] = self.ends(edge);
        if edge.u == edge.v {
            model.fitting_loops(a, from)
        } else {
            model.fitting([a, b], from)
        }
    }

    /// Put the ends of `edge` into the states `before` with label `l` added.
    #[inline(always)]
    fn add(&mut self, model: &Model, edge: Edge, l: usize, before: [VertexState; 2]) {
        let after = if edge.u == edge.v {
            [model.add_loop(before[0], l); 2]
        } else {
            model.add_label(before, l)
        };
        self.set_ends(edge, after);
    }

    /// Every later edge at the ends of edge `e`, just labelled, still has a label that fits its
    /// ends and does not decrease along the bundle of `e`. The edges of a bundle share their ends
    /// and lower bound, so the next edge of the bundle of `e` (from the label of `e` on) and the
    /// first edge of every later bundle at these ends (from any label) stand for all.
    #[inline(always)]
    fn forward_check(&self, model: &Model, e: usize) -> bool {
        let b = self.edges[e].bundle;
        if let Some(&next) = self.edges.get(e + 1) {
            let from = self.labels[e] as usize;
            if next.bundle == b && self.fitting(model, next, from).next().is_none() {
                return false;
            }
        }
        for &f in &self.later[self.later_at[b as usize]..self.later_at[b as usize + 1]] {
            if self.fitting(model, f, 0).next().is_none() {
                return false;
            }
        }
        true
    }

    /// Whether `filter` keeps the labelling up to edge `e`. While checking heredity, a rejection
    /// is recorded instead of pruning.
    fn keeps(&mut self, t: &Topology, filter: &dyn Prune, e: usize) -> bool {
        // a rejection recorded at `e` or later was undone by backtracking
        if self.vetoed.is_some_and(|v| v >= e) {
            self.vetoed = None;
        }
        if filter.keeps(State::Decoration(self.labelled(t, e + 1))) {
            return true;
        }
        if filter.checks_heredity() {
            self.vetoed.get_or_insert(e);
            return true;
        }
        false
    }

    /// Compare the image of the labelling under non-identity element `k` with the labelling.
    fn compare_image(&mut self, model: &Model, k: usize) -> Ordering {
        let nb = self.bundles.len();
        for j in 0..nb {
            let i = self.preimage[k * nb + j] as usize;
            let (start, end) = self.bundles[i];
            let loop_bundle = self.edges[start].u == self.edges[start].v;
            self.image.clear();
            self.image.extend(self.labels[start..end].iter().map(|&l| {
                if self.flips[k * nb + j] && !loop_bundle {
                    model.labels[l as usize].flip
                } else {
                    l
                }
            }));
            self.image.sort_unstable();
            let (start, end) = self.bundles[j];
            match self.image[..].cmp(&self.labels[start..end]) {
                Ordering::Equal => {}
                other => return other,
            }
        }
        Ordering::Equal
    }

    /// The complete labelling is minimal among its images; records its stabilizer.
    fn labelling_minimal(&mut self, model: &Model) -> bool {
        self.stabilizer.clear();
        for k in 0..self.group.len() - 1 {
            match self.compare_image(model, k) {
                Ordering::Less => return false,
                Ordering::Equal => self.stabilizer.push(k + 1),
                Ordering::Greater => {}
            }
        }
        true
    }

    /// Advance to the next rule tuple (the first if `first`) that is minimal under the
    /// stabilizer of the labelling; `false` after the last.
    fn next_rules(&mut self, model: &Model, mut first: bool) -> bool {
        let n = self.choices.len();
        loop {
            if !first {
                // next rule tuple
                let mut v = 0;
                loop {
                    if v == n {
                        return false;
                    }
                    self.choices[v] += 1;
                    if self.choices[v] < model.rules(self.states[v]).len() {
                        break;
                    }
                    self.choices[v] = 0;
                    v += 1;
                }
            }
            first = false;
            self.rules.clear();
            self.rules
                .extend((0..n).map(|v| model.rules(self.states[v])[self.choices[v]]));
            if self.rules_minimal() {
                return true;
            }
        }
    }

    /// The rule tuple is minimal among its images under the stabilizer of the labelling.
    fn rules_minimal(&mut self) -> bool {
        let n = self.rules.len();
        self.rule_image.resize(n, 0);
        for &k in &self.stabilizer {
            let g = self.group.elements().nth(k).expect("stabilizer element");
            for (v, &r) in self.rules.iter().enumerate() {
                self.rule_image[g[v] as usize] = r;
            }
            if self.rule_image[..] < self.rules[..] {
                return false;
            }
        }
        true
    }
}
