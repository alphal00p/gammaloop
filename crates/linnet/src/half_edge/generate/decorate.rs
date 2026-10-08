//! Stage 2: fields, orientations and vertex rules on a complete topology, one diagram per orbit of
//! the topology's node automorphism group `H`.
//!
//! Edges are grouped into bundles (all edges between one vertex pair, or all self-loops at one
//! vertex) and labelled edge by edge with non-decreasing labels inside a bundle, so permutations
//! of parallel edges and self-loop flips never produce duplicates. Vertices are checked after
//! every label against the sub-multisets of the signatures of their degree, and bundles are
//! ordered breadth-first from the legs so that vertices complete early. A complete labelling is
//! kept if it is lexicographically minimal among its images under `H`; the vertex rules are then
//! chosen as minimal tuples under the stabilizer of the labelling. The generator's filter sees every
//! partial labelling that fits the signatures and prunes the labellings it rejects with all their
//! extensions.

use std::cmp::Ordering;

use super::{
    model::{Model, StubSet},
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
    /// Per vertex: the stubs of its legs and labelled edges.
    pub(crate) stubs: &'d [StubSet],
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
    bundle_at: Vec<u16>,
    rank: Vec<usize>,
    /// Bundle sort keys: (later rank, earlier rank, u, v).
    pairs: Vec<(usize, usize, usize, usize)>,
    queue: Vec<usize>,
    stubs: Vec<StubSet>,
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
        self.stubs.clear();
        for v in 0..n {
            let vertex = &topology.vertices[v];
            let stubs = model.leg_stubs(vertex.degree, vertex.legs);
            if !model.feasible(&stubs) {
                return;
            }
            self.stubs.push(stubs);
        }
        self.labels.clear();
        self.labels.resize(self.edges.len(), 0);
        if single {
            if self.single(model) {
                self.phase = Phase::Single;
            }
            return;
        }
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
            stubs: &self.stubs,
            rules: &[],
            vetoed: false,
        }
    }

    /// Label every edge with the only label and every vertex with its only rule; `false` if a
    /// vertex has no rule.
    fn single(&mut self, model: &Model) -> bool {
        for edge in &self.edges {
            let (u, v) = (edge.u as usize, edge.v as usize);
            if u == v {
                let Some(label) = model.loop_labels.first() else {
                    return false;
                };
                self.stubs[u].insert(label.stubs[0]);
                self.stubs[u].insert(label.stubs[1]);
            } else {
                let Some(label) = model.labels.first() else {
                    return false;
                };
                self.stubs[u].insert(label.at_u);
                self.stubs[v].insert(label.at_v);
            }
        }
        self.rules.clear();
        for stubs in &self.stubs {
            let Some(&rule) = model.rules(stubs).first() else {
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
        self.bundle_at.clear();
        self.bundle_at.resize(n * n, u16::MAX);
        for (b, &(_, _, u, v)) in self.pairs.iter().enumerate() {
            let m = if u == v {
                t.vertices[u].loops
            } else {
                t.multiplicity(u, v)
            };
            let start = self.edges.len();
            for _ in 0..m {
                self.edges.push(Edge {
                    u: u as u8,
                    v: v as u8,
                    bundle: b as u16,
                });
            }
            self.bundles.push((start, self.edges.len()));
            self.bundle_at[u * n + v] = b as u16;
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
            for (i, &(first, _)) in self.bundles.iter().enumerate() {
                let Edge { u, v, .. } = self.edges[first];
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
                let l = self.labels[e];
                self.unplace(model, e, l);
                self.descending = self.place(model, e, l as usize + 1, keeps);
            }
        }
    }

    /// Give edge `e` the first label from `first` on that keeps its vertices feasible and that
    /// `keeps`.
    fn place(
        &mut self,
        model: &Model,
        e: usize,
        first: usize,
        keeps: &mut impl FnMut(&mut Self, usize) -> bool,
    ) -> bool {
        let Edge { u, v, .. } = self.edges[e];
        let (u, v) = (u as usize, v as usize);
        if u == v {
            for (l, label) in model.loop_labels.iter().enumerate().skip(first) {
                let [a, b] = label.stubs;
                self.stubs[u].insert(a);
                self.stubs[u].insert(b);
                if model.feasible(&self.stubs[u]) {
                    self.labels[e] = l as u8;
                    if keeps(self, e) {
                        self.assigned = e + 1;
                        return true;
                    }
                }
                self.stubs[u].remove(a);
                self.stubs[u].remove(b);
            }
        } else {
            for (l, label) in model.labels.iter().enumerate().skip(first) {
                self.stubs[u].insert(label.at_u);
                self.stubs[v].insert(label.at_v);
                if model.feasible(&self.stubs[u]) && model.feasible(&self.stubs[v]) {
                    self.labels[e] = l as u8;
                    if keeps(self, e) {
                        self.assigned = e + 1;
                        return true;
                    }
                }
                self.stubs[u].remove(label.at_u);
                self.stubs[v].remove(label.at_v);
            }
        }
        false
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

    /// Remove the stubs of label `l` on edge `e`.
    fn unplace(&mut self, model: &Model, e: usize, l: u8) {
        let Edge { u, v, .. } = self.edges[e];
        let (u, v) = (u as usize, v as usize);
        if u == v {
            let [a, b] = model.loop_labels[l as usize].stubs;
            self.stubs[u].remove(a);
            self.stubs[u].remove(b);
        } else {
            let label = model.labels[l as usize];
            self.stubs[u].remove(label.at_u);
            self.stubs[v].remove(label.at_v);
        }
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
                    if self.choices[v] < model.rules(&self.stubs[v]).len() {
                        break;
                    }
                    self.choices[v] = 0;
                    v += 1;
                }
            }
            first = false;
            self.rules.clear();
            self.rules
                .extend((0..n).map(|v| model.rules(&self.stubs[v])[self.choices[v]]));
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
