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
//! The budgets of the request (see [`super::budget`]) bound their values on every partial
//! labelling, one row of intervals per labelled edge, and a labelling whose intervals miss a range
//! is dropped like one that fails forward checking; before the first label, the legs alone must
//! admit the budgets.
//!
//! A complete labelling is kept if it is lexicographically minimal among its images under `H`; the
//! vertex rules are then chosen as minimal tuples under the stabilizer of the labelling whose
//! weights meet the budgets exactly. The weights are invariant under the stabilizer, so the
//! minimal tuple of an orbit meets them if any tuple does. The generator's filter sees every
//! partial labelling that passes forward checking and the budgets, and prunes the labellings it
//! rejects with all their extensions.

use std::cmp::Ordering;

use crate::half_edge::canon::{factorial, scaled};

use super::{
    budget::Lines,
    model::{Fitting, Model, VertexState},
    suspend::{Step, Work},
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
    /// Fine particle labels may be assigned by whole line rather than by edge prefix. MAX is
    /// unresolved; the complete coarse labels still expose its class and orientation.
    pub(crate) particles: bool,
    pub(crate) coarse_labels: Option<&'d [u8]>,
    /// Per vertex: the automaton state of its legs and labelled edges.
    pub(crate) states: &'d [VertexState],
    /// Per vertex: the signature index; empty until the rules are chosen.
    pub(crate) rules: &'d [u32],
    /// Full half-edge group order, meaningful only after choosing the rules.
    pub(crate) automorphism_count: Option<u128>,
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
    /// Per labelled prefix and budget bound: `tallies[e * width + b]` is the interval of the values
    /// of bound `b` over the completions of the labelling of the first `e` edges.
    tallies: Vec<[i64; 2]>,
    /// Per edge: it is a bridge of the topology (for chord counts).
    bridges: Vec<bool>,
    /// Per line count: the lines of the labelled edges.
    lines: Vec<Lines>,
    labels: Vec<u8>,
    group: Group,
    /// Per non-identity group element and bundle: the preimage bundle and whether the pair flips.
    preimage: Vec<u16>,
    flips: Vec<bool>,
    stabilizer: Vec<usize>,
    /// The node stabilizer of the final rule tuple, and the kernel fixing all nodes.
    node_stabilizer: Option<u128>,
    leg_kernel: Option<u128>,
    kernel: Option<u128>,
    image: Vec<u8>,
    choices: Vec<usize>,
    rules: Vec<u32>,
    rule_image: Vec<u32>,
    phase: Phase,
    refining: bool,
    coarse_labels: Vec<u8>,
    coarse_states: Vec<VertexState>,
    coarse_stabilizer: Vec<usize>,
    refine_lines: Vec<RefinementLine>,
    refine_next: Vec<usize>,
    refine_trails: Vec<Vec<(usize, [VertexState; 2])>>,
    refine_depth: usize,
    refine_descending: bool,
    /// Incremental pin propagation is needed only for external nontrivial-class legs.
    pinned_lines: Option<PinnedLines>,
    pin_checkpoints: Vec<usize>,
    refine_image: Vec<(u8, u8)>,
    /// Number of edges labelled so far, and whether the labelling extends (or backtracks next).
    assigned: usize,
    descending: bool,
    /// An unfinished label attempt retains its next candidate and the states before this edge.
    placement: Option<Placement>,
    /// The current candidate passed native feasibility and awaits exactly one pruning verdict.
    parked: bool,
    /// While checking the filter's heredity: the first edge whose partial labelling the filter
    /// rejected.
    vetoed: Option<usize>,
}

#[derive(Clone, Copy, Debug)]
struct Placement {
    edge: usize,
    next_label: usize,
    before: [VertexState; 2],
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
    /// Particle choices on the directed lines of a minimal coarse assignment.
    Refinements,
    /// Enumerating the rule tuples of a minimal labelling; `true` before the first.
    Rules(bool),
}

/// A nontrivial class occupies one diagonal directed line. The original external legs may
/// fix its particle; otherwise every class member is tried once, before the rule fibres.
#[derive(Debug)]
struct RefinementLine {
    class: usize,
    edges: Vec<usize>,
    pin: Option<u8>,
}

/// Rollback connectivity for externally pinned lines. No instance is created for purely
/// internal flavour sums; those components are built once per complete coarse assignment.
#[derive(Debug)]
struct PinnedLines {
    parent: Vec<usize>,
    size: Vec<usize>,
    pin: Vec<Option<u8>>,
    trail: Vec<(usize, usize, usize, Option<u8>)>,
}

impl PinnedLines {
    fn new(pin: Vec<Option<u8>>) -> Self {
        Self {
            parent: (0..pin.len()).collect(),
            size: vec![1; pin.len()],
            pin,
            trail: Vec::new(),
        }
    }

    fn root(&self, mut vertex: usize) -> usize {
        while self.parent[vertex] != vertex {
            vertex = self.parent[vertex];
        }
        vertex
    }

    fn join(&mut self, u: usize, v: usize) -> bool {
        let (mut u, mut v) = (self.root(u), self.root(v));
        if u == v {
            return true;
        }
        if self.pin[u].zip(self.pin[v]).is_some_and(|(a, b)| a != b) {
            return false;
        }
        if self.size[u] < self.size[v] {
            std::mem::swap(&mut u, &mut v);
        }
        for vertex in [u, v] {
            self.trail.push((
                vertex,
                self.parent[vertex],
                self.size[vertex],
                self.pin[vertex],
            ));
        }
        self.parent[v] = u;
        self.size[u] += self.size[v];
        self.pin[u] = self.pin[u].or(self.pin[v]);
        true
    }

    fn undo(&mut self, checkpoint: usize) {
        while self.trail.len() > checkpoint {
            let (vertex, parent, size, pin) = self.trail.pop().unwrap();
            self.parent[vertex] = parent;
            self.size[vertex] = size;
            self.pin[vertex] = pin;
        }
    }
}

impl Decorator {
    /// Start decorating `topology`, whose node automorphism generators (stride `n`) are given
    /// unless the model allows a single diagram per topology. [`Self::next`] yields one diagram
    /// per isomorphism class.
    pub(crate) fn start(
        &mut self,
        model: &Model,
        topology: &Topology,
        generators: &[u32],
        node_group_order: Option<u128>,
    ) {
        self.phase = Phase::Done;
        self.refining = false;
        self.pinned_lines = None;
        self.placement = None;
        self.parked = false;
        self.vetoed = None;
        self.node_stabilizer = Some(1);
        self.leg_kernel = topology.vertices[..topology.n]
            .iter()
            .try_fold(1, |order, v| {
                (0..model.legs.n_classes()).try_fold(order, |order, class| {
                    scaled(
                        Some(order),
                        factorial(model.legs.count(v.labelled, class) as usize),
                    )
                })
            });
        self.kernel = self.leg_kernel;
        let n = topology.n;
        let single = model.refinement.is_none() && model.single_diagram();
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
        self.trail.resize(self.edges.len(), [0; 2]);
        if !self.start_pins(model, topology) {
            return;
        }
        if !self.start_budgets(model, topology) {
            return;
        }
        if single {
            if self.single(model) {
                // All labels and rules are unique: every topology automorphism survives.
                // Do not close the group here (a 13-leaf star already has 13! elements).
                self.node_stabilizer = node_group_order;
                self.kernel = self.kernel_order(model);
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
        self.group.close(n, generators);
        self.lift();
        self.assigned = 0;
        self.descending = true;
        self.phase = Phase::Labels;
    }

    /// The next diagram on the topology given to [`Self::start`], with its partial labellings
    /// pruned by `filter`. Scalar iteration resolves the same parked prefixes as batching,
    /// checking cancellation around each verdict and between bounded work allowances.
    pub(crate) fn next<'d>(
        &'d mut self,
        model: &Model,
        t: &'d Topology,
        filter: Option<&dyn Prune>,
        stopped: impl Fn() -> bool,
    ) -> Option<Decoration<'d>> {
        while !stopped() {
            let step = self.poll(model, t, filter, filter.is_some(), &mut Work::bounded(128));
            if stopped() {
                return None;
            }
            match step {
                Step::Ready => return Some(self.decoration(t)),
                Step::Done => return None,
                Step::Pending => {
                    let filter = filter.expect("a deferred dressing filter");
                    let keep = filter.keeps(State::Decoration(self.pending(t)));
                    if stopped() {
                        return None;
                    }
                    self.resume(model, keep, filter.checks_heredity());
                }
                Step::Yield => {}
            }
        }
        None
    }

    /// Advance within the work allowance. Deferred pruning parks a feasible field assignment
    /// before invoking any hook; the caller supplies its combined verdict through [`Self::resume`].
    pub(crate) fn poll(
        &mut self,
        model: &Model,
        t: &Topology,
        filter: Option<&dyn Prune>,
        defer: bool,
        work: &mut Work,
    ) -> Step {
        if self.parked {
            return Step::Pending;
        }
        loop {
            match self.phase {
                Phase::Done => return Step::Done,
                Phase::Single => {
                    if !work.spend() {
                        return Step::Yield;
                    }
                    self.phase = Phase::Done;
                    return Step::Ready;
                }
                Phase::Rules(first) => match self.poll_rules(model.particles(), first, work) {
                    Step::Ready => return Step::Ready,
                    Step::Yield => return Step::Yield,
                    Step::Done => {
                        self.phase = if self.refining {
                            Phase::Refinements
                        } else {
                            Phase::Labels
                        }
                    }
                    Step::Pending => unreachable!("rule tuples are judged after dressing"),
                },
                Phase::Labels => match self.poll_labelling(model, t, filter, defer, work) {
                    Step::Done => self.phase = Phase::Done,
                    Step::Pending => return Step::Pending,
                    Step::Yield => return Step::Yield,
                    Step::Ready => {
                        if self.labelling_minimal(model) {
                            if model.refinement.is_some() {
                                if self.start_refinement(model, t) {
                                    self.phase = Phase::Refinements;
                                }
                            } else {
                                self.choices.clear();
                                self.choices.resize(t.n, 0);
                                self.phase = Phase::Rules(true);
                            }
                        }
                    }
                },
                Phase::Refinements => match self.poll_refinement(model, t, filter, defer, work) {
                    Step::Pending => return Step::Pending,
                    Step::Yield => return Step::Yield,
                    Step::Ready => {
                        if self.refinement_minimal(model) {
                            self.choices.clear();
                            self.choices.resize(t.n, 0);
                            self.phase = Phase::Rules(true);
                        }
                    }
                    Step::Done => {
                        self.restore_coarse();
                        self.phase = Phase::Labels;
                    }
                },
            }
        }
    }

    /// The complete assignment from the last [`Step::Ready`].
    pub(crate) fn decoration<'d>(&'d self, t: &'d Topology) -> Decoration<'d> {
        Decoration {
            rules: &self.rules,
            automorphism_count: scaled(self.node_stabilizer, self.kernel),
            vetoed: self.vetoed.is_some(),
            ..self.labelled(t, self.edges.len())
        }
    }

    /// The field prefix from the last [`Step::Pending`], before its line counts are committed.
    pub(crate) fn pending<'d>(&'d self, t: &'d Topology) -> Decoration<'d> {
        assert!(self.parked, "a dressing verdict must be pending");
        self.labelled(
            t,
            if self.refining {
                self.edges.len()
            } else {
                self.placement.expect("parked edge").edge + 1
            },
        )
    }

    /// Resolve one parked edge. Rejection resumes at its next label; accepted prefixes commit
    /// line counts once, including rejected prefixes retained for the heredity diagnostic.
    pub(crate) fn resume(&mut self, model: &Model, keep: bool, check_heredity: bool) {
        assert!(self.parked, "a dressing verdict must be pending");
        self.parked = false;
        if self.refining {
            if self.verdict(self.edges.len() + self.refine_depth, keep, check_heredity) {
                self.refine_depth += 1;
                self.refine_descending = true;
            } else {
                self.undo_refinement_line(self.refine_depth);
            }
            return;
        }
        let placement = self.placement.expect("parked edge");
        if self.verdict(placement.edge, keep, check_heredity) {
            self.commit(model, placement.edge);
        } else {
            self.set_ends(self.edges[placement.edge], placement.before);
            self.undo_pins(placement.edge);
        }
    }

    fn commit(&mut self, model: &Model, edge: usize) {
        self.follow_lines(model, edge, false);
        self.assigned = edge + 1;
        self.descending = true;
        self.placement = None;
    }

    /// The diagram with its first `n` edges labelled and no rules.
    fn labelled<'d>(&'d self, t: &'d Topology, n: usize) -> Decoration<'d> {
        Decoration {
            topology: t,
            edges: &self.edges,
            labels: &self.labels[..n],
            particles: self.refining,
            coarse_labels: self.refining.then_some(&self.coarse_labels),
            states: &self.states,
            rules: &[],
            automorphism_count: None,
            vetoed: false,
        }
    }

    /// Label every edge with the only label and every vertex with its only rule; `false` if the
    /// label does not fit an edge, a vertex has no rule or the budgets reject the diagram.
    fn single(&mut self, model: &Model) -> bool {
        for e in 0..self.edges.len() {
            let edge = self.edges[e];
            let Some(l) = self.fitting(model, edge, 0).next() else {
                return false;
            };
            let before = self.ends(edge);
            self.trail[e] = before;
            self.add(model, edge, l, before);
            self.labels[e] = l as u8;
            if !self.within_budgets(model, e) {
                return false;
            }
            self.follow_lines(model, e, false);
        }
        self.rules.clear();
        for &state in &self.states {
            let Some(&rule) = model.rules(state).first() else {
                return false;
            };
            self.rules.push(rule);
        }
        model.budgets.rules_admit(&self.rules)
    }

    /// Start the budget bounds from the legs: the first row of tallies, the bridges for chord
    /// counts and the lines of the line counts. `false` if the budgets admit no labelling.
    fn start_budgets(&mut self, model: &Model, t: &Topology) -> bool {
        let budgets = &model.budgets;
        // the decorator may come from a generation with other budgets
        self.lines
            .resize_with(budgets.n_line_counts(), Lines::default);
        let width = budgets.len();
        if width == 0 {
            return true;
        }
        let n_edges = self.edges.len();
        self.tallies.clear();
        self.tallies.resize((n_edges + 1) * width, [0; 2]);
        let mut chords = n_edges;
        if budgets.counts_chords() {
            self.bridges.clear();
            self.bridges.resize(n_edges, false);
            let (edges, bridges) = (&self.edges, &mut self.bridges);
            t.bridges(|u, v| {
                let (u, v) = (u.min(v) as u8, u.max(v) as u8);
                // a bridge has no parallel edge
                let e = edges.iter().position(|e| (e.u, e.v) == (u, v));
                bridges[e.expect("a bridge is an edge")] = true;
                chords -= 1;
            });
        }
        let k = budgets.n_weights;
        for w in 0..k {
            self.tallies[w] = self.states.iter().fold([0; 2], |[low, high], &state| {
                let [least, most] = model.weight(state, w);
                [low + least as i64, high + most as i64]
            });
        }
        let n_counts = budgets.n_edge_counts();
        for c in 0..n_counts {
            let eligible = if budgets.chords_only(c) {
                chords
            } else {
                n_edges
            };
            self.tallies[k + c] = [0, eligible as i64];
        }
        for (j, lines) in self.lines.iter_mut().enumerate() {
            let legs = |v: usize| {
                (0..model.legs.n_classes())
                    .filter(|&c| budgets.leg_lines(c) >> j & 1 == 1)
                    .map(|c| model.legs.count(t.vertices[v].labelled, c))
                    .sum()
            };
            lines.reset((0..t.n).map(legs));
            self.tallies[k + n_counts + j] = Lines::interval(
                lines.closed,
                lines.branching,
                n_edges == 0,
                budgets.branching(j),
            );
        }
        budgets.admit(&self.tallies[..width])
    }

    /// Whether the budgets admit the labelling up to edge `e`, just labelled.
    #[inline(always)]
    fn within_budgets(&mut self, model: &Model, e: usize) -> bool {
        model.budgets.is_empty() || self.tally(model, e)
    }

    /// Whether budgets admit the labelling up to edge `e`, just labelled. Computes the bounds
    /// after `e` (row `e + 1` of the tallies) from those before it: the weight bounds of the two
    /// ends move with their states, an edge count counts the edge or no longer can, and the line
    /// counts follow the lines; on the last edge the counts are exact. Kept out of the labelling
    /// loop, which runs without budgets as fast as before.
    #[inline(never)]
    fn tally(&mut self, model: &Model, e: usize) -> bool {
        let budgets = &model.budgets;
        let width = budgets.len();
        let edge = self.edges[e];
        let (u, v) = (edge.u as usize, edge.v as usize);
        let (before, after) = (self.trail[e], self.ends(edge));
        let (done, rest) = self.tallies.split_at_mut((e + 1) * width);
        let (row, next) = (&done[e * width..], &mut rest[..width]);
        let k = budgets.n_weights;
        for w in 0..k {
            let shift = |[low, high]: [i64; 2], from, to| {
                let ([a, b], [c, d]) = (model.weight(from, w), model.weight(to, w));
                [low + (c - a) as i64, high + (d - b) as i64]
            };
            next[w] = shift(row[w], before[0], after[0]);
            if u != v {
                next[w] = shift(next[w], before[1], after[1]);
            }
        }
        let counted = budgets.counted(self.labels[e] as usize, u == v);
        let n_counts = budgets.n_edge_counts();
        for c in 0..n_counts {
            let [low, high] = row[k + c];
            next[k + c] = if budgets.chords_only(c) && self.bridges[e] {
                [low, high]
            } else if counted >> c & 1 == 1 {
                [low + 1, high]
            } else {
                [low, high - 1]
            };
        }
        let complete = e + 1 == self.edges.len();
        for (j, lines) in self.lines.iter().enumerate() {
            let (closed, branching) = if counted >> (n_counts + j) & 1 == 1 {
                lines.after(u, v)
            } else {
                (lines.closed, lines.branching)
            };
            next[k + n_counts + j] =
                Lines::interval(closed, branching, complete, budgets.branching(j));
        }
        budgets.admit(next)
    }

    /// Add edge `e`, labelled and accepted, to the lines of the line counts that its label adds
    /// to, or with `undo` remove it again, the last edge added.
    #[inline(always)]
    fn follow_lines(&mut self, model: &Model, e: usize, undo: bool) {
        if !self.lines.is_empty() {
            self.follow_lines_of(model, e, undo);
        }
    }

    #[inline(never)]
    fn follow_lines_of(&mut self, model: &Model, e: usize, undo: bool) {
        let Edge { u, v, .. } = self.edges[e];
        let budgets = &model.budgets;
        let counted = budgets.counted(self.labels[e] as usize, u == v) >> budgets.n_edge_counts();
        for (j, lines) in self.lines.iter_mut().enumerate() {
            if counted >> j & 1 == 1 {
                if undo {
                    lines.undo(u as usize, v as usize);
                } else {
                    lines.add(u as usize, v as usize);
                }
            }
        }
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

    /// Depth-first field assignment with one work checkpoint per candidate. The placement
    /// cursor survives both a work yield and a deferred verdict, so neither repeats an edge.
    fn poll_labelling(
        &mut self,
        model: &Model,
        t: &Topology,
        filter: Option<&dyn Prune>,
        defer: bool,
        work: &mut Work,
    ) -> Step {
        loop {
            if self.placement.is_none() && self.descending && self.assigned == self.edges.len() {
                self.descending = false;
                return Step::Ready;
            }
            if !work.spend() {
                return Step::Yield;
            }
            let placement = match self.placement {
                Some(placement) => placement,
                None => {
                    let (e, first) = if self.descending {
                        let e = self.assigned;
                        let first = if e > 0 && self.edges[e - 1].bundle == self.edges[e].bundle {
                            self.labels[e - 1] as usize
                        } else {
                            0
                        };
                        (e, first)
                    } else {
                        let Some(e) = self.assigned.checked_sub(1) else {
                            return Step::Done;
                        };
                        self.assigned = e;
                        self.set_ends(self.edges[e], self.trail[e]);
                        self.follow_lines(model, e, true);
                        self.undo_pins(e);
                        (e, self.labels[e] as usize + 1)
                    };
                    let before = self.ends(self.edges[e]);
                    self.trail[e] = before;
                    if let Some(pins) = &self.pinned_lines {
                        self.pin_checkpoints[e] = pins.trail.len();
                    }
                    Placement {
                        edge: e,
                        next_label: first,
                        before,
                    }
                }
            };
            let e = placement.edge;
            let edge = self.edges[e];
            let Some(label) = self.fitting(model, edge, placement.next_label).next() else {
                self.set_ends(edge, placement.before);
                self.placement = None;
                self.descending = false;
                continue;
            };
            self.placement = Some(Placement {
                next_label: label + 1,
                ..placement
            });
            self.add(model, edge, label, placement.before);
            self.labels[e] = label as u8;
            if self.forward_check(model, e)
                && self.within_budgets(model, e)
                && self.join_pins(model, e)
            {
                if defer {
                    self.parked = true;
                    return Step::Pending;
                }
                if filter.is_none_or(|filter| self.keeps(t, filter, e)) {
                    self.commit(model, e);
                    continue;
                }
            }
            self.set_ends(edge, placement.before);
            self.undo_pins(e);
        }
    }

    fn class_field(model: &Model, edge: Edge, label: u8) -> usize {
        if edge.u == edge.v {
            model.theory.loop_labels[label as usize].field as usize
        } else {
            model.theory.labels[label as usize].field as usize
        }
    }

    fn start_pins(&mut self, model: &Model, topology: &Topology) -> bool {
        let Some(refinement) = &model.refinement else {
            return true;
        };
        let fine = &refinement.fine;
        let mut pins = vec![None; topology.n];
        for (vertex, pin) in pins.iter_mut().enumerate() {
            for class in 0..fine.legs.n_classes() {
                if fine.legs.count(topology.vertices[vertex].labelled, class) == 0 {
                    continue;
                }
                let field = fine.legs.stub[class] / 3;
                if refinement.members[refinement.field_class[field as usize] as usize].len() == 1 {
                    continue;
                }
                if pin.is_some_and(|previous| previous != field) {
                    return false;
                }
                *pin = Some(field);
            }
        }
        if pins.iter().any(Option::is_some) {
            self.pinned_lines = Some(PinnedLines::new(pins));
            self.pin_checkpoints.clear();
            self.pin_checkpoints.resize(self.edges.len(), 0);
        }
        true
    }

    fn join_pins(&mut self, model: &Model, edge: usize) -> bool {
        let Some(pins) = &mut self.pinned_lines else {
            return true;
        };
        let e = self.edges[edge];
        let class = Self::class_field(model, e, self.labels[edge]);
        if model.refinement.as_ref().unwrap().members[class].len() == 1 {
            return true;
        }
        pins.join(e.u as usize, e.v as usize)
    }

    fn undo_pins(&mut self, edge: usize) {
        if let Some(pins) = &mut self.pinned_lines {
            pins.undo(self.pin_checkpoints[edge]);
        }
    }

    /// Map a label between the coarse and particle alphabet without changing its orientation.
    fn corresponding_label(
        from: &Model,
        to: &Model,
        edge: Edge,
        label: u8,
        field: u8,
    ) -> Option<u8> {
        if edge.u == edge.v {
            let directed = from.theory.loop_labels[label as usize].directed;
            to.theory
                .loop_labels
                .iter()
                .position(|l| l.field == field && l.directed == directed)
        } else {
            let direction = from.theory.labels[label as usize].direction;
            to.theory
                .labels
                .iter()
                .position(|l| l.field == field && l.direction == direction)
        }
        .map(|index| index as u8)
    }

    fn start_refinement(&mut self, model: &Model, topology: &Topology) -> bool {
        let refinement = model.refinement.as_ref().unwrap();
        let fine = &refinement.fine;
        self.coarse_labels.clone_from(&self.labels);
        self.coarse_states.clone_from(&self.states);
        self.coarse_stabilizer.clone_from(&self.stabilizer);
        self.refining = true;
        self.labels.fill(u8::MAX);
        self.states.clear();
        for vertex in &topology.vertices[..topology.n] {
            let Some(state) = fine.labelled_state(vertex.degree, vertex.labelled) else {
                self.restore_coarse();
                return false;
            };
            self.states.push(state);
        }
        // Unpinned line connectivity is deliberately constructed only here, once per complete
        // coarse assignment. Incremental rollback connectivity exists solely for external pins.
        let mut internal = None;
        if self.pinned_lines.is_none() {
            let mut lines = PinnedLines::new(vec![None; topology.n]);
            for (index, &edge) in self.edges.iter().enumerate() {
                let class = Self::class_field(model, edge, self.coarse_labels[index]);
                if refinement.members[class].len() > 1 {
                    assert!(lines.join(edge.u as usize, edge.v as usize));
                }
            }
            internal = Some(lines);
        }
        let lines = self.pinned_lines.as_ref().or(internal.as_ref()).unwrap();
        self.refine_lines.clear();
        for (index, &edge) in self.edges.iter().enumerate() {
            let class = Self::class_field(model, edge, self.coarse_labels[index]);
            if refinement.members[class].len() == 1 {
                continue;
            }
            let root = lines.root(edge.u as usize);
            if let Some(line) = self
                .refine_lines
                .iter_mut()
                .find(|line| lines.root(self.edges[line.edges[0]].u as usize) == root)
            {
                assert_eq!(
                    line.class, class,
                    "a safe directed line has one field class"
                );
                line.edges.push(index);
            } else {
                self.refine_lines.push(RefinementLine {
                    class,
                    edges: vec![index],
                    pin: lines.pin[root],
                });
            }
        }
        // Singleton edges are already fully known in the coarse view. Their fine automaton
        // transitions are fixed, not new branching decisions or extra hook calls.
        for index in 0..self.edges.len() {
            let edge = self.edges[index];
            let class = Self::class_field(model, edge, self.coarse_labels[index]);
            if refinement.members[class].len() != 1 {
                continue;
            }
            let label = Self::corresponding_label(
                model,
                fine,
                edge,
                self.coarse_labels[index],
                refinement.members[class][0],
            )
            .unwrap();
            let Some(after) = fine.try_label(self.ends(edge), label as usize, edge.u == edge.v)
            else {
                self.restore_coarse();
                return false;
            };
            self.set_ends(edge, after);
            self.labels[index] = label;
        }
        self.refine_next.clear();
        self.refine_next.resize(self.refine_lines.len(), 0);
        self.refine_trails.clear();
        self.refine_trails
            .resize_with(self.refine_lines.len(), Vec::new);
        self.refine_depth = 0;
        self.refine_descending = true;
        true
    }

    fn restore_coarse(&mut self) {
        self.labels.clone_from(&self.coarse_labels);
        self.states.clone_from(&self.coarse_states);
        self.refining = false;
    }

    fn undo_refinement_line(&mut self, depth: usize) {
        while let Some((index, before)) = self.refine_trails[depth].pop() {
            self.set_ends(self.edges[index], before);
            self.labels[index] = u8::MAX;
        }
    }

    fn refine_line(&mut self, model: &Model, depth: usize, field: u8) -> bool {
        let fine = model.particles();
        for position in 0..self.refine_lines[depth].edges.len() {
            let index = self.refine_lines[depth].edges[position];
            let edge = self.edges[index];
            let label =
                Self::corresponding_label(model, fine, edge, self.coarse_labels[index], field)
                    .unwrap();
            let before = self.ends(edge);
            let Some(after) = fine.try_label(before, label as usize, edge.u == edge.v) else {
                return false;
            };
            self.refine_trails[depth].push((index, before));
            self.set_ends(edge, after);
            self.labels[index] = label;
        }
        // Count membership is constant throughout each class. Only rule-weight intervals can
        // narrow during particle refinement; exact rule weights are still checked by Rules.
        let width = model.budgets.len();
        if width != 0 {
            let start = self.edges.len() * width;
            let mut row = self.tallies[start..start + width].to_vec();
            for (weight, interval) in row[..fine.budgets.n_weights].iter_mut().enumerate() {
                *interval = self.states.iter().fold([0; 2], |[a, b], &state| {
                    let [low, high] = fine.weight(state, weight);
                    [a + low as i64, b + high as i64]
                });
            }
            if !fine.budgets.admit(&row) {
                return false;
            }
        }
        true
    }

    fn poll_refinement(
        &mut self,
        model: &Model,
        topology: &Topology,
        filter: Option<&dyn Prune>,
        defer: bool,
        work: &mut Work,
    ) -> Step {
        loop {
            if self.refine_descending && self.refine_depth == self.refine_lines.len() {
                self.refine_descending = false;
                return Step::Ready;
            }
            if !work.spend() {
                return Step::Yield;
            }
            if !self.refine_descending {
                let Some(depth) = self.refine_depth.checked_sub(1) else {
                    return Step::Done;
                };
                self.refine_depth = depth;
                self.undo_refinement_line(depth);
                self.refine_descending = true;
            }
            let depth = self.refine_depth;
            let line = &self.refine_lines[depth];
            let members = &model.refinement.as_ref().unwrap().members[line.class];
            let field = match line.pin {
                Some(pin) => (self.refine_next[depth] == 0).then_some(pin),
                None => members.get(self.refine_next[depth]).copied(),
            };
            let Some(field) = field else {
                self.refine_next[depth] = 0;
                self.refine_descending = false;
                continue;
            };
            self.refine_next[depth] += 1;
            if self.refine_line(model, depth, field) {
                if defer {
                    self.parked = true;
                    return Step::Pending;
                }
                let keep = filter.is_none_or(|filter| {
                    filter.keeps(State::Decoration(self.labelled(topology, self.edges.len())))
                });
                if self.verdict(
                    self.edges.len() + depth,
                    keep,
                    filter.is_some_and(|filter| filter.checks_heredity()),
                ) {
                    self.refine_depth += 1;
                    continue;
                }
            }
            self.undo_refinement_line(depth);
        }
    }

    /// The second quotient uses only the stabilizer of this coarse assignment. Bundle order
    /// remains coarse-first even when a class's original particle ranks are interleaved.
    fn refinement_minimal(&mut self, model: &Model) -> bool {
        let fine = model.particles();
        let field_class = &model.refinement.as_ref().unwrap().field_class;
        self.stabilizer.clear();
        for position in 0..self.coarse_stabilizer.len() {
            let k = self.coarse_stabilizer[position] - 1;
            let mut comparison = Ordering::Equal;
            for bundle in 0..self.bundles.len() {
                let preimage = self.preimage[k * self.bundles.len() + bundle] as usize;
                let (start, end) = self.bundles[preimage];
                let edge = self.edges[start];
                self.refine_image.clear();
                for &label in &self.labels[start..end] {
                    let label = if self.flips[k * self.bundles.len() + bundle] && edge.u != edge.v {
                        fine.theory.labels[label as usize].flip
                    } else {
                        label
                    };
                    let field = Self::class_field(fine, edge, label);
                    let coarse =
                        Self::corresponding_label(fine, model, edge, label, field_class[field])
                            .unwrap();
                    self.refine_image.push((coarse, label));
                }
                self.refine_image.sort_unstable();
                let (start, end) = self.bundles[bundle];
                comparison = self
                    .refine_image
                    .iter()
                    .map(|&(_, label)| label)
                    .cmp(self.labels[start..end].iter().copied());
                if comparison != Ordering::Equal {
                    break;
                }
            }
            match comparison {
                Ordering::Less => return false,
                Ordering::Equal => self.stabilizer.push(k + 1),
                Ordering::Greater => {}
            }
        }
        self.kernel = self.kernel_order(fine);
        true
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
        let keep = filter.keeps(State::Decoration(self.labelled(t, e + 1)));
        self.verdict(e, keep, filter.checks_heredity())
    }

    fn verdict(&mut self, e: usize, keep: bool, check_heredity: bool) -> bool {
        // a rejection recorded at `e` or later was undone by backtracking
        if self.vetoed.is_some_and(|v| v >= e) {
            self.vetoed = None;
        }
        if keep {
            return true;
        }
        if check_heredity {
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
                    model.theory.labels[l as usize].flip
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
        self.kernel = self.kernel_order(model);
        true
    }

    /// The kernel fixing every node: equal edges and legs permute, and undirected loops flip.
    fn kernel_order(&self, model: &Model) -> Option<u128> {
        let mut order = self.leg_kernel;
        let mut start = 0;
        while start < self.edges.len() {
            let edge = self.edges[start];
            let label = self.labels[start];
            let mut end = start + 1;
            // Single uses bundle zero for every edge, so compare actual endpoints and labels.
            while end < self.edges.len()
                && (self.edges[end].u, self.edges[end].v, self.labels[end])
                    == (edge.u, edge.v, label)
            {
                end += 1;
            }
            order = scaled(order, factorial(end - start));
            if edge.u == edge.v && !model.theory.loop_labels[label as usize].directed {
                for _ in start..end {
                    order = scaled(order, Some(2));
                }
            }
            start = end;
        }
        order
    }

    /// Advance to the next minimal rule tuple. A rejected tuple still spends work; preserving
    /// the first-tuple flag prevents its replay when every tuple in a quantum is rejected.
    fn poll_rules(&mut self, model: &Model, mut first: bool, work: &mut Work) -> Step {
        let n = self.choices.len();
        loop {
            if !work.spend() {
                return Step::Yield;
            }
            if !first {
                // next rule tuple
                let mut v = 0;
                loop {
                    if v == n {
                        return Step::Done;
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
            self.phase = Phase::Rules(false);
            self.rules.clear();
            self.rules
                .extend((0..n).map(|v| model.rules(self.states[v])[self.choices[v]]));
            if self.rules_minimal() && model.budgets.rules_admit(&self.rules) {
                return Step::Ready;
            }
        }
    }

    /// The rule tuple is minimal among its images under the stabilizer of the labelling.
    fn rules_minimal(&mut self) -> bool {
        let n = self.rules.len();
        self.rule_image.resize(n, 0);
        self.node_stabilizer = Some(1); // identity is omitted from `stabilizer`
        for &k in &self.stabilizer {
            let g = self.group.elements().nth(k).expect("stabilizer element");
            for (v, &r) in self.rules.iter().enumerate() {
                self.rule_image[g[v] as usize] = r;
            }
            match self.rule_image[..].cmp(&self.rules[..]) {
                Ordering::Less => return false,
                Ordering::Equal => {
                    self.node_stabilizer = self.node_stabilizer.and_then(|n| n.checked_add(1));
                }
                Ordering::Greater => {}
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use std::sync::{atomic::AtomicBool, Mutex};

    use super::*;
    use crate::half_edge::generate::{Countable, Legs, Plan, Query, Structure, Stub, Theory};

    #[derive(Debug, PartialEq, Eq)]
    struct Assignment {
        labels: Vec<u8>,
        rules: Vec<u32>,
        automorphisms: Option<u128>,
        vetoed: bool,
    }

    struct Prefixes {
        reject: bool,
        check: bool,
        seen: Mutex<Vec<Vec<u8>>>,
    }

    impl Prune for Prefixes {
        fn keeps(&self, state: State<'_>) -> bool {
            let State::Decoration(d) = state else {
                panic!("only field prefixes belong to decorator pruning");
            };
            assert!(d.rules.is_empty());
            self.seen.lock().unwrap().push(d.labels.to_vec());
            !self.reject || d.labels.last() == Some(&0)
        }

        fn checks_heredity(&self) -> bool {
            self.check
        }
    }

    fn plan() -> Plan<u8, u8> {
        let (a, b) = (Stub::undirected(0), Stub::undirected(1));
        let theory = Theory::builder()
            .weights(["order"])
            .vertex([a.clone(), a.clone(), a.clone()], [0])
            .vertex([a.clone(), a.clone(), a.clone()], [2])
            .vertex([a.clone(), a.clone(), b.clone()], [1])
            .vertex([a.clone(), b.clone(), b.clone()], [2])
            .vertex([b.clone(), b.clone(), b], [3])
            .build()
            .unwrap();
        let legs = Legs::new().leg(1, a.clone()).leg(2, a);
        let query = Query::new(&legs)
            .loops(1..=1)
            .structure(Structure {
                self_loops: true,
                zero_flow_edges: true,
                ..Structure::default()
            })
            .weight("order", 0..=3)
            .count(Countable::Propagators(vec![1]), ..=1)
            .count(Countable::ClosedLines(vec![1]), ..=1);
        Plan::new(&theory, &query).unwrap()
    }

    /// Compare owned assignments and the exact callback sequence across every topology.
    fn collect(
        plan: &Plan<u8, u8>,
        filter: Option<&Prefixes>,
        deferred: bool,
        quantum: Option<usize>,
    ) -> (Vec<Assignment>, usize) {
        use crate::half_edge::generate::topology::{Context, Search};
        let model = &plan.compiled.model;
        let limits = &plan.orders[0].1;
        let stop = AtomicBool::new(false);
        let cx = Context {
            model,
            limits,
            stop: &stop,
            abort_check: None,
            filter: None,
        };
        let mut search = Search::new(model, Topology::empty(model, limits), Default::default());
        let (mut assignments, mut yields) = (Vec::new(), 0);
        while search.next(&cx) {
            let t = search.topology();
            let mut decorator = Decorator::default();
            decorator.start(model, t, search.generators(), search.node_group_order());
            loop {
                let mut work = quantum.map_or_else(Work::unlimited, Work::bounded);
                match decorator.poll(
                    model,
                    t,
                    filter.map(|f| f as &dyn Prune),
                    deferred,
                    &mut work,
                ) {
                    Step::Ready => {
                        let d = decorator.decoration(t);
                        assignments.push(Assignment {
                            labels: d.labels.to_vec(),
                            rules: d.rules.to_vec(),
                            automorphisms: d.automorphism_count,
                            vetoed: d.vetoed,
                        });
                    }
                    Step::Pending => {
                        let filter = filter.expect("deferred hook");
                        let labels = decorator.pending(t).labels.to_vec();
                        // Polling a parked cursor spends nothing and never invokes the hook.
                        let seen = filter.seen.lock().unwrap().len();
                        assert_eq!(
                            decorator.poll(model, t, Some(filter), true, &mut Work::bounded(0)),
                            Step::Pending
                        );
                        assert_eq!(filter.seen.lock().unwrap().len(), seen);
                        assert_eq!(decorator.pending(t).labels, labels);
                        let keep = filter.keeps(State::Decoration(decorator.pending(t)));
                        decorator.resume(model, keep, filter.check);
                    }
                    Step::Yield => yields += 1,
                    Step::Done => break,
                }
            }
        }
        (assignments, yields)
    }

    #[test]
    fn tiny_work_quanta_preserve_assignments_and_rule_tuple_order() {
        let plan = plan();
        let (expected, _) = collect(&plan, None, false, None);
        assert!(!expected.is_empty());
        for quantum in [1, 2, 3] {
            let (actual, yields) = collect(&plan, None, false, Some(quantum));
            assert_eq!(actual, expected);
            assert!(yields > 0);
        }
    }

    #[test]
    fn native_budget_rejections_yield_without_any_hook_or_output() {
        // A connected diagram with A legs can only use AAA vertices: BBB is disconnected.
        // Degree/weight partitions admit one BBB vertex, but dressing must reject every graph.
        let (a, b) = (Stub::undirected(0), Stub::undirected(1));
        let theory = Theory::builder()
            .weights(["b_vertices"])
            .vertex([a.clone(), a.clone(), a.clone()], [0])
            .vertex([b.clone(), b.clone(), b], [1])
            .build()
            .unwrap();
        let legs = Legs::new().leg(1, a.clone()).leg(2, a);
        let query = Query::new(&legs)
            .loops(2..=2)
            .structure(Structure {
                self_loops: true,
                zero_flow_edges: true,
                ..Structure::default()
            })
            .weight("b_vertices", 1..=1);
        let plan = Plan::new(&theory, &query).unwrap();
        assert!(collect(&plan, None, false, None).0.is_empty());
        let (assignments, yields) = collect(&plan, None, false, Some(1));
        assert!(assignments.is_empty());
        assert!(yields > 0);
    }

    #[test]
    fn deferred_rejections_resume_once_and_preserve_heredity_lineage() {
        let plan = plan();
        for check in [false, true] {
            let direct = Prefixes {
                reject: true,
                check,
                seen: Mutex::default(),
            };
            let (expected, _) = collect(&plan, Some(&direct), false, None);
            assert!(!expected.is_empty());
            assert_eq!(expected.iter().any(|d| d.vetoed), check);
            for quantum in [1, 2, 3] {
                let deferred = Prefixes {
                    reject: true,
                    check,
                    seen: Mutex::default(),
                };
                let (actual, _) = collect(&plan, Some(&deferred), true, Some(quantum));
                assert_eq!(actual, expected);
                assert_eq!(*deferred.seen.lock().unwrap(), *direct.seen.lock().unwrap());
                assert!(!deferred.seen.lock().unwrap().is_empty());
            }
        }
    }
}

#[cfg(test)]
mod refinement_tests {
    use std::{collections::BTreeMap, sync::Mutex};

    use super::*;
    use crate::half_edge::{
        generate::{
            batch::{Event, HookMask, Output},
            Classes, Countable, Dressed, GenerationError, Legs, Plan, Query, Run, Structure, Stub,
            Theory,
        },
        involution::Flow,
    };

    fn theory() -> Theory<char> {
        let mut builder = Theory::builder().weights(["order"]);
        // The two members straddle the fixed boson in field-rank order. Equal signatures are
        // distinct rule fibres, and selecting a weight does not merge their original IDs.
        for field in ['a', 'c'] {
            for weight in [0, 1] {
                builder = builder.vertex(
                    [
                        Stub::incoming(field),
                        Stub::outgoing(field),
                        Stub::undirected('b'),
                    ],
                    [weight],
                );
            }
        }
        builder
            .vertex(vec![Stub::undirected('b'); 3], [0])
            .build()
            .unwrap()
    }

    fn query(legs: &Legs<char, u8>, loops: usize) -> Query<char, u8> {
        Query::new(legs).loops(loops..=loops).structure(Structure {
            self_loops: true,
            zero_flow_edges: true,
            ..Structure::default()
        })
    }

    fn exact(
        theory: &Theory<char>,
        request: &Query<char, u8>,
        mode: Classes,
    ) -> BTreeMap<crate::half_edge::canon::CanonicalForm<usize, char, Option<u8>>, u128> {
        let plan = Plan::new(theory, &request.clone().classes(mode)).unwrap();
        let mut reference = None;
        let executions: &[(bool, usize)] =
            if cfg!(all(feature = "rayon", not(target_arch = "wasm32"))) {
                &[(false, 1), (true, 1), (true, 8), (true, 32)]
            } else {
                &[(false, 1)]
            };
        for &(parallel, threads) in executions {
            let generate = || {
                let run = Run::default().parallel(parallel);
                let entry = |diagram: crate::half_edge::generate::Diagram<char, u8>| {
                    (
                        diagram.canonical_form().clone(),
                        diagram.symmetry_factor().unwrap(),
                    )
                };
                let ordered = plan.diagrams(&run).unwrap().map(entry).collect::<Vec<_>>();
                let diagrams = ordered.iter().cloned().collect::<BTreeMap<_, _>>();
                assert_eq!(diagrams.len(), ordered.len(), "duplicate coloured key");
                let counts = plan.count(&run).unwrap();
                assert_eq!(counts.diagrams, num_bigint::BigUint::from(diagrams.len()));
                let weight: num_rational::BigRational = diagrams
                    .values()
                    .map(|&aut| num_rational::BigRational::new(1.into(), aut.into()))
                    .sum();
                assert_eq!(counts.sum_inv_aut, weight);
                assert!(counts.complete);
                let snapshots = plan
                    .skeletons(&run)
                    .unwrap()
                    .map(|skeleton| skeleton.snapshot())
                    .collect::<Vec<_>>();
                let retained = snapshots
                    .iter()
                    .flat_map(|skeleton| plan.dressings(skeleton, &run).unwrap().map(entry))
                    .collect::<Vec<_>>();
                let retained_map = retained.iter().cloned().collect::<BTreeMap<_, _>>();
                assert_eq!(retained_map.len(), retained.len(), "duplicate retained key");
                assert_eq!(retained_map, diagrams, "retained skeleton dressing");
                ordered
            };
            #[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
            let actual = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(generate);
            #[cfg(not(all(feature = "rayon", not(target_arch = "wasm32"))))]
            let actual = generate();
            if let Some(expected) = &reference {
                assert_eq!(
                    &actual, expected,
                    "{mode:?}: order with parallel={parallel} at {threads} threads"
                );
            } else {
                reference = Some(actual);
            }
        }
        reference.unwrap().into_iter().collect()
    }

    #[test]
    fn class_refinement_preserves_fine_keys_rules_weights_and_external_pins() {
        let theory = theory();
        let photon = Legs::new()
            .leg(1, Stub::undirected('b'))
            .leg(2, Stub::undirected('b'));
        let four_photons = (1..=4).fold(Legs::new(), |legs, label| {
            legs.leg(label, Stub::undirected('b'))
        });
        let pinned = Legs::new()
            .leg(1, Stub::incoming('a'))
            .leg(2, Stub::outgoing('a'))
            .leg(3, Stub::incoming('c'))
            .leg(4, Stub::outgoing('c'));
        for request in [
            query(&photon, 1),
            query(&photon, 2).weight("order", 1..=2),
            query(&four_photons, 1),
            query(&pinned, 1),
            query(&Legs::new(), 2),
        ] {
            let plan = Plan::new(&theory, &request).unwrap();
            assert_eq!(
                plan.compiled
                    .tables
                    .class_members
                    .iter()
                    .map(|m| m.to_vec())
                    .collect::<Vec<_>>(),
                vec![vec!['a', 'c'], vec!['b']]
            );
            let particle = exact(&theory, &request, Classes::None);
            assert!(!particle.is_empty());
            assert_eq!(exact(&theory, &request, Classes::Auto), particle);
        }
    }

    #[test]
    fn original_rule_exclusions_and_active_colours_split_before_quotienting() {
        let base = theory();
        let legs = Legs::new()
            .leg(1, Stub::undirected('b'))
            .leg(2, Stub::undirected('b'));
        let plain = query(&legs, 1);
        for request in [
            plain.clone().class_distinctions([('a', 1)], []),
            plain.clone().class_distinctions([], [(0, 1)]),
            plain.clone().count(Countable::Propagators(vec!['a']), ..=1),
        ] {
            let plan = Plan::new(&base, &request).unwrap();
            assert!(plan.compiled.model.refinement.is_none());
            assert_eq!(
                exact(&base, &request, Classes::Auto),
                exact(&base, &request, Classes::None)
            );
        }
        for mixed in [
            vec![
                Stub::incoming('a'),
                Stub::outgoing('c'),
                Stub::undirected('b'),
            ],
            vec![
                Stub::incoming('a'),
                Stub::outgoing('a'),
                Stub::incoming('c'),
                Stub::outgoing('c'),
            ],
        ] {
            let mut builder = Theory::builder();
            for field in ['a', 'c'] {
                builder = builder.vertex(
                    [
                        Stub::incoming(field),
                        Stub::outgoing(field),
                        Stub::undirected('b'),
                    ],
                    [],
                );
            }
            let unsafe_theory = builder.vertex(mixed, []).build().unwrap();
            let plan = Plan::new(&unsafe_theory, &plain).unwrap();
            assert!(plan.compiled.model.refinement.is_none());
            // The particle-only fallback needs no allocated class-output cache.
            assert!(plan.compiled.tables.class_members.is_empty());
            let particle = exact(&unsafe_theory, &plain, Classes::None);
            assert!(!particle.is_empty());
            assert_eq!(exact(&unsafe_theory, &plain, Classes::Auto), particle);
        }
        let undirected = Theory::builder()
            .vertex(vec![Stub::undirected('a'); 3], [])
            .vertex(vec![Stub::undirected('c'); 3], [])
            .build()
            .unwrap();
        assert!(Plan::new(&undirected, &query(&Legs::new(), 2))
            .unwrap()
            .compiled
            .model
            .refinement
            .is_none());
    }

    #[test]
    fn class_refinement_preserves_thirteen_identical_and_distinct_legs() {
        for identical in [true, false] {
            let mut builder = Theory::builder();
            for field in ['a', 'c'] {
                let mut signature = vec![Stub::incoming(field), Stub::outgoing(field)];
                signature.extend(vec![Stub::undirected('b'); 13]);
                builder = builder.vertex(signature, []);
            }
            let theory = builder.build().unwrap();
            let mut legs = Legs::new()
                .leg(20, Stub::incoming('a'))
                .leg(21, Stub::outgoing('a'));
            for i in 0..13 {
                legs = legs.leg(if identical { 0 } else { i }, Stub::undirected('b'));
            }
            let request = query(&legs, 0).structure(Structure {
                max_vertices: Some(1),
                ..Structure::default()
            });
            assert!(Plan::new(&theory, &request)
                .unwrap()
                .compiled
                .model
                .refinement
                .is_some());
            let particle = exact(&theory, &request, Classes::None);
            assert_eq!(particle.len(), 1);
            assert_eq!(
                particle.values().copied().next(),
                Some(if identical { (1..=13).product() } else { 1 })
            );
            assert_eq!(exact(&theory, &request, Classes::Auto), particle);
        }
        // The same leg kernel also applies when there is no nontrivial flavour class.
        let scalar = Theory::builder()
            .vertex(vec![Stub::undirected('b'); 13], [])
            .build()
            .unwrap();
        for identical in [true, false] {
            let legs = (0..13).fold(Legs::new(), |legs, i| {
                legs.leg(if identical { 0 } else { i }, Stub::undirected('b'))
            });
            let request = query(&legs, 0).structure(Structure {
                max_vertices: Some(1),
                ..Structure::default()
            });
            let particle = exact(&scalar, &request, Classes::None);
            assert_eq!(
                particle.values().copied().collect::<Vec<_>>(),
                [if identical { (1..=13).product() } else { 1 }]
            );
            assert_eq!(exact(&scalar, &request, Classes::Auto), particle);
        }
        // The unrelated-field singleton path must retain the cached 13! node-group order.
        let theory = Theory::builder()
            .vertex([Stub::undirected('b')], [])
            .vertex(vec![Stub::undirected('b'); 13], [])
            .build()
            .unwrap();
        let request = query(&Legs::new(), 0).structure(Structure {
            max_vertices: Some(14),
            zero_flow_edges: true,
            ..Structure::default()
        });
        let particle = exact(&theory, &request, Classes::None);
        let mut orders = particle.values().copied().collect::<Vec<_>>();
        orders.sort();
        assert_eq!(orders, [2, (1..=13).product()]);
        assert_eq!(exact(&theory, &request, Classes::Auto), particle);
    }

    #[test]
    fn anonymous_flavour_pins_and_active_weight_bounds_preserve_complete_maps() {
        let theory = theory();
        let legs = Legs::new()
            .leg(0, Stub::incoming('a'))
            .leg(0, Stub::outgoing('a'))
            .leg(0, Stub::incoming('c'))
            .leg(0, Stub::outgoing('c'));
        let request = query(&legs, 1);
        for bounded in [request.clone(), request.clone().weight("order", 1..=2)] {
            let particle = exact(&theory, &bounded, Classes::None);
            assert!(!particle.is_empty());
            assert_eq!(exact(&theory, &bounded, Classes::Auto), particle);
        }
        // A fine rule's order cannot be hidden by the coarse fibre or a loose partial bound.
        let impossible = request.weight("order", 100..=100);
        for mode in [Classes::None, Classes::Auto] {
            match Plan::new(&theory, &impossible.clone().classes(mode)) {
                Err(GenerationError::Infeasible(report)) => {
                    assert_eq!(
                        report.infeasible.as_deref(),
                        Some("no loop partition satisfies the request")
                    );
                    assert_eq!(
                        report.empty_partitions,
                        vec![(1, "no degree/weight partition satisfies the request".into())]
                    );
                }
                _ => panic!("{mode:?}: impossible order must be rejected during planning"),
            }
        }
    }

    #[test]
    fn class_views_keep_flow_ranked_members_and_owned_refinement_snapshots() {
        let theory = theory();
        let legs = Legs::new()
            .leg(1, Stub::undirected('b'))
            .leg(2, Stub::undirected('b'));
        let request = query(&legs, 1).deferred_hooks(HookMask {
            dressings: true,
            ..HookMask::default()
        });
        let plan = Plan::new(&theory, &request).unwrap();
        let run = Run::default();
        let mut driver = plan.batched(&run, Output::Diagrams).unwrap();
        let mut retained = None;
        let mut expected = BTreeMap::new();
        let mut saw_particle = false;
        loop {
            match driver.poll().unwrap() {
                Event::Batch(requests) => {
                    let replies = requests
                        .iter()
                        .map(|request| request.reply(true))
                        .collect::<Vec<_>>();
                    for request in requests {
                        let ranked = request.into_ranked();
                        let retain = {
                            let view = ranked.dressing().unwrap();
                            let mut retain = None;
                            for (index, edge) in view.edges().enumerate() {
                                match view.field(index) {
                                    Some(Dressed::Class(fields)) => {
                                        assert_eq!(fields, &[0, 2]);
                                        assert!(matches!(
                                            view.flow(index),
                                            Some(Flow::Source | Flow::Sink)
                                        ));
                                        assert_eq!(edge.flow, view.flow(index));
                                        assert!(edge.stub.is_none());
                                        if retained.is_none() {
                                            retain = Some(index);
                                        }
                                    }
                                    Some(Dressed::Particle(_))
                                        if view.edges().all(|edge| edge.stub.is_some()) =>
                                    {
                                        saw_particle = true
                                    }
                                    _ => {}
                                }
                            }
                            retain
                        };
                        if let Some(index) = retain {
                            retained = Some((ranked, index));
                        }
                    }
                    driver.resolve(&replies).unwrap();
                }
                Event::Diagram(diagram) => {
                    assert!(expected
                        .insert(
                            diagram.canonical_form().clone(),
                            diagram.symmetry_factor().unwrap()
                        )
                        .is_none());
                }
                Event::Done => break,
                Event::Yield => {}
                _ => panic!("unexpected output"),
            }
        }
        assert!(saw_particle);
        let (retained, index) = retained.expect("a coarse callback snapshot");
        assert!(matches!(
            retained.dressing().unwrap().field(index),
            Some(Dressed::Class(&[0, 2]))
        ));
        assert_eq!(expected, exact(&theory, &query(&legs, 1), Classes::None));
    }

    #[test]
    fn class_veto_survives_refinement_for_the_existing_heredity_diagnostic() {
        let theory = theory();
        let legs = Legs::new()
            .leg(1, Stub::undirected('b'))
            .leg(2, Stub::undirected('b'));
        let request = query(&legs, 1).prune_dressings(|view| {
            !(0..view.edges().count())
                .any(|edge| matches!(view.field(edge), Some(Dressed::Class(_))))
        });
        let plan = Plan::new(&theory, &request).unwrap();
        for parallel in [false, true] {
            let run = Run::default().parallel(parallel).check_heredity(true);
            assert!(
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| plan.count(&run)))
                    .is_err()
            );
        }
    }

    #[test]
    fn refinement_tiny_work_and_deferred_verdicts_preserve_every_assignment_and_hook() {
        use crate::half_edge::generate::topology::{Context, Search};
        use std::sync::atomic::AtomicBool;
        struct Trace(Mutex<Vec<(Vec<u8>, bool)>>);
        impl Prune for Trace {
            fn keeps(&self, state: State<'_>) -> bool {
                let State::Decoration(d) = state else {
                    unreachable!()
                };
                self.0
                    .lock()
                    .unwrap()
                    .push((d.labels.to_vec(), d.particles));
                true
            }
            fn checks_heredity(&self) -> bool {
                false
            }
        }
        let theory = theory();
        let legs = Legs::new()
            .leg(1, Stub::undirected('b'))
            .leg(2, Stub::undirected('b'));
        let plan = Plan::new(&theory, &query(&legs, 1)).unwrap();
        let model = &plan.compiled.model;
        let limits = &plan.orders[0].1;
        let stop = AtomicBool::new(false);
        let cx = Context {
            model,
            limits,
            stop: &stop,
            abort_check: None,
            filter: None,
        };
        let mut search = Search::new(model, Topology::empty(model, limits), Default::default());
        while search.next(&cx) {
            let t = search.topology();
            let mut reference = None;
            for (defer, quantum) in [(false, None), (false, Some(1)), (true, Some(1))] {
                let trace = Trace(Mutex::new(Vec::new()));
                let mut decorator = Decorator::default();
                decorator.start(model, t, search.generators(), search.node_group_order());
                let mut outputs = Vec::new();
                loop {
                    match decorator.poll(
                        model,
                        t,
                        Some(&trace),
                        defer,
                        &mut quantum.map_or_else(Work::unlimited, Work::bounded),
                    ) {
                        Step::Ready => {
                            let d = decorator.decoration(t);
                            outputs.push((
                                d.labels.to_vec(),
                                d.rules.to_vec(),
                                d.automorphism_count,
                            ));
                            // No nontrivial same-label parallel kernel can exchange two distinct
                            // refinement components under the single diagonal line contract.
                            for &(start, end) in &decorator.bundles {
                                for pair in (start..end).collect::<Vec<_>>().windows(2) {
                                    let a = pair[0];
                                    let b = pair[1];
                                    if decorator.coarse_labels.get(a)
                                        == decorator.coarse_labels.get(b)
                                    {
                                        assert_eq!(d.labels[a], d.labels[b]);
                                    }
                                }
                            }
                        }
                        Step::Pending => {
                            let keep = trace.keeps(State::Decoration(decorator.pending(t)));
                            decorator.resume(model, keep, false);
                        }
                        Step::Yield => {}
                        Step::Done => break,
                    }
                }
                let result = (outputs, trace.0.into_inner().unwrap());
                if let Some(expected) = &reference {
                    assert_eq!(&result, expected);
                } else {
                    reference = Some(result);
                }
            }
        }
    }
    #[test]
    fn low_degree_rule_fibres_preserve_every_coloured_assignment() {
        let mut builder = Theory::builder().weights(["order"]);
        for field in ['a', 'c'] {
            builder = builder.vertex(
                [
                    Stub::incoming(field),
                    Stub::outgoing(field),
                    Stub::undirected('b'),
                ],
                [1],
            );
            for weight in [0, 1] {
                builder = builder.vertex([Stub::incoming(field), Stub::outgoing(field)], [weight]);
            }
        }
        let theory = builder.build().unwrap();
        let legs = Legs::new()
            .leg(1, Stub::undirected('b'))
            .leg(2, Stub::undirected('b'));
        let request = query(&legs, 1)
            .structure(Structure {
                self_loops: true,
                zero_flow_edges: true,
                max_vertices: Some(4),
                ..Structure::default()
            })
            .weight("order", 2..=3);
        let expected = exact(&theory, &request, Classes::None);
        assert!(!expected.is_empty());
        assert_eq!(exact(&theory, &request, Classes::Auto), expected);
    }

    #[test]
    fn cancellation_at_class_and_particle_hooks_stops_before_any_output() {
        use std::sync::{
            atomic::{AtomicBool, AtomicUsize, Ordering},
            Arc,
        };
        for particle in [false, true] {
            let mut builder = Theory::builder();
            for field in ['a', 'c'] {
                builder = builder.vertex(
                    [
                        Stub::incoming(field),
                        Stub::outgoing(field),
                        Stub::undirected('b'),
                    ],
                    [],
                );
            }
            let theory = builder.build().unwrap();
            let legs = Legs::new()
                .leg(1, Stub::undirected('b'))
                .leg(2, Stub::undirected('b'));
            for mode in 0..3 {
                let stop = Arc::new(AtomicBool::new(false));
                let triggered = Arc::new(AtomicUsize::new(0));
                let callback_stop = Arc::clone(&stop);
                let callback_triggered = Arc::clone(&triggered);
                let request = query(&legs, 1).prune_dressings(move |view| {
                    assert!(
                        !callback_stop.load(Ordering::Relaxed),
                        "no callback after cancellation"
                    );
                    let boundary = (0..view.edges().count()).any(|edge| match view.field(edge) {
                        Some(Dressed::Class(_)) => !particle,
                        Some(Dressed::Particle(&('a' | 'c'))) => particle,
                        _ => false,
                    });
                    if boundary {
                        callback_triggered.fetch_add(1, Ordering::Relaxed);
                        callback_stop.store(true, Ordering::Relaxed);
                    }
                    true
                });
                let plan = Plan::new(&theory, &request).unwrap();
                let run = Run::default().abort(&stop);
                match mode {
                    0 => {
                        let counts = plan.count(&run).unwrap();
                        assert_eq!(counts.diagrams, num_bigint::BigUint::from(0u8));
                        assert!(!counts.complete);
                    }
                    1 => assert!(plan.diagrams(&run).unwrap().next().is_none()),
                    _ => {
                        let skeleton_run = Run::default();
                        let skeleton = plan
                            .skeletons(&skeleton_run)
                            .unwrap()
                            .next()
                            .unwrap()
                            .snapshot();
                        assert!(plan.dressings(&skeleton, &run).unwrap().next().is_none());
                    }
                }
                assert_eq!(triggered.load(Ordering::Relaxed), 1);
            }
        }
    }
}
