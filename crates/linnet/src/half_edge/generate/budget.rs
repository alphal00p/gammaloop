//! Budgets: sums of signature weights (coupling orders, counterterm orders, any integer tag) and
//! counts of vertices, propagators, chords and closed lines that the generated graphs keep within
//! ranges.
//!
//! Each budget is a bound with an inclusive range. Stage 2 keeps, per bound, the interval of the
//! values that the completions of a partial field assignment can still take, one row of intervals
//! per labelled edge, and drops the assignment as soon as an interval misses its range:
//!
//! - A weight is summed over the vertices. A partial vertex is a state of the vertex automaton, and
//!   the least and the most weight of the signatures containing the state bound the weight of every
//!   vertex completing it; the interval is the sum of these bounds. Labelling an edge moves its two
//!   ends only, and the bounds only narrow as states grow, so pruning is hereditary and costs a few
//!   additions per label. The vertex rules are checked exactly once they are chosen.
//! - Vertices of some signatures are a weight of one on those signatures.
//! - Propagators and chords of some fields: every labelled edge either counts or can no longer
//!   count, so the count so far and the count plus the edges still to label bound the value. The
//!   topology decides which edges are chords (no bridges).
//! - Closed lines of some fields: a line closes with its last edge and stays closed, so the lines
//!   closed so far bound the value from below as long as no vertex can take three stubs of the
//!   fields. If a signature can, lines may branch, and they are counted on complete labellings.
//!
//! Before any search, [`Budgets::partition_exists`] proves a loop order empty when no multiset of
//! signatures has the vertex degrees that the loop order needs (`sum (degree - 2) = 2 L - 2 +
//! legs`) and weights within the ranges. It follows the reachable sums of degrees and weights
//! rather than the multisets themselves (the partitions of QGRAF and grcc), which multiply with the
//! weight classes: on SMEFTsim with NP at most two, 2 -> 3 at three loops has 46,697 partitions but
//! 2,533 sums, while the Standard Model requests measured have at most 30. The search carries the
//! budgets as the running bounds above instead of running once per partition.

use std::ops::{Bound, RangeBounds};

use ahash::AHashSet;

use super::{
    model::{Label, LoopLabel},
    GenerationError,
};

/// A count over a graph, which [`super::Query::count`] keeps within a range.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Countable<F> {
    /// Vertices whose signature index is listed.
    Vertices(Vec<usize>),
    /// Paired edges (propagators) carrying one of the fields; legs are not propagators.
    Propagators(Vec<F>),
    /// Propagators carrying one of the fields that lie on a cycle, that is, are no bridges.
    Chords(Vec<F>),
    /// Closed lines of the fields: cycles of edges carrying them through vertices without other
    /// stubs of them, so a line through a leg of one of the fields is open. Lines are defined
    /// while every vertex has at most two stubs of the fields; a graph with a vertex of three or
    /// more is not constrained.
    ClosedLines(Vec<F>),
}

/// The inclusive range `[min, max]` of `range`, with values converted by `value`; empty
/// (`min > max`) if `range` is.
pub(crate) fn inclusive<T: Copy>(
    range: &impl RangeBounds<T>,
    value: impl Fn(T) -> i64,
) -> [i64; 2] {
    let min = match range.start_bound() {
        Bound::Included(&t) => value(t),
        Bound::Excluded(&t) => value(t).saturating_add(1),
        Bound::Unbounded => i64::MIN,
    };
    let max = match range.end_bound() {
        Bound::Included(&t) => value(t),
        Bound::Excluded(&t) => value(t).saturating_sub(1),
        Bound::Unbounded => i64::MAX,
    };
    [min, max]
}

/// The budgets of a request with fields and signatures as dense indices, as [`Budgets::new`]
/// compiles them.
#[derive(Clone, Debug, Default)]
pub(crate) struct Request {
    /// Per weight: the weight of every signature and the range of the sum.
    pub(crate) weights: Vec<(Vec<i32>, [i64; 2])>,
    /// Per edge count: whether it counts chords only, its fields (a bit per field id) and range.
    pub(crate) edges: Vec<(bool, u128, [i64; 2])>,
    /// Per line count: its fields and range.
    pub(crate) lines: Vec<(u128, [i64; 2])>,
}

impl Request {
    /// The weights of all signatures, `weights.len()` per signature.
    pub(crate) fn signature_weights(&self, signatures: usize) -> Vec<i32> {
        (0..signatures)
            .flat_map(|s| self.weights.iter().map(move |(weights, _)| weights[s]))
            .collect()
    }
}

/// The budgets of a request compiled against its model: what stage 2 bounds and checks.
#[derive(Debug, Default)]
pub(crate) struct Budgets {
    /// Per bound, its range: the weights, then the edge counts, then the line counts.
    ranges: Vec<[i64; 2]>,
    pub(crate) n_weights: usize,
    /// Per signature, its weights.
    weights: Vec<i32>,
    /// Per edge count: whether it counts chords only.
    chords: Vec<bool>,
    /// Per line count: whether some signature has three or more of its stubs.
    branching: Vec<bool>,
    /// Per paired label and per self-loop label: a bit per edge count, then per line count, whose
    /// fields contain its field.
    paired: Vec<u64>,
    looped: Vec<u64>,
    /// Per leg class: a bit per line count whose fields contain its field.
    legs: Vec<u64>,
    /// The distinct `[degree - 2, weights..]` of the signatures.
    classes: Vec<Vec<i32>>,
}

/// Most edge and line counts: a label marks those it adds to in a `u64`.
const MAX_COUNTS: usize = u64::BITS as usize;

/// Most sums of degrees and weights that [`Budgets::partition_exists`] follows before it gives up.
const MAX_SUMS: usize = 1 << 16;

impl Budgets {
    /// `weights` are those of [`Request::signature_weights`]; `signatures` and leg classes
    /// `(stub, multiplicity)` are given as stub types `3 * field + direction`.
    pub(crate) fn new(
        request: &Request,
        weights: Vec<i32>,
        signatures: &[Vec<u8>],
        labels: &[Label],
        loop_labels: &[LoopLabel],
        legs: &[(u8, u8)],
    ) -> Result<Self, GenerationError> {
        let n_counts = request.edges.len() + request.lines.len();
        if n_counts > MAX_COUNTS {
            return Err(GenerationError::TooManyCounts);
        }
        let n_weights = request.weights.len();
        let has = |fields: u128, field: u8| fields >> field & 1 == 1;
        let counted = |field: u8| -> u64 {
            let edges = request.edges.iter().map(|&(_, fields, _)| fields);
            let lines = request.lines.iter().map(|&(fields, _)| fields);
            edges
                .chain(lines)
                .enumerate()
                .filter(|&(_, fields)| has(fields, field))
                .fold(0, |bits, (c, _)| bits | 1 << c)
        };
        // 64 edge counts leave no bit to line counts
        let line_bits = |field: u8| {
            counted(field)
                .checked_shr(request.edges.len() as u32)
                .unwrap_or(0)
        };
        let mut classes: Vec<Vec<i32>> = match n_weights {
            0 => vec![],
            k => signatures
                .iter()
                .zip(weights.chunks(k))
                .map(|(signature, w)| [&[signature.len() as i32 - 2], w].concat())
                .collect(),
        };
        classes.sort();
        classes.dedup();
        Ok(Budgets {
            ranges: (request.weights.iter().map(|w| w.1))
                .chain(request.edges.iter().map(|e| e.2))
                .chain(request.lines.iter().map(|l| l.1))
                .collect(),
            n_weights,
            weights,
            chords: request.edges.iter().map(|e| e.0).collect(),
            branching: (request.lines.iter())
                .map(|&(fields, _)| {
                    signatures.iter().any(|signature| {
                        signature.iter().filter(|&&s| has(fields, s / 3)).count() > 2
                    })
                })
                .collect(),
            paired: labels.iter().map(|label| counted(label.field)).collect(),
            looped: loop_labels
                .iter()
                .map(|label| counted(label.field))
                .collect(),
            legs: legs.iter().map(|&(stub, _)| line_bits(stub / 3)).collect(),
            classes,
        })
    }

    /// Number of bounds.
    #[inline(always)]
    pub(crate) fn len(&self) -> usize {
        self.ranges.len()
    }

    /// The request has no budgets.
    #[inline(always)]
    pub(crate) fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }

    pub(crate) fn n_edge_counts(&self) -> usize {
        self.chords.len()
    }

    pub(crate) fn n_line_counts(&self) -> usize {
        self.branching.len()
    }

    /// Some range is empty, so no graph meets the budgets.
    pub(crate) fn unsatisfiable(&self) -> bool {
        self.ranges.iter().any(|&[min, max]| min > max)
    }

    /// Some edge count counts chords only.
    pub(crate) fn counts_chords(&self) -> bool {
        self.chords.contains(&true)
    }

    /// Edge count `c` counts chords only.
    #[inline(always)]
    pub(crate) fn chords_only(&self, c: usize) -> bool {
        self.chords[c]
    }

    /// Lines of line count `j` may branch: some signature has three or more of its stubs.
    #[inline(always)]
    pub(crate) fn branching(&self, j: usize) -> bool {
        self.branching[j]
    }

    /// The edge counts, then the line counts, that a paired label (a self-loop label if
    /// `self_loop`) adds to, a bit each.
    #[inline(always)]
    pub(crate) fn counted(&self, label: usize, self_loop: bool) -> u64 {
        if self_loop {
            self.looped[label]
        } else {
            self.paired[label]
        }
    }

    /// The line counts that the legs of class `c` carry stubs of, a bit each.
    pub(crate) fn leg_lines(&self, c: usize) -> u64 {
        self.legs[c]
    }

    /// Every interval of `row`, one per bound, meets its range.
    #[inline(always)]
    pub(crate) fn admit(&self, row: &[[i64; 2]]) -> bool {
        row.iter()
            .zip(&self.ranges)
            .all(|(&[low, high], &[min, max])| low <= max && high >= min)
    }

    /// The weights of the vertex rules `rules` sum into their ranges.
    #[inline]
    pub(crate) fn rules_admit(&self, rules: &[u32]) -> bool {
        let k = self.n_weights;
        (0..k).all(|w| {
            let sum: i64 = rules
                .iter()
                .map(|&r| self.weights[r as usize * k + w] as i64)
                .sum();
            let [min, max] = self.ranges[w];
            (min..=max).contains(&sum)
        })
    }

    /// Some multiset of between one and `max_vertices` signatures has `sum (degree - 2) =
    /// excess` and weights summing into their ranges: a necessary condition for a graph of the
    /// loop order with that `excess`. Sums that only grow (or only shrink) are dropped once past
    /// their range, so the sums followed stay few on coupling orders; past [`MAX_SUMS`] of them the
    /// test gives up and admits the loop order.
    pub(crate) fn partition_exists(&self, excess: i32, max_vertices: usize) -> bool {
        if self.n_weights == 0 {
            return true;
        }
        // the sum of `degree - 2` is one more sum, whose range is the excess
        let ranges: Vec<[i64; 2]> = std::iter::once([excess as i64; 2])
            .chain(self.ranges[..self.n_weights].iter().copied())
            .collect();
        let monotone =
            |d: usize, sign: fn(&i32) -> bool| self.classes.iter().all(|class| sign(&class[d]));
        let grows: Vec<bool> = (0..ranges.len())
            .map(|d| monotone(d, |&x| x >= 0))
            .collect();
        let shrinks: Vec<bool> = (0..ranges.len())
            .map(|d| monotone(d, |&x| x <= 0))
            .collect();
        let dead = |sums: &[i64]| {
            (0..ranges.len()).any(|d| {
                let [min, max] = ranges[d];
                (grows[d] && sums[d] > max) || (shrinks[d] && sums[d] < min)
            })
        };
        let done = |sums: &[i64]| {
            sums.iter()
                .zip(&ranges)
                .all(|(sum, [min, max])| (min..=max).contains(&sum))
        };
        // a sum reached with fewer vertices has every continuation of a later one
        let mut seen: AHashSet<Vec<i64>> = AHashSet::new();
        let mut frontier = vec![vec![0i64; ranges.len()]];
        for _ in 0..max_vertices {
            let mut next = vec![];
            for sums in &frontier {
                for class in &self.classes {
                    let sums: Vec<i64> =
                        sums.iter().zip(class).map(|(s, &c)| s + c as i64).collect();
                    if done(&sums) {
                        return true;
                    }
                    if !dead(&sums) && seen.insert(sums.clone()) {
                        next.push(sums);
                    }
                }
            }
            if seen.len() > MAX_SUMS {
                return true;
            }
            frontier = next;
        }
        false
    }
}

/// The lines of one closed-line count on a topology: a union-find over the vertices whose sets
/// are joined by the labelled edges carrying its fields, with a log to undo them in reverse order.
#[derive(Debug, Default)]
pub(crate) struct Lines {
    parent: Vec<u8>,
    /// Per root: the vertices and the edges of its component.
    vertices: Vec<u16>,
    edges: Vec<u16>,
    /// Per vertex: its stubs of the fields so far, legs included.
    stubs: Vec<u8>,
    /// Vertices with three or more stubs.
    pub(crate) branching: u32,
    /// Components with as many edges as vertices: cycles, while no vertex branches.
    pub(crate) closed: u32,
    /// Per added edge: the root that the edge joined (twice), or the root merged and the root
    /// merged into.
    log: Vec<[u8; 2]>,
}

impl Lines {
    /// No edge labelled yet; vertex `v` has `legs[v]` legs of the fields.
    pub(crate) fn reset(&mut self, legs: impl ExactSizeIterator<Item = u8>) {
        let n = legs.len();
        self.parent.clear();
        self.parent.extend(0..n as u8);
        self.vertices.clear();
        self.vertices.resize(n, 1);
        self.edges.clear();
        self.edges.resize(n, 0);
        self.stubs.clear();
        self.stubs.extend(legs);
        self.branching = self.stubs.iter().filter(|&&s| s > 2).count() as u32;
        self.closed = 0;
        self.log.clear();
    }

    /// The interval of a line count with `closed` closed components and `branching` branching
    /// vertices: exact once every edge is labelled (any value with a branching vertex, where
    /// lines are not defined), and from `closed` on while labelling, unless lines may branch.
    #[inline(always)]
    pub(crate) fn interval(
        closed: u32,
        branching: u32,
        complete: bool,
        may_branch: bool,
    ) -> [i64; 2] {
        match (complete, branching > 0, may_branch) {
            (true, false, _) => [closed as i64; 2],
            (true, true, _) | (false, _, true) => [i64::MIN, i64::MAX],
            (false, _, false) => [closed as i64, i64::MAX],
        }
    }

    fn find(&self, mut v: usize) -> usize {
        while self.parent[v] as usize != v {
            v = self.parent[v] as usize;
        }
        v
    }

    /// The closed components and the branching vertices once edge `(u, v)` is added.
    #[inline(always)]
    pub(crate) fn after(&self, u: usize, v: usize) -> (u32, u32) {
        let branches = |stubs: u8, added: u8| (stubs <= 2 && stubs + added > 2) as u32;
        let branching = self.branching
            + if u == v {
                branches(self.stubs[u], 2)
            } else {
                branches(self.stubs[u], 1) + branches(self.stubs[v], 1)
            };
        let (ru, rv) = (self.find(u), self.find(v));
        let (vertices, edges) = if ru == rv {
            (self.vertices[ru], self.edges[ru] + 1)
        } else {
            (
                self.vertices[ru] + self.vertices[rv],
                self.edges[ru] + self.edges[rv] + 1,
            )
        };
        (self.closed + (edges == vertices) as u32, branching)
    }

    /// Add edge `(u, v)`.
    pub(crate) fn add(&mut self, u: usize, v: usize) {
        (self.closed, self.branching) = self.after(u, v);
        self.stubs[u] += 1;
        self.stubs[v] += 1;
        let (ru, rv) = (self.find(u), self.find(v));
        if ru == rv {
            self.edges[ru] += 1;
            self.log.push([ru as u8; 2]);
            return;
        }
        let (child, root) = if self.vertices[ru] < self.vertices[rv] {
            (ru, rv)
        } else {
            (rv, ru)
        };
        self.parent[child] = root as u8;
        self.vertices[root] += self.vertices[child];
        self.edges[root] += self.edges[child] + 1;
        self.log.push([child as u8, root as u8]);
    }

    /// Undo the last added edge, `(u, v)`.
    pub(crate) fn undo(&mut self, u: usize, v: usize) {
        let [child, root] = self.log.pop().expect("an added edge").map(usize::from);
        self.closed -= (self.edges[root] == self.vertices[root]) as u32;
        if child == root {
            self.edges[root] -= 1;
        } else {
            self.parent[child] = child as u8;
            self.vertices[root] -= self.vertices[child];
            self.edges[root] -= self.edges[child] + 1;
        }
        self.stubs[u] -= 1;
        self.stubs[v] -= 1;
        let branched = |stubs: u8, removed: u8| (stubs > 2 && stubs - removed <= 2) as u32;
        // the stubs are already decremented: compare with the counts before
        self.branching -= if u == v {
            branched(self.stubs[u] + 2, 2)
        } else {
            branched(self.stubs[u] + 1, 1) + branched(self.stubs[v] + 1, 1)
        };
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::{inclusive, Budgets, Lines, Request};

    #[test]
    fn ranges_are_inclusive() {
        assert_eq!(inclusive(&(2..5), i64::from), [2, 4]);
        assert_eq!(inclusive(&(..=3), i64::from), [i64::MIN, 3]);
        assert_eq!(inclusive(&(1..), i64::from), [1, i64::MAX]);
        assert_eq!(inclusive(&(3..3), i64::from), [3, 2]);
        let n = |n: usize| i64::try_from(n).unwrap_or(i64::MAX);
        assert_eq!(inclusive(&(0..=usize::MAX), n), [0, i64::MAX]);
    }

    /// Budgets with weights only, on signatures of the given degrees.
    fn weighted(degrees: &[usize], weights: &[(Vec<i32>, [i64; 2])]) -> Budgets {
        let request = Request {
            weights: weights.to_vec(),
            ..Request::default()
        };
        let signatures: Vec<Vec<u8>> = degrees.iter().map(|&d| vec![0; d]).collect();
        let flat = request.signature_weights(signatures.len());
        Budgets::new(&request, flat, &signatures, &[], &[], &[]).unwrap()
    }

    #[test]
    fn partitions_meet_the_euler_relation_and_the_weights() {
        // electroweak-like: cubic vertices of QED order 0 or 1, quartic ones of order 0 or 2
        let ew = |min, max| weighted(&[3, 3, 4, 4], &[(vec![0, 1, 0, 2], [min, max])]);
        // a 2 -> 2 tree (excess 2): two cubic vertices or one quartic one
        assert!(ew(2, 2).partition_exists(2, 2));
        assert!(ew(0, 0).partition_exists(2, 2));
        assert!(!ew(3, 3).partition_exists(2, 2));
        // one loop (excess 4) reaches order 3 only with a cubic vertex of order one, so with
        // three vertices or more
        assert!(!ew(3, 3).partition_exists(4, 2));
        assert!(ew(3, 3).partition_exists(4, 3));
        assert!(!ew(5, 9).partition_exists(4, 4));
        // two orders that both count every vertex
        let both = weighted(&[3, 4], &[(vec![1, 2], [4, 4]), (vec![1, 1], [3, 3])]);
        assert!(both.partition_exists(4, 4));
        assert!(!both.partition_exists(2, 4));
        // signed weights and vertices of degree one, whose sums are not monotone
        let signed = weighted(&[1, 3], &[(vec![-1, 1], [0, 0])]);
        assert!(signed.partition_exists(0, 2));
        assert!(!signed.partition_exists(-1, 1));
        assert!(!signed.partition_exists(1, 1));
        // without weights every loop order may have graphs
        assert!(weighted(&[3], &[]).partition_exists(5, 1));
    }

    /// Closed components (as many edges as vertices) and branching vertices (three or more
    /// stubs) of the edges `edges` on vertices with `legs` legs, counted directly.
    fn count_lines(legs: &[u8], edges: &[(usize, usize)]) -> (u32, u32) {
        let n = legs.len();
        let mut stubs = legs.to_vec();
        let mut component: Vec<usize> = (0..n).collect();
        let find = |c: &[usize], mut v: usize| {
            while c[v] != v {
                v = c[v];
            }
            v
        };
        for &(u, v) in edges {
            stubs[u] += 1;
            stubs[v] += 1;
            let (a, b) = (find(&component, u), find(&component, v));
            component[a] = b;
        }
        let closed = (0..n)
            .filter(|&r| find(&component, r) == r)
            .filter(|&r| {
                let in_component = |v: &usize| find(&component, *v) == r;
                let vertices = (0..n).filter(in_component).count();
                let lines = edges.iter().filter(|e| in_component(&e.0)).count();
                lines == vertices
            })
            .count();
        (
            closed as u32,
            stubs.iter().filter(|&&s| s > 2).count() as u32,
        )
    }

    proptest! {
        /// Adding edges one at a time, as `after` predicts them, counts the branching vertices and,
        /// while none branches, the closed lines as counting them directly does; undoing the edges
        /// in reverse order restores every earlier count.
        #[test]
        fn lines_follow_added_and_undone_edges(
            legs in proptest::collection::vec(0u8..2, 1..7),
            ends in proptest::collection::vec((0usize..7, 0usize..7), 0..10),
        ) {
            let n = legs.len();
            let edges: Vec<_> = ends.into_iter().map(|(u, v)| (u % n, v % n)).collect();
            let mut lines = Lines::default();
            lines.reset(legs.iter().copied());
            prop_assert_eq!((lines.closed, lines.branching), count_lines(&legs, &[]));
            let mut counts = vec![(lines.closed, lines.branching)];
            for (i, &(u, v)) in edges.iter().enumerate() {
                let predicted = lines.after(u, v);
                lines.add(u, v);
                prop_assert_eq!(predicted, (lines.closed, lines.branching));
                let (closed, branching) = count_lines(&legs, &edges[..=i]);
                prop_assert_eq!(lines.branching, branching);
                if branching == 0 {
                    prop_assert_eq!(lines.closed, closed);
                }
                counts.push(predicted);
            }
            for i in (0..edges.len()).rev() {
                lines.undo(edges[i].0, edges[i].1);
                prop_assert_eq!((lines.closed, lines.branching), counts[i]);
            }
        }
    }

    #[test]
    fn line_intervals_are_exact_on_complete_labellings() {
        assert_eq!(Lines::interval(2, 0, true, true), [2, 2]);
        assert_eq!(Lines::interval(2, 1, true, true), [i64::MIN, i64::MAX]);
        assert_eq!(Lines::interval(2, 0, false, false), [2, i64::MAX]);
        assert_eq!(Lines::interval(2, 0, false, true), [i64::MIN, i64::MAX]);
    }
}
