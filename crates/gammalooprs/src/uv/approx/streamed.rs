//! Memory-bounded forms of the whole-expression maps of the Vakint bridge.
//!
//! The R* outer integrand of a five-loop diagram can be a nested product of
//! sums of almost a gigabyte. Tensor materialization, expansion, metric
//! contraction and the propagator rewrites act on it as a whole, and each of
//! them transiently holds many times its size. Every map here acts on pieces
//! of the expression instead and recombines them into the atom the whole map
//! returns: a normalized Symbolica sum or product is canonical, so adding or
//! multiplying the normalized pieces gives the same atom.
//!
//! Each streamed map falls back to the whole map when a piece introduces a
//! dummy index, because fresh dummies are numbered over the whole expression.
//! Expressions below [`piece_limit`] always take the whole map.
//!
//! Validation only: `GAMMALOOP_STREAM_PIECE_BYTES=<n>` sets the piece limit,
//! and `GAMMALOOP_STREAM_CHECK_BYTES=<n>` also runs the whole map on every
//! streamed input of at most n bytes and reports on stderr whether both
//! results are identical (`STREAM_CHECK stage=... equal=...`).

use std::{collections::HashSet, sync::LazyLock};

use color_eyre::Result;
use spenso::{
    network::{
        parsing::{AtomStructureExt, StrictTensorFilter},
        tags::SPENSO_TAG,
    },
    structure::slot::{SlotMatch, SlotMatcher},
};
use symbolica::atom::{Atom, AtomCore, AtomView};

use crate::utils::GS;

/// Expressions and pieces smaller than this are mapped as a whole.
const PIECE_BYTES: usize = 4 << 20;

static PIECE_LIMIT: LazyLock<usize> = LazyLock::new(|| {
    std::env::var("GAMMALOOP_STREAM_PIECE_BYTES")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(PIECE_BYTES)
});

static CHECK_LIMIT: LazyLock<Option<usize>> = LazyLock::new(|| {
    std::env::var("GAMMALOOP_STREAM_CHECK_BYTES")
        .ok()
        .and_then(|value| value.parse().ok())
});

/// Tests stream every piece, so that the streamed maps are compared with the
/// whole maps on every expression a test reaches.
pub(crate) fn piece_limit() -> usize {
    if cfg!(test) { 1 } else { *PIECE_LIMIT }
}

/// Whether the whole map should also be run on this streamed input to check
/// the streamed result. Always in tests.
pub(crate) fn check_streamed(input: AtomView<'_>) -> bool {
    cfg!(test) || CHECK_LIMIT.is_some_and(|limit| input.get_byte_size() <= limit)
}

/// Report a streamed result against the whole map. Tests fail on a difference.
pub(crate) fn report_check(stage: &str, input: AtomView<'_>, streamed: &Atom, whole: &Atom) {
    let equal = streamed == whole;
    if cfg!(test) {
        assert!(
            equal,
            "streamed {stage} differs from the whole map\ninput: {input}\nstreamed: {streamed}\nwhole: {whole}"
        );
    } else {
        eprintln!(
            "STREAM_CHECK stage={stage} bytes={} equal={equal}",
            input.get_byte_size()
        );
    }
}

/// The labels of the `dummy(i)` indices of an expression, or `None` if one of
/// them is not a single integer.
pub(crate) fn dummy_indices(expression: AtomView<'_>) -> Option<HashSet<i64>> {
    let mut indices = HashSet::new();
    let mut opaque = false;
    expression.visitor(&mut |part| {
        if let AtomView::Fun(function) = part
            && function.get_symbol() == GS.dummyaind
        {
            match (function.get_nargs(), i64::try_from(function.get(0))) {
                (1, Ok(index)) => {
                    indices.insert(index);
                }
                _ => opaque = true,
            }
            return false;
        }
        true
    });
    (!opaque).then_some(indices)
}

/// Whether `output` has a dummy index that `input` lacks.
pub(crate) fn introduces_dummies(input: AtomView<'_>, output: AtomView<'_>) -> bool {
    if !output.contains_symbol(GS.dummyaind) {
        return false;
    }
    match (dummy_indices(input), dummy_indices(output)) {
        (Some(before), Some(after)) => !after.is_subset(&before),
        _ => true,
    }
}

/// Whether tensor materialization leaves `expression` unchanged: it has no
/// trace, chain or scalar marker for the shorthand rewrites, and the parser
/// classifies it as a pure scalar, which it retains verbatim.
pub(crate) fn nothing_to_materialize(expression: AtomView<'_>) -> bool {
    !expression.contains_symbol(SPENSO_TAG.trace)
        && !expression.contains_symbol(SPENSO_TAG.chain)
        && !expression.contains_symbol(SPENSO_TAG.pure_scalar)
        && !expression.is_tensorial(StrictTensorFilter::ContainsReps)
}

/// `whole(expression)` for a map that is additive over summands and
/// multiplicative over index-disjoint factors (tensor materialization).
///
/// Sums are mapped summand by summand. Products are split into groups of
/// factors that share no explicit index (all factors with compact slots
/// outside a `dot` form one group), and each group is mapped alone. Pieces
/// below `limit`, and products that do not split, are mapped by `whole`;
/// pieces for which `skip` holds are kept. Returns `None` if a piece
/// introduces a dummy index.
pub(crate) fn map_multilinear(
    expression: AtomView<'_>,
    limit: usize,
    skip: &dyn Fn(AtomView<'_>) -> bool,
    whole: &dyn Fn(AtomView<'_>) -> Result<Atom>,
) -> Result<Option<Atom>> {
    let leaf = |piece: AtomView<'_>| -> Result<Option<Atom>> {
        if skip(piece) {
            return Ok(Some(piece.to_owned()));
        }
        let mapped = whole(piece)?;
        Ok((!introduces_dummies(piece, mapped.as_view())).then_some(mapped))
    };
    if expression.get_byte_size() < limit {
        return leaf(expression);
    }
    match expression {
        AtomView::Add(sum) => {
            let mut terms = Vec::with_capacity(sum.get_nargs());
            for term in sum.iter() {
                let Some(mapped) = map_multilinear(term, limit, skip, whole)? else {
                    return Ok(None);
                };
                terms.push(mapped);
            }
            Ok(Some(Atom::add_many(terms.iter())))
        }
        AtomView::Mul(product) => {
            let factors = product.iter().collect::<Vec<_>>();
            let groups = index_disjoint_groups(&factors);
            if groups.len() < 2 {
                return leaf(expression);
            }
            let mut mapped = Vec::with_capacity(groups.len());
            for group in groups {
                let owned;
                let piece = if let [factor] = group[..] {
                    factors[factor]
                } else {
                    owned = Atom::mul_many(group.iter().map(|&factor| factors[factor]));
                    owned.as_view()
                };
                let Some(result) = map_multilinear(piece, limit, skip, whole)? else {
                    return Ok(None);
                };
                mapped.push(result);
            }
            Ok(Some(Atom::mul_many(mapped.iter())))
        }
        _ => leaf(expression),
    }
}

/// Partition the factors of a product (by position, in order) into groups
/// that share no explicit index. Factors with a compact or malformed slot
/// outside a `dot` are put in one group, since such a slot can hide an index.
fn index_disjoint_groups(factors: &[AtomView<'_>]) -> Vec<Vec<usize>> {
    let mut matcher = SlotMatcher::default();
    let interfaces = factors
        .iter()
        .map(|factor| factor_interface(*factor, &mut matcher))
        .collect::<Vec<_>>();
    let mut parent = (0..factors.len()).collect::<Vec<_>>();
    fn root(parent: &mut [usize], mut node: usize) -> usize {
        while parent[node] != node {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        node
    }
    for left in 0..factors.len() {
        for right in left + 1..factors.len() {
            let (a, b) = (&interfaces[left], &interfaces[right]);
            let linked =
                (a.opaque && b.opaque) || a.indices.iter().any(|index| b.indices.contains(index));
            if linked {
                let (left_root, right_root) = (root(&mut parent, left), root(&mut parent, right));
                parent[left_root.max(right_root)] = left_root.min(right_root);
            }
        }
    }
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut group_of_root = vec![usize::MAX; factors.len()];
    for factor in 0..factors.len() {
        let factor_root = root(&mut parent, factor);
        if group_of_root[factor_root] == usize::MAX {
            group_of_root[factor_root] = groups.len();
            groups.push(Vec::new());
        }
        groups[group_of_root[factor_root]].push(factor);
    }
    groups
}

struct FactorInterface {
    indices: Vec<Atom>,
    opaque: bool,
}

fn factor_interface(factor: AtomView<'_>, matcher: &mut SlotMatcher) -> FactorInterface {
    let mut interface = FactorInterface {
        indices: Vec::new(),
        opaque: false,
    };
    collect_interface(factor, matcher, false, &mut interface);
    interface
}

fn collect_interface(
    expression: AtomView<'_>,
    matcher: &mut SlotMatcher,
    inside_dot: bool,
    interface: &mut FactorInterface,
) {
    expression.visitor(&mut |part| {
        let AtomView::Fun(function) = part else {
            return true;
        };
        match matcher.classify(part) {
            SlotMatch::Explicit(slot) => {
                interface.indices.push(slot.index().to_owned());
                return false;
            }
            SlotMatch::Opaque => {
                interface.opaque |= !inside_dot;
                return false;
            }
            SlotMatch::Other => {}
        }
        if !inside_dot && function.get_symbol() == SPENSO_TAG.dot {
            for argument in function.iter() {
                collect_interface(argument, matcher, true, interface);
            }
            return false;
        }
        true
    });
}

/// Whether a power of a composite base with an exponent other than a positive
/// integer occurs outside function arguments. `expand` keeps such a power as
/// an opaque factor, which a distributed part could cancel against
/// (`(x+y)*(1/(x+y)+z)`) in a product that `expand` never forms.
fn has_opaque_composite_power(expression: AtomView<'_>) -> bool {
    let mut found = false;
    expression.visitor(&mut |part| {
        if found {
            return false;
        }
        match part {
            AtomView::Fun(_) => false,
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                if matches!(base, AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_))
                    && !i64::try_from(exponent).is_ok_and(|exponent| exponent > 0)
                {
                    found = true;
                }
                !found
            }
            _ => true,
        }
    });
    found
}

/// Pass the expansion of `expression` to `sink` in pieces whose sum is
/// `expression.expand()`. Each piece is an expanded, normalized atom.
///
/// Symbolica expands a product of sums by concatenating every distributed
/// term and normalizing once, so a product whose expansion merges heavily
/// holds many times its result. Here a sum is expanded in chunks of
/// consecutive summands and a product is distributed over parts of its
/// largest sum factor, recursively, down to pieces below `limit`. The
/// expansion of a normalized sum is canonical (sorted monomials, like terms
/// merged with exact coefficients, zeros dropped), so adding the pieces gives
/// the atom `expand` returns. An expression with an opaque composite power
/// is expanded as a whole.
pub(crate) fn expand_in_pieces(
    expression: AtomView<'_>,
    limit: usize,
    sink: &mut dyn FnMut(Atom) -> Result<()>,
) -> Result<()> {
    if expression.get_byte_size() < limit || has_opaque_composite_power(expression) {
        return sink(expression.expand());
    }
    expand_pieces_unchecked(expression, limit, sink)
}

fn expand_pieces_unchecked(
    expression: AtomView<'_>,
    limit: usize,
    sink: &mut dyn FnMut(Atom) -> Result<()>,
) -> Result<()> {
    const PARTS: usize = 16;
    if expression.get_byte_size() < limit {
        return sink(expression.expand());
    }
    match expression {
        AtomView::Add(sum) => {
            let mut chunk = Vec::new();
            let mut chunk_bytes = 0;
            let flush = |chunk: &mut Vec<AtomView<'_>>,
                         sink: &mut dyn FnMut(Atom) -> Result<()>|
             -> Result<()> {
                match chunk[..] {
                    [] => Ok(()),
                    [term] => expand_pieces_unchecked(term, limit, sink),
                    _ => sink(Atom::add_many(chunk.iter().copied()).expand()),
                }
                .map(|()| chunk.clear())
            };
            for term in sum.iter() {
                let bytes = term.get_byte_size();
                if chunk_bytes + bytes >= limit {
                    flush(&mut chunk, sink)?;
                    chunk_bytes = 0;
                }
                chunk.push(term);
                chunk_bytes += bytes;
            }
            flush(&mut chunk, sink)
        }
        AtomView::Mul(product) => {
            let factors = product.iter().collect::<Vec<_>>();
            let Some((split, sum)) = factors
                .iter()
                .enumerate()
                .filter_map(|(index, factor)| match factor {
                    AtomView::Add(sum) => Some((index, *sum)),
                    _ => None,
                })
                .max_by_key(|(_, sum)| sum.as_view().get_byte_size())
            else {
                return sink(expression.expand());
            };
            let terms = sum.iter().collect::<Vec<_>>();
            for part in terms.chunks(terms.len().div_ceil(PARTS)) {
                let part = Atom::add_many(part.iter().copied());
                let distributed =
                    Atom::mul_many(factors.iter().enumerate().map(|(index, factor)| {
                        if index == split {
                            part.as_view()
                        } else {
                            *factor
                        }
                    }));
                drop(part);
                expand_pieces_unchecked(distributed.as_view(), limit, sink)?;
            }
            Ok(())
        }
        _ => sink(expression.expand()),
    }
}

/// Visit the summands of `expression` (itself, if it is not a sum) in chunks
/// of consecutive summands of about `limit` bytes. A chunk of several
/// summands is passed as their normalized sum.
pub(crate) fn for_each_summand_chunk(
    expression: AtomView<'_>,
    limit: usize,
    visit: &mut dyn FnMut(AtomView<'_>) -> Result<()>,
) -> Result<()> {
    let AtomView::Add(sum) = expression else {
        return visit(expression);
    };
    if expression.get_byte_size() < limit {
        return visit(expression);
    }
    let mut chunk = Vec::new();
    let mut chunk_bytes = 0;
    let flush = |chunk: &mut Vec<AtomView<'_>>,
                 visit: &mut dyn FnMut(AtomView<'_>) -> Result<()>|
     -> Result<()> {
        match chunk[..] {
            [] => Ok(()),
            [term] => visit(term),
            _ => visit(Atom::add_many(chunk.iter().copied()).as_view()),
        }
        .map(|()| chunk.clear())
    };
    for term in sum.iter() {
        let bytes = term.get_byte_size();
        if chunk_bytes + bytes >= limit {
            flush(&mut chunk, visit)?;
            chunk_bytes = 0;
        }
        chunk.push(term);
        chunk_bytes += bytes;
    }
    flush(&mut chunk, visit)
}

/// A running sum of normalized atoms. Pending summands are merged into the
/// accumulated sum whenever they outgrow it, so that the unmerged part never
/// holds much more than the result. A normalized sum is canonical, so the
/// result is `Atom::add_many` of all pushed atoms.
#[derive(Default)]
pub(crate) struct SumAccumulator {
    accumulated: Atom,
    pending: Vec<Atom>,
    pending_bytes: usize,
}

impl SumAccumulator {
    pub(crate) fn push(&mut self, summand: Atom) {
        self.pending_bytes += summand.as_view().get_byte_size();
        self.pending.push(summand);
        if self.pending_bytes > self.accumulated.as_view().get_byte_size().max(64 << 20) {
            self.merge();
        }
    }

    fn merge(&mut self) {
        let merged = Atom::add_many(std::iter::once(&self.accumulated).chain(self.pending.iter()));
        self.pending.clear();
        self.pending_bytes = 0;
        self.accumulated = merged;
    }

    pub(crate) fn finish(mut self) -> Atom {
        if !self.pending.is_empty() {
            self.merge();
        }
        self.accumulated
    }
}

/// Names listed in `GAMMALOOP_ALGEBRA_LEGACY` (comma separated, or `all`)
/// select the previous implementation of that algebra fast path, so that a
/// regression can be bisected without a rebuild.
pub(crate) fn legacy(name: &str) -> bool {
    static LEGACY: LazyLock<Vec<String>> =
        LazyLock::new(|| switch_list("GAMMALOOP_ALGEBRA_LEGACY"));
    LEGACY.iter().any(|entry| entry == name || entry == "all")
}

/// Names listed in `GAMMALOOP_ALGEBRA_CHECK` (or `all`) also run the previous
/// implementation and report on stderr whether both results are identical
/// (`STREAM_CHECK stage=<name> ... equal=...`). Tests always check.
pub(crate) fn check(name: &str) -> bool {
    static CHECK: LazyLock<Vec<String>> = LazyLock::new(|| switch_list("GAMMALOOP_ALGEBRA_CHECK"));
    cfg!(test) || CHECK.iter().any(|entry| entry == name || entry == "all")
}

fn switch_list(variable: &str) -> Vec<String> {
    std::env::var(variable)
        .map(|value| {
            value
                .split(',')
                .map(|entry| entry.trim().to_owned())
                .collect()
        })
        .unwrap_or_default()
}
