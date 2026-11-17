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
//!
//! `GAMMALOOP_STREAM_BY_ESTIMATE=1` (opt-in, default off): the Vakint bridge
//! expands in pieces when the estimated size of the expansion
//! ([`expansion_estimate`], bytes before like terms merge) reaches the piece
//! limit, not only when its input does, and [`expand_in_pieces`] sizes its
//! pieces by that estimate (see [`stream_by_estimate`]). A small input with a
//! large expansion (a product of sums) is then distributed in pieces too. The
//! pieces add up to the same atom, so results do not change.

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

#[cfg(test)]
thread_local! {
    static TEST_PIECE_LIMIT: std::cell::Cell<Option<usize>> = const { std::cell::Cell::new(None) };
}

static BY_ESTIMATE: LazyLock<bool> =
    LazyLock::new(|| std::env::var("GAMMALOOP_STREAM_BY_ESTIMATE").as_deref() == Ok("1"));

/// `GAMMALOOP_STREAM_BY_ESTIMATE=1`: decide and size the pieces of the
/// to_vakint expansion by the estimated size of the expansion (see the module
/// documentation). Default off: by the input bytes, as before.
pub(crate) fn stream_by_estimate() -> bool {
    *BY_ESTIMATE
}

/// The size by which expansions are cut into pieces: the bytes of
/// `expression`, or with `by_estimate` the larger of those and the estimated
/// bytes of its expansion ([`expansion_estimate`]).
pub(crate) fn expansion_size(expression: AtomView<'_>, by_estimate: bool) -> usize {
    let bytes = expression.get_byte_size();
    if !by_estimate {
        return bytes;
    }
    let estimate = expansion_estimate(expression).1;
    if estimate >= usize::MAX as f64 {
        usize::MAX
    } else {
        bytes.max(estimate as usize)
    }
}

/// Tests only: run `f` on this thread with the production piece limit, for
/// tests that compare results rather than the streamed maps themselves.
#[cfg(test)]
pub(crate) fn with_production_piece_limit<T>(f: impl FnOnce() -> T) -> T {
    with_piece_limit(*PIECE_LIMIT, f)
}

/// Tests only: run `f` on this thread with the given piece limit.
#[cfg(test)]
pub(crate) fn with_piece_limit<T>(limit: usize, f: impl FnOnce() -> T) -> T {
    let previous = TEST_PIECE_LIMIT.with(|cell| cell.replace(Some(limit)));
    let result = f();
    TEST_PIECE_LIMIT.with(|cell| cell.set(previous));
    result
}

/// Tests stream every piece, so that the streamed maps are compared with the
/// whole maps on every expression a test reaches.
pub(crate) fn piece_limit() -> usize {
    #[cfg(test)]
    if let Some(limit) = TEST_PIECE_LIMIT.with(|limit| limit.get()) {
        return limit;
    }
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

/// The size of `expression.expand()` before like terms merge, as
/// `(summands, bytes)`: a sum adds the sizes of its terms; a product
/// multiplies the summand counts of its factors (every summand of the
/// expansion concatenates one summand of each factor, so its bytes are the
/// sum of the factors' average summand bytes); a positive integer power of a
/// sum or product is that many factors. Functions and other powers are not
/// distributed and count as one summand of their own bytes. Merging like
/// terms only shrinks the expansion, so this bounds what `expand` holds
/// before it normalizes. In floating point: it never overflows.
pub(crate) fn expansion_estimate(expression: AtomView<'_>) -> (f64, f64) {
    match expression {
        AtomView::Add(sum) => sum
            .iter()
            .map(expansion_estimate)
            .fold((0.0, 0.0), |(terms, bytes), (more_terms, more_bytes)| {
                (terms + more_terms, bytes + more_bytes)
            }),
        AtomView::Mul(product) => {
            let (mut terms, mut summand_bytes) = (1.0, 0.0);
            for factor in product.iter() {
                let (factor_terms, factor_bytes) = expansion_estimate(factor);
                if factor_terms == 0.0 {
                    return (0.0, 0.0);
                }
                terms *= factor_terms;
                summand_bytes += factor_bytes / factor_terms;
            }
            (terms, terms * summand_bytes)
        }
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            match (base, i64::try_from(exponent)) {
                (AtomView::Add(_) | AtomView::Mul(_), Ok(n)) if n > 0 => {
                    let (base_terms, base_bytes) = expansion_estimate(base);
                    if base_terms == 0.0 {
                        return (0.0, 0.0);
                    }
                    let terms = base_terms.powf(n as f64);
                    (terms, terms * n as f64 * base_bytes / base_terms)
                }
                _ => (1.0, expression.get_byte_size() as f64),
            }
        }
        _ => (1.0, expression.get_byte_size() as f64),
    }
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
/// is expanded as a whole. Sizes are bytes, or with
/// `GAMMALOOP_STREAM_BY_ESTIMATE=1` estimated expansion bytes
/// ([`expansion_size`]).
pub(crate) fn expand_in_pieces(
    expression: AtomView<'_>,
    limit: usize,
    sink: &mut dyn FnMut(Atom) -> Result<()>,
) -> Result<()> {
    expand_in_pieces_by(expression, limit, stream_by_estimate(), sink)
}

/// [`expand_in_pieces`] with an explicit size measure: `by_estimate` cuts by
/// [`expansion_size`] instead of input bytes.
pub(crate) fn expand_in_pieces_by(
    expression: AtomView<'_>,
    limit: usize,
    by_estimate: bool,
    sink: &mut dyn FnMut(Atom) -> Result<()>,
) -> Result<()> {
    if expansion_size(expression, by_estimate) < limit {
        return sink(expression.expand());
    }
    let blocked = if legacy("opaque_pieces") {
        has_opaque_composite_power(expression)
    } else {
        opaque_power_blocks_pieces(expression)
    };
    if blocked {
        return sink(expression.expand());
    }
    expand_pieces_unchecked(expression, limit, by_estimate, sink)
}

/// Whether distributing `expression` in pieces could cancel an opaque
/// composite power (see [`has_opaque_composite_power`]) that `expand` keeps.
///
/// The pieces multiply chunks of consecutive summands of a sum (or the sum
/// itself) with the other factors of its product; a normalized product
/// merges only factors with the same base, and the only merge that `expand`
/// never makes is a chunk `B` with an opaque power `B^e` (`e` not a positive
/// integer). Merges of a chunk with a positive integer power give the same
/// expansion. A chunk is a sum of terms of one sum of the expression, so it
/// can equal `B` only if some sum of the expression holds every term of
/// `B`. So the pieces are blocked iff an opaque power has a base that is not
/// a sum, or a sum outside function arguments and outside the opaque bases
/// holds every term of an opaque base. Without opaque powers nothing is
/// blocked. `GAMMALOOP_ALGEBRA_LEGACY=opaque_pieces` blocks every
/// expression with an opaque power, as before.
fn opaque_power_blocks_pieces(expression: AtomView<'_>) -> bool {
    let mut bases = Vec::<Vec<&[u8]>>::new();
    let mut blocked = false;
    expression.visitor(&mut |part| {
        if blocked {
            return false;
        }
        match part {
            AtomView::Fun(_) => false,
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                if matches!(base, AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_))
                    && !i64::try_from(exponent).is_ok_and(|exponent| exponent > 0)
                {
                    match base {
                        AtomView::Add(sum) => {
                            let mut terms =
                                sum.iter().map(|term| term.get_data()).collect::<Vec<_>>();
                            terms.sort_unstable();
                            terms.dedup();
                            bases.push(terms);
                        }
                        _ => blocked = true,
                    }
                    // The sums inside an opaque power are not distributed.
                    return false;
                }
                true
            }
            _ => true,
        }
    });
    if blocked {
        return true;
    }
    if bases.is_empty() {
        return false;
    }
    let all_terms = bases.iter().flatten().copied().collect::<HashSet<&[u8]>>();
    expression.visitor(&mut |part| {
        if blocked {
            return false;
        }
        match part {
            AtomView::Fun(_) => false,
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                // Opaque powers were collected above; their sums stay whole.
                !(matches!(base, AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_))
                    && !i64::try_from(exponent).is_ok_and(|exponent| exponent > 0))
            }
            AtomView::Add(sum) => {
                let present = sum
                    .iter()
                    .map(|term| term.get_data())
                    .filter(|term| all_terms.contains(term))
                    .collect::<HashSet<_>>();
                if !present.is_empty()
                    && bases
                        .iter()
                        .any(|base| base.iter().all(|term| present.contains(term)))
                {
                    blocked = true;
                    return false;
                }
                true
            }
            _ => true,
        }
    });
    blocked
}

fn expand_pieces_unchecked(
    expression: AtomView<'_>,
    limit: usize,
    by_estimate: bool,
    sink: &mut dyn FnMut(Atom) -> Result<()>,
) -> Result<()> {
    const PARTS: usize = 16;
    if expansion_size(expression, by_estimate) < limit {
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
                    [term] => expand_pieces_unchecked(term, limit, by_estimate, sink),
                    _ => sink(Atom::add_many(chunk.iter().copied()).expand()),
                }
                .map(|()| chunk.clear())
            };
            for term in sum.iter() {
                let bytes = expansion_size(term, by_estimate);
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
                expand_pieces_unchecked(distributed.as_view(), limit, by_estimate, sink)?;
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

/// Pieces of [`map_term_local`]: `GAMMALOOP_UV_MAP_PIECE_BYTES` (default 64 MiB).
fn term_local_piece_bytes() -> usize {
    static BYTES: LazyLock<usize> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_UV_MAP_PIECE_BYTES")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(64 << 20)
    });
    *BYTES
}

/// A resource checkpoint of every budget active on this thread (the process RSS
/// guard of the beta driver), as an error naming `stage`.
pub(crate) fn guard_checkpoint(stage: &str) -> Result<()> {
    symbolica::utils::resources::checkpoint()
        .map_err(|error| color_eyre::eyre::eyre!("resource budget exhausted in {stage}: {error}"))
}

/// Apply a term-local map `f` to `expression` in pieces, with a resource
/// checkpoint before every piece.
///
/// `f` must commute with sums and products: `f(a + b) = f(a) + f(b)` and
/// `f(a * b) = f(a) * f(b)` after normalization. A replacement of function or
/// symbol patterns whose images do not depend on the surrounding term (momentum
/// routing, edge scaling, the edge-combination substitution) is such a map.
/// A sum is mapped in chunks of consecutive summands (a summand larger than a
/// piece is mapped recursively), a product factor by factor, anything else
/// whole; normalized sums and products are canonical, so the result is the
/// atom `f(expression)`. The whole-expression map holds the input, its copy
/// under construction and the normalization buffers at once (the r5 kernel
/// OOMs); here only one piece is in flight, and a process over its RSS guard
/// fails this diagram at the next piece instead of growing until the kernel
/// kills the process.
pub(crate) fn map_term_local(
    stage: &str,
    expression: AtomView<'_>,
    f: &dyn Fn(AtomView<'_>) -> Atom,
) -> Result<Atom> {
    map_term_local_with(stage, expression, term_local_piece_bytes(), f)
}

pub(crate) fn map_term_local_with(
    stage: &str,
    expression: AtomView<'_>,
    limit: usize,
    f: &dyn Fn(AtomView<'_>) -> Atom,
) -> Result<Atom> {
    guard_checkpoint(stage)?;
    if expression.get_byte_size() < limit {
        return Ok(f(expression));
    }
    match expression {
        AtomView::Add(sum) => {
            let mut accumulated = SumAccumulator::default();
            let mut chunk = Vec::new();
            let mut chunk_bytes = 0;
            let flush =
                |chunk: &mut Vec<AtomView<'_>>, accumulated: &mut SumAccumulator| -> Result<()> {
                    if chunk.is_empty() {
                        return Ok(());
                    }
                    guard_checkpoint(stage)?;
                    let piece = Atom::add_many(chunk.iter().copied());
                    accumulated.push(f(piece.as_view()));
                    chunk.clear();
                    Ok(())
                };
            for term in sum.iter() {
                let bytes = term.get_byte_size();
                if bytes >= limit {
                    flush(&mut chunk, &mut accumulated)?;
                    chunk_bytes = 0;
                    accumulated.push(map_term_local_with(stage, term, limit, f)?);
                    continue;
                }
                if chunk_bytes + bytes >= limit {
                    flush(&mut chunk, &mut accumulated)?;
                    chunk_bytes = 0;
                }
                chunk.push(term);
                chunk_bytes += bytes;
            }
            flush(&mut chunk, &mut accumulated)?;
            guard_checkpoint(stage)?;
            Ok(accumulated.finish())
        }
        AtomView::Mul(product) => {
            let mut factors = Vec::with_capacity(product.get_nargs());
            for factor in product.iter() {
                factors.push(map_term_local_with(stage, factor, limit, f)?);
            }
            guard_checkpoint(stage)?;
            Ok(Atom::mul_many(factors.iter()))
        }
        _ => Ok(f(expression)),
    }
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

/// Diagnostics: with `GAMMALOOP_ALGEBRA_DUMP=<directory>`, write the first
/// `GAMMALOOP_ALGEBRA_DUMP_COUNT` (default 40) inputs that are smaller than
/// `GAMMALOOP_ALGEBRA_DUMP_MAX_BYTES` (default 2 MB) to
/// `<directory>/<stage>-<n>.txt`, only of `GAMMALOOP_ALGEBRA_DUMP_STAGE` if set.
pub(crate) fn dump(stage: &str, expression: AtomView<'_>) {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static DIRECTORY: LazyLock<Option<String>> =
        LazyLock::new(|| std::env::var("GAMMALOOP_ALGEBRA_DUMP").ok());
    static COUNT: LazyLock<usize> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_ALGEBRA_DUMP_COUNT")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(40)
    });
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    static STAGE: LazyLock<Option<String>> =
        LazyLock::new(|| std::env::var("GAMMALOOP_ALGEBRA_DUMP_STAGE").ok());
    static MAX_BYTES: LazyLock<usize> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_ALGEBRA_DUMP_MAX_BYTES")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(2 << 20)
    });
    let Some(directory) = DIRECTORY.as_ref() else {
        return;
    };
    if STAGE
        .as_ref()
        .is_some_and(|only| !only.split(',').any(|entry| entry == stage))
    {
        return;
    }
    if expression.get_byte_size() > *MAX_BYTES {
        return;
    }
    let n = NEXT.fetch_add(1, Ordering::Relaxed);
    if n >= *COUNT {
        return;
    }
    let _ = std::fs::write(
        format!("{directory}/{stage}-{n}.txt"),
        expression.to_canonical_string(),
    );
}

/// `whole(expression)` for a map that acts summand by summand and, inside a
/// summand, only on its factors with explicit (or compact) indices: tensor
/// contraction and the dot-notation pass. A summand is a product `S * T` of
/// index-free factors `S` (numbers, symbols, `dot`s, denominators) and
/// indexed factors `T`; its image is `S * whole(T)`. The products `T` of a
/// five-loop integrand take a handful of shapes (the open external vectors
/// and metrics), so `whole` runs once per distinct `T`, cached in `images`,
/// instead of building a tensor network of the whole sum. Normalized sums and
/// products are canonical, so the rebuilt sum is the atom `whole` returns.
/// Returns `None` if an image introduces a dummy index (numbered over the
/// whole expression), so that the caller maps the whole expression.
pub(crate) fn map_indexed_factors<'a>(
    expression: AtomView<'a>,
    images: &mut std::collections::HashMap<Atom, Atom>,
    whole: &dyn Fn(AtomView<'_>) -> Result<Atom>,
) -> Result<Option<Atom>> {
    let mut matcher = SlotMatcher::default();
    let mut scalars = Vec::new();
    let mut indexed = Vec::new();
    let mut map_term = |term: AtomView<'a>| -> Result<Option<Atom>> {
        scalars.clear();
        indexed.clear();
        let factors = match term {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            factor => vec![factor],
        };
        for factor in factors {
            let interface = factor_interface(factor, &mut matcher);
            if interface.indices.is_empty() && !interface.opaque {
                scalars.push(factor);
            } else {
                indexed.push(factor);
            }
        }
        let key = match indexed[..] {
            [] => Atom::num(1),
            [single] => single.to_owned(),
            _ => Atom::mul_many(indexed.iter().copied()),
        };
        let image = match images.entry(key) {
            std::collections::hash_map::Entry::Occupied(entry) => entry.into_mut(),
            std::collections::hash_map::Entry::Vacant(entry) => {
                let image = whole(entry.key().as_view())?;
                if introduces_dummies(entry.key().as_view(), image.as_view()) {
                    return Ok(None);
                }
                entry.insert(image)
            }
        };
        if scalars.is_empty() {
            return Ok(Some(image.clone()));
        }
        Ok(Some(Atom::mul_many(
            scalars
                .iter()
                .copied()
                .chain(std::iter::once(image.as_view())),
        )))
    };
    match expression {
        AtomView::Add(sum) => {
            let mut terms = Vec::with_capacity(sum.get_nargs());
            for term in sum.iter() {
                let Some(mapped) = map_term(term)? else {
                    return Ok(None);
                };
                terms.push(mapped);
            }
            Ok(Some(Atom::add_many(terms.iter())))
        }
        term => map_term(term),
    }
}

/// Diagnostics: with `GAMMALOOP_ALGEBRA_SIZES=<bytes>`, report on stderr the
/// shape of every input of a stage of at least that many bytes
/// (`ALGEBRA_SIZE stage=... bytes=... kind=... args=...`).
pub(crate) fn report_size(stage: &str, expression: AtomView<'_>) {
    static THRESHOLD: LazyLock<Option<usize>> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_ALGEBRA_SIZES")
            .ok()
            .and_then(|value| value.parse().ok())
    });
    let Some(threshold) = *THRESHOLD else {
        return;
    };
    let bytes = expression.get_byte_size();
    if bytes < threshold {
        return;
    }
    let (kind, args) = match expression {
        AtomView::Add(sum) => ("add", sum.get_nargs()),
        AtomView::Mul(product) => ("mul", product.get_nargs()),
        AtomView::Pow(_) => ("pow", 2),
        AtomView::Fun(function) => ("fun", function.get_nargs()),
        _ => ("leaf", 0),
    };
    let largest = match expression {
        AtomView::Mul(product) => product.iter().map(|f| f.get_byte_size()).max().unwrap_or(0),
        AtomView::Add(sum) => sum.iter().map(|f| f.get_byte_size()).max().unwrap_or(0),
        _ => 0,
    };
    eprintln!(
        "ALGEBRA_SIZE stage={stage} bytes={bytes} kind={kind} args={args} largest_arg_bytes={largest}"
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::parse;

    #[test]
    fn expansion_estimate_counts_distributed_summands() {
        let expression = parse!("(a+b+c)*(d+e)*f + (g+h)^2 + k(x+y)");
        let (terms, bytes) = expansion_estimate(expression.as_view());
        assert_eq!(terms, 3.0 * 2.0 + 2.0 * 2.0 + 1.0);
        assert!(bytes > 0.0);
        let unexpanded = parse!("(a+b)*(c+d)");
        assert_eq!(
            expansion_size(unexpanded.as_view(), false),
            unexpanded.as_view().get_byte_size()
        );
        assert!(expansion_size(unexpanded.as_view(), true) >= unexpanded.as_view().get_byte_size());
    }

    #[test]
    fn term_local_map_in_pieces_is_the_whole_map() {
        use symbolica::id::Replacement;
        let expression = parse!(
            "(f(1,x)+2*f(2,y)+3)*(f(1,z)+g(4)+f(1,x)*f(2,x))*h(f(1,x)) + 5*f(1,y)^2 + f(3,x) + (f(1,x)+f(2,x))^3 + f(1,x)*f(1,y)"
        );
        let rules = vec![Replacement::new(
            parse!("f(1,x___)").to_pattern(),
            parse!("a*f(1,x___)+f(9,x___)").to_pattern(),
        )];
        let whole = expression.replace_multiple(&rules);
        for limit in [1, 8, 32, 64, 128, 1 << 20] {
            let pieces = map_term_local_with("test", expression.as_view(), limit, &|piece| {
                piece.replace_multiple(&rules)
            })
            .unwrap();
            assert_eq!(pieces, whole, "limit {limit}");
        }
    }

    #[test]
    fn pieces_by_estimate_add_up_to_the_whole_expansion() {
        let expression =
            parse!("(a+b+c+d)*(x1+x2+x3+x4+x5)*(y1+y2)^2*z + (p+q)*(r+s+t) + 3*u^2*(v+w)^(-1)");
        let whole = expression.expand();
        for limit in [1, 64, 256, 1024, 1 << 20] {
            for by_estimate in [false, true] {
                let mut pieces = Vec::new();
                expand_in_pieces_by(expression.as_view(), limit, by_estimate, &mut |piece| {
                    pieces.push(piece);
                    Ok(())
                })
                .unwrap();
                assert_eq!(
                    Atom::add_many(pieces.iter()),
                    whole,
                    "limit {limit} by_estimate {by_estimate}"
                );
                if by_estimate && limit <= 256 {
                    assert!(pieces.len() > 1, "limit {limit}: {} piece(s)", pieces.len());
                }
            }
        }
    }
}
