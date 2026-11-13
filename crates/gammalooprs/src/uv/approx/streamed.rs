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
