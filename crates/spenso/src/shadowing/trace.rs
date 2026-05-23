use crate::network::tags::SPENSO_TAG;

use super::{
    IntoAtom,
    projectors::{CYCLIC, cyclic, sym},
};
use symbolica::atom::{Atom, AtomView, FunctionBuilder, Symbol, representation::FunView};

pub fn trace<R, F>(rep: R, factors: impl IntoIterator<Item = F>) -> Atom
where
    R: IntoAtom,
    F: IntoAtom,
{
    let rep = rep.into_atom();
    let factors = factors
        .into_iter()
        .map(IntoAtom::into_atom)
        .collect::<Vec<_>>();

    if factors.is_empty() {
        return chain_like_with_factors(SPENSO_TAG.trace, &[rep], std::iter::empty());
    }

    let trace_factor = if let [factor] = factors.as_slice() {
        if is_cyclic_projector(factor.as_view()) {
            factor.clone()
        } else {
            cyclic(factors)
        }
    } else {
        cyclic(factors)
    };

    chain_like_with_factors(SPENSO_TAG.trace, &[rep], std::iter::once(trace_factor))
}

/// Builds a trace whose payload is a single symmetric projector.
///
/// This is the compact form for fully symmetric trace invariants such as color
/// `d` tensors: `trace(rep, sym(factors...))`. Ordinary ordered closed traces
/// should use [`trace`], which wraps factors in `cyclic(...)` instead.
pub fn trace_sym<R, F>(rep: R, factors: impl IntoIterator<Item = F>) -> Atom
where
    R: IntoAtom,
    F: IntoAtom,
{
    let rep = rep.into_atom();
    let projector = sym(factors);
    chain_like_with_factors(SPENSO_TAG.trace, &[rep], std::iter::once(projector))
}

/// Return the representation and actual factor sequence of a trace.
///
/// Canonical traces store their factors under a single `cyclic(...)` argument;
/// this helper unwraps that representation for algorithms that need the
/// ordered factor list. The non-canonical raw shape `trace(rep, factors...)` is
/// accepted on input so readers can normalize externally constructed atoms, but
/// all builders emit the cyclic form.
pub fn trace_parts<'a>(fun: FunView<'a>) -> Option<(AtomView<'a>, Vec<AtomView<'a>>)> {
    if fun.get_symbol() != SPENSO_TAG.trace {
        return None;
    }

    let args = fun.iter().collect::<Vec<_>>();
    let [rep, factors @ ..] = args.as_slice() else {
        return None;
    };

    Some((*rep, trace_factor_views(factors)))
}

/// Flatten the canonical cyclic trace factor wrapper.
pub fn trace_factor_views<'a>(factors: &[AtomView<'a>]) -> Vec<AtomView<'a>> {
    if let [factor] = factors
        && let AtomView::Fun(fun) = *factor
        && fun.get_symbol() == *CYCLIC
    {
        return fun.iter().collect();
    }

    factors.to_vec()
}

fn is_cyclic_projector(factor: AtomView<'_>) -> bool {
    matches!(factor, AtomView::Fun(fun) if fun.get_symbol() == *CYCLIC)
}

pub(crate) fn chain_like_with_factors(
    symbol: Symbol,
    prefix: &[Atom],
    factors: impl IntoIterator<Item = Atom>,
) -> Atom {
    prefix
        .iter()
        .cloned()
        .chain(factors)
        .fold(FunctionBuilder::new(symbol), |builder, arg| {
            builder.add_arg(arg)
        })
        .finish()
}
