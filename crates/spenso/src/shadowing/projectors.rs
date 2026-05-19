use std::sync::LazyLock;

use crate::network::tags::SPENSO_TAG;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    symbol,
    utils::Settable,
};

use super::{
    IntoAtom, TensorCollectExt,
    trace::{chain_like_with_factors, trace, trace_factor_views},
};

/// Normalized symmetrizer over an ordered list of chain/trace factors.
///
/// This is an inert projector head. Its arguments are canonicalized by
/// Symbolica's `Symmetric` attribute, but it is only expanded by explicit
/// calls to [`ProjectorExpander::expand_projectors`].
pub static SYM: LazyLock<Symbol> = LazyLock::new(|| symbol!("spenso::sym"; Symmetric));

/// Normalized antisymmetrizer over an ordered list of chain/trace factors.
///
/// This is an inert projector head. Its arguments are canonicalized by
/// Symbolica's `Antisymmetric` attribute, but it is only expanded by explicit
/// calls to [`ProjectorExpander::expand_projectors`].
pub static ANTISYM: LazyLock<Symbol> = LazyLock::new(|| symbol!("spenso::antisym"; Antisymmetric));

/// Normalized cyclic symmetrizer over an ordered list of chain/trace factors.
///
/// This is an inert projector head. Its arguments are canonicalized by
/// Symbolica's `Cyclesymmetric` attribute, but it is only expanded by explicit
/// calls to [`ProjectorExpander::expand_projectors`]. In traces it is also the
/// compact representation of cyclic trace equivalence:
/// `trace(rep, cyclic(a,b,c))`.
pub static CYCLIC: LazyLock<Symbol> =
    LazyLock::new(|| symbol!("spenso::cyclic"; Cyclesymmetric; norm = normalize_cyclic_projector));

/// Builds an inert normalized symmetrizer over chain/trace factors.
pub fn sym<F>(factors: impl IntoIterator<Item = F>) -> Atom
where
    F: IntoAtom,
{
    projector(*SYM, factors)
}

/// Builds an inert normalized antisymmetrizer over chain/trace factors.
pub fn antisym<F>(factors: impl IntoIterator<Item = F>) -> Atom
where
    F: IntoAtom,
{
    projector(*ANTISYM, factors)
}

/// Builds an inert normalized cyclic symmetrizer over chain/trace factors.
pub fn cyclic<F>(factors: impl IntoIterator<Item = F>) -> Atom
where
    F: IntoAtom,
{
    projector(*CYCLIC, factors)
}

fn normalize_cyclic_projector(view: AtomView<'_>, out: &mut Settable<Atom>) {
    let AtomView::Fun(fun) = view else {
        return;
    };
    let args = fun.iter().collect::<Vec<_>>();
    if let [AtomView::Fun(inner)] = args.as_slice()
        && inner.get_symbol() == *SYM
    {
        out.set_from_view(&args[0]);
    }
}

fn projector<F>(symbol: Symbol, factors: impl IntoIterator<Item = F>) -> Atom
where
    F: IntoAtom,
{
    factors
        .into_iter()
        .fold(FunctionBuilder::new(symbol), |builder, factor| {
            builder.add_arg(factor.into_atom())
        })
        .finish()
}

pub trait ProjectorExpander {
    fn expand_projectors(&self) -> Atom;
}

impl ProjectorExpander for Atom {
    fn expand_projectors(&self) -> Atom {
        expand_projectors_impl(self.as_atom_view())
    }
}

impl ProjectorExpander for AtomView<'_> {
    fn expand_projectors(&self) -> Atom {
        expand_projectors_impl(*self)
    }
}

fn expand_projectors_impl(expression: AtomView) -> Atom {
    let mut current = expression.to_owned();
    loop {
        let next = current
            .replace_map(|arg, _context, out| {
                if let Some(expanded) = expand_chain_like_projector(arg) {
                    **out = expanded;
                }
            })
            .collect_tensors();

        if next == current {
            return next;
        }
        current = next;
    }
}

fn expand_chain_like_projector(arg: AtomView) -> Option<Atom> {
    let AtomView::Fun(f) = arg else {
        return None;
    };
    let offset = if f.get_symbol() == SPENSO_TAG.chain {
        2
    } else if f.get_symbol() == SPENSO_TAG.trace {
        1
    } else {
        return None;
    };
    if f.get_nargs() < offset {
        return None;
    }

    let args = f.iter().collect::<Vec<_>>();
    let (prefix, raw_factors) = args.split_at(offset);
    let factors = if f.get_symbol() == SPENSO_TAG.trace {
        trace_factor_views(raw_factors)
    } else {
        raw_factors.to_vec()
    };
    for (position, factor) in factors.iter().enumerate() {
        let Some((prefactor, projector_symbol, projector_factors)) = projector_factor(*factor)
        else {
            continue;
        };

        let mut sum = Atom::Zero;
        let is_trace = f.get_symbol() == SPENSO_TAG.trace;
        for (coefficient, expanded_factors) in
            projector_expansion_terms(&projector_factors, projector_symbol)?
        {
            let new_factors = factors[..position]
                .iter()
                .map(|factor| factor.to_owned())
                .chain(expanded_factors)
                .chain(
                    factors[position + 1..]
                        .iter()
                        .map(|factor| factor.to_owned()),
                );
            let rebuilt = if is_trace {
                trace_with_cyclic_factors(prefix, new_factors)
            } else {
                let prefix = prefix.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
                chain_like_with_factors(f.get_symbol(), &prefix, new_factors)
            };
            sum += prefactor.clone() * coefficient * rebuilt;
        }
        return Some(sum);
    }

    None
}

fn projector_factor(factor: AtomView) -> Option<(Atom, Symbol, Vec<Atom>)> {
    if let Some((projector_symbol, factors)) = projector_parts(factor) {
        return Some((Atom::num(1), projector_symbol, factors));
    }

    let AtomView::Mul(product) = factor else {
        return None;
    };

    let mut prefactor = Atom::num(1);
    let mut projector = None;
    for factor in product.iter() {
        if let Some(parts) = projector_parts(factor) {
            if projector.is_some() {
                return None;
            }
            projector = Some(parts);
        } else if matches!(factor, AtomView::Num(_)) {
            prefactor *= factor.to_owned();
        } else {
            return None;
        }
    }

    let (projector_symbol, factors) = projector?;
    Some((prefactor, projector_symbol, factors))
}

fn projector_parts(projector: AtomView) -> Option<(Symbol, Vec<Atom>)> {
    let AtomView::Fun(f) = projector else {
        return None;
    };
    if f.get_symbol() != *SYM && f.get_symbol() != *ANTISYM && f.get_symbol() != *CYCLIC {
        return None;
    }

    Some((f.get_symbol(), f.iter().map(|arg| arg.to_owned()).collect()))
}

fn projector_expansion_terms(
    factors: &[Atom],
    projector_symbol: Symbol,
) -> Option<Vec<(Atom, Vec<Atom>)>> {
    if projector_symbol == *CYCLIC {
        return cyclic_projector_terms(factors);
    }

    let denominator = Atom::num(factorial_i64(factors.len())?);
    let mut permutation = (0..factors.len()).collect::<Vec<_>>();
    let mut terms = Vec::new();
    collect_projector_terms(
        factors,
        projector_symbol == *ANTISYM,
        0,
        1,
        &mut permutation,
        &denominator,
        &mut terms,
    );
    Some(terms)
}

fn cyclic_projector_terms(factors: &[Atom]) -> Option<Vec<(Atom, Vec<Atom>)>> {
    if factors.is_empty() {
        return Some(vec![(Atom::num(1), Vec::new())]);
    }

    let denominator = Atom::num(i64::try_from(factors.len()).ok()?);
    Some(
        (0..factors.len())
            .map(|shift| {
                let rotated = factors
                    .iter()
                    .cycle()
                    .skip(shift)
                    .take(factors.len())
                    .cloned()
                    .collect();
                (Atom::num(1) / denominator.clone(), rotated)
            })
            .collect(),
    )
}

fn collect_projector_terms(
    factors: &[Atom],
    antisymmetric: bool,
    position: usize,
    sign: i64,
    permutation: &mut [usize],
    denominator: &Atom,
    terms: &mut Vec<(Atom, Vec<Atom>)>,
) {
    if position == permutation.len() {
        let coefficient = if antisymmetric {
            Atom::num(sign) / denominator.clone()
        } else {
            Atom::num(1) / denominator.clone()
        };
        let permuted = permutation
            .iter()
            .map(|index| factors[*index].clone())
            .collect();
        terms.push((coefficient, permuted));
        return;
    }

    for i in position..permutation.len() {
        permutation.swap(position, i);
        let next_sign = if i == position { sign } else { -sign };
        collect_projector_terms(
            factors,
            antisymmetric,
            position + 1,
            next_sign,
            permutation,
            denominator,
            terms,
        );
        permutation.swap(position, i);
    }
}

fn factorial_i64(n: usize) -> Option<i64> {
    (1..=n).try_fold(1_i64, |acc, value| {
        acc.checked_mul(i64::try_from(value).ok()?)
    })
}

fn trace_with_cyclic_factors(
    prefix: &[AtomView<'_>],
    factors: impl IntoIterator<Item = Atom>,
) -> Atom {
    debug_assert_eq!(prefix.len(), 1);
    trace(prefix[0], factors)
}
