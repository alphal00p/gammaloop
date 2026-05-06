use std::sync::LazyLock;

use crate::{
    network::tags::SPENSO_TAG,
    structure::{
        representation::{RepName, Representation},
        slot::{AbsInd, IsAbstractSlot, ParseableAind, Slot},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    symbol,
};

/// Converts common symbolic Spenso values into owned Symbolica atoms.
///
/// This is used by the variadic chain and trace macros so callers can pass
/// representation slots, stripped representations, symbols, atoms, or atom
/// views without spelling out the conversion at every call site.
pub trait IntoAtom {
    fn into_atom(self) -> Atom;
}

impl IntoAtom for Atom {
    fn into_atom(self) -> Atom {
        self
    }
}

impl IntoAtom for &Atom {
    fn into_atom(self) -> Atom {
        self.clone()
    }
}

impl IntoAtom for AtomView<'_> {
    fn into_atom(self) -> Atom {
        self.to_owned()
    }
}

impl IntoAtom for AtomOrView<'_> {
    fn into_atom(self) -> Atom {
        self.into_owned()
    }
}

impl IntoAtom for Symbol {
    fn into_atom(self) -> Atom {
        Atom::var(self)
    }
}

impl<R, A> IntoAtom for Slot<R, A>
where
    R: RepName,
    A: AbsInd + ParseableAind,
{
    fn into_atom(self) -> Atom {
        self.to_atom()
    }
}

impl<R> IntoAtom for Representation<R>
where
    R: RepName,
{
    fn into_atom(self) -> Atom {
        self.to_symbolic(std::iter::empty::<Atom>())
    }
}

impl<R> IntoAtom for &Representation<R>
where
    R: RepName,
{
    fn into_atom(self) -> Atom {
        self.to_symbolic(std::iter::empty::<Atom>())
    }
}

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

/// Expands `sym(...)` and `antisym(...)` factors inside `chain(...)` and
/// `trace(...)`.
///
/// The projectors are interpreted as normalized sums over permutations of the
/// factor list. They are not expanded when they appear outside a chain-like
/// container, because the ordering only has a noncommutative meaning in that
/// context.
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
            .expand();

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

    let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
    let (prefix, factors) = args.split_at(offset);
    for (position, factor) in factors.iter().enumerate() {
        let Some((prefactor, projector_symbol, projector_factors)) =
            projector_factor(factor.as_view())
        else {
            continue;
        };

        let mut sum = Atom::Zero;
        for (coefficient, expanded_factors) in
            projector_expansion_terms(&projector_factors, projector_symbol == *ANTISYM)?
        {
            let new_factors = factors[..position]
                .iter()
                .cloned()
                .chain(expanded_factors)
                .chain(factors[position + 1..].iter().cloned());
            let rebuilt = chain_like_with_factors(f.get_symbol(), prefix, new_factors);
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
    if f.get_symbol() != *SYM && f.get_symbol() != *ANTISYM {
        return None;
    }

    Some((f.get_symbol(), f.iter().map(|arg| arg.to_owned()).collect()))
}

fn projector_expansion_terms(
    factors: &[Atom],
    antisymmetric: bool,
) -> Option<Vec<(Atom, Vec<Atom>)>> {
    let denominator = Atom::num(factorial_i64(factors.len())?);
    let mut permutation = (0..factors.len()).collect::<Vec<_>>();
    let mut terms = Vec::new();
    collect_projector_terms(
        factors,
        antisymmetric,
        0,
        1,
        &mut permutation,
        &denominator,
        &mut terms,
    );
    Some(terms)
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

fn chain_like_with_factors(
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

/// Creates a Symbolica symbol from an identifier.
///
/// This is a small shorthand for `symbolica::symbol!(stringify!(name))`.
/// It is useful in expression-building code where many symbolic abstract
/// indices are needed.
///
/// # Examples
///
/// ```ignore
/// use spenso::s;
///
/// let mu = s!(mu);
/// let dim = s!(D);
/// ```
#[macro_export]
macro_rules! s {
    ($name:ident) => {
        symbolica::symbol!(stringify!($name))
    };
}

/// Creates an abstract-index slot from a representation and an index.
///
/// The slot index is fixed to Spenso's symbolic `AbstractIndex`, which avoids
/// type-inference ambiguity when passing `Symbol`s. The second argument can be
/// an identifier, expanded through [`s!`], an integer index, or an arbitrary
/// expression convertible into `AbstractIndex`.
///
/// # Examples
///
/// ```ignore
/// use spenso::{s, slot};
/// use spenso::structure::representation::{Minkowski, RepName};
///
/// let mink4 = Minkowski {}.new_rep(4);
/// let mu = slot!(mink4, mu);
/// let nu = slot!(mink4, s!(nu));
/// let one = slot!(mink4, 1);
/// ```
#[macro_export]
macro_rules! slot {
    ($rep:expr, $index:ident) => {
        ($rep).slot::<$crate::structure::abstract_index::AbstractIndex, _>($crate::s!($index))
    };
    ($rep:expr, $index:expr) => {
        ($rep).slot::<$crate::structure::abstract_index::AbstractIndex, _>($index)
    };
}

/// Builds an inert normalized symmetrizer over chain/trace factors.
///
/// The resulting `sym(...)` is canonicalized through Symbolica's symmetric
/// function attribute. It expands only when explicitly passed through
/// [`ProjectorExpander::expand_projectors`].
///
/// # Examples
///
/// ```ignore
/// use spenso::{chain, sym};
///
/// let projected = chain!(i, j, sym!(a, b, c));
/// let from_vec = sym!(; vec![a, b, c]);
/// ```
#[macro_export]
macro_rules! sym {
    (; $factors:expr $(,)?) => {
        $crate::symbolica_atom::sym(($factors).into_iter())
    };
    ($($factor:expr),* $(,)?) => {
        $crate::symbolica_atom::sym(vec![$($crate::symbolica_atom::IntoAtom::into_atom($factor)),*])
    };
}

/// Builds an inert normalized antisymmetrizer over chain/trace factors.
///
/// The resulting `antisym(...)` is canonicalized through Symbolica's
/// antisymmetric function attribute. It expands only when explicitly passed
/// through [`ProjectorExpander::expand_projectors`].
///
/// # Examples
///
/// ```ignore
/// use spenso::{antisym, chain};
///
/// let commutator_part = chain!(i, j, antisym!(a, b));
/// let from_vec = antisym!(; vec![a, b, c]);
/// ```
#[macro_export]
macro_rules! antisym {
    (; $factors:expr $(,)?) => {
        $crate::symbolica_atom::antisym(($factors).into_iter())
    };
    ($($factor:expr),* $(,)?) => {
        $crate::symbolica_atom::antisym(vec![$($crate::symbolica_atom::IntoAtom::into_atom($factor)),*])
    };
}

/// Builds a symbolic open chain with a start slot, an end slot, and any number
/// of factors.
///
/// Arguments are converted through [`IntoAtom`], so callers can pass Spenso
/// slots, atoms, atom views, or symbols. This macro is only the variadic surface
/// for [`SPENSO_TAG.chain`](crate::network::tags::SPENSO_TAG); it does not
/// simplify or normalize the expression.
///
/// # Examples
///
/// ```ignore
/// use spenso::{chain, slot};
///
/// let expr = chain!(
///     slot!(bis4, a),
///     slot!(bis4, b),
///     gamma_mu,
///     gamma_nu,
/// );
///
/// let factors = vec![gamma_mu, gamma_nu];
/// let expr = chain!(slot!(bis4, a), slot!(bis4, b); factors);
/// ```
#[macro_export]
macro_rules! chain {
    ($start:expr, $end:expr; $factors:expr $(,)?) => {
        $crate::network::tags::SPENSO_TAG.chain(
            $crate::symbolica_atom::IntoAtom::into_atom($start),
            $crate::symbolica_atom::IntoAtom::into_atom($end),
            ($factors)
                .into_iter()
                .map($crate::symbolica_atom::IntoAtom::into_atom),
        )
    };
    ($start:expr, $end:expr $(, $factor:expr)+; $factors:expr $(,)?) => {
        $crate::network::tags::SPENSO_TAG.chain(
            $crate::symbolica_atom::IntoAtom::into_atom($start),
            $crate::symbolica_atom::IntoAtom::into_atom($end),
            vec![$($crate::symbolica_atom::IntoAtom::into_atom($factor)),*]
                .into_iter()
                .chain(($factors).into_iter().map($crate::symbolica_atom::IntoAtom::into_atom)),
        )
    };
    ($start:expr, $end:expr $(, $factor:expr)* $(,)?) => {
        $crate::network::tags::SPENSO_TAG.chain(
            $crate::symbolica_atom::IntoAtom::into_atom($start),
            $crate::symbolica_atom::IntoAtom::into_atom($end),
            vec![$($crate::symbolica_atom::IntoAtom::into_atom($factor)),*],
        )
    };
}

/// Builds a symbolic trace with a representation and any number of factors.
///
/// The representation and factors are converted through [`IntoAtom`]. This is
/// the variadic surface for [`SPENSO_TAG.trace`](crate::network::tags::SPENSO_TAG)
/// and intentionally leaves any algebraic simplification to the caller.
///
/// # Examples
///
/// ```ignore
/// use spenso::trace;
///
/// let expr = trace!(
///     &cof_nc,
///     color_t_a,
///     color_t_b,
/// );
///
/// let factors = vec![color_t_a, color_t_b];
/// let expr = trace!(&cof_nc; factors);
/// ```
#[macro_export]
macro_rules! trace {
    ($rep:expr; $factors:expr $(,)?) => {
        $crate::network::tags::SPENSO_TAG.trace(
            $crate::symbolica_atom::IntoAtom::into_atom($rep),
            ($factors)
                .into_iter()
                .map($crate::symbolica_atom::IntoAtom::into_atom),
        )
    };
    ($rep:expr $(, $factor:expr)+; $factors:expr $(,)?) => {
        $crate::network::tags::SPENSO_TAG.trace(
            $crate::symbolica_atom::IntoAtom::into_atom($rep),
            vec![$($crate::symbolica_atom::IntoAtom::into_atom($factor)),*]
                .into_iter()
                .chain(($factors).into_iter().map($crate::symbolica_atom::IntoAtom::into_atom)),
        )
    };
    ($rep:expr $(, $factor:expr)* $(,)?) => {
        $crate::network::tags::SPENSO_TAG.trace(
            $crate::symbolica_atom::IntoAtom::into_atom($rep),
            vec![$($crate::symbolica_atom::IntoAtom::into_atom($factor)),*],
        )
    };
}

#[cfg(test)]
mod tests {
    use crate::symbolica_atom::ProjectorExpander;
    use symbolica::{atom::Atom, symbol};

    #[test]
    fn chain_macro_accepts_iterable_factors() {
        let start = Atom::var(symbol!("start"));
        let end = Atom::var(symbol!("end"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        assert_eq!(
            crate::chain!(start.clone(), end.clone(); vec![first.clone(), second.clone()]),
            crate::chain!(start, end, first, second)
        );
    }

    #[test]
    fn trace_macro_accepts_iterable_factors() {
        let rep = Atom::var(symbol!("rep"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        assert_eq!(
            crate::trace!(rep.clone(); vec![first.clone(), second.clone()]),
            crate::trace!(rep, first, second)
        );
    }

    #[test]
    fn sym_macro_canonicalizes_arguments() {
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        assert_eq!(
            crate::sym!(second.clone(), first.clone()),
            crate::sym!(first, second)
        );
    }

    #[test]
    fn expand_sym_in_chain_preserves_endpoints() {
        let start = Atom::var(symbol!("start"));
        let end = Atom::var(symbol!("end"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        let expanded = crate::chain!(
            start.clone(),
            end.clone(),
            crate::sym!(first.clone(), second.clone())
        )
        .expand_projectors();
        let expected = Atom::num(1) / Atom::num(2)
            * crate::chain!(start.clone(), end.clone(), first.clone(), second.clone())
            + Atom::num(1) / Atom::num(2) * crate::chain!(start, end, second, first);

        assert_eq!(expanded, expected);
    }

    #[test]
    fn expand_antisym_in_trace_preserves_sign() {
        let rep = Atom::var(symbol!("rep"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        let expanded = crate::trace!(rep.clone(), crate::antisym!(first.clone(), second.clone()))
            .expand_projectors();
        let expected = Atom::num(1) / Atom::num(2)
            * crate::trace!(rep.clone(), first.clone(), second.clone())
            - Atom::num(1) / Atom::num(2) * crate::trace!(rep, second, first);

        assert_eq!(expanded, expected);
    }

    #[test]
    fn expand_antisym_in_chain_handles_canonicalization_sign() {
        let start = Atom::var(symbol!("start"));
        let end = Atom::var(symbol!("end"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        let expanded = crate::chain!(
            start.clone(),
            end.clone(),
            crate::antisym!(second.clone(), first.clone())
        )
        .expand_projectors();
        let expected = Atom::num(-1) / Atom::num(2)
            * crate::chain!(start.clone(), end.clone(), first.clone(), second.clone())
            + Atom::num(1) / Atom::num(2) * crate::chain!(start, end, second, first);

        assert_eq!(expanded, expected);
    }
}
