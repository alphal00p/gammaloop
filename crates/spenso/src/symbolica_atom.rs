use std::{collections::BTreeSet, sync::LazyLock};

use crate::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        abstract_index::AIND_SYMBOLS,
        representation::{RepName, Representation},
        slot::{AbsInd, IsAbstractSlot, ParseableAind, Slot},
    },
};
use symbolica::{
    atom::{
        Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol, representation::FunView,
    },
    symbol,
    utils::Settable,
};

/// Converts common symbolic Spenso values into owned Symbolica atoms.
///
/// This is used by the variadic chain and trace macros so callers can pass
/// representation slots, stripped representations, symbols, atoms, or atom
/// views without spelling out the conversion at every call site.
pub trait IntoAtom {
    fn into_atom(self) -> Atom;
}

macro_rules! impl_into_atom_via_atom_or_view {
    (<$lt:lifetime> $($ty:ty),+ $(,)?) => {
        $(
            impl<$lt> IntoAtom for $ty {
                fn into_atom(self) -> Atom {
                    let atom: AtomOrView<$lt> = self.into();
                    atom.into_owned()
                }
            }
        )+
    };
    ($($ty:ty),+ $(,)?) => {
        $(
            impl IntoAtom for $ty {
                fn into_atom(self) -> Atom {
                    let atom: AtomOrView<'_> = self.into();
                    atom.into_owned()
                }
            }
        )+
    };
}

impl_into_atom_via_atom_or_view!(Atom, Symbol);
impl_into_atom_via_atom_or_view!(<'a> &'a Atom, AtomView<'a>, AtomOrView<'a>);

macro_rules! impl_into_atom_integer {
    ($($ty:ty),+ $(,)?) => {
        $(
            impl IntoAtom for $ty {
                fn into_atom(self) -> Atom {
                    Atom::num(self)
                }
            }
        )+
    };
}

macro_rules! impl_into_atom_integer_via {
    ($target:ty; $($ty:ty),+ $(,)?) => {
        $(
            impl IntoAtom for $ty {
                fn into_atom(self) -> Atom {
                    Atom::num(<$target>::from(self))
                }
            }
        )+
    };
}

impl_into_atom_integer!(i32, i64, isize, u32, u64, usize);
impl_into_atom_integer_via!(i64; i8, i16);
impl_into_atom_integer_via!(u64; u8, u16);

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

/// Normalized cyclic symmetrizer over an ordered list of chain/trace factors.
///
/// This is an inert projector head. Its arguments are canonicalized by
/// Symbolica's `Cyclesymmetric` attribute, but it is only expanded by explicit
/// calls to [`ProjectorExpander::expand_projectors`]. In traces it is also the
/// compact representation of cyclic trace equivalence:
/// `trace(rep, cyclic(a,b,c))`.
pub static CYCLIC: LazyLock<Symbol> =
    LazyLock::new(|| symbol!("spenso::cyclic"; Cyclesymmetric; norm = normalize_cyclic_projector));

/// Temporary head used to collect tensor factors without expanding products.
///
/// The collection helpers wrap maximal matching tensorial subexpressions as
/// `spenso::collect(...)`, call Symbolica's `collect_symbol` on this inert
/// head, then unwrap it again. This gives simplifiers a way to combine terms
/// with common tensor factors without distributing scalar sums through the
/// whole expression.
pub static COLLECT: LazyLock<Symbol> = LazyLock::new(|| symbol!("spenso::collect"));

/// Selects which tensorial subexpressions are protected during collection.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TensorCollectFilter {
    /// Collect every syntactically tensorial subexpression.
    Tensors,
    /// Collect tensorial subexpressions that contain this representation.
    Rep(crate::structure::representation::LibraryRep),
    /// Collect only metric tensors.
    Metrics,
}

impl TensorCollectFilter {
    fn collect(self, expression: AtomView<'_>) -> Atom {
        let wrapped = expression.replace_map(|arg, _context, out| {
            if self.matches(arg) {
                **out = self.wrap(arg);
            }
        });

        if wrapped == expression.to_owned() {
            return wrapped;
        }

        let collected = wrapped.collect_symbol::<i16>(*COLLECT, None, None);
        self.unwrap(collected.as_view())
    }

    fn matches(self, arg: AtomView<'_>) -> bool {
        let AtomView::Fun(fun) = arg else {
            return false;
        };

        match self {
            Self::Metrics => fun.get_symbol() == ETS.metric,
            Self::Tensors => {
                !Self::is_non_tensorial_wrapper(fun.get_symbol()) && Self::contains_rep(fun)
            }
            Self::Rep(rep) => {
                !Self::is_non_tensorial_wrapper(fun.get_symbol())
                    && fun.iter().any(|arg| Self::arg_contains_rep(arg, rep))
            }
        }
    }

    fn wrap(self, arg: AtomView<'_>) -> Atom {
        FunctionBuilder::new(*COLLECT).add_arg(arg).finish()
    }

    fn unwrap(self, expression: AtomView<'_>) -> Atom {
        expression.replace_map(|arg, _context, out| {
            if let Some(inner) = Self::collect_inner(arg) {
                **out = inner.to_owned();
            }
        })
    }

    fn collect_inner(arg: AtomView<'_>) -> Option<AtomView<'_>> {
        let AtomView::Fun(fun) = arg else {
            return None;
        };
        if fun.get_symbol() != *COLLECT || fun.get_nargs() != 1 {
            return None;
        }
        fun.iter().next()
    }

    fn contains_rep(fun: FunView<'_>) -> bool {
        fun.iter().any(Self::arg_is_or_contains_rep)
    }

    fn arg_is_or_contains_rep(arg: AtomView<'_>) -> bool {
        if Self::is_rep(arg) {
            return true;
        }

        match arg {
            AtomView::Fun(fun) if !Self::is_non_tensorial_wrapper(fun.get_symbol()) => {
                fun.iter().any(Self::arg_is_or_contains_rep)
            }
            AtomView::Mul(mul) => mul.iter().any(Self::arg_is_or_contains_rep),
            AtomView::Add(add) => add.iter().any(Self::arg_is_or_contains_rep),
            AtomView::Pow(pow) => Self::arg_is_or_contains_rep(pow.get_base_exp().0),
            _ => false,
        }
    }

    fn arg_contains_rep(
        arg: AtomView<'_>,
        rep: crate::structure::representation::LibraryRep,
    ) -> bool {
        if Self::is_exact_rep(arg, rep) {
            return true;
        }

        match arg {
            AtomView::Fun(fun) if !Self::is_non_tensorial_wrapper(fun.get_symbol()) => {
                fun.iter().any(|arg| Self::arg_contains_rep(arg, rep))
            }
            AtomView::Mul(mul) => mul.iter().any(|arg| Self::arg_contains_rep(arg, rep)),
            AtomView::Add(add) => add.iter().any(|arg| Self::arg_contains_rep(arg, rep)),
            AtomView::Pow(pow) => Self::arg_contains_rep(pow.get_base_exp().0, rep),
            _ => false,
        }
    }

    fn is_exact_rep(arg: AtomView<'_>, rep: crate::structure::representation::LibraryRep) -> bool {
        Representation::<crate::structure::representation::LibraryRep>::try_from(arg)
            .is_ok_and(|found| found.rep.base() == rep.base())
    }

    fn is_rep(arg: AtomView<'_>) -> bool {
        matches!(arg, AtomView::Fun(fun) if fun.get_symbol().has_tag(&SPENSO_TAG.representation))
            || arg.has_attributes_of(SPENSO_TAG.rep_)
    }

    fn is_non_tensorial_wrapper(symbol: Symbol) -> bool {
        symbol == *COLLECT
            || symbol == AIND_SYMBOLS.aind
            || symbol == AIND_SYMBOLS.dind
            || symbol == AIND_SYMBOLS.uind
            || symbol == AIND_SYMBOLS.selfdualind
            || symbol.has_tag(&SPENSO_TAG.representation)
    }
}

#[derive(Default)]
struct ScalarSymbolCollector {
    symbols: BTreeSet<Symbol>,
}

impl ScalarSymbolCollector {
    fn collect(expression: AtomView<'_>) -> Atom {
        let mut collector = Self::default();
        collector.visit(expression);

        collector
            .symbols
            .into_iter()
            .fold(expression.to_owned(), |expr, symbol| {
                expr.collect_symbol::<i16>(symbol, None, None)
            })
    }

    fn visit(&mut self, expression: AtomView<'_>) {
        match expression {
            AtomView::Var(var) => self.push(var.get_symbol()),
            AtomView::Add(add) => add.iter().for_each(|arg| self.visit(arg)),
            AtomView::Mul(mul) => mul.iter().for_each(|arg| self.visit(arg)),
            AtomView::Pow(pow) => self.visit(pow.get_base_exp().0),
            AtomView::Fun(_) | AtomView::Num(_) => {}
        }
    }

    fn push(&mut self, symbol: Symbol) {
        if symbol.get_wildcard_level() == 0 {
            self.symbols.insert(symbol);
        }
    }
}

/// Collect tensor factors without full expression expansion.
pub trait TensorCollectExt {
    /// Collect common tensor leaves by temporarily wrapping them in `spenso::collect(...)`.
    fn collect_tensors(&self) -> Atom;

    /// Collect common tensor leaves that contain `rep` as one of their slot representations.
    fn collect_rep(&self, rep: crate::structure::representation::LibraryRep) -> Atom;

    /// Collect common metric tensors.
    fn collect_metrics(&self) -> Atom;
}

/// Collect visible scalar variables without distributing products.
pub trait ScalarCollectExt {
    /// Repeatedly call Symbolica's `collect_symbol` on scalar variable factors.
    fn collect_scalar_symbols(&self) -> Atom;
}

impl ScalarCollectExt for Atom {
    fn collect_scalar_symbols(&self) -> Atom {
        self.as_view().collect_scalar_symbols()
    }
}

impl ScalarCollectExt for AtomView<'_> {
    fn collect_scalar_symbols(&self) -> Atom {
        ScalarSymbolCollector::collect(*self)
    }
}

impl TensorCollectExt for Atom {
    fn collect_tensors(&self) -> Atom {
        self.as_view().collect_tensors()
    }

    fn collect_rep(&self, rep: crate::structure::representation::LibraryRep) -> Atom {
        self.as_view().collect_rep(rep)
    }

    fn collect_metrics(&self) -> Atom {
        self.as_view().collect_metrics()
    }
}

impl TensorCollectExt for AtomView<'_> {
    fn collect_tensors(&self) -> Atom {
        TensorCollectFilter::Tensors.collect(*self)
    }

    fn collect_rep(&self, rep: crate::structure::representation::LibraryRep) -> Atom {
        TensorCollectFilter::Rep(rep).collect(*self)
    }

    fn collect_metrics(&self) -> Atom {
        TensorCollectFilter::Metrics.collect(*self)
    }
}

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

/// Builds a trace in canonical projector form.
///
/// Empty traces are represented as `trace(rep)`. Non-empty traces are
/// represented as `trace(rep, cyclic(factors...))`. Passing a single existing
/// `cyclic(...)` factor is accepted and not wrapped again. The `cyclic`
/// constructor itself normalizes `cyclic(sym(...))` to `sym(...)`.
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

/// Expands `sym(...)`, `antisym(...)`, and `cyclic(...)` factors inside
/// `chain(...)` and `trace(...)`.
///
/// `sym` and `antisym` are interpreted as normalized sums over all
/// permutations of the factor list. `cyclic` is interpreted as a normalized sum
/// over cyclic rotations. Inside traces, expanded terms are emitted as
/// `trace(rep, cyclic(...))`, so trace cyclicity stays explicit in the
/// expression instead of being handled by a side-channel grouping pass.
/// Projectors are not expanded outside chain-like containers because the
/// ordering only has a noncommutative meaning in that context.
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

fn trace_with_cyclic_factors(
    prefix: &[AtomView<'_>],
    factors: impl IntoIterator<Item = Atom>,
) -> Atom {
    debug_assert_eq!(prefix.len(), 1);
    trace(prefix[0], factors)
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

#[doc(hidden)]
#[macro_export]
macro_rules! spenso_rep_atom {
    ($rep:expr; $dim:ident $(,)?) => {{
        let rep = $crate::structure::representation::RepName::new_rep(&$rep, $crate::s!($dim));
        rep.to_symbolic(std::iter::empty::<symbolica::atom::Atom>())
    }};
    ($rep:expr; $dim:expr $(,)?) => {{
        let rep = $crate::structure::representation::RepName::new_rep(&$rep, $dim);
        rep.to_symbolic(std::iter::empty::<symbolica::atom::Atom>())
    }};
    ($rep:expr; $dim:ident, $index:ident $(,)?) => {{
        let rep = $crate::structure::representation::RepName::new_rep(&$rep, $crate::s!($dim));
        rep.to_symbolic([$crate::symbolica_atom::IntoAtom::into_atom($crate::s!(
            $index
        ))])
    }};
    ($rep:expr; $dim:ident, $index:expr $(,)?) => {{
        let rep = $crate::structure::representation::RepName::new_rep(&$rep, $crate::s!($dim));
        rep.to_symbolic([$crate::symbolica_atom::IntoAtom::into_atom($index)])
    }};
    ($rep:expr; $dim:expr, $index:ident $(,)?) => {{
        let rep = $crate::structure::representation::RepName::new_rep(&$rep, $dim);
        rep.to_symbolic([$crate::symbolica_atom::IntoAtom::into_atom($crate::s!(
            $index
        ))])
    }};
    ($rep:expr; $dim:expr, $index:expr $(,)?) => {{
        let rep = $crate::structure::representation::RepName::new_rep(&$rep, $dim);
        rep.to_symbolic([$crate::symbolica_atom::IntoAtom::into_atom($index)])
    }};
}

/// Builds a stripped or indexed Minkowski representation atom.
#[macro_export]
macro_rules! mink {
    ($($args:tt)*) => {
        $crate::spenso_rep_atom!($crate::structure::representation::Minkowski {}; $($args)*)
    };
}

/// Builds a stripped or indexed Euclidean representation atom.
#[macro_export]
macro_rules! euc {
    ($($args:tt)*) => {
        $crate::spenso_rep_atom!($crate::structure::representation::Euclidean {}; $($args)*)
    };
}

/// Builds a stripped or indexed Lorentz representation atom.
#[macro_export]
macro_rules! lor {
    ($($args:tt)*) => {
        $crate::spenso_rep_atom!($crate::structure::representation::Lorentz {}; $($args)*)
    };
}

/// Builds a symbolic tensor function from a name and tensor arguments.
///
/// Identifier heads are created with [`tensor_symbol!`], so they keep the
/// caller's symbol namespace and carry the generic Spenso tensor tag. Pass an
/// explicit `Symbol` expression when the head was built elsewhere.
///
/// Arguments are converted through [`IntoAtom`], so callers can mix scalar
/// literal arguments, atoms, slots, and stripped representations. Passing a
/// representation, for example `tensor!(p, rep)`, emits compact Schoonschip
/// syntax; passing a slot, for example `tensor!(p, slot!(rep, i))`, emits an
/// explicitly indexed tensor.
///
/// # Examples
///
/// ```ignore
/// use spenso::{p, q, slot, tensor};
///
/// let compact = tensor!(k, 1, mink4);
/// let indexed = p!(1, slot!(mink4, mu)) * q!(2, slot!(mink4, mu));
/// ```
#[macro_export]
macro_rules! tensor {
    ($name:ident $(,)?) => {
        symbolica::atom::FunctionBuilder::new($crate::tensor_symbol!($name)).finish()
    };
    ($name:ident, $($arg:expr),+ $(,)?) => {{
        let mut tensor = symbolica::atom::FunctionBuilder::new($crate::tensor_symbol!($name));
        $(
            tensor = tensor.add_arg($crate::symbolica_atom::IntoAtom::into_atom($arg));
        )*
        tensor.finish()
    }};
    ($name:expr $(,)?) => {
        symbolica::atom::FunctionBuilder::new($name).finish()
    };
    ($name:expr, $($arg:expr),+ $(,)?) => {{
        let mut tensor = symbolica::atom::FunctionBuilder::new($name);
        $(
            tensor = tensor.add_arg($crate::symbolica_atom::IntoAtom::into_atom($arg));
        )*
        tensor.finish()
    }};
}

/// Builds a symbolic vector function from a name and tensor arguments.
///
/// The head is created with [`vector_symbol!`], so `vector!(p, ...)` is a
/// rank-one tensor in the caller's symbol namespace.
#[macro_export]
macro_rules! vector {
    ($name:ident $(,)?) => {
        symbolica::atom::FunctionBuilder::new($crate::vector_symbol!($name)).finish()
    };
    ($name:ident, $($arg:expr),+ $(,)?) => {{
        let mut tensor = symbolica::atom::FunctionBuilder::new($crate::vector_symbol!($name));
        $(
            tensor = tensor.add_arg($crate::symbolica_atom::IntoAtom::into_atom($arg));
        )*
        tensor.finish()
    }};
    ($name:expr $(,)?) => {
        symbolica::atom::FunctionBuilder::new($name).finish()
    };
    ($name:expr, $($arg:expr),+ $(,)?) => {{
        let mut tensor = symbolica::atom::FunctionBuilder::new($name);
        $(
            tensor = tensor.add_arg($crate::symbolica_atom::IntoAtom::into_atom($arg));
        )*
        tensor.finish()
    }};
}

/// Builds a symbolic vector `p(...)`.
///
/// This is a convenience wrapper around [`vector!`] with a rank-one tensor head.
#[macro_export]
macro_rules! p {
    () => {
        $crate::vector!(p)
    };
    ($($arg:expr),+ $(,)?) => {
        $crate::vector!(p, $($arg),+)
    };
}

/// Builds a symbolic vector `q(...)`.
///
/// This is a convenience wrapper around [`vector!`] with a rank-one tensor head.
#[macro_export]
macro_rules! q {
    () => {
        $crate::vector!(q)
    };
    ($($arg:expr),+ $(,)?) => {
        $crate::vector!(q, $($arg),+)
    };
}

/// Builds a metric tensor atom `g(a,b)`.
#[macro_export]
macro_rules! metric {
    ($a:expr, $b:expr $(,)?) => {
        symbolica::function!(
            $crate::network::library::symbolic::ETS.metric,
            $crate::symbolica_atom::IntoAtom::into_atom($a),
            $crate::symbolica_atom::IntoAtom::into_atom($b)
        )
    };
}

/// Builds a metric tensor atom `g(a,b)`.
#[macro_export]
macro_rules! g {
    ($a:expr, $b:expr $(,)?) => {
        $crate::metric!($a, $b)
    };
}

/// Builds a two-argument compact dot-product atom.
#[macro_export]
macro_rules! dot {
    ($a:expr, $b:expr $(,)?) => {
        symbolica::function!(
            $crate::network::tags::SPENSO_TAG.dot,
            $crate::symbolica_atom::IntoAtom::into_atom($a),
            $crate::symbolica_atom::IntoAtom::into_atom($b)
        )
    };
}

/// Forces an expression through the parser's pure-scalar boundary.
#[macro_export]
macro_rules! pure_scalar {
    ($expr:expr $(,)?) => {
        symbolica::function!(
            $crate::network::tags::SPENSO_TAG.pure_scalar,
            $crate::symbolica_atom::IntoAtom::into_atom($expr)
        )
    };
}

/// Builds a parser grouping over network factors.
#[macro_export]
macro_rules! bracket {
    ($($expr:expr),* $(,)?) => {
        symbolica::function!(
            $crate::network::tags::SPENSO_TAG.bracket,
            $($crate::symbolica_atom::IntoAtom::into_atom($expr)),*
        )
    };
}

/// Bundles structural slots inside one tensor argument.
#[macro_export]
macro_rules! aind {
    ($($slot:expr),* $(,)?) => {
        symbolica::function!(
            $crate::structure::abstract_index::AIND_SYMBOLS.aind,
            $($crate::symbolica_atom::IntoAtom::into_atom($slot)),*
        )
    };
}

/// Wraps a slot or representation atom in the dual-index marker.
#[macro_export]
macro_rules! dind {
    ($slot:expr $(,)?) => {
        symbolica::function!(
            $crate::structure::abstract_index::AIND_SYMBOLS.dind,
            $crate::symbolica_atom::IntoAtom::into_atom($slot)
        )
    };
}

/// Builds a generic chain or trace factor with explicit `in` and `out` markers.
///
/// Use the identifier tokens `in` and `out` in the argument list to insert the
/// chain placeholder symbols. Other arguments are converted through `IntoAtom`.
#[macro_export]
macro_rules! chain_factor {
    ($name:ident $(,)?) => {{
        let factor = symbolica::atom::FunctionBuilder::new($crate::tensor_symbol!($name));
        $crate::chain_factor!(@finish factor)
    }};
    ($name:ident, $($args:tt)+) => {{
        let factor = symbolica::atom::FunctionBuilder::new($crate::tensor_symbol!($name));
        $crate::chain_factor!(@finish factor, $($args)+)
    }};
    ($name:expr; $($args:tt)*) => {{
        let factor = symbolica::atom::FunctionBuilder::new($name);
        $crate::chain_factor!(@finish factor $(, $args)*)
    }};
    (@finish $factor:ident) => {
        $factor.finish()
    };
    (@finish $factor:ident, in $(, $rest:tt)*) => {{
        let $factor = $factor.add_arg(symbolica::atom::Atom::var($crate::network::tags::SPENSO_TAG.chain_in));
        $crate::chain_factor!(@finish $factor $(, $rest)*)
    }};
    (@finish $factor:ident, out $(, $rest:tt)*) => {{
        let $factor = $factor.add_arg(symbolica::atom::Atom::var($crate::network::tags::SPENSO_TAG.chain_out));
        $crate::chain_factor!(@finish $factor $(, $rest)*)
    }};
    (@finish $factor:ident, $arg:expr $(, $rest:tt)*) => {{
        let $factor = $factor.add_arg($crate::symbolica_atom::IntoAtom::into_atom($arg));
        $crate::chain_factor!(@finish $factor $(, $rest)*)
    }};
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

/// Builds an inert normalized cyclic symmetrizer over chain/trace factors.
///
/// The resulting `cyclic(...)` is canonicalized through Symbolica's
/// cyclesymmetric function attribute. Inside traces it is the compact syntax
/// for cyclic trace equivalence, for example `trace(rep, cyclic(a,b,c))`.
///
/// # Examples
///
/// ```ignore
/// use spenso::{cyclic, trace};
///
/// let cyclic_trace = trace!(rep, cyclic!(a, b, c));
/// let from_vec = cyclic!(; vec![a, b, c]);
/// ```
#[macro_export]
macro_rules! cyclic {
    (; $factors:expr $(,)?) => {
        $crate::symbolica_atom::cyclic(($factors).into_iter())
    };
    ($($factor:expr),* $(,)?) => {
        $crate::symbolica_atom::cyclic(vec![$($crate::symbolica_atom::IntoAtom::into_atom($factor)),*])
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
    ($start:expr, $end:expr $(,)?) => {
        $crate::network::tags::SPENSO_TAG.chain(
            $crate::symbolica_atom::IntoAtom::into_atom($start),
            $crate::symbolica_atom::IntoAtom::into_atom($end),
            std::iter::empty::<symbolica::atom::Atom>(),
        )
    };
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
/// The representation and factors are converted through [`IntoAtom`]. Empty
/// traces are emitted as `trace(rep)`. Non-empty traces are emitted in the
/// canonical cyclic form `trace(rep, cyclic(factors...))`.
///
/// # Examples
///
/// ```ignore
/// use spenso::{cyclic, trace};
///
/// let expr = trace!(
///     &cof_nc,
///     color_t_a,
///     color_t_b,
/// );
/// assert_eq!(expr, trace!(&cof_nc, cyclic!(color_t_a, color_t_b)));
///
/// let factors = vec![color_t_a, color_t_b];
/// let expr = trace!(&cof_nc; factors);
/// ```
#[macro_export]
macro_rules! trace {
    ($rep:expr $(,)?) => {
        $crate::symbolica_atom::trace(
            $crate::symbolica_atom::IntoAtom::into_atom($rep),
            std::iter::empty::<symbolica::atom::Atom>(),
        )
    };
    ($rep:expr; $factors:expr $(,)?) => {
        $crate::symbolica_atom::trace(
            $crate::symbolica_atom::IntoAtom::into_atom($rep),
            ($factors)
                .into_iter()
                .map($crate::symbolica_atom::IntoAtom::into_atom),
        )
    };
    ($rep:expr $(, $factor:expr)+; $factors:expr $(,)?) => {
        $crate::symbolica_atom::trace(
            $crate::symbolica_atom::IntoAtom::into_atom($rep),
            vec![$($crate::symbolica_atom::IntoAtom::into_atom($factor)),*]
                .into_iter()
                .chain(($factors).into_iter().map($crate::symbolica_atom::IntoAtom::into_atom)),
        )
    };
    ($rep:expr $(, $factor:expr)* $(,)?) => {
        $crate::symbolica_atom::trace(
            $crate::symbolica_atom::IntoAtom::into_atom($rep),
            vec![$($crate::symbolica_atom::IntoAtom::into_atom($factor)),*],
        )
    };
}

/// Builds a symbolic trace with a single symmetric projector payload.
///
/// This emits `trace(rep, sym(factors...))`. Use `trace!` for ordinary cyclic
/// closed traces.
///
/// # Examples
///
/// ```ignore
/// use spenso::trace_sym;
///
/// let invariant = trace_sym!(&cof_nc, color_t_a, color_t_b, color_t_c);
/// let from_vec = trace_sym!(&cof_nc; factors);
/// ```
#[macro_export]
macro_rules! trace_sym {
    ($rep:expr; $factors:expr $(,)?) => {
        $crate::symbolica_atom::trace_sym(
            $crate::symbolica_atom::IntoAtom::into_atom($rep),
            ($factors)
                .into_iter()
                .map($crate::symbolica_atom::IntoAtom::into_atom),
        )
    };
    ($rep:expr $(, $factor:expr)+; $factors:expr $(,)?) => {
        $crate::symbolica_atom::trace_sym(
            $crate::symbolica_atom::IntoAtom::into_atom($rep),
            vec![$($crate::symbolica_atom::IntoAtom::into_atom($factor)),*]
                .into_iter()
                .chain(($factors).into_iter().map($crate::symbolica_atom::IntoAtom::into_atom)),
        )
    };
    ($rep:expr $(, $factor:expr)+ $(,)?) => {
        $crate::symbolica_atom::trace_sym(
            $crate::symbolica_atom::IntoAtom::into_atom($rep),
            vec![$($crate::symbolica_atom::IntoAtom::into_atom($factor)),*],
        )
    };
}

#[cfg(test)]
mod tests {
    use crate::structure::representation::{LibraryRep, Minkowski};
    #[allow(unused_imports)]
    use crate::{
        aind, antisym, bracket, chain, chain_factor, cyclic, dind, dot, euc, g, lor, mink,
        network::tags::SPENSO_TAG,
        p, pure_scalar, q,
        shadowing::symbolica_utils::AtomCoreExt,
        sym,
        symbolica_atom::{ProjectorExpander, ScalarCollectExt, TensorCollectExt},
        trace, trace_sym,
    };
    use symbolica::{
        atom::{Atom, AtomView},
        symbol,
    };

    #[test]
    fn tensor_macros_create_tagged_heads() {
        #[allow(unused_imports)]
        use crate::{tensor, vector};

        let tensor = tensor!(tensor_macro_test_head);
        let vector = vector!(vector_macro_test_head, 1);

        let AtomView::Fun(tensor) = tensor.as_view() else {
            panic!("tensor macro should produce a function");
        };
        assert!(tensor.get_symbol().has_tag(&SPENSO_TAG.tensor));
        assert!(!tensor.get_symbol().has_tag(&SPENSO_TAG.rank1));

        let AtomView::Fun(vector) = vector.as_view() else {
            panic!("vector macro should produce a function");
        };
        assert!(vector.get_symbol().has_tag(&SPENSO_TAG.tensor));
        assert!(vector.get_symbol().has_tag(&SPENSO_TAG.rank1));
    }

    #[test]
    fn surface_macros_build_parser_syntax() {
        let mu = mink!(4, mu);
        let nu = mink!(4, nu);
        let p = p!(mink!(4));
        let q = q!(mink!(4));

        insta::assert_snapshot!(g!(mu.clone(), nu.clone()).to_bare_ordered_string(), @"g(mink(4,mu),mink(4,nu))");
        insta::assert_snapshot!(dot!(p.clone(), q.clone()).to_bare_ordered_string(), @"dot(p(mink(4)),q(mink(4)))");
        insta::assert_snapshot!(aind!(mu.clone(), nu.clone()).to_bare_ordered_string(), @"aind(mink(4,mu),mink(4,nu))");
        let dual = dind!(lor!(4, rho));
        let AtomView::Fun(dual) = dual.as_view() else {
            panic!("dind macro should produce a function");
        };
        assert_eq!(
            dual.get_symbol(),
            crate::structure::abstract_index::AIND_SYMBOLS.dind
        );
        insta::assert_snapshot!(pure_scalar!(Atom::num(1)).to_bare_ordered_string(), @"pure_scalar(1)");
        insta::assert_snapshot!(bracket!(p, q).to_bare_ordered_string(), @"bracket(p(mink(4)),q(mink(4)))");
        insta::assert_snapshot!(
            chain_factor!(factor, in, mu, out).to_bare_ordered_string(),
            @"factor(in,mink(4,mu),out)"
        );
    }

    #[test]
    fn collect_tensors_keeps_scalar_products_factored() {
        let (a, b, c, d) = symbol!("a", "b", "c", "d");
        let p = p!(mink!(4));
        let expr = (Atom::var(a) + Atom::var(b)) * (Atom::var(c) + Atom::var(d)) * p.clone()
            + Atom::var(a) * p;

        insta::assert_snapshot!(
            expr.collect_tensors().to_bare_ordered_string(),
            @"((a+b)*(c+d)+a)*p(mink(4))"
        );
    }

    #[test]
    fn collect_tensors_marks_chain_like_forms_as_maximal_factors() {
        let (a, b) = symbol!("a", "b");
        let mu = mink!(4, mu);
        let gamma_mu = chain_factor!(collect_test_gamma, in, mu, out);

        let chain_expr = chain!(mink!(4, i), mink!(4, j), gamma_mu.clone());
        insta::assert_snapshot!(
            (Atom::var(a) * chain_expr.clone() + Atom::var(b) * chain_expr)
                .collect_tensors()
                .to_bare_ordered_string(),
            @"(a+b)*chain(mink(4,i),mink(4,j),collect_test_gamma(in,mink(4,mu),out))"
        );

        let trace_expr = trace!(mink!(4), gamma_mu.clone());
        insta::assert_snapshot!(
            (Atom::var(a) * trace_expr.clone() + Atom::var(b) * trace_expr)
                .collect_tensors()
                .to_bare_ordered_string(),
            @"(a+b)*trace(mink(4),cyclic(collect_test_gamma(in,mink(4,mu),out)))"
        );

        let sym_expr = sym!(gamma_mu.clone());
        insta::assert_snapshot!(
            (Atom::var(a) * sym_expr.clone() + Atom::var(b) * sym_expr)
                .collect_tensors()
                .to_bare_ordered_string(),
            @"(a+b)*sym(collect_test_gamma(in,mink(4,mu),out))"
        );

        let cyclic_expr = cyclic!(gamma_mu);
        insta::assert_snapshot!(
            (Atom::var(a) * cyclic_expr.clone() + Atom::var(b) * cyclic_expr)
                .collect_tensors()
                .to_bare_ordered_string(),
            @"(a+b)*cyclic(collect_test_gamma(in,mink(4,mu),out))"
        );
    }

    #[test]
    fn collect_rep_only_wraps_matching_representations() {
        let (a, b) = symbol!("a", "b");
        let mink_p = p!(mink!(4));
        let euc_q = q!(euc!(3));
        let expr = Atom::var(a) * mink_p.clone()
            + Atom::var(b) * mink_p
            + Atom::var(a) * euc_q.clone()
            + Atom::var(b) * euc_q;

        insta::assert_snapshot!(
            expr.collect_rep(LibraryRep::from(Minkowski {}))
                .to_bare_ordered_string(),
            @"(a+b)*p(mink(4))+a*q(euc(3))+b*q(euc(3))"
        );
    }

    #[test]
    fn collect_metrics_only_wraps_metric_heads() {
        let (a, b) = symbol!("a", "b");
        let metric = g!(mink!(4, mu), mink!(4, nu));
        let vector = p!(mink!(4));
        let expr = Atom::var(a) * metric.clone()
            + Atom::var(b) * metric
            + Atom::var(a) * vector.clone()
            + Atom::var(b) * vector;

        insta::assert_snapshot!(
            expr.collect_metrics().to_bare_ordered_string(),
            @"(a+b)*g(mink(4,mu),mink(4,nu))+a*p(mink(4))+b*p(mink(4))"
        );
    }

    #[test]
    fn collect_scalar_symbols_combines_nested_coefficients() {
        let (a, b, c) = symbol!("a", "b", "c");
        let expr = ((-Atom::var(a) + Atom::var(b)) * Atom::var(c) - Atom::var(b) * Atom::var(c))
            * p!(mink!(4));

        insta::assert_snapshot!(
            expr.collect_tensors()
                .collect_scalar_symbols()
                .to_bare_ordered_string(),
            @"-1*a*c*p(mink(4))"
        );
    }

    #[test]
    fn chain_macro_accepts_iterable_factors() {
        let start = Atom::var(symbol!("start"));
        let end = Atom::var(symbol!("end"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        assert_eq!(
            chain!(start.clone(), end.clone(); vec![first.clone(), second.clone()]),
            chain!(start, end, first, second)
        );
    }

    #[test]
    fn trace_macro_accepts_iterable_factors() {
        let rep = Atom::var(symbol!("rep"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        assert_eq!(
            trace!(rep.clone(); vec![first.clone(), second.clone()]),
            trace!(rep.clone(), first.clone(), second.clone())
        );
        assert_eq!(
            trace!(rep, first.clone(), second.clone()),
            trace!(Atom::var(symbol!("rep")), cyclic!(first, second))
        );
    }

    #[test]
    fn sym_macro_canonicalizes_arguments() {
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        assert_eq!(sym!(second.clone(), first.clone()), sym!(first, second));
    }

    #[test]
    fn cyclic_macro_canonicalizes_rotations() {
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));
        let third = Atom::var(symbol!("third"));

        assert_eq!(
            cyclic!(second.clone(), third.clone(), first.clone()),
            cyclic!(first.clone(), second.clone(), third.clone())
        );
        assert_ne!(
            cyclic!(first.clone(), third.clone(), second.clone()),
            cyclic!(second, third, first)
        );
    }

    #[test]
    fn cyclic_unwraps_single_symmetric_projector() {
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        assert_eq!(
            cyclic!(sym!(first.clone(), second.clone())),
            sym!(first, second)
        );
    }

    #[test]
    fn trace_sym_macro_builds_symmetric_trace() {
        let rep = Atom::var(symbol!("rep"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        assert_eq!(
            trace_sym!(rep.clone(), first.clone(), second.clone()),
            trace!(rep, cyclic!(sym!(first, second)))
        );
    }

    #[test]
    fn expand_sym_in_chain_preserves_endpoints() {
        let start = Atom::var(symbol!("start"));
        let end = Atom::var(symbol!("end"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        let expanded = chain!(
            start.clone(),
            end.clone(),
            sym!(first.clone(), second.clone())
        )
        .expand_projectors();
        let expected = Atom::num(1) / Atom::num(2)
            * chain!(start.clone(), end.clone(), first.clone(), second.clone())
            + Atom::num(1) / Atom::num(2) * chain!(start, end, second, first);

        assert_eq!(expanded, expected);
    }

    #[test]
    fn expand_antisym_in_trace_uses_cyclicity() {
        let rep = Atom::var(symbol!("rep"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        let expanded =
            trace!(rep.clone(), antisym!(first.clone(), second.clone())).expand_projectors();

        assert_eq!(expanded, Atom::Zero);
    }

    #[test]
    fn expand_sym_in_trace_emits_cyclic_projectors() {
        let rep = Atom::var(symbol!("rep"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));
        let third = Atom::var(symbol!("third"));

        let expanded = trace!(
            rep.clone(),
            sym!(first.clone(), second.clone(), third.clone())
        )
        .expand_projectors();
        let expected = Atom::num(1) / Atom::num(2)
            * trace!(
                rep.clone(),
                cyclic!(first.clone(), second.clone(), third.clone())
            )
            + Atom::num(1) / Atom::num(2) * trace!(rep, cyclic!(first, third, second));

        assert_eq!(expanded, expected);
    }

    #[test]
    fn expand_cyclic_in_chain_rotates_factors() {
        let start = Atom::var(symbol!("start"));
        let end = Atom::var(symbol!("end"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));
        let third = Atom::var(symbol!("third"));

        let expanded = chain!(
            start.clone(),
            end.clone(),
            cyclic!(first.clone(), second.clone(), third.clone())
        )
        .expand_projectors();
        let expected = Atom::num(1) / Atom::num(3)
            * chain!(
                start.clone(),
                end.clone(),
                first.clone(),
                second.clone(),
                third.clone()
            )
            + Atom::num(1) / Atom::num(3)
                * chain!(
                    start.clone(),
                    end.clone(),
                    second.clone(),
                    third.clone(),
                    first.clone()
                )
            + Atom::num(1) / Atom::num(3) * chain!(start, end, third, first, second);

        assert_eq!(expanded, expected);
    }

    #[test]
    fn expand_cyclic_in_trace_keeps_compact_cyclic_trace() {
        let rep = Atom::var(symbol!("rep"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));
        let third = Atom::var(symbol!("third"));

        let expanded = trace!(
            rep.clone(),
            cyclic!(first.clone(), second.clone(), third.clone())
        )
        .expand_projectors();
        let expected = trace!(rep, cyclic!(first, second, third));

        assert_eq!(expanded, expected);
    }

    #[test]
    fn expand_antisym_in_chain_handles_canonicalization_sign() {
        let start = Atom::var(symbol!("start"));
        let end = Atom::var(symbol!("end"));
        let first = Atom::var(symbol!("first"));
        let second = Atom::var(symbol!("second"));

        let expanded = chain!(
            start.clone(),
            end.clone(),
            antisym!(second.clone(), first.clone())
        )
        .expand_projectors();
        let expected = Atom::num(-1) / Atom::num(2)
            * chain!(start.clone(), end.clone(), first.clone(), second.clone())
            + Atom::num(1) / Atom::num(2) * chain!(start, end, second, first);

        assert_eq!(expanded, expected);
    }
}
