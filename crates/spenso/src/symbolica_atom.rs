use crate::structure::{
    representation::{RepName, Representation},
    slot::{AbsInd, IsAbstractSlot, ParseableAind, Slot},
};
use symbolica::atom::{Atom, AtomOrView, AtomView, Symbol};

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
}
