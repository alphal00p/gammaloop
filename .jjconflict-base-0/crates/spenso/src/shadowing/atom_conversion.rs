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
