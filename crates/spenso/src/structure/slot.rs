use super::{
    abstract_index::{AbstractIndex, AbstractIndexError},
    dimension::DimensionError,
    representation::{
        BaseRepName, LibraryRep, LibrarySlot, RepName, Representation, RepresentationError,
    },
};
use crate::structure::dimension::Dimension;
use bincode::Encode;
// #[cfg(feature = "shadowing")]
use serde::{Deserialize, Serialize};
use std::hash::Hash;
use std::{
    cmp::Ordering,
    fmt::{Debug, Display},
};
#[cfg(feature = "shadowing")]
use symbolica::{
    atom::{Atom, AtomView, Symbol},
    symbol,
};

use thiserror::Error;

#[derive(
    Debug,
    Copy,
    Clone,
    Ord,
    PartialOrd,
    Eq,
    PartialEq,
    Hash,
    Serialize,
    Deserialize,
    Encode,
    bincode_trait_derive::Decode,
    // bincode_trait_derive::BorrowDecodeFromDecode,
)]
#[cfg_attr(
    feature = "shadowing",
    trait_decode(trait = symbolica::state::HasStateMap),
)]
/// A [`Slot`] is an index, identified by a `usize` and a [`Representation`].
///
/// A vector of slots thus identifies the shape and type of the tensor.
/// Two indices are considered matching if *both* the `Dimension` and the [`Representation`] matches.
///
/// # Example
///
/// It can be built from a `Representation` calling one of the built in representations e.g.
/// Or one can define custom representations{}
pub struct Slot<T: RepName, Aind = AbstractIndex> {
    pub(crate) rep: Representation<T>,
    pub aind: Aind,
}

impl<T: RepName, Aind> Slot<T, Aind> {
    pub fn cast<U: RepName + From<T>>(self) -> Slot<U, Aind> {
        Slot {
            aind: self.aind,
            rep: Representation {
                dim: self.rep.dim,
                rep: U::from(self.rep.rep),
            },
        }
    }
}

#[derive(Error, Debug)]
pub enum SlotError {
    #[error("Dimension is not concrete")]
    NotConcrete,
    #[error("Empty structure")]
    EmptyStructure,
    #[error("Argument is not a natural number")]
    NotNatural,
    #[error("Abstract index error :{0}")]
    AindError(#[from] AbstractIndexError),
    #[error("Representation error :{0}")]
    RepError(#[from] RepresentationError),
    #[error("Argument is not a number")]
    NotNumber,
    #[error("No more arguments")]
    NoMoreArguments,
    #[error("Too many arguments")]
    TooManyArguments,
    #[error("Not a slot, isn't a representation")]
    NotRepresentation,
    #[error("Not a slot, is composite")]
    Composite,
    #[error("{0}")]
    DimErr(#[from] DimensionError),
    #[error("{0}")]
    Any(#[from] eyre::Error),
}

#[cfg(feature = "shadowing")]
/// Can possibly constuct a Slot from an `AtomView`, if it is of the form: <representation>(<dimension>,<index>)
///
impl<'a, T: RepName, Aind> TryFrom<AtomView<'a>> for Slot<T, Aind>
where
    Aind: ParseableAind,
{
    type Error = SlotError;

    fn try_from(value: AtomView<'a>) -> Result<Self, Self::Error> {
        let (rep, mut iter) = if let AtomView::Fun(f) = value {
            let name = f.get_symbol();

            let innerf = f.iter().next().ok_or(SlotError::Composite)?;

            if let AtomView::Fun(innerf) = innerf {
                let rep =
                    T::try_from_symbol(innerf.get_symbol(), name).map_err(SlotError::RepError)?;

                (rep, innerf.iter())
            } else {
                let rep = T::try_from_symbol_coerced(name).map_err(SlotError::RepError)?;
                (rep, f.iter())
            }
        } else {
            return Err(SlotError::Composite);
        };

        let dim: Dimension = if let Some(a) = iter.next() {
            Dimension::try_from(a).map_err(SlotError::DimErr)?
        } else {
            return Err(SlotError::NoMoreArguments);
        };

        let index: Aind = if let Some(a) = iter.next() {
            Ok(Aind::from_view(a).map_err(|a| a.into())?)
        } else {
            Err(SlotError::NoMoreArguments)
        }?;

        if if let Some(a) = iter.next() {
            Ok(Aind::from_view(a).map_err(|a| a.into())?)
        } else {
            Err(SlotError::NoMoreArguments)
        }
        .is_ok()
        {
            return Err(SlotError::TooManyArguments);
        }

        Ok(Slot {
            rep: Representation { dim, rep },
            aind: index,
        })
    }
}

// pub trait SlotFromRep<S:IsAbstractSlot>:Rep{

// }

pub trait ConstructibleSlot<T: RepName, Aind> {
    fn new(rep: T, dim: Dimension, aind: Aind) -> Self;
}

impl<T: BaseRepName, Aind> ConstructibleSlot<T, Aind> for Slot<T, Aind> {
    fn new(_: T, dim: Dimension, aind: Aind) -> Self {
        Slot {
            aind,
            rep: Representation {
                dim,
                rep: T::default(),
            },
        }
    }
}

pub trait AbsInd:
    Copy + PartialEq + Eq + Debug + Clone + Hash + Ord + Display + Send + Sync + 'static
{
}

#[cfg(feature = "shadowing")]
pub trait ParseableAind: Sized {
    type Error: Into<SlotError>;
    fn from_view(view: AtomView<'_>) -> Result<Self, Self::Error>;

    fn to_atom(&self) -> Atom;
}

pub trait DummyAind {
    fn new_dummy() -> Self;
    fn new_dummy_at(i: usize) -> Self;
    fn is_dummy(&self) -> bool;
}

pub trait IsAbstractSlot: Copy + PartialEq + Eq + Debug + Clone + Hash + Ord + Display {
    type Aind: AbsInd;
    type R: RepName;

    fn reindex(self, id: Self::Aind) -> Self;
    fn dim(&self) -> Dimension;
    fn to_dummy_rep(&self) -> LibrarySlot<Self::Aind> {
        let rep = self.rep().to_dummy().to_lib();
        let aind = self.aind();
        Slot { rep, aind }
    }

    fn to_dummy_ind(&self) -> LibrarySlot<Self::Aind>
    where
        Self::Aind: DummyAind,
    {
        Slot {
            rep: self.rep().to_lib(),
            aind: Self::Aind::new_dummy(),
        }
    }

    fn to_lib(&self) -> LibrarySlot<Self::Aind> {
        let rep: LibraryRep = self.rep_name().into();
        rep.new_slot(self.dim(), self.aind())
    }
    fn aind(&self) -> Self::Aind;
    fn set_aind(&mut self, aind: Self::Aind);
    fn rep_name(&self) -> Self::R;
    fn rep(&self) -> Representation<Self::R> {
        Representation {
            dim: self.dim(),
            rep: self.rep_name(),
        }
    }

    #[cfg(feature = "shadowing")]
    /// using the function builder of the representation add the abstract index as an argument, and finish it to an Atom.
    fn to_atom(&self) -> Atom
    where
        Self::Aind: ParseableAind;
    #[cfg(feature = "shadowing")]
    fn to_symbolic_wrapped(&self) -> Atom
    where
        Self::Aind: ParseableAind;
    // #[cfg(feature = "shadowing")]
    // fn try_from_view<'a>(v: AtomView<'a>) -> Result<Self, SlotError>
    // where
    //     Self::Aind: TryFrom<AtomView<'a>>,
    //     SlotError: From<<Self::Aind as TryFrom<AtomView<'a>>>::Error>;
}

pub trait DualSlotTo: IsAbstractSlot {
    type Dual: IsAbstractSlot;
    fn dual(&self) -> Self::Dual;
    fn matches(&self, other: &Self::Dual) -> bool;

    fn match_cmp(&self, other: &Self::Dual) -> Ordering;
}

impl<T: RepName, A: AbsInd> IsAbstractSlot for Slot<T, A> {
    type Aind = A;
    type R = T;
    // type Dual = GenSlot<T::Dual>;
    fn dim(&self) -> Dimension {
        self.rep.dim
    }

    fn reindex(mut self, id: Self::Aind) -> Self {
        self.aind = id;
        self
    }
    fn aind(&self) -> Self::Aind {
        self.aind
    }
    fn rep_name(&self) -> Self::R {
        self.rep.rep
    }

    fn set_aind(&mut self, aind: Self::Aind) {
        self.aind = aind;
    }
    #[cfg(feature = "shadowing")]
    fn to_atom(&self) -> Atom
    where
        Self::Aind: ParseableAind,
    {
        self.rep.to_symbolic([self.aind.to_atom()])
    }
    #[cfg(feature = "shadowing")]
    fn to_symbolic_wrapped(&self) -> Atom
    where
        Self::Aind: ParseableAind,
    {
        use symbolica::function;

        self.rep
            .to_symbolic([function!(symbol!("indexid"), self.aind.to_atom())])
    }
    // #[cfg(feature = "shadowing")]
    // fn try_from_view<'a>(v: AtomView<'a>) -> Result<Self, SlotError>
    // where
    //     Self::Aind: TryFrom<AtomView<'a>>,
    //     SlotError: From<<Self::Aind as TryFrom<AtomView<'a>>>::Error>,
    // {
    //     Slot::try_from(v)
    // }
}

impl<T: RepName, Aind: AbsInd> std::fmt::Display for Slot<T, Aind> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.rep.rep.is_self_dual() {
            write!(f, "{}|{}", self.rep, self.aind)
        } else if self.rep.rep.is_dual() {
            write!(f, "{}|{:-}", self.rep, self.aind)
        } else {
            write!(f, "{}|{:+}", self.rep, self.aind)
        }
    }
}

#[cfg(feature = "shadowing")]
impl<T: RepName, Aind: AbsInd> Slot<T, Aind>
where
    Atom: From<Aind>,
{
    pub fn to_pattern(&self, dimension: Symbol) -> Atom {
        self.rep
            .rep
            .to_symbolic([Atom::var(dimension), Atom::from(self.aind)])
    }
}

impl<T: RepName, Aind: AbsInd> DualSlotTo for Slot<T, Aind> {
    type Dual = Slot<T::Dual, Aind>;
    fn dual(&self) -> Slot<T::Dual, Aind> {
        Slot {
            rep: self.rep.dual(),
            aind: self.aind,
        }
    }
    fn matches(&self, other: &Self::Dual) -> bool {
        self.rep.matches(&other.rep) && self.aind() == other.aind()
    }

    fn match_cmp(&self, other: &Self::Dual) -> Ordering {
        self.rep
            .match_cmp(&other.rep)
            .then(self.aind.cmp(&other.aind))
    }
}

#[cfg(test)]
#[cfg(feature = "shadowing")]
mod shadowing_tests {
    use insta::assert_snapshot;
    use symbolica::{atom::AtomCore, parse, symbol};

    use crate::structure::{
        representation::{DualLorentz, LibraryRep, Lorentz, RepName, Representation},
        slot::{DualSlotTo, IsAbstractSlot, Slot},
    };

    #[test]
    fn doc_slot() {
        let mink: Representation<Lorentz> = Lorentz {}.new_rep(4);

        let mud: Slot<Lorentz> = mink.slot(0);
        let muu: Slot<DualLorentz> = mink.slot(0).dual();

        assert!(mud.matches(&muu));
        assert_eq!("lor🠓4|₀", format!("{muu}"));

        let custom_mink = LibraryRep::new_dual("custom_lor").unwrap();

        let nud: Slot<LibraryRep> = custom_mink.new_slot(4, 0);
        let nuu: Slot<LibraryRep> = nud.dual();

        assert!(nuu.matches(&nud));
        assert_eq!("custom_lor🠓4|₀", format!("{nuu}"));
    }

    #[test]
    fn to_symbolic() {
        let mink = Lorentz {}.new_rep(4);
        let mu: Slot<Lorentz> = mink.slot(0);
        println!("{}", mu.to_atom());
        assert_snapshot!( mu.to_atom().to_canonical_string(),@"spenso::{spenso::upper}::lor(4,0)");
        // assert_eq!("lor🠑4|₀", mu.dual().to_string());

        let mink = Lorentz {}.new_rep(4);
        let mu: Slot<Lorentz> = mink.slot(0);
        let atom = mu.to_atom();
        let slot = Slot::try_from(atom.as_view()).unwrap();
        assert_eq!(slot, mu);
    }

    #[test]
    fn slot_from_atom_view() {
        let mink = Lorentz {}.new_rep(4);
        let mu = mink.slot(0);
        let atom = mu.to_atom();
        assert_eq!(Slot::try_from(atom.as_view()).unwrap(), mu);
        assert_eq!(
            Slot::<Lorentz>::try_from(atom.as_view()).unwrap().dual(),
            mu.dual()
        );
        assert_eq!(
            Slot::try_from(mu.dual().to_atom().as_view()).unwrap(),
            mu.dual()
        );

        let expr = parse!("dind(lor(4,-1))");

        let _slot: Slot<LibraryRep> = Slot::try_from(expr.as_view()).unwrap();
        let _slot: Slot<DualLorentz> = Slot::try_from(expr.as_view()).unwrap();

        println!("{}", _slot.to_symbolic_wrapped());
        println!("{}", _slot.to_pattern(symbol!("d_")));
    }
}
