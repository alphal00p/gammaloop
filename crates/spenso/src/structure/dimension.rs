use derive_more::Display;
use serde::{Deserialize, Serialize};

#[cfg(feature = "shadowing")]
use symbolica::{atom::AtomView, coefficient::CoefficientView};
use thiserror::Error;

#[cfg(feature = "shadowing")]
use symbolica::atom::{Atom, Symbol};

#[cfg(feature = "shadowing")]
use crate::shadowing::symbolica_utils::SerializableSymbol;

/// A Dimension
#[derive(
    Debug,
    Copy,
    Clone,
    Eq,
    PartialEq,
    Hash,
    Display,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    // bincode_trait_derive::BorrowDecodeFromDecode,
)]
#[cfg_attr(
    feature = "shadowing",
    trait_decode(trait = symbolica::state::HasStateMap),
)]
pub enum Dimension {
    Concrete(usize),
    #[cfg(feature = "shadowing")]
    Symbolic(SerializableSymbol),
}

impl Dimension {
    pub fn new_concrete(value: usize) -> Self {
        Self::Concrete(value)
    }

    #[cfg(feature = "shadowing")]
    pub fn to_symbolic(&self) -> Atom {
        match self {
            Self::Concrete(c) => Atom::num(*c as i64),
            Self::Symbolic(s) => Atom::var((*s).into()),
        }
    }
}

impl Ord for Dimension {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self, other) {
            (Dimension::Concrete(s), Dimension::Concrete(o)) => s.cmp(o),
            #[cfg(feature = "shadowing")]
            _ => std::cmp::Ordering::Equal,
        }
    }
}

impl PartialOrd for Dimension {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Error, Debug)]
pub enum DimensionError {
    #[error("Dimension is not concrete")]
    NotConcrete,
    #[error("Dimension too large")]
    TooLarge,
    #[error("Dimension negative")]
    Negative,
    #[error("Dimension is not natural")]
    NotNatural,

    #[cfg(feature = "shadowing")]
    #[error("Dimension is not a var or int: {0}")]
    NotVarOrInt(Atom),
}

#[allow(unreachable_patterns)]
impl TryFrom<Dimension> for usize {
    type Error = DimensionError;

    fn try_from(value: Dimension) -> Result<Self, Self::Error> {
        match value {
            Dimension::Concrete(c) => Ok(c),
            _ => Err(DimensionError::NotConcrete),
        }
    }
}

impl From<usize> for Dimension {
    fn from(value: usize) -> Self {
        Dimension::Concrete(value)
    }
}

#[cfg(feature = "shadowing")]
impl From<Symbol> for Dimension {
    fn from(value: Symbol) -> Self {
        Dimension::Symbolic(value.into())
    }
}

#[cfg(feature = "shadowing")]
impl<'a> TryFrom<AtomView<'a>> for Dimension {
    type Error = DimensionError;
    fn try_from(value: AtomView<'a>) -> Result<Self, Self::Error> {
        match value {
            AtomView::Var(a) => Ok(Dimension::Symbolic(a.get_symbol().into())),
            AtomView::Num(n) => match n.get_coeff_view() {
                CoefficientView::Natural(n, 1, _, _) => {
                    if n < 0 {
                        return Err(DimensionError::Negative);
                    }
                    Ok(Dimension::Concrete(n as usize))
                }
                _ => {
                    let mut err = Atom::new();
                    err.set_from_view(&value);
                    Err(DimensionError::NotVarOrInt(err))
                }
            },
            _ => {
                let mut err = Atom::new();
                err.set_from_view(&value);
                Err(DimensionError::NotVarOrInt(err))
            }
        }
    }
}

#[allow(unreachable_patterns)]
impl PartialEq<usize> for Dimension {
    fn eq(&self, other: &usize) -> bool {
        match self {
            Self::Concrete(c) => c == other,
            _ => false,
        }
    }
}

#[allow(unreachable_patterns)]
impl PartialEq<Dimension> for usize {
    fn eq(&self, other: &Dimension) -> bool {
        match other {
            Dimension::Concrete(c) => c == self,
            _ => false,
        }
    }
}

#[cfg(test)]
#[cfg(feature = "shadowing")]
mod shadowing_tests {
    use symbolica::{atom::Atom, function, symbol};

    use super::Dimension;

    #[test]
    fn dimension_from_view() {
        let a = Atom::num(5);
        let b = Atom::var(symbol!("b"));
        let c = function!(symbol!("a"), symbol!("b"));
        let d = Atom::num(-1);
        let e = Atom::num((1, 2));

        let dima = Dimension::try_from(a.as_view()).unwrap();

        assert_eq!(dima, Dimension::Concrete(5));
        let dimb = Dimension::try_from(b.as_view()).unwrap();
        assert_eq!(dimb, Dimension::Symbolic(symbol!("b").into()));
        let dimc = Dimension::try_from(c.as_view());
        assert!(dimc.is_err());
        let dimd = Dimension::try_from(d.as_view());
        assert!(dimd.is_err());
        let dime = Dimension::try_from(e.as_view());
        assert!(dime.is_err());
    }
}
