use derive_more::Display;
use serde::{Deserialize, Serialize};

#[cfg(feature = "symbolica")]
use symbolica::atom::{Atom, AtomView, Symbol};
#[cfg(feature = "symbolica")]
use symbolica::coefficient::CoefficientView;
#[cfg(feature = "symbolica")]
use symbolica_utils::{InternedAtomId, SerializableSymbol};
use thiserror::Error;

/// The dimension of a representation.
///
/// Variants are ordered `Concrete < Symbolic < Expr`, comparing payloads
/// within a variant, which gives a total order consistent with equality.
///
/// ```
/// # #[cfg(feature = "symbolica")]
/// # {
/// use spenso::structure::dimension::Dimension;
/// use symbolica::parse;
///
/// let adjoint = parse!("dimension_doc::nc^2 - 1");
/// let dimension = Dimension::try_from(adjoint.as_view()).unwrap();
///
/// assert!(matches!(dimension, Dimension::Expr(_)));
/// assert_eq!(dimension.to_symbolic(), adjoint);
/// assert!(Dimension::Concrete(4) < dimension);
/// # }
/// ```
#[derive(
    Debug,
    Copy,
    Clone,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
    Hash,
    Display,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    // bincode_trait_derive::BorrowDecodeFromDecode,
)]
#[cfg_attr(
    feature = "symbolica",
    trait_decode(trait = symbolica::state::HasStateMap),
)]
pub enum Dimension {
    Concrete(usize),
    #[cfg(feature = "symbolica")]
    Symbolic(SerializableSymbol),
    /// A composite symbolic dimension such as `Nc^2 - 1` or `4 - 2*eps`.
    #[cfg(feature = "symbolica")]
    Expr(InternedAtomId),
}

impl Dimension {
    pub fn new_concrete(value: usize) -> Self {
        Self::Concrete(value)
    }

    #[cfg(feature = "symbolica")]
    pub fn to_symbolic(&self) -> Atom {
        match self {
            Self::Concrete(c) => Atom::num(*c as i64),
            Self::Symbolic(s) => Atom::var((*s).into()),
            Self::Expr(expression) => expression.to_atom(),
        }
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

    #[cfg(feature = "symbolica")]
    #[error("Dimension is not a natural number, variable or composite expression: {0}")]
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

#[cfg(feature = "symbolica")]
impl From<Symbol> for Dimension {
    fn from(value: Symbol) -> Self {
        Dimension::Symbolic(value.into())
    }
}

#[cfg(feature = "symbolica")]
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
            AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_) => {
                Ok(Dimension::Expr(InternedAtomId::intern(value)))
            }
            // A function argument in the dimension position marks a wrapped
            // representation instead, such as `dind(cof(3, i))`.
            AtomView::Fun(_) => Err(DimensionError::NotVarOrInt(value.to_owned())),
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
#[cfg(feature = "symbolica")]
mod shadowing_tests {
    use std::cmp::Ordering;
    use std::collections::BTreeSet;

    use symbolica::atom::Atom;
    use symbolica::{function, parse, symbol};

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

    #[test]
    fn composite_dimensions_round_trip() {
        let adjoint = parse!("dimension_test::nc^2 - 1");
        let regulated = parse!("4 - 2*dimension_test::eps");

        let dimension = Dimension::try_from(adjoint.as_view()).unwrap();
        assert!(matches!(dimension, Dimension::Expr(_)));
        assert_eq!(dimension.to_symbolic(), adjoint);
        assert_eq!(dimension, Dimension::try_from(adjoint.as_view()).unwrap());
        assert_ne!(dimension, Dimension::try_from(regulated.as_view()).unwrap());
    }

    #[test]
    fn ordering_is_total_and_consistent_with_equality() {
        let dimensions = [
            Dimension::Concrete(3),
            Dimension::Concrete(4),
            Dimension::Symbolic(symbol!("dimension_test::a").into()),
            Dimension::Symbolic(symbol!("dimension_test::b").into()),
            Dimension::try_from(parse!("dimension_test::nc^2 - 1").as_view()).unwrap(),
            Dimension::try_from(parse!("4 - 2*dimension_test::eps").as_view()).unwrap(),
        ];

        for left in &dimensions {
            for right in &dimensions {
                assert_eq!(left.cmp(right) == Ordering::Equal, left == right);
                assert_eq!(left.cmp(right), right.cmp(left).reverse());
            }
        }
        assert!(dimensions[1] < dimensions[2]);
        assert!(dimensions[3] < dimensions[4]);
        assert_eq!(BTreeSet::from(dimensions).len(), dimensions.len());
    }
}
