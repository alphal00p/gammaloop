use std::sync::atomic::{AtomicUsize, Ordering};

use linnet::half_edge::{
    involution::{EdgeIndex, Hedge},
    NodeIndex,
};
use serde::{Deserialize, Serialize};
use spenso::{
    structure::slot::{AbsInd, DummyAind, SlotError},
    utils::{to_subscript, to_superscript},
};
use symbolica::{
    atom::{Atom, AtomView},
    coefficient::CoefficientView,
    function, symbol,
};
use thiserror::Error;

static DUMMYCOUNTER: AtomicUsize = AtomicUsize::new(0);
/// A type that represents the name of an index in a tensor.
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
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    // bincode_trait_derive::BorrowDecodeFromDecode,
)]
pub enum Aind {
    Normal(usize),
    Hedge(u16, u16),
    Edge(u16, u16),
    Vertex(u16, u16),
    Dummy(usize),
}

pub trait NewAind {
    fn aind(self, local: u16) -> Aind;
}

impl NewAind for EdgeIndex {
    fn aind(self, local: u16) -> Aind {
        Aind::Edge(usize::from(self) as u16, local)
    }
}

impl NewAind for NodeIndex {
    fn aind(self, local: u16) -> Aind {
        Aind::Vertex(usize::from(self) as u16, local)
    }
}
impl NewAind for Hedge {
    fn aind(self, local: u16) -> Aind {
        Aind::Hedge(self.0 as u16, local)
    }
}

impl AbsInd for Aind {}

impl DummyAind for Aind {
    fn new_dummy() -> Self {
        Aind::Dummy(DUMMYCOUNTER.fetch_add(1, Ordering::Relaxed))
    }
}

impl std::fmt::Display for Aind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Aind::Normal(v) => {
                if f.sign_minus() {
                    write!(f, "{}", to_subscript(*v as isize))
                } else if f.sign_plus() {
                    write!(f, "{}", to_superscript(*v as isize))
                } else {
                    write!(f, "{}", v)
                }
            }
            Aind::Hedge(i, j) => {
                if f.sign_minus() {
                    write!(
                        f,
                        "h{}.{}",
                        to_subscript(*i as isize),
                        to_subscript(*j as isize)
                    )
                } else if f.sign_plus() {
                    write!(
                        f,
                        "h{}'{}",
                        to_superscript(*i as isize),
                        to_superscript(*j as isize)
                    )
                } else {
                    write!(f, "h{}-{}", i, j)
                }
            }
            Aind::Edge(i, j) => {
                if f.sign_minus() {
                    write!(
                        f,
                        "e{}.{}",
                        to_subscript(*i as isize),
                        to_subscript(*j as isize)
                    )
                } else if f.sign_plus() {
                    write!(
                        f,
                        "e{}'{}",
                        to_superscript(*i as isize),
                        to_superscript(*j as isize)
                    )
                } else {
                    write!(f, "e{}-{}", i, j)
                }
            }
            Aind::Vertex(i, j) => {
                if f.sign_minus() {
                    write!(
                        f,
                        "v{}.{}",
                        to_subscript(*i as isize),
                        to_subscript(*j as isize)
                    )
                } else if f.sign_plus() {
                    write!(
                        f,
                        "v{}'{}",
                        to_superscript(*i as isize),
                        to_superscript(*j as isize)
                    )
                } else {
                    write!(f, "v{}-{}", i, j)
                }
            }

            Aind::Dummy(v) => {
                if f.sign_minus() {
                    write!(f, "d{}", to_subscript(*v as isize))
                } else if f.sign_plus() {
                    write!(f, "d{}", to_superscript(*v as isize))
                } else {
                    write!(f, "d{}", v)
                }
            }
        }
    }
}

impl From<usize> for Aind {
    fn from(value: usize) -> Self {
        Aind::Normal(value)
    }
}

impl From<Aind> for Atom {
    fn from(value: Aind) -> Self {
        match value {
            Aind::Dummy(i) => function!(symbol!("dummy"), i as i64),
            Aind::Normal(i) => Atom::num(i as i64),
            Aind::Edge(i, j) => {
                if j != 0 {
                    function!(symbol!("edge"), i as i64, j as i64)
                } else {
                    function!(symbol!("edge"), i as i64)
                }
            }
            Aind::Hedge(i, j) => {
                if j != 0 {
                    function!(symbol!("hedge"), i as i64, j as i64)
                } else {
                    function!(symbol!("hedge"), i as i64)
                }
            }
            Aind::Vertex(i, j) => {
                if j != 0 {
                    function!(symbol!("vertex"), i as i64, j as i64)
                } else {
                    function!(symbol!("vertex"), i as i64)
                }
            }
        }
    }
}
impl<'a> From<Aind> for symbolica::atom::AtomOrView<'a> {
    fn from(value: Aind) -> Self {
        symbolica::atom::AtomOrView::Atom(Atom::from(value))
    }
}

#[derive(Error, Debug)]
pub enum AindError {
    #[error("Argument is not a natural number")]
    NotNatural,
    #[error("Argument  {0} is not a valid index")]
    NotIndex(String),
    #[error("parsing error")]
    ParsingError(String),
}

impl From<AindError> for SlotError {
    fn from(value: AindError) -> Self {
        SlotError::Any(value.into())
    }
}

impl TryFrom<AtomView<'_>> for Aind {
    type Error = AindError;

    fn try_from(view: AtomView<'_>) -> Result<Self, Self::Error> {
        match view {
            AtomView::Num(n) => match n.get_coeff_view() {
                CoefficientView::Natural(n, 1, _, _) => Ok(Aind::Normal(
                    usize::try_from(n).map_err(|e| AindError::NotIndex(e.to_string()))?,
                )),
                _ => Err(AindError::NotNatural),
            },
            AtomView::Fun(f) => {
                if f.get_symbol() == symbol!("dummy") {
                    if f.get_nargs() == 1 {
                        let arg = f.iter().next().unwrap();

                        if let Ok(a) = i64::try_from(arg) {
                            if a >= 0 {
                                DUMMYCOUNTER.fetch_max(a as usize + 1, Ordering::Relaxed);
                                Ok(Aind::Dummy(a as usize))
                            } else {
                                Err(AindError::NotIndex(format!("Negative index {}", a)))
                            }
                        } else {
                            Err(AindError::NotIndex(format!("Invalid index {}", arg)))
                        }
                    } else {
                        Err(AindError::ParsingError(format!(
                            "Too many arguments to dummy:{}",
                            f.as_view()
                        )))
                    }
                } else if f.get_symbol() == symbol!("edge") {
                    if f.get_nargs() <= 2 {
                        let mut iter = f.iter();
                        let i = iter.next().unwrap();
                        let j = iter.next().map(i64::try_from).unwrap_or(Ok(0));

                        if let (Ok(a), Ok(b)) = (i64::try_from(i), j) {
                            if a >= 0 && b >= 0 {
                                Ok(Aind::Edge(a as u16, b as u16))
                            } else {
                                Err(AindError::NotIndex(format!("Negative index {}", a)))
                            }
                        } else {
                            Err(AindError::NotIndex(format!(
                                "Invalid index {}",
                                f.as_view()
                            )))
                        }
                    } else {
                        Err(AindError::ParsingError(format!(
                            "Incorrect number of arguments to edge:{}",
                            f.as_view()
                        )))
                    }
                } else if f.get_symbol() == symbol!("hedge") {
                    if f.get_nargs() <= 2 {
                        let mut iter = f.iter();
                        let i = iter.next().unwrap();
                        let j = iter.next().map(i64::try_from).unwrap_or(Ok(0));

                        if let (Ok(a), Ok(b)) = (i64::try_from(i), j) {
                            if a >= 0 && b >= 0 {
                                Ok(Aind::Hedge(a as u16, b as u16))
                            } else {
                                Err(AindError::NotIndex(format!("Negative index {}", a)))
                            }
                        } else {
                            Err(AindError::NotIndex(format!(
                                "Invalid index {}",
                                f.as_view()
                            )))
                        }
                    } else {
                        Err(AindError::ParsingError(format!(
                            "Incorrect number of arguments to hedge:{}",
                            f.as_view()
                        )))
                    }
                } else if f.get_symbol() == symbol!("vertex") {
                    if f.get_nargs() <= 2 {
                        let mut iter = f.iter();
                        let i = iter.next().unwrap();
                        let j = iter.next().map(i64::try_from).unwrap_or(Ok(0));

                        if let (Ok(a), Ok(b)) = (i64::try_from(i), j) {
                            if a >= 0 && b >= 0 {
                                Ok(Aind::Vertex(a as u16, b as u16))
                            } else {
                                Err(AindError::NotIndex(format!("Negative index {}", a)))
                            }
                        } else {
                            Err(AindError::NotIndex(format!(
                                "Invalid index {}",
                                f.as_view()
                            )))
                        }
                    } else {
                        Err(AindError::ParsingError(format!(
                            "Incorrect number of arguments to vertex:{}",
                            f.as_view()
                        )))
                    }
                } else {
                    Err(AindError::NotIndex(format!(
                        "Invalid index {}",
                        f.as_view()
                    )))
                }
            }

            _ => Err(AindError::NotIndex(view.to_string())),
        }
    }
}

#[cfg(test)]
mod tests {
    use idenso::metric::MetricSimplifier;
    use spenso::{
        network::parsing::ShadowedStructure,
        structure::{permuted::Perm, PermutedStructure, TensorStructure},
        tensors::symbolic::SymbolicTensor,
    };
    use symbolica::{atom::AtomCore, parse_lit};

    use super::*;

    #[test]
    fn test_structure_parsing() {
        let expr = parse_lit!(gamma(
            spenso::mink(4, edge(1, 1)),
            spenso::mink(4, edge(1)),
            spenso::mink(4, hedge(1, 1)),
            spenso::mink(4, hedge(1)),
            spenso::mink(4, vertex(2, 1)),
            spenso::mink(4, vertex(2)),
            spenso::mink(4, dummy(111)),
            spenso::mink(4, 1)
        ));
        let structure = PermutedStructure::<ShadowedStructure<Aind>>::try_from(expr.clone());

        match structure {
            Ok(s) => {
                let pexpr = s
                    .clone()
                    .map_structure(|a| SymbolicTensor::from_named(&a).unwrap())
                    .permute_inds()
                    .expression
                    .simplify_metrics();
                assert_eq!(s.structure.order(), 8);
                assert_eq!(
                    expr,
                    pexpr,
                    "{}\n not equal to\n{}",
                    expr.to_canonical_string(),
                    pexpr.to_canonical_string()
                );
            }
            Err(e) => println!("Error parsing structure: {}", e),
        }
        // assert_eq!(structure.to_string(), "s₁₂");
    }
}
