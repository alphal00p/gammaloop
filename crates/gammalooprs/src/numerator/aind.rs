use std::sync::atomic::{AtomicUsize, Ordering};

use linnet::half_edge::{
    NodeIndex,
    involution::{EdgeIndex, Hedge},
};
use spenso::{
    structure::slot::{AbsInd, DummyAind, ParseableAind, SlotError},
    utils::{to_subscript, to_superscript},
};
use symbolica::{
    atom::{Atom, AtomView, Symbol, representation::FunView},
    coefficient::CoefficientView,
    function, symbol,
};
use thiserror::Error;

use crate::utils::GS;

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
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    // bincode_trait_derive::BorrowDecodeFromDecode,
)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub enum Aind {
    Normal(usize),
    Hedge(u16, u16),
    Edge(u16, u16),
    Vertex(u16, u16),
    UVTerm(u16, u16),
    Symbol(Symbol),
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

    fn is_dummy(&self) -> bool {
        matches!(self, Aind::Dummy(_))
    }

    fn new_dummy_at(i: usize) -> Self {
        Aind::Dummy(i)
    }
}

impl ParseableAind for Aind {
    type Error = AindError;

    fn from_view(view: AtomView<'_>) -> Result<Self, Self::Error> {
        view.try_into()
    }

    fn to_atom(&self) -> Atom {
        (*self).into()
    }
}

impl std::fmt::Display for Aind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Aind::Symbol(v) => {
                if f.sign_minus() {
                    write!(f, "_{}", v)
                } else if f.sign_plus() {
                    write!(f, "^{}", v)
                } else {
                    write!(f, "{}", v)
                }
            }
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
            Aind::UVTerm(i, j) => {
                if f.sign_minus() {
                    write!(
                        f,
                        "u{}.{}",
                        to_subscript(*i as isize),
                        to_subscript(*j as isize)
                    )
                } else if f.sign_plus() {
                    write!(
                        f,
                        "u{}'{}",
                        to_superscript(*i as isize),
                        to_superscript(*j as isize)
                    )
                } else {
                    write!(f, "u{}-{}", i, j)
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
            Aind::Symbol(s) => Atom::var(s),
            Aind::Dummy(i) => function!(GS.dummyaind, i as i64),
            Aind::Normal(i) => Atom::num(i as i64),
            Aind::UVTerm(i, j) => {
                if j != 0 {
                    function!(GS.uvaind, i as i64, j as i64)
                } else {
                    function!(GS.uvaind, i as i64)
                }
            }
            Aind::Edge(i, j) => {
                if j != 0 {
                    function!(GS.edgeaind, i as i64, j as i64)
                } else {
                    function!(GS.edgeaind, i as i64)
                }
            }
            Aind::Hedge(i, j) => {
                if j != 0 {
                    function!(GS.hedgeaind, i as i64, j as i64)
                } else {
                    function!(GS.hedgeaind, i as i64)
                }
            }
            Aind::Vertex(i, j) => {
                if j != 0 {
                    function!(GS.vertexaind, i as i64, j as i64)
                } else {
                    function!(GS.vertexaind, i as i64)
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

fn parse_natural_i64(arg: AtomView<'_>) -> Result<i64, AindError> {
    let index =
        i64::try_from(arg).map_err(|_| AindError::NotIndex(format!("Invalid index {arg}")))?;
    if index >= 0 {
        Ok(index)
    } else {
        Err(AindError::NotIndex(format!("Negative index {index}")))
    }
}

fn parse_dummy_aind(f: FunView<'_>) -> Result<Aind, AindError> {
    if f.get_nargs() != 1 {
        return Err(AindError::ParsingError(format!(
            "Incorrect number of arguments to dummy:{}",
            f.as_view()
        )));
    }

    let i = parse_natural_i64(f.iter().next().unwrap())?;
    let i = usize::try_from(i).map_err(|e| AindError::NotIndex(e.to_string()))?;
    DUMMYCOUNTER.fetch_max(i + 1, Ordering::Relaxed);
    Ok(Aind::Dummy(i))
}

fn parse_u16_aind_pair(f: FunView<'_>, name: &str) -> Result<(u16, u16), AindError> {
    if !(1..=2).contains(&f.get_nargs()) {
        return Err(AindError::ParsingError(format!(
            "Incorrect number of arguments to {name}:{}",
            f.as_view()
        )));
    }

    let mut iter = f.iter();
    let i = parse_natural_i64(iter.next().unwrap())?;
    let j = iter.next().map(parse_natural_i64).unwrap_or(Ok(0))?;

    let i = u16::try_from(i).map_err(|e| AindError::NotIndex(e.to_string()))?;
    let j = u16::try_from(j).map_err(|e| AindError::NotIndex(e.to_string()))?;

    Ok((i, j))
}

impl TryFrom<AtomView<'_>> for Aind {
    type Error = AindError;

    fn try_from(view: AtomView<'_>) -> Result<Self, Self::Error> {
        match view {
            AtomView::Var(v) => Ok(Aind::Symbol(v.get_symbol())),
            AtomView::Num(n) => match n.get_coeff_view() {
                CoefficientView::Natural(n, 1, _, _) => Ok(Aind::Normal(
                    usize::try_from(n).map_err(|e| AindError::NotIndex(e.to_string()))?,
                )),
                _ => Err(AindError::NotNatural),
            },
            AtomView::Fun(f) => {
                let symbol = f.get_symbol();
                if symbol == GS.dummyaind || symbol == symbol!("dummy") {
                    parse_dummy_aind(f)
                } else if symbol == GS.edgeaind || symbol == symbol!("edge") {
                    let (i, j) = parse_u16_aind_pair(f, "edge")?;
                    Ok(Aind::Edge(i, j))
                } else if symbol == GS.hedgeaind || symbol == symbol!("hedge") {
                    let (i, j) = parse_u16_aind_pair(f, "hedge")?;
                    Ok(Aind::Hedge(i, j))
                } else if symbol == GS.vertexaind {
                    let (i, j) = parse_u16_aind_pair(f, "vertex")?;
                    Ok(Aind::Vertex(i, j))
                } else if symbol == GS.uvaind {
                    let (i, j) = parse_u16_aind_pair(f, "uv")?;
                    Ok(Aind::UVTerm(i, j))
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
    use idenso::{
        dirac::AGS,
        metric::MetricSimplifier,
        tensor::{SymbolicNetParse, SymbolicTensor},
    };
    use spenso::{
        network::parsing::{ParseSettings, ShadowedStructure, StructureFromAtom},
        structure::{TensorStructure, permuted::Perm},
    };
    use symbolica::{
        atom::{Atom, AtomCore},
        function, parse_lit, symbol,
    };

    use crate::initialisation::initialise;

    use super::*;

    #[test]
    fn parses_gammaloop_owned_aind_symbols() {
        initialise().unwrap();

        let cases = [
            (Aind::Dummy(111), function!(GS.dummyaind, 111i64)),
            (Aind::Edge(1, 0), function!(GS.edgeaind, 1i64)),
            (Aind::Edge(1, 2), function!(GS.edgeaind, 1i64, 2i64)),
            (Aind::Hedge(3, 0), function!(GS.hedgeaind, 3i64)),
            (Aind::Hedge(3, 4), function!(GS.hedgeaind, 3i64, 4i64)),
            (Aind::Vertex(5, 0), function!(GS.vertexaind, 5i64)),
            (Aind::Vertex(5, 6), function!(GS.vertexaind, 5i64, 6i64)),
            (Aind::UVTerm(7, 0), function!(GS.uvaind, 7i64)),
            (Aind::UVTerm(7, 8), function!(GS.uvaind, 7i64, 8i64)),
        ];

        for (expected, atom) in cases {
            assert_eq!(Aind::try_from(atom.as_view()).unwrap(), expected);
        }
    }

    #[test]
    fn parses_indexed_gammaloop_tensors_without_compact_rewrite() {
        initialise().unwrap();

        let bis0 = function!(symbol!("spenso::bis"), 4i64, Atom::from(Aind::Hedge(0, 0)));
        let bis1 = function!(symbol!("spenso::bis"), 4i64, Atom::from(Aind::Hedge(1, 0)));
        let mink = function!(symbol!("spenso::mink"), 4i64, Atom::from(Aind::Edge(2, 1)));

        let expr = function!(GS.ubar, 1i64, bis0.clone())
            * function!(GS.u, 1i64, bis1.clone())
            * function!(GS.emr_mom, 2i64, mink.clone())
            * function!(AGS.gamma, bis0, bis1, mink);

        let net = expr
            .parse_to_symbolic_net::<Aind>(&ParseSettings::default())
            .unwrap();
        let dangling = net.graph.dangling_indices();

        assert!(
            dangling.is_empty(),
            "indexed GammaLoop tensors were parsed with dangling slots: {dangling:?}"
        );
    }

    #[test]
    fn test_structure_parsing() {
        initialise().unwrap();
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
        let structure = ShadowedStructure::<Aind>::parse(expr.as_view());

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
