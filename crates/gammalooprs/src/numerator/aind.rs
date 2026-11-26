use std::collections::BTreeSet;
use std::sync::atomic::{AtomicUsize, Ordering};

use linnet::half_edge::{
    NodeIndex,
    involution::{EdgeIndex, Hedge},
};
use spenso::structure::abstract_index::AbstractIndex;
use spenso::{
    structure::slot::{AbsInd, DummyAind, ParseableAind, SlotError},
    utils::{to_subscript, to_superscript},
};
use symbolica::id::Replacement;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol, representation::FunView},
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
        let index = DUMMYCOUNTER.fetch_add(1, Ordering::Relaxed);
        crate::debug_tags!(#generation, #inspect;
            aind_dummy = true,
            stage = "aind_new_dummy",
            dummy_index = index,
            "Aind dummy allocated"
        );
        Aind::Dummy(index)
    }

    fn is_dummy(&self) -> bool {
        matches!(self, Aind::Dummy(_))
    }

    fn new_dummy_at(i: usize) -> Self {
        Aind::Dummy(i)
    }
}

/// A contracted index name that one derivation allocates for itself.
///
/// Each allocator restarts its numbering for every expression it produces, so
/// independently derived expressions reuse these names. None of them ever
/// connects an expression to the rest of the graph.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum PrivateIndex {
    /// Integrated-UV tensor-reduction index, keyed only by topology order.
    UvTerm(u16, u16),
    /// GammaLoop dummy, numbered from one offset by every tensor-network parse.
    Dummy(usize),
    /// Spenso dummy `d_n`, numbered from zero by the color simplifier.
    Spenso(usize),
}

impl PrivateIndex {
    /// The private index that `view` names, if any. `spenso_namespace` is the
    /// namespace of Spenso's dummy variables.
    fn from_view(view: AtomView<'_>, spenso_namespace: &str) -> Option<Self> {
        match view {
            AtomView::Fun(function)
                if function.get_symbol() == GS.dummyaind || function.get_symbol() == GS.uvaind =>
            {
                // A malformed index is not a slot name that a product could
                // contract, so only valid indices are renamed.
                match Aind::try_from(view).ok()? {
                    Aind::Dummy(number) => Some(Self::Dummy(number)),
                    Aind::UVTerm(order, local) => Some(Self::UvTerm(order, local)),
                    _ => None,
                }
            }
            AtomView::Var(variable) => {
                let symbol = variable.get_symbol();
                if symbol.get_namespace() != spenso_namespace {
                    return None;
                }
                let number = symbol
                    .get_stripped_name()
                    .strip_prefix("d_")?
                    .parse()
                    .ok()?;
                let index = Self::Spenso(number);
                // A name with leading zeros is a different variable.
                (index.to_atom().as_view() == view).then_some(index)
            }
            _ => None,
        }
    }

    fn to_atom(self) -> Atom {
        match self {
            Self::UvTerm(order, local) => Aind::UVTerm(order, local).into(),
            Self::Dummy(number) => Aind::Dummy(number).into(),
            Self::Spenso(number) => AbstractIndex::Dummy(number).into(),
        }
    }
}

impl Aind {
    /// Renamings that keep the private indices of independently derived factors
    /// apart before the factors are multiplied.
    ///
    /// Each item of `factors` lists the expressions of one factor. Private
    /// indices are the dummies of tensor-network parses and of the color
    /// simplifier, and integrated-UV tensor-reduction indices: a derivation
    /// allocates them only for contractions inside its own expression. Every
    /// allocator restarts its numbering per expression, so separately derived
    /// factors reuse these names. Their unchanged product would contract
    /// unrelated slots of different factors with each other.
    ///
    /// The first factor to use a name keeps it. Later factors receive fresh
    /// names, renamed consistently across all of a factor's expressions: Spenso
    /// dummies beyond every Spenso dummy, and other private indices as dummies
    /// beyond every dummy of every factor. Open graph indices are never renamed.
    /// The result has one, possibly empty, replacement list per factor.
    pub(crate) fn disjoint_private_renamings<'a, F>(
        factors: impl IntoIterator<Item = F>,
    ) -> Vec<Vec<Replacement>>
    where
        F: IntoIterator<Item = AtomView<'a>>,
    {
        let spenso_namespace = Atom::from(AbstractIndex::Dummy(0))
            .get_symbol()
            .expect("Spenso names its dummies with variables")
            .get_namespace();
        let private_indices = factors
            .into_iter()
            .map(|atoms| {
                let mut indices = BTreeSet::new();
                for atom in atoms {
                    atom.visitor(&mut |view| {
                        let Some(index) = PrivateIndex::from_view(view, spenso_namespace) else {
                            return true;
                        };
                        indices.insert(index);
                        false
                    });
                }
                indices
            })
            .collect::<Vec<_>>();

        let (mut next_dummy, mut next_spenso) = (0, 0);
        for index in private_indices.iter().flatten() {
            match *index {
                PrivateIndex::Dummy(number) => next_dummy = next_dummy.max(number + 1),
                PrivateIndex::Spenso(number) => next_spenso = next_spenso.max(number + 1),
                PrivateIndex::UvTerm(..) => {}
            }
        }
        let mut used = BTreeSet::new();
        private_indices
            .into_iter()
            .map(|indices| {
                indices
                    .into_iter()
                    .filter(|index| !used.insert(*index))
                    .map(|index| {
                        // Keep color dummies in the form the color simplifier parses.
                        let fresh = match index {
                            PrivateIndex::Spenso(_) => {
                                next_spenso += 1;
                                PrivateIndex::Spenso(next_spenso - 1)
                            }
                            PrivateIndex::UvTerm(..) | PrivateIndex::Dummy(_) => {
                                next_dummy += 1;
                                PrivateIndex::Dummy(next_dummy - 1)
                            }
                        };
                        Replacement::new(index.to_atom().to_pattern(), fresh.to_atom())
                    })
                    .collect()
            })
            .collect()
    }

    /// `atom` with its private indices renamed apart from those of `others`.
    ///
    /// Use this before multiplying `atom` with expressions derived
    /// independently of it; see [`Self::disjoint_private_renamings`].
    pub(crate) fn private_indices_apart<'a>(
        atom: Atom,
        others: impl IntoIterator<Item = AtomView<'a>>,
    ) -> Atom {
        let mut renamings = Self::disjoint_private_renamings([
            others.into_iter().collect::<Vec<_>>(),
            vec![atom.as_view()],
        ]);
        let renaming = renamings.pop().expect("one renaming list per factor");
        if renaming.is_empty() {
            atom
        } else {
            atom.replace_multiple(&renaming)
        }
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
    crate::debug_tags!(#generation, #inspect;
        aind_dummy = true,
        stage = "aind_parse_dummy",
        dummy_index = i,
        "Aind dummy parsed"
    );
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
        shorthands::metric::MetricSimplifier,
        tensor::{SymbolicNetParse, SymbolicTensor},
    };
    use spenso::{
        network::parsing::{ParseSettings, ShadowedStructure, StructureFromAtom},
        structure::TensorStructure,
    };
    use symbolica::{
        atom::{Atom, AtomCore, FunctionBuilder},
        function, parse_lit, symbol,
    };

    use crate::initialisation::initialise;

    use super::*;

    #[test]
    fn disjoint_private_renamings_rename_only_reused_private_indices() {
        initialise().unwrap();

        let mink = |index: Atom| function!(symbol!("spenso::mink"), 4i64, index);
        let coad = |index: Atom| function!(symbol!("spenso::coad"), 8i64, index);
        let tensor = |name: &str, slots: &[Atom]| {
            let mut tensor = FunctionBuilder::new(symbol!(name));
            for slot in slots {
                tensor = tensor.add_arg(slot.as_view());
            }
            tensor.finish()
        };
        let dummy = |number| mink(Aind::Dummy(number).into());
        let color = |number| coad(AbstractIndex::Dummy(number).into());
        let open = mink(Aind::Hedge(3, 0).into());
        let first = tensor(
            "first",
            &[
                dummy(1_000_000),
                mink(Aind::UVTerm(1, 1).into()),
                color(0),
                open.clone(),
            ],
        );
        // The shared open index belongs to the product; the private names only
        // coincide because both factors were derived separately.
        let second = tensor(
            "second",
            &[
                dummy(1_000_000),
                mink(Aind::UVTerm(1, 1).into()),
                dummy(1_000_003),
                color(0),
                color(2),
                open.clone(),
            ],
        );
        let separated = tensor(
            "second",
            &[
                dummy(1_000_005),
                dummy(1_000_004),
                dummy(1_000_003),
                color(3),
                color(2),
                open,
            ],
        );

        let renamings = Aind::disjoint_private_renamings([[first.as_view()], [second.as_view()]]);

        assert!(renamings[0].is_empty());
        assert_eq!(second.replace_multiple(&renamings[1]), separated);
        assert_eq!(
            Aind::private_indices_apart(second, [first.as_view()]),
            separated
        );
    }

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
                let pexpr = SymbolicTensor::from_canonicalized(&s)
                    .unwrap()
                    .expression
                    .simplify_metrics();
                assert_eq!(s.canonical().order(), 8);
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
