use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::atomic::{AtomicUsize, Ordering};

use linnet::half_edge::NodeIndex;
use linnet::half_edge::involution::{EdgeIndex, Hedge};
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::slot::{AbsInd, DummyAind, ParseableAind, SlotError};
use spenso::utils::{to_subscript, to_superscript};
use symbolica::atom::representation::FunView;
use symbolica::atom::{Atom, AtomCore, AtomOrView, AtomView, Symbol};
use symbolica::coefficient::CoefficientView;
use symbolica::{function, symbol};
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
    Symbol(Symbol),
    Dummy(usize),
    /// A private dummy owned by the producer that `Symbol` names.
    ///
    /// See [`Aind::scope_private_indices`].
    ScopedDummy(Symbol, usize),
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

fn allocate_dummy_from(counter: &AtomicUsize) -> Aind {
    let index = counter
        .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |current| {
            i64::try_from(current).ok()?;
            current.checked_add(1)
        })
        .expect("dummy index counter exhausted representable range");
    crate::debug_tags!(#generation, #inspect;
        aind_dummy = true,
        stage = "aind_new_dummy",
        dummy_index = index,
        "Aind dummy allocated"
    );
    Aind::Dummy(index)
}

fn record_parsed_dummy(counter: &AtomicUsize, index: usize) {
    counter.fetch_max(index.saturating_add(1), Ordering::Relaxed);
}

impl DummyAind for Aind {
    fn new_dummy() -> Self {
        allocate_dummy_from(&DUMMYCOUNTER)
    }

    fn is_dummy(&self) -> bool {
        matches!(self, Aind::Dummy(_) | Aind::ScopedDummy(..))
    }

    fn new_dummy_at(i: usize) -> Self {
        Aind::Dummy(i)
    }
}

impl Aind {
    /// `atom` with every unscoped private index owned by `scope`.
    ///
    /// Private indices are contractions that one derivation allocates inside
    /// its own expression: GammaLoop dummies, renumbered by every
    /// tensor-network parse, and Spenso dummies `d_n` of the color simplifier.
    /// Each allocator restarts its numbering, so independently derived
    /// expressions reuse these names, and multiplying them would contract
    /// unrelated slots. A producer whose output is multiplied with other
    /// graph pieces stamps it once with a scope symbol naming that producer:
    /// the indices become [`Aind::ScopedDummy`] values numbered in first
    /// occurrence order, skipping the numbers that `atom` already uses in
    /// `scope`. Indices already scoped, and open graph indices, are kept;
    /// `atom` is moved back unchanged and `scope` is not evaluated when `atom`
    /// has no unscoped private index.
    pub(crate) fn scope_private_indices(atom: Atom, scope: impl FnOnce() -> Symbol) -> Atom {
        match Self::scope_private_indices_cow(&atom, scope) {
            Cow::Owned(scoped) => scoped,
            Cow::Borrowed(_) => atom,
        }
    }

    /// Borrowing form of [`Aind::scope_private_indices`], which borrows
    /// `atom` without allocating or evaluating `scope` when `atom` has no
    /// unscoped private index.
    #[must_use]
    pub(crate) fn scope_private_indices_cow(
        atom: &Atom,
        scope: impl FnOnce() -> Symbol,
    ) -> Cow<'_, Atom> {
        let mut first_occurrences: Vec<Atom> = Vec::new();
        atom.visitor(&mut |view| {
            if is_unscoped_private_index(view) {
                if first_occurrences
                    .iter()
                    .all(|index| index.as_view() != view)
                {
                    first_occurrences.push(view.to_owned());
                }
                return false;
            }
            !(matches!(view, AtomView::Fun(function) if function.get_symbol() == GS.dummyaind))
        });
        if first_occurrences.is_empty() {
            return Cow::Borrowed(atom);
        }
        let scope = scope();
        let mut occupied = BTreeSet::new();
        atom.visitor(&mut |view| {
            let AtomView::Fun(function) = view else {
                return true;
            };
            if function.get_symbol() != GS.dummyaind {
                return true;
            }
            if let Ok(Aind::ScopedDummy(owner, number)) = peek_dummy_aind(function)
                && owner == scope
            {
                occupied.insert(number);
            }
            false
        });
        let free = (0..).filter(|number| !occupied.contains(number));
        let renaming: BTreeMap<Atom, Atom> = first_occurrences
            .into_iter()
            .zip(free)
            .map(|(index, number)| (index, Aind::ScopedDummy(scope, number).to_atom()))
            .collect();
        Cow::Owned(atom.replace_map(|view, _, out| {
            if is_unscoped_private_index(view)
                && let Some(scoped) = renaming.get(&view.to_owned())
            {
                **out = scoped.clone();
            }
        }))
    }

    /// Whether no private index, scoped or not, occurs in two of `factors`.
    ///
    /// Each factor lists its expressions. Products of independently derived
    /// factors rely on this; it holds when each producer scoped its output.
    pub(crate) fn private_indices_are_disjoint<'a, F>(factors: impl IntoIterator<Item = F>) -> bool
    where
        F: IntoIterator<Item = AtomView<'a>>,
    {
        let mut seen = BTreeSet::new();
        factors.into_iter().all(|factor| {
            let mut indices = BTreeSet::new();
            for atom in factor {
                atom.visitor(&mut |view| {
                    let private = match view {
                        AtomView::Fun(function) if function.get_symbol() == GS.dummyaind => {
                            peek_dummy_aind(function).is_ok_and(|index| index.is_dummy())
                        }
                        _ => is_unscoped_private_index(view),
                    };
                    if private {
                        indices.insert(view.to_owned());
                    }
                    !private
                });
            }
            indices.into_iter().all(|index| seen.insert(index))
        })
    }

    /// One past the largest unscoped dummy number in `atom`, or zero if none
    /// occurs.
    ///
    /// Every number from the result upwards is free, so callers may offset it
    /// by any natural number. Inspection does not advance the process-wide
    /// dummy counter.
    pub(crate) fn next_free_dummy(atom: &Atom) -> usize {
        let mut next = 0;
        atom.visitor(&mut |view| {
            if let AtomView::Fun(function) = view
                && function.get_symbol() == GS.dummyaind
                && let Ok(Aind::Dummy(number)) = peek_dummy_aind(function)
            {
                next = next.max(number + 1);
            }
            true
        });
        next
    }
}

/// Whether `view` names a private index that no producer has scoped yet,
/// without advancing the process-wide dummy counter.
fn is_unscoped_private_index(view: AtomView<'_>) -> bool {
    match view {
        AtomView::Fun(function) if function.get_symbol() == GS.dummyaind => {
            // A malformed index is not a slot name that a product could
            // contract, so only valid dummies are scoped.
            matches!(peek_dummy_aind(function), Ok(Aind::Dummy(_)))
        }
        AtomView::Var(variable) => variable
            .get_symbol()
            .get_stripped_name()
            .strip_prefix("d_")
            .and_then(|number| number.parse().ok())
            // A name with leading zeros, or outside Spenso, is another variable.
            .is_some_and(|number| Atom::from(AbstractIndex::Dummy(number)).as_view() == view),
        _ => false,
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
            Aind::ScopedDummy(scope, i) => {
                if f.sign_minus() {
                    write!(f, "d{scope}{}", to_subscript(*i as isize))
                } else if f.sign_plus() {
                    write!(f, "d{scope}{}", to_superscript(*i as isize))
                } else {
                    write!(f, "d{scope}-{i}")
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
            Aind::ScopedDummy(scope, i) => function!(GS.dummyaind, Atom::var(scope), i as i64),
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
impl<'a> From<Aind> for AtomOrView<'a> {
    fn from(value: Aind) -> Self {
        AtomOrView::Atom(Atom::from(value))
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

/// Parse a `dummy(...)` index view without updating [`DUMMYCOUNTER`].
fn peek_dummy_aind(f: FunView<'_>) -> Result<Aind, AindError> {
    let mut args = f.iter();
    let (scope, index) = match (args.next(), args.next(), args.next()) {
        (Some(index), None, _) => (None, index),
        (Some(AtomView::Var(scope)), Some(index), None) => (Some(scope.get_symbol()), index),
        _ => {
            return Err(AindError::ParsingError(format!(
                "Incorrect number of arguments to dummy:{}",
                f.as_view()
            )));
        }
    };
    let i = usize::try_from(parse_natural_i64(index)?)
        .map_err(|e| AindError::NotIndex(e.to_string()))?;
    Ok(match scope {
        Some(scope) => Aind::ScopedDummy(scope, i),
        None => Aind::Dummy(i),
    })
}

fn parse_dummy_aind(f: FunView<'_>) -> Result<Aind, AindError> {
    let parsed = peek_dummy_aind(f)?;
    if let Aind::Dummy(i) = parsed {
        record_parsed_dummy(&DUMMYCOUNTER, i);
        crate::debug_tags!(#generation, #inspect;
            aind_dummy = true,
            stage = "aind_parse_dummy",
            dummy_index = i,
            "Aind dummy parsed"
        );
    }
    Ok(parsed)
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
    use idenso::dirac::AGS;
    use idenso::shorthands::metric::MetricSimplifier;
    use idenso::tensor::{SymbolicNetParse, SymbolicTensor};
    use spenso::network::parsing::{ParseSettings, ShadowedStructure, StructureFromAtom};
    use spenso::structure::TensorStructure;
    use symbolica::atom::{Atom, AtomCore};
    use symbolica::{function, parse_lit, symbol};

    use idenso::IndexTooling;
    use idenso::representations::{Bispinor, ColorAdjoint};
    use spenso::structure::representation::{Minkowski, RepName};

    use crate::initialisation::initialise;

    use super::*;

    #[test]
    fn scope_private_indices_stamps_only_unscoped_private_indices() {
        initialise().unwrap();

        let scope = symbol!("aind_test::producer");
        let other = symbol!("aind_test::other_producer");
        let mink = |index: Aind| Minkowski {}.new_rep(4).to_symbolic([index.to_atom()]);
        let coad = |index: AbstractIndex| ColorAdjoint {}.new_rep(8).to_symbolic([index.into()]);
        let open = mink(Aind::Hedge(3, 0));
        let already_scoped = mink(Aind::ScopedDummy(other, 0));
        let expression = |momentum: Atom, color: Atom| {
            function!(
                symbol!("aind_test::first"),
                &momentum,
                &color,
                &already_scoped,
                &open
            ) * function!(symbol!("aind_test::second"), momentum, color)
        };
        let input = expression(mink(Aind::Dummy(1_000_000)), coad(AbstractIndex::Dummy(0)));

        let stamped = Aind::scope_private_indices(input, || scope);

        let scoped = [0, 1].map(|index| mink(Aind::ScopedDummy(scope, index)));
        let color_scoped = [0, 1].map(|index| {
            ColorAdjoint {}
                .new_rep(8)
                .to_symbolic([Aind::ScopedDummy(scope, index).to_atom()])
        });
        // First-occurrence order fixes which private index is numbered first.
        assert!(
            [
                expression(scoped[0].clone(), color_scoped[1].clone()),
                expression(scoped[1].clone(), color_scoped[0].clone()),
            ]
            .contains(&stamped),
            "{stamped}"
        );
        assert!(matches!(
            Aind::scope_private_indices_cow(&stamped, || {
                panic!("already-scoped expressions must not evaluate the scope builder")
            }),
            Cow::Borrowed(_)
        ));
        assert!(matches!(
            Aind::scope_private_indices_cow(&open, || {
                panic!("open indices must not evaluate the scope builder")
            }),
            Cow::Borrowed(_)
        ));
    }

    #[test]
    fn scope_private_indices_skips_numbers_the_scope_already_owns() {
        initialise().unwrap();

        let scope = symbol!("aind_test::occupied_producer");
        let mink = |index: Aind| Minkowski {}.new_rep(4).to_symbolic([index.to_atom()]);
        let coad = |index: Atom| ColorAdjoint {}.new_rep(8).to_symbolic([index]);
        let contraction = |name: &str, slot: Atom| {
            function!(symbol!(&format!("aind_test::{name}_left")), &slot)
                * function!(symbol!(&format!("aind_test::{name}_right")), slot)
        };
        let input = contraction("owned", mink(Aind::ScopedDummy(scope, 0)))
            * contraction("gammaloop", mink(Aind::Dummy(9)))
            * contraction("spenso", coad(AbstractIndex::Dummy(0).into()));

        let stamped = Aind::scope_private_indices(input, || scope);

        let reference = |gammaloop: usize, spenso: usize| {
            contraction("owned", mink(Aind::ScopedDummy(scope, 0)))
                * contraction("gammaloop", mink(Aind::ScopedDummy(scope, gammaloop)))
                * contraction("spenso", coad(Aind::ScopedDummy(scope, spenso).to_atom()))
        };
        // The two unscoped indices take the free numbers in first occurrence
        // order, never the number the scope already owns.
        assert!(
            [reference(1, 2), reference(2, 1)].contains(&stamped),
            "{stamped}"
        );
    }

    #[test]
    #[should_panic(expected = "dummy index counter exhausted representable range")]
    fn parsed_maximum_dummy_exhausts_counter_instead_of_reusing_it() {
        let counter = AtomicUsize::new(0);
        let penultimate = (i64::MAX - 1) as usize;
        let maximum = i64::MAX as usize;
        record_parsed_dummy(&counter, penultimate);
        record_parsed_dummy(&counter, maximum);
        let _ = allocate_dummy_from(&counter);
    }

    #[test]
    fn scoped_dummies_round_trip_and_canonize_as_dummies() {
        initialise().unwrap();

        let scoped = Aind::ScopedDummy(symbol!("aind_test::round_trip"), 3);
        assert!(scoped.is_dummy());
        assert_eq!(Aind::try_from(scoped.to_atom().as_view()).unwrap(), scoped);
        assert_eq!(
            scoped,
            Aind::ScopedDummy(symbol!("aind_test::round_trip"), 3)
        );
        assert_ne!(scoped, Aind::ScopedDummy(symbol!("aind_test::other"), 3));

        let slashed = |index: Aind| {
            let mu = Minkowski {}.new_rep(4).to_symbolic([index.to_atom()]);
            function!(GS.emr_mom, 1, &mu)
                * function!(
                    AGS.gamma,
                    Bispinor {}
                        .new_rep(4)
                        .to_symbolic([Aind::Hedge(0, 0).to_atom()]),
                    Bispinor {}
                        .new_rep(4)
                        .to_symbolic([Aind::Hedge(1, 0).to_atom()]),
                    mu
                )
        };
        assert_eq!(
            slashed(scoped).canonize(Aind::Dummy).unwrap(),
            slashed(Aind::Dummy(5)).canonize(Aind::Dummy).unwrap()
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
            (
                Aind::ScopedDummy(symbol!("aind_test::scope"), 8),
                function!(GS.dummyaind, symbol!("aind_test::scope"), 8i64),
            ),
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
