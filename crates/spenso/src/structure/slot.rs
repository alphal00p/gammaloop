use super::{
    abstract_index::{AbstractIndex, AbstractIndexError},
    dimension::DimensionError,
    representation::{
        BaseRepName, LibraryRep, LibrarySlot, RepName, Representation, RepresentationError,
    },
};
use crate::structure::dimension::Dimension;
#[cfg(feature = "shadowing")]
use crate::{network::tags::SPENSO_TAG, structure::abstract_index::AIND_SYMBOLS};
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
    atom::{
        Atom, AtomView, Symbol,
        representation::{FunView, ListIterator},
    },
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
    #[error("Not a slot, isn't a representation")]
    NotRepresentation,
    #[error("Not a slot, is composite")]
    Composite,
    #[error(
        "Expected rep(dimension, index), optionally inside one single-argument duality wrapper"
    )]
    Malformed,
    #[error("{0}")]
    DimErr(#[from] DimensionError),
    #[error("{0}")]
    Any(#[from] eyre::Error),
}

#[cfg(feature = "shadowing")]
#[derive(Clone, Copy)]
enum SlotHead {
    Representation,
    Wrapper,
    Other,
}

#[cfg(feature = "shadowing")]
impl SlotHead {
    fn from_symbol(head: Symbol) -> Self {
        if head.has_tag(&SPENSO_TAG.representation) {
            Self::Representation
        } else if [
            AIND_SYMBOLS.dind,
            AIND_SYMBOLS.uind,
            AIND_SYMBOLS.selfdualind,
        ]
        .contains(&head)
        {
            Self::Wrapper
        } else {
            Self::Other
        }
    }
}

#[cfg(feature = "shadowing")]
/// Borrowed slot syntax. Recognition leaves dimension and index payloads unvalidated.
#[derive(Clone, Copy, Debug)]
pub struct SlotView<'a> {
    function: FunView<'a>,
    wrapper: Option<Symbol>,
    arguments: ListIterator<'a>,
}

#[cfg(feature = "shadowing")]
/// Classification of an expression for a structural tensor-index walk.
pub enum SlotMatch<'a> {
    Explicit(SlotView<'a>),
    /// A compact or malformed slot: its payload is opaque to tensor-index scans.
    Opaque,
    Other,
}

#[cfg(feature = "shadowing")]
impl<'a> SlotMatch<'a> {
    fn into_slot(self, value: AtomView<'a>) -> Result<SlotView<'a>, SlotError> {
        match self {
            Self::Explicit(slot) => Ok(slot),
            Self::Opaque => Err(SlotError::Malformed),
            Self::Other if matches!(value, AtomView::Fun(_)) => Err(SlotError::NotRepresentation),
            Self::Other => Err(SlotError::Composite),
        }
    }
}

#[cfg(feature = "shadowing")]
impl<'a> SlotView<'a> {
    /// Decode arguments once after recognizing a representation or duality wrapper.
    #[inline]
    fn classify(
        value: AtomView<'a>,
        mut classify_head: impl FnMut(FunView<'a>) -> SlotHead,
    ) -> SlotMatch<'a> {
        let AtomView::Fun(mut function) = value else {
            return SlotMatch::Other;
        };
        let mut wrapper = None;
        match classify_head(function) {
            SlotHead::Other => return SlotMatch::Other,
            SlotHead::Representation => {}
            SlotHead::Wrapper => {
                let mut arguments = function.iter();
                if arguments.len() != 1 {
                    return SlotMatch::Opaque;
                }
                let Some(AtomView::Fun(inner)) = arguments.next() else {
                    return SlotMatch::Opaque;
                };
                if !matches!(classify_head(inner), SlotHead::Representation) {
                    return SlotMatch::Opaque;
                }
                wrapper = Some(function.get_symbol());
                function = inner;
            }
        }
        let arguments = function.iter();
        if arguments.len() != 2 {
            return SlotMatch::Opaque;
        }
        SlotMatch::Explicit(Self {
            function,
            wrapper,
            arguments,
        })
    }

    /// Borrow the dimension expression without requiring a number or symbol.
    #[inline]
    pub fn dimension(mut self) -> AtomView<'a> {
        self.arguments.next().unwrap()
    }

    #[inline]
    /// Borrow the original explicit index payload without coercing its value.
    pub fn index(mut self) -> AtomView<'a> {
        self.arguments.next();
        self.arguments.next().unwrap()
    }

    fn parse<T: RepName, Aind: ParseableAind>(
        mut self,
        rep: T,
    ) -> Result<Slot<T, Aind>, SlotError> {
        let dim = Dimension::try_from(self.arguments.next().unwrap())?;
        let aind = Aind::from_view(self.arguments.next().unwrap()).map_err(Into::into)?;
        Ok(Slot {
            rep: Representation { dim, rep },
            aind,
        })
    }
}

#[cfg(feature = "shadowing")]
/// Recognition state shared by slot scans and typed structure inference.
///
/// The small, bounded caches avoid repeated tag lookups and representation
/// resolution without making scans quadratic in the number of distinct heads.
#[derive(Default)]
pub struct SlotMatcher {
    heads: [Option<(u32, SlotHead)>; 16],
    resolved: Vec<(Symbol, Option<Symbol>, LibraryRep)>,
}

#[cfg(feature = "shadowing")]
impl SlotMatcher {
    #[inline]
    fn cache_index(id: u32) -> usize {
        // Symbol IDs can share their low bits; mix before selecting one of 16 buckets.
        (id.wrapping_mul(0x9e37_79b9) >> 28) as usize
    }

    #[inline]
    fn classify_head(&mut self, function: FunView<'_>) -> SlotHead {
        let id = function.get_symbol_id();
        let entry = &mut self.heads[Self::cache_index(id)];
        if let Some((cached_id, kind)) = *entry
            && cached_id == id
        {
            return kind;
        }
        let kind = SlotHead::from_symbol(function.get_symbol());
        *entry = Some((id, kind));
        kind
    }

    /// Identify a slot or an opaque slot payload in one classification pass.
    #[inline]
    pub fn classify<'a>(&mut self, value: AtomView<'a>) -> SlotMatch<'a> {
        SlotView::classify(value, |function| self.classify_head(function))
    }

    /// Validate a recognized slot using the requested representation and index types.
    /// Use [`SlotView::index`] when exact symbolic index identity must be preserved.
    pub fn parse<T: RepName, Aind: ParseableAind>(
        &mut self,
        value: AtomView<'_>,
    ) -> Result<Slot<T, Aind>, SlotError> {
        let slot = self.classify(value).into_slot(value)?;
        let rep = self.representation(slot)?;
        slot.parse(T::from_library_rep(rep)?)
    }

    /// Resolve the representation and duality while retaining arbitrary dimension
    /// and index expressions in the borrowed view.
    pub fn representation(
        &mut self,
        slot: SlotView<'_>,
    ) -> Result<LibraryRep, RepresentationError> {
        let head = slot.function.get_symbol();
        let rep = if let Some((_, _, rep)) = self
            .resolved
            .iter()
            .find(|(symbol, wrapper, _)| *symbol == head && *wrapper == slot.wrapper)
        {
            *rep
        } else {
            let rep = match slot.wrapper {
                Some(wrapper) => LibraryRep::try_from_symbol(head, wrapper)?,
                None => LibraryRep::try_from_symbol_coerced(head)?,
            };
            let entry = (head, slot.wrapper, rep);
            if self.resolved.len() < 16 {
                self.resolved.push(entry);
            } else {
                self.resolved[Self::cache_index(head.get_id())] = entry;
            }
            rep
        };
        Ok(rep)
    }
}

#[cfg(feature = "shadowing")]
/// Parse `rep(d, i)` or one single-argument duality wrapper in the normalized Atom tree.
impl<'a, T: RepName, Aind: ParseableAind> TryFrom<AtomView<'a>> for Slot<T, Aind> {
    type Error = SlotError;

    fn try_from(value: AtomView<'a>) -> Result<Self, Self::Error> {
        let slot = SlotView::classify(value, |function| {
            SlotHead::from_symbol(function.get_symbol())
        })
        .into_slot(value)?;
        let head = slot.function.get_symbol();
        let rep = match slot.wrapper {
            Some(wrapper) => T::try_from_symbol(head, wrapper)?,
            None => T::try_from_symbol_coerced(head)?,
        };
        slot.parse(rep)
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

/// One external tensor axis: representation, dimension, and abstract index.
///
/// Slot equality includes all three fields. Dummy abstract indices remain
/// ordinary identity-bearing values; their dummy classification is not a
/// wildcard for contraction.
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

/// Duality and contraction matching for tensor slots.
pub trait DualSlotTo: IsAbstractSlot {
    /// Slot type obtained by reversing the representation orientation.
    type Dual: IsAbstractSlot;
    /// Returns the same dimension and abstract index in the dual
    /// representation.
    fn dual(&self) -> Self::Dual;
    /// Reports whether `other` has the same dimension and abstract index and a
    /// representation dual-compatible with `self`.
    fn matches(&self, other: &Self::Dual) -> bool;

    /// Ordering counterpart of [`Self::matches`] used by canonical structure
    /// merge algorithms.
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
    use symbolica::{
        atom::{Atom, AtomCore, NamespacedSymbol, SymbolBuilder},
        function, parse, symbol,
    };

    use crate::{
        network::tags::SPENSO_TAG,
        structure::{
            abstract_index::{AIND_SYMBOLS, AbstractIndex},
            representation::{
                DualLorentz, IndexDisplay, IndexPalette, IndexRow, LibraryRep, Lorentz, Minkowski,
                RepName, Representation, RepresentationClass, RepresentationError,
                RepresentationMetadata,
            },
            slot::{DualSlotTo, IsAbstractSlot, Slot, SlotError, SlotMatch, SlotMatcher},
        },
    };

    #[test]
    fn borrowed_slot_recognition_reads_only_explicit_index_positions() {
        let rep = LibraryRep::from(Lorentz {}).symbol();
        let index = Atom::var(symbol!("slot_match_tests::mu"));
        let slot = function!(rep, 4, &index);
        let dimension_function = function!(symbol!("slot_match_tests::dimension"), 4);
        let compound_index = function!(symbol!("slot_match_tests::label"), 7);
        let mut matcher = SlotMatcher::default();

        for atom in [
            slot.clone(),
            function!(AIND_SYMBOLS.dind, &slot),
            function!(AIND_SYMBOLS.uind, &slot),
            function!(AIND_SYMBOLS.selfdualind, &slot),
            function!(rep, &dimension_function, &index),
        ] {
            assert_eq!(
                matcher
                    .classify(atom.as_view())
                    .into_slot(atom.as_view())
                    .unwrap()
                    .index(),
                index.as_view(),
                "{atom}"
            );
        }
        let atom = function!(rep, 4, &compound_index);
        assert_eq!(
            matcher
                .classify(atom.as_view())
                .into_slot(atom.as_view())
                .unwrap()
                .index(),
            compound_index.as_view()
        );

        for atom in [
            Atom::num(4),
            index.clone(),
            function!(rep),
            function!(rep, 4),
            function!(rep, 4, &index, 0),
            function!(AIND_SYMBOLS.dind, function!(rep, 4)),
            function!(AIND_SYMBOLS.dind, &slot, 0),
            function!(symbol!("slot_match_tests::foreign"), 4, &index),
        ] {
            assert!(
                !matches!(matcher.classify(atom.as_view()), SlotMatch::Explicit(_)),
                "{atom}"
            );
        }
    }

    #[test]
    fn slot_head_cache_spreads_aligned_trace_heads() {
        // Observed mink, epsilon, and metric IDs in the axial trace benchmark.
        let buckets = [32, 128, 144].map(SlotMatcher::cache_index);
        assert_ne!(buckets[0], buckets[1]);
        assert_ne!(buckets[0], buckets[2]);
        assert_ne!(buckets[1], buckets[2]);
    }

    #[test]
    fn slot_head_cache_checks_keys_after_collision() {
        let rep = LibraryRep::from(Lorentz {}).symbol();
        let colliding = (0..64)
            .map(|index| {
                let name = format!("slot_match_tests::cache_collision_{index}");
                SymbolBuilder::new(NamespacedSymbol::parse(&name))
                    .build()
                    .unwrap()
            })
            .find(|head| {
                SlotMatcher::cache_index(head.get_id()) == SlotMatcher::cache_index(rep.get_id())
            })
            .unwrap();
        let slot = function!(rep, 4, 1);
        let ordinary = function!(colliding, 4, 1);
        let mut matcher = SlotMatcher::default();
        for _ in 0..3 {
            assert!(matches!(
                matcher.classify(slot.as_view()),
                SlotMatch::Explicit(_)
            ));
            assert!(matches!(
                matcher.classify(ordinary.as_view()),
                SlotMatch::Other
            ));
        }
        assert_eq!(
            matcher
                .parse::<LibraryRep, AbstractIndex>(slot.as_view())
                .unwrap(),
            Slot::<LibraryRep>::try_from(slot.as_view()).unwrap()
        );
    }

    #[test]
    fn cached_slot_parsing_preserves_validation_and_duality() {
        let rep = LibraryRep::from(Lorentz {});
        let head = rep.symbol();
        let index = symbol!("slot_match_tests::typed_mu");
        let slot = function!(head, 4, index);
        let dual = function!(AIND_SYMBOLS.dind, &slot);
        let mut matcher = SlotMatcher::default();

        for (atom, accepted) in [
            (slot.clone(), true),
            (dual.clone(), true),
            (function!(AIND_SYMBOLS.dind, &slot, 99), false),
            (function!(head, symbol!("slot_match_tests::D"), index), true),
            (function!(head), false),
            (function!(head, 4), false),
            (function!(head, -1, index), false),
            (function!(head, 4, index, 1), false),
            (function!(head, 4, index, function!(head, 4)), false),
            (function!(head, 4, function!(head, 4)), false),
            (
                function!(symbol!("slot_match_tests::unknown"), 4, index),
                false,
            ),
            (
                function!(
                    AIND_SYMBOLS.dind,
                    function!(LibraryRep::from(Minkowski {}).symbol(), 4, index)
                ),
                true,
            ),
        ] {
            let expected = Slot::<LibraryRep>::try_from(atom.as_view());
            assert_eq!(expected.is_ok(), accepted, "{atom}: {expected:?}");
            let actual = matcher.parse::<LibraryRep, AbstractIndex>(atom.as_view());
            match (expected, actual) {
                (Ok(expected), Ok(actual)) => assert_eq!(actual, expected, "{atom}"),
                (Err(expected), Err(actual)) => {
                    assert_eq!(actual.to_string(), expected.to_string(), "{atom}")
                }
                (expected, actual) => panic!("{atom}: expected {expected:?}, got {actual:?}"),
            }
        }
        let functional_dimension = function!(
            head,
            function!(symbol!("slot_match_tests::dimension"), 4),
            index
        );
        assert!(matches!(
            matcher.classify(functional_dimension.as_view()),
            SlotMatch::Explicit(_)
        ));
        assert!(matches!(
            matcher.parse::<LibraryRep, AbstractIndex>(functional_dimension.as_view()),
            Err(SlotError::DimErr(_))
        ));
        let actual = matcher
            .parse::<LibraryRep, AbstractIndex>(dual.as_view())
            .unwrap();
        assert_eq!(actual.rep_name(), rep.dual());
        assert!(
            matcher
                .parse::<Lorentz, AbstractIndex>(dual.as_view())
                .is_err()
        );
        assert!(
            matcher
                .parse::<DualLorentz, AbstractIndex>(dual.as_view())
                .is_ok()
        );
    }

    #[test]
    fn cached_slot_parsing_hydrates_portable_representations() {
        let mut matcher = SlotMatcher::default();
        for (name, class) in [
            ("ImportedSelfDualSlot", RepresentationClass::SelfDual),
            ("ImportedDualSlot", RepresentationClass::Dualizable),
            ("ImportedInlineSlot", RepresentationClass::InlineMetric),
        ] {
            let class_tag = if class == RepresentationClass::Dualizable {
                SPENSO_TAG.dualizable.clone()
            } else {
                SPENSO_TAG.self_dual.clone()
            };
            let name = format!("slot_match_tests::{name}");
            let head = SymbolBuilder::new(NamespacedSymbol::parse(&name))
                .with_tags(vec![SPENSO_TAG.representation.clone(), class_tag])
                .with_user_data(
                    RepresentationMetadata {
                        class,
                        label: IndexDisplay::symbol(name.rsplit("::").next().unwrap()).unwrap(),
                        index_palette: IndexPalette::Numeric,
                        index_row: IndexRow::Top,
                    }
                    .to_user_data(),
                )
                .build()
                .unwrap();
            let atom = function!(head, 4, 17);
            assert!(matches!(
                matcher.classify(atom.as_view()),
                SlotMatch::Explicit(_)
            ));
            let parsed = matcher.parse::<LibraryRep, AbstractIndex>(atom.as_view());
            if class == RepresentationClass::InlineMetric {
                assert!(matches!(parsed, Err(SlotError::RepError(
                    RepresentationError::ImportedInlineMetricRequiresLocalRegistration(symbol)
                )) if symbol == head));
            } else {
                let parsed = parsed.unwrap();
                assert_eq!(parsed.rep_name().symbol(), head);
                assert_eq!(
                    matcher
                        .parse::<LibraryRep, AbstractIndex>(atom.as_view())
                        .unwrap(),
                    parsed
                );
            }
        }
    }

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
        assert_snapshot!(mu.to_atom().to_canonical_string(), @"spenso::{spenso::representation,spenso::dualizable}::lor(4,0)");
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
