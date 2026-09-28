use super::{
    abstract_index::{AbstractIndex, AbstractIndexError},
    dimension::DimensionError,
    representation::{
        BaseRepName, LibraryRep, LibrarySlot, RepName, Representation, RepresentationError,
    },
};
use crate::structure::dimension::Dimension;
#[cfg(feature = "shadowing")]
use crate::{
    network::tags::{SPENSO_TAG, SpensoTags},
    structure::abstract_index::AIND_SYMBOLS,
};
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
        Atom, AtomView, FunctionBuilder, Symbol,
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
    fn from_symbol(head: Symbol, representation_tag: &str, wrappers: &[u32; 3]) -> Self {
        if head.has_tag(representation_tag) {
            Self::Representation
        } else if wrappers.contains(&head.get_id()) {
            Self::Wrapper
        } else {
            Self::Other
        }
    }
}

#[cfg(feature = "shadowing")]
/// Borrowed representation syntax with an optional variance wrapper.
/// Dimension and index payloads remain symbolic rather than being coerced.
#[derive(Clone, Copy, Debug)]
pub struct RepresentationView<'a> {
    function: FunView<'a>,
    wrapper: Option<Symbol>,
    arguments: ListIterator<'a>,
}

#[cfg(feature = "shadowing")]
enum RepresentationMatch<'a> {
    Recognized(RepresentationView<'a>),
    Opaque,
    Other,
}

#[cfg(feature = "shadowing")]
/// Borrowed explicit slot syntax, containing exactly a dimension and an index.
#[derive(Clone, Copy, Debug)]
pub struct SlotView<'a> {
    representation: RepresentationView<'a>,
}

#[cfg(feature = "shadowing")]
/// Classification of an expression for a structural tensor-index walk.
#[derive(Clone, Copy)]
pub enum SlotMatch<'a> {
    Explicit(SlotView<'a>),
    /// A compact or malformed slot: its payload is opaque to tensor-index scans.
    Opaque,
    Other,
}

#[cfg(feature = "shadowing")]
impl<'a> SlotMatch<'a> {
    /// Validate this classification without repeating the structural scan.
    /// Index conversion still runs for every occurrence and requested index type.
    pub(crate) fn parse<T: RepName, Aind: ParseableAind>(
        self,
        value: AtomView<'a>,
        matcher: &mut SlotMatcher,
    ) -> Result<Slot<T, Aind>, SlotError> {
        let slot = self.into_slot(value)?;
        let rep = matcher.representation(slot)?;
        slot.parse(T::from_library_rep(rep)?)
    }

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
impl<'a> RepresentationView<'a> {
    /// Decode arguments once after recognizing a representation or duality wrapper.
    #[inline]
    fn classify(
        value: AtomView<'a>,
        mut classify_head: impl FnMut(FunView<'a>) -> SlotHead,
    ) -> RepresentationMatch<'a> {
        let AtomView::Fun(mut function) = value else {
            return RepresentationMatch::Other;
        };
        let mut wrapper = None;
        match classify_head(function) {
            SlotHead::Other => return RepresentationMatch::Other,
            SlotHead::Representation => {}
            SlotHead::Wrapper => {
                let mut arguments = function.iter();
                if arguments.len() != 1 {
                    return RepresentationMatch::Opaque;
                }
                let Some(AtomView::Fun(inner)) = arguments.next() else {
                    return RepresentationMatch::Opaque;
                };
                if !matches!(classify_head(inner), SlotHead::Representation) {
                    return RepresentationMatch::Opaque;
                }
                wrapper = Some(function.get_symbol());
                function = inner;
            }
        }
        let arguments = function.iter();
        if !matches!(arguments.len(), 1 | 2) {
            return RepresentationMatch::Opaque;
        }
        RepresentationMatch::Recognized(Self {
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

    /// The representation head, without any variance wrapper.
    pub fn head(self) -> Symbol {
        self.function.get_symbol()
    }

    pub fn wrapper(self) -> Option<Symbol> {
        self.wrapper
    }

    /// Whether the syntax denotes the base orientation of its representation.
    pub fn is_base(self) -> bool {
        self.wrapper.is_none()
            || (self.wrapper == Some(AIND_SYMBOLS.uind)
                && self.head().has_tag(&SPENSO_TAG.dualizable))
            || (self.wrapper == Some(AIND_SYMBOLS.selfdualind)
                && self.head().has_tag(&SPENSO_TAG.self_dual))
    }

    pub fn is_self_dual(self) -> bool {
        self.is_base() && self.head().has_tag(&SPENSO_TAG.self_dual)
    }

    pub fn is_dual(self) -> bool {
        self.wrapper == Some(AIND_SYMBOLS.dind) && self.head().has_tag(&SPENSO_TAG.dualizable)
    }

    /// Compare exact symbolic dimensions and complementary variance using tags.
    /// This also supports tagged representations that have no library registration.
    pub fn matches(self, other: Self) -> bool {
        self.head() == other.head()
            && self.dimension() == other.dimension()
            && ((self.is_self_dual() && other.is_self_dual())
                || (self.head().has_tag(&SPENSO_TAG.dualizable)
                    && ((self.is_base() && other.is_dual())
                        || (self.is_dual() && other.is_base()))))
    }

    /// Build the compact base representation, retaining its exact dimension.
    pub fn compact(self) -> Atom {
        FunctionBuilder::new(self.head())
            .add_arg(self.dimension())
            .finish()
    }
}

#[cfg(feature = "shadowing")]
impl<'a> SlotView<'a> {
    #[inline]
    fn classify(
        value: AtomView<'a>,
        classify_head: impl FnMut(FunView<'a>) -> SlotHead,
    ) -> SlotMatch<'a> {
        match RepresentationView::classify(value, classify_head) {
            RepresentationMatch::Recognized(representation)
                if representation.arguments.len() == 2 =>
            {
                SlotMatch::Explicit(Self { representation })
            }
            RepresentationMatch::Other => SlotMatch::Other,
            _ => SlotMatch::Opaque,
        }
    }

    pub fn representation(self) -> RepresentationView<'a> {
        self.representation
    }

    /// Borrow the dimension expression without requiring a number or symbol.
    #[inline]
    pub fn dimension(self) -> AtomView<'a> {
        self.representation.dimension()
    }

    #[inline]
    /// Borrow the original explicit index payload without coercing its value.
    pub fn index(mut self) -> AtomView<'a> {
        self.representation.arguments.next();
        self.representation.arguments.next().unwrap()
    }

    /// Concrete component markers are not summed abstract index labels.
    pub fn is_concrete_index(self) -> bool {
        matches!(self.index(), AtomView::Fun(index)
            if [AIND_SYMBOLS.cind, AIND_SYMBOLS.find].contains(&index.get_symbol()))
    }

    fn parse<T: RepName, Aind: ParseableAind>(
        mut self,
        rep: T,
    ) -> Result<Slot<T, Aind>, SlotError> {
        let dim = Dimension::try_from(self.representation.arguments.next().unwrap())?;
        let aind =
            Aind::from_view(self.representation.arguments.next().unwrap()).map_err(Into::into)?;
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
pub struct SlotMatcher {
    heads: [Option<(u32, SlotHead)>; 16],
    resolved: Vec<(Symbol, Option<Symbol>, LibraryRep)>,
    tags: &'static SpensoTags,
    aind: Symbol,
    wrappers: [u32; 3],
    metric: std::sync::OnceLock<Symbol>,
    projectors: [std::sync::OnceLock<Symbol>; 3],
}

#[cfg(feature = "shadowing")]
impl Default for SlotMatcher {
    fn default() -> Self {
        // Force each symbol bundle once per scan, including when many distinct
        // tensor heads miss the small recognition cache.
        let wrappers = &*AIND_SYMBOLS;
        Self {
            metric: std::sync::OnceLock::new(),
            projectors: std::array::from_fn(|_| std::sync::OnceLock::new()),
            heads: [None; 16],
            resolved: Vec::new(),
            tags: &SPENSO_TAG,
            aind: wrappers.aind,
            wrappers: [wrappers.dind, wrappers.uind, wrappers.selfdualind]
                .map(|symbol| symbol.get_id()),
        }
    }
}

#[cfg(feature = "shadowing")]
impl SlotMatcher {
    pub(crate) fn index_bundle(&self) -> Symbol {
        self.aind
    }

    /// Bundle values stay local to one observation; the first access retains
    /// the ordinary initialization boundary and a new matcher probes again.
    pub(crate) fn metric(&self) -> Symbol {
        *self
            .metric
            .get_or_init(|| crate::network::library::symbolic::ETS.metric)
    }

    fn projector(&self, index: usize) -> Symbol {
        *self.projectors[index].get_or_init(|| match index {
            0 => *crate::shadowing::SYM,
            1 => *crate::shadowing::ANTISYM,
            2 => *crate::shadowing::CYCLIC,
            _ => unreachable!("three projector symbols"),
        })
    }

    pub(crate) fn is_projector(&self, symbol: Symbol) -> bool {
        (0..3).any(|index| self.projector(index) == symbol)
    }

    pub(crate) fn projectors(&self) -> [Symbol; 3] {
        std::array::from_fn(|index| self.projector(index))
    }

    pub(crate) fn is_variance_wrapper(&self, symbol: Symbol) -> bool {
        self.wrappers.contains(&symbol.get_id())
    }

    pub(crate) fn tags(&self) -> &'static SpensoTags {
        self.tags
    }

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
        let kind = SlotHead::from_symbol(
            function.get_symbol(),
            &self.tags.representation,
            &self.wrappers,
        );
        *entry = Some((id, kind));
        kind
    }

    /// Identify a slot or an opaque slot payload in one classification pass.
    #[inline]
    pub fn classify<'a>(&mut self, value: AtomView<'a>) -> SlotMatch<'a> {
        SlotView::classify(value, |function| self.classify_head(function))
    }

    /// Recognize exactly `rep(dimension)`, optionally inside one variance wrapper.
    pub fn compact_representation<'a>(
        &mut self,
        value: AtomView<'a>,
    ) -> Option<RepresentationView<'a>> {
        match RepresentationView::classify(value, |function| self.classify_head(function)) {
            RepresentationMatch::Recognized(representation)
                if representation.arguments.len() == 1 =>
            {
                Some(representation)
            }
            _ => None,
        }
    }

    /// Borrow the final argument of a canonical vector-shaped function.
    /// Earlier scalar metadata is opaque, but direct representation arguments
    /// are rejected. Callers decide whether the head must carry a rank-one tag.
    pub fn vector_argument<'a>(&mut self, function: FunView<'a>) -> Option<AtomView<'a>> {
        if function.get_symbol().is_scalar()
            || !matches!(self.classify(function.as_view()), SlotMatch::Other)
        {
            return None;
        }
        let mut arguments = function.iter();
        let mut last = arguments.next()?;
        for argument in arguments {
            if !matches!(self.classify(last), SlotMatch::Other) {
                return None;
            }
            last = argument;
        }
        Some(last)
    }

    /// Read a canonical singleton concrete component marker. It consumes no
    /// abstract slot. A recognized malformed marker is distinguished from
    /// ordinary scalar metadata so typed rank-one admission cannot hide it.
    pub fn concrete_component(&self, value: AtomView<'_>) -> Option<Result<usize, SlotError>> {
        let AtomView::Fun(function) = value else {
            return None;
        };
        if ![AIND_SYMBOLS.cind, AIND_SYMBOLS.find].contains(&function.get_symbol()) {
            return None;
        }
        Some(if function.get_nargs() == 1 {
            usize::try_from(function.get(0)).map_err(|_| SlotError::NotNatural)
        } else {
            Err(SlotError::NotNatural)
        })
    }

    /// Validate a recognized slot using the requested representation and index types.
    /// Use [`SlotView::index`] when exact symbolic index identity must be preserved.
    pub fn parse<T: RepName, Aind: ParseableAind>(
        &mut self,
        value: AtomView<'_>,
    ) -> Result<Slot<T, Aind>, SlotError> {
        self.classify(value).parse(value, self)
    }

    /// Parse a compact representation with the same recognition and resolution
    /// caches used for explicit slots. Indexed or malformed arguments are not
    /// compact ports.
    pub fn parse_representation<T: RepName>(
        &mut self,
        value: AtomView<'_>,
    ) -> Result<Representation<T>, SlotError> {
        let representation = self
            .compact_representation(value)
            .ok_or(SlotError::NotRepresentation)?;
        let rep = T::from_library_rep(
            self.resolve_representation(representation.head(), representation.wrapper())?,
        )?;
        let dim = Dimension::try_from(representation.dimension())?;
        Ok(Representation { dim, rep })
    }

    /// Resolve the representation and duality while retaining arbitrary dimension
    /// and index expressions in the borrowed view.
    pub fn representation(
        &mut self,
        slot: SlotView<'_>,
    ) -> Result<LibraryRep, RepresentationError> {
        self.resolve_representation(slot.representation.head(), slot.representation.wrapper())
    }

    /// Test the existing representation-prefix grammar without resolving heads
    /// which cannot denote either a representation or a variance wrapper.
    /// The full parser still owns dimensions, portable payloads and trailing args.
    pub(crate) fn is_representation(&mut self, value: AtomView<'_>) -> bool {
        matches!(value, AtomView::Fun(function)
            if !matches!(self.classify_head(function), SlotHead::Other))
            && self.representation_from_atom(value).is_ok()
    }

    /// Parse the representation prefix accepted by `Representation::try_from`,
    /// sharing this operation's successful head and variance resolutions.
    pub(crate) fn representation_from_atom(
        &mut self,
        value: AtomView<'_>,
    ) -> Result<Representation<LibraryRep>, SlotError> {
        Representation::parse_with(value, |head, wrapper| {
            self.resolve_representation(head, wrapper)
        })
    }

    fn resolve_representation(
        &mut self,
        head: Symbol,
        wrapper: Option<Symbol>,
    ) -> Result<LibraryRep, RepresentationError> {
        let rep = if let Some((_, _, rep)) = self
            .resolved
            .iter()
            .find(|(symbol, cached_wrapper, _)| *symbol == head && *cached_wrapper == wrapper)
        {
            *rep
        } else {
            let rep = match wrapper {
                Some(wrapper) => LibraryRep::try_from_symbol(head, wrapper)?,
                None => LibraryRep::try_from_symbol_coerced(head)?,
            };
            let entry = (head, wrapper, rep);
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
        let slot = SlotMatcher::default().classify(value).into_slot(value)?;
        let head = slot.representation.head();
        let rep = match slot.representation.wrapper() {
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
    fn cached_representation_prefix_preserves_direct_conversion_grammar() {
        let rep = LibraryRep::from(Lorentz {}).symbol();
        let unknown = symbol!("slot_match_tests::unknown_representation_prefix");
        let dimension = symbol!("slot_match_tests::prefix_dimension");
        let mut matcher = SlotMatcher::default();
        // Prefix conversion intentionally accepts trailing arguments. The strict
        // compact-port recognizer has a different contract and cannot replace it.
        for (atom, accepted) in [
            (function!(rep, 4), true),
            (function!(rep, dimension), true),
            (function!(rep, 4, 7), true),
            (function!(rep, 4, 7, 8), true),
            (function!(AIND_SYMBOLS.dind, function!(rep, 4)), true),
            (function!(AIND_SYMBOLS.dind, function!(rep, 4, 7), 8), true),
            (function!(rep), false),
            (function!(rep, -1), false),
            (function!(rep, function!(unknown, 4)), false),
            (function!(unknown, 4), false),
            (function!(unknown, function!(rep, 4)), false),
            (Atom::num(4), false),
        ] {
            for _ in 0..2 {
                let expected = Representation::<LibraryRep>::try_from(atom.as_view());
                let actual = matcher.representation_from_atom(atom.as_view());
                assert_eq!(expected.is_ok(), accepted, "{atom}: {expected:?}");
                assert_eq!(
                    matcher.is_representation(atom.as_view()),
                    accepted,
                    "{atom}"
                );
                match (expected, actual) {
                    (Ok(expected), Ok(actual)) => assert_eq!(actual, expected, "{atom}"),
                    (Err(expected), Err(actual)) => {
                        assert_eq!(actual.to_string(), expected.to_string(), "{atom}")
                    }
                    (expected, actual) => panic!("{atom}: expected {expected:?}, got {actual:?}"),
                }
            }
        }
    }

    #[test]
    fn representation_predicate_skips_impossible_vector_heads() {
        let vector = crate::vector_symbol!("slot_match_tests::representation_probe_vector");
        let custom = LibraryRep::new_dual("slot_match_tests::predicate_custom_rep").unwrap();
        let mut matcher = SlotMatcher::default();
        for rep in [
            LibraryRep::from(Minkowski {}),
            LibraryRep::from(Lorentz {}),
            custom,
        ] {
            let compact = rep.to_symbolic([Atom::num(4)]);
            let value = function!(vector, &compact);
            assert!(Representation::<LibraryRep>::try_from(value.as_view()).is_err());
            assert!(!matcher.is_representation(value.as_view()));
            // This boolean rejection must not resolve a compact vector's port
            // as though the vector head were a variance wrapper.
            assert!(matcher.resolved.is_empty());
        }
        for rep in [
            LibraryRep::from(Minkowski {}),
            LibraryRep::from(Lorentz {}),
            custom,
        ] {
            for atom in [
                rep.to_symbolic([Atom::num(4)]),
                rep.to_symbolic([Atom::num(4), Atom::num(7), Atom::num(8)]),
                rep.dual().to_symbolic([Atom::num(4)]),
            ] {
                assert_eq!(
                    matcher.is_representation(atom.as_view()),
                    Representation::<LibraryRep>::try_from(atom.as_view()).is_ok(),
                    "{atom}"
                );
            }
        }
    }

    #[test]
    fn reused_classification_still_parses_each_custom_index_occurrence() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        static CALLS: AtomicUsize = AtomicUsize::new(0);
        #[derive(Clone, Copy, Debug, PartialEq)]
        struct Index(usize);
        impl super::ParseableAind for Index {
            type Error = SlotError;
            fn from_view(value: symbolica::atom::AtomView<'_>) -> Result<Self, Self::Error> {
                CALLS.fetch_add(1, Ordering::SeqCst);
                usize::try_from(value)
                    .map(Self)
                    .map_err(|_| SlotError::NotNatural)
            }
            fn to_atom(&self) -> Atom {
                Atom::num(self.0)
            }
        }
        let atom = function!(LibraryRep::from(Minkowski {}).symbol(), 4, 7);
        let mut matcher = SlotMatcher::default();
        let recognized = matcher.classify(atom.as_view());
        CALLS.store(0, Ordering::SeqCst);
        for _ in 0..2 {
            let slot = recognized
                .parse::<LibraryRep, Index>(atom.as_view(), &mut matcher)
                .unwrap();
            assert_eq!(slot.aind, Index(7));
        }
        assert_eq!(CALLS.load(Ordering::SeqCst), 2);
    }

    #[test]
    fn slot_matcher_remains_send_and_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<SlotMatcher>();
    }

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
    fn cached_compact_representations_preserve_variance_and_reject_extra_arguments() {
        let head = LibraryRep::from(Lorentz {}).symbol();
        let compact = function!(head, 4);
        let explicit = function!(head, 4, symbol!("slot_match_tests::compact_mu"));
        let mut matcher = SlotMatcher::default();
        for atom in [
            compact.clone(),
            function!(AIND_SYMBOLS.dind, &compact),
            function!(AIND_SYMBOLS.uind, &compact),
            // Symbolica removes a double dual before the matcher sees it.
            function!(AIND_SYMBOLS.dind, function!(AIND_SYMBOLS.dind, &compact)),
            function!(head, symbol!("slot_match_tests::compact_D")),
            function!(LibraryRep::from(Minkowski {}).symbol(), 4),
        ] {
            let expected = Representation::<LibraryRep>::try_from(atom.as_view()).unwrap();
            for _ in 0..2 {
                assert_eq!(
                    matcher
                        .parse_representation::<LibraryRep>(atom.as_view())
                        .unwrap(),
                    expected,
                    "{atom}"
                );
                // Explicit and compact parsing share the cache; alternating them
                // must not confuse the variance or the argument count.
                assert_eq!(
                    matcher
                        .parse::<LibraryRep, AbstractIndex>(explicit.as_view())
                        .unwrap()
                        .rep(),
                    Representation::<LibraryRep>::try_from(compact.as_view()).unwrap()
                );
            }
        }
        for atom in [
            explicit,
            function!(head),
            function!(head, -1),
            function!(
                head,
                function!(symbol!("slot_match_tests::compact_dimension"), 4)
            ),
            function!(AIND_SYMBOLS.dind, &compact, 1),
            function!(symbol!("slot_match_tests::compact_unknown"), 4),
        ] {
            assert!(
                matcher
                    .parse_representation::<LibraryRep>(atom.as_view())
                    .is_err(),
                "{atom}"
            );
        }
        assert!(
            matcher
                .parse_representation::<Lorentz>(function!(AIND_SYMBOLS.dind, &compact).as_view())
                .is_err()
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
            // The prefix predicate must hydrate genuine portable heads through
            // the existing resolver, before an explicit-slot parse uses them.
            assert_eq!(
                matcher.is_representation(atom.as_view()),
                class != RepresentationClass::InlineMetric
            );
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
    #[test]
    fn canonical_singleton_components_are_checked_without_abstract_slots() {
        let matcher = SlotMatcher::default();
        for marker in [AIND_SYMBOLS.cind, AIND_SYMBOLS.find] {
            for index in [0, 3] {
                let component = function!(marker, Atom::num(index));
                assert_eq!(
                    matcher
                        .concrete_component(component.as_view())
                        .unwrap()
                        .unwrap(),
                    index as usize
                );
            }
            for component in [
                symbolica::atom::FunctionBuilder::new(marker).finish(),
                function!(marker, -1),
                function!(marker, Atom::num(1) / Atom::num(2)),
                function!(marker, symbol!("component_slot::unknown")),
                function!(marker, 0, 1),
            ] {
                assert!(
                    matcher
                        .concrete_component(component.as_view())
                        .unwrap()
                        .is_err(),
                    "{component}"
                );
            }
        }
        assert!(
            matcher
                .concrete_component(function!(symbol!("component_slot::cind"), 0).as_view())
                .is_none()
        );
        assert!(matcher.concrete_component(Atom::num(0).as_view()).is_none());
    }
}
