use super::{
    abstract_index::{AbstractIndex, AbstractIndexError},
    concrete_index::ConcreteIndex,
    dimension::{Dimension, DimensionError},
    slot::Slot,
};
use ahash::AHashMap;
use append_only_vec::AppendOnlyVec;
use linnet::half_edge::involution::Orientation;
use serde::{Deserialize, Serialize};
use spenso_macros::SimpleRepresentation;
use std::ops::Index;
use std::{
    cmp::Ordering,
    convert::Infallible,
    fmt::{Debug, Display},
    hash::{Hash, Hasher},
    sync::{LazyLock, RwLock},
};

use bincode::{Decode, Encode};

#[cfg(feature = "shadowing")]
use crate::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    self_dual_symbol,
    structure::slot::SlotError,
};

#[cfg(feature = "shadowing")]
use symbolica::{
    atom::{
        Atom, AtomOrView, AtomView, FunctionBuilder, NamespacedSymbol, Symbol, SymbolBuilder,
        UserData,
    },
    coefficient::CoefficientView,
    function,
    printer::PrintUserData,
    symbol,
};

use thiserror::Error;

use eyre::Result;

pub trait BaseRepName: RepName<Dual: RepName> + Default {
    const NAME: &'static str;
    // fn selfless_name() -> String;
    fn selfless_base() -> Self::Base;

    #[cfg(feature = "shadowing")]
    fn selfless_symbol() -> Symbol {
        symbol!(Self::NAME)
    }
    fn selfless_dual() -> Self::Dual;
    fn selfless_rep<D: Into<Dimension>>(dim: D) -> Representation<Self>
    where
        Self: Sized,
    {
        Representation {
            dim: dim.into(),
            rep: Self::default(),
        }
    }

    #[cfg(feature = "shadowing")]
    fn pattern(symbol: Symbol) -> Atom {
        Self::default().to_symbolic([Atom::var(symbol)])
    }

    fn slot<D: Into<Dimension>, A: Into<AbstractIndex>>(dim: D, aind: A) -> Slot<Self>
    where
        Self: Sized,
    {
        let aind: AbstractIndex = aind.into();
        Slot {
            rep: Self::selfless_rep(dim),
            aind,
        }
    }

    #[cfg(feature = "shadowing")]
    fn new_slot_from(
        sym: Symbol,
        dim: Dimension,
        aind: AbstractIndex,
    ) -> Result<Slot<Self>, SlotError> {
        if sym == Self::selfless_symbol() {
            ::std::result::Result::Ok(Slot {
                rep: Representation {
                    dim,
                    rep: Self::default(),
                },
                aind,
            })
        } else {
            Err(SlotError::NotRepresentation)
        }
    }
}

#[derive(Error, Debug)]
pub enum RepresentationError {
    #[cfg(feature = "shadowing")]
    #[error("Symbol {0} isn't one of [sind,uind,dind]")]
    SymbolError(Symbol),
    #[cfg(feature = "shadowing")]
    #[error("Expected dual state: {0} but got {1}")]
    ExpectedDualStateError(Symbol, Symbol),
    #[cfg(feature = "shadowing")]
    #[error("{0} is not a possible Representation")]
    NotRepresentationError(Symbol),
    #[cfg(feature = "shadowing")]
    #[error(
        "inline-metric representation {0} must be registered locally because its metric function is not portable"
    )]
    ImportedInlineMetricRequiresLocalRegistration(Symbol),
    #[error("representation library error: {0}")]
    RepLibrary(#[from] RepLibraryError),
    #[error("Wrong representation, expected {0},got {1}")]
    WrongRepresentationError(String, String),
    #[error("Abstract index error :{0}")]
    AindError(#[from] AbstractIndexError),
    #[error("{0}")]
    DimErr(#[from] DimensionError),
    #[error("{0}")]
    Any(#[from] eyre::Error),
    #[error("infallible")]
    Infallible(#[from] Infallible),
}

pub trait RepName:
    Copy + Clone + Debug + PartialEq + Eq + Hash + Display + Ord + Into<LibraryRep>
{
    type Dual: RepName<Dual = Self, Base = Self::Base>;
    type Base: RepName;

    fn from_library_rep(rep: LibraryRep) -> Result<Self, RepresentationError>;
    fn is_dummy(self) -> bool {
        false
    }
    fn orientation(self) -> Orientation;
    fn dual(self) -> Self::Dual;
    fn is_dual(self) -> bool;
    fn base(&self) -> Self::Base;
    fn is_base(&self) -> bool;
    fn is_self_dual(&self) -> bool {
        self.is_base() && self.is_dual()
    }

    fn matches(&self, other: &Self::Dual) -> bool;

    fn match_cmp(&self, _other: &Self::Dual) -> Ordering {
        Ordering::Equal
    }

    #[cfg(feature = "shadowing")]
    fn try_from_symbol(sym: Symbol, aind: Symbol) -> Result<Self, RepresentationError> {
        Self::from_library_rep(LibraryRep::try_from_symbol(sym, aind)?)
    }
    #[cfg(feature = "shadowing")]
    fn try_from_symbol_coerced(sym: Symbol) -> Result<Self, RepresentationError> {
        Self::from_library_rep(LibraryRep::try_from_symbol_coerced(sym)?)
    }

    // fn try_from<B: BaseRepName>(b: B) -> Result<B, SlotError>;

    /// for the given concrete index, says whether it should have a minus sign during contraction
    ///
    /// for example see [`Self::negative`]
    #[must_use]
    fn is_neg(self, _i: usize) -> bool {
        false
    }

    #[allow(clippy::cast_possible_wrap)]
    #[cfg(feature = "shadowing")]
    /// yields a function builder for the representation, adding a first variable: the dimension.
    ///
    fn to_symbolic<'a, It: Into<AtomOrView<'a>>>(
        &self,
        args: impl IntoIterator<Item = It>,
    ) -> Atom {
        let librep: LibraryRep = (*self).into();

        librep.to_symbolic(args)
    }

    #[cfg(feature = "shadowing")]
    /// An atom representing the identity function for that representation.
    /// a is dualized
    /// b is not
    fn id_atom<'a, It: Into<AtomOrView<'a>>>(
        &self,
        a: impl IntoIterator<Item = It>,
        b: impl IntoIterator<Item = It>,
    ) -> Atom {
        let librep: LibraryRep = (*self).into();
        function!(
            ETS.metric,
            librep.dual().to_symbolic(a),
            librep.to_symbolic(b)
        )
    }

    #[allow(clippy::cast_possible_wrap)]
    #[cfg(feature = "shadowing")]
    /// yields a function builder for the representation, adding a first variable: the dimension.
    ///
    fn metric_atom<'a, It: Into<AtomOrView<'a>>>(
        &self,
        a: impl IntoIterator<Item = It>,
        b: impl IntoIterator<Item = It>,
    ) -> Atom {
        let librep: LibraryRep = (*self).into();
        function!(ETS.metric, librep.to_symbolic(a), librep.to_symbolic(b))
    }

    fn new_slot<Aind, D: Into<Dimension>, A: Into<Aind>>(self, dim: D, aind: A) -> Slot<Self, Aind>
    where
        Self: Sized,
    {
        Slot {
            rep: self.new_rep(dim),
            aind: aind.into(),
        }
    }

    fn new_rep<D: Into<Dimension>>(&self, dim: D) -> Representation<Self>
    where
        Self: Sized,
    {
        Representation {
            dim: dim.into(),
            rep: *self,
        }
    }
}

#[rustfmt::skip]
#[derive(SimpleRepresentation)]
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Hash,
    PartialOrd,
    Ord,
    Default,
    Serialize,
    Deserialize,
    Encode,
    Decode,
)]
#[representation(name = "euc", self_dual)] // Specify the dual name
pub struct Euclidean {}

#[rustfmt::skip]
#[derive(SimpleRepresentation)]
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Hash,
    PartialOrd,
    Ord,
    Default,
    Serialize,
    Encode,
    Decode,
    Deserialize,
)]
#[representation(name = "lor")] // Specify the dual name
pub struct Lorentz {}

#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Hash,
    PartialOrd,
    Ord,
    Serialize,
    Deserialize,
    Default,
    Encode,
    Decode,
)]
pub struct Minkowski {}

impl From<Minkowski> for LibraryRep {
    fn from(_value: Minkowski) -> Self {
        ExtendibleReps::MINKOWSKI
    }
}

impl RepName for Minkowski {
    type Base = Minkowski;
    type Dual = Minkowski;

    fn from_library_rep(rep: LibraryRep) -> Result<Self, RepresentationError> {
        rep.try_into()
    }

    fn orientation(self) -> ::linnet::half_edge::involution::Orientation {
        ::linnet::half_edge::involution::Orientation::Undirected
    }

    fn base(&self) -> Self::Base {
        Minkowski::selfless_base()
    }

    fn is_base(&self) -> bool {
        true
    }

    fn dual(self) -> Self::Dual {
        Minkowski::selfless_dual()
    }

    fn is_dual(self) -> bool {
        true
    }

    fn matches(&self, _: &Self::Dual) -> bool {
        true
    }

    fn is_neg(self, i: usize) -> bool {
        i > 0
    }
}

impl Display for Minkowski {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "mink")
    }
}

impl TryFrom<LibraryRep> for Minkowski {
    type Error = RepresentationError;

    fn try_from(value: LibraryRep) -> std::result::Result<Self, Self::Error> {
        if value == ExtendibleReps::MINKOWSKI {
            std::result::Result::Ok(Minkowski {})
        } else {
            Err(RepresentationError::WrongRepresentationError(
                "mink".to_string(),
                value.to_string(),
            ))
        }
    }
}

impl BaseRepName for Minkowski {
    const NAME: &'static str = "mink";

    fn selfless_base() -> Self::Base {
        Self::default()
    }

    fn selfless_dual() -> Self::Dual {
        Self::default()
    }
}
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Hash,
    PartialOrd,
    Ord,
    Serialize,
    Deserialize,
    Default,
    Encode,
    Decode,
)]
pub struct Dummy {}

impl From<Dummy> for LibraryRep {
    fn from(_value: Dummy) -> Self {
        LibraryRep::Dummy
    }
}

impl RepName for Dummy {
    type Base = Dummy;
    type Dual = Dummy;

    fn from_library_rep(rep: LibraryRep) -> Result<Self, RepresentationError> {
        rep.try_into()
    }

    fn orientation(self) -> ::linnet::half_edge::involution::Orientation {
        ::linnet::half_edge::involution::Orientation::Undirected
    }

    fn base(&self) -> Self::Base {
        Dummy::selfless_base()
    }

    fn is_base(&self) -> bool {
        true
    }

    fn dual(self) -> Self::Dual {
        Dummy::selfless_dual()
    }

    fn is_dual(self) -> bool {
        true
    }

    fn matches(&self, _: &Self::Dual) -> bool {
        true
    }
}

impl Display for Dummy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "dummy")
    }
}

impl TryFrom<LibraryRep> for Dummy {
    type Error = RepresentationError;

    fn try_from(value: LibraryRep) -> std::result::Result<Self, Self::Error> {
        if value == LibraryRep::Dummy {
            std::result::Result::Ok(Dummy {})
        } else {
            Err(RepresentationError::WrongRepresentationError(
                "dummy".to_string(),
                value.to_string(),
            ))
        }
    }
}

impl BaseRepName for Dummy {
    const NAME: &'static str = "dummy";

    fn selfless_base() -> Self::Base {
        Self::default()
    }

    fn selfless_dual() -> Self::Dual {
        Self::default()
    }
}

#[derive(
    Debug,
    Copy,
    Clone,
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
pub struct Representation<T: RepName> {
    pub rep: T,
    pub dim: Dimension,
}

impl<T: RepName> Ord for Representation<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.rep.is_dummy() && other.rep.is_dummy() {
            Ordering::Equal
        } else {
            self.rep.cmp(&other.rep).then(self.dim.cmp(&other.dim))
        }
    }
}

impl<T: RepName> PartialOrd for Representation<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T: Hash + RepName> Hash for Representation<T> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.rep.hash(state);
        self.dim.hash(state);
    }
}

impl<T: RepName> PartialEq for Representation<T> {
    fn eq(&self, other: &Self) -> bool {
        self.rep == other.rep && self.dim == other.dim
    }
}

impl<T: RepName> Eq for Representation<T> {}

impl<T: RepName> Representation<T> {
    #[allow(clippy::cast_possible_wrap)]
    #[cfg(feature = "shadowing")]
    /// An atom representing the identity tensor with aind a, and b.
    /// a is dualized, b is not.
    ///
    pub fn id<'a, It: Into<AtomOrView<'a>>>(&self, a: It, b: It) -> Atom {
        let a: AtomOrView<'a> = a.into();
        let b: AtomOrView<'a> = b.into();
        function!(ETS.metric, self.dual().pattern(a), self.pattern(b))
    }

    // [allow(clippy::cast_possible_wrap)]
    #[cfg(feature = "shadowing")]
    /// An atom representing the identity tensor with aind a, and b.
    /// a is dualized, b is not.
    ///
    pub fn inner_product<'a, It: Into<AtomOrView<'a>>>(&self, a: It, b: It) -> Atom {
        fn with_rep(value: AtomView<'_>, rep: &Atom) -> Atom {
            match value {
                AtomView::Fun(fun) => {
                    let mut rebuilt = FunctionBuilder::new(fun.get_symbol());
                    for arg in fun.iter() {
                        rebuilt = rebuilt.add_arg(arg);
                    }
                    rebuilt.add_arg(rep).finish()
                }
                AtomView::Var(var) => FunctionBuilder::new(var.get_symbol()).add_arg(rep).finish(),
                _ => value.to_owned(),
            }
        }

        let a: AtomOrView<'a> = a.into();
        let b: AtomOrView<'a> = b.into();
        let rep = self.to_symbolic([]);
        function!(
            SPENSO_TAG.dot,
            with_rep(a.as_view(), &rep),
            with_rep(b.as_view(), &rep)
        )
    }

    pub fn base(self) -> Representation<T::Base> {
        Representation {
            dim: self.dim,
            rep: self.rep.base(),
        }
    }

    #[allow(clippy::cast_possible_wrap)]
    #[cfg(feature = "shadowing")]
    /// An atom representing the metric tensor with aind a, and b.
    pub fn g<'a, It: Into<AtomOrView<'a>>>(&self, a: It, b: It) -> Atom {
        let a: AtomOrView<'a> = a.into();
        let b: AtomOrView<'a> = b.into();
        function!(ETS.metric, self.pattern(a), self.pattern(b))
    }
    #[cfg(feature = "shadowing")]
    /// An atom representing the musical isomorphism tensor with aind a, and b.
    pub fn flat<'a, It: Into<AtomOrView<'a>>>(&self, a: It, b: It) -> Atom {
        let a: AtomOrView<'a> = a.into();
        let b: AtomOrView<'a> = b.into();
        function!(ETS.flat, self.pattern(a), self.pattern(b))
    }

    pub fn to_lib(self) -> Representation<LibraryRep> {
        let rep: LibraryRep = self.rep.into();
        Representation { dim: self.dim, rep }
    }

    pub fn to_dummy(self) -> Representation<Dummy> {
        Representation {
            dim: self.dim,
            rep: Dummy {},
        }
    }

    pub fn dot(&self) -> String {
        format!(
            "<<TABLE><TR><TD>{}</TD><TD>{}</TD></TR></TABLE>>",
            self.rep, self.dim
        )
    }

    #[cfg(feature = "shadowing")]
    pub fn pattern<'a, A: Into<AtomOrView<'a>>>(&self, aind: A) -> Atom {
        let dim = AtomOrView::Atom(self.dim.to_symbolic());
        let a = aind.into();
        self.rep.to_symbolic([dim, a])
    }

    #[cfg(feature = "shadowing")]
    pub fn to_pattern_wrapped(&self, aind: Symbol) -> Atom {
        self.rep.to_symbolic([
            self.dim.to_symbolic(),
            function!(symbol!("indexid"), Atom::var(aind)),
        ])
    }
}

#[derive(PartialEq, Eq, Clone, Copy, Debug, Hash, Serialize, Deserialize, Encode, Decode)]
pub enum LibraryRep {
    SelfDual(u16),
    InlineMetric(u16),
    Dualizable(i16),
    Dummy,
}

impl Ord for LibraryRep {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self, other) {
            (LibraryRep::SelfDual(a), LibraryRep::SelfDual(b)) => a.cmp(b),
            (LibraryRep::InlineMetric(a), LibraryRep::InlineMetric(b)) => a.cmp(b),
            (LibraryRep::Dualizable(a), LibraryRep::Dualizable(b)) => {
                if *a < 0 {
                    if *b < 0 {
                        a.abs().cmp(&b.abs())
                    } else {
                        Ordering::Greater
                    }
                } else if *b < 0 {
                    Ordering::Less
                } else {
                    a.cmp(b)
                }
            }
            (LibraryRep::SelfDual(_), LibraryRep::Dualizable(_))
            | (LibraryRep::SelfDual(_), LibraryRep::InlineMetric(_))
            | (LibraryRep::InlineMetric(_), LibraryRep::Dualizable(_)) => Ordering::Less,
            (LibraryRep::Dualizable(_), LibraryRep::SelfDual(_))
            | (LibraryRep::InlineMetric(_), LibraryRep::SelfDual(_))
            | (LibraryRep::Dualizable(_), LibraryRep::InlineMetric(_)) => Ordering::Greater,
            (LibraryRep::Dummy, LibraryRep::Dummy) => Ordering::Equal,
            (LibraryRep::Dummy, _) => Ordering::Less,
            (_, LibraryRep::Dummy) => Ordering::Greater,
        }
    }
}

#[test]
fn sorting_reps() {
    use linnet::permutation::Permutation;
    let mut a = [
        Euclidean {}.new_rep(4).cast(),
        Euclidean {}.new_rep(4).cast(),
        LibraryRep::from(Minkowski {}).new_rep(4),
    ];

    let perm = Permutation::sort(a);
    perm.apply_slice_in_place(&mut a);

    let mut b = [
        Euclidean {}.new_rep(4).cast(),
        LibraryRep::from(Minkowski {}).new_rep(4),
        Euclidean {}.new_rep(4).cast(),
    ];

    let perm = Permutation::sort(b);
    perm.apply_slice_in_place(&mut b);

    let mut c = [
        LibraryRep::from(Minkowski {}).new_rep(4),
        Euclidean {}.new_rep(4).cast(),
        Euclidean {}.new_rep(4).cast(),
    ];

    let perm = Permutation::sort(c);
    perm.apply_slice_in_place(&mut c);

    assert_eq!(a, b);
    assert_eq!(a, c);
    // assert_eq!()
}

impl PartialOrd for LibraryRep {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub type LibrarySlot<Aind> = Slot<LibraryRep, Aind>;

#[cfg(feature = "shadowing")]
crate::symbolica_init_lazy_static! {
    pub(crate) static REPS, REPS_INNER: RwLock<ExtendibleReps> =
        || RwLock::new(ExtendibleReps::new());
}
#[cfg(not(feature = "shadowing"))]
pub(crate) static REPS: LazyLock<RwLock<ExtendibleReps>> =
    LazyLock::new(|| RwLock::new(ExtendibleReps::new()));
pub(crate) static SELF_DUAL: AppendOnlyVec<(LibraryRep, RepData)> = AppendOnlyVec::new();
pub(crate) static INLINE_METRIC: AppendOnlyVec<(LibraryRep, MetricRepData)> = AppendOnlyVec::new();
pub(crate) static DUALIZABLE: AppendOnlyVec<(LibraryRep, RepData)> = AppendOnlyVec::new();

pub const LATIN: [char; 26] = [
    'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's',
    't', 'u', 'v', 'w', 'x', 'y', 'z',
];

pub const GREEK: [char; 24] = [
    'α', 'β', 'γ', 'δ', 'ε', 'ζ', 'η', 'θ', 'ι', 'κ', 'λ', 'μ', 'ν', 'ξ', 'ο', 'π', 'ρ', 'σ', 'τ',
    'υ', 'φ', 'χ', 'ψ', 'ω',
];

pub const CYRILLIC: [char; 33] = [
    'а', 'б', 'в', 'г', 'д', 'е', 'ё', 'ж', 'з', 'и', 'й', 'к', 'л', 'м', 'н', 'о', 'п', 'р', 'с',
    'т', 'у', 'ф', 'х', 'ц', 'ч', 'ш', 'щ', 'ъ', 'ы', 'ь', 'э', 'ю', 'я',
];

pub const CIRCLED: [char; 26] = [
    'ⓐ', 'ⓑ', 'ⓒ', 'ⓓ', 'ⓔ', 'ⓕ', 'ⓖ', 'ⓗ', 'ⓘ', 'ⓙ', 'ⓚ', 'ⓛ', 'ⓜ', 'ⓝ', 'ⓞ', 'ⓟ', 'ⓠ', 'ⓡ', 'ⓢ',
    'ⓣ', 'ⓤ', 'ⓥ', 'ⓦ', 'ⓧ', 'ⓨ', 'ⓩ',
];

pub const MATH_BOLD: [char; 26] = [
    '𝐚', '𝐛', '𝐜', '𝐝', '𝐞', '𝐟', '𝐠', '𝐡', '𝐢', '𝐣', '𝐤', '𝐥', '𝐦', '𝐧', '𝐨', '𝐩', '𝐪', '𝐫', '𝐬',
    '𝐭', '𝐮', '𝐯', '𝐰', '𝐱', '𝐲', '𝐳',
];

pub fn encode_base(mut n: usize, alphabet: &[char]) -> String {
    let base = alphabet.len();
    assert!(base > 0);
    let mut parts: Vec<char> = Vec::new();
    loop {
        parts.push(alphabet[n % base]);
        n /= base;
        if n == 0 {
            break;
        }
    }
    parts.into_iter().rev().collect::<String>()
}

fn canonical_representation_name(name: &str) -> String {
    if name.contains("::") {
        name.to_owned()
    } else {
        format!("spenso::{name}")
    }
}

fn representation_label_name(name: &str) -> &str {
    name.rsplit("::").next().unwrap_or(name)
}

fn default_index_row(name: &str) -> IndexRow {
    // Idenso's canonical self-dual bispinor indices conventionally sit below
    // the tensor head. A same-named representation in another namespace is a
    // separate declaration and retains the ordinary top-row default.
    if canonical_representation_name(name) == "spenso::bis" {
        IndexRow::Bottom
    } else {
        IndexRow::Top
    }
}

fn default_index_palette(name: &str) -> IndexPalette {
    let labels: &[&str] = match canonical_representation_name(name).as_str() {
        "spenso::mink" | "spenso::lor" => &["mu", "nu", "rho", "sigma"],
        "spenso::euc" | "spenso::cof" => &["i", "j", "k", "l"],
        "spenso::bis" | "spenso::coad" => &["a", "b", "c", "d"],
        "spenso::spf" => &["alpha", "beta", "gamma", "delta"],
        "spenso::cos" => &["I", "J", "K", "L"],
        _ => return IndexPalette::Numeric,
    };

    IndexPalette::cyclic(
        1,
        labels.iter().map(|label| {
            IndexDisplay::symbol(*label)
                .expect("the built-in representation palettes contain valid symbols")
        }),
    )
    .expect("the built-in representation palettes are non-empty and valid")
}

fn representation_typst_body(label: &IndexDisplay) -> String {
    format!(
        "(dim, ind ) = (content: $ {}^#dim _#ind $, upper:true)",
        label.to_typst_source()
    )
}

/// A safe, printer-neutral description of one mathematical index label.
///
/// The tree is deliberately small: plugins may deserialize it from Symbolica
/// symbol metadata and turn it into Typst source without evaluating arbitrary
/// source supplied by an Atom payload.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IndexDisplay {
    Symbol(String),
    /// Upright text embedded in a mathematical label.
    Text(String),
    Number(i64),
    Sequence(Vec<IndexDisplay>),
    /// A nested argument/row list. Unlike `Sequence`, entries retain their
    /// comma/row boundary when used by a safe mathematical display node.
    List(Vec<IndexDisplay>),
    /// A validated mathematical display node.
    ///
    /// `head` is never emitted directly. Rendering dispatches through a fixed
    /// allowlist, so imported symbol metadata cannot inject Typst source.
    Math {
        head: String,
        arguments: Vec<IndexDisplay>,
    },
    Attach {
        base: Box<IndexDisplay>,
        top: Option<Box<IndexDisplay>>,
        bottom: Option<Box<IndexDisplay>>,
    },
}

impl IndexDisplay {
    pub fn symbol(name: impl Into<String>) -> Result<Self, RepLibraryError> {
        let name = name.into();
        if name.is_empty() {
            return Err(RepLibraryError::InvalidIndexDisplay(
                "an index label cannot be empty".to_owned(),
            ));
        }
        if name.chars().count() > 128 || name.chars().any(char::is_control) {
            return Err(RepLibraryError::InvalidIndexDisplay(
                "an index label must contain at most 128 non-control characters".to_owned(),
            ));
        }
        Ok(Self::Symbol(name))
    }

    pub fn with_bottom(self, bottom: IndexDisplay) -> Self {
        match self {
            Self::Attach {
                base,
                top,
                bottom: None,
            } => Self::Attach {
                base,
                top,
                bottom: Some(Box::new(bottom)),
            },
            display => Self::Attach {
                base: Box::new(display),
                top: None,
                bottom: Some(Box::new(bottom)),
            },
        }
    }

    pub fn with_top(self, top: IndexDisplay) -> Self {
        match self {
            Self::Attach {
                base,
                top: None,
                bottom,
            } => Self::Attach {
                base,
                top: Some(Box::new(top)),
                bottom,
            },
            display => Self::Attach {
                base: Box::new(display),
                top: Some(Box::new(top)),
                bottom: None,
            },
        }
    }

    pub fn text(value: impl Into<String>) -> Result<Self, RepLibraryError> {
        let value = value.into();
        if value.chars().count() > 1024 || value.chars().any(char::is_control) {
            return Err(RepLibraryError::InvalidIndexDisplay(
                "display text must contain at most 1024 non-control characters".to_owned(),
            ));
        }
        Ok(Self::Text(value))
    }

    pub fn math(
        head: impl Into<String>,
        arguments: Vec<IndexDisplay>,
    ) -> Result<Self, RepLibraryError> {
        let head = head.into();
        if !matches!(
            head.as_str(),
            "arg"
                | "add"
                | "sub"
                | "plus"
                | "neg"
                | "times"
                | "dot"
                | "factorial"
                | "mul"
                | "()"
                | "group"
                | "lr"
                | "pow"
                | "frac"
                | "root"
                | "call"
                | "op-call"
                | "op"
                | "accent"
                | "underline"
                | "overline"
                | "underbrace"
                | "overbrace"
                | "underbracket"
                | "overbracket"
                | "underparen"
                | "overparen"
                | "undershell"
                | "overshell"
                | "cancel"
                | "vec"
                | "mat"
                | "cases"
                | "class"
                | "mid"
                | "scripts"
                | "limits"
                | "stretch"
        ) {
            return Err(RepLibraryError::InvalidIndexDisplay(format!(
                "unsupported mathematical display node {head:?}"
            )));
        }
        if arguments.len() > 64 {
            return Err(RepLibraryError::InvalidIndexDisplay(
                "a mathematical display node cannot contain more than 64 arguments".to_owned(),
            ));
        }
        Ok(Self::Math { head, arguments })
    }

    fn escaped_typst_string(value: &str) -> String {
        value
            .replace('\\', "\\\\")
            .replace('"', "\\\"")
            .replace('\n', "\\n")
            .replace('\r', "\\r")
    }

    fn literal_text(&self) -> Option<&str> {
        match self {
            Self::Symbol(value) | Self::Text(value) => Some(value),
            _ => None,
        }
    }

    fn typst_math_source(head: &str, arguments: &[Self]) -> String {
        let sources = arguments
            .iter()
            .map(Self::to_typst_source)
            .collect::<Vec<_>>();
        let joined = |separator: &str| sources.join(separator);
        let unary = |prefix: &str, suffix: &str| {
            sources
                .first()
                .map(|body| format!("{prefix}{body}{suffix}"))
                .unwrap_or_else(|| r#"upright("?")"#.to_owned())
        };
        let call = |name: &str| format!("{name}({})", joined(","));

        match head {
            "arg" => joined(","),
            "add" => joined(" + "),
            "sub" => joined(" - "),
            "plus" => unary("+", ""),
            "neg" => unary("-", ""),
            "times" => joined(" times "),
            "dot" => joined(" dot "),
            "factorial" => unary("", "!"),
            "mul" => joined(" "),
            "()" | "group" | "lr" => unary("lr((", "))"),
            "pow" if sources.len() == 2 => {
                format!("attach({},t:{})", sources[0], sources[1])
            }
            "frac" if sources.len() == 2 => format!("frac({},{})", sources[0], sources[1]),
            "root" if sources.len() == 1 => format!("root({})", sources[0]),
            "root" if sources.len() == 2 => {
                format!("root({},{})", sources[0], sources[1])
            }
            "call" | "op-call" if !sources.is_empty() => {
                let body = sources[1..].join(",");
                format!("{}({body})", sources[0])
            }
            "op" if arguments.len() == 1 => {
                let text = arguments[0].literal_text().unwrap_or("?");
                format!(r#"op("{}")"#, Self::escaped_typst_string(text))
            }
            "accent" if arguments.len() == 2 => {
                let accent = arguments[1].literal_text().unwrap_or("?");
                format!(
                    r#"accent({},"{}")"#,
                    sources[0],
                    Self::escaped_typst_string(accent)
                )
            }
            "underline" | "overline" | "cancel" | "mid" | "scripts" | "limits" | "stretch"
                if sources.len() == 1 =>
            {
                call(head)
            }
            "underbrace" | "overbrace" | "underbracket" | "overbracket" | "underparen"
            | "overparen" | "undershell" | "overshell"
                if (1..=2).contains(&sources.len()) =>
            {
                call(head)
            }
            "vec" => call("vec"),
            "mat" => {
                let rows = arguments
                    .iter()
                    .map(|row| match row {
                        Self::List(entries) => entries
                            .iter()
                            .map(Self::to_typst_source)
                            .collect::<Vec<_>>()
                            .join(","),
                        value => value.to_typst_source(),
                    })
                    .collect::<Vec<_>>()
                    .join(";");
                format!("mat({rows})")
            }
            "cases" => call("cases"),
            "class" if arguments.len() == 2 => {
                let class = arguments[0].literal_text().unwrap_or("normal");
                format!(
                    r#"class("{}",{})"#,
                    Self::escaped_typst_string(class),
                    sources[1]
                )
            }
            // The constructor rejects unknown nodes. This fixed fallback also
            // keeps forged/imported user data from becoming executable source.
            _ => r#"upright("?")"#.to_owned(),
        }
    }

    fn is_typst_math_name(name: &str) -> bool {
        matches!(
            name,
            "alpha"
                | "beta"
                | "gamma"
                | "delta"
                | "epsilon"
                | "zeta"
                | "eta"
                | "theta"
                | "iota"
                | "kappa"
                | "lambda"
                | "mu"
                | "nu"
                | "xi"
                | "omicron"
                | "pi"
                | "rho"
                | "sigma"
                | "tau"
                | "upsilon"
                | "phi"
                | "chi"
                | "psi"
                | "omega"
        )
    }

    pub fn to_typst_source(&self) -> String {
        match self {
            Self::Symbol(name)
                if name.chars().count() == 1
                    && name
                        .chars()
                        .next()
                        .is_some_and(|character| character.is_alphanumeric()) =>
            {
                name.clone()
            }
            Self::Symbol(name) if Self::is_typst_math_name(name) => name.clone(),
            Self::Symbol(name) => {
                format!(r#"italic("{}")"#, Self::escaped_typst_string(name))
            }
            Self::Text(value) => {
                format!(r#"upright("{}")"#, Self::escaped_typst_string(value))
            }
            Self::Number(number) => number.to_string(),
            Self::Sequence(items) => items
                .iter()
                .map(Self::to_typst_source)
                .collect::<Vec<_>>()
                .join(" "),
            Self::List(items) => items
                .iter()
                .map(Self::to_typst_source)
                .collect::<Vec<_>>()
                .join(","),
            Self::Math { head, arguments } => Self::typst_math_source(head, arguments),
            Self::Attach { base, top, bottom } => {
                let mut source = format!("attach({}", base.to_typst_source());
                if let Some(top) = top {
                    source.push_str(",t:");
                    source.push_str(&top.to_typst_source());
                }
                if let Some(bottom) = bottom {
                    source.push_str(",b:");
                    source.push_str(&bottom.to_typst_source());
                }
                source.push(')');
                source
            }
        }
    }

    pub fn to_native_string(&self) -> String {
        match self {
            Self::Symbol(name) => name.clone(),
            Self::Text(value) => value.clone(),
            Self::Number(number) => number.to_string(),
            Self::Sequence(items) => items
                .iter()
                .map(Self::to_native_string)
                .collect::<Vec<_>>()
                .join(" "),
            Self::List(items) => items
                .iter()
                .map(Self::to_native_string)
                .collect::<Vec<_>>()
                .join(","),
            Self::Math { .. } => self.to_typst_source(),
            Self::Attach { base, top, bottom } => {
                let mut output = base.to_native_string();
                if let Some(top) = top {
                    output.push_str("^(");
                    output.push_str(&top.to_native_string());
                    output.push(')');
                }
                if let Some(bottom) = bottom {
                    output.push_str("_(");
                    output.push_str(&bottom.to_native_string());
                    output.push(')');
                }
                output
            }
        }
    }

    #[cfg(feature = "shadowing")]
    fn node_user_data(&self) -> UserData {
        match self {
            Self::Symbol(name) => UserData::List(vec![
                UserData::String("symbol".to_owned()),
                UserData::String(name.clone()),
            ]),
            Self::Text(value) => UserData::List(vec![
                UserData::String("text".to_owned()),
                UserData::String(value.clone()),
            ]),
            Self::Number(number) => UserData::List(vec![
                UserData::String("number".to_owned()),
                UserData::Integer(*number),
            ]),
            Self::Sequence(items) => UserData::List(vec![
                UserData::String("sequence".to_owned()),
                UserData::List(items.iter().map(Self::node_user_data).collect()),
            ]),
            Self::List(items) => UserData::List(vec![
                UserData::String("list".to_owned()),
                UserData::List(items.iter().map(Self::node_user_data).collect()),
            ]),
            Self::Math { head, arguments } => UserData::List(vec![
                UserData::String("math".to_owned()),
                UserData::String(head.clone()),
                UserData::List(arguments.iter().map(Self::node_user_data).collect()),
            ]),
            Self::Attach { base, top, bottom } => UserData::List(vec![
                UserData::String("attach".to_owned()),
                base.node_user_data(),
                top.as_deref()
                    .map(Self::node_user_data)
                    .unwrap_or(UserData::None),
                bottom
                    .as_deref()
                    .map(Self::node_user_data)
                    .unwrap_or(UserData::None),
            ]),
        }
    }

    #[cfg(feature = "shadowing")]
    fn from_node_user_data(data: &UserData, depth: usize) -> Option<Self> {
        if depth > 16 {
            return None;
        }
        let UserData::List(fields) = data else {
            return None;
        };
        match fields.as_slice() {
            [UserData::String(kind), UserData::String(name)] if kind == "symbol" => {
                Self::symbol(name.clone()).ok()
            }
            [UserData::String(kind), UserData::String(value)] if kind == "text" => {
                Self::text(value.clone()).ok()
            }
            [UserData::String(kind), UserData::Integer(number)] if kind == "number" => {
                Some(Self::Number(*number))
            }
            [UserData::String(kind), UserData::List(items)] if kind == "sequence" => {
                if items.len() > 64 {
                    return None;
                }
                Some(Self::Sequence(
                    items
                        .iter()
                        .map(|item| Self::from_node_user_data(item, depth + 1))
                        .collect::<Option<Vec<_>>>()?,
                ))
            }
            [UserData::String(kind), UserData::List(items)] if kind == "list" => {
                if items.len() > 64 {
                    return None;
                }
                Some(Self::List(
                    items
                        .iter()
                        .map(|item| Self::from_node_user_data(item, depth + 1))
                        .collect::<Option<Vec<_>>>()?,
                ))
            }
            [
                UserData::String(kind),
                UserData::String(head),
                UserData::List(arguments),
            ] if kind == "math" => {
                let arguments = arguments
                    .iter()
                    .map(|argument| Self::from_node_user_data(argument, depth + 1))
                    .collect::<Option<Vec<_>>>()?;
                Self::math(head.clone(), arguments).ok()
            }
            [UserData::String(kind), base, top, bottom] if kind == "attach" => Some(Self::Attach {
                base: Box::new(Self::from_node_user_data(base, depth + 1)?),
                top: match top {
                    UserData::None => None,
                    value => Some(Box::new(Self::from_node_user_data(value, depth + 1)?)),
                },
                bottom: match bottom {
                    UserData::None => None,
                    value => Some(Box::new(Self::from_node_user_data(value, depth + 1)?)),
                },
            }),
            _ => None,
        }
    }

    #[cfg(feature = "shadowing")]
    pub fn symbol_user_data(&self) -> UserData {
        UserData::List(vec![
            UserData::String("spenso::index-display-v1".to_owned()),
            self.node_user_data(),
        ])
    }

    #[cfg(feature = "shadowing")]
    pub fn from_symbol(symbol: Symbol) -> Option<Self> {
        match symbol.get_data() {
            UserData::List(fields) => match fields.as_slice() {
                [UserData::String(kind), display] if kind == "spenso::index-display-v1" => {
                    Self::from_node_user_data(display, 0)
                }
                _ => None,
            },
            _ => None,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum IndexPalette {
    #[default]
    Numeric,
    Cyclic {
        start: usize,
        labels: Box<[IndexDisplay]>,
    },
}

impl IndexPalette {
    pub fn cyclic(
        start: usize,
        labels: impl IntoIterator<Item = IndexDisplay>,
    ) -> Result<Self, RepLibraryError> {
        if i64::try_from(start).is_err() {
            return Err(RepLibraryError::InvalidIndexPalette(
                "the palette start must fit in a signed 64-bit integer".to_owned(),
            ));
        }
        let labels = labels.into_iter().collect::<Vec<_>>().into_boxed_slice();
        if labels.is_empty() {
            return Err(RepLibraryError::InvalidIndexPalette(
                "a cyclic index palette needs at least one label".to_owned(),
            ));
        }
        if labels.len() > 64 {
            return Err(RepLibraryError::InvalidIndexPalette(
                "a cyclic index palette may contain at most 64 labels".to_owned(),
            ));
        }
        Ok(Self::Cyclic { start, labels })
    }

    pub fn resolve(&self, index: usize) -> Option<IndexDisplay> {
        let Self::Cyclic { start, labels } = self else {
            return None;
        };
        let offset = index.checked_sub(*start)?;
        let cycle = offset / labels.len();
        let display = labels[offset % labels.len()].clone();
        if cycle == 0 {
            Some(display)
        } else {
            Some(display.with_bottom(IndexDisplay::Number(i64::try_from(cycle).ok()?)))
        }
    }

    #[cfg(feature = "shadowing")]
    fn to_user_data(&self) -> UserData {
        match self {
            Self::Numeric => UserData::List(vec![UserData::String("numeric".to_owned())]),
            Self::Cyclic { start, labels } => UserData::List(vec![
                UserData::String("cyclic".to_owned()),
                UserData::Integer(*start as i64),
                UserData::List(labels.iter().map(IndexDisplay::node_user_data).collect()),
            ]),
        }
    }

    #[cfg(feature = "shadowing")]
    fn from_user_data(data: &UserData) -> Option<Self> {
        let UserData::List(fields) = data else {
            return None;
        };
        match fields.as_slice() {
            [UserData::String(kind)] if kind == "numeric" => Some(Self::Numeric),
            [
                UserData::String(kind),
                UserData::Integer(start),
                UserData::List(labels),
            ] if kind == "cyclic" && *start >= 0 && !labels.is_empty() && labels.len() <= 64 => {
                Self::cyclic(
                    usize::try_from(*start).ok()?,
                    labels
                        .iter()
                        .map(|label| IndexDisplay::from_node_user_data(label, 0))
                        .collect::<Option<Vec<_>>>()?,
                )
                .ok()
            }
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RepresentationClass {
    SelfDual,
    InlineMetric,
    Dualizable,
}

/// Preferred Typst script row for the base orientation of a representation.
///
/// Dualizable representations use the opposite row for their `dind(...)`
/// orientation. Self-dual representations keep this row in either spelling.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum IndexRow {
    #[default]
    Top,
    Bottom,
}

impl IndexRow {
    pub fn opposite(self) -> Self {
        match self {
            Self::Top => Self::Bottom,
            Self::Bottom => Self::Top,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Top => "top",
            Self::Bottom => "bottom",
        }
    }

    fn from_str(value: &str) -> Option<Self> {
        match value {
            "top" => Some(Self::Top),
            "bottom" => Some(Self::Bottom),
            _ => None,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RepresentationMetadata {
    pub class: RepresentationClass,
    pub label: IndexDisplay,
    pub index_palette: IndexPalette,
    pub index_row: IndexRow,
}

#[cfg(feature = "shadowing")]
impl RepresentationMetadata {
    pub fn to_user_data(&self) -> UserData {
        let class = match self.class {
            RepresentationClass::SelfDual => "self-dual",
            RepresentationClass::InlineMetric => "inline-metric",
            RepresentationClass::Dualizable => "dualizable",
        };
        UserData::List(vec![
            UserData::String("spenso::representation-v2".to_owned()),
            UserData::String(class.to_owned()),
            self.label.node_user_data(),
            self.index_palette.to_user_data(),
            UserData::String(self.index_row.as_str().to_owned()),
        ])
    }

    pub fn from_symbol(symbol: Symbol) -> Option<Self> {
        let UserData::List(fields) = symbol.get_data() else {
            return None;
        };
        let (version, class, label, palette, index_row) = match fields.as_slice() {
            [
                UserData::String(version),
                UserData::String(class),
                label,
                palette,
            ] if version == "spenso::representation-v1" => (
                version,
                class,
                label,
                palette,
                default_index_row(symbol.get_name()),
            ),
            [
                UserData::String(version),
                UserData::String(class),
                label,
                palette,
                UserData::String(index_row),
            ] if version == "spenso::representation-v2" => (
                version,
                class,
                label,
                palette,
                IndexRow::from_str(index_row)?,
            ),
            _ => return None,
        };
        debug_assert!(
            version == "spenso::representation-v1" || version == "spenso::representation-v2"
        );
        Some(Self {
            class: match class.as_str() {
                "self-dual" => RepresentationClass::SelfDual,
                "inline-metric" => RepresentationClass::InlineMetric,
                "dualizable" => RepresentationClass::Dualizable,
                _ => return None,
            },
            label: IndexDisplay::from_node_user_data(label, 0)?,
            index_palette: IndexPalette::from_user_data(palette)?,
            index_row,
        })
    }
}

impl LibraryRep {
    #[cfg(feature = "shadowing")]
    pub fn new_symbol(
        &self,
        name: &str,
        index_palette: IndexPalette,
    ) -> Result<Symbol, RepLibraryError> {
        self.new_symbol_with_index_row(name, index_palette, default_index_row(name))
    }

    #[cfg(feature = "shadowing")]
    fn new_symbol_with_index_row(
        &self,
        name: &str,
        index_palette: IndexPalette,
        index_row: IndexRow,
    ) -> Result<Symbol, RepLibraryError> {
        use symbolica::{atom::AtomCore, printer::PrintState};

        use crate::shadowing::symbolica_utils::SpensoPrintSettings;

        let (class, rep_name, tags) = match self {
            LibraryRep::SelfDual(a) => (
                RepresentationClass::SelfDual,
                encode_base(*a as usize, &LATIN),
                vec![
                    SPENSO_TAG.representation.clone(),
                    SPENSO_TAG.self_dual.clone(),
                ],
            ),
            LibraryRep::Dualizable(a) => (
                RepresentationClass::Dualizable,
                encode_base(a.unsigned_abs() as usize, &CYRILLIC),
                vec![
                    SPENSO_TAG.representation.clone(),
                    SPENSO_TAG.dualizable.clone(),
                ],
            ),
            LibraryRep::InlineMetric(a) => (
                RepresentationClass::InlineMetric,
                encode_base(*a as usize, &GREEK),
                vec![
                    SPENSO_TAG.representation.clone(),
                    SPENSO_TAG.self_dual.clone(),
                ],
            ),
            LibraryRep::Dummy => (
                RepresentationClass::SelfDual,
                String::new(),
                vec![
                    SPENSO_TAG.representation.clone(),
                    SPENSO_TAG.self_dual.clone(),
                ],
            ),
        };
        let label_name = representation_label_name(name);
        let metadata = RepresentationMetadata {
            class,
            label: IndexDisplay::symbol(label_name)?,
            index_palette,
            index_row,
        };
        let body = representation_typst_body(&metadata.label);
        let qualified_name = canonical_representation_name(name);
        let namespaced = NamespacedSymbol::parse(&qualified_name);

        if let Some(existing) = Symbol::get_symbol(namespaced.clone()) {
            let correct_tags = tags.iter().all(|tag| existing.has_tag(tag))
                && match class {
                    RepresentationClass::SelfDual | RepresentationClass::InlineMetric => {
                        !existing.has_tag(&SPENSO_TAG.dualizable)
                    }
                    RepresentationClass::Dualizable => !existing.has_tag(&SPENSO_TAG.self_dual),
                };
            if !correct_tags || RepresentationMetadata::from_symbol(existing) != Some(metadata) {
                return Err(RepLibraryError::AlreadyExistsDifferentMetadata(
                    name.to_owned(),
                ));
            }
            return Ok(existing);
        }

        let print_name = name.to_owned();
        SymbolBuilder::new(namespaced)
            .with_tags(tags)
            .with_user_data(metadata.to_user_data())
            .with_print_function(move |a, opt, _state| {
                if matches!(
                    opt.custom_print_mode.get("typst"),
                    Some(PrintUserData::Integer(1))
                ) {
                    return Some(body.clone());
                }

                if let Some(resolved) = SpensoPrintSettings::resolve(opt) {
                    let (script_open, script_close) = resolved.script_delimiters();
                    let SpensoPrintSettings {
                        with_dim,
                        commas,
                        parens,
                        index_subscripts,
                        ..
                    } = resolved.presentation;
                    let AtomView::Fun(f) = a else {
                        return None;
                    };

                    let mut out = if opt.color_builtin_symbols {
                        nu_ansi_term::Color::DarkGray.paint(&rep_name).to_string()
                    } else {
                        rep_name.clone()
                    };

                    let mut arg_iter = f.iter();
                    let dim = arg_iter.next()?;
                    if f.get_nargs() == 1 {
                        if with_dim {
                            out.push('(');
                            dim.format(&mut out, opt, PrintState::new()).ok()?;
                            out.push(')');
                        }
                        return Some(out);
                    }

                    if f.get_nargs() == 2 {
                        if index_subscripts {
                            out.push('_');
                        }
                        if parens && index_subscripts {
                            out.push(script_open);
                        }
                        if with_dim {
                            dim.format(&mut out, opt, PrintState::new()).ok()?;
                            if commas {
                                out.push(',');
                            } else {
                                out.push(' ');
                            }
                        }
                        let ind = arg_iter.next()?;
                        let palette_index = if let AtomView::Num(number) = ind {
                            match number.get_coeff_view() {
                                CoefficientView::Natural(index, 1, 0, 1) => {
                                    usize::try_from(index).ok()
                                }
                                _ => None,
                            }
                        } else {
                            None
                        };
                        let palette = RepresentationMetadata::from_symbol(f.get_symbol())
                            .map(|metadata| metadata.index_palette);
                        if let Some(display) =
                            palette_index.and_then(|index| palette.as_ref()?.resolve(index))
                        {
                            out.push_str(&display.to_native_string());
                        } else {
                            ind.format(&mut out, opt, PrintState::new()).ok()?;
                        }
                        if parens && index_subscripts {
                            out.push(script_close);
                        }

                        return Some(out);
                    }

                    return None;
                }

                let AtomView::Fun(f) = a else {
                    return None;
                };

                let mut out = if opt.color_builtin_symbols {
                    nu_ansi_term::Color::DarkGray.paint(&print_name).to_string()
                } else {
                    return None;
                };

                out.push('(');
                let mut first = true;
                for arg in f.iter() {
                    if !first {
                        out.push_str(", ");
                    } else {
                        first = false;
                    }
                    out.push_str(&arg.to_string());
                }
                out.push(')');
                Some(out)
            })
            .build()
            .map_err(|error| RepLibraryError::SymbolRegistration {
                name: name.to_owned(),
                reason: error.to_string(),
            })
    }

    pub fn new_dual(name: &str) -> Result<Self, RepLibraryError> {
        REPS.write().unwrap().new_dual_impl(name)
    }

    #[cfg(feature = "shadowing")]
    pub fn new_dual_with_index_palette(
        name: &str,
        index_palette: IndexPalette,
    ) -> Result<Self, RepLibraryError> {
        REPS.write()
            .unwrap()
            .new_dual_impl_with_index_palette(name, index_palette)
    }

    #[cfg(feature = "shadowing")]
    pub fn new_dual_with_index_palette_and_row(
        name: &str,
        index_palette: IndexPalette,
        index_row: IndexRow,
    ) -> Result<Self, RepLibraryError> {
        REPS.write()
            .unwrap()
            .new_dual_impl_with_index_palette_and_row(name, index_palette, index_row)
    }

    #[cfg(feature = "shadowing")]
    pub fn symbol(&self) -> Symbol {
        REPS.read().unwrap()[*self].symbol
    }

    pub fn name(&self) -> String {
        REPS.read().unwrap()[*self].name.clone()
    }

    #[cfg(feature = "shadowing")]
    pub fn metadata(&self) -> Option<RepresentationMetadata> {
        RepresentationMetadata::from_symbol(self.symbol())
    }

    #[cfg(feature = "shadowing")]
    fn from_registered_or_portable_symbol(symbol: Symbol) -> Result<Self, RepresentationError> {
        let registered = { REPS.read().unwrap().find_symbol(symbol) };
        if let Some(representation) = registered {
            return Ok(representation);
        }

        if !symbol.has_tag(&SPENSO_TAG.representation) {
            return Err(RepresentationError::NotRepresentationError(symbol));
        }
        let metadata = RepresentationMetadata::from_symbol(symbol)
            .ok_or(RepresentationError::NotRepresentationError(symbol))?;
        let self_dual = symbol.has_tag(&SPENSO_TAG.self_dual);
        let dualizable = symbol.has_tag(&SPENSO_TAG.dualizable);
        let name = symbol.get_name();

        match metadata.class {
            RepresentationClass::SelfDual if self_dual && !dualizable => {
                Self::new_self_dual_with_index_palette_and_row(
                    name,
                    metadata.index_palette,
                    metadata.index_row,
                )
                .map_err(Into::into)
            }
            RepresentationClass::Dualizable if dualizable && !self_dual => {
                Self::new_dual_with_index_palette_and_row(
                    name,
                    metadata.index_palette,
                    metadata.index_row,
                )
                .map_err(Into::into)
            }
            RepresentationClass::InlineMetric if self_dual && !dualizable => {
                Err(RepresentationError::ImportedInlineMetricRequiresLocalRegistration(symbol))
            }
            _ => Err(RepresentationError::NotRepresentationError(symbol)),
        }
    }

    pub fn new_self_dual(name: &str) -> Result<Self, RepLibraryError> {
        REPS.write().unwrap().new_self_dual(name)
    }

    #[cfg(feature = "shadowing")]
    pub fn new_self_dual_with_index_palette(
        name: &str,
        index_palette: IndexPalette,
    ) -> Result<Self, RepLibraryError> {
        REPS.write()
            .unwrap()
            .new_self_dual_with_index_palette(name, index_palette)
    }

    #[cfg(feature = "shadowing")]
    pub fn new_self_dual_with_index_palette_and_row(
        name: &str,
        index_palette: IndexPalette,
        index_row: IndexRow,
    ) -> Result<Self, RepLibraryError> {
        REPS.write()
            .unwrap()
            .new_self_dual_with_index_palette_and_row(name, index_palette, index_row)
    }

    pub fn all_self_duals() -> impl Iterator<Item = &'static LibraryRep> {
        SELF_DUAL.iter().map(|(rep, _)| rep)
    }

    pub fn all_dualizables() -> impl Iterator<Item = &'static LibraryRep> {
        DUALIZABLE.iter().map(|(rep, _)| rep)
    }

    pub fn all_inline_metrics() -> impl Iterator<Item = &'static LibraryRep> {
        INLINE_METRIC.iter().map(|(rep, _)| rep)
    }

    pub fn all_representations() -> impl Iterator<Item = &'static LibraryRep> {
        Self::all_self_duals()
            .chain(Self::all_dualizables())
            .chain(Self::all_inline_metrics())
    }
}

pub struct MetricRepData {
    metric_data: fn(ConcreteIndex) -> bool,
    rep_data: RepData,
}

pub struct RepData {
    // metric_data: Fn(Dimension)->SparseTensor<i8,IndexLess>
    name: String,
    #[cfg(feature = "shadowing")]
    symbol: Symbol,
}

static DUMMY_REP_DATA: LazyLock<RepData> = LazyLock::new(|| RepData {
    name: "Dummy".to_string(),
    #[cfg(feature = "shadowing")]
    symbol: self_dual_symbol!("Dummy"),
});

pub struct ExtendibleReps {
    name_map: AHashMap<String, LibraryRep>,
    #[cfg(feature = "shadowing")]
    symbol_map: AHashMap<Symbol, LibraryRep>,
}

#[derive(Debug, Error)]
pub enum RepLibraryError {
    #[error("{0} Already exists and is of different type")]
    AlreadyExistsDifferentType(String),
    #[error("{0} Already exists and has different metric function")]
    AlreadyExistsDifferentMetric(String),
    #[error("invalid index display: {0}")]
    InvalidIndexDisplay(String),
    #[error("invalid index palette: {0}")]
    InvalidIndexPalette(String),
    #[error("{0} already exists with different representation metadata")]
    AlreadyExistsDifferentMetadata(String),
    #[error("could not register representation symbol {name}: {reason}")]
    SymbolRegistration { name: String, reason: String },
}

impl ExtendibleReps {
    pub fn reps(&self) -> impl Iterator<Item = &LibraryRep> {
        self.name_map.values()
    }

    pub fn new_dual_impl(&mut self, name: &str) -> Result<LibraryRep, RepLibraryError> {
        self.new_dual_impl_with_metadata_request(name, None, None)
    }

    #[cfg(feature = "shadowing")]
    pub fn new_dual_impl_with_index_palette(
        &mut self,
        name: &str,
        index_palette: IndexPalette,
    ) -> Result<LibraryRep, RepLibraryError> {
        self.new_dual_impl_with_metadata_request(name, Some(index_palette), None)
    }

    #[cfg(feature = "shadowing")]
    pub fn new_dual_impl_with_index_palette_and_row(
        &mut self,
        name: &str,
        index_palette: IndexPalette,
        index_row: IndexRow,
    ) -> Result<LibraryRep, RepLibraryError> {
        self.new_dual_impl_with_metadata_request(name, Some(index_palette), Some(index_row))
    }

    #[cfg(feature = "shadowing")]
    fn new_dual_impl_with_metadata_request(
        &mut self,
        name: &str,
        index_palette: Option<IndexPalette>,
        index_row: Option<IndexRow>,
    ) -> Result<LibraryRep, RepLibraryError> {
        let canonical_name = canonical_representation_name(name);
        if let Some(rep) = self.name_map.get(&canonical_name) {
            if !matches!(rep, LibraryRep::Dualizable(_)) {
                return Err(RepLibraryError::AlreadyExistsDifferentType(name.into()));
            }
            if let Some(index_palette) = &index_palette
                && RepresentationMetadata::from_symbol(self[*rep].symbol)
                    .is_none_or(|metadata| metadata.index_palette != *index_palette)
            {
                return Err(RepLibraryError::AlreadyExistsDifferentMetadata(name.into()));
            }
            if let Some(index_row) = index_row
                && RepresentationMetadata::from_symbol(self[*rep].symbol)
                    .is_none_or(|metadata| metadata.index_row != index_row)
            {
                return Err(RepLibraryError::AlreadyExistsDifferentMetadata(name.into()));
            }
            return Ok(*rep);
        }

        let rep = LibraryRep::Dualizable(DUALIZABLE.len() as i16 + 1);
        let symbol = rep.new_symbol_with_index_row(
            &canonical_name,
            index_palette.unwrap_or_else(|| default_index_palette(&canonical_name)),
            index_row.unwrap_or_else(|| default_index_row(&canonical_name)),
        )?;
        self.name_map.insert(canonical_name, rep);
        self.symbol_map.insert(symbol, rep);

        DUALIZABLE.push((
            rep,
            RepData {
                name: name.to_string(),
                #[cfg(feature = "shadowing")]
                symbol,
            },
        ));
        Ok(rep)
    }

    #[cfg(not(feature = "shadowing"))]
    fn new_dual_impl_with_metadata_request(
        &mut self,
        name: &str,
        _index_palette: Option<IndexPalette>,
        _index_row: Option<IndexRow>,
    ) -> Result<LibraryRep, RepLibraryError> {
        let canonical_name = canonical_representation_name(name);
        if let Some(rep) = self.name_map.get(&canonical_name) {
            if !matches!(rep, LibraryRep::Dualizable(_)) {
                return Err(RepLibraryError::AlreadyExistsDifferentType(name.into()));
            }
            return Ok(*rep);
        }

        let rep = LibraryRep::Dualizable(DUALIZABLE.len() as i16 + 1);
        self.name_map.insert(canonical_name, rep);
        DUALIZABLE.push((
            rep,
            RepData {
                name: name.to_string(),
            },
        ));
        Ok(rep)
    }

    pub fn new_dual(name: &str) -> Result<LibraryRep, RepLibraryError> {
        REPS.write().unwrap().new_dual_impl(name)
    }

    pub fn new_self_dual(&mut self, name: &str) -> Result<LibraryRep, RepLibraryError> {
        self.new_self_dual_impl(name, None, None)
    }

    #[cfg(feature = "shadowing")]
    pub fn new_self_dual_with_index_palette(
        &mut self,
        name: &str,
        index_palette: IndexPalette,
    ) -> Result<LibraryRep, RepLibraryError> {
        self.new_self_dual_impl(name, Some(index_palette), None)
    }

    #[cfg(feature = "shadowing")]
    pub fn new_self_dual_with_index_palette_and_row(
        &mut self,
        name: &str,
        index_palette: IndexPalette,
        index_row: IndexRow,
    ) -> Result<LibraryRep, RepLibraryError> {
        self.new_self_dual_impl(name, Some(index_palette), Some(index_row))
    }

    fn new_self_dual_impl(
        &mut self,
        name: &str,
        index_palette: Option<IndexPalette>,
        index_row: Option<IndexRow>,
    ) -> Result<LibraryRep, RepLibraryError> {
        let canonical_name = canonical_representation_name(name);
        if let Some(rep) = self.name_map.get(&canonical_name) {
            if let LibraryRep::Dualizable(_) = rep {
                return Err(RepLibraryError::AlreadyExistsDifferentType(name.into()));
            }
            #[cfg(feature = "shadowing")]
            if let Some(index_palette) = &index_palette
                && RepresentationMetadata::from_symbol(self[*rep].symbol)
                    .is_none_or(|metadata| metadata.index_palette != *index_palette)
            {
                return Err(RepLibraryError::AlreadyExistsDifferentMetadata(name.into()));
            }
            #[cfg(feature = "shadowing")]
            if let Some(index_row) = index_row
                && RepresentationMetadata::from_symbol(self[*rep].symbol)
                    .is_none_or(|metadata| metadata.index_row != index_row)
            {
                return Err(RepLibraryError::AlreadyExistsDifferentMetadata(name.into()));
            }
            return Ok(*rep);
        }

        let rep = LibraryRep::SelfDual(SELF_DUAL.len() as u16);

        #[cfg(feature = "shadowing")]
        let symbol = rep.new_symbol_with_index_row(
            &canonical_name,
            index_palette.unwrap_or_else(|| default_index_palette(&canonical_name)),
            index_row.unwrap_or_else(|| default_index_row(&canonical_name)),
        )?;
        #[cfg(not(feature = "shadowing"))]
        let _ = (index_palette, index_row);
        self.name_map.insert(canonical_name, rep);
        #[cfg(feature = "shadowing")]
        self.symbol_map.insert(symbol, rep);

        SELF_DUAL.push((
            rep,
            RepData {
                name: name.to_string(),
                #[cfg(feature = "shadowing")]
                symbol,
            },
        ));
        Ok(rep)
    }

    #[allow(unpredictable_function_pointer_comparisons)]
    pub fn new_inline_metric(
        &mut self,
        name: &str,
        metric_fn: fn(ConcreteIndex) -> bool,
    ) -> Result<LibraryRep, RepLibraryError> {
        let canonical_name = canonical_representation_name(name);
        if let Some(rep) = self.name_map.get(&canonical_name) {
            match rep {
                LibraryRep::SelfDual(_) | LibraryRep::Dualizable(_) | LibraryRep::Dummy => {
                    return Err(RepLibraryError::AlreadyExistsDifferentType(name.into()));
                }
                LibraryRep::InlineMetric(a) => {
                    if INLINE_METRIC[*a as usize].1.metric_data == metric_fn {
                        #[cfg(feature = "shadowing")]
                        if RepresentationMetadata::from_symbol(self[*rep].symbol).is_none_or(
                            |metadata| {
                                metadata.index_palette != default_index_palette(&canonical_name)
                            },
                        ) {
                            return Err(RepLibraryError::AlreadyExistsDifferentMetadata(
                                name.into(),
                            ));
                        }
                        return Ok(*rep);
                    } else {
                        return Err(RepLibraryError::AlreadyExistsDifferentMetric(
                            name.to_string(),
                        ));
                    }
                }
            }
        }

        let rep = LibraryRep::InlineMetric(INLINE_METRIC.len() as u16);
        #[cfg(feature = "shadowing")]
        let symbol = rep.new_symbol(&canonical_name, default_index_palette(&canonical_name))?;
        self.name_map.insert(canonical_name, rep);
        #[cfg(feature = "shadowing")]
        self.symbol_map.insert(symbol, rep);

        INLINE_METRIC.push((
            rep,
            MetricRepData {
                metric_data: metric_fn,
                rep_data: RepData {
                    name: name.to_string(),
                    #[cfg(feature = "shadowing")]
                    symbol,
                },
            },
        ));
        Ok(rep)
    }
}

impl Index<LibraryRep> for ExtendibleReps {
    type Output = RepData;

    fn index(&self, index: LibraryRep) -> &Self::Output {
        match index {
            LibraryRep::Dummy => &DUMMY_REP_DATA,
            LibraryRep::SelfDual(l) => &SELF_DUAL[l as usize].1,
            LibraryRep::InlineMetric(l) => &INLINE_METRIC[l as usize].1.rep_data,
            LibraryRep::Dualizable(l) => &DUALIZABLE[l.unsigned_abs() as usize - 1].1,
        }
    }
}

impl ExtendibleReps {
    pub const EUCLIDEAN: LibraryRep = LibraryRep::SelfDual(0);
    pub const MINKOWSKI: LibraryRep = LibraryRep::InlineMetric(0);
    pub const LORENTZ_UP: LibraryRep = LibraryRep::Dualizable(1);
    pub const LORENTZ_DOWN: LibraryRep = LibraryRep::Dualizable(-1);

    pub fn new() -> Self {
        let mut new = Self {
            name_map: AHashMap::new(),
            #[cfg(feature = "shadowing")]
            symbol_map: AHashMap::new(),
        };

        // #[cfg(feature = "shadowing")]
        // let _ = AIND_SYMBOLS.aind;
        new.new_self_dual(Euclidean::NAME).unwrap();
        fn mink_is_neg(id: ConcreteIndex) -> bool {
            Minkowski {}.is_neg(id)
        }

        new.new_inline_metric(Minkowski::NAME, mink_is_neg).unwrap();
        new.new_dual_impl(Lorentz::NAME).unwrap();

        new
    }

    #[cfg(feature = "shadowing")]
    pub fn find_symbol(&self, symbol: Symbol) -> Option<LibraryRep> {
        self.symbol_map.get(&symbol).cloned()
    }
}

impl Default for ExtendibleReps {
    fn default() -> Self {
        Self::new()
    }
}

impl Display for LibraryRep {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Dummy => write!(f, "Dummy"),
            Self::SelfDual(_) => write!(f, "{}", REPS.read().unwrap()[*self].name),
            Self::InlineMetric(_) => write!(f, "{}", REPS.read().unwrap()[*self].name),
            Self::Dualizable(l) => {
                if *l < 0 {
                    write!(f, "{}🠓", REPS.read().unwrap()[*self].name)
                } else {
                    write!(f, "{}🠑", REPS.read().unwrap()[*self].name)
                }
            }
        }
    }
}

pub fn initialize() {
    let _ = LibraryRep::from(Minkowski {}).to_string();
    #[cfg(feature = "shadowing")]
    let _ = ETS.metric;
}

impl RepName for LibraryRep {
    type Dual = LibraryRep;
    type Base = LibraryRep;

    fn from_library_rep(rep: LibraryRep) -> Result<Self, RepresentationError> {
        Ok(rep)
    }

    fn orientation(self) -> Orientation {
        match self {
            Self::Dummy => Orientation::Undirected,
            Self::SelfDual(_) => Orientation::Undirected,
            Self::InlineMetric(_) => Orientation::Undirected,
            Self::Dualizable(l) => match l.cmp(&0) {
                Ordering::Greater => Orientation::Default,
                Ordering::Less => Orientation::Reversed,
                Ordering::Equal => panic!("dualizable with 0"),
            },
        }
    }

    #[inline]
    fn dual(self) -> Self::Dual {
        match self {
            Self::Dummy => Self::Dummy,
            Self::SelfDual(l) => Self::SelfDual(l),
            Self::InlineMetric(l) => Self::InlineMetric(l),
            Self::Dualizable(l) => Self::Dualizable(-l),
        }
    }

    #[inline]
    fn is_base(&self) -> bool {
        match self {
            Self::Dualizable(l) => *l > 0,
            _ => true,
        }
    }

    #[inline]
    fn is_dual(self) -> bool {
        match self {
            Self::Dualizable(l) => l < 0,
            _ => true,
        }
    }

    fn is_self_dual(&self) -> bool {
        !matches!(self, Self::Dualizable(_))
    }

    #[inline]
    fn base(&self) -> Self::Base {
        match self {
            Self::Dualizable(l) => Self::Dualizable(l.abs()),
            x => *x,
        }
    }

    #[inline]
    fn matches(&self, other: &Self::Dual) -> bool {
        match (self, other) {
            (Self::SelfDual(s), Self::SelfDual(o)) => s == o,
            (Self::Dualizable(s), Self::Dualizable(o)) => *s == -o,
            (Self::InlineMetric(s), Self::InlineMetric(o)) => s == o,
            _ => false,
        }
    }

    fn match_cmp(&self, other: &Self::Dual) -> Ordering {
        match (self, other) {
            (Self::SelfDual(s), Self::SelfDual(o))
            | (Self::InlineMetric(s), Self::InlineMetric(o)) => s.cmp(o),
            (Self::Dualizable(s), Self::Dualizable(o)) => s.abs().cmp(&o.abs()),
            _ => self.cmp(other),
        }
    }

    #[cfg(feature = "shadowing")]
    fn try_from_symbol(sym: Symbol, aind: Symbol) -> Result<Self, RepresentationError> {
        use super::abstract_index::AIND_SYMBOLS;

        let rep = Self::from_registered_or_portable_symbol(sym)?;

        match rep {
            LibraryRep::Dualizable(_) => {
                if aind == AIND_SYMBOLS.dind {
                    Ok(rep.dual())
                } else if aind == AIND_SYMBOLS.uind {
                    Ok(rep)
                } else if aind == AIND_SYMBOLS.selfdualind {
                    Err(RepresentationError::ExpectedDualStateError(
                        AIND_SYMBOLS.uind,
                        aind,
                    ))
                } else {
                    Err(RepresentationError::SymbolError(aind))
                }
            }
            LibraryRep::SelfDual(_) | LibraryRep::InlineMetric(_) | LibraryRep::Dummy => {
                if aind == AIND_SYMBOLS.selfdualind {
                    Ok(rep)
                } else if aind == AIND_SYMBOLS.dind || aind == AIND_SYMBOLS.uind {
                    Err(RepresentationError::ExpectedDualStateError(
                        AIND_SYMBOLS.selfdualind,
                        aind,
                    ))
                } else {
                    Err(RepresentationError::SymbolError(aind))
                }
            }
        }
    }

    #[cfg(feature = "shadowing")]
    fn try_from_symbol_coerced(sym: Symbol) -> Result<Self, RepresentationError> {
        Self::from_registered_or_portable_symbol(sym)
    }

    fn is_neg(self, i: usize) -> bool {
        if let LibraryRep::InlineMetric(a) = self {
            (INLINE_METRIC[a as usize].1.metric_data)(i)
        } else {
            false
        }
    }

    fn is_dummy(self) -> bool {
        matches!(self, LibraryRep::Dummy)
    }

    #[cfg(feature = "shadowing")]
    /// yields a function builder for the representation
    fn to_symbolic<'a, It: Into<AtomOrView<'a>>>(
        &self,
        args: impl IntoIterator<Item = It>,
    ) -> Atom {
        use crate::structure::abstract_index::AIND_SYMBOLS;

        let mut fun = FunctionBuilder::new(self.symbol());
        for a in args {
            fun = fun.add_arg(a);
        }
        let inner = fun.finish();

        match self {
            Self::SelfDual(_) | Self::Dummy | Self::InlineMetric(_) => inner,
            Self::Dualizable(l) => {
                if *l < 0 {
                    function!(AIND_SYMBOLS.dind, &inner)
                } else {
                    inner
                }
            }
        }
    }
}

#[test]
fn extendible_reps() {
    let r = LibraryRep::new_dual("lor").unwrap();
    let rd = r.dual();
    let e = LibraryRep::new_self_dual("euc").unwrap();

    println!(
        "{r}{r:?}, {e}{e:?},{rd}{rd:?}",
        // ExtendibleReps::BISPINOR.base()
    );

    assert!(ExtendibleReps::LORENTZ_UP.matches(&ExtendibleReps::LORENTZ_DOWN));
    assert!(!ExtendibleReps::LORENTZ_UP.matches(&ExtendibleReps::LORENTZ_UP));
    // assert!(ExtendibleReps::BISPINOR.matches(&ExtendibleReps::BISPINOR));

    // let rs = r.new_slot(10, 1);
    // let rr = r.new_dimed_rep(1);

    // // println!("{}", rs.to_symbolic());
    // println!("{}", rs.dual());
    // println!("{}", rr)
}

// struct UserDefRep{
//     dual: usize,
//     name: String,
// }

// struct UserDefReps

// }

// impl RepName for usize{
//     type
// }

// pub trait HasDimension: RepName {
//     fn dim(&self) -> Dimension;

//     fn to_fnbuilder(&self) -> FunctionBuilder {
//         ::to_fnbuilder().add_arg(self.dim().to_symbolic().as_atom_view())
//     }
// }

// impl<T: BaseRepName<Base: BaseRepName, Dual: BaseRepName>> Representation<T> {
//     pub fn dual_pair(self) -> Representation<DualPair<T::Base>>
//     where
//         <T::Base as RepName>::Dual: RepName<Dual = T::Base, Base = T::Base>,
//         T::Base: RepName<Dual = T::Dual, Base = T::Base>,
//     {
//         Representation {
//             dim: self.dim,
//             rep: T::selfless_pair(),
//         }
//     }
// }

impl<T: RepName> Representation<T> {
    pub fn matches(&self, other: &Representation<T::Dual>) -> bool {
        self.dim == other.dim && self.rep.matches(&other.rep)
    }

    pub fn match_cmp(&self, other: &Representation<T::Dual>) -> Ordering {
        self.dim
            .cmp(&other.dim)
            .then(self.rep.match_cmp(&other.rep))
    }

    #[cfg(feature = "shadowing")]
    /// yields a function builder for the representation, adding a first variable: the dimension.
    ///
    pub fn to_symbolic(&self, args: impl IntoIterator<Item = Atom>) -> Atom {
        self.rep
            .to_symbolic([self.dim.to_symbolic()].into_iter().chain(args))
    }
    pub fn dual(self) -> Representation<T::Dual> {
        Representation {
            dim: self.dim,
            rep: self.rep.dual(),
        }
    }

    pub fn cast<U: RepName + From<T>>(self) -> Representation<U> {
        Representation {
            dim: self.dim,
            rep: U::from(self.rep),
        }
    }

    pub fn is_neg(&self, i: usize) -> bool {
        self.rep.is_neg(i)
    }

    pub fn slot<Aind, A: Into<Aind>>(&self, aind: A) -> Slot<T, Aind> {
        Slot {
            aind: aind.into(),
            rep: *self,
        }
    }

    #[inline]
    // this could be implemented directly in the fiberiterator.
    /// gives the vector of booleans, saying which concrete index along a Dimension/Abstract Index should have a minus sign during contraction.
    ///
    pub fn negative(&self) -> Result<Vec<bool>> {
        Ok((0..usize::try_from(self.dim)?)
            .map(|i| self.is_neg(i))
            .collect())
    }
}

#[test]
fn test_negative() {
    let spin: Representation<Euclidean> = Euclidean {}.new_rep(5);

    let metric_diag: Vec<bool> = spin.negative().unwrap();

    let mut agree = true;

    for (i, r) in metric_diag.iter().enumerate() {
        if r ^ spin.is_neg(i) {
            agree = false;
        }
    }

    assert!(agree);
}

#[cfg(feature = "shadowing")]
/// Can possibly constuct a Representation from an `AtomView`, if it is of the form: <representation>(<dimension>)
///
impl<'a, T: RepName> TryFrom<AtomView<'a>> for Representation<T> {
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

        Ok(Representation { dim, rep })
    }
}

impl<T: RepName> Display for Representation<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.rep, self.dim)
    }
}

impl From<Dimension> for Representation<Euclidean> {
    fn from(value: Dimension) -> Self {
        Representation {
            dim: value,
            rep: Euclidean {},
        }
    }
}

impl<T: RepName> From<Representation<T>> for Dimension {
    fn from(value: Representation<T>) -> Self {
        value.dim
    }
}

impl<T: RepName> TryFrom<Representation<T>> for usize {
    type Error = DimensionError;
    fn try_from(value: Representation<T>) -> std::result::Result<Self, Self::Error> {
        usize::try_from(value.dim)
    }
}

impl<'a, T: RepName> FromIterator<&'a Representation<T>> for Vec<Dimension> {
    fn from_iter<I: IntoIterator<Item = &'a Representation<T>>>(iter: I) -> Self {
        iter.into_iter().map(|rep| rep.dim).collect()
    }
}

#[cfg(test)]
mod test {
    use linnet::permutation::Permutation;

    use crate::structure::representation::{RepName, Representation};

    use super::{Euclidean, LibraryRep, Lorentz, Minkowski};

    #[test]
    fn ordering() {
        let reps: Vec<LibraryRep> = vec![
            Euclidean {}.into(),
            Minkowski {}.into(),
            Lorentz {}.into(),
            Lorentz {}.dual().into(),
        ];

        assert!(Permutation::sort(&reps).is_identity());
        let reps: Vec<Representation<LibraryRep>> = vec![
            Euclidean {}.new_rep(3).cast(),
            Minkowski {}.new_rep(2).cast(),
            Lorentz {}.new_rep(1).cast(),
        ];

        assert!(Permutation::sort(&reps).is_identity())
    }
}

#[cfg(test)]
#[cfg(feature = "shadowing")]
mod shadowing_tests {
    use symbolica::atom::{Atom, AtomView, NamespacedSymbol, Symbol, SymbolBuilder};

    use crate::network::tags::SPENSO_TAG;

    use super::{
        IndexDisplay, IndexPalette, IndexRow, LibraryRep, Minkowski, RepLibraryError, RepName,
        RepresentationClass, RepresentationError, RepresentationMetadata,
        representation_typst_body,
    };

    fn greek_palette() -> IndexPalette {
        IndexPalette::cyclic(
            1,
            [
                IndexDisplay::symbol("mu").unwrap(),
                IndexDisplay::symbol("nu").unwrap(),
            ],
        )
        .unwrap()
    }

    fn portable_representation_symbol(
        name: &str,
        class: RepresentationClass,
        palette: IndexPalette,
    ) -> Symbol {
        let label = IndexDisplay::symbol(name.rsplit("::").next().unwrap()).unwrap();
        let tags = match class {
            RepresentationClass::SelfDual | RepresentationClass::InlineMetric => vec![
                SPENSO_TAG.representation.clone(),
                SPENSO_TAG.self_dual.clone(),
            ],
            RepresentationClass::Dualizable => vec![
                SPENSO_TAG.representation.clone(),
                SPENSO_TAG.dualizable.clone(),
            ],
        };
        SymbolBuilder::new(NamespacedSymbol::parse(name))
            .with_tags(tags)
            .with_user_data(
                RepresentationMetadata {
                    class,
                    label,
                    index_palette: palette,
                    index_row: IndexRow::Top,
                }
                .to_user_data(),
            )
            .build()
            .unwrap()
    }

    #[test]
    fn cyclic_palette_wraps_with_one_based_subscripts() {
        let palette = greek_palette();

        assert_eq!(palette.resolve(0), None);
        assert_eq!(
            palette.resolve(1),
            Some(IndexDisplay::symbol("mu").unwrap())
        );
        assert_eq!(
            palette.resolve(2),
            Some(IndexDisplay::symbol("nu").unwrap())
        );
        assert_eq!(
            palette.resolve(3),
            Some(
                IndexDisplay::symbol("mu")
                    .unwrap()
                    .with_bottom(IndexDisplay::Number(1))
            )
        );
        assert_eq!(
            palette.resolve(5),
            Some(
                IndexDisplay::symbol("mu")
                    .unwrap()
                    .with_bottom(IndexDisplay::Number(2))
            )
        );
    }

    #[test]
    fn canonical_representations_own_conventional_index_palettes() {
        let cases = [
            (
                LibraryRep::from(Minkowski {}),
                RepresentationClass::InlineMetric,
                &["mu", "nu", "rho", "sigma"][..],
                IndexRow::Top,
            ),
            (
                LibraryRep::new_self_dual("euc").unwrap(),
                RepresentationClass::SelfDual,
                &["i", "j", "k", "l"][..],
                IndexRow::Top,
            ),
            (
                LibraryRep::new_dual("lor").unwrap(),
                RepresentationClass::Dualizable,
                &["mu", "nu", "rho", "sigma"][..],
                IndexRow::Top,
            ),
            (
                LibraryRep::new_self_dual("bis").unwrap(),
                RepresentationClass::SelfDual,
                &["a", "b", "c", "d"][..],
                IndexRow::Bottom,
            ),
            (
                LibraryRep::new_dual("spf").unwrap(),
                RepresentationClass::Dualizable,
                &["alpha", "beta", "gamma", "delta"][..],
                IndexRow::Top,
            ),
            (
                LibraryRep::new_dual("cof").unwrap(),
                RepresentationClass::Dualizable,
                &["i", "j", "k", "l"][..],
                IndexRow::Top,
            ),
            (
                LibraryRep::new_self_dual("coad").unwrap(),
                RepresentationClass::SelfDual,
                &["a", "b", "c", "d"][..],
                IndexRow::Top,
            ),
            (
                LibraryRep::new_dual("cos").unwrap(),
                RepresentationClass::Dualizable,
                &["I", "J", "K", "L"][..],
                IndexRow::Top,
            ),
        ];

        for (representation, class, labels, row) in cases {
            let metadata = representation.metadata().unwrap();
            assert_eq!(metadata.class, class);
            assert_eq!(metadata.index_row, row);
            let expected = IndexPalette::cyclic(
                1,
                labels
                    .iter()
                    .map(|label| IndexDisplay::symbol(*label).unwrap()),
            )
            .unwrap();
            assert_eq!(metadata.index_palette, expected);
            assert_eq!(
                metadata.index_palette.resolve(5),
                Some(
                    IndexDisplay::symbol(labels[0])
                        .unwrap()
                        .with_bottom(IndexDisplay::Number(1))
                )
            );
        }
    }

    #[test]
    fn representation_symbol_is_the_palette_authority() {
        let name = "spenso_rep_metadata_tests::M";
        let palette = greek_palette();
        let representation =
            LibraryRep::new_self_dual_with_index_palette(name, palette.clone()).unwrap();
        let repeated = LibraryRep::new_self_dual_with_index_palette(name, palette.clone()).unwrap();

        assert_eq!(representation, repeated);
        let metadata = representation.metadata().unwrap();
        assert_eq!(metadata.label, IndexDisplay::symbol("M").unwrap());
        assert_eq!(metadata.index_palette, palette);

        let conflicting = IndexPalette::cyclic(1, [IndexDisplay::symbol("rho").unwrap()]).unwrap();
        assert!(matches!(
            LibraryRep::new_self_dual_with_index_palette(name, conflicting),
            Err(RepLibraryError::AlreadyExistsDifferentMetadata(_))
        ));
    }

    #[test]
    fn qualified_and_unqualified_names_share_one_registry_entry() {
        let name = "CanonicalPaletteAliasForSpensoTests";
        let qualified_name = format!("spenso::{name}");
        let palette = greek_palette();
        let representation =
            LibraryRep::new_self_dual_with_index_palette(name, palette.clone()).unwrap();

        assert_eq!(
            LibraryRep::new_self_dual_with_index_palette(&qualified_name, palette).unwrap(),
            representation
        );
        assert_eq!(
            LibraryRep::new_self_dual(name).unwrap(),
            representation,
            "the legacy constructor must look up a palette-bearing representation without requesting a numeric palette"
        );

        let dual_name = "CanonicalDualPaletteAliasForSpensoTests";
        let qualified_dual_name = format!("spenso::{dual_name}");
        let dual = LibraryRep::new_dual_with_index_palette(dual_name, greek_palette()).unwrap();
        assert_eq!(LibraryRep::new_dual(dual_name).unwrap(), dual);
        assert_eq!(
            LibraryRep::new_dual_with_index_palette(&qualified_dual_name, greek_palette()).unwrap(),
            dual
        );
    }

    #[test]
    fn dual_constructor_rejects_an_existing_inline_metric() {
        assert!(matches!(
            LibraryRep::new_dual_with_index_palette("mink", IndexPalette::Numeric),
            Err(RepLibraryError::AlreadyExistsDifferentType(_))
        ));
    }

    #[test]
    fn portable_metadata_hydrates_unregistered_self_dual_and_dualizable_reps() {
        let self_dual_symbol = portable_representation_symbol(
            "spenso_rep_metadata_tests::ImportedSelfDual",
            RepresentationClass::SelfDual,
            greek_palette(),
        );
        let self_dual = LibraryRep::try_from_symbol_coerced(self_dual_symbol).unwrap();
        assert!(matches!(self_dual, LibraryRep::SelfDual(_)));
        assert_eq!(self_dual.symbol(), self_dual_symbol);

        let dualizable_symbol = portable_representation_symbol(
            "spenso_rep_metadata_tests::ImportedDualizable",
            RepresentationClass::Dualizable,
            greek_palette(),
        );
        let dualizable = LibraryRep::try_from_symbol_coerced(dualizable_symbol).unwrap();
        assert!(matches!(dualizable, LibraryRep::Dualizable(index) if index > 0));
        assert_eq!(dualizable.symbol(), dualizable_symbol);
    }

    #[test]
    fn imported_inline_metric_requires_a_local_metric_function() {
        let symbol = portable_representation_symbol(
            "spenso_rep_metadata_tests::ImportedInlineMetric",
            RepresentationClass::InlineMetric,
            IndexPalette::Numeric,
        );

        assert!(matches!(
            LibraryRep::try_from_symbol_coerced(symbol),
            Err(RepresentationError::ImportedInlineMetricRequiresLocalRegistration(
                error_symbol
            )) if error_symbol == symbol
        ));
    }

    #[test]
    fn representation_typst_body_escapes_the_symbol_label() {
        let body = representation_typst_body(
            &IndexDisplay::symbol("M\"; raw-source").expect("valid display label"),
        );

        assert!(body.contains(r#"italic("M\"; raw-source")"#));
        assert!(!body.contains(r#""M"; raw-source"#));
    }

    #[test]
    fn mathematical_display_nodes_render_and_round_trip_as_safe_user_data() {
        let display = IndexDisplay::math(
            "add",
            vec![
                IndexDisplay::math(
                    "accent",
                    vec![
                        IndexDisplay::symbol("p").unwrap(),
                        IndexDisplay::symbol("⃗").unwrap(),
                    ],
                )
                .unwrap(),
                IndexDisplay::math(
                    "frac",
                    vec![IndexDisplay::symbol("q").unwrap(), IndexDisplay::Number(2)],
                )
                .unwrap(),
                IndexDisplay::text("soft").unwrap(),
            ],
        )
        .unwrap();

        assert_eq!(
            display.to_typst_source(),
            r#"accent(p,"⃗") + frac(q,2) + upright("soft")"#
        );
        let symbol = SymbolBuilder::new(NamespacedSymbol::parse(
            "spenso_rep_metadata_tests::MathDisplayRoundTrip",
        ))
        .with_user_data(display.symbol_user_data())
        .build()
        .unwrap();
        assert_eq!(IndexDisplay::from_symbol(symbol), Some(display));
    }

    #[test]
    fn unknown_math_display_heads_cannot_be_constructed() {
        assert!(IndexDisplay::math("raw-code", vec![]).is_err());
    }

    #[test]
    fn palette_changes_only_printing_not_numeric_index_identity() {
        let representation = LibraryRep::new_self_dual_with_index_palette(
            "spenso_rep_metadata_tests::NumericIdentity",
            greek_palette(),
        )
        .unwrap();
        let atom = representation.to_symbolic([Atom::num(4), Atom::num(3)]);
        let AtomView::Fun(function) = atom.as_view() else {
            panic!("representation should remain a function atom");
        };
        let mut arguments = function.iter();

        assert_eq!(arguments.next().unwrap(), Atom::num(4).as_view());
        assert_eq!(arguments.next().unwrap(), Atom::num(3).as_view());
    }
}
