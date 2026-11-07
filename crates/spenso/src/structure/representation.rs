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
use symbolica_notation::{IndexDisplay, IndexDisplayError};

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

#[cfg(feature = "shadowing")]
use crate::network::tags::mathml;

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
    /// For example, see [`Representation::negative`].
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
        self.rep.cmp(&other.rep).then(self.dim.cmp(&other.dim))
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
    /// Construct compact `dot(a(rep), b(rep))` notation from unrepresented
    /// vector function heads (or vector symbols).
    ///
    /// This appends the representation to both operands. Inputs which already
    /// carry a representation or explicit index do not satisfy the precondition;
    /// applying this method to them would append a second structure argument.
    pub fn inner_product<'a, It: Into<AtomOrView<'a>>>(&self, a: It, b: It) -> Atom {
        let a: AtomOrView<'a> = a.into();
        let b: AtomOrView<'a> = b.into();
        function!(
            SPENSO_TAG.dot,
            self.vector(a.as_view(), []),
            self.vector(b.as_view(), [])
        )
    }

    /// Attach this representation, optionally with an index, to a vector name.
    ///
    /// For example, `p` becomes `p(mink(D,mu))`, while `Q(1)` becomes
    /// `Q(1,mink(D,mu))`. An empty index list produces the compact vector used
    /// by scalar products and Dirac slashes. Non-function expressions are
    /// retained, matching the scalar-product convention.
    #[cfg(feature = "shadowing")]
    pub fn vector(&self, value: AtomView<'_>, indices: impl IntoIterator<Item = Atom>) -> Atom {
        let rep = self.to_symbolic(indices);
        match value {
            AtomView::Fun(fun) => FunctionBuilder::new(fun.get_symbol())
                .add_args(fun.iter())
                .add_arg(rep)
                .finish(),
            AtomView::Var(var) => FunctionBuilder::new(var.get_symbol()).add_arg(rep).finish(),
            _ => value.to_owned(),
        }
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

#[cfg(feature = "shadowing")]
fn representation_label_name(name: &str) -> &str {
    name.rsplit("::").next().unwrap_or(name)
}

#[cfg(feature = "shadowing")]
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

#[cfg(feature = "shadowing")]
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

#[cfg(feature = "shadowing")]
fn representation_typst_body(label: &IndexDisplay) -> String {
    format!(
        "(dim, ind ) = (content: $ {}^#dim _#ind $, upper:true)",
        label.to_typst_source()
    )
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

    #[cfg(feature = "shadowing")]
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

        SymbolBuilder::new(namespaced)
            .with_tags(tags)
            .with_user_data(metadata.to_user_data())
            .with_print_function(move |a, opt, state| {
                if opt.mode.is_mathml() {
                    return mathml::render(a, opt, state, &|_, _| None);
                }
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
                        if let Some(display) = palette_index
                            .and_then(|index| palette.as_ref()?.resolve(index))
                            .or_else(|| match ind {
                                AtomView::Var(variable) => {
                                    IndexDisplay::from_symbol(variable.get_symbol())
                                }
                                _ => None,
                            })
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

                // Ordinary expressions use Symbolica's namespace, bracket,
                // and color handling, including the current nested print state.
                None
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
    #[error(transparent)]
    IndexDisplay(#[from] IndexDisplayError),
    #[error("{0} Already exists and is of different type")]
    AlreadyExistsDifferentType(String),
    #[error("{0} Already exists and has different metric function")]
    AlreadyExistsDifferentMetric(String),
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
/// Can construct a `Representation` from an `AtomView` of the form
/// `<representation>(<dimension>)`.
///
impl<'a, T: RepName> TryFrom<AtomView<'a>> for Representation<T> {
    type Error = SlotError;

    fn try_from(value: AtomView<'a>) -> Result<Self, Self::Error> {
        Self::parse_with(value, |head, wrapper| match wrapper {
            Some(wrapper) => T::try_from_symbol(head, wrapper),
            None => T::try_from_symbol_coerced(head),
        })
    }
}

#[cfg(feature = "shadowing")]
impl<T: RepName> Representation<T> {
    /// Read the existing representation-prefix grammar with a caller-owned
    /// resolver. Dimension validation and acceptance of trailing arguments stay
    /// identical for direct conversion and operation-scoped cached recognition.
    pub(crate) fn parse_with(
        value: AtomView<'_>,
        mut resolve: impl FnMut(Symbol, Option<Symbol>) -> Result<T, RepresentationError>,
    ) -> Result<Self, SlotError> {
        let (rep, mut iter) = if let AtomView::Fun(f) = value {
            let name = f.get_symbol();

            let innerf = f.iter().next().ok_or(SlotError::Composite)?;

            if let AtomView::Fun(innerf) = innerf {
                let rep = resolve(innerf.get_symbol(), Some(name)).map_err(SlotError::RepError)?;

                (rep, innerf.iter())
            } else {
                let rep = resolve(name, None).map_err(SlotError::RepError)?;
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
    use std::cmp::Ordering;

    #[cfg(feature = "shadowing")]
    use symbolica::symbol;

    use crate::structure::dimension::Dimension;
    use crate::structure::representation::{RepName, Representation};

    use super::{Dummy, Euclidean, LibraryRep, Lorentz, Minkowski};

    #[test]
    fn representation_order_is_total_and_consistent_with_equality() {
        let dimensions = [
            Dimension::Concrete(0),
            Dimension::Concrete(1),
            Dimension::Concrete(4),
            #[cfg(feature = "shadowing")]
            Dimension::from(symbol!("representation_order_a")),
            #[cfg(feature = "shadowing")]
            Dimension::from(symbol!("representation_order_b")),
        ];
        let kinds = [
            LibraryRep::Dummy,
            Euclidean {}.into(),
            Minkowski {}.into(),
            Lorentz {}.into(),
            Lorentz {}.dual().into(),
        ];
        let representations: Vec<_> = kinds
            .into_iter()
            .flat_map(|rep| dimensions.map(|dim| Representation { rep, dim }))
            .collect();

        for left in &representations {
            for right in &representations {
                assert_eq!(left.cmp(right) == Ordering::Equal, left == right);
                assert_eq!(left.partial_cmp(right), Some(left.cmp(right)));
                assert_eq!(left.cmp(right), right.cmp(left).reverse());
                for third in &representations {
                    if left <= right && right <= third {
                        assert!(left <= third, "{left:?} <= {right:?} <= {third:?}");
                    }
                }
            }
        }

        let small = Dummy {}.new_rep(1);
        let large = Dummy {}.new_rep(4);
        assert!(small < large);
        assert!(small.cast::<LibraryRep>() < large.cast::<LibraryRep>());
    }

    #[test]
    fn sorting_respects_kind_orientation_and_dimension_including_dummy() {
        let expected = [
            LibraryRep::Dummy.new_rep(0),
            LibraryRep::Dummy.new_rep(4),
            Euclidean {}.new_rep(0).cast(),
            Euclidean {}.new_rep(4).cast(),
            Minkowski {}.new_rep(2).cast(),
            Lorentz {}.new_rep(1).cast(),
            Lorentz {}.new_rep(4).cast(),
            Lorentz {}.dual().new_rep(1).cast(),
            Lorentz {}.dual().new_rep(4).cast(),
        ];
        assert!(Permutation::sort(expected).is_identity());
        for rotation in 0..expected.len() {
            let mut input = expected;
            input.reverse();
            input.rotate_left(rotation);
            input.sort();
            assert_eq!(input, expected);
        }
    }

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
    use symbolica_notation::IndexDisplay;

    use crate::network::tags::SPENSO_TAG;

    use super::{
        IndexPalette, IndexRow, LibraryRep, Minkowski, RepLibraryError, RepName,
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
    fn ordinary_representation_calls_reuse_symbolica_namespace_and_color_handling() {
        use symbolica::{
            atom::AtomCore,
            function,
            printer::{ColorMode, PrintOptions},
            symbol,
        };

        let dimension = function!(
            symbol!("representation_print_test::dim"),
            symbol!("representation_print_test::D")
        );
        let representation = LibraryRep::from(Minkowski {}).to_symbolic([dimension.clone()]);
        let ordinary = function!(symbol!("representation_print_control::mink"), dimension);
        let outer = symbol!("representation_print_test::P");
        let actual = function!(outer, representation);
        let control = function!(outer, ordinary);
        let canonical = actual.to_canonical_string();

        for brackets in [('(', ')'), ('[', ']')] {
            let options = PrintOptions {
                function_brackets: brackets,
                hide_all_namespaces: true,
                color_mode: ColorMode::Always,
                ..PrintOptions::new()
            };
            let actual_text = actual.printer(options.clone()).to_string();
            let expected_text = control.printer(options).to_string();
            assert!(!actual_text.contains("::"));
            assert!(actual_text.contains("\x1b["));
            assert_eq!(actual_text, expected_text);
        }
        let explicit = actual
            .printer(PrintOptions {
                hide_all_namespaces: false,
                color_mode: ColorMode::Never,
                ..PrintOptions::new()
            })
            .to_string();
        assert!(explicit.contains("spenso::mink("));
        assert!(explicit.contains("representation_print_test::D"));
        assert_eq!(actual.to_canonical_string(), canonical);
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
