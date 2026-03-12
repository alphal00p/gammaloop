use super::{
    abstract_index::{AbstractIndex, AbstractIndexError},
    concrete_index::ConcreteIndex,
    dimension::{Dimension, DimensionError},
    slot::Slot,
};
use ahash::AHashMap;
use append_only_vec::AppendOnlyVec;
use linnet::half_edge::involution::Orientation;
use once_cell::sync::Lazy;
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
    network::{library::symbolic::ETS, parsing::SPENSO_TAG},
    structure::{abstract_index::AIND_SYMBOLS, slot::SlotError},
};

#[cfg(feature = "shadowing")]
use symbolica::{
    atom::{Atom, AtomOrView, AtomView, FunctionBuilder, Symbol},
    function, get_symbol, symbol,
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
        let a: AtomOrView<'a> = a.into();
        let b: AtomOrView<'a> = b.into();
        function!(ETS.metric, self.to_symbolic([]), a.as_view(), b.as_view())
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

pub(crate) static REPS: Lazy<RwLock<ExtendibleReps>> =
    Lazy::new(|| RwLock::new(ExtendibleReps::new()));
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

impl LibraryRep {
    #[cfg(feature = "shadowing")]
    pub fn new_symbol(&self, name: &str) -> Symbol {
        if let Some(s) = get_symbol!(name) {
            s
        } else {
            use symbolica::{atom::AtomCore, printer::PrintState};

            use crate::shadowing::symbolica_utils::SpensoPrintSettings;

            let body = format!(
                "(dim, ind ) = (content: $ \"{}\"^#dim _#ind $, upper:true)",
                name
            );

            let (rep_name, tages) = match self {
                LibraryRep::SelfDual(a) => (
                    encode_base(*a as usize, &LATIN),
                    &[
                        // &SPENSO_TAG.upper,
                        &SPENSO_TAG.representation,
                        &SPENSO_TAG.self_dual,
                    ],
                ),
                LibraryRep::Dualizable(a) => (
                    encode_base(a.unsigned_abs() as usize, &CYRILLIC),
                    &[
                        // &SPENSO_TAG.upper,
                        &SPENSO_TAG.representation,
                        &SPENSO_TAG.dualizable,
                    ],
                ),
                LibraryRep::InlineMetric(a) => (
                    encode_base(*a as usize, &GREEK),
                    &[
                        // &SPENSO_TAG.upper,
                        &SPENSO_TAG.representation,
                        &SPENSO_TAG.self_dual,
                    ],
                ),
                LibraryRep::Dummy => (
                    String::new(),
                    &[
                        // &SPENSO_TAG.upper,
                        &SPENSO_TAG.representation,
                        &SPENSO_TAG.self_dual,
                    ],
                ),
            };

            let name = name.to_string();

            symbol!(
                &name,
                print = move |a, opt| {
                    match opt.custom_print_mode {
                        Some(("typst", 1)) => Some(body.clone()),

                        Some(("spenso", i)) => {
                            let SpensoPrintSettings {
                                with_dim,
                                commas,
                                parens,
                                index_subscripts,
                                ..
                            } = SpensoPrintSettings::from(i);
                            let AtomView::Fun(f) = a else {
                                return None;
                            };

                            let mut out = if opt.color_builtin_symbols {
                                nu_ansi_term::Color::DarkGray.paint(&rep_name).to_string()
                            } else {
                                rep_name.clone()
                            };

                            if index_subscripts {
                                out.push('_');
                            }
                            if parens && index_subscripts {
                                out.push('(');
                            }

                            if f.get_nargs() == 2 {
                                let mut arg_iter = f.iter();
                                let dim = arg_iter.next()?;

                                if with_dim {
                                    dim.format(&mut out, opt, PrintState::new()).unwrap();
                                    if commas {
                                        out.push(',');
                                    } else {
                                        out.push(' ');
                                    }
                                }
                                let ind = arg_iter.next()?;

                                ind.format(&mut out, opt, PrintState::new()).unwrap();
                                if parens && index_subscripts {
                                    out.push(')');
                                }

                                return Some(out);
                            }

                            None
                        }
                        _ => {
                            let AtomView::Fun(f) = a else {
                                return None;
                            };

                            let mut out = if opt.color_builtin_symbols {
                                nu_ansi_term::Color::DarkGray.paint(&name).to_string()
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
                        }
                    }
                },
                norm = |ain, out| {
                    let AtomView::Fun(f) = ain else {
                        return;
                    };

                    // println!("Normalizing representation function: {:#}", ain);

                    let name = f.get_symbol();

                    if f.get_nargs() == 3 {
                        let mut args = f.iter();
                        let dim = args.next().unwrap();
                        let a = args.next().unwrap();
                        let b = args.next().unwrap();

                        // A rep is linear in each index slot;
                        match a {
                            AtomView::Add(a) => {
                                let builder = FunctionBuilder::new(name);

                                let mut new_out = Atom::Zero;
                                for t in a.iter() {
                                    new_out += builder.clone().add_args(&[dim, t, b]).finish();
                                }
                                **out = new_out;
                                return;
                            }

                            AtomView::Mul(a) => {
                                let builder = FunctionBuilder::new(name);

                                let mut new_out = Atom::num(1);
                                let mut coef = Atom::num(1);
                                for t in a.iter() {
                                    if let AtomView::Num(_) = t {
                                        coef *= t;
                                    } else {
                                        new_out *= t;
                                    }
                                }

                                if coef != Atom::num(1) {
                                    **out = coef
                                        * builder.add_args(&[dim, new_out.as_view(), b]).finish();
                                    return;
                                }
                            }

                            _ => {}
                        }

                        match b {
                            AtomView::Add(b) => {
                                let builder = FunctionBuilder::new(name);

                                let mut new_out = Atom::Zero;
                                for t in b.iter() {
                                    new_out += builder.clone().add_args(&[dim, a, t]).finish();
                                }
                                **out = new_out;
                                return;
                            }

                            AtomView::Mul(b) => {
                                let builder = FunctionBuilder::new(name);

                                let mut new_out = Atom::num(1);
                                let mut coef = Atom::num(1);
                                for t in b.iter() {
                                    if let AtomView::Num(_) = t {
                                        coef *= t;
                                    } else {
                                        new_out *= t;
                                    }
                                }

                                if coef != Atom::num(1) {
                                    **out = coef
                                        * builder.add_args(&[dim, a, new_out.as_view()]).finish();
                                    return;
                                }
                            }

                            AtomView::Pow(_) => {
                                panic!("Powers are not supported in index slots {ain}");
                            }
                            _ => {}
                        }

                        let a_build = match a {
                            AtomView::Fun(f) => {
                                let sym = f.get_symbol();
                                if sym.get_wildcard_level() > 0 {
                                    return;
                                }

                                if sym.has_tag(&SPENSO_TAG.index) {
                                    None
                                } else if sym.has_tag(&SPENSO_TAG.representation) {
                                    return;
                                } else {
                                    Some((sym, f.iter().collect::<Vec<_>>()))
                                }
                            }
                            AtomView::Var(v) => {
                                let sym = v.get_symbol();
                                if sym.get_wildcard_level() > 0 {
                                    return;
                                }
                                if sym.has_tag(&SPENSO_TAG.index) {
                                    None
                                } else if sym.has_tag(&SPENSO_TAG.representation) {
                                    return;
                                } else {
                                    Some((sym, vec![]))
                                }
                            }
                            _ => None,
                        };

                        let b_build = match b {
                            AtomView::Fun(f) => {
                                let sym = f.get_symbol();
                                if sym.get_wildcard_level() > 0 {
                                    return;
                                }
                                if sym.has_tag(&SPENSO_TAG.index) {
                                    None
                                } else if sym.has_tag(&SPENSO_TAG.representation) {
                                    return;
                                } else {
                                    Some((sym, f.iter().collect::<Vec<_>>()))
                                }
                            }
                            AtomView::Var(v) => {
                                let sym = v.get_symbol();
                                if sym.get_wildcard_level() > 0 {
                                    return;
                                }
                                // sym
                                if sym.has_tag(&SPENSO_TAG.index) {
                                    None
                                } else if sym.has_tag(&SPENSO_TAG.representation) {
                                    return;
                                } else {
                                    Some((sym, vec![]))
                                }
                            }
                            _ => None,
                        };

                        match (a_build, b_build) {
                            (Some((a, args)), Some((b, args_b))) => {
                                // println!(
                                //     "Both are functions or variables, so this is a dot product;"
                                // );
                                let builder = FunctionBuilder::new(ETS.metric);
                                let rep = FunctionBuilder::new(name).add_arg(dim).finish();

                                let a = if args.is_empty() {
                                    Atom::var(a)
                                } else {
                                    FunctionBuilder::new(a).add_args(&args).finish()
                                };

                                let b = if args_b.is_empty() {
                                    Atom::var(b)
                                } else {
                                    FunctionBuilder::new(b).add_args(&args_b).finish()
                                };

                                let atom = builder.add_args(&[rep, a, b]).finish();
                                // println!("{atom}");
                                **out = atom;
                            }
                            (None, Some((b, args))) => {
                                // println!("a is an index and b is a tensor, we just append to b");
                                let index =
                                    FunctionBuilder::new(name).add_arg(dim).add_arg(a).finish();
                                // println!("index: {index}");
                                let atom = FunctionBuilder::new(b)
                                    .add_args(&args)
                                    .add_arg(index)
                                    .finish();
                                // println!("atom:{atom}");
                                **out = atom;
                            }
                            (Some((a, args)), None) => {
                                // // println!(
                                //     "b is a dual index (in second position) and a is a tensor, we just append to a"
                                // );
                                let index =
                                    // FunctionBuilder::new(AIND_SYMBOLS.dind)
                                    // .add_arg(
                                        FunctionBuilder::new(name).add_arg(dim).add_arg(b).finish();
                                // );
                                // .finish();
                                // println!("index: {index}");
                                let atom = FunctionBuilder::new(a)
                                    .add_args(&args)
                                    .add_arg(index)
                                    .finish();
                                // println!("atom: {atom}");
                                **out = atom;
                            }
                            (None, None) => {
                                // println!("both are indices so this is a metric shorthand");
                                let index_builder = FunctionBuilder::new(name).add_arg(dim);

                                let atom = FunctionBuilder::new(ETS.metric)
                                    .add_args(&[
                                        index_builder.clone().add_arg(a).finish(),
                                        index_builder.add_arg(b).finish(),
                                    ])
                                    .finish();

                                // println!("atom:{atom}");
                                **out = atom;
                            }
                        }
                    }
                    // println!("DOne");
                    // else {
                    //     // println!("Does not have 3 arguments: {}", f.get_nargs());
                    // }
                },
                tags = tages
            )
        }
    }

    pub fn new_dual(name: &str) -> Result<Self, RepLibraryError> {
        REPS.write().unwrap().new_dual_impl(name)
    }

    #[cfg(feature = "shadowing")]
    pub fn symbol(&self) -> Symbol {
        REPS.read().unwrap()[*self].symbol
    }

    pub fn name(&self) -> String {
        REPS.read().unwrap()[*self].name.clone()
    }

    pub fn new_self_dual(name: &str) -> Result<Self, RepLibraryError> {
        REPS.write().unwrap().new_self_dual(name)
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
    symbol: symbol!("Dummy"),
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
}

impl ExtendibleReps {
    pub fn reps(&self) -> impl Iterator<Item = &LibraryRep> {
        self.name_map.values()
    }

    pub fn new_dual_impl(&mut self, name: &str) -> Result<LibraryRep, RepLibraryError> {
        if let Some(rep) = self.name_map.get(name) {
            if let LibraryRep::SelfDual(_) = rep {
                return Err(RepLibraryError::AlreadyExistsDifferentType(name.into()));
            } else {
                return Ok(*rep);
            }
        }

        let rep = LibraryRep::Dualizable(DUALIZABLE.len() as i16 + 1);

        self.name_map.insert(name.into(), rep);
        #[cfg(feature = "shadowing")]
        let symbol = rep.new_symbol(name);
        #[cfg(feature = "shadowing")]
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

    pub fn new_dual(name: &str) -> Result<LibraryRep, RepLibraryError> {
        REPS.write().unwrap().new_dual_impl(name)
    }

    pub fn new_self_dual(&mut self, name: &str) -> Result<LibraryRep, RepLibraryError> {
        if let Some(rep) = self.name_map.get(name) {
            if let LibraryRep::Dualizable(_) = rep {
                return Err(RepLibraryError::AlreadyExistsDifferentType(name.into()));
            } else {
                return Ok(*rep);
            }
        }

        let rep = LibraryRep::SelfDual(SELF_DUAL.len() as u16);
        self.name_map.insert(name.into(), rep);
        #[cfg(feature = "shadowing")]
        let symbol = rep.new_symbol(name);
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
        if let Some(rep) = self.name_map.get(name) {
            match rep {
                LibraryRep::SelfDual(_) | LibraryRep::Dualizable(_) | LibraryRep::Dummy => {
                    return Err(RepLibraryError::AlreadyExistsDifferentType(name.into()));
                }
                LibraryRep::InlineMetric(a) => {
                    if INLINE_METRIC[*a as usize].1.metric_data == metric_fn {
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
        self.name_map.insert(name.into(), rep);
        #[cfg(feature = "shadowing")]
        let symbol = rep.new_symbol(name);
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

        #[cfg(feature = "shadowing")]
        let _ = ETS.metric;

        #[cfg(feature = "shadowing")]
        let _ = AIND_SYMBOLS.aind;
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

        let rep = REPS
            .read()
            .unwrap()
            .find_symbol(sym)
            .ok_or(RepresentationError::NotRepresentationError(sym))?;

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
        REPS.read()
            .unwrap()
            .find_symbol(sym)
            .ok_or(RepresentationError::NotRepresentationError(sym))
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
    /// yields a function builder for the representation, adding a first variable: the dimension.
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
    use crate::network::parsing::SPENSO_TAG;
    use crate::structure::representation::initialize;
    use symbolica::{atom::AtomCore, parse, parse_lit, symbol};

    #[test]
    fn normalization_shortcuts() {
        initialize();
        let mu = symbol!("mu", tags = [&SPENSO_TAG.index]);
        let a = parse!("mink(4,a,mu)");
        let b = parse!("mink(4,mu,a)");
        assert_eq!(a, b);
        let a = parse!("lor(4,a,mu)");
        let b = parse!("lor(4,mu,a)");

        let a = parse!("mink(4,mu,mink(4,mu))");
        // assert_ne!(a, b);
    }
    // use symbolica::symbol;

    // use crate::structure::representation::BaseRepName;

    // use super::Lorentz;

    // #[test]
    // fn rep_pattern() {
    //     println!("{}", Dual::<Lorentz>::pattern(symbol!("d_")));
    //     println!(
    //         "{}",
    //         Dual::<Lorentz>::rep(3).to_pattern_wrapped(symbol!("d_"))
    //     );
    // }
}
