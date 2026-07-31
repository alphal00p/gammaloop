use crate::{network::tags::SPENSO_TAG, structure::concrete_index::ConcreteIndex};
use ahash::HashMap;
use derive_more::Display;
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    id::ReplaceBuilder,
    printer::{CanonicalOrderingSettings, PrintOptions, PrintUserData},
    symbol,
};

extern crate derive_more;

use std::fmt::{Debug, Display};

use eyre::Result;

pub trait ReplaceBuilderExt {
    fn has_matches(&self) -> bool;
    // Returns true if the pattern matches entirely the given atom
    fn matches(self) -> bool;
}

impl<'a, 'b> ReplaceBuilderExt for ReplaceBuilder<'a, 'b> {
    fn has_matches(&self) -> bool {
        self.match_iter().next().is_some()
    }

    fn matches(self) -> bool {
        self.partial(false).has_matches()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct SpensoPrintSettings {
    pub with_dim: bool,
    pub parens: bool,
    pub commas: bool,
    pub index_subscripts: bool,
    pub symbol_scripts: bool,
}

impl From<SpensoPrintSettings> for HashMap<String, PrintUserData> {
    fn from(settings: SpensoPrintSettings) -> Self {
        HashMap::from_iter([(
            "spenso".to_string(),
            PrintUserData::Integer(settings.to_usize() as i64),
        )])
    }
}

impl From<&SpensoPrintSettings> for HashMap<String, PrintUserData> {
    fn from(settings: &SpensoPrintSettings) -> Self {
        (*settings).into()
    }
}

impl From<usize> for SpensoPrintSettings {
    fn from(flag: usize) -> Self {
        Self::from_usize(flag)
    }
}

impl SpensoPrintSettings {
    fn from_usize(x: usize) -> Self {
        Self {
            parens: (x & 0b00001) != 0,
            commas: (x & 0b00010) != 0,
            with_dim: (x & 0b00100) != 0,
            symbol_scripts: (x & 0b01000) != 0,
            index_subscripts: (x & 0b10000) != 0,
        }
    }

    fn to_usize(self) -> usize {
        (self.parens as usize)
            | ((self.commas as usize) << 1)
            | ((self.with_dim as usize) << 2)
            | ((self.symbol_scripts as usize) << 3)
            | ((self.index_subscripts as usize) << 4)
    }

    pub fn typst() -> Self {
        Self {
            parens: true,
            commas: false,
            with_dim: false,
            symbol_scripts: true,
            index_subscripts: true,
        }
    }

    pub fn is_typst(&self) -> bool {
        self == &Self::typst()
    }

    // x^-2*a^-2*b^-2 -> 1/(x * a * b)^2
    // x^-1*a^-1*b^-1 -> ((1/x)/a)/b
    pub fn compact() -> Self {
        Self {
            parens: true,
            commas: false,
            with_dim: false,
            symbol_scripts: false,
            index_subscripts: false,
        }
    }

    pub fn nice_symbolica(&self) -> PrintOptions {
        PrintOptions {
            custom_print_mode: self.into(),
            color_builtin_symbols: true,
            terms_on_new_line: true,
            max_line_length: Some(120),
            color_namespace: false,
            multiplication_operator: '·',
            hide_all_namespaces: true,
            color_top_level_sum: true,
            num_exp_as_superscript: true,
            ..Default::default()
        }
    }

    /// Use Symbolica's Typst arithmetic while retaining Spenso's tensor and index notation.
    pub fn typst_options() -> PrintOptions {
        PrintOptions {
            custom_print_mode: Self::typst().into(),
            ..PrintOptions::typst()
        }
    }
}

pub trait AtomCoreExt {
    fn to_bare_ordered_string(&self) -> String;

    fn is_upper(&self) -> bool;
    fn is_lower(&self) -> bool;
}

impl<A: AtomCore> AtomCoreExt for A {
    fn to_bare_ordered_string(&self) -> String {
        self.to_canonically_ordered_string(
            CanonicalOrderingSettings::new()
                .include_namespace(false)
                .include_attributes(false),
        )
    }

    fn is_upper(&self) -> bool {
        match self.as_atom_view() {
            AtomView::Fun(a) => a.get_symbol().has_tag(&SPENSO_TAG.upper),
            AtomView::Var(a) => a.get_symbol().has_tag(&SPENSO_TAG.upper),
            _ => false,
        }
    }

    fn is_lower(&self) -> bool {
        match self.as_atom_view() {
            AtomView::Fun(a) => a.get_symbol().has_tag(&SPENSO_TAG.lower),
            AtomView::Var(a) => a.get_symbol().has_tag(&SPENSO_TAG.lower),
            _ => false,
        }
    }
}
#[derive(
    Debug,
    Copy,
    Clone,
    Ord,
    PartialOrd,
    Eq,
    PartialEq,
    Hash,
    Display,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    bincode_trait_derive::BorrowDecodeFromDecode,
)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct SerializableSymbol {
    symbol: Symbol,
}

impl SerializableSymbol {
    pub fn get_id(&self) -> u32 {
        self.symbol.get_id()
    }

    pub fn get_name(&self) -> &str {
        self.symbol.get_name()
    }
}

impl Serialize for SerializableSymbol {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.symbol.get_name().serialize(serializer)
    }
}

impl<'d> Deserialize<'d> for SerializableSymbol {
    fn deserialize<D>(deserializer: D) -> Result<SerializableSymbol, D::Error>
    where
        D: serde::Deserializer<'d>,
    {
        let value = String::deserialize(deserializer)?;
        Ok(SerializableSymbol {
            symbol: symbol!(&value),
        })
    }
}

impl From<Symbol> for SerializableSymbol {
    fn from(value: Symbol) -> Self {
        Self { symbol: value }
    }
}

impl From<SerializableSymbol> for Symbol {
    fn from(value: SerializableSymbol) -> Self {
        value.symbol
    }
}

impl From<SerializableSymbol> for u32 {
    fn from(value: SerializableSymbol) -> Self {
        value.symbol.get_id()
    }
}

pub trait LogPrint {
    fn log_print(&self, max_line_length: Option<usize>) -> String;
}

impl<A: AtomCore> LogPrint for A {
    fn log_print(&self, max: Option<usize>) -> String {
        let mut settings = SpensoPrintSettings::compact().nice_symbolica();
        settings.max_line_length = max;
        self.printer(settings).to_string()
    }
}

pub fn atomic_expanded_label<I: IntoSymbol>(indices: &[ConcreteIndex], name: I) -> Atom {
    let id = name.ref_into_symbol();
    atomic_expanded_label_id(indices, id, &[])
}
#[cfg(feature = "shadowing")]
pub fn atomic_flat_label<I: IntoSymbol>(index: usize, name: I) -> Atom {
    let id = name.ref_into_symbol();
    atomic_flat_label_id(index, id)
}

#[allow(clippy::cast_possible_wrap)]
#[cfg(feature = "shadowing")]
pub fn atomic_flat_label_id(index: usize, id: Symbol) -> Atom {
    let mut value_builder = FunctionBuilder::new(id);
    value_builder = value_builder.add_arg(Atom::num(index as i64).as_atom_view());
    value_builder.finish()
}
#[cfg(feature = "shadowing")]
#[allow(clippy::cast_possible_wrap)]
pub fn atomic_expanded_label_id(indices: &[ConcreteIndex], name: Symbol, args: &[Atom]) -> Atom {
    let mut value_builder = FunctionBuilder::new(name);
    let mut index_func = FunctionBuilder::new(symbol!("cind"));
    for arg in args {
        value_builder = value_builder.add_arg(arg);
    }
    for &index in indices {
        index_func = index_func.add_arg(Atom::num(index as i64).as_atom_view());
    }

    let indices = index_func.finish();
    value_builder.add_arg(&indices).finish()
}

#[cfg(feature = "shadowing")]
pub trait IntoSymbol {
    fn ref_into_symbol(&self) -> Symbol;

    fn from_str(s: &str) -> Self;
}

#[cfg(feature = "shadowing")]
pub trait IntoArgs {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom>;
    fn args(&self) -> Vec<Atom> {
        self.ref_into_args().collect()
    }
    fn cooked_name(&self) -> std::string::String;
}

#[cfg(feature = "shadowing")]
impl IntoArgs for usize {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        std::iter::once(Atom::num(*self as i64))
    }
    fn cooked_name(&self) -> std::string::String {
        format!("{self}")
    }
}

#[cfg(feature = "shadowing")]
impl IntoArgs for () {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        std::iter::empty()
    }
    fn cooked_name(&self) -> std::string::String {
        "".into()
    }
}

#[derive(Debug, Default, Clone, Copy, Serialize, Deserialize)]
pub struct NoArgs;

impl Display for NoArgs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "")
    }
}

#[cfg(feature = "shadowing")]
impl IntoArgs for NoArgs {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        std::iter::empty()
    }
    fn cooked_name(&self) -> std::string::String {
        "".into()
    }
}

#[cfg(feature = "shadowing")]
impl IntoArgs for Atom {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        std::iter::once(self.clone())
    }

    fn cooked_name(&self) -> std::string::String {
        self.to_string()
    }
}

#[cfg(feature = "shadowing")]
impl IntoArgs for Vec<Atom> {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        self.iter().cloned()
    }

    fn cooked_name(&self) -> std::string::String {
        let init = "".into();
        self.iter()
            .fold(init, |acc, x| acc + x.to_string().as_str())
    }
}

#[cfg(feature = "shadowing")]
impl<const N: usize> IntoArgs for [Atom; N] {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        self.iter().cloned()
    }

    fn cooked_name(&self) -> std::string::String {
        let init = "".into();
        self.iter()
            .fold(init, |acc, x| acc + x.to_string().as_str())
    }
}

// #[cfg(feature = "shadowing")]
// impl IntoSymbol for String {
//     fn ref_into_symbol(&self) -> Symbol {
//         symbol!(self)
//     }

//     fn from_str(s: &str) -> Self {
//         s.into()
//     }
// }

#[cfg(feature = "shadowing")]
impl IntoSymbol for Symbol {
    fn ref_into_symbol(&self) -> Symbol {
        *self
    }

    fn from_str(s: &str) -> Self {
        symbol!(s)
    }
}

#[cfg(feature = "shadowing")]
impl IntoSymbol for SerializableSymbol {
    fn ref_into_symbol(&self) -> Symbol {
        self.symbol
    }

    fn from_str(s: &str) -> Self {
        Self { symbol: symbol!(s) }
    }
}

#[cfg(feature = "shadowing")]
impl IntoSymbol for std::string::String {
    fn ref_into_symbol(&self) -> Symbol {
        symbol!(self)
    }
    fn from_str(s: &str) -> Self {
        s.into()
    }
}

#[cfg(test)]
mod tests {
    use super::SpensoPrintSettings;
    use symbolica::printer::PrintUserData;

    #[test]
    fn typst_options_use_typst_mode() {
        assert!(SpensoPrintSettings::typst_options().mode.is_typst());
    }

    #[test]
    fn typst_options_include_spenso_typst_settings() {
        let options = SpensoPrintSettings::typst_options();
        let Some(PrintUserData::Integer(encoded)) = options.custom_print_mode.get("spenso") else {
            panic!("missing Spenso print settings");
        };

        assert_eq!(
            SpensoPrintSettings::from(usize::try_from(*encoded).unwrap()),
            SpensoPrintSettings::typst()
        );
    }
}
