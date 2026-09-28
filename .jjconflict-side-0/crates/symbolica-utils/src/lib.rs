use std::fmt::{Debug, Display};

use itertools::Itertools;
use serde::{Deserialize, Serialize};
use symbolica::{
    domains::finite_field::PrimeIteratorU64, id::Context, prelude::*, printer::PrintUserData,
    utils::Settable,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TypstMode {
    Math,
    Markup,
    Code,
}

pub struct TypstModeSettings {
    mode: TypstMode,
}

impl TypstModeSettings {
    pub fn new(mode: &PrintUserData) -> Self {
        let mode = match mode {
            PrintUserData::Integer(0) => TypstMode::Math,
            PrintUserData::Integer(1) => TypstMode::Markup,
            PrintUserData::Integer(2) => TypstMode::Code,
            PrintUserData::String(mode) if mode.eq_ignore_ascii_case("math") => TypstMode::Math,
            PrintUserData::String(mode) if mode.eq_ignore_ascii_case("markup") => TypstMode::Markup,
            PrintUserData::String(mode) if mode.eq_ignore_ascii_case("code") => TypstMode::Code,
            _ => TypstMode::Math,
        };

        Self { mode }
    }

    pub fn to_user_data(&self) -> PrintUserData {
        PrintUserData::Integer(self.mode as i64)
    }
}

pub trait PrintSettingsExt {
    fn typst_mode(&self) -> Option<TypstMode>;
    fn set_typst_mode(&mut self, mode: TypstMode);
}

impl PrintSettingsExt for PrintOptions {
    fn typst_mode(&self) -> Option<TypstMode> {
        if self.mode.is_typst() {
            if let Some(a) = self.custom_print_mode.get("typst") {
                Some(TypstModeSettings::new(a).mode)
            } else {
                Some(TypstMode::Math)
            }
        } else {
            None
        }
    }

    fn set_typst_mode(&mut self, mode: TypstMode) {
        self.custom_print_mode.insert(
            "typst".to_string(),
            TypstModeSettings { mode }.to_user_data(),
        );
    }
}

pub trait ReplaceBuilderExt {
    fn has_matches(&self) -> bool;

    /// Returns true if the pattern matches the whole atom.
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

pub trait AtomPrintExt {
    fn to_bare_ordered_string(&self) -> String;
}

impl<A: AtomCore> AtomPrintExt for A {
    fn to_bare_ordered_string(&self) -> String {
        self.to_canonically_ordered_string(
            CanonicalOrderingSettings::new()
                .include_namespace(false)
                .include_attributes(false),
        )
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

impl Display for SerializableSymbol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.symbol, f)
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

pub trait IntoSymbol {
    fn ref_into_symbol(&self) -> Symbol;

    fn from_str(s: &str) -> Self;
}

pub trait IntoArgs {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom>;

    fn args(&self) -> Vec<Atom> {
        self.ref_into_args().collect()
    }

    fn cooked_name(&self) -> String;
}

impl IntoArgs for usize {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        std::iter::once(Atom::num(*self as i64))
    }

    fn cooked_name(&self) -> String {
        format!("{self}")
    }
}

impl IntoArgs for () {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        std::iter::empty()
    }

    fn cooked_name(&self) -> String {
        String::new()
    }
}

#[derive(Debug, Default, Clone, Copy, Serialize, Deserialize)]
pub struct NoArgs;

impl Display for NoArgs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "")
    }
}

impl IntoArgs for NoArgs {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        std::iter::empty()
    }

    fn cooked_name(&self) -> String {
        String::new()
    }
}

impl IntoArgs for Atom {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        std::iter::once(self.clone())
    }

    fn cooked_name(&self) -> String {
        self.to_string()
    }
}

impl IntoArgs for Vec<Atom> {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        self.iter().cloned()
    }

    fn cooked_name(&self) -> String {
        self.iter().map(ToString::to_string).join("")
    }
}

impl<const N: usize> IntoArgs for [Atom; N] {
    fn ref_into_args(&self) -> impl Iterator<Item = Atom> {
        self.iter().cloned()
    }

    fn cooked_name(&self) -> String {
        self.iter().map(ToString::to_string).join("")
    }
}

impl IntoSymbol for Symbol {
    fn ref_into_symbol(&self) -> Symbol {
        *self
    }

    fn from_str(s: &str) -> Self {
        symbol!(s)
    }
}

impl IntoSymbol for SerializableSymbol {
    fn ref_into_symbol(&self) -> Symbol {
        self.symbol
    }

    fn from_str(s: &str) -> Self {
        Self { symbol: symbol!(s) }
    }
}

impl IntoSymbol for String {
    fn ref_into_symbol(&self) -> Symbol {
        symbol!(self)
    }

    fn from_str(s: &str) -> Self {
        s.into()
    }
}

pub trait IsNeg {
    fn is_negative(&self) -> bool;
}

impl<A: AtomCore> IsNeg for A {
    fn is_negative(&self) -> bool {
        match self.as_atom_view() {
            AtomView::Num(a) => match a.get_coeff_view() {
                CoefficientView::FiniteField(_, _) => false,
                CoefficientView::Float(re, im) => im.is_zero() && re.to_float().is_negative(),
                CoefficientView::Indeterminate => false,
                CoefficientView::Natural(n_re, _d_re, n_im, _de_im) => {
                    n_im == 0 && n_re.is_negative()
                }
                CoefficientView::Infinity(_) => false,
                CoefficientView::Large(re, im) => {
                    let re = re.to_rat();
                    let im = im.to_rat();
                    im.is_zero() && re.is_negative()
                }
                CoefficientView::RationalPolynomial(_) => false,
            },
            AtomView::Mul(a) => a.iter().any(|term| term.is_negative()),
            _ => false,
        }
    }
}

pub trait PrimeGenerate {
    fn prime_generate_int(sample_iterator: &mut PrimeIteratorU64) -> Atom;
    fn prime_generate_ff_complex(
        sample_iterator: &mut PrimeIteratorU64,
        finite_field: &Zp64,
    ) -> Atom;
    fn prime_generate_ff(sample_iterator: &mut PrimeIteratorU64, finite_field: &Zp64) -> Atom;
    fn prime_generate_rat_complex(sample_iterator: &mut PrimeIteratorU64) -> Atom;
    fn prime_generate_rat(sample_iterator: &mut PrimeIteratorU64) -> Atom;
}

impl PrimeGenerate for Atom {
    fn prime_generate_int(sample_iterator: &mut PrimeIteratorU64) -> Atom {
        Atom::num(sample_iterator.next().unwrap())
    }

    fn prime_generate_rat(sample_iterator: &mut PrimeIteratorU64) -> Atom {
        let a = Rational::from((
            sample_iterator.next().unwrap(),
            sample_iterator.next().unwrap(),
        ));
        Atom::num(a)
    }

    fn prime_generate_ff_complex(
        sample_iterator: &mut PrimeIteratorU64,
        finite_field: &Zp64,
    ) -> Atom {
        let _a = finite_field.to_element(sample_iterator.next().unwrap());
        let _b = finite_field.to_element(sample_iterator.next().unwrap());
        Atom::num(1)
    }

    fn prime_generate_ff(sample_iterator: &mut PrimeIteratorU64, _finite_field: &Zp64) -> Atom {
        let a = Rational::from((
            sample_iterator.next().unwrap(),
            sample_iterator.next().unwrap(),
        ));
        let b = Rational::from((
            sample_iterator.next().unwrap(),
            sample_iterator.next().unwrap(),
        ));
        Atom::num(Complex::new(a, b))
    }

    fn prime_generate_rat_complex(sample_iterator: &mut PrimeIteratorU64) -> Atom {
        let a = Rational::from((
            sample_iterator.next().unwrap(),
            sample_iterator.next().unwrap(),
        ));
        let b = Rational::from((
            sample_iterator.next().unwrap(),
            sample_iterator.next().unwrap(),
        ));
        Atom::num(Complex::new(a, b))
    }
}

pub trait Replaces {
    fn replace_with<R: Into<ReplaceWith<'static>>>(&self, rhs: R) -> Replacement;
}

impl<A: AtomCore> Replaces for A {
    fn replace_with<R: Into<ReplaceWith<'static>>>(&self, rhs: R) -> Replacement {
        Replacement::new(self.to_pattern(), rhs)
    }
}

pub trait AtomFloatExt {
    fn floatify(&self, prec: u32) -> Atom;
}

impl<A: AtomCore<Output = Atom>> AtomFloatExt for A {
    fn floatify(&self, prec: u32) -> Atom {
        self.map_coefficient(|c| match c {
            CoefficientView::Natural(r, d, ri, di) => {
                if (ri == 0 || di == 1) && (r == 0 || d == 1) {
                    Coefficient::Complex(Complex::new(Rational::new(r, d), Rational::new(ri, di)))
                } else {
                    Coefficient::Float(Complex::new(
                        Rational::from((r, d)).to_multi_prec_float(prec),
                        Rational::from((ri, di)).to_multi_prec_float(prec),
                    ))
                }
            }
            CoefficientView::Large(r, ri) => Coefficient::Float(Complex::new(
                r.to_rat().to_multi_prec_float(prec),
                ri.to_rat().to_multi_prec_float(prec),
            )),
            _ => c.to_owned(),
        })
    }
}

pub trait PatternReplacement {
    fn replace_multiple_repeat<T: BorrowReplacement>(&self, replacements: &[T]) -> Self;
    fn replace_multiple_mut<T: BorrowReplacement>(&mut self, replacements: &[T]);
    fn replace_multiple_repeat_mut<T: BorrowReplacement>(&mut self, replacements: &[T]);
    fn replace_map_mut<F: Fn(AtomView, &Context, &mut Settable<'_, Atom>)>(&mut self, m: &F);
}

impl PatternReplacement for Atom {
    fn replace_multiple_mut<T: BorrowReplacement>(&mut self, replacements: &[T]) {
        *self = self
            .as_atom_view()
            .replace_multiple(replacements.iter().map(BorrowReplacement::borrow));
    }

    fn replace_multiple_repeat<T: BorrowReplacement>(&self, replacements: &[T]) -> Atom {
        let mut out = self.clone();
        let mut out_mut = out.clone();

        while out.replace_multiple_into(
            replacements.iter().map(BorrowReplacement::borrow),
            &mut out_mut,
        ) {
            if out == out_mut {
                break;
            }
            std::mem::swap(&mut out, &mut out_mut)
        }
        out
    }

    fn replace_multiple_repeat_mut<T: BorrowReplacement>(&mut self, replacements: &[T]) {
        let atom = self.replace_multiple_repeat(replacements);
        *self = atom;
    }

    fn replace_map_mut<F: Fn(AtomView, &Context, &mut Settable<'_, Atom>)>(&mut self, m: &F) {
        *self = self.replace_map(m);
    }
}
