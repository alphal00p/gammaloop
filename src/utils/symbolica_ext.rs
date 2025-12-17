use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use schemars::{JsonSchema, json_schema};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Display,
    ops::{Deref, DerefMut},
    sync::LazyLock,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    coefficient::CoefficientView,
    domains::{
        algebraic_number::AlgebraicExtension,
        finite_field::{FiniteFieldCore, PrimeIteratorU64, Zp64},
        float::{Complex, FloatLike},
        integer::IntegerRing,
        rational::{FractionField, Q, Rational},
    },
    parse,
    poly::polynomial::PolynomialRing,
    printer::{PrintMode, PrintOptions, PrintState},
};

use crate::GammaLoopContext;

use super::{GS, W_};

pub static Q_I: LazyLock<AlgebraicExtension<FractionField<IntegerRing>>> =
    LazyLock::new(|| AlgebraicExtension::new_complex(Q));

pub static COMPLEXRATPOLYFIELD: LazyLock<
    FractionField<PolynomialRing<AlgebraicExtension<FractionField<IntegerRing>>, u16>>,
> = LazyLock::new(|| FractionField::new(PolynomialRing::<_, u16>::new(Q_I.clone())));

pub static LOGPRINTOPTS: PrintOptions = PrintOptions {
    hide_all_namespaces: false,
    color_namespace: false,
    color_builtin_symbols: false,
    color_top_level_sum: false,
    terms_on_new_line: false,
    print_ring: false,
    include_attributes: false,
    symmetric_representation_for_finite_field: false,
    explicit_rational_polynomial: false,
    number_thousands_separator: None,
    multiplication_operator: '*',
    double_star_for_exponentiation: false,
    square_brackets_for_function: false,
    num_exp_as_superscript: false,
    mode: PrintMode::Symbolica,
    precision: None,
    pretty_matrix: false,
    max_terms: None,
    custom_print_mode: None,
    hide_namespace: Some("gammalooprs"),
};

pub trait IsNeg {
    fn is_negative(&self) -> bool;
}

impl<A: AtomCore> IsNeg for A {
    fn is_negative(&self) -> bool {
        match self.as_atom_view() {
            AtomView::Num(a) => match a.get_coeff_view() {
                CoefficientView::FiniteField(_, _) => false,
                CoefficientView::Float(re, im) => {
                    if im.is_zero() {
                        re.to_float().is_negative()
                    } else {
                        false
                    }
                }
                CoefficientView::Indeterminate => false,
                CoefficientView::Natural(n_re, d_re, n_im, de_im) => {
                    if n_im == 0 {
                        n_re.is_negative()
                    } else {
                        false
                    }
                }
                CoefficientView::Infinity(_) => false,
                CoefficientView::Large(re, im) => {
                    let re = re.to_rat();
                    let im = im.to_rat();
                    if im.is_zero() {
                        re.is_negative()
                    } else {
                        false
                    }
                }
                CoefficientView::RationalPolynomial(_) => false,
            },
            AtomView::Mul(a) => {
                if let Some(first) = a.iter().next() {
                    first.is_negative()
                } else {
                    false
                }
            }
            _ => false,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TypstState {
    pub inproduct: bool,
    pub without_minus: bool,
    pub in_power: bool,
}
pub trait TypstFormat {
    fn typst_string(&self) -> String {
        let mut s = String::new();

        let symbols = self.preable(&mut s).unwrap();
        s.push_str("$");
        self.fmt_output(
            &mut s,
            &symbols,
            TypstState {
                inproduct: false,
                without_minus: false,
                in_power: false,
            },
        )
        .unwrap();
        s.push_str("$");
        s
    }

    fn fmt_output<W: std::fmt::Write>(
        &self,
        f: &mut W,
        symbols: &BTreeMap<Symbol, bool>,
        print_state: TypstState,
    ) -> Result<bool>;

    fn preable<W: std::fmt::Write>(&self, f: &mut W) -> Result<BTreeMap<Symbol, bool>>;
}

impl<A: AtomCore> TypstFormat for A {
    fn preable<W: std::fmt::Write>(&self, fmt: &mut W) -> Result<BTreeMap<Symbol, bool>> {
        let mut symbols = BTreeMap::new();
        self.visitor(&mut |a| match a {
            AtomView::Add(a) => true,
            AtomView::Mul(m) => true,
            AtomView::Fun(f) => {
                *(symbols.entry(f.get_symbol()).or_insert_with(|| true)) = true;
                true
            }
            AtomView::Num(a) => false,
            AtomView::Pow(a) => true,
            AtomView::Var(v) => {
                symbols.entry(v.get_symbol()).or_insert_with(|| false);
                false
            }
        });

        for (s, isfun) in &symbols {
            if let Some(p) = s.get_print_function() {
                if let Some(s) = p(
                    if *isfun {
                        Atom::var(GS.is_function)
                    } else {
                        Atom::var(GS.is_symbol)
                    }
                    .as_view(),
                    &PrintOptions {
                        custom_print_mode: Some(("typst", 2)),
                        ..Default::default()
                    },
                ) {
                    writeln!(fmt, "{}", s).unwrap();
                    continue;
                }
            }

            if *isfun {
                writeln!(
                    fmt,
                    "#let {}-{}(..args)= $op(\"{}\")$",
                    s.get_namespace(),
                    s.get_stripped_name(),
                    s.get_stripped_name()
                )
                .unwrap();
            } else {
                writeln!(
                    fmt,
                    "#let {}-{}= $\"{}\"$",
                    s.get_namespace(),
                    s.get_stripped_name(),
                    s.get_stripped_name()
                )
                .unwrap();
            }
        }

        Ok(symbols)
    }
    fn fmt_output<W: std::fmt::Write>(
        &self,
        fmt: &mut W,
        symbols: &BTreeMap<Symbol, bool>,
        mut print_state: TypstState,
    ) -> Result<bool> {
        match self.as_atom_view() {
            AtomView::Num(n) => match n.get_coeff_view() {
                CoefficientView::FiniteField(e, i) => Ok(true),
                CoefficientView::Float(re, im) => {
                    let re = re.to_float();
                    let im = im.to_float();

                    if im.is_fully_zero() {
                        if print_state.without_minus && re.is_negative() {
                            write!(fmt, "{}", -re)?;
                        } else {
                            write!(fmt, "{}", re)?;
                        }
                    } else {
                        if print_state.inproduct || print_state.in_power {
                            write!(fmt, "(")?;
                        }
                        write!(fmt, "{} + {} i", re, im)?;
                        if print_state.inproduct || print_state.in_power {
                            write!(fmt, ")")?;
                        }
                    }
                    Ok(true)
                }
                CoefficientView::Indeterminate => {
                    write!(fmt, "NA")?;
                    Ok(true)
                }
                CoefficientView::Natural(n_re, d_re, n_im, de_im) => {
                    if d_re == 1 && de_im == 1 {
                        if n_im == 0 {
                            write!(fmt, "{}", n_re)?;
                        } else {
                            if print_state.inproduct || print_state.in_power {
                                write!(fmt, "(")?;
                            }
                            write!(fmt, "{} + {} i", n_re, n_im)?;
                            if print_state.inproduct || print_state.in_power {
                                write!(fmt, ")")?;
                            }
                        }
                    } else {
                        if n_im == 0 {
                            if d_re == 1 {
                                write!(fmt, "{}", n_re)?;
                                return Ok(true);
                            }
                            write!(fmt, "({})/({})", n_re, d_re)?;
                        } else {
                            if de_im == 1 {
                                if print_state.inproduct || print_state.in_power {
                                    write!(fmt, "(")?;
                                }
                                write!(fmt, "({})/({}) + {} i", n_re, d_re, n_im)?;
                                if print_state.inproduct || print_state.in_power {
                                    write!(fmt, ")")?;
                                }
                                return Ok(true);
                            }
                            if print_state.inproduct || print_state.in_power {
                                write!(fmt, "(")?;
                            }
                            write!(fmt, "({})/({})+ ({})/({})i", n_re, d_re, n_im, de_im)?;
                            if print_state.inproduct || print_state.in_power {
                                write!(fmt, ")")?;
                            }
                        }
                    }
                    Ok(true)
                }
                CoefficientView::Infinity(a) => {
                    write!(fmt, "oo")?;
                    Ok(true)
                }
                CoefficientView::Large(re, im) => {
                    let re = re.to_rat();
                    let im = im.to_rat();
                    if im.is_zero() {
                        if re.is_integer() {
                            write!(fmt, "{}", re.numerator())?;
                        } else {
                            write!(fmt, "({})/({})", re.numerator(), re.denominator())?;
                        }
                    } else if re.is_zero() {
                        if im.is_integer() {
                            write!(fmt, "{} i", im.numerator())?;
                        } else {
                            write!(fmt, "({})/({}) i", im.numerator(), im.denominator())?;
                        }
                    } else {
                        if print_state.inproduct || print_state.in_power {
                            write!(fmt, "(")?;
                        }
                        if re.is_integer() && im.is_integer() {
                            write!(fmt, "{} + {} i", re.numerator(), im.numerator())?;
                        } else if re.is_integer() {
                            write!(
                                fmt,
                                "{} + ({})/({}) i",
                                re.numerator(),
                                im.numerator(),
                                im.denominator()
                            )?;
                        } else if im.is_integer() {
                            write!(
                                fmt,
                                "({})/({}) + {} i",
                                re.numerator(),
                                re.denominator(),
                                im.numerator()
                            )?;
                        } else {
                            write!(
                                fmt,
                                "({})/({}) + ({})/({}) i",
                                re.numerator(),
                                re.denominator(),
                                im.numerator(),
                                im.denominator()
                            )?;
                        }
                        if print_state.inproduct || print_state.in_power {
                            write!(fmt, ")")?;
                        }
                    }
                    Ok(true)
                }
                CoefficientView::RationalPolynomial(a) => {
                    write!(fmt, "{}", a.deserialize())?;
                    Ok(true)
                }
            },

            AtomView::Var(v) => {
                write!(
                    fmt,
                    "#{}-{}",
                    v.get_symbol().get_namespace(),
                    v.get_symbol().get_stripped_name()
                )?;
                Ok(true)
            }
            AtomView::Fun(f) => {
                write!(
                    fmt,
                    "#{}-{}(",
                    f.get_symbol().get_namespace(),
                    f.get_symbol().get_stripped_name()
                )?;

                for i in f.iter() {
                    i.fmt_output(
                        fmt,
                        symbols,
                        TypstState {
                            inproduct: false,
                            without_minus: false,
                            in_power: false,
                        },
                    )?;
                    write!(fmt, ", ")?;
                }
                write!(fmt, ")")?;

                Ok(true)
            }
            AtomView::Pow(p) => {
                let (base, exp) = p.get_base_exp();

                if exp.is_negative() {
                    write!(fmt, "1/")?;
                    base.fmt_output(
                        fmt,
                        symbols,
                        TypstState {
                            inproduct: false,
                            without_minus: false,
                            in_power: true,
                        },
                    )?;
                    write!(fmt, "^(")?;
                    exp.fmt_output(
                        fmt,
                        symbols,
                        TypstState {
                            inproduct: false,
                            without_minus: true,
                            in_power: false,
                        },
                    )?;
                    write!(fmt, ")")?;
                } else {
                    base.fmt_output(
                        fmt,
                        symbols,
                        TypstState {
                            inproduct: false,
                            without_minus: false,
                            in_power: true,
                        },
                    )?;
                    write!(fmt, "^(")?;
                    exp.fmt_output(
                        fmt,
                        symbols,
                        TypstState {
                            inproduct: false,
                            without_minus: false,
                            in_power: false,
                        },
                    )?;

                    write!(fmt, ")")?;
                }
                Ok(true)
            }

            AtomView::Mul(t) => {
                print_state.inproduct = true;
                for term in t.iter() {
                    term.fmt_output(fmt, symbols, print_state)?;
                    write!(fmt, " ")?;
                }

                Ok(true)
            }

            AtomView::Add(e) => {
                let mut first = true;
                for term in e.iter() {
                    if term.is_negative() {
                        write!(fmt, " - ")?;
                    } else if !first {
                        write!(fmt, " + ")?;
                    }
                    first = false;
                    print_state.without_minus = true;
                    term.fmt_output(fmt, symbols, print_state)?;
                }
                Ok(true)
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct StringSerializedAtom(pub Atom);

impl JsonSchema for StringSerializedAtom {
    fn schema_name() -> std::borrow::Cow<'static, str> {
        "ParseableAtom".into()
    }
    fn json_schema(_generator: &mut schemars::SchemaGenerator) -> schemars::Schema {
        json_schema!({
            "description": "An atom that is serialized as a string. Do not make any assumptions about the state",
            "type": ["string"]
        })
    }
}

impl Deref for StringSerializedAtom {
    type Target = Atom;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl DerefMut for StringSerializedAtom {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
impl Display for StringSerializedAtom {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Serialize for StringSerializedAtom {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.0.to_canonical_string().serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for StringSerializedAtom {
    fn deserialize<D>(deserializer: D) -> std::result::Result<StringSerializedAtom, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        Ok(StringSerializedAtom(parse!(String::deserialize(
            deserializer
        )?)))
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
        let a = finite_field.to_element(sample_iterator.next().unwrap());
        let b = finite_field.to_element(sample_iterator.next().unwrap());
        Atom::num(1)
        // Atom::i() * Atom::num(b) + Atom::num(a)
    }

    fn prime_generate_ff(sample_iterator: &mut PrimeIteratorU64, finite_field: &Zp64) -> Atom {
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

pub trait DOD {
    fn dod(&self) -> i32;
}

impl DOD for Atom {
    fn dod(&self) -> i32 {
        self.as_view().dod()
    }
}
impl DOD for AtomView<'_> {
    fn dod(&self) -> i32 {
        let rescaled = self
            .replace(GS.emr_mom.f(&[W_.a__]))
            .with(GS.emr_mom.f(&[W_.a__]) * GS.rescale);

        let series = rescaled
            .series(GS.rescale, Atom::Zero, (4, 1).into(), true)
            .unwrap();

        let dod = series.degree();

        if dod.is_integer() {
            dod.numerator().to_i64().unwrap() as i32
        } else {
            panic!("{dod} for {self}")
        }
    }
}

pub trait CallSymbol<T> {
    fn f(&self, args: T) -> Atom;
}

impl<'a, I> CallSymbol<&'a [I]> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
{
    fn f(&self, args: &'a [I]) -> Atom {
        FunctionBuilder::new(*self).add_args(args).finish()
    }
}
impl<'a, I, const N: usize> CallSymbol<&'a [I; N]> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
{
    fn f(&self, args: &'a [I; N]) -> Atom {
        FunctionBuilder::new(*self).add_args(args).finish()
    }
}

impl<'a, 'b, I, J> CallSymbol<(&'a [I], &'b [J])> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
    &'b J: Into<AtomOrView<'b>>,
{
    fn f(&self, args: (&'a [I], &'b [J])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0)
            .add_args(args.1)
            .finish()
    }
}

impl<'a, 'b, I, J, const N: usize> CallSymbol<(&'a [I; N], &'b [J])> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
    &'b J: Into<AtomOrView<'b>>,
{
    fn f(&self, args: (&'a [I; N], &'b [J])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0.as_slice())
            .add_args(args.1)
            .finish()
    }
}

impl<'a, 'b, I, J, const N: usize, const M: usize> CallSymbol<(&'a [I; N], &'b [J; M])> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
    &'b J: Into<AtomOrView<'b>>,
{
    fn f(&self, args: (&'a [I; N], &'b [J; M])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0.as_slice())
            .add_args(args.1)
            .finish()
    }
}
impl<'a, 'b, 'c, I, J, K, const N: usize, const M: usize, const O: usize>
    CallSymbol<(&'a [I; N], &'b [J; M], &'c [K; O])> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
    &'b J: Into<AtomOrView<'b>>,
    &'c K: Into<AtomOrView<'c>>,
{
    fn f(&self, args: (&'a [I; N], &'b [J; M], &'c [K; O])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0.as_slice())
            .add_args(args.1)
            .add_args(args.2)
            .finish()
    }
}

impl<I, const N: usize> CallSymbol<[I; N]> for Symbol
where
    for<'a> &'a I: Into<AtomOrView<'a>>,
{
    fn f(&self, args: [I; N]) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.as_slice())
            .finish()
    }
}

impl<'b, I, J, const N: usize> CallSymbol<([I; N], &'b [J])> for Symbol
where
    for<'a> &'a I: Into<AtomOrView<'a>>,
    &'b J: Into<AtomOrView<'b>>,
{
    fn f(&self, args: ([I; N], &'b [J])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0.as_slice())
            .add_args(args.1)
            .finish()
    }
}

impl<'a, I, J, K> CallSymbol<(&'a [I], &'a [J], &'a [K])> for Symbol
where
    &'a I: Into<AtomOrView<'a>>,
    &'a J: Into<AtomOrView<'a>>,
    &'a K: Into<AtomOrView<'a>>,
{
    fn f(&self, args: (&'a [I], &'a [J], &'a [K])) -> Atom {
        FunctionBuilder::new(*self)
            .add_args(args.0)
            .add_args(args.1)
            .add_args(args.2)
            .finish()
    }
}
