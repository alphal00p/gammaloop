use bincode_trait_derive::{Decode, Encode};
use schemars::{json_schema, JsonSchema};
use serde::{Deserialize, Serialize};
use std::{
    fmt::Display,
    num,
    ops::{Deref, DerefMut},
    sync::LazyLock,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    domains::{
        algebraic_number::AlgebraicExtension,
        finite_field::PrimeIteratorU64,
        float::Complex,
        integer::IntegerRing,
        rational::{FractionField, Rational, Q},
    },
    parse,
    poly::polynomial::PolynomialRing,
    printer::{PrintMode, PrintOptions},
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
    print_finite_field: true,
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
