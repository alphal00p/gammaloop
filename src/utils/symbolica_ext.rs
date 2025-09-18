use bincode_trait_derive::{Decode, Encode};
use schemars::{json_schema, JsonSchema};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::{
        algebraic_number::AlgebraicExtension,
        integer::IntegerRing,
        rational::{FractionField, Q},
    },
    parse,
    poly::polynomial::PolynomialRing,
    printer::{PrintMode, PrintOptions},
};

use std::{
    fmt::Display,
    ops::{Deref, DerefMut},
    sync::LazyLock,
};

use crate::GammaLoopContext;

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
