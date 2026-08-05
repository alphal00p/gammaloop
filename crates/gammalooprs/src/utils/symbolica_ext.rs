use std::{
    fmt::Display,
    ops::{Deref, DerefMut},
    sync::LazyLock,
};

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::involution::EdgeIndex;
use schemars::{JsonSchema, json_schema};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    domains::{
        algebraic_number::AlgebraicExtension,
        integer::IntegerRing,
        rational::{FractionField, Q},
    },
    function, parse,
    poly::polynomial::PolynomialRing,
    printer::{PrintMode, PrintOptions},
    symbol,
};

use crate::GammaLoopContext;

use super::{GS, W_};

pub static Q_I: LazyLock<AlgebraicExtension<FractionField<IntegerRing>>> =
    LazyLock::new(|| AlgebraicExtension::new_complex(Q));
static RAW_UFO_MOMENTUM: LazyLock<Symbol> = LazyLock::new(|| symbol!("UFO::P"));
static RAW_UFO_PSLASH: LazyLock<Symbol> = LazyLock::new(|| symbol!("UFO::PSlash"));

pub static COMPLEXRATPOLYFIELD: LazyLock<
    FractionField<PolynomialRing<AlgebraicExtension<FractionField<IntegerRing>>, u16>>,
> = LazyLock::new(|| FractionField::new(PolynomialRing::<_, u16>::new(Q_I.clone())));

pub static LOGPRINTOPTS: LazyLock<PrintOptions> = LazyLock::new(|| PrintOptions {
    hide_all_namespaces: true,
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
    num_exp_as_superscript: false,
    mode: PrintMode::Symbolica,
    precision: None,
    pretty_matrix: false,
    max_terms: None,
    custom_print_mode: Default::default(),
    hide_namespace: Some(std::borrow::Cow::Borrowed("gammalooprs")),
    ..PrintOptions::new()
});

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

pub trait DOD {
    /// Rescales momentum of edge `eid`, and computes the leading scaling.
    fn edge_dod(&self, eid: EdgeIndex) -> i32;

    /// Rescales all momenta, and computes the leading scaling.
    fn all_dod(&self) -> i32;

    fn trailing_exponent(&self) -> i32;
}

impl DOD for Atom {
    fn edge_dod(&self, eid: EdgeIndex) -> i32 {
        self.as_view().edge_dod(eid)
    }

    fn all_dod(&self) -> i32 {
        self.as_view().all_dod()
    }

    fn trailing_exponent(&self) -> i32 {
        self.as_view().trailing_exponent()
    }
}

impl DOD for AtomView<'_> {
    fn edge_dod(&self, eid: EdgeIndex) -> i32 {
        self.replace(GS.emr_mom(eid, W_.a___))
            .with(GS.emr_mom(eid, W_.a___) / GS.rescale)
            .trailing_exponent()
    }

    fn all_dod(&self) -> i32 {
        self.replace(function!(GS.emr_mom, W_.a___))
            .with(function!(GS.emr_mom, W_.a___) / GS.rescale)
            .replace(function!(*RAW_UFO_MOMENTUM, W_.a___))
            .with(function!(*RAW_UFO_MOMENTUM, W_.a___) / GS.rescale)
            .replace(function!(*RAW_UFO_PSLASH, W_.a___))
            .with(function!(*RAW_UFO_PSLASH, W_.a___) / GS.rescale)
            .trailing_exponent()
    }

    fn trailing_exponent(&self) -> i32 {
        let series = self.series(GS.rescale, Atom::Zero, 1).unwrap();
        let dod = series.get_trailing_exponent();

        if dod.is_integer() {
            -(dod.numerator().to_i64().unwrap() as i32)
        } else {
            panic!("{dod} for {self}")
        }
    }
}

#[test]
fn test_dod() {
    let (e1, e2) = (EdgeIndex(1), EdgeIndex(2));

    let atom = (GS.emr_mom(e1, Atom::Zero) * GS.emr_mom(e2, Atom::Zero)
        + GS.emr_mom(e2, Atom::Zero))
        / (GS.emr_mom(e1, Atom::Zero) * GS.emr_mom(e1, Atom::Zero));

    let atom2 =
        Atom::num(1) / (GS.emr_mom(e1, Atom::Zero) * GS.emr_mom(e1, Atom::Zero) + parse!("m"));
    let atom3 =
        GS.emr_mom(e1, Atom::Zero) * GS.emr_mom(e2, Atom::Zero) + GS.emr_mom(e2, Atom::Zero);

    assert_eq!(-1, atom.edge_dod(e1));
    assert_eq!(1, atom.edge_dod(e2));
    assert_eq!(-2, atom2.edge_dod(e1));
    assert_eq!(1, atom3.edge_dod(e1));
    assert_eq!(1, atom3.edge_dod(e2));
}
