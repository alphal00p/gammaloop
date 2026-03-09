use spenso::{
    network::{library::symbolic::ETS, parsing::SPENSO_TAG},
    structure::{
        abstract_index::AIND_SYMBOLS,
        representation::{Euclidean, Lorentz, Minkowski, RepName},
    },
};
use spenso_macros::SimpleRepresentation;
use symbolica::atom::Atom;

use super::{color::CS, gamma::AGS, metric::MS, rep_symbols::RS};

#[cfg(feature = "python")]
use pyo3::pyfunction;
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::gen_stub_pyfunction;
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
)]
#[cfg_attr(
    feature = "bincode",
    derive(bincode_trait_derive::Encode),
    derive(bincode_trait_derive::Decode),
    derive(bincode_trait_derive::BorrowDecodeFromDecode),
)]
#[representation(name = "spf", dual_name = "SpinAntiFundamental")] // Specify the dual name
pub struct SpinFundamental {}

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
)]
#[cfg_attr(
    feature = "bincode",
    derive(bincode_trait_derive::Encode),
    derive(bincode_trait_derive::Decode),
    derive(bincode_trait_derive::BorrowDecodeFromDecode),
)]
#[representation(name = "cof", dual_name = "ColorAntiFundamental")] // Specify the dual name
pub struct ColorFundamental {}

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
)]
#[cfg_attr(
    feature = "bincode",
    derive(bincode_trait_derive::Encode),
    derive(bincode_trait_derive::Decode),
    derive(bincode_trait_derive::BorrowDecodeFromDecode),
)]
#[representation(name = "cos", dual_name = "ColorAntiSextet")] // Specify the dual name
pub struct ColorSextet {}

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
)]
#[cfg_attr(
    feature = "bincode",
    derive(bincode_trait_derive::Encode),
    derive(bincode_trait_derive::Decode),
    derive(bincode_trait_derive::BorrowDecodeFromDecode),
)]
#[representation(name = "bis", self_dual)] // Specify the dual name
pub struct Bispinor {}

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
)]
#[cfg_attr(
    feature = "bincode",
    derive(bincode_trait_derive::Encode),
    derive(bincode_trait_derive::Decode),
    derive(bincode_trait_derive::BorrowDecodeFromDecode),
)]
#[representation(name = "coad", self_dual)] // Specify the dual name
pub struct ColorAdjoint {}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.spenso")
)]
#[cfg_attr(feature = "python", pyfunction)]
pub fn initialize() {
    let _ = AIND_SYMBOLS.dind;
    let _ = Minkowski {}.to_symbolic([Atom::Zero]);
    let _ = Euclidean {}.to_symbolic([Atom::Zero]);
    let _ = Lorentz {}.to_symbolic([Atom::Zero]);
    let _ = SpinFundamental {}.to_symbolic([Atom::Zero]);
    let _ = Bispinor {}.to_symbolic([Atom::Zero]);
    let _ = ColorAdjoint {}.to_symbolic([Atom::Zero]);
    let _ = ColorFundamental {}.to_symbolic([Atom::Zero]);
    let _ = ColorSextet {}.to_symbolic([Atom::Zero]);
    let _ = RS.a_;
    let _ = MS.dummy;
    let _ = AGS.gamma;
    let _ = ETS.metric;
    let _ = CS.f;
    let _ = SPENSO_TAG.bracket;
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_spin_fundamental() {
        let spin_fundamental = SpinFundamental::default();

        // let dual = spin_fundamental.dual();

        assert_eq!(spin_fundamental.dual(), SpinAntiFundamental {});
        // assert_eq!(spin_fundamental.dual_name(), "SpinAntiFundamental");
    }
}
