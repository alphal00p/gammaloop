//! Shared tensor heads for edge momenta and their stored routing coordinates.

use spenso::{
    structure::representation::{Minkowski, RepName},
    symbolica_init::in_symbolica_initializer,
};
use symbolica::{
    atom::{Atom, Symbol},
    initialize, symbol,
};

initialize!(|| in_symbolica_initializer(|| {
    let _ = Minkowski {}.to_symbolic([Atom::Zero]);
    let _ = dimension();
    let _ = momentum();
    let _ = loop_momentum();
    let _ = external_momentum();
}));

/// Symbolic Lorentz dimension used by graph-derived Python integral families.
pub fn dimension() -> Symbol {
    symbol!("gammalooprs::dim")
}

/// Edge-indexed momentum, matching the diagram generator's numerator convention.
pub fn momentum() -> Symbol {
    symbol!(
        "FeynKit::Momentum",
        tags = ["spenso::tensor", "spenso::rank1"]
    )
}

/// Independent loop coordinate, indexed by its position in the stored basis.
pub fn loop_momentum() -> Symbol {
    symbol!(
        "FeynKit::LoopMomentum",
        tags = ["spenso::tensor", "spenso::rank1"]
    )
}

/// External coordinate, indexed before elimination of dependent external columns.
pub fn external_momentum() -> Symbol {
    symbol!(
        "FeynKit::ExternalMomentum",
        tags = ["spenso::tensor", "spenso::rank1"]
    )
}
