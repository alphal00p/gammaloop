use std::hash::Hash;

use crate::{numerator::aind::Aind, utils::GS};
use spenso::{
    network::parsing::ShadowedStructure,
    structure::{
        dimension::Dimension,
        representation::{Minkowski, RepName},
        NamedStructure, ToSymbolic,
    },
};
use symbolica::atom::Atom;

use linnet::half_edge::involution::HedgePair;

// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

pub(crate) fn spenso_lor(
    tag: i32,
    ind: impl Into<Aind>,
    dim: impl Into<Dimension>,
) -> ShadowedStructure<Aind> {
    let mink = Minkowski {}.new_slot(dim, ind);
    NamedStructure::from_iter([mink], GS.emr_mom, Some(vec![Atom::num(tag)])).structure
}

pub(crate) fn spenso_lor_atom(tag: i32, ind: impl Into<Aind>, dim: impl Into<Dimension>) -> Atom {
    spenso_lor(tag, ind, dim).to_symbolic(None).unwrap()
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct IntegrandExpr {
    integrand: Atom,
    // add_arg: Option<Atom>,
}

pub(crate) fn is_not_paired(pair: &HedgePair) -> bool {
    !pair.is_paired()
}

pub mod uv_graph;
pub use uv_graph::UltravioletGraph;

pub mod poset;
pub use poset::Poset;

pub mod wood;
pub use wood::Wood;

pub mod approx;
pub use approx::ApproxOp;

pub mod forest;
pub use forest::Forest;

#[cfg(test)]
mod tests;
