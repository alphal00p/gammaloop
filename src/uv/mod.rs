use std::hash::Hash;

use crate::{numerator::aind::Aind, utils::GS};
use spenso::{
    network::parsing::ShadowedStructure,
    structure::{
        NamedStructure, ToSymbolic,
        dimension::Dimension,
        representation::{Minkowski, RepName},
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

#[allow(dead_code)]
pub(crate) fn is_not_paired(pair: &HedgePair) -> bool {
    !pair.is_paired()
}

pub mod hedge_poset;
pub mod settings;
pub use settings::UVgenerationSettings;
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

pub mod profile;
pub use profile::{
    UVProfileConfig, UVProfileResult, display_results_summary, run_uv_profile, write_results,
};

#[cfg(test)]
mod tests;
