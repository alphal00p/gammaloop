//! Graph-first R* preparation with the IR counterterms expressed as UV maps.
//!
//! The plan preserves the separate counterterm scopes in the unfolded form of
//! Chetyrkin's Eq. (4.49). Taylor operations act on explicit expressions only;
//! the remaining single-mass integrals are evaluated through Vakint p-integrals.
//! Hard vacuum tensors are projected onto local covariants before the quotient
//! UV maps. Their scalar hard integrals stay outside those maps, which reuse
//! the ordinary massive R operation.

mod evaluate;
mod expression;
mod integrate;
mod massive;
#[cfg(test)]
mod massive_tests;
#[cfg(test)]
mod master_tests;
#[cfg(test)]
mod matter_tests;
mod moments;
mod plan;
#[cfg(test)]
mod pole_tests;
mod project;
#[cfg(test)]
mod scalar;
mod taylor;
#[cfg(test)]
mod tests;
mod uv;

pub use integrate::integrate_single_mass;
pub use massive::{
    CountertermCache, DistributionTerm, MassiveVacuumSubdivergences, MassiveVacuumTerm,
    NumeratorMonomial, RecursionContext, VacuumTopology, massive_vacuum_direct,
    massive_vacuum_settings,
};
pub use plan::*;
pub use taylor::TaylorExpansion;
pub use uv::MassiveUvEvaluator;
