//! Graph-first R* preparation with the IR counterterms expressed as UV maps.
//!
//! The plan preserves the separate counterterm scopes in the unfolded form of
//! Chetyrkin's Eq. (4.49). Taylor operations act on explicit expressions only;
//! the remaining single-mass integrals are evaluated through Vakint p-integrals.
//! Hard vacuum tensors are projected onto local covariants before the quotient
//! UV maps. Their scalar hard integrals stay outside those maps, which reuse
//! the ordinary massive R operation.

mod boundary;
#[cfg(test)]
mod boundary_tests;
#[cfg(test)]
mod cache_tests;
mod closed;
#[cfg(test)]
mod closed_skip_tests;
mod edge_experiment;
mod evaluate;
mod expression;
mod insertion_identity;
mod insertions;
mod integrate;
mod ir_dirac;
#[cfg(test)]
mod ir_dirac_tests;
mod linearize;
mod massive;
#[cfg(test)]
mod massive_regression_tests;
#[cfg(test)]
mod massive_tests;
#[cfg(test)]
mod master_tests;
#[cfg(test)]
mod matter_tests;
mod moments;
mod physical;
pub mod pint6_dump;
mod plan;
#[cfg(test)]
mod pole_tests;
mod routing;
mod routing_num;

/// [`linearize::linearize`] (`GAMMALOOP_LINEARIZE_MOMENTA`) for the converters
/// outside R*.
pub(crate) fn linearize_momenta(atom: symbolica::atom::Atom) -> symbolica::atom::Atom {
    linearize::linearize(atom)
}
mod project;
mod propagator;
#[cfg(test)]
mod recursive_tests;
#[cfg(test)]
mod scalar;
mod taylor;
#[cfg(test)]
mod tests;
mod top_mass;
mod uncontracted;
#[cfg(test)]
mod uncontracted_init_refusal_tests;
mod uncontracted_ledger;
#[cfg(test)]
mod uncontracted_ledger_tests;
#[cfg(test)]
mod uncontracted_oracle_k2_tests;
#[cfg(test)]
mod uncontracted_oracle_tests;
#[cfg(test)]
mod uncontracted_prereview_tests;
#[cfg(test)]
mod uncontracted_refusal_tests;
#[cfg(test)]
mod uncontracted_review_tests;
#[cfg(test)]
mod uncontracted_tests;
mod uv;

pub use crate::uv::counterterm_cache::{
    CacheTimingStats, CountertermCache, CountertermCacheStats, CountertermPrewarmPlan,
    CountertermRequest,
};

pub(super) static REFERENCE_MASS_SQUARED: std::sync::LazyLock<symbolica::atom::Symbol> =
    std::sync::LazyLock::new(|| symbolica::symbol!("rstar_reference_mass_squared"; Scalar));

/// A thread-local `RefCell<Option<T>>` set for one scope and restored when the
/// scope ends, also when it unwinds: beta catches a panicking diagram (an
/// exhausted budget inside idenso panics) and its worker thread goes on with
/// the next diagram, which must not inherit the scope.
struct ScopedSlot<T: 'static> {
    key: &'static std::thread::LocalKey<std::cell::RefCell<Option<T>>>,
    previous: Option<Option<T>>,
}

impl<T: 'static> ScopedSlot<T> {
    /// Set the slot to `value` until the returned scope ends.
    fn enter(key: &'static std::thread::LocalKey<std::cell::RefCell<Option<T>>>, value: T) -> Self {
        let previous = key.with(|slot| slot.replace(Some(value)));
        Self {
            key,
            previous: Some(previous),
        }
    }

    /// End the scope normally: the previous value is back, and this scope's
    /// value is returned.
    fn exit(mut self) -> Option<T> {
        let previous = self.previous.take().flatten();
        self.key.with(|slot| slot.replace(previous))
    }
}

impl<T: 'static> Drop for ScopedSlot<T> {
    fn drop(&mut self) {
        // Only when the scope did not end through `exit` (unwinding). Inner
        // borrows are released before this runs; never panic here.
        if let Some(previous) = self.previous.take() {
            let _ = self.key.try_with(|slot| {
                if let Ok(mut slot) = slot.try_borrow_mut() {
                    *slot = previous;
                }
            });
        }
    }
}
pub use boundary::{
    BoundaryMode, ChildBoundary, ClosedRootNormalization, child_boundary_mode,
    remainder_boundary_mode, remainder_boundary_provenance, startup_line,
    validate_remainder_boundary,
};
pub use expression::{ClosedSkipMode, closed_dot_scalar, closed_skip_mode, set_memory_release};
pub use insertion_identity::{InsertionSkeletonIdentity, InsertionSkeletonKey};
pub use insertions::{
    GroupedInsertionPoleResult, InsertionPoleResult, PreparedInsertion, PreparedInsertionHole,
    PropagatorInsertionEvaluator, surrogate_closed_root,
};
pub use integrate::integrate_single_mass;
pub use ir_dirac::{IrDiracMode, ir_dirac_mode};
pub(crate) use ir_dirac::{diagnostic_failure, is_diagnostic_failure};
pub use massive::{
    DistributionTerm, MassiveVacuumSubdivergences, MassiveVacuumTerm, NumeratorMonomial,
    VacuumTopology, massive_vacuum_direct, massive_vacuum_settings,
};
pub use physical::{PhysicalUvGroup, PhysicalUvKey, PreparedPhysicalUvTerm};
pub use plan::*;
pub use propagator::integrate_propagator;
pub use taylor::TaylorExpansion;
pub use top_mass::{
    RStarTerms, TopMassTerms, auxiliary_mass_squared, rstar_terms, split_auxiliary_mass,
    split_enabled as top_mass_split_enabled, top_level_switches_line,
    top_level_switches_provenance, top_mass_terms, validate_top_level_switches,
};
pub use uncontracted_ledger::{DiagramLedger, Resumed, RunTally, resumed as uncontracted_resumed};
pub use uv::MassiveUvEvaluator;

/// `GAMMALOOP_UNCONTRACTED_NUMERATOR`'s initialisation (call at start-up for
/// an early error): a bad value, `GAMMALOOP_UNCONTRACTED_STRICT=0` (on and
/// check), a bad `GAMMALOOP_UNCONTRACTED_TIMING` (every mode), check mode
/// with `PINT6_CLOSURE_DUMP`, a kind list without a kind this build emits,
/// an enabled kind whose route is off (`child` with
/// `GAMMALOOP_UV_CHILDREN=massive`, `ir-hard` off the IR hard child route),
/// a quotient route other than `GAMMALOOP_IR_QUOTIENT=rstar` (on and check;
/// helper VH's C7 settings exception E2) or a bad exception list
/// (`GAMMALOOP_UNCONTRACTED_EXCEPTIONS`) is refused.
pub fn validate_uncontracted() -> eyre::Result<()> {
    uncontracted::validate()
}

/// Whether a diagram ledger opens ([`DiagramLedger::open`] is `Some`): the
/// switch on or check, or `GAMMALOOP_UNCONTRACTED_TIMING=1` (off included).
pub fn uncontracted_ledger_active() -> eyre::Result<bool> {
    uncontracted::validate()?;
    Ok(uncontracted_ledger::active())
}

/// The identity of `GAMMALOOP_UNCONTRACTED_NUMERATOR` for a checkpoint
/// (U7a delta, U7a.1 (i)): `None` off (an off checkpoint's identity is
/// unchanged); else the mode, the enabled kinds, strictness and the
/// exception list's digest (`check:root,child,ir-hard:strict`).
pub fn uncontracted_identity() -> eyre::Result<Option<String>> {
    uncontracted::validate()?;
    Ok(uncontracted_ledger::current_identity())
}

/// The run-level tally of `GAMMALOOP_UNCONTRACTED_NUMERATOR` on / check
/// ([`RunTally`]; `None` off): the beta adds every diagram record of a loop
/// order (computed and reused), prints its totals on `UNCONTRACTED_RUN` and
/// fails the run on its failures (pre-review of the U7a delta, finding 6).
pub fn uncontracted_run_tally() -> eyre::Result<Option<RunTally>> {
    uncontracted::validate()?;
    Ok(RunTally::new())
}
