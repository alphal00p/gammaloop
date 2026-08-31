//! FORM-independent scalar reduction through sealed RustRed artifacts.

mod artifact;
mod matching;
mod materialize;
mod terminal;

use std::fmt;

use symbolica::atom::{Atom, AtomView};
use thiserror::Error;

use crate::{ReplacementRules, Topology, Vakint, VakintError, VakintSettings};

use matching::MatchedScalarFamily;

/// Options for the FORM-independent RustRed scalar-integral evaluator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RustRedEvaluationOptions {
    /// Replace RustRed's terminal master integrals by Vakint's known evaluations.
    pub substitute_masters: bool,
}

impl Default for RustRedEvaluationOptions {
    fn default() -> Self {
        Self {
            substitute_masters: true,
        }
    }
}

impl fmt::Display for RustRedEvaluationOptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "substitute_masters={}", self.substitute_masters)
    }
}

/// Errors owned by the RustRed scalar-evaluation adapter.
#[derive(Debug, Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum RustRedEvaluationError {
    #[deprecated(
        since = "0.1.2",
        note = "retained for source compatibility with the pre-artifact RustRed seam"
    )]
    #[error(
        "RustRed scalar reduction is not yet available for the matched {loop_count}-loop topology"
    )]
    ReducerUnavailable { loop_count: usize },
    #[error("RustRed has no sealed scalar artifact for this matched family: {detail}")]
    UnsupportedMatchedFamily { detail: String },
    #[error("RustRed could not read propagator {propagator} power as an integer: {power}")]
    InvalidPower { propagator: usize, power: String },
    #[error("RustRed rejected the matched single-scale vacuum family: {detail}")]
    InvalidMatchedFamily { detail: String },
    #[error("RustRed could not load the shipped {family} artifact: {detail}")]
    ArtifactLoad {
        family: &'static str,
        detail: String,
    },
    #[error("RustRed could not construct the integral key: {detail}")]
    IntegralKey { detail: String },
    #[error("RustRed could not reduce the matched integral: {detail}")]
    Reduction { detail: String },
    #[error("RustRed could not lower the scalar numerator: {detail}")]
    ScalarNumerator { detail: String },
    #[error("RustRed returned an unsupported {family} terminal master {powers:?}")]
    UnsupportedMaster {
        family: &'static str,
        powers: Vec<i64>,
    },
    #[error("RustRed common-mass exponent {exponent} does not fit Vakint's exact exponent type")]
    MassExponentOverflow { exponent: i128 },
    #[error(
        "RustRed common-mass exponent overflowed while adding reduction power {reduction} and scalar-numerator power {numerator}"
    )]
    MassExponentAdditionOverflow { reduction: i128, numerator: u32 },
}

/// Whether the opt-in scalar backend owns a sealed artifact for this matcher class.
pub(crate) fn supports(
    settings: &VakintSettings,
    topology: &Topology,
    options: &RustRedEvaluationOptions,
) -> bool {
    (!options.substitute_masters || settings.number_of_terms_in_epsilon_expansion <= 5)
        && MatchedScalarFamily::try_from_topology(topology).is_ok()
}

impl Vakint {
    /// Reduce one already-matched and simultaneously routed scalar integral.
    ///
    /// `integral_specs` is the sole topology/routing witness. This adapter does
    /// not rematch a graph and never dispatches on a topology name.
    pub(crate) fn rustred_evaluate(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        integral_specs: &ReplacementRules,
        options: &RustRedEvaluationOptions,
    ) -> Result<Atom, VakintError> {
        MatchedScalarFamily::try_new(integral_specs)?
            .evaluate(settings, numerator, options)
            .map_err(Into::into)
    }
}
