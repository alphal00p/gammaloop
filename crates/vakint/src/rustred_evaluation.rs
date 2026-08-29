use std::fmt;

use symbolica::atom::{Atom, AtomView};
use thiserror::Error;

use crate::{ReplacementRules, Vakint, VakintError, VakintSettings};

/// Reserved options for the FORM-independent RustRed scalar-integral evaluator.
///
/// The public boundary is available so callers can prepare their evaluation stacks,
/// but it does not claim support for any topology until closing IBP artifacts land.
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
pub enum RustRedEvaluationError {
    #[error(
        "RustRed scalar reduction is not yet available for the matched {loop_count}-loop topology"
    )]
    ReducerUnavailable { loop_count: usize },
}

impl Vakint {
    /// Reserved scalar-evaluation seam for RustRed.
    ///
    /// `integral_specs` is the routing and topology witness already produced by Vakint's
    /// matcher. Once closing artifacts are available, the reducer must consume this
    /// witness directly rather than matching the graph a second time. Until then, direct
    /// dispatch returns [`RustRedEvaluationError::ReducerUnavailable`].
    pub(crate) fn rustred_evaluate(
        &self,
        _settings: &VakintSettings,
        _numerator: AtomView,
        integral_specs: &ReplacementRules,
        _options: &RustRedEvaluationOptions,
    ) -> Result<Atom, VakintError> {
        Err(RustRedEvaluationError::ReducerUnavailable {
            loop_count: integral_specs.canonical_topology.get_integral().n_loops,
        }
        .into())
    }
}
