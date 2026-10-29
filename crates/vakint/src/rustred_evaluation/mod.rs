//! FORM-independent scalar reduction through shipped RustRed programs.
//!
//! One through three loops use sealed artifacts. Four loops use finite-tested
//! candidate programs with pointwise guard, descent, and terminal checks;
//! successful evaluation is not an unrestricted family-closure claim.

mod artifact;
pub(crate) mod data;
pub(crate) mod masters;
pub(crate) use experimental::timing::MathSpan;
pub(crate) use masters::combine_exact;
#[cfg(all(unix, not(target_arch = "wasm32")))]
mod official;
#[cfg(all(unix, not(target_arch = "wasm32")))]
pub use official::{ExternalDataLoadReport, ExternalDataLoopReport, OfficialExternalDataOptions};
pub mod cache;
#[cfg(test)]
mod cache_four_loop_tests;
#[cfg(test)]
mod cache_lower_tests;
#[cfg(test)]
mod cache_tests;
mod factor_projection;
mod factor_provider;
pub(crate) mod five_loop;
mod four_loop;
mod mass_channels;
mod matching;
mod materialize;
pub(crate) mod monitor;
mod numerator;
mod ordered;
#[cfg(test)]
mod terminal;

#[cfg(feature = "experimental-rustred")]
pub mod experimental;
#[cfg(not(feature = "experimental-rustred"))]
// Ordinary 4L evaluation shares the implementation; optional public diagnostics
// and offline preparation APIs remain unused when this module is private.
#[allow(dead_code)]
mod experimental;

use monitor::MonitorContext;
use rustred::solver::CandidateRoutedApplicationObserver;
use std::fmt;

use symbolica::atom::{Atom, AtomView};
use thiserror::Error;

use crate::{ReplacementRules, Topology, Vakint, VakintError, VakintSettings};

use matching::MatchedScalarFamily;

/// One contribution after the shared topology witness and numerator routing.
/// Ordered reduction batches all contributions before expanding any integral.
pub(crate) struct MatchedContribution {
    pub(crate) numerator: Atom,
    pub(crate) matched: ReplacementRules,
}

/// How the final RustRed masters are represented.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord)]
pub enum RustRedMasterMode {
    Retain,
    Symbolic,
    #[default]
    Numerical,
}

impl fmt::Display for RustRedMasterMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Retain => "retain",
            Self::Symbolic => "symbolic",
            Self::Numerical => "numerical",
        })
    }
}

/// Options for the FORM-independent RustRed scalar-integral evaluator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RustRedEvaluationOptions {
    pub master_mode: RustRedMasterMode,
}

impl Default for RustRedEvaluationOptions {
    fn default() -> Self {
        Self {
            master_mode: RustRedMasterMode::Numerical,
        }
    }
}

impl fmt::Display for RustRedEvaluationOptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "master_mode={}", self.master_mode)
    }
}

/// Errors owned by the RustRed scalar-evaluation adapter.
#[derive(Debug, Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum RustRedEvaluationError {
    #[error("RustRed resource limit: {0}")]
    ResourceLimit(#[from] symbolica::utils::resources::ResourceLimitError),
    #[error("RustRed external data: {detail}")]
    ExternalData { detail: String },
    #[error("RustRed master data capability: {detail}")]
    MasterCapability { detail: String },
    #[error("RustRed cache backend is unavailable for {loop_count} loops")]
    CacheBackendUnavailable { loop_count: usize },
    #[error("RustRed lower-loop cache lock is poisoned")]
    CacheLockPoisoned,
    #[error("RustRed cache operation failed: {source}")]
    CacheOperation { source: rustred_app::AppError },
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
    #[error(
        "RustRed has no shipped scalar artifact for the {loop_count}-loop family yet: {detail}"
    )]
    FourLoopArtifactUnavailable { loop_count: usize, detail: String },
    #[error("five-loop RustRed evaluation requires a loaded saved application package")]
    SavedApplicationUnavailable,
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

impl RustRedEvaluationError {
    pub(super) fn from_native(error: impl std::error::Error + 'static) -> Self {
        match Self::resource_cause(&error) {
            Some(resource) => Self::ResourceLimit(resource.clone()),
            None => Self::Reduction {
                detail: error.to_string(),
            },
        }
    }

    pub(super) fn resource_cause<'a>(
        error: &'a (dyn std::error::Error + 'static),
    ) -> Option<&'a symbolica::utils::resources::ResourceLimitError> {
        let mut source = Some(error);
        while let Some(error) = source {
            if let Some(error) =
                error.downcast_ref::<symbolica::utils::resources::ResourceLimitError>()
            {
                return Some(error);
            }
            // io::Error's source() delegates to the wrapped error's source(),
            // which can skip a ResourceLimitError leaf entirely.
            source = match error
                .downcast_ref::<std::io::Error>()
                .and_then(std::io::Error::get_ref)
            {
                Some(inner) => Some(inner as &(dyn std::error::Error + 'static)),
                None => error.source(),
            };
        }
        None
    }
}

/// Whether the opt-in scalar backend owns a shipped program for this matcher class.
pub(crate) fn supports(
    settings: &VakintSettings,
    topology: &Topology,
    options: &RustRedEvaluationOptions,
) -> bool {
    let _ = (settings, options);
    MatchedScalarFamily::try_from_topology(topology).is_ok()
        || four_loop::supports(topology)
        || five_loop::supports(topology)
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
        let raw = self.rustred_reduce(settings, numerator, integral_specs)?;
        self.finalize_rustred(
            raw,
            settings,
            integral_specs.canonical_topology.get_integral().n_loops,
            options.master_mode,
            None,
        )
    }

    pub(crate) fn finalize_rustred(
        &self,
        raw: Atom,
        settings: &VakintSettings,
        loop_count: usize,
        mode: RustRedMasterMode,
        observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    ) -> Result<Atom, VakintError> {
        // Ordinary reduction simplifies each normalized contribution before
        // this boundary. Ordered reduction combines those contributions first,
        // so apply the same convention once here, before master substitution.
        let raw = if settings.rustred_application_options.reduction_strategy
            == crate::RustRedReductionStrategy::Ordered
        {
            crate::utils::simplify_real(raw.as_view())
        } else {
            raw
        };
        let data = self.rustred_data()?;
        let announced = self.announced_bare_masters()?;
        let result = if mode == RustRedMasterMode::Retain {
            data.masters.retain(raw, &announced, observation)
        } else {
            if !announced.is_empty() {
                data.masters
                    .require_coefficients(&combine_exact(raw.clone()), &announced)?;
            }
            data.masters.finalize(raw, settings, loop_count, mode)?
        };
        if settings.use_dot_product_notation {
            Ok(result)
        } else {
            let indexed = Vakint::convert_from_dot_notation(result.as_view());
            Ok(crate::VakintTerm::repair_expanded_spenso_vector_slots(
                indexed.as_view(),
            ))
        }
    }

    /// Reduce and normalize without inspecting coefficient availability. The
    /// expression owner combines all contributions before materializing masters.
    pub(crate) fn rustred_reduce(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        integral_specs: &ReplacementRules,
    ) -> Result<Atom, VakintError> {
        self.rustred_reduce_observed(settings, numerator, integral_specs, None)
    }

    pub(crate) fn rustred_reduce_observed(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        integral_specs: &ReplacementRules,
        mut observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    ) -> Result<Atom, VakintError> {
        if settings.rustred_application_options.reduction_strategy
            == crate::RustRedReductionStrategy::Ordered
        {
            return self.rustred_reduce_ordered(
                settings,
                vec![MatchedContribution {
                    numerator: numerator.to_owned(),
                    matched: integral_specs.clone(),
                }],
                observation,
                None,
            );
        }
        // Raw contributions must combine before the public notation is restored.
        let mut scalar_settings = settings.clone();
        scalar_settings.use_dot_product_notation = true;
        let settings = &scalar_settings;
        let options = &RustRedEvaluationOptions {
            master_mode: RustRedMasterMode::Retain,
        };
        crate::citations::USED_CITATIONS.record(crate::citations::CitationSource::RustRed);
        if integral_specs.canonical_topology.get_integral().n_loops == 5 {
            MonitorContext::native_phase(
                observation.as_deref_mut(),
                "saved_application_preparation",
                None,
            );
            return self.saved_application(settings)?.evaluate(
                settings,
                numerator,
                integral_specs,
                observation.as_deref_mut(),
            );
        }
        MonitorContext::native_phase(observation, "lower_loop_reduction_unobserved", None);
        if integral_specs.canonical_topology.get_integral().n_loops == 4 {
            return four_loop::evaluate(self, settings, numerator, integral_specs, options);
        }
        MatchedScalarFamily::try_new(integral_specs)?
            .evaluate(self, settings, numerator, options)
            .map_err(Into::into)
    }
}

#[cfg(test)]
mod bundled_cache_tests;
