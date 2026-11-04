//! FORM-independent scalar reduction through shipped RustRed programs.
//!
//! Every request reduces with RustRed's ordered scheduler: all matched
//! contributions of a loop count share one exact frontier per compatible
//! coefficient channel (see `ordered.rs`).
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
pub(crate) mod cache;
pub use cache::{
    DEFAULT_ORDERED_CACHE_MAX_RESIDENT_BYTES, OrderedCacheAutosave, OrderedCacheMode,
    OrderedCacheProgramSaveReport, OrderedCacheProgramStatistics, OrderedCacheSaveReport,
    OrderedCacheSettings, OrderedCacheStatistics, OrderedCacheWarmReport,
    RustRedOrderedStepCacheSaveReport, RustRedOrderedStepCacheShardStatistics,
    RustRedOrderedStepCacheStatistics,
};
pub(crate) mod shipped_cache;
pub use shipped_cache::{
    METADATA_SCHEMA as SHIPPED_ORDERED_CACHE_SCHEMA, OrderedCacheExportReport,
    ShippedOrderedCacheCheck, ShippedOrderedCacheMetadata, ShippedOrderedCacheProgram,
    ShippedOrderedCacheStatus, ShippedOrderedCacheUse, ordered_cache_producer,
};
mod factor_projection;
mod factor_provider;
pub(crate) mod five_loop;
pub(crate) use five_loop::rustred_evaluates_loop;
mod four_loop;
mod mass_channels;
mod matching;
mod materialize;
pub(crate) mod monitor;
mod numerator;
mod ordered;
mod parallel;
mod request_dump;
#[cfg(test)]
mod terminal;
#[cfg(test)]
mod test_data;

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
    #[error("RustRed application lock is poisoned")]
    CacheLockPoisoned,
    #[error("RustRed cache operation failed: {source}")]
    CacheOperation { source: rustred_app::AppError },
    #[error("RustRed Ordered step cache: {detail}")]
    OrderedStepCache { detail: String },
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

/// Worker threads for the independent caller-side work around a native
/// reduction: guard and retained-coefficient materialization and the retained
/// catalogue combination. Each item is computed exactly as on one thread and
/// the results are consumed in input order, so the output never depends on
/// this count; `1` runs everything inline on the calling thread. Follows
/// [`VakintSettings::n_cores`].
pub(crate) fn postprocessing_threads(settings: &VakintSettings) -> usize {
    settings.n_cores.max(1)
}

/// The Ordered scheduler of four- and five-loop requests: worker threads
/// follow [`VakintSettings::n_cores`] (`1` is RustRed's sequential scheduler
/// on the calling thread, with no pool), with RustRed's default window. The
/// result is identical for every setting. Closed one- to three-loop programs
/// always use the sequential scheduler.
pub(crate) fn ordered_scheduler(
    settings: &VakintSettings,
) -> rustred::reduction::OrderedSchedulerOptions {
    rustred::reduction::OrderedSchedulerOptions::with_threads(settings.n_cores.max(1))
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

/// Native application options of one request. Vakint always reduces with
/// RustRed's ordered scheduler and preserves exact uncovered leaves; the
/// `reduction_strategy` of [`VakintSettings::rustred_application_options`] is
/// ignored. The Ordered weight arithmetic comes from
/// [`VakintSettings::ordered_weight_arithmetic`]; the scheduler threads and
/// the five-loop route-verification workers follow [`VakintSettings::n_cores`].
pub(crate) fn application_options(
    settings: &VakintSettings,
) -> rustred_app::SavedApplicationOptions {
    rustred_app::SavedApplicationOptions {
        reduction_strategy: rustred::reduction::ReductionStrategy::Ordered,
        unresolved_policy: rustred::solver::UnresolvedPolicy::Preserve,
        ordered_weight_arithmetic: five_loop_ordered_weight_arithmetic(settings),
        ordered_scheduler: ordered_scheduler(settings),
        preparation_workers: settings.n_cores,
        ..settings.rustred_application_options
    }
}

/// Ordered weight arithmetic of five-loop saved-application requests.
///
/// RustRed composes the raw terminals of a modular request with the package's
/// published terminal rows modulo every prime, before the rational
/// reconstruction (compose before lift), so only the coefficients on the
/// masters are reconstructed. Their degree in `d` is far below the degree of
/// the raw terminal coefficients: `modular_sampled` replays the recorded walk
/// at as many points as that degree needs, while `modular_polynomial` carries
/// weights of the raw-terminal size through every replay. With the default
/// selection ([`crate::OrderedWeightArithmetic::is_default_selection`]),
/// five-loop requests therefore use `modular_sampled` on the same workers; an
/// explicit selection, including an explicit `modular_polynomial`, is kept.
/// Every arithmetic returns the same exact result.
pub(crate) fn five_loop_ordered_weight_arithmetic(
    settings: &VakintSettings,
) -> crate::RustRedOrderedWeightArithmetic {
    let selection = settings.ordered_weight_arithmetic;
    if !selection.is_default_selection() {
        return settings.rustred_ordered_weight_arithmetic();
    }
    let sampled = crate::OrderedWeightArithmetic::modular_sampled();
    match selection.threads() {
        Some(threads) => sampled.with_threads(threads),
        None => sampled,
    }
    .resolve(settings.n_cores)
}

/// Whether the data needed to reduce a supported topology is loaded. One to
/// three loops use the vendored artifacts; four and five loops need external
/// data or, at five loops, a configured application package.
pub(crate) fn has_data(vakint: &Vakint, settings: &VakintSettings, topology: &Topology) -> bool {
    match topology.get_integral().n_loops {
        4 => vakint
            .rustred_data
            .as_ref()
            .is_some_and(|data| data.four_loop.is_some()),
        5 => settings.rustred_artifact_path.is_some() || vakint.has_rustred_loop(5),
        _ => true,
    }
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
        // Ordered reduction combines all normalized contributions first, so
        // simplify real logarithms and zero powers once here, before master
        // substitution, rather than per contribution.
        let raw = crate::utils::simplify_real(raw.as_view());
        let data = self.rustred_data()?;
        let announced = self.announced_bare_masters()?;
        let result = if mode == RustRedMasterMode::Retain {
            data.masters.retain(
                raw,
                &announced,
                observation,
                postprocessing_threads(settings),
            )
        } else {
            if !announced.is_empty() {
                data.masters
                    .require_coefficients(&combine_exact(raw.clone()), &announced)?;
            }
            data.masters.finalize(raw, settings, loop_count, mode)?
        };
        // The catalogue checks its own substitution; every RustRed entry
        // point (single integral, ordered batch, term) passes here as well.
        crate::truncation_guard::ensure_known_orders(
            &result,
            crate::utils::vakint_macros::vk_symbol!(&settings.epsilon_symbol),
            settings
                .number_of_terms_in_epsilon_expansion
                .checked_sub(i64::try_from(loop_count).map_err(data::failure)?)
                .and_then(|value| value.checked_sub(1)),
            "RustRed",
            crate::truncation_guard::plain_name,
        )?;
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
        observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    ) -> Result<Atom, VakintError> {
        self.rustred_reduce_ordered(
            settings,
            vec![MatchedContribution {
                numerator: numerator.to_owned(),
                matched: integral_specs.clone(),
            }],
            observation,
            None,
        )
    }
}

#[cfg(test)]
mod selection_tests {
    use crate::{
        EvaluationMethod, EvaluationOrder, FMFTOptions, MATADOptions, RustRedEvaluationOptions,
        Vakint, VakintSettings, topologies::Topology, utils::vakint_macros::vk_parse,
    };

    const FOUR_LOOP_H: &str = "topo(prop(1,edge(5,1),k(1),muvsq,1)*prop(2,edge(2,6),k(2),muvsq,1)\
        *prop(3,edge(6,5),k(3),muvsq,1)*prop(4,edge(3,4),k(4),muvsq,1)\
        *prop(5,edge(4,5),k(1)-k(3),muvsq,1)*prop(6,edge(6,3),k(2)-k(3),muvsq,1)\
        *prop(7,edge(4,1),k(3)-k(1)+k(4),muvsq,1)*prop(8,edge(2,3),k(3)-k(2)+k(4),muvsq,1)\
        *prop(9,edge(1,2),k(3)+k(4),muvsq,1))";

    fn topology(engine: &Vakint, integral: &str) -> Topology {
        let mut matched = engine
            .topologies
            .match_topologies_to_user_input(vk_parse!(integral).unwrap().as_view(), false)
            .unwrap()
            .unwrap();
        matched.apply_replacement_rules().unwrap();
        matched.canonical_topology
    }

    fn order(methods: Vec<EvaluationMethod>) -> VakintSettings {
        VakintSettings {
            evaluation_order: EvaluationOrder(methods),
            ..Default::default()
        }
    }

    #[test]
    fn rustred_yields_to_a_later_backend_only_without_its_data() {
        let rustred = EvaluationMethod::RustRed(RustRedEvaluationOptions::default());
        let with_fmft = order(vec![
            rustred.clone(),
            EvaluationMethod::FMFT(FMFTOptions::default()),
        ]);
        let rustred_only = order(vec![rustred.clone()]);
        let bare = Vakint::new().unwrap();
        let four_loop = topology(&bare, FOUR_LOOP_H);
        assert!(matches!(
            bare.evaluation_method(&with_fmft, &four_loop),
            Some(EvaluationMethod::FMFT(_))
        ));
        // Without a supporting alternative RustRed reports its missing data.
        assert!(matches!(
            bare.evaluation_method(&rustred_only, &four_loop),
            Some(EvaluationMethod::RustRed(_))
        ));
        let loaded = super::test_data::external_engine();
        assert!(matches!(
            loaded.evaluation_method(&with_fmft, &four_loop),
            Some(EvaluationMethod::RustRed(_))
        ));
        // One to three loops use the vendored programs.
        let with_matad = order(vec![
            rustred,
            EvaluationMethod::MATAD(MATADOptions::default()),
        ]);
        let three_loop = topology(&bare, "topo(I3L(muvsq,1,1,1,1,1,1))");
        assert!(matches!(
            bare.evaluation_method(&with_matad, &three_loop),
            Some(EvaluationMethod::RustRed(_))
        ));
    }
}

#[cfg(test)]
mod five_loop_arithmetic_tests {
    use super::five_loop_ordered_weight_arithmetic;
    use crate::{OrderedWeightArithmetic, VakintSettings};

    /// Only the default selection switches five-loop requests to
    /// `modular_sampled` (on `n_cores` workers, or on the workers it fixes);
    /// an explicit selection is kept, also an explicit `modular_polynomial`.
    #[test]
    fn five_loop_requests_use_modular_sampled_only_for_the_default_selection() {
        let base = VakintSettings {
            n_cores: 3,
            ..Default::default()
        };
        let resolved = five_loop_ordered_weight_arithmetic(&base);
        assert_eq!(
            (resolved.name(), resolved.threads()),
            ("modular_sampled", 3)
        );
        let fixed = VakintSettings {
            ordered_weight_arithmetic: OrderedWeightArithmetic::default().with_threads(2),
            ..base.clone()
        };
        let resolved = five_loop_ordered_weight_arithmetic(&fixed);
        assert_eq!(
            (resolved.name(), resolved.threads()),
            ("modular_sampled", 2)
        );
        for explicit in [
            OrderedWeightArithmetic::modular_polynomial(),
            "modular_polynomial".parse().unwrap(),
            OrderedWeightArithmetic::exact(),
        ] {
            let settings = VakintSettings {
                ordered_weight_arithmetic: explicit,
                ..base.clone()
            };
            assert_eq!(
                five_loop_ordered_weight_arithmetic(&settings),
                settings.rustred_ordered_weight_arithmetic()
            );
            assert_ne!(
                five_loop_ordered_weight_arithmetic(&settings).name(),
                "modular_sampled"
            );
        }
    }
}
