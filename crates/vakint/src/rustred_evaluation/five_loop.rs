//! Exact symbolic application of a separately loaded immutable publication.
//!
//! Masters denote unit-mass integrals in the FMFT/MATAD convention, qualified by
//! the authenticated output-family fingerprint. The requested Minkowski measure
//! and mass scaling are restored once. Coefficients are exact in d=4-2*epsilon.
//! Saved-program guards depend only on d; caller numerator and common-mass pole
//! guards retain their scalar parameters. RustRed owns scalar lowering and rule
//! application; this adapter performs no IBP search or master substitution.

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
    path::{Path, PathBuf},
    sync::{Arc, LazyLock, Mutex},
};

use super::experimental::timing::MathSpan;
use super::monitor::{MonitorContext, RustRedMonitorOptions, RustRedMonitorScope};
use rustred::algebra::Coefficient;
use rustred::family::IntegralKey;
use rustred::input::{Compiler, Limits, LoweringLimits, TextProject, TextPropagator};
use rustred::reduction::terminal_normalization::{
    VacuumFactor, VacuumProductCacheStatistics, VacuumProductChartLimits, VacuumProductProjection,
    VacuumProductReducer, VacuumProductReductionLimits, VacuumTensorProjection,
    VacuumTensorProjectionInput, VerifiedVacuumProductChart,
};
use rustred::scalar_numerator::{FamilyScalarNumeratorService, ScalarNumeratorLimits};
use rustred::solver::CandidateRoutedApplicationObserver;
use rustred_app::{
    PreparedSavedApplication, SavedApplicationOptions, SavedApplicationProductChartStatistics,
    SavedApplicationProductStatistics, load_saved_application,
};
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};

use super::{
    MatchedContribution, RustRedEvaluationError, factor_projection,
    factor_provider::ShippedFactorProvider,
    numerator::{MatchedScalarInput, ScalarNumeratorForm},
};
use crate::topologies::five_loop::{ParentDescriptor, common_momentum_rows, descriptors};
use crate::utils::vakint_macros::vk_symbol;
use crate::{
    ReplacementRules, Topology, Vakint, VakintError, VakintExpression, VakintSettings, symbols::S,
};

static PARENTS: LazyLock<Vec<ParentDescriptor>> = LazyLock::new(descriptors);

#[path = "five_loop/witness_cache.rs"]
mod witness_cache;
pub use witness_cache::RustRedWitnessCacheStatistics;
use witness_cache::{Policy as WitnessPolicy, WitnessCache};

#[cfg(all(test, feature = "feynkit-ingress"))]
#[path = "five_loop/r2_acceptance.rs"]
mod r2_acceptance;

#[cfg(test)]
#[path = "five_loop/admission_census_tests.rs"]
mod admission_census_tests;

#[cfg(test)]
#[path = "five_loop/entry_probe_diagnostics.rs"]
mod entry_probe_diagnostics;

#[cfg(test)]
#[path = "five_loop/test_profiles.rs"]
mod test_profiles;

#[cfg(test)]
#[path = "five_loop/partial_acceptance.rs"]
mod partial_acceptance;

#[cfg(all(
    test,
    feature = "weighted-lowering-diagnostics",
    feature = "feynkit-ingress"
))]
#[path = "five_loop/weighted_diagnostic.rs"]
mod weighted_diagnostic;

#[cfg(all(test, feature = "feynkit-ingress"))]
#[path = "../../tests/fixtures/qcd_ladder/slow.rs"]
mod slow;

#[cfg(all(test, feature = "feynkit-ingress"))]
#[path = "../../tests/fixtures/qcd_ladder/diagnostics.rs"]
mod qcd_diagnostics;

/// Exact symbolic value with the saved publication's generic-d assumptions.
#[derive(Debug, Clone)]
pub struct RustRedSymbolicEvaluation {
    pub value: Atom,
    pub nonzero_conditions: Vec<Atom>,
    /// Master family labels map back to their full authenticated semantic identities.
    pub families: BTreeMap<Atom, String>,
}

/// Declared before batch locals, so its diagnostic runs after their natural drop.
/// It retains no algebraic value and is armed only on completed batch paths.
#[derive(Default)]
struct RetainedBatchCleanup(Option<std::time::Instant>);

impl RetainedBatchCleanup {
    fn arm(&mut self) {
        self.0 = Some(std::time::Instant::now());
    }
}

impl Drop for RetainedBatchCleanup {
    fn drop(&mut self) {
        if let Some(started) = self.0 {
            log::debug!(target: "vakint::saved_application", "phase=retained_batch_cleanup elapsed_seconds={:.6}", started.elapsed().as_secs_f64());
        }
    }
}

pub(crate) struct SavedApplication {
    pub(super) prepared: Mutex<PreparedSavedApplication<15>>,
    pub(super) cache_digest: [u8; 32],
    pub(super) cache: Mutex<Option<super::cache::SavedFiveLoopCache>>,
    fingerprint: String,
    data_version: String,
    products: SavedApplicationProductStatistics,
    product_charts: SavedApplicationProductChartStatistics,
    product_reducer: Mutex<VacuumProductReducer<ShippedFactorProvider>>,
    data: Arc<super::data::RustRedData>,
    witnesses: Mutex<WitnessCache>,
}

#[derive(Debug, Clone)]
pub(crate) struct LoadedApplication {
    pub(super) path: PathBuf,
    pub(super) options: SavedApplicationOptions,
    pub(super) product_limits: VacuumProductReductionLimits,
    pub(super) application: Arc<SavedApplication>,
    bare: bool,
    catalogue: Arc<crate::topologies::FiveLoopCatalogue>,
}

impl fmt::Debug for SavedApplication {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SavedApplication")
            .field("family", &self.fingerprint)
            .finish_non_exhaustive()
    }
}

fn failure(error: impl fmt::Display) -> RustRedEvaluationError {
    RustRedEvaluationError::Reduction {
        detail: error.to_string(),
    }
}

fn native_failure(error: impl std::error::Error + 'static) -> RustRedEvaluationError {
    if let Some(resource) = RustRedEvaluationError::resource_cause(&error) {
        RustRedEvaluationError::ResourceLimit(resource.clone())
    } else {
        failure(error)
    }
}

fn project_chart(
    chart: &VerifiedVacuumProductChart,
    dot_head: Symbol,
    loop_momenta: Vec<Atom>,
    limits: ScalarNumeratorLimits,
    mut projector: impl FnMut(
        &VacuumTensorProjectionInput,
    ) -> Result<VacuumTensorProjection, RustRedEvaluationError>,
) -> Result<VacuumProductProjection, RustRedEvaluationError> {
    let mut callback_resource = None;
    let projected = chart.project(dot_head, loop_momenta, limits, |input| {
        let result = projector(input);
        if let Err(error) = &result {
            callback_resource = RustRedEvaluationError::resource_cause(error).cloned();
        }
        result
    });
    if let Some(resource) = callback_resource {
        return Err(RustRedEvaluationError::ResourceLimit(resource));
    }
    projected.map_err(native_failure)
}

fn select(topology: &Topology) -> Result<&'static ParentDescriptor, RustRedEvaluationError> {
    if !matches!(topology, Topology::FiveLoop(_)) {
        return Err(failure(
            "saved five-loop application requires a registered five-loop topology",
        ));
    }
    PARENTS
        .iter()
        .find(|parent| {
            topology
                .get_integral()
                .validate_parent_routing(&parent.momenta)
                .is_ok()
        })
        .ok_or_else(|| {
            failure("no saved application descriptor matches the retained parent routing")
        })
}

pub(super) fn supports(topology: &Topology) -> bool {
    select(topology).is_ok()
}

impl Vakint {
    pub(crate) fn match_rustred_input(
        &self,
        settings: &VakintSettings,
        integral: AtomView,
        numerator: AtomView,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        // Select and validate the active immutable package before any memo hit.
        // The ordinary no-package matching path retains its previous behavior.
        if settings.rustred_artifact_path.is_none() && !self.has_rustred_loop(5) {
            return self.match_rustred_input_uncached(settings, integral, numerator);
        }
        let application = self.saved_application(settings)?;
        // Loading also installs the catalogue in a clone created before load.
        let registry = self.topologies.registry_identity();
        let policy = WitnessPolicy::from(settings);
        let cached = application
            .witnesses
            .lock()
            .ok()
            .and_then(|mut cache| cache.get(&registry, &policy, integral, numerator));
        if let Some(cached) = cached {
            // The memo lock is released before native lowering/package locking.
            // Typed policy errors propagate; a scope miss invalidates the entry.
            match self.rustred_matched_input_key(settings, numerator, &cached.witness)? {
                Some(key) if key == cached.accepted => {
                    log::debug!(target: "vakint::saved_application", "phase=frontend_witness_memo_hit");
                    return Ok(Some(cached.witness.clone()));
                }
                _ => {
                    if let Ok(mut cache) = application.witnesses.lock() {
                        cache.invalidate(&cached);
                    }
                }
            }
        }
        let matched = self.match_rustred_input_uncached(settings, integral, numerator)?;
        if let Some(witness) = &matched {
            // A verified fallback outside recorded coverage remains usable by
            // the partial reducer. Only covered witnesses enter this memo.
            if matches!(witness.canonical_topology, Topology::FiveLoop(_))
                && let Some(accepted) =
                    self.rustred_matched_input_key(settings, numerator, witness)?
                && let Ok(mut cache) = application.witnesses.lock()
            {
                cache.insert(registry, policy, integral, numerator, witness, accepted);
            }
        }
        Ok(matched)
    }

    fn match_rustred_input_uncached(
        &self,
        settings: &VakintSettings,
        integral: AtomView,
        numerator: AtomView,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        #[cfg(all(
            test,
            feature = "weighted-lowering-diagnostics",
            feature = "feynkit-ingress"
        ))]
        let diagnostic_thread = weighted_diagnostic::begin_search(integral, numerator)?;
        let result = self
            .topologies
            .match_topologies_to_user_input_with_acceptance(
                integral,
                settings.allow_unknown_integrals,
                settings.rustred_momentum_match_limits,
                |matched| {
                    #[cfg(all(
                        test,
                        feature = "weighted-lowering-diagnostics",
                        feature = "feynkit-ingress"
                    ))]
                    weighted_diagnostic::assert_thread(diagnostic_thread)?;
                    self.rustred_matched_input_key(settings, numerator, matched)
                },
            );
        #[cfg(all(
            test,
            feature = "weighted-lowering-diagnostics",
            feature = "feynkit-ingress"
        ))]
        if let Ok(selected) = &result {
            weighted_diagnostic::end_search(selected.as_ref())?;
        }
        result
    }

    /// Order already verified matcher witnesses by their exact surviving
    /// common-family integrand. This checks the saved entry union without
    /// applying any rule or preparing an owner.
    pub(crate) fn rustred_matched_input_key(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        matched: &ReplacementRules,
    ) -> Result<Option<Vec<(IntegralKey, Atom)>>, VakintError> {
        #[cfg(all(
            test,
            feature = "weighted-lowering-diagnostics",
            feature = "feynkit-ingress"
        ))]
        let diagnostic_original = matched;
        let mut matched = matched.clone();
        matched.apply_replacement_rules()?;
        let mut term = crate::VakintTerm {
            integral: Atom::one(),
            numerator: numerator.to_owned(),
            vectors: Vec::new(),
        };
        term.apply_scalar_numerator_replacement_rules(&matched, settings)?;
        let descriptor = select(&matched.canonical_topology)?;
        let input = MatchedScalarInput::prepare(
            term.numerator.as_view(),
            &matched,
            &descriptor.momenta,
            descriptor.physical_slots,
            15,
            settings,
            ScalarNumeratorForm::PreserveProducts,
        )?;
        let application = self.saved_application(settings)?;
        let prepared = application.prepared.lock().map_err(failure)?;
        let scalar = FamilyScalarNumeratorService::try_new(
            prepared.family(),
            S.dot,
            (1..=5).map(|axis| S.k.call(axis)).collect(),
            settings.rustred_scalar_numerator_limits,
        )
        .map_err(native_failure)?;
        #[cfg(all(
            test,
            feature = "weighted-lowering-diagnostics",
            feature = "feynkit-ingress"
        ))]
        let diagnostic_before = weighted_diagnostic::before_lowering()?;
        let lowered = scalar
            .lower_weighted(&input.numerator, &input.target, &input.mass_squared)
            .map_err(native_failure)?;
        let mut key = Vec::with_capacity(lowered.terms().len());
        for term in lowered.terms() {
            if !prepared
                .routed_entry_is_admitted(term.integral())
                .map_err(native_failure)?
            {
                log::debug!(target: "vakint::five_loop_match", "phase=candidate_rejected lowered_count={} rejected_key={:?} lowered_keys_prefix={:?}", lowered.terms().len(), term.integral().powers(), lowered.terms().iter().take(4).map(|term| term.integral().powers()).collect::<Vec<_>>());
                #[cfg(all(
                    test,
                    feature = "weighted-lowering-diagnostics",
                    feature = "feynkit-ingress"
                ))]
                weighted_diagnostic::candidate(
                    diagnostic_original,
                    &input,
                    &lowered,
                    Some(key.len()),
                    None,
                    diagnostic_before,
                )?;
                return Ok(None);
            }
            key.push((
                term.integral().clone(),
                super::combine_exact(term.coefficient().to_expression() * term.scalar_spectator()),
            ));
        }
        key.sort();
        #[cfg(all(
            test,
            feature = "weighted-lowering-diagnostics",
            feature = "feynkit-ingress"
        ))]
        weighted_diagnostic::candidate(
            diagnostic_original,
            &input,
            &lowered,
            None,
            Some(&key),
            diagnostic_before,
        )?;
        Ok(Some(key))
    }

    /// Load one authenticated application package and reuse its owned memoization.
    ///
    /// Master substitution uses the active external-data catalog. Uncomputed
    /// terminals remain symbolic or produce an explicit numerical-data error.
    /// Clones share this immutable package and synchronized application cache.
    /// Only successful loading prepares and publishes the five-loop catalogue;
    /// `Vakint::new` prepares just the ordinary lower-loop registry.
    pub fn with_rustred_application(
        mut self,
        settings: &mut VakintSettings,
        path: impl AsRef<Path>,
    ) -> Result<Self, VakintError> {
        let mut configured = settings.clone();
        let path = path.as_ref();
        configured.rustred_artifact_path = Some(path.canonicalize().map_err(failure)?);
        let data = self.rustred_data()?;
        let root = super::data::RustRedData::release_root(path)?
            .or_else(|| data.directory_root_for_ibps(5, path));
        if let Some(root) = root
            && !data.matches_directory(&root, 5)?
        {
            self.load_external_data_with_prepared_caches(
                root,
                5,
                configured.rustred_application_options.reduction_strategy
                    == crate::RustRedReductionStrategy::DepthFirst,
            )?;
        }
        self.saved_application(&configured)?;
        *settings = configured;
        Ok(self)
    }

    pub(crate) fn saved_application(
        &self,
        settings: &VakintSettings,
    ) -> Result<Arc<SavedApplication>, VakintError> {
        let data = self.rustred_data()?;
        let official_path = data
            .loops
            .get(&5)
            .map(|pack| pack.path("ibps"))
            .transpose()?
            .flatten();
        let path = settings
            .rustred_artifact_path
            .as_ref()
            .or(official_path.as_ref())
            .ok_or(RustRedEvaluationError::SavedApplicationUnavailable)?;
        let path = path.canonicalize().map_err(failure)?;
        let mut options = settings.rustred_application_options;
        // Strategy changes do not change immutable tables or the DFS cache
        // namespace. Select the request engine under the prepared-owner lock.
        options.reduction_strategy = crate::RustRedReductionStrategy::DepthFirst;
        options.unresolved_policy = rustred::solver::UnresolvedPolicy::Preserve;
        let mut loaded = self.rustred_application.lock().map_err(failure)?;
        if loaded.as_ref().is_none_or(|current| {
            current.path != path
                || current.options != options
                || current.product_limits != settings.rustred_product_limits
        }) {
            let bare = super::data::RustRedData::release_root(&path)?.is_none();
            let started = std::time::Instant::now();
            let cache_digest = super::cache::saved_application_digest(&data, &path)?;
            log::debug!(target: "vakint::saved_application", "phase=native_load_begin");
            let prepared = load_saved_application::<15>(&path, options).map_err(native_failure)?;
            if super::cache::saved_application_digest(&data, &path)? != cache_digest {
                return Err(failure(
                    "saved application changed while its native payload was loading",
                )
                .into());
            }
            log::debug!(target: "vakint::saved_application", "phase=native_load_end elapsed_seconds={:.6}", started.elapsed().as_secs_f64());
            let application = Arc::new(SavedApplication::new(
                prepared,
                settings.rustred_product_limits,
                Arc::clone(&data),
                cache_digest,
            )?);
            if let Some(previous) = loaded.as_ref()
                && previous.path == path
                && previous.application.cache_digest == cache_digest
            {
                application.inherit_compatible_cache(&previous.application)?;
            }
            // Only a fully authenticated, prepared publication can activate
            // the optional catalogue. Cache/options reloads reuse its Arc.
            let catalogue = match self
                .topologies
                .five_loop_catalogue()
                .or_else(|| loaded.as_ref().map(|current| &current.catalogue))
            {
                Some(catalogue) => Arc::clone(catalogue),
                None => crate::topologies::FiveLoopCatalogue::prepare()?,
            };
            *loaded = Some(LoadedApplication {
                path,
                options,
                product_limits: settings.rustred_product_limits,
                application,
                bare,
                catalogue,
            });
        }
        let loaded = loaded
            .as_ref()
            .expect("successful load installed a package");
        // A Vakint cloned before the shared package was loaded has its own
        // empty registry view; cache hits must publish to that view as well.
        self.topologies
            .install_five_loop_catalogue(Arc::clone(&loaded.catalogue));
        Ok(Arc::clone(&loaded.application))
    }

    pub(super) fn loaded_application(&self) -> Result<Arc<SavedApplication>, VakintError> {
        let loaded = self.rustred_application.lock().map_err(failure)?;
        loaded
            .as_ref()
            .map(|current| Arc::clone(&current.application))
            .ok_or_else(|| RustRedEvaluationError::SavedApplicationUnavailable.into())
    }

    /// Inspect the already authenticated native inventory without loading data
    /// or manufacturing coefficient/graph metadata for a bare application.
    pub(crate) fn announced_bare_masters(
        &self,
    ) -> Result<super::masters::BareMasterDeclarations, VakintError> {
        let (path, application) = {
            let loaded = self.rustred_application.lock().map_err(failure)?;
            let Some(source) = loaded.as_ref().filter(|source| source.bare) else {
                return Ok(BTreeMap::new());
            };
            (source.path.clone(), Arc::clone(&source.application))
        };
        // Release the session lock before acquiring the native application.
        // All declarations belong to this immutable prepared snapshot.
        let data = self.rustred_data()?;
        let prepared = application.prepared.lock().map_err(failure)?;
        Ok(prepared
            .remaining_terminals()
            .iter()
            .map(|terminal| {
                (
                    data.masters
                        .symbol(terminal.family_fingerprint(), terminal.integral().powers()),
                    super::masters::BareMasterDeclaration {
                        application_path: path.clone(),
                    },
                )
            })
            .collect())
    }

    pub(crate) fn require_bare_master_coefficients(
        &self,
        expression: &Atom,
    ) -> Result<(), VakintError> {
        let announced = self.announced_bare_masters()?;
        if !announced.is_empty() {
            self.rustred_data()?
                .masters
                .require_coefficients(&super::combine_exact(expression.clone()), &announced)?;
        }
        Ok(())
    }

    /// Clear the explicitly loaded package's cache without changing its entry scope.
    pub fn clear_rustred_application_cache(&self) -> Result<(), VakintError> {
        self.loaded_application()?.clear_cache_state(false)
    }

    /// Observe the optional frontend witness memo, not the IBP coefficient cache.
    /// A poisoned memo is bypassed during matching and reported as unavailable.
    pub fn rustred_application_witness_cache_statistics(
        &self,
    ) -> Result<RustRedWitnessCacheStatistics, VakintError> {
        Ok(self
            .loaded_application()?
            .witnesses
            .lock()
            .map_err(failure)?
            .statistics())
    }

    pub fn rustred_application_cache_statistics(
        &self,
    ) -> Result<rustred::solver::CandidateStatistics, VakintError> {
        let application = self.loaded_application()?;
        Ok(application
            .prepared
            .lock()
            .map_err(failure)?
            .cache_statistics())
    }

    /// Cumulative eviction/recomputation work, including failed requests.
    /// Clearing retained entries does not reset these counters. Recomputed
    /// nodes are repeated misses within one request, not cross-request history.
    pub fn rustred_application_cache_work_statistics(
        &self,
    ) -> Result<crate::RustRedCacheWorkStatistics, VakintError> {
        let application = self.loaded_application()?;
        Ok(application
            .prepared
            .lock()
            .map_err(failure)?
            .cache_work_statistics())
    }

    /// One-time exact factor-product preparation, excluded from warm application.
    pub fn rustred_application_product_statistics(
        &self,
    ) -> Result<SavedApplicationProductStatistics, VakintError> {
        Ok(self.loaded_application()?.products.clone())
    }

    /// Numerator-aware geometry only; projections run when their terminals occur.
    pub fn rustred_application_product_chart_statistics(
        &self,
    ) -> Result<SavedApplicationProductChartStatistics, VakintError> {
        Ok(self.loaded_application()?.product_charts.clone())
    }

    pub fn rustred_application_owner_preparation_statistics(
        &self,
    ) -> Result<rustred::solver::DeferredCandidateOwnerStatistics, VakintError> {
        Ok(self
            .loaded_application()?
            .prepared
            .lock()
            .map_err(failure)?
            .owner_preparation_statistics())
    }

    /// Qualified lower-factor memoization, separate from saved-entry memoization.
    pub fn rustred_application_factor_cache_statistics(
        &self,
    ) -> Result<VacuumProductCacheStatistics, VakintError> {
        Ok(self
            .loaded_application()?
            .product_reducer
            .lock()
            .map_err(failure)?
            .cache_statistics())
    }

    /// Evaluate five-loop input while retaining each exact nonzero condition.
    /// This uses the same matcher, tensor prepass and scalar tail as `evaluate`.
    pub fn evaluate_rustred_symbolic(
        &self,
        settings: &VakintSettings,
        input: AtomView,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        self.evaluate_rustred_symbolic_observed(settings, input, None, None)
    }

    /// Experimental five-loop diagnostics on the ordinary public evaluation.
    /// The borrowed callback observes Native boundaries without presentation;
    /// it must not block or reenter the engine. Observation faults preserve the
    /// mathematical result; callers must inspect the observer's final status.
    #[cfg(feature = "experimental-rustred")]
    pub fn evaluate_rustred_symbolic_with_observer(
        &self,
        settings: &VakintSettings,
        input: AtomView,
        observer: &mut CandidateRoutedApplicationObserver<'_>,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        self.evaluate_rustred_symbolic_observed(settings, input, Some(observer), None)
    }

    /// Retain exact five-loop masters and guards with optional per-call diagnostics.
    pub fn evaluate_rustred_symbolic_with_monitor(
        &self,
        settings: &VakintSettings,
        input: AtomView,
        monitor: &RustRedMonitorOptions,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        MonitorContext::run(monitor, |observation, progress| {
            self.evaluate_rustred_symbolic_observed(settings, input, observation, progress)
        })
    }

    /// Retain exact masters inside an existing caller-owned job monitor.
    /// Returning completes this API call; the outer job scope remains active.
    pub fn evaluate_rustred_symbolic_in_scope(
        &self,
        settings: &VakintSettings,
        input: AtomView,
        scope: &mut RustRedMonitorScope<'_, '_>,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        scope.evaluate(|observation, progress| {
            self.evaluate_rustred_symbolic_observed(settings, input, observation, progress)
        })
    }

    fn evaluate_rustred_symbolic_observed(
        &self,
        settings: &VakintSettings,
        input: AtomView,
        mut observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
        progress: Option<&super::monitor::CallerProgress>,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        settings.with_resource_budget(|| {
            let ordered = settings.rustred_application_options.reduction_strategy
                == crate::RustRedReductionStrategy::Ordered;
            if ordered && let Some(progress) = progress {
                progress.ordered_batch(0, 1, "5L / Retain");
            }
            let setup_started = std::time::Instant::now();
            MonitorContext::native_phase(observation.as_deref_mut(), "saved_application_preparation", Some(input));
            self.validate_settings(settings)?;
            let application = self.saved_application(settings)?;
            let mut expression = VakintExpression::try_from(input)?;
            log::debug!(target: "vakint::saved_application", "phase=frontend_setup elapsed_seconds={:.6}", setup_started.elapsed().as_secs_f64());
            log::debug!(target: "vakint::five_loop_match", "phase=symbolic_input input={input}");
            let canonical_started = std::time::Instant::now();
            MonitorContext::native_phase(observation.as_deref_mut(), "canonicalization", None);
            expression.canonicalize(settings, &self.topologies, false)?;
            log::debug!(target: "vakint::saved_application", "phase=frontend_canonicalization elapsed_seconds={:.6}", canonical_started.elapsed().as_secs_f64());
            log::debug!(target: "vakint::five_loop_match", "phase=symbolic_canonical input={expression}");
            // This is an internal scalar handoff, not the public tensor-output
            // presentation boundary. Reopening scalar powers here can manufacture
            // a rank22 indexed expression before a rank0 scalar application.
            let mut tensor_settings = settings.clone();
            tensor_settings.use_dot_product_notation = true;
            let tensor_started = std::time::Instant::now();
            MonitorContext::native_phase(observation.as_deref_mut(), "tensor_prepass", None);
            expression.tensor_reduce(self, &tensor_settings)?;
            log::debug!(target: "vakint::saved_application", "phase=frontend_tensor_prepass elapsed_seconds={:.6}", tensor_started.elapsed().as_secs_f64());
            log::debug!(target: "vakint::five_loop_match", "phase=symbolic_tensor_reduced input={expression}");
            let mut value = Atom::Zero;
            let mut conditions = BTreeSet::new();
            let mut families = BTreeMap::new();
            let mut ordered_contributions = Vec::new();
            let split_started = std::time::Instant::now();
            let expression = Atom::from(expression);
            MonitorContext::native_phase(observation.as_deref_mut(), "scalar_input", Some(expression.as_view()));
            let terms = VakintExpression::split_integrals(expression.as_view())?;
            drop(expression);
            log::debug!(target: "vakint::saved_application", "phase=scalar_input_split elapsed_seconds={:.6}", split_started.elapsed().as_secs_f64());
            for mut term in terms {
                symbolica::utils::resources::checkpoint()?;
                let matching_started = std::time::Instant::now();
                MonitorContext::native_phase(observation.as_deref_mut(), "routing_verification", None);
                let mut matched = self
                    .match_rustred_input(settings, term.integral.as_view(), term.numerator.as_view())?
                    .ok_or_else(|| failure("input has no registered topology match"))?;
                log::debug!(target: "vakint::saved_application", "phase=frontend_admitted_witness_search elapsed_seconds={:.6}", matching_started.elapsed().as_secs_f64());
                let replacement_started = std::time::Instant::now();
                MonitorContext::native_phase(observation.as_deref_mut(), "numerator_replacement", None);
                matched.apply_replacement_rules()?;
                term.apply_scalar_numerator_replacement_rules(&matched, settings)?;
                log::debug!(target: "vakint::saved_application", "phase=numerator_replacement elapsed_seconds={:.6}", replacement_started.elapsed().as_secs_f64());
                log::debug!(target: "vakint::five_loop_match", "phase=symbolic_candidate_numerator numerator={}", term.numerator);
                if settings.rustred_application_options.reduction_strategy
                    == crate::RustRedReductionStrategy::Ordered
                {
                    ordered_contributions.push(MatchedContribution {
                        numerator: term.numerator,
                        matched,
                    });
                    continue;
                }
                // This parent span includes callee cleanup. Its recorded
                // subphases must not be added to it as disjoint durations.
                let detailed_started = std::time::Instant::now();
                let result =
                    application.evaluate_detailed(&tensor_settings, term.numerator.as_view(), &matched, observation.as_deref_mut(), progress);
                log::debug!(target: "vakint::saved_application", "phase=detailed_evaluation_total elapsed_seconds={:.6}", detailed_started.elapsed().as_secs_f64());
                let result = result?;
                let merge_started = std::time::Instant::now();
                MonitorContext::native_phase(observation.as_deref_mut(), "caller_result_merge", Some(result.value.as_view()));
                {
                    let _math = MathSpan::new();
                    value += result.value;
                    conditions.extend(result.nonzero_conditions);
                    families.extend(result.families);
                }
                log::debug!(target: "vakint::saved_application", "phase=caller_result_merge elapsed_seconds={:.6}", merge_started.elapsed().as_secs_f64());
            }
            if !ordered_contributions.is_empty() {
                let result = application.evaluate_ordered_batch(
                    &tensor_settings,
                    &ordered_contributions,
                    observation.as_deref_mut(),
                    progress,
                )?;
                {
                    let _math = MathSpan::new();
                    value += result.value;
                    conditions.extend(result.nonzero_conditions);
                    families.extend(result.families);
                }
            } else if ordered && let Some(progress) = progress {
                progress.cohorts(0, 0);
            }
            if ordered && let Some(progress) = progress {
                progress.ordered_batch(1, 1, "5L / Retain");
            }
            // The catalogue's own exact-combination and classification spans
            // are children of this complete finalization call.
            let finalization_started = std::time::Instant::now();
            let value = self.finalize_rustred(
                value,
                settings,
                5,
                super::RustRedMasterMode::Retain,
                observation.as_deref_mut(),
            );
            log::debug!(target: "vakint::saved_application", "phase=retained_finalization_total elapsed_seconds={:.6}", finalization_started.elapsed().as_secs_f64());
            let value = value?;
            let packing_started = std::time::Instant::now();
            MonitorContext::native_phase(observation.as_deref_mut(), "retained_public_result_packing", Some(value.as_view()));
            let result = RustRedSymbolicEvaluation {
                value,
                nonzero_conditions: conditions.into_iter().collect(),
                families,
            };
            log::debug!(target: "vakint::saved_application", "phase=public_result_packing elapsed_seconds={:.6}", packing_started.elapsed().as_secs_f64());
            MonitorContext::native_phase(observation, "retained_output", Some(result.value.as_view()));
            Ok(result)
        })
    }
}

impl SavedApplication {
    pub(in crate::rustred_evaluation) fn clear_cache_state(
        &self,
        detach_backing: bool,
    ) -> Result<(), VakintError> {
        // Evaluation, attachment and statistics all acquire Prepared before Cache.
        // Wait for every owner before checking the resource budget, so a refusal
        // cannot leave only some resident caches cleared or the backing detached.
        let mut prepared = self.prepared.lock().map_err(failure)?;
        let mut attached = self.cache.lock().map_err(failure)?;
        let mut products = self.product_reducer.lock().map_err(failure)?;
        let mut witnesses = self.witnesses.lock().map_err(failure)?;
        super::cache::resource_checkpoint()?;
        if detach_backing && let Some(cache) = attached.as_ref() {
            cache.flush()?;
        }
        // These native clears release retained entries without allocation or
        // resource polls. The commit therefore cannot introduce a ResourceError.
        prepared.clear_cache().map_err(super::cache::operation)?;
        products.clear_cache().map_err(native_failure)?;
        witnesses.clear();
        if detach_backing {
            *attached = None;
        }
        Ok(())
    }

    fn new(
        mut prepared: PreparedSavedApplication<15>,
        product_limits: VacuumProductReductionLimits,
        data: Arc<super::data::RustRedData>,
        cache_digest: [u8; 32],
    ) -> Result<Self, RustRedEvaluationError> {
        let started = std::time::Instant::now();
        // Compile the same complete common basis through RustRed's existing
        // input boundary. No denominator inversion is implemented in Vakint.
        let expected = Compiler::new(Limits::default())
            .map_err(native_failure)?
            .compile_text(TextProject {
                name: Some(prepared.family().name().to_owned()),
                parameters: None,
                loop_momenta: (1..=5).map(|axis| format!("k{axis}")).collect(),
                external_momenta: Vec::new(),
                dimension: "d".into(),
                propagators: common_momentum_rows()
                    .iter()
                    .enumerate()
                    .map(|(slot, row)| {
                        let momentum = row
                            .iter()
                            .enumerate()
                            .filter(|(_, coefficient)| **coefficient != 0)
                            .map(|(axis, coefficient)| format!("({coefficient})*k{}", axis + 1))
                            .collect::<Vec<_>>()
                            .join("+");
                        TextPropagator {
                            id: format!("D{}", slot + 1),
                            expression: format!("({momentum})^2-1"),
                            target_power: 0,
                            power_shift: None,
                        }
                    })
                    .collect(),
                external_gram: Vec::new(),
                numerator: None,
            })
            .map_err(native_failure)?
            .into_lowered(LoweringLimits::default())
            .map_err(native_failure)?
            .into_family();
        if prepared.family().fingerprint() != expected.fingerprint() {
            return Err(failure(
                "saved package does not use the registered common momentum basis",
            ));
        }
        log::debug!(target: "vakint::saved_application", "phase=family_verified elapsed_seconds={:.6}", started.elapsed().as_secs_f64());
        let products = prepared
            .prepare_terminal_products(product_limits.products)
            .map_err(native_failure)?;
        let product_charts = prepared
            .prepare_terminal_product_charts(VacuumProductChartLimits {
                products: product_limits.products,
                ..Default::default()
            })
            .map_err(native_failure)?;
        log::debug!(target: "vakint::saved_application", "phase=products_prepared elapsed_seconds={:.6}", started.elapsed().as_secs_f64());
        let anchor = prepared
            .output_family_shared(prepared.family().fingerprint())
            .ok_or_else(|| failure("saved application does not expose its authenticated family"))?;
        let data_version = data.masters.version.clone();
        let product_reducer = VacuumProductReducer::try_new(
            ShippedFactorProvider::with_data(Arc::clone(&data)),
            anchor,
            product_limits,
        )
        .map_err(native_failure)?;
        Ok(Self {
            fingerprint: prepared.family().fingerprint().to_owned(),
            data_version,
            prepared: Mutex::new(prepared),
            cache_digest,
            cache: Mutex::new(None),
            products,
            product_charts,
            product_reducer: Mutex::new(product_reducer),
            data,
            witnesses: Mutex::new(WitnessCache::default()),
        })
    }

    fn symbolic_integral(
        family_fingerprint: &str,
        integral: &IntegralKey,
        head: Symbol,
        families: &mut BTreeMap<Atom, String>,
    ) -> Atom {
        let identity = Atom::var(vk_symbol!(format!(
            "rustred_family_{:x}",
            Sha256::digest(family_fingerprint.as_bytes())
        )));
        families.insert(identity.clone(), family_fingerprint.to_owned());
        let mut arguments = vec![identity];
        arguments.extend(integral.powers().iter().map(|power| Atom::num(*power)));
        head.call(&arguments)
    }

    fn normalization(settings: &VakintSettings, loop_count: usize) -> Result<Atom, VakintError> {
        let normalization = settings
            .get_staged_integral_normalization_factor_atom()?
            .replace(S.n_loops.to_pattern())
            .with(Atom::num(loop_count).to_pattern())
            .replace(Atom::var(vk_symbol!("log_mu_sq")).to_pattern())
            .with(
                Symbol::LOG
                    .call(Atom::var(vk_symbol!(settings.mu_r_sq_symbol.as_str())))
                    .to_pattern(),
            );
        // Public masters use the unit-mass FMFT/MATAD convention, including
        // lower-loop product factors. Normalize both sides of the measure
        // conversion together so the two EulerGamma spellings cancel exactly.
        let measure = crate::utils::vakint_macros::vk_parse!(format!(
            "(𝑖*𝜋^(2-({eps}))*exp(-EulerGamma*({eps})))^{loop_count}",
            eps = settings.epsilon_symbol
        ))
        .map_err(failure)?;
        let factor =
            Vakint::normalize_gamma_functions((normalization * measure).as_view()).expand();
        Ok(crate::utils::simplify_real(factor.as_view()))
    }

    pub(super) fn evaluate(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        matched: &ReplacementRules,
        observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    ) -> Result<Atom, VakintError> {
        Ok(self
            .evaluate_detailed(settings, numerator, matched, observation, None)?
            .value)
    }

    fn evaluate_detailed(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        matched: &ReplacementRules,
        observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
        progress: Option<&super::monitor::CallerProgress>,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        self.evaluate_matched_batch(
            settings,
            &[MatchedContribution {
                numerator: numerator.to_owned(),
                matched: matched.clone(),
            }],
            observation,
            progress,
        )
    }

    pub(super) fn evaluate_ordered_batch(
        &self,
        settings: &VakintSettings,
        contributions: &[MatchedContribution],
        observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
        progress: Option<&super::monitor::CallerProgress>,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        self.evaluate_matched_batch(settings, contributions, observation, progress)
    }

    fn evaluate_matched_batch(
        &self,
        settings: &VakintSettings,
        contributions: &[MatchedContribution],
        mut observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
        progress: Option<&super::monitor::CallerProgress>,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        MonitorContext::native_phase(observation.as_deref_mut(), "scalar_lowering", None);
        let mut application = self.prepared.lock().map_err(failure)?;
        application.set_reduction_strategy(settings.rustred_application_options.reduction_strategy);
        let family = application.family();
        let loop_count = family.loop_count();
        let dimension = family.dimension().to_expression();
        let epsilon = Atom::var(vk_symbol!(settings.epsilon_symbol.as_str()));
        let physical_dimension = Atom::num(4) - Atom::num(2) * &epsilon;
        let scalar = FamilyScalarNumeratorService::try_new(
            family,
            S.dot,
            (1..=loop_count).map(|axis| S.k.call(axis)).collect(),
            settings.rustred_scalar_numerator_limits,
        )
        .map_err(native_failure)?;
        let mut batches = BTreeMap::<(Atom, Atom, i64), Vec<(IntegralKey, Coefficient)>>::new();
        let mut uncovered = BTreeSet::new();
        let mut conditions = BTreeSet::new();
        let mut mass_exponent_origins = BTreeMap::<Atom, i64>::new();
        let mut symbolic_masses = BTreeSet::new();
        let mut input_progress = progress.map(|owner| {
            owner.phase(
                "matched_input_batching",
                contributions.len(),
                "matched contributions lowered",
            )
        });
        for (input_index, contribution) in contributions.iter().enumerate() {
            let lowering_started = std::time::Instant::now();
            let descriptor = select(&contribution.matched.canonical_topology)?;
            let mut input = MatchedScalarInput::prepare(
                contribution.numerator.as_view(),
                &contribution.matched,
                &descriptor.momenta,
                descriptor.physical_slots,
                15,
                settings,
                ScalarNumeratorForm::PreserveProducts,
            )?;
            if settings.rustred_application_options.reduction_strategy
                == crate::RustRedReductionStrategy::Ordered
            {
                // Preserve the received rational domain before identifying
                // epsilon with (4-d)/2 can remove a denominator algebraically.
                conditions.extend(input.input_conditions_after_mapping(&scalar, |condition| {
                    condition
                        .replace(dimension.to_pattern())
                        .with(physical_dimension.to_pattern())
                })?);
                input.numerator = MatchedScalarInput::numerator_in_family_dimension(
                    input.numerator.as_view(),
                    settings,
                    &dimension,
                );
            }
            let lowering = scalar
                .lower_weighted(&input.numerator, &input.target, &input.mass_squared)
                .map_err(native_failure)?;
            log::debug!(target: "vakint::saved_application",
                "phase=parent_transport_and_scalar_lowering elapsed_seconds={:.6} lowered_count={}",
                lowering_started.elapsed().as_secs_f64(), lowering.terms().len());
            MonitorContext::native_phase(
                observation.as_deref_mut(),
                "coverage_mass_batching",
                Some(input.numerator.as_view()),
            );
            let batching_started = std::time::Instant::now();
            let mut coverage_progress = progress.map(|owner| {
                owner.phase(
                    "coverage_mass_batching",
                    lowering.terms().len(),
                    "lowered keys checked",
                )
            });
            for (index, term) in lowering.terms().iter().enumerate() {
                if !application
                    .routed_entry_is_admitted(term.integral())
                    .map_err(native_failure)?
                {
                    uncovered.insert(term.integral().clone());
                }
                let exponent = term
                    .integral()
                    .powers()
                    .iter()
                    .try_fold(0i128, |sum, power| sum.checked_sub(i128::from(*power)))
                    .ok_or_else(|| failure("saved application mass exponent overflow"))?;
                let exponent = i64::try_from(exponent)
                    .map_err(|_| RustRedEvaluationError::MassExponentOverflow { exponent })?;
                let anchor = *mass_exponent_origins
                    .entry(input.mass_squared.clone())
                    .or_insert(exponent);
                let (spectator, exponent) = MatchedScalarInput::mass_channel(
                    term.scalar_spectator(),
                    &input.mass_squared,
                    exponent,
                    anchor,
                    settings
                        .rustred_scalar_numerator_limits
                        .diagnostics_enabled
                        .then_some(
                            settings
                                .rustred_scalar_numerator_limits
                                .max_scalar_product_degree,
                        ),
                );
                let (spectator, coefficient) =
                    if settings.rustred_application_options.reduction_strategy
                        == crate::RustRedReductionStrategy::Ordered
                    {
                        MatchedScalarInput::channel_coefficient(
                            spectator,
                            term.coefficient(),
                            family.coefficient_context(),
                            settings
                                .rustred_application_options
                                .reduction_limits
                                .exact_algebra,
                        )?
                    } else {
                        (spectator, term.coefficient().clone())
                    };
                batches
                    .entry((input.mass_squared.clone(), spectator, exponent))
                    .or_default()
                    .push((term.integral().clone(), coefficient));
                if let Some(progress) = &mut coverage_progress {
                    progress.completed(index + 1);
                }
            }
            log::debug!(target: "vakint::saved_application", "phase=coverage_mass_batching elapsed_seconds={:.6}", batching_started.elapsed().as_secs_f64());
            let lowering_guards_started = std::time::Instant::now();
            MonitorContext::native_phase(
                observation.as_deref_mut(),
                "lowering_guard_materialization",
                None,
            );
            conditions.extend(lowering.nonzero_conditions().iter().map(|condition| {
                condition
                    .replace(dimension.to_pattern())
                    .with(physical_dimension.to_pattern())
            }));
            log::debug!(target: "vakint::saved_application", "phase=lowering_guard_materialization elapsed_seconds={:.6}", lowering_guards_started.elapsed().as_secs_f64());
            if !matches!(input.mass_squared.as_view(), AtomView::Num(_)) {
                symbolic_masses.insert(input.mass_squared);
            }
            if let Some(progress) = &mut input_progress {
                progress.completed(input_index + 1);
            }
        }
        drop(scalar);
        if let Some(first) = uncovered.first() {
            log::warn!(
                "RustRed 5L: {} scalar key(s) produced by numerator lowering are outside the recorded shifted-index entry predicates (data {}); reduction will be attempted, but complete reduction to the known masters is not guaranteed. First common-family powers: {:?}",
                uncovered.len(),
                self.data_version,
                first.powers()
            );
            log::debug!(target: "vakint::five_loop_match", "phase=outside_coverage data={} keys={uncovered:?}; index convention: product_i D_i^(-a_i), D_i=q_i^2-M^2; negative a_i encode numerator powers, not the original physical dot-product rank", self.data_version);
        }
        let mut value = Atom::Zero;
        let mut families = BTreeMap::new();
        // Factorization is an existing complete checked action. Its internal
        // memo belongs to this ordered request, never the DFS session.
        let mut ordered_products = if settings.rustred_application_options.reduction_strategy
            == crate::RustRedReductionStrategy::Ordered
        {
            let _math = MathSpan::new();
            let anchor = application
                .output_family_shared(application.family().fingerprint())
                .ok_or_else(|| failure("saved application has no output family"))?;
            Some(
                VacuumProductReducer::try_new(
                    ShippedFactorProvider::with_data(Arc::clone(&self.data)),
                    anchor,
                    settings.rustred_product_limits,
                )
                .map_err(native_failure)?,
            )
        } else {
            None
        };
        let cohort_total = batches.len();
        let cohort_progress = progress.filter(|_| {
            settings.rustred_application_options.reduction_strategy
                == crate::RustRedReductionStrategy::Ordered
        });
        for (cohort_index, ((mass_squared, spectator, exponent), terms)) in
            batches.into_iter().enumerate()
        {
            if let Some(progress) = cohort_progress {
                progress.cohorts(cohort_index, cohort_total);
            }
            // Declare this before the returned native values: the complete
            // exact reduction/assembly interval includes their natural cleanup.
            let math = MathSpan::new();
            if settings.rustred_application_options.preparation_workers != 1 {
                math.invalidate("net timing requires one owner-preparation worker");
            }
            let spectator = spectator
                .replace(dimension.to_pattern())
                .with(physical_dimension.to_pattern());
            let mut cleanup = RetainedBatchCleanup::default();
            let application_started = std::time::Instant::now();
            MonitorContext::native_phase(
                observation.as_deref_mut(),
                "saved_rule_application",
                None,
            );
            let preparation_before = math
                .is_active()
                .then(|| application.owner_preparation_statistics().preparation_time);
            let reduced = self.reduce_saved(&mut application, terms, observation.as_deref_mut())?;
            if let Some(before) = preparation_before {
                // Native reports elapsed WALL time for each deferred owner.
                // The prepared-owner mutex and one-worker qualification keep
                // these intervals serial and enclosed by this math span.
                math.owner_preparation(
                    before,
                    application.owner_preparation_statistics().preparation_time,
                );
            }
            log::debug!(target: "vakint::saved_application", "phase=saved_rule_application elapsed_seconds={:.6}", application_started.elapsed().as_secs_f64());
            MonitorContext::native_phase(observation.as_deref_mut(), "terminal_geometry", None);
            let products_started = std::time::Instant::now();
            let mut contribution = Atom::Zero;
            let mut products = Vec::new();
            let mut product_conditions = reduced.nonzero_conditions().to_vec();
            for (terminal, coefficient) in reduced.terms() {
                match application.terminal_product(terminal) {
                    Some(product) => {
                        let factors = product
                            .factors()
                            .iter()
                            .map(|factor| {
                                Ok(VacuumFactor {
                                    family: application
                                        .output_family_shared(factor.family_fingerprint())
                                        .ok_or_else(|| {
                                            failure("product factor family is unavailable")
                                        })?,
                                    integral: factor.integral().clone(),
                                })
                            })
                            .collect::<Result<Vec<_>, RustRedEvaluationError>>()?;
                        let coefficient = application
                            .family()
                            .coefficient_context()
                            .try_mul(
                                coefficient,
                                product.coefficient(),
                                settings
                                    .rustred_application_options
                                    .reduction_limits
                                    .exact_algebra,
                            )
                            .map_err(native_failure)?;
                        products.push((factors, coefficient));
                    }
                    None => {
                        if let Some(chart) = application.terminal_product_chart(terminal) {
                            let projection = project_chart(
                                chart,
                                S.dot,
                                (1..=chart.combined_family().loop_count())
                                    .map(|axis| S.k.call(axis))
                                    .collect(),
                                ScalarNumeratorLimits {
                                    exact_algebra: settings
                                        .rustred_product_limits
                                        .reduction
                                        .exact_algebra,
                                    ..Default::default()
                                },
                                factor_projection::project,
                            )?;
                            product_conditions.extend(projection.nonzero_conditions);
                            for (factors, projected_coefficient) in projection.terms {
                                let coefficient = application
                                    .family()
                                    .coefficient_context()
                                    .try_mul(
                                        coefficient,
                                        &projected_coefficient,
                                        settings.rustred_product_limits.reduction.exact_algebra,
                                    )
                                    .map_err(native_failure)?;
                                products.push((factors, coefficient));
                            }
                            continue;
                        }
                        contribution += coefficient
                            .to_expression()
                            .replace(dimension.to_pattern())
                            .with(physical_dimension.to_pattern())
                            * Self::symbolic_integral(
                                terminal.family_fingerprint(),
                                terminal.integral(),
                                vk_symbol!("RustRedMaster"),
                                &mut families,
                            );
                    }
                }
            }
            MonitorContext::native_phase(
                observation.as_deref_mut(),
                "lower_factor_normalization_unobserved",
                None,
            );
            let normalized = match ordered_products.as_mut() {
                Some(reducer) => reducer.reduce_products(products, product_conditions),
                None => self
                    .product_reducer
                    .lock()
                    .map_err(failure)?
                    .reduce_products(products, product_conditions),
            };
            let normalized = match normalized {
                Ok(normalized) => normalized,
                Err(rustred::reduction::terminal_normalization::VacuumProductReductionError::UnresolvedFactor { family, target }) => {
                    log::debug!(target: "vakint::saved_application", "phase=terminal_geometry_and_lower_factors elapsed_seconds={:.6}", products_started.elapsed().as_secs_f64());
                    let fallback_started = std::time::Instant::now();
                    MonitorContext::native_phase(observation.as_deref_mut(), "unprojected_terminal_materialization", None);
                    log::warn!("RustRed lower-factor coverage is incomplete for {family} {:?}; retaining the exact parent terminals", target.powers());
                    // Preserve the original reduced sum and its complete guards;
                    // no failed or partially normalized product enters a cache.
                    let original = reduced.terms().iter().map(|(terminal, coefficient)| {
                        coefficient.to_expression().replace(dimension.to_pattern()).with(physical_dimension.to_pattern())
                            * Self::symbolic_integral(terminal.family_fingerprint(), terminal.integral(), vk_symbol!("RustRedMaster"), &mut families)
                    }).sum::<Atom>();
                    for condition in reduced.nonzero_conditions() {
                        let condition = condition.to_expression();
                        if condition.is_zero() || condition.get_all_symbols(true).iter().any(|symbol| !dimension.get_all_symbols(true).contains(symbol)) {
                            return Err(failure("saved application has an invalid non-generic-d guard").into());
                        }
                        conditions.insert(condition.replace(dimension.to_pattern()).with(physical_dimension.to_pattern()));
                    }
                    let mass_power = Atom::num((loop_count as i64, 2)) * &physical_dimension + Atom::num(exponent);
                    value += original * spectator * mass_squared.clone().pow(mass_power);
                    log::debug!(target: "vakint::saved_application", "phase=unprojected_terminal_materialization elapsed_seconds={:.6}", fallback_started.elapsed().as_secs_f64());
                    MonitorContext::native_phase(observation.as_deref_mut(), "retained_batch_cleanup", None);
                    cleanup.arm();
                    continue;
                }
                Err(error) => return Err(native_failure(error).into()),
            };
            log::debug!(target: "vakint::saved_application", "phase=terminal_geometry_and_lower_factors elapsed_seconds={:.6}", products_started.elapsed().as_secs_f64());
            let guards_started = std::time::Instant::now();
            MonitorContext::native_phase(
                observation.as_deref_mut(),
                "retained_guard_materialization",
                None,
            );
            let mut guard_progress = progress.map(|owner| {
                owner.phase(
                    "retained_guard_materialization",
                    normalized.nonzero_conditions().len(),
                    "guards processed",
                )
            });
            for (index, condition) in normalized.nonzero_conditions().iter().enumerate() {
                let condition = condition.to_expression();
                if condition.is_zero()
                    || condition
                        .get_all_symbols(true)
                        .iter()
                        .any(|symbol| !dimension.get_all_symbols(true).contains(symbol))
                {
                    return Err(failure(
                        "saved application has a non-generic-d guard; the symbolic Vakint result cannot discharge it",
                    ).into());
                }
                conditions.insert(
                    condition
                        .replace(dimension.to_pattern())
                        .with(physical_dimension.to_pattern()),
                );
                if let Some(progress) = &mut guard_progress {
                    progress.completed(index + 1);
                }
            }
            log::debug!(target: "vakint::saved_application", "phase=retained_guard_materialization elapsed_seconds={:.6}", guards_started.elapsed().as_secs_f64());
            let coefficients_started = std::time::Instant::now();
            MonitorContext::native_phase(
                observation.as_deref_mut(),
                "retained_coefficient_materialization",
                Some(contribution.as_view()),
            );
            let mut coefficient_progress = progress.map(|owner| {
                owner.phase(
                    "retained_coefficient_materialization",
                    normalized.terms().len(),
                    "master products materialized",
                )
            });
            for (index, (product, coefficient)) in normalized.terms().iter().enumerate() {
                let master = product
                    .iter()
                    .map(|factor| {
                        Self::symbolic_integral(
                            factor.family_fingerprint(),
                            factor.integral(),
                            vk_symbol!("RustRedMaster"),
                            &mut families,
                        )
                    })
                    .product::<Atom>();
                contribution += coefficient
                    .to_expression()
                    .replace(dimension.to_pattern())
                    .with(physical_dimension.to_pattern())
                    * master;
                if let Some(progress) = &mut coefficient_progress {
                    progress.completed_expression(
                        index + 1,
                        contribution.as_view(),
                        observation.as_deref_mut(),
                    );
                }
            }
            log::debug!(target: "vakint::saved_application", "phase=retained_coefficient_materialization elapsed_seconds={:.6}", coefficients_started.elapsed().as_secs_f64());
            let mass_started = std::time::Instant::now();
            MonitorContext::native_phase(
                observation.as_deref_mut(),
                "retained_batch_mass_assembly",
                Some(contribution.as_view()),
            );
            let mass_power =
                Atom::num((loop_count as i64, 2)) * &physical_dimension + Atom::num(exponent);
            value += contribution * spectator * mass_squared.clone().pow(mass_power);
            log::debug!(target: "vakint::saved_application", "phase=retained_batch_mass_assembly elapsed_seconds={:.6}", mass_started.elapsed().as_secs_f64());
            MonitorContext::native_phase(
                observation.as_deref_mut(),
                "retained_batch_cleanup",
                None,
            );
            cleanup.arm();
        }
        if let Some(progress) = cohort_progress {
            progress.cohorts(cohort_total, cohort_total);
        }
        conditions.extend(symbolic_masses);
        let normalization_started = std::time::Instant::now();
        MonitorContext::native_phase(
            observation.as_deref_mut(),
            "retained_normalization_and_notation",
            Some(value.as_view()),
        );
        let value = value * Self::normalization(settings, loop_count)?;
        let value = if settings.use_dot_product_notation {
            value
        } else {
            // All integrated scalar products have been eliminated. Apply the
            // requested presentation only to surviving external spectators.
            Vakint::convert_from_dot_notation(value.as_view())
        };
        log::debug!(target: "vakint::saved_application", "phase=retained_normalization_and_notation elapsed_seconds={:.6}", normalization_started.elapsed().as_secs_f64());
        let packing_started = std::time::Instant::now();
        MonitorContext::native_phase(
            observation.as_deref_mut(),
            "retained_detailed_result_packing",
            Some(value.as_view()),
        );
        let result = RustRedSymbolicEvaluation {
            value,
            nonzero_conditions: conditions.into_iter().collect(),
            families,
        };
        log::debug!(target: "vakint::saved_application", "phase=detailed_result_packing elapsed_seconds={:.6}", packing_started.elapsed().as_secs_f64());
        MonitorContext::native_phase(
            observation,
            "normalized_retained_output",
            Some(result.value.as_view()),
        );
        // The ordered product memo is request-owned; the DFS product memo
        // belongs to the engine and keeps its existing lifetime and policy.
        if let Some(products) = ordered_products {
            let _math = MathSpan::new();
            drop(products);
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vakint_macros::vk_parse;
    use crate::{EvaluationMethod, EvaluationOrder, RustRedEvaluationOptions};

    #[test]
    fn fmft_measure_conversion_cancels_exactly_in_each_epsilon_spelling() {
        for epsilon_symbol in ["ep", "eps", "eta"] {
            let settings = VakintSettings {
                epsilon_symbol: epsilon_symbol.into(),
                integral_normalization_factor: crate::LoopNormalizationFactor::FMFTandMATAD,
                ..Default::default()
            };
            for loops in 1..=5 {
                assert_eq!(
                    SavedApplication::normalization(&settings, loops).unwrap(),
                    Atom::one()
                );
            }
        }
    }

    #[test]
    fn chart_projector_resource_refusal_preserves_the_original_cause() {
        use rustred::algebra::CoefficientContext;
        use rustred::family::{AffineDenominator, IntegralFamily};
        use rustred::reduction::terminal_normalization::VacuumProductChartPlan;
        use symbolica::utils::resources::ResourceLimitError;

        let context = CoefficientContext::try_new(["d"]).unwrap();
        let family = IntegralFamily::new(
            "projector-resource-test",
            vec!["k1".into(), "k2".into()],
            vec![],
            context.clone(),
            context.parameter("d").unwrap(),
            [(-1, [1, 0, 0]), (-1, [0, 0, 1]), (0, [0, 1, 0])]
                .into_iter()
                .map(|(constant, row)| {
                    AffineDenominator::new(
                        context.integer(constant),
                        row.into_iter().map(|n| context.integer(n)).collect(),
                    )
                })
                .collect(),
            vec![],
            vec![context.zero(); 3],
        )
        .unwrap();
        let target = IntegralKey::try_new([1, 1, -1]).unwrap();
        let plan = VacuumProductChartPlan::prepare(
            &family,
            &BTreeSet::from([target.clone()]),
            Default::default(),
        )
        .unwrap();
        let chart = &plan.charts()[&target];
        let cause = ResourceLimitError::Memory {
            rss_bytes: 0,
            additional_bytes: 2,
            limit_bytes: 1,
        };
        let result = project_chart(
            chart,
            S.dot,
            vec![S.k.call(1), S.k.call(2)],
            Default::default(),
            |_| Err(RustRedEvaluationError::ResourceLimit(cause.clone())),
        );
        assert!(
            matches!(result, Err(RustRedEvaluationError::ResourceLimit(error)) if error == cause)
        );

        let healthy = project_chart(
            chart,
            S.dot,
            vec![S.k.call(1), S.k.call(2)],
            Default::default(),
            |_| {
                Ok(VacuumTensorProjection {
                    numerator: Atom::Zero,
                    nonzero_conditions: vec![context.parameter("d").unwrap().to_expression()],
                })
            },
        )
        .unwrap();
        assert!(healthy.terms.is_empty());
        assert!(healthy.nonzero_conditions.contains(family.dimension()));
    }

    #[test]
    fn saved_cache_policy_defaults_to_lru_and_participates_in_session_identity() {
        let defaults = VakintSettings::default();
        assert_eq!(
            defaults.rustred_scalar_numerator_limits,
            ScalarNumeratorLimits::default()
        );
        assert_eq!(
            defaults
                .rustred_scalar_numerator_limits
                .exact_algebra
                .max_term_operations,
            16_000_000
        );
        assert_eq!(
            defaults.rustred_application_options.cache_policy,
            crate::RustRedCachePolicy::EvictLeastRecentlyUsed
        );
        let mut fail = defaults.rustred_application_options;
        fail.cache_policy = crate::RustRedCachePolicy::FailOnLimit;
        assert_ne!(fail, defaults.rustred_application_options);
        assert!(matches!(
            Vakint::new()
                .unwrap()
                .rustred_application_cache_work_statistics(),
            Err(VakintError::RustRedEvaluation(
                RustRedEvaluationError::SavedApplicationUnavailable
            ))
        ));
    }

    fn root() -> Atom {
        vk_parse!("topo(I5L_32745(muvsq,1,1,1,1,1,1,1,1,1,1,1,1))").unwrap()
    }

    fn settings() -> VakintSettings {
        VakintSettings {
            form_exe_path: "/five-loop-must-not-use-form".into(),
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                RustRedEvaluationOptions {
                    master_mode: crate::RustRedMasterMode::Retain,
                },
            )]),
            ..Default::default()
        }
    }

    #[test]
    fn absent_artifact_path_does_not_enable_five_loop_application() {
        let vakint = Vakint::new().unwrap();
        let mut settings = settings();
        settings.allow_unknown_integrals = false;
        assert!(settings.rustred_artifact_path.is_none());
        assert!(vakint.topologies.five_loop_catalogue_statistics().is_none());
        let error = vakint
            .evaluate_integral(&settings, root().as_view())
            .unwrap_err();
        assert!(error.to_string().contains("5L RustRed data is not loaded"));
        assert!(matches!(
            vakint.evaluate_rustred_symbolic(&settings, root().as_view()),
            Err(VakintError::RustRedEvaluation(
                RustRedEvaluationError::SavedApplicationUnavailable
            ))
        ));
    }

    #[test]
    #[cfg(feature = "experimental-rustred")]
    fn raw_observer_preserves_public_admission_errors_and_isolates_faults() {
        use rustred::solver::{
            CandidateRoutedObservationOptions, CandidateRoutedObserverFault,
            CandidateRoutedObserverStatus,
        };
        let vakint = Vakint::new().unwrap();
        let settings = settings();
        let input = root();
        let ordinary = vakint
            .evaluate_rustred_symbolic(&settings, input.as_view())
            .unwrap_err();
        let mut callback = |_: &rustred::solver::CandidateRoutedApplicationEnvelope<'_>| Err(());
        let mut observer = CandidateRoutedApplicationObserver::new(
            CandidateRoutedObservationOptions::default(),
            &mut callback,
        );
        let observed = vakint
            .evaluate_rustred_symbolic_with_observer(&settings, input.as_view(), &mut observer)
            .unwrap_err();
        assert_eq!(ordinary.to_string(), observed.to_string());
        assert!(matches!(
            observed,
            VakintError::RustRedEvaluation(RustRedEvaluationError::SavedApplicationUnavailable)
        ));
        assert_eq!(
            observer.status(),
            CandidateRoutedObserverStatus::Faulted(CandidateRoutedObserverFault::CallbackError)
        );
    }

    #[test]
    fn wrong_artifact_path_fails_before_evaluation_without_form() {
        let vakint = Vakint::new().unwrap();
        let mut settings = settings();
        settings.rustred_artifact_path =
            Some("/nonexistent-vakint-five-loop-application/application.json".into());
        let error = vakint
            .evaluate_integral(&settings, root().as_view())
            .unwrap_err();
        assert!(matches!(
            error,
            VakintError::RustRedEvaluation(RustRedEvaluationError::Reduction { .. })
        ));
        assert!(vakint.topologies.five_loop_catalogue_statistics().is_none());
        // A numerical request must still authenticate its data; it no longer
        // fails early merely because the matched input has five loops.
        settings.evaluation_order = EvaluationOrder::rustred_only();
        let error = vakint
            .evaluate_integral(&settings, root().as_view())
            .unwrap_err();
        assert!(matches!(
            error,
            VakintError::RustRedEvaluation(RustRedEvaluationError::Reduction { .. })
        ));
        assert!(vakint.topologies.five_loop_catalogue_statistics().is_none());
    }

    #[test_log::test]
    #[ignore = "requires VAKINT_SAVED_APPLICATION pointing to an external generation-54 package"]
    fn external_package_bounded_lru_cache_preserves_admitted_repeats() {
        // A measured 1 MiB allowance cannot hold this root's protected live
        // children. 4 MiB still forces eviction below its 7.8 MB full cache.
        const CACHE_BYTES: usize = 4 * 1024 * 1024;
        let path = std::env::var_os("VAKINT_SAVED_APPLICATION")
            .expect("set VAKINT_SAVED_APPLICATION to the uncommitted external package");
        let mut settings = settings();
        settings.rustred_application_options.cache_policy =
            crate::RustRedCachePolicy::EvictLeastRecentlyUsed;
        settings
            .rustred_application_options
            .reduction_limits
            .max_cached_coefficient_bytes = CACHE_BYTES;
        settings
            .rustred_application_options
            .reduction_limits
            .max_cached_coefficient_terms = usize::MAX;
        let vakint = Vakint::new()
            .unwrap()
            .with_rustred_application(&mut settings, PathBuf::from(path))
            .unwrap();
        let application = vakint.loaded_application().unwrap();
        let mut controls = Vec::new();
        for powers in [
            // Generation-54 points162 and163. The second requires more than
            // this retained allowance with fail-on-limit memoization.
            [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 5],
            [0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 4],
        ] {
            let key = IntegralKey::try_new(powers).unwrap();
            assert!(
                application
                    .prepared
                    .lock()
                    .unwrap()
                    .routed_entry_is_admitted(&key)
                    .unwrap()
            );
            let input = super::admission_census_tests::input_for_key(&key)
                .unwrap()
                .expression();
            let first = vakint
                .evaluate_rustred_symbolic(&settings, input.as_view())
                .unwrap();
            assert!(
                vakint
                    .rustred_application_cache_statistics()
                    .unwrap()
                    .cached_coefficient_bytes()
                    <= CACHE_BYTES
            );
            controls.push((input, first));
        }
        assert!(
            controls.len() > 1,
            "the pinned package must provide distinct admitted cache controls"
        );
        for (input, first) in &controls {
            let repeated = vakint
                .evaluate_rustred_symbolic(&settings, input.as_view())
                .unwrap();
            assert_eq!(first.value, repeated.value);
            assert_eq!(first.nonzero_conditions, repeated.nonzero_conditions);
            assert_eq!(first.families, repeated.families);
            assert!(
                vakint
                    .rustred_application_cache_statistics()
                    .unwrap()
                    .cached_coefficient_bytes()
                    <= CACHE_BYTES
            );
        }
        let work = vakint.rustred_application_cache_work_statistics().unwrap();
        assert!(
            work.evicted_entries > 0,
            "the real workload must actually exercise eviction"
        );
        assert!(work.evicted_coefficient_bytes > 0);
        // A warmed rule cache must not bypass a tightened caller-numerator
        // policy, in either witness selection or the final weighted handoff.
        let mut strict = settings.clone();
        strict.rustred_scalar_numerator_limits.max_input_nodes = 0;
        strict.rustred_scalar_numerator_limits =
            strict.rustred_scalar_numerator_limits.with_diagnostics();
        assert!(Arc::ptr_eq(
            &application,
            &vakint.saved_application(&strict).unwrap()
        ));
        let expression = VakintExpression::try_from(controls[0].0.as_view()).unwrap();
        let mut matched = vakint
            .topologies
            .match_topologies_to_user_input(expression.0[0].integral.as_view(), false)
            .unwrap()
            .unwrap();
        let expected = rustred::scalar_numerator::ScalarNumeratorError::ResourceLimit {
            resource: "weighted lowering input nodes",
            requested: 1,
            limit: 0,
        }
        .to_string();
        let assert_scalar_limit = |error: VakintError| match error {
            VakintError::RustRedEvaluation(RustRedEvaluationError::Reduction { detail }) => {
                assert_eq!(detail, expected)
            }
            other => panic!("unexpected tightened scalar-policy failure: {other}"),
        };
        let before = vakint.rustred_application_cache_statistics().unwrap();
        assert_scalar_limit(
            vakint
                .rustred_matched_input_key(&strict, expression.0[0].numerator.as_view(), &matched)
                .unwrap_err(),
        );
        matched.apply_replacement_rules().unwrap();
        let mut term = expression.0[0].clone();
        term.apply_scalar_numerator_replacement_rules(&matched, &strict)
            .unwrap();
        assert_scalar_limit(
            application
                .evaluate_detailed(&strict, term.numerator.as_view(), &matched, None, None)
                .unwrap_err(),
        );
        let after = vakint.rustred_application_cache_statistics().unwrap();
        assert_eq!(before.cached_integrals(), after.cached_integrals());
        assert_eq!(
            before.cached_coefficient_bytes(),
            after.cached_coefficient_bytes()
        );
        assert_eq!(
            work,
            vakint.rustred_application_cache_work_statistics().unwrap()
        );
        assert!(
            application
                .prepared
                .lock()
                .unwrap()
                .routed_entry_is_admitted(
                    &IntegralKey::try_new([0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 5]).unwrap()
                )
                .unwrap()
        );
        println!(
            "tightened_scalar_policy_before_warm_cache=refused_without_reload_or_scope_change"
        );
        println!("admitted_lru_controls={} work={work:?}", controls.len());
        vakint.clear_rustred_application_cache().unwrap();
        assert_eq!(
            vakint
                .rustred_application_cache_statistics()
                .unwrap()
                .cached_integrals(),
            0
        );
        assert_eq!(
            vakint.rustred_application_cache_work_statistics().unwrap(),
            work
        );
    }

    #[test_log::test]
    #[ignore = "requires VAKINT_SAVED_APPLICATION pointing to <data root>/5L/ibps of the external 640-terminal package"]
    fn external_package_all_factorized_terminals_reduce_through_existing_lower_programs() {
        use rustred::reduction::terminal_normalization::VacuumProductProjection;
        let path = PathBuf::from(
            std::env::var_os("VAKINT_SAVED_APPLICATION")
                .expect("set VAKINT_SAVED_APPLICATION to the uncommitted external package"),
        );
        let mut settings = settings();
        let vakint = Vakint::new().unwrap();
        let before_load_clone = vakint.clone();
        assert!(vakint.topologies.five_loop_catalogue_statistics().is_none());
        assert!(
            before_load_clone
                .topologies
                .five_loop_catalogue_statistics()
                .is_none()
        );
        let vakint = vakint
            .with_rustred_application(&mut settings, path)
            .unwrap();
        // Activating a complete release replaces this instance's registry and
        // application cache. Earlier clones retain their previous snapshot.
        assert!(before_load_clone.loaded_application().is_err());
        assert!(
            before_load_clone
                .topologies
                .five_loop_catalogue_statistics()
                .is_none()
        );
        assert!(
            !before_load_clone
                .rustred_data()
                .unwrap()
                .loops
                .contains_key(&5)
        );
        let after_load_clone = vakint.clone();
        let data = vakint.rustred_data().unwrap();
        assert!((1..=5).all(|loops| data.loops.contains_key(&loops)));
        assert!(data.four_loop.is_some());
        let application = vakint.loaded_application().unwrap();
        let cloned_application = after_load_clone.saved_application(&settings).unwrap();
        assert!(Arc::ptr_eq(&application, &cloned_application));
        assert!(Arc::ptr_eq(
            vakint.topologies.five_loop_catalogue().unwrap(),
            after_load_clone.topologies.five_loop_catalogue().unwrap(),
        ));
        assert_eq!(
            vakint
                .topologies
                .five_loop_catalogue_statistics()
                .unwrap()
                .topology_count,
            67
        );
        let cache_work = vakint.rustred_application_cache_work_statistics().unwrap();
        assert_eq!(cache_work, crate::RustRedCacheWorkStatistics::default());
        vakint.clear_rustred_application_cache().unwrap();
        assert_eq!(
            vakint.rustred_application_cache_work_statistics().unwrap(),
            cache_work
        );
        let prepared = application.prepared.lock().unwrap();
        let mut reducer = application.product_reducer.lock().unwrap();
        let mut projected_count = 0;
        let mut cross_count = 0;
        let mut positive_count = 0;
        let mut positive_supports = BTreeSet::new();
        for terminal in prepared.remaining_terminals() {
            println!(
                "factor_terminal_candidate={:?}",
                terminal.integral().powers()
            );
            let mut crosses = false;
            let projection = if let Some(chart) = prepared.terminal_product_chart(terminal) {
                projected_count += 1;
                project_chart(
                    chart,
                    S.dot,
                    (1..=chart.combined_family().loop_count())
                        .map(|axis| S.k.call(axis))
                        .collect(),
                    Default::default(),
                    |input| {
                        crosses = true;
                        factor_projection::project(input)
                    },
                )
                .unwrap()
            } else if let Some(product) = prepared.terminal_product(terminal) {
                positive_count += 1;
                positive_supports.insert(
                    terminal
                        .integral()
                        .powers()
                        .iter()
                        .map(|power| *power > 0)
                        .collect::<Vec<_>>(),
                );
                VacuumProductProjection {
                    terms: vec![(
                        product
                            .factors()
                            .iter()
                            .map(|factor| VacuumFactor {
                                family: prepared
                                    .output_family_shared(factor.family_fingerprint())
                                    .unwrap(),
                                integral: factor.integral().clone(),
                            })
                            .collect(),
                        product.coefficient().clone(),
                    )],
                    nonzero_conditions: Vec::new(),
                }
            } else {
                continue;
            };
            let first = reducer
                .reduce_products(
                    projection.terms.clone(),
                    projection.nonzero_conditions.clone(),
                )
                .unwrap();
            let warm = reducer
                .reduce_products(projection.terms, projection.nonzero_conditions)
                .unwrap();
            assert_eq!(first.terms(), warm.terms());
            assert_eq!(first.nonzero_conditions(), warm.nonzero_conditions());
            assert!(
                first
                    .terms()
                    .keys()
                    .flatten()
                    .all(|factor| { factor.family_fingerprint() != terminal.family_fingerprint() })
            );
            println!(
                "cross_projected={crosses} product_terms={} nonzero_conditions={} cache={:?}",
                first.terms().len(),
                first.nonzero_conditions().len(),
                reducer.cache_statistics()
            );
            reducer.clear_cache().unwrap();
            cross_count += usize::from(crosses);
        }
        // Pinned to the published 640-terminal package (generation 54 had
        // 38 scalar products, 14 numerator charts and 3 cross projections).
        assert_eq!(positive_count, 25, "pinned 640-pack scalar products");
        assert_eq!(
            positive_supports.len(),
            19,
            "all factorized support classes"
        );
        assert_eq!(
            projected_count, 7,
            "pinned 640-pack numerator product inventory"
        );
        assert_eq!(cross_count, 2, "native cross-ISP transport inventory");
    }
}
