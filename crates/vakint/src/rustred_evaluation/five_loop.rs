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
use super::parallel;
use super::request_dump::{self, DumpChannel, DumpProgram, DumpRequest, RequestDump};
use rustred::algebra::Coefficient;
use rustred::family::IntegralKey;
use rustred::input::{Compiler, Limits, LoweringLimits, TextProject, TextPropagator};
use rustred::reduction::terminal_normalization::{
    VacuumFactor, VacuumProductChartLimits, VacuumProductProjection, VacuumProductReducer,
    VacuumProductReductionLimits, VacuumTensorProjection, VacuumTensorProjectionInput,
    VerifiedVacuumProductChart,
};
use rustred::scalar_numerator::{FamilyScalarNumeratorService, ScalarNumeratorLimits};
use rustred::solver::CandidateRoutedApplicationObserver;
use rustred_app::{
    PreparedSavedApplication, SAVED_ROUTE_PREPARATION_SCHEMA, SavedApplicationOptions,
    SavedApplicationProductChartStatistics, SavedApplicationProductStatistics, SavedRouteAddon,
    SavedRoutePreparationSource, load_saved_application_with_route_addon,
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
    ReplacementRules, Topology, Vakint, VakintError, VakintExpression, VakintSettings, VakintTerm,
    symbols::S,
};

static PARENTS: LazyLock<Vec<ParentDescriptor>> = LazyLock::new(descriptors);
/// Guards materialized per parallel chunk between progress updates.
const GUARD_CHUNK: usize = 4096;

/// Worker threads for the independent five-loop frontend items (witness
/// screening, full witness keys, and per-contribution scalar lowering):
/// [`VakintSettings::n_cores`]. `1` is strictly sequential on the calling
/// thread, with no pool; results never depend on this value.
fn frontend_threads(settings: &VakintSettings) -> usize {
    settings.n_cores.max(1)
}

/// A selected five-loop witness and, when a selector ran, whether it proved
/// the witness admitted (with its exact lowered key support) or not.
type SelectedWitness = (ReplacementRules, Option<Option<Vec<IntegralKey>>>);

#[path = "five_loop/master_transport.rs"]
mod master_transport;
#[path = "five_loop/witness_cache.rs"]
mod witness_cache;
pub(crate) use master_transport::rustred_evaluates_loop;
#[path = "five_loop/witness_screen.rs"]
mod witness_screen;
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
    /// Canonical package directory, recorded by native request dumps.
    path: PathBuf,
    fingerprint: String,
    data_version: String,
    products: SavedApplicationProductStatistics,
    product_charts: SavedApplicationProductChartStatistics,
    data: Arc<super::data::RustRedData>,
    witnesses: Mutex<WitnessCache>,
    /// The package masters in their own contracted long forms, prepared on
    /// the first master-transport attempt (see
    /// [`Vakint::transport_rustred_master_term`]).
    master_forms: std::sync::OnceLock<Vec<master_transport::MasterForm>>,
    /// Contracted scalar inputs (digests) the master transport found not to
    /// be masters, so a repeated input skips its map enumeration; cleared with
    /// the witness memo.
    master_transport_misses: Mutex<BTreeSet<[u8; 32]>>,
    /// Ordered step-cache identity of the loaded package: its
    /// saved-application identity as rule tables and RustRed's binding
    /// (family, variables, policy).
    step_cache: super::cache::OrderedCacheProgram,
    /// A request stopped by a cancellation or a resource limit left a failed
    /// one-time owner preparation behind: the native application keeps that
    /// failure for good, so the next request reloads the package.
    poisoned: std::sync::atomic::AtomicBool,
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

/// Whether `error` reports an engine time or memory limit (also when RustRed
/// carried it as text).
fn stopped_by_resource_limit(error: &VakintError) -> bool {
    if matches!(
        error,
        VakintError::ResourceLimit(_)
            | VakintError::RustRedEvaluation(RustRedEvaluationError::ResourceLimit(_))
    ) {
        return true;
    }
    let text = error.to_string();
    // Symbolica's `ResourceLimitError` texts.
    [
        "exceeds RAM budget",
        "configured time budget",
        "cannot enforce ",
    ]
    .iter()
    .any(|marker| text.contains(marker))
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

/// The registered complete common momentum basis as a RustRed family.
fn common_basis_family(
    name: &str,
) -> Result<rustred::family::IntegralFamily, RustRedEvaluationError> {
    Ok(Compiler::new(Limits::default())
        .map_err(native_failure)?
        .compile_text(TextProject {
            name: Some(name.to_owned()),
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
        .into_family())
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
            return Ok(self
                .match_rustred_input_selected(settings, integral, numerator)?
                .map(|(witness, _)| witness));
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
            match self.rustred_matched_input_support(settings, numerator, &cached.witness)? {
                Some(support) if support == cached.accepted => {
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
        let Some((witness, support)) =
            self.match_rustred_input_selected(settings, integral, numerator)?
        else {
            return Ok(None);
        };
        // A verified fallback outside recorded coverage remains usable by
        // the partial reducer. Only covered witnesses enter this memo, with
        // the exact support their selection proved admitted.
        if matches!(witness.canonical_topology, Topology::FiveLoop(_))
            && let Some(support) = support.flatten()
            && let Ok(mut cache) = application.witnesses.lock()
        {
            cache.insert(registry, policy, integral, numerator, &witness, support);
        }
        Ok(Some(witness))
    }

    #[cfg(test)]
    fn match_rustred_input_uncached(
        &self,
        settings: &VakintSettings,
        integral: AtomView,
        numerator: AtomView,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        Ok(self
            .match_rustred_input_selected(settings, integral, numerator)?
            .map(|(witness, _)| witness))
    }

    /// Select the verified witness whose exact surviving common-family
    /// integrand is minimal (see [`Self::select_rustred_witness`]). The second
    /// value is `Some(Some(support))` when the selector proved the returned
    /// witness admitted with that exact lowered key support, `Some(None)` when
    /// it proved it outside the saved entry union, and `None` when no
    /// selector ran (no five-loop witness). The selector also runs on a single
    /// parametric/short-form five-loop fallback, so its support is recorded.
    fn match_rustred_input_selected(
        &self,
        settings: &VakintSettings,
        integral: AtomView,
        numerator: AtomView,
    ) -> Result<Option<SelectedWitness>, VakintError> {
        #[cfg(all(
            test,
            feature = "weighted-lowering-diagnostics",
            feature = "feynkit-ingress"
        ))]
        let diagnostic_thread = weighted_diagnostic::begin_search(integral, numerator)?;
        let mut selection = None;
        let result = self
            .topologies
            .match_topologies_to_user_input_with_selection(
                integral,
                settings.allow_unknown_integrals,
                settings.rustred_momentum_match_limits,
                |witnesses| {
                    #[cfg(all(
                        test,
                        feature = "weighted-lowering-diagnostics",
                        feature = "feynkit-ingress"
                    ))]
                    weighted_diagnostic::assert_thread(diagnostic_thread)?;
                    let (index, support) =
                        self.select_rustred_witness(settings, numerator, witnesses)?;
                    selection = Some(support);
                    Ok(index)
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
        let Some(witness) = result? else {
            return Ok(None);
        };
        // The selector also screens the single short-form fallback witness;
        // either way its support belongs to the returned witness.
        Ok(Some((witness, selection)))
    }

    /// Exact selection among verified witnesses, identical to taking the first
    /// witness (in enumeration order) with the minimal admitted key of
    /// [`Self::rustred_matched_input_key`], or `None` if no witness is admitted,
    /// except that admitted witnesses whose whole exact support consists of
    /// the package's remaining terminals (masters) are preferred over all
    /// others. Every witness is the same integral; a terminal-only witness
    /// needs no IBP walk, while the lexicographically minimal one can be a
    /// symmetric image that the package reduces by a long walk (the 5L master
    /// with powers [1,1,0,1,2,-1,1,1,1,1,0,1,0,0,1] was sent to its image
    /// [-1,1,0,1,1,1,1,1,1,1,0,2,0,0,1]: an Ordered walk of 1.57 M keys and a
    /// modular lift of hours, instead of the identity).
    ///
    /// The full key needs a weighted lowering over symbolic coefficients per
    /// witness. Instead, admission is first decided from the exact lowered key
    /// support computed from rational groups ([`witness_screen`]) where the
    /// numerator allows it: a certainly nonzero key outside the saved entry
    /// union proves rejection, as the full key rejects exactly such witnesses
    /// (screening stops at the first one). Among admitted witnesses,
    /// keys are compared lexicographically and the first entries' integral keys
    /// are distinct whenever their smallest support keys are, so only witnesses
    /// sharing the smallest first key can be minimal. A unique one is selected
    /// without its full key; otherwise the full keys of exactly those witnesses
    /// decide as before. Witnesses outside the grouped grammar always use the
    /// full key. Returns the selected index and its exact support, if admitted.
    ///
    /// Behavioural difference, errors only: a witness that this exact
    /// screening rejects or outranks is never fully lowered, so an error that
    /// only its full lowering would raise (e.g. an exact-algebra size limit)
    /// no longer aborts the match, and admission is probed in discovery rather
    /// than key order. Count diagnostics disable the screening.
    fn select_rustred_witness(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        witnesses: &[ReplacementRules],
    ) -> Result<(Option<usize>, Option<Vec<IntegralKey>>), VakintError> {
        let started = std::time::Instant::now();
        let threads = frontend_threads(settings);
        let application = self.saved_application(settings)?;
        let prepared = application.prepared.lock().map_err(failure)?;
        let scalar = FamilyScalarNumeratorService::try_new(
            prepared.family(),
            S.dot,
            (1..=5).map(|axis| S.k.call(axis)).collect(),
            settings.rustred_scalar_numerator_limits,
        )
        .map_err(native_failure)?;
        let prepared = &*prepared;
        let groups = if settings.verify_numerator_identification
            || settings.rustred_scalar_numerator_limits.diagnostics_enabled
        {
            None
        } else {
            // The full key starts with this witness-independent scalar form.
            let scalar_form = VakintTerm::scalar_product_form(settings, numerator)?;
            witness_screen::NumeratorGroups::new(scalar_form.as_view())
        };
        log::debug!(target: "vakint::five_loop_match", "phase=witness_screen_groups witnesses={} threads={threads} groups={:?}", witnesses.len(), groups.as_ref().map(witness_screen::NumeratorGroups::shape));
        // Admission of a key is a pure function of the locked package; share
        // it between witnesses. Errors are never cached.
        let admission = Mutex::new(BTreeMap::<IntegralKey, bool>::new());
        let admitted = |key: &IntegralKey| -> Result<bool, VakintError> {
            if let Some(known) = admission.lock().map_err(failure)?.get(key) {
                return Ok(*known);
            }
            let known = prepared
                .routed_entry_is_admitted(key)
                .map_err(native_failure)?;
            admission
                .lock()
                .map_err(failure)?
                .insert(key.clone(), known);
            Ok(known)
        };
        enum Status {
            Rejected,
            Admitted(Vec<IntegralKey>),
            Unknown,
        }
        let status = match &groups {
            Some(groups) => super::parallel::map_ordered(threads, witnesses, |_, witness| {
                groups.screen(settings, witness, &scalar, &admitted)
            }),
            None => witnesses.iter().map(|_| Ok(None)).collect(),
        }
        .into_iter()
        .map(|screen| {
            Ok(match screen? {
                None => Status::Unknown,
                Some(witness_screen::Screen::Rejected) => Status::Rejected,
                Some(witness_screen::Screen::Admitted(support)) => Status::Admitted(support),
            })
        })
        .collect::<Result<Vec<_>, VakintError>>()?;
        let full_key = |index: usize| {
            self.matched_input_key_with(settings, numerator, &witnesses[index], prepared, &scalar)
        };
        let unknown = status
            .iter()
            .enumerate()
            .filter_map(|(index, status)| matches!(status, Status::Unknown).then_some(index))
            .collect::<Vec<_>>();
        let mut keys = BTreeMap::new();
        for (index, key) in unknown.iter().zip(super::parallel::map_ordered(
            threads,
            &unknown,
            |_, index| full_key(*index),
        )) {
            keys.insert(*index, key?);
        }
        // Smallest first integral key of every admitted witness (`None` for
        // an empty support, which precedes every nonempty key).
        let heads = status
            .iter()
            .enumerate()
            .filter_map(|(index, status)| match status {
                Status::Rejected => None,
                Status::Admitted(support) => Some((index, support.first().cloned())),
                Status::Unknown => keys[&index]
                    .as_ref()
                    .map(|key| (index, key.first().map(|(key, _)| key.clone()))),
            })
            .collect::<Vec<_>>();
        let support_of =
            |index: usize, keys: &BTreeMap<usize, Option<Vec<(IntegralKey, Atom)>>>| match &status
                [index]
            {
                Status::Admitted(support) => support.clone(),
                _ => keys[&index]
                    .as_ref()
                    .expect("admitted full key")
                    .iter()
                    .map(|(key, _)| key.clone())
                    .collect(),
            };
        // Prefer witnesses whose exact support is all remaining terminals: the
        // input is then already in masters and its reduction is the identity.
        // An empty support (a vanishing integrand) qualifies vacuously and still
        // precedes every nonempty key below.
        let terminals = prepared
            .remaining_terminals()
            .iter()
            .filter(|terminal| terminal.family_fingerprint() == prepared.family().fingerprint())
            .map(|terminal| terminal.integral().clone())
            .collect::<BTreeSet<IntegralKey>>();
        let terminal_only = heads
            .iter()
            .filter(|(index, _)| {
                support_of(*index, &keys)
                    .iter()
                    .all(|key| terminals.contains(key))
            })
            .cloned()
            .collect::<Vec<_>>();
        let preferred_terminals = !terminal_only.is_empty();
        let heads = if preferred_terminals {
            terminal_only
        } else {
            heads
        };
        let Some(minimum) = heads.iter().map(|(_, head)| head).min().cloned() else {
            log::debug!(target: "vakint::five_loop_match", "phase=witness_selection admitted=0 screened={} full_keys={} elapsed_seconds={:.6}", witnesses.len() - unknown.len(), keys.len(), started.elapsed().as_secs_f64());
            return Ok((None, None));
        };
        let tied = heads
            .iter()
            .filter(|(_, head)| *head == minimum)
            .map(|(index, _)| *index)
            .collect::<Vec<_>>();
        let selected = if let [only] = tied.as_slice() {
            *only
        } else {
            let missing = tied
                .iter()
                .copied()
                .filter(|index| !keys.contains_key(index))
                .collect::<Vec<_>>();
            for (index, key) in missing.iter().zip(super::parallel::map_ordered(
                threads,
                &missing,
                |_, index| full_key(*index),
            )) {
                keys.insert(*index, key?);
            }
            // The first minimal full key in enumeration order, as before.
            let mut best: Option<(&Vec<(IntegralKey, Atom)>, usize)> = None;
            for index in &tied {
                let key = keys[index].as_ref();
                // Screened witnesses are cross-checked against their full key.
                if let Status::Admitted(support) = &status[*index]
                    && key.is_none_or(|key| !key.iter().map(|(key, _)| key).eq(support.iter()))
                {
                    return Err(failure(
                        "exact witness screening disagrees with the full weighted lowering",
                    )
                    .into());
                }
                let key = key.ok_or_else(|| failure("tied witness has no admitted key"))?;
                if best.is_none_or(|(old, _)| key < old) {
                    best = Some((key, *index));
                }
            }
            best.expect("tied witnesses are nonempty").1
        };
        log::debug!(target: "vakint::five_loop_match", "phase=witness_selection admitted={} tied={} selected={selected} terminal_only={preferred_terminals} screened={} full_keys={} elapsed_seconds={:.6}", heads.len(), tied.len(), witnesses.len() - unknown.len(), keys.len(), started.elapsed().as_secs_f64());
        Ok((Some(selected), Some(support_of(selected, &keys))))
    }

    /// The exact lowered key support of an admitted witness, or `None` if it
    /// is outside the saved entry union. This is the memo's admission token.
    fn rustred_matched_input_support(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        matched: &ReplacementRules,
    ) -> Result<Option<Vec<IntegralKey>>, VakintError> {
        Ok(self
            .select_rustred_witness(settings, numerator, std::slice::from_ref(matched))?
            .1)
    }

    /// Order already verified matcher witnesses by their exact surviving
    /// common-family integrand. This checks the saved entry union without
    /// applying any rule or preparing an owner. Selection computes the same
    /// key where it needs one; tests and diagnostics call this directly.
    #[cfg(test)]
    pub(crate) fn rustred_matched_input_key(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        matched: &ReplacementRules,
    ) -> Result<Option<Vec<(IntegralKey, Atom)>>, VakintError> {
        let application = self.saved_application(settings)?;
        let prepared = application.prepared.lock().map_err(failure)?;
        let scalar = FamilyScalarNumeratorService::try_new(
            prepared.family(),
            S.dot,
            (1..=5).map(|axis| S.k.call(axis)).collect(),
            settings.rustred_scalar_numerator_limits,
        )
        .map_err(native_failure)?;
        self.matched_input_key_with(settings, numerator, matched, &prepared, &scalar)
    }

    fn matched_input_key_with(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        matched: &ReplacementRules,
        prepared: &PreparedSavedApplication<15>,
        scalar: &FamilyScalarNumeratorService<'_>,
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
            self.load_external_data_on(root, 5, configured.n_cores)?;
        }
        self.saved_application(&configured)?;
        *settings = configured;
        Ok(self)
    }

    pub(crate) fn saved_application(
        &self,
        settings: &VakintSettings,
    ) -> Result<Arc<SavedApplication>, VakintError> {
        /// Route preparation add-on choice of one native load.
        enum RouteAddonChoice {
            /// The official 5L pack declares it and its sha256 matched.
            Declared(PathBuf),
            /// Declared but unusable here (missing, unlisted, sha256
            /// mismatch, other format; reported on stderr): verify every
            /// route.
            Rejected,
            /// Not declared: RustRed's default (a `route-preparation.rrbin`
            /// next to `application.json`, when present).
            Undeclared,
        }
        /// The add-on applies only to the official pack's own `ibps/` saved
        /// application; its file is checked against the pack manifest here
        /// (the data load skips it, like a shipped Ordered cache). A refusal
        /// prints one line on stderr, as RustRed does for an add-on it
        /// rejects; resource limits stopping the check are the load's error.
        fn route_preparation_addon(
            data: &super::data::RustRedData,
            official: Option<&Path>,
            path: &Path,
        ) -> Result<RouteAddonChoice, RustRedEvaluationError> {
            let Some(pack) = data.loops.get(&5) else {
                return Ok(RouteAddonChoice::Undeclared);
            };
            let Some(declared) = pack.manifest.route_preparation.as_ref() else {
                return Ok(RouteAddonChoice::Undeclared);
            };
            if official
                .and_then(|official| official.canonicalize().ok())
                .as_deref()
                != Some(path)
            {
                return Ok(RouteAddonChoice::Undeclared);
            }
            let reason = if declared.format != SAVED_ROUTE_PREPARATION_SCHEMA {
                format!(
                    "format {} is not {SAVED_ROUTE_PREPARATION_SCHEMA}",
                    declared.format
                )
            } else {
                match pack.read(&declared.path) {
                    Err(error @ RustRedEvaluationError::ResourceLimit(_)) => return Err(error),
                    Err(error) => error.to_string(),
                    Ok(_) => match pack.path(&declared.path) {
                        Ok(Some(addon)) => return Ok(RouteAddonChoice::Declared(addon)),
                        Ok(None) => "the pack has no file system path".to_owned(),
                        Err(error) => error.to_string(),
                    },
                }
            };
            let shown = pack
                .path(&declared.path)
                .ok()
                .flatten()
                .unwrap_or_else(|| PathBuf::from(format!("5L/{}", declared.path)));
            eprintln!(
                "vakint: 5L route preparation add-on {} not used ({reason}); verifying every \
                 route",
                shown.display()
            );
            Ok(RouteAddonChoice::Rejected)
        }

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
        let options = super::application_options(settings);
        // Requests select their weight arithmetic and scheduler under the
        // prepared-owner lock before reducing, and preparation workers only
        // parallelize route verification during the load: none of them
        // changes the prepared application, so none may reload and re-prepare
        // the package.
        let prepared_identity = |options: SavedApplicationOptions| SavedApplicationOptions {
            ordered_weight_arithmetic: Default::default(),
            ordered_scheduler: Default::default(),
            preparation_workers: 1,
            ..options
        };
        let mut loaded = self.rustred_application.lock().map_err(failure)?;
        // A package poisoned by a stopped request is dropped before its
        // replacement loads (its memory is not held twice); the catalogue is
        // kept. Cache records survive: the reloaded package has the same
        // rule-table identity.
        let mut kept_catalogue = None;
        if loaded
            .as_ref()
            .is_some_and(|current| current.application.is_poisoned())
            && let Some(poisoned) = loaded.take()
        {
            eprintln!(
                "vakint: reloading the 5L application: a stopped request left a failed \
                 one-time owner preparation in it"
            );
            kept_catalogue = Some(poisoned.catalogue);
        }
        if loaded.as_ref().is_none_or(|current| {
            current.path != path
                || prepared_identity(current.options) != prepared_identity(options)
                || current.product_limits != settings.rustred_product_limits
        }) {
            let bare = super::data::RustRedData::release_root(&path)?.is_none();
            let started = std::time::Instant::now();
            let identity = super::cache::saved_application_digest(&data, &path)?;
            log::debug!(target: "vakint::saved_application", "phase=native_load_begin");
            let route_addon = route_preparation_addon(&data, official_path.as_deref(), &path)?;
            crate::execution::checkpoint().map_err(VakintError::from)?;
            let prepared = load_saved_application_with_route_addon::<15>(
                &path,
                options,
                match &route_addon {
                    RouteAddonChoice::Declared(addon) => SavedRouteAddon::File(addon),
                    RouteAddonChoice::Rejected => SavedRouteAddon::Disabled,
                    RouteAddonChoice::Undeclared => SavedRouteAddon::Auto,
                },
            )
            .map_err(native_failure)?;
            // RustRed prints its own line when it rejects the add-on.
            let routes = prepared.route_preparation();
            log::debug!(
                target: "vakint::saved_application",
                "phase=native_routes source={}",
                match routes.source {
                    SavedRoutePreparationSource::Addon => "addon",
                    SavedRoutePreparationSource::Verified => "verified",
                }
            );
            if super::cache::saved_application_digest(&data, &path)? != identity {
                return Err(failure(
                    "saved application changed while its native payload was loading",
                )
                .into());
            }
            log::debug!(target: "vakint::saved_application", "phase=native_load_end elapsed_seconds={:.6}", started.elapsed().as_secs_f64());
            crate::execution::checkpoint().map_err(VakintError::from)?;
            let application = Arc::new(SavedApplication::new(
                prepared,
                path.clone(),
                identity,
                settings.rustred_product_limits,
                Arc::clone(&data),
            )?);
            // Only a fully authenticated, prepared publication can activate
            // the optional catalogue. Cache/options reloads reuse its Arc.
            let catalogue = match self
                .topologies
                .five_loop_catalogue()
                .or_else(|| loaded.as_ref().map(|current| &current.catalogue))
                .or(kept_catalogue.as_ref())
            {
                Some(catalogue) => Arc::clone(catalogue),
                None => crate::topologies::FiveLoopCatalogue::prepare()?,
            };
            let replaced = loaded.replace(LoadedApplication {
                path,
                options,
                product_limits: settings.rustred_product_limits,
                application: Arc::clone(&application),
                bare,
                catalogue,
            });
            // The replaced package's cache is unreachable unless both share
            // their rule tables (a copy, or an options reload).
            if replaced.is_some() {
                let mut programs = data.ordered_cache_programs()?;
                programs.push(application.ordered_cache_program().clone());
                self.ordered_cache.detach_except(&programs);
            }
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

    /// Clear the optional frontend witness memo of the loaded package and the
    /// master-transport memo of non-master inputs. Ordered reduction retains
    /// no other result state between requests.
    pub fn clear_rustred_application_witness_cache(&self) -> Result<(), VakintError> {
        let application = self.loaded_application()?;
        application.witnesses.lock().map_err(failure)?.clear();
        application
            .master_transport_misses
            .lock()
            .map_err(failure)?
            .clear();
        Ok(())
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

    /// Evaluate five-loop input while retaining each exact nonzero condition.
    /// This uses the same matcher, tensor prepass and scalar tail as `evaluate`.
    pub fn evaluate_rustred_symbolic(
        &self,
        settings: &VakintSettings,
        input: AtomView,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        self.autosave_ordered_cache(
            self.evaluate_rustred_symbolic_observed(settings, input, None, None),
        )
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
        self.autosave_ordered_cache(self.evaluate_rustred_symbolic_observed(
            settings,
            input,
            Some(observer),
            None,
        ))
    }

    /// Retain exact five-loop masters and guards with optional per-call diagnostics.
    pub fn evaluate_rustred_symbolic_with_monitor(
        &self,
        settings: &VakintSettings,
        input: AtomView,
        monitor: &RustRedMonitorOptions,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        self.autosave_ordered_cache(MonitorContext::run(monitor, |observation, progress| {
            self.evaluate_rustred_symbolic_observed(settings, input, observation, progress)
        }))
    }

    /// Retain exact masters inside an existing caller-owned job monitor.
    /// Returning completes this API call; the outer job scope remains active.
    pub fn evaluate_rustred_symbolic_in_scope(
        &self,
        settings: &VakintSettings,
        input: AtomView,
        scope: &mut RustRedMonitorScope<'_, '_>,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        self.autosave_ordered_cache(scope.evaluate(|observation, progress| {
            self.evaluate_rustred_symbolic_observed(settings, input, observation, progress)
        }))
    }

    fn evaluate_rustred_symbolic_observed(
        &self,
        settings: &VakintSettings,
        input: AtomView,
        mut observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
        progress: Option<&super::monitor::CallerProgress>,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        settings.with_resource_budget(|| {
            if let Some(progress) = progress {
                progress.configure(settings);
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
            // A term equal to a package master is first written as the master's
            // own parent form: canonicalization may move a contracted input onto
            // another (possibly 2-isomorphic) census embedding of its graph.
            for term in &mut expression.0 {
                self.transport_rustred_master_term(settings, &mut term.integral, &mut term.numerator)?;
            }
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
            let mut contributions = Vec::new();
            let split_started = std::time::Instant::now();
            let expression = Atom::from(expression);
            MonitorContext::native_phase(observation.as_deref_mut(), "scalar_input", Some(expression.as_view()));
            let terms = VakintExpression::split_integrals(expression.as_view())?;
            drop(expression);
            log::debug!(target: "vakint::saved_application", "phase=scalar_input_split elapsed_seconds={:.6}", split_started.elapsed().as_secs_f64());
            for mut term in terms {
                crate::execution::checkpoint()?;
                let matching_started = std::time::Instant::now();
                MonitorContext::native_phase(observation.as_deref_mut(), "routing_verification", None);
                // A term equal to a package master is rewritten as the master's
                // own representation and keeps its identity witness (no walk).
                let mut matched = match self.transport_rustred_master_term(
                    settings,
                    &mut term.integral,
                    &mut term.numerator,
                )? {
                    Some(witness) => witness,
                    None => self
                        .match_rustred_input(settings, term.integral.as_view(), term.numerator.as_view())?
                        .ok_or_else(|| failure("input has no registered topology match"))?,
                };
                log::debug!(target: "vakint::saved_application", "phase=frontend_admitted_witness_search elapsed_seconds={:.6}", matching_started.elapsed().as_secs_f64());
                let replacement_started = std::time::Instant::now();
                MonitorContext::native_phase(observation.as_deref_mut(), "numerator_replacement", None);
                matched.apply_replacement_rules()?;
                term.apply_scalar_numerator_replacement_rules(&matched, settings)?;
                log::debug!(target: "vakint::saved_application", "phase=numerator_replacement elapsed_seconds={:.6}", replacement_started.elapsed().as_secs_f64());
                log::debug!(target: "vakint::five_loop_match", "phase=symbolic_candidate_numerator numerator={}", term.numerator);
                contributions.push(MatchedContribution {
                    numerator: term.numerator,
                    matched,
                });
            }
            let (value, conditions, families) = if contributions.is_empty() {
                if let Some(progress) = progress {
                    progress.cohorts(0, 0);
                }
                (Atom::Zero, Vec::new(), BTreeMap::new())
            } else {
                // This parent span includes callee cleanup. Its recorded
                // subphases must not be added to it as disjoint durations.
                let batch_started = std::time::Instant::now();
                let step_cache =
                    self.ordered_step_cache(settings, application.ordered_cache_program())?;
                let result = application.evaluate_ordered_batch(
                    &tensor_settings,
                    &contributions,
                    step_cache,
                    observation.as_deref_mut(),
                    progress,
                );
                log::debug!(target: "vakint::saved_application", "phase=ordered_batch_total elapsed_seconds={:.6}", batch_started.elapsed().as_secs_f64());
                let result = result?;
                (result.value, result.nonzero_conditions, result.families)
            };
            if let Some(progress) = progress {
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
                nonzero_conditions: conditions,
                families,
            };
            log::debug!(target: "vakint::saved_application", "phase=public_result_packing elapsed_seconds={:.6}", packing_started.elapsed().as_secs_f64());
            MonitorContext::native_phase(observation, "retained_output", Some(result.value.as_view()));
            Ok(result)
        })
    }
}

impl SavedApplication {
    /// `identity` is the package's [`super::cache::saved_application_digest`],
    /// the rule-table identity of its Ordered step cache.
    fn new(
        mut prepared: PreparedSavedApplication<15>,
        path: PathBuf,
        identity: [u8; 32],
        product_limits: VacuumProductReductionLimits,
        data: Arc<super::data::RustRedData>,
    ) -> Result<Self, RustRedEvaluationError> {
        let started = std::time::Instant::now();
        // Compile the same complete common basis through RustRed's existing
        // input boundary. No denominator inversion is implemented in Vakint.
        let expected = common_basis_family(prepared.family().name())?;
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
        if prepared
            .output_family_shared(prepared.family().fingerprint())
            .is_none()
        {
            return Err(failure(
                "saved application does not expose its authenticated family",
            ));
        }
        let data_version = data.masters.version.clone();
        let step_cache = super::cache::OrderedCacheProgram::new(
            5,
            prepared.family().name(),
            &data.version,
            identity,
            prepared.ordered_step_cache_binding(),
        );
        Ok(Self {
            step_cache,
            fingerprint: prepared.family().fingerprint().to_owned(),
            data_version,
            prepared: Mutex::new(prepared),
            path,
            products,
            product_charts,
            data,
            witnesses: Mutex::new(WitnessCache::default()),
            master_forms: std::sync::OnceLock::new(),
            master_transport_misses: Mutex::new(BTreeSet::new()),
            poisoned: std::sync::atomic::AtomicBool::new(false),
        })
    }

    /// Whether a stopped request left a failed owner preparation behind.
    pub(super) fn is_poisoned(&self) -> bool {
        self.poisoned.load(std::sync::atomic::Ordering::Acquire)
    }

    /// Deferred owners whose one-time preparation failed (RustRed keeps the
    /// failure; only a reload retries them).
    fn failed_owners(&self) -> Option<usize> {
        self.prepared
            .lock()
            .ok()
            .map(|prepared| prepared.owner_preparation_statistics().failed_owners)
    }

    /// RustRed retries an owner preparation refused by a typed resource
    /// error, but a cancellation or limit reached inside its exact algebra
    /// can arrive as an ordinary failure, which the owner then keeps. After
    /// a request that was cancelled or stopped by a limit and left a new
    /// failed owner, the package is reloaded before its next use.
    fn note_stopped_request(&self, failed_before: Option<usize>, error: &VakintError) {
        let stopped = crate::Cancellation::current().is_some_and(|scope| scope.is_cancelled())
            || stopped_by_resource_limit(error);
        let failed_after = self.failed_owners();
        if stopped
            && let (Some(before), Some(after)) = (failed_before, failed_after)
            && after > before
        {
            self.poisoned
                .store(true, std::sync::atomic::Ordering::Release);
        }
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

    pub(super) fn normalization(
        settings: &VakintSettings,
        loop_count: usize,
    ) -> Result<Atom, VakintError> {
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

    /// Ordered step-cache identity of this package.
    pub(super) fn ordered_cache_program(&self) -> &super::cache::OrderedCacheProgram {
        &self.step_cache
    }

    /// The reduction of a request whose keys are all remaining terminals
    /// (masters) of the package's own family: those terminals with their
    /// coalesced coefficients (each reduces to itself under the published
    /// substitutions). `None` when some key is not such a terminal.
    fn remaining_terminal_request(
        application: &PreparedSavedApplication<15>,
        own_terminals: &BTreeMap<IntegralKey, rustred_app::SavedApplicationIntegral>,
        terms: &[(IntegralKey, Coefficient)],
        settings: &VakintSettings,
    ) -> Result<
        Option<BTreeMap<rustred_app::SavedApplicationIntegral, Coefficient>>,
        RustRedEvaluationError,
    > {
        if terms.is_empty()
            || terms
                .iter()
                .any(|(key, _)| !own_terminals.contains_key(key))
        {
            return Ok(None);
        }
        let context = application.family().coefficient_context();
        let limits = settings
            .rustred_application_options
            .reduction_limits
            .exact_algebra;
        let mut reduced = BTreeMap::new();
        for (key, coefficient) in terms {
            let terminal = own_terminals[key].clone();
            let sum = match reduced.remove(&terminal) {
                Some(previous) => context
                    .try_add(&previous, coefficient, limits)
                    .map_err(native_failure)?,
                None => coefficient.clone(),
            };
            if !sum.is_zero() {
                reduced.insert(terminal, sum);
            }
        }
        Ok(Some(reduced))
    }

    /// Reduce all matched contributions of one request with the ordered
    /// scheduler, one native frontier per (mass, spectator, exponent) channel.
    /// `step_cache` (bound to [`Self::ordered_cache_program`]) is attached to
    /// the native application for this batch, or detached when absent.
    pub(super) fn evaluate_ordered_batch(
        &self,
        settings: &VakintSettings,
        contributions: &[MatchedContribution],
        step_cache: Option<Arc<rustred::reduction::OrderedStepCache>>,
        observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
        progress: Option<&super::monitor::CallerProgress>,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        let failed_before = self.failed_owners();
        let result = self.evaluate_ordered_batch_locked(
            settings,
            contributions,
            step_cache,
            observation,
            progress,
        );
        if let Err(error) = &result {
            self.note_stopped_request(failed_before, error);
        }
        result
    }

    fn evaluate_ordered_batch_locked(
        &self,
        settings: &VakintSettings,
        contributions: &[MatchedContribution],
        step_cache: Option<Arc<rustred::reduction::OrderedStepCache>>,
        mut observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
        progress: Option<&super::monitor::CallerProgress>,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        MonitorContext::native_phase(observation.as_deref_mut(), "scalar_lowering", None);
        let mut application = self.prepared.lock().map_err(failure)?;
        application.set_reduction_strategy(rustred::reduction::ReductionStrategy::Ordered);
        application
            .set_ordered_weight_arithmetic(super::five_loop_ordered_weight_arithmetic(settings));
        let scheduler = super::ordered_scheduler(settings);
        application.set_ordered_scheduler_options(scheduler);
        // The cache is the one of this package's own program (its rule tables
        // and binding); RustRed checks family, variables and policy again
        // here and at every request.
        application
            .set_ordered_step_cache(step_cache)
            .map_err(super::cache::operation)?;
        let family = application.family();
        let loop_count = family.loop_count();
        let family_coefficient_variables = family.coefficient_context().parameter_names().to_vec();
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
        // The remaining terminals (masters) of the package's own family: such a
        // key reduces to itself, inside or outside the recorded entry union.
        let own_terminals = application
            .remaining_terminals()
            .iter()
            .filter(|terminal| terminal.family_fingerprint() == family.fingerprint())
            .map(|terminal| (terminal.integral().clone(), terminal.clone()))
            .collect::<BTreeMap<_, _>>();
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
        // Transport and weighted lowering are independent per contribution
        // and run first (in parallel when configured); everything below then
        // consumes them strictly in contribution order, as before.
        let batch_started = std::time::Instant::now();
        let lowered = super::parallel::map_ordered_by_cost(
            frontend_threads(settings),
            contributions,
            // Lowering cost grows with the numerator; start the largest first.
            |contribution| contribution.numerator.as_view().get_byte_size(),
            |_, contribution| -> Result<_, VakintError> {
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
                // Preserve the received rational domain before identifying
                // epsilon with (4-d)/2 can remove a denominator algebraically.
                let input_conditions =
                    input.input_conditions_after_mapping(&scalar, |condition| {
                        condition
                            .replace(dimension.to_pattern())
                            .with(physical_dimension.to_pattern())
                    })?;
                input.numerator = MatchedScalarInput::numerator_in_family_dimension(
                    input.numerator.as_view(),
                    settings,
                    &dimension,
                );
                let lowering = scalar
                    .lower_weighted(&input.numerator, &input.target, &input.mass_squared)
                    .map_err(native_failure)?;
                log::debug!(target: "vakint::saved_application",
                    "phase=parent_transport_and_scalar_lowering elapsed_seconds={:.6} lowered_count={}",
                    lowering_started.elapsed().as_secs_f64(), lowering.terms().len());
                Ok((input, input_conditions, lowering))
            },
        );
        log::debug!(target: "vakint::saved_application",
            "phase=scalar_lowering contributions={} threads={} elapsed_seconds={:.6}",
            contributions.len(), frontend_threads(settings), batch_started.elapsed().as_secs_f64());
        for (input_index, lowered) in lowered.into_iter().enumerate() {
            let (input, input_conditions, lowering) = lowered?;
            conditions.extend(input_conditions);
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
                if !own_terminals.contains_key(term.integral())
                    && !application
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
                for (spectator, coefficient) in MatchedScalarInput::channel_coefficients(
                    spectator,
                    term.coefficient(),
                    family.coefficient_context(),
                    settings
                        .rustred_application_options
                        .reduction_limits
                        .exact_algebra,
                )? {
                    batches
                        .entry((input.mass_squared.clone(), spectator, exponent))
                        .or_default()
                        .push((term.integral().clone(), coefficient));
                }
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
        // Lower-loop factor products are reduced by a request-local checked
        // factorization primitive; no product result outlives the request.
        let mut products_reducer = {
            let _math = MathSpan::new();
            let anchor = application
                .output_family_shared(application.family().fingerprint())
                .ok_or_else(|| failure("saved application has no output family"))?;
            VacuumProductReducer::try_new(
                ShippedFactorProvider::with_data(Arc::clone(&self.data)),
                anchor,
                settings.rustred_product_limits,
            )
            .map_err(native_failure)?
        };
        let cohort_total = batches.len();
        let mut skipped = 0;
        for (cohort_index, ((mass_squared, spectator, exponent), terms)) in
            batches.into_iter().enumerate()
        {
            if let Some(progress) = progress {
                progress.cohorts(cohort_index, cohort_total);
            }
            // Declare this before the returned native values: the complete
            // exact reduction/assembly interval includes their natural cleanup.
            // With one scheduler thread, deferred owners are prepared on this
            // request's thread (preparation workers only verify routes while
            // the application loads), so their summed preparation time is a
            // part of this span. Concurrent Apply steps may prepare owners at
            // the same time; their summed time then no longer nests here.
            let math = MathSpan::new();
            if scheduler.threads > 1 {
                math.invalidate("ordered scheduler threads > 1");
            }
            let spectator = spectator
                .replace(dimension.to_pattern())
                .with(physical_dimension.to_pattern());
            let mut cleanup = RetainedBatchCleanup::default();
            // A request made of remaining terminals (masters) only is its own
            // reduction: it is built without the native request. A terminal of
            // a composed pack that is not a terminal of the original rule
            // tables would otherwise be walked through them and composed back
            // onto itself, a first touch prepares owners, and the published
            // terminal rows attach the package-wide guard vector, which the
            // identity does not need.
            let direct =
                Self::remaining_terminal_request(&application, &own_terminals, &terms, settings)?;
            let native;
            let direct_terms;
            let (reduced_terms, reduced_conditions): (
                &BTreeMap<rustred_app::SavedApplicationIntegral, Coefficient>,
                &[Coefficient],
            ) = match direct {
                Some(terms) => {
                    MonitorContext::native_phase(
                        observation.as_deref_mut(),
                        "terminal_substitution",
                        None,
                    );
                    log::debug!(target: "vakint::saved_application", "phase=saved_rule_application outcome=remaining_terminals_only terms={}", terms.len());
                    direct_terms = terms;
                    (&direct_terms, &[])
                }
                None => {
                    let options = super::application_options(settings);
                    let mut dump = RequestDump::begin(DumpRequest {
                        loop_count,
                        program: DumpProgram::SavedApplication {
                            path: &self.path,
                            family_fingerprint: application.family().fingerprint(),
                            options: &options,
                        },
                        weight_arithmetic: options.ordered_weight_arithmetic,
                        master_mode: request_dump::master_mode(settings),
                        channel: DumpChannel {
                            index: cohort_index,
                            total: cohort_total,
                            mass_squared: &mass_squared,
                            spectator: &spectator,
                            exponent,
                        },
                        limits: settings.rustred_application_options.reduction_limits,
                        coefficient_variables: family_coefficient_variables.clone(),
                        inputs: &terms,
                        caller_conditions: conditions
                            .iter()
                            .map(|condition| condition.to_canonical_string())
                            .collect(),
                    })?;
                    if dump.as_ref().is_some_and(RequestDump::only) {
                        skipped += 1;
                        continue;
                    }
                    let application_started = std::time::Instant::now();
                    MonitorContext::native_phase(
                        observation.as_deref_mut(),
                        "saved_rule_application",
                        None,
                    );
                    let preparation_before = math
                        .is_active()
                        .then(|| application.owner_preparation_statistics().preparation_time);
                    let reduced = request_dump::observe_native(
                        dump.as_mut(),
                        observation.as_deref_mut(),
                        |observation| {
                            match observation {
                                Some(observation) => application
                                    .reduce_routed_linear_combination_with_observation(
                                        terms,
                                        observation,
                                        None,
                                    ),
                                None => application.reduce_routed_linear_combination(terms),
                            }
                            .map_err(super::cache::operation)
                        },
                    )??;
                    if let Some(dump) = dump {
                        dump.finish(
                            reduced.terms().iter().map(|(terminal, coefficient)| {
                                (
                                    terminal.family_fingerprint(),
                                    terminal.integral(),
                                    coefficient,
                                )
                            }),
                            reduced.nonzero_conditions(),
                            observation
                                .as_deref()
                                .and_then(|observer| observer.last_ordered_progress()),
                        )?;
                    }
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
                    native = reduced;
                    (native.terms(), native.nonzero_conditions())
                }
            };
            MonitorContext::native_phase(observation.as_deref_mut(), "terminal_geometry", None);
            let products_started = std::time::Instant::now();
            let mut contribution = Atom::Zero;
            let mut products = Vec::new();
            let mut product_conditions = reduced_conditions.to_vec();
            for (terminal, coefficient) in reduced_terms {
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
            let normalized = products_reducer.reduce_products(products, product_conditions);
            let normalized = match normalized {
                Ok(normalized) => normalized,
                Err(rustred::reduction::terminal_normalization::VacuumProductReductionError::UnresolvedFactor { family, target }) => {
                    log::debug!(target: "vakint::saved_application", "phase=terminal_geometry_and_lower_factors elapsed_seconds={:.6}", products_started.elapsed().as_secs_f64());
                    let fallback_started = std::time::Instant::now();
                    MonitorContext::native_phase(observation.as_deref_mut(), "unprojected_terminal_materialization", None);
                    log::warn!("RustRed lower-factor coverage is incomplete for {family} {:?}; retaining the exact parent terminals", target.powers());
                    // Preserve the original reduced sum and its complete guards;
                    // no failed or partially normalized product enters a cache.
                    let original = reduced_terms.iter().map(|(terminal, coefficient)| {
                        coefficient.to_expression().replace(dimension.to_pattern()).with(physical_dimension.to_pattern())
                            * Self::symbolic_integral(terminal.family_fingerprint(), terminal.integral(), vk_symbol!("RustRedMaster"), &mut families)
                    }).sum::<Atom>();
                    for condition in reduced_conditions {
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
            // Each guard is materialized independently; chunks are consumed in
            // input order, so the first invalid guard reports as sequentially.
            let threads = super::postprocessing_threads(settings);
            let dimension_symbols = dimension.get_all_symbols(true);
            let dimension_pattern = dimension.to_pattern();
            let physical_pattern = physical_dimension.to_pattern();
            let guard_chunk = GUARD_CHUNK.max(threads);
            for (chunk_index, chunk) in normalized
                .nonzero_conditions()
                .chunks(guard_chunk)
                .enumerate()
            {
                let materialized = parallel::map_ordered(threads, chunk, |_, condition| {
                    let condition = condition.to_expression();
                    (!condition.is_zero()
                        && condition
                            .get_all_symbols(true)
                            .iter()
                            .all(|symbol| dimension_symbols.contains(symbol)))
                    .then(|| {
                        condition
                            .replace(&dimension_pattern)
                            .with(&physical_pattern)
                    })
                });
                for condition in materialized {
                    conditions.insert(condition.ok_or_else(|| failure(
                        "saved application has a non-generic-d guard; the symbolic Vakint result cannot discharge it",
                    ))?);
                }
                if let Some(progress) = &mut guard_progress {
                    progress.completed(chunk_index * guard_chunk + chunk.len());
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
        if let Some(progress) = progress {
            progress.cohorts(cohort_total, cohort_total);
        }
        if skipped != 0 {
            return Err(request_dump::skipped(skipped).into());
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
        // The request-local product state is part of the timed reduction work.
        {
            let _math = MathSpan::new();
            drop(products_reducer);
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
    #[ignore = "loads the official five-loop package twice over (minutes); needs VAKINT_TEST_EXTERNAL_DATA"]
    fn weight_arithmetic_changes_reuse_the_loaded_application() {
        let artifacts = PathBuf::from(
            std::env::var_os("VAKINT_TEST_EXTERNAL_DATA")
                .expect("set VAKINT_TEST_EXTERNAL_DATA to the official external-data checkout"),
        );
        let mut vakint = Vakint::new().unwrap();
        vakint
            .load_official_external_data(
                5,
                &crate::OfficialExternalDataOptions {
                    local_artifacts_path: Some(artifacts),
                    ..Default::default()
                },
            )
            .unwrap();
        let loaded = vakint.loaded_application().unwrap();
        let mut settings = VakintSettings::default();
        for arithmetic in [
            "modular:threads=1",
            "modular:threads=2,batch=8",
            "modular_polynomial",
            "exact",
        ] {
            settings.ordered_weight_arithmetic = arithmetic.parse().unwrap();
            let reused = vakint.saved_application(&settings).unwrap();
            assert!(
                Arc::ptr_eq(&loaded, &reused),
                "{arithmetic} reloaded the package"
            );
        }
        // n_cores sets the scheduler threads and the route-verification
        // workers, neither of which is part of the prepared application.
        for n_cores in [1, 2] {
            settings.n_cores = n_cores;
            let reused = vakint.saved_application(&settings).unwrap();
            assert!(
                Arc::ptr_eq(&loaded, &reused),
                "n_cores = {n_cores} reloaded the package"
            );
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
    fn scalar_numerator_defaults_and_unloaded_application_statistics() {
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
        assert!(matches!(
            Vakint::new()
                .unwrap()
                .rustred_application_owner_preparation_statistics(),
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
        // application. Earlier clones retain their previous snapshot.
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
        let prepared = application.prepared.lock().unwrap();
        // The same request-local checked factorization as ordered reduction.
        let mut reducer = VacuumProductReducer::try_new(
            ShippedFactorProvider::with_data(Arc::clone(&application.data)),
            prepared
                .output_family_shared(prepared.family().fingerprint())
                .unwrap(),
            settings.rustred_product_limits,
        )
        .unwrap();
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
