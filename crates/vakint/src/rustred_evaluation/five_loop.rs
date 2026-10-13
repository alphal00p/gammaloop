//! Exact symbolic application of a separately loaded immutable publication.
//!
//! Masters denote unit-mass Minkowski integrals with unnormalized loop measures,
//! qualified by the authenticated output-family fingerprint. Coefficients are
//! exact in generic d=4-2*epsilon. Saved-program guards depend only on d;
//! caller numerator and common-mass pole guards retain their scalar parameters.
//! This adapter performs no scalar reduction, master substitution or IBP search.

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
    path::{Path, PathBuf},
    sync::{Arc, LazyLock, Mutex},
};

use rustred::algebra::Coefficient;
use rustred::family::IntegralKey;
use rustred::input::{Compiler, Limits, LoweringLimits, TextProject, TextPropagator};
use rustred::reduction::terminal_normalization::{
    VacuumFactor, VacuumProductCacheStatistics, VacuumProductChartLimits, VacuumProductReducer,
    VacuumProductReductionLimits,
};
use rustred::scalar_numerator::{FamilyScalarNumeratorService, ScalarNumeratorLimits};
use rustred_app::{
    PreparedSavedApplication, SavedApplicationOptions, SavedApplicationProductChartStatistics,
    SavedApplicationProductStatistics, load_saved_application,
};
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};

use super::{
    RustRedEvaluationError, factor_projection,
    factor_provider::ShippedFactorProvider,
    numerator::{MatchedScalarInput, ScalarNumeratorForm},
};
use crate::topologies::five_loop::{ParentDescriptor, common_momentum_rows, descriptors};
use crate::utils::vakint_macros::vk_symbol;
use crate::{
    ReplacementRules, Topology, Vakint, VakintError, VakintExpression, VakintSettings, symbols::S,
};

static PARENTS: LazyLock<Vec<ParentDescriptor>> = LazyLock::new(descriptors);

#[cfg(test)]
#[path = "five_loop/admission_census_tests.rs"]
mod admission_census_tests;

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

pub(crate) struct SavedApplication {
    prepared: Mutex<PreparedSavedApplication<15>>,
    fingerprint: String,
    products: SavedApplicationProductStatistics,
    product_charts: SavedApplicationProductChartStatistics,
    product_reducer: Mutex<VacuumProductReducer<ShippedFactorProvider>>,
}

#[derive(Debug)]
pub(crate) struct LoadedApplication {
    path: PathBuf,
    options: SavedApplicationOptions,
    product_limits: VacuumProductReductionLimits,
    application: Arc<SavedApplication>,
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
        self.topologies
            .match_topologies_to_user_input_with_acceptance(
                integral,
                settings.allow_unknown_integrals,
                settings.rustred_momentum_match_limits,
                |matched| self.rustred_matched_input_key(settings, numerator, matched),
            )
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
        .map_err(failure)?;
        let lowered = scalar
            .lower_weighted(&input.numerator, &input.target, &input.mass_squared)
            .map_err(failure)?;
        let mut key = Vec::with_capacity(lowered.terms().len());
        for term in lowered.terms() {
            if !prepared
                .routed_entry_is_admitted(term.integral())
                .map_err(failure)?
            {
                log::debug!(target: "vakint::five_loop_match", "phase=candidate_rejected lowered_count={} rejected_key={:?} lowered_keys_prefix={:?}", lowered.terms().len(), term.integral().powers(), lowered.terms().iter().take(4).map(|term| term.integral().powers()).collect::<Vec<_>>());
                return Ok(None);
            }
            key.push((
                term.integral().clone(),
                (term.coefficient().to_expression() * term.scalar_spectator()).together(),
            ));
        }
        key.sort();
        Ok(Some(key))
    }

    /// Load one authenticated application package and reuse its owned memoization.
    ///
    /// The five-loop backend is symbolic only (`substitute_masters=false`).
    /// Clones share this immutable package and synchronized application cache.
    /// Only successful loading prepares and publishes the five-loop catalogue;
    /// `Vakint::new` prepares just the ordinary lower-loop registry.
    pub fn with_rustred_application(
        self,
        settings: &mut VakintSettings,
        path: impl AsRef<Path>,
    ) -> Result<Self, VakintError> {
        let mut configured = settings.clone();
        configured.rustred_artifact_path = Some(path.as_ref().canonicalize().map_err(failure)?);
        self.saved_application(&configured)?;
        *settings = configured;
        Ok(self)
    }

    pub(crate) fn saved_application(
        &self,
        settings: &VakintSettings,
    ) -> Result<Arc<SavedApplication>, VakintError> {
        let path = settings
            .rustred_artifact_path
            .as_ref()
            .ok_or(RustRedEvaluationError::SavedApplicationUnavailable)?;
        let path = path.canonicalize().map_err(failure)?;
        let options = settings.rustred_application_options;
        let mut loaded = self.rustred_application.lock().map_err(failure)?;
        if loaded.as_ref().is_none_or(|current| {
            current.path != path
                || current.options != options
                || current.product_limits != settings.rustred_product_limits
        }) {
            let started = std::time::Instant::now();
            log::debug!(target: "vakint::saved_application", "phase=native_load_begin");
            let prepared = load_saved_application::<15>(&path, options).map_err(failure)?;
            log::debug!(target: "vakint::saved_application", "phase=native_load_end elapsed_seconds={:.6}", started.elapsed().as_secs_f64());
            let application = Arc::new(SavedApplication::new(
                prepared,
                settings.rustred_product_limits,
            )?);
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

    fn loaded_application(&self) -> Result<Arc<SavedApplication>, VakintError> {
        let loaded = self.rustred_application.lock().map_err(failure)?;
        loaded
            .as_ref()
            .map(|current| Arc::clone(&current.application))
            .ok_or_else(|| RustRedEvaluationError::SavedApplicationUnavailable.into())
    }

    /// Clear the explicitly loaded package's cache without changing its entry scope.
    pub fn clear_rustred_application_cache(&self) -> Result<(), VakintError> {
        let application = self.loaded_application()?;
        application
            .prepared
            .lock()
            .map_err(failure)?
            .clear_cache()
            .map_err(failure)?;
        application
            .product_reducer
            .lock()
            .map_err(failure)?
            .clear_cache()
            .map_err(failure)?;
        Ok(())
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
        self.validate_settings(settings)?;
        let application = self.saved_application(settings)?;
        let mut expression = VakintExpression::try_from(input)?;
        log::debug!(target: "vakint::five_loop_match", "phase=symbolic_input input={input}");
        let canonical_started = std::time::Instant::now();
        expression.canonicalize(settings, &self.topologies, false)?;
        log::debug!(target: "vakint::saved_application", "phase=frontend_canonicalization elapsed_seconds={:.6}", canonical_started.elapsed().as_secs_f64());
        log::debug!(target: "vakint::five_loop_match", "phase=symbolic_canonical input={expression}");
        // This is an internal scalar handoff, not the public tensor-output
        // presentation boundary. Reopening scalar powers here can manufacture
        // a rank22 indexed expression before a rank0 scalar application.
        let mut tensor_settings = settings.clone();
        tensor_settings.use_dot_product_notation = true;
        let tensor_started = std::time::Instant::now();
        expression.tensor_reduce(self, &tensor_settings)?;
        log::debug!(target: "vakint::saved_application", "phase=frontend_tensor_prepass elapsed_seconds={:.6}", tensor_started.elapsed().as_secs_f64());
        log::debug!(target: "vakint::five_loop_match", "phase=symbolic_tensor_reduced input={expression}");
        let mut value = Atom::Zero;
        let mut conditions = BTreeSet::new();
        let mut families = BTreeMap::new();
        for mut term in VakintExpression::split_integrals(Atom::from(expression).as_view())? {
            let matching_started = std::time::Instant::now();
            let mut matched = self
                .match_rustred_input(settings, term.integral.as_view(), term.numerator.as_view())?
                .ok_or_else(|| failure("input has no registered topology match"))?;
            log::debug!(target: "vakint::saved_application", "phase=frontend_admitted_witness_search elapsed_seconds={:.6}", matching_started.elapsed().as_secs_f64());
            matched.apply_replacement_rules()?;
            term.apply_scalar_numerator_replacement_rules(&matched, settings)?;
            log::debug!(target: "vakint::five_loop_match", "phase=symbolic_candidate_numerator numerator={}", term.numerator);
            let result =
                application.evaluate_detailed(settings, term.numerator.as_view(), &matched)?;
            value += result.value;
            conditions.extend(result.nonzero_conditions);
            families.extend(result.families);
        }
        Ok(RustRedSymbolicEvaluation {
            value,
            nonzero_conditions: conditions.into_iter().collect(),
            families,
        })
    }
}

impl SavedApplication {
    fn new(
        mut prepared: PreparedSavedApplication<15>,
        product_limits: VacuumProductReductionLimits,
    ) -> Result<Self, RustRedEvaluationError> {
        let started = std::time::Instant::now();
        // Compile the same complete common basis through RustRed's existing
        // input boundary. No denominator inversion is implemented in Vakint.
        let expected = Compiler::new(Limits::default())
            .map_err(failure)?
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
            .map_err(failure)?
            .into_lowered(LoweringLimits::default())
            .map_err(failure)?
            .into_family();
        if prepared.family().fingerprint() != expected.fingerprint() {
            return Err(failure(
                "saved package does not use the registered common momentum basis",
            ));
        }
        log::debug!(target: "vakint::saved_application", "phase=family_verified elapsed_seconds={:.6}", started.elapsed().as_secs_f64());
        let products = prepared
            .prepare_terminal_products(product_limits.products)
            .map_err(failure)?;
        let product_charts = prepared
            .prepare_terminal_product_charts(VacuumProductChartLimits {
                products: product_limits.products,
                ..Default::default()
            })
            .map_err(failure)?;
        log::debug!(target: "vakint::saved_application", "phase=products_prepared elapsed_seconds={:.6}", started.elapsed().as_secs_f64());
        let anchor = prepared
            .output_family_shared(prepared.family().fingerprint())
            .ok_or_else(|| failure("saved application does not expose its authenticated family"))?;
        let product_reducer =
            VacuumProductReducer::try_new(ShippedFactorProvider::default(), anchor, product_limits)
                .map_err(failure)?;
        Ok(Self {
            fingerprint: prepared.family().fingerprint().to_owned(),
            prepared: Mutex::new(prepared),
            products,
            product_charts,
            product_reducer: Mutex::new(product_reducer),
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

    pub(super) fn evaluate(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        matched: &ReplacementRules,
    ) -> Result<Atom, VakintError> {
        Ok(self.evaluate_detailed(settings, numerator, matched)?.value)
    }

    fn evaluate_detailed(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        matched: &ReplacementRules,
    ) -> Result<RustRedSymbolicEvaluation, VakintError> {
        let lowering_started = std::time::Instant::now();
        let descriptor = select(&matched.canonical_topology)?;
        let input = MatchedScalarInput::prepare(
            numerator,
            matched,
            &descriptor.momenta,
            descriptor.physical_slots,
            15,
            settings,
            ScalarNumeratorForm::PreserveProducts,
        )?;
        log::debug!(target: "vakint::five_loop_match", "phase=symbolic_parent_numerator numerator={}", input.numerator);
        let mut application = self.prepared.lock().map_err(failure)?;
        let family = application.family();
        let loop_count = family.loop_count();
        let dimension = family.dimension().to_expression();
        let scalar = FamilyScalarNumeratorService::try_new(
            family,
            S.dot,
            (1..=loop_count).map(|axis| S.k.call(axis)).collect(),
            settings.rustred_scalar_numerator_limits,
        )
        .map_err(failure)?;
        let lowering = scalar
            .lower_weighted(&input.numerator, &input.target, &input.mass_squared)
            .map_err(failure)?;
        log::debug!(target: "vakint::saved_application", "phase=parent_transport_and_scalar_lowering elapsed_seconds={:.6}", lowering_started.elapsed().as_secs_f64());
        log::debug!(target: "vakint::saved_application",
            "phase=scalar_lowered target={:?} lowered_count={} lowered_keys_prefix={:?}",
            input.target.powers(), lowering.terms().len(),
            lowering.terms().iter().take(8).map(|term| term.integral().powers()).collect::<Vec<_>>());
        let mut batches = BTreeMap::<(Atom, i64), Vec<(IntegralKey, Coefficient)>>::new();
        for term in lowering.terms() {
            if !application
                .routed_entry_is_admitted(term.integral())
                .map_err(failure)?
            {
                return Err(RustRedEvaluationError::SavedApplicationEntryNotAdmitted {
                    powers: term.integral().powers().to_vec(),
                }
                .into());
            }
            // Unit-mass output masters include their full integration measure.
            // RustRed has already combined numerator mass factors and scalar
            // weights before admission. Each remaining key scales as
            // (m²)^(Ld/2-sum(a)); no expanded zero endpoint is admitted alone.
            let exponent = term
                .integral()
                .powers()
                .iter()
                .try_fold(0i128, |sum, power| sum.checked_sub(i128::from(*power)))
                .ok_or_else(|| failure("saved application mass exponent overflow"))?;
            let exponent = i64::try_from(exponent)
                .map_err(|_| RustRedEvaluationError::MassExponentOverflow { exponent })?;
            batches
                .entry((term.scalar_spectator().clone(), exponent))
                .or_default()
                .push((term.integral().clone(), term.coefficient().clone()));
        }
        let epsilon = Atom::var(vk_symbol!(settings.epsilon_symbol.as_str()));
        let physical_dimension = Atom::num(4) - Atom::num(2) * &epsilon;
        let mut value = Atom::Zero;
        // Caller-owned scalar and mass poles remain conditions even when their
        // weighted integral cancels. They are separate from saved-program
        // guards, whose authenticated coefficient context is d-only.
        let mut conditions = lowering
            .nonzero_conditions()
            .iter()
            .map(|condition| {
                condition
                    .replace(dimension.to_pattern())
                    .with(physical_dimension.to_pattern())
            })
            .collect::<BTreeSet<_>>();
        let mut families = BTreeMap::new();
        for ((spectator, exponent), terms) in batches {
            let application_started = std::time::Instant::now();
            let reduced = application
                .reduce_routed_linear_combination(terms)
                .map_err(failure)?;
            log::debug!(target: "vakint::saved_application", "phase=saved_rule_application elapsed_seconds={:.6}", application_started.elapsed().as_secs_f64());
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
                            .map_err(failure)?;
                        products.push((factors, coefficient));
                    }
                    None => {
                        if let Some(chart) = application.terminal_product_chart(terminal) {
                            let projection = chart
                                .project(
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
                                )
                                .map_err(failure)?;
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
                                    .map_err(failure)?;
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
            let normalized = self
                .product_reducer
                .lock()
                .map_err(failure)?
                .reduce_products(products, product_conditions)
                .map_err(failure)?;
            log::debug!(target: "vakint::saved_application", "phase=terminal_geometry_and_lower_factors elapsed_seconds={:.6}", products_started.elapsed().as_secs_f64());
            for condition in normalized.nonzero_conditions() {
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
            }
            for (product, coefficient) in normalized.terms() {
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
            }
            let mass_power =
                Atom::num((loop_count as i64, 2)) * &physical_dimension + Atom::num(exponent);
            value += contribution * spectator * input.mass_squared.clone().pow(mass_power);
        }
        let normalization = settings
            .get_integral_normalization_factor_atom()?
            .replace(S.n_loops.to_pattern())
            .with(Atom::num(loop_count).to_pattern())
            .replace(Atom::var(vk_symbol!("log_mu_sq")).to_pattern())
            .with(
                Symbol::LOG
                    .call(Atom::var(vk_symbol!(settings.mu_r_sq_symbol.as_str())))
                    .to_pattern(),
            );
        if !matches!(input.mass_squared.as_view(), AtomView::Num(_)) {
            conditions.insert(input.mass_squared);
        }
        let value = value * normalization;
        let value = if settings.use_dot_product_notation {
            value
        } else {
            // All integrated scalar products have been eliminated. Apply the
            // requested presentation only to surviving external spectators.
            Vakint::convert_from_dot_notation(value.as_view())
        };
        Ok(RustRedSymbolicEvaluation {
            value,
            nonzero_conditions: conditions.into_iter().collect(),
            families,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vakint_macros::vk_parse;
    use crate::{EvaluationMethod, EvaluationOrder, RustRedEvaluationOptions};

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
                    substitute_masters: false,
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
        assert!(matches!(
            vakint.evaluate_integral(&settings, root().as_view()),
            Err(VakintError::UnreckognizedIntegral(..))
        ));
        assert!(matches!(
            vakint.evaluate_rustred_symbolic(&settings, root().as_view()),
            Err(VakintError::RustRedEvaluation(
                RustRedEvaluationError::SavedApplicationUnavailable
            ))
        ));
    }

    #[test]
    fn wrong_artifact_path_and_numerical_substitution_fail_without_form() {
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
        // Independently exercise numerical-master rejection after explicit
        // test-only matching activation; no failed load publishes a catalogue.
        vakint
            .topologies
            .prepare_five_loop_catalogue_for_tests()
            .unwrap();
        settings.evaluation_order = EvaluationOrder::rustred_only();
        let error = vakint
            .evaluate_integral(&settings, root().as_view())
            .unwrap_err();
        assert!(matches!(
            error,
            VakintError::RustRedEvaluation(
                RustRedEvaluationError::FiveLoopNumericalMastersUnsupported
            )
        ));
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
                .evaluate_detailed(&strict, term.numerator.as_view(), &matched)
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
    #[ignore = "requires VAKINT_SAVED_APPLICATION pointing to an external generation-54 package"]
    fn external_package_all_factorized_terminals_reduce_through_existing_lower_programs() {
        use rustred::reduction::terminal_normalization::VacuumProductProjection;
        let path = std::env::var_os("VAKINT_SAVED_APPLICATION")
            .expect("set VAKINT_SAVED_APPLICATION to the uncommitted external package");
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
            .with_rustred_application(&mut settings, PathBuf::from(path))
            .unwrap();
        let application = vakint.loaded_application().unwrap();
        let cloned_application = before_load_clone.saved_application(&settings).unwrap();
        assert!(Arc::ptr_eq(&application, &cloned_application));
        assert!(Arc::ptr_eq(
            vakint.topologies.five_loop_catalogue().unwrap(),
            before_load_clone.topologies.five_loop_catalogue().unwrap(),
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
                chart
                    .project(
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
        assert_eq!(positive_count, 38, "pinned generation-54 scalar products");
        assert_eq!(
            positive_supports.len(),
            19,
            "all factorized support classes"
        );
        assert_eq!(
            projected_count, 14,
            "pinned generation-54 numerator product inventory"
        );
        assert_eq!(cross_count, 3, "native cross-ISP transport inventory");
    }
}
