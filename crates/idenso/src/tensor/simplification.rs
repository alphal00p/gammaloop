//! One domain planner for structural contraction and selected tensor identities.
//! Definitions stay shared; only a changed domain or dependent call site is
//! reconsidered. Materialization and explicit notation undo are separate APIs.

use super::{
    ContractSettings, SymbolicTensor,
    aliases::{AliasInterfaces, Definition, TENSOR_ALIAS_SYMBOL},
    inference::TensorInferenceError,
};
use crate::{
    color::ColorSimplifySettings, dirac::GammaSimplifySettings, epsilon::EpsilonSimplifierPass,
};
use spenso::structure::partial::PartialStructure;
use std::{
    collections::{HashMap, HashSet, VecDeque},
    sync::Arc,
};
use symbolica::atom::{AliasedAtom, Atom, AtomView};

pub(crate) mod observation;
use observation::{DomainObservations, SettledRegions};

/// Completion is relative to the requested capabilities and output forms.
/// A deferred or capped result retains the exact expression and all definitions.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub enum ReductionStatus {
    #[default]
    Complete,
    Deferred,
    Capped,
}

/// Independent algebra capabilities. Enabling a family authorizes its required
/// structural contractions, but never enables a different identity family.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AlgebraSettings {
    pub gamma: Option<GammaSimplifySettings>,
    pub color: Option<ColorSimplifySettings>,
    pub epsilon: bool,
    /// Maximum successful transformations of each domain in this invocation.
    pub max_passes: usize,
}

impl Default for AlgebraSettings {
    fn default() -> Self {
        Self {
            gamma: None,
            color: None,
            epsilon: false,
            max_passes: 16,
        }
    }
}

impl AlgebraSettings {
    pub fn hep() -> Self {
        Self {
            gamma: Some(GammaSimplifySettings::default()),
            color: Some(ColorSimplifySettings::default()),
            epsilon: true,
            ..Self::default()
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ReductionPolicy {
    Algebra(AlgebraSettings),
    Contract {
        representations: Option<Vec<spenso::structure::representation::LibraryRep>>,
        metrics: bool,
        rank_one: bool,
        collect_chains: bool,
        collect_traces: bool,
        order: Option<Vec<usize>>,
        max_passes: usize,
    },
}

#[derive(Clone, Debug)]
pub(crate) struct ReductionRecord {
    // Every entry describes this exact immutable root and alias DAG. A no-op
    // under another policy preserves those facts; any domain change drops them.
    completed: Vec<ReductionPolicy>,
    pub(crate) status: ReductionStatus,
}

#[derive(Clone, Copy)]
pub(crate) enum ReductionRequest<'a> {
    Algebra(&'a AlgebraSettings),
    Contract(ContractSettings<'a>),
}

const STRUCTURE: u8 = 1;
const GAMMA: u8 = 2;
const COLOR: u8 = 4;
const EPSILON: u8 = 8;

impl ReductionRequest<'_> {
    fn policy(self) -> ReductionPolicy {
        match self {
            Self::Algebra(settings) => ReductionPolicy::Algebra(*settings),
            Self::Contract(settings) => ReductionPolicy::Contract {
                representations: settings.representations.map(<[_]>::to_vec),
                metrics: settings.metrics,
                rank_one: settings.rank_one,
                collect_chains: settings.collect_chains,
                collect_traces: settings.collect_traces,
                order: settings.order.map(<[_]>::to_vec),
                max_passes: settings.max_passes,
            },
        }
    }

    fn budget(self) -> usize {
        match self {
            Self::Algebra(settings) => settings.max_passes,
            Self::Contract(settings) => settings.max_passes,
        }
    }

    fn enabled(self) -> u8 {
        match self {
            Self::Contract(_) => STRUCTURE,
            Self::Algebra(settings) => {
                (u8::from(settings.gamma.is_some()) * GAMMA)
                    | (u8::from(settings.color.is_some()) * COLOR)
                    | (u8::from(settings.epsilon) * EPSILON)
            }
        }
    }
}

/// Scratch for one domain; the payload remains the existing symbolic tensor.
struct DomainPlan {
    value: SymbolicTensor<PartialStructure>,
    observed: Arc<DomainObservations>,
    dependencies: HashSet<Atom>,
    pending: u8,
    effective: u8,
    prerequisites: bool,
    changes: usize,
    deferred: u8,
    settled: HashMap<u8, SettledRegions>,
    status: ReductionStatus,
}

impl DomainPlan {
    fn new(
        value: SymbolicTensor<PartialStructure>,
        registry: &HashMap<Atom, Definition>,
        pending: u8,
    ) -> Self {
        let observed = Arc::clone(value.reduction_observations());
        let dependencies = observed
            .literals
            .iter()
            .filter(|literal| registry.contains_key(*literal))
            .cloned()
            .collect();
        let prerequisites = value.proofs.prerequisites && pending & GAMMA != 0;
        let mut plan = Self {
            value,
            observed,
            dependencies,
            pending: pending | if prerequisites { STRUCTURE } else { 0 },
            effective: 0,
            prerequisites,
            changes: 0,
            deferred: 0,
            settled: HashMap::new(),
            status: ReductionStatus::Complete,
        };
        plan.effective = plan.local_candidates();
        plan
    }

    fn local_candidates(&self) -> u8 {
        self.observed.identity_candidates() | STRUCTURE
    }
}

impl SymbolicTensor<PartialStructure> {
    /// Reduce only the selected identity families, keeping scalar sums factored.
    pub fn simplify_algebra(
        &self,
        settings: &AlgebraSettings,
    ) -> Result<Arc<SymbolicTensor<AliasInterfaces, AliasedAtom>>, TensorInferenceError> {
        self.ensure_validated()?;
        self.reduction_observations();
        Arc::new(self.clone().with_aliases([])?).plan_reduction(ReductionRequest::Algebra(settings))
    }
}

impl SymbolicTensor<AliasInterfaces, AliasedAtom> {
    pub fn simplify_algebra(
        self: &Arc<Self>,
        settings: &AlgebraSettings,
    ) -> Result<Arc<Self>, TensorInferenceError> {
        self.plan_reduction(ReductionRequest::Algebra(settings))
    }

    pub fn reduction_status(&self) -> ReductionStatus {
        self.proofs
            .reduction
            .as_ref()
            .map_or(ReductionStatus::Deferred, |record| record.status)
    }

    /// Both public settings types lower into this queue. A head is an admission
    /// candidate; an unchanged kernel discharges it until its inputs change.
    pub(crate) fn plan_reduction(
        self: &Arc<Self>,
        request: ReductionRequest<'_>,
    ) -> Result<Arc<Self>, TensorInferenceError> {
        #[cfg(feature = "reference-cases")]
        let _phase =
            crate::reference_cases::timing::scope(crate::reference_cases::timing::Phase::Planning);
        #[cfg(feature = "reference-cases")]
        let diagnostic_operation = match request {
            ReductionRequest::Contract(_) => "contract",
            ReductionRequest::Algebra(_) => "simplify_algebra",
        };
        let policy = request.policy();
        if let Some(record) = &self.proofs.reduction
            && record.completed.contains(&policy)
        {
            #[cfg(feature = "reference-cases")]
            crate::reference_cases::timing::outcome(
                diagnostic_operation,
                ReductionStatus::Complete,
            );
            let mut result = Arc::clone(self);
            if record.status != ReductionStatus::Complete {
                Arc::make_mut(&mut result)
                    .proofs
                    .reduction
                    .as_mut()
                    .unwrap()
                    .status = ReductionStatus::Complete;
            }
            return Ok(result);
        }
        let enabled = request.enabled();
        let structural = match request {
            ReductionRequest::Contract(settings) => settings,
            ReductionRequest::Algebra(_) => ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            },
        };
        let original = self.aliases()?;
        let mut registry: HashMap<_, _> = original
            .iter()
            .cloned()
            .map(|mut pair| {
                if pair.1.proofs.literal_origin.is_none()
                    && matches!(pair.0.expression.as_view(), AtomView::Fun(function) if function.get_symbol() == *TENSOR_ALIAS_SYMBOL)
                {
                    pair.1.proofs.literal_origin = Some(Arc::new(pair.1.clone()));
                }
                (pair.0.expression.clone(), pair)
            })
            .collect();
        let mut domains = HashMap::<Option<Atom>, DomainPlan>::new();
        domains.insert(None, DomainPlan::new(self.root(), &registry, enabled));
        let reachable = self.reachable_definitions();
        for (handle, _) in original {
            if reachable.contains(&handle.expression) {
                let body = registry[&handle.expression].1.clone();
                let mut plan = DomainPlan::new(body, &registry, enabled);
                for dependency in &plan.dependencies {
                    if let Some(previous) = domains.get(&Some(dependency.clone())) {
                        plan.effective |= previous.effective;
                    }
                }
                domains.insert(Some(handle.expression), plan);
            }
        }
        let mut dependents = HashMap::<Atom, HashSet<Option<Atom>>>::new();
        for (key, plan) in &domains {
            for dependency in &plan.dependencies {
                dependents
                    .entry(dependency.clone())
                    .or_default()
                    .insert(key.clone());
            }
        }
        let mut queue = VecDeque::from([None]);
        let mut initial = domains.keys().filter_map(Clone::clone).collect::<Vec<_>>();
        initial.sort_unstable_by(|a, b| AtomView::cmp(&a.as_view(), &b.as_view()));
        queue.extend(initial.into_iter().map(Some));
        let mut queued = queue.iter().cloned().collect::<HashSet<_>>();
        let mut alias_candidates = domains
            .iter()
            .filter_map(|(key, plan)| key.as_ref().map(|key| (key.clone(), plan.effective)))
            .collect::<HashMap<_, _>>();
        let mut available = None;
        let mut changed = false;
        let mut status = ReductionStatus::Complete;
        while let Some(key) = queue.pop_front() {
            queued.remove(&key);
            // Local occurrence information crosses an alias boundary through its
            // dependency edge, never through a registry-wide family flag.
            let candidates = domains[&key].dependencies.iter().fold(
                domains[&key].local_candidates(),
                |mask, dependency| {
                    mask | domains
                        .get(&Some(dependency.clone()))
                        .map_or(0, |plan| plan.effective)
                },
            );
            if let Some(handle) = &key {
                alias_candidates.insert(handle.clone(), candidates);
                // A dependency can expose a candidate without rewriting this
                // body's syntax. Propagate that changed observation through
                // call sites as well as propagating actual body replacements.
                if domains[&key].effective != candidates
                    && let Some(users) = dependents.get(handle)
                {
                    for user in users {
                        let dependent = domains.get_mut(user).unwrap();
                        dependent.settled.clear();
                        dependent.pending = enabled
                            | if dependent.prerequisites {
                                STRUCTURE
                            } else {
                                0
                            };
                        if queued.insert(user.clone()) {
                            queue.push_back(user.clone());
                        }
                    }
                }
            }
            let plan = domains.get_mut(&key).unwrap();
            plan.effective = candidates;
            let mut eligible = plan.pending & candidates;
            plan.pending &= candidates;
            // An identity authorizes cleanup of its emitted connections, but
            // many emitted definitions contain only free ports. Discharge those
            // prerequisites from their existing observations before allocating
            // contraction state or charging the domain's work budget.
            if eligible & STRUCTURE != 0
                && (structural.order.is_none() || key.is_some() || plan.changes != 0)
                && !plan.value.structural_work_permitted(structural)
            {
                plan.pending &= !STRUCTURE;
                eligible &= !STRUCTURE;
                plan.deferred &= !STRUCTURE;
                if plan.deferred == 0 {
                    plan.status = ReductionStatus::Complete;
                }
            }
            if eligible == 0 {
                continue;
            }
            if plan.changes >= request.budget() {
                plan.status = ReductionStatus::Capped;
                status = status.max(plan.status);
                continue;
            }
            // Prefer cheap contractions and short words to a determinant or
            // trace expansion. Recompute this choice only for changed inputs.
            let operation = [STRUCTURE, GAMMA, COLOR, EPSILON]
                .into_iter()
                .filter(|operation| eligible & operation != 0)
                .min_by_key(|operation| {
                    let counts = &plan.observed.candidates.counts;
                    match *operation {
                        STRUCTURE => 0,
                        GAMMA => 1 + counts[8..17].iter().sum::<usize>(),
                        COLOR => 2 + counts[1..4].iter().sum::<usize>(),
                        _ => 8 + counts[4].saturating_mul(counts[4]),
                    }
                })
                .unwrap();
            plan.pending &= !operation;
            let mut before = plan.value.clone();
            // Frontier facts describe one kernel invocation, not a permanent
            // limitation of this expression or a subsequent resume.
            before.proofs.frontier = ReductionStatus::Complete;
            let mut emitted;
            let (mut value, operation_status) = match request {
                ReductionRequest::Contract(settings) => {
                    let settings = ContractSettings {
                        order: if key.is_none() && plan.changes == 0 {
                            settings.order
                        } else {
                            None
                        },
                        ..settings
                    };
                    let (value, definitions, status) =
                        before.contract_domain(settings, &registry)?;
                    emitted = definitions;
                    (value, status)
                }
                ReductionRequest::Algebra(settings) => {
                    if operation == STRUCTURE {
                        let (value, definitions, status) =
                            before.contract_domain(structural, &registry)?;
                        emitted = definitions;
                        (value, status)
                    } else {
                        let definitions = available
                            .get_or_insert_with(|| registry.values().cloned().collect::<Vec<_>>());
                        let observed = &plan.observed;
                        let (prerequisites, additions) = before.contract_prerequisites(structural, definitions, |value| {
                            if alias_candidates.get(value.get_data()).is_some_and(|mask| mask & operation != 0) { return true; }
                            // A metric is a prerequisite, not the identity's seed.
                            // Nested leaves may lack individual regional records;
                            // the representation fallback must not select their
                            // metrics and hide the incident connection from the
                            // prerequisite contractor.
                            if matches!(value, AtomView::Fun(function) if function.get_symbol() == spenso::network::library::symbolic::ETS.metric) { return false; }
                            if let Some(region) = observed.region(value) {
                                return region.identity_candidates() & operation != 0;
                            }
                            match operation {
                                GAMMA => spenso::shadowing::TensorCollectFilter::Reps([crate::representations::Bispinor {}.into()]).matches(value),
                                COLOR => spenso::shadowing::TensorCollectFilter::Reps([crate::representations::ColorFundamental {}.into(), crate::representations::ColorAdjoint {}.into(), crate::representations::ColorSextet {}.into()]).matches(value),
                                _ => matches!(value, AtomView::Fun(function) if function.get_symbol() == *crate::epsilon::EPSILON_SYMBOL),
                            }
                        })?;
                        if prerequisites.root != before
                            || additions
                                .iter()
                                .any(|pair| registry.get(&pair.0.expression) != Some(pair))
                        {
                            let mut retained = additions
                                .into_iter()
                                .map(|pair| (pair.0.expression.clone(), pair))
                                .collect();
                            Self::retain_contracted_definitions(
                                &prerequisites,
                                &registry,
                                &mut retained,
                            )?;
                            emitted = retained.into_values().collect();
                            (prerequisites.root, prerequisites.status)
                        } else {
                            let (value, definitions) = match operation {
                                GAMMA => before.simplify_gamma_parts(definitions, settings.gamma.unwrap(), plan.settled.entry(GAMMA).or_default())?,
                                COLOR => before.simplify_color_parts(definitions, settings.color.unwrap(), plan.settled.entry(COLOR).or_default())?,
                                EPSILON if before.reduction_observations().epsilon_degree < 2
                                    && !before.reduction_observations().literals.iter().any(|literal| registry.get(literal).is_some_and(|(_, body)| !body.is_scalar())) => (before.clone(), Vec::new()),
                                EPSILON => before.collect_with_map(None, definitions,
                                    |value| matches!(value, AtomView::Fun(function) if function.get_symbol() == *crate::epsilon::EPSILON_SYMBOL),
                                    |selected, registry, complete| plan.settled.entry(EPSILON).or_default().run(selected, complete, |selected| {
                                        let expression = EpsilonSimplifierPass::step(selected.expression.as_view());
                                        let rewritten = selected.with_identity_result(expression)?;
                                        if rewritten == selected { return Ok((rewritten, Vec::new())); }
                                        let registry = registry.iter().cloned().map(|pair| (pair.0.expression.clone(), pair)).collect();
                                        let (mut value, definitions, status) = rewritten.contract_domain(structural, &registry)?;
                                        value.proofs.frontier = status;
                                        Ok((value, definitions))
                                    }))?,
                                _ => unreachable!("only enabled identities enter the algebra queue"),
                            };
                            emitted = definitions;
                            let status = value.proofs.frontier;
                            (value, status)
                        }
                    }
                }
            };
            // Kernels may return their available registry along with additions.
            // Unchanged associations are references, not replacement definitions.
            emitted.retain(|pair| {
                registry
                    .get(&pair.0.expression)
                    .is_none_or(|previous| previous != pair)
            });
            if value.structure != before.structure {
                return Err(TensorInferenceError::invalid(
                    "reduction changes its established logical interface",
                ));
            }
            let actual_change = value != before;
            if actual_change {
                value.proofs.literal_origin = before.proofs.literal_origin.clone();
                value.observe_replacement(&before);
                plan.observed = Arc::clone(value.reduction_observations());
                for dependency in &plan.dependencies {
                    if let Some(users) = dependents.get_mut(dependency) {
                        users.remove(&key);
                    }
                }
                plan.dependencies = plan
                    .observed
                    .literals
                    .iter()
                    .filter(|literal| registry.contains_key(*literal))
                    .cloned()
                    .collect();
                for dependency in &plan.dependencies {
                    dependents
                        .entry(dependency.clone())
                        .or_default()
                        .insert(key.clone());
                }
                plan.value = value;
                plan.changes += 1;
                plan.pending = enabled | if plan.prerequisites { STRUCTURE } else { 0 };
                changed = true;
                if let Some(handle) = &key {
                    registry.get_mut(handle).unwrap().1 = plan.value.clone();
                    available = None;
                }
            }
            // A completed family cannot discharge another family's exact
            // unfinished frontier. Only revisiting that operation can do so.
            if operation_status == ReductionStatus::Deferred {
                plan.deferred |= operation;
            } else {
                plan.deferred &= !operation;
            }
            plan.status = if operation_status == ReductionStatus::Capped {
                ReductionStatus::Capped
            } else if plan.deferred != 0 {
                ReductionStatus::Deferred
            } else {
                ReductionStatus::Complete
            };
            status = status.max(plan.status);
            if operation_status == ReductionStatus::Capped {
                plan.pending = 0;
            }
            if plan.pending != 0 && queued.insert(key.clone()) {
                queue.push_back(key.clone());
            }
            // Register all additions before deriving edges: a new definition
            // may refer to another addition later in the kernel's result.
            for pair in &mut emitted {
                // Seed the new domain before cloning it into the registry and
                // plan. Both owners then share its single initial observation.
                pair.1.reduction_observations();
                if pair.1.proofs.literal_origin.is_none()
                    && matches!(pair.0.expression.as_view(), AtomView::Fun(function) if function.get_symbol() == *TENSOR_ALIAS_SYMBOL)
                {
                    pair.1.proofs.literal_origin = Some(Arc::new(pair.1.clone()));
                }
                let handle = pair.0.expression.clone();
                if registry
                    .get(&handle)
                    .is_some_and(|previous| previous != pair)
                {
                    return Err(TensorInferenceError::invalid(
                        "conflicting reduction alias definitions",
                    ));
                }
                registry.insert(handle, pair.clone());
                available = None;
            }
            for pair in emitted {
                let handle = pair.0.expression.clone();
                if let std::collections::hash_map::Entry::Vacant(entry) =
                    domains.entry(Some(handle.clone()))
                {
                    let mut plan = DomainPlan::new(pair.1, &registry, enabled);
                    if operation_status == ReductionStatus::Capped {
                        plan.pending = 0;
                        plan.status = ReductionStatus::Capped;
                    }
                    for dependency in &plan.dependencies {
                        dependents
                            .entry(dependency.clone())
                            .or_default()
                            .insert(Some(handle.clone()));
                    }
                    alias_candidates.insert(handle.clone(), plan.effective);
                    entry.insert(plan);
                    if queued.insert(Some(handle.clone())) {
                        queue.push_back(Some(handle));
                    }
                    changed = true;
                }
            }
            {
                let plan = domains.get_mut(&key).unwrap();
                for literal in &plan.observed.literals {
                    if registry.contains_key(literal) && plan.dependencies.insert(literal.clone()) {
                        dependents
                            .entry(literal.clone())
                            .or_default()
                            .insert(key.clone());
                    }
                }
            }
            if actual_change
                && let Some(handle) = &key
                && let Some(users) = dependents.get(handle)
            {
                for user in users {
                    let plan = domains.get_mut(user).unwrap();
                    plan.settled.clear();
                    plan.pending = enabled | if plan.prerequisites { STRUCTURE } else { 0 };
                    if queued.insert(user.clone()) {
                        queue.push_back(user.clone());
                    }
                }
            }
        }
        status = domains
            .values()
            .map(|plan| plan.status)
            .max()
            .unwrap_or(status);
        let mut result = if changed {
            Arc::new(
                domains
                    .remove(&None)
                    .unwrap()
                    .value
                    .attach_aliases(registry.into_values(), Some(self))?,
            )
        } else {
            Arc::clone(self)
        };
        let proofs = &mut Arc::make_mut(&mut result).proofs;
        proofs.contracted = (!changed && self.proofs.contracted)
            || (matches!(request, ReductionRequest::Contract(settings) if settings.unrestricted())
                && status == ReductionStatus::Complete);
        let mut completed = if changed {
            Vec::new()
        } else {
            self.proofs
                .reduction
                .as_ref()
                .map_or_else(Vec::new, |record| record.completed.clone())
        };
        if status == ReductionStatus::Complete {
            completed.push(policy);
        }
        proofs.reduction = Some(ReductionRecord { completed, status });
        #[cfg(feature = "reference-cases")]
        crate::reference_cases::timing::outcome(diagnostic_operation, status);
        Ok(result)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::dirac::GammaOutput;
    use crate::shorthands::chain::Chain;
    use spenso::structure::representation::LibraryRep;
    use spenso::{g, mink, p};
    use symbolica::{
        atom::{Atom, AtomCore},
        symbol,
    };

    #[test]
    fn admitted_gamma_traces_do_not_revalidate_generated_results() {
        use super::super::inference::tests::{INFERENCE_CALLS, SCOPE_VALIDATIONS};
        use crate::gamma;
        use spenso::shadowing;
        use spenso::structure::{
            dimension::Dimension,
            representation::{Minkowski, RepName},
        };

        let reps = crate::test_support::test_initialize();
        let spin = reps.bis4.to_symbolic([]);
        for dimension in [
            Dimension::Concrete(4),
            Dimension::from(symbol!("trusted_trace_dimension")),
        ] {
            let representation = Minkowski {}.new_rep(dimension);
            for length in [4, 6, 8] {
                let indices = (0..length)
                    .map(|i| representation.to_symbolic([Atom::num(98400 + i)]))
                    .collect::<Vec<_>>();
                let source = SymbolicTensor::infer(shadowing::trace(
                    &spin,
                    indices.iter().map(|index| gamma!(index)),
                ))
                .unwrap();
                assert_eq!(source.proofs.validated.get(), Some(&true));
                INFERENCE_CALLS.with(|count| count.set(0));
                SCOPE_VALIDATIONS.with(|count| count.set(0));
                let result = source
                    .simplify_algebra(&AlgebraSettings {
                        gamma: Some(GammaSimplifySettings::default()),
                        epsilon: true,
                        ..Default::default()
                    })
                    .unwrap();
                assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
                assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
                assert_eq!(result.root().structure, source.structure);
                let repeated = result
                    .simplify_algebra(&AlgebraSettings {
                        gamma: Some(GammaSimplifySettings::default()),
                        epsilon: true,
                        ..Default::default()
                    })
                    .unwrap();
                assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
                assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
                assert!(Arc::ptr_eq(&result, &repeated));

                // Independent validation is a test oracle, outside the algebra.
                let resolved = result.resolved().unwrap();
                SymbolicTensor::validate_interface(&resolved.expression, &source.structure)
                    .unwrap();
                SymbolicTensor::validate_atom(&resolved.expression).unwrap();
            }
        }
    }

    #[test]
    fn compact_inner_product_gamma_reduces_both_raw_and_chain_operands() {
        use crate::gamma;
        use spenso::{chain, network::tags::SPENSO_TAG};
        use symbolica::function;
        let reps = crate::test_support::test_initialize();
        let left = reps.bis4.to_symbolic([Atom::num(91301)]);
        let inner = reps.bis4.to_symbolic([Atom::num(91302)]);
        let right = reps.bis4.to_symbolic([Atom::num(91303)]);
        let scalar_dot = function!(
            SPENSO_TAG.dot,
            p!(reps.mink4.to_symbolic([])),
            spenso::q!(reps.mink4.to_symbolic([]))
        );
        for dimension in [
            Atom::num(4),
            Atom::var(symbol!("compact_inner_product_dimension")),
        ] {
            let compact = function!(
                LibraryRep::from(spenso::structure::representation::Minkowski {}).symbol(),
                &dimension
            );
            let direct = function!(
                SPENSO_TAG.dot,
                gamma!(&left, &inner, &compact),
                gamma!(&inner, &right, &compact)
            );
            let chained = function!(
                SPENSO_TAG.dot,
                chain!(&left, &inner, gamma!(&compact)),
                chain!(&inner, &right, gamma!(&compact))
            );
            for expression in [direct, chained] {
                let source =
                    SymbolicTensor::<PartialStructure>::infer(&scalar_dot * expression).unwrap();
                let collected = source
                    .simplify_algebra(&AlgebraSettings {
                        gamma: Some(GammaSimplifySettings {
                            output: GammaOutput::Chains,
                            ..GammaSimplifySettings::default()
                        }),
                        ..Default::default()
                    })
                    .unwrap();
                assert_eq!(collected.root().structure, source.structure);
                let result = source
                    .simplify_algebra(&AlgebraSettings {
                        gamma: Some(GammaSimplifySettings::default()),
                        epsilon: true,
                        ..Default::default()
                    })
                    .unwrap();
                assert_eq!(result.root().structure, source.structure);
                assert_eq!(
                    result.resolved().unwrap().expression,
                    &dimension * &scalar_dot * g!(&left, &right)
                );
                let handle = source.alias_handle().unwrap();
                let aliased = Arc::new(
                    handle
                        .clone()
                        .with_aliases([(handle, source.clone())])
                        .unwrap(),
                );
                let alias_result = aliased
                    .simplify_algebra(&AlgebraSettings {
                        gamma: Some(GammaSimplifySettings::default()),
                        epsilon: true,
                        ..Default::default()
                    })
                    .unwrap();
                assert_eq!(alias_result.root().structure, source.structure);
                assert_eq!(
                    alias_result.resolved().unwrap().expression,
                    result.resolved().unwrap().expression
                );
                let repeated = collected
                    .simplify_algebra(&AlgebraSettings {
                        gamma: Some(GammaSimplifySettings::default()),
                        epsilon: true,
                        ..Default::default()
                    })
                    .unwrap();
                assert_eq!(
                    repeated.resolved().unwrap().expression,
                    result.resolved().unwrap().expression
                );
            }
        }
    }

    #[test]
    fn shared_pipeline_preserves_factored_spectators_and_typed_zero() {
        crate::test_support::test_initialize();
        let x = Atom::var(symbol!("shared_pipeline_x"));
        let y = Atom::var(symbol!("shared_pipeline_y"));
        let source = SymbolicTensor::<PartialStructure>::infer(
            (x + Atom::one()) * (y + Atom::one()) * p!(mink!(4, mu)),
        )
        .unwrap();
        let settings = AlgebraSettings::hep();
        assert_eq!(source.simplify_algebra(&settings).unwrap().root(), source);
        let zero = source.with_rewritten_expression(Atom::Zero).unwrap();
        assert_eq!(zero.simplify_algebra(&settings).unwrap().root(), zero);
    }

    #[test]
    fn shared_pipeline_preserves_exact_results_at_its_fixed_point_budget() {
        crate::test_support::test_initialize();
        let source = SymbolicTensor::<PartialStructure>::infer(
            g!(mink!(4, mu), mink!(4, nu)) * p!(mink!(4, nu)),
        )
        .unwrap();
        let disabled = source
            .simplify_algebra(&AlgebraSettings::default())
            .unwrap();
        assert_eq!(disabled.root(), source);
        assert_eq!(disabled.reduction_status(), ReductionStatus::Complete);
        let mut settings = ContractSettings::default();
        let result = Arc::new(source.contract(settings).unwrap());
        assert_eq!(result.root().expression, p!(mink!(4, mu)));
        assert_eq!(result.root().structure, source.structure);
        let rerun = result.contract(settings).unwrap();
        assert!(Arc::ptr_eq(&result, &rerun));
        assert_eq!(rerun.root(), result.root());
        assert_eq!(
            rerun.expression.get_aliases(),
            result.expression.get_aliases()
        );
        settings.max_passes = 1;
        let partial = source.contract(settings).unwrap();
        assert_eq!(partial.reduction_status(), ReductionStatus::Complete);
        assert_eq!(partial.root().expression, p!(mink!(4, mu)));
        assert_eq!(partial.root().structure, source.structure);
        settings.max_passes = 0;
        let untouched = source.contract(settings).unwrap();
        assert_eq!(untouched.reduction_status(), ReductionStatus::Capped);
        assert_eq!(untouched.root(), source);
    }

    #[test]
    fn shared_pipeline_registers_trace_alias_port_changes() {
        use crate::{gamma, test_support::test_initialize};
        use spenso::trace;

        let reps = test_initialize();
        let slots = (0..6)
            .map(|i| reps.mink_d.to_symbolic([Atom::num(93700 + i)]))
            .collect::<Vec<_>>();
        let trace = trace!(reps.bis4.to_symbolic([]);
            slots[..4].iter().map(|slot| gamma!(slot)));
        let traced = SymbolicTensor::<PartialStructure>::infer(trace.clone())
            .unwrap()
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::default()),
                epsilon: true,
                ..Default::default()
            })
            .unwrap();
        let metric = g!(&slots[0], &slots[4]);
        let input =
            SymbolicTensor::<PartialStructure>::infer(&metric * &traced.root().expression).unwrap();
        let source = Arc::new(
            input
                .clone()
                .with_aliases(traced.aliases().unwrap())
                .unwrap(),
        );
        let result = source
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::default()),
                epsilon: true,
                ..Default::default()
            })
            .unwrap()
            .contract(ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap();
        let reference_input = SymbolicTensor::<PartialStructure>::infer(metric * trace).unwrap();
        let reference = reference_input
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::default()),
                epsilon: true,
                ..Default::default()
            })
            .unwrap()
            .contract(ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(result.root().structure, input.structure);
        assert_eq!(reference.root().structure, reference_input.structure);
        // The two syntactic inputs establish different logical port orders.
        // Each result retains its own order; compare their algebra by labels.
        assert_eq!(
            result.expanded().unwrap().expression,
            reference.expanded().unwrap().expression
        );
        assert!(result.aliases().unwrap().len() > traced.aliases().unwrap().len());
    }

    #[test]
    fn shared_pipeline_checks_callback_rank_loss_inside_a_relabelled_alias() {
        use symbolica::atom::{AtomView, FunctionBuilder};

        crate::test_support::test_initialize();
        let a = mink!(4, 93801);
        let b = mink!(4, 93803);
        let target = b.clone();
        let head = spenso::tensor_symbol!(
            "shared_pipeline_alias_callback",
            norm = move |node, out| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **out = Atom::one();
                }
            }
        );
        let body = SymbolicTensor::<PartialStructure>::infer(
            FunctionBuilder::new(head).add_arg(&a).finish(),
        )
        .unwrap();
        let handle = body.alias_handle().unwrap();
        let root =
            SymbolicTensor::<PartialStructure>::infer(g!(&a, &b) * &handle.expression).unwrap();
        let source = Arc::new(root.with_aliases([(handle, body)]).unwrap());
        assert!(source.contract(ContractSettings::default()).is_err());
    }
    #[test]
    fn nested_compact_odd_trace_reduces_without_expanding_spectators() {
        use crate::gamma;
        use spenso::chain;

        let reps = crate::test_support::test_initialize();
        let left = reps.bis4.to_symbolic([Atom::num(93911)]);
        let right = reps.bis4.to_symbolic([Atom::num(93912)]);
        let momentum = p!(reps.mink4.to_symbolic([]));
        let spectator = (Atom::var(symbol!("nested_chain_spectator")) + Atom::one()).pow(20);
        let source = SymbolicTensor::<PartialStructure>::infer(
            &spectator
                * (Atom::num(3) + chain!(&left, &right, gamma!(&momentum)) * g!(&left, &right)),
        )
        .unwrap();
        let settings = AlgebraSettings {
            gamma: Some(GammaSimplifySettings::default()),
            ..Default::default()
        };
        let result = source.simplify_algebra(&settings).unwrap();
        // The trace of one gamma matrix vanishes independently of its momentum.
        assert_eq!(
            result.resolved().unwrap().expression,
            Atom::num(3) * spectator
        );
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert!(Arc::ptr_eq(
            &result,
            &result.simplify_algebra(&settings).unwrap()
        ));
    }

    #[test]
    fn shared_gamma_collection_combines_an_owned_alias_with_an_adjacent_word() {
        use crate::gamma;
        use spenso::{chain, q};

        let reps = crate::test_support::test_initialize();
        let left = reps.bis4.to_symbolic([Atom::num(93901)]);
        let right = reps.bis4.to_symbolic([Atom::num(93902)]);
        let momentum = reps.mink4.to_symbolic([]);
        let p = p!(&momentum);
        let q = q!(&momentum);
        let r = p!(1, &momentum);
        let body = SymbolicTensor::<PartialStructure>::infer(
            chain!(&left, &right, gamma!(&p)) + chain!(&left, &right, gamma!(&q)),
        )
        .unwrap();
        let handle = body.alias_handle().unwrap();
        let spectator = Atom::var(symbol!("shared_gamma_alias_spectator")) + Atom::one();
        let root = SymbolicTensor::<PartialStructure>::infer(
            &spectator * &handle.expression * chain!(&right, &left, gamma!(&r)),
        )
        .unwrap();
        let source = Arc::new(root.clone().with_aliases([(handle, body)]).unwrap());
        let collected = source
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings {
                    output: GammaOutput::Chains,
                    ..GammaSimplifySettings::default()
                }),
                ..Default::default()
            })
            .unwrap();
        assert_eq!(collected.root().structure, root.structure);
        let result = collected
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::default()),
                epsilon: true,
                ..Default::default()
            })
            .unwrap();
        let expected = Atom::num(4) * spectator * (g!(&p, &r) + g!(&q, &r));
        assert_eq!(result.expanded().unwrap().expression, expected.expand());
    }

    #[test]
    fn shared_color_collection_combines_an_owned_alias_with_an_adjacent_generator() {
        use crate::color_t;

        crate::test_support::test_initialize();
        let a_ik = color_t!([8, 94001], [3, 94003], [3, 94005]);
        let b_kj = color_t!([8, 94002], [3, 94005], [3, 94004]);
        let b_ik = color_t!([8, 94002], [3, 94003], [3, 94005]);
        let a_kj = color_t!([8, 94001], [3, 94005], [3, 94004]);
        let body = SymbolicTensor::<PartialStructure>::infer(a_ik * b_kj + b_ik * a_kj).unwrap();
        let handle = body.alias_handle().unwrap();
        let root = SymbolicTensor::<PartialStructure>::infer(
            &handle.expression * color_t!([8, 94002], [3, 94004], [3, 94006]),
        )
        .unwrap();
        let source = Arc::new(root.clone().with_aliases([(handle, body)]).unwrap());
        let result = source
            .simplify_algebra(&AlgebraSettings {
                color: Some(ColorSimplifySettings {
                    substitute_cof_dimension_invariants: true,
                    ..ColorSimplifySettings::default()
                }),
                ..Default::default()
            })
            .unwrap();
        assert_eq!(result.root().structure, root.structure);
        // Sum_b {T^a,T^b} T^b = (2 C_F - C_A/2) T^a = 7 T^a / 6.
        let expected = Atom::num(7) * color_t!([8, 94001], [3, 94003], [3, 94006]) / Atom::num(6);
        // The color output may retain its one-generator chain. Compare in the
        // existing explicit tensor representation, preserving the exact coefficient.
        assert_eq!(
            result.expanded().unwrap().expression.undo_single_length(),
            expected
        );
    }
    #[test]
    fn massive_chiral_current_product_registers_all_literal_aliases() {
        use spenso::structure::partial::PartialStructureExt;

        crate::test_support::test_initialize();
        let _ = (
            spenso::tensor_symbol!("python::p1"),
            spenso::tensor_symbol!("python::p2"),
            spenso::tensor_symbol!("python::k1"),
            spenso::tensor_symbol!("python::k2"),
        );
        // Exact factorized spin-sum inputs from the installed FeynCalc tree test.
        let incoming = Atom::parse(
            r#"1/4*(UFO::Me*spenso::g(spenso::bis(4,python::i2),spenso::bis(4,python::i3))
    +spenso::gamma(spenso::bis(4,python::i2),spenso::bis(4,python::i3),
        python::p1(spenso::mink(4))
    )
)*(-UFO::Me*spenso::g(spenso::bis(4,python::i0),spenso::bis(4,python::i1))
        +spenso::gamma(spenso::bis(4,python::i0),spenso::bis(4,python::i1),
            python::p2(spenso::mink(4))
        )
    )*spenso::chain(spenso::bis(4,python::i1),spenso::bis(4,python::i2),
        spenso::gamma(spenso::in,spenso::out,spenso::mink(4,python::mu)),
        spenso::projp(spenso::in,spenso::out)
    )*spenso::chain(spenso::bis(4,python::i3),spenso::bis(4,python::i0),
        spenso::gamma(spenso::in,spenso::out,spenso::mink(4,python::nu)),
        spenso::projp(spenso::in,spenso::out)
    )"#,
            "python",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap();
        let outgoing = Atom::parse(
            r#"(UFO::MM*spenso::g(spenso::bis(4,python::j0),spenso::bis(4,python::j1))
    +spenso::gamma(spenso::bis(4,python::j0),spenso::bis(4,python::j1),
        python::k1(spenso::mink(4))
    )
)*(-UFO::MM*spenso::g(spenso::bis(4,python::j2),spenso::bis(4,python::j3))
        +spenso::gamma(spenso::bis(4,python::j2),spenso::bis(4,python::j3),
            python::k2(spenso::mink(4))
        )
    )*spenso::chain(spenso::bis(4,python::j1),spenso::bis(4,python::j2),
        spenso::gamma(spenso::in,spenso::out,spenso::mink(4,python::mu)),
        spenso::projp(spenso::in,spenso::out)
    )*spenso::chain(spenso::bis(4,python::j3),spenso::bis(4,python::j0),
        spenso::gamma(spenso::in,spenso::out,spenso::mink(4,python::nu)),
        spenso::projp(spenso::in,spenso::out)
    )"#,
            "python",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap();
        let product = &incoming * &outgoing;
        let masses = [symbol!("UFO::Me"), symbol!("UFO::MM")];
        let massless = product.replace_map(|value, _, output| {
            if matches!(value, symbolica::atom::AtomView::Var(var) if masses.contains(&var.get_symbol())) {
                **output = Atom::Zero;
            }
        });
        for (name, rank, input) in [
            ("incoming", 2, incoming),
            ("outgoing", 2, outgoing),
            ("product", 0, product),
            ("massless", 0, massless),
        ] {
            let source = SymbolicTensor::<PartialStructure>::infer(input).unwrap();
            assert_eq!(source.structure.logical_slots().len(), rank, "{name}");
            let result = source
                .simplify_algebra(&AlgebraSettings {
                    gamma: Some(GammaSimplifySettings::default()),
                    epsilon: true,
                    ..Default::default()
                })
                .unwrap_or_else(|error| panic!("{name}: {error}"));
            assert_eq!(result.root().structure, source.structure, "{name}");
            assert!(
                result.resolved().is_ok(),
                "{name}: all literal aliases must resolve"
            );
        }
    }
}

#[cfg(test)]
mod planner_tests;

#[cfg(test)]
mod completion_tests;
