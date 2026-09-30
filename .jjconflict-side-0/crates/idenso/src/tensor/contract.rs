use spenso::structure::{
    partial::{PartialStructure, PartialStructureExt},
    representation::{LibraryRep, RepName},
};
use std::{collections::HashMap, sync::Arc};
#[cfg(test)]
use symbolica::atom::AtomCore;
use symbolica::atom::{AliasedAtom, Atom, AtomView};

use super::{SymbolicTensor, aliases::AliasInterfaces, inference::TensorInferenceError};
use crate::shorthands::schoonschip::{
    DotNormalizer, ReductionStatus, SchoonschipSettings, SchoonschipWithSettings,
    SimplificationCandidates, SlotContraction,
};

#[cfg(test)]
thread_local! {
    pub(crate) static CONTRACT_DOMAIN_CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Contraction scheduling and whether rank-one tensors may form compact products.
/// Metric-only contraction retains explicit vector endpoints, including Vakint inputs.
#[derive(Clone, Copy, Debug)]
pub struct ContractSettings<'a> {
    /// Optional permutation of normalized factors for the initial contraction.
    pub order: Option<&'a [usize]>,
    /// Select connections by representation; `None` permits all and `Some(&[])` none.
    /// Construction's intrinsic normalization, including symmetric scalar
    /// products, has already happened and is not undone by this filter.
    pub representations: Option<&'a [LibraryRep]>,
    /// Substitute compatible metric/identity ports and form dimension factors.
    pub metrics: bool,
    /// Substitute vector ports and form normalized scalar products.
    pub rank_one: bool,
    /// Collect connected matrix factors while preserving order and orientation.
    pub collect_chains: bool,
    /// Close compatible matrix chains into trace notation without evaluating it.
    pub collect_traces: bool,
    /// Maximum successful planning rounds per domain; unfinished results stay exact.
    pub max_passes: usize,
}

impl Default for ContractSettings<'_> {
    fn default() -> Self {
        Self {
            order: None,
            representations: None,
            metrics: true,
            rank_one: true,
            collect_chains: true,
            collect_traces: true,
            max_passes: 16,
        }
    }
}

impl<'a> ContractSettings<'a> {
    pub(crate) fn permits(self, representation: LibraryRep) -> bool {
        self.representations.is_none_or(|allowed| {
            allowed
                .iter()
                .any(|candidate| candidate.base() == representation.base())
        })
    }

    pub(crate) fn unrestricted(self) -> bool {
        self.representations.is_none() && self.metrics && self.rank_one
    }

    /// Visit normalized top-level factors in this explicit permutation.
    pub fn with_order(mut self, order: &'a [usize]) -> Self {
        self.order = Some(order);
        self
    }

    /// Contract metrics while retaining explicit rank-one tensor endpoints.
    pub fn without_rank_one_tensors(mut self) -> Self {
        self.rank_one = false;
        self
    }
}

impl SymbolicTensor<PartialStructure> {
    /// A cheap admission predicate from the existing domain observation. A true
    /// result is only a candidate; the graph still establishes actual eligibility.
    pub(crate) fn structural_work_permitted(&self, settings: ContractSettings<'_>) -> bool {
        let observed = self.reduction_observations();
        let candidates = &observed.candidates;
        if (candidates.dots || (settings.metrics && candidates.symbols[5]))
            && self.normalized_contract_notation(settings) != self.expression
        {
            return true;
        }
        if settings
            .representations
            .is_some_and(<[LibraryRep]>::is_empty)
            || (!settings.metrics
                && !settings.rank_one
                && !settings.collect_chains
                && !settings.collect_traces)
        {
            return false;
        }
        if !settings.collect_chains
            && !settings.collect_traces
            && !candidates.brackets
            && observed.excludes_contraction_sources(settings)
        {
            return false;
        }
        if candidates.complete && candidates.traversal_complete {
            if !settings.collect_chains
                && !settings.collect_traces
                && !(settings.metrics && (candidates.symbols[17] || candidates.symbols[5]))
                && !(settings.rank_one && candidates.dots)
                && !candidates.brackets
                && !observed.literals.iter().any(|literal| {
                    matches!(literal.as_view(), AtomView::Fun(function)
                        if function.get_symbol() == *super::aliases::TENSOR_ALIAS_SYMBOL)
                })
            {
                return false;
            }
            if !observed
                .representations
                .iter()
                .any(|rep| settings.permits(*rep))
            {
                return false;
            }
            return candidates.repeated_indices || candidates.brackets;
        }
        true
    }

    fn contraction_is_noop(&self) -> bool {
        let observed = self.reduction_observations();
        self.contraction_is_noop_observed(Default::default(), &observed.candidates)
    }

    fn contraction_is_noop_observed<const N: usize>(
        &self,
        settings: ContractSettings<'_>,
        observed: &SimplificationCandidates<N>,
    ) -> bool {
        observed.complete
            && observed.intrinsic
            && observed.traversal_complete
            && !observed.repeated_indices
            && !observed.brackets
            // A no-op leaf needs its exact validated interface, not the
            // stronger self-dual proof required for arbitrary algebraic rewrites.
            && (self.rewrites_preserve_interface()
                || (matches!(self.expression.as_view(), AtomView::Fun(_))
                    && self.structure.open_positions().is_empty()
                    && !self.expression.as_view().needs_normalization()
                    && self.established_interface_is_valid()))
            // The observer conservatively marks explicit vector leaves as dot
            // candidates. Only the existing normalizer can discharge that flag;
            // intrinsic observation above excludes user callbacks from this check.
            && (!observed.dots || DotNormalizer::with_settings(self.expression.as_view(), settings) == self.expression)
    }

    /// Contract the metric/vector graph while keeping generated sums in aliases.
    /// An explicit order permutes the top-level factors after dot normalization.
    /// Materializing the result is a separate operation on the returned tensor.
    pub fn contract(
        &self,
        settings: ContractSettings<'_>,
    ) -> Result<SymbolicTensor<AliasInterfaces, AliasedAtom>, TensorInferenceError> {
        self.ensure_validated()?;
        self.reduction_observations();
        Arc::new(self.clone().with_aliases([])?)
            .plan_reduction(super::simplification::ReductionRequest::Contract(settings))
            .map(Arc::unwrap_or_clone)
    }

    /// Finish one domain, exposing only literal boundaries selected by the
    /// existing collector. Registry traversal and scheduling belong to the planner.
    pub(crate) fn contract_domain(
        &self,
        settings: ContractSettings<'_>,
        registry: &HashMap<Atom, super::aliases::Definition>,
    ) -> Result<(Self, Vec<super::aliases::Definition>, ReductionStatus), TensorInferenceError>
    {
        #[cfg(test)]
        CONTRACT_DOMAIN_CALLS.with(|count| count.set(count.get() + 1));
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::Contraction,
        );
        let contracted = self.contract_parts(settings)?;
        let mut status = contracted.status;
        let mut definitions = HashMap::new();
        SymbolicTensor::<AliasInterfaces, AliasedAtom>::retain_contracted_definitions(
            &contracted,
            registry,
            &mut definitions,
        )?;
        let mut root = contracted.root;
        root.observe_replacement(self);
        if status != ReductionStatus::Capped
            && root.reduction_observations().candidates.repeated_indices
            && root
                .reduction_observations()
                .literals
                .iter()
                .any(|literal| {
                    definitions
                        .get(literal)
                        .or_else(|| registry.get(literal))
                        .is_some_and(|(_, body)| !body.is_scalar())
                })
        {
            use spenso::{
                network::{library::symbolic::ETS, tags::SPENSO_TAG},
                structure::{
                    abstract_index::AbstractIndex,
                    slot::{IsAbstractSlot, SlotMatcher},
                },
            };
            // Only selected exposure needs the collector's mutable literal
            // registry. Ordinary domain contraction borrows the planner's
            // shared registry and reports only new associations.
            let available = registry
                .values()
                .chain(definitions.values())
                .cloned()
                .collect::<Vec<_>>();
            let mut slots = SlotMatcher::default();
            let mut completed = false;
            let (rewritten, emitted) = root.collect_with_map(
                Some(&mut |_, _, _| { completed = true; Ok(()) }),
                &available,
                |value| {
                    let AtomView::Fun(function) = value else { return false; };
                    let mut ports = HashMap::<LibraryRep, usize>::new();
                    for argument in function.iter() {
                        if let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(argument)
                            && settings.permits(slot.rep().rep)
                        {
                            *ports.entry(slot.rep().rep.base()).or_default() += 1;
                        }
                    }
                    !ports.is_empty() && ((settings.metrics && function.get_symbol() == ETS.metric)
                        || (settings.rank_one && function.get_symbol().has_tag(&SPENSO_TAG.rank1))
                        || ((settings.collect_chains || settings.collect_traces)
                            && (function.get_symbol().has_tag(&SPENSO_TAG.tensor)
                                || function.get_symbol() == SPENSO_TAG.chain)
                            && ports.values().any(|&count| count >= 2)))
                },
                |selected, available, _| {
                    selected.observe_replacement(&root);
                    let contracted = selected.contract_parts(ContractSettings { order: None, ..settings })?;
                    status = status.max(contracted.status);
                    let mut definitions = HashMap::new();
                    if !contracted.aliases.is_empty() || !contracted.literal_relabellings.is_empty() {
                        let inherited = available.iter().cloned().map(|pair| (pair.0.expression.clone(), pair)).collect::<HashMap<_, _>>();
                        SymbolicTensor::<AliasInterfaces, AliasedAtom>::retain_contracted_definitions(&contracted, &inherited, &mut definitions)?;
                    }
                    Ok((contracted.root, definitions.into_values().collect()))
                },
            )?;
            root = rewritten;
            for pair in emitted {
                definitions.insert(pair.0.expression.clone(), pair);
            }
            if !completed {
                status = status.max(ReductionStatus::Deferred);
            } else if status != ReductionStatus::Capped
                && root.proofs.frontier == ReductionStatus::Complete
            {
                root.observe_replacement(self);
                if !root.structural_work_permitted(settings) {
                    // Exposure has discharged the formerly opaque boundary.
                    // Alias bodies retain independent planner completion.
                    status = ReductionStatus::Complete;
                }
            }
        }
        Ok((root, definitions.into_values().collect(), status))
    }

    /// Contract only metric/vector paths incident to selected identity factors.
    /// Unselected branches retain their exact payload and normalization state.
    pub(crate) fn contract_prerequisites(
        &self,
        settings: ContractSettings<'_>,
        definitions: &[super::aliases::Definition],
        mut selected: impl FnMut(AtomView<'_>) -> bool,
    ) -> Result<
        (
            crate::shorthands::schoonschip::FactorizedContraction<Self>,
            Vec<super::aliases::Definition>,
        ),
        TensorInferenceError,
    > {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::Contraction,
        );
        if !self.reduction_observations().candidates.repeated_indices
            && !self.structural_work_permitted(settings)
        {
            // Free-index results need no prerequisite graph. A repeated alias
            // connection remains a candidate even when its hidden definition
            // supplies the metric/vector head that authorizes local exposure.
            return Ok((
                crate::shorthands::schoonschip::FactorizedContraction {
                    root: self.clone(),
                    aliases: Vec::new(),
                    status: ReductionStatus::Complete,
                    literal_relabellings: Vec::new(),
                },
                Vec::new(),
            ));
        }
        let mut emitted = Vec::new();
        let contracted = self.contract_prerequisite_scopes(
            settings,
            definitions,
            &mut emitted,
            &mut selected,
            0,
        )?;
        Ok((contracted, emitted))
    }

    fn contract_prerequisite_scopes(
        &self,
        settings: ContractSettings<'_>,
        definitions: &[super::aliases::Definition],
        emitted: &mut Vec<super::aliases::Definition>,
        selected: &mut dyn FnMut(AtomView<'_>) -> bool,
        depth: usize,
    ) -> Result<crate::shorthands::schoonschip::FactorizedContraction<Self>, TensorInferenceError>
    {
        use super::inference::InterfaceInference;
        use crate::shorthands::schoonschip::FactorizedContraction;
        use spenso::structure::TensorStructure;
        let mut result = FactorizedContraction {
            root: self.clone(),
            aliases: Vec::new(),
            status: ReductionStatus::Complete,
            literal_relabellings: Vec::new(),
        };
        if depth >= spenso::shadowing::TermTape::<usize>::MAX_EXPANSION_DEPTH {
            result.status = ReductionStatus::Deferred;
            return Ok(result);
        }
        let emitted_before = emitted.len();
        result = self.contract_prerequisite_boundary(settings, definitions, emitted, selected)?;
        if result.root != *self
            || !result.aliases.is_empty()
            || !result.literal_relabellings.is_empty()
            || emitted.len() != emitted_before
            || result.status == ReductionStatus::Capped
        {
            // Let the shared planner reconsider the changed outer boundary
            // before any selected child is opened.
            return Ok(result);
        }
        // Only an unchanged outer boundary admits work inside a selected
        // arithmetic child. Additive branches and scalar power bases own their
        // independent dummy scopes and remain factorized.
        let power_selected = match self.expression.as_view() {
            AtomView::Pow(power) => selected(power.get_base()),
            _ => true,
        };
        let mut visit = |value: AtomView<'_>,
                         interface: Option<PartialStructure>|
         -> Result<Atom, TensorInferenceError> {
            if !selected(value) {
                return Ok(value.to_owned());
            }
            let interface = match interface {
                Some(interface) => interface,
                None => InterfaceInference::replacement_interface(value)?,
            };
            let domain = if self.normalization_is_intrinsic() {
                Self::from_validated_parts(value.to_owned(), interface)
            } else {
                Self::checked_parts(value.to_owned(), interface)?
            };
            domain.observe_replacement(self);
            let nested = domain.contract_prerequisite_scopes(
                settings,
                definitions,
                emitted,
                selected,
                depth + 1,
            )?;
            result.aliases.extend(nested.aliases);
            result
                .literal_relabellings
                .extend(nested.literal_relabellings);
            result.status = result.status.max(nested.status);
            Ok(nested.root.expression)
        };
        let expression = match self.expression.as_view() {
            AtomView::Add(sum) => Atom::add_many(
                sum.iter()
                    .map(|term| visit(term, Some(self.structure.clone())))
                    .collect::<Result<Vec<_>, _>>()?,
            ),
            AtomView::Pow(power) if power_selected => {
                let (base, exponent) = power.get_base_exp();
                let interface = InterfaceInference::replacement_interface(base)?;
                if interface.canonical().is_scalar() {
                    use symbolica::atom::AtomCore;
                    visit(base, Some(interface))?.pow(exponent)
                } else {
                    self.expression.clone()
                }
            }
            AtomView::Mul(product) => {
                let mut factors = Vec::new();
                for factor in product.iter() {
                    factors.push(
                        if matches!(
                            factor,
                            AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_)
                        ) {
                            visit(factor, None)?
                        } else {
                            factor.to_owned()
                        },
                    );
                }
                Atom::mul_many(factors)
            }
            _ => self.expression.clone(),
        };
        result.root = self.with_identity_result(expression)?;
        Ok(result)
    }

    fn contract_prerequisite_boundary(
        &self,
        settings: ContractSettings<'_>,
        definitions: &[super::aliases::Definition],
        emitted: &mut Vec<super::aliases::Definition>,
        selected: &mut dyn FnMut(AtomView<'_>) -> bool,
    ) -> Result<crate::shorthands::schoonschip::FactorizedContraction<Self>, TensorInferenceError>
    {
        let mut result = crate::shorthands::schoonschip::FactorizedContraction {
            root: self.clone(),
            aliases: Vec::new(),
            status: ReductionStatus::Complete,
            literal_relabellings: Vec::new(),
        };
        let contractor = SlotContraction::configured(settings.metrics, settings.representations)
            .observed(result.root.reduction_observations());
        let contractor = match result.root.shallow_graph() {
            Ok(graph) => contractor.planned(
                result.root.expression.clone(),
                graph,
                result.root.proofs.leaf_interfaces.get().cloned(),
            ),
            Err(_) => contractor,
        };
        fn supports_structural(
            literal: &Atom,
            definitions: &HashMap<Atom, &SymbolicTensor<PartialStructure>>,
            cache: &mut HashMap<Atom, bool>,
            depth: usize,
        ) -> bool {
            if let Some(&candidate) = cache.get(literal) {
                return candidate;
            }
            if depth >= spenso::shadowing::TermTape::<usize>::MAX_EXPANSION_DEPTH {
                return false;
            }
            let Some(body) = definitions.get(literal) else {
                return false;
            };
            let observed = body.reduction_observations();
            let candidate = observed.candidates.dots
                || observed.candidates.symbols[17]
                || observed.literals.iter().any(|dependency| {
                    supports_structural(dependency, definitions, cache, depth + 1)
                });
            cache.insert(literal.clone(), candidate);
            candidate
        }
        let definitions_by_handle = definitions
            .iter()
            .map(|(handle, body)| (handle.expression.clone(), body))
            .collect::<HashMap<_, _>>();
        let mut cache = HashMap::new();
        let candidates = result
            .root
            .reduction_observations()
            .literals
            .iter()
            .filter(|literal| supports_structural(literal, &definitions_by_handle, &mut cache, 0))
            .cloned()
            .collect::<std::collections::HashSet<_>>();
        let Some((expression, interface, mut spectators, expose, incident_ports)) = contractor
            .prerequisite_factors(result.root.expression.as_view(), &mut *selected, |value| {
                candidates.contains(value.get_data())
            })
        else {
            return Ok(result);
        };
        let domain = if self.normalization_is_intrinsic() {
            Self::from_validated_parts(expression, interface)
        } else {
            Self::checked_parts(expression, interface)?
        };
        domain.observe_replacement(self);
        let policy = ContractSettings {
            collect_chains: false,
            collect_traces: false,
            order: None,
            ..settings
        };
        if expose {
            use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
            let mut available = definitions.to_vec();
            available.extend(emitted.iter().cloned());
            let incident_ports = incident_ports
                .into_iter()
                .collect::<std::collections::HashSet<_>>();
            let (root, additional) = domain.collect_with_map(None, &available,
                |value| selected(value) || matches!(value, AtomView::Fun(function)
                    if (function.get_symbol() == ETS.metric || function.get_symbol().has_tag(&SPENSO_TAG.rank1))
                        && function.iter().any(|argument| incident_ports.contains(argument.get_data()))),
                |part, available, _| {
                    part.observe_replacement(&domain);
                    let contracted = part.contract_parts(policy)?;
                    result.status = result.status.max(contracted.status);
                    let mut emitted = HashMap::new();
                    if !contracted.aliases.is_empty() || !contracted.literal_relabellings.is_empty() {
                        let inherited = available.iter().cloned().map(|pair| (pair.0.expression.clone(), pair)).collect::<HashMap<_,_>>();
                        SymbolicTensor::<AliasInterfaces, AliasedAtom>::retain_contracted_definitions(&contracted, &inherited, &mut emitted)?;
                    }
                    Ok((contracted.root, emitted.into_values().collect()))
                })?;
            result.status = result.status.max(root.proofs.frontier);
            emitted.extend(additional);
            spectators.push(root.expression);
            result.root = self.with_identity_result(Atom::mul_many(spectators))?;
            return Ok(result);
        }
        let contracted = domain.contract_parts(ContractSettings {
            collect_chains: false,
            collect_traces: false,
            order: None,
            ..settings
        })?;
        result.aliases.extend(contracted.aliases);
        result
            .literal_relabellings
            .extend(contracted.literal_relabellings);
        result.status = result.status.max(contracted.status);
        spectators.push(contracted.root.expression);
        result.root = self.with_identity_result(Atom::mul_many(spectators))?;
        Ok(result)
    }

    pub(crate) fn contract_parts(
        &self,
        settings: ContractSettings<'_>,
    ) -> Result<crate::shorthands::schoonschip::FactorizedContraction<Self>, TensorInferenceError>
    {
        let observed = self.reduction_observations();
        self.contract_parts_observed(settings, &observed.candidates)
    }

    pub(crate) fn contract_parts_observed<const N: usize>(
        &self,
        settings: ContractSettings<'_>,
        observed: &SimplificationCandidates<N>,
    ) -> Result<crate::shorthands::schoonschip::FactorizedContraction<Self>, TensorInferenceError>
    {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::Contraction,
        );
        let mut result = self.contract_kernel(settings, observed)?;
        if settings.rank_one && observed.dots {
            let expression =
                SlotContraction::configured(settings.metrics, settings.representations)
                    .contract_inner_products(
                        result.root.expression.as_view(),
                        &mut result.literal_relabellings,
                    );
            result.root = result.root.with_identity_result(expression)?;
        }
        if settings.collect_chains || settings.collect_traces {
            result.root.observe_replacement(self);
            let mut representations = self
                .reduction_observations()
                .representations
                .iter()
                .copied()
                .filter(|rep| settings.permits(*rep))
                .map(|rep| rep.base())
                .collect::<Vec<_>>();
            representations.sort_unstable();
            representations.dedup();
            if representations.is_empty() {
                return Ok(result);
            }
            let contractor =
                SlotContraction::configured(settings.metrics, settings.representations)
                    .observed(result.root.reduction_observations());
            let contractor = match result.root.shallow_graph() {
                Ok(graph) => contractor.planned(
                    result.root.expression.clone(),
                    graph,
                    result.root.proofs.leaf_interfaces.get().cloned(),
                ),
                Err(_) => contractor,
            };
            let expression = contractor
                .collect_matrix_notation(
                    result.root.expression.as_view(),
                    &representations,
                    settings.collect_chains,
                    settings.collect_traces,
                )
                .unwrap_or_else(|| result.root.expression.clone());
            result.root = result.root.with_identity_result(expression)?;
        }
        Ok(result)
    }

    fn normalized_contract_notation(&self, settings: ContractSettings<'_>) -> Atom {
        let mut source = self.expression.clone();
        let observations = self.reduction_observations();
        if settings.metrics && observations.candidates.symbols[5] {
            use crate::shorthands::chain::Chain;
            for representation in observations
                .representations
                .iter()
                .copied()
                .filter(|rep| settings.permits(*rep))
            {
                source = source.normalize_chain_identities(Some(representation.base()));
            }
        }
        if observations.candidates.dots {
            DotNormalizer::with_settings(source.as_view(), settings)
        } else {
            source
        }
    }

    fn contract_kernel<const N: usize>(
        &self,
        settings: ContractSettings<'_>,
        observed: &SimplificationCandidates<N>,
    ) -> Result<crate::shorthands::schoonschip::FactorizedContraction<Self>, TensorInferenceError>
    {
        let order = settings.order;
        let source = self.normalized_contract_notation(settings);
        if let Some(order) = order {
            let count = match source.as_view() {
                AtomView::Mul(product) => product.iter().len(),
                _ => 1,
            };
            let mut positions = order.to_vec();
            positions.sort_unstable();
            if positions != (0..count).collect::<Vec<_>>() {
                return Err(TensorInferenceError::Invalid(
                    "contraction order must list every normalized top-level factor exactly once"
                        .into(),
                ));
            }
        }
        if (settings
            .representations
            .is_some_and(<[LibraryRep]>::is_empty)
            || (!settings.metrics && !settings.rank_one))
            || (order.is_none()
                && source == self.expression
                && self.contraction_is_noop_observed(settings, observed))
        {
            // Complete intrinsic observation and the carried interface proof
            // exclude callbacks, opaque payloads and AUTO. The alias owner
            // still visits every original reachable definition.
            return Ok(crate::shorthands::schoonschip::FactorizedContraction {
                root: self.with_identity_result(source)?,
                aliases: Vec::new(),
                status: ReductionStatus::Complete,
                literal_relabellings: Vec::new(),
            });
        }
        let contractor = SlotContraction::configured(settings.metrics, settings.representations)
            .observed(self.reduction_observations());
        let contractor = match self.shallow_graph() {
            Ok(graph) => contractor.planned(
                self.expression.clone(),
                graph,
                self.proofs.leaf_interfaces.get().cloned(),
            ),
            Err(_) => contractor,
        };
        if let Some(contracted) =
            contractor.contract_factorized(source.as_view(), order, settings.rank_one)
        {
            return Ok(crate::shorthands::schoonschip::FactorizedContraction {
                root: self.with_identity_result(contracted.root)?,
                aliases: contracted.aliases,
                status: contracted.status,
                literal_relabellings: contracted.literal_relabellings,
            });
        }
        // Callback-sensitive leaves retain the existing ordered rewrite and
        // checked finisher. A changed result can expose intrinsic work (for
        // example alpha-equivalent dummy labels after metadata cleanup); retry
        // the same graph owner only on that changed input.
        // Preserve the established network sum schedule: each already-present
        // term is rewritten once, then collected. A nested whole-sum discovery
        // probe would normalize its first callback-sensitive term twice.
        // The typed interface includes ports inside chain/trace words. Their
        // builtin projectors can decline intrinsic graph admission; the ordered
        // reducer must still apply substitutions to those established ports.
        let mut rewrite = SchoonschipSettings::default().with_chain_like_functions();
        rewrite.schoonschip_rank1_tensors = settings.rank_one;
        rewrite.metrics = settings.metrics;
        rewrite.representations = settings.representations.map(<[LibraryRep]>::to_vec);
        let mut literal_relabellings = Vec::new();
        let rewritten = match source.as_view() {
            AtomView::Add(sum) => Atom::add_many(sum.iter().map(|term| {
                SchoonschipWithSettings { settings: &rewrite }.run_observed(
                    term,
                    &mut literal_relabellings,
                    Some(observed),
                )
            })),
            _ => SchoonschipWithSettings { settings: &rewrite }.run_observed(
                source.as_view(),
                &mut literal_relabellings,
                Some(observed),
            ),
        };
        if rewritten != source
            && let Some(contracted) =
                contractor.contract_factorized(rewritten.as_view(), None, settings.rank_one)
        {
            literal_relabellings.extend(contracted.literal_relabellings);
            return Ok(crate::shorthands::schoonschip::FactorizedContraction {
                root: self.with_identity_result(contracted.root)?,
                aliases: contracted.aliases,
                status: contracted.status,
                literal_relabellings,
            });
        }
        let root = self.with_identity_result(rewritten)?;
        root.observe_replacement(self);
        // The ordered fallback supports spaces declined by the coefficient
        // collector. Repeated gamma ports or compact slash vectors are not
        // structural source endpoints. Their existing regional observations can
        // certify completion without a new walk; explicit metric/vector sources
        // inside unresolved sums remain pending. An incomplete collector frontier
        // returned above retains its budget.
        let complete = !root.structural_work_permitted(settings)
            || root
                .reduction_observations()
                .excludes_contraction_sources(settings)
            || (root.expression != self.expression && root.contraction_is_noop());
        Ok(crate::shorthands::schoonschip::FactorizedContraction {
            root,
            aliases: Vec::new(),
            status: if complete {
                ReductionStatus::Complete
            } else {
                ReductionStatus::Deferred
            },
            literal_relabellings,
        })
    }
}

impl SymbolicTensor<AliasInterfaces, AliasedAtom> {
    /// Contract the root and its reachable definition templates without resolving
    /// their DAG. A completed default result reuses the same sealed allocation.
    /// Disconnected definitions are preserved without invoking their callbacks.
    pub fn contract(
        self: &Arc<Self>,
        settings: ContractSettings<'_>,
    ) -> Result<Arc<Self>, TensorInferenceError> {
        self.plan_reduction(super::simplification::ReductionRequest::Contract(settings))
    }

    pub(crate) fn retain_contracted_definitions(
        contracted: &crate::shorthands::schoonschip::FactorizedContraction<
            SymbolicTensor<PartialStructure>,
        >,
        inherited: &HashMap<Atom, super::aliases::Definition>,
        definitions: &mut HashMap<Atom, super::aliases::Definition>,
    ) -> Result<(), TensorInferenceError> {
        for (source, target) in &contracted.literal_relabellings {
            Self::register_literal_use(
                source,
                target,
                &|key| inherited.get(key).cloned(),
                definitions,
            )?;
        }
        for (handle, body) in &contracted.aliases {
            let scalar = PartialStructure::from_logical_slots([]);
            // The factorized kernel emits only proven scalar coefficient
            // definitions; callback-sensitive fallback emits no new aliases.
            let pair = (
                SymbolicTensor::from_validated_parts(handle.clone(), scalar.clone()),
                SymbolicTensor::from_validated_parts(body.clone(), scalar),
            );
            if let Some(previous) = definitions.get(handle).or_else(|| inherited.get(handle)) {
                if previous != &pair {
                    return Err(TensorInferenceError::invalid(
                        "conflicting contraction alias definitions",
                    ));
                }
                // An unchanged association retains its interface and domain
                // observations, even when the kernel reports it again.
                continue;
            }
            definitions.insert(handle.clone(), pair);
        }
        Ok(())
    }
}

#[cfg(test)]
mod alias_tests;
#[cfg(test)]
mod factorized_tests;
#[cfg(test)]
mod network_tests;

#[cfg(test)]
pub mod test {
    use super::super::AbstractIndex;
    use super::*;
    use crate::shorthands::schoonschip::Schoonschip;
    use crate::test_support::test_initialize;
    use spenso::{g, mink, p};

    #[test]
    fn compact_vector_dot_exposes_closed_spinor_connection() {
        use crate::{bis, gamma, tensor::AlgebraSettings};
        use spenso::{dot, q, trace};
        use symbolica::atom::AtomCore;
        let reps = test_initialize();
        let p = p!(mink!(4));
        let q = q!(mink!(4));
        // Multiplying a rank-one composite by a vector already constructs this
        // dot, before contract is called. Its spinor connection remains internal.
        let source = SymbolicTensor::infer(dot!(
            gamma!(bis!(4, 831), bis!(4, 832), mink!(4, 833))
                * p!(mink!(4, 833))
                * gamma!(bis!(4, 832), bis!(4, 831), mink!(4)),
            &q
        ))
        .unwrap();
        let expected = trace!(bis!(4), gamma!(&p), gamma!(&q));
        let result = Arc::new(source.contract(Default::default()).unwrap());
        assert_eq!(result.resolved().unwrap().expression, expected);
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert!(Arc::ptr_eq(
            &result,
            &result.contract(Default::default()).unwrap()
        ));
        assert_eq!(
            source
                .simplify_algebra(&AlgebraSettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .expression,
            expected
        );
        let no_trace = source
            .contract(ContractSettings {
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .resolved()
            .unwrap();
        assert!(
            !no_trace
                .expression
                .contains_symbol(spenso::network::tags::SPENSO_TAG.trace)
        );
        assert!(
            no_trace
                .expression
                .contains_symbol(spenso::network::tags::SPENSO_TAG.chain)
        );
        let filtered = source
            .contract(ContractSettings {
                representations: Some(&[reps.mink4.rep.into()]),
                ..Default::default()
            })
            .unwrap()
            .resolved()
            .unwrap();
        assert!(
            !filtered
                .expression
                .contains_symbol(spenso::network::tags::SPENSO_TAG.trace)
        );
        assert!(
            !filtered
                .expression
                .contains_symbol(spenso::network::tags::SPENSO_TAG.dot)
        );
        assert_eq!(
            source
                .contract(ContractSettings {
                    representations: Some(&[]),
                    ..Default::default()
                })
                .unwrap()
                .resolved()
                .unwrap()
                .expression,
            source.expression
        );
    }

    #[test]
    fn schoonschip_simplifies_metrics() {
        test_initialize();
        let expr = g!(mink!(4, mu), mink!(4, nu)) * p!(mink!(4, nu));

        assert_eq!(
            crate::test_support::contracted_atom(expr.as_view()).unwrap(),
            p!(mink!(4, mu))
        );
    }

    #[test]
    fn metric_only_retains_vector_pairs_and_default_finishes_them() {
        test_initialize();
        let mu = mink!(4, 94401);
        let nu = mink!(4, 94402);
        let spectator = symbolica::parse_lit!((metric_policy_x + metric_policy_y) ^ 30);
        let source = &spectator * g!(&mu, &nu) * p!(&nu) * spenso::q!(&mu);
        let value = SymbolicTensor::infer(source.clone()).unwrap();
        let policy = ContractSettings::default().without_rank_one_tensors();
        let metric_only = Arc::new(value.contract(policy).unwrap());
        let explicit = metric_only.resolved().unwrap();
        assert_eq!(explicit.expression, &spectator * p!(&nu) * spenso::q!(&nu));
        assert!(!metric_only.contraction_complete());
        assert_eq!(explicit.structure, value.structure);
        assert_eq!(
            metric_only.contract(policy).unwrap().resolved().unwrap(),
            explicit
        );

        let full = metric_only.contract(Default::default()).unwrap();
        assert_eq!(
            full.resolved().unwrap().expression,
            spectator * g!(p!(mink!(4)), spenso::q!(mink!(4)))
        );
        assert!(full.contraction_complete());
        assert!(Arc::ptr_eq(
            &full,
            &full.contract(Default::default()).unwrap()
        ));
    }

    #[test]
    fn metric_only_preserves_factored_branches_and_scoped_alias_powers() {
        test_initialize();
        let mu = mink!(4, 94411);
        let nu = mink!(4, 94412);
        let spectator = symbolica::parse_lit!((metric_alias_x + metric_alias_y) ^ 30);
        let body =
            SymbolicTensor::infer((g!(&mu, &nu) * p!(&nu) * spenso::q!(&mu) + Atom::one()).pow(2))
                .unwrap();
        let handle = body.alias_handle().unwrap();
        let root = SymbolicTensor::infer(spectator.clone())
            .unwrap()
            .multiply(&handle)
            .unwrap();
        let value = Arc::new(root.with_aliases([(handle, body)]).unwrap());
        let metric_only = value
            .contract(ContractSettings::default().without_rank_one_tensors())
            .unwrap();
        assert!(!metric_only.contraction_complete());
        assert_eq!(
            metric_only.resolved().unwrap().expression,
            &spectator * (p!(&nu) * spenso::q!(&nu) + Atom::one()).pow(2)
        );
        let full = metric_only.contract(Default::default()).unwrap();
        assert_eq!(
            full.resolved().unwrap().expression,
            spectator * (g!(p!(mink!(4)), spenso::q!(mink!(4))) + Atom::one()).pow(2)
        );
        assert!(full.contraction_complete());
    }

    #[test]
    fn metric_only_relabels_final_vector_port_and_preserves_metadata() {
        use symbolica::atom::FunctionBuilder;
        test_initialize();
        let mu = mink!(4, 94421);
        let nu = mink!(4, 94422);
        let metadata = symbolica::function!(
            symbolica::symbol!("metric_policy_scalar_metadata"; Scalar),
            mink!(4, 94423)
        );
        let vector = spenso::vector_symbol!("metric_policy_metadata");
        let vector_at = |port: &Atom| {
            FunctionBuilder::new(vector)
                .add_arg(&metadata)
                .add_arg(port)
                .finish()
        };
        let invalid = FunctionBuilder::new(vector)
            .add_arg(mink!(4, 94423))
            .add_arg(&nu)
            .finish();
        assert!(SymbolicTensor::infer(g!(&mu, &nu) * invalid).is_err());
        let value = SymbolicTensor::infer(g!(&mu, &nu) * vector_at(&nu)).unwrap();
        assert!(
            SlotContraction::new()
                .contract_factorized(value.expression.as_view(), None, false)
                .is_some()
        );
        let result = value
            .contract(ContractSettings::default().without_rank_one_tensors())
            .unwrap()
            .resolved()
            .unwrap();
        assert_eq!(result.expression, vector_at(&mu));
        assert_eq!(result.structure, value.structure);
    }

    #[test]
    fn contraction_returns_checked_aliases_and_preserves_typed_zero() {
        use spenso::structure::{
            partial::{PartialIndex, PartialStructureExt},
            representation::{LibraryRep, Minkowski, RepName},
            slot::IsAbstractSlot,
        };
        use symbolica::atom::{Atom, AtomCore, FunctionBuilder};

        test_initialize();
        let p = spenso::vector_symbol!("contract_alias_p");
        let q = spenso::vector_symbol!("contract_alias_q");
        let r = spenso::vector_symbol!("contract_alias_r");
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let slot = rep.slot::<AbstractIndex, _>(92101).to_atom();
        let vector = |head| FunctionBuilder::new(head).add_arg(&slot).finish();
        let expression = (vector(p) + vector(q)) * vector(r);
        let source = SymbolicTensor::infer(expression).unwrap();
        let contracted = source.contract(Default::default()).unwrap();
        assert!(!contracted.expression.get_aliases().is_empty());
        assert_eq!(
            contracted.resolved().unwrap().expression.expand(),
            source.expression.expand().schoonschip().expand()
        );
        let zero = SymbolicTensor::checked_parts(
            Atom::Zero,
            PartialStructure::from_logical_slots([
                rep.slot(PartialIndex::Explicit(AbstractIndex::from(92103)))
            ]),
        )
        .unwrap();
        assert_eq!(zero.contract(Default::default()).unwrap().root(), zero);
    }

    #[test]
    fn generated_contraction_weights_are_not_reentered_as_new_templates() {
        test_initialize();
        let left = p!(mink!(4, 92801));
        let right = spenso::q!(mink!(4, 92801));
        let spectator = p!(mink!(4, 92801));
        let body = SymbolicTensor::infer((&left + &right) * &spectator).unwrap();
        let handle = body.alias_handle().unwrap();
        let source = Arc::new(handle.clone().with_aliases([(handle, body)]).unwrap());
        let result = source.contract(Default::default()).unwrap();
        assert!(
            result.expression.get_aliases().len() <= 8,
            "one completed frontier must retain a finite definition registry"
        );
        let expected = g!(p!(mink!(4)), p!(mink!(4))) + g!(spenso::q!(mink!(4)), p!(mink!(4)));
        assert_eq!(
            result.resolved().unwrap().expression.expand(),
            expected.expand()
        );
        assert!(Arc::ptr_eq(
            &result,
            &result.contract(Default::default()).unwrap()
        ));
    }

    #[test]
    fn completed_noop_contraction_does_not_wrap_existing_coefficients() {
        test_initialize();
        let a = mink!(4, 92803);
        let b = mink!(4, 92804);
        let weight = SymbolicTensor::infer(
            Atom::var(symbolica::symbol!("noop_contraction_weight"))
                * g!(p!(mink!(4)), p!(mink!(4))),
        )
        .unwrap();
        let handle = weight.alias_handle().unwrap();
        let metric = SymbolicTensor::infer(g!(&a, &b)).unwrap();
        let root = metric.multiply(&handle).unwrap();
        let value = Arc::new(root.with_aliases([(handle, weight)]).unwrap());
        assert!(!value.contraction_complete());
        let contracted = value.contract(Default::default()).unwrap();
        assert!(contracted.contraction_complete());
        assert_eq!(contracted.root(), value.root());
        assert_eq!(contracted.aliases().unwrap(), value.aliases().unwrap());

        // Domain rewrites legitimately clear completion. Reconstructing the
        // same checked registry must still be stable without that shortcut.
        let unattested = Arc::new(
            contracted
                .root()
                .with_aliases(contracted.aliases().unwrap())
                .unwrap(),
        );
        assert!(!unattested.contraction_complete());
        let repeated = unattested.contract(Default::default()).unwrap();
        assert_eq!(repeated.root(), contracted.root());
        assert_eq!(repeated.aliases().unwrap(), contracted.aliases().unwrap());
        assert!(repeated.contraction_complete());
        assert!(Arc::ptr_eq(
            &repeated,
            &repeated.contract(Default::default()).unwrap()
        ));
    }

    #[test]
    fn contraction_noop_observation_keeps_free_metric_polynomials_factored() {
        use spenso::structure::{representation::RepName, slot::IsAbstractSlot};
        test_initialize();
        let [a, b, c, d] = [92811, 92812, 92813, 92814].map(|index| {
            spenso::structure::representation::LibraryRep::from(
                spenso::structure::representation::Minkowski {},
            )
            .new_rep(4)
            .slot::<AbstractIndex, _>(index)
            .to_atom()
        });
        let spectator = symbolica::parse_lit!((noop_x + noop_y) ^ 8);
        let expression = &spectator * (g!(&a, &b) * g!(&c, &d) + g!(&a, &c) * g!(&b, &d));
        let source = SymbolicTensor::infer(expression.clone()).unwrap();
        let result = source.contract(Default::default()).unwrap();
        assert!(result.contraction_complete());
        assert_eq!(result.root(), source);
        assert!(result.expression.get_aliases().is_empty());
        assert_eq!(result.resolved().unwrap().expression, expression);
    }

    #[test]
    fn contraction_noop_root_still_reduces_scoped_definition_powers() {
        use spenso::structure::{
            abstract_index::AbstractIndex,
            representation::{LibraryRep, Minkowski, RepName},
            slot::IsAbstractSlot,
        };
        use symbolica::atom::FunctionBuilder;
        test_initialize();
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let index = AbstractIndex::Normal(92821).scoped(symbolica::symbol!("noop_scope"));
        let slot = rep.slot::<AbstractIndex, _>(index).to_atom();
        let p = spenso::vector_symbol!("noop_scoped_p");
        let q = spenso::vector_symbol!("noop_scoped_q");
        let vector = |head, arg: &Atom| FunctionBuilder::new(head).add_arg(arg).finish();
        let x = Atom::var(symbolica::symbol!("noop_scoped_weight"));
        for power in [2, 3] {
            let body = SymbolicTensor::infer((vector(p, &slot) * vector(q, &slot) + &x).pow(power))
                .unwrap();
            let handle = body.alias_handle().unwrap();
            let source = Arc::new(handle.clone().with_aliases([(handle, body)]).unwrap());
            let result = source.contract(Default::default()).unwrap();
            let compact = rep.to_symbolic([]);
            let expected = (g!(vector(p, &compact), vector(q, &compact)) + &x).pow(power);
            assert_eq!(
                result.resolved().unwrap().expression.expand(),
                expected.expand()
            );
            assert!(result.contraction_complete());
            assert!(Arc::ptr_eq(
                &result,
                &result.contract(Default::default()).unwrap()
            ));
        }
    }

    #[test]
    fn contraction_checks_callback_induced_rank_loss() {
        use spenso::network::library::symbolic::ETS;
        use symbolica::atom::{Atom, AtomView, FunctionBuilder};

        test_initialize();
        let a = mink!(4, 92201);
        let b = mink!(4, 92203);
        let target = b.clone();
        let head = spenso::tensor_symbol!(
            "contract_alias_callback",
            norm = move |node, output| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let metric = FunctionBuilder::new(ETS.metric)
            .add_arg(&a)
            .add_arg(&b)
            .finish();
        let tensor = FunctionBuilder::new(head).add_arg(&a).finish();
        let source = SymbolicTensor::infer(metric * tensor).unwrap();
        assert!(source.contract(Default::default()).is_err());
        assert!(
            source
                .contract(ContractSettings::default().without_rank_one_tensors())
                .is_err()
        );
    }
}

#[cfg(test)]
mod completion_tests {
    use super::*;
    use spenso::{
        network::library::symbolic::ETS,
        structure::{
            abstract_index::AbstractIndex,
            representation::{LibraryRep, RepName},
            slot::IsAbstractSlot,
        },
    };
    use std::sync::atomic::{AtomicUsize, Ordering};
    use symbolica::atom::FunctionBuilder;

    #[test]
    fn fallback_completes_unevaluated_trace_with_compact_vectors() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 98901);
        let b = spenso::mink!(4, 98902);
        let gamma = |port: Atom| {
            symbolica::function!(
                crate::dirac::AGS.gamma,
                spenso::network::tags::SPENSO_TAG.chain_in,
                spenso::network::tags::SPENSO_TAG.chain_out,
                port
            )
        };
        for rank_one in [false, true] {
            let settings = ContractSettings {
                rank_one,
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            };
            for axial in [false, true] {
                let mut word = vec![
                    gamma(spenso::p!(spenso::mink!(4))),
                    gamma(a.clone()),
                    gamma(b.clone()),
                    gamma(spenso::q!(spenso::mink!(4))),
                ];
                if axial {
                    word.push(symbolica::function!(
                        crate::dirac::AGS.gamma5,
                        spenso::network::tags::SPENSO_TAG.chain_in,
                        spenso::network::tags::SPENSO_TAG.chain_out
                    ));
                }
                let trace = spenso::trace!(crate::bis!(4); word);
                let source = SymbolicTensor::infer(spenso::g!(&a, &b) * trace).unwrap();
                let result = Arc::new(source.contract(settings).unwrap());
                assert_eq!(result.reduction_status(), ReductionStatus::Complete);
                let explicit = result.resolved().unwrap();
                assert!(matches!(explicit.expression.as_view(), AtomView::Fun(fun)
                if fun.get_symbol() == spenso::network::tags::SPENSO_TAG.trace));
                assert!(explicit.expression.contains_symbol(crate::dirac::AGS.gamma));
                assert_eq!(
                    result
                        .contract(settings)
                        .unwrap()
                        .resolved()
                        .unwrap()
                        .expression,
                    explicit.expression,
                );
                assert!(
                    explicit
                        .reduction_observations()
                        .excludes_contraction_sources(settings)
                );
                let zero_budget = explicit
                    .contract(ContractSettings {
                        max_passes: 0,
                        ..settings
                    })
                    .unwrap();
                assert_eq!(zero_budget.reduction_status(), ReductionStatus::Complete);
                assert_eq!(
                    zero_budget.resolved().unwrap().expression,
                    explicit.expression
                );
            }
        }
    }

    #[test]
    fn metric_substitution_into_trace_sum_discharges_repeated_free_ports() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 98931);
        let b = spenso::mink!(4, 98932);
        let c = spenso::mink!(4, 98933);
        let d = spenso::mink!(4, 98934);
        let e = spenso::mink!(4, 98935);
        let gamma = |port: &Atom| {
            symbolica::function!(
                crate::dirac::AGS.gamma,
                spenso::network::tags::SPENSO_TAG.chain_in,
                spenso::network::tags::SPENSO_TAG.chain_out,
                port
            )
        };
        let first = spenso::trace!(crate::bis!(4); [&a, &c, &d, &e].map(gamma));
        let second = spenso::trace!(crate::bis!(4); [&a, &d, &c, &e].map(gamma));
        let source = SymbolicTensor::infer(spenso::g!(&a, &b) * (first + second)).unwrap();
        let settings = ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        };
        let result = Arc::new(source.contract(settings).unwrap());
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        let explicit = result.resolved().unwrap();
        assert!(
            explicit
                .expression
                .contains_symbol(spenso::network::tags::SPENSO_TAG.trace)
        );
        assert!(
            !explicit
                .reduction_observations()
                .candidates
                .repeated_indices
        );
        assert_eq!(
            result
                .contract(settings)
                .unwrap()
                .resolved()
                .unwrap()
                .expression,
            explicit.expression
        );
        assert_eq!(
            explicit.contract(settings).unwrap().reduction_status(),
            ReductionStatus::Complete
        );
    }

    #[test]
    fn fallback_retains_unfinished_sum_connections_with_inexact_coefficients() {
        crate::test_support::test_initialize();
        let port = spenso::mink!(4, 98911);
        let r = spenso::vector_symbol!("completion_inexact::r");
        let s = spenso::vector_symbol!("completion_inexact::s");
        let rounded = Atom::num(symbolica::domains::float::Float::parse("0.1", Some(11)).unwrap());
        // A branch-local rounded coefficient cannot be retained as one opaque
        // scalar spectator of the whole contracted product.
        let expression = (rounded * spenso::p!(&port) + spenso::q!(&port))
            * (symbolica::function!(r, &port) + symbolica::function!(s, &port));
        let source = SymbolicTensor::infer(expression.clone()).unwrap();
        let settings = ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        };
        assert!(
            !source
                .reduction_observations()
                .excludes_contraction_sources(settings)
        );
        assert!(
            SlotContraction::new()
                .contract_factorized(expression.as_view(), None, true)
                .is_none()
        );
        let result = source.contract(settings).unwrap();
        assert_eq!(result.reduction_status(), ReductionStatus::Deferred);
        assert_eq!(result.resolved().unwrap().expression, expression);
    }

    #[test]
    fn dual_fallback_checks_callback_before_certifying_completion() {
        crate::test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let vector = spenso::vector_symbol!(
            "completion_callback_q",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let rep = LibraryRep::new_dual("completion_callback::R")
            .unwrap()
            .new_rep(4);
        let i = rep
            .slot::<AbstractIndex, _>(AbstractIndex::Symbol(
                symbolica::symbol!("completion_callback::i").into(),
            ))
            .to_atom();
        let j = AbstractIndex::Symbol(symbolica::symbol!("completion_callback::j").into());
        let expected = FunctionBuilder::new(vector).add_arg(&i).finish();
        let expression = FunctionBuilder::new(ETS.metric)
            .add_arg(&i)
            .add_arg(rep.dual().slot::<AbstractIndex, _>(j).to_atom())
            .finish()
            * FunctionBuilder::new(vector)
                .add_arg(rep.slot::<AbstractIndex, _>(j).to_atom())
                .finish();
        let source = SymbolicTensor::infer(expression).unwrap();
        calls.store(0, Ordering::Relaxed);
        let result = source.contract(Default::default()).unwrap();
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert!(calls.load(Ordering::Relaxed) > 0);
        assert_eq!(result.resolved().unwrap().expression(), &expected);
    }
}

#[cfg(test)]
mod policy_tests {
    use super::*;
    use spenso::{
        g, mink,
        network::tags::SPENSO_TAG,
        structure::{
            abstract_index::AbstractIndex,
            representation::{Minkowski, RepName},
            slot::IsAbstractSlot,
        },
    };
    use symbolica::atom::{AtomCore, FunctionBuilder};

    #[test]
    fn incompatible_repeated_ports_are_complete_without_spending_budget() {
        crate::test_support::test_initialize();
        let dual = LibraryRep::new_dual("completion_external_ports::R")
            .unwrap()
            .new_rep(4);
        let index =
            AbstractIndex::Symbol(symbolica::symbol!("completion_external_ports::i").into());
        let same_variance = dual.slot::<AbstractIndex, _>(index).to_atom();
        for expression in [
            g!(&same_variance, &same_variance),
            spenso::p!(mink!(4, 98591)) * spenso::q!(mink!(5, 98591)),
        ] {
            let source = SymbolicTensor::infer(expression.clone()).unwrap();
            assert_eq!(source.structure.logical_slots().len(), 2);
            for max_passes in [0, 16] {
                let result = source
                    .contract(ContractSettings {
                        max_passes,
                        ..Default::default()
                    })
                    .unwrap();
                assert_eq!(result.reduction_status(), ReductionStatus::Complete);
                assert_eq!(result.resolved().unwrap().expression, expression);
                assert_eq!(result.root().structure, source.structure);
            }
        }
    }

    #[test]
    fn filters_select_connections_without_discarding_other_ports() {
        crate::test_support::test_initialize();
        let spin = LibraryRep::from(crate::representations::Bispinor {}).new_rep(4);
        let lorentz = LibraryRep::from(Minkowski {});
        let a = mink!(4, 98201);
        let b = mink!(4, 98202);
        let i = spin.slot::<AbstractIndex, _>(98211).to_atom();
        let j = spin.slot::<AbstractIndex, _>(98212).to_atom();
        let head = spenso::tensor_symbol!("contract_filter_mixed_matrix");
        let matrix = |port: &Atom| {
            FunctionBuilder::new(head)
                .add_arg(&i)
                .add_arg(&j)
                .add_arg(port)
                .finish()
        };
        let source = SymbolicTensor::infer(g!(&a, &b) * matrix(&a)).unwrap();
        let empty = source
            .contract(ContractSettings {
                representations: Some(&[]),
                ..Default::default()
            })
            .unwrap();
        assert_eq!(empty.root(), source);
        assert_eq!(
            empty.reduction_status(),
            super::super::simplification::ReductionStatus::Complete
        );
        let result = source
            .contract(ContractSettings {
                representations: Some(&[lorentz]),
                ..Default::default()
            })
            .unwrap();
        assert_eq!(result.resolved().unwrap().expression, matrix(&b));
        assert_eq!(result.structure.root(), source.structure());
        assert!(
            !result
                .expression
                .get_root()
                .contains_symbol(SPENSO_TAG.chain)
        );
    }

    #[test]
    fn collecting_closed_matrix_products_does_not_evaluate_traces() {
        crate::test_support::test_initialize();
        let rep = LibraryRep::new_dual("contract_collect_policy::R").unwrap();
        let r = rep.new_rep(3);
        let i = r.slot::<AbstractIndex, _>(98301).to_atom();
        let j = r.slot::<AbstractIndex, _>(98302).to_atom();
        let di = r.dual().slot::<AbstractIndex, _>(98301).to_atom();
        let dj = r.dual().slot::<AbstractIndex, _>(98302).to_atom();
        let a = spenso::tensor_symbol!("contract_collect_policy::A");
        let b = spenso::tensor_symbol!("contract_collect_policy::B");
        let first = FunctionBuilder::new(a).add_arg(&i).add_arg(&dj).finish();
        let second = FunctionBuilder::new(b).add_arg(&j).add_arg(&di).finish();
        let source = SymbolicTensor::infer(&first * &second).unwrap();
        let open = source
            .contract(ContractSettings {
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .resolved()
            .unwrap();
        use crate::shorthands::chain::Chain;
        assert!(
            open.expression.contains_symbol(SPENSO_TAG.chain),
            "source={} chainify={} primitive={} result={}",
            source.expression,
            source.expression.chainify(rep),
            source.expression.collect_chains(rep, true, false),
            open.expression
        );
        assert!(!open.expression.contains_symbol(SPENSO_TAG.trace));
        assert!(open.expression.contains_symbol(a) && open.expression.contains_symbol(b));
        let traced = source
            .contract(ContractSettings::default())
            .unwrap()
            .resolved()
            .unwrap();
        assert!(traced.expression.contains_symbol(SPENSO_TAG.trace));
        assert!(traced.expression.contains_symbol(a) && traced.expression.contains_symbol(b));
        assert!(!matches!(traced.expression.as_view(), AtomView::Num(_)));

        let body = SymbolicTensor::infer(first).unwrap();
        let handle = body.alias_handle().unwrap();
        let root = SymbolicTensor::infer(&handle.expression * second).unwrap();
        let aliased = Arc::new(root.with_aliases([(handle, body)]).unwrap());
        let collected = aliased
            .contract(ContractSettings {
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .resolved()
            .unwrap();
        assert!(collected.expression.contains_symbol(SPENSO_TAG.chain));
        assert!(!collected.expression.contains_symbol(SPENSO_TAG.trace));
        assert!(collected.expression.contains_symbol(a) && collected.expression.contains_symbol(b));
    }

    #[test]
    fn collects_chains_inside_opaque_sum_scopes() {
        crate::test_support::test_initialize();
        let rep = LibraryRep::new_dual("contract_sum_chain::R")
            .unwrap()
            .new_rep(3);
        let i = rep.slot::<AbstractIndex, _>(98431).to_atom();
        let j = rep.slot::<AbstractIndex, _>(98432).to_atom();
        let dj = rep.dual().slot::<AbstractIndex, _>(98432).to_atom();
        let dk = rep.dual().slot::<AbstractIndex, _>(98433).to_atom();
        let heads = [
            "contract_sum_chain::A",
            "contract_sum_chain::B",
            "contract_sum_chain::C",
            "contract_sum_chain::D",
        ]
        .map(|name| SPENSO_TAG.tensor_symbol(name));
        let matrix = |head, first: &Atom, second: &Atom| {
            FunctionBuilder::new(head)
                .add_arg(first)
                .add_arg(second)
                .finish()
        };
        let spectator = symbolica::parse_lit!((chain_scope_x + chain_scope_y) ^ 40);
        let source = SymbolicTensor::infer(
            &spectator
                * (matrix(heads[0], &i, &dj) * matrix(heads[1], &j, &dk)
                    + matrix(heads[2], &i, &dj) * matrix(heads[3], &j, &dk)),
        )
        .unwrap();
        let settings = ContractSettings {
            collect_traces: false,
            ..Default::default()
        };
        let result = source.contract(settings).unwrap().resolved().unwrap();
        let mut chains = 0;
        result.expression.visitor(&mut |node| {
            chains += usize::from(matches!(node, AtomView::Fun(function)
                if function.get_symbol() == SPENSO_TAG.chain));
            true
        });
        assert_eq!(chains, 2, "{}", result.expression);
        assert_eq!(result.structure(), source.structure());
        assert!(matches!(result.expression.as_view(), AtomView::Mul(product)
            if product.iter().any(|factor| factor == spectator.as_view())));
        assert_eq!(
            result.contract(settings).unwrap().resolved().unwrap(),
            result
        );
    }

    #[test]
    fn metric_alias_relabelling_preserves_disabled_vector_powers() {
        crate::test_support::test_initialize();
        let [a, b, c] = [98421, 98422, 98423].map(|index| {
            LibraryRep::from(Minkowski {})
                .new_rep(4)
                .slot::<AbstractIndex, _>(index)
                .to_atom()
        });
        let tensor = spenso::tensor_symbol!("metric_alias_power_tensor");
        let vector = spenso::vector_symbol!("metric_alias_power_vector");
        let at = |port: &Atom| FunctionBuilder::new(tensor).add_arg(port).finish();
        let square = FunctionBuilder::new(vector).add_arg(c).finish().pow(2);
        let body = SymbolicTensor::infer(at(&a) * &square).unwrap();
        let handle = body.alias_handle().unwrap();
        let root = SymbolicTensor::infer(g!(&a, &b) * &handle.expression).unwrap();
        let source = Arc::new(root.with_aliases([(handle, body)]).unwrap());
        let result = source
            .contract(ContractSettings {
                rank_one: false,
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .resolved()
            .unwrap();
        assert_eq!(result.expression, at(&b) * square);
    }

    #[test]
    fn prerequisite_outer_annihilation_keeps_selected_scalar_child_opaque() {
        crate::test_support::test_initialize();
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let slots = (98501..98509)
            .map(|index| rep.slot::<AbstractIndex, _>(index).to_atom())
            .collect::<Vec<_>>();
        let epsilon = |ports: &[Atom]| {
            FunctionBuilder::new(*crate::epsilon::EPSILON_SYMBOL)
                .add_args(ports)
                .finish()
        };
        let selected_child = (Atom::one() + epsilon(&slots[4..]).pow(2)).pow(20);
        let source =
            SymbolicTensor::infer(g!(&slots[0], &slots[1]) * epsilon(&slots[..4]) * selected_child)
                .unwrap();
        assert!(!source.expression.is_zero());
        super::super::collection::SHALLOW_BUILDS.with(|count| count.set(0));
        super::super::collection::SELECTED_OPENS.with(|count| count.set(0));
        let (result, definitions) = source
            .contract_prerequisites(
                ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                },
                &[],
                |value| value.contains_symbol(*crate::epsilon::EPSILON_SYMBOL),
            )
            .unwrap();
        assert!(result.root.expression.is_zero());
        assert_eq!(result.root.structure(), source.structure());
        assert!(definitions.is_empty());
        super::super::collection::SELECTED_OPENS.with(|count| assert_eq!(count.get(), 0));
        super::super::collection::SHALLOW_BUILDS.with(|count| {
            assert!(
                count.get() <= 2,
                "opened a selected child before outer annihilation: {} views",
                count.get()
            );
        });
    }

    #[test]
    fn prerequisite_contraction_leaves_unrelated_metric_work_exact() {
        crate::test_support::test_initialize();
        let [a, b, c, d] = [98401, 98402, 98403, 98404].map(|index| {
            LibraryRep::from(Minkowski {})
                .new_rep(4)
                .slot::<AbstractIndex, _>(index)
                .to_atom()
        });
        let head = spenso::tensor_symbol!("contract_prerequisite_selected");
        let vector = spenso::vector_symbol!("contract_prerequisite_spectator");
        let selected = |port: &Atom| FunctionBuilder::new(head).add_arg(port).finish();
        let spectator = g!(&c, &d) * FunctionBuilder::new(vector).add_arg(&c).finish();
        let source = SymbolicTensor::infer(g!(&a, &b) * selected(&a) * &spectator).unwrap();
        let result = source
            .contract_prerequisites(
                Default::default(),
                &[],
                |value| matches!(value, AtomView::Fun(function) if function.get_symbol() == head),
            )
            .unwrap()
            .0;
        assert_eq!(result.root.expression, selected(&b) * spectator);
        assert!(result.aliases.is_empty());
    }
}
