use spenso::structure::{
    partial::{PartialStructure, PartialStructureExt},
    representation::{LibraryRep, RepName},
};
use spenso::{
    network::library::symbolic::ETS,
    structure::slot::{DualSlotTo, IsAbstractSlot},
};
use symbolica::atom::{Atom, AtomCore, AtomView};

use super::{
    SymbolicTensor, inference::TensorInferenceError,
    simplification::observation::DomainObservations,
};
use crate::shorthands::schoonschip::{
    DotNormalizer, ReductionStatus, SimplificationCandidates, SlotContraction,
};

#[cfg(test)]
thread_local! {
    pub(crate) static CONTRACT_DOMAIN_CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    pub(crate) static CONTRACT_PARTS_CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
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
    /// Permit local distribution required to contract independent sum alternatives.
    /// False retains such products and powers, while allowing port substitutions
    /// through existing sums and contractions within their branches.
    pub expand: bool,
    /// Optional limit on successful planning rounds per domain. `None` is unlimited.
    pub max_passes: Option<usize>,
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
            expand: true,
            max_passes: None,
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
        self.representations.is_none() && self.metrics && self.rank_one && self.expand
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
        if observed.excludes_internal_connections() {
            return false;
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
            if !(candidates.repeated_indices || candidates.brackets) {
                return false;
            }
            if !settings.expand && settings.max_passes.is_some() {
                // A deliberate budget counts permitted transformations, not
                // products whose remaining contractions require distribution.
                // Reuse admission's graph; probing must not perform rewrites.
                let contractor =
                    SlotContraction::configured(settings.metrics, settings.representations)
                        .observed(observed);
                let contractor = match self.shallow_graph() {
                    Ok(graph) => contractor.planned(
                        self.expression.clone(),
                        graph,
                        self.proofs.leaf_interfaces.get().cloned(),
                    ),
                    Err(_) => contractor,
                };
                return contractor
                    .minimal_work_permitted(self.expression.as_view(), settings)
                    .unwrap_or(true);
            }
            return true;
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

    /// Contract compatible structural connections without evaluating identities
    /// or traces. [`ContractSettings`] controls representations and notation.
    /// Setting `expand = false` retains independent sum alternatives; completion
    /// is relative to that policy. Both policies reuse the outside-in planner and
    /// shallow graph, opening only selected tensor scopes.
    pub fn contract(&self, settings: ContractSettings<'_>) -> Result<Self, TensorInferenceError> {
        self.plan_reduction(super::simplification::ReductionRequest::Contract(settings))
    }

    /// Also reports a notation-only change: nothing was contracted and only
    /// chain and trace collection changed the expression. Collection is
    /// idempotent and adds no contraction source, so such a result is its
    /// own fixed point and needs no confirming round.
    pub(crate) fn contract_domain(
        &self,
        settings: ContractSettings<'_>,
    ) -> Result<(Self, ReductionStatus, bool), TensorInferenceError> {
        #[cfg(test)]
        CONTRACT_DOMAIN_CALLS.with(|count| count.set(count.get() + 1));
        let observed = self.reduction_observations();
        let (contracted, notation_only) =
            self.contract_parts_observed(settings, &observed.candidates)?;
        contracted.root.observe_replacement(self);
        Ok((contracted.root, contracted.status, notation_only))
    }

    /// Contract only metric/vector paths incident to selected identity factors.
    /// Unselected branches retain their exact payload and normalization state.
    pub(crate) fn contract_prerequisites(
        &self,
        settings: ContractSettings<'_>,
        mut selected: impl FnMut(AtomView<'_>) -> bool,
        allowed_vector: impl Fn(LibraryRep) -> bool,
    ) -> Result<crate::shorthands::schoonschip::FactorizedContraction<Self>, TensorInferenceError>
    {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::Contraction,
        );
        self.contract_prerequisite_scopes(settings, &mut selected, &allowed_vector)
    }

    fn contract_prerequisite_scopes(
        &self,
        settings: ContractSettings<'_>,
        selected: &mut dyn FnMut(AtomView<'_>) -> bool,
        allowed_vector: &dyn Fn(LibraryRep) -> bool,
    ) -> Result<crate::shorthands::schoonschip::FactorizedContraction<Self>, TensorInferenceError>
    {
        use super::inference::InterfaceInference;
        use spenso::structure::TensorStructure;
        let observed = self.reduction_observations();
        if observed.excludes_contraction_sources(settings)
            || (!observed.candidates.repeated_indices && !self.structural_work_permitted(settings))
            || self.only_external_metric_sources(settings)
        {
            return Ok(crate::shorthands::schoonschip::FactorizedContraction {
                root: self.clone(),
                status: ReductionStatus::Complete,
                observations: None,
            });
        }
        let mut result = self.contract_prerequisite_boundary(settings, selected, allowed_vector)?;
        if result.root != *self || result.status == ReductionStatus::Capped {
            // Let the shared planner reconsider the changed outer boundary
            // before any selected child is opened.
            return Ok(result);
        }
        // Only an unchanged outer boundary admits work inside a selected
        // arithmetic child. Additive branches and scalar power bases own their
        // independent dummy scopes and remain factorized.
        let needs_sources = |value: AtomView<'_>| {
            observed
                .region(value)
                .is_none_or(|region| !region.excludes_contraction_sources(settings))
        };
        let power_selected = match self.expression.as_view() {
            AtomView::Pow(power) => needs_sources(power.get_base()) && selected(power.get_base()),
            _ => true,
        };
        let mut visit = |value: AtomView<'_>,
                         interface: Option<PartialStructure>|
         -> Result<Atom, TensorInferenceError> {
            // Repeated indices in an algebraic word do not themselves permit
            // structural work. Inspect retained regional facts before inferring
            // a child interface or building another shallow graph.
            if !needs_sources(value) || !selected(value) {
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
            let nested = domain.contract_prerequisite_scopes(settings, selected, allowed_vector)?;
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
        result.root = self.with_identity_result(expression, None)?;
        Ok(result)
    }

    fn contract_prerequisite_boundary(
        &self,
        settings: ContractSettings<'_>,
        selected: &mut dyn FnMut(AtomView<'_>) -> bool,
        allowed_vector: &dyn Fn(LibraryRep) -> bool,
    ) -> Result<crate::shorthands::schoonschip::FactorizedContraction<Self>, TensorInferenceError>
    {
        let mut result = crate::shorthands::schoonschip::FactorizedContraction {
            root: self.clone(),
            status: ReductionStatus::Complete,
            observations: None,
        };
        let contractor = SlotContraction::configured(settings.metrics, settings.representations)
            .expanding(settings.expand)
            .observed(result.root.reduction_observations());
        let contractor = match result.root.shallow_graph() {
            Ok(graph) => contractor.planned(
                result.root.expression.clone(),
                graph,
                result.root.proofs.leaf_interfaces.get().cloned(),
            ),
            Err(_) => contractor,
        };
        let Some((expression, interface, mut spectators)) = contractor.prerequisite_factors(
            result.root.expression.as_view(),
            &mut *selected,
            allowed_vector,
        ) else {
            return Ok(result);
        };
        let domain = if self.normalization_is_intrinsic() {
            Self::from_validated_parts(expression, interface)
        } else {
            Self::checked_parts(expression, interface)?
        };
        domain.observe_replacement(self);
        let contracted = domain.contract_parts(ContractSettings {
            collect_chains: false,
            collect_traces: false,
            order: None,
            ..settings
        })?;
        result.status = result.status.max(contracted.status);
        spectators.push(contracted.root.expression);
        result.root = self.with_identity_result(Atom::mul_many(spectators), None)?;
        Ok(result)
    }

    pub(crate) fn contract_parts(
        &self,
        settings: ContractSettings<'_>,
    ) -> Result<crate::shorthands::schoonschip::FactorizedContraction<Self>, TensorInferenceError>
    {
        #[cfg(test)]
        CONTRACT_PARTS_CALLS.with(|count| count.set(count.get() + 1));
        let observed = self.reduction_observations();
        Ok(self
            .contract_parts_observed(settings, &observed.candidates)?
            .0)
    }

    /// See [`Self::contract_domain`] for the notation-only flag.
    fn contract_parts_observed<const N: usize>(
        &self,
        settings: ContractSettings<'_>,
        observed: &SimplificationCandidates<N>,
    ) -> Result<
        (
            crate::shorthands::schoonschip::FactorizedContraction<Self>,
            bool,
        ),
        TensorInferenceError,
    > {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::Contraction,
        );
        let mut result = self.contract_kernel(settings, observed)?;
        if result
            .root
            .proofs
            .observations
            .get()
            .is_some_and(|observed| observed.excludes_internal_connections())
        {
            return Ok((result, false));
        }
        if settings.rank_one && observed.dots {
            let expression =
                SlotContraction::configured(settings.metrics, settings.representations)
                    .expanding(settings.expand)
                    .contract_inner_products(result.root.expression.as_view());
            result.root = result.root.with_identity_result(expression, None)?;
        }
        let uncontracted =
            result.status == ReductionStatus::Complete && result.root.expression == self.expression;
        if settings.collect_chains || settings.collect_traces {
            result.root.observe_replacement(self);
            let updated = result.root.reduction_observations();
            let mut representations = self
                .reduction_observations()
                .representations
                .iter()
                .copied()
                .filter(|rep| settings.permits(*rep))
                .filter(|rep| updated.has_chain_channel(*rep))
                .map(|rep| rep.base())
                .collect::<Vec<_>>();
            representations.sort_unstable();
            representations.dedup();
            if representations.is_empty() {
                return Ok((result, false));
            }
            let contractor =
                SlotContraction::configured(settings.metrics, settings.representations)
                    .expanding(settings.expand)
                    .observed(updated);
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
            result.root = result.root.with_identity_result(expression, None)?;
        }
        let notation_only = uncontracted && result.root.expression != self.expression;
        Ok((result, notation_only))
    }

    /// A metric whose two ports are both external to the whole domain is
    /// terminal: no factor can absorb it, even when other factors carry
    /// internal dummies. With no vector source either, nothing is left to
    /// contract.
    pub(crate) fn only_external_metric_sources(&self, settings: ContractSettings<'_>) -> bool {
        if !self
            .reduction_observations()
            .excludes_contraction_sources(ContractSettings {
                metrics: false,
                ..settings
            })
        {
            return false;
        }
        let Ok(slots) = self.structure.slots() else {
            return false;
        };
        let external = slots
            .iter()
            .flat_map(|slot| [slot.to_atom(), slot.dual().to_atom()])
            .collect::<std::collections::HashSet<_>>();
        let mut terminal = true;
        self.expression.visitor(&mut |node| {
            if let AtomView::Fun(metric) = node
                && metric.get_symbol() == ETS.metric
            {
                terminal &= metric
                    .iter()
                    .all(|port| external.contains(&port.to_owned()));
                return false;
            }
            terminal
        });
        terminal
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
        let observations = self.reduction_observations();
        let terminal_polynomial = observations.is_terminal_polynomial();
        if observations.excludes_internal_connections()
            || (settings
                .representations
                .is_some_and(<[LibraryRep]>::is_empty)
                || (!settings.metrics && !settings.rank_one))
            || (order.is_none()
                && source == self.expression
                && self.contraction_is_noop_observed(settings, observed))
            || (!settings.expand
                && source == self.expression
                && !source.as_view().needs_normalization()
                && observations.candidates.intrinsic
                && observations.excludes_contraction_sources(settings))
        {
            // Repeated indices in algebraic invariants do not themselves
            // authorize metric/vector work. Reuse the source inventory before
            // constructing its network; unknown and bracket regions still need
            // planning. Chain/trace collection remains in contract_parts.
            return Ok(crate::shorthands::schoonschip::FactorizedContraction {
                root: self.with_identity_result(
                    source,
                    terminal_polynomial.then(|| {
                        // Structural normalization can retain compact metrics;
                        // only the notation owner can discharge their dot work.
                        DomainObservations::terminal_polynomial(
                            &self.structure,
                            observations.needs_dot_normalization(),
                        )
                    }),
                )?,
                status: ReductionStatus::Complete,
                observations: None,
            });
        }
        // Notation normalization can replace whole arithmetic scopes. Keep the
        // payload, regional observations and shallow graph on the same owner.
        let domain = self.with_identity_result(source, None)?;
        let contractor = SlotContraction::configured(settings.metrics, settings.representations)
            .expanding(settings.expand)
            .observed(domain.reduction_observations());
        let contractor = match domain.shallow_graph() {
            Ok(graph) => contractor.planned(
                domain.expression.clone(),
                graph,
                domain.proofs.leaf_interfaces.get().cloned(),
            ),
            Err(_) => contractor,
        };
        if let Some(contracted) =
            contractor.contract_factorized(domain.expression.as_view(), order, settings.rank_one)
        {
            return Ok(crate::shorthands::schoonschip::FactorizedContraction {
                root: self.with_identity_result(contracted.root, contracted.observations)?,
                status: contracted.status,
                observations: None,
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
        let rewritten = contractor.contract_ordered(
            domain.expression.as_view(),
            settings.rank_one,
            Some(&domain.reduction_observations().candidates),
        );
        let root = domain.with_identity_result(rewritten, None)?;
        if root.expression.is_zero() {
            // The ordered rewrite already annihilated the whole domain. Keep
            // its typed boundary without opening any now-irrelevant scopes.
            return Ok(crate::shorthands::schoonschip::FactorizedContraction {
                root,
                status: ReductionStatus::Complete,
                observations: None,
            });
        }
        if root.expression != domain.expression {
            let contractor = contractor.observed(root.reduction_observations());
            let contractor = match root.shallow_graph() {
                Ok(graph) => contractor.planned(
                    root.expression.clone(),
                    graph,
                    root.proofs.leaf_interfaces.get().cloned(),
                ),
                Err(_) => contractor,
            };
            if let Some(contracted) =
                contractor.contract_factorized(root.expression.as_view(), None, settings.rank_one)
            {
                return Ok(crate::shorthands::schoonschip::FactorizedContraction {
                    root: root.with_identity_result(contracted.root, contracted.observations)?,
                    status: contracted.status,
                    observations: None,
                });
            }
        }
        // The ordered fallback supports spaces declined by the coefficient
        // collector. Repeated gamma ports or compact slash vectors are not
        // structural source endpoints. Their existing regional observations can
        // certify completion without a new walk; explicit metric/vector sources
        // inside unsupported sums remain pending. An incomplete collector frontier
        // returned above retains its exact remaining factors.
        let complete = !root.structural_work_permitted(settings)
            || root
                .reduction_observations()
                .excludes_contraction_sources(settings)
            || (root.expression != self.expression && root.contraction_is_noop())
            || root.only_external_metric_sources(settings);
        Ok(crate::shorthands::schoonschip::FactorizedContraction {
            root,
            observations: None,
            status: if complete {
                ReductionStatus::Complete
            } else {
                ReductionStatus::Deferred
            },
        })
    }
}

#[cfg(test)]
mod factorized_tests;
#[cfg(test)]
mod minimal_tests;
#[cfg(test)]
mod network_tests;

#[cfg(test)]
pub mod test {
    use super::super::AbstractIndex;
    use super::*;
    use crate::shorthands::schoonschip::Schoonschip;
    use crate::test_support::test_initialize;
    use spenso::{g, mink, p};
    use std::sync::Arc;

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
        let result = source.contract(Default::default()).unwrap();
        assert_eq!(result.expression, expected);
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert!(Arc::ptr_eq(
            result.reduction_observations(),
            result
                .contract(Default::default())
                .unwrap()
                .reduction_observations()
        ));
        assert_eq!(
            source
                .simplify_algebra(&AlgebraSettings::default())
                .unwrap()
                .expression,
            expected
        );
        let no_trace = source
            .contract(ContractSettings {
                collect_traces: false,
                ..Default::default()
            })
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
        let metric_only = value.contract(policy).unwrap();
        let explicit = metric_only.clone();
        assert_eq!(explicit.expression, &spectator * p!(&nu) * spenso::q!(&nu));
        assert!(!metric_only.contraction_complete());
        assert_eq!(explicit.structure, value.structure);
        assert_eq!(metric_only.contract(policy).unwrap(), explicit);

        let full = metric_only.contract(Default::default()).unwrap();
        assert_eq!(
            full.expression,
            spectator * g!(p!(mink!(4)), spenso::q!(mink!(4)))
        );
        assert!(full.contraction_complete());
        assert!(Arc::ptr_eq(
            full.reduction_observations(),
            full.contract(Default::default())
                .unwrap()
                .reduction_observations()
        ));
    }

    #[test]
    fn metric_only_preserves_factored_branches_and_scoped_powers() {
        test_initialize();
        let mu = mink!(4, 94411);
        let nu = mink!(4, 94412);
        let spectator = symbolica::parse_lit!((metric_scope_x + metric_scope_y) ^ 30);
        let body =
            SymbolicTensor::infer((g!(&mu, &nu) * p!(&nu) * spenso::q!(&mu) + Atom::one()).pow(2))
                .unwrap();
        let value = SymbolicTensor::infer(spectator.clone())
            .unwrap()
            .multiply(&body)
            .unwrap();
        let metric_only = value
            .contract(ContractSettings::default().without_rank_one_tensors())
            .unwrap();
        assert!(!metric_only.contraction_complete());
        assert_eq!(
            metric_only.expression,
            &spectator * (p!(&nu) * spenso::q!(&nu) + Atom::one()).pow(2)
        );
        let full = metric_only.contract(Default::default()).unwrap();
        assert_eq!(
            full.expression,
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
            .unwrap();
        assert_eq!(result.expression, vector_at(&mu));
        assert_eq!(result.structure, value.structure);
    }

    #[test]
    fn contraction_returns_a_plain_tensor_and_preserves_typed_zero() {
        use spenso::structure::{
            partial::{PartialIndex, PartialStructureExt},
            representation::{LibraryRep, Minkowski, RepName},
            slot::IsAbstractSlot,
        };
        use symbolica::atom::{Atom, AtomCore, FunctionBuilder};

        test_initialize();
        let p = spenso::vector_symbol!("contract_plain_p");
        let q = spenso::vector_symbol!("contract_plain_q");
        let r = spenso::vector_symbol!("contract_plain_r");
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let slot = rep.slot::<AbstractIndex, _>(92101).to_atom();
        let vector = |head| FunctionBuilder::new(head).add_arg(&slot).finish();
        let expression = (vector(p) + vector(q)) * vector(r);
        let source = SymbolicTensor::infer(expression).unwrap();
        let contracted = source.contract(Default::default()).unwrap();
        assert_eq!(contracted.structure, source.structure);
        assert_eq!(
            contracted.expression.expand(),
            source.expression.expand().schoonschip().expand()
        );
        let zero = SymbolicTensor::checked_parts(
            Atom::Zero,
            PartialStructure::from_logical_slots([
                rep.slot(PartialIndex::Explicit(AbstractIndex::from(92103)))
            ]),
        )
        .unwrap();
        assert_eq!(zero.contract(Default::default()).unwrap(), zero);
    }

    #[test]
    fn generated_contraction_weights_are_stable_on_reruns() {
        test_initialize();
        let left = p!(mink!(4, 92801));
        let right = spenso::q!(mink!(4, 92801));
        let spectator = p!(mink!(4, 92801));
        let source = SymbolicTensor::infer((&left + &right) * &spectator).unwrap();
        let result = source.contract(Default::default()).unwrap();
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        let expected = g!(p!(mink!(4)), p!(mink!(4))) + g!(spenso::q!(mink!(4)), p!(mink!(4)));
        assert_eq!(result.expression.expand(), expected.expand());
        assert!(Arc::ptr_eq(
            result.reduction_observations(),
            result
                .contract(Default::default())
                .unwrap()
                .reduction_observations()
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
        let metric = SymbolicTensor::infer(g!(&a, &b)).unwrap();
        let value = metric.multiply(&weight).unwrap();
        assert!(!value.contraction_complete());
        let contracted = value.contract(Default::default()).unwrap();
        assert!(contracted.contraction_complete());
        assert_eq!(contracted, value);

        // Re-admission clears completion. The same expression must still be
        // stable without that shortcut.
        let unattested = SymbolicTensor::checked_parts(
            contracted.expression.clone(),
            contracted.structure.clone(),
        )
        .unwrap();
        assert!(!unattested.contraction_complete());
        let repeated = unattested.contract(Default::default()).unwrap();
        assert_eq!(repeated, contracted);
        assert!(repeated.contraction_complete());
        assert!(Arc::ptr_eq(
            repeated.reduction_observations(),
            repeated
                .contract(Default::default())
                .unwrap()
                .reduction_observations()
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
        assert_eq!(result, source);
        assert_eq!(result.expression, expression);
    }

    #[test]
    fn contraction_reduces_scoped_scalar_powers() {
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
            let source =
                SymbolicTensor::infer((vector(p, &slot) * vector(q, &slot) + &x).pow(power))
                    .unwrap();
            let result = source.contract(Default::default()).unwrap();
            let compact = rep.to_symbolic([]);
            let expected = (g!(vector(p, &compact), vector(q, &compact)) + &x).pow(power);
            assert_eq!(result.expression.expand(), expected.expand());
            assert!(result.contraction_complete());
            assert!(Arc::ptr_eq(
                result.reduction_observations(),
                result
                    .contract(Default::default())
                    .unwrap()
                    .reduction_observations()
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
            "contract_rank_callback",
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
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
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
                let result = source.contract(settings).unwrap();
                assert_eq!(result.reduction_status(), ReductionStatus::Complete);
                let explicit = result.clone();
                assert!(matches!(explicit.expression.as_view(), AtomView::Fun(fun)
                if fun.get_symbol() == spenso::network::tags::SPENSO_TAG.trace));
                assert!(explicit.expression.contains_symbol(crate::dirac::AGS.gamma));
                assert_eq!(
                    result.contract(settings).unwrap().expression,
                    explicit.expression,
                );
                assert!(
                    explicit
                        .reduction_observations()
                        .excludes_contraction_sources(settings)
                );
                let zero_budget = explicit
                    .contract(ContractSettings {
                        max_passes: Some(0),
                        ..settings
                    })
                    .unwrap();
                assert_eq!(zero_budget.reduction_status(), ReductionStatus::Complete);
                assert_eq!(zero_budget.expression, explicit.expression);
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
        let result = source.contract(settings).unwrap();
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        let explicit = result.clone();
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
            result.contract(settings).unwrap().expression,
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
        assert_eq!(result.expression, expression);
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
        assert_eq!(result.expression(), &expected);
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
            for max_passes in [Some(0), Some(16), None] {
                let result = source
                    .contract(ContractSettings {
                        max_passes,
                        ..Default::default()
                    })
                    .unwrap();
                assert_eq!(result.reduction_status(), ReductionStatus::Complete);
                assert_eq!(result.expression, expression);
                assert_eq!(result.structure, source.structure);
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
        assert_eq!(empty, source);
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
        assert_eq!(result.expression, matrix(&b));
        assert_eq!(result.structure(), source.structure());
        assert!(!result.expression.contains_symbol(SPENSO_TAG.chain));
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
            .unwrap();
        use crate::shorthands::chain::Chain;
        assert!(
            open.expression.contains_symbol(SPENSO_TAG.chain),
            "source={} chainify={} primitive={} result={}",
            source.expression,
            source.expression.chainify(rep),
            source.expression.collect_chains(rep, true, false, false),
            open.expression
        );
        assert!(!open.expression.contains_symbol(SPENSO_TAG.trace));
        assert!(open.expression.contains_symbol(a) && open.expression.contains_symbol(b));
        let traced = source.contract(ContractSettings::default()).unwrap();
        assert!(traced.expression.contains_symbol(SPENSO_TAG.trace));
        assert!(traced.expression.contains_symbol(a) && traced.expression.contains_symbol(b));
        assert!(!matches!(traced.expression.as_view(), AtomView::Num(_)));
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
        let result = source.contract(settings).unwrap();
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
        assert_eq!(result.contract(settings).unwrap(), result);
    }

    #[test]
    fn metric_relabelling_preserves_disabled_vector_powers() {
        crate::test_support::test_initialize();
        let [a, b, c] = [98421, 98422, 98423].map(|index| {
            LibraryRep::from(Minkowski {})
                .new_rep(4)
                .slot::<AbstractIndex, _>(index)
                .to_atom()
        });
        let tensor = spenso::tensor_symbol!("metric_scoped_power_tensor");
        let vector = spenso::vector_symbol!("metric_scoped_power_vector");
        let at = |port: &Atom| FunctionBuilder::new(tensor).add_arg(port).finish();
        let square = FunctionBuilder::new(vector).add_arg(c).finish().pow(2);
        let source = SymbolicTensor::infer(g!(&a, &b) * at(&a) * &square).unwrap();
        let result = source
            .contract(ContractSettings {
                rank_one: false,
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
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
        let result = source
            .contract_prerequisites(
                ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                },
                |value| value.contains_symbol(*crate::epsilon::EPSILON_SYMBOL),
                |_| true,
            )
            .unwrap();
        assert!(result.root.expression.is_zero());
        assert_eq!(result.root.structure(), source.structure());
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
    fn prerequisites_reuse_source_free_factored_color_regions() {
        use super::super::{collection::SHALLOW_BUILDS, inference::tests::SCOPE_VALIDATIONS};
        use crate::{color_f, color_t, test_support::test_initialize};
        use std::sync::Arc;
        let reps = test_initialize();
        let ports = [98521, 98522, 98523].map(|index| reps.coad_da.to_symbolic([Atom::num(index)]));
        let rep = reps.cof_nc.to_symbolic([] as [Atom; 0]);
        let word = spenso::trace!(&rep; ports.iter().map(|port| color_t!(port)));
        let symmetric = spenso::trace_sym!(&rep; ports.iter().map(|port| color_t!(port)));
        let closed = color_f!(&ports[0], &ports[1], &ports[2]) * (word + symmetric);
        let source = SymbolicTensor::infer((Atom::one() + closed).pow(3)).unwrap();
        let settings = ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        };
        let selected = |value: AtomView<'_>| value.contains_symbol(crate::color::CS.t);
        let observed = source.reduction_observations();
        assert!(observed.candidates.repeated_indices);
        SHALLOW_BUILDS.with(|count| count.set(0));
        SCOPE_VALIDATIONS.with(|count| count.set(0));
        let result = source
            .contract_prerequisites(settings, selected, |_| true)
            .unwrap();
        assert_eq!(result.root, source);
        assert_eq!(result.status, ReductionStatus::Complete);
        assert!(Arc::ptr_eq(result.root.reduction_observations(), observed));
        SHALLOW_BUILDS.with(|count| assert_eq!(count.get(), 0));
        SCOPE_VALIDATIONS.with(|count| assert_eq!(count.get(), 0));

        // A source elsewhere in the domain must not reopen the colour scope.
        // Warm the actual boundary graph so the counter isolates child intake.
        let mu = mink!(4, 98524);
        let nu = mink!(4, 98525);
        let mixed =
            SymbolicTensor::infer(&source.expression * g!(&mu, &nu) * spenso::p!(&nu)).unwrap();
        mixed.shallow_graph().unwrap();
        mixed.reduction_observations();
        SCOPE_VALIDATIONS.with(|count| count.set(0));
        let result = mixed
            .contract_prerequisites(settings, selected, |_| true)
            .unwrap();
        assert_eq!(result.root, mixed);
        assert_eq!(result.status, ReductionStatus::Complete);
        SCOPE_VALIDATIONS.with(|count| assert_eq!(count.get(), 0));
    }

    #[test]
    fn prerequisites_contract_metrics_inside_selected_scalar_powers() {
        use crate::{color_t, test_support::test_initialize};
        let reps = test_initialize();
        let [a, b] = [98531, 98532].map(|index| reps.coad_da.to_symbolic([Atom::num(index)]));
        let rep = reps.cof_nc.to_symbolic([] as [Atom; 0]);
        let trace = spenso::trace!(&rep, color_t!(&a), color_t!(&b));
        let source = SymbolicTensor::infer((Atom::one() + g!(&a, &b) * trace).pow(2)).unwrap();
        let settings = ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        };
        let selected = |value: AtomView<'_>| value.contains_symbol(crate::color::CS.t);
        let filtered = source
            .contract_prerequisites(
                ContractSettings {
                    representations: Some(&[]),
                    ..settings
                },
                selected,
                |_| true,
            )
            .unwrap();
        assert_eq!(filtered.root, source);
        let result = source
            .contract_prerequisites(settings, selected, |_| true)
            .unwrap();
        let expected = (Atom::one() + spenso::trace!(&rep, color_t!(&b), color_t!(&b))).pow(2);
        assert_eq!(result.root.expression, expected);
        assert_eq!(result.root.structure(), source.structure());
        assert_eq!(result.status, ReductionStatus::Complete);
        let rerun = result
            .root
            .contract_prerequisites(settings, selected, |_| true)
            .unwrap();
        assert_eq!(rerun.root, result.root);
        assert_eq!(rerun.status, ReductionStatus::Complete);
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
                |value| matches!(value, AtomView::Fun(function) if function.get_symbol() == head),
                |_| true,
            )
            .unwrap();
        assert_eq!(result.root.expression, selected(&b) * spectator);
    }
}
