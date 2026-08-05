use spenso::{
    contraction::Contract,
    network::{
        ExecutionResult, Sequential, SmallestDegree,
        library::{DummyLibrary, function_lib::Wrap},
    },
    structure::slot::{AbsInd, DummyAind, ParseableAind},
};
use symbolica::atom::Atom;

use super::{
    CanonicalPolicyNet, CanonicalizationError,
    projection::{DEFAULT_GRAPH_BUDGET, ProblemIdentity, project},
    reconstruct::{
        DummyAllocator, RetryReason, SignedProblemIdentity, prepare_reconstruction, reconstruct,
    },
};
use crate::tensor::SymbolicTensor;

// See the Phase 6 budget table in the signed-canonicalization architecture plan.
pub(super) const DEFAULT_ITERATION_LIMIT: usize = 8;

#[cfg(test)]
#[derive(Clone, Copy, PartialEq, Eq)]
enum InjectedFailureStage {
    Execution,
    ResultExtraction,
}

#[cfg(test)]
thread_local! {
    static EXECUTION_CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static TEMPORARY_SCOPE_EXECUTIONS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static TEMPORARY_SCOPE_ACTIVE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static INJECTED_EXECUTION_FAILURE: std::cell::Cell<Option<(bool, InjectedFailureStage)>> = const { std::cell::Cell::new(None) };
}

#[cfg(test)]
fn reset_execution_calls() {
    EXECUTION_CALLS.with(|calls| calls.set(0));
    TEMPORARY_SCOPE_EXECUTIONS.with(|calls| calls.set(0));
    TEMPORARY_SCOPE_ACTIVE.with(|active| active.set(false));
    INJECTED_EXECUTION_FAILURE.with(|failure| failure.set(None));
}

#[cfg(test)]
fn execution_calls() -> (usize, usize) {
    let total = EXECUTION_CALLS.with(std::cell::Cell::get);
    let temporary = TEMPORARY_SCOPE_EXECUTIONS.with(std::cell::Cell::get);
    (total - temporary, temporary)
}

#[cfg(test)]
pub(super) fn record_temporary_scope_execution() {
    TEMPORARY_SCOPE_EXECUTIONS.with(|calls| calls.set(calls.get() + 1));
    TEMPORARY_SCOPE_ACTIVE.with(|active| active.set(true));
}

#[cfg(test)]
fn inject_execution_failure(temporary: bool) {
    INJECTED_EXECUTION_FAILURE
        .with(|failure| failure.set(Some((temporary, InjectedFailureStage::Execution))));
}

#[cfg(test)]
fn inject_result_extraction_failure(temporary: bool) {
    INJECTED_EXECUTION_FAILURE
        .with(|failure| failure.set(Some((temporary, InjectedFailureStage::ResultExtraction))));
}

#[cfg(test)]
fn take_injected_failure(temporary: bool, stage: InjectedFailureStage) -> bool {
    INJECTED_EXECUTION_FAILURE.with(|failure| {
        if failure.get() == Some((temporary, stage)) {
            failure.set(None);
            true
        } else {
            false
        }
    })
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CanonicalProblemIdentity {
    graph: ProblemIdentity,
    signed: SignedProblemIdentity,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum IdentityObservation {
    New,
    Consecutive,
    Cycle {
        first_iteration: usize,
        repeated_iteration: usize,
        cycle_length: usize,
        retry_reasons: Vec<RetryReason>,
    },
}

impl IdentityObservation {
    fn into_cycle_error(self) -> CanonicalizationError {
        let Self::Cycle {
            first_iteration,
            repeated_iteration,
            cycle_length,
            retry_reasons,
        } = self
        else {
            unreachable!("only a cycle observation can produce a cycle error")
        };
        CanonicalizationError::ConvergenceCycle {
            first_iteration,
            repeated_iteration,
            cycle_length,
            retry_reasons: retry_reasons
                .into_iter()
                .map(|reason| reason.to_string())
                .collect(),
        }
    }
}

struct IdentityHistory<I> {
    identities: Vec<I>,
    retry_reasons: Vec<RetryReason>,
}

impl<I: Eq> IdentityHistory<I> {
    fn new() -> Self {
        Self {
            identities: Vec::new(),
            retry_reasons: Vec::new(),
        }
    }

    fn observe(&mut self, identity: I) -> IdentityObservation {
        let repeated_iteration = self.identities.len();
        if self.identities.last() == Some(&identity) {
            return IdentityObservation::Consecutive;
        }
        if let Some(first_iteration) = self
            .identities
            .iter()
            .position(|candidate| candidate == &identity)
        {
            return IdentityObservation::Cycle {
                first_iteration,
                repeated_iteration,
                cycle_length: repeated_iteration - first_iteration,
                retry_reasons: self.retry_reasons[first_iteration..repeated_iteration].to_vec(),
            };
        }
        self.identities.push(identity);
        IdentityObservation::New
    }

    fn record_retry(&mut self, reason: RetryReason) {
        debug_assert_eq!(self.retry_reasons.len() + 1, self.identities.len());
        self.retry_reasons.push(reason);
    }
}

impl<Aind> CanonicalPolicyNet<Aind>
where
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
{
    pub(crate) fn canonize<F>(self, mut new_dummy: F) -> Result<Self, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
    {
        self.canonize_with(
            &mut new_dummy,
            DEFAULT_ITERATION_LIMIT,
            DEFAULT_GRAPH_BUDGET,
        )
    }

    fn canonize_with<F>(
        mut self,
        new_dummy: &mut F,
        iteration_limit: usize,
        graph_budget: super::projection::GraphBudget,
    ) -> Result<Self, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
    {
        let mut identities = IdentityHistory::<CanonicalProblemIdentity>::new();
        let mut dummy_allocator = DummyAllocator::<Aind>::new();
        let mut last_reason = "initial canonical projection".to_owned();

        for _ in 0..iteration_limit {
            let projection = project(&self, graph_budget)?;
            let prepared = prepare_reconstruction(&projection)?;
            let identity = CanonicalProblemIdentity {
                graph: projection.identity.clone(),
                signed: prepared.identity.clone(),
            };
            match identities.observe(identity) {
                IdentityObservation::Consecutive => return Ok(self),
                cycle @ IdentityObservation::Cycle { .. } => {
                    return Err(cycle.into_cycle_error());
                }
                IdentityObservation::New => {}
            }

            let reconstruction =
                reconstruct(&projection, prepared, &mut dummy_allocator, new_dummy)?;
            let atom = execute_atom(reconstruction.network.clone())?;
            let terminal = atom.as_view().is_zero() || atom.as_view().is_one();
            // The parser and projection are deterministic, so exact equality
            // after mandatory execution is a complete extensional fixed point.
            let execution_unchanged = atom == *self.normalized_atom();
            let reparsed = CanonicalPolicyNet::parse(atom)?;
            if terminal || execution_unchanged {
                return Ok(reparsed);
            }
            let mut retry_reason = reconstruction.retry_reason;
            if retry_reason == Some(RetryReason::IncompleteStabilityCertificate) {
                // This coarse profile only makes retry diagnostics more specific.
                // Fixed-point correctness is still decided by exact executed atoms.
                let operation_profile = |network: &crate::tensor::SymbolicNet<Aind>| {
                    let mut profile = [0_usize; 5];
                    for (_, _, node) in network.graph.graph.iter_nodes() {
                        let spenso::network::graph::NetworkNode::Op(operation) = node else {
                            continue;
                        };
                        let kind = match operation {
                            spenso::network::graph::NetworkOp::Product => 0,
                            spenso::network::graph::NetworkOp::Sum => 1,
                            spenso::network::graph::NetworkOp::Neg => 2,
                            spenso::network::graph::NetworkOp::Function(_) => 3,
                            spenso::network::graph::NetworkOp::Power(_) => 4,
                        };
                        profile[kind] += 1;
                    }
                    profile
                };
                if operation_profile(&reconstruction.network)
                    != operation_profile(reparsed.network())
                {
                    retry_reason = Some(RetryReason::ExecutionTopologyChange);
                }
            }
            let Some(retry_reason) = retry_reason else {
                return Ok(reparsed);
            };
            last_reason = retry_reason.to_string();
            identities.record_retry(retry_reason);
            self = reparsed;
        }

        Err(CanonicalizationError::IterationLimit {
            limit: iteration_limit,
            last_reason,
        })
    }
}

pub(super) fn execute_atom<Aind>(
    mut network: crate::tensor::SymbolicNet<Aind>,
) -> Result<Atom, CanonicalizationError>
where
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
    SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
{
    #[cfg(test)]
    let temporary = {
        EXECUTION_CALLS.with(|calls| calls.set(calls.get() + 1));
        TEMPORARY_SCOPE_ACTIVE.with(|active| active.replace(false))
    };
    #[cfg(test)]
    if take_injected_failure(temporary, InjectedFailureStage::Execution) {
        return Err(CanonicalizationError::Execution {
            scope: if temporary {
                "temporary canonical scope".to_owned()
            } else {
                "complete network".to_owned()
            },
            error: "injected execution failure".to_owned(),
        });
    }
    let library = DummyLibrary::<SymbolicTensor<Aind>>::new();
    network
        .execute::<Sequential, SmallestDegree<()>, _, _, _>(&library, &Wrap {})
        .map_err(|error| CanonicalizationError::Execution {
            scope: "complete network".to_owned(),
            error: error.to_string(),
        })?;
    #[cfg(test)]
    if take_injected_failure(temporary, InjectedFailureStage::ResultExtraction) {
        return Err(CanonicalizationError::Execution {
            scope: if temporary {
                "temporary canonical scope result".to_owned()
            } else {
                "complete network result".to_owned()
            },
            error: "injected result extraction failure".to_owned(),
        });
    }
    match network
        .result_tensor(&library)
        .map_err(|error| CanonicalizationError::Execution {
            scope: "complete network result".to_owned(),
            error: error.to_string(),
        })? {
        ExecutionResult::One => Ok(Atom::num(1)),
        ExecutionResult::Zero => Ok(Atom::Zero),
        ExecutionResult::Val(tensor) => Ok(tensor.expression.clone()),
    }
}

#[cfg(test)]
mod tests {
    use std::{
        collections::HashSet,
        hash::{Hash, Hasher},
    };

    use spenso::{
        network::{
            graph::{NetworkLeaf, NetworkNode},
            tags::SPENSO_TAG,
        },
        slot,
        structure::{
            abstract_index::AbstractIndex, representation::LibrarySlot, slot::IsAbstractSlot,
        },
        tensor_symbol,
    };
    use symbolica::{
        atom::{Atom, AtomCore, AtomView, NamespacedSymbol, SymbolBuilder},
        function,
    };

    use super::*;
    use crate::{
        CookMode, CookSettings, IndexTooling, tensor::SymbolicNetExt, test_support::test_initialize,
    };

    #[test]
    fn consecutive_identity_is_a_fixed_point() {
        let mut history = IdentityHistory::new();

        assert_eq!(history.observe('a'), IdentityObservation::New);
        history.record_retry(RetryReason::IncompleteStabilityCertificate);
        assert_eq!(history.observe('a'), IdentityObservation::Consecutive);
    }

    #[test]
    fn identity_cycles_report_their_exact_span() {
        let mut two_state = IdentityHistory::new();
        assert_eq!(two_state.observe('a'), IdentityObservation::New);
        two_state.record_retry(RetryReason::SignedPayloadEdit);
        assert_eq!(two_state.observe('b'), IdentityObservation::New);
        two_state.record_retry(RetryReason::VisibleClassMergeSplit);
        let two_state_cycle = two_state.observe('a');
        assert_eq!(
            two_state_cycle,
            IdentityObservation::Cycle {
                first_iteration: 0,
                repeated_iteration: 2,
                cycle_length: 2,
                retry_reasons: vec![
                    RetryReason::SignedPayloadEdit,
                    RetryReason::VisibleClassMergeSplit,
                ],
            }
        );
        assert_eq!(
            two_state_cycle.into_cycle_error(),
            CanonicalizationError::ConvergenceCycle {
                first_iteration: 0,
                repeated_iteration: 2,
                cycle_length: 2,
                retry_reasons: vec![
                    "signed payload edit".to_owned(),
                    "visible-class merge/split".to_owned(),
                ],
            }
        );

        let mut three_state = IdentityHistory::new();
        assert_eq!(three_state.observe('a'), IdentityObservation::New);
        three_state.record_retry(RetryReason::ZeroSubstitution);
        assert_eq!(three_state.observe('b'), IdentityObservation::New);
        three_state.record_retry(RetryReason::ResidualPowerPhase);
        assert_eq!(three_state.observe('c'), IdentityObservation::New);
        three_state.record_retry(RetryReason::ExecutionTopologyChange);
        let three_state_cycle = three_state.observe('a');
        assert_eq!(
            three_state_cycle,
            IdentityObservation::Cycle {
                first_iteration: 0,
                repeated_iteration: 3,
                cycle_length: 3,
                retry_reasons: vec![
                    RetryReason::ZeroSubstitution,
                    RetryReason::ResidualPowerPhase,
                    RetryReason::ExecutionTopologyChange,
                ],
            }
        );
        assert_eq!(
            three_state_cycle.into_cycle_error(),
            CanonicalizationError::ConvergenceCycle {
                first_iteration: 0,
                repeated_iteration: 3,
                cycle_length: 3,
                retry_reasons: vec![
                    "zero substitution".to_owned(),
                    "residual Power phase".to_owned(),
                    "execution-topology change".to_owned(),
                ],
            }
        );
    }

    #[test]
    fn identity_cycles_return_the_typed_error_with_retry_reasons() {
        let error = IdentityObservation::Cycle {
            first_iteration: 1,
            repeated_iteration: 4,
            cycle_length: 3,
            retry_reasons: vec![
                RetryReason::SignedPayloadEdit,
                RetryReason::VisibleClassMergeSplit,
                RetryReason::ExecutionTopologyChange,
            ],
        }
        .into_cycle_error();

        assert_eq!(
            error,
            CanonicalizationError::ConvergenceCycle {
                first_iteration: 1,
                repeated_iteration: 4,
                cycle_length: 3,
                retry_reasons: vec![
                    "signed payload edit".to_owned(),
                    "visible-class merge/split".to_owned(),
                    "execution-topology change".to_owned(),
                ],
            }
        );
    }

    #[test]
    fn iteration_limit_reports_the_last_retry_reason() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, canonical_driver_limit_first);
        let second = slot!(rep, canonical_driver_limit_second);
        let tensor = tensor_symbol!(canonical_driver_limit_tensor; Symmetric);
        let expression = function!(tensor, first.to_atom(), second.to_atom()).pow(Atom::num(2));
        let mut new_dummy = AbstractIndex::Dummy;

        let result = CanonicalPolicyNet::<AbstractIndex>::parse(expression)
            .unwrap()
            .canonize_with(
                &mut new_dummy,
                1,
                super::super::projection::DEFAULT_GRAPH_BUDGET,
            );
        let Err(error) = result else {
            panic!("one uncertified reconstruction must exhaust a one-iteration budget");
        };

        assert_eq!(
            error,
            CanonicalizationError::IterationLimit {
                limit: 1,
                last_reason: "incomplete stability certificate".to_owned(),
            }
        );
    }

    #[test]
    fn distinct_synthetic_identities_reach_the_explicit_budget_without_false_repetition() {
        let limit = 3;
        let states = [
            ('a', RetryReason::SignedPayloadEdit),
            ('b', RetryReason::ZeroSubstitution),
            ('c', RetryReason::ExecutionTopologyChange),
        ];
        let mut history = IdentityHistory::new();
        let mut last_reason = "initial canonical projection".to_owned();

        for (identity, reason) in states {
            assert_eq!(history.observe(identity), IdentityObservation::New);
            last_reason = reason.to_string();
            history.record_retry(reason);
        }

        assert_eq!(history.identities, vec!['a', 'b', 'c']);
        assert_eq!(history.retry_reasons.len(), limit);
        assert_eq!(
            CanonicalizationError::IterationLimit { limit, last_reason },
            CanonicalizationError::IterationLimit {
                limit: 3,
                last_reason: "execution-topology change".to_owned(),
            }
        );
    }

    #[test]
    fn colliding_candidate_hashes_do_not_replace_full_problem_identity_equality() {
        #[derive(Clone, Debug, PartialEq, Eq)]
        struct CollidingIdentity(CanonicalProblemIdentity);

        impl Hash for CollidingIdentity {
            fn hash<H: Hasher>(&self, state: &mut H) {
                // Model an optional candidate index whose hash gives no
                // information; equality still traverses the complete identity.
                0_u8.hash(state);
            }
        }

        let identity_for = |atom| {
            let policy = CanonicalPolicyNet::<AbstractIndex>::parse(atom).unwrap();
            let projection = super::super::projection::project(
                &policy,
                super::super::projection::DEFAULT_GRAPH_BUDGET,
            )
            .unwrap();
            let prepared = super::super::reconstruct::prepare_reconstruction(&projection).unwrap();
            CanonicalProblemIdentity {
                graph: projection.identity,
                signed: prepared.identity,
            }
        };
        let left = CollidingIdentity(identity_for(Atom::num(2)));
        let right = CollidingIdentity(identity_for(Atom::num(3)));
        assert_ne!(left.0, right.0);

        let mut candidates = HashSet::new();
        assert!(candidates.insert(left.clone()));
        assert!(candidates.insert(right.clone()));
        assert_eq!(candidates.len(), 2);

        let mut history = IdentityHistory::new();
        assert_eq!(history.observe(left.clone()), IdentityObservation::New);
        history.record_retry(RetryReason::IncompleteStabilityCertificate);
        assert_eq!(history.observe(right), IdentityObservation::New);
        history.record_retry(RetryReason::IncompleteStabilityCertificate);
        assert!(matches!(
            history.observe(left),
            IdentityObservation::Cycle {
                first_iteration: 0,
                repeated_iteration: 2,
                cycle_length: 2,
                ..
            }
        ));
    }

    #[test]
    fn canonicalization_orders_semantic_colors_independently_of_symbol_ids() {
        let rep = test_initialize().mink4;
        let tagged_symbol = |name: &str| {
            SymbolBuilder::new(NamespacedSymbol::parse(name))
                .with_tags([SPENSO_TAG.tensor.as_str()])
                .build()
                .unwrap()
        };
        let canonical_line_order = |namespace: &str, reverse_interning_order: bool| {
            let first_name = format!("{namespace}::a");
            let last_name = format!("{namespace}::z");
            let (first, last) = if reverse_interning_order {
                // Symbolica orders Symbols by their process-local numeric ID.
                let last = tagged_symbol(&last_name);
                let first = tagged_symbol(&first_name);
                assert!(last < first);
                (first, last)
            } else {
                let first = tagged_symbol(&first_name);
                let last = tagged_symbol(&last_name);
                assert!(first < last);
                (first, last)
            };
            assert!(
                super::super::semantic::SemanticSymbolKey::from(first)
                    < super::super::semantic::SemanticSymbolKey::from(last)
            );
            let common = tagged_symbol(&format!("{namespace}::common"));
            let first_line = slot!(rep, AbstractIndex::Dummy(91));
            let last_line = slot!(rep, AbstractIndex::Dummy(92));
            let expression = function!(first, first_line.to_atom())
                * function!(common, first_line.to_atom())
                * function!(last, last_line.to_atom())
                * function!(common, last_line.to_atom());
            let canonical = CanonicalPolicyNet::<AbstractIndex>::parse(expression)
                .unwrap()
                .canonize(AbstractIndex::Dummy)
                .unwrap()
                .into_atom();
            let slot_for = |head| {
                let AtomView::Mul(product) = canonical.as_view() else {
                    panic!("the two closed components must remain one Product")
                };
                let argument = product
                    .iter()
                    .find_map(|factor| match factor {
                        AtomView::Fun(function) if function.get_symbol() == head => {
                            function.iter().next()
                        }
                        _ => None,
                    })
                    .expect("the marker tensor must remain in the canonical Product");
                LibrarySlot::<AbstractIndex>::try_from(argument)
                    .unwrap()
                    .aind()
            };
            let first_slot = slot_for(first);
            let last_slot = slot_for(last);
            assert_ne!(first_slot, last_slot);
            first_slot < last_slot
        };

        let reversed = canonical_line_order("canonical_interner_reversed", true);
        let forward = canonical_line_order("canonical_interner_forward", false);
        assert_eq!(reversed, forward);
    }

    #[test]
    fn exact_executed_atom_ends_an_incomplete_power_certificate_in_one_call() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, canonical_driver_exact_power_first);
        let second = slot!(rep, canonical_driver_exact_power_second);
        let tensor = tensor_symbol!(canonical_driver_exact_power_tensor; Symmetric);
        let expression = function!(tensor, first.to_atom(), second.to_atom()).pow(Atom::num(2));
        let retained = CanonicalPolicyNet::<AbstractIndex>::parse(expression)
            .unwrap()
            .canonize(AbstractIndex::Dummy)
            .unwrap()
            .into_atom();

        let probe = CanonicalPolicyNet::<AbstractIndex>::parse(retained.clone()).unwrap();
        let projection = super::super::projection::project(
            &probe,
            super::super::projection::DEFAULT_GRAPH_BUDGET,
        )
        .unwrap();
        let prepared = super::super::reconstruct::prepare_reconstruction(&projection).unwrap();
        let mut allocator = super::super::reconstruct::DummyAllocator::new();
        let reconstructed = super::super::reconstruct::reconstruct(
            &projection,
            prepared,
            &mut allocator,
            &mut AbstractIndex::Dummy,
        )
        .unwrap();
        assert_eq!(
            reconstructed.retry_reason,
            Some(RetryReason::IncompleteStabilityCertificate)
        );
        assert_eq!(execute_atom(reconstructed.network).unwrap(), retained);

        super::super::projection::reset_graphica_calls();
        reset_execution_calls();
        let canonical = CanonicalPolicyNet::<AbstractIndex>::parse(retained.clone())
            .unwrap()
            .canonize(AbstractIndex::Dummy)
            .unwrap();
        assert_eq!(canonical.normalized_atom(), &retained);
        assert_eq!(super::super::projection::graphica_calls(), 1);
        let (full, temporary) = execution_calls();
        assert_eq!(full, 1);
        assert!(temporary > 0);
    }

    #[test]
    fn unsigned_transport_certifies_after_one_graphica_call() {
        let rep = test_initialize().mink4;
        let index = slot!(rep, canonical_driver_unsigned_index);
        let left = tensor_symbol!(canonical_driver_unsigned_left);
        let right = tensor_symbol!(canonical_driver_unsigned_right);
        let expression = function!(left, index.to_atom()) * function!(right, index.to_atom());

        super::super::projection::reset_graphica_calls();
        reset_execution_calls();
        CanonicalPolicyNet::<AbstractIndex>::parse(expression)
            .unwrap()
            .canonize(AbstractIndex::Dummy)
            .unwrap();
        assert_eq!(super::super::projection::graphica_calls(), 1);
        let (full, temporary) = execution_calls();
        assert_eq!(full, 1);
        assert!(temporary > 0);
    }

    #[test]
    fn unchanged_sum_partition_certifies_after_one_graphica_call() {
        let rep = test_initialize().mink4;
        let index = slot!(rep, canonical_driver_sum_index);
        let left = tensor_symbol!(canonical_driver_sum_left);
        let right = tensor_symbol!(canonical_driver_sum_right);
        let expression = function!(left, index.to_atom()) + function!(right, index.to_atom());

        super::super::projection::reset_graphica_calls();
        CanonicalPolicyNet::<AbstractIndex>::parse(expression)
            .unwrap()
            .canonize(AbstractIndex::Dummy)
            .unwrap();
        assert_eq!(super::super::projection::graphica_calls(), 1);
    }

    #[test]
    fn local_zero_removal_exposes_a_global_odd_stabilizer_on_retry() {
        let rep = test_initialize().mink4;
        let i = slot!(rep, canonical_driver_cascade_i);
        let j = slot!(rep, canonical_driver_cascade_j);
        let a = slot!(rep, canonical_driver_cascade_a);
        let b = slot!(rep, canonical_driver_cascade_b);
        let c = slot!(rep, canonical_driver_cascade_c);
        let outer = tensor_symbol!(canonical_driver_cascade_outer; Antisymmetric);
        let cycle = tensor_symbol!(canonical_driver_cascade_cycle; Antisymmetric);
        let left = tensor_symbol!(canonical_driver_cascade_left);
        let right = tensor_symbol!(canonical_driver_cascade_right);
        let breaker_left = tensor_symbol!(canonical_driver_cascade_breaker_left);
        let breaker_right = tensor_symbol!(canonical_driver_cascade_breaker_right);
        let symmetric_branches = function!(left, i.to_atom()) * function!(right, j.to_atom())
            + function!(left, j.to_atom()) * function!(right, i.to_atom());
        let local_zero = function!(cycle, a.to_atom(), b.to_atom())
            * function!(cycle, b.to_atom(), c.to_atom())
            * function!(cycle, c.to_atom(), a.to_atom())
            * function!(breaker_left, i.to_atom())
            * function!(breaker_right, j.to_atom());
        let expression =
            function!(outer, i.to_atom(), j.to_atom()) * (symmetric_branches + local_zero);

        super::super::projection::reset_graphica_calls();
        reset_execution_calls();
        let canonical = CanonicalPolicyNet::<AbstractIndex>::parse(expression)
            .unwrap()
            .canonize(AbstractIndex::Dummy)
            .unwrap();
        assert!(canonical.normalized_atom().is_zero());
        assert_eq!(super::super::projection::graphica_calls(), 2);
        let (full, temporary) = execution_calls();
        assert_eq!(full, 2);
        assert!(temporary > 0);
    }

    #[test]
    fn terminal_zero_and_one_do_not_trigger_a_verification_iteration() {
        let rep = test_initialize().mink4;
        let index = slot!(rep, canonical_driver_terminal_zero_index);
        let antisymmetric = tensor_symbol!(canonical_driver_terminal_zero_tensor; Antisymmetric);
        let cases = [
            function!(antisymmetric, index.to_atom(), index.to_atom()),
            Atom::num(1),
        ];

        for expression in cases {
            super::super::projection::reset_graphica_calls();
            reset_execution_calls();
            let canonical = CanonicalPolicyNet::<AbstractIndex>::parse(expression)
                .unwrap()
                .canonize(AbstractIndex::Dummy)
                .unwrap();
            assert!(canonical.normalized_atom().is_zero() || canonical.normalized_atom().is_one());
            assert_eq!(super::super::projection::graphica_calls(), 1);
            let (full, temporary) = execution_calls();
            assert_eq!(full, 1);
            assert!(temporary > 0);
        }
    }

    #[test]
    fn execution_failures_retain_full_scope_and_result_context() {
        let rep = test_initialize().mink4;
        let index = slot!(rep, canonical_driver_execution_failure_index);
        let left = tensor_symbol!(canonical_driver_execution_failure_left);
        let right = tensor_symbol!(canonical_driver_execution_failure_right);
        let expression = function!(left, index.to_atom()) * function!(right, index.to_atom());

        reset_execution_calls();
        inject_execution_failure(true);
        assert!(matches!(
            CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone())
                .unwrap()
                .canonize(AbstractIndex::Dummy),
            Err(CanonicalizationError::Execution { scope, error })
                if scope.starts_with("canonical scope ")
                    && error == "injected execution failure"
        ));

        reset_execution_calls();
        inject_execution_failure(false);
        assert!(matches!(
            CanonicalPolicyNet::<AbstractIndex>::parse(expression)
                .unwrap()
                .canonize(AbstractIndex::Dummy),
            Err(CanonicalizationError::Execution { scope, error })
                if scope == "complete network" && error == "injected execution failure"
        ));

        let mut malformed = crate::tensor::SymbolicNet::<AbstractIndex>::from_scalar(Atom::num(1));
        let (_, _, node) = malformed.graph.graph.iter_nodes_mut().next().unwrap();
        *node = NetworkNode::Leaf(NetworkLeaf::TensorSum(Vec::new()));
        assert_eq!(
            execute_atom(malformed),
            Err(CanonicalizationError::Execution {
                scope: "complete network result".to_owned(),
                error: "empty tensor sum leaf".to_owned(),
            })
        );

        reset_execution_calls();
        inject_result_extraction_failure(true);
        let power = tensor_symbol!(canonical_driver_result_failure_power; Symmetric);
        let power_expression = function!(power, index.to_atom(), index.to_atom()).pow(Atom::num(2));
        assert!(matches!(
            CanonicalPolicyNet::<AbstractIndex>::parse(power_expression)
                .unwrap()
                .canonize(AbstractIndex::Dummy),
            Err(CanonicalizationError::Execution { scope, error })
                if (scope.starts_with("canonical scope ") || scope.starts_with("Power "))
                    && error == "injected result extraction failure"
        ));
    }

    #[test]
    fn final_policy_network_and_atom_projections_agree_and_are_idempotent() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, AbstractIndex::Dummy(0));
        let second = slot!(rep, AbstractIndex::Dummy(1));
        let antisymmetric = tensor_symbol!(canonical_driver_projection_tensor; Antisymmetric);
        let left = tensor_symbol!(canonical_driver_projection_left);
        let right = tensor_symbol!(canonical_driver_projection_right);
        let expression = Atom::num(2)
            * function!(antisymmetric, second.to_atom(), first.to_atom())
            * function!(left, first.to_atom())
            * function!(right, second.to_atom());
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let cooked = cooking.try_cook_indices(expression.as_view()).unwrap();

        let canonical = CanonicalPolicyNet::<AbstractIndex>::parse(cooked)
            .unwrap()
            .canonize(AbstractIndex::Dummy)
            .unwrap();
        let canonical_network = canonical.network().snapshot_dot();
        let retained = canonical.normalized_atom().clone();
        assert_eq!(execute_atom(canonical.into_network()).unwrap(), retained);

        let facade = expression
            .as_view()
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();
        assert_eq!(facade, cooking.uncook(retained.as_view()));
        assert_eq!(
            facade
                .as_view()
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap(),
            facade
        );

        let reparsed = CanonicalPolicyNet::<AbstractIndex>::parse(retained.clone())
            .unwrap()
            .canonize(AbstractIndex::Dummy)
            .unwrap();
        assert_eq!(reparsed.network().snapshot_dot(), canonical_network);
        assert_eq!(reparsed.normalized_atom(), &retained);
        assert_eq!(execute_atom(reparsed.into_network()).unwrap(), retained);
    }

    #[test]
    fn branch_local_dummy_reuse_conservatively_declines_the_certificate() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, AbstractIndex::Dummy(0));
        let second = slot!(rep, AbstractIndex::Dummy(1));
        let left = tensor_symbol!(canonical_driver_certificate_left);
        let right = tensor_symbol!(canonical_driver_certificate_right);
        let expression = function!(left, first.to_atom(), first.to_atom())
            + function!(right, second.to_atom(), second.to_atom());
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();

        super::super::projection::reset_graphica_calls();
        reset_execution_calls();
        let projection = super::super::projection::project(
            &policy,
            super::super::projection::DEFAULT_GRAPH_BUDGET,
        )
        .unwrap();
        assert_eq!(execution_calls(), (0, 0));
        let prepared = super::super::reconstruct::prepare_reconstruction(&projection).unwrap();
        let (full, temporary) = execution_calls();
        assert_eq!(full, 0);
        assert!(temporary > 0);
        let mut allocator = super::super::reconstruct::DummyAllocator::new();
        let mut new_dummy = AbstractIndex::Dummy;
        let reconstructed = super::super::reconstruct::reconstruct(
            &projection,
            prepared,
            &mut allocator,
            &mut new_dummy,
        )
        .unwrap();
        assert_eq!(
            reconstructed.retry_reason,
            Some(RetryReason::VisibleClassMergeSplit)
        );
        assert_eq!(super::super::projection::graphica_calls(), 1);
        assert_ne!(
            execute_atom(reconstructed.network).unwrap(),
            *policy.normalized_atom()
        );

        super::super::projection::reset_graphica_calls();
        reset_execution_calls();
        policy.canonize(AbstractIndex::Dummy).unwrap();
        assert_eq!(super::super::projection::graphica_calls(), 2);
        let (full, temporary) = execution_calls();
        assert_eq!(full, 1);
        assert!(temporary > 0);
    }

    #[test]
    fn visible_sign_edit_uses_conditional_second_iteration() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, AbstractIndex::Dummy(0));
        let second = slot!(rep, AbstractIndex::Dummy(1));
        let tensor = tensor_symbol!(canonical_driver_signed_tensor; Antisymmetric);
        let left = tensor_symbol!(canonical_driver_signed_left);
        let right = tensor_symbol!(canonical_driver_signed_right);
        let expression = Atom::num(2)
            * function!(tensor, first.to_atom(), second.to_atom())
            * function!(left, second.to_atom())
            * function!(right, first.to_atom());
        let expected = Atom::num(-2)
            * function!(tensor, first.to_atom(), second.to_atom())
            * function!(left, first.to_atom())
            * function!(right, second.to_atom());

        super::super::projection::reset_graphica_calls();
        reset_execution_calls();
        let canonical = CanonicalPolicyNet::<AbstractIndex>::parse(expression)
            .unwrap()
            .canonize(AbstractIndex::Dummy)
            .unwrap();
        assert_eq!(canonical.normalized_atom(), &expected);
        assert_eq!(super::super::projection::graphica_calls(), 2);
        let (full, temporary) = execution_calls();
        assert_eq!(full, 2);
        assert!(temporary > 0);
    }

    #[test]
    fn request_local_dummy_allocator_reuses_values_across_retry() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, AbstractIndex::Dummy(0));
        let second = slot!(rep, AbstractIndex::Dummy(1));
        let tensor = tensor_symbol!(canonical_driver_retry_allocator_tensor; Antisymmetric);
        let left = tensor_symbol!(canonical_driver_retry_allocator_left);
        let right = tensor_symbol!(canonical_driver_retry_allocator_right);
        let expression = Atom::num(2)
            * function!(tensor, first.to_atom(), second.to_atom())
            * function!(left, second.to_atom())
            * function!(right, first.to_atom());
        let mut allocated_positions = Vec::new();

        super::super::projection::reset_graphica_calls();
        CanonicalPolicyNet::<AbstractIndex>::parse(expression)
            .unwrap()
            .canonize(|position| {
                allocated_positions.push(position);
                AbstractIndex::Dummy(position + 10)
            })
            .unwrap();

        assert_eq!(super::super::projection::graphica_calls(), 2);
        assert_eq!(allocated_positions, vec![0, 1]);
    }
}
