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
    projection::{
        CanonicalProjection, DEFAULT_GRAPH_BUDGET, ExternalLineMode, GraphBudget, ProblemIdentity,
        project,
    },
    reconstruct::{
        DummyAllocator, PreparedReconstruction, RetryReason, SignedProblemIdentity,
        prepare_reconstruction, reconstruct,
    },
    semantic::SemanticAtomKey,
};
use crate::tensor::SymbolicTensor;

// See the Phase 6 budget table in the signed-canonicalization architecture plan.
pub(super) const DEFAULT_ITERATION_LIMIT: usize = 8;

pub(super) struct CanonicalRequest<'a, Aind: AbsInd, F> {
    new_dummy: &'a mut F,
    dummy_allocator: DummyAllocator<Aind>,
    temporary_dummy_next: usize,
    iteration_limit: usize,
    graph_budget: GraphBudget,
    external_line_mode: ExternalLineMode,
}

impl<'a, Aind: AbsInd, F> CanonicalRequest<'a, Aind, F> {
    fn new(new_dummy: &'a mut F, iteration_limit: usize, graph_budget: GraphBudget) -> Self {
        Self {
            new_dummy,
            dummy_allocator: DummyAllocator::new(),
            temporary_dummy_next: 0,
            iteration_limit,
            graph_budget,
            external_line_mode: ExternalLineMode::Preserve,
        }
    }
}

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

struct PreparedGraphStep<Aind: AbsInd> {
    projection: CanonicalProjection<Aind>,
    prepared: PreparedReconstruction,
}

impl<Aind: AbsInd> PreparedGraphStep<Aind> {
    fn identity(&self) -> CanonicalProblemIdentity {
        CanonicalProblemIdentity {
            graph: self.projection.identity.clone(),
            signed: self.prepared.identity.clone(),
        }
    }
}

struct CompletedGraphStep<Aind: AbsInd> {
    /// Reparsed graph result, absent when execution reproduced the input Atom
    /// exactly and the caller can retain its parser-proven policy.
    policy: Option<CanonicalPolicyNet<Aind>>,
    retry_reason: Option<RetryReason>,
    terminal: bool,
    execution_unchanged: bool,
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

// Boxing the policy variant would allocate on the ordinary canonicalization
// path merely to shrink this short-lived result enum.
#[allow(clippy::large_enum_variant)]
enum Canonicalized<Aind: AbsInd> {
    Policy(CanonicalPolicyNet<Aind>),
    Terminal(Atom),
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
    #[cfg(test)]
    pub(crate) fn canonize<F>(self, mut new_dummy: F) -> Result<Self, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
    {
        match self.canonize_output(&mut new_dummy)? {
            Canonicalized::Policy(policy) => Ok(policy),
            Canonicalized::Terminal(atom) => CanonicalPolicyNet::parse(atom),
        }
    }

    pub(crate) fn canonize_atom<F>(self, mut new_dummy: F) -> Result<Atom, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
    {
        Ok(match self.canonize_output(&mut new_dummy)? {
            Canonicalized::Policy(policy) => policy.into_atom(),
            Canonicalized::Terminal(atom) => atom,
        })
    }

    fn canonize_output<F>(
        self,
        new_dummy: &mut F,
    ) -> Result<Canonicalized<Aind>, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
    {
        let mut request =
            CanonicalRequest::new(new_dummy, DEFAULT_ITERATION_LIMIT, DEFAULT_GRAPH_BUDGET);
        match super::young::straighten(self, &mut request)? {
            super::young::YoungStraightened::Terminal(atom) => Ok(Canonicalized::Terminal(atom)),
            super::young::YoungStraightened::Policy {
                policy,
                graph_canonical,
            } => {
                if graph_canonical {
                    Ok(Canonicalized::Policy(policy))
                } else {
                    request.canonize_graph(policy).map(Canonicalized::Policy)
                }
            }
        }
    }

    #[cfg(test)]
    fn canonize_with<F>(
        self,
        new_dummy: &mut F,
        iteration_limit: usize,
        graph_budget: super::projection::GraphBudget,
    ) -> Result<Self, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
    {
        CanonicalRequest::new(new_dummy, iteration_limit, graph_budget).canonize_graph(self)
    }
}

impl<Aind, F> CanonicalRequest<'_, Aind, F>
where
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
    F: FnMut(usize) -> Aind,
    SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
{
    /// Reserve one deterministic source-disjoint namespace for a candidate batch.
    pub(super) fn reserve_temporary_namespace(
        &mut self,
        source_indices: usize,
    ) -> Result<usize, CanonicalizationError> {
        // One graph can skip at most one currently visible name per line and
        // then allocate at most one replacement per line. Source-wide names
        // are an additional conservative forbidden set.
        let width = self
            .graph_budget
            .vertices
            .checked_mul(2)
            .and_then(|width| width.checked_add(source_indices))
            .and_then(|width| width.checked_add(1))
            .ok_or_else(|| {
                CanonicalizationError::Projection(
                    "temporary canonical namespace width overflowed usize".to_owned(),
                )
            })?;
        let start = self.temporary_dummy_next;
        self.temporary_dummy_next = start.checked_add(width).ok_or_else(|| {
            CanonicalizationError::Projection(
                "temporary canonical namespace allocation overflowed usize".to_owned(),
            )
        })?;
        Ok(start)
    }

    pub(super) fn canonize_graph(
        &mut self,
        mut policy: CanonicalPolicyNet<Aind>,
    ) -> Result<CanonicalPolicyNet<Aind>, CanonicalizationError> {
        let mut identities = IdentityHistory::<CanonicalProblemIdentity>::new();
        let mut last_reason = "initial canonical projection".to_owned();

        for _ in 0..self.iteration_limit {
            let step = self.prepare_graph_step(&policy)?;
            match identities.observe(step.identity()) {
                IdentityObservation::Consecutive => return Ok(policy),
                cycle @ IdentityObservation::Cycle { .. } => {
                    return Err(cycle.into_cycle_error());
                }
                IdentityObservation::New => {}
            }
            let completed = self.complete_graph_step(&policy, step)?;
            if completed.execution_unchanged {
                return Ok(policy);
            }
            let next = completed
                .policy
                .expect("changed graph execution is reparsed");
            if completed.terminal {
                return Ok(next);
            }
            let Some(retry_reason) = completed.retry_reason else {
                return Ok(next);
            };
            last_reason = retry_reason.to_string();
            identities.record_retry(retry_reason);
            policy = next;
        }

        Err(CanonicalizationError::IterationLimit {
            limit: self.iteration_limit,
            last_reason,
        })
    }

    /// Canonize modulo an exact deterministic normalizer.
    ///
    /// `policy` is already normalized. Each transition applies one complete
    /// graph step and then re-applies `normalize`. Only consecutive exact
    /// equality of the post-normalizer states certifies the composite fixed
    /// point. The retained middle graph result is returned, so a fresh request
    /// reproduces that graph-after-normalizer representative exactly: if
    /// `x = normalize(input)`, `y = graph(x)`, and `normalize(y) = x`, then the
    /// next request repeats the same `x -> y` transition.
    pub(super) fn canonize_graph_modulo<N>(
        &mut self,
        mut policy: CanonicalPolicyNet<Aind>,
        mut normalize: N,
    ) -> Result<CanonicalPolicyNet<Aind>, CanonicalizationError>
    where
        N: FnMut(
            &CanonicalPolicyNet<Aind>,
        ) -> Result<CanonicalPolicyNet<Aind>, CanonicalizationError>,
    {
        let mut prepared = self.prepare_graph_step(&policy)?;
        let mut atoms = vec![policy.normalized_atom().clone()];
        let mut retry_reasons = Vec::<String>::new();
        let mut last_reason = "initial canonical projection".to_owned();

        for iteration in 0..self.iteration_limit {
            let completed = self.complete_graph_step(&policy, prepared)?;
            let graph = if completed.execution_unchanged {
                policy
            } else {
                completed
                    .policy
                    .expect("changed graph execution is reparsed")
            };
            let normalized = normalize(&graph)?;
            if normalized.normalized_atom().as_view().is_zero()
                || normalized.normalized_atom().as_view().is_one()
            {
                return Ok(normalized);
            }

            last_reason = completed.retry_reason.map_or_else(
                || "composite normalizer reapplication".to_owned(),
                |reason| reason.to_string(),
            );
            retry_reasons.push(last_reason.clone());
            let repeated_iteration = atoms.len();
            let atom = normalized.normalized_atom();
            if atoms.last() == Some(atom) {
                return Ok(graph);
            }
            let repeated_atom = atoms.iter().position(|candidate| candidate == atom);
            if let Some(first_iteration) = repeated_atom {
                return Err(CanonicalizationError::ConvergenceCycle {
                    first_iteration,
                    repeated_iteration,
                    cycle_length: repeated_iteration - first_iteration,
                    retry_reasons: retry_reasons[first_iteration..repeated_iteration].to_vec(),
                });
            }
            atoms.push(atom.clone());
            policy = normalized;
            if iteration + 1 == self.iteration_limit {
                break;
            }
            prepared = self.prepare_graph_step(&policy)?;
        }

        Err(CanonicalizationError::IterationLimit {
            limit: self.iteration_limit,
            last_reason,
        })
    }

    fn prepare_graph_step(
        &self,
        policy: &CanonicalPolicyNet<Aind>,
    ) -> Result<PreparedGraphStep<Aind>, CanonicalizationError> {
        let projection = project(policy, self.graph_budget, self.external_line_mode)?;
        let prepared = prepare_reconstruction(&projection)?;
        Ok(PreparedGraphStep {
            projection,
            prepared,
        })
    }

    fn complete_graph_step(
        &mut self,
        policy: &CanonicalPolicyNet<Aind>,
        step: PreparedGraphStep<Aind>,
    ) -> Result<CompletedGraphStep<Aind>, CanonicalizationError> {
        let reconstruction = reconstruct(
            &step.projection,
            step.prepared,
            &mut self.dummy_allocator,
            self.new_dummy,
        )?;
        let mut retry_reason = reconstruction.retry_reason;
        let operation_profile = (retry_reason == Some(RetryReason::IncompleteStabilityCertificate))
            .then(|| Self::operation_profile(&reconstruction.network));
        let atom = execute_atom(reconstruction.network)?;
        let atom = if self.external_line_mode == ExternalLineMode::AnonymousEvenPower {
            // This mode is used only for a base bound by an even Power, so
            // its two global orientations represent the same final value.
            let negative = -atom.clone();
            if SemanticAtomKey::new(negative.as_view()) < SemanticAtomKey::new(atom.as_view()) {
                negative
            } else {
                atom
            }
        } else {
            atom
        };
        let terminal = atom.as_view().is_zero() || atom.as_view().is_one();
        // The parser and projection are deterministic, so exact equality after
        // mandatory execution is a complete extensional graph fixed point.
        let execution_unchanged = atom == *policy.normalized_atom();
        if execution_unchanged {
            return Ok(CompletedGraphStep {
                policy: None,
                retry_reason,
                terminal,
                execution_unchanged,
            });
        }
        let reparsed = CanonicalPolicyNet::parse(atom)?;
        if operation_profile
            .is_some_and(|profile| profile != Self::operation_profile(reparsed.network()))
        {
            retry_reason = Some(RetryReason::ExecutionTopologyChange);
        }
        Ok(CompletedGraphStep {
            policy: Some(reparsed),
            retry_reason,
            terminal,
            execution_unchanged,
        })
    }

    fn operation_profile(network: &crate::tensor::SymbolicNet<Aind>) -> [usize; 5] {
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
    }

    /// Canonize a staged graph in a source-disjoint temporary dummy namespace.
    ///
    /// The external-line mode either retains the currently exposed interface
    /// or chooses the anonymous frame of an even-Power base.
    pub(super) fn canonize_temporary_graph(
        &mut self,
        policy: CanonicalPolicyNet<Aind>,
        source_indices: &[Atom],
        namespace_start: usize,
        external_line_mode: ExternalLineMode,
    ) -> Result<CanonicalPolicyNet<Aind>, CanonicalizationError> {
        let mut forbidden = source_indices.to_vec();
        forbidden.extend(policy.network().store.tensors.iter().flat_map(|tensor| {
            (&tensor.structure)
                .into_iter()
                .map(|slot| slot.aind.to_atom())
        }));
        let mut next = namespace_start;
        let mut temporary_dummy = |_: usize| loop {
            let dummy = Aind::new_dummy_at(next);
            next += 1;
            let atom = dummy.to_atom();
            if !forbidden.contains(&atom) {
                forbidden.push(atom);
                break dummy;
            }
        };
        let mut request = CanonicalRequest::new(
            &mut temporary_dummy,
            self.iteration_limit,
            self.graph_budget,
        );
        request.external_line_mode = external_line_mode;
        request.canonize_graph(policy)
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
                super::super::projection::ExternalLineMode::Preserve,
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
            super::super::projection::ExternalLineMode::Preserve,
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
    fn anonymous_boundary_frame_is_alpha_invariant_without_adding_symmetry() {
        let rep = test_initialize().mink4;
        let [first, second, renamed_first, renamed_second] = [
            slot!(rep, canonical_driver_anonymous_first),
            slot!(rep, canonical_driver_anonymous_second),
            slot!(rep, canonical_driver_anonymous_renamed_first),
            slot!(rep, canonical_driver_anonymous_renamed_second),
        ]
        .map(|slot| slot.to_lib());
        let tensor = tensor_symbol!(canonical_driver_anonymous_tensor; Antisymmetric);
        let canonize = |left: LibrarySlot<AbstractIndex>, right: LibrarySlot<AbstractIndex>| {
            let expression = function!(tensor, left.to_atom(), right.to_atom());
            let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
            let source_indices = [left.aind().to_atom(), right.aind().to_atom()];
            let mut new_dummy = AbstractIndex::Dummy;
            let mut request = CanonicalRequest::new(
                &mut new_dummy,
                DEFAULT_ITERATION_LIMIT,
                super::super::projection::DEFAULT_GRAPH_BUDGET,
            );
            let namespace = request
                .reserve_temporary_namespace(source_indices.len())
                .unwrap();
            request
                .canonize_temporary_graph(
                    policy,
                    &source_indices,
                    namespace,
                    ExternalLineMode::AnonymousEvenPower,
                )
                .unwrap()
                .into_atom()
        };

        let canonical = canonize(first, second);
        assert!(!canonical.is_zero());
        assert_eq!(canonical, canonize(renamed_first, renamed_second));
        assert!(canonize(first, first).is_zero());
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
        assert_eq!(temporary, 1);
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
            super::super::projection::ExternalLineMode::Preserve,
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
