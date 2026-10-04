//! Optional diagnostic wall clocks for the shared implementation. Captures are
//! thread-local and inactive outside `measure`; worker-thread CPU is not summed.
//! Inclusive phase totals overlap. `outermost_ns` avoids counting recursion of
//! the same phase twice, but still includes other phases nested within it.

use std::{cell::RefCell, fmt, time::Instant};

use crate::tensor::ReductionStatus;

#[derive(Clone, Copy)]
pub(crate) enum Phase {
    Admission,
    InterfaceDiscovery,
    Validation,
    CandidateScan,
    Planning,
    ShallowGraph,
    GammaKernel,
    ColorKernel,
    EpsilonKernel,
    Contraction,
    DummyReservation,
    Materialization,
    OutputNormalization,
    ContractionStates,
    CoefficientMerge,
    ContractionEmission,
}

const PHASES: [&str; 16] = [
    "admission",
    "interface_discovery",
    "validation",
    "candidate_scan",
    "planning",
    "shallow_graph",
    "gamma_kernel",
    "color_kernel",
    "epsilon_kernel",
    "contraction",
    "dummy_reservation",
    "materialization",
    "output_normalization",
    "contraction_states",
    "coefficient_merge",
    "contraction_emission",
];

/// Storage and branching at one selected contraction factor. Byte sizes describe
/// stored coefficient buffers and estimated workspace, not resident memory.
#[derive(Clone, Copy)]
pub(crate) struct Frontier {
    pub selected: usize,
    pub states: usize,
    pub local_terms: usize,
    pub generated: usize,
    pub predicted_bytes: usize,
    pub coefficient_bytes: usize,
    pub max_coefficient_bytes: usize,
}

#[derive(Clone, Copy)]
pub(crate) enum ScopeOpenReason {
    Unobserved,
    OpaqueOrBracket,
    MetricOrVector,
}

/// Requests may reuse an already opened scope. Indexed powers overlap the
/// exclusive reasons; metric/vector presence takes precedence over opacity.
#[derive(Default)]
struct ScopeOpenRequests {
    reasons: [usize; 3],
    indexed_powers: usize,
}

#[derive(Clone, Copy, Default)]
struct Entry {
    calls: usize,
    inclusive_ns: u128,
    outermost_ns: u128,
    depth: usize,
}

#[derive(Default)]
struct Capture {
    entries: [Entry; PHASES.len()],
    outcomes: Vec<(&'static str, ReductionStatus)>,
    frontiers: Vec<Frontier>,
    scope_open_requests: ScopeOpenRequests,
}

thread_local! {
    static CAPTURE: RefCell<Option<Capture>> = const { RefCell::new(None) };
}

pub(crate) fn outcome(operation: &'static str, status: ReductionStatus) {
    CAPTURE.with_borrow_mut(|capture| {
        if let Some(capture) = capture {
            capture.outcomes.push((operation, status));
        }
    });
}

pub(crate) fn frontier(observation: Frontier) {
    CAPTURE.with_borrow_mut(|capture| {
        if let Some(capture) = capture {
            capture.frontiers.push(observation);
        }
    });
}

pub(crate) fn frontier_workspace(bytes: usize) {
    CAPTURE.with_borrow_mut(|capture| {
        if let Some(frontier) = capture
            .as_mut()
            .and_then(|capture| capture.frontiers.last_mut())
        {
            frontier.predicted_bytes = frontier.predicted_bytes.max(bytes);
        }
    });
}

pub(crate) fn scope_open(reason: ScopeOpenReason, indexed_power: bool) {
    CAPTURE.with_borrow_mut(|capture| {
        if let Some(capture) = capture {
            capture.scope_open_requests.reasons[reason as usize] += 1;
            capture.scope_open_requests.indexed_powers += usize::from(indexed_power);
        }
    });
}

pub(crate) struct Scope {
    phase: usize,
    start: Instant,
    outermost: bool,
}

pub(crate) fn scope(phase: Phase) -> Option<Scope> {
    CAPTURE.with_borrow_mut(|capture| {
        let entry = &mut capture.as_mut()?.entries[phase as usize];
        let outermost = entry.depth == 0;
        entry.calls += 1;
        entry.depth += 1;
        Some(Scope {
            phase: phase as usize,
            start: Instant::now(),
            outermost,
        })
    })
}

impl Drop for Scope {
    fn drop(&mut self) {
        let elapsed = self.start.elapsed().as_nanos();
        CAPTURE.with_borrow_mut(|capture| {
            if let Some(capture) = capture {
                let entry = &mut capture.entries[self.phase];
                entry.depth -= 1;
                entry.inclusive_ns += elapsed;
                if self.outermost {
                    entry.outermost_ns += elapsed;
                }
            }
        });
    }
}

/// One instrumented diagnostic invocation, separate from primary benchmark clocks.
pub struct Measurements {
    elapsed_ns: u128,
    entries: [Entry; PHASES.len()],
    outcomes: Vec<(&'static str, ReductionStatus)>,
    frontiers: Vec<Frontier>,
    scope_open_requests: ScopeOpenRequests,
}

impl fmt::Display for Measurements {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{{\"clock\":\"same_thread_wall_nested_inclusive\",\"elapsed_ns\":{},\"phases\":{{",
            self.elapsed_ns
        )?;
        for (index, (name, entry)) in PHASES.iter().zip(&self.entries).enumerate() {
            if index != 0 {
                f.write_str(",")?;
            }
            write!(
                f,
                "\"{name}\":{{\"calls\":{},\"inclusive_ns\":{},\"outermost_ns\":{}}}",
                entry.calls, entry.inclusive_ns, entry.outermost_ns
            )?;
        }
        f.write_str("},\"outcomes\":[")?;
        for (index, (operation, status)) in self.outcomes.iter().enumerate() {
            if index != 0 {
                f.write_str(",")?;
            }
            write!(
                f,
                "{{\"operation\":\"{operation}\",\"status\":\"{status:?}\"}}"
            )?;
        }
        f.write_str("],\"frontiers\":[")?;
        for (index, frontier) in self.frontiers.iter().enumerate() {
            if index != 0 {
                f.write_str(",")?;
            }
            write!(
                f,
                "{{\"selected\":{},\"states\":{},\"local_terms\":{},\"generated\":{},\"predicted_bytes\":{},\"coefficient_bytes\":{},\"max_coefficient_bytes\":{} }}",
                frontier.selected,
                frontier.states,
                frontier.local_terms,
                frontier.generated,
                frontier.predicted_bytes,
                frontier.coefficient_bytes,
                frontier.max_coefficient_bytes,
            )?;
        }
        let [unobserved, opaque_or_bracket, metric_or_vector] = self.scope_open_requests.reasons;
        write!(
            f,
            "],\"scope_open_requests\":{{\"unobserved\":{unobserved},\"opaque_or_bracket\":{opaque_or_bracket},\"metric_or_vector\":{metric_or_vector},\"indexed_powers\":{}}}}}",
            self.scope_open_requests.indexed_powers,
        )
    }
}

/// Capture nested phase clocks on this thread. The capture is cleared on panic.
pub fn measure<T>(label: &str, operation: impl FnOnce() -> T) -> (T, Measurements) {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            CAPTURE.with_borrow_mut(|capture| *capture = None);
        }
    }
    CAPTURE.with_borrow_mut(|capture| {
        assert!(capture.is_none(), "diagnostic captures cannot nest");
        *capture = Some(Capture::default());
    });
    let reset = Reset;
    // The existing process-wide Spenso profiler is only enabled by its own
    // environment switch. Run those captures in a separate diagnostic process;
    // its atomics and reporting are outside the primary benchmark protocol.
    spenso::network::profile::reset();
    let start = Instant::now();
    let result = operation();
    let elapsed_ns = start.elapsed().as_nanos();
    let capture = CAPTURE.with_borrow_mut(|capture| capture.take().unwrap());
    drop(reset);
    spenso::network::profile::report(&format!("idenso_diagnostic {label}"));
    (
        result,
        Measurements {
            elapsed_ns,
            entries: capture.entries,
            outcomes: capture.outcomes,
            frontiers: capture.frontiers,
            scope_open_requests: capture.scope_open_requests,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tensor::{ContractSettings, SymbolicTensor};
    use spenso::{g, mink, p};

    #[test]
    fn diagnostic_capture_preserves_the_contracted_tensor_and_closes_nested_scopes() {
        crate::test_support::test_initialize();
        let expression = g!(mink!(4, 921001), mink!(4, 921003)) * p!(mink!(4, 921003));
        let (result, clocks) = measure("complete contraction", || {
            SymbolicTensor::infer(expression)
                .unwrap()
                .contract(ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })
                .unwrap()
        });
        assert_eq!(result.expression(), &p!(mink!(4, 921001)));
        assert_eq!(clocks.entries[Phase::Admission as usize].calls, 1);
        assert!(clocks.entries[Phase::InterfaceDiscovery as usize].calls > 0);
        assert!(clocks.entries[Phase::Contraction as usize].calls > 0);
        assert_eq!(clocks.entries[Phase::Materialization as usize].calls, 0);
        assert_eq!(clocks.outcomes, [("contract", ReductionStatus::Complete)]);
        for entry in clocks.entries {
            assert_eq!(entry.depth, 0);
            assert!(entry.outermost_ns <= entry.inclusive_ns);
        }
        assert!(scope(Phase::Admission).is_none());
    }

    #[test]
    fn diagnostic_outcomes_distinguish_cached_completion_from_exact_capping() {
        crate::test_support::test_initialize();
        let source =
            SymbolicTensor::infer(g!(mink!(4, 921011), mink!(4, 921013)) * p!(mink!(4, 921013)))
                .unwrap();
        let settings = ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        };
        let owner = source.contract(settings).unwrap();
        let ((rerun, capped), clocks) = measure("completion outcomes", || {
            (
                owner.contract(settings).unwrap(),
                source
                    .contract(ContractSettings {
                        max_passes: Some(0),
                        ..settings
                    })
                    .unwrap(),
            )
        });
        assert!(std::sync::Arc::ptr_eq(
            owner.reduction_observations(),
            rerun.reduction_observations()
        ));
        assert_eq!(capped.expression(), source.expression());
        assert_eq!(
            clocks.outcomes,
            [
                ("contract", ReductionStatus::Complete),
                ("contract", ReductionStatus::Capped),
            ]
        );
    }

    #[test]
    fn diagnostic_capture_is_cleared_after_an_uncontrolled_panic() {
        assert!(
            std::panic::catch_unwind(|| measure("callback panic", || {
                let _scope = scope(Phase::GammaKernel);
                panic!("uncontrolled diagnostic callback");
            }))
            .is_err()
        );
        let (value, clocks) = measure("after callback panic", || 7);
        assert_eq!(value, 7);
        assert!(clocks.entries.iter().all(|entry| entry.calls == 0));
        assert!(clocks.outcomes.is_empty());
    }
}
