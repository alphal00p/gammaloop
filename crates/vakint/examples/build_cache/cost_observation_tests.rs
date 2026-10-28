use super::*;
use rustred::{
    family::IntegralKey,
    solver::{CandidateRoutedApplicationSnapshot, CandidateRoutedApplicationStatistics},
};
use std::{collections::HashSet, time::Duration};

struct Harness {
    observation: CostObservation,
    sequence: u64,
    local_sequence: u64,
    completed: usize,
    unique: HashSet<(bool, Vec<i64>)>,
}

impl Harness {
    fn new() -> Self {
        Self {
            observation: CostObservation::new(json!({"test_binding": true})),
            sequence: 0,
            local_sequence: 0,
            completed: 0,
            unique: HashSet::new(),
        }
    }

    fn start(&mut self, now: u64) {
        self.local_sequence = 0;
        self.completed = 0;
        self.unique.clear();
        self.send(Event::Snapshot, Phase::Admission, now, 0, None)
            .unwrap();
    }

    fn send(
        &mut self,
        kind: Event,
        phase: Phase,
        now: u64,
        cost: usize,
        current: Option<Work<'_>>,
    ) -> Result<Control, ()> {
        self.send_family(kind, phase, now, cost, current, "test-family")
    }

    fn send_family(
        &mut self,
        kind: Event,
        phase: Phase,
        now: u64,
        cost: usize,
        current: Option<Work<'_>>,
        family: &str,
    ) -> Result<Control, ()> {
        self.sequence += 1;
        self.local_sequence += 1;
        if matches!(kind, Event::NodeCompleted(_)) {
            self.completed += 1;
            let work = current.unwrap();
            self.unique.insert((
                matches!(work, Work::Apply { .. }),
                work.target().powers().to_vec(),
            ));
        }
        let mut envelope = Envelope::new(kind, family);
        envelope.sequence = self.sequence;
        envelope.elapsed = Duration::from_nanos(now);
        envelope.activity = match phase {
            Phase::Expansion => Activity::Expansion,
            Phase::Combination => Activity::CoefficientCombination,
            _ => Activity::Admission,
        };
        envelope.current = current;
        envelope.status = if kind == Event::CoreFinished {
            Status::Succeeded
        } else {
            Status::Running
        };
        envelope.metrics.completion_attempts = self.completed;
        envelope.metrics.completed_unique_operational_nodes = self.unique.len();
        envelope.snapshot = Some(CandidateRoutedApplicationSnapshot {
            sequence: self.local_sequence,
            elapsed: Duration::from_nanos(now),
            phase,
            status: envelope.status,
            statistics: CandidateRoutedApplicationStatistics {
                transport_calls: cost,
                ..Default::default()
            },
            unique_nodes: self.unique.len(),
            active_nodes: 0,
            pending_frames: 0,
            cached: Default::default(),
            pinned: None,
            eligible: None,
        });
        self.observation.observe(&envelope)
    }

    fn finish(self) -> Value {
        self.observation.finish(true, ObserverStatus::Active)
    }
}

fn transport(value: &Value) -> u64 {
    value[3].as_u64().unwrap()
}

#[test]
fn nested_and_shared_visits_have_disjoint_self_costs_and_intervals() {
    let root = IntegralKey::try_new([2, 1]).unwrap();
    let shared = IntegralKey::try_new([1, 1]).unwrap();
    let mut h = Harness::new();
    h.start(0);
    h.send(
        Event::Snapshot,
        Phase::Expansion,
        1,
        1,
        Some(Work::Route(&root)),
    )
    .unwrap();
    h.send(
        Event::Snapshot,
        Phase::Expansion,
        3,
        3,
        Some(Work::Route(&shared)),
    )
    .unwrap();
    h.send(
        Event::NodeCompleted(Source::Computed),
        Phase::Expansion,
        8,
        8,
        Some(Work::Route(&shared)),
    )
    .unwrap();
    // Same key is attempted again after residency has been established. Its
    // own hit cost is observed independently of the first computed interval.
    h.send(
        Event::Snapshot,
        Phase::Expansion,
        10,
        10,
        Some(Work::Route(&shared)),
    )
    .unwrap();
    h.send(
        Event::NodeCompleted(Source::RamHit),
        Phase::Expansion,
        11,
        11,
        Some(Work::Route(&shared)),
    )
    .unwrap();
    h.send(
        Event::Snapshot,
        Phase::Combination,
        12,
        12,
        Some(Work::Route(&root)),
    )
    .unwrap();
    h.send(
        Event::NodeCompleted(Source::Computed),
        Phase::Combination,
        15,
        15,
        Some(Work::Route(&root)),
    )
    .unwrap();
    h.send(Event::CoreFinished, Phase::RootCollection, 18, 18, None)
        .unwrap();
    let document = h.finish();
    assert_eq!(document["complete"], true);
    let attempts = document["attempts"].as_array().unwrap();
    assert_eq!(attempts.len(), 3);
    assert_eq!(attempts[1]["parent_id"], 0);
    assert_eq!(attempts[2]["parent_id"], 0);
    assert_eq!(transport(&attempts[0]["inclusive_counters"]), 14);
    assert_eq!(transport(&attempts[0]["self_counters"]), 8);
    assert_eq!(transport(&attempts[1]["self_counters"]), 5);
    assert_eq!(transport(&attempts[2]["self_counters"]), 1);
    assert_eq!(attempts[1]["source"], "computed");
    assert_eq!(attempts[2]["source"], "ram_hit");
    assert_eq!(document["keys"][1]["visits"], 2);
    let dispatch = &document["dispatches"][0];
    assert_eq!(transport(&dispatch["total_counters"]), 18);
    assert_eq!(transport(&dispatch["attributed_self_counters"]), 14);
    assert_eq!(transport(&dispatch["unscoped_counters"]), 4);
    assert_eq!(dispatch["attributed_self_wall_ns"], 14);
    assert_eq!(dispatch["unscoped_wall_ns"], 4);
    assert_eq!(dispatch["reconciled"], true);
}

#[test]
fn request_counters_reset_and_route_apply_identities_remain_distinct() {
    let key = IntegralKey::try_new([1, 1]).unwrap();
    let owners = [true, true];
    let mut h = Harness::new();
    h.start(0);
    h.send(
        Event::Snapshot,
        Phase::Expansion,
        1,
        0,
        Some(Work::Route(&key)),
    )
    .unwrap();
    h.send(
        Event::Snapshot,
        Phase::Expansion,
        2,
        1,
        Some(Work::Apply {
            owner_sector: &owners,
            target: &key,
        }),
    )
    .unwrap();
    h.send(
        Event::NodeCompleted(Source::Terminal),
        Phase::Expansion,
        3,
        2,
        Some(Work::Apply {
            owner_sector: &owners,
            target: &key,
        }),
    )
    .unwrap();
    h.send(
        Event::NodeCompleted(Source::Computed),
        Phase::Combination,
        4,
        3,
        Some(Work::Route(&key)),
    )
    .unwrap();
    h.send(Event::CoreFinished, Phase::RootCollection, 5, 4, None)
        .unwrap();
    h.start(10);
    h.send(
        Event::Snapshot,
        Phase::Expansion,
        11,
        0,
        Some(Work::Route(&key)),
    )
    .unwrap();
    h.send(
        Event::NodeCompleted(Source::BackingHit),
        Phase::Expansion,
        12,
        1,
        Some(Work::Route(&key)),
    )
    .unwrap();
    h.send(Event::CoreFinished, Phase::RootCollection, 13, 1, None)
        .unwrap();
    let document = h.finish();
    assert_eq!(document["complete"], true);
    assert_eq!(document["families"].as_array().unwrap().len(), 1);
    assert_eq!(document["keys"].as_array().unwrap().len(), 2);
    assert_eq!(
        document["keys"][1]["work"]["owner_sector"],
        json!([true, true])
    );
    assert_eq!(transport(&document["dispatches"][1]["total_counters"]), 1);
    assert!(document["attempts"][2]["parent_id"].is_null());
}

#[test]
fn caller_and_forwarded_saved_summaries_are_ignored() {
    let mut h = Harness::new();
    let mut ignored = Envelope::new(Event::RequestFinished, "stale-family");
    ignored.caller_phase = Some("terminal_substitution");
    assert_eq!(h.observation.observe(&ignored), Ok(Control::Continue));
    ignored.caller_phase = None;
    ignored.activity = Activity::Caller;
    assert_eq!(h.observation.observe(&ignored), Ok(Control::Continue));
    h.start(1);
    h.send(Event::CoreFinished, Phase::RootCollection, 2, 0, None)
        .unwrap();
    let document = h.finish();
    assert_eq!(document["complete"], true);
    assert_eq!(document["metadata"]["native_events"], 2);
}

#[test]
fn mismatched_completion_faults_only_the_diagnostic() {
    let first = IntegralKey::try_new([1]).unwrap();
    let other = IntegralKey::try_new([2]).unwrap();
    let mut h = Harness::new();
    h.start(0);
    h.send(
        Event::Snapshot,
        Phase::Expansion,
        1,
        0,
        Some(Work::Route(&first)),
    )
    .unwrap();
    assert_eq!(
        h.send(
            Event::NodeCompleted(Source::Computed),
            Phase::Expansion,
            2,
            1,
            Some(Work::Route(&other))
        ),
        Err(())
    );
    let document = h.finish();
    assert_eq!(document["complete"], false);
    assert_eq!(document["error"], "completion key was never entered");
    assert_eq!(document["metadata"]["uncompleted_attempts"], 1);
}

#[test]
fn truncated_or_detached_observation_never_claims_completion() {
    let key = IntegralKey::try_new([1]).unwrap();
    let mut h = Harness::new();
    h.start(0);
    h.send(
        Event::Snapshot,
        Phase::Expansion,
        1,
        0,
        Some(Work::Route(&key)),
    )
    .unwrap();
    let document = h.finish();
    assert_eq!(document["complete"], false);
    assert_eq!(
        document["error"],
        "observation ended before core completion"
    );
    let mut h = Harness::new();
    h.start(0);
    h.send(Event::CoreFinished, Phase::RootCollection, 1, 0, None)
        .unwrap();
    assert_eq!(
        h.observation.finish(true, ObserverStatus::Detached)["complete"],
        false
    );
}

#[test]
fn decreasing_and_saturated_counters_and_overflow_are_rejected() {
    let mut h = Harness::new();
    h.start(0);
    h.send(Event::ActivityStarted, Phase::Admission, 1, 2, None)
        .unwrap();
    assert_eq!(
        h.send(Event::CoreFinished, Phase::RootCollection, 2, 1, None),
        Err(())
    );
    assert_eq!(h.finish()["error"], "native counter decreased");
    let mut h = Harness::new();
    h.start(0);
    assert_eq!(
        h.send(
            Event::CoreFinished,
            Phase::RootCollection,
            1,
            usize::MAX,
            None
        ),
        Err(())
    );
    assert_eq!(h.finish()["error"], "native counter saturated");
    let mut counters = Counters([u64::MAX; COUNTER_NAMES.len()]);
    assert_eq!(
        counters.add(Counters([1; COUNTER_NAMES.len()])),
        Err("diagnostic counter overflow")
    );
}

#[test]
fn unmatched_dispatch_end_and_completion_census_are_rejected() {
    let mut h = Harness::new();
    assert_eq!(
        h.send(Event::CoreFinished, Phase::RootCollection, 1, 0, None),
        Err(())
    );
    assert_eq!(h.finish()["error"], "native event without dispatch start");
    let mut h = Harness::new();
    h.start(0);
    h.completed = 1;
    assert_eq!(
        h.send(Event::CoreFinished, Phase::RootCollection, 1, 0, None),
        Err(())
    );
    assert_eq!(
        h.finish()["error"],
        "native completion census differs from observed attempts"
    );
}

#[test]
fn changed_family_borrow_is_authenticated_without_requiring_the_same_address() {
    let mut h = Harness::new();
    h.start(0);
    let same_family = "test-family".to_owned();
    h.send_family(
        Event::ActivityStarted,
        Phase::Admission,
        1,
        0,
        None,
        &same_family,
    )
    .unwrap();
    h.send(Event::CoreFinished, Phase::RootCollection, 2, 0, None)
        .unwrap();
    assert_eq!(h.finish()["complete"], true);
    let mut h = Harness::new();
    h.start(0);
    assert_eq!(
        h.send_family(
            Event::CoreFinished,
            Phase::RootCollection,
            1,
            0,
            None,
            "other-family"
        ),
        Err(())
    );
    assert_eq!(h.finish()["error"], "family changed inside dispatch");
}

#[test]
fn final_self_sum_is_independent_and_observer_fault_status_is_preserved() {
    let key = IntegralKey::try_new([1]).unwrap();
    let mut h = Harness::new();
    h.start(0);
    h.send(
        Event::Snapshot,
        Phase::Expansion,
        1,
        0,
        Some(Work::Route(&key)),
    )
    .unwrap();
    h.send(
        Event::NodeCompleted(Source::Computed),
        Phase::Expansion,
        2,
        1,
        Some(Work::Route(&key)),
    )
    .unwrap();
    h.send(Event::CoreFinished, Phase::RootCollection, 3, 2, None)
        .unwrap();
    h.observation.attempts[0].self_counters.0[3] += 1;
    let document = h.finish();
    assert_eq!(document["complete"], false);
    assert_eq!(
        document["error"],
        "completed attempt self sum did not reconcile"
    );
    let mut h = Harness::new();
    h.start(0);
    h.send(Event::CoreFinished, Phase::RootCollection, 1, 0, None)
        .unwrap();
    let status =
        ObserverStatus::Faulted(rustred::solver::CandidateRoutedObserverFault::CallbackError);
    let document = h.observation.finish(true, status);
    assert_eq!(document["operation_succeeded"], true);
    assert_eq!(document["complete"], false);
    assert_eq!(document["observer_status"], "Faulted(CallbackError)");
}

#[test]
fn missing_native_boundary_sequence_and_failed_operation_are_incomplete() {
    let mut h = Harness::new();
    h.start(0);
    h.local_sequence += 1;
    assert_eq!(
        h.send(Event::CoreFinished, Phase::RootCollection, 1, 0, None),
        Err(())
    );
    assert_eq!(
        h.finish()["error"],
        "native snapshot sequence was not contiguous"
    );
    let mut h = Harness::new();
    h.start(0);
    h.send(Event::CoreFinished, Phase::RootCollection, 1, 0, None)
        .unwrap();
    assert_eq!(
        h.observation.finish(false, ObserverStatus::Active)["complete"],
        false
    );
}
