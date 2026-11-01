use super::*;
use rustred::solver::CandidateOrderedPhaseKind;

fn snapshot(seconds: u64, entered: u64, completed: u64, roots: u64) -> DisplaySnapshot {
    DisplaySnapshot {
        sequence: seconds + 1,
        elapsed: Duration::from_secs(seconds),
        quantitative_sample_at: Some(Duration::from_secs(seconds)),
        quantitative_sample_sequence: seconds + 1,
        phase: "native CoefficientCombination".into(),
        stage: NATIVE_STAGE.into(),
        native_work_active: true,
        unique_entered: entered,
        completed_unique: Some(completed),
        completion_attempts: Some(completed),
        completed_roots: Some(roots),
        total_roots: Some(5),
        pending_unique: Some(FrontierCount {
            count: 12,
            exact: false,
            operational_scope: true,
        }),
        ..Default::default()
    }
}

fn estimates(state: &MonitorState, now: u64) -> DfsEstimates {
    state
        .live_model(Duration::from_secs(now))
        .dfs
        .expect("depth-first census")
        .estimates
}

fn census(
    completed_phases: usize,
    processed_keys: usize,
    hits: usize,
    misses: usize,
) -> CandidateOrderedReductionProgress {
    CandidateOrderedReductionProgress {
        completed_phases,
        total_phases: 4,
        phase_kind: CandidateOrderedPhaseKind::Reduction,
        active_support: Some(completed_phases + 1),
        processed_keys,
        frontier_keys: 10 - completed_phases,
        frontier_coefficient_monomials: 100,
        step_cache_hits: hits,
        step_cache_misses: misses,
        ..Default::default()
    }
}

#[test]
fn ordered_census_shows_phase_position_grouping_and_step_cache_without_time_percent() {
    let mut state = MonitorState::default();
    let mut sample = snapshot(10, 9, 8, 1);
    sample.native_epoch = 2;
    sample.ordered = Some(OrderedReductionSample {
        counts: CandidateOrderedReductionProgress {
            completed_phases: 3,
            total_phases: 12,
            phase_kind: CandidateOrderedPhaseKind::Reduction,
            active_support: Some(6),
            owner_index: Some(2),
            input_keys: 7,
            frontier_keys: 13,
            peak_frontier_keys: 29,
            terminal_keys: 4,
            frontier_coefficient_monomials: 41,
            peak_frontier_coefficient_monomials: 127,
            terminal_coefficient_monomials: 17,
            peak_terminal_coefficient_monomials: 53,
            emitted_terms: 73,
            merged_terms: 17,
            cancelled_terms: 11,
            terminal_merges: 2,
            terminal_cancellations: 1,
            processed_keys: 19,
            step_cache_hits: 15,
            step_cache_misses: 4,
        },
        native_epoch: 2,
        family_fingerprint: "ordered-family".into(),
        sampled_at: Duration::from_secs(8),
        phase_start: PhaseStart {
            at: Duration::from_secs(4),
            processed: 15,
            hits: 12,
            misses: 3,
        },
        phases: Vec::new(),
        first: (Duration::ZERO, 0),
        window: VecDeque::from([(Duration::from_secs(4), 15)]),
    });
    state.observe(sample, Duration::from_secs(10));
    let model = state.live_model(Duration::from_secs(13));
    let ordered = model.ordered.as_ref().unwrap();
    assert_eq!(ordered.request, 2);
    assert_eq!((ordered.phase.index, ordered.phase.total), (3, 12));
    assert_eq!(ordered.phase.kind, PhaseKind::Apply);
    assert_eq!(
        (ordered.phase.support, ordered.phase.owner),
        (Some(6), Some(2))
    );
    // In phase since the census at t = 4 s; the census itself is 5 s old.
    assert!(ordered.phase.active);
    assert_eq!(ordered.phase.seconds, Some(9.0));
    assert_eq!(ordered.phase.keys, 4);
    assert_eq!(ordered.census_age_seconds, 5.0);
    // Rates include the 5 s since the census: a silent core slows down.
    assert_eq!(ordered.recent_rate.map(Rate::per_second), Some(4.0 / 9.0));
    assert_eq!(
        ordered.request_rate.map(Rate::per_second),
        Some(19.0 / 13.0)
    );
    assert_eq!(
        (ordered.frontier.keys, ordered.frontier.peak_keys),
        (13, 29)
    );
    assert_eq!(
        (ordered.frontier.monomials, ordered.frontier.peak_monomials),
        (41, 127)
    );
    assert_eq!(
        (ordered.output.keys, ordered.output.peak_monomials),
        (4, 53)
    );
    assert_eq!(
        (ordered.grouping.merged, ordered.grouping.cancelled),
        (17, 11)
    );
    assert_eq!(
        (
            ordered.grouping.output_merged,
            ordered.grouping.output_cancelled
        ),
        (2, 1)
    );
    assert_eq!(
        (ordered.step_cache_hits, ordered.step_cache_misses),
        (15, 4)
    );
    assert_eq!(
        (
            ordered.phase.step_cache_hits,
            ordered.phase.step_cache_misses
        ),
        (3, 1)
    );
    assert!(model.dfs.is_none());
    assert_eq!(model.details.family.as_deref(), Some("ordered-family"));
    let document = state.presentation(Duration::from_secs(13));
    assert_eq!(document["schema"], model::MODEL_SCHEMA);
    assert_eq!(document["ordered_reduction"]["native_request"], 2);
    assert_eq!(document["ordered_reduction"]["sample_age_seconds"], 5.0);
    assert_eq!(document["ordered_reduction"]["frontier_keys"], 13);
    assert_eq!(
        document["ordered_reduction"]["peak_frontier_coefficient_monomials"],
        127
    );
    assert_eq!(
        document["ordered_reduction"]["phase_fraction_is_time_fraction"],
        false
    );
    let cache = &document["ordered_reduction"]["step_cache"];
    assert_eq!(
        (cache["hits"].as_u64(), cache["misses"].as_u64()),
        (Some(15), Some(4))
    );
    assert_eq!(
        (cache["phase_hits"].as_u64(), cache["phase_misses"].as_u64()),
        (Some(3), Some(1))
    );
    assert_eq!(document["model"]["ordered"]["phase"]["kind"], "apply");
    assert_eq!(document["model"]["ordered"]["frontier"]["keys"], 13);
}

#[test]
fn every_phase_end_closes_one_phase_record_with_its_time_keys_and_cache() {
    let at = Duration::from_secs;
    let first = OrderedReductionSample::next(None, census(0, 0, 0, 0), 1, "family", at(1));
    let same = OrderedReductionSample::next(Some(&first), census(0, 6, 5, 1), 1, "family", at(3));
    assert_eq!(same.phase_step_cache(), (5, 1));
    assert_eq!(same.phase_start.at, at(1));
    // The census RustRed delivers at the end of phase 0 closes it exactly.
    let boundary =
        OrderedReductionSample::next(Some(&same), census(1, 9, 7, 2), 1, "family", at(4));
    assert_eq!(boundary.phase_step_cache(), (0, 0));
    assert_eq!(
        boundary.phases,
        [PhaseRecord {
            index: 0,
            kind: PhaseKind::Reduce,
            support: Some(1),
            owner: None,
            seconds: 3.0,
            keys: 9,
            step_cache_hits: 7,
            step_cache_misses: 2,
            frontier_keys: 9,
            frontier_monomials: 100,
        }]
    );
    // Empty phases are not logged; their time belongs to the next record.
    let empty =
        OrderedReductionSample::next(Some(&boundary), census(2, 9, 7, 2), 1, "family", at(5));
    assert_eq!(empty.phases.len(), 1);
    assert_eq!(empty.phase_start.at, at(5));
    let busy = OrderedReductionSample::next(Some(&empty), census(2, 29, 7, 2), 1, "family", at(15));
    let (recent, request) = busy.rates(0.0);
    assert_eq!(recent.map(Rate::per_second), Some(29.0 / 14.0));
    assert_eq!(request.map(Rate::per_second), Some(29.0 / 14.0));
    let done = OrderedReductionSample::next(Some(&busy), census(4, 29, 7, 2), 1, "family", at(16));
    assert_eq!(done.phases.len(), 2);
    assert_eq!((done.phases[1].index, done.phases[1].keys), (2, 20));
    assert_eq!(done.phases[1].seconds, 11.0);
    let model = done.model(None, false, 0.0);
    assert_eq!(model.phase.kind, PhaseKind::Done);
    // A new native request restarts the log, the phase and the rates.
    let next = OrderedReductionSample::next(Some(&done), census(0, 1, 1, 0), 2, "family", at(17));
    assert!(next.phases.is_empty());
    assert_eq!(next.phase_step_cache(), (0, 0));
    assert_eq!(next.rates(0.0), (None, None));
    let (recent, request) = busy.rates(14.0);
    assert_eq!(recent.map(Rate::per_second), Some(29.0 / 28.0));
    assert_eq!(request.map(Rate::per_second), Some(29.0 / 28.0));
}

#[test]
fn a_phase_interrupted_by_a_failure_keeps_its_time_and_claims_no_current_rate() {
    let at = Duration::from_secs;
    let mut state = MonitorState::default();
    // Phase 1 opens at t = 10 s; the core works on it until it fails at
    // t = 112 s, whose final census RustRed reports with the failure.
    let start = OrderedReductionSample::next(None, census(0, 0, 0, 0), 1, "family", at(0));
    let opened =
        OrderedReductionSample::next(Some(&start), census(1, 10, 0, 0), 1, "family", at(10));
    let working =
        OrderedReductionSample::next(Some(&opened), census(1, 1010, 0, 0), 1, "family", at(110));
    let mut native = snapshot(110, 0, 0, 0);
    native.native_epoch = 1;
    native.ordered = Some(working.clone());
    state.observe(native.clone(), at(110));
    let running = state.live_model(at(111));
    let phase = &running.ordered.as_ref().unwrap().phase;
    assert!(phase.active);
    assert_eq!(phase.seconds, Some(101.0));
    let mut failed = native;
    failed.sequence += 1;
    failed.elapsed = at(112);
    failed.native_work_active = false;
    failed.status = MonitorStatus::Failed;
    failed.request_finished = true;
    failed.ordered = Some(OrderedReductionSample::next(
        Some(&working),
        census(1, 1010, 0, 0),
        1,
        "family",
        at(112),
    ));
    state.observe(failed, at(112));
    let model = state.live_model(at(400));
    let ordered = model.ordered.as_ref().unwrap();
    assert_eq!(model.status, MonitorStatus::Failed);
    assert!(!ordered.phase.active);
    // The phase ran from its opening census to the failure, not 0 s.
    assert_eq!(ordered.phase.seconds, Some(102.0));
    assert_eq!(ordered.phase.keys, 1000);
    assert!(ordered.recent_rate.is_none());
    assert_eq!(
        ordered.request_rate,
        Some(Rate {
            keys: 1010,
            seconds: 112.0
        })
    );
    let document = state.presentation(at(400));
    assert_eq!(document["model"]["ordered"]["phase"]["seconds"], 102.0);
    assert_eq!(document["model"]["ordered"]["phase"]["active"], false);
    #[cfg(feature = "progress-tui")]
    {
        use super::tui::tests::{render, text};
        let screen = text(&render(&model, UiView::default(), 80, 24, true));
        let row = screen
            .lines()
            .find(|line| line.contains("Reduce s2") && !line.contains("Phase"))
            .unwrap_or_else(|| panic!("no row for the interrupted phase:\n{screen}"));
        assert!(row.contains("1:42") && row.contains("9.8/s") && row.contains("idle"));
        assert!(!row.contains("0.0s"), "{screen}");
        assert!(screen.contains("1:42 in phase") && screen.contains("FAILED"));
        assert!(!screen.contains(" now"), "{screen}");
    }
}

#[test]
fn step_cache_phase_entries_keep_the_native_kind_and_their_cache_activity() {
    let at = Duration::from_secs;
    let routing = |completed, processed, hits, misses| CandidateOrderedReductionProgress {
        phase_kind: CandidateOrderedPhaseKind::Routing,
        ..census(completed, processed, hits, misses)
    };
    let mut state = MonitorState::default();
    let mut native = snapshot(10, 0, 0, 0);
    let first = OrderedReductionSample::next(None, routing(0, 0, 0, 0), 1, "family", at(1));
    // Phase 0 (Route) replays 4 steps; phase 1 (Reduce) has no cache activity.
    let routed = OrderedReductionSample::next(Some(&first), census(1, 5, 4, 1), 1, "family", at(4));
    let reduced =
        OrderedReductionSample::next(Some(&routed), census(2, 9, 4, 1), 1, "family", at(6));
    assert_eq!(reduced.phases.len(), 2);
    native.ordered = Some(reduced);
    state.observe(native, at(10));
    let document = state.presentation(at(10));
    let phases = document["ordered_reduction"]["step_cache"]["completed_phases"]
        .as_array()
        .unwrap();
    assert_eq!(phases.len(), 1, "{phases:?}");
    assert_eq!(phases[0]["phase"], 0);
    assert_eq!(phases[0]["phase_kind"], "Routing");
    assert_eq!(phases[0]["kind"], "route");
    assert_eq!(
        (phases[0]["hits"].as_u64(), phases[0]["misses"].as_u64()),
        (Some(4), Some(1))
    );
    assert_eq!(document["ordered_reduction"]["phase_kind"], "Reduction");
}

#[test]
fn the_phase_log_and_rate_window_stay_bounded() {
    let mut sample = OrderedReductionSample::next(
        None,
        CandidateOrderedReductionProgress {
            total_phases: 1000,
            ..Default::default()
        },
        1,
        "family",
        Duration::ZERO,
    );
    for step in 1..=500_usize {
        let counts = CandidateOrderedReductionProgress {
            completed_phases: step / 10,
            total_phases: 1000,
            processed_keys: step * 3,
            ..Default::default()
        };
        sample = OrderedReductionSample::next(
            Some(&sample),
            counts,
            1,
            "family",
            Duration::from_millis(250 * step as u64),
        );
    }
    assert_eq!(sample.phases.len(), PHASE_LOG);
    assert!(sample.window.len() <= 62, "{}", sample.window.len());
    let (recent, _) = sample.rates(0.0);
    assert_eq!(recent.map(Rate::per_second), Some(12.0));
}

#[test]
fn batch_and_cohort_counters_stay_live_while_native_panels_are_frozen() {
    let mut state = MonitorState::default();
    let mut initial = snapshot(10, 8, 5, 0);
    initial.ordered_batches = Some(CallerWork {
        phase: "1L / Retain".into(),
        completed: 0,
        total: 2,
        unit: "ordered batches completed".into(),
        sampled_at: Duration::from_secs(10),
    });
    initial.reduction_cohorts = Some(CallerWork {
        phase: "current collected ordered batch".into(),
        completed: 1,
        total: 2,
        unit: "compatible scalar cohorts completed".into(),
        sampled_at: Duration::from_secs(10),
    });
    state.observe(initial, Duration::from_secs(10));
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(10));
    state.latest.ordered_batches.as_mut().unwrap().completed = 1;
    state.latest.ordered_batches.as_mut().unwrap().phase = "2L / Retain".into();
    state.latest.reduction_cohorts.as_mut().unwrap().total = 3;
    state.latest.completed_unique = Some(7);
    let model = state.model(Duration::from_secs(15));
    assert!(state.ui().frozen);
    let batch = model.batch.as_ref().unwrap();
    assert_eq!((batch.completed, batch.total), (1, 2));
    assert_eq!(batch.label, "2L / Retain");
    let cohort = model.cohort.as_ref().unwrap();
    assert_eq!((cohort.completed, cohort.total), (1, 3));
    // The native census is the frozen one.
    assert_eq!(model.dfs.as_ref().unwrap().completed, Some(5));
    let decoded: serde_json::Value = serde_json::from_str(
        &serde_json::to_string(&state.presentation(Duration::from_secs(15))).unwrap(),
    )
    .unwrap();
    assert_eq!(decoded["ordered_batches"]["completed"], 1);
    assert_eq!(decoded["ordered_batches"]["label"], "2L / Retain");
    assert_eq!(decoded["reduction_cohorts"]["total"], 3);
    assert_eq!(
        decoded["reduction_cohorts"]["fraction_is_time_fraction"],
        false
    );
    // The published model is live, never the frozen display.
    let published: MonitorModel = serde_json::from_value(decoded["model"].clone()).unwrap();
    assert_eq!(published.dfs.unwrap().completed, Some(7));
    assert_eq!(published.cohort.unwrap().total, 3);
    state.latest.reduction_cohorts = None;
    assert!(state.model(Duration::from_secs(16)).cohort.is_none());
}

#[test]
fn proc_memory_preserves_process_peak_and_host_usable_unused_scopes() {
    let resources = ResourceUsage::from_proc(
        "Name:\tfixture\nVmRSS:\t17 kB\nVmHWM:\t29 kB\n",
        "MemTotal: 1000 kB\nMemAvailable: 700 kB\nMemFree: 200 kB\n",
    );
    assert_eq!(resources.rss_bytes, Some(17 * 1024));
    assert_eq!(resources.peak_rss_bytes, Some(29 * 1024));
    assert_eq!(resources.mem_total_bytes, Some(1000 * 1024));
    assert_eq!(resources.mem_available_bytes, Some(700 * 1024));
    assert_eq!(resources.mem_free_bytes, Some(200 * 1024));
    assert_ne!(resources.mem_available_bytes, resources.mem_free_bytes);
    assert!(resources.sample_age_seconds.is_none());
}

#[test]
fn readable_memory_uses_iec_units_without_changing_exact_raw_observations() {
    for (bytes, expected) in [
        (0, "0 B"),
        (1023, "1023 B"),
        (1024, "1.0 KiB"),
        (1024 * 1024, "1.0 MiB"),
        (50_000_000_000, "46.6 GiB"),
        (u64::MAX, "16.0 EiB"),
    ] {
        assert_eq!(model::format::bytes(bytes), expected);
    }
    let mut state = MonitorState::default();
    state.latest.resources = ResourceUsage {
        process_identity_verified: true,
        rss_bytes: Some(50_000_000_000),
        peak_rss_bytes: Some(64 * 1024 * 1024 * 1024),
        mem_total_bytes: Some(u64::MAX),
        mem_available_bytes: Some(1024 * 1024 * 1024),
        mem_free_bytes: Some(0),
        guard: ResourceGuardView {
            threshold_bytes: Some(50_000_000_000),
            ..Default::default()
        },
        ..Default::default()
    };
    let published = state.presentation(Duration::ZERO);
    let memory = state.live_model(Duration::ZERO).memory;
    assert_eq!(memory.rss_bytes, Some(50_000_000_000));
    assert_eq!(memory.total_bytes, Some(u64::MAX));
    let guard = memory.guard.unwrap();
    assert_eq!(guard.status, "unknown");
    assert_eq!(guard.threshold_bytes, Some(50_000_000_000));
    assert!(!guard.covered);
    assert_eq!(published["resources"]["rss_bytes"], 50_000_000_000_u64);
    assert_eq!(published["resources"]["mem_free_bytes"], 0);
    assert_eq!(published["resources"]["mem_total_bytes"], u64::MAX);
    assert_eq!(
        published["resources"]["guard"]["threshold_bytes"],
        50_000_000_000_u64
    );
    assert_eq!(published["resources"]["guard"]["status"], "unknown");
    // Process values of an unverified worker are never shown.
    state.latest.resources.process_identity_verified = false;
    let memory = state.live_model(Duration::ZERO).memory;
    assert!(memory.rss_bytes.is_none() && memory.peak_rss_bytes.is_none());
    assert_eq!(memory.available_bytes, Some(1024 * 1024 * 1024));
}

#[test]
fn proc_memory_missing_malformed_units_duplicates_and_overflow_stay_unknown() {
    assert!(ResourceUsage::from_proc("", "").source.is_none());
    for status in [
        "",
        "VmRSS: missing kB",
        "VmRSS: -1 kB",
        "VmRSS: 12 MB",
        "VmRSS: 12 kB extra",
        "VmRSS: 1 kB\nVmRSS: 2 kB",
        "VmRSS: 18446744073709551615 kB",
    ] {
        let resources = ResourceUsage::from_proc(status, "MemAvailable: 9 kB\n");
        assert_eq!(resources.rss_bytes, None, "{status}");
        assert_eq!(resources.mem_available_bytes, Some(9 * 1024));
        assert_eq!(resources.mem_free_bytes, None);
        assert_eq!(resources.peak_rss_bytes, None);
    }
}

#[test]
fn display_freeze_keeps_resources_and_status_live_in_model_and_json() {
    let mut state = MonitorState::default();
    state.observe(snapshot(0, 3, 1, 0), Duration::ZERO);
    state.latest.resources = ResourceUsage::from_proc("VmRSS: 1 kB", "MemAvailable: 2 kB");
    state.latest.resource_sample_at = Some(Duration::ZERO);
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(1));
    state.latest.resources = ResourceUsage::from_proc("VmRSS: 3 kB", "MemAvailable: 4 kB");
    state.latest.resources.process_identity_verified = true;
    state.latest.resources.guard = ResourceGuardView {
        status: "latched_abort".into(),
        sample_seconds: Some(1),
        report_age_seconds: Some(0.0),
        reported_latched_abort: Some(true),
        ..Default::default()
    };
    state.latest.resource_sample_at = Some(Duration::from_secs(5));
    // A synthetic terminal failure after native silence retains the old census.
    let mut failed = state.latest.clone();
    failed.elapsed = Duration::from_secs(6);
    failed.status = MonitorStatus::Failed;
    failed.request_finished = true;
    failed.phase = "synthetic_failed_finalization".into();
    failed
        .events
        .push("cycle detected · synthetic fixture".into());
    state.observe(failed, Duration::from_secs(6));
    let frozen = state.model(Duration::from_secs(7));
    assert!(state.ui().frozen);
    assert_eq!(frozen.status, MonitorStatus::Failed);
    assert_eq!(frozen.memory.rss_bytes, Some(3 * 1024));
    assert_eq!(
        frozen.memory.guard.as_ref().unwrap().status,
        "latched_abort"
    );
    assert_eq!(frozen.events.len(), 1);
    let live = serde_json::to_value(state.live_resources(Duration::from_secs(7))).unwrap();
    assert_eq!(live["rss_bytes"], 3 * 1024);
    assert_eq!(live["mem_available_bytes"], 4 * 1024);
    assert_eq!(live["sample_age_seconds"], 2.0);
    assert_eq!(live["guard"]["status"], "latched_abort");
    assert_eq!(live["guard"]["report_age_seconds"], 2.0);
    let stale = state.live_resources(Duration::from_secs(9));
    assert_eq!(stale.guard.status, "unknown");
    assert_eq!(stale.guard.reported_latched_abort, Some(true));
    let published = state.presentation(Duration::from_secs(7));
    assert_eq!(published["status"], "Failed");
    assert_eq!(published["phase"], "synthetic_failed_finalization");
    assert!(published.get("view").is_none());
    assert_eq!(published["model"]["status"], "Failed");
    // A finished run has no live estimate.
    assert_eq!(published["model"]["dfs"]["estimates"]["state"], "inactive");
    assert_eq!(published["resources"]["guard"]["status"], "latched_abort");
}

#[test]
fn fresh_os_memory_does_not_refresh_native_progress_eta_or_infer_cycle() {
    let mut state = MonitorState::default();
    let mut native = snapshot(0, 3, 1, 0);
    native.native_observed_at = Some(Duration::ZERO);
    native.current_key = Some("Route [2,1]".into());
    native.current_degrees = Some((0, 3, 3));
    state.observe(native, Duration::ZERO);
    state.latest.resources = ResourceUsage::from_proc(
        "VmRSS: 10 kB\nVmHWM: 20 kB",
        "MemTotal: 100 kB\nMemAvailable: 80 kB\nMemFree: 30 kB",
    );
    state.latest.resource_sample_at = Some(Duration::from_secs(5));
    let model = state.live_model(Duration::from_secs(6));
    assert_eq!(model.memory.sample_age_seconds, Some(1.0));
    assert_eq!(model.details.native_observation_age_seconds, Some(6.0));
    let estimates = model.dfs.as_ref().unwrap().estimates.clone();
    assert_eq!(estimates.sample_age_seconds, Some(6.0));
    assert!(estimates.frontier_eta_seconds.is_none());
    assert_eq!(
        model.details.current_degrees,
        Some(["0".into(), "3".into(), "3".into()])
    );
    assert!(model.details.caller_expression.is_none());
    assert_eq!(model.status, MonitorStatus::Running);
    assert!(model.events.is_empty());
    let json = serde_json::to_value(state.live_resources(Duration::from_secs(6))).unwrap();
    assert_eq!(json["mem_available_bytes"], 80 * 1024);
    assert_eq!(json["mem_free_bytes"], 30 * 1024);
    assert_eq!(json["sample_age_seconds"], 1.0);
}

#[test]
fn stored_symbolic_partitions_and_caller_sizes_keep_independent_scopes_and_ages() {
    let mut state = MonitorState::default();
    let mut observed = snapshot(10, 8, 3, 1);
    observed.stored_state = Some(StoredStateSample {
        counts: CandidateRoutedStoredState {
            unstarted_roots: Some(2),
            deferred_rule_maps: Some(4),
            deferred_child_integral_occurrences: Some(11),
            collected_output_integrals: Some(3),
        },
        family_fingerprint: "native-family".into(),
        sampled_at: Duration::from_secs(8),
    });
    observed.caller_expression = Some(CallerExpressionSample {
        encoded_bytes: 937,
        stored_summands: Some(7),
        phase: "scalar_input".into(),
        sampled_at: Duration::from_secs(2),
    });
    // Reception/rendering and native samples intentionally use different epochs.
    state.observe(observed, Duration::from_secs(100));
    let model = state.live_model(Duration::from_secs(103));
    let expression = model.details.caller_expression.unwrap();
    assert_eq!(
        (
            expression.bytes,
            expression.summands,
            expression.stage.as_str()
        ),
        (937, Some(7), "scalar_input")
    );
    assert_eq!(model.details.family.as_deref(), Some("native-family"));
    let published = state.presentation(Duration::from_secs(103));
    assert_eq!(
        published["stored_symbolic_state"]["sample_age_seconds"],
        5.0
    );
    assert_eq!(
        published["stored_symbolic_state"]["deferred_child_integral_occurrences"],
        11
    );
    assert_eq!(published["caller_expression"]["sample_age_seconds"], 11.0);
    assert_eq!(published["caller_expression"]["encoded_bytes"], 937);
}

#[test]
fn default_is_unobserved_but_snapshot_path_is_an_explicit_opt_in() {
    let mut options = RustRedMonitorOptions::default();
    assert!(!options.enabled());
    options.validate().unwrap();
    options.progress_snapshot = Some(PathBuf::from("call-progress.jsonl"));
    assert!(options.enabled());
    options.validate().unwrap();
    options.display = RustRedMonitorDisplay::Json;
    options.progress_snapshot = None;
    assert!(options.validate().is_err());
}

#[test]
fn writer_preserves_the_native_projection_without_a_consumer_schema() {
    let native = serde_json::json!({
        "schema": "rustred.saved-application-progress.v1",
        "closure_claim": false,
        "core": {"pinned": null, "pending_frames": 19}
    });
    let mut writer = SnapshotWriter::new(Vec::new());
    writer.write(&native).unwrap();
    assert_eq!(writer.0.last(), Some(&b'\n'));
    assert_eq!(
        serde_json::from_slice::<serde_json::Value>(&writer.0).unwrap(),
        native
    );
}

#[cfg(all(unix, not(target_arch = "wasm32")))]
mod atomic_file {
    use super::*;
    use std::{cell::Cell, collections::BTreeSet, fs, path::Path};

    fn entries(path: &Path) -> BTreeSet<PathBuf> {
        fs::read_dir(path)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .collect()
    }

    fn sink(path: &Path) -> LatestJsonFile {
        LatestJsonFile::from_options(&RustRedMonitorOptions {
            progress_snapshot: Some(path.to_owned()),
            ..Default::default()
        })
    }

    #[test]
    fn default_off_does_not_create_a_temporary_or_call_the_writer() {
        let writes = Cell::new(0);
        let file = LatestJsonFile::from_options(&RustRedMonitorOptions::default());
        file.write_with(|_| {
            writes.set(writes.get() + 1);
            Err(io::Error::other("disabled output must not reach I/O"))
        })
        .unwrap();
        assert_eq!(writes.get(), 0);
    }

    #[test]
    fn path_alone_publishes_one_complete_latest_json_and_leaves_no_scratch() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("latest.json");
        let file = sink(&path);
        assert!(entries(directory.path()).is_empty());
        for sequence in [1, 2] {
            let native = serde_json::json!({
                "schema": "rustred.saved-application-progress.v1",
                "sequence": sequence,
                "unknown_future_field": [null, {"scope": "native cache payload"}],
                "text": "x".repeat(16384)
            });
            file.write(&native).unwrap();
            let actual = fs::read(&path).unwrap();
            assert_eq!(
                serde_json::from_slice::<serde_json::Value>(&actual).unwrap(),
                native
            );
            assert_eq!(entries(directory.path()), BTreeSet::from([path.clone()]));
        }
    }

    #[test]
    fn staged_partial_json_does_not_replace_the_old_snapshot() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("latest.json");
        let previous = b"{\"sequence\":1}\n";
        fs::write(&path, previous).unwrap();
        sink(&path)
            .write_with(|file| {
                file.write_all(b"{\"sequence\":")?;
                assert_eq!(fs::read(&path)?, previous);
                file.write_all(b"2}\n")
            })
            .unwrap();
        assert_eq!(fs::read(&path).unwrap(), b"{\"sequence\":2}\n");
        assert_eq!(entries(directory.path()), BTreeSet::from([path]));
    }

    struct FailAfter<W> {
        writer: W,
        remaining: usize,
    }

    impl<W: Write> Write for FailAfter<W> {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            if self.remaining == 0 {
                return Err(io::Error::other("intentional serializer output failure"));
            }
            let count = bytes.len().min(self.remaining);
            let written = self.writer.write(&bytes[..count])?;
            self.remaining -= written;
            Ok(written)
        }

        fn flush(&mut self) -> io::Result<()> {
            self.writer.flush()
        }
    }

    #[test]
    fn serializer_write_failure_preserves_old_bytes_and_removes_scratch() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("latest.json");
        let previous = b"{\"sequence\":1}\n";
        fs::write(&path, previous).unwrap();
        let error = sink(&path)
            .write_with(|file| {
                serde_json::to_writer(
                    FailAfter {
                        writer: file,
                        remaining: 3,
                    },
                    &serde_json::json!({"sequence": 2, "payload": [1, 2, 3]}),
                )
                .map_err(io::Error::from)
            })
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("intentional serializer output failure")
        );
        assert_eq!(fs::read(&path).unwrap(), previous);
        assert_eq!(entries(directory.path()), BTreeSet::from([path]));
    }

    #[test]
    fn real_persist_failure_preserves_destination_and_removes_scratch() {
        let directory = tempfile::tempdir().unwrap();
        let destination = directory.path().join("existing-directory");
        fs::create_dir(&destination).unwrap();
        let sentinel = destination.join("sentinel");
        fs::write(&sentinel, b"original").unwrap();
        sink(&destination)
            .write(&serde_json::json!({"sequence": 2}))
            .unwrap_err();
        assert_eq!(fs::read(&sentinel).unwrap(), b"original");
        assert_eq!(entries(directory.path()), BTreeSet::from([destination]));
    }

    #[test]
    fn missing_parent_is_not_created_by_snapshot_configuration_or_failed_write() {
        let directory = tempfile::tempdir().unwrap();
        let missing = directory.path().join("missing");
        let file = sink(&missing.join("latest.json"));
        assert!(!missing.exists());
        file.write(&serde_json::json!({"sequence": 1})).unwrap_err();
        assert!(!missing.exists());
        assert!(entries(directory.path()).is_empty());
    }
}

#[test]
fn observational_etas_report_method_scope_growth_and_age() {
    let mut state = MonitorState::default();
    state.observe(snapshot(0, 10, 0, 0), Duration::ZERO);
    state.observe(snapshot(10, 20, 5, 1), Duration::from_secs(10));
    let estimates = estimates(&state, 13);
    assert_eq!(estimates.state, EstimateState::Live);
    assert_eq!(estimates.frontier_eta_seconds, Some(24.0));
    assert_eq!(estimates.completion_rate, Some(0.5));
    assert_eq!(estimates.root_eta_seconds, Some(40.0));
    assert_eq!(estimates.growth_rate, Some(0.5));
    assert_eq!(estimates.sample_age_seconds, Some(3.0));
    assert_eq!(
        (estimates.window_seconds, estimates.samples),
        (Some(10.0), 2)
    );
}

#[test]
fn unknown_or_incomparable_frontier_never_becomes_a_node_eta() {
    let mut state = MonitorState::default();
    state.observe(snapshot(0, 10, 0, 0), Duration::ZERO);
    let mut next = snapshot(10, 20, 5, 0);
    next.pending_unique.as_mut().unwrap().operational_scope = false;
    state.observe(next, Duration::from_secs(10));
    let estimates = estimates(&state, 10);
    assert!(estimates.frontier_eta_seconds.is_none());
    assert!(estimates.root_eta_seconds.is_none());
    let mut next = snapshot(20, 20, 5, 0);
    next.pending_unique = None;
    state.observe(next, Duration::from_secs(20));
    let model = state.live_model(Duration::from_secs(20));
    assert!(model.dfs.unwrap().pending_unique.is_none());
}

#[test]
fn pause_and_detach_only_change_the_display() {
    let mut state = MonitorState::default();
    state.observe(snapshot(1, 10, 1, 0), Duration::from_secs(1));
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(1));
    let frozen = state.model(Duration::from_secs(1)).dfs;
    state.observe(snapshot(2, 50, 20, 1), Duration::from_secs(2));
    assert_eq!(state.model(Duration::from_secs(2)).dfs, frozen);
    assert_eq!(state.latest.unique_entered, 50);
    let now = Duration::from_secs(2);
    let live = serde_json::to_value(state.live_model(now)).unwrap();
    for (action, scroll, help, details) in [
        (UiAction::ScrollUp, 0, false, false),
        (UiAction::PageUp, 0, false, false),
        (UiAction::ToggleHelp, 0, true, false),
        (UiAction::ScrollDown, 1, true, false),
        (UiAction::PageDown, 11, true, false),
        (UiAction::ToggleDetails, 0, true, true),
        (UiAction::PageDown, 10, true, true),
        (UiAction::ScrollUp, 9, true, true),
        (UiAction::ToggleHelp, 9, false, true),
    ] {
        state.apply(action, now);
        let ui = state.ui();
        assert_eq!(
            (ui.scroll, ui.show_help, ui.show_details),
            (scroll, help, details)
        );
        assert_eq!(state.model(now).dfs, frozen);
        assert_eq!(serde_json::to_value(state.live_model(now)).unwrap(), live);
    }
    state.clamp_scroll(3);
    assert_eq!(state.ui().scroll, 3);
    state.terminal.ui.scroll = usize::MAX - 1;
    for action in [UiAction::PageDown, UiAction::ScrollDown] {
        state.apply(action, now);
        assert_eq!(state.ui().scroll, usize::MAX);
    }
    state.apply(UiAction::Detach, Duration::from_secs(2));
    assert!(state.detached() && !state.interrupted());
    assert_eq!(state.latest.status, MonitorStatus::Running);
    assert_eq!(state.latest.completed_unique, Some(20));
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(2));
    assert!(!state.ui().frozen);
    state.apply(UiAction::Interrupt, Duration::from_secs(2));
    assert!(state.detached() && state.interrupted());
}

#[test]
fn core_completion_is_not_final_request_success_or_mathematical_closure() {
    let mut state = MonitorState::default();
    let mut finished_core = snapshot(1, 10, 10, 1);
    finished_core.status = MonitorStatus::Succeeded;
    finished_core.request_finished = false;
    state.observe(finished_core, Duration::from_secs(1));
    let model = state.live_model(Duration::from_secs(1));
    assert_eq!(model.status, MonitorStatus::Running);
    assert_eq!(
        state.presentation(Duration::from_secs(1))["status"],
        "Running"
    );
}

#[test]
fn silent_work_advances_observation_age_and_completion_age_without_inventing_work() {
    let mut state = MonitorState::default();
    state.observe(snapshot(0, 10, 0, 0), Duration::ZERO);
    state.observe(snapshot(10, 20, 5, 1), Duration::from_secs(10));
    let model = state.live_model(Duration::from_secs(30));
    assert_eq!(model.details.event_age_seconds, Some(20.0));
    let dfs = model.dfs.unwrap();
    assert_eq!(dfs.since_completion_seconds, Some(20.0));
    assert_eq!(dfs.completed, Some(5));
    assert_eq!(dfs.entered, 20);
}

#[test]
fn request_reset_discards_previous_rate_history() {
    let mut state = MonitorState::default();
    state.observe(snapshot(10, 10, 5, 1), Duration::from_secs(10));
    state.observe(snapshot(20, 20, 15, 2), Duration::from_secs(20));
    state.observe(snapshot(0, 1, 0, 0), Duration::from_secs(30));
    assert_eq!(state.samples.len(), 1);
    assert!(state.last_completion.is_none());
    let estimates = estimates(&state, 30);
    assert!(estimates.root_eta_seconds.is_none());
    assert_eq!(estimates.state, EstimateState::Warmup);
}

#[test]
fn carried_samples_do_not_refresh_eta_or_evict_the_rate_window() {
    let mut state = MonitorState::default();
    state.observe(snapshot(0, 10, 0, 0), Duration::ZERO);
    let quantitative = snapshot(10, 20, 5, 1);
    state.observe(quantitative.clone(), Duration::from_secs(10));
    // Fake caller/activity frames carry the same accepted quantitative marker.
    for seconds in 11..=90 {
        let mut carried = quantitative.clone();
        carried.sequence = seconds + 1;
        carried.elapsed = Duration::from_secs(seconds);
        state.observe(carried, Duration::from_secs(seconds));
    }
    let model = state.live_model(Duration::from_secs(93));
    assert_eq!(state.samples.len(), 2);
    assert_eq!(model.details.event_age_seconds, Some(3.0));
    let dfs = model.dfs.unwrap();
    assert_eq!(dfs.since_completion_seconds, Some(83.0));
    let estimates = dfs.estimates;
    assert_eq!(estimates.state, EstimateState::Stale);
    assert_eq!(estimates.sample_age_seconds, Some(83.0));
    assert_eq!(estimates.window_seconds, Some(10.0));
    assert!(estimates.frontier_eta_seconds.is_none() && estimates.root_eta_seconds.is_none());
    assert!(estimates.growth_rate.is_none());
}

#[test]
fn genuine_unchanged_sample_extends_the_observed_rate_window() {
    let mut state = MonitorState::default();
    state.observe(snapshot(0, 10, 0, 0), Duration::ZERO);
    state.observe(snapshot(10, 20, 5, 1), Duration::from_secs(10));
    state.observe(snapshot(30, 20, 5, 1), Duration::from_secs(30));
    let model = state.live_model(Duration::from_secs(30));
    assert_eq!(state.samples.len(), 3);
    let dfs = model.dfs.unwrap();
    assert_eq!(dfs.estimates.sample_age_seconds, Some(0.0));
    assert_eq!(dfs.estimates.frontier_eta_seconds, Some(72.0));
    assert_eq!(dfs.estimates.root_eta_seconds, Some(120.0));
    assert_eq!(dfs.since_completion_seconds, Some(20.0));
}

#[test]
fn unobserved_quantitative_data_never_creates_a_rate_sample() {
    let mut state = MonitorState::default();
    let mut notice = snapshot(20, 20, 5, 1);
    notice.quantitative_sample_at = None;
    notice.quantitative_sample_sequence = 0;
    state.observe(notice, Duration::from_secs(20));
    assert!(state.samples.is_empty());
    let estimates = state.dfs_estimates(Duration::from_secs(5));
    assert_eq!(estimates.state, EstimateState::Warmup);
    assert!(estimates.sample_age_seconds.is_none());
    // Without a census there is no depth-first panel at all.
    assert!(state.live_model(Duration::from_secs(25)).dfs.is_none());
}

#[test]
fn new_epoch_with_higher_counters_does_not_invent_a_completion() {
    let mut state = MonitorState::default();
    state.observe(snapshot(0, 10, 0, 0), Duration::ZERO);
    state.observe(snapshot(10, 20, 5, 1), Duration::from_secs(10));
    let mut next = snapshot(30, 100, 80, 4);
    next.native_epoch = 1;
    state.observe(next, Duration::from_secs(30));
    assert_eq!(state.samples.len(), 1);
    assert!(state.last_completion.is_none());
    assert!(estimates(&state, 30).frontier_eta_seconds.is_none());
}

#[test]
fn current_canonical_degrees_are_shown_only_for_the_current_key() {
    let mut data = snapshot(1, 10, 2, 0);
    data.current_key = Some("Route [-1,2,8]".into());
    data.current_degrees = Some((1, 10, 9));
    let mut state = MonitorState::default();
    state.observe(data.clone(), Duration::from_secs(1));
    let details = state.live_model(Duration::from_secs(1)).details;
    assert_eq!(details.current_key.as_deref(), Some("Route [-1,2,8]"));
    assert_eq!(
        details.current_degrees,
        Some(["1".into(), "10".into(), "9".into()])
    );
    data.current_key = None;
    data.current_degrees = None;
    data.elapsed = Duration::from_secs(2);
    state.observe(data, Duration::from_secs(2));
    let details = state.live_model(Duration::from_secs(2)).details;
    assert!(details.current_key.is_none() && details.current_degrees.is_none());
}

#[test]
fn stage_clock_and_resources_remain_live_while_native_panels_are_frozen() {
    let mut state = MonitorState::default();
    let mut event = snapshot(10, 7, 3, 1);
    event.phase = "retained_catalogue_exact_combination".into();
    event.stage = "retained_catalogue_exact_combination".into();
    event.native_work_active = false;
    event.phase_started_at = Duration::from_secs(10);
    event.stage_started_at = Duration::from_secs(10);
    state.observe(event, Duration::from_secs(20));
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(20));
    state.latest.resources.process_identity_verified = true;
    state.latest.resources.rss_bytes = Some(12345);
    state.latest.resource_sample_at = Some(Duration::from_secs(56));
    let model = state.model(Duration::from_secs(56));
    assert!(state.ui().frozen);
    assert_eq!(model.stage.kind, StageKind::PostProcessing);
    assert_eq!(model.stage.seconds, 36.0);
    assert_eq!(model.details.event_age_seconds, Some(36.0));
    assert_eq!(model.memory.rss_bytes, Some(12345));
    let document = state.presentation(Duration::from_secs(56));
    assert!(document["caller_work"].is_null());
    assert_eq!(document["phase_elapsed_seconds"], 36.0);
    assert_eq!(document["event_age_seconds"], 36.0);
    assert_eq!(state.latest.completed_unique, Some(3));
}

#[test]
fn completed_native_frontier_never_finishes_or_zeroes_the_public_postprocessing_eta() {
    let mut state = MonitorState::default();
    let mut first = snapshot(0, 7, 0, 0);
    first.total_roots = Some(1);
    state.observe(first, Duration::from_secs(90));
    let mut core = snapshot(10, 7, 7, 1);
    core.total_roots = Some(1);
    core.pending_unique.as_mut().unwrap().count = 0;
    core.native_observed_at = Some(Duration::from_secs(10));
    core.native_work_active = false;
    core.status = MonitorStatus::Succeeded;
    state.observe(core, Duration::from_secs(100));
    let just_finished = state.live_model(Duration::from_secs(100));
    assert_eq!(just_finished.status, MonitorStatus::Running);
    let estimates = &just_finished.dfs.as_ref().unwrap().estimates;
    assert_eq!(estimates.state, EstimateState::Inactive);
    assert!(estimates.frontier_eta_seconds.is_none() && estimates.root_eta_seconds.is_none());
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(100));
    let mut caller = state.latest.clone();
    caller.phase = "retained_coefficient_materialization".into();
    caller.stage = "retained_coefficient_materialization".into();
    caller.phase_started_at = Duration::from_secs(11);
    caller.stage_started_at = Duration::from_secs(11);
    caller.elapsed = Duration::from_secs(11);
    caller.sequence += 1;
    caller.status = MonitorStatus::Running;
    caller.caller_work = Some(CallerWork {
        phase: caller.phase.clone(),
        completed: 1,
        total: 58,
        unit: "master products materialized".into(),
        sampled_at: Duration::from_secs(101),
    });
    state.observe(caller, Duration::from_secs(101));
    let model = state.model(Duration::from_secs(137));
    assert_eq!(model.status, MonitorStatus::Running);
    assert_eq!(model.stage.kind, StageKind::PostProcessing);
    assert_eq!(model.stage.id, "retained_coefficient_materialization");
    assert_eq!(model.stage.seconds, 36.0);
    let work = model.stage.work.as_ref().unwrap();
    assert_eq!((work.completed, work.total), (1, 58));
    assert_eq!(work.label, "master products materialized");
    let dfs = model.dfs.as_ref().unwrap();
    assert_eq!((dfs.roots_completed, dfs.roots_total), (Some(1), Some(1)));
    assert_eq!(dfs.estimates.state, EstimateState::Inactive);
    assert!(
        !state.presentation(Duration::from_secs(137))["native_work_active"]
            .as_bool()
            .unwrap()
    );
}

#[test]
fn frozen_job_keeps_api_and_whole_job_completion_distinct() {
    let mut state = MonitorState::default();
    let mut event = snapshot(10, 7, 7, 1);
    event.scope_kind = ScopeKind::Job;
    event.api_status = ApiStatus::Running;
    event.native_work_active = false;
    state.observe(event, Duration::from_secs(10));
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(10));
    state.latest.api_status = ApiStatus::Succeeded;
    state.latest.phase = "result_expression_serialization".into();
    let serializing = state.model(Duration::from_secs(46));
    assert_eq!(serializing.scope, ScopeKind::Job);
    assert_eq!(serializing.status, MonitorStatus::Running);
    assert_eq!(serializing.api_status, ApiStatus::Succeeded);
    assert_eq!(
        state.presentation(Duration::from_secs(46))["job_finished"],
        false
    );
    state.latest.request_finished = true;
    state.latest.status = MonitorStatus::Succeeded;
    let finished = state.model(Duration::from_secs(47));
    assert_eq!(finished.status, MonitorStatus::Succeeded);
    assert_eq!(
        state.presentation(Duration::from_secs(47))["job_finished"],
        true
    );
}

#[test]
fn frozen_frontier_keeps_collection_work_and_its_independent_age_live() {
    let mut state = MonitorState::default();
    state.observe(snapshot(10, 20, 5, 1), Duration::from_secs(100));
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(100));
    let frozen = state.model(Duration::from_secs(100)).dfs.unwrap();
    let mut notice = state.latest.clone();
    notice.elapsed = Duration::from_secs(30);
    notice.sequence += 1;
    notice.completed_unique = Some(9);
    notice.native_collection = Some(NativeCollectionSample {
        completed_roots: Some(177),
        total_roots: Some(1750),
        contribution_index: Some(3),
        contribution_total: Some(153),
        source: "native root collection".into(),
        sequence: 42,
        sample_age_seconds: Some(0.0),
    });
    notice.native_collection_sample_at = Some(Duration::from_secs(30));
    state.observe(notice, Duration::from_secs(120));
    let dfs = state.model(Duration::from_secs(156)).dfs.unwrap();
    // The frontier census is frozen; collection work and its age are live.
    assert_eq!(dfs.completed, frozen.completed);
    assert_eq!((dfs.roots_completed, dfs.roots_total), (Some(1), Some(5)));
    let collection = dfs.collection.unwrap();
    assert_eq!(
        (collection.completed_roots, collection.total_roots),
        (Some(177), Some(1750))
    );
    assert_eq!(
        (collection.contribution, collection.contributions),
        (Some(3), Some(153))
    );
    assert_eq!(collection.age_seconds, Some(36.0));
    let document = state.presentation(Duration::from_secs(156));
    assert_eq!(document["native_collection"]["sample_age_seconds"], 36.0);
    assert_eq!(
        state.latest.quantitative_sample_at,
        Some(Duration::from_secs(10))
    );
    assert!(!state.latest.request_finished);
}
