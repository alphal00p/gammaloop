use super::*;

fn snapshot(seconds: u64, entered: u64, completed: u64, roots: u64) -> DisplaySnapshot {
    DisplaySnapshot {
        sequence: seconds + 1,
        elapsed: Duration::from_secs(seconds),
        quantitative_sample_at: Some(Duration::from_secs(seconds)),
        quantitative_sample_sequence: seconds + 1,
        phase: "native coefficient combination".into(),
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
            scope: "sampled Route/Apply identities".into(),
        }),
        ..Default::default()
    }
}

fn row<'a>(rows: &'a [(String, String)], label: &str) -> &'a str {
    &rows.iter().find(|(key, _)| key == label).unwrap().1
}

#[test]
fn ordered_progress_exposes_phase_and_cancellation_scopes_without_time_percent() {
    use rustred::solver::CandidateOrderedPhaseKind;
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
        },
        native_epoch: 2,
        family_fingerprint: "ordered-family".into(),
        sampled_at: Duration::from_secs(8),
    });
    state.observe(sample, Duration::from_secs(10));
    let view = state.view(Duration::from_secs(13));
    assert_eq!(
        row(&view.overview, "Ordered phases completed / total"),
        "3 / 12 · fixed schedule; not time fraction"
    );
    assert_eq!(
        row(&view.overview, "Ordered topology stage"),
        "Reduction · support 6 · owner 2"
    );
    assert!(row(&view.overview, "Ordered phase scope").contains("request 2"));
    assert!(
        !view
            .overview
            .iter()
            .any(|(label, _)| label == "Native roots at frontier census")
    );
    assert_eq!(
        row(&view.frontier, "Input / current / peak frontier keys"),
        "7 / 13 / 29"
    );
    assert_eq!(
        row(
            &view.frontier,
            "Frontier nonzero merges / exact cancellations"
        ),
        "17 / 11"
    );
    assert_eq!(
        row(
            &view.frontier,
            "Terminal nonzero merges / exact cancellations"
        ),
        "2 / 1"
    );
    assert_eq!(row(&view.complexity, "Ordered frontier census age"), "5.0s");
    assert_eq!(
        row(
            &view.complexity,
            "Frontier coefficient monomials (current / peak)"
        ),
        "41 / 127 · num+den proxy; not bytes"
    );
    assert_eq!(
        row(
            &view.complexity,
            "Terminal coefficient monomials (current / peak)"
        ),
        "17 / 53 · num+den proxy; not bytes"
    );
    assert!(!view.estimates.iter().any(|(name, _)| name.ends_with("ETA")));
    let document = state.presentation(Duration::from_secs(13));
    assert_eq!(document["ordered_reduction"]["native_request"], 2);
    assert_eq!(document["ordered_reduction"]["sample_age_seconds"], 5.0);
    assert_eq!(document["ordered_reduction"]["frontier_keys"], 13);
    assert_eq!(
        document["ordered_reduction"]["frontier_coefficient_monomials"],
        41
    );
    assert_eq!(
        document["ordered_reduction"]["peak_frontier_coefficient_monomials"],
        127
    );
    assert_eq!(
        document["ordered_reduction"]["terminal_coefficient_monomials"],
        17
    );
    assert_eq!(
        document["ordered_reduction"]["peak_terminal_coefficient_monomials"],
        53
    );
    assert_eq!(
        document["ordered_reduction"]["phase_fraction_is_time_fraction"],
        false
    );
}

#[test]
fn ordered_batch_metadata_round_trips_and_stays_live_across_frozen_native_phases() {
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
    state.latest.reduction_cohorts.as_mut().unwrap().completed = 1;
    state.latest.reduction_cohorts.as_mut().unwrap().total = 3;
    let view = state.view(Duration::from_secs(15));
    assert!(view.frozen);
    assert_eq!(
        row(&view.overview, "Ordered batches completed / total"),
        "1 / 2 · 2L / Retain"
    );
    assert_eq!(
        row(&view.overview, "Compatible cohorts completed / total"),
        "1 / 3 · current collected batch"
    );
    assert_eq!(
        row(&view.frontier, "Completed unique identities (not closure)"),
        "5"
    );
    let encoded = serde_json::to_string(&state.presentation(Duration::from_secs(15))).unwrap();
    let decoded: serde_json::Value = serde_json::from_str(&encoded).unwrap();
    assert_eq!(decoded["ordered_batches"]["completed"], 1);
    assert_eq!(decoded["ordered_batches"]["total"], 2);
    assert_eq!(decoded["ordered_batches"]["label"], "2L / Retain");
    assert_eq!(decoded["reduction_cohorts"]["completed"], 1);
    assert_eq!(decoded["reduction_cohorts"]["total"], 3);
    assert_eq!(
        decoded["reduction_cohorts"]["fraction_is_time_fraction"],
        false
    );
    let decoded_view: MonitorView = serde_json::from_value(decoded["view"].clone()).unwrap();
    assert_eq!(
        row(
            &decoded_view.overview,
            "Compatible cohorts completed / total"
        ),
        "1 / 3 · current collected batch"
    );
    state.latest.reduction_cohorts = None;
    assert!(
        !state
            .view(Duration::from_secs(16))
            .overview
            .iter()
            .any(|(label, _)| label == "Compatible cohorts completed / total")
    );
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
        (None, "unknown"),
        (Some(0), "0 B"),
        (Some(1023), "1023 B"),
        (Some(1024), "1.0 KiB"),
        (Some(1024 * 1024), "1.0 MiB"),
        (Some(50_000_000_000), "46.6 GiB"),
        (Some(u64::MAX), "16.0 EiB"),
    ] {
        assert_eq!(ResourceUsage::display_bytes(bytes), expected);
    }
    let mut state = MonitorState::default();
    state.latest.resources = ResourceUsage {
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
    let view = state.view(Duration::ZERO);
    assert_eq!(
        row(&view.cache, "Process RSS / peak"),
        "46.6 GiB / 64.0 GiB (current / process peak; not payload)"
    );
    assert_eq!(
        row(&view.cache, "Host MemAvailable / MemFree"),
        "1.0 GiB / 0 B (usable / unused)"
    );
    assert!(row(&view.cache, "Memory units / raw values").contains("resources JSON"));
    assert!(
        row(&view.cache, "Guard MemAvailable floor / report age")
            .starts_with("46.6 GiB / unknown (explicit local policy only)")
    );
    assert_eq!(published["resources"]["rss_bytes"], 50_000_000_000_u64);
    assert_eq!(published["resources"]["mem_free_bytes"], 0);
    assert_eq!(published["resources"]["mem_total_bytes"], u64::MAX);
    assert_eq!(
        published["resources"]["guard"]["threshold_bytes"],
        50_000_000_000_u64
    );
    assert_eq!(published["resources"]["guard"]["status"], "unknown");
    assert_eq!(
        published["resources"]["guard"]["current_process_covered"],
        false
    );
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
fn display_freeze_keeps_resources_and_status_live_in_view_and_json() {
    let mut state = MonitorState::default();
    state.observe(snapshot(0, 3, 1, 0), Duration::ZERO);
    state.latest.resources = ResourceUsage::from_proc("VmRSS: 1 kB", "MemAvailable: 2 kB");
    state.latest.resource_sample_at = Some(Duration::ZERO);
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(1));
    state.latest.resources = ResourceUsage::from_proc("VmRSS: 3 kB", "MemAvailable: 4 kB");
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
        .push("CycleDetected / synthetic fixture".into());
    state.observe(failed, Duration::from_secs(6));
    let frozen = state.view(Duration::from_secs(7));
    assert!(frozen.frozen);
    assert!(row(&frozen.cache, "Process RSS / peak").starts_with("3.0 KiB /"));
    let live = serde_json::to_value(state.live_resources(Duration::from_secs(7))).unwrap();
    assert_eq!(live["rss_bytes"], 3 * 1024);
    assert_eq!(live["mem_available_bytes"], 4 * 1024);
    assert_eq!(live["sample_age_seconds"], 2.0);
    assert_eq!(live["guard"]["status"], "latched_abort");
    assert_eq!(live["guard"]["report_age_seconds"], 2.0);
    assert!(
        row(&frozen.cache, "Local guard state / current-PID coverage")
            .starts_with("latched_abort /")
    );
    let stale = state.live_resources(Duration::from_secs(9));
    assert_eq!(stale.guard.status, "unknown");
    assert_eq!(stale.guard.reported_latched_abort, Some(true));
    let published = state.presentation(Duration::from_secs(7));
    assert_eq!(published["status"], "Failed");
    assert_eq!(published["phase"], "synthetic_failed_finalization");
    assert_eq!(published["view"]["frozen"], true);
    assert_eq!(published["view"]["status"], "Failed");
    let estimates: Vec<(String, String)> =
        serde_json::from_value(published["estimates"].clone()).unwrap();
    assert_eq!(row(&estimates, "Estimate sample age"), "7.0s");
    assert!(row(&estimates, "Discovered-frontier ETA").contains("unknown"));
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
    let view = state.view(Duration::from_secs(6));
    assert_eq!(
        state
            .live_resources(Duration::from_secs(6))
            .sample_age_seconds,
        Some(1.0)
    );
    assert_eq!(row(&view.overview, "Last native observation age"), "6.0s");
    assert_eq!(row(&view.estimates, "Estimate sample age"), "6.0s");
    assert!(row(&view.estimates, "Discovered-frontier ETA").contains("unknown"));
    assert_eq!(
        row(&view.complexity, "Current canonical R / A / D"),
        "0 / 3 / 3"
    );
    assert_eq!(
        row(&view.complexity, "Caller Atom bytes / stored summands"),
        "unobserved"
    );
    assert_eq!(view.status, MonitorStatus::Running);
    assert!(view.events.is_empty());
    assert!(row(&view.cache, "Host MemAvailable / MemFree").contains("usable / unused"));
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
    let view = state.view(Duration::from_secs(103));
    assert_eq!(
        row(&view.complexity, "Native unstarted roots / output keys"),
        "2 / 3"
    );
    assert_eq!(
        row(&view.complexity, "Native deferred maps / child occurrences"),
        "4 / 11"
    );
    assert!(
        row(&view.complexity, "Caller Atom bytes / stored summands")
            .contains("937 / 7 · scalar_input · age 11.0s")
    );
    assert!(
        row(&view.complexity, "Stored-state scope / sample age").contains("native-family · 5.0s")
    );
    assert!(row(&view.complexity, "Symbolic size conventions").contains("partitions overlap"));
    assert!(
        !view
            .complexity
            .iter()
            .any(|(label, _)| label == "Full symbolic expression size")
    );
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
    let view = state.view(Duration::from_secs(13));
    assert!(row(&view.estimates, "Discovered-frontier ETA").contains("24.0s"));
    assert!(row(&view.estimates, "Discovered-frontier ETA").contains("sampled"));
    assert!(row(&view.estimates, "Discovered-frontier ETA").contains("completions/s"));
    assert!(row(&view.estimates, "Root-rate ETA").contains("40.0s"));
    assert!(row(&view.estimates, "Root-rate ETA").contains("known selected roots"));
    assert!(row(&view.estimates, "Net frontier growth").contains("+0.50"));
    assert!(row(&view.estimates, "Net frontier growth").contains("growing"));
    assert_eq!(row(&view.estimates, "Estimate sample age"), "3.0s");
    assert!(row(&view.estimates, "Estimate scope").contains("future discovery excluded"));
    assert!(row(&view.estimates, "Confidence").contains("no remaining-work certificate"));
}

#[test]
fn unknown_or_incomparable_frontier_never_becomes_a_node_eta() {
    let mut state = MonitorState::default();
    state.observe(snapshot(0, 10, 0, 0), Duration::ZERO);
    let mut next = snapshot(10, 20, 5, 0);
    next.pending_unique.as_mut().unwrap().operational_scope = false;
    state.observe(next, Duration::from_secs(10));
    let view = state.view(Duration::from_secs(10));
    assert_eq!(
        row(&view.estimates, "Discovered-frontier ETA"),
        "warmup / unknown"
    );
    assert_eq!(row(&view.estimates, "Root-rate ETA"), "warmup / unknown");
    let mut next = snapshot(20, 20, 5, 0);
    next.pending_unique = None;
    state.observe(next, Duration::from_secs(20));
    assert!(
        row(
            &state.view(Duration::from_secs(20)).frontier,
            "Pending unique identities"
        )
        .contains("unknown")
    );
}

#[test]
fn pause_and_detach_only_change_the_display() {
    let mut state = MonitorState::default();
    state.observe(snapshot(1, 10, 1, 0), Duration::from_secs(1));
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(1));
    let frozen = state.view(Duration::from_secs(1));
    state.observe(snapshot(2, 50, 20, 1), Duration::from_secs(2));
    assert_eq!(state.view(Duration::from_secs(2)).frontier, frozen.frontier);
    assert_eq!(state.latest.unique_entered, 50);
    let now = Duration::from_secs(2);
    state.latest.resources = ResourceUsage::from_proc(
        "VmRSS: 10 kB\nVmHWM: 20 kB",
        "MemTotal: 100 kB\nMemAvailable: 80 kB\nMemFree: 30 kB",
    );
    state.latest.resource_sample_at = Some(now);
    let live = serde_json::to_value(state.live_view(now)).unwrap();
    for (action, scroll, help) in [
        (UiAction::ScrollUp, 0, false),
        (UiAction::PageUp, 0, false),
        (UiAction::ToggleHelp, 0, true),
        (UiAction::ScrollDown, 1, true),
        (UiAction::PageDown, 11, true),
        (UiAction::PageUp, 1, true),
        (UiAction::ScrollUp, 0, true),
        (UiAction::ToggleHelp, 0, false),
    ] {
        state.apply(action, now);
        let view = state.view(now);
        assert_eq!(view.scroll, scroll);
        assert_eq!(view.show_help, help);
        assert_eq!(view.frontier, frozen.frontier);
        assert_eq!(serde_json::to_value(state.live_view(now)).unwrap(), live);
    }
    state.scroll = usize::MAX - 1;
    for action in [UiAction::PageDown, UiAction::ScrollDown] {
        state.apply(action, now);
        assert_eq!(state.view(now).scroll, usize::MAX);
        assert_eq!(serde_json::to_value(state.live_view(now)).unwrap(), live);
    }
    state.apply(UiAction::Detach, Duration::from_secs(2));
    assert!(state.detached());
    assert_eq!(state.latest.status, MonitorStatus::Running);
    assert_eq!(state.latest.completed_unique, Some(20));
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(2));
    assert!(!state.view(Duration::from_secs(2)).frozen);
}

#[test]
fn core_completion_is_not_final_request_success_or_mathematical_closure() {
    let mut state = MonitorState::default();
    let mut finished_core = snapshot(1, 10, 10, 1);
    finished_core.status = MonitorStatus::Succeeded;
    finished_core.request_finished = false;
    state.observe(finished_core, Duration::from_secs(1));
    let view = state.view(Duration::from_secs(1));
    assert_eq!(row(&view.overview, "Request finished"), "false");
    assert_eq!(view.status, MonitorStatus::Running);
    assert!(
        view.frontier
            .iter()
            .any(|(label, _)| label.contains("not closure"))
    );
}

#[test]
fn silent_work_advances_observation_age_and_completion_age_without_inventing_work() {
    let mut state = MonitorState::default();
    state.observe(snapshot(0, 10, 0, 0), Duration::ZERO);
    state.observe(snapshot(10, 20, 5, 1), Duration::from_secs(10));
    let view = state.view(Duration::from_secs(30));
    assert_eq!(row(&view.overview, "Snapshot age"), "20.0s");
    assert_eq!(
        row(&view.frontier, "Since last observed new completion"),
        "20.0s"
    );
    assert_eq!(state.latest.completed_unique, Some(5));
    assert_eq!(state.latest.unique_entered, 20);
}

#[test]
fn request_reset_discards_previous_rate_history() {
    let mut state = MonitorState::default();
    state.observe(snapshot(10, 10, 5, 1), Duration::from_secs(10));
    state.observe(snapshot(20, 20, 15, 2), Duration::from_secs(20));
    state.observe(snapshot(0, 1, 0, 0), Duration::from_secs(30));
    assert_eq!(state.samples.len(), 1);
    assert!(state.last_completion.is_none());
    assert_eq!(
        row(
            &state.view(Duration::from_secs(30)).estimates,
            "Root-rate ETA"
        ),
        "warmup / unknown"
    );
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
    let view = state.view(Duration::from_secs(93));
    assert_eq!(state.samples.len(), 2);
    assert_eq!(row(&view.overview, "Snapshot age"), "3.0s");
    assert_eq!(row(&view.estimates, "Estimate sample age"), "83.0s");
    assert_eq!(
        row(&view.frontier, "Since last observed new completion"),
        "83.0s"
    );
    assert_eq!(
        row(&view.estimates, "Discovered-frontier ETA"),
        "stale / unknown"
    );
    assert_eq!(row(&view.estimates, "Root-rate ETA"), "stale / unknown");
    assert_eq!(
        row(&view.estimates, "Net frontier growth"),
        "stale / unknown"
    );
    assert_eq!(
        row(&view.estimates, "Observation window"),
        "10.0s / 2 snapshots"
    );
}

#[test]
fn genuine_unchanged_sample_extends_the_observed_rate_window() {
    let mut state = MonitorState::default();
    state.observe(snapshot(0, 10, 0, 0), Duration::ZERO);
    state.observe(snapshot(10, 20, 5, 1), Duration::from_secs(10));
    state.observe(snapshot(30, 20, 5, 1), Duration::from_secs(30));
    let view = state.view(Duration::from_secs(30));
    assert_eq!(state.samples.len(), 3);
    assert_eq!(row(&view.estimates, "Estimate sample age"), "0.0s");
    assert!(row(&view.estimates, "Discovered-frontier ETA").contains("72.0s"));
    assert!(row(&view.estimates, "Root-rate ETA").contains("120.0s"));
    assert_eq!(
        row(&view.frontier, "Since last observed new completion"),
        "20.0s"
    );
}

#[test]
fn unobserved_quantitative_data_never_creates_a_rate_sample() {
    let mut state = MonitorState::default();
    let mut notice = snapshot(20, 20, 5, 1);
    notice.quantitative_sample_at = None;
    notice.quantitative_sample_sequence = 0;
    state.observe(notice, Duration::from_secs(20));
    assert!(state.samples.is_empty());
    let view = state.view(Duration::from_secs(25));
    assert_eq!(row(&view.estimates, "Estimate sample age"), "unknown");
    assert_eq!(
        row(&view.estimates, "Discovered-frontier ETA"),
        "warmup / unknown"
    );
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
    assert_eq!(
        row(
            &state.view(Duration::from_secs(30)).estimates,
            "Discovered-frontier ETA"
        ),
        "warmup / unknown"
    );
}

#[test]
fn current_canonical_degrees_are_not_maxima_or_physical_input_ranks() {
    let mut data = snapshot(1, 10, 2, 0);
    data.current_key = Some("Route [-1,2,8]".into());
    data.current_degrees = Some((1, 10, 9));
    data.complexity = vec![
        ("Max seen canonical R".into(), "99".into()),
        (
            "Physical input rank".into(),
            "unobserved; not inferred".into(),
        ),
    ];
    let mut state = MonitorState::default();
    state.observe(data.clone(), Duration::from_secs(1));
    let view = state.view(Duration::from_secs(1));
    assert_eq!(
        row(&view.complexity, "Current canonical R / A / D"),
        "1 / 10 / 9"
    );
    assert_eq!(
        row(&view.complexity, "Current canonical native key"),
        "Route [-1,2,8]"
    );
    assert_eq!(row(&view.complexity, "Max seen canonical R"), "99");
    assert_eq!(
        row(&view.complexity, "Physical input rank"),
        "unobserved; not inferred"
    );
    data.current_key = None;
    data.current_degrees = None;
    data.elapsed = Duration::from_secs(2);
    state.observe(data, Duration::from_secs(2));
    let view = state.view(Duration::from_secs(2));
    assert_eq!(
        row(&view.complexity, "Current canonical R / A / D"),
        "unknown"
    );
    assert!(row(&view.complexity, "Current canonical native key").contains("unknown"));
    assert_eq!(row(&view.complexity, "Max seen canonical R"), "99");
}

#[test]
fn phase_clock_and_resources_remain_live_while_native_panels_are_frozen() {
    let mut state = MonitorState::default();
    let mut event = snapshot(10, 7, 3, 1);
    event.phase = "retained_catalogue_exact_combination".into();
    event.native_work_active = false;
    event.phase_started_at = Duration::from_secs(10);
    state.observe(event, Duration::from_secs(20));
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(20));
    state.latest.resources.rss_bytes = Some(12345);
    state.latest.resource_sample_at = Some(Duration::from_secs(56));
    let view = state.view(Duration::from_secs(56));
    assert!(view.frozen);
    assert!(view.subtitle.contains("active 36.0s"));
    assert!(view.subtitle.contains("event age 36.0s"));
    assert!(row(&view.cache, "Process RSS / peak").contains("12.1 KiB"));
    assert!(state.presentation(Duration::from_secs(56))["caller_work"].is_null());
    assert_eq!(
        state.presentation(Duration::from_secs(56))["phase_elapsed_seconds"],
        36.0
    );
    assert_eq!(
        state.presentation(Duration::from_secs(56))["event_age_seconds"],
        36.0
    );
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
    let just_finished = state.view(Duration::from_secs(100));
    assert_eq!(just_finished.status, MonitorStatus::Running);
    assert!(
        just_finished
            .estimates
            .iter()
            .filter(|(label, _)| label.ends_with("ETA"))
            .all(|(_, value)| value.contains("indeterminate"))
    );
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(100));
    let mut caller = state.latest.clone();
    caller.phase = "retained_coefficient_materialization".into();
    caller.phase_started_at = Duration::from_secs(11);
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
    let view = state.view(Duration::from_secs(137));
    assert!(view.frozen);
    assert_eq!(view.status, MonitorStatus::Running);
    assert!(
        view.subtitle
            .contains("retained coefficient materialization")
            && view.subtitle.contains("active 36.0s")
    );
    assert_eq!(
        row(&view.overview, "Phase"),
        "retained_coefficient_materialization"
    );
    assert!(
        row(&view.overview, "Caller phase progress")
            .contains("1 / 58 master products materialized")
    );
    assert_eq!(
        row(&view.overview, "Native roots at frontier census"),
        "1 / 1"
    );
    assert!(
        view.estimates
            .iter()
            .filter(|(label, _)| label.ends_with("ETA"))
            .all(|(_, value)| value.contains("indeterminate"))
    );
    assert!(row(&view.estimates, "Estimate scope").contains("excludes caller postprocessing"));
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
    let serializing = state.view(Duration::from_secs(46));
    assert!(serializing.frozen);
    assert_eq!(serializing.status, MonitorStatus::Running);
    assert_eq!(
        row(&serializing.overview, "Reduction API status"),
        "Succeeded"
    );
    assert_eq!(row(&serializing.overview, "Whole job finished"), "false");
    state.latest.request_finished = true;
    state.latest.status = MonitorStatus::Succeeded;
    let finished = state.view(Duration::from_secs(47));
    assert_eq!(finished.status, MonitorStatus::Succeeded);
    assert_eq!(row(&finished.overview, "Whole job finished"), "true");
}

#[test]
fn frozen_frontier_keeps_collection_work_and_its_independent_age_live() {
    let mut state = MonitorState::default();
    state.observe(snapshot(10, 20, 5, 1), Duration::from_secs(100));
    state.apply(UiAction::ToggleFreeze, Duration::from_secs(100));
    let frozen_frontier = state.view(Duration::from_secs(100)).frontier;
    let mut notice = state.latest.clone();
    notice.elapsed = Duration::from_secs(30);
    notice.sequence += 1;
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
    let view = state.view(Duration::from_secs(156));
    assert!(view.frozen);
    assert_eq!(view.frontier, frozen_frontier);
    assert_eq!(
        row(&view.overview, "Native roots at frontier census"),
        "1 / 5"
    );
    assert_eq!(
        row(&view.overview, "Current collection completed roots / total"),
        "177 / 1750"
    );
    assert!(
        row(
            &view.overview,
            "Contribution currently assembled (index / total)"
        )
        .contains("3 / 153 · current item, not completed work")
    );
    assert_eq!(
        row(&view.overview, "Collection source / sample age"),
        "native root collection / 36.0s"
    );
    let document = state.presentation(Duration::from_secs(156));
    assert_eq!(document["native_collection"]["sample_age_seconds"], 36.0);
    assert_eq!(
        state.latest.quantitative_sample_at,
        Some(Duration::from_secs(10))
    );
    assert!(!state.latest.request_finished);
    assert_eq!(view.status, MonitorStatus::Running);
}
