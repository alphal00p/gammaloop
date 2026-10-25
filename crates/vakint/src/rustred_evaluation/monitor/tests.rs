use super::*;

fn snapshot(seconds: u64, entered: u64, completed: u64, roots: u64) -> DisplaySnapshot {
    DisplaySnapshot {
        sequence: seconds + 1,
        elapsed: Duration::from_secs(seconds),
        quantitative_sample_at: Some(Duration::from_secs(seconds)),
        quantitative_sample_sequence: seconds + 1,
        phase: "coefficient_combination".into(),
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
