use gammaloop_api::{
    state::State,
    tracing::{set_file_log_filter_override, set_stderr_log_filter_override},
};

#[test]
fn progress_layer_does_not_construct_filtered_event_payloads() {
    let log_dir = tempfile::tempdir().unwrap();
    let _state = State::new(log_dir.path(), None);

    let mut observed = Vec::new();
    for spec in [
        "off",
        "[{#generation,#profile,#compile,#summary}]=info",
        "info",
    ] {
        set_stderr_log_filter_override(Some(spec.to_string())).unwrap();
        set_file_log_filter_override(Some(spec.to_string())).unwrap();
        let evaluations = std::cell::Cell::new(0);
        // INFO keeps the regression active in builds that compile out DEBUG.
        tracing::info!(
            generation = true,
            compile = true,
            term = true,
            dump = true,
            file.network = %{
                evaluations.set(evaluations.get() + 1);
                "network payload"
            },
            "Event payload evaluation regression"
        );
        observed.push(evaluations.get());
    }

    assert_eq!(observed, [0, 0, 1]);
}
