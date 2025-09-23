#![allow(dead_code)]

use color_eyre::Result;
use criterion::{criterion_group, criterion_main, Criterion};
use gammaloop_api::{
    commands::inspect::Inspect,
    state::{RunHistory, State},
    GlobalCliSettings, OneShot, Repl,
};
use gammalooprs::settings::{GlobalSettings, RuntimeSettings};
use gammalooprs::utils::serde_utils::SmartSerde;
use gammalooprs::utils::test_utils::load_generic_model;
use pprof::criterion::{Output, PProfProfiler};
use std::{env, path::Path, time::Duration};
use tracing::{debug, warn};

fn new_cli_for_bench(test_path: impl AsRef<Path>) -> (OneShot, State) {
    let state_path = if let Ok(user_specified_state_path) = env::var("TESTS_GAMMALOOP_STATE_PATH") {
        if user_specified_state_path.to_ascii_uppercase() == "AUTO" {
            env::temp_dir().join("tests_gammaloop_state")
        } else {
            std::path::PathBuf::from(user_specified_state_path)
        }
    } else {
        test_path.as_ref().join("tests_gammaloop_state")
    };
    debug!("Using gammaloop state path: {}", state_path.display());

    let mut state = State::new_bench(state_path.clone());
    state.model = load_generic_model("sm");

    (state.new_test_cli(), state)
}

fn clean_test(state: &State) {
    if let Err(_) = env::var("GAMMALOOP_TESTS_NO_CLEAN_STATE") {
        match std::fs::remove_dir_all(state.save_path.clone()) {
            Ok(()) => debug!(
                "Gammaloop state folder '{}' deleted successfully",
                state.save_path.display()
            ),
            Err(e) => debug!(
                "Error deleting Gammaloop state folder '{}'. Error: {}",
                state.save_path.display(),
                e
            ),
        }
    } else {
        warn!("Environment variable 'GAMMALOOP_TESTS_NO_CLEAN_STATE' is set so that the gammaloop state test folder '{}' is not cleaned.",state.save_path.display());
    }
}

fn criterion_benchmark(c: &mut Criterion) -> Result<()> {
    let mut inspect_group = c.benchmark_group("inspect");
    inspect_group.measurement_time(Duration::from_secs(10));

    let inspect = Inspect {
        process_id: Some(0),
        integrand_name: Some("qqx_aaa_subtracted".into()),
        point: vec![0.1, 0.2, 0.3],
        use_f128: false,
        force_radius: false,
        momentum_space: true,
        discrete_dim: vec![0],
    };
    let mut qqx_aaa_run: RunHistory = RunHistory::from_file(
        "./benches/qqx_aaa_amplitude/qqx_aaa_subtracted_nlo_amplitude_iterative_compiled.yaml",
        "run history",
    )?;
    let (mut cli, mut state) = new_cli_for_bench("./benches/qqx_aaa_amplitude");
    let mut global_settings = GlobalCliSettings::default();
    let mut default_runtime_settings = RuntimeSettings::default();
    let _ = qqx_aaa_run.run(
        &mut state,
        &mut global_settings,
        &mut default_runtime_settings,
    )?;

    // state
    //

    inspect_group.bench_function("qqx_aaa_amp", |b| b.iter(|| inspect.run(&mut state)));
    clean_test(&state);
    Ok(())
}

criterion_group! {
    name = benches;
    config = Criterion::default().with_profiler(PProfProfiler::new(100, Output::Flamegraph(None)));
    targets = criterion_benchmark
}
criterion_main!(benches);
