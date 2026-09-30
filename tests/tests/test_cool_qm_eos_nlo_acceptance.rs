use std::{collections::BTreeSet, time::Instant};

use color_eyre::{Result, eyre::eyre};
use gammaloop_api::{
    commands::{Integrate, integrate::RendererOption},
    state::ProcessRef,
};
use gammaloop_integration_tests::{clean_test, get_example_cli, get_tests_workspace_path};
use gammalooprs::settings::global::Parallelisation;

#[test]
fn cool_qm_eos_nlo_matches_example_target() -> Result<()> {
    const PROCESS: &str = "cool_qm_eos";
    const INTEGRAND: &str = "NLO";
    // Benchmark recorded in examples/cli/eos/cool_qm/NLO/cool_qm_eos_NLO.toml:
    // beta=10, muB=3, MB=0.2, aS=1, m_uv=mu_r=localization_scale=2.
    const TARGET: f64 = -0.06107718321357327;
    const N_SAMPLES: usize = 750_000;
    const SEED: u64 = 1337;
    // One-worker dev-optim calibration on macOS arm64 (2026-09-30), ten
    // 75,000-sample iterations: 65 s excluding generation/compilation,
    // Im(I)=-0.061206761125173606 +/- 0.000610891261316893 (0.21 sigma).
    const BASELINE_RELATIVE_ERROR: f64 = 0.009980780719103282;
    const RELATIVE_ERROR_LIMIT: f64 = 1.5 * BASELINE_RELATIVE_ERROR;

    let test_root = get_tests_workspace_path().join("cool_qm_eos_nlo_matches_example_target");
    clean_test(&test_root);
    let mut cli = get_example_cli(
        "eos/cool_qm/NLO/cool_qm_eos_NLO.toml",
        &[],
        Some(test_root.join("state")),
        None,
        true,
    )?;
    // Override the example's ten-worker settings before executing any stage.
    cli.cli_settings.global.n_cores = Parallelisation {
        feyngen: 1,
        generate: 1,
        compile: 1,
        integrate: 1,
    };
    let integrator = &mut cli.default_runtime_settings.integrator;
    integrator.n_start = N_SAMPLES / 10;
    integrator.n_increase = 0;
    integrator.n_max = N_SAMPLES;
    integrator.seed = SEED;
    integrator.target_relative_accuracy = None;
    integrator.target_absolute_accuracy = None;

    cli.run_command("run generate")?;
    cli.run_command("set process -p cool_qm_eos -i NLO defaults")?;

    let process = ProcessRef::Unqualified(PROCESS.to_string());
    let integrand_name = INTEGRAND.to_string();
    let info = cli
        .state
        .get_integrand_info(Some(&process), Some(&integrand_name))?;
    let masters = info
        .graph_groups
        .iter()
        .flat_map(|group| &group.graphs)
        .filter(|graph| graph.is_master)
        .map(|graph| graph.name.as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        masters,
        BTreeSet::from(["GL0", "GL1", "GL2", "GL3", "GL5"]),
        "the acceptance must generate the complete cool-QM NLO example",
    );

    let started = Instant::now();
    let output = Integrate {
        process: vec![process],
        integrand_name: vec![integrand_name],
        n_cores: Some(1),
        workspace_path: Some(test_root.join("integration_workspace")),
        target: vec![format!("{PROCESS}@{INTEGRAND}=0.0,{TARGET}")],
        restart: true,
        renderer: RendererOption::Tabled,
        show_max_weight_info: false,
        no_stream_iterations: true,
        no_stream_updates: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;
    let elapsed = started.elapsed();
    let slot = output
        .single_slot()
        .ok_or_else(|| eyre!("expected one cool-QM NLO integration slot"))?;
    let estimate = &slot.integral;
    let value = estimate.result.im.0;
    let error = estimate.error.im.0;
    let relative_error = error / value.abs();
    let delta = (value - TARGET).abs();
    tracing::info!(
        ?elapsed,
        samples = estimate.neval,
        value,
        error,
        relative_error,
        target = TARGET,
        delta_sigma = delta / error,
        "Cool-QM NLO acceptance"
    );
    assert_eq!(estimate.neval, N_SAMPLES, "the sample budget must be fixed");
    assert!(
        [value, error, estimate.result.re.0, estimate.error.re.0]
            .into_iter()
            .all(f64::is_finite)
            && value < 0.0
            && error > 0.0
            && estimate.error.re.0 >= 0.0,
        "expected a finite negative imaginary integral and finite uncertainties: {estimate}",
    );
    assert!(
        estimate.result.re.0.abs() <= 2.0 * estimate.error.re.0 + 1.0e-12
            && estimate.error.re.0 <= 1.0e-12,
        "the real component must vanish: {estimate}",
    );
    assert!(
        relative_error <= RELATIVE_ERROR_LIMIT,
        "cool-QM NLO uncertainty exceeds 1.5 times the calibrated baseline: {value:e} +/- {error:e}, target={TARGET:e}, samples={}, relative error={relative_error:e}, limit={RELATIVE_ERROR_LIMIT:e}, delta/sigma={:e}",
        estimate.neval,
        delta / error,
    );
    assert!(
        delta <= 2.0 * error,
        "cool-QM NLO target mismatch: {value:e} +/- {error:e}, target={TARGET:e}, samples={}, delta/sigma={:e}",
        estimate.neval,
        delta / error,
    );

    let breakdown = slot
        .grid_breakdown
        .im
        .as_ref()
        .filter(|breakdown| breakdown.axis_label == "graph")
        .ok_or_else(|| eyre!("expected an imaginary graph-MC breakdown"))?;
    let sampled_graphs = breakdown
        .entries
        .iter()
        .filter(|entry| entry.processed_samples > 0)
        .filter_map(|entry| entry.bin_label.as_deref())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        sampled_graphs, masters,
        "every master graph must be sampled"
    );

    clean_test(&test_root);
    Ok(())
}
