use std::{
    collections::{BTreeMap, BTreeSet},
    f64::consts::PI,
};

use color_eyre::{Result, eyre::eyre};
use gammaloop_api::{
    commands::{
        Inspect, Integrate,
        evaluate_samples::{EvaluateSamplesPrecise, evaluate_samples_precise},
        integrate::RendererOption,
    },
    state::ProcessRef,
};
use gammaloop_integration_tests::{clean_test, get_example_cli, get_tests_workspace_path};
use gammalooprs::{
    integrands::evaluation::PreciseEvaluationResultOutput,
    settings::global::VectorPolarizationSumGauge,
    utils::{ArbPrec, F, FloatLike},
    uv::settings::{ApproximationType, FinalIntegrandDimension},
};
use ndarray::Array2;
use serial_test::serial;
use spenso::algebra::algebraic_traits::IsZero;
use symbolica::domains::float::Real;
use which::which;

#[test]
#[serial]
fn gamma_star_ddx_msbar_nlo_matches_the_published_graph_cross_sections() -> Result<()> {
    const LO_PROCESS: &str = "gamma_star_ddx_lo";
    const NLO_PROCESS: &str = "gamma_star_ddx_explicit_local_3d";
    const ROUTES: [(&str, bool, bool); 3] = [
        ("gamma_star_ddx_orientation_local_3d", false, false),
        (NLO_PROCESS, true, false),
        ("gamma_star_ddx_local_4d_then_cff", true, true),
    ];
    const LO: &str = "LO";
    const NLO: &str = "NLO";
    const E_CM: f64 = 400.0;
    const GEV_SQUARED_TO_PICOBARN: f64 = 0.389379304e9;
    const PUBLISHED_LO: f64 = 5.031049e-1;
    const PUBLISHED_B1_VERTEX: f64 = 5.03926e-2;
    const PUBLISHED_B2_SELF_ENERGY: f64 = -3.14956e-2;
    const PUBLISHED_NLO: f64 = PUBLISHED_B1_VERTEX + PUBLISHED_B2_SELF_ENERGY;
    let test_root = get_tests_workspace_path()
        .join("gamma_star_ddx_msbar_nlo_matches_the_published_graph_cross_sections");
    clean_test(&test_root);

    let mut cli = get_example_cli(
        "epem_a_ddx/NLO/epem_a_ddx_NLO.toml",
        &[],
        Some(test_root.join("state")),
        None,
        true,
    )?;
    {
        let generation = &mut cli.cli_settings.global.generation;
        generation.orientation_pattern = Default::default();
        generation.explicit_orientation_sum_only = true;
        generation
            .tropical_subgraph_table
            .disable_tropical_generation = true;
        generation.threshold_subtraction.enable_thresholds = true;
        generation.threshold_subtraction.disable_integrated_ct = false;
        // A virtual incoming photon is not an on-shell external state. In
        // Feynman gauge its `summed` projector is exactly -g^{mu nu}.
        generation.vector_polarization_sum_gauge = VectorPolarizationSumGauge::Feynman;
        let uv = &mut generation.uv;
        uv.softct = false;
        uv.generate_integrated = true;
        uv.subtract_uv = true;
        uv.final_integrand = FinalIntegrandDimension::ThreeD;
        uv.local_uv_cts_from_expanded_4d_integrands = false;
        uv.renormalization_prescription.log_divergent = ApproximationType::MUV;
        uv.renormalization_prescription.massive_power_divergent = ApproximationType::MUV;
        uv.renormalization_prescription.massless_power_divergent = ApproximationType::MUV;
        uv.renormalization_prescription.overrides.clear();
        uv.vakint.normalization = "MSbar".to_string();
        uv.vakint.additional_normalization = "-1".to_string();
        uv.vakint.form_exe_path = which("form")?.display().to_string();
    }
    cli.default_runtime_settings
        .subtraction
        .disable_threshold_subtraction = false;

    // Spin sums enter during integrand preprocessing, so install the direct-
    // current external state before calling `generate existing`.
    cli.run_command(&format!(
        r#"set default-runtime string '
            [kinematics]
            e_cm = {E_CM}

            [kinematics.externals]
            type = "constant"

            [kinematics.externals.data]
            momenta = [[{E_CM}, 0.0, 0.0, 0.0]]
            helicities = ["summed"]

            [general]
            mu_r = {E_CM}
            m_uv = {E_CM}
            renormalization_localization_scale = {E_CM}
            integral_unit = "none"
            disable_flux_factor = false

            [sampling]
            graphs = "monte_carlo"
            graph_names = []
            orientations = "summed"
            lmb_multichanneling = true
            lmb_channels = "summed"
            lmb_channel_weight = "ose"
        '"#,
    ))?;

    cli.run_command("import model sm-default.json")?;
    cli.run_command("run set_model_parameters")?;
    cli.run_command(&format!(
        r#"generate xs a > d d~ | a d d~ g ghG ghG~ QCD^2==0 QED^2==2 [{{{{1}}}} QCD=0]
            --numerator-grouping group_identical_graphs_up_to_scalar_rescaling
            --symmetrize-left-right-states true
            -p {LO_PROCESS} -i {LO} --global-prefactor-num "1/2" --only-diagrams"#,
    ))?;
    cli.run_command(&format!("generate existing -p {LO_PROCESS} -i {LO}"))?;

    // The paper averages over the two incoming-current polarizations. The
    // factor 1/2 is kept outside `summed`, which supplies -g^{mu nu}.
    for (process, explicit_orientation_sum_only, project_local_4d_to_cff) in ROUTES {
        {
            let generation = &mut cli.cli_settings.global.generation;
            generation.explicit_orientation_sum_only = explicit_orientation_sum_only;
            generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d_to_cff;
        }
        cli.run_command(&format!(
            r#"generate xs a > d d~ | a d d~ g ghG ghG~ QCD^2==2 QED^2==2 [{{{{2}}}} QCD=1]
                --numerator-grouping group_identical_graphs_up_to_scalar_rescaling
                --symmetrize-left-right-states true
                -p {process} -i {NLO} --global-prefactor-num "-1𝑖/2" --only-diagrams"#,
        ))?;
        cli.run_command(&format!("generate existing -p {process} -i {NLO}"))?;

        let info = cli.state.get_integrand_info(
            Some(&ProcessRef::Unqualified(process.to_string())),
            Some(&NLO.to_string()),
        )?;
        let masters = info
            .graph_groups
            .iter()
            .flat_map(|group| group.graphs.iter())
            .filter(|graph| graph.is_master)
            .map(|graph| graph.name.clone())
            .collect::<BTreeSet<_>>();
        assert_eq!(
            masters,
            BTreeSet::from(["GL0".to_string(), "GL2".to_string()]),
            "the {process} direct-current acceptance must contain exactly the self-energy and vertex supergraphs",
        );
    }
    {
        let generation = &mut cli.cli_settings.global.generation;
        generation.explicit_orientation_sum_only = true;
        generation.uv.local_uv_cts_from_expanded_4d_integrands = false;
    }

    let first_process_ref = ProcessRef::Unqualified(NLO_PROCESS.to_string());
    let (process_id, resolved_nlo_name) = cli
        .state
        .find_integrand_ref(Some(&first_process_ref), Some(&NLO.to_string()))?;
    let parameters = cli
        .state
        .resolve_effective_model_parameter_card_for_integrand(process_id, &resolved_nlo_name)?;
    let (alpha_s_re, alpha_s_im) = parameters
        .data
        .get("aS")
        .ok_or_else(|| eyre!("the effective NLO model parameters do not contain aS"))?;
    let (alpha_qed_inverse_re, alpha_qed_inverse_im) = parameters
        .data
        .get("aEWM1")
        .ok_or_else(|| eyre!("the effective NLO model parameters do not contain aEWM1"))?;
    assert_eq!(alpha_s_im.0, 0.0, "the acceptance expects real aS");
    assert_eq!(
        alpha_qed_inverse_im.0, 0.0,
        "the acceptance expects real aEWM1"
    );
    let alpha_s = alpha_s_re.0;
    let alpha_qed = alpha_qed_inverse_re.0.recip();

    // Contracting the direct current with the massless e+e- tensor and adding
    // the photon propagator gives the normalization used by the e+e- setup:
    //
    // sigma(e+e- -> gamma* -> X)[pb]
    //   = sigma(gamma* -> X) * 2 e^2/(3 E_cm^3) * (GeV^-2 -> pb).
    let direct_current_to_epem_pb =
        2.0 * (4.0 * PI * alpha_qed) / (3.0 * E_CM.powi(3)) * GEV_SQUARED_TO_PICOBARN;
    let published_epem_lo_pb = PUBLISHED_LO * direct_current_to_epem_pb;
    let published_epem_nlo_pb = PUBLISHED_NLO * direct_current_to_epem_pb;
    let analytic_epem_lo_pb = 4.0 * PI * alpha_qed.powi(2) * 3.0 * (1.0 / 3.0_f64).powi(2)
        / (3.0 * E_CM.powi(2))
        * GEV_SQUARED_TO_PICOBARN;
    assert!(
        (published_epem_lo_pb - analytic_epem_lo_pb).abs() <= 5.0e-5 * analytic_epem_lo_pb,
        "the published direct-current LO and analytic e+e- normalization disagree: converted={published_epem_lo_pb:e} pb, analytic={analytic_epem_lo_pb:e} pb",
    );
    assert!(
        (PUBLISHED_NLO / PUBLISHED_LO - alpha_s / PI).abs() <= 2.0e-5 * alpha_s / PI,
        "the published direct-current inclusive ratio is not alpha_s/pi: published={:e}, alpha_s/pi={:e}",
        PUBLISHED_NLO / PUBLISHED_LO,
        alpha_s / PI,
    );

    cli.run_command(&format!(
        r#"set process -p {LO_PROCESS} -i {LO} kv
            sampling.graphs="summed"
            sampling.orientations="summed"
            sampling.lmb_multichanneling=true
            sampling.lmb_channels="summed"
            sampling.lmb_channel_weight="ose""#,
    ))?;
    let (_, lo_probe) = Inspect {
        process: Some(ProcessRef::Unqualified(LO_PROCESS.to_string())),
        integrand_name: Some(LO.to_string()),
        point: vec![0.19, 0.43, 0.71],
        ..Default::default()
    }
    .run(&mut cli.state)?;
    assert!(
        lo_probe.im < 0.0 && lo_probe.re.abs() <= 1.0e-12 * lo_probe.im.abs(),
        "the direct-current LO convention must be purely negative imaginary, got {lo_probe:e}",
    );
    cli.run_command(&format!(
        r#"set process -p {LO_PROCESS} -i {LO} kv
            integrator.integrated_phase="imag"
            integrator.min_samples_for_update=5000
            integrator.n_start=5000
            integrator.n_increase=5000
            integrator.n_max=10000
            integrator.target_relative_accuracy=1.0e-12
            integrator.seed=7331
            integrator.discrete_dim_learning_rate=0.0
            integrator.continuous_dim_learning_rate=0.1"#,
    ))?;
    let lo_output = Integrate {
        process: vec![ProcessRef::Unqualified(LO_PROCESS.to_string())],
        integrand_name: vec![LO.to_string()],
        n_cores: Some(10),
        workspace_path: Some(test_root.join("lo_integration_workspace")),
        restart: true,
        renderer: RendererOption::Tabled,
        show_max_weight_info: false,
        no_stream_iterations: true,
        no_stream_updates: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;
    let lo_estimate = lo_output
        .single_slot_integral()
        .ok_or_else(|| eyre!("expected one direct-current LO integration slot"))?;
    let lo_value = -lo_estimate.result.im.0;
    let lo_error = lo_estimate.error.im.0.abs();
    assert!(
        lo_error / lo_value <= 0.05,
        "direct-current LO uncertainty is too large: {lo_value:e} ± {lo_error:e}",
    );
    assert!(
        (lo_value - PUBLISHED_LO).abs() <= 3.0 * lo_error,
        "direct-current LO mismatch: {lo_value:e} ± {lo_error:e}, published={PUBLISHED_LO:e}",
    );
    assert!(
        (lo_value * direct_current_to_epem_pb - published_epem_lo_pb).abs()
            <= 3.0 * lo_error * direct_current_to_epem_pb,
        "converted direct-current LO does not reproduce the e+e- normalization",
    );

    for (process, _, _) in ROUTES {
        cli.run_command(&format!(
            r#"set process -p {process} -i {NLO} kv
                sampling.graphs="summed"
                sampling.orientations="summed"
                sampling.lmb_multichanneling=true
                sampling.lmb_channels="summed"
                sampling.lmb_channel_weight="ose""#,
        ))?;
    }
    let reference_process = NLO_PROCESS;
    let (_, nlo_probe) = Inspect {
        process: Some(ProcessRef::Unqualified(reference_process.to_string())),
        integrand_name: Some(NLO.to_string()),
        point: vec![17.0, -31.0, 43.0, -29.0, 53.0, 71.0],
        momentum_space: true,
        ..Default::default()
    }
    .run(&mut cli.state)?;
    assert!(
        nlo_probe.re.is_finite()
            && nlo_probe.im.is_finite()
            && nlo_probe.re.abs() > 0.0
            && nlo_probe.im.abs() <= 1.0e-8 * nlo_probe.re.abs(),
        "the {reference_process} NLO pointwise probe must be finite, nonzero and real, got {nlo_probe:e}",
    );

    let points = Array2::from_shape_vec(
        (2, 6),
        vec![
            17.0, -31.0, 43.0, -29.0, 53.0, 71.0, 17.0, -31.0, 43.0, -29.0, 53.0, 71.0,
        ],
    )?;
    let graph_names = ["GL0", "GL2"];
    // Retaining one eighth of ArbPrec's requested precision permits substantial
    // expression-dependent loss while still testing the precision-increase path.
    let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
    let mut route_values = Vec::new();
    for (process, _, _) in ROUTES {
        let process_id = cli
            .state
            .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
        let result = evaluate_samples_precise(
            &mut cli.state,
            &EvaluateSamplesPrecise {
                process_id: Some(process_id),
                integrand_name: Some(NLO.to_string()),
                use_arb_prec: true,
                minimal_output: true,
                return_generated_events: Some(false),
                momentum_space: true,
                points: points.view(),
                integrator_weights: None,
                discrete_dims: None,
                graph_names: Some(graph_names.map(|graph| Some(graph.to_string())).to_vec()),
                orientations: Some(vec![None; graph_names.len()]),
            },
        )?;
        let values = result
            .samples
            .into_iter()
            .zip(graph_names)
            .map(|(sample, graph_name)| match sample.evaluation {
                PreciseEvaluationResultOutput::Arb(result) => Ok(result.integrand_result),
                evaluation => Err(eyre!(
                    "{process} {graph_name} comparison requested Arb but got {:?}",
                    evaluation.precision()
                )),
            })
            .collect::<Result<Vec<_>>>()?;
        route_values.push((process, values));
    }
    let reference = &route_values[1];
    for (graph_index, graph_name) in graph_names.iter().copied().enumerate() {
        let reference_value = &reference.1[graph_index];
        for actual in [&route_values[0], &route_values[2]] {
            let actual_value = &actual.1[graph_index];
            let distance = (actual_value.clone() - reference_value.clone()).norm().re;
            let actual_norm = actual_value.norm().re;
            let reference_norm = reference_value.norm().re;
            let scale = if actual_norm > reference_norm {
                actual_norm
            } else {
                reference_norm
            };
            let relative_distance = if scale.is_zero() {
                distance.clone()
            } else {
                distance.clone() / scale
            };
            assert!(
                relative_distance <= tolerance.clone(),
                "{graph_name} differs between precise routes {} and {}: actual={:e}, reference={:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}",
                actual.0,
                reference.0,
                actual_value,
                reference_value,
            );
        }
    }

    // The Arb checks transfer the explicit-local-3D published oracle to both
    // alternative routes, so only the reference needs the expensive graph MC.
    cli.run_command(&format!(
        r#"set process -p {reference_process} -i {NLO} kv
            sampling.graphs="monte_carlo"
            sampling.graph_names=[]
            sampling.orientations="summed"
            sampling.lmb_multichanneling=true
            sampling.lmb_channels="summed"
            sampling.lmb_channel_weight="ose"
            integrator.integrated_phase="real"
            integrator.min_samples_for_update=20000
            integrator.n_start=20000
            integrator.n_increase=20000
            integrator.n_max=40000
            integrator.target_relative_accuracy=1.0e-12
            integrator.seed=220311038
            integrator.discrete_dim_learning_rate=0.0
            integrator.continuous_dim_learning_rate=0.1"#,
    ))?;
    let output = Integrate {
        process: vec![ProcessRef::Unqualified(reference_process.to_string())],
        integrand_name: vec![NLO.to_string()],
        n_cores: Some(10),
        workspace_path: Some(test_root.join("nlo_graph_mc_integration_workspace")),
        restart: true,
        renderer: RendererOption::Tabled,
        show_max_weight_info: false,
        no_stream_iterations: true,
        no_stream_updates: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;
    let slot = output
        .single_slot()
        .ok_or_else(|| eyre!("expected one direct-current NLO integration slot"))?;
    let nlo_value = slot.integral.result.re.0;
    let nlo_error = slot.integral.error.re.0.abs();
    assert!(
        nlo_value > 0.0 && nlo_error / nlo_value <= 0.15,
        "direct-current NLO uncertainty is too large: {nlo_value:e} ± {nlo_error:e}",
    );
    assert!(
        (nlo_value - PUBLISHED_NLO).abs() <= 3.0 * nlo_error,
        "direct-current NLO mismatch: {nlo_value:e} ± {nlo_error:e}, published={PUBLISHED_NLO:e}",
    );
    assert!(
        (nlo_value * direct_current_to_epem_pb - published_epem_nlo_pb).abs()
            <= 3.0 * nlo_error * direct_current_to_epem_pb,
        "converted direct-current NLO does not reproduce the e+e- normalization",
    );

    let graph_breakdown = slot
        .grid_breakdown
        .re
        .as_ref()
        .filter(|breakdown| breakdown.axis_label == "graph")
        .ok_or_else(|| eyre!("direct-current graph-MC result has no real graph breakdown"))?;
    let graph_estimates = graph_breakdown
        .entries
        .iter()
        .map(|entry| {
            let label = entry
                .bin_label
                .clone()
                .unwrap_or_else(|| entry.bin_index.to_string());
            (label, (entry.value.0, entry.error.0.abs()))
        })
        .collect::<BTreeMap<_, _>>();
    assert_eq!(
        graph_estimates.keys().cloned().collect::<BTreeSet<_>>(),
        BTreeSet::from(["GL0".to_string(), "GL2".to_string()]),
        "the direct-current graph-MC breakdown must contain exactly GL0 and GL2",
    );
    // In this two-master generation GL2 is the B.1 vertex graph and GL0 is
    // the B.2 self-energy graph, so the direct-current rows can be tested
    // separately without invoking a graphwise e+e- Ward decomposition.
    for (graph, target) in [
        ("GL0", PUBLISHED_B2_SELF_ENERGY),
        ("GL2", PUBLISHED_B1_VERTEX),
    ] {
        let (value, error) = graph_estimates
            .get(graph)
            .copied()
            .ok_or_else(|| eyre!("graph-MC result is missing {graph}"))?;
        assert!(
            error / target.abs() <= 0.15,
            "{graph} uncertainty is too large: {value:e} ± {error:e}, published={target:e}",
        );
        assert!(
            (value - target).abs() <= 3.0 * error,
            "{graph} mismatch: {value:e} ± {error:e}, published={target:e}",
        );
    }
    let graph_sum = graph_estimates
        .values()
        .map(|(value, _)| value)
        .sum::<f64>();
    let graph_sum_error = graph_estimates
        .values()
        .map(|(_, error)| error.powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(
        (graph_sum - nlo_value).abs() <= 4.0 * graph_sum_error.hypot(nlo_error),
        "sum of graph-MC rows does not reproduce the slot integral: graphs={graph_sum:e} ± {graph_sum_error:e}, slot={nlo_value:e} ± {nlo_error:e}",
    );
    assert!(
        (graph_sum - PUBLISHED_NLO).abs() <= 3.0 * graph_sum_error,
        "sum of graph-MC NLO rows mismatches the paper: {graph_sum:e} ± {graph_sum_error:e}, published={PUBLISHED_NLO:e}",
    );

    clean_test(&test_root);
    Ok(())
}
