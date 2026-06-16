use std::{collections::BTreeSet, f64::consts::PI};

use color_eyre::{Result, eyre::eyre};
use gammaloop_api::{
    commands::{
        Inspect, Integrate, Profile,
        evaluate_samples::{EvaluateSamplesPrecise, evaluate_samples_precise},
        integrate::RendererOption,
        profile::UltraVioletProfile,
    },
    state::ProcessRef,
};
use gammaloop_integration_tests::{clean_test, get_example_cli, get_tests_workspace_path};
use gammalooprs::{
    integrands::evaluation::PreciseEvaluationResultOutput,
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
fn epem_a_ttx_msbar_nlo_matches_the_published_inclusive_ratio_in_all_local_uv_routes() -> Result<()>
{
    const LO_PROCESS: &str = "epem_a_ttx_lo";
    const EXPLICIT_3D_PROCESS: &str = "epem_a_ttx_explicit_local_3d";
    const LOCAL_ONLY_ROUTES: [(&str, bool, bool); 3] = [
        ("epem_a_ttx_local_only_3d", false, false),
        ("epem_a_ttx_local_only_explicit_3d", true, false),
        ("epem_a_ttx_local_only_4d_then_cff", true, true),
    ];
    const LO: &str = "LO";
    const NLO: &str = "NLO";
    const E_CM: f64 = 600.0;
    const BEAM_ENERGY: f64 = E_CM / 2.0;
    const GEV_SQUARED_TO_PICOBARN: f64 = 0.389379304e9;
    const ROUTES: [(&str, bool); 2] = [
        ("epem_a_ttx_local_3d", false),
        ("epem_a_ttx_local_4d_then_cff", true),
    ];

    // The published absolute values are for gamma* -> t t~, whereas this
    // acceptance computes e+ e- -> gamma* -> t t~. A direct off-shell
    // gamma* run would have to generate the summed Feynman-gauge vector
    // projector -g^{mu nu}; here the valid external-lepton spin sums and
    // common photon factors cancel only in the inclusive NLO/LO ratio.
    const PUBLISHED_GAMMA_STAR_LO: f64 = 2.876302;
    const PUBLISHED_GAMMA_STAR_NLO: f64 = 0.201520;
    const PUBLISHED_NLO_OVER_LO: f64 = PUBLISHED_GAMMA_STAR_NLO / PUBLISHED_GAMMA_STAR_LO;

    let test_root = get_tests_workspace_path()
        .join("epem_a_ttx_msbar_nlo_matches_the_published_inclusive_ratio_in_all_local_uv_routes");
    clean_test(&test_root);

    // Reuse the production LU settings, while generating the two-graph top
    // process explicitly because the shipped example is the massless case.
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
        generation.explicit_orientation_sum_only = false;
        generation
            .tropical_subgraph_table
            .disable_tropical_generation = true;
        generation.threshold_subtraction.enable_thresholds = true;
        generation.threshold_subtraction.disable_integrated_ct = false;
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

    cli.run_command("import model sm-default.json")?;
    cli.run_command("run set_model_parameters")?;
    cli.run_command(
        r#"generate xs e+ e- > t t~ | e+ e- t t~ g ghG ghG~ a QCD^2==0 QED^2==4 [{{1}} QCD=0]
                --numerator-grouping group_identical_graphs_up_to_scalar_rescaling
                --symmetrize-left-right-states true
                -p epem_a_ttx_lo -i LO --only-diagrams"#,
    )?;
    cli.run_command("generate existing -p epem_a_ttx_lo -i LO")?;

    for (process, project_local_4d_to_cff) in ROUTES {
        {
            let generation = &mut cli.cli_settings.global.generation;
            generation.explicit_orientation_sum_only = project_local_4d_to_cff;
            generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d_to_cff;
        }
        cli.run_command(&format!(
                r#"generate xs e+ e- > t t~ | e+ e- t t~ g ghG ghG~ a QCD^2==2 QED^2==4 [{{{{2}}}} QCD=1]
                    --numerator-grouping group_identical_graphs_up_to_scalar_rescaling
                    --symmetrize-left-right-states true
                    -p {process} -i NLO --global-prefactor-num "1𝑖" --only-diagrams"#,
            ))?;
        cli.run_command(&format!("generate existing -p {process} -i NLO"))?;

        let process_ref = ProcessRef::Unqualified(process.to_string());
        let info = cli
            .state
            .get_integrand_info(Some(&process_ref), Some(&NLO.to_string()))?;
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
            "the {process} acceptance must contain exactly the self-energy and vertex supergraphs",
        );
    }

    let lo_process_ref = ProcessRef::Unqualified(LO_PROCESS.to_string());
    let (lo_process_id, resolved_lo_name) = cli
        .state
        .find_integrand_ref(Some(&lo_process_ref), Some(&LO.to_string()))?;
    let model_parameters = cli
        .state
        .resolve_effective_model_parameter_card_for_integrand(lo_process_id, &resolved_lo_name)?;
    for (parameter, expected) in [
        ("aS", 0.118),
        ("aEWM1", 132.507),
        ("MT", 173.0),
        ("ymt", 173.0),
        ("WT", 0.0),
    ] {
        let (re, im) = model_parameters
            .data
            .get(parameter)
            .ok_or_else(|| eyre!("the effective LO model parameters do not contain {parameter}"))?;
        assert_eq!(
            im.0, 0.0,
            "the published benchmark expects a real {parameter}",
        );
        assert_eq!(
            re.0, expected,
            "the effective {parameter} does not match the published benchmark",
        );
    }
    let alpha_qed = model_parameters.data["aEWM1"].0.0.recip();
    // Eq. (7.1) of arXiv:2203.11038 closes the external lepton trace;
    // Q_e^2=1 for this process and the generated cross-section is in pb.
    let gamma_star_to_epem_pb =
        2.0 * (4.0 * PI * alpha_qed) / (3.0 * E_CM.powi(3)) * GEV_SQUARED_TO_PICOBARN;
    let published_lo_pb = PUBLISHED_GAMMA_STAR_LO * gamma_star_to_epem_pb;
    let published_nlo_pb = PUBLISHED_GAMMA_STAR_NLO * gamma_star_to_epem_pb;

    // Keep this representation in a separate process. Integrand generation
    // annotates the processed graph, so regenerating the orientation-local
    // process in place would not be a representation-independent comparison.
    {
        let generation = &mut cli.cli_settings.global.generation;
        generation.explicit_orientation_sum_only = true;
        generation.uv.local_uv_cts_from_expanded_4d_integrands = false;
    }
    cli.run_command(&format!(
        r#"generate xs e+ e- > t t~ | e+ e- t t~ g ghG ghG~ a QCD^2==2 QED^2==4 [{{{{2}}}} QCD=1]
            --numerator-grouping group_identical_graphs_up_to_scalar_rescaling
            --symmetrize-left-right-states true
            -p {EXPLICIT_3D_PROCESS} -i NLO --global-prefactor-num "1𝑖" --only-diagrams"#,
    ))?;
    cli.run_command(&format!(
        "generate existing -p {EXPLICIT_3D_PROCESS} -i NLO"
    ))?;

    // Isolate the original graph plus unintegrated local UV subtraction. These
    // processes deliberately omit integrated UV counterterms, so equality is
    // a local test of the two 3D routes and the local-4D-to-CFF projection.
    cli.cli_settings.global.generation.uv.generate_integrated = false;
    for (process, explicit_orientation_sum_only, project_local_4d_to_cff) in LOCAL_ONLY_ROUTES {
        {
            let generation = &mut cli.cli_settings.global.generation;
            generation.explicit_orientation_sum_only = explicit_orientation_sum_only;
            generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d_to_cff;
        }
        cli.run_command(&format!(
            r#"generate xs e+ e- > t t~ | e+ e- t t~ g ghG ghG~ a QCD^2==2 QED^2==4 [{{{{2}}}} QCD=1]
                --numerator-grouping group_identical_graphs_up_to_scalar_rescaling
                --symmetrize-left-right-states true
                -p {process} -i NLO --global-prefactor-num "1𝑖" --only-diagrams"#,
        ))?;
        cli.run_command(&format!("generate existing -p {process} -i NLO"))?;

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
            "the local-only {process} acceptance must contain exactly the self-energy and vertex supergraphs",
        );
    }
    {
        let generation = &mut cli.cli_settings.global.generation;
        generation.explicit_orientation_sum_only = true;
        generation.uv.local_uv_cts_from_expanded_4d_integrands = false;
        generation.uv.generate_integrated = true;
    }

    for process in std::iter::once(LO_PROCESS)
        .chain(ROUTES.map(|(process, _)| process))
        .chain(std::iter::once(EXPLICIT_3D_PROCESS))
        .chain(LOCAL_ONLY_ROUTES.map(|(process, _, _)| process))
    {
        cli.run_command(&format!(
            r#"set process -p {process} string '
                    [kinematics]
                    e_cm = {E_CM}

                    [kinematics.externals]
                    type = "constant"

                    [kinematics.externals.data]
                    momenta = [
                        [{BEAM_ENERGY}, 0.0, 0.0, {BEAM_ENERGY}],
                        [{BEAM_ENERGY}, 0.0, 0.0, -{BEAM_ENERGY}],
                    ]
                    helicities = ["summed_averaged", "summed_averaged"]

                    [general]
                    mu_r = {E_CM}
                    m_uv = {E_CM}
                    renormalization_localization_scale = {E_CM}
                    integral_unit = "picobarn"

                    [sampling]
                    graphs = "monte_carlo"
                    graph_names = []
                    orientations = "summed"
                    lmb_multichanneling = true
                    lmb_channels = "summed"
                    lmb_channel_weight = "ose"
                '"#,
        ))?;
    }

    let mut reference_profile_signatures = None;
    let profile_routes = [ROUTES[0].0, EXPLICIT_3D_PROCESS, ROUTES[1].0];
    for (route_index, process) in profile_routes.into_iter().enumerate() {
        let analysis = Profile::UltraViolet(UltraVioletProfile {
            process: Some(ProcessRef::Unqualified(process.to_string())),
            integrand_name: Some(NLO.to_string()),
            n_points: 6,
            seed: Some(9200 + route_index as u64),
            ..Default::default()
        })
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
        assert_eq!(
            analysis.graphs.len(),
            2,
            "the {process} default divergent-only profile must cover exactly two NLO supergraphs",
        );
        assert_eq!(
            analysis
                .graphs
                .iter()
                .map(|graph| graph.graph_name.clone())
                .collect::<BTreeSet<_>>(),
            BTreeSet::from(["GL0".to_string(), "GL2".to_string()]),
            "the {process} default divergent-only profile must cover both NLO supergraphs",
        );
        let profile_signatures = analysis
            .graphs
            .iter()
            .flat_map(|graph| {
                graph.lmbs.iter().flat_map(|lmb| {
                    lmb.subsets.iter().map(|subset| {
                        (
                            graph.graph_name.clone(),
                            subset
                                .free
                                .iter()
                                .copied()
                                .map(usize::from)
                                .collect::<Vec<_>>(),
                            subset.initial_dod,
                        )
                    })
                })
            })
            .collect::<BTreeSet<_>>();
        assert!(
            profile_signatures
                .iter()
                .all(|(_, _, initial_dod)| *initial_dod >= 0),
            "the {process} default divergent-only profile selected a convergent limit: {profile_signatures:?}",
        );
        if let Some(reference) = &reference_profile_signatures {
            assert_eq!(
                &profile_signatures, reference,
                "the {process} default divergent-only profile selected a different deterministic limit set",
            );
        } else {
            reference_profile_signatures = Some(profile_signatures);
        }
        for graph in &analysis.graphs {
            assert!(
                graph
                    .lmbs
                    .iter()
                    .flat_map(|lmb| &lmb.subsets)
                    .any(|subset| subset.estimated_dod().is_some()),
                "the {process} default divergent-only profile must fit a nonvanishing limit for {}",
                graph.graph_name,
            );
        }
        let profile = analysis.pass_fail(-0.9);
        assert!(
            profile.total > 0,
            "the {process} default divergent-only UV profile must test at least one limit",
        );
        assert_eq!(
            profile.failed, 0,
            "the {process} default divergent-only UV profile failed:\n{profile}",
        );
    }

    let route_processes = [ROUTES[0].0, EXPLICIT_3D_PROCESS, ROUTES[1].0];
    let local_only_processes = LOCAL_ONLY_ROUTES.map(|(process, _, _)| process);
    for process in route_processes.into_iter().chain(local_only_processes) {
        cli.run_command(&format!(
            r#"set process -p {process} -i {NLO} kv
                sampling.graphs="summed"
                sampling.orientations="summed"
                sampling.lmb_multichanneling=true
                sampling.lmb_channels="summed"
                sampling.lmb_channel_weight="ose""#,
        ))?;
    }

    // A cheap f64 probe selects the integrated phase. Route equivalence itself
    // is checked below in Arb so that precision-escalation losses remain visible.
    let (_, nlo_probe) = Inspect {
        process: Some(ProcessRef::Unqualified(EXPLICIT_3D_PROCESS.to_string())),
        integrand_name: Some(NLO.to_string()),
        point: vec![17.0, -31.0, 43.0, -29.0, 53.0, 71.0],
        momentum_space: true,
        ..Default::default()
    }
    .run(&mut cli.state)?;
    assert!(
        nlo_probe.re.is_finite()
            && nlo_probe.im.is_finite()
            && nlo_probe.re.abs().max(nlo_probe.im.abs()) > 0.0,
        "the explicit-local-3D NLO pointwise probe must be finite and nonzero, got {nlo_probe:e}",
    );
    assert!(
        nlo_probe.re.abs().min(nlo_probe.im.abs())
            <= 1.0e-8 * nlo_probe.re.abs().max(nlo_probe.im.abs()),
        "the explicit-local-3D NLO probe is not phase-pure enough for a magnitude acceptance: {nlo_probe:e}",
    );
    let nlo_phase = if nlo_probe.re.abs() >= nlo_probe.im.abs() {
        "real"
    } else {
        "imag"
    };

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
    for (comparison, processes) in [
        ("integrated", route_processes),
        ("local-only", local_only_processes),
    ] {
        let mut route_values = Vec::new();
        for process in processes {
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
                        "{comparison} {process} {graph_name} comparison requested Arb but got {:?}",
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
                    "{comparison} {graph_name} differs between precise routes {} and {}: actual={:e}, reference={:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}",
                    actual.0,
                    reference.0,
                    actual_value,
                    reference_value,
                );
            }
        }
    }

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
        lo_probe.re.is_finite()
            && lo_probe.im.is_finite()
            && lo_probe.re.abs().max(lo_probe.im.abs()) > 0.0,
        "the LO phase probe must be finite and nonzero, got {lo_probe:e}",
    );
    assert!(
        lo_probe.re.abs().min(lo_probe.im.abs())
            <= 1.0e-8 * lo_probe.re.abs().max(lo_probe.im.abs()),
        "the LO probe is not phase-pure enough for a magnitude acceptance: {lo_probe:e}",
    );
    let lo_phase = if lo_probe.re.abs() >= lo_probe.im.abs() {
        "real"
    } else {
        "imag"
    };
    cli.run_command(&format!(
        r#"set process -p {LO_PROCESS} -i {LO} kv
                integrator.integrated_phase="{lo_phase}"
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
        .ok_or_else(|| eyre!("expected one LO integration slot"))?;
    let lo_value = lo_estimate.result.re.0.hypot(lo_estimate.result.im.0);
    let lo_error = lo_estimate.error.re.0.hypot(lo_estimate.error.im.0);
    assert!(lo_value > 0.0, "expected a nonzero LO cross-section");
    assert!(
        lo_error / lo_value <= 0.05,
        "LO uncertainty is too large for the normalization acceptance: {lo_value:e} ± {lo_error:e}",
    );
    let lo_absolute_delta = (lo_value - published_lo_pb).abs();
    let lo_absolute_tolerance = (3.0 * lo_error).max(0.05 * published_lo_pb);
    assert!(
        lo_absolute_delta <= lo_absolute_tolerance,
        "LO absolute normalization mismatch: |LO|={lo_value:e} ± {lo_error:e} pb, converted published value={published_lo_pb:e} pb, |delta|={lo_absolute_delta:e}, tolerance={lo_absolute_tolerance:e}",
    );

    // The Arb checks establish route equivalence locally, so only the explicit
    // local-3D route needs the comparatively expensive inclusive Monte Carlo.
    cli.run_command(&format!(
        r#"set process -p {EXPLICIT_3D_PROCESS} -i {NLO} kv
            sampling.graphs="summed"
            sampling.orientations="summed"
            sampling.lmb_multichanneling=true
            sampling.lmb_channels="summed"
            sampling.lmb_channel_weight="ose"
            integrator.integrated_phase="{nlo_phase}"
            integrator.min_samples_for_update=10000
            integrator.n_start=10000
            integrator.n_increase=10000
            integrator.n_max=20000
            integrator.target_relative_accuracy=1.0e-12
            integrator.seed=1337
            integrator.discrete_dim_learning_rate=0.0
            integrator.continuous_dim_learning_rate=0.1"#,
    ))?;
    let nlo_output = Integrate {
        process: vec![ProcessRef::Unqualified(EXPLICIT_3D_PROCESS.to_string())],
        integrand_name: vec![NLO.to_string()],
        n_cores: Some(10),
        workspace_path: Some(test_root.join("explicit_local_3d_nlo_integration_workspace")),
        restart: true,
        renderer: RendererOption::Tabled,
        show_max_weight_info: false,
        no_stream_iterations: true,
        no_stream_updates: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;
    let nlo_estimate = nlo_output
        .single_slot_integral()
        .ok_or_else(|| eyre!("expected one explicit-local-3D NLO integration slot"))?;
    let nlo_value = nlo_estimate.result.re.0.hypot(nlo_estimate.result.im.0);
    let nlo_error = nlo_estimate.error.re.0.hypot(nlo_estimate.error.im.0);
    let nlo_relative_error = nlo_error / nlo_value;
    assert!(
        nlo_relative_error.is_finite() && nlo_relative_error <= 0.15,
        "explicit-local-3D NLO uncertainty is too large for the normalization acceptance: {nlo_value:e} ± {nlo_error:e} ({:.2}%)",
        100.0 * nlo_relative_error,
    );

    let expected_nlo = PUBLISHED_NLO_OVER_LO * lo_value;
    let combined_error = nlo_error.hypot(PUBLISHED_NLO_OVER_LO * lo_error);
    let delta = (nlo_value - expected_nlo).abs();
    assert!(
        combined_error.is_finite() && delta <= 3.0 * combined_error,
        "explicit-local-3D inclusive NLO normalization mismatch: |NLO|={nlo_value:e} ± {nlo_error:e}, published-ratio*LO={expected_nlo:e} ± {:e}, |delta|={delta:e}, delta/sigma={:e}",
        PUBLISHED_NLO_OVER_LO * lo_error,
        delta / combined_error,
    );
    let absolute_delta = (nlo_value - published_nlo_pb).abs();
    let absolute_tolerance = (3.0 * nlo_error).max(0.05 * published_nlo_pb);
    assert!(
        absolute_delta <= absolute_tolerance,
        "explicit-local-3D absolute NLO normalization mismatch: |NLO|={nlo_value:e} ± {nlo_error:e} pb, converted published value={published_nlo_pb:e} pb, |delta|={absolute_delta:e}, tolerance={absolute_tolerance:e}",
    );

    clean_test(&test_root);
    Ok(())
}
