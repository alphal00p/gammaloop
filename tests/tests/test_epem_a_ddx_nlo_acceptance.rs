use std::{
    collections::{BTreeMap, BTreeSet},
    f64::consts::PI,
};

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
    cff::CffEnergyBoundSourceKind,
    integrands::evaluation::PreciseEvaluationResultOutput,
    utils::{ArbPrec, F, FloatLike},
    uv::settings::{ApproximationType, FinalIntegrandDimension},
};
use ndarray::Array2;
use serial_test::serial;
use spenso::algebra::algebraic_traits::IsZero;
use symbolica::domains::float::Real;

#[test]
#[serial]
fn epem_a_ddx_nlo_is_alpha_s_over_pi_times_lo_in_all_local_uv_routes() -> Result<()> {
    const LO_PROCESS: &str = "epem_a_ddx_lo";
    const EXPLICIT_3D_PROCESS: &str = "epem_a_ddx_explicit_local_3d";
    const LOCAL_ONLY_ROUTES: [(&str, bool, bool); 3] = [
        ("epem_a_ddx_local_only_3d", false, false),
        ("epem_a_ddx_local_only_explicit_3d", true, false),
        ("epem_a_ddx_local_only_4d_then_cff", true, true),
    ];
    const NLO: &str = "NLO";
    const LO: &str = "LO";
    const E_CM: f64 = 400.0;
    const GEV_SQUARED_TO_PICOBARN: f64 = 0.389379304e9;
    const PUBLISHED_GAMMA_STAR_LO: f64 = 5.031049e-1;
    const PUBLISHED_GAMMA_STAR_NLO: f64 = 5.03926e-2 - 3.14956e-2;
    const COMMON_UV_SCALE: f64 = 1000.0;
    const SHIFTED_SCALE: f64 = 400.0;
    const MUV_SHIFTED: &str = "NLO_muv_shifted";
    const MUR_SHIFTED: &str = "NLO_mur_shifted";
    const RHO_SHIFTED: &str = "NLO_rho_shifted";
    const ROUTES: [(&str, bool); 2] = [
        ("epem_a_ddx_local_3d", false),
        ("epem_a_ddx_local_4d_then_cff", true),
    ];

    let test_root = get_tests_workspace_path()
        .join("epem_a_ddx_nlo_is_alpha_s_over_pi_times_lo_in_all_local_uv_routes");
    clean_test(&test_root);

    // Use the production NLO example settings through integrand generation, but keep this
    // acceptance inclusive and independent from its semi-inclusive integration card.
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
        uv.subtract_uv = true;
        uv.generate_integrated = true;
        uv.final_integrand = FinalIntegrandDimension::ThreeD;
        uv.local_uv_cts_from_expanded_4d_integrands = false;
        uv.renormalization_prescription.log_divergent = ApproximationType::MUV;
        uv.renormalization_prescription.massive_power_divergent = ApproximationType::MUV;
        uv.renormalization_prescription.massless_power_divergent = ApproximationType::MUV;
        uv.renormalization_prescription.overrides.clear();
    }
    cli.default_runtime_settings
        .subtraction
        .disable_threshold_subtraction = false;
    cli.run_command("import model sm-default.json")?;
    cli.run_command("run set_model_parameters")?;
    cli.run_command(
        r#"generate xs e+ e- > d d~ | e+ e- d d~ g ghG ghG~ a QCD^2==0 QED^2==4 [{{1}} QCD=0]
            --numerator-grouping group_identical_graphs_up_to_scalar_rescaling
            --symmetrize-left-right-states true
            -p epem_a_ddx_lo -i LO --only-diagrams"#,
    )?;
    cli.run_command("generate existing -p epem_a_ddx_lo -i LO")?;

    for (process, project_local_4d_to_cff) in ROUTES {
        {
            let generation = &mut cli.cli_settings.global.generation;
            generation.explicit_orientation_sum_only = project_local_4d_to_cff;
            generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d_to_cff;
        }
        cli.run_command(
            &format!(
                r#"generate xs e+ e- > d d~ | e+ e- d d~ g ghG ghG~ a QCD^2==2 QED^2==4 [{{{{2}}}} QCD=1]
                --numerator-grouping group_identical_graphs_up_to_scalar_rescaling
                --symmetrize-left-right-states true
                -p {process} -i NLO --global-prefactor-num "-1𝑖" --only-diagrams"#,
            ),
        )?;
        cli.run_command(&format!("generate existing -p {process} -i NLO"))?;

        let nlo_process_ref = ProcessRef::Unqualified(process.to_string());
        let nlo_name = NLO.to_string();
        let nlo_info = cli
            .state
            .get_integrand_info(Some(&nlo_process_ref), Some(&nlo_name))?;
        let nlo_masters = nlo_info
            .graph_groups
            .iter()
            .flat_map(|group| group.graphs.iter())
            .filter(|graph| graph.is_master)
            .map(|graph| graph.name.clone())
            .collect::<BTreeSet<_>>();
        assert_eq!(
            nlo_masters,
            BTreeSet::from(["GL0".to_string(), "GL2".to_string()]),
            "the {process} acceptance must integrate the complete GL0+GL2 NLO example",
        );

        let (process_id, resolved_nlo_name) = cli
            .state
            .find_integrand_ref(Some(&nlo_process_ref), Some(&nlo_name))?;
        let reports = &cli
            .state
            .generation_summary(process_id, &resolved_nlo_name)
            .ok_or_else(|| eyre!("missing generation summary for {process}"))?
            .reports;
        let gl0_reports = &reports
            .iter()
            .find(|report| report.graph_name == "GL0")
            .ok_or_else(|| eyre!("missing GL0 generation report for {process}"))?
            .stats
            .cff_energy_degree_bound_reports;
        let expected_contracted_gl0_bounds = [(2, 1), (3, 2), (5, 1)];
        assert!(
            gl0_reports.iter().any(|report| {
                report.source_kind == CffEnergyBoundSourceKind::PhysicalGraph
                    && report.physical_parent_bounds == expected_contracted_gl0_bounds
                    && report.assigned_cff_source_bounds == expected_contracted_gl0_bounds
            }),
            "{process} must retain the contracted physical-graph GL0 source with exactly its degree-two edge-3 bound: {gl0_reports:?}",
        );
        if project_local_4d_to_cff {
            let expected_exact_gl0_parents = [(7, 1), (8, 2)];
            assert!(
                gl0_reports.iter().any(|report| {
                    report.source_kind == CffEnergyBoundSourceKind::ExactFourD
                        && report.physical_parent_bounds == expected_exact_gl0_parents
                        && report.assigned_cff_source_bounds.len() == 3
                        && report
                            .assigned_cff_source_bounds
                            .iter()
                            .all(|(_, degree)| *degree == 1)
                }),
                "{process} must dispatch the projected exact-4D GL0 parent bounds over exactly three degree-one occurrence-local energies: {gl0_reports:?}",
            );
        }
        assert!(
            gl0_reports.iter().all(|report| report.source_kind
                != CffEnergyBoundSourceKind::PhysicalGraph
                || report.physical_parent_bounds == report.assigned_cff_source_bounds),
            "{process} physical-graph sources must use their physical edge bounds unchanged: {gl0_reports:?}",
        );
        let gl2_reports = &reports
            .iter()
            .find(|report| report.graph_name == "GL2")
            .ok_or_else(|| eyre!("missing GL2 generation report for {process}"))?
            .stats
            .cff_energy_degree_bound_reports;
        let physical_gl2_max = gl2_reports
            .iter()
            .filter(|report| report.source_kind == CffEnergyBoundSourceKind::PhysicalGraph)
            .flat_map(|report| report.assigned_cff_source_bounds.iter())
            .map(|(_, degree)| *degree)
            .max();
        let exact_gl2_max = gl2_reports
            .iter()
            .filter(|report| report.source_kind == CffEnergyBoundSourceKind::ExactFourD)
            .flat_map(|report| report.assigned_cff_source_bounds.iter())
            .map(|(_, degree)| *degree)
            .max();
        assert_eq!(
            physical_gl2_max,
            Some(1),
            "{process} GL2 physical-graph CFF-source bounds must have maximum one: {gl2_reports:?}",
        );
        assert!(
            exact_gl2_max.is_none_or(|maximum| maximum == 1),
            "{process} any exact-4D GL2 source, including an integrated counterterm source in the direct-3D route, must have maximum one: {gl2_reports:?}",
        );
        if project_local_4d_to_cff {
            assert_eq!(
                exact_gl2_max,
                Some(1),
                "{process} projected local-4D GL2 must report an exact source with maximum bound one: {gl2_reports:?}",
            );
        }
    }

    let first_process_ref = ProcessRef::Unqualified(ROUTES[0].0.to_string());
    let nlo_name = NLO.to_string();
    let (process_id, resolved_nlo_name) = cli
        .state
        .find_integrand_ref(Some(&first_process_ref), Some(&nlo_name))?;
    let model_parameters = cli
        .state
        .resolve_effective_model_parameter_card_for_integrand(process_id, &resolved_nlo_name)?;
    let (alpha_s_re, alpha_s_im) = model_parameters
        .data
        .get("aS")
        .ok_or_else(|| eyre!("the effective NLO model parameters do not contain aS"))?;
    assert_eq!(
        alpha_s_im.0, 0.0,
        "the acceptance expects a real strong coupling"
    );
    let alpha_s = alpha_s_re.0;
    let (alpha_qed_inverse_re, alpha_qed_inverse_im) = model_parameters
        .data
        .get("aEWM1")
        .ok_or_else(|| eyre!("the effective NLO model parameters do not contain aEWM1"))?;
    assert_eq!(
        alpha_qed_inverse_im.0, 0.0,
        "the acceptance expects a real electromagnetic coupling",
    );
    assert!(
        alpha_qed_inverse_re.0 > 0.0,
        "the acceptance expects a positive inverse electromagnetic coupling",
    );
    let gamma_star_to_epem_pb =
        2.0 * (4.0 * PI / alpha_qed_inverse_re.0) / (3.0 * E_CM.powi(3)) * GEV_SQUARED_TO_PICOBARN;
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
        r#"generate xs e+ e- > d d~ | e+ e- d d~ g ghG ghG~ a QCD^2==2 QED^2==4 [{{{{2}}}} QCD=1]
            --numerator-grouping group_identical_graphs_up_to_scalar_rescaling
            --symmetrize-left-right-states true
            -p {EXPLICIT_3D_PROCESS} -i NLO --global-prefactor-num "-1𝑖" --only-diagrams"#,
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
            r#"generate xs e+ e- > d d~ | e+ e- d d~ g ghG ghG~ a QCD^2==2 QED^2==4 [{{{{2}}}} QCD=1]
                --numerator-grouping group_identical_graphs_up_to_scalar_rescaling
                --symmetrize-left-right-states true
                -p {process} -i NLO --global-prefactor-num "-1𝑖" --only-diagrams"#,
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
                    [200.0, 0.0, 0.0, 200.0],
                    [200.0, 0.0, 0.0, -200.0],
                ]
                helicities = ["summed_averaged", "summed_averaged"]

                [general]
                mu_r = {COMMON_UV_SCALE}
                m_uv = {COMMON_UV_SCALE}
                renormalization_localization_scale = 400.0

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
    let mut failed_profiles = Vec::new();
    let profile_routes = [ROUTES[0].0, EXPLICIT_3D_PROCESS, ROUTES[1].0];
    for (route_index, process) in profile_routes.into_iter().enumerate() {
        let analysis = Profile::UltraViolet(UltraVioletProfile {
            process: Some(ProcessRef::Unqualified(process.to_string())),
            integrand_name: Some(NLO.to_string()),
            n_points: 6,
            seed: Some(9100 + route_index as u64),
            ..Default::default()
        })
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
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
        if profile.failed != 0 {
            failed_profiles.push(format!(
                "the {process} default divergent-only UV profile failed:\n{profile}"
            ));
        }

        let selected_cut_analysis = Profile::UltraViolet(UltraVioletProfile {
            process: Some(ProcessRef::Unqualified(process.to_string())),
            integrand_name: Some(NLO.to_string()),
            graph: Some("GL0".to_string()),
            cutkosky_cut: vec![2, 5],
            n_points: 6,
            seed: Some(9200 + route_index as u64),
            ..Default::default()
        })
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
        let selected_cut_profile = selected_cut_analysis.pass_fail(-0.9);
        assert!(
            selected_cut_profile.total > 0,
            "the {process} selected GL0 cut [2,5] UV profile must test at least one limit",
        );
        if selected_cut_profile.failed != 0 {
            failed_profiles.push(format!(
                "the {process} selected GL0 cut [2,5] UV profile failed:\n{selected_cut_profile}"
            ));
        }
    }
    assert!(
        failed_profiles.is_empty(),
        "{}",
        failed_profiles.join("\n\n"),
    );

    let local_uv_processes = [ROUTES[0].0, EXPLICIT_3D_PROCESS, ROUTES[1].0];
    let local_only_processes = LOCAL_ONLY_ROUTES.map(|(process, _, _)| process);
    for process in local_uv_processes.into_iter().chain(local_only_processes) {
        cli.run_command(&format!(
            r#"set process -p {process} -i {NLO} kv
                sampling.graphs="summed"
                sampling.orientations="summed"
                sampling.lmb_multichanneling=true
                sampling.lmb_channels="summed"
                sampling.lmb_channel_weight="ose""#,
        ))?;
    }

    // A cheap f64 probe selects the integrated phase. Route equivalence and
    // every scale variant are compared graph by graph in Arb below.
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
        "the explicit-local-3D NLO probe is not phase-pure enough for a single-phase magnitude acceptance: {nlo_probe:e}",
    );
    let nlo_phase = if nlo_probe.re.abs() >= nlo_probe.im.abs() {
        "real"
    } else {
        "imag"
    };

    let scale_settings = [
        (MUV_SHIFTED, format!("general.m_uv={SHIFTED_SCALE}")),
        (MUR_SHIFTED, format!("general.mu_r={SHIFTED_SCALE}")),
        (
            RHO_SHIFTED,
            format!("general.renormalization_localization_scale={COMMON_UV_SCALE}"),
        ),
    ];

    // Reuse the explicit-local-3D integrand with three runtime-only scale changes. The common
    // M_UV=mu_r point is deliberately unequal to E_cm, so the normalization check below cannot
    // accidentally pass only at E_cm=M_UV=mu_r. The Arb comparison below independently verifies
    // all three representations at the central point before the integrated scale test.
    let scale_process = EXPLICIT_3D_PROCESS;
    for (variant, setting) in &scale_settings {
        cli.run_command(&format!(
            "duplicate integrand -p {scale_process} -i {NLO} --output_process_name {scale_process} --output_integrand_name {variant}",
        ))?;
        cli.run_command(&format!(
            "set process -p {scale_process} -i {variant} kv {setting}"
        ))?;
    }

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
    for (comparison, processes, disable_threshold_subtraction) in [
        ("integrated", local_uv_processes, false),
        ("local-only", local_only_processes, false),
        ("local-only-no-threshold", local_only_processes, true),
    ] {
        if disable_threshold_subtraction {
            for process in processes {
                cli.run_command(&format!(
                    "set process -p {process} -i {NLO} kv subtraction.disable_threshold_subtraction=true"
                ))?;
            }
        }
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
        "the LO probe is not phase-pure enough for a single-phase magnitude acceptance: {lo_probe:e}",
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
    let published_lo_delta = (lo_value - published_lo_pb).abs();
    let published_lo_tolerance = (3.0 * lo_error).max(0.05 * published_lo_pb);
    assert!(
        published_lo_delta <= published_lo_tolerance,
        "LO absolute normalization mismatch: |LO|={lo_value:e} ± {lo_error:e} pb, converted published gamma-star value={published_lo_pb:e} pb, |delta|={published_lo_delta:e}, tolerance={published_lo_tolerance:e}",
    );

    let normalization = alpha_s / PI;

    let scale_variants = [NLO, MUV_SHIFTED, MUR_SHIFTED, RHO_SHIFTED];
    for variant in scale_variants {
        cli.run_command(&format!(
            r#"set process -p {scale_process} -i {variant} kv
                sampling.graphs="monte_carlo"
                sampling.graph_names=[]
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
                integrator.seed=220311038
                integrator.discrete_dim_learning_rate=0.0
                integrator.continuous_dim_learning_rate=0.1"#,
        ))?;
    }

    let scale_process_ref = ProcessRef::Unqualified(scale_process.to_string());
    let output = Integrate {
        process: scale_variants.map(|_| scale_process_ref.clone()).to_vec(),
        integrand_name: scale_variants.map(str::to_string).to_vec(),
        n_cores: Some(10),
        workspace_path: Some(test_root.join("explicit_local_3d_scale_integration_workspace")),
        restart: true,
        renderer: RendererOption::Tabled,
        show_max_weight_info: false,
        no_stream_iterations: true,
        no_stream_updates: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    // Resolve every individual scale law from the graph-MC breakdown of one
    // correlated four-slot run. The orientations remain explicitly summed and
    // each LMB is sampled through the OSE multichannel configured above.
    let mut scale_estimates: BTreeMap<String, BTreeMap<String, (f64, f64)>> = BTreeMap::new();
    for variant in scale_variants {
        let slot = output
            .slot(&format!("{scale_process}@{variant}"))
            .ok_or_else(|| eyre!("expected {variant} integration slot"))?;
        let breakdown = match nlo_phase {
            "real" => slot.grid_breakdown.re.as_ref(),
            "imag" => slot.grid_breakdown.im.as_ref(),
            _ => unreachable!("the phase probe returns real or imag"),
        }
        .filter(|breakdown| breakdown.axis_label == "graph")
        .ok_or_else(|| eyre!("{variant} graph-MC result has no {nlo_phase} graph breakdown"))?;
        for entry in &breakdown.entries {
            let graph = entry
                .bin_label
                .clone()
                .unwrap_or_else(|| entry.bin_index.to_string());
            scale_estimates
                .entry(graph)
                .or_default()
                .insert(variant.to_string(), (entry.value.0, entry.error.0.abs()));
        }
    }
    let estimate = |graph: &str, variant: &str| -> Result<(f64, f64)> {
        scale_estimates
            .get(graph)
            .and_then(|estimates| estimates.get(variant))
            .copied()
            .ok_or_else(|| eyre!("missing {graph} {variant} scale estimate"))
    };
    let base_graphs = ["GL0", "GL2"]
        .into_iter()
        .map(|graph| Ok((graph.to_string(), estimate(graph, NLO)?)))
        .collect::<Result<BTreeMap<_, _>>>()?;

    let base_slot = output
        .slot(&format!("{scale_process}@{NLO}"))
        .ok_or_else(|| eyre!("expected the base NLO integration slot"))?;
    let (base_integral, base_integral_error) = match nlo_phase {
        "real" => (
            base_slot.integral.result.re.0,
            base_slot.integral.error.re.0.abs(),
        ),
        "imag" => (
            base_slot.integral.result.im.0,
            base_slot.integral.error.im.0.abs(),
        ),
        _ => unreachable!("the phase probe returns real or imag"),
    };
    let nlo_value = base_integral;
    let nlo_error = base_integral_error;
    assert!(
        nlo_value > 0.0,
        "the shipped -i NLO prefactor must produce a positive inclusive correction, got {nlo_value:e} ± {nlo_error:e}",
    );
    let nlo_relative_error = nlo_error / nlo_value;
    assert!(
        nlo_relative_error.is_finite() && nlo_relative_error <= 0.15,
        "explicit-local-3D non-E-equal-scale NLO uncertainty is too large for the normalization acceptance: {nlo_value:e} ± {nlo_error:e} ({:.2}%)",
        100.0 * nlo_relative_error,
    );
    let expected_nlo = normalization * lo_value;
    let combined_error = nlo_error.hypot(normalization * lo_error);
    let delta = (nlo_value - expected_nlo).abs();
    assert!(
        combined_error.is_finite() && delta <= 3.0 * combined_error,
        "explicit-local-3D inclusive NLO normalization mismatch at E_cm={E_CM}, M_UV=mu_r={COMMON_UV_SCALE}: NLO={nlo_value:e} ± {nlo_error:e}, (aS/pi)*LO={expected_nlo:e} ± {:e}, |delta|={delta:e}, delta/sigma={:e}",
        normalization * lo_error,
        delta / combined_error,
    );
    // The paper's individual gamma-star rows cannot be converted graph by
    // graph because GL0 and GL2 are not separately Ward-transverse. Their
    // inclusive sum is transverse and therefore obeys this closure factor.
    let published_delta = (nlo_value - published_nlo_pb).abs();
    let published_tolerance = (3.0 * nlo_error).max(0.05 * published_nlo_pb);
    assert!(
        published_delta <= published_tolerance,
        "explicit-local-3D absolute NLO normalization mismatch at E_cm={E_CM}, M_UV=mu_r={COMMON_UV_SCALE}: NLO={nlo_value:e} ± {nlo_error:e} pb, converted published gamma-star value={published_nlo_pb:e} pb, |delta|={published_delta:e}, tolerance={published_tolerance:e}",
    );
    let graph_sum = base_graphs.values().map(|(value, _)| value).sum::<f64>();
    let graph_sum_error = base_graphs
        .values()
        .map(|(_, error)| error * error)
        .sum::<f64>()
        .sqrt();
    assert!(
        (graph_sum - base_integral).abs() <= 4.0 * graph_sum_error.hypot(base_integral_error),
        "base graph-MC breakdown does not reproduce its slot integral: graphs={graph_sum:e} ± {graph_sum_error:e}, slot={base_integral:e} ± {base_integral_error:e}",
    );

    for (variant, scale_name) in [(MUV_SHIFTED, "M_UV"), (RHO_SHIFTED, "localization rho")] {
        for graph in ["GL0", "GL2"] {
            let (base_value, base_error) = base_graphs[graph];
            let (variant_value, variant_error) = estimate(graph, variant)?;
            let error = base_error.hypot(variant_error);
            let value_scale = base_value
                .abs()
                .max(variant_value.abs())
                .max(f64::MIN_POSITIVE);
            assert!(
                error.is_finite() && error / value_scale <= 0.20,
                "{graph} uncertainty is too large for the {scale_name}-invariance acceptance: base={base_value:e} ± {base_error:e}, shifted={variant_value:e} ± {variant_error:e}",
            );
            let delta = (base_value - variant_value).abs();
            assert!(
                delta <= 4.0 * error,
                "{graph} is not {scale_name} independent: base={base_value:e} ± {base_error:e}, shifted={variant_value:e} ± {variant_error:e}, delta/sigma={:e}",
                delta / error,
            );
        }
    }

    let mur_graphs = ["GL0", "GL2"]
        .into_iter()
        .map(|graph| Ok((graph.to_string(), estimate(graph, MUR_SHIFTED)?)))
        .collect::<Result<BTreeMap<_, _>>>()?;
    let base_total = base_graphs.values().map(|(value, _)| value).sum::<f64>();
    let mur_total = mur_graphs.values().map(|(value, _)| value).sum::<f64>();
    let mur_total_delta = (base_total - mur_total).abs();
    let mur_total_error = base_graphs
        .values()
        .chain(mur_graphs.values())
        .map(|(_, error)| error * error)
        .sum::<f64>()
        .sqrt();
    // Do not normalize this uncertainty by either total: both are small
    // differences of the individually resolved GL0 and GL2 contributions.
    // The graphwise logarithms below are the well-conditioned scale-law oracle.
    assert!(
        mur_total_delta <= 4.0 * mur_total_error,
        "GL0+GL2 is not mu_r independent across graphs: base={base_total:e}, shifted={mur_total:e}, delta/sigma={:e}",
        mur_total_delta / mur_total_error,
    );

    let expected_mur_shift = 2.0 / 3.0
        * normalization
        * lo_value
        * (COMMON_UV_SCALE.powi(2) / SHIFTED_SCALE.powi(2)).ln();
    let expected_mur_error = 2.0 / 3.0
        * normalization
        * lo_error
        * (COMMON_UV_SCALE.powi(2) / SHIFTED_SCALE.powi(2)).ln();
    for (graph, sign) in [("GL0", -1.0), ("GL2", 1.0)] {
        let observed_shift = base_graphs[graph].0 - mur_graphs[graph].0;
        let observed_error = base_graphs[graph].1.hypot(mur_graphs[graph].1);
        let expected = sign * expected_mur_shift;
        let combined_error = observed_error.hypot(expected_mur_error);
        let delta = (observed_shift - expected).abs();
        assert!(
            combined_error.is_finite() && combined_error / expected_mur_shift <= 0.55,
            "{graph} uncertainty is too large to resolve its individual mu_r logarithm: observed central-minus-shifted={observed_shift:e} ± {observed_error:e}, expected magnitude={expected_mur_shift:e}",
        );
        assert!(
            combined_error.is_finite() && delta <= 4.0 * combined_error,
            "{graph} has the wrong individual mu_r logarithm: observed central-minus-shifted={observed_shift:e} ± {observed_error:e}, expected={expected:e} ± {expected_mur_error:e}, delta/sigma={:e}",
            delta / combined_error,
        );
    }
    let graph_delta_sum = ["GL0", "GL2"]
        .iter()
        .map(|graph| base_graphs[*graph].0 - mur_graphs[*graph].0)
        .sum::<f64>();
    let graph_delta_error = ["GL0", "GL2"]
        .iter()
        .flat_map(|graph| [base_graphs[*graph].1, mur_graphs[*graph].1])
        .map(|error| error * error)
        .sum::<f64>()
        .sqrt();
    assert!(
        graph_delta_sum.abs() <= 4.0 * graph_delta_error,
        "the GL0 and GL2 mu_r shifts do not cancel: GL0 shift={:e}, GL2 shift={:e}, sum/sigma={:e}",
        base_graphs["GL0"].0 - mur_graphs["GL0"].0,
        base_graphs["GL2"].0 - mur_graphs["GL2"].0,
        graph_delta_sum.abs() / graph_delta_error,
    );

    clean_test(&test_root);
    Ok(())
}
