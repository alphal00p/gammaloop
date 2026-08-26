use std::{collections::BTreeSet, f64::consts::PI};

use color_eyre::{Result, eyre::eyre};
use gammaloop_api::{
    commands::{
        Inspect, Integrate, Profile, integrate::RendererOption, profile::UltraVioletProfile,
    },
    state::ProcessRef,
};
use gammaloop_integration_tests::{clean_test, get_example_cli, get_tests_workspace_path};
use gammalooprs::uv::settings::{ApproximationType, FinalIntegrandDimension};
use serial_test::serial;

#[test]
#[serial]
fn epem_a_ttx_msbar_nlo_matches_the_published_inclusive_ratio_in_both_local_uv_routes() -> Result<()>
{
    const LO_PROCESS: &str = "epem_a_ttx_lo";
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
        .join("epem_a_ttx_msbar_nlo_matches_the_published_inclusive_ratio_in_both_local_uv_routes");
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
    }

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
    let (alpha_qed_inverse_re, alpha_qed_inverse_im) = model_parameters
        .data
        .get("aEWM1")
        .ok_or_else(|| eyre!("the effective LO model parameters do not contain aEWM1"))?;
    assert_eq!(
        alpha_qed_inverse_im.0, 0.0,
        "the acceptance expects a real electromagnetic coupling",
    );
    assert!(
        alpha_qed_inverse_re.0 > 0.0,
        "the acceptance expects a positive inverse electromagnetic coupling",
    );
    let alpha_qed = alpha_qed_inverse_re.0.recip();
    // Eq. (7.1) of arXiv:2203.11038 closes the external lepton trace;
    // Q_e^2=1 for this process and the generated cross-section is in pb.
    let gamma_star_to_epem_pb =
        2.0 * (4.0 * PI * alpha_qed) / (3.0 * E_CM.powi(3)) * GEV_SQUARED_TO_PICOBARN;
    let published_lo_pb = PUBLISHED_GAMMA_STAR_LO * gamma_star_to_epem_pb;
    let published_nlo_pb = PUBLISHED_GAMMA_STAR_NLO * gamma_star_to_epem_pb;

    for process in std::iter::once(LO_PROCESS).chain(ROUTES.map(|(process, _)| process)) {
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
                    lmb_channels = "monte_carlo"
                    lmb_channel_weight = "ose"
                '"#,
        ))?;
    }

    for (route_index, (process, _)) in ROUTES.into_iter().enumerate() {
        let analysis = Profile::UltraViolet(UltraVioletProfile {
            process: Some(ProcessRef::Unqualified(process.to_string())),
            integrand_name: Some(NLO.to_string()),
            n_points: 8,
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

    let mut route_probe_phases = Vec::new();
    let mut route_probes = Vec::new();
    for (process, _) in ROUTES {
        // Compare the complete orientation-local and explicit-orientation-sum
        // representations before spending time on either numerical integral.
        cli.run_command(&format!(
            r#"set process -p {process} -i {NLO} kv
                sampling.graphs="summed"
                sampling.orientations="summed"
                sampling.lmb_channels="summed""#,
        ))?;
        let (_, nlo_probe) = Inspect {
            process: Some(ProcessRef::Unqualified(process.to_string())),
            integrand_name: Some(NLO.to_string()),
            point: vec![0.11, 0.23, 0.37, 0.41, 0.53, 0.67],
            ..Default::default()
        }
        .run(&mut cli.state)?;
        assert!(
            nlo_probe.re.is_finite()
                && nlo_probe.im.is_finite()
                && nlo_probe.re.abs().max(nlo_probe.im.abs()) > 0.0,
            "the {process} NLO pointwise probe must be finite and nonzero, got {nlo_probe:e}",
        );
        assert!(
            nlo_probe.re.abs().min(nlo_probe.im.abs())
                <= 1.0e-8 * nlo_probe.re.abs().max(nlo_probe.im.abs()),
            "the {process} NLO probe is not phase-pure enough for a magnitude acceptance: {nlo_probe:e}",
        );
        route_probe_phases.push(if nlo_probe.re.abs() >= nlo_probe.im.abs() {
            "real"
        } else {
            "imag"
        });
        route_probes.push((process, nlo_probe.re, nlo_probe.im));
    }
    let [
        (first_probe_process, first_probe_re, first_probe_im),
        (second_probe_process, second_probe_re, second_probe_im),
    ] = route_probes.as_slice()
    else {
        unreachable!("both local UV routes are always probed")
    };
    let pointwise_scale = first_probe_re
        .hypot(*first_probe_im)
        .max(second_probe_re.hypot(*second_probe_im));
    let pointwise_delta =
        (first_probe_re - second_probe_re).hypot(first_probe_im - second_probe_im);
    assert!(
        pointwise_delta <= 1.0e-9 * pointwise_scale,
        "summed-orientation local UV route mismatch at the identical point: {first_probe_process}=({first_probe_re:e},{first_probe_im:e}i), {second_probe_process}=({second_probe_re:e},{second_probe_im:e}i), relative delta={:e}",
        pointwise_delta / pointwise_scale,
    );

    // Regenerate the direct local-3D route in the same selector-free mode
    // required by the projected local-4D route. This reuses its later Monte
    // Carlo slot while the probe above retains coverage of orientation-local
    // direct-3D generation.
    {
        let generation = &mut cli.cli_settings.global.generation;
        generation.explicit_orientation_sum_only = true;
        generation.uv.local_uv_cts_from_expanded_4d_integrands = false;
    }
    cli.run_command(&format!("generate existing -p {} -i {NLO}", ROUTES[0].0))?;
    cli.run_command(&format!(
        r#"set process -p {} -i {NLO} kv
            sampling.graphs="summed"
            sampling.orientations="summed"
            sampling.lmb_channels="summed""#,
        ROUTES[0].0,
    ))?;
    let (_, explicit_direct_probe) = Inspect {
        process: Some(ProcessRef::Unqualified(ROUTES[0].0.to_string())),
        integrand_name: Some(NLO.to_string()),
        point: vec![0.11, 0.23, 0.37, 0.41, 0.53, 0.67],
        ..Default::default()
    }
    .run(&mut cli.state)?;
    let explicit_direct_delta = (first_probe_re - explicit_direct_probe.re)
        .hypot(first_probe_im - explicit_direct_probe.im);
    let explicit_direct_scale = first_probe_re
        .hypot(*first_probe_im)
        .max(explicit_direct_probe.re.hypot(explicit_direct_probe.im));
    assert!(
        explicit_direct_scale > 0.0 && explicit_direct_delta <= 1.0e-9 * explicit_direct_scale,
        "direct local-3D orientation-local versus explicit-sum mismatch at the identical point: orientation-local=({first_probe_re:e},{first_probe_im:e}i), explicit-sum=({:e},{:e}i), relative delta={:e}",
        explicit_direct_probe.re,
        explicit_direct_probe.im,
        explicit_direct_delta / explicit_direct_scale,
    );

    for graph_id in 0..2 {
        let mut graph_probes = Vec::new();
        for (process, _) in ROUTES {
            let (_, probe) = Inspect {
                process: Some(ProcessRef::Unqualified(process.to_string())),
                integrand_name: Some(NLO.to_string()),
                point: vec![17.0, -31.0, 43.0, -29.0, 53.0, 71.0],
                momentum_space: true,
                graph_id: Some(graph_id),
                ..Default::default()
            }
            .run(&mut cli.state)?;
            assert!(
                probe.re.is_finite() && probe.im.is_finite(),
                "the {process} graph #{graph_id} momentum-space probe must be finite, got {probe:e}",
            );
            graph_probes.push((process, probe.re, probe.im));
        }
        let [
            (first_process, first_re, first_im),
            (second_process, second_re, second_im),
        ] = graph_probes.as_slice()
        else {
            unreachable!("both local UV routes are always probed")
        };
        let scale = first_re.hypot(*first_im).max(second_re.hypot(*second_im));
        let delta = (first_re - second_re).hypot(first_im - second_im);
        assert!(
            scale > 0.0 && delta <= 1.0e-9 * scale,
            "summed-orientation graph #{graph_id} local UV route mismatch at the identical momentum-space point: {first_process}=({first_re:e},{first_im:e}i), {second_process}=({second_re:e},{second_im:e}i), relative delta={:e}",
            delta / scale,
        );
    }
    for (process, _) in ROUTES {
        cli.run_command(&format!(
            r#"set process -p {process} -i {NLO} kv
                sampling.graphs="monte_carlo"
                sampling.orientations="summed"
                sampling.lmb_channels="monte_carlo""#,
        ))?;
    }

    let (_, lo_probe) = Inspect {
        process: Some(ProcessRef::Unqualified(LO_PROCESS.to_string())),
        integrand_name: Some(LO.to_string()),
        point: vec![0.19, 0.43, 0.71],
        discrete_dim: vec![0, 0],
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
                integrator.min_samples_for_update=10000
                integrator.n_start=10000
                integrator.n_increase=10000
                integrator.n_max=30000
                integrator.target_relative_accuracy=0.01
                integrator.seed=7331
                integrator.discrete_dim_learning_rate=0.0
                integrator.continuous_dim_learning_rate=0.25"#,
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

    let mut route_estimates = Vec::new();
    for (route_index, (process, _)) in ROUTES.into_iter().enumerate() {
        let process_ref = ProcessRef::Unqualified(process.to_string());
        let nlo_phase = route_probe_phases[route_index];
        cli.run_command(&format!(
            r#"set process -p {process} -i {NLO} kv
                    integrator.integrated_phase="{nlo_phase}"
                    integrator.min_samples_for_update=100000
                    integrator.n_start=100000
                    integrator.n_increase=100000
                    integrator.n_max=200000
                    integrator.target_relative_accuracy=1.0e-12
                    integrator.seed={}
                    integrator.discrete_dim_learning_rate=0.0
                    integrator.continuous_dim_learning_rate=0.25"#,
            1337 + route_index,
        ))?;

        let nlo_output = Integrate {
            process: vec![process_ref],
            integrand_name: vec![NLO.to_string()],
            n_cores: Some(10),
            workspace_path: Some(test_root.join(format!("{process}_integration_workspace"))),
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
            .ok_or_else(|| eyre!("expected one {process} NLO integration slot"))?;
        let nlo_value = nlo_estimate.result.re.0.hypot(nlo_estimate.result.im.0);
        let nlo_error = nlo_estimate.error.re.0.hypot(nlo_estimate.error.im.0);
        let nlo_relative_error = nlo_error / nlo_value;
        assert!(
            nlo_relative_error.is_finite() && nlo_relative_error <= 0.12,
            "{process} NLO uncertainty is too large for the normalization acceptance: {nlo_value:e} ± {nlo_error:e} ({:.2}%)",
            100.0 * nlo_relative_error,
        );

        let expected_nlo = PUBLISHED_NLO_OVER_LO * lo_value;
        let combined_error = nlo_error.hypot(PUBLISHED_NLO_OVER_LO * lo_error);
        let delta = (nlo_value - expected_nlo).abs();
        assert!(
            combined_error.is_finite() && delta <= 3.0 * combined_error,
            "{process} inclusive NLO normalization mismatch: |NLO|={nlo_value:e} ± {nlo_error:e}, published-ratio*LO={expected_nlo:e} ± {:e}, |delta|={delta:e}, delta/sigma={:e}",
            PUBLISHED_NLO_OVER_LO * lo_error,
            delta / combined_error,
        );
        let absolute_delta = (nlo_value - published_nlo_pb).abs();
        let absolute_tolerance = (3.0 * nlo_error).max(0.05 * published_nlo_pb);
        assert!(
            absolute_delta <= absolute_tolerance,
            "{process} absolute NLO normalization mismatch: |NLO|={nlo_value:e} ± {nlo_error:e} pb, converted published value={published_nlo_pb:e} pb, |delta|={absolute_delta:e}, tolerance={absolute_tolerance:e}",
        );
        route_estimates.push((process, nlo_value, nlo_error));
    }

    let [
        (first_process, first_value, first_error),
        (second_process, second_value, second_error),
    ] = route_estimates.as_slice()
    else {
        unreachable!("the acceptance defines exactly two local-UV routes")
    };
    let route_delta = (first_value - second_value).abs();
    let route_error = first_error.hypot(*second_error);
    assert!(
        route_delta <= 3.0 * route_error,
        "local-UV route mismatch: {first_process}={first_value:e} ± {first_error:e} pb, {second_process}={second_value:e} ± {second_error:e} pb, |delta|={route_delta:e}, delta/sigma={:e}",
        route_delta / route_error,
    );

    clean_test(&test_root);
    Ok(())
}
