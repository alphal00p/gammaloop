use super::*;
use gammalooprs::{
    cff::surface::HybridSurfaceID,
    integrands::evaluation::StabilityStatus,
    processes::{CutId, ProcessCollection},
};
use std::{collections::BTreeSet, f64::consts::PI};

#[test]
fn thermal_vector_current_lo_matches_fermi_phase_space() -> Result<()> {
    let _serial_guard = gammaloop_integration_tests::SERIAL_TEST_LOCK.lock();
    let root = get_tests_workspace_path().join(format!("thermal_vector_lo_{}", std::process::id()));
    let mut cli = get_example_cli(
        "thermal_rates/qcd_vector/run.toml",
        &["model", "lo"],
        Some(root),
        None,
        true,
    )?;
    cli.run_command("set process -p vector_lo -i LO kv sampling.sampling_multichanneling=false")?;
    let (process_id, integrand_name) = cli.state.find_integrand_ref(
        Some(&ProcessRef::Unqualified("vector_lo".into())),
        Some(&"LO".into()),
    )?;
    let ProcessCollection::CrossSections(cross_sections) =
        &cli.state.process_list.processes[process_id].collection
    else {
        unreachable!()
    };
    let graph = &cross_sections["LO"].supergraphs[0];
    let routes = graph
        .graph
        .iter_edges_of(&graph.cuts[CutId(0)].cut)
        .map(|(_, edge, _)| {
            let signature = &graph.graph.loop_momentum_basis.edge_signatures[edge];
            (
                *signature.internal.iter().next().unwrap() as i8 as f64,
                signature
                    .external
                    .iter()
                    .map(|sign| *sign as i8 as f64)
                    .sum::<f64>(),
            )
        })
        .collect_vec();
    assert_eq!(routes.len(), 2);
    let electric = cli.state.model.get_parameter("ee").value.unwrap().re.0;
    let strong = cli.state.model.get_parameter("G").value.unwrap().re.0;
    assert!((electric - 1.0).abs() < 1e-12 && (strong - 1.0).abs() < 1e-12);
    let points = arr2(&[[0.6_f64, 0.0, 0.8], [0.3, -0.4, -0.7]]);
    let discrete = ndarray::Array2::<usize>::zeros((points.nrows(), 0));
    let batch = EvaluateSamples {
        process_id: Some(process_id),
        integrand_name: Some(integrand_name.clone()),
        use_arb_prec: false,
        minimal_output: false,
        return_generated_events: Some(true),
        momentum_space: true,
        points: points.view(),
        integrator_weights: None,
        discrete_dims: Some(discrete.view()),
        graph_names: None,
        orientations: None,
    }
    .run(&mut cli.state)?;
    for (point, sample) in points.rows().into_iter().zip(batch.samples) {
        let energies = |t: f64| {
            routes
                .iter()
                .map(|(sign, shift)| {
                    ((t * sign * point[0]).powi(2)
                        + (t * sign * point[1]).powi(2)
                        + (t * sign * point[2] + shift).powi(2))
                    .sqrt()
                })
                .collect_vec()
        };
        let (mut low, mut high) = (0.0, 1.0);
        while energies(high).iter().sum::<f64>() < 4.0 {
            high *= 2.0;
        }
        for _ in 0..80 {
            let mid = (low + high) / 2.0;
            if energies(mid).iter().sum::<f64>() < 4.0 {
                low = mid;
            } else {
                high = mid;
            }
        }
        let t = (low + high) / 2.0;
        let es = energies(t);
        let derivative = routes
            .iter()
            .zip(&es)
            .map(|((sign, shift), energy)| {
                (t * point.iter().map(|x| x * x).sum::<f64>() + sign * shift * point[2]) / energy
            })
            .sum::<f64>();
        // -g_mu_nu times the native spin/color-summed current gives
        // |M|²=4 Nc e² Qd² K²=20. No initial polarization average.
        let weight = 1.0 - 1.0 / (es[0].exp() + 1.0) - 1.0 / (es[1].exp() + 1.0);
        let h = 2.0 / PI.sqrt() * (2.0 - t * t - 1.0 / (t * t)).exp();
        let expected = 20.0 * 2.0 * PI * t.powi(3) * h * weight
            / (derivative * 4.0 * es[0] * es[1] * (2.0 * PI).powi(3) * 8.0);
        let evaluation = sample.evaluation;
        assert!(!evaluation.evaluation_metadata.as_ref().unwrap().is_nan);
        let actual = evaluation.integrand_result;
        assert!(
            actual.im.0.abs() < 1e-12 && (actual.re.0 / expected - 1.0).abs() < 1e-9,
            "vector LO {point:?}: {actual:?}, independent phase space={expected:e}"
        );
        let events = evaluation
            .event_groups
            .iter()
            .flat_map(|group| group.iter())
            .collect_vec();
        assert_eq!(events.len(), 1);
        assert!((events[0].weight.re.0 / expected - 1.0).abs() < 1e-9);
        let mut pdgs = events[0]
            .cut_info
            .particle_pdgs
            .1
            .iter()
            .copied()
            .collect_vec();
        pdgs.sort_unstable();
        assert_eq!(pdgs, [-1, 1]);
    }

    // Jackson arXiv:1910.07552v2, vector trace: the analytic one-loop
    // rho includes Nc=3 but not e²Qd². Native flux gives rho=36 Gamma_rho.
    let target_rho = 45.0 / (4.0 * PI)
        * (1.0 + 2.0 * ((1.0 + (-2.5_f64).exp()) / (1.0 + (-1.5_f64).exp())).ln());
    let mut estimates = Vec::new();
    for count in [16, 32] {
        let mut coordinates = Vec::new();
        let mut weights = Vec::new();
        for i in 0..2 * count {
            let x = (i as f64 + 0.5) / (2 * count) as f64;
            let radius = 2.0 * x / (1.0 - x);
            for j in 0..count {
                let cosine = 2.0 * (j as f64 + 0.5) / count as f64 - 1.0;
                coordinates.extend([
                    radius * (1.0 - cosine * cosine).sqrt(),
                    0.0,
                    radius * cosine,
                ]);
                weights.push(
                    8.0 * PI * radius.powi(2) / ((1.0 - x).powi(2) * (2 * count * count) as f64),
                );
            }
        }
        let points = ndarray::Array2::from_shape_vec((weights.len(), 3), coordinates)?;
        let discrete = ndarray::Array2::<usize>::zeros((points.nrows(), 0));
        let batch = EvaluateSamples {
            process_id: Some(process_id),
            integrand_name: Some(integrand_name.clone()),
            use_arb_prec: false,
            minimal_output: false,
            return_generated_events: Some(false),
            momentum_space: true,
            points: points.view(),
            integrator_weights: None,
            discrete_dims: Some(discrete.view()),
            graph_names: None,
            orientations: None,
        }
        .run(&mut cli.state)?;
        estimates.push(
            batch
                .samples
                .iter()
                .zip(weights)
                .map(|(sample, weight)| {
                    assert!(
                        !sample
                            .evaluation
                            .evaluation_metadata
                            .as_ref()
                            .unwrap()
                            .is_nan
                    );
                    sample.evaluation.integrand_result.re.0 * weight * 36.0
                })
                .sum::<f64>(),
        );
    }
    let extrapolated = (4.0 * estimates[1] - estimates[0]) / 3.0;
    gammalooprs::debug_tags!(#integration,#thermal,#benchmark,#summary;
        stage="thermal_rate_reference_comparison",fixture="qcd_vector",observable="rho_lo",
        coarse=estimates[0],fine=estimates[1],extrapolated,target=target_rho,
        relative_error=(extrapolated/target_rho-1.0).abs(),
        "Native vector-current LO compared with the analytic Fermi phase-space integral"
    );
    assert!(
        (extrapolated / target_rho - 1.0).abs() < 1e-5,
        "thermal vector LO quadrature {estimates:?}, extrapolated {extrapolated:e}, target {target_rho:e}"
    );
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn thermal_vector_current_nlo_retains_real_virtual_and_raised_cuts() -> Result<()> {
    let _serial_guard = gammaloop_integration_tests::SERIAL_TEST_LOCK.lock();
    let root =
        get_tests_workspace_path().join(format!("thermal_vector_nlo_{}", std::process::id()));
    let mut cli = get_example_cli(
        "thermal_rates/qcd_vector/run.toml",
        &["model", "nlo"],
        Some(root),
        None,
        true,
    )?;
    let (process_id, integrand_name) = cli.state.find_integrand_ref(
        Some(&ProcessRef::Unqualified("vector_nlo".into())),
        Some(&"NLO".into()),
    )?;
    let info = cli
        .state
        .get_integrand_info(Some(&ProcessRef::Id(process_id)), Some(&integrand_name))?;
    let masters = info
        .graph_groups
        .iter()
        .flat_map(|group| group.graphs.iter())
        .filter(|graph| graph.is_master)
        .map(|graph| graph.name.as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(masters, BTreeSet::from(["GL0", "GL2"]));
    let ProcessCollection::CrossSections(cross_sections) =
        &cli.state.process_list.processes[process_id].collection
    else {
        unreachable!()
    };
    let mut cut_multiplicities = BTreeSet::new();
    let mut raised = false;
    let mut crossed = false;
    for graph in &cross_sections["NLO"].supergraphs {
        for cut in &graph.cuts {
            cut_multiplicities.insert(graph.graph.iter_edges_of(&cut.cut).count());
        }
        for group in &graph.derived_data.cut_group_data.cut_groups {
            raised |= group.residue.max_occurrence() > 1;
            crossed |= group
                .residue
                .surface_ids
                .iter()
                .any(|surface| matches!(surface, HybridSurfaceID::Hsurface(_)));
        }
    }
    assert!(
        cut_multiplicities.contains(&2) && cut_multiplicities.contains(&3),
        "inclusive NLO requires two-particle virtual and three-particle real cuts"
    );
    assert!(
        raised && crossed,
        "thermal self-energy poles and crossed bath channels must remain present"
    );
    assert!(cli.cli_settings.global.generation.uv.subtract_uv);
    assert!(cli.cli_settings.global.generation.uv.generate_integrated);
    assert!(
        !cli.cli_settings
            .global
            .generation
            .uv
            .local_uv_cts_from_expanded_4d_integrands
    );

    // One ordinary point exercises the complete generated graph sum and each
    // master separately. Scale changes are finite point checks: their physical
    // invariance is an integrated property tested by the example command card.
    let point = [0.17, -0.31, 0.43, -0.29, 0.53, 0.71];
    let points = ndarray::Array2::from_shape_vec((3, 6), point.repeat(3))?;
    let mut baseline = None;
    for scale in [4.0_f64, 2.0, 4.0] {
        cli.run_command(&format!("set process -p vector_nlo -i NLO kv general.m_uv={scale} general.renormalization_localization_scale={scale}"))?;
        let batch = EvaluateSamples {
            process_id: Some(process_id),
            integrand_name: Some(integrand_name.clone()),
            use_arb_prec: false,
            minimal_output: false,
            return_generated_events: Some(true),
            momentum_space: true,
            points: points.view(),
            integrator_weights: None,
            discrete_dims: None,
            graph_names: Some(vec![None, Some("GL0".into()), Some("GL2".into())]),
            orientations: None,
        }
        .run(&mut cli.state)?;
        let values = batch
            .samples
            .iter()
            .map(|sample| sample.evaluation.integrand_result)
            .collect_vec();
        for (graph, sample) in ["inclusive", "GL0", "GL2"].into_iter().zip(&batch.samples) {
            let evaluation = &sample.evaluation;
            let metadata = evaluation.evaluation_metadata.as_ref().unwrap();
            assert!(
                !metadata.is_nan,
                "{graph}/scale{scale}: discarded nonfinite point"
            );
            assert!(
                !matches!(
                    metadata
                        .stability_results
                        .last()
                        .map(|result| &result.status),
                    Some(StabilityStatus::Unstable(_))
                ),
                "{graph}/scale{scale}: unstable point"
            );
            let value = evaluation.integrand_result;
            assert!(value.re.0.is_finite() && value.im.0.is_finite() && value.re.0.abs() > 0.0);
            let mut event_sum = Complex::new(0.0, 0.0);
            let mut event_count = 0;
            for event in evaluation
                .event_groups
                .iter()
                .flat_map(|group| group.iter())
            {
                event_count += 1;
                event_sum += Complex::new(event.weight.re.0, event.weight.im.0);
                let mut residual = [0.0_f64; 4];
                for (sign, momenta) in [
                    (1.0, event.kinematic_configuration.0.as_slice()),
                    (-1.0, event.kinematic_configuration.1.as_slice()),
                ] {
                    for p in momenta {
                        let components = [
                            p.temporal.value.0,
                            p.spatial.px.0,
                            p.spatial.py.0,
                            p.spatial.pz.0,
                        ];
                        assert!(components[0] > 0.0);
                        for (sum, component) in residual.iter_mut().zip(components) {
                            *sum += sign * component;
                        }
                    }
                }
                assert!(residual.iter().all(|component| component.abs() < 1e-9));
            }
            assert!(event_count > 0);
            let factor = evaluation.integrator_weight.0
                * evaluation.parameterization_jacobian.map_or(1.0, |j| j.0);
            let expected = Complex::new(value.re.0 * factor, value.im.0 * factor);
            assert!(
                (event_sum.re - expected.re).hypot(event_sum.im - expected.im)
                    < 1e-11 + 1e-9 * expected.re.hypot(expected.im),
                "{graph}/scale{scale}: event closure"
            );
            gammalooprs::debug_tags!(#evaluation,#thermal,#benchmark,#summary;
                stage="thermal_qcd_nlo_point",graph,auxiliary_scale=scale,
                value_re=value.re.0,value_im=value.im.0,event_count,
                "Native thermal QCD point with shared local and integrated ultraviolet subtraction"
            );
        }
        let graph_sum = values[1] + values[2];
        assert!(
            (values[0].re.0 - graph_sum.re.0).hypot(values[0].im.0 - graph_sum.im.0)
                < 1e-10 + 1e-8 * values[0].re.0.hypot(values[0].im.0)
        );
        if scale == 4.0 {
            if let Some(reference) = baseline {
                let reference: Complex<F<f64>> = reference;
                assert!(
                    (values[0].re.0 - reference.re.0).hypot(values[0].im.0 - reference.im.0)
                        < 1e-12 + 1e-10 * reference.re.0.hypot(reference.im.0)
                );
            } else {
                baseline = Some(values[0]);
                // The helper already saved the baseline. Verify rebuilding
                // its native evaluator backend preserves the same point.
                cli.state = gammaloop_api::state::State::load(
                    cli.cli_settings.state.folder.clone(),
                    None,
                    None,
                )?;
                cli.state.activate_loaded_integrand_backends(false)?;
            }
        }
    }
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
