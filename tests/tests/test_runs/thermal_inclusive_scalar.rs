use super::*;
use gammalooprs::processes::ProcessCollection;
use std::collections::BTreeMap;

#[test]
fn thermal_virtual_scalar_cuts_preserve_channel_and_event_closure() -> Result<()> {
    let _serial_guard = gammaloop_integration_tests::SERIAL_TEST_LOCK.lock();
    let state_path =
        get_tests_workspace_path().join(format!("thermal_double_triangle_{}", std::process::id()));
    let mut cli = get_example_cli(
        "thermal_rates/massive_double_triangle/run.toml",
        &["generate"],
        Some(state_path),
        None,
        true,
    )?;
    cli.state =
        gammaloop_api::state::State::load(cli.cli_settings.state.folder.clone(), None, None)?;
    cli.state.activate_loaded_integrand_backends(false)?;
    for (block, historical) in [
        (
            "below_threshold",
            [
                5.137413496697824e-5,
                5.558131236202565e-5,
                -4.207177395047416e-6,
            ],
        ),
        (
            "above_threshold",
            [
                -9.58821293108844e-6,
                -6.330517695597147e-6,
                -3.2576952354912934e-6,
            ],
        ),
    ] {
        cli.run_command(&format!("run {block}"))?;
        let points = arr2(&[[0.4, 0.2, 0.1, -0.3, 0.5, 0.2]]);
        let discrete = ndarray::Array2::<usize>::zeros((1, 0));
        let mut values = BTreeMap::new();
        for name in ["thermal", "virtual", "real", "vacuum", "zero_temperature"] {
            let (process_id, integrand_name) = cli.state.find_integrand_ref(
                Some(&ProcessRef::Unqualified(name.into())),
                Some(&"NLO".into()),
            )?;
            let process = &cli.state.process_list.processes[process_id];
            let ProcessCollection::CrossSections(cross_sections) = &process.collection else {
                unreachable!()
            };
            let graph = &cross_sections["NLO"].supergraphs[0];
            assert!(
                graph
                    .derived_data
                    .cut_group_data
                    .cut_groups
                    .iter()
                    .all(|group| group.residue.max_occurrence() <= 1)
            );
            assert!(
                !process
                    .settings_history
                    .as_ref()
                    .unwrap()
                    .generation
                    .threshold_subtraction
                    .enable_thresholds
            );
            let result = EvaluateSamples {
                process_id: Some(process_id),
                integrand_name: Some(integrand_name),
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
            let evaluation = &result.samples[0].evaluation;
            let value = evaluation.integrand_result;
            assert!(value.re.0.is_finite() && value.im.0.abs() < 1e-12);
            let mut event_sum = 0.0;
            let mut event_count = 0;
            for event in evaluation
                .event_groups
                .iter()
                .flat_map(|group| group.iter())
            {
                event_count += 1;
                event_sum += event.weight.re.0;
                let mut residual = [0.0_f64; 4];
                for (sign, momenta) in [
                    (1.0, event.kinematic_configuration.0.as_slice()),
                    (-1.0, event.kinematic_configuration.1.as_slice()),
                ] {
                    for p in momenta {
                        assert!(p.temporal.value.0 > 0.0);
                        for (sum, component) in residual.iter_mut().zip([
                            p.temporal.value.0,
                            p.spatial.px.0,
                            p.spatial.py.0,
                            p.spatial.pz.0,
                        ]) {
                            *sum += sign * component;
                        }
                    }
                }
                assert!(residual.iter().all(|component| component.abs() < 1e-10));
            }
            assert!(event_count > 0);
            assert!((event_sum - value.re.0).abs() < 1e-12 + 1e-10 * value.re.0.abs());
            values.insert(name, value.re.0);
        }
        assert!((values["thermal"] - values["virtual"] - values["real"]).abs() < 1e-12);
        assert!((values["zero_temperature"] - values["vacuum"]).abs() < 1e-12);
        // These are historical pointwise baselines for the same graph, routing,
        // and LU profile. Independent integrated literature targets are checked
        // by the separately recorded native acceptance runs, not by these points.
        for (name, expected) in ["thermal", "virtual", "real"].into_iter().zip(historical) {
            gammalooprs::debug_tags!(#integration, #thermal, #benchmark, #summary;
                stage="thermal_inclusive_point_comparison", fixture="massive_double_triangle",
                block, process=name, value=values[name], historical=expected,
                relative_difference=(values[name]-expected).abs()/expected.abs(),
                "Restored scalar cut compared with its retained pointwise baseline"
            );
            assert!(
                (values[name] / expected - 1.0).abs() < 1e-8,
                "{block}/{name}: {} != {expected}",
                values[name]
            );
        }
    }
    clean_test(cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn thermal_raised_scalar_uv_addback_preserves_channel_closure() -> Result<()> {
    let _serial_guard = gammaloop_integration_tests::SERIAL_TEST_LOCK.lock();
    let state_path =
        get_tests_workspace_path().join(format!("thermal_self_energy_uv_{}", std::process::id()));
    let mut cli = get_example_cli(
        "thermal_rates/massive_self_energy/run.toml",
        &["generate"],
        Some(state_path),
        None,
        true,
    )?;
    cli.state =
        gammaloop_api::state::State::load(cli.cli_settings.state.folder.clone(), None, None)?;
    cli.state.activate_loaded_integrand_backends(false)?;
    let points = arr2(&[
        [0.4, 0.2, 0.1, -0.3, 0.5, 0.2],
        [0.7, -0.3, 0.2, 0.1, 0.6, -0.4],
    ]);
    let discrete = ndarray::Array2::<usize>::zeros((points.nrows(), 0));
    let mut values = BTreeMap::new();
    for name in [
        "thermal",
        "vacuum",
        "zero_temperature",
        "thermal_hh",
        "thermal_hhl",
        "thermal_aux",
        "thermal_bare",
    ] {
        let (process_id, integrand_name) = cli.state.find_integrand_ref(
            Some(&ProcessRef::Unqualified(name.into())),
            Some(&"NLO".into()),
        )?;
        let process = &cli.state.process_list.processes[process_id];
        let generation = &process.settings_history.as_ref().unwrap().generation;
        assert!(!generation.medium.vacuum_subtraction);
        assert!(!generation.threshold_subtraction.enable_thresholds);
        assert_eq!(generation.uv.subtract_uv, name != "thermal_bare");
        assert_eq!(generation.uv.generate_integrated, name != "thermal_bare");
        let ProcessCollection::CrossSections(cross_sections) = &process.collection else {
            unreachable!()
        };
        let graph = &cross_sections["NLO"].supergraphs[0];
        assert_eq!(
            graph
                .derived_data
                .cut_group_data
                .cut_groups
                .iter()
                .any(|group| group.residue.max_occurrence() == 2),
            name != "thermal_hhl"
        );
        let result = EvaluateSamples {
            process_id: Some(process_id),
            integrand_name: Some(integrand_name),
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
        let mut actual = Vec::new();
        for sample in result.samples {
            let evaluation = sample.evaluation;
            let value = evaluation.integrand_result;
            assert!(value.re.0.is_finite() && value.im.0.abs() < 1e-11);
            let event_sum = evaluation
                .event_groups
                .iter()
                .flat_map(|group| group.iter())
                .map(|event| event.weight.re.0)
                .sum::<f64>();
            assert!(
                (event_sum - value.re.0).abs() < 1e-11 + 1e-9 * value.re.0.abs(),
                "{name}: event closure"
            );
            actual.push(value.re.0);
        }
        values.insert(name, actual);
    }
    for ((inclusive, hh), hhl) in values["thermal"]
        .iter()
        .zip(&values["thermal_hh"])
        .zip(&values["thermal_hhl"])
    {
        let parts = hh + hhl;
        assert!((inclusive - parts).abs() < 1e-11 + 1e-9 * inclusive.abs());
    }
    for (vacuum, cold) in values["vacuum"].iter().zip(&values["zero_temperature"]) {
        assert!((vacuum - cold).abs() < 1e-11);
    }
    assert_ne!(
        values["thermal"], values["thermal_bare"],
        "UV subtraction must alter the bare point"
    );
    assert_ne!(
        values["thermal"], values["thermal_aux"],
        "Auxiliary UV changes alter local terms before integration"
    );
    // Bare point is unchanged by the UV normalization choice and preserves the
    // exact raised-line continuation's earlier split-mass check.
    assert!((values["thermal_bare"][0] / 2.0918100443146576e-5 - 1.0).abs() < 1e-8);
    // Probe the integrated addback separately: mu_r changes only its finite
    // logarithm; the ordinary local subtraction is independent of that scale.
    cli.run_command("set process -p thermal -i NLO kv general.mu_r=2.0")?;
    let (process_id, integrand_name) = cli.state.find_integrand_ref(
        Some(&ProcessRef::Unqualified("thermal".into())),
        Some(&"NLO".into()),
    )?;
    let changed = EvaluateSamples {
        process_id: Some(process_id),
        integrand_name: Some(integrand_name),
        use_arb_prec: false,
        minimal_output: true,
        return_generated_events: Some(false),
        momentum_space: true,
        points: points.view(),
        integrator_weights: None,
        discrete_dims: Some(discrete.view()),
        graph_names: None,
        orientations: None,
    }
    .run(&mut cli.state)?;
    assert!(
        changed
            .samples
            .iter()
            .zip(&values["thermal"])
            .any(
                |(sample, original)| (sample.evaluation.integrand_result.re.0 - original).abs()
                    > 1e-9
            )
    );
    clean_test(cli.cli_settings.state.folder);
    Ok(())
}
