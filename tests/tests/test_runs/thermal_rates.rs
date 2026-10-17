use super::*;
use std::f64::consts::PI;

#[test]
fn thermal_scalar_sunset_matches_phase_space_and_conserves_events() -> Result<()> {
    let _serial_guard = gammaloop_integration_tests::SERIAL_TEST_LOCK.lock();
    let state_path =
        get_tests_workspace_path().join(format!("thermal_scalar_sunset_{}", std::process::id()));
    let mut cli = get_example_cli(
        "thermal_rates/scalar_sunset/run.toml",
        &["generate", "rest"],
        Some(state_path),
        None,
        true,
    )?;
    // The example helper saves the generated state. Exercise the persisted
    // signed cut catalogue as well as freshly generated expressions.
    cli.state =
        gammaloop_api::state::State::load(cli.cli_settings.state.folder.clone(), None, None)?;
    cli.state.activate_loaded_integrand_backends(false)?;

    let evaluate =
        |cli: &mut gammaloop_integration_tests::CLIState, process: &str, retain_events: bool| {
            let (process_id, integrand_name) = cli.state.find_integrand_ref(
                Some(&ProcessRef::Unqualified(process.to_owned())),
                Some(&"LO".to_owned()),
            )?;
            let points = arr2(&[[0.6, 0.0, 0.0, 0.0, 0.8, 0.0]]);
            let discrete = ndarray::Array2::<usize>::zeros((1, 0));
            evaluate_sample(
                &mut cli.state,
                &EvaluateSamples {
                    process_id: Some(process_id),
                    integrand_name: Some(integrand_name),
                    use_arb_prec: false,
                    minimal_output: false,
                    return_generated_events: Some(retain_events),
                    momentum_space: true,
                    points: points.view(),
                    integrator_weights: None,
                    discrete_dims: Some(discrete.view()),
                    graph_names: None,
                    orientations: None,
                },
            )
        };

    let mut rest_values = Vec::new();
    for block in ["rest", "boosted", "fast"] {
        cli.run_command(&format!("run {block}"))?;
        for process in ["equilibrium", "vacuum"] {
            let result = evaluate(&mut cli, process, true)?;
            let evaluation = &result.sample.evaluation;
            let events = evaluation
                .event_groups
                .iter()
                .flat_map(|group| group.iter())
                .collect_vec();
            assert!(!events.is_empty(), "{block}/{process}: missing events");
            if block == "rest" {
                assert_eq!(events.len(), if process == "equilibrium" { 4 } else { 1 });
                rest_values.push(evaluation.integrand_result.re.0);
            }
            assert!(evaluation.integrand_result.re.0.is_finite());
            assert!(evaluation.integrand_result.re.0 > 0.0);
            assert!(evaluation.integrand_result.im.0.abs() < 1e-12);
            let mut weight = Complex::new(0.0, 0.0);
            for event in events {
                let incoming = &event.kinematic_configuration.0;
                let outgoing = &event.kinematic_configuration.1;
                assert_eq!(event.cut_info.particle_pdgs.0[0], 9000001);
                assert!(matches!((incoming.len(), outgoing.len()), (1, 3) | (2, 2)));
                let mut residual = [0.0_f64; 4];
                for (sign, momenta) in [(1.0, incoming.as_slice()), (-1.0, outgoing.as_slice())] {
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
                assert!(
                    residual.iter().all(|x| x.abs() < 1e-10),
                    "{block}/{process}: nonconserving event {event:?}"
                );
                weight += Complex::new(event.weight.re.0, event.weight.im.0);
            }
            let factor = evaluation.integrator_weight.0
                * evaluation.parameterization_jacobian.map_or(1.0, |j| j.0);
            let expected = evaluation.integrand_result.re.0 * factor;
            assert!(
                (weight.re - expected).abs() <= 1e-10 * expected.abs(),
                "{block}/{process}: events do not reproduce the integrand"
            );
        }
    }

    // Independent momentum-space phase-space density. At rest the three
    // massless lines have lengths 0.6, 0.8, 1.0 in any sunset LMB. Each physical
    // channel has t = 1 / sum(sigma_i |q_i|). No generated residue enters here.
    let mut spectral_density = 0.0;
    let mut vacuum_density = 0.0;
    for absorbed in [None, Some(0), Some(1), Some(2)] {
        let raw_energies = [0.6_f64, 0.8, 1.0];
        let derivative: f64 = raw_energies
            .iter()
            .enumerate()
            .map(|(i, e)| if absorbed == Some(i) { -e } else { *e })
            .sum();
        let radius = 1.0 / derivative;
        let energies = raw_energies.map(|e| radius * e);
        let h = 2.0 / PI.sqrt() * (2.0 - radius * radius - 1.0 / (radius * radius)).exp();
        let measure = 2.0 * PI * radius.powi(6) * h
            / (derivative.abs()
                * energies.iter().map(|e| 2.0 * e).product::<f64>()
                * (2.0 * PI).powi(6)
                * 2.0);
        let occupations = energies.map(|e| 1.0 / e.exp_m1());
        let mut forward = 1.0;
        let mut reverse = 1.0;
        for (i, n) in occupations.into_iter().enumerate() {
            if absorbed == Some(i) {
                forward *= n;
                reverse *= 1.0 + n;
            } else {
                forward *= 1.0 + n;
                reverse *= n;
            }
        }
        spectral_density += measure * (forward - reverse);
        if absorbed.is_none() {
            vacuum_density = measure;
        }
    }
    for (actual, expected) in rest_values
        .into_iter()
        .zip([spectral_density, vacuum_density])
    {
        assert!(
            (actual - expected).abs() <= 1e-9 * expected,
            "sunset phase-space oracle: {actual:e}, expected {expected:e}"
        );
    }

    // The supported finite-temperature mode must approach the vacuum result.
    cli.run_command("run rest")?;
    cli.run_command("set process -p equilibrium -i LO kv general.inverse_temperature=1000.0")?;
    let cold = evaluate(&mut cli, "equilibrium", true)?
        .sample
        .evaluation
        .integrand_result;
    let vacuum = evaluate(&mut cli, "vacuum", true)?
        .sample
        .evaluation
        .integrand_result;
    assert!((cold.re.0 - vacuum.re.0).abs() <= 1e-9 * vacuum.re.0.abs());
    // The shared selector must see crossed incoming particles before counting
    // outgoing B1s. It removes exactly the B1-absorption channel, while the
    // vacuum decay survives. Event retention must not change that selection.
    cli.run_command("set process -p equilibrium -i LO kv general.inverse_temperature=1.0")?;
    for process in ["equilibrium", "vacuum"] {
        let baseline = evaluate(&mut cli, process, true)?.sample.evaluation;
        assert_eq!(baseline.event_groups.len(), 1);
        let expected = baseline
            .event_groups
            .iter()
            .flat_map(|group| group.iter())
            .filter(|event| event.cut_info.particle_pdgs.1.contains(&9000002))
            .map(|event| event.weight.re.0)
            .sum::<f64>();
        assert!(expected > 0.0);
        cli.run_command(&format!(
            r#"set process -p {process} -i LO string '
[quantities.outgoing_b1_count]
type = "particle"
computation = "count"
pdgs = [9000002]

[selectors.require_outgoing_b1]
quantity = "outgoing_b1_count"
selector = "discrete_range"
min = 1
max = 1
'"#
        ))?;
        for retain_events in [true, false] {
            let selected = evaluate(&mut cli, process, retain_events)?
                .sample
                .evaluation;
            let factor = selected.integrator_weight.0
                * selected.parameterization_jacobian.map_or(1.0, |j| j.0);
            assert!(
                (selected.integrand_result.re.0 * factor - expected).abs()
                    <= 1e-10 * expected.abs(),
                "{process}: crossing-aware selector changed the selected event weight"
            );
            if retain_events {
                assert_eq!(selected.event_groups.len(), 1);
                let events = selected
                    .event_groups
                    .iter()
                    .flat_map(|group| group.iter())
                    .collect_vec();
                assert_eq!(events.len(), if process == "equilibrium" { 3 } else { 1 });
                assert!(
                    events
                        .iter()
                        .all(|event| event.cut_info.particle_pdgs.1.contains(&9000002))
                );
            } else {
                assert!(selected.event_groups.is_empty());
                let metadata = selected
                    .evaluation_metadata
                    .as_ref()
                    .expect("selectors retain evaluation diagnostics");
                assert!(metadata.accepted_event_count > 0);
                if process == "equilibrium" {
                    assert!(metadata.generated_event_count > metadata.accepted_event_count);
                }
            }
        }
    }
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
